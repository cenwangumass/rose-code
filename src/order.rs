//! Barrier setup structure: the `A A'` pattern, the ordering levels, the minimum-degree ordering and the supernodal symbolic
//! factorization of the reference (`reverse-engineer/specs/barrier.md`; port of `tools/barrier_order.py`).
//!
//! Conventions: `perm[original] = position`. A [`Pattern`] is the strictly off-diagonal symmetric pattern of `A A'` over the
//! model's rows in file order, with each dense column an extra node after the rows.
//!
//! Known gaps, reported as errors (never approximated):
//! - the partition routine `0x325850` for connected components larger than `max(600, n/100)` (PEROLD, 80BAU3B, DEGEN3):
//!   [`levels`] returns [`StructError::Partition`] unless the caller supplies the labels;
//! - the quicksort partition of the sorter `0x8c8b70` for lists of 100 or more: [`StructError::SortLarge`] (not reached by the 12 models);
//! - the dense-column threshold is not recovered: [`dense_columns`] takes the minimum column count as a parameter.
use crate::model::LpModel;
use std::collections::VecDeque;

#[derive(Debug, Clone, PartialEq)]
pub enum StructError {
    /// A connected component of this many nodes needs the partition routine `0x325850`, which is not decoded.
    Partition(usize),
    /// The sorter `0x8c8b70` meets a list of 100 or more keys (quicksort partition, not decoded).
    SortLarge(usize),
}

impl std::fmt::Display for StructError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            StructError::Partition(n) => write!(f, "partition routine 0x325850 needed for a component of {n} nodes (not decoded)"),
            StructError::SortLarge(n) => write!(f, "sorter 0x8c8b70 quicksort partition needed for {n} keys (not decoded)"),
        }
    }
}

impl std::error::Error for StructError {}

#[derive(Debug, Clone, PartialEq)]
pub struct Pattern {
    pub n: usize,
    pub colptr: Vec<usize>,
    pub rowind: Vec<usize>,
}

impl Pattern {
    pub fn nnz(&self) -> usize {
        self.colptr[self.n]
    }
    fn adj(&self, i: usize) -> &[usize] {
        &self.rowind[self.colptr[i]..self.colptr[i + 1]]
    }
}

#[derive(Debug, Clone, PartialEq)]
pub struct OrderResult {
    pub perm: Vec<i64>,
    pub nz: i64,
    pub ops: f64,
}

#[derive(Debug, Clone, PartialEq)]
pub struct Symbolic {
    /// Exact nnz(L) with the diagonal in the ordering's numbering (the value of `0x31ebb0`).
    pub exact_nz: i64,
    pub perm: Vec<i64>,
    pub nz: i64,
    pub ops: f64,
    pub width: Vec<i64>,
    pub count: Vec<i64>,
    pub rowoff: Vec<i64>,
    pub idx: Vec<i64>,
}

/// Constraint columns' row lists (explicit zero coefficients dropped), in the model's column order.
fn columns(m: &LpModel) -> Vec<Vec<usize>> {
    (0..m.num_cols()).map(|j| m.col(j).filter(|&(_, v)| v != 0.0).map(|(i, _)| i).collect()).collect()
}

/// The columns taken as dense, ascending. The reference's threshold (`P+0x218`) is not recovered; the data only bound it
/// (49 nonzeros is not dense, 389 is), so `min_nnz` is a parameter and the default [`DENSE_COL_MIN`] is a placeholder inside that range.
pub fn dense_columns(m: &LpModel, min_nnz: usize) -> Vec<usize> {
    columns(m).iter().enumerate().filter(|(_, c)| c.len() >= min_nnz).map(|(j, _)| j).collect()
}

/// Placeholder for the unrecovered dense-column threshold (between 50 and 389 nonzeros on the captured models).
pub const DENSE_COL_MIN: usize = 100;

/// The `AA'` pattern the ordering gets (`PAT` of the captures): rows in file order, dense columns left out of the `AA'` part and
/// added as nodes `rows + k` adjacent to the rows of their column. Returns the pattern and `n1` (= rows).
pub fn aat_pattern(m: &LpModel, dense: &[usize]) -> (Pattern, usize) {
    let nr = m.num_rows();
    let cols = columns(m);
    let mut is_dense = vec![false; cols.len()];
    for &j in dense {
        is_dense[j] = true;
    }
    let mut adj: Vec<Vec<usize>> = vec![Vec::new(); nr];
    for (j, c) in cols.iter().enumerate() {
        if is_dense[j] {
            continue;
        }
        for &a in c {
            adj[a].extend(c.iter().copied());
        }
    }
    for (i, a) in adj.iter_mut().enumerate() {
        a.sort_unstable();
        a.dedup();
        a.retain(|&x| x != i);
    }
    let mut ext = Vec::new();
    for (k, &j) in dense.iter().enumerate() {
        for &r in &cols[j] {
            adj[r].push(nr + k);
        }
        let mut e = cols[j].clone();
        e.sort_unstable();
        ext.push(e);
    }
    for a in adj.iter_mut() {
        a.sort_unstable();
    }
    adj.extend(ext);
    let mut colptr = vec![0];
    let mut rowind = Vec::new();
    for a in &adj {
        rowind.extend_from_slice(a);
        colptr.push(rowind.len());
    }
    (Pattern { n: nr + dense.len(), colptr, rowind }, nr)
}

/// `max(100, int(min(int(50 nnz / n), 10 sqrt(n))))`: the degree above which a node is dense (driver `0x324330` and core `0x320fa0`).
pub fn dense_threshold(n: usize, nnz: usize) -> i64 {
    let a = (50.0 * nnz as f64 / n as f64) as i64 as f64;
    let b = 10.0 * (n as f64).sqrt();
    100.max(a.min(b) as i64)
}

/// Labels of the partition routine `0x325850` for a component: `bisect(size, colptr, rowind)`, label 2 = vertex separator.
pub type Bisect<'a> = &'a mut dyn FnMut(usize, &[usize], &[usize]) -> Vec<i64>;

/// Driver `0x324330` (mode 3): the level flags handed to the core (-1 normal, -2 later, -3 dense-column nodes). A component of
/// the -1 nodes larger than `max(600, n/100)` needs the partition routine: without `bisect` this is [`StructError::Partition`].
pub fn levels(p: &Pattern, n1: usize, bisect: Option<Bisect>) -> Result<Vec<i64>, StructError> {
    let n = p.n;
    let n1 = n1.min(n);
    let thr = dense_threshold(n, p.nnz());
    let mut lv: Vec<i64> = (0..n)
        .map(|i| if i < n1 { if (p.adj(i).len() as i64) > thr { -2 } else { -1 } } else { -3 })
        .collect();
    let mut comp = vec![usize::MAX; n];
    let mut sizes = Vec::new();
    for st in 0..n {
        if lv[st] == -1 && comp[st] == usize::MAX {
            let c = sizes.len();
            comp[st] = c;
            let mut stack = vec![st];
            let mut cnt = 1;
            while let Some(u) = stack.pop() {
                for &x in p.adj(u) {
                    if lv[x] == -1 && comp[x] == usize::MAX {
                        comp[x] = c;
                        stack.push(x);
                        cnt += 1;
                    }
                }
            }
            sizes.push(cnt);
        }
    }
    let big = 600.max(n / 100);
    let mut bisect = bisect;
    for (c, &sz) in sizes.iter().enumerate() {
        if sz > big {
            let Some(b) = bisect.as_mut() else { return Err(StructError::Partition(sz)) };
            let nodes: Vec<usize> = (0..n).filter(|&i| comp[i] == c).collect();
            let mut loc = vec![usize::MAX; n];
            for (k, &v) in nodes.iter().enumerate() {
                loc[v] = k;
            }
            let mut cp = vec![0];
            let mut ri = Vec::new();
            for &v in &nodes {
                ri.extend(p.adj(v).iter().filter(|&&x| loc[x] != usize::MAX).map(|&x| loc[x]));
                cp.push(ri.len());
            }
            let lab = b(sz, &cp, &ri);
            for (k, &v) in nodes.iter().enumerate() {
                if lab[k] == 2 {
                    lv[v] = -2;
                }
            }
        }
    }
    Ok(lv)
}

/// `0x8c8b70` for fewer than 100 keys: shell sort, gap `n/2` then `int(gap/2.2)` (a gap of 2 is followed by 1), insertion with a strict `>` shift.
fn shellsort(idx: &mut [i64], key: &[f64]) -> Result<(), StructError> {
    let n = idx.len();
    if n >= 100 {
        return Err(StructError::SortLarge(n));
    }
    let mut gap = n / 2;
    while gap > 0 {
        let gaps: &[usize] = if gap == 2 { &[2, 1] } else { &[gap][..] };
        for &g in gaps {
            for i in g..n {
                let x = idx[i];
                let kx = key[x as usize];
                let mut j = i;
                while j >= g {
                    let y = idx[j - g];
                    if key[y as usize] <= kx {
                        break;
                    }
                    idx[j] = y;
                    j -= g;
                }
                idx[j] = x;
            }
        }
        gap = (gap as f64 / 2.2) as usize;
    }
    Ok(())
}

/// The ordering core `0x320fa0`: an approximate-minimum-degree variant on the quotient graph. `levels` are the input flags
/// (negative; -1 first). Returns `perm[original] = elimination position`, the core's off-diagonal count and Factor Ops.
pub fn order(p: &Pattern, levels: &[i64]) -> Result<OrderResult, StructError> {
    let n = p.n;
    let ni = n as i64;
    let mut adj: Vec<Vec<i64>> = (0..n).map(|i| p.adj(i).iter().map(|&x| x as i64).collect()).collect();
    let mut ln: Vec<i64> = adj.iter().map(|a| a.len() as i64).collect();
    let mut elen = vec![0i64; n];
    let mut deg = ln.clone();
    let mut nv = vec![1i64; n];
    let mut parent: Vec<usize> = (0..n).collect();
    let mut perm: Vec<i64> = levels.to_vec();
    let maxlev = perm.iter().filter(|&&x| x < 0).map(|&x| -x).max().unwrap_or(0);
    let mut head = vec![-1i64; n + 1];
    let mut prev = vec![-1i64; n];
    let mut nxt = vec![-1i64; n];
    let mut w = vec![-1i64; n];
    let mut hsh = vec![0f64; n];
    let mut nz: i64 = 0;
    let mut ops = 0f64;
    let mut nelim: i64 = 0;
    let thr = if n > 0 { dense_threshold(n, p.nnz()) } else { 100 };
    let mut ndense = 0;
    for i in 0..n {
        if ln[i] > thr {
            ndense += 1;
            perm[i] = ni - perm[i];
            ln[i] = -1;
        }
    }
    if ndense > 0 {
        for i in 0..n {
            let l = ln[i];
            if l > 0 {
                adj[i] = adj[i][..l as usize].iter().copied().filter(|&j| perm[j as usize] < 0).collect();
                ln[i] = adj[i].len() as i64;
            } else {
                ln[i] = 0;
                adj[i] = Vec::new();
            }
        }
    }
    for i in 0..n {
        if perm[i] == -1 && ln[i] == 0 {
            perm[i] = nelim;
            nelim += 1;
            ln[i] = -1;
            ops += 1.0;
        }
    }
    macro_rules! link {
        ($v:expr) => {{
            let v = $v as usize;
            let d = deg[v] as usize;
            nxt[v] = head[d];
            if head[d] != -1 {
                prev[head[d] as usize] = v as i64;
            }
            prev[v] = -1;
            head[d] = v as i64;
        }};
    }
    macro_rules! unlink {
        ($v:expr) => {{
            let v = $v as usize;
            let (pp, qq) = (prev[v], nxt[v]);
            if pp == -1 {
                head[deg[v] as usize] = qq;
            } else {
                nxt[pp as usize] = qq;
            }
            if qq != -1 {
                prev[qq as usize] = pp;
            }
        }};
    }
    while nelim < ni {
        let mut level = maxlev;
        for i in 0..n {
            if perm[i] < 0 && nv[i] > 0 {
                level = level.min(-perm[i]);
            } else if perm[i] > ni && nv[i] > 0 {
                level = level.min(perm[i] - ni);
            }
        }
        let mut mind = ni;
        for v in (0..n).rev() {
            if perm[v] == -level && nv[v] > 0 {
                link!(v);
                mind = mind.min(deg[v]);
            }
        }
        loop {
            let mut d = mind;
            while d < ni && head[d as usize] == -1 {
                d += 1;
            }
            if d >= ni {
                break;
            }
            let v = head[d as usize] as usize;
            unlink!(v);
            perm[v] = nelim;
            w[v] = v as i64;
            let mut lme: Vec<usize> = Vec::new();
            let mut degme: i64 = 0;
            let elv = elen[v] as usize;
            let a = adj[v].clone();
            for &e in &a[..elv] {
                let e = e as usize;
                for &j in &adj[e][..ln[e] as usize] {
                    let j = j as usize;
                    if w[j] != v as i64 && nv[j] > 0 && perm[j] < 0 {
                        lme.push(j);
                        w[j] = v as i64;
                        degme += nv[j];
                    }
                }
                elen[e] = -1;
                ln[e] = -1;
            }
            for &j in &a[elv..ln[v] as usize] {
                let j = j as usize;
                if w[j] != v as i64 && nv[j] > 0 && perm[j] < 0 {
                    lme.push(j);
                    w[j] = v as i64;
                    degme += nv[j];
                }
            }
            for &i in &lme {
                for &e in &adj[i][..elen[i] as usize] {
                    deg[e as usize] = nv[e as usize];
                }
            }
            for &i in &lme {
                for &e in &adj[i][..elen[i] as usize] {
                    let e = e as usize;
                    deg[e] -= nv[i];
                    if deg[e] == 0 {
                        ln[e] = 0;
                    }
                }
            }
            for &i in &lme {
                if perm[i] == -level {
                    unlink!(i);
                }
            }
            for &i in &lme {
                let mut s: i64 = 0;
                let mut h = (perm[i] * ni) as f64;
                let mut na: Vec<i64> = Vec::new();
                let ei = elen[i] as usize;
                for &e in &adj[i][..ei] {
                    if ln[e as usize] > 0 {
                        s += deg[e as usize];
                        h += e as f64;
                        na.push(e);
                    }
                }
                na.push(v as i64);
                let ne = na.len() as i64;
                for &j in &adj[i][ei..ln[i] as usize] {
                    if v as i64 != w[j as usize] && nv[j as usize] > 0 {
                        s += nv[j as usize];
                        na.push(j);
                        h += j as f64;
                    }
                }
                ln[i] = na.len() as i64;
                adj[i] = na;
                elen[i] = ne;
                let dd = degme + deg[i] - nv[i];
                deg[i] = s.min(dd);
                hsh[i] = h;
                if 5 * ni - 5 * nelim < s + nv[i] {
                    perm[i] = ni - perm[i];
                    ln[i] = -1;
                }
            }
            let mut keep: Vec<i64> = Vec::new();
            for &i in &lme {
                if perm[i] == -level && ln[i] == 1 {
                    nv[v] += nv[i];
                    degme -= nv[i];
                    nv[i] = 0;
                    parent[i] = v;
                    ln[i] = -1;
                } else if perm[i] < ni {
                    keep.push(i as i64);
                }
            }
            let k = nv[v];
            nz += k * degme + k * (k - 1) / 2;
            let s = degme + k;
            let kk = (k * (k - 1)) as f64;
            ops = (2 * k - 1) as f64 * kk / 6.0 + ((k as f64 * s as f64 * s as f64 + ops) - s as f64 * kk);
            nelim += k;
            nv[v] = degme;
            for &i in &keep {
                let i = i as usize;
                deg[i] = (ni - 1).min(degme - nv[i] + deg[i]);
            }
            let mut lst = keep;
            shellsort(&mut lst, &hsh)?;
            let m = lst.len();
            let mut pos = 0;
            while pos < m {
                let start = pos;
                pos = start + 1;
                while pos < m && hsh[lst[start] as usize] == hsh[lst[pos] as usize] {
                    pos += 1;
                }
                for a0 in start..pos {
                    for b0 in a0 + 1..pos {
                        let (i, j) = (lst[a0], lst[b0]);
                        if i >= 0 && j >= 0 {
                            let (i, j) = (i as usize, j as usize);
                            if perm[i] == perm[j]
                                && ln[i] == ln[j]
                                && elen[i] == elen[j]
                                && adj[i][..ln[i] as usize] == adj[j][..ln[j] as usize]
                            {
                                nv[i] += nv[j];
                                deg[i] -= nv[j];
                                nv[j] = 0;
                                parent[j] = i;
                                ln[j] = -1;
                                lst[b0] = -1;
                            }
                        }
                    }
                }
            }
            lst.retain(|&x| x >= 0);
            ln[v] = lst.len() as i64;
            elen[v] = 0;
            mind = d;
            for &i in &lst {
                if perm[i as usize] == -level {
                    link!(i);
                    mind = mind.min(deg[i as usize]);
                }
            }
            adj[v] = lst;
        }
        for i in 0..n {
            if perm[i] == ni + level {
                perm[i] = nelim;
                ln[i] = -1;
                nelim += nv[i];
                nz += ni - 1;
                nv[i] = 0;
            }
        }
    }
    for i in 0..n {
        if perm[i] < 0 {
            let mut r = parent[i];
            while parent[r] != r {
                r = parent[r];
            }
            perm[i] = perm[r];
            perm[r] += 1;
        }
    }
    Ok(OrderResult { perm, nz, ops })
}

/// Symbolic factorization `0x2fedb0`: relaxed supernodes, the final permutation, Factor NZ/Ops and the row lists.
/// `perm1` is the core's permutation; `n1` = rows (columns `>= n1` are dense-column nodes, never merged with rows).
pub fn symbolic(p: &Pattern, perm1: &[i64], n1: usize) -> Symbolic {
    let n = p.n;
    let mut inv = vec![0usize; n];
    for i in 0..n {
        inv[perm1[i] as usize] = i;
    }
    // 1. exact structure of L in the perm1 numbering
    let mut rows: Vec<Vec<usize>> = vec![Vec::new(); n];
    let mut kids: Vec<Vec<usize>> = vec![Vec::new(); n];
    for j in 0..n {
        let mut s: Vec<usize> = p.adj(inv[j]).iter().map(|&i| perm1[i] as usize).filter(|&x| x > j).collect();
        for &c in &kids[j] {
            s.extend(rows[c][1..].iter().copied().filter(|&x| x > j));
        }
        s.sort_unstable();
        s.dedup();
        if let Some(&mn) = s.first() {
            kids[mn].push(j);
        }
        let mut r = vec![j];
        r.extend(s);
        rows[j] = r;
    }
    let cnt: Vec<usize> = rows.iter().map(|r| r.len()).collect();
    let exact_nz: usize = cnt.iter().sum();
    // 2. fundamental supernodes
    let mut first = Vec::new();
    let mut width: Vec<i64> = Vec::new();
    let mut sn_of = vec![0usize; n];
    let mut j = 0;
    while j < n {
        let mut k = 1;
        while j + k < n && cnt[j + k] + k == cnt[j] && k < rows[j].len() && rows[j][k] == j + k && (n1 <= j) == (n1 <= j + k) {
            k += 1;
        }
        for t in 0..k {
            sn_of[j + t] = first.len();
        }
        first.push(j);
        width.push(k as i64);
        j += k;
    }
    let nsn = first.len();
    let mut scnt: Vec<i64> = first.iter().map(|&f| cnt[f] as i64).collect();
    let mut members: Vec<VecDeque<usize>> = vec![VecDeque::new(); nsn];
    for c in 0..n {
        members[sn_of[c]].push_front(c);
    }
    let mut par = vec![nsn; nsn + 1];
    for k in 0..nsn {
        if width[k] != scnt[k] {
            par[k] = sn_of[rows[first[k]][width[k] as usize]];
        }
    }
    let mut ch: Vec<Vec<usize>> = vec![Vec::new(); nsn + 1];
    for k in 0..nsn {
        ch[par[k]].insert(0, k);
    }
    let f = |x: i64| -> f64 {
        let x = x as f64;
        x * (x + 1.0) * (2.0 * x + 1.0) / 6.0
    };
    // 3. relaxed merging
    for c in 0..nsn {
        let pp = par[c];
        if pp >= nsn {
            continue;
        }
        let (wc, cc, cp) = (width[c], scnt[c], scnt[pp]);
        if cc <= 1 || width[c] <= 0 {
            continue;
        }
        if (n1 <= members[c][0]) != (n1 <= members[pp][0]) {
            continue;
        }
        let r = (cc - wc) as f64;
        let t = if wc < 192 { 3.0 } else { 0.1 };
        let extra = (f(cp + wc) - f(cp)) - (f(cc) - f(cc - wc));
        if extra < t * r * r + 1000.0 {
            let mem = std::mem::take(&mut members[c]);
            for x in mem {
                members[pp].push_front(x);
                sn_of[x] = pp;
            }
            let at = ch[pp].iter().position(|&x| x == c).unwrap();
            ch[pp].remove(at);
            let cc_children = std::mem::take(&mut ch[c]);
            for x in cc_children {
                par[x] = pp;
                ch[pp].insert(0, x);
            }
            scnt[pp] += wc;
            width[pp] += wc;
            width[c] = -1;
        }
    }
    // 4. pass 2: singleton siblings with identical structure (transcribed from the decompile; never exercised by the 12 models)
    for k in 0..nsn {
        if width[k] <= 0 {
            continue;
        }
        let mut lst = std::mem::take(&mut ch[k]);
        let mut i = 0;
        while i < lst.len() {
            let c = lst[i];
            if width[c] == 1 && scnt[c] > 24 && i + 1 < lst.len() {
                let a = members[c][0];
                let mut j2 = i + 1;
                while j2 < lst.len() {
                    if width[c] >= 24 {
                        break;
                    }
                    let s2 = lst[j2];
                    if width[s2] == 1 {
                        let b = members[s2][0];
                        if cnt[b] == cnt[a] && (n1 <= b) == (n1 <= a) && rows[a][1..] == rows[b][1..] {
                            scnt[c] += 1;
                            width[c] += 1;
                            members[c].push_front(b);
                            members[s2].clear();
                            let kids2 = std::mem::take(&mut ch[s2]);
                            for x in kids2 {
                                par[x] = c;
                                ch[c].insert(0, x);
                            }
                            width[s2] = -1;
                            lst.remove(j2);
                            continue;
                        }
                    }
                    j2 += 1;
                }
            }
            i += 1;
        }
        ch[k] = lst;
    }
    // 5. postorder: preorder from the virtual root with a stack, reversed
    let mut stack: Vec<usize> = ch[nsn].clone();
    let mut pre = Vec::new();
    while let Some(x) = stack.pop() {
        pre.push(x);
        stack.extend(ch[x].iter().copied());
    }
    pre.reverse();
    let mut post: Vec<usize> = pre.iter().copied().filter(|&x| members[x][0] < n1).collect();
    post.extend(pre.iter().copied().filter(|&x| members[x][0] >= n1));
    let mut order = Vec::new();
    let mut xmap = vec![0usize; n];
    for &s in &post {
        order.extend(members[s].iter().copied());
    }
    for (pos, &c) in order.iter().enumerate() {
        xmap[c] = pos;
    }
    let perm: Vec<i64> = (0..n).map(|i| xmap[perm1[i] as usize] as i64).collect();
    let (mut width_out, mut count, mut rowoff) = (vec![0i64; n], vec![0i64; n], vec![0i64; n]);
    let mut idx: Vec<i64> = Vec::new();
    let (mut nz, mut ops, mut pos) = (0i64, 0f64, 0usize);
    for &s in &post {
        let m = &members[s];
        let mx = *m.iter().max().unwrap();
        let w_ = m.len();
        let mut rl: Vec<usize> = m.iter().map(|&c| xmap[c]).collect();
        rl.extend(rows[mx][1..].iter().map(|&r| xmap[r]));
        rl.sort_unstable();
        let base = idx.len() as i64;
        idx.extend(rl.iter().map(|&x| x as i64));
        let rr = rl.len() as i64;
        for t in 0..w_ {
            let tt = t as i64;
            count[pos + t] = rr - tt;
            rowoff[pos + t] = base + tt;
            width_out[pos + t] = if t == 0 { w_ as i64 } else { -tt };
            nz += rr - tt;
            ops += ((rr - tt) as f64) * ((rr - tt) as f64);
        }
        pos += w_;
    }
    Symbolic { exact_nz: exact_nz as i64, perm, nz, ops, width: width_out, count, rowoff, idx }
}

/// The whole structure stage for a model: pattern, levels, ordering, symbolic factorization.
#[derive(Debug, Clone)]
pub struct Structure {
    pub pattern: Pattern,
    pub n1: usize,
    pub dense: Vec<usize>,
    pub levels: Vec<i64>,
    pub order: OrderResult,
    pub sym: Symbolic,
}

/// Runs the stage with the dense columns chosen by [`dense_columns`]`(m, dense_min)`. Errors name the undecoded routine reached.
pub fn structure(m: &LpModel, dense_min: usize) -> Result<Structure, (Pattern, usize, Vec<usize>, StructError)> {
    let dense = dense_columns(m, dense_min);
    let (pattern, n1) = aat_pattern(m, &dense);
    let run = || -> Result<(Vec<i64>, OrderResult, Symbolic), StructError> {
        let lv = levels(&pattern, n1, None)?;
        let o = order(&pattern, &lv)?;
        let s = symbolic(&pattern, &o.perm, n1);
        Ok((lv, o, s))
    };
    match run() {
        Ok((levels, order, sym)) => Ok(Structure { pattern, n1, dense, levels, order, sym }),
        Err(e) => Err((pattern, n1, dense, e)),
    }
}
