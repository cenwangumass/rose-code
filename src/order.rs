//! Barrier setup structure: the `A A'` pattern, the ordering levels, the minimum-degree ordering and the supernodal symbolic
//! factorization of the reference (`reverse-engineer/specs/barrier.md`; port of `tools/barrier_order.py`).
//!
//! Conventions: `perm[original] = position`. A [`Pattern`] is the strictly off-diagonal symmetric pattern of `A A'` over the
//! model's rows in file order, with each dense column an extra node after the rows.
//!
use crate::model::LpModel;
use std::collections::VecDeque;

#[derive(Debug, Clone, PartialEq)]
pub enum StructError {
    /// Reserved for compatibility; the dense-column scan is implemented.
    DenseScan,
}

impl std::fmt::Display for StructError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            StructError::DenseScan => {
                write!(f, "dense-column scan error")
            }
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
    /// Work units the core adds to the manager's counter (`0x323507`, `0x3231f5`, `0x323827`): `5.25 elen_i` and `3.75 len_i` per variable of every pivot
    /// list, `15 |kept|` per pivot (the ordering manager `0x31da80` compares the spent work with its estimate).
    pub work: f64,
}

#[derive(Debug, Clone, PartialEq)]
pub struct Symbolic {
    /// Exact nnz(L) with the diagonal in the ordering's numbering (the value of `0x31ebb0`).
    pub exact_nz: i64,
    /// Sum of the squared exact column counts (`0x31f620`: the ordering driver replaces the core's Ops by it; the manager `0x31da80` budgets against it).
    pub exact_ops: f64,
    pub perm: Vec<i64>,
    pub nz: i64,
    pub ops: f64,
    pub width: Vec<i64>,
    pub count: Vec<i64>,
    pub rowoff: Vec<i64>,
    pub idx: Vec<i64>,
    /// State of the reference's generator (`ctx+0xb0`, `0x88a870`) after the ordering: the early-centering solve draws its start vector from it.
    pub rng: u64,
}

/// Constraint columns' row lists (explicit zero coefficients dropped), in the model's column order.
fn columns(m: &LpModel) -> Vec<Vec<usize>> {
    (0..m.num_cols())
        .map(|j| {
            m.col(j)
                .filter(|&(_, v)| v != 0.0)
                .map(|(i, _)| i)
                .collect()
        })
        .collect()
}

/// The dense columns of the model conversion `0x31af70` (`P+0x218`; port of `dense_columns` of `tools/barrier_order.py`, `specs/barrier.md` "Dense columns").
/// With `nrows < 1000` there are none. Otherwise: the histogram of the column counts is walked from `d = 0` with `S` = sum of `d'^2 hist[d']` (double) over
/// `d' < d` and `T` = number of columns with count `>= d`; the first `d` with `hist[d] > 0`, `d^2 > S` and `d > 5 T` is the cut (none: `int(max(50 nnz /
/// max(ncols, 1), 500))`). Columns at or above the cut are dense. A remaining column with at least `int(max(10 nnz / max(ncols, 1), 100))` entries starts
/// the cost scan `0x31bc1a`, which can select further columns.
pub fn dense_rule(m: &LpModel) -> Result<Vec<usize>, StructError> {
    let nrows = m.num_rows();
    if nrows < 1000 {
        return Ok(Vec::new());
    }
    let cols = columns(m);
    let cnt: Vec<usize> = cols.iter().map(|c| c.len()).collect();
    let (nc, nnz) = (cnt.len(), cnt.iter().sum::<usize>());
    let mut hist = vec![0usize; nrows + 1];
    for &c in &cnt {
        hist[c] += 1;
    }
    let (mut cut, mut s, mut t) = (None, 0.0f64, nc);
    for d in 0..=nrows {
        let h = hist[d];
        let sq = d as f64 * d as f64;
        if h > 0 && sq > s && d as f64 > t as f64 * 5.0 {
            cut = Some(d);
            break;
        }
        s += sq * h as f64;
        t -= h;
    }
    let avg = nnz as f64 / nc.max(1) as f64;
    let cut = cut.unwrap_or_else(|| (50.0 * avg).max(500.0) as usize);
    let thr2 = (10.0 * avg).max(100.0) as usize;
    let mut dense: Vec<usize> = cnt
        .iter()
        .enumerate()
        .filter(|&(_, &c)| c >= cut)
        .map(|(j, _)| j)
        .collect();
    if cnt.iter().any(|&c| c >= thr2 && c < cut) {
        let mut ids: Vec<i64> = (0..nc)
            .filter(|&j| cnt[j] < cut)
            .map(|j| j as i64)
            .collect();
        let keys: Vec<f64> = cnt.iter().map(|&c| c as f64).collect();
        shellsort(&mut ids, &keys);
        let mut incident = vec![Vec::new(); nrows];
        for &j in &ids {
            for &r in &cols[j as usize] {
                incident[r].push(j as usize);
            }
        }
        // For each row, charge each column for its newly exposed neighbors;
        // earlier columns in sorted order own the neighbors they share.
        let mut seen = vec![usize::MAX; nrows];
        let mut cost = vec![0.0; nc];
        for r in 0..nrows {
            seen[r] = r;
            for &j in &incident[r] {
                let mut k = 0;
                for &s in &cols[j] {
                    if seen[s] != r {
                        seen[s] = r;
                        k += 1;
                        cost[j] += k as f64;
                    }
                }
            }
        }
        // Walk backward, resetting both totals whenever the cost clears the threshold.
        let (mut start, mut c, mut s) = (ids.len(), 0.0, 0.0);
        for i in (0..ids.len()).rev() {
            c += cost[ids[i] as usize];
            s += (nc + nrows - i) as f64;
            if c > 10.0 * s + 1000.0 {
                start = i;
                c = 0.0;
                s = 0.0;
            }
        }
        dense.extend(ids[start..].iter().map(|&j| j as usize));
    }
    Ok(dense)
}

/// The columns taken as dense, in dense-node order: the reference's rule ([`dense_rule`]) when `min_nnz` is [`DENSE_AUTO`], else every column with at least
/// `min_nnz` nonzeros (tests).
pub fn dense_columns(m: &LpModel, min_nnz: usize) -> Vec<usize> {
    if min_nnz == DENSE_AUTO {
        return dense_rule(m).unwrap_or_else(|e| panic!("{e}"));
    }
    columns(m)
        .iter()
        .enumerate()
        .filter(|(_, c)| c.len() >= min_nnz)
        .map(|(j, _)| j)
        .collect()
}

/// `min_nnz` value that selects the reference's dense-column rule.
pub const DENSE_AUTO: usize = 0;
/// Kept name of [`DENSE_AUTO`] for the callers that pass the default.
pub const DENSE_COL_MIN: usize = DENSE_AUTO;

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
    // Row-to-column incidence lets each row form its set before sorting,
    // avoiding duplicate neighbour storage from overlapping sparse columns.
    let mut offsets = vec![0usize; nr + 1];
    for (j, c) in cols.iter().enumerate() {
        if !is_dense[j] {
            for &r in c {
                offsets[r + 1] += 1;
            }
        }
    }
    for r in 0..nr {
        offsets[r + 1] += offsets[r];
    }
    let mut next = offsets[..nr].to_vec();
    let mut rowcols = vec![0usize; offsets[nr]];
    for (j, c) in cols.iter().enumerate() {
        if !is_dense[j] {
            for &r in c {
                rowcols[next[r]] = j;
                next[r] += 1;
            }
        }
    }
    // Reuse the completed fill cursors as generation markers.
    next.fill(usize::MAX);
    for r in 0..nr {
        for &j in &rowcols[offsets[r]..offsets[r + 1]] {
            for &other in &cols[j] {
                if other != r && next[other] != r {
                    next[other] = r;
                    adj[r].push(other);
                }
            }
        }
        adj[r].sort_unstable();
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
    (
        Pattern {
            n: nr + dense.len(),
            colptr,
            rowind,
        },
        nr,
    )
}

#[cfg(test)]
fn frozen_aat_pattern(m: &LpModel, dense: &[usize]) -> (Pattern, usize) {
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
    (
        Pattern {
            n: nr + dense.len(),
            colptr,
            rowind,
        },
        nr,
    )
}

/// `max(100, int(min(int(50 nnz / n), 10 sqrt(n))))`: the degree above which a node is dense (driver `0x324330` and core `0x320fa0`).
pub fn dense_threshold(n: usize, nnz: usize) -> i64 {
    let a = (50.0 * nnz as f64 / n as f64) as i64 as f64;
    let b = 10.0 * (n as f64).sqrt();
    100.max(a.min(b) as i64)
}

/// Labels of the partition routine `0x325850` for a component: `bisect(size, colptr, rowind)`, label 2 = vertex separator.
pub type Bisect<'a> = &'a mut dyn FnMut(usize, &[usize], &[usize]) -> Vec<i64>;

thread_local! {
    /// Work units the partition code adds to the ordering manager's counter since the last [`take_partition_work`].
    static PARTITION_WORK: std::cell::Cell<f64> = const { std::cell::Cell::new(0.0) };
}

thread_local! {
    /// Matching scheme of the partition routine (its last argument, `0xf0(rsp)` of `0x325850`): 3 = the first unmatched neighbour in adjacency order (`0x326eb3`), every other
    /// value (-1 for the first ordering, the task index 0..4 for the `BarOrder = -1` extra orderings) = the unmatched neighbour of smallest degree (`0x326e45`).
    static MATCH_SCHEME: std::cell::Cell<i32> = const { std::cell::Cell::new(-1) };
}

/// Sets the matching scheme of the partitions run on this thread (see [`partition`]); -1 is the default.
pub fn set_match_scheme(s: i32) {
    MATCH_SCHEME.with(|c| c.set(s));
}

/// Adds `(a deg + b) 0.75` work units (all values are multiples of 1/4, so the sum does not depend on the order).
fn add_work(deg: usize, a: f64, b: f64) {
    PARTITION_WORK.with(|w| w.set(w.get() + (a * deg as f64 + b) * 0.75));
}

/// Returns and clears the work of the partitions run on this thread: per fine vertex of the coarse-adjacency build `(3 deg + 2) 0.75`
/// (`0x326d88`), per vertex pulled into the separator by a move `(6 deg + 10) 0.75` (`0x327b1b`) and per one restored by the undo `(3 deg + 4) 0.75` (`0x327f57`).
pub fn take_partition_work() -> f64 {
    PARTITION_WORK.with(|w| w.replace(0.0))
}

/// The generator `0x88a870` (31 bits): `state = state * 0x5851f42d4c957f2d + 0x14057b7ef767814f`, tempering of the high word.
fn rand31(rng: &mut u64) -> u32 {
    *rng = rng
        .wrapping_mul(0x5851f42d4c957f2d)
        .wrapping_add(0x14057b7ef767814f);
    let mut e = (*rng >> 32) as u32;
    e ^= e >> 11;
    e ^= (e << 7) & 0x9d2c5680;
    e ^= (e << 15) & 0xefc60000;
    e ^= e >> 18;
    e & 0x7fffffff
}

/// Max-heap of items with integer keys (`0x887c20` ... `0x889290`): sift up while the key is greater than the parent's, sift down to the left child
/// unless `key[left] < key[right]`, until the key is not smaller than the child's. `pos[item]` is the heap position or -1.
struct Heap {
    h: Vec<usize>,
    k: Vec<i64>,
    pos: Vec<i64>,
}

impl Heap {
    fn new(n: usize) -> Self {
        Heap {
            h: Vec::new(),
            k: Vec::new(),
            pos: vec![-1; n],
        }
    }

    fn clear(&mut self) {
        for &it in &self.h {
            self.pos[it] = -1;
        }
        self.h.clear();
        self.k.clear();
    }

    fn top(&self) -> i64 {
        self.h.first().map_or(-1, |&v| v as i64)
    }

    fn contains(&self, it: usize) -> bool {
        self.pos[it] >= 0
    }

    fn swap(&mut self, a: usize, b: usize) {
        self.h.swap(a, b);
        self.k.swap(a, b);
        self.pos[self.h[a]] = a as i64;
        self.pos[self.h[b]] = b as i64;
    }

    fn up(&mut self, mut p: usize) {
        while p > 0 {
            let q = (p - 1) / 2;
            if self.k[p] <= self.k[q] {
                break;
            }
            self.swap(p, q);
            p = q;
        }
    }

    fn down(&mut self, mut p: usize) {
        let n = self.h.len();
        while 2 * p + 1 < n {
            let mut c = 2 * p + 1;
            if 2 * p + 2 < n && self.k[2 * p + 1] < self.k[2 * p + 2] {
                c = 2 * p + 2;
            }
            if self.k[p] >= self.k[c] {
                break;
            }
            self.swap(p, c);
            p = c;
        }
    }

    fn insert(&mut self, it: usize, key: i64) {
        self.h.push(it);
        self.k.push(key);
        self.pos[it] = self.h.len() as i64 - 1;
        self.up(self.h.len() - 1);
    }

    fn remove(&mut self, it: usize) {
        let p = self.pos[it];
        if p < 0 {
            return;
        }
        let p = p as usize;
        self.pos[it] = -1;
        let n = self.h.len();
        if p == n - 1 {
            self.h.pop();
            self.k.pop();
            return;
        }
        self.h[p] = self.h[n - 1];
        self.k[p] = self.k[n - 1];
        self.h.pop();
        self.k.pop();
        self.pos[self.h[p]] = p as i64;
        if p != 0 && self.k[p] > self.k[(p - 1) / 2] {
            self.up(p);
        } else {
            self.down(p);
        }
    }

    fn pop(&mut self) -> usize {
        let it = self.h[0];
        self.remove(it);
        it
    }

    fn update(&mut self, it: usize, key: i64) {
        let p = self.pos[it] as usize;
        self.k[p] = key;
        if p != 0 && key > self.k[(p - 1) / 2] {
            self.up(p);
        } else {
            self.down(p);
        }
    }
}

/// One matching level of `0x325850`: vertices sorted by degree (counting sort, ascending index), every run of equal degree shuffled (Fisher-Yates from the end, one
/// draw per element but the first); in that order each unmatched vertex is paired with its unmatched neighbour of smallest degree (first in adjacency order on ties)
/// or stays alone. Returns `(cptr, members, cxadj, cadj, cvw)`; coarse vertices are numbered in order of creation.
fn coarsen(
    n: usize,
    xadj: &[usize],
    adj: &[usize],
    vw: &[i64],
    rng: &mut u64,
) -> (Vec<usize>, Vec<usize>, Vec<usize>, Vec<usize>, Vec<i64>) {
    let deg: Vec<usize> = (0..n).map(|v| xadj[v + 1] - xadj[v]).collect();
    let mut cnt = vec![0usize; n + 1];
    for &d in &deg {
        cnt[d] += 1;
    }
    let mut s = 0;
    for c in cnt.iter_mut() {
        let t = *c;
        *c = s;
        s += t;
    }
    let mut order = vec![0usize; n];
    for v in 0..n {
        order[cnt[deg[v]]] = v;
        cnt[deg[v]] += 1;
    }
    let mut a = 0;
    while a < n {
        let mut b = a + 1;
        while b < n && deg[order[b]] == deg[order[a]] {
            b += 1;
        }
        let (mut last, mut size) = (b - 1, b - a);
        while last > a {
            let j = rand31(rng) as usize % size;
            order.swap(last, a + j);
            last -= 1;
            size -= 1;
        }
        a = b;
    }
    let mut matched = vec![usize::MAX; n];
    let mut cptr = vec![0];
    let mut members = Vec::new();
    for &v in &order {
        if matched[v] != usize::MAX {
            continue;
        }
        let (mut best, mut bd) = (usize::MAX, n + 1);
        let first = MATCH_SCHEME.with(|c| c.get()) == 3;
        for &u in &adj[xadj[v]..xadj[v + 1]] {
            if matched[u] == usize::MAX && (first || deg[u] < bd) {
                bd = deg[u];
                best = u;
                if first {
                    break;
                }
            }
        }
        let c = cptr.len() - 1;
        matched[v] = c;
        members.push(v);
        if best != usize::MAX {
            matched[best] = c;
            members.push(best);
        }
        cptr.push(members.len());
    }
    let cn = cptr.len() - 1;
    let mut mark = vec![usize::MAX; cn];
    let (mut cxadj, mut cadj, mut cvw) = (vec![0], Vec::new(), Vec::new());
    for c in 0..cn {
        mark[c] = c;
        let mut w = 0;
        for &v in &members[cptr[c]..cptr[c + 1]] {
            add_work(xadj[v + 1] - xadj[v], 3.0, 2.0);
            w += vw[v];
            for &u in &adj[xadj[v]..xadj[v + 1]] {
                let cu = matched[u];
                if mark[cu] != c {
                    cadj.push(cu);
                    mark[cu] = c;
                }
            }
        }
        cxadj.push(cadj.len());
        cvw.push(w);
    }
    (cptr, members, cxadj, cadj, cvw)
}

/// Separator quality `p2 / (p0 p1)`; `100 (p0 + p1 + p2)` when a side is empty.
fn separator_score(p0: i64, p1: i64, p2: i64) -> f64 {
    if p0 != 0 && p1 != 0 {
        p2 as f64 / (p0 as f64 * p1 as f64)
    } else {
        (p0 + p1 + p2) as f64 * 100.0
    }
}

/// `0x327190`: Fiduccia-Mattheyses refinement of a vertex separator (labels 0 and 1 sides, 2 separator), in place. At most 10 passes. A pass shuffles the separator
/// vertices (one draw per element), puts them in two heaps (key `w(v) - ed[other side][v]`: moving `v` to a side pulls its neighbours on the other side into the
/// separator), and moves the best vertex of the side chosen by the balance rule (the lighter side if the other weighs more than `param` times as much, else the
/// side with the smaller `p2 / (p0 p1)` after the move, ties to side 1). The best state of the pass (strictly smaller score) is kept and the later moves undone; a
/// pass stops `min(500, 5 |separator|)` moves after its best, and no pass follows one without improvement. The reference tests `contains` with the moved vertex
/// instead of the neighbour when it updates the key of a separator neighbour of the moved vertex, so that update never happens (checked against the generator state).
fn refine(
    n: usize,
    xadj: &[usize],
    adj: &[usize],
    vw: &[i64],
    lab: &mut [i64],
    param: f64,
    rng: &mut u64,
) {
    let mut pw = [0i64; 3];
    let mut ed = [vec![0i64; n], vec![0i64; n]];
    let mut nsep = 0usize;
    for v in 0..n {
        pw[lab[v] as usize] += vw[v];
        if lab[v] == 2 {
            nsep += 1;
            for &u in &adj[xadj[v]..xadj[v + 1]] {
                if lab[u] != 2 {
                    ed[lab[u] as usize][v] += vw[u];
                }
            }
        }
    }
    let ntry = (5 * nsep).min(500) as i64;
    let mut heaps = [Heap::new(n), Heap::new(n)];
    for _pass in 0..10 {
        heaps[0].clear();
        heaps[1].clear();
        let mut sep: Vec<usize> = (0..n).filter(|&v| lab[v] == 2).collect();
        for i in (0..sep.len()).rev() {
            let j = rand31(rng) as usize % (i + 1);
            sep.swap(i, j);
        }
        for &v in &sep {
            heaps[0].insert(v, vw[v] - ed[1][v]);
            heaps[1].insert(v, vw[v] - ed[0][v]);
        }
        let mut where_ = vec![-1i64; n];
        let mut best_sc = separator_score(pw[0], pw[1], pw[2]);
        let mut best = -1i64;
        let (mut moves, mut pulled, mut poff): (Vec<usize>, Vec<usize>, Vec<usize>) =
            (Vec::new(), Vec::new(), Vec::new());
        let mut k = 0usize;
        while k < n {
            if k as i64 - best > ntry {
                break;
            }
            let (c0, c1) = (heaps[0].top(), heaps[1].top());
            let side = if c0 < 0 {
                1
            } else if c1 < 0 {
                0
            } else if pw[1] as f64 > param * pw[0] as f64 {
                0
            } else if pw[0] as f64 > param * pw[1] as f64 {
                1
            } else {
                let (a, b) = (c0 as usize, c1 as usize);
                let s0 = separator_score(pw[0] + vw[a], pw[1] - ed[1][a], pw[2] - vw[a] + ed[1][a]);
                let s1 = separator_score(pw[0] - ed[0][b], pw[1] + vw[b], pw[2] - vw[b] + ed[0][b]);
                if s0 >= s1 {
                    1
                } else {
                    0
                }
            };
            if [c0, c1][side] < 0 {
                break;
            }
            let other = 1 - side;
            let v = heaps[side].pop();
            heaps[other].remove(v);
            moves.push(v);
            poff.push(pulled.len());
            lab[v] = side as i64;
            where_[v] = k as i64;
            pw[2] -= vw[v];
            pw[side] += vw[v];
            for &u in &adj[xadj[v]..xadj[v + 1]] {
                if lab[u] == 2 {
                    ed[side][u] += vw[v];
                    if heaps[other].contains(v) {
                        heaps[other].update(u, vw[u] - ed[side][u]);
                    }
                } else if lab[u] == other as i64 {
                    lab[u] = 2;
                    pw[other] -= vw[u];
                    pulled.push(u);
                    add_work(xadj[u + 1] - xadj[u], 6.0, 10.0);
                    ed[0][u] = 0;
                    ed[1][u] = 0;
                    pw[2] += vw[u];
                    for &x in &adj[xadj[u]..xadj[u + 1]] {
                        if lab[x] == 2 {
                            ed[other][x] -= vw[u];
                            if heaps[side].contains(x) {
                                heaps[side].update(x, vw[x] - ed[other][x]);
                            }
                        } else {
                            ed[lab[x] as usize][u] += vw[x];
                        }
                    }
                    if where_[u] == -1 {
                        heaps[side].insert(u, vw[u] - ed[other][u]);
                        where_[u] = -2;
                    }
                }
            }
            let sc = separator_score(pw[0], pw[1], pw[2]);
            if best_sc > sc {
                best_sc = sc;
                best = k as i64;
            }
            k += 1;
        }
        poff.push(pulled.len());
        for j in ((best + 1) as usize..k).rev() {
            let v = moves[j];
            let s = lab[v] as usize;
            lab[v] = 2;
            pw[s] -= vw[v];
            ed[0][v] = 0;
            ed[1][v] = 0;
            pw[2] += vw[v];
            for &x in &adj[xadj[v]..xadj[v + 1]] {
                if lab[x] == 2 {
                    ed[s][x] -= vw[v];
                } else {
                    ed[lab[x] as usize][v] += vw[x];
                }
            }
            for &u in &pulled[poff[j]..poff[j + 1]] {
                add_work(xadj[u + 1] - xadj[u], 3.0, 4.0);
                pw[2] -= vw[u];
                pw[1 - s] += vw[u];
                lab[u] = 1 - s as i64;
                for &x in &adj[xadj[u]..xadj[u + 1]] {
                    if lab[x] == 2 {
                        ed[1 - s][x] += vw[u];
                    }
                }
            }
        }
        if best == -1 {
            break;
        }
    }
}

/// The partition `0x325850` of one component (`size` nodes, `colptr`/`rowind` its induced graph): a multilevel vertex separator. The graph is coarsened by
/// [`coarsen`] until fewer than 100 vertices remain, the size stops shrinking or 19 levels exist; one draw gives `param = 5 u + 5`; the coarsest level starts with every
/// vertex in the separator and each level is refined ([`refine`]) and projected to the next finer one. Label 2 is the separator. `rng` is the state of the
/// generator `0x88a870` (0 at the start of a solve), advanced by every draw.
pub fn partition(size: usize, colptr: &[usize], rowind: &[usize], rng: &mut u64) -> Vec<i64> {
    let mut levels: Vec<(usize, Vec<usize>, Vec<usize>, Vec<i64>)> =
        vec![(size, colptr.to_vec(), rowind.to_vec(), vec![1; size])];
    let mut maps: Vec<(Vec<usize>, Vec<usize>)> = Vec::new();
    loop {
        let (n, cptr, members, cxadj, cadj, cvw) = {
            let (n, xa, ad, vw) = levels.last().unwrap();
            let (cptr, members, cxadj, cadj, cvw) = coarsen(*n, xa, ad, vw, rng);
            (*n, cptr, members, cxadj, cadj, cvw)
        };
        let cn = cvw.len();
        maps.push((cptr, members));
        levels.push((cn, cxadj, cadj, cvw));
        if cn < 100 || cn == n || levels.len() - 1 >= 19 {
            break;
        }
    }
    let top = levels.len() - 1;
    let param = (rand31(rng) as f64 * 4.656612873077393e-10) * 5.0 + 5.0;
    let mut lab = vec![2i64; levels[top].0];
    for l in (0..=top).rev() {
        let (n, xa, ad, vw) = &levels[l];
        refine(*n, xa, ad, vw, &mut lab, param, rng);
        if l > 0 {
            let (cptr, members) = &maps[l - 1];
            let mut fine = vec![0i64; levels[l - 1].0];
            for c in 0..*n {
                for &m in &members[cptr[c]..cptr[c + 1]] {
                    fine[m] = lab[c];
                }
            }
            lab = fine;
        }
    }
    lab
}

/// Driver `0x324330` (mode 3): the level flags handed to the core (-1 normal, -2 later, -3 dense-column nodes). A component of
/// the -1 nodes larger than `max(600, n/100)` goes through the partition routine [`partition`] (one generator, state 0, shared by the components);
/// `bisect` replaces it.
pub fn levels(p: &Pattern, n1: usize, bisect: Option<Bisect>) -> Result<Vec<i64>, StructError> {
    levels_rng(p, n1, bisect, &mut 0)
}

/// [`levels`] with the generator state `*rng_out` at the start (0 for the first run), which also returns the state after the partitions in `rng_out`.
pub fn levels_rng(
    p: &Pattern,
    n1: usize,
    bisect: Option<Bisect>,
    rng_out: &mut u64,
) -> Result<Vec<i64>, StructError> {
    levels_mode(p, n1, bisect, rng_out, false)
}

/// [`levels_rng`] in mode 1 of `0x324330` (`repeat`): the component pass is repeated until no component of -1 nodes exceeds the size limit.
pub fn levels_mode(
    p: &Pattern,
    n1: usize,
    bisect: Option<Bisect>,
    rng_out: &mut u64,
    repeat: bool,
) -> Result<Vec<i64>, StructError> {
    let n = p.n;
    let n1 = n1.min(n);
    let thr = dense_threshold(n, p.nnz());
    let mut lv: Vec<i64> = (0..n)
        .map(|i| {
            if i < n1 {
                if (p.adj(i).len() as i64) > thr {
                    -2
                } else {
                    -1
                }
            } else {
                -3
            }
        })
        .collect();
    let big = 600.max(n / 100);
    let mut bisect = bisect;
    let mut rng = *rng_out;
    let mut default = |sz: usize, cp: &[usize], ri: &[usize]| partition(sz, cp, ri, &mut rng);
    loop {
        let sep0 = lv.iter().filter(|&&x| x == -2).count();
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
        for (c, &sz) in sizes.iter().enumerate() {
            if sz > big {
                let b: Bisect = match bisect.as_mut() {
                    Some(b) => &mut **b,
                    None => &mut default,
                };
                let nodes: Vec<usize> = (0..n).filter(|&i| comp[i] == c).collect();
                let mut loc = vec![usize::MAX; n];
                for (k, &v) in nodes.iter().enumerate() {
                    loc[v] = k;
                }
                let mut cp = vec![0];
                let mut ri = Vec::new();
                for &v in &nodes {
                    ri.extend(
                        p.adj(v)
                            .iter()
                            .filter(|&&x| loc[x] != usize::MAX)
                            .map(|&x| loc[x]),
                    );
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
        if !repeat
            || !sizes.iter().any(|&z| z > big)
            || lv.iter().filter(|&&x| x == -2).count() == sep0
        {
            break;
        }
    }
    drop(default);
    *rng_out = rng;
    Ok(lv)
}

/// Shell sort of `[lo, hi)` of the paired lists (keys `k`, payload `id`): gap `size/2`, then `int(gap/2.2)` (a gap of 2 is followed by 1), insertion with a strict `>` shift.
fn shell_range(k: &mut [f64], id: &mut [i64], lo: usize, hi: usize) {
    let mut gap = (hi - lo) / 2;
    while gap > 0 {
        for i in lo + gap..hi {
            let (kx, x) = (k[i], id[i]);
            let mut j = i;
            while j >= lo + gap && k[j - gap] > kx {
                k[j] = k[j - gap];
                id[j] = id[j - gap];
                j -= gap;
            }
            k[j] = kx;
            id[j] = x;
        }
        gap = if gap == 2 {
            1
        } else {
            (gap as f64 / 2.2) as usize
        };
    }
}

/// The quicksort loop of `0x8c8b70` (`0x8bf1f0`): while the range holds 100 or more entries and fewer than 50 partitions were made on
/// this call chain, a three-way partition around the median-of-three key (first, middle, last) puts `< pivot`, `== pivot` and
/// `> pivot` in that order; the smaller side is sorted by a recursive call (depth + 1), the loop continues on the larger one; the
/// remaining range gets the shell sort.
fn quick_range(k: &mut [f64], id: &mut [i64], mut lo: usize, mut hi: usize, mut depth: u32) {
    while depth < 50 && hi - lo >= 100 {
        let (a, b, c) = (k[lo], k[lo + (hi - lo) / 2], k[hi - 1]);
        let p = if a > b {
            if c > a {
                a
            } else if a == c {
                a
            } else if b > c {
                b
            } else {
                c
            }
        } else if c > b {
            b
        } else if c > a {
            c
        } else {
            a
        };
        let (mut i, mut eq, mut j) = (lo, lo, hi);
        if p > k[lo] {
            loop {
                i += 1;
                if i >= hi {
                    eq = i;
                    break;
                }
                if p > k[i] {
                    continue;
                }
                eq = i;
                break;
            }
        }
        while i < j {
            let x = k[i];
            if p > x {
                k[i] = k[eq];
                k[eq] = x;
                id.swap(i, eq);
                i += 1;
                eq += 1;
            } else if x == p {
                i += 1;
            } else {
                j -= 1;
                let mut y = k[j];
                loop {
                    if p > y {
                        if i < j {
                            let xi = id[i];
                            k[i] = k[eq];
                            id[i] = id[eq];
                            k[eq] = k[j];
                            id[eq] = id[j];
                            k[j] = x;
                            id[j] = xi;
                            i += 1;
                            eq += 1;
                        }
                        break;
                    } else if y == p {
                        if i < j {
                            let xi = id[i];
                            k[i] = y;
                            id[i] = id[j];
                            k[j] = x;
                            id[j] = xi;
                            i += 1;
                        }
                        break;
                    } else {
                        if i >= j {
                            break;
                        }
                        j -= 1;
                        y = k[j];
                    }
                }
            }
        }
        depth += 1;
        if eq - lo <= hi - j {
            quick_range(k, id, lo, eq, depth);
            lo = j;
        } else {
            quick_range(k, id, j, hi, depth);
            hi = eq;
        }
    }
    shell_range(k, id, lo, hi);
}

/// The sorter `0x8c8b70`: ascending by `key[idx]`, not stable. The keys are copied in list order and sorted together with the list
/// (the same result as the indirect variant used below 32 entries); below 100 entries it is the shell sort alone.
fn shellsort(idx: &mut [i64], key: &[f64]) {
    let mut k: Vec<f64> = idx.iter().map(|&x| key[x as usize]).collect();
    let n = idx.len();
    quick_range(&mut k, idx, 0, n, 0);
}

/// The ordering core `0x320fa0`: an approximate-minimum-degree variant on the quotient graph. `levels` are the input flags
/// (negative; -1 first). Returns `perm[original] = elimination position`, the core's off-diagonal count and Factor Ops.
pub fn order(p: &Pattern, levels: &[i64]) -> Result<OrderResult, StructError> {
    let n = p.n;
    let ni = n as i64;
    let mut adj: Vec<Vec<i64>> = (0..n)
        .map(|i| p.adj(i).iter().map(|&x| x as i64).collect())
        .collect();
    let mut ln: Vec<i64> = adj.iter().map(|a| a.len() as i64).collect();
    let mut elen = vec![0i64; n];
    let mut deg = ln.clone();
    let mut nv = vec![1i64; n];
    let mut parent: Vec<usize> = (0..n).collect();
    let mut perm: Vec<i64> = levels.to_vec();
    let maxlev = perm
        .iter()
        .filter(|&&x| x < 0)
        .map(|&x| -x)
        .max()
        .unwrap_or(0);
    let mut head = vec![-1i64; n + 1];
    let mut prev = vec![-1i64; n];
    let mut nxt = vec![-1i64; n];
    let mut w = vec![-1i64; n];
    let mut hsh = vec![0f64; n];
    let mut nz: i64 = 0;
    let mut ops = 0f64;
    let mut work = 0f64;
    let mut nelim: i64 = 0;
    let thr = if n > 0 {
        dense_threshold(n, p.nnz())
    } else {
        100
    };
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
                adj[i] = adj[i][..l as usize]
                    .iter()
                    .copied()
                    .filter(|&j| perm[j as usize] < 0)
                    .collect();
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
                work += 5.25 * elen[i] as f64;
            }
            for &i in &lme {
                if perm[i] == -level {
                    unlink!(i);
                }
            }
            for &i in &lme {
                let mut s: i64 = 0;
                let mut h = (perm[i] * ni) as f64;
                let mut na = std::mem::take(&mut adj[i]);
                let old_len = ln[i] as usize;
                let ei = elen[i] as usize;
                let mut used = 0usize;
                for pos in 0..ei {
                    let e = na[pos];
                    if ln[e as usize] > 0 {
                        s += deg[e as usize];
                        h += e as f64;
                        na[used] = e;
                        used += 1;
                    }
                }
                let pivot_pos = used;
                let ne = (pivot_pos + 1) as i64;
                for pos in ei..old_len {
                    let j = na[pos];
                    if v as i64 != w[j as usize] && nv[j as usize] > 0 {
                        s += nv[j as usize];
                        na[used] = j;
                        used += 1;
                        h += j as f64;
                    }
                }
                // Compact both segments before inserting the pivot, so it cannot
                // overwrite an unread entry at the start of the tail.
                na.truncate(used);
                na.insert(pivot_pos, v as i64);
                ln[i] = na.len() as i64;
                work += 3.75 * ln[i] as f64;
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
            ops = (2 * k - 1) as f64 * kk / 6.0
                + ((k as f64 * s as f64 * s as f64 + ops) - s as f64 * kk);
            nelim += k;
            nv[v] = degme;
            work += 15.0 * keep.len() as f64;
            for &i in &keep {
                let i = i as usize;
                deg[i] = (ni - 1).min(degme - nv[i] + deg[i]);
            }
            let mut lst = keep;
            shellsort(&mut lst, &hsh);
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
    Ok(OrderResult {
        perm,
        nz,
        ops,
        work,
    })
}

#[cfg(test)]
fn frozen_allocating_order(p: &Pattern, levels: &[i64]) -> Result<OrderResult, StructError> {
    let n = p.n;
    let ni = n as i64;
    let mut adj: Vec<Vec<i64>> = (0..n)
        .map(|i| p.adj(i).iter().map(|&x| x as i64).collect())
        .collect();
    let mut ln: Vec<i64> = adj.iter().map(|a| a.len() as i64).collect();
    let mut elen = vec![0i64; n];
    let mut deg = ln.clone();
    let mut nv = vec![1i64; n];
    let mut parent: Vec<usize> = (0..n).collect();
    let mut perm: Vec<i64> = levels.to_vec();
    let maxlev = perm
        .iter()
        .filter(|&&x| x < 0)
        .map(|&x| -x)
        .max()
        .unwrap_or(0);
    let mut head = vec![-1i64; n + 1];
    let mut prev = vec![-1i64; n];
    let mut nxt = vec![-1i64; n];
    let mut w = vec![-1i64; n];
    let mut hsh = vec![0f64; n];
    let mut nz: i64 = 0;
    let mut ops = 0f64;
    let mut work = 0f64;
    let mut nelim: i64 = 0;
    let thr = if n > 0 {
        dense_threshold(n, p.nnz())
    } else {
        100
    };
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
                adj[i] = adj[i][..l as usize]
                    .iter()
                    .copied()
                    .filter(|&j| perm[j as usize] < 0)
                    .collect();
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
                work += 5.25 * elen[i] as f64;
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
                work += 3.75 * ln[i] as f64;
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
            ops = (2 * k - 1) as f64 * kk / 6.0
                + ((k as f64 * s as f64 * s as f64 + ops) - s as f64 * kk);
            nelim += k;
            nv[v] = degme;
            work += 15.0 * keep.len() as f64;
            for &i in &keep {
                let i = i as usize;
                deg[i] = (ni - 1).min(degme - nv[i] + deg[i]);
            }
            let mut lst = keep;
            shellsort(&mut lst, &hsh);
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
    Ok(OrderResult {
        perm,
        nz,
        ops,
        work,
    })
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
        let mut s: Vec<usize> = p
            .adj(inv[j])
            .iter()
            .map(|&i| perm1[i] as usize)
            .filter(|&x| x > j)
            .collect();
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
    let exact_ops: f64 = cnt.iter().map(|&c| (c * c) as f64).sum();
    // 2. fundamental supernodes
    let mut first = Vec::new();
    let mut width: Vec<i64> = Vec::new();
    let mut sn_of = vec![0usize; n];
    let mut j = 0;
    while j < n {
        let mut k = 1;
        while j + k < n
            && cnt[j + k] + k == cnt[j]
            && k < rows[j].len()
            && rows[j][k] == j + k
            && (n1 <= j) == (n1 <= j + k)
        {
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
                        if cnt[b] == cnt[a]
                            && (n1 <= b) == (n1 <= a)
                            && rows[a][1..] == rows[b][1..]
                        {
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
    let mut post: Vec<usize> = pre
        .iter()
        .copied()
        .filter(|&x| members[x][0] < n1)
        .collect();
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
    Symbolic {
        exact_nz: exact_nz as i64,
        exact_ops,
        perm,
        nz,
        ops,
        width: width_out,
        count,
        rowoff,
        idx,
        rng: 0,
    }
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
    /// `BarOrder = -1`: the first ordering's work over its budget (`xmm1` of `0x31dd14`) and the number of extra orderings run.
    pub ratio: f64,
    pub extra: usize,
}

/// Runs the stage with the dense columns chosen by [`dense_columns`]`(m, dense_min)`. Errors name the undecoded routine reached.
pub fn structure(
    m: &LpModel,
    dense_min: usize,
) -> Result<Structure, (Pattern, usize, Vec<usize>, StructError)> {
    structure_order(m, dense_min, false)
}

/// [`structure`] with `BarOrder = -1` (`auto`: the release profile) or `0`.
// GUROBI: the manager `0x31da80` reads BarOrder from `P+0x2bc4`; values other than -1, 0 and 1 (1 = always run the extra orderings) are not
// modelled; no witness (`plans/PLAN-barrier.md`, "Open leads").
pub fn structure_order(
    m: &LpModel,
    dense_min: usize,
    auto: bool,
) -> Result<Structure, (Pattern, usize, Vec<usize>, StructError)> {
    let dense = if dense_min == DENSE_AUTO {
        match dense_rule(m) {
            Ok(d) => d,
            Err(e) => {
                let (pattern, n1) = aat_pattern(m, &[]);
                return Err((pattern, n1, Vec::new(), e));
            }
        }
    } else {
        dense_columns(m, dense_min)
    };
    let (pattern, n1) = aat_pattern(m, &dense);
    let run = |seed: u64,
               repeat: bool,
               scheme: i32|
     -> Result<(Vec<i64>, OrderResult, Symbolic), StructError> {
        let mut rng = seed;
        take_partition_work();
        set_match_scheme(scheme);
        let lv = levels_mode(&pattern, n1, None, &mut rng, repeat);
        set_match_scheme(-1);
        let lv = lv?;
        let pw = take_partition_work();
        let mut o = order(&pattern, &lv)?;
        o.work += pw;
        let mut s = symbolic(&pattern, &o.perm, n1);
        s.rng = rng;
        Ok((lv, o, s))
    };
    let first = run(0, false, -1);
    let mut ratio = 0.0;
    let mut extra = 0usize;
    let best = match first {
        // `BarOrder = -1` (the release default), ordering manager `0x31da80` with one task per round (Threads = 1): `spent` is the work the first ordering added to
        // the manager's counter (`OrderResult::work`) and `budget` the sum of the squared exact column counts. Extra orderings (generator seeds
        // `AUTO_ORDER_SEED + t`, mode 1 of `0x324330`) run when `spent / budget < 0.1`, or when `n > 1000`, `500 nnz <= budget`, `nnz <= NZ` and `spent / budget <= 2`.
        // They run one per round, at most five, the counter growing by each one's work; the rounds stop once `(counter - start) / budget > 0.05`. The cheapest by
        // exact factor operations wins (the first on ties).
        Ok(f) if auto => {
            let budget = f.2.exact_ops;
            let nnz = pattern.nnz() as f64;
            let spent = f.1.work;
            ratio = if budget > 0.0 { spent / budget } else { 0.0 };
            let go = ratio < 0.1
                || (n_nodes(&pattern) > 1000
                    && 500.0 * nnz <= budget
                    && nnz <= f.2.exact_nz as f64
                    && ratio <= 2.0);
            if go {
                // The generator the iterations draw from (early centering) keeps the state of the first ordering: the tasks run on their own seeded states.
                let rng0 = f.2.rng;
                let mut total = spent;
                let mut best = f;
                for t in 0..5 {
                    if let Ok(c) = run(AUTO_ORDER_SEED + t, true, t as i32) {
                        extra += 1;
                        total += c.1.work;
                        if c.2.exact_ops < best.2.exact_ops {
                            best = c;
                        }
                    }
                    // GUROBI: the round loop also checks the `0x80(rsp)` budget at `0x31e47a`, not modelled (it never fires with `1e100`); no witness
                    // (`plans/PLAN-barrier.md`, "Open leads").
                    if budget <= 0.0 || total / budget > 0.05 {
                        break;
                    }
                }
                best.2.rng = rng0;
                Ok(best)
            } else {
                Ok(f)
            }
        }
        other => other,
    };
    match best {
        Ok((levels, order, sym)) => Ok(Structure {
            pattern,
            n1,
            dense,
            levels,
            order,
            sym,
            ratio,
            extra,
        }),
        Err(e) => Err((pattern, n1, dense, e)),
    }
}

fn n_nodes(p: &Pattern) -> usize {
    p.n
}

/// Generator state of the first extra ordering (`BarOrder = -1`): `Seed + task + 11` (`0x31e07d..0x31e09c`, `Seed` = `P+0x2cd0` = 0). The orderings of tasks 0, 1, 2, 4 of `fx_bar_dense3400`
/// match the reference bit for bit; task 3 does not, and none of the five do on `fx_bar_mixed12` (mode 1 of `0x324330` differs from mode 3 in more than the seed), so the
/// release profile does not use it yet.
pub const AUTO_ORDER_SEED: u64 = 11;

#[cfg(test)]
mod tests {
    use super::*;

    /// The sorter `0x8c8b70` on 1000 keys with many ties: ascending, a permutation of the input list (the exact tie order is
    /// checked against the reference captures by `compare.py barrier-struct`, `fx_bar_dense3400`).
    #[test]
    fn sorter_large_ascending_permutation() {
        let key: Vec<f64> = (0..1000).map(|i| ((i * 7919) % 37) as f64).collect();
        let mut idx: Vec<i64> = (0..1000).collect();
        shellsort(&mut idx, &key);
        assert!(idx
            .windows(2)
            .all(|w| key[w[0] as usize] <= key[w[1] as usize]));
        let mut s = idx.clone();
        s.sort();
        assert_eq!(s, (0..1000).collect::<Vec<i64>>());
    }
}

#[cfg(test)]
mod aat_unique_tests {
    use super::*;

    #[test]
    fn marked_pattern_matches_original_expansion_and_dense_nodes() {
        let base = crate::mps::parse(include_str!("../tests/fixtures/afiro.mps")).unwrap();
        let mut seed = 0xc25491abe3u64;
        let mut next = || {
            seed ^= seed << 13;
            seed ^= seed >> 7;
            seed ^= seed << 17;
            seed
        };
        for case in 0..512 {
            let nr = case % 37;
            let nc = case % 67;
            let mut model = base.clone();
            model.row_lo = vec![0.0; nr];
            model.row_up = vec![1.0; nr];
            model.obj = vec![0.0; nc];
            model.col_lo = vec![0.0; nc];
            model.col_up = vec![1.0; nc];
            model.col_start = vec![0];
            model.row_idx.clear();
            model.val.clear();
            for _ in 0..nc {
                let len = if nr == 0 {
                    0
                } else {
                    next() as usize % (nr + 3)
                };
                for _ in 0..len {
                    model.row_idx.push(next() as usize % nr);
                    model
                        .val
                        .push([0.0, -0.0, 1.0, -2.0, f64::NAN][next() as usize % 5]);
                }
                model.col_start.push(model.val.len());
            }
            let mut choices = vec![
                vec![],
                (0..nc).collect(),
                (0..nc).rev().filter(|j| j % 3 == 0).collect(),
            ];
            if nc > 0 {
                choices.push(vec![0, nc - 1, 0]);
            }
            for dense in choices {
                assert_eq!(
                    frozen_aat_pattern(&model, &dense),
                    aat_pattern(&model, &dense),
                    "case {case}, dense {dense:?}"
                );
            }
        }
    }

    #[test]
    fn marked_pattern_matches_original_small_models() {
        for source in [
            include_str!("../tests/fixtures/afiro.mps"),
            include_str!("../tests/fixtures/sc50a.mps"),
            include_str!("../tests/fixtures/blend.mps"),
            include_str!("../tests/fixtures/adlittle.mps"),
            include_str!("../tests/fixtures/bore3d.mps"),
            include_str!("../tests/fixtures/degen2.mps"),
        ] {
            let model = crate::mps::parse(source).unwrap();
            for dense in [
                vec![],
                dense_columns(&model, DENSE_AUTO),
                dense_columns(&model, 3),
            ] {
                assert_eq!(
                    frozen_aat_pattern(&model, &dense),
                    aat_pattern(&model, &dense)
                );
            }
        }
    }
}

#[cfg(test)]
mod order_row_buffer_tests {
    use super::*;

    fn check(pattern: &Pattern, levels: &[i64]) {
        let old = frozen_allocating_order(pattern, levels).unwrap();
        let new = order(pattern, levels).unwrap();
        assert_eq!(old.perm, new.perm);
        assert_eq!(old.nz, new.nz);
        assert_eq!(old.ops.to_bits(), new.ops.to_bits());
        assert_eq!(old.work.to_bits(), new.work.to_bits());
    }

    #[test]
    fn row_buffers_preserve_original_ordering_and_float_hash_effects() {
        let mut seed = 0xc1a958b62u64;
        for case in 0..512 {
            let n = case % 130;
            let mut adj = vec![Vec::new(); n];
            for i in 0..n {
                for j in i + 1..n {
                    seed ^= seed << 13;
                    seed ^= seed >> 7;
                    seed ^= seed << 17;
                    let edge = match case % 7 {
                        0 => false,
                        1 => j == i + 1,
                        2 => i % 2 != j % 2,
                        3 => true,
                        _ => seed % 17 < (case % 9) as u64,
                    };
                    if edge {
                        adj[i].push(j);
                        adj[j].push(i);
                    }
                }
            }
            let mut offsets = vec![0];
            let mut rows = Vec::new();
            for a in &mut adj {
                a.sort_unstable();
                rows.extend_from_slice(a);
                offsets.push(rows.len());
            }
            let p = Pattern {
                n,
                colptr: offsets,
                rowind: rows,
            };
            check(&p, &vec![-1; n]);
            check(&p, &(0..n).map(|i| -1 - (i % 3) as i64).collect::<Vec<_>>());
            // Large legal level labels make the floating hash additions round,
            // exercising their original sequential order without a huge graph.
            check(&p, &vec![-((1i64 << 56) / (n as i64).max(1)); n]);
        }
    }

    #[test]
    fn row_buffers_preserve_original_small_model_orderings() {
        for source in [
            include_str!("../tests/fixtures/afiro.mps"),
            include_str!("../tests/fixtures/sc50a.mps"),
            include_str!("../tests/fixtures/blend.mps"),
            include_str!("../tests/fixtures/adlittle.mps"),
            include_str!("../tests/fixtures/bore3d.mps"),
            include_str!("../tests/fixtures/degen2.mps"),
        ] {
            let model = crate::mps::parse(source).unwrap();
            for dense in [
                vec![],
                dense_columns(&model, DENSE_AUTO),
                dense_columns(&model, 3),
            ] {
                let (p, n1) = aat_pattern(&model, &dense);
                check(&p, &vec![-1; p.n]);
                check(&p, &levels(&p, n1, None).unwrap());
            }
        }
    }
}
