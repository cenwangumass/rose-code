//! Frozen pre-flat-array nucleus (Vec<Vec> chains, per-pivot clones), kept as the reference of the differential tests.

use super::*;

struct Node<R> {
    r: usize,
    c: usize,
    v: R,
    alive: bool,
}


struct Nucleus<R> {
    nodes: Vec<Node<R>>,
    colchain: Vec<Vec<usize>>,
    rowchain: Vec<Vec<usize>>,
    ccnt: Vec<usize>,
    rcnt: Vec<usize>,
    cb: Buckets,
    rb: Buckets,
    /// Upper bounds on occupied count buckets; search tightens them after scanning.
    cmax: usize,
    rmax: usize,
    cst: Vec<usize>,
    rst: Vec<usize>,
}

impl<R: Real> Nucleus<R> {
    fn cins(&mut self, c: usize) {
        let k = self.ccnt[c] as isize - 1;
        if k == 0 {
            self.cst.push(c);
        } else if k > 0 {
            self.cb.insert(k as usize, c);
            self.cmax = self.cmax.max(k as usize);
        }
    }

    fn rins(&mut self, r: usize) {
        let k = self.rcnt[r] as isize - 1;
        if k == 0 {
            self.rst.push(r);
        } else if k > 0 {
            self.rb.insert(k as usize, r);
            self.rmax = self.rmax.max(k as usize);
        }
    }

    fn cdel(&mut self, c: usize) {
        if self.ccnt[c] >= 2 {
            let k = self.ccnt[c] - 1;
            self.cb.remove(k, c);
        }
    }

    fn rdel(&mut self, r: usize) {
        if self.rcnt[r] >= 2 {
            let k = self.rcnt[r] - 1;
            self.rb.remove(k, r);
        }
    }

    fn clean_col(&mut self, c: usize) {
        let nodes = &self.nodes;
        self.colchain[c].retain(|&i| nodes[i].alive);
    }

    fn clean_row(&mut self, r: usize) {
        let nodes = &self.nodes;
        self.rowchain[r].retain(|&i| nodes[i].alive);
    }

    fn add_node(&mut self, r: usize, c: usize, v: R) {
        let id = self.nodes.len();
        self.nodes.push(Node {
            r,
            c,
            v,
            alive: true,
        });
        self.rowchain[r].insert(0, id);
        self.colchain[c].insert(0, id);
        self.rcnt[r] += 1;
        self.ccnt[c] += 1;
    }
}

/// `cols`: (position, entries (row, value) in column-chain order). Rows and positions are the global indices of the basis.
type NucleusOut<R> = (
    Vec<(usize, usize)>,
    Vec<Vec<usize>>,
    Vec<Vec<R>>,
    Vec<R>,
    Vec<Vec<(usize, R)>>,
    f64,
);

pub(super) fn nucleus<R: Real>(m: usize, nnz_b: usize, cols: &[(usize, Vec<(usize, R)>)]) -> NucleusOut<R> {
    DENSE_N.with(|c| c.set(0));
    let npos = cols.iter().map(|c| c.0 + 1).max().unwrap_or(0);
    let mut s = Nucleus {
        nodes: Vec::new(),
        colchain: vec![Vec::new(); npos],
        rowchain: vec![Vec::new(); m],
        ccnt: vec![0; npos],
        rcnt: vec![0; m],
        cb: Buckets::new(m.max(npos) + 2),
        rb: Buckets::new(m.max(npos) + 2),
        cmax: 0,
        rmax: 0,
        cst: Vec::new(),
        rst: Vec::new(),
    };
    let ckeys: Vec<usize> = cols.iter().map(|c| c.0).collect();
    for (pos, ents) in cols {
        for &(i, x) in ents {
            let id = s.nodes.len();
            s.nodes.push(Node {
                r: i,
                c: *pos,
                v: x,
                alive: true,
            });
            s.colchain[*pos].push(id);
        }
    }
    for &c in &ckeys {
        let ids = s.colchain[c].clone();
        for id in ids {
            let r = s.nodes[id].r;
            s.rowchain[r].insert(0, id);
        }
    }
    let rkeys: Vec<usize> = (0..m).filter(|&r| !s.rowchain[r].is_empty()).collect();
    for &c in &ckeys {
        s.ccnt[c] = s.colchain[c].len();
    }
    for &r in &rkeys {
        s.rcnt[r] = s.rowchain[r].len();
    }
    for &c in &ckeys {
        s.cins(c);
    }
    for &r in &rkeys {
        s.rins(r);
    }
    let mut out = Vec::new();
    let mut lout: Vec<Vec<usize>> = Vec::new();
    let mut lvout: Vec<Vec<R>> = Vec::new();
    let mut pout: Vec<R> = Vec::new();
    let mut uout: Vec<Vec<(usize, R)>> = Vec::new();
    let mut prows = vec![false; m];
    let mut pcols = vec![false; npos];
    // Each temporary seen set gets a fresh generation without clearing or allocating per row/column.
    let mut seen = vec![0usize; m.max(npos)];
    let mut generation = 0usize;
    let mut remaining = ckeys.len();
    let n0 = remaining;
    // The reference's hand-off test (`0xd6736`): a nonzero estimate that starts at the nonzeros of the basis matrix and grows by
    // `(r - 2)(c - 2) - 2` per pivot, with row and column counts `r`, `c` of the pivot; past the gate it is replaced by the exact count.
    let mut z = nnz_b as f64;
    let mut work = 0.0;
    while remaining > 0 {
        let mut sing: Option<usize> = None;
        // Whether the singleton is a column singleton (true) or a row singleton.
        let mut sing_kind: Option<bool> = None;
        while let Some(c) = s.cst.pop() {
            if s.ccnt[c] == 1 {
                if let Some(&id) = s.colchain[c].iter().find(|&&i| s.nodes[i].alive) {
                    sing = Some(id);
                    sing_kind = Some(true);
                    break;
                }
            }
        }
        if sing.is_none() {
            while let Some(r) = s.rst.pop() {
                if s.rcnt[r] == 1 {
                    if let Some(&id) = s.rowchain[r].iter().find(|&&i| s.nodes[i].alive) {
                        sing = Some(id);
                        sing_kind = Some(false);
                        break;
                    }
                }
            }
        }
        if sing.is_none() {
            let n = remaining as f64;
            if gate_log() {
                eprintln!("RGATE n={} z={}", remaining, z);
            }
            if remaining > 50
                && z > 0.36 * (n * n)
                && n * n > 4.0 * m as f64
                && m as f64 + n * n < 2e9
            {
                let exact: usize = ckeys
                    .iter()
                    .filter(|c| !pcols[**c])
                    .map(|&c| s.ccnt[c])
                    .sum();
                z = exact as f64;
                if gate_log() {
                    eprintln!(
                        "RGATE2 n={} exact={} thr={}",
                        remaining,
                        exact,
                        0.325 * n * n
                    );
                }
                if z > 0.325 * (n * n) {
                    // Hand the active matrix to the dense LU: rows by index, columns in ascending bucket order. (Rows in bucket order
                    // reproduce the reference's pivots but not the parity of its vectorised dot products.)
                    let mut rows: Vec<usize> = Vec::new();
                    for k in 0..s.rb.heads.len() {
                        if k == 0 {
                            for &r in s.rst.iter().rev() {
                                if !prows[r] && s.rcnt[r] == 1 && !rows.contains(&r) {
                                    rows.push(r);
                                }
                            }
                        } else {
                            rows.extend(s.rb.list(k));
                        }
                    }
                    let mut cols_a: Vec<usize> = Vec::new();
                    for k in 0..s.cb.heads.len() {
                        if k == 0 {
                            for &c in s.cst.iter().rev() {
                                if !pcols[c] && s.ccnt[c] == 1 && !cols_a.contains(&c) {
                                    cols_a.push(c);
                                }
                            }
                        } else {
                            cols_a.extend(s.cb.list(k));
                        }
                    }
                    rows.sort_unstable();
                    let nn = remaining;
                    if rows.len() == nn && cols_a.len() == nn {
                        let mut a = vec![vec![R::ZERO; nn]; nn];
                        let rpos: std::collections::HashMap<usize, usize> =
                            rows.iter().enumerate().map(|(i, &r)| (r, i)).collect();
                        for (cj, &c) in cols_a.iter().enumerate() {
                            for &id in &s.colchain[c] {
                                let nd = &s.nodes[id];
                                if nd.alive {
                                    a[rpos[&nd.r]][cj] = nd.v;
                                }
                            }
                        }
                        let base = std::env::var("ROSE_DBASE")
                            .ok()
                            .and_then(|v| v.parse().ok())
                            .unwrap_or(pool_fill(m, nnz_b) & 1);
                        if std::env::var_os("ROSE_DENSE_INFO").is_some() {
                            eprintln!("DINFO n0={} n={}", n0, nn);
                        }
                        if std::env::var_os("ROSE_DENSE_DUMP").is_some() {
                            let flat: Vec<f64> = (0..nn)
                                .flat_map(|c| (0..nn).map(move |r| (r, c)))
                                .map(|(r, c)| a[r][c].to_f64())
                                .collect();
                            println!(
                                "DENSE {}",
                                serde_json::to_string(&(nn, &rows, &cols_a, flat, base)).unwrap()
                            );
                        }
                        let perm = dense_lu(&mut a, base);
                        for k in 0..nn {
                            let piv = a[k][k];
                            out.push((rows[perm[k]], cols_a[k]));
                            let mut lr = Vec::new();
                            let mut lv = Vec::new();
                            for i in k + 1..nn {
                                if a[i][k].abs() > R::c(ABS_TOL) {
                                    lr.push(rows[perm[i]]);
                                    lv.push(a[i][k]);
                                }
                            }
                            lout.push(lr);
                            lvout.push(lv);
                            pout.push(piv);
                            uout.push(
                                (k + 1..nn)
                                    .filter(|&j| a[k][j] != R::ZERO)
                                    .map(|j| (cols_a[j], a[k][j]))
                                    .collect(),
                            );
                        }
                        DENSE_N.with(|c| c.set(nn));
                        return (out, lout, lvout, pout, uout, work);
                    }
                }
            }
        }
        let mut search_cost = 0.0;
        let best = match sing {
            Some(id) => id,
            None => match search(&mut s, remaining) {
                Some((id, c)) => {
                    search_cost = c;
                    id
                }
                None => break,
            },
        };
        let (pr, pc, pv) = (s.nodes[best].r, s.nodes[best].c, s.nodes[best].v);
        if gate_log() {
            eprintln!(
                "RPIV n={} r={} c={} sing={:?} row={} col={}",
                remaining, s.rcnt[pr], s.ccnt[pc], sing_kind, pr, pc
            );
        }
        if sing.is_none() {
            z += (s.rcnt[pr] as f64 - 2.0) * (s.ccnt[pc] as f64 - 2.0) - 2.0;
        }
        out.push((pr, pc));
        remaining -= 1;
        prows[pr] = true;
        pcols[pc] = true;
        s.clean_row(pr);
        s.clean_col(pc);
        let s_col: usize = s.colchain[pc]
            .iter()
            .map(|&i| s.rcnt[s.nodes[i].r] - 1)
            .sum();
        let s_row: usize = s.rowchain[pr]
            .iter()
            .map(|&i| s.ccnt[s.nodes[i].c] - 1)
            .sum();
        s.rdel(pr);
        s.cdel(pc);
        let prow: Vec<usize> = s.rowchain[pr]
            .iter()
            .copied()
            .filter(|&i| s.nodes[i].c != pc)
            .collect();
        let pcol: Vec<usize> = s.colchain[pc]
            .iter()
            .copied()
            .filter(|&i| s.nodes[i].r != pr)
            .collect();
        let u: Vec<(usize, R)> = prow.iter().map(|&i| (s.nodes[i].c, s.nodes[i].v)).collect();
        let l: Vec<(usize, R)> = pcol
            .iter()
            .map(|&i| (s.nodes[i].r, s.nodes[i].v / pv))
            .collect();
        lout.push(l.iter().map(|x| x.0).collect());
        lvout.push(l.iter().map(|x| x.1).collect());
        pout.push(pv);
        uout.push(u.clone());
        let col_ids: Vec<usize> = s.colchain[pc].clone();
        let row_ids: Vec<usize> = s.rowchain[pr].clone();
        // The search's work (site `0xd85e6`): a singleton costs its remaining row (column singleton) or column (row
        // singleton); a Markowitz step its search and the cheaper of the two elimination sweeps.
        work += match sing_kind {
            Some(true) => 10.0 * (prow.len() + 1) as f64,
            Some(false) => 10.0 * (pcol.len() + 1) as f64,
            None => 6.0 * search_cost + 10.0 * (s_row - pcol.len()).min(s_col - prow.len()) as f64,
        };
        if s_col < s_row {
            for &i in &prow {
                let c = s.nodes[i].c;
                s.cdel(c);
                s.ccnt[c] -= 1;
                s.nodes[i].alive = false;
            }
            for &i in &pcol {
                let r = s.nodes[i].r;
                s.rdel(r);
            }
            for &i in &col_ids {
                s.nodes[i].alive = false;
            }
            for &i in &row_ids {
                s.nodes[i].alive = false;
            }
            for &i in &pcol {
                let r2 = s.nodes[i].r;
                s.nodes[i].alive = false;
                s.rcnt[r2] -= 1;
                let lr = l.iter().find(|x| x.0 == r2).unwrap().1;
                generation = generation.wrapping_add(1);
                if generation == 0 {
                    seen.fill(0);
                    generation = 1;
                }
                let chain = s.rowchain[r2].clone();
                for mid in chain {
                    if !s.nodes[mid].alive || s.nodes[mid].c == pc {
                        continue;
                    }
                    let c = s.nodes[mid].c;
                    if let Some(&(_, uc)) = u.iter().find(|x| x.0 == c) {
                        seen[c] = generation;
                        s.nodes[mid].v -= lr * uc;
                        if s.nodes[mid].v.abs() < R::c(ABS_TOL) {
                            s.nodes[mid].alive = false;
                            s.rcnt[r2] -= 1;
                            s.ccnt[c] -= 1;
                        }
                    }
                }
                for &(c, uv) in &u {
                    if seen[c] == generation {
                        continue;
                    }
                    let val = -lr * uv;
                    if val.abs() >= R::c(ABS_TOL) {
                        s.add_node(r2, c, val);
                    }
                }
                s.rins(r2);
            }
            for &(c, _) in &u {
                s.cins(c);
            }
        } else {
            for &i in &pcol {
                let r = s.nodes[i].r;
                s.rdel(r);
                s.rcnt[r] -= 1;
                s.nodes[i].alive = false;
            }
            for &i in &prow {
                let c = s.nodes[i].c;
                s.cdel(c);
            }
            for &i in &col_ids {
                s.nodes[i].alive = false;
            }
            for &i in &row_ids {
                s.nodes[i].alive = false;
            }
            for &(c2, uc) in &u {
                s.ccnt[c2] -= 1;
                generation = generation.wrapping_add(1);
                if generation == 0 {
                    seen.fill(0);
                    generation = 1;
                }
                let chain = s.colchain[c2].clone();
                for mid in chain {
                    if !s.nodes[mid].alive || s.nodes[mid].r == pr {
                        continue;
                    }
                    let r = s.nodes[mid].r;
                    if let Some(&(_, lv)) = l.iter().find(|x| x.0 == r) {
                        seen[r] = generation;
                        s.nodes[mid].v -= lv * uc;
                        if s.nodes[mid].v.abs() < R::c(ABS_TOL) {
                            s.nodes[mid].alive = false;
                            s.rcnt[r] -= 1;
                            s.ccnt[c2] -= 1;
                        }
                    }
                }
                for &(r, lv) in &l {
                    if seen[r] == generation {
                        continue;
                    }
                    let val = -lv * uc;
                    if val.abs() >= R::c(ABS_TOL) {
                        s.add_node(r, c2, val);
                    }
                }
                s.cins(c2);
            }
            for &(r, _) in &l {
                s.rins(r);
            }
        }
    }
    (out, lout, lvout, pout, uout, work)
}

/// Markowitz search over the count buckets; returns the node of the chosen pivot (rotating bucket heads as the reference does)
/// and the search's work factor: the candidates examined times the cost bound of the last bucket.
fn search<R: Real>(s: &mut Nucleus<R>, remaining: usize) -> Option<(usize, f64)> {
    let thr = THR.with(|c| c.get());
    let first = |b: &Buckets| b.heads.iter().position(|&v| v != usize::MAX);
    let cmin = first(&s.cb);
    let rmin = first(&s.rb);
    if cmin.is_none() && rmin.is_none() {
        return None;
    }
    let big = usize::MAX / 4;
    let (u74, u69) = (cmin.unwrap_or(big), rmin.unwrap_or(big));
    let mut k = u74.min(u69);
    // The loop bound (double `0xd6a77`, quad `0x229a3c..0x229a6d`): max(4, k / 2) examined lists, at most the number of columns left (seeds 2448, 2324; bh50093: three left).
    let limit = ((k as f64 * 0.5) as usize).max(4).min(remaining);
    let cmax = s.cb.heads[..=s.cmax].iter().rposition(|&v| v != usize::MAX);
    let rmax = s.rb.heads[..=s.rmax].iter().rposition(|&v| v != usize::MAX);
    s.cmax = cmax.unwrap_or(0);
    s.rmax = rmax.unwrap_or(0);
    let maxk = cmax.into_iter().chain(rmax).max().unwrap();
    let mut best: Option<usize> = None;
    let mut bcost = f64::INFINITY;
    let mut bratio = R::ZERO;
    let mut count = 0usize;
    let mut done = false;
    let mut last_bound = 0.0;
    while k <= maxk && !done {
        let bound = if u74 < big && u69 < big {
            ((k * k) as f64).max((u74 * u69) as f64)
        } else {
            (k * k) as f64
        };
        last_bound = bound;
        for col in [true, false] {
            let head = if col { s.cb.heads[k] } else { s.rb.heads[k] };
            if head != usize::MAX {
                let mut last: usize;
                let mut x = head;
                loop {
                    if col {
                        s.clean_col(x);
                    } else {
                        s.clean_row(x);
                    }
                    last = x;
                    let ch: Vec<usize> = if col {
                        s.colchain[x].clone()
                    } else {
                        s.rowchain[x].clone()
                    };
                    let mx = ch.iter().map(|&n| s.nodes[n].v.abs()).fold(R::ZERO, R::max);
                    if mx >= R::c(ABS_TOL) {
                        for &n in &ch {
                            let cost = (((if col {
                                s.rcnt[s.nodes[n].r]
                            } else {
                                s.ccnt[s.nodes[n].c]
                            }) - 1)
                                * k) as f64;
                            let a = s.nodes[n].v.abs();
                            if (cost < bcost && a > R::c(thr) * mx)
                                || (cost == bcost && a > bratio * mx)
                            {
                                bcost = cost;
                                bratio = a / mx;
                                best = Some(n);
                            }
                        }
                        count += 1;
                        if count >= limit {
                            best?;
                            done = true;
                            break;
                        }
                        if !(bound < bcost) {
                            done = true;
                            break;
                        }
                    }
                    x = if col { s.cb.next[x] } else { s.rb.next[x] };
                    if x == head {
                        break;
                    }
                }
                if let (true, Some(b)) = (done, best) {
                    let li = last;
                    let owner = if col { s.nodes[b].c } else { s.nodes[b].r };
                    let next = if col { s.cb.next[li] } else { s.rb.next[li] };
                    let nh = if owner == li { li } else { next };
                    if col {
                        s.cb.heads[k] = nh;
                    } else {
                        s.rb.heads[k] = nh;
                    }
                }
            }
            if done {
                break;
            }
            if best.is_some() && bound >= bcost {
                done = true;
                break;
            }
        }
        k += 1;
    }
    best.map(|b| (b, count as f64 * last_bound))
}
