//! Pivot order of the reference LU factorization, recovered from the Gurobi 13.0.3 primal simplex (`reverse-engineer/specs/lu-factor.md`):
//! four singleton phases followed by a bucketed Markowitz search on the remaining nucleus. Only the order is reproduced here.

use std::collections::HashSet;

const THRESHOLD: f64 = 0.0078125;
const ABS_TOL: f64 = 1e-13;

/// `cols[pos] = (is_slack, rows, values)` in incoming basis-position order. Returns `(row, pos)` per pivot step.
pub fn pivot_order(m: usize, cols: &[(bool, Vec<usize>, Vec<f64>)]) -> Vec<(usize, usize)> {
    let mut seq = phases(m, cols);
    let done_r: HashSet<usize> = seq.iter().map(|&(r, _)| r).collect();
    let done_c: HashSet<usize> = seq.iter().map(|&(_, c)| c).collect();
    let mut nuc: Vec<(usize, Vec<(usize, f64)>)> = Vec::new();
    for (pos, (_, idx, val)) in cols.iter().enumerate() {
        if done_c.contains(&pos) {
            continue;
        }
        let mut d: Vec<(usize, f64)> = idx.iter().copied().zip(val.iter().copied()).filter(|(i, _)| !done_r.contains(i)).collect();
        d.sort_by(|a, b| b.0.cmp(&a.0));
        nuc.push((pos, d));
    }
    seq.extend(nucleus(m, &nuc));
    seq
}

fn phases(m: usize, cols: &[(bool, Vec<usize>, Vec<f64>)]) -> Vec<(usize, usize)> {
    let n = cols.len();
    let mut rc = vec![0usize; m];
    for (is_slack, idx, _) in cols {
        if *is_slack {
            rc[idx[0]] += 1 + m;
        } else {
            for &i in idx {
                rc[i] += 1;
            }
        }
    }
    let mut row_done = vec![false; m];
    let mut col_done = vec![false; n];
    let mut seq = Vec::with_capacity(m);
    macro_rules! pivot {
        ($i:expr, $p:expr) => {{
            row_done[$i] = true;
            col_done[$p] = true;
            seq.push(($i, $p));
        }};
    }
    // Phase 1: rows with one structural column and no slack, LIFO stack.
    let mut st: Vec<usize> = (0..m).filter(|&i| rc[i] == 1).collect();
    while let Some(i) = st.pop() {
        if rc[i] == 0 || row_done[i] {
            continue;
        }
        rc[i] = 0;
        let Some(pos) = (0..n).find(|&p| !col_done[p] && !cols[p].0 && cols[p].1.contains(&i)) else { continue };
        pivot!(i, pos);
        for &r2 in &cols[pos].1 {
            if !row_done[r2] {
                rc[r2] -= 1;
                if rc[r2] == 1 {
                    st.push(r2);
                }
            }
        }
    }
    // Phase 2: slack columns whose row has only the slack left.
    for pos in 0..n {
        if col_done[pos] || !cols[pos].0 {
            continue;
        }
        let i = cols[pos].1[0];
        if rc[i] == m + 1 {
            rc[i] = 1;
            pivot!(i, pos);
        }
    }
    // Phase 3: one ascending pass.
    for pos in 0..n {
        if col_done[pos] {
            continue;
        }
        let (is_slack, idx, _) = &cols[pos];
        if !*is_slack {
            if idx.len() == 1 && rc[idx[0]] <= m && !row_done[idx[0]] {
                let i = idx[0];
                rc[i] = 1;
                pivot!(i, pos);
            }
        } else {
            let i = idx[0];
            if !row_done[i] {
                rc[i] = rc[i].wrapping_sub(m);
                pivot!(i, pos);
            }
        }
    }
    // Phase 4: column singletons over the unassigned rows.
    let mut cnt = vec![0i64; n];
    let mut rowcols: Vec<Vec<usize>> = vec![Vec::new(); m];
    for (pos, (_, idx, _)) in cols.iter().enumerate() {
        for &i in idx {
            rowcols[i].push(pos);
        }
        if !col_done[pos] {
            cnt[pos] = idx.iter().filter(|&&i| !row_done[i]).count() as i64;
        }
    }
    let mut st: Vec<usize> = (0..n).filter(|&p| !col_done[p] && cnt[p] == 1).collect();
    while let Some(pos) = st.pop() {
        if cnt[pos] <= 0 || col_done[pos] {
            continue;
        }
        let i = *cols[pos].1.iter().find(|&&i| !row_done[i]).unwrap();
        pivot!(i, pos);
        cnt[pos] = 0;
        for &c2 in &rowcols[i] {
            if !col_done[c2] && c2 != pos {
                cnt[c2] -= 1;
                if cnt[c2] == 1 {
                    st.push(c2);
                }
            }
        }
    }
    seq
}

struct Node {
    r: usize,
    c: usize,
    v: f64,
    alive: bool,
}

struct Nucleus {
    nodes: Vec<Node>,
    colchain: Vec<Vec<usize>>,
    rowchain: Vec<Vec<usize>>,
    ccnt: Vec<usize>,
    rcnt: Vec<usize>,
    cb: Vec<Vec<usize>>,
    rb: Vec<Vec<usize>>,
    cst: Vec<usize>,
    rst: Vec<usize>,
    building: bool,
}

impl Nucleus {
    fn cins(&mut self, c: usize) {
        let k = self.ccnt[c] as isize - 1;
        if k == 0 {
            self.cst.push(c);
        } else if k > 0 && !self.building {
            self.cb[k as usize].insert(0, c);
        }
    }
    fn rins(&mut self, r: usize) {
        let k = self.rcnt[r] as isize - 1;
        if k == 0 {
            self.rst.push(r);
        } else if k > 0 && !self.building {
            self.rb[k as usize].insert(0, r);
        }
    }
    fn cdel(&mut self, c: usize) {
        if self.ccnt[c] >= 2 {
            let k = self.ccnt[c] - 1;
            if let Some(p) = self.cb[k].iter().position(|&x| x == c) {
                self.cb[k].remove(p);
            }
        }
    }
    fn rdel(&mut self, r: usize) {
        if self.rcnt[r] >= 2 {
            let k = self.rcnt[r] - 1;
            if let Some(p) = self.rb[k].iter().position(|&x| x == r) {
                self.rb[k].remove(p);
            }
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
    fn add_node(&mut self, r: usize, c: usize, v: f64) {
        let id = self.nodes.len();
        self.nodes.push(Node { r, c, v, alive: true });
        self.rowchain[r].insert(0, id);
        self.colchain[c].insert(0, id);
        self.rcnt[r] += 1;
        self.ccnt[c] += 1;
    }
}

/// `cols`: (position, entries (row, value) in column-chain order). Rows and positions are the global indices of the basis.
fn nucleus(m: usize, cols: &[(usize, Vec<(usize, f64)>)]) -> Vec<(usize, usize)> {
    let npos = cols.iter().map(|c| c.0 + 1).max().unwrap_or(0);
    let mut s = Nucleus {
        nodes: Vec::new(),
        colchain: vec![Vec::new(); npos],
        rowchain: vec![Vec::new(); m],
        ccnt: vec![0; npos],
        rcnt: vec![0; m],
        cb: vec![Vec::new(); m.max(npos) + 2],
        rb: vec![Vec::new(); m.max(npos) + 2],
        cst: Vec::new(),
        rst: Vec::new(),
        building: true,
    };
    let ckeys: Vec<usize> = cols.iter().map(|c| c.0).collect();
    for (pos, ents) in cols {
        for &(i, x) in ents {
            let id = s.nodes.len();
            s.nodes.push(Node { r: i, c: *pos, v: x, alive: true });
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
    let mut prows: HashSet<usize> = HashSet::new();
    let mut pcols: HashSet<usize> = HashSet::new();
    let mut remaining = ckeys.len();
    while remaining > 0 {
        let mut sing: Option<usize> = None;
        while let Some(c) = s.cst.pop() {
            if s.ccnt[c] == 1 {
                if let Some(&id) = s.colchain[c].iter().find(|&&i| s.nodes[i].alive) {
                    sing = Some(id);
                    break;
                }
            }
        }
        if sing.is_none() {
            while let Some(r) = s.rst.pop() {
                if s.rcnt[r] == 1 {
                    if let Some(&id) = s.rowchain[r].iter().find(|&&i| s.nodes[i].alive) {
                        sing = Some(id);
                        break;
                    }
                }
            }
        }
        if sing.is_none() && s.building {
            s.building = false;
            for &c in &ckeys {
                if !pcols.contains(&c) && s.ccnt[c] >= 2 {
                    s.cins(c);
                }
            }
            for &r in &rkeys {
                if !prows.contains(&r) && s.rcnt[r] >= 2 {
                    s.rins(r);
                }
            }
            continue;
        }
        let best = match sing {
            Some(id) => id,
            None => match search(&mut s) {
                Some(id) => id,
                None => break,
            },
        };
        let (pr, pc, pv) = (s.nodes[best].r, s.nodes[best].c, s.nodes[best].v);
        out.push((pr, pc));
        remaining -= 1;
        prows.insert(pr);
        pcols.insert(pc);
        s.clean_row(pr);
        s.clean_col(pc);
        let s_col: usize = s.colchain[pc].iter().map(|&i| s.rcnt[s.nodes[i].r] - 1).sum();
        let s_row: usize = s.rowchain[pr].iter().map(|&i| s.ccnt[s.nodes[i].c] - 1).sum();
        s.rdel(pr);
        s.cdel(pc);
        let prow: Vec<usize> = s.rowchain[pr].iter().copied().filter(|&i| s.nodes[i].c != pc).collect();
        let pcol: Vec<usize> = s.colchain[pc].iter().copied().filter(|&i| s.nodes[i].r != pr).collect();
        let u: Vec<(usize, f64)> = prow.iter().map(|&i| (s.nodes[i].c, s.nodes[i].v)).collect();
        let l: Vec<(usize, f64)> = pcol.iter().map(|&i| (s.nodes[i].r, s.nodes[i].v / pv)).collect();
        let col_ids: Vec<usize> = s.colchain[pc].clone();
        let row_ids: Vec<usize> = s.rowchain[pr].clone();
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
                let mut seen: HashSet<usize> = HashSet::new();
                let chain = s.rowchain[r2].clone();
                for mid in chain {
                    if !s.nodes[mid].alive || s.nodes[mid].c == pc {
                        continue;
                    }
                    let c = s.nodes[mid].c;
                    if let Some(&(_, uc)) = u.iter().find(|x| x.0 == c) {
                        seen.insert(c);
                        s.nodes[mid].v -= lr * uc;
                        if s.nodes[mid].v.abs() < ABS_TOL {
                            s.nodes[mid].alive = false;
                            s.rcnt[r2] -= 1;
                            s.ccnt[c] -= 1;
                        }
                    }
                }
                for &(c, uv) in &u {
                    if seen.contains(&c) {
                        continue;
                    }
                    let val = -lr * uv;
                    if val.abs() >= ABS_TOL {
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
                let mut seen: HashSet<usize> = HashSet::new();
                let chain = s.colchain[c2].clone();
                for mid in chain {
                    if !s.nodes[mid].alive || s.nodes[mid].r == pr {
                        continue;
                    }
                    let r = s.nodes[mid].r;
                    if let Some(&(_, lv)) = l.iter().find(|x| x.0 == r) {
                        seen.insert(r);
                        s.nodes[mid].v -= lv * uc;
                        if s.nodes[mid].v.abs() < ABS_TOL {
                            s.nodes[mid].alive = false;
                            s.rcnt[r] -= 1;
                            s.ccnt[c2] -= 1;
                        }
                    }
                }
                for &(r, lv) in &l {
                    if seen.contains(&r) {
                        continue;
                    }
                    let val = -lv * uc;
                    if val.abs() >= ABS_TOL {
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
    out
}

/// Markowitz search over the count buckets; returns the node of the chosen pivot (rotating bucket heads as the reference does).
fn search(s: &mut Nucleus) -> Option<usize> {
    let first = |b: &Vec<Vec<usize>>| b.iter().position(|v| !v.is_empty());
    let cmin = first(&s.cb);
    let rmin = first(&s.rb);
    if cmin.is_none() && rmin.is_none() {
        return None;
    }
    let big = usize::MAX / 4;
    let (u74, u69) = (cmin.unwrap_or(big), rmin.unwrap_or(big));
    let mut k = u74.min(u69);
    let limit = ((k as f64 * 0.5) as usize).max(4);
    let maxk = s.cb.iter().rposition(|v| !v.is_empty()).into_iter().chain(s.rb.iter().rposition(|v| !v.is_empty())).max().unwrap();
    let mut best: Option<usize> = None;
    let mut bcost = f64::INFINITY;
    let mut bratio = 0.0f64;
    let mut count = 0usize;
    let mut done = false;
    while k <= maxk && !done {
        let bound = if u74 < big && u69 < big { ((k * k) as f64).max((u74 * u69) as f64) } else { (k * k) as f64 };
        for col in [true, false] {
            let lst: Vec<usize> = if col { s.cb[k].clone() } else { s.rb[k].clone() };
            if !lst.is_empty() {
                let l = lst.len();
                let mut last: Option<usize> = None;
                let mut i = 0;
                while i < l {
                    let x = lst[i];
                    if col {
                        s.clean_col(x);
                    } else {
                        s.clean_row(x);
                    }
                    last = Some(i);
                    let ch: Vec<usize> = if col { s.colchain[x].clone() } else { s.rowchain[x].clone() };
                    let mx = ch.iter().map(|&n| s.nodes[n].v.abs()).fold(0.0f64, f64::max);
                    if mx >= ABS_TOL {
                        for &n in &ch {
                            let cost = (((if col { s.rcnt[s.nodes[n].r] } else { s.ccnt[s.nodes[n].c] }) - 1) * k) as f64;
                            let a = s.nodes[n].v.abs();
                            if (cost < bcost && a > THRESHOLD * mx) || (cost == bcost && a > bratio * mx) {
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
                    i += 1;
                }
                if let (Some(li), true, Some(b)) = (last, done, best) {
                    let owner = if col { s.nodes[b].c } else { s.nodes[b].r };
                    let nh = if owner == lst[li] { li } else { (li + 1) % l };
                    let rotated: Vec<usize> = lst[nh..].iter().chain(lst[..nh].iter()).copied().collect();
                    if col {
                        s.cb[k] = rotated;
                    } else {
                        s.rb[k] = rotated;
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
    best
}
