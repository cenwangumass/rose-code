//! Pivot order of the reference LU factorization, recovered from the Gurobi 13.0.3 primal simplex (`reverse-engineer/specs/lu-factor.md`):
//! four singleton phases followed by a bucketed Markowitz search on the remaining nucleus. Only the order is reproduced here.

use std::collections::HashSet;

const THRESHOLD: f64 = 0.0078125;
const ABS_TOL: f64 = 1e-13;

/// `cols[pos] = (is_slack, rows, values)` in incoming basis-position order. Returns `(row, pos)` per pivot step.
pub fn pivot_order(m: usize, cols: &[(bool, Vec<usize>, Vec<f64>)]) -> Vec<(usize, usize)> {
    factor_structure(m, cols).pivots
}

/// The pivot sequence together with what the sparse FTRAN model needs (`crate::ftran_order`): the singleton phase of each step
/// (1..=4, 0 for the Markowitz nucleus) and the L column of each step as row indices in chain order.
pub struct Structure {
    pub pivots: Vec<(usize, usize)>,
    pub phase: Vec<u8>,
    pub l_rows: Vec<Vec<usize>>,
    /// L multipliers matching `l_rows` and the pivot value of each step.
    pub l_vals: Vec<Vec<f64>>,
    /// The column entries before the division by the pivot, for the singleton-phase steps; NaN where not known (pair
    /// overrides, nucleus steps).
    pub l_raw: Vec<Vec<f64>>,
    pub diag: Vec<f64>,
    /// Pivot-row entries `(basis position, value)` of the Markowitz steps (empty for singleton steps).
    pub u_rows: Vec<Vec<(usize, f64)>>,
    /// Entries of the nucleus (the columns left after the singleton phases, in the rows left) and the Markowitz search's
    /// work estimate.
    pub nuc_nnz: usize,
    pub nuc_work: f64,
    /// Steps done when the doubleton pass starts (`L+0xc4`): the pass's own steps count in `phase` but not here.
    pub loff: usize,
    /// U-row nonzeros of the column pairs' second pivots, by step (`crate::lu_pairs::Overrides::unz`).
    pub pair_unz: std::collections::HashMap<usize, usize>,
    /// Column order of the U rows of the column pairs' pivots (`crate::lu_pairs::Overrides::urow`).
    pub pair_urow: std::collections::HashMap<usize, Vec<usize>>,
    pub dense_n: usize,
    /// Rows of the factored matrix (the pivot count is smaller for a singular basis).
    pub m: usize,
}

pub fn factor_structure(m: usize, cols: &[(bool, Vec<usize>, Vec<f64>)]) -> Structure {
    let (mut seq, mut phase, ov, loff) = phases(m, cols);
    let nsing = seq.len();
    let done_r: HashSet<usize> = seq.iter().map(|&(r, _)| r).collect();
    let done_c: HashSet<usize> = seq.iter().map(|&(_, c)| c).collect();
    let mut nuc: Vec<(usize, Vec<(usize, f64)>)> = Vec::new();
    for (pos, (_, idx, val)) in cols.iter().enumerate() {
        if done_c.contains(&pos) {
            continue;
        }
        let mut d: Vec<(usize, f64)> = idx.iter().copied().zip(val.iter().copied()).filter(|(i, _)| !done_r.contains(i)).collect();
        d.reverse();
        nuc.push((pos, d));
    }
    let nuc_nnz = nuc.iter().map(|c| c.1.len()).sum();
    let (np, nl, nv, npiv, nu, nuc_work) = nucleus(m, cols.iter().map(|c| c.1.len()).sum(), &nuc);
    seq.extend(np);
    phase.resize(seq.len(), 0);
    let mut step = vec![usize::MAX; m];
    for (s, &(r, _)) in seq.iter().enumerate() {
        step[r] = s;
    }
    let mut l_rows: Vec<Vec<usize>> = Vec::with_capacity(seq.len());
    let mut l_vals: Vec<Vec<f64>> = Vec::with_capacity(seq.len());
    let mut l_raw: Vec<Vec<f64>> = Vec::with_capacity(seq.len());
    let mut diag: Vec<f64> = Vec::with_capacity(seq.len());
    let mut u_rows: Vec<Vec<(usize, f64)>> = Vec::with_capacity(seq.len());
    for (s, &(prow, c)) in seq.iter().enumerate().take(nsing) {
        let (_, idx, val) = &cols[c];
        let pv = val[idx.iter().position(|&i| i == prow).unwrap()];
        l_rows.push(idx.iter().copied().filter(|&i| step[i] > s).collect());
        l_vals.push(idx.iter().zip(val).filter(|(&i, _)| step[i] > s).map(|(_, &x)| x / pv).collect());
        l_raw.push(idx.iter().zip(val).filter(|(&i, _)| step[i] > s).map(|(_, &x)| x).collect());
        diag.push(pv);
        u_rows.push(Vec::new());
    }
    let pair_unz = ov.unz.clone();
    let pair_urow = ov.urow.clone();
    for (&s, &d) in &ov.diag {
        diag[s] = d;
    }
    for (s, (r, v)) in ov.l {
        l_rows[s] = r;
        l_raw[s] = vec![f64::NAN; v.len()];
        l_vals[s] = v;
    }
    l_rows.extend(nl);
    l_raw.extend(nv.iter().map(|v| vec![f64::NAN; v.len()]));
    l_vals.extend(nv);
    diag.extend(npiv);
    u_rows.extend(nu);
    Structure { pivots: seq, phase, l_rows, l_vals, l_raw, diag, u_rows, nuc_nnz, nuc_work, loff, pair_unz, pair_urow, dense_n: DENSE_N.with(|c| c.get()), m }
}

impl Structure {
    /// The factorization's increments of the reference's work counter (`0xd4190`, `crate::work`) as `(site, estimate)`, given
    /// the factored columns and whether the L solve levels are built (`crate::ftran_order`).
    pub fn work(&self, cols: &[(bool, Vec<usize>, Vec<f64>)], levels: bool) -> Vec<(u32, f64)> {
        let m = self.m as f64;
        let nnz_b: usize = cols.iter().map(|c| if c.0 { 1 } else { c.1.len() }).sum();
        let nsing = self.phase.iter().filter(|&&p| p != 0).count();
        let nnz_l: usize = self.l_rows.iter().map(Vec::len).sum();
        // The entries of the singleton steps' columns that go to U (rows pivoted before them), and with the Markowitz pivot
        // rows the off-diagonal U.
        let nnz_l_sing: usize = (0..nsing).map(|k| self.l_rows[k].len()).sum();
        // The U entries of the singleton steps' rows: basis entries in a row pivoted in a singleton step, in a column pivoted
        // later (the reference counts them directly; `nnz(B) - nsing - nnz_l_sing - nuc` differs when a pair step drops or
        // creates L entries).
        let u_sing = {
            let mut step_r = vec![usize::MAX; self.m];
            let mut step_c = vec![usize::MAX; cols.len()];
            for (k, &(r, c)) in self.pivots.iter().enumerate() {
                step_r[r] = k;
                step_c[c] = k;
            }
            let mut per = vec![0usize; nsing];
            for (pos, c) in cols.iter().enumerate() {
                for &i in &c.1 {
                    if step_r[i] < nsing && step_c[pos] > step_r[i] {
                        per[step_r[i]] += 1;
                    }
                }
            }
            // A column pair's second pivot row is updated by the first pivot's row: the reference counts the entries that
            // result (fill included, cancelled ones left out).
            for (&s, &n) in &self.pair_unz {
                per[s] = n;
            }
            per.iter().sum::<usize>()
        };
        let mut out = vec![
            (0xd4ded, 5.0 * m + 100.0 + 3.0 * (nnz_l_sing + u_sing) as f64),
            (0xd5300, 10.0 * m + 3.0 * nnz_b as f64 + 7.0 * self.nuc_nnz.max(1) as f64),
        ];
        if self.nuc_work > 0.0 {
            out.push((0xd85e6, self.nuc_work));
        }
        if self.dense_n > 0 {
            let n = self.dense_n as f64;
            out.push((0xd8973, 0.15 * n * (n + 1.0)));
        }
        // The reference's U count leaves out the entries that cancelled to rounding noise (`|v| <= 1e-13`).
        let u_rows_nnz: usize = self.u_rows.iter().flatten().filter(|e| e.1.abs() > 1e-13).count();
        let nnz_lu = (nnz_l + u_sing + u_rows_nnz) as f64;
        // A size factor that ramps from 1 to 10 between 2e6 and 4e6 for 24 bytes per entry.
        let a = 24.0 * nnz_lu;
        let g = if a <= 2e6 { 1.0 } else if a < 4e6 { 1.0 + (a - 2e6) * 9.0 / 2e6 } else { 10.0 };
        out.push((0xdc7b2, (2.0 * m + nnz_lu) * 2.0 * g));
        if levels {
            out.push((0xdcb7b, 2.0 * (nnz_l as f64 + m)));
        }
        out
    }
}

fn phases(m: usize, cols: &[(bool, Vec<usize>, Vec<f64>)]) -> (Vec<(usize, usize)>, Vec<u8>, crate::lu_pairs::Overrides, usize) {
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
    let mut phase: Vec<u8> = Vec::with_capacity(m);
    let mut ph = 1u8;
    macro_rules! pivot {
        ($i:expr, $p:expr) => {{
            row_done[$i] = true;
            col_done[$p] = true;
            seq.push(($i, $p));
            phase.push(ph);
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
    ph = 2;
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
    ph = 3;
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
    let (row_done3, col_done3) = (row_done.clone(), col_done.clone());
    // Phase 4: column singletons over the unassigned rows.
    ph = 4;
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
    let mut ov = crate::lu_pairs::Overrides::default();
    let loff = seq.len();
    if pair_gate(m - seq.len()) && std::env::var_os("ROSE_NO_PAIRS").is_none() {
        let before = seq.len();
        ov = crate::lu_pairs::run(m, cols, &mut seq, &mut phase, &row_done3, &col_done3, &row_done, &col_done);
        pair_gate_end(seq.len() - before);
    }
    (seq, phase, ov, loff)
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
type NucleusOut = (Vec<(usize, usize)>, Vec<Vec<usize>>, Vec<Vec<f64>>, Vec<f64>, Vec<Vec<(usize, f64)>>, f64);

thread_local! {
    /// Doubleton pass gate state kept across factorizations: `(ctx+0x4e0 attempts, ctx+0x4e4 gate word)`, and `ctx+0x2c0`, the number of
    /// factorizations including the running one.
    static PAIR_GATE: std::cell::Cell<(i64, i64)> = std::cell::Cell::new((0, 0));
    static FAC_NO: std::cell::Cell<i64> = std::cell::Cell::new(1);
}

/// Start of a solve: the gate word starts at 0 (parameter `[env+0x2450]` neither 0 nor 1; `0x2f7e1+0x4e7`), no attempts.
pub fn reset_pair_gate() {
    PAIR_GATE.with(|c| c.set((0, 0)));
    FAC_NO.with(|c| c.set(1));
}

/// `ctx+0x2c0` while the next factorization runs.
pub fn set_factorization_number(n: usize) {
    FAC_NO.with(|c| c.set(n as i64));
}

const GATE_ALWAYS: i64 = 2_000_000_000;

/// The doubleton pass gate (`0xe08f0..0xe0990`): at least 4 rows remain and the gate word is `>= 2e9`, or `>= 0` with at least 50 rows remaining.
/// A run counts as an attempt.
fn pair_gate(remaining: usize) -> bool {
    let (att, gate) = PAIR_GATE.with(|c| c.get());
    let run = remaining >= 4 && remaining - 4 <= 0x1dcd64fc && gate >= 0 && (gate >= GATE_ALWAYS || remaining >= 50);
    if run {
        PAIR_GATE.with(|c| c.set((att + 1, gate)));
    }
    run
}

/// End of a pass that did `steps` pivots (`0xe277e`): a productive pass (>= 20 steps) raises the gate word; from the 21st factorization on the word
/// then becomes 2e9 when `2 * word >= attempts` and -1 (pass off) otherwise.
fn pair_gate_end(steps: usize) {
    let (att, mut gate) = PAIR_GATE.with(|c| c.get());
    if gate < GATE_ALWAYS {
        if steps >= 20 {
            gate += 1;
        }
        if FAC_NO.with(|c| c.get()) > 20 {
            gate = if 2 * gate >= att { GATE_ALWAYS } else { -1 };
        }
    }
    PAIR_GATE.with(|c| c.set((att, gate)));
}

thread_local! { pub static DENSE_N: std::cell::Cell<usize> = std::cell::Cell::new(0); }
fn nucleus(m: usize, nnz_b: usize, cols: &[(usize, Vec<(usize, f64)>)]) -> NucleusOut {
    DENSE_N.with(|c| c.set(0));
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
    let mut lout: Vec<Vec<usize>> = Vec::new();
    let mut lvout: Vec<Vec<f64>> = Vec::new();
    let mut pout: Vec<f64> = Vec::new();
    let mut uout: Vec<Vec<(usize, f64)>> = Vec::new();
    let mut prows: HashSet<usize> = HashSet::new();
    let mut pcols: HashSet<usize> = HashSet::new();
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
        if sing.is_none() && !s.building {
            let n = remaining as f64;
            if std::env::var_os("ROSE_GATE_LOG").is_some() && !s.building {
                eprintln!("RGATE n={} z={}", remaining, z);
            }
            if remaining > 50 && z > 0.36 * (n * n) && n * n > 4.0 * m as f64 && m as f64 + n * n < 2e9 {
                let exact: usize = ckeys.iter().filter(|c| !pcols.contains(c)).map(|&c| s.ccnt[c]).sum();
                z = exact as f64;
                if std::env::var_os("ROSE_GATE_LOG").is_some() {
                    eprintln!("RGATE2 n={} exact={} thr={}", remaining, exact, 0.325 * n * n);
                }
                if z > 0.325 * (n * n) && !s.building {
                    // Hand the active matrix to the dense LU: rows by index, columns in ascending bucket order. (Rows in bucket order
                    // reproduce the reference's pivots but not the parity of its vectorised dot products.)
                    let mut rows: Vec<usize> = Vec::new();
                    for k in 0..s.rb.len() {
                        if k == 0 {
                            for &r in s.rst.iter().rev() {
                                if !prows.contains(&r) && s.rcnt[r] == 1 && !rows.contains(&r) {
                                    rows.push(r);
                                }
                            }
                        } else {
                            rows.extend(s.rb[k].iter().copied());
                        }
                    }
                    let mut cols_a: Vec<usize> = Vec::new();
                    for k in 0..s.cb.len() {
                        if k == 0 {
                            for &c in s.cst.iter().rev() {
                                if !pcols.contains(&c) && s.ccnt[c] == 1 && !cols_a.contains(&c) {
                                    cols_a.push(c);
                                }
                            }
                        } else {
                            cols_a.extend(s.cb[k].iter().copied());
                        }
                    }
                    rows.sort_unstable();
                    let nn = remaining;
                    if rows.len() == nn && cols_a.len() == nn {
                        let mut a = vec![vec![0.0f64; nn]; nn];
                        let rpos: std::collections::HashMap<usize, usize> = rows.iter().enumerate().map(|(i, &r)| (r, i)).collect();
                        for (cj, &c) in cols_a.iter().enumerate() {
                            for &id in &s.colchain[c] {
                                let nd = &s.nodes[id];
                                if nd.alive {
                                    a[rpos[&nd.r]][cj] = nd.v;
                                }
                            }
                        }
                        let base = std::env::var("ROSE_DBASE").ok().and_then(|v| v.parse().ok()).unwrap_or(pool_fill(m, nnz_b) & 1);
                        if std::env::var_os("ROSE_DENSE_INFO").is_some() {
                            eprintln!("DINFO n0={} n={}", n0, nn);
                        }
                        if std::env::var_os("ROSE_DENSE_DUMP").is_some() {
                            let flat: Vec<f64> = (0..nn).flat_map(|c| (0..nn).map(move |r| (r, c))).map(|(r, c)| a[r][c]).collect();
                            println!("DENSE {}", serde_json::to_string(&(nn, &rows, &cols_a, flat, base)).unwrap());
                        }
                        let perm = dense_lu(&mut a, base);
                        for k in 0..nn {
                            let piv = a[k][k];
                            out.push((rows[perm[k]], cols_a[k]));
                            let mut lr = Vec::new();
                            let mut lv = Vec::new();
                            for i in k + 1..nn {
                                if a[i][k].abs() > ABS_TOL {
                                    lr.push(rows[perm[i]]);
                                    lv.push(a[i][k]);
                                }
                            }
                            lout.push(lr);
                            lvout.push(lv);
                            pout.push(piv);
                            uout.push((k + 1..nn).filter(|&j| a[k][j] != 0.0).map(|j| (cols_a[j], a[k][j])).collect());
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
            None => match search(&mut s) {
                Some((id, c)) => {
                    search_cost = c;
                    id
                }
                None => break,
            },
        };
        let (pr, pc, pv) = (s.nodes[best].r, s.nodes[best].c, s.nodes[best].v);
        if std::env::var_os("ROSE_GATE_LOG").is_some() && !s.building {
            eprintln!("RPIV n={} r={} c={} sing={:?}", remaining, s.rcnt[pr], s.ccnt[pc], sing_kind);
        }
        if sing.is_none() {
            z += (s.rcnt[pr] as f64 - 2.0) * (s.ccnt[pc] as f64 - 2.0) - 2.0;
        }
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
    (out, lout, lvout, pout, uout, work)
}

/// Markowitz search over the count buckets; returns the node of the chosen pivot (rotating bucket heads as the reference does)
/// and the search's work factor: the candidates examined times the cost bound of the last bucket.
fn search(s: &mut Nucleus) -> Option<(usize, f64)> {
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
    let mut last_bound = 0.0;
    while k <= maxk && !done {
        let bound = if u74 < big && u69 < big { ((k * k) as f64).max((u74 * u69) as f64) } else { (k * k) as f64 };
        last_bound = bound;
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
    best.map(|b| (b, count as f64 * last_bound))
}

/// The fill of the reference's double pool at the start of the factorization (`0xd4190`): twice the row count times the smaller of
/// `m/2 + 2` and a power of the average column length. The nucleus nodes take 4 doubles per entry after it, so the dense work array
/// starts at an odd or even double by this parity, which fixes how its vectorised dot products peel (row strides are in doubles).
fn pool_fill(m: usize, nnz_b: usize) -> usize {
    let avg = nnz_b as f64 / m as f64;
    let g = if m < 500 {
        avg * avg
    } else if m < 3000 {
        avg * avg.sqrt()
    } else if m < 10000 {
        avg * avg.sqrt().sqrt()
    } else {
        avg
    };
    (2.0 * ((m / 2 + 2) as f64).min(g) * m as f64) as usize
}

fn dot_lanes<const W: usize>(mut acc: f64, a: &[f64], b: &[f64], peel: usize) -> f64 {
    let len = a.len();
    if len < peel + W || len < 2 {
        for t in 0..len {
            acc -= a[t] * b[t];
        }
        return acc;
    }
    let end = len - ((len - peel) % W);
    for t in 0..peel {
        acc -= a[t] * b[t];
    }
    let mut ln = [0.0f64; W];
    ln[0] = acc;
    let mut t = peel;
    while t < end {
        for q in 0..W {
            ln[q] -= a[t + q] * b[t + q];
        }
        t += W;
    }
    // Eight lanes are four xmm accumulators: (p0 + p1) + (p2 + p3) per lane parity, then even + odd.
    let r = if W == 2 {
        ln[0] + ln[1]
    } else {
        let ev = (ln[0] + ln[2]) + (ln[4] + ln[6]);
        let od = (ln[1] + ln[3]) + (ln[5] + ln[7]);
        ev + od
    };
    let mut acc = r;
    for t in end..len {
        acc -= a[t] * b[t];
    }
    acc
}

/// The reference's dense LU on `a[row][col]` (`0xdd1b2`): column blocks of four with partial pivoting over the rows, the dot
/// products of the block updates split in two lanes and those of the leftover columns in eight. Returns the row permutation.
fn dense_lu(a: &mut [Vec<f64>], base: usize) -> Vec<usize> {
    let n = a.len();
    let mut perm: Vec<usize> = (0..n).collect();
    let nb = n / 4 * 4;
    let peel = |i: usize| (base + i * n) & 1;
    let pick = |a: &[Vec<f64>], k: usize| {
        let mut p = k;
        for i in k + 1..n {
            if a[i][k].abs() > a[p][k].abs() {
                p = i;
            }
        }
        p
    };
    for j0 in (0..nb).step_by(4) {
        for i in 1..j0 {
            for c in j0..j0 + 4 {
                let col: Vec<f64> = (0..i).map(|t| a[t][c]).collect();
                a[i][c] = dot_lanes::<2>(a[i][c], &a[i][..i], &col, peel(i));
            }
        }
        for i in j0..n {
            for c in j0..j0 + 4 {
                let col: Vec<f64> = (0..j0).map(|t| a[t][c]).collect();
                a[i][c] = dot_lanes::<2>(a[i][c], &a[i][..j0], &col, peel(i));
            }
        }
        for k in j0..j0 + 4 {
            let p = pick(a, k);
            if p != k {
                a.swap(p, k);
                perm.swap(p, k);
            }
            let piv = a[k][k];
            for i in k + 1..n {
                let l = a[i][k] / piv;
                a[i][k] = l;
                for c in k + 1..j0 + 4 {
                    let u = a[k][c];
                    a[i][c] -= l * u;
                }
            }
        }
    }
    for k in nb..n {
        for i in 1..k {
            let col: Vec<f64> = (0..i).map(|t| a[t][k]).collect();
            a[i][k] = dot_lanes::<8>(a[i][k], &a[i][..i], &col, 0);
        }
        for i in k..n {
            let col: Vec<f64> = (0..k).map(|t| a[t][k]).collect();
            a[i][k] = dot_lanes::<8>(a[i][k], &a[i][..k], &col, 0);
        }
        let p = pick(a, k);
        if p != k {
            a.swap(p, k);
            perm.swap(p, k);
        }
        let piv = a[k][k];
        for i in k + 1..n {
            a[i][k] /= piv;
        }
    }
    perm
}
