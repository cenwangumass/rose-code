//! Pivot order of the reference LU factorization, recovered from the Gurobi 13.0.3 primal simplex (`reverse-engineer/specs/lu-factor.md`):
//! four singleton phases followed by a bucketed Markowitz search on the remaining nucleus. Only the order is reproduced here.

use crate::real::Real;

/// The relative pivot threshold at the start of a solve (`ctx+0x390` = 2^-7).
pub const THRESHOLD: f64 = 0.0078125;
const ABS_TOL: f64 = 1e-13;

/// `cols[pos] = (is_slack, rows, values)` in incoming basis-position order. Returns `(row, pos)` per pivot step.
pub fn pivot_order(m: usize, cols: &[(bool, Vec<usize>, Vec<f64>)]) -> Vec<(usize, usize)> {
    factor_structure(m, cols).pivots
}

/// The pivot sequence together with what the sparse FTRAN model needs (`crate::ftran_order`): the singleton phase of each step
/// (1..=4, 0 for the Markowitz nucleus) and the L column of each step as row indices in chain order.
#[derive(Clone)]
pub struct Structure<R = f64> {
    pub pivots: Vec<(usize, usize)>,
    pub phase: Vec<u8>,
    pub l_rows: Vec<Vec<usize>>,
    /// L multipliers matching `l_rows` and the pivot value of each step.
    pub l_vals: Vec<Vec<R>>,
    /// The column entries before the division by the pivot, for the singleton-phase steps; NaN where not known (pair
    /// overrides, nucleus steps).
    pub l_raw: Vec<Vec<R>>,
    pub diag: Vec<R>,
    /// Pivot-row entries `(basis position, value)` of the Markowitz steps (empty for singleton steps).
    pub u_rows: Vec<Vec<(usize, R)>>,
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
    /// The structural pass emits at most one gate record; retain it when a fresh factor is copied.
    pair_gate_record: Option<PairGateRecord>,
}

pub fn factor_structure(m: usize, cols: &[(bool, Vec<usize>, Vec<f64>)]) -> Structure {
    factor_structure_r::<f64>(m, cols)
}

/// `factor_structure` with the numbers of the factor in the number type `R` (the quad mode runs the same LU in extended precision).
pub fn factor_structure_r<R: Real>(
    m: usize,
    cols: &[(bool, Vec<usize>, Vec<f64>)],
) -> Structure<R> {
    let (mut seq, mut phase, ov, loff, pair_gate_record) = phases(m, cols);
    let nsing = seq.len();
    let mut done_r = vec![false; m];
    let mut done_c = vec![false; cols.len()];
    for &(r, c) in &seq {
        done_r[r] = true;
        done_c[c] = true;
    }
    let mut nuc: Vec<(usize, Vec<(usize, R)>)> = Vec::new();
    for (pos, (_, idx, val)) in cols.iter().enumerate() {
        if done_c[pos] {
            continue;
        }
        let mut d: Vec<(usize, R)> = idx
            .iter()
            .copied()
            .zip(val.iter().copied())
            // The reference's matrix leaves out entries below the absolute tolerance (fx_e2_r25: a scaled entry 2^-44 in a row of the nucleus changed its count).
            .filter(|(i, v)| !done_r[*i] && v.abs() >= ABS_TOL)
            .map(|(i, v)| (i, R::c(v)))
            .collect();
        d.reverse();
        nuc.push((pos, d));
    }
    // The work estimate counts every input entry in active rows and columns
    // before the numeric nucleus applies its magnitude filter (0xd4e4b,
    // 0xdf521..0xdf5ed); scaled entries below ABS_TOL still contribute.
    let nuc_nnz = cols
        .iter()
        .enumerate()
        .filter(|(pos, _)| !done_c[*pos])
        .map(|(_, (slack, rows, _))| {
            if *slack {
                1
            } else {
                rows.iter().filter(|r| !done_r[**r]).count()
            }
        })
        .sum();
    let (np, nl, nv, npiv, nu, nuc_work) = nucleus(m, cols.iter().map(|c| c.1.len()).sum(), &nuc);
    seq.extend(np);
    phase.resize(seq.len(), 0);
    let mut step = vec![usize::MAX; m];
    for (s, &(r, _)) in seq.iter().enumerate() {
        step[r] = s;
    }
    let mut l_rows: Vec<Vec<usize>> = Vec::with_capacity(seq.len());
    let mut l_vals: Vec<Vec<R>> = Vec::with_capacity(seq.len());
    let mut l_raw: Vec<Vec<R>> = Vec::with_capacity(seq.len());
    let mut diag: Vec<R> = Vec::with_capacity(seq.len());
    let mut u_rows: Vec<Vec<(usize, R)>> = Vec::with_capacity(seq.len());
    for (s, &(prow, c)) in seq.iter().enumerate().take(nsing) {
        let (_, idx, val) = &cols[c];
        let pv = R::c(val[idx.iter().position(|&i| i == prow).unwrap()]);
        l_rows.push(idx.iter().copied().filter(|&i| step[i] > s).collect());
        l_vals.push(
            idx.iter()
                .zip(val)
                .filter(|(&i, _)| step[i] > s)
                .map(|(_, &x)| R::c(x) / pv)
                .collect(),
        );
        l_raw.push(
            idx.iter()
                .zip(val)
                .filter(|(&i, _)| step[i] > s)
                .map(|(_, &x)| R::c(x))
                .collect(),
        );
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
        l_raw[s] = vec![R::c(f64::NAN); v.len()];
        l_vals[s] = v;
    }
    l_rows.extend(nl);
    l_raw.extend(nv.iter().map(|v| vec![R::c(f64::NAN); v.len()]));
    l_vals.extend(nv);
    diag.extend(npiv);
    u_rows.extend(nu);
    Structure {
        pivots: seq,
        phase,
        l_rows,
        l_vals,
        l_raw,
        diag,
        u_rows,
        nuc_nnz,
        nuc_work,
        loff,
        pair_unz,
        pair_urow,
        dense_n: DENSE_N.with(|c| c.get()),
        m,
        pair_gate_record,
    }
}

impl<R: Real> Structure<R> {
    /// Repeat the records of an identical structure build without repeating its numeric work.
    pub fn repeat_dumps(&self) {
        if let Some(record) = self.pair_gate_record {
            record.emit();
        }
    }

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
            (
                0xd4ded,
                5.0 * m + 100.0 + 3.0 * (nnz_l_sing + u_sing) as f64,
            ),
            (
                0xd5300,
                10.0 * m + 3.0 * nnz_b as f64 + 7.0 * self.nuc_nnz.max(1) as f64,
            ),
        ];
        if self.nuc_work > 0.0 {
            out.push((0xd85e6, self.nuc_work));
        }
        if self.dense_n > 0 {
            let n = self.dense_n as f64;
            out.push((0xd8973, 0.15 * n * (n + 1.0)));
        }
        // The reference's U count leaves out the entries that cancelled to rounding noise (`|v| <= 1e-13`).
        let u_rows_nnz: usize = self
            .u_rows
            .iter()
            .flatten()
            .filter(|e| e.1.abs() > R::c(1e-13))
            .count();
        let nnz_lu = (nnz_l + u_sing + u_rows_nnz) as f64;
        // A size factor that ramps from 1 to 10 between 2e6 and 4e6 for 24 bytes per entry.
        let a = 24.0 * nnz_lu;
        let g = if a <= 2e6 {
            1.0
        } else if a < 4e6 {
            1.0 + (a - 2e6) * 9.0 / 2e6
        } else {
            10.0
        };
        out.push((0xdc7b2, (2.0 * m + nnz_lu) * 2.0 * g));
        if levels {
            out.push((0xdcb7b, 2.0 * (nnz_l as f64 + m)));
        }
        out
    }
}

fn phases<R: Real>(
    m: usize,
    cols: &[(bool, Vec<usize>, Vec<f64>)],
) -> (
    Vec<(usize, usize)>,
    Vec<u8>,
    crate::lu_pairs::Overrides<R>,
    usize,
    Option<PairGateRecord>,
) {
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
        let Some(pos) = (0..n).find(|&p| !col_done[p] && !cols[p].0 && cols[p].1.contains(&i))
        else {
            continue;
        };
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
    let mut ov = crate::lu_pairs::Overrides::<R>::default();
    let loff = seq.len();
    let pair_gate_record = pair_gate(m - seq.len());
    if pair_gate_record.is_some() && std::env::var_os("ROSE_NO_PAIRS").is_none() {
        let before = seq.len();
        ov = crate::lu_pairs::run(
            m, cols, &mut seq, &mut phase, &row_done3, &col_done3, &row_done, &col_done,
        );
        pair_gate_end(seq.len() - before);
    }
    (seq, phase, ov, loff, pair_gate_record)
}

struct Node<R> {
    r: usize,
    c: usize,
    v: R,
    alive: bool,
}

/// Circular lists preserve bucket order while allowing constant-time membership updates.
struct Buckets {
    heads: Vec<usize>,
    next: Vec<usize>,
    prev: Vec<usize>,
    owner: Vec<usize>,
}

impl Buckets {
    fn new(n: usize) -> Self {
        Self {
            heads: vec![usize::MAX; n],
            next: vec![usize::MAX; n],
            prev: vec![usize::MAX; n],
            owner: vec![usize::MAX; n],
        }
    }

    fn insert(&mut self, k: usize, x: usize) {
        debug_assert_eq!(self.owner[x], usize::MAX);
        let h = self.heads[k];
        if h == usize::MAX {
            self.next[x] = x;
            self.prev[x] = x;
        } else {
            let tail = self.prev[h];
            self.next[tail] = x;
            self.prev[x] = tail;
            self.next[x] = h;
            self.prev[h] = x;
        }
        self.heads[k] = x;
        self.owner[x] = k;
    }

    fn remove(&mut self, k: usize, x: usize) {
        if self.owner[x] != k {
            return;
        }
        let next = self.next[x];
        if next == x {
            self.heads[k] = usize::MAX;
        } else {
            self.next[self.prev[x]] = next;
            self.prev[next] = self.prev[x];
            if self.heads[k] == x {
                self.heads[k] = next;
            }
        }
        self.owner[x] = usize::MAX;
    }

    fn list(&self, k: usize) -> Vec<usize> {
        let h = self.heads[k];
        let mut out = Vec::new();
        if h != usize::MAX {
            let mut x = h;
            loop {
                out.push(x);
                x = self.next[x];
                if x == h {
                    break;
                }
            }
        }
        out
    }
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

thread_local! {
    /// Doubleton pass gate state kept across factorizations: `(ctx+0x4e0 attempts, ctx+0x4e4 gate word)`, and `ctx+0x2c0`, the number of
    /// factorizations including the running one.
    static PAIR_GATE: std::cell::Cell<(i64, i64)> = std::cell::Cell::new((0, 0));
    static FAC_NO: std::cell::Cell<i64> = std::cell::Cell::new(1);
    /// The relative pivot threshold `ctx+0x390` for the next factorization; the dual raises it (`set_threshold`), every solve starts at `THRESHOLD`.
    static THR: std::cell::Cell<f64> = std::cell::Cell::new(THRESHOLD);
    /// The model's columns and constraint nonzeros, which size the reference's pool (`pool_cap`).
    static POOL_DIMS: std::cell::Cell<(usize, usize)> = std::cell::Cell::new((0, 0));
}

/// Sets the model's column count and constraint nonzeros for the pool capacity; call it with `reset_pair_gate`.
pub fn set_pool_dims(n: usize, nnz_a: usize) {
    POOL_DIMS.with(|c| c.set((n, nnz_a)));
}

/// Start of a solve: the gate word starts at 0 (parameter `[env+0x2450]` neither 0 nor 1; `0x2f7e1+0x4e7`), no attempts.
pub fn reset_pair_gate() {
    PAIR_GATE.with(|c| c.set((0, 0)));
    FAC_NO.with(|c| c.set(1));
    THR.with(|c| c.set(THRESHOLD));
}

/// The relative pivot threshold (`ctx+0x390`) the next factorizations use.
pub fn set_threshold(t: f64) {
    THR.with(|c| c.set(t));
}

/// Raises the pivot threshold to `t` when it is lower.
pub fn raise_threshold(t: f64) {
    THR.with(|c| {
        if c.get() < t {
            c.set(t)
        }
    });
}

/// The doubleton gate state, saved and restored around a second structure build of the same factorization (the quad twin of the primal's mirror).
pub fn gate_get() -> (i64, i64) {
    PAIR_GATE.with(|c| c.get())
}

pub fn gate_set(g: (i64, i64)) {
    PAIR_GATE.with(|c| c.set(g));
}

/// The pivot threshold `ctx+0x390` in force.
pub fn threshold() -> f64 {
    THR.with(|c| c.get())
}

/// `ctx+0x2c0` while the next factorization runs.
pub fn set_factorization_number(n: usize) {
    FAC_NO.with(|c| c.set(n as i64));
}

const GATE_ALWAYS: i64 = 2_000_000_000;

#[derive(Clone, Copy)]
struct PairGateRecord {
    remaining: usize,
    att: i64,
    gate: i64,
    fac_no: i64,
}

impl PairGateRecord {
    fn emit(self) {
        if crate::dump::on("pairgate") {
            crate::dump::emit(
                "pairgate",
                serde_json::json!({"rem": self.remaining, "att": self.att, "gate": self.gate, "c2c0": self.fac_no}),
            );
        }
    }
}

/// The doubleton pass gate (`0xe08f0..0xe0990`): at least 4 rows remain and the gate word is `>= 2e9`, or `>= 0` with at least 50 rows remaining.
/// A run counts as an attempt.
fn pair_gate(remaining: usize) -> Option<PairGateRecord> {
    let (att, gate) = PAIR_GATE.with(|c| c.get());
    let run = remaining >= 4
        && remaining - 4 <= 0x1dcd64fc
        && gate >= 0
        && (gate >= GATE_ALWAYS || remaining >= 50);
    if run {
        let record = PairGateRecord {
            remaining,
            att,
            gate,
            fac_no: FAC_NO.with(|c| c.get()),
        };
        record.emit();
        PAIR_GATE.with(|c| c.set((att + 1, gate)));
        Some(record)
    } else {
        None
    }
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

fn nucleus<R: Real>(m: usize, nnz_b: usize, cols: &[(usize, Vec<(usize, R)>)]) -> NucleusOut<R> {
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
            if std::env::var_os("ROSE_GATE_LOG").is_some() {
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
                if std::env::var_os("ROSE_GATE_LOG").is_some() {
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
        if std::env::var_os("ROSE_GATE_LOG").is_some() {
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

/// The reference's `g` of the pool formulas: `d` raised to a power that grows with the row or column count `cnt`.
fn pool_g(d: f64, cnt: usize) -> f64 {
    if cnt < 500 {
        d * d
    } else if cnt < 3000 {
        d * d.sqrt()
    } else if cnt < 10000 {
        d * d.sqrt().sqrt()
    } else {
        d
    }
}

/// Capacity in doubles of the reference's pool (`0xe3130`, the factor object's constructor, stored at `+0x158`), from the model's rows `a`,
/// columns `b` and constraint nonzeros `nz` (parameter `[+0x50] = 5`, parameter `[+0x44]` = 0): `trunc(0.75 E) + 1`. Verified on 10 models
/// (a < 3000 and b < 3000 only; larger sizes follow `pool_g` and are untested).
fn pool_cap(a: usize, b: usize, nz: usize) -> usize {
    let (af, bf, nzf) = (a as f64, b as f64, nz as f64);
    let k = 5.0;
    let x6 = if a < b {
        ((a / 2 + 2) as f64).min(pool_g(nzf / bf, a)) * af
    } else {
        0.5 * nzf + ((b / 2) as f64).min(pool_g(nzf / af, b)) * bf
    };
    let r8 = (x6 as i64 + 1) as f64;
    let e = if nzf > 10.0 * bf {
        r8 * 6.0 + af * 30.0 * k
    } else if a < b {
        r8 * 6.0 + af * 3.0 * nzf * k / bf
    } else {
        r8 * 6.0 + af * 3.0 * (af + nzf) * k / (af + bf + 1.0)
    };
    (e * 0.75) as usize + 1
}

/// The fill of the reference's double pool at the start of the factorization (`0xd4190`, `0xd4508`): `max(cap / 5, 2 min(m/2 + 2, g) m)`
/// with `cap` the pool capacity (`pool_cap`). The nucleus nodes take 4 doubles per entry after it, so the dense work array
/// starts at an odd or even double by this parity, which fixes how its vectorised dot products peel (row strides are in doubles).
fn pool_fill(m: usize, nnz_b: usize) -> usize {
    let avg = nnz_b as f64 / m as f64;
    let g = pool_g(avg, m);
    let (n, nz) = POOL_DIMS.with(|c| c.get());
    let from_cap = (pool_cap(m, n, nz) / 5) as f64;
    from_cap.max(2.0 * ((m / 2 + 2) as f64).min(g) * m as f64) as usize
}

fn dot_lanes<R: Real, const W: usize>(mut acc: R, a: &[R], b: &[R], peel: usize) -> R {
    let len = a.len();
    if !R::LANES || len < peel + W || len < 2 {
        for t in 0..len {
            acc -= a[t] * b[t];
        }
        return acc;
    }
    let end = len - ((len - peel) % W);
    for t in 0..peel {
        acc -= a[t] * b[t];
    }
    let mut ln = [R::ZERO; W];
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
fn dense_lu<R: Real>(a: &mut [Vec<R>], base: usize) -> Vec<usize> {
    // Different column storage can change the compiler's NaN operand priority.
    // Preserve the original input for replay, including NaNs from finite overflow.
    let n = a.len();
    let original = if R::LANES {
        let mut saved = Vec::with_capacity(n * n);
        for row in a.iter() {
            saved.extend_from_slice(row);
        }
        Some(saved)
    } else {
        None
    };
    let mut perm: Vec<usize> = (0..n).collect();
    let nb = n / 4 * 4;
    let peel = |i: usize| (base + i * n) & 1;
    let pick = |a: &[Vec<R>], k: usize| {
        let mut p = k;
        for i in k + 1..n {
            if a[i][k].abs() > a[p][k].abs() {
                p = i;
            }
        }
        p
    };
    // A panel's column prefixes change only in the first row sweep. Mirror each
    // updated entry there, then reuse the completed prefixes for all remaining rows.
    let mut columns: [Vec<R>; 4] = std::array::from_fn(|_| Vec::with_capacity(n));
    for j0 in (0..nb).step_by(4) {
        for (q, col) in columns.iter_mut().enumerate() {
            col.clear();
            col.extend((0..j0).map(|t| a[t][j0 + q]));
        }
        for i in 1..j0 {
            for (q, col) in columns.iter_mut().enumerate() {
                let c = j0 + q;
                a[i][c] = dot_lanes::<R, 2>(a[i][c], &a[i][..i], &col[..i], peel(i));
                col[i] = a[i][c];
            }
        }
        for i in j0..n {
            for (q, col) in columns.iter().enumerate() {
                let c = j0 + q;
                a[i][c] = dot_lanes::<R, 2>(a[i][c], &a[i][..j0], col, peel(i));
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
        let col = &mut columns[0];
        col.clear();
        col.extend((0..k).map(|t| a[t][k]));
        for i in 1..k {
            a[i][k] = dot_lanes::<R, 8>(a[i][k], &a[i][..i], &col[..i], 0);
            col[i] = a[i][k];
        }
        for i in k..n {
            a[i][k] = dot_lanes::<R, 8>(a[i][k], &a[i][..k], col, 0);
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
    if let Some(original) = original {
        if a.iter().flatten().any(|v| v.is_nan()) {
            // Every updated entry includes its previous value, so an intermediate
            // NaN remains in the final matrix even if subsequent pivots move it.
            for (row, saved) in a.iter_mut().zip(original.chunks_exact(n)) {
                row.copy_from_slice(saved);
            }
            return dense_lu_uncached(a, base);
        }
    }
    perm
}

// Preserve the original arithmetic lowering when a factor contains a NaN.
fn dense_lu_uncached<R: Real>(a: &mut [Vec<R>], base: usize) -> Vec<usize> {
    let n = a.len();
    let mut perm: Vec<usize> = (0..n).collect();
    let nb = n / 4 * 4;
    let peel = |i: usize| (base + i * n) & 1;
    let pick = |a: &[Vec<R>], k: usize| {
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
                let col: Vec<R> = (0..i).map(|t| a[t][c]).collect();
                a[i][c] = dot_lanes::<R, 2>(a[i][c], &a[i][..i], &col, peel(i));
            }
        }
        for i in j0..n {
            for c in j0..j0 + 4 {
                let col: Vec<R> = (0..j0).map(|t| a[t][c]).collect();
                a[i][c] = dot_lanes::<R, 2>(a[i][c], &a[i][..j0], &col, peel(i));
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
            let col: Vec<R> = (0..i).map(|t| a[t][k]).collect();
            a[i][k] = dot_lanes::<R, 8>(a[i][k], &a[i][..i], &col, 0);
        }
        for i in k..n {
            let col: Vec<R> = (0..k).map(|t| a[t][k]).collect();
            a[i][k] = dot_lanes::<R, 8>(a[i][k], &a[i][..k], &col, 0);
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

#[cfg(test)]
mod bucket_tests {
    use super::Buckets;

    #[test]
    fn circular_buckets_preserve_vector_order() {
        let n = 37;
        let mut buckets = Buckets::new(n);
        let mut old = vec![Vec::new(); n];
        let mut owner: Vec<Option<usize>> = vec![None; n];
        let mut seed = 0x371829abdu64;
        for step in 0..10000 {
            seed ^= seed << 13;
            seed ^= seed >> 7;
            seed ^= seed << 17;
            let x = seed as usize % n;
            let k = (seed >> 16) as usize % n;
            if step % 3 == 0 {
                if let Some(i) = owner[x] {
                    let pos = old[i].iter().position(|&v| v == x).unwrap();
                    old[i].remove(pos);
                    buckets.remove(i, x);
                    owner[x] = None;
                    // Missing-member deletion must be harmless.
                    buckets.remove(i, x);
                }
            } else if owner[x].is_none() {
                old[k].insert(0, x);
                buckets.insert(k, x);
                owner[x] = Some(k);
            } else if !old[k].is_empty() {
                let offset = x % old[k].len();
                old[k].rotate_left(offset);
                buckets.heads[k] = old[k][0];
            }
            for (i, expected) in old.iter().enumerate() {
                assert_eq!(buckets.list(i), *expected);
            }
        }
    }
}

#[cfg(test)]
mod dense_lu_tests {
    use crate::real::Real;

    fn dot_lanes<R: Real, const W: usize>(mut acc: R, a: &[R], b: &[R], peel: usize) -> R {
        let len = a.len();
        if !R::LANES || len < peel + W || len < 2 {
            for t in 0..len {
                acc -= a[t] * b[t];
            }
            return acc;
        }
        let end = len - ((len - peel) % W);
        for t in 0..peel {
            acc -= a[t] * b[t];
        }
        let mut ln = [R::ZERO; W];
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
    fn dense_lu<R: Real>(a: &mut [Vec<R>], base: usize) -> Vec<usize> {
        let n = a.len();
        let mut perm: Vec<usize> = (0..n).collect();
        let nb = n / 4 * 4;
        let peel = |i: usize| (base + i * n) & 1;
        let pick = |a: &[Vec<R>], k: usize| {
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
                    let col: Vec<R> = (0..i).map(|t| a[t][c]).collect();
                    a[i][c] = dot_lanes::<R, 2>(a[i][c], &a[i][..i], &col, peel(i));
                }
            }
            for i in j0..n {
                for c in j0..j0 + 4 {
                    let col: Vec<R> = (0..j0).map(|t| a[t][c]).collect();
                    a[i][c] = dot_lanes::<R, 2>(a[i][c], &a[i][..j0], &col, peel(i));
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
                let col: Vec<R> = (0..i).map(|t| a[t][k]).collect();
                a[i][k] = dot_lanes::<R, 8>(a[i][k], &a[i][..i], &col, 0);
            }
            for i in k..n {
                let col: Vec<R> = (0..k).map(|t| a[t][k]).collect();
                a[i][k] = dot_lanes::<R, 8>(a[i][k], &a[i][..k], &col, 0);
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

    fn check<R: Real>(input: &[Vec<R>], base: usize, bits: impl Fn(R) -> (u64, u16)) {
        let mut expected = input.to_vec();
        let mut actual = input.to_vec();
        let old_perm = dense_lu(&mut expected, base);
        let new_perm = super::dense_lu(&mut actual, base);
        assert_eq!(
            new_perm,
            old_perm,
            "pivot permutation, n={}, base={base}",
            input.len()
        );
        for (i, (old, new)) in expected.iter().zip(&actual).enumerate() {
            for (j, (&old, &new)) in old.iter().zip(new).enumerate() {
                assert_eq!(
                    bits(new),
                    bits(old),
                    "entry ({i},{j}), n={}, base={base}",
                    input.len()
                );
            }
        }
    }

    fn cases<R: Real>(bits: impl Fn(R) -> (u64, u16) + Copy) {
        let mut seed = 0x9a3f701d284bc561u64;
        let mut next = || {
            seed ^= seed << 13;
            seed ^= seed >> 7;
            seed ^= seed << 17;
            seed
        };
        let extremes = [
            0.0,
            -0.0,
            1.0,
            -1.0,
            f64::from_bits(1),
            -f64::MIN_POSITIVE,
            1e-300,
            -1e300,
            f64::MAX,
            f64::INFINITY,
            f64::NEG_INFINITY,
            f64::from_bits(0x7ff8_0000_0000_1234),
            f64::from_bits(0xfff8_0000_0000_5678),
        ];
        for n in (0..=18).chain([23, 24, 25, 31, 32, 33, 47]) {
            for kind in 0..4 {
                let input: Vec<Vec<R>> = (0..n)
                    .map(|i| {
                        (0..n)
                            .map(|j| {
                                let r = next();
                                let x = if kind == 0 {
                                    ((r >> 11) as f64 / ((1u64 << 53) as f64) - 0.5) * 8.0
                                        + if i == j { n as f64 } else { 0.0 }
                                } else if kind == 1 {
                                    extremes[r as usize % extremes.len()]
                                } else if kind == 3 {
                                    // Finite inputs can overflow and create NaNs internally.
                                    extremes[r as usize % 9]
                                } else if i == j {
                                    1.0
                                } else if r & 1 == 0 {
                                    0.0
                                } else {
                                    -0.0
                                };
                                R::c(x)
                            })
                            .collect()
                    })
                    .collect();
                for base in 0..=1 {
                    check(&input, base, bits);
                }
            }
        }
        // First dense handoff of the parent PGO CLI on PILOT4 primal-release.
        let captured: serde_json::Value =
            serde_json::from_str(include_str!("../tests/fixtures/dense_lu_captured.json")).unwrap();
        let input: Vec<Vec<R>> = captured["matrix"]
            .as_array()
            .unwrap()
            .iter()
            .map(|row| {
                row.as_array()
                    .unwrap()
                    .iter()
                    .map(|x| R::c(x.as_f64().unwrap()))
                    .collect()
            })
            .collect();
        for base in 0..=1 {
            check(&input, base, bits);
        }
    }

    #[test]
    fn panel_columns_match_frozen_dense_lu_f64() {
        cases::<f64>(|v| (v.to_bits(), 0));
    }

    #[test]
    fn panel_columns_match_frozen_dense_lu_f80() {
        use crate::f80::F80;
        cases::<F80>(F80::repr);
        let values = [
            F80::from_repr(0x8000_0000_0000_0001, 0x3fff),
            F80::from_repr(0xffff_ffff_ffff_ffff, 0x7ffe),
            F80::from_repr(0x8000_0000_0000_0000, 1),
            F80::from_repr(1, 0),
            F80::from_repr(0x8000_0000_0000_0001, 0xbfff),
        ];
        for n in [3, 4, 5, 8, 9, 12] {
            let input: Vec<Vec<F80>> = (0..n)
                .map(|i| (0..n).map(|j| values[(i * n + j) % values.len()]).collect())
                .collect();
            for base in 0..=1 {
                check(&input, base, F80::repr);
            }
        }
    }
}
