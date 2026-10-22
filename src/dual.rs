//! Dual simplex of the reference's `dual-controlled` profile (`Method=1`, `Presolve=0`, `ScaleFlag=0`, `SimplexPricing=1`, `PerturbValue=0`).
//!
//! Rules: `reverse-engineer/specs/dual-simplex.md`. The reference keeps the reduced costs and the pivot row in a frame in which
//! the dual-feasible sign is positive: `d'_j = d_j` for a variable at its lower bound (status -1) and `-d_j` at its upper bound
//! (status -2), and the pivot row is negated for the same variables and when the leaving variable is above its bound, so every
//! ratio-test candidate has `alpha' < 0`. Basic variables are indexed by position (the row of the pivot-row BTRAN).

use crate::basis::Basis;
use crate::ftran_order::RefStruct;
use crate::lu::{Column, Factor};
use crate::model::{LpModel, Sense, INF};
use crate::simplex::{build, ref_sort, Internal, SolveError, SolveOptions, SolveResult, Status, TraceEvent};

/// Status of a nonbasic variable: at its lower bound, upper bound, free, fixed (`ctx+0x168`); a basic variable holds its position.
const LB: i32 = -1;
const UB: i32 = -2;
const FREE: i32 = -3;
const FIXED: i32 = -4;
/// Tolerances (`ctx+0x370..0x3b0`).
const BIG: f64 = 1e30;
const PRIMAL_TOL: f64 = 1e-6;
const DUAL_TOL: f64 = 1e-6;
const PIVOT_FLOOR: f64 = 3e-8;
const DROP_TOL: f64 = 1e-13;
const ALPHA_TOL: f64 = 1e-10;
const REL_PIVOT_TOL: f64 = 3e-6;
const W_EPS: f64 = 1e-8;

fn is_inf(v: f64) -> bool {
    v.abs() >= INF
}

fn minsd(a: f64, b: f64) -> f64 {
    if a < b { a } else { b }
}

fn maxsd(a: f64, b: f64) -> f64 {
    if a > b { a } else { b }
}

/// Outcome of the ratio test (`0xc74a0`).
struct Ratio {
    /// Index into the row of the entering variable (`None`: no candidate, the dual is unbounded).
    q: Option<usize>,
    /// Row entries of the bound flips, in the order of the reference's candidate list.
    flips: Vec<usize>,
}

struct Dual<'a> {
    p: &'a Internal,
    n: usize,
    m: usize,
    head: Vec<usize>,
    /// Status per variable: the basis position when basic, else `LB`, `UB`, `FREE`, `FIXED`.
    st: Vec<i32>,
    /// Basic values by position.
    xb: Vec<f64>,
    /// Reduced costs in the reference's frame.
    d: Vec<f64>,
    /// Dual steepest-edge weights by position.
    w: Vec<f64>,
    mirror: RefStruct,
    factor: Factor,
    /// Row-wise copy of the nonbasic, non-fixed structurals in the reference's storage order.
    rrows: Vec<Vec<usize>>,
    refactorizations: usize,
    /// Tolerance level (`ctx+0x2cc`): the number of nonbasic variables at an artificial upper bound.
    level: usize,
    /// Pending refactorization request (`ctx+0x294`).
    request: u32,
    /// Iterations made (`ctx+0x2a0`).
    iters: usize,
    // Pricing object P (`ctx+0x478`).
    prio: Vec<f64>,
    flags: Vec<u8>,
    list: Vec<usize>,
    thr: f64,
    value_mode: bool,
    cursor: usize,
    walk_len: usize,
    /// Slacks the row product offers, by row: those nonbasic at the last factorization and those that have left the basis since. A slack that enters the
    /// basis stays on offer until the next factorization.
    slack_nb: Vec<bool>,
    /// Rows the list walk may still see become infeasible before the list is rebuilt (`P+0x28` in list-walk mode): the candidate capacity
    /// `max(20, m^0.7)` after a rebuild, one less per touched row that is infeasible (negative: rebuild).
    pend: i64,
    cap: i64,
    /// Rebuilds of the walk list left (`P+0x38`).
    rebuilds: i64,
    first_recompute: bool,
    /// Running average of the number of rows touched per iteration.
    avg_touch: f64,
    /// Entry count of the last pivot column (`*(ctx+0x210)`).
    last_cap: usize,
    /// Updates of the factor since the last refactorization.
    updates: usize,
    /// The pivot column of the last iteration by position (for the factor's update).
    pending_alpha: Option<Vec<f64>>,
    /// Final status once the loop ends.
    status: Option<Status>,
}

/// The start basis of the dual path (setup `0xe7531`; assignments at `0xed2e8` and `0xeb69f`, `replay/gdb_dual_crash.py`). The equality rows (fixed slack)
/// are taken in ascending order of their number of structural entries (the reference's unstable sort); a row is skipped once it is marked. Two passes
/// over that order, sharing the marks and the cost array:
/// 1. (`0xeb3c5`) among the leading columns of the row (a flagged column or a coefficient below 0.1 among them cancels the row), the one with the smallest `cost / a` (ties within
///    1e-10: the column with fewer entries); it becomes basic in the row, the costs of the other columns of the row drop by `a_k c / a`, all columns of the row
///    are flagged, and every row of the chosen column is marked.
/// 2. (`0xeb5c0`) among the columns of the row that are not basic, have no cost left, at most five entries and a coefficient of at least 1e-6, the one with
///    the fewest entries (ties: the larger coefficient, then the first); every row of the chosen column is marked.
/// Position `r` of `basic` holds the column of row `r` (the slack otherwise).
pub(crate) fn dual_crash_basic(p: &Internal) -> Vec<usize> {
    let (n, m) = (p.n, p.m);
    let mut basic: Vec<usize> = (n..n + m).collect();
    let is_eq: Vec<bool> = (0..m).map(|i| p.up[n + i] - p.lo[n + i] <= 0.0).collect();
    let eq: Vec<usize> = (0..m).filter(|&i| is_eq[i]).collect();
    let mut row_cols: Vec<Vec<(usize, f64)>> = vec![Vec::new(); m];
    for j in 0..n {
        for &(r, a) in &p.cols[j] {
            row_cols[r].push((j, a));
        }
    }
    let mut keys: Vec<f64> = eq.iter().map(|&r| row_cols[r].len() as f64).collect();
    let mut order = eq.clone();
    ref_sort(&mut keys, &mut order);
    if crate::dump::on("dcr") {
        crate::dump::emit("dcr", serde_json::json!({"order": order, "counts": order.iter().map(|&r| row_cols[r].len()).collect::<Vec<_>>()}));
    }
    let mut cost: Vec<f64> = p.cost[..n].to_vec();
    let mut marked = vec![false; m];
    // Columns that start at an infinite bound (negative cost and infinite upper bound, or positive cost and infinite lower bound) are flagged from the start.
    let mut flagged: Vec<bool> = (0..n).map(|j| (p.cost[j] < 0.0 && p.up[j] >= BIG) || (p.cost[j] >= 0.0 && p.lo[j] <= -BIG)).collect();
    let mut taken = vec![false; n];
    let emit = |j: usize, r: usize| {
        if crate::dump::on("dcr") {
            crate::dump::emit("dcr", serde_json::json!({"col": j, "row": r}));
        }
    };
    for &r in &order {
        if marked[r] {
            continue;
        }
        let (mut best_ratio, mut best_count, mut best): (f64, usize, Option<(usize, f64)>) = (1e100, m, None);
        for &(j, a) in &row_cols[r] {
            if flagged[j] || a < 0.1 {
                best = None;
                break;
            }
            let ratio = cost[j] / a;
            let cnt = p.cols[j].len();
            if best_ratio - 1e-10 > ratio || (best_ratio + 1e-10 > ratio && cnt < best_count) {
                best_ratio = ratio;
                best_count = cnt;
                best = Some((j, a));
            }
        }
        if let Some((j, aj)) = best {
            emit(j, r);
            basic[r] = j;
            taken[j] = true;
            let c = cost[j];
            cost[j] = 0.0;
            for &(k, a) in &row_cols[r] {
                flagged[k] = true;
                if k != j {
                    cost[k] -= a * c / aj;
                }
            }
            for &(i, _) in &p.cols[j] {
                marked[i] = true;
            }
        }
    }
    for &r in &order {
        if marked[r] {
            continue;
        }
        let mut best: Option<usize> = None;
        for &(j, a) in &row_cols[r] {
            if taken[j] || cost[j] != 0.0 || p.up[j] - p.lo[j] <= 0.0 || p.cols[j].len() > 5 || a.abs() < 1e-6 {
                continue;
            }
            let a = a.abs();
            let aa = |c: usize| p.cols[c].iter().find(|e| e.0 == r).map_or(0.0, |e| e.1.abs());
            if best.map_or(true, |b| p.cols[j].len() < p.cols[b].len() || (p.cols[j].len() == p.cols[b].len() && a > aa(b))) {
                best = Some(j);
            }
        }
        if let Some(j) = best {
            emit(j, r);
            basic[r] = j;
            taken[j] = true;
            for &(i, _) in &p.cols[j] {
                marked[i] = true;
            }
        }
    }
    basic
}

/// The reference's pivot order of the basis columns (`lu_order`), as in the primal path.
fn reference_order(p: &Internal, basic: &[usize], fac_no: usize) -> (Vec<usize>, crate::lu_order::Structure, Vec<(bool, Vec<usize>, Vec<f64>)>) {
    let n = p.n;
    let cols: Vec<(bool, Vec<usize>, Vec<f64>)> = basic
        .iter()
        .map(|&v| (v >= n, p.fcols[v].iter().map(|&(r, _)| r).collect(), p.fcols[v].iter().map(|&(_, a)| a).collect()))
        .collect();
    crate::lu_order::set_factorization_number(fac_no);
    let st = crate::lu_order::factor_structure(p.m, &cols);
    let order = st.pivots.iter().map(|&(_, pos)| basic[pos]).collect();
    (order, st, cols)
}

impl<'a> Dual<'a> {
    fn nonbasic_value(&self, j: usize) -> f64 {
        match self.st[j] {
            LB | FIXED => self.p.lo[j],
            UB => self.p.up[j],
            _ => 0.0,
        }
    }

    /// Primal feasibility tolerance in force (`0x380`, or `0x370 * (1e-5 * 0x380)` while the level is not zero).
    fn ptol(&self) -> f64 {
        if self.level != 0 { BIG * (1e-5 * PRIMAL_TOL) } else { PRIMAL_TOL }
    }

    fn infeas(&self, r: usize) -> f64 {
        let v = self.head[r];
        let x = self.xb[r];
        let (lo, up) = (self.p.lo[v], self.p.up[v]);
        if x < lo {
            lo - x
        } else if up < BIG && x > up {
            x - up
        } else {
            0.0
        }
    }

    // ------------------------------------------------------------------ factorization

    /// Factorize the basis (`0xd4190`): the positions follow the pivot order of the reference's LU; the weights move with their variables.
    fn refactor(&mut self) -> Result<(), SolveError> {
        let (order, st, rows) = reference_order(self.p, &self.head, self.refactorizations + 1);
        if crate::dump::on("luin") {
            crate::dump::emit("luin", serde_json::json!({"n": self.refactorizations, "cols": self.head, "out": order, "prow": st.pivots.iter().map(|p| p.0).collect::<Vec<_>>(), "phase": st.phase}));
        }
        let mut wv = vec![0.0; self.n + self.m];
        for (pos, &v) in self.head.iter().enumerate() {
            wv[v] = self.w[pos];
        }
        for (pos, &v) in order.iter().enumerate() {
            self.st[v] = pos as i32;
            self.w[pos] = wv[v];
        }
        self.head = order;
        self.mirror = RefStruct::fresh(self.m, self.n, &st, &rows);
        let cols: Vec<&Column> = self.head.iter().map(|&j| &self.p.cols[j]).collect();
        self.factor = Factor::new(self.m, &cols).map_err(|_| SolveError::SingularBasis)?;
        self.refactorizations += 1;
        self.updates = 0;
        self.slack_nb = (0..self.m).map(|i| self.st[self.n + i] < 0).collect();
        Ok(())
    }

    /// The row-wise copy of the nonbasic, non-fixed structurals (`ctx+0x1d0`), built once at the start, column by column in index order; the pivots update
    /// it in place and no factorization rebuilds it.
    fn build_rrows(&mut self) {
        for r in self.rrows.iter_mut() {
            r.clear();
        }
        for j in 0..self.n {
            if self.st[j] < 0 && self.p.lo[j] != self.p.up[j] {
                for &(i, _) in &self.p.cols[j] {
                    self.rrows[i].push(j);
                }
            }
        }
    }

    // ------------------------------------------------------------------ basic solution

    /// x_B from scratch (`0xc0c29` path of `0xc0050`): FTRAN of `b - N x_N`.
    fn recompute_xb(&mut self) {
        let p = self.p;
        let mut r = p.rhs.clone();
        for j in 0..p.n + p.m {
            if self.st[j] < 0 {
                let v = self.nonbasic_value(j);
                if v != 0.0 {
                    for &(i, a) in &p.cols[j] {
                        r[i] -= a * v;
                    }
                }
            }
        }
        let rows: Vec<usize> = (0..p.m).filter(|&i| r[i] != 0.0).collect();
        let vals: Vec<f64> = rows.iter().map(|&i| r[i]).collect();
        let out = self.mirror.ftran_basic(&rows, &vals);
        let mut xb = vec![0.0; p.m];
        for (i, v) in out {
            xb[i] = v;
        }
        let mut res = r.clone();
        for (pos, &j) in self.head.iter().enumerate() {
            for &(i, a) in &p.cols[j] {
                res[i] -= a * xb[pos];
            }
        }
        let scale = 1.0 + r.iter().fold(0.0f64, |a, v| a.max(v.abs()));
        if res.iter().all(|v| v.abs() <= 1e-7 * scale) {
            self.xb = xb;
        } else {
            self.factor.ftran(&mut r);
            self.xb = r;
        }
    }

    /// d recomputed from the costs after a factorization, in the reference frame.
    fn recompute_d(&mut self) {
        let dt = self.true_d();
        for j in 0..self.n + self.m {
            if self.st[j] < 0 {
                self.d[j] = if self.st[j] == UB { -dt[j] } else { dt[j] };
            } else {
                self.d[j] = 0.0;
            }
        }
    }

    /// d from the costs (`d_j = c_j - y a_j`, in the true frame) with y = B^-T c_B from the reference BTRAN; the terms are subtracted in the order of the
    /// rows in the BTRAN output list (the clean-up pass of the primal path does the same, `Solver::set_d_from_y`).
    fn true_d(&mut self) -> Vec<f64> {
        let p = self.p;
        let inp: Vec<(usize, f64)> = self.head.iter().enumerate().map(|(pos, &j)| (pos, p.cost[j])).filter(|e| e.1 != 0.0).collect();
        let out = self.mirror.btran(&inp, false, true);
        let mut y = vec![0.0; p.m];
        let mut rank = vec![usize::MAX; p.m];
        for (k, &(i, v)) in out.iter().enumerate() {
            y[i] = v;
            rank[i] = k;
        }
        let mut d = vec![0.0; p.n + p.m];
        for j in 0..p.n + p.m {
            if self.st[j] < 0 {
                let mut e: Vec<(usize, f64)> = p.fcols[j].clone();
                e.sort_by_key(|x| rank[x.0]);
                let mut acc = p.cost[j];
                for &(i, a) in &e {
                    acc -= y[i] * a;
                }
                d[j] = acc;
            }
        }
        d
    }

    // ------------------------------------------------------------------ pricing object

    /// `0xbdc10`: flags of all rows and the list of the infeasible ones sorted by `infeas^2 / W` (descending).
    fn build_walk_list(&mut self) {
        let ptol = self.ptol();
        let mut rows = Vec::new();
        let mut keys = Vec::new();
        for r in 0..self.m {
            let f = self.infeas(r);
            if f > ptol {
                self.flags[r] = 0;
                rows.push(r);
                keys.push(-(f * f / self.w[r]));
            } else {
                self.flags[r] = 2;
            }
        }
        ref_sort(&mut keys, &mut rows);
        self.walk_len = rows.len();
        self.list = rows;
        self.cursor = 0;
        self.pend = self.cap;
    }

    /// The rebuild of the walk list inside the update (`0xc2291..0xc2900`): the rows that became infeasible outside the list (flag 1), sorted by priority,
    /// are merged into the entries still waiting (flag 0, kept in order); ties go to the waiting entries.
    fn merge_walk_list(&mut self) {
        let mut fresh: Vec<usize> = (0..self.m).filter(|&r| self.flags[r] == 1).collect();
        let key = |s: &Self, r: usize| {
            let f = s.infeas(r);
            -(f * f / s.w[r])
        };
        let mut keys: Vec<f64> = fresh.iter().map(|&r| key(self, r)).collect();
        ref_sort(&mut keys, &mut fresh);
        let old: Vec<usize> = self.list[..self.walk_len.min(self.list.len())].iter().copied().filter(|&r| self.flags[r] == 0).collect();
        let mut merged = Vec::with_capacity(old.len() + fresh.len());
        let (mut a, mut b) = (0, 0);
        while a < old.len() && b < fresh.len() {
            if key(self, old[a]) <= keys[b] {
                merged.push(old[a]);
                a += 1;
            } else {
                merged.push(fresh[b]);
                b += 1;
            }
        }
        merged.extend_from_slice(&old[a..]);
        merged.extend_from_slice(&fresh[b..]);
        for &r in &fresh {
            self.flags[r] = 0;
        }
        self.walk_len = merged.len();
        self.list = merged;
        self.cursor = 0;
        self.pend = self.cap;
    }

    /// The recompute path of `0xc0050` after a factorization: priorities, and (after the first one) the value mode with its threshold and list. An
    /// infeasible row below the threshold gets flag 0 and stays out of the list.
    fn recompute_flags(&mut self, big: bool, after_pivot: bool) {
        let ptol = self.ptol();
        let mut maxp = 0.0f64;
        let mut inf = vec![false; self.m];
        for r in 0..self.m {
            let f = self.infeas(r);
            if f > ptol {
                inf[r] = true;
                self.prio[r] = f * f / self.w[r];
                maxp = maxsd(maxp, self.prio[r]);
            }
        }
        if self.first_recompute {
            self.first_recompute = false;
            self.build_walk_list();
            return;
        }
        self.value_mode = true;
        self.cursor = 0;
        // New threshold (`0xc1e92..0xc1ee9`): max priority / 1e5, or / 1e3 when the level was zero before this factorization (no huge infeasibilities)
        // and the list is not shorter than the pivot column's entry count (always, after a no-row call); an existing threshold below the maximum
        // priority is only ever raised.
        let thr_old = self.thr;
        let new5 = maxp / 1e5;
        let thr5 = if thr_old != 0.0 && maxp > thr_old { maxsd(new5, thr_old) } else { new5 };
        let count = (0..self.m).filter(|&r| inf[r] && self.prio[r] >= thr5).count();
        let div = if after_pivot && count < self.last_cap || big { 1e5 } else { 1e3 };
        let new = maxp / div;
        self.thr = if thr_old != 0.0 && maxp > thr_old { maxsd(new, thr_old) } else { new };
        self.list.clear();
        for r in 0..self.m {
            if !inf[r] {
                self.flags[r] = 2;
            } else if self.prio[r] >= self.thr {
                self.flags[r] = 1;
                self.list.push(r);
            } else {
                self.flags[r] = 0;
            }
        }
    }

    /// The leaving row (`0xccbc0`), `None` when no row is infeasible.
    fn price(&mut self) -> Option<usize> {
        if !self.value_mode {
            loop {
                while self.cursor < self.walk_len {
                    let r = self.list[self.cursor];
                    self.cursor += 1;
                    if self.flags[r] == 0 {
                        return Some(r);
                    }
                }
                // The list is exhausted: rebuild it (`0xbdc10`), the first entry is the argmax.
                self.rebuilds -= 1;
                self.build_walk_list();
                if self.walk_len == 0 {
                    return None;
                }
                let r = self.list[0];
                self.cursor = 1;
                return Some(r);
            }
        }
        // Value mode (`0xccd77`): the list is compacted; a row that is feasible (flag 2..4) leaves it with flag 2, one below the threshold leaves it with
        // flag 0; the leaving row is the first maximum of the priority among the rest.
        let mut best = 0.5 * self.thr;
        let mut row = None;
        let mut kept = Vec::with_capacity(self.list.len());
        for &e in &self.list {
            if self.flags[e] >= 2 {
                if self.flags[e] < 5 {
                    self.flags[e] = 2;
                }
                continue;
            }
            if self.thr > self.prio[e] {
                self.flags[e] = 0;
                continue;
            }
            kept.push(e);
            if self.prio[e] > best {
                best = self.prio[e];
                row = Some(e);
            }
        }
        self.list = kept;
        if row.is_none() {
            // Refresh (`0xcd031`): argmax over the rows with flag 0, threshold divided by 10, those at or above it join the list.
            self.thr /= 10.0;
            let mut b = 0.0;
            for r in 0..self.m {
                if self.flags[r] == 0 {
                    let p = self.prio[r];
                    if p > b {
                        b = p;
                        row = Some(r);
                    }
                    if p >= self.thr {
                        self.list.push(r);
                        self.flags[r] = 1;
                    }
                }
            }
            if self.list.is_empty() {
                if let Some(r) = row {
                    self.list.push(r);
                    self.flags[r] = 1;
                    self.thr = b;
                } else {
                    self.thr = 1e-13;
                }
            }
        }
        row
    }

    /// Rows whose x_B changed in a column or a bound flip (`0xc0050` loops). List-walk mode: infeasible rows get flag 1 and use up the pending count,
    /// feasible ones flag 2. Value mode: a feasible row gets bit 1 of its flag set (flag 1 becomes 3, to leave the list at the next pricing call); an
    /// infeasible one has its priority recomputed and, unless it is in the list (flag 1), joins it when the priority reaches the threshold (flag 1) or gets
    /// flag 0 otherwise; flag 3 becomes 1.
    fn touch_rows(&mut self, rows: &[usize]) {
        let ptol = self.ptol();
        for &r in rows {
            let f = self.infeas(r);
            if !self.value_mode {
                if f > ptol {
                    self.flags[r] = 1;
                    self.pend -= 1;
                } else {
                    self.flags[r] = 2;
                }
            } else if f <= ptol {
                self.flags[r] |= 2;
            } else {
                self.prio[r] = f * f / self.w[r];
                match self.flags[r] {
                    3 => self.flags[r] = 1,
                    1 => {}
                    _ => {
                        if self.prio[r] >= self.thr {
                            self.list.push(r);
                            self.flags[r] = 1;
                        } else {
                            self.flags[r] = 0;
                        }
                    }
                }
            }
        }
    }

    // ------------------------------------------------------------------ iteration

    /// The row product (`0x134d90`): the pivot row `rho^T A_N` over the nonbasic columns (and the leaving slack), in the reference frame.
    /// Returns the entries `(variable, value)` with the structurals first, and their count.
    fn pivot_row(&self, rho: &[(usize, f64)], leave: usize, above: bool) -> (Vec<(usize, f64)>, usize) {
        let (n, m) = (self.n, self.m);
        let p = self.p;
        // The row product runs over the rows of `rho` unless their entries (`0x134aa0` counts them) outnumber the structurals n by 1.8 (`0xc5355`, `ctx+0x68`);
        // then it is columnwise, in index order.
        let est: usize = rho.iter().map(|e| self.rrows[e.0].len()).sum();
        let mut acc = vec![0.0f64; n + m];
        let mut touched = vec![false; n + m];
        let mut order: Vec<usize> = Vec::new();
        let sparse = !((n as f64) < 1.8 * est as f64);
        for &(i, rv) in rho {
            for &j in &self.rrows[i] {
                let a = p.cols[j].iter().find(|e| e.0 == i).map_or(0.0, |e| e.1);
                acc[j] = acc[j] + rv * a;
                if !touched[j] {
                    touched[j] = true;
                    if sparse {
                        order.push(j);
                    }
                }
            }
        }
        if !sparse {
            for j in 0..n {
                if touched[j] {
                    order.push(j);
                }
            }
        }
        // Entries that cancelled to (nearly) nothing are dropped (the smallest captured entry is 1.04e-13).
        order.retain(|&j| acc[j].abs() > DROP_TOL);
        let nst = order.len();
        for &(i, rv) in rho {
            let j = n + i;
            if self.slack_nb[i] || j == leave {
                acc[j] = rv * p.cols[j][0].1;
                order.push(j);
            }
        }
        let sdir = if above { -1.0 } else { 1.0 };
        let out = order
            .into_iter()
            .map(|j| {
                let s = if self.st[j] == UB { -1.0 } else { 1.0 };
                (j, acc[j] * (s * sdir))
            })
            .collect();
        (out, nst)
    }

    /// The bound-flipping Harris ratio test (`0xc74a0`, `ctx+0x468 != 0`), `specs/dual-simplex.md`.
    fn ratio_test(&self, row: &[(usize, f64)], esi: usize, leave: usize, above: bool, xbr: f64) -> Ratio {
        let p = self.p;
        let delta = DUAL_TOL;
        let (vlb, vub) = (p.lo[leave], p.up[leave]);
        let xmm10 = if above { vub - xbr } else { xbr - vlb };
        let ptol = self.ptol();
        let xmm12 = -2.0 * xmm10;
        let xmm7 = ptol - xmm10;
        let cnt = row.len();
        let mut bound = -1e100f64;
        let (mut best, mut best_i, mut bestc) = (-1e100f64, -1i64, 0.0f64);
        let (mut xm13, mut xm9) = (0.0f64, 0.0f64);
        let mut ratio_arr: Vec<f64> = Vec::new();
        let mut con: Vec<f64> = Vec::new();
        let mut lst: Vec<usize> = Vec::new();
        macro_rules! track {
            ($i:expr, $ratio:expr, $c:expr) => {{
                let (ratio, c) = ($ratio, $c);
                let k = lst.len() as i64;
                if ratio > 1e-10 + best {
                    best_i = k;
                    best = ratio;
                    bestc = c;
                } else if ratio >= best - 1e-10 {
                    let old = best;
                    bestc = minsd(bestc, c);
                    if ratio > old {
                        best_i = k;
                    }
                    best = maxsd(ratio, old);
                }
                con.push(c);
                lst.push($i);
                ratio_arr.push(ratio);
            }};
        }
        for i in esi..cnt {
            let (j, v) = row[i];
            if v < -ALPHA_TOL && (self.st[j] == LB || self.st[j] == UB) {
                let dj = self.d[j];
                let ratio = dj / v;
                if bound > ratio {
                    xm9 = minsd(v, xm9);
                } else {
                    xm13 = minsd(v, xm13);
                    if dj >= 0.0 {
                        bound = maxsd(PIVOT_DELTA(delta, v, ratio), bound);
                    }
                    track!(i, ratio, xmm12);
                }
            }
        }
        let variant = self.level != 0 && self.iters > self.m;
        for i in 0..esi.min(cnt) {
            let (j, v) = row[i];
            if v < -ALPHA_TOL {
                let dj = self.d[j];
                let ratio = dj / v;
                if bound > ratio {
                    xm9 = minsd(v, xm9);
                    continue;
                }
                xm13 = minsd(v, xm13);
                let c;
                if variant && self.st[j] == LB && p.up[j] >= BIG {
                    if dj >= 0.0 {
                        bound = maxsd(PIVOT_DELTA(delta, v, ratio), bound);
                    }
                    c = xmm12;
                } else {
                    let c0 = (p.lo[j] - p.up[j]) * v;
                    if c0 > xmm7 {
                        if dj >= 0.0 {
                            bound = maxsd(PIVOT_DELTA(delta, v, ratio), bound);
                        }
                        c = minsd(xmm12, c0);
                    } else {
                        c = c0;
                    }
                }
                track!(i, ratio, c);
            }
        }
        let nn = lst.len();
        if nn == 0 {
            return Ratio { q: None, flips: Vec::new() };
        }
        if nn == 1 {
            return Ratio { q: Some(lst[0]), flips: Vec::new() };
        }
        let maxa = -minsd(xm13, xm9);
        let mut nflip = 0usize;
        if bestc + xmm10 < 0.0 {
            let mut s = xmm10;
            let (mut lo, mut hi) = (0usize, nn - 1);
            while lo < hi {
                let (pr, pc) = (ratio_arr[lo], con[lo]);
                let mut h = hi + 1;
                let mut l = lo as isize - 1;
                loop {
                    h -= 1;
                    while pr > ratio_arr[h] {
                        h -= 1;
                    }
                    while ratio_arr[h] == pr && con[h] > pc {
                        h -= 1;
                        while pr > ratio_arr[h] {
                            h -= 1;
                        }
                    }
                    l += 1;
                    while ratio_arr[l as usize] > pr || (ratio_arr[l as usize] == pr && pc > con[l as usize]) {
                        l += 1;
                    }
                    if l as usize >= h {
                        break;
                    }
                    let lu = l as usize;
                    if ratio_arr[h] > ratio_arr[lu] || con[lu] > con[h] {
                        ratio_arr.swap(lu, h);
                        con.swap(lu, h);
                        lst.swap(lu, h);
                    }
                }
                if s < 0.0 {
                    s += con[lo..=h].iter().fold(0.0, |a, c| a + c);
                } else {
                    s -= con[h + 1..=hi].iter().fold(0.0, |a, c| a + c);
                }
                if 0.0 > s {
                    lo = h + 1;
                } else {
                    hi = h;
                }
            }
            nflip = lo;
        } else if best_i > 0 {
            let b = best_i as usize;
            ratio_arr.swap(0, b);
            con.swap(0, b);
            lst.swap(0, b);
        }
        if nflip > 0 && !(self.d[row[lst[nflip]].0] > delta) {
            nflip = 0;
        }
        let pivtol = maxsd(0.0, maxsd(1e-13 * maxa, PIVOT_FLOOR));
        if !(pivtol <= maxa) {
            return Ratio { q: None, flips: Vec::new() };
        }
        let mut b = nflip;
        let mut kept: Vec<usize> = Vec::new();
        let mut xm3 = -1e100f64;
        let mut xm15 = bound;
        for t in b..nn {
            let e = lst[t];
            let (j, av) = row[e];
            let dj = self.d[j];
            if pivtol > -av {
                continue;
            }
            if xm15 * av > dj {
                let r = (delta + dj) / av;
                xm3 = maxsd(r, xm3);
                kept.push(e);
                if dj >= 0.0 {
                    xm15 = maxsd(r, xm15);
                }
            }
        }
        if kept.is_empty() && b > 0 {
            for t in 0..b {
                let e = lst[t];
                let (j, av) = row[e];
                let dj = self.d[j];
                if pivtol > -av {
                    continue;
                }
                if xm15 * av > dj {
                    let r = (delta + dj) / av;
                    xm3 = maxsd(r, xm3);
                    kept.push(e);
                    if dj >= 0.0 {
                        xm15 = maxsd(r, xm15);
                    }
                }
            }
            b = 0;
        }
        if kept.is_empty() {
            return Ratio { q: None, flips: Vec::new() };
        }
        let mut chosen = kept[0];
        let mut xm5 = 0.0f64;
        for &e in &kept {
            let (j, av) = row[e];
            let wgt = if j < self.n { p.cols[j].len() as f64 } else { 0.5 };
            if xm3 * av < self.d[j] {
                continue;
            }
            if xm5 * wgt > av {
                xm5 = av / wgt;
                chosen = e;
            }
        }
        Ratio { q: Some(chosen), flips: lst[..b].to_vec() }
    }
}

/// `delta/a + ratio`: the reference adds two terms (not `(delta + d)/a`).
#[allow(non_snake_case)]
fn PIVOT_DELTA(delta: f64, v: f64, ratio: f64) -> f64 {
    delta / v + ratio
}

/// Result of one call of the iteration function.
enum Step {
    /// A basis change: leaving row, leaving and entering variable, the pivot row (reference frame) and the data of the post-iteration passes.
    Pivot(Box<Pivot>),
    /// No leaving row.
    NoRow,
    /// The ratio test found no entering variable: the dual is unbounded (the primal infeasible).
    Unbounded,
}

struct Pivot {
    r: usize,
    leave: usize,
    q: usize,
    row: Vec<(usize, f64)>,
    /// Flipped variables, with their statuses before the flip.
    flips: Vec<usize>,
    /// The pivot column without the pivot entry, tau at the same positions, the pivot alpha and tau_r, true frame.
    col: Vec<(usize, f64)>,
    tau: Vec<f64>,
    alpha: f64,
    w: f64,
    /// Value of the entering variable before the basis change (`ctx+0x318`).
    t0: f64,
}

impl<'a> Dual<'a> {
    /// One call of the iteration function `0xc4fe0`.
    fn iteration(&mut self, k: usize, dump: &mut Option<serde_json::Value>) -> Step {
        let p = self.p;
        let (n, m) = (self.n, self.m);
        let snap = crate::dump::on("dsx").then(|| self.snapshot(k));
        let row_sel = self.price();
        let Some(r) = row_sel else {
            if let Some(mut s) = snap {
                s["r"] = serde_json::json!(-1);
                *dump = Some(s);
            }
            return Step::NoRow;
        };
        let leave = self.head[r];
        let xbr = self.xb[r];
        let above = xbr > p.up[leave];
        // Pivot-row BTRAN, drop of the small entries.
        let rho_full = self.mirror.btran(&[(r, 1.0)], true, false);
        let rho: Vec<(usize, f64)> = rho_full.iter().copied().filter(|e| e.1.abs() > DROP_TOL).collect();
        let (row, esi) = self.pivot_row(&rho, leave, above);
        let rt = self.ratio_test(&row, esi, leave, above, xbr);
        let mut rec = snap;
        if let Some(s) = rec.as_mut() {
            s["r"] = serde_json::json!(r);
            s["var"] = serde_json::json!(leave);
            s["rho"] = crate::dump::sparse(&rho);
            s["row"] = crate::dump::sparse(&row);
        }
        let Some(qe) = rt.q else {
            if let Some(s) = rec.as_mut() {
                s["q"] = serde_json::json!(-1);
            }
            *dump = rec;
            return Step::Unbounded;
        };
        let (q, aq_row) = row[qe];
        let flips: Vec<usize> = rt.flips.iter().map(|&e| row[e].0).collect();
        if let Some(s) = rec.as_mut() {
            s["q"] = serde_json::json!(q);
            s["alpha_row"] = serde_json::json!(aq_row);
            s["flips"] = serde_json::json!(flips);
        }
        // tau = B^-1 rho first, then the pivot column (the update keeps the spike of the last FTRAN).
        let rho_rows: Vec<usize> = rho.iter().map(|e| e.0).collect();
        let rho_vals: Vec<f64> = rho.iter().map(|e| e.1).collect();
        let tau_out = self.mirror.ftran(&rho_rows, &rho_vals);
        let mut tau_dense = vec![0.0; m];
        for &(i, v) in &tau_out {
            tau_dense[i] = v;
        }
        let rows_q: Vec<usize> = p.fcols[q].iter().map(|&(i, _)| i).collect();
        let vals_q: Vec<f64> = p.fcols[q].iter().map(|&(_, a)| a).collect();
        let col_out = self.mirror.ftran(&rows_q, &vals_q);
        let mut alpha_dense = vec![0.0; m];
        for &(i, v) in &col_out {
            alpha_dense[i] = v;
        }
        let alpha = alpha_dense[r];
        // Pivot check: the row's value in the true frame against the column's.
        let sq = if self.st[q] == UB { -1.0 } else { 1.0 };
        let sdir = if above { -1.0 } else { 1.0 };
        let alpha_row_true = aq_row * (sq * sdir);
        let rel = (alpha_row_true - alpha).abs() / (alpha.abs() + 1.0);
        if rel > REL_PIVOT_TOL {
            self.request = 6;
        }
        if let Some(s) = rec.as_mut() {
            s["alpha_col"] = serde_json::json!(alpha);
            s["tau_r"] = serde_json::json!(tau_dense[r]);
            s["w"] = serde_json::json!(rho_norm(&rho));
            s["col"] = crate::dump::sparse(&col_out);
            s["tau"] = crate::dump::sparse(&tau_out);
        }
        // Column without the pivot entry: the last entry moves into the pivot's slot.
        let mut col: Vec<(usize, f64)> = col_out.clone();
        if let Some(k) = col.iter().position(|e| e.0 == r) {
            col.swap_remove(k);
        }
        let tau: Vec<f64> = col.iter().map(|e| tau_dense[e.0]).collect();
        // Bound flips: x_B and the statuses (`0xcd720`).
        if !flips.is_empty() {
            let mut rhs = vec![0.0; m];
            for &j in &flips {
                let delta = if self.st[j] == LB { p.up[j] - p.lo[j] } else { p.lo[j] - p.up[j] };
                for &(i, a) in &p.cols[j] {
                    rhs[i] += a * delta;
                }
            }
            self.factor.ftran(&mut rhs);
            let mut changed = Vec::new();
            for pos in 0..m {
                if rhs[pos] != 0.0 {
                    self.xb[pos] -= rhs[pos];
                    changed.push(pos);
                }
            }
            for &j in &flips {
                self.st[j] = if self.st[j] == LB { UB } else { LB };
                self.d[j] = -self.d[j];
            }
            self.touch_rows(&changed);
        }
        let t0 = self.nonbasic_value(q);
        // Basis change.
        let new_st_leave = if p.lo[leave] == p.up[leave] { FIXED } else if above { UB } else { LB };
        self.st[leave] = new_st_leave;
        if leave >= self.n {
            self.slack_nb[leave - self.n] = true;
        }
        self.st[q] = r as i32;
        self.head[r] = q;
        // Row-wise copy: in a row holding both columns the leaving one takes the entering one's place; elsewhere the entering one is deleted (the last
        // entry moves into its slot) and the leaving one is appended.
        let leave_in = leave < n && p.lo[leave] != p.up[leave];
        if q < n {
            for &(i, _) in &p.cols[q] {
                if let Some(kk) = self.rrows[i].iter().position(|&j| j == q) {
                    if leave_in && p.cols[leave].iter().any(|e| e.0 == i) {
                        self.rrows[i][kk] = leave;
                    } else {
                        self.rrows[i].swap_remove(kk);
                    }
                }
            }
        }
        if leave_in {
            for &(i, _) in &p.cols[leave] {
                if !self.rrows[i].contains(&leave) {
                    self.rrows[i].push(leave);
                }
            }
        }
        self.iters += 1;
        *dump = rec;
        let _ = flips.len();
        // The factor's update is made by the post-iteration pass (it refactorizes instead when one is due).
        let alpha_vec = alpha_dense;
        self.pending_alpha = Some(alpha_vec);
        Step::Pivot(Box::new(Pivot { r, leave, q, row, flips, col, tau, alpha, w: rho_norm(&rho), t0 }))
    }
}

impl<'a> Dual<'a> {
    fn snapshot(&self, k: usize) -> serde_json::Value {
        serde_json::json!({
            "k": k, "lvl": self.level, "req": self.request, "it": self.iters as f64, "head": self.head, "xB": self.xb, "d": self.d, "W": self.w,
            "st": self.st, "flags": self.flags, "prio": self.prio,
            "list": if self.value_mode { self.list.clone() } else { self.list[..self.walk_len.min(self.list.len())].to_vec() },
            "thr": self.thr, "p40": self.value_mode as u8, "cur": self.cursor,
        })
    }

    /// The no-leaving-row handler (`0x1aa280`): a refactorization request, optimality or infeasibility.
    fn no_row(&mut self) {
        if self.updates != 0 {
            // `ctx+0x3f8 != 0` in all captures: confirm on a fresh factor.
            self.request = 4;
        } else if self.level == 0 {
            self.status = Some(Status::Optimal);
        } else {
            let p = self.p;
            for j in 0..self.n + self.m {
                let bad = (self.st[j] == UB && p.up[j] >= BIG && self.d[j] > DUAL_TOL) || (self.st[j] == LB && p.lo[j] <= -BIG && self.d[j] > DUAL_TOL);
                if bad {
                    self.status = Some(Status::Infeasible);
                    return;
                }
            }
            self.level = 0;
        }
    }

    /// Level of the tolerance (`0xe45c0`, after a refactorization): with request 4 the nonbasic variables are first put on the dual-feasible bound.
    fn set_level(&mut self) {
        let p = self.p;
        if self.request == 4 {
            for j in 0..self.n + self.m {
                if self.st[j] == LB && self.d[j] < -DUAL_TOL {
                    self.st[j] = UB;
                    self.d[j] = -self.d[j];
                } else if self.st[j] == UB && p.up[j] >= BIG && self.d[j] <= DUAL_TOL {
                    self.st[j] = LB;
                    self.d[j] = -self.d[j];
                }
            }
        }
        self.level = (0..self.n + self.m).filter(|&j| self.st[j] == UB && p.up[j] >= BIG).count();
    }

    /// The post-iteration passes `0xf6c82..0xf6d0a`: factor update or refactorization, the dual update, the primal and weight update.
    fn post(&mut self, pv: Option<&Pivot>, schedule_due: bool) -> Result<(), SolveError> {
        let refactor = self.request != 0 || schedule_due;
        if let Some(pv) = pv {
            self.last_cap = pv.col.len();
            let alpha = self.pending_alpha.take().unwrap();
            if refactor {
                // The weights of the pivot are updated before the factorization renumbers the positions; the basic solution is recomputed.
                self.update_weights(pv);
                self.refactor()?;
            } else {
                self.factor.update(pv.r, alpha);
                self.mirror.update(pv.r);
                self.updates += 1;
            }
            // d update.
            let aq = pv.row.iter().find(|e| e.0 == pv.q).unwrap().1;
            let theta = self.d[pv.q] / aq;
            for &(j, a) in &pv.row {
                let a = if pv.flips.contains(&j) { -a } else { a };
                self.d[j] = self.d[j] - theta * a;
            }
            self.d[pv.leave] = -theta;
        } else if refactor {
            self.refactor()?;
        }
        if refactor {
            if self.request == 0 {
                self.request = 10;
            }
            self.recompute_d();
            let level_before = self.level;
            self.set_level();
            self.request = 0;
            self.recompute_xb();
            self.recompute_flags(level_before != 0, pv.is_some());
        } else if let Some(pv) = pv {
            self.update_primal_and_weights(pv);
        }
        Ok(())
    }

    /// x_B and DSE weight update (`0xc0050`, fast path), bit for bit as in the spec.
    fn update_weights(&mut self, pv: &Pivot) {
        let (r, alpha, w) = (pv.r, pv.alpha, pv.w);
        if w == 0.0 {
            self.w[r] = W_EPS;
        } else {
            let a2 = if alpha >= 0.0 { maxsd(alpha, W_EPS) } else { minsd(alpha, -W_EPS) };
            self.w[r] = (w / a2) / a2;
            for (j, &(i, cv)) in pv.col.iter().enumerate() {
                let a = cv / a2;
                let t = pv.tau[j];
                let x = ((w * a) * a - (2.0 * a) * t) + self.w[i];
                self.w[i] = maxsd(W_EPS, x);
            }
        }
    }

    fn update_primal_and_weights(&mut self, pv: &Pivot) {
        self.update_weights(pv);
        let (r, alpha) = (pv.r, pv.alpha);
        let mut changed: Vec<usize> = Vec::new();
        let v = pv.leave;
        let b = if self.st[v] == LB || self.st[v] == FIXED { self.p.lo[v] } else { self.p.up[v] };
        let theta = (self.xb[r] - b) / alpha;
        self.xb[r] = pv.t0 + theta;
        for &(i, cv) in &pv.col {
            self.xb[i] = self.xb[i] - cv * theta;
            changed.push(i);
        }
        changed.push(r);
        // Running average of the changed-row count and the growth of the list threshold (`0xc0532`, `0xc06d7`): while more rows are on the list than
        // that average and than m^0.45 the threshold grows by 10 percent per iteration (before the rows are touched).
        let ninf = self.list.len() as f64;
        let it = self.iters as f64;
        self.avg_touch = (2.0 * changed.len() as f64 + (it - 2.0) * self.avg_touch) / it;
        if self.value_mode {
            if ninf > self.avg_touch.max((self.m as f64).powf(0.45).floor()) {
                self.thr *= 1.1;
            }
            // The pivot row leaves the list at the next pricing call whatever its value (`0xc0732`); its priority is refreshed for the entering variable.
            let q = self.head[r];
            let d = self.xb[r] - self.p.lo[q];
            self.prio[r] = d * d / self.w[r];
            self.flags[r] = 3;
            changed.pop();
            self.touch_rows(&changed);
        } else {
            self.touch_rows(&changed);
            if self.pend < 0 {
                // Too many rows outside the walk became infeasible (`0xc2248`): the list is rebuilt from the current values.
                self.rebuilds -= 1;
                self.merge_walk_list();
            }
        }
    }
}

/// The squared norm of rho as the reference sums it (`0x11f9..0x11fa5`, the value `ctx+0x310`): sequentially up to 6 entries, else in the vectorised loop
/// of eight entries per pass (four accumulators of two lanes, combined pairwise, then lane 0 + lane 1, then the tail).
fn rho_norm(rho: &[(usize, f64)]) -> f64 {
    let sq: Vec<f64> = rho.iter().map(|e| e.1 * e.1).collect();
    let n = sq.len();
    if n <= 6 {
        return sq.iter().fold(0.0, |a, x| a + x);
    }
    let end = n - (n & 7);
    let mut acc = [[0.0f64; 2]; 4];
    let mut i = 0;
    while i < end {
        for q in 0..4 {
            acc[q][0] += sq[i + 2 * q];
            acc[q][1] += sq[i + 2 * q + 1];
        }
        i += 8;
    }
    let lane = |l: usize| (acc[0][l] + acc[1][l]) + (acc[2][l] + acc[3][l]);
    let mut t = lane(0) + lane(1);
    for x in &sq[end..] {
        t += x;
    }
    t
}

/// Dual simplex from a cold start (the crash basis) or a supplied basis.
pub fn solve(model: &LpModel, start: Option<&Basis>, opts: &SolveOptions) -> Result<SolveResult, SolveError> {
    let p = build(model);
    let (n, m) = (p.n, p.m);
    let basis = match start {
        Some(b) => {
            let mut b = b.clone();
            b.validate(n, m).map_err(SolveError::InvalidBasis)?;
            b.basic.sort_unstable();
            b
        }
        None => Basis { basic: dual_crash_basic(&p), at_upper: vec![false; n + m] },
    };
    basis.validate(n, m).map_err(SolveError::InvalidBasis)?;
    crate::lu_order::reset_pair_gate();
    let (order, stc, rows) = reference_order(&p, &basis.basic, 1);
    let mirror = RefStruct::fresh(m, n, &stc, &rows);
    let cols: Vec<&Column> = order.iter().map(|&j| &p.cols[j]).collect();
    let factor = Factor::new(m, &cols).map_err(|_| SolveError::SingularBasis)?;
    let mut st = vec![LB; n + m];
    for (pos, &j) in order.iter().enumerate() {
        st[j] = pos as i32;
    }
    let cap = ((m as f64).powf(0.7).floor() as i64).max(20);
    let mut s = Dual {
        p: &p,
        n,
        m,
        head: order,
        st,
        xb: vec![0.0; m],
        d: vec![0.0; n + m],
        w: vec![1.0; m],
        mirror,
        factor,
        rrows: vec![Vec::new(); m],
        refactorizations: 1,
        level: 0,
        request: 0,
        iters: 0,
        prio: vec![0.0; m],
        flags: vec![2; m],
        list: Vec::new(),
        thr: 0.0,
        value_mode: false,
        cursor: 0,
        walk_len: 0,
        slack_nb: vec![false; m],
        pend: cap,
        cap,
        rebuilds: 40,
        first_recompute: true,
        avg_touch: 0.0,
        last_cap: 0,
        updates: 0,
        pending_alpha: None,
        status: None,
    };
    // Initial statuses of the nonbasic variables: the dual-feasible bound (an infinite bound is the artificial one, value +-1e30).
    for j in 0..n + m {
        if s.st[j] >= 0 {
            continue;
        }
        let (lo, up) = (p.lo[j], p.up[j]);
        s.st[j] = if lo == up {
            FIXED
        } else if is_inf(lo) && is_inf(up) {
            FREE
        } else {
            LB
        };
    }
    let dt = s.true_d();
    for j in 0..n + m {
        if s.st[j] >= 0 || s.st[j] == FIXED || s.st[j] == FREE {
            continue;
        }
        if dt[j] < -DUAL_TOL {
            s.st[j] = UB;
        } else if dt[j] > DUAL_TOL {
            s.st[j] = LB;
        }
        s.d[j] = if s.st[j] == UB { -dt[j] } else { dt[j] };
    }
    s.build_rrows();
    s.slack_nb = (0..m).map(|i| s.st[n + i] < 0).collect();
    s.level = (0..n + m).filter(|&j| s.st[j] == UB && p.up[j] >= BIG).count();
    // Initial weights: the squared norms of the rows of B^-1.
    for r in 0..m {
        let out = s.mirror.btran(&[(r, 1.0)], false, false);
        s.w[r] = out.iter().fold(0.0, |a, e| a + e.1 * e.1);
    }
    s.recompute_xb();
    s.recompute_flags(false, false);
    let schedule = opts.refactor_schedule.clone();
    let mut trace = Vec::new();
    let mut k = 0usize;
    let status;
    loop {
        k += 1;
        if s.iters >= opts.max_iterations {
            status = Status::IterationLimit;
            break;
        }
        let mut dump = None;
        let step = s.iteration(k, &mut dump);
        if let Some(rec) = dump {
            crate::dump::emit("dsx", rec);
        }
        match step {
            Step::Pivot(pv) => {
                if opts.record_trace {
                    trace.push(TraceEvent { iteration: s.iters, entering: pv.q, leaving: Some(pv.leave), position: Some(pv.r), step: 0.0 });
                }
                let due = match &schedule {
                    Some(sch) => sch.contains(&s.iters),
                    None => s.mirror.refactor_due(),
                };
                s.post(Some(&pv), due)?;
            }
            Step::NoRow => {
                s.no_row();
                if let Some(stt) = s.status {
                    status = stt;
                    break;
                }
                s.post(None, false)?;
            }
            Step::Unbounded => {
                status = Status::Infeasible;
                break;
            }
        }
    }
    // Read-out.
    let sgn = if model.sense == Sense::Maximize { -1.0 } else { 1.0 };
    let mut x = vec![0.0; n + m];
    for j in 0..n + m {
        x[j] = if s.st[j] >= 0 { s.xb[s.st[j] as usize] } else { s.nonbasic_value(j) };
    }
    let mut y: Vec<f64> = s.head.iter().map(|&j| p.cost[j]).collect();
    s.factor.btran(&mut y);
    let dtrue = s.true_d();
    let nm = model.num_cols();
    let xs: Vec<f64> = x[..nm].to_vec();
    let at_upper: Vec<bool> = s.st.iter().map(|&t| t == UB).collect();
    Ok(SolveResult {
        status,
        objective: model.objective(&xs),
        x: xs,
        dual: y.iter().map(|v| sgn * v).collect(),
        reduced_cost: dtrue[..nm].iter().map(|v| sgn * v).collect(),
        basis: Basis { basic: s.head.clone(), at_upper },
        iterations: s.iters,
        refactorizations: s.refactorizations,
        trace,
        farkas: None,
        ray: None,
    })
}
