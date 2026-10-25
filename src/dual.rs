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
use crate::simplex::{
    build, ref_sort, Internal, SolveError, SolveOptions, SolveResult, Status, TraceEvent,
};

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

/// The reference's integer generator (`0x88a870`): a 64-bit congruential step, the high word tempered, 31 bits.
fn lcg_int(state: &mut u64) -> u32 {
    *state = state
        .wrapping_mul(0x5851f42d4c957f2d)
        .wrapping_add(0x14057b7ef767814f);
    let mut e = (*state >> 32) as u32;
    e ^= e >> 11;
    e ^= (e << 7) & 0x9d2c5680;
    e ^= (e << 15) & 0xefc60000;
    e ^= e >> 18;
    e & 0x7fff_ffff
}

fn is_inf(v: f64) -> bool {
    v.abs() >= INF
}

/// The pivot threshold asked for by the accuracy of the x_B recompute (`0xc1880`, `0xc3e42..0xc3ed8`, `0xc18a5`): `m0` = max |x_B| / R and `e2` = max |rhs - B x_B| / R,
/// R = max(1e-9, max |rhs|). The first rung whose `m0` or `e2` bound is exceeded gives the candidate; the driver keeps the larger of it and `ctx+0x390`.
fn accuracy_threshold(m0: f64, e2: f64) -> f64 {
    const RUNGS: [(f64, f64, f64); 5] = [
        (1e10, 1e-8, 0.5),
        (1e9, 4e-9, 0.25),
        (1e8, 1e-9, 0.125),
        (1e7, 4e-10, 0.0625),
        (1e6, 1e-10, 0.03125),
    ];
    for (a, b, t) in RUNGS {
        if m0 > a || e2 > b {
            return t;
        }
    }
    1e-4
}

fn minsd(a: f64, b: f64) -> f64 {
    if a < b {
        a
    } else {
        b
    }
}

fn maxsd(a: f64, b: f64) -> f64 {
    if a > b {
        a
    } else {
        b
    }
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
    /// The relative pivot threshold of the LU (`ctx+0x390`, 2^-7 at the start): raised after a recompute of x_B by `accuracy_threshold`, never lowered.
    lu_thr: f64,
    /// The pricing mode `ctx+0x3ec`: 1 for dual steepest edge, 3 after the switch of `0xf770d` (weights restarted at 1).
    pmode: u8,
    /// The number of infeasible rows at the last value-mode recompute (`P+0x3c`).
    n_infeas: usize,
    /// The reference framework of mode 3 (`ctx+0x418`): the variables nonbasic when the weights were restarted.
    dflag: Vec<bool>,
    /// `SimplexPricing` left at its default (-1): the release profile, which may switch the pricing mode.
    switch_enabled: bool,
    /// Tolerance level (`ctx+0x2cc`): the number of nonbasic variables at an artificial upper bound.
    level: usize,
    /// Pending refactorization request (`ctx+0x294`).
    request: u32,
    /// Iterations made (`ctx+0x2a0`).
    iters: usize,
    // Pricing object P (`ctx+0x478`).
    prio: Vec<f64>,
    /// Pricing type `P+0x44`: 0 plain dual steepest edge priorities, 1 scaled by `P70[var] / P68[row]`, 2 by `1 / P68[row]`.
    ptype: u8,
    /// Per-row weights `P+0x68` (by position) and per-variable constants `P+0x70` of the scaled priorities.
    p68: Vec<f64>,
    p70: Vec<f64>,
    /// Touch increment `P+0x60`: the number of entries of the pivot row's BTRAN result.
    p60: f64,
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
    /// A supplied start basis: the reference runs it in the pricing mode 3 (`ctx+0x3ec`, `ctx+0x298` = 2) with unit weights, and the first
    /// recompute starts the value mode (`ctx+0x3fc` = -1).
    warm: bool,
    /// Running average of the number of rows touched per iteration.
    avg_touch: f64,
    /// Entry count of the last pivot column (`*(ctx+0x210)`).
    last_cap: usize,
    /// Updates of the factor since the last refactorization.
    updates: usize,
    /// The update limit of the factor object once the value mode has started (`F+0x188`).
    upd_limit: Option<usize>,
    /// The pivot column of the last iteration by position (for the factor's update).
    pending_alpha: Option<Vec<f64>>,
    /// Final status once the loop ends.
    status: Option<Status>,
    /// Number of structural columns of the model (the length of a returned ray).
    nm: usize,
    ray: Option<Vec<f64>>,
    /// The read-out of a solve that stopped inside the artificial-bound phase takes the nonbasic variables at artificial bounds as 0 instead of +-1e30.
    zero_artificial: bool,
    /// Costs in force: the problem's costs plus the shifts of the ratio test (`*(P+0)`).
    cost: Vec<f64>,
    /// The accumulated cost shifts per variable: the reference keeps them apart and forms `cost = base + shift` (the sum of the shifts is rounded once).
    cshift: Vec<f64>,
    /// State of the reference's random generator (`0x88a870`, `0x88a8d0`), seeded 101.
    rng: u64,
    /// Perturbation flag `ctx+0x338`: 1 after a small cost shift, 2 after a large one or one against an infinite bound.
    pert: u8,
    /// Objective offset of the cost shifts (`ctx+0x130`).
    obj_off: f64,
    /// `ctx+0x3fc` (-1 once set): the second selection among the bound flips of the ratio test. Set at a refactorization (`0xe6dcc`) when the
    /// perturbation flag is 3 or more or 3 iterations exceed the row count (`ctx+0x178`, which also sets it, is 0 in every capture).
    flip_select: bool,
    /// The reference's work counter (`ctx+0x118`); in the dual path every estimate is on the scale 1.3.
    work: crate::work::Work,
    /// Entry count of the pivot column when the pending factorization follows a pivot (its weight update and spike discard are booked then).
    pivot_cols: Option<usize>,
    /// `nnz(A)` of the structural columns (`ctx+0x3d0`), for the time clause.
    nnz_a: usize,
    /// Cost perturbation state (`reverse-engineer/specs/dual-simplex.md`, "Cost perturbation"): the scale `ctx+0x348` (0: off), the count of iterations
    /// with an entering `d` of at most delta (`ctx+0x330`) with its value (`0x334`) and the iteration count (`0x2b0`) at the last stall check, the progress
    /// flag `0x110`, the best objective `0x140`, the perturbations left `0x364`, the parameter value, the objective scale `ctx+0xd8` and the cost
    /// base that perturbations accumulate in.
    pstep: f64,
    stall: i64,
    stall_snap: i64,
    stall_it: f64,
    prog_on: bool,
    best_obj: f64,
    pert_left: i64,
    pert_param: f64,
    obj_scale: f64,
    pcost: Vec<f64>,
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
        crate::dump::emit(
            "dcr",
            serde_json::json!({"order": order, "counts": order.iter().map(|&r| row_cols[r].len()).collect::<Vec<_>>()}),
        );
    }
    let mut cost: Vec<f64> = p.cost[..n].to_vec();
    let mut marked = vec![false; m];
    // Columns that start at an infinite bound (negative cost and infinite upper bound, or positive cost and infinite lower bound) are flagged from the start.
    let mut flagged: Vec<bool> = (0..n)
        .map(|j| (p.cost[j] < 0.0 && p.up[j] >= BIG) || (p.cost[j] >= 0.0 && p.lo[j] <= -BIG))
        .collect();
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
        let (mut best_ratio, mut best_count, mut best): (f64, usize, Option<(usize, f64)>) =
            (1e100, m, None);
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
        // Pass 2 (`0xeb5c0`, entry test `0xeb61b`, choice `0xed270`): a column with a cost left (above 1e-10) is passed over unless its lower bound is
        // infinite; it must be nonbasic, |a| >= 1e-6 and have at most five entries. A column with an infinite lower bound has the key `count - 2m` (it
        // wins over every other); the best key kept is the count of the best column, the key compared is the candidate's.
        let mut best: Option<usize> = None;
        let (mut best_cnt, mut best_abs) = (5i64, 0.0f64);
        for &(j, a) in &row_cols[r] {
            let free_lb = p.lo[j] <= -BIG;
            if (cost[j].abs() > 1e-10 && !free_lb)
                || taken[j]
                || a.abs() < 1e-6
                || p.cols[j].len() > 5
            {
                continue;
            }
            let cnt = p.cols[j].len() as i64;
            let key = if free_lb { cnt - 2 * m as i64 } else { cnt };
            if key < best_cnt || (key == best_cnt && a.abs() > best_abs) {
                best = Some(j);
                best_cnt = cnt;
                best_abs = a.abs();
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
    crash_pass3(p, &mut basic, &mut cost, &taken, &is_eq, &row_cols);
    basic
}

/// Pass 3 of the start basis (`0xeb7ee..0xeba97`, assignment `0xeb9db`; `replay/gdb_dual_crash3.py`): every row that is not an equality and still has its slack
/// is looked at in ascending order, in the reference's row sign (a `>=` row is negated). The nonbasic status of a column is `at lower` when its cost is at
/// least 1e-6, `at upper` otherwise (a fixed column keeps its bound). The row is abandoned when it has a basic column, a bound of 1e30 or more on a
/// column, or a column with a zero cost and a zero bound that points the wrong way. Otherwise `sum = a x - rhs` at those bounds must exceed 1e-6.
/// Candidates are the columns with `a cost < 0` and ratio `-cost / a`: a column with a bound above 1e-10 in absolute value competes for the largest ratio
/// (`jmax`), the others for the smallest (`jmin`); the largest wins when it is smaller than the smallest, else (ties included) the smallest, and nothing is taken when
/// the winner is missing (a smallest without a largest is dropped). The chosen column replaces the slack (`cost = 0`), the costs of the other columns of the row
/// grow by `a ratio` and their statuses follow the sign of the cost.
fn crash_pass3(
    p: &Internal,
    basic: &mut [usize],
    cost: &mut [f64],
    taken: &[bool],
    is_eq: &[bool],
    row_cols: &[Vec<(usize, f64)>],
) {
    let (n, m) = (p.n, p.m);
    // -4 fixed, -1 at lower, -2 at upper, >= 0 basic in that row.
    let mut st: Vec<i64> = (0..n + m)
        .map(|j| {
            if p.up[j] - p.lo[j] <= 0.0 {
                -4
            } else if p.cost[j] >= 1e-6 {
                -1
            } else {
                -2
            }
        })
        .collect();
    for (r, &j) in basic.iter().enumerate() {
        if j < n {
            st[j] = r as i64;
            st[n + r] = -4;
        }
    }
    let _ = taken;
    for r in 0..m {
        // A free row (free slack) keeps its slack: a nonbasic free slack with a d would be dual infeasible for good (no evidence for such a row).
        if is_eq[r] || basic[r] < n || (p.lo[n + r] <= -BIG && p.up[n + r] >= BIG) {
            continue;
        }
        let sg = p.cols[n + r][0].1;
        let (mut sum, mut max_r, mut min_r) = (0.0f64, 0.0f64, 1e100f64);
        let (mut jmax, mut jmin): (Option<usize>, Option<usize>) = (None, None);
        let mut abandoned = false;
        for &(j, a0) in &row_cols[r] {
            let a = sg * a0;
            let bound = match st[j] {
                -4 | -1 => p.lo[j],
                -2 => p.up[j],
                _ => {
                    abandoned = true;
                    break;
                }
            };
            if bound.abs() >= BIG {
                abandoned = true;
                break;
            }
            sum += bound * a;
            if st[j] == -4 {
                continue;
            }
            let c = cost[j];
            if bound.abs() <= 1e-10 && c == 0.0 {
                if (st[j] == -1 && a < 0.0) || (st[j] == -2 && a > 0.0) {
                    abandoned = true;
                    break;
                }
                continue;
            }
            if a * c < 0.0 {
                let ratio = -c / a;
                if bound.abs() <= 1e-10 {
                    if min_r > ratio {
                        jmin = Some(j);
                    }
                    min_r = ratio.min(min_r);
                } else {
                    if ratio > max_r {
                        jmax = Some(j);
                    }
                    max_r = ratio.max(max_r);
                }
            }
        }
        if abandoned || sum - sg * p.rhs[r] <= 1e-6 {
            continue;
        }
        let (chosen, ratio) = if max_r >= min_r {
            (jmin, min_r)
        } else {
            (jmax, max_r)
        };
        let Some(c) = chosen else { continue };
        if crate::dump::on("dcr") {
            crate::dump::emit("dcr", serde_json::json!({"col": c, "row": r, "pass": 3}));
        }
        for &(k, a0) in &row_cols[r] {
            if k == c {
                st[k] = r as i64;
                cost[k] = 0.0;
                continue;
            }
            cost[k] += sg * a0 * ratio;
            if st[k] != -4 {
                if cost[k] < 0.0 {
                    if st[k] == -1 {
                        st[k] = -2;
                    }
                } else if st[k] == -2 {
                    st[k] = -1;
                }
            }
        }
        st[basic[r]] = -1;
        basic[r] = c;
    }
}

/// The reference's pivot order of the basis columns (`lu_order`), as in the primal path.
fn reference_order(
    p: &Internal,
    basic: &[usize],
    fac_no: usize,
) -> (
    Vec<usize>,
    crate::lu_order::Structure,
    Vec<(bool, Vec<usize>, Vec<f64>)>,
) {
    let n = p.n;
    let cols: Vec<(bool, Vec<usize>, Vec<f64>)> = basic
        .iter()
        .map(|&v| {
            (
                v >= n,
                p.fcols[v].iter().map(|&(r, _)| r).collect(),
                p.fcols[v].iter().map(|&(_, a)| a).collect(),
            )
        })
        .collect();
    crate::lu_order::set_factorization_number(fac_no);
    let st = crate::lu_order::factor_structure(p.m, &cols);
    let order = st.pivots.iter().map(|&(_, pos)| basic[pos]).collect();
    (order, st, cols)
}

impl<'a> Dual<'a> {
    fn nonbasic_value(&self, j: usize) -> f64 {
        if self.zero_artificial && self.at_artificial(j) {
            return 0.0;
        }
        match self.st[j] {
            LB | FIXED => self.p.lo[j],
            UB => self.p.up[j],
            _ => 0.0,
        }
    }

    /// Primal feasibility tolerance in force (`0x380`, or `0x370 * (1e-5 * 0x380)` while the level is not zero).
    fn ptol(&self) -> f64 {
        if self.level != 0 {
            BIG * (1e-5 * PRIMAL_TOL)
        } else {
            PRIMAL_TOL
        }
    }

    fn infeas(&self, r: usize) -> f64 {
        let v = self.head[r];
        let x = self.xb[r];
        let (lo, up) = (self.p.lo[v], self.p.up[v]);
        if lo > -BIG && x < lo {
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
        crate::lu_order::set_threshold(self.lu_thr);
        let (order, st, rows) = reference_order(self.p, &self.head, self.refactorizations + 1);
        let discard = if self.pivot_cols.is_some() {
            self.mirror.discard_work()
        } else {
            0.0
        };
        if crate::dump::on("luin") {
            crate::dump::emit(
                "luin",
                serde_json::json!({"n": self.refactorizations, "cols": self.head, "out": order, "prow": st.pivots.iter().map(|p| p.0).collect::<Vec<_>>(), "phase": st.phase, "thr": self.lu_thr}),
            );
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
        self.reset_pricing(&st, &rows);
        self.mirror = RefStruct::fresh(self.m, self.n, &st, &rows);
        self.mirror.upd_limit = self.upd_limit;
        self.factor_work(&st, &rows, discard);
        let cols: Vec<&Column> = self.head.iter().map(|&j| &self.p.cols[j]).collect();
        self.factor = Factor::new(self.m, &cols).map_err(|_| SolveError::SingularBasis)?;
        self.refactorizations += 1;
        self.updates = 0;
        self.slack_nb = (0..self.m).map(|i| self.st[self.n + i] < 0).collect();
        Ok(())
    }

    /// The free-variable pivots of the setup (`0xe88xx`, before the iteration loop; `replay/gdb_dual_free.py`): every free nonbasic structural, in
    /// ascending index, is pivoted into the basis. Its column is transformed; among the rows whose basic variable has a finite lower bound the entry of
    /// largest magnitude wins (first in the result's order on ties; below `PIVOT_FLOOR` the variable is deferred, unexercised). The leaving variable
    /// goes to its lower bound (fixed when the range is zero). The factor is updated, or refactorized when the work-based trigger is due.
    fn free_pivots(&mut self) -> Result<(), SolveError> {
        let (n, m) = (self.n, self.m);
        let free: Vec<usize> = (0..n).filter(|&j| self.st[j] == FREE).collect();
        for j in free {
            let rows_q: Vec<usize> = self.p.fcols[j].iter().map(|&(i, _)| i).collect();
            let vals_q: Vec<f64> = self.p.fcols[j].iter().map(|&(_, a)| a).collect();
            let mut col_out = self.mirror.ftran(&rows_q, &vals_q);
            if self.mirror.last_dense() {
                col_out.sort_by_key(|e| e.0);
            }
            let mut best: Option<(usize, f64)> = None;
            let mut mx = 0.0f64;
            for &(pos, v) in &col_out {
                if self.p.lo[self.head[pos]] > -BIG && v.abs() > mx {
                    mx = v.abs();
                    best = Some((pos, v));
                }
            }
            let Some((r, _)) = best.filter(|_| mx >= PIVOT_FLOOR) else {
                continue;
            };
            let mut alpha = vec![0.0; m];
            for &(i, v) in &col_out {
                alpha[i] = v;
            }
            let leave = self.head[r];
            self.st[leave] = if self.p.lo[leave] == self.p.up[leave] {
                FIXED
            } else {
                LB
            };
            self.st[j] = r as i32;
            self.head[r] = j;
            self.iters += 1;
            if self.mirror.refactor_due() {
                self.refactor()?;
            } else {
                self.factor.update(r, alpha);
                self.mirror.update(r);
                self.mirror.upd_work.clear();
                self.updates += 1;
            }
        }
        Ok(())
    }

    /// The work of a factorization and the installation of its factor (`0xd4190` from `0xf34f2`, the update routine `0xf15c0`), after the fresh mirror is built:
    /// the driver's sites, the spike discard and the weight update of the pivot that preceded it, then the factor object's own passes.
    fn factor_work(
        &mut self,
        st: &crate::lu_order::Structure,
        rows: &[(bool, Vec<usize>, Vec<f64>)],
        discard: f64,
    ) {
        use crate::work::{CTX, FAC};
        let lu = st.work(rows, self.mirror.has_levels());
        self.work.factor_driver(&lu, CTX);
        for (site, est) in lu {
            self.work.add(site, est, CTX);
        }
        if let Some(c) = self.pivot_cols.take() {
            if discard > 0.0 {
                self.work.add(0xf461f, discard, CTX);
            }
            self.work.add(0xbe793, 3.0 * c as f64, CTX);
        }
        let m3 = 3.0 * self.m as f64;
        self.work.add(0xf48ae, m3, FAC);
        // The scaled-priority pass `0x132460` (only with a pricing type above zero): `nnz(U)` and `3 m`.
        if self.ptype > 0 {
            let u = self.mirror.u_nnz();
            if u > 0 {
                self.work.add(0x1327dc, u as f64, FAC);
            }
            self.work.add(0x1329c1, m3, FAC);
        }
        for (site, est) in self.mirror.fresh_work() {
            self.work.add(site, est, CTX);
        }
        self.work.factor_done();
    }

    /// The row-wise copy of the nonbasic, non-fixed structurals (`ctx+0x1d0`), built once at the start, column by column in index order; the pivots update
    /// it in place and no factorization rebuilds it.
    fn build_rrows(&mut self) {
        for r in self.rrows.iter_mut() {
            r.clear();
        }
        for j in 0..self.n {
            if self.st[j] < 0 && self.st[j] != FIXED && self.p.lo[j] != self.p.up[j] {
                for &(i, _) in &self.p.cols[j] {
                    self.rrows[i].push(j);
                }
            }
        }
    }

    // ------------------------------------------------------------------ basic solution

    /// x_B from scratch (`0xc0c29` path of `0xc0050`): FTRAN of `b - N x_N`.
    ///
    /// `req` is the refactorization request in force (`ctx+0x294`) when the recompute follows a factorization, 0 otherwise. The accuracy of the solve then
    /// raises the LU pivot threshold for the next factorizations (`0xc1861`, `accuracy_threshold`).
    fn recompute_xb(&mut self, req: u32) {
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
        let est = self.mirror.ftran_dense_work(&r);
        self.work
            .add_from(0x1219d6, 0x120673, est, crate::work::FAC);
        // `0xc142b`: the solution read-out sums the columns of the nonbasic variables at a nonzero value.
        let cols: usize = (0..p.n)
            .filter(|&j| self.st[j] < 0 && self.nonbasic_value(j) != 0.0)
            .map(|j| p.cols[j].len())
            .sum();
        self.work.add(
            0xc142b,
            4.0 * p.m as f64 + 2.0 * cols as f64,
            crate::work::CTX,
        );
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
        // `0xc10ab..0xc1129`, `0xc1480..0xc1517`: R = max(1e-9, max|rhs|), the largest |x_B| and the residual of B x_B against the rhs, both relative to R;
        // the residual is read only when the factorization count is 4 mod 5 or the request is 4, 6 or 8 (`0xc1137..0xc11ce`, `0xc165b..0xc1699`).
        let rmax = r.iter().fold(1e-9f64, |a, v| a.max(v.abs()));
        let xmax = xb.iter().fold(0.0f64, |a, v| a.max(v.abs()));
        let e2 = if self.refactorizations % 5 == 4 || matches!(req, 4 | 6 | 8) {
            res.iter().fold(0.0f64, |a, v| a.max(v.abs())) / rmax
        } else {
            0.0
        };
        let t = accuracy_threshold(xmax / rmax, e2);
        if crate::dump::on("dualacc") {
            crate::dump::emit(
                "dualacc",
                serde_json::json!({"n": self.refactorizations, "req": req, "R": rmax, "xmax": xmax, "e2": e2, "cand": t, "thr": self.lu_thr}),
            );
        }
        if t > self.lu_thr {
            self.lu_thr = t;
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
        self.work
            .add_from(0x12de50, 0xe53a7, self.mirror.work, crate::work::FAC);
        for j in 0..self.n + self.m {
            if self.st[j] == FIXED {
                // A fixed nonbasic variable that has a d (it left the basis, or the start gave it one) gets its true d like one at its lower bound
                // (DEGEN2 d[700] after the refactorization at iteration 16); the fixed variables whose d is 0 keep it (PILOT4: 183 of 214).
                if self.d[j] != 0.0 {
                    self.d[j] = if j < self.n { self.cost[j] } else { dt[j] };
                }
            } else if self.st[j] < 0 {
                self.d[j] = if self.st[j] == UB { -dt[j] } else { dt[j] };
            } else {
                self.d[j] = 0.0;
            }
        }
    }

    /// d from the costs (`d_j = c_j - y a_j`, in the true frame) with y = B^-T c_B from the reference BTRAN; the terms are subtracted row by row in
    /// ascending row order, as the reference's row-wise product `0x136b00` does (rows with |y| <= 1e-13 are skipped).
    fn true_d(&mut self) -> Vec<f64> {
        let p = self.p;
        let inp: Vec<(usize, f64)> = self
            .head
            .iter()
            .enumerate()
            .map(|(pos, &j)| (pos, self.cost[j]))
            .filter(|e| e.1 != 0.0)
            .collect();
        let out = self.mirror.btran(&inp, false, true);
        if crate::dump::on("dual") {
            crate::dump::emit(
                "dual",
                serde_json::json!({"refactor": self.refactorizations, "iter": self.iters, "out": crate::dump::sparse(&out)}),
            );
        }
        let mut y = vec![0.0; p.m];
        for &(i, v) in out.iter() {
            y[i] = v;
        }
        let mut d = vec![0.0; p.n + p.m];
        for j in 0..p.n + p.m {
            if self.st[j] < 0 {
                // The reference accumulates row-wise (`0x136b00`): rows in ascending order, skipping |y_i| <= 1e-13 (`ctx+0x3a8`).
                let mut e: Vec<(usize, f64)> = p.fcols[j].clone();
                e.sort_by_key(|x| x.0);
                let mut acc = self.cost[j];
                for &(i, a) in &e {
                    if y[i].abs() > DROP_TOL {
                        acc -= y[i] * a;
                    }
                }
                d[j] = acc;
            }
        }
        d
    }

    // ------------------------------------------------------------------ pricing object

    /// Priority of row `r` with primal infeasibility `f` (`0xc0861`, `0xc23fe`, `0xcdc9`): `f^2 / W[r]`, scaled by the pricing type
    /// (the factor 1 of the reference's variable-class multiplier holds for continuous variables).
    fn pri(&self, r: usize, f: f64) -> f64 {
        let base = f * f / self.w[r];
        match self.ptype {
            1 => base * (self.p70[self.head[r]] / self.p68[r]),
            2 => base / self.p68[r],
            _ => base,
        }
    }

    /// Setup of the pricing type (`0xe7f7c..0xe7fe5`, `SimplexPricing` -1): type 1 for m >= 100 and n >= 100 unless the matrix is dense
    /// (more than 20 entries per column and m * nnz > 10 n^2); `P+0x70` is the column count plus the average count.
    fn init_pricing_type(&mut self) {
        let (n, m) = (self.n, self.m);
        let nnz: usize = self.p.cols[..n].iter().map(|c| c.len()).sum();
        let (nf, mf, z) = (n as f64, m as f64, nnz as f64);
        self.ptype = if m >= 100 && n >= 100 && (z <= 20.0 * nf || mf * z <= 10.0 * nf * nf) {
            1
        } else {
            0
        };
        if self.ptype > 0 {
            let avg = z / nf;
            self.p70 = (0..n + m)
                .map(|j| (if j < n { self.p.cols[j].len() } else { 1 }) as f64 + avg)
                .collect();
        }
    }

    /// Reset of the row weights at a factorization (`0x132460`, from `0xf48fd`): in pivot order, `w_k = 1` combined with the U entries of column `k`
    /// (rows above, descending) as `w + w_j - w w_j / m`, then backwards with the L entries of column `k`; the total decides between the
    /// capped weights `min(0.1 m + 100, w + 100)` and (total above `m (0.1 m + 100)`, or type 2 already) the constant cap with type 2.
    fn reset_pricing(
        &mut self,
        st: &crate::lu_order::Structure,
        cols: &[(bool, Vec<usize>, Vec<f64>)],
    ) {
        if self.ptype == 0 {
            return;
        }
        let m = self.m;
        let mut step_of_row = vec![usize::MAX; m];
        let mut step_of_pos = vec![usize::MAX; m];
        for (s, &(r, c)) in st.pivots.iter().enumerate() {
            step_of_row[r] = s;
            step_of_pos[c] = s;
        }
        // Forward lists (`0xf48fd` structure +0x70/+0x80): the earlier steps whose L column holds the row of step k, descending.
        let mut fwd: Vec<Vec<usize>> = vec![Vec::new(); m];
        for e in 0..m {
            for &r in &st.l_rows[e] {
                fwd[step_of_row[r]].push(e);
            }
        }
        // Backward lists (+0xb0/+0xd0): the later steps in the pivot row of step k: nucleus rows with their fill (`u_rows`), the others with their
        // original entries, by ascending basis position.
        let mut row_pos: Vec<Vec<usize>> = vec![Vec::new(); m];
        for (c, col) in cols.iter().enumerate() {
            for &i in &col.1 {
                row_pos[i].push(c);
            }
        }
        let bwd: Vec<Vec<usize>> = (0..m)
            .map(|k| {
                if st.phase[k] != 0 {
                    row_pos[st.pivots[k].0]
                        .iter()
                        .map(|&c| step_of_pos[c])
                        .filter(|&t| t > k)
                        .collect()
                } else {
                    st.u_rows[k]
                        .iter()
                        .map(|&(c, _)| step_of_pos[c])
                        .filter(|&t| t > k)
                        .collect()
                }
            })
            .collect();
        let n = m as f64;
        let mut w = vec![1.0f64; m];
        for k in 0..m {
            fwd[k].sort_unstable_by(|a, b| b.cmp(a));
            fwd[k].dedup();
            for &j in &fwd[k] {
                w[k] = (w[k] + w[j]) - (w[j] * w[k]) / n;
            }
        }
        for k in (0..m).rev() {
            for &j in &bwd[k] {
                w[k] = (w[k] + w[j]) - (w[j] * w[k]) / n;
            }
        }
        let total: f64 = w.iter().sum();
        let cap = 0.1 * n + 100.0;
        if self.ptype == 2 || total > n * cap {
            self.ptype = 2;
            w.iter_mut().for_each(|x| *x = cap);
        } else {
            for x in w.iter_mut() {
                let v = *x + 100.0;
                *x = if cap < v { cap } else { v };
            }
        }
        if crate::dump::on("p68") {
            crate::dump::emit(
                "p68",
                serde_json::json!({"n": self.refactorizations, "it": self.iters, "type": self.ptype, "ucol": fwd, "lcol": bwd, "p68": w}),
            );
        }
        self.p68 = w;
    }

    /// Touch of the row weights at a factor update (`0x132a10`, from `0xf2eb3`): the rows of the pivot column grow by the pivot row's entry count,
    /// up to m / 2, and the pivot row is at least that count plus 100.
    fn touch_pricing(&mut self, r: usize, col: &[(usize, f64)]) {
        if self.ptype == 0 {
            return;
        }
        let delta = self.p60;
        let cap = self.m as f64 * 0.5;
        for &(i, _) in col {
            self.p68[i] = minsd(cap, self.p68[i] + delta);
        }
        self.p68[r] = maxsd(self.p68[r], 100.0 + delta);
    }

    /// `0xbdc10`: flags of all rows and the list of the infeasible ones sorted by priority (descending).
    fn build_walk_list(&mut self) {
        let ptol = self.ptol();
        let mut rows = Vec::new();
        let mut keys = Vec::new();
        for r in 0..self.m {
            let f = self.infeas(r);
            if f > ptol {
                self.flags[r] = 0;
                rows.push(r);
                keys.push(-self.pri(r, f));
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
            -s.pri(r, f)
        };
        let mut keys: Vec<f64> = fresh.iter().map(|&r| key(self, r)).collect();
        ref_sort(&mut keys, &mut fresh);
        let old: Vec<usize> = self.list[..self.walk_len.min(self.list.len())]
            .iter()
            .copied()
            .filter(|&r| self.flags[r] == 0)
            .collect();
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
        // The work of the rebuild (`0xc26ec`: the rows by class, `0xc2761`: the waiting entries, `0xc2842`: the new ones, `0xc2ad4`: the merged list).
        let (m, n1) = (self.m as f64, merged.len() as f64);
        let old_len = self.walk_len.min(self.list.len());
        // The loop over the rows (`0xc22f3` / `0xc24b7`) counts every row as feasible when `ctx+0x2c8` is 0 (the largest infeasibility sum at the last recompute is
        // at most 1e20, i.e. the tolerance level is 0), so that the site books `3 m`.
        let n_inf = if self.level != 0 { n1 } else { 0.0 };
        self.work.add(
            0xc26ec,
            3.0 * (m - n_inf) + (5.0 + if self.ptype == 1 { 1.0 } else { 0.0 }) * n_inf,
            crate::work::CTX,
        );
        self.work.add(
            0xc2761,
            2.0 * old.len() as f64 + old_len as f64,
            crate::work::CTX,
        );
        self.work
            .add(0xc2842, 2.0 * fresh.len() as f64, crate::work::CTX);
        self.work.add(0xc2ad4, n1, crate::work::CTX);
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
        self.rebuilds = 40;
        let mut maxp = 0.0f64;
        let mut inf = vec![false; self.m];
        for r in 0..self.m {
            let f = self.infeas(r);
            if f > ptol {
                inf[r] = true;
                self.prio[r] = self.pri(r, f);
                maxp = maxsd(maxp, self.prio[r]);
            }
        }
        if self.first_recompute {
            self.first_recompute = false;
            if !self.warm {
                self.build_walk_list();
                return;
            }
        }
        // The mode switch (`0xc0c29`): the value mode starts when 15 times the average of the changed-row counts exceeds m, or when the iteration
        // count does; until then the walk list is rebuilt.
        if !self.value_mode
            && !self.warm
            && !(15.0 * self.avg_touch > self.m as f64 || self.iters as f64 > self.m as f64)
        {
            self.build_walk_list();
            return;
        }
        if !self.value_mode {
            // The value mode lowers the factor's update limit to max(20, m^0.65) (`0xc0cad`).
            let lim = ((self.m as f64).powf(0.65) as usize).max(20);
            let cur = self
                .mirror
                .upd_limit
                .unwrap_or(((self.m as f64).powf(0.7) as usize).max(20));
            self.upd_limit = Some(lim.min(cur));
            self.mirror.upd_limit = self.upd_limit;
        }
        self.value_mode = true;
        self.cursor = 0;
        // New threshold (`0xc1e92..0xc1ee9`): max priority / 1e5, or / 1e3 when the level was zero before this factorization (no huge infeasibilities)
        // and the list is not shorter than the pivot column's entry count (always, after a no-row call); an existing threshold below the maximum
        // priority is only ever raised.
        let thr_old = self.thr;
        // The count compared with the pivot column's entries is the number of infeasible rows (`[P+0x3c]`, `0xc1e8e`).
        let count = inf.iter().filter(|&&f| f).count();
        self.n_infeas = count;
        let div = if after_pivot && count < self.last_cap || big {
            1e5
        } else {
            1e3
        };
        let new = maxp / div;
        self.thr = if thr_old != 0.0 && maxp > thr_old {
            maxsd(new, thr_old)
        } else {
            new
        };
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
                // The list is exhausted: rebuild it (`0xbdc10`), the first entry is the argmax. When the rebuild budget `P+0x38` (40 after every recompute) runs
                // out, there is no row and a refactorization is requested instead (`0xcce9f..0xcceb5`, request 3; 80BAU3B call 1336).
                self.rebuilds -= 1;
                self.cursor = 0;
                if self.rebuilds < 0 {
                    self.request = 3;
                    return None;
                }
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
        let (walked, mut removed) = (self.list.len(), 0usize);
        for &e in &self.list {
            if self.flags[e] >= 2 {
                if self.flags[e] < 5 {
                    self.flags[e] = 2;
                }
                removed += 1;
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
        // The walk's work (`0xcce72`): 0.5 per entry kept, 1 per entry that left as feasible, 2 per other entry walked.
        if walked > 0 {
            let est = 0.5 * kept.len() as f64 + removed as f64 + 2.0 * (walked - removed) as f64;
            self.work.add(0xcce72, est, crate::work::CTX);
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
                self.prio[r] = self.pri(r, f);
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
    fn pivot_row(
        &self,
        rho: &[(usize, f64)],
        leave: usize,
        above: bool,
    ) -> (Vec<(usize, f64)>, usize) {
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

    /// The second FTRAN vector of the pricing mode 3 (`0xc57d8..0xc58c4`, `0x11faf9..0x11fcd0`): the entries of the pivot row whose variable was basic when the
    /// weights restarted, in the true frame (the leaving variable with 1), combined column by column into a dense vector in the order of the list; the
    /// entries above 1e-13 in magnitude, in the order of their first touch, are the vector. Returns it with the squared norm of the list and the work of
    /// the wrapper without the column copy (`0x120114`: 4 per list entry, 2 per column entry of the structurals, 2 per touched row, 2 per kept row).
    fn framework_rhs(
        &self,
        row: &[(usize, f64)],
        leave: usize,
        above: bool,
    ) -> (Vec<usize>, Vec<f64>, f64, f64) {
        let p = self.p;
        let sdir = if above { -1.0 } else { 1.0 };
        let mut list: Vec<(usize, f64)> = Vec::new();
        for &(j, a) in row {
            if self.dflag[j] {
                continue;
            }
            let a = if a.abs() < ALPHA_TOL { 0.0 } else { a };
            let s = if self.st[j] == UB { -1.0 } else { 1.0 };
            list.push((j, a * (s * sdir)));
        }
        if leave < self.n && !self.dflag[leave] {
            list.push((leave, 1.0));
        }
        let w = rho_norm(&list);
        let mut dense = vec![0.0f64; self.m];
        let mut seen = vec![false; self.m];
        let mut order: Vec<usize> = Vec::new();
        let mut work = 4.0 * list.len() as f64;
        for &(j, v) in &list {
            if j >= self.n {
                let i = j - self.n;
                if !seen[i] {
                    seen[i] = true;
                    order.push(i);
                }
                dense[i] += p.fcols[j][0].1 * v;
            } else {
                work += 2.0 * p.fcols[j].len() as f64;
                for &(i, a) in &p.fcols[j] {
                    if !seen[i] {
                        seen[i] = true;
                        order.push(i);
                    }
                    dense[i] = a * v + dense[i];
                }
            }
        }
        let mut rows = Vec::new();
        let mut vals = Vec::new();
        let order_len = order.len();
        for i in order {
            if dense[i].abs() > 1e-13 {
                rows.push(i);
                vals.push(dense[i]);
            }
        }
        work += 2.0 * order_len as f64 + 2.0 * rows.len() as f64;
        (rows, vals, w, work)
    }

    /// Work of the row product (`0x134d90`, site `0x135e92`, same estimate as the primal path's `row_product_work`).
    fn row_product_work(&mut self, rho: &[(usize, f64)], nslack: usize) {
        let p5: usize = rho.iter().map(|e| self.rrows[e.0].len()).sum();
        let (n, m) = (self.n, self.m);
        let mut seen = vec![false; n + m];
        let mut d = 0usize;
        for &(i, _) in rho {
            for &j in &self.rrows[i] {
                if !seen[j] {
                    seen[j] = true;
                    d += 1;
                }
            }
        }
        let nr = rho.len() as f64;
        let est = if rho.len() == 1 {
            if rho[0].1.abs() < 1.0 {
                10.5 + 5.0 * p5 as f64
            } else {
                10.5 + 4.0 * p5 as f64
            }
        } else if (n as f64) < 1.8 * p5 as f64 {
            // Columnwise in index order (the same test as `pivot_row`).
            8.5 * nr + 2.0 * p5 as f64 + 2.5 * n as f64
        } else {
            11.5 * nr + 3.5 * p5 as f64 + 3.0 * d as f64
        };
        // Every slack entry of the result (the leaving one and the offered nonbasic ones) adds four.
        let slack = if self.level != 0 {
            4.0 * nslack as f64
        } else {
            0.0
        };
        self.work
            .add_from(0x135e92, 0xc5381, est + slack, crate::work::CTX);
    }

    /// The bound-flipping Harris ratio test (`0xc74a0`, `ctx+0x468 != 0`), `specs/dual-simplex.md`.
    fn ratio_test(
        &mut self,
        row: &[(usize, f64)],
        esi: usize,
        leave: usize,
        above: bool,
        xbr: f64,
    ) -> Ratio {
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
        self.work.add(0xc7bcd, 3.0 * cnt as f64, crate::work::CTX);
        if nn == 0 {
            self.work.add(0xc7c4f, 0.0, crate::work::CTX);
            return Ratio {
                q: None,
                flips: Vec::new(),
            };
        }
        // The partition passes of the selection book `3 (hi - lo + 4)` each; the store is separate (`0xcbcae`) when the final sum is below minus the primal tolerance (`xmm9 = -1e19` at a nonzero level, else `-1e-6`, against `xmm10` at `0xcbc8d`).
        let mut sel = 0.0f64;
        let mut sel_sep = false;
        if nn == 1 {
            self.work.add(0xc7c4f, 5.0, crate::work::CTX);
            let (q, _) = self.shift_entering(row, &[lst[0]], lst[0], 0.0, 0.0, 0);
            return Ratio {
                q: Some(q),
                flips: Vec::new(),
            };
        }
        let maxa = -minsd(xm13, xm9);
        let mut nflip = 0usize;
        if bestc + xmm10 < 0.0 {
            let mut s = xmm10;
            let (mut lo, mut hi) = (0usize, nn - 1);
            while lo < hi {
                sel += (hi - lo + 4) as f64;
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
                    while ratio_arr[l as usize] > pr
                        || (ratio_arr[l as usize] == pr && pc > con[l as usize])
                    {
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
                    s = accumulate(s, &con[lo..=h], lo % 2, false);
                } else {
                    s = accumulate(s, &con[h + 1..=hi], (h + 1) % 2, true);
                }
                if 0.0 > s {
                    lo = h + 1;
                } else {
                    hi = h;
                }
            }
            nflip = lo;
            sel_sep = -self.ptol() > s;
        } else if best_i > 0 {
            let b = best_i as usize;
            ratio_arr.swap(0, b);
            con.swap(0, b);
            lst.swap(0, b);
        }
        let w5 = 5.0 * nn as f64;
        if sel == 0.0 {
            self.work.add(0xc7c4f, w5, crate::work::CTX);
        } else if sel_sep {
            self.work.add(0xcbcae, 3.0 * sel, crate::work::CTX);
            self.work.add(0xc7c4f, w5, crate::work::CTX);
        } else {
            self.work
                .add_merged(0xc7c4f, &[3.0 * sel, w5], crate::work::CTX);
        }
        if nflip > 0 && !(self.d[row[lst[nflip]].0] > delta) {
            nflip = 0;
        }
        let pivtol = maxsd(0.0, maxsd(1e-13 * maxa, PIVOT_FLOOR));
        if !(pivtol <= maxa) {
            return Ratio {
                q: None,
                flips: Vec::new(),
            };
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
        self.work.add(
            0xc82ee,
            3.0 * (nn - b) as f64 + kept.len() as f64,
            crate::work::CTX,
        );
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
            return Ratio {
                q: None,
                flips: Vec::new(),
            };
        }
        // The entries of the list that the second loop `0xc8a05` walks: the flips and the kept candidates (`r8d` at `0xc8b2e`).
        self.work
            .add(0xc8b67, 4.0 * (b + kept.len()) as f64, crate::work::CTX);
        let (chosen, xm5) = self.select_entering(row, &kept, xm3);
        let (mut chosen, xm5) = self.shift_entering(row, &kept, chosen, xm5, xm15, b);
        let mut flips = lst[..b].to_vec();
        if self.flip_select && !flips.is_empty() {
            (chosen, flips) = self.select_among_flips(row, chosen, flips, xm5);
        }
        Ratio {
            q: Some(chosen),
            flips,
        }
    }
}

impl<'a> Dual<'a> {
    /// Choice of the entering variable among the Harris-filtered candidates (`0xc89b3`): the smallest `a / w`, `w` the column count of a structural, 0.5 of a slack.
    fn select_entering(&self, row: &[(usize, f64)], kept: &[usize], xm3: f64) -> (usize, f64) {
        let mut chosen = kept[0];
        let mut xm5 = 0.0f64;
        for &e in kept {
            let (j, av) = row[e];
            let wgt = if j < self.n {
                self.p.cols[j].len() as f64
            } else {
                0.5
            };
            if xm3 * av < self.d[j] {
                continue;
            }
            if xm5 * wgt > av {
                xm5 = av / wgt;
                chosen = e;
            }
        }
        (chosen, xm5)
    }

    /// The next value of the reference's random generator (`0x88a870`: 31 bits; `0x88a8d0`: that times 2^-31).
    fn rng_int(&mut self) -> u32 {
        lcg_int(&mut self.rng)
    }

    fn rng_double(&mut self) -> f64 {
        self.rng_int() as f64 * 4.656612873077393e-10
    }

    fn raise_pert(&mut self, to: u8) {
        self.pert = if to == 1 {
            if self.pert == 0 {
                1
            } else {
                self.pert
            }
        } else {
            self.pert.max(2)
        };
    }

    /// Cost shift of a nonbasic variable at a bound whose reduced cost has the wrong sign (`d' < -1e-13`; `0xcac52..0xcaeac`, `0xca1a0..0xca4c6`): with
    /// `m = -d'`, `d'` grows by `m (1 + u)` with `u` uniform in [0, 1) when `m` is tiny (below 5e-4 delta, or below 0.05 delta with `0.1 delta` above the
    /// range times `m`), else by `m ((r mod 500 + 111) / 1000 + 1)` with `r` an integer draw, so that `d'` becomes positive (0.1 to 0.6 of `m`). The cost
    /// takes the shift (plus at the lower bound, minus at the upper), the objective offset the shift times the bound.
    fn shift_cost(&mut self, j: usize) {
        let delta = DUAL_TOL;
        let m = -self.d[j];
        let (lo, up) = (self.p.lo[j], self.p.up[j]);
        let small = delta * 0.0005 > m || (delta * 0.05 > m && delta * 0.1 > (up - lo) * m);
        let s = if small {
            self.raise_pert(1);
            let u = self.rng_double();
            m * (u + 1.0)
        } else {
            self.raise_pert(2);
            let r = self.rng_int() as i32;
            m * ((((r % 500) + 111) as f64) / 1000.0 + 1.0)
        };
        self.d[j] += s;
        if self.st[j] == LB {
            self.cshift[j] += s;
            self.cost[j] = self.pcost[j] + self.cshift[j];
            self.obj_off -= lo * s;
            if -BIG >= lo {
                self.raise_pert(2);
            }
        } else {
            self.cshift[j] -= s;
            self.cost[j] = self.pcost[j] + self.cshift[j];
            self.obj_off += up * s;
            if up >= BIG {
                self.raise_pert(2);
            }
        }
    }

    /// Treatment of the chosen entering variable's reduced cost (`0xcabc2`, `0xc8c12..0xc8e27`). A variable at a bound with `d' < -1e-13` is shifted; the other
    /// candidates of the list with a negative `d'` are then shifted or zeroed, the Harris filter runs again on the new values (continuing from the bound
    /// `xm15` reached, with `xm3 = xm15`) and the entering variable is chosen again. A `d'` in [-1e-13, 0) is zeroed; a free variable's is zeroed.
    fn shift_entering(
        &mut self,
        row: &[(usize, f64)],
        kept: &[usize],
        chosen: usize,
        xm5: f64,
        mut xm15: f64,
        nflips: usize,
    ) -> (usize, f64) {
        let mut xm3;
        let mut xm5 = xm5;
        let delta = DUAL_TOL;
        let mut chosen = chosen;
        let mut kept: Vec<usize> = kept.to_vec();
        for _ in 0..8 {
            let j = row[chosen].0;
            let dj = self.d[j];
            if self.st[j] == LB || self.st[j] == UB {
                if dj < -1e-13 {
                    self.shift_cost(j);
                    if self.d[j] == 0.0 {
                        break;
                    }
                } else {
                    if dj < 0.0 {
                        self.d[j] = 0.0;
                    }
                    break;
                }
            } else if self.st[j] == FREE {
                if dj.abs() > 1e-13 {
                    self.raise_pert(if delta * 0.0005 > dj.abs() { 1 } else { 2 });
                }
                self.d[j] = 0.0;
                break;
            } else {
                break;
            }
            xm3 = xm15;
            let mut next: Vec<usize> = Vec::new();
            for &e in &kept {
                let (jj, av) = row[e];
                if self.d[jj] < 0.0 {
                    if self.st[jj] == LB || self.st[jj] == UB {
                        if self.d[jj] < -1e-13 {
                            self.shift_cost(jj);
                        } else {
                            self.d[jj] = 0.0;
                        }
                    } else if self.st[jj] == FREE {
                        if self.d[jj].abs() > 1e-13 {
                            self.raise_pert(if delta * 0.0005 > self.d[jj].abs() {
                                1
                            } else {
                                2
                            });
                        }
                        self.d[jj] = 0.0;
                    }
                }
                let dj2 = self.d[jj];
                if xm15 * av > dj2 {
                    let r = (delta + dj2) / av;
                    xm3 = maxsd(r, xm3);
                    next.push(e);
                    if dj2 >= 0.0 {
                        xm15 = maxsd(r, xm15);
                    }
                }
            }
            if next.is_empty() {
                chosen = kept[0];
                xm5 = 0.0;
            } else {
                // The selection loop runs again over the compacted list and books its count (`0xc8b67`) like the first one.
                self.work.add(
                    0xc8b67,
                    4.0 * (nflips + next.len()) as f64,
                    crate::work::CTX,
                );
                (chosen, xm5) = self.select_entering(row, &next, xm3);
                kept = next;
            }
        }
        (chosen, xm5)
    }

    /// Second selection among the bound-flip candidates once `ctx+0x3fc` is set (`0xc8ecf..0xc9439`): with `theta = d_q / a_q` of the chosen variable, a flipped
    /// entry with `0 <= d <= theta a - 1000 delta` has the weight `(count + ((theta a - d) / delta)^0.25) (0.3 f + 0.7)` (count 0.5 for a slack, factor `f` 1
    /// for continuous variables; the power term while `ctx+0x178 < 10`) and replaces the chosen one when `xm5 w > a` (`xm5` the best `a / w` so far, that of
    /// the first selection at the start). The new entering variable keeps as flips the old flips (other than itself) with `d <= a theta'`.
    fn select_among_flips(
        &mut self,
        row: &[(usize, f64)],
        chosen: usize,
        flips: Vec<usize>,
        xm5: f64,
    ) -> (usize, Vec<usize>) {
        let delta = DUAL_TOL;
        let (qj, aq) = row[chosen];
        let theta = self.d[qj] / aq;
        let tol = 1000.0 * delta;
        let mut xm5 = xm5;
        let mut best: Option<usize> = None;
        // The sites `0xc92b8` and `0xc9431` add to another accumulator than the driver's W (no write to W in any capture), so nothing is booked here.
        for (k, &e) in flips.iter().enumerate() {
            let (j, aj) = row[e];
            let dj = self.d[j];
            if 0.0 > dj {
                continue;
            }
            let t = theta * aj;
            if dj > t - tol {
                continue;
            }
            let pw = ((t - dj) / delta).powf(0.25);
            let w0 = if j < self.n {
                self.p.cols[j].len() as f64 + pw
            } else {
                pw + 0.5
            };
            let w = w0 * (1.0 * 0.3 + 0.7);
            if xm5 * w > aj {
                xm5 = aj / w;
                best = Some(k);
            }
        }
        let Some(k) = best else {
            return (chosen, flips);
        };
        let qe = flips[k];
        let (nj, na) = row[qe];
        let theta2 = self.d[nj] / na;
        let kept: Vec<usize> = flips
            .iter()
            .copied()
            .filter(|&e| e != qe && self.d[row[e].0] <= row[e].1 * theta2)
            .collect();
        (qe, kept)
    }
}

/// The running sum of the ratio test's selection (`0xcba58..0xcbc3a`): `s` plus (or minus) the entries of `c` in the reference's vectorized order. Fewer
/// than `peel + 8` entries are added one by one. Otherwise `peel` entries (one when the first is at an odd index of the 16-byte aligned array, so that the
/// pairs are aligned) go in one by one, then blocks of eight enter four two-lane accumulators (lane 0 of the first one starts at `s`), which are combined
/// as `(a0 + a1) + (a2 + a3)` and then lane 0 + lane 1, and the remaining entries follow one by one.
fn accumulate(s: f64, c: &[f64], peel: usize, sub: bool) -> f64 {
    let step = |a: f64, x: f64| if sub { a - x } else { a + x };
    let n = c.len();
    let mut s = s;
    if n < peel + 8 {
        for &x in c {
            s = step(s, x);
        }
        return s;
    }
    for &x in &c[..peel] {
        s = step(s, x);
    }
    let vend = n - ((n - peel) & 7);
    let mut acc = [[0.0f64; 2]; 4];
    acc[0][0] = s;
    let mut i = peel;
    while i < vend {
        for (q, a) in acc.iter_mut().enumerate() {
            a[0] = step(a[0], c[i + 2 * q]);
            a[1] = step(a[1], c[i + 2 * q + 1]);
        }
        i += 8;
    }
    let l0 = [acc[0][0] + acc[1][0], acc[0][1] + acc[1][1]];
    let l2 = [acc[2][0] + acc[3][0], acc[2][1] + acc[3][1]];
    let l = [l0[0] + l2[0], l0[1] + l2[1]];
    s = l[0] + l[1];
    for &x in &c[vend..] {
        s = step(s, x);
    }
    s
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
    /// The ratio test found no entering variable: the dual is unbounded (the primal infeasible). Carries the BTRAN row `e_r' B^-1` (dense, by
    /// row), which is the Farkas certificate.
    Unbounded(Vec<f64>),
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
    above: bool,
}

impl<'a> Dual<'a> {
    /// One call of the iteration function `0xc4fe0`.
    fn iteration(&mut self, k: usize, dump: &mut Option<serde_json::Value>) -> Step {
        let p = self.p;
        let (n, m) = (self.n, self.m);
        let snap = if crate::dump::on_from("dsx", k) {
            Some(self.snapshot(k))
        } else if crate::dump::on("dsxl") {
            Some(serde_json::json!({"k": k}))
        } else {
            None
        };
        self.work.trig();
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
        let rho: Vec<(usize, f64)> = rho_full
            .iter()
            .copied()
            .filter(|e| e.1.abs() > DROP_TOL)
            .collect();
        self.p60 = rho.len() as f64;
        self.work
            .add_from(0x13174d, 0xc51dd, self.mirror.work, crate::work::FAC);
        let (row, esi) = self.pivot_row(&rho, leave, above);
        self.row_product_work(&rho, row.len() - esi);
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
            let mut y = vec![0.0; m];
            for &(i, v) in &rho_full {
                y[i] = if above { -v } else { v };
            }
            // Roundoff in a row with an infinite bound would void the certificate: entries at noise level are dropped.
            let big = y.iter().fold(0.0f64, |a, v| a.max(v.abs()));
            y.iter_mut()
                .filter(|v| v.abs() <= 1e-11 * big)
                .for_each(|v| *v = 0.0);
            return Step::Unbounded(y);
        };
        let (q, aq_row) = row[qe];
        let flips: Vec<usize> = rt.flips.iter().map(|&e| row[e].0).collect();
        if !flips.is_empty() {
            self.work
                .add(0xc561a, 7.0 * flips.len() as f64, crate::work::CTX);
        }
        if let Some(s) = rec.as_mut() {
            s["q"] = serde_json::json!(q);
            s["alpha_row"] = serde_json::json!(aq_row);
            s["flips"] = serde_json::json!(flips);
        }
        // The reference's two-vector FTRAN: the pivot column first, then tau = B^-1 rho at the column's positions only.
        let rows_q: Vec<usize> = p.fcols[q].iter().map(|&(i, _)| i).collect();
        let vals_q: Vec<f64> = p.fcols[q].iter().map(|&(_, a)| a).collect();
        let col_out = self.mirror.ftran(&rows_q, &vals_q);
        // Mode 3: the second vector is not rho but the combination of the columns of the row's entries that were basic when the weights restarted
        // (`0xc57d8`, `0x11f620`), and its squared norm is the weight of the pivot row.
        let (rho_rows, rho_vals, w_row, vec_work) = if self.pmode == 3 {
            self.framework_rhs(&row, leave, above)
        } else {
            (
                rho.iter().map(|e| e.0).collect::<Vec<usize>>(),
                rho.iter().map(|e| e.1).collect::<Vec<f64>>(),
                rho_norm(&rho),
                rho.len() as f64,
            )
        };
        let needed: Vec<usize> = col_out.iter().map(|e| e.0).collect();
        // Mode 3: the framework vector is a dense vector (count -1) once it has more than m/8 entries (SC50A calls 11, 14..; AFIRO call 12: 4 of m = 27;
        // 6 of 50 stay sparse); its L stage then runs dense.
        let dense_in = self.pmode == 3 && rho_rows.len() * 8 > m;
        let tau_out = self
            .mirror
            .ftran_needed(&rho_rows, &rho_vals, &needed, dense_in);
        {
            // The two-vector FTRAN books the L stages of both vectors and then the U stage of the column at the factor's scale
            // (`0x129671`, call `0x120200`); the wrapper `0x11f620` books `4 nnz(column) + nnz(rho)` (`0x120114`).
            let (e1, u1) = self.mirror.work_parts;
            let est = match self.mirror.joint_u(&needed) {
                Some(j) => (e1 + self.mirror.tau_l.get()) + j,
                None => ((e1 + self.mirror.tau_l.get()) + u1) + self.mirror.tau_u.get(),
            };
            self.work
                .add_from(0x129671, 0x120200, est, crate::work::FAC);
            // A slack column (`ctx+0x2fc >= n`, `0x11f788`) is not copied from the matrix: its `4 nnz` term is absent.
            let copied = if q < self.n { rows_q.len() } else { 0 };
            self.work
                .add(0x120114, 4.0 * copied as f64 + vec_work, crate::work::CTX);
        }
        let mut tau_dense = vec![0.0; m];
        for &(i, v) in &tau_out {
            tau_dense[i] = v;
        }
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
            s["w"] = serde_json::json!(w_row);
            s["col"] = crate::dump::sparse(&col_out);
            s["tau"] = crate::dump::sparse(&tau_out);
        }
        // Column without the pivot entry: the last entry moves into the pivot's slot.
        // A dense FTRAN result (count -1) is compacted by an ascending scan that skips the pivot and the entries up to the drop tolerance (`0x12013c`).
        let mut col: Vec<(usize, f64)> = col_out.clone();
        if self.mirror.last_dense() {
            col.retain(|e| e.0 != r && e.1.abs() > DROP_TOL);
        } else if let Some(k) = col.iter().position(|e| e.0 == r) {
            col.swap_remove(k);
        }
        let tau: Vec<f64> = col.iter().map(|e| tau_dense[e.0]).collect();
        // Bound flips: x_B and the statuses (`0xcd720`).
        if !flips.is_empty() {
            let mut rhs = vec![0.0; m];
            let mut touched: Vec<usize> = Vec::new();
            let mut est = 0usize;
            for &j in &flips {
                let delta = if self.st[j] == LB {
                    p.up[j] - p.lo[j]
                } else {
                    p.lo[j] - p.up[j]
                };
                est += p.cols[j].len();
                for &(i, a) in &p.fcols[j] {
                    if !touched.contains(&i) {
                        touched.push(i);
                    }
                    rhs[i] += a * delta;
                }
            }
            // The reference's input vector (`0xcd720`): a dense one when twice the summed column counts exceed m, else the touched
            // rows in order of first touch. Then the plain FTRAN of the factor.
            let dense_in = 2 * est > m;
            // Work of `0xcd720`: the input vector (`4 sum len` for a list, `2 m` for a dense one, doubled) before the FTRAN.
            self.work.add(
                0xcdbb5,
                2.0 * if dense_in {
                    2.0 * m as f64
                } else {
                    4.0 * est as f64
                },
                crate::work::CTX,
            );
            let rows: Vec<usize> = if dense_in {
                (0..m).filter(|&i| rhs[i] != 0.0).collect()
            } else {
                touched
            };
            let vals: Vec<f64> = rows.iter().map(|&i| rhs[i]).collect();
            let dense_est = if dense_in {
                Some(self.mirror.ftran_dense_work(&rhs))
            } else {
                None
            };
            let sol = self.mirror.ftran_keep(&rows, &vals, dense_in);
            self.work.add_from(
                0x1219d6,
                0x120673,
                dense_est.unwrap_or(self.mirror.keep_work),
                crate::work::FAC,
            );
            if dense_in || self.mirror.keep_dense {
                self.work
                    .add(0xce316, (m + sol.len()) as f64, crate::work::CTX);
            }
            // The count of the result vector times 7 on the context's scale (`(7 * scale) * count`).
            self.work
                .add(0xce1ef, sol.len() as f64, 7.0 * crate::work::FAC);
            // The x_B update and the touching follow the result list of the FTRAN in its own order (FIT2P call 44: the rows join the pricing list in that order).
            let mut changed = Vec::new();
            for (i, v) in sol {
                if v != 0.0 {
                    self.xb[i] -= v;
                    changed.push(i);
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
        let new_st_leave = if p.lo[leave] == p.up[leave] {
            FIXED
        } else if above {
            UB
        } else {
            LB
        };
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
        if self.d[q] <= DUAL_TOL {
            self.stall += 1;
        }
        self.iters += 1;
        *dump = rec;
        let _ = flips.len();
        // The factor's update is made by the post-iteration pass (it refactorizes instead when one is due).
        let alpha_vec = alpha_dense;
        self.pending_alpha = Some(alpha_vec);
        // The ratio test leaves the entries of the row list below 1e-10 in magnitude as exact zeros (RIN against ROUT of `0xc74a0`: 41376 entries on PEROLD,
        // largest 9.9997e-11, smallest kept 1.0000751e-10), so that the d update does not see them.
        let row: Vec<(usize, f64)> = row
            .into_iter()
            .map(|(j, v)| (j, if v.abs() < ALPHA_TOL { 0.0 } else { v }))
            .collect();
        Step::Pivot(Box::new(Pivot {
            r,
            leave,
            q,
            row,
            flips,
            col,
            tau,
            alpha,
            w: w_row,
            t0,
            above,
        }))
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
            // `ctx+0x3f8 != 0` in all captures: confirm on a fresh factor (`0x1aa280`). The cost shifts stay in force through the confirming refactorization
            // (`ctx+0x290` is 0, `0xe4891` is not taken) and are removed at the end of the solve.
            self.request = 4;
        } else if self.free_infeasible() {
            // A free nonbasic variable (never pivoted in: no basic variable with a finite lower bound) with a nonzero d.
            self.dual_unbounded();
        } else if self.level == 0 {
            self.status = Some(Status::Optimal);
        } else {
            let p = self.p;
            for j in 0..self.n + self.m {
                let bad = (self.st[j] == UB && p.up[j] >= BIG && self.d[j] > DUAL_TOL)
                    || (self.st[j] == LB && p.lo[j] <= -BIG && self.d[j] > DUAL_TOL);
                if bad {
                    self.dual_unbounded();
                    return;
                }
            }
            // No artificial bound is required by its d (`ctx+0x290 = 1`, read in the reference but not observed): a refactorization puts those variables
            // on their finite bound (`set_level`) and recomputes x_B without the sentinels.
            self.request = 4;
        }
    }

    /// Nonbasic at an artificial (infinite, 1e30) bound: at the upper bound of a variable without one, or at the lower bound of a variable without one.
    /// The recovered code of the reference counts the first kind only (`0xe5d5d`); the second occurs only with a lower bound of -1e30 and a finite upper
    /// bound, which no corpus model has.
    fn at_artificial(&self, j: usize) -> bool {
        (self.st[j] == UB && self.p.up[j] >= BIG) || (self.st[j] == LB && self.p.lo[j] <= -BIG)
    }

    /// A free nonbasic variable has a reduced cost beyond the dual tolerance.
    fn free_infeasible(&self) -> bool {
        (0..self.n + self.m).any(|j| self.st[j] == FREE && self.d[j].abs() > DUAL_TOL)
    }

    /// The dual is infeasible (unbounded dual ray in the primal's terms): a nonbasic variable sits at an artificial bound with the wrong-signed d, or is free
    /// with a nonzero d. The basic solution is primal feasible at the scale of the artificial bounds, so the point it describes is a finite point plus 1e30 times
    /// a recession direction. That direction moves every nonbasic variable at an artificial bound away from its finite bound, every free one with a nonzero d
    /// against its d, and the basic variables with them (the factor is fresh here). Whether the primal is feasible is for the caller to settle.
    fn dual_unbounded(&mut self) {
        let p = self.p;
        let mut rhs = vec![0.0; self.m];
        let mut ray = vec![0.0; self.n + self.m];
        for k in 0..self.n + self.m {
            let dir = if self.st[k] == UB && p.up[k] >= BIG {
                1.0
            } else if self.st[k] == LB && p.lo[k] <= -BIG {
                -1.0
            } else if self.st[k] == FREE && self.d[k].abs() > DUAL_TOL {
                -self.d[k].signum()
            } else {
                continue;
            };
            ray[k] = dir;
            for &(i, a) in &p.cols[k] {
                rhs[i] -= a * dir;
            }
        }
        self.factor.ftran(&mut rhs);
        for (pos, &b) in self.head.iter().enumerate() {
            ray[b] = rhs[pos];
        }
        ray.truncate(self.nm);
        self.ray = Some(ray);
        self.status = Some(Status::Unbounded);
    }

    /// Work of `0xe45c0` after a refactorization (site `0xe7000`): two copies of x_B and d's rows (`2 m`), the status pass (`m + 2 n`), a pass over all variables
    /// when the level was not zero, and the tail that is longer for the confirming refactorization (request 4 or 6). With the perturbation flag set (`ctx+0x338`
    /// != 0, kept to the end of the solve) the request-4 tail is the plain `1.5 n` and every such factorization adds `2 m + n` (the cost recompute with the shift
    /// array, `0xe555c`; ADLITTLE, BOEING2, PILOT4).
    fn set_level_work(&mut self, level_was_set: bool) {
        let (n, m) = (self.n as f64, self.m as f64);
        let mut est = 2.0 * m + (m + 2.0 * n);
        if level_was_set {
            est += 0.5 * (n + m);
        }
        est += if (self.request == 4 || self.request == 6) && self.pert == 0 {
            2.0 * n + (n + m)
        } else {
            1.5 * n
        };
        if self.pert != 0 {
            est += 2.0 * m + n;
        }
        self.work.add(0xe7000, est, crate::work::CTX);
    }

    /// Work of the pricing setup after a refactorization: the list of infeasible rows (`5 m + 4 len`, len the walk list or, in the value mode, the number of infeasible rows `P+0x3c`; `0xc4329` in the walk mode, `0xc1fe3` in the value mode)
    /// and the closing pass over all variables.
    fn flags_work(&mut self) {
        let (n, m) = (self.n as f64, self.m as f64);
        let site = if self.value_mode { 0xc1fe3 } else { 0xc4329 };
        let listed = if self.value_mode {
            self.n_infeas
        } else {
            self.list.len()
        };
        self.work
            .add(site, 5.0 * m + 4.0 * listed as f64, crate::work::CTX);
        // The closing pass (`0xbf5b0`): the variant for a nonzero tolerance level books `bfed2`, the other `bf981`; the perturbation flag adds one unit per variable.
        let pert = if self.pert != 0 { 1.0 } else { 0.0 };
        self.work.add(
            if self.level != 0 { 0xbfed2 } else { 0xbf981 },
            (2.0 + pert) * (n + m),
            crate::work::CTX,
        );
    }

    /// Level of the tolerance (`0xe45c0`, after a refactorization): with request 4 the nonbasic variables are first put on the dual-feasible bound.
    fn set_level(&mut self, old_level: usize) {
        let p = self.p;
        if self.request == 4 {
            for j in 0..self.n + self.m {
                let (lo_inf, up_inf) = (p.lo[j] <= -BIG, p.up[j] >= BIG);
                if (self.st[j] == LB || self.st[j] == UB)
                    && lo_inf
                    && up_inf
                    && self.d[j].abs() <= DUAL_TOL
                {
                    // A free variable no longer tied to an artificial bound by its d goes back to the free status (value 0).
                    self.st[j] = FREE;
                    self.d[j] = 0.0;
                } else if self.st[j] == LB && self.d[j] < -DUAL_TOL {
                    self.st[j] = UB;
                    self.d[j] = -self.d[j];
                } else if self.st[j] == UB && up_inf && self.d[j] <= DUAL_TOL {
                    self.st[j] = LB;
                    self.d[j] = -self.d[j];
                } else if self.st[j] == LB && lo_inf && !up_inf && self.d[j] <= DUAL_TOL {
                    // The mirror of the previous rule for a variable without a lower bound.
                    self.st[j] = UB;
                    self.d[j] = -self.d[j];
                }
            }
        }
        self.level = (0..self.n + self.m)
            .filter(|&j| self.at_artificial(j))
            .count();
        self.stall_check(old_level);
        if !self.flip_select && (self.pert >= 3 || 3.0 * self.iters as f64 > self.m as f64) {
            self.flip_select = true;
        }
    }

    /// The dual objective as the reference keeps it (`ctx+0x118`): the working costs times the basic solution and the nonbasic values, plus the shift
    /// offset, divided by the objective scale.
    fn fresh_objective(&self) -> f64 {
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
        self.factor.ftran(&mut r);
        let mut s = 0.0;
        for (pos, &j) in self.head.iter().enumerate() {
            s += self.cost[j] * r[pos];
        }
        for j in 0..p.n + p.m {
            if self.st[j] < 0 {
                s += self.cost[j] * self.nonbasic_value(j);
            }
        }
        (s + self.obj_off) / self.obj_scale
    }

    /// The pricing switch after a factorization (`0xf770d..0xf82c5`, release profile only): at every third factorization, with the steepest-edge pricing
    /// still on, more than m + 1000 iterations done and fewer than 100 m structurals, the weights restart at 1 and the mode becomes 3 ("simplex pricing
    /// changes from 1 to 3"). The pricing update formulas stay those of the steepest edge.
    fn pricing_switch(&mut self) {
        if !self.switch_enabled || self.pmode != 1 || self.refactorizations % 3 != 0 {
            return;
        }
        let m = self.m as f64;
        if self.iters as f64 > m + 1000.0 && (self.n as f64) < 100.0 * m {
            self.pmode = 3;
            self.dflag = self.st.iter().map(|&t| t < 0).collect();
            for w in self.w.iter_mut() {
                *w = 1.0;
            }
        }
    }

    /// The stall test of the refactorization post pass (`0xe60cc..0xe6237`) and the cost perturbation it triggers (`0xe6af8..0xe6cd7`, `0xe74b1`).
    fn stall_check(&mut self, old_level: usize) {
        if self.pstep <= 0.0 && self.pert_param <= 0.0 {
            return;
        }
        // The test runs at every second refactorization (every fifth once the perturbation flag has reached 3).
        let period = if self.pert > 2 { 5 } else { 2 };
        if self.refactorizations % period != 1 {
            return;
        }
        let mut trig = false;
        if self.pstep > 0.0 && old_level <= self.level {
            // More than 90 percent of the iterations since the last check had an entering d of at most delta.
            trig =
                (self.stall - self.stall_snap) as f64 > 0.9 * (self.iters as f64 - self.stall_it);
            if !trig && self.prog_on {
                // No objective progress since the last check.
                let obj = self.fresh_objective();
                trig = 1e-13 * (obj.abs().min(1e10) + 1.0) + self.best_obj > obj;
                if !trig {
                    self.best_obj = maxsd(obj, self.best_obj);
                }
            }
        } else if self.prog_on {
            let obj = self.fresh_objective();
            self.best_obj = maxsd(obj, self.best_obj);
        }
        if trig && self.pert >= 2 && !(self.pert == 3 && self.pstep >= 0.01) {
            self.perturb();
        } else if trig && self.pert < 2 {
            if self.obj_scale >= 1.0 {
                self.rescale_objective();
            }
            self.perturb();
        }
        self.stall_snap = self.stall;
        self.stall_it = self.iters as f64;
    }

    /// The objective rescale of the stall trigger with the perturbation flag below 2 and `ctx+0xd8 >= 1` (`0xe6266..0xe6a88`, then the perturbation): the
    /// factor is a power of two from the largest `|d|` among the free nonbasic variables, `M` (`ctx+0x388`, the tolerance, is the reference value `a`):
    /// 128 when `a > 128 M`; else the smallest power of two `x < 128` with `a > M x` (or 128), halved when above 2, else 1 (no rescale). With a factor above 1
    /// the costs, the reduced costs and the objective scale `ctx+0xd8` are multiplied by it.
    fn rescale_objective(&mut self) {
        let a = DUAL_TOL;
        let mut mx = 0.0f64;
        for j in 0..self.n + self.m {
            if self.st[j] == FREE {
                mx = maxsd(self.d[j].abs(), mx);
            }
        }
        let f = if a > 128.0 * mx {
            128.0
        } else {
            let mut x = 1.0f64;
            loop {
                if a > mx * x {
                    break;
                }
                x += x;
                if !(128.0 > x) {
                    break;
                }
            }
            if 2.0 < x {
                x * 0.5
            } else {
                1.0
            }
        };
        if f <= 1.0 {
            return;
        }
        for v in self.cost.iter_mut() {
            *v *= f;
        }
        for v in self.pcost.iter_mut() {
            *v *= f;
        }
        for v in self.d.iter_mut() {
            *v *= f;
        }
        self.obj_scale *= f;
    }

    /// The cost perturbation (`0xe6af8..0xe6cd7`, `0xe74b1..0xe7585`): every nonbasic variable draws an integer `r` (also those that stay unchanged), and
    /// a variable at a bound moves its cost by `(r mod 250 + 1) (c348 / 500) log2(1 + count)` (count: the column count of a structural, 1 for a slack) in the
    /// direction that increases `d'` (the cost goes up at the lower bound and down at a finite upper bound; a variable at an infinite upper bound with
    /// `d' <= delta` is moved to its lower bound first). The scale grows threefold for the next perturbation.
    fn perturb(&mut self) {
        let p = self.p;
        let sc = self.pstep / 500.0;
        // The perturbation has its own generator, seeded with the iteration count (`0x115220`, state = it + the `Seed` parameter); the generator of the
        // cost shifts is not touched.
        let mut g = self.iters as u64;
        let mut acc = 0.0f64;
        for j in 0..self.n + self.m {
            if self.st[j] >= 0 {
                continue;
            }
            let r = lcg_int(&mut g) as i32;
            let count = if j < self.n {
                p.cols[j].len() as f64
            } else {
                1.0
            };
            let f = 1.0f64;
            let x = count * (f * 0.3 + 0.7);
            let xi = (((r % 250 + 1) as f64) * sc) * (x + 1.0).log2();
            if self.st[j] == LB {
                self.pcost[j] += xi;
                self.cost[j] += xi;
                self.d[j] += xi;
                acc -= p.lo[j] * xi;
            } else if self.st[j] == UB {
                if p.up[j] < BIG {
                    self.pcost[j] -= xi;
                    self.cost[j] -= xi;
                    self.d[j] += xi;
                    acc += p.up[j] * xi;
                } else if self.d[j] <= DUAL_TOL {
                    self.d[j] = -self.d[j];
                    self.pcost[j] += xi;
                    self.cost[j] += xi;
                    acc -= p.lo[j] * xi;
                    self.d[j] += xi;
                    self.st[j] = LB;
                }
            }
        }
        self.obj_off = acc + self.obj_off;
        self.prog_on = false;
        self.pert = 3;
        self.best_obj = -1e100;
        self.pert_left -= 1;
        let x = 3.0 * self.pstep;
        self.pstep = if self.pert_param == 0.0002 {
            minsd(0.01, x)
        } else {
            minsd(1000.0 * self.pert_param, x)
        };
    }

    /// The post-iteration passes `0xf6c82..0xf6d0a`: factor update or refactorization, the dual update, the primal and weight update.
    fn post(&mut self, pv: Option<&Pivot>, schedule_due: bool) -> Result<(), SolveError> {
        let refactor = self.request != 0 || schedule_due;
        if refactor {
            self.work.decided(self.work.w);
        }
        if let Some(pv) = pv {
            self.last_cap = pv.col.len();
            let alpha = self.pending_alpha.take().unwrap();
            if refactor {
                // The weights of the pivot are updated before the factorization renumbers the positions; the basic solution is recomputed.
                self.update_weights(pv);
                self.pivot_cols = Some(pv.col.len());
                self.refactor()?;
            } else {
                self.factor.update(pv.r, alpha);
                self.mirror.update(pv.r);
                for (site, est, scale) in std::mem::take(&mut self.mirror.upd_work) {
                    self.work.add(site, est, scale);
                }
                self.touch_pricing(pv.r, &pv.col);
                self.updates += 1;
            }
            // d update (`0xe45c0`).
            // The reference's row list also holds the leaving variable (rose's only when it is a slack).
            let aq = pv.row.iter().find(|e| e.0 == pv.q).unwrap().1;
            let theta = self.d[pv.q] / aq;
            // A zero step skips the loop (`0xe4719`): no work.
            if !refactor && theta != 0.0 {
                let row_cnt = pv.row.len() + usize::from(pv.leave < self.n);
                self.work
                    .add(0xe7000, 3.0 * row_cnt as f64, crate::work::CTX);
            }
            let old_leave = self.d[pv.leave];
            // The reference's row holds a slack leaving variable with its computed coefficient (about 1, 0.9999999999999997 in PILOT.JA call 30657), the step on it uses that.
            let lcoef = pv.row.iter().find(|e| e.0 == pv.leave).map_or(1.0, |e| {
                if pv.above {
                    -e.1
                } else {
                    e.1
                }
            });
            for &(j, a) in &pv.row {
                let a = if pv.flips.contains(&j) { -a } else { a };
                self.d[j] = self.d[j] - theta * a;
            }
            // A leaving variable that becomes fixed keeps the true reduced cost, not the frame value.
            // The leaving variable is in the reference's row with the coefficient 1, so its stored d (0 up to the noise of the last recompute) takes the step too.
            // `0xe45c0` (`0xe47c7..0xe47ff`): the leaving variable is the row entry with coefficient +1 (leaving below) or -1 (above, the row is negated);
            // a status -2 (at the upper bound) variable then has its d negated, so for it d = -(old + theta * coef) (GROW7 release call 311; PILOT.JA calls 30572, 30580 with a noise `old` of 7e-18).
            self.d[pv.leave] = if self.st[pv.leave] == FIXED && pv.above {
                theta
            } else if pv.above {
                -(old_leave + theta * lcoef)
            } else {
                old_leave - theta * lcoef
            };
        } else if refactor {
            self.refactor()?;
        }
        if refactor {
            self.pricing_switch();
            if self.request == 0 {
                self.request = 10;
            }
            self.recompute_d();
            let level_before = self.level;
            self.set_level_work(level_before != 0);
            self.set_level(level_before);
            let req = self.request;
            self.request = 0;
            self.recompute_xb(req);
            self.prog_on = self.level == 0;
            self.recompute_flags(level_before != 0, pv.is_some());
            self.flags_work();
            // The recompute pass `0xbf5b0` (called from `0xf6d0f`) runs the no-leaving-row handler itself when the pricing is in the list-walk mode
            // (`P+0x40 == 0`) and the list is empty (`P+0x3c == 0`, `0xbffbe`): no pricing call follows (RECIPE call 46).
            if !self.value_mode && self.list.is_empty() && self.status.is_none() {
                self.no_row();
            }
        } else if let Some(pv) = pv {
            self.work
                .add(0xc0410, 3.0 * pv.col.len() as f64, crate::work::CTX);
            // With `ctx+0x2c8` set (huge infeasibilities: level != 0) the value-mode walk is the 8-per-entry loop `0xc09bd`.
            if self.value_mode && self.level != 0 {
                self.work
                    .add(0xc09bd, 8.0 * pv.col.len() as f64, crate::work::CTX);
            } else {
                let site = if self.value_mode { 0xc0c24 } else { 0xc06c7 };
                self.work
                    .add(site, 7.0 * pv.col.len() as f64, crate::work::CTX);
            }
            self.update_primal_and_weights(pv);
        }
        Ok(())
    }

    /// x_B and DSE weight update (`0xc0050`, fast path), bit for bit as in the spec.
    fn update_weights(&mut self, pv: &Pivot) {
        let (r, alpha, w) = (pv.r, pv.alpha, pv.w);
        if crate::dump::on("dtau") {
            crate::dump::emit(
                "dtau",
                serde_json::json!({"r": r, "alpha": alpha, "w": w, "col": crate::dump::sparse(&pv.col), "tau": pv.tau}),
            );
        }
        if w == 0.0 {
            self.w[r] = W_EPS;
        } else {
            let a2 = if alpha >= 0.0 {
                maxsd(alpha, W_EPS)
            } else {
                minsd(alpha, -W_EPS)
            };
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
        let b = if self.st[v] == LB || self.st[v] == FIXED {
            self.p.lo[v]
        } else {
            self.p.up[v]
        };
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
        self.avg_touch = (2.0 * pv.col.len() as f64 + (it - 2.0) * self.avg_touch) / it;
        if self.value_mode {
            if ninf > self.avg_touch.max((self.m as f64).powf(0.45).floor()) {
                self.thr *= 1.1;
            }
            // The pivot row leaves the list at the next pricing call whatever its value (`0xc0732`); its priority is refreshed for the entering variable.
            // An entering variable that lands beyond its bound stays on the list with its infeasibility (BOEING2 call 57); a feasible one gets the distance
            // from its lower bound and flag 3.
            let q = self.head[r];
            let fi = self.infeas(r);
            if fi > self.ptol() {
                self.prio[r] = self.pri(r, fi);
            } else {
                let d = self.xb[r] - self.p.lo[q];
                self.prio[r] = self.pri(r, d);
                self.flags[r] = 3;
            }
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
/// The initial weight of a dense BTRAN result (`0xc392b..0xc3b08`): the squares are summed over all `m` slots of the dense vector (zeros included), so
/// each nonzero falls in the lane of its row index; the alignment peel was never needed in the captures (16-byte aligned buffers).
fn rho_norm_dense(rho: &[(usize, f64)], m: usize) -> f64 {
    let mut sq = vec![0.0f64; m];
    for &(i, v) in rho {
        sq[i] = v * v;
    }
    if m <= 6 {
        return sq.iter().fold(0.0, |a, x| a + x);
    }
    let end = m - (m & 7);
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
pub fn solve(
    model: &LpModel,
    start: Option<&Basis>,
    opts: &SolveOptions,
) -> Result<SolveResult, SolveError> {
    let p = build(model);
    let (n, m) = (p.n, p.m);
    let mut dropped = None;
    let mut supplied = None;
    let basis = match start {
        Some(b) => {
            let mut b = b.clone();
            b.validate(n, m).map_err(SolveError::InvalidBasis)?;
            b.basic.sort_unstable();
            supplied = Some(b.basic.clone());
            // A singular basis loses the columns that get no pivot in the reference's LU order and takes the slacks of the uncovered rows instead.
            dropped = crate::simplex::repair_basis(&p, &mut b).map(|(j, _)| j);
            b
        }
        None => Basis {
            basic: dual_crash_basic(&p),
            at_upper: vec![false; n + m],
        },
    };
    basis.validate(n, m).map_err(SolveError::InvalidBasis)?;
    crate::lu_order::reset_pair_gate();
    crate::lu_order::set_pool_dims(p.n, p.cols[..p.n].iter().map(Vec::len).sum());
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
        lu_thr: crate::lu_order::THRESHOLD,
        pmode: if start.is_some() { 3 } else { 1 },
        n_infeas: 0,
        dflag: Vec::new(),
        switch_enabled: opts.perturb > 0.0,
        level: 0,
        request: 0,
        iters: 0,
        prio: vec![0.0; m],
        ptype: 0,
        p68: Vec::new(),
        p70: Vec::new(),
        p60: 0.0,
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
        warm: start.is_some(),
        avg_touch: 0.0,
        last_cap: 0,
        updates: 0,
        upd_limit: None,
        pending_alpha: None,
        status: None,
        nm: model.num_cols(),
        ray: None,
        zero_artificial: false,
        cost: p.cost.clone(),
        cshift: vec![0.0; p.cost.len()],
        rng: 101,
        pert: 0,
        obj_off: 0.0,
        flip_select: start.is_some(),
        work: crate::work::Work::dual(),
        pivot_cols: None,
        nnz_a: p.cols[..n].iter().map(Vec::len).sum(),
        pstep: opts.perturb,
        stall: 0,
        stall_snap: 0,
        stall_it: 0.0,
        prog_on: false,
        best_obj: -1e100,
        pert_left: 1_000_000_000,
        pert_param: opts.perturb,
        obj_scale: opts.obj_scale,
        pcost: p.cost.clone(),
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
    if let Some(j) = dropped {
        // The variable the singular start basis lost stays out of the dual (status -4 like a fixed one, at its lower bound, d = its cost); the primal
        // pass after the dual optimum brings it back (`0x155139`).
        s.st[j] = FIXED;
    }
    s.init_pricing_type();
    s.work.decided(s.work.w);
    s.factor_work(&stc, &rows, 0.0);
    let free_phase = (0..n).any(|j| s.st[j] == FREE);
    if free_phase {
        s.free_pivots()?;
        s.refactor()?;
    }
    s.build_rrows();
    let dt = s.true_d();
    s.work
        .add_from(0x12de50, 0xe53a7, s.mirror.work, crate::work::FAC);
    s.request = 1;
    s.set_level_work(false);
    s.request = 0;
    for j in 0..n + m {
        if s.st[j] >= 0 {
            continue;
        }
        if s.st[j] == FREE && dt[j].abs() > DUAL_TOL {
            // A free variable the setup pivots could not bring into the basis (no basic variable with a finite lower bound in its column) and that is dual
            // infeasible: it takes the artificial bound its d asks for, like any other variable without the bound it would need.
            s.st[j] = if dt[j] > 0.0 { LB } else { UB };
            s.d[j] = dt[j].abs();
            continue;
        }
        if s.st[j] == FIXED || s.st[j] == FREE {
            s.d[j] = if s.st[j] == FIXED && j < n {
                s.cost[j]
            } else {
                dt[j]
            };
            continue;
        }
        if dt[j] < -DUAL_TOL {
            s.st[j] = UB;
        } else if dt[j] > DUAL_TOL {
            s.st[j] = LB;
        } else if p.lo[j] <= -BIG {
            // A d within the tolerance does not tie the variable to the artificial bound: the finite upper bound (the mirror of the default LB for a
            // variable with a finite lower bound; no corpus model has a variable without a lower bound and with an upper one).
            s.st[j] = UB;
        }
        s.d[j] = if s.st[j] == UB { -dt[j] } else { dt[j] };
    }
    if !free_phase {
        s.reset_pricing(&stc, &rows);
    }
    s.slack_nb = (0..m).map(|i| s.st[n + i] < 0).collect();
    s.level = (0..n + m).filter(|&j| s.at_artificial(j)).count();
    if let Some(sup) = &supplied {
        // The flags of the pricing mode 3 follow the supplied basis, before the repair of a singular one (a slack that replaced a dropped column is flagged
        // nonbasic, the dropped column basic).
        s.dflag = vec![true; n + m];
        for &j in sup {
            s.dflag[j] = false;
        }
    }
    // Initial weights: the squared norms of the rows of B^-1 (a supplied basis starts with unit weights).
    for r in 0..if s.warm { 0 } else { m } {
        let out = s.mirror.btran(&[(r, 1.0)], false, false);
        s.work
            .add_from(0x12de50, 0xc36f1, s.mirror.work, crate::work::FAC);
        s.w[r] = if s.mirror.bt_dense {
            rho_norm_dense(&out, m)
        } else {
            rho_norm(&out)
        };
    }
    s.stall_check(0);
    s.recompute_xb(0);
    s.prog_on = s.level == 0;
    s.work.add(0xc391f, m as f64, crate::work::CTX);
    s.recompute_flags(s.warm && s.level != 0, false);
    s.flags_work();
    let schedule = opts.refactor_schedule.clone();
    let mut trace = Vec::new();
    let mut k = 0usize;
    let mut farkas = None;
    let status;
    loop {
        k += 1;
        if s.iters >= opts.max_iterations {
            status = Status::IterationLimit;
            break;
        }
        let mut dump = None;
        let step = s.iteration(k, &mut dump);
        if let Some(mut rec) = dump {
            if crate::dump::on("dsxl") {
                let lean = ["k", "r", "var", "q", "flips"]
                    .iter()
                    .filter_map(|f| rec.get(*f).map(|v| (f.to_string(), v.clone())))
                    .collect::<serde_json::Map<_, _>>();
                rec = serde_json::Value::Object(lean);
            }
            crate::dump::emit(
                if crate::dump::on("dsxl") {
                    "dsxl"
                } else {
                    "dsx"
                },
                rec,
            );
        }
        match step {
            Step::Pivot(pv) => {
                if opts.record_trace {
                    trace.push(TraceEvent {
                        iteration: s.iters,
                        entering: pv.q,
                        leaving: Some(pv.leave),
                        position: Some(pv.r),
                        step: 0.0,
                        flips: pv.flips.clone(),
                    });
                }
                let due = match &schedule {
                    Some(sch) => sch.contains(&s.iters),
                    None => {
                        let w_now = s.work.w;
                        s.mirror.refactor_due() || s.work.clause_dual(s.updates, s.nnz_a, w_now)
                    }
                };
                s.post(Some(&pv), due)?;
            }
            Step::NoRow => {
                if s.request != 3 {
                    s.no_row();
                }
                if s.status.is_none() {
                    s.post(None, false)?;
                }
            }
            Step::Unbounded(y) => {
                status = Status::Infeasible;
                farkas = Some(y);
                break;
            }
        }
        if let Some(stt) = s.status {
            status = stt;
            if stt == Status::Optimal && s.pert != 0 {
                // The end of the dual with the cost shifts in force (`ctx+0x338` != 0, `0xf7033`): the shifts are removed (status -1004) and the primal
                // clean-up starts with a factorization (`0x142350` from `0xe3ee1`), the reference's last one.
                s.cost.clone_from(&p.cost);
                s.obj_off = 0.0;
                s.pert = 0;
                s.refactor()?;
                s.recompute_xb(0);
            }
            break;
        }
    }
    if status == Status::Optimal && dropped.is_some() {
        // The primal pass of the reference after the dual optimum: the dropped variable is priced again (a bound flip or a pivot, each one an iteration).
        let mut popts = opts.clone();
        popts.method = crate::simplex::Method::Primal;
        let b = Basis {
            basic: s.head.clone(),
            at_upper: s.st.iter().map(|&t| t == UB).collect(),
        };
        let mut r = crate::simplex::solve(model, Some(&b), &popts)?;
        r.iterations += s.iters;
        r.refactorizations += s.refactorizations;
        return Ok(r);
    }
    if status == Status::Unbounded {
        // The dual is infeasible (a nonbasic variable stayed at an artificial bound, or a free one has a nonzero d). That proves the primal unbounded only
        // if it is feasible, and the point the dual loop holds is a finite point plus 1e30 times the ray, not a feasible point. A zero objective is dual
        // feasible, so the dual simplex on it decides: a feasible point (unbounded) or a Farkas certificate (primal and dual infeasible).
        let zero = LpModel {
            obj: vec![0.0; model.num_cols()],
            offset: 0.0,
            ..model.clone()
        };
        // No cost perturbation here: the stall trigger would perturb the zero costs and the loop could end unbounded again, without end.
        let mut zopts = opts.clone();
        zopts.perturb = 0.0;
        let mut r = solve(&zero, start, &zopts)?;
        r.iterations += s.iters;
        r.refactorizations += s.refactorizations;
        if r.status == Status::Optimal {
            r.status = Status::Unbounded;
            r.objective = model.objective(&r.x);
            r.ray = s.ray.take().filter(|ray| {
                crate::verify::ray_improvement(model, ray, 1e-9).is_some_and(|g| g > 0.0)
            });
            if r.ray.is_none() {
                // The ray of the artificial-bound point did not verify; the primal simplex may still find one.
                let mut popts = opts.clone();
                popts.method = crate::simplex::Method::Primal;
                r.ray = crate::simplex::solve(model, None, &popts)?.ray;
            }
        }
        return Ok(r);
    }
    if status == Status::IterationLimit && s.level != 0 {
        s.zero_artificial = true;
        s.recompute_xb(0);
    }
    // Read-out.
    let sgn = if model.sense == Sense::Maximize {
        -1.0
    } else {
        1.0
    };
    let mut x = vec![0.0; n + m];
    for j in 0..n + m {
        x[j] = if s.st[j] >= 0 {
            s.xb[s.st[j] as usize]
        } else {
            s.nonbasic_value(j)
        };
    }
    s.cost.clone_from(&p.cost);
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
        basis: Basis {
            basic: s.head.clone(),
            at_upper,
        },
        iterations: s.iters,
        refactorizations: s.refactorizations,
        trace,
        farkas,
        ray: s.ray.take(),
    })
}
