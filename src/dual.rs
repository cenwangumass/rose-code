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
use crate::real::Real;
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
const DUAL_TOL0: f64 = 1e-6;

thread_local! {
    /// The dual feasibility tolerance (`ctx+0x388`) and the counter of failed dual residual checks (`ctx+0x174`); both outlive the restarts of one solve.
    static DUAL_TOL_CELL: std::cell::Cell<f64> = const { std::cell::Cell::new(DUAL_TOL0) };
    /// The primal feasibility tolerance in force (`ctx+0x380`); the primal's drift test raises it (up to 1e-3), the dual passes after it read the raised value.
    static PRIMAL_TOL_CELL: std::cell::Cell<f64> = const { std::cell::Cell::new(PRIMAL_TOL) };
    static RNG_CELL: std::cell::Cell<u64> = const { std::cell::Cell::new(101) };
    static PERT_RERUNS: std::cell::Cell<u32> = const { std::cell::Cell::new(0) };
    static RESID_FAILS: std::cell::Cell<u32> = const { std::cell::Cell::new(0) };
    /// Sum of dual infeasibilities from the last full recomputation (`ctx+0x2d8`).
    static DUAL_INFEAS_SUM: std::cell::Cell<f64> = const { std::cell::Cell::new(0.0) };
}

thread_local! {
    /// The number of passes started in this solve (selects the entry of a per-pass refactorization schedule).
    static LOOP_TOL: std::cell::Cell<f64> = const { std::cell::Cell::new(1e-6) };
    /// Dual passes ended with a primal handoff in this solve (`ctx+0x35c` - 1 outside the quad passes).
    static DUAL_HANDOFFS: std::cell::Cell<u32> = const { std::cell::Cell::new(0) };
    /// `ctx+0x35c`: 1 at the start, plus one for every run that ends with status -1004 (writers `0xf8691` double dual, `0xe404a` double primal, `0x24b98b` quad dual,
    /// `0x2371bd` x87 primal; `replay/gdb_c35c_trace.py`). The quad primals' method (`ctx+0x3ec`) and the switch of the first ones read it.
    pub(crate) static C35C: std::cell::Cell<u32> = const { std::cell::Cell::new(1) };
    static PASS_NO: std::cell::Cell<usize> = const { std::cell::Cell::new(0) };
    /// The stability counter `ctx+0x178`: one per pivot-check mismatch (`0xc5a35`, quad `0x21da43`) and per singular quad restart (`0x248f99`); never reset in a solve.
    pub(crate) static C178: std::cell::Cell<u32> = const { std::cell::Cell::new(0) };
    /// Large basic-value flag (`ctx+0x2c8`), retained across dual passes.
    static LARGE_XB: std::cell::Cell<bool> = const { std::cell::Cell::new(false) };
    /// The basis snapshot of the reference's log timer (`ctx+0x190`) and the iteration count of the last one (`ctx+0x40`).
    static SNAP: std::cell::RefCell<Option<(Vec<usize>, Vec<bool>)>> = const { std::cell::RefCell::new(None) };
    static INITIAL_SNAP: std::cell::RefCell<Option<(Vec<usize>, Vec<bool>)>> = const { std::cell::RefCell::new(None) };
    /// Snapshot quality (`ctx+0x188`).
    static SNAP_QUALITY: std::cell::Cell<u32> = const { std::cell::Cell::new(0) };
    /// The relative tolerance of the primal's pivot-pair check (`ctx+0x3a0`): 3e-6 from the parameters, 1e-7 from the first restore of a dropped variable on (`0x1551d4`).
    static PAIR_TOL: std::cell::Cell<f64> = const { std::cell::Cell::new(3e-6) };
    /// The primal's start basis when its drift test ended it (status -1002): the restart takes it (`0xcedb0`, seed 2717).
    static DRIFT_START: std::cell::RefCell<Option<Basis>> = const { std::cell::RefCell::new(None) };
    static SNAP_C40: std::cell::Cell<f64> = const { std::cell::Cell::new(-100.0) };
}

/// The log timer of the primal (`ctx+0x40`): true once `total` is 100 past the last check, which it then sets (`(total / 100) * 100`, `0x27e1a1..0x27e1a9`).
pub(crate) fn timer_check(total: f64) -> bool {
    let c40 = SNAP_C40.with(|c| c.get());
    if total - c40 < 100.0 {
        return false;
    }
    SNAP_C40.with(|c| c.set((total / 100.0) * 100.0));
    true
}

pub(crate) fn resid_fails() -> u32 {
    RESID_FAILS.with(|c| c.get())
}

pub(crate) fn set_resid_fails(v: u32) {
    RESID_FAILS.with(|c| c.set(v));
}

pub(crate) fn dual_infeas_sum() -> f64 {
    DUAL_INFEAS_SUM.with(|c| c.get())
}

pub(crate) fn set_dual_infeas_sum(v: f64) {
    DUAL_INFEAS_SUM.with(|c| c.set(v));
}

pub(crate) fn primal_tol() -> f64 {
    PRIMAL_TOL_CELL.with(|c| c.get())
}

pub(crate) fn set_primal_tol(v: f64) {
    PRIMAL_TOL_CELL.with(|c| c.set(v));
}

pub(crate) fn dual_tol() -> f64 {
    DUAL_TOL_CELL.with(|c| c.get())
}

pub(crate) fn set_dual_tol(v: f64) {
    DUAL_TOL_CELL.with(|c| c.set(v));
}

const PIVOT_FLOOR: f64 = 3e-8;
const DROP_TOL: f64 = 1e-13;
const ALPHA_TOL: f64 = 1e-10;
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
pub(crate) fn accuracy_threshold(m0: f64, e2: f64) -> f64 {
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

fn minr<R: Real>(a: R, b: R) -> R {
    if a < b {
        a
    } else {
        b
    }
}

fn maxr<R: Real>(a: R, b: R) -> R {
    if a > b {
        a
    } else {
        b
    }
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
    /// The row is rejected (`0xc7ca4..0xc7e4f`): its largest entry is below the pivot tolerance; the measure stored for the worst rejected row.
    reject: Option<f64>,
}

struct Dual<'a, R: Real = f64> {
    p: &'a Internal,
    n: usize,
    m: usize,
    head: Vec<usize>,
    /// Status per variable: the basis position when basic, else `LB`, `UB`, `FREE`, `FIXED`.
    st: Vec<i32>,
    /// Basic values by position.
    xb: Vec<R>,
    /// Reduced costs in the reference's frame.
    d: Vec<R>,
    /// Dual steepest-edge weights by position.
    w: Vec<R>,
    mirror: RefStruct<R>,
    factor: Factor,
    /// Row-wise copy of the nonbasic, non-fixed structurals in the reference's storage order.
    rrows: Vec<Vec<usize>>,
    /// The ratio test asked for a refactorization (request 5) instead of ending the dual: the iteration returns without a pivot.
    absorb_request: bool,
    refactorizations: usize,
    /// The relative pivot threshold of the LU (`ctx+0x390`, 2^-7 at the start): raised after a recompute of x_B by `accuracy_threshold`, never lowered.
    lu_thr: f64,
    /// Whether either side of the last basic-value recompute has large values.
    big_pricing: bool,
    /// The pricing mode `ctx+0x3ec`: 1 for dual steepest edge, 3 after the switch of `0xf770d` (weights restarted at 1).
    pmode: u8,
    /// The number of infeasible rows at the last value-mode recompute (`P+0x3c`).
    n_infeas: usize,
    /// Rows rejected by the ratio test since the last factorization (`ctx+0x420` count), the worst measure (`ctx+0x438`) and its row (`ctx+0x430`).
    bad_rows: Vec<usize>,
    /// Rejected rows at the entry of the current driver call (`0xf66fb`: `ctx+0x420` copied into the local that `0xc50e4` reads): the pricing that finds no row
    /// releases that many of the oldest ones (`0xcd250`) before it gives up.
    old_rej: usize,
    /// The pass ended because every infeasible row was rejected (hand over to the primal like the repair case).
    rejected_end: bool,
    /// The pass ended at the request-4 refactorization with variables at an artificial bound (status -1004): the primal decides the outcome.
    art_exit: bool,
    worst: f64,
    worst_row: i64,
    /// The retry of the worst rejected row (`ctx+0x420` negated, `0x21cf50` / `0xc4fe0`): the ratio test rejects nothing.
    retry: bool,
    /// The end of a retried row whose largest candidate is below a tenth of the pivot floor (`0xcc3bc`, quad `0x223f49`): 0 none, 1 the dual ends infeasible (status 3),
    /// 2 a refactorization is requested (5).
    retry_end: u8,
    /// `ctx+0x4d4`: the retried rows of that kind so far (at most 10 before the dual ends infeasible).
    retry_tiny: u32,
    /// The infeasible-row count of the last flag recompute, whatever the mode (`[P+0x3c]`); a restart after a singular factorization compares the new count with it.
    inf_count: usize,
    /// Restarts so far (the reference's `ctx+0x184`).
    retries: u32,
    /// The variable the in-place repair after the retry limit left out of the dual (`[ctx+0x454]` list, status -4); the controller restores it and runs the primal.
    repaired: Option<usize>,
    snapshot_excluded: bool,
    /// Set at a restart: the next flag recompute compares the stale `inf_count`.
    stale_count: bool,
    /// The last `post` call followed a pivot (the leaving row `ctx+0x2f8` is set).
    last_pivoted: bool,
    /// `Carry::it_base`.
    it_base: usize,
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
    prio: Vec<R>,
    /// Pricing type `P+0x44`: 0 plain dual steepest edge priorities, 1 scaled by `P70[var] / P68[row]`, 2 by `1 / P68[row]`.
    ptype: u8,
    /// Per-row weights `P+0x68` (by position) and per-variable constants `P+0x70` of the scaled priorities.
    p68: Vec<f64>,
    p70: Vec<f64>,
    /// Touch increment `P+0x60`: the number of entries of the pivot row's BTRAN result.
    p60: f64,
    /// `P+0x58`: the entries the row products have visited since the pass (or the last weight reset) started (the row lengths of the nonzero rho entries).
    p58: f64,
    flags: Vec<u8>,
    list: Vec<usize>,
    thr: R,
    /// Relative residual of the last full basic-solution recomputation.
    xb_residual: R,
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
    /// The duals `y = B^-T c_B` of the last recompute.
    last_y: Vec<f64>,
    /// The factorization count the pass started from (`ctx+0x2c0` counts from there).
    fac0: usize,
    /// `ctx+0x2c0` at the start of the pass (`Carry::nref0`), and the factorizations that `refactorizations` counts but `ctx+0x2c0` does not (the failed one of the in-place repair:
    /// fx_midsing_1031 and 1390, `dual_acc` captures).
    nref0: usize,
    nref_skip: usize,
    /// The iterations of the earlier passes (`ctx+0x20`).
    it0: usize,
    /// Updates of the factor since the last refactorization.
    updates: usize,
    /// The update limit of the factor object once the value mode has started (`F+0x188`).
    upd_limit: Option<usize>,
    /// The pivot column of the last iteration by position (for the factor's update).
    pending_alpha: Option<Vec<R>>,
    /// The pivot-check pair of the last iteration that reached the check (`ctx+0x300` = the column's pivot, `ctx+0x308` = the row's): the primal started by the
    /// controller reads it at its first pricing call.
    pair: (f64, f64),
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
    /// The part of a shifted cost the double does not hold (the quad pass keeps the cost shifted in long double).
    cost_lo: Vec<f64>,
    /// The quad instance's shift array (`0xc0(rsp)` of `0x237750`): the accumulated shifts and perturbations per variable in long double; the cost used is
    /// the double base cost plus this, rounded once to long double.
    shift_r: Vec<R>,
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
    /// The progress accumulator `ctx+0x100` (objective gain since the last x_B recompute, scaled by `ctx+0xf8`).
    prog_acc: f64,
    best_obj: f64,
    pert_left: i64,
    pert_param: f64,
    obj_scale: f64,
    /// The model was scaled (`ctx+0xe0` non-null): the growth check of the factorization uses 1e12.
    scaled: bool,
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
pub(crate) fn dual_crash_basic(p: &Internal) -> (Vec<usize>, Vec<bool>) {
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
    let upper = crash_pass3(p, &mut basic, &mut cost, &taken, &is_eq, &row_cols);
    (basic, upper)
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
) -> Vec<bool> {
    let (n, m) = (p.n, p.m);
    // -4 fixed, -1 at lower, -2 at upper, >= 0 basic in that row.
    let mut st: Vec<i64> = (0..n + m)
        .map(|j| {
            if p.up[j] - p.lo[j] <= 0.0 {
                -4
            } else if p.cost[j] < 0.0 {
                -2
            } else {
                -1
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
    if crate::dump::on("dcr") {
        crate::dump::emit("dcr", serde_json::json!({"crash_st": st}));
    }
    st.iter().map(|&v| v == -2).collect()
}

/// The reference's pivot order of the basis columns (`lu_order`), as in the primal path.
/// The `lustart` record of a pass's starting factorization.
fn dump_start<R: Real>(
    cols: &[usize],
    order: &[usize],
    st: &crate::lu_order::Structure<R>,
    thr: f64,
) {
    if crate::dump::on("lustart") {
        crate::dump::emit(
            "lustart",
            serde_json::json!({"cols": cols, "out": order, "prow": st.pivots.iter().map(|p| p.0).collect::<Vec<_>>(), "thr": thr}),
        );
    }
}

fn reference_order<R: Real>(
    p: &Internal,
    basic: &[usize],
    fac_no: usize,
) -> (
    Vec<usize>,
    crate::lu_order::Structure<R>,
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
    let st = crate::lu_order::factor_structure_r::<R>(p.m, &cols);
    let order = st.pivots.iter().map(|&(_, pos)| basic[pos]).collect();
    (order, st, cols)
}

thread_local! {
    static QUAD_RESULT: std::cell::Cell<bool> = const { std::cell::Cell::new(false) };
    /// The scaled long-double solution of a quad final pass as (high, low) double pairs, for the reference's objective read-out.
    static QUAD_XR: std::cell::RefCell<Option<Vec<(f64, f64)>>> = const { std::cell::RefCell::new(None) };
}

/// The objective read-out accumulates in the working precision of the final pass (`long double` in the quad pass: seed 2150, the double sum is one
/// bit off), the offset added last.
fn objective_in<R: Real>(model: &crate::model::LpModel, x: &[f64]) -> f64 {
    QUAD_RESULT.with(|c| c.set(R::QUAD));
    let mut acc = R::ZERO;
    for (c, v) in model.obj.iter().zip(x) {
        acc = acc + R::c(*c) * R::c(*v);
    }
    (acc + R::c(model.offset)).to_f64()
}

/// Preserve the x87 primal's freshly recomputed basic values for the unscaled objective.
pub(crate) fn set_x87_result(x: Vec<(f64, f64)>) {
    QUAD_XR.with(|c| *c.borrow_mut() = Some(x));
    QUAD_RESULT.with(|c| c.set(true));
}

/// The objective of the unscaled solution of the dual's last result: `long double` accumulation when its final pass was the quad pass. Clears the mark.
pub fn unscaled_objective(model: &crate::model::LpModel, x: &[f64], col: &[f64]) -> f64 {
    if QUAD_RESULT.with(|c| c.replace(false)) {
        use crate::f80::F80;
        if let Some(xr) = QUAD_XR.with(|c| c.borrow_mut().take()) {
            if xr.len() >= model.obj.len() {
                // The scale factors are powers of two: the scaled long-double value times the column scale is exact.
                let mut acc = F80::ZERO;
                for (j, c) in model.obj.iter().enumerate() {
                    let v =
                        (F80::from_f64(xr[j].0) + F80::from_f64(xr[j].1)) * F80::from_f64(col[j]);
                    acc = acc + F80::from_f64(*c) * v;
                }
                return (acc + F80::from_f64(model.offset)).to_f64();
            }
        }
        objective_in::<crate::f80::F80>(model, x)
    } else {
        model.objective(x)
    }
}

impl<'a, R: Real> Dual<'a, R> {
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

    /// The reference's factorization count `ctx+0x2c0`.
    fn nref(&self) -> usize {
        self.refactorizations - self.fac0 - self.nref_skip + self.nref0
    }

    /// Primal feasibility tolerance in force (`0x380`, or `(0x370 * 1e-5) * 0x380` while the level is not zero).
    fn ptol(&self) -> f64 {
        if self.level != 0 {
            BIG * 1e-5 * primal_tol()
        } else {
            primal_tol()
        }
    }

    fn infeas(&self, r: usize) -> f64 {
        self.infeas_real(r).to_f64()
    }

    fn infeas_real(&self, r: usize) -> R {
        let v = self.head[r];
        let x = self.xb[r];
        let (lo, up) = (self.p.lo[v], self.p.up[v]);
        if lo > -BIG && x < R::c(lo) {
            R::c(lo) - x
        } else if up < BIG && x > R::c(up) {
            x - R::c(up)
        } else {
            R::ZERO
        }
    }

    // ------------------------------------------------------------------ factorization

    /// One call of the reference's LU (`0xd4190`) on the current head: the pivot order, the structure and, when every column found a pivot, the fresh mirror. Each call is
    /// dumped (`luin`, `luout`), the failed ones and the retries inside the driver included.
    fn lu_attempt(
        &mut self,
    ) -> (
        Vec<usize>,
        crate::lu_order::Structure<R>,
        Vec<(bool, Vec<usize>, Vec<f64>)>,
        Option<RefStruct<R>>,
    ) {
        crate::lu_order::set_threshold(self.lu_thr);
        let (order, st, rows) = reference_order::<R>(self.p, &self.head, self.refactorizations + 1);
        let mi = (order.len() == self.m).then(|| RefStruct::fresh(self.m, self.n, &st, &rows));
        if crate::dump::on("luin") {
            crate::dump::emit(
                "luin",
                serde_json::json!({"n": self.refactorizations, "cols": self.head, "out": order, "prow": st.pivots.iter().map(|p| p.0).collect::<Vec<_>>(), "phase": st.phase, "thr": self.lu_thr}),
            );
        }
        if let (Some(mi), true) = (mi.as_ref(), crate::dump::on("luout")) {
            let (diag, l) = mi.numerics();
            let diag: Vec<f64> = diag.iter().map(|v| v.to_f64()).collect();
            let l: Vec<Vec<(usize, f64)>> = l
                .iter()
                .map(|c| c.iter().map(|e| (e.0, e.1.to_f64())).collect())
                .collect();
            let u: Vec<Vec<(usize, f64)>> = mi
                .u_numerics()
                .iter()
                .map(|c| c.iter().map(|e| (e.0, e.1.to_f64())).collect())
                .collect();
            crate::dump::emit(
                "luout",
                serde_json::json!({"n": self.refactorizations, "diag": diag, "l": l, "u": u}),
            );
        }
        (order, st, rows, mi)
    }

    /// The payload of a pass that ends in a restart (`SINGULAR`): the state the next pass inherits.
    fn singular_payload(&self) {
        // The double driver's singular handler increments the stability counter while
        // its retry count is below five (`0xf58f7`, default environment).
        if R::QUAD || self.retries < 5 {
            C178.with(|c| c.set(c.get() + 1));
        }
        SINGULAR_BASIS.with(|c| {
            *c.borrow_mut() = Some(Basis {
                basic: self.head.clone(),
                at_upper: self.st.iter().map(|&t| t == UB).collect(),
            })
        });
        SINGULAR.with(|c| {
            c.set(Some((
                self.iters,
                self.refactorizations + 1,
                self.lu_thr,
                self.avg_touch,
                // The pivot-column length carries over from any pass (`[P+0x?]` read at `0xc1e9e`: 1494 compares with 6 after a list-walk pass).
                self.last_cap,
                self.inf_count,
                self.last_pivoted,
                self.retries,
                self.it_base + self.iters,
            )))
        });
    }

    /// Factorize the basis (`0xd4190`): the positions follow the pivot order of the reference's LU; the weights move with their variables.
    fn refactor(&mut self, cause: &str) -> Result<(), SolveError> {
        if crate::dump::on("facctx") {
            crate::dump::emit(
                "facctx",
                serde_json::json!({"it": self.iters, "it_base": self.it_base, "total_it": self.it_base + self.iters,
                    "n": self.refactorizations, "req": self.request, "cause": cause, "quad": R::QUAD,
                    "level": self.level, "updates": self.updates, "pert": self.pert}),
            );
        }
        // LU driver (`0xf32a4..0xf32ec`): a factorization requested as 5, 6 or 7 (no entering variable, pivot mismatch) first raises `ctx+0x390` to 0.5.
        if (5..=7).contains(&self.request) && self.lu_thr < 0.5 {
            self.lu_thr = 0.5;
        }
        let (mut order, mut st, mut rows, mut mi) = self.lu_attempt();
        if order.len() < self.m && self.retries >= 5 {
            // After five restarts the driver refuses another one (`0xf3aba`) and repairs in place: the columns without a pivot leave the basis and the
            // slacks of the rows left without a pivot replace them (the supplied-basis repair); the last variable dropped stays out of pricing (status -4).
            // The failed factorization counts as one (`ctx+0x2c0`).
            let n = self.n;
            let (mut row_done, mut pos_done) = (vec![false; self.m], vec![false; self.m]);
            for &(r, pos) in &st.pivots {
                row_done[r] = true;
                pos_done[pos] = true;
            }
            let free_rows = (0..self.m).filter(|&r| !row_done[r]);
            let gone = (0..self.m).filter(|&pos| !pos_done[pos]);
            let mut last = None;
            for (pos, r) in gone.zip(free_rows) {
                let out = self.head[pos];
                self.st[out] = LB;
                self.head[pos] = n + r;
                last = Some(out);
            }
            if let Some(j) = last {
                self.st[j] = FIXED;
                self.d[j] = R::c(if j < n { self.cost[j] } else { 0.0 });
                self.repaired = Some(j);
            }
            // The post pass after the repair runs as at a normal start (shifts, no flip to the artificial bound): the request is not 4.
            self.request = 10;
            self.refactorizations += 1;
            self.nref_skip += 1;
            // The failure counts as one more retry (`ctx+0x184` 5 -> 6 at `0xf5921`, fx_midsing_1390): later pivot mismatches on a fresh factor end the run.
            self.retries += 1;
            // The in-place repair also counts for the stability counter (`ctx+0x178`; quad seeds 2081, 2448).
            C178.with(|c| c.set(c.get() + 1));
            (order, st, rows, mi) = self.lu_attempt();
        }
        if order.len() < self.m {
            // A basis column without a pivot: the reference restarts the dual from the slack basis (spec, "Singular factorization in the middle of a dual run").
            self.singular_payload();
            return Err(SolveError::SingularBasis);
        }
        let discard = if self.pivot_cols.is_some() {
            self.mirror.discard_work()
        } else {
            0.0
        };
        // Growth check after the factorization (`0xf351d..0xf3629`): while the pivot threshold is below 0.9 and max |diag| / min |diag| of the new factor exceeds 1e15
        // (1e12 when the model is scaled, `ctx+0xe0`), the threshold rises (0.5, then 0.9) and the same columns are factorized again (a retry inside the driver: the
        // factorization count `ctx+0x2c0` counts once; fx_midsing_1390 factorizations 3 and 4 both have a ratio of 3.7e16).
        loop {
            if self.lu_thr >= 0.9 {
                break;
            }
            let (mx, mn) = mi
                .as_ref()
                .unwrap()
                .numerics()
                .0
                .iter()
                .fold((0.0f64, f64::INFINITY), |(a, b), v| {
                    (a.max(v.abs().to_f64()), b.min(v.abs().to_f64()))
                });
            if !(mx / mn > if self.scaled { 1e12 } else { 1e15 }) {
                break;
            }
            self.factor_work(&st, &rows, discard);
            self.lu_thr = if self.lu_thr < 0.5 { 0.5 } else { 0.9 };
            (order, st, rows, mi) = self.lu_attempt();
            if order.len() < self.m {
                self.singular_payload();
                return Err(SolveError::SingularBasis);
            }
        }
        let mi = mi.unwrap();
        let mut wv = vec![R::ZERO; self.n + self.m];
        for (pos, &v) in self.head.iter().enumerate() {
            wv[v] = self.w[pos];
        }
        for (pos, &v) in order.iter().enumerate() {
            self.st[v] = pos as i32;
            self.w[pos] = wv[v];
        }
        self.head = order;
        self.reset_pricing(&st, &rows);
        self.mirror = mi;
        self.mirror.upd_limit = self.upd_limit;
        self.factor_work(&st, &rows, discard);
        let cols: Vec<&Column> = self.head.iter().map(|&j| &self.p.cols[j]).collect();
        self.factor =
            Factor::with_floor(self.m, &cols, 0.0).map_err(|_| SolveError::SingularBasis)?;
        self.refactorizations += 1;
        self.updates = 0;
        self.bad_rows.clear();
        self.old_rej = 0;
        self.worst = 0.0;
        self.worst_row = -2;
        self.slack_nb = (0..self.m).map(|i| self.st[self.n + i] < 0).collect();
        Ok(())
    }

    /// The free-variable pivots of the setup (`0xe88xx`, before the iteration loop; `replay/gdb_dual_free.py`): every free nonbasic structural, in
    /// ascending index, is pivoted into the basis. Its column is transformed; among the rows whose basic variable has a finite lower bound the entry of
    /// largest magnitude wins (first in the result's order on ties; below `PIVOT_FLOOR` the variable is deferred, unexercised). The leaving variable
    /// goes to its lower bound (fixed when the range is zero). The factor is updated, or refactorized when the work-based trigger is due.
    fn free_pivots(&mut self) -> Result<Option<(usize, f64)>, SolveError> {
        let n = self.n;
        let free: Vec<usize> = (0..n).filter(|&j| self.st[j] == FREE).collect();
        let mut deferred = Vec::new();
        for j in free {
            match self.free_candidate(j, PIVOT_FLOOR) {
                Some((r, alpha)) => self.pivot_free(j, r, alpha)?,
                None => deferred.push(j),
            }
        }
        if deferred.is_empty() {
            return Ok(None);
        }
        // Second pass (`0xe8fc6..0xe98e1`): request 10, a refactorization; then for each deferred variable in turn the reduced cost `d = c - y a` with
        // `y = B^-T c_B` (recomputed after every pivot of this pass). Within `dual_tol()` the variable stays free and nonbasic; beyond it the column is
        // transformed again and pivoted in on the largest entry among the rows whose basic variable has a finite lower bound (floor `0.01 PIVOT_FLOOR`);
        // without such an entry the dual is infeasible at once (`0xea0b7`: status 5, `ctx+0x2fc` the variable), before any iteration. `s70`, `s143`.
        self.request = 10;
        self.refactor("free_deferred")?;
        self.request = 0;
        let mut k = 0;
        while k < deferred.len() {
            let mut y: Vec<f64> = self.head.iter().map(|&v| self.cost[v]).collect();
            self.factor.btran(&mut y);
            while k < deferred.len() {
                let j = deferred[k];
                k += 1;
                let mut d = self.cost[j];
                for &(i, a) in &self.p.fcols[j] {
                    d -= y[i] * a;
                }
                if d.abs() <= dual_tol() {
                    continue;
                }
                match self.free_candidate(j, 0.01 * PIVOT_FLOOR) {
                    Some((r, alpha)) => self.pivot_free(j, r, alpha)?,
                    None => return Ok(Some((j, d))),
                }
                break;
            }
        }
        Ok(None)
    }

    /// The pivot row of a free column: the column is transformed and the entry of largest magnitude among the rows whose basic variable has a finite lower
    /// bound wins (first in the result's order on ties); `None` below `floor` (or without such a row).
    fn free_candidate(&mut self, j: usize, floor: f64) -> Option<(usize, Vec<R>)> {
        let m = self.m;
        let rows_q: Vec<usize> = self.p.fcols[j].iter().map(|&(i, _)| i).collect();
        let vals_q: Vec<R> = self.p.fcols[j].iter().map(|&(_, a)| R::c(a)).collect();
        let mut col_out = self.mirror.ftran(&rows_q, &vals_q);
        if self.mirror.last_dense() {
            col_out.sort_by_key(|e| e.0);
        }
        let mut best: Option<(usize, R)> = None;
        let mut mx = R::ZERO;
        for &(pos, v) in &col_out {
            if self.p.lo[self.head[pos]] > -BIG && v.abs() > mx {
                mx = v.abs();
                best = Some((pos, v));
            }
        }
        let (r, _) = best.filter(|_| mx >= R::c(floor))?;
        let mut alpha = vec![R::ZERO; m];
        for &(i, v) in &col_out {
            alpha[i] = v;
        }
        Some((r, alpha))
    }

    /// One pivot of the setup: the free variable `j` replaces the basic variable of position `r`, which goes to its lower bound (fixed when the range is
    /// zero); the factor is updated, or refactorized when the work-based trigger is due.
    fn pivot_free(&mut self, j: usize, r: usize, alpha: Vec<R>) -> Result<(), SolveError> {
        if crate::dump::on("dcr") {
            crate::dump::emit(
                "dcr",
                serde_json::json!({"free": j, "row": r, "val": alpha[r].to_f64(), "it": self.iters}),
            );
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
            self.refactor("free_pivot")?;
        } else {
            self.factor
                .update(r, alpha.iter().map(|v| v.to_f64()).collect());
            self.mirror.update(r);
            self.mirror.upd_work.clear();
            self.updates += 1;
        }
        Ok(())
    }

    /// The work of a factorization and the installation of its factor (`0xd4190` from `0xf34f2`, the update routine `0xf15c0`), after the fresh mirror is built:
    /// the driver's sites, the spike discard and the weight update of the pivot that preceded it, then the factor object's own passes.
    fn factor_work(
        &mut self,
        st: &crate::lu_order::Structure<R>,
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
        let mut r: Vec<R> = p.rhs.iter().map(|&v| R::c(v)).collect();
        for j in 0..p.n + p.m {
            if self.st[j] < 0 {
                let v = self.nonbasic_value(j);
                if v != 0.0 {
                    for &(i, a) in &p.cols[j] {
                        r[i] -= R::c(a) * R::c(v);
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
        let rows: Vec<usize> = (0..p.m).filter(|&i| r[i] != R::ZERO).collect();
        let vals: Vec<R> = rows.iter().map(|&i| r[i]).collect();
        let out = self.mirror.ftran_basic(&rows, &vals);
        let mut xb = vec![R::ZERO; p.m];
        for (i, v) in out {
            xb[i] = v;
        }
        // The reference forms B*x_B from zero first, then subtracts it from the
        // right-hand side (`0x115d60`, `0xc1480`); subtracting each term from
        // the right-hand side directly can change the residual and its gate.
        let mut product = vec![R::ZERO; p.m];
        for (pos, &j) in self.head.iter().enumerate() {
            for &(i, a) in &p.cols[j] {
                product[i] += R::c(a) * xb[pos];
            }
        }
        // `0xc10ab..0xc1129`, `0xc1480..0xc1517`: R = max(1e-9, max|rhs|), the largest |x_B| and the residual of B x_B against the rhs, both relative to R;
        // the residual is read only when the factorization count is 4 mod 5 or the request is 4, 6 or 8 (`0xc1137..0xc11ce`, `0xc165b..0xc1699`).
        let rmax = r.iter().fold(R::c(1e-9), |a, v| a.max(v.abs()));
        let xmax = xb.iter().fold(R::ZERO, |a, v| a.max(v.abs()));
        let e2 = if self.nref() % 5 == 4 || matches!(req, 4 | 6 | 8) {
            let residual = r
                .iter()
                .zip(&product)
                .fold(R::ZERO, |a, (b, v)| a.max((*b - *v).abs()));
            if self.level == 0 && residual > R::c(primal_tol()) {
                self.work.add(
                    0xc177e,
                    4.0 * p.m as f64 + 2.0 * cols as f64,
                    crate::work::CTX,
                );
            }
            residual / rmax
        } else {
            R::ZERO
        };
        self.xb_residual = e2;
        let (rmax, xmax, e2) = (rmax.to_f64(), xmax.to_f64(), e2.to_f64());
        if !R::QUAD {
            let large = rmax + xmax > 1e20;
            self.big_pricing = LARGE_XB.with(|c| {
                let old = c.get();
                c.set(large);
                old || large
            });
        }
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
        // The reference keeps the solution of its own factor whatever its residual (no refinement): with an ill-conditioned factor x_B is garbage of the size
        // 1e12 and more, and the next passes work on it (fx_midsing_1390 calls 7, 46).
        self.xb = xb;
    }

    /// d recomputed from the costs after a factorization, in the reference frame.
    fn recompute_d(&mut self) {
        let dt = self.true_d();
        self.dual_residual_check();
        self.work
            .add_from(0x12de50, 0xe53a7, self.mirror.work, crate::work::FAC);
        for j in 0..self.n + self.m {
            if self.st[j] == FIXED {
                if j < self.n {
                    // The cost in force: the double cost and, in the quad pass, its long-double remainder (seeds 8531, 7084).
                    self.d[j] = R::c(self.cost[j]) + R::c(self.cost_lo[j]);
                } else {
                    // The row product subtracts y from every slack irrespective of its status (`0x136b5d..0x136b85`, quad `0x273cd0`).
                    self.d[j] = dt[j];
                }
            } else if self.st[j] < 0 {
                self.d[j] = if self.st[j] == UB { -dt[j] } else { dt[j] };
            } else {
                self.d[j] = R::ZERO;
            }
        }
    }

    /// The log timer of the wrapper (`0xf9306..0xf6e13`): on a fresh factor with no request pending, the first time the total iteration count is 100 past the last
    /// snapshot, the basis is saved (`0xd40b0`, only while the restart counter is below 5). The first one (at the start) is the basis a restart from scratch builds.
    fn log_snapshot(&mut self) {
        if self.updates != 0 || self.request != 0 || self.retries >= 5 || self.status.is_some() {
            return;
        }
        let total = (self.it0 + self.iters) as f64;
        let c40 = SNAP_C40.with(|c| c.get());
        if total - c40 < 100.0 {
            return;
        }
        SNAP_C40.with(|c| c.set(total));
        let quality = if dual_infeas_sum() > 0.0 { 1 } else { 2 };
        // The dual snapshot (`0xd40b0`) skips a dropped variable and keeps a
        // previous snapshot of higher quality. The forced-row stack gate
        // remains unmodelled (no witness in the release profile).
        if self.snapshot_excluded
            || self.repaired.is_some()
            || SNAP_QUALITY.with(|c| quality < c.get())
        {
            return;
        }
        SNAP_QUALITY.with(|c| c.set(quality));
        let basic = self.head.clone();
        let at_upper = self.st.iter().map(|&t| t == UB).collect();
        if c40 == -100.0 {
            INITIAL_SNAP.with(|c| *c.borrow_mut() = Some((basic, at_upper)));
        } else {
            SNAP.with(|c| *c.borrow_mut() = Some((basic, at_upper)));
        }
    }

    /// Full dual recompute of `ctx+0x2d8` (`0xbf5b0`, quad `0x217600`).
    /// Its terms are positive reduced costs at artificial upper bounds.
    fn snapshot_recomputed(&self) {
        if R::QUAD && !self.value_mode {
            return;
        }
        let tol = R::c(dual_tol());
        let mut sum = R::ZERO;
        if self.level != 0 {
            for j in 0..self.n + self.m {
                if self.st[j] == UB && self.p.up[j] >= BIG && self.d[j] > tol {
                    sum = sum + self.d[j];
                }
            }
        }
        set_dual_infeas_sum(sum.to_f64());
    }

    /// The dual residual check of the post-refactor pass (`0xe5965..0xe5d4c`): `max |c_j - y a_j|` over the basic variables, read when the factorization count is 4 mod 5
    /// or the request is 4 or 6. A residual above the tolerance counts a failure; at the sixth the tolerance grows tenfold (up to 1e-3) and the counter restarts.
    fn dual_residual_check(&mut self) {
        let req = self.request;
        if !(self.nref() % 5 == 4 || req == 4 || req == 6) {
            return;
        }
        let mut meas = 0.0f64;
        for &j in &self.head {
            // The kernel of `0xe59a1..0xe5bd9`: the cost (with the shift), minus the column products, sequentially below 8 entries past the alignment peel, else the eight-lane
            // body of the read-out kernel (fx_midsing_1127: 165.177734375 against 142.407).
            let acc = crate::simplex::dot_lanes_from(
                &self.p.fcols[j],
                &self.last_y,
                self.p.offs[j] & 1,
                self.cost[j],
            );
            meas = meas.max(acc.abs());
        }
        let tol = dual_tol();
        if crate::dump::on("dualtol") {
            crate::dump::emit(
                "dualtol",
                serde_json::json!({"meas": meas, "tol": tol, "cnt": RESID_FAILS.with(|c| c.get()), "req": req, "n": self.nref()}),
            );
        }
        let mut cnt = RESID_FAILS.with(|c| c.get());
        if meas > tol {
            cnt += 1;
        }
        if cnt > 5 && tol < 1e-3 {
            set_dual_tol(tol * 10.0);
            cnt = 0;
        }
        RESID_FAILS.with(|c| c.set(cnt));
    }

    /// d from the costs (`d_j = c_j - y a_j`, in the true frame) with y = B^-T c_B from the reference BTRAN; the terms are subtracted row by row in
    /// ascending row order, as the reference's row-wise product `0x136b00` does (rows with |y| <= 1e-13 are skipped).
    fn true_d(&mut self) -> Vec<R> {
        let p = self.p;
        let inp: Vec<(usize, R)> = self
            .head
            .iter()
            .enumerate()
            .map(|(pos, &j)| (pos, R::c(self.cost[j]) + R::c(self.cost_lo[j])))
            .filter(|e| e.1 != R::ZERO)
            .collect();
        let out = self.mirror.btran(&inp, false, true);
        if crate::dump::on("dual") {
            crate::dump::emit(
                "dual",
                serde_json::json!({"refactor": self.refactorizations, "iter": self.iters, "out": crate::dump::sparse_r(&out)}),
            );
        }
        let mut y = vec![R::ZERO; p.m];
        for &(i, v) in out.iter() {
            y[i] = v;
        }
        self.last_y = y.iter().map(|v| v.to_f64()).collect();
        let mut d = vec![R::ZERO; p.n + p.m];
        for j in 0..p.n + p.m {
            if self.st[j] < 0 {
                // The reference accumulates row-wise (`0x136b00`): rows in ascending order, skipping |y_i| <= 1e-13 (`ctx+0x3a8`).
                let mut e: Vec<(usize, f64)> = p.fcols[j].clone();
                e.sort_by_key(|x| x.0);
                let mut acc = R::c(self.cost[j]) + R::c(self.cost_lo[j]);
                for &(i, a) in &e {
                    if y[i].abs() > R::c(DROP_TOL) {
                        acc -= y[i] * R::c(a);
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
    fn pri(&self, r: usize, f: R) -> R {
        let base = f * f / self.w[r];
        match self.ptype {
            1 => base * R::c(self.p70[self.head[r]] / self.p68[r]),
            2 => base / R::c(self.p68[r]),
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
        st: &crate::lu_order::Structure<R>,
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
    fn touch_pricing(&mut self, r: usize, col: &[(usize, R)]) {
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
                keys.push((-self.pri(r, self.infeas_real(r))).to_f64());
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
        let key = |s: &Self, r: usize| (-s.pri(r, s.infeas_real(r))).to_f64();
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
        let mut maxp = R::ZERO;
        let mut inf = vec![false; self.m];
        for r in 0..self.m {
            let f = self.infeas(r);
            if f > ptol {
                inf[r] = true;
                self.prio[r] = self.pri(r, self.infeas_real(r));
                maxp = maxr(maxp, self.prio[r]);
            }
        }
        let count_now = inf.iter().filter(|&&f| f).count();
        let count_cmp = if self.stale_count {
            self.inf_count
        } else {
            count_now
        };
        self.stale_count = false;
        self.inf_count = count_now;
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
        // New threshold (`0xc1e92..0xc1ee9`): max priority / 1e5, or / 1e3 when neither the old nor the new basic values are huge
        // and the list is not shorter than the pivot column's entry count (always, after a no-row call); an existing threshold below the maximum
        // priority is only ever raised.
        let thr_old = self.thr;
        // The count compared with the pivot column's entries is the number of infeasible rows (`[P+0x3c]`, `0xc1e8e`).
        let count = count_cmp;
        self.n_infeas = count_now;
        let div = if after_pivot && count < self.last_cap || big {
            1e5
        } else {
            1e3
        };
        let new = maxp / R::c(div);
        self.thr = if thr_old != R::ZERO && maxp > thr_old {
            maxr(new, thr_old)
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

    /// The acceptance test of the recompute (`0x21a383..0x21a4b6`): with `v` the infeasibility of a row beyond the tolerance, `rel` the largest `v / (|x| + 1)`, `a` the largest `v`
    /// and `xm` the largest `|x|`, the comparison threshold is the greater of the recompute's relative residual and `1e-4` times the primal tolerance.
    fn tiny_infeasibility(&self) -> bool {
        let tol = self.ptol();
        let (mut rel, mut a, mut xm) = (R::ZERO, R::ZERO, R::ZERO);
        let mut any = false;
        for r in 0..self.m {
            let x = self.xb[r];
            xm = xm.max(x.abs());
            let f = self.infeas_real(r);
            if f > R::c(tol) {
                any = true;
                a = a.max(f);
                rel = rel.max(f / (x.abs() + R::ONE));
            }
        }
        let s3 = self.xb_residual.max(R::c(tol * 1e-4));
        any && (rel < s3 || (R::ONE / R::c(10.0)) * xm * s3 > a)
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
                // Work of the pricing-side walk rebuild (`0xbdc10`): classify
                // all rows, retain waiting entries, and merge the newly infeasible rows.
                let ptol = self.ptol();
                let infeas = (0..self.m).filter(|&i| self.infeas(i) > ptol).count();
                let old_len = self.walk_len.min(self.list.len());
                let old = self.list[..old_len]
                    .iter()
                    .filter(|&&i| self.flags[i] == 0 && self.infeas(i) > ptol)
                    .count();
                let fresh = (0..self.m)
                    .filter(|&i| self.flags[i] != 0 && self.infeas(i) > ptol)
                    .count();
                let m = self.m as f64;
                let scan = if self.level == 0 {
                    3.0 * m
                } else {
                    3.0 * (self.m - infeas) as f64
                        + (5.0 + if self.ptype == 1 { 1.0 } else { 0.0 }) * infeas as f64
                };
                self.work.add(0xbe039, scan, crate::work::CTX);
                self.work
                    .add(0xbe0a1, (old_len + 2 * old) as f64, crate::work::CTX);
                self.work.add(0xbe17e, (2 * fresh) as f64, crate::work::CTX);
                self.work
                    .add(0xbe34d, (old + fresh) as f64, crate::work::CTX);
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
        let mut best = R::c(0.5) * self.thr;
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
            self.thr /= R::c(10.0);
            let mut b = R::ZERO;
            let before = self.list.len();
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
            // The refresh scan books two units per row and one per newly listed row,
            // before the fallback that can add the argmax (`0xcd18c`).
            self.work.add(
                0xcd191,
                2.0 * self.m as f64 + (self.list.len() - before) as f64,
                crate::work::CTX,
            );
            if self.list.is_empty() {
                if let Some(r) = row {
                    self.list.push(r);
                    self.flags[r] = 1;
                    self.thr = b;
                } else {
                    self.thr = R::c(1e-13);
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
                self.prio[r] = self.pri(r, self.infeas_real(r));
                match self.flags[r] {
                    3 => self.flags[r] = 1,
                    1 | 5 | 7 => {}
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
    fn pivot_row(&self, rho: &[(usize, R)], leave: usize, above: bool) -> (Vec<(usize, R)>, usize) {
        let (n, m) = (self.n, self.m);
        let p = self.p;
        // The row product runs over the rows of `rho` unless their entries (`0x134aa0` counts them) outnumber the structurals n by 1.8 (`0xc5355`, `ctx+0x68`);
        // then it is columnwise, in index order.
        let est: usize = rho.iter().map(|e| self.rrows[e.0].len()).sum();
        let mut acc = vec![R::ZERO; n + m];
        let mut touched = vec![false; n + m];
        let mut order: Vec<usize> = Vec::new();
        let sparse = rho.len() == 1 || !((n as f64) < 1.8 * est as f64);
        for &(i, rv) in rho {
            for &j in &self.rrows[i] {
                let a = p.cols[j].iter().find(|e| e.0 == i).map_or(0.0, |e| e.1);
                acc[j] = acc[j] + rv * R::c(a);
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
        // A single entry in rho takes the scaled-copy path of `0x134d90` (`0x1360ac`), which has no tolerance test for a scale of 1 or more (fx_e2_r25 call 1: entries of 6e-14 stay).
        // The copy keeps its entries only for a scale of magnitude 1 or more; below 1 the same tolerance filter as the general path follows (`0x136629`: `|scale| < 1.0`; seed 50165 call 5: rho -5e-5, a product of 1.25e-17 is dropped).
        if rho.len() != 1 || rho[0].1.abs() < R::ONE {
            order.retain(|&j| acc[j].abs() > R::c(DROP_TOL));
        }
        let nst = order.len();
        for &(i, rv) in rho {
            let j = n + i;
            if self.slack_nb[i] || j == leave {
                acc[j] = rv * R::c(p.cols[j][0].1);
                order.push(j);
            }
        }
        let sdir = if above { -1.0 } else { 1.0 };
        let out = order
            .into_iter()
            .map(|j| {
                let s = if self.st[j] == UB { -1.0 } else { 1.0 };
                (j, acc[j] * R::c(s * sdir))
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
        row: &[(usize, R)],
        leave: usize,
        above: bool,
    ) -> (Vec<usize>, Vec<R>, R, f64, bool) {
        let p = self.p;
        let sdir = if above { -1.0 } else { 1.0 };
        let mut list: Vec<(usize, R)> = Vec::new();
        for &(j, a) in row {
            if self.dflag[j] {
                continue;
            }
            let a = if a.abs() < R::c(ALPHA_TOL) {
                R::ZERO
            } else {
                a
            };
            let s = if self.st[j] == UB { -1.0 } else { 1.0 };
            list.push((j, a * R::c(s * sdir)));
        }
        if leave < self.n && !self.dflag[leave] {
            list.push((leave, R::ONE));
        }
        let w = rho_norm(&list);
        let mut dense = vec![R::ZERO; self.m];
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
                dense[i] += R::c(p.fcols[j][0].1) * v;
            } else {
                work += 2.0 * p.fcols[j].len() as f64;
                for &(i, a) in &p.fcols[j] {
                    if !seen[i] {
                        seen[i] = true;
                        order.push(i);
                    }
                    dense[i] = R::c(a) * v + dense[i];
                }
            }
        }
        // The quad update (`0x25d220`, `0x25d510`) builds the vector as a sparse list (entries above 1e-13 in magnitude) only when
        // `n m > 3 count nnz(A)` (`count` the list length, `n` = `ctx+0x68`, one above rose's structural count: 28 for the 27 columns of seed 2772); otherwise it
        // is a dense vector (count -1) whose entries are all kept, however small (`0x25d7ac`; seed 2772 call 93: a row of 1.9e-14).
        let dense_route = R::QUAD
            && ((self.n + 1) * self.m) as f64 <= list.len() as f64 * 3.0 * self.nnz_a as f64;
        let mut rows = Vec::new();
        let mut vals = Vec::new();
        let order_len = order.len();
        for i in order {
            if dense[i].abs() > R::c(1e-13) || (dense_route && dense[i] != R::ZERO) {
                rows.push(i);
                vals.push(dense[i]);
            }
        }
        work += 2.0 * order_len as f64 + 2.0 * rows.len() as f64;
        (rows, vals, w, work, dense_route)
    }

    /// Work of the row product (`0x134d90`, site `0x135e92`, same estimate as the primal path's `row_product_work`).
    fn row_product_work(&mut self, rho: &[(usize, R)], nslack: usize) {
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
            if rho[0].1.abs() < R::ONE {
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
        row: &[(usize, R)],
        esi: usize,
        leave: usize,
        above: bool,
        xbr: R,
    ) -> Ratio {
        let p = self.p;
        let delta = R::c(dual_tol());
        let (vlb, vub) = (R::c(p.lo[leave]), R::c(p.up[leave]));
        let xmm10 = if above { vub - xbr } else { xbr - vlb };
        let ptol = R::c(self.ptol());
        let xmm12 = R::c(-2.0) * xmm10;
        let xmm7 = ptol - xmm10;
        let cnt = row.len();
        let mut bound = R::c(-1e100);
        let (mut best, mut best_i, mut bestc) = (R::c(-1e100), -1i64, R::ZERO);
        let (mut xm13, mut xm9) = (R::ZERO, R::ZERO);
        let mut ratio_arr: Vec<R> = Vec::new();
        let mut con: Vec<R> = Vec::new();
        let mut lst: Vec<usize> = Vec::new();
        macro_rules! track {
            ($i:expr, $ratio:expr, $c:expr) => {{
                let (ratio, c) = ($ratio, $c);
                let k = lst.len() as i64;
                if ratio > R::c(1e-10) + best {
                    best_i = k;
                    best = ratio;
                    bestc = c;
                } else if ratio >= best - R::c(1e-10) {
                    let old = best;
                    bestc = minr(bestc, c);
                    if ratio > old {
                        best_i = k;
                    }
                    best = maxr(ratio, old);
                }
                con.push(c);
                lst.push($i);
                ratio_arr.push(ratio);
            }};
        }
        for i in esi..cnt {
            let (j, v) = row[i];
            if v < R::c(-ALPHA_TOL) && (self.st[j] == LB || self.st[j] == UB) {
                let dj = self.d[j];
                let ratio = dj / v;
                if bound > ratio {
                    xm9 = minr(v, xm9);
                } else {
                    xm13 = minr(v, xm13);
                    if dj >= R::ZERO {
                        bound = maxr(PIVOT_DELTA(delta, v, ratio), bound);
                    }
                    track!(i, ratio, xmm12);
                }
            }
        }
        let variant = self.level != 0 && self.iters + self.it_base > self.m;
        for i in 0..esi.min(cnt) {
            let (j, v) = row[i];
            if v < R::c(-ALPHA_TOL) {
                let dj = self.d[j];
                let ratio = dj / v;
                if bound > ratio {
                    xm9 = minr(v, xm9);
                    continue;
                }
                xm13 = minr(v, xm13);
                let c;
                if variant && self.st[j] == LB && p.up[j] >= BIG {
                    if dj >= R::ZERO {
                        bound = maxr(PIVOT_DELTA(delta, v, ratio), bound);
                    }
                    c = xmm12;
                } else {
                    let c0 = R::c(p.lo[j] - p.up[j]) * v;
                    if c0 > xmm7 {
                        if dj >= R::ZERO {
                            bound = maxr(PIVOT_DELTA(delta, v, ratio), bound);
                        }
                        c = minr(xmm12, c0);
                    } else {
                        c = c0;
                    }
                }
                track!(i, ratio, c);
            }
        }
        let nn = lst.len();
        self.work.add(0xc7bcd, 3.0 * cnt as f64, crate::work::CTX);
        self.retry_end = 0;
        if self.retry && self.level == 0 && R::c(PIVOT_FLOOR / 10.0) > -minr(xm13, xm9) {
            // The retried row whose largest candidate is below a tenth of the pivot floor (`0xcc3bc..0xcc583`, quad `0x223f49..0x224183`; the factor is fresh and
            // `ctx+0x3f8` is 0 here): the dual ends infeasible (status 3) when the infeasibility exceeds 100 primal tolerances, after ten such rows, or when the finite
            // ranges of the candidates cannot absorb it (a tiny infeasibility under ten tolerances first raises the LU threshold to 0.5 and refactorizes);
            // otherwise a refactorization (5) is requested after the bound flips the loop compacted.
            let pt = self.ptol();
            let end = if R::c(pt * -100.0) > xmm10 || self.retry_tiny >= 10 {
                1
            } else {
                self.retry_tiny += 1;
                let mut acc = xmm10;
                let mut done = 0usize;
                for &i in &lst {
                    let (j, av) = row[i];
                    if 1e30 > p.up[j] {
                        acc = acc + (R::c(p.up[j] - p.lo[j]) * av).abs();
                        if acc >= R::c(-pt) {
                            break;
                        }
                    }
                    done += 1;
                }
                if done == lst.len() {
                    if R::c(10.0 * pt) <= -xmm10 {
                        1
                    } else if 0.3 > self.lu_thr {
                        self.lu_thr = 0.5;
                        2
                    } else {
                        1
                    }
                } else {
                    2
                }
            };
            self.retry_end = end;
            return Ratio {
                q: None,
                flips: Vec::new(),
                reject: None,
            };
        }
        if self.retry
            && R::c(PIVOT_FLOOR / 10.0) <= -minr(xm13, xm9)
            && R::c(PIVOT_FLOOR * 10.0) > -minr(xm13, xm9)
        {
            // The retry of the worst rejected row (`ctx+0x420` < 0, `0xc7bcd` -> `0xcc39e`, `0xcc59d`): a row whose largest candidate is below ten pivot floors ends the
            // pass with status -1004 and no row or entering variable, before the selection (fx_midsing_1031 and 1496: the dual hands over to the primal).
            return Ratio {
                q: None,
                flips: Vec::new(),
                reject: None,
            };
        }
        if nn == 0 {
            self.work.add(0xc7c4f, 0.0, crate::work::CTX);
            return Ratio {
                q: None,
                flips: Vec::new(),
                reject: None,
            };
        }
        // The partition passes of the selection book `3 (hi - lo + 4)` each; the store is separate (`0xcbcae`) when the final sum is below minus the primal tolerance (`xmm9 = -1e19` at a nonzero level, else `-1e-6`, against `xmm10` at `0xcbc8d`).
        let mut sel = 0.0f64;
        let mut sel_sep = false;
        // The candidates cannot absorb the infeasibility (`0xcbcb3..0xcbcc8`: the sum of the contributions with the last one added is below minus the primal
        // tolerance, `xmm9`; only at level 0: `0xcbcc8` tests `ctx+0x2cc`, the level != 0 path was not read): no entering variable, as for an empty list (row of a primal infeasible model; found on random LPs, e.g. 3 rows and one structural
        // candidate with range * |a| = 2000 against an infeasibility of 2003).
        if nn == 1
            && bestc + xmm10 < R::ZERO
            && R::c(-self.ptol()) > bestc + xmm10
            && (self.level == 0 || self.updates > 0)
            // The same level-nonzero re-selection used for larger lists first
            // gives an artificial upper-bound candidate 2 * infeasibility
            // (`0xcbcc8`, quad `0x223793`), so it can absorb the row.
            && !(self.level != 0 && self.st[row[lst[0]].0] == UB && p.up[row[lst[0]].0] >= BIG)
        {
            self.work.add(0xc7c4f, 5.0, crate::work::CTX);
            return self.no_absorb();
        }
        if nn == 1 {
            self.work.add(0xc7c4f, 5.0, crate::work::CTX);
            let (q, _) = self.shift_entering(row, &[lst[0]], lst[0], R::ZERO, R::ZERO, 0);
            return Ratio {
                q: Some(q),
                flips: Vec::new(),
                reject: None,
            };
        }
        let maxa = -minr(xm13, xm9);
        // A row whose largest entry is below the pivot tolerance (`max(|worst|, 1e-13 maxa, 3e-8)`, `0xc7cd9`) is rejected: no entering variable is chosen, the
        // row gets flag 5 and the pricing goes on with the next one (fx_midsing_m6 calls 10 and 13). The measure kept for the worst rejected row is
        // `0.01 maxa` when the infeasibility is at least 0.01, else `infeasibility * maxa`, negated.
        if !self.retry
            && maxr(
                R::c(self.worst.abs()),
                maxr(R::c(1e-13) * maxa, R::c(PIVOT_FLOOR)),
            ) > maxa
        {
            let m = if R::c(0.01) <= -xmm10 {
                maxa * R::c(0.01)
            } else {
                xmm10 * maxa
            };
            return Ratio {
                q: None,
                flips: Vec::new(),
                reject: Some(-m.abs().to_f64()),
            };
        }
        let mut nflip = 0usize;
        if bestc + xmm10 < R::ZERO {
            let mut second = false;
            loop {
                let mut s = xmm10;
                sel = 0.0;
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
                    if s < R::ZERO {
                        s = accumulate(s, &con[lo..=h], lo % 2, false);
                    } else {
                        s = accumulate(s, &con[h + 1..=hi], (h + 1) % 2, true);
                    }
                    if R::ZERO > s {
                        lo = h + 1;
                    } else {
                        hi = h;
                    }
                }
                nflip = lo;
                sel_sep = R::c(-self.ptol()) > s;
                if self.level == 0 && R::c(-self.ptol()) > s + con[lo] {
                    self.work.add(0xc7c4f, 5.0 * nn as f64, crate::work::CTX);
                    return self.no_absorb();
                }
                // Level != 0 and the candidates cannot absorb the infeasibility (`0xcbcc8..0xcbdd4`, first pass only): a candidate at its upper bound
                // (status -2) whose bound is infinite takes the contribution `2 infeas` (the collection does it for a status -1 one only, and only from
                // iteration m on) and the selection runs again from the start on the arrays as the first pass left them (`0xc96c5`; random LP seed 305,
                // call 1: three such candidates against an infeasibility of 6e36).
                if self.level != 0 && !second && sel_sep && R::c(-self.ptol()) > s + con[lo] {
                    let mut changed = false;
                    for t in 0..nn {
                        let j = row[lst[t]].0;
                        if self.st[j] == UB && p.up[j] >= BIG {
                            con[t] = xmm12;
                            changed = true;
                        }
                    }
                    if changed {
                        // The first pass's partition work is booked at its own site (`0xcbcae`).
                        self.work.add(0xcbcae, 3.0 * sel, crate::work::CTX);
                        second = true;
                        continue;
                    }
                }
                if self.level != 0
                    && self.updates > 0
                    && sel_sep
                    && R::c(-self.ptol()) > s + con[lo]
                {
                    self.work.add(0xc7c4f, 5.0 * nn as f64, crate::work::CTX);
                    return self.no_absorb();
                }
                break;
            }
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
        // The pivot floor is `max(|ctx+0x438|, 1e-13 maxa, 3e-8)`: the measure of the worst rejected row counts (random LP seed 278: 0.0667 left by an earlier rejection).
        let pivtol = maxr(
            R::c(self.worst.abs()),
            maxr(R::c(1e-13) * maxa, R::c(PIVOT_FLOOR)),
        );
        if !(pivtol <= maxa) {
            return Ratio {
                q: None,
                flips: Vec::new(),
                reject: None,
            };
        }
        let mut b = nflip;
        let mut kept: Vec<usize> = Vec::new();
        let mut xm3 = R::c(-1e100);
        let mut xm15 = bound;
        for t in b..nn {
            let e = lst[t];
            let (j, av) = row[e];
            let dj = self.shift_below_delta(j, delta);
            if pivtol > -av {
                continue;
            }
            if xm15 * av > dj {
                let r = (delta + dj) / av;
                xm3 = maxr(r, xm3);
                kept.push(e);
                if dj >= R::ZERO {
                    xm15 = maxr(r, xm15);
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
                let dj = self.shift_below_delta(j, delta);
                if pivtol > -av {
                    continue;
                }
                if xm15 * av > dj {
                    let r = (delta + dj) / av;
                    xm3 = maxr(r, xm3);
                    kept.push(e);
                    if dj >= R::ZERO {
                        xm15 = maxr(r, xm15);
                    }
                }
            }
            b = 0;
        }
        if kept.is_empty() {
            // Relaxed fallback (`0xc8585..0xc8887`): the whole row is scanned again, without the `bound > ratio` filter. Every entry with `a < -1e-10` (a slack only
            // at a bound) joins the list; `a <= -pivtol` entries with `d >= -delta` are kept when `r = (delta + d)/a` exceeds the running threshold, which starts at
            // `bound - 100 delta / pivtol` and follows the kept `r`. `ctx+0x338` becomes 2. An entry with `d < -delta` takes the cost-shift path (`0xcb0be`: `shift_cost`; a free
            // variable's d is zeroed).
            self.pert = 2;
            let mut thr = xm15 - R::c(100.0) * delta / pivtol;
            let mut xm10 = R::c(-1e100);
            lst.clear();
            for (i, &(j, av)) in row.iter().enumerate() {
                if !(R::c(-1e-10) > av) {
                    continue;
                }
                if i >= esi && self.st[j] != LB && self.st[j] != UB {
                    continue;
                }
                lst.push(i);
                if av > -pivtol {
                    continue;
                }
                let mut dj = self.d[j];
                if -delta > dj {
                    if self.st[j] == LB || self.st[j] == UB {
                        // The cost shift of `0xcb0be` (the formulas of `shift_cost`): the entry's d becomes positive and joins the list.
                        self.shift_cost(j);
                        dj = self.d[j];
                    } else if self.st[j] == FREE {
                        // A free variable's d is zeroed with the perturbation flag raised (`0xc8739`).
                        self.raise_pert(if delta * R::c(0.0005) > dj.abs() {
                            1
                        } else {
                            2
                        });
                        self.d[j] = R::ZERO;
                        dj = R::ZERO;
                    }
                }
                let r = (delta + dj) / av;
                if r > thr {
                    xm10 = maxr(xm10, r);
                    thr = r;
                    kept.push(i);
                }
            }
            if !kept.is_empty() {
                xm3 = xm10;
                xm15 = thr;
            }
        }
        if kept.is_empty() && self.retry && !lst.is_empty() {
            // The retry of the worst rejected row (`ctx+0x420` < 0, `0x320cb6`; the quad copy of the ratio test only, the double copy rejects: fx_midsing_1031 call 26) chooses among the whole list of the fallback instead of rejecting: the
            // entries whose ratio is not below the largest `(delta + d)/a` are admissible (seed 2238: the entry of 1e-10 wins over the one of 1e-3).
            let mut mx = R::c(-1e100);
            for &i in &lst {
                let (j, av) = row[i];
                mx = maxr(mx, (delta + self.d[j]) / av);
            }
            xm3 = mx;
            xm15 = mx;
            kept = lst.clone();
        }
        if kept.is_empty() {
            // Nothing is kept after the relaxed fallback: with the rejected count `ctx+0x420` not negative (`0xcb094`: jump to the rejection `0xc7cfd`) the row is
            // rejected like a row with a small pivot row (flag 5, measure of the same formula) and the pricing goes on; the dual is unbounded only in the retry of the
            // worst rejected row (`ctx+0x420` < 0), which is not modelled.
            let m = if R::c(0.01) <= -xmm10 {
                maxa * R::c(0.01)
            } else {
                xmm10 * maxa
            };
            return Ratio {
                q: None,
                flips: Vec::new(),
                reject: Some(-m.abs().to_f64()),
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
            reject: None,
        }
    }
}

impl<'a, R: Real> Dual<'a, R> {
    /// Choice of the entering variable among the Harris-filtered candidates (`0xc89b3`): the smallest `a / w`, `w` the column count of a structural, 0.5 of a slack.
    fn select_entering(&self, row: &[(usize, R)], kept: &[usize], xm3: R) -> (usize, R) {
        let mut chosen = kept[0];
        let mut xm5 = R::ZERO;
        for &e in kept {
            let (j, av) = row[e];
            let wgt = if j < self.n {
                R::from_usize(self.p.cols[j].len())
            } else {
                R::c(0.5)
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
        let v = lcg_int(&mut self.rng);
        RNG_CELL.with(|c| c.set(self.rng));
        v
    }

    fn rng_double(&mut self) -> f64 {
        self.rng_int() as f64 * 4.656612873077393e-10
    }

    /// The quad instance adds a shift to its long double shift array; the cost in force is the double base plus that array, rounded once to long double.
    fn add_shift_q(&mut self, j: usize, delta: R) {
        self.shift_r[j] = self.shift_r[j] + delta;
        let total = R::c(self.pcost[j]) + self.shift_r[j];
        self.cost[j] = total.to_f64();
        self.cost_lo[j] = (total - R::c(self.cost[j])).to_f64();
        self.cshift[j] = self.shift_r[j].to_f64();
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
    /// The Harris filter's read of `d_j` (`0xc8138..0xc8194`): a candidate with `d_j < -delta` goes through `0xc4ca0` first (before the pivot-floor skip) and the
    /// filter uses the new value. A variable at a bound gets its cost shifted (`shift_cost`), a free one has `d` zeroed with the perturbation flag raised.
    /// Found on random near-dependent LPs (`fx_dual_dneg_1934`: a candidate with `a = -1.4e-10`, `d = -3.3e-5`).
    fn shift_below_delta(&mut self, j: usize, delta: R) -> R {
        if -delta > self.d[j] {
            if self.st[j] == LB || self.st[j] == UB {
                self.shift_cost(j);
            } else if self.st[j] == FREE {
                let a = self.d[j].abs();
                if a > R::c(1e-13) {
                    self.raise_pert(if delta * R::c(0.0005) > a { 1 } else { 2 });
                }
                self.d[j] = R::ZERO;
            }
        }
        self.d[j]
    }

    fn shift_cost(&mut self, j: usize) {
        let delta = R::c(dual_tol());
        let m = -self.d[j];
        let (lo, up) = (self.p.lo[j], self.p.up[j]);
        let small = delta * R::c(0.0005) > m
            || (delta * R::c(0.05) > m && delta * R::c(0.1) > R::c(up - lo) * m);
        let s = if small {
            self.raise_pert(1);
            let u = self.rng_double();
            m * R::c(u + 1.0)
        } else {
            self.raise_pert(2);
            let r = self.rng_int() as i32;
            m * R::c((((r % 500) + 111) as f64) / 1000.0 + 1.0)
        };
        self.d[j] += s;
        let sr = s;
        let s = s.to_f64();
        if self.st[j] == LB {
            if R::QUAD {
                self.add_shift_q(j, sr);
            } else {
                self.cshift[j] += s;
                self.cost[j] = self.pcost[j] + self.cshift[j];
            }
            self.obj_off -= lo * s;
            if -BIG >= lo {
                self.raise_pert(2);
            }
        } else {
            if R::QUAD {
                self.add_shift_q(j, -sr);
            } else {
                self.cshift[j] -= s;
                self.cost[j] = self.pcost[j] + self.cshift[j];
            }
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
        row: &[(usize, R)],
        kept: &[usize],
        chosen: usize,
        xm5: R,
        mut xm15: R,
        nflips: usize,
    ) -> (usize, R) {
        let mut xm3;
        let mut xm5 = xm5;
        let delta = R::c(dual_tol());
        let mut chosen = chosen;
        let mut kept: Vec<usize> = kept.to_vec();
        for _ in 0..8 {
            let j = row[chosen].0;
            let dj = self.d[j];
            if self.st[j] == LB || self.st[j] == UB {
                if dj < R::c(-1e-13) {
                    self.shift_cost(j);
                    if self.d[j] == R::ZERO {
                        break;
                    }
                } else {
                    if dj < R::ZERO {
                        self.d[j] = R::ZERO;
                    }
                    break;
                }
            } else if self.st[j] == FREE {
                if dj.abs() > R::c(1e-13) {
                    self.raise_pert(if delta * R::c(0.0005) > dj.abs() {
                        1
                    } else {
                        2
                    });
                }
                self.d[j] = R::ZERO;
                break;
            } else {
                break;
            }
            xm3 = xm15;
            let mut next: Vec<usize> = Vec::new();
            for &e in &kept {
                let (jj, av) = row[e];
                if self.d[jj] < R::ZERO {
                    if self.st[jj] == LB || self.st[jj] == UB {
                        if self.d[jj] < R::c(-1e-13) {
                            self.shift_cost(jj);
                        } else {
                            self.d[jj] = R::ZERO;
                        }
                    } else if self.st[jj] == FREE {
                        if self.d[jj].abs() > R::c(1e-13) {
                            self.raise_pert(if delta * R::c(0.0005) > self.d[jj].abs() {
                                1
                            } else {
                                2
                            });
                        }
                        self.d[jj] = R::ZERO;
                    }
                }
                let dj2 = self.d[jj];
                if xm15 * av > dj2 {
                    let r = (delta + dj2) / av;
                    xm3 = maxr(r, xm3);
                    next.push(e);
                    if dj2 >= R::ZERO {
                        xm15 = maxr(r, xm15);
                    }
                }
            }
            if next.is_empty() {
                // The empty-list branch combines the re-pass and selection work in one store (`0xc8b62`).
                self.work.add(
                    0xc8b67,
                    (3 * kept.len() + 4 * nflips) as f64,
                    crate::work::CTX,
                );
                chosen = kept[0];
                xm5 = R::ZERO;
            } else {
                // The re-pass scan and compaction book separately before the next selection (`0xc8e73`).
                self.work.add(
                    0xc8e78,
                    (3 * kept.len() + next.len()) as f64,
                    crate::work::CTX,
                );
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
        row: &[(usize, R)],
        chosen: usize,
        flips: Vec<usize>,
        xm5: R,
    ) -> (usize, Vec<usize>) {
        let delta = R::c(dual_tol());
        let (qj, aq) = row[chosen];
        let theta = self.d[qj] / aq;
        let tol = R::c(1000.0) * delta;
        let mut xm5 = xm5;
        let mut best: Option<usize> = None;
        // The flip-selection scan books four units per original flip (`0xc92b8`).
        self.work
            .add(0xc92bd, 4.0 * flips.len() as f64, crate::work::CTX);
        for (k, &e) in flips.iter().enumerate() {
            let (j, aj) = row[e];
            let dj = self.d[j];
            if R::ZERO > dj {
                continue;
            }
            let t = theta * aj;
            if dj > t - tol {
                continue;
            }
            // The power term only while `ctx+0x178` < 10 (quad copy `0x22159e`: past it the weight is the column count alone).
            let pw = if C178.with(|c| c.get()) >= 10 {
                R::ZERO
            } else {
                R::c(((t - dj) / delta).to_f64().powf(0.25))
            };
            let w0 = if j < self.n {
                R::from_usize(self.p.cols[j].len()) + pw
            } else {
                pw + R::c(0.5)
            };
            let w = w0 * R::c(1.0 * 0.3 + 0.7);
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
        // A replacement also books the scan and the rebuilt flip count (`0xc9431`).
        self.work.add(
            0xc9436,
            (4 * flips.len() + kept.len()) as f64,
            crate::work::CTX,
        );
        (qe, kept)
    }
}

/// The running sum of the ratio test's selection (`0xcba58..0xcbc3a`): `s` plus (or minus) the entries of `c` in the reference's vectorized order. Fewer
/// than `peel + 8` entries are added one by one. Otherwise `peel` entries (one when the first is at an odd index of the 16-byte aligned array, so that the
/// pairs are aligned) go in one by one, then blocks of eight enter four two-lane accumulators (lane 0 of the first one starts at `s`), which are combined
/// as `(a0 + a1) + (a2 + a3)` and then lane 0 + lane 1, and the remaining entries follow one by one.
fn accumulate<R: Real>(s: R, c: &[R], peel: usize, sub: bool) -> R {
    let step = |a: R, x: R| if sub { a - x } else { a + x };
    let n = c.len();
    let mut s = s;
    if !R::LANES || n < peel + 8 {
        for &x in c {
            s = step(s, x);
        }
        return s;
    }
    for &x in &c[..peel] {
        s = step(s, x);
    }
    let vend = n - ((n - peel) & 7);
    let mut acc = [[R::ZERO; 2]; 4];
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
fn PIVOT_DELTA<R: Real>(delta: R, v: R, ratio: R) -> R {
    delta / v + ratio
}

thread_local! {
    /// The leaving variable of the last ratio test without an entering variable and the number of non-slack columns (read by `scale::solve`).
    pub(crate) static FARKAS_LEAVE: std::cell::Cell<Option<(usize, usize)>> = const { std::cell::Cell::new(None) };
}

/// Result of one call of the iteration function.
enum Step<R> {
    /// A basis change: leaving row, leaving and entering variable, the pivot row (reference frame) and the data of the post-iteration passes.
    Pivot(Box<Pivot<R>>),
    /// No leaving row.
    NoRow,
    /// The row was rejected by the ratio test (pivot row too small); the pricing goes on with another row, no iteration is done.
    Reject,
    /// The pivot check rejected the pivot (row and column values disagree): no iteration; the refactorization in `request` (or the status) follows.
    Mismatch,
    /// The pivot check on a fresh factor ended the pass with status -1002: the controller restarts the dual (`0xc5b9a`).
    Restart,
    /// The pivot check on a fresh factor ended the pass with status -1005 (`0xc5f43`, "Switch to quad due to numerics"): the controller reruns the dual in
    /// quad precision from the current basis (`0xcec20`).
    SwitchQuad,
    /// The ratio test found no entering variable: the dual is unbounded (the primal infeasible). Carries the BTRAN row `e_r' B^-1` (dense, by
    /// row), which is the Farkas certificate.
    Unbounded(Vec<f64>),
}

struct Pivot<R> {
    r: usize,
    leave: usize,
    q: usize,
    row: Vec<(usize, R)>,
    /// Flipped variables, with their statuses before the flip.
    flips: Vec<usize>,
    /// The pivot column without the pivot entry, tau at the same positions, the pivot alpha and tau_r, true frame.
    col: Vec<(usize, R)>,
    tau: Vec<R>,
    alpha: R,
    w: R,
    /// Value of the entering variable before the basis change (`ctx+0x318`).
    t0: R,
    above: bool,
}

impl<'a, R: Real> Dual<'a, R> {
    /// One call of the iteration function `0xc4fe0`.
    fn iteration(&mut self, k: usize, dump: &mut Option<serde_json::Value>) -> Step<R> {
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
        let mut row_sel = self.price();
        // No row while rows were rejected before this driver call: they are released and the pricing runs again inside the same call (the loop of
        // `0xc4fe0`/`0x21cf50` around `0x2896a3` and `0xcd250`, quad seed 2081 call 19).
        while R::QUAD
            && row_sel.is_none()
            && self.request == 0
            && !self.bad_rows.is_empty()
            && self.value_mode
            && self.old_rej > 0
        {
            self.release_rejected();
            row_sel = self.price();
        }
        if row_sel.is_none()
            && self.request == 0
            && !self.bad_rows.is_empty()
            && self.value_mode
            && self.old_rej == 0
            && self.updates == 0
            && self.worst_row >= 0
        {
            // No row to price while rows were rejected, on a fresh factor and without rows released: the worst rejected row is tried again with the
            // rejected count negated (`0xc50fb`, `0x21d0xx`), its measure cleared; the ratio test then rejects nothing.
            let r = self.worst_row as usize;
            self.retry = true;
            self.worst = 0.0;
            self.flags[r] = 1;
            self.list.push(r);
            row_sel = Some(r);
        }
        let Some(r) = row_sel else {
            if let Some(mut s) = snap {
                s["r"] = serde_json::json!(-1);
                *dump = Some(s);
            }
            // No row and rows were rejected (`0xc50d4..0xc5105`, `0xc67f6..0xc683b`): the rows rejected before this driver call are released (`0xcd250`: flag 0 and
            // a new priority if still infeasible, else flag 2) and the pricing runs again; with none of them left, a factor with updates is refactorized with
            // a higher pivot threshold (request 6), a fresh one is the end handled in `no_row` (fx_midsing_1106 call 21).
            if self.request == 0 && !self.bad_rows.is_empty() && self.value_mode {
                if self.old_rej > 0 {
                    self.release_rejected();
                    return Step::Reject;
                }
                if self.updates != 0 {
                    self.lu_thr = if 0.125 > self.lu_thr {
                        0.125
                    } else {
                        self.lu_thr.max(0.5)
                    };
                    self.request = 6;
                    return Step::Mismatch;
                }
            }
            return Step::NoRow;
        };
        let leave = self.head[r];
        let xbr = self.xb[r];
        let above = xbr > R::c(p.up[leave]);
        // Pivot-row BTRAN, drop of the small entries.
        let rho_full = self.mirror.btran(&[(r, R::ONE)], true, false);
        let rho: Vec<(usize, R)> = rho_full
            .iter()
            .copied()
            .filter(|e| e.1.abs() > R::c(DROP_TOL))
            .collect();
        self.p60 = rho.len() as f64;
        self.work
            .add_from(0x13174d, 0xc51dd, self.mirror.work, crate::work::FAC);
        let (row, esi) = self.pivot_row(&rho, leave, above);
        self.row_product_work(&rho, row.len() - esi);
        self.p58 += rho.iter().map(|e| self.rrows[e.0].len()).sum::<usize>() as f64;
        let rt = self.ratio_test(&row, esi, leave, above, xbr);
        if std::mem::take(&mut self.absorb_request) {
            if let Some(mut sn) = snap {
                sn["r"] = serde_json::json!(r);
                sn["var"] = serde_json::json!(leave);
                sn["rho"] = crate::dump::sparse_r(&rho);
                sn["row"] = crate::dump::sparse_r(&row);
                sn["q"] = serde_json::json!(-1);
                *dump = Some(sn);
            }
            return Step::Mismatch;
        }
        let was_retry = self.retry;
        if self.retry {
            self.retry = false;
            if rt.q.is_none() && self.retry_end == 2 {
                self.retry_end = 0;
                self.request = 5;
                if let Some(mut sn) = snap {
                    sn["r"] = serde_json::json!(r);
                    sn["var"] = serde_json::json!(leave);
                    sn["rho"] = crate::dump::sparse_r(&rho);
                    sn["row"] = crate::dump::sparse_r(&row);
                    sn["q"] = serde_json::json!(-1);
                    *dump = Some(sn);
                }
                return Step::Mismatch;
            }
            if rt.q.is_some() {
                // The wrapper releases after the call returns; the worst row of the call (the mismatch ladder reads it) stays.
                let (w, wr) = (self.worst, self.worst_row);
                self.release_all();
                self.worst = w;
                self.worst_row = wr;
            }
            if rt.q.is_none() && self.retry_end == 0 {
                // The retried row has no usable entering variable: the end of the rejected rows (status -1004 on the fresh factor).
                if let Some(mut s) = snap {
                    s["r"] = serde_json::json!(-1);
                    s["q"] = serde_json::json!(-1);
                    *dump = Some(s);
                }
                return Step::NoRow;
            }
            // The retried row is taken: the rejected count (`ctx+0x420`, negated for the retry) is back to 0 at the next call (seed 2717 call 18).
            self.bad_rows.clear();
            self.old_rej = 0;
        }
        let mut rec = snap;
        if let Some(s) = rec.as_mut() {
            s["r"] = serde_json::json!(r);
            s["var"] = serde_json::json!(leave);
            s["rho"] = crate::dump::sparse_r(&rho);
            s["row"] = crate::dump::sparse_r(&row);
        }
        // The row's value of the chosen entering variable is stored (`ctx+0x308`) before the single-candidate check rejects it (fx_midsing_1496 call 5: -4.5e-9).
        if let Some(qe) = rt.q {
            self.pair.1 = row[qe].1.to_f64();
        }
        // A row rejected for its small entries stores the measure of the rejection there too (`ctx+0x308` = `ctx+0x438`: fx_midsing_1106 pairs of -6.2287e-11).
        if let (Some(rj), true) = (rt.reject, true) {
            self.pair.1 = rj;
        }
        // The chosen pivot is checked once more (`0xc94a3..0xc9528`): a value above minus the pivot floor on a non-free variable rejects the row like the
        // test on the whole row (a single candidate bypasses the tolerance of the selection; rnd1496 call 5).
        let rt = match rt.q {
            Some(qe)
                if !was_retry && self.st[row[qe].0] != FREE && row[qe].1 > R::c(-PIVOT_FLOOR) =>
            {
                Ratio {
                    q: None,
                    flips: Vec::new(),
                    reject: Some(-row[qe].1.abs().to_f64()),
                }
            }
            _ => rt,
        };
        if let Some(meas) = rt.reject {
            if let Some(s) = rec.as_mut() {
                s["q"] = serde_json::json!(-1);
                s["reject"] = serde_json::json!(true);
            }
            *dump = rec;
            if self.value_mode {
                self.bad_rows.push(r);
                self.flags[r] = 5;
                if self.worst > meas {
                    self.worst = meas;
                    self.worst_row = r as i64;
                }
            } else {
                // The list-walk mode (`0xc94fb` -> `0xc9577`): the row is dropped without a flag; at the end of the walk list the factor is refactorized or the value mode starts.
                self.walk_end_request();
            }
            return Step::Reject;
        }
        let Some(qe) = rt.q else {
            if let Some(s) = rec.as_mut() {
                s["q"] = serde_json::json!(-1);
            }
            *dump = rec;
            let mut y = vec![0.0; m];
            for &(i, v) in &rho_full {
                let v = v.to_f64();
                y[i] = if above { -v } else { v };
            }
            // Roundoff in a row with an infinite bound would void the certificate: entries at noise level are dropped.
            let big = y.iter().fold(0.0f64, |a, v| a.max(v.abs()));
            y.iter_mut()
                .filter(|v| v.abs() <= 1e-11 * big)
                .for_each(|v| *v = 0.0);
            FARKAS_LEAVE.with(|c| c.set(Some((leave, self.n))));
            return Step::Unbounded(y);
        };
        let (q, aq_row) = row[qe];
        let flips: Vec<usize> = rt.flips.iter().map(|&e| row[e].0).collect();
        // The objective progress flag `ctx+0x110` (`0x21cf50` at `0x31dc95`, `0xc56f0` in the double code): cleared by a pivot below 1e-6, a leaving
        // infeasibility above 1e8, or when the gain accumulated since the recompute (times `ctx+0xf8`, 1e-13) exceeds 100 times the dual tolerance.
        {
            let infeas = if above {
                R::c(p.up[leave]) - xbr
            } else {
                xbr - R::c(p.lo[leave])
            };
            let mut ob = R::ZERO;
            let mut inf2 = infeas;
            if !rt.flips.is_empty() {
                let mut sm = R::ZERO;
                for &e in &rt.flips {
                    let (j, a) = row[e];
                    if j != q {
                        let g = R::c(p.up[j] - p.lo[j]);
                        sm = sm - a * g;
                        ob = self.d[j] * g + ob;
                    }
                }
                inf2 = infeas + sm;
            }
            let a308 = aq_row.to_f64();
            ob = ob + (self.d[q] / R::c(a308)) * inf2;
            ob = ob / R::c(self.obj_scale);
            if a308.abs() < 1e-6 || 1e8 < -infeas.to_f64() {
                self.prog_on = false;
            } else {
                let tol = 100.0 * dual_tol();
                let acc = (R::c(self.obj_scale) * ob * R::c(1e-13) + R::c(self.prog_acc)).to_f64();
                self.prog_acc = acc;
                if tol < acc {
                    self.prog_on = false;
                }
            }
        }
        if !flips.is_empty() {
            self.work
                .add(0xc561a, 7.0 * flips.len() as f64, crate::work::CTX);
        }
        if let Some(s) = rec.as_mut() {
            s["q"] = serde_json::json!(q);
            s["alpha_row"] = serde_json::json!(aq_row.to_f64());
            s["flips"] = serde_json::json!(flips);
        }
        // The reference's two-vector FTRAN: the pivot column first, then tau = B^-1 rho at the column's positions only.
        let rows_q: Vec<usize> = p.fcols[q].iter().map(|&(i, _)| i).collect();
        let vals_q: Vec<R> = p.fcols[q].iter().map(|&(_, a)| R::c(a)).collect();
        let col_out = self.mirror.ftran(&rows_q, &vals_q);
        // Mode 3: the second vector is not rho but the combination of the columns of the row's entries that were basic when the weights restarted
        // (`0xc57d8`, `0x11f620`), and its squared norm is the weight of the pivot row.
        let (rho_rows, rho_vals, w_row, vec_work, dense_route) = if self.pmode == 3 {
            self.framework_rhs(&row, leave, above)
        } else {
            (
                rho.iter().map(|e| e.0).collect::<Vec<usize>>(),
                rho.iter().map(|e| e.1).collect::<Vec<R>>(),
                rho_norm(&rho),
                rho.len() as f64,
                false,
            )
        };
        let needed: Vec<usize> = col_out.iter().map(|e| e.0).collect();
        // Mode 3: the framework vector is a dense vector (count -1) once it has more than m/8 entries (SC50A calls 11, 14..; AFIRO call 12: 4 of m = 27;
        // 6 of 50 stay sparse); its L stage then runs dense.
        let dense_in = self.pmode == 3 && (dense_route || rho_rows.len() * 8 > m);
        let tau_out = self
            .mirror
            .ftran_needed(&rho_rows, &rho_vals, &needed, dense_in);
        if R::QUAD && crate::dump::on("qtau") {
            let v: Vec<(usize, R)> = rho_rows
                .iter()
                .copied()
                .zip(rho_vals.iter().copied())
                .collect();
            crate::dump::emit_tag(
                "qtau",
                "QTAU",
                serde_json::json!({"r": r, "pmode": self.pmode, "dense_in": dense_in, "in": crate::dump::sparse_r(&v), "col": crate::dump::sparse_r(&col_out), "tau": crate::dump::sparse_r(&tau_out)}),
            );
        }
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
        let mut tau_dense = vec![R::ZERO; m];
        for &(i, v) in &tau_out {
            tau_dense[i] = v;
        }
        let mut alpha_dense = vec![R::ZERO; m];
        for &(i, v) in &col_out {
            alpha_dense[i] = v;
        }
        let alpha = alpha_dense[r];
        // Pivot check: the row's value in the true frame against the column's.
        let sq = if self.st[q] == UB { -1.0 } else { 1.0 };
        let sdir = if above { -1.0 } else { 1.0 };
        let alpha_row_true = aq_row * R::c(sq * sdir);
        self.pair.0 = alpha.to_f64();
        let rel = (alpha_row_true - alpha).abs() / (alpha.abs() + R::ONE);
        // Mismatch (`0xc5a0b..0xc5a19`, quad `0x21da1f..0x21da2d`): the two values differ in sign (`(min - 1e-13)(max + 1e-13) <= 0`) or by more than `ctx+0x3a0` (3e-6, 1e-7 from the first restore of a dropped variable on: it persists into later dual passes, seed 2948).
        let (amin, amax) = (alpha_row_true.min(alpha), alpha_row_true.max(alpha));
        let mismatch = (amin - R::c(DROP_TOL)) * (amax + R::c(DROP_TOL)) <= R::ZERO
            || rel > R::c(PAIR_TOL.with(|c| c.get()));
        if let Some(s) = rec.as_mut() {
            s["alpha_col"] = serde_json::json!(alpha.to_f64());
            s["tau_r"] = serde_json::json!(tau_dense[r].to_f64());
            s["w"] = serde_json::json!(w_row.to_f64());
            s["col"] = crate::dump::sparse_r(&col_out);
            s["tau"] = crate::dump::sparse_r(&tau_out);
        }
        if mismatch {
            C178.with(|c| c.set(c.get() + 1));
            if let Some(s) = rec.as_mut() {
                // The ladder ends the call without an entering variable (`ctx+0x2fc` = -1, `0xc5abd`).
                s["mismatch"] = serde_json::json!(true);
            }
            // The ladder of `0xc5a2f..0xc5c93`: a pivot whose row and column values disagree is rejected (no iteration) while the factorization can still be made more
            // stable: below 0.5 the pivot threshold becomes 0.5 and the basis is refactorized (request 6); after a restart (`ctx+0x184` != 0) below 0.9 it becomes 0.9
            // (and the retry counter grows). Past that, in the value mode: a factor with updates is refactorized, a fresh one ends the pass with status -1002 (restart,
            // counter + 1) up to five times and later accepts a pivot of agreeing sign or ends with status 12. With `Quad` at its default (-1, the release profile) and `ctx+0x28` = 0 a fresh
            // factor ends with status -1005 instead (`0xc5f43`, "Switch to quad due to numerics": the controller reruns the dual in quad precision, `0xcec20`), which is not modelled.
            // The list-walk mode is handled in the branch below.
            if self.lu_thr < 0.5 {
                self.lu_thr = 0.5;
                self.request = 6;
                *dump = rec;
                return Step::Mismatch;
            }
            if self.lu_thr < 0.9 && self.retries != 0 {
                self.retries += 1;
                self.lu_thr = 0.9;
                self.request = 6;
                *dump = rec;
                return Step::Mismatch;
            }
            if self.value_mode {
                if self.worst_row < 0 {
                    // The first rejection since the last factorization (`ctx+0x430` < 0, `0xc6007`): the row joins the rejected list with flag 5 and becomes
                    // the worst one (`ctx+0x430`, with minus the magnitude of 0.001 of its pivot row value in `ctx+0x438` when that is lower), and the pricing
                    // goes on with another row; the second one climbs the ladder (fx_midsing_1390 calls 7 and 8).
                    let meas = -(alpha_row_true * R::c(0.001)).abs().to_f64();
                    self.bad_rows.push(r);
                    self.flags[r] = 5;
                    if self.worst > meas {
                        self.worst = meas;
                        self.worst_row = r as i64;
                    }
                    if let Some(s) = rec.as_mut() {
                        s["reject"] = serde_json::json!(true);
                    }
                    *dump = rec;
                    return Step::Reject;
                }
                if self.updates != 0 {
                    self.request = 6;
                    *dump = rec;
                    return Step::Mismatch;
                }
                if !R::QUAD && self.switch_enabled {
                    *dump = rec;
                    return Step::SwitchQuad;
                }
                if self.retries <= 5 {
                    *dump = rec;
                    return Step::Restart;
                }
                if alpha * alpha_row_true <= R::ZERO {
                    self.status = Some(Status::NumericalFailure);
                    *dump = rec;
                    return Step::Mismatch;
                }
            } else {
                // The list-walk mode (`0xc5ab2` with `P+0x40` = 0, `0xc5bd0`): the row is rejected without an iteration (`ctx+0x2f8`, `0x2fc` = -1, no flag
                // set; the cursor has passed it). When the walk has reached the end of its list (`ctx+0x2c4` = `P+0x3c`): more updates than half the
                // cursor ask for a refactorization (request 6, none for a fresh factor); otherwise the pricing goes to the value mode (`P+0x40` = -1) and
                // asks for one.
                self.walk_end_request();
                *dump = rec;
                return Step::Reject;
            }
        }
        if let Some(o) = rec.as_mut().and_then(|s| s.as_object_mut()) {
            // A mismatch the ladder accepted (agreeing signs past the sixth retry) is an iteration: no mismatch mark in the record.
            o.remove("mismatch");
        }
        // Column without the pivot entry: the last entry moves into the pivot's slot.
        // A dense FTRAN result (count -1) is compacted by an ascending scan that skips the pivot and the entries up to the drop tolerance (`0x12013c`).
        let mut col: Vec<(usize, R)> = col_out.clone();
        if self.mirror.last_dense() {
            col.retain(|e| e.0 != r && e.1.abs() > R::c(DROP_TOL));
        } else if let Some(k) = col.iter().position(|e| e.0 == r) {
            col.swap_remove(k);
        }
        let tau: Vec<R> = col.iter().map(|e| tau_dense[e.0]).collect();
        // Bound flips: x_B and the statuses (`0xcd720`).
        if !flips.is_empty() {
            let mut rhs = vec![R::ZERO; m];
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
                    rhs[i] += R::c(a) * R::c(delta);
                }
            }
            // The reference's input vector (`0xcd720`): a dense one when twice the summed column counts exceed m, else the touched
            // rows above the drop tolerance in order of first touch (`0xcdb53..0xcdb8a`). Then the plain FTRAN of the factor.
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
                (0..m).filter(|&i| rhs[i] != R::ZERO).collect()
            } else {
                touched
                    .into_iter()
                    .filter(|&i| rhs[i].abs() > R::c(DROP_TOL))
                    .collect()
            };
            let vals: Vec<R> = rows.iter().map(|&i| rhs[i]).collect();
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
                if v != R::ZERO {
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
        let t0 = R::c(self.nonbasic_value(q));
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
        if self.d[q] <= R::c(dual_tol()) {
            self.stall += 1;
        }
        self.iters += 1;
        *dump = rec;
        let _ = flips.len();
        // The factor's update is made by the post-iteration pass (it refactorizes instead when one is due).
        let alpha_vec = alpha_dense;
        self.pending_alpha = Some(alpha_vec);
        // The ratio test zeroes row entries of magnitude at most 1e-10 (`0xc7747..0xc7769`, `0xcc606`; bh50022 has an entry exactly -1e-10).
        let row: Vec<(usize, R)> = row
            .into_iter()
            .map(|(j, v)| {
                (
                    j,
                    if v.abs() <= R::c(ALPHA_TOL) {
                        R::ZERO
                    } else {
                        v
                    },
                )
            })
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

impl<'a, R: Real> Dual<'a, R> {
    fn snapshot(&self, k: usize) -> serde_json::Value {
        serde_json::json!({
            "c2d8": dual_infeas_sum(), "k": k, "lvl": self.level, "req": self.request, "it": self.iters as f64, "head": self.head, "xB": self.xb.iter().map(|v| v.json()).collect::<Vec<_>>(), "d": self.d.iter().map(|v| v.json()).collect::<Vec<_>>(), "W": self.w.iter().map(|v| v.json()).collect::<Vec<_>>(),
            "st": self.st, "flags": self.flags, "prio": self.prio.iter().map(|v| v.to_f64()).collect::<Vec<_>>(), "diag": self.mirror.numerics().0.iter().map(|v| v.json()).collect::<Vec<_>>(),
            "list": if self.value_mode { self.list.clone() } else { self.list[..self.walk_len.min(self.list.len())].to_vec() },
            "thr": self.thr.to_f64(), "p40": self.value_mode as u8, "cur": self.cursor,
        })
    }

    /// The driver around the pricing returns (`0xf6745..0xf6768`, `0xf8b51..0xf8b77`): the worst rejected row is cleared and the rejected rows so far are the
    /// ones a later call may release.
    fn driver_return(&mut self) {
        self.worst = 0.0;
        self.worst_row = -1;
        self.old_rej = self.bad_rows.len();
    }

    /// `0xcd250` on the rows rejected before the current driver call: a row still infeasible gets flag 0 and its priority, any other flag 2; the rest of the
    /// list moves up.
    fn release_rejected(&mut self) {
        let tol = self.ptol();
        let first: Vec<usize> = self
            .bad_rows
            .drain(..self.old_rej.min(self.bad_rows.len()))
            .collect();
        self.old_rej = 0;
        for r in first {
            if self.flags[r] < 5 {
                continue;
            }
            let v = self.head[r];
            let (x, lo, up) = (self.xb[r].to_f64(), self.p.lo[v], self.p.up[v]);
            if (lo > -BIG && lo - tol > x) || (up < BIG && x > tol + up) {
                self.prio[r] = self.pri(r, self.infeas_real(r));
                self.flags[r] = 0;
            } else {
                self.flags[r] = 2;
            }
        }
    }

    /// `0x2252b0` with a count of 0 (quad; the wrapper `0x24bdc2` calls it after a retried worst row, `0x24c2fb`): the list of rejected rows is emptied, the worst row cleared, and every row
    /// with flag 5 or more is released: still infeasible it gets flag 0 and its priority, and joins the pricing list (flag 1) when the priority reaches the threshold; any other row
    /// gets flag 2.
    fn release_all(&mut self) {
        let tol = self.ptol();
        self.bad_rows.clear();
        self.old_rej = 0;
        self.worst = 0.0;
        self.worst_row = -2;
        for r in 0..self.m {
            if self.flags[r] < 5 {
                continue;
            }
            let v = self.head[r];
            let (x, lo, up) = (self.xb[r].to_f64(), self.p.lo[v], self.p.up[v]);
            if (lo > -BIG && lo - tol > x) || (up < BIG && x > tol + up) {
                self.prio[r] = self.pri(r, self.infeas_real(r));
                self.flags[r] = 0;
                if self.value_mode && self.prio[r] >= self.thr {
                    self.list.push(r);
                    self.flags[r] = 1;
                }
            } else {
                self.flags[r] = 2;
            }
        }
    }

    /// The candidates cannot absorb the infeasibility (`0xcc061`: the sum of all contributions with the infeasibility is below minus the primal tolerance): on a factor with
    /// updates the ratio test asks for a refactorization (request 5, `0xcc1d6`, `ctx+0x3f8` != 0) and returns no entering variable without ending the dual (random LP seed
    /// 31249 call 3: one candidate of range 10 against an infeasibility of 2.9e19 at level 1); on a fresh factor the dual ends as with an empty row.
    fn no_absorb(&mut self) -> Ratio {
        if self.updates > 0 {
            self.request = 5;
            self.absorb_request = true;
        }
        Ratio {
            q: None,
            flips: Vec::new(),
            reject: None,
        }
    }

    /// The end of the list walk after a rejected row (`0xc9577..0xc95b3`, `0xc5bd0`): at the end of the list (`ctx+0x2c4` = `P+0x3c`) more updates than half the
    /// cursor ask for a refactorization (request 6, none for a fresh factor); otherwise the pricing goes to the value mode (`P+0x40` = -1) and asks for one.
    fn walk_end_request(&mut self) {
        if self.cursor == self.walk_len {
            if self.updates > self.cursor / 2 {
                if self.updates != 0 {
                    self.request = 6;
                }
            } else {
                self.value_mode = true;
                self.request = 6;
            }
        }
    }

    /// The no-leaving-row handler (`0x1aa280`): a refactorization request, optimality or infeasibility.
    fn no_row(&mut self) {
        let bad_recomputed =
            R::QUAD && self.updates != 0 && self.level != 0 && self.artificial_count() != 0;
        if bad_recomputed {
            // The quad driver's post-iteration function (`0x237750`) with the confirming refactorization pending (`ctx+0x294` = 4) and updates on the
            // factor: it recomputes d from the costs with the updated factor (`0x26ce00` BTRAN, `0x273cd0`) and runs the status count (`0x238c38..0x238d1f`);
            // variables at the bound that d turns wrong, with an infinite range on the far side, end the pass with status -1004 before the refactorization.
            self.art_exit = true;
            self.dual_unbounded();
            return;
        }
        if self.updates != 0 {
            // `ctx+0x3f8 != 0` in all captures: confirm on a fresh factor (`0x1aa280`). The cost shifts stay in force through the confirming refactorization
            // (`ctx+0x290` is 0, `0xe4891` is not taken) and are removed at the end of the solve.
            self.request = 4;
        } else if !self.bad_rows.is_empty() {
            // Rows were rejected by the ratio test (pivot row too small) and the fresh factor leaves none to price: the dual ends with the cost shifts in
            // force (status -1004 in the reference); the controller hands over to the primal, which removes them, and reruns the dual (rnd1496 pass 1).
            self.status = Some(Status::Infeasible);
            self.rejected_end = true;
        } else if self.free_infeasible() {
            // A free nonbasic variable (never pivoted in: no basic variable with a finite lower bound) with a nonzero d.
            self.dual_unbounded();
        } else if self.level == 0 {
            self.status = Some(Status::Optimal);
        } else {
            let p = self.p;
            for j in 0..self.n + self.m {
                let bad = (self.st[j] == UB && p.up[j] >= BIG && self.d[j] > R::c(dual_tol()))
                    || (self.st[j] == LB && p.lo[j] <= -BIG && self.d[j] > R::c(dual_tol()));
                if bad {
                    self.dual_unbounded();
                    return;
                }
            }
            // No artificial bound is required by its d (`ctx+0x290 = 1`, read in the reference but not observed): a refactorization puts those variables
            // on their finite bound (`set_level`) and recomputes x_B without the sentinels. The quad driver ends with status 2 on the fresh factor whatever the level.
            if R::QUAD {
                self.status = Some(Status::Optimal);
            } else {
                self.request = 4;
            }
        }
    }

    /// The status count of the quad post-iteration function (`0x238c38..0x238d1f`) on d recomputed with the current factor: a variable at its lower bound
    /// (status -1) with `d < -tol` and no finite upper bound, or at its upper bound (status -2) with `d > tol` and no finite upper bound, counts.
    fn artificial_count(&mut self) -> usize {
        let y = self.last_y.clone();
        let dt = self.true_d();
        self.last_y = y;
        let tol = R::c(dual_tol());
        let mut count = 0;
        for j in 0..self.n + self.m {
            let inf = self.p.up[j] >= BIG;
            if self.st[j] == LB {
                if dt[j] < -tol && inf {
                    count += 1;
                }
            } else if self.st[j] == UB && inf && (-dt[j]) > tol {
                count += 1;
            }
        }
        count
    }

    /// Nonbasic at an artificial (infinite, 1e30) bound: at the upper bound of a variable without one, or at the lower bound of a variable without one.
    /// The recovered code of the reference counts the first kind only (`0xe5d5d`); the second occurs only with a lower bound of -1e30 and a finite upper
    /// bound, which no corpus model has.
    fn at_artificial(&self, j: usize) -> bool {
        (self.st[j] == UB && self.p.up[j] >= BIG) || (self.st[j] == LB && self.p.lo[j] <= -BIG)
    }

    /// A free nonbasic variable has a reduced cost beyond the dual tolerance.
    fn free_infeasible(&self) -> bool {
        (0..self.n + self.m).any(|j| self.st[j] == FREE && self.d[j].abs() > R::c(dual_tol()))
    }

    /// The dual is infeasible (unbounded dual ray in the primal's terms): a nonbasic variable sits at an artificial bound with the wrong-signed d, or is free
    /// with a nonzero d. The basic solution is primal feasible at the scale of the artificial bounds, so the point it describes is a finite point plus 1e30 times
    /// a recession direction. That direction moves every nonbasic variable at an artificial bound away from its finite bound, every free one with a nonzero d
    /// against its d, and the basic variables with them (the factor is fresh here). Whether the primal is feasible is for the caller to settle.
    fn dual_unbounded(&mut self) {
        // The reference reads the ray from the point the loop holds (`0x1b4870`, branch without an entering variable): the primal solution of the
        // artificial-bound problem (nonbasic variables at `+-1e30`, basics from the FTRAN) divided by `1e30` (`ctx+0x370`), entries of 1e-13 or less dropped.
        let mut ray = vec![0.0; self.n + self.m];
        for j in 0..self.n + self.m {
            let v = if self.st[j] >= 0 {
                self.xb[self.st[j] as usize]
            } else {
                R::c(self.nonbasic_value(j))
            };
            let r = (v / R::c(BIG)).to_f64();
            if r.abs() > 1e-13 {
                ray[j] = r;
            }
        }
        ray.truncate(self.nm);
        self.ray = Some(ray);
        self.status = Some(Status::Unbounded);
    }

    /// Work of `0xe45c0` after a refactorization (site `0xe7000`): two copies of x_B and d's rows (`2 m`), the status pass (`m + 2 n`), a pass over all variables
    /// when the level was not zero, and the tail that is longer for a request-4 confirmation with a nonzero incoming level (or request 6). With the perturbation flag set (`ctx+0x338`
    /// != 0, kept to the end of the solve) the request-4 tail is the plain `1.5 n` and every such factorization adds `2 m + n` (the cost recompute with the shift
    /// array, `0xe555c`; ADLITTLE, BOEING2, PILOT4).
    fn set_level_work(&mut self, level_was_set: bool) {
        let (n, m) = (self.n as f64, self.m as f64);
        let mut est = 2.0 * m + (m + 2.0 * n);
        if level_was_set {
            est += 0.5 * (n + m);
        }
        est += if ((self.request == 4 && level_was_set) || self.request == 6) && self.pert == 0 {
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

    /// The cost shift of the post-refactorization pass (`0xe45c0`, loop `0xe5f15..0xe5f3d`, branches `0xe730b`, `0xe72bb`, `0xe7263`; every request but 1 and 4):
    /// a nonbasic variable at its lower bound (status -1) with `d < -tol`, or at its upper bound (status -2, d in the negated frame) with `d < -tol`, takes the
    /// exact shift that zeroes its true d (`cost -= d` at the lower bound, `cost += d` at the upper), and its d becomes 0.0; the perturbation flag
    /// `ctx+0x338` becomes at least 2. A free nonbasic variable (status -3) with `|d| > tol` is shifted like the lower-bound case. Fixed ones (status -4) stay.
    /// In a quad pass the request-4 refactorization shifts too (the reference keeps st and sets d = 0 where the double dual moves the variable to its other bound; seed 5247).
    fn shift_wrong_d(&mut self) {
        // The double dual's confirming refactorization (request 4) shifts only with the cost shifts already in force (`ctx+0x338` != 0; seed 3981: d = -29 and -58 at
        // a lower bound with an infinite range have d = 0 at the next call, not a level of artificial bounds).
        if self.request == 1 || (self.request == 4 && !R::QUAD && self.pert == 0) {
            return;
        }
        for j in 0..self.n + self.m {
            let d = self.d[j];
            let tol = R::c(dual_tol());
            let shift = if self.st[j] == LB {
                (d < -tol).then_some(-d)
            } else if self.st[j] == UB {
                (d < -tol).then_some(d)
            } else if self.st[j] == FREE {
                (d.abs() > tol).then_some(-d)
            } else {
                None
            };
            if let Some(s) = shift {
                if R::QUAD {
                    // `0x239ff7..0x23a02e`: the shift array takes `-d` (lower bound, free) or `+d` (upper bound) in long double.
                    self.add_shift_q(j, s);
                } else {
                    let s = s.to_f64();
                    self.cshift[j] += s;
                    self.cost[j] += s;
                }
                self.d[j] = R::ZERO;
                self.raise_pert(2);
            }
        }
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
                    && self.d[j].abs() <= R::c(dual_tol())
                {
                    // A free variable no longer tied to an artificial bound by its d goes back to the free status (value 0).
                    self.st[j] = FREE;
                    self.d[j] = R::ZERO;
                } else if self.st[j] == LB && self.d[j] < R::c(-dual_tol()) {
                    self.st[j] = UB;
                    self.d[j] = -self.d[j];
                } else if self.st[j] == UB && up_inf && self.d[j] <= R::c(dual_tol()) {
                    self.st[j] = LB;
                    self.d[j] = -self.d[j];
                } else if self.st[j] == UB && !up_inf && self.d[j] < R::c(-dual_tol()) {
                    // At a finite upper bound with the wrong-signed d (`0xe7820`, status -2: d' < -tol): the variable goes to its lower bound (fx_badscale
                    // release, x with a scaled cost of -2e-9 that the pivot turns into +2e-6).
                    self.st[j] = LB;
                    self.d[j] = -self.d[j];
                } else if self.st[j] == LB && lo_inf && !up_inf && self.d[j] <= R::c(dual_tol()) {
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
                *w = R::ONE;
            }
        }
    }

    /// The weight-overflow reset of the quad driver (`0x24ac2c..0x24b51a`, after each iteration call): with a pending request other than 1, the factorization count
    /// a multiple of 3 and the pricing mode 1 or 3, when fewer than `10 + m/10` entries were visited since the pass began or the iterations exceed 3 m or n exceeds 3 m,
    /// any weight above 1e20 (1e30 outside the mode 3) or not a number restarts the weights at 1 in mode 3.
    fn weight_overflow_gate(&mut self, req: u32) {
        if !R::QUAD || req == 1 || self.nref() % 3 != 0 || !matches!(self.pmode, 1 | 3) {
            return;
        }
        let m = self.m as f64;
        let few = self.p58 < (10 + self.m / 10) as f64;
        if !(few || 3.0 * m < (self.iters + self.it_base) as f64 || 3.0 * m < self.n as f64) {
            return;
        }
        let thr = R::c(if self.pmode == 3 { 1e20 } else { 1e30 });
        if self.w.iter().any(|&w| w > thr || w != w) {
            self.pmode = 3;
            self.p58 = 0.0;
            self.dflag = self.st.iter().map(|&t| t < 0).collect();
            for w in self.w.iter_mut() {
                *w = R::ONE;
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
        if self.nref() % period != 1 {
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
        let a = dual_tol();
        let mut mx = 0.0f64;
        for j in 0..self.n + self.m {
            if self.st[j] == FREE {
                mx = maxsd(self.d[j].abs().to_f64(), mx);
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
        for v in self.cost_lo.iter_mut() {
            *v *= f;
        }
        for v in self.pcost.iter_mut() {
            *v *= f;
        }
        for v in self.d.iter_mut() {
            *v *= R::c(f);
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
        // The seed is the pass's iteration counter including the iterations before the quad switch (`ctx+0x2a8`), as the stall and the ratio test variant count.
        let mut g = (self.iters + if R::QUAD { self.it_base } else { 0 }) as u64;
        let mut acc = R::ZERO;
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
            // The quad instance keeps xi in long double (`0x23a293`: `(r mod 250 + 1) * sc` in x87, times the double log2, not rounded to a double).
            let xi = R::c((r % 250 + 1) as f64) * R::c(sc) * R::c((x + 1.0).log2());
            let xd = xi.to_f64();
            let bump = |s: &mut Self, up: bool| {
                if R::QUAD {
                    s.add_shift_q(j, if up { xi } else { -xi });
                } else {
                    s.pcost[j] += if up { xd } else { -xd };
                    s.cost[j] += if up { xd } else { -xd };
                }
            };
            if self.st[j] == LB {
                bump(self, true);
                self.d[j] += xi;
                acc -= R::c(p.lo[j]) * xi;
            } else if self.st[j] == UB {
                if p.up[j] < BIG {
                    bump(self, false);
                    self.d[j] += xi;
                    acc += R::c(p.up[j]) * xi;
                } else if self.d[j] <= R::c(dual_tol()) {
                    self.d[j] = -self.d[j];
                    bump(self, true);
                    acc -= R::c(p.lo[j]) * xi;
                    self.d[j] += xi;
                    self.st[j] = LB;
                }
            }
        }
        self.obj_off = acc.to_f64() + self.obj_off;
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
    fn post(&mut self, pv: Option<&Pivot<R>>, schedule_due: bool) -> Result<(), SolveError> {
        self.post_inner(pv, schedule_due)?;
        self.log_snapshot();
        Ok(())
    }

    fn post_inner(&mut self, pv: Option<&Pivot<R>>, schedule_due: bool) -> Result<(), SolveError> {
        let refactor = self.request != 0 || schedule_due;
        if refactor {
            self.work.decided(self.work.w);
        }
        self.last_pivoted = pv.is_some();
        if let Some(pv) = pv {
            self.last_cap = pv.col.len();
            let alpha = self.pending_alpha.take().unwrap();
            if refactor {
                // The weights of the pivot are updated before the factorization renumbers the positions; the basic solution is recomputed.
                self.update_weights(pv);
                self.pivot_cols = Some(pv.col.len());
                self.refactor(if schedule_due { "scheduled" } else { "request" })?;
            } else {
                self.factor
                    .update(pv.r, alpha.iter().map(|v| v.to_f64()).collect());
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
            // The step divides by the pivot as the double `ctx+0x308` (the row's pivot, stored as a double); the same for the primal step and `ctx+0x300`.
            let theta = self.d[pv.q] / R::c(aq.to_f64());
            // A zero step skips the loop (`0xe4719`): no work.
            if !refactor && theta != R::ZERO {
                let row_cnt = pv.row.len() + usize::from(pv.leave < self.n);
                self.work
                    .add(0xe7000, 3.0 * row_cnt as f64, crate::work::CTX);
            }
            let old_leave = self.d[pv.leave];
            // The reference's row holds a slack leaving variable with its computed coefficient (about 1, 0.9999999999999997 in PILOT.JA call 30657), the step on it uses that.
            let lcoef = pv.row.iter().find(|e| e.0 == pv.leave).map_or(R::ONE, |e| {
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
                theta + R::ZERO
            } else if pv.above {
                -(old_leave + theta * lcoef)
            } else {
                old_leave - theta * lcoef
            };
        } else if refactor {
            self.refactor(if schedule_due { "scheduled" } else { "request" })?;
        }
        if refactor {
            self.pricing_switch();
            if self.request == 0 {
                self.request = 10;
            }
            self.recompute_d();
            self.shift_wrong_d();
            let level_before = self.level;
            self.set_level_work(level_before != 0);
            self.set_level(level_before);
            let req = self.request;
            self.request = 0;
            self.recompute_xb(req);
            if req == 4
                && (!R::QUAD && self.level != 0)
                && self.status.is_none()
                && !NO_ART_EXIT.with(|c| c.get())
            {
                // The confirming refactorization leaves variables at an artificial bound: the post-refactor pass (`0xe45c0`, `0xe5e33..0xe5e63`) counts
                // them and, with `ctx+0x35c` != 0 (1 from the start), ends the dual at once with status -1004 (s16, s198 `dual_exit`: the exit comes
                // with `ctx+0x294` = 4 pending, before any further pricing call). The controller continues with the primal from this basis.
                self.art_exit = true;
                self.dual_unbounded();
                return Ok(());
            }
            self.weight_overflow_gate(req);
            self.prog_on = self.level == 0;
            self.prog_acc = 0.0;
            self.recompute_flags(
                if R::QUAD {
                    level_before != 0
                } else {
                    self.big_pricing
                },
                pv.is_some(),
            );
            self.snapshot_recomputed();
            self.flags_work();
            if R::QUAD
                && self.status.is_none()
                && self.level == 0
                && self.lu_thr >= 0.5
                && (req == 4 || crate::simplex::QUAD_PRIMALS.with(|c| c.get()) >= 10)
                && self.tiny_infeasibility()
            {
                // The recompute after a factorization of a pass that follows the tenth quad primal (`ctx+0x35c` above 10, `0x21a315..0x21a4b6`) takes the basic solution as
                // feasible when the infeasibilities are small against the size of x_B: the flags are all 2 and the dual ends optimal (seed 2950, 3e11 against 7e-3).
                self.flags.iter_mut().for_each(|f| *f = 2);
                self.list.clear();
                self.status = Some(Status::Optimal);
            }
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
    fn update_weights(&mut self, pv: &Pivot<R>) {
        let (r, alpha, w) = (pv.r, pv.alpha, pv.w);
        // The quad weight update works on the double-rounded rho norm and pivot (verified against the reference W bits, `0x225820` callers).
        let (alpha, w) = if R::QUAD {
            (R::c(alpha.to_f64()), R::c(w.to_f64()))
        } else {
            (alpha, w)
        };
        if crate::dump::on("dtau") {
            crate::dump::emit(
                "dtau",
                serde_json::json!({"r": r, "alpha": alpha.to_f64(), "w": w.to_f64(), "col": crate::dump::sparse_r(&pv.col), "tau": pv.tau.iter().map(|v| v.to_f64()).collect::<Vec<_>>()}),
            );
        }
        if w == R::ZERO {
            self.w[r] = R::c(W_EPS);
        } else {
            let a2 = if alpha >= R::ZERO {
                maxr(alpha, R::c(W_EPS))
            } else {
                minr(alpha, R::c(-W_EPS))
            };
            self.w[r] = (w / a2) / a2;
            for (j, &(i, cv)) in pv.col.iter().enumerate() {
                let a = cv / a2;
                let t = pv.tau[j];
                // Quad kernel `0x2165c0..0x21668d`: a * ((w a) - 2 tau) + W (the double one is (w a) a - (2 a) tau + W).
                let x = if R::QUAD {
                    (a * ((w * a) - t * R::c(2.0))) + self.w[i]
                } else {
                    ((w * a) * a - (R::c(2.0) * a) * t) + self.w[i]
                };
                self.w[i] = maxr(R::c(W_EPS), x);
            }
        }
    }

    fn update_primal_and_weights(&mut self, pv: &Pivot<R>) {
        self.update_weights(pv);
        let (r, alpha) = (pv.r, pv.alpha);
        let mut changed: Vec<usize> = Vec::new();
        let v = pv.leave;
        let b = if self.st[v] == LB || self.st[v] == FIXED {
            self.p.lo[v]
        } else {
            self.p.up[v]
        };
        let theta = (self.xb[r] - R::c(b)) / R::c(alpha.to_f64());
        self.xb[r] = pv.t0 + theta;
        for &(i, cv) in &pv.col {
            self.xb[i] = self.xb[i] - cv * theta;
            changed.push(i);
        }
        changed.push(r);
        // Running average of the changed-row count and the growth of the list threshold (`0xc0532`, `0xc06d7`): while more rows are on the list than
        // that average and than m^0.45 the threshold grows by 10 percent per iteration (before the rows are touched).
        let ninf = self.list.len() as f64;
        let it = (self.iters + self.it_base) as f64;
        self.avg_touch = (2.0 * pv.col.len() as f64 + (it - 2.0) * self.avg_touch) / it;
        if self.value_mode {
            if ninf > self.avg_touch.max((self.m as f64).powf(0.45).floor()) {
                self.thr *= R::c(1.1);
            }
            // The pivot row leaves the list at the next pricing call whatever its value (`0xc0732`); its priority is refreshed for the entering variable.
            // An entering variable that lands beyond its bound stays on the list with its infeasibility (BOEING2 call 57); a feasible one gets the distance
            // from its lower bound and flag 3.
            let q = self.head[r];
            let fi = self.infeas(r);
            if fi > self.ptol() {
                self.prio[r] = self.pri(r, self.infeas_real(r));
            } else {
                let d = self.xb[r] - R::c(self.p.lo[q]);
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
fn rho_norm_dense<R: Real>(rho: &[(usize, R)], m: usize) -> R {
    let mut sq = vec![R::ZERO; m];
    for &(i, v) in rho {
        sq[i] = v * v;
    }
    sum_squares(&sq)
}

fn rho_norm<R: Real>(rho: &[(usize, R)]) -> R {
    if R::QUAD {
        // The quad kernel (`0x25d3ff..0x25d49a`) squares in extended precision and keeps the running sum in a double: each step is
        // rounded to extended by the add and then to double by the store.
        return rho
            .iter()
            .fold(R::ZERO, |a, e| R::c((a + e.1 * e.1).to_f64()));
    }
    let sq: Vec<R> = rho.iter().map(|e| e.1 * e.1).collect();
    sum_squares(&sq)
}

fn sum_squares<R: Real>(sq: &[R]) -> R {
    let n = sq.len();
    if !R::LANES || n <= 6 {
        return sq.iter().fold(R::ZERO, |a, &x| a + x);
    }
    let end = n - (n & 7);
    let mut acc = [[R::ZERO; 2]; 4];
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
    for &x in &sq[end..] {
        t += x;
    }
    t
}

/// State the reference's driver keeps across a singular factorization (`0xf398e..0xf3b84`): the retry counter `ctx+0x184`, the LU pivot threshold `ctx+0x390`,
/// the factorization count `ctx+0x2c0` (the failed one included) and the iterations done.
#[derive(Clone, Copy)]
struct Carry {
    retries: u32,
    thr: f64,
    fac: usize,
    iters: usize,
    avg_touch: f64,
    last_cap: usize,
    inf_count: usize,
    /// The stale leaving row `ctx+0x2f8` was set when the pass failed (the failing factorization followed a pivot, not a no-row request): the first recompute of the next
    /// pass compares the infeasible-row count with the carried pivot-column length (`0xc1e92`: `edx >= 0`).
    after_pivot: bool,
    /// The iteration counter `ctx+0x20` at the start of the pass: 0 after a restart, the iterations so far when the dual runs again after the primal.
    it_base: usize,
    /// The factorization count `ctx+0x2c0` at the start of the pass: 0 after a restart (a fresh context), the count so far plus the primal's factorizations when the dual runs again
    /// after the primal.
    nref0: usize,
    /// The perturbation scale `ctx+0x348` the pass starts with: the default, or after the primal's clean-up max(1e-7, 0.05 times the previous pass's scale) (reference quad runs, seeds 2950, 2840).
    pstep: Option<f64>,
}

thread_local! {
    /// Set where a mid-run factorization found no pivot for a basic column: `(iterations done, factorizations including the failed one, LU threshold)`.
    /// The pass ended with the quad switch (`Step::SwitchQuad`): the counters and the basis the quad pass continues from.
    static SWITCH: std::cell::Cell<Option<(usize, usize, f64, f64, usize, usize, bool, u32, usize, usize)>> = const { std::cell::Cell::new(None) };
    static FROM_SWITCH: std::cell::Cell<bool> = const { std::cell::Cell::new(false) };
    static SWITCH_BASIS: std::cell::RefCell<Option<Basis>> = const { std::cell::RefCell::new(None) };
    /// Set by the primal's clean-up pass: the dual runs inside the primal's context, with the LU threshold the primal's accuracy ladder left in `ctx+0x390` (gen3 seeds 247, 680: raised to 0.5 at the first x_B recompute; with 2^-7 the pivot order of the final factorization differs).
    pub static CLEANUP_PASS: std::cell::Cell<Option<f64>> = const { std::cell::Cell::new(None) };
    /// The perturbation flag the quad pass inherits from the double pass (`ctx+0x338`).
    /// Out of a clean-up pass dual: (the dual ended with cost shifts, so its last factorization was the sorted primal clean-up one, the LU threshold at the end).
    pub static CLEANUP_OUT: std::cell::Cell<(bool, f64)> = const { std::cell::Cell::new((false, 0.0)) };
    static SWITCH_PERT: std::cell::Cell<Option<u8>> = const { std::cell::Cell::new(None) };
    static SINGULAR: std::cell::Cell<Option<(usize, usize, f64, f64, usize, usize, bool, u32, usize)>> = const { std::cell::Cell::new(None) };
    static SINGULAR_BASIS: std::cell::RefCell<Option<Basis>> = const { std::cell::RefCell::new(None) };
}

/// Dual simplex from a cold start (the crash basis) or a supplied basis. A singular factorization in the middle of the run restarts the dual from
/// the slack basis (the reference's driver, `0xf398e`): the retry counter picks the pivot threshold of the next factorizations, 0.5 for the first
/// retry and 0.9 for the later ones, and the iteration count goes on.
pub fn solve(
    model: &LpModel,
    start: Option<&Basis>,
    opts: &SolveOptions,
) -> Result<SolveResult, SolveError> {
    // The reference works on the mirror image of a column with no lower bound and a finite upper bound (`x' = -x`, bounds `[-up, inf)`, `0xeb61b`
    // sees `lb = -5` for `x in (-inf, 5]`): every rule below (crash tests on the coefficient sign, the status by the cost sign, the d signs) runs on that frame.
    let flip: Vec<bool> = (0..model.num_cols())
        .map(|j| model.col_lo[j] <= -crate::model::INF && model.col_up[j] < crate::model::INF)
        .collect();
    if !flip.iter().any(|&f| f) {
        return solve_frame(model, start, opts);
    }
    let mut mirrored = model.clone();
    for j in 0..model.num_cols() {
        if flip[j] {
            for k in model.col_start[j]..model.col_start[j + 1] {
                mirrored.val[k] = -model.val[k];
            }
            mirrored.obj[j] = -model.obj[j];
            mirrored.col_lo[j] = -model.col_up[j];
            mirrored.col_up[j] = crate::model::INF;
        }
    }
    let start = start.map(|b| {
        let mut b = b.clone();
        for j in (0..model.num_cols()).filter(|&j| flip[j]) {
            // The mirror of a variable at its upper bound (the only finite one) is at its lower bound.
            b.at_upper[j] = !b.at_upper[j];
        }
        b
    });
    let mut r = solve_frame(&mirrored, start.as_ref(), opts)?;
    for j in (0..model.num_cols()).filter(|&j| flip[j]) {
        r.x[j] = -r.x[j];
        r.reduced_cost[j] = -r.reduced_cost[j];
        r.basis.at_upper[j] = !r.basis.at_upper[j];
        if let Some(ray) = r.ray.as_mut() {
            ray[j] = -ray[j];
        }
    }
    Ok(r)
}

thread_local! {
    /// Quad resumes after the primal gave up at an artificial-bound exit (one per solve).
    static QUAD_RESUMES: std::cell::Cell<u32> = const { std::cell::Cell::new(0) };
    static NO_ART_EXIT: std::cell::Cell<bool> = const { std::cell::Cell::new(false) };
}

fn solve_frame(
    model: &LpModel,
    start: Option<&Basis>,
    opts: &SolveOptions,
) -> Result<SolveResult, SolveError> {
    let mut carry = Carry {
        retries: 0,
        thr: CLEANUP_PASS
            .with(|c| c.get())
            .unwrap_or(crate::lu_order::THRESHOLD),
        fac: 0,
        iters: 0,
        avg_touch: 0.0,
        last_cap: 0,
        inf_count: 0,
        after_pivot: false,
        it_base: 0,
        nref0: 0,
        pstep: None,
    };
    set_dual_tol(DUAL_TOL0);
    set_dual_infeas_sum(0.0);
    set_primal_tol(PRIMAL_TOL);
    RESID_FAILS.with(|c| c.set(0));
    PERT_RERUNS.with(|c| c.set(0));
    QUAD_RESUMES.with(|c| c.set(0));
    RNG_CELL.with(|c| c.set(101));
    SNAP.with(|c| *c.borrow_mut() = None);
    INITIAL_SNAP.with(|c| *c.borrow_mut() = None);
    SNAP_QUALITY.with(|c| c.set(0));
    DRIFT_START.with(|c| *c.borrow_mut() = None);
    SNAP_C40.with(|c| c.set(-100.0));
    PAIR_TOL.with(|c| c.set(3e-6));
    PASS_NO.with(|c| c.set(0));
    DUAL_HANDOFFS.with(|c| c.set(0));
    C35C.with(|c| c.set(1));
    C178.with(|c| c.set(0));
    crate::simplex::HANDOFF_ROWS.with(|c| *c.borrow_mut() = None);
    LARGE_XB.with(|c| c.set(false));
    crate::simplex::QUAD_PRIMALS.with(|c| c.set(0));
    // The quad mode (`ctx+0x28`): once the controller switched, every later pass runs the extended-precision twin of the dual.
    let mut quad = false;
    let mut switch_start: Option<Basis> = None;
    let mut restart_start: Option<Basis> = None;
    // The factor's `+0x0c`: zero supplies the current statuses to setup,
    // one asks setup for a fresh crash basis.
    let mut factor_crash = start.is_none();
    // The controller snapshots the tolerance once (`0x154ccc`), before any
    // singular retries; its primal handoffs scale this original value.
    LOOP_TOL.with(|c| c.set(dual_tol()));
    loop {
        SINGULAR.with(|c| c.set(None));
        SINGULAR_BASIS.with(|c| *c.borrow_mut() = None);
        SWITCH.with(|c| c.set(None));
        let sw = switch_start.take();
        let restart_basis = restart_start.take();
        let drift_start = DRIFT_START.with(|c| c.borrow_mut().take());
        let start_now = sw
            .as_ref()
            .or(drift_start.as_ref())
            .or(restart_basis.as_ref())
            .or_else(|| if carry.retries == 0 { start } else { None });
        FROM_SWITCH.with(|c| c.set(sw.is_some()));
        let res = if quad {
            solve_pass::<crate::f80::F80>(model, start_now, opts, &carry)
        } else {
            solve_pass::<f64>(model, start_now, opts, &carry)
        };
        match res {
            Err(SolveError::SingularBasis) if SWITCH.with(|c| c.get()).is_some() => {
                let (
                    iters,
                    fac,
                    thr,
                    avg_touch,
                    last_cap,
                    inf_count,
                    after_pivot,
                    pass_retries,
                    it_base,
                    nref,
                ) = SWITCH.with(|c| c.take()).unwrap();
                let b = SWITCH_BASIS.with(|c| c.borrow_mut().take()).unwrap();
                if crate::dump::on("dcr") {
                    crate::dump::emit(
                        "dcr",
                        serde_json::json!({"switch": true, "iters": carry.iters + iters}),
                    );
                }
                quad = true;
                factor_crash = false;
                carry = Carry {
                    retries: pass_retries,
                    thr,
                    fac,
                    iters: carry.iters + iters,
                    avg_touch,
                    last_cap,
                    inf_count,
                    after_pivot,
                    // The precision switch keeps the iteration count in the
                    // same context; a singular restart below starts it anew.
                    it_base,
                    // A direct switch keeps the factor count as well.
                    nref0: nref,
                    pstep: None,
                };
                switch_start = Some(b);
            }
            Err(SolveError::SingularBasis) if carry.retries < 100 => {
                crate::simplex::HANDOFF_ROWS.with(|c| *c.borrow_mut() = None);
                let Some((
                    iters,
                    fac,
                    thr,
                    avg_touch,
                    last_cap,
                    inf_count,
                    after_pivot,
                    pass_retries,
                    context_iters,
                )) = SINGULAR.with(|c| c.take())
                else {
                    return Err(SolveError::SingularBasis);
                };
                // The restart handler (`0xcedb0`) switches to quad when iterations were made (`ctx+0x20 + ctx+0x2a0`), with `Quad` at its default (-1, the release
                // profile) and `ctx+0x28` still 0 (the checks of `0x2c8d0` and `ctx+0x188` hold in a pure dual run).
                if crate::dump::on("dcr") {
                    crate::dump::emit(
                        "dcr",
                        serde_json::json!({"restart": true, "quad": quad, "perturb": opts.perturb, "iters": carry.iters + iters}),
                    );
                }
                if !quad && opts.perturb > 0.0 && carry.iters + iters > 0 {
                    quad = true;
                }
                // The third restart starts the pricing from scratch (list-walk mode, empty running averages): the reference's p40/thr schedule over the
                // restarts is value (thr maxp/1e3), value (thr maxp/1e5), list-walk, then value again (fx_midsing_m5, seeds 1494, 1390, 1496, 50010).
                let fresh = pass_retries + 1 == 3;
                if fresh {
                    factor_crash = true;
                } else if context_iters > 0 && SNAP_QUALITY.with(|c| c.get() > 0) {
                    factor_crash = false;
                    let snap = SNAP
                        .with(|c| c.borrow().clone())
                        .or_else(|| INITIAL_SNAP.with(|c| c.borrow().clone()));
                    restart_start = snap.map(|(basic, at_upper)| Basis { basic, at_upper });
                } else if !factor_crash {
                    restart_start = SINGULAR_BASIS.with(|c| c.borrow_mut().take());
                }
                let (avg_touch, last_cap, inf_count) = if fresh {
                    (0.0, 0, 0)
                } else {
                    (avg_touch, last_cap, inf_count)
                };
                carry = Carry {
                    retries: pass_retries + 1,
                    thr: thr.max(if pass_retries == 0 { 0.5 } else { 0.9 }),
                    fac,
                    iters: carry.iters + iters,
                    avg_touch,
                    last_cap,
                    inf_count,
                    after_pivot,
                    it_base: 0,
                    nref0: 0,
                    pstep: None,
                };
            }
            Ok(mut r) => {
                r.iterations += carry.iters;
                return Ok(r);
            }
            e => return e,
        }
    }
}

/// A restart pass that starts the pricing from scratch: only the third (`ctx+0x184 == 3`, `0xcedb0`), in the double and in the quad pass (the passes at 4, 5
/// and 6 restarts price in the value mode: random LP seed 50010, restarts at calls 9 .. 54).
fn fresh_pass<R: Real>(retries: u32) -> bool {
    retries > 0 && retries == 3
}

fn solve_pass<R: Real>(
    model: &LpModel,
    start: Option<&Basis>,
    opts: &SolveOptions,
    carry: &Carry,
) -> Result<SolveResult, SolveError> {
    if fresh_pass::<R>(carry.retries) {
        // The fresh start of the third restart initialises the context again (`0xcf29b`): the tolerance goes back to the parameter, the counter stays.
        set_dual_tol(DUAL_TOL0);
        set_primal_tol(PRIMAL_TOL);
        LARGE_XB.with(|c| c.set(false));
    }
    let p = build(model);
    let (n, m) = (p.n, p.m);
    let mut dropped = None;
    let mut cleaned = false;
    let mut supplied = None;
    let mut crash_upper: Option<Vec<bool>> = None;
    let basis = match start {
        Some(b) => {
            let mut b = b.clone();
            b.validate(n, m).map_err(SolveError::InvalidBasis)?;
            b.basic.sort_unstable();
            supplied = Some(b.basic.clone());
            // A singular basis loses the columns that get no pivot in the reference's LU order and takes the slacks of the uncovered rows instead.
            let mut probe = b.clone();
            if FROM_SWITCH.with(|c| c.get())
                && crate::simplex::repair_basis(&p, &mut probe).is_some()
            {
                // The basis the quad switch carries over is not a supplied one: when its quad factorization is singular (the double structure finds the same columns
                // dependent), the pass ends before its first iteration and the restart takes the snapshot basis (seed 2706: first quad factorization has `p13` = 1,
                // the next one is the slack basis, `ctx+0x184` = 1).
                SINGULAR_BASIS.with(|c| *c.borrow_mut() = Some(b.clone()));
                SINGULAR.with(|c| {
                    c.set(Some((
                        0,
                        carry.fac + 1,
                        carry.thr,
                        carry.avg_touch,
                        carry.last_cap,
                        carry.inf_count,
                        carry.after_pivot,
                        carry.retries,
                        carry.it_base,
                    )))
                });
                return Err(SolveError::SingularBasis);
            }
            dropped = crate::simplex::repair_basis(&p, &mut b).map(|(j, _)| j);
            b
        }
        None => {
            let (basic, upper) = dual_crash_basic(&p);
            crash_upper = Some(upper);
            Basis {
                basic,
                at_upper: vec![false; n + m],
            }
        }
    };
    basis.validate(n, m).map_err(SolveError::InvalidBasis)?;
    if carry.retries == 0 {
        crate::lu_order::reset_pair_gate();
    }
    crate::lu_order::set_threshold(carry.thr);
    crate::lu_order::set_pool_dims(p.n, p.cols[..p.n].iter().map(Vec::len).sum());
    let (mut order, mut stc, mut rows) = reference_order::<R>(&p, &basis.basic, carry.fac + 1);
    dump_start(&basis.basic, &order, &stc, carry.thr);
    let mut mirror = RefStruct::fresh(m, n, &stc, &rows);
    let mut carry_thr = carry.thr;
    // Growth check of the setup factorization (quad setup `0x246d39..0x246df1`, double driver `0xf351d..0xf3629`): max|diag| / min|diag| above 1e15 (1e12 scaled)
    // raises the threshold to 0.5, then 0.9, and the same columns are factorized again.
    if order.len() == m {
        while carry_thr < 0.9 {
            let (mx, mn) = mirror
                .numerics()
                .0
                .iter()
                .fold((0.0f64, f64::INFINITY), |(a, b), v| {
                    (a.max(v.abs().to_f64()), b.min(v.abs().to_f64()))
                });
            if !(mx / mn > if opts.scaled { 1e12 } else { 1e15 }) {
                break;
            }
            carry_thr = if carry_thr < 0.5 { 0.5 } else { 0.9 };
            crate::lu_order::set_threshold(carry_thr);
            (order, stc, rows) = reference_order::<R>(&p, &basis.basic, carry.fac + 1);
            dump_start(&basis.basic, &order, &stc, carry_thr);
            if order.len() != m {
                break;
            }
            mirror = RefStruct::fresh(m, n, &stc, &rows);
        }
    }
    let cols: Vec<&Column> = order.iter().map(|&j| &p.cols[j]).collect();
    let factor = if R::QUAD {
        // The quad pass starts from an ill-conditioned basis the double factor may find singular; the quad LU decides (the double factor is only for read-out).
        Factor::with_floor(m, &cols, 0.0)
    } else {
        Factor::new(m, &cols)
    }
    .map_err(|_| SolveError::SingularBasis)?;
    let mut st = vec![LB; n + m];
    for (pos, &j) in order.iter().enumerate() {
        st[j] = pos as i32;
    }
    let cap = ((m as f64).powf(0.7).floor() as i64).max(20);
    let mut s: Dual<R> = Dual {
        p: &p,
        n,
        m,
        head: order,
        st,
        xb: vec![R::ZERO; m],
        d: vec![R::ZERO; n + m],
        w: vec![R::ONE; m],
        mirror,
        factor,
        rrows: vec![Vec::new(); m],
        absorb_request: false,
        refactorizations: carry.fac + 1,
        lu_thr: carry_thr,
        big_pricing: false,
        pmode: if start.is_some() || (carry.retries > 0 && !fresh_pass::<R>(carry.retries)) {
            3
        } else {
            1
        },
        n_infeas: 0,
        inf_count: carry.inf_count,
        retries: carry.retries,
        repaired: None,
        snapshot_excluded: dropped.is_some(),
        bad_rows: Vec::new(),
        old_rej: 0,
        rejected_end: false,
        art_exit: false,
        worst: 0.0,
        worst_row: -1,
        retry: false,
        retry_end: 0,
        retry_tiny: 0,
        stale_count: false,
        last_pivoted: false,
        it_base: carry.it_base,
        dflag: Vec::new(),
        switch_enabled: opts.perturb > 0.0,
        level: 0,
        request: 0,
        iters: 0,
        prio: vec![R::ZERO; m],
        ptype: 0,
        p68: Vec::new(),
        p70: Vec::new(),
        p60: 0.0,
        p58: 0.0,
        flags: vec![2; m],
        list: Vec::new(),
        thr: R::ZERO,
        xb_residual: R::ZERO,
        value_mode: false,
        cursor: 0,
        walk_len: 0,
        slack_nb: vec![false; m],
        pend: cap,
        cap,
        rebuilds: 40,
        first_recompute: carry.retries == 0 || fresh_pass::<R>(carry.retries),
        warm: start.is_some(),
        avg_touch: carry.avg_touch,
        last_cap: carry.last_cap,
        last_y: Vec::new(),
        fac0: carry.fac,
        nref0: carry.nref0,
        nref_skip: 0,
        it0: carry.iters,
        updates: 0,
        upd_limit: None,
        pending_alpha: None,
        pair: (0.0, 0.0),
        status: None,
        nm: model.num_cols(),
        ray: None,
        zero_artificial: false,
        cost: p.cost.clone(),
        cshift: vec![0.0; p.cost.len()],
        cost_lo: vec![0.0; p.cost.len()],
        shift_r: vec![R::ZERO; p.cost.len()],
        rng: RNG_CELL.with(|c| c.get()),
        pert: 0,
        obj_off: 0.0,
        flip_select: start.is_some(),
        work: if R::QUAD {
            crate::work::Work::dual_quad()
        } else {
            crate::work::Work::dual()
        },
        pivot_cols: None,
        nnz_a: p.cols[..n].iter().map(Vec::len).sum(),
        pstep: carry.pstep.unwrap_or(opts.perturb),
        stall: 0,
        stall_snap: 0,
        stall_it: 0.0,
        prog_on: false,
        prog_acc: 0.0,
        best_obj: -1e100,
        pert_left: 1_000_000_000,
        pert_param: opts.perturb,
        obj_scale: opts.obj_scale,
        scaled: opts.scaled,
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
    let mut setup_free = None;
    if free_phase {
        setup_free = s.free_pivots()?;
        if setup_free.is_none() {
            s.refactor("free_setup")?;
        }
    }
    s.build_rrows();
    if let Some(rows) = crate::simplex::HANDOFF_ROWS.with(|c| c.borrow_mut().take()) {
        s.rrows = rows;
    }
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
        if s.st[j] == FREE && dt[j].abs() > R::c(dual_tol()) {
            // A free variable the setup pivots could not bring into the basis (no basic variable with a finite lower bound in its column) and that is dual
            // infeasible: it takes the artificial bound its d asks for, like any other variable without the bound it would need.
            s.st[j] = if dt[j] > R::ZERO { LB } else { UB };
            s.d[j] = dt[j].abs();
            continue;
        }
        if s.st[j] == FIXED || s.st[j] == FREE {
            s.d[j] = if s.st[j] == FIXED && j < n {
                R::c(s.cost[j])
            } else {
                dt[j]
            };
            continue;
        }
        if dt[j] < R::c(-dual_tol()) {
            s.st[j] = UB;
        } else if dt[j] > R::c(dual_tol()) {
            s.st[j] = LB;
        } else if crash_upper.as_ref().is_some_and(|u| u[j]) {
            // A d within the tolerance keeps the status the crash gave (`0xead36`: at the upper bound for a negative cost, else at the lower bound, then
            // the sign of the cost after the crash's third pass; f94 release: x120 with cost -1 and d = 0 starts at the artificial upper bound).
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
    if supplied.is_none() && carry.retries > 0 && !fresh_pass::<R>(carry.retries) {
        // A pass after a singular factorization runs in the pricing mode 3 from the slack basis (fx_midsing_1106: the second vector of the two-vector FTRAN is
        // dense from the first iteration of the second pass, the weights differ from the first pass in the last bits).
        s.dflag = s.st.iter().map(|&t| t < 0).collect();
    }
    // Initial weights: the squared norms of the rows of B^-1 (a supplied basis starts with unit weights).
    for r in 0..if s.warm { 0 } else { m } {
        let out = s.mirror.btran(&[(r, R::ONE)], false, false);
        s.work
            .add_from(0x12de50, 0xc36f1, s.mirror.work, crate::work::FAC);
        s.w[r] = if s.mirror.bt_dense {
            rho_norm_dense(&out, m)
        } else {
            rho_norm(&out)
        };
    }
    s.stall_check(0);
    // The setup's post-refactor pass (`0xe6dcc`) sees the iterations of the free-column pivots.
    if !s.flip_select && 3.0 * s.iters as f64 > s.m as f64 {
        s.flip_select = true;
    }
    if R::QUAD {
        // The quad pass keeps the perturbation flag of the double pass but not its cost shifts: d is recomputed from the true costs (seed 2240: the shifted
        // d of the double pass are 10^3..10^4 times the reference's at the first quad call, and so are the d of every quad refactorization).
        if let Some(pert) = SWITCH_PERT.with(|c| c.take()) {
            s.pert = pert;
        }
    }
    s.recompute_xb(0);
    s.prog_on = s.level == 0;
    s.prog_acc = 0.0;
    s.work.add(0xc391f, m as f64, crate::work::CTX);
    s.recompute_flags(
        if R::QUAD {
            s.warm && s.level != 0
        } else {
            s.big_pricing
        },
        carry.retries > 0 && carry.after_pivot,
    );
    s.snapshot_recomputed();
    // The carried average decides the pricing mode of the first recompute; the running average of the pass (`P+0x48`) starts at 0 (m6, m5 captures).
    s.avg_touch = 0.0;
    s.flags_work();
    s.log_snapshot();
    // A schedule with `usize::MAX` separators is one list of per-pass iteration counts for each pass of a run with restarts (the counter restarts at 0).
    let pass_no = PASS_NO.with(|c| {
        let k = c.get();
        c.set(k + 1);
        k
    });
    let schedule = opts.refactor_schedule.as_ref().map(|sch| {
        if sch.contains(&usize::MAX) {
            sch.split(|&x| x == usize::MAX)
                .nth(pass_no)
                .map_or_else(Vec::new, |v| v.to_vec())
        } else {
            sch.clone()
        }
    });
    let mut trace = Vec::new();
    let mut k = 0usize;
    let mut farkas = None;
    let mut status;
    // A deferred free column that has a reduced cost and no row to enter on ends the dual before its first iteration (`0xea0b7`, status 5).
    if let Some((j, d)) = setup_free {
        // The deferred-free exit records the entering column in ctx+0x2fc. Its ray is the unit step against d and the corresponding basic response,
        // not the aggregate artificial-bound point used by the no-row handler.
        let dir = -d.signum();
        let mut col = vec![0.0; m];
        for &(i, a) in &p.fcols[j] {
            col[i] = a * dir;
        }
        s.factor.ftran(&mut col);
        let mut ray = vec![0.0; n + m];
        ray[j] = dir;
        for (i, &v) in col.iter().enumerate() {
            ray[s.head[i]] = -v;
        }
        ray.truncate(s.nm);
        s.ray = Some(ray);
        s.status = Some(Status::Unbounded);
    }
    loop {
        k += 1;
        if setup_free.is_some() {
            status = Status::Unbounded;
            break;
        }
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
                // The driver around the pricing (`0xf6745..0xf6768`, `0xf8b51..0xf8b77`) clears the worst rejected row (`ctx+0x430` = -1, `ctx+0x438` = 0)
                // each time it returns, so the "first rejection" is per iteration; the rejected count `ctx+0x420` stays (fx_midsing_1106 call 20).
                s.driver_return();
            }
            Step::Reject => {
                // A row dropped in the list-walk mode returns from the iteration function with the request it set (`0xc9577`): the post-iteration pass refactorizes.
                if s.request != 0 && s.status.is_none() {
                    s.post(None, false)?;
                    s.driver_return();
                }
            }
            Step::Mismatch => {
                if s.status.is_none() {
                    s.post(None, false)?;
                }
                s.driver_return();
            }
            Step::Restart => {
                SINGULAR_BASIS.with(|c| {
                    *c.borrow_mut() = Some(Basis {
                        basic: s.head.clone(),
                        at_upper: s.st.iter().map(|&t| t == UB).collect(),
                    })
                });
                SINGULAR.with(|c| {
                    c.set(Some((
                        s.iters,
                        s.refactorizations,
                        s.lu_thr,
                        if R::QUAD && s.iters == 0 {
                            carry.avg_touch
                        } else {
                            s.avg_touch
                        },
                        s.last_cap,
                        s.inf_count,
                        // The stale leaving row is cleared on this path (`ctx+0x2f8` = -1 in the recompute of the next pass, fx_midsing_1390).
                        false,
                        s.retries,
                        s.it_base + s.iters,
                    )))
                });
                return Err(SolveError::SingularBasis);
            }
            Step::SwitchQuad => {
                SWITCH.with(|c| {
                    c.set(Some((
                        s.iters,
                        s.refactorizations,
                        s.lu_thr,
                        s.avg_touch,
                        s.last_cap,
                        s.inf_count,
                        false,
                        s.retries,
                        s.it_base + s.iters,
                        s.nref(),
                    )))
                });
                SWITCH_PERT.with(|c| c.set(Some(s.pert)));
                SWITCH_BASIS.with(|c| {
                    *c.borrow_mut() = Some(Basis {
                        basic: s.head.clone(),
                        at_upper: (0..n + m).map(|j| s.st[j] == UB && p.up[j] < BIG).collect(),
                    })
                });
                return Err(SolveError::SingularBasis);
            }
            Step::NoRow => {
                if s.request != 3 {
                    s.no_row();
                }
                if s.status.is_none() {
                    s.post(None, false)?;
                }
                s.driver_return();
            }
            Step::Unbounded(y) => {
                if s.updates != 0 {
                    // No entering variable on a factor with updates (ctx+0x294 = 5, RECIPE-like bgdbg1 calls 29 and 31): refactorize and price again; the
                    // certificate is the row of the first such failure on a fresh factor.
                    s.request = 5;
                    s.post(None, false)?;
                    if s.status.is_none() {
                        continue;
                    }
                } else {
                    status = Status::Infeasible;
                    farkas = Some(y);
                    break;
                }
            }
        }
        if let Some(stt) = s.status {
            status = stt;
            if R::QUAD
                && matches!(stt, Status::Optimal | Status::Unbounded)
                && s.pert != 0
                && !NO_ART_EXIT.with(|c| c.get())
            {
                // The quad driver's end (`0x24a4f6..0x24a548`): a final status 2 (dual optimal) or 5 with the cost shifts in force (`ctx+0x338` != 0) clears them and returns
                // status -1004; the controller hands over to the x87 primal (seeds 5247, 6972, 6610, 8179). It replaces the double dual's exit at the request-4 refactor.
                s.art_exit = true;
                s.dual_unbounded();
                status = Status::Unbounded;
                break;
            }
            if !R::QUAD
                && stt == Status::Optimal
                && s.pert != 0
                && !NO_ART_EXIT.with(|c| c.get())
                && CLEANUP_PASS.with(|c| c.get()).is_none()
            {
                // The double dual ends optimal with cost shifts in force: restore the costs and start the primal (-1004).

                // The primal can return optimal, restart the dual, or report a singular starting factor.
                s.rejected_end = true;
                status = Status::Infeasible;
                break;
            }
            if stt == Status::Optimal && s.pert != 0 {
                // The end of the dual with the cost shifts in force (`ctx+0x338` != 0, `0xf7033`): the shifts are removed (status -1004) and the primal
                // clean-up starts with a factorization (`0x142350` from `0xe3ee1`), the reference's last one.
                s.cost.clone_from(&p.cost);
                s.obj_off = 0.0;
                s.pert = 0;
                // The clean-up factorization is handed the basic columns in ascending variable order (PILOT.JA `lu_in`: the last record is sorted, the dual's own are in
                // basis order); the final x_B comes from that factor.
                let mut wv = vec![R::ZERO; s.n + s.m];
                for (pos, &v) in s.head.iter().enumerate() {
                    wv[v] = s.w[pos];
                }
                s.head.sort_unstable();
                for pos in 0..s.m {
                    let v = s.head[pos];
                    s.st[v] = pos as i32;
                    s.w[pos] = wv[v];
                }
                if let Err(e) = s.refactor("optimal_cleanup") {
                    // The double driver's shifted optimum returns -1004 and
                    // increments ctx+0x35c before the primal cleanup. If our
                    // substitute cleanup factor is singular, preserve that
                    // increment across its restart (seed 2594).
                    if !R::QUAD {
                        C35C.with(|c| c.set(c.get() + 1));
                    }
                    return Err(e);
                }
                s.recompute_xb(0);
                cleaned = true;
                CLEANUP_OUT.with(|c| c.set((true, s.lu_thr)));
            }
            break;
        }
    }
    if !cleaned {
        CLEANUP_OUT.with(|c| c.set((false, s.lu_thr)));
    }
    let restored = dropped.or(s.repaired);
    // The eleventh restore (`0x15553e`, `PERT_RERUNS` past 10) ends the solve with status 12 before the primal runs: the reference's iteration count has no primal pivot
    // in that last cycle, and its x is that of the dual's last basis (seed 2147: 189 iterations).
    if status == Status::Infeasible && s.pert != 0 && PERT_RERUNS.with(|c| c.get()) >= 10 {
        status = Status::NumericalFailure;
    }
    if (restored.is_some()
        && (matches!(status, Status::Optimal | Status::Infeasible)
            || (R::QUAD && s.art_exit && status == Status::Unbounded)))
        || (status == Status::Infeasible && s.rejected_end)
    {
        if restored.is_some() {
            PAIR_TOL.with(|c| c.set(1e-7));
        }
        // The controller (`0x154f43..0x155190`) after a dual that left a variable out (a singular supplied basis, or the repair after the retry limit):
        // the variable goes back to a bound (status -4 + 3) and the primal runs from the dual's basis; the dropped variable is priced again (a bound flip
        // or a pivot, each one an iteration).
        let mut popts = opts.clone();
        popts.method = crate::simplex::Method::Primal;
        // The dual's teacher-forced factorization schedule does not apply to the primal, which decides its own (work based).
        popts.refactor_schedule = None;
        if R::QUAD && crate::simplex::QUAD_PRIMALS.with(|c| c.get()) == 0 && s.iters > 0 {
            // The x87 start selects method 0, 2, or 1 for ctx+0x35c equal
            // to 1, 2, or at least 3. The post-factor decision can also move
            // method 0 to 1 when the dual made iterations since its snapshot.
            let increments = !(s.pert != 0 && !s.rejected_end && !s.art_exit);
            let c35c_at_start = C35C.with(|c| c.get()) + u32::from(increments);
            let fires = C178.with(|c| c.get()) > 0 && s.stall_it != s.iters as f64;
            popts.pricing0 = c35c_at_start == 1 && !fires;
        }
        if opts.perturb != 0.0 {
            popts.sw_carry = Some(crate::simplex::SwCarry {
                iters: s.iters as f64,
                degenerate: s.stall as i32,
                snap_degenerate: s.stall_snap as i32,
                snap_iters: s.stall_it,
                factorizations: s.refactorizations as u32,
                level: 0,
                pert: s.pstep * 0.05,
            });
        }
        let b = Basis {
            basic: s.head.clone(),
            at_upper: s.st.iter().map(|&t| t == UB).collect(),
        };
        crate::simplex::HANDOFF.with(|c| {
            c.set(Some(crate::simplex::Handoff {
                stale: Some(s.pair),
                thr: s.lu_thr,
                retries: s.retries as usize,
                pert: s.pert >= 2 && !s.rejected_end && !s.art_exit,
                pair_tol: PAIR_TOL.with(|c| c.get()),
                quad: R::QUAD,
                dual_level: s.level,
                shifted_exit: !R::QUAD && s.rejected_end && s.pert != 0,
                iters: (s.it_base + s.iters) as f64,
                snap: (s.it_base + s.stall_it as usize) as f64,
                c178: C178.with(|c| c.get()),
                excl: if R::QUAD && (s.pert == 0 || s.art_exit) {
                    s.repaired
                } else {
                    None
                },
            }))
        });
        // The controller scales the tolerance by a factor of `ctx+0x35c` (`0x154cf0..0x154d13`, store `0x15647a`): 1 up to 5, 0.5 for 6 and 7, 0.2 above (seed 2950:
        // `ctx+0x388` 1e-6, 5e-7, 2e-7), and 1 when the tolerance now in force is above the parameter's 1e-6 (`0x156464`; fx_midsing_1106 keeps 1e-3). `ctx+0x35c` is 1 at the
        // start and counts the dual passes ended with status -1004 (`0xe404a`, `0xf8691`): the primals the controller starts after a dual (fx_midsing_1496: 2, 3, 4, 5, 6 -> 5e-7,
        // 7); a primal rerun after the restore does not count. In a quad solve rose counts the quad primals started.
        let handoffs = DUAL_HANDOFFS.with(|c| {
            c.set(c.get() + 1);
            c.get()
        });
        if !(R::QUAD && s.pert != 0 && !s.rejected_end && !s.art_exit) {
            C35C.with(|c| c.set(c.get() + 1));
        }
        let c35c = if R::QUAD {
            C35C.with(|c| c.get())
        } else {
            handoffs + 1
        };
        let tol_factor = if dual_tol() > DUAL_TOL0 {
            1.0
        } else {
            match c35c {
                0..=5 => 1.0,
                6 | 7 => 0.5,
                _ => 0.2,
            }
        };
        set_dual_tol(LOOP_TOL.with(|c| c.get()) * tol_factor);
        crate::simplex::PRIMAL_SINGULAR.with(|c| c.set(None));
        crate::simplex::PRIMAL_START_SINGULAR.with(|c| c.set(false));
        crate::simplex::PRIMAL_FARKAS_SKIPPED.with(|c| c.set(0));
        crate::simplex::PRIMAL_REPAIRED.with(|c| c.set(false));
        crate::simplex::PRIMAL_SHIFT_END.with(|c| c.set(false));
        let primal_result = crate::simplex::solve(model, Some(&b), &popts);
        // The controller saves the tolerance left by the primal, undoing the
        // handoff scale (`0x1564a4`), for the next handoff.
        LOOP_TOL.with(|c| c.set(dual_tol() / tol_factor));
        let mut r = match primal_result {
            Err(SolveError::SingularBasis) => {
                // The primal's factorization has a column without a pivot (status -1002): the controller restarts the dual from the slack basis, like the dual's
                // own failed factorization, with the primal's iterations counted.
                let Some((pit, pfac, pthr)) = crate::simplex::PRIMAL_SINGULAR.with(|c| c.take())
                else {
                    return Err(SolveError::SingularBasis);
                };
                if s.retries < 5 {
                    // The primal LU singular handler increments the same counter
                    // used by the next dual ratio test (0xf58f7 / 0x248f99).
                    C178.with(|c| c.set(c.get() + 1));
                }
                if !R::QUAD
                    && s.rejected_end
                    && !crate::simplex::PRIMAL_START_SINGULAR.with(|c| c.take())
                {
                    // The double primal after a rejected end made its pivots and ended with status -1004 (`0xe404a`, seed 3981), then the dual restarts in quad.
                    C35C.with(|c| c.set(c.get() + 1));
                }
                if crate::simplex::PRIMAL_DRIFT.with(|c| c.take()) {
                    DRIFT_START.with(|c| *c.borrow_mut() = Some(b.clone()));
                }
                // The reference's iteration counter (`ctx+0x2a8`, the log timer's clock) runs five above the primal's pivots here (fx_midsing_1106: 9 per run).

                SINGULAR_BASIS.with(|c| *c.borrow_mut() = Some(b.clone()));
                SINGULAR.with(|c| {
                    c.set(Some((
                        s.iters + pit,
                        s.refactorizations + pfac,
                        s.lu_thr.max(pthr),
                        s.avg_touch,
                        s.last_cap,
                        s.inf_count,
                        false,
                        s.retries,
                        s.it_base + s.iters + pit,
                    )))
                });
                return Err(SolveError::SingularBasis);
            }
            other => other?,
        };
        let mut repaired =
            crate::simplex::PRIMAL_REPAIRED.with(|c| c.take()) && r.status == Status::Infeasible;
        let had_repair = repaired;
        let shift_end = crate::simplex::PRIMAL_SHIFT_END.with(|c| c.take());
        let mut restore_dual = false;
        let mut repair_count = 0usize;
        while repaired {
            repair_count += 1;
            if R::QUAD {
                // Each dropped-variable restore increments the controller's
                // counter. The first two retry primal; from the third it
                // switches back to dual, and the eleventh ends numerically
                // (`0x1551ac..0x15527c`, seed 4737).
                let count = PERT_RERUNS.with(|c| {
                    c.set(c.get() + 1);
                    c.get()
                });
                if count > 10 {
                    r.status = Status::NumericalFailure;
                    break;
                }
                if count >= 3 {
                    restore_dual = true;
                    break;
                }
            }
            // Restore the dropped variable and restart primal from its basis.
            crate::simplex::HANDOFF.with(|c| {
                c.set(Some(crate::simplex::Handoff {
                    stale: None,
                    thr: s.lu_thr,
                    retries: s.retries as usize + repair_count,
                    pert: false,
                    pair_tol: PAIR_TOL.with(|c| c.get()),
                    quad: R::QUAD,
                    dual_level: s.level,
                    shifted_exit: false,
                    iters: (s.it_base + s.iters) as f64,
                    snap: (s.it_base + s.stall_it as usize) as f64,
                    c178: C178.with(|c| c.get()),
                    excl: None,
                }))
            });
            let b2 = r.basis.clone();
            let r2 = crate::simplex::solve(model, Some(&b2), &popts)?;
            r.iterations += r2.iterations;
            r.refactorizations += r2.refactorizations;
            r.status = r2.status;
            r.farkas = r2.farkas;
            r.ray = r2.ray;
            r.objective = r2.objective;
            r.x = r2.x;
            r.dual = r2.dual;
            r.reduced_cost = r2.reduced_cost;
            r.basis = r2.basis;
            repaired = R::QUAD
                && crate::simplex::PRIMAL_REPAIRED.with(|c| c.take())
                && r.status == Status::Infeasible;
        }
        if opts.record_trace {
            trace.append(&mut r.trace);
            r.trace = trace;
        }
        r.iterations += s.iters + crate::simplex::EARLIER_ITERS.with(|c| c.take());
        r.refactorizations += s.refactorizations;
        if restore_dual {
            let carry2 = Carry {
                retries: s.retries + repair_count as u32,
                thr: s.lu_thr.max(0.5),
                fac: r.refactorizations,
                iters: r.iterations,
                avg_touch: 0.0,
                last_cap: 0,
                inf_count: 0,
                after_pivot: true,
                it_base: s.it_base + r.iterations,
                nref0: s.nref() + (r.refactorizations - s.refactorizations).max(1),
                pstep: None,
            };
            let mut r2 = solve_pass::<R>(model, Some(&r.basis), opts, &carry2)?;
            r2.iterations += carry2.iters;
            return Ok(r2);
        }
        if R::QUAD && crate::simplex::PRIMAL_STALL.with(|c| c.take()) {
            C35C.with(|c| c.set(c.get() + 1));
            // The x87 primal ended on its progress check (status -1004): the controller runs the quad dual again from the primal's basis, and the primal after it.
            let l = r.refactorizations - s.refactorizations;
            let carry2 = Carry {
                retries: s.retries,
                thr: s.lu_thr,
                fac: r.refactorizations,
                iters: r.iterations,
                avg_touch: 0.0,
                last_cap: 0,
                inf_count: 0,
                after_pivot: true,
                it_base: s.it_base + r.iterations,
                nref0: s.nref() + l.max(1),
                pstep: if s.pstep > 0.0 {
                    Some(maxsd(1e-7, s.pstep * 0.05))
                } else {
                    None
                },
            };
            let mut r2 = solve_pass::<R>(model, Some(&r.basis), opts, &carry2)?;
            r2.iterations += carry2.iters;
            return Ok(r2);
        }
        if r.status == Status::NumericalFailure {
            // The primal's phase 1 ratio test finds no row on a fresh factor: the reference ends the solve with status 12 (fx_midsing_1031), the controller does not
            // restart the dual.
            return Ok(r);
        }
        if status == Status::Infeasible
            && s.pert != 0
            && !had_repair
            && PERT_RERUNS.with(|c| c.get()) >= 10
        {
            // The controller's restore counter (`ebx` of `0x1551ac`) passes 10 at the eleventh restore: status 12 (`0x15553e`).
            r.status = Status::NumericalFailure;
            return Ok(r);
        }
        // A double primal that reaches an optimum ends the controller: do not restart the dual just because its earlier perturbation flag was set
        // (fx_e2_r03125 and fx_e2_r0625: the reference has no additional dual pass or factorization).
        if status == Status::Infeasible
            && s.pert != 0
            && !had_repair
            && (R::QUAD || r.status != Status::Optimal)
            && (!(R::QUAD && s.rejected_end) || shift_end)
            && PERT_RERUNS.with(|c| {
                // A rejected-end handoff is not a restore; the controller's
                // restore counter increments on the later status-3 exits.
                if !(R::QUAD && s.rejected_end) {
                    c.set(c.get() + 1);
                }
                c.get() <= 10
            })
        {
            C35C.with(|c| c.set(c.get() + 1));
            // The dual ended with cost shifts in force: the primal removes them and ends with status -1004 (its price scan finds nothing), and the
            // controller runs the dual again from the primal's basis (m5: one primal pivot, then the dual proves the infeasibility on the true costs).
            // Status -1004 switches methods without a restart (`0x155024`);
            // only -1002 calls the snapshot restart handler (seed 4918, 1106).
            let carry2 = Carry {
                retries: s.retries,
                thr: s.lu_thr,
                fac: r.refactorizations,
                iters: r.iterations,
                avg_touch: 0.0,
                last_cap: 0,
                inf_count: 0,
                after_pivot: true,
                it_base: s.it_base + r.iterations,
                nref0: {
                    // `ctx+0x2c0` after the primal: its start factorization counts one (the primal's own refactorizations at its end, two in rose, do not appear in the
                    // reference's `lu_in`: fx_midsing_m5, 1494, 1390, `dual_acc` captures). The x87 primal with shifts in force ends at its first failed scan with
                    // status -1004 and makes no such factorizations: every one counts (seed 2950: four in the reference, `ctx+0x2c0` 24 -> 28).
                    let l = r.refactorizations - s.refactorizations;
                    s.nref()
                        + if R::QUAD {
                            l.max(1)
                        } else {
                            l.saturating_sub(
                                2 - crate::simplex::PRIMAL_FARKAS_SKIPPED.with(|c| c.take()),
                            )
                        }
                },
                pstep: if R::QUAD && s.pstep > 0.0 {
                    Some(maxsd(1e-7, s.pstep * 0.05))
                } else {
                    None
                },
            };
            let mut r2 = solve_pass::<R>(model, Some(&r.basis), opts, &carry2)?;
            r2.iterations += carry2.iters;
            return Ok(r2);
        }
        return Ok(r);
    }
    if status == Status::Unbounded && s.art_exit {
        if R::QUAD && s.art_exit {
            C35C.with(|c| c.set(c.get() + 1));
        }
        // A request-4 artificial-bound exit returns -1004 to the controller, which runs the primal from this basis. Its completed status is final.
        let mut popts = opts.clone();
        popts.method = crate::simplex::Method::Primal;
        if opts.perturb != 0.0 {
            // The release profile: the primal goes on in method 2 with the dual's counters (`compare.py switch` call 0 of `dual-release`: iterations, the stall
            // counters, `ctx+0x2c0` one above the dual's factorizations, `ctx+0x340` = 1e-5 for the default `PerturbValue`).
            popts.sw_carry = Some(crate::simplex::SwCarry {
                iters: s.iters as f64,
                degenerate: s.stall as i32,
                snap_degenerate: s.stall_snap as i32,
                snap_iters: s.stall_it,
                factorizations: s.refactorizations as u32,

                level: 0,
                pert: s.pstep * 0.05,
            });
        }
        let from_dual = Basis {
            basic: s.head.clone(),
            at_upper: (0..n + m).map(|j| s.st[j] == UB && p.up[j] < BIG).collect(),
        };
        if R::QUAD && s.art_exit {
            // After a quad pass the reference's primal is the x87 one (`ctx+0x28` stays 1), whatever ended the dual (random LP seed 5247: the request-4
            // exit of a quad pass; the double primal found the ill-conditioned basis singular and the run ended with an error).
            if crate::simplex::QUAD_PRIMALS.with(|c| c.get()) == 0 && s.iters > 0 {
                popts.pricing0 = C35C.with(|c| c.get()) == 1;
            }
            crate::simplex::HANDOFF.with(|c| {
                c.set(Some(crate::simplex::Handoff {
                    stale: Some(s.pair),
                    thr: s.lu_thr,
                    retries: s.retries as usize,
                    // The quad dual driver ends with status 2 or 5 (an unbounded verdict here) clear the cost shifts `ctx+0x338` (`0x24a539`), so the x87 primal runs without them and
                    // ends with a refactorization and an optimal verdict (seed 6610) instead of the request-9 end that restores the dual (seeds 2133, 2950).
                    pert: false,
                    pair_tol: PAIR_TOL.with(|c| c.get()),
                    quad: true,
                    dual_level: s.level,
                    shifted_exit: false,
                    iters: (s.it_base + s.iters) as f64,
                    snap: (s.it_base + s.stall_it as usize) as f64,
                    c178: C178.with(|c| c.get()),
                    excl: if s.pert == 0 { s.repaired } else { None },
                }))
            });
        }
        let prres = crate::simplex::solve(model, Some(&from_dual), &popts).ok();
        if R::QUAD && s.art_exit {
            if let Some(pr) = prres
                .as_ref()
                .filter(|pr| pr.status == Status::Infeasible && pr.farkas.is_some())
            {
                // A completed phase-1 infeasibility verdict ends the controller;
                // a failed scan requesting a restart resumes below. The reference's
                // raw Farkas vector need not verify on ill-conditioned models.
                let mut r = pr.clone();
                r.iterations += s.iters;
                r.refactorizations += s.refactorizations;
                FARKAS_LEAVE.with(|c| c.set(None));
                return Ok(r);
            }
        }
        if R::QUAD && s.art_exit {
            if let Some(pr) = prres
                .as_ref()
                .filter(|pr| pr.status == Status::NumericalFailure)
            {
                // The x87 primal ended with status 12 (`0x280dd9`: no leaving row for the entering variable, again after the request-7 refactorization, with updates on
                // the factor): the controller ends the solve with it (seed 8179).
                let mut r = pr.clone();
                r.iterations += s.iters;
                r.refactorizations += s.refactorizations;
                return Ok(r);
            }
        }
        // The primal's unit step of the entering variable belongs to its own ray only, not to the dual's artificial-bound point.
        let prres_enter = crate::simplex::RAY_ENTER.with(|c| c.take());
        if s.art_exit
            && !prres
                .as_ref()
                .is_some_and(|pr| matches!(pr.status, Status::Optimal | Status::Unbounded))
        {
            // The primal gives up from the dual's basis: the dual runs on past the request-4 refactorization instead.
            NO_ART_EXIT.with(|c| c.set(true));
            let resume = prres.as_ref().map(|pr| pr.basis.clone());
            if R::QUAD
                && QUAD_RESUMES.with(|c| {
                    c.set(c.get() + 1);
                    c.get() <= 1
                })
            {
                // A quad pass resumes as a quad pass from the primal's basis (the x87 primal gave up at its first scan), the counters going on (random LP seed 5247).
                let pit = prres.as_ref().map_or(0, |pr| pr.iterations);
                let pfac = prres.as_ref().map_or(0, |pr| pr.refactorizations);
                let carry3 = Carry {
                    retries: s.retries,
                    thr: s.lu_thr,
                    fac: s.refactorizations + pfac,
                    iters: s.iters + pit,
                    avg_touch: 0.0,
                    last_cap: 0,
                    inf_count: 0,
                    after_pivot: false,
                    it_base: s.it_base + s.iters + pit,
                    nref0: s.nref() + pfac.max(1),
                    pstep: None,
                };
                let again = solve_pass::<R>(model, resume.as_ref().or(start), opts, &carry3);
                NO_ART_EXIT.with(|c| c.set(false));
                let mut again = again?;
                again.iterations += carry3.iters;
                return Ok(again);
            }
            let again = solve_frame(model, resume.as_ref().or(start), opts);
            NO_ART_EXIT.with(|c| c.set(false));
            return again;
        }
        if s.art_exit {
            if let Some(pr) = prres.as_ref().filter(|pr| pr.status == Status::Optimal) {
                // The variables left at an artificial bound were not a proof: the primal from the dual's basis found the optimum.
                let mut r = pr.clone();
                r.iterations += s.iters;
                r.refactorizations += s.refactorizations;
                return Ok(r);
            }
        }
        if let Some(mut r) = prres {
            r.iterations += s.iters;
            r.refactorizations += s.refactorizations;
            crate::simplex::RAY_ENTER.with(|c| c.set(prres_enter));
            return Ok(r);
        }
        return Err(SolveError::SingularBasis);
    }
    if status == Status::Unbounded {
        // Status 5 on the dual itself goes directly to the controller's final status store. Keep its raw artificial-bound ray and return the finite
        // point obtained by clearing the artificial bounds, without solving a second feasibility problem.
        s.zero_artificial = true;
        s.recompute_xb(0);
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
    let mut xr: Vec<(f64, f64)> = Vec::new();
    for j in 0..n + m {
        if R::QUAD {
            // The long-double value the quad pass holds (exact as two doubles): the reference's objective is computed from it, not from the rounded x.
            let v = if s.st[j] >= 0 {
                s.xb[s.st[j] as usize]
            } else {
                R::c(s.nonbasic_value(j))
            };
            let hi = v.to_f64();
            xr.push((hi, (v - R::c(hi)).to_f64()));
        }
        x[j] = if s.st[j] >= 0 {
            let v = s.xb[s.st[j] as usize];
            let (lo, up) = (s.p.lo[j], s.p.up[j]);
            // After the clean-up factorization the reference shifts the basic solution to zero lower bounds (`x - lo`, loop `0x10731f..0x1073d2`, run on
            // the FTRAN result; the threshold `[rsp+0x78]` is `ctx+0x370` = 1e30) and adds the bound back at the read-out: a rounding round trip for a
            // finite lower bound (not on the other exits: 80BAU3B is bit-equal without it). A variable with only an upper bound is mirrored in the
            // reference (`x' = -x`, lower bound `-up`), so its round trip uses the upper bound: `(-x + up) - up` negated is `(x - up) + up`. A free basic
            // variable stays as it is (bit-equal on fx_cleanrt_free_105, _1331, _2122: 19 to 36 nonzero free basic values each).
            if cleaned && lo > -1e30 && lo.is_finite() {
                ((v - R::c(lo)) + R::c(lo)).to_f64()
            } else if cleaned && up < 1e30 && up.is_finite() {
                ((v - R::c(up)) + R::c(up)).to_f64()
            } else {
                v.to_f64()
            }
        } else {
            s.nonbasic_value(j)
        };
    }
    s.cost.clone_from(&p.cost);
    let mut y: Vec<f64> = s.head.iter().map(|&j| p.cost[j]).collect();
    s.factor.btran(&mut y);
    let dtrue: Vec<f64> = s.true_d().iter().map(|v| v.to_f64()).collect();
    let nm = model.num_cols();
    let xs: Vec<f64> = x[..nm].to_vec();
    let at_upper: Vec<bool> = s.st.iter().map(|&t| t == UB).collect();
    QUAD_XR.with(|c| *c.borrow_mut() = if R::QUAD { Some(xr) } else { None });
    Ok(SolveResult {
        status,
        objective: objective_in::<R>(model, &xs),
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
