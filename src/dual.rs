//! Dual simplex of the reference's `dual-controlled` profile (`Method=1`, `Presolve=0`, `ScaleFlag=0`, `SimplexPricing=1`, `PerturbValue=0`).
//!
//! Rules: `reverse-engineer/specs/dual-simplex.md`. The reference keeps the reduced costs and the pivot row in a frame in which
//! the dual-feasible sign is positive: `d'_j = d_j` for a variable at its lower bound (status -1) and `-d_j` at its upper bound
//! (status -2), and the pivot row is negated for the same variables and when the leaving variable is above its bound, so every
//! ratio-test candidate has `alpha' < 0`. Basic variables are indexed by position (the row of the pivot-row BTRAN).

use crate::basis::Basis;
#[cfg(not(feature = "factor-shadow"))]
pub(crate) type Mirror<R> = crate::dual_factor::DualFactor<R>;
#[cfg(feature = "factor-shadow")]
pub(crate) type Mirror<R> = crate::dual_factor::shadow::Shadow<R>;
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

// These independent context fields survive method handoffs and ordinary restarts.
thread_local! {
    static PRIMAL_PERT_SIZE: std::cell::Cell<f64> = const { std::cell::Cell::new(0.0) };
    static DUAL_PERT_SIZE: std::cell::Cell<f64> = const { std::cell::Cell::new(0.0) };
}

pub(crate) fn reset_pert_sizes(value: f64) {
    PRIMAL_PERT_SIZE.with(|c| c.set(value));
    DUAL_PERT_SIZE.with(|c| c.set(value));
}

fn enter_pert_size(
    cell: &'static std::thread::LocalKey<std::cell::Cell<f64>>,
    primal: bool,
) -> f64 {
    let count = C35C.with(|c| c.get());
    let size = cell.with(|c| {
        let old = c.get();
        let size = if count > 1 && old > 1e-7 {
            maxsd(1e-7, old / 20.0)
        } else {
            old
        };
        c.set(size);
        size
    });
    C35C.with(|c| {
        c.set(if primal {
            count.wrapping_abs()
        } else if count == -1 {
            0
        } else {
            count
        })
    });
    size
}

fn enter_dual_pert_size() -> (f64, i64) {
    let left = if C35C.with(|c| c.get()) == -1 {
        1
    } else {
        1_000_000_000
    };
    (enter_pert_size(&DUAL_PERT_SIZE, false), left)
}

pub(crate) fn enter_primal_pert_size() -> f64 {
    enter_pert_size(&PRIMAL_PERT_SIZE, true)
}

pub(crate) fn primal_pert_size() -> f64 {
    PRIMAL_PERT_SIZE.with(|c| c.get())
}

pub(crate) fn set_primal_pert_size(value: f64) {
    PRIMAL_PERT_SIZE.with(|c| c.set(value));
}

// State held by the native simplex controller, outside either callback.
#[derive(Clone, Copy)]
struct CallbackController {
    selector: bool,
    result: i32,
    skip: bool,
    primal_carry: f64,
    dual_carry: f64,
}

thread_local! {
    static CALLBACK: std::cell::Cell<CallbackController> = const { std::cell::Cell::new(CallbackController {
        selector: true, result: 0, skip: false, primal_carry: PRIMAL_TOL, dual_carry: DUAL_TOL0,
    }) };
}

fn callback_reset() {
    CALLBACK.with(|c| {
        c.set(CallbackController {
            selector: true,
            result: 0,
            skip: false,
            primal_carry: primal_tol(),
            dual_carry: dual_tol(),
        })
    });
}

struct CallbackGuard {
    dual: bool,
    factor: f64,
    active: bool,
}

impl CallbackGuard {
    fn enter(dual: bool, opts: &SolveOptions) -> Self {
        let count = C35C.with(|c| c.get());
        let mut factor = if count <= 5 {
            1.0
        } else if count <= 7 {
            0.5
        } else {
            0.2
        };
        CALLBACK.with(|c| {
            let mut st = c.get();
            st.selector = dual;
            // Only ordinary callback dispatches are wired up here. Native
            // control-result and one-shot skip routing is not yet recovered.
            assert!(st.result == 0 && !st.skip, "unimplemented callback control route");
            if crate::dump::on("callback") {
                crate::dump::emit("callback", serde_json::json!({"event":"enter", "selector":i32::from(dual), "result":st.result, "skip":st.skip, "count":count, "c380":primal_tol(), "c388":dual_tol(), "carry":st.primal_carry, "total":crate::simplex::TOTAL_ITERS.with(|c| c.get())}));
            }
            let live = if dual { primal_tol() } else { dual_tol() };
            let parameter = if dual { opts.feasibility_tol } else { crate::simplex::reference_optimality_tol(opts) };
            if live > parameter { factor = 1.0; }
            if dual { set_primal_tol(st.primal_carry * factor); }
            else { set_dual_tol(st.dual_carry * factor); }
            c.set(st);
        });
        Self {
            dual,
            factor,
            active: true,
        }
    }

    fn finish(&mut self) {
        if !self.active {
            return;
        }
        CALLBACK.with(|c| {
            let mut st = c.get();
            if self.dual { st.primal_carry = primal_tol() / self.factor; }
            else { st.dual_carry = dual_tol() / self.factor; }
            if crate::dump::on("callback") {
                crate::dump::emit("callback", serde_json::json!({"event":"exit", "selector":i32::from(self.dual), "count":C35C.with(|c| c.get()), "c380":primal_tol(), "c388":dual_tol(), "carry":st.primal_carry, "total":crate::simplex::TOTAL_ITERS.with(|c| c.get())}));
            }
            c.set(st);
        });
        self.active = false;
    }
}

impl Drop for CallbackGuard {
    fn drop(&mut self) {
        self.finish();
    }
}

fn primal_callback(
    model: &LpModel,
    start: Option<&Basis>,
    opts: &SolveOptions,
) -> Result<SolveResult, SolveError> {
    let _callback = CallbackGuard::enter(false, opts);
    crate::simplex::solve(model, start, opts)
}

thread_local! {
    /// The dual feasibility tolerance (`ctx+0x388`) and the counter of failed dual residual checks (`ctx+0x174`); both outlive the restarts of one solve.
    static DUAL_TOL_CELL: std::cell::Cell<f64> = const { std::cell::Cell::new(DUAL_TOL0) };
    /// The primal feasibility tolerance in force (`ctx+0x380`); the primal's drift test raises it (up to 1e-3), the dual passes after it read the raised value.
    pub(crate) static PRIMAL_DRIFT_COUNT: std::cell::Cell<u32> = const { std::cell::Cell::new(0) };
    static PRIMAL_TOL_CELL: std::cell::Cell<f64> = const { std::cell::Cell::new(PRIMAL_TOL) };
    static RNG_CELL: std::cell::Cell<u64> = const { std::cell::Cell::new(101) };
    static PERT_RERUNS: std::cell::Cell<u32> = const { std::cell::Cell::new(0) };
    static RESID_FAILS: std::cell::Cell<u32> = const { std::cell::Cell::new(0) };
    /// Sum of dual infeasibilities from the last full recomputation (`ctx+0x2d8`).
    static DUAL_INFEAS_SUM: std::cell::Cell<f64> = const { std::cell::Cell::new(0.0) };
}

thread_local! {
    /// `ctx+0x35c`: 1 at the start, plus one for every run that ends with status -1004 (writers `0xf8691` double dual, `0xe404a` double primal, `0x24b98b` quad dual,
    /// `0x2371bd` x87 primal; `replay/gdb_c35c_trace.py`). The quad primals' method (`ctx+0x3ec`) and the switch of the first ones read it.
    pub(crate) static C35C: std::cell::Cell<i32> = const { std::cell::Cell::new(1) };
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
    static SNAP_PRIMAL: std::cell::Cell<bool> = const { std::cell::Cell::new(false) };
    static SINGULAR_FROM_PRIMAL: std::cell::Cell<bool> = const { std::cell::Cell::new(false) };
    static SINGULAR_TOL_ABORT: std::cell::Cell<bool> = const { std::cell::Cell::new(false) };
    pub(crate) static RESTORED_PRIMAL: std::cell::Cell<bool> = const { std::cell::Cell::new(false) };
    /// The pending flag `ctx+0x180` (`0xcef91..0xcef9b`): a quality-2 snapshot restore asks the restarted x87 primal for the infeasible-basic repair at its first basic-value recompute.
    pub(crate) static RESTORE_REPAIR: std::cell::Cell<bool> = const { std::cell::Cell::new(false) };
    /// The relative tolerance of the primal's pivot-pair check (`ctx+0x3a0`): 3e-6 from the parameters, 1e-7 from the first restore of a dropped variable on (`0x1551d4`).
    static PAIR_TOL: std::cell::Cell<f64> = const { std::cell::Cell::new(3e-6) };
    /// The primal's start basis when its drift test ended it (status -1002): the restart takes it (`0xcedb0`, seed 2717).
    static DRIFT_START: std::cell::RefCell<Option<Basis>> = const { std::cell::RefCell::new(None) };
    static SNAP_C40: std::cell::Cell<f64> = const { std::cell::Cell::new(-100.0) };
    /// The tolerance reset of the third restart (`0xcf29b`, once per restart) already happened in the primal run that the restart started first; the dual pass that
    /// follows it must not reset again (witness_primal_r31_double4: the reference's dual pass after the primal runs at 0.005, not 1e-6).
    static TOL_RESET_DONE: std::cell::Cell<bool> = const { std::cell::Cell::new(false) };
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

/// A standalone primal starts with a fresh snapshot and progress clock.
pub(crate) fn reset_primal_snapshot() {
    SNAP.with(|c| *c.borrow_mut() = None);
    SNAP_QUALITY.with(|c| c.set(0));
    SNAP_C40.with(|c| c.set(-100.0));
    set_primal_tol(PRIMAL_TOL);
    PRIMAL_DRIFT_COUNT.with(|c| c.set(0));
}

/// The primal run of a double dual's artificial-bound exit (no `Handoff`) starts without clearing the controller's snapshot and its timer (`ctx+0x40`): the reference
/// takes no snapshot in it and the restart after its singular start restores the dual's own (witness_primal_r31_double4).
pub(crate) fn reset_primal_snapshot_keep() {
    set_primal_tol(PRIMAL_TOL);
    PRIMAL_DRIFT_COUNT.with(|c| c.set(0));
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

pub(crate) fn has_restart_snapshot() -> bool {
    SNAP_QUALITY.with(|c| c.get() > 0)
}

/// Save the primal's current status at its progress check (`0x1422a5`, `0x27e75a`,
/// `0xd40b0` mode 0). Quality uses primal infeasibility, not the dual sum.
pub(crate) fn primal_snapshot(
    basic: &[usize],
    at_upper: &[bool],
    infeas: f64,
    retries: usize,
    excluded: bool,
) {
    let quality = if infeas > 0.0 { 1 } else { 2 };
    if retries >= 5 || excluded || SNAP_QUALITY.with(|c| quality < c.get()) {
        return;
    }
    SNAP_QUALITY.with(|c| c.set(quality));
    SNAP_PRIMAL.with(|c| c.set(true));
    SNAP.with(|c| *c.borrow_mut() = Some((basic.to_vec(), at_upper.to_vec())));
    #[cfg(test)]
    PRIMAL_SNAPSHOT_TEST.with(|c| {
        if let Some(records) = c.borrow_mut().as_mut() {
            records.push((basic.to_vec(), at_upper.to_vec(), quality));
        }
    });
}

#[cfg(test)]
thread_local! {
    pub(crate) static PRIMAL_SUM_TEST: std::cell::RefCell<Option<Vec<(bool, usize, usize, f64)>>> = const { std::cell::RefCell::new(None) };
    pub(crate) static PRIMAL_SNAPSHOT_TEST: std::cell::RefCell<Option<Vec<(Vec<usize>, Vec<bool>, u32)>>> = const { std::cell::RefCell::new(None) };
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
pub(crate) fn lcg_int(state: &mut u64) -> u32 {
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

/// The row weights of `Dual::reset_pricing` before the cap, without materialized lists: the forward lists are the factor's L transpose read
/// backwards (the reference's descending order), the backward lists the factor's numeric U rows in their stored order. Both are those of the
/// fresh factor `f` (no update since `fresh`); the lists of the previous implementation (CSRs built from the structure and the basis columns)
/// are rebuilt and compared in debug and shadow builds.
fn pricing_weights<R: Real>(
    m: usize,
    st: &crate::lu_order::Structure<R>,
    cols: crate::lu_order::Cols<'_>,
    f: &Mirror<R>,
) -> Vec<f64> {
    #[cfg(any(debug_assertions, feature = "factor-shadow"))]
    check_reset_lists(m, st, cols, f);
    #[cfg(not(any(debug_assertions, feature = "factor-shadow")))]
    let _ = (st, cols);
    let n = m as f64;
    let mut w = vec![1.0f64; m];
    for k in 0..m {
        for &j in f.reset_forward(k).iter().rev() {
            let j = j as usize;
            w[k] = (w[k] + w[j]) - (w[j] * w[k]) / n;
        }
    }
    for k in (0..m).rev() {
        for &j in f.reset_backward(k) {
            let j = j as usize;
            w[k] = (w[k] + w[j]) - (w[j] * w[k]) / n;
        }
    }
    w
}

/// The lists the reset used to build (forward CSR deduplicated, backward filtered to later steps) against the factor's: same entries in the
/// same order, so the sequential sums of `pricing_weights` see identical terms.
#[cfg(any(debug_assertions, feature = "factor-shadow"))]
fn check_reset_lists<R: Real>(
    m: usize,
    st: &crate::lu_order::Structure<R>,
    cols: crate::lu_order::Cols<'_>,
    f: &Mirror<R>,
) {
    let (fwd, bwd, _) = pricing_weights_lists(m, st, cols);
    for k in 0..m {
        let got: Vec<usize> = f.reset_forward(k).iter().rev().map(|&j| j as usize).collect();
        assert_eq!(got, fwd[k], "reset forward list of step {k}");
        let got: Vec<usize> = f.reset_backward(k).iter().map(|&j| j as usize).collect();
        assert_eq!(got, bwd[k], "reset backward list of step {k}");
    }
}

/// The lists and weights of the reference form of `Dual::reset_pricing` (the `p68` dump carries the lists).
fn pricing_weights_lists<R: Real>(
    m: usize,
    st: &crate::lu_order::Structure<R>,
    cols: crate::lu_order::Cols<'_>,
) -> (Vec<Vec<usize>>, Vec<Vec<usize>>, Vec<f64>) {
    let mut step_of_row = vec![usize::MAX; m];
    let mut step_of_pos = vec![usize::MAX; m];
    for (s, &(r, c)) in st.pivots.iter().enumerate() {
        step_of_row[r] = s;
        step_of_pos[c] = s;
    }
    // Forward lists (`0xf48fd` structure +0x70/+0x80): the earlier steps whose L column holds the row of step k, descending.
    let mut fwd: Vec<Vec<usize>> = vec![Vec::new(); m];
    for e in 0..m {
        for &r in st.l_rows(e) {
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
                st.u_row(k)
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
    (fwd, bwd, w)
}

/// The cap of `Dual::reset_pricing`: the capped weights `min(0.1 m + 100, w + 100)`, or (total above `m (0.1 m + 100)`, or type 2 already)
/// the constant cap with type 2. Returns the new pricing type.
fn cap_pricing_weights(ptype: u8, m: usize, w: &mut [f64]) -> u8 {
    let n = m as f64;
    let total: f64 = w.iter().sum();
    let cap = 0.1 * n + 100.0;
    if ptype == 2 || total > n * cap {
        w.iter_mut().for_each(|x| *x = cap);
        2
    } else {
        for x in w.iter_mut() {
            let v = *x + 100.0;
            *x = if cap < v { cap } else { v };
        }
        ptype
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

// The multiply and comparison stay in the working precision, including x87.
fn primal_drift_excess<R: Real>(residual: R, violation: R) -> bool {
    residual * R::c(100.0) > violation
}

fn primal_drift_action(
    count: u32,
    tol: f64,
    quad: bool,
    automatic_quad: bool,
) -> (f64, bool, bool) {
    if count <= 2 {
        (tol, false, false)
    } else if tol < 1e-3 {
        (tol * 10.0, true, false)
    } else {
        (tol, false, !quad && automatic_quad)
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

/// The independent LU is only read by the stall check and deferred-free setup.
/// Its base build is lazy: `Dual::refactor` stores only the head (the columns the
/// eager code factored) and the first reading call builds `Factor::with_floor`.
/// Its pivot order and dropping differ from the mirror, so mirror success does
/// not prove the build succeeds; a failed build is reported by the reading call
/// as `SingularBasis`. Updates cannot fail and are replayed, with their original
/// representation, on demand.
struct AuxiliaryFactor {
    base: Option<Factor>,
    head: Vec<usize>,
    /// Absolute pivot floor of the lazy build (`Factor::with_floor`).
    floor: f64,
    pending: Vec<AuxiliaryUpdate>,
    /// Shadow build: `base` is the eager build of `Dual::refactor`; the first read builds the lazy one and compares.
    #[cfg(feature = "factor-shadow")]
    verify: bool,
}

enum AuxiliaryUpdate {
    Column(usize, Column),
    Dense(usize, Vec<f64>),
}

impl AuxiliaryFactor {
    #[cfg(test)]
    fn new(base: Factor) -> Self {
        Self {
            base: Some(base),
            head: Vec::new(),
            floor: 0.0,
            pending: Vec::new(),
            #[cfg(feature = "factor-shadow")]
            verify: false,
        }
    }

    fn lazy(head: Vec<usize>) -> Self {
        Self::lazy_floor(head, 0.0)
    }

    fn lazy_floor(head: Vec<usize>, floor: f64) -> Self {
        Self {
            base: None,
            head,
            floor,
            pending: Vec::new(),
            #[cfg(feature = "factor-shadow")]
            verify: false,
        }
    }

    fn build(m: usize, head: &[usize], p: &Internal, floor: f64) -> Result<Factor, usize> {
        let cols: Vec<&Column> = head.iter().map(|&j| &p.cols[j]).collect();
        Factor::with_floor(m, &cols, floor)
    }

    fn update_column<R: Real>(&mut self, r: usize, column: Vec<(usize, R)>) {
        self.pending.push(AuxiliaryUpdate::Column(
            r, column.into_iter().map(|(i, v)| (i, v.to_f64())).collect(),
        ));
    }

    fn update(&mut self, r: usize, alpha: Vec<f64>) {
        self.pending.push(AuxiliaryUpdate::Dense(r, alpha));
    }

    fn materialize(&mut self, p: &Internal) -> Result<&Factor, SolveError> {
        let floor = self.floor;
        self.materialize_with(|head| Self::build(p.m, head, p, floor))
    }

    fn materialize_with(&mut self, build: impl Fn(&[usize]) -> Result<Factor, usize>) -> Result<&Factor, SolveError> {
        if self.base.is_none() {
            self.base = Some(build(&self.head).map_err(|_| SolveError::SingularBasis)?);
        }
        #[cfg(feature = "factor-shadow")]
        if std::mem::take(&mut self.verify) {
            let lazy = build(&self.head).expect("shadow: lazy auxiliary build failed where the eager one succeeded");
            let eager = self.base.as_ref().unwrap();
            let m = self.head.len();
            let probe: Vec<f64> = (0..m).map(|i| 1.0 + 0.37 * i as f64).collect();
            let bits = |v: &[f64]| v.iter().map(|x| x.to_bits()).collect::<Vec<_>>();
            let (mut a, mut b) = (probe.clone(), probe.clone());
            eager.ftran(&mut a);
            lazy.ftran(&mut b);
            assert_eq!(bits(&a), bits(&b), "shadow: lazy auxiliary ftran differs from the eager build");
            let (mut a, mut b) = (probe.clone(), probe);
            eager.btran(&mut a);
            lazy.btran(&mut b);
            assert_eq!(bits(&a), bits(&b), "shadow: lazy auxiliary btran differs from the eager build");
        }
        let base = self.base.as_mut().unwrap();
        for update in self.pending.drain(..) {
            match update {
                AuxiliaryUpdate::Column(r, column) => base.update_column(r, column),
                AuxiliaryUpdate::Dense(r, alpha) => base.update(r, alpha),
            }
        }
        Ok(self.base.as_ref().unwrap())
    }

    fn ftran(&mut self, p: &Internal, rhs: &mut Vec<f64>) -> Result<(), SolveError> {
        self.materialize(p)?.ftran(rhs);
        Ok(())
    }

    fn btran(&mut self, p: &Internal, rhs: &mut Vec<f64>) -> Result<(), SolveError> {
        self.materialize(p)?.btran(rhs);
        Ok(())
    }
}

#[cfg(test)]
#[test]
fn auxiliary_updates_are_lazy_and_bit_identical() {
    fn check<R: Real>(values: &[f64]) {
        let m = 7;
        let cols: Vec<Column> = (0..m).map(|i| vec![(i, 2.0), ((i + 1) % m, 0.25)]).collect();
        let mut eager = Factor::new(m, &cols).unwrap();
        let mut lazy = AuxiliaryFactor::new(Factor::new(m, &cols).unwrap());
        for round in 0..32 {
            let r = round % m;
            let mut column: Vec<_> = (0..round % (m + 2))
                .map(|i| ((m - i % m - 1), R::c(values[(round + i) % values.len()])))
                .collect();
            column.push((r, R::c(1.25)));
            eager.update_column(r, column.clone());
            lazy.update_column(r, column);
            let mut alpha: Vec<_> = (0..m).map(|i| values[(round + i + 2) % values.len()]).collect();
            alpha[r] = 1.25;
            eager.update(r, alpha.clone());
            lazy.update(r, alpha);
            assert_eq!(lazy.base.as_ref().unwrap().num_updates(), (round / 3 * 3) * 2);
            if round % 3 == 2 || round == 31 {
                let bits = |v: Vec<f64>| v.into_iter().map(f64::to_bits).collect::<Vec<_>>();
                let rhs: Vec<_> = (0..m).map(|i| values[(round + i) % values.len()]).collect();
                let (mut a, mut b) = (rhs.clone(), rhs.clone());
                eager.ftran(&mut a);
                lazy.materialize_with(|_| unreachable!()).unwrap().ftran(&mut b);
                assert_eq!(bits(a), bits(b));
                let (mut a, mut b) = (rhs.clone(), rhs);
                eager.btran(&mut a);
                lazy.materialize_with(|_| unreachable!()).unwrap().btran(&mut b);
                assert_eq!(bits(a), bits(b));
                assert!(lazy.pending.is_empty());
                assert_eq!(lazy.base.as_ref().unwrap().nnz_eta, eager.nnz_eta);
            }
        }
        lazy.update(0, vec![1.0; m]);
        // Installing a fresh base discards all updates, including unread ones.
        lazy = AuxiliaryFactor::new(Factor::new(m, &cols).unwrap());
        assert!(lazy.pending.is_empty());
        assert_eq!(lazy.base.as_ref().unwrap().num_updates(), 0);
    }
    let values = [0.0, -0.0, 1.25, -0.75, f64::from_bits(1), f64::INFINITY,
        f64::from_bits(0x7ff8000000000123), f64::from_bits(0xfff8000000000456)];
    for values in [&values[..5], &values[..]] {
        check::<f64>(values);
        check::<crate::f80::F80>(values);
    }
}

#[cfg(test)]
#[test]
fn auxiliary_lazy_build_matches_eager_build() {
    fn check<R: Real>(values: &[f64]) {
        let m = 7;
        let all: Vec<Column> = (0..2 * m).map(|j| vec![(j % m, 2.0), ((j + 3) % m, 0.25), ((j + 1) % m, -0.5)]).collect();
        let head: Vec<usize> = (0..m).map(|i| (i * 2 + 1) % (2 * m)).collect();
        let build = |h: &[usize]| {
            let cols: Vec<&Column> = h.iter().map(|&j| &all[j]).collect();
            Factor::with_floor(m, &cols, 0.0)
        };
        let mut eager = build(&head).unwrap();
        let mut lazy = AuxiliaryFactor::lazy(head.clone());
        assert!(lazy.base.is_none());
        let bits = |v: Vec<f64>| v.into_iter().map(f64::to_bits).collect::<Vec<_>>();
        for round in 0..32 {
            let r = round % m;
            let mut column: Vec<_> = (0..round % (m + 2))
                .map(|i| ((m - i % m - 1), R::c(values[(round + i) % values.len()])))
                .collect();
            column.push((r, R::c(1.25)));
            eager.update_column(r, column.clone().into_iter().map(|(i, v)| (i, v.to_f64())).collect());
            lazy.update_column(r, column);
            let mut alpha: Vec<_> = (0..m).map(|i| values[(round + i + 2) % values.len()]).collect();
            alpha[r] = 1.25;
            eager.update(r, alpha.clone());
            lazy.update(r, alpha);
        }
        let rhs: Vec<_> = (0..m).map(|i| values[i % values.len()] + 1.0).collect();
        let (mut a, mut b) = (rhs.clone(), rhs.clone());
        eager.ftran(&mut a);
        lazy.materialize_with(&build).unwrap().ftran(&mut b);
        assert_eq!(bits(a), bits(b));
        let (mut a, mut b) = (rhs.clone(), rhs);
        eager.btran(&mut a);
        lazy.materialize_with(&build).unwrap().btran(&mut b);
        assert_eq!(bits(a), bits(b));
        // A failed lazy build is reported by the reading call as `SingularBasis`.
        let mut bad = AuxiliaryFactor::lazy(head);
        assert!(matches!(bad.materialize_with(|_| Err(0)), Err(SolveError::SingularBasis)));
    }
    let values = [0.0, -0.0, 1.25, -0.75, 3.5, f64::from_bits(1)];
    check::<f64>(&values);
    check::<crate::f80::F80>(&values);
}

/// The row-wise copy of the nonbasic structurals as flat pools: row `i` is `idx/val[start[i]..start[i] + len[i]]` with `cap[i] - len[i]` spare slots.
/// The per-pivot edits keep the in-row order of a `Vec<(usize, f64)>` (`swap_remove` moves the last entry into the hole; `push` appends): the order of
/// the scatter is the order of the floating-point sums. A full row moves to the end of the pool; the pool is compacted when most of it is dead.
struct FlatRows {
    start: Vec<u32>,
    len: Vec<u32>,
    cap: Vec<u32>,
    idx: Vec<u32>,
    val: Vec<f64>,
    /// Slots of the pool owned by no row (left behind by moved rows).
    dead: usize,
    /// Exclusive bound of every stored index (`n + m`).
    nv: usize,
    n: usize,
    /// Debug check: the same rows as `Vec<Vec>`, edited with the same operations.
    #[cfg(debug_assertions)]
    mirror: Vec<Vec<(usize, f64)>>,
}

impl FlatRows {
    fn new(m: usize, n: usize) -> Self {
        FlatRows {
            start: vec![0; m],
            len: vec![0; m],
            cap: vec![0; m],
            idx: Vec::new(),
            val: Vec::new(),
            dead: 0,
            nv: n + m,
            n,
            #[cfg(debug_assertions)]
            mirror: vec![Vec::new(); m],
        }
    }

    /// The rows of `rows` (indices below `n + m`) with spare capacity.
    fn from_rows(rows: &[Vec<(usize, f64)>], n: usize) -> Self {
        let m = rows.len();
        let mut f = FlatRows::new(m, n);
        assert!(n + m < u32::MAX as usize / 2);
        let total: usize = rows.iter().map(|r| r.len() + (r.len() / 2).max(4)).sum();
        assert!(total < u32::MAX as usize / 2);
        f.idx.reserve(total);
        f.val.reserve(total);
        for (i, r) in rows.iter().enumerate() {
            let cap = r.len() + (r.len() / 2).max(4);
            f.start[i] = f.idx.len() as u32;
            f.len[i] = r.len() as u32;
            f.cap[i] = cap as u32;
            for &(j, a) in r {
                assert!(j < n, "slack index {j} in a pivot row (n = {n})");
                f.idx.push(j as u32);
                f.val.push(a);
            }
            f.idx.resize(f.idx.len() + cap - r.len(), 0);
            f.val.resize(f.val.len() + cap - r.len(), 0.0);
        }
        #[cfg(debug_assertions)]
        {
            f.mirror = rows.to_vec();
        }
        f.check();
        f
    }

    #[inline(always)]
    fn row_len(&self, i: usize) -> usize {
        self.len[i] as usize
    }

    #[inline(always)]
    fn row(&self, i: usize) -> (&[u32], &[f64]) {
        let (s, l) = (self.start[i] as usize, self.len[i] as usize);
        (&self.idx[s..s + l], &self.val[s..s + l])
    }

    /// The first position of variable `j` in row `i`.
    fn position(&self, i: usize, j: usize) -> Option<usize> {
        self.row(i).0.iter().position(|&x| x as usize == j)
    }

    /// Replace entry `k` of row `i`.
    fn set(&mut self, i: usize, k: usize, j: usize, a: f64) {
        assert!(k < self.row_len(i) && j < self.n, "slack index {j} in a pivot row (n = {})", self.n);
        let s = self.start[i] as usize + k;
        self.idx[s] = j as u32;
        self.val[s] = a;
        #[cfg(debug_assertions)]
        {
            self.mirror[i][k] = (j, a);
        }
        self.check();
    }

    /// Remove entry `k` of row `i`; the last entry takes its place.
    fn swap_remove(&mut self, i: usize, k: usize) {
        let l = self.row_len(i);
        assert!(k < l);
        let s = self.start[i] as usize;
        self.idx[s + k] = self.idx[s + l - 1];
        self.val[s + k] = self.val[s + l - 1];
        self.len[i] -= 1;
        #[cfg(debug_assertions)]
        {
            self.mirror[i].swap_remove(k);
        }
        self.check();
    }

    /// Append `(j, a)` to row `i`.
    fn push(&mut self, i: usize, j: usize, a: f64) {
        assert!(j < self.n, "slack index {j} in a pivot row (n = {})", self.n);
        let l = self.row_len(i);
        if l == self.cap[i] as usize {
            self.grow(i);
        }
        let s = self.start[i] as usize + l;
        self.idx[s] = j as u32;
        self.val[s] = a;
        self.len[i] += 1;
        #[cfg(debug_assertions)]
        {
            self.mirror[i].push((j, a));
        }
        self.check();
    }

    /// Move the full row `i` to the end of the pool with doubled capacity (the pool is compacted first when most of it is dead).
    fn grow(&mut self, i: usize) {
        if self.dead > 4096 && self.dead * 2 > self.idx.len() {
            self.compact();
        }
        let l = self.row_len(i);
        let cap = (2 * l).max(l + 4);
        let (s, new_s) = (self.start[i] as usize, self.idx.len());
        assert!(new_s + cap < u32::MAX as usize / 2);
        self.idx.extend_from_within(s..s + l);
        self.val.extend_from_within(s..s + l);
        self.idx.resize(new_s + cap, 0);
        self.val.resize(new_s + cap, 0.0);
        self.dead += self.cap[i] as usize;
        self.start[i] = new_s as u32;
        self.cap[i] = cap as u32;
    }

    /// Rewrite the pool with the live rows only (same order, same capacities).
    fn compact(&mut self) {
        let m = self.start.len();
        let live: usize = self.cap.iter().map(|&c| c as usize).sum();
        let (mut idx, mut val) = (Vec::with_capacity(live), Vec::with_capacity(live));
        for i in 0..m {
            let (s, l, c) = (self.start[i] as usize, self.len[i] as usize, self.cap[i] as usize);
            self.start[i] = idx.len() as u32;
            idx.extend_from_slice(&self.idx[s..s + l]);
            val.extend_from_slice(&self.val[s..s + l]);
            idx.resize(idx.len() + c - l, 0);
            val.resize(val.len() + c - l, 0.0);
        }
        self.idx = idx;
        self.val = val;
        self.dead = 0;
    }

    /// Debug builds: the flat rows equal the `Vec<Vec>` mirror, entry by entry and bit by bit.
    #[inline(always)]
    fn check(&self) {
        #[cfg(debug_assertions)]
        {
            assert_eq!(self.mirror.len(), self.start.len());
            for (i, r) in self.mirror.iter().enumerate() {
                let (ix, vl) = self.row(i);
                assert!(self.len[i] <= self.cap[i]);
                assert_eq!(ix.len(), r.len(), "flat row {i} length");
                for (k, &(j, a)) in r.iter().enumerate() {
                    assert!(ix[k] as usize == j && vl[k].to_bits() == a.to_bits(), "flat row {i} entry {k}");
                }
            }
        }
    }
}

/// `d_j = c_j - y a_j` for the nonbasic `j` (`st[j] < 0`). The reference accumulates row-wise (`0x136b00`): rows in ascending order, skipping
/// |y_i| <= 1e-13 (`ctx+0x3a8`). `p.rwise[i]` lists the entries of row `i` by ascending column in the entry order of `p.cols`, so every column
/// sees its rows in ascending order.
fn true_d_rowwise<R: Real>(p: &Internal, st: &[i32], cost: &[f64], cost_lo: &[f64], y: &[R]) -> Vec<R> {
    let mut d = vec![R::ZERO; p.n + p.m];
    for j in 0..p.n + p.m {
        if st[j] < 0 {
            d[j] = R::c(cost[j]) + R::c(cost_lo[j]);
        }
    }
    for (i, row) in p.rwise.iter().enumerate() {
        let yi = y[i];
        if yi.abs() > R::c(DROP_TOL) {
            for &(j, a) in row {
                if st[j] < 0 {
                    d[j] -= yi * R::c(a);
                }
            }
        }
    }
    d
}

/// The retained workspace of the row product: `acc` is +0.0 and `bits` zero everywhere between calls.
struct PivotWs<R: Real> {
    acc: Vec<R>,
    /// The touched flags, one bit per variable.
    bits: Vec<u64>,
    /// Every index touched by the scatter, in first-touch order (capacity `n + m`).
    all: Vec<u32>,
    /// The coefficient of the first entry of each slack column.
    slack_coef: Vec<f64>,
    /// Result of the last call: the position in the returned row of the first entry of the leaving variable (it is only ever a slack entry).
    leave_pos: Option<usize>,
    /// Result of the last call: some variable occurs twice in the returned row (a repeated slack index in `rho`), so positions are not unique.
    slack_dup: bool,
    /// The `(start, len)` of each row of `rho` in the pool of `rrows`, found by the size pass and read again by the scatter.
    spans: Vec<(u32, u32)>,
}

impl<R: Real> Default for PivotWs<R> {
    fn default() -> Self {
        PivotWs {
            acc: Vec::new(),
            bits: Vec::new(),
            all: Vec::new(),
            slack_coef: Vec::new(),
            leave_pos: None,
            slack_dup: false,
            spans: Vec::new(),
        }
    }
}

/// The row product (see `Dual::pivot_row`) on a retained workspace; every operation and its order as in `pivot_row_ref`.
#[allow(clippy::too_many_arguments)]
fn pivot_row_ws<R: Real>(
    ws: &mut PivotWs<R>,
    rrows: &FlatRows,
    slack_nb: &[bool],
    st: &[i32],
    cols: &[Vec<(usize, f64)>],
    n: usize,
    m: usize,
    rho: &[(usize, R)],
    leave: usize,
    above: bool,
) -> (Vec<(usize, R)>, usize, usize, usize) {
    if ws.acc.len() != n + m {
        ws.acc = vec![R::ZERO; n + m];
        ws.bits = vec![0; (n + m).div_ceil(64)];
        ws.all = Vec::with_capacity(n + m);
        ws.slack_coef = (0..m).map(|i| cols[n + i][0].1).collect();
    }
    assert!(rrows.nv == n + m && rrows.start.len() == m && st.len() >= n + m && slack_nb.len() >= m && ws.slack_coef.len() == m);
    // Local owners keep the pointers and lengths in registers across the scatter (stores through `ws` could alias them).
    let mut acc_v = std::mem::take(&mut ws.acc);
    let mut bits_v = std::mem::take(&mut ws.bits);
    let mut all_v = std::mem::take(&mut ws.all);
    let acc = &mut acc_v[..];
    let bits = &mut bits_v[..];
    let slack_coef = &ws.slack_coef[..];
    // The row product runs over the rows of `rho` unless their entries (`0x134aa0` counts them) outnumber the structurals n by 1.8 (`0xc5355`, `ctx+0x68`);
    // then it is columnwise, in index order.
    // The size pass also checks every index of `rho` against `m` (the unchecked reads below rely on it) and keeps the row spans for the scatter.
    let mut spans_v = std::mem::take(&mut ws.spans);
    spans_v.clear();
    spans_v.reserve(rho.len());
    assert!(rrows.len.len() == m && rrows.idx.len() == rrows.val.len());
    let mut est = 0usize;
    {
        let spans_p = spans_v.as_mut_ptr();
        for (k, e) in rho.iter().enumerate() {
            // One check per entry; the reads of `start`/`len` (both of length m) need no more.
            assert!(e.0 < m, "pivot_row_ws: row index {} out of range (m = {m})", e.0);
            let (s, l) = (rrows.start[e.0], rrows.len[e.0]);
            debug_assert!(s as usize + l as usize <= rrows.idx.len());
            est += l as usize;
            // SAFETY: `k < rho.len()`, the reserved capacity.
            unsafe { spans_p.add(k).write((s, l)) };
        }
        // SAFETY: `rho.len()` spans were written.
        unsafe { spans_v.set_len(rho.len()) };
    }
    let sparse = rho.len() == 1 || !((n as f64) < 1.8 * est as f64);
    all_v.clear();
    all_v.reserve(n + m);
    let all_p = all_v.as_mut_ptr();
    let mut cnt = 0usize;
    // `acc`, `bits` and `all` hold `n + m`, `ceil((n + m) / 64)` and at least `n + m` slots; every stored index is below `n + m` (checked when stored),
    // and the first-touch bit makes the entries of `all` distinct.
    macro_rules! step {
        ($rv:expr, $j:expr, $a:expr) => {{
            let j = $j as usize;
            debug_assert!(j < acc.len());
            // SAFETY: j < n + m, see above.
            unsafe {
                let p = acc.get_unchecked_mut(j);
                *p = *p + $rv * R::c($a);
                let w = bits.get_unchecked_mut(j >> 6);
                let b = 1u64 << (j & 63);
                if *w & b == 0 {
                    *w |= b;
                    *all_p.add(cnt) = j as u32;
                    cnt += 1;
                }
            }
        }};
    }
    // The same update without the first-touch tracking (the columnwise route finds the touched slots by scanning `acc`).
    macro_rules! step_untracked {
        ($rv:expr, $j:expr, $a:expr) => {{
            let j = $j as usize;
            debug_assert!(j < n);
            // SAFETY: j < n (every stored index of `rrows` is, asserted above).
            unsafe {
                let p = acc.get_unchecked_mut(j);
                *p = *p + $rv * R::c($a);
            }
        }};
    }
    // SAFETY (the spans): `start + len <= idx.len() = val.len()` is an invariant of `FlatRows` (`check`).
    let (idx_p, val_p) = (rrows.idx.as_ptr(), rrows.val.as_ptr());
    macro_rules! scatter {
        ($step:ident) => {
            for (&(_, rv), &(s, l)) in rho.iter().zip(spans_v.iter()) {
                // SAFETY: the span is inside the pools (invariant above), so every read below is.
                let (ix, vl) = unsafe { (idx_p.add(s as usize), val_p.add(s as usize)) };
                let k = l as usize;
                let mut t = 0;
                while t + 4 <= k {
                    // The loads of four entries; the updates follow in order.
                    let (j0, j1, j2, j3) = unsafe { (*ix.add(t), *ix.add(t + 1), *ix.add(t + 2), *ix.add(t + 3)) };
                    let (a0, a1, a2, a3) = unsafe { (*vl.add(t), *vl.add(t + 1), *vl.add(t + 2), *vl.add(t + 3)) };
                    $step!(rv, j0, a0);
                    $step!(rv, j1, a1);
                    $step!(rv, j2, a2);
                    $step!(rv, j3, a3);
                    t += 4;
                }
                // The last 0-3 entries, in order.
                let r = k - t;
                if r >= 1 {
                    let (j0, a0) = unsafe { (*ix.add(t), *vl.add(t)) };
                    $step!(rv, j0, a0);
                    if r >= 2 {
                        let (j1, a1) = unsafe { (*ix.add(t + 1), *vl.add(t + 1)) };
                        $step!(rv, j1, a1);
                        if r == 3 {
                            let (j2, a2) = unsafe { (*ix.add(t + 2), *vl.add(t + 2)) };
                            $step!(rv, j2, a2);
                        }
                    }
                }
            }
        };
    }
    // The columnwise route (`!sparse`) outputs the touched structurals in index order and needs no first-touch list: the scatter is untracked and
    // the scan below walks `acc[0..n]`. This is exact because (1) `rho.len() >= 2` there, so `filter` holds and an entry is output only when
    // |v| > DROP_TOL: a slot that cancelled to zero (either sign) or was never touched (+0.0) is not output, and a NaN is not either (`!(NaN > tol)`),
    // exactly as the touched list would have it; (2) `rrows` stores structural indices only (`FlatRows::push/set/from_rows` assert j < n), so the
    // scatter writes `acc[0..n]` only and the scan restores all of it; (3) `row_product_work` books `8.5 nr + 2 p5 + 2.5 n` for this route and does
    // not read the distinct count (the returned `distinct` is 0 there). Every slot is restored to +0.0 (an untouched or cancelled slot is already
    // +0.0: `+0.0 + x` is never -0.0 in round-to-nearest, so `acc == ZERO` skips exactly the slots that need no reset).
    if sparse {
        scatter!(step);
    } else {
        assert!(rrows.n == n && rho.len() >= 2);
        scatter!(step_untracked);
    }
    // SAFETY: `cnt` slots were written, and cnt <= n + m (distinct indices).
    unsafe { all_v.set_len(cnt) };
    let all = &all_v[..];
    // Entries that cancelled to (nearly) nothing are dropped (the smallest captured entry is 1.04e-13).
    // A single entry in rho takes the scaled-copy path of `0x134d90` (`0x1360ac`), which has no tolerance test for a scale of 1 or more (fx_e2_r25 call 1: entries of 6e-14 stay);
    // below 1 the same tolerance filter as the general path follows (`0x136629`: `|scale| < 1.0`; seed 50165 call 5: rho -5e-5, a product of 1.25e-17 is dropped).
    let filter = rho.len() != 1 || rho[0].1.abs() < R::ONE;
    let sdir = if above { -1.0 } else { 1.0 };
    // `s * sdir` for s = +-1 is exactly +-sdir.
    let (c_pos, c_neg) = (R::c(sdir), R::c(-sdir));
    let tol = R::c(DROP_TOL);
    // Every entry is written through `out_p`: at most one per touched slot plus one per slack entry of `rho`, which is the reserved capacity
    // (the touched slots are `all.len()` on the sparse route and at most `min(n, est)` on the columnwise one).
    let mut out: Vec<(usize, R)> = Vec::with_capacity(if sparse { all.len() } else { est.min(n) } + rho.len());
    let out_p = out.as_mut_ptr();
    let mut no = 0usize;
    // The count of distinct touched columns: booked only by the sparse multi-row branch of `row_product_work`, 0 on the columnwise route.
    let mut distinct = 0usize;
    // SAFETY: every `j` is below `n + m` (see above), `acc` has `n + m` slots and `st` at least that many (asserted); `no <= all.len()` here.
    if sparse {
        distinct = all.len();
        // One pass per touched slot: the kept entry is appended and the accumulator is restored to +0.0 (the slack entries below do not read it).
        for &j in all {
            let j = j as usize;
            unsafe {
                let v = *acc.get_unchecked(j);
                *acc.get_unchecked_mut(j) = R::ZERO;
                if !filter || v.abs() > tol {
                    out_p.add(no).write((j, v * if *st.get_unchecked(j) == UB { c_neg } else { c_pos }));
                    no += 1;
                }
            }
        }
    } else {
        // Ascending index order over the structurals, four slots at a time: a block of slots equal to zero (all +0.0, see above) is skipped.
        // SAFETY: `acc[0..n]` exists (n <= n + m), `st[j]` for j < n, and `no` counts distinct kept slots below `min(n, est)` (the capacity).
        let mut j = 0usize;
        // Double mode: eight slots at a time are tested as raw words first (an all-zero block has every bit clear; a block holding -0.0 or anything
        // else takes the per-slot path below, which skips the zeros itself), so the scan of a mostly empty `acc` costs about one instruction per slot.
        if !R::QUAD {
            assert!(std::mem::size_of::<R>() == 8 && acc.len() >= n);
            let w = acc.as_ptr() as *const u64;
            while j + 8 <= n {
                // SAFETY: `R` is `f64` (not `QUAD`, size 8), `acc` holds at least `n` slots, and `j + 8 <= n`.
                let z = unsafe {
                    w.add(j).read() | w.add(j + 1).read() | w.add(j + 2).read() | w.add(j + 3).read() | w.add(j + 4).read() | w.add(j + 5).read()
                        | w.add(j + 6).read() | w.add(j + 7).read()
                };
                if z == 0 {
                    j += 8;
                    continue;
                }
                let end = j + 8;
                while j < end {
                    unsafe {
                        let v = *acc.get_unchecked(j);
                        if v != R::ZERO {
                            *acc.get_unchecked_mut(j) = R::ZERO;
                            if v.abs() > tol {
                                out_p.add(no).write((j, v * if *st.get_unchecked(j) == UB { c_neg } else { c_pos }));
                                no += 1;
                            }
                        }
                    }
                    j += 1;
                }
            }
        }
        while j + 4 <= n {
            unsafe {
                let (v0, v1, v2, v3) = (*acc.get_unchecked(j), *acc.get_unchecked(j + 1), *acc.get_unchecked(j + 2), *acc.get_unchecked(j + 3));
                if !((v0 == R::ZERO) & (v1 == R::ZERO) & (v2 == R::ZERO) & (v3 == R::ZERO)) {
                    macro_rules! slot {
                        ($j:expr, $v:expr) => {{
                            let (j, v) = ($j, $v);
                            if v != R::ZERO {
                                *acc.get_unchecked_mut(j) = R::ZERO;
                                if v.abs() > tol {
                                    out_p.add(no).write((j, v * if *st.get_unchecked(j) == UB { c_neg } else { c_pos }));
                                    no += 1;
                                }
                            }
                        }};
                    }
                    slot!(j, v0);
                    slot!(j + 1, v1);
                    slot!(j + 2, v2);
                    slot!(j + 3, v3);
                }
            }
            j += 4;
        }
        while j < n {
            unsafe {
                let v = *acc.get_unchecked(j);
                if v != R::ZERO {
                    *acc.get_unchecked_mut(j) = R::ZERO;
                    if v.abs() > tol {
                        out_p.add(no).write((j, v * if *st.get_unchecked(j) == UB { c_neg } else { c_pos }));
                        no += 1;
                    }
                }
            }
            j += 1;
        }
    }
    let nst = no;
    // Slack entries (their accumulators are not read: a slack touched by the scatter has been restored above).
    // The scatter touches structurals only (`FlatRows` stores indices below `n`), and the rows of `rho` are distinct (a BTRAN result has one entry per row),
    // so a slack occurs at most once in the row: `slack_dup` stays false and the leaving variable's position is that of its single entry.
    let leave_i = leave.wrapping_sub(n);
    let mut lp = usize::MAX;
    // SAFETY: every `i` of `rho` is below `m` (size pass), `slack_nb.len() >= m`, `slack_coef.len() == m`, `st.len() >= n + m`;
    // at most `rho.len()` slack entries follow the `all.len()` structurals.
    for &(i, rv) in rho {
        unsafe {
            let is_leave = i == leave_i;
            if *slack_nb.get_unchecked(i) | is_leave {
                let j = n + i;
                let v = rv * R::c(*slack_coef.get_unchecked(i));
                let s = if *st.get_unchecked(j) == UB { c_neg } else { c_pos };
                // The rows of `rho` are distinct, so `is_leave` holds at most once: a conditional move, no first-occurrence test.
                lp = if is_leave { no } else { lp };
                out_p.add(no).write((j, v * s));
                no += 1;
            }
        }
    }
    let leave_pos = if lp == usize::MAX { None } else { Some(lp) };
    // `slack_dup` is never set: BTRAN's result holds each row once. Checked in debug and shadow builds (shadow parity runs it).
    #[cfg(any(debug_assertions, feature = "factor-shadow"))]
    assert!(
        {
            let mut seen = std::collections::HashSet::new();
            rho.iter().all(|e| seen.insert(e.0))
        },
        "pivot_row_ws: repeated row index in rho"
    );
    // SAFETY: `no` entries were written, no <= capacity.
    unsafe { out.set_len(no) };
    if sparse {
        bits.fill(0);
    }
    ws.acc = acc_v;
    ws.bits = bits_v;
    ws.all = all_v;
    ws.spans = spans_v;
    ws.leave_pos = leave_pos;
    ws.slack_dup = false;
    (out, nst, est, distinct)
}

/// The dual step over the pivot row: `d[j] -= theta * a` for every entry, with the entry negated when `j` is one of the flipped variables.
///
/// # Safety
/// Every index of `row` must be below `d.len()` (`flips` are checked). The rows come from `pivot_row_ws`, which stores only indices below `n + m` (checked when
/// stored), and `d` has `n + m` entries (`Dual::iteration` asserts that before the call).
#[inline]
unsafe fn d_row_step<R: Real>(d: &mut [R], row: &[(usize, R)], flips: &[usize], theta: R, mark: &mut Vec<bool>) {
    if cfg!(test) {
        assert!(row.iter().all(|e| e.0 < d.len()), "d_row_step: index out of range");
    }
    let dp = d.as_mut_ptr();
    if flips.is_empty() {
        for &(j, a) in row {
            // SAFETY: j < d.len(), see above. Repeated indices in a row update sequentially, as with a checked loop.
            unsafe {
                let p = dp.add(j);
                *p = *p - theta * a;
            }
        }
        return;
    }
    if mark.len() < d.len() {
        mark.resize(d.len(), false);
    }
    for &j in flips {
        mark[j] = true;
    }
    let mp = mark.as_ptr();
    for &(j, a) in row {
        // SAFETY: j < d.len() <= mark.len(), see above.
        unsafe {
            let a = if *mp.add(j) { -a } else { a };
            let p = dp.add(j);
            *p = *p - theta * a;
        }
    }
    for &j in flips {
        mark[j] = false;
    }
}

/// What the candidate loops of the ratio test leave besides the lists: the running `bound`, the position and contribution of the best ratio and the
/// largest negative entries of the kept (`xm13`) and rejected (`xm9`) candidates.
struct Cand<R> {
    bound: R,
    best_i: i64,
    bestc: R,
    xm13: R,
    xm9: R,
    /// Some entry has magnitude exactly `ALPHA_TOL` (`ZERO` loops only): it is still raw in the row and the zeroing pass has to run.
    edge: bool,
}

/// The two candidate loops of `ratio_test_with` (entries from `esi` on that sit at a bound, then the earlier ones), filling the retained lists.
///
/// With `ZERO` the loops also write `+0.0` over every entry of magnitude below `ALPHA_TOL` (the zeroing of the ratio test, `0xc7741`, `0xc7870`; the entries that
/// are candidates have `v < -ALPHA_TOL` and are not touched, NaN is not touched); entries of magnitude exactly `ALPHA_TOL` stay raw and set `edge`, so that the
/// reads of the row before the zeroing pass (`framework_rhs` zeroes below the tolerance only) see what they saw before the fusion.
///
/// # Safety
/// `st` and `d` must be longer than every index of `row`.
#[allow(clippy::too_many_arguments)]
#[inline]
unsafe fn ratio_candidates<R: Real, const ZERO: bool>(
    lo: &[f64],
    up: &[f64],
    st: &[i32],
    d: &[R],
    row: &mut [(usize, R)],
    esi: usize,
    delta: R,
    xmm7: R,
    xmm12: R,
    variant: bool,
    ratio_arr: &mut Vec<R>,
    con: &mut Vec<R>,
    lst: &mut Vec<usize>,
) -> Cand<R> {
    let cnt = row.len();
    if cfg!(test) {
        assert!(row.iter().all(|e| e.0 < st.len() && e.0 < d.len()), "ratio_candidates: index out of range");
    }
    let mut bound = R::c(-1e100);
    let (mut best, mut best_i, mut bestc) = (R::c(-1e100), -1i64, R::ZERO);
    let (mut xm13, mut xm9) = (R::ZERO, R::ZERO);
    let (neg_tol, tol10) = (R::c(-ALPHA_TOL), R::c(1e-10));
    let mut edge = false;
    ratio_arr.clear();
    con.clear();
    lst.clear();
    macro_rules! track {
        ($i:expr, $ratio:expr, $c:expr) => {{
            let (ratio, c) = ($ratio, $c);
            let k = lst.len() as i64;
            if ratio > tol10 + best {
                best_i = k;
                best = ratio;
                bestc = c;
            } else if ratio >= best - tol10 {
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
        // SAFETY: j < n + m <= st.len(), d.len() (see the contract).
        if v < neg_tol && unsafe { *st.get_unchecked(j) == LB || *st.get_unchecked(j) == UB } {
            let dj = unsafe { *d.get_unchecked(j) };
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
        } else if ZERO {
            let a = v.abs();
            if a <= tol10 {
                if a < tol10 {
                    row[i].1 = R::ZERO;
                } else {
                    edge = true;
                }
            }
        }
    }
    for i in 0..esi.min(cnt) {
        let (j, v) = row[i];
        if v < neg_tol {
            // SAFETY: as above.
            let dj = unsafe { *d.get_unchecked(j) };
            let ratio = dj / v;
            if bound > ratio {
                xm9 = minr(v, xm9);
                continue;
            }
            xm13 = minr(v, xm13);
            let c;
            if variant && unsafe { *st.get_unchecked(j) } == LB && up[j] >= BIG {
                if dj >= R::ZERO {
                    bound = maxr(PIVOT_DELTA(delta, v, ratio), bound);
                }
                c = xmm12;
            } else {
                let c0 = R::c(lo[j] - up[j]) * v;
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
        } else if ZERO {
            let a = v.abs();
            if a <= tol10 {
                if a < tol10 {
                    row[i].1 = R::ZERO;
                } else {
                    edge = true;
                }
            }
        }
    }
    Cand { bound, best_i, bestc, xm13, xm9, edge }
}

/// The candidate list of the ratio test (ratios, contributions, row positions) and the kept list; retained across calls, cleared at the start of each.
struct RatioBufs<R> {
    ratio_arr: Vec<R>,
    con: Vec<R>,
    lst: Vec<usize>,
    kept: Vec<usize>,
    /// `Cand::edge` of the last call.
    edge: bool,
}

impl<R> Default for RatioBufs<R> {
    fn default() -> Self {
        RatioBufs { ratio_arr: Vec::new(), con: Vec::new(), lst: Vec::new(), kept: Vec::new(), edge: false }
    }
}

/// Primal infeasibility of a basic variable with value `x` and bounds `lo`, `up`.
#[cfg(test)]
#[inline(always)]
fn infeas_of<R: Real>(x: R, lo: f64, up: f64) -> R {
    if lo > -BIG && x < R::c(lo) {
        R::c(lo) - x
    } else if up < BIG && x > R::c(up) {
        x - R::c(up)
    } else {
        R::ZERO
    }
}

/// The bounds of a basic variable as the comparisons of `infeas_of` use them: `R::c(lo)`, or `-inf` when `lo <= -BIG` (the test `x < -inf` is never true, as the
/// skipped branch), and `R::c(up)` or `+inf` likewise.
#[inline(always)]
fn basic_bound_pair<R: Real>(lo: f64, up: f64) -> (R, R) {
    (
        if lo > -BIG { R::c(lo) } else { R::c(f64::NEG_INFINITY) },
        if up < BIG { R::c(up) } else { R::c(f64::INFINITY) },
    )
}

/// `infeas_of` with the bounds in `basic_bound_pair` form: the same comparisons on the same values.
#[inline(always)]
fn infeas_cached<R: Real>(x: R, blo: R, bup: R) -> R {
    if x < blo {
        blo - x
    } else if x > bup {
        x - bup
    } else {
        R::ZERO
    }
}

/// The classification pass of the walk-list rebuild over rows `0..xb.len()`: the infeasible rows (infeasibility above `ptol`) go to `rows` in ascending order with
/// `key(row, infeasibility)` in `keys`, get flag 0, and every other row gets flag 2. Returns the number of infeasible rows and of those whose flag was not 0 before.
#[inline(always)]
fn rebuild_rows<R: Real>(
    xb: &[R],
    blo: &[R],
    bup: &[R],
    ptol: f64,
    flags: &mut [u8],
    rows: &mut Vec<usize>,
    keys: &mut Vec<f64>,
    mut key: impl FnMut(usize, R) -> f64,
) -> (usize, usize) {
    rows.clear();
    keys.clear();
    let m = xb.len();
    let (blo, bup, flags) = (&blo[..m], &bup[..m], &mut flags[..m]);
    let mut n_fresh = 0usize;
    for (r, ((&x, &lo), &up)) in xb.iter().zip(blo).zip(bup).enumerate() {
        let fr = infeas_cached::<R>(x, lo, up);
        if fr.to_f64() > ptol {
            n_fresh += (flags[r] != 0) as usize;
            rows.push(r);
            keys.push(key(r, fr));
        }
    }
    flags.fill(2);
    for &r in rows.iter() {
        flags[r] = 0;
    }
    (rows.len(), n_fresh)
}

/// The per-row cache of the bounds of the basic variables for the basis `head`.
fn basic_bounds<R: Real>(p: &Internal, head: &[usize]) -> (Vec<R>, Vec<R>) {
    head.iter().map(|&v| basic_bound_pair::<R>(p.lo[v], p.up[v])).unzip()
}

/// Retained per-iteration buffers of `Dual::iteration` and `update_primal_and_weights`. `tau` and `rhs` are either empty or of length `m` with every entry
/// `R::ZERO`; the user clears the entries it set and puts the buffer back (a path that returns early just drops it).
struct GlueBufs<R: Real> {
    tau: Vec<R>,
    rhs: Vec<R>,
    changed: Vec<usize>,
    /// Per-iteration lists of `Dual::iteration`, taken (cleared) at the start and put back when the pivot is built; capacities only.
    rho: Vec<(usize, R)>,
    rows_q: Vec<usize>,
    vals_q: Vec<R>,
    rho_rows: Vec<usize>,
    rho_vals: Vec<R>,
    needed: Vec<usize>,
    sq: Vec<R>,
    /// The `col`/`tau` of the last `Pivot` and its box, recovered by `Dual::reclaim` after the post-iteration pass.
    col: Vec<(usize, R)>,
    tau_col: Vec<R>,
    pivot: Option<Box<Pivot<R>>>,
    /// The walk list before its last rebuild (swapped with `Dual::list`), and the sort keys of the rebuild.
    spare_list: Vec<usize>,
    keys: Vec<f64>,
}

impl<R: Real> Default for GlueBufs<R> {
    fn default() -> Self {
        GlueBufs {
            tau: Vec::new(),
            rhs: Vec::new(),
            changed: Vec::new(),
            rho: Vec::new(),
            rows_q: Vec::new(),
            vals_q: Vec::new(),
            rho_rows: Vec::new(),
            rho_vals: Vec::new(),
            needed: Vec::new(),
            sq: Vec::new(),
            col: Vec::new(),
            tau_col: Vec::new(),
            pivot: None,
            spare_list: Vec::new(),
            keys: Vec::new(),
        }
    }
}

/// Retained workspace of `Dual::framework_rhs`: `dense` (all zero) and `seen` (all false) of length `m` between calls, the other lists empty.
struct FwWs<R: Real> {
    dense: Vec<R>,
    seen: Vec<bool>,
    order: Vec<usize>,
    list: Vec<(usize, R)>,
    /// Squares of `list` for its norm.
    sq: Vec<R>,
}

impl<R: Real> Default for FwWs<R> {
    fn default() -> Self {
        FwWs { dense: Vec::new(), seen: Vec::new(), order: Vec::new(), list: Vec::new(), sq: Vec::new() }
    }
}

struct Dual<'a, R: Real = f64> {
    p: &'a Internal,
    n: usize,
    m: usize,
    head: Vec<usize>,
    /// Bounds of the basic variable of each row (`basic_bound_pair`), kept current at every change of `head` (`set_head`, `sync_bounds`).
    blo: Vec<R>,
    bup: Vec<R>,
    /// Status per variable: the basis position when basic, else `LB`, `UB`, `FREE`, `FIXED`.
    st: Vec<i32>,
    /// Basic values by position.
    xb: Vec<R>,
    /// Reduced costs in the reference's frame.
    d: Vec<R>,
    /// Dual steepest-edge weights by position.
    w: Vec<R>,
    mirror: Mirror<R>,
    factor: AuxiliaryFactor,
    /// Row-wise copy of the nonbasic, non-fixed structurals in the reference's storage order, with their immutable first-match coefficients.
    rrows: FlatRows,
    /// The ratio test asked for a refactorization (request 5) instead of ending the dual: the iteration returns without a pivot.
    absorb_request: bool,
    /// The candidates of a row on a fresh factor at a nonzero level cannot absorb the infeasibility (`0xcc1ca`: `ctx+0x290` = 1, no entering variable): the pass ends with status -1004 at the post-iteration pass.
    fresh_absorb: bool,
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
    /// The primal gave up on the artificial-bound handoff of the enclosing pass: this pass runs past the request-4 refactorization. Only the first pass after the give-up (the enclosing controller's later passes end as usual: 5959 dual from iteration 5371, status -1004).
    no_art_exit: bool,
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
    automatic_quad: bool,
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
    /// The double driver's work counters `P+0x50` (the nonbasic row lengths of the nonzero rho entries) and `P+0x58` (the result counts of the pivot columns' two-vector
    /// FTRANs, `0x120238`; m/2 for a dense result), both since the start or the last switch to the mode 3 (`0xf7f22`); the gate of the switch to the mode 2 (`0xf7a54`) compares them.
    p50d: f64,
    p58d: f64,
    /// The mode 2 weight update's counters (`ctx+0x400`, `0x408`, `0x410`, `0x414`, `0xc2f8e`): the framework size, the sum of the accuracy ratios, their count,
    /// and the number of times the size was cut.
    dv_f: f64,
    dv_g: f64,
    dv_k: i32,
    dv_c: i32,
    flags: Vec<u8>,
    list: Vec<usize>,
    thr: R,
    /// Relative residual of the last full basic-solution recomputation.
    xb_residual: R,
    value_mode: bool,
    cursor: usize,
    walk_len: usize,
    /// Scratch of the pricing passes: the infeasibility of every row.
    /// Slacks the row product offers, by row: those nonbasic at the last factorization and those that have left the basis since. A slack that enters the
    /// basis stays on offer until the next factorization.
    slack_nb: Vec<bool>,
    /// The retained workspace of the row product.
    prow_ws: std::cell::RefCell<PivotWs<R>>,
    /// Retained buffers of the per-iteration glue: dense scratch vectors of length `m` (all zero between uses, cleared through the touched lists),
    /// the changed-row list, and the workspace of `framework_rhs`.
    glue: GlueBufs<R>,
    fw_ws: std::cell::RefCell<FwWs<R>>,
    /// Retained buffers of the ratio test and the flip marks of the d update (all empty between calls).
    rbuf: RatioBufs<R>,
    flip_mark: Vec<bool>,
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
    last_y: Vec<R>,
    /// Retained scratch of the refactor-time passes (no per-refactor allocation): the weights by variable (all zero between uses), the nonzero rows and values of
    /// the x_B right-hand side, the double copy of `last_y`, the nonbasic structural values of the last x_B recompute.
    scr_wv: Vec<R>,
    scr_rows: Vec<usize>,
    scr_vals: Vec<R>,
    scr_ydbl: Vec<f64>,
    scr_nbv: Vec<f64>,
    /// The variables at their upper bound found by `assign_d_single_pass`'s initialisation (negated at its end), and the right-hand side of `btran_y`.
    scr_ub: Vec<u32>,
    scr_binp: Vec<(usize, R)>,
    /// `rrows` is the set of the nonbasic, non-fixed structurals as the pivots keep it (false after a handoff list or an in-place repair, which change the
    /// set without the row copy).
    rrows_exact: bool,
    /// A structural column repeats a row (MPS input is not merged): the reset then keeps the deduplicating list path.
    dup_entries: bool,
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
    pending_alpha: Option<Vec<(usize, R)>>,
    /// The pivot-check pair of the last iteration that reached the check (`ctx+0x300` = the column's pivot, `ctx+0x308` = the row's): the primal started by the
    /// controller reads it at its first pricing call.
    pair: (f64, f64),
    /// Final status once the loop ends.
    status: Option<Status>,
    // ctx+0xf0 = -1002, pending until the recompute finishes: the
    // residual-scale acceptance gate may replace it with optimal status.
    drift_status: bool,
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
    /// Quad working objective, absent after a writer whose arithmetic is not yet modelled.
    working_obj: Option<f64>,
    /// The double dual's maintained objective (`ctx+0x118`): reset at every x_B recompute, advanced by each accepted pivot, read by the stall check.
    obj: f64,
    /// The model's objective constant (`ctx+0x128`, raw, not divided by the objective scale; -0.0 when the model has none).
    obj_const: f64,
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
    /// Set by `perturb` inside the post-factorization pass, which the reference books into the same `0xe7000` total.
    perturbed: bool,
    pert_param: f64,
    obj_scale: f64,
    /// The model was scaled (`ctx+0xe0` non-null): the growth check of the factorization uses 1e12.
    scaled: bool,
    pcost: Vec<f64>,
    /// The forced-row stack of the public `DStart` start (`ctx+0x530`, count `ctx+0x52c`): the basic variables with |d| > 1e-7 at the start, the last entry on top.
    /// An empty stack is the reference's count -1 or 0.
    forced: Vec<usize>,
    /// The stack was filled (count > 0 at the start): `ctx+0x52c` is never 0 again, so the leaving variable's d is not clamped (`0xc5d85`).
    forced_ever: bool,
    /// The last pricing call popped a stale entry (a nonbasic variable or a basic one with d = 0): the iteration ends without a row.
    forced_stale: bool,
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

fn reference_order<'p, R: Real>(
    p: &'p Internal,
    basic: &[usize],
    fac_no: usize,
) -> (Vec<usize>, crate::lu_order::Structure<R>, crate::lu_order::BasisCols<'p>) {
    let cols = crate::lu_order::BasisCols::new(p.split_cols(), basic);
    crate::lu_order::set_factorization_number(fac_no);
    let st = crate::lu_order::factor_structure_r::<R>(p.m, &cols);
    let order = st.pivots.iter().map(|&(_, pos)| basic[pos]).collect();
    (order, st, cols)
}

thread_local! {
    static QUAD_RESULT: std::cell::Cell<bool> = const { std::cell::Cell::new(false) };
    /// The scaled long-double solution of a quad final pass as (high, low) double pairs, for the reference's objective read-out.
    static QUAD_XR: std::cell::RefCell<Option<Vec<(f64, f64)>>> = const { std::cell::RefCell::new(None) };
    static DUAL_XR_AT_HANDOFF: std::cell::RefCell<Option<Vec<(f64, f64)>>> = const { std::cell::RefCell::new(None) };
    /// The last dual pass ended with status 12 from its own pivot check (`0x21db9a`, double copy `0xc5a..`): the controller writes the status and the iteration count
    /// at `0x15637b` and skips the objective store of `0x153eba`, so the result keeps objective 0 (witness_primal_r31_double4; the `0x15553e` end stores `ctx+0x118`).
    static PASS_END12: std::cell::Cell<bool> = const { std::cell::Cell::new(false) };
    /// The double dual's x (scaled frame) at its last -1004 handoff to the primal: `ctx+0x118` keeps the objective of that x, and the controller's status-12 end
    /// (the primal's own, `0x144e6f`) stores it as ObjVal (`fx_midsing_1031`: 197979259134.57623, not the 0.0 of the `0x15637b` path).
    static DUAL_X_AT_HANDOFF: std::cell::RefCell<Option<Vec<f64>>> = const { std::cell::RefCell::new(None) };
    /// The same x, kept for the unscaled read-out of a status-12 result.
    static END12_X: std::cell::RefCell<Option<Vec<f64>>> = const { std::cell::RefCell::new(None) };
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

/// The read-out state of the last result (long-double solution and mark), so that a pass whose result is discarded can leave it as it was.
pub(crate) fn quad_result_state() -> (bool, Option<Vec<(f64, f64)>>) {
    (
        QUAD_RESULT.with(|c| c.get()),
        QUAD_XR.with(|c| c.borrow().clone()),
    )
}

pub(crate) fn set_quad_result_state(st: (bool, Option<Vec<(f64, f64)>>)) {
    QUAD_RESULT.with(|c| c.set(st.0));
    QUAD_XR.with(|c| *c.borrow_mut() = st.1);
}

/// The objective of the unscaled solution of the dual's last result: `long double` accumulation when its final pass was the quad pass. Clears the mark.
pub fn unscaled_objective(model: &crate::model::LpModel, x: &[f64], col: &[f64]) -> f64 {
    if PASS_END12.with(|c| c.get()) {
        QUAD_RESULT.with(|c| c.set(false));
        return 0.0;
    }
    if let Some(xs) = END12_X.with(|c| c.borrow_mut().take()) {
        QUAD_RESULT.with(|c| c.set(false));
        if xs.len() >= model.obj.len() {
            let mut acc = 0.0f64;
            for (j, c) in model.obj.iter().enumerate() {
                acc += *c * (xs[j] * col[j]);
            }
            return acc + model.offset;
        }
    }
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

/// The value of a nonbasic variable of status `s` with bounds `lo`, `up`: zero at an artificial bound when `zero_art` (the first kind of `at_artificial`
/// for the upper bound, the second for the lower), else the bound it sits at, 0 for a free one.
#[inline(always)]
fn nb_value(zero_art: bool, s: i32, lo: f64, up: f64) -> f64 {
    match s {
        LB => {
            if zero_art && lo <= -BIG {
                0.0
            } else {
                lo
            }
        }
        FIXED => lo,
        UB => {
            if zero_art && up >= BIG {
                0.0
            } else {
                up
            }
        }
        _ => 0.0,
    }
}

impl<'a, R: Real> Dual<'a, R> {
    fn nonbasic_value(&self, j: usize) -> f64 {
        nb_value(self.zero_artificial, self.st[j], self.p.lo[j], self.p.up[j])
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

    #[inline(always)]
    fn infeas(&self, r: usize) -> f64 {
        self.infeas_real(r).to_f64()
    }

    #[inline(always)]
    fn infeas_real(&self, r: usize) -> R {
        infeas_cached::<R>(self.xb[r], self.blo[r], self.bup[r])
    }

    /// `head[r] = j` with the bound cache of the row.
    #[inline(always)]
    fn set_head(&mut self, r: usize, j: usize) {
        self.head[r] = j;
        let (lo, up) = basic_bound_pair::<R>(self.p.lo[j], self.p.up[j]);
        self.blo[r] = lo;
        self.bup[r] = up;
    }

    /// Rebuild of the bound cache after a wholesale change of `head`.
    fn sync_bounds(&mut self) {
        (self.blo, self.bup) = basic_bounds::<R>(self.p, &self.head);
    }

    /// The cache equals the tables (debug and shadow builds).
    #[cfg(any(debug_assertions, feature = "factor-shadow"))]
    fn check_bounds(&self) {
        assert_eq!(self.blo.len(), self.head.len());
        assert_eq!(self.bup.len(), self.head.len());
        for (r, &v) in self.head.iter().enumerate() {
            let (lo, up) = basic_bound_pair::<R>(self.p.lo[v], self.p.up[v]);
            assert!(self.blo[r] == lo && self.bup[r] == up, "stale bound cache at row {r}");
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
        crate::lu_order::BasisCols<'a>,
        Option<Mirror<R>>,
    ) {
        crate::lu_order::set_threshold(self.lu_thr);
        let (order, st, rows) = reference_order::<R>(self.p, &self.head, self.refactorizations + 1);
        let mi = (order.len() == self.m).then(|| Mirror::fresh(self.m, self.n, &st, rows.view()));
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
                self.set_head(pos, n + r);
                last = Some(out);
            }
            if let Some(j) = last {
                self.st[j] = FIXED;
                self.d[j] = R::c(if j < n { self.cost[j] } else { 0.0 });
                self.repaired = Some(j);
            }
            // The repair retains the pending request: in particular, a request-4
            // confirmation still runs its residual-scale feasibility test.
            self.rrows_exact = false;
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
                .diag()
                .iter()
                .fold((0.0f64, f64::INFINITY), |(a, b), v| {
                    (a.max(v.abs().to_f64()), b.min(v.abs().to_f64()))
                });
            if !(mx / mn > if self.scaled { 1e12 } else { 1e15 }) {
                break;
            }
            self.factor_work(&st, rows.view(), discard);
            self.lu_thr = if self.lu_thr < 0.5 { 0.5 } else { 0.9 };
            if let Some(retired) = mi.take() {
                retired.recycle();
            }
            (order, st, rows, mi) = self.lu_attempt();
            if order.len() < self.m {
                self.singular_payload();
                return Err(SolveError::SingularBasis);
            }
        }
        let mi = mi.unwrap();
        // The weights by variable on a retained buffer (all zero between uses).
        let mut wv = std::mem::take(&mut self.scr_wv);
        wv.resize(self.n + self.m, R::ZERO);
        for (&v, &w) in self.head.iter().zip(&self.w[..self.m]) {
            wv[v] = w;
        }
        // `order` holds each basic variable once (one pivot per position), so zeroing `wv[v]` as it is read leaves the buffer all zero, as the separate
        // pass over `head` did.
        for (pos, (&v, w)) in order.iter().zip(&mut self.w[..self.m]).enumerate() {
            self.st[v] = pos as i32;
            *w = std::mem::replace(&mut wv[v], R::ZERO);
        }
        self.scr_wv = wv;
        self.head = order;
        self.sync_bounds();
        self.reset_pricing(&st, rows.view(), &mi);
        std::mem::replace(&mut self.mirror, mi).recycle();
        self.mirror.set_upd_limit(self.upd_limit);
        self.factor_work(&st, rows.view(), discard);
        self.factor = AuxiliaryFactor::lazy(self.head.clone());
        // Shadow build: keep the eager build and its error return; the lazy build must agree on success and failure.
        #[cfg(feature = "factor-shadow")]
        match AuxiliaryFactor::build(self.m, &self.head, self.p, 0.0) {
            Ok(eager) => {
                self.factor.base = Some(eager);
                self.factor.verify = true;
            }
            Err(_) => {
                assert!(
                    AuxiliaryFactor::build(self.m, &self.head, self.p, 0.0).is_err(),
                    "shadow: lazy auxiliary build succeeded where the eager one failed"
                );
                return Err(SolveError::SingularBasis);
            }
        }
        self.refactorizations += 1;
        self.updates = 0;
        self.bad_rows.clear();
        self.old_rej = 0;
        self.worst = 0.0;
        self.worst_row = -2;
        let mut slack_nb = std::mem::take(&mut self.slack_nb);
        slack_nb.clear();
        slack_nb.extend(self.st[self.n..self.n + self.m].iter().map(|&t| t < 0));
        self.slack_nb = slack_nb;
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
            self.factor.btran(self.p, &mut y)?;
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
        self.set_head(r, j);
        self.iters += 1;
        crate::simplex::TOTAL_ITERS.with(|c| c.set(c.get() + 1.0));
        if self.mirror.refactor_due() {
            self.refactor("free_pivot")?;
        } else {
            self.factor
                .update(r, alpha.iter().map(|v| v.to_f64()).collect());
            self.mirror.update(r);
            self.mirror.clear_update_work();
            self.updates += 1;
        }
        Ok(())
    }

    /// The work of a factorization and the installation of its factor (`0xd4190` from `0xf34f2`, the update routine `0xf15c0`), after the fresh mirror is built:
    /// the driver's sites, the spike discard and the weight update of the pivot that preceded it, then the factor object's own passes.
    fn factor_work(
        &mut self,
        st: &crate::lu_order::Structure<R>,
        rows: crate::lu_order::Cols<'_>,
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
        let mut rows: Vec<Vec<(usize, f64)>> = vec![Vec::new(); self.m];
        for j in 0..self.n {
            if self.st[j] < 0 && self.st[j] != FIXED && self.p.lo[j] != self.p.up[j] {
                for &(i, _) in &self.p.cols[j] {
                    let a = self.p.cols[j].iter().find(|e| e.0 == i).unwrap().1;
                    rows[i].push((j, a));
                }
            }
        }
        self.rrows = FlatRows::from_rows(&rows, self.n);
        // `rrows` takes the first entry of a repeated (row, column) pair for every copy, `p.rwise` keeps each entry's own
        // value: with a repeated row in any structural column the single-pass d could differ, so keep the full path.
        let mut seen = vec![usize::MAX; self.m];
        'cols: for j in 0..self.n {
            for &(i, _) in &self.p.cols[j] {
                if seen[i] == j {
                    self.rrows_exact = false;
                    self.dup_entries = true;
                    break 'cols;
                }
                seen[i] = j;
            }
        }
    }

    // ------------------------------------------------------------------ basic solution

    /// x_B from scratch (`0xc0c29` path of `0xc0050`): FTRAN of `b - N x_N`.
    ///
    /// `req` is the refactorization request in force (`ctx+0x294`) when the recompute follows a factorization, 0 otherwise. The accuracy of the solve then
    /// raises the LU pivot threshold for the next factorizations (`0xc1861`, `accuracy_threshold`).
    fn recompute_xb(&mut self, req: u32) -> Result<(), SolveError> {
        let mut tol_abort = false;
        let p = self.p;
        let mut r: Vec<R> = p.rhs.iter().map(|&v| R::c(v)).collect();
        // Entries of the structural nonbasic columns at a nonzero value (the read-out below, `0xc142b`), counted in the same pass.
        let mut cols = 0usize;
        // The nonbasic structural values, computed once here and read again by the objective sum below.
        let mut nbv = std::mem::take(&mut self.scr_nbv);
        nbv.clear();
        nbv.resize(p.n, 0.0);
        let zero_art = self.zero_artificial;
        let nv = p.n + p.m;
        for (j, ((&s, &lo), &up)) in self.st[..nv].iter().zip(&p.lo[..nv]).zip(&p.up[..nv]).enumerate() {
            if s < 0 {
                let v = nb_value(zero_art, s, lo, up);
                if j < p.n {
                    nbv[j] = v;
                }
                if v != 0.0 {
                    if j < p.n {
                        cols += p.cols[j].len();
                    }
                    for &(i, a) in &p.cols[j] {
                        r[i] -= R::c(a) * R::c(v);
                    }
                }
            }
        }
        let mut rows = std::mem::take(&mut self.scr_rows);
        let mut vals = std::mem::take(&mut self.scr_vals);
        rows.clear();
        vals.clear();
        for i in 0..p.m {
            if r[i] != R::ZERO {
                rows.push(i);
                vals.push(r[i]);
            }
        }
        // The estimate of the dense FTRAN is booked before the read-out term. On the dense route it is the work of the real solve below
        // (booked right after it, the solve books nothing here); on the others it is the replay, booked first.
        let from_solve = self.mirror.dense_route(rows.len());
        #[cfg(any(test, debug_assertions, feature = "factor-shadow"))]
        let replay = self.mirror.ftran_dense_work(&r);
        if !from_solve {
            let est = self.mirror.ftran_dense_work(&r);
            self.work
                .add_from(0x1219d6, 0x120673, est, crate::work::FAC);
            // `0xc142b`: the solution read-out sums the columns of the nonbasic variables at a nonzero value.
            self.work.add(
                0xc142b,
                4.0 * p.m as f64 + 2.0 * cols as f64,
                crate::work::CTX,
            );
        }
        let out = self.mirror.ftran_basic(&rows, &vals);
        if from_solve {
            let est = self.mirror.dense_work_of_solve(self.mirror.work(), rows.len());
            #[cfg(any(test, debug_assertions, feature = "factor-shadow"))]
            assert_eq!(est.to_bits(), replay.to_bits(), "recompute_xb: dense work of the solve differs from the replay");
            self.work
                .add_from(0x1219d6, 0x120673, est, crate::work::FAC);
            self.work.add(
                0xc142b,
                4.0 * p.m as f64 + 2.0 * cols as f64,
                crate::work::CTX,
            );
        }
        self.scr_rows = rows;
        self.scr_vals = vals;
        let mut xb = std::mem::take(&mut self.xb);
        xb.clear();
        xb.resize(p.m, R::ZERO);
        for (i, v) in out {
            xb[i] = v;
        }
        // `0xc10ab..0xc1129`, `0xc1480..0xc1517`: R = max(1e-9, max|rhs|), the largest |x_B| and the residual of B x_B against the rhs, both relative to R;
        // the residual is read only when the factorization count is 4 mod 5 or the request is 4, 6 or 8 (`0xc1137..0xc11ce`, `0xc165b..0xc1699`).
        let rmax = r.iter().fold(R::c(1e-9), |a, v| a.max(v.abs()));
        let xmax = xb.iter().fold(R::ZERO, |a, v| a.max(v.abs()));
        let e2 = if self.nref() % 5 == 4 || matches!(req, 4 | 6 | 8) {
            // The reference forms B*x_B from zero first, then subtracts it from the
            // right-hand side (`0x115d60`, `0xc1480`); subtracting each term from
            // the right-hand side directly can change the residual and its gate.
            let mut product = vec![R::ZERO; p.m];
            for (pos, &j) in self.head.iter().enumerate() {
                for &(i, a) in &p.cols[j] {
                    product[i] += R::c(a) * xb[pos];
                }
            }
            let residual = r
                .iter()
                .zip(&product)
                .fold(R::ZERO, |a, (b, v)| a.max((*b - *v).abs()));
            if self.level == 0 && residual > R::c(primal_tol()) {
                #[cfg(test)]
                CEILING_TEST.with(|c| {
                    if !R::QUAD
                        && req == 4
                        && self.iters == 16
                        && c.borrow().as_ref().is_some_and(|v| v.is_empty())
                    {
                        PRIMAL_DRIFT_COUNT.with(|n| n.set(2));
                        set_primal_tol(1e-3);
                    }
                });
                self.work.add(
                    0xc177e,
                    4.0 * p.m as f64 + 2.0 * cols as f64,
                    crate::work::CTX,
                );
                {
                    let violation = self
                        .head
                        .iter()
                        .enumerate()
                        .fold(R::ZERO, |max, (pos, &j)| {
                            let x = xb[pos];
                            let lo = R::c(p.lo[j]);
                            let up = R::c(p.up[j]);
                            if p.lo[j] > -BIG && lo > x {
                                max.max(lo - x)
                            } else if p.up[j] < BIG && x > up {
                                max.max(x - up)
                            } else {
                                max
                            }
                        });
                    if primal_drift_excess(residual, violation) {
                        PRIMAL_DRIFT_COUNT.with(|c| c.set(c.get() + 1));
                    }
                    if crate::dump::on("pdrift") {
                        crate::dump::emit(
                            "pdrift",
                            serde_json::json!({"kind":"dual", "total": crate::simplex::TOTAL_ITERS.with(|c| c.get()), "nref":self.nref(), "req":req, "residual":residual.to_f64(), "violation":violation.to_f64(), "tol":primal_tol(), "cnt":PRIMAL_DRIFT_COUNT.with(|c| c.get())}),
                        );
                    }
                    let (tol, clear, abort) = primal_drift_action(
                        PRIMAL_DRIFT_COUNT.with(|c| c.get()),
                        primal_tol(),
                        R::QUAD,
                        self.automatic_quad,
                    );
                    if clear {
                        set_primal_tol(tol);
                        PRIMAL_DRIFT_COUNT.with(|c| c.set(0));
                    }
                    tol_abort = abort;
                    #[cfg(test)]
                    if abort {
                        CEILING_TEST.with(|c| {
                            if let Some(v) = c.borrow_mut().as_mut() {
                                v.push(("pending", self.iters, req, self.head.clone()));
                            }
                        });
                    }
                }
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
        if R::QUAD {
            let mut sum = R::ZERO;
            for j in 0..p.n {
                let x = if self.st[j] >= 0 {
                    self.xb[self.st[j] as usize]
                } else {
                    R::c(nbv[j])
                };
                sum += x * (R::c(self.cost[j]) + R::c(self.cost_lo[j]));
            }
            self.working_obj = Some(((sum + R::c(self.obj_off)) / R::c(self.obj_scale)).to_f64());
            #[cfg(test)]
            record_quad_obj(self.working_obj.unwrap());
            if crate::dump::on("qobj") {
                crate::dump::emit(
                    "qobj",
                    serde_json::json!({"kind":"recompute", "it":self.iters, "value":self.working_obj}),
                );
            }
        }
        if !R::QUAD {
            // `0xbfc80..0xbfd96`: the maintained objective restarts from the structurals' sum (index order) over the x_B just computed.
            let mut sum = 0.0f64;
            for j in 0..p.n {
                let x = if self.st[j] >= 0 {
                    self.xb[self.st[j] as usize].to_f64()
                } else {
                    nbv[j]
                };
                sum += x * self.cost[j];
            }
            self.obj = ((sum + self.obj_off) / self.obj_scale) + self.obj_const;
            if crate::dump::on("stallobj") {
                crate::dump::emit(
                    "stallobj",
                    serde_json::json!({"ev":"recompute", "it":self.iters, "nref":self.nref(), "sum":format!("{:016x}", sum.to_bits()),
                        "obj":format!("{:016x}", self.obj.to_bits())}),
                );
            }
        }
        self.scr_nbv = nbv;
        self.drift_status |= tol_abort;
        Ok(())
    }

    fn finish_drift_status(&mut self) -> Result<(), SolveError> {
        if std::mem::take(&mut self.drift_status) {
            #[cfg(test)]
            CEILING_TEST.with(|c| {
                if let Some(v) = c.borrow_mut().as_mut() {
                    v.push(("restart", self.iters, self.request, self.head.clone()));
                }
            });
            // A tolerance abort writes status -1002 without a failed LU.
            SINGULAR_TOL_ABORT.with(|c| c.set(true));
            SINGULAR_BASIS.with(|c| {
                *c.borrow_mut() = Some(Basis {
                    basic: self.head.clone(),
                    at_upper: self.st.iter().map(|&t| t == UB).collect(),
                })
            });
            SINGULAR.with(|c| {
                c.set(Some((
                    self.iters,
                    self.refactorizations,
                    self.lu_thr,
                    self.avg_touch,
                    self.last_cap,
                    self.inf_count,
                    self.last_pivoted,
                    self.retries,
                    self.it_base + self.iters,
                )))
            });
            return Err(SolveError::SingularBasis);
        }
        Ok(())
    }

    /// d recomputed from the costs after a factorization, in the reference frame.
    fn recompute_d(&mut self) {
        // The d of the forced-row stack entries survive the recompute (`0xe53a7..0xe5400` save, `0xe57a8..` restore).
        let saved: Vec<(usize, R)> = self.forced.iter().map(|&v| (v, self.d[v])).collect();
        self.btran_y();
        self.dual_residual_check();
        self.work
            .add_from(0x12de50, 0xe53a7, self.mirror.work(), crate::work::FAC);
        if self.rrows_exact {
            self.assign_d_single_pass();
            #[cfg(any(debug_assertions, feature = "factor-shadow"))]
            {
                let want = self.d_from_full_product(&saved);
                for (v, dv) in saved {
                    self.d[v] = dv;
                }
                for j in 0..self.n + self.m {
                    assert!(
                        want[j].json().to_string() == self.d[j].json().to_string(),
                        "single-pass d differs at {j}: {:?} against {:?}",
                        self.d[j],
                        want[j]
                    );
                }
            }
            #[cfg(not(any(debug_assertions, feature = "factor-shadow")))]
            for (v, dv) in saved {
                self.d[v] = dv;
            }
            return;
        }
        let dt = true_d_rowwise(self.p, &self.st, &self.cost, &self.cost_lo, &self.last_y);
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
        for (v, dv) in saved {
            self.d[v] = dv;
        }
    }

    /// The final d of `recompute_d` in one pass over the nonbasic non-fixed structurals' row copy `rrows` (the reference's list `ctx+0x1d0`): every d_j starts
    /// from its cost (nonbasic) or zero (basic), receives `- y_i a_ij` row by row in ascending row order as `true_d_rowwise` does, and a variable at its upper
    /// bound ends as `-(c - y a)`. The entries the full product visits and `recompute_d` then discards (the basic and fixed structurals) are not visited: the
    /// fixed ones end at their cost, the basic ones at zero. The slack of row `i` is its single entry, subtracted in its row.
    fn assign_d_single_pass(&mut self) {
        let (n, m) = (self.n, self.m);
        let nv = n + m;
        // The initialisation also lists the variables at their upper bound, so the closing negation visits only those (each slot is negated once, after
        // its last subtraction, as the full sweep did).
        let mut ubs = std::mem::take(&mut self.scr_ub);
        ubs.clear();
        let st = &self.st[..nv];
        let d = &mut self.d[..nv];
        for (j, (((dj, &s), &c), &cl)) in d.iter_mut().zip(st).zip(&self.cost[..nv]).zip(&self.cost_lo[..nv]).enumerate() {
            *dj = if s < 0 {
                if s == UB {
                    ubs.push(j as u32);
                }
                R::c(c) + R::c(cl)
            } else {
                R::ZERO
            };
        }
        let drop = R::c(DROP_TOL);
        for (i, &yi) in self.last_y[..m].iter().enumerate() {
            if yi.abs() > drop {
                let (idx, val) = self.rrows.row(i);
                for (&j, &a) in idx.iter().zip(val) {
                    let j = j as usize;
                    let s = st[j];
                    if s < 0 && s != FIXED {
                        d[j] -= yi * R::c(a);
                    }
                }
                if st[n + i] < 0 {
                    let a = self.p.cols[n + i][0].1;
                    d[n + i] -= yi * R::c(a);
                }
            }
        }
        for &j in &ubs {
            let dj = &mut d[j as usize];
            *dj = -*dj;
        }
        self.scr_ub = ubs;
    }

    /// The full product and the assignment loop of `recompute_d` as they were before the single pass (reference for the debug and shadow check); the d of the forced rows
    /// is not restored.
    #[cfg(any(debug_assertions, feature = "factor-shadow"))]
    fn d_from_full_product(&self, saved: &[(usize, R)]) -> Vec<R> {
        let dt = true_d_rowwise(self.p, &self.st, &self.cost, &self.cost_lo, &self.last_y);
        let mut d = vec![R::ZERO; self.n + self.m];
        for j in 0..self.n + self.m {
            d[j] = if self.st[j] == FIXED {
                if j < self.n {
                    R::c(self.cost[j]) + R::c(self.cost_lo[j])
                } else {
                    dt[j]
                }
            } else if self.st[j] < 0 {
                if self.st[j] == UB { -dt[j] } else { dt[j] }
            } else {
                R::ZERO
            };
        }
        for &(v, dv) in saved {
            d[v] = dv;
        }
        d
    }

    /// The log timer of the wrapper (`0xf9306..0xf6e13`): on a fresh factor with no request pending, the first time the total iteration count is 100 past the last
    /// snapshot, the basis is saved (`0xd40b0`, only while the restart counter is below 5). The first one (at the start) is the basis a restart from scratch builds.
    fn log_snapshot(&mut self) {
        if self.updates != 0 || self.request != 0 || self.status.is_some() {
            return;
        }
        // `ctx+0x2a8`: every iteration and bound flip of the controller, which `it0 + iters` misses after a restore (5959: the dual from iteration 130 reaches 401).
        let total = if R::QUAD {
            crate::simplex::TOTAL_ITERS.with(|c| c.get())
        } else {
            (self.it0 + self.iters) as f64
        };
        let c40 = SNAP_C40.with(|c| c.get());
        if total - c40 < 100.0 {
            return;
        }
        SNAP_C40.with(|c| c.set(total));
        // The timer advances at every restart counter; only the save is skipped from five (standalone_restart_context 5959: `ctx+0x40` = 401 after the dual that ended at iteration 147 with `ctx+0x184` = 9).
        if self.retries >= 5 {
            return;
        }
        let quality = if dual_infeas_sum() > 0.0 { 1 } else { 2 };
        // The dual snapshot (`0xd40b0`) skips a dropped variable and keeps a
        // previous snapshot of higher quality.
        // GUROBI: `0xd40b0` also skips the save while the forced-row stack (`ctx+0x52c`) holds rows; not modelled, no witness (none in the release profile).
        if self.snapshot_excluded
            || self.repaired.is_some()
            || SNAP_QUALITY.with(|c| quality < c.get())
        {
            return;
        }
        SNAP_QUALITY.with(|c| c.set(quality));
        SNAP_PRIMAL.with(|c| c.set(false));
        let basic = self.head.clone();
        let at_upper = self.st.iter().map(|&t| t == UB).collect();
        if c40 == -100.0 {
            INITIAL_SNAP.with(|c| *c.borrow_mut() = Some((basic, at_upper)));
        } else {
            SNAP.with(|c| *c.borrow_mut() = Some((basic, at_upper)));
        }
    }

    /// Full dual recompute (`0xbf5b0`, quad `0x217600`). The basic primal
    /// infeasibility sum updates `ctx+0x120`; `ctx+0x2d8` separately sums
    /// positive reduced costs at artificial upper bounds.
    fn snapshot_recomputed(&self) {
        if R::QUAD && !self.value_mode {
            return;
        }
        let mut primal_sum = 0.0;
        let ptol = R::c(primal_tol());
        // The second sum (artificial upper bounds with a positive d, level != 0) shares the pass: the two conditions are disjoint (basic against nonbasic), the
        // accumulators independent, each in ascending j order.
        let tol = R::c(dual_tol());
        let mut sum = R::ZERO;
        let lvl = self.level != 0;
        let nv = self.n + self.m;
        for (j, &s) in self.st[..nv].iter().enumerate() {
            if s < 0 {
                if lvl && s == UB && self.p.up[j] >= BIG && self.d[j] > tol {
                    sum = sum + self.d[j];
                }
                continue;
            }
            let x = self.xb[s as usize];
            let lo = R::c(self.p.lo[j]);
            let up = R::c(self.p.up[j]);
            if lo - x > ptol {
                primal_sum = if R::QUAD {
                    (R::c(primal_sum) + lo - x).to_f64()
                } else {
                    (R::c(primal_sum) + (lo - x)).to_f64()
                };
            } else if self.p.up[j] < BIG && x - up > ptol {
                primal_sum = if R::QUAD {
                    (R::c(primal_sum) - up + x).to_f64()
                } else {
                    (R::c(primal_sum) + (x - up)).to_f64()
                };
            }
        }
        crate::simplex::STALL_PREV.with(|c| c.set(primal_sum));
        #[cfg(test)]
        PRIMAL_SUM_TEST.with(|c| {
            if let Some(records) = c.borrow_mut().as_mut() {
                records.push((R::QUAD, self.iters, self.it0, primal_sum));
            }
        });
        set_dual_infeas_sum(sum.to_f64());
    }

    /// The dual residual check of the post-refactor pass (`0xe5965..0xe5d4c`): `max |c_j - y a_j|` over the basic variables, read when the factorization count is 4 mod 5
    /// or the request is 4 or 6. A residual above the tolerance counts a failure; at the sixth the tolerance grows tenfold (up to 1e-3) and the counter restarts.
    fn dual_residual_check(&mut self) {
        let req = self.request;
        if !(self.nref() % 5 == 4 || req == 4 || req == 6) {
            return;
        }
        let mut meas = R::ZERO;
        let mut y_double = std::mem::take(&mut self.scr_ydbl);
        y_double.clear();
        if !R::QUAD {
            y_double.extend(self.last_y.iter().map(|v| v.to_f64()));
        }
        for &j in &self.head {
            // The kernel of `0xe59a1..0xe5bd9`: the cost (with the shift), minus the column products, sequentially below 8 entries past the alignment peel, else the eight-lane
            // body of the read-out kernel (fx_midsing_1127: 165.177734375 against 142.407).
            let acc = if R::QUAD {
                // The x87 residual check (`0x2387b6..0x238a4d`) accumulates eight
                // long-double lanes, folds them left to right, then walks the tail backwards.
                let col = &self.p.fcols[j];
                let n8 = col.len() / 8 * 8;
                let mut lanes = [R::ZERO; 8];
                lanes[0] = R::c(self.cost[j]) + R::c(self.cost_lo[j]);
                for i in 0..n8 {
                    let (row, a) = col[i];
                    lanes[i % 8] -= R::c(a) * self.last_y[row];
                }
                let mut sum = lanes[0];
                for lane in &lanes[1..] {
                    sum += *lane;
                }
                for &(row, a) in col[n8..].iter().rev() {
                    sum -= R::c(a) * self.last_y[row];
                }
                sum
            } else {
                R::c(crate::simplex::dot_lanes_from(
                    &self.p.fcols[j],
                    &y_double,
                    self.p.offs[j] & 1,
                    self.cost[j],
                ))
            };
            meas = maxr(meas, acc.abs());
        }
        self.scr_ydbl = y_double;
        let tol = dual_tol();
        if crate::dump::on("dualtol") {
            crate::dump::emit(
                "dualtol",
                serde_json::json!({"meas": meas.to_f64(), "tol": tol, "cnt": RESID_FAILS.with(|c| c.get()), "req": req, "n": self.nref()}),
            );
        }
        let mut cnt = RESID_FAILS.with(|c| c.get());
        if meas > R::c(tol) {
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
        self.btran_y();
        true_d_rowwise(self.p, &self.st, &self.cost, &self.cost_lo, &self.last_y)
    }

    /// `y = B^-T c_B` of the reference BTRAN into `last_y` (the right-hand side takes the forced-row stack into account).
    fn btran_y(&mut self) {
        let p = self.p;
        let mut inp = std::mem::take(&mut self.scr_binp);
        inp.clear();
        inp.reserve(p.m);
        for (pos, &j) in self.head.iter().enumerate() {
            let v = R::c(self.cost[j]) + R::c(self.cost_lo[j]);
            if v != R::ZERO {
                inp.push((pos, v));
            }
        }
        if !self.forced.is_empty() {
            // The forced-row stack (`0xe52c3..0xe537f`): the basic stack variables keep their d, so the right-hand side is c_B minus their d.
            let mut rhs = vec![R::ZERO; p.m];
            for &(pos, v) in &inp {
                rhs[pos] = v;
            }
            for &v in &self.forced {
                if self.st[v] >= 0 {
                    rhs[self.st[v] as usize] = rhs[self.st[v] as usize] - self.d[v];
                }
            }
            inp.clear();
            inp.extend(rhs.into_iter().enumerate().filter(|e| e.1 != R::ZERO));
        }
        let out = self.mirror.btran(&inp, false, true);
        self.scr_binp = inp;
        if crate::dump::on("dual") {
            crate::dump::emit(
                "dual",
                serde_json::json!({"refactor": self.refactorizations, "iter": self.iters, "out": crate::dump::sparse_r(&out)}),
            );
        }
        let mut y = std::mem::take(&mut self.last_y);
        y.clear();
        y.resize(p.m, R::ZERO);
        for &(i, v) in out.iter() {
            y[i] = v;
        }
        self.last_y = y;
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
        cols: crate::lu_order::Cols<'_>,
        f: &Mirror<R>,
    ) {
        if self.ptype == 0 {
            return;
        }
        // A repeated (row, column) entry can repeat a row in a singleton L column; the old lists deduplicate it, the factor's do not.
        if crate::dump::on("p68") || self.dup_entries {
            return self.reset_pricing_ref(st, cols);
        }
        // `f` is the factor just built by `lu_attempt`, before it replaces `self.mirror`: no update since `fresh`.
        let w = pricing_weights(self.m, st, cols, f);
        self.set_pricing_weights(w);
    }

    /// The cap and store of a reset's weights.
    fn set_pricing_weights(&mut self, mut w: Vec<f64>) {
        self.ptype = cap_pricing_weights(self.ptype, self.m, &mut w);
        self.p68 = w;
    }

    /// `reset_pricing` with the forward/backward lists materialized, for the `p68` dump (which carries them).
    fn reset_pricing_ref(
        &mut self,
        st: &crate::lu_order::Structure<R>,
        cols: crate::lu_order::Cols<'_>,
    ) {
        let (fwd, bwd, mut w) = pricing_weights_lists(self.m, st, cols);
        self.ptype = cap_pricing_weights(self.ptype, self.m, &mut w);
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

    /// `0xbdc10`: flags of all rows and the list of the infeasible ones sorted by priority (descending), in one pass over the rows.
    ///
    /// Returns the number of infeasible rows and the number of those whose flag was not 0 before (the rows that became infeasible outside the list).
    fn build_walk_list(&mut self) -> (usize, usize) {
        #[cfg(any(debug_assertions, feature = "factor-shadow"))]
        self.check_bounds();
        let ptol = self.ptol();
        let m = self.m;
        let mut rows = std::mem::take(&mut self.glue.spare_list);
        let mut keys = std::mem::take(&mut self.glue.keys);
        let mut flags = std::mem::take(&mut self.flags);
        let counts = rebuild_rows::<R>(&self.xb[..m], &self.blo[..m], &self.bup[..m], ptol, &mut flags, &mut rows, &mut keys, |r, fr| {
            (-self.pri(r, fr)).to_f64()
        });
        self.flags = flags;
        ref_sort(&mut keys, &mut rows);
        self.walk_len = rows.len();
        self.glue.spare_list = std::mem::replace(&mut self.list, rows);
        self.glue.keys = keys;
        self.cursor = 0;
        self.pend = self.cap;
        counts
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
                .upd_limit()
                .unwrap_or(((self.m as f64).powf(0.7) as usize).max(20));
            self.upd_limit = Some(lim.min(cur));
            self.mirror.set_upd_limit(self.upd_limit);
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
                if R::QUAD {
                    rel = rel.max(f / (x.abs() + R::ONE));
                } else {
                    let denominator = x.abs() + R::ONE;
                    if f > rel * denominator {
                        rel = f / denominator;
                    }
                }
            }
        }
        let s3 = self.xb_residual.max(R::c(tol * 1e-4));
        any && (rel < s3 || (R::ONE / R::c(10.0)) * xm * s3 > a)
    }

    /// The leaving row (`0xccbc0`), `None` when no row is infeasible.
    fn price(&mut self) -> Option<usize> {
        if let Some(&j) = self.forced.last() {
            // The forced-row stack (`0xccbc0..0xccd1f`): the top entry is the leaving variable while it is basic with d != 0; a stale entry is popped, and the
            // pop of the last one requests a refactorization (10).
            if self.st[j] >= 0 && self.d[j] != R::ZERO {
                return Some(self.st[j] as usize);
            }
            self.forced.pop();
            self.forced_stale = true;
            if self.forced.is_empty() {
                self.request = 10;
            }
            return None;
        }
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
                let old_len = self.walk_len.min(self.list.len());
                let old = self.list[..old_len]
                    .iter()
                    .filter(|&&i| self.flags[i] == 0 && self.infeas(i) > ptol)
                    .count();
                // The row counts (`infeas`, `fresh`) come from the pass of the build, which reads each flag before it sets it; `old` reads the flags of the old
                // list, so it is counted before the build.
                let (infeas, fresh) = self.build_walk_list();
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
        let mut kept = std::mem::take(&mut self.glue.spare_list);
        kept.clear();
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
        self.glue.spare_list = std::mem::replace(&mut self.list, kept);
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
            let fr = self.infeas_real(r);
            let f = fr.to_f64();
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
                self.prio[r] = self.pri(r, fr);
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
    /// Returns entries `(variable, value)` with structurals first, the structural result count, and the occurrence/distinct-column counts before filtering.
    fn pivot_row(
        &self,
        rho: &[(usize, R)],
        leave: usize,
        above: bool,
    ) -> (Vec<(usize, R)>, usize, usize, usize) {
        // No reentrancy: the workspace is borrowed only for the duration of this call, which calls nothing else of `self`.
        let mut ws = self.prow_ws.borrow_mut();
        pivot_row_ws(&mut ws, &self.rrows, &self.slack_nb, &self.st, &self.p.cols, self.n, self.m, rho, leave, above)
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
        mut rows: Vec<usize>,
        mut vals: Vec<R>,
    ) -> (Vec<usize>, Vec<R>, R, f64, bool) {
        let p = self.p;
        let sdir = if above { -1.0 } else { 1.0 };
        let mut guard = self.fw_ws.borrow_mut();
        let FwWs { dense, seen, order, list, sq } = &mut *guard;
        if dense.len() != self.m {
            *dense = vec![R::ZERO; self.m];
            *seen = vec![false; self.m];
        }
        list.clear();
        list.reserve(row.len() + 1);
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
        let w = rho_norm_with(list, sq);
        order.clear();
        let mut work = 4.0 * list.len() as f64;
        for &(j, v) in list.iter() {
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
        // The result lists are the caller's retained buffers (capacities only).
        rows.clear();
        vals.clear();
        rows.reserve(order.len());
        vals.reserve(order.len());
        let order_len = order.len();
        // The scratch goes back to all zero / all false through the touched list.
        for &i in order.iter() {
            if dense[i].abs() > R::c(1e-13) || (dense_route && dense[i] != R::ZERO) {
                rows.push(i);
                vals.push(dense[i]);
            }
            dense[i] = R::ZERO;
            seen[i] = false;
        }
        work += 2.0 * order_len as f64 + 2.0 * rows.len() as f64;
        (rows, vals, w, work, dense_route)
    }

    /// Work of the row product (`0x134d90`, site `0x135e92`, same estimate as the primal path's `row_product_work`).
    fn row_product_work(&mut self, rho: &[(usize, R)], nslack: usize, p5: usize, d: usize) {
        let n = self.n;
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
    // GUROBI: with `ctx+0x468 == 0` the reference takes `0xc6900` instead; not modelled, no witness (never hit by `replay/gdb_dual_branchhits.py`).
    ///
    /// # Safety
    /// Every index of `row` must be below `n + m` (rows of `pivot_row_ws` are; the candidate loops read `st` and `d` unchecked).
    /// With `zero` the candidate loops also zero the entries below `ALPHA_TOL` in magnitude (`self.rbuf.edge` tells that the zeroing is not complete).
    unsafe fn ratio_test(&mut self, row: &mut [(usize, R)], zero: bool, esi: usize, leave: usize, above: bool, xbr: R) -> Ratio {
        let mut bufs = std::mem::take(&mut self.rbuf);
        let result = unsafe { self.ratio_test_with(&mut bufs, row, zero, esi, leave, above, xbr) };
        self.rbuf = bufs;
        result
    }

    /// # Safety
    /// As `ratio_test`.
    unsafe fn ratio_test_with(
        &mut self,
        bufs: &mut RatioBufs<R>,
        row: &mut [(usize, R)],
        zero: bool,
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
        let RatioBufs { ratio_arr, con, lst, kept, edge } = bufs;
        assert!(self.st.len() >= self.n + self.m && self.d.len() >= self.n + self.m, "st or d shorter than n + m");
        let variant = self.level != 0 && self.iters + self.it_base > self.m;
        // SAFETY: the lengths were checked above (n + m); the row indices are below n + m (the contract of this function).
        let Cand { bound, best_i, bestc, xm13, xm9, edge: e } = unsafe {
            if zero {
                ratio_candidates::<R, true>(&p.lo, &p.up, &self.st, &self.d, row, esi, delta, xmm7, xmm12, variant, ratio_arr, con, lst)
            } else {
                ratio_candidates::<R, false>(&p.lo, &p.up, &self.st, &self.d, row, esi, delta, xmm7, xmm12, variant, ratio_arr, con, lst)
            }
        };
        *edge = e;
        let row: &[(usize, R)] = row;
        kept.clear();
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
                for &i in lst.iter() {
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
            // An empty list at a nonzero level on a fresh factor, the infeasibility beyond the (artificial-bound) primal tolerance: `0xcc061` -> `0xcc1ca`.
            if self.level != 0 && self.updates == 0 && R::c(-self.ptol()) > xmm10 {
                self.fresh_absorb = true;
            }
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
                    // GUROBI: both arms (the cost shift `0xcb0be`, the free-variable zeroing `0xc8739`) follow the disassembly; no witness (never hit by
                    // `replay/gdb_dual_branchhits.py`, PLAN-dual.md "Unexercised in-scope branches").
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
            for &i in lst.iter() {
                let (j, av) = row[i];
                mx = maxr(mx, (delta + self.d[j]) / av);
            }
            xm3 = mx;
            xm15 = mx;
            kept.clear();
            kept.extend_from_slice(lst);
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
    /// The choice between the forced ratio `0xcc6a0` and the normal one (`0xc53a7..0xc54d3`): the normal ratio needs the leaving variable infeasible on the side its
    /// d sign wants (d > 0: below its lower bound by more than the primal tolerance; d <= 0: above its upper bound by more); an infinite bound there, or no such
    /// infeasibility, takes the forced one.
    // GUROBI: `ctx+0x198 != 0` also takes the forced ratio (`0xc53a7..0xc54d3`); not modelled, no witness (0 in every capture).
    fn forced_ratio_chosen(&self, leave: usize, xbr: R) -> bool {
        let tol = R::c(self.ptol());
        if self.d[leave] > R::ZERO {
            let lo = self.p.lo[leave];
            -BIG >= lo || xbr >= R::c(lo) - tol
        } else {
            let up = self.p.up[leave];
            up >= BIG || R::c(up) + tol >= xbr
        }
    }

    /// The forced ratio test `0xcc6a0` on the pivot row (reference frame, entries > 0 qualify): a Harris pass with the bound |d_p| (candidates with
    /// bound * alpha >= d_j; those above the pivot floor tighten the bound to (d_j + delta) / alpha, the others go to a small list), then the entry of the largest alpha with
    /// d_j / alpha <= bound. `None`: the bound stays at |d_p|, the entering variable is p itself. The small candidates that the step would drive negative have their cost shifted by
    /// (theta alpha - d_j) (1000 + draw % 200) / 1000 (`0xcc9b5`).
    fn forced_ratio(&mut self, row: &[(usize, R)], esi: usize, leave: usize) -> Option<usize> {
        let _ = leave;
        let dp = self.d[leave].abs();
        let delta = R::c(dual_tol());
        let floor = R::c(PIVOT_FLOOR);
        let mut bound = dp;
        let (mut big, mut small): (Vec<usize>, Vec<usize>) = (Vec::new(), Vec::new());
        for (i, &(j, v)) in row.iter().enumerate() {
            let ok = if i < esi {
                self.st[j] < 0
            } else {
                self.st[j] < 0 && self.st[j] >= UB
            };
            if !(v > R::ZERO && ok) {
                continue;
            }
            let dj = self.d[j];
            if bound * v < dj {
                continue;
            }
            if v > floor {
                big.push(i);
                let nb = (dj + delta) / v;
                bound = if nb < bound { nb } else { bound };
            } else {
                small.push(i);
            }
        }
        let (qe, theta) = if bound < dp {
            let mut best = floor;
            let mut q = None;
            for &i in &big {
                let (j, v) = row[i];
                if bound >= self.d[j] / v && v > best {
                    best = v;
                    q = Some(i);
                }
            }
            let qe = q?;
            (Some(qe), self.d[row[qe].0] / row[qe].1)
        } else {
            (None, dp)
        };
        // GUROBI: the small-candidate shift loop `0xcc9b5` follows the disassembly; no witness (STATE.md, "Known unmodelled edges").
        for &i in &small {
            let (j, v) = row[i];
            if theta > self.d[j] / v {
                let rem = (self.rng_int() as i32) % 200;
                let sh = (theta * v - self.d[j]) * R::c((rem + 1000) as f64) / R::c(1000.0);
                if R::c(0.001) * delta > sh {
                    self.raise_pert(1);
                } else {
                    self.raise_pert(2);
                }
                let s64 = sh.to_f64();
                if self.st[j] == LB {
                    self.cshift[j] += s64;
                    self.cost[j] += s64;
                } else {
                    self.cshift[j] -= s64;
                    self.cost[j] -= s64;
                }
                self.d[j] += sh;
            }
        }
        qe
    }

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
        kept0: &[usize],
        chosen: usize,
        xm5: R,
        mut xm15: R,
        nflips: usize,
    ) -> (usize, R) {
        let mut xm3;
        let mut xm5 = xm5;
        let delta = R::c(dual_tol());
        let mut chosen = chosen;
        // The list is copied only when a pass replaces it.
        let mut kept_own: Option<Vec<usize>> = None;
        for _ in 0..8 {
            let kept: &[usize] = kept_own.as_deref().unwrap_or(kept0);
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
            for &e in kept {
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
                // The empty-list branch combines the re-pass and selection work in one store (`0xc8b62`). A list of one entry never reaches it
                // (GREENBEA dual-controlled iteration 150: no `0xc8dd3` hit and no booking; the captured hits all have two or more entries).
                if kept.len() > 1 {
                    self.work.add(
                        0xc8b67,
                        (3 * kept.len() + 4 * nflips) as f64,
                        crate::work::CTX,
                    );
                }
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
                kept_own = Some(next);
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

/// The running sum of the ratio test's selection. The double path (`0xcba58..0xcbc3a`) uses vectorized order: fewer
/// than `peel + 8` entries are added one by one. Otherwise `peel` entries (one when the first is at an odd index of the 16-byte aligned array, so that the
/// pairs are aligned) go in one by one, then blocks of eight enter four two-lane accumulators (lane 0 of the first one starts at `s`), which are combined
/// as `(a0 + a1) + (a2 + a3)` and then lane 0 + lane 1, and the remaining entries follow one by one. The quad path (`0x22342b..0x2236d4`) uses eight scalar lanes without a peel, reduces them in order, then walks the remainder backward.
fn accumulate<R: Real>(s: R, c: &[R], peel: usize, sub: bool) -> R {
    let step = |a: R, x: R| if sub { a - x } else { a + x };
    let n = c.len();
    let mut s = s;
    if !R::LANES {
        let vend = n & !7;
        if vend != 0 {
            // Seed the first x87 lane with the running sum.
            let mut a = [R::ZERO; 8];
            a[0] = s;
            for block in c[..vend].chunks_exact(8) {
                for i in 0..8 {
                    a[i] = step(a[i], block[i]);
                }
            }
            s = a[0];
            for &lane in &a[1..] {
                s = s + lane;
            }
        }
        for &x in c[vend..].iter().rev() {
            s = step(s, x);
        }
        return s;
    }
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
    /// A forced row of the DStart start was resolved without a basis change (`0xc660e`): the iteration counts, the request (10 after the last entry) follows.
    Shifted,
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
    /// The position in `row` of the first entry of `q`, and of the first entry of the leaving variable (`None`: not in the row).
    q_pos: usize,
    leave_pos: Option<usize>,
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
    /// Keeps the pivot's `col` and `tau` buffers and its box for the next iteration (contents are cleared when taken).
    fn reclaim(&mut self, mut pv: Box<Pivot<R>>) {
        self.glue.col = std::mem::take(&mut pv.col);
        self.glue.tau_col = std::mem::take(&mut pv.tau);
        self.glue.pivot = Some(pv);
    }

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
        if std::mem::take(&mut self.forced_stale) {
            return Step::Reject;
        }
        let forced_mode = !self.forced.is_empty();
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
        // With the forced-row stack filled the side comes from the sign of the leaving variable's d (`0xc516c`: d <= 0 leaves at the upper bound).
        let above = if forced_mode {
            self.d[leave] <= R::ZERO
        } else {
            xbr > R::c(p.up[leave])
        };
        // Pivot-row BTRAN, drop of the small entries.
        let rho_full = self.mirror.btran(&[(r, R::ONE)], true, false);
        let drop_tol = R::c(DROP_TOL);
        let mut rho: Vec<(usize, R)> = std::mem::take(&mut self.glue.rho);
        rho.clear();
        rho.reserve(rho_full.len());
        for &e in &rho_full {
            if e.1.abs() > drop_tol {
                rho.push(e);
            }
        }
        self.p60 = rho.len() as f64;
        self.work
            .add_from(0x13174d, 0xc51dd, self.mirror.work(), crate::work::FAC);
        let (mut row, esi, occurrences, distinct) = self.pivot_row(&rho, leave, above);
        self.row_product_work(&rho, row.len() - esi, occurrences, distinct);
        self.p58 += occurrences as f64;
        self.p50d += occurrences as f64;
        let forced_route = forced_mode && self.forced_ratio_chosen(leave, xbr);
        // The zeroing of the small row entries is made by the candidate loops of the ratio test (`zeroed`), unless a dump of the row or the forced route
        // needs the raw values.
        let mut zeroed = false;
        let rt = if forced_route {
            match self.forced_ratio(&row, esi, leave) {
                Some(qe) => Ratio {
                    q: Some(qe),
                    flips: Vec::new(),
                    reject: None,
                },
                None => {
                    // q = p (`0xc660e`): no basis change; every d of the row moves by -|d_p| times its entry, d_p becomes 0, the iteration counts.
                    self.forced.pop();
                    let last = self.forced.is_empty();
                    let dp = self.d[leave].abs();
                    for &(j, a) in &row {
                        self.d[j] = self.d[j] + a * (-dp);
                    }
                    self.d[leave] = R::ZERO;
                    self.iters += 1;
                    if last {
                        self.request = 10;
                    }
                    return Step::Shifted;
                }
            }
        } else {
            // SAFETY: the row comes from `pivot_row_ws` (indices below n + m).
            zeroed = snap.is_none();
            unsafe { self.ratio_test(&mut row, zeroed, esi, leave, above, xbr) }
        };
        if std::mem::take(&mut self.fresh_absorb) && !self.no_art_exit {
            if let Some(mut sn) = snap {
                sn["r"] = serde_json::json!(r);
                sn["var"] = serde_json::json!(leave);
                sn["rho"] = crate::dump::sparse_r(&rho);
                sn["row"] = crate::dump::sparse_r(&row);
                sn["q"] = serde_json::json!(-1);
                *dump = Some(sn);
            }
            // `ctx+0x290` = 1: the post-iteration pass counts the variables at an artificial bound and, with `ctx+0x35c` != 0, ends the pass with status -1004.
            self.art_exit = true;
            self.dual_unbounded();
            return Step::Mismatch;
        }
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
                if !forced_route
                    && !was_retry
                    && self.st[row[qe].0] != FREE
                    && row[qe].1 > R::c(-PIVOT_FLOOR) =>
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
            // GUROBI: the pop of the forced-row stack at the ratio-reject returns (`0xc67a1`; here and at the no-entering return below) follows the
            // disassembly; no witness (STATE.md, "Known unmodelled edges").
            if forced_mode {
                self.forced.pop();
            }
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
            if forced_mode {
                self.forced.pop();
            }
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
        // GUROBI: the reference has a loop raising `ctx+0xf8` from 1e-13; rose keeps 1e-13, no witness (PLAN-dual.md, "Remaining extended-precision rules").
        let obj_increment = {
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
            ob
        };
        if !flips.is_empty() {
            self.work
                .add(0xc561a, 7.0 * flips.len() as f64, crate::work::CTX);
        }
        if let Some(s) = rec.as_mut() {
            s["q"] = serde_json::json!(q);
            s["alpha_row"] = serde_json::json!(aq_row.to_f64());
            s["flips"] = serde_json::json!(flips);
        }
        if forced_mode {
            // The pop of the forced-row stack (`0xc5738`) precedes the FTRAN and the pivot check.
            self.forced.pop();
        }
        // The reference's two-vector FTRAN: the pivot column first, then tau = B^-1 rho at the column's positions only.
        let mut rows_q = std::mem::take(&mut self.glue.rows_q);
        let mut vals_q = std::mem::take(&mut self.glue.vals_q);
        rows_q.clear();
        rows_q.extend(p.fcols[q].iter().map(|&(i, _)| i));
        vals_q.clear();
        vals_q.extend(p.fcols[q].iter().map(|&(_, a)| R::c(a)));
        let col_out = self.mirror.ftran(&rows_q, &vals_q);
        // Mode 3: the second vector is not rho but the combination of the columns of the row's entries that were basic when the weights restarted
        // (`0xc57d8`, `0x11f620`), and its squared norm is the weight of the pivot row.
        let (rho_rows, rho_vals, w_row, vec_work, dense_route) = if self.pmode >= 2 {
            let rr = std::mem::take(&mut self.glue.rho_rows);
            let rv = std::mem::take(&mut self.glue.rho_vals);
            self.framework_rhs(&row, leave, above, rr, rv)
        } else {
            let mut rr = std::mem::take(&mut self.glue.rho_rows);
            let mut rv = std::mem::take(&mut self.glue.rho_vals);
            rr.clear();
            rr.extend(rho.iter().map(|e| e.0));
            rv.clear();
            rv.extend(rho.iter().map(|e| e.1));
            (
                rr,
                rv,
                rho_norm_with(&rho, &mut self.glue.sq),
                rho.len() as f64,
                false,
            )
        };
        let mut needed = std::mem::take(&mut self.glue.needed);
        needed.clear();
        needed.extend(col_out.iter().map(|e| e.0));
        // Mode 3: the framework vector is a dense vector (count -1) once it has more than m/8 entries (SC50A calls 11, 14..; AFIRO call 12: 4 of m = 27;
        // 6 of 50 stay sparse); its L stage then runs dense.
        let dense_in = self.pmode == 3 && (dense_route || rho_rows.len() * 8 > m);
        // The mode 2 solves the pivot column alone (`0x11f620` calls the plain FTRAN `0x1205d0` unless the mode is 1 or 3); the two-vector FTRAN adds the result count
        // of the column to the double's counter `P+0x58` (`0x120238`).
        let tau_out = if self.pmode == 2 {
            Vec::new()
        } else {
            let t = self
                .mirror
                .ftran_needed(&rho_rows, &rho_vals, &needed, dense_in);
            let first = if self.mirror.last_dense() {
                1 + col_out
                    .iter()
                    .filter(|e| e.0 != r && e.1.abs() > R::c(DROP_TOL))
                    .count()
            } else {
                col_out.len()
            };
            self.p58d += (first + self.mirror.tau_fin()) as f64;
            t
        };
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
            let (e1, u1) = self.mirror.work_parts();
            if self.pmode == 2 {
                // The plain FTRAN of the column alone (`0x1205d0` -> `0x120673`).
                self.work
                    .add_from(0x1219d6, 0x120673, e1 + u1, crate::work::FAC);
            } else {
                let est = match self.mirror.joint_u(&needed) {
                    Some(j) => (e1 + self.mirror.tau_l()) + j,
                    None => ((e1 + self.mirror.tau_l()) + u1) + self.mirror.tau_u(),
                };
                self.work
                    .add_from(0x129671, 0x120200, est, crate::work::FAC);
            }
            // A slack column (`ctx+0x2fc >= n`, `0x11f788`) is not copied from the matrix: its `4 nnz` term is absent.
            let copied = if q < self.n { rows_q.len() } else { 0 };
            self.work
                .add(0x120114, 4.0 * copied as f64 + vec_work, crate::work::CTX);
        }
        let mut tau_dense = std::mem::take(&mut self.glue.tau);
        if tau_dense.len() != m {
            tau_dense = vec![R::ZERO; m];
        }
        for &(i, v) in &tau_out {
            tau_dense[i] = v;
        }
        let alpha = col_out.iter().rfind(|e| e.0 == r).map_or(R::ZERO, |e| e.1);
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
            // factor ends with status -1005 instead (`0xc5f43`, "Switch to quad due to numerics": the controller reruns the dual in quad precision, `0xcec20`; `Step::SwitchQuad` below).
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
        if !R::QUAD {
            // `0xc5d46`: the accepted pivot adds its objective step to the maintained objective.
            self.obj += obj_increment.to_f64();
            if crate::dump::on("stallobj") {
                crate::dump::emit(
                    "stallobj",
                    serde_json::json!({"ev":"pivot", "it":self.iters, "step":format!("{:016x}", obj_increment.to_f64().to_bits()),
                        "obj":format!("{:016x}", self.obj.to_bits())}),
                );
            }
        }
        if R::QUAD {
            if let Some(old) = self.working_obj {
                let new = (R::c(old) + obj_increment).to_f64();
                self.working_obj = Some(new);
                #[cfg(test)]
                record_quad_obj(new);
                if crate::dump::on("qobj") {
                    crate::dump::emit(
                        "qobj",
                        serde_json::json!({"kind":"pivot", "it":self.iters, "old":old, "delta":obj_increment.to_f64(), "value":new}),
                    );
                }
            }
        }
        // Column without the pivot entry: the last entry moves into the pivot's slot.
        // A dense FTRAN result (count -1) is compacted by an ascending scan that skips the pivot and the entries up to the drop tolerance (`0x12013c`).
        let mut col: Vec<(usize, R)> = std::mem::take(&mut self.glue.col);
        col.clear();
        col.extend_from_slice(&col_out);
        if self.mirror.last_dense() {
            col.retain(|e| e.0 != r && e.1.abs() > R::c(DROP_TOL));
        } else if let Some(k) = col.iter().position(|e| e.0 == r) {
            col.swap_remove(k);
        }
        let mut tau: Vec<R> = std::mem::take(&mut self.glue.tau_col);
        tau.clear();
        tau.extend(col.iter().map(|e| tau_dense[e.0]));
        for &(i, _) in &tau_out {
            tau_dense[i] = R::ZERO;
        }
        self.mirror.give_list(tau_out);
        self.glue.tau = tau_dense;
        // Each accepted pivot and each bound flip advances the shared clock.
        crate::simplex::TOTAL_ITERS.with(|c| c.set(c.get() + 1.0 + flips.len() as f64));
        // Bound flips: x_B and the statuses (`0xcd720`).
        if !flips.is_empty() {
            let mut rhs = std::mem::take(&mut self.glue.rhs);
            if rhs.len() != m {
                rhs = vec![R::ZERO; m];
            }
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
                    .iter()
                    .copied()
                    .filter(|&i| rhs[i].abs() > R::c(DROP_TOL))
                    .collect()
            };
            let vals: Vec<R> = rows.iter().map(|&i| rhs[i]).collect();
            // On the dense route the estimate is the work of the real solve below; otherwise the replay of the dense FTRAN.
            let from_solve = dense_in && self.mirror.dense_route(rows.len());
            #[cfg(any(test, debug_assertions, feature = "factor-shadow"))]
            let replay = dense_in.then(|| self.mirror.ftran_dense_work(&rhs));
            let dense_est = if dense_in && !from_solve {
                Some(self.mirror.ftran_dense_work(&rhs))
            } else {
                None
            };
            for &i in &touched {
                rhs[i] = R::ZERO;
            }
            self.glue.rhs = rhs;
            let sol = self.mirror.ftran_keep(&rows, &vals, dense_in);
            let dense_est = if from_solve {
                Some(self.mirror.dense_work_of_solve(self.mirror.keep_work(), rows.len()))
            } else {
                dense_est
            };
            #[cfg(any(test, debug_assertions, feature = "factor-shadow"))]
            if from_solve {
                assert_eq!(dense_est.map(f64::to_bits), replay.map(f64::to_bits), "dense work of the solve differs from the replay");
            }
            self.work.add_from(
                0x1219d6,
                0x120673,
                dense_est.unwrap_or(self.mirror.keep_work()),
                crate::work::FAC,
            );
            if dense_in || self.mirror.keep_dense() {
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
        self.set_head(r, q);
        // Row-wise copy: in a row holding both columns the leaving one takes the entering one's place; elsewhere the entering one is deleted (the last
        // entry moves into its slot) and the leaving one is appended.
        let leave_in = leave < n && p.lo[leave] != p.up[leave];
        if q < n {
            for &(i, _) in &p.cols[q] {
                if let Some(kk) = self.rrows.position(i, q) {
                    if let Some(a) = leave_in
                        .then(|| p.cols[leave].iter().find(|e| e.0 == i))
                        .flatten()
                        .map(|e| e.1)
                    {
                        self.rrows.set(i, kk, leave, a);
                    } else {
                        self.rrows.swap_remove(i, kk);
                    }
                }
            }
        }
        if leave_in {
            for &(i, _) in &p.cols[leave] {
                if self.rrows.position(i, leave).is_none() {
                    let a = p.cols[leave].iter().find(|e| e.0 == i).unwrap().1;
                    self.rrows.push(i, leave, a);
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
        self.pending_alpha = Some(col_out);
        // The ratio test zeroes row entries of magnitude at most 1e-10 (`0xc7747..0xc7769`, `0xcc606`; bh50022 has an entry exactly -1e-10).
        // The candidate loops have zeroed the entries below the tolerance in magnitude; the entries of exactly the tolerance (`rbuf.edge`) were left raw for
        // `framework_rhs` and are zeroed here.
        if !zeroed || self.rbuf.edge {
            for e in row.iter_mut() {
                if e.1.abs() <= R::c(ALPHA_TOL) {
                    e.1 = R::ZERO;
                }
            }
        }
        // First-match positions of `q` and the leaving variable: the row's entries are unique by variable unless `rho` repeats a slack index.
        let (q_pos, leave_pos) = {
            let ws = self.prow_ws.borrow();
            if ws.slack_dup {
                (row.iter().position(|e| e.0 == q).unwrap(), row.iter().position(|e| e.0 == leave))
            } else {
                debug_assert_eq!(row[qe].0, q);
                (qe, ws.leave_pos)
            }
        };
        // The lists of this iteration go back to the glue (capacities only), the pivot reuses its box.
        self.mirror.give_list(rho_full);
        self.glue.rho = rho;
        self.glue.rows_q = rows_q;
        self.glue.vals_q = vals_q;
        self.glue.rho_rows = rho_rows;
        self.glue.rho_vals = rho_vals;
        self.glue.needed = needed;
        let pivot = Pivot { r, leave, q, row, q_pos, leave_pos, flips, col, tau, alpha, w: w_row, t0, above };
        Step::Pivot(match self.glue.pivot.take() {
            Some(mut b) => {
                *b = pivot;
                b
            }
            None => Box::new(pivot),
        })
    }
}

impl<'a, R: Real> Dual<'a, R> {
    fn snapshot(&self, k: usize) -> serde_json::Value {
        serde_json::json!({
            "c2d8": dual_infeas_sum(), "k": k, "lvl": self.level, "req": self.request, "it": self.iters as f64, "context_it": (self.it_base + self.iters) as f64, "clock": crate::simplex::TOTAL_ITERS.with(|c| c.get()), "head": self.head, "xB": self.xb.iter().map(|v| v.json()).collect::<Vec<_>>(), "d": self.d.iter().map(|v| v.json()).collect::<Vec<_>>(), "W": self.w.iter().map(|v| v.json()).collect::<Vec<_>>(),
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
    fn set_level_work(&mut self, level_was_set: bool) -> f64 {
        let (n, m) = (self.n as f64, self.m as f64);
        let mut est = 2.0 * m + (m + 2.0 * n);
        if level_was_set {
            est += 0.5 * (n + m);
        }
        est += if (self.request == 4 && level_was_set) || (self.request == 6 && self.pert == 0) {
            2.0 * n + (n + m)
        } else {
            1.5 * n
        };
        if self.pert != 0 {
            est += 2.0 * m + n;
        }
        est
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
        // The double post-refactor pass (`0xe5d54..0xe5ead`) chooses the bound-count path when the previous level (`ctx+0x2cc`) is nonzero;
        // otherwise request 4 shifts wrong-signed reduced costs even when `ctx+0x338` was zero (randlp5802).
        if self.request == 1 || (self.request == 4 && !R::QUAD && self.level != 0) {
            return;
        }
        // `add_shift_q` and `raise_pert` leave the dual tolerance alone.
        let tol = R::c(dual_tol());
        for j in 0..self.n + self.m {
            // Basic and fixed variables (and any other status) never shift.
            let s = self.st[j];
            if s >= 0 || s < FREE {
                continue;
            }
            let d = self.d[j];
            let shift = if s == LB {
                (d < -tol).then_some(-d)
            } else if s == UB {
                (d < -tol).then_some(d)
            } else {
                (d.abs() > tol).then_some(-d)
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
    fn set_level(&mut self, old_level: usize) -> Result<(), SolveError> {
        let p = self.p;
        let nv = self.n + self.m;
        // The count of the variables at an artificial bound (`at_artificial` on the final status), one pass: with request 4 inside the loop that moves the
        // statuses (a status changes only at its own `j`), else alone.
        let mut level = 0usize;
        if self.request == 4 {
            for j in 0..nv {
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
                level += self.at_artificial(j) as usize;
            }
        } else {
            for ((&s, &lo), &up) in self.st[..nv].iter().zip(&p.lo[..nv]).zip(&p.up[..nv]) {
                level += ((s == UB && up >= BIG) || (s == LB && lo <= -BIG)) as usize;
            }
        }
        self.level = level;
        self.stall_check(old_level)?;
        if !self.flip_select && (self.pert >= 3 || 3.0 * self.iters as f64 > self.m as f64) {
            self.flip_select = true;
        }
        Ok(())
    }

    /// The objective the stall check reads (`ctx+0x118`, read without a solve): the maintained `obj` in the double dual. The quad dual still solves it
    /// afresh through the auxiliary factor: the working costs times the basic solution and the nonbasic values, plus the shift offset, divided by the
    /// objective scale.
    fn stall_objective(&mut self) -> Result<f64, SolveError> {
        if R::QUAD {
            self.fresh_objective()
        } else {
            Ok(self.obj)
        }
    }

    fn fresh_objective(&mut self) -> Result<f64, SolveError> {
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
        self.factor.ftran(p, &mut r)?;
        let mut s = 0.0;
        for (pos, &j) in self.head.iter().enumerate() {
            s += self.cost[j] * r[pos];
        }
        for j in 0..p.n + p.m {
            if self.st[j] < 0 {
                s += self.cost[j] * self.nonbasic_value(j);
            }
        }
        Ok((s + self.obj_off) / self.obj_scale)
    }

    /// The pricing switch after a factorization (`0xf770d..0xf82c5`, release profile only): at every third factorization, with the steepest-edge pricing
    /// still on, more than m + 1000 iterations done and fewer than 100 m structurals, the weights restart at 1 and the mode becomes 3 ("simplex pricing
    /// changes from 1 to 3"). The pricing update formulas stay those of the steepest edge.
    fn pricing_switch(&mut self) {
        // GUROBI: in mode 2 the reference can switch to mode 3 (`0xf7780..0xf7a21`: iterations > 3m, or `ctx+0x35c` > 2 with `ctx+0x178` > 0); rose
        // leaves mode 2 unchanged, no witness (`specs/dual-simplex.md`, pricing modes).
        if !self.switch_enabled || !matches!(self.pmode, 1 | 3) || self.refactorizations % 3 != 0 {
            return;
        }
        let m = self.m as f64;
        if self.pmode == 1 && self.iters as f64 > m + 1000.0 && (self.n as f64) < 100.0 * m {
            self.pmode = 3;
            self.p50d = 0.0;
            self.p58d = 0.0;
            self.dflag = self.st.iter().map(|&t| t < 0).collect();
            for w in self.w.iter_mut() {
                *w = R::ONE;
            }
            return;
        }
        if R::QUAD {
            return;
        }
        // The switch to the mode 2 (`0xf7a54..0xf7bb2`, not nested): enough row products since the start (`P+0x50 >= m/10 + 10`), at most 3 m iterations and
        // n <= 3 m, and the result counts of the pivot columns' FTRANs (`P+0x58`) exceed `M P+0x50`: `M` = 1e100 below 15000 rows; from 15000 rows in the mode 1 it is 1e100 below m/4
        // iterations, 4 from m/4 and 3 from m/3 iterations (the mode 3 always 3), plus 1 after a perturbation (`ctx+0x338 == 3`). The weights restart at 1.
        // GUROBI: the nested case (`ctx+0xb8` != 0) and the `ctx+0x3f0` condition of `0xf7a54..0xf7bb2` are read but not modelled; no witness (only the
        // standalone `witness_dual_r2_mode_switch` reaches the switch).
        if self.p50d < ((self.m / 10) + 10) as f64
            || 3.0 * m < self.iters as f64
            || 3.0 * m < self.n as f64
        {
            return;
        }
        let it = self.iters as f64;
        let base = if self.m < 15000 {
            1e100
        } else if self.pmode == 1 && it < m / 3.0 {
            if it >= m * 0.25 {
                4.0
            } else {
                1e100
            }
        } else {
            3.0
        };
        let mult = if self.m >= 15000 && self.pert == 3 && base < 1e50 {
            base + 1.0
        } else {
            base
        };
        if self.p58d > mult * self.p50d {
            self.pmode = 2;
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

    fn dump_stall_check(&self, obj: f64) {
        if crate::dump::on("stallobj") {
            crate::dump::emit(
                "stallobj",
                serde_json::json!({"ev":"check", "it":self.iters, "nref":self.nref(), "obj":format!("{:016x}", obj.to_bits()),
                    "best":format!("{:016x}", self.best_obj.to_bits()), "prog_on":self.prog_on}),
            );
        }
    }

    /// The stall test of the refactorization post pass (`0xe60cc..0xe6237`) and the cost perturbation it triggers (`0xe6af8..0xe6cd7`, `0xe74b1`).
    fn stall_check(&mut self, old_level: usize) -> Result<(), SolveError> {
        if self.pstep <= 0.0 && self.pert_param <= 0.0 {
            return Ok(());
        }
        // Quad 0x23904b..0x239085: the period is 2, 5, or 10000 according to
        // whether neither, one, or both of the handoff counter and perturbation flag exceed 2.
        let handoffs = R::QUAD && C35C.with(|c| c.get()) > 2;
        let period = match (handoffs, self.pert > 2) {
            (true, true) => 10_000,
            (true, false) | (false, true) => 5,
            (false, false) => 2,
        };
        if self.nref() % period != 1 {
            return Ok(());
        }
        let mut trig = false;
        if self.pstep > 0.0 && old_level <= self.level {
            // More than 90 percent of the iterations since the last check had an entering d of at most delta.
            trig =
                (self.stall - self.stall_snap) as f64 > 0.9 * (self.iters as f64 - self.stall_it);
            if !trig && self.prog_on {
                // No objective progress since the last check.
                let obj = self.stall_objective()?;
                self.dump_stall_check(obj);
                trig = 1e-13 * (obj.abs().min(1e10) + 1.0) + self.best_obj > obj;
                if !trig {
                    self.best_obj = maxsd(obj, self.best_obj);
                }
            }
        } else if self.prog_on {
            let obj = self.stall_objective()?;
            self.dump_stall_check(obj);
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
        Ok(())
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
        self.working_obj = None;
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
        self.perturbed = true;
        self.working_obj = None;
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
        DUAL_PERT_SIZE.with(|c| c.set(self.pstep));
    }

    #[cfg(test)]
    fn lowgate_record(&self, tag: &'static str, req: u32) {
        LOWGATE_TEST.with(|c| {
            if let Some(v) = c.borrow_mut().as_mut() {
                v.push((tag, self.iters, req, self.lu_thr, self.head.clone()));
            }
        });
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
                self.update_weights(pv, false);
                self.pivot_cols = Some(pv.col.len());
                self.refactor(if schedule_due { "scheduled" } else { "request" })?;
            } else {
                self.factor.update_column(pv.r, alpha);
                self.mirror.update(pv.r);
                let uw = self.mirror.take_update_work();
                for &(site, est, scale) in &uw {
                    self.work.add(site, est, scale);
                }
                self.mirror.give_update_work(uw);
                self.touch_pricing(pv.r, &pv.col);
                self.updates += 1;
            }
            // d update (`0xe45c0`).
            // The reference's row list also holds the leaving variable (rose's only when it is a slack).
            let aq = pv.row[pv.q_pos].1;
            // The step divides by the pivot as the double `ctx+0x308` (the row's pivot, stored as a double); the same for the primal step and `ctx+0x300`.
            let theta = self.d[pv.q] / R::c(aq.to_f64());
            // A zero step skips the loop (`0xe4719`): no work.
            if !refactor && theta != R::ZERO {
                let row_cnt = pv.row.len() + usize::from(pv.leave < self.n);
                self.work
                    .add(0xe7000, 3.0 * row_cnt as f64, crate::work::CTX);
            }
            // The iteration clamps the leaving variable's stored d from below at -dual_tol before the step when the forced-row stack is empty (`0x21dd8e..0x21ddcd`, `fcomip`/`jbe`: d <= -tol gives d = -tol; the double iteration has the same `maxsd` at `0xc5d85..0xc5dab`; 6762 call 43: -1.7136e-6 -> -1e-6).
            if !self.forced_ever && self.d[pv.leave] <= R::c(-dual_tol()) {
                self.d[pv.leave] = R::c(-dual_tol());
            }
            let old_leave = self.d[pv.leave];
            // The reference's row holds a slack leaving variable with its computed coefficient (about 1, 0.9999999999999997 in PILOT.JA call 30657), the step on it uses that.
            let lcoef = pv.leave_pos.map_or(R::ONE, |p| {
                let e = &pv.row[p];
                if pv.above {
                    -e.1
                } else {
                    e.1
                }
            });
            assert!(self.d.len() >= self.n + self.m, "d shorter than n + m");
            // SAFETY: the row's indices are below n + m (`pivot_row_ws`), d has n + m entries (the flips are checked).
            unsafe { d_row_step(&mut self.d, &pv.row, &pv.flips, theta, &mut self.flip_mark) };
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
            #[cfg(test)]
            if self.iters == 15 && self.request == 4 {
                self.lowgate_record("postfactor", self.request);
            }
            self.recompute_d();
            self.shift_wrong_d();
            let level_before = self.level;
            let level_work = self.set_level_work(level_before != 0);
            self.perturbed = false;
            self.set_level(level_before)?;
            // A cost perturbation made by the stall check inside the same pass adds `2 m + 4 n` to the `0xe6fdd` total (kb2 release, iteration 15:
            // 522.5 = 272.5 + 250 in extended units of 1.3; sctap1, grow7 the same).
            let perturb_work = if self.perturbed && !R::QUAD {
                2.0 * self.m as f64 + 4.0 * self.n as f64
            } else {
                0.0
            };
            self.work
                .add(0xe7000, level_work + perturb_work, crate::work::CTX);
            #[allow(unused_mut)]
            let mut req = self.request;
            self.request = 0;
            self.recompute_xb(req)?;
            #[cfg(test)]
            if self.drift_status && CEILING_NO_ACCEPT.with(|c| c.get()) {
                req = 0;
            }
            if req == 4
                && (!R::QUAD && self.level != 0)
                && self.status.is_none()
                && !self.no_art_exit
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
            #[cfg(test)]
            if !R::QUAD
                && self.iters == 15
                && req == 4
                && self.status.is_none()
                && level_before == 0
                && self.tiny_infeasibility()
            {
                LOWGATE_TEST.with(|c| {
                    if let Some(v) = c.borrow_mut().as_mut() {
                        v.push(("eligible", self.iters, req, self.lu_thr, self.head.clone()));
                        if v.iter().filter(|x| x.0 == "eligible").count() == 1 {
                            self.lu_thr = 0.03125;
                        }
                    }
                });
            }
            if self.status.is_none()
                && if R::QUAD { self.level == 0 } else { level_before == 0 }
                // Both recomputations read `ctx+0x35c` (`0x21a315..0x21a329` quad, `0xc1fe3` double).
                && (req == 4 || C35C.with(|c| c.get()) > 10)
                && (
                    // The low arm (`0xc2025..0xc2044`) needs only a nonempty infeasible list (`[P+0x3c]`, `0xc1ffb`) in the value mode; the acceptance
                    // test that follows it (threshold at least 0.5) is the tiny-infeasibility one (witness_dual_r3_residual_low_small).
                    // GUROBI: the low arm in the list-walk mode is not modelled (rose keeps the tiny-infeasibility test there); no witness, no reference data.
                    (!R::QUAD && self.value_mode && self.n_infeas > 0 && self.lu_thr < 0.5)
                        || self.tiny_infeasibility()
                )
            {
                // At a low double LU threshold the recompute returns the pending request
                // to the driver for another factorization before accepting infeasibilities.
                if !R::QUAD && self.lu_thr < 0.5 {
                    self.lu_thr = (4.0 * self.lu_thr).min(0.5);
                    self.request = req;
                    #[cfg(test)]
                    self.lowgate_record("after", self.request);
                    self.finish_drift_status()?;
                    return self.post_inner(None, false);
                }
                if self.lu_thr >= 0.5 {
                    self.flags.iter_mut().for_each(|f| *f = 2);
                    self.list.clear();
                    self.status = Some(Status::Optimal);
                    #[cfg(test)]
                    if self.drift_status {
                        CEILING_TEST.with(|c| {
                            if let Some(v) = c.borrow_mut().as_mut() {
                                v.push(("accepted", self.iters, req, self.head.clone()));
                            }
                        });
                    }
                    self.drift_status = false;
                    #[cfg(test)]
                    self.lowgate_record("optimal", 0);
                }
            }
            self.finish_drift_status()?;
            // The recompute pass `0xbf5b0` (called from `0xf6d0f`) runs the no-leaving-row handler itself when the pricing is in the list-walk mode
            // (`P+0x40 == 0`) and the list is empty (`P+0x3c == 0`, `0xbffbe`): no pricing call follows (RECIPE call 46).
            // GUROBI: the second identical check `0xc000e` (the `ctx+0x290 == 1` path) is not modelled; its witness `witness_dual_norow` (public DStart)
            // agrees in result without it (`tests/dual_start.rs`), so the branch itself is not compared.
            if !self.value_mode && self.list.is_empty() && self.status.is_none() {
                self.no_row();
            }
        } else if let Some(pv) = pv {
            if self.pmode != 2 {
                self.work
                    .add(0xc0410, 3.0 * pv.col.len() as f64, crate::work::CTX);
            }
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
    fn update_weights(&mut self, pv: &Pivot<R>, fast: bool) {
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
        if self.pmode == 2 {
            self.devex_update(pv, fast);
            return;
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

    /// The mode 2 weight update (`0xc2f8e`, `0xbe7b0`): the weights of the entries of the pivot column grow to `(T / alpha^2) a_i^2` when that is larger, with the
    /// framework counters `ctx+0x400..0x414`: the size is cut to a quarter when it exceeds 30 m, and the framework restarts (weights 1, flags of the nonbasic
    /// variables) when the mean accuracy ratio of the leaving row's weight exceeds 5 or after ten cuts. No update when the row's squared norm `T` is zero.
    fn devex_update(&mut self, pv: &Pivot<R>, fast: bool) {
        let (r, alpha, t) = (pv.r, pv.alpha.to_f64(), pv.w.to_f64());
        let m = self.m;
        if self.dv_f > 30.0 * m as f64 {
            self.dv_k = (self.dv_k as f64 * 0.25) as i32;
            self.dv_f *= 0.25;
            self.dv_g *= 0.25;
            self.dv_c += 1;
        }
        if t == 0.0 {
            return;
        }
        let wr = self.w[r].to_f64();
        let g = if t > wr {
            t / wr + self.dv_g
        } else {
            self.dv_g + wr / t
        };
        let k = self.dv_k + 1;
        // GUROBI: the framework restart branch (`0xc2f8e`, `0xbe7b0`) follows the disassembly; no witness (not known to fire on `witness_dual_r2_mode_switch`),
        // nor for the nested case (`ctx+0xb8` != 0), which is not modelled.
        if g / k as f64 > 5.0 || self.dv_c >= 10 {
            self.dv_k = 0;
            self.dv_f = 0.0;
            self.dv_g = 0.0;
            self.dv_c = 0;
            self.dflag = self.st.iter().map(|&t| t < 0).collect();
            for w in self.w.iter_mut() {
                *w = R::ONE;
            }
            if fast {
                self.work
                    .add(0xc34f5, (2 * m + self.n) as f64, crate::work::CTX);
            }
            return;
        }
        self.dv_k = k;
        let mut f = (m / 50) as f64 + self.dv_f;
        self.dv_g = g;
        let a2 = alpha * alpha;
        if t > a2 {
            self.dv_f = f;
            self.w[r] = R::c((t / alpha) / alpha);
        }
        f += pv.col.len() as f64;
        self.dv_f = f;
        if fast {
            self.work
                .add(0xc0410, 3.0 * pv.col.len() as f64, crate::work::CTX);
        }
        let factor = t / a2;
        for &(i, cv) in &pv.col {
            let v = cv.to_f64();
            let x = (factor * v) * v;
            let cur = self.w[i].to_f64();
            self.w[i] = R::c(if x > cur { x } else { cur });
        }
    }

    fn update_primal_and_weights(&mut self, pv: &Pivot<R>) {
        self.update_weights(pv, true);
        let (r, alpha) = (pv.r, pv.alpha);
        let mut changed = std::mem::take(&mut self.glue.changed);
        changed.clear();
        changed.reserve(pv.col.len() + 1);
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
        self.glue.changed = changed;
    }
}

/// The squared norm of rho as the reference sums it (`0x11f9..0x11fa5`, the value `ctx+0x310`): sequentially up to 6 entries, else in the vectorised loop
/// of eight entries per pass (four accumulators of two lanes, combined pairwise, then lane 0 + lane 1, then the tail).
/// The initial weight of a dense BTRAN result (`0xc392b..0xc3b08`): the squares are summed over all `m` slots of the dense vector (zeros included), so
/// each nonzero falls in the lane of its row index; the alignment peel was never needed in the captures (16-byte aligned buffers).
fn rho_norm_dense_with<R: Real>(rho: &[(usize, R)], m: usize, sq: &mut Vec<R>) -> R {
    sq.clear();
    sq.resize(m, R::ZERO);
    for &(i, v) in rho {
        sq[i] = v * v;
    }
    sum_squares(sq)
}

/// `rho_norm` with the squares built in the caller's buffer.
fn rho_norm_with<R: Real>(rho: &[(usize, R)], sq: &mut Vec<R>) -> R {
    if R::QUAD {
        // The quad kernel (`0x25d3ff..0x25d49a`) squares in extended precision and keeps the running sum in a double: each step is
        // rounded to extended by the add and then to double by the store.
        return rho
            .iter()
            .fold(R::ZERO, |a, e| R::c((a + e.1 * e.1).to_f64()));
    }
    sq.clear();
    sq.extend(rho.iter().map(|e| e.1 * e.1));
    sum_squares(sq)
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
    /// The third restart of a quad pass starts from the slack basis without a crash (the restart handler's `factor+0xc` = 1 store is the double handler's: `fx_witness_dual_c6007`, one crash in the whole run).
    static NO_CRASH: std::cell::Cell<bool> = const { std::cell::Cell::new(false) };
    static SWITCH_PERT: std::cell::Cell<Option<u8>> = const { std::cell::Cell::new(None) };
    /// The stall counters the quad pass inherits from the double pass in the same context (`ctx+0x330`, `ctx+0x334`, and `ctx+0x2b0` relative to the new pass's iteration base):
    /// the quad setup's stall test reads them (`fx_witness_dual_c6007`: 2, 0, 0 at iteration 2).
    static SWITCH_STALL: std::cell::Cell<Option<(i64, i64, f64)>> = const { std::cell::Cell::new(None) };
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
    if CLEANUP_PASS.with(|c| c.get()).is_none() {
        reset_pert_sizes(opts.perturb);
        PRIMAL_DRIFT_COUNT.with(|c| c.set(0));
        C35C.with(|c| c.set(1));
    }
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

thread_local! {
    /// The standalone primal's tolerance-ceiling exit hands its state to the controller (`0xcedb0`): iterations of the aborted run, the LU threshold, the inherited pivot pair.
    static STANDALONE_RESUME: std::cell::Cell<Option<(usize, f64, Option<(f64, f64)>)>> = const { std::cell::Cell::new(None) };
}

/// The controller of a standalone primal after its tolerance-ceiling exit: the first pass is the x87 primal from the crash basis (the pricing0 start of the aborted run kept),
/// the counters of the aborted run survive; what follows is the dual controller's loop (primal -1004 -> x87 dual -> restarts).
pub(crate) fn resume_standalone_controller(
    model: &LpModel,
    opts: &SolveOptions,
    iters: usize,
    thr: f64,
    pair: Option<(f64, f64)>,
) -> Result<SolveResult, SolveError> {
    STANDALONE_RESUME.with(|c| c.set(Some((iters, thr, pair))));
    solve_frame(model, None, opts)
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
    };
    let resume = STANDALONE_RESUME.with(|c| c.take());
    if let Some((iters, thr, _)) = resume {
        carry.iters = iters;
        carry.thr = thr;
    }
    if resume.is_some() {
    } else {
        if CLEANUP_PASS.with(|c| c.get()).is_none() {
            crate::simplex::STALL_PREV.with(|c| c.set(0.0));
            crate::simplex::TOTAL_ITERS.with(|c| c.set(0.0));
            crate::simplex::EARLIER_ITERS.with(|c| c.set(0));
        }
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
        SNAP_PRIMAL.with(|c| c.set(false));
        DRIFT_START.with(|c| *c.borrow_mut() = None);
        SNAP_C40.with(|c| c.set(-100.0));
        PAIR_TOL.with(|c| c.set(3e-6));
        PASS_NO.with(|c| c.set(0));
        C178.with(|c| c.set(0));
        crate::simplex::HANDOFF_ROWS.with(|c| *c.borrow_mut() = None);
        LARGE_XB.with(|c| c.set(false));
        crate::simplex::QUAD_PRIMALS.with(|c| c.set(0));
    }
    // The quad mode (`ctx+0x28`): once the controller switched, every later pass runs the extended-precision twin of the dual.
    let mut quad = resume.is_some();
    let mut switch_start: Option<Basis> = None;
    let mut restart_start: Option<Basis> = None;
    let mut restart_primal = resume.is_some();
    let mut first_resume = resume.is_some();
    // The factor's `+0x0c`: zero supplies the current statuses to setup,
    // one asks setup for a fresh crash basis.
    let mut factor_crash = start.is_none();
    // The controller snapshots the tolerance once (`0x154ccc`), before any
    // singular retries; its primal handoffs scale this original value.
    if resume.is_none() {
        callback_reset();
    }
    loop {
        SINGULAR.with(|c| c.set(None));
        SINGULAR_FROM_PRIMAL.with(|c| c.set(false));
        SINGULAR_TOL_ABORT.with(|c| c.set(false));
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
        // Each pass shares the controller's iteration allowance.
        let mut pass_opts = opts.clone();
        pass_opts.max_iterations = opts.max_iterations.saturating_sub(carry.iters);
        let res = if restart_primal {
            restart_primal = false;
            let resumed = std::mem::take(&mut first_resume);
            let fresh = restart_basis.is_none() && !resumed;
            if fresh {
                set_dual_tol(DUAL_TOL0);
                set_primal_tol(PRIMAL_TOL);
                LARGE_XB.with(|c| c.set(false));
                TOL_RESET_DONE.with(|c| c.set(fresh_pass::<f64>(carry.retries)));
            }
            let mut popts = pass_opts.clone();
            popts.method = crate::simplex::Method::Primal;
            popts.refactor_schedule = None;
            // The tolerance-ceiling restart keeps the aborted run's `pricing0` start (x87 method 0, `ctx+0x3ec` re-initialized at the entry).
            if !resumed {
                popts.pricing0 = false;
                popts.sw_carry = Some(crate::simplex::SwCarry {
                    iters: 0.0,
                    degenerate: 0,
                    snap_degenerate: 0,
                    snap_iters: 0.0,
                    factorizations: 0,
                    level: 0,
                    pert: 0.0,
                });
            } else {
                popts.sw_carry = None;
            }
            crate::simplex::HANDOFF.with(|c| {
                c.set(Some(crate::simplex::Handoff {
                    stale: if resumed {
                        resume.and_then(|r| r.2)
                    } else {
                        None
                    },
                    working_obj: None,
                    nref0: 0,
                    thr: carry.thr,
                    retries: carry.retries as usize,
                    pert: !fresh && SNAP_QUALITY.with(|c| c.get() == 2),
                    pair_tol: PAIR_TOL.with(|c| c.get()),
                    quad: true,
                    dual_level: 0,
                    shifted_exit: false,
                    iters: 0.0,
                    pricing_iters: 0.0,
                    snap: 0.0,
                    c178: C178.with(|c| c.get()),
                    preserve_ratio: resumed,
                    excl: None,
                }))
            });
            crate::simplex::PRIMAL_STALL.with(|c| c.set(false));
            crate::simplex::PRIMAL_SHIFT_END.with(|c| c.set(false));
            RESTORED_PRIMAL.with(|c| c.set(true));
            RESTORE_REPAIR
                .with(|c| c.set(restart_basis.is_some() && SNAP_QUALITY.with(|q| q.get() == 2)));
            let pr = primal_callback(model, restart_basis.as_ref(), &popts);
            RESTORED_PRIMAL.with(|c| c.set(false));
            match pr {
                // GUROBI: the reference installs the quad dual's x on this `PRIMAL_STALL` route; rose does not, no witness (STATE.md, "Known unmodelled edges").
                Ok(r)
                    if crate::simplex::PRIMAL_STALL.with(|c| c.take())
                        || crate::simplex::PRIMAL_SHIFT_END.with(|c| c.take()) =>
                {
                    C35C.with(|c| c.set(c.get() + 1));
                    carry.fac += r.refactorizations;
                    carry.it_base = r.iterations;
                    carry.nref0 = r.refactorizations.max(1);
                    carry.after_pivot = true;
                    let mut dual_opts = pass_opts.clone();
                    dual_opts.max_iterations =
                        pass_opts.max_iterations.saturating_sub(r.iterations);
                    let result = solve_pass::<crate::f80::F80>(
                        model,
                        Some(&r.basis),
                        &dual_opts,
                        &carry,
                        false,
                    );
                    if result.is_ok() {
                        carry.iters += r.iterations;
                    }
                    result
                }
                other => other,
            }
        } else if quad {
            solve_pass::<crate::f80::F80>(
                model,
                start_now,
                &pass_opts,
                &carry,
                restart_basis.is_some(),
            )
        } else {
            solve_pass::<f64>(
                model,
                start_now,
                &pass_opts,
                &carry,
                restart_basis.is_some(),
            )
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
                };
                switch_start = Some(b);
            }
            Err(SolveError::SingularBasis) if carry.retries < 100 => {
                let failed_quad_dual = quad && !SINGULAR_FROM_PRIMAL.with(|c| c.get());
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
                // GUROBI: the other switch conditions of `0x2c8d0` (`ctx+0x188/0x380/0x388`, `env+0x2728 < 2`) are assumed to hold; no witness
                // (PLAN-dual.md, "Remaining extended-precision rules").
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
                // The primal tolerance-exhaustion exit writes -1002 without
                // incrementing ctx+0x184; a singular factor handler increments
                // it before entering the common restart handler.
                let next_retries = pass_retries + u32::from(!SINGULAR_TOL_ABORT.with(|c| c.get()));
                let fresh = next_retries == 3;
                TOL_RESET_DONE.with(|c| c.set(false));
                if fresh {
                    reset_pert_sizes(opts.perturb);
                    factor_crash = true;
                    NO_CRASH.with(|c| c.set(quad));
                    restart_primal = failed_quad_dual && SNAP_PRIMAL.with(|c| c.get());
                } else if context_iters > 0 && SNAP_QUALITY.with(|c| c.get() > 0) {
                    factor_crash = false;
                    let snap = SNAP
                        .with(|c| c.borrow().clone())
                        .or_else(|| INITIAL_SNAP.with(|c| c.borrow().clone()));
                    restart_start = snap.map(|(basic, at_upper)| Basis { basic, at_upper });
                    restart_primal = failed_quad_dual && SNAP_PRIMAL.with(|c| c.get());
                } else if !factor_crash {
                    restart_start = SINGULAR_BASIS.with(|c| c.borrow_mut().take());
                }
                if !quad {
                    // The double dual setup clears ctx+0x2c8 (`0xe7e1b`) before the restarted pass's first pricing rebuild (randlp50154).
                    LARGE_XB.with(|c| c.set(false));
                }
                let (avg_touch, last_cap, inf_count) = if fresh {
                    (0.0, 0, 0)
                } else {
                    (avg_touch, last_cap, inf_count)
                };
                carry = Carry {
                    retries: next_retries,
                    thr: thr.max(if pass_retries == 0 { 0.5 } else { 0.9 }),
                    fac,
                    iters: carry.iters + if restart_primal { context_iters } else { iters },
                    avg_touch,
                    last_cap,
                    inf_count,
                    after_pivot,
                    it_base: 0,
                    nref0: 0,
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

// These two handoffs return ownership to the controller; other historical
// handoff paths still run through their existing recursive routes.
struct HandoffSnapshot {
    iters: usize,
    refactorizations: usize,
    lu_thr: f64,
    avg_touch: f64,
    last_cap: usize,
    inf_count: usize,
    retries: u32,
    it_base: usize,
    nref: usize,
    pert: u8,
    rejected_end: bool,
    level: usize,
    stall_it: f64,
    art_exit: bool,
}

struct PrimalTransfer {
    basis: Basis,
    opts: SolveOptions,
    snapshot: HandoffSnapshot,
    status: Status,
    trace: Vec<TraceEvent>,
    control_status: Option<i32>,
    route: PrimalRoute,
}

enum PrimalRoute {
    Restored,
    Artificial { start: Option<Basis> },
}

enum DriverOutcome {
    Finished(SolveResult),
    ToPrimal(PrimalTransfer),
}

#[derive(Clone, Copy)]
enum ControlTransfer {
    MethodSwitch,
    Restart,
}

impl ControlTransfer {
    fn status(self) -> i32 {
        match self {
            Self::MethodSwitch => -1004,
            Self::Restart => -1002,
        }
    }
}

struct PrimalOutcome {
    result: Result<SolveResult, SolveError>,
    control: Option<ControlTransfer>,
}

fn dispatch_primal(model: &LpModel, basis: &Basis, opts: &SolveOptions) -> PrimalOutcome {
    let result = primal_callback(model, Some(basis), opts);
    // These are the control routes represented by the current primal driver;
    // other native callback statuses have not yet been surfaced.
    let control = match &result {
        Err(SolveError::SingularBasis) => Some(ControlTransfer::Restart),
        Ok(_) if crate::simplex::PRIMAL_STALL.with(|c| c.get()) => {
            Some(ControlTransfer::MethodSwitch)
        }
        _ => None,
    };
    PrimalOutcome { result, control }
}

fn solve_pass<R: Real>(
    model: &LpModel,
    start: Option<&Basis>,
    opts: &SolveOptions,
    carry: &Carry,
    snapshot_restart: bool,
) -> Result<SolveResult, SolveError> {
    match solve_pass_callback::<R>(model, start, opts, carry, snapshot_restart)? {
        DriverOutcome::Finished(result) => Ok(result),
        DriverOutcome::ToPrimal(transfer) => {
            CALLBACK.with(|c| {
                let mut state = c.get();
                state.selector = false;
                c.set(state);
            });
            if crate::dump::on("callback") {
                let route = match &transfer.route {
                    PrimalRoute::Restored => "restored",
                    PrimalRoute::Artificial { .. } => "artificial",
                };
                crate::dump::emit(
                    "callback",
                    serde_json::json!({"event":"transfer", "route":route, "selector":0, "result":0, "status":transfer.control_status, "total":crate::simplex::TOTAL_ITERS.with(|c| c.get())}),
                );
            }
            match transfer.route {
                PrimalRoute::Restored => resume_primal_transfer::<R>(model, opts, transfer),
                PrimalRoute::Artificial { .. } => {
                    resume_artificial_transfer::<R>(model, opts, transfer)
                }
            }
        }
    }
}

fn solve_pass_callback<R: Real>(
    model: &LpModel,
    start: Option<&Basis>,
    opts: &SolveOptions,
    carry: &Carry,
    snapshot_restart: bool,
) -> Result<DriverOutcome, SolveError> {
    if fresh_pass::<R>(carry.retries) {
        // The fresh start of the third restart initialises the context again (`0xcf29b`): the tolerance goes back to the parameter, the counter stays.
        if !TOL_RESET_DONE.with(|c| c.replace(false)) {
            set_dual_tol(DUAL_TOL0);
        }
        set_primal_tol(PRIMAL_TOL);
        LARGE_XB.with(|c| c.set(false));
    }
    let mut callback = CallbackGuard::enter(true, opts);
    let (dual_pert_size, dual_pert_left) = enter_dual_pert_size();
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
            if (FROM_SWITCH.with(|c| c.get()) || (has_restart_snapshot() && carry.retries < 5))
                && crate::simplex::repair_basis(&p, &mut probe).is_some()
            {
                // The basis the quad switch carries over is not a supplied one: when its quad factorization is singular (the double structure finds the same columns
                // dependent), the pass ends before its first iteration and the restart takes the snapshot basis (seed 2706: first quad factorization has `p13` = 1,
                // the next one is the slack basis, `ctx+0x184` = 1).
                if !FROM_SWITCH.with(|c| c.get()) {
                    // The x87 dual run after a primal's -1004 fails its start factorization (`0x248f99`): the stability counter counts it.
                    C178.with(|c| c.set(c.get() + 1));
                }
                SINGULAR_BASIS.with(|c| *c.borrow_mut() = Some(b.clone()));
                SINGULAR.with(|c| {
                    c.set(Some((
                        if R::QUAD { 0 } else { carry.it_base },
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
            if R::QUAD && carry.after_pivot && carry.retries >= 5 && dropped.is_some() {
                // The singular handler counts the in-place repair after the retry limit too (`0x248f99`, `standalone_restart_context` dual-out with status 3: `ctx+0x178` 8 -> 9).
                C178.with(|c| c.set(c.get() + 1));
            }
            b
        }
        None if NO_CRASH.with(|c| c.take()) => Basis {
            basic: (n..n + m).collect(),
            at_upper: vec![false; n + m],
        },
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
    let mut mirror = Mirror::fresh(m, n, &stc, rows.view());
    let mut carry_thr = carry.thr;
    // Growth check of the setup factorization (quad setup `0x246d39..0x246df1`, double driver `0xf351d..0xf3629`): max|diag| / min|diag| above 1e15 (1e12 scaled)
    // raises the threshold to 0.5, then 0.9, and the same columns are factorized again.
    if order.len() == m {
        while carry_thr < 0.9 {
            let (mx, mn) = mirror
                .diag()
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
            mirror = Mirror::fresh(m, n, &stc, rows.view());
        }
    }
    // The auxiliary LU is built on its first read (the deferred-free pass). The quad pass starts from an ill-conditioned basis the double factor may find
    // singular; the quad LU decides (the double factor is only for read-out).
    let factor = AuxiliaryFactor::lazy_floor(order.clone(), if R::QUAD { 0.0 } else { 1e-11 });
    let mut st = vec![LB; n + m];
    for (pos, &j) in order.iter().enumerate() {
        st[j] = pos as i32;
    }
    let cap = ((m as f64).powf(0.7).floor() as i64).max(20);
    let (blo, bup) = basic_bounds::<R>(&p, &order);
    let mut s: Dual<R> = Dual {
        p: &p,
        n,
        m,
        head: order,
        blo,
        bup,
        st,
        xb: vec![R::ZERO; m],
        d: vec![R::ZERO; n + m],
        w: vec![R::ONE; m],
        mirror,
        factor,
        rrows: FlatRows::new(m, n),
        absorb_request: false,
        fresh_absorb: false,
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
        no_art_exit: NO_ART_EXIT.with(|c| c.replace(false)),
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
        automatic_quad: opts.automatic_quad,
        level: 0,
        request: 0,
        iters: 0,
        prio: vec![R::ZERO; m],
        ptype: 0,
        p68: Vec::new(),
        p70: Vec::new(),
        p60: 0.0,
        p58: 0.0,
        p50d: 0.0,
        p58d: 0.0,
        dv_f: 0.0,
        dv_g: 0.0,
        dv_k: 0,
        dv_c: 0,
        flags: vec![2; m],
        list: Vec::new(),
        thr: R::ZERO,
        xb_residual: R::ZERO,
        value_mode: false,
        cursor: 0,
        walk_len: 0,
        slack_nb: vec![false; m],
        prow_ws: Default::default(),
        glue: Default::default(),
        fw_ws: Default::default(),
        rbuf: Default::default(),
        flip_mark: Vec::new(),
        pend: cap,
        cap,
        rebuilds: 40,
        first_recompute: carry.retries == 0 || fresh_pass::<R>(carry.retries),
        warm: start.is_some(),
        avg_touch: carry.avg_touch,
        last_cap: carry.last_cap,
        last_y: Vec::new(),
        scr_wv: Vec::new(),
        scr_rows: Vec::new(),
        scr_vals: Vec::new(),
        scr_ydbl: Vec::new(),
        scr_nbv: Vec::new(),
        scr_ub: Vec::new(),
        scr_binp: Vec::new(),
        rrows_exact: true,
        dup_entries: false,
        fac0: carry.fac,
        nref0: carry.nref0,
        nref_skip: 0,
        it0: carry.iters,
        updates: 0,
        upd_limit: None,
        pending_alpha: None,
        pair: (0.0, 0.0),
        status: None,
        drift_status: false,
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
        working_obj: None,
        obj: 0.0,
        // GUROBI: `ctx+0x128` holds -0.0 when the model has no objective constant; the term `[ctx+0x490]+0x10` of the recompute is null in every run seen
        // and taken as 0.
        obj_const: if model.offset == 0.0 { -0.0 } else { model.offset },
        flip_select: start.is_some(),
        work: if R::QUAD {
            crate::work::Work::dual_quad()
        } else {
            crate::work::Work::dual()
        },
        pivot_cols: None,
        nnz_a: p.cols[..n].iter().map(Vec::len).sum(),
        pstep: dual_pert_size,
        stall: 0,
        stall_snap: 0,
        stall_it: 0.0,
        prog_on: false,
        prog_acc: 0.0,
        best_obj: -1e100,
        pert_left: dual_pert_left,
        perturbed: false,
        pert_param: opts.perturb,
        obj_scale: opts.obj_scale,
        scaled: opts.scaled,
        pcost: p.cost.clone(),
        forced: Vec::new(),
        forced_ever: false,
        forced_stale: false,
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
        } else if R::QUAD && snapshot_restart && basis.at_upper[j] {
            UB
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
    // The public DStart (`reverse-engineer/PLAN.md`): d = c - A^T y for every variable, basic ones included; the basic variables of the start basis with |d| > 1e-7
    // form the forced-row stack (`0xeef09`, ascending), taken before the setup pivots of the free variables.
    // GUROBI: witnessed only in the first double pass with one y value for every row (`tests/dual_start.rs`: forced1, norow, AFIRO); DStart in the
    // scaled frame, max-sense models, non-uniform or mixed-sign y, and DStart in quad and restart passes have no witness (STATE.md, "Known unmodelled edges").
    let dstart_d: Option<Vec<R>> = match &opts.dstart {
        Some(y)
            if carry.retries == 0
                && carry.iters == 0
                && !R::QUAD
                && supplied.is_none()
                && y.len() == m =>
        {
            let mut d = vec![R::ZERO; n + m];
            for j in 0..n + m {
                let mut e: Vec<(usize, f64)> = p.fcols[j].clone();
                e.sort_by_key(|x| x.0);
                let mut acc = R::c(p.cost[j]);
                for &(i, a) in &e {
                    if R::c(y[i]).abs() > R::c(DROP_TOL) {
                        acc -= R::c(y[i]) * R::c(a);
                    }
                }
                d[j] = acc;
            }
            let mut stack: Vec<usize> = s
                .head
                .iter()
                .copied()
                .filter(|&j| d[j].abs() > R::c(1e-7))
                .collect();
            stack.sort_unstable();
            s.forced_ever = !stack.is_empty();
            s.forced = stack;
            Some(d)
        }
        _ => None,
    };
    s.work.decided(s.work.w);
    s.factor_work(&stc, rows.view(), 0.0);
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
        let rows: Vec<Vec<(usize, f64)>> = rows
            .into_iter()
            .enumerate()
            .map(|(i, row)| {
                row.into_iter()
                    .map(|j| {
                        let a = p.cols[j].iter().find(|e| e.0 == i).map_or(0.0, |e| e.1);
                        (j, a)
                    })
                    .collect()
            })
            .collect();
        s.rrows = FlatRows::from_rows(&rows, n);
        s.rrows_exact = false;
    }
    if let Some(dy) = &dstart_d {
        // Basic variables keep their d (the forced rows read it); the recompute below takes y from c_B minus the d of the stack entries.
        for j in 0..n + m {
            if s.st[j] >= 0 {
                s.d[j] = dy[j];
            }
        }
    }
    let dt = s.true_d();
    // The start pass runs the residual check of the post-refactor pass with the request 1 (`0xe5cda`; witness_primal_r31_double4: factorization count 9, residual 0.452).
    s.request = 1;
    s.dual_residual_check();
    s.work
        .add_from(0x12de50, 0xe53a7, s.mirror.work(), crate::work::FAC);
    let level_work = s.set_level_work(false);
    s.work.add(0xe7000, level_work, crate::work::CTX);
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
        if R::QUAD && snapshot_restart {
            s.d[j] = if s.st[j] == UB { -dt[j] } else { dt[j] };
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
    // A restarted quad setup with a supplied saved basis retains the saved
    // bound sides (`0x23a9b0`). Its request-1 recompute (`0x23a6d5`) shifts
    // negative reduced costs by only the random fraction of their magnitude;
    // they may consequently remain negative at the first pricing call.
    // Only a reduced cost beyond the dual tolerance is shifted: 4250's restart keeps a wrong-signed d of -3.6e-10 unshifted (8631 shifts -0.36).
    if R::QUAD && snapshot_restart {
        for j in 0..n + m {
            if (s.st[j] != LB && s.st[j] != UB) || s.d[j] >= R::c(-dual_tol()) {
                continue;
            }
            let magnitude = -s.d[j];
            let delta = R::c(dual_tol());
            let small = delta * R::c(0.0005) > magnitude
                || (delta * R::c(0.05) > magnitude
                    && delta * R::c(0.1) > R::c(p.up[j] - p.lo[j]) * magnitude);
            let factor = if small {
                s.raise_pert(1);
                s.rng_double()
            } else {
                s.raise_pert(2);
                ((s.rng_int() as i32 % 500 + 111) as f64) / 1000.0
            };
            let shift = magnitude * R::c(factor);
            s.d[j] += shift;
            if s.st[j] == UB {
                s.add_shift_q(j, -shift);
                s.obj_off += p.up[j] * shift.to_f64();
            } else {
                s.add_shift_q(j, shift);
                s.obj_off -= p.lo[j] * shift.to_f64();
            }
        }
    }
    // Not in the free phase, `s.mirror` is the setup factor with no update yet (only read by `true_d` and the residual check above).
    if !free_phase && s.ptype > 0 {
        if crate::dump::on("p68") || s.dup_entries {
            s.reset_pricing_ref(&stc, rows.view());
        } else {
            let w = pricing_weights(m, &stc, rows.view(), &s.mirror);
            s.set_pricing_weights(w);
        }
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
    let mut sq_buf: Vec<R> = Vec::new();
    for r in 0..if s.warm { 0 } else { m } {
        let out = s.mirror.btran(&[(r, R::ONE)], false, false);
        s.work
            .add_from(0x12de50, 0xc36f1, s.mirror.work(), crate::work::FAC);
        s.w[r] = if s.mirror.bt_dense() {
            rho_norm_dense_with(&out, m, &mut sq_buf)
        } else {
            rho_norm_with(&out, &mut sq_buf)
        };
        s.mirror.give_list(out);
    }
    if R::QUAD {
        // The quad pass keeps the perturbation flag and the stall counters of the double pass (same context), so its setup stall test can perturb (`0x239085`, c6007 quad pass 1).
        if let Some((a, b, c)) = SWITCH_STALL.with(|c| c.take()) {
            s.stall = a;
            s.stall_snap = b;
            s.stall_it = c;
        }
        if let Some(pert) = SWITCH_PERT.with(|c| c.take()) {
            s.pert = pert;
        }
    }
    s.stall_check(0)?;
    // The setup's post-refactor pass (`0xe6dcc`) sees the iterations of the free-column pivots.
    if !s.flip_select && 3.0 * s.iters as f64 > s.m as f64 {
        s.flip_select = true;
    }
    // The quad pass keeps the perturbation flag of the double pass (taken above) but not its cost shifts: d is recomputed from the true costs (seed 2240: the shifted
    // d of the double pass are 10^3..10^4 times the reference's at the first quad call, and so are the d of every quad refactorization).
    s.recompute_xb(0)?;
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
    s.finish_drift_status()?;
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
        s.factor.ftran(s.p, &mut col)?;
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
            if let Some(o) = rec.as_object_mut() {
                // The dual feasibility tolerance `ctx+0x388` in force at this call.
                o.insert("tol".to_string(), serde_json::json!(dual_tol()));
            }
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
                s.reclaim(pv);
                // The driver around the pricing (`0xf6745..0xf6768`, `0xf8b51..0xf8b77`) clears the worst rejected row (`ctx+0x430` = -1, `ctx+0x438` = 0)
                // each time it returns, so the "first rejection" is per iteration; the rejected count `ctx+0x420` stays (fx_midsing_1106 call 20).
                s.driver_return();
            }
            Step::Reject => {
                // The native iteration has non-pivot objective stores not yet modelled.
                s.working_obj = None;
                // A row dropped in the list-walk mode returns from the iteration function with the request it set (`0xc9577`): the post-iteration pass refactorizes.
                if s.request != 0 && s.status.is_none() {
                    s.post(None, false)?;
                    s.driver_return();
                }
            }
            Step::Mismatch => {
                s.working_obj = None;
                if s.status.is_none() {
                    s.post(None, false)?;
                }
                s.driver_return();
            }
            Step::Shifted => {
                s.working_obj = None;
                if s.request != 0 && s.status.is_none() {
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
                SWITCH_STALL
                    .with(|c| c.set(Some((s.stall, s.stall_snap, s.stall_it - s.iters as f64))));
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
                s.working_obj = None;
                if s.updates != 0 {
                    // No entering variable on a factor with updates (ctx+0x294 = 5, RECIPE-like bgdbg1 calls 29 and 31): refactorize and price again; the
                    // certificate is the row of the first such failure on a fresh factor.
                    s.request = 5;
                    s.post(None, false)?;
                    if s.status.is_none() {
                        continue;
                    }
                } else {
                    // GUROBI: the reference's combined status `INF_OR_UNBD` (4) is never produced; no witness (none in the dual searches,
                    // `evidence/witness-dual-search.md` item 14).
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
                && !s.no_art_exit
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
                && !s.no_art_exit
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
                s.sync_bounds();
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
                s.recompute_xb(0)?;
                s.finish_drift_status()?;
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
        popts.max_iterations = opts.max_iterations.saturating_sub(s.iters);
        popts.method = crate::simplex::Method::Primal;
        // The dual's teacher-forced factorization schedule does not apply to the primal, which decides its own (work based).
        popts.refactor_schedule = None;
        if R::QUAD && s.iters > 0 {
            // The x87 start selects method 0, 2, or 1 for ctx+0x35c equal
            // to 1, 2, or at least 3. The post-factor decision can also move
            // method 0 to 1 when the dual made iterations since its snapshot.
            let increments = !(s.pert != 0 && !s.rejected_end && !s.art_exit)
                && !(s.retries >= 5 && dropped.is_some() && s.pert == 0 && carry.after_pivot);
            let c35c_at_start = C35C.with(|c| c.get()) + i32::from(increments);
            // GUROBI: the decision's `ctx+0x178 <= 0` cost comparison and `ctx+0x4d8` are not modelled; no witness (PLAN-dual.md, "Remaining
            // extended-precision rules").
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
                pert: 0.0,
            });
        }
        let b = Basis {
            basic: s.head.clone(),
            at_upper: s.st.iter().map(|&t| t == UB).collect(),
        };
        crate::simplex::HANDOFF.with(|c| {
            c.set(Some(crate::simplex::Handoff {
                stale: Some(s.pair),
                // Shift and offset accounting in the primal progress writer is not yet recovered.
                working_obj: (R::QUAD && s.pert == 0 && s.obj_off == 0.0 && model.offset == 0.0)
                    .then_some(s.working_obj)
                    .flatten(),
                nref0: s.nref(),
                thr: s.lu_thr,
                retries: s.retries as usize,
                pert: s.pert >= 2 && !s.rejected_end && !s.art_exit,
                pair_tol: PAIR_TOL.with(|c| c.get()),
                quad: R::QUAD,
                dual_level: s.level,
                shifted_exit: !R::QUAD && s.rejected_end && s.pert != 0,
                iters: (s.it_base + s.iters) as f64,
                pricing_iters: (s.it_base + s.iters) as f64,
                snap: (s.it_base + s.stall_it as usize) as f64,
                c178: C178.with(|c| c.get()),
                preserve_ratio: false,
                // A shifted artificial-bound exit restores the repaired variable
                // before starting primal pricing (seed 6575, status -1 at entry).
                excl: if R::QUAD && s.pert == 0 {
                    s.repaired
                } else {
                    None
                },
            }))
        });
        // The in-place repair after the retry limit ends the dual with status 3, not -1004: `ctx+0x35c` stays (`standalone_restart_context`, 5959 dual-out 3).
        let repaired_end =
            R::QUAD && s.retries >= 5 && dropped.is_some() && s.pert == 0 && carry.after_pivot;
        if !(R::QUAD && s.pert != 0 && !s.rejected_end && !s.art_exit) && !repaired_end {
            C35C.with(|c| c.set(c.get() + 1));
        }
        // The status-3 end of the in-place repair is a restore of the controller (`ebx` of `0x1551ac`, `controller_restore_count`: 1 at iteration 16, 2 at 24, 3 at 122, 4 at 130
        // in 5959): from the third the controller alternates the methods, so the next primal's own restore switches to the dual.
        if repaired_end {
            PERT_RERUNS.with(|c| c.set(c.get() + 1));
        }
        crate::simplex::PRIMAL_SINGULAR.with(|c| c.set(None));
        crate::simplex::PRIMAL_TOL_ABORT.with(|c| c.set(false));
        crate::simplex::PRIMAL_START_SINGULAR.with(|c| c.set(false));
        crate::simplex::PRIMAL_FARKAS_SKIPPED.with(|c| c.set(0));
        crate::simplex::PRIMAL_FARKAS_MADE.with(|c| c.set(0));
        crate::simplex::PRIMAL_REPAIRED.with(|c| c.set(false));
        crate::simplex::PRIMAL_SHIFT_END.with(|c| c.set(false));
        if !R::QUAD {
            let nn = s.p.lo.len();
            let xd: Vec<f64> = (0..nn)
                .map(|j| {
                    if s.st[j] >= 0 {
                        s.xb[s.st[j] as usize].to_f64()
                    } else {
                        s.nonbasic_value(j)
                    }
                })
                .collect();
            DUAL_X_AT_HANDOFF.with(|c| *c.borrow_mut() = Some(xd));
        }
        if R::QUAD {
            let nn = s.p.lo.len();
            let xr: Vec<(f64, f64)> = (0..nn)
                .map(|j| {
                    let v = if s.st[j] >= 0 {
                        s.xb[s.st[j] as usize]
                    } else {
                        R::c(s.nonbasic_value(j))
                    };
                    let hi = v.to_f64();
                    (hi, (v - R::c(hi)).to_f64())
                })
                .collect();
            DUAL_XR_AT_HANDOFF.with(|c| *c.borrow_mut() = Some(xr));
        }
        let transfer = PrimalTransfer {
            basis: b,
            opts: popts,
            status,
            trace,
            route: PrimalRoute::Restored,
            control_status: (!(R::QUAD && s.pert != 0 && !s.rejected_end && !s.art_exit))
                .then_some(-1004),
            snapshot: HandoffSnapshot {
                iters: s.iters,
                refactorizations: s.refactorizations,
                lu_thr: s.lu_thr,
                avg_touch: s.avg_touch,
                last_cap: s.last_cap,
                inf_count: s.inf_count,
                retries: s.retries,
                it_base: s.it_base,
                nref: s.nref(),
                pert: s.pert,
                rejected_end: s.rejected_end,
                level: s.level,
                stall_it: s.stall_it,
                art_exit: s.art_exit,
            },
        };
        callback.finish();
        return Ok(DriverOutcome::ToPrimal(transfer));
    }
    if status == Status::Unbounded && s.art_exit {
        // The double dual's request-4 exit is a run that ends with -1004 too (`0xf8691`).
        C35C.with(|c| c.set(c.get() + 1));
        // A request-4 artificial-bound exit returns -1004 to the controller, which runs the primal from this basis. Its completed status is final.
        let mut popts = opts.clone();
        popts.max_iterations = opts.max_iterations.saturating_sub(s.iters);
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
                pert: 0.0,
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
                    working_obj: (s.obj_off == 0.0 && model.offset == 0.0)
                        .then_some(s.working_obj)
                        .flatten(),
                    nref0: s.nref(),
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
                    pricing_iters: (s.it_base + s.iters) as f64,
                    snap: (s.it_base + s.stall_it as usize) as f64,
                    c178: C178.with(|c| c.get()),
                    preserve_ratio: false,
                    excl: if s.pert == 0 { s.repaired } else { None },
                }))
            });
        }
        if !R::QUAD {
            crate::simplex::PRIMAL_START_THR.with(|c| c.set(Some(s.lu_thr)));
        }
        crate::simplex::PRIMAL_SINGULAR.with(|c| c.set(None));
        crate::simplex::PRIMAL_START_SINGULAR.with(|c| c.set(false));
        crate::simplex::PRIMAL_REPAIRED.with(|c| c.set(false));
        let transfer = PrimalTransfer {
            basis: from_dual,
            opts: popts,
            status,
            trace,
            control_status: Some(-1004),
            route: PrimalRoute::Artificial {
                start: start.cloned(),
            },
            snapshot: HandoffSnapshot {
                iters: s.iters,
                refactorizations: s.refactorizations,
                lu_thr: s.lu_thr,
                avg_touch: s.avg_touch,
                last_cap: s.last_cap,
                inf_count: s.inf_count,
                retries: s.retries,
                it_base: s.it_base,
                nref: s.nref(),
                pert: s.pert,
                rejected_end: s.rejected_end,
                level: s.level,
                stall_it: s.stall_it,
                art_exit: s.art_exit,
            },
        };
        callback.finish();
        return Ok(DriverOutcome::ToPrimal(transfer));
    }
    if status == Status::Unbounded {
        // Status 5 on the dual itself goes directly to the controller's final status store. Keep its raw artificial-bound ray and return the finite
        // point obtained by clearing the artificial bounds, without solving a second feasibility problem.
        s.zero_artificial = true;
        s.recompute_xb(0)?;
        s.finish_drift_status()?;
    }
    if status == Status::IterationLimit && s.level != 0 {
        s.zero_artificial = true;
        s.recompute_xb(0)?;
        s.finish_drift_status()?;
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
    let dtrue: Vec<f64> = s.true_d().iter().map(|v| v.to_f64()).collect();
    // The multipliers are those of the same reference btran as the reduced costs (`true_d`; long double in the quad pass), rounded once.
    let y: Vec<f64> = s.last_y.iter().map(|v| v.to_f64()).collect();
    let nm = model.num_cols();
    let xs: Vec<f64> = x[..nm].to_vec();
    let at_upper: Vec<bool> = s.st.iter().map(|&t| t == UB).collect();
    QUAD_XR.with(|c| *c.borrow_mut() = if R::QUAD { Some(xr) } else { None });
    PASS_END12.with(|c| c.set(s.status == Some(Status::NumericalFailure)));
    Ok(DriverOutcome::Finished(SolveResult {
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
    }))
}

// Resume the artificial-bound continuation after the dual callback has returned.
fn resume_artificial_transfer<R: Real>(
    model: &LpModel,
    opts: &SolveOptions,
    transfer: PrimalTransfer,
) -> Result<SolveResult, SolveError> {
    let PrimalTransfer {
        basis: from_dual,
        opts: popts,
        snapshot: s,
        route,
        ..
    } = transfer;
    let PrimalRoute::Artificial { start } = route else {
        unreachable!()
    };
    let PrimalOutcome { result, control } = dispatch_primal(model, &from_dual, &popts);
    crate::simplex::PRIMAL_START_THR.with(|c| c.set(None));
    if crate::dump::on("callback") {
        crate::dump::emit(
            "callback",
            serde_json::json!({"event":"primal_outcome", "result":if result.is_ok() { Some(0) } else { None }, "status":control.map(ControlTransfer::status), "error":result.is_err(), "total":crate::simplex::TOTAL_ITERS.with(|c| c.get())}),
        );
    }
    let mut prres = match result {
        Err(SolveError::SingularBasis) if s.art_exit => {
            SINGULAR_FROM_PRIMAL.with(|c| c.set(true));
            // A failed x87 handoff factorization returns -1002 to the
            // controller, which restores the saved snapshot (seed 8825).
            let Some((pit, pfac, pthr)) = crate::simplex::PRIMAL_SINGULAR.with(|c| c.take()) else {
                return Err(SolveError::SingularBasis);
            };
            if s.retries < 5 {
                C178.with(|c| c.set(c.get() + 1));
            }
            SINGULAR_BASIS.with(|c| *c.borrow_mut() = Some(from_dual.clone()));
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
        other => other.ok(),
    };
    // A progress exit after an artificial-bound handoff resumes dual from
    // the primal basis and increments the shared -1004 counter.
    let primal_stalled = R::QUAD && crate::simplex::PRIMAL_STALL.with(|c| c.take());
    // GUROBI: the reference installs the quad dual's x on the `primal_stalled` routes; rose does not, no witness (STATE.md, "Known unmodelled edges").
    if primal_stalled {
        C35C.with(|c| c.set(c.get() + 1));
        // The controller ends a -1004 cycle at status 12 once the
        // shared handoff counter exceeds 20 (0x15557d..0x15558c).
        if C35C.with(|c| c.get() > 20) {
            if let Some(mut r) = prres {
                r.status = Status::NumericalFailure;
                r.iterations += s.iters;
                r.refactorizations += s.refactorizations;
                return Ok(r);
            }
        }
    }
    // An artificial-bound dual exit can also hand off to a primal whose
    // start factorization drops a singular column. Restore that column and
    // retry the primal when its failed scan ends with the repair still out.
    if R::QUAD
        && s.art_exit
        && !primal_stalled
        && crate::simplex::PRIMAL_REPAIRED.with(|c| c.take())
    {
        if let Some(pr) = prres.as_mut().filter(|pr| pr.status == Status::Infeasible) {
            crate::simplex::HANDOFF.with(|c| {
                c.set(Some(crate::simplex::Handoff {
                    stale: None,
                    // Restored-variable objective carry is not yet recovered.
                    working_obj: None,
                    nref0: s.nref + pr.refactorizations,
                    thr: s.lu_thr,
                    retries: s.retries as usize + 1,
                    pert: false,
                    pair_tol: PAIR_TOL.with(|c| c.get()),
                    quad: true,
                    dual_level: s.level,
                    shifted_exit: false,
                    iters: (s.it_base + s.iters) as f64,
                    pricing_iters: (s.it_base + s.iters + pr.iterations) as f64,
                    snap: (s.it_base + s.stall_it as usize) as f64,
                    c178: C178.with(|c| c.get()),
                    preserve_ratio: true,
                    excl: None,
                }))
            });
            let mut retry_opts = popts.clone();
            retry_opts.max_iterations = popts.max_iterations.saturating_sub(pr.iterations);
            let mut retry = primal_callback(model, Some(&pr.basis), &retry_opts)?;
            retry.iterations += pr.iterations;
            retry.refactorizations += pr.refactorizations;
            *pr = retry;
        }
    }
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
                c.get() <= 1 || primal_stalled
            })
        {
            // A stalled x87 primal can resume quad repeatedly, retaining
            // the context and its clocks (5959). Other routes retain the
            // existing single-resume behavior pending their recovery.
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
                after_pivot: primal_stalled,
                it_base: s.it_base + s.iters + pit,
                nref0: s.nref + pfac.max(1),
            };
            let mut retry_opts = opts.clone();
            retry_opts.max_iterations = opts.max_iterations.saturating_sub(carry3.iters);
            let again = solve_pass::<R>(
                model,
                resume.as_ref().or(start.as_ref()),
                &retry_opts,
                &carry3,
                false,
            );
            NO_ART_EXIT.with(|c| c.set(false));
            let mut again = again?;
            again.iterations += carry3.iters;
            return Ok(again);
        }
        let again = solve_frame(model, resume.as_ref().or(start.as_ref()), opts);
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

// The callback has returned before this controller dispatch runs the primal.
fn resume_primal_transfer<R: Real>(
    model: &LpModel,
    opts: &SolveOptions,
    transfer: PrimalTransfer,
) -> Result<SolveResult, SolveError> {
    let PrimalTransfer {
        basis: b,
        opts: popts,
        snapshot: s,
        status,
        mut trace,
        control_status: _,
        route: _,
    } = transfer;
    let PrimalOutcome {
        result: primal_result,
        control,
    } = dispatch_primal(model, &b, &popts);
    // The inferred control status remains separate from the solve result
    // and error so the legacy continuation can distinguish all three.
    let control_status = control.map(ControlTransfer::status);
    if crate::dump::on("callback") {
        crate::dump::emit(
            "callback",
            serde_json::json!({"event":"primal_outcome", "result":if primal_result.is_ok() { Some(0) } else { None }, "status":control_status, "error":primal_result.is_err(), "total":crate::simplex::TOTAL_ITERS.with(|c| c.get())}),
        );
    }
    let mut r = match primal_result {
        Err(SolveError::SingularBasis) => {
            SINGULAR_FROM_PRIMAL.with(|c| c.set(true));
            // The primal's factorization has a column without a pivot (status -1002): the controller restarts the dual from the slack basis, like the dual's
            // own failed factorization, with the primal's iterations counted.
            let Some((pit, pfac, pthr)) = crate::simplex::PRIMAL_SINGULAR.with(|c| c.take()) else {
                return Err(SolveError::SingularBasis);
            };
            let tol_abort = crate::simplex::PRIMAL_TOL_ABORT.with(|c| c.take());
            SINGULAR_TOL_ABORT.with(|c| c.set(tol_abort));
            if s.retries < 5 && !tol_abort {
                // The primal LU singular handler increments the same counter
                // used by the next dual ratio test (0xf58f7 / 0x248f99).
                C178.with(|c| c.set(c.get() + 1));
            }
            if !tol_abort
                && !R::QUAD
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
                working_obj: None,
                nref0: s.nref + r.refactorizations,
                thr: s.lu_thr,
                retries: s.retries as usize + repair_count,
                pert: false,
                pair_tol: PAIR_TOL.with(|c| c.get()),
                quad: R::QUAD,
                dual_level: s.level,
                shifted_exit: false,
                iters: (s.it_base + s.iters) as f64,
                // The native context keeps its iteration counter across this retry;
                // `r.iterations` includes every completed primal run here.
                pricing_iters: (s.it_base + s.iters + r.iterations) as f64,
                snap: (s.it_base + s.stall_it as usize) as f64,
                c178: C178.with(|c| c.get()),
                preserve_ratio: R::QUAD,
                excl: None,
            }))
        });
        let b2 = r.basis.clone();
        let mut retry_opts = popts.clone();
        retry_opts.max_iterations = popts.max_iterations.saturating_sub(r.iterations);
        let r2 = primal_callback(model, Some(&b2), &retry_opts)?;
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
            nref0: s.nref + (r.refactorizations - s.refactorizations).max(1),
        };
        let mut retry_opts = opts.clone();
        retry_opts.max_iterations = opts.max_iterations.saturating_sub(carry2.iters);
        let mut r2 = solve_pass::<R>(model, Some(&r.basis), &retry_opts, &carry2, false)?;
        r2.iterations += carry2.iters;
        return Ok(r2);
    }
    // GUROBI: the reference installs the quad dual's x on this `PRIMAL_STALL` route; rose does not, no witness (STATE.md, "Known unmodelled edges").
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
            nref0: s.nref + l.max(1),
        };
        let mut retry_opts = opts.clone();
        retry_opts.max_iterations = opts.max_iterations.saturating_sub(carry2.iters);
        let mut r2 = solve_pass::<R>(model, Some(&r.basis), &retry_opts, &carry2, false)?;
        r2.iterations += carry2.iters;
        return Ok(r2);
    }
    if r.status == Status::NumericalFailure {
        // The primal's phase 1 ratio test finds no row on a fresh factor: the reference ends the solve with status 12 (fx_midsing_1031), the controller does not
        // restart the dual. ObjVal is `ctx+0x118`, still the objective of the dual's x at its last handoff.
        if !R::QUAD {
            if let Some(xd) = DUAL_X_AT_HANDOFF.with(|c| c.borrow_mut().take()) {
                let nm = model.num_cols();
                if xd.len() >= model.obj.len() && nm == model.obj.len() {
                    let mut acc = 0.0f64;
                    for (c, v) in model.obj.iter().zip(&xd) {
                        acc += *c * *v;
                    }
                    r.objective = acc + model.offset;
                    // GUROBI: the scaled-model path of this primal-ended status-12 objective (`0x144e6f`, read out in `unscaled_objective`) is
                    // unverified; no witness (STATE.md, "Known unmodelled edges").
                    if opts.scaled {
                        END12_X.with(|c| *c.borrow_mut() = Some(xd));
                    }
                    PASS_END12.with(|c| c.set(false));
                }
            }
        }
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
            // A rejected-row dual handoff can end with an ordinary primal
            // infeasibility result. The controller exits on status 3 when no
            // dropped variable remains; only a shifted primal exit reruns it.
            && (!s.rejected_end || shift_end)
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
        // The -1004 return of the primal ends the solve with status 12 once the shared counter exceeds 20 (`0x15557d` -> `0x15553e`; witness_primal_r25_quad20).
        if R::QUAD && C35C.with(|c| c.get() > 20) {
            if let Some(xr) = DUAL_XR_AT_HANDOFF.with(|c| c.borrow_mut().take()) {
                set_x87_result(xr);
            }
            r.status = Status::NumericalFailure;
            return Ok(r);
        }
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
                s.nref
                    + if R::QUAD {
                        l.max(1)
                    } else {
                        {
                            crate::simplex::PRIMAL_FARKAS_SKIPPED.with(|c| c.take());
                            l.saturating_sub(crate::simplex::PRIMAL_FARKAS_MADE.with(|c| c.take()))
                        }
                    }
            },
        };
        let mut retry_opts = opts.clone();
        retry_opts.max_iterations = opts.max_iterations.saturating_sub(carry2.iters);
        let mut r2 = match solve_pass::<R>(model, Some(&r.basis), &retry_opts, &carry2, false) {
            // A singular end of the rerun reaches the controller loop; a record made by the primal holds only the iterations since the last dual start, so the earlier ones are added (`ctx+0x20`).
            Err(SolveError::SingularBasis) => {
                SINGULAR.with(|c| {
                    if let Some(mut t) = c.take() {
                        if SINGULAR_FROM_PRIMAL.with(|c| c.get()) {
                            t.0 += carry2.iters;
                        }
                        c.set(Some(t));
                    }
                });
                return Err(SolveError::SingularBasis);
            }
            other => other?,
        };
        r2.iterations += carry2.iters;
        return Ok(r2);
    }
    return Ok(r);
}

#[cfg(test)]
mod captured_quad_partition {
    use super::accumulate;
    use crate::f80::F80;
    use serde_json::Value;

    fn value(v: &Value) -> F80 {
        let part =
            |i| u64::from_str_radix(v[i].as_str().unwrap().trim_start_matches("0x"), 16).unwrap();
        F80::from_repr(part(0), part(1) as u16)
    }

    #[test]
    fn captured_add_and_subtract() {
        // qwide8360, call 70: reference partition operands and x87 results.
        let cases: Value = serde_json::from_str(include_str!(
            "../tests/fixtures/quad_partition_capture.json"
        ))
        .unwrap();
        for case in cases.as_array().unwrap() {
            let c: Vec<_> = case["con"].as_array().unwrap().iter().map(value).collect();
            assert_eq!(
                accumulate(value(&case["start"]), &c, 0, case["sub"].as_bool().unwrap()).repr(),
                value(&case["end"]).repr()
            );
        }
    }
}

#[cfg(test)]
mod captured_bh305_ratio {
    use super::*;
    use serde_json::Value;

    #[test]
    fn fresh_factor_artificial_bound_retry() {
        // Captured reference call 1, not a reconstruction of the unknown original model.
        let case: Value =
            serde_json::from_str(include_str!("../tests/fixtures/bh305_ratio_capture.json"))
                .unwrap();
        let input = &case["input"];
        let n = input["n"].as_u64().unwrap() as usize;
        let m = input["m"].as_u64().unwrap() as usize;
        let counts = input["colcnt"].as_array().unwrap();
        // Only column lengths are read by ratio_test; the placeholder matrix is used
        // to construct an otherwise unused factor and mirror.
        let mut starts = vec![0];
        let mut rows = Vec::new();
        for count in counts {
            for i in 0..count.as_u64().unwrap() as usize {
                rows.push(i);
            }
            starts.push(rows.len());
        }
        let model = LpModel {
            name: "captured-ratio-input".into(),
            sense: Sense::Minimize,
            offset: 0.0,
            obj: vec![0.0; n],
            col_lo: vec![0.0; n],
            col_up: vec![INF; n],
            row_lo: vec![-INF; m],
            row_up: vec![0.0; m],
            row_range: vec![],
            row_kind: vec![],
            range_order: vec![],
            range_entries: vec![],
            range_dup: false,
            col_start: starts,
            val: vec![1.0; rows.len()],
            row_idx: rows,
            col_names: vec![],
            row_names: vec![],
        };
        let mut p = build(&model);
        let idx = input["row"]["idx"].as_array().unwrap();
        let vals = input["row"]["val"].as_array().unwrap();
        let row: Vec<(usize, f64)> = idx
            .iter()
            .zip(vals)
            .map(|(j, a)| (j.as_u64().unwrap() as usize, a.as_f64().unwrap()))
            .collect();
        let mut st = vec![FIXED; n + m];
        let mut d = vec![0.0; n + m];
        for (i, &(j, _)) in row.iter().enumerate() {
            st[j] = input["st"][i].as_i64().unwrap() as i32;
            d[j] = input["d"][i].as_f64().unwrap();
            p.lo[j] = input["lb"][i].as_f64().unwrap();
            p.up[j] = input["ub"][i].as_f64().unwrap();
        }
        let basic: Vec<_> = (n..n + m).collect();
        let (order, structure, cols) = reference_order::<f64>(&p, &basic, 1);
        let mirror = Mirror::fresh(m, n, &structure, &cols);
        let basis_cols: Vec<_> = order.iter().map(|&j| &p.cols[j]).collect();
        let factor = Factor::new(m, &basis_cols).unwrap();
        let cap = 20;
        let (blo, bup) = basic_bounds::<f64>(&p, &order);
        let mut s: Dual<f64> = Dual {
            p: &p,
            n,
            m,
            head: order,
            blo,
            bup,
            st,
            xb: vec![0.0; m],
            d: vec![0.0; n + m],
            w: vec![1.0; m],
            mirror,
            factor: AuxiliaryFactor::new(factor),
            rrows: FlatRows::new(m, n),
            absorb_request: false,
            fresh_absorb: false,
            refactorizations: 1,
            lu_thr: 0.0078125,
            big_pricing: false,
            pmode: 1,
            n_infeas: 0,
            inf_count: 0,
            retries: 0,
            repaired: None,
            snapshot_excluded: false,
            bad_rows: Vec::new(),
            old_rej: 0,
            rejected_end: false,
            art_exit: false,
            no_art_exit: false,
            worst: 0.0,
            worst_row: -1,
            retry: false,
            retry_end: 0,
            retry_tiny: 0,
            stale_count: false,
            last_pivoted: false,
            it_base: 0,
            dflag: Vec::new(),
            switch_enabled: false,
            automatic_quad: false,
            level: 0,
            request: 0,
            iters: 0,
            prio: vec![0.0; m],
            ptype: 0,
            p68: Vec::new(),
            p70: Vec::new(),
            p60: 0.0,
            p58: 0.0,
            p50d: 0.0,
            p58d: 0.0,
            dv_f: 0.0,
            dv_g: 0.0,
            dv_k: 0,
            dv_c: 0,
            flags: vec![2; m],
            list: Vec::new(),
            thr: 0.0,
            xb_residual: 0.0,
            value_mode: false,
            cursor: 0,
            walk_len: 0,
                slack_nb: vec![false; m],
            prow_ws: Default::default(),
            glue: Default::default(),
            fw_ws: Default::default(),
            rbuf: Default::default(),
            flip_mark: Vec::new(),
            pend: cap,
            cap,
            rebuilds: 40,
            first_recompute: true,
            warm: false,
            avg_touch: 0.0,
            last_cap: 0,
            last_y: Vec::new(),
            scr_wv: Vec::new(),
            scr_rows: Vec::new(),
            scr_vals: Vec::new(),
            scr_ydbl: Vec::new(),
            scr_nbv: Vec::new(),
            scr_ub: Vec::new(),
            scr_binp: Vec::new(),
            rrows_exact: true,
            dup_entries: false,
            fac0: 0,
            nref0: 0,
            nref_skip: 0,
            it0: 0,
            updates: 0,
            upd_limit: None,
            pending_alpha: None,
            pair: (0.0, 0.0),
            status: None,
            drift_status: false,
            nm: n,
            ray: None,
            zero_artificial: false,
            cost: p.cost.clone(),
            cshift: vec![0.0; p.cost.len()],
            cost_lo: vec![0.0; p.cost.len()],
            shift_r: vec![0.0; p.cost.len()],
            rng: 101,
            pert: 0,
            obj_off: 0.0,
            working_obj: None,
            obj: 0.0,
            obj_const: -0.0,
            flip_select: false,
            work: crate::work::Work::dual(),
            pivot_cols: None,
            nnz_a: p.cols[..n].iter().map(Vec::len).sum(),
            pstep: 0.0,
            stall: 0,
            stall_snap: 0,
            stall_it: 0.0,
            prog_on: false,
            prog_acc: 0.0,
            best_obj: -1e100,
            pert_left: 1_000_000_000,
            perturbed: false,
            pert_param: 0.0,
            obj_scale: 1.0,
            scaled: false,
            pcost: p.cost.clone(),
            forced: Vec::new(),
            forced_ever: false,
            forced_stale: false,
        };
        s.d = d;
        s.level = input["lvl"].as_u64().unwrap() as usize;
        s.iters = input["it"].as_f64().unwrap() as usize;
        let rc = &case["candidate_collection"];
        assert_eq!(rc["n"].as_u64().unwrap(), 7);
        assert_eq!(rc["pos"], serde_json::json!([0, 2, 3, 4, 5, 6, 9]));
        assert_eq!(
            rc["x12"].as_f64().unwrap().to_bits(),
            (-2.0 * input["xB"].as_f64().unwrap()).to_bits()
        );
        for pos in [0, 2, 4, 9] {
            let j = row[pos].0;
            assert_eq!(s.st[j], UB);
            assert!(p.up[j] >= BIG);
        }
        // The capture does not record the factor's update count: exercise a
        // fresh factor explicitly (updates = 0).
        let leave = input["var"].as_u64().unwrap() as usize;
        // Stale contents of the retained buffers must not matter.
        s.rbuf.lst.extend([7, 3, 9]);
        s.rbuf.kept.extend([1, 2]);
        s.rbuf.con.extend([f64::NAN; 4]);
        s.rbuf.ratio_arr.extend([f64::NAN; 2]);
        let mut row = row;
        let result = unsafe { s.ratio_test(
            &mut row,
            false,
            input["esi"].as_u64().unwrap() as usize,
            leave,
            false,
            input["xB"].as_f64().unwrap(),
        ) };
        let q = result.q.expect("reference selected a candidate");
        assert_eq!(row[q].0 as i64, case["expected"]["q"].as_i64().unwrap());
        assert_eq!(
            result.flips.len(),
            case["expected"]["nflip"].as_u64().unwrap() as usize
        );
        assert_eq!(
            row[q].1.to_bits(),
            case["expected"]["a308"].as_f64().unwrap().to_bits()
        );
        assert!(result.reject.is_none());
    }
}

#[cfg(test)]
thread_local! {
    static QUAD_OBJ_TEST: std::cell::RefCell<Option<Vec<f64>>> = const { std::cell::RefCell::new(None) };
}

#[cfg(test)]
fn record_quad_obj(value: f64) {
    QUAD_OBJ_TEST.with(|c| {
        if let Some(values) = c.borrow_mut().as_mut() {
            values.push(value);
        }
    });
}

#[cfg(test)]
#[test]
fn quad_working_objective_captured_stores() {
    // Reference stores at 0x217d6c and 0x21dd67 through the first quad dual.
    let expected = [
        0.0,
        0.25,
        0.33347262103511405,
        0.5165277613548308,
        0.5165277613548308,
        0.8565894421544799,
        1.4295000733793137,
        3.5351275804953075,
        3.5357673398085137,
        3.5955687150231457,
        3.595568715023162,
        4.517244820478645,
        4.517246676550257,
        349.14162441750926,
        349.14813281525466,
        6352.899623308976,
        6352.899623308975,
        95526943.50521018,
        221699539.5595171,
        221699549.95117602,
    ];
    let model = crate::mps::parse(include_str!("../tests/fixtures/fx_quad_obj_5371.mps")).unwrap();
    let mut opts = crate::SolveOptions::dual_release();
    opts.max_iterations = 200;
    QUAD_OBJ_TEST.with(|c| *c.borrow_mut() = Some(Vec::new()));
    crate::solve(&model, None, &opts).unwrap();
    let actual = QUAD_OBJ_TEST.with(|c| c.borrow_mut().take().unwrap());
    assert!(actual.len() >= expected.len());
    for (index, (actual, expected)) in actual.iter().zip(expected).enumerate() {
        assert_eq!(actual.to_bits(), f64::to_bits(expected), "store {index}");
    }
}

#[cfg(test)]
thread_local! {
    static LOWGATE_TEST: std::cell::RefCell<Option<Vec<(&'static str, usize, u32, f64, Vec<usize>)>>> = const { std::cell::RefCell::new(None) };
    static CEILING_TEST: std::cell::RefCell<Option<Vec<(&'static str, usize, u32, Vec<usize>)>>> = const { std::cell::RefCell::new(None) };
    static CEILING_NO_ACCEPT: std::cell::Cell<bool> = const { std::cell::Cell::new(false) };
}

#[cfg(test)]
#[test]
fn drift_ceiling_without_acceptance_reaches_restart() {
    // Native noaccept capture changes request 4 to 0 at c1fe3.
    let model = crate::mps::parse(include_str!("../tests/fixtures/fx_quad_obj_5371.mps")).unwrap();
    let mut opts = crate::SolveOptions::dual_release();
    opts.max_iterations = 200;
    CEILING_TEST.with(|c| *c.borrow_mut() = Some(Vec::new()));
    CEILING_NO_ACCEPT.with(|c| c.set(true));
    let _ = crate::solve(&model, None, &opts);
    CEILING_NO_ACCEPT.with(|c| c.set(false));
    let events = CEILING_TEST.with(|c| c.borrow_mut().take().unwrap());
    let head = vec![23, 31, 32, 6, 1, 4, 5, 10, 21, 18, 20];
    assert_eq!(
        events,
        vec![("pending", 16, 4, head.clone()), ("restart", 16, 0, head)]
    );
}

#[cfg(test)]
#[test]
fn drift_ceiling_acceptance_overwrites_pending_status() {
    // Counterfactual native 5371 capture: count 2 and tolerance 1e-3 at
    // c177e; c17e7 writes -1002, then c21f6 overwrites it with 2.
    let model = crate::mps::parse(include_str!("../tests/fixtures/fx_quad_obj_5371.mps")).unwrap();
    let mut opts = crate::SolveOptions::dual_release();
    opts.max_iterations = 200;
    CEILING_TEST.with(|c| *c.borrow_mut() = Some(Vec::new()));
    let result = crate::solve(&model, None, &opts);
    let events = CEILING_TEST.with(|c| c.borrow_mut().take().unwrap());
    let head = vec![23, 31, 32, 6, 1, 4, 5, 10, 21, 18, 20];
    assert_eq!(
        events,
        vec![("pending", 16, 4, head.clone()), ("accepted", 16, 4, head)]
    );
    assert!(result.is_ok());
}

#[cfg(test)]
#[test]
fn lowgate_counterfactual_refactor_replay() {
    let model = crate::mps::parse(include_str!("../tests/fixtures/bh50028.mps")).unwrap();
    let mut opts = crate::SolveOptions::dual();
    opts.max_iterations = 30;
    opts.refactor_schedule = Some(vec![3, 8, 12]);
    LOWGATE_TEST.with(|c| *c.borrow_mut() = Some(Vec::new()));
    let result = crate::solve(&model, None, &opts).unwrap();
    let records = LOWGATE_TEST.with(|c| c.borrow_mut().take().unwrap());
    // Reference gdb_dual_lowgate_head_replay, CF=1: the intervention changes
    // only the LU threshold at the first eligible gate, after 16 matching calls.
    let head = vec![20, 19, 18, 26, 2, 6, 13, 17, 12];
    let expected = [
        ("postfactor", 15, 4, 0.5),
        ("eligible", 15, 4, 0.5),
        ("after", 15, 4, 0.125),
        ("postfactor", 15, 4, 0.125),
        ("eligible", 15, 4, 0.5),
        ("optimal", 15, 0, 0.5),
    ];
    assert_eq!(records.len(), expected.len());
    for ((tag, it, req, thr, actual_head), expected) in records.iter().zip(expected) {
        assert_eq!((*tag, *it, *req, *thr), expected);
        assert_eq!(*actual_head, head);
    }
    assert_eq!(result.status, Status::Optimal);
}

#[cfg(test)]
mod perturb_size_lifecycle_tests {
    use super::*;

    #[test]
    fn captured_sizes_and_signed_entry_rules() {
        reset_pert_sizes(0.0002);
        C35C.with(|c| c.set(1));
        assert_eq!(
            enter_pert_size(&DUAL_PERT_SIZE, false).to_bits(),
            0.0002_f64.to_bits()
        );
        C35C.with(|c| c.set(2));
        assert_eq!(enter_primal_pert_size().to_bits(), 1e-5_f64.to_bits());
        assert_eq!(
            DUAL_PERT_SIZE.with(|c| c.get()).to_bits(),
            0.0002_f64.to_bits()
        );
        assert_eq!(
            enter_pert_size(&DUAL_PERT_SIZE, false).to_bits(),
            1e-5_f64.to_bits()
        );
        C35C.with(|c| c.set(3));
        assert_eq!(
            enter_primal_pert_size().to_bits(),
            5.000000000000001e-7_f64.to_bits()
        );
        C35C.with(|c| c.set(4));
        assert_eq!(
            enter_pert_size(&DUAL_PERT_SIZE, false).to_bits(),
            5.000000000000001e-7_f64.to_bits()
        );
        assert_eq!(enter_primal_pert_size().to_bits(), 1e-7_f64.to_bits());
        C35C.with(|c| c.set(-4));
        assert_eq!(
            enter_pert_size(&DUAL_PERT_SIZE, false).to_bits(),
            5.000000000000001e-7_f64.to_bits()
        );
        assert_eq!(C35C.with(|c| c.get()), -4);
        enter_primal_pert_size();
        assert_eq!(C35C.with(|c| c.get()), 4);
        C35C.with(|c| c.set(-1));
        assert_eq!(enter_dual_pert_size().1, 1);
        assert_eq!(C35C.with(|c| c.get()), 0);
        C35C.with(|c| c.set(-4));
        assert_eq!(enter_dual_pert_size().1, 1_000_000_000);
        assert_eq!(C35C.with(|c| c.get()), -4);
        reset_pert_sizes(0.0002);
        assert_eq!(PRIMAL_PERT_SIZE.with(|c| c.get()), 0.0002);
        assert_eq!(DUAL_PERT_SIZE.with(|c| c.get()), 0.0002);
    }
}

#[cfg(test)]
#[test]
fn method_handoffs_share_iteration_limit() {
    for (source, limits) in [
        (
            include_str!("../tests/fixtures/fx_quad_obj_5371.mps"),
            &[17, 32][..],
        ),
        (
            include_str!("../tests/fixtures/primal_5959.mps"),
            &[12, 32][..],
        ),
    ] {
        let model = crate::mps::parse(source).unwrap();
        for &limit in limits {
            let mut opts = SolveOptions::dual_release();
            opts.max_iterations = limit;
            let result = crate::solve(&model, None, &opts).unwrap();
            assert!(
                result.iterations <= limit,
                "limit {limit}, got {}",
                result.iterations
            );
            let clock = crate::simplex::TOTAL_ITERS.with(|c| c.get());
            assert!(clock <= limit as f64, "limit {limit}, clock {clock}");
        }
    }
}

#[cfg(test)]
#[test]
fn iteration_limit_keeps_the_last_allowed_pivot() {
    let model = crate::mps::parse(include_str!("../tests/fixtures/afiro.mps")).unwrap();
    let opts = SolveOptions::dual_release();
    let complete = crate::solve(&model, None, &opts).unwrap();
    assert!(complete.iterations > 0);
    let mut bounded = opts;
    bounded.max_iterations = complete.iterations;
    let result = crate::solve(&model, None, &bounded).unwrap();
    assert_eq!(result.iterations, complete.iterations);
    let mut bounded_basic = result.basis.basic.clone();
    let mut complete_basic = complete.basis.basic.clone();
    bounded_basic.sort_unstable();
    complete_basic.sort_unstable();
    assert_eq!(bounded_basic, complete_basic);
    assert_eq!(result.basis.at_upper, complete.basis.at_upper);
}

#[test]
fn callback_carry_survives_primal_and_live_override() {
    let opts = SolveOptions::dual_release();
    set_primal_tol(1e-6);
    set_dual_tol(1e-6);
    callback_reset();
    C35C.with(|c| c.set(8));
    let mut dual = CallbackGuard::enter(true, &opts);
    assert_eq!(primal_tol().to_bits(), (1e-6_f64 * 0.2).to_bits());
    set_primal_tol(4e-7);
    dual.finish();
    let expected_carry = 4e-7_f64 / 0.2;
    {
        let _primal = CallbackGuard::enter(false, &opts);
        set_primal_tol(1e-5);
    }
    let _dual = CallbackGuard::enter(true, &opts);
    assert_eq!(primal_tol().to_bits(), expected_carry.to_bits());
}

#[cfg(test)]
mod captured_primal_drift_gate {
    use super::*;
    use crate::f80::F80;

    #[test]
    fn reference_5371_dual_and_ceiling_actions() {
        // Native double event at total 16: residual and bound violation.
        assert!(primal_drift_excess(
            0.010148395773596355_f64,
            0.002116363008820097_f64
        ));
        assert!(!primal_drift_excess(0.001_f64, 0.1_f64));
        assert_eq!(
            primal_drift_action(2, 1e-6, false, true),
            (1e-6, false, false)
        );
        assert_eq!(
            primal_drift_action(3, 1e-6, false, true),
            (1e-6 * 10.0, true, false)
        );
        assert_eq!(
            primal_drift_action(3, 1e-3, false, true),
            (1e-3, false, true)
        );
        assert_eq!(
            primal_drift_action(3, 1e-3, false, false),
            (1e-3, false, false)
        );
        assert_eq!(
            primal_drift_action(3, 1e-3, true, true),
            (1e-3, false, false)
        );
        // A sub-double-ULP distinction must survive until the x87 comparison.
        let one = F80::from_f64(1.0);
        let epsilon = F80::from_f64(crate::real::pow2(-60));
        assert!(primal_drift_excess(
            (one + epsilon) / F80::from_f64(100.0),
            one
        ));
    }
}

#[cfg(test)]
mod pivot_row_diff {
    use super::*;

    /// The previous implementation, kept as the reference of the differential test.
    #[allow(clippy::too_many_arguments)]
    fn pivot_row_ref<R: Real>(
        rrows: &[Vec<(usize, f64)>],
        slack_nb: &[bool],
        st: &[i32],
        cols: &[Vec<(usize, f64)>],
        n: usize,
        m: usize,
        rho: &[(usize, R)],
        leave: usize,
        above: bool,
    ) -> (Vec<(usize, R)>, usize, usize, usize) {
        let est: usize = rho.iter().map(|e| rrows[e.0].len()).sum();
        let mut acc = vec![R::ZERO; n + m];
        let mut touched = vec![false; n + m];
        let mut order: Vec<usize> = Vec::new();
        let sparse = rho.len() == 1 || !((n as f64) < 1.8 * est as f64);
        for &(i, rv) in rho {
            for &(j, a) in &rrows[i] {
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
        let distinct = order.len();
        if rho.len() != 1 || rho[0].1.abs() < R::ONE {
            order.retain(|&j| acc[j].abs() > R::c(DROP_TOL));
        }
        let nst = order.len();
        for &(i, rv) in rho {
            let j = n + i;
            if slack_nb[i] || j == leave {
                acc[j] = rv * R::c(cols[j][0].1);
                order.push(j);
            }
        }
        let sdir = if above { -1.0 } else { 1.0 };
        let out = order
            .into_iter()
            .map(|j| {
                let s = if st[j] == UB { -1.0 } else { 1.0 };
                (j, acc[j] * R::c(s * sdir))
            })
            .collect();
        (out, nst, est, distinct)
    }

    struct Rng(u64);
    impl Rng {
        fn next(&mut self) -> u64 {
            self.0 ^= self.0 << 13;
            self.0 ^= self.0 >> 7;
            self.0 ^= self.0 << 17;
            self.0
        }
        fn below(&mut self, k: usize) -> usize {
            (self.next() % k as u64) as usize
        }
    }

    fn value(r: &mut Rng) -> f64 {
        match r.below(9) {
            0 => 0.0,
            1 => -0.0,
            2 => 1.0,
            3 => -1.0,
            4 => 1e-13,
            5 => 3e-14,
            6 => 5e-5,
            7 => (r.below(2001) as f64 - 1000.0) / 8.0,
            _ => (r.below(2001) as f64 - 1000.0) / 1000.0,
        }
    }

    /// `distinct` is booked (and compared) only on the sparse route; the columnwise route returns 0 (`row_product_work` does not read it there).
    fn same<R: Real>(
        a: &(Vec<(usize, R)>, usize, usize, usize),
        b: &(Vec<(usize, R)>, usize, usize, usize),
        n: usize,
        nrho: usize,
    ) -> bool {
        let sparse = nrho == 1 || !((n as f64) < 1.8 * a.2 as f64);
        a.1 == b.1
            && a.2 == b.2
            && (if sparse { a.3 == b.3 } else { b.3 == 0 })
            && a.0.len() == b.0.len()
            && a.0.iter().zip(&b.0).all(|(x, y)| x.0 == y.0 && format!("{:?}", x.1).as_bytes() == format!("{:?}", y.1).as_bytes()
                && x.1.to_f64().to_bits() == y.1.to_f64().to_bits())
    }

    fn run<R: Real>() {
        let mut r = Rng(0x9e3779b97f4a7c15);
        for case in 0..300 {
            let n = 1 + r.below(40);
            let m = 1 + r.below(20);
            let mut rrows: Vec<Vec<(usize, f64)>> = vec![Vec::new(); m];
            for j in 0..n {
                for i in 0..m {
                    if r.below(4) == 0 {
                        rrows[i].push((j, value(&mut r)));
                    }
                }
            }
            // Per-row entry order is shuffled-ish to exercise first-touch order.
            for row in rrows.iter_mut() {
                if row.len() > 1 && r.below(2) == 0 {
                    let k = row.len();
                    row.swap(0, k - 1);
                }
            }
            let slack_nb: Vec<bool> = (0..m).map(|_| r.below(2) == 0).collect();
            let st: Vec<i32> = (0..n + m).map(|_| if r.below(2) == 0 { UB } else { LB }).collect();
            let cols: Vec<Vec<(usize, f64)>> = (0..n + m)
                .map(|j| if j < n { Vec::new() } else { vec![(j - n, if r.below(2) == 0 { 1.0 } else { -1.0 })] })
                .collect();
            let mut flat = FlatRows::from_rows(&rrows, n);
            let mut ws = PivotWs::<R>::default();
            for call in 0..8 {
                // Pivot-like edits of both copies (position, in-place replace, swap_remove, push, also duplicates and emptied rows): same in-row order.
                for _ in 0..r.below(6) * (call.min(1)) {
                    let i = r.below(m);
                    let nj = r.below(n);
                    match r.below(3) {
                        0 if !rrows[i].is_empty() => {
                            let k = r.below(rrows[i].len());
                            let j = rrows[i][k].0;
                            assert_eq!(flat.position(i, j), rrows[i].iter().position(|e| e.0 == j));
                            let kk = flat.position(i, j).unwrap();
                            let a = value(&mut r);
                            rrows[i][kk] = (nj, a);
                            flat.set(i, kk, nj, a);
                        }
                        1 if !rrows[i].is_empty() => {
                            let k = r.below(rrows[i].len());
                            rrows[i].swap_remove(k);
                            flat.swap_remove(i, k);
                        }
                        _ => {
                            let a = value(&mut r);
                            rrows[i].push((nj, a));
                            flat.push(i, nj, a);
                        }
                    }
                }
                if r.below(8) == 0 {
                    flat.compact();
                }
                for i in 0..m {
                    let (ix, vl) = flat.row(i);
                    assert!(
                        ix.len() == rrows[i].len()
                            && rrows[i].iter().enumerate().all(|(k, e)| e.0 == ix[k] as usize && e.1.to_bits() == vl[k].to_bits()),
                        "case {case} call {call} row {i}"
                    );
                }
                let k = if r.below(4) == 0 { 1 } else { 1 + r.below(m) };
                let mut rows: Vec<usize> = (0..m).collect();
                for t in 0..m {
                    let u = r.below(m);
                    rows.swap(t, u);
                }
                let rho: Vec<(usize, R)> = rows[..k.min(m)]
                    .iter()
                    .map(|&i| (i, R::c(if r.below(3) == 0 { value(&mut r) * 4.0 } else { value(&mut r) })))
                    .collect();
                // The rows of `rho` are distinct (a BTRAN result has one entry per row), see `pivot_row_ws`.
                let leave = n + r.below(m);
                let above = r.below(2) == 0;
                let want = pivot_row_ref(&rrows, &slack_nb, &st, &cols, n, m, &rho, leave, above);
                let got = pivot_row_ws(&mut ws, &flat, &slack_nb, &st, &cols, n, m, &rho, leave, above);
                assert!(same(&want, &got, n, rho.len()), "case {case} call {call}");
                // The first-match position of the leaving variable and the duplicate flag.
                assert_eq!(ws.leave_pos, got.0.iter().position(|e| e.0 == leave), "case {case} call {call}");
                let dup = (0..got.0.len()).any(|x| got.0[..x].iter().any(|e| e.0 == got.0[x].0));
                assert_eq!(ws.slack_dup, dup, "case {case} call {call}");
                assert!(ws.acc.iter().all(|x| x.to_f64().to_bits() == 0), "acc not restored");
                assert!(ws.bits.iter().all(|&w| w == 0));
            }
        }
    }

    /// Rows of exactly the lengths 0..=5, 7, 8 (and 6, 9, 12), columns shared between rows, `rho` over all or some rows in a random order,
    /// the leaving slack inside and outside `rho`, on the columnwise (small n) and the row-wise (large n) route.
    fn run_lengths<R: Real>() {
        let mut r = Rng(0x2545f4914f6cdd1d);
        let lens = [0usize, 1, 2, 3, 4, 5, 7, 8, 6, 9, 12, 4, 0, 8];
        let m = lens.len();
        for n in [12usize, 16, 200] {
            for case in 0..60 {
                let mut rrows: Vec<Vec<(usize, f64)>> = vec![Vec::new(); m];
                for (i, &l) in lens.iter().enumerate() {
                    let mut cs: Vec<usize> = (0..n).collect();
                    for t in 0..l.min(n) {
                        let u = t + r.below(n - t);
                        cs.swap(t, u);
                        rrows[i].push((cs[t], value(&mut r)));
                    }
                }
                let slack_nb: Vec<bool> = (0..m).map(|_| r.below(2) == 0).collect();
                let st: Vec<i32> = (0..n + m).map(|_| if r.below(2) == 0 { UB } else { LB }).collect();
                let cols: Vec<Vec<(usize, f64)>> = (0..n + m)
                    .map(|j| if j < n { Vec::new() } else { vec![(j - n, if r.below(2) == 0 { 1.0 } else { -1.0 })] })
                    .collect();
                let flat = FlatRows::from_rows(&rrows, n);
                let mut ws = PivotWs::<R>::default();
                for call in 0..12 {
                    let mut rows: Vec<usize> = (0..m).collect();
                    for t in 0..m {
                        let u = r.below(m);
                        rows.swap(t, u);
                    }
                    let k = match call % 4 {
                        0 => m,
                        1 => 1,
                        _ => 1 + r.below(m),
                    };
                    let rho: Vec<(usize, R)> = rows[..k].iter().map(|&i| (i, R::c(value(&mut r)))).collect();
                    // Inside `rho` on the even calls, outside it (when some row is left out) on the odd ones.
                    let leave = if call % 2 == 0 { n + rows[r.below(k)] } else { n + rows[(k + r.below(m - k + 1)).min(m - 1)] };
                    let above = r.below(2) == 0;
                    let want = pivot_row_ref(&rrows, &slack_nb, &st, &cols, n, m, &rho, leave, above);
                    let got = pivot_row_ws(&mut ws, &flat, &slack_nb, &st, &cols, n, m, &rho, leave, above);
                    assert!(same(&want, &got, n, rho.len()), "n {n} case {case} call {call}");
                    assert_eq!(ws.leave_pos, got.0.iter().position(|e| e.0 == leave), "n {n} case {case} call {call}");
                    assert!(!ws.slack_dup);
                    assert!(ws.acc.iter().all(|x| x.to_f64().to_bits() == 0) && ws.bits.iter().all(|&w| w == 0));
                }
            }
        }
    }

    /// The columnwise route with many cancellations: values from a small grid so that sums hit exactly 0.0 (also from -0.0 and -x + x products),
    /// tiny values around DROP_TOL, and one or two huge ones; columns shared by many rows.
    fn run_cancel<R: Real>() {
        let mut r = Rng(0x1234_5678_9abc_def1);
        let grid = [0.0, -0.0, 1.0, -1.0, 0.5, -0.5, 1e-13, -1e-13, 1.0000000000000002e-13, 9.999999999999e-14, 3e-14, -3e-14, 1e300, -1e300, 1e-300];
        let mut sparse_calls = 0;
        let mut col_calls = 0;
        for case in 0..400 {
            let n = 1 + r.below(70);
            let m = 2 + r.below(24);
            let mut rrows: Vec<Vec<(usize, f64)>> = vec![Vec::new(); m];
            for row in rrows.iter_mut() {
                let l = r.below(n.min(12) + 1);
                let mut cs: Vec<usize> = (0..n).collect();
                for t in 0..l {
                    let u = t + r.below(n - t);
                    cs.swap(t, u);
                    row.push((cs[t], grid[r.below(grid.len())]));
                }
            }
            let slack_nb: Vec<bool> = (0..m).map(|_| r.below(2) == 0).collect();
            let st: Vec<i32> = (0..n + m).map(|_| if r.below(2) == 0 { UB } else { LB }).collect();
            let cols: Vec<Vec<(usize, f64)>> = (0..n + m)
                .map(|j| if j < n { Vec::new() } else { vec![(j - n, if r.below(2) == 0 { 1.0 } else { -1.0 })] })
                .collect();
            let flat = FlatRows::from_rows(&rrows, n);
            let mut ws = PivotWs::<R>::default();
            for call in 0..10 {
                let mut rows: Vec<usize> = (0..m).collect();
                for t in 0..m {
                    let u = r.below(m);
                    rows.swap(t, u);
                }
                let k = if call % 5 == 0 { 1 } else { 2 + r.below(m - 1) };
                let rho: Vec<(usize, R)> = rows[..k].iter().map(|&i| (i, R::c(grid[r.below(grid.len())]))).collect();
                let leave = n + r.below(m);
                let above = r.below(2) == 0;
                let want = pivot_row_ref(&rrows, &slack_nb, &st, &cols, n, m, &rho, leave, above);
                let got = pivot_row_ws(&mut ws, &flat, &slack_nb, &st, &cols, n, m, &rho, leave, above);
                assert!(same(&want, &got, n, rho.len()), "cancel case {case} call {call}");
                if rho.len() == 1 || !((n as f64) < 1.8 * want.2 as f64) {
                    sparse_calls += 1;
                } else {
                    col_calls += 1;
                }
                assert_eq!(ws.leave_pos, got.0.iter().position(|e| e.0 == leave), "cancel case {case} call {call}");
                assert!(ws.acc.iter().all(|x| x.to_f64().to_bits() == 0) && ws.bits.iter().all(|&w| w == 0), "workspace not restored");
            }
        }
        assert!(sparse_calls > 100 && col_calls > 500, "{sparse_calls} sparse, {col_calls} columnwise");
    }

    #[test]
    fn matches_reference_cancel_f64() {
        run_cancel::<f64>();
    }

    #[test]
    fn matches_reference_cancel_f80() {
        run_cancel::<crate::f80::F80>();
    }

    #[test]
    fn matches_reference_lengths_f64() {
        run_lengths::<f64>();
    }

    #[test]
    fn matches_reference_lengths_f80() {
        run_lengths::<crate::f80::F80>();
    }

    #[test]
    fn matches_reference_f64() {
        run::<f64>();
    }

    #[test]
    fn matches_reference_f80() {
        run::<crate::f80::F80>();
    }
}

#[cfg(test)]
mod d_step_diff {
    use super::*;

    /// The previous d update: a `contains` test per entry.
    fn d_step_ref<R: Real>(d: &mut [R], row: &[(usize, R)], flips: &[usize], theta: R) {
        for &(j, a) in row {
            let a = if flips.contains(&j) { -a } else { a };
            d[j] = d[j] - theta * a;
        }
    }

    fn run<R: Real>() {
        let mut seed = 0x2545_f491_4f6c_dd1du64;
        let mut next = move || {
            seed ^= seed << 13;
            seed ^= seed >> 7;
            seed ^= seed << 17;
            seed
        };
        let values = [0.0, -0.0, 1.0, -1.0, 0.1, 1e-14, -3.5, 1e300, -1e300, f64::MIN_POSITIVE, 0.3, f64::NAN, f64::INFINITY, f64::NEG_INFINITY, 1e-10, -1e-10];
        let mut mark = Vec::new();
        for _ in 0..2000 {
            let len = 1 + next() as usize % 12;
            let d0: Vec<R> = (0..len).map(|_| R::c(values[next() as usize % values.len()])).collect();
            let rl = next() as usize % 20;
            // Repeated variables in the row are allowed: the update then applies once per entry.
            let row: Vec<(usize, R)> =
                (0..rl).map(|_| (next() as usize % len, R::c(values[next() as usize % values.len()]))).collect();
            let nf = if next() % 2 == 0 { 0 } else { next() as usize % 6 };
            let flips: Vec<usize> = (0..nf).map(|_| next() as usize % len).collect();
            let theta = R::c(values[next() as usize % values.len()]);
            let (mut want, mut got) = (d0.clone(), d0.clone());
            d_step_ref(&mut want, &row, &flips, theta);
            unsafe { d_row_step(&mut got, &row, &flips, theta, &mut mark) };
            assert!(want.iter().zip(&got).all(|(x, y)| x.to_f64().to_bits() == y.to_f64().to_bits()));
            assert!(mark.iter().all(|&b| !b), "marks not cleared");
        }
    }

    #[test]
    fn matches_reference_f64() {
        run::<f64>();
    }

    #[test]
    fn matches_reference_f80() {
        run::<crate::f80::F80>();
    }
}

#[cfg(test)]
mod row_loops_diff {
    use super::*;

    /// The candidate loops before the unchecked rewrite.
    #[allow(clippy::too_many_arguments)]
    fn ratio_candidates_ref<R: Real>(
        lo: &[f64],
        up: &[f64],
        st: &[i32],
        d: &[R],
        row: &[(usize, R)],
        esi: usize,
        delta: R,
        xmm7: R,
        xmm12: R,
        variant: bool,
        ratio_arr: &mut Vec<R>,
        con: &mut Vec<R>,
        lst: &mut Vec<usize>,
    ) -> Cand<R> {
        let cnt = row.len();
        let mut bound = R::c(-1e100);
        let (mut best, mut best_i, mut bestc) = (R::c(-1e100), -1i64, R::ZERO);
        let (mut xm13, mut xm9) = (R::ZERO, R::ZERO);
        ratio_arr.clear();
        con.clear();
        lst.clear();
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
            if v < R::c(-ALPHA_TOL) && (st[j] == LB || st[j] == UB) {
                let dj = d[j];
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
        for i in 0..esi.min(cnt) {
            let (j, v) = row[i];
            if v < R::c(-ALPHA_TOL) {
                let dj = d[j];
                let ratio = dj / v;
                if bound > ratio {
                    xm9 = minr(v, xm9);
                    continue;
                }
                xm13 = minr(v, xm13);
                let c;
                if variant && st[j] == LB && up[j] >= BIG {
                    if dj >= R::ZERO {
                        bound = maxr(PIVOT_DELTA(delta, v, ratio), bound);
                    }
                    c = xmm12;
                } else {
                    let c0 = R::c(lo[j] - up[j]) * v;
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
        Cand { bound, best_i, bestc, xm13, xm9, edge: false }
    }

    fn same<R: Real>(a: R, b: R) -> bool {
        a.to_f64().to_bits() == b.to_f64().to_bits() && format!("{a:?}") == format!("{b:?}")
    }

    fn run<R: Real>() {
        let mut seed = 0x9e37_79b9_7f4a_7c15u64;
        let mut next = move || {
            seed ^= seed << 13;
            seed ^= seed >> 7;
            seed ^= seed << 17;
            seed
        };
        let tol = ALPHA_TOL;
        // Values at and around the tolerances, signed zeros, infinities, NaN.
        let vals = [
            0.0, -0.0, tol, -tol, -tol * (1.0 + 1e-15), -tol * (1.0 - 1e-15), 1e-10 + 1e-25, -1.0, -2.5, -3e-9, 1.0, 0.5, -1e-300, -1e300, f64::NAN,
            f64::INFINITY, f64::NEG_INFINITY, 1e-5, -7.25,
        ];
        let sts = [LB, UB, FREE, FIXED, 0, 3];
        let (mut nc, mut nq) = (0usize, 0usize);
        for it in 0..20000 {
            let nv = 1 + next() as usize % 12;
            let mut st: Vec<i32> = (0..nv).map(|_| sts[next() as usize % sts.len()]).collect();
            if it % 7 == 0 {
                st.iter_mut().for_each(|s| *s = LB);
            }
            let d: Vec<R> = (0..nv).map(|_| R::c(vals[next() as usize % vals.len()])).collect();
            let lo: Vec<f64> = (0..nv).map(|_| [0.0, -1.0, -BIG, 2.0, -0.0][next() as usize % 5]).collect();
            let up: Vec<f64> = (0..nv).map(|_| [1.0, BIG, 0.5, 1e30, 3.0][next() as usize % 5]).collect();
            let cnt = next() as usize % 20;
            // Repeated variables are allowed in a row.
            let row: Vec<(usize, R)> = (0..cnt).map(|_| (next() as usize % nv, R::c(vals[next() as usize % vals.len()]))).collect();
            let esi = (next() as usize % (cnt + 3)).min(cnt + 2);
            let delta = R::c(1e-7);
            let xmm12 = R::c([-2.0, -0.0, 0.0, -1e-3, f64::NAN][next() as usize % 5]);
            let xmm7 = R::c([-1.0, 0.0, 1e-6, 5.0][next() as usize % 4]);
            let variant = next() % 2 == 0;
            let (mut a, mut b, mut c) = (vec![R::c(f64::NAN); 2], vec![R::ZERO; 3], vec![7usize]);
            let (mut a2, mut b2, mut c2) = (vec![], vec![R::c(5.0)], vec![]);
            let want = ratio_candidates_ref(&lo, &up, &st, &d, &row, esi, delta, xmm7, xmm12, variant, &mut a, &mut b, &mut c);
            let mut row2 = row.clone();
            let got = unsafe { ratio_candidates::<R, false>(&lo, &up, &st, &d, &mut row2, esi, delta, xmm7, xmm12, variant, &mut a2, &mut b2, &mut c2) };
            assert!(row2.iter().zip(&row).all(|(x, y)| x.0 == y.0 && same(x.1, y.1)), "the plain loops must not write the row");
            // The fused zeroing: the same lists, and the row after the loops (plus the old pass when an entry sits exactly at the tolerance) is the old pass's.
            let mut row3 = row.clone();
            let (mut a3, mut b3, mut c3) = (vec![], vec![], vec![]);
            let fused = unsafe { ratio_candidates::<R, true>(&lo, &up, &st, &d, &mut row3, esi, delta, xmm7, xmm12, variant, &mut a3, &mut b3, &mut c3) };
            assert!(same(fused.bound, got.bound) && same(fused.bestc, got.bestc) && same(fused.xm13, got.xm13) && same(fused.xm9, got.xm9));
            assert_eq!((fused.best_i, &c3), (got.best_i, &c2));
            assert!(a3.iter().zip(&a2).all(|(x, y)| same(*x, *y)) && b3.iter().zip(&b2).all(|(x, y)| same(*x, *y)));
            let mut old = row.clone();
            for e in old.iter_mut() {
                if e.1.abs() <= R::c(ALPHA_TOL) {
                    e.1 = R::ZERO;
                }
            }
            let has_edge = row.iter().any(|e| e.1.abs() == R::c(ALPHA_TOL));
            assert_eq!(fused.edge, has_edge);
            // Before the old pass the row differs from the raw one only by zeroed entries below the tolerance (the entries read by `framework_rhs`).
            for (x, y) in row3.iter().zip(&row) {
                if y.1.abs() < R::c(ALPHA_TOL) {
                    assert!(same(x.1, R::ZERO) && format!("{:?}", x.1) == format!("{:?}", R::ZERO));
                } else {
                    assert!(same(x.1, y.1));
                }
            }
            if fused.edge {
                for e in row3.iter_mut() {
                    if e.1.abs() <= R::c(ALPHA_TOL) {
                        e.1 = R::ZERO;
                    }
                }
            }
            assert!(row3.iter().zip(&old).all(|(x, y)| x.0 == y.0 && same(x.1, y.1)), "fused zeroing differs from the pass");
            assert!(same(want.bound, got.bound) && same(want.bestc, got.bestc) && same(want.xm13, got.xm13) && same(want.xm9, got.xm9));
            assert_eq!(want.best_i, got.best_i);
            assert_eq!(c, c2);
            assert!(a.len() == a2.len() && a.iter().zip(&a2).all(|(x, y)| same(*x, *y)));
            assert!(b.len() == b2.len() && b.iter().zip(&b2).all(|(x, y)| same(*x, *y)));
            nc += c.len();
            nq += usize::from(!c.is_empty());
        }
        assert!(nc > 5000 && nq > 3000, "the random rows should produce candidates ({nc}, {nq})");
    }

    #[test]
    fn candidates_match_reference_f64() {
        run::<f64>();
    }

    #[test]
    fn candidates_match_reference_f80() {
        run::<crate::f80::F80>();
    }
}

#[cfg(test)]
mod true_d_column_order {
    use super::*;

    /// `true_d` sums over `p.cols[j]`; the old code sorted a clone of `p.fcols[j]` by row. Both must be the same list bit for bit,
    /// ties (repeated rows) and signed zeros included, so the accumulation order is unchanged.
    #[test]
    fn sorted_clone_of_fcols_is_cols() {
        let mut seed = 0x9e37_79b9_7f4a_7c15u64;
        let mut next = move || {
            seed ^= seed << 13;
            seed ^= seed >> 7;
            seed ^= seed << 17;
            seed
        };
        let values = [0.0, -0.0, 1.0, -1.0, 0.1, 1e-14, -3.5, f64::MIN_POSITIVE];
        for round in 0..200 {
            let m = 1 + next() as usize % 9;
            let nc = 1 + next() as usize % 12;
            let mut col_start = vec![0];
            let (mut row_idx, mut val) = (Vec::new(), Vec::new());
            for _ in 0..nc {
                for _ in 0..next() % 8 {
                    row_idx.push(next() as usize % m);
                    val.push(values[next() as usize % values.len()]);
                }
                col_start.push(row_idx.len());
            }
            let kinds = [(-INF, 1.0), (1.0, INF), (1.0, 1.0), (-1.0, 2.0), (-INF, INF)];
            let rows: Vec<(f64, f64)> = (0..m).map(|_| kinds[next() as usize % kinds.len()]).collect();
            let model = LpModel {
                name: "order".into(),
                sense: if round % 2 == 0 { Sense::Minimize } else { Sense::Maximize },
                offset: 0.0,
                obj: vec![1.0; nc],
                col_lo: vec![0.0; nc],
                col_up: vec![INF; nc],
                row_lo: rows.iter().map(|r| r.0).collect(),
                row_up: rows.iter().map(|r| r.1).collect(),
                row_range: vec![],
                row_kind: vec![],
                range_order: vec![],
                range_entries: vec![],
                range_dup: false,
                col_start,
                val,
                row_idx,
                col_names: vec![],
                row_names: vec![],
            };
            let p = build(&model);
            let y: Vec<f64> = (0..p.m).map(|i| 0.5 + i as f64 * 0.37).collect();
            for j in 0..p.n + p.m {
                let mut e: Vec<(usize, f64)> = p.fcols[j].clone();
                e.sort_by_key(|x| x.0);
                assert_eq!(e.len(), p.cols[j].len());
                for (a, b) in e.iter().zip(&p.cols[j]) {
                    assert_eq!((a.0, a.1.to_bits()), (b.0, b.1.to_bits()), "round {round} col {j}");
                }
                let (mut old, mut new) = (0.1f64, 0.1f64);
                for &(i, a) in &e {
                    old -= y[i] * a;
                }
                for &(i, a) in &p.cols[j] {
                    new -= y[i] * a;
                }
                assert_eq!(old.to_bits(), new.to_bits());
            }
        }
    }
}

/// `true_d_rowwise` against the version it replaced: per column, a stably sorted clone of `fcols[j]`, subtracting `y[i] a` in row order for
/// `|y[i]| > DROP_TOL`. Columns with repeated rows, slack and range columns, signed zeros in the values, costs and `y`.
#[cfg(test)]
mod true_d_differential {
    use super::*;

    fn old_true_d(p: &Internal, st: &[i32], cost: &[f64], cost_lo: &[f64], y: &[f64]) -> Vec<f64> {
        let mut d = vec![0.0; p.n + p.m];
        for j in 0..p.n + p.m {
            if st[j] < 0 {
                let mut e = p.fcols[j].clone();
                e.sort_by_key(|x| x.0);
                let mut acc = cost[j] + cost_lo[j];
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

    #[test]
    fn rowwise_matches_per_column_clone_and_sort() {
        let mut seed = 0x1b87_3593_cc9e_2d51u64;
        let mut next = move || {
            seed ^= seed << 13;
            seed ^= seed >> 7;
            seed ^= seed << 17;
            seed
        };
        let values = [0.0, -0.0, 1.0, -1.0, 0.1, 1e-14, -3.5, f64::MIN_POSITIVE, 1e16];
        let ys = [0.0, -0.0, 1e-13, -1e-13, 2e-13, 0.37, -1.25, 1e8];
        let costs = [0.0, -0.0, 1.0, -2.5, 1e-9];
        let (mut dup_cols, mut range_rows, mut compared) = (0, 0, 0);
        for round in 0..300 {
            let m = 1 + next() as usize % 9;
            let nc = 1 + next() as usize % 12;
            let mut col_start = vec![0];
            let (mut row_idx, mut val) = (Vec::new(), Vec::new());
            for _ in 0..nc {
                let first = row_idx.len();
                for _ in 0..next() % 8 {
                    row_idx.push(next() as usize % m);
                    val.push(values[next() as usize % values.len()]);
                }
                let mut rs = row_idx[first..].to_vec();
                rs.sort_unstable();
                rs.dedup();
                dup_cols += usize::from(rs.len() < row_idx.len() - first);
                col_start.push(row_idx.len());
            }
            let kinds = [(-INF, 1.0), (1.0, INF), (1.0, 1.0), (-1.0, 2.0), (-INF, INF), (0.5, 4.0)];
            let rows: Vec<(f64, f64)> = (0..m).map(|_| kinds[next() as usize % kinds.len()]).collect();
            range_rows += rows.iter().filter(|r| r.0 > -INF && r.1 < INF && r.0 != r.1).count();
            let model = LpModel {
                name: "true-d".into(),
                sense: if round % 2 == 0 { Sense::Minimize } else { Sense::Maximize },
                offset: 0.0,
                obj: (0..nc).map(|_| costs[next() as usize % costs.len()]).collect(),
                col_lo: vec![0.0; nc],
                col_up: vec![INF; nc],
                row_lo: rows.iter().map(|r| r.0).collect(),
                row_up: rows.iter().map(|r| r.1).collect(),
                row_range: vec![],
                row_kind: vec![],
                range_order: vec![],
                range_entries: vec![],
                range_dup: false,
                col_start,
                val,
                row_idx,
                col_names: vec![],
                row_names: vec![],
            };
            let p = build(&model);
            let nv = p.n + p.m;
            let st: Vec<i32> = (0..nv).map(|_| if next() % 3 == 0 { 0 } else { -1 }).collect();
            let cost: Vec<f64> = (0..nv).map(|_| costs[next() as usize % costs.len()]).collect();
            let cost_lo: Vec<f64> = (0..nv).map(|_| costs[next() as usize % costs.len()]).collect();
            let y: Vec<f64> = (0..p.m).map(|_| ys[next() as usize % ys.len()]).collect();
            let want = old_true_d(&p, &st, &cost, &cost_lo, &y);
            let got = true_d_rowwise::<f64>(&p, &st, &cost, &cost_lo, &y);
            assert_eq!(got.len(), want.len());
            for (j, (g, w)) in got.iter().zip(&want).enumerate() {
                assert_eq!(g.to_bits(), w.to_bits(), "round {round} column {j}");
            }
            compared += nv;
        }
        assert!(dup_cols > 0 && range_rows > 0 && compared > 0, "dup {dup_cols} range {range_rows}");
    }
}

#[cfg(test)]
mod reset_pricing_flat {
    use super::{cap_pricing_weights, pricing_weights, pricing_weights_lists, Mirror};

    fn bits(v: &[f64]) -> Vec<u64> {
        v.iter().map(|x| x.to_bits()).collect()
    }

    /// The weights from the factor's lists equal the list form bit for bit on random factor structures.
    #[test]
    fn flat_weights_match_the_list_form() {
        let mut seed = 0x2545f4914f6cdd1du64;
        let mut rnd = move |n: usize| {
            seed ^= seed << 13;
            seed ^= seed >> 7;
            seed ^= seed << 17;
            (seed % n as u64) as usize
        };
        let mut checked = 0;
        for case in 0..400 {
            let m = 1 + rnd(40);
            let extra = rnd(4 * m + 1);
            let mut cols: Vec<(bool, Vec<usize>, Vec<f64>)> = Vec::new();
            for c in 0..m {
                let mut rows = vec![c];
                for _ in 0..rnd(1 + extra / m.max(1) + 1) {
                    let r = rnd(m);
                    if !rows.contains(&r) {
                        rows.push(r);
                    }
                }
                rows.sort_unstable();
                let vals = rows.iter().map(|&r| if r == c { 4.0 } else { 0.5 + rnd(7) as f64 }).collect();
                cols.push((false, rows, vals));
            }
            // Permute the columns so the pivot order is not the identity.
            for i in (1..m).rev() {
                cols.swap(i, rnd(i + 1));
            }
            let st = crate::lu_order::factor_structure(m, &cols);
            if st.pivots.len() != m {
                continue;
            }
            let f = Mirror::fresh(m, m, &st, &cols);
            let flat = pricing_weights(m, &st, (&cols).into(), &f);
            let (_, _, lists) = pricing_weights_lists(m, &st, (&cols).into());
            assert_eq!(bits(&flat), bits(&lists), "case {case}, m {m}");
            checked += 1;
        }
        assert!(checked > 200);
    }

    /// The cap, including the type-2 and the over-total branches, against the original inline code.
    #[test]
    fn cap_branches() {
        fn old(ptype: u8, m: usize, w: &mut Vec<f64>) -> u8 {
            let mut ptype = ptype;
            let n = m as f64;
            let total: f64 = w.iter().sum();
            let cap = 0.1 * n + 100.0;
            if ptype == 2 || total > n * cap {
                ptype = 2;
                w.iter_mut().for_each(|x| *x = cap);
            } else {
                for x in w.iter_mut() {
                    let v = *x + 100.0;
                    *x = if cap < v { cap } else { v };
                }
            }
            ptype
        }
        for (ptype, m, scale) in [(1u8, 10usize, 1.0), (2, 10, 1.0), (1, 10, 1e3), (1, 200, 0.1), (1, 200, 1e6), (2, 200, 1e6)] {
            let w0: Vec<f64> = (0..m).map(|i| 1.0 + i as f64 * scale).collect();
            let (mut a, mut b) = (w0.clone(), w0);
            let (pa, pb) = (cap_pricing_weights(ptype, m, &mut a), old(ptype, m, &mut b));
            assert_eq!((pa, bits(&a)), (pb, bits(&b)));
        }
        // Both type outcomes occur among the cases above.
        let mut w = vec![1e9; 10];
        assert_eq!(cap_pricing_weights(1, 10, &mut w), 2);
        let mut w = vec![1.0; 10];
        assert_eq!(cap_pricing_weights(1, 10, &mut w), 1);
    }
}

#[cfg(test)]
mod infeas_tests {
    use super::*;

    /// The former `infeas_real` body, as reference.
    fn old_infeas_real<R: Real>(x: R, lo: f64, up: f64) -> R {
        if lo > -BIG && x < R::c(lo) {
            R::c(lo) - x
        } else if up < BIG && x > R::c(up) {
            x - R::c(up)
        } else {
            R::ZERO
        }
    }

    fn same<R: Real>(a: R, b: R) -> bool {
        let (a, b) = (a.to_f64(), b.to_f64());
        a.to_bits() == b.to_bits() || (a.is_nan() && b.is_nan())
    }

    fn check<R: Real>() {
        let vals = [
            0.0, -0.0, 1.0, -1.0, 1e-9, -1e-9, 0.5, 2.0, 1e29, -1e29, 1e30, -1e30, 1e31, -1e31, f64::INFINITY, f64::NEG_INFINITY, f64::NAN,
        ];
        let mut s = 0x9e3779b97f4a7c15u64;
        let mut next = || {
            s ^= s << 13;
            s ^= s >> 7;
            s ^= s << 17;
            s
        };
        let mut rnd = |next: &mut dyn FnMut() -> u64| {
            let u = next();
            if u & 1 == 0 {
                vals[(u >> 8) as usize % vals.len()]
            } else {
                ((next() >> 11) as f64 / (1u64 << 53) as f64 - 0.5) * 4.0
            }
        };
        for _ in 0..200_000 {
            let (x, lo, up) = (rnd(&mut next), rnd(&mut next), rnd(&mut next));
            // ties: x equal to a bound
            let x = match next() % 5 {
                0 => lo,
                1 => up,
                _ => x,
            };
            let (a, b) = (infeas_of::<R>(R::c(x), lo, up), old_infeas_real::<R>(R::c(x), lo, up));
            assert!(same(a, b), "x {x} lo {lo} up {up}");
            let (blo, bup) = basic_bound_pair::<R>(lo, up);
            let c = infeas_cached::<R>(R::c(x), blo, bup);
            assert!(same(c, b), "cached: x {x} lo {lo} up {up}");
        }
    }

    /// The former two-pass rebuild (`fill_infeas` + `build_walk_list_from`) on slices.
    fn old_rebuild<R: Real>(
        xb: &[R],
        blo: &[R],
        bup: &[R],
        ptol: f64,
        flags: &mut [u8],
        key: &dyn Fn(usize, R) -> f64,
    ) -> (Vec<usize>, Vec<f64>, usize, usize) {
        let infs: Vec<R> = xb.iter().zip(blo).zip(bup).map(|((&x, &l), &u)| infeas_cached::<R>(x, l, u)).collect();
        let (mut n_inf, mut n_fresh) = (0usize, 0usize);
        let (mut rows, mut keys) = (Vec::new(), Vec::new());
        for (r, &fr) in infs.iter().enumerate() {
            if fr.to_f64() > ptol {
                n_inf += 1;
                n_fresh += (flags[r] != 0) as usize;
                flags[r] = 0;
                rows.push(r);
                keys.push(key(r, fr));
            } else {
                flags[r] = 2;
            }
        }
        (rows, keys, n_inf, n_fresh)
    }

    fn rebuild_check<R: Real>() {
        let vals = [
            0.0, -0.0, 1.0, -1.0, 1e-9, 0.5, 2.0, 1e30, -1e30, f64::INFINITY, f64::NEG_INFINITY, f64::NAN,
        ];
        let mut s = 0x2545f4914f6cdd1du64;
        let mut next = move || {
            s ^= s << 13;
            s ^= s >> 7;
            s ^= s << 17;
            s
        };
        for round in 0..3000 {
            let m = match round % 4 {
                0 => 0,
                1 => 1 + (next() % 5) as usize,
                _ => (next() % 300) as usize,
            };
            let gen = |next: &mut dyn FnMut() -> u64| -> f64 {
                let u = next();
                if u % 3 == 0 {
                    vals[(u >> 8) as usize % vals.len()]
                } else {
                    // few distinct values: ties between keys
                    ((u >> 16) % 7) as f64 * 0.25 - 0.5
                }
            };
            let all_bad = round % 11 == 0;
            let xb: Vec<R> = (0..m).map(|_| R::c(gen(&mut next))).collect();
            let (blo, bup): (Vec<R>, Vec<R>) = (0..m)
                .map(|_| {
                    let (l, u) = if all_bad { (1.0, -1.0) } else { (gen(&mut next), gen(&mut next)) };
                    basic_bound_pair::<R>(l, u)
                })
                .unzip();
            let w: Vec<f64> = (0..m).map(|_| gen(&mut next)).collect();
            let ptol = [1e-6, 0.0, -1.0, 0.25, 1e-9][round % 5];
            let flags0: Vec<u8> = (0..m + 3).map(|_| (next() % 6) as u8).collect();
            let key = |r: usize, fr: R| (-(fr * fr / R::c(w[r]))).to_f64();
            let (mut fa, mut fb) = (flags0.clone(), flags0.clone());
            let (rows_a, keys_a, ni_a, nf_a) = old_rebuild::<R>(&xb, &blo, &bup, ptol, &mut fa, &key);
            let (mut rows_b, mut keys_b) = (vec![9usize; 3], vec![f64::NAN; 2]);
            let (ni_b, nf_b) = rebuild_rows::<R>(&xb, &blo, &bup, ptol, &mut fb, &mut rows_b, &mut keys_b, key);
            assert_eq!((ni_a, nf_a), (ni_b, nf_b));
            assert_eq!(fa, fb);
            assert_eq!(rows_a, rows_b);
            assert!(keys_a.iter().zip(&keys_b).all(|(a, b)| a.to_bits() == b.to_bits()));
            assert_eq!(keys_a.len(), keys_b.len());
        }
    }

    #[test]
    fn walk_rebuild_matches_two_pass_f64() {
        rebuild_check::<f64>();
    }

    #[test]
    fn walk_rebuild_matches_two_pass_quad() {
        rebuild_check::<crate::f80::F80>();
    }

    #[test]
    fn infeas_of_matches_reference_f64() {
        check::<f64>();
    }

    #[test]
    fn infeas_of_matches_reference_quad() {
        check::<crate::f80::F80>();
    }
}
