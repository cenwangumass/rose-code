//! Primal simplex from a supplied feasible basis.
//!
//! Rules follow the recovered reference behaviour (see `reverse-engineer/specs/primal-simplex.md`):
//! Devex pricing over a fixed reference framework, Harris two-pass ratio test with bound flips, and
//! reduced-cost updates from the pivot row. Internal form: `A x + s = rhs`, columns `0..n` structural then
//! one slack per row (`+e_i` for `=`/`<=`, `-e_i` for `>=`).

use crate::basis::Basis;
use crate::lu::{Column, Factor};
use crate::model::{LpModel, Sense, INF};

/// Algorithm behind `solve`.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Default)]
pub enum Method {
    /// The primal simplex (the default).
    #[default]
    Primal,
    /// The dual simplex (`crate::dual`). It reads `SolveOptions::refactor_schedule` as the iteration counts after which the reference's
    /// work-based trigger refactorizes (the refactorizations of the no-leaving-row handler are modelled).
    Dual,
    /// The restarted reflected Halpern PDHG (`crate::pdhg`), a first-order method: the result is approximate (tolerances in `SolveOptions::pdhg`),
    /// has no basis (`SolveResult::basis` is the all-slack placeholder; the start basis is ignored), no Farkas certificate and no ray.
    Pdhg,
    /// The primal-dual interior point method of the reference's `Method=2`, `Crossover=0` (`crate::barrier`, `crate::chol`), the controlled profile
    /// (`BarCorrectors=1`, no scaling, `BarOrder=0`), bit for bit on the models of `reverse-engineer/specs/barrier*.md`. The point is the barrier's own (no
    /// crossover): `SolveResult::basis` is the all-slack placeholder, there is no Farkas certificate and no ray. Features outside the recovered LP case
    /// (maximization, free or lower-unbounded columns, ranged or free rows, dense columns, partitioned orderings) return `SolveError::Unsupported`.
    Barrier,
}

#[derive(Clone, Debug)]
pub struct SolveOptions {
    pub method: Method,
    /// Options of `Method::Pdhg`.
    pub pdhg: crate::pdhg::PdhgOptions,
    pub feasibility_tol: f64,
    pub optimality_tol: f64,
    pub pivot_tol: f64,
    pub max_iterations: usize,
    /// Refactor after this many product-form updates.
    pub max_updates: usize,
    pub record_trace: bool,
    /// Reproduce the reference's Devex weights exactly (weights outside the framework may fall below 1). The reference's
    /// anti-degeneracy perturbation is not implemented; no stall was seen on the benchmark corpus. On by default; `SolveOptions::plain()` turns it off.
    pub faithful_weights: bool,
    /// Renumber basis positions after each refactorization the way the reference does (only the singleton part of the
    /// order is recovered; the nucleus is ordered by remaining entries and previous position). Affects tie-breaking only.
    pub faithful_positions: bool,
    /// Refactorize exactly before these iteration counts instead of after `max_updates` updates (taken from a reference trace;
    /// the reference's own trigger depends on its LU fill, which is not reproduced).
    pub refactor_schedule: Option<Vec<usize>>,
    /// List the nonzeros of the entering column in the order the reference's sparse FTRAN produces them (needs `faithful_positions`).
    /// The Phase 1 ratio test depends on that order (`crate::ftran_order`).
    pub faithful_ftran: bool,
    /// The release profile's initial pricing (`SimplexPricing=-1` starts as `SimplexPricing=0`, `reverse-engineer/specs/primal-simplex.md`):
    /// no Devex weights and no maintained reduced costs; the duals are updated by one BTRAN per pivot and the entering variable is
    /// the best of a window of a persistent, bucket-sorted candidate list. Needs `faithful_ftran`. Phase 1 and the later switch
    /// to Devex pricing (`0x142500`) are not implemented. Off by default.
    pub pricing0: bool,
    /// Scale the model the way the reference does by default (`ScaleFlag`, power-of-two factors, `crate::scale`) and unscale the
    /// result. Off by default; on in `SolveOptions::release()`.
    pub scaling: bool,
    /// `PerturbValue` of the dual's cost perturbation (`crate::dual`, `reverse-engineer/specs/dual-simplex.md`): 0 turns it off (the controlled
    /// profile), 0.0002 is the library default.
    pub perturb: f64,
    /// The objective scale the model's costs carry (set by `crate::scale::solve`); the dual's objective value divides by it.
    pub obj_scale: f64,
    /// `BarConvTol` of `Method::Barrier`.
    pub barrier_tol: f64,
}

impl Default for SolveOptions {
    fn default() -> Self {
        SolveOptions {
            method: Method::Primal,
            pdhg: crate::pdhg::PdhgOptions::default(),
            feasibility_tol: 1e-6,
            optimality_tol: 1e-7,
            pivot_tol: 1e-9,
            max_iterations: 1_000_000,
            max_updates: 20,
            record_trace: true,
            faithful_weights: true,
            faithful_positions: true,
            refactor_schedule: None,
            faithful_ftran: true,
            pricing0: false,
            scaling: false,
            perturb: 0.0,
            obj_scale: 1.0,
            barrier_tol: 1e-8,
        }
    }
}

impl SolveOptions {
    /// The barrier at the reference's `Method=2`, `Crossover=0` defaults as far as recovered: the controlled profile (`BarCorrectors=1`, `ScaleFlag=0`,
    /// `BarOrder=0`), `BarConvTol=1e-8`, `BarIterLimit=1000`.
    pub fn barrier() -> Self {
        SolveOptions {
            method: Method::Barrier,
            max_iterations: 1000,
            ..SolveOptions::default()
        }
    }

    /// The barrier at the reference's release defaults (`ScaleFlag` default scaling, adaptive `BarCorrectors`), `specs/barrier-release.md`.
    pub fn barrier_release() -> Self {
        SolveOptions {
            scaling: true,
            ..SolveOptions::barrier()
        }
    }

    /// The PDHG at the reference's `Method=6` defaults (tolerances 1e-6, no iteration limit). The reference's controlled and release PDHG profiles are identical.
    pub fn pdhg() -> Self {
        SolveOptions {
            method: Method::Pdhg,
            ..SolveOptions::default()
        }
    }

    /// The solver without the reference-matching tie-breaking (own Devex floor, own basis position order, no FTRAN mirror).
    pub fn plain() -> Self {
        SolveOptions {
            faithful_weights: false,
            faithful_positions: false,
            faithful_ftran: false,
            ..SolveOptions::default()
        }
    }

    /// The reference's `dual-controlled` profile (`Method=1`, no scaling): the dual simplex of `crate::dual`.
    pub fn dual() -> Self {
        SolveOptions {
            method: Method::Dual,
            ..SolveOptions::default()
        }
    }

    /// The reference's `dual-release` profile (`Method=1`, library defaults) as far as it is recovered: the dual loop on the scaled model.
    pub fn dual_release() -> Self {
        SolveOptions {
            scaling: true,
            perturb: 0.0002,
            ..SolveOptions::dual()
        }
    }

    /// The reference's release profile as far as it is recovered: the faithful options with the initial `SimplexPricing=0` pricing.
    pub fn release() -> Self {
        SolveOptions {
            pricing0: true,
            scaling: true,
            ..SolveOptions::default()
        }
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Status {
    Optimal,
    Infeasible,
    Unbounded,
    IterationLimit,
    NumericalFailure,
}

/// One simplex iteration as seen by the diagnostic stream.
#[derive(Clone, Debug, PartialEq)]
pub struct TraceEvent {
    pub iteration: usize,
    /// Entering variable (internal index: columns first, then slacks).
    pub entering: usize,
    /// Leaving variable; `None` for a bound flip of the entering variable.
    pub leaving: Option<usize>,
    pub position: Option<usize>,
    pub step: f64,
    /// Variables flipped between their bounds in the same iteration (dual ratio test); empty for the primal.
    pub flips: Vec<usize>,
}

#[derive(Clone, Debug)]
pub struct SolveResult {
    pub status: Status,
    pub objective: f64,
    pub x: Vec<f64>,
    /// Row duals in the sign convention of the original objective sense.
    pub dual: Vec<f64>,
    pub reduced_cost: Vec<f64>,
    pub basis: Basis,
    pub iterations: usize,
    pub refactorizations: usize,
    pub trace: Vec<TraceEvent>,
    /// Row multipliers `y` proving infeasibility (see `verify::check_farkas`) when the status is `Infeasible`.
    pub farkas: Option<Vec<f64>>,
    /// Recession direction over the original columns proving unboundedness when the status is `Unbounded`.
    pub ray: Option<Vec<f64>>,
}

#[derive(Debug)]
pub enum SolveError {
    InvalidBasis(String),
    SingularBasis,
    /// A model feature or option outside the recovered case of the method (named); never answered approximately.
    Unsupported(String),
}

impl std::fmt::Display for SolveError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            SolveError::InvalidBasis(s) => write!(f, "invalid basis: {s}"),
            SolveError::SingularBasis => write!(f, "starting basis is singular"),
            SolveError::Unsupported(s) => write!(f, "unsupported: {s}"),
        }
    }
}

impl std::error::Error for SolveError {}

const NONE: usize = usize::MAX;

pub(crate) struct Internal {
    pub(crate) n: usize,
    pub(crate) m: usize,
    pub(crate) cols: Vec<Column>,
    /// The same columns in the model's entry order (the reference's LU structure and FTRAN depend on it).
    pub(crate) fcols: Vec<Column>,
    pub(crate) lo: Vec<f64>,
    pub(crate) up: Vec<f64>,
    pub(crate) rhs: Vec<f64>,
    pub(crate) cost: Vec<f64>,
    /// Row-wise copy of `cols` (all columns, slacks included): row `i` holds `(j, a)` in ascending `j`, a column's entries in its order.
    rwise: Vec<Vec<(usize, f64)>>,
    /// Offset of each column in the concatenated value array (prefix sums of the column lengths).
    offs: Vec<usize>,
}

impl Internal {
    /// `rank_dot(&cols[j], v, rank)` for every column at once, iterating only the rows with a rank (in rank order), so each
    /// column's sum is accumulated in the same order from 0.0 as `rank_dot` does.
    fn rank_products(&self, v: &[f64], rank: &[usize]) -> Vec<f64> {
        let mut order: Vec<(usize, usize)> = (0..self.m)
            .filter(|&i| rank[i] != usize::MAX)
            .map(|i| (rank[i], i))
            .collect();
        order.sort_unstable();
        let mut acc = vec![0.0f64; self.n + self.m];
        for &(_, i) in &order {
            let vi = v[i];
            for &(j, a) in &self.rwise[i] {
                acc[j] = acc[j] + a * vi;
            }
        }
        acc
    }
}

fn is_inf(v: f64) -> bool {
    v.abs() >= INF
}

/// Internal layout: structurals, then one range column per ranged inequality row (as the reference does), then
/// one slack per row. A ranged row `rl <= a x <= ru` becomes `a x + s + g = ru` with `s` fixed at 0 and
/// `g in [0, ru - rl]`; a ranged `>=` row uses `-1` for both unit columns.
pub(crate) fn build(model: &LpModel) -> Internal {
    let (nm, m) = (model.num_cols(), model.num_rows());
    let mut cols: Vec<Column> = (0..nm).map(|j| model.col(j).collect()).collect();
    // The model keeps the file's entry order (the crash tie-break needs it); the solver works on columns sorted by row.
    let mut fcols = cols.clone();
    for c in &mut cols {
        c.sort_by_key(|e| e.0);
    }
    let mut lo = model.col_lo.clone();
    let mut up = model.col_up.clone();
    let mut rhs = vec![0.0; m];
    let mut slack: Vec<(f64, f64, f64)> = Vec::with_capacity(m);
    for i in 0..m {
        let (rl, ru) = (model.row_lo[i], model.row_up[i]);
        let (sign, r, slo, sup) = match (is_inf(rl), is_inf(ru)) {
            (false, false) if rl == ru => (1.0, ru, 0.0, 0.0),
            (false, false) => {
                cols.push(vec![(i, 1.0)]);
                lo.push(0.0);
                up.push(model.range_width(i));
                (1.0, ru, 0.0, 0.0)
            }
            (true, false) => (1.0, ru, 0.0, INF),
            (false, true) => (-1.0, rl, 0.0, INF),
            (true, true) => (1.0, 0.0, -INF, INF),
        };
        rhs[i] = r;
        slack.push((sign, slo, sup));
    }
    let n = cols.len();
    for (i, &(sign, slo, sup)) in slack.iter().enumerate() {
        cols.push(vec![(i, sign)]);
        lo.push(slo);
        up.push(sup);
    }
    let sgn = if model.sense == Sense::Maximize {
        -1.0
    } else {
        1.0
    };
    let mut cost: Vec<f64> = model.obj.iter().map(|c| sgn * c).collect();
    cost.resize(n + m, 0.0);
    fcols.extend(cols[fcols.len()..].iter().cloned());
    let mut rwise: Vec<Vec<(usize, f64)>> = vec![Vec::new(); m];
    for (j, c) in cols.iter().enumerate() {
        for &(i, a) in c {
            rwise[i].push((j, a));
        }
    }
    let mut offs = vec![0usize; n + m];
    for j in 1..n + m {
        offs[j] = offs[j - 1] + cols[j - 1].len();
    }
    Internal {
        n,
        m,
        cols,
        fcols,
        lo,
        up,
        rhs,
        cost,
        rwise,
        offs,
    }
}

/// Crash basis of the reference (primal simplex setup, `0x1387d0`): each equality row (fixed slack) may receive one structural
/// column in place of its slack. Singleton columns with an infinite upper bound go first; then columns with fewer than 10
/// entries are taken bucket by bucket (by static entry count, buckets scanned from the highest index first) while their dynamic
/// count of still-unassigned equality rows is within the pass limit; a column takes the unassigned equality row with the
/// largest `|a|` (first on ties) and marks all of its equality rows used. Recovered from the assignments of the reference on
/// ADLITTLE, BLEND, SC50A and BORE3D (the set of crash columns matches on all four).
pub fn crash_basis(model: &LpModel) -> Basis {
    let p = build(model);
    let (n, m) = (p.n, p.m);
    let file_cols: Vec<Column> = (0..n)
        .map(|j| {
            if j < model.num_cols() {
                model.col(j).collect()
            } else {
                p.cols[j].clone()
            }
        })
        .collect();
    let mut b = Basis::slack(n, m);
    for j in 0..n {
        b.at_upper[j] = is_inf(p.lo[j]) && !is_inf(p.up[j]);
    }
    // used[r]: row r is not an equality row, or already has a column.
    let mut used: Vec<bool> = (0..m).map(|i| p.up[n + i] - p.lo[n + i] > 0.0).collect();
    let mut cols_of_row: Vec<Vec<usize>> = vec![Vec::new(); m];
    for j in 0..n {
        for &(r, _) in &p.cols[j] {
            cols_of_row[r].push(j);
        }
    }
    // Fixed columns never enter the crash.
    let mut eligible: Vec<bool> = (0..n).map(|j| p.up[j] - p.lo[j] > 0.0).collect();
    let assign =
        |b: &mut Basis, used: &mut Vec<bool>, eligible: &mut Vec<bool>, j: usize, r: usize| {
            b.basic[r] = j;
            used[r] = true;
            eligible[j] = false;
        };
    // Singleton columns take their row in two sweeps over the columns: first those that start at a value inside their bounds
    // (`rhs / a >= lo`; D2Q06C: of `-x1 = 6` and `x2 = 6` the second takes the row), then the rest (PILOT4 takes `x = -2.171`).
    for sweep in 0..2 {
        for j in 0..n {
            if eligible[j] && p.cols[j].len() == 1 && is_inf(p.up[j]) {
                let (r, a) = p.cols[j][0];
                if !used[r] && (sweep == 1 || p.rhs[r] / a >= p.lo[j]) {
                    assign(&mut b, &mut used, &mut eligible, j, r);
                }
            }
        }
    }
    // The reference's per-column counter starts at -1, so it holds the number of unassigned equality rows minus one.
    let mut count = vec![-1i64; n];
    for r in 0..m {
        if !used[r] {
            for &j in &cols_of_row[r] {
                count[j] += 1;
            }
        }
    }
    let mut buckets: Vec<Vec<usize>> = vec![Vec::new(); 11];
    for j in (0..n).rev() {
        if eligible[j] && p.cols[j].len() < 10 {
            buckets[p.cols[j].len()].push(j);
        }
    }
    let mut marked = used.clone();
    for bk in 1..10usize {
        for t in 0..bk as i64 {
            for &j in &buckets[bk] {
                if !(0..=t).contains(&count[j]) {
                    continue;
                }
                let mut best: Option<(usize, f64)> = None;
                for &(r, a) in &file_cols[j] {
                    if !marked[r] && best.map_or(a.abs() > 0.0, |(_, v)| a.abs() > v) {
                        best = Some((r, a.abs()));
                    }
                }
                if let Some((r, _)) = best {
                    assign(&mut b, &mut used, &mut eligible, j, r);
                }
                for &(r, _) in &p.cols[j] {
                    if !marked[r] {
                        marked[r] = true;
                        for &k in &cols_of_row[r] {
                            if (0..10).contains(&count[k]) {
                                count[k] -= 1;
                            }
                        }
                    }
                }
            }
        }
    }
    flip_to_upper(&p, &mut b);
    b
}

/// The reference's sort of a key array with a companion index array (`0x8bee60`, `0x8bf1f0`): ascending, not stable. Ranges of 100
/// or more are partitioned around the median of three (the smaller side is sorted by a recursive call, the larger one goes on in
/// the loop, at most 50 rounds); what is left, and every smaller range, is shell-sorted with gaps `n/2`, then `gap/2.2`, ending
/// with 1 (a gap of 2 goes to 1). Ties therefore keep the input order only in the insertion-sort stage.
pub(crate) fn ref_sort<T: Copy>(a: &mut [f64], idx: &mut [T]) {
    fn quick<T: Copy>(a: &mut [f64], idx: &mut [T], mut lo: usize, mut hi: usize, mut depth: u32) {
        while hi - lo >= 100 && depth < 50 {
            let n = hi - lo;
            let (x, y, z) = (a[lo], a[lo + n / 2], a[hi - 1]);
            let p = if (y <= x && x <= z) || (z <= x && x <= y) {
                x
            } else if (x <= y && y <= z) || (z <= y && y <= x) {
                y
            } else {
                z
            };
            let (mut i, mut r, mut lp, mut rp) = (lo, hi, lo, 0usize);
            if lo < hi && p > a[lo] {
                while i < hi && p > a[i] {
                    i += 1;
                }
                lp = i;
            }
            while i < r {
                let v = a[i];
                if p > v {
                    a.swap(i, lp);
                    idx.swap(i, lp);
                    i += 1;
                    lp += 1;
                } else if v == p {
                    i += 1;
                } else {
                    loop {
                        r -= 1;
                        rp += 1;
                        let w = a[r];
                        if p > w {
                            if i < r {
                                // a[i] <- a[lp], a[lp] <- a[r], a[r] <- v
                                let (il, ir) = (idx[lp], idx[r]);
                                a[i] = a[lp];
                                a[lp] = w;
                                a[r] = v;
                                idx[r] = idx[i];
                                idx[i] = il;
                                idx[lp] = ir;
                                i += 1;
                                lp += 1;
                            }
                            break;
                        } else if w == p {
                            if i < r {
                                a[i] = w;
                                a[r] = v;
                                idx.swap(i, r);
                                i += 1;
                            }
                            break;
                        } else if i >= r {
                            break;
                        }
                    }
                }
            }
            depth += 1;
            let (left, right) = (lp - lo, rp);
            if left > right {
                quick(a, idx, hi - rp, hi, depth);
                hi = lp;
            } else {
                quick(a, idx, lo, lp, depth);
                lo = hi - rp;
            }
        }
        shell(a, idx, lo, hi);
    }

    fn shell<T: Copy>(a: &mut [f64], idx: &mut [T], lo: usize, hi: usize) {
        let mut gap = (hi - lo) / 2;
        while gap > 0 {
            for s in lo + gap..hi {
                let (key, ki) = (a[s], idx[s]);
                let mut k = s;
                while k >= lo + gap && a[k - gap] > key {
                    a[k] = a[k - gap];
                    idx[k] = idx[k - gap];
                    k -= gap;
                }
                a[k] = key;
                idx[k] = ki;
            }
            gap = if gap == 2 {
                1
            } else {
                (gap as f64 / 2.2) as usize
            };
        }
    }
    quick(a, idx, 0, a.len(), 0);
}

/// The reference's bound placement after the crash (`0x1387d0`, decompiled lines 405-560): each nonbasic column with a finite
/// range, taken in ascending order of `(up - lo) * cost` (ties by descending index, which reproduces BOEING2), moves to its upper
/// bound when that lowers the total infeasibility of the row slacks (a slack is feasible in `[0, up]`), or leaves it unchanged
/// while the cost is negative. The row activities follow the placements already made.
fn flip_to_upper(p: &Internal, b: &mut Basis) {
    let (n, m) = (p.n, p.m);
    let mut is_basic = vec![false; n + m];
    for &j in &b.basic {
        is_basic[j] = true;
    }
    // Slack-space activity (every structural at its bound, the crash columns too): `(rhs - A x_N) / sign` for every row.
    let mut act = p.rhs.clone();
    for j in 0..n {
        let v = if b.at_upper[j] {
            p.up[j]
        } else if is_inf(p.lo[j]) {
            0.0
        } else {
            p.lo[j]
        };
        if v != 0.0 {
            for &(i, a) in &p.cols[j] {
                act[i] -= a * v;
            }
        }
    }
    let sign: Vec<f64> = (0..m).map(|i| p.cols[n + i][0].1).collect();
    let mut sact: Vec<f64> = (0..m).map(|i| act[i] / sign[i]).collect();
    let mut cand: Vec<usize> = (0..n)
        .filter(|&j| !b.at_upper[j] && !is_inf(p.lo[j]) && !is_inf(p.up[j]) && p.up[j] > p.lo[j])
        .collect();
    let key = |j: usize| (p.up[j] - p.lo[j]) * p.cost[j];
    let mut keys: Vec<f64> = cand.iter().map(|&j| key(j)).collect();
    ref_sort(&mut keys, &mut cand);
    for j in cand {
        let range = p.up[j] - p.lo[j];
        let mut gain = 0.0;
        for &(i, a) in &p.cols[j] {
            let hi = p.up[n + i];
            let infeas = |x: f64| {
                if x >= 0.0 {
                    if hi < x {
                        x - hi
                    } else {
                        0.0
                    }
                } else {
                    -x
                }
            };
            let x = sact[i];
            gain += infeas(x) - infeas(x - a / sign[i] * range);
        }
        if gain > 0.0 || (gain == 0.0 && p.cost[j] < 0.0) {
            // A crash column keeps its basic status (the reference writes the positions afterwards) but its placement moves the row activities.
            if !is_basic[j] {
                b.at_upper[j] = true;
            }
            for &(i, a) in &p.cols[j] {
                sact[i] -= a / sign[i] * range;
            }
        }
    }
}

struct Solver<'a> {
    p: &'a Internal,
    opts: &'a SolveOptions,
    /// Variable dropped last from a singular supplied basis: not priced until the basis is primal feasible.
    excluded: Option<usize>,
    restart: bool,
    pre_restart: bool,
    /// The duals were recomputed at the start of this iteration (the end of phase 1), so the scan that finds nothing does not recompute them again.
    fresh_duals: bool,
    /// The driver restarted inside the last pricing call: the refactorization check read before it is stale.
    restarted: bool,
    basic: Vec<usize>,
    pos_of: Vec<usize>,
    at_upper: Vec<bool>,
    x: Vec<f64>,
    d: Vec<f64>,
    w: Vec<f64>,
    in_r0: Vec<bool>,
    factor: Factor,
    refactorizations: usize,
    /// Persistent pricing candidate list and its bucket level (see `bucket`).
    level: usize,
    list: Vec<usize>,
    /// Feasibility class (-1 below, 1 above, 0 inside) of each basic variable at the last phase 1 cost evaluation.
    p1_prev: Vec<(usize, i8)>,
    /// Variable that left the basis in the last pivot.
    leave_prev: Option<usize>,
    /// Structural columns of each row, ascending.
    rows: Vec<Vec<usize>>,
    /// The reference's row-wise copy of the nonbasic structurals, in its storage order: a column that enters the basis is
    /// removed by moving the row's last entry into its place, a column that leaves is appended.
    rrows: Vec<Vec<usize>>,
    /// The last `price` call scanned a non-empty list and found nothing.
    list_failed: bool,
    /// The mirror's penalised dual multipliers `y` by row, from the last `dual_btran_work` (they carry the reference's rounding).
    y_ref: Option<Vec<f64>>,
    /// Position of each row in the output list of the last dual BTRAN (`usize::MAX` for rows without an entry).
    y_rank: Vec<usize>,
    /// Nonbasic variables that left the basis in phase 2: the reference keeps their computed value (not the bound) for the next
    /// recomputation of `x_B`.
    unsnapped: std::collections::HashSet<usize>,
    /// Distance of each basic variable in the frame of the ratio test (`x - lo`, or `up - x` without a lower bound).
    dst: Vec<f64>,
    /// The reference's persistent rhs shift by row: each time a variable leaves outside its bound, `a_ij * (bound - x)` is added
    /// for its column; the recomputation of `x_B` starts from `b + shift` and the shift is never reset (PEROLD: it stays after the
    /// variable re-enters).
    rhs_shift: Vec<f64>,
    /// Mirror of the reference's FTRAN structure (`faithful_ftran`).
    mirror: Option<crate::ftran_order::RefStruct>,
    /// The reference's work counter (`crate::work`).
    work: crate::work::Work,
    /// The last pivot, whose primal update the reference runs at the start of the next iteration (`0x105f20`).
    last_pivot: Option<LastPivot>,
    /// Whether basic variables are infeasible (`ctx+0x2d0 != 0`): the list update and the pricing scan then take the dual
    /// tolerance `P[0x3c0] * P[0]` = `1e-7 * 1e6` instead of the plain one.
    p1_tol_on: bool,
    /// Whether the last pricing call recomputed the duals (the list scan failed).
    recomputed: bool,
    /// The scan that follows a recomputation of the duals: the reference's recomputation already covered it, so a scan that finds nothing is not booked.
    rescan: bool,
    /// The final factorization of a solve that ended with right-hand-side shifts: the reference removes the shifts in a dual clean-up pass
    /// that factorizes the basis in ascending index order (`0xf6a00`), instead of refactorizing in pivot order.
    cleanup: bool,
    /// The pivot row and `extra` of the last list update from the pivot step (`0x1442ab`), which `0x107e14` repeats.
    last_upd: Option<(Vec<usize>, usize)>,
    /// The pivot row's values (`arj`) of that update, and the row being built for the next one.
    last_row: Vec<(usize, f64)>,
    pending_row: Option<Vec<(usize, f64)>>,
    /// Work of the last column replacement, which the reference does after the next refactorization check (and not at all
    /// when it refactors instead).
    upd_work: Vec<(u32, f64, f64)>,
    /// State of the `pricing0` option (`None` when it is off or after the switch to Devex pricing).
    p0: Option<P0>,
    /// Set by `price0` when the scan found nothing after many updates: refactorize and price again.
    want_refactor: bool,
    /// Counters of the method decision (`crate::method`).
    sw: crate::method::Method,
}

/// State of the release profile's initial pricing (`0xfbe60`, list build `0xfecb0`, dual update `0xfb790`).
struct P0 {
    /// All variables, sorted at each rebuild by bucket key (descending), then by index (descending); the scan permutes it in place.
    list: Vec<usize>,
    /// Scan position (`P+0x54`).
    pos: usize,
    /// The reference's stored `d` (sign-adjusted so that attractive is negative), `ctx+0x1e0`.
    dstore: Vec<f64>,
    /// Dual multipliers, kept by one BTRAN per pivot.
    y: Vec<f64>,
    /// The BTRAN input of the next dual update: the leaving position with the entering `d`, then (phase 1) the cost shifts of the basic
    /// variables that changed feasibility class, in the reference's order.
    pending: Option<Vec<(usize, f64)>>,
    /// Reduced-cost tolerance of the scan and the rebuild: `ctx+0x388`, or in phase 1 `P[0x3c0] * P[0]` (`P0_DTOL_P1`).
    tol: f64,
    /// The duals are recomputed and the list rebuilt at the next pricing call.
    rebuild: bool,
    /// Pricing calls so far (the dump's call index).
    calls: usize,
    /// Scan window (`P+0x50`): set at each refactorization to `max(10, trunc(avg / 2))`, `avg = trunc(S / U)` the mean FTRAN result count of the `U` updates since the last one.
    window: usize,
    /// The bucket-chain variant of the candidate list (`P+0x80 == 2`, see `Buckets`), for models with at least twice as many columns as
    /// rows per column entry (`0x137fc9`); `list` is unused then.
    bk: Option<Buckets>,
}

/// The reference's bucket candidate list when `P+0x80 == 2` (`0xfecb0` second loop, scan `0xfbe60` from `0xfbfcc`): `next` is the array
/// `P+0x40` (the successor in the chain of the same bucket), `head[k]` / `tail[k]` the chains `P+0x88 + 4k` / `P+0x21c + 4k` for the keys
/// `1..=100` (`-1` empty), `visited` the entries scanned since the rebuild (`P+0x3c8`).
struct Buckets {
    head: [i64; 101],
    tail: [i64; 101],
    next: Vec<i64>,
    visited: usize,
    /// The chains were built once; later refactorizations only recompute the duals (`0xe4092`: `0x145520` is skipped while `P+0x78 > 0`).
    built: bool,
    /// `P+0x78 = r8d` of the gate (`p0_bucket_mode`); the sign flip at each factorization (`0x1428ce`..`0x1428d3`) makes `0xe4092` rebuild when the
    /// factorization is a clean-up (`ctx+0x294` in 3..8) or its number (`ctx+0x2c0`, the first is 1) is 1 modulo `r8d`.
    r8: usize,
    /// The refactorization just made calls the rebuild (the first one always does).
    full: bool,
}

impl Buckets {
    /// Appends `j` to the chain of key `k`.
    fn append(&mut self, k: usize, j: usize) {
        if self.tail[k] != -1 {
            self.next[self.tail[k] as usize] = j as i64;
        } else {
            self.head[k] = j as i64;
        }
        self.tail[k] = j as i64;
        self.next[j] = -1;
    }
}

struct LastPivot {
    /// Rows of the FTRAN result in its output order, without the pivot row.
    rows: Vec<usize>,
    /// Whether the reference held the FTRAN result as a dense vector.
    dense: bool,
    q: usize,
    leave: usize,
    alpha_rq: f64,
}

fn dot(col: &Column, v: &[f64]) -> f64 {
    col.iter().map(|&(r, a)| a * v[r]).sum()
}

/// Column dot product as the reference's vectorised loop evaluates it (`0x102d20`): below `peel + 8` terms it is sequential,
/// otherwise the peeled head goes to lane 0, eight lane-wise partial sums follow, the lanes reduce pairwise and the tail is
/// sequential. `peel` is the parity of the column's offset in the value array.
/// Dot product summing the terms in the rank order of the vector's support (the order of the reference's BTRAN output list).
fn rank_dot(col: &Column, v: &[f64], rank: &[usize]) -> f64 {
    let mut terms: Vec<(usize, f64)> = col
        .iter()
        .filter(|&&(i, _)| rank[i] != usize::MAX)
        .map(|&(i, a)| (rank[i], a * v[i]))
        .collect();
    terms.sort_by_key(|t| t.0);
    terms.iter().fold(0.0, |acc, t| acc + t.1)
}

/// `init - sum a_i y_i` accumulated as the reference's read-out kernel does: sequentially below `peel + 8` entries, else the peeled head
/// and the eight-lane body, the lanes folded pairwise, then the tail (as `dot_lanes`, with the start value in lane 0).
fn dot_lanes_from(col: &Column, y: &[f64], peel: usize, init: f64) -> f64 {
    let n = col.len();
    let t = |i: usize| -(col[i].1 * y[col[i].0]);
    if n < peel + 8 {
        let mut acc = init;
        for i in 0..n {
            acc += t(i);
        }
        return acc;
    }
    let n8 = n - ((n - peel) & 7);
    let mut lanes = [0.0f64; 8];
    lanes[0] = init;
    for i in 0..peel {
        lanes[0] += t(i);
    }
    for i in peel..n8 {
        lanes[(i - peel) & 7] += t(i);
    }
    let mut tot = ((lanes[0] + lanes[2]) + (lanes[4] + lanes[6]))
        + ((lanes[1] + lanes[3]) + (lanes[5] + lanes[7]));
    for i in n8..n {
        tot += t(i);
    }
    tot
}

fn dot_lanes(col: &Column, v: &[f64], peel: usize) -> f64 {
    let n = col.len();
    let t = |i: usize| col[i].1 * v[col[i].0];
    if n < peel + 8 {
        let mut acc = 0.0;
        for i in 0..n {
            acc += t(i);
        }
        return acc;
    }
    let n8 = n - ((n - peel) & 7);
    let mut lanes = [0.0f64; 8];
    for i in 0..peel {
        lanes[0] += t(i);
    }
    for i in peel..n8 {
        lanes[(i - peel) & 7] += t(i);
    }
    let mut tot = ((lanes[0] + lanes[2]) + (lanes[4] + lanes[6]))
        + ((lanes[1] + lanes[3]) + (lanes[5] + lanes[7]));
    for i in n8..n {
        tot += t(i);
    }
    tot
}

impl<'a> Solver<'a> {
    fn nonbasic_value(&self, j: usize) -> f64 {
        let (l, u) = (self.p.lo[j], self.p.up[j]);
        if self.at_upper[j] && !is_inf(u) {
            u
        } else if !is_inf(l) {
            l
        } else if !is_inf(u) {
            u
        } else {
            0.0
        }
    }

    /// Basis positions after a refactorization: the order of the reference LU's pivot steps (`crate::lu_order`, recovered in
    /// `reverse-engineer/specs/lu-factor.md`), applied to the basis in its previous position order.
    fn reference_order(
        &self,
        fac_no: usize,
    ) -> (
        Vec<usize>,
        crate::lu_order::Structure,
        Vec<(bool, Vec<usize>, Vec<f64>)>,
    ) {
        let n = self.p.n;
        let cols: Vec<(bool, Vec<usize>, Vec<f64>)> = self
            .basic
            .iter()
            .map(|&v| {
                (
                    v >= n,
                    self.p.fcols[v].iter().map(|&(r, _)| r).collect(),
                    self.p.fcols[v].iter().map(|&(_, a)| a).collect(),
                )
            })
            .collect();
        crate::lu_order::set_factorization_number(fac_no);
        let st = crate::lu_order::factor_structure(self.p.m, &cols);
        let order = st.pivots.iter().map(|&(_, pos)| self.basic[pos]).collect();
        (order, st, cols)
    }

    fn refactor(&mut self) -> Result<(), SolveError> {
        let mut sw_avg = 0usize;
        if self.opts.faithful_positions {
            if self.cleanup || self.restart {
                self.basic.sort_unstable();
            }
            // The driver restart after a repaired basis (`0xcda5b`) factors like the starting factorization.
            let (order, st, rows) = self.reference_order(if self.restart {
                1
            } else {
                self.refactorizations + 1
            });
            if crate::dump::on("luin") {
                crate::dump::emit(
                    "luin",
                    serde_json::json!({"n": self.refactorizations, "cols": self.basic, "out": order, "prow": st.pivots.iter().map(|p| p.0).collect::<Vec<_>>(), "phase": st.phase}),
                );
            }
            debug_assert_eq!(order.len(), self.basic.len());
            for (pos, &j) in order.iter().enumerate() {
                self.pos_of[j] = pos;
            }
            self.basic = order;
            if self.opts.faithful_ftran {
                if let Some(old) = self.mirror.as_ref() {
                    let (u, sum) = old.window_stats();
                    sw_avg = if u > 0 { (sum / u as f64) as usize } else { 0 };
                }
                if let (Some(old), Some(p0)) = (self.mirror.as_ref(), self.p0.as_mut()) {
                    // `0xf15c0` (`0xf4562`): the scan window follows the mean FTRAN result count of the updates made before the last one
                    // (the refactorization is decided before that update is counted).
                    let (u, sum) = old.window_stats();
                    let avg = if u > 0 { (sum / u as f64) as usize } else { 0 };
                    p0.window = (avg / 2).max(P0_WINDOW);
                }
                if let Some(old) = self.mirror.as_ref() {
                    let est = if self.restart {
                        0.0
                    } else if self.cleanup || self.pre_restart {
                        old.discard_work_unused()
                    } else {
                        old.discard_work()
                    };
                    if est > 0.0 {
                        self.work.add(0xf461f, est, crate::work::CTX);
                    }
                }
                let mi = crate::ftran_order::RefStruct::fresh(self.p.m, self.p.n, &st, &rows);
                if crate::dump::on("luout") {
                    let (diag, l) = mi.numerics();
                    crate::dump::emit(
                        "luout",
                        serde_json::json!({"n": self.refactorizations, "diag": diag, "l": l}),
                    );
                }
                let lu = st.work(&rows, mi.has_levels());
                if self.cleanup {
                    // The clean-up pass's setup (`0xe84b9`): two units per nonzero of the constraint matrix; the dual simplex install of the
                    // factor then books `3 m` (`0xf48ae`), `nnz(U)` (`0x1327dc`) and `3 m` (`0x1329c1`) at the factor object's scale.
                    let nnz_a: usize = (0..self.p.n).map(|j| self.p.cols[j].len()).sum();
                    self.work.add(0xe84b9, 2.0 * nnz_a as f64, crate::work::CTX);
                }
                self.work.factor_driver(&lu, crate::work::CTX);
                for (site, est) in lu.into_iter().chain(mi.fresh_work()) {
                    self.work.add(site, est, crate::work::CTX);
                }
                self.work.factor_done();
                if self.cleanup {
                    let m = self.p.m as f64;
                    self.work.add(0xf48ae, 3.0 * m, crate::work::FAC);
                    self.work.add(0x1327dc, mi.u_nnz() as f64, crate::work::FAC);
                    self.work.add(0x1329c1, 3.0 * m, crate::work::FAC);
                }
                self.mirror = Some(mi);
            }
        }
        let cols: Vec<&Column> = self.basic.iter().map(|&j| &self.p.cols[j]).collect();
        self.factor = Factor::new(self.p.m, &cols).map_err(|_| SolveError::SingularBasis)?;
        self.refactorizations += 1;
        self.list.clear();
        if self.cleanup {
            // The clean-up pass removes the shifts before it solves for the basic values.
            self.rhs_shift.clear();
        }
        // The bucket-chain variant recomputes x_B and the distances only at its full refactorizations; at the others the
        // reference keeps its updated values (and runs no dense FTRAN).
        let keep_x = !self.cleanup
            && self.p0.as_ref().is_some_and(|p0| {
                p0.bk
                    .as_ref()
                    .is_some_and(|bk| self.refactorizations % bk.r8 != 1)
            });
        if !keep_x {
            self.recompute_primal();
        }
        if self.cleanup {
            self.cleanup_work();
        } else {
            self.primal_full_work();
            self.dual_btran_work();
        }
        self.recompute_dual();
        if let Some(p0) = self.p0.as_mut() {
            p0.rebuild = true;
            p0.pending = None;
            if let Some(bk) = p0.bk.as_mut() {
                bk.full = self.cleanup || self.refactorizations % bk.r8 == 1;
            }
        }
        if !self.cleanup
            && self
                .p0
                .as_ref()
                .is_some_and(|p0| p0.bk.as_ref().is_none_or(|bk| bk.full))
        {
            // `0xe4092` runs `0x145520` (duals, residual check, rebuild) after the factorization.
            self.residual_check_work();
        }
        if self.sw.active && !self.cleanup {
            self.method_decision(sw_avg);
            if self.sw.perturb_due {
                self.perturb();
            }
        }
        Ok(())
    }

    /// The bound perturbation (`0x141a60`, called from the pricing driver after a factorization): every basic variable at a bound (type
    /// `L`, or `U`) moves off it by `(r % m + m) * size / m`, `r` the next output of a 64-bit LCG seeded with the iteration count, in
    /// position order. The reference changes its distance array and the persistent right-hand side shift (the basic values follow at the
    /// next recomputation), not `x_B`.
    fn perturb(&mut self) {
        let (p, m) = (self.p, self.p.m);
        let tol = self.opts.feasibility_tol;
        let size = self.sw.pert;
        let mut state = self.sw.iters as u64;
        for pos in 0..m {
            let j = self.basic[pos];
            let Some((dist, range, typ, mirrored)) =
                classify(self.x[j], p.lo[j], p.up[j], tol, self.dst[j])
            else {
                continue;
            };
            if typ != b'L' && typ != b'U' {
                continue;
            }
            state = state
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            let mut t = (state >> 32) as u32;
            t ^= t >> 11;
            t ^= (t << 7) & 0x9d2c5680;
            t ^= (t << 15) & 0xefc60000;
            t ^= t >> 18;
            let r = (t & 0x7fff_ffff) as i64;
            let k = ((r % m as i64) + m as i64) as f64 * size / m as f64;
            let amt = if typ == b'L' { k } else { -k };
            let nd = dist + amt;
            if (typ == b'L' && nd > range) || (typ == b'U' && 0.0 > nd) {
                continue;
            }
            self.dst[j] = nd;
            if self.rhs_shift.is_empty() {
                self.rhs_shift = vec![0.0; m];
            }
            let a = if mirrored { -amt } else { amt };
            for &(i, v) in &p.cols[j] {
                self.rhs_shift[i] += v * a;
            }
        }
        self.sw.perturbed();
    }

    /// The method decision after a factorization (`0x142500`): switches from `pricing0` to Devex pricing, restarting the weights and the list.
    fn method_decision(&mut self, avg: usize) {
        let mi = self.mirror.as_ref().unwrap();
        let inp = crate::method::Inputs {
            rows: self.p.m,
            cols: self.p.n,
            avg,
            l_total: mi.l_nnz() + self.p.m,
            u_total: mi.u_nnz() + self.p.m,
            infeasible: self.has_infeasible_basic(),
        };
        let before = self.sw.clone();
        let switch = self.sw.decide(&inp);
        if crate::dump::on("switch") {
            crate::dump::emit(
                "switch",
                serde_json::json!({"call": self.refactorizations - 1, "iters": before.iters, "degenerate": before.degenerate,
                "snap_degenerate": before.snap_degenerate, "snap_iters": before.snap_iters, "avg": avg, "l_total": inp.l_total, "u_total": inp.u_total, "switch": switch, "mode": self.sw.mode}),
            );
        }
        if switch {
            // The switch books the Devex initialization (`0x142c4c`: `3 (n + m)`; `0x132d88`: `2 nnz(A)`), then the reference runs the Devex
            // factorization handling too: the primal update's `x_B` pass (`0x107b57`) and the reduced-cost pass (`0x1021b0`), on top of the
            // `pricing0` handling already booked.
            let nnz_a: usize = (0..self.p.n).map(|j| self.p.cols[j].len()).sum();
            self.work.add(
                0x142c4c,
                3.0 * (self.p.n + self.p.m) as f64,
                crate::work::CTX,
            );
            if before.mode == 0 {
                self.work
                    .add(0x132d88, 2.0 * nnz_a as f64, crate::work::CTX);
            }
            self.p0 = None;
            // The Devex initialization builds the row-wise copy afresh (`0x132d88`, `2 nnz(A)`): nonbasic, non-fixed structurals in index order.
            if before.mode == 0 {
                self.rrows = vec![Vec::new(); self.p.m];
                for j in 0..self.p.n {
                    if self.pos_of[j] == NONE && !self.is_fixed(j) {
                        for &(i, _) in &self.p.cols[j] {
                            self.rrows[i].push(j);
                        }
                    }
                }
            }
            if before.mode == 0 {
                self.primal_full_work();
                self.reduced_cost_work();
                // The reduced-cost pass (`0x1021b0`) recomputes the nonbasic d sequentially in stored column order, not with the
                // `pricing0` eight-lane kernel.
                self.recompute_dual();
            }
            self.reset_weights();
            self.sw.ratio_sum = 0.0;
            self.sw.ratio_cnt = 1;
            self.list.clear();
            self.level = 0;
        }
    }

    /// The dual clean-up pass that follows a primal optimum with shifts (`0xf6a00`), after its factorization and the recomputation of the duals:
    /// the pricing setup (`0xe7000`, `3 m + 3.5 n`) and the final solution handling (`0xc362f`, `0xc1fe3`, `0xbf981`, `0xcd191`).
    fn cleanup_work(&mut self) {
        let (m, n) = (self.p.m as f64, self.p.n as f64);
        self.dual_btran_work();
        self.work.add(0xe7000, 3.0 * m + 3.5 * n, crate::work::CTX);
        // `0xc142b`: the nonbasic structurals at a nonzero value, whose columns the solution read-out sums (`0xc0ef0`).
        let r: usize = (0..self.p.n)
            .filter(|&j| self.pos_of[j] == NONE && self.nonbasic_value(j) != 0.0)
            .map(|j| self.p.cols[j].len())
            .sum();
        self.work
            .add(0xc142b, 4.0 * m + 2.0 * r as f64, crate::work::CTX);
        self.work.add(0xc362f, m, crate::work::CTX);
        self.work.add(0xc1fe3, 5.0 * m, crate::work::CTX);
        self.work.add(0xbf981, 2.0 * (n + m), crate::work::CTX);
        self.work.add(0xcd191, 2.0 * m, crate::work::CTX);
    }

    /// The reference's residual check (`0x145a36` in `0x145520`): after the duals are recomputed, when no basic variable is
    /// infeasible and the factorization count is 4 modulo 5, it books `(rows + 2 * sum(len + 1 over basic structurals)) * c` for the pass that forms their reduced costs.
    fn residual_check_work(&mut self) {
        if self.mirror.is_none() || self.refactorizations % 5 != 4 || self.has_infeasible_basic() {
            return;
        }
        let nnz: usize = self
            .basic
            .iter()
            .filter(|&&j| j < self.p.n)
            .map(|&j| self.p.cols[j].len() + 1)
            .sum();
        self.work
            .add(0x145a36, (self.p.m + 2 * nnz) as f64, crate::work::CTX);
    }

    /// Work of the reference's primal update (`0x105f20`, site `0x107b57`) when it recomputes `x_B` after a factorization:
    /// copies and scans of length `m`, the columns of the nonbasic structurals away from zero, and every fifth time an extra copy.
    fn primal_full_work(&mut self) {
        if self.mirror.is_none() || self.p0.is_some() {
            return;
        }
        let (n, m) = (self.p.n, self.p.m);
        let mut est = (10 * m + 2 * n) as f64;
        for j in 0..n {
            if self.pos_of[j] == NONE && self.nonbasic_value(j) != 0.0 {
                est += 2.0 * (self.p.cols[j].len() + 1) as f64;
            }
        }
        if self.refactorizations % 5 == 4 {
            est += 2.0 * m as f64;
        }
        // The right-hand side is built as a copy (2 units per row), or as `b + shift` (3 units) once the shift array exists.
        if !self.rhs_shift.is_empty() {
            est += m as f64;
        }
        self.work.add(0x107b57, est, crate::work::CTX);
    }

    /// The reference's handling of phase 1 cost changes after a pivot (`0x105f20`): a BTRAN of the cost changes `dc`
    /// (positions and values), the row product of the result, and a list update over that row, which also carries the
    /// entering column when it is structural.
    fn cost_change_work(
        &mut self,
        dc: &[(usize, f64)],
        q: usize,
        leave: usize,
        alpha_rq: f64,
        q_only: bool,
    ) {
        if q_only {
            // Only the entering variable changed class: no BTRAN and no row product. The reference applies `d_j -= dc * arj` over the
            // previous pivot row (`0x1084a0` again, with `dc` as the step) and repeats its list update (`0x107e14`).
            let step = dc[0].1 / alpha_rq;
            for &(j, arj) in &self.last_row.clone() {
                if j == q {
                    // The entering variable carries its own cost change: `(d - y) + dc` keeps the noise of `d`.
                    self.d[j] = (self.d[j] - step * arj) + dc[0].1;
                } else if self.pos_of[j] == NONE && !self.is_fixed(j)
                    || j >= self.p.n && self.pos_of[j] != NONE
                {
                    self.d[j] -= step * arj;
                }
            }
            if let Some((c, e)) = self.last_upd.clone() {
                self.update_list(&c, e, 0x107e14);
            }
            return;
        }
        // The entering variable's own change next to others: the reference first repeats the previous pivot row for it alone
        // (`0x107e14` pass), then solves for the remaining changes.
        if q != NONE && self.pos_of[q] != NONE && dc.len() > 1 {
            let r = self.pos_of[q];
            if let Some(k) = dc.iter().position(|e| e.0 == r) {
                self.cost_change_work(&[dc[k]], q, leave, alpha_rq, true);
                let rest: Vec<(usize, f64)> = dc.iter().copied().filter(|e| e.0 != r).collect();
                self.cost_change_work(&rest, q, leave, alpha_rq, false);
                return;
            }
        }
        let Some(mi) = self.mirror.as_mut() else {
            return;
        };
        let (n, m) = (self.p.n, self.p.m);
        let dc_sorted = dc.to_vec();
        let out = mi.btran(&dc_sorted, false, false);
        self.work
            .add_from(0x12de50, 0x10804b, mi.work, crate::work::FAC);
        // The mirror's values, which carry the reference's cancellation noise; a dense result is compacted with the row
        // product's drop tolerance.
        let bt_dense = mi.bt_dense;
        let mut rho = vec![0.0; m];
        for &(i, v) in &out {
            if !bt_dense || v.abs() > ROW_DROP_TOL {
                rho[i] = v;
            }
        }
        self.row_product_work(&rho, q, leave, 0x10820b);
        // The row's slack part is `rho` itself, basic slacks included (they are the positions whose cost changed). The
        // reference subtracts the row product from `d` in place, summing in the order of the BTRAN output list.
        let p = self.p;
        let mut rank = vec![usize::MAX; m];
        for (k, &(i, _)) in out.iter().enumerate() {
            rank[i] = k;
        }
        // The reference's pass also reaches basic entries of the stored d. The entering variable (basic by now) loses its
        // row product when it is a slack, `alpha_rq` when structural (it is the entering column's pivot element, as in
        // the pivot-row update). Basic slacks outside the cost-change positions lose their raw row entry (rounding noise);
        // basic structurals are not touched.
        let prods = p.rank_products(&rho, &rank);
        if q != NONE && !q_only {
            self.d[q] = if q < n {
                self.d[q] - alpha_rq
            } else {
                self.d[q] - prods[q]
            };
        }
        for j in n..n + m {
            if self.pos_of[j] != NONE && j != q {
                let ab = dot(&p.cols[j], &rho);
                // A basic slack at a cost-change position gets the row entry and its own cost change: `(d - y) + dc`, whose
                // rounding keeps the stored noise of `d` (`0x1084a0` on the row vector, then the cost change).
                if let Some(e) = dc.iter().find(|e| self.basic[e.0] == j) {
                    self.d[j] = (self.d[j] - ab) + e.1;
                } else if ab != 0.0 {
                    self.d[j] -= ab;
                }
            }
        }
        let mut changed: Vec<usize> = Vec::new();
        for j in 0..n + m {
            let a = if self.pos_of[j] == NONE {
                prods[j]
            } else {
                0.0
            };
            if self.pos_of[j] == NONE && a.abs() > ROW_DROP_TOL {
                self.d[j] -= a;
            }
            if ((j >= n && rho[j - n] != 0.0) || (self.pos_of[j] == NONE && a.abs() > ROW_DROP_TOL))
                && self.excluded != Some(j)
            {
                changed.push(j);
            }
        }
        let nz: Vec<usize> = out.iter().map(|e| e.0).filter(|&i| rho[i] != 0.0).collect();
        if crate::dump::on("ord") {
            eprintln!("CORD {}", serde_json::json!({"nz": nz, "changed": changed}));
        }
        let changed = self.pivot_row_order_by(&nz, &changed, nz.len());
        self.update_list(&changed, usize::from(q < n), 0x108280);
    }

    /// Work of the reference's dual recomputation after a factorization (`0x108d70`): a BTRAN of the dense vector of basic costs,
    /// with the phase 1 penalty on the infeasible ones.
    fn dual_btran_work(&mut self) {
        if self.mirror.is_none() {
            return;
        }
        let p = self.p;
        // The classes come from the distance array (`dst`), as in the reference's pricing0 rebuild (`0xfb790`).
        let sg = self.sig();
        let mi = self.mirror.as_mut().unwrap();
        let inp: Vec<(usize, f64)> = self
            .basic
            .iter()
            .enumerate()
            .map(|(pos, &j)| (pos, p.cost[j] + PHASE1_M * f64::from(sg[pos])))
            .filter(|e| e.1 != 0.0)
            .collect();
        let out = mi.btran(&inp, false, true);
        if crate::dump::on("dual") {
            crate::dump::emit(
                "dual",
                serde_json::json!({"refactor": self.refactorizations, "out": crate::dump::sparse(&out)}),
            );
        }
        let caller = if self.cleanup { 0xe53a7 } else { 0x109289 };
        self.work
            .add_from(0x12de50, caller, mi.work, crate::work::FAC);
        let mut y = vec![0.0; p.m];
        for &(i, v) in &out {
            y[i] = if self.cleanup && v.abs() < 1e-13 {
                0.0
            } else {
                v
            };
        }
        let mut rk = vec![usize::MAX; p.m];
        for (k, &(i, _)) in out.iter().enumerate() {
            rk[i] = k;
        }
        self.y_rank = rk;
        self.y_ref = Some(y);
        // The reduced costs that follow (`0x1020a0`, site `0x1021b0`): a pass over all variables and the columns of the
        // nonbasic structurals that are not fixed.
        if self.cleanup || self.p0.is_some() {
            return;
        }
        self.reduced_cost_work();
        self.residual_check_work();
    }

    /// Work of the reduced-cost pass after the duals (`0x1020a0`, site `0x1021b0`).
    fn reduced_cost_work(&mut self) {
        let p = self.p;
        let cols: usize = (0..p.n)
            .filter(|&j| self.pos_of[j] == NONE && !self.is_fixed(j) && self.excluded != Some(j))
            .map(|j| p.cols[j].len())
            .sum();
        self.work
            .add(0x1021b0, 3.0 * (p.n + p.m + cols) as f64, crate::work::CTX);
    }

    fn recompute_primal(&mut self) {
        self.recompute_primal_x();
        for pos in 0..self.basic.len() {
            let j = self.basic[pos];
            self.dst[j] = self.dist_of(j, self.x[j]);
        }
    }

    fn dist_of(&self, j: usize, x: f64) -> f64 {
        let p = self.p;
        if is_inf(p.lo[j]) {
            -x + p.up[j]
        } else {
            x - p.lo[j]
        }
    }

    fn recompute_primal_x(&mut self) {
        let p = self.p;
        let mut r = p.rhs.clone();
        for (ri, sh) in r.iter_mut().zip(&self.rhs_shift) {
            *ri += sh;
        }
        for j in 0..p.n + p.m {
            if self.pos_of[j] == NONE {
                let v = self.nonbasic_value(j);
                self.x[j] = v;
                if v != 0.0 {
                    for &(i, a) in &p.cols[j] {
                        r[i] -= a * v;
                    }
                }
            }
        }
        if let Some(mi) = self.mirror.as_ref() {
            // The clean-up pass solves for the unshifted right-hand side.
            let mut r0 = r.clone();
            if self.cleanup {
                for (ri, sh) in r0.iter_mut().zip(&self.rhs_shift) {
                    *ri -= sh;
                }
            }
            let est = mi.ftran_dense_work(&r0);
            self.work
                .add_from(0x1219d6, 0x120673, est, crate::work::FAC);
        }
        if let Some(mi) = self.mirror.as_mut() {
            // The reference's FTRAN of the basic solution, value for value. The mirror's numeric factor is exact only where the
            // reference's elimination arithmetic is modelled, so a result that fails the residual test falls back to the solve
            // with the working factor.
            let rows: Vec<usize> = (0..p.m).filter(|&i| r[i] != 0.0).collect();
            let vals: Vec<f64> = rows.iter().map(|&i| r[i]).collect();
            let out = mi.ftran_basic(&rows, &vals);
            let mut xb = vec![0.0; p.m];
            for (i, v) in out {
                xb[i] = v;
            }
            if crate::dump::on("xb") {
                crate::dump::emit("xb", serde_json::json!({"rhs": r, "xb": xb}));
            }
            let mut res = r.clone();
            for (pos, &j) in self.basic.iter().enumerate() {
                for &(i, a) in &p.cols[j] {
                    res[i] -= a * xb[pos];
                }
            }
            let scale = 1.0 + r.iter().fold(0.0f64, |a, v| a.max(v.abs()));
            if res.iter().all(|v| v.abs() <= 1e-7 * scale) {
                for (pos, &j) in self.basic.iter().enumerate() {
                    self.x[j] = xb[pos];
                }
                return;
            }
        }
        self.factor.ftran(&mut r);
        for (pos, &j) in self.basic.iter().enumerate() {
            self.x[j] = r[pos];
        }
    }

    /// `d_j = c_j - sum y_i a_ij`, subtracting the terms in the column's stored order (`0x1020a0`).
    fn set_d_from_y(&mut self, y: &[f64]) {
        let p = self.p;
        let offs = &p.offs;
        let col_off = |j: usize| offs[j] & 1;
        for j in 0..p.n + p.m {
            self.d[j] = if self.pos_of[j] == NONE {
                if self.cleanup && p.lo[j] == p.up[j] {
                    // The API read-out (`0x1a7120..0x1a73c0`) recomputes the reduced cost of a fixed variable from the matrix: the cost plus
                    // the negated column products in stored column order, with the eight-lane kernel (peel by the value pointer's alignment).
                    dot_lanes_from(&p.fcols[j], y, col_off(j), p.cost[j])
                } else if self.cleanup {
                    // Every other nonbasic variable returns the solver's d, here the sum in the order of the rows in the dual BTRAN output list.
                    let mut e: Vec<(usize, f64)> = p.fcols[j].clone();
                    e.sort_by_key(|x| self.y_rank[x.0]);
                    let mut acc = p.cost[j];
                    for &(i, a) in &e {
                        acc -= y[i] * a;
                    }
                    acc
                } else if self.p0.is_some() {
                    // `pricing0`: the returned reduced costs are the scan's scores (`0xfecb0`), the eight-lane kernel.
                    dot_lanes_from(&p.fcols[j], y, col_off(j), p.cost[j])
                } else {
                    let mut acc = p.cost[j];
                    for &(i, a) in &p.fcols[j] {
                        acc -= y[i] * a;
                    }
                    acc
                }
            } else {
                0.0
            };
        }
    }

    fn recompute_dual(&mut self) {
        let p = self.p;
        if let Some(y) = self.y_ref.take() {
            self.set_d_from_y(&y);
            self.y_ref = Some(y);
            return;
        }
        let mut y: Vec<f64> = self.basic.iter().map(|&j| p.cost[j]).collect();
        self.factor.btran(&mut y);
        for j in 0..p.n + p.m {
            self.d[j] = if self.pos_of[j] == NONE {
                p.cost[j] - dot(&p.cols[j], &y)
            } else {
                0.0
            };
        }
    }

    fn is_fixed(&self, j: usize) -> bool {
        self.p.lo[j] == self.p.up[j]
    }

    /// Order in which the reference stores the pivot-row entries: with a sparse `rho`, structurals in first-touch
    /// order over the rows of `rho`; otherwise structurals in index order; then slacks in row order.
    fn pivot_row_order(&self, rho: &[f64], changed: &[usize]) -> Vec<usize> {
        let nz: Vec<usize> = (0..self.p.m).rev().filter(|&i| rho[i] != 0.0).collect();
        self.pivot_row_order_by(&nz, changed, nz.len())
    }

    /// As `pivot_row_order`, with the rows of `rho` in the order the reference stores them.
    fn pivot_row_order_by(&self, nz: &[usize], changed: &[usize], cnt: usize) -> Vec<usize> {
        let (n, m) = (self.p.n, self.p.m);
        let nnz_a: usize = self.p.cols[..n].iter().map(Vec::len).sum();
        let mut keep = vec![false; n + m];
        for &j in changed {
            keep[j] = true;
        }
        let mut out = Vec::with_capacity(changed.len());
        if cnt * 3 * nnz_a < m * n {
            for &i in nz {
                for &j in &self.rrows[i] {
                    if keep[j] {
                        keep[j] = false;
                        out.push(j);
                    }
                }
            }
        } else {
            for j in 0..n {
                if keep[j] {
                    keep[j] = false;
                    out.push(j);
                }
            }
        }
        for &i in nz {
            if keep[n + i] {
                keep[n + i] = false;
                out.push(n + i);
            }
        }
        out.extend((0..n + m).filter(|&j| keep[j]));
        out
    }

    /// Work of the reference's pivot-row product (`0x134d90`) for the row combination `rho` (`reverse-engineer/specs/primal-simplex.md`,
    /// "Work counter"): its rows hold the nonbasic structural columns that are not fixed, where `q` already counts as basic and
    /// `leave` as nonbasic.
    fn row_product_work(&mut self, rho: &[f64], q: usize, leave: usize, caller: u32) {
        let (n, m) = (self.p.n, self.p.m);
        let nz: Vec<usize> = (0..m).filter(|&i| rho[i] != 0.0).collect();
        let nnz_a: usize = self.p.cols[..n].iter().map(Vec::len).sum();
        let mut seen = vec![false; n];
        let (mut p5, mut d) = (0usize, 0usize);
        for &i in &nz {
            for &j in &self.rows[i] {
                if ((self.pos_of[j] == NONE && j != q) || j == leave)
                    && self.p.lo[j] != self.p.up[j]
                    && self.excluded != Some(j)
                {
                    p5 += 1;
                    if !seen[j] {
                        seen[j] = true;
                        d += 1;
                    }
                }
            }
        }
        let nr = nz.len() as f64;
        let est = if nz.len() == 1 {
            // A multiplier below one in magnitude goes through the filtering loop, one more pass over the row.
            if rho[nz[0]].abs() < 1.0 {
                10.5 + 5.0 * p5 as f64
            } else {
                10.5 + 4.0 * p5 as f64
            }
        } else if m * n <= nz.len() * 3 * nnz_a {
            8.5 * nr + 2.0 * p5 as f64 + 2.5 * n as f64
        } else {
            11.5 * nr + 3.5 * p5 as f64 + 3.0 * d as f64
        };
        self.work.add_from(0x135e92, caller, est, crate::work::CTX);
    }

    /// Eligibility and score `d^2 / w` of a nonbasic variable, or `None` if it cannot enter.
    fn score(&self, j: usize) -> Option<f64> {
        if self.pos_of[j] != NONE || self.is_fixed(j) || self.excluded == Some(j) {
            return None;
        }
        // The reference's candidate scan accepts |d| above 1e-6 (it rejected a d of 1.5e-7 that the default 1e-7 takes).
        let tol = if self.opts.faithful_weights && self.p1_tol_on {
            1e-7 * PHASE1_M
        } else if self.opts.faithful_weights {
            self.opts.optimality_tol.max(1e-6)
        } else {
            self.opts.optimality_tol
        };
        let free = is_inf(self.p.lo[j]) && is_inf(self.p.up[j]);
        let dj = self.d[j];
        if free && self.opts.faithful_weights {
            // The reference lists nonbasic free variables (status -3) in the top bucket and takes the first one whatever its reduced cost.
            return Some(f64::INFINITY);
        }
        let eligible = if free {
            dj.abs() > tol
        } else if self.at_upper[j] {
            dj > tol
        } else {
            dj < -tol
        };
        eligible.then(|| dj * dj / self.w[j])
    }

    /// Entering variable: the first maximum-score entry of the persistent candidate list, rebuilt by a full
    /// pass over all variables when the list holds no eligible entry.
    fn price(&mut self) -> Option<usize> {
        if !self.list.is_empty() {
            let best = self.best_in_list();
            self.work
                .add(0x105d03, self.scanned(best) as f64 * 4.0, crate::work::CTX);
            if best.is_some() {
                return best;
            }
            if self.mirror.is_some() {
                // The reference does not rebuild the list from here: the caller recomputes the duals first.
                self.list_failed = true;
                return None;
            }
        }
        self.level = 0;
        self.list.clear();
        let nm = self.p.n + self.p.m;
        let full_pass = nm as f64 + nm as f64;
        if !self.rescan {
            self.work.add(0x105b2d, full_pass, crate::work::CTX);
        }
        for j in 0..nm {
            if let Some(sc) = self.score(j) {
                if thr(self.level) < sc {
                    let b = bucket(sc);
                    if b > self.level {
                        self.level = b;
                        self.list.clear();
                    }
                    self.list.push(j);
                }
            }
        }
        let best = self.best_in_list();
        if self.rescan && best.is_some() {
            self.work.add(0x105b2d, full_pass, crate::work::CTX);
        }
        if !self.list.is_empty() {
            self.work
                .add(0x105c0c, self.scanned(best) as f64 * 4.0, crate::work::CTX);
        }
        best
    }

    /// `price`, with the reference's recovery when the list scan finds nothing: the duals are recomputed from scratch
    /// (no pivot) and the list is rebuilt by the next pass.
    fn price_faithful(&mut self) -> Option<usize> {
        self.list_failed = false;
        self.recomputed = false;
        let q = self.price();
        if q.is_none() && self.list_failed {
            self.list_failed = false;
            self.recomputed = true;
            self.dual_btran_work();
            self.recompute_dual();
            self.list.clear();
            self.level = 0;
            if crate::dump::on("pricing") {
                crate::dump::emit(
                    "pricing",
                    serde_json::json!({"n": PRC.load(std::sync::atomic::Ordering::Relaxed), "iter": 0, "infeasible": true, "list": self.list, "basic": self.basic, "d": self.d, "w": self.w, "x": self.x}),
                );
            }
            PRC.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
            self.rescan = true;
            let q = self.price();
            self.rescan = false;
            return q;
        }
        q
    }

    /// `price_faithful`, except that the variable dropped from a singular supplied basis stays out of the scan until the scan finds no
    /// other candidate; it is then admitted and the list rebuilt (the reference adds 3 to its status, -4 -> -1, in `0xcda5b`).
    fn price_admitting(&mut self) -> Option<usize> {
        let w0 = self.work.w;
        let q = self.price_faithful();
        if q.is_none() && self.excluded.is_some() {
            // Without updates the scan books nothing (the driver's own list build, `0x138532` below, takes its place).
            if self.factor.num_updates() == 0 {
                self.work.rollback(w0);
            }
            // The scan found nothing: like the optimality exit, the reference first refactors when updates were made (`0xd4190`,
            // the factorization before the status add of `0x155139`, the excluded variable still out of the row copy); only then
            // does the driver restart below.
            if self.factor.num_updates() > 0 {
                // As at the optimality exit: the duals are recomputed once more when the scan finds nothing (`0x145520`).
                if self.mirror.is_some() && !self.recomputed && !self.fresh_duals {
                    self.dual_btran_work();
                }
                self.work.trig();
                self.work.decided(self.work.w);
                self.pre_restart = true;
                let r = self.refactor();
                self.pre_restart = false;
                r.ok()?;
            }
            self.excluded = None;
            // Admitting the variable rebuilds the candidate list (`0x138532`: one pass over all variables) and the row-wise copy
            // (`0x132d88`: `2 nnz(A)`), before the driver restarts.
            let nnz_a: usize = (0..self.p.n).map(|j| self.p.cols[j].len()).sum();
            let nm = self.p.n + self.p.m;
            self.work.add(0x138532, 2.0 * nm as f64, crate::work::CTX);
            self.work
                .add(0x132d88, 2.0 * nnz_a as f64, crate::work::CTX);
            // The rebuilt row-wise copy holds every nonbasic, non-fixed structural (the admitted one included), in index order.
            self.rrows = vec![Vec::new(); self.p.m];
            for j in 0..self.p.n {
                if self.pos_of[j] == NONE && !self.is_fixed(j) {
                    for &(i, _) in &self.p.cols[j] {
                        self.rrows[i].push(j);
                    }
                }
            }
            self.work.restart();
            // The reference re-enters its driver, which rebuilds the basis in ascending column order and refactors.
            self.restart = true;
            let r = self.refactor();
            self.restart = false;
            self.restarted = true;
            if r.is_err() {
                return None;
            }
            self.list.clear();
            self.level = 0;
            let w1 = self.work.w;
            let q = self.price_faithful();
            if q.is_none() {
                self.work.rollback(w1);
            }
            return q;
        }
        q
    }

    /// The reference's "Resetting Norm" (`0x102d20`, after the basis change): w = 1 and framework membership for every nonbasic
    /// variable, w = 0 and no membership for the basic ones. The candidate list keeps its entries and takes the normal update.
    fn reset_weights(&mut self) {
        let nm = self.p.n + self.p.m;
        for j in 0..nm {
            let nb = self.pos_of[j] == NONE;
            self.w[j] = if nb { 1.0 } else { 0.0 };
            self.in_r0[j] = nb;
        }
    }

    /// Entries of the list the reference's scan looks at: all of them, except that it stops at a free variable (top bucket,
    /// nothing beats it), after the entries before it.
    fn scanned(&self, best: Option<usize>) -> usize {
        match best {
            Some(j) if self.score(j) == Some(f64::INFINITY) => {
                self.list.iter().position(|&x| x == j).unwrap()
            }
            _ => self.list.len(),
        }
    }

    fn best_in_list(&self) -> Option<usize> {
        let (mut best, mut best_score) = (None, 0.0);
        for &j in &self.list {
            if let Some(sc) = self.score(j) {
                if best_score < sc {
                    best = Some(j);
                    best_score = sc;
                }
            }
        }
        best
    }

    /// After a pivot: score the changed entries (in index order), raise the level if a higher bucket appears,
    /// keep still-qualifying old entries first, then append the new ones.
    fn update_list(&mut self, changed: &[usize], extra: usize, caller: u32) {
        if caller == 0x1442ab {
            self.last_upd = Some((changed.to_vec(), extra));
            self.last_row = self.pending_row.take().unwrap_or_default();
        }
        let mut fresh: Vec<usize> = Vec::new();
        // The reference empties the old list when the level rises, so it scans none of it.
        let mut raised = false;
        for &j in changed {
            // The reference's update scores only the entries with status -1 or -2; a free variable (-3) keeps its place.
            if self.opts.faithful_weights && is_inf(self.p.lo[j]) && is_inf(self.p.up[j]) {
                continue;
            }
            if let Some(sc) = self.score(j) {
                if thr(self.level) < sc {
                    let b = bucket(sc);
                    if self.level < b {
                        self.level = b;
                        fresh.clear();
                        raised = true;
                    }
                    fresh.push(j);
                }
            }
        }
        let old = std::mem::take(&mut self.list);
        let old_len = if raised { 0 } else { old.len() };
        let mut kept: Vec<usize> = Vec::new();
        for j in old {
            if fresh.contains(&j) {
                continue;
            }
            if let Some(sc) = self.score(j) {
                if thr(self.level) < sc {
                    kept.push(j);
                }
            }
        }
        kept.extend(fresh);
        self.list = kept;
        // The reference's pivot row has no fixed structural columns (its row storage leaves them out).
        let n = self.p.n;
        let rowlen = changed
            .iter()
            .filter(|&&j| j >= n || !self.is_fixed(j))
            .count();
        let est = (rowlen + extra) as f64 * 4.0
            + old_len as f64
            + old_len as f64
            + self.list.len() as f64 * 2.0;
        if crate::dump::on("work") {
            crate::dump::emit(
                "work",
                serde_json::json!({"k": "upd", "caller": caller, "changed": changed, "old": old_len, "new": self.list, "level": self.level, "est": est, "extra": extra}),
            );
        }
        self.work.add_from(0x1087b1, caller, est, crate::work::CTX);
    }
}

/// Bucket key of the `pricing0` candidate list (the sort routine `0xfffc0`, called by `0xfecb0` when `P+0x80 == 1`, as on every
/// corpus model): the binary exponent `k` of `|s|` (`frexp`), clamped to -25..=25 (25 at `|s| >= 2^25`, -25 below 2^-25 plus an ulp);
/// the key is `0x4b + k` for an attractive (negative) score and `0x19 - k` for the others, so a score of 2^24 or more ranks with
/// the basic and fixed variables (key 0). With `m >= 500` (`P+0x1c = n`, `0x137ce0`) a structural variable's `|s|` is first divided by
/// its column count plus 0.01 (`div`, `ctx+0x78` = the column's entry count).
fn p0_key(s: f64, div: Option<f64>) -> usize {
    const HI: f64 = 33554432.0; // 2^25
    const LO: f64 = f64::from_bits(0x3e60000000000001);
    let a = match div {
        Some(c) => s.abs() / (c + 0.01),
        None => s.abs(),
    };
    let k: i32 = if a >= HI {
        25
    } else if a < LO {
        -25
    } else {
        ((a.to_bits() >> 52) & 0x7ff) as i32 - 1022
    };
    (if s < 0.0 { k + 0x4b } else { 0x19 - k }) as usize
}

/// Whether the reference runs the bucket-chain candidate list (`P+0x80 = 2`, `P+0x78 = r8d` at `0x1385d1`): with `SimplexPricing` at its release
/// default and no warm start, when `trunc(n * (n / m) / max(n, nnz)) >= 2` (`0x137f86`..`0x137fce`, `nnz` = `ctx+0x3d0`). That is a model of
/// one column count per row at most: `fx_unbounded`, `fx_tied` (1 row, 2 columns); no Netlib model of the corpus below 80BAU3B.
fn p0_bucket_mode(n: usize, m: usize, nnz: usize) -> bool {
    p0_bucket_r8(n, m, nnz) >= 2
}

/// `r8d` of `0x137fc9`: `trunc(n * (n / m) / max(n, nnz))`.
fn p0_bucket_r8(n: usize, m: usize, nnz: usize) -> usize {
    if m == 0 {
        return 0;
    }
    ((n as f64 / m as f64 * n as f64) / n.max(nnz) as f64) as usize
}

/// Bucket key of the bucket-chain list (`0xfecb0`, `0xfbe60`): like `p0_key` with the exponent clamped to -15..=34 (`17179500000` and
/// `3.051758e-05`), `0x41 + k` for an attractive score and `0x23 - k` for the others. Keys run from 1 to 99.
fn p0_key2(s: f64, div: Option<f64>) -> usize {
    const HI: f64 = 17179500000.0;
    const LO: f64 = 3.051758e-05;
    let a = match div {
        Some(c) => s.abs() / (c + 0.01),
        None => s.abs(),
    };
    let k: i32 = if a >= HI {
        34
    } else if a < LO {
        -15
    } else {
        ((a.to_bits() >> 52) & 0x7ff) as i32 - 1022
    };
    (if s < 0.0 { 0x41 + k } else { 0x23 - k }) as usize
}

/// `PerturbValue` default (`ctx+0x340` at the start).
const PERTURB_VALUE: f64 = 0.0002;

/// Initial window of the `pricing0` scan: it ends at the first position more than this many past its start once a candidate was counted (`P+0x50`).
const P0_WINDOW: usize = 10;
/// Reduced-cost tolerance of the `pricing0` scan (`ctx+0x388`).
const P0_DTOL: f64 = 1e-6;
/// Phase 1 form of the tolerance (`ctx+0x2d0 != 0`): `P+0x3c0` (1e-7) times the infeasibility cost shift `P[0]` (1e6).
const P0_DTOL_P1: f64 = 1e-7 * PHASE1_M;

impl<'a> Solver<'a> {
    /// Nonbasic, not fixed: the variables the `pricing0` list ranks (status -1, -2 or -3 of the reference).
    fn p0_eligible(&self, j: usize) -> bool {
        self.pos_of[j] == NONE && !self.is_fixed(j)
    }

    /// The reference's score `s` of a nonbasic variable (`d`, negated at the upper bound) and whether it is free, from the duals `y`.
    fn p0_score(&self, j: usize, y: &[f64], col_off: &dyn Fn(usize) -> usize) -> (f64, bool) {
        let p = self.p;
        let d = dot_lanes_from(&p.fcols[j], y, col_off(j), p.cost[j]);
        let free = is_inf(p.lo[j]) && is_inf(p.up[j]);
        if free {
            (d, true)
        } else if self.at_upper[j] {
            (-d, false)
        } else {
            (d, false)
        }
    }

    fn p0_col_offsets(&self) -> &'a [usize] {
        &self.p.offs
    }

    /// Recomputes the duals (a dense BTRAN of the basic costs), then scores every variable and sorts the list (`0xfb790`, `0xfecb0`).
    fn p0_rebuild(&mut self) {
        if self.y_ref.is_none() {
            self.dual_btran_work();
            // `0x145520` runs its residual check after the duals (the refactorization's own run is booked in `refactor`).
            self.residual_check_work();
        }
        let y = self.y_ref.clone().unwrap();
        let nm = self.p.n + self.p.m;
        let offs = self.p0_col_offsets();
        let col_off = |j: usize| offs[j] & 1;
        let mut dstore = vec![0.0; nm];
        let mut keys = vec![0usize; nm];
        // Work of the sort routine `0xfffc0`: the bucket table (101), 2 per column entry of every eligible structural, then 4 per
        // structural, 5 per row and 3 per variable (the count and placement passes).
        let mut sort_work: Vec<(u32, f64)> = vec![(0x100fa5, 101.0)];
        for j in 0..nm {
            if !self.p0_eligible(j) {
                continue;
            }
            let (s, free) = self.p0_score(j, &y, &col_off);
            if j < self.p.n {
                sort_work.push((0x1010f3, 2.0 * self.p.cols[j].len() as f64));
            }
            dstore[j] = s;
            let div = if self.p.m >= 500 && j < self.p.n {
                Some(self.p.cols[j].len() as f64)
            } else {
                None
            };
            keys[j] = if free { 100 } else { p0_key(s, div) };
        }
        if self.p0.as_ref().unwrap().bk.is_some() {
            // Bucket chains: the variables of each key in index order (`0xfecb0`, the loop that appends to `P+0x21c + 4 key`).
            // Work of the bucket rebuild (`0xfecb0`): the score of every eligible structural books `3 len` (`0xffd5b`, in index
            // order), then `7 n` (`0xffab4`), `4 m` (`0xff07e`) and `5 (n + m)` (`0xff207`).
            for j in 0..self.p.n {
                if self.p0_eligible(j) && !self.p.cols[j].is_empty() {
                    self.work
                        .add(0xffd5b, 3.0 * self.p.cols[j].len() as f64, crate::work::CTX);
                }
            }
            self.work
                .add(0xffab4, 7.0 * self.p.n as f64, crate::work::CTX);
            self.work
                .add(0xff07e, 4.0 * self.p.m as f64, crate::work::CTX);
            self.work.add(0xff207, 5.0 * nm as f64, crate::work::CTX);
            let keys2: Vec<usize> = (0..nm)
                .map(|j| {
                    if !self.p0_eligible(j) {
                        return 0;
                    }
                    let div = if self.p.m >= 500 && j < self.p.n {
                        Some(self.p.cols[j].len() as f64)
                    } else {
                        None
                    };
                    p0_key2(dstore[j], div)
                })
                .collect();
            let p0 = self.p0.as_mut().unwrap();
            let bk = p0.bk.as_mut().unwrap();
            bk.head = [-1; 101];
            bk.tail = [-1; 101];
            bk.visited = 0;
            bk.built = true;
            bk.full = false;
            for j in 0..nm {
                let k = keys2[j];
                if k != 0 {
                    bk.append(k, j);
                }
            }
            p0.dstore = dstore;
            p0.y = y;
            p0.pos = 0;
            p0.rebuild = false;
            p0.pending = None;
            return;
        }
        sort_work.push((0x1010a3, 4.0 * self.p.n as f64));
        sort_work.push((0x1002f8, 5.0 * self.p.m as f64));
        sort_work.push((0x100404, 3.0 * nm as f64));
        for (site, est) in sort_work {
            self.work.add(site, est, crate::work::CTX);
        }
        let mut list: Vec<usize> = (0..nm).collect();
        list.sort_by(|&a, &b| keys[b].cmp(&keys[a]).then(b.cmp(&a)));
        let p0 = self.p0.as_mut().unwrap();
        p0.list = list;
        p0.dstore = dstore;
        p0.y = y;
        p0.pos = 0;
        p0.rebuild = false;
        p0.pending = None;
    }

    /// Phase 1 of `pricing0`: the cost shifts of the basic variables that changed feasibility class since the last pricing call join the
    /// BTRAN input of the next dual update (after the entering `d`). Without a valid dual (`update` false: a fresh start, a refactorization
    /// or the first infeasible pass) the duals are rebuilt from the costs and the shifts instead.
    fn p0_cost_changes(&mut self, update: bool, dc: &[(usize, f64)]) {
        let p0 = self.p0.as_mut().unwrap();
        if !update {
            p0.rebuild = true;
            return;
        }
        if !dc.is_empty() {
            p0.pending
                .get_or_insert_with(Vec::new)
                .extend_from_slice(dc);
        }
    }

    /// Applies the pending dual update: the BTRAN of the entering `d` at the leaving position, in the updated basis (`0xfb790`, mode 2).
    fn p0_update_duals(&mut self) {
        let Some(inp) = self.p0.as_mut().unwrap().pending.take() else {
            return;
        };
        let mi = self.mirror.as_mut().unwrap();
        let out = mi.btran(&inp, false, false);
        let dense = mi.bt_dense;
        if crate::dump::on("dual") {
            crate::dump::emit(
                "dual",
                serde_json::json!({"update": inp, "out": crate::dump::sparse(&out)}),
            );
        }
        self.work
            .add_from(0x12de50, 0xfb7ee, mi.work, crate::work::FAC);
        let y = &mut self.p0.as_mut().unwrap().y;
        for &(i, v) in &out {
            if !dense || v < -ROW_DROP_TOL || ROW_DROP_TOL < v {
                y[i] += v;
            }
        }
    }

    /// One scan of the list window from the current position (`0xfbe60`). Returns the entering variable and its `d`.
    fn p0_scan(&mut self) -> Option<(usize, f64)> {
        if self.p0.as_ref().unwrap().bk.is_some() {
            return self.p0_scan2();
        }
        let m = self.p.m;
        let nm = self.p.n + self.p.m;
        let offs = self.p0_col_offsets();
        let col_off = |j: usize| offs[j] & 1;
        let mut p0 = self.p0.take().unwrap();
        let start = p0.pos;
        let (mut i, mut cnt, mut best, mut q): (usize, usize, f64, Option<(usize, f64)>) =
            (start, 0, 0.0, None);
        let mut marked = vec![false; nm];
        // Work of the reference's scan (`0xfdef9`): 2 per entry visited, and per candidate 6 plus 2 per column entry of a structural.
        let mut scan_est = 0.0f64;
        if start < m {
            loop {
                let j = p0.list[i];
                scan_est += 2.0;
                if self.p0_eligible(j) {
                    scan_est += 6.0;
                    if j < self.p.n {
                        scan_est += 2.0 * self.p.cols[j].len() as f64;
                    }
                    let (s, free) = self.p0_score(j, &p0.y, &col_off);
                    if free {
                        // A free variable is taken at once, whatever its reduced cost.
                        p0.dstore[j] = if s <= 0.0 { s } else { -s };
                        q = Some((j, s));
                        break;
                    }
                    if s < -p0.tol {
                        if s + s < p0.dstore[j] {
                            cnt += 1;
                            marked[i] = true;
                            // With `m >= 500` (`P+0x1c = n`) a structural's score is divided by its column count plus 0.01 for the comparison.
                            let sc = if m >= 500 && j < self.p.n {
                                s / (self.p.cols[j].len() as f64 + 0.01)
                            } else {
                                s
                            };
                            if best > sc {
                                best = sc;
                                q = Some((j, if self.at_upper[j] { -s } else { s }));
                            }
                        }
                        p0.dstore[j] = s;
                    }
                }
                if p0.window < i - start && cnt != 0 {
                    break;
                }
                i += 1;
                if i >= m {
                    break;
                }
            }
        }
        let newpos = i + 1 - cnt;
        if self.mirror.is_some() {
            if cnt != 0 {
                scan_est += 4.0 * cnt as f64;
            }
            if scan_est != 0.0 {
                self.work.add(0xfdef9, scan_est, crate::work::CTX);
            }
        }
        // Gather the counted entries into the last `cnt` positions of the window, in place.
        let mut h = start;
        let mut t = newpos;
        while t <= i {
            if marked[t] {
                marked[t] = false;
            } else {
                while !marked[h] {
                    h += 1;
                }
                p0.list.swap(h, t);
                marked[h] = false;
                h += 1;
            }
            t += 1;
        }
        p0.pos = newpos;
        self.p0 = Some(p0);
        q
    }

    /// The scan of the bucket-chain list (`0xfbe60` from `0xfbfcc`, non-sparse loop `0xfce50`): walks the chains from the highest non-empty
    /// bucket down through the visited entries (basic ones are only linked), scoring every nonbasic non-fixed one (the stored `d` takes the
    /// score whether attractive or not); the window ends after `window` entries once an attractive one was seen. The visited entries are
    /// then unlinked from their buckets (the buckets above the end are emptied, the partly scanned one keeps the rest; every `window`-th
    /// update the rest and the lower buckets move up one key) and re-inserted by their new score, except the entering variable.
    fn p0_scan2(&mut self) -> Option<(usize, f64)> {
        let nm = self.p.n + self.p.m;
        let offs = self.p0_col_offsets();
        let col_off = |j: usize| offs[j] & 1;
        let upd = self.mirror.as_ref().map_or(0, |mi| mi.updates());
        let mut p0 = self.p0.take().unwrap();
        let mut bk = p0.bk.take().unwrap();
        let m500 = self.p.m >= 500;
        let n = self.p.n;
        if bk.visited as f64 > 3.0 * nm as f64 {
            // `P+0x3c8` above three times the variable count: nothing is returned, the caller recomputes and rebuilds.
            p0.bk = Some(bk);
            self.p0 = Some(p0);
            return None;
        }
        let mut b0 = 0usize;
        if bk.head[100] >= 0 {
            b0 = 100;
        } else {
            for k in (1..100).rev() {
                if bk.head[k] >= 0 {
                    b0 = k;
                    break;
                }
            }
        }
        let (mut prev, mut first): (i64, i64) = (-1, -1);
        let mut last: i64 = -1;
        let mut q: Option<(usize, f64)> = None;
        let mut best = 0.0f64;
        let mut nattr = 0usize;
        let mut visited = 0usize;
        let mut b = b0;
        // Work of the scan (`0xfc693`): 100, 5 per step of the walk (an entry or the end of a chain), `2 len` per eligible structural and 6
        // per entry of the visited chain in the relink pass.
        let mut scan_est = 100.0f64;
        if b0 > 0 {
            let mut cur = bk.head[b];
            'scan: loop {
                scan_est += 5.0;
                if cur < 0 {
                    loop {
                        b -= 1;
                        if b == 0 {
                            last = -1;
                            break 'scan;
                        }
                        cur = bk.head[b];
                        if cur >= 0 {
                            break;
                        }
                    }
                }
                let j = cur as usize;
                visited += 1;
                if self.p0_eligible(j) {
                    if j < n {
                        scan_est += 2.0 * self.p.cols[j].len() as f64;
                    }
                    let (s, free) = self.p0_score(j, &p0.y, &col_off);
                    p0.dstore[j] = s;
                    if free {
                        // A free variable is taken at once and is not linked into the visited chain.
                        q = Some((j, s));
                        last = cur;
                        break 'scan;
                    }
                    if s < -p0.tol {
                        nattr += 1;
                        let sc = if m500 && j < n {
                            s / (self.p.cols[j].len() as f64 + 0.01)
                        } else {
                            s
                        };
                        if best > sc {
                            best = sc;
                            q = Some((j, if self.at_upper[j] { -s } else { s }));
                        }
                    }
                }
                if prev != -1 {
                    bk.next[prev as usize] = cur;
                } else {
                    first = cur;
                }
                prev = cur;
                last = cur;
                if visited > p0.window && nattr != 0 {
                    break 'scan;
                }
                cur = bk.next[j];
            }
        }
        bk.visited += visited;
        // Unlink: the buckets above the end of the scan are empty, the one it ended in keeps the rest of its chain.
        if b == 100 {
            let rest = if last >= 0 {
                bk.next[last as usize]
            } else {
                -1
            };
            bk.head[100] = rest;
            if rest < 0 {
                bk.tail[100] = -1;
            }
        } else {
            for k in b + 1..=b0.min(99) {
                bk.head[k] = -1;
                bk.tail[k] = -1;
            }
            if b != 0 {
                let rest = if last >= 0 {
                    bk.next[last as usize]
                } else {
                    -1
                };
                let w = p0.window.max(1);
                if upd % w != w - 1 {
                    bk.head[b + 1] = -1;
                    bk.tail[b + 1] = -1;
                    bk.head[b] = rest;
                    if rest < 0 {
                        bk.tail[b] = -1;
                    }
                } else {
                    if rest < 0 {
                        bk.head[b + 1] = -1;
                        bk.tail[b + 1] = -1;
                    } else {
                        bk.head[b + 1] = rest;
                        bk.tail[b + 1] = bk.tail[b];
                    }
                    for k in (2..=b).rev() {
                        bk.head[k] = bk.head[k - 1];
                        bk.tail[k] = bk.tail[k - 1];
                    }
                    bk.head[1] = -1;
                    bk.tail[1] = -1;
                }
            } else {
                bk.head[1] = -1;
                bk.tail[1] = -1;
            }
        }
        if prev >= 0 {
            bk.next[prev as usize] = -1;
        }
        // Re-insert the visited entries by their new score.
        let mut cur = first;
        while cur >= 0 {
            scan_est += 6.0;
            let j = cur as usize;
            let nxt = bk.next[j];
            bk.next[j] = -1;
            if q.map(|e| e.0) != Some(j) && self.p0_eligible(j) {
                let free = is_inf(self.p.lo[j]) && is_inf(self.p.up[j]);
                if !free {
                    let div = if m500 && j < n {
                        Some(self.p.cols[j].len() as f64)
                    } else {
                        None
                    };
                    bk.append(p0_key2(p0.dstore[j], div), j);
                }
            }
            cur = nxt;
        }
        p0.bk = Some(bk);
        self.p0 = Some(p0);
        if self.mirror.is_some() {
            self.work.add(0xfc693, scan_est, crate::work::CTX);
        }
        q
    }

    /// Entering variable of the `pricing0` option, or `None` at optimality (after the recomputation of the duals the reference makes
    /// when the list holds nothing).
    fn price0(&mut self) -> Option<usize> {
        self.recomputed = false;
        if self.p0.as_ref().unwrap().rebuild {
            if self
                .p0
                .as_ref()
                .unwrap()
                .bk
                .as_ref()
                .is_some_and(|b| b.built && !b.full)
            {
                // Bucket chains keep their state across a refactorization: only the duals are recomputed (`P+0x10 = 0`, `0xfb790` in the scan).
                if self.y_ref.is_none() {
                    self.dual_btran_work();
                }
                let y = self.y_ref.clone().unwrap();
                let p0 = self.p0.as_mut().unwrap();
                p0.y = y;
                p0.rebuild = false;
                p0.pending = None;
            } else {
                self.p0_rebuild();
            }
        } else {
            self.p0_update_duals();
        }
        let dump = |s: &Self, tag: &str, q: Option<(usize, f64)>| {
            if crate::dump::on("pricing0") {
                let p0 = s.p0.as_ref().unwrap();
                crate::dump::emit(
                    "pricing0",
                    serde_json::json!({"n": p0.calls, "tag": tag, "list": p0.bk.as_ref().map_or(p0.list.iter().map(|&x| x as i64).collect::<Vec<_>>(), |b| b.next.clone()), "head": p0.bk.as_ref().map(|b| b.head.to_vec()), "tail": p0.bk.as_ref().map(|b| b.tail.to_vec()), "pos": p0.pos, "y": p0.y, "dstore": p0.dstore, "q": q.map(|e| e.0), "d": q.map(|e| e.1)}),
                );
            }
        };
        let before = |s: &Self| {
            if crate::dump::on("pricing0") {
                let p0 = s.p0.as_ref().unwrap();
                crate::dump::emit(
                    "pricing0",
                    serde_json::json!({"n": p0.calls, "tag": "in", "list": p0.bk.as_ref().map_or(p0.list.iter().map(|&x| x as i64).collect::<Vec<_>>(), |b| b.next.clone()), "head": p0.bk.as_ref().map(|b| b.head.to_vec()), "tail": p0.bk.as_ref().map(|b| b.tail.to_vec()), "pos": p0.pos, "y": p0.y, "dstore": p0.dstore}),
                );
            }
        };
        before(self);
        let mut q = self.p0_scan();
        dump(self, "out", q);
        self.p0.as_mut().unwrap().calls += 1;
        if q.is_none()
            && self
                .mirror
                .as_ref()
                .is_some_and(|mi| 5 * mi.updates() > 4 * mi.max_updates())
        {
            // The reference does not recompute the duals when the scan finds nothing after many updates (`0x1445ad`: five times
            // the updates above four times the limit): it sets the refactorization state 9 and prices again on the new factor.
            self.want_refactor = true;
            return None;
        }
        if q.is_none() {
            // Nothing in the window and the list: recompute the duals and rebuild, then scan once more.
            self.recomputed = true;
            self.y_ref = None;
            self.p0_rebuild();
            // The reference scans again only if the rebuild counted a dual infeasibility (a free variable counts always).
            let p0 = self.p0.as_ref().unwrap();
            let any = (0..self.p.n + self.p.m).any(|j| {
                self.p0_eligible(j)
                    && (p0.dstore[j] < -p0.tol || (is_inf(self.p.lo[j]) && is_inf(self.p.up[j])))
            });
            if !any {
                return None;
            }
            before(self);
            q = self.p0_scan();
            dump(self, "out", q);
            self.p0.as_mut().unwrap().calls += 1;
        }
        let (j, d) = q?;
        self.d[j] = d;
        Some(j)
    }
}

/// Drop tolerance of the pivot-row product (`0x134d90`'s first argument).
const ROW_DROP_TOL: f64 = 1e-13;
/// Pivot-row entries at or below this are left out of the Devex weight update (`0x103d0d`).
const WEIGHT_DROP_TOL: f64 = 1e-12;

/// Pricing calls of the reference's `0x105a50` (the dump's `n`): the retry after a failed scan of the candidate list is one too.
static PRC: std::sync::atomic::AtomicUsize = std::sync::atomic::AtomicUsize::new(0);

/// Penalty on the sum of bound violations in the composite phase-1 cost (recovered from the reference's `d`).
const PHASE1_M: f64 = 1e6;

/// Rounding guard for Devex weights of variables outside the reference framework (the reference's smallest observed weight is 3e-4).
const W_MIN: f64 = 1e-12;
/// Devex weight below which the reference resets the framework (`0x102d20`: an updated nonbasic weight below it, or the leaving weight at or below it).
const NORM_RESET_TOL: f64 = 1e-10;

impl<'a> Solver<'a> {
    fn sig(&self) -> Vec<i8> {
        let (p, tol) = (self.p, self.opts.feasibility_tol);
        self.basic
            .iter()
            .map(|&j| {
                let (li, ui) = (is_inf(p.lo[j]), is_inf(p.up[j]));
                if li && ui {
                    return 0;
                }
                let d = self.dst[j];
                if li {
                    return if d < -tol { 1 } else { 0 };
                }
                if d < -tol {
                    -1
                } else if !ui && d - (p.up[j] - p.lo[j]) > tol {
                    1
                } else {
                    0
                }
            })
            .collect()
    }

    /// The reference's dual clean-up pass after a shifted primal optimum (`0xf6a00`): the unshifted basic solution may leave a bound by
    /// more than the tolerance (FIT2P: one variable, 4.15e-4), and the reference then pivots with the dual simplex instead of re-entering
    /// the primal phase 1 (2 pivots, counted in its iteration count; rose's phase 1 took 6 to reach the same basis). Textbook dual
    /// simplex on the true costs: the leaving row is the largest bound violation, the entering variable the one with the least
    /// `|d_j / alpha_rj|` (ties: the largest `|alpha_rj|`); every pivot is followed by the clean-up refactorization (ascending order),
    /// which is the state the final solution is read from. Returns the number of pivots, or `None` when a row has no entering variable
    /// (the caller then falls back to the primal phase 1).
    fn dual_cleanup(&mut self) -> Result<Option<usize>, SolveError> {
        let p = self.p;
        let tol = self.opts.feasibility_tol;
        let mut pivots = 0;
        while pivots < 10 * p.m {
            let sg = self.sig();
            // Leaving row: the largest violation.
            let mut best: Option<(usize, f64, bool)> = None;
            for (r, &c) in sg.iter().enumerate() {
                if c == 0 {
                    continue;
                }
                let j = self.basic[r];
                let (li, d) = (is_inf(p.lo[j]), self.dst[j]);
                // `viol`: how far outside; `upper`: the bound it goes to.
                let (viol, upper) = if li {
                    (-d, true)
                } else if c < 0 {
                    (-d, false)
                } else {
                    (d - (p.up[j] - p.lo[j]), true)
                };
                if best.is_none_or(|b| viol > b.1) {
                    best = Some((r, viol, upper));
                }
            }
            let Some((r, _, to_upper)) = best else {
                return Ok(Some(pivots));
            };
            let _ = tol;
            // x_Br must increase when it goes to its lower bound, decrease when it goes to its upper one.
            let sr = if to_upper { -1.0 } else { 1.0 };
            let mut y: Vec<f64> = self.basic.iter().map(|&j| p.cost[j]).collect();
            self.factor.btran(&mut y);
            let mut rho = vec![0.0; p.m];
            rho[r] = 1.0;
            self.factor.btran(&mut rho);
            let mut enter: Option<(usize, f64, f64)> = None;
            for j in 0..p.n + p.m {
                if self.pos_of[j] != NONE || self.is_fixed(j) {
                    continue;
                }
                let alpha = dot(&p.cols[j], &rho);
                if alpha.abs() <= 1e-9 {
                    continue;
                }
                let (li, ui) = (is_inf(p.lo[j]), is_inf(p.up[j]));
                let ok = if li && ui {
                    true
                } else if self.at_upper[j] && !ui || li {
                    alpha * sr > 0.0
                } else {
                    alpha * sr < 0.0
                };
                if !ok {
                    continue;
                }
                let ratio = (p.cost[j] - dot(&p.cols[j], &y)).abs() / alpha.abs();
                let better = match enter {
                    None => true,
                    Some((_, br, ba)) => ratio < br || (ratio == br && alpha.abs() > ba),
                };
                if better {
                    enter = Some((j, ratio, alpha.abs()));
                }
            }
            let Some((q, _, _)) = enter else {
                return Ok(None);
            };
            let leave = self.basic[r];
            self.basic[r] = q;
            self.pos_of[q] = r;
            self.pos_of[leave] = NONE;
            self.at_upper[leave] = to_upper;
            self.refactor()?;
            pivots += 1;
        }
        Ok(None)
    }

    fn has_infeasible_basic(&self) -> bool {
        self.sig().iter().any(|&c| c != 0)
    }

    /// Reduced costs of the reference's infeasible-start mode: true cost plus `PHASE1_M` times the cost vector
    /// `-1` (basic below its lower bound) / `+1` (above its upper bound). `pure` drops the true cost and the
    /// penalty factor, which gives the pure sum-of-infeasibilities duals used for the infeasibility certificate.
    /// Returns the basic cost vector's multipliers `y`.
    fn composite_dual(&mut self, pure: bool) -> Vec<f64> {
        let p = self.p;
        if !pure {
            if let Some(y) = self.y_ref.take() {
                self.set_d_from_y(&y);
                return y;
            }
        }
        let sg = self.sig();
        let mut y: Vec<f64> = self
            .basic
            .iter()
            .enumerate()
            .map(|(i, &j)| {
                let inf = f64::from(sg[i]);
                if pure {
                    inf
                } else {
                    p.cost[j] + PHASE1_M * inf
                }
            })
            .collect();
        self.factor.btran(&mut y);
        for j in 0..p.n + p.m {
            self.d[j] = if self.pos_of[j] == NONE {
                (if pure { 0.0 } else { p.cost[j] }) - dot(&p.cols[j], &y)
            } else {
                0.0
            };
        }
        y
    }
}

/// One breakpoint-based ratio test result: entry index, step, `|pivot|`, and whether the leaving variable ends at its upper bound.
struct Phase1Pick {
    entry: Option<usize>,
    step: f64,
    pivot: f64,
    to_upper: bool,
}

/// Type byte of a basic variable for the phase-1 ratio test, in a frame where the variable has a finite lower bound
/// (variables bounded only above are mirrored): returns `(dist, range, type, mirrored)`; `None` for free variables.
fn classify(x: f64, lo: f64, up: f64, tol: f64, kept: f64) -> Option<(f64, f64, u8, bool)> {
    let (lo_inf, up_inf) = (is_inf(lo), is_inf(up));
    if lo_inf && up_inf {
        return None;
    }
    let _ = x;
    let (dist, range, mirrored) = if lo_inf {
        (kept, INF, true)
    } else {
        (kept, if up_inf { INF } else { up - lo }, false)
    };
    let t = if dist < -tol {
        b'B'
    } else if range < INF && dist - range > tol {
        b'A'
    } else if range < INF && dist - range >= -tol {
        if range <= tol {
            b'X'
        } else {
            b'U'
        }
    } else if dist <= tol {
        b'L'
    } else {
        b'N'
    };
    Some((dist, range, t, mirrored))
}

/// `s` plus (`sign` = 1) or minus (`sign` = -1) the weights of `part`, in the reference's order (`0xf0a58`..`0xf0cc0`):
/// fewer than 8 terms sequentially; otherwise four two-lane accumulators over blocks of 8 (`s` in lane 0 of the first),
/// combined as `((a0 + a1) + (a2 + a3))`, then lane 0 + lane 1, then the remaining terms one by one.
fn lane_sum(s: f64, part: &[(f64, usize, bool)], sign: f64, weight: &dyn Fn(usize) -> f64) -> f64 {
    let n = part.len();
    let nv = n & !7;
    if nv == 0 {
        return part.iter().fold(s, |s, b| s + sign * weight(b.1));
    }
    let mut acc = [[s, 0.0f64], [0.0; 2], [0.0; 2], [0.0; 2]];
    for blk in part[..nv].chunks(8) {
        for a in 0..4 {
            for l in 0..2 {
                acc[a][l] += sign * weight(blk[2 * a + l].1);
            }
        }
    }
    let t0 = acc[0][0] + acc[1][0];
    let t1 = acc[0][1] + acc[1][1];
    let u0 = acc[2][0] + acc[3][0];
    let u1 = acc[2][1] + acc[3][1];
    let mut r = (t0 + u0) + (t1 + u1);
    for b in &part[nv..] {
        r += sign * weight(b.1);
    }
    r
}

/// Mid-point search of the breakpoint where the cumulative slope reaches zero: a Hoare partition around the first element
/// that adds the lower part / subtracts the upper part of the weights `|a|` (transcribed from the reference).
fn find_theta(s0: f64, bp: &mut [(f64, usize, bool)], weight: &dyn Fn(usize) -> f64) -> f64 {
    let n = bp.len();
    let (mut s, mut lo, mut hi) = (s0, 0usize, n - 1);
    if hi > 0 {
        loop {
            let piv = bp[lo].0;
            let mut i = lo as isize - 1;
            let mut j = hi as isize + 1;
            loop {
                j -= 1;
                while piv < bp[j as usize].0 {
                    j -= 1;
                }
                i += 1;
                while bp[i as usize].0 < piv {
                    i += 1;
                }
                if j <= i {
                    break;
                }
                bp.swap(i as usize, j as usize);
            }
            let j = j as usize;
            if s >= 0.0 {
                s = lane_sum(s, &bp[j + 1..=hi], -1.0, weight);
            } else if lo <= j {
                s = lane_sum(s, &bp[lo..=j], 1.0, weight);
            }
            if s < 0.0 {
                lo = j + 1;
            } else {
                hi = j;
            }
            if lo >= hi {
                break;
            }
        }
    }
    bp[lo].0
}

/// The reference's infeasible-start ratio test (`0xf0120`, see `reverse-engineer/replay/phase1_ratio.py`).
/// `ent[k] = (position, alpha)` with `alpha > 0` meaning the basic variable decreases; `info(position)` gives the
/// typed frame from `classify` (already mirrored, so `alpha` is flipped for mirrored variables here).
/// `passes` receives the work of each pass: `Some((eligible, breakpoints))` for a breakpoint pass, `None` for a selection pass.
fn phase1_ratio(
    tol_h: f64,
    tol_p0: f64,
    own: f64,
    ent: &[(usize, f64)],
    info: &dyn Fn(usize) -> Option<(f64, f64, u8, bool)>,
    passes: &mut Vec<Option<(usize, usize)>>,
    dense: bool,
) -> Phase1Pick {
    let inf = INF;
    let frames: Vec<Option<(f64, f64, u8, bool)>> = ent.iter().map(|&(pos, _)| info(pos)).collect();
    let val = |k: usize| -> f64 {
        let a = ent[k].1;
        match frames[k] {
            Some((_, _, _, true)) => -a,
            _ => a,
        }
    };
    let mut tol_p = tol_p0;
    let mut first = true;
    loop {
        let (mut slope_a, mut slope_b) = (0.0f64, 0.0f64);
        let mut bp: Vec<(f64, usize, bool)> = Vec::new();
        for k in 0..ent.len() {
            let Some((d, r, c, _)) = frames[k] else {
                continue;
            };
            let a = val(k);
            if a <= tol_p {
                if !(a < -tol_p) {
                    continue;
                }
                match c {
                    b'A' => slope_a -= a,
                    b'B' => {
                        if d / a < own {
                            bp.push((d / a, k, false));
                        }
                        if r < inf && (d - r) / a < own {
                            bp.push(((d - r) / a, k, true));
                        }
                        slope_a += a;
                    }
                    b'L' | b'N' => {
                        if r < inf && (d - r) / a < own {
                            bp.push(((d - r) / a, k, true));
                        }
                    }
                    b'U' | b'X' => slope_b -= a,
                    _ => {}
                }
            } else {
                match c {
                    b'A' => {
                        if (d - r) / a < own {
                            bp.push(((d - r) / a, k, true));
                        }
                        if d / a < own {
                            bp.push((d / a, k, false));
                        }
                        slope_a -= a;
                    }
                    b'B' => slope_a += a,
                    b'L' | b'X' => slope_b += a,
                    // The reference's dense-result variant records only the breakpoints below the entering variable's range.
                    b'N' | b'U' => {
                        if !dense || d / a < own {
                            bp.push((d / a, k, false));
                        }
                    }
                    _ => {}
                }
            }
        }
        passes.push(Some((
            (0..ent.len()).filter(|&k| ent[k].1.abs() > tol_p).count(),
            bp.len(),
        )));
        let s0 = slope_b + slope_a;
        if s0 < 0.0 && bp.is_empty() {
            return Phase1Pick {
                entry: None,
                step: own.min(inf),
                pivot: 0.0,
                to_upper: false,
            };
        }
        let mut theta = 0.0;
        if s0 < 0.0 {
            let mut work = bp.clone();
            theta = find_theta(s0, &mut work, &|k| val(k).abs());
            if own < inf && own <= theta {
                return Phase1Pick {
                    entry: None,
                    step: own,
                    pivot: 0.0,
                    to_upper: false,
                };
            }
        }
        let (mut best, mut best_abs, mut step, mut to_upper) = (NONE, 0.0f64, 0.0f64, false);
        passes.push(None);
        for k in 0..ent.len() {
            let Some((d, r, c, _)) = frames[k] else {
                continue;
            };
            let a = val(k);
            let aa = a.abs();
            if a <= tol_p {
                if !(a < -tol_p) {
                    continue;
                }
                match c {
                    b'L' | b'N' => {
                        if r < inf && ((d - r) - theta * a).abs() * 2.0 <= tol_h && best_abs < aa {
                            (step, best, best_abs, to_upper) = ((d - r) / a, k, aa, true);
                        }
                    }
                    b'B' => {
                        if (d - theta * a).abs() * 2.0 <= tol_h && best_abs < aa {
                            (step, best, best_abs, to_upper) = (d / a, k, aa, false);
                        } else if r < inf
                            && ((d - r) - theta * a).abs() * 2.0 <= tol_h
                            && best_abs < aa
                        {
                            (step, best, best_abs, to_upper) = ((d - r) / a, k, aa, true);
                        }
                    }
                    b'U' | b'X' => {
                        if ((d - r) - theta * a).abs() <= tol_h && best_abs < aa {
                            (step, best, best_abs, to_upper) = ((d - r) / a, k, aa, true);
                        }
                    }
                    _ => {}
                }
            } else {
                match c {
                    b'N' | b'U' => {
                        if (d - theta * a).abs() * 2.0 <= tol_h && best_abs < aa {
                            (step, best, best_abs, to_upper) = (d / a, k, aa, false);
                        }
                    }
                    b'A' => {
                        if tol_h < ((d - r) - theta * a).abs() * 2.0 || a <= best_abs {
                            if (d - theta * a).abs() * 2.0 <= tol_h && best_abs < a {
                                (step, best, best_abs, to_upper) = (d / a, k, aa, false);
                            }
                        } else {
                            (step, best, best_abs, to_upper) = ((d - r) / a, k, aa, true);
                        }
                    }
                    b'L' | b'X' => {
                        if (d - theta * a).abs() <= tol_h && best_abs < a {
                            (step, best, best_abs, to_upper) = (d / a, k, aa, false);
                        }
                    }
                    _ => {}
                }
            }
        }
        if best == NONE {
            step = inf;
            for &(rr, k, tu) in &bp {
                if rr < step {
                    (step, best, to_upper) = (rr, k, tu);
                }
            }
            if best == NONE {
                return Phase1Pick {
                    entry: None,
                    step,
                    pivot: 0.0,
                    to_upper: false,
                };
            }
            // The routine itself does not compare this step with the entering variable's range; a bound flip is needed
            // when the entering variable reaches its opposite bound first.
            if own < inf && own <= step {
                return Phase1Pick {
                    entry: None,
                    step: own,
                    pivot: 0.0,
                    to_upper: false,
                };
            }
        }
        // The same flip when the second pass found a blocking entry beyond the entering variable's own range.
        if own < inf && own <= step {
            return Phase1Pick {
                entry: None,
                step: own,
                pivot: 0.0,
                to_upper: false,
            };
        }
        let pa = val(best).abs();
        // The iteration step (`0x143b8c`..`0x143cf8`) picks the bound the leaving variable ends at from the step: a type B variable
        // (below its lower bound) ends at the lower bound when `|d / alpha - step| < 2 tolH (step + 1)`, a type A variable (above
        // its upper bound) at the upper one when `|(d - r) / alpha - step| < 2 tolH (step + 1)`, and the update moves by that
        // bound's step (PILOT.JA iteration 13: both breakpoints of a range of 1e-5 are within the tolerance, the search picked
        // the upper one and the reference moves to the lower).
        match frames[best] {
            Some((d, _, b'B', _)) if (d / val(best) - step).abs() < 2.0 * tol_h * (step + 1.0) => {
                step = d / val(best);
                to_upper = false;
            }
            Some((d, r, b'A', _))
                if (r < inf && ((d - r) / val(best) - step).abs() < 2.0 * tol_h * (step + 1.0)) =>
            {
                step = (d - r) / val(best);
                to_upper = true;
            }
            _ => {}
        }
        // A mirrored variable reaches the opposite bound of the original frame.
        let mirrored = matches!(frames[best], Some((_, _, _, true)));
        let pick = Phase1Pick {
            entry: Some(best),
            step,
            pivot: pa,
            to_upper: to_upper != mirrored,
        };
        if !first || pa >= 0.1 {
            return pick;
        }
        let mx = ent
            .iter()
            .enumerate()
            .fold(0.0f64, |m, (k, _)| m.max(val(k).abs()));
        let mult = if mx >= 1e13 {
            1e-6
        } else if mx >= 1e10 {
            1e-5
        } else if mx >= 1e7 {
            1e-4
        } else if mx >= 1e4 {
            1e-3
        } else if mx >= 10.0 {
            1e-2
        } else {
            0.1
        };
        if tol_p * mx * mult <= pa {
            return pick;
        }
        tol_p *= mx * mult;
        first = false;
    }
}

/// The phase 1 ratio test used once the method has been switched (`ctx+0x468 == 0`, `0xeefa0`; `reverse-engineer/replay/ratio_eefa0.py`,
/// spec "Phase 1 ratio test with ctx+0x468 == 0"): a Harris test over the typed entries. Arguments as `phase1_ratio`.
fn phase1_eefa0(
    tol_h: f64,
    tol_p0: f64,
    own: f64,
    ent: &[(usize, f64)],
    info: &dyn Fn(usize) -> Option<(f64, f64, u8, bool)>,
    passes: &mut Vec<(usize, usize)>,
) -> Phase1Pick {
    let inf = INF;
    let frames: Vec<Option<(f64, f64, u8, bool)>> = ent.iter().map(|&(pos, _)| info(pos)).collect();
    let val = |k: usize| -> f64 {
        let a = ent[k].1;
        match frames[k] {
            Some((_, _, _, true)) => -a,
            _ => a,
        }
    };
    let none = |step: f64| Phase1Pick {
        entry: None,
        step,
        pivot: 0.0,
        to_upper: false,
    };
    let mut tol_p = tol_p0;
    let mut first = true;
    loop {
        let mut bound = own;
        let mut rec: Vec<(usize, f64, bool)> = Vec::new();
        for k in 0..ent.len() {
            let Some((d, r, c, _)) = frames[k] else {
                continue;
            };
            let a = val(k);
            if a <= tol_p {
                if !(a < -tol_p) {
                    continue;
                }
                let (dd, up) = match c {
                    b'L' | b'N' => {
                        if !(r < inf) {
                            continue;
                        }
                        (d - r, true)
                    }
                    b'B' => (d, false),
                    b'U' | b'X' => (d - r, true),
                    _ => continue,
                };
                let ratio = dd / a;
                if ratio <= bound {
                    rec.push((k, ratio, up));
                    bound = bound.min((dd - tol_h) / a);
                }
            } else {
                let (dd, up) = match c {
                    b'N' | b'U' | b'L' | b'X' => (d, false),
                    b'A' => (d - r, true),
                    _ => continue,
                };
                let ratio = dd / a;
                if ratio <= bound {
                    rec.push((k, ratio, up));
                    bound = bound.min((dd + tol_h) / a);
                }
            }
        }
        // Work of the pass: recorded and eligible entries (`0xef226` for a list, `0xef8c3` for a dense result; the capture names the next pc).
        passes.push((
            (0..ent.len()).filter(|&k| ent[k].1.abs() > tol_p).count(),
            rec.len(),
        ));
        if own < inf && own <= bound {
            return none(own);
        }
        if bound == inf || rec.is_empty() {
            return none(if first { inf } else { 0.0 });
        }
        let (mut best, mut best_abs, mut step, mut up) =
            (rec[rec.len() - 1].0, 0.0f64, 0.0f64, false);
        for &(k, ratio, u) in &rec {
            if ratio <= bound {
                let aa = val(k).abs();
                if best_abs < aa {
                    (step, best, best_abs, up) = (ratio, k, aa, u);
                }
            }
        }
        let mirrored = matches!(frames[best], Some((_, _, _, true)));
        let pick = Phase1Pick {
            entry: Some(best),
            step,
            pivot: best_abs,
            to_upper: up != mirrored,
        };
        if !first || best_abs >= 0.1 {
            return pick;
        }
        let mx = (0..ent.len()).fold(0.0f64, |m, k| m.max(val(k).abs()));
        let mult = if mx >= 1e13 {
            1e-6
        } else if mx >= 1e10 {
            1e-5
        } else if mx >= 1e7 {
            1e-4
        } else if mx >= 1e4 {
            1e-3
        } else if mx >= 10.0 {
            1e-2
        } else {
            0.1
        };
        let nt = tol_p * mx * mult;
        if nt <= best_abs {
            return pick;
        }
        tol_p = nt;
        first = false;
    }
}

/// Bucket thresholds of the candidate list: `T[0] = 0`, `T[k] = 2^(k-20)`.
fn thr(k: usize) -> f64 {
    if k == 0 {
        0.0
    } else {
        2f64.powi(k as i32 - 20)
    }
}

/// Largest `k` with `T[k] < score`, capped at the reference's top level 58.
fn bucket(score: f64) -> usize {
    let mut k = 0;
    while k < 58 && thr(k + 1) < score {
        k += 1;
    }
    k
}

/// Harris two-pass ratio test over `(position, a)` entries where `a > 0` means the basic variable decreases.
/// Returns `(entry index, step, |pivot|)`; `None` entry means no basic variable blocks (bound flip or unbounded).
#[allow(clippy::too_many_arguments)]
fn ratio_test(
    tol_h: f64,
    mut tol_p: f64,
    own: f64,
    entries: &[(usize, f64)],
    dist_down: &dyn Fn(usize) -> f64,
    dist_up: &dyn Fn(usize) -> f64,
    col_count: &dyn Fn(usize) -> Option<usize>,
    is_free: &dyn Fn(usize) -> bool,
    passes: &mut Vec<(usize, usize)>,
) -> (Option<usize>, f64, f64) {
    let mut first = true;
    loop {
        let mut theta = own;
        let mut cands: Vec<(usize, f64)> = Vec::new();
        // Entries past the pivot tolerance with a bound in their direction (the reference counts those moving down
        // unless the variable is free), for the pass's work.
        let mut eligible = 0usize;
        for (k, &(pos, a)) in entries.iter().enumerate() {
            if a > tol_p {
                let d = dist_down(pos);
                if !is_free(pos) {
                    eligible += 1;
                }
                if d < INF {
                    let r = d / a;
                    if r <= theta {
                        cands.push((k, r));
                        theta = theta.min((d + tol_h) / a);
                    }
                }
            } else if a < -tol_p {
                let d = dist_up(pos);
                if d < INF {
                    eligible += 1;
                    let r = d / -a;
                    if r <= theta {
                        cands.push((k, r));
                        theta = theta.min((d + tol_h) / -a);
                    }
                }
            }
        }
        passes.push((eligible, cands.len()));
        if own < INF && own <= theta {
            return (None, own, 0.0);
        }
        if theta >= INF || cands.is_empty() {
            return (None, INF, 0.0);
        }
        let (mut best_w, mut chosen, mut step) = (0.0, NONE, 0.0);
        for &(k, r) in &cands {
            if r <= theta {
                // Structural basic variables are favoured by their column count (the reference's tie-break among equal pivots).
                let mut w = entries[k].1.abs();
                if let Some(c) = col_count(entries[k].0) {
                    w *= (c as f64 - 1.0).min(10.0) / 10.0 + 1.0;
                }
                let better = best_w < w;
                if better {
                    best_w = w;
                    chosen = k;
                    step = r;
                }
            }
        }
        if chosen == NONE {
            return (None, INF, 0.0);
        }
        let pa = entries[chosen].1.abs();
        if !first || pa >= 0.1 {
            return (Some(chosen), step, pa);
        }
        let mx = entries.iter().fold(0.0f64, |m, e| m.max(e.1.abs()));
        let mult = if mx >= 1e13 {
            1e-6
        } else if mx >= 1e10 {
            1e-5
        } else if mx >= 1e7 {
            1e-4
        } else if mx >= 1e4 {
            1e-3
        } else if mx >= 10.0 {
            1e-2
        } else {
            0.1
        };
        if tol_p * mx * mult <= pa {
            return (Some(chosen), step, pa);
        }
        tol_p *= mx * mult;
        first = false;
    }
}

/// The reference's repair of a singular supplied basis ("variables dropped from basis"): the columns that get no pivot in the
/// reference's LU order leave the basis (nonbasic at their lower bound) and the slacks of the rows left without a pivot replace them.
/// Returns the last variable dropped (highest basis position), which the reference keeps out of pricing until the basis is
/// primal feasible (status -4 until `0x155139` adds 3 back).
pub(crate) fn repair_basis(p: &Internal, b: &mut Basis) -> Option<(usize, Vec<(u32, f64)>)> {
    let (n, m) = (p.n, p.m);
    let cols: Vec<(bool, Vec<usize>, Vec<f64>)> = b
        .basic
        .iter()
        .map(|&v| {
            (
                v >= n,
                p.fcols[v].iter().map(|&(r, _)| r).collect(),
                p.fcols[v].iter().map(|&(_, a)| a).collect(),
            )
        })
        .collect();
    crate::lu_order::reset_pair_gate();
    crate::lu_order::set_pool_dims(p.n, p.cols[..p.n].iter().map(Vec::len).sum());
    crate::lu_order::set_factorization_number(1);
    let st = crate::lu_order::factor_structure(m, &cols);
    let piv = st.pivots.clone();
    if piv.len() == m {
        return None;
    }
    // The reference's first factorization (of the singular basis) counts its work like any other (`Structure::work` with the full row count).
    let fwork = st.work(&cols, false);
    let mut row_done = vec![false; m];
    let mut pos_done = vec![false; m];
    for &(r, pos) in &piv {
        row_done[r] = true;
        pos_done[pos] = true;
    }
    let free_rows = (0..m).filter(|&r| !row_done[r]);
    let dropped = (0..m).filter(|&pos| !pos_done[pos]);
    let mut last = None;
    for (pos, r) in dropped.zip(free_rows) {
        let out = b.basic[pos];
        b.at_upper[out] = false;
        b.basic[pos] = n + r;
        b.at_upper[n + r] = false;
        last = Some(out);
    }
    last.map(|l| (l, fwork))
}

/// Solve starting from `start` (a primal feasible basis). `None` uses the all-slack basis with columns at their lower bounds.
pub fn solve(
    model: &LpModel,
    start: Option<&Basis>,
    opts: &SolveOptions,
) -> Result<SolveResult, SolveError> {
    if opts.method == Method::Pdhg {
        let r = crate::pdhg::solve(model, &opts.pdhg);
        return Ok(SolveResult {
            status: r.status,
            objective: r.objective,
            basis: Basis::slack(model.num_cols(), model.num_rows()),
            x: r.x,
            dual: r.dual,
            reduced_cost: r.reduced_cost,
            iterations: r.iterations as usize,
            refactorizations: 0,
            trace: Vec::new(),
            farkas: None,
            ray: None,
        });
    }
    if opts.method == Method::Barrier {
        return crate::barrier::solve_method(model, opts);
    }
    if opts.scaling {
        return crate::scale::solve(model, start, opts);
    }
    if opts.method == Method::Dual {
        return crate::dual::solve(model, start, opts);
    }
    let p = build(model);
    let (n, m) = (p.n, p.m);
    let mut excluded = None;
    let mut singular_work: Vec<(u32, f64)> = Vec::new();
    let basis = match start {
        Some(b) => {
            let mut b = b.clone();
            b.validate(n, m).map_err(SolveError::InvalidBasis)?;
            // The reference gets a basis as variable statuses, so its incoming position order is the variable order.
            b.basic.sort_unstable();
            if let Some((e, w)) = repair_basis(&p, &mut b) {
                excluded = Some(e);
                singular_work = w;
            }
            b
        }
        // Cold start: the reference's crash basis (`crash_basis`).
        None => crash_basis(model),
    };
    basis.validate(n, m).map_err(SolveError::InvalidBasis)?;
    let mut pos_of = vec![NONE; n + m];
    for (pos, &j) in basis.basic.iter().enumerate() {
        pos_of[j] = pos;
    }
    let in_r0: Vec<bool> = (0..n + m).map(|j| pos_of[j] == NONE).collect();
    let cols: Vec<Column> = basis.basic.iter().map(|&j| p.cols[j].clone()).collect();
    let factor = Factor::new(m, &cols).map_err(|_| SolveError::SingularBasis)?;
    let mut s = Solver {
        p: &p,
        opts,
        excluded,
        restart: false,
        pre_restart: false,
        fresh_duals: false,
        restarted: false,
        basic: basis.basic.clone(),
        pos_of,
        at_upper: basis.at_upper.clone(),
        x: vec![0.0; n + m],
        d: vec![0.0; n + m],
        w: vec![1.0; n + m],
        in_r0,
        factor,
        refactorizations: 1,
        level: 0,
        list: Vec::new(),
        p1_prev: Vec::new(),
        leave_prev: None,
        rows: Vec::new(),
        rrows: Vec::new(),
        list_failed: false,
        y_ref: None,
        y_rank: Vec::new(),
        unsnapped: Default::default(),
        want_refactor: false,
        dst: vec![0.0; p.n + p.m],
        rhs_shift: Vec::new(),
        mirror: None,
        work: crate::work::Work::default(),
        last_pivot: None,
        p1_tol_on: false,
        recomputed: false,
        rescan: false,
        cleanup: false,
        last_upd: None,
        last_row: Vec::new(),
        pending_row: None,
        upd_work: Vec::new(),
        p0: None,
        sw: crate::method::Method {
            factorizations: 1,
            ..Default::default()
        },
    };
    s.rows = vec![Vec::new(); m];
    s.rrows = vec![Vec::new(); m];
    for j in 0..n {
        for &(i, _) in &s.p.cols[j] {
            s.rows[i].push(j);
            // The reference's row-wise copy leaves out fixed columns.
            if s.pos_of[j] == NONE && !s.is_fixed(j) && excluded != Some(j) {
                s.rrows[i].push(j);
            }
        }
    }
    // The starting factorization's work, which the reference counts after its first refactorization check.
    let mut lu_work = Vec::new();
    if opts.faithful_ftran && opts.faithful_positions {
        // The reference factorizes the starting basis before the first iteration and numbers the positions by pivot step.
        crate::lu_order::reset_pair_gate();
        crate::lu_order::set_pool_dims(p.n, p.cols[..p.n].iter().map(Vec::len).sum());
        let (order, st, rows) = s.reference_order(1);
        for (pos, &j) in order.iter().enumerate() {
            s.pos_of[j] = pos;
        }
        s.basic = order;
        s.mirror = Some(crate::ftran_order::RefStruct::fresh(m, n, &st, &rows));
        let mi = s.mirror.as_ref().unwrap();
        lu_work = st.work(&rows, mi.has_levels());
        lu_work.extend(mi.fresh_work());

        let cols: Vec<&Column> = s.basic.iter().map(|&j| &s.p.cols[j]).collect();
        s.factor = Factor::new(m, &cols).map_err(|_| SolveError::SingularBasis)?;
    }
    s.work.trig();
    s.work.decided(s.work.w);
    if !singular_work.is_empty() && s.mirror.is_some() {
        // The singular start basis: the driver factors it first (its cost enters the smoothed cost `ctx+0x158` and the factorization
        // count), then the repaired basis.
        s.work.factor_driver(&singular_work, crate::work::CTX);
        for (site, est) in singular_work {
            s.work.add(site, est, crate::work::CTX);
        }
        s.work.factor_done();
    }
    s.work.factor_driver(&lu_work, crate::work::CTX);
    for (site, est) in lu_work {
        s.work.add(site, est, crate::work::CTX);
    }
    if opts.pricing0 {
        assert!(
            s.mirror.is_some(),
            "pricing0 needs faithful_ftran and faithful_positions"
        );
        s.sw.active = true;
        s.sw.pert = PERTURB_VALUE;
        let nnz: usize = (0..n).map(|j| p.cols[j].len()).sum();
        s.p0 = Some(P0 {
            list: Vec::new(),
            pos: 0,
            dstore: vec![0.0; n + m],
            y: Vec::new(),
            pending: None,
            tol: P0_DTOL,
            rebuild: true,
            calls: 0,
            window: P0_WINDOW,
            bk: p0_bucket_mode(n, m, nnz).then(|| Buckets {
                head: [-1; 101],
                tail: [-1; 101],
                next: vec![-1; n + m],
                visited: 0,
                built: false,
                r8: p0_bucket_r8(n, m, nnz),
                full: true,
            }),
        });
    }
    s.work.factor_done();
    s.recompute_primal();
    s.primal_full_work();
    s.dual_btran_work();
    s.recompute_dual();
    let mut trace = Vec::new();
    let mut iterations = 0usize;
    // Pricing passes of the main loop, the call index of the `pricing` dump.
    if crate::dump::on("header") {
        let fixed: Vec<usize> = (0..p.n + p.m).filter(|&j| s.is_fixed(j)).collect();
        let upper: Vec<usize> = (0..p.n + p.m)
            .filter(|&j| s.pos_of[j] == NONE && s.at_upper[j])
            .collect();
        crate::dump::emit(
            "header",
            serde_json::json!({"n": p.n, "m": p.m, "fixed": fixed, "basic": s.basic, "upper": upper}),
        );
    }
    let mut phase1 = false;
    let mut refactor_next = false;
    let mut farkas: Option<Vec<f64>> = None;
    let mut ray: Option<Vec<f64>> = None;
    let status = loop {
        if iterations >= opts.max_iterations {
            break Status::IterationLimit;
        }
        let due = match &opts.refactor_schedule {
            Some(sched) => sched.contains(&iterations),
            None if opts.faithful_ftran && s.mirror.is_some() => {
                std::mem::replace(&mut refactor_next, false)
            }
            None => s.factor.num_updates() >= opts.max_updates,
        };
        // The reference's check reads `W` at the end of the iteration (right after `trig()`), before this update's own work.
        let mut pending_upd: Vec<(u32, f64, f64)> = Vec::new();
        let mut chk: Option<(bool, usize, usize)> = None;
        if due {
            s.refactor()?;
        }
        // The reference makes this check after the ratio test, and not at all for a bound flip: `check_state` takes it back.
        let mut check_state = None;
        if opts.refactor_schedule.is_none() {
            if let Some(mi) = s.mirror.as_ref() {
                let nnz_a: usize = p.cols[..p.n].iter().map(Vec::len).sum();
                check_state = Some(s.work.decision_state());
                chk = Some((mi.refactor_due(), mi.updates(), nnz_a));
            }
        }
        // The reference's primal update of the last pivot (`0x105f20`, site `0x107b57`): a pass over the FTRAN result, which
        // also collects the phase 1 cost changes of the basic variables that changed feasibility class.
        let last = s.last_pivot.take().filter(|_| !due);
        let mut dc: Vec<(usize, f64)> = Vec::new();
        // The entering variable alone changes class: the reference neither solves nor takes a row product for it (`0x107e14`).
        let mut q_only = false;
        if let Some(lp) = &last {
            if phase1 {
                let cls = s.sig();
                for &i in &lp.rows {
                    let (j, c) = s.p1_prev[i];
                    if j == s.basic[i] && cls[i] != c {
                        dc.push((i, PHASE1_M * f64::from(cls[i] - c)));
                    }
                }
                // The variable that just became basic may be infeasible (the infeasible clamp below lets it): its position
                // takes the penalty too.
                if lp.q != NONE && s.pos_of[lp.q] != NONE {
                    let r = s.pos_of[lp.q];
                    if cls[r] != 0 {
                        dc.push((r, PHASE1_M * f64::from(cls[r])));
                        q_only = dc.len() == 1;
                    }
                }
            }
            if s.mirror.is_some() && s.p0.is_none() {
                let est = if lp.dense {
                    m as f64 * 3.0
                } else {
                    lp.rows.len() as f64 * 6.0
                } + if q_only { 0.0 } else { dc.len() as f64 * 4.0 };
                s.work.add(0x107b57, est, crate::work::CTX);
            }
        }
        let infeasible = s.has_infeasible_basic();
        s.p1_tol_on = infeasible;
        s.fresh_duals = false;

        if infeasible {
            let was = phase1;
            phase1 = true;
            if s.p0.is_some() {
                // `pricing0`: no maintained `d`; the cost shifts join the BTRAN input of the dual update (a refactorization or a
                // fresh start rebuilds the duals instead).
                s.p0_cost_changes(was && !due, &dc);
            } else if was && !due {
                // The reference corrects `d` in place when basic variables change feasibility class and re-scores the
                // entries of the row of the cost changes.
                if s.mirror.is_some() {
                    // The reference keeps `d` and applies the cost changes to it in place.
                    if let (Some(lp), false) = (&last, dc.is_empty()) {
                        s.cost_change_work(&dc, lp.q, lp.leave, lp.alpha_rq, q_only);
                    }
                } else {
                    // Without the mirror: re-score the columns whose `d` moved.
                    let old_d = s.d.clone();
                    let cls = s.sig();
                    let cost_changed = s
                        .p1_prev
                        .iter()
                        .any(|&(j, c)| s.pos_of[j] != NONE && cls[s.pos_of[j]] != c);
                    s.composite_dual(false);
                    let tol = s.opts.optimality_tol;
                    let delta: Vec<usize> = (0..s.d.len())
                        .filter(|&j| (s.d[j] - old_d[j]).abs() > tol)
                        .collect();
                    if cost_changed && !delta.is_empty() {
                        s.update_list(&delta, 0, 0x108280);
                    }
                }
            } else {
                s.composite_dual(false);
                s.list.clear();
                s.level = 0;
            }
            s.p1_prev = s.basic.iter().zip(s.sig()).map(|(&j, c)| (j, c)).collect();
        } else if phase1 {
            // Back in a feasible region: continue with the true costs (weights and framework are kept).
            phase1 = false;
            if s.p0.is_some() && !due {
                // The end of phase 1 recomputes the duals from the true costs (one dense BTRAN) and keeps the list.
                s.dual_btran_work();
                let p0 = s.p0.as_mut().unwrap();
                if !p0.rebuild {
                    p0.y = s.y_ref.clone().unwrap();
                    p0.pending = None;
                }
            }
            // The reference's primal update of the last pivot still runs the phase 1 cost-change pass on the iteration that ends
            // phase 1 (its work counts toward the refactorization clause).
            if s.mirror.is_some() && !due && s.p0.is_none() {
                if let (Some(lp), false) = (&last, dc.is_empty()) {
                    s.cost_change_work(&dc, lp.q, lp.leave, lp.alpha_rq, q_only);
                }
            }
            s.list.clear();
            s.level = 0;
            // After a refactorization the duals were just recomputed with the same (now feasible) basic values.
            if !due && !s.p0.is_some() {
                s.dual_btran_work();
                s.recompute_dual();
                s.fresh_duals = true;
            }
        }
        if false {
            let fr: Vec<usize> = (0..p.n + p.m)
                .filter(|&j| s.pos_of[j] == NONE && is_inf(p.lo[j]) && is_inf(p.up[j]))
                .collect();
            eprintln!(
                "FREE it={iterations} nfree={} list={} faithful={} first={:?}",
                fr.len(),
                s.list.len(),
                opts.faithful_weights,
                &fr[..fr.len().min(5)]
            );
        }
        if crate::dump::on("pricing") {
            crate::dump::emit(
                "pricing",
                serde_json::json!({"n": PRC.load(std::sync::atomic::Ordering::Relaxed), "iter": iterations, "infeasible": infeasible, "list": s.list, "basic": s.basic, "d": s.d, "w": s.w, "x": s.x}),
            );
        }
        PRC.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
        if let Some(p0) = s.p0.as_mut() {
            p0.tol = if infeasible { P0_DTOL_P1 } else { P0_DTOL };
        }
        let q = if infeasible && !s.p0.is_some() {
            match s.price_admitting() {
                Some(q) => q,
                None => {
                    // No improving column for the penalised cost: try the pure sum of infeasibilities, which either
                    // continues the search or proves infeasibility with the multipliers as certificate.
                    let y = s.composite_dual(true);
                    s.list.clear();
                    s.level = 0;
                    match s.price_admitting() {
                        Some(q) => q,
                        None => {
                            farkas = Some(y.iter().map(|v| -v).collect());
                            break Status::Infeasible;
                        }
                    }
                }
            }
        } else {
            match if s.p0.is_some() {
                s.price0()
            } else {
                s.price_admitting()
            } {
                Some(q) => q,
                None if std::mem::take(&mut s.want_refactor) => {
                    s.work.trig();
                    s.work.decided(s.work.w);
                    if let Some(mi) = s.mirror.as_mut() {
                        mi.count_all_updates();
                    }
                    refactor_next = true;
                    continue;
                }
                None => {
                    // The reference recomputes the duals once more when the scan finds nothing (`0x145520`, with the residual check).
                    if s.mirror.is_some() && !s.recomputed && !s.fresh_duals {
                        s.dual_btran_work();
                    }
                    // Before declaring optimality, rebuild the basic solution from a fresh factorization: drift in the updated
                    // values can hide a bound violation just above the tolerance.
                    if s.factor.num_updates() > 0 {
                        s.cleanup = s.mirror.is_some() && !s.rhs_shift.is_empty();
                        s.work.trig();
                        s.refactor()?;
                        if s.cleanup && s.has_infeasible_basic() {
                            // The dual clean-up pass; the primal phase 1 only when it finds no entering variable.
                            if let Some(k) = s.dual_cleanup()? {
                                iterations += k;
                            }
                        }
                        if s.has_infeasible_basic() {
                            continue;
                        }
                    }
                    if infeasible && s.has_infeasible_basic() {
                        // The scan found no improving column in phase 1 on a basis that is still infeasible: the pure sum of infeasibilities
                        // either continues the search or proves infeasibility with the multipliers as certificate (as without the scan).
                        let y = s.composite_dual(true);
                        s.list.clear();
                        s.level = 0;
                        match s.price_admitting() {
                            Some(q) => q,
                            None => {
                                farkas = Some(y.iter().map(|v| -v).collect());
                                break Status::Infeasible;
                            }
                        }
                    } else {
                        break Status::Optimal;
                    }
                }
            }
        };
        if std::mem::take(&mut s.restarted) && opts.refactor_schedule.is_none() {
            // The check of this iteration reads the factor that the restart replaced: take it again from the new one.
            if let Some(mi) = s.mirror.as_ref() {
                let nnz_a: usize = p.cols[..p.n].iter().map(Vec::len).sum();
                check_state = Some(s.work.decision_state());
                chk = Some((mi.refactor_due(), mi.updates(), nnz_a));
            }
        }
        let mut alpha = vec![0.0; m];
        let dir = if s.d[q] < 0.0 { 1.0 } else { -1.0 };
        let mut ftran_out: Vec<(usize, f64)> = Vec::new();
        let mut ft_dense = false;
        let entries: Vec<(usize, f64)> = if let Some(mi) = s.mirror.as_mut() {
            // The reference's FTRAN, value for value and in its output order.
            let rows_q: Vec<usize> = p.fcols[q].iter().map(|&(i, _)| i).collect();
            let vals_q: Vec<f64> = p.fcols[q].iter().map(|&(_, a)| a).collect();
            let out = mi.ftran(&rows_q, &vals_q);
            if crate::dump::on("ftran") {
                crate::dump::emit(
                    "ftran",
                    serde_json::json!({"iter": iterations, "q": q, "in": rows_q, "dense": mi.last_dense(), "out": crate::dump::sparse(&out)}),
                );
            }
            if crate::dump::on("lstage") {
                let (l, e) = mi.lstage_snap();
                crate::dump::emit(
                    "lstage",
                    serde_json::json!({"iter": iterations, "q": q, "l": l, "e": e}),
                );
            }
            if crate::dump::on("ustage") {
                let (spike, diag, uc, slot, kind, i21) = mi.ustage_state();
                let tail: Vec<serde_json::Value> = (0..diag.len())
                    .filter(|&p| kind[p] == -2)
                    .map(|p| serde_json::json!([p, slot[p], crate::dump::sparse(&uc[p])]))
                    .collect();
                crate::dump::emit(
                    "ustage",
                    serde_json::json!({"iter": iterations, "q": q, "i21": i21, "spike": crate::dump::sparse(spike), "diag": diag, "tail": tail}),
                );
            }
            s.work
                .add_from(0x1219d6, 0x120673, mi.work, crate::work::FAC);
            for &(i, a) in &out {
                alpha[i] = a;
            }
            ft_dense = mi.last_dense();
            let e = out.iter().map(|&(i, a)| (i, dir * a)).collect();
            ftran_out = out;
            e
        } else {
            for &(i, a) in &p.cols[q] {
                alpha[i] = a;
            }
            s.factor.ftran(&mut alpha);
            (0..m)
                .filter(|&i| alpha[i].abs() > 1e-14)
                .map(|i| (i, dir * alpha[i]))
                .collect()
        };
        let own = if !is_inf(p.lo[q]) && !is_inf(p.up[q]) {
            p.up[q] - p.lo[q]
        } else {
            INF
        };
        let mut leave_to_upper: Option<bool> = None;
        let (chosen, mut step, _pivot) = if infeasible {
            let (xs, basic, lo, up) = (&s.x, &s.basic, &p.lo, &p.up);
            let tol = opts.feasibility_tol;
            let info = |pos: usize| {
                let j = basic[pos];
                classify(xs[j], lo[j], up[j], tol, s.dst[j])
            };
            if crate::dump::on("ent") {
                let dims = s.mirror.as_ref().map(|mi| mi.dims());
                let kinds = s.mirror.as_ref().map(|mi| mi.kinds().to_vec());
                let bl = s
                    .mirror
                    .as_ref()
                    .map(|mi| mi.blists().to_vec())
                    .unwrap_or_default();
                // Per basis position, as the reference's ratio-test arrays: distance, range (null when infinite) and type
                // letter. A free variable has no class; the dump gives it type 'F' and its value as the distance.
                let cls: Vec<_> = (0..m).map(&info).collect();
                let dist: Vec<f64> = cls
                    .iter()
                    .enumerate()
                    .map(|(pos, c)| c.map_or(xs[basic[pos]], |c| c.0))
                    .collect();
                let rng: Vec<Option<f64>> = cls
                    .iter()
                    .map(|c| c.and_then(|c| (c.1 < INF).then_some(c.1)))
                    .collect();
                let typ: String = cls.iter().map(|c| c.map_or('F', |c| c.2 as char)).collect();
                crate::dump::emit(
                    "ent",
                    serde_json::json!({"iter": iterations, "q": q, "own": (own < INF).then_some(own),
                    "idx": entries.iter().map(|e| e.0).collect::<Vec<_>>(), "val": entries.iter().map(|e| e.1).collect::<Vec<_>>(),
                    "d": dist, "rng": rng, "typ": typ, "dense": ft_dense, "bl": bl, "kinds": kinds, "dims": dims}),
                );
            }
            let mut passes = Vec::new();
            let eefa = s.sw.mode != 0;
            let mut passes2: Vec<(usize, usize)> = Vec::new();
            let pick = if eefa {
                phase1_eefa0(tol, opts.pivot_tol, own, &entries, &info, &mut passes2)
            } else {
                phase1_ratio(
                    tol,
                    opts.pivot_tol,
                    own,
                    &entries,
                    &info,
                    &mut passes,
                    ft_dense,
                )
            };
            if let Some(mi) = s.mirror.as_ref() {
                let dense = mi.last_dense();
                for (eligible, rec) in passes2 {
                    let (site, est) = if dense {
                        (0xef8d5, (2 * rec + 2 * eligible + m) as f64)
                    } else {
                        (0xef238, (2 * rec + 4 * entries.len()) as f64)
                    };
                    s.work.add(site, est, crate::work::CTX);
                }
            }
            if let Some(mi) = s.mirror.as_ref().filter(|_| !eefa) {
                // The reference's phase 1 ratio test (`0xf0120`): a breakpoint pass and a selection pass over the FTRAN result.
                let dense = mi.last_dense();
                for pass in passes {
                    let (site, est) = match (pass, dense) {
                        (Some((eligible, stored)), true) => {
                            (0xf15aa, (4 * stored + eligible + m) as f64)
                        }
                        (Some((_, stored)), false) => {
                            (0xf0427, (4 * stored + 2 * entries.len()) as f64)
                        }
                        (None, true) => (0xf125c, 3.0 * m as f64),
                        (None, false) => (0xf07e3, 4.0 * entries.len() as f64),
                    };
                    s.work.add(site, est, crate::work::CTX);
                }
            }
            let mut step_used = pick.step;
            if let Some(k) = pick.entry {
                leave_to_upper = Some(pick.to_upper);
                // To the upper bound the iteration step moves by the distance of the leaving variable's position, rebuilt in the original
                // frame, to the bound: `((dist + lo) - up) / alpha` (`(up - dist)` for a variable without lower bound), not the ratio
                // test's `(dist - range) / alpha`, which differs in the last bits when `lo != 0` (PILOT.JA iterations 1096, 1123,
                // 1219: the reference's next distances equal `dist - step alpha_i` with this step for every entry).
                let leave = s.basic[entries[k].0];
                if let Some(c) = info(entries[k].0).filter(|_| pick.to_upper) {
                    let bound = p.up[leave];
                    let xr = if is_inf(p.lo[leave]) {
                        p.up[leave] - c.0
                    } else {
                        c.0 + p.lo[leave]
                    };
                    if !is_inf(bound) && entries[k].1 != 0.0 {
                        step_used = (xr - bound) / entries[k].1;
                    }
                }
            }
            (pick.entry, step_used, pick.pivot)
        } else {
            let (xs, basic, lo, up) = (&s.x, &s.basic, &p.lo, &p.up);
            // After a perturbation the reference's Harris test (`0xef960`) reads the shifted distance array (`dst`), not `x_B`.
            let perturbed = s.sw.level > 0;
            let dst = &s.dst;
            let down = |pos: usize| {
                let j = basic[pos];
                if is_inf(lo[j]) {
                    INF
                } else {
                    dst[j]
                }
            };
            let upd = |pos: usize| {
                let j = basic[pos];
                if is_inf(up[j]) {
                    INF
                } else if !is_inf(lo[j]) {
                    (up[j] - lo[j]) - dst[j]
                } else if perturbed {
                    dst[j]
                } else {
                    up[j] - xs[j]
                }
            };
            if crate::dump::on("ent2") {
                let dist: Vec<f64> = (0..m)
                    .map(|pos| {
                        let j = basic[pos];
                        if perturbed || !is_inf(lo[j]) {
                            dst[j]
                        } else {
                            xs[j]
                        }
                    })
                    .collect();
                let rng: Vec<Option<f64>> = (0..m)
                    .map(|pos| {
                        let j = basic[pos];
                        (!is_inf(lo[j]) && !is_inf(up[j])).then_some(up[j] - lo[j])
                    })
                    .collect();
                let typ: String = (0..m)
                    .map(|pos| {
                        let j = basic[pos];
                        if is_inf(lo[j]) && is_inf(up[j]) {
                            'F'
                        } else {
                            'N'
                        }
                    })
                    .collect();
                crate::dump::emit(
                    "ent2",
                    serde_json::json!({"iter": iterations, "q": q, "own": (own < INF).then_some(own),
                    "idx": entries.iter().map(|e| e.0).collect::<Vec<_>>(), "val": entries.iter().map(|e| e.1).collect::<Vec<_>>(),
                    "d": dist, "rng": rng, "typ": typ}),
                );
            }
            let mut passes = Vec::new();
            // Set up at `0x138722`: the weight applies from 500 rows (the field is 0 below that).
            let col_count = |pos: usize| {
                let j = basic[pos];
                (m >= 500 && j < n).then(|| p.cols[j].len())
            };
            let free = |pos: usize| {
                let j = basic[pos];
                is_inf(lo[j]) && is_inf(up[j])
            };
            let res = ratio_test(
                opts.feasibility_tol,
                opts.pivot_tol,
                own,
                &entries,
                &down,
                &upd,
                &col_count,
                &free,
                &mut passes,
            );
            if let Some(mi) = s.mirror.as_ref() {
                // The reference's ratio test (`0xef960`): per pass, a scan of the FTRAN result and its stored candidates.
                let dense = mi.last_dense();
                for (eligible, stored) in passes {
                    let (site, est) = if dense {
                        (0xf010c, (2 * eligible + m + 3 * stored) as f64)
                    } else {
                        (0xefb7d, (4 * entries.len() + 3 * stored) as f64)
                    };
                    s.work.add(site, est, crate::work::CTX);
                }
            }
            let mut res = res;
            if let Some(k) = res.0 {
                // To the upper bound the step is the distance of the position rebuilt in the original frame, `((dist + lo) - up) / alpha`,
                // as in phase 1 (PILOT.JA phase 2 ratio call 17).
                let (pos, a) = entries[k];
                let j = basic[pos];
                if a < 0.0 && !is_inf(lo[j]) && !is_inf(up[j]) {
                    res.1 = ((dst[j] + lo[j]) - up[j]) / a;
                }
            }
            res
        };
        if chosen.is_none() && step >= INF && infeasible && own < INF {
            // No admissible pivot and no improving slope in phase 1: the entering variable goes to its opposite bound.
            step = own;
        }
        if chosen.is_none() && step >= INF {
            if infeasible {
                // Rebuild the factorization once and price again.
                if s.factor.num_updates() > 0 {
                    s.refactor()?;
                    continue;
                }
                break Status::NumericalFailure;
            }
            let mut dir_x = vec![0.0; n + m];
            dir_x[q] = dir;
            for &(i, a) in &entries {
                dir_x[s.basic[i]] -= a;
            }
            dir_x.truncate(model.num_cols());
            ray = Some(dir_x);
            break Status::Unbounded;
        }
        // A phase 1 step may be negative (a basic variable inside its tolerance leaves from the wrong side of its bound).
        iterations += 1;
        if s.sw.active {
            let counted = chosen.is_some_and(|k| !s.is_fixed(s.basic[entries[k].0]))
                && !(is_inf(p.lo[q]) && is_inf(p.up[q]));
            s.sw.iteration(chosen.is_none(), counted, step);
        }
        // After a perturbation the leaving variable's position is judged in the shifted frame of the distance array.
        let shifted_frame = s.sw.level > 0 && !s.rhs_shift.is_empty();
        let x_leave_before = chosen.map(|k| {
            let j = s.basic[entries[k].0];
            if shifted_frame && !infeasible && !(is_inf(p.lo[j]) && is_inf(p.up[j])) {
                if is_inf(p.lo[j]) {
                    p.up[j] - s.dst[j]
                } else {
                    p.lo[j] + s.dst[j]
                }
            } else {
                s.x[j]
            }
        });
        // Move the basic variables and the entering variable.
        // Phase 2 does not move the point by a negative Harris step (reference, PEROLD ratio calls 50, 82, 354, ...: the
        // basics stay as they are, the entering variable stays at its bound and the leaving one keeps its computed value).
        let move_step = if infeasible { step } else { step.max(0.0) };
        for &(i, a) in &entries {
            let j = s.basic[i];
            s.x[j] -= move_step * a;
            s.dst[j] = if is_inf(p.lo[j]) {
                s.dst[j] + move_step * a
            } else {
                s.dst[j] - move_step * a
            };
        }
        let dq_base = s.dist_of(q, s.nonbasic_value(q));
        if s.unsnapped.remove(&q) {
            // A variable that left outside its bound enters from the bound, not from the value it kept.
            s.x[q] = s.nonbasic_value(q);
        }
        s.x[q] += dir * move_step;
        s.dst[q] = if is_inf(p.lo[q]) {
            dq_base - dir * move_step
        } else if dir < 0.0 {
            // An entering variable that leaves its upper bound has its distance rebuilt from the moved value (PILOT.JA
            // iteration 1187: (3.5 - step) - 2.0, not 1.5 - step).
            s.dist_of(q, s.x[q])
        } else {
            dq_base + dir * move_step
        };
        let Some(k) = chosen else {
            // Bound flip: the entering variable travels to its opposite bound.
            for (site, est, scale) in std::mem::take(&mut pending_upd) {
                s.work.add(site, est, scale);
            }
            if let Some(st) = check_state {
                refactor_next = false;
                s.work.restore_decision(st);
            }
            s.at_upper[q] = !s.at_upper[q];
            s.x[q] = s.nonbasic_value(q);
            s.y_ref = None;
            // The reference still runs the primal update of the move at the next iteration (`0x105f20`): the phase 1 cost
            // changes of the basic variables whose class changed re-score the columns they touch.
            s.last_pivot = Some(LastPivot {
                rows: ftran_out.iter().map(|e| e.0).collect(),
                dense: s.mirror.as_ref().is_some_and(|mi| mi.last_dense()),
                q: NONE,
                leave: NONE,
                alpha_rq: 0.0,
            });
            if opts.record_trace {
                trace.push(TraceEvent {
                    iteration: iterations,
                    entering: q,
                    leaving: None,
                    position: None,
                    step,
                    flips: Vec::new(),
                });
            }
            continue;
        };
        let (r, a_r) = entries[k];
        let leave = s.basic[r];
        let to_upper = leave_to_upper.unwrap_or(a_r < 0.0);
        // Snap the leaving variable to the bound it reached. In phase 2 a variable that ends outside its bound (a negative
        // Harris step) keeps its computed value, which the reference's next recomputation of `x_B` uses.
        let x_before = [x_leave_before.unwrap_or(0.0)];
        let bound = if to_upper { p.up[leave] } else { p.lo[leave] };
        let outside = if to_upper {
            x_before[0] > bound
        } else {
            x_before[0] < bound
        };
        if infeasible || !outside {
            s.x[leave] = bound;
            s.unsnapped.remove(&leave);
        } else {
            s.x[leave] = x_before[0];
            s.unsnapped.insert(leave);
            // `0x143f38`: the shift of a leaving variable that ends outside its bound also marks the shifted frame (`ctx+0x338` 0 -> 1).
            if s.sw.active && s.sw.level == 0 {
                s.sw.level = 1;
            }
            let delta = bound - x_before[0];
            if s.rhs_shift.is_empty() {
                s.rhs_shift = vec![0.0; p.m];
            }
            for &(i, a) in &p.cols[leave] {
                s.rhs_shift[i] += a * delta;
            }
        }
        let alpha_rq = alpha[r];

        if s.p0.is_some() {
            // `pricing0`: no pivot row, no Devex weights, no maintained `d`; the next pricing call updates the duals by one BTRAN.
            let dq = s.d[q];
            s.y_ref = None;
            s.pos_of[leave] = NONE;
            s.at_upper[leave] = to_upper;
            if s.p0_eligible(leave) {
                // The bucket-chain list takes the leaving variable at the tail of its lowest bucket (key 1), where the next scan reaches it last.
                if let Some(bk) = s.p0.as_mut().unwrap().bk.as_mut() {
                    bk.append(1, leave);
                }
            }
            s.basic[r] = q;
            s.pos_of[q] = r;
            s.factor.update(r, alpha);
            if let Some(mi) = s.mirror.as_mut() {
                mi.update(r);
                s.upd_work = std::mem::take(&mut mi.upd_work);
            }
            s.leave_prev = Some(leave);
            s.last_pivot = Some(LastPivot {
                rows: {
                    let mut rows: Vec<usize> = ftran_out.iter().map(|e| e.0).collect();
                    if let Some(k) = rows.iter().position(|&i| i == r) {
                        rows.swap_remove(k);
                    }
                    rows
                },
                dense: s.mirror.as_ref().is_some_and(|mi| mi.last_dense()),
                q,
                leave,
                alpha_rq,
            });
            s.p0.as_mut().unwrap().pending = Some(vec![(r, dq)]);
            s.work.trig();
            if let Some((d, u, na)) = chk {
                let w_now = s.work.w;
                refactor_next = d || s.work.clause(u, na, w_now);
                if refactor_next {
                    s.work.decided(w_now);
                }
            }
            let upd_work = std::mem::take(&mut s.upd_work);
            if !refactor_next {
                for (site, est, scale) in upd_work {
                    s.work.add(site, est, scale);
                }
            }
            if opts.record_trace {
                trace.push(TraceEvent {
                    iteration: iterations,
                    entering: q,
                    leaving: Some(leave),
                    position: Some(r),
                    step,
                    flips: Vec::new(),
                });
            }
            continue;
        }

        // Pivot row alpha_r = e_r' B^-1 [A I] and the Devex cross-term vector tau = B^-T u.
        let mut rho_ref: Option<Vec<(usize, f64)>> = None;
        // The pivot-row update's own estimate (`0x102d20`, site `0x1040c7`) apart from its pivot-row term.
        let mut devex_est = 0.0f64;
        // Method 2 (`0x104a03`) has no cross-term BTRAN and a simpler work estimate.
        let m2 = s.sw.mode == 2;
        let mut alpha_dense = false;
        let mut alpha_ordered = false;
        let mut tau_ref: Option<Vec<(usize, f64)>> = None;
        if let Some(mi) = s.mirror.as_mut() {
            rho_ref = Some(mi.btran(&[(r, 1.0)], true, false));
            if crate::dump::on("rho") {
                crate::dump::emit(
                    "rho",
                    serde_json::json!({"iter": iterations, "r": r, "out": crate::dump::sparse(rho_ref.as_deref().unwrap())}),
                );
            }
            s.work
                .add_from(0x13174d, 0x102e5b, mi.work, crate::work::FAC);
            // The reference's input list for the cross-term BTRAN (`0x102d20`): the pivot entry first (when its variable is in `R0`),
            // then the pivot column's list with the pivot entry dropped by moving the last entry into its slot (when the column
            // came out as a list rather than ascending).
            let inp: Vec<(usize, f64)> = match ftran_out.iter().position(|e| e.0 == r) {
                Some(k) if !alpha_ordered => {
                    let mut rest = ftran_out.clone();
                    let piv = rest.swap_remove(k);
                    let mut v: Vec<(usize, f64)> = Vec::with_capacity(rest.len() + 1);
                    if s.in_r0[s.basic[piv.0]] {
                        v.push(piv);
                    }
                    v.extend(rest.into_iter().filter(|&(i, _)| s.in_r0[s.basic[i]]));
                    v
                }
                _ => ftran_out
                    .iter()
                    .copied()
                    .filter(|&(i, _)| s.in_r0[s.basic[i]])
                    .collect(),
            };
            // The cross-term vector goes in dense when the FTRAN result was dense or too long for a list (`0x102d20`).
            // By then the list no longer holds the pivot entry.
            let dense_in = mi.last_dense() || m <= mi.dthr() * (ftran_out.len() - 1);
            alpha_dense = mi.last_dense();
            alpha_ordered = mi.last_ordered();
            let alen = ftran_out.len() - 1;
            if dense_in && !m2 {
                devex_est = m as f64;
            }
            devex_est = if mi.last_dense() {
                m as f64 * 1.5
            } else {
                alen as f64 * 3.0
            } + devex_est;
            let mut tau_len = Some(0usize);
            tau_ref = Some(Vec::new());
            if !m2 && (dense_in || !inp.is_empty()) {
                let tau = mi.btran(&inp, false, dense_in);
                if crate::dump::on("tau") {
                    crate::dump::emit(
                        "tau",
                        serde_json::json!({"iter": iterations, "out": crate::dump::sparse(&tau)}),
                    );
                }
                tau_ref = Some(tau.clone());
                s.work
                    .add_from(0x12de50, 0x103ac7, mi.work, crate::work::FAC);
                tau_len = if mi.bt_dense { None } else { Some(tau.len()) };
            }
            if let (Some(t), false) = (tau_len, m2) {
                devex_est = devex_est + m as f64 + t as f64 * 2.0;
            }
        }
        let mut rho = vec![0.0; m];
        if let Some(rr) = &rho_ref {
            // The reference's BTRAN values and support (`0x131710`).
            for &(i, v) in rr {
                rho[i] = v;
            }
        } else {
            rho[r] = 1.0;
            s.factor.btran(&mut rho);
        }
        let mut u = vec![0.0; m];
        if let Some(t) = &tau_ref {
            for &(i, v) in t {
                u[i] = v;
            }
        } else {
            for i in 0..m {
                if s.in_r0[s.basic[i]] {
                    u[i] = alpha[i];
                }
            }
            s.factor.btran(&mut u);
        }
        if s.mirror.is_some() {
            let rho_work: Vec<f64> = rho
                .iter()
                .map(|&v| if v.abs() > ROW_DROP_TOL { v } else { 0.0 })
                .collect();
            s.row_product_work(&rho_work, q, leave, 0x105498);
        }
        let theta_d = s.d[q] / alpha_rq;
        // The reference recomputes the entering column's weight from the pivot column (`0x102d20`): the entering variable's own 1,
        // the pivot square, then the squares of the other entries in list order, over the reference framework.
        let wq = if ftran_out.is_empty() {
            s.w[q]
        } else {
            // Mode 2 (`0x104a03`) takes the entering variable's 1 always and counts the variables whose weight is positive instead of the flags.
            let m2 = s.sw.mode == 2;
            let inr = |j: usize| if m2 { s.w[j] > 0.0 } else { s.in_r0[j] };
            let mut t = if m2 || s.in_r0[q] { 1.0 } else { 0.0 };
            if inr(s.basic[r]) {
                t += alpha_rq * alpha_rq;
            }
            // The reference drops the pivot entry from the list by moving the last entry into its slot, then sums.
            let mut order = ftran_out.clone();
            if let (false, Some(k)) = (alpha_ordered, order.iter().position(|e| e.0 == r)) {
                order.swap_remove(k);
            }
            for &(i, v) in &order {
                if i != r && inr(s.basic[i]) {
                    t += v * v;
                }
            }
            t
        };
        let col_off = &p.offs;
        let rho_rank: Option<Vec<usize>> = rho_ref.as_ref().map(|rr| {
            let mut rank = vec![usize::MAX; m];
            for (k, &(i, _)) in rr.iter().enumerate() {
                // Entries below the drop tolerance are not part of the reference's row-product input.
                if rho[i].abs() > ROW_DROP_TOL {
                    rank[i] = k;
                }
            }
            rank
        });
        // The leaving variable takes part in the row update like any nonbasic: its row product is 1 up to rounding.
        let prods: Option<Vec<f64>> = rho_rank.as_ref().map(|rank| p.rank_products(&rho, rank));
        let leave_arj = match prods.as_ref() {
            Some(pr) => pr[leave],
            None => 1.0,
        };
        let theta_leave = match rho_rank.as_ref() {
            Some(_) => theta_d * leave_arj,
            None => theta_d,
        };
        // The entering variable takes part in the row update before it turns basic: a slack with its row product, a
        // structural with the pivot element (`alpha_rq`), so `d` ends as rounding noise (exactly 0.0 in most pivots).
        let arq_x = match prods.as_ref() {
            Some(pr) => pr[q],
            None => dot(&p.cols[q], &rho),
        };
        let mut norm_reset = false;
        let mut changed: Vec<usize> = Vec::new();
        let mut rowvals: Vec<(usize, f64)> = Vec::new();
        for j in 0..n + m {
            if j >= n && j != q && j != leave && s.pos_of[j] != NONE && s.mirror.is_some() {
                let arj = match prods.as_ref() {
                    Some(pr) => pr[j],
                    None => dot(&p.cols[j], &rho),
                };
                if arj.abs() > ROW_DROP_TOL {
                    if !s.is_fixed(j) {
                        s.d[j] -= theta_d * arj;
                    }
                    // The reference keeps a basic slack's rounding-noise row product in the pivot row like any entry.
                    changed.push(j);
                    rowvals.push((j, arj));
                }
                continue;
            }
            if s.pos_of[j] != NONE || j == q {
                continue;
            }
            let arj = match prods.as_ref() {
                Some(pr) => pr[j],
                None => dot(&p.cols[j], &rho),
            };
            // The reference's row product keeps the entries above its drop tolerance (`0x134d90`).
            if arj.abs() <= ROW_DROP_TOL {
                continue;
            }
            // The Devex update has its own, larger drop tolerance (`0x103d0d`); d and the list keep the row product's.
            if m2 {
                // Mode 2 (`0x104c2d`): `w_j = max(w_j, wq / alpha^2 * arj^2 + 1)`, no cross term.
                if arj.abs() > WEIGHT_DROP_TOL {
                    s.w[j] = s.w[j].max(((wq / alpha_rq) / alpha_rq) * arj * arj + 1.0);
                }
            } else if arj.abs() > WEIGHT_DROP_TOL {
                let ratio = arj / alpha_rq;
                let vj = dot_lanes(&p.fcols[j], &u, col_off[j] & 1);
                s.w[j] = ratio * ((wq * ratio - vj) - vj) + s.w[j];
                // A fixed structural never resets the framework (PILOT.JA iteration 3856: the fixed EMROP06 went to -3e4); a fixed
                // slack does (PEROLD iteration 1626, PILOTNOV iteration 1026).
                if opts.faithful_weights && s.w[j] < NORM_RESET_TOL && !(j < n && s.is_fixed(j)) {
                    norm_reset = true;
                }
                let floor = if opts.faithful_weights { W_MIN } else { 1.0 };
                s.w[j] = s.w[j].max(floor);
            }
            if !s.is_fixed(j) {
                // The reference sums the row product in the order of the BTRAN output list.
                s.d[j] -= theta_d * arj;
            }
            // The variable dropped from a singular start basis is not in the reference's row copy, so not in its pivot row.
            if s.excluded != Some(j) {
                changed.push(j);
                rowvals.push((j, arj));
            }
        }
        // The penalty leaves `d` of the leaving variable after the list update, which scores it with the penalty still in.
        let mut leave_penalty = 0.0;
        if s.mirror.is_some() && phase1 {
            // The reference keeps the stored value of a basic variable: the leaving one moves by the pivot-row update like any
            // other entry and loses its phase 1 penalty; the entering one is left holding the negated pivot element (zero for a
            // degenerate step; this matches 277 of 330 DEGEN2 pivots).
            s.d[leave] -= theta_leave;
            if let Some(&(j, c)) = s.p1_prev.get(r) {
                if j == leave && c != 0 {
                    leave_penalty = PHASE1_M * f64::from(c);
                }
            }
            s.d[q] = s.d[q] - theta_d * if q < n { alpha_rq } else { arq_x };
        } else {
            // The reference subtracts from the stored value (rounding noise included) in phase 2 as well.
            if s.mirror.is_some() {
                s.d[leave] -= theta_leave;
            } else {
                s.d[leave] = -theta_leave;
            }
            s.d[q] = if s.mirror.is_some() {
                s.d[q] - theta_d * if q < n { alpha_rq } else { arq_x }
            } else {
                0.0
            };
        }
        s.w[leave] = (wq / alpha_rq) / alpha_rq;
        if m2 {
            // `0x104b4d`: the error of the stored weight of the entering variable against the recomputed one, larger over smaller.
            let ws = s.w[q];
            s.sw.ratio_sum += if wq > ws { wq / ws } else { ws / wq };
            s.sw.ratio_cnt += 1;
        }
        if m2 {
            // Mode 2: the leaving weight is floored at 1 and the pivot-row loop then reaches the leaving variable with its row product.
            let wl = s.w[leave].max(1.0);
            s.w[leave] = wl;
            if leave_arj.abs() > WEIGHT_DROP_TOL {
                s.w[leave] = wl.max(((wq / alpha_rq) / alpha_rq) * leave_arj * leave_arj + 1.0);
            }
            s.w[q] = 1.0;
        } else {
            if opts.faithful_weights && s.w[leave] <= NORM_RESET_TOL {
                norm_reset = true;
            }
            s.w[leave] = s.w[leave].max(if opts.faithful_weights { W_MIN } else { 1.0 });
        }

        // Basis change.
        s.y_ref = None;
        s.pos_of[leave] = NONE;
        s.at_upper[leave] = to_upper;
        s.basic[r] = q;
        s.pos_of[q] = r;
        if q < n {
            for &(i, _) in &p.cols[q] {
                if let Some(k) = s.rrows[i].iter().position(|&j| j == q) {
                    s.rrows[i].swap_remove(k);
                }
            }
        }
        if leave < n && !s.is_fixed(leave) {
            for &(i, _) in &p.cols[leave] {
                s.rrows[i].push(leave);
            }
        }
        s.factor.update(r, alpha);
        if let Some(mi) = s.mirror.as_mut() {
            mi.update(r);
            s.upd_work = std::mem::take(&mut mi.upd_work);
        }
        changed.push(leave);
        s.leave_prev = Some(leave);
        s.last_pivot = Some(LastPivot {
            rows: {
                // The reference removes the pivot entry from the list by moving the last entry into its place.
                let mut rows: Vec<usize> = ftran_out.iter().map(|e| e.0).collect();
                if let Some(k) = rows.iter().position(|&i| i == r) {
                    rows.swap_remove(k);
                }
                rows
            },
            dense: alpha_dense,
            q,
            leave,
            alpha_rq,
        });
        if crate::dump::on("ord") {
            let rho_nz = rho_ref.as_ref().map(|r| {
                r.iter()
                    .map(|e| (e.0, e.1 != 0.0 && rho[e.0] != 0.0))
                    .collect::<Vec<_>>()
            });
            crate::dump::emit(
                "ord",
                serde_json::json!({"iter": iterations, "rho": rho_nz, "changed": changed}),
            );
        }
        let changed = match (true, rho_ref.as_ref()) {
            (true, Some(rr)) => {
                let nz: Vec<usize> = rr
                    .iter()
                    .map(|e| e.0)
                    .filter(|&i| rho[i].abs() > ROW_DROP_TOL)
                    .collect();
                {
                    // The reference's density test counts the entries that survive the drop tolerance.
                    let cnt = nz.iter().filter(|&&i| rho[i].abs() > ROW_DROP_TOL).count();
                    s.pivot_row_order_by(&nz, &changed, cnt)
                }
            }
            _ => s.pivot_row_order(&rho, &changed),
        };
        if s.mirror.is_some() {
            let nn = p.n;
            let nnz_a: usize = p.cols[..nn].iter().map(Vec::len).sum();
            let rowlen = changed
                .iter()
                .filter(|&&j| j >= nn || !s.is_fixed(j))
                .count()
                + 1;
            let est = if m2 {
                3.0 * rowlen as f64 + devex_est
            } else {
                2.0 * rowlen as f64
                    * ((m as f64 + m as f64 + nnz_a as f64) / (m + 1 + nn) as f64 + 1.5)
                    + devex_est
            };
            s.work.add_from(0x1040c7, 0x143ff9, est, crate::work::CTX);
        }
        // The leaving variable's entry of the pivot row is `rho_r . a_leave` (1 up to rounding).
        let arl = match rho_rank.as_ref() {
            Some(rank) => rank_dot(&p.cols[leave], &rho, rank),
            None => dot(&p.cols[leave], &rho),
        };
        rowvals.push((leave, arl));
        rowvals.push((q, if q < n { alpha_rq } else { arq_x }));
        s.pending_row = Some(rowvals);
        if s.mirror.is_some() && opts.faithful_weights {
            s.p1_tol_on = s.has_infeasible_basic();
        }
        // The reference's pivot row also carries the entering column (`0x102d20`).
        // Mode 2 restarts the weights when the average error exceeds 3 (`0x10882c`, before the list update).
        let m2_reset = m2 && s.sw.ratio_sum / f64::from(s.sw.ratio_cnt) > 3.0;
        if m2_reset {
            s.sw.ratio_sum = 0.0;
            s.sw.ratio_cnt = 1;
        }
        if (norm_reset && opts.faithful_weights) || m2_reset {
            // The reset sets the weights; the candidate list keeps its entries and level and takes the normal update.
            s.reset_weights();
            s.update_list(&changed, 1, 0x1442ab);
        } else {
            s.update_list(&changed, 1, 0x1442ab);
        }
        s.d[leave] -= leave_penalty;
        s.work.trig();
        if let Some((d, u, na)) = chk {
            let w_now = s.work.w;
            refactor_next = d || s.work.clause(u, na, w_now);
            if refactor_next {
                s.work.decided(w_now);
            }
        }
        // The reference updates the factor after the check; when the check calls for a refactorization the update is not made.
        let upd_work = std::mem::take(&mut s.upd_work);
        if !refactor_next {
            for (site, est, scale) in upd_work {
                s.work.add(site, est, scale);
            }
        }
        if opts.record_trace {
            trace.push(TraceEvent {
                iteration: iterations,
                entering: q,
                leaving: Some(leave),
                position: Some(r),
                step,
                flips: Vec::new(),
            });
        }
    };

    if status == Status::Optimal {
        // The solution is read from the factorization that was just built; its work was booked then.
        s.work.mute = s.cleanup;
        s.recompute_primal();
        s.work.mute = false;
        s.recompute_dual();
    }
    let sgn = if model.sense == Sense::Maximize {
        -1.0
    } else {
        1.0
    };
    let mut y: Vec<f64> = s.basic.iter().map(|&j| p.cost[j]).collect();
    s.factor.btran(&mut y);
    if status == Status::Optimal {
        if let Some(yr) = s.y_ref.as_ref() {
            y = yr.clone();
        }
    }
    let nm = model.num_cols();
    let x: Vec<f64> = s.x[..nm].to_vec();
    // The reference's primal (controlled profile) returns a certificate that uses an inequality row in its own convention (`verify::farkas_margins`:
    // `y >= 0` on rows with only a lower side, `<= 0` on rows with only an upper side); with equality rows only, as in the release profile and
    // the dual method, the opposite signs.
    let farkas = farkas.map(|y| {
        let one_sided = (0..model.num_rows())
            .any(|i| y[i] != 0.0 && (is_inf(model.row_lo[i]) || is_inf(model.row_up[i])));
        let neg: Vec<f64> = y.iter().map(|v| -v).collect();
        let (a, a_neg) = (
            crate::verify::farkas_margins(model, &y).0,
            crate::verify::farkas_margins(model, &neg).0,
        );
        if !opts.pricing0 && one_sided && a_neg > 0.0 && a_neg > a {
            neg
        } else {
            y
        }
    });
    Ok(SolveResult {
        status,
        objective: model.objective(&x),
        x,
        dual: y.iter().map(|v| sgn * v).collect(),
        reduced_cost: s.d[..nm].iter().map(|v| sgn * v).collect(),
        basis: Basis {
            basic: s.basic.clone(),
            at_upper: s.at_upper.clone(),
        },
        iterations,
        refactorizations: s.refactorizations,
        trace,
        farkas,
        ray,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::mps;

    #[test]
    fn solves_a_small_lp_with_slack_start() {
        // max 3x + 2y  s.t. x + y <= 4, x + 3y <= 6, x <= 3, x,y >= 0   -> x=3, y=1, z=11
        let text = "NAME T\nOBJSENSE\n MAX\nROWS\n N obj\n L c1\n L c2\nCOLUMNS\n x obj 3 c1 1\n x c2 1\n y obj 2 c1 1\n y c2 3\nRHS\n rhs c1 4 c2 6\nBOUNDS\n UP bnd x 3\nENDATA\n";
        let model = mps::parse(text).unwrap();
        let r = solve(&model, None, &SolveOptions::default()).unwrap();
        assert_eq!(r.status, Status::Optimal);
        assert!(
            (r.objective - 11.0).abs() < 1e-9,
            "objective {}",
            r.objective
        );
        assert!((r.x[0] - 3.0).abs() < 1e-9 && (r.x[1] - 1.0).abs() < 1e-9);
    }

    #[test]
    fn detects_unbounded() {
        let text = "NAME T\nROWS\n N obj\n G c1\nCOLUMNS\n x obj -1 c1 1\nRHS\n rhs c1 1\nENDATA\n";
        let model = mps::parse(text).unwrap();
        // x >= 1 is infeasible at the lower bound 0, so start from the basis where the slack is basic at x=1... use x basic.
        let start = Basis {
            basic: vec![0],
            at_upper: vec![false, false],
        };
        let r = solve(&model, Some(&start), &SolveOptions::default()).unwrap();
        assert_eq!(r.status, Status::Unbounded);
    }
}
