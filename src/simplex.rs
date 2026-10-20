//! Primal simplex from a supplied feasible basis.
//!
//! Rules follow the recovered reference behaviour (see `reverse-engineer/specs/primal-simplex.md`):
//! Devex pricing over a fixed reference framework, Harris two-pass ratio test with bound flips, and
//! reduced-cost updates from the pivot row. Internal form: `A x + s = rhs`, columns `0..n` structural then
//! one slack per row (`+e_i` for `=`/`<=`, `-e_i` for `>=`).

use crate::basis::Basis;
use crate::lu::{Column, Factor};
use crate::model::{LpModel, Sense, INF};

#[derive(Clone, Debug)]
pub struct SolveOptions {
    pub feasibility_tol: f64,
    pub optimality_tol: f64,
    pub pivot_tol: f64,
    pub max_iterations: usize,
    /// Refactor after this many product-form updates.
    pub max_updates: usize,
    pub record_trace: bool,
    /// Reproduce the reference's Devex weights exactly (weights outside the framework may fall below 1). The reference's
    /// anti-degeneracy perturbation is not implemented, so this mode can stall on degenerate models; off by default.
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
}

impl Default for SolveOptions {
    fn default() -> Self {
        SolveOptions {
            feasibility_tol: 1e-6,
            optimality_tol: 1e-6,
            pivot_tol: 1e-9,
            max_iterations: 1_000_000,
            max_updates: 20,
            record_trace: true,
            faithful_weights: false,
            faithful_positions: false,
            refactor_schedule: None,
            faithful_ftran: false,
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
}

impl std::fmt::Display for SolveError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            SolveError::InvalidBasis(s) => write!(f, "invalid basis: {s}"),
            SolveError::SingularBasis => write!(f, "starting basis is singular"),
        }
    }
}

impl std::error::Error for SolveError {}

const NONE: usize = usize::MAX;

struct Internal {
    n: usize,
    m: usize,
    cols: Vec<Column>,
    lo: Vec<f64>,
    up: Vec<f64>,
    rhs: Vec<f64>,
    cost: Vec<f64>,
}

fn is_inf(v: f64) -> bool {
    v.abs() >= INF
}

/// Internal layout: structurals, then one range column per ranged inequality row (as the reference does), then
/// one slack per row. A ranged row `rl <= a x <= ru` becomes `a x + s + g = ru` with `s` fixed at 0 and
/// `g in [0, ru - rl]`; a ranged `>=` row uses `-1` for both unit columns.
fn build(model: &LpModel) -> Internal {
    let (nm, m) = (model.num_cols(), model.num_rows());
    let mut cols: Vec<Column> = (0..nm).map(|j| model.col(j).collect()).collect();
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
                up.push(ru - rl);
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
    let sgn = if model.sense == Sense::Maximize { -1.0 } else { 1.0 };
    let mut cost: Vec<f64> = model.obj.iter().map(|c| sgn * c).collect();
    cost.resize(n + m, 0.0);
    Internal { n, m, cols, lo, up, rhs, cost }
}

struct Solver<'a> {
    p: &'a Internal,
    opts: &'a SolveOptions,
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
    /// Structural columns of each row, ascending.
    rows: Vec<Vec<usize>>,
    /// Mirror of the reference's FTRAN structure (`faithful_ftran`).
    mirror: Option<crate::ftran_order::RefStruct>,
}

fn dot(col: &Column, v: &[f64]) -> f64 {
    col.iter().map(|&(r, a)| a * v[r]).sum()
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
    fn reference_order(&self) -> (Vec<usize>, crate::lu_order::Structure, Vec<Vec<usize>>) {
        let n = self.p.n;
        let cols: Vec<(bool, Vec<usize>, Vec<f64>)> = self
            .basic
            .iter()
            .map(|&v| (v >= n, self.p.cols[v].iter().map(|&(r, _)| r).collect(), self.p.cols[v].iter().map(|&(_, a)| a).collect()))
            .collect();
        let st = crate::lu_order::factor_structure(self.p.m, &cols);
        let order = st.pivots.iter().map(|&(_, pos)| self.basic[pos]).collect();
        let rows = cols.into_iter().map(|c| c.1).collect();
        (order, st, rows)
    }

    fn refactor(&mut self) -> Result<(), SolveError> {
        if self.opts.faithful_positions {
            let (order, st, rows) = self.reference_order();
            debug_assert_eq!(order.len(), self.basic.len());
            for (pos, &j) in order.iter().enumerate() {
                self.pos_of[j] = pos;
            }
            self.basic = order;
            if self.opts.faithful_ftran {
                self.mirror = Some(crate::ftran_order::RefStruct::fresh(self.p.m, &st, &rows));
            }
        }
        let cols: Vec<Column> = self.basic.iter().map(|&j| self.p.cols[j].clone()).collect();
        self.factor = Factor::new(self.p.m, &cols).map_err(|_| SolveError::SingularBasis)?;
        self.refactorizations += 1;
        self.list.clear();
        self.recompute_primal();
        self.recompute_dual();
        Ok(())
    }

    fn recompute_primal(&mut self) {
        let p = self.p;
        let mut r = p.rhs.clone();
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
        self.factor.ftran(&mut r);
        for (pos, &j) in self.basic.iter().enumerate() {
            self.x[j] = r[pos];
        }
    }

    fn recompute_dual(&mut self) {
        let p = self.p;
        let mut y: Vec<f64> = self.basic.iter().map(|&j| p.cost[j]).collect();
        self.factor.btran(&mut y);
        for j in 0..p.n + p.m {
            self.d[j] = if self.pos_of[j] == NONE { p.cost[j] - dot(&p.cols[j], &y) } else { 0.0 };
        }
    }

    fn is_fixed(&self, j: usize) -> bool {
        self.p.lo[j] == self.p.up[j]
    }

    /// Order in which the reference stores the pivot-row entries: with a sparse `rho`, structurals in first-touch
    /// order over the rows of `rho`; otherwise structurals in index order; then slacks in row order.
    fn pivot_row_order(&self, rho: &[f64], changed: &[usize]) -> Vec<usize> {
        let (n, m) = (self.p.n, self.p.m);
        let nz: Vec<usize> = (0..m).rev().filter(|&i| rho[i] != 0.0).collect();
        let nnz_a: usize = self.p.cols[..n].iter().map(Vec::len).sum();
        let mut keep = vec![false; n + m];
        for &j in changed {
            keep[j] = true;
        }
        let mut out = Vec::with_capacity(changed.len());
        if nz.len() * 3 * nnz_a < m * n {
            for &i in &nz {
                for &j in &self.rows[i] {
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
        for &i in &nz {
            if keep[n + i] {
                keep[n + i] = false;
                out.push(n + i);
            }
        }
        out.extend((0..n + m).filter(|&j| keep[j]));
        out
    }

    /// Eligibility and score `d^2 / w` of a nonbasic variable, or `None` if it cannot enter.
    fn score(&self, j: usize) -> Option<f64> {
        if self.pos_of[j] != NONE || self.is_fixed(j) {
            return None;
        }
        let tol = self.opts.optimality_tol;
        let free = is_inf(self.p.lo[j]) && is_inf(self.p.up[j]);
        let dj = self.d[j];
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
        if let Some(j) = self.best_in_list() {
            return Some(j);
        }
        self.level = 0;
        self.list.clear();
        for j in 0..self.p.n + self.p.m {
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
        self.best_in_list()
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
    fn update_list(&mut self, changed: &[usize]) {
        let mut fresh: Vec<usize> = Vec::new();
        for &j in changed {
            if let Some(sc) = self.score(j) {
                if thr(self.level) < sc {
                    let b = bucket(sc);
                    if self.level < b {
                        self.level = b;
                        fresh.clear();
                    }
                    fresh.push(j);
                }
            }
        }
        let old = std::mem::take(&mut self.list);
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
    }
}

/// Penalty on the sum of bound violations in the composite phase-1 cost (recovered from the reference's `d`).
const PHASE1_M: f64 = 1e6;

/// Rounding guard for Devex weights of variables outside the reference framework (the reference's smallest observed weight is 3e-4).
const W_MIN: f64 = 1e-6;

impl<'a> Solver<'a> {
    fn has_infeasible_basic(&self) -> bool {
        let tol = self.opts.feasibility_tol;
        self.basic.iter().any(|&j| {
            (!is_inf(self.p.lo[j]) && self.x[j] < self.p.lo[j] - tol) || (!is_inf(self.p.up[j]) && self.x[j] > self.p.up[j] + tol)
        })
    }

    /// Reduced costs of the reference's infeasible-start mode: true cost plus `PHASE1_M` times the cost vector
    /// `-1` (basic below its lower bound) / `+1` (above its upper bound). `pure` drops the true cost and the
    /// penalty factor, which gives the pure sum-of-infeasibilities duals used for the infeasibility certificate.
    /// Returns the basic cost vector's multipliers `y`.
    fn composite_dual(&mut self, pure: bool) -> Vec<f64> {
        let p = self.p;
        let tol = self.opts.feasibility_tol;
        let mut y: Vec<f64> = self
            .basic
            .iter()
            .map(|&j| {
                let inf = if !is_inf(p.lo[j]) && self.x[j] < p.lo[j] - tol {
                    -1.0
                } else if !is_inf(p.up[j]) && self.x[j] > p.up[j] + tol {
                    1.0
                } else {
                    0.0
                };
                if pure { inf } else { p.cost[j] + PHASE1_M * inf }
            })
            .collect();
        self.factor.btran(&mut y);
        for j in 0..p.n + p.m {
            self.d[j] = if self.pos_of[j] == NONE { (if pure { 0.0 } else { p.cost[j] }) - dot(&p.cols[j], &y) } else { 0.0 };
        }
        self.list.clear();
        self.level = 0;
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
fn classify(x: f64, lo: f64, up: f64, tol: f64) -> Option<(f64, f64, u8, bool)> {
    let (lo_inf, up_inf) = (is_inf(lo), is_inf(up));
    if lo_inf && up_inf {
        return None;
    }
    let (dist, range, mirrored) = if lo_inf { (-x + up, INF, true) } else { (x - lo, if up_inf { INF } else { up - lo }, false) };
    let t = if dist < -tol {
        b'B'
    } else if range < INF && dist - range > tol {
        b'A'
    } else if range < INF && dist - range >= -tol {
        if range <= tol { b'X' } else { b'U' }
    } else if dist <= tol {
        b'L'
    } else {
        b'N'
    };
    Some((dist, range, t, mirrored))
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
                for t in j + 1..=hi {
                    s -= weight(bp[t].1);
                }
            } else if lo <= j {
                for t in lo..=j {
                    s += weight(bp[t].1);
                }
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
fn phase1_ratio(tol_h: f64, tol_p0: f64, own: f64, ent: &[(usize, f64)], info: &dyn Fn(usize) -> Option<(f64, f64, u8, bool)>) -> Phase1Pick {
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
            let Some((d, r, c, _)) = frames[k] else { continue };
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
                    b'N' | b'U' => bp.push((d / a, k, false)),
                    _ => {}
                }
            }
        }
        let s0 = slope_b + slope_a;
        if s0 < 0.0 && bp.is_empty() {
            return Phase1Pick { entry: None, step: own.min(inf), pivot: 0.0, to_upper: false };
        }
        let mut theta = 0.0;
        if s0 < 0.0 {
            let mut work = bp.clone();
            theta = find_theta(s0, &mut work, &|k| val(k).abs());
            if own < inf && own <= theta {
                return Phase1Pick { entry: None, step: own, pivot: 0.0, to_upper: false };
            }
        }
        let (mut best, mut best_abs, mut step, mut to_upper) = (NONE, 0.0f64, 0.0f64, false);
        for k in 0..ent.len() {
            let Some((d, r, c, _)) = frames[k] else { continue };
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
                        } else if r < inf && ((d - r) - theta * a).abs() * 2.0 <= tol_h && best_abs < aa {
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
                return Phase1Pick { entry: None, step, pivot: 0.0, to_upper: false };
            }
            // The routine itself does not compare this step with the entering variable's range; a bound flip is needed
            // when the entering variable reaches its opposite bound first.
            if own < inf && own <= step {
                return Phase1Pick { entry: None, step: own, pivot: 0.0, to_upper: false };
            }
        }
        let pa = val(best).abs();
        // A mirrored variable reaches the opposite bound of the original frame.
        let mirrored = matches!(frames[best], Some((_, _, _, true)));
        let pick = Phase1Pick { entry: Some(best), step, pivot: pa, to_upper: to_upper != mirrored };
        if !first || pa >= 0.1 {
            return pick;
        }
        let mx = ent.iter().enumerate().fold(0.0f64, |m, (k, _)| m.max(val(k).abs()));
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

/// Bucket thresholds of the candidate list: `T[0] = 0`, `T[k] = 2^(k-20)`.
fn thr(k: usize) -> f64 {
    if k == 0 {
        0.0
    } else {
        2f64.powi(k as i32 - 20)
    }
}

/// Largest `k` with `T[k] < score`.
fn bucket(score: f64) -> usize {
    let mut k = 0;
    while thr(k + 1) < score {
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
) -> (Option<usize>, f64, f64) {
    let mut first = true;
    loop {
        let mut theta = own;
        let mut cands: Vec<(usize, f64)> = Vec::new();
        for (k, &(pos, a)) in entries.iter().enumerate() {
            if a > tol_p {
                let d = dist_down(pos);
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
                    let r = d / -a;
                    if r <= theta {
                        cands.push((k, r));
                        theta = theta.min((d + tol_h) / -a);
                    }
                }
            }
        }
        if own < INF && own <= theta {
            return (None, own, 0.0);
        }
        if theta >= INF || cands.is_empty() {
            return (None, INF, 0.0);
        }
        let (mut best_w, mut chosen, mut step) = (0.0, NONE, 0.0);
        for &(k, r) in &cands {
            if r <= theta {
                let w = entries[k].1.abs();
                if best_w < w {
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

/// Solve starting from `start` (a primal feasible basis). `None` uses the all-slack basis with columns at their lower bounds.
pub fn solve(model: &LpModel, start: Option<&Basis>, opts: &SolveOptions) -> Result<SolveResult, SolveError> {
    let p = build(model);
    let (n, m) = (p.n, p.m);
    let basis = match start {
        Some(b) => b.clone(),
        None => {
            let mut b = Basis::slack(n, m);
            for j in 0..n {
                b.at_upper[j] = is_inf(p.lo[j]) && !is_inf(p.up[j]);
            }
            b
        }
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
        rows: Vec::new(),
        mirror: None,
    };
    s.rows = vec![Vec::new(); m];
    for j in 0..n {
        for &(i, _) in &s.p.cols[j] {
            s.rows[i].push(j);
        }
    }
    if opts.faithful_ftran && opts.faithful_positions {
        let (order, st, rows) = s.reference_order();
        let identity = order.iter().zip(&s.basic).all(|(a, b)| a == b);
        if std::env::var_os("ROSE_STATE").is_some() {
            eprintln!("INITIAL position order is the pivot order: {identity}");
        }
        if identity {
            s.mirror = Some(crate::ftran_order::RefStruct::fresh(m, &st, &rows));
        }
    }
    s.recompute_primal();
    s.recompute_dual();
    let mut trace = Vec::new();
    let mut iterations = 0usize;
    let mut phase1 = false;
    let mut farkas: Option<Vec<f64>> = None;
    let mut ray: Option<Vec<f64>> = None;
    let status = loop {
        if iterations >= opts.max_iterations {
            break Status::IterationLimit;
        }
        let due = match &opts.refactor_schedule {
            Some(sched) => sched.contains(&iterations),
            None => s.factor.num_updates() >= opts.max_updates,
        };
        if due {
            s.refactor()?;
        }
        let infeasible = s.has_infeasible_basic();
        if infeasible {
            phase1 = true;
            s.composite_dual(false);
        } else if phase1 {
            // Back in a feasible region: continue with the true costs (weights and framework are kept).
            phase1 = false;
            s.list.clear();
            s.level = 0;
            s.recompute_dual();
        }
        if std::env::var_os("ROSE_STATE").is_some() {
            let f = |v: &[f64]| v.iter().map(|x| format!("{x:e}")).collect::<Vec<_>>().join(",");
            eprintln!("STATE {{\"iter\":{iterations},\"infeasible\":{infeasible},\"basic\":{:?},\"d\":[{}],\"w\":[{}]}}", s.basic, f(&s.d), f(&s.w));
        }
        let q = if infeasible {
            match s.price() {
                Some(q) => q,
                None => {
                    // No improving column for the penalised cost: try the pure sum of infeasibilities, which either
                    // continues the search or proves infeasibility with the multipliers as certificate.
                    let y = s.composite_dual(true);
                    match s.price() {
                        Some(q) => q,
                        None => {
                            farkas = Some(y.iter().map(|v| -v).collect());
                            break Status::Infeasible;
                        }
                    }
                }
            }
        } else {
            match s.price() {
                Some(q) => q,
                None => break Status::Optimal,
            }
        };
        let mut alpha = vec![0.0; m];
        for &(i, a) in &p.cols[q] {
            alpha[i] = a;
        }
        s.factor.ftran(&mut alpha);
        let dir = if s.d[q] < 0.0 { 1.0 } else { -1.0 };
        let rows_q: Vec<usize> = p.cols[q].iter().map(|&(i, _)| i).collect();
        let spike: Option<Vec<usize>> = s.mirror.as_ref().map(|mi| mi.lstage(&rows_q));
        let entries: Vec<(usize, f64)> = match (&s.mirror, &spike) {
            (Some(mi), Some(sp)) => {
                let ord = mi.uorder(sp);
                let mut seen = vec![false; m];
                let mut e: Vec<(usize, f64)> = Vec::new();
                for &i in &ord {
                    if !seen[i] {
                        seen[i] = true;
                        if alpha[i].abs() > 1e-14 {
                            e.push((i, dir * alpha[i]));
                        }
                    }
                }
                for i in 0..m {
                    if !seen[i] && alpha[i].abs() > 1e-14 {
                        e.push((i, dir * alpha[i]));
                    }
                }
                e
            }
            _ => (0..m).filter(|&i| alpha[i].abs() > 1e-14).map(|i| (i, dir * alpha[i])).collect(),
        };
        let own = if !is_inf(p.lo[q]) && !is_inf(p.up[q]) { p.up[q] - p.lo[q] } else { INF };
        let mut leave_to_upper: Option<bool> = None;
        let (chosen, mut step, _pivot) = if infeasible {
            let (xs, basic, lo, up) = (&s.x, &s.basic, &p.lo, &p.up);
            let tol = opts.feasibility_tol;
            let info = |pos: usize| {
                let j = basic[pos];
                classify(xs[j], lo[j], up[j], tol)
            };
            let pick = phase1_ratio(tol, opts.pivot_tol, own, &entries, &info);
            if pick.entry.is_some() {
                leave_to_upper = Some(pick.to_upper);
            }
            (pick.entry, pick.step, pick.pivot)
        } else {
            let (xs, basic, lo, up) = (&s.x, &s.basic, &p.lo, &p.up);
            let down = |pos: usize| {
                let j = basic[pos];
                if is_inf(lo[j]) { INF } else { xs[j] - lo[j] }
            };
            let upd = |pos: usize| {
                let j = basic[pos];
                if is_inf(up[j]) { INF } else { up[j] - xs[j] }
            };
            ratio_test(opts.feasibility_tol, opts.pivot_tol, own, &entries, &down, &upd)
        };
        if chosen.is_none() && step >= INF {
            if infeasible {
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
        step = step.max(0.0);
        iterations += 1;
        // Move the basic variables and the entering variable.
        for &(i, a) in &entries {
            let j = s.basic[i];
            s.x[j] -= step * a;
        }
        s.x[q] += dir * step;
        let Some(k) = chosen else {
            // Bound flip: the entering variable travels to its opposite bound.
            s.at_upper[q] = !s.at_upper[q];
            s.x[q] = s.nonbasic_value(q);
            if opts.record_trace {
                trace.push(TraceEvent { iteration: iterations, entering: q, leaving: None, position: None, step });
            }
            continue;
        };
        let (r, a_r) = entries[k];
        let leave = s.basic[r];
        let to_upper = leave_to_upper.unwrap_or(a_r < 0.0);
        // Snap the leaving variable to the bound it reached.
        s.x[leave] = if to_upper { p.up[leave] } else { p.lo[leave] };
        let alpha_rq = alpha[r];

        // Pivot row alpha_r = e_r' B^-1 [A I] and the Devex cross-term vector tau = B^-T u.
        let mut rho = vec![0.0; m];
        rho[r] = 1.0;
        s.factor.btran(&mut rho);
        let mut u = vec![0.0; m];
        for i in 0..m {
            if s.in_r0[s.basic[i]] {
                u[i] = alpha[i];
            }
        }
        s.factor.btran(&mut u);
        let theta_d = s.d[q] / alpha_rq;
        let wq = s.w[q];
        let mut changed: Vec<usize> = Vec::new();
        for j in 0..n + m {
            if s.pos_of[j] != NONE || j == q {
                continue;
            }
            let arj = dot(&p.cols[j], &rho);
            if arj == 0.0 {
                continue;
            }
            let ratio = arj / alpha_rq;
            let vj = dot(&p.cols[j], &u);
            s.w[j] = s.w[j] - 2.0 * ratio * vj + ratio * ratio * wq;
            let floor = if opts.faithful_weights && !s.in_r0[j] { W_MIN } else { 1.0 };
            s.w[j] = s.w[j].max(floor);
            if !s.is_fixed(j) {
                s.d[j] -= theta_d * arj;
            }
            changed.push(j);
        }
        s.d[leave] = -theta_d;
        s.d[q] = 0.0;
        s.w[leave] = wq / (alpha_rq * alpha_rq);
        s.w[leave] = s.w[leave].max(if opts.faithful_weights && !s.in_r0[leave] { W_MIN } else { 1.0 });

        // Basis change.
        s.pos_of[leave] = NONE;
        s.at_upper[leave] = to_upper;
        s.basic[r] = q;
        s.pos_of[q] = r;
        s.factor.update(r, alpha);
        if let (Some(mi), Some(sp)) = (s.mirror.as_mut(), spike.as_ref()) {
            mi.update(r, sp);
        }
        changed.push(leave);
        let changed = s.pivot_row_order(&rho, &changed);
        s.update_list(&changed);
        if opts.record_trace {
            trace.push(TraceEvent { iteration: iterations, entering: q, leaving: Some(leave), position: Some(r), step });
        }
    };

    if status == Status::Optimal {
        s.recompute_primal();
        s.recompute_dual();
    }
    let sgn = if model.sense == Sense::Maximize { -1.0 } else { 1.0 };
    let mut y: Vec<f64> = s.basic.iter().map(|&j| p.cost[j]).collect();
    s.factor.btran(&mut y);
    let nm = model.num_cols();
    let x: Vec<f64> = s.x[..nm].to_vec();
    Ok(SolveResult {
        status,
        objective: model.objective(&x),
        x,
        dual: y.iter().map(|v| sgn * v).collect(),
        reduced_cost: s.d[..nm].iter().map(|v| sgn * v).collect(),
        basis: Basis { basic: s.basic.clone(), at_upper: s.at_upper.clone() },
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
        assert!((r.objective - 11.0).abs() < 1e-9, "objective {}", r.objective);
        assert!((r.x[0] - 3.0).abs() < 1e-9 && (r.x[1] - 1.0).abs() < 1e-9);
    }

    #[test]
    fn detects_unbounded() {
        let text = "NAME T\nROWS\n N obj\n G c1\nCOLUMNS\n x obj -1 c1 1\nRHS\n rhs c1 1\nENDATA\n";
        let model = mps::parse(text).unwrap();
        // x >= 1 is infeasible at the lower bound 0, so start from the basis where the slack is basic at x=1... use x basic.
        let start = Basis { basic: vec![0], at_upper: vec![false, false] };
        let r = solve(&model, Some(&start), &SolveOptions::default()).unwrap();
        assert_eq!(r.status, Status::Unbounded);
    }
}
