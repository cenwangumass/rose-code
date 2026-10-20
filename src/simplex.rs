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
        }
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Status {
    Optimal,
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
}

#[derive(Debug)]
pub enum SolveError {
    InvalidBasis(String),
    /// The starting point is primal infeasible; Phase 1 is not implemented yet.
    InfeasibleStart { max_violation: f64 },
    SingularBasis,
}

impl std::fmt::Display for SolveError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            SolveError::InvalidBasis(s) => write!(f, "invalid basis: {s}"),
            SolveError::InfeasibleStart { max_violation } => write!(f, "starting basis is infeasible (violation {max_violation:.3e})"),
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

fn build(model: &LpModel) -> Internal {
    let (n, m) = (model.num_cols(), model.num_rows());
    let mut cols: Vec<Column> = (0..n).map(|j| model.col(j).collect()).collect();
    let mut lo = model.col_lo.clone();
    let mut up = model.col_up.clone();
    let mut rhs = vec![0.0; m];
    for i in 0..m {
        let (rl, ru) = (model.row_lo[i], model.row_up[i]);
        let (sign, r, slo, sup) = match (is_inf(rl), is_inf(ru)) {
            (false, false) if rl == ru => (1.0, ru, 0.0, 0.0),
            (false, false) => (1.0, ru, 0.0, ru - rl),
            (true, false) => (1.0, ru, 0.0, INF),
            (false, true) => (-1.0, rl, 0.0, INF),
            (true, true) => (1.0, 0.0, -INF, INF),
        };
        cols.push(vec![(i, sign)]);
        rhs[i] = r;
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

    fn refactor(&mut self) -> Result<(), SolveError> {
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
    };
    s.recompute_primal();
    s.recompute_dual();
    let mut worst = 0.0f64;
    for &j in &s.basic {
        if !is_inf(p.lo[j]) {
            worst = worst.max(p.lo[j] - s.x[j]);
        }
        if !is_inf(p.up[j]) {
            worst = worst.max(s.x[j] - p.up[j]);
        }
    }
    if worst > opts.feasibility_tol {
        return Err(SolveError::InfeasibleStart { max_violation: worst });
    }

    let mut trace = Vec::new();
    let mut iterations = 0usize;
    let status = loop {
        if iterations >= opts.max_iterations {
            break Status::IterationLimit;
        }
        if s.factor.num_updates() >= opts.max_updates {
            s.refactor()?;
        }
        let Some(q) = s.price() else {
            break Status::Optimal;
        };
        let mut alpha = vec![0.0; m];
        for &(i, a) in &p.cols[q] {
            alpha[i] = a;
        }
        s.factor.ftran(&mut alpha);
        let dir = if s.d[q] < 0.0 { 1.0 } else { -1.0 };
        let entries: Vec<(usize, f64)> = (0..m).filter(|&i| alpha[i].abs() > 1e-14).map(|i| (i, dir * alpha[i])).collect();
        let own = if !is_inf(p.lo[q]) && !is_inf(p.up[q]) { p.up[q] - p.lo[q] } else { INF };
        let (chosen, mut step, _pivot) = {
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
        let to_upper = a_r < 0.0;
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
            s.w[j] = (s.w[j] - 2.0 * ratio * vj + ratio * ratio * wq).max(1.0);
            if !s.is_fixed(j) {
                s.d[j] -= theta_d * arj;
            }
            changed.push(j);
        }
        s.d[leave] = -theta_d;
        s.d[q] = 0.0;
        s.w[leave] = (wq / (alpha_rq * alpha_rq)).max(1.0);

        // Basis change.
        s.pos_of[leave] = NONE;
        s.at_upper[leave] = to_upper;
        s.basic[r] = q;
        s.pos_of[q] = r;
        s.factor.update(r, alpha);
        changed.push(leave);
        changed.sort_unstable();
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
    let x: Vec<f64> = s.x[..n].to_vec();
    Ok(SolveResult {
        status,
        objective: model.objective(&x),
        x,
        dual: y.iter().map(|v| sgn * v).collect(),
        reduced_cost: s.d[..n].iter().map(|v| sgn * v).collect(),
        basis: Basis { basic: s.basic.clone(), at_upper: s.at_upper.clone() },
        iterations,
        refactorizations: s.refactorizations,
        trace,
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
