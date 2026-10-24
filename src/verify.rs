//! Independent checks of solutions and certificates against the original model. Nothing here uses the simplex
//! internals; every quantity is recomputed from `LpModel` and the returned vectors.

use crate::model::{LpModel, Sense, INF};
use crate::simplex::{SolveResult, Status};

fn fin(v: f64) -> bool {
    v.abs() < INF
}

#[derive(Clone, Debug, Default)]
pub struct Report {
    /// Largest violation of a column or row bound.
    pub primal_violation: f64,
    /// Largest dual sign violation (row duals and reduced costs).
    pub dual_violation: f64,
    /// Largest `|multiplier| * slack` over rows and columns.
    pub complementarity: f64,
    /// `|reported objective - c'x - offset|`.
    pub objective_error: f64,
    /// `|primal objective - dual objective|`, with the dual objective rebuilt from `y`, `d` and the bounds.
    pub duality_gap: f64,
    /// Largest `|d_j - (c_j - a_j'y)|`.
    pub reduced_cost_error: f64,
}

impl Report {
    pub fn max_error(&self) -> f64 {
        self.primal_violation
            .max(self.dual_violation)
            .max(self.complementarity)
            .max(self.objective_error)
            .max(self.duality_gap)
            .max(self.reduced_cost_error)
    }
}

/// Check an `Optimal` result: primal feasibility, dual feasibility, complementarity, objective and gap.
/// Duals and reduced costs are in the original sense (for maximization the signs are flipped internally).
pub fn check_optimal(model: &LpModel, r: &SolveResult) -> Report {
    let (n, m) = (model.num_cols(), model.num_rows());
    let sgn = if model.sense == Sense::Maximize {
        -1.0
    } else {
        1.0
    };
    let act = model.activities(&r.x);
    let mut rep = Report::default();
    for j in 0..n {
        if fin(model.col_lo[j]) {
            rep.primal_violation = rep.primal_violation.max(model.col_lo[j] - r.x[j]);
        }
        if fin(model.col_up[j]) {
            rep.primal_violation = rep.primal_violation.max(r.x[j] - model.col_up[j]);
        }
    }
    for i in 0..m {
        if fin(model.row_lo[i]) {
            rep.primal_violation = rep.primal_violation.max(model.row_lo[i] - act[i]);
        }
        if fin(model.row_up[i]) {
            rep.primal_violation = rep.primal_violation.max(act[i] - model.row_up[i]);
        }
    }
    rep.objective_error = (r.objective - model.objective(&r.x)).abs();

    // Work in minimization form: c' = sgn c, y' = sgn y, d' = sgn d.
    let y: Vec<f64> = r.dual.iter().map(|v| sgn * v).collect();
    let d: Vec<f64> = r.reduced_cost.iter().map(|v| sgn * v).collect();
    let mut dual_obj = sgn * model.offset;
    for j in 0..n {
        let aty: f64 = model.col(j).map(|(i, a)| a * y[i]).sum();
        rep.reduced_cost_error = rep
            .reduced_cost_error
            .max((d[j] - (sgn * model.obj[j] - aty)).abs());
        // d > 0: pushes against the lower bound; d < 0: against the upper bound.
        let (side, slack) = if d[j] > 0.0 {
            (model.col_lo[j], r.x[j] - model.col_lo[j])
        } else {
            (model.col_up[j], model.col_up[j] - r.x[j])
        };
        if fin(side) {
            rep.complementarity = rep.complementarity.max(d[j].abs() * slack.abs());
            dual_obj += d[j] * side;
        } else if d[j] != 0.0 {
            rep.dual_violation = rep.dual_violation.max(d[j].abs());
        }
    }
    for i in 0..m {
        let (side, slack) = if y[i] > 0.0 {
            (model.row_lo[i], act[i] - model.row_lo[i])
        } else {
            (model.row_up[i], model.row_up[i] - act[i])
        };
        if fin(side) {
            rep.complementarity = rep.complementarity.max(y[i].abs() * slack.abs());
            dual_obj += y[i] * side;
        } else if y[i] != 0.0 {
            rep.dual_violation = rep.dual_violation.max(y[i].abs());
        }
    }
    rep.duality_gap = (sgn * model.objective(&r.x) - dual_obj).abs();
    rep
}

/// Margin by which `y` proves infeasibility: positive if `y'A x` over the column box cannot meet `y'r` over the
/// row ranges. The multipliers are scaled to unit max-norm first.
pub fn farkas_margin(model: &LpModel, y: &[f64]) -> f64 {
    let (a, b) = farkas_margins(model, y);
    a.max(b)
}

/// The two forms of `farkas_margin`: `(a, b)` with `a` = the row side from below exceeds the column side from above (`y'r > y'Ax`, the reference's
/// `FarkasDual` convention: `y >= 0` on rows with only a lower side, `<= 0` on rows with only an upper side) and `b` = the column side from below
/// exceeds the row side from above (the opposite signs, which rose's primal and the dual method return).
pub fn farkas_margins(model: &LpModel, y: &[f64]) -> (f64, f64) {
    let scale = y.iter().fold(0.0f64, |a, v| a.max(v.abs()));
    if scale == 0.0 {
        return (0.0, 0.0);
    }
    let y: Vec<f64> = y.iter().map(|v| v / scale).collect();
    let (mut gl, mut gu) = (0.0f64, 0.0f64);
    for j in 0..model.num_cols() {
        let g: f64 = model.col(j).map(|(i, a)| a * y[i]).sum();
        // Rounding noise left by cancelling terms is not a coefficient: with an infinite bound it would make the combined range unbounded.
        if g.abs() <= 1e-12 {
            continue;
        }
        // Bound attaining the minimum / maximum of g*x_j.
        let (lo, hi) = if g > 0.0 {
            (model.col_lo[j], model.col_up[j])
        } else {
            (model.col_up[j], model.col_lo[j])
        };
        if fin(lo) {
            gl += g * lo
        } else {
            gl = -INF
        }
        if fin(hi) {
            gu += g * hi
        } else {
            gu = INF
        }
    }
    let (mut rl, mut ru) = (0.0f64, 0.0f64);
    for i in 0..model.num_rows() {
        let (lo, up) = (model.row_lo[i], model.row_up[i]);
        if y[i] == 0.0 {
            continue;
        }
        let (a, b) = if y[i] > 0.0 { (lo, up) } else { (up, lo) };
        if fin(a) {
            rl += y[i] * a
        } else {
            rl = -INF
        }
        if fin(b) {
            ru += y[i] * b
        } else {
            ru = INF
        }
    }
    (rl - gu, gl - ru)
}

/// Margin of an unbounded-ray certificate: `Some(-c'r)` (positive means the objective improves along the ray) if
/// `r` is a recession direction of the feasible region, `None` if it violates a bound by more than `tol`.
pub fn ray_improvement(model: &LpModel, ray: &[f64], tol: f64) -> Option<f64> {
    let scale = ray.iter().fold(0.0f64, |a, v| a.max(v.abs()));
    if scale == 0.0 {
        return None;
    }
    let r: Vec<f64> = ray.iter().map(|v| v / scale).collect();
    for j in 0..model.num_cols() {
        if (r[j] < -tol && fin(model.col_lo[j])) || (r[j] > tol && fin(model.col_up[j])) {
            return None;
        }
    }
    let act = model.activities(&r);
    for i in 0..model.num_rows() {
        if (act[i] < -tol && fin(model.row_lo[i])) || (act[i] > tol && fin(model.row_up[i])) {
            return None;
        }
    }
    let sgn = if model.sense == Sense::Maximize {
        -1.0
    } else {
        1.0
    };
    Some(-sgn * model.obj.iter().zip(&r).map(|(c, v)| c * v).sum::<f64>())
}

/// Full check of a result of any status against the model; `Err` describes the first failed condition.
pub fn check(model: &LpModel, r: &SolveResult, tol: f64) -> Result<(), String> {
    match r.status {
        Status::Optimal => {
            let rep = check_optimal(model, r);
            let scale = 1.0 + r.objective.abs();
            if rep.primal_violation > tol
                || rep.dual_violation > tol
                || rep.reduced_cost_error > tol
                || rep.complementarity > tol * scale
                || rep.duality_gap > tol * scale
                || rep.objective_error > 1e-9 * scale
            {
                return Err(format!("{rep:?}"));
            }
            Ok(())
        }
        Status::Infeasible => match &r.farkas {
            Some(y) if farkas_margin(model, y) > 1e-9 => Ok(()),
            _ => Err("missing or invalid Farkas certificate".into()),
        },
        Status::Unbounded => {
            let rep = check_optimal_primal(model, &r.x);
            if rep > tol {
                return Err(format!(
                    "unbounded result point violates bounds by {rep:.3e}"
                ));
            }
            match r
                .ray
                .as_ref()
                .and_then(|ray| ray_improvement(model, ray, tol))
            {
                Some(g) if g > 1e-9 => Ok(()),
                _ => Err("missing or invalid unbounded ray".into()),
            }
        }
        s => Err(format!("status {s:?}")),
    }
}

fn check_optimal_primal(model: &LpModel, x: &[f64]) -> f64 {
    let act = model.activities(x);
    let mut v = 0.0f64;
    for j in 0..model.num_cols() {
        if fin(model.col_lo[j]) {
            v = v.max(model.col_lo[j] - x[j]);
        }
        if fin(model.col_up[j]) {
            v = v.max(x[j] - model.col_up[j]);
        }
    }
    for i in 0..model.num_rows() {
        if fin(model.row_lo[i]) {
            v = v.max(model.row_lo[i] - act[i]);
        }
        if fin(model.row_up[i]) {
            v = v.max(act[i] - model.row_up[i]);
        }
    }
    v
}
