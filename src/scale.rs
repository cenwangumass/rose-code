//! Power-of-two scaling of the reference's default `ScaleFlag` (`reverse-engineer/specs/primal-simplex.md`, "Scaling").
//!
//! The reference scales the model before the simplex starts: `A' = R A C`, `c' = C c`, `b' = R b`, bounds `l' = l / C`
//! (the factors are powers of two, so the scaled data are exact), solves it, and unscales the returned solution.
//! The routine is the setup routine `0x116520` (mode 3); the factors come from an alternating column/row max-norm
//! equilibration with a damping of 2 per pass. [`compute`] reproduces it, [`solve`] runs the simplex on the scaled model.

use crate::basis::Basis;
use crate::model::{LpModel, INF};
use crate::simplex::{SolveError, SolveOptions, SolveResult};

/// Scale factors: `col` has one entry per structural column followed by one per range column (a row with
/// two distinct finite bounds, in row order; the reference adds these columns to the matrix and scales them like the others),
/// `row` one per row.
#[derive(Clone, Debug, PartialEq)]
pub struct Scaling {
    pub col: Vec<f64>,
    pub row: Vec<f64>,
    /// Objective scale (a power of two): the scaled costs are multiplied by it (see [`obj_scale`]).
    pub obj: f64,
}

const LOWEST: f64 = 7.105427357601002e-15; // 2^-47
const HIGHEST: f64 = 140737488355328.0; // 2^47
const PASSES: usize = 10;

fn is_inf(v: f64) -> bool {
    v.abs() >= INF
}

/// `scalbn(1, ilogb(v))` for `v > 0`: the largest power of two not above `v`.
fn pow2_floor(v: f64) -> f64 {
    if v > 0.5 && v < 2.0 {
        return 1.0;
    }
    if !v.is_finite() {
        return v;
    }
    let bits = v.to_bits();
    let e = ((bits >> 52) & 0x7ff) as i32;
    let exp = if e == 0 {
        63 - (bits & ((1u64 << 52) - 1)).leading_zeros() as i32 - 1074
    } else {
        e - 1023
    };
    f64::from_bits(((exp + 1023).max(1) as u64) << 52)
}

/// The target `1 / m` (1 for `m == 0`), rounded down to a power of two, moved from `old` by at most a factor `thr`.
/// Returns the new value and whether the move was limited.
fn damped(target: f64, old: f64, thr: f64) -> (f64, bool) {
    let p = if target > 0.5 && target < 2.0 {
        1.0
    } else {
        pow2_floor(target)
    };
    if p > thr * old {
        return (thr * old, true);
    }
    let lowest = old / thr;
    (lowest.max(p), lowest > p)
}

/// The reference's scale factors of `model` (ScaleFlag default, no objective scaling).
pub fn compute(model: &LpModel) -> Scaling {
    let ranges: Vec<(usize, f64)> = (0..model.num_rows())
        .filter(|&i| {
            !is_inf(model.row_lo[i])
                && !is_inf(model.row_up[i])
                && model.row_lo[i] != model.row_up[i]
        })
        .map(|i| (i, model.row_up[i] - model.row_lo[i]))
        .collect();
    compute_with_ranges(model, &ranges)
}

/// [`compute`] with the range columns given as `(row, upper bound)` in column order (the barrier's: RANGES order, first entry per row, an `E` row with
/// `R = 0` included, bound `|R|`).
pub fn compute_with_ranges(model: &LpModel, ranges: &[(usize, f64)]) -> Scaling {
    let (nm, m) = (model.num_cols(), model.num_rows());
    // Columns: structurals, then one unit column per range column.
    let mut cols: Vec<Vec<(usize, f64)>> = (0..nm)
        .map(|j| model.col(j).map(|(i, a)| (i, a.abs())).collect())
        .collect();
    let mut lo = model.col_lo.clone();
    let mut up = model.col_up.clone();
    for &(i, ub) in ranges {
        cols.push(vec![(i, 1.0)]);
        lo.push(0.0);
        up.push(ub);
    }
    let n = cols.len();
    let mut col = vec![1.0f64; n];
    let mut row = vec![1.0f64; m];

    // One pass over the columns. With `bounds` the target is `sqrt(1/max * max(1, |l|, |u|))` (finite bounds only).
    let col_pass = |col: &mut Vec<f64>, row: &[f64], thr: f64, bounds: bool| -> bool {
        let mut changed = false;
        for j in 0..n {
            let mut cm = 0.0f64;
            for &(i, a) in &cols[j] {
                cm = cm.max((a * row[i]).abs());
            }
            let mut v = if cm != 0.0 { 1.0 / cm } else { 1.0 };
            if bounds {
                let mut b = 1.0f64;
                if !is_inf(lo[j]) {
                    b = b.max(lo[j].abs());
                }
                if !is_inf(up[j]) {
                    b = b.max(up[j].abs());
                }
                v = (v * b).sqrt();
            }
            let (new, lim) = damped(v, col[j], thr);
            col[j] = new;
            changed |= lim;
        }
        changed
    };
    // One pass over the rows: the largest scaled entry of each row (a column with an infinite upper bound compresses an
    // entry above 1e8 to its square root, not below 1e8).
    let row_pass = |row: &mut Vec<f64>, col: &[f64], thr: f64| -> bool {
        let mut rm = vec![0.0f64; m];
        for j in 0..n {
            for &(i, a) in &cols[j] {
                let t = (a * col[j]).abs();
                if t > rm[i] {
                    rm[i] = t;
                    if up[j] >= INF && t > 1e8 {
                        rm[i] = t.sqrt().max(1e8);
                    }
                }
            }
        }
        let mut changed = false;
        for i in 0..m {
            if rm[i] > 0.0 {
                let (new, lim) = damped(1.0 / rm[i], row[i], thr);
                row[i] = new;
                changed |= lim;
            }
        }
        changed
    };

    let (mut pass, mut changed, mut thr) = (0usize, false, 2.0f64);
    loop {
        // The first three column passes use the bound-aware target.
        changed |= col_pass(&mut col, &row, thr, pass < 3);
        if pass > 0 && !changed {
            if pass < 3 {
                // A last plain column pass, after which the pass count is 3.
                changed = col_pass(&mut col, &row, thr, false);
                pass = 3;
                if !changed {
                    break;
                }
            } else {
                break;
            }
        }
        pass += 1;
        if pass >= PASSES {
            break;
        }
        thr = if pass == PASSES - 1 { 1e100 } else { 2.0 };
        changed = row_pass(&mut row, &col, thr);
    }
    let clamp = |x: f64| if x < LOWEST { LOWEST } else { x.min(HIGHEST) };
    let col: Vec<f64> = col.into_iter().map(clamp).collect();
    let obj = obj_scale(model, &col);
    Scaling {
        col,
        row: row.into_iter().map(clamp).collect(),
        obj,
    }
}

/// The objective scale of the setup routine (`0x11c803..0x11cbad`, `ObjScale` 0): from the largest and smallest nonzero
/// scaled cost `|c_j C_j|` it takes `t = min(100 / max, 1 / (10 g))` with `g = sqrt(max(1e-6, min) max)`; `t < 2` leaves the
/// costs alone, else they are multiplied by `2^ilogb(t)` (at most `2^9`). Only PEROLD and PILOT4 of the corpus get a factor (2).
fn obj_scale(model: &LpModel, col: &[f64]) -> f64 {
    let (mut mx, mut mn) = (0.0f64, f64::INFINITY);
    for j in 0..model.num_cols() {
        let c = (model.obj[j] * col[j]).abs();
        if c != 0.0 {
            mx = mx.max(c);
            mn = mn.min(c);
        }
    }
    if mx == 0.0 {
        return 1.0;
    }
    let g = (mn.max(1e-6) * mx).sqrt();
    let t = (100.0 / mx).min(1.0 / (10.0 * g));
    if t < 2.0 {
        1.0
    } else if t > 1000.0 {
        512.0
    } else {
        pow2_floor(t)
    }
}

fn scale_bound(v: f64, f: f64) -> f64 {
    if is_inf(v) {
        v
    } else {
        v * f
    }
}

impl Scaling {
    /// The scaled model, or `None` when a range column's scale is not the inverse of its row's (the solver builds that
    /// column as a unit column; the reference keeps it as scaled data, which this model cannot hold).
    pub fn apply(&self, model: &LpModel) -> Option<LpModel> {
        let (nm, m) = (model.num_cols(), model.num_rows());
        let mut k = nm;
        for i in 0..m {
            let (rl, ru) = (model.row_lo[i], model.row_up[i]);
            if !is_inf(rl) && !is_inf(ru) && rl != ru {
                if self.col[k] * self.row[i] != 1.0 {
                    return None;
                }
                k += 1;
            }
        }
        let mut s = model.clone();
        for j in 0..nm {
            let c = self.col[j];
            s.obj[j] *= c;
            s.obj[j] *= self.obj;
            s.col_lo[j] = scale_bound(model.col_lo[j], 1.0 / c);
            s.col_up[j] = scale_bound(model.col_up[j], 1.0 / c);
            for p in model.col_start[j]..model.col_start[j + 1] {
                s.val[p] = model.val[p] * self.row[model.row_idx[p]] * c;
            }
        }
        for i in 0..m {
            s.row_lo[i] = scale_bound(model.row_lo[i], self.row[i]);
            s.row_up[i] = scale_bound(model.row_up[i], self.row[i]);
            if let Some(r) = s.row_range.get_mut(i) {
                *r *= self.row[i];
            }
        }
        Some(s)
    }

    /// The scaled data of every row and structural column without the check of [`Scaling::apply`]: the barrier builds a range column itself
    /// with the scaled coefficient `R_i C_r` and bound `range / C_r` (`range_entries` hold the entries multiplied by `R_i`).
    pub fn apply_loose(&self, model: &LpModel) -> LpModel {
        let (nm, m) = (model.num_cols(), model.num_rows());
        let mut s = model.clone();
        let mut ranged = vec![false; m];
        for &(i, _) in &model.range_entries {
            ranged[i] = true;
        }
        for j in 0..nm {
            let c = self.col[j];
            s.obj[j] *= c;
            s.obj[j] *= self.obj;
            s.col_lo[j] = scale_bound(model.col_lo[j], 1.0 / c);
            s.col_up[j] = scale_bound(model.col_up[j], 1.0 / c);
            for p in model.col_start[j]..model.col_start[j + 1] {
                s.val[p] = model.val[p] * self.row[model.row_idx[p]] * c;
            }
        }
        for i in 0..m {
            // A file equality's RHS remains finite beyond the variable-bound infinity sentinel.
            let equality = model.row_kind.get(i) == Some(&b'E')
                && !ranged[i]
                && model.row_lo[i] == model.row_up[i]
                && model.row_lo[i].is_finite();
            s.row_lo[i] = if equality {
                model.row_lo[i] * self.row[i]
            } else {
                scale_bound(model.row_lo[i], self.row[i])
            };
            s.row_up[i] = if equality {
                model.row_up[i] * self.row[i]
            } else {
                scale_bound(model.row_up[i], self.row[i])
            };
            if let Some(r) = s.row_range.get_mut(i) {
                *r *= self.row[i];
            }
        }
        for e in s.range_entries.iter_mut() {
            e.1 *= self.row[e.0];
        }
        s
    }

    /// Maps a result of the scaled model back: `x = C x'`, `y = R y'`, `d = d' / C`; the objective, the basis and the
    /// iteration data are unchanged.
    pub fn unscale(&self, model: &LpModel, mut r: SolveResult) -> SolveResult {
        let nm = model.num_cols();
        for j in 0..nm {
            r.x[j] *= self.col[j];
            r.reduced_cost[j] /= self.col[j] * self.obj;
        }
        for (i, y) in r.dual.iter_mut().enumerate() {
            *y *= self.row[i] / self.obj;
        }
        if let Some(f) = r.farkas.as_mut() {
            for (i, y) in f.iter_mut().enumerate() {
                *y *= self.row[i];
            }
        }
        if let Some(ray) = r.ray.as_mut() {
            for j in 0..nm {
                ray[j] *= self.col[j];
            }
        }
        r.objective = crate::dual::unscaled_objective(model, &r.x, &self.col);
        r
    }
}

thread_local! {
    /// The unscaled model of the running scaled solve with its column and row factors: the primal's end-of-solve test and its dual clean-up pass use it
    /// (the reference's `0xcff40` measures the primal infeasibility of the unscaled problem and `0x1bbd80` runs the dual on it).
    pub static ORIG: std::cell::RefCell<Option<(LpModel, Vec<f64>, Vec<f64>)>> = const { std::cell::RefCell::new(None) };
}

thread_local! {
    /// The dual clean-up pass ran inside the scaled primal: the reference makes no unscaled second primal pass after it.
    pub static PASS_RAN: std::cell::Cell<bool> = const { std::cell::Cell::new(false) };
    /// The unscaled (primal violation, dual infeasibility) of the scaled primal's optimum (`0xcff40`, `r9+0x90`, `r9+0xb0`).
    pub static END_INF: std::cell::Cell<Option<(f64, f64)>> = const { std::cell::Cell::new(None) };
}

/// Solves the scaled model and returns the solution of `model`.
pub fn solve(
    model: &LpModel,
    start: Option<&Basis>,
    opts: &SolveOptions,
) -> Result<SolveResult, SolveError> {
    let sc = compute(model);
    if crate::dump::on("scale") {
        crate::dump::emit("scale", serde_json::json!({"col": sc.col, "row": sc.row}));
    }
    let inner = SolveOptions {
        scaling: false,
        obj_scale: sc.obj,
        ..opts.clone()
    };
    match sc.apply(model) {
        Some(scaled) => {
            let inner = SolveOptions {
                scaled: true,
                ..inner
            };
            crate::simplex::RAY_ENTER.with(|c| c.set(None));
            ORIG.with(|c| *c.borrow_mut() = Some((model.clone(), sc.col.clone(), sc.row.clone())));
            END_INF.with(|c| c.set(None));
            PASS_RAN.with(|c| c.set(false));
            let mut r = crate::simplex::solve(&scaled, start, &inner);
            ORIG.with(|c| *c.borrow_mut() = None);
            let mut limited_pass: Option<(usize, Vec<_>)> = None;
            if opts.pricing0 || END_INF.with(|c| c.get()).is_some() {
                if let (Ok(r1), Some((pv, dv)), false) = (
                    r.as_ref(),
                    END_INF.with(|c| c.take()),
                    PASS_RAN.with(|c| c.get()),
                ) {
                    if r1.status == crate::simplex::Status::Optimal && pv <= 1e-6 && dv >= 1e-6 {
                        let n1 = r1.iterations;
                        // The second pass runs under an iteration limit (`0xd0233`, mode `ctx`+`0x2514` = -1): total iterations `T` so far plus a fifth of
                        // `max(T, cols / 3 + rows)`; the check after each iteration (`0x11f430`) ends the pass with status 7.
                        let t = n1 as f64;
                        let limit =
                            t + 0.2 * t.max((model.num_cols() / 3 + model.num_rows()) as f64);
                        let allowed = (1..).find(|&k| t + k as f64 >= limit).unwrap_or(1);
                        let second = SolveOptions {
                            scaled: false,
                            obj_scale: 1.0,
                            max_iterations: inner.max_iterations.min(allowed),
                            ..inner.clone()
                        };
                        let basis = r1.basis.clone();
                        let readout = crate::dual::quad_result_state();
                        crate::simplex::SECOND_PASS.with(|c| c.set(true));
                        let r2 = crate::simplex::solve(model, Some(&basis), &second);
                        crate::simplex::SECOND_PASS.with(|c| c.set(false));
                        if let Ok(mut r2) = r2 {
                            r2.iterations += n1;
                            // The diagnostic stream keeps the pivots of the scaled pass.
                            r2.trace.splice(0..0, r1.trace.iter().cloned());
                            // Status 7 (the iteration limit): the driver restores the scaled problem and the basis saved before the pass (`0xcf8d0`, `0xcfa8e`); the pass's
                            // pivot is undone but stays counted, and the dual pass that follows makes no pivot (free-shift witnesses): the result is the scaled pass's.
                            if r2.status == crate::simplex::Status::IterationLimit
                                && allowed < inner.max_iterations
                            {
                                limited_pass = Some((r2.iterations, std::mem::take(&mut r2.trace)));
                                crate::dual::set_quad_result_state(readout);
                            } else if r2.status == crate::simplex::Status::Optimal {
                                let dopts = SolveOptions {
                                    method: crate::simplex::Method::Dual,
                                    record_trace: true,
                                    ..SolveOptions::default()
                                };
                                crate::dual::CLEANUP_PASS.with(|c| {
                                    c.set(Some(crate::simplex::ACC_THR_PASS.with(|a| a.get())))
                                });
                                let via = crate::dual::solve(model, Some(&r2.basis), &dopts);
                                crate::dual::CLEANUP_PASS.with(|c| c.set(None));
                                if let Ok(r3) = via {
                                    if r3.iterations > 0
                                        || r3.status == crate::simplex::Status::Optimal
                                    {
                                        let mut r3 = r3;
                                        r3.iterations += r2.iterations;
                                        r3.trace.splice(0..0, r2.trace.iter().cloned());
                                        return Ok(r3);
                                    }
                                }
                            }
                            if limited_pass.is_none() {
                                return Ok(r2);
                            }
                        }
                    }
                }
            }
            if let (Some((iterations, trace)), Ok(r1)) = (limited_pass, r.as_mut()) {
                r1.iterations = iterations;
                r1.trace = trace;
            }
            let r = r?;
            let mut out = sc.unscale(model, r);
            if let (Some(ray), Some((q, n))) = (
                out.ray.as_mut(),
                crate::simplex::RAY_ENTER.with(|c| c.take()),
            ) {
                // Unit step of the entering variable in the user's frame: a slack's scaled value is `R s`, a structural's is `x / C`.
                let g = if q >= n {
                    sc.row[q - n]
                } else {
                    1.0 / sc.col[q]
                };
                ray.iter_mut().for_each(|v| *v *= g);
            }
            // The dual's FarkasDual: the reference enters the unit vector of the leaving position in the user's frame, so the wrapper `0x1b32e0`
            // divides it by the row factor of a leaving slack (multiplies it by the column factor of a leaving structural) before the BTRAN;
            // the result is `R e_r' B'^-1` times that factor.
            if opts.method == crate::simplex::Method::Dual {
                if let (Some(f), Some((leave, n))) = (
                    out.farkas.as_mut(),
                    crate::dual::FARKAS_LEAVE.with(|c| c.take()),
                ) {
                    let g = if leave >= n {
                        1.0 / sc.row[leave - n]
                    } else {
                        sc.col[leave]
                    };
                    f.iter_mut().for_each(|y| *y *= g);
                }
            }
            Ok(out)
        }
        None => crate::simplex::solve(model, start, &inner),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn pow2_floor_matches_ilogb() {
        for (v, p) in [
            (0.3, 0.25),
            (0.5, 0.5),
            (0.49, 0.25),
            (2.0, 2.0),
            (3.9, 2.0),
            (4.0, 4.0),
            (1e-3, crate::real::pow2(-10)),
            (1.5, 1.0),
            (1e5, 65536.0),
        ] {
            assert_eq!(pow2_floor(v), p, "{v}");
        }
    }
}
