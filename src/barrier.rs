//! Barrier iteration of the reference (`reverse-engineer/specs/barrier-iteration.md`; port of `tools/barrier_iter.py`): conversion to the
//! solved standard form, the starting point, residuals and scalars, the predictor / corrector decisions, the step-back update and the
//! termination test `0x2e2fe0`. LP only: free and lower-unbounded columns are converted (implied bounds, regularised D, refinement); maximization, dense columns and cones (reported as [`BarrierError::Unsupported`]).
//!
//! The linear algebra behind the iteration (factor of `A D A'`, the solves of `0x305f10`) is the [`NormalSolver`] trait. The reference's
//! supernodal Cholesky and its solves are not recovered bit for bit yet, so two implementations exist:
//! - [`ReplaySolver`]: returns captured solve outputs, so that everything else can be compared bitwise per iteration
//!   (`examples/barrier_replay.rs`, `tools/compare.py trace --profile barrier-controlled --rose`);
//! - [`DenseSolver`]: a plain f64 dense `LDL'` for tests at tolerance. **Not bit-exact**: its results differ from the reference in the
//!   last bits, so the iteration it drives is only as close to the reference as the solves are.
//!
//! No fused multiply-add anywhere: the reference's kernels use separate `mulsd/addsd` (SSE2 pairs, lane order in the helpers below).
//! Confidence labels as in the spec: C = verified bit for bit against captures, H = read from code only.
use crate::dump;
use crate::model::{LpModel, Sense, INF};
use serde_json::{json, Value};

#[derive(Debug, Clone, PartialEq)]
pub enum BarrierError {
    /// A model feature outside the recovered LP case (named).
    Unsupported(String),
    /// A solver failure (singular system, replay exhausted, ...).
    Solver(String),
}

impl std::fmt::Display for BarrierError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            BarrierError::Unsupported(s) => write!(f, "unsupported: {s}"),
            BarrierError::Solver(s) => write!(f, "solver: {s}"),
        }
    }
}

impl std::error::Error for BarrierError {}

fn unsupported<T>(s: impl Into<String>) -> Result<T, BarrierError> {
    Err(BarrierError::Unsupported(s.into()))
}

// ---------------------------------------------------------------------------------------------------------------------------------
// Standard form (conversion `0x31af70`, `0x3548a0`).

/// The solved problem `min c'x, A x + S s = b, 0 <= x (<= ub)` in the reference's layout. Variables: the model's columns in file order
/// (`n_struct`), the range columns of the ranged rows (appended in RANGES order, `n2 = n_struct + ranges`), then one slack per row (`m`, flag 0 and no
/// complementarity pair for an '=' row). Rows are in the supernodal elimination order (`SYM` perm), `A` columns keep their entries in ascending permuted row
/// order. A ranged row is an equality row plus its range column (`MPS_Rg<row>`, cost 0, bounds `[0, |R|]`, coefficient +1 for `L` and `E` rows, -1 for `G`);
/// a '>' row is negated (row and right-hand side, sense kept). Lower bounds are shifted to 0 (`b -= A lo`, `const = c lo`, `ub -= lo`); upper bounds
/// implied by the rows are dropped (`classify_bounds`).
#[derive(Debug, Clone)]
pub struct StdForm {
    pub m: usize,
    /// Columns of the solved problem: the model's, then the range columns.
    pub n2: usize,
    /// The model's own columns (the first `n_struct` of `n2`).
    pub n_struct: usize,
    pub col_start: Vec<usize>,
    pub row_idx: Vec<usize>,
    pub val: Vec<f64>,
    /// `b'='`, `b'<'` or `b'>'` per row (permuted order).
    pub sense: Vec<u8>,
    /// Flag per variable (`n2 + m`): bit 0 = lower bound 0 (pair x, z), bit 1 = upper bound (pair u = ub - x, v), `Q+0x640`.
    pub flags: Vec<u8>,
    /// Right-hand side `Q+0x610` (permuted rows).
    pub b: Vec<f64>,
    /// Right-hand side in file row order before the lower-bound shift, with the row sign (input of the setup scalars `0x2e4f10`: all zero on BORE3D).
    pub b_file: Vec<f64>,
    /// Costs of the columns, `Q+0x608`.
    pub c: Vec<f64>,
    /// Upper bounds after the shift (`Q+0x600`; [`INF`] where there is none).
    pub ub: Vec<f64>,
    /// `Q+0x630`, `Q+0x634`: the range of columns that holds the upper-bounded ones.
    pub ub_lo: usize,
    pub ub_hi: usize,
    pub offset: f64,
    /// -1 for a maximization model (solved as the minimization of the negated objective), else 1.
    pub max_sign: f64,
    /// `c . lower` of the shifted columns.
    pub cshift: f64,
    /// Lower bound removed from each column.
    pub lower: Vec<f64>,
    /// File row of each position (`inv[k]`).
    pub row_of: Vec<usize>,
    /// Number of complementarity pairs (`Q+0x210`): lower bounds, the upper bounds present before the implied ones were dropped, and the inequality slacks.
    pub ncomp: usize,
    /// `Q+0x648` ([`uniform_scale`]): per column the common value of its entries, 0 for the others; `None` when the reference frees the vector.
    pub scale: Option<Vec<f64>>,
    /// Dense columns of the model conversion (`P+0x218`, [`crate::order::dense_rule`]), ascending. A dense column is a node of the factor (`m + dense_k`),
    /// not a term of `A D A'`.
    pub dense: Vec<usize>,
    /// `k` of each dense column: its factor node is `m + k` (`SYM perm[m + r] - m`).
    pub dense_k: Vec<usize>,
    /// The release profile's scale factors (`reverse-engineer/specs/barrier-release.md`, "Scaling") the data were scaled with, if any.
    pub scaling: Option<crate::scale::Scaling>,
    /// File row of the model for each scaled row factor and column of each range column in `scaling.col`: the range column `k` (RANGES order) is `scaling.col[n_struct + range_rank[k]]`.
    pub range_rank: Vec<usize>,
}

/// `Q+0x648` (allocated at `0x2e5533`): for every column with an entry whose entries all lie within `1e-10` of the first, `s_j` = the first value
/// (0 otherwise). The vector is freed (`None`) when `2 nc > sum over the uniform columns of count^2` (`0x2e56e8`; DEGEN3 keeps it).
pub fn uniform_scale(n2: usize, col_start: &[usize], val: &[f64]) -> Option<Vec<f64>> {
    let mut s = vec![0.0; n2];
    let mut tot = 0.0f64;
    for j in 0..n2 {
        let (a, e) = (col_start[j], col_start[j + 1]);
        if e == a {
            continue;
        }
        let v0 = val[a];
        if (a + 1..e).all(|k| (v0 - val[k]).abs() <= 1e-10) {
            s[j] = v0;
            let c = (e - a) as f64;
            tot = tot + c * c;
        }
    }
    if 2.0 * n2 as f64 > tot {
        None
    } else {
        Some(s)
    }
}

impl StdForm {
    pub fn n(&self) -> usize {
        self.n2 + self.m
    }

    /// Number of complementarity pairs (`Q+0x210`).
    pub fn ncomp(&self) -> usize {
        self.ncomp
    }

    pub fn has_ub(&self) -> bool {
        self.ub_lo < self.ub_hi
    }

    /// `(A x)_i` over the columns only, accumulated column by column (`a x_j + out_i`).
    pub(crate) fn ax_cols(&self, x: &[f64]) -> Vec<f64> {
        let mut out = vec![0.0; self.m];
        for j in 0..self.n2 {
            for k in self.col_start[j]..self.col_start[j + 1] {
                out[self.row_idx[k]] = self.val[k] * x[j] + out[self.row_idx[k]];
            }
        }
        out
    }

    /// `Af x = A x + S x_slack` (slack added for the non-'=' rows).
    fn af(&self, x: &[f64]) -> Vec<f64> {
        let mut out = self.ax_cols(x);
        for i in 0..self.m {
            if self.sense[i] != b'=' {
                out[i] = out[i] + x[self.n2 + i];
            }
        }
        out
    }

    /// `(A' y)_j` of column `j`, entries summed sequentially from 0. A scaled column (`Q+0x648[j] != 0`) sums `y_i s_j` in two accumulators (even and
    /// odd positions, added at the end), then the last entry of an odd count.
    pub(crate) fn aty_col(&self, j: usize, y: &[f64]) -> f64 {
        if let Some(sc) = self.scale.as_ref().map(|v| v[j]).filter(|&v| v != 0.0) {
            let (a, cn) = (self.col_start[j], self.col_start[j + 1] - self.col_start[j]);
            let (mut ev, mut od) = (0.0, 0.0);
            for k in 0..cn - (cn & 1) {
                if k % 2 == 0 {
                    ev = ev + y[self.row_idx[a + k]] * sc;
                } else {
                    od = od + y[self.row_idx[a + k]] * sc;
                }
            }
            let mut t = ev + od;
            if cn & 1 == 1 {
                t = t + y[self.row_idx[a + cn - 1]] * sc;
            }
            return t;
        }
        let mut t = 0.0;
        for k in self.col_start[j]..self.col_start[j + 1] {
            t = t + y[self.row_idx[k]] * self.val[k];
        }
        t
    }
}

fn is_inf(v: f64) -> bool {
    v.abs() >= INF
}

/// Implied bounds of the free columns (`0x3283f0`, port of the first part of `classify_bounds` of `tools/barrier_model.py`; only when a column has neither bound).
/// Two sweeps over the rows into work copies `lo`, `hi` of the bounds. For a row and each `s` in (-1 ('=' rows only), +1): `t_k = s a_k`;
/// `acc` = sum in entry order of `t_k lo_k` (`t_k > 0`, `lo_k > -INF`) or `t_k hi_k` (`t_k <= 0`, `hi_k < INF`), the other entries counted as infinite;
/// `slack = s rhs - acc`, skipped if `slack < -1e-6`. No infinite entry: per entry `t > 1e-6`: `hi_k = min(hi_k, slack/t + lo_k)`; `t < -1e-6`:
/// `lo_k = max(lo_k, slack/t + hi_k)`. One infinite entry (the first with `t > 0` and `lo <= -INF`, or `t < 0` and `hi >= INF`): `t > 1e-6`: `hi_k = slack/t`;
/// `t < -1e-6`: `lo_k = slack/t`. Then every column still free whose smaller implied `|bound|` is below `1e6` gets `lb = lo - 0.01` if `|hi| > |lo|`, else
/// `ub = hi + 0.01`.
fn free_column_bounds(
    m: usize,
    n: usize,
    sense: &[u8],
    col_start: &[usize],
    row_idx: &[usize],
    val: &[f64],
    b: &[f64],
    lb: &mut [f64],
    ub: &mut [f64],
) {
    if !(0..n).any(|j| lb[j] <= -INF && ub[j] >= INF) {
        return;
    }
    let mut rows: Vec<Vec<(usize, f64)>> = vec![Vec::new(); m];
    for j in 0..n {
        for q in col_start[j]..col_start[j + 1] {
            rows[row_idx[q]].push((j, val[q]));
        }
    }
    let (mut lo, mut hi) = (lb.to_vec(), ub.to_vec());
    for _sweep in 0..2 {
        for i in 0..m {
            for s in [-1.0f64, 1.0] {
                if s == -1.0 && sense[i] != b'=' {
                    continue;
                }
                let (mut acc, mut ninf) = (0.0f64, 0usize);
                for &(j, a) in &rows[i] {
                    let t = a * s;
                    if t > 0.0 {
                        if lo[j] > -INF {
                            acc = acc + lo[j] * t;
                        } else {
                            ninf += 1;
                        }
                    } else if hi[j] < INF {
                        acc = acc + t * hi[j];
                    } else {
                        ninf += 1;
                    }
                }
                let slack = s * b[i] - acc;
                if -1e-6 > slack {
                    continue;
                }
                if ninf == 0 {
                    for &(j, a) in &rows[i] {
                        let t = s * a;
                        if t > 1e-6 {
                            hi[j] = hi[j].min(slack / t + lo[j]);
                        } else if -1e-6 > t {
                            lo[j] = lo[j].max(slack / t + hi[j]);
                        }
                    }
                } else if ninf == 1 {
                    let mut found = None;
                    for &(j, a) in &rows[i] {
                        let t = s * a;
                        if (t > 0.0 && -INF >= lo[j]) || (0.0 > t && hi[j] >= INF) {
                            found = Some((j, t));
                            break;
                        }
                    }
                    if let Some((j, t)) = found {
                        if t > 1e-6 {
                            hi[j] = slack / t;
                        } else if -1e-6 > t {
                            lo[j] = slack / t;
                        }
                    }
                }
            }
        }
    }
    for j in 0..n {
        if lb[j] <= -INF && ub[j] >= INF && lo[j].abs().min(hi[j].abs()) < 1e6 {
            if hi[j].abs() > lo[j].abs() {
                lb[j] = lo[j] - 0.01;
            } else {
                ub[j] = hi[j] + 0.01;
            }
        }
    }
}

/// Bound classification `0x3283f0` (port of the second part of `classify_bounds` of `tools/barrier_model.py`: the implied upper bounds; the free columns come before, [`free_column_bounds`]). Inputs: the converted model (rows permuted, '>' rows negated, so every non-'=' row is a
/// '<=' row), the right-hand side before the lower-bound shift, bounds with `INF` as infinity. Per row, in order: `mn = sum a>0: a lb + a<=0: a ub`,
/// `mx` the reverse (column order, an infinite bound enters as the number `1e30`), `up = rhs - mn`, `lo = mx - rhs` (`1e100` for non-'=' rows). For every
/// entry whose bounds are both finite, `w = ub - lb`: `a > 1e-4`: the upper bound is dropped (set to `1e100`, which also sets `lo = 1e100`) if
/// `0.9999 a w > up` or the row has no entry with `a <= 0` and `a w > up - 1e-10`; `a < -1e-4`: dropped (sets `up = 1e100`) if `0.9999 a w > lo` or the row
/// has no entry with `a > 0` and `a w > lo - 1e-10`. Later rows see the dropped bounds.
fn drop_implied_upper_bounds(
    m: usize,
    n: usize,
    sense: &[u8],
    col_start: &[usize],
    row_idx: &[usize],
    val: &[f64],
    b: &[f64],
    lb: &[f64],
    ub: &mut [f64],
) {
    let mut rows: Vec<Vec<(usize, f64)>> = vec![Vec::new(); m];
    for j in 0..n {
        for q in col_start[j]..col_start[j + 1] {
            rows[row_idx[q]].push((j, val[q]));
        }
    }
    let big = 1e100;
    for i in 0..m {
        let (mut npos, mut nneg) = (0usize, 0usize);
        let (mut mn, mut mx) = (0.0f64, 0.0f64);
        for &(j, a) in &rows[i] {
            if a > 0.0 {
                npos += 1;
                mn += a * lb[j];
                mx += a * ub[j];
            } else {
                nneg += 1;
                mn += a * ub[j];
                mx += a * lb[j];
            }
        }
        let mut up = b[i] - mn;
        let mut lo = mx - b[i];
        if sense[i] != b'=' {
            lo = big;
        }
        for &(j, a) in &rows[i] {
            let w = ub[j] - lb[j];
            if lb[j] > -INF && ub[j] < INF {
                if a > 1e-4 {
                    if 0.9999 * a * w > up || (nneg == 0 && a * w > up - 1e-10) {
                        ub[j] = big;
                        lo = big;
                    }
                } else if a < -1e-4 && (0.9999 * a * w > lo || (npos == 0 && a * w > lo - 1e-10)) {
                    ub[j] = big;
                    up = big;
                }
            }
        }
    }
}

/// The range columns of the barrier conversion as `(row, upper bound)` in column order: RANGES order, the first entry of a row, bound `R` for an `E` row
/// with `R > 0`, else `|R|` (`-R` for an `E` row with `R <= 0`).
pub fn range_columns(model: &LpModel) -> Vec<(usize, f64)> {
    let mut v: Vec<(usize, f64)> = Vec::new();
    for &(i, r) in &model.range_entries {
        if !v.iter().any(|&(k, _)| k == i) {
            let kind = model.row_kind.get(i).copied().unwrap_or(b'L');
            v.push((
                i,
                if kind == b'E' {
                    if r > 0.0 {
                        r
                    } else {
                        -r
                    }
                } else {
                    r.abs()
                },
            ));
        }
    }
    v
}

/// Conversion of a model to the solved form, with `perm[original row] = position` (the final `SYM` permutation of [`crate::order`]).
pub fn std_form(model: &LpModel, perm: &[i64]) -> Result<StdForm, BarrierError> {
    std_form_scaled(model, perm, None)
}

/// [`std_form`] of the model scaled with `sc` (release profile): `A' = R A C`, `c' = f C c`, `b' = R b`, bounds `/ C`; a range column has coefficient `R_i C_r` and bound `range / C_r`.
pub fn std_form_scaled(
    model: &LpModel,
    perm: &[i64],
    sc: Option<&crate::scale::Scaling>,
) -> Result<StdForm, BarrierError> {
    let orig = model;
    let scaled_model;
    let model = match sc {
        Some(s) => {
            scaled_model = s.apply_loose(model);
            &scaled_model
        }
        None => model,
    };
    let m = model.num_rows();
    let n_struct = model.num_cols();
    // The final permutation covers the rows, then the dense columns as nodes after them.
    let dense =
        crate::order::dense_rule(model).map_err(|e| BarrierError::Unsupported(e.to_string()))?;
    if perm.len() != m + dense.len() {
        return unsupported(format!(
            "permutation of {} entries for {m} rows and {} dense columns",
            perm.len(),
            dense.len()
        ));
    }
    let dense_k: Vec<usize> = (0..dense.len()).map(|r| perm[m + r] as usize - m).collect();
    // A row given twice in RANGES keeps its first entry (the later ones are ignored by the reference); `r_last` is the entry the model's row bounds were built from.
    let mut ranges: Vec<(usize, f64)> = Vec::new();
    let mut r_last = vec![0.0; m];
    for &(i, r) in &model.range_entries {
        if !ranges.iter().any(|&(k, _)| k == i) {
            ranges.push((i, r));
        }
        r_last[i] = r;
    }
    // `Scaling::col` lists the range columns after the structurals in the order of `ranges` (see [`range_columns`]).
    let range_rank: Vec<usize> = (0..ranges.len()).collect();
    let mut row_of = vec![0usize; m];
    for (i, &p) in perm[..m].iter().enumerate() {
        row_of[p as usize] = i;
    }
    // Ranged rows (RANGES order): an equality row plus a range column with upper bound |R| (coefficient +1; -1 for a G row).
    let mut is_range = vec![false; m];
    for &(i, _) in &ranges {
        if model.row_kind.get(i).is_none() {
            return unsupported("ranged row of unknown type");
        }
        is_range[i] = true;
    }
    // Rows: sense, sign and right-hand side in file order.
    let mut sense_file = vec![b'='; m];
    let mut sign = vec![1.0; m];
    // Right-hand side with the sign of the row (a '>' row gives `-lo`).
    let mut b_file = vec![0.0; m];
    for i in 0..m {
        let (lo, up) = (model.row_lo[i], model.row_up[i]);
        if is_range[i] {
            // L: [r - |R|, r], rhs r = row_up; G: [r, r + |R|], rhs r = row_lo; E: [r, r + R] or [r + R, r], rhs row_up.
            b_file[i] = if model.row_kind[i] == b'G' {
                lo
            } else if model.row_kind[i] == b'E' {
                // E: [r, r + R] (R > 0) or [r + R, r] with the last entry R_last; the rhs is r, plus the first entry R when it is positive.
                let r = if r_last[i] > 0.0 { lo } else { up };
                let r1 = ranges.iter().find(|&&(k, _)| k == i).unwrap().1;
                if r1 > 0.0 {
                    r + r1
                } else {
                    r
                }
            } else {
                up
            };
            if is_inf(lo) || is_inf(up) {
                return unsupported("ranged row with an infinite side");
            }
            continue;
        }
        if is_inf(lo) && is_inf(up) {
            return unsupported("free row");
        }
        if !is_inf(lo) && !is_inf(up) {
            if lo != up {
                return unsupported("ranged row without a RANGES entry");
            }
            b_file[i] = lo;
        } else if is_inf(lo) {
            sense_file[i] = b'<';
            b_file[i] = up;
        } else {
            sense_file[i] = b'>';
            sign[i] = -1.0;
            b_file[i] = -lo;
        }
    }
    // Columns: the model's, then the range columns; bounds as the file gives them.
    let nr = ranges.len();
    let n2 = n_struct + nr;
    let mut lb_in = vec![0.0; n2];
    let mut ub_in = vec![INF; n2];
    for j in 0..n_struct {
        let (lo, up) = (model.col_lo[j], model.col_up[j]);
        if !is_inf(lo) && !is_inf(up) && up < lo {
            return unsupported("inconsistent column bounds");
        }
        lb_in[j] = if is_inf(lo) { -INF } else { lo };
        if !is_inf(up) {
            ub_in[j] = up;
        }
    }
    for (k, &(i, r)) in ranges.iter().enumerate() {
        // L and G: |R|; E: R, or -R for R <= 0 (R = 0 gives -0.0).
        let r = match sc {
            // The bound of a range column is the unscaled range over its column factor.
            Some(_) => orig.range_entries.iter().find(|&&(t, _)| t == i).unwrap().1,
            None => r,
        };
        let b = if model.row_kind[i] == b'E' {
            if r > 0.0 {
                r
            } else {
                -r
            }
        } else {
            r.abs()
        };
        ub_in[n_struct + k] = match sc {
            Some(s) => b / s.col[n_struct + range_rank[k]],
            None => b,
        };
    }
    // Matrix in permuted rows, entries ascending by position (range column: its row, coefficient +-1).
    let mut col_start = vec![0usize];
    let mut row_idx = Vec::new();
    let mut val = Vec::new();
    for j in 0..n2 {
        let mut e: Vec<(usize, f64)> = if j < n_struct {
            model
                .col(j)
                .filter(|&(_, v)| v != 0.0)
                .map(|(i, v)| (perm[i] as usize, sign[i] * v))
                .collect()
        } else {
            let i = ranges[j - n_struct].0;
            let u = if model.row_kind[i] == b'G' { -1.0 } else { 1.0 };
            let u = match sc {
                Some(s) => u * s.row[i] * s.col[n_struct + range_rank[j - n_struct]],
                None => u,
            };
            vec![(perm[i] as usize, u)]
        };
        e.sort_by_key(|t| t.0);
        for (i, v) in e {
            row_idx.push(i);
            val.push(v);
        }
        col_start.push(row_idx.len());
    }
    let sense: Vec<u8> = (0..m).map(|k| sense_file[row_of[k]]).collect();
    let mut b = vec![0.0; m];
    for i in 0..m {
        b[perm[i] as usize] = b_file[i];
    }
    // Bound classification: implied bounds of the free columns, then the upper bounds implied by the rows are dropped (the complementarity pairs count the
    // bounds present before).
    let mut lb_cls = lb_in.clone();
    let mut ub_cls = ub_in.clone();
    free_column_bounds(
        m,
        n2,
        &sense,
        &col_start,
        &row_idx,
        &val,
        &b,
        &mut lb_cls,
        &mut ub_cls,
    );
    drop_implied_upper_bounds(
        m,
        n2,
        &sense,
        &col_start,
        &row_idx,
        &val,
        &b,
        &lb_cls,
        &mut ub_cls,
    );
    let mut lower = vec![0.0; n2];
    let mut flags = vec![0u8; n2 + m];
    let mut ub = vec![INF; n2];
    for j in 0..n2 {
        let has_lb = lb_cls[j] > -INF;
        if has_lb {
            lower[j] = lb_cls[j];
        }
        flags[j] = u8::from(has_lb) | if ub_cls[j] < INF { 2 } else { 0 };
        if ub_cls[j] < INF {
            ub[j] = if has_lb {
                ub_cls[j] - lb_cls[j]
            } else {
                ub_cls[j]
            };
        }
    }
    for k in 0..m {
        if sense[k] != b'=' {
            flags[n2 + k] = 1;
        }
    }
    let ncomp = lb_in.iter().filter(|&&l| l > -INF).count()
        + ub_in.iter().filter(|&&u| u < INF).count()
        + sense.iter().filter(|&&s| s != b'=').count();
    let ub_cols: Vec<usize> = (0..n2).filter(|&j| flags[j] & 2 != 0).collect();
    let (ub_lo, ub_hi) = if ub_cols.is_empty() {
        (n2, n2)
    } else {
        (ub_cols[0], ub_cols[ub_cols.len() - 1] + 1)
    };
    // Shift b by the lower bounds: one entry at a time in column order, `b_i = b_i - a_ij lo_j`.
    let mut cshift = 0.0;
    // Maximization: the reference minimises the negated objective (c and the constant negated); the read-out negates the objective, the duals and the reduced costs back.
    let max_sign = if model.sense == Sense::Maximize {
        -1.0
    } else {
        1.0
    };
    let mut c: Vec<f64> = if max_sign < 0.0 {
        model.obj.iter().map(|v| -v).collect()
    } else {
        model.obj.clone()
    };
    c.resize(n2, 0.0);
    for j in 0..n2 {
        let l = lower[j];
        if l != 0.0 {
            for q in col_start[j]..col_start[j + 1] {
                b[row_idx[q]] = b[row_idx[q]] - val[q] * l;
            }
            cshift = cshift + c[j] * l;
        }
    }
    let scale = uniform_scale(n2, &col_start, &val);
    Ok(StdForm {
        b_file,
        m,
        n2,
        n_struct,
        col_start,
        row_idx,
        val,
        sense,
        flags,
        b,
        c,
        ub,
        ub_lo,
        ub_hi,
        offset: if max_sign < 0.0 {
            -model.offset
        } else {
            model.offset
        },
        max_sign,
        cshift,
        lower,
        row_of,
        dense,
        dense_k,
        ncomp,
        scale,
        scaling: sc.cloned(),
        range_rank,
    })
}

// ---------------------------------------------------------------------------------------------------------------------------------
// x87 extended-precision emulation for the compensated dot of `0x35adf0` (64-bit significand, round to nearest even; no overflow,
// no denormals).

#[derive(Clone, Copy, Debug)]
struct Ext {
    neg: bool,
    /// Significand with the top bit set, or 0 for zero.
    m: u64,
    /// Value = m * 2^e.
    e: i32,
}

impl Ext {
    const ZERO: Ext = Ext {
        neg: false,
        m: 0,
        e: 0,
    };

    fn from_f64(d: f64) -> Ext {
        if d == 0.0 {
            return Ext::ZERO;
        }
        let bits = d.to_bits();
        let ex = ((bits >> 52) & 0x7ff) as i32;
        let frac = bits & ((1u64 << 52) - 1);
        let (mant, e) = if ex == 0 {
            (frac, -1074)
        } else {
            (frac | (1u64 << 52), ex - 1075)
        };
        let lz = mant.leading_zeros();
        Ext {
            neg: d < 0.0,
            m: mant << lz,
            e: e - lz as i32,
        }
    }

    /// Rounds `m * 2^e` (plus a sticky tail below the lowest bit when `sticky`) to 64 significant bits.
    fn round(neg: bool, m: u128, e: i32, sticky: bool) -> Ext {
        if m == 0 {
            return Ext::ZERO;
        }
        let bits = 128 - m.leading_zeros() as i32;
        if bits <= 64 {
            let sh = 64 - bits;
            return Ext {
                neg,
                m: (m << sh) as u64,
                e: e - sh,
            };
        }
        let sh = (bits - 64) as u32;
        let rem = m & ((1u128 << sh) - 1);
        let half = 1u128 << (sh - 1);
        let mut q = (m >> sh) as u64;
        let mut e = e + sh as i32;
        if rem > half || (rem == half && (sticky || q & 1 == 1)) {
            let (nq, over) = q.overflowing_add(1);
            q = nq;
            if over {
                q = 1u64 << 63;
                e += 1;
            }
        }
        Ext { neg, m: q, e }
    }

    fn mul_f64(a: f64, b: f64) -> Ext {
        let (x, y) = (Ext::from_f64(a), Ext::from_f64(b));
        if x.m == 0 || y.m == 0 {
            return Ext::ZERO;
        }
        Ext::round(x.neg != y.neg, x.m as u128 * y.m as u128, x.e + y.e, false)
    }

    fn neg(self) -> Ext {
        Ext {
            neg: !self.neg,
            ..self
        }
    }

    fn add(self, other: Ext) -> Ext {
        if self.m == 0 {
            return other;
        }
        if other.m == 0 {
            return self;
        }
        // a = the larger magnitude (significands are normalized, so the exponent decides, then the significand).
        let (a, b) = if (self.e, self.m) >= (other.e, other.m) {
            (self, other)
        } else {
            (other, self)
        };
        let shift = (a.e - b.e) as u32;
        let aa = (a.m as u128) << 62;
        let (mut bb, mut sticky) = (0u128, true);
        if shift < 126 {
            let full = (b.m as u128) << 62;
            bb = full >> shift;
            sticky = shift > 0 && (full & ((1u128 << shift) - 1)) != 0;
        }
        if sticky {
            bb |= 1;
        }
        let e = a.e - 62;
        if a.neg == b.neg {
            Ext::round(a.neg, aa + bb, e, false)
        } else {
            Ext::round(a.neg, aa - bb, e, false)
        }
    }

    fn sub(self, other: Ext) -> Ext {
        self.add(other.neg())
    }

    fn to_f64(self) -> f64 {
        if self.m == 0 {
            return 0.0;
        }
        let rem = self.m & 0x7ff;
        let mut q = self.m >> 11;
        let mut e = self.e + 11;
        if rem > 0x400 || (rem == 0x400 && q & 1 == 1) {
            q += 1;
            if q == 1u64 << 53 {
                q >>= 1;
                e += 1;
            }
        }
        let v = ldexp(q as f64, e);
        if self.neg {
            -v
        } else {
            v
        }
    }
}

fn ldexp(v: f64, e: i32) -> f64 {
    let h = e / 2;
    v * 2f64.powi(h) * 2f64.powi(e - h)
}

/// `0x35adf0`: Kahan summation in extended precision of the extended products `a_i b_i`:
/// `t = a b; y = t - c; s' = s + y; c = (s' - s) - y` with every operation rounded to 64 bits (C).
fn kahan_dot(a: &[f64], b: &[f64]) -> Ext {
    let (mut s, mut c) = (Ext::ZERO, Ext::ZERO);
    for i in 0..a.len() {
        let t = Ext::mul_f64(a[i], b[i]);
        let y = t.sub(c);
        let sn = s.add(y);
        c = sn.sub(s).sub(y);
        s = sn;
    }
    s
}

/// `Q+0x1f0` (mu) of `0x2ecde0` (C): `r2 = x.z` over all entries, `r3 = u.v` over `[lo, hi)`, `s = double(r3 + double(r2 + 0))`,
/// `mu = (tau kappa + s) / ncomp`.
pub fn complementarity(
    x: &[f64],
    z: &[f64],
    u: &[f64],
    v: &[f64],
    lo: usize,
    hi: usize,
    tau: f64,
    kappa: f64,
    ncomp: usize,
) -> f64 {
    let r2 = kahan_dot(x, z);
    let r3 = kahan_dot(&u[lo..hi], &v[lo..hi]);
    let c = r2.add(Ext::from_f64(0.0)).to_f64();
    let d = r3.add(Ext::from_f64(c)).to_f64();
    (tau * kappa + d) / ncomp as f64
}

// ---------------------------------------------------------------------------------------------------------------------------------
// Double kernels.

/// `0x311b80` (C): `sum a_i b_i`. `n < 8` (or `n < peel + 8`, `peel = parity & 1`): sequential; otherwise `peel` scalar products, blocks of 8 into
/// four 2-lane accumulators (acc0 starts at the running sum), `(acc0 + acc1) + (acc2 + acc3)` per lane, lane 0 + lane 1, then the scalar tail.
pub fn dot_lanes(a: &[f64], b: &[f64], parity: usize) -> f64 {
    let n = a.len();
    if n == 0 {
        return 0.0;
    }
    let peel = parity & 1;
    if n >= 8 && n >= peel + 8 {
        let mut s = 0.0;
        for i in 0..peel {
            s = s + a[i] * b[i];
        }
        let mut acc = [[s, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0]];
        let end = n - ((n - peel) & 7);
        let mut k = peel;
        while k < end {
            for g in 0..4 {
                for l in 0..2 {
                    acc[g][l] = acc[g][l] + a[k + 2 * g + l] * b[k + 2 * g + l];
                }
            }
            k += 8;
        }
        let t0 = [acc[0][0] + acc[1][0], acc[0][1] + acc[1][1]];
        let t1 = [acc[2][0] + acc[3][0], acc[2][1] + acc[3][1]];
        let t = [t0[0] + t1[0], t0[1] + t1[1]];
        let mut s = t[0] + t[1];
        for i in end..n {
            s = s + a[i] * b[i];
        }
        return s;
    }
    let mut s = 0.0;
    for i in 0..n {
        s = s + a[i] * b[i];
    }
    s
}

/// `0x311d00`: `max |a_i|`.
pub fn max_abs(a: &[f64]) -> f64 {
    let mut m = 0.0;
    for &v in a {
        if v.abs() >= m {
            m = v.abs();
        }
    }
    m
}

/// The scalars of `0x2e4f10` (C): norms of `b` (rows in file order) and `c` (columns in file order) and the feasibility thresholds.
#[derive(Debug, Clone, Copy)]
pub struct Setup {
    pub q160: f64,
    pub q168: f64,
    pub q170: f64,
    pub q178: f64,
    pub q260: f64,
    pub q268: f64,
    pub q270: f64,
}

pub fn setup_scalars(b_file: &[f64], c_file: &[f64], tol: f64) -> Setup {
    let q160 = dot_lanes(b_file, b_file, 0).sqrt();
    let q170 = dot_lanes(c_file, c_file, 0).sqrt();
    Setup {
        q160,
        q168: max_abs(b_file),
        q170,
        q178: max_abs(c_file),
        q260: ((b_file.len() as f64 + 1.0).sqrt() + q160) * tol,
        q268: ((c_file.len() as f64 + 1.0).sqrt() + q170) * tol,
        q270: 10.0 * tol,
    }
}

/// Sum of `t(base + k)`, `k < count`, as `0x30c670` adds a block whose 16-byte alignment is the parity of `base`: fewer than 8 (or `8 + peel`)
/// elements one by one onto `a0`; otherwise `peel` scalars first, blocks of 8 into four lane accumulators (lane `l` gets `t[k+l]` then
/// `t[k+l+4]`, lane 0 starts at `a0`), `(a0 + a2) + (a1 + a3)`, then the scalar tail (C).
fn lane_block(t: &dyn Fn(usize) -> f64, base: usize, count: usize, a0: f64) -> f64 {
    let peel = base & 1;
    if count >= 8 && count >= peel + 8 {
        let mut a = [a0, 0.0, 0.0, 0.0];
        for i in 0..peel {
            a[0] = a[0] + t(base + i);
        }
        let end = count - ((count - peel) & 7);
        let mut k = peel;
        while k < end {
            for l in 0..4 {
                a[l] = a[l] + t(base + k + l) + t(base + k + l + 4);
            }
            k += 8;
        }
        let mut tot = (a[0] + a[2]) + (a[1] + a[3]);
        for i in end..count {
            tot = tot + t(base + i);
        }
        return tot;
    }
    let mut tot = a0;
    for i in 0..count {
        tot = tot + t(base + i);
    }
    tot
}

/// The current iterate (`Q+0x658` x, `0x668` z, `0x660` y, `0x670` u, `0x678` v). x, z have `n2 + m` entries, u, v `n2`, y `m`.
#[derive(Debug, Clone)]
pub struct State {
    pub x: Vec<f64>,
    pub z: Vec<f64>,
    pub y: Vec<f64>,
    pub u: Vec<f64>,
    pub v: Vec<f64>,
}

/// A search direction: `dx`, `dz` over all variables, `dy` over the rows, `du`, `dv` over the columns.
#[derive(Debug, Clone, Default)]
pub struct Direction {
    pub dx: Vec<f64>,
    pub dy: Vec<f64>,
    pub dz: Vec<f64>,
    pub du: Vec<f64>,
    pub dv: Vec<f64>,
}

/// `0x30c670` (C): mean complementarity after the trial step.
fn mu_aff(
    f: &StdForm,
    s: &State,
    d: &Direction,
    ap: f64,
    ad: f64,
    ncomp: usize,
    hk: Option<(f64, f64, f64, f64)>,
) -> f64 {
    let t = |i: usize| (d.dx[i] * ap + s.x[i]) * (d.dz[i] * ad + s.z[i]);
    let mut acc = 0.0;
    for j in 0..f.n2 {
        if f.flags[j] & 1 != 0 {
            acc = acc + t(j);
        }
    }
    let mut tot = lane_block(&t, f.n2, f.m, acc);
    if f.has_ub() {
        let tu = |j: usize| (d.du[j] * ap + s.u[j]) * (d.dv[j] * ad + s.v[j]);
        tot = lane_block(&tu, f.ub_lo, f.ub_hi - f.ub_lo, tot);
    }
    if let Some((dtau, dkap, tau, kappa)) = hk {
        // Homogeneous term: `(dkappa ad + kappa)(dtau ap + tau)`.
        tot = (dkap * ad + kappa) * (dtau * ap + tau) + tot;
    }
    tot / ncomp as f64
}

/// `0x2e0410` after `0x30c670` (C): `r = aff / max(mu, 1e-30)`, `sigma = r r min(r, 0.1)`, `target = mu sigma`.
fn sigma_target(mu: f64, aff: f64) -> f64 {
    let r = aff / mu.max(1e-30);
    mu * (r * r * r.min(0.1))
}

/// `0x30c050` / `0x30ed50` (C): step to the boundary of `x + a d >= 0` over the flagged entries (and `u + a e` over `[lo, hi)`).
fn ratio(
    x: &[f64],
    d: &dyn Fn(usize) -> f64,
    fl: &[u8],
    up: Option<(usize, usize, &[f64], &dyn Fn(usize) -> f64)>,
    tk: &[(f64, f64, f64)],
    w: f64,
) -> f64 {
    let mut a = 1.001;
    for j in 0..x.len() {
        if fl[j] & 1 != 0 {
            let dj = d(j);
            if a * dj + x[j] < 0.0 {
                a = -(x[j] + (-1e-30)) / dj - 1e-13;
            }
        }
    }
    if let Some((lo, hi, u, e)) = up {
        for j in lo..hi {
            let dj = e(j);
            if a * dj + u[j] < 0.0 {
                a = -(u[j] + (-1e-30)) / dj - 1e-13;
            }
        }
    }
    // The homogeneous scalars tau and kappa (value, direction, corrector direction) are tested after the columns, like further entries.
    for &(t, e0, e1) in tk {
        let dj = e1 * w + e0;
        if a * dj + t < 0.0 {
            a = -(t + (-1e-30)) / dj - 1e-13;
        }
    }
    if a <= 0.0 {
        0.0
    } else {
        a
    }
}

fn step_pair(
    f: &StdForm,
    s: &State,
    d: &Direction,
    d2: Option<(&Direction, f64)>,
    primal: bool,
    tk: &[(f64, f64, f64)],
) -> f64 {
    let w = d2.map_or(0.0, |(_, w)| w);
    // Direction d, or d + w d2 as `d2 w + d`.
    if primal {
        let dir = |j: usize| match d2 {
            Some((e, w)) => e.dx[j] * w + d.dx[j],
            None => d.dx[j],
        };
        let du = |j: usize| match d2 {
            Some((e, w)) => e.du[j] * w + d.du[j],
            None => d.du[j],
        };
        let up: Option<(usize, usize, &[f64], &dyn Fn(usize) -> f64)> = if f.has_ub() {
            Some((f.ub_lo, f.ub_hi, &s.u, &du))
        } else {
            None
        };
        ratio(&s.x, &dir, &f.flags, up, tk, w)
    } else {
        let dir = |j: usize| match d2 {
            Some((e, w)) => e.dz[j] * w + d.dz[j],
            None => d.dz[j],
        };
        let dv = |j: usize| match d2 {
            Some((e, w)) => e.dv[j] * w + d.dv[j],
            None => d.dv[j],
        };
        let up: Option<(usize, usize, &[f64], &dyn Fn(usize) -> f64)> = if f.has_ub() {
            Some((f.ub_lo, f.ub_hi, &s.v, &dv))
        } else {
            None
        };
        ratio(&s.z, &dir, &f.flags, up, tk, w)
    }
}

/// `0x30ed50` pass 1 (C): weights `w_i = ((ntry-1-i)(1-wlo))/(ntry-1) + wlo`, `wlo = clamp(ap0 ad0, 0.01, 1)`; a weight replaces the best
/// step only when `best * 1.0001 < step` (primal and dual separately). Returns `(best_p, best_d, wp, wd)`.
fn select(
    f: &StdForm,
    s: &State,
    d: &Direction,
    d2: &Direction,
    ap0: f64,
    ad0: f64,
    ntry: usize,
    tk: &[(f64, f64, f64)],
    k: usize,
) -> (f64, f64, f64, f64) {
    // `0x30ed95`: for a corrector k > 1 the smallest weight is 0 (the last trial keeps the steps so far).
    let wlo = if k > 1 {
        0.0
    } else {
        (ap0 * ad0).min(1.0).max(0.01)
    };
    let (mut bp, mut bd) = (-1e100, -1e100);
    let (mut wp, mut wd) = (0.0, 0.0);
    for i in 0..ntry {
        let w = if ntry > 1 {
            (((ntry - 1 - i) as f64) * (1.0 - wlo)) / ((ntry - 1) as f64) + wlo
        } else {
            wlo
        };
        let (ap, ad) = if w != 0.0 {
            (
                step_pair(f, s, d, Some((d2, w)), true, tk),
                step_pair(f, s, d, Some((d2, w)), false, tk),
            )
        } else {
            (ap0, ad0)
        };
        if bp * 1.0001 < ap {
            bp = ap;
            wp = w;
        }
        if bd * 1.0001 < ad {
            bd = ad;
            wd = w;
        }
    }
    (bp, bd, wp, wd)
}

fn combine(d0: &[f64], d1: &[f64], w: f64) -> Vec<f64> {
    (0..d0.len()).map(|j| d1[j] * w + d0[j]).collect()
}

/// The homogeneous scalars of an update: `(tau, kappa, dtau, dkappa, dtau of the dual step)`.
type HsdStep = (f64, f64, f64, f64, f64);

/// New homogeneous scalars after an update: `tau` after the primal step, `kappa`, and `tau` after the dual step.
#[derive(Debug, Clone, Copy)]
struct HsdOut {
    tau: f64,
    kappa: f64,
    tau_d: f64,
}

/// `0x30d3e0` (C): step-back gamma and the update of x, z, y (u, v). `g0 = 1e-6` (no cones).
/// Returns the new state, `(g, sp, sd)` and, with the homogeneous scalars `hk`, their new values.
fn update(
    f: &StdForm,
    s: &State,
    d: &Direction,
    ap: f64,
    ad: f64,
    ncomp: usize,
    hk: Option<HsdStep>,
) -> (State, (f64, f64, f64), Option<HsdOut>) {
    let g0 = 1e-6;
    let steps = |g: f64| ((ap * (1.0 - g)).min(1.0), ((1.0 - g) * ad).min(1.0));
    let (sp, sd) = steps(g0);
    let (mut tot, mut minp) = (0.0, 1e100);
    for j in 0..f.n() {
        if f.flags[j] & 1 != 0 {
            let p = (d.dx[j] * sp + s.x[j]) * (d.dz[j] * sd + s.z[j]);
            minp = f64::min(minp, p);
            tot = tot + p;
        }
    }
    for j in f.ub_lo..f.ub_hi {
        if f.flags[j] & 2 != 0 {
            let p = (d.du[j] * sp + s.u[j]) * (d.dv[j] * sd + s.v[j]);
            minp = f64::min(minp, p);
            tot = tot + p;
        }
    }
    // Trial term of the homogeneous scalars `(dkappa sd + kappa)(dtau sp + tau)` (`0 * 1` without them).
    let (tau0, kap0, dtau, dkap) = hk.map_or((1.0, 0.0, 0.0, 0.0), |h| (h.0, h.1, h.2, h.3));
    tot = tot + (dkap * sd + kap0) * (dtau * sp + tau0);
    tot = tot / ncomp as f64;
    let mut g = g0;
    if minp < 0.01 * tot {
        let gg = (0.01 * g0 * tot) / minp;
        if g0 <= gg {
            g = gg;
        }
    }
    g = g.min(0.5);
    let (mut sp, mut sd) = steps(g);
    if kap0 > 0.0 {
        // `0x30dc2d`: while the trial `tau kappa` product is below 0.01 of the mean, g grows tenfold (at most 0.5).
        while g < 0.5 && 0.01 * tot > (dkap * sd + kap0) * (dtau * sp + tau0) {
            g = (g * 10.0).min(0.5);
            (sp, sd) = steps(g);
        }
    }
    let hout = hk.map(|h| HsdOut {
        tau: dtau * sp + tau0,
        kappa: dkap * sp + kap0,
        tau_d: h.4 * sd + tau0,
    });
    let n = f.n();
    let (mut x, mut z) = (vec![0.0; n], vec![0.0; n]);
    for j in 0..n {
        x[j] = sp * d.dx[j] + s.x[j];
        z[j] = sd * d.dz[j] + s.z[j];
        if f.flags[j] & 1 != 0 {
            x[j] = x[j].max(1e-30);
            z[j] = z[j].max(1e-30);
        }
    }
    let y: Vec<f64> = (0..f.m).map(|i| d.dy[i] * sd + s.y[i]).collect();
    let (mut u, mut v) = (s.u.clone(), s.v.clone());
    for j in f.ub_lo..f.ub_hi {
        if f.flags[j] & 2 != 0 {
            u[j] = (sp * d.du[j] + s.u[j]).max(1e-30);
            v[j] = (sd * d.dv[j] + s.v[j]).max(1e-30);
        }
    }
    (State { x, z, y, u, v }, (g, sp, sd), hout)
}

// ---------------------------------------------------------------------------------------------------------------------------------
// Residuals and scalars of `0x2ecde0` (C).

/// Residual vectors and scalars of one evaluation; scalar names are the reference's `Q` offsets.
#[derive(Debug, Clone)]
pub struct Residuals {
    /// `tau b - Af x` (`Q+0x7a0`, m entries).
    pub rp: Vec<f64>,
    /// Dual residual (`Q+0x7a8`, n entries).
    pub rd: Vec<f64>,
    /// `ub tau - x - u` over `[lo, hi)` (zero elsewhere; n2 entries).
    pub ru: Vec<f64>,
    pub q130: f64,
    pub q138: f64,
    pub q140: f64,
    pub q148: f64,
    pub q198: f64,
    pub q1a0: f64,
    pub q1a8: f64,
    pub q1b0: f64,
    pub q1c0: f64,
    pub q1d8: f64,
    pub q1e0: f64,
    pub q1f0: f64,
    pub q1f8: f64,
}

fn dual_residual(f: &StdForm, s: &State, tau: f64) -> Vec<f64> {
    let (n2, m) = (f.n2, f.m);
    let mut rd = vec![0.0; n2 + m];
    let end = if n2 >= 8 { n2 - (n2 & 7) } else { 0 };
    for j in 0..n2 {
        let t = f.aty_col(j, &s.y);
        let zj = if f.flags[j] & 1 != 0 { s.z[j] } else { 0.0 };
        let vj = if f.flags[j] & 2 != 0 { s.v[j] } else { 0.0 };
        rd[j] = if j < end {
            tau * f.c[j] - ((zj + t) - vj)
        } else {
            ((tau * f.c[j] - t) - zj) + vj
        };
    }
    // Slack part of `0x2ecde0` (C, DEGEN2 and the five models): four loops with different zero signs for '=' rows. m > 6, m >= 8: `peel` (parity of
    // n2) elements `(-z) - t`, blocks of 8 `-(z + t)`, pairs `0 - (t + z)`, a last single `(0 - z) - t`; m = 7 or m + peel < 8: pairs and tail only;
    // m <= 6: pairs and tail of `(-z) - t`.
    let ts: Vec<f64> = (0..m)
        .map(|i| if f.sense[i] != b'=' { s.y[i] } else { 0.0 })
        .collect();
    let zs = &s.z[n2..n2 + m];
    if m > 6 {
        let mut peel = if m >= 8 { n2 & 1 } else { 0 };
        if m >= 8 && m < peel + 8 {
            peel = 0;
        }
        let end = if m >= 8 && m >= peel + 8 {
            m - ((m - peel) & 7)
        } else {
            0
        };
        for i in 0..peel {
            rd[n2 + i] = (-zs[i]) - ts[i];
        }
        for i in (if end != 0 { peel } else { 0 })..end {
            rd[n2 + i] = -(zs[i] + ts[i]);
        }
        let mut k = end;
        while m - k >= 2 {
            rd[n2 + k] = 0.0 - (ts[k] + zs[k]);
            rd[n2 + k + 1] = 0.0 - (ts[k + 1] + zs[k + 1]);
            k += 2;
        }
        for i in k..m {
            rd[n2 + i] = (0.0 - zs[i]) - ts[i];
        }
    } else {
        let mut k = 0;
        while m - k >= 2 {
            rd[n2 + k] = (-zs[k]) - ts[k];
            rd[n2 + k + 1] = (-zs[k + 1]) - ts[k + 1];
            k += 2;
        }
        for i in k..m {
            rd[n2 + i] = (-zs[i]) - ts[i];
        }
    }
    rd
}

/// `scale` is the objective scale `Q+0x158` (1 for the controlled profile).
pub fn residuals(
    f: &StdForm,
    s: &State,
    tau: f64,
    kappa: f64,
    ncomp: usize,
    scale: f64,
) -> Residuals {
    let (m, n2) = (f.m, f.n2);
    let (lo, hi) = (f.ub_lo, f.ub_hi);
    let ax = f.af(&s.x);
    let rp: Vec<f64> = (0..m).map(|i| tau * f.b[i] - ax[i]).collect();
    let rd = dual_residual(f, s, tau);
    let mut ru = vec![0.0; n2];
    for j in lo..hi {
        if f.flags[j] & 2 != 0 {
            ru[j] = (f.ub[j] * tau - s.x[j]) - s.u[j];
        }
    }
    let q1a0 = dot_lanes(&rp, &rp, 0).sqrt();
    let q1c0 = if hi > lo {
        dot_lanes(&ru[lo..hi], &ru[lo..hi], lo).sqrt()
    } else {
        0.0
    };
    let q1b0 = dot_lanes(&rd, &rd, 0).sqrt();
    let mut mx = 0.0;
    for &t in rp
        .iter()
        .chain((lo..hi).filter(|&j| f.flags[j] & 2 != 0).map(|j| &ru[j]))
    {
        if mx <= t.abs() {
            mx = t.abs();
        }
    }
    let q1d8 = mx / tau;
    let mut dm = 0.0;
    for j in 0..n2 {
        let mut d = rd[j];
        if f.flags[j] & 1 == 0 || 0.0 <= d {
            d = d.abs();
        } else {
            d = -(d + s.z[j]);
            if d <= 0.0 {
                d = 0.0;
            }
        }
        if dm <= d {
            dm = d;
        }
    }
    for i in 0..m {
        if f.sense[i] != b'=' && dm <= s.y[i] {
            dm = s.y[i];
        }
    }
    let q138 = dot_lanes(&s.x[..n2], &f.c, 0) + 0.0 + f.cshift * tau;
    let q130 = (q138 / scale) / tau + f.offset;
    let g = dot_lanes(&s.y, &f.b, 0)
        - if hi > lo {
            dot_lanes(&s.v[lo..hi], &f.ub[lo..hi], lo)
        } else {
            0.0
        };
    let q148 = ((0.0 + g) - 0.0) + tau * f.cshift;
    let q140 = (q148 / scale) / tau + f.offset;
    let q1f0 = complementarity(&s.x, &s.z, &s.u, &s.v, lo, hi, tau, kappa, ncomp);
    Residuals {
        rp,
        rd,
        ru,
        q130,
        q138,
        q140,
        q148,
        q198: q1a0 / tau,
        q1a0,
        q1a8: q1b0 / tau,
        q1b0,
        q1c0,
        q1d8,
        q1e0: dm / tau,
        q1f0,
        q1f8: (q1f0 / tau) / tau,
    }
}

// ---------------------------------------------------------------------------------------------------------------------------------
// Starting point `0x2e8ef0` (C, five models).

/// The diagonal `D` of the starting-point factor: 1 for a column with a lower bound only (or an upper bound only), 0.5 with both bounds, `1e-6` for a free
/// column, 0 for the slack of an '=' row, 1 for the other slacks.
pub fn start_diag(f: &StdForm) -> Vec<f64> {
    (0..f.n())
        .map(|j| match (j < f.n2, f.flags[j]) {
            (true, 3) => 0.5,
            (true, 0) => 1e-6,
            (true, _) => 1.0,
            (false, 0) => 0.0,
            (false, _) => 1.0,
        })
        .collect()
}

/// `zx = Q+0x800` and the diagonal `D` of `A D A'` at an iterate (`0x2ecde0`, C): `zx = z/x (+ v/u for a boxed column)`, `D = min(1 / zx, 1e12)`. A column
/// with `zx <= 0` (a free column) gets `zx = 0` and `D = 1 / min(max((|pobj| + 1) 1e-8 / max(|x|, 1e-10)^2, 1e-12), 1e-6)`; the slack of an '=' row 0 and 0.
pub fn normal_diag(f: &StdForm, s: &State, pobj: f64) -> (Vec<f64>, Vec<f64>) {
    let n = f.n();
    let (mut zx, mut d) = (vec![0.0; n], vec![0.0; n]);
    let reg = (pobj.abs() + 1.0) * 1e-8;
    for j in 0..n {
        let fl = f.flags[j];
        if j < f.n2 {
            let mut t = 0.0;
            if fl & 1 != 0 {
                t = s.z[j] / s.x[j] + 0.0;
            }
            if fl & 2 != 0 {
                t = t + s.v[j] / s.u[j];
            }
            if t > 0.0 {
                zx[j] = t;
                d[j] = (1.0 / t).min(1e12);
            } else {
                let a = s.x[j].abs().max(1e-10);
                d[j] = 1.0 / (reg / (a * a)).max(1e-12).min(1e-6);
            }
        } else if fl != 0 {
            zx[j] = s.z[j] / s.x[j];
            d[j] = (1.0 / zx[j]).min(1e12);
        }
    }
    (zx, d)
}

/// `0x2eb560` (C): the shifts `(p3, p4)` of the starting point from the minima and sums over the flagged entries.
fn start_scalars(
    f: &StdForm,
    x: &[f64],
    z: &[f64],
    u: &[f64],
    v: &[f64],
    q168: f64,
    q178: f64,
    q30: f64,
) -> (f64, f64) {
    let (n2, n) = (f.n2, f.n());
    let fl = &f.flags;
    let (mut mx, mut mz) = (1e100, 1e100);
    for j in 0..n {
        if fl[j] & 1 != 0 {
            if x[j] <= mx {
                mx = x[j];
            }
            if z[j] <= mz {
                mz = z[j];
            }
        }
    }
    for j in 0..n2 {
        if fl[j] & 2 != 0 {
            if u[j] <= mx {
                mx = u[j];
            }
            if v[j] <= mz {
                mz = v[j];
            }
        }
    }
    let mut dx0 = mx * -1.5;
    if dx0 <= 1e-4 {
        dx0 = 1e-4;
    }
    let mut dz0 = -1.5 * mz;
    if dz0 <= 1e-4 {
        dz0 = 1e-4;
    }
    let (mut sx, mut sz, mut sxz) = (0.0, 0.0, 0.0);
    let mut acc = |a: f64, b: f64| {
        sx = sx + a;
        sz = sz + b;
        sxz = sxz + a * b;
    };
    for j in 0..n2 {
        if fl[j] & 1 != 0 {
            acc(x[j] + dx0, z[j] + dz0);
        }
        if fl[j] & 2 != 0 {
            acc(u[j] + dx0, v[j] + dz0);
        }
    }
    for j in n2..n {
        if fl[j] & 1 != 0 {
            acc(x[j] + dx0, z[j] + dz0);
        }
    }
    let (mut p3, mut p4) = (100.0, 100.0);
    if 0.0 < sz {
        p3 = (sxz * 0.5) / sz + dx0;
    }
    if 0.0 < sx {
        p4 = dz0 + (sxz * 0.5) / sx;
    }
    let t = (q168 * 1.0001 * p3).sqrt();
    if p3 <= t {
        p3 = t;
    }
    let t = (q178 * 1.0001 * p4).sqrt();
    if p4 <= t {
        p4 = t;
    }
    // Floor Q+0x30 (0 -> 0.1 for an LP; 1000 at the restart).
    let floor = if q30 == 0.0 { 0.1 } else { q30 };
    if p3 <= floor {
        p3 = floor;
    }
    if p4 <= floor {
        p4 = floor;
    }
    (p3, p4)
}

fn start_dual_stage1(f: &StdForm) -> (Vec<f64>, Vec<f64>) {
    let mut z = vec![0.0; f.n()];
    let mut v = vec![0.0; f.n2];
    for j in 0..f.n2 {
        match f.flags[j] {
            0 => {}
            3 => {
                let t = f.c[j] * 0.5;
                z[j] = t;
                v[j] = -t;
            }
            fl if fl & 1 == 0 => v[j] = -f.c[j],
            _ => z[j] = f.c[j],
        }
    }
    (z, v)
}

fn start_dual_stage2(f: &StdForm, r: &[f64]) -> (Vec<f64>, Vec<f64>) {
    let mut z = vec![0.0; f.n()];
    let mut v = vec![0.0; f.n2];
    for j in 0..f.n2 {
        let fl = f.flags[j];
        if fl == 3 {
            z[j] = -r[j];
            v[j] = r[j];
        } else if fl & 1 == 0 {
            if fl & 2 != 0 {
                v[j] = r[j];
            }
        } else {
            z[j] = -r[j];
        }
    }
    for i in 0..f.m {
        z[f.n2 + i] = -r[f.n2 + i];
    }
    (z, v)
}

/// The starting point from the two solves of `0x305f10` (stage 1: `x~`, `y~`; stage 2: `r`, `y`), after the clamps and before the z reduction.
/// `mode` is the parameter block word `+0x2bf0` (0 by default). Returns `(state, p3, p4, stage)`.
pub fn start_point(
    f: &StdForm,
    s1: (Vec<f64>, Vec<f64>),
    s2: (Vec<f64>, Vec<f64>),
    sc: &Setup,
    mode: i32,
    q30: f64,
) -> (State, f64, f64, u8) {
    let (n2, n) = (f.n2, f.n());
    let mut x = s1.0;
    let mut u = vec![0.0; n2];
    for j in 0..n2 {
        if f.flags[j] & 2 != 0 {
            u[j] = f.ub[j] - x[j];
        }
    }
    let (mut z, mut v) = start_dual_stage1(f);
    let (mut p3, mut p4) = start_scalars(f, &x, &z, &u, &v, sc.q168, sc.q178, q30);
    let (mut stage, mut y) = (1u8, vec![0.0; f.m]);
    if mode != 1 {
        let (z2, v2) = start_dual_stage2(f, &s2.0);
        let (q3, q4) = start_scalars(f, &x, &z2, &u, &v2, sc.q168, sc.q178, q30);
        if mode == 2 || q4 < p4 {
            z = z2;
            v = v2;
            p3 = q3;
            p4 = q4;
            stage = 2;
            y = s2.1[..f.m].to_vec();
        }
    }
    for j in 0..n {
        if f.flags[j] & 1 != 0 {
            x[j] = (p3 + x[j]).min(1e10);
            z[j] = (z[j] + p4).min(1e10);
        } else {
            x[j] = x[j].max(-1e10).min(1e10);
        }
    }
    for j in 0..n2 {
        if f.flags[j] & 2 != 0 {
            u[j] = (p3 + u[j]).min(1e10);
            v[j] = (v[j] + p4).min(1e10);
        }
    }
    let y = y.into_iter().map(|t| t.max(-1e10).min(1e10)).collect();
    (State { x, z, y, u, v }, p3, p4, stage)
}

/// `0x2f1560` (C): pulls z (and v) towards `mu / x` within the dual residual: for every flagged entry `d = mu/x - z`,
/// `d <= 0: d = max(d, min(-1e-20, r))`, else `d = min(d, max(1e-20, r))`; `z += d; r -= d`. Slacks only when the flag is exactly 1.
pub fn z_reduce(f: &StdForm, s: &mut State, r: &mut [f64], mu: f64) {
    let clamp = |d: f64, rr: f64| {
        if d <= 0.0 {
            d.max((-1e-20f64).min(rr))
        } else {
            d.min(1e-20f64.max(rr))
        }
    };
    for j in 0..f.n2 {
        if f.flags[j] & 1 != 0 {
            let d = clamp(mu / s.x[j] - s.z[j], r[j]);
            s.z[j] = s.z[j] + d;
            r[j] = r[j] - d;
        }
        if f.flags[j] & 2 != 0 {
            let d = clamp(mu / s.u[j] - s.v[j], -r[j]);
            s.v[j] = s.v[j] + d;
            r[j] = r[j] + d;
        }
    }
    for j in f.n2..f.n() {
        if f.flags[j] == 1 {
            let d = clamp(mu / s.x[j] - s.z[j], r[j]);
            s.z[j] = s.z[j] + d;
            r[j] = r[j] - d;
        }
    }
}

// ---------------------------------------------------------------------------------------------------------------------------------
// Termination `0x2e2fe0` (C for the continue and optimal returns of the corpus, H for the others).

/// Inputs of the status test: the reference's `Q` and parameter-block words under their offsets.
#[derive(Debug, Clone)]
pub struct TermIn {
    pub it: i32,
    pub q130: f64,
    pub q138: f64,
    pub q140: f64,
    pub q148: f64,
    pub q128: f64,
    pub q1e8: f64,
    pub off: f64,
    pub q1f0: f64,
    pub q1f8: f64,
    pub tol: f64,
    pub ncomp: f64,
    pub tau: f64,
    pub kappa: f64,
    pub q198: f64,
    pub q1a8: f64,
    pub q1d8: f64,
    pub q1e0: f64,
    pub q1a0: f64,
    pub q1b0: f64,
    pub q1c0: f64,
    pub q1d0: f64,
    pub q220: f64,
    pub q228: f64,
    pub q160: f64,
    pub q170: f64,
    pub q260: f64,
    pub q268: f64,
    pub q270: f64,
    pub q278: f64,
    pub q2a0: i32,
    pub q5e0: f64,
    pub q5e8: f64,
    /// `Q+0x88`, `0xf8`, `0xa8`, `0xf4`: iteration of the last improvement of the progress measures (written by `0x2e39b0`, undecoded).
    pub i88: i32,
    pub if8: i32,
    pub ia8: i32,
    pub if4: i32,
    /// `P+0x2dd4` (stall flag), `P+0x2460`, `P+0x2bbc` (iteration limit), `P+0x2c10`; time limit exceeded.
    pub p2dd4: bool,
    pub p2460: f64,
    pub p2bbc: i32,
    pub p2c10: i32,
    pub time_up: bool,
    /// Histories by iteration: `0x50` (Q198), `0x68` (Q1a8), `0x78` (Q1f8), `0x80` (Q1f0), each including the current iteration.
    pub h50: Vec<f64>,
    pub h68: Vec<f64>,
    pub h78: Vec<f64>,
    pub h80: Vec<f64>,
}

/// `0x2e2fe0`: 1 continue, 2 optimal, 3 / 5 infeasible variants, 6, 7 iteration limit, 0xb time, 0xd numerical. Returns `(code, deciding return number)`.
pub fn status(t: &TermIn) -> (i32, u32) {
    let it = t.it;
    let iu = it as usize;
    let nx = t.q220.sqrt();
    let nz = t.q228.sqrt();
    let gap = (t.q130 - t.q140).abs();
    let mut relgap = gap / ((t.q130 - t.off).abs() + 1.0);
    if t.q2a0 == 0 {
        let alt = gap / (100.0 * t.q130.abs() + 1.0);
        if relgap <= alt {
            relgap = alt;
        }
    }
    let c17 = (t.q1f0 * t.ncomp) / (t.q138.abs() + 1.0);
    let c23 = (t.ncomp * t.q1f8) / (t.q130.abs() + 1.0);
    let stall = (it > 0 && t.tol < t.q1f8 && t.h78[iu] < 0.5 * t.h78[iu - 1]) || t.p2dd4;
    let n18 = (t.q198 * t.q198 + t.q1c0 * t.q1c0 + t.q1d0 * t.q1d0).sqrt();
    let b6 = n18 < t.q260 || t.q1d8 < t.q270;
    let b7 = t.q1a8 < t.q268 || t.q1e0 < t.q278;
    let tol = t.tol;
    if it > 5 && 1e20 < t.q1f0 && 0.99 * t.h80[iu - 1] < t.q1f0 {
        return (0xd, 1);
    }
    if 1e100 < t.q1f8 || 1e100 < t.q220 || 1e100 < t.q228 {
        return (0xd, 2);
    }
    for v in [t.q1f8, t.q220, t.q228] {
        if v.is_nan() || v.is_infinite() {
            return (0xd, 3);
        }
    }
    if t.q128 < t.q1e8 {
        return (6, 4);
    }
    let kap_small = t.kappa < 100.0 * tol * t.tau;
    let a = relgap < tol && c23 < 100.0 * tol && b6 && b7 && kap_small && !stall;
    let b = c17 < tol && c23 < 100.0 * tol && relgap < 0.1 && b6 && b7 && !stall;
    let c = relgap < tol
        && c17 < tol
        && b6
        && b7
        && it > 1
        && !t.p2dd4
        && 0.5 * t.h78[iu - 1] < t.h78[iu];
    let d = c17 < tol && t.p2dd4 && t.q1d8 < t.p2460 && kap_small;
    let e = t.q1f8 < 1e-16 && b6 && b7;
    if a || b || c || d || e {
        return (2, 5);
    }
    if t.q1f0 < tol && t.tau < tol && 1e6 * t.tau < t.kappa && t.q1a0 < t.q260 && t.q1b0 < t.q268 {
        if 1e-10 < t.q148 && t.q138 < -1e-10 {
            return (5 - 2 * (-t.q138 < t.q148) as i32, 6);
        }
        if 1e-10 < t.q148 {
            return (3, 7);
        }
        if t.q138 < -1e-10 {
            return (5, 8);
        }
    }
    if t.tau < tol && t.q1d8 < 0.01 * t.p2460 * nx && t.q138 < -1000.0 * tol * nx {
        return (5, 9);
    }
    if t.tau < tol
        && 1e6 * t.tau < t.kappa
        && t.q1e0 < 0.01 * t.p2460 * nz
        && 1000.0 * tol * nz < t.q148
    {
        return (3, 10);
    }
    if it > 10
        && t.kappa < 100.0 * tol * t.tau
        && relgap < 10.0 * tol
        && b6
        && b7
        && 0.8 * t.h80[iu - 3] < t.q1f0
    {
        return (2, 11);
    }
    if t.p2bbc <= it {
        return (7, 12);
    }
    if t.time_up {
        return (0xb, 13);
    }
    if t.p2bbc != 1000 && t.p2c10 == 0 {
        return (1, 14);
    }
    if t.if8.max(t.i88) + 200 < it {
        return (0xd, 15);
    }
    let d1 = t.q160 + 1.0;
    let d2 = 1.0 + t.q170;
    let m = |j: usize| (t.h50[j] / d1).max(t.h68[j] / d2);
    let late;
    if it < 11 {
        if it > 5 {
            if m(iu) < 0.9 * m(iu - 5) && t.q1f0 <= 0.99 * t.h78[iu - 5] {
                return (1, 16);
            }
            if m(iu) < 0.9 * m(iu - 1) && t.q1f0 <= 0.99 * t.h78[iu - 1] {
                return (1, 17);
            }
            late = t.if8 + 0x14 < it;
        } else if it > 1 {
            if m(iu) < 0.9 * m(iu - 1) && t.q1f0 <= 0.99 * t.h78[iu - 1] {
                return (1, 18);
            }
            late = t.if8 + 0x14 < it;
        } else {
            late = it > t.if8 + 0x14;
        }
    } else {
        if 1e-10 < t.q1f0 && t.q1f0 < t.h80[iu - 10] * 0.9 {
            return (1, 19);
        }
        if m(iu) < 0.9 * m(iu - 5) && t.q1f0 <= 0.99 * t.h78[iu - 5] {
            return (1, 20);
        }
        if m(iu) < 0.9 * m(iu - 1) && t.q1f0 <= 0.99 * t.h78[iu - 1] {
            return (1, 21);
        }
        late = t.if8 + 0x14 < it;
    }
    if late {
        let lim = if t.q2a0 != 0 { t.if4 + 4 } else { t.if4 };
        if t.q2a0 != 0 && lim < it {
            return (0xd, 22);
        }
        if t.if4 + 0x13 < it {
            return (0xd, 23);
        }
    }
    if it > 0x1e {
        if t.q2a0 != 0 && t.i88 + 9 < it && t.ia8 + 9 < it && t.if8 <= t.i88 && t.if8 <= t.ia8 {
            return (0xd, 24);
        }
        if it > 0x32 && t.i88 + 9 < it && t.ia8 + 9 < it && t.if8 <= t.i88 && t.if8 <= t.ia8 {
            return (0xd, 25);
        }
    }
    if (t.ncomp < 1.0 || 1e-16 <= t.q1f8) && (t.q5e0 != 0.0 || t.q5e8 != 0.0 || it <= t.if8 + 2) {
        return (1, 26);
    }
    (0xd, 27)
}

// ---------------------------------------------------------------------------------------------------------------------------------
// Linear algebra behind the iteration.

/// Right-hand sides of one direction solve (`Af dx = rp`, `Af'dy + dz - dv = rd`, `dx + du = ru`, `z dx + x dz = rcx`, `v du + u dv = rcu`).
#[derive(Debug, Clone)]
pub struct Rhs {
    pub rp: Vec<f64>,
    pub rd: Vec<f64>,
    pub ru: Vec<f64>,
    pub rcx: Vec<f64>,
    pub rcu: Vec<f64>,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum DirKind {
    Predictor,
    Corrector,
}

/// The two solves of the starting point (`0x2e8ef0`, factor with `D = start_diag`).
#[derive(Debug, Clone, Copy)]
pub enum StartKind<'a> {
    /// `A D A' y = b - Af w`, returns `(x = w + D Af' y, y)` (`w` = half the upper bound of a boxed column).
    Primal { w: &'a [f64] },
    /// `A D A' y = Af D c`, returns `(r = Af' y - c, y)` over all `n` variables.
    Dual,
}

/// Factor of `Af D Af'` and the solves with it (`0x2f1720` + `0x30fc90` / `0x305f10`).
pub trait NormalSolver {
    fn factor(&mut self, f: &StdForm, d: &[f64]) -> Result<(), BarrierError>;
    /// The factor of an iteration: `zx` (`Q+0x800`) is also the diagonal `-zx` of the dense-column nodes (the starting point's factors use `-1/D`).
    fn factor_iter(&mut self, f: &StdForm, d: &[f64], _zx: &[f64]) -> Result<(), BarrierError> {
        self.factor(f, d)
    }
    fn start(&mut self, f: &StdForm, kind: StartKind)
        -> Result<(Vec<f64>, Vec<f64>), BarrierError>;
    fn direction(
        &mut self,
        f: &StdForm,
        s: &State,
        kind: DirKind,
        rhs: &Rhs,
    ) -> Result<Direction, BarrierError>;
    /// Per iteration, before the directions: `Q+0x800` (`zx`) and the refinement threshold `Q+0x280` for the solvers that refine.
    fn prepare(&mut self, _zx: &[f64], _thr: f64) {}
    /// `Q+0x770 / 0x780` of the last direction: the norms of its solve residual (0 where the solver does not refine).
    fn norms(&self) -> (f64, f64) {
        (0.0, 0.0)
    }
    /// Number of trial weights of the corrector selection at iteration `it` (`Q+0x29c`). The reference's adaptation rule is not decoded
    /// (5 at iteration 0 and 3 later on AFIRO, SC50A, BLEND, BORE3D, DEGEN2; 3 at iteration 0 on ADLITTLE), so a replaying solver supplies it.
    fn ntry(&mut self, _it: i32, default: usize) -> usize {
        default
    }
    /// Inputs of the release profile's work model (`NZ = SYM+0x38`, `S = SYM+0x168`; `None` for solvers without a symbolic factor).
    fn work_info(&self, _f: &StdForm) -> Option<(f64, f64)> {
        None
    }
    /// `SYM+0x178 == 1` after the last factorization (the solves refine their result).
    fn refine_mode(&self, _f: &StdForm) -> bool {
        false
    }
    /// Number of correction solves of the refinement since the last call.
    fn take_extra_solves(&mut self) -> usize {
        0
    }
    /// The plain solve `(dx, dy)` for `rx` over the variables and `ry` over the rows (`0x305f10`; the homogeneous algorithm's solve with `(c, b)`).
    fn solve_raw(
        &mut self,
        _f: &StdForm,
        _rx: &[f64],
        _ry: &[f64],
    ) -> Result<(Vec<f64>, Vec<f64>), BarrierError> {
        unsupported("the homogeneous algorithm needs the exact solver")
    }
}

/// Plain f64 dense `LDL'` of `Af D Af'` in the permuted row order. **Not bit-exact** with the reference (its supernodal accumulation order
/// is not recovered); for tests at tolerance only. A pivot below `1e-13` times the largest diagonal is treated as a dependent row (its
/// unknown is set to 0).
#[derive(Default)]
pub struct DenseSolver {
    d: Vec<f64>,
    l: Vec<f64>,
    dinv: Vec<f64>,
    m: usize,
}

impl DenseSolver {
    fn solve_normal(&self, g: &[f64]) -> Vec<f64> {
        let m = self.m;
        let mut w = g.to_vec();
        for k in 0..m {
            for i in k + 1..m {
                w[i] -= self.l[i * m + k] * w[k];
            }
        }
        for k in 0..m {
            w[k] *= self.dinv[k];
        }
        for k in (0..m).rev() {
            for i in k + 1..m {
                w[k] -= self.l[i * m + k] * w[i];
            }
        }
        w
    }

    /// `Af theta v` for a vector over the variables.
    fn af_theta(&self, f: &StdForm, v: &[f64]) -> Vec<f64> {
        let mut out = vec![0.0; f.m];
        for j in 0..f.n2 {
            for k in f.col_start[j]..f.col_start[j + 1] {
                out[f.row_idx[k]] += f.val[k] * self.d[j] * v[j];
            }
        }
        for i in 0..f.m {
            if f.sense[i] != b'=' {
                out[i] += self.d[f.n2 + i] * v[f.n2 + i];
            }
        }
        out
    }

    fn at_y(f: &StdForm, y: &[f64]) -> Vec<f64> {
        let mut t = vec![0.0; f.n()];
        for j in 0..f.n2 {
            t[j] = f.aty_col(j, y);
        }
        for i in 0..f.m {
            if f.sense[i] != b'=' {
                t[f.n2 + i] = y[i];
            }
        }
        t
    }
}

impl NormalSolver for DenseSolver {
    fn factor(&mut self, f: &StdForm, d: &[f64]) -> Result<(), BarrierError> {
        let m = f.m;
        self.m = m;
        self.d = d.to_vec();
        let mut a = vec![0.0; m * m];
        for j in 0..f.n2 {
            let r = f.col_start[j]..f.col_start[j + 1];
            for p in r.clone() {
                for q in r.clone() {
                    a[f.row_idx[p] * m + f.row_idx[q]] += f.val[p] * d[j] * f.val[q];
                }
            }
        }
        for i in 0..m {
            if f.sense[i] != b'=' {
                a[i * m + i] += d[f.n2 + i];
            }
        }
        let maxd = (0..m).map(|i| a[i * m + i]).fold(0.0, f64::max);
        let mut l = vec![0.0; m * m];
        let mut dinv = vec![0.0; m];
        let mut dd = vec![0.0; m];
        for k in 0..m {
            let mut piv = a[k * m + k];
            for t in 0..k {
                piv -= l[k * m + t] * l[k * m + t] * dd[t];
            }
            if piv <= 1e-13 * maxd {
                continue; // dependent row: dinv = 0, column of L = 0
            }
            dd[k] = piv;
            dinv[k] = 1.0 / piv;
            for i in k + 1..m {
                let mut s = a[i * m + k];
                for t in 0..k {
                    s -= l[i * m + t] * l[k * m + t] * dd[t];
                }
                l[i * m + k] = s / piv;
            }
        }
        self.l = l;
        self.dinv = dinv;
        Ok(())
    }

    fn start(
        &mut self,
        f: &StdForm,
        kind: StartKind,
    ) -> Result<(Vec<f64>, Vec<f64>), BarrierError> {
        let n = f.n();
        match kind {
            StartKind::Primal { w } => {
                let aw = f.af(w);
                let g: Vec<f64> = (0..f.m).map(|i| f.b[i] - aw[i]).collect();
                let y = self.solve_normal(&g);
                let t = Self::at_y(f, &y);
                let x = (0..n).map(|j| w[j] + self.d[j] * t[j]).collect();
                Ok((x, y))
            }
            StartKind::Dual => {
                let mut cf = f.c.clone();
                cf.resize(n, 0.0);
                let g = self.af_theta(f, &cf);
                let y = self.solve_normal(&g);
                let t = Self::at_y(f, &y);
                let r = (0..n).map(|j| t[j] - cf[j]).collect();
                Ok((r, y))
            }
        }
    }

    fn direction(
        &mut self,
        f: &StdForm,
        s: &State,
        _kind: DirKind,
        rhs: &Rhs,
    ) -> Result<Direction, BarrierError> {
        let (n, n2, m) = (f.n(), f.n2, f.m);
        let mut w = vec![0.0; n];
        for j in 0..n {
            if f.flags[j] & 1 != 0 {
                w[j] = rhs.rd[j] - rhs.rcx[j] / s.x[j];
                if j < n2 && f.flags[j] & 2 != 0 {
                    w[j] += rhs.rcu[j] / s.u[j] - s.v[j] * rhs.ru[j] / s.u[j];
                }
            }
        }
        let aw = self.af_theta(f, &w);
        let g: Vec<f64> = (0..m).map(|i| rhs.rp[i] + aw[i]).collect();
        let dy = self.solve_normal(&g);
        let t = Self::at_y(f, &dy);
        let mut dir = Direction {
            dx: vec![0.0; n],
            dy,
            dz: vec![0.0; n],
            du: vec![0.0; n2],
            dv: vec![0.0; n2],
        };
        for j in 0..n {
            if f.flags[j] & 1 != 0 {
                dir.dx[j] = self.d[j] * (t[j] - w[j]);
                dir.dz[j] = (rhs.rcx[j] - s.z[j] * dir.dx[j]) / s.x[j];
                if j < n2 && f.flags[j] & 2 != 0 {
                    dir.du[j] = rhs.ru[j] - dir.dx[j];
                    dir.dv[j] = (rhs.rcu[j] - s.v[j] * dir.du[j]) / s.u[j];
                }
            }
        }
        Ok(dir)
    }
}

/// Captured solve outputs fed back in the order the iteration asks for them (starting point: primal, dual; then predictor and corrector per
/// iteration). Test-only input of `examples/barrier_replay.rs`; the factor call is ignored.
#[derive(Default)]
pub struct ReplaySolver {
    starts: Vec<(Vec<f64>, Vec<f64>)>,
    dirs: Vec<Direction>,
    ntry: Vec<Option<usize>>,
    ns: usize,
    nd: usize,
}

fn floats(v: &Value) -> Vec<f64> {
    v.as_array()
        .map(|a| a.iter().map(|x| x.as_f64().unwrap()).collect())
        .unwrap_or_default()
}

fn direction_json(v: &Value, sfx: &str) -> Direction {
    Direction {
        dx: floats(&v[format!("dx{sfx}")]),
        dy: floats(&v[format!("dy{sfx}")]),
        dz: floats(&v[format!("dz{sfx}")]),
        du: floats(&v[format!("du{sfx}")]),
        dv: floats(&v[format!("dv{sfx}")]),
    }
}

impl ReplaySolver {
    /// `{"start": {"x1", "y1", "r2", "y2"}, "iters": [{"pred": {dx, dy, dz, du, dv}, "corr": {dx2, dy2, dz2, du2, dv2}, "ntry": n}, ...]}`
    /// (`dy` over the rows, `du`, `dv` over the columns, the others over all variables; `ntry` optional).
    pub fn from_json(v: &Value) -> ReplaySolver {
        let s = &v["start"];
        let starts = vec![
            (floats(&s["x1"]), floats(&s["y1"])),
            (floats(&s["r2"]), floats(&s["y2"])),
        ];
        let mut dirs = Vec::new();
        let mut ntry = Vec::new();
        for it in v["iters"].as_array().map(|a| a.as_slice()).unwrap_or(&[]) {
            dirs.push(direction_json(&it["pred"], ""));
            dirs.push(direction_json(&it["corr"], "2"));
            ntry.push(it["ntry"].as_u64().map(|n| n as usize));
        }
        ReplaySolver {
            starts,
            dirs,
            ntry,
            ns: 0,
            nd: 0,
        }
    }

    /// Whether every recorded solve was used.
    pub fn consumed(&self) -> bool {
        self.ns == self.starts.len() && self.nd == self.dirs.len()
    }
}

impl NormalSolver for ReplaySolver {
    fn factor(&mut self, _f: &StdForm, _d: &[f64]) -> Result<(), BarrierError> {
        Ok(())
    }

    fn start(
        &mut self,
        _f: &StdForm,
        _kind: StartKind,
    ) -> Result<(Vec<f64>, Vec<f64>), BarrierError> {
        let r =
            self.starts.get(self.ns).cloned().ok_or_else(|| {
                BarrierError::Solver("replay: no more starting-point solves".into())
            })?;
        self.ns += 1;
        Ok(r)
    }

    fn direction(
        &mut self,
        _f: &StdForm,
        _s: &State,
        _kind: DirKind,
        _rhs: &Rhs,
    ) -> Result<Direction, BarrierError> {
        let r = self
            .dirs
            .get(self.nd)
            .cloned()
            .ok_or_else(|| BarrierError::Solver("replay: no more direction solves".into()))?;
        self.nd += 1;
        Ok(r)
    }

    fn ntry(&mut self, it: i32, default: usize) -> usize {
        self.ntry
            .get(it as usize)
            .copied()
            .flatten()
            .unwrap_or(default)
    }
}

// ---------------------------------------------------------------------------------------------------------------------------------
// The iteration.

#[derive(Debug, Clone, Copy)]
pub struct Options {
    /// `BarConvTol` (`Q+0x258`).
    pub tol: f64,
    /// Iteration limit (`P+0x2bbc`).
    pub max_iter: i32,
    /// Trial weights of the first corrector (`Q+0x29c`): 5 at iteration 0, 3 later (the adaptation rule is not decoded).
    pub ntry_first: usize,
    pub ntry_later: usize,
    /// Release profile (`BarCorrectors = -1`): the corrector limit and the trial weights adapt (`reverse-engineer/specs/barrier-release.md`).
    pub release: bool,
}

impl Default for Options {
    fn default() -> Options {
        Options {
            tol: 1e-8,
            max_iter: 1000,
            ntry_first: 5,
            ntry_later: 3,
            release: false,
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum BarrierStatus {
    Optimal,
    IterationLimit,
    /// A status code of `0x2e2fe0` that this port does not interpret further (3 / 5 infeasible variants, 6, 0xb, 0xd).
    Other(i32),
}

/// Per-iteration decisions (`BPRED` / `BAFF` / `BCORR` / `BSEL` / `BUPD` of the reference captures).
#[derive(Debug, Clone, Default)]
pub struct Decisions {
    pub ap: f64,
    pub ad: f64,
    pub aff: f64,
    pub target: f64,
    pub tp: f64,
    pub td: f64,
    pub wp: f64,
    pub wd: f64,
    pub ntry: usize,
    pub gamma: f64,
}

#[derive(Debug, Clone)]
pub struct IterRecord {
    pub it: i32,
    pub state: State,
    pub res: Residuals,
    pub status: i32,
    /// Absent at the starting point (it = -1) and at the final iteration.
    pub dec: Option<Decisions>,
}

#[derive(Debug, Clone)]
pub struct BarrierRun {
    pub status: BarrierStatus,
    pub iterations: i32,
    pub state: State,
    pub records: Vec<IterRecord>,
}

fn vec_json(v: &[f64]) -> Value {
    json!(v)
}

fn emit_record(r: &IterRecord) {
    let k = "barrier";
    if !dump::on(k) {
        return;
    }
    let mut rec = json!({
        "it": r.it, "x": vec_json(&r.state.x), "z": vec_json(&r.state.z), "y": vec_json(&r.state.y),
        "u": vec_json(&r.state.u), "v": vec_json(&r.state.v), "mu": r.res.q1f0, "pobj": r.res.q130, "dobj": r.res.q140,
        "rp": r.res.q1d8, "rd": r.res.q1e0, "status": r.status,
    });
    if let Some(d) = &r.dec {
        let o = rec.as_object_mut().unwrap();
        for (name, v) in [
            ("ap", d.ap),
            ("ad", d.ad),
            ("aff", d.aff),
            ("target", d.target),
            ("tp", d.tp),
            ("td", d.td),
            ("wp", d.wp),
            ("wd", d.wd),
            ("gamma", d.gamma),
        ] {
            o.insert(name.into(), json!(v));
        }
        o.insert("ntry".into(), json!(d.ntry));
    }
    dump::emit(k, rec);
}

fn term_input(
    f: &StdForm,
    sc: &Setup,
    o: &Options,
    it: i32,
    res: &Residuals,
    ncomp: usize,
    h: &Hist,
    best_dual: f64,
    q5: (f64, f64),
    tk: (f64, f64),
) -> TermIn {
    TermIn {
        it,
        q130: res.q130,
        q138: res.q138,
        q140: res.q140,
        q148: res.q148,
        q128: 1e100,
        q1e8: best_dual,
        off: f.offset,
        q1f0: res.q1f0,
        q1f8: res.q1f8,
        tol: o.tol,
        ncomp: ncomp as f64,
        tau: tk.0,
        kappa: tk.1,
        q198: res.q198,
        q1a8: res.q1a8,
        q1d8: res.q1d8,
        q1e0: res.q1e0,
        q1a0: res.q1a0,
        q1b0: res.q1b0,
        q1c0: res.q1c0,
        q1d0: 0.0,
        q220: h.x2,
        q228: h.z2,
        q160: sc.q160,
        q170: sc.q170,
        q260: sc.q260,
        q268: sc.q268,
        q270: sc.q270,
        q278: sc.q270,
        q2a0: 0,
        q5e0: q5.0,
        q5e8: q5.1,
        i88: it,
        if8: 0,
        ia8: it,
        if4: it,
        p2dd4: false,
        p2460: 1e-6,
        p2bbc: o.max_iter,
        p2c10: 0,
        time_up: false,
        h50: h.h50.clone(),
        h68: h.h68.clone(),
        h78: h.h78.clone(),
        h80: h.h80.clone(),
    }
}

#[derive(Default)]
struct Hist {
    h50: Vec<f64>,
    h68: Vec<f64>,
    h78: Vec<f64>,
    h80: Vec<f64>,
    /// `Q+0x220` (|x|^2): the reference's value lags the iterate (undecoded); only the 1e100 / NaN tests read it.
    x2: f64,
    z2: f64,
}

/// The second solve of the homogeneous algorithm (`0x3117e0`: `(c, b)`) and `den = kappa / tau + b'dy_s - c'dx_s`.
struct Hsd {
    cx: Vec<f64>,
    dxs: Vec<f64>,
    dys: Vec<f64>,
    den: f64,
}

impl Hsd {
    /// `0x30fc90` with kappa > 0: `num = ((G - b'dy) + c'dx) + Rk / tau`, `dtau = num / den`, `dx += dtau dx_s`, `dy += dtau dy_s`, `dz = (rcx - z dx) / x`
    /// again; returns `(dtau, dkappa = (Rk - dtau kappa) / tau)`.
    fn apply(
        &self,
        f: &StdForm,
        s: &State,
        rhs: &Rhs,
        d: &mut Direction,
        g: f64,
        rk: f64,
        tau: f64,
        kappa: f64,
    ) -> (f64, f64) {
        let (m, n) = (f.m, f.n());
        let num =
            ((g - dot_lanes(&f.b[..m], &d.dy[..m], 0)) + dot_lanes(&self.cx, &d.dx, 0)) + rk / tau;
        let dtau = num / self.den;
        let dkap = (rk - dtau * kappa) / tau;
        for j in 0..n {
            d.dx[j] = d.dx[j] + dtau * self.dxs[j];
        }
        for i in 0..m {
            d.dy[i] = d.dy[i] + dtau * self.dys[i];
        }
        for j in 0..n {
            if f.flags[j] & 1 != 0 {
                d.dz[j] = (rhs.rcx[j] - s.z[j] * d.dx[j]) / s.x[j];
            }
        }
        (dtau, dkap)
    }
}

// ---------------------------------------------------------------------------------------------------------------------------------
// Release profile: Gondzio correctors, the work counter and the adaptation of the corrector count (`reverse-engineer/specs/barrier-release.md`).

/// `0x2f0c86` (corrector `k >= 2`): the projection of a trial product `p` on `[target / 10, 10 target]`.
fn centrality(target: f64, p: f64) -> f64 {
    let (lo, hi) = (target / 10.0, 10.0 * target);
    if lo > p {
        lo - p
    } else if p > hi {
        (hi - p).max(-10.0 * target)
    } else {
        0.0
    }
}

/// Right-hand side of corrector `k` for the direction `d` so far: `target - dx dz` (`steps = None`, `k = 1`) or, for `k >= 2` with the steps of the direction so far,
/// the centrality correction of the stretched trial point `(dx ap~ + x)(dz ad~ + z)`, `ap~ = min(1.5 ap + 0.3, 1.001)`, `ad~ = min(0.3 + 1.5 ad, 1.001)`.
fn corrector_rhs(
    f: &StdForm,
    s: &State,
    d: &Direction,
    target: f64,
    steps: Option<(f64, f64)>,
) -> Rhs {
    let (n, n2) = (f.n(), f.n2);
    let (rcx, rcu): (Vec<f64>, Vec<f64>) = match steps {
        None => (
            (0..n)
                .map(|j| {
                    if f.flags[j] & 1 != 0 {
                        target - d.dx[j] * d.dz[j]
                    } else {
                        0.0
                    }
                })
                .collect(),
            (0..n2)
                .map(|j| {
                    if f.flags[j] & 2 != 0 {
                        target - d.du[j] * d.dv[j]
                    } else {
                        0.0
                    }
                })
                .collect(),
        ),
        Some((tp, td)) => {
            let (apt, adt) = ((tp * 1.5 + 0.3).min(1.001), (0.3 + 1.5 * td).min(1.001));
            (
                (0..n)
                    .map(|j| {
                        if f.flags[j] & 1 != 0 {
                            centrality(target, (d.dx[j] * apt + s.x[j]) * (d.dz[j] * adt + s.z[j]))
                        } else {
                            0.0
                        }
                    })
                    .collect(),
                (0..n2)
                    .map(|j| {
                        if f.flags[j] & 2 != 0 {
                            centrality(target, (d.du[j] * apt + s.u[j]) * (d.dv[j] * adt + s.v[j]))
                        } else {
                            0.0
                        }
                    })
                    .collect(),
            )
        }
    };
    Rhs {
        rp: vec![0.0; f.m],
        rd: vec![0.0; n],
        ru: vec![0.0; n2],
        rcx,
        rcu,
    }
}

/// The work counter the driver reads through its timers: every site adds 0.75 times a count, a window adds 1.0 on top.
struct Work {
    m: f64,
    nc: f64,
    nub: f64,
    nnz: f64,
    nz: f64,
    s: f64,
    mode1: bool,
    /// Totals of the windows: factor (`W138`), predictor (`W130`), first corrector (`W148`), later correctors (`W140`), update (`W128`).
    t: [f64; 5],
}

const WORK_SCALE: f64 = 0.75;

impl Work {
    fn new(f: &StdForm, nz: f64, s: f64) -> Work {
        Work {
            m: f.m as f64,
            nc: f.n2 as f64,
            nub: (f.ub_hi - f.ub_lo) as f64,
            nnz: f.row_idx.len() as f64,
            nz,
            s,
            mode1: false,
            t: [0.0; 5],
        }
    }

    fn add(&mut self, k: usize, sites: &[f64]) {
        let mut w = 0.0;
        for x in sites {
            w = w + x;
        }
        self.t[k] = (self.t[k] + 1.0) + w;
    }

    /// The solve `0x30fc90` with `extra` correction solves of the refinement.
    fn solve_sites(&self, extra: usize) -> Vec<f64> {
        let blk = [
            WORK_SCALE * (2.0 * self.nnz),
            WORK_SCALE * (2.0 * self.nz),
            WORK_SCALE * (2.0 * self.nnz),
            (self.m * 11.0 + self.nc * 4.0) * WORK_SCALE,
        ];
        let res = WORK_SCALE * (4.0 * self.nnz);
        let mut v = Vec::new();
        for _ in 0..=extra {
            v.extend_from_slice(&blk);
            if self.mode1 {
                v.push(res);
            }
        }
        v.push(((self.m + self.nc) * 5.0 + self.nub * 9.0) * WORK_SCALE);
        v
    }

    fn step_site(&self) -> f64 {
        (((self.m + self.nc) * 7.0) + self.nub * 6.0) * WORK_SCALE
    }

    fn rhs_site(&self) -> f64 {
        (((self.m + self.nub) * 3.0) + self.nc * 13.0) * WORK_SCALE
    }

    fn factor(&mut self) {
        let s = self.s * WORK_SCALE;
        self.add(0, &[s]);
    }

    fn predictor(&mut self, extra: usize) {
        let mut v = self.solve_sites(extra);
        v.push(self.step_site());
        self.add(1, &v);
    }

    fn corrector(&mut self, k: usize, trials: usize, extra: usize) {
        let mut v = vec![self.rhs_site()];
        v.extend(self.solve_sites(extra));
        v.extend(std::iter::repeat(self.step_site()).take(trials));
        self.add(if k == 1 { 2 } else { 3 }, &v);
    }

    fn update(&mut self) {
        self.add(4, &[]);
    }

    fn ratio(&self) -> f64 {
        let t = &self.t;
        t[3] / (((t[0] + t[1]) + t[2]) + t[4])
    }
}

/// `0x2e22bf`: `max(1, floor(log10(max(columns, 1))) - 2)`.
fn corrector_cap(ncol: usize) -> usize {
    let e = (ncol.max(1) as f64).log10().floor() as i64 - 2;
    if e > 0 {
        e as usize
    } else {
        1
    }
}

/// `0x2e22bf..0x2e24d2`: the corrector limit `Q+0x298` and the trial weights `Q+0x29c` of the next iteration from the work ratio and the mean step product of the last three iterations.
fn adapt_correctors(q298: usize, cap: usize, ratio: f64, avg: f64) -> (usize, usize) {
    let (mut inc, mut dec) = (false, false);
    if q298 < 20 && 0.25 > ratio {
        inc = true;
    } else if q298 < cap {
        if 0.5 > ratio && 0.1 > avg {
            inc = true;
        } else if q298 > 1 && ratio > 0.25 && avg > 0.1 {
            dec = true;
        }
    } else if q298 > cap {
        dec = ratio > 0.25;
    } else {
        dec = q298 > 1 && ratio > 0.25 && avg > 0.1;
    }
    let r = if inc {
        q298 + 1
    } else if dec {
        q298 - 1
    } else {
        q298
    };
    (r, (r + 1).min(10).max(3))
}

/// Runs the barrier iteration on `f` with `solver` behind the factor and the solves.
///
/// Gaps (H): the progress counters `Q+0x88/0xa8/0xf4` are modelled as "improved at every iteration" (true on the five models), `Q+0x5e0/0x5e8`
/// as nonzero, `Q+0x1e8` as the best dual objective, `Q+0x220/0x228` as the norms of the iterate. They matter only in branches of
/// `0x2e2fe0` that the corpus does not reach.
pub fn run<S: NormalSolver>(
    f: &StdForm,
    solver: &mut S,
    o: &Options,
) -> Result<BarrierRun, BarrierError> {
    // The objective scale `Q+0x158` (the release profile's `f`).
    let scale = f.scaling.as_ref().map_or(1.0, |s| s.obj);
    let release = o.release;
    // Work model of the corrector adaptation (release only).
    let mut work = if release {
        let (nz, sv) = solver.work_info(f).ok_or_else(|| {
            BarrierError::Unsupported("the release profile needs the exact solver".into())
        })?;
        Some(Work::new(f, nz, sv))
    } else {
        None
    };
    // The homogeneous self-dual scalars `Q+0x690 / 0x698` (kappa 0: not started; started after restart A) and the number of complementarity terms.
    let (mut tau, mut kappa) = (1.0f64, 0.0f64);
    let mut ncomp = f.ncomp();
    let mut q2c = false;
    let c_file = f.c.clone();
    let sc = setup_scalars(&f.b_file, &c_file, o.tol);
    // Starting point.
    solver.factor(f, &start_diag(f))?;
    if let Some(wk) = work.as_mut() {
        wk.factor();
    }
    let w: Vec<f64> = (0..f.n())
        .map(|j| {
            if j < f.n2 && f.flags[j] == 3 {
                0.5 * f.ub[j]
            } else {
                0.0
            }
        })
        .collect();
    let s1 = solver.start(f, StartKind::Primal { w: &w })?;
    let s2 = solver.start(f, StartKind::Dual)?;
    let (mut st, p3_0, p4_0, _stage) = start_point(f, s1.clone(), s2.clone(), &sc, 0, 0.0);
    // `Q+0x38` of the starting point: the floor of a type B restart is `max(1000, 10 Q+0x38)` (`0x2e2983`).
    let q38 = p3_0.min(p4_0);
    let mut records = Vec::new();
    let res = residuals(f, &st, tau, kappa, ncomp, scale);
    let rec = IterRecord {
        it: -1,
        state: st.clone(),
        res: res.clone(),
        status: 1,
        dec: None,
    };
    emit_record(&rec);
    records.push(rec);
    let mut rd = res.rd.clone();
    z_reduce(f, &mut st, &mut rd, res.q1f0);

    let mut h = Hist::default();
    let mut best_dual = -1e100;
    // `Q+0x5e0 / 0x5e8 / 0x5f0`: the step products of the last three iterations.
    let mut w5 = [1e100f64; 3];
    // `Q+0xf4`: the iteration after the last update with a step product above 5e-5; `Q+0xf8`: the iteration after the restart.
    let (mut qf4, mut qf8) = (0i32, 0i32);
    // `Q+0x298` (corrector limit), `Q+0x29c` of the next iteration.
    let (mut q298, mut ntry_next) = (1usize, 3usize);
    let mut it = 0;
    loop {
        let mut res = residuals(f, &st, tau, kappa, ncomp, scale);
        if it > 10 && h.h78[it as usize - 1] > 0.99 * h.h78[it as usize - 2] {
            // `0x2ecde0` at `0x2eed95`: complementarity decreased by less than 1% over the last iteration: the z reduction of the starting point with the
            // current mu, after the scalars were computed.
            z_reduce(f, &mut st, &mut res.rd, res.q1f0);
        }
        h.h50.push(res.q198);
        h.h68.push(res.q1a8);
        h.h78.push(res.q1f8);
        h.h80.push(res.q1f0);
        h.x2 = dot_lanes(&st.x, &st.x, 0);
        h.z2 = dot_lanes(&st.z, &st.z, 0);
        if res.q148 > best_dual && it > 0 {
            best_dual = res.q148;
        }
        let ti = term_input(
            f,
            &sc,
            o,
            it,
            &res,
            ncomp,
            &h,
            best_dual,
            (w5[0], w5[1]),
            (tau, kappa),
        );
        let (code, _why) = status(&ti);
        if code != 1 || it >= o.max_iter + 1 {
            let rec = IterRecord {
                it,
                state: st.clone(),
                res,
                status: code,
                dec: None,
            };
            emit_record(&rec);
            records.push(rec);
            let status = match code {
                2 => BarrierStatus::Optimal,
                7 => BarrierStatus::IterationLimit,
                c => BarrierStatus::Other(c),
            };
            return Ok(BarrierRun {
                status,
                iterations: it,
                state: st,
                records,
            });
        }
        let (n, n2, m) = (f.n(), f.n2, f.m);
        // The reference's iteration record (`BIT`) is taken before the homogeneous start recomputes the scalars.
        let res_bit = res.clone();
        if q2c && kappa == 0.0 {
            // `0x2ebf20` (homogeneous start after restart A): tau = 1, kappa = max(c'x - b'y + v.ub, 1) and sqrt(kappa mu) when mu > 0, one more
            // complementarity term (`Q+0x210`), `0x2ecde0` again.
            if f.has_ub() {
                return unsupported("the homogeneous algorithm with upper-bounded columns");
            }
            tau = 1.0;
            let g = dot_lanes(&f.c[..n2], &st.x[..n2], 0) - dot_lanes(&f.b[..m], &st.y[..m], 0);
            kappa = g.max(1.0);
            if res.q1f0 > 0.0 {
                kappa = (kappa * res.q1f0).sqrt();
            }
            ncomp += 1;
            res = residuals(f, &st, tau, kappa, ncomp, scale);
            let i = it as usize;
            h.h50[i] = res.q198;
            h.h68[i] = res.q1a8;
            h.h78[i] = res.q1f8;
            h.h80[i] = res.q1f0;
        }
        let hsd = kappa > 0.0;
        // Predictor.
        let (zx, dd) = normal_diag(f, &st, res.q130);
        solver.factor_iter(f, &dd, &zx)?;
        if let Some(wk) = work.as_mut() {
            wk.factor();
            wk.mode1 = solver.refine_mode(f);
        }
        solver.prepare(
            &zx,
            ((res.q198 * res.q198 + res.q1a8 * res.q1a8).sqrt() / 10.0).max(1e-8),
        );
        let rcx: Vec<f64> = (0..n)
            .map(|j| {
                if f.flags[j] & 1 != 0 {
                    -(st.x[j] * st.z[j])
                } else {
                    0.0
                }
            })
            .collect();
        let rcu: Vec<f64> = (0..n2)
            .map(|j| {
                if f.flags[j] & 2 != 0 {
                    -(st.u[j] * st.v[j])
                } else {
                    0.0
                }
            })
            .collect();
        let rhs = Rhs {
            rp: res.rp.clone(),
            rd: res.rd.clone(),
            ru: res.ru.clone(),
            rcx,
            rcu,
        };
        let mut pred = solver.direction(f, &st, DirKind::Predictor, &rhs)?;
        let (mut q770, mut q780) = solver.norms();
        if let Some(wk) = work.as_mut() {
            let extra = solver.take_extra_solves();
            wk.predictor(extra);
        }
        // Homogeneous step (`0x30fc90` with kappa > 0): the second solve with `(c, b)` (`0x3117e0`), `dtau = num / den`, `dx += dtau dx_s`, `dy += dtau dy_s`, then dz again.
        let mut hq = (0.0f64, 0.0f64);
        let mut hs: Option<Hsd> = None;
        if hsd {
            let mut cx = f.c.clone();
            cx.resize(n, 0.0);
            let (dxs, dys) = solver.solve_raw(f, &cx, &f.b)?;
            let den =
                ((kappa / tau) + dot_lanes(&f.b[..m], &dys[..m], 0)) - dot_lanes(&cx, &dxs, 0);
            let hd = Hsd { cx, dxs, dys, den };
            hq = hd.apply(
                f,
                &st,
                &rhs,
                &mut pred,
                (kappa - res.q148) + res.q138,
                -tau * kappa,
                tau,
                kappa,
            );
            hs = Some(hd);
        }
        let tkp: Vec<(f64, f64, f64)> = if hsd {
            vec![(tau, hq.0, 0.0), (kappa, hq.1, 0.0)]
        } else {
            Vec::new()
        };
        let ap0 = step_pair(f, &st, &pred, None, true, &tkp);
        let ad0 = step_pair(f, &st, &pred, None, false, &tkp);
        let aff = mu_aff(
            f,
            &st,
            &pred,
            ap0,
            ad0,
            ncomp,
            if hsd {
                Some((hq.0, hq.1, tau, kappa))
            } else {
                None
            },
        );
        let mut target = sigma_target(res.q1f0, aff);
        if it > qf4 {
            // `0x2e2047`: the last update with a step product above 5e-5 was not the previous iteration.
            target = res.q1f0 * 0.9;
        }
        // Corrector loop (`0x2e20e3..0x2e2276`): the direction so far starts as the predictor and becomes the combined direction after each selection.
        if !release {
            q298 = 1; // `0x2e20b5`: `Q+0x298 = max(BarCorrectors, 1)`
        }
        let ntry = solver.ntry(
            it,
            if it == 0 {
                o.ntry_first
            } else if release {
                ntry_next
            } else {
                o.ntry_later
            },
        );
        let (mut tp0, mut td0) = (ap0, ad0);
        let mut dirk = pred.clone();
        let mut k = 1usize;
        let mut hc = (0.0f64, 0.0f64);
        let (mut tp, mut td, mut wp, mut wd);
        let comb = loop {
            let rhs2 = corrector_rhs(
                f,
                &st,
                &dirk,
                target,
                if k > 1 { Some((tp0, td0)) } else { None },
            );
            let mut corr = solver.direction(f, &st, DirKind::Corrector, &rhs2)?;
            let (q778, q788) = solver.norms();
            let mut tks: Vec<(f64, f64, f64)> = Vec::new();
            if let Some(hd) = &hs {
                hc = hd.apply(
                    f,
                    &st,
                    &rhs2,
                    &mut corr,
                    0.0,
                    target - hq.0 * hq.1,
                    tau,
                    kappa,
                );
                tks = vec![(tau, hq.0, hc.0), (kappa, hq.1, hc.1)];
            }
            (tp, td, wp, wd) = select(f, &st, &dirk, &corr, tp0, td0, ntry, &tks, k);
            // `0x30f5xx..0x30fc3e`: the residual norms of the combined step.
            (q770, q780) = (q770.max(q778 * wp), q780.max(q788 * wd));
            let c = Direction {
                dx: combine(&dirk.dx, &corr.dx, wp),
                dz: combine(&dirk.dz, &corr.dz, wd),
                dy: combine(&dirk.dy[..f.m], &corr.dy[..f.m], wd),
                du: combine(&dirk.du, &corr.du, wp),
                dv: combine(&dirk.dv, &corr.dv, wd),
            };
            // `0x2e2176..0x2e2276`: leave when both steps reach 1.001, when the step product of a corrector k > 1 is below 1e-100, when (k > max(Q+0x298 / 10, 1))
            // the previous product exceeds 0.99 times the new one, or when a solve residual exceeds 1e30; the steps become `Q+0x5b8 / 0x5c8` in any case.
            let mut stop = tp >= 1.001 && td >= 1.001;
            if k > 1 && 1e-100 > td * tp {
                stop = true;
            }
            if k > (q298 / 10).max(1) && tp0 * td0 > (0.99 * tp) * td {
                stop = true;
            }
            if q770 > 1e30 || q780 > 1e30 {
                stop = true;
            }
            (tp0, td0) = (tp, td);
            if let Some(wk) = work.as_mut() {
                let extra = solver.take_extra_solves();
                wk.corrector(k, ntry - usize::from(k > 1 && ntry > 1), extra);
            }
            k += 1;
            if stop || hsd || k > q298 {
                break c;
            }
            dirk = c;
        };
        // `0x2e22bf..0x2e24d2`: adaptation of the corrector limit and of the trial weights of the next iteration (release).
        if let Some(wk) = work.as_mut() {
            let avg = if it > 2 {
                (w5[0] + w5[1] + w5[2]) * (1.0 / 3.0)
            } else {
                1.0
            };
            (q298, ntry_next) = adapt_correctors(q298, corrector_cap(f.n2), wk.ratio(), avg);
            wk.update();
        }
        // `0x2e24f0..0x2e2647`: with a small complementarity measure the steps are capped by the solve residuals against the iterate's.
        let (mut cp, mut cd) = (tp, td);
        if res.q1f8 < 1e-4 {
            let t1 = (res.q198 * 1.2).max(sc.q260);
            if q770 > t1 {
                cp = cp.min(t1 / q770);
            }
            let t2 = (res.q1a8 * 1.2).max(sc.q268);
            if q780 > t2 {
                cd = cd.min(t2 / q780);
            }
        }
        w5 = [cp * cd, w5[0], w5[1]];
        // Combination words `Q+0x720..0x750` (`0x30fb8d`): dtau and dkappa of the primal (wp) and dual (wd) candidates.
        let hk = if hsd {
            Some((
                tau,
                kappa,
                hq.0 + wp * hc.0,
                hq.1 + wp * hc.1,
                hq.0 + wd * hc.0,
            ))
        } else {
            None
        };
        let (mut next, (g, _sp, _sd), hout) = update(f, &st, &comb, cp, cd, ncomp, hk);
        if let Some(h) = hout {
            // `0x30e619`: the dual iterate is scaled by tau_p / tau_d when the two differ.
            if h.tau != h.tau_d {
                let sc_ = h.tau / h.tau_d;
                next.z.iter_mut().for_each(|t| *t *= sc_);
                next.y.iter_mut().for_each(|t| *t *= sc_);
                next.v.iter_mut().for_each(|t| *t *= sc_);
            }
            tau = h.tau;
            kappa = h.kappa;
        }
        let dec = Decisions {
            ap: ap0,
            ad: ad0,
            aff,
            target,
            tp,
            td,
            wp,
            wd,
            ntry,
            gamma: g,
        };
        let (q1d8, q1e0, q138, q1f0) = (res.q1d8, res.q1e0, res.q138, res.q1f0);
        let rec = IterRecord {
            it,
            state: st.clone(),
            res: res_bit,
            status: code,
            dec: Some(dec),
        };
        emit_record(&rec);
        records.push(rec);
        st = next;
        // `0x30d3e0` prologue: Q+0xf4 = it + 1 when the step product is above 5e-5.
        if (cp * (1.0 - 1e-6)).min(1.0) * ((1.0 - 1e-6) * cd).min(1.0) > 5e-5 {
            qf4 = it + 1;
        }
        // `0x2e2711..0x2e29d5`: restart of the starting point (only before the first restart): the last three step products below 1e-3, complementarity
        // stalled (the log gap history over 3 iterations), both solve residuals below 1e-3 and a relative gap above 1e-4.
        let mut relgap = ncomp as f64 * q1f0 / (q138.abs() + 1.0);
        if 1e-7 > q1d8 && 1e-7 > q1e0 {
            relgap = relgap * 0.1;
        }
        let ratio3 = if it - qf8 >= 3 {
            let den = h.h78[it as usize - 3];
            if den != 0.0 {
                h.h78[it as usize] / den
            } else {
                1.0
            }
        } else {
            0.0
        };
        let mut restart = qf8 == 0
            && 1e-3 > w5[2]
            && 1e-3 > w5[1]
            && 1e-3 > w5[0]
            && ratio3 > 0.98
            && 1e-3 > q770
            && 1e-3 > q780
            && relgap > 1e-4;
        // Type B restart floor `Q+0x30` (`0x2e2983`): `max(1000, 10 Q+0x38)`.
        let mut q30r = (q38 * 10.0).max(1000.0);
        if !restart && !q2c {
            // Restart A (`0x2e2aa1`, once; reached whenever B does not fire): |x| or |y| of the updated iterate huge (1e15 times the data norms), or from
            // iteration 10 on above 1e8 times them while the primal (dual) residual history did not fall. The starting point's floor stays 0.
            let (nx, ny) = (
                dot_lanes(&st.x, &st.x, 0).sqrt(),
                dot_lanes(&st.y, &st.y, 0).sqrt(),
            );
            let (t1, t2) = (sc.q168.max(1.0), sc.q178.max(1.0));
            let (mut fp, mut fd) = (false, false);
            if it >= 2 {
                let i = it as usize;
                fp = sc.q270 > q1d8 || 0.9 * h.h50[i - 2] > h.h50[i];
                fd = sc.q270 > q1e0 || 0.9 * h.h68[i - 2] > h.h68[i];
            }
            if nx > 1e15 * t1 || ny > 1e15 * t2 {
                restart = true;
            } else if it >= 10 {
                restart = (nx > 1e8 * t1 && !fp) || (ny > 1e8 * t2 && !fd);
            }
            if restart {
                q2c = true;
                q30r = 0.0;
            }
        }
        if release && !restart && it > qf4 {
            // `0x2e28ef..0x2e291f` (`BarCorrectors = -1`, no restart): the last update had a step product of at most 5e-5: at least 10 correctors and 10 trial weights next.
            q298 = q298.max(10);
            ntry_next = 10;
        }
        if restart {
            (st, _, _, _) = start_point(f, s1.clone(), s2.clone(), &sc, 0, q30r);
            if let Some(wk) = work.as_mut() {
                wk.factor();
            }
            qf8 = it + 1;
            let r = residuals(f, &st, tau, kappa, ncomp, scale);
            let mut rd = r.rd.clone();
            z_reduce(f, &mut st, &mut rd, r.q1f0);
        }
        it += 1;
    }
}

/// The solution of the solved columns (the model's, then the range columns; the lower-bound shift undone), the objective (`Q+0x130`) and the row duals in file order.
#[derive(Debug, Clone)]
pub struct Solution {
    pub run: BarrierRun,
    pub x: Vec<f64>,
    pub objective: f64,
    /// Row duals in file row order (the row sign of a '>' row undone).
    pub y: Vec<f64>,
    /// Reduced costs `z - v` of the model's columns.
    pub rc: Vec<f64>,
}

pub fn solution(f: &StdForm, run: BarrierRun) -> Solution {
    let x: Vec<f64> = (0..f.n2)
        .map(|j| {
            if f.lower[j] != 0.0 {
                run.state.x[j] + f.lower[j]
            } else {
                run.state.x[j]
            }
        })
        .collect();
    let last = run.records.last().unwrap();
    let mut y = vec![0.0; f.m];
    for k in 0..f.m {
        y[f.row_of[k]] = run.state.y[k] * if f.sense[k] == b'>' { -1.0 } else { 1.0 };
    }
    let mut rc: Vec<f64> = (0..f.n2)
        .map(|j| {
            (if f.flags[j] & 1 != 0 {
                run.state.z[j]
            } else {
                0.0
            }) - (if f.flags[j] & 2 != 0 {
                run.state.v[j]
            } else {
                0.0
            })
        })
        .collect();
    let mut objective = last.res.q130;
    let mut x = x;
    if let Some(sc) = &f.scaling {
        // Release profile: `x = C x'`, `y = R y' / f`, `rc = rc' / (C f)` (powers of two; the objective is already in the units of the model).
        for j in 0..f.n_struct {
            x[j] *= sc.col[j];
        }
        for k in 0..f.range_rank.len() {
            x[f.n_struct + k] *= sc.col[f.n_struct + f.range_rank[k]];
        }
        for (i, t) in y.iter_mut().enumerate() {
            *t = *t * sc.row[i] / sc.obj;
        }
        for (j, t) in rc.iter_mut().enumerate() {
            let c = if j < f.n_struct {
                sc.col[j]
            } else {
                sc.col[f.n_struct + f.range_rank[j - f.n_struct]]
            };
            *t = *t / (c * sc.obj);
        }
    }
    if f.max_sign < 0.0 {
        y.iter_mut().for_each(|t| *t = -*t);
        rc.iter_mut().for_each(|t| *t = -*t);
        objective = -objective;
    }
    Solution {
        objective,
        x,
        y,
        rc,
        run,
    }
}

/// Convenience driver: structure (ordering) by [`crate::order::structure`], conversion, iteration with the given solver.
pub fn solve_with<S: NormalSolver>(
    model: &LpModel,
    solver: &mut S,
    o: &Options,
) -> Result<Solution, BarrierError> {
    let s = crate::order::structure(model, crate::order::DENSE_COL_MIN)
        .map_err(|e| BarrierError::Unsupported(e.3.to_string()))?;
    let f = std_form(model, &s.sym.perm)?;
    let run = run(&f, solver, o)?;
    Ok(solution(&f, run))
}

/// [`solve_with`] with the bit-exact [`crate::chol::CholSolver`] (the controlled profile of the reference: `BarCorrectors = 1`, no scaling, `BarOrder = 0`).
pub fn solve_exact(model: &LpModel, o: &Options) -> Result<Solution, BarrierError> {
    let s = crate::order::structure(model, crate::order::DENSE_COL_MIN)
        .map_err(|e| BarrierError::Unsupported(e.3.to_string()))?;
    let f = std_form(model, &s.sym.perm)?;
    let mut solver = crate::chol::CholSolver::new(&f, &s.sym);
    let run = run(&f, &mut solver, o)?;
    Ok(solution(&f, run))
}

/// The release profile of the reference (`reverse-engineer/specs/barrier-release.md`): the model scaled by [`crate::scale::compute_with_ranges`] over [`range_columns`], adaptive correctors.
pub fn solve_release(model: &LpModel, o: &Options) -> Result<Solution, BarrierError> {
    let sc = crate::scale::compute_with_ranges(model, &range_columns(model));
    let s = crate::order::structure(model, crate::order::DENSE_COL_MIN)
        .map_err(|e| BarrierError::Unsupported(e.3.to_string()))?;
    let f = std_form_scaled(model, &s.sym.perm, Some(&sc))?;
    let mut solver = crate::chol::CholSolver::new(&f, &s.sym);
    let o = Options {
        release: true,
        ..o.clone()
    };
    let run = run(&f, &mut solver, &o)?;
    Ok(solution(&f, run))
}

/// `Method::Barrier` of [`crate::solve`]: [`solve_exact`] with the options mapped (`barrier_tol`, `max_iterations`), the result in the public types.
/// Status: optimal, the iteration limit (API 7), unbounded (code 5, API 5) and infeasible (code 3, never observed); infeasible-or-unbounded (code 0xd, API 12)
/// and the other codes (6, 0xb) are reported as `NumericalFailure` (no certificate is produced).
pub fn solve_method(
    model: &LpModel,
    opts: &crate::SolveOptions,
) -> Result<crate::SolveResult, crate::simplex::SolveError> {
    use crate::simplex::{SolveError, Status};
    let o = Options {
        tol: opts.barrier_tol,
        max_iter: opts.max_iterations.min(i32::MAX as usize) as i32,
        ..Options::default()
    };
    let sol = (if opts.scaling {
        solve_release(model, &o)
    } else {
        solve_exact(model, &o)
    })
    .map_err(|e| match e {
        BarrierError::Unsupported(s) | BarrierError::Solver(s) => SolveError::Unsupported(s),
    })?;
    let status = match sol.run.status {
        BarrierStatus::Optimal => Status::Optimal,
        BarrierStatus::IterationLimit => Status::IterationLimit,
        // 5: unbounded (API 5, fx_unbounded); 3: infeasible (never observed); 0xd: infeasible or unbounded (API 12, fx_infeasible), reported like PDHG's 12.
        BarrierStatus::Other(5) => Status::Unbounded,
        BarrierStatus::Other(3) => Status::Infeasible,
        BarrierStatus::Other(_) => Status::NumericalFailure,
    };
    Ok(crate::SolveResult {
        status,
        objective: sol.objective,
        x: sol.x[..model.num_cols()].to_vec(),
        dual: sol.y,
        reduced_cost: sol.rc[..model.num_cols()].to_vec(),
        basis: crate::Basis::slack(model.num_cols(), model.num_rows()),
        iterations: sol.run.iterations.max(0) as usize,
        refactorizations: 0,
        trace: Vec::new(),
        farkas: None,
        ray: None,
    })
}

/// [`solve_with`] on the plain f64 [`DenseSolver`] (not bit-exact, tests at tolerance).
pub fn solve_dense(model: &LpModel) -> Result<Solution, BarrierError> {
    solve_with(model, &mut DenseSolver::default(), &Options::default())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn ext_rounds_to_64_bits() {
        // 1 + 2^-64 is a tie in 64 bits and rounds to even (1); 1 + 3 * 2^-64 rounds up to 1 + 2^-62.
        let one = Ext::from_f64(1.0);
        let t = Ext::from_f64(2f64.powi(-64));
        assert_eq!(one.add(t).to_f64(), 1.0);
        let e = one.add(Ext::from_f64(3.0 * 2f64.powi(-64)));
        assert_eq!(e.sub(one).to_f64(), 2f64.powi(-62));
        // Product of two doubles is exact in 106 bits and rounded once to 64.
        let a = 1.0 + 2f64.powi(-52);
        let p = Ext::mul_f64(a, a);
        assert_eq!(p.sub(Ext::from_f64(1.0)).to_f64(), 2f64.powi(-51));
        // Cancellation is exact.
        assert_eq!(Ext::from_f64(3.5).sub(Ext::from_f64(3.5)).to_f64(), 0.0);
        assert_eq!(Ext::from_f64(-2.0).add(Ext::from_f64(0.75)).to_f64(), -1.25);
    }

    #[test]
    fn kahan_dot_is_more_accurate_than_double() {
        let a = [1.0, 1e-17, 1e-17, -1.0];
        let b = [1.0, 1.0, 1.0, 1.0];
        let seq = a.iter().zip(&b).fold(0.0, |s, (x, y)| s + x * y);
        assert_eq!(seq, 0.0);
        assert!((kahan_dot(&a, &b).to_f64() - 2e-17).abs() < 1e-19);
    }

    #[test]
    fn dot_lanes_matches_sequential_below_eight() {
        let a: Vec<f64> = (0..7).map(|i| 0.1 * i as f64 + 0.3).collect();
        let s: f64 = a.iter().fold(0.0, |s, v| s + v * v);
        assert_eq!(dot_lanes(&a, &a, 0), s);
    }
}
