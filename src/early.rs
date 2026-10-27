//! The early-centering solve `0x34a580` of the reference (`reverse-engineer/specs/barrier-solve.md`, "Early centering"): in place of the Cholesky factor
//! of `A D A'` it solves the normal equations with a diagonally preconditioned conjugate residual method, started from a random vector.
//! Same arguments and results as `0x305f10` (`rx` over the variables, `ry` over the rows; `(dx, dy)`), Threads = 1, mode 1 of the matrix product
//! (`2 threads nrow < nnz(A)`: the column scatter of `0x351160`), no dense-column nodes (`Q+0x28 = 0`).
//!
//! `w = rx D`, `b = A w + ry` (`0x2ec200`), `M = diag(A D A')` (column order), `x0 = 1e-8 u` with `u` the next `nrow` draws of the generator,
//! `r = b - A D A' x0`, `z = M^-1 r`, `p = z`, `Az = A D A' z`, `num = z.Az`, `Ap = Az`; then per step: `den = sum (Ap_i Ap_i) / M_i`,
//! `alpha = num / den`, `x += alpha p`, `r -= alpha Ap`, stop when `step > 2` and `tol > |r|`, `z = M^-1 r`, `Az`, `num' = z.Az`,
//! `beta = num' / num`, `p = z + beta p`, `Ap = Az + beta Ap`. The results are `dy = x`, `dx = (A' dy - rx) D` (`0x2ec480`).
use crate::barrier::{dot_lanes, StdForm};
use crate::chol::{dot8, Amat};

/// The reference's generator `0x88a8d0` (31 bits times 2^-31), state `ctx+0xb0`.
pub fn draw(state: &mut u64) -> f64 {
    *state = state
        .wrapping_mul(0x5851f42d4c957f2d)
        .wrapping_add(0x14057b7ef767814f);
    let mut e = (*state >> 32) as u32;
    e ^= e >> 11;
    e ^= (e << 7) & 0x9d2c5680;
    e ^= (e << 15) & 0xefc60000;
    e ^= e >> 18;
    f64::from((e & 0x7fffffff) as i32) * 4.656612873077393e-10
}

/// `0x351160` with one thread: `y = A D A' v` over the structural columns and the slack of every row other than '=' (`D[nc + i] v_i`, then the
/// columns in order: `t_j = (sum_k a_kj v_k) D_j`, `y_k += a_kj t_j`; a uniform column (`Q+0x648`) sums `v_k s_j` without the alignment peel and adds `s_j t_j`).
fn matvec(f: &StdForm, d: &[f64], v: &[f64]) -> Vec<f64> {
    let (m, nc) = (f.m, f.n2);
    let mut y: Vec<f64> = (0..m)
        .map(|i| {
            if f.sense[i] == b'=' {
                0.0
            } else {
                v[i] * d[nc + i]
            }
        })
        .collect();
    for j in 0..nc {
        let (s, e) = (f.col_start[j], f.col_start[j + 1]);
        let cnt = e - s;
        if cnt == 0 {
            continue;
        }
        let sj = f.scale.as_ref().map_or(0.0, |u| u[j]);
        if sj != 0.0 {
            let t = dot8(cnt, 0, &|q| v[f.row_idx[s + q]] * sj) * d[j];
            let u = sj * t;
            for q in s..e {
                let r = f.row_idx[q];
                y[r] = u + y[r];
            }
        } else {
            let t = dot8(cnt, s & 1, &|q| f.val[s + q] * v[f.row_idx[s + q]]) * d[j];
            for q in s..e {
                let r = f.row_idx[q];
                y[r] = f.val[q] * t + y[r];
            }
        }
    }
    y
}

/// The result of one solve.
pub struct Solve {
    pub dx: Vec<f64>,
    pub dy: Vec<f64>,
    /// `[rsp+0x38]` of the reference at the end: the index of the last step started (added to the work counter `Q+0x290`).
    pub steps: i64,
    /// The stopping tolerance.
    pub tol: f64,
}

/// `0x34a580` (Threads = 1, `Q+0x28 = 0`). `it` is `Q+0x40` (-1 for the starting point: the tolerance comes from `|b|`), `q190` is `Q+0x190`
/// (the primal residual norm of the iterate), `rng` the generator state.
pub fn solve(
    a: &Amat,
    f: &StdForm,
    d: &[f64],
    rx: &[f64],
    ry: &[f64],
    it: i32,
    q190: f64,
    rng: &mut u64,
) -> Solve {
    let (m, nc, n) = (f.m, f.n2, f.n());
    // `w = rx D`; `b = A w + ry`.
    let w: Vec<f64> = (0..n).map(|i| rx[i] * d[i]).collect();
    let t = a.aw(f, &w);
    let b: Vec<f64> = (0..m).map(|i| t[i] + ry[i]).collect();
    // The start vector.
    let mut x: Vec<f64> = (0..m).map(|_| draw(rng) * 1e-8).collect();
    // The preconditioner: the diagonal of `A D A'`.
    let mut inv = vec![0.0; m];
    for j in 0..nc {
        for q in f.col_start[j]..f.col_start[j + 1] {
            let (r, av) = (f.row_idx[q], f.val[q]);
            inv[r] = (d[j] * av) * av + inv[r];
        }
    }
    for i in 0..m {
        if f.sense[i] != b'=' {
            inv[i] = inv[i] + d[nc + i];
        }
        inv[i] = 1.0 / inv[i];
    }
    // The tolerance (`0x34b60c`; at the starting point `0x34c6c9`).
    let tol = if it == -1 {
        let nb = dot_lanes(&b, &b, 0).sqrt() / 100.0;
        if nb > 1.0 {
            nb
        } else {
            1.0
        }
    } else {
        let t = q190 / 10.0;
        if t > 1e-8 {
            t
        } else {
            1e-8
        }
    };
    let ax = matvec(f, d, &x);
    let mut r: Vec<f64> = (0..m).map(|i| b[i] - ax[i]).collect();
    let mut z: Vec<f64> = (0..m).map(|i| r[i] * inv[i]).collect();
    let mut az = matvec(f, d, &z);
    let mut num = 0.0 + dot_lanes(&az, &z, 0);
    let mut p = z.clone();
    let mut ap = az.clone();
    let cap = 100 * m as i64;
    let mut k: i64 = 0;
    let mut last = 0i64;
    while k < cap {
        last = k;
        let den = 0.0 + dot8(m, 0, &|i| (ap[i] * ap[i]) * inv[i]);
        if den == 0.0 {
            break;
        }
        let alpha = num / den;
        for i in 0..m {
            x[i] = p[i] * alpha + x[i];
        }
        let na = -alpha;
        for i in 0..m {
            r[i] = ap[i] * na + r[i];
        }
        let rn = (0.0 + dot_lanes(&r, &r, 0)).sqrt();
        if k > 2 && tol > rn {
            break;
        }
        for i in 0..m {
            z[i] = r[i] * inv[i];
        }
        az = matvec(f, d, &z);
        let num2 = 0.0 + dot_lanes(&az, &z, 0);
        if num2 == 0.0 {
            break;
        }
        let beta = num2 / num;
        num = num2;
        for i in 0..m {
            p[i] = p[i] * beta + z[i];
        }
        for i in 0..m {
            ap[i] = ap[i] * beta + az[i];
        }
        k += 1;
    }
    let s = a.aty(f, &x);
    let dx: Vec<f64> = (0..n).map(|j| (s[j] - rx[j]) * d[j]).collect();
    Solve {
        dx,
        dy: x,
        steps: last,
        tol,
    }
}
