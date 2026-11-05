//! Restarted reflected Halpern PDHG, the algorithm of the reference's `Method=6` on the CPU (`Threads=1`, `Crossover=0`).
//!
//! Every rule here is taken from `reverse-engineer/specs/pdhg*.md` (setup, iteration, residuals and termination, restart, exit); the
//! sums follow the reference's accumulation orders (`blocksum`, `gather_sum`), which assume 16-byte aligned arrays. Arrays are
//! scaled by the reference's Ruiz + Pock-Chambolle scaling before the iteration and unscaled on return.
//! Diagnostic records: `ROSE_DUMP=pdhg` (`crate::dump`).

mod ordered;
use crate::dump;
use crate::model::{LpModel, Sense};
use crate::simplex::Status;
use ordered::Ordered;

/// The reference's infinity inside the scaled model.
const BIG: f64 = 1e30;
/// Interval of the residual, termination and restart tests (hard coded in the reference).
const CHECK: u64 = 100;

#[derive(Clone, Debug)]
pub struct PdhgOptions {
    /// `PDHGIterLimit`; the count includes the power iterations done before the loop.
    pub iter_limit: f64,
    /// `PDHGRelTol`, `PDHGAbsTol`, `PDHGConvTol`.
    pub rel_tol: f64,
    pub abs_tol: f64,
    pub conv_tol: f64,
    /// Keep one record per residual test in `PdhgResult::checks` (the fields of the `ROSE_DUMP=pdhg` record).
    pub record_checks: bool,
    /// The `Seed` parameter: any nonzero value starts the primal iterate at `clamp(1e-6 * u_j, lb, ub)` with `u_j` the first draws of the generator (the power
    /// iteration's start vector follows in the same stream); the value itself does not matter (reference: Seed 1 and 7 give the same run).
    pub seed: i32,
    /// Differential-test hook (`reverse-engineer/replay/gdb_pdhg_checks.py`, `INJECT`): at the residual test with `n8 == at`, set the merit R (`Merit`: restart
    /// rule 5, the R == 0 power-iteration budget) or the relative dual residual (`RelDual`: ratio 0, the primal-weight fallback) to 0 before the rule runs.
    pub inject: Option<(u64, Inject)>,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Inject {
    Merit,
    RelDual,
    /// The primal objective set to 1e101 before the termination test: status 12.
    Pobj,
}

impl Default for PdhgOptions {
    fn default() -> Self {
        PdhgOptions {
            iter_limit: f64::INFINITY,
            rel_tol: 1e-6,
            abs_tol: 1e-6,
            conv_tol: 1e-6,
            record_checks: false,
            seed: 0,
            inject: None,
        }
    }
}

#[derive(Clone, Debug)]
pub struct PdhgResult {
    pub status: Status,
    /// Objective in the model's own sense, offset included.
    pub objective: f64,
    pub x: Vec<f64>,
    /// Row duals in the sign convention of the original objective sense (`Pi`).
    pub dual: Vec<f64>,
    pub reduced_cost: Vec<f64>,
    /// `PDHGIterCount`: the power iterations of the setup and of the restarts plus the loop iterations.
    pub iterations: u64,
    pub restarts: u32,
    /// Residuals at the end: absolute primal and dual (max norm), relative primal, dual and gap.
    pub primal_residual: f64,
    pub dual_residual: f64,
    pub relative: [f64; 3],
    /// With `PdhgOptions::record_checks`: the state at each residual test after the one at the starting point.
    pub checks: Vec<serde_json::Value>,
}

// ---- sums in the reference's order --------------------------------------------------------------

/// Sum of `n` terms in the order of the compiler's SSE loop (`0x311b80`, the squared-norm loops, the matvec rows): left to right below 8 terms,
/// else `peel` terms (0 or 1) summed first, blocks of 8 in four two-lane accumulators, `(a0+a1)+(a2+a3)`, the two lanes added, the tail left to right.
fn blocksum(n: usize, peel: usize, term: impl Fn(usize) -> f64) -> f64 {
    if n < 8 {
        let mut s = 0.0;
        for i in 0..n {
            s += term(i);
        }
        return s;
    }
    let mut s = 0.0;
    for i in 0..peel {
        s += term(i);
    }
    let main = n - ((n - peel) & 7);
    let mut a = [[s, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0]];
    let mut k = peel;
    while k < main {
        for q in 0..4 {
            for l in 0..2 {
                a[q][l] += term(k + 2 * q + l);
            }
        }
        k += 8;
    }
    let t0 = a[0][0] + a[1][0];
    let t1 = a[0][1] + a[1][1];
    let u0 = a[2][0] + a[3][0];
    let u1 = a[2][1] + a[3][1];
    let mut r = (t0 + u0) + (t1 + u1);
    for i in main..n {
        r += term(i);
    }
    r
}

fn dot(a: &[f64], b: &[f64], peel: usize) -> f64 {
    blocksum(a.len(), peel, |i| a[i] * b[i])
}

fn sqdist(a: &[f64], b: &[f64], peel: usize) -> f64 {
    blocksum(a.len(), peel, |i| {
        let d = a[i] - b[i];
        d * d
    })
}

fn seqsum(it: impl Iterator<Item = f64>) -> f64 {
    let mut s = 0.0;
    for v in it {
        s += v;
    }
    s
}

/// `sqrt(sum (a_i * s_i)^2)` as the setup and the residual routine compute their norms.
fn scaled_norm(a: &[f64], s: &[f64]) -> f64 {
    blocksum(a.len(), 0, |i| {
        let t = a[i] * s[i];
        t * t
    })
    .sqrt()
}

/// One row (or column) of a matvec: `vals[b..e] * x[idx[b..e]]` (`0x36ae60`, `0x36b160`). The peel is 1 when `b` is odd (16-byte aligned values).
fn gather_sum(vals: &[f64], idx: &[u32], x: &[f64], b: usize, e: usize) -> f64 {
    let l = e - b;
    let seq = |from: usize, mut s: f64| {
        for q in from..l {
            s += vals[b + q] * x[idx[b + q] as usize];
        }
        s
    };
    if l < 8 {
        return seq(0, 0.0);
    }
    let peel = b & 1;
    if l < peel + 8 {
        return seq(0, 0.0);
    }
    let mut s = 0.0;
    for q in 0..peel {
        s += vals[b + q] * x[idx[b + q] as usize];
    }
    let mut a = [[s, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0]];
    let nb = (l - peel) / 8;
    for t in 0..nb {
        let q0 = peel + 8 * t;
        for k in 0..4 {
            for ln in 0..2 {
                let q = q0 + 2 * k + ln;
                a[k][ln] += x[idx[b + q] as usize] * vals[b + q];
            }
        }
    }
    let t0 = a[0][0] + a[1][0];
    let t1 = a[0][1] + a[1][1];
    let u0 = a[2][0] + a[3][0];
    let u1 = a[2][1] + a[3][1];
    let s = (t0 + u0) + (t1 + u1);
    seq(peel + 8 * nb, s)
}

/// The reference's pseudo random generator (`0x88a8d0`): a 64-bit LCG with the PCG multiplier, tempered like MT19937.
struct Rng(u64);

impl Rng {
    fn next(&mut self) -> f64 {
        self.0 = self
            .0
            .wrapping_mul(0x5851f42d4c957f2d)
            .wrapping_add(0x14057b7ef767814f);
        let mut y = (self.0 >> 32) as u32;
        y ^= y >> 11;
        y ^= (y << 7) & 0x9d2c5680;
        y ^= (y << 15) & 0xefc60000;
        y ^= y >> 18;
        f64::from(y & 0x7fffffff) * 4.656612873077393e-10
    }
}

// ---- the model the reference iterates on --------------------------------------------------------

struct Col {
    ent: Vec<(usize, f64)>,
    c: f64,
    lb: f64,
    ub: f64,
}

/// The internal LP before scaling: `min c'x`, rows `A x (=,<=) b`, `lb <= x <= ub` (infinite bounds are +-infinity here).
struct Internal {
    m: usize,
    cols: Vec<Col>,
    /// Row types: true for an equality; every other row is `A x <= b` (a G row is stored negated).
    eq: Vec<bool>,
    neg: Vec<bool>,
    b: Vec<f64>,
    /// Internal row of each original row (`None`: a free row, dropped).
    row_map: Vec<Option<usize>>,
    flip: Vec<bool>,
    norig: usize,
    sense: f64,
    offset: f64,
}

fn finite_lo(v: f64) -> f64 {
    if v <= -1e20 {
        f64::NEG_INFINITY
    } else {
        v
    }
}

fn finite_up(v: f64) -> f64 {
    if v >= 1e20 {
        f64::INFINITY
    } else {
        v
    }
}

/// A row side that is not the file's right-hand side: infinite from the model's own infinity (1e30).
fn row_side(v: f64) -> f64 {
    if v.abs() >= crate::model::INF {
        v.signum() * f64::INFINITY
    } else {
        v
    }
}

/// The sides of row `i` as the reference uses them. The right-hand side is used as written, even 1e20 or 1e30 or 1e100 (it enters ||b|| and the
/// residuals; the reference never treats a row side as infinite); the other side of an L or G row is infinite unless a range makes it finite.
fn row_sides(model: &LpModel, i: usize) -> (f64, f64) {
    let ranged = model.row_range.get(i).is_some_and(|&r| r > 0.0);
    // Gurobi's own infinity: a side beyond 1e100 is stored as +-1e100.
    let c = |v: f64| {
        if v.is_finite() {
            v.clamp(-1e100, 1e100)
        } else {
            v
        }
    };
    let (lo, up) = match model.row_kind.get(i) {
        Some(&b'L') => (
            if ranged {
                model.row_lo[i]
            } else {
                row_side(model.row_lo[i])
            },
            model.row_up[i],
        ),
        Some(&b'G') => (
            model.row_lo[i],
            if ranged {
                model.row_up[i]
            } else {
                row_side(model.row_up[i])
            },
        ),
        Some(_) => (model.row_lo[i], model.row_up[i]),
        None => (row_side(model.row_lo[i]), row_side(model.row_up[i])),
    };
    (c(lo), c(up))
}

fn build(model: &LpModel) -> Internal {
    let n0 = model.num_cols();
    let m0 = model.num_rows();
    let sense = if model.sense == Sense::Maximize {
        -1.0
    } else {
        1.0
    };
    let mut row_map = vec![None; m0];
    let mut m = 0;
    let mut eq = Vec::new();
    let mut neg = Vec::new();
    let mut b = Vec::new();
    // (internal row, sign, width) of the range columns, by original row
    let mut slack: Vec<Option<(usize, f64, f64)>> = vec![None; m0];
    for i in 0..m0 {
        let (lo, up) = row_sides(model, i);
        if lo == f64::NEG_INFINITY && up == f64::INFINITY {
            continue;
        }
        row_map[i] = Some(m);
        let g = model.row_kind.get(i) == Some(&b'G');
        if lo == up {
            eq.push(true);
            neg.push(false);
            b.push(lo);
        } else if lo > f64::NEG_INFINITY && up < f64::INFINITY {
            // ranged: an equality row plus a range column [0, |R|]
            let w = match model.row_range.get(i) {
                Some(&r) if r > 0.0 => r,
                _ => up - lo,
            };
            eq.push(true);
            neg.push(false);
            if g {
                b.push(lo);
                slack[i] = Some((m, -1.0, w));
            } else {
                b.push(up);
                slack[i] = Some((m, 1.0, w));
            }
        } else if up < f64::INFINITY {
            eq.push(false);
            neg.push(false);
            b.push(up);
        } else {
            eq.push(false);
            neg.push(true);
            b.push(-lo);
        }
        m += 1;
    }
    let mut cols = Vec::with_capacity(n0);
    let mut flip = vec![false; n0];
    for j in 0..n0 {
        let (mut lb, mut ub) = (finite_lo(model.col_lo[j]), finite_up(model.col_up[j]));
        let mut c = sense * model.obj[j];
        let mut ent: Vec<(usize, f64)> = Vec::new();
        for (i, v) in model.col(j) {
            if let Some(r) = row_map[i] {
                if !ent.iter().any(|e| e.0 == r) {
                    ent.push((r, if neg[r] { -v } else { v }));
                }
            }
        }
        if lb == f64::NEG_INFINITY && ub != f64::INFINITY {
            flip[j] = true;
            for e in ent.iter_mut() {
                e.1 = -e.1;
            }
            c = -c;
            (lb, ub) = (-ub, f64::INFINITY);
        }
        cols.push(Col { ent, c, lb, ub });
    }
    let order: Vec<usize> = if model
        .range_order
        .iter()
        .filter(|&&i| slack[i].is_some())
        .count()
        == slack.iter().filter(|s| s.is_some()).count()
    {
        model
            .range_order
            .iter()
            .copied()
            .filter(|&i| slack[i].is_some())
            .collect()
    } else {
        (0..m0).filter(|&i| slack[i].is_some()).collect()
    };
    for i in order {
        let (r, sg, w) = slack[i].unwrap();
        cols.push(Col {
            ent: vec![(r, sg)],
            c: 0.0,
            lb: 0.0,
            ub: finite_up(w),
        });
    }
    Internal {
        m,
        cols,
        eq,
        neg,
        b,
        row_map,
        flip,
        norig: n0,
        sense,
        offset: sense * model.offset,
    }
}

/// Ruiz x10 (simultaneous) and one Pock-Chambolle pass; returns the row and column scale factors.
fn scale_factors(p: &Internal) -> (Vec<f64>, Vec<f64>) {
    let (m, n) = (p.m, p.cols.len());
    let inv_sqrt = |x: f64| if x == 0.0 { 1.0 } else { 1.0 / x.sqrt() };
    let mut r = vec![1.0; m];
    let mut c = vec![1.0; n];
    for _ in 0..10 {
        let mut rmax = vec![0.0f64; m];
        let mut cmax = vec![0.0f64; n];
        for (j, col) in p.cols.iter().enumerate() {
            for &(i, v) in &col.ent {
                let val = (r[i] * c[j]) * v.abs();
                if val > rmax[i] {
                    rmax[i] = val;
                }
                if val > cmax[j] {
                    cmax[j] = val;
                }
            }
        }
        for j in 0..n {
            c[j] *= inv_sqrt(cmax[j]);
        }
        for i in 0..m {
            r[i] *= inv_sqrt(rmax[i]);
        }
    }
    let mut rs = vec![0.0f64; m];
    let mut cs = vec![0.0f64; n];
    for (j, col) in p.cols.iter().enumerate() {
        for &(i, v) in &col.ent {
            let val = (c[j] * r[i]) * v.abs();
            rs[i] += val;
            cs[j] += val;
        }
    }
    for j in 0..n {
        c[j] *= inv_sqrt(cs[j]);
    }
    for i in 0..m {
        r[i] *= inv_sqrt(rs[i]);
    }
    (r, c)
}

/// The objective scale: a power of two (`0x11c8e0`).
fn objective_factor(c: &[f64]) -> f64 {
    let nz: Vec<f64> = c.iter().map(|v| v.abs()).filter(|&v| v != 0.0).collect();
    if nz.is_empty() {
        return 1.0;
    }
    let mx = nz.iter().cloned().fold(f64::MIN, f64::max);
    let mn = nz.iter().cloned().fold(f64::MAX, f64::min);
    let x = mn.max(1e-6) * mx;
    let x4 = (100.0 / mx).min(1.0 / (10.0 * x.sqrt()));
    if 2.0 > x4 {
        return 1.0;
    }
    let v = if x4 > 1000.0 { 1000.0 } else { x4 };
    crate::real::pow2(v.log2().floor() as i32)
}

/// The scaled model the iteration runs on.
struct Scaled {
    m: usize,
    n: usize,
    rp: Vec<usize>,
    rcol: Vec<u32>,
    rval: Vec<f64>,
    cp: Vec<usize>,
    crow: Vec<u32>,
    cval: Vec<f64>,
    c: Vec<f64>,
    b: Vec<f64>,
    lb: Vec<f64>,
    ub: Vec<f64>,
    eq: Vec<bool>,
    /// Column and row scale factors, and their inverses (`state+0x158`, `+0x160`).
    cs: Vec<f64>,
    rs: Vec<f64>,
    cinv: Vec<f64>,
    rinv: Vec<f64>,
    objf: f64,
}

impl Scaled {
    fn new(p: &Internal) -> Scaled {
        let (m, n) = (p.m, p.cols.len());
        let (r, c) = scale_factors(p);
        let mut cp = vec![0usize];
        let mut crow = Vec::new();
        let mut cval = Vec::new();
        let mut rows: Vec<Vec<(u32, f64)>> = vec![Vec::new(); m];
        for (j, col) in p.cols.iter().enumerate() {
            for &(i, v) in &col.ent {
                let sv = (r[i] * c[j]) * v;
                crow.push(i as u32);
                cval.push(sv);
                rows[i].push((j as u32, sv));
            }
            cp.push(crow.len());
        }
        let mut rp = vec![0usize];
        let mut rcol = Vec::new();
        let mut rval = Vec::new();
        for row in &rows {
            for &(j, v) in row {
                rcol.push(j);
                rval.push(v);
            }
            rp.push(rcol.len());
        }
        let mut cvec: Vec<f64> = (0..n).map(|j| c[j] * p.cols[j].c).collect();
        let f = objective_factor(&cvec);
        if f != 1.0 {
            for v in cvec.iter_mut() {
                *v *= f;
            }
        }
        let b: Vec<f64> = (0..m)
            .map(|i| {
                if p.b[i].abs() >= 1e30 {
                    p.b[i]
                } else {
                    r[i] * p.b[i]
                }
            })
            .collect();
        let lb = (0..n)
            .map(|j| {
                if p.cols[j].lb > f64::NEG_INFINITY {
                    p.cols[j].lb / c[j]
                } else {
                    -BIG
                }
            })
            .collect();
        let ub = (0..n)
            .map(|j| {
                if p.cols[j].ub < f64::INFINITY {
                    p.cols[j].ub / c[j]
                } else {
                    BIG
                }
            })
            .collect();
        let cinv = c.iter().map(|v| 1.0 / v).collect();
        let rinv = r.iter().map(|v| 1.0 / v).collect();
        Scaled {
            m,
            n,
            rp,
            rcol,
            rval,
            cp,
            crow,
            cval,
            c: cvec,
            b,
            lb,
            ub,
            eq: p.eq.clone(),
            cs: c,
            rs: r,
            cinv,
            rinv,
            objf: f,
        }
    }

    /// `out = A x`.
    fn mul(&self, x: &[f64], out: &mut [f64]) {
        for i in 0..self.m {
            out[i] = gather_sum(&self.rval, &self.rcol, x, self.rp[i], self.rp[i + 1]);
        }
    }

    /// `out = A' y`.
    fn mul_t(&self, y: &[f64], out: &mut [f64]) {
        for j in 0..self.n {
            out[j] = gather_sum(&self.cval, &self.crow, y, self.cp[j], self.cp[j + 1]);
        }
    }
}

// ---- the iteration ------------------------------------------------------------------------------

/// Scalars the residual routine (`0x3685a0`) writes.
#[derive(Clone, Copy, Default, Debug)]
struct Resid {
    pobj: f64,
    dobj: f64,
    compl: f64,
    p_inf: f64,
    d_inf: f64,
    p2: f64,
    d2: f64,
    rel_p: f64,
    rel_d: f64,
    gap: f64,
}

#[cfg_attr(test, derive(Clone))]
struct State<'a> {
    a: &'a Scaled,
    objoff: f64,
    x: Vec<f64>,
    y: Vec<f64>,
    x0: Vec<f64>,
    y0: Vec<f64>,
    xh: Vec<f64>,
    yh: Vec<f64>,
    zc: Vec<f64>,
    xb: Vec<f64>,
    yb: Vec<f64>,
    ax: Vec<f64>,
    aty: Vec<f64>,
    eta: f64,
    omega: f64,
    omega0: f64,
    integral: f64,
    /// `state[0]`, `[8]`, `[0x10]`, `[0x18]`.
    total: u64,
    n8: u64,
    k: u64,
    power_iters: u64,
    r: f64,
    r_prev: f64,
    r0: f64,
    nb: f64,
    nc: f64,
    ncons: usize,
    res: Resid,
    /// The last power iteration vector (`state+0x180`), its residual (`+0x188`) and tolerance (`+0x190`).
    pv: Vec<f64>,
    perr: f64,
    ptol: f64,
}

// Matching slices establish all bounds before entering either dense kernel.
// Keeping these passes separate from sparse gathering also gives the compiler
// non-aliasing input/output references and independent coordinates to vectorize.
#[inline(never)]
fn primal_dense(
    x: &mut [f64],
    xb: &mut [f64],
    anchor: &[f64],
    aty: &[f64],
    a: &Scaled,
    tau: f64,
    w: f64,
) {
    let n = x.len();
    let xb = &mut xb[..n];
    let anchor = &anchor[..n];
    let aty = &aty[..n];
    let c = &a.c[..n];
    let lb = &a.lb[..n];
    let ub = &a.ub[..n];
    let w1 = 1.0 - w;
    for j in 0..n {
        let old = x[j];
        let g = (c[j] - aty[j]) * tau;
        let t = old - g;
        let xp = t.max(lb[j]).min(ub[j]);
        let reflected = (xp + xp) - old;
        xb[j] = reflected;
        x[j] = reflected * w + anchor[j] * w1;
    }
}
#[inline(never)]
fn dual_dense(y: &mut [f64], anchor: &[f64], ax: &[f64], a: &Scaled, sigma: f64, w: f64) {
    let n = y.len();
    let anchor = &anchor[..n];
    let ax = &ax[..n];
    let rhs = &a.b[..n];
    let eq = &a.eq[..n];
    let w1 = 1.0 - w;
    for i in 0..n {
        let old = y[i];
        let t3 = old - ax[i] * sigma;
        let nsb = (-sigma) * rhs[i];
        let yp = if eq[i] { t3 - nsb } else { t3 - t3.max(nsb) };
        let reflected = (yp - old) + yp;
        y[i] = reflected * w + anchor[i] * w1;
    }
}

impl<'a> State<'a> {
    /// Outside residual, restart and termination boundaries, only the reflected
    /// primal vector and the next x/y feed subsequent work. Keep both sparse
    /// matvecs intact and evaluate the dense state transforms in independent lanes.
    fn ordinary_step(&mut self, tau: f64, sigma: f64) {
        let a = self.a;
        let w = (self.k + 1) as f64 / (self.k + 2) as f64;
        a.mul_t(&self.y, &mut self.aty);
        primal_dense(&mut self.x, &mut self.xb, &self.x0, &self.aty, a, tau, w);
        a.mul(&self.xb, &mut self.ax);
        dual_dense(&mut self.y, &self.y0, &self.ax, a, sigma, w);
    }

    fn ordinary_ordered_step(&mut self, tau: f64, sigma: f64, ordered: &mut Ordered) {
        let w = (self.k + 1) as f64 / (self.k + 2) as f64;
        if ordered.step(self.a, &self.y, &self.x, &self.x0, &mut self.ax, tau, w) {
            dual_dense(&mut self.y, &self.y0, &self.ax, self.a, sigma, w);
        } else {
            self.ordinary_step(tau, sigma);
        }
    }
    /// `0x3685a0`: the residual scalars at `(xh, yh)`. `flag1`: the first pass and the final pass.
    fn residual(&mut self, flag1: bool) -> Resid {
        let a = self.a;
        let (m, n) = (a.m, a.n);
        a.mul(&self.xh, &mut self.ax);
        a.mul_t(&self.yh, &mut self.aty);
        let mut o = Resid::default();
        let v: Vec<f64> = (0..m)
            .map(|i| {
                let s = (a.b[i] - self.ax[i]) * a.rinv[i];
                if a.eq[i] {
                    s.abs().max(0.0)
                } else {
                    (-s).max(0.0)
                }
            })
            .collect();
        o.p_inf = v.iter().cloned().fold(0.0, f64::max);
        o.p2 = seqsum(v.iter().map(|t| t * t)).sqrt();
        let zz: Vec<f64>;
        if !flag1 {
            let t: Vec<f64> = (0..n)
                .map(|j| (a.c[j] - self.aty[j] - self.zc[j]) * a.cinv[j])
                .collect();
            o.d_inf = t.iter().map(|v| v.abs()).fold(0.0, f64::max);
            o.d2 = seqsum(t.iter().map(|v| v * v)).sqrt();
            zz = self.zc.clone();
        } else {
            let mut inf = f64::MIN;
            let mut ss = 0.0;
            for j in 0..n {
                let t = (a.c[j] - self.aty[j]) * a.cinv[j];
                let v = if a.lb[j] > -BIG {
                    if a.ub[j] >= BIG {
                        (-t).max(0.0)
                    } else {
                        0.0
                    }
                } else if a.ub[j] < BIG {
                    t.max(0.0)
                } else {
                    t.max(-t)
                };
                inf = inf.max(v);
                ss += v * v;
            }
            o.d_inf = inf;
            o.d2 = ss.sqrt();
            zz = (0..n).map(|j| a.c[j] - self.aty[j]).collect();
        }
        o.pobj = dot(&self.xh, &a.c, 0) / a.objf + self.objoff;
        let mut acc = 0.0;
        for i in 0..m {
            acc += self.yh[i] * a.b[i];
        }
        for j in 0..n {
            if zz[j] > 0.0 && a.lb[j] > -BIG {
                acc += zz[j] * a.lb[j];
            } else if zz[j] < 0.0 && a.ub[j] < BIG {
                acc += zz[j] * a.ub[j];
            }
        }
        o.dobj = acc / a.objf + self.objoff;
        let mut acc = 0.0;
        for j in 0..n {
            let r = a.c[j] - self.aty[j];
            if a.lb[j] > -BIG && r > 1e-10 {
                acc += (self.xh[j] - a.lb[j]) * r;
            } else if a.ub[j] < BIG && r < -1e-10 {
                acc += (self.xh[j] - a.ub[j]) * r;
            }
        }
        for i in 0..m {
            let sl = a.b[i] - self.ax[i];
            if !a.eq[i] && sl > 1e-10 {
                acc += -sl * self.yh[i];
            }
        }
        o.compl = if self.ncons > 0 {
            (acc / self.ncons as f64).max(0.0)
        } else {
            acc.max(0.0)
        };
        o.rel_p = o.p2 / (self.nb + 1.0);
        o.rel_d = o.d2 / (self.nc + 1.0);
        o.gap = (o.pobj - o.dobj).abs() / ((1.0 + o.pobj.abs()) + o.dobj.abs());
        o
    }

    /// `0x36a7d0` with the previous status 1: 2 optimal, 7 iteration limit, 12 numerics, else 1.
    fn status(&self, o: &Resid, opts: &PdhgOptions) -> u8 {
        let (a, r, c) = (opts.abs_tol, opts.rel_tol, opts.conv_tol);
        if a > o.p_inf && a > o.d_inf && c > o.gap {
            return 2;
        }
        if r > 0.0 && r > o.rel_p && r > o.rel_d && c > o.gap {
            return 2;
        }
        if self.total as f64 >= opts.iter_limit {
            return 7;
        }
        if o.pobj.abs() > 1e100 || o.dobj.abs() > 1e100 || o.compl.abs() > 1e100 {
            return 12;
        }
        1
    }

    /// `0x364d50`.
    fn restart_rule(&self) -> u8 {
        let n = self.n8;
        if n < CHECK {
            0
        } else if n == CHECK {
            1
        } else if self.r == 0.0 {
            5
        } else if 0.2 * self.r0 >= self.r {
            2
        } else if 0.5 * self.r0 >= self.r && self.r > self.r_prev {
            3
        } else if n as f64 * 0.36 <= self.k as f64 {
            4
        } else {
            0
        }
    }

    /// `0x3678d0`: the merit R at the current point and its reflection.
    fn merit(&mut self) -> f64 {
        let a = self.a;
        let dx: Vec<f64> = (0..a.n).map(|j| self.x[j] - self.xb[j]).collect();
        let dy: Vec<f64> = (0..a.m).map(|i| self.y[i] - self.yb[i]).collect();
        a.mul_t(&dy, &mut self.aty);
        let ddx = dot(&dx, &dx, 0);
        let ddy = dot(&dy, &dy, a.n % 2);
        let cr = dot(&dx, &self.aty, 0);
        let t = cr * 2.0 * self.eta;
        let s = ddx * self.omega + ddy / self.omega;
        let r2 = t + s;
        // `0x368395`: a negative (or NaN) radicand is not rooted and R is stored as 0.
        if r2 >= 0.0 {
            r2.sqrt()
        } else {
            0.0
        }
    }

    /// `0x364df0`: restart to the step output, with the primal weight update.
    fn apply_restart(&mut self) {
        let dxn = sqdist(&self.xh, &self.x0, 0).sqrt();
        let dyn_ = sqdist(&self.yh, &self.y0, 0).sqrt();
        let ratio = if self.res.rel_p != 0.0 {
            self.res.rel_d / self.res.rel_p
        } else {
            1e100
        };
        if dxn > 1e-16 && dyn_ > 1e-16 && 1e12 > dxn && 1e12 > dyn_ && ratio > 1e-8 && 1e8 > ratio {
            let d = (dyn_.ln() - dxn.ln()) - self.omega.ln();
            let i2 = 0.3 * self.integral + d;
            self.omega *= (0.99 * d + 0.01 * i2).exp();
            self.integral = i2;
        } else {
            self.omega = self.omega0;
            self.integral = 0.0;
        }
        self.k = 0;
        self.r_prev = 1e100;
        self.x.copy_from_slice(&self.xh);
        self.x0.copy_from_slice(&self.xh);
        self.y.copy_from_slice(&self.yh);
        self.y0.copy_from_slice(&self.yh);
    }

    /// `0x365e60`: the power iteration on `A A'` that sets the step size `eta`. `budget` iterations at most; `first`: the setup call.
    fn power(&mut self, budget: u64, first: bool) {
        let a = self.a;
        let m = a.m;
        let mut s = self.pv.clone();
        let mut w = vec![0.0; a.n];
        let mut u = vec![0.0; m];
        let mut done = 0u64;
        let mut lam = 0.0;
        let mut err = self.perr;
        let mut prev: Option<f64> = None;
        while done < budget {
            if let Some(p) = prev {
                let chunk_end = first && done % 100 == 0;
                if self.ptol > p || (chunk_end && self.ptol >= p) {
                    break;
                }
            }
            let nrm = dot(&s, &s, 0).sqrt();
            let v: Vec<f64> = s.iter().map(|t| t / nrm).collect();
            a.mul_t(&v, &mut w);
            a.mul(&w, &mut u);
            lam = dot(&u, &v, 0);
            let r: Vec<f64> = (0..m).map(|i| u[i] - lam * v[i]).collect();
            err = dot(&r, &r, 0).sqrt();
            s.copy_from_slice(&u);
            prev = Some(err);
            done += 1;
        }
        if done == 0 {
            return;
        }
        let sigma = lam.sqrt();
        self.eta = if err < 0.5 {
            (1.0 - err - err) / sigma
        } else {
            0.01 / sigma
        };
        self.perr = err;
        self.pv = s;
        self.total += done;
        self.power_iters += done;
    }
}

/// A model without rows or without nonzeros never reaches the PDHG driver (the reference answers it without iterating): every column goes to
/// the bound its cost prefers (to the point of smallest magnitude in its bounds when the cost is 0), the duals are 0 and the reduced costs are the costs.
fn trivial(model: &LpModel) -> PdhgResult {
    let sense = if model.sense == Sense::Maximize {
        -1.0
    } else {
        1.0
    };
    let mut status = Status::Optimal;
    let mut x = Vec::with_capacity(model.num_cols());
    for j in 0..model.num_cols() {
        let (lb, ub) = (finite_lo(model.col_lo[j]), finite_up(model.col_up[j]));
        let c = sense * model.obj[j];
        let v = if c > 0.0 {
            lb
        } else if c < 0.0 {
            ub
        } else {
            0f64.max(lb).min(ub)
        };
        if !v.is_finite() {
            status = Status::Unbounded;
        }
        x.push(if v.is_finite() { v } else { 0.0 });
    }
    for i in 0..model.num_rows() {
        let (lo, up) = (row_side(model.row_lo[i]), row_side(model.row_up[i]));
        if model.val.is_empty() && (0.0 < lo || 0.0 > up) {
            status = Status::Infeasible;
        }
    }
    let objective = model.objective(&x);
    PdhgResult {
        status,
        objective,
        reduced_cost: model.obj.clone(),
        dual: vec![0.0; model.num_rows()],
        x,
        iterations: 0,
        restarts: 0,
        primal_residual: 0.0,
        dual_residual: 0.0,
        relative: [0.0; 3],
        checks: Vec::new(),
    }
}

/// Solves `model` with the reflected restarted Halpern PDHG.
pub fn solve(model: &LpModel, opts: &PdhgOptions) -> PdhgResult {
    solve_inner::<true>(model, opts)
}
fn solve_inner<const BATCHED: bool>(model: &LpModel, opts: &PdhgOptions) -> PdhgResult {
    if model.num_rows() == 0 || model.val.is_empty() {
        return trivial(model);
    }
    let p = build(model);
    let a = Scaled::new(&p);
    let (m, n) = (a.m, a.n);
    let nb = scaled_norm(&a.b, &a.rinv);
    let nc = scaled_norm(&a.c, &a.cinv);
    let omega = if nb > 0.0 && nc > 0.0 { nc / nb } else { 1.0 };
    let mut rng = Rng(0);
    let x0: Vec<f64> = (0..n)
        .map(|j| {
            let v = if opts.seed != 0 {
                rng.next() * 1e-6
            } else {
                0.0
            };
            v.max(a.lb[j]).min(a.ub[j])
        })
        .collect();
    let ncons = (0..n).filter(|&j| a.lb[j] > -BIG || a.ub[j] < BIG).count()
        + a.eq.iter().filter(|&&e| !e).count();
    let pv: Vec<f64> = (0..m).map(|_| rng.next()).collect();
    let mut st = State {
        a: &a,
        objoff: p.offset,
        x: x0.clone(),
        y: vec![0.0; m],
        x0: x0.clone(),
        y0: vec![0.0; m],
        xh: x0,
        yh: vec![0.0; m],
        zc: vec![0.0; n],
        xb: vec![0.0; n],
        yb: vec![0.0; m],
        ax: vec![0.0; m],
        aty: vec![0.0; n],
        eta: 0.0,
        omega,
        omega0: omega,
        integral: 0.0,
        total: 0,
        n8: 0,
        k: 0,
        power_iters: 0,
        r: 0.0,
        r_prev: 0.0,
        r0: 0.0,
        nb,
        nc,
        ncons,
        res: Resid::default(),
        pv,
        perr: 1e100,
        ptol: 1e-4,
    };
    let mut restarts = 0u32;
    let mut checks = Vec::new();
    let mut status;
    // iteration 0: the residual at the starting point; the solve ends here if it is already optimal
    let o = st.residual(true);
    st.res = o;
    status = st.status(&o, opts);
    if status == 1 {
        st.r_prev = st.r;
        st.power(250, true);
        let mut ordered = Ordered::new(&a);
        let mut just_restarted = false;
        let mut tau;
        let mut sigma;
        loop {
            let chk = (st.n8 + 1) % CHECK == 0;
            tau = st.eta / st.omega;
            sigma = st.eta * st.omega;
            if BATCHED && !chk && !just_restarted && ((st.total + 1) as f64) < opts.iter_limit {
                st.ordinary_ordered_step(tau, sigma, &mut ordered);
                st.total += 1;
                st.n8 += 1;
                st.k += 1;
                continue;
            }
            ordered.observe(&mut st.x, &mut st.xb);
            // primal step
            a.mul_t(&st.y, &mut st.aty);
            for j in 0..n {
                let g = (a.c[j] - st.aty[j]) * tau;
                let t = st.x[j] - g;
                let xp = t.max(a.lb[j]).min(a.ub[j]);
                if chk {
                    st.xh[j] = xp;
                    st.zc[j] = (xp - t) / tau;
                }
                st.xb[j] = (xp + xp) - st.x[j];
            }
            // dual step
            a.mul(&st.xb, &mut st.ax);
            for i in 0..m {
                let t3 = st.y[i] - st.ax[i] * sigma;
                let nsb = (-sigma) * a.b[i];
                let yp = if a.eq[i] { t3 - nsb } else { t3 - t3.max(nsb) };
                st.yh[i] = yp;
                st.yb[i] = (yp - st.y[i]) + yp;
            }
            if chk || just_restarted {
                st.r = st.merit();
                if just_restarted {
                    st.r0 = st.r;
                }
            }
            // Halpern average
            let w = (st.k + 1) as f64 / (st.k + 2) as f64;
            let w1 = 1.0 - w;
            for j in 0..n {
                st.x[j] = st.xb[j] * w + st.x0[j] * w1;
            }
            for i in 0..m {
                st.y[i] = st.yb[i] * w + st.y0[i] * w1;
            }
            st.total += 1;
            st.n8 += 1;
            st.k += 1;
            just_restarted = false;
            if chk || st.total as f64 >= opts.iter_limit {
                let mut o = st.residual(false);
                if let Some((at, Inject::Pobj)) = opts.inject {
                    if st.n8 == at {
                        o.pobj = 1e101;
                    }
                }
                status = st.status(&o, opts);
                if let (Some((at, what)), 1) = (opts.inject, status) {
                    if st.n8 == at {
                        match what {
                            Inject::Merit => st.r = 0.0,
                            Inject::RelDual => o.rel_d = 0.0,
                            Inject::Pobj => {}
                        }
                    }
                }
                st.res = o;
                if opts.record_checks || dump::on("pdhg") {
                    let rule = if chk && status == 1 {
                        st.restart_rule()
                    } else {
                        0
                    };
                    let rec = serde_json::json!({"total": st.total, "n8": st.n8, "k": st.k, "status": status, "rule": rule, "eta": st.eta, "omega": st.omega,
                        "r": st.r, "r_prev": st.r_prev, "r0": st.r0, "pobj": o.pobj, "dobj": o.dobj, "compl": o.compl,
                        "p_inf": o.p_inf, "d_inf": o.d_inf, "rel_p": o.rel_p, "rel_d": o.rel_d, "gap": o.gap});
                    if dump::on("pdhg") {
                        dump::emit("pdhg", rec.clone());
                    }
                    if opts.record_checks {
                        checks.push(rec);
                    }
                }
                if status != 1 {
                    break;
                }
                if chk {
                    let rule = st.restart_rule();
                    st.r_prev = st.r;
                    if rule != 0 {
                        st.apply_restart();
                        restarts += 1;
                        just_restarted = true;
                        // a continued power iteration while the step size estimate has not converged
                        // (`0x365ee6`: a merit of 0 gives a budget of 1000 and first tightens a converged tolerance, `0x366b7e`)
                        let b = if st.r == 0.0 {
                            if st.ptol >= st.perr {
                                st.ptol = (0.1 * st.ptol).max(1e-16);
                            }
                            1000
                        } else {
                            (0.1 * st.total as f64 - st.power_iters as f64) as i64
                        };
                        {
                            if st.ptol < st.perr && b > 100 {
                                let left = opts.iter_limit - st.total as f64;
                                let b = if b as f64 >= left {
                                    left.max(0.0) as u64
                                } else {
                                    b as u64
                                };
                                st.power(b, false);
                                // the termination test after the chunk (`0x36a7d0`) can only stop on the iteration limit here
                                if st.total as f64 >= opts.iter_limit {
                                    status = 7;
                                    break;
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    // the final residual pass (flag 1) and the read-out
    let o = st.residual(true);
    st.res = o;
    readout(&p, &a, &st, status, restarts, checks)
}

fn readout(
    p: &Internal,
    a: &Scaled,
    st: &State,
    status: u8,
    restarts: u32,
    checks: Vec<serde_json::Value>,
) -> PdhgResult {
    let sense = p.sense;
    let os = a.objf;
    let (m, n) = (a.m, a.n);
    // unscale (0x10b0f0): x = x_s cs, rc = rc_s / (cs os), y = (rs / os) y_s
    let mut x = vec![0.0; p.norig];
    let mut rc = vec![0.0; p.norig];
    for j in 0..p.norig {
        let xj = st.xh[j] * a.cs[j];
        let rcj = (a.c[j] - st.aty[j]) / (a.cs[j] * os);
        let flip = if p.flip[j] { -1.0 } else { 1.0 };
        x[j] = flip * xj;
        rc[j] = sense * flip * rcj;
    }
    let _ = n;
    let mut dual = vec![0.0; p.row_map.len()];
    for (i, r) in p.row_map.iter().enumerate() {
        if let Some(r) = *r {
            let y = (a.rs[r] / os) * st.yh[r];
            let sg = if p.neg[r] { -1.0 } else { 1.0 };
            dual[i] = sg * sense * y;
        }
    }
    let _ = m;
    let st_status = match status {
        2 => Status::Optimal,
        7 => Status::IterationLimit,
        _ => Status::NumericalFailure,
    };
    PdhgResult {
        status: st_status,
        objective: st.res.pobj * sense,
        x,
        dual,
        reduced_cost: rc,
        iterations: st.total,
        restarts,
        primal_residual: st.res.p_inf,
        dual_residual: st.res.d_inf,
        relative: [st.res.rel_p, st.res.rel_d, st.res.gap],
        checks,
    }
}

#[cfg(test)]
mod ordinary_step_tests {
    use super::*;

    // The original complete ordinary iteration, including both staging matvecs
    // and the two separate Halpern passes. Observations remain in solve_inner.
    fn original(st: &mut State<'_>, tau: f64, sigma: f64) {
        let a = st.a;
        a.mul_t(&st.y, &mut st.aty);
        for j in 0..a.n {
            let g = (a.c[j] - st.aty[j]) * tau;
            let t = st.x[j] - g;
            let xp = t.max(a.lb[j]).min(a.ub[j]);
            st.xb[j] = (xp + xp) - st.x[j];
        }
        a.mul(&st.xb, &mut st.ax);
        for i in 0..a.m {
            let t3 = st.y[i] - st.ax[i] * sigma;
            let nsb = (-sigma) * a.b[i];
            let yp = if a.eq[i] { t3 - nsb } else { t3 - t3.max(nsb) };
            st.yh[i] = yp;
            st.yb[i] = (yp - st.y[i]) + yp;
        }
        let w = (st.k + 1) as f64 / (st.k + 2) as f64;
        let w1 = 1.0 - w;
        for j in 0..a.n {
            st.x[j] = st.xb[j] * w + st.x0[j] * w1;
        }
        for i in 0..a.m {
            st.y[i] = st.yb[i] * w + st.y0[i] * w1;
        }
    }

    fn initial(a: &Scaled, mut next: impl FnMut() -> f64) -> State<'_> {
        let mut fill = |n| (0..n).map(|_| next()).collect();
        State {
            a,
            objoff: 0.0,
            x: fill(a.n),
            y: fill(a.m),
            x0: fill(a.n),
            y0: fill(a.m),
            xh: fill(a.n),
            yh: fill(a.m),
            zc: fill(a.n),
            xb: fill(a.n),
            yb: fill(a.m),
            ax: fill(a.m),
            aty: fill(a.n),
            eta: 0.5,
            omega: 2.0,
            omega0: 2.0,
            integral: 0.0,
            total: 143,
            n8: 21,
            k: 12,
            power_iters: 122,
            r: 1.0,
            r_prev: 2.0,
            r0: 3.0,
            nb: 1.0,
            nc: 1.0,
            ncons: a.n,
            res: Resid::default(),
            pv: fill(a.m),
            perr: 0.0,
            ptol: 1e-4,
        }
    }

    fn bits(values: &[f64]) -> Vec<u64> {
        values.iter().map(|x| x.to_bits()).collect()
    }

    #[test]
    fn ordered_step_matches_materialized_vectors() {
        let model = crate::mps::parse(include_str!("../tests/fixtures/afiro.mps")).unwrap();
        let a = Scaled::new(&build(&model));
        let mut rng = 19_u64;
        for case in 0..32 {
            let mut next = || {
                rng = rng.wrapping_mul(6364136223846793005).wrapping_add(1);
                match case % 4 {
                    0 => (rng as i64 as f64) * 1e-16,
                    1 => f64::from_bits(rng),
                    2 => [-0.0, f64::MAX, f64::MIN_POSITIVE, -1.0, 2.0, f64::INFINITY]
                        [(rng >> 32) as usize % 6],
                    _ => [
                        f64::from_bits(0x7ff8_0000_0000_0001),
                        -0.0,
                        -f64::INFINITY,
                        -4.0,
                    ][(rng >> 32) as usize % 4],
                }
            };
            let mut fused = initial(&a, &mut next);
            let mut frozen = fused.clone();
            for _ in 0..3 {
                let (tau, sigma) = (next(), next());
                original(&mut frozen, tau, sigma);
                fused.ordinary_step(tau, sigma);
                assert_eq!(bits(&frozen.x), bits(&fused.x), "x case {case}");
                assert_eq!(bits(&frozen.y), bits(&fused.y), "y case {case}");
                assert_eq!(bits(&frozen.xb), bits(&fused.xb), "xb case {case}");
                frozen.k += 1;
                fused.k += 1;
            }
        }
    }

    #[test]
    fn complete_results_around_checks_and_restarts() {
        let m = crate::mps::parse(include_str!("../tests/fixtures/adlittle.mps")).unwrap();
        let mut opts = PdhgOptions {
            iter_limit: 4000.0,
            abs_tol: 0.0,
            rel_tol: 0.0,
            conv_tol: 0.0,
            record_checks: true,
            ..Default::default()
        };
        let first = solve_inner::<false>(&m, &opts);
        let at = first.checks[0]["total"].as_u64().unwrap();
        for limit in [1, at - 1, at, at + 1, at + 99, at + 100, at + 101, at + 850] {
            for inject in [None, Some((100, Inject::Merit)), Some((100, Inject::Pobj))] {
                opts.iter_limit = limit as f64;
                opts.inject = inject;
                let a = solve_inner::<false>(&m, &opts);
                let b = solve_inner::<true>(&m, &opts);
                assert_eq!(
                    (a.status, a.iterations, a.restarts),
                    (b.status, b.iterations, b.restarts)
                );
                assert_eq!(a.checks, b.checks);
                assert_eq!(a.objective.to_bits(), b.objective.to_bits());
                assert_eq!(bits(&a.x), bits(&b.x));
                assert_eq!(bits(&a.dual), bits(&b.dual));
                assert_eq!(bits(&a.reduced_cost), bits(&b.reduced_cost));
                assert_eq!(bits(&a.relative), bits(&b.relative));
                assert_eq!(a.primal_residual.to_bits(), b.primal_residual.to_bits());
                assert_eq!(a.dual_residual.to_bits(), b.dual_residual.to_bits());
            }
        }
    }
    #[test]
    fn ordered_recurrence_matches_full_observations_and_restart() {
        let model = crate::mps::parse(include_str!("../tests/fixtures/afiro.mps")).unwrap();
        let mut a = Scaled::new(&build(&model));
        for j in 0..a.n {
            a.c[j] = if j % 2 == 0 { 100. } else { -100. };
            a.lb[j] = -0.0;
            a.ub[j] = 1.0;
        }
        let (tau, sigma) = (0.1, 0.001);
        for inject in [
            0.0,
            -0.0,
            f64::INFINITY,
            f64::NEG_INFINITY,
            f64::from_bits(0x7ff8_0000_0000_1001),
            f64::from_bits(0xfff8_0000_0000_0042),
        ] {
            let mut frozen = initial(&a, || 0.0);
            let mut current = frozen.clone();
            let mut ordered = Ordered::new(&a);
            for step in 0..83 {
                if step == 38 {
                    current.y[0] = inject;
                    frozen.y[0] = inject;
                }
                if step == 48 || step == 82 {
                    ordered.observe(&mut current.x, &mut current.xb);
                    for st in [&mut frozen, &mut current] {
                        a.mul_t(&st.y, &mut st.aty);
                        for j in 0..a.n {
                            let t = st.x[j] - (a.c[j] - st.aty[j]) * tau;
                            st.xh[j] = t.max(a.lb[j]).min(a.ub[j]);
                            st.zc[j] = (st.xh[j] - t) / tau;
                        }
                        original(st, tau, sigma);
                        st.r = st.merit();
                        st.res = st.residual(false);
                    }
                    assert_eq!(
                        bits(&frozen.aty),
                        bits(&current.aty),
                        "raw observed transpose {step}"
                    );
                    assert_eq!(bits(&frozen.xh), bits(&current.xh));
                    assert_eq!(bits(&frozen.zc), bits(&current.zc));
                    assert_eq!(frozen.res.pobj.to_bits(), current.res.pobj.to_bits());
                    assert_eq!(frozen.res.dobj.to_bits(), current.res.dobj.to_bits());
                    assert_eq!(frozen.r.to_bits(), current.r.to_bits());
                    frozen.apply_restart();
                    current.apply_restart();
                } else {
                    original(&mut frozen, tau, sigma);
                    current.ordinary_ordered_step(tau, sigma, &mut ordered);
                }
                let (cx, cxb) = ordered.snapshot(&current.x, &current.xb);
                assert_eq!(
                    bits(&frozen.xb),
                    bits(&cxb),
                    "xb {step}, injection {:x}",
                    inject.to_bits()
                );
                assert_eq!(bits(&frozen.x), bits(&cx));
                assert_eq!(bits(&frozen.y), bits(&current.y));
                frozen.k += 1;
                current.k += 1;
            }
        }
    }
    #[test]
    fn randomized_permuted_state_matches_all_original_outputs() {
        let model = crate::mps::parse(include_str!("../tests/fixtures/greenbea.mps")).unwrap();
        let a = Scaled::new(&build(&model));
        let mut rng = 2613921u64;
        for case in 0..10 {
            let mut next = || {
                rng ^= rng << 13;
                rng ^= rng >> 7;
                rng ^= rng << 17;
                if case % 3 == 0 {
                    f64::from_bits(rng)
                } else if case % 3 == 1 {
                    [
                        -0.0,
                        0.0,
                        f64::from_bits(1),
                        f64::from_bits(0x7ff8_0000_0000_0132),
                        f64::INFINITY,
                        -1.0,
                        3.0,
                    ][rng as usize % 7]
                } else {
                    (rng as i64 as f64) * 1e-18
                }
            };
            let mut old = initial(&a, &mut next);
            let mut curr = old.clone();
            let mut ord = Ordered::new(&a);
            for k in 0..7 {
                let tau = next();
                let sigma = next();
                original(&mut old, tau, sigma);
                curr.ordinary_ordered_step(tau, sigma, &mut ord);
                let (x, xb) = ord.snapshot(&curr.x, &curr.xb);
                for (v, actual) in [(&old.x, &x), (&old.xb, &xb), (&old.y, &curr.y)] {
                    for (j, (a, b)) in v.iter().zip(actual).enumerate() {
                        assert_eq!(a.to_bits(), b.to_bits(), "case {case} it {k} pos {j}");
                    }
                }
                if k == 3 {
                    ord.observe(&mut curr.x, &mut curr.xb);
                }
                old.k += 1;
                curr.k += 1;
            }
        }
    }
}
