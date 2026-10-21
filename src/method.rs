//! The release profile's method decision (`0x142500`): it starts as `SimplexPricing=0` (`crate::simplex::SolveOptions::pricing0`) and, at a
//! refactorization, switches to Devex pricing when a work estimate favours it. Rules: `reverse-engineer/specs/primal-simplex.md`,
//! "Method switch". Captured inputs: `reverse-engineer/replay/gdb_switch.py`, `compare.py switch`.

/// Counters of the iteration step (`ctx+0x2a0`, `0x330`, `0x334`, `0x2b0`) and the decision made from them.
#[derive(Clone, Debug, Default)]
pub struct Method {
    /// Iterations so far, bound flips included (`ctx+0x2a0`).
    pub iters: f64,
    /// Iterations with a step of at most the tolerance (`ctx+0x330`), not counting a free entering variable, a fixed leaving variable or a flip.
    pub degenerate: i32,
    /// `degenerate` and `iters` at the last snapshot (`ctx+0x334`, `ctx+0x2b0`).
    pub snap_degenerate: i32,
    pub snap_iters: f64,
    /// Factorizations so far, the initial one included (`ctx+0x2c0`).
    pub factorizations: u32,
    /// The method is Devex pricing (`ctx+0x3ec == 1`).
    pub devex: bool,
}

/// What the decision reads from the factor and the model.
pub struct Inputs {
    /// Rows and columns (`ctx+0x64`, `ctx+0x68`).
    pub rows: usize,
    pub cols: usize,
    /// Mean FTRAN result count of the updates since the last factorization, truncated (`F+0x1f0`).
    pub avg: usize,
    /// `nnz(L) + rows` and `nnz(U) + rows` of the new factor (`F+0x1e0`, `F+0x1e8`).
    pub l_total: usize,
    pub u_total: usize,
    /// Basic variables violate their bounds (`ctx+0x2d0 != 0`, `P+0xc == 1`).
    pub infeasible: bool,
}

/// Tolerance of the degenerate-step count (`ctx+0x380`).
const STEP_TOL: f64 = 1e-6;

impl Method {
    /// An iteration was made (`0x143a22`); `step` is its step length and `counted` whether it has a leaving variable that is not fixed
    /// and an entering variable that is not free.
    pub fn iteration(&mut self, flip: bool, counted: bool, step: f64) {
        self.iters += 1.0;
        if !flip && counted && step.abs() <= STEP_TOL {
            self.degenerate += 1;
        }
    }

    /// The decision after a factorization (`0x142500`, the work-based branch for the default parameters; state 10 or 2): `true` switches to
    /// Devex pricing. Afterwards the snapshot of the pricing driver (`0xfb207`) is taken on odd factorization counts.
    pub fn decide(&mut self, inp: &Inputs) -> bool {
        self.factorizations += 1;
        let (a, b) = (self.iters, self.snap_iters);
        let mut switch = false;
        if !self.devex && a != 0.0 && a != b {
            let rows = inp.rows as f64;
            let cols = inp.cols as f64;
            let mut big = a / (rows + cols)
                + 0.3
                + 1.2 * inp.avg as f64 / rows
                + 0.6 * (f64::from(self.degenerate - self.snap_degenerate) / (a - b));
            if inp.infeasible {
                big *= 0.5;
            }
            // `P+0x3b0`: rows * (rows^2 / cols)^0.35 / 20.
            let p3b0 = if inp.cols == 0 { 1.0 } else { rows * (rows * rows / cols).powf(0.35) / 20.0 };
            let denom = inp.u_total as f64 - 2.0 * rows + 3.0 * inp.l_total as f64;
            let rhs = p3b0 / (inp.avg as f64 + 1.0) + 1.0 + 9.0 * cols / denom;
            switch = big + big > rhs;
        }
        if switch {
            self.devex = true;
        }
        if self.factorizations % 2 == 1 {
            self.snap_degenerate = self.degenerate;
            self.snap_iters = self.iters;
        }
        switch
    }
}
