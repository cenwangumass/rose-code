//! The release profile's method decision (`0x142500`): it starts as `SimplexPricing=0` (`crate::simplex::SolveOptions::pricing0`) and, at a
//! refactorization, switches to Devex pricing when a work estimate favours it. Rules: `reverse-engineer/specs/primal-simplex.md`,
//! "Method switch" (method 2 is the second test of `decide`; see `crate::simplex::phase1_eefa0`). Captured inputs: `reverse-engineer/replay/gdb_switch.py`, `compare.py switch`.

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
    /// The method (`ctx+0x3ec`): 0 `pricing0`, 1 Devex pricing, 2 Devex pricing reached by the second test of `decide`. A non-zero mode also
    /// clears `ctx+0x468`, which selects the phase 1 ratio test `0xeefa0` (`crate::simplex::phase1_eefa0`).
    pub mode: u8,
    /// The release profile's bookkeeping runs (counters, snapshot, perturbation); off in the controlled profile.
    pub active: bool,
    /// Perturbation level (`ctx+0x338`, +2 per perturbation) and size (`ctx+0x340`, tripled by each one; 0 disables).
    pub level: u32,
    pub pert: f64,
    /// The last `decide` found the perturbation due (`0xfb207..0xfb283`): the caller runs `Solver::perturb`.
    pub perturb_due: bool,
    /// Mode 2's weight-error average (`ctx+0x408`, `ctx+0x410`): the sum of `max(wq / w_q, w_q / wq)` over the iterations since the weights
    /// were initialised, and the count (starting at 1). The weights restart when the average exceeds 3 (`0x10882c`).
    pub ratio_sum: f64,
    pub ratio_cnt: u32,
    /// `ctx+0x468` is still set although the mode is not 0 (the primal that continues a dual run starts in mode 2 without the clearing): the phase 1 ratio test
    /// stays `0xf0120` until the decision switches the method.
    pub keep468: bool,
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
    /// The perturbation routine (`0x141a60`) ran: the level rises by 2 and the size triples (while below `1e6 * PerturbValue`).
    pub fn perturbed(&mut self) {
        self.level += 2;
        if 1e6 * 0.0002 > self.pert {
            self.pert *= 3.0;
        }
    }

    pub fn iteration(&mut self, flip: bool, counted: bool, step: f64) {
        self.iters += 1.0;
        if !flip && counted && step.abs() <= STEP_TOL {
            self.degenerate += 1;
        }
    }

    /// The decision after a factorization (`0x142500`, the work-based branch for the default parameters; state 10 or 2): `true` enters Devex
    /// pricing (mode 1 or 2; the weights are initialised again). Afterwards the snapshot of the pricing driver (`0xfb207`) is taken when the factorization count mod (level + 2) is 1.
    pub fn decide(&mut self, inp: &Inputs) -> bool {
        self.factorizations += 1;
        let (a, b) = (self.iters, self.snap_iters);
        let mut switch = false;
        if self.mode != 1 && a != 0.0 && a != b {
            let rows = inp.rows as f64;
            let cols = inp.cols as f64;
            // `ctx+0x338 >> 1` (the perturbation level, even) joins the degenerate share.
            let mut big = a / (rows + cols)
                + 0.3
                + 1.2 * inp.avg as f64 / rows
                + 0.6
                    * (f64::from(self.level as i32 >> 1)
                        + f64::from(self.degenerate - self.snap_degenerate) / (a - b));
            // Halved in phase 1 while `ctx+0x468 != 0`, that is before any switch.
            if inp.infeasible && (self.mode == 0 || self.keep468) {
                big *= 0.5;
            }
            // `P+0x3b0`: rows * (rows^2 / cols)^0.35 / 20.
            let p3b0 = if inp.cols == 0 {
                1.0
            } else {
                rows * (rows * rows / cols).powf(0.35) / 20.0
            };
            let denom = inp.u_total as f64 - 2.0 * rows + 3.0 * inp.l_total as f64;
            let rhs = p3b0 / (inp.avg as f64 + 1.0) + 1.0 + 9.0 * cols / denom;
            if big + big > rhs {
                self.mode = 1;
                self.keep468 = false;
                switch = true;
            } else if self.mode == 0 && 1.2 * big > 1.0 + 6.0 * cols / denom {
                self.mode = 2;
                self.keep468 = false;
                switch = true;
            }
        }
        self.perturb_due = false;
        if self.factorizations % (self.level + 2) == 1 {
            // The pricing driver (`0xfb207`) perturbs when the degenerate iterations since the snapshot exceed 0.9 of all of them.
            if self.pert > 0.0 && f64::from(self.degenerate - self.snap_degenerate) > (a - b) * 0.9
            {
                self.perturb_due = true;
            }
            self.snap_degenerate = self.degenerate;
            self.snap_iters = self.iters;
        }
        switch
    }
}
