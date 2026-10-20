//! The reference's deterministic work counter (`W`, `reverse-engineer/specs/primal-simplex.md`, "Refactor time clause").
//!
//! Each routine of the reference adds an estimate of its own loop trip counts, scaled by a per-object factor: 1.43 for the
//! solver context (`ctx+0x148`) and 1.3 for the factor object (`F+0x1f8`). The refactorization trigger compares the work since
//! the last factorization with the factorization's own cost. Sites are named by the module-relative address of the increment in
//! the reference, so a log (`ROSE_WORK`) can be compared with a GDB capture (`replay/gdb_wwatch.py`).

/// Scale of the solver context's estimates (`ctx+0x148`).
pub const CTX: f64 = 1.43;
/// Scale of the factor object's estimates (`F+0x1f8`).
pub const FAC: f64 = 1.3;

pub struct Work {
    pub w: f64,
    /// `W` at the refactorization check that last decided to refactorize (`ctx+0x150`); zero before the first.
    w150: f64,
    /// Whether a check has decided to refactorize yet. The reference's `w150 > 0` holds after the setup work (not modelled in
    /// `W`) of the first check, so it is a flag here.
    decided_once: bool,
    /// Smoothed cost of the factorization driver `0xd4190` (`ctx+0x158`): half the last cost plus half the previous value.
    w158: f64,
    /// Cost of the factorization in progress: the sum of the driver's own increments.
    fac_cost: f64,
    log: bool,
}

impl Default for Work {
    fn default() -> Self {
        Work { w: 0.0, w150: 0.0, decided_once: false, w158: 0.0, fac_cost: 0.0, log: std::env::var_os("ROSE_WORK").is_some() }
    }
}

impl Work {
    /// `W = est * scale + W`, in the reference's operation order.
    pub fn add(&mut self, site: u32, est: f64, scale: f64) {
        let d = est * scale;
        self.w = d + self.w;
        if self.log {
            eprintln!("WINC {d:.4} {site:#x}");
        }
    }

    /// As `add`, for a shared routine (solve wrapper, row product) whose increment is told apart by its call site.
    pub fn add_from(&mut self, site: u32, caller: u32, est: f64, scale: f64) {
        let d = est * scale;
        self.w = d + self.w;
        if self.log {
            eprintln!("WINC {d:.4} {site:#x} {caller:#x}");
        }
    }

    /// The time clause of the refactorization check: every tenth update, when the work since the last decision exceeds
    /// `8 * w158 + 100 * nnz(A)`. `w` is `W` at the check.
    pub fn clause(&self, upd: usize, nnz_a: usize, w: f64) -> bool {
        upd % 10 == 9 && self.decided_once && w - self.w150 > 8.0 * self.w158 + 100.0 * nnz_a as f64
    }

    /// Records a check that decided to refactorize, at `W = w`.
    pub fn decided(&mut self, w: f64) {
        self.w150 = w;
        self.decided_once = true;
    }

    /// Adds the increments of the factorization driver `0xd4190` (the sites that make up its cost) and, at its end, smooths
    /// the cost into `w158`.
    pub fn factor_driver(&mut self, sites: &[(u32, f64)], scale: f64) {
        self.fac_cost = 0.0;
        for &(site, est) in sites {
            if matches!(site, 0xd4ded | 0xd5300 | 0xd85e6 | 0xdc7b2 | 0xdcb7b) {
                self.fac_cost += est * scale;
            }
        }
    }

    /// Smooths the cost of the factorization just done into `w158`.
    pub fn factor_done(&mut self) {
        // The first factorization stores its cost as it is; later ones average it with the previous value.
        self.w158 = if self.w158 == 0.0 { self.fac_cost } else { 0.5 * (self.fac_cost + self.w158) };
        if self.log {
            eprintln!("W158 {} cost {} w150 {}", self.w158, self.fac_cost, self.w150);
        }
    }

    /// Marks the reference's refactorization check (`0x142350`), which separates the iterations in the logs.
    pub fn trig(&self) {
        if self.log {
            eprintln!("TRIG {}", self.w);
        }
    }
}
