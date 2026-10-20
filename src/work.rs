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
    log: bool,
}

impl Default for Work {
    fn default() -> Self {
        Work { w: 0.0, log: std::env::var_os("ROSE_WORK").is_some() }
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

    /// Marks the reference's refactorization check (`0x142350`), which separates the iterations in the logs.
    pub fn trig(&self) {
        if self.log {
            eprintln!("TRIG {}", self.w);
        }
    }
}
