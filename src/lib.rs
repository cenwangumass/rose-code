//! Standalone LP solver. The numerical core is our own implementation; no external solver is used.

pub mod barrier;
pub mod basis;
pub mod chol;
pub mod dual;
pub mod dump;
pub mod early;
pub mod f80;
pub mod ftran_order;
pub mod lu;
pub mod lu_order;
pub mod lu_pairs;
pub mod method;
pub mod model;
pub mod mps;
pub mod order;
pub mod pdhg;
mod primal_rhs_perturb;
pub mod real;
pub mod scale;
pub mod sgemm32;
pub mod simplex;
pub mod work;

pub use basis::Basis;
pub use model::{LpModel, Sense};

pub mod verify;

pub use simplex::{crash_basis, solve, Method, SolveOptions, SolveResult, Status, TraceEvent};
