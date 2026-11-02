//! `witness_primal_r3_fresh_high` (`reverse-engineer/corpus/fixtures`): the first sixteen structurals plus `s0`, `s1` as the start basis, one release solve. The
//! reference factors the start basis at thresholds 1/128, 0.5, 0.9 (near-dependent `s0`, `s1`), so at iteration 8 the fresh factor (zero updates) fails the
//! previous-column residual (1.52587890625e-5 > 3e-6 times 2.83) at threshold 0.9 with `ctx+0x3f8` = 1 and takes `0x1453c7 -> 0x1453d0 -> 0x145403`, the ordinary
//! rejection (an empty trace row, refactorization). The termination `0x1453d2` needs `ctx+0x3f8` = 0 (nested MIP) and is not reached.
//! Reference: `compare.py seq --rounds 1 --basic 18` (11 pivots) and the runner's final objective, x, pi, rc.
use rose::{mps, solve, Basis, SolveOptions, Status};

#[test]
fn fresh_factor_high_threshold_rejection() {
    let model = mps::parse(include_str!("fixtures/witness_primal_r3_fresh_high.mps")).unwrap();
    let (n, m) = (model.num_cols(), model.num_rows());
    assert_eq!(m, 18);
    let basis = Basis {
        basic: (0..m).collect(),
        at_upper: vec![false; n + m],
    };
    let r = solve(&model, Some(&basis), &SolveOptions::release()).unwrap();
    assert_eq!(
        (r.status, r.iterations, r.objective.to_bits()),
        (Status::Optimal, 11, 0xc00a661bb72c889e)
    );
}
