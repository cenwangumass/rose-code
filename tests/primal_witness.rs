use rose::{mps, solve, SolveOptions, Status};

fn release(model: &str) -> (Status, usize, f64) {
    let r = solve(&mps::parse(model).unwrap(), None, &SolveOptions::release()).unwrap();
    (r.status, r.iterations, r.objective)
}

#[test]
fn residual_above_ten_tolerances_ends_the_primal_without_a_follow_up_pass() {
    // `witness_primal_residual_loop` (`0xd0382`, `r14d == -1`): the unscaled primal residual is 1.4e-4 with the gate open, so the end test returns without the dual
    // clean-up pass: the reference ends Optimal after 7 iterations at 44166500.00015778 (rose ran a clean-up pivot and ended with 8 and 44166499.99586613).
    let (s, n, obj) = release(include_str!("fixtures/witness_primal_residual_loop.mps"));
    assert_eq!((s, n, obj), (Status::Optimal, 7, 44166500.00015778));
}

#[test]
fn small_unscaled_violation_runs_the_second_primal_pass() {
    // `witness_primal_small_bound_second`: unscaled violation 9.9e-7 with the gate open and a small residual; the reference runs the unscaled second primal pass and
    // ends Optimal after 9 iterations.
    let (s, n, obj) = release(include_str!("fixtures/witness_primal_small_bound_second.mps"));
    assert_eq!((s, n, obj), (Status::Optimal, 9, 414273.1355091953));
}

#[test]
fn accuracy_ladder_witnesses_agree_with_the_reference_iterations() {
    // Natural low-threshold rejection (`witness_primal_accuracy_low`) and the five `e2` rungs of the residual ladder (thresholds 0.03125 .. 0.5): the release traces
    // agree with the reference on every pivot (`compare.py trace --decisions`).
    let counts = [
        (include_str!("fixtures/witness_primal_accuracy_low.mps"), 9),
        (include_str!("fixtures/witness_primal_e2_03125.mps"), 19),
        (include_str!("fixtures/witness_primal_e2_0625.mps"), 23),
        (include_str!("fixtures/witness_primal_e2_125.mps"), 16),
        (include_str!("fixtures/witness_primal_e2_25.mps"), 13),
        (include_str!("fixtures/witness_primal_e2_5.mps"), 12),
    ];
    for (model, want) in counts {
        let (s, n, _) = release(model);
        assert_eq!((s, n), (Status::Optimal, want));
    }
}
