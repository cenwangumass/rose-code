use rose::{mps, solve, SolveOptions, Status};

#[test]
fn shifted_quad_handoff_checks_recomputed_reduced_costs() {
    let model = mps::parse(include_str!("fixtures/quad_8179.mps")).unwrap();
    let result = solve(&model, None, &SolveOptions::dual_release()).unwrap();
    assert_eq!(
        (result.status, result.iterations),
        (Status::NumericalFailure, 152)
    );
}
