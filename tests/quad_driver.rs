use rose::{mps, solve, verify, SolveOptions, Status};

#[test]
fn primal_stale_pivot_rejection_carries_stability_to_quad_dual() {
    let model = mps::parse(include_str!("fixtures/quad_3260.mps")).unwrap();
    let result = solve(&model, None, &SolveOptions::dual_release()).unwrap();
    assert_eq!((result.status, result.iterations), (Status::Infeasible, 96));
}

#[test]
fn shifted_quad_handoff_checks_recomputed_reduced_costs() {
    let model = mps::parse(include_str!("fixtures/quad_8179.mps")).unwrap();
    let result = solve(&model, None, &SolveOptions::dual_release()).unwrap();
    assert_eq!(
        (result.status, result.iterations),
        (Status::NumericalFailure, 152)
    );
}

#[test]
fn shifted_quad_handoff_preserves_reference_infeasible_readout() {
    let model = mps::parse(include_str!("fixtures/quad_5885.mps")).unwrap();
    let result = solve(&model, None, &SolveOptions::dual_release()).unwrap();
    assert_eq!(
        (result.status, result.iterations),
        (Status::Infeasible, 116)
    );
    let expected = [
        0.0,
        -6.507161132832055,
        -0.0006507161132571568,
        0.0,
        -3.2555340950683644e-09,
        0.0,
        0.00390625,
        0.0,
        -1.0417317773488966e-08,
        1.0417317779753883e-10,
    ];
    assert_eq!(result.farkas.as_deref(), Some(expected.as_slice()));
    // The reference returns this raw vector, but it is not a verified proof.
    assert!(verify::farkas_margin(&model, &expected) <= 0.0);
}

#[test]
fn shifted_quad_handoff_repairs_basis_at_the_carried_threshold() {
    let model = mps::parse(include_str!("fixtures/quad_2594.mps")).unwrap();
    let result = solve(&model, None, &SolveOptions::dual_release()).unwrap();
    assert_eq!(
        (result.status, result.iterations),
        (Status::Infeasible, 110)
    );
}

#[test]
fn x87_pure_penalty_rescans_updated_duals_before_ending() {
    let model = mps::parse(include_str!("fixtures/quad_2717.mps")).unwrap();
    let result = solve(&model, None, &SolveOptions::dual_release()).unwrap();
    assert_eq!((result.status, result.iterations), (Status::Optimal, 16010));
}
