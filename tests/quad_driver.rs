use rose::{mps, solve, verify, SolveOptions, Status};

#[test]
fn quad_primal_repeated_repairs_resume_dual() {
    let model = mps::parse(include_str!("fixtures/quad_4737.mps")).unwrap();
    let result = solve(&model, None, &SolveOptions::dual_release()).unwrap();
    assert_eq!(
        (result.status, result.iterations),
        (Status::NumericalFailure, 228)
    );
}

#[test]
fn x87_tiny_pivot_retries_the_active_pricing_bucket() {
    for (text, status, iterations) in [
        (include_str!("fixtures/quad_6230.mps"), Status::Infeasible, 40),
        (include_str!("fixtures/quad_4564.mps"), Status::Infeasible, 117),
    ] {
        let model = mps::parse(text).unwrap();
        let result = solve(&model, None, &SolveOptions::dual_release()).unwrap();
        assert_eq!((result.status, result.iterations), (status, iterations));
    }
}

#[test]
fn shifted_primal_handoff_continues_from_its_basis_with_a_snapshot() {
    for (input, status, iterations) in [
        (include_str!("fixtures/quad_4918.mps"), Status::Optimal, 182),
        (include_str!("fixtures/quad_5147.mps"), Status::Infeasible, 194),
    ] {
        let model = mps::parse(input).unwrap();
        let result = solve(&model, None, &SolveOptions::dual_release()).unwrap();
        assert_eq!((result.status, result.iterations), (status, iterations));
    }
}

#[test]
fn quad_residual_uses_extended_precision() {
    let model = mps::parse(include_str!("fixtures/quad_7374.mps")).unwrap();
    let result = solve(&model, None, &SolveOptions::dual_release()).unwrap();
    assert_eq!((result.status, result.iterations), (Status::NumericalFailure, 79));
}

#[test]
fn shifted_quad_handoff_retries_until_restore_limit() {
    let model = mps::parse(include_str!("fixtures/quad_3076.mps")).unwrap();
    let result = solve(&model, None, &SolveOptions::dual_release()).unwrap();
    assert_eq!(
        (result.status, result.iterations),
        (Status::NumericalFailure, 189)
    );
}

#[test]
fn primal_stale_pivot_rejection_carries_stability_to_quad_dual() {
    let model = mps::parse(include_str!("fixtures/quad_3260.mps")).unwrap();
    let result = solve(&model, None, &SolveOptions::dual_release()).unwrap();
    assert_eq!((result.status, result.iterations), (Status::Infeasible, 96));
}

#[test]
fn quad_switch_keeps_iteration_count_for_infinite_bound_ratio() {
    let model = mps::parse(include_str!("fixtures/quad_8151.mps")).unwrap();
    let result = solve(&model, None, &SolveOptions::dual_release()).unwrap();
    assert_eq!((result.status, result.iterations), (Status::Infeasible, 31));
}

#[test]
fn quad_switch_carries_iteration_count_to_ratio_selection() {
    let model = mps::parse(include_str!("fixtures/quad_5103.mps")).unwrap();
    let result = solve(&model, None, &SolveOptions::dual_release()).unwrap();
    assert_eq!((result.status, result.iterations), (Status::Optimal, 27));
}

#[test]
fn quad_singular_restart_restores_initial_snapshot() {
    let model = mps::parse(include_str!("fixtures/quad_7973.mps")).unwrap();
    let result = solve(&model, None, &SolveOptions::dual_release()).unwrap();
    assert_eq!((result.status, result.iterations), (Status::Optimal, 54));
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

#[test]
fn shifted_quad_retry_selects_pricing_from_handoff_count() {
    let model = mps::parse(include_str!("fixtures/quad_3196.mps")).unwrap();
    let result = solve(&model, None, &SolveOptions::dual_release()).unwrap();
    assert_eq!((result.status, result.iterations), (Status::Infeasible, 142));
}
