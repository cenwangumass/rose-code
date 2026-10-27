use rose::{mps, solve, verify, SolveOptions, Status};

#[test]
fn rejected_dual_handoff_returns_the_primal_infeasibility() {
    let model = mps::parse(include_str!("fixtures/bh50000.mps")).unwrap();
    let controlled = solve(&model, None, &SolveOptions::dual()).unwrap();
    assert_eq!(
        (controlled.status, controlled.iterations),
        (Status::Infeasible, 14)
    );
    assert_eq!(
        controlled.farkas,
        Some(vec![
            -0.002180437468045592,
            0.517123627222175,
            -1.0,
            0.3066188460023535,
            -37.70483043996373,
            0.0,
            -0.00517123627171923,
            -0.12943650937193896,
            0.0,
            -0.019999999961840684,
            0.0,
            0.06284140406668667,
        ])
    );
    // The reference vector is reproduced, but it is not a valid certificate.
    assert!(verify::check(&model, &controlled, 1e-6).is_err());
    let release = solve(&model, None, &SolveOptions::dual_release()).unwrap();
    assert_eq!(
        (release.status, release.iterations),
        (Status::Infeasible, 13)
    );
}
