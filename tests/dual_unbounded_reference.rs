use rose::{mps, solve, verify, SolveOptions, Status};

#[test]
fn reference_dual_status_on_contradictions_and_infinite_lower_bound() {
    let independent = mps::parse(include_str!("fixtures/fx_infu_independent.mps")).unwrap();
    let nonempty = mps::parse(include_str!("fixtures/fx_infu_nonempty.mps")).unwrap();
    let upper = mps::parse(include_str!("fixtures/fx_unbounded_p1_mi.mps")).unwrap();
    for opts in [SolveOptions::dual(), SolveOptions::dual_release()] {
        let r = solve(&independent, None, &opts).unwrap();
        assert_eq!((r.status, r.iterations), (Status::Unbounded, 0));
        assert_eq!(r.ray, Some(vec![0.0, 1.0]));
        assert!(verify::check(&independent, &r, 1e-6).is_err());
        let r = solve(&nonempty, None, &opts).unwrap();
        assert_eq!((r.status, r.iterations), (Status::Infeasible, 1));
        assert!(verify::check(&nonempty, &r, 1e-6).is_ok());
        let r = solve(&upper, None, &opts).unwrap();
        assert_eq!(r.status, Status::Optimal);
        assert_eq!(r.objective, -42.0);
        assert!(verify::check(&upper, &r, 1e-6).is_ok());
    }
}

#[test]
fn valid_unbounded_reference_fixtures_verify() {
    for (name, source) in [
        ("unbounded", include_str!("fixtures/fx_unbounded.mps")),
        ("unbounded2", include_str!("fixtures/fx_unbounded2.mps")),
        ("unbounded_p1", include_str!("fixtures/fx_unbounded_p1.mps")),
        ("eq", include_str!("fixtures/fx_unb_eq.mps")),
        ("free", include_str!("fixtures/fx_unb_free.mps")),
        ("scaled", include_str!("fixtures/fx_unb_scaled.mps")),
        ("col", include_str!("fixtures/fx_unb_col.mps")),
        ("mi", include_str!("fixtures/fx_unb_mi.mps")),
    ] {
        let m = mps::parse(source).unwrap();
        for opts in [SolveOptions::dual(), SolveOptions::dual_release()] {
            let r = solve(&m, None, &opts).unwrap();
            assert_eq!(r.status, Status::Unbounded, "{name}");
            // Status 5 does not certify a primal point. Obtain a separate feasible witness for the strict checker, in tests only.
            let mut zero = m.clone();
            zero.obj.fill(0.0);
            let point = solve(&zero, None, &opts).unwrap();
            assert_eq!(point.status, Status::Optimal, "{name}");
            let mut certified = r;
            certified.x = point.x;
            assert!(
                verify::check(&m, &certified, 1e-6).is_ok(),
                "{name}: {:?}",
                verify::check(&m, &certified, 1e-6)
            );
        }
    }
}
