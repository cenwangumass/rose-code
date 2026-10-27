use rose::{mps, solve, verify, SolveOptions, Status};

#[test]
fn deferred_free_column_ray_excludes_other_artificial_bounds() {
    let m = mps::parse(include_str!("fixtures/fx_dual_free_ray_122.mps")).unwrap();
    for opts in [SolveOptions::dual(), SolveOptions::dual_release()] {
        let r = solve(&m, None, &opts).unwrap();
        assert_eq!((r.status, r.iterations), (Status::Unbounded, 0));
        assert_eq!(r.ray, Some(vec![-1.0, 0.0, 0.0, 0.0, 0.0]));
        assert!(verify::ray_improvement(&m, r.ray.as_ref().unwrap(), 1e-6).is_some_and(|g| g > 0.0));
    }
}
