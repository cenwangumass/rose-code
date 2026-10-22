//! General-LP behaviour: cold starts, infeasible starts (Phase 1), certificates, warm starts, and the small
//! Netlib models against published objectives. Every result is validated by `rose::verify`, which recomputes
//! feasibility, duality and certificates from the original model.

use rose::verify;
use rose::{mps, solve, Basis, LpModel, SolveOptions, Status};

/// The methods every test runs under: the primal simplex and the dual simplex (the controlled profiles; pivot sequences differ).
fn opts() -> Vec<SolveOptions> {
    vec![SolveOptions::default(), SolveOptions::dual()]
}

fn fixture(name: &str) -> LpModel {
    let text = std::fs::read_to_string(format!("{}/tests/fixtures/{name}", env!("CARGO_MANIFEST_DIR"))).unwrap();
    mps::parse(&text).unwrap()
}

fn run(model: &LpModel, o: &SolveOptions) -> rose::SolveResult {
    let r = solve(model, None, o).unwrap();
    if let Err(e) = verify::check(model, &r, 1e-6) {
        panic!("independent check failed: {e}");
    }
    r
}

fn close(a: f64, b: f64) -> bool {
    (a - b).abs() <= 1e-7 + 1e-8 * a.abs().max(b.abs())
}

#[test]
fn netlib_objectives_from_cold_start() {
    for o in opts() {
        for (name, z) in [
            ("afiro.mps", -464.75314285714285),
            ("sc50a.mps", -64.57507705856452),
            ("blend.mps", -30.812149845828216),
            ("adlittle.mps", 225494.96316238036),
            ("boeing2.mps", -315.0187280152027),
            ("bore3d.mps", 1373.0803942084924),
        ] {
            let model = fixture(name);
            let r = run(&model, &o);
            assert_eq!(r.status, Status::Optimal, "{name}");
            assert!(close(r.objective, z), "{name}: {} vs {z}", r.objective);
        }
    }
}

#[test]
fn phase1_from_infeasible_start() {
    for o in opts() {
        // min x + y  s.t. x + y >= 2, x - y = 0: the slack basis violates both rows.
        let text = "NAME T\nROWS\n N obj\n G c1\n E c2\nCOLUMNS\n x obj 1 c1 1\n x c2 1\n y obj 1 c1 1\n y c2 -1\nRHS\n rhs c1 2\nENDATA\n";
        let r = run(&mps::parse(text).unwrap(), &o);
        assert_eq!(r.status, Status::Optimal);
        assert!(close(r.objective, 2.0));
    }
}

#[test]
fn detects_infeasible_with_certificate() {
    for o in opts() {
        let text = "NAME T\nROWS\n N obj\n L c1\n G c2\nCOLUMNS\n x obj 1 c1 1\n x c2 1\nRHS\n rhs c1 1 c2 3\nENDATA\n";
        let model = mps::parse(text).unwrap();
        let r = run(&model, &o);
        assert_eq!(r.status, Status::Infeasible);
        assert!(verify::farkas_margin(&model, r.farkas.as_ref().unwrap()) > 1e-6);
    }
}

#[test]
fn detects_unbounded_with_ray() {
    for o in opts() {
        let text = "NAME T\nROWS\n N obj\n G c1\nCOLUMNS\n x obj -1 c1 1\n y obj 0 c1 -1\nRHS\n rhs c1 1\nENDATA\n";
        let model = mps::parse(text).unwrap();
        let r = run(&model, &o);
        assert_eq!(r.status, Status::Unbounded);
    }
}

#[test]
fn ranged_rows_and_free_variables() {
    for o in opts() {
        // min x - y, 1 <= x + y <= 3 (ranged), y free, x in [0, 2], row2: y <= 5
        let text = "NAME T\nROWS\n N obj\n G c1\n L c2\nCOLUMNS\n x obj 1 c1 1\n y obj -1 c1 1\n y c2 1\nRHS\n rhs c1 1 c2 5\nRANGES\n rng c1 2\nBOUNDS\n UP bnd x 2\n FR bnd y\nENDATA\n";
        let model = mps::parse(text).unwrap();
        let r = run(&model, &o);
        assert_eq!(r.status, Status::Optimal);
        assert!(close(r.objective, -3.0), "{}", r.objective); // x = 0, y = 3 (upper end of the range)
    }
}

#[test]
fn warm_start_from_optimal_basis_takes_no_iterations() {
    for o in opts() {
        let model = fixture("blend.mps");
        let cold = run(&model, &o);
        let warm = solve(&model, Some(&cold.basis), &o).unwrap();
        assert_eq!(warm.status, Status::Optimal);
        assert_eq!(warm.iterations, 0);
        assert!(close(warm.objective, cold.objective));
    }
}

#[test]
fn invalid_basis_is_rejected() {
    for o in opts() {
        let model = fixture("afiro.mps");
        let bad = Basis { basic: vec![0; model.num_rows()], at_upper: vec![false; model.num_rows() + model.num_cols()] };
        assert!(solve(&model, Some(&bad), &o).is_err());
    }
}

#[test]
fn iteration_limit_is_reported_with_a_finite_point() {
    let model = fixture("blend.mps");
    for o in opts() {
        let r = solve(&model, None, &SolveOptions { max_iterations: 5, record_trace: false, ..o }).unwrap();
        assert_eq!(r.status, Status::IterationLimit);
        assert!(r.x.iter().all(|v| v.abs() < 1e20), "{:?}", r.x);
    }
}

#[test]
fn singular_supplied_basis_is_repaired() {
    // Two identical columns both basic: the second leaves for the slack of the uncovered row.
    let text = "NAME T\nROWS\n N obj\n L c1\n L c2\nCOLUMNS\n x obj -1 c1 1\n x c2 1\n y obj -2 c1 1\n y c2 1\nRHS\n rhs c1 4 c2 5\nENDATA\n";
    let model = mps::parse(text).unwrap();
    let singular = Basis { basic: vec![0, 1], at_upper: vec![false; 4] };
    for o in opts() {
        let r = solve(&model, Some(&singular), &o).unwrap();
        assert_eq!(r.status, Status::Optimal);
        assert!(close(r.objective, -8.0), "{}", r.objective);
        verify::check(&model, &r, 1e-6).unwrap();
    }
}
