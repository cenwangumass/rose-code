//! NEW FILE (reference-corpus cases C11/C12/C14). Only cases where rose agrees with the controlled-profile reference
//! (statuses, x, ray, entering/leaving). Fixtures: reverse-engineer/corpus/fixtures/fx_*.mps; skipped if absent.
//! Known disagreements are NOT tested here: fx_tied (alternate optimum), fx_singular (supplied singular basis).

use rose::{mps, solve, Basis, LpModel, SolveOptions, Status};

fn fixture(name: &str) -> Option<LpModel> {
    let p = format!("{}/../reverse-engineer/corpus/fixtures/{name}.mps", env!("CARGO_MANIFEST_DIR"));
    std::fs::read_to_string(p).ok().map(|t| mps::parse(&t).unwrap())
}

#[test]
fn infeasible_and_unbounded_statuses_match_reference() {
    if let Some(m) = fixture("fx_infeasible") {
        let r = solve(&m, None, &SolveOptions::default()).unwrap();
        assert_eq!(r.status, Status::Infeasible);
        assert!(r.farkas.is_some());
    }
    if let Some(m) = fixture("fx_unbounded") {
        let r = solve(&m, None, &SolveOptions::default()).unwrap();
        assert_eq!(r.status, Status::Unbounded);
        assert_eq!(r.ray, Some(vec![1.0, 1.0])); // reference UnbdRay is [1, 1] too
    }
}

#[test]
fn neardep_and_badscale_solutions_match_reference_bits() {
    for (name, x, z) in [("fx_neardep", vec![0.0, 4.0], -8.0), ("fx_badscale", vec![0.0, 5e-6], -5.0)] {
        if let Some(m) = fixture(name) {
            let r = solve(&m, None, &SolveOptions::default()).unwrap();
            assert_eq!(r.status, Status::Optimal, "{name}");
            assert_eq!(r.x, x, "{name}");
            assert_eq!(r.objective, z, "{name}");
        }
    }
}

#[test]
fn warm_starts_from_nonoptimal_basis_match_reference() {
    let Some(m) = fixture("fx_warm") else { return };
    // y basic, slack c1 basic, x at lower: reference takes 1 iteration (a bound flip of x to its upper bound 3).
    let b = Basis { basic: vec![2, 1], at_upper: vec![false; 4] };
    let r = solve(&m, Some(&b), &SolveOptions::default()).unwrap();
    assert_eq!((r.status, r.iterations, r.x.clone()), (Status::Optimal, 1, vec![3.0, 1.0]));
    // same with x already at its upper bound: reference takes 0 iterations.
    let b = Basis { basic: vec![2, 1], at_upper: vec![true, false, false, false] };
    let r = solve(&m, Some(&b), &SolveOptions::default()).unwrap();
    assert_eq!((r.status, r.iterations, r.x.clone()), (Status::Optimal, 0, vec![3.0, 1.0]));
}
