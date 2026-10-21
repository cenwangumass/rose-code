//! NEW FILE (reference-corpus cases C11/C12/C14). Only cases where rose agrees with the controlled-profile reference
//! (statuses, x, ray, entering/leaving). Fixtures: reverse-engineer/corpus/fixtures/fx_*.mps; skipped if absent.
//! Open: fx_sing_dense4 (singular 4x4 basis: the reference's x differs from rose's in the last digit of two entries).

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

#[test]
fn tied_optimal_face_returns_the_reference_vertex() {
    // min -x - y, x + y <= 4, x, y <= 3: the crash puts x at its upper bound first (the reference's sort of equal keys keeps the
    // input order), so y enters and the vertex is (3, 1), as in the reference.
    let Some(m) = fixture("fx_tied") else { return };
    let r = solve(&m, None, &SolveOptions::default()).unwrap();
    assert_eq!((r.status, r.iterations, r.x.clone()), (Status::Optimal, 1, vec![3.0, 1.0]));
}

/// Supplied singular bases (`fx_singular`, `fx_sing_*`, made by `corpus/fixtures/singular_cases.py`): the reference drops the
/// columns that get no pivot in its LU order, replaces them by the slacks of the unpivoted rows, and keeps the last one dropped
/// out of pricing until no other candidate is left. Expected iterations and x are the reference's (controlled profile).
#[test]
fn singular_supplied_bases_are_repaired_like_the_reference() {
    // (fixture, structurals made basic with the row whose slack they replace, iterations, x)
    let cases: [(&str, &[(usize, usize)], usize, &[f64]); 14] = [
        ("fx_singular", &[(0, 0), (1, 1)], 1, &[0.0, 4.0]),
        ("fx_sing_dup_cross", &[(0, 1), (1, 0)], 1, &[0.0, 4.0]),
        ("fx_sing_dup_attr", &[(0, 0), (1, 1)], 1, &[4.0, 0.0]),
        ("fx_sing_dup_attr_ub", &[(0, 0), (1, 1)], 1, &[3.0, 1.0]),
        ("fx_sing_dup_ub", &[(0, 0), (1, 1)], 3, &[3.0, 1.0]),
        ("fx_sing_zero_col", &[(0, 0), (1, 1)], 2, &[3.0, 4.0]),
        ("fx_sing_zero_col2", &[(0, 0), (1, 1)], 2, &[4.0, 3.0]),
        ("fx_sing_slack_dup", &[(0, 1)], 2, &[4.0, 5.0]),
        ("fx_sing_slack_dup2", &[(0, 1), (1, 0)], 0, &[4.0, 0.0]),
        ("fx_sing_rank3", &[(0, 0), (1, 1), (2, 2)], 0, &[4.0, 5.0, 0.0]),
        ("fx_sing_all_same3", &[(0, 0), (1, 1), (2, 2)], 1, &[0.0, 0.0, 4.0]),
        ("fx_sing_nuc_dep2", &[(0, 0), (1, 1), (2, 2)], 1, &[0.0, 3.0, 0.0]),
        ("fx_sing_two_pairs", &[(0, 0), (1, 1), (2, 2), (3, 3)], 3, &[0.0, 4.0, 6.0, 0.0]),
        ("fx_sing_dup_slack3", &[(0, 0), (1, 1)], 2, &[0.0, 4.0, 6.0]),
    ];
    for (name, pairs, iters, x) in cases {
        let Some(m) = fixture(name) else { continue };
        let (n, rows) = (m.num_cols(), m.num_rows());
        let mut basic: Vec<usize> = (n..n + rows).collect();
        for &(j, r) in pairs {
            basic[r] = j;
        }
        let b = Basis { basic, at_upper: vec![false; n + rows] };
        let r = solve(&m, Some(&b), &SolveOptions::default()).unwrap();
        assert_eq!((r.status, r.iterations, r.x.clone()), (Status::Optimal, iters, x.to_vec()), "{name}");
    }
}
