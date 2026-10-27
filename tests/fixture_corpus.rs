//! NEW FILE (reference-corpus cases C11/C12/C14). Only cases where rose agrees with the controlled-profile reference
//! (statuses, x, ray, entering/leaving). Fixtures: reverse-engineer/corpus/fixtures/fx_*.mps; skipped if absent.

use rose::{mps, solve, Basis, LpModel, SolveOptions, Status};

fn fixture(name: &str) -> Option<LpModel> {
    let p = format!(
        "{}/../reverse-engineer/corpus/fixtures/{name}.mps",
        env!("CARGO_MANIFEST_DIR")
    );
    std::fs::read_to_string(p)
        .ok()
        .map(|t| mps::parse(&t).unwrap())
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
    for (name, x, z) in [
        ("fx_neardep", vec![0.0, 4.0], -8.0),
        ("fx_badscale", vec![0.0, 5e-6], -5.0),
    ] {
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
    let b = Basis {
        basic: vec![2, 1],
        at_upper: vec![false; 4],
    };
    let r = solve(&m, Some(&b), &SolveOptions::default()).unwrap();
    assert_eq!(
        (r.status, r.iterations, r.x.clone()),
        (Status::Optimal, 1, vec![3.0, 1.0])
    );
    // same with x already at its upper bound: reference takes 0 iterations.
    let b = Basis {
        basic: vec![2, 1],
        at_upper: vec![true, false, false, false],
    };
    let r = solve(&m, Some(&b), &SolveOptions::default()).unwrap();
    assert_eq!(
        (r.status, r.iterations, r.x.clone()),
        (Status::Optimal, 0, vec![3.0, 1.0])
    );
}

#[test]
fn tied_optimal_face_returns_the_reference_vertex() {
    // min -x - y, x + y <= 4, x, y <= 3: the crash puts x at its upper bound first (the reference's sort of equal keys keeps the
    // input order), so y enters and the vertex is (3, 1), as in the reference.
    let Some(m) = fixture("fx_tied") else { return };
    let r = solve(&m, None, &SolveOptions::default()).unwrap();
    assert_eq!(
        (r.status, r.iterations, r.x.clone()),
        (Status::Optimal, 1, vec![3.0, 1.0])
    );
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
        (
            "fx_sing_rank3",
            &[(0, 0), (1, 1), (2, 2)],
            0,
            &[4.0, 5.0, 0.0],
        ),
        (
            "fx_sing_all_same3",
            &[(0, 0), (1, 1), (2, 2)],
            1,
            &[0.0, 0.0, 4.0],
        ),
        (
            "fx_sing_nuc_dep2",
            &[(0, 0), (1, 1), (2, 2)],
            1,
            &[0.0, 3.0, 0.0],
        ),
        (
            "fx_sing_two_pairs",
            &[(0, 0), (1, 1), (2, 2), (3, 3)],
            3,
            &[0.0, 4.0, 6.0, 0.0],
        ),
        ("fx_sing_dup_slack3", &[(0, 0), (1, 1)], 2, &[0.0, 4.0, 6.0]),
    ];
    for (name, pairs, iters, x) in cases {
        let Some(m) = fixture(name) else { continue };
        let (n, rows) = (m.num_cols(), m.num_rows());
        let mut basic: Vec<usize> = (n..n + rows).collect();
        for &(j, r) in pairs {
            basic[r] = j;
        }
        let b = Basis {
            basic,
            at_upper: vec![false; n + rows],
        };
        let r = solve(&m, Some(&b), &SolveOptions::default()).unwrap();
        assert_eq!(
            (r.status, r.iterations, r.x.clone()),
            (Status::Optimal, iters, x.to_vec()),
            "{name}"
        );
    }
}

/// The same fixtures under the dual method, both profiles (controlled `dual()` and `dual_release()`): statuses, certificates, objective and x to
/// 1e-9, and `verify::check` of every result. Iteration counts are not asserted: the singular-basis and warm-start fixtures differ from the
/// reference's counts (documented in the README, "Recovery gaps").
#[test]
fn fixtures_under_the_dual_method_are_correct() {
    for o in [SolveOptions::dual(), SolveOptions::dual_release()] {
        if let Some(m) = fixture("fx_infeasible") {
            let r = solve(&m, None, &o).unwrap();
            assert_eq!(r.status, Status::Infeasible);
            assert!(r.farkas.is_some() && rose::verify::check(&m, &r, 1e-6).is_ok());
        }
        if let Some(m) = fixture("fx_unbounded") {
            let r = solve(&m, None, &o).unwrap();
            assert_eq!(r.status, Status::Unbounded);
            assert!(r.ray.is_some() && rose::verify::check(&m, &r, 1e-6).is_ok());
        }
        for (name, x, z) in [
            ("fx_neardep", vec![0.0, 4.0], -8.0),
            ("fx_badscale", vec![0.0, 5e-6], -5.0),
            ("fx_tied", vec![3.0, 1.0], -4.0),
        ] {
            if let Some(m) = fixture(name) {
                let r = solve(&m, None, &o).unwrap();
                assert_eq!(r.status, Status::Optimal, "{name}");
                assert!(
                    (r.objective - z).abs() <= 1e-9 * z.abs().max(1.0),
                    "{name}: {}",
                    r.objective
                );
                assert!(rose::verify::check(&m, &r, 1e-6).is_ok(), "{name}");
                if name != "fx_tied" {
                    assert!(
                        r.x.iter().zip(&x).all(|(a, b)| (a - b).abs() <= 1e-9),
                        "{name}: {:?}",
                        r.x
                    );
                }
            }
        }
        if let Some(m) = fixture("fx_warm") {
            for up in [false, true] {
                let b = Basis {
                    basic: vec![2, 1],
                    at_upper: vec![up, false, false, false],
                };
                let r = solve(&m, Some(&b), &o).unwrap();
                assert_eq!((r.status, r.x.clone()), (Status::Optimal, vec![3.0, 1.0]));
                assert!(rose::verify::check(&m, &r, 1e-6).is_ok());
            }
        }
        for (name, pairs) in [
            ("fx_singular", vec![(0, 0), (1, 1)]),
            ("fx_sing_dup_cross", vec![(0, 1), (1, 0)]),
            ("fx_sing_dup_attr", vec![(0, 0), (1, 1)]),
            ("fx_sing_dup_attr_ub", vec![(0, 0), (1, 1)]),
            ("fx_sing_dup_ub", vec![(0, 0), (1, 1)]),
            ("fx_sing_zero_col", vec![(0, 0), (1, 1)]),
            ("fx_sing_zero_col2", vec![(0, 0), (1, 1)]),
            ("fx_sing_slack_dup2", vec![(0, 1), (1, 0)]),
            ("fx_sing_rank3", vec![(0, 0), (1, 1), (2, 2)]),
            ("fx_sing_all_same3", vec![(0, 0), (1, 1), (2, 2)]),
            ("fx_sing_dense4", vec![(0, 0), (1, 1), (2, 2), (3, 3)]),
            ("fx_sing_two_pairs", vec![(0, 0), (1, 1), (2, 2), (3, 3)]),
        ] {
            let Some(m) = fixture(name) else { continue };
            let (n, rows) = (m.num_cols(), m.num_rows());
            let mut basic: Vec<usize> = (n..n + rows).collect();
            for &(j, r) in &pairs {
                basic[r] = j;
            }
            let b = Basis {
                basic,
                at_upper: vec![false; n + rows],
            };
            let r = solve(&m, Some(&b), &o).unwrap();
            assert_eq!(r.status, Status::Optimal, "{name}");
            assert!(rose::verify::check(&m, &r, 1e-6).is_ok(), "{name}");
            // The optimum is the same one the primal reaches from the same supplied basis (objective, not vertex: ties are legal).
            let p = solve(&m, Some(&b), &SolveOptions::default()).unwrap();
            assert!(
                (r.objective - p.objective).abs() <= 1e-9 * p.objective.abs().max(1.0),
                "{name}: {} vs {}",
                r.objective,
                p.objective
            );
        }
    }
}

/// Infeasible Netlib models (`netlib.org/lp/infeas`): every simplex profile ends `Infeasible` with a certificate that `verify::check` accepts
/// (the primal release profile once ended `Optimal` here). In the primal controlled profile a certificate that uses an inequality row has the reference's
/// signs (`y >= 0` on a row with only a lower side), as the reference's `FarkasDual` of `inf_galenet` ([0, 0, 0, 0, 1, 0, 1, 1]) does.
#[test]
fn infeasible_models_end_infeasible_in_every_simplex_profile() {
    for name in ["inf_galenet", "inf_itest6", "fx_infeasible"] {
        let path = format!("{}/tests/fixtures/{name}.mps", env!("CARGO_MANIFEST_DIR"));
        let m = match std::fs::read_to_string(&path) {
            Ok(t) => mps::parse(&t).unwrap(),
            Err(_) => fixture(name).unwrap(),
        };
        for (profile, o) in [
            ("primal", SolveOptions::default()),
            ("primal release", SolveOptions::release()),
            ("dual", SolveOptions::dual()),
            ("dual release", SolveOptions::dual_release()),
        ] {
            let r = solve(&m, None, &o).unwrap();
            assert_eq!(r.status, Status::Infeasible, "{name} {profile}");
            if name == "inf_itest6" && profile == "primal" {
                // The reference's controlled primal returns the cost-vector BTRAN as it is: `+1` on every violated slack whatever the row's
                // sense, which is not a proof on this model (a lower-only and an upper-only row both get `+1`) but is what rose reproduces.
                let want = [1.0, 0.0, 0.0, 1.0, 1.0, -2.15, 0.0, 0.0, -39.0, 40.0, 1.0];
                let y = r.farkas.as_ref().unwrap();
                assert!(
                    y.iter().zip(want).all(|(a, b)| (a - b).abs() < 1e-9),
                    "{y:?}"
                );
                continue;
            }
            assert!(
                rose::verify::check(&m, &r, 1e-6).is_ok(),
                "{name} {profile}"
            );
        }
    }
    let m = mps::parse(
        &std::fs::read_to_string(format!(
            "{}/tests/fixtures/inf_galenet.mps",
            env!("CARGO_MANIFEST_DIR")
        ))
        .unwrap(),
    )
    .unwrap();
    let y = solve(&m, None, &SolveOptions::default())
        .unwrap()
        .farkas
        .unwrap();
    assert_eq!(y.len(), 8);
    assert!(y[4] > 0.0 && y[6] > 0.0 && y[7] > 0.0, "{y:?}");
}

/// The dual clean-up pass after a shifted optimum (random LPs of `corpus/randlp/gen3.py`, controlled profile; reference run through `tools/ref_run.py`).
/// The pass is the dual simplex run from the primal basis (`crate::dual`) with the LU threshold of the primal's accuracy ladder; the primal then prices
/// again (`fx_cleanup_one`: one dual pivot, then a primal pivot, 16 iterations). x, reduced costs and duals of the returned solution are compared
/// bit for bit with the reference by `tools/compare.py`-style runs (`gen3` seeds 243, 247, 566, 680, 1296); here the iteration count and the objective.
/// `fx_cleanup_infeasible` (`gen` seed 2020): the dual proves the infeasibility after four pivots (22 iterations).
#[test]
fn dual_cleanup_pass_matches_the_reference() {
    for (name, iters, obj) in [
        ("fx_cleanup_one", 16, Some(9.999962571150556)),
        ("fx_cleanup_two", 11, Some(56.38623931869912)),
        ("fx_cleanup_four", 9, Some(0.0004169818070081059)),
        ("fx_cleanup_five", 25, Some(0.0005515042229220582)),
        ("fx_cleanup_six", 34, Some(1.2259691964239767)),
        ("fx_cleanup_infeasible", 22, None),
    ] {
        let Some(m) = fixture(name) else { continue };
        let opts = SolveOptions {
            faithful_weights: true,
            faithful_positions: true,
            faithful_ftran: true,
            ..SolveOptions::default()
        };
        let r = solve(&m, None, &opts).unwrap();
        assert_eq!(r.iterations, iters, "{name}");
        match obj {
            Some(o) => {
                assert_eq!(r.status, Status::Optimal, "{name}");
                assert_eq!(r.objective.to_bits(), f64::to_bits(o), "{name}");
            }
            None => assert_eq!(r.status, Status::Infeasible, "{name}"),
        }
    }
}

#[test]
fn infeasible_dual_cleanup_returns_its_reference_certificate() {
    let m = mps::parse(include_str!("fixtures/fx_cleanup_infeasible.mps")).unwrap();
    for (opts, want) in [
        (
            SolveOptions::default(),
            [0.0003322259136212624, 1.661129568106312],
        ),
        (SolveOptions::release(), [3.125e-06, 0.015625]),
    ] {
        let r = solve(&m, None, &opts).unwrap();
        assert_eq!(r.status, Status::Infeasible);
        let mut expected = vec![0.0; 12];
        expected[5] = want[0];
        expected[10] = want[1];
        assert_eq!(r.farkas, Some(expected));
    }
}

#[test]
fn phase_one_retries_preserve_pricing_after_bound_flips() {
    let model = mps::parse(include_str!("fixtures/fx_ladder_386.mps")).unwrap();
    let opts = SolveOptions {
        max_iterations: 24,
        ..SolveOptions::default()
    };
    let r = solve(&model, None, &opts).unwrap();
    let tail: Vec<_> = r
        .trace
        .iter()
        .skip(16)
        .map(|e| (e.entering, e.leaving))
        .collect();
    assert_eq!(
        tail,
        [
            (3, None),
            (40, None),
            (40, None),
            (17, Some(50)),
            (28, None),
            (40, None),
            (40, None),
            (36, Some(57))
        ]
    );
}
