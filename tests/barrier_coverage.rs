//! Diagnostic coverage of `Method::Barrier` (`SolveOptions::barrier()`, the controlled profile): the barrier's own interior point is checked by the
//! independent report of `rose::verify` (relative to the objective's magnitude), the optimum against the primal simplex, and the statuses of the
//! infeasible and unbounded fixtures. Models: reverse-engineer/corpus/fixtures/fx_*.mps (skipped if absent) and inline models.
//!
//! Unsupported or different from the primal result, on purpose:
//! - `fx_infeasible`: the barrier ends with its status 0xd (API status 12, infeasible or unbounded), reported as `NumericalFailure`; no certificate.
//! - `fx_badscale`: the reference's barrier stops at -4.99001868824 after 7 iterations (optimum -5); rose agrees with the reference, not with the optimum.
//! - an extra N row is ignored like the reference does (`unsupported_inputs_are_explicit`).

use rose::verify;
use rose::{mps, solve, LpModel, SolveOptions, Status};

fn load(path: String) -> Option<LpModel> {
    std::fs::read_to_string(path).ok().map(|t| mps::parse(&t).unwrap())
}

fn fixture(name: &str) -> Option<LpModel> {
    load(format!("{}/../reverse-engineer/corpus/fixtures/{name}.mps", env!("CARGO_MANIFEST_DIR")))
}

fn inline(text: &str) -> LpModel {
    mps::parse(text).unwrap()
}

fn barrier(m: &LpModel) -> rose::SolveResult {
    solve(m, None, &SolveOptions::barrier()).unwrap()
}

/// Optimal point: the independent report is small relative to the objective, and the objective equals the primal simplex's to the barrier's gap.
fn optimal_agrees(name: &str, m: &LpModel) -> rose::SolveResult {
    let r = barrier(m);
    assert_eq!(r.status, Status::Optimal, "{name}");
    let rep = verify::check_optimal(m, &r);
    let scale = 1.0 + r.objective.abs();
    assert!(rep.max_error() <= 1e-6 * scale.max(1e3), "{name}: {rep:?}");
    let p = solve(m, None, &SolveOptions::default()).unwrap();
    if p.status == Status::Optimal {
        assert!((r.objective - p.objective).abs() <= 1e-6 * scale, "{name}: barrier {} vs primal {}", r.objective, p.objective);
    }
    r
}

#[test]
fn fixtures_with_an_optimum() {
    // Degenerate, near-dependent, rank-deficient (`fx_sing_*`: duplicated rows and columns, empty columns), ranged, maximized, tied and warm-start models.
    let names = [
        "fx_bar_max", "fx_bar_nrow", "fx_bar_rngdup", "fx_bar_rngGE", "fx_bar_rngL", "fx_neardep", "fx_pdhg_max", "fx_singular", "fx_tied", "fx_warm",
        "fx_sing_all_same3", "fx_sing_dense4", "fx_sing_dense4n", "fx_sing_dup3_attr", "fx_sing_dup_attr", "fx_sing_dup_attr_ub", "fx_sing_dup_costs",
        "fx_sing_dup_cross", "fx_sing_dup_slack3", "fx_sing_dup_swapped", "fx_sing_dup_ub", "fx_sing_nuc_dep", "fx_sing_nuc_dep2", "fx_sing_rank3", "fx_sing_rank3b",
        "fx_sing_slack_dup", "fx_sing_slack_dup2", "fx_sing_two_deps", "fx_sing_two_pairs", "fx_sing_zero_col", "fx_sing_zero_col2",
    ];
    for n in names {
        if let Some(m) = fixture(n) {
            optimal_agrees(n, &m);
        }
    }
}

#[test]
fn badscale_matches_the_reference_stop() {
    let Some(m) = fixture("fx_badscale") else { return };
    let r = barrier(&m);
    assert_eq!(r.status, Status::Optimal);
    assert_eq!(r.iterations, 7);
    assert!((r.objective - -4.9900186882395445).abs() < 1e-12, "{}", r.objective);
}

#[test]
fn unbounded_fixtures() {
    for n in ["fx_unbounded", "fx_unbounded2", "fx_unbounded_p1"] {
        if let Some(m) = fixture(n) {
            let r = barrier(&m);
            assert_eq!(r.status, Status::Unbounded, "{n}");
            assert!(r.farkas.is_none() && r.ray.is_none(), "{n}: the barrier produces no certificate or ray");
        }
    }
}

#[test]
fn infeasible_fixture_is_reported_without_certificate() {
    if let Some(m) = fixture("fx_infeasible") {
        let r = barrier(&m);
        assert_eq!(r.status, Status::NumericalFailure);
        assert!(r.farkas.is_none());
    }
}

#[test]
fn alternate_optima_give_the_analytic_center() {
    // min x + y s.t. x + y = 2, 0 <= x, y <= 2: the whole segment is optimal; the barrier converges to its analytic center (1, 1), not a vertex.
    let m = inline("NAME T\nROWS\n N obj\n E c1\nCOLUMNS\n x obj 1 c1 1\n y obj 1 c1 1\nRHS\n rhs c1 2\nBOUNDS\n UP b x 2\n UP b y 2\nENDATA\n");
    let r = optimal_agrees("alt", &m);
    assert!((r.x[0] - 1.0).abs() < 1e-6 && (r.x[1] - 1.0).abs() < 1e-6, "{:?}", r.x);
}

#[test]
fn boxed_fixed_free_and_empty_columns_and_rows() {
    // fixed x (FX), boxed y, free z (FR) with a row that bounds it, an empty column w (cost 1, lower bound 2), a row without entries.
    let m = inline(
        "NAME T\nROWS\n N obj\n L c1\n G c2\n E c3\n L c4\nCOLUMNS\n x obj 1 c1 1\n x c2 1\n y obj 2 c1 1\n y c3 1\n z obj 1 c2 1\n z c3 -1\n w obj 1\nRHS\n rhs c1 10 c2 1\n rhs c3 0 c4 5\nBOUNDS\n FX b x 3\n UP b y 4\n FR b z\n LO b w 2\nENDATA\n",
    );
    let r = optimal_agrees("mixed", &m);
    assert!((r.x[0] - 3.0).abs() < 1e-6 && (r.x[3] - 2.0).abs() < 1e-6, "{:?}", r.x);
}

#[test]
fn iteration_limit() {
    let Some(m) = load(format!("{}/tests/fixtures/afiro.mps", env!("CARGO_MANIFEST_DIR"))) else { return };
    let r = solve(&m, None, &SolveOptions { max_iterations: 3, ..SolveOptions::barrier() }).unwrap();
    assert_eq!(r.status, Status::IterationLimit);
    assert_eq!(r.iterations, 3);
}

#[test]
fn unsupported_inputs_are_explicit() {
    let free_row = inline("NAME T\nROWS\n N obj\n N c1\n L c2\nCOLUMNS\n x obj 1 c1 1\n x c2 1\nRHS\n rhs c2 3\nENDATA\n");
    // an extra N row is ignored (the reference ignores it as well), so this one must solve
    assert_eq!(barrier(&free_row).status, Status::Optimal);
    let scaled = SolveOptions { scaling: true, ..SolveOptions::barrier() };
    // the release profile (scaling on) is supported
    let r = solve(&free_row, None, &scaled).unwrap();
    assert_eq!(r.status, Status::Optimal);
    assert!((r.objective - barrier(&free_row).objective).abs() < 1e-6);
}
