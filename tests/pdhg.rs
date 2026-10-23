//! PDHG (`rose::pdhg`) against the reference's recorded residual tests (`tests/fixtures/MODEL-pdhg-release.checks.json`, made by
//! `reverse-engineer/tools/make_pdhg_fixture.py`): counters, step size, primal weight, merit, objectives, residuals, status and
//! restart rule at every test must equal the reference's bit for bit; the final iteration count and objective too.

use rose::pdhg::{solve, PdhgOptions};
use rose::{mps, Status};
use serde_json::Value;

fn read(name: &str) -> String {
    std::fs::read_to_string(format!("{}/tests/fixtures/{name}", env!("CARGO_MANIFEST_DIR"))).unwrap()
}

fn check_model(name: &str) {
    check_with(name, PdhgOptions::default());
}

fn check_with(name: &str, opts: PdhgOptions) {
    let model = mps::parse(&read(&format!("{name}.mps"))).unwrap();
    let fix: Value = serde_json::from_str(&read(&format!("{name}-pdhg-release.checks.json"))).unwrap();
    let r = solve(&model, &PdhgOptions { record_checks: true, ..opts });
    let want_status = match fix["status"].as_f64() {
        Some(7.0) => Status::IterationLimit,
        Some(12.0) => Status::NumericalFailure,
        _ => Status::Optimal,
    };
    assert_eq!(r.status, want_status, "{name}");
    assert_eq!(r.iterations, fix["iterations"].as_u64().unwrap(), "{name}: iterations");
    if let Some(z) = fix["objective"].as_f64() {
        assert_eq!(r.objective, z, "{name}: objective");
    }
    let want = fix["checks"].as_array().unwrap();
    assert_eq!(r.checks.len(), want.len(), "{name}: number of tests");
    for (i, (got, exp)) in r.checks.iter().zip(want).enumerate() {
        for (key, e) in exp.as_object().unwrap() {
            let g = &got[key];
            let same = match (g.as_f64(), e.as_f64()) {
                (Some(a), Some(b)) => a == b,
                _ => g == e,
            };
            assert!(same, "{name}: test {} field {key}: rose {g}, reference {e}", i + 1);
        }
    }
}

#[test]
fn afiro_tests_equal_the_reference() {
    check_model("afiro");
}

#[test]
fn sc50a_tests_equal_the_reference() {
    check_model("sc50a");
}

#[test]
fn blend_tests_equal_the_reference() {
    check_model("blend");
}

#[test]
fn adlittle_tests_equal_the_reference() {
    check_model("adlittle");
}

#[test]
fn iteration_limit_stops_at_the_first_count_reaching_it() {
    // The limit counts the setup's power iterations too: AFIRO does 93 of them, so a limit of 300 stops at 300 (reference: status 7, 300).
    let model = mps::parse(&read("afiro.mps")).unwrap();
    let r = solve(&model, &PdhgOptions { iter_limit: 300.0, ..PdhgOptions::default() });
    assert_eq!(r.status, Status::IterationLimit);
    assert_eq!(r.iterations, 300);
}

#[test]
fn a_model_without_rows_is_solved_at_its_bounds() {
    let model = mps::parse("NAME T\nROWS\n N OBJ\nCOLUMNS\n X1 OBJ 1.0\n X2 OBJ -1.0\nBOUNDS\n UP BND X1 5.0\n UP BND X2 5.0\nENDATA\n").unwrap();
    let r = solve(&model, &PdhgOptions::default());
    assert_eq!(r.status, Status::Optimal);
    assert_eq!((r.x.clone(), r.objective, r.iterations), (vec![0.0, 5.0], -5.0, 0));
}

#[test]
fn nonzero_seed_perturbs_the_start_point() {
    // Reference (Seed=1, pdhg-release): AFIRO 602 iterations, SC50A 1562; the seed's value does not matter.
    for (name, iters, seed) in [("afiro", 602, 1), ("afiro", 602, 7), ("sc50a", 1562, 1)] {
        let model = mps::parse(&read(&format!("{name}.mps"))).unwrap();
        let r = solve(&model, &PdhgOptions { seed, ..PdhgOptions::default() });
        assert_eq!((r.status, r.iterations), (Status::Optimal, iters), "{name} seed {seed}");
    }
}

#[test]
fn solve_options_select_the_method() {
    let model = mps::parse(&read("afiro.mps")).unwrap();
    let r = rose::solve(&model, None, &rose::SolveOptions::pdhg()).unwrap();
    assert_eq!((r.status, r.iterations), (Status::Optimal, 593));
    assert!(rose::verify::check_optimal(&model, &r).max_error() < 1e-3);
}

#[test]
fn tight_tolerances_and_limit_exit_inside_the_power_iteration() {
    // pdhg_e3 at PDHGAbsTol = PDHGRelTol = PDHGConvTol = 1e-9, PDHGIterLimit 8000: restarts whose continued power iteration runs into the limit.
    check_with("pdhg_e3_tight", PdhgOptions { abs_tol: 1e-9, rel_tol: 1e-9, conv_tol: 1e-9, iter_limit: 8000.0, ..PdhgOptions::default() });
}

#[test]
fn row_sides_at_and_above_1e20_are_used_as_written() {
    // Right-hand sides 1e30 (pdhg_e17, never infinite for a row) and a range of 1e30 on an L row (pdhg_e22, a range column with no upper bound), limit 3000.
    for name in ["pdhg_e17", "pdhg_e22"] {
        check_with(name, PdhgOptions { iter_limit: 3000.0, ..PdhgOptions::default() });
    }
}

#[test]
fn zero_merit_takes_rule_5_and_the_zero_dual_residual_takes_the_weight_fallback() {
    // Reference with the state written by a debugger at one residual test (replay/gdb_pdhg_checks.py INJECT): merit R = 0 at n8 = 500 (restart rule 5, power-iteration
    // budget 1000) and relative dual residual = 0 at n8 = 700 (a restart with the primal weight reset to its initial value).
    check_with("adlittle_merit0", PdhgOptions { inject: Some((500, rose::pdhg::Inject::Merit)), ..PdhgOptions::default() });
    check_with("adlittle_reld0", PdhgOptions { inject: Some((700, rose::pdhg::Inject::RelDual)), ..PdhgOptions::default() });
}

#[test]
fn a_huge_objective_value_stops_with_status_12() {
    // Reference with the primal objective written to 1e101 at the residual test n8 = 300 (INJECT=pobj@300): the termination test returns 12.
    check_with("sc50a_pobj_inf", PdhgOptions { inject: Some((300, rose::pdhg::Inject::Pobj)), ..PdhgOptions::default() });
}

#[test]
fn a_long_run_takes_rule_5_and_the_weight_fallback_without_injection() {
    // pdhg_e3 to the limit 300000 (PDHGIterLimit): 271 restarts by rule 5 (merit 0 from a negative radicand) and 3 resets of the primal weight, all as in the reference.
    check_with("pdhg_e3_long", PdhgOptions { iter_limit: 300000.0, ..PdhgOptions::default() });
}
