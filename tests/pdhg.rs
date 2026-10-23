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
    let model = mps::parse(&read(&format!("{name}.mps"))).unwrap();
    let fix: Value = serde_json::from_str(&read(&format!("{name}-pdhg-release.checks.json"))).unwrap();
    let r = solve(&model, &PdhgOptions { record_checks: true, ..PdhgOptions::default() });
    assert_eq!(r.status, Status::Optimal, "{name}");
    assert_eq!(r.iterations, fix["iterations"].as_u64().unwrap(), "{name}: iterations");
    assert_eq!(r.objective, fix["objective"].as_f64().unwrap(), "{name}: objective");
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
