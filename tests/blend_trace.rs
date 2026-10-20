//! Partial decision-level agreement with the reference on BLEND (controlled profile) using the options that reproduce the
//! reference's basis position order after refactorization (`faithful_positions`, `refactor_schedule`) and its Devex weights.
//! The first 37 of 43 decisions agree; the remainder differs at a Phase 2 tie whose order depends on the reference's sparse
//! FTRAN output, which is not reproduced. This test pins the agreed prefix so regressions are visible.

use rose::{mps, solve, Basis, SolveOptions, Status};

const AGREED_PREFIX: usize = 37;

fn fixture(path: &str) -> String {
    std::fs::read_to_string(format!("{}/tests/fixtures/{path}", env!("CARGO_MANIFEST_DIR"))).unwrap()
}

#[test]
fn blend_agrees_on_the_recovered_prefix() {
    let model = mps::parse(&fixture("blend.mps")).unwrap();
    let (n, m) = (model.num_cols(), model.num_rows());
    let fx: serde_json::Value = serde_json::from_str(&fixture("blend-controlled.decisions.json")).unwrap();
    let mut basic = vec![usize::MAX; m];
    let mut at_upper = vec![false; n + m];
    for (j, s) in fx["status0"].as_array().unwrap().iter().enumerate() {
        let s = s.as_i64().unwrap();
        if s >= 0 {
            basic[s as usize] = j;
        } else if s == -2 || s == -5 {
            at_upper[j] = true;
        }
    }
    let schedule: Vec<usize> = fx["refactor_schedule"].as_array().unwrap().iter().map(|v| v.as_u64().unwrap() as usize).collect();
    let opts = SolveOptions { faithful_weights: true, faithful_positions: true, refactor_schedule: Some(schedule), ..SolveOptions::default() };
    let r = solve(&model, Some(&Basis { basic, at_upper }), &opts).unwrap();
    assert_eq!(r.status, Status::Optimal);
    assert!((r.objective - (-30.812149845828216)).abs() < 1e-7 + 1e-8 * 30.81, "objective {}", r.objective);
    let expected: Vec<(usize, Option<usize>)> =
        fx["decisions"].as_array().unwrap().iter().map(|d| (d[0].as_u64().unwrap() as usize, d[1].as_u64().map(|v| v as usize))).collect();
    let got: Vec<(usize, Option<usize>)> = r.trace.iter().map(|e| (e.entering, e.leaving)).collect();
    assert_eq!(&got[..AGREED_PREFIX], &expected[..AGREED_PREFIX], "decisions differ from the reference within the agreed prefix");
}
