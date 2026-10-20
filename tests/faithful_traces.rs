//! BLEND, SC50A and BORE3D (controlled profile) reproduce every reference decision when the reference's basis position order,
//! refactor iterations, Devex weights and sparse FTRAN output order are all mirrored (`faithful_*`).

use rose::{mps, solve, Basis, SolveOptions, Status};

fn fixture(path: &str) -> String {
    std::fs::read_to_string(format!("{}/tests/fixtures/{path}", env!("CARGO_MANIFEST_DIR"))).unwrap()
}

fn check(mps_name: &str, decisions: &str, objective: f64) {
    let model = mps::parse(&fixture(mps_name)).unwrap();
    let (n, m) = (model.num_cols(), model.num_rows());
    let fx: serde_json::Value = serde_json::from_str(&fixture(decisions)).unwrap();
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
    let opts = SolveOptions {
        faithful_weights: true,
        faithful_positions: true,
        faithful_ftran: true,
        refactor_schedule: Some(schedule),
        ..SolveOptions::default()
    };
    let r = solve(&model, Some(&Basis { basic, at_upper }), &opts).unwrap();
    assert_eq!(r.status, Status::Optimal);
    assert!((r.objective - objective).abs() < 1e-7 + 1e-8 * objective.abs(), "objective {}", r.objective);
    let expected: Vec<(usize, Option<usize>)> =
        fx["decisions"].as_array().unwrap().iter().map(|d| (d[0].as_u64().unwrap() as usize, d[1].as_u64().map(|v| v as usize))).collect();
    let got: Vec<(usize, Option<usize>)> = r.trace.iter().map(|e| (e.entering, e.leaving)).collect();
    assert_eq!(got, expected, "entering/leaving variables differ from the reference");
}

#[test]
fn blend_matches_all_reference_decisions() {
    check("blend.mps", "blend-controlled.decisions.json", -30.812149845828216);
}

#[test]
fn bore3d_matches_all_reference_decisions() {
    check("bore3d.mps", "bore3d-controlled.decisions.json", 1373.0803942085);
}

#[test]
fn sc50a_matches_all_reference_decisions() {
    check("sc50a.mps", "sc50a-controlled.decisions.json", -64.57507705856449);
}
