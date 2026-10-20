//! Decision-level agreement with the reference on ADLITTLE (controlled profile), which starts from an infeasible
//! basis and so exercises the recovered Phase 1 (big-M composite cost, typed breakpoint ratio test) and the
//! Phase 2 Devex/Harris rules. The start basis is the reference's basis at its first iteration.
//! The fixture stores the reference's `(entering, leaving)` variables per iteration (`leaving` is null for a bound flip).

use rose::{mps, solve, Basis, SolveOptions, Status};

fn fixture(path: &str) -> String {
    std::fs::read_to_string(format!("{}/tests/fixtures/{path}", env!("CARGO_MANIFEST_DIR"))).unwrap()
}

#[test]
fn adlittle_matches_all_reference_decisions() {
    let model = mps::parse(&fixture("adlittle.mps")).unwrap();
    let (n, m) = (model.num_cols(), model.num_rows());
    let fx: serde_json::Value = serde_json::from_str(&fixture("adlittle-controlled.decisions.json")).unwrap();
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
    let r = solve(&model, Some(&Basis { basic, at_upper }), &SolveOptions { faithful_weights: true, ..SolveOptions::default() }).unwrap();
    assert_eq!(r.status, Status::Optimal);
    assert!((r.objective - 225494.96316238036).abs() < 1e-7 + 1e-8 * 225494.96, "objective {}", r.objective);
    let expected: Vec<(usize, Option<usize>)> =
        fx["decisions"].as_array().unwrap().iter().map(|d| (d[0].as_u64().unwrap() as usize, d[1].as_u64().map(|v| v as usize))).collect();
    let got: Vec<(usize, Option<usize>)> = r.trace.iter().map(|e| (e.entering, e.leaving)).collect();
    assert_eq!(got, expected, "entering/leaving variables differ from the reference");
}

/// Same trace with the reference's basis position order after refactorization, its refactor iterations and the order of the
/// nonzeros its sparse FTRAN returns (`faithful_ftran`): the Phase 1 ratio test is sensitive to that order.
#[test]
fn adlittle_matches_with_faithful_positions_and_ftran_order() {
    let model = mps::parse(&fixture("adlittle.mps")).unwrap();
    let (n, m) = (model.num_cols(), model.num_rows());
    let fx: serde_json::Value = serde_json::from_str(&fixture("adlittle-controlled.decisions.json")).unwrap();
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
    let expected: Vec<(usize, Option<usize>)> =
        fx["decisions"].as_array().unwrap().iter().map(|d| (d[0].as_u64().unwrap() as usize, d[1].as_u64().map(|v| v as usize))).collect();
    let got: Vec<(usize, Option<usize>)> = r.trace.iter().map(|e| (e.entering, e.leaving)).collect();
    assert_eq!(got, expected);
}
