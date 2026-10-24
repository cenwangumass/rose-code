//! The crash basis picks the same structural columns as the reference's start basis, and a cold start with the faithful
//! options then reproduces the reference decisions.

use rose::{crash_basis, mps, solve, SolveOptions, Status};

fn fixture(path: &str) -> String {
    std::fs::read_to_string(format!(
        "{}/tests/fixtures/{path}",
        env!("CARGO_MANIFEST_DIR")
    ))
    .unwrap()
}

fn check(mps_name: &str, decisions: &str, objective: f64) {
    let model = mps::parse(&fixture(mps_name)).unwrap();
    let n = model.num_cols();
    let fx: serde_json::Value = serde_json::from_str(&fixture(decisions)).unwrap();
    let mut want: Vec<usize> = fx["status0"]
        .as_array()
        .unwrap()
        .iter()
        .enumerate()
        .filter(|(_, s)| s.as_i64().unwrap() >= 0)
        .map(|(j, _)| j)
        .filter(|&j| j < n)
        .collect();
    let basis = crash_basis(&model);
    let mut got: Vec<usize> = basis.basic.iter().copied().filter(|&j| j < n).collect();
    want.sort();
    got.sort();
    assert_eq!(
        got, want,
        "crash columns differ from the reference start basis"
    );
    let opts = SolveOptions::default();
    let r = solve(&model, Some(&basis), &opts).unwrap();
    assert_eq!(r.status, Status::Optimal);
    assert!((r.objective - objective).abs() < 1e-7 + 1e-8 * objective.abs());
    let expected: Vec<(usize, Option<usize>)> = fx["decisions"]
        .as_array()
        .unwrap()
        .iter()
        .map(|d| {
            (
                d[0].as_u64().unwrap() as usize,
                d[1].as_u64().map(|v| v as usize),
            )
        })
        .collect();
    let got: Vec<(usize, Option<usize>)> =
        r.trace.iter().map(|e| (e.entering, e.leaving)).collect();
    assert_eq!(got, expected);
}

#[test]
fn sc50a() {
    check(
        "sc50a.mps",
        "sc50a-controlled.decisions.json",
        -64.57507705856449,
    );
}

#[test]
fn blend() {
    check(
        "blend.mps",
        "blend-controlled.decisions.json",
        -30.812149845828216,
    );
}

#[test]
fn adlittle() {
    check(
        "adlittle.mps",
        "adlittle-controlled.decisions.json",
        225494.96316238036,
    );
}

#[test]
fn bore3d() {
    check(
        "bore3d.mps",
        "bore3d-controlled.decisions.json",
        1373.0803942085,
    );
}
