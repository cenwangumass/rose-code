//! The dual simplex reproduces the reference's decisions: the leaving row (basis position) and the entering variable of every pivot of the
//! reference's dual simplex (`Method=1`, the controlled profile), from fixtures in `tests/fixtures/*-dual-controlled.decisions.json`
//! (written by `reverse-engineer/tools/make_dual_fixture.py` from GDB captures; the files are self-contained). The refactorizations follow the
//! reference's schedule (`refactor_schedule`), as its work-based trigger is not exact on every model. The bound flips of the fixture are not
//! compared (a `TraceEvent` carries none).

use rose::{mps, solve, SolveOptions, Status};

fn fixture(path: &str) -> String {
    std::fs::read_to_string(format!(
        "{}/tests/fixtures/{path}",
        env!("CARGO_MANIFEST_DIR")
    ))
    .unwrap()
}

fn check(name: &str) {
    let fx: serde_json::Value =
        serde_json::from_str(&fixture(&format!("{name}-dual-controlled.decisions.json"))).unwrap();
    assert_eq!(fx["version"], 1);
    let model = mps::parse(&fixture(fx["model"].as_str().unwrap())).unwrap();
    let schedule: Vec<usize> = fx["refactor_schedule"]
        .as_array()
        .unwrap()
        .iter()
        .map(|v| v.as_u64().unwrap() as usize)
        .collect();
    let opts = SolveOptions {
        refactor_schedule: Some(schedule),
        ..SolveOptions::dual()
    };
    let r = solve(&model, None, &opts).unwrap();
    assert_eq!(r.status, Status::Optimal);
    let z = fx["objective"].as_f64().unwrap();
    assert!(
        (r.objective - z).abs() <= 1e-7 + 1e-8 * z.abs(),
        "objective {} vs {z}",
        r.objective
    );
    let want: Vec<(usize, usize)> = fx["decisions"]
        .as_array()
        .unwrap()
        .iter()
        .map(|d| {
            (
                d[0].as_u64().unwrap() as usize,
                d[1].as_u64().unwrap() as usize,
            )
        })
        .collect();
    let got: Vec<(usize, usize)> = r
        .trace
        .iter()
        .map(|e| (e.position.unwrap(), e.entering))
        .collect();
    if let Some(k) = (0..want.len().min(got.len())).find(|&k| want[k] != got[k]) {
        panic!(
            "{name}: first difference at pivot {}: reference (row, entering) {:?}, rose {:?}",
            k + 1,
            want[k],
            got[k]
        );
    }
    assert_eq!(got.len(), want.len(), "{name}: number of pivots");
}

#[test]
fn afiro() {
    check("afiro");
}

#[test]
fn sc50a() {
    check("sc50a");
}

#[test]
fn blend() {
    check("blend");
}

#[test]
fn adlittle() {
    check("adlittle");
}
