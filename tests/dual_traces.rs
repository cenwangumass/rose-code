//! The dual simplex reproduces the reference's decisions: the leaving row (basis position) and the entering variable of every pivot of the
//! reference's dual simplex (`Method=1`; the controlled profile and, for the four small models, the release profile), from fixtures in `tests/fixtures/*-dual-{controlled,release}.decisions.json`
//! (written by `reverse-engineer/tools/make_dual_fixture.py` from GDB captures; the files are self-contained). The refactorizations follow the
//! reference's schedule (`refactor_schedule`), as its work-based trigger is not exact on every model. The bound flips of the fixture are
//! compared too, as a set per pivot.

use rose::{mps, solve, SolveOptions, Status};

fn fixture(path: &str) -> String {
    std::fs::read_to_string(format!(
        "{}/tests/fixtures/{path}",
        env!("CARGO_MANIFEST_DIR")
    ))
    .unwrap()
}

fn check(name: &str) {
    check_profile(name, "dual-controlled");
}

fn check_release(name: &str) {
    check_profile(name, "dual-release");
}

fn check_profile(name: &str, profile: &str) {
    let fx: serde_json::Value =
        serde_json::from_str(&fixture(&format!("{name}-{profile}.decisions.json"))).unwrap();
    assert_eq!(fx["version"], 1);
    let file = fx["model"].as_str().unwrap();
    let text = fixture(file);
    let model = mps::parse(&text).unwrap();
    let schedule: Vec<usize> = fx["refactor_schedule"]
        .as_array()
        .unwrap()
        .iter()
        .map(|v| v.as_u64().unwrap() as usize)
        .collect();
    let base = if profile == "dual-release" {
        SolveOptions::dual_release()
    } else {
        SolveOptions::dual()
    };
    let opts = SolveOptions {
        refactor_schedule: Some(schedule),
        ..base
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
    let sorted = |v: &Vec<usize>| {
        let mut v = v.clone();
        v.sort_unstable();
        v
    };
    let want_flips: Vec<Vec<usize>> = fx["decisions"]
        .as_array()
        .unwrap()
        .iter()
        .map(|d| {
            sorted(
                &d[2]
                    .as_array()
                    .unwrap()
                    .iter()
                    .map(|v| v.as_u64().unwrap() as usize)
                    .collect(),
            )
        })
        .collect();
    let got_flips: Vec<Vec<usize>> = r.trace.iter().map(|e| sorted(&e.flips)).collect();
    if let Some(k) = (0..want.len().min(got.len())).find(|&k| want[k] != got[k]) {
        panic!(
            "{name}: first difference at pivot {}: reference (row, entering) {:?}, rose {:?}",
            k + 1,
            want[k],
            got[k]
        );
    }
    assert_eq!(got.len(), want.len(), "{name}: number of pivots");
    if let Some(k) = (0..want.len()).find(|&k| want_flips[k] != got_flips[k]) {
        panic!(
            "{name}: flips differ at pivot {}: reference {:?}, rose {:?}",
            k + 1,
            want_flips[k],
            got_flips[k]
        );
    }
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

#[test]
fn afiro_release() {
    check_release("afiro");
}

#[test]
fn sc50a_release() {
    check_release("sc50a");
}

#[test]
fn blend_release() {
    check_release("blend");
}

#[test]
fn adlittle_release() {
    check_release("adlittle");
}

// The larger models: `cargo test --release --test dual_traces -- --ignored`.

#[test]
#[ignore]
fn large_bore3d() {
    check("bore3d");
}

#[test]
#[ignore]
fn large_bore3d_release() {
    check_release("bore3d");
}

#[test]
#[ignore]
fn large_degen2() {
    check("degen2");
}

#[test]
#[ignore]
fn large_degen2_release() {
    check_release("degen2");
}

#[test]
#[ignore]
fn large_boeing2() {
    check("boeing2");
}

#[test]
#[ignore]
fn large_boeing2_release() {
    check_release("boeing2");
}

#[test]
#[ignore]
fn large_pilot4() {
    check("pilot4");
}

#[test]
#[ignore]
fn large_pilot4_release() {
    check_release("pilot4");
}

#[test]
#[ignore]
fn large_degen3() {
    check("degen3");
}

#[test]
#[ignore]
fn large_degen3_release() {
    check_release("degen3");
}

#[test]
#[ignore]
fn large_perold() {
    check("perold");
}

#[test]
#[ignore]
fn large_perold_release() {
    check_release("perold");
}

#[test]
#[ignore]
fn large_80bau3b() {
    check("80bau3b");
}

#[test]
#[ignore]
fn large_80bau3b_release() {
    check_release("80bau3b");
}

#[test]
#[ignore]
fn large_fit2p() {
    check("fit2p");
}

#[test]
#[ignore]
fn large_fit2p_release() {
    check_release("fit2p");
}
