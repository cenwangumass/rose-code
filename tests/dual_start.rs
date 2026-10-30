//! The public `DStart` start of the dual (`reverse-engineer/PLAN.md`, "Public PStart/DStart start"): `witness_dual_forced1`, `witness_dual_norow` (dual-release) and AFIRO
//! (dual-controlled) with a uniform DStart. `PStart` has no effect on the dual. The reference rows of `tests/fixtures/dual_start_reference.json` come from
//! `reverse-engineer/tools/witness_dual_start_runner.c` (status, IterCount, ObjVal, X, Pi, RC at 17 digits).
use rose::{mps, solve, SolveOptions, Status};

fn bits(v: &[f64]) -> Vec<u64> {
    v.iter().map(|x| x.to_bits()).collect()
}

fn floats(v: &serde_json::Value) -> Vec<f64> {
    v.as_array().unwrap().iter().map(|x| x.as_f64().unwrap()).collect()
}

#[test]
fn dstart_runs_agree_with_the_reference() {
    let reference: serde_json::Value = serde_json::from_str(include_str!("fixtures/dual_start_reference.json")).unwrap();
    let mut checked = 0;
    for case in reference.as_array().unwrap() {
        let name = case["model"].as_str().unwrap();
        let text = match name {
            "witness_dual_forced1" => include_str!("fixtures/witness_dual_forced1.mps"),
            "witness_dual_norow" => include_str!("fixtures/witness_dual_norow.mps"),
            "afiro" => include_str!("fixtures/afiro.mps"),
            other => panic!("{other}"),
        };
        let model = mps::parse(text).unwrap();
        let mut o = if case["profile"] == "release" { SolveOptions::dual_release() } else { SolveOptions::dual() };
        let y = case["dstart"].as_f64().unwrap();
        o.dstart = Some(vec![y; model.num_rows()]);
        let r = solve(&model, None, &o).unwrap();
        let tag = format!("{name} DStart {y}");
        assert_eq!(r.status, Status::Optimal, "{tag}");
        assert_eq!(r.iterations as u64, case["iters"].as_u64().unwrap(), "{tag}: iterations");
        assert_eq!(r.objective.to_bits(), case["obj"].as_f64().unwrap().to_bits(), "{tag}: objective");
        assert_eq!(bits(&r.x), bits(&floats(&case["x"])), "{tag}: x");
        assert_eq!(bits(&r.dual), bits(&floats(&case["pi"])), "{tag}: pi");
        assert_eq!(bits(&r.reduced_cost), bits(&floats(&case["rc"])), "{tag}: rc");
        checked += 1;
    }
    assert_eq!(checked, 12);
}

#[test]
fn zero_dstart_is_the_cold_start() {
    // A zero DStart gives d = c: no basic variable has a d, the stack is empty and the run is the cold one (AFIRO dual-controlled: 10 iterations).
    let model = mps::parse(include_str!("fixtures/afiro.mps")).unwrap();
    let cold = solve(&model, None, &SolveOptions::dual()).unwrap();
    let mut o = SolveOptions::dual();
    o.dstart = Some(vec![0.0; model.num_rows()]);
    let r = solve(&model, None, &o).unwrap();
    assert_eq!((r.status, r.iterations), (cold.status, cold.iterations));
    assert_eq!(bits(&r.x), bits(&cold.x));
}
