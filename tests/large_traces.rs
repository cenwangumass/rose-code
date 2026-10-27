//! FIT2P and 80BAU3B (controlled profile) reproduce every reference decision and the returned solution bit for bit.
//! The models and compact reference fixtures are in `tests/fixtures/`; the tests are ignored by default (about
//! 30 s each in the release build): `cargo test --release --test large_traces -- --ignored`.

use rose::{mps, solve, SolveOptions, Status};

fn check(name: &str) {
    let root = format!("{}/tests/fixtures", env!("CARGO_MANIFEST_DIR"));
    let model =
        mps::parse(&std::fs::read_to_string(format!("{root}/{name}.mps")).unwrap()).unwrap();
    let opts = SolveOptions {
        faithful_weights: true,
        faithful_positions: true,
        faithful_ftran: true,
        ..SolveOptions::default()
    };
    let r = solve(&model, None, &opts).unwrap();
    assert_eq!(r.status, Status::Optimal);

    let fx: serde_json::Value = serde_json::from_str(
        &std::fs::read_to_string(format!("{root}/{name}-primal-controlled.json")).unwrap(),
    )
    .unwrap();
    assert_eq!(fx["version"], 1);
    let rows = fx["decisions"].as_array().unwrap();
    assert_eq!(r.trace.len(), rows.len(), "pivot count");
    for (k, (row, e)) in rows.iter().zip(&r.trace).enumerate() {
        assert_eq!(
            e.entering,
            row[0].as_u64().unwrap() as usize,
            "row {k}: entering"
        );
        if let Some(leaving) = row[1].as_u64() {
            assert_eq!(e.leaving, Some(leaving as usize), "row {k}: leaving");
        }
    }
    let res = &fx["result"];
    let bits = |v: &serde_json::Value| {
        v.as_array()
            .unwrap()
            .iter()
            .map(|x| x.as_f64().unwrap().to_bits())
            .collect::<Vec<_>>()
    };
    let mine = |v: &[f64]| v.iter().map(|x| x.to_bits()).collect::<Vec<_>>();
    assert_eq!(
        r.objective.to_bits(),
        res["obj"].as_f64().unwrap().to_bits(),
        "objective"
    );
    assert_eq!(mine(&r.x), bits(&res["x"]), "x");
    assert_eq!(mine(&r.reduced_cost), bits(&res["rc"]), "reduced costs");
    assert_eq!(mine(&r.dual), bits(&res["pi"]), "duals");
}

#[test]
#[ignore = "about 30 s"]
fn fit2p_matches_the_reference() {
    check("fit2p");
}

#[test]
#[ignore = "about 30 s"]
fn bau3b_matches_the_reference() {
    check("80bau3b");
}
