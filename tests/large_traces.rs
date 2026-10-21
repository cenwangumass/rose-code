//! FIT2P and 80BAU3B (controlled profile) reproduce every reference decision and the returned solution bit for bit.
//! The models, reference traces and reference runs live in `reverse-engineer/`; the tests are ignored by default (about
//! 30 s each in the release build): `cargo test --release --test large_traces -- --ignored`.

use rose::{mps, solve, SolveOptions, Status};

fn check(name: &str) {
    let root = format!("{}/../reverse-engineer", env!("CARGO_MANIFEST_DIR"));
    let model = mps::parse(&std::fs::read_to_string(format!("{root}/corpus/netlib/mps/{name}.mps")).unwrap()).unwrap();
    let opts = SolveOptions { faithful_weights: true, faithful_positions: true, faithful_ftran: true, ..SolveOptions::default() };
    let r = solve(&model, None, &opts).unwrap();
    assert_eq!(r.status, Status::Optimal);

    // Reference decisions: rows with a pivot (`f2fc` is the entering variable); the leaving variable is the one whose new
    // position is -1 (absent for a bound flip or a slack leaving from its own row).
    let rows: Vec<serde_json::Value> = std::fs::read_to_string(format!("{root}/traces/{name}-controlled.jsonl"))
        .unwrap()
        .lines()
        .map(|l| serde_json::from_str(l).unwrap())
        .filter(|r: &serde_json::Value| r["f2fc"].as_i64().unwrap_or(-1) != -1)
        .collect();
    assert_eq!(r.trace.len(), rows.len(), "pivot count");
    for (k, (row, e)) in rows.iter().zip(&r.trace).enumerate() {
        let entering = row["f2fc"].as_u64().unwrap() as usize;
        assert_eq!(e.entering, entering, "row {k}: entering");
        let leaving = row["changed"].as_array().unwrap().iter().find_map(|c| {
            (c[2].as_i64() == Some(-1) && c[0].as_u64() != Some(entering as u64)).then(|| c[0].as_u64().unwrap() as usize)
        });
        if leaving.is_some() {
            assert_eq!(e.leaving, leaving, "row {k}: leaving");
        }
    }

    let run: serde_json::Value = serde_json::from_str(&std::fs::read_to_string(format!("{root}/runs/{name}-controlled.json")).unwrap()).unwrap();
    let res = &run["result"];
    let bits = |v: &serde_json::Value| v.as_array().unwrap().iter().map(|x| x.as_f64().unwrap().to_bits()).collect::<Vec<_>>();
    let mine = |v: &[f64]| v.iter().map(|x| x.to_bits()).collect::<Vec<_>>();
    assert_eq!(r.objective.to_bits(), res["obj"].as_f64().unwrap().to_bits(), "objective");
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
