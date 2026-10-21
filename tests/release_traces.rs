//! AFIRO (release profile, `SimplexPricing=-1` before its switch) reproduces every reference decision and the returned solution bit
//! for bit with the `pricing0` and `scaling` options; SC50A and BLEND also reproduce the method switch (`0x142500`)
//! and the Devex pricing after it. The models, reference traces and reference runs live in `reverse-engineer/`.

use rose::{mps, solve, SolveOptions, Status};

fn check(name: &str) {
    check_rows(name, usize::MAX);
}

/// The first `limit` decisions (all of them for `usize::MAX`, which also checks the returned solution).
fn check_rows(name: &str, limit: usize) {
    let root = format!("{}/../reverse-engineer", env!("CARGO_MANIFEST_DIR"));
    let model = mps::parse(&std::fs::read_to_string(format!("{root}/corpus/netlib/mps/{name}.mps")).unwrap()).unwrap();
    let r = solve(&model, None, &SolveOptions::release()).unwrap();
    if limit == usize::MAX {
        assert_eq!(r.status, Status::Optimal);
    }

    // Reference decisions: rows with a pivot (`f2fc` is the entering variable); the leaving variable is the one whose new
    // position is -1 (absent for a bound flip or a slack leaving from its own row).
    let rows: Vec<serde_json::Value> = std::fs::read_to_string(format!("{root}/traces/{name}-release.jsonl"))
        .unwrap()
        .lines()
        .map(|l| serde_json::from_str(l).unwrap())
        .filter(|r: &serde_json::Value| r["f2fc"].as_i64().unwrap_or(-1) != -1)
        .collect();
    if limit == usize::MAX {
        assert_eq!(r.trace.len(), rows.len(), "pivot count");
    }
    for (k, (row, e)) in rows.iter().zip(&r.trace).enumerate().take(limit) {
        let entering = row["f2fc"].as_u64().unwrap() as usize;
        assert_eq!(e.entering, entering, "row {k}: entering");
        let leaving = row["changed"].as_array().unwrap().iter().find_map(|c| {
            (c[2].as_i64() == Some(-1) && c[0].as_u64() != Some(entering as u64)).then(|| c[0].as_u64().unwrap() as usize)
        });
        if leaving.is_some() {
            assert_eq!(e.leaving, leaving, "row {k}: leaving");
        }
    }

    if limit != usize::MAX {
        return;
    }
    let run: serde_json::Value = serde_json::from_str(&std::fs::read_to_string(format!("{root}/runs/{name}-release.json")).unwrap()).unwrap();
    let res = &run["result"];
    let bits = |v: &serde_json::Value| v.as_array().unwrap().iter().map(|x| x.as_f64().unwrap().to_bits()).collect::<Vec<_>>();
    let mine = |v: &[f64]| v.iter().map(|x| x.to_bits()).collect::<Vec<_>>();
    assert_eq!(r.iterations as f64, res["iter_count"].as_f64().unwrap(), "iterations");
    assert_eq!(r.objective.to_bits(), res["obj"].as_f64().unwrap().to_bits(), "objective");
    assert_eq!(mine(&r.x), bits(&res["x"]), "x");
    assert_eq!(mine(&r.reduced_cost), bits(&res["rc"]), "reduced costs");
    assert_eq!(mine(&r.dual), bits(&res["pi"]), "duals");
}

#[test]
fn afiro_release_matches_the_reference() {
    check("afiro");
}

/// ADLITTLE starts in phase 1 (tolerance `P[0x3c0] * P[0]`, cost shifts in the dual update, dense dual recomputation at the end of phase 1)
/// and never switches the method: all 104 decisions and the returned solution agree.
#[test]
fn adlittle_release_matches_the_reference() {
    check("adlittle");
}

/// DEGEN2: 690 decisions (phase 1 first), the method switch at the twelfth factorization (step 526) and Devex pricing after it.
#[test]
fn degen2_release_matches_the_reference() {
    check("degen2");
}

/// SC50A: `SimplexPricing=0` for 25 decisions, the method switch at the sixth factorization (`0x142500`, `crate::method`) and Devex
/// pricing afterwards; all 32 decisions and the returned solution agree.
#[test]
fn sc50a_release_matches_the_reference() {
    check("sc50a");
}

/// BLEND: the switch at the third factorization; all 64 decisions and the returned solution agree.
#[test]
fn blend_release_matches_the_reference() {
    check("blend");
}

/// BORE3D: 101 decisions (phase 1 first, no switch). The candidate list's key is the one of the sort routine `0xfffc0`
/// (exponent clamped to -25..=25, key `0x19 - k` for the non-attractive), which puts scores of 2^24 and more with the basic variables.
#[test]
fn bore3d_release_matches_the_reference() {
    check("bore3d");
}

/// BOEING2 (ranges, scaling, shifted bounds): the start point agrees (bound placement of the crash with the crash columns, exact range
/// widths) and the first 120 decisions. A numeric difference in the last bit of a basic variable's distance at the 16th ratio
/// test (x_B accumulation) makes the next decision differ from row 120 on; not modelled.
#[test]
fn boeing2_release_matches_the_reference_at_the_start() {
    check_rows("boeing2", 120);
}

/// The scale factors against the reference's (`replay/gdb_scale_factors.py`, captured by `compare.py scale`; skipped where a
/// capture is missing, as the captures are not under version control).
#[test]
fn scale_factors_match_the_reference() {
    let root = format!("{}/../reverse-engineer", env!("CARGO_MANIFEST_DIR"));
    let mut checked = 0;
    for name in ["afiro", "sc50a", "blend", "adlittle", "bore3d", "boeing2", "degen2", "pilot4", "perold", "degen3"] {
        let Ok(text) = std::fs::read_to_string(format!("{root}/captures/{name}-release/scale_factors.LIM=1.jsonl")) else { continue };
        let line = text.lines().find(|l| l.starts_with("SCALE ")).unwrap();
        let rec: serde_json::Value = serde_json::from_str(&line[6..]).unwrap();
        let model = mps::parse(&std::fs::read_to_string(format!("{root}/corpus/netlib/mps/{name}.mps")).unwrap()).unwrap();
        let sc = rose::scale::compute(&model);
        let f = |v: &serde_json::Value| v.as_array().unwrap().iter().map(|x| x.as_f64().unwrap()).collect::<Vec<_>>();
        assert_eq!(sc.col, f(&rec["col"]), "{name}: column factors");
        assert_eq!(sc.row, f(&rec["row"]), "{name}: row factors");
        checked += 1;
    }
    eprintln!("{checked} models checked");
}
