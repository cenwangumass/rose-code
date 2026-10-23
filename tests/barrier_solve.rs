//! The whole barrier solve from the MPS alone (`rose::barrier` with the bit-exact `rose::chol::CholSolver`) against the reference
//! (`tests/fixtures/barrier-solve-controlled.json`, written by `reverse-engineer/tools/make_barrier_solve_fixture.py`): every iteration's state and
//! decisions, and the returned x, duals, reduced costs and objective, agree bit for bit on AFIRO, SC50A, BLEND, ADLITTLE, BORE3D, DEGEN2, BOEING2 (ranged
//! rows, implied upper bounds dropped), PILOT4 (free columns, refinement, step caps, restart), DEGEN3 (scaled columns `Q+0x648`, supernodes cut to 192 columns, generic assembly of the root) and the two ranged-row fixtures (`fx_bar_rngL`, `fx_bar_rngGE`: L, G and E rows, R of both signs and 0).
use rose::barrier::{self, BarrierStatus, Options};
use rose::{mps, solve, SolveOptions, Status};
use serde_json::Value;

const MODELS: [(&str, &str); 11] = [
    ("afiro", include_str!("fixtures/afiro.mps")),
    ("sc50a", include_str!("fixtures/sc50a.mps")),
    ("blend", include_str!("fixtures/blend.mps")),
    ("adlittle", include_str!("fixtures/adlittle.mps")),
    ("bore3d", include_str!("fixtures/bore3d.mps")),
    ("degen2", include_str!("fixtures/degen2.mps")),
    ("boeing2", include_str!("fixtures/boeing2.mps")),
    ("fx_bar_rngL", include_str!("fixtures/fx_bar_rngL.mps")),
    ("fx_bar_rngGE", include_str!("fixtures/fx_bar_rngGE.mps")),
    ("pilot4", include_str!("fixtures/pilot4.mps")),
    ("degen3", include_str!("fixtures/degen3.mps")),
];

fn floats(v: &Value) -> Vec<f64> {
    v.as_array().unwrap().iter().map(|x| x.as_f64().unwrap()).collect()
}

/// FNV-1a 64 over the little-endian bytes (`hash_f64` of `make_barrier_solve_fixture.py`).
fn hash_f64(v: &[f64]) -> String {
    let mut h: u64 = 0xcbf29ce484222325;
    for x in v {
        for b in x.to_bits().to_le_bytes() {
            h = (h ^ b as u64).wrapping_mul(0x100000001b3);
        }
    }
    h.to_string()
}

fn same(name: &str, field: &str, it: i64, want: &[f64], got: &[f64]) {
    assert_eq!(want.len(), got.len(), "{name} it {it} {field}: length");
    for (j, (w, g)) in want.iter().zip(got).enumerate() {
        assert_eq!(w.to_bits(), g.to_bits(), "{name} it {it} {field}[{j}]: reference {w:e} rose {g:e}");
    }
}

fn scalar(name: &str, field: &str, it: i64, want: &Value, got: f64) {
    if let Some(w) = want.as_f64() {
        assert_eq!(w.to_bits(), got.to_bits(), "{name} it {it} {field}: reference {w:e} rose {got:e}");
    }
}

#[test]
fn whole_solves_match_reference_bitwise() {
    let fx: Value = serde_json::from_str(include_str!("fixtures/barrier-solve-controlled.json")).unwrap();
    for (name, text) in MODELS {
        let model = mps::parse(text).unwrap();
        let sol = barrier::solve_exact(&model, &Options::default()).unwrap_or_else(|e| panic!("{name}: {e}"));
        assert_eq!(sol.run.status, BarrierStatus::Optimal, "{name}");
        for want in fx[name]["expect"].as_array().unwrap() {
            let it = want["it"].as_i64().unwrap();
            let got = sol.run.records.iter().find(|r| r.it as i64 == it).unwrap_or_else(|| panic!("{name}: no record for iteration {it}"));
            for (f, v) in [("x", &got.state.x), ("z", &got.state.z), ("y", &got.state.y), ("u", &got.state.u), ("v", &got.state.v)] {
                if let Some(h) = want.get(f) {
                    assert_eq!(h.as_str().unwrap(), hash_f64(v), "{name} it {it} {f}: differs from the reference (compare.py trace --rose --solve shows where)");
                }
            }
            scalar(name, "mu", it, &want["mu"], got.res.q1f0);
            scalar(name, "pobj", it, &want["pobj"], got.res.q130);
            scalar(name, "dobj", it, &want["dobj"], got.res.q140);
            if let Some(s) = want["status"].as_i64() {
                assert_eq!(s as i32, got.status, "{name} it {it}: status");
            }
            if want.get("tp").is_some() {
                let d = got.dec.as_ref().unwrap_or_else(|| panic!("{name} it {it}: no decisions"));
                for (f, v) in [("ap", d.ap), ("ad", d.ad), ("aff", d.aff), ("target", d.target), ("tp", d.tp), ("td", d.td)] {
                    scalar(name, f, it, &want[f], v);
                }
                assert_eq!(want["ntry"].as_u64().unwrap() as usize, d.ntry, "{name} it {it}: ntry");
            }
        }
        let run = &fx[name]["run"];
        assert_eq!(run["iters"].as_i64().unwrap() as i32, sol.run.iterations, "{name} iterations");
        // x and rc cover the model's columns and the range columns of the ranged rows.
        same(name, "x (read-out)", 0, &floats(&run["x"]), &sol.x);
        same(name, "pi (read-out)", 0, &floats(&run["pi"]), &sol.y);
        same(name, "rc (read-out)", 0, &floats(&run["rc"]), &sol.rc);
        scalar(name, "objective", 0, &run["obj"], sol.objective);
    }
}

/// The public entry point: `Method::Barrier` returns the same solution in the `SolveResult` fields.
#[test]
fn method_barrier_returns_the_reference_solution() {
    let fx: Value = serde_json::from_str(include_str!("fixtures/barrier-solve-controlled.json")).unwrap();
    for (name, text) in MODELS {
        let model = mps::parse(text).unwrap();
        let r = solve(&model, None, &SolveOptions::barrier()).unwrap_or_else(|e| panic!("{name}: {e}"));
        assert_eq!(r.status, Status::Optimal, "{name}");
        let run = &fx[name]["run"];
        assert_eq!(run["iters"].as_u64().unwrap() as usize, r.iterations, "{name}");
        // The public vectors cover the model's own columns: the reference's x and rc also hold the range columns, after them.
        let (x, rc) = (floats(&run["x"]), floats(&run["rc"]));
        same(name, "x", 0, &x[..model.num_cols()], &r.x);
        same(name, "dual", 0, &floats(&run["pi"]), &r.dual);
        same(name, "reduced_cost", 0, &rc[..model.num_cols()], &r.reduced_cost);
        scalar(name, "objective", 0, &run["obj"], r.objective);
    }
}

#[test]
fn iteration_limit_is_reported() {
    let model = mps::parse(MODELS[0].1).unwrap();
    let r = solve(&model, None, &SolveOptions { max_iterations: 3, ..SolveOptions::barrier() }).unwrap();
    assert_eq!(r.status, Status::IterationLimit);
    assert_eq!(r.iterations, 3);
}

/// Features outside the recovered case are explicit errors, never a wrong answer.
#[test]
fn unsupported_models_are_explicit_errors() {
    let unsupported = |model: &rose::LpModel| matches!(solve(model, None, &SolveOptions::barrier()), Err(rose::simplex::SolveError::Unsupported(_)));
    let afiro = mps::parse(MODELS[0].1).unwrap();
    let mut max = afiro.clone();
    max.sense = rose::Sense::Maximize;
    assert!(unsupported(&max));
    // A row ranged twice (the later entry wins in the bounds, the reference's columns are not recovered for it).
    let twice = mps::parse(&include_str!("fixtures/fx_bar_rngL.mps").replace("RANGES\n", "RANGES\n    rng2      r1                 1\n")).unwrap();
    assert!(twice.range_dup);
    assert!(unsupported(&twice));
    // The release profile (scaling) is not recovered.
    assert!(matches!(solve(&afiro, None, &SolveOptions { scaling: true, ..SolveOptions::barrier() }), Err(rose::simplex::SolveError::Unsupported(_))));
}
