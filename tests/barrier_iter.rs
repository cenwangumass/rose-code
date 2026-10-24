//! Barrier iteration (`rose::barrier`) against the reference captures (`tests/fixtures/barrier-iter-controlled.json`, written by
//! `reverse-engineer/tools/make_barrier_fixture.py`): with the captured solve outputs replayed, every iteration's state and decisions
//! agree bit for bit (AFIRO and SC50A whole solves; the first iterations of BLEND, ADLITTLE and BORE3D, whose upper-bound pairs, lower-bound
//! shifts and fixed column it exercises). The plain f64 `DenseSolver` is checked at tolerance only (it is not bit-exact).
use rose::barrier::{self, BarrierStatus, Options, ReplaySolver, StdForm};
use rose::{mps, order};
use serde_json::Value;

const MODELS: [(&str, &str); 5] = [
    ("afiro", include_str!("fixtures/afiro.mps")),
    ("sc50a", include_str!("fixtures/sc50a.mps")),
    ("blend", include_str!("fixtures/blend.mps")),
    ("adlittle", include_str!("fixtures/adlittle.mps")),
    ("bore3d", include_str!("fixtures/bore3d.mps")),
];

fn floats(v: &Value) -> Vec<f64> {
    v.as_array()
        .unwrap()
        .iter()
        .map(|x| x.as_f64().unwrap())
        .collect()
}

fn std_form(text: &str) -> StdForm {
    let model = mps::parse(text).unwrap();
    let s = order::structure(&model, order::DENSE_COL_MIN).unwrap();
    barrier::std_form(&model, &s.sym.perm).unwrap()
}

fn same(name: &str, field: &str, it: i64, want: &[f64], got: &[f64]) {
    assert_eq!(want.len(), got.len(), "{name} it {it} {field}: length");
    for (j, (w, g)) in want.iter().zip(got).enumerate() {
        assert_eq!(
            w.to_bits(),
            g.to_bits(),
            "{name} it {it} {field}[{j}]: reference {w:e} rose {g:e}"
        );
    }
}

fn scalar(name: &str, field: &str, it: i64, want: &Value, got: f64) {
    if let Some(w) = want.as_f64() {
        assert_eq!(
            w.to_bits(),
            got.to_bits(),
            "{name} it {it} {field}: reference {w:e} rose {got:e}"
        );
    }
}

#[test]
fn replayed_iterations_match_reference_bitwise() {
    let fx: Value =
        serde_json::from_str(include_str!("fixtures/barrier-iter-controlled.json")).unwrap();
    for (name, text) in MODELS {
        let f = std_form(text);
        let expect = fx[name]["expect"].as_array().unwrap();
        let whole = fx[name]["whole"].as_bool().unwrap();
        let mut rs = ReplaySolver::from_json(&fx[name]["solves"]);
        // A truncated fixture stops at the iteration limit (the last expected iteration + 1); the iterations before it are compared.
        let niter = fx[name]["solves"]["iters"].as_array().unwrap().len() as i32;
        let opts = Options {
            max_iter: if whole { 1000 } else { niter },
            ..Options::default()
        };
        let run = barrier::run(&f, &mut rs, &opts).unwrap_or_else(|e| panic!("{name}: {e}"));
        assert!(rs.consumed(), "{name}: recorded solves left over");
        let mut checked = 0;
        for want in expect {
            let it = want["it"].as_i64().unwrap();
            let got = run
                .records
                .iter()
                .find(|r| r.it as i64 == it)
                .unwrap_or_else(|| panic!("{name}: no record for iteration {it}"));
            same(name, "x", it, &floats(&want["x"]), &got.state.x);
            same(name, "z", it, &floats(&want["z"]), &got.state.z);
            same(name, "y", it, &floats(&want["y"]), &got.state.y);
            if want.get("u").is_some() {
                same(name, "u", it, &floats(&want["u"]), &got.state.u);
                same(name, "v", it, &floats(&want["v"]), &got.state.v);
            }
            scalar(name, "mu", it, &want["mu"], got.res.q1f0);
            scalar(name, "pobj", it, &want["pobj"], got.res.q130);
            scalar(name, "dobj", it, &want["dobj"], got.res.q140);
            if let Some(s) = want["status"].as_i64() {
                assert_eq!(s as i32, got.status, "{name} it {it}: status");
            }
            if want.get("tp").is_some() {
                let d = got
                    .dec
                    .as_ref()
                    .unwrap_or_else(|| panic!("{name} it {it}: no decisions"));
                scalar(name, "ap", it, &want["ap"], d.ap);
                scalar(name, "ad", it, &want["ad"], d.ad);
                scalar(name, "aff", it, &want["aff"], d.aff);
                scalar(name, "target", it, &want["target"], d.target);
                scalar(name, "tp", it, &want["tp"], d.tp);
                scalar(name, "td", it, &want["td"], d.td);
                assert_eq!(
                    want["ntry"].as_u64().unwrap() as usize,
                    d.ntry,
                    "{name} it {it}: ntry"
                );
            }
            checked += 1;
        }
        assert_eq!(checked, expect.len());
        if whole {
            assert_eq!(run.status, BarrierStatus::Optimal, "{name}");
        }
    }
}

/// The reference's log: iterations and final objective of `runs/<model>-barrier-controlled.json`.
#[test]
fn dense_solver_reaches_the_reference_optimum() {
    for (name, text, iters, obj) in [
        ("afiro", MODELS[0].1, 8, -464.75314285176614),
        ("sc50a", MODELS[1].1, 8, -64.57507705583878),
        ("blend", MODELS[2].1, 10, -30.81214984582793),
        ("adlittle", MODELS[3].1, 12, 225494.96316238563),
        ("bore3d", MODELS[4].1, 20, 1373.0803942085195),
    ] {
        let model = mps::parse(text).unwrap();
        let sol = barrier::solve_dense(&model).unwrap_or_else(|e| panic!("{name}: {e}"));
        assert_eq!(sol.run.status, BarrierStatus::Optimal, "{name}");
        assert_eq!(sol.run.iterations, iters, "{name} iterations");
        assert!(
            (sol.objective - obj).abs() <= 1e-9 * obj.abs().max(1.0),
            "{name}: {} vs {obj}",
            sol.objective
        );
        // The returned point is feasible for the rows (file order, original columns).
        let ax = model.activities(&sol.x);
        for i in 0..model.num_rows() {
            let tol = 1e-5 * (1.0 + ax[i].abs());
            assert!(
                ax[i] >= model.row_lo[i] - tol && ax[i] <= model.row_up[i] + tol,
                "{name} row {i}: {} not in [{}, {}]",
                ax[i],
                model.row_lo[i],
                model.row_up[i]
            );
        }
        for j in 0..model.num_cols() {
            assert!(
                sol.x[j] >= model.col_lo[j] - 1e-6 && sol.x[j] <= model.col_up[j] + 1e-6,
                "{name} col {j}"
            );
        }
    }
}

#[test]
fn converted_models_are_accepted() {
    let model = mps::parse(include_str!("fixtures/afiro.mps")).unwrap();
    let mut max = model.clone();
    max.sense = rose::Sense::Maximize;
    let s = order::structure(&model, order::DENSE_COL_MIN).unwrap();
    // Maximization is converted (c and the offset negated; the solution negates the objective, y and rc back).
    assert!(barrier::std_form(&max, &s.sym.perm).is_ok());
    // A free column is converted (implied bounds, flag 0 when none is tight enough).
    let mut free = model.clone();
    free.col_lo[0] = -1e30;
    assert!(barrier::std_form(&free, &s.sym.perm).is_ok());
}
