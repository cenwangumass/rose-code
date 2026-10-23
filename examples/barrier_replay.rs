//! `ROSE_DUMP=barrier[,barrier_factor] barrier_replay model.mps [solves.json|dense]`: the barrier iteration on a model (`rose::barrier`), with the factor and
//! the solves replayed from captured outputs (`solves.json`, written by `reverse-engineer/tools/compare.py trace --profile barrier-controlled --rose`;
//! see `ReplaySolver::from_json`), by the plain f64 `DenseSolver` (`dense`, not bit-exact, tolerance only) or, without a second argument, by the
//! bit-exact `CholSolver` (`rose::chol`: the whole solve from the MPS alone). Records (stderr): BMODEL (the converted form), one BARRIER per iteration,
//! and BRES (final x, dual, reduced cost, objective, in the model's own coordinates) for the bit-exact solver.
use rose::barrier::{self, BarrierStatus, DenseSolver, Options, ReplaySolver};
use rose::dump::{emit_tag, on};
use rose::{mps, order};
use serde_json::json;

fn main() {
    let a: Vec<String> = std::env::args().collect();
    let model = mps::parse(&std::fs::read_to_string(&a[1]).unwrap()).unwrap();
    let s = match order::structure(&model, order::DENSE_COL_MIN) {
        Ok(s) => s,
        Err(e) => {
            println!("unsupported: {}", e.3);
            std::process::exit(2);
        }
    };
    // BARRIER_PROFILE=release: the release profile (scaled model, adaptive correctors; the whole solve from the MPS alone).
    let release = std::env::var("BARRIER_PROFILE").map_or(false, |v| v == "release");
    let sc = if release { Some(rose::scale::compute_with_ranges(&model, &rose::barrier::range_columns(&model))) } else { None };
    let f = match barrier::std_form_scaled(&model, &s.sym.perm, sc.as_ref()) {
        Ok(f) => f,
        Err(e) => {
            println!("{e}");
            std::process::exit(2);
        }
    };
    if on("barrier") {
        emit_tag("barrier", "BMODEL", json!({
            "nrow": f.m, "ncol": f.n2, "start": f.col_start, "idx": f.row_idx, "val": f.val, "sense": f.sense, "flags": f.flags,
            "b610": f.b, "b_file": f.b_file, "c": f.c, "ub": f.ub, "lo": f.ub_lo, "hi": f.ub_hi, "obj": f.offset, "const": f.cshift,
        }));
    }
    let mut o = Options { release, ..Options::default() };
    // PARAMS="K=V,...": the reference parameters the run is compared under (BarIterLimit).
    for kv in std::env::var("PARAMS").unwrap_or_default().split(',').filter(|t| !t.is_empty()) {
        if let Some(("BarIterLimit", v)) = kv.split_once('=') {
            o.max_iter = v.parse().expect("BarIterLimit");
        }
    }
    let mut solver = rose::chol::CholSolver::new(&f, &s.sym);
    let run = if a.get(2).map(String::as_str) == Some("dense") {
        barrier::run(&f, &mut DenseSolver::default(), &o)
    } else if let Some(p) = a.get(2) {
        let v: serde_json::Value = serde_json::from_str(&std::fs::read_to_string(p).unwrap()).unwrap();
        let mut rs = ReplaySolver::from_json(&v);
        let r = barrier::run(&f, &mut rs, &o);
        if r.is_ok() && !rs.consumed() {
            println!("replay: recorded solves left over");
        }
        r
    } else {
        barrier::run(&f, &mut solver, &o)
    };
    match run {
        Ok(r) => {
            if on("barrier") && a.get(2).is_none() {
                let sol = barrier::solution(&f, r.clone());
                let code = match r.status {
                    BarrierStatus::Optimal => 2,
                    BarrierStatus::IterationLimit => 7,
                    BarrierStatus::Other(c) => c,
                };
                emit_tag("barrier", "BRES", json!({"x": sol.x, "pi": sol.y, "rc": sol.rc, "obj": sol.objective, "iters": r.iterations, "status": code}));
            }
            let last = r.records.last().unwrap();
            let st = match r.status {
                BarrierStatus::Optimal => "optimal".to_string(),
                BarrierStatus::IterationLimit => "iteration limit".to_string(),
                BarrierStatus::Other(c) => format!("status {c}"),
            };
            println!("{st} after {} iterations, objective {:.9e}", r.iterations, last.res.q130);
        }
        Err(e) => {
            println!("{e}");
            std::process::exit(2);
        }
    }
}
