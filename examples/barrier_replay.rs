//! `ROSE_DUMP=barrier barrier_replay model.mps [solves.json]`: the barrier iteration on a model (`rose::barrier`), with the factor and the solves
//! replayed from captured outputs (`solves.json`, written by `reverse-engineer/tools/compare.py trace --profile barrier-controlled --rose`;
//! see `ReplaySolver::from_json`) or, without it, solved by the plain f64 `DenseSolver` (not bit-exact, tolerance only).
//! Test-only: the replay file is the only coupling to the analysis tooling. Records (stderr): BMODEL (the converted form), then one BARRIER per iteration.
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
    let f = match barrier::std_form(&model, &s.sym.perm) {
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
    let o = Options::default();
    let run = if let Some(p) = a.get(2) {
        let v: serde_json::Value = serde_json::from_str(&std::fs::read_to_string(p).unwrap()).unwrap();
        let mut rs = ReplaySolver::from_json(&v);
        let r = barrier::run(&f, &mut rs, &o);
        if r.is_ok() && !rs.consumed() {
            println!("replay: recorded solves left over");
        }
        r
    } else {
        barrier::run(&f, &mut DenseSolver::default(), &o)
    };
    match run {
        Ok(r) => {
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
