//! Compare solver decisions with a reference state file (JSON lines from reverse-engineer/tools/gdb_state.py).
//! Usage: cargo run --release --example compare -- model.mps states.jsonl
//! Reports the first divergence in (entering, leaving) variables, and whether the optimal objective agrees.

use rose::{mps, solve, Basis, SolveOptions};

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let model = mps::parse(&std::fs::read_to_string(&args[1]).unwrap()).unwrap();
    let (n, m) = (model.num_cols(), model.num_rows());
    let states: Vec<serde_json::Value> =
        std::fs::read_to_string(&args[2])
            .unwrap()
            .lines()
            .map(|l| serde_json::from_str(&l.replace("NaN", "null").replace("Infinity", "null")).unwrap())
            .collect();
    eprintln!("model {} x {}, state length {}", m, n, states[0]["status"].as_array().unwrap().len());
    let status: Vec<i64> = states[0]["status"].as_array().unwrap().iter().map(|s| s.as_i64().unwrap()).collect();
    let mut basic = vec![usize::MAX; m];
    let mut at_upper = vec![false; n + m];
    for (j, &s) in status.iter().enumerate() {
        if s >= 0 {
            basic[s as usize] = j;
        } else if s == -2 || s == -5 {
            at_upper[j] = true;
        }
    }
    let r = match solve(&model, Some(&Basis { basic, at_upper }), &SolveOptions::default()) {
        Ok(r) => r,
        Err(e) => {
            println!("{}: {e}", args[1]);
            return;
        }
    };
    let expected: Vec<(usize, Option<usize>)> = states
        .iter()
        .filter(|st| st["after"]["enter"].as_i64().unwrap() >= 0)
        .map(|st| {
            let enter = st["after"]["enter"].as_u64().unwrap() as usize;
            let pos = st["after"]["pos"].as_i64().unwrap();
            let sv: Vec<i64> = st["status"].as_array().unwrap().iter().map(|s| s.as_i64().unwrap()).collect();
            (enter, sv.iter().position(|&s| s == pos && pos >= 0))
        })
        .collect();
    let got: Vec<(usize, Option<usize>)> = r.trace.iter().map(|e| (e.entering, e.leaving)).collect();
    let agree = got.iter().zip(&expected).take_while(|(a, b)| a == b).count();
    println!(
        "{}: status {:?}, objective {:.10}, iterations {} (reference {}), agree on first {} decisions{}",
        args[1], r.status, r.objective, r.iterations, expected.len(), agree,
        if agree < got.len().min(expected.len()) {
            format!("; first difference: got {:?}, reference {:?}", got[agree], expected[agree])
        } else {
            String::new()
        }
    );
}
