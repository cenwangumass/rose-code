//! `final_dump model.mps`: prints the faithful solve's objective, x, reduced costs, duals and iteration count as JSON (for comparison with `runs/*.json`).
use rose::{mps, solve, SolveOptions};
fn main() {
    let a: Vec<String> = std::env::args().collect();
    let model = mps::parse(&std::fs::read_to_string(&a[1]).unwrap()).unwrap();
    let opts = SolveOptions { faithful_weights: true, faithful_positions: true, faithful_ftran: true, ..SolveOptions::default() };
    let r = solve(&model, None, &opts).unwrap();
    println!("{}", serde_json::json!({"obj": r.objective, "x": r.x, "rc": r.reduced_cost, "pi": r.dual, "iters": r.iterations, "pivots": r.trace.len()}));
}
