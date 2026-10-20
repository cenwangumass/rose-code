//! Print the first iterations' decisions: cargo run --release --example trace -- model.mps [max_iterations]
use rose::{mps, solve, SolveOptions};
fn main() {
    let a: Vec<String> = std::env::args().collect();
    let model = mps::parse(&std::fs::read_to_string(&a[1]).unwrap()).unwrap();
    let opts = SolveOptions { max_iterations: a.get(2).map(|v| v.parse().unwrap()).unwrap_or(100), ..SolveOptions::default() };
    let r = solve(&model, None, &opts).unwrap();
    for e in &r.trace {
        println!("{} enter {} leave {:?} pos {:?} step {:e}", e.iteration, e.entering, e.leaving, e.position, e.step);
    }
    println!("{:?} obj {}", r.status, r.objective);
}
