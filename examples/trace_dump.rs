//! `trace_dump model.mps faithful|release|default`: prints entering/leaving per iteration from a cold (crash) start.
use rose::{mps, solve, SolveOptions};
fn main() {
    let a: Vec<String> = std::env::args().collect();
    let model = mps::parse(&std::fs::read_to_string(&a[1]).unwrap()).unwrap();
    let opts = if a[2] == "release" { SolveOptions::release() } else if a[2] == "faithful" { SolveOptions { faithful_weights: true, faithful_positions: true, faithful_ftran: true, ..SolveOptions::default() } } else { SolveOptions::plain() };
    let r = solve(&model, None, &opts).unwrap();
    for e in &r.trace.iter().take(100000).collect::<Vec<_>>() {
        println!("{} {} {:e}", e.entering, e.leaving.map_or(-1, |v| v as i64), e.step);
    }
}
