//! Times a cold-start solve with and without the FTRAN mirror: `faithful_time model.mps`.
use rose::{crash_basis, mps, solve, SolveOptions};
fn main() {
    let model = mps::parse(&std::fs::read_to_string(std::env::args().nth(1).unwrap()).unwrap()).unwrap();
    for (name, opts) in [
        ("default", SolveOptions::default()),
        ("faithful", SolveOptions { faithful_weights: true, faithful_positions: true, faithful_ftran: true, ..SolveOptions::default() }),
        ("faithful+crash", SolveOptions { faithful_weights: true, faithful_positions: true, faithful_ftran: true, ..SolveOptions::default() }),
        ("default+crash", SolveOptions::default()),
    ] {
        let crash = crash_basis(&model);
        let t = std::time::Instant::now();
        match solve(&model, if name.ends_with("+crash") { Some(&crash) } else { None }, &opts) {
            Ok(r) => println!("{name}: {:?} obj {:.6} iters {} in {:.2?}", r.status, r.objective, r.iterations, t.elapsed()),
            Err(e) => println!("{name}: error {e}"),
        }
    }
}
