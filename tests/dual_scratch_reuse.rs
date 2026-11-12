//! Retained scratch (result lists, L-stage buffers, level workspace, per-iteration glue lists) must not leak state between solves: a dual solve
//! that runs after solves of other dimensions on the same thread is bit-identical to the same solve on a fresh thread.

use rose::{mps, solve, SolveOptions};

type Print = (u64, Vec<u64>, Vec<u64>, Vec<u64>, usize, usize);

fn print(name: &str) -> Print {
    let text = std::fs::read_to_string(format!("{}/tests/fixtures/{name}.mps", env!("CARGO_MANIFEST_DIR"))).unwrap();
    let model = mps::parse(&text).unwrap();
    let r = solve(&model, None, &SolveOptions::dual_release()).unwrap();
    let bits = |v: &[f64]| v.iter().map(|x| x.to_bits()).collect::<Vec<_>>();
    (r.objective.to_bits(), bits(&r.x), bits(&r.dual), bits(&r.reduced_cost), r.iterations, r.refactorizations)
}

#[test]
fn consecutive_solves_of_different_dimensions_are_bit_identical() {
    let names = ["afiro", "sc50a", "adlittle", "boeing2", "degen2"];
    // Reference: each model on a thread of its own (empty retained scratch).
    let fresh: Vec<Print> = names
        .iter()
        .map(|&n| std::thread::spawn(move || print(n)).join().unwrap())
        .collect();
    // One thread, models of growing, shrinking and repeated sizes, so every buffer is reused at another dimension.
    let order = [4usize, 0, 2, 4, 3, 1, 0, 4, 1, 3, 2, 0];
    let seen = std::thread::spawn(move || order.iter().map(|&i| (i, print(names[i]))).collect::<Vec<_>>())
        .join()
        .unwrap();
    for (i, p) in seen {
        assert!(p == fresh[i], "{}: a solve after other solves differs from a fresh-thread solve", names[i]);
    }
}
