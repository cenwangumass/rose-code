//! `resolve_seq model.mps faithful|release|dual-release|default [basic_cols] [rounds]`: the public reoptimization sequence of the reference runner
//! `tools/witness_primal_r2_unbuffered.c`. Start basis: structurals `0..basic_cols` basic and the rest at the lower bound (0 = cold start);
//! solve with `c`, then with `-c`, `c`, ... each from the previous final basis. One `RESULT {json}` line per solve (stdout) after its
//! trace rows; the dump kinds (`ROSE_DUMP`) carry the records of every solve in order.
use rose::{mps, solve, Basis, SolveOptions};

fn main() {
    let a: Vec<String> = std::env::args().collect();
    let mut model = mps::parse(&std::fs::read_to_string(&a[1]).unwrap()).unwrap();
    let opts = if a[2] == "release" {
        SolveOptions::release()
    } else if a[2] == "dual-release" {
        SolveOptions::dual_release()
    } else if a[2] == "faithful" {
        SolveOptions {
            faithful_weights: true,
            faithful_positions: true,
            faithful_ftran: true,
            ..SolveOptions::default()
        }
    } else {
        SolveOptions::plain()
    };
    let (n, m) = (model.num_cols(), model.num_rows());
    let k: usize = a.get(3).map_or(0, |s| s.parse().unwrap());
    let rounds: usize = a.get(4).map_or(3, |s| s.parse().unwrap());
    let mut basis = if k == 0 {
        None
    } else {
        // Public VBasis 0 = basic, -1 = at lower; CBasis -1 = slack nonbasic at lower. The basis is the k structurals; the slacks
        // complete it only when k < m (never here: k == m).
        // `HEAD=j,...` (as for the reference runner) lists the basic variables instead: structurals below n, then slacks by row.
        let basic: Vec<usize> = match std::env::var("HEAD") {
            Ok(h) => h.split(',').map(|t| t.parse().unwrap()).collect(),
            Err(_) => {
                assert_eq!(k, m);
                (0..k).collect()
            }
        };
        Some(Basis {
            basic,
            at_upper: vec![false; n + m],
        })
    };
    for r in 0..rounds {
        if r > 0 {
            for c in model.obj.iter_mut() {
                *c = -*c;
            }
        }
        let res = solve(&model, basis.as_ref(), &opts).unwrap();
        for e in &res.trace {
            println!(
                "{} {} {:e}",
                e.entering,
                e.leaving.map_or(-1, |v| v as i64),
                e.step
            );
        }
        println!(
            "RESULT {}",
            serde_json::json!({"round": r, "status": format!("{:?}", res.status), "iterations": res.iterations, "objective": res.objective, "x": res.x, "pi": res.dual, "rc": res.reduced_cost})
        );
        basis = Some(res.basis);
    }
}
