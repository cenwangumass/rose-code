//! Compare the solver's decisions with the reference recorded inside the pinned binary (AFIRO, controlled profile).
//! `afiro-controlled.states.jsonl` holds the reference's state before each iteration (status per variable:
//! `>= 0` basis position, `-1`/`-4` at lower, `-2`/`-5` at upper) and its choice (`after.enter`, `after.pos`).
//! The starting basis is the reference's basis at its first iteration; decisions are compared as variables
//! (entering, leaving), because the reference renumbers basis positions at each refactorization.

use rose::{mps, solve, Basis, SolveOptions, Status};

fn fixture(path: &str) -> String {
    std::fs::read_to_string(format!(
        "{}/tests/fixtures/{path}",
        env!("CARGO_MANIFEST_DIR")
    ))
    .unwrap()
}

fn states() -> Vec<serde_json::Value> {
    fixture("afiro-controlled.states.jsonl")
        .lines()
        .map(|l| serde_json::from_str(l).unwrap())
        .collect()
}

fn status_of(v: &serde_json::Value) -> Vec<i64> {
    v["status"]
        .as_array()
        .unwrap()
        .iter()
        .map(|s| s.as_i64().unwrap())
        .collect()
}

fn start_basis(status: &[i64], n: usize, m: usize) -> Basis {
    let mut basic = vec![usize::MAX; m];
    let mut at_upper = vec![false; n + m];
    for (j, &s) in status.iter().enumerate() {
        if s >= 0 {
            basic[s as usize] = j;
        } else if s == -2 || s == -5 {
            at_upper[j] = true;
        }
    }
    Basis { basic, at_upper }
}

#[test]
fn afiro_matches_reference_decisions() {
    let model = mps::parse(&fixture("afiro.mps")).unwrap();
    let (n, m) = (model.num_cols(), model.num_rows());
    let states = states();
    let basis = start_basis(&status_of(&states[0]), n, m);
    let r = solve(&model, Some(&basis), &SolveOptions::default()).unwrap();
    assert_eq!(r.status, Status::Optimal);
    assert!(
        (r.objective - (-464.75314285714285)).abs() < 1e-7 + 1e-8 * 464.75,
        "objective {}",
        r.objective
    );

    let expected: Vec<(usize, usize)> = states
        .iter()
        .map(|st| {
            let enter = st["after"]["enter"].as_u64().unwrap() as usize;
            let pos = st["after"]["pos"].as_i64().unwrap();
            let leave = status_of(st).iter().position(|&s| s == pos).unwrap();
            (enter, leave)
        })
        .collect();
    let got: Vec<(usize, usize)> = r
        .trace
        .iter()
        .map(|e| (e.entering, e.leaving.unwrap()))
        .collect();
    assert_eq!(
        got, expected,
        "entering/leaving variables differ from the reference"
    );
}
