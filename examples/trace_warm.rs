//! `trace_warm model.mps faithful|release|default basis.json`: like trace_dump, from a supplied start basis
//! ({"basic": [...], "at_upper": [...]}); then one `RESULT {json}` line with status, iterations, objective and x.
use rose::{mps, solve, Basis, SolveOptions};
fn main() {
    let a: Vec<String> = std::env::args().collect();
    let model = mps::parse(&std::fs::read_to_string(&a[1]).unwrap()).unwrap();
    let opts = if a[2] == "release" {
        SolveOptions::release()
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
    let v: serde_json::Value =
        serde_json::from_str(&std::fs::read_to_string(&a[3]).unwrap()).unwrap();
    let basic = v["basic"]
        .as_array()
        .unwrap()
        .iter()
        .map(|x| x.as_u64().unwrap() as usize)
        .collect();
    let at_upper = v["at_upper"]
        .as_array()
        .unwrap()
        .iter()
        .map(|x| x.as_bool().unwrap())
        .collect();
    let r = solve(&model, Some(&Basis { basic, at_upper }), &opts).unwrap();
    for e in &r.trace {
        println!(
            "{} {} {:e}",
            e.entering,
            e.leaving.map_or(-1, |v| v as i64),
            e.step
        );
    }
    println!(
        "RESULT {}",
        serde_json::json!({"status": format!("{:?}", r.status), "iterations": r.iterations, "objective": r.objective, "x": r.x})
    );
}
