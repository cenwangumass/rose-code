//! `trace_dump model.mps faithful|release|dual|default [schedule] [basis.json]`: prints entering/leaving per iteration from a cold
//! (crash) start, or from the supplied start basis ({"basic": [...], "at_upper": [...]}); the final `STATUS {json}` record (dump
//! kind `status`, stderr) holds the status, iterations and objective.
//! `dual` runs the dual simplex (`SolveOptions::dual`); `schedule` is a comma list of iteration counts after which the work-based refactorization fires.
use rose::{mps, solve, Basis, SolveOptions};
fn main() {
    let a: Vec<String> = std::env::args().collect();
    let model = mps::parse(&std::fs::read_to_string(&a[1]).unwrap()).unwrap();
    let basis_arg = a.iter().skip(3).find(|x| x.ends_with(".json"));
    let schedule_arg = a.iter().skip(3).find(|x| !x.ends_with(".json"));
    let opts = if a[2] == "dual" || a[2] == "dual-controlled" || a[2] == "dual-release" {
        let mut o = if a[2] != "dual-release" { SolveOptions::dual() } else { SolveOptions::dual_release() };
        if let Some(sch) = schedule_arg {
            o.refactor_schedule = Some(sch.split(',').filter(|x| !x.is_empty()).map(|x| x.parse().unwrap()).collect());
        }
        if let Ok(v) = std::env::var("ROSE_MAXIT") {
            o.max_iterations = v.parse().unwrap();
        }
        o
    } else if a[2] == "release" { SolveOptions::release() } else if a[2] == "faithful" { SolveOptions { faithful_weights: true, faithful_positions: true, faithful_ftran: true, ..SolveOptions::default() } } else { SolveOptions::plain() };
    let basis = basis_arg.map(|f| {
        let v: serde_json::Value = serde_json::from_str(&std::fs::read_to_string(f).unwrap()).unwrap();
        let basic = v["basic"].as_array().unwrap().iter().map(|x| x.as_u64().unwrap() as usize).collect();
        let at_upper = v["at_upper"].as_array().unwrap().iter().map(|x| x.as_bool().unwrap()).collect();
        Basis { basic, at_upper }
    });
    let r = solve(&model, basis.as_ref(), &opts).unwrap();
    for e in &r.trace.iter().take(100000).collect::<Vec<_>>() {
        println!("{} {} {:e}", e.entering, e.leaving.map_or(-1, |v| v as i64), e.step);
    }
    rose::dump::emit("status", serde_json::json!({"status": format!("{:?}", r.status), "iterations": r.iterations, "objective": r.objective}));
}
