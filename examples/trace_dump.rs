//! `trace_dump model.mps faithful|release|dual|default [schedule]`: prints entering/leaving per iteration from a cold (crash) start.
//! `dual` runs the dual simplex (`SolveOptions::dual`); `schedule` is a comma list of iteration counts after which the work-based refactorization fires.
use rose::{mps, solve, SolveOptions};
fn main() {
    let a: Vec<String> = std::env::args().collect();
    let model = mps::parse(&std::fs::read_to_string(&a[1]).unwrap()).unwrap();
    let opts = if a[2] == "dual" || a[2] == "dual-controlled" || a[2] == "dual-release" {
        let mut o = if a[2] != "dual-release" { SolveOptions::dual() } else { SolveOptions::dual_release() };
        if let Some(sch) = a.get(3) {
            o.refactor_schedule = Some(sch.split(',').filter(|x| !x.is_empty()).map(|x| x.parse().unwrap()).collect());
        }
        if let Ok(v) = std::env::var("ROSE_MAXIT") {
            o.max_iterations = v.parse().unwrap();
        }
        o
    } else if a[2] == "release" { SolveOptions::release() } else if a[2] == "faithful" { SolveOptions { faithful_weights: true, faithful_positions: true, faithful_ftran: true, ..SolveOptions::default() } } else { SolveOptions::plain() };
    let r = solve(&model, None, &opts).unwrap();
    for e in &r.trace.iter().take(100000).collect::<Vec<_>>() {
        println!("{} {} {:e}", e.entering, e.leaving.map_or(-1, |v| v as i64), e.step);
    }
}
