//! `rose solve model.mps --output result.json [--basis-in b.json] [--basis-out b.json] [--method primal|dual|pdhg] [--profile controlled|release]`
//! Progress goes to stderr; stdout stays empty so machine-readable output is only ever in the files.

use rose::verify;
use rose::{mps, solve, Basis, SolveOptions, Status};
use serde_json::{json, Value};

fn usage() -> ! {
    eprintln!("usage: rose solve <model.mps> [--output result.json] [--basis-in basis.json] [--basis-out basis.json] [--max-iterations N] [--method primal|dual|pdhg|barrier] [--profile controlled|release] [--pdhg-iter-limit N] [--pdhg-rel-tol T] [--pdhg-abs-tol T] [--pdhg-conv-tol T] [--pdhg-seed N] [--pdhg-inject merit@N|reld@N|pobj@N]");
    std::process::exit(2)
}

fn basis_to_json(b: &Basis) -> Value {
    json!({ "version": 1, "basic": b.basic, "at_upper": b.at_upper })
}

fn basis_from_json(v: &Value) -> Option<Basis> {
    let basic = v["basic"].as_array()?.iter().map(|x| x.as_u64().map(|u| u as usize)).collect::<Option<Vec<_>>>()?;
    let at_upper = v["at_upper"].as_array()?.iter().map(|x| x.as_bool()).collect::<Option<Vec<_>>>()?;
    Some(Basis { basic, at_upper })
}

fn main() {
    let args: Vec<String> = std::env::args().skip(1).collect();
    if args.first().map(String::as_str) != Some("solve") || args.len() < 2 {
        usage();
    }
    let path = &args[1];
    let (mut output, mut basis_in, mut basis_out) = (None, None, None);
    let (mut release, mut max_it) = (false, None);
    let mut method = String::from("primal");
    let mut popts = rose::pdhg::PdhgOptions::default();
    let mut it = args[2..].iter();
    while let Some(a) = it.next() {
        let val = it.next().unwrap_or_else(|| usage());
        match a.as_str() {
            "--output" => output = Some(val.clone()),
            "--basis-in" => basis_in = Some(val.clone()),
            "--basis-out" => basis_out = Some(val.clone()),
            "--max-iterations" => max_it = Some(val.parse().unwrap_or_else(|_| usage())),
            "--method" => method = val.clone(),
            "--profile" => {
                release = match val.as_str() {
                    "controlled" => false,
                    "release" => true,
                    _ => usage(),
                }
            }
            "--pdhg-iter-limit" => popts.iter_limit = val.parse().unwrap_or_else(|_| usage()),
            "--pdhg-rel-tol" => popts.rel_tol = val.parse().unwrap_or_else(|_| usage()),
            "--pdhg-abs-tol" => popts.abs_tol = val.parse().unwrap_or_else(|_| usage()),
            "--pdhg-inject" => {
                let (w, at) = val.split_once('@').unwrap_or_else(|| usage());
                let at = at.parse().unwrap_or_else(|_| usage());
                popts.inject = Some((at, match w {
                    "merit" => rose::pdhg::Inject::Merit,
                    "reld" => rose::pdhg::Inject::RelDual,
                    "pobj" => rose::pdhg::Inject::Pobj,
                    _ => usage(),
                }));
            }
            "--pdhg-seed" => popts.seed = val.parse().unwrap_or_else(|_| usage()),
            "--pdhg-conv-tol" => popts.conv_tol = val.parse().unwrap_or_else(|_| usage()),
            _ => usage(),
        }
    }
    let dual = match method.as_str() {
        "primal" | "pdhg" | "barrier" => false,
        "dual" => true,
        _ => usage(),
    };
    let base = match (dual, release) {
        (false, false) => SolveOptions::default(),
        (false, true) => SolveOptions::release(),
        (true, false) => SolveOptions::dual(),
        (true, true) => SolveOptions::dual_release(),
    };
    let mut opts = SolveOptions { record_trace: false, ..base };
    if let Some(n) = max_it {
        opts.max_iterations = n;
    }
    let text = std::fs::read_to_string(path).unwrap_or_else(|e| fail(&format!("cannot read {path}: {e}")));
    let model = mps::parse(&text).unwrap_or_else(|e| fail(&format!("invalid input: {e}")));
    if method == "pdhg" {
        return run_pdhg(&model, &popts, output);
    } else if method == "barrier" {
        // The barrier's own point (no crossover): the controlled profile, BarIterLimit 1000 unless --max-iterations is given.
        let max_iterations = if args.iter().any(|a| a == "--max-iterations") { opts.max_iterations } else { 1000 };
        opts = SolveOptions { record_trace: false, max_iterations, ..SolveOptions::barrier() };
        if basis_in.is_some() || basis_out.is_some() {
            fail("--method barrier has no basis (Crossover=0): --basis-in and --basis-out are not supported");
        }
    }
    let start = basis_in.map(|p| {
        let v: Value = serde_json::from_str(&std::fs::read_to_string(&p).unwrap_or_else(|e| fail(&format!("cannot read {p}: {e}")))).unwrap_or_else(|e| fail(&format!("invalid basis file: {e}")));
        basis_from_json(&v).unwrap_or_else(|| fail("invalid basis file"))
    });
    eprintln!("{}: {} rows, {} columns", model.name, model.num_rows(), model.num_cols());
    let t = std::time::Instant::now();
    let r = solve(&model, start.as_ref(), &opts).unwrap_or_else(|e| match e {
        rose::simplex::SolveError::Unsupported(_) => fail(&e.to_string()),
        _ => fail(&format!("invalid input: {e}")),
    });
    let secs = t.elapsed().as_secs_f64();
    let status = match r.status {
        Status::Optimal => "optimal",
        Status::Infeasible => "infeasible",
        Status::Unbounded => "unbounded",
        Status::IterationLimit => "iteration_limit",
        Status::NumericalFailure => "numerical_failure",
    };
    let check = verify::check(&model, &r, 1e-6);
    eprintln!("status {status}, objective {:.12e}, {} iterations, {} refactorizations, {secs:.3}s", r.objective, r.iterations, r.refactorizations);
    if let Err(e) = &check {
        eprintln!("independent check failed: {e}");
    }
    let out = json!({
        "version": 1,
        "status": status,
        "objective": if r.status == Status::Optimal { json!(r.objective) } else { Value::Null },
        // like the reference, which returns no solution vectors unless the status is optimal; `PdhgResult` keeps the last iterate for library users
        "x": if r.status == Status::Optimal { json!(r.x) } else { Value::Null },
        "dual": if r.status == Status::Optimal { json!(r.dual) } else { Value::Null },
        "reduced_cost": if r.status == Status::Optimal { json!(r.reduced_cost) } else { Value::Null },
        "basis": basis_to_json(&r.basis),
        "farkas": r.farkas, "ray": r.ray,
        "diagnostics": { "iterations": r.iterations, "refactorizations": r.refactorizations, "seconds": secs, "check": check.err() },
    });
    if let Some(p) = basis_out {
        std::fs::write(&p, serde_json::to_string(&basis_to_json(&r.basis)).unwrap()).unwrap_or_else(|e| fail(&format!("cannot write {p}: {e}")));
    }
    let text = serde_json::to_string_pretty(&out).unwrap();
    match output {
        Some(p) => std::fs::write(&p, text).unwrap_or_else(|e| fail(&format!("cannot write {p}: {e}"))),
        None => println!("{text}"),
    }
    std::process::exit(match r.status {
        Status::Optimal => 0,
        Status::Infeasible | Status::Unbounded => 0,
        _ => 3,
    });
}

fn fail(msg: &str) -> ! {
    eprintln!("error: {msg}");
    std::process::exit(1)
}

/// `--method pdhg`: the reflected restarted Halpern PDHG (`rose::pdhg`); no basis.
fn run_pdhg(model: &rose::LpModel, popts: &rose::pdhg::PdhgOptions, output: Option<String>) {
    eprintln!("{}: {} rows, {} columns", model.name, model.num_rows(), model.num_cols());
    let t = std::time::Instant::now();
    let r = rose::pdhg::solve(model, popts);
    let secs = t.elapsed().as_secs_f64();
    let status = match r.status {
        Status::Optimal => "optimal",
        Status::IterationLimit => "iteration_limit",
        Status::Infeasible => "infeasible",
        Status::Unbounded => "unbounded",
        _ => "numerical_failure",
    };
    eprintln!("status {status}, objective {:.12e}, {} iterations, {} restarts, {secs:.3}s", r.objective, r.iterations, r.restarts);
    let out = json!({
        "version": 1,
        "status": status,
        "objective": if r.status == Status::Optimal { json!(r.objective) } else { Value::Null },
        "x": r.x, "dual": r.dual, "reduced_cost": r.reduced_cost,
        "diagnostics": { "iterations": r.iterations, "restarts": r.restarts, "seconds": secs,
            "primal_residual": r.primal_residual, "dual_residual": r.dual_residual, "relative": r.relative },
    });
    let text = serde_json::to_string_pretty(&out).unwrap();
    match output {
        Some(p) => std::fs::write(&p, text).unwrap_or_else(|e| fail(&format!("cannot write {p}: {e}"))),
        None => println!("{text}"),
    }
    std::process::exit(if r.status == Status::Optimal { 0 } else { 3 });
}
