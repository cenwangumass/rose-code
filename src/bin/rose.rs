//! `rose solve model.mps --output result.json [--basis-in b.json] [--basis-out b.json] [--method primal|dual]`
//! Progress goes to stderr; stdout stays empty so machine-readable output is only ever in the files.

use rose::verify;
use rose::{mps, solve, Basis, SolveOptions, Status};
use serde_json::{json, Value};

fn usage() -> ! {
    eprintln!("usage: rose solve <model.mps> [--output result.json] [--basis-in basis.json] [--basis-out basis.json] [--max-iterations N] [--method primal|dual]");
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
    let mut opts = SolveOptions { record_trace: false, ..SolveOptions::default() };
    let mut it = args[2..].iter();
    while let Some(a) = it.next() {
        let val = it.next().unwrap_or_else(|| usage());
        match a.as_str() {
            "--output" => output = Some(val.clone()),
            "--basis-in" => basis_in = Some(val.clone()),
            "--basis-out" => basis_out = Some(val.clone()),
            "--max-iterations" => opts.max_iterations = val.parse().unwrap_or_else(|_| usage()),
            "--method" => {
                opts.method = match val.as_str() {
                    "primal" => rose::method::SolveMethod::Primal,
                    "dual" => rose::method::SolveMethod::Dual,
                    _ => usage(),
                }
            }
            _ => usage(),
        }
    }
    let text = std::fs::read_to_string(path).unwrap_or_else(|e| fail(&format!("cannot read {path}: {e}")));
    let model = mps::parse(&text).unwrap_or_else(|e| fail(&format!("invalid input: {e}")));
    let start = basis_in.map(|p| {
        let v: Value = serde_json::from_str(&std::fs::read_to_string(&p).unwrap_or_else(|e| fail(&format!("cannot read {p}: {e}")))).unwrap_or_else(|e| fail(&format!("invalid basis file: {e}")));
        basis_from_json(&v).unwrap_or_else(|| fail("invalid basis file"))
    });
    eprintln!("{}: {} rows, {} columns", model.name, model.num_rows(), model.num_cols());
    let t = std::time::Instant::now();
    let r = solve(&model, start.as_ref(), &opts).unwrap_or_else(|e| fail(&format!("invalid input: {e}")));
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
        "x": r.x, "dual": r.dual, "reduced_cost": r.reduced_cost,
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
