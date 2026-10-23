//! Metamorphic tests: transformed models must give the same optimal objective (up to the
//! known effect of the transformation) and pass the independent check in `verify`. Pivot sequences may differ.
//! Models: tests/fixtures/*.mps and reverse-engineer/corpus/fixtures/fx_*.mps (skipped if absent).

use rose::model::{Sense, INF};
use rose::verify;
use rose::{mps, solve, LpModel, SolveOptions, Status};

fn load(path: String) -> Option<LpModel> {
    std::fs::read_to_string(path).ok().map(|t| mps::parse(&t).unwrap())
}

fn small() -> Vec<(String, LpModel)> {
    let d = env!("CARGO_MANIFEST_DIR");
    let mut v = vec![];
    for n in ["afiro", "sc50a", "blend", "adlittle", "bore3d", "boeing2"] {
        if let Some(m) = load(format!("{d}/tests/fixtures/{n}.mps")) {
            v.push((n.to_string(), m));
        }
    }
    for n in ["fx_tied", "fx_neardep", "fx_badscale", "fx_singular", "fx_warm", "fx_infeasible", "fx_unbounded"] {
        if let Some(m) = load(format!("{d}/../reverse-engineer/corpus/fixtures/{n}.mps")) {
            v.push((n.to_string(), m));
        }
    }
    v
}

/// The methods every case runs under: the primal simplex and the dual simplex (the controlled profiles; pivot sequences differ).
fn opts() -> Vec<SolveOptions> {
    vec![SolveOptions::default(), SolveOptions::dual(), SolveOptions::dual_release()]
}

fn fin(b: f64) -> bool {
    b.abs() < INF
}

fn sc(b: f64, s: f64) -> f64 {
    if fin(b) { b * s } else { b }
}

/// Deterministic pseudo-random permutation of 0..n.
fn perm(n: usize, seed: u64) -> Vec<usize> {
    let mut p: Vec<usize> = (0..n).collect();
    let mut s = seed.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
    for i in (1..n).rev() {
        s = s.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
        p.swap(i, (s >> 33) as usize % (i + 1));
    }
    p
}

/// New row `k` is old row `p[k]`; new column `k` is old column `q[k]`.
fn permute(m: &LpModel, p: &[usize], q: &[usize]) -> LpModel {
    let mut inv = vec![0; p.len()];
    for (k, &o) in p.iter().enumerate() {
        inv[o] = k;
    }
    let mut r = m.clone();
    r.row_lo = p.iter().map(|&i| m.row_lo[i]).collect();
    r.row_up = p.iter().map(|&i| m.row_up[i]).collect();
    r.row_names = p.iter().map(|&i| m.row_names[i].clone()).collect();
    r.obj = q.iter().map(|&j| m.obj[j]).collect();
    r.col_lo = q.iter().map(|&j| m.col_lo[j]).collect();
    r.col_up = q.iter().map(|&j| m.col_up[j]).collect();
    r.col_names = q.iter().map(|&j| m.col_names[j].clone()).collect();
    r.col_start = vec![0];
    r.row_idx.clear();
    r.val.clear();
    for &j in q {
        let mut e: Vec<(usize, f64)> = m.col(j).map(|(i, a)| (inv[i], a)).collect();
        e.sort_by_key(|x| x.0);
        for (i, a) in e {
            r.row_idx.push(i);
            r.val.push(a);
        }
        r.col_start.push(r.row_idx.len());
    }
    r
}

/// Multiply row `i` by `s > 0` (every `i % every == 0`).
fn scale_rows(m: &LpModel, s: f64, every: usize) -> LpModel {
    let mut r = m.clone();
    for i in (0..m.num_rows()).step_by(every) {
        r.row_lo[i] = sc(m.row_lo[i], s);
        r.row_up[i] = sc(m.row_up[i], s);
    }
    for k in 0..r.val.len() {
        if r.row_idx[k] % every == 0 {
            r.val[k] *= s;
        }
    }
    r
}

/// Substitute `x_j = s * x'_j` for every `j % every == 0`, `s > 0`.
fn scale_cols(m: &LpModel, s: f64, every: usize) -> LpModel {
    let mut r = m.clone();
    for j in (0..m.num_cols()).step_by(every) {
        r.obj[j] *= s;
        r.col_lo[j] = sc(m.col_lo[j], 1.0 / s);
        r.col_up[j] = sc(m.col_up[j], 1.0 / s);
        for k in r.col_start[j]..r.col_start[j + 1] {
            r.val[k] *= s;
        }
    }
    r
}

fn negate(m: &LpModel) -> LpModel {
    let mut r = m.clone();
    r.sense = if m.sense == Sense::Minimize { Sense::Maximize } else { Sense::Minimize };
    r.obj.iter_mut().for_each(|c| *c = -*c);
    r.offset = -m.offset;
    r
}

fn run(m: &LpModel, name: &str, what: &str, o: &SolveOptions) -> rose::SolveResult {
    let r = solve(m, None, o).unwrap();
    if r.status == Status::Optimal {
        // scaled models have larger magnitudes, so use a mildly loose tolerance on the independent check
        if let Err(e) = verify::check(m, &r, 1e-6) {
            panic!("{name} [{what}]: check failed: {e}");
        }
    }
    r
}

/// Solve the original and `t(model)`; objective of the transformed model must equal `f(original objective)`.
fn case(what: &str, t: impl Fn(&LpModel) -> LpModel, f: impl Fn(f64) -> f64) {
    for (name, m) in small() {
        for o in opts() {
            let base = run(&m, &name, "base", &o);
            let r = run(&t(&m), &name, what, &o);
            assert_eq!(r.status, base.status, "{name} [{what}]: status");
            if base.status == Status::Optimal {
                let want = f(base.objective);
                let err = (r.objective - want).abs() / (1.0 + want.abs());
                assert!(err <= 1e-8, "{name} [{what}]: objective {} vs {} (rel {err:.2e})", r.objective, want);
            }
        }
    }
}

#[test]
fn row_and_column_permutation() {
    case("perm", |m| permute(m, &perm(m.num_rows(), 1), &perm(m.num_cols(), 2)), |z| z);
    case("rows reversed", |m| permute(m, &(0..m.num_rows()).rev().collect::<Vec<_>>(), &(0..m.num_cols()).collect::<Vec<_>>()), |z| z);
    case("cols reversed", |m| permute(m, &(0..m.num_rows()).collect::<Vec<_>>(), &(0..m.num_cols()).rev().collect::<Vec<_>>()), |z| z);
}

#[test]
fn row_scaling() {
    case("rows x16", |m| scale_rows(m, 16.0, 2), |z| z);
    case("rows x1/8", |m| scale_rows(m, 0.125, 3), |z| z);
    case("rows x7", |m| scale_rows(m, 7.0, 2), |z| z);
    case("rows x0.3", |m| scale_rows(m, 0.3, 3), |z| z);
}

#[test]
fn column_scaling() {
    case("cols x8", |m| scale_cols(m, 8.0, 2), |z| z);
    case("cols x1/4", |m| scale_cols(m, 0.25, 3), |z| z);
    case("cols x3", |m| scale_cols(m, 3.0, 2), |z| z);
    case("cols x0.7", |m| scale_cols(m, 0.7, 3), |z| z);
}

#[test]
fn negated_objective_with_sense_switch() {
    case("negate", negate, |z| -z);
}

#[test]
fn objective_offset() {
    case("offset +123.5", |m| LpModel { offset: m.offset + 123.5, ..m.clone() }, |z| z + 123.5);
    case("offset -7 then negate", |m| negate(&LpModel { offset: m.offset - 7.0, ..m.clone() }), |z| -(z - 7.0));
}
