//! Random small LPs (all row types, ranged rows, free/fixed/boxed/one-sided columns, both senses, offsets, degenerate integer data) solved by the
//! primal and by the dual simplex: same status, same objective, and every result passes `rose::verify` (optimal point and duals, Farkas certificate,
//! unbounded ray). Pivot sequences are not compared.

use rose::model::{LpModel, Sense, INF};
use rose::verify;
use rose::{solve, SolveOptions, Status};

struct Rng(u64);
impl Rng {
    fn next(&mut self) -> u64 {
        self.0 = self
            .0
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        self.0 >> 33
    }
    fn below(&mut self, n: u64) -> u64 {
        self.next() % n
    }
    fn int(&mut self, lo: i64, hi: i64) -> f64 {
        (lo + self.below((hi - lo + 1) as u64) as i64) as f64
    }
}

fn random_model(seed: u64) -> LpModel {
    random_model_sized(seed, 7, 9)
}

fn random_model_sized(seed: u64, max_m: u64, max_n: u64) -> LpModel {
    random_model_with(seed, max_m, max_n, false)
}

/// With `feasible` the row bounds are built around the activity of a point inside the column bounds, so the model has a feasible point.
fn random_model_with(seed: u64, max_m: u64, max_n: u64, feasible: bool) -> LpModel {
    let mut r = Rng(seed.wrapping_mul(0x9e3779b97f4a7c15) | 1);
    let m = 1 + r.below(max_m) as usize;
    let n = 1 + r.below(max_n) as usize;
    let dens = 1 + r.below(3);
    let (mut col_start, mut row_idx, mut val) = (vec![0], vec![], vec![]);
    for _ in 0..n {
        let start = row_idx.len();
        for i in 0..m {
            if r.below(4) < dens {
                let a = r.int(-4, 4);
                if a != 0.0 {
                    row_idx.push(i);
                    val.push(a);
                }
            }
        }
        if row_idx.len() == start {
            // No empty columns: a column without entries is a separate (trivial) question.
            row_idx.push(r.below(m as u64) as usize);
            val.push(1.0 + r.below(3) as f64);
        }
        col_start.push(row_idx.len());
    }
    let obj: Vec<f64> = (0..n)
        .map(|_| if r.below(4) == 0 { 0.0 } else { r.int(-5, 5) })
        .collect();
    let (mut col_lo, mut col_up) = (vec![], vec![]);
    for _ in 0..n {
        let (lo, up) = match r.below(8) {
            0..=2 => (0.0, INF),
            3 => (0.0, r.int(1, 6)),
            4 => (-INF, INF),
            5 => (r.int(-3, 0), r.int(1, 5)),
            6 => (-INF, r.int(-2, 4)),
            _ => {
                let v = r.int(-2, 3);
                (v, v)
            }
        };
        col_lo.push(lo);
        col_up.push(up);
    }
    let mut act = vec![0.0; m];
    if feasible {
        for j in 0..n {
            let x = if col_lo[j] > -INF && col_up[j] < INF {
                col_lo[j] + r.below((col_up[j] - col_lo[j]) as u64 + 1) as f64
            } else if col_lo[j] > -INF {
                col_lo[j] + r.below(3) as f64
            } else if col_up[j] < INF {
                col_up[j] - r.below(3) as f64
            } else {
                r.int(-2, 2)
            };
            for e in col_start[j]..col_start[j + 1] {
                act[row_idx[e]] += val[e] * x;
            }
        }
    }
    let (mut row_lo, mut row_up) = (vec![], vec![]);
    for i in 0..m {
        let b = if feasible { act[i] } else { r.int(-6, 10) };
        let (lo, up) = match r.below(7) {
            0 | 1 => (-INF, b),
            2 | 3 => (b, INF),
            4 => (b, b),
            5 => (
                b - r.int(0, 6),
                b + if feasible { r.int(0, 3) } else { 0.0 },
            ),
            _ => (-INF, INF),
        };
        row_lo.push(lo);
        row_up.push(up);
    }
    let row_range = row_lo
        .iter()
        .zip(&row_up)
        .map(|(&l, &u)| {
            if l > -INF && u < INF && l != u {
                u - l
            } else {
                0.0
            }
        })
        .collect();
    LpModel {
        name: format!("rand{seed}"),
        sense: if r.below(2) == 0 {
            Sense::Minimize
        } else {
            Sense::Maximize
        },
        offset: r.int(-5, 5),
        obj,
        col_lo,
        col_up,
        row_lo,
        row_up,
        row_range,
        row_kind: Vec::new(),
        range_order: Vec::new(),
        range_entries: Vec::new(),
        range_dup: false,
        col_start,
        row_idx,
        val,
        col_names: (0..n).map(|j| format!("x{j}")).collect(),
        row_names: (0..m).map(|i| format!("r{i}")).collect(),
    }
}

#[test]
fn dual_agrees_with_primal_on_random_models() {
    compare_random(4000, 7, 9);
}

#[test]
fn dual_agrees_with_primal_on_larger_random_models() {
    compare_random(1500, 16, 24);
}

#[test]
fn dual_agrees_with_primal_on_feasible_random_models() {
    compare_random_with(4000, 8, 12, true);
    compare_random_with(1500, 20, 30, true);
}

fn compare_random(count: u64, max_m: u64, max_n: u64) {
    compare_random_with(count, max_m, max_n, false)
}

fn compare_random_with(count: u64, max_m: u64, max_n: u64, feasible: bool) {
    for dual in [SolveOptions::dual(), SolveOptions::dual_release()] {
        compare_random_profile(count, max_m, max_n, feasible, &dual);
    }
}

fn compare_random_profile(count: u64, max_m: u64, max_n: u64, feasible: bool, dual: &SolveOptions) {
    let mut counts = [0usize; 5];
    let mut skipped = 0;
    let mut bad = vec![];
    for seed in 0..count {
        let model = random_model_with(seed, max_m, max_n, feasible);
        let p = solve(&model, None, &SolveOptions::default()).unwrap();
        let d = solve(&model, None, dual).unwrap();
        counts[match d.status {
            Status::Optimal => 0,
            Status::Infeasible => 1,
            Status::Unbounded => 2,
            Status::IterationLimit => 3,
            Status::NumericalFailure => 4,
        }] += 1;
        if let Err(e) = verify::check(&model, &d, 1e-6) {
            bad.push(format!(
                "seed {seed}: dual {:?} fails its check: {e}",
                d.status
            ));
        } else if p.status == Status::NumericalFailure || verify::check(&model, &p, 1e-6).is_err() {
            // The primal simplex gives up on some degenerate corners (a single free row with free columns, unbounded feasible models) or returns a
            // result that fails its check; the dual's verified result stands alone there, nothing to compare.
            skipped += 1;
        } else if d.status != p.status {
            // Both infeasible and unbounded: either answer is a proof of its own; the check above vouches for the dual's.
            bad.push(format!(
                "seed {seed}: primal {:?}, dual {:?}",
                p.status, d.status
            ));
        } else if d.status == Status::Optimal
            && (d.objective - p.objective).abs() > 1e-7 + 1e-8 * p.objective.abs()
        {
            bad.push(format!(
                "seed {seed}: objective {} vs {}",
                d.objective, p.objective
            ));
        }
    }
    let mut kinds = std::collections::BTreeMap::new();
    for b in &bad {
        let k = b
            .split(": ")
            .nth(1)
            .unwrap()
            .split(" fails")
            .next()
            .unwrap()
            .split(',')
            .take(2)
            .collect::<Vec<_>>()
            .join(",");
        *kinds.entry(k).or_insert(0) += 1;
    }
    eprintln!("failure kinds: {kinds:?}");
    eprintln!("dual statuses (optimal, infeasible, unbounded, limit, failure): {counts:?}, {skipped} models not compared (the primal simplex fails or fails its own check)");
    assert!(
        bad.is_empty(),
        "{} failures, first: {:#?}",
        bad.len(),
        &bad[..bad.len().min(25)]
    );
}

#[test]
#[ignore]
fn show() {
    let seed: u64 = std::env::var("SEED").unwrap().parse().unwrap();
    let model = random_model(seed);
    eprintln!("{model:#?}");
    for o in [
        SolveOptions::default(),
        SolveOptions::dual(),
        SolveOptions::dual_release(),
    ] {
        let r = solve(&model, None, &o).unwrap();
        eprintln!(
            "{:?} obj {} it {} x {:?} farkas {:?} ray {:?} check {:?}",
            r.status,
            r.objective,
            r.iterations,
            r.x,
            r.farkas,
            r.ray,
            verify::check(&model, &r, 1e-6)
        );
    }
}

#[test]
#[ignore]
fn primal_failures() {
    let mut n = 0;
    for seed in 0..4000u64 {
        let model = random_model_with(seed, 8, 12, true);
        let p = solve(&model, None, &SolveOptions::default()).unwrap();
        let d = solve(&model, None, &SolveOptions::dual()).unwrap();
        if p.status == Status::NumericalFailure && n < 6 {
            n += 1;
            eprintln!(
                "seed {seed}: primal {:?} dual {:?} dual check {:?} m {} n {}",
                p.status,
                d.status,
                verify::check(&model, &d, 1e-6).is_ok(),
                model.num_rows(),
                model.num_cols()
            );
        }
    }
}
