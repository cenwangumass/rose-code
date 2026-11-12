//! The full U kernels of F3b (`ftran_with`, `ftran_with_dense`: one lookup per column, direct block marks, fused final pass and output
//! scratch, no normalizing pass over the dense L-stage vector) against the frozen ones of the stack before it (`ftran_with_old`): output
//! order and bits, WORK, the spike, the F2 state (column snapshot, dirty blocks, survivors) and the second vector of the pair, over random
//! factors and update sequences. Repeated calls on the same factors catch stale scratch.
use super::handoff_tests::{all_bits, assert_same_with, bits};

/// The dirty blocks of the sparse U stage are a subset of the older ones (P6-4), so they are compared by coverage.
fn assert_same<R: Real>(a: &mut DualFactor<R>, b: &mut DualFactor<R>, what: &str) {
    assert_same_with(a, b, what, false);
}
use super::*;

fn run<R: Real>(seeds: &[u64], vals_pool: &[f64]) -> (usize, usize) {
    let (mut checks, mut dense_out) = (0usize, 0usize);
    for &seed0 in seeds {
        let mut seed = seed0;
        let mut next = || {
            seed ^= seed << 13;
            seed ^= seed >> 7;
            seed ^= seed << 17;
            seed as usize
        };
        for m in [9, 40, 130, 300] {
            let cols: Vec<_> = (0..m)
                .map(|i| match next() % 6 {
                    0 => (true, vec![i], vec![1.0]),
                    1 => (false, vec![i], vec![2.0]),
                    _ => {
                        let mut rows = vec![i];
                        let mut vals = vec![4.0];
                        for _ in 0..(1 + next() % 5) {
                            let r = next() % m;
                            if !rows.contains(&r) {
                                rows.push(r);
                                vals.push(0.25 * vals_pool[next() % vals_pool.len()]);
                            }
                        }
                        (false, rows, vals)
                    }
                })
                .collect();
            let st = crate::lu_order::factor_structure_r::<R>(m, &cols);
            let mut a = DualFactor::<R>::fresh(m, m, &st, &cols);
            let mut b = DualFactor::<R>::fresh(m, m, &st, &cols);
            for step in 0..m.min(40) {
                // A (nearly) singular draw breeds NaNs, for which the work estimates of both kernels are not compared with their oracle.
                if a.diag.iter().any(|d| !d.to_f64().is_normal() || !(1e-9..1e9).contains(&d.to_f64().abs())) {
                    a = DualFactor::<R>::fresh(m, m, &st, &cols);
                    b = DualFactor::<R>::fresh(m, m, &st, &cols);
                }
                let n = [1, 2, 3, 5, m / 8 + 1, m / 4 + 1, m / 2 + 1, m][next() % 8].min(m);
                let mut rows: Vec<usize> = Vec::new();
                for _ in 0..n {
                    let r = next() % m;
                    if !rows.contains(&r) {
                        rows.push(r);
                    }
                }
                let vals: Vec<R> = rows.iter().map(|_| R::c([0.25, -1.5, 1e-15, 3.0, -0.0, 0.0, TOL, -TOL, 1e-14, 2e-14][next() % 10])).collect();
                let force = next() % 3 == 0;
                let y = b.ftran_with_old(&rows, &vals, force);
                let x = a.ftran_with(&rows, &vals, force);
                assert_eq!(bits(&x), bits(&y), "m {m} step {step}: output");
                assert_same(&mut a, &mut b, "ftran_with");
                dense_out += usize::from(a.last_dense);
                checks += 1;
                // The second vector of the pair reads the column's F2 state.
                if next() % 2 == 0 {
                    let needed: Vec<usize> = x.iter().map(|e| e.0).collect();
                    let dense_in = next() % 3 == 0;
                    let p = a.ftran_needed(&rows, &vals, &needed, dense_in);
                    let q = b.ftran_needed_old(&rows, &vals, &needed, dense_in);
                    assert_eq!(bits(&p), bits(&q), "m {m} step {step}: needed output");
                    assert_same(&mut a, &mut b, "ftran_needed");
                    assert_eq!(all_bits(&a.f2.get_mut().x_tau), all_bits(&b.f2.get_mut().x_tau), "x_tau");
                    assert_eq!(a.joint_u(&needed).map(f64::to_bits), b.joint_u(&needed).map(f64::to_bits), "joint_u");
                }
                // The kept and the basis variants wrap the same kernel.
                if next() % 5 == 0 {
                    let kx = a.ftran_keep(&rows, &vals, force);
                    let ky = b.ftran_keep(&rows, &vals, force);
                    assert_eq!(bits(&kx), bits(&ky));
                }
                if step % 3 == 2 {
                    let r = next() % m;
                    let mut rows = vec![r];
                    let mut vals = vec![R::c(20.0 + step as f64)];
                    for _ in 0..(1 + next() % 6) {
                        let q = next() % m;
                        if !rows.contains(&q) {
                            rows.push(q);
                            vals.push(R::c([0.25, -0.5, 1.5][next() % 3]));
                        }
                    }
                    let x = a.ftran(&rows, &vals);
                    let y = b.ftran_with_old(&rows, &vals, false);
                    assert_eq!(bits(&x), bits(&y));
                    assert_same(&mut a, &mut b, "pivot column");
                    a.update(r);
                    b.update(r);
                    assert_eq!(a.upd_work, b.upd_work);
                }
            }
        }
    }
    (checks, dense_out)
}

const SMOOTH: &[f64] = &[0.5, -0.75, 1.25, -0.0, 1e-15];
const WILD: &[f64] = &[0.5, -0.75, 1.25, 0.0, 1e-15, -2.0, 1e-13, 0.125, -1e-14];

#[test]
fn full_kernels_match_the_frozen_ones_f64() {
    let seeds = [0x2545f4914f6cdd1d, 0x9e3779b97f4a7c15, 0x43e941132abcu64, 0xdeadbeefcafe1234, 0x1234567887654321, 0xfeedfacecafebeef];
    for pool in [SMOOTH, WILD] {
        let handoffs = DENSE_HANDOFFS.with(|k| k.get());
        let (checks, dense) = run::<f64>(&seeds, pool);
        assert!(DENSE_HANDOFFS.with(|k| k.get()) > handoffs + 20, "the dense L-stage handoff must run");
        assert!(checks > 200);
        assert!(dense > 20 && dense < checks, "both the dense and the sparse routes must run: {dense} of {checks}");
    }
}

#[test]
fn full_kernels_match_the_frozen_ones_f80() {
    let seeds = [0x2545f4914f6cdd1d, 0x43e941132abcu64, 0x1234567887654321];
    for pool in [SMOOTH, WILD] {
        assert!(run::<crate::f80::F80>(&seeds, pool).0 > 100);
    }
}
