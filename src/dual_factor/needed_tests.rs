//! The needed U kernel of F3a (`ftran_needed`: dense tail slots walked in place, counts booked by `TailCounts`, one-pass tail count)
//! against the frozen kernel of the stack before it (`ftran_needed_old`): output order and bits, WORK, the F2 state the joint stages
//! read, and the tau vector, over random factors and update sequences (repeated calls on the same factor catch stale scratch).
use super::handoff_tests::{all_bits, assert_same, bits};
use super::*;

fn run<R: Real>(seeds: &[u64], vals_pool: &[f64]) -> usize {
    let mut checks = 0usize;
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
                                vals.push(vals_pool[next() % vals_pool.len()]);
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
                let n = [1, 2, 3, 5, m / 8 + 1, m / 4 + 1, m / 2 + 1, m][next() % 8].min(m);
                let mut rows: Vec<usize> = Vec::new();
                for _ in 0..n {
                    let r = next() % m;
                    if !rows.contains(&r) {
                        rows.push(r);
                    }
                }
                let vals: Vec<R> = rows.iter().map(|_| R::c([0.25, -1.5, 1e-15, 3.0, -0.0, 0.0, TOL, -TOL, 1e-14, 2e-14][next() % 10])).collect();
                // The first vector of the pair, then the second one through the new and the frozen kernel.
                let x = a.ftran(&rows, &vals);
                let y = b.ftran(&rows, &vals);
                assert_eq!(bits(&x), bits(&y), "m {m} step {step}: ftran output");
                let mut needed: Vec<usize> = x.iter().map(|e| e.0).collect();
                if next() % 3 == 0 {
                    // Arbitrary positions (singletons, nucleus rows and tail slots), also repeated ones.
                    needed = (0..1 + next() % m.min(30)).map(|_| next() % m).collect();
                }
                let (rows2, vals2): (Vec<usize>, Vec<R>) = if next() % 4 == 0 {
                    (rows.clone(), vals.iter().map(|v| *v).collect())
                } else {
                    (rows.iter().rev().copied().collect(), vals.iter().rev().copied().collect())
                };
                let dense_in = next() % 3 == 0;
                let p = a.ftran_needed(&rows2, &vals2, &needed, dense_in);
                let q = b.ftran_needed_old(&rows2, &vals2, &needed, dense_in);
                assert_eq!(bits(&p), bits(&q), "m {m} step {step}: needed output");
                assert_same(&mut a, &mut b, "ftran_needed");
                assert_eq!(all_bits(&a.f2.get_mut().x_tau), all_bits(&b.f2.get_mut().x_tau), "x_tau");
                assert_eq!(a.joint_u(&needed).map(f64::to_bits), b.joint_u(&needed).map(f64::to_bits), "joint_u");
                checks += 1;
                tail_pair(&a, &mut next);
                if step % 3 == 2 {
                    let r = next() % m;
                    let mut rows = vec![r];
                    let mut vals = vec![R::c(5.0 + step as f64)];
                    for _ in 0..(1 + next() % 6) {
                        let q = next() % m;
                        if !rows.contains(&q) {
                            rows.push(q);
                            vals.push(R::c([0.25, -0.5, 1.5][next() % 3]));
                        }
                    }
                    let x = a.ftran(&rows, &vals);
                    let y = b.ftran(&rows, &vals);
                    assert_eq!(bits(&x), bits(&y));
                    a.update(r);
                    b.update(r);
                    assert_eq!(a.upd_work, b.upd_work);
                }
            }
        }
    }
    checks
}

/// The stream loop against the `at`/`diag`/`ab` loop it replaces (`needed_tail::<true, false, TD>` over the descending tail slots):
/// the vector bits, the raw flags and every booked count, on vectors with zeros of both signs, tiny and NaN values.
fn tail_pair<R: Real>(f: &DualFactor<R>, next: &mut impl FnMut() -> usize) {
    let m = f.m;
    let pool = [0.0, -0.0, 0.5, -1.5, 1e-15, -1e-14, 3.0, TOL, -TOL, 2e-14, 7.0, f64::NAN];
    for round in 0..4 {
        let x0: Vec<R> = (0..m)
            .map(|_| match next() % (2 + round * 2) {
                0 | 1 => R::ZERO,
                _ => R::c(pool[next() % pool.len()]),
            })
            .collect();
        let x_col: Vec<R> = (0..m).map(|_| R::c(pool[next() % 4])).collect();
        macro_rules! both {
            ($td:literal) => {{
                let (mut xa, mut xb) = (x0.clone(), x0.clone());
                let (mut ra, mut rb) = (vec![false; m], vec![false; m]);
                let order = f.at[f.i21..m].iter().rev().copied();
                let a = f.needed_tail::<true, false, $td>(&mut xa, order, &x_col, &mut ra);
                let b = f.needed_tail_stream::<$td>(&mut xb, &x_col, &mut rb);
                assert_eq!(all_bits(&xa), all_bits(&xb), "stream tail x");
                assert_eq!(ra, rb, "stream tail raw flags");
                let key = |t: &TailCounts| (t.nz, t.abs_a, t.tp_a2, t.tp_only_a, t.surv_nz, t.s, t.n, t.n_tiny, t.nan_a);
                assert_eq!(key(&a), key(&b), "stream tail counts");
            }};
        }
        both!(false);
        both!(true);
    }
}

const SMOOTH: &[f64] = &[0.5, -0.75, 1.25, -0.0, 1e-15];
const WILD: &[f64] = &[0.5, -0.75, 1.25, 0.0, 1e-15, -2.0, 1e-13, 0.125, -1e-14];

#[test]
fn needed_kernel_matches_the_frozen_one_f64() {
    let seeds = [0x2545f4914f6cdd1d, 0x9e3779b97f4a7c15, 0x43e941132abcu64, 0xdeadbeefcafe1234, 0x1234567887654321, 0xfeedfacecafebeef];
    assert!(run::<f64>(&seeds, SMOOTH) > 200);
    assert!(run::<f64>(&seeds, WILD) > 200);
}

#[test]
fn needed_kernel_matches_the_frozen_one_f80() {
    let seeds = [0x2545f4914f6cdd1d, 0x43e941132abcu64, 0x1234567887654321];
    assert!(run::<crate::f80::F80>(&seeds, SMOOTH) > 100);
    assert!(run::<crate::f80::F80>(&seeds, WILD) > 100);
}

#[test]
fn tail_nonzero_matches_the_kind_filter() {
    let m = 60;
    let cols: Vec<_> = (0..m)
        .map(|i| if i % 4 == 0 { (false, vec![i], vec![2.0]) } else { (false, vec![i, (i * 7 + 3) % m], vec![4.0, 0.5]) })
        .collect();
    let st = crate::lu_order::factor_structure_r::<f64>(m, &cols);
    let f = DualFactor::<f64>::fresh(m, m, &st, &cols);
    for k in 0..m {
        let x: Vec<f64> = (0..m).map(|p| if (p * 31 + k) % 5 < 2 { 0.0 } else { -0.0 + p as f64 }).collect();
        let want = (0..m).filter(|&p| f.kind[p] == -2 && x[p] != 0.0).count();
        assert_eq!(f.tail_nonzero(&x), want);
    }
}
