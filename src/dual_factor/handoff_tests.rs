//! The dense L-stage vector feeds the U stage directly; the old list round trip (`dense_spike`, re-scatter into a zero-filled
//! vector) is kept behind the `FORCE_LIST_HANDOFF` test hook and must give identical bits, orders and WORK.
use super::*;

pub(super) fn bits<R: Real>(v: &[(usize, R)]) -> Vec<(usize, u64)> {
    v.iter().map(|&(i, x)| (i, x.to_f64().to_bits())).collect()
}

pub(super) fn all_bits<R: Real>(v: &[R]) -> Vec<u64> {
    v.iter().map(|x| x.to_f64().to_bits()).collect()
}

fn list_handoff<T>(on: bool, f: impl FnOnce() -> T) -> T {
    FORCE_LIST_HANDOFF.with(|k| k.set(on));
    let out = f();
    FORCE_LIST_HANDOFF.with(|k| k.set(false));
    out
}

/// Everything either path leaves behind that a later consumer reads.
pub(super) fn assert_same<R: Real>(a: &mut DualFactor<R>, b: &mut DualFactor<R>, what: &str) {
    assert_same_with(a, b, what, true);
}

/// `exact_dirty`: the lists of dirty blocks are equal; otherwise `a`'s are a subset of `b`'s that still covers every nonzero of `a`'s column.
pub(super) fn assert_same_with<R: Real>(a: &mut DualFactor<R>, b: &mut DualFactor<R>, what: &str, exact_dirty: bool) {
    assert_eq!(a.work.to_bits(), b.work.to_bits(), "{what}: work");
    assert_eq!((a.work_parts.0.to_bits(), a.work_parts.1.to_bits()), (b.work_parts.0.to_bits(), b.work_parts.1.to_bits()), "{what}: parts");
    assert_eq!((a.last_dense, a.last_ordered, a.last_alpha), (b.last_dense, b.last_ordered, b.last_alpha), "{what}: flags");
    assert_eq!((a.tau_l.get().to_bits(), a.tau_u.get().to_bits(), a.tau_fin.get()), (b.tau_l.get().to_bits(), b.tau_u.get().to_bits(), b.tau_fin.get()), "{what}: tau");
    assert_eq!(bits(&a.spike), bits(&b.spike), "{what}: spike");
    assert_eq!(a.spike_used, b.spike_used, "{what}: spike_used");
    let (fa, fb) = (a.f2.get_mut(), b.f2.get_mut());
    assert_eq!((fa.col_n, &fa.col_in, &fa.col_surv, fa.col_sa, fa.col_sb), (fb.col_n, &fb.col_in, &fb.col_surv, fb.col_sa, fb.col_sb), "{what}: f2 column");
    assert_eq!(all_bits(&fa.x_col), all_bits(&fb.x_col), "{what}: x_col");
    if exact_dirty {
        assert_eq!(fa.dirty_blocks, fb.dirty_blocks, "{what}: dirty");
    } else {
        let (da, db) = (fa.dirty_blocks.as_deref(), fb.dirty_blocks.as_deref());
        if let (Some(da), Some(db)) = (da, db) {
            assert!(da.iter().all(|k| db.contains(k)), "{what}: dirty blocks of the new route are not a subset");
            assert!(da.windows(2).all(|w| w[0] < w[1]), "{what}: dirty blocks not ascending");
        }
        if let Some(da) = da {
            for (i, &v) in fa.x_col.iter().enumerate() {
                assert!(v == R::ZERO || da.contains(&(i / 64)), "{what}: nonzero outside the dirty blocks at {i}");
            }
        }
    }
    assert_eq!((fa.tau_n, fa.tp_a2, fa.tp_only_a), (fb.tau_n, fb.tp_a2, fb.tp_only_a), "{what}: f2 tau");
    assert_eq!(all_bits(&fa.x_tau), all_bits(&fb.x_tau), "{what}: x_tau");
    assert_eq!(fa.tau_raw, fb.tau_raw, "{what}: tau_raw");
}

fn run<R: Real>(seeds: &[u64]) -> (usize, usize) {
    let (mut dense_calls, mut checks) = (0usize, 0usize);
    for &seed0 in seeds {
        let mut seed = seed0;
        let mut next = || {
            seed ^= seed << 13;
            seed ^= seed >> 7;
            seed ^= seed << 17;
            seed as usize
        };
        for m in [9, 40, 130, 300] {
            let cols: Vec<_> = (0..m).map(|i| match next() % 5 {
                0 => (true, vec![i], vec![1.0]),
                1 => (false, vec![i], vec![2.0]),
                _ => {
                    let mut rows = vec![i];
                    let mut vals = vec![4.0];
                    for _ in 0..(1 + next() % 4) {
                        let r = next() % m;
                        if !rows.contains(&r) {
                            rows.push(r);
                            vals.push([0.5, -0.75, 1.25, -0.0, 1e-15][next() % 5]);
                        }
                    }
                    (false, rows, vals)
                }
            }).collect();
            let st = crate::lu_order::factor_structure_r::<R>(m, &cols);
            let mut a = DualFactor::<R>::fresh(m, m, &st, &cols);
            let mut b = DualFactor::<R>::fresh(m, m, &st, &cols);
            for step in 0..m.min(36) {
                // Sparse, medium and dense right-hand sides: a dense one makes the L stage dense, the rest depends on the counts.
                let n = [1, 2, 5, m / 8 + 1, m / 2 + 1, m][next() % 6].min(m);
                let mut rows: Vec<usize> = Vec::new();
                for _ in 0..n {
                    let r = next() % m;
                    if !rows.contains(&r) {
                        rows.push(r);
                    }
                }
                let vals: Vec<R> = rows.iter().map(|_| R::c([0.25, -1.5, 1e-15, 3.0, -0.0, 0.0, TOL, -TOL, 1e-14][next() % 9])).collect();
                let force = next() % 3 == 0;
                let before = DENSE_HANDOFFS.with(|k| k.get());
                let (x, y) = match next() % 3 {
                    0 => (a.ftran(&rows, &vals), list_handoff(true, || b.ftran(&rows, &vals))),
                    1 => (a.ftran_basic(&rows, &vals), list_handoff(true, || b.ftran_basic(&rows, &vals))),
                    _ => (a.ftran_keep(&rows, &vals, force), list_handoff(true, || b.ftran_keep(&rows, &vals, force))),
                };
                dense_calls += DENSE_HANDOFFS.with(|k| k.get()) - before;
                assert_eq!(bits(&x), bits(&y), "m {m} step {step}: ftran output");
                assert_same(&mut a, &mut b, "ftran");
                // The second vector of the pair, then the joint U work.
                let needed: Vec<usize> = x.iter().map(|e| e.0).collect();
                let rows2: Vec<usize> = rows.iter().rev().copied().collect();
                let vals2: Vec<R> = vals.iter().rev().copied().collect();
                let dense_in = next() % 2 == 0;
                let before = DENSE_HANDOFFS.with(|k| k.get());
                let p = a.ftran_needed(&rows2, &vals2, &needed, dense_in);
                let q = list_handoff(true, || b.ftran_needed(&rows2, &vals2, &needed, dense_in));
                dense_calls += DENSE_HANDOFFS.with(|k| k.get()) - before;
                assert_eq!(bits(&p), bits(&q), "m {m} step {step}: needed output");
                assert_same(&mut a, &mut b, "ftran_needed");
                let ju = a.joint_u(&needed).map(f64::to_bits);
                assert_eq!(ju, b.joint_u(&needed).map(f64::to_bits), "joint_u");
                checks += 1;
                if step % 3 == 2 {
                    // A pivot: the spike the dense route builds feeds the update, and the etas make later L stages differ.
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
                    let y = list_handoff(true, || b.ftran(&rows, &vals));
                    assert_eq!(bits(&x), bits(&y));
                    a.update(r);
                    b.update(r);
                    assert_eq!(a.upd_work, b.upd_work);
                }
            }
        }
    }
    (dense_calls, checks)
}

#[test]
fn dense_handoff_matches_the_list_round_trip_f64() {
    let (dense, checks) = run::<f64>(&[0x2545f4914f6cdd1d, 0x9e3779b97f4a7c15, 0x43e941132abcu64, 0xdeadbeefcafe1234]);
    assert!(dense > 100 && checks > 100, "dense handoffs {dense}, checks {checks}");
}

#[test]
fn dense_handoff_matches_the_list_round_trip_f80() {
    // `LANES` is false for F80: the list stays, and must stay unchanged.
    let (dense, checks) = run::<crate::f80::F80>(&[0x2545f4914f6cdd1d, 0x43e941132abcu64]);
    assert_eq!(dense, 0);
    assert!(checks > 20);
}

#[test]
fn normalize_count_counts_nonzero_entries() {
    let nan = f64::from_bits(0x7ff8_0000_dead_beef);
    let sub = f64::from_bits(1);
    let mut x = vec![0.0, 0.0, 1e-300, -sub, nan, -1.5, f64::INFINITY, 0.0, 0.0, 1e-15];
    let want = all_bits(&x);
    let nz = DualFactor::<f64>::normalize_count(&mut x);
    assert_eq!(nz, 6);
    assert_eq!(all_bits(&x), want);
}

/// The no -0.0 contract of the dense L-stage vector is checked in the test, debug and shadow builds.
#[test]
#[should_panic(expected = "-0.0")]
fn normalize_count_rejects_a_negative_zero() {
    DualFactor::<f64>::normalize_count(&mut [1.0, -0.0]);
}

#[test]
fn dense_lists_match_the_filters() {
    let mut inp = vec![99usize; 3];
    let mut sp = vec![(7usize, 1.0f64); 2];
    for x in [
        vec![],
        vec![0.0; 5],
        vec![1.0, 0.0, 1e-15, -2.0, TOL, -TOL, f64::NAN, 5e-14, 0.0, 0.0],
        (0..70).map(|i| if i % 3 == 0 { 0.0 } else { i as f64 * 1e-14 }).collect(),
    ] {
        let nz = DualFactor::<f64>::dense_lists(&x, &mut inp, &mut sp);
        let want_i: Vec<usize> = (0..x.len()).filter(|&i| x[i] != 0.0).collect();
        let want_s: Vec<(usize, f64)> = want_i.iter().map(|&i| (i, x[i])).filter(|e| e.1.abs() > TOL).collect();
        assert_eq!(nz, want_i.len());
        assert_eq!(inp, want_i);
        assert_eq!(bits(&sp), bits(&want_s));
    }
}
