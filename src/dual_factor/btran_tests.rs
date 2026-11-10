//! Differential oracle: `ftran_order::RefStruct` retains the pre-workspace kernels.
use super::*;
use crate::ftran_order::RefStruct;

fn bits(v: Vec<(usize, f64)>) -> Vec<(usize, u64)> {
    v.into_iter().map(|(i, x)| (i, x.to_bits())).collect()
}

#[test]
fn retained_slot_counts_keep_duplicates_and_clear_between_calls() {
    let m = 321;
    let cols: Vec<_> = (0..m).map(|i| (true, vec![i], vec![1.0])).collect();
    let st = crate::lu_order::factor_structure(m, &cols);
    let mut f = DualFactor::fresh(m, m, &st, &cols);
    f.at.rotate_left(71);
    for (s, &p) in f.at.iter().enumerate() { f.slot[p] = s; }
    for pass in 0..16 {
        for n in [0, 1, 63, 64, 65, 200, 500] {
            let mut v: Vec<_> = (0..n).map(|i| (i * 47 + pass * 13) % m).collect();
            let mut expected = v.clone();
            expected.sort_unstable_by_key(|&p| std::cmp::Reverse(f.slot[p]));
            f.sort_slot_desc(&mut v);
            assert_eq!(v, expected);
            assert!(f.slot_counts.borrow().iter().all(|&n| n == 0));
        }
    }
}

#[test]
fn retained_btran_work_and_updates_match_original_repeatedly() {
    let mut seed = 0x43e941132abcu64;
    let mut next = || {
        seed ^= seed << 13;
        seed ^= seed >> 7;
        seed ^= seed << 17;
        seed as usize
    };
    for m in [17, 65, 321, 4096] {
        let cols: Vec<_> = (0..m).map(|i| {
            let mut rows = vec![i];
            let mut vals = vec![2.0];
            for _ in 0..3 {
                let r = next() % m;
                if r != i && !rows.contains(&r) {
                    rows.push(r);
                    vals.push(if next() % 2 == 0 { 0.125 } else { -0.25 });
                }
            }
            (true, rows, vals)
        }).collect();
        let st = crate::lu_order::factor_structure(m, &cols);
        let mut new = DualFactor::fresh(m, m, &st, &cols);
        let mut old = RefStruct::fresh(m, m, &st, &cols);
        for pass in 0..48 {
            let n = [0, 1, 2, m / 3, m][pass % 5];
            let inp: Vec<_> = (0..n).map(|_| {
                (next() % m, [0.0, -0.0, 1.0, -0.25, TOL, -TOL][next() % 6])
            }).collect();
            let direct = pass % 2 == 0;
            let dense = pass % 3 == 0;
            assert_eq!(bits(new.btran(&inp, direct, dense)), bits(old.btran(&inp, direct, dense)));
            assert_eq!(new.work.to_bits(), old.work.to_bits());
            assert_eq!(new.bt_dense, old.bt_dense);
            let b: Vec<_> = (0..m).map(|_| [0.0, -0.0, 0.5, -2.0][next() % 4]).collect();
            assert_eq!(new.ftran_dense_work(&b).to_bits(), old.ftran_dense_work(&b).to_bits());
            // Perturb the current basis column, keeping a safe nonzero replacement pivot.
            if pass % 4 == 0 {
                let r = next() % m;
                let mut rows = cols[r].1.clone();
                let mut vals = cols[r].2.clone();
                vals[0] += (pass + 1) as f64 * 0.001;
                rows.push((r + 1) % m);
                vals.push(-0.0);
                assert_eq!(bits(new.ftran(&rows, &vals)), bits(old.ftran(&rows, &vals)));
                assert_eq!(new.work.to_bits(), old.work.to_bits());
                new.update(r);
                old.update(r);
                assert_eq!(new.upd_work, old.upd_work);
                assert_eq!(bits(new.diag.iter().copied().enumerate().collect()),
                    bits(old.numerics().0.iter().copied().enumerate().collect()));
                for (a, b) in new.uc.iter().zip(old.u_numerics()) {
                    assert_eq!(bits(a.clone()), bits(b.clone()));
                }
            }
        }
        let mut clone = new.clone();
        assert_eq!(bits(new.btran(&[(0, 1.0)], true, false)), bits(clone.btran(&[(0, 1.0)], true, false)));
    }
}

#[test]
fn captured_transpose_retained_tail_matches_dense_scratch() {
    let records: serde_json::Value =
        serde_json::from_str(include_str!("../../tests/fixtures/ftran_dense_tail.json")).unwrap();
    for rec in records.as_array().unwrap() {
        let diag: Vec<f64> = serde_json::from_value(rec["diag"].clone()).unwrap();
        let m = diag.len();
        let cols: Vec<_> = (0..m).map(|i| (true, vec![i], vec![1.0])).collect();
        let st = crate::lu_order::factor_structure(m, &cols);
        let mut f = DualFactor::fresh(m, m, &st, &cols);
        f.diag = diag;
        f.kind = serde_json::from_value(rec["kind"].clone()).unwrap();
        f.slot = serde_json::from_value(rec["slot"].clone()).unwrap();
        f.i21 = f.kind.iter().filter(|&&k| k != -2).count();
        for (p, &s) in f.slot.iter().enumerate() { f.at[s] = p; }
        let uc: Vec<Vec<(usize, f64)>> = serde_json::from_value(rec["uc"].clone()).unwrap();
        for (c, col) in uc.iter().enumerate() {
            for &(p, v) in col { f.ut[p].push((c, v)); }
        }
        let mut tail = SparseBtranTail::default();
        for h in [2.0, 7.0] {
            f.h33 = h;
            for p in 0..m {
                let inp = [(p, 1.0), ((p + 1) % m, -0.0), (p, -0.5)];
                let mut old_x = vec![0.0; m];
                let mut new_x = vec![0.0; m];
                let expected = f.sparse_btran_u(&inp, &mut old_x, &mut vec![0.0; m], Vec::new());
                tail.reset(m);
                let actual = f.sparse_btran_u(&inp, &mut new_x, &mut tail, Vec::new());
                assert_eq!(actual.0, expected.0);
                assert_eq!(actual.1.to_bits(), expected.1.to_bits());
                assert_eq!(bits(new_x.into_iter().enumerate().collect()), bits(old_x.into_iter().enumerate().collect()));
            }
        }
    }
}

/// The walk's column lengths and the count scan of `update` (O(1) forms) against the retained original on many mixed replacements.
#[test]
fn update_walk_lengths_and_count_scan_match_original() {
    let mut seed = 0x9e3779b97f4a7c15u64;
    let mut next = || {
        seed ^= seed << 13;
        seed ^= seed >> 7;
        seed ^= seed << 17;
        seed as usize
    };
    let mut walks = 0;
    for m in [24, 70, 200] {
        let mut cols: Vec<_> = (0..m).map(|i| {
            let mut rows = vec![i];
            let mut vals = vec![4.0];
            for _ in 0..4 {
                let r = next() % m;
                if !rows.contains(&r) {
                    rows.push(r);
                    vals.push(if next() % 2 == 0 { 0.5 } else { -0.75 });
                }
            }
            (true, rows, vals)
        }).collect();
        let st = crate::lu_order::factor_structure(m, &cols);
        let mut new = DualFactor::fresh(m, m, &st, &cols);
        let mut old = RefStruct::fresh(m, m, &st, &cols);
        for step in 0..(m.min(60)) {
            let r = next() % m;
            let mut rows = vec![r];
            let mut vals = vec![5.0 + step as f64];
            for _ in 0..(1 + next() % 8) {
                let q = next() % m;
                if !rows.contains(&q) {
                    rows.push(q);
                    vals.push([0.25, -0.5, 1.5][next() % 3]);
                }
            }
            assert_eq!(bits(new.ftran(&rows, &vals)), bits(old.ftran(&rows, &vals)));
            new.update(r);
            old.update(r);
            assert_eq!(new.upd_work.len(), old.upd_work.len());
            for (a, b) in new.upd_work.iter().zip(&old.upd_work) {
                assert_eq!((a.0, a.1.to_bits()), (b.0, b.1.to_bits()), "m {m} step {step}");
            }
            assert_eq!(bits(new.diag.iter().copied().enumerate().collect()),
                bits(old.numerics().0.iter().copied().enumerate().collect()));
            walks += new.upd_work.iter().filter(|w| w.0 == 0xf3ca4).count();
            cols[r] = (true, rows, vals);
        }
    }
    assert!(walks > 0, "the walk's jumps were never exercised");
}

/// Every field `update` writes, bitwise (floats as bits, the unordered `empty_tail` as its sorted members).
fn update_state(f: &DualFactor<f64>) -> String {
    let ub = |v: &Vec<Vec<(usize, f64)>>| v.iter().map(|c| c.iter().map(|e| (e.0, e.1.to_bits())).collect::<Vec<_>>()).collect::<Vec<_>>();
    let etas = f.etas.iter().map(|e| (e.idx, e.ind.clone(), e.val.iter().map(|x| x.to_bits()).collect::<Vec<_>>(), e.peel)).collect::<Vec<_>>();
    let work = f.upd_work.iter().map(|w| (w.0, w.1.to_bits(), w.2.to_bits())).collect::<Vec<_>>();
    let diag = f.diag.iter().map(|x| x.to_bits()).collect::<Vec<_>>();
    let et = f.empty_tail.iter().copied().collect::<Vec<_>>();
    format!(
        "{:?}",
        (
            (&f.kind, &f.nuc, f.i21, f.mark, &f.a, &f.b, &f.slot, &f.at, &f.rowl),
            (ub(&f.uc), ub(&f.ut), diag, etas, &f.first_eta, f.fill, &f.row_cap, &f.row_par, et),
            (work, f.sum_eta.to_bits(), f.sum_alpha.to_bits(), f.upd, (f.before_last.0, f.before_last.1.to_bits())),
            (f.bt_occupied, f.dense_est.get().to_bits(), f.spike_used, f.row_cached),
        )
    )
}

/// `update` without per-call temporaries against the retained original (`update_old`) over random replacement sequences, including bare
/// replacements (the entering column is a multiple of the current one), demotions of empty tail columns and nucleus replacements.
#[test]
fn update_without_temporaries_matches_original() {
    let mut seed = 0x1f83d9abfb41bd6bu64;
    let mut next = || {
        seed ^= seed << 13;
        seed ^= seed >> 7;
        seed ^= seed << 17;
        seed as usize
    };
    let (mut bare, mut demoted, mut tail_repl, mut et_seen) = (0, 0, 0, 0);
    for m in [24, 70, 200] {
        for round in 0..16 {
            let mut cols: Vec<_> = (0..m).map(|i| {
                let mut rows = vec![i];
                let mut vals = vec![4.0];
                for _ in 0..(1 + round % 5) {
                    let r = next() % m;
                    if !rows.contains(&r) {
                        rows.push(r);
                        vals.push(if next() % 2 == 0 { 0.5 } else { -0.75 });
                    }
                }
                (true, rows, vals)
            }).collect();
            let st = crate::lu_order::factor_structure(m, &cols);
            let mut new = DualFactor::fresh(m, m, &st, &cols);
            let mut old = new.clone();
            assert_eq!(update_state(&new), update_state(&old));
            for step in 0..(m.min(100)) {
                let r = next() % m;
                let (rows, vals) = if next() % 6 == 0 {
                    (cols[r].1.clone(), cols[r].2.iter().map(|v| 2.0 * v).collect())
                } else {
                    let mut rows = vec![r];
                    let mut vals = vec![5.0 + step as f64];
                    for _ in 0..(1 + next() % 8) {
                        let q = next() % m;
                        if !rows.contains(&q) {
                            rows.push(q);
                            vals.push([0.25, -0.5, 1.5][next() % 3]);
                        }
                    }
                    (rows, vals)
                };
                let a = new.ftran(&rows, &vals);
                let b = old.ftran(&rows, &vals);
                assert_eq!(bits(a), bits(b));
                bare += usize::from(new.spike.iter().all(|e| e.0 == r));
                tail_repl += usize::from(new.kind[r] == -2);
                et_seen += usize::from(!new.empty_tail.list.is_empty() && new.kind[r] != -2);
                let kinds = new.kind.clone();
                new.update(r);
                old.update_old(r);
                demoted += (0..m).filter(|&x| x != r && kinds[x] == -2 && new.kind[x] != -2).count();
                assert_eq!(update_state(&new), update_state(&old), "m {m} round {round} step {step}");
                cols[r] = (true, rows, vals);
            }
        }
    }
    assert!(bare > 0 && demoted > 0 && tail_repl > 0, "bare {bare} demoted {demoted} tail {tail_repl} et {et_seen}");
}

/// `uorder` (bit-ordered lists, slot heap, retained lists) against the previous implementation after fresh factors and random update
/// sequences, on the solve's own L-stage lists and on random lists with repeats and random zero patterns.
#[test]
fn uorder_matches_original_after_updates() {
    let mut seed = 0x9e3779b97f4a7c15u64;
    let mut next = || {
        seed ^= seed << 13;
        seed ^= seed >> 7;
        seed ^= seed << 17;
        seed as usize
    };
    let (mut dense_seen, mut sparse_seen, mut sorted_seen, mut checks) = (0, 0, 0, 0);
    for m in [24, 70, 200, 700] {
        for round in 0..10 {
            let mut cols: Vec<_> = (0..m).map(|i| {
                let mut rows = vec![i];
                let mut vals = vec![4.0];
                for _ in 0..(1 + round % 5) {
                    let r = next() % m;
                    if !rows.contains(&r) {
                        rows.push(r);
                        vals.push(if next() % 2 == 0 { 0.5 } else { -0.75 });
                    }
                }
                (true, rows, vals)
            }).collect();
            let st = crate::lu_order::factor_structure(m, &cols);
            let mut f = DualFactor::fresh(m, m, &st, &cols);
            for step in 0..m.min(60) + 1 {
                let mut lists: Vec<(Vec<usize>, Vec<f64>)> = Vec::new();
                for k in 0..6 {
                    let n = [1, 3, 9, m / 8 + 1, m / 3 + 1, m][(k + step) % 6].min(3 * m);
                    let inp: Vec<usize> = (0..n).map(|_| next() % m).collect();
                    let x: Vec<f64> = (0..m).map(|_| [0.0, -0.0, 1.0, 0.5][next() % 4]).collect();
                    lists.push((inp, x));
                }
                lists.push((f.f2.get_mut().col_in.clone(), f.f2.get_mut().x_col.clone()));
                for (inp, x) in lists {
                    if x.len() != m {
                        continue;
                    }
                    let want = f.uorder_ref(&inp, &x);
                    let got = f.uorder(&inp, &x);
                    assert_eq!(got, want, "m {m} round {round} step {step}");
                    f.uorder_done(got);
                    let n = {
                        let mut s = inp.clone();
                        s.sort_unstable();
                        s.dedup();
                        s.len()
                    };
                    if (m as f64) <= (n as f64) * f.dens { dense_seen += 1 } else { sparse_seen += 1 }
                    sorted_seen += usize::from(want.len() >= 24);
                    checks += 1;
                }
                let r = next() % m;
                let mut rows = vec![r];
                let mut vals = vec![5.0 + step as f64];
                for _ in 0..(1 + next() % 8) {
                    let q = next() % m;
                    if !rows.contains(&q) {
                        rows.push(q);
                        vals.push([0.25, -0.5, 1.5][next() % 3]);
                    }
                }
                f.ftran(&rows, &vals);
                f.update(r);
                cols[r] = (true, rows, vals);
            }
        }
    }
    assert!(dense_seen > 0 && sparse_seen > 0 && sorted_seen > 0, "{dense_seen} {sparse_seen} {sorted_seen} {checks}");
}

/// The FTRAN work estimates booked inside the U loops (`u_work`, `tau_u_finish`, `third_work`, `dense_pair_work`) against the retained scans
/// (`*_ref`, asserted inline under cfg(test)) and against the original backend, over random factors with updates, sparse and dense inputs,
/// tiny and zero values, both vectors sparse or dense.
#[test]
fn ftran_work_booked_in_u_loops_matches_scans() {
    let mut seed = 0x2545f4914f6cdd1du64;
    let mut next = || {
        seed ^= seed << 13;
        seed ^= seed >> 7;
        seed ^= seed << 17;
        seed as usize
    };
    // joint selections seen: (third, dense pair, usual pair), tau dense/sparse, column dense/sparse
    let (mut third, mut pair, mut usual, mut checks) = (0, 0, 0, 0);
    let (mut tau_dense_seen, mut col_dense_seen, mut tiny_seen) = (0, 0, 0);
    for m in [24, 70, 200, 700, 1500] {
        for round in 0..10 {
            let cols: Vec<_> = (0..m).map(|i| {
                let mut rows = vec![i];
                let mut vals = vec![4.0];
                for _ in 0..(1 + round % 5) {
                    let r = next() % m;
                    if !rows.contains(&r) {
                        rows.push(r);
                        vals.push(if next() % 2 == 0 { 0.5 } else { -0.75 });
                    }
                }
                (true, rows, vals)
            }).collect();
            let st = crate::lu_order::factor_structure(m, &cols);
            let mut new = DualFactor::fresh(m, m, &st, &cols);
            let mut old = RefStruct::fresh(m, m, &st, &cols);
            for step in 0..m.min(40) {
                let kind = next() % 6;
                let n = [1, 2, 4, m / 6 + 1, m / 2 + 1, m][(step + kind) % 6].min(m);
                let mut rows: Vec<usize> = Vec::new();
                for _ in 0..n {
                    let r = next() % m;
                    if !rows.contains(&r) {
                        rows.push(r);
                    }
                }
                let val = |k: usize| match k % 9 {
                    0 => 1e-15,
                    1 => -1e-14,
                    2 => 0.0,
                    3 => -0.0,
                    4 => 1e-9,
                    5 => 0.25,
                    6 => -1.5,
                    7 => 3.0,
                    _ => 1e-12,
                };
                let vals: Vec<f64> = rows.iter().map(|_| val(next())).collect();
                let rho: Vec<usize> = {
                    let mut r: Vec<usize> = Vec::new();
                    for _ in 0..[1, 3, m / 5 + 1, m / 2 + 1, m][next() % 5].min(m) {
                        let q = next() % m;
                        if !r.contains(&q) {
                            r.push(q);
                        }
                    }
                    r
                };
                let rvals: Vec<f64> = rho.iter().map(|_| val(next())).collect();
                let force_dense = next() % 7 == 0;
                let (a, b) = if force_dense {
                    (new.ftran_basic(&rows, &vals), old.ftran_basic(&rows, &vals))
                } else {
                    (new.ftran(&rows, &vals), old.ftran(&rows, &vals))
                };
                assert_eq!(bits(a.clone()), bits(b));
                assert_eq!(new.work.to_bits(), old.work.to_bits());
                assert_eq!((new.work_parts.0.to_bits(), new.work_parts.1.to_bits()), (old.work_parts.0.to_bits(), old.work_parts.1.to_bits()));
                assert_eq!(new.last_dense, old.last_dense());
                let needed: Vec<usize> = a.iter().map(|e| e.0).collect();
                let dense_in = next() % 5 == 0;
                let ta = new.ftran_needed(&rho, &rvals, &needed, dense_in);
                let tb = old.ftran_needed(&rho, &rvals, &needed, dense_in);
                assert_eq!(bits(ta), bits(tb));
                assert_eq!(new.tau_l.get().to_bits(), old.tau_l.get().to_bits());
                assert_eq!(new.tau_u.get().to_bits(), old.tau_u.get().to_bits());
                assert_eq!(new.tau_fin.get(), old.tau_fin.get());
                let (ja, jb) = (new.joint_u(&needed), old.joint_u(&needed));
                assert_eq!(ja.map(f64::to_bits), jb.map(f64::to_bits), "m {m} round {round} step {step}");
                {
                    let f = new.f2.borrow();
                    tau_dense_seen += usize::from(f.tau_n.is_none());
                    col_dense_seen += usize::from(f.col_n.is_none());
                }
                tiny_seen += usize::from(new.f2.borrow().x_tau.iter().any(|v| *v != 0.0 && v.abs() <= TOL));
                match ja {
                    None => usual += 1,
                    Some(_) if new.f2.borrow().col_n.is_some() => third += 1,
                    Some(_) => pair += 1,
                }
                checks += 1;
                let b: Vec<_> = (0..m).map(|_| [0.0, -0.0, 0.5, -2.0][next() % 4]).collect();
                assert_eq!(new.ftran_dense_work(&b).to_bits(), old.ftran_dense_work(&b).to_bits());
                if next() % 2 == 0 {
                    let r = next() % m;
                    let mut rows = vec![r];
                    let mut vals = vec![5.0 + step as f64];
                    for _ in 0..(1 + next() % 8) {
                        let q = next() % m;
                        if !rows.contains(&q) {
                            rows.push(q);
                            vals.push([0.25, -0.5, 1.5][next() % 3]);
                        }
                    }
                    assert_eq!(bits(new.ftran(&rows, &vals)), bits(old.ftran(&rows, &vals)));
                    new.update(r);
                    old.update(r);
                }
            }
        }
    }
    assert!(third > 0 && pair > 0 && usual > 0, "third {third} pair {pair} usual {usual}");
    assert!(tau_dense_seen > 0 && col_dense_seen > 0 && tiny_seen > 0, "tau dense {tau_dense_seen} col dense {col_dense_seen} tiny {tiny_seen} checks {checks}");
}
