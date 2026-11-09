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
