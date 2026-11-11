//! Differential tests of the A/B lists with coefficients: the in-place edits against `Vec` semantics, random update sequences against the
//! retained original `update_old` (lists, order and coefficient bits after every step, every branch of `update`), and the needed tail
//! scatter over A against the scan of `uc` it replaced.
use super::ab::{same_bits, AbLists};
use super::*;

fn rng(mut seed: u64) -> impl FnMut() -> usize {
    move || {
        seed ^= seed << 13;
        seed ^= seed >> 7;
        seed ^= seed << 17;
        seed as usize
    }
}

fn coeffs() -> [f64; 5] {
    [0.0, -0.0, 0.5, -3.0, 1e-300]
}

/// Random edits through every primitive against the `Vec` operations they stand for (A and B as nested vectors).
fn ab_ops<R: Real>() {
    let mut next = rng(0x2545f4914f6cdd1d);
    let m = 12;
    let mut ab = AbLists::<R>::new();
    ab.reset(m);
    let mut model: Vec<(Vec<(u32, R)>, Vec<(u32, R)>)> = vec![(Vec::new(), Vec::new()); m];
    let mut seen = [0usize; 5];
    for step in 0..4000 {
        let c = next() % m;
        let (na, nb) = (model[c].0.len(), model[c].1.len());
        match next() % 5 {
            0 => {
                let mk = |n: usize, next: &mut dyn FnMut() -> usize, base: u32| -> Vec<(u32, R)> {
                    (0..n).map(|i| (base + i as u32 * 3 + (next() % 2) as u32, R::c(coeffs()[next() % 5]))).collect()
                };
                let (a, b) = (mk(next() % 6, &mut next, 0), mk(next() % 6, &mut next, 100));
                ab.set(c, &a, &b);
                model[c] = (a, b);
                seen[0] += 1;
            }
            1 if na > 0 => {
                // A.swap_remove(p) then B.insert(0, removed)
                let p = next() % na;
                let (a, b) = &mut model[c];
                let e = a.swap_remove(p);
                b.insert(0, e);
                ab.demote(c, p);
                seen[1] += 1;
            }
            2 if na > 0 => {
                // A.swap_remove(p), B's last entry to its front
                let p = next() % na;
                let (a, b) = &mut model[c];
                a.swap_remove(p);
                if let Some(last) = b.pop() {
                    b.insert(0, last);
                }
                ab.drop_a(c, p);
                seen[2] += 1;
            }
            3 if nb > 0 => {
                let p = next() % nb;
                model[c].1.swap_remove(p);
                ab.drop_b(c, p);
                seen[3] += 1;
            }
            4 if nb > 0 => {
                // remove_front_fill from B, push to A
                let p = next() % nb;
                let (a, b) = &mut model[c];
                let e = b[p];
                if p != 0 {
                    b[p] = b[0];
                }
                b.remove(0);
                a.push(e);
                ab.b_to_a(c, p);
                seen[4] += 1;
            }
            _ => {}
        }
        for (c, (a, b)) in model.iter().enumerate() {
            assert_eq!(ab.a(c), a.iter().map(|e| e.0).collect::<Vec<_>>(), "A of {c} at step {step}");
            assert_eq!(ab.b(c), b.iter().map(|e| e.0).collect::<Vec<_>>(), "B of {c} at step {step}");
            assert!(ab.a_val(c).iter().zip(a).all(|(&x, e)| same_bits(x, e.1)), "A coefficients of {c} at step {step}");
            assert!(ab.b_val(c).iter().zip(b).all(|(&x, e)| same_bits(x, e.1)), "B coefficients of {c} at step {step}");
            assert_eq!((ab.na(c), ab.nb(c)), (a.len(), b.len()));
        }
    }
    assert!(seen.iter().all(|&n| n > 100), "{seen:?}");
}

#[test]
fn ab_edits_match_vec_semantics() {
    ab_ops::<f64>();
    ab_ops::<crate::f80::F80>();
}

/// Lists, order and coefficient bits of every column of two factors.
fn assert_same_lists<R: Real>(new: &DualFactor<R>, old: &DualFactor<R>, what: &str) {
    assert_eq!(new.ab.to_vecs(), old.ab.to_vecs(), "{what}: lists");
    for c in 0..new.m {
        for (x, y) in [(new.ab.a_val(c), old.ab.a_val(c)), (new.ab.b_val(c), old.ab.b_val(c))] {
            assert!(x.iter().zip(y).all(|(&p, &q)| same_bits(p, q)), "{what}: coefficients of column {c}");
        }
    }
}

/// Random replacement sequences (the cases of `update_without_temporaries_matches_original`, with a triangular start so that columns
/// empty and demote, and bare replacements) on f64 and F80, with the branch of each step counted from its work bookings.
fn update_sequences<R: Real>(seen: &mut std::collections::BTreeMap<&'static str, usize>) {
    let mut next = rng(0x6a09e667f3bcc908);
    for m in [18, 60, 150] {
        for round in 0..14 {
            let mut cols: Vec<(bool, Vec<usize>, Vec<f64>)> = (0..m)
                .map(|i| {
                    let mut rows = vec![i];
                    let mut vals = vec![4.0];
                    for _ in 0..(1 + round % 5) {
                        let r = if round % 2 == 1 {
                            if i == 0 {
                                continue;
                            }
                            next() % i
                        } else {
                            next() % m
                        };
                        if !rows.contains(&r) {
                            rows.push(r);
                            vals.push(coeffs()[2 + next() % 3] + 0.25);
                        }
                    }
                    (true, rows, vals)
                })
                .collect();
            let st = crate::lu_order::factor_structure_r::<R>(m, &cols);
            let mut new = DualFactor::<R>::fresh(m, m, &st, &cols);
            let mut old = new.clone();
            assert_same_lists(&new, &old, "fresh");
            for step in 0..m.min(90) {
                let r = next() % m;
                let (rows, vals): (Vec<usize>, Vec<f64>) = if next() % 6 == 0 {
                    (cols[r].1.clone(), cols[r].2.iter().map(|v| 2.0 * v).collect())
                } else {
                    let mut rows = vec![r];
                    let mut vals = vec![5.0 + step as f64];
                    for _ in 0..(1 + next() % 8) {
                        let q = next() % m;
                        if !rows.contains(&q) {
                            rows.push(q);
                            vals.push([0.25, -0.5, 1.5, -0.0][next() % 4]);
                        }
                    }
                    (rows, vals)
                };
                let rv: Vec<R> = vals.iter().map(|&v| R::c(v)).collect();
                new.ftran(&rows, &rv);
                old.ftran(&rows, &rv);
                new.update(r);
                old.update_old(r);
                for w in &new.upd_work {
                    let name = match w.0 {
                        0xf3eac => "demote",
                        0xf4073 => "keep",
                        0xf2969 => "eliminate from A",
                        0xf3f71 => "eliminate from B",
                        0xf3ca4 => "walk jumps",
                        0xf3d83 => "walk shift",
                        _ => continue,
                    };
                    *seen.entry(name).or_default() += 1;
                }
                *seen.entry("bare").or_default() += usize::from(new.spike.iter().all(|e| e.0 == r));
                assert_same_lists(&new, &old, &format!("m {m} round {round} step {step}"));
                new.ab.check(&new.uc, &new.kind);
                cols[r] = (true, rows, vals);
            }
            new.recycle();
        }
    }
}

#[test]
fn ab_lists_follow_every_update_branch() {
    let mut seen = std::collections::BTreeMap::new();
    update_sequences::<f64>(&mut seen);
    update_sequences::<crate::f80::F80>(&mut seen);
    for name in ["demote", "keep", "eliminate from A", "eliminate from B", "walk jumps", "walk shift", "bare"] {
        assert!(seen.get(name).copied().unwrap_or(0) > 0, "branch {name} never taken: {seen:?}");
    }
}

/// The needed tail scatter over A against the `uc` scan with the `kind` filter it replaced: same visited rows, same bits of `x`
/// (signed zeros included), on fresh and updated factors.
fn scatter<R: Real>() {
    let mut next = rng(0xbb67ae8584caa73b);
    let mut visited = 0usize;
    for m in [20, 90] {
        for round in 0..8 {
            let cols: Vec<(bool, Vec<usize>, Vec<f64>)> = (0..m)
                .map(|i| {
                    let mut rows = vec![i];
                    let mut vals = vec![3.0];
                    for _ in 0..(1 + round % 4) {
                        let r = if round % 2 == 1 { if i == 0 { continue } else { next() % i } } else { next() % m };
                        if !rows.contains(&r) {
                            rows.push(r);
                            vals.push(coeffs()[2 + next() % 3]);
                        }
                    }
                    (true, rows, vals)
                })
                .collect();
            let st = crate::lu_order::factor_structure_r::<R>(m, &cols);
            let mut f = DualFactor::<R>::fresh(m, m, &st, &cols);
            for step in 0..30 {
                if step > 0 {
                    let r = next() % m;
                    let rows = vec![r, next() % m, next() % m];
                    let rows: Vec<usize> = rows.iter().enumerate().filter(|&(i, x)| rows[..i].iter().all(|y| y != x)).map(|(_, &x)| x).collect();
                    let vals: Vec<R> = rows.iter().map(|&q| R::c(if q == r { 6.0 } else { -0.5 })).collect();
                    f.ftran(&rows, &vals);
                    f.update(r);
                }
                let x0: Vec<R> = (0..m).map(|_| R::c(coeffs()[next() % 5] + (next() % 3) as f64)).collect();
                for c in (0..m).filter(|&c| f.kind[c] == -2) {
                    let xc = R::c(coeffs()[2 + next() % 3]);
                    let (mut old_x, mut new_x) = (x0.clone(), x0.clone());
                    let mut old_rows = Vec::new();
                    for &(j, v) in &f.uc[c] {
                        if f.kind[j] == -2 {
                            old_x[j] -= v * xc;
                            old_rows.push(j);
                        }
                    }
                    let mut new_rows = Vec::new();
                    for (&j, &v) in f.ab.a(c).iter().zip(f.ab.a_val(c)) {
                        new_x[j as usize] -= v * xc;
                        new_rows.push(j as usize);
                    }
                    old_rows.sort_unstable();
                    new_rows.sort_unstable();
                    assert_eq!(old_rows, new_rows, "visited rows of column {c}");
                    assert!(old_x.iter().zip(&new_x).all(|(&p, &q)| same_bits(p, q)), "x after column {c}");
                    visited += new_rows.len();
                }
            }
        }
    }
    assert!(visited > 100);
}

#[test]
fn needed_scatter_over_a_matches_the_filtered_uc_scan() {
    scatter::<f64>();
    scatter::<crate::f80::F80>();
}
