//! Differential tests of the BTRAN output packers against their previous filter/collect forms.
use super::*;

struct Rng(u64);
impl Rng {
    fn next(&mut self) -> u64 {
        self.0 ^= self.0 << 13;
        self.0 ^= self.0 >> 7;
        self.0 ^= self.0 << 17;
        self.0
    }
    fn below(&mut self, n: usize) -> usize {
        (self.next() % n as u64) as usize
    }
}

fn value(rng: &mut Rng) -> f64 {
    match rng.below(8) {
        0 | 1 | 2 => 0.0,
        3 => -0.0,
        4 => f64::NAN,
        5 => f64::from_bits(0x7ff8_0000_0000_0000 | (rng.next() & 0xffff) | 1 << 63),
        6 => f64::MIN_POSITIVE / 4.0,
        _ => (rng.next() % 1000) as f64 / 7.0 - 70.0,
    }
}

fn bits(v: &[(usize, f64)]) -> Vec<(usize, u64)> {
    v.iter().map(|&(i, x)| (i, x.to_bits())).collect()
}

#[test]
fn dense_output_matches_filter_collect() {
    let mut rng = Rng(0x9e37_79b9_7f4a_7c15);
    for m in [0usize, 1, 2, 3, 7, 64, 333] {
        for _ in 0..40 {
            let mut perm: Vec<usize> = (0..m).collect();
            for k in (1..m).rev() {
                perm.swap(k, rng.below(k + 1));
            }
            let x: Vec<f64> = (0..m + rng.below(3)).map(|_| value(&mut rng)).collect();
            let m = perm.len();
            let want: Vec<(usize, f64)> =
                (0..m).filter_map(|r| (x[perm[r]] != 0.0).then_some((r, x[perm[r]]))).collect();
            let bound = perm.iter().map(|&p| p + 1).max().unwrap_or(0);
            let got = DualFactor::<f64>::dense_btran_output(&x, &perm, bound, vec![(9, 9.0); 3]);
            assert_eq!(bits(&got), bits(&want));
        }
    }
}

#[test]
fn sparse_output_matches_filter_collect() {
    let mut rng = Rng(0x1234_5678_9abc_def1);
    for m in [1usize, 2, 5, 64, 300] {
        for _ in 0..60 {
            let x: Vec<f64> = (0..m).map(|_| value(&mut rng)).collect();
            let row_of: Vec<usize> = (0..m).map(|_| rng.below(m)).collect();
            // duplicates and any order, as the level walk may seed
            let res: Vec<usize> = (0..rng.below(2 * m + 1)).map(|_| rng.below(m)).collect();
            let want: Vec<(usize, f64)> = res.iter().copied().filter(|&i| x[i] != 0.0).map(|i| (row_of[i], x[i])).collect();
            let got = DualFactor::<f64>::sparse_btran_output(&x, &res, &row_of, vec![(7, 7.0); 5]);
            assert_eq!(bits(&got), bits(&want));
            let xf: Vec<crate::f80::F80> = x.iter().map(|&v| crate::f80::F80::from_f64(v)).collect();
            let gotf = DualFactor::<crate::f80::F80>::sparse_btran_output(&xf, &res, &row_of, Vec::new());
            let wantf: Vec<(usize, crate::f80::F80)> =
                res.iter().copied().filter(|&i| xf[i] != crate::f80::F80::ZERO).map(|i| (row_of[i], xf[i])).collect();
            assert_eq!(gotf.len(), wantf.len());
            for (a, b) in gotf.iter().zip(&wantf) {
                assert!(a.0 == b.0 && a.1.repr() == b.1.repr());
            }
        }
    }
}

/// `btran_rho`'s compaction against the previous pipeline (`dense_btran_output`/`sparse_btran_output`, then the `|v| > tol` copy), and `rho_full`
/// against the packers' unfiltered result.
fn rho_output_matches<R: Real>(key: impl Fn(R) -> u128) {
    let mut rng = Rng(0x0dd5_eed5_f00d_1234);
    let tol = R::c(1e-13);
    // The tolerance edge, signed zeros, NaN (both signs), and a normal range.
    let edge = [1e-13, -1e-13, 1e-13 * (1.0 + 1e-15), 1e-13 * (1.0 - 1e-15), 5e-14, -5e-14, 0.0, -0.0, f64::NAN, -f64::NAN, 0.5, -3.0];
    let mut dense_calls = 0;
    for m in [0usize, 1, 2, 3, 7, 64, 333] {
        for call in 0..60 {
            let x: Vec<R> = (0..m)
                .map(|_| R::c(if rng.below(3) == 0 { edge[rng.below(edge.len())] } else { value(&mut rng) }))
                .collect();
            let mut perm: Vec<usize> = (0..m).collect();
            for k in (1..m).rev() {
                perm.swap(k, rng.below(k + 1));
            }
            let row_of: Vec<usize> = (0..m).map(|_| rng.below(m)).collect();
            let res: Vec<usize> = (0..if m == 0 { 0 } else { rng.below(2 * m + 1) }).map(|_| rng.below(m)).collect();
            let dense = call % 2 == 0;
            let (full, mut ws) = (
                if dense {
                    DualFactor::<R>::dense_btran_output(&x, &perm, m, vec![(9, R::ONE); 2])
                } else {
                    DualFactor::<R>::sparse_btran_output(&x, &res, &row_of, vec![(9, R::ONE); 2])
                },
                BtranWork::<R>::default(),
            );
            let want: Vec<(usize, R)> = full.iter().copied().filter(|e| e.1.abs() > tol).collect();
            let (mut rows, mut vals) = (vec![5usize; 3], vec![R::ONE; 3]);
            rows.clear();
            vals.clear();
            let o = RhoOut { rows: &mut rows, vals: &mut vals, tol };
            if dense {
                DualFactor::<R>::rho_output(m, |k| (k, x[perm[k]]), o, &mut ws);
                dense_calls += 1;
            } else {
                DualFactor::<R>::rho_output(res.len(), |k| (row_of[res[k]], x[res[k]]), o, &mut ws);
            }
            assert_eq!(rows.len(), vals.len());
            assert_eq!(rows.len(), want.len(), "m {m} call {call}");
            for (k, w) in want.iter().enumerate() {
                assert!(rows[k] == w.0 && key(vals[k]) == key(w.1), "m {m} call {call} entry {k}");
            }
            // The unfiltered list rebuilt from the kept entries and the dropped ones.
            let rebuilt = DualFactor::<R>::merge_dropped(&rows, &vals, &ws.dropped);
            assert_eq!(rebuilt.len(), full.len(), "m {m} call {call}");
            for (a, b) in rebuilt.iter().zip(&full) {
                assert!(a.0 == b.0 && key(a.1) == key(b.1), "m {m} call {call}");
            }
        }
    }
    assert!(dense_calls > 100);
}

#[test]
fn rho_output_matches_packer_and_filter() {
    rho_output_matches::<f64>(|v| v.to_bits() as u128);
    rho_output_matches::<crate::f80::F80>(|v| {
        let (a, b) = v.repr();
        ((b as u128) << 64) | a as u128
    });
}
