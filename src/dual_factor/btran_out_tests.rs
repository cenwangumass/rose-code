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
