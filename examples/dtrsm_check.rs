//! Reads `m n diag` lines followed by hex words (64 bit) of A (m x m), B (m x n) and the library result of DTRSM('L','U','T',diag) (`tools/mkl_probe dtrsm`);
//! reports the shapes where `rose::chol::dtrsm_llu_col` differs bit for bit (diag U only).
use std::io::Read;

fn main() {
    let mut s = String::new();
    std::io::stdin().read_to_string(&mut s).unwrap();
    let mut it = s.split_whitespace();
    let (mut total, mut bad_shapes) = (0usize, 0usize);
    while let Some(t) = it.next() {
        let m: usize = t.parse().unwrap();
        let n: usize = it.next().unwrap().parse().unwrap();
        let _diag = it.next().unwrap();
        let mut rd = |len: usize| -> Vec<f64> {
            (0..len)
                .map(|_| f64::from_bits(u64::from_str_radix(it.next().unwrap(), 16).unwrap()))
                .collect()
        };
        let a = rd(m * m);
        let b = rd(m * n);
        let r = rd(m * n);
        let mut l = vec![0.0f64; m * m];
        for i in 0..m {
            for j in 0..i {
                l[i * m + j] = a[j + m * i];
            }
        }
        let (mut bad, mut first) = (0, None);
        for c in 0..n {
            let x = rose::chol::dtrsm_llu_col(&l, m, &b[c * m..(c + 1) * m], c, n);
            for i in 0..m {
                if x[i].to_bits() != r[c * m + i].to_bits() {
                    bad += 1;
                    if first.is_none() {
                        first = Some((i, c));
                    }
                }
            }
        }
        total += 1;
        if bad > 0 {
            bad_shapes += 1;
            println!("({m},{n}): {bad} of {} differ, first {:?}", m * n, first);
        }
    }
    println!("shapes {total}, differing {bad_shapes}");
}
