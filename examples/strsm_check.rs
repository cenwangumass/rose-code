//! Reads `m n diag` lines from stdin followed by the hex words of A (m x m column-major, lda = m), B (m x n) and the library result of
//! STRSM('L','U','T',diag), as written by `tools/mkl_probe trsm`; reports the bit mismatches of `rose::sgemm32::strsm_llu_col`.
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
        let mut rd = |len: usize| -> Vec<f32> {
            (0..len)
                .map(|_| f32::from_bits(u32::from_str_radix(it.next().unwrap(), 16).unwrap()))
                .collect()
        };
        let a = rd(m * m);
        let b = rd(m * n);
        let r = rd(m * n);
        // column-major A(j, i) (j < i) = a[j + m i] is the lower entry (i, j) of the row-major l
        let mut l = vec![0.0f32; m * m];
        for i in 0..m {
            for j in 0..i {
                l[i * m + j] = a[j + m * i];
            }
        }
        let mut bad = 0;
        let mut first = None;
        for c in 0..n {
            let mut x = b[c * m..(c + 1) * m].to_vec();
            if true {
                rose::sgemm32::strsm_llu_col(&l, m, &mut x, c, n);
            }
            for i in 0..m {
                if x[i].to_bits() != r[c * m + i].to_bits() {
                    bad += 1;
                    if std::env::var("VERB").is_ok() {
                        eprintln!(
                            "diff row {i} col {c}: ref {:e} rose {:e}",
                            r[c * m + i],
                            x[i]
                        );
                    }
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
