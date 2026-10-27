//! Reads `m n k` lines from stdin followed by the hex words of A (k x m), B (k x n) and the library result C (m x n), column-major with lda = ldb = k,
//! ldc = m; reports the bit mismatches of `rose::sgemm32`. Input written by `tools/mkl_probe`.
use std::io::Read;

fn main() {
    let mut s = String::new();
    std::io::stdin().read_to_string(&mut s).unwrap();
    let mut it = s.split_whitespace();
    let (mut total, mut bad_shapes) = (0usize, 0usize);
    while let Some(t) = it.next() {
        let m: usize = t.parse().unwrap();
        let n: usize = it.next().unwrap().parse().unwrap();
        let k: usize = it.next().unwrap().parse().unwrap();
        let mut rd = |len: usize| -> Vec<f32> {
            (0..len)
                .map(|_| f32::from_bits(u32::from_str_radix(it.next().unwrap(), 16).unwrap()))
                .collect()
        };
        let a = rd(k * m);
        let b = rd(k * n);
        let r = rd(m * n);
        let mut c = vec![0.0f32; m * n];
        rose::sgemm32::sgemm_t_n(m, n, k, &a, k, &b, k, &mut c, m);
        let mut bad = 0;
        let mut first = None;
        for j in 0..n {
            for i in 0..m {
                if c[j * m + i].to_bits() != r[j * m + i].to_bits() {
                    bad += 1;
                    if first.is_none() {
                        first = Some((i, j));
                    }
                }
            }
        }
        total += 1;
        if bad > 0 {
            bad_shapes += 1;
            println!(
                "({m},{n},{k}): {bad} of {} differ, first {:?}",
                m * n,
                first
            );
        }
    }
    println!("shapes {total}, differing {bad_shapes}");
}
