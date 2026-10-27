//! Reads records `m n k lda ldb ldc` followed by hex words of A (lda x m), B (ldb x n), C before and C after (ldc x n) of captured SGEMM('T','N') calls;
//! reports the bit mismatches of `rose::sgemm32::sgemm_t_n` (only the m x n part of C).
use std::io::Read;

fn main() {
    let mut s = String::new();
    std::io::stdin().read_to_string(&mut s).unwrap();
    let mut it = s.split_whitespace();
    let (mut total, mut bad_shapes) = (0usize, 0usize);
    while let Some(t) = it.next() {
        let m: usize = t.parse().unwrap();
        let mut nx = || -> usize { it.next().unwrap().parse().unwrap() };
        let (n, k, lda, ldb, ldc) = (nx(), nx(), nx(), nx(), nx());
        let mut rd = |len: usize| -> Vec<f32> {
            (0..len)
                .map(|_| f32::from_bits(u32::from_str_radix(it.next().unwrap(), 16).unwrap()))
                .collect()
        };
        let a = rd(lda * m);
        let b = rd(ldb * n);
        let mut c = rd(ldc * n);
        let r = rd(ldc * n);
        rose::sgemm32::sgemm_t_n(m, n, k, &a, lda, &b, ldb, &mut c, ldc);
        let mut bad = 0;
        let mut first = None;
        for j in 0..n {
            for i in 0..m {
                if c[j * ldc + i].to_bits() != r[j * ldc + i].to_bits() {
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
                "({m},{n},{k}) ldc {ldc}: {bad} of {} differ, first {:?}",
                m * n,
                first
            );
        }
    }
    println!("calls {total}, differing {bad_shapes}");
}
