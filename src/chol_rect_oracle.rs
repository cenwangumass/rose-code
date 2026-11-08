// Complete original factor_block frozen before RHS batching.
use super::*;
pub(super) fn factor_block(
    mm: &mut [f64],
    r: &mut [f64],
    w: usize,
    c: usize,
    orig: &[f64],
    reg: Option<(usize, usize)>,
) -> (Vec<f64>, Vec<f64>) {
    let dense = (1..w).any(|i| (0..i).any(|j| mm[i * w + j] != 0.0));
    let (d, rcp) = if dense {
        dense_ldl(mm, w, orig, reg)
    } else {
        diag_only(mm, w, orig, reg)
    };
    let rows = c - w;
    if rows != 0 {
        if dense {
            let big32 = SINGLE.with(|c| c.get()) && w >= 2;
            let l32: Vec<f32> = if big32 {
                mm.iter().map(|&v| v as f32).collect()
            } else {
                Vec::new()
            };
            for rr in 0..rows {
                let x = if big32 {
                    let mut x: Vec<f32> =
                        r[rr * w..(rr + 1) * w].iter().map(|&v| v as f32).collect();
                    crate::sgemm32::strsm_llu_col(&l32, w, &mut x, rr, rows);
                    x.iter().map(|&v| v as f64).collect()
                } else if SINGLE.with(|c| c.get()) {
                    strsm_seq(mm, w, &r[rr * w..(rr + 1) * w])
                } else {
                    dtrsm_llu_col(mm, w, &r[rr * w..(rr + 1) * w], rr, rows)
                };
                r[rr * w..(rr + 1) * w].copy_from_slice(&x);
            }
        }
        for rr in 0..rows {
            for k in 0..w {
                r[rr * w + k] = rf(r[rr * w + k] * if rcp[k] != 0.0 { rf(rcp[k]) } else { 1.0 });
            }
        }
    }
    (d, rcp)
}

