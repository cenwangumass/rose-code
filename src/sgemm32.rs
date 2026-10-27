//! Bit model of single-precision `SGEMM('T', 'N', alpha = 1, beta = 1)` (`0xa009a0`, MKL `avx2`): `C(i,j) += sum_q A(q,i) B(q,j)`, column-major,
//! `A(q,i) = a[i lda + q]`, `B(q,j) = b[j ldb + q]`. Per-element accumulation order, recovered by running the library on random data
//! (`tools/mkl_probe.c`) and checked bit for bit against it for every shape the factorisation produces.
//!
//! Paths (`mkl_blas_avx2_xsgemm`): `m == 1` or `n == 1`: SGEMV; `m >= 56 && n >= 56`: copy-based `sgemm_kernel_0` (A in panels of 24 rows,
//! then the bits 16, 8, 4, 2, 1 of `m % 24`; B in columns kernels of 4, then a 2 and a 1); otherwise `sgemm_kernel_nocopy_TN_b1` (rows in blocks of 16,
//! a last block of 9..15 rows is a masked 16-row block, 1..8 rows an 8-row block; columns in blocks of 6 per driver column block with remainders).

#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum Mode {
    /// One accumulator, fma over all terms from 0.
    Chain,
    /// Accumulators of the even and the odd terms over the first `2 (k / 2)` terms (fma), their sum, then the remaining terms fma'd into the sum.
    Par,
    /// As `Par` with rounded products and rounded sums.
    ParU,
    /// As `Par` over the first `4 (k / 4)` terms (the copy-based kernel).
    Par4,
}

fn rd(x: f32) -> f32 {
    std::hint::black_box(x)
}

/// Sum of the products `a[q] b[q]`, `q < k`, in the order of `mode`.
pub fn acc(mode: Mode, k: usize, a: &[f32], b: &[f32]) -> f32 {
    match mode {
        Mode::Chain => {
            let mut s = 0.0f32;
            for q in 0..k {
                s = a[q].mul_add(b[q], s);
            }
            s
        }
        Mode::Par | Mode::Par4 => {
            let m4 = if mode == Mode::Par {
                k - k % 2
            } else {
                k - k % 4
            };
            let (mut e, mut o) = (0.0f32, 0.0f32);
            for q in 0..m4 {
                if q % 2 == 0 {
                    e = a[q].mul_add(b[q], e);
                } else {
                    o = a[q].mul_add(b[q], o);
                }
            }
            let mut s = e + o;
            for q in m4..k {
                s = a[q].mul_add(b[q], s);
            }
            s
        }
        Mode::ParU => {
            let m4 = k - k % 2;
            let (mut e, mut o) = (0.0f32, 0.0f32);
            for q in 0..m4 {
                let p = rd(a[q] * b[q]);
                if q % 2 == 0 {
                    e = rd(e + p);
                } else {
                    o = rd(o + p);
                }
            }
            let mut s = rd(e + o);
            for q in m4..k {
                s = rd(s + rd(a[q] * b[q]));
            }
            s
        }
    }
}

/// Column blocks of the nocopy driver (`0xa7be50`): remaining >= 192: 96; > 144: ceil(rem / 2); else rem.
pub fn n_blocks(n: usize) -> Vec<usize> {
    let mut out = Vec::new();
    let mut rem = n;
    while rem > 0 {
        let b = if rem >= 192 {
            96
        } else if rem > 144 {
            (rem + 1) / 2
        } else {
            rem
        };
        out.push(b);
        rem -= b;
    }
    out
}

/// Row tile heights of the copy-based kernel: panels of 24, then the bits 16, 8, 4, 2, 1 of `m % 24`.
fn tiles24(m: usize) -> Vec<usize> {
    let mut t = vec![24; m / 24];
    for g in [16usize, 8, 4, 2, 1] {
        if (m % 24) & g != 0 {
            t.push(g);
        }
    }
    t
}

#[derive(Clone, Copy, PartialEq, Eq)]
enum ColK {
    G4,
    C2,
    C1,
}

fn col_kernel(n: usize, j: usize) -> ColK {
    let (q, r) = (n / 4, n % 4);
    if j < 4 * q {
        return ColK::G4;
    }
    match r {
        1 => ColK::C1,
        2 => ColK::C2,
        _ => {
            if j < 4 * q + 2 {
                ColK::C2
            } else {
                ColK::C1
            }
        }
    }
}

fn copy_mode(tile: usize, ck: ColK) -> Mode {
    match (tile, ck) {
        (24, _) | (16, _) => Mode::Chain,
        (1, ColK::C1) => Mode::Par4,
        (1, _) => Mode::Chain,
        (_, ColK::G4) | (_, ColK::C1) => Mode::Par4,
        (_, ColK::C2) => Mode::Chain,
    }
}

/// `true`: copy-based kernel, `false`: nocopy driver.
pub fn copy_path(m: usize, n: usize) -> bool {
    m >= 56 && n >= 56
}

/// Mode of element (i, j) of the nocopy path; `cg` is the width of the column group (6 for a full group of 6).
fn nocopy_mode(m: usize, i: usize, cg: usize) -> Mode {
    let r = m % 16;
    let big = i < m - r || r >= 9;
    if cg >= 4 {
        if big {
            Mode::Chain
        } else {
            Mode::Par
        }
    } else if cg == 3 {
        Mode::Par
    } else {
        Mode::ParU
    }
}

/// Accumulation mode of element (i, j) (also for tests), `m, n >= 2`.
pub fn element_mode(m: usize, n: usize, i: usize, j: usize) -> Mode {
    if copy_path(m, n) {
        let mut r0 = 0;
        for t in tiles24(m) {
            if i < r0 + t {
                return copy_mode(t, col_kernel(n, j));
            }
            r0 += t;
        }
        unreachable!()
    }
    let mut j0 = 0;
    for nb in n_blocks(n) {
        if j < j0 + nb {
            let jj = j - j0;
            let r6 = nb % 6;
            let cg = if jj < nb - r6 { 6 } else { r6 };
            return nocopy_mode(m, i, cg);
        }
        j0 += nb;
    }
    unreachable!()
}

/// One dot product of the SGEMV('T') kernel (`mkl_blas_avx2_xsgemv_t`), the float analogue of `dgemv_t_col`: 8 lanes; blocks of 16 terms
/// `t = fma(x_c, a_c, x_{c+8} a_{c+8})`, `lane += t`; then a block of 8 (`fused`: fma, else rounded product and sum) and the remaining `k % 8` terms fma'd
/// into lanes `0..`; the lanes are reduced as `((l0 + l1) + (l4 + l5)) + ((l2 + l3) + (l6 + l7))`.
fn gemv_dot(a: &[f32], x: &[f32], fused: bool) -> f32 {
    let k = a.len();
    let mut l = [0.0f32; 8];
    let main = k - k % 16;
    let tail8 = (k % 16) & !7;
    let r = k % 8;
    let mut c = 0;
    while c < main {
        for q in 0..8 {
            let t = x[c + q].mul_add(a[c + q], rd(x[c + 8 + q] * a[c + 8 + q]));
            l[q] = rd(t + l[q]);
        }
        c += 16;
    }
    let mut p = main;
    if tail8 != 0 {
        for q in 0..8 {
            l[q] = if fused {
                x[p + q].mul_add(a[p + q], l[q])
            } else {
                rd(rd(x[p + q] * a[p + q]) + l[q])
            };
        }
        p += 8;
    }
    for q in 0..r {
        l[q] = x[p + q].mul_add(a[p + q], l[q]);
    }
    rd(rd(rd(l[0] + l[1]) + rd(l[4] + l[5])) + rd(rd(l[2] + l[3]) + rd(l[6] + l[7])))
}

/// `C += A^T B`, in place (`ldc` row stride of the column-major `C`).
pub fn sgemm_t_n(
    m: usize,
    n: usize,
    k: usize,
    a: &[f32],
    lda: usize,
    b: &[f32],
    ldb: usize,
    c: &mut [f32],
    ldc: usize,
) {
    for j in 0..n {
        for i in 0..m {
            let s = if m == 1 || n == 1 {
                let var = if m == 1 {
                    if j < n - n % 4 {
                        4
                    } else {
                        1
                    }
                } else if i < m - m % 4 {
                    4
                } else {
                    1
                };
                gemv_dot(&a[i * lda..i * lda + k], &b[j * ldb..j * ldb + k], var == 4)
            } else {
                // The nocopy driver (`mkl_blas_avx2_xsgemm_nocopy_driver`) calls the kernel per piece of `k`: blocks of 256 while more than 512 remain, then
                // `ceil(r / 2)` and the rest when `r > 256` (traced for k = 257 .. 1100); each piece is a `beta = 1` update of `C`.
                let mode = element_mode(m, n, i, j);
                if !copy_path(m, n) && k > 256 {
                    let (mut k0, mut r) = (0, k);
                    let mut ps = Vec::new();
                    while r > 512 {
                        ps.push(256);
                        r -= 256;
                    }
                    if r > 256 {
                        ps.push((r + 1) / 2);
                        r -= (r + 1) / 2;
                    }
                    ps.push(r);
                    let mut cv = c[j * ldc + i];
                    for kk in ps {
                        let t = acc(
                            mode,
                            kk,
                            &a[i * lda + k0..i * lda + k0 + kk],
                            &b[j * ldb + k0..j * ldb + k0 + kk],
                        );
                        cv = rd(cv + t);
                        k0 += kk;
                    }
                    c[j * ldc + i] = cv;
                    continue;
                }
                acc(mode, k, &a[i * lda..i * lda + k], &b[j * ldb..j * ldb + k])
            };
            c[j * ldc + i] = rd(c[j * ldc + i] + s);
        }
    }
}

// ---------------------------------------------------------------------------------------------------------------------------------
// STRSM('L', 'U', 'T', 'U')

/// STRSM('L','U','T','U', alpha = 1) with `n` right-hand sides (`0xa00350`); `l[i * w + j]` (row-major, `j < i`) is the unit lower triangle of the
/// `w x w` block (the transposed upper triangle of the call), `x` the column `col` of `n` right-hand sides, solved in place.
/// `w >= 25`: `strsm_kernel_ll_0`: rows 0..23 an fnma chain from `b`; then tiles of 24 rows and the bits 16, 8, 4, 2, 1 of the rest, each with the
/// update from all rows solved before it, `x_i = b_i + acc` (accumulator structure of [`copy_mode`] per tile and column kernel), then the chain inside
/// the tile. Smaller `w` is not modelled here.
pub fn strsm_llu_col(l: &[f32], w: usize, x: &mut [f32], col: usize, n: usize) {
    if n == 1 {
        strsv_utu_small(l, w, x);
        return;
    }
    if w < 25 {
        if w >= 8 && col < 8 * (n / 8) {
            strsm_a8_col(l, w, x);
        } else {
            strsm_left_ker(l, w, x);
        }
        return;
    }
    let ck = col_kernel(n, col);
    let blk = 24;
    for i in 0..blk.min(w) {
        let mut v = x[i];
        for j in 0..i {
            v = (-l[i * w + j]).mul_add(x[j], v);
        }
        x[i] = v;
    }
    let mut r0 = blk;
    let mut tiles = vec![24usize; (w - blk) / 24];
    for g in [16usize, 8, 4, 2, 1] {
        if ((w - blk) % 24) & g != 0 {
            tiles.push(g);
        }
    }
    for t in tiles {
        let mode = copy_mode(t, ck);
        for i in r0..r0 + t {
            let neg: Vec<f32> = (0..r0).map(|j| -l[i * w + j]).collect();
            let acc = acc(mode, r0, &neg, &x[..r0]);
            x[i] = rd(x[i] + acc);
        }
        for i in r0..r0 + t {
            let mut v = x[i];
            for j in r0..i {
                v = (-l[i * w + j]).mul_add(x[j], v);
            }
            x[i] = v;
        }
        r0 += t;
    }
}

/// `strsm_left_ker` (`m <= 24` with fewer than 8 right-hand sides, and the remainder columns of the 8-row kernel): `x_i = b_i - sum_{j<i} l_ij x_j` with
/// rounded products and sums; rows `i < 16` sequentially, rows `i >= 16` in 8 lanes over blocks of 8 terms (lane 0 starts at `b_i`),
/// `((l0 + l4) + (l2 + l6)) + ((l1 + l5) + (l3 + l7))`, then the remaining terms sequentially.
pub fn strsm_left_ker(l: &[f32], w: usize, x: &mut [f32]) {
    for i in 0..w {
        let mut v = x[i];
        let mut j0 = 0;
        if i >= 16 {
            let nch = i / 8;
            let mut lane = [0.0f32; 8];
            lane[0] = v;
            for c in 0..nch {
                for k in 0..8 {
                    let j = 8 * c + k;
                    // The first block of lanes is fused (`b - l x` in one rounding, visible when `l_i0 x_0` cancels `b_i`), the later blocks unfused.
                    lane[k] = if c == 0 {
                        (-l[i * w + j]).mul_add(x[j], lane[k])
                    } else {
                        rd(lane[k] - rd(l[i * w + j] * x[j]))
                    };
                }
            }
            v = rd(rd(rd(lane[0] + lane[4]) + rd(lane[2] + lane[6]))
                + rd(rd(lane[1] + lane[5]) + rd(lane[3] + lane[7])));
            j0 = 8 * nch;
        }
        for j in j0..i {
            v = rd(v - rd(l[i * w + j] * x[j]));
        }
        x[i] = v;
    }
}

/// `strsm_ker_llu_a8_b8` for the columns of full blocks of 8 right-hand sides (`8 <= w <= 24`; the other columns use [`strsm_left_ker`]): the first `8 (w / 8)` rows an fma chain from `b`, the
/// remaining rows `x_i = b_i - s` with `s` the 8-lane sum over the solved rows (see [`strsm_left_ker`] for the lanes), then the chain over the rows of the tail.
pub fn strsm_a8_col(l: &[f32], w: usize, x: &mut [f32]) {
    let r0 = 8 * (w / 8);
    for i in 0..r0 {
        let mut v = x[i];
        for j in 0..i {
            v = (-l[i * w + j]).mul_add(x[j], v);
        }
        x[i] = v;
    }
    for i in r0..w {
        let mut lane = [0.0f32; 8];
        for c in 0..r0 / 8 {
            for k in 0..8 {
                let j = 8 * c + k;
                lane[k] = rd(rd(l[i * w + j] * x[j]) + lane[k]);
            }
        }
        let s = rd(rd(rd(lane[0] + lane[4]) + rd(lane[2] + lane[6]))
            + rd(rd(lane[1] + lane[5]) + rd(lane[3] + lane[7])));
        let mut v = rd(x[i] + (-s));
        for j in r0..i {
            v = rd(v - rd(l[i * w + j] * x[j]));
        }
        x[i] = v;
    }
}

/// `strsv_utu` (`n == 1`): rows in pairs `(p, p + 1)`, `p` even (a last odd row alone). `reduce(p)` is `b_p` minus the products of the first `p` solved rows: `p < 8`
/// sequentially, else 8 lanes over blocks of 8 terms (lane 0 starts at `b`) reduced as in [`strsm_left_ker`], then the remaining terms sequentially; all with rounded
/// products and sums. A lone last row uses the threshold 16 of [`strsm_left_ker`] instead of 8. Row `p` takes `reduce(p)`; row `p + 1` takes `fma(-x_p, l_{p+1,p}, reduce_{p+1}(p))`. Verified for `w <= 24` only.
pub fn strsv_utu_small(l: &[f32], w: usize, x: &mut [f32]) {
    let reduce = |x: &[f32], i: usize, p: usize, thr: usize| -> f32 {
        let mut v = x[i];
        let mut j0 = 0;
        if p >= thr {
            let nch = p / 8;
            let mut lane = [0.0f32; 8];
            lane[0] = v;
            for c in 0..nch {
                for k in 0..8 {
                    let j = 8 * c + k;
                    lane[k] = rd(lane[k] - rd(l[i * w + j] * x[j]));
                }
            }
            v = rd(rd(rd(lane[0] + lane[4]) + rd(lane[2] + lane[6]))
                + rd(rd(lane[1] + lane[5]) + rd(lane[3] + lane[7])));
            j0 = 8 * nch;
        }
        for j in j0..p {
            v = rd(v - rd(l[i * w + j] * x[j]));
        }
        v
    };
    let mut p = 0;
    while p + 1 < w {
        let v0 = reduce(x, p, p, 8);
        let v1 = reduce(x, p + 1, p, 8);
        x[p] = v0;
        x[p + 1] = (-v0).mul_add(l[(p + 1) * w + p], v1);
        p += 2;
    }
    if p < w {
        x[p] = reduce(x, p, p, 16);
    }
}
