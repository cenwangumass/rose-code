//! Logical destinations for original grouped double GEMM. Public GEMM remains the oracle.
use super::*;
#[derive(Clone, Copy)]
pub(super) struct Index {
    pub row: usize,
    pub col: usize,
}
pub(super) trait Destination {
    // bounded() validates every physical location accessed within shape().
    fn shape(&self) -> (usize, usize);
    /// Caller supplies a row/column strictly inside shape(), after bounded().
    unsafe fn load(&self, index: Index) -> f64;
    /// Same contract as load; each cell is written only once by this kernel.
    unsafe fn store(&mut self, index: Index, value: f64);
    fn bounded(&self) -> bool;
}
// Zero is allowed; normal operands have quantum >=2^-502 and magnitude <2^451.
// Products/partial sums for k<=240 have quantum >=2^-1004 and magnitude <2^910.
// Gathered initial C values in [-950,950] exponent range preserve normal-or-zero
// intermediates and finite results in every rounding/flush mode. Zero-C scatter
// admits any parent: each final addition has at most one nonfinite operand.
// Out-of-range A/B or gathered C use the materialized path before any mutation.
pub(super) fn operand_bound(x: f64) -> bool {
    let bits = x.to_bits() & 0x7fff_ffff_ffff_ffff;
    bits == 0 || (573..=1473).contains(&((bits >> 52) as u16))
}
pub(super) fn parent_bound(x: f64) -> bool {
    let bits = x.to_bits() & 0x7fff_ffff_ffff_ffff;
    bits == 0 || (73..=1973).contains(&((bits >> 52) as u16))
}
pub(super) fn try_update<C: Destination>(
    m: usize,
    n: usize,
    k: usize,
    a: &[f64],
    lda: usize,
    b: &[f64],
    ldb: usize,
    c: &mut C,
    ldc: usize,
) -> bool {
    #[cfg(target_arch = "x86_64")]
    if !SINGLE.with(|s| s.get())
        && k > 0
        && k <= 240
        && (k <= 192 || gemm_path2(m, n, k))
        && std::arch::is_x86_feature_detected!("fma")
    {
        let mut controls = 0u32;
        unsafe {
            std::arch::asm!("stmxcsr [{p}]",p=in(reg)&mut controls,options(nostack,preserves_flags));
        }
        if controls & 0x1f80 == 0x1f80
            && a.iter().copied().all(operand_bound)
            && b.iter().copied().all(operand_bound)
            && m <= c.shape().0
            && n <= c.shape().1
            && c.bounded()
        {
            unsafe {
                native(m, n, k, a, lda, b, ldb, c, ldc);
            }
            return true;
        }
    }
    false
}
#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "fma")]
unsafe fn native<C: Destination>(
    m: usize,
    n: usize,
    k: usize,
    a: &[f64],
    lda: usize,
    b: &[f64],
    ldb: usize,
    c: &mut C,
    ldc: usize,
) {
    body(m, n, k, a, lda, b, ldb, c, ldc);
}
#[inline(always)]
unsafe fn body<C: Destination>(
    m: usize,
    n: usize,
    k: usize,
    a: &[f64],
    lda: usize,
    b: &[f64],
    ldb: usize,
    c: &mut C,
    _ldc: usize,
) {
    const GROUP_S: bool = true;
    if m == 1 {
        for j in 0..n {
            let variant = if j < n - n % 4 {
                4
            } else if n % 4 & 2 != 0 && j < n - n % 4 + 2 {
                2
            } else {
                1
            };
            let index = Index { row: 0, col: j };
            let value = dgemv_t_col(
                &b[j * ldb..j * ldb + k],
                &a[..k],
                variant,
                1.0,
                c.load(index),
            );
            c.store(index, value);
        }
        return;
    }
    if n == 1 {
        for i in 0..m {
            let variant = if i < m - m % 4 {
                4
            } else if m % 4 & 2 != 0 && i < m - m % 4 + 2 {
                2
            } else {
                1
            };
            let index = Index { row: i, col: 0 };
            let value = dgemv_t_col(
                &a[i * lda..i * lda + k],
                &b[..k],
                variant,
                1.0,
                c.load(index),
            );
            c.store(index, value);
        }
        return;
    }
    if gemm_path2(m, n, k) {
        for j in 0..n {
            let ck = col_kernel(n, j);
            for i in 0..m {
                let mode = tile_mode(k0_row_tile(m, i), ck);
                let (ar, br) = (&a[i * lda..i * lda + k], &b[j * ldb..j * ldb + k]);
                let s = acc_terms(mode, k, &|t| (1.0 * ar[t], br[t]));
                let index = Index { row: i, col: j };
                c.store(index, c.load(index) + s);
            }
        }
        return;
    }
    // The nocopy driver (`0xa7be50`) cuts `k > 192` into calls of the kernel, each a `beta = 1` update of `C` (the dispatch above sees the whole `k`;
    // the copy-based kernel takes `k` up to 240 in one call).
    // The pieces: while more than 384 remain a block of 192; then up to 384 `ceil(r / 2)` and the rest; up to 192 one call (traced for k = 385, 400, 512, 700).
    let mut ps = Vec::new();
    let mut r = k;
    while r > 384 {
        ps.push(192);
        r -= 192;
    }
    if r > 192 {
        ps.push((r + 1) / 2);
        r -= (r + 1) / 2;
    }
    if r > 0 {
        ps.push(r);
    }
    let mut k0 = 0;
    for &kk in &ps {
        let mut j0 = 0;
        for nb in nocopy_n_blocks(n) {
            let mut j = 0;
            while j < nb {
                #[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
                if GROUP_S
                    && j + 3 < nb
                    && gemm_kind(m, nb, 0, j) == Kind::S
                    && gemm_kind(m, nb, 0, j + 3) == Kind::S
                {
                    let ar = |row: usize| &a[row * lda + k0..row * lda + k0 + kk];
                    let br = |col: usize| &b[(j0 + col) * ldb + k0..(j0 + col) * ldb + k0 + kk];
                    let mut i = 0;
                    while i < m {
                        if i + 4 <= m
                            && gemm_kind(m, nb, i, j) == Kind::S
                            && gemm_kind(m, nb, i + 3, j) == Kind::S
                        {
                            // SAFETY: GROUP_S is enabled only by the runtime-detected FMA backend.
                            let sums = unsafe {
                                dgemm_s4x4(
                                    [ar(i), ar(i + 1), ar(i + 2), ar(i + 3)],
                                    [br(j), br(j + 1), br(j + 2), br(j + 3)],
                                )
                            };
                            for col in 0..4 {
                                for row in 0..4 {
                                    let ci = Index {
                                        row: i + row,
                                        col: j0 + j + col,
                                    };
                                    c.store(ci, c.load(ci) + 1.0 * sums[col][row]);
                                }
                            }
                            i += 4;
                        } else {
                            for col in j..j + 4 {
                                let ci = Index {
                                    row: i,
                                    col: j0 + col,
                                };
                                c.store(
                                    ci,
                                    dgemm_tn_elem(
                                        gemm_kind(m, nb, i, col),
                                        ar(i),
                                        br(col),
                                        c.load(ci),
                                    ),
                                );
                            }
                            i += 1;
                        }
                    }
                    j += 4;
                    continue;
                }
                if GROUP_S
                    && j + 1 < nb
                    && gemm_kind(m, nb, 0, j) == Kind::S
                    && gemm_kind(m, nb, 0, j + 1) == Kind::S
                {
                    let br = |col: usize| &b[(j0 + col) * ldb + k0..(j0 + col) * ldb + k0 + kk];
                    let ar = |row: usize| &a[row * lda + k0..row * lda + k0 + kk];
                    let mut i = 0;
                    while i < m {
                        if i + 4 <= m
                            && gemm_kind(m, nb, i, j) == Kind::S
                            && gemm_kind(m, nb, i + 3, j) == Kind::S
                        {
                            let sums = dgemm_s4x2(
                                ar(i),
                                ar(i + 1),
                                ar(i + 2),
                                ar(i + 3),
                                br(j),
                                br(j + 1),
                            );
                            for col in 0..2 {
                                for row in 0..4 {
                                    let ci = Index {
                                        row: i + row,
                                        col: j0 + j + col,
                                    };
                                    c.store(ci, c.load(ci) + 1.0 * sums[col][row]);
                                }
                            }
                            i += 4;
                        } else {
                            for col in j..j + 2 {
                                let ci = Index {
                                    row: i,
                                    col: j0 + col,
                                };
                                c.store(
                                    ci,
                                    dgemm_tn_elem(
                                        gemm_kind(m, nb, i, col),
                                        ar(i),
                                        br(col),
                                        c.load(ci),
                                    ),
                                );
                            }
                            i += 1;
                        }
                    }
                    j += 2;
                    continue;
                }
                if GROUP_S {
                    let br = &b[(j0 + j) * ldb + k0..(j0 + j) * ldb + k0 + kk];
                    let mut i = 0;
                    while i < m {
                        let kind = gemm_kind(m, nb, i, j);
                        if kind == Kind::S && i + 4 <= m && gemm_kind(m, nb, i + 3, j) == Kind::S {
                            let ar = |row: usize| &a[row * lda + k0..row * lda + k0 + kk];
                            let sums = dgemm_s4(ar(i), ar(i + 1), ar(i + 2), ar(i + 3), br);
                            for t in 0..4 {
                                let ci = Index {
                                    row: i + t,
                                    col: j0 + j,
                                };
                                c.store(ci, c.load(ci) + 1.0 * sums[t]);
                            }
                            i += 4;
                        } else if kind == Kind::S
                            && i + 2 <= m
                            && gemm_kind(m, nb, i + 1, j) == Kind::S
                        {
                            let ar = |row: usize| &a[row * lda + k0..row * lda + k0 + kk];
                            let sums = dgemm_s2(ar(i), ar(i + 1), br);
                            for t in 0..2 {
                                let ci = Index {
                                    row: i + t,
                                    col: j0 + j,
                                };
                                c.store(ci, c.load(ci) + 1.0 * sums[t]);
                            }
                            i += 2;
                        } else {
                            let ci = Index {
                                row: i,
                                col: j0 + j,
                            };
                            c.store(
                                ci,
                                dgemm_tn_elem(
                                    kind,
                                    &a[i * lda + k0..i * lda + k0 + kk],
                                    br,
                                    c.load(ci),
                                ),
                            );
                            i += 1;
                        }
                    }
                } else {
                    for i in 0..m {
                        let kind = gemm_kind(m, nb, i, j);
                        let ci = Index {
                            row: i,
                            col: j0 + j,
                        };
                        c.store(
                            ci,
                            dgemm_tn_elem(
                                kind,
                                &a[i * lda + k0..i * lda + k0 + kk],
                                &b[(j0 + j) * ldb + k0..(j0 + j) * ldb + k0 + kk],
                                c.load(ci),
                            ),
                        );
                    }
                }
                j += 1;
            }
            j0 += nb;
        }
        k0 += kk;
    }
}

#[cfg(all(test, target_arch = "x86_64"))]
mod tests {
    use super::*;
    struct Matrix {
        values: Vec<f64>,
        leading: usize,
    }
    impl Destination for Matrix {
        fn shape(&self) -> (usize, usize) {
            (self.leading, self.values.len() / self.leading)
        }
        unsafe fn load(&self, i: Index) -> f64 {
            self.values[i.col * self.leading + i.row]
        }
        unsafe fn store(&mut self, i: Index, x: f64) {
            self.values[i.col * self.leading + i.row] = x;
        }
        fn bounded(&self) -> bool {
            self.values.iter().copied().all(parent_bound)
        }
    }
    fn csr() -> u32 {
        let mut x = 0;
        unsafe {
            std::arch::asm!("stmxcsr [{p}]", p=in(reg)&mut x, options(nostack,preserves_flags));
        }
        x
    }
    fn set_csr(x: u32) {
        unsafe {
            std::arch::asm!("ldmxcsr [{p}]", p=in(reg)&x, options(nostack,preserves_flags));
        }
    }
    struct Restore(u32);
    impl Drop for Restore {
        fn drop(&mut self) {
            set_csr(self.0);
        }
    }
    #[test]
    fn range_edges_cancellation_and_unmasked_decline() {
        let _restore = Restore(csr());
        let _single = SingleGuard::new(false);
        if !std::arch::is_x86_feature_detected!("fma") {
            return;
        }
        let lo = f64::from_bits((573u64 << 52) | 0x000f_ffff_ffff_ffff);
        let hi = f64::from_bits((1473u64 << 52) | 0x000a_51fc_3461_8673);
        let a = [lo, -lo, hi, -hi, lo, lo, -lo, -lo];
        let b = [hi, hi, lo, lo, lo, -lo, lo, -lo];
        for initial in [
            0.0,
            -0.0,
            f64::from_bits(73u64 << 52),
            f64::from_bits(1973u64 << 52),
        ] {
            for mode in 0..16 {
                let controls = 0x1f80
                    | ((mode & 3) << 13)
                    | (((mode >> 2) & 1) << 15)
                    | (((mode >> 3) & 1) << 6);
                let mut old = vec![initial; 4];
                set_csr(controls);
                std::hint::black_box(dgemm_t_n)(2, 2, 4, &a, 4, &b, 4, &mut old, 2);
                let flags = csr();
                let mut dst = Matrix {
                    values: vec![initial; 4],
                    leading: 2,
                };
                set_csr(controls);
                assert!(try_update(2, 2, 4, &a, 4, &b, 4, &mut dst, 2));
                assert_eq!(flags, csr());
                assert_eq!(
                    old.iter().map(|x| x.to_bits()).collect::<Vec<_>>(),
                    dst.values.iter().map(|x| x.to_bits()).collect::<Vec<_>>()
                );
            }
        }
        let mut dst = Matrix {
            values: vec![1.0; 4],
            leading: 2,
        };
        set_csr(0x1f00);
        assert!(!try_update(2, 2, 4, &a, 4, &b, 4, &mut dst, 2));
        assert_eq!(dst.values, vec![1.0; 4]);
        assert_eq!(csr(), 0x1f00);
    }

    #[test]
    fn logical_destination_matches_public_gemm() {
        let _restore = Restore(csr());
        let _single = SingleGuard::new(false);
        let mut seed = 0x68da1b734f36u64;
        let mut value = || {
            seed ^= seed << 13;
            seed ^= seed >> 7;
            seed ^= seed << 17;
            f64::from_bits((seed & 0x800f_ffff_ffff_ffff) | ((1019 + (seed >> 52) % 9) << 52))
        };
        let mut count = 0;
        for (m, n) in [
            (1, 1),
            (1, 7),
            (7, 1),
            (2, 3),
            (4, 5),
            (11, 9),
            (25, 25),
            (101, 4),
        ] {
            for k in [1, 2, 3, 4, 7, 8, 9, 191, 192, 193, 240, 385] {
                let a: Vec<_> = (0..m * k).map(|_| value()).collect();
                let b: Vec<_> = (0..n * k).map(|_| value()).collect();
                let c: Vec<_> = (0..(m + 3) * n).map(|_| value()).collect();
                for mode in 0..16 {
                    let controls = 0x1f80
                        | ((mode & 3) << 13)
                        | (((mode >> 2) & 1) << 15)
                        | (((mode >> 3) & 1) << 6);
                    let mut expected = c.clone();
                    set_csr(controls);
                    std::hint::black_box(dgemm_t_n)(m, n, k, &a, k, &b, k, &mut expected, m + 3);
                    let flags = csr();
                    let mut actual = Matrix {
                        values: c.clone(),
                        leading: m + 3,
                    };
                    set_csr(controls);
                    let used = try_update(m, n, k, &a, k, &b, k, &mut actual, m + 3);
                    if used {
                        assert_eq!(
                            expected.iter().map(|x| x.to_bits()).collect::<Vec<_>>(),
                            actual
                                .values
                                .iter()
                                .map(|x| x.to_bits())
                                .collect::<Vec<_>>(),
                            "m={m} n={n} k={k} mode={mode}"
                        );
                        assert_eq!(flags, csr(), "flags m={m} n={n} k={k} mode={mode}");
                        count += 1;
                    } else {
                        assert_eq!(actual.values, c);
                        assert!(
                            k > 240
                                || (k > 192 && !gemm_path2(m, n, k))
                                || !std::arch::is_x86_feature_detected!("fma")
                        );
                    }
                }
            }
        }
        assert!(count > 1000);
        for exceptional in [f64::NAN, f64::INFINITY, f64::from_bits(1), 1e300, 1e-300] {
            let mut dst = Matrix {
                values: vec![1.0],
                leading: 1,
            };
            assert!(!try_update(
                1,
                1,
                1,
                &[exceptional],
                1,
                &[1.0],
                1,
                &mut dst,
                1
            ));
            assert_eq!(dst.values, vec![1.0]);
        }
    }
}
