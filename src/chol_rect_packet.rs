//! Complete rectangular substitution with one independent RHS per SIMD lane.
//! LDL and reciprocal scaling remain in the original factor_block.
use super::*;

pub(super) fn try_solve(l: &[f64], w: usize, r: &mut [f64]) -> bool {
    #[cfg(target_arch = "x86_64")]
    if w >= 2
        && r.len() / w >= 4
        && std::arch::is_x86_feature_detected!("avx2")
        && std::arch::is_x86_feature_detected!("fma")
    {
        let mut controls = 0u32;
        unsafe {
            std::arch::asm!("stmxcsr [{p}]", p = in(reg) &mut controls, options(nostack, preserves_flags));
        }
        // Independent lanes can reorder exceptions between RHSs. Unmasked
        // environments retain the exact original ordered column execution.
        if controls & 0x1f80 == 0x1f80 {
            unsafe {
                native::solve(l, w, r);
            }
            return true;
        }
    }
    false
}

#[cfg(target_arch = "x86_64")]
mod native {
    use super::*;
    use std::arch::x86_64::*;

    #[target_feature(enable = "avx2,fma")]
    pub(super) unsafe fn solve(l: &[f64], w: usize, r: &mut [f64]) {
        let n = r.len() / w;
        let n4 = n / 4 * 4;
        let mut x = vec![_mm256_setzero_pd(); w];
        let ll0 = w >= TRSM_LL0_MIN && !(w % 4 == 0 && w <= 88);
        let mut tiles = Vec::new();
        if ll0 {
            let remaining = w - TRSM_LL0_BLOCK;
            tiles.resize(remaining / 12, 12usize);
            for t in [8usize, 4, 2, 1] {
                if remaining % 12 & t != 0 {
                    tiles.push(t);
                }
            }
        }
        let zero = _mm256_setzero_pd();
        let sign = _mm256_set1_pd(-0.0);
        let exponent = _mm256_set1_epi64x(0x7ff0_0000_0000_0000);
        for first in (0..n4).step_by(4) {
            for i in 0..w {
                x[i] = _mm256_set_pd(
                    r[(first + 3) * w + i],
                    r[(first + 2) * w + i],
                    r[(first + 1) * w + i],
                    r[first * w + i],
                );
            }
            if ll0 {
                for i in 0..TRSM_LL0_BLOCK {
                    let mut v = x[i];
                    for j in 0..i {
                        v = _mm256_fmadd_pd(_mm256_set1_pd(-l[i * w + j]), x[j], v);
                    }
                    x[i] = v;
                }
                let mut r0 = TRSM_LL0_BLOCK;
                for &t in &tiles {
                    // Every packed column belongs to the original full G4 class.
                    // Its tile-4 and tile-1 updates use even/odd accumulators;
                    // tile-12/8/2 use one FMA chain. No C1 padded lanes execute.
                    for i in r0..r0 + t {
                        let acc = if t == 4 || t == 1 {
                            let m4 = r0 - r0 % 4;
                            let (mut even, mut odd) = (zero, zero);
                            for j in (0..m4).step_by(2) {
                                even = _mm256_fmadd_pd(_mm256_set1_pd(-l[i * w + j]), x[j], even);
                                odd = _mm256_fmadd_pd(
                                    _mm256_set1_pd(-l[i * w + j + 1]),
                                    x[j + 1],
                                    odd,
                                );
                            }
                            let mut a = _mm256_add_pd(even, odd);
                            for j in m4..r0 {
                                a = _mm256_fmadd_pd(_mm256_set1_pd(-l[i * w + j]), x[j], a);
                            }
                            a
                        } else {
                            let mut a = zero;
                            for j in 0..r0 {
                                a = _mm256_fmadd_pd(_mm256_set1_pd(-l[i * w + j]), x[j], a);
                            }
                            a
                        };
                        x[i] = _mm256_add_pd(x[i], acc);
                    }
                    // Original LL0 boundary: external updates of the WHOLE tile
                    // precede its within-tile substitution.
                    for i in r0..r0 + t {
                        let mut v = x[i];
                        for j in r0..i {
                            v = _mm256_fmadd_pd(_mm256_set1_pd(-l[i * w + j]), x[j], v);
                        }
                        x[i] = v;
                    }
                    r0 += t;
                }
            } else {
                let r0 = w / 4 * 4;
                for i in 0..r0 {
                    let mut v = x[i];
                    for j in 0..i {
                        v = _mm256_fmadd_pd(_mm256_set1_pd(-l[i * w + j]), x[j], v);
                    }
                    x[i] = v;
                }
                for i in r0..w {
                    if r0 != 0 {
                        let mut acc = [zero; 4];
                        for block in 0..r0 / 4 {
                            for k in 0..4 {
                                let j = block * 4 + k;
                                let product = _mm256_mul_pd(_mm256_set1_pd(l[i * w + j]), x[j]);
                                acc[k] = _mm256_add_pd(product, acc[k]);
                            }
                        }
                        let s = _mm256_add_pd(
                            _mm256_add_pd(acc[0], acc[2]),
                            _mm256_add_pd(acc[1], acc[3]),
                        );
                        x[i] = _mm256_add_pd(x[i], _mm256_xor_pd(s, sign));
                    }
                }
                for i in r0..w {
                    let mut v = x[i];
                    for j in r0..i {
                        let product = _mm256_mul_pd(_mm256_set1_pd(l[i * w + j]), x[j]);
                        v = _mm256_sub_pd(v, product);
                    }
                    x[i] = v;
                }
            }
            let mut bad = 0;
            for &v in &x {
                let fields = _mm256_and_si256(_mm256_castpd_si256(v), exponent);
                let nonfinite = _mm256_cmpeq_epi64(fields, exponent);
                bad |= _mm256_movemask_pd(_mm256_castsi256_pd(nonfinite));
            }
            // __m256d holds exactly four contiguous f64 values. All values were
            // initialized; this immutable view cannot outlive or mutate x.
            let packed = std::slice::from_raw_parts(x.as_ptr().cast::<f64>(), 4 * w);
            for lane in 0..4 {
                let col = first + lane;
                if bad & (1 << lane) != 0 {
                    // Inputs in r have not been overwritten. Retain the original
                    // complete native attempt + nonfinite fallback for this RHS.
                    let old = dtrsm_llu_col(l, w, &r[col * w..(col + 1) * w], col, n);
                    r[col * w..(col + 1) * w].copy_from_slice(&old);
                } else {
                    for i in 0..w {
                        r[col * w + i] = packed[i * 4 + lane];
                    }
                }
            }
        }
        for col in n4..n {
            let old = dtrsm_llu_col(l, w, &r[col * w..(col + 1) * w], col, n);
            r[col * w..(col + 1) * w].copy_from_slice(&old);
        }
    }
}

#[cfg(all(test, target_arch = "x86_64"))]
mod tests {
    use super::*;
    use std::hint::black_box;

    fn csr() -> u32 {
        let mut x = 0;
        unsafe {
            std::arch::asm!("stmxcsr [{p}]", p = in(reg) &mut x, options(nostack,preserves_flags));
        }
        x
    }
    fn set_csr(x: u32) {
        unsafe {
            std::arch::asm!("ldmxcsr [{p}]", p = in(reg) &x, options(nostack,preserves_flags));
        }
    }
    struct Restore(u32);
    impl Drop for Restore {
        fn drop(&mut self) {
            set_csr(self.0);
        }
    }
    fn bits(v: &[f64]) -> Vec<u64> {
        v.iter().map(|v| v.to_bits()).collect()
    }

    #[inline(never)]
    fn original_rect(l: &[f64], w: usize, r: &mut [f64]) {
        let n = r.len() / w;
        for col in 0..n {
            let x = dtrsm_llu_col(l, w, &r[col * w..(col + 1) * w], col, n);
            r[col * w..(col + 1) * w].copy_from_slice(&x);
        }
    }
    #[inline(never)]
    fn packet_rect(l: &[f64], w: usize, r: &mut [f64]) {
        if !try_solve(l, w, r) {
            original_rect(l, w, r);
        }
    }
    fn compare_rect(l: &[f64], w: usize, r: &[f64]) {
        let before = csr() & !63;
        let (mut a, mut b) = (r.to_vec(), r.to_vec());
        set_csr(before);
        black_box(original_rect)(l, w, &mut a);
        let old = csr();
        set_csr(before);
        black_box(packet_rect)(l, w, &mut b);
        let new = csr();
        assert_eq!(bits(&a), bits(&b), "rectangle w={w} n={}", r.len() / w);
        assert_eq!(old, new, "rectangle flags w={w} n={}", r.len() / w);
    }

    #[test]
    fn rectangular_packets_match_original_column_programs() {
        let _restore = Restore(csr());
        let mut seed = 0x3729a17b_u64;
        let mut next = || {
            seed ^= seed << 13;
            seed ^= seed >> 7;
            seed ^= seed << 17;
            seed
        };
        if std::arch::is_x86_feature_detected!("avx2") && std::arch::is_x86_feature_detected!("fma")
        {
            set_csr(0x1f80);
            assert!(try_solve(&[0.125; 4], 2, &mut [0.25; 8]));
        }
        let mut checks = 0;
        for w in [
            1, 2, 3, 4, 7, 8, 12, 13, 14, 15, 16, 20, 24, 25, 28, 32, 40, 64, 88, 89, 96, 101, 192,
        ] {
            for n in [1, 2, 3, 4, 5, 7, 8, 9, 12, 17] {
                for pass in 0..4 {
                    let l: Vec<_> = (0..w * w)
                        .map(|i| match pass {
                            0 | 2 => (next() % 32) as f64 * 0.015625 - 0.25,
                            1 => f64::from_bits(next()),
                            _ => {
                                if i % 3 == 0 {
                                    1e200
                                } else {
                                    -1e200
                                }
                            }
                        })
                        .collect();
                    let r: Vec<_> = (0..w * n)
                        .map(|i| match pass {
                            0 => (next() % 128) as f64 * 0.125 - 7.0,
                            1 => f64::from_bits(next()),
                            2 => match i % 13 {
                                0 => -0.0,
                                1 => f64::from_bits(1),
                                2 => f64::INFINITY,
                                3 => f64::from_bits(0x7ff8_0000_0000_0123),
                                4 => f64::from_bits(0x7ff8_0000_0000_0456),
                                _ => 0.125,
                            },
                            _ => {
                                if i % 2 == 0 {
                                    1e200
                                } else {
                                    -1e200
                                }
                            }
                        })
                        .collect();
                    for mode in 0..if w <= 14 && n <= 9 { 16 } else { 1 } {
                        set_csr(
                            0x1f80
                                | ((mode & 3) << 13)
                                | (((mode >> 2) & 1) << 15)
                                | (((mode >> 3) & 1) << 6),
                        );
                        compare_rect(&l, w, &r);
                        checks += 1;
                    }
                }
            }
        }
        // These cases stay finite, so zero signs and denormal behavior are
        // checked on the packet path rather than hidden by exceptional replay.
        for w in [2, 3, 12, 13, 16, 25, 88, 89] {
            for pattern in 0..3 {
                let l: Vec<_> = (0..w * w)
                    .map(|i| match pattern {
                        0 => {
                            if i % 2 == 0 {
                                -0.0
                            } else {
                                0.0
                            }
                        }
                        _ => {
                            if i % 3 == 0 {
                                0.5
                            } else {
                                -0.5
                            }
                        }
                    })
                    .collect();
                let r: Vec<_> = (0..w * 4)
                    .map(|i| match pattern {
                        0 => {
                            if i % 2 == 0 {
                                -0.0
                            } else {
                                0.0
                            }
                        }
                        1 => {
                            if i % 2 == 0 {
                                f64::from_bits(1)
                            } else {
                                -f64::MIN_POSITIVE
                            }
                        }
                        _ => {
                            if i % 3 == 0 {
                                1e-300
                            } else {
                                -1e-310
                            }
                        }
                    })
                    .collect();
                for mode in 0..16 {
                    set_csr(
                        0x1f80
                            | ((mode & 3) << 13)
                            | (((mode >> 2) & 1) << 15)
                            | (((mode >> 3) & 1) << 6),
                    );
                    compare_rect(&l, w, &r);
                    checks += 1;
                }
            }
        }
        // Unmasked invalid/denormal inputs are safe finite values here. The
        // packet entry must decline before changing any data or executing math.
        let l = vec![0.125; 64];
        let mut r = vec![0.25; 32];
        let original = r.clone();
        for controls in [0x1f00, 0x1e80] {
            set_csr(controls);
            assert!(!try_solve(&l, 8, &mut r));
            assert_eq!(bits(&r), bits(&original));
            assert_eq!(csr(), controls);
        }
        println!("rectangular_packet_checks={checks}");
    }

    fn compare_block(
        mm: &[f64],
        r: &[f64],
        w: usize,
        c: usize,
        orig: &[f64],
        reg: Option<(usize, usize)>,
    ) {
        let before = csr() & !63;
        let (mut ma, mut mb) = (mm.to_vec(), mm.to_vec());
        let (mut ra, mut rb) = (r.to_vec(), r.to_vec());
        set_csr(before);
        let a = black_box(rect_oracle::factor_block)(&mut ma, &mut ra, w, c, orig, reg);
        let fa = csr();
        set_csr(before);
        let b = black_box(factor_block)(&mut mb, &mut rb, w, c, orig, reg);
        let fb = csr();
        assert_eq!(bits(&ma), bits(&mb), "factor block M w={w} n={}", c - w);
        assert_eq!(bits(&ra), bits(&rb), "factor block R w={w} n={}", c - w);
        assert_eq!(bits(&a.0), bits(&b.0), "factor diagonal");
        assert_eq!(bits(&a.1), bits(&b.1), "factor reciprocal");
        assert_eq!(fa, fb, "factor block flags w={w} n={}", c - w);
    }
    #[test]
    fn rectangular_packets_preserve_complete_factor_blocks() {
        let _restore = Restore(csr());
        let mut checks = 0;
        for w in [1, 2, 4, 8, 12, 13, 16, 24, 25, 32, 64, 88, 89] {
            for n in [0, 1, 3, 4, 7, 8, 9, 17] {
                let mut mm = vec![0.0; w * w];
                for i in 0..w {
                    for j in 0..=i {
                        mm[i * w + j] = if i == j {
                            32.0
                        } else {
                            ((i * 7 + j * 3) % 13) as f64 * 0.015625 - 0.125
                        };
                    }
                }
                let orig: Vec<_> = (0..w).map(|i| mm[i * w + i]).collect();
                let r: Vec<_> = (0..w * n).map(|i| (i % 31) as f64 * 0.125 - 2.0).collect();
                for single in [true, false] {
                    let _single = SingleGuard::new(single);
                    for mode in 0..if w <= 13 && n <= 8 { 16 } else { 1 } {
                        set_csr(
                            0x1f80
                                | ((mode & 3) << 13)
                                | (((mode >> 2) & 1) << 15)
                                | (((mode >> 3) & 1) << 6),
                        );
                        compare_block(
                            &mm,
                            &r,
                            w,
                            w + n,
                            &orig,
                            if n % 2 == 0 { Some((0, w)) } else { None },
                        );
                        checks += 1;
                    }
                }
            }
        }
        let records: serde_json::Value =
            serde_json::from_str(include_str!("../tests/fixtures/barrier_rect_blocks.json"))
                .unwrap();
        for rec in records.as_array().unwrap() {
            let values = |key: &str| {
                rec[key]
                    .as_array()
                    .unwrap()
                    .iter()
                    .map(|v| f64::from_bits(v.as_u64().unwrap()))
                    .collect::<Vec<_>>()
            };
            let (mm, r, orig) = (values("mm"), values("r"), values("orig"));
            let w = rec["w"].as_u64().unwrap() as usize;
            let c = rec["c"].as_u64().unwrap() as usize;
            let reg: Option<(usize, usize)> = serde_json::from_value(rec["reg"].clone()).unwrap();
            for single in [true, false] {
                let _single = SingleGuard::new(single);
                set_csr(0x1f80);
                compare_block(&mm, &r, w, c, &orig, reg);
                checks += 1;
            }
        }
        println!("complete_factor_block_checks={checks}");
    }
}
