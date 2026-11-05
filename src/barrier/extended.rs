//! Exact 64-bit extended arithmetic for the barrier Kahan recurrence. Each
//! product, correction, sum and subtraction retains the reference order and
//! returns the same extended sign, exponent and significand as software Ext.
use super::Ext;

pub(super) fn default_environment() -> bool {
    #[cfg(target_arch = "x86_64")]
    {
        let mut csr = 0u32;
        let mut cw = 0u16;
        unsafe {
            std::arch::asm!("stmxcsr [{p}]",p=in(reg) &mut csr,options(nostack,preserves_flags));
            std::arch::asm!("fnstcw [{p}]",p=in(reg) &mut cw,options(nostack,preserves_flags));
        }
        csr & 0xffc0 == 0x1f80 && cw & 0x0f3f == 0x033f
    }
    #[cfg(not(target_arch = "x86_64"))]
    {
        false
    }
}

// A 128-bit product or compact extended sum with one final nearest-even rounding.
#[inline(always)]
fn finish_significand(neg: bool, v: u128, e: i32, shift: u32) -> Ext {
    let head = (v >> shift) as u64;
    let rem = (v - ((head as u128) << shift)) as u64;
    let half = 1u64 << (shift - 1);
    let (q, carry) = head.overflowing_add(u64::from(rem > half || rem == half && head & 1 != 0));
    if carry {
        Ext {
            neg,
            m: 1u64 << 63,
            e: e + shift as i32 + 1,
        }
    } else {
        Ext {
            neg,
            m: q,
            e: e + shift as i32,
        }
    }
}
#[inline(always)]
fn checked_product(a: f64, b: f64) -> Option<Ext> {
    let (x, y) = (a.to_bits(), b.to_bits());
    let (ex, ey) = (((x >> 52) & 2047) as i32, ((y >> 52) & 2047) as i32);
    let nx = (1..2047).contains(&ex);
    let ny = (1..2047).contains(&ey);
    if !nx || !ny {
        if !nx && x << 1 != 0 || !ny && y << 1 != 0 {
            return None;
        }
        return Some(Ext::ZERO);
    }
    let m = |v: u64| (v & ((1u64 << 52) - 1)) | (1 << 52);
    let raw = u128::from(m(x)) * u128::from(m(y));
    let shift = 41 + u32::from(raw >= 1u128 << 105);
    Some(finish_significand(
        (x ^ y) >> 63 != 0,
        raw,
        ex + ey - 2150,
        shift,
    ))
}
#[cfg(test)]
fn exact_product(a: f64, b: f64) -> Ext {
    checked_product(a, b).unwrap_or_else(|| Ext::mul_f64(a, b))
}

#[inline(always)]
fn normalize_difference(neg: bool, v: u64, e: i32) -> Ext {
    if v == 0 {
        return Ext::ZERO;
    }
    let lz = v.leading_zeros();
    Ext {
        neg,
        m: v << lz,
        e: e - lz as i32,
    }
}
#[inline(always)]
// Align three low guard bits, using the same low-bit jam as Ext::add. The
// legacy alignment uses 62 bits. Shifting its result right 59 bits with a jam
// gives this compact result: the large operand has 3 trailing zeros, so the
// low-bit jam commutes with both adding and subtracting the smaller operand.
// For a gap above three, at most one leading bit cancels; at least two guard
// bits remain. Rounding the compact result to nearest-even then gives exactly
// the original significand. Gaps up to three are aligned without a lost bit.
fn exact_add(x: Ext, y: Ext) -> Ext {
    if x.m == 0 {
        return y;
    }
    if y.m == 0 {
        return x;
    }
    let (a, b) = if (x.e, x.m) >= (y.e, y.m) {
        (x, y)
    } else {
        (y, x)
    };
    let gap = (a.e - b.e) as u32;
    // Below a quarter ulp of a, including the smaller predecessor spacing
    // when a is a power of two. The original jammed alignment is also below it.
    if gap >= 66 {
        return a;
    }
    let same = a.neg == b.neg;
    if gap == 0 {
        if !same {
            return normalize_difference(a.neg, a.m - b.m, a.e);
        }
        return finish_significand(a.neg, u128::from(a.m) + u128::from(b.m), a.e, 1);
    }
    if gap == 1 && !same {
        let v = (u128::from(a.m) << 1) - u128::from(b.m);
        if v <= u128::from(u64::MAX) {
            return normalize_difference(a.neg, v as u64, b.e);
        }
        return finish_significand(a.neg, v, b.e, 1);
    }
    let large = u128::from(a.m) << 3;
    let small = u128::from(b.m) << 3;
    let shifted = (small >> gap) | u128::from(small & ((1u128 << gap) - 1) != 0);
    let v = if same {
        large + shifted
    } else {
        large - shifted
    };
    let shift = if v >= 1u128 << 67 {
        4
    } else if v >= 1u128 << 66 {
        3
    } else {
        2
    };
    finish_significand(a.neg, v, a.e - 3, shift)
}
#[inline(always)]
fn exact_sub(x: Ext, y: Ext) -> Ext {
    exact_add(x, y.neg())
}
/// Computes both original feedback operations. With a common sign and
/// exponent the difference is an exact signed integer at that exponent.
#[inline(always)]
fn feedback(next: Ext, sum: Ext, y: Ext) -> (Ext, Ext) {
    if next.m != 0 && sum.m != 0 && next.e == sum.e && next.neg == sum.neg {
        let (m, neg) = if next.m >= sum.m {
            (next.m - sum.m, next.neg)
        } else {
            (sum.m - next.m, !next.neg)
        };
        let h = normalize_difference(neg, m, sum.e);
        if m == 0 || y.m == 0 {
            return (h, h.sub(y));
        }
        let gap = sum.e - y.e;
        if (0..64).contains(&gap) {
            let lhs = u128::from(m) << gap;
            let rhs = u128::from(y.m);
            let (v, sign) = if neg != y.neg {
                (lhs + rhs, neg)
            } else if lhs >= rhs {
                (lhs - rhs, neg)
            } else {
                (rhs - lhs, !neg)
            };
            if v <= u128::from(u64::MAX) {
                // h-y fits the original extended precision on this common
                // grid. No tail is present to round away.
                return (h, normalize_difference(sign, v as u64, y.e));
            }
        }
        return (h, exact_sub(h, y));
    }
    let h = exact_sub(next, sum);
    (h, exact_sub(h, y))
}
#[inline(never)]
pub(super) fn stream(a: &[f64], b: &[f64]) -> Ext {
    let (mut sum, mut c) = (Ext::ZERO, Ext::ZERO);
    for i in 0..a.len() {
        let t = match checked_product(a[i], b[i]) {
            Some(v) => v,
            None => return super::software_kahan_tail(&a[i..], &b[i..], sum, c),
        };
        let y = exact_sub(t, c);
        let next = exact_add(sum, y);
        let (_, new_c) = feedback(next, sum, y);
        c = new_c;
        sum = next;
    }
    sum
}

#[cfg(test)]
mod tests {
    use super::*;
    fn same(x: Ext, y: Ext) {
        assert_eq!((x.neg, x.m, x.e), (y.neg, y.m, y.e));
    }
    fn random(s: &mut u64) -> u64 {
        *s ^= *s << 13;
        *s ^= *s >> 7;
        *s ^= *s << 17;
        *s
    }
    fn shadow(a: &[f64], b: &[f64]) {
        let (mut s, mut c, mut sr, mut cr) = (Ext::ZERO, Ext::ZERO, Ext::ZERO, Ext::ZERO);
        for (&a, &b) in a.iter().zip(b) {
            let t = exact_product(a, b);
            let tr = Ext::mul_f64(a, b);
            same(t, tr);
            let y = exact_sub(t, c);
            let yr = tr.sub(cr);
            same(y, yr);
            let n = exact_add(s, y);
            let nr = sr.add(yr);
            same(n, nr);
            let (h, next_c) = feedback(n, s, y);
            let hr = nr.sub(sr);
            same(h, hr);
            c = next_c;
            cr = hr.sub(yr);
            same(c, cr);
            s = n;
            sr = nr;
        }
        same(
            stream(a, b),
            super::super::software_kahan_tail(a, b, Ext::ZERO, Ext::ZERO),
        );
    }
    #[test]
    fn extended_add_matches_full_alignment() {
        let mant = [
            0,
            1 << 63,
            (1 << 63) + 1,
            (1 << 63) + 2,
            (1 << 63) + 3,
            u64::MAX - 2,
            u64::MAX - 1,
            u64::MAX,
        ];
        for gap in (0..=140).chain([300, 2000]) {
            for &m in &mant {
                for &n in &mant {
                    for sign in 0..4 {
                        let a = Ext {
                            neg: sign & 1 != 0,
                            m,
                            e: 17,
                        };
                        let b = Ext {
                            neg: sign & 2 != 0,
                            m: n,
                            e: 17 - gap,
                        };
                        same(exact_add(a, b), a.add(b));
                        same(exact_sub(a, b), a.sub(b));
                        same(exact_add(b, a), b.add(a));
                    }
                }
            }
        }
        let mut seed = 0x5a179f33472;
        for _ in 0..32_768 {
            let p = random(&mut seed);
            let q = random(&mut seed);
            let a = Ext {
                neg: p & 1 != 0,
                m: p | 1 << 63,
                e: (p % 143) as i32,
            };
            let b = Ext {
                neg: q & 2 != 0,
                m: q | 1 << 63,
                e: (q % 143) as i32,
            };
            same(exact_add(a, b), a.add(b));
            same(exact_sub(a, b), a.sub(b));
        }
    }
    #[test]
    fn all_five_kahan_values_match_software() {
        let special = [
            0.0,
            -0.0,
            1.0,
            -1.0,
            f64::MIN_POSITIVE,
            f64::from_bits(1),
            f64::MAX,
            f64::INFINITY,
            f64::NEG_INFINITY,
            f64::from_bits(0x7ff0_0000_0000_0012),
            f64::from_bits(0xfff8_0000_0000_ff45),
        ];
        let (mut a, mut b) = (Vec::new(), Vec::new());
        for &x in &special {
            for &y in &special {
                a.push(x);
                b.push(y);
            }
        }
        shadow(&a, &b);
        a.clear();
        b.clear();
        let mut seed = 0x1ce7f09a4117;
        for _ in 0..32_768 {
            a.push(f64::from_bits(random(&mut seed)));
            b.push(f64::from_bits(random(&mut seed)));
        }
        shadow(&a, &b);
        for (i, x) in a.iter_mut().enumerate() {
            *x = (1. + (i as f64 % 23.)) / 7.;
        }
        for (i, x) in b.iter_mut().enumerate() {
            *x = if i % 3 == 0 {
                -0.0
            } else {
                (i as f64 % 31.) / 2.
            };
        }
        shadow(&a, &b);
    }
}
