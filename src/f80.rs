//! The x87 extended precision type ("long double", 64-bit significand) in software, round to nearest even.
//!
//! The reference's quad mode (`Quad = -1` in the release profile switches to it on a numerical failure of the dual) runs the whole
//! dual simplex on `long double` arrays with the x87 unit, so its values are reproduced bit for bit by this type: every operation
//! is rounded to 64 bits (the processor's default precision control on Linux), including gradual underflow.
//!
//! `repr` / `from_repr` give the 10 stored bytes as the explicit significand and the sign/exponent word.

use std::cmp::Ordering;
use std::ops::{Add, AddAssign, Div, DivAssign, Mul, MulAssign, Neg, Sub, SubAssign};

const E_INF: i32 = i32::MAX;
const E_NAN: i32 = i32::MAX - 1;
const E_MIN: i32 = -16382;
const E_MAX: i32 = 16383;

#[derive(Clone, Copy, Debug)]
pub struct F80 {
    neg: bool,
    /// Unbiased exponent: the value is `m * 2^(e - 63)` with `m` in `[2^63, 2^64)`. Zero has `m == 0`; infinity and NaN carry a sentinel `e`.
    e: i32,
    m: u64,
}

impl F80 {
    pub const ZERO: F80 = F80 {
        neg: false,
        e: 0,
        m: 0,
    };
    pub const ONE: F80 = F80 {
        neg: false,
        e: 0,
        m: 1 << 63,
    };
    pub const INFINITY: F80 = F80 {
        neg: false,
        e: E_INF,
        m: 1 << 63,
    };
    pub const NAN: F80 = F80 {
        neg: false,
        e: E_NAN,
        m: 3 << 62,
    };

    pub fn is_zero(self) -> bool {
        self.m == 0
    }

    pub fn is_nan(self) -> bool {
        self.e == E_NAN
    }

    pub fn is_infinite(self) -> bool {
        self.e == E_INF
    }

    pub fn is_sign_negative(self) -> bool {
        self.neg
    }

    /// The significand and the sign/exponent word as stored by `fstp tbyte` (biased exponent, sign in bit 15).
    pub fn repr(self) -> (u64, u16) {
        if self.is_nan() {
            return (3 << 62, 0x7fff);
        }
        if self.is_infinite() {
            return (1 << 63, 0x7fff | if self.neg { 0x8000 } else { 0 });
        }
        let s = if self.neg { 0x8000u16 } else { 0 };
        if self.m == 0 {
            return (0, s);
        }
        if self.e < E_MIN {
            return (self.m >> (E_MIN - self.e), s);
        }
        (self.m, s | (self.e + 16383) as u16)
    }

    pub fn from_repr(m: u64, se: u16) -> F80 {
        let neg = se & 0x8000 != 0;
        let be = (se & 0x7fff) as i32;
        if be == 0x7fff {
            return if m << 1 == 0 {
                F80 {
                    neg,
                    e: E_INF,
                    m: 1 << 63,
                }
            } else {
                F80::NAN
            };
        }
        if m == 0 {
            return F80 { neg, ..F80::ZERO };
        }
        // Denormals and pseudo-denormals: normalize (exact).
        let lz = m.leading_zeros() as i32;
        let e = if be == 0 { -16382 } else { be - 16383 } - lz;
        F80 { neg, e, m: m << lz }
    }

    /// Packs `mant * 2^exp` (plus a positive fraction of the last place when `sticky`) rounded to 64 bits.
    fn pack(neg: bool, mant: u128, exp: i32, sticky: bool) -> F80 {
        if mant == 0 {
            return F80 { neg, ..F80::ZERO };
        }
        let lz = mant.leading_zeros() as i32;
        let mant = mant << lz;
        let exp = exp - lz;
        let mut e = exp + 127;
        // Round once at the destination's precision. Subnormals have fewer
        // significand bits; rounding first to 64 bits would double-round.
        let (mut hi, round_up) = if e >= E_MIN {
            // Normal destinations always drop exactly 64 bits of the normalized significand.
            let hi = (mant >> 64) as u64;
            let rem = mant as u64;
            let half = 1u64 << 63;
            (hi, rem > half || (rem == half && (sticky || hi & 1 != 0)))
        } else {
            let shift = 64 + (E_MIN - e);
            if shift < 128 {
                let hi = (mant >> shift) as u64;
                let rem = mant & ((1u128 << shift) - 1);
                let half = 1u128 << (shift - 1);
                (hi, rem > half || (rem == half && (sticky || hi & 1 != 0)))
            } else if shift == 128 {
                (
                    0,
                    mant > (1u128 << 127) || (mant == (1u128 << 127) && sticky),
                )
            } else {
                (0, false)
            }
        };
        if round_up {
            let (h, carry) = hi.overflowing_add(1);
            if carry {
                hi = 1 << 63;
                e += 1;
            } else {
                hi = h;
            }
        }
        if e > E_MAX {
            return F80 {
                neg,
                e: E_INF,
                m: 1 << 63,
            };
        }
        if hi == 0 {
            return F80 { neg, ..F80::ZERO };
        }
        if e < E_MIN {
            let lz = hi.leading_zeros() as i32;
            return F80 {
                neg,
                e: E_MIN - lz,
                m: hi << lz,
            };
        }
        F80 { neg, e, m: hi }
    }

    pub fn from_f64(x: f64) -> F80 {
        let b = x.to_bits();
        let neg = b >> 63 != 0;
        let be = ((b >> 52) & 0x7ff) as i32;
        let f = b & ((1 << 52) - 1);
        if be == 0x7ff {
            return if f == 0 {
                F80 {
                    neg,
                    e: E_INF,
                    m: 1 << 63,
                }
            } else {
                F80::NAN
            };
        }
        if be == 0 {
            if f == 0 {
                return F80 { neg, ..F80::ZERO };
            }
            let lz = f.leading_zeros() as i32;
            // value = f * 2^-1074, f < 2^52
            return F80 {
                neg,
                e: 63 - lz - 1074,
                m: f << lz,
            };
        }
        F80 {
            neg,
            e: be - 1023,
            m: ((1u64 << 52) | f) << 11,
        }
    }

    pub fn from_i64(x: i64) -> F80 {
        if x == 0 {
            return F80::ZERO;
        }
        let neg = x < 0;
        let u = x.unsigned_abs();
        let lz = u.leading_zeros() as i32;
        F80 {
            neg,
            e: 63 - lz,
            m: u << lz,
        }
    }

    /// Round to the nearest double (ties to even); overflow gives infinity.
    pub fn to_f64(self) -> f64 {
        if self.is_nan() {
            return f64::NAN;
        }
        let sgn = if self.neg { -1.0 } else { 1.0 };
        if self.is_infinite() {
            return sgn * f64::INFINITY;
        }
        if self.m == 0 {
            return sgn * 0.0;
        }
        // Number of significand bits kept: 53 for normal doubles, fewer below 2^-1022.
        let keep = if self.e >= -1022 {
            53
        } else {
            53 - (-1022 - self.e)
        };
        if keep < 0 {
            return sgn * 0.0;
        }
        let drop = 64 - keep as u32;
        let (mut q, rem, half) = if drop >= 64 {
            (0u64, self.m, 1u64 << 63)
        } else {
            (
                self.m >> drop,
                self.m & ((1u64 << drop) - 1),
                1u64 << (drop - 1),
            )
        };
        if rem > half || (rem == half && q & 1 == 1) {
            q += 1;
        }
        // value = q * 2^(e - 63 + drop)
        let r = (q as f64) * pow2(self.e - 63 + drop as i32);
        sgn * r
    }

    pub fn abs(self) -> F80 {
        F80 { neg: false, ..self }
    }

    pub fn max(self, o: F80) -> F80 {
        if o > self {
            o
        } else {
            self
        }
    }

    pub fn min(self, o: F80) -> F80 {
        if o < self {
            o
        } else {
            self
        }
    }

    fn add_signed(self, o: F80, o_neg: bool) -> F80 {
        let (a, b) = (self, F80 { neg: o_neg, ..o });
        if a.is_nan() || b.is_nan() {
            return F80::NAN;
        }
        if a.is_infinite() || b.is_infinite() {
            if a.is_infinite() && b.is_infinite() && a.neg != b.neg {
                return F80::NAN;
            }
            return if a.is_infinite() { a } else { b };
        }
        if a.m == 0 && b.m == 0 {
            // -0 + -0 = -0, otherwise +0 (round to nearest).
            return F80 {
                neg: a.neg && b.neg,
                ..F80::ZERO
            };
        }
        if a.m == 0 {
            return b;
        }
        if b.m == 0 {
            return a;
        }
        // Order by magnitude.
        let (big, small) = if (a.e, a.m) >= (b.e, b.m) {
            (a, b)
        } else {
            (b, a)
        };
        let d = (big.e - small.e) as u32;
        let ba = (big.m as u128) << 63;
        let sm = (small.m as u128) << 63;
        let (bb, sticky) = if d == 0 {
            (sm, false)
        } else if d >= 128 {
            (0, true)
        } else {
            (sm >> d, sm & ((1u128 << d) - 1) != 0)
        };
        let exp = big.e - 126;
        if big.neg == small.neg {
            // Sum below 2^128: both are below 2^127.
            F80::pack(big.neg, ba + bb, exp, sticky)
        } else {
            let mut r = ba - bb;
            if sticky {
                // The true subtrahend is slightly larger than bb.
                r -= 1;
            }
            if r == 0 && !sticky {
                return F80::ZERO;
            }
            F80::pack(big.neg, r, exp, sticky)
        }
    }

    /// Fused-free multiply-subtract helpers are deliberately absent: the x87 rounds every operation.
    pub fn sqrt(self) -> F80 {
        if self.is_nan() || (self.neg && self.m != 0) {
            return F80::NAN;
        }
        if self.m == 0 || self.is_infinite() {
            return self;
        }
        // value = m * 2^t, t = e - 63. Make the exponent of the integer even and the integer 126..127 bits wide.
        let t = self.e - 63;
        let (big, tt) = if t & 1 == 0 {
            ((self.m as u128) << 62, t - 62)
        } else {
            ((self.m as u128) << 63, t - 63)
        };
        // big in [2^125, 2^127); isqrt in [2^62.5, 2^63.5) has 63 or 64 bits: take more bits by scaling 4x twice.
        // Compute the root of big << 64 (a 192-bit number) with the digit-by-digit method on two u128 halves.
        let (root, exact) = isqrt192(big);
        // value = (root + frac) * 2^(tt/2 - 32)
        F80::pack(false, root as u128, tt / 2 - 32, !exact)
    }
}

/// Integer square root of `x * 2^64` (a 192-bit radicand, 96 root bits): the floor and whether the square is exact.
fn isqrt192(x: u128) -> (u128, bool) {
    let words: [u64; 3] = [0, x as u64, (x >> 64) as u64];
    let bit = |p: usize| -> u128 { ((words[p / 64] >> (p % 64)) & 1) as u128 };
    let (mut rem, mut root) = (0u128, 0u128);
    for i in 0..96 {
        let p = 190 - 2 * i;
        rem = (rem << 2) | (bit(p + 1) << 1) | bit(p);
        let trial = (root << 2) | 1;
        root <<= 1;
        if rem >= trial {
            rem -= trial;
            root |= 1;
        }
    }
    (root, rem == 0)
}

fn pow2(e: i32) -> f64 {
    // Exact power of two for the exponents that matter here (the double exponent range).
    if e > 1023 {
        return f64::INFINITY;
    }
    if e >= -1022 {
        return f64::from_bits(((e + 1023) as u64) << 52);
    }
    if e >= -1074 {
        return f64::from_bits(1u64 << (e + 1074));
    }
    0.0
}

impl From<f64> for F80 {
    fn from(x: f64) -> F80 {
        F80::from_f64(x)
    }
}

impl Add for F80 {
    type Output = F80;

    fn add(self, o: F80) -> F80 {
        self.add_signed(o, o.neg)
    }
}

impl Sub for F80 {
    type Output = F80;

    fn sub(self, o: F80) -> F80 {
        self.add_signed(o, !o.neg)
    }
}

impl Mul for F80 {
    type Output = F80;

    fn mul(self, o: F80) -> F80 {
        let neg = self.neg != o.neg;
        if self.is_nan() || o.is_nan() {
            return F80::NAN;
        }
        if self.is_infinite() || o.is_infinite() {
            if self.m == 0 || o.m == 0 {
                return F80::NAN;
            }
            return F80 {
                neg,
                e: E_INF,
                m: 1 << 63,
            };
        }
        if self.m == 0 || o.m == 0 {
            return F80 { neg, ..F80::ZERO };
        }
        let product = (self.m as u128) * (o.m as u128);
        // Nonzero finite significands are normalized: their product has 127 or 128 bits.
        let top = product >> 127 != 0;
        let mut e = self.e + o.e + i32::from(top);
        // The general packer rounds directly to subnormal precision and handles overflow at the boundary.
        if e < E_MIN || e >= E_MAX {
            return F80::pack(neg, product, self.e + o.e - 126, false);
        }
        let (mut hi, lo) = if top {
            ((product >> 64) as u64, product as u64)
        } else {
            ((product >> 63) as u64, (product as u64) << 1)
        };
        let half = 1u64 << 63;
        if lo > half || (lo == half && hi & 1 != 0) {
            let (rounded, carry) = hi.overflowing_add(1);
            hi = if carry { 1 << 63 } else { rounded };
            e += i32::from(carry);
        }
        F80 { neg, e, m: hi }
    }
}

impl Div for F80 {
    type Output = F80;

    fn div(self, o: F80) -> F80 {
        let neg = self.neg != o.neg;
        if self.is_nan() || o.is_nan() {
            return F80::NAN;
        }
        if self.is_infinite() {
            return if o.is_infinite() {
                F80::NAN
            } else {
                F80 {
                    neg,
                    e: E_INF,
                    m: 1 << 63,
                }
            };
        }
        if o.is_infinite() {
            return F80 { neg, ..F80::ZERO };
        }
        if o.m == 0 {
            return if self.m == 0 {
                F80::NAN
            } else {
                F80 {
                    neg,
                    e: E_INF,
                    m: 1 << 63,
                }
            };
        }
        if self.m == 0 {
            return F80 { neg, ..F80::ZERO };
        }
        let (a, b) = (self.m as u128, o.m as u128);
        let n1 = a << 63;
        let q1 = n1 / b;
        let r1 = n1 % b;
        let n2 = r1 << 64;
        let q2 = n2 / b;
        let r2 = n2 % b;
        F80::pack(neg, (q1 << 64) | q2, self.e - o.e - 127, r2 != 0)
    }
}

impl Neg for F80 {
    type Output = F80;

    fn neg(self) -> F80 {
        F80 {
            neg: !self.neg,
            ..self
        }
    }
}

macro_rules! assign_op {
    ($tr:ident, $f:ident, $op:tt) => {
        impl $tr for F80 {
            fn $f(&mut self, o: F80) {
                *self = *self $op o;
            }
        }
    };
}

assign_op!(AddAssign, add_assign, +);
assign_op!(SubAssign, sub_assign, -);
assign_op!(MulAssign, mul_assign, *);
assign_op!(DivAssign, div_assign, /);

impl PartialEq for F80 {
    fn eq(&self, o: &F80) -> bool {
        self.partial_cmp(o) == Some(Ordering::Equal)
    }
}

impl PartialOrd for F80 {
    fn partial_cmp(&self, o: &F80) -> Option<Ordering> {
        if self.is_nan() || o.is_nan() {
            return None;
        }
        if self.m == 0 && o.m == 0 {
            return Some(Ordering::Equal);
        }
        let key = |x: &F80| -> (i32, i32, u64) {
            // Sign-aware key: (class, exponent, mantissa) with the sign applied by the caller.
            if x.m == 0 {
                (0, 0, 0)
            } else {
                (1, x.e, x.m)
            }
        };
        let (ka, kb) = (key(self), key(o));
        match (self.neg && self.m != 0, o.neg && o.m != 0) {
            (false, true) => Some(Ordering::Greater),
            (true, false) => Some(Ordering::Less),
            (false, false) => Some(ka.cmp(&kb)),
            (true, true) => Some(kb.cmp(&ka)),
        }
    }
}

impl Default for F80 {
    fn default() -> F80 {
        F80::ZERO
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn previous_pack(neg: bool, mant: u128, exp: i32, sticky: bool) -> F80 {
        if mant == 0 {
            return F80 { neg, ..F80::ZERO };
        }
        let lz = mant.leading_zeros() as i32;
        let mant = mant << lz;
        let exp = exp - lz;
        let mut e = exp + 127;
        // Round once at the destination's precision. Subnormals have fewer
        // significand bits; rounding first to 64 bits would double-round.
        let shift = 64 + (E_MIN - e).max(0);
        let (mut hi, round_up) = if shift < 128 {
            let hi = (mant >> shift) as u64;
            let rem = mant & ((1u128 << shift) - 1);
            let half = 1u128 << (shift - 1);
            (hi, rem > half || (rem == half && (sticky || hi & 1 != 0)))
        } else if shift == 128 {
            (
                0,
                mant > (1u128 << 127) || (mant == (1u128 << 127) && sticky),
            )
        } else {
            (0, false)
        };
        if round_up {
            let (h, carry) = hi.overflowing_add(1);
            if carry {
                hi = 1 << 63;
                e += 1;
            } else {
                hi = h;
            }
        }
        if e > E_MAX {
            return F80 {
                neg,
                e: E_INF,
                m: 1 << 63,
            };
        }
        if hi == 0 {
            return F80 { neg, ..F80::ZERO };
        }
        if e < E_MIN {
            let lz = hi.leading_zeros() as i32;
            return F80 {
                neg,
                e: E_MIN - lz,
                m: hi << lz,
            };
        }
        F80 { neg, e, m: hi }
    }

    // Frozen multiply before the normal-product specialization; keep the general pack path as the oracle.
    fn previous_mul(a: F80, o: F80) -> F80 {
        let neg = a.neg != o.neg;
        if a.is_nan() || o.is_nan() {
            return F80::NAN;
        }
        if a.is_infinite() || o.is_infinite() {
            if a.m == 0 || o.m == 0 {
                return F80::NAN;
            }
            return F80 {
                neg,
                e: E_INF,
                m: 1 << 63,
            };
        }
        if a.m == 0 || o.m == 0 {
            return F80 { neg, ..F80::ZERO };
        }
        previous_pack(neg, (a.m as u128) * (o.m as u128), a.e + o.e - 126, false)
    }

    #[test]
    fn pack_matches_previous() {
        let check = |neg, mant, exp, sticky| {
            assert_eq!(
                F80::pack(neg, mant, exp, sticky).repr(),
                previous_pack(neg, mant, exp, sticky).repr(),
                "{neg} {mant:x} {exp} {sticky}"
            );
        };
        for mant in [
            0,
            1,
            2,
            1 << 63,
            1 << 64,
            (1 << 127) - 1,
            1 << 127,
            (1 << 127) + (1 << 63) - 1,
            (1 << 127) + (1 << 63),
            (1 << 127) + (1 << 63) + 1,
            u128::MAX,
        ] {
            for exp in [
                E_MIN - 255,
                E_MIN - 192,
                E_MIN - 128,
                E_MIN - 127,
                E_MIN - 126,
                E_MIN - 64,
                -127,
                -126,
                0,
                E_MAX - 128,
                E_MAX - 127,
                E_MAX - 126,
                E_MAX,
                E_MAX + 128,
            ] {
                for neg in [false, true] {
                    for sticky in [false, true] {
                        check(neg, mant, exp, sticky);
                    }
                }
            }
        }
        // Column values from the first fx_quad_obj_5371 QTAU dump, in the add/pack significand scale.
        for (mant, se) in [
            (0xf424_0000_0000_0000, 0xbffe),
            (0xb71b_0000_0000_0000, 0xbffd),
            (0xea60_0000_0000_0000, 0x3ff7),
            (0x8000_0000_0000_0000, 0xbff1),
            (0x9c40_0000_0000_0000, 0xbff6),
            (0xc000_0000_0000_0000, 0x3ffe),
        ] {
            let v = F80::from_repr(mant, se);
            check(v.neg, (v.m as u128) << 63, v.e - 126, false);
        }
        let mut seed = 0x198d_ac52_6214_be9du64;
        let mut next = || {
            seed ^= seed << 13;
            seed ^= seed >> 7;
            seed ^= seed << 17;
            seed
        };
        for _ in 0..100_000 {
            let mant = (((next() as u128) << 64) | next() as u128) >> (next() & 127);
            let exp = (next() % 66_001) as i32 - 33_000;
            let flags = next();
            check(flags & 1 != 0, mant, exp, flags & 2 != 0);
        }
    }

    #[test]
    fn multiply_matches_previous() {
        let mut values = vec![
            F80::ZERO,
            -F80::ZERO,
            F80::ONE,
            F80::INFINITY,
            -F80::INFINITY,
            F80::NAN,
        ];
        for se in [0, 1, 2, 0x3ffe, 0x3fff, 0x4000, 0x7ffd, 0x7ffe] {
            for m in [
                1,
                0x7fff_ffff_ffff_ffff,
                0x8000_0000_0000_0000,
                0x8000_0000_0000_0001,
                0xffff_ffff_ffff_ffff,
            ] {
                for sign in [0, 0x8000] {
                    values.push(F80::from_repr(m, se | sign));
                }
            }
        }
        // First QTAU column from fx_quad_obj_5371.dual-r, ROSE_DUMP=qtau, max-iterations=5400.
        for (m, se) in [
            (0xf424_0000_0000_0000, 0xbffe),
            (0xb71b_0000_0000_0000, 0xbffd),
            (0xea60_0000_0000_0000, 0x3ff7),
            (0x8000_0000_0000_0000, 0xbff1),
            (0x8000_0000_0000_0000, 0x3ffe),
            (0x9c40_0000_0000_0000, 0xbff6),
            (0xf424_0000_0000_0000, 0x3ffe),
            (0x8000_0000_0000_0000, 0xbfeb),
            (0xc000_0000_0000_0000, 0x3ffe),
            (0x8000_0000_0000_0000, 0x3ffd),
        ] {
            values.push(F80::from_repr(m, se));
        }
        let check = |a: F80, b: F80| {
            assert_eq!(
                (a * b).repr(),
                previous_mul(a, b).repr(),
                "{:?} * {:?}",
                a.repr(),
                b.repr()
            )
        };
        for &a in &values {
            for &b in &values {
                check(a, b);
            }
        }
        let mut seed = 0xe5a1_d839_1208_759bu64;
        let mut next = || {
            seed ^= seed << 13;
            seed ^= seed >> 7;
            seed ^= seed << 17;
            seed
        };
        for _ in 0..100_000 {
            let a = F80::from_repr(next(), next() as u16);
            let b = F80::from_repr(next(), next() as u16);
            check(a, b);
        }
    }

    #[test]
    fn roundtrips_doubles() {
        for &x in &[
            0.0,
            1.0,
            -1.0,
            0.1,
            1e-300,
            3.5e300,
            -2.5e-310,
            f64::MIN_POSITIVE,
            123456.789,
        ] {
            assert_eq!(F80::from_f64(x).to_f64().to_bits(), x.to_bits(), "{x}");
        }
    }

    #[test]
    fn small_arithmetic() {
        let a = F80::from_f64(1.0) / F80::from_f64(3.0);
        let (m, se) = a.repr();
        assert_eq!((m, se), (0xaaaa_aaaa_aaaa_aaab, 0x3ffd));
        let s = F80::from_f64(2.0).sqrt();
        assert_eq!(s.repr(), (0xb504_f333_f9de_6484, 0x3fff));
        let t = F80::from_f64(0.1) + F80::from_f64(0.2);
        assert_eq!(t.to_f64(), 0.30000000000000004);
    }

    #[cfg(all(target_arch = "x86_64", target_os = "linux"))]
    mod hw {
        use super::*;
        use std::arch::asm;

        fn bytes(x: F80) -> [u8; 16] {
            let (m, se) = x.repr();
            let mut b = [0u8; 16];
            b[..8].copy_from_slice(&m.to_le_bytes());
            b[8..10].copy_from_slice(&se.to_le_bytes());
            b
        }

        fn unbytes(b: [u8; 16]) -> F80 {
            F80::from_repr(
                u64::from_le_bytes(b[..8].try_into().unwrap()),
                u16::from_le_bytes(b[8..10].try_into().unwrap()),
            )
        }

        fn op(name: u8, a: F80, b: F80) -> F80 {
            let (ba, bb) = (bytes(a), bytes(b));
            let mut out = [0u8; 16];
            unsafe {
                match name {
                    b'+' => {
                        asm!("fld tbyte ptr [{a}]", "fld tbyte ptr [{b}]", "faddp st(1), st", "fstp tbyte ptr [{o}]",
                        a = in(reg) ba.as_ptr(), b = in(reg) bb.as_ptr(), o = in(reg) out.as_mut_ptr(), out("st(0)") _, out("st(1)") _, options(nostack))
                    }
                    b'-' => {
                        asm!("fld tbyte ptr [{a}]", "fld tbyte ptr [{b}]", "fsubp st(1), st", "fstp tbyte ptr [{o}]",
                        a = in(reg) ba.as_ptr(), b = in(reg) bb.as_ptr(), o = in(reg) out.as_mut_ptr(), out("st(0)") _, out("st(1)") _, options(nostack))
                    }
                    b'*' => {
                        asm!("fld tbyte ptr [{a}]", "fld tbyte ptr [{b}]", "fmulp st(1), st", "fstp tbyte ptr [{o}]",
                        a = in(reg) ba.as_ptr(), b = in(reg) bb.as_ptr(), o = in(reg) out.as_mut_ptr(), out("st(0)") _, out("st(1)") _, options(nostack))
                    }
                    b'/' => {
                        asm!("fld tbyte ptr [{a}]", "fld tbyte ptr [{b}]", "fdivp st(1), st", "fstp tbyte ptr [{o}]",
                        a = in(reg) ba.as_ptr(), b = in(reg) bb.as_ptr(), o = in(reg) out.as_mut_ptr(), out("st(0)") _, out("st(1)") _, options(nostack))
                    }
                    b'q' => asm!("fld tbyte ptr [{a}]", "fsqrt", "fstp tbyte ptr [{o}]",
                        a = in(reg) ba.as_ptr(), o = in(reg) out.as_mut_ptr(), out("st(0)") _, out("st(1)") _, options(nostack)),
                    _ => unreachable!(),
                }
            }
            unbytes(out)
        }

        fn rnd(s: &mut u64) -> F80 {
            *s = s
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            let m = (*s >> 11) | (1 << 63);
            *s = s
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            let e = 16383 - 40 + ((*s >> 40) % 80) as u16;
            let neg = (*s >> 20) & 1 == 1;
            // Sometimes use short significands (many trailing zeros) to hit ties and exact cancellations.
            let m = match (*s >> 10) % 4 {
                0 => m & !((1u64 << 40) - 1),
                1 => m & !((1u64 << 11) - 1),
                _ => m,
            };
            F80::from_repr(m, e | if neg { 0x8000 } else { 0 })
        }

        // Set the same 64-bit precision and round-to-nearest mode as the
        // reference; restore the calling thread's control word on exit.
        struct ControlWord(u16);

        impl ControlWord {
            fn nearest_extended() -> Self {
                let mut old = 0u16;
                let cw = 0x037fu16;
                unsafe {
                    asm!("fnstcw [{old}]", "fldcw [{cw}]", old = in(reg) &mut old, cw = in(reg) &cw, options(nostack));
                }
                Self(old)
            }
        }

        impl Drop for ControlWord {
            fn drop(&mut self) {
                unsafe {
                    asm!("fldcw [{cw}]", cw = in(reg) &self.0, options(nostack));
                }
            }
        }

        fn hw_to_f64(x: F80) -> f64 {
            let a = bytes(x);
            let mut out = 0f64;
            unsafe {
                asm!("fld tbyte ptr [{a}]", "fstp qword ptr [{o}]", a = in(reg) a.as_ptr(), o = in(reg) &mut out, out("st(0)") _, options(nostack));
            }
            out
        }

        #[test]
        fn conversion_boundaries_match_x87() {
            let _cw = ControlWord::nearest_extended();
            for bits in [
                0,
                1,
                2,
                (1u64 << 52) - 1,
                1u64 << 52,
                (1u64 << 52) + 1,
                0x7fef_ffff_ffff_ffff,
                0x7ff0_0000_0000_0000,
            ] {
                for sign in [0, 1u64 << 63] {
                    let x = f64::from_bits(bits | sign);
                    assert_eq!(F80::from_f64(x).to_f64().to_bits(), x.to_bits());
                }
            }
            // Double normal/subnormal transition, half-way ties and overflow.
            for e in [
                -16445, -16383, -1076, -1075, -1074, -1023, -1022, 1023, 1024,
            ] {
                for m in [
                    1u64 << 63,
                    (1u64 << 63) + 1,
                    (1u64 << 63) + 1024,
                    (1u64 << 63) + 1025,
                    u64::MAX,
                ] {
                    let x = F80::pack(false, m as u128, e - 63, false);
                    for x in [x, -x] {
                        assert_eq!(
                            x.to_f64().to_bits(),
                            hw_to_f64(x).to_bits(),
                            "{:?}",
                            x.repr()
                        );
                    }
                }
            }
        }

        #[test]
        fn subnormal_boundaries_match_x87() {
            let _cw = ControlWord::nearest_extended();
            let values = [
                (0, 0),
                (0, 0x8000),
                (1, 0),
                (2, 0),
                (3, 0),
                (0x3fff_ffff_ffff_ffff, 0),
                (0x4000_0000_0000_0000, 0),
                (0x7fff_ffff_ffff_fffe, 0),
                (0x7fff_ffff_ffff_ffff, 0),
                (0x8000_0000_0000_0000, 1),
                (0x8000_0000_0000_0001, 1),
                (0xffff_ffff_ffff_ffff, 1),
                (0x8000_0000_0000_0000, 2),
                (0x8000_0000_0000_0000, 0x3ffe),
                (0xc000_0000_0000_0000, 0x3ffe),
                (0x8000_0000_0000_0000, 0x3fff),
                (0xc000_0000_0000_0000, 0x3fff),
                (0x8000_0000_0000_0000, 0x4000),
                (0xffff_ffff_ffff_ffff, 0x7ffe),
                (0x8000_0000_0000_0000, 0x7fff),
            ];
            for &(m, se) in &values {
                for sign in [0, 0x8000] {
                    let a = F80::from_repr(m, se ^ sign);
                    assert_eq!(a.repr(), (m, se ^ sign));
                    if !a.is_sign_negative() {
                        let h = op(b'q', a, a);
                        assert_eq!(h.repr(), a.sqrt().repr(), "sqrt {:?}", a.repr());
                    }
                    for &(n, te) in &values {
                        for sign2 in [0, 0x8000] {
                            let b = F80::from_repr(n, te ^ sign2);
                            for (c, soft) in
                                [(b'+', a + b), (b'-', a - b), (b'*', a * b), (b'/', a / b)]
                            {
                                let h = op(c, a, b);
                                if h.is_nan() {
                                    assert!(soft.is_nan());
                                } else {
                                    assert_eq!(
                                        h.repr(),
                                        soft.repr(),
                                        "{} {:?} {:?}",
                                        c as char,
                                        a.repr(),
                                        b.repr()
                                    );
                                }
                            }
                        }
                    }
                }
            }
            // Exact, half-way and sticky products on both sides of the
            // subnormal-to-normal boundary, plus very small round-to-zero.
            let mut seed = 11u64;
            for _ in 0..20_000 {
                seed = seed.wrapping_mul(6364136223846793005).wrapping_add(1);
                let a = F80::from_repr((seed & 0x7fff_ffff_ffff_ffff) | 1, 0);
                seed = seed.wrapping_mul(6364136223846793005).wrapping_add(1);
                let b = F80::from_repr(seed | (1 << 63), 0x3ffd + (seed % 4) as u16);
                for (c, soft) in [(b'+', a + b), (b'-', a - b), (b'*', a * b), (b'/', a / b)] {
                    assert_eq!(
                        op(c, a, b).repr(),
                        soft.repr(),
                        "{} {:?} {:?}",
                        c as char,
                        a.repr(),
                        b.repr()
                    );
                }
            }
        }

        #[test]
        fn matches_the_x87_unit() {
            let _cw = ControlWord::nearest_extended();
            let mut s = 12345u64;
            for i in 0..200_000 {
                let (a, b) = (rnd(&mut s), rnd(&mut s));
                // Exact cancellation and near-equal operands.
                let b = if i % 7 == 0 {
                    a
                } else if i % 11 == 0 {
                    F80::from_repr(a.repr().0 ^ 1, a.repr().1)
                } else {
                    b
                };
                for (c, soft) in [(b'+', a + b), (b'-', a - b), (b'*', a * b), (b'/', a / b)] {
                    let h = op(c, a, b);
                    assert_eq!(
                        h.repr(),
                        soft.repr(),
                        "{} {:?} {:?}",
                        c as char,
                        a.repr(),
                        b.repr()
                    );
                }
                let h = op(b'q', a.abs(), a);
                assert_eq!(h.repr(), a.abs().sqrt().repr(), "sqrt {:?}", a.repr());
            }
        }
    }
}
