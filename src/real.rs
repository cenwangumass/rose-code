//! The number type of the LU, FTRAN/BTRAN and dual-simplex arithmetic: `f64`, or `F80` for the reference's quad mode (x87 `long double`).
//!
//! The reference's dual is one template instantiated for `double` and for `long double`; the same source gives the same operations in the
//! same order, so the code that models it is generic over `Real`. A double constant meeting a `Real` value is converted exactly (`c`), as the
//! x87 loads of a `double` operand do; a stored `Real` is rounded back to a double only where the reference stores into a `double` (`to_f64`).

use crate::f80::F80;
use std::fmt::Debug;
use std::ops::{Add, AddAssign, Div, DivAssign, Mul, MulAssign, Neg, Sub, SubAssign};

pub trait Real:
    Copy
    + Debug
    + Default
    + PartialOrd
    + PartialEq
    + Add<Output = Self>
    + Sub<Output = Self>
    + Mul<Output = Self>
    + Div<Output = Self>
    + Neg<Output = Self>
    + AddAssign
    + SubAssign
    + MulAssign
    + DivAssign
    + 'static
{
    const ZERO: Self;
    const ONE: Self;
    /// The SSE2 kernels of the double code add in vector lanes; the x87 code of the quad mode adds one term after the other.
    const LANES: bool;

    /// Exact conversion of a double.
    fn c(x: f64) -> Self;

    /// Rounding to the nearest double.
    fn to_f64(self) -> f64;

    fn abs(self) -> Self;

    fn max(self, o: Self) -> Self;

    fn min(self, o: Self) -> Self;

    fn is_nan(self) -> bool;

    fn from_usize(n: usize) -> Self {
        Self::c(n as f64)
    }
}

impl Real for f64 {
    const ZERO: f64 = 0.0;
    const ONE: f64 = 1.0;
    const LANES: bool = true;

    #[inline(always)]
    fn c(x: f64) -> f64 {
        x
    }

    #[inline(always)]
    fn to_f64(self) -> f64 {
        self
    }

    #[inline(always)]
    fn abs(self) -> f64 {
        f64::abs(self)
    }

    #[inline(always)]
    fn max(self, o: f64) -> f64 {
        f64::max(self, o)
    }

    #[inline(always)]
    fn min(self, o: f64) -> f64 {
        f64::min(self, o)
    }

    #[inline(always)]
    fn is_nan(self) -> bool {
        f64::is_nan(self)
    }
}

impl Real for F80 {
    const ZERO: F80 = F80::ZERO;
    const ONE: F80 = F80::ONE;
    const LANES: bool = false;

    fn c(x: f64) -> F80 {
        F80::from_f64(x)
    }

    fn to_f64(self) -> f64 {
        F80::to_f64(self)
    }

    fn abs(self) -> F80 {
        F80::abs(self)
    }

    fn max(self, o: F80) -> F80 {
        if self.is_nan() {
            o
        } else if o.is_nan() || self >= o {
            self
        } else {
            o
        }
    }

    fn min(self, o: F80) -> F80 {
        if self.is_nan() {
            o
        } else if o.is_nan() || self <= o {
            self
        } else {
            o
        }
    }

    fn is_nan(self) -> bool {
        F80::is_nan(self)
    }
}
