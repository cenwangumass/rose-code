//! The A and B lists of the tail columns with their coefficients, in one shared pair of arrays.
//!
//! A tail column's numeric U entries split into those whose row is a tail position (`kind == -2`, list A) and those whose row is in
//! the nucleus list (`kind >= 0`, list B). The reference keeps both in one array, A first and B after it (`H+0x68/0x70/0x78`),
//! and every entry carries its coefficient: the exact bits of the matching `uc[c]` entry. List order is behaviour (the work-estimate
//! loops and the dependency walk read it), so the in-place edits below reproduce `Vec::swap_remove`, `insert(0, _)` and `remove(0)`
//! on the two lists through the slot the boundary between them gives up or takes.

use crate::real::Real;

/// Where a column's lists live: `na` entries of A and then `nb` of B from `start`, with room for `cap` entries.
#[derive(Clone, Copy, Default)]
struct Span {
    start: u32,
    na: u32,
    nb: u32,
    cap: u32,
}

#[derive(Clone)]
pub(super) struct AbLists<R> {
    span: Vec<Span>,
    idx: Vec<u32>,
    val: Vec<R>,
}

impl<R: Real> AbLists<R> {
    pub fn new() -> Self {
        Self { span: Vec::new(), idx: Vec::new(), val: Vec::new() }
    }

    /// `m` empty columns; the arrays keep their allocations.
    pub fn reset(&mut self, m: usize) {
        self.span.clear();
        self.span.resize(m, Span::default());
        self.idx.clear();
        self.val.clear();
    }

    #[inline(always)]
    pub fn na(&self, c: usize) -> usize {
        self.span[c].na as usize
    }

    #[inline(always)]
    pub fn nb(&self, c: usize) -> usize {
        self.span[c].nb as usize
    }

    /// `(na, nb)` from one lookup of the column's span.
    #[inline(always)]
    pub fn counts(&self, c: usize) -> (usize, usize) {
        let s = &self.span[c];
        (s.na as usize, s.nb as usize)
    }

    #[inline(always)]
    pub fn is_empty(&self, c: usize) -> bool {
        let s = &self.span[c];
        s.na == 0 && s.nb == 0
    }

    /// A list of column `c` (row positions).
    #[inline(always)]
    pub fn a(&self, c: usize) -> &[u32] {
        let s = &self.span[c];
        &self.idx[s.start as usize..(s.start + s.na) as usize]
    }

    /// The A list and its coefficients, from one lookup of the column's span.
    #[inline(always)]
    pub fn a_pair(&self, c: usize) -> (&[u32], &[R]) {
        let s = &self.span[c];
        let (lo, hi) = (s.start as usize, (s.start + s.na) as usize);
        #[cfg(any(test, debug_assertions, feature = "factor-shadow"))]
        {
            (&self.idx[lo..hi], &self.val[lo..hi])
        }
        #[cfg(not(any(test, debug_assertions, feature = "factor-shadow")))]
        {
            // Every span lies inside both arrays (`idx.len() == val.len()`, `start + na + nb <= start + cap <= len`): every edit above
            // writes through checked indices into room it has sized, and the test, debug and shadow builds check it on each read.
            // SAFETY: that invariant.
            unsafe { (self.idx.get_unchecked(lo..hi), self.val.get_unchecked(lo..hi)) }
        }
    }

    /// A then B and their coefficients, from one lookup of the column's span: the whole numeric U column, entry for entry.
    #[inline(always)]
    pub fn both_pair(&self, c: usize) -> (&[u32], &[R]) {
        let s = &self.span[c];
        let (lo, hi) = (s.start as usize, (s.start + s.na + s.nb) as usize);
        #[cfg(any(test, debug_assertions, feature = "factor-shadow"))]
        {
            (&self.idx[lo..hi], &self.val[lo..hi])
        }
        #[cfg(not(any(test, debug_assertions, feature = "factor-shadow")))]
        {
            // SAFETY: as `a_pair` (`start + na + nb <= start + cap <= len` of both arrays).
            unsafe { (self.idx.get_unchecked(lo..hi), self.val.get_unchecked(lo..hi)) }
        }
    }

    /// `(na, nb)` and the whole numeric U column (A then B with the coefficients) from one lookup of the column's span.
    #[inline(always)]
    pub fn counts_both_pair(&self, c: usize) -> (usize, usize, &[u32], &[R]) {
        let s = &self.span[c];
        let (lo, hi) = (s.start as usize, (s.start + s.na + s.nb) as usize);
        #[cfg(any(test, debug_assertions, feature = "factor-shadow"))]
        {
            (s.na as usize, s.nb as usize, &self.idx[lo..hi], &self.val[lo..hi])
        }
        #[cfg(not(any(test, debug_assertions, feature = "factor-shadow")))]
        {
            // SAFETY: as `both_pair`.
            unsafe { (s.na as usize, s.nb as usize, self.idx.get_unchecked(lo..hi), self.val.get_unchecked(lo..hi)) }
        }
    }

    /// `x[j] -= v * xv` over the A list of column `c`, entry by entry in list order (4-way unrolled, same operations and order as the
    /// plain loop even if an index repeats). Every A entry is a position below `x.len()` (asserted by the caller once per call;
    /// checked per entry in the test, debug and shadow builds).
    #[inline(always)]
    pub fn scatter_a(&self, c: usize, x: &mut [R], xv: R) {
        let (idx, val) = self.a_pair(c);
        #[cfg(any(test, debug_assertions, feature = "factor-shadow"))]
        {
            for (&j, &v) in idx.iter().zip(val) {
                x[j as usize] -= v * xv;
            }
        }
        #[cfg(not(any(test, debug_assertions, feature = "factor-shadow")))]
        {
            // SAFETY: the caller's invariant: every A entry is below `x.len()`.
            unsafe { scatter_unrolled(idx, val, x, xv) }
        }
    }

    /// Coefficients of the A list, entry for entry.
    #[cfg(any(test, debug_assertions, feature = "factor-shadow"))]
    #[inline(always)]
    pub fn a_val(&self, c: usize) -> &[R] {
        let s = &self.span[c];
        &self.val[s.start as usize..(s.start + s.na) as usize]
    }

    #[inline(always)]
    pub fn b(&self, c: usize) -> &[u32] {
        let s = &self.span[c];
        &self.idx[(s.start + s.na) as usize..(s.start + s.na + s.nb) as usize]
    }

    #[cfg(any(test, debug_assertions, feature = "factor-shadow"))]
    #[inline(always)]
    pub fn b_val(&self, c: usize) -> &[R] {
        let s = &self.span[c];
        &self.val[(s.start + s.na) as usize..(s.start + s.na + s.nb) as usize]
    }

    /// The coefficients of A then B, entry for entry with `both`.
    #[inline(always)]
    pub fn both_val(&self, c: usize) -> &[R] {
        let s = &self.span[c];
        &self.val[s.start as usize..(s.start + s.na + s.nb) as usize]
    }

    /// A then B.
    #[inline(always)]
    pub fn both(&self, c: usize) -> &[u32] {
        let s = &self.span[c];
        &self.idx[s.start as usize..(s.start + s.na + s.nb) as usize]
    }

    pub fn a_pos(&self, c: usize, x: usize) -> Option<usize> {
        self.a(c).iter().position(|&y| y as usize == x)
    }

    pub fn b_pos(&self, c: usize, x: usize) -> Option<usize> {
        self.b(c).iter().position(|&y| y as usize == x)
    }

    /// The lists of column `c`, which is empty and has no room yet (a column `fresh` has not touched), appended at the end of the
    /// arrays with exactly their length as room: the layout `set` gives such a column, without zero-filling the new room first.
    pub fn append_fresh(&mut self, c: usize, a: impl Iterator<Item = (u32, R)>, b: &[(u32, R)]) {
        debug_assert!(self.span[c].cap == 0 && self.span[c].na == 0 && self.span[c].nb == 0);
        let start = self.idx.len();
        for (x, v) in a {
            self.idx.push(x);
            self.val.push(v);
        }
        let na = self.idx.len() - start;
        for &(x, v) in b {
            self.idx.push(x);
            self.val.push(v);
        }
        let n = self.idx.len() - start;
        if n > 0 {
            self.span[c] = Span { start: start as u32, na: na as u32, nb: (n - na) as u32, cap: n as u32 };
        }
    }

    /// Replace both lists of column `c`: in place when they fit its room, else at the end of the arrays (Rose-internal; not projected).
    pub fn set(&mut self, c: usize, a: &[(u32, R)], b: &[(u32, R)]) {
        let n = a.len() + b.len();
        let s = &mut self.span[c];
        if n > s.cap as usize {
            s.start = self.idx.len() as u32;
            s.cap = n as u32;
            let end = self.idx.len() + n;
            self.idx.resize(end, 0);
            self.val.resize(end, R::ZERO);
        }
        s.na = a.len() as u32;
        s.nb = b.len() as u32;
        let base = s.start as usize;
        for (i, &(x, v)) in a.iter().chain(b).enumerate() {
            self.idx[base + i] = x;
            self.val[base + i] = v;
        }
    }

    /// Both lists of `c` become empty (the room stays).
    pub fn clear(&mut self, c: usize) {
        let s = &mut self.span[c];
        s.na = 0;
        s.nb = 0;
    }

    #[inline(always)]
    fn copy(&mut self, from: usize, to: usize) {
        self.idx[to] = self.idx[from];
        self.val[to] = self.val[from];
    }

    /// `A.swap_remove(p)` followed by `B.insert(0, removed)`: the removed entry and the last A entry trade places, and the boundary moves
    /// below the removed entry, which becomes the front of B.
    pub fn demote(&mut self, c: usize, p: usize) {
        let s = &mut self.span[c];
        debug_assert!(p < s.na as usize);
        let (base, last) = (s.start as usize, (s.start + s.na - 1) as usize);
        s.na -= 1;
        s.nb += 1;
        self.idx.swap(base + p, last);
        self.val.swap(base + p, last);
    }

    /// `A.swap_remove(p)` followed by moving B's last entry to its front (`pop` then `insert(0, _)`), dropping the A entry: the last A entry
    /// fills the hole, and B's last entry takes the slot A gave up.
    pub fn drop_a(&mut self, c: usize, p: usize) {
        let s = &mut self.span[c];
        debug_assert!(p < s.na as usize);
        let (base, last_a, last_b) = (s.start as usize, (s.start + s.na - 1) as usize, (s.start + s.na + s.nb).wrapping_sub(1) as usize);
        let nb = s.nb;
        s.na -= 1;
        self.copy(last_a, base + p);
        if nb > 0 {
            self.copy(last_b, last_a);
        }
    }

    /// `B.swap_remove(p)`.
    pub fn drop_b(&mut self, c: usize, p: usize) {
        let s = &mut self.span[c];
        debug_assert!(p < s.nb as usize);
        let (first, last) = ((s.start + s.na) as usize, (s.start + s.na + s.nb - 1) as usize);
        s.nb -= 1;
        self.copy(last, first + p);
    }

    /// The entry at `B[p]` leaves B the way `remove_front_fill` removes it (B's front entry fills the hole, the front slot is dropped) and
    /// is pushed to the end of A, which is that same front slot.
    pub fn b_to_a(&mut self, c: usize, p: usize) {
        let s = &mut self.span[c];
        debug_assert!(p < s.nb as usize);
        let front = (s.start + s.na) as usize;
        s.na += 1;
        s.nb -= 1;
        if p != 0 {
            self.idx.swap(front, front + p);
            self.val.swap(front, front + p);
        }
    }

    /// The lists as nested vectors `(A, B)`, for the oracle projection.
    #[cfg(any(test, debug_assertions, feature = "factor-shadow"))]
    pub fn to_vecs(&self) -> (Vec<Vec<usize>>, Vec<Vec<usize>>) {
        let m = self.span.len();
        let a = (0..m).map(|c| self.a(c).iter().map(|&x| x as usize).collect()).collect();
        let b = (0..m).map(|c| self.b(c).iter().map(|&x| x as usize).collect()).collect();
        (a, b)
    }

    /// Every list and coefficient against the numeric columns `uc` (the audited invariants): A holds the rows of `uc[c]` with `kind == -2`,
    /// B those with `kind >= 0`, no duplicates, and the coefficients carry the bits of the matching `uc` entry.
    #[cfg(any(test, debug_assertions, feature = "factor-shadow"))]
    pub fn check(&self, uc: &[Vec<(usize, R)>], kind: &[i64]) {
        assert_eq!(self.span.len(), uc.len());
        for c in 0..uc.len() {
            let (a, b) = (self.a(c), self.b(c));
            let (av, bv) = (self.a_val(c), self.b_val(c));
            let cols = |sel: &dyn Fn(i64) -> bool| uc[c].iter().filter(|e| sel(kind[e.0])).map(|e| e.0).collect::<std::collections::BTreeSet<_>>();
            let in_a = cols(&|k| k == -2);
            let in_b = cols(&|k| k >= 0);
            assert_eq!(a.iter().map(|&x| x as usize).collect::<std::collections::BTreeSet<_>>(), in_a, "A set of column {c}");
            assert_eq!(b.iter().map(|&x| x as usize).collect::<std::collections::BTreeSet<_>>(), in_b, "B set of column {c}");
            assert_eq!(a.len(), in_a.len(), "duplicate in A of column {c}");
            assert_eq!(b.len(), in_b.len(), "duplicate in B of column {c}");
            for (list, vals) in [(a, av), (b, bv)] {
                for (&x, &v) in list.iter().zip(vals) {
                    let u = uc[c].iter().find(|e| e.0 == x as usize).expect("entry in uc").1;
                    assert!(same_bits(u, v), "coefficient of row {x} in column {c}");
                }
            }
        }
    }
}

/// Bitwise equality of two values (signed zeros and NaN payloads included; the x87 format by its representation).
#[cfg(any(test, debug_assertions, feature = "factor-shadow"))]
pub(super) fn same_bits<R: Real>(a: R, b: R) -> bool {
    use std::any::Any;
    let (a, b) = (&a as &dyn Any, &b as &dyn Any);
    if let Some(a) = a.downcast_ref::<f64>() {
        a.to_bits() == b.downcast_ref::<f64>().unwrap().to_bits()
    } else {
        let (a, b) = (a.downcast_ref::<crate::f80::F80>().unwrap(), b.downcast_ref::<crate::f80::F80>().unwrap());
        a.repr() == b.repr()
    }
}

/// The representation of a value (signed zeros, NaN payloads and the x87 format included), for multiset comparisons.
#[cfg(feature = "factor-shadow")]
pub(super) fn bits_of<R: Real>(a: R) -> Vec<u8> {
    use std::any::Any;
    let a = &a as &dyn Any;
    if let Some(a) = a.downcast_ref::<f64>() {
        a.to_bits().to_le_bytes().to_vec()
    } else {
        format!("{:?}", a.downcast_ref::<crate::f80::F80>().unwrap().repr()).into_bytes()
    }
}

/// `x[idx[k]] -= val[k] * xv` for `k` in order, 4-way unrolled with the entries in textual order (an index may repeat).
///
/// # Safety
/// Every `idx` entry is below `x.len()`.
#[inline(always)]
#[cfg_attr(any(test, debug_assertions, feature = "factor-shadow"), allow(dead_code))]
pub(super) unsafe fn scatter_unrolled<R: Real>(idx: &[u32], val: &[R], x: &mut [R], xv: R) {
    let n = idx.len().min(val.len());
    let mut k = 0;
    // SAFETY: `k + r < n <= idx.len(), val.len()`; destinations are in range (caller).
    unsafe {
        while k + 4 <= n {
            *x.get_unchecked_mut(*idx.get_unchecked(k) as usize) -= *val.get_unchecked(k) * xv;
            *x.get_unchecked_mut(*idx.get_unchecked(k + 1) as usize) -= *val.get_unchecked(k + 1) * xv;
            *x.get_unchecked_mut(*idx.get_unchecked(k + 2) as usize) -= *val.get_unchecked(k + 2) * xv;
            *x.get_unchecked_mut(*idx.get_unchecked(k + 3) as usize) -= *val.get_unchecked(k + 3) * xv;
            k += 4;
        }
        while k < n {
            *x.get_unchecked_mut(*idx.get_unchecked(k) as usize) -= *val.get_unchecked(k) * xv;
            k += 1;
        }
    }
}
