//! The numeric U rows in one pool: the row-wise copy of U (`H+0x90..0xb0` in the reference).
//!
//! Row `y` lists the tail columns `t` with a numeric U entry at `y` together with its value, in the reference's physical row order
//! (the order in which `fresh` laid the row out, then `update` appends at the end and removes by moving the last entry into the hole).
//! Rose used three structures for this: the transpose `ut` (sorted by column), the member lists `rowl` (physical order) and a stamped
//! gather (`TauRow`) to read the values in physical order. A row has no duplicate column; every value carries the bits of the matching
//! `uc` entry. The capacity of a row in this pool is Rose-internal (the reference's `row_cap` is a separate integer rule).

use crate::real::Real;

/// Where a row lives: `len` entries from `start`, with room for `cap`.
#[derive(Clone, Copy, Default)]
struct Span {
    start: u32,
    len: u32,
    cap: u32,
}

#[derive(Clone)]
pub(super) struct URows<R> {
    span: Vec<Span>,
    col: Vec<u32>,
    val: Vec<R>,
}

impl<R: Real> URows<R> {
    pub fn new() -> Self {
        Self { span: Vec::new(), col: Vec::new(), val: Vec::new() }
    }

    /// `m` empty rows; the arrays keep their allocations.
    pub fn reset(&mut self, m: usize) {
        self.span.clear();
        self.span.resize(m, Span::default());
        self.col.clear();
        self.val.clear();
    }

    #[inline(always)]
    pub fn len(&self, y: usize) -> usize {
        self.span[y].len as usize
    }

    #[inline(always)]
    pub fn is_empty(&self, y: usize) -> bool {
        self.span[y].len == 0
    }

    /// The columns of row `y`, in physical order.
    #[inline(always)]
    pub fn cols(&self, y: usize) -> &[u32] {
        let s = &self.span[y];
        &self.col[s.start as usize..(s.start + s.len) as usize]
    }

    /// The columns and values of row `y`, entry for entry.
    #[inline(always)]
    pub fn row(&self, y: usize) -> (&[u32], &[R]) {
        let s = &self.span[y];
        let r = s.start as usize..(s.start + s.len) as usize;
        (&self.col[r.clone()], &self.val[r])
    }

    /// Mutable view of row `y` (the entries only).
    pub fn row_mut(&mut self, y: usize) -> (&mut [u32], &mut [R]) {
        let s = &self.span[y];
        let r = s.start as usize..(s.start + s.len) as usize;
        (&mut self.col[r.clone()], &mut self.val[r])
    }

    /// Make room for `n` entries in total in the empty row `y` (the row must be empty).
    #[cfg(test)]
    pub fn reserve_row(&mut self, y: usize, n: usize) {
        debug_assert!(self.span[y].len == 0);
        if n > self.span[y].cap as usize {
            let start = self.col.len();
            self.col.resize(start + n, 0);
            self.val.resize(start + n, R::ZERO);
            self.span[y] = Span { start: start as u32, len: 0, cap: n as u32 };
        }
    }

    /// Room for `caps[y]` entries in every (empty) row `y`, laid out in row order in one go: the layout `reserve_row` per row gives,
    /// with one growth of the arrays instead of one per row. Rows with room 0 get none.
    pub fn reserve_all(&mut self, caps: impl Iterator<Item = usize>) {
        debug_assert!(self.col.is_empty() && self.span.iter().all(|s| s.len == 0 && s.cap == 0));
        let mut at = 0usize;
        for (s, n) in self.span.iter_mut().zip(caps) {
            if n > 0 {
                *s = Span { start: at as u32, len: 0, cap: n as u32 };
                at += n;
            }
        }
        self.col.resize(at, 0);
        self.val.resize(at, R::ZERO);
    }

    /// Append `(c, v)` to row `y`; a full row moves to the end of the arrays with room to grow.
    pub fn push(&mut self, y: usize, c: usize, v: R) {
        debug_assert!(!self.cols(y).contains(&(c as u32)), "column {c} twice in row {y}");
        let s = self.span[y];
        if s.len == s.cap {
            let cap = (s.cap + (s.cap >> 1)).max(s.len + 4);
            let start = self.col.len();
            self.col.extend_from_within(s.start as usize..(s.start + s.len) as usize);
            self.val.extend_from_within(s.start as usize..(s.start + s.len) as usize);
            self.col.resize(start + cap as usize, 0);
            self.val.resize(start + cap as usize, R::ZERO);
            self.span[y] = Span { start: start as u32, len: s.len, cap };
        }
        let s = &mut self.span[y];
        let at = (s.start + s.len) as usize;
        s.len += 1;
        self.col[at] = c as u32;
        self.val[at] = v;
    }

    /// The position of column `c` in row `y`.
    #[inline(always)]
    pub fn position(&self, y: usize, c: usize) -> Option<usize> {
        self.cols(y).iter().position(|&t| t as usize == c)
    }

    /// `Vec::swap_remove` on the entry at `p` of row `y`: the last entry fills the hole.
    pub fn swap_remove(&mut self, y: usize, p: usize) {
        let s = &mut self.span[y];
        debug_assert!(p < s.len as usize);
        let (hole, last) = ((s.start as usize) + p, (s.start + s.len - 1) as usize);
        s.len -= 1;
        self.col[hole] = self.col[last];
        self.val[hole] = self.val[last];
    }

    /// Row `y` becomes empty (the room stays).
    #[inline(always)]
    pub fn clear_row(&mut self, y: usize) {
        self.span[y].len = 0;
    }

    /// The rows as `rowl` was: the columns of each row, in physical order.
    #[cfg(any(test, feature = "factor-shadow"))]
    pub fn to_rowl(&self) -> Vec<Vec<usize>> {
        (0..self.span.len()).map(|y| self.cols(y).iter().map(|&c| c as usize).collect()).collect()
    }

    /// The rows as `ut` was: `(column, value)` sorted by column.
    #[cfg(any(test, feature = "factor-shadow"))]
    pub fn to_ut(&self) -> Vec<Vec<(usize, R)>> {
        (0..self.span.len())
            .map(|y| {
                let (c, v) = self.row(y);
                let mut e: Vec<(usize, R)> = c.iter().zip(v).map(|(&c, &v)| (c as usize, v)).collect();
                e.sort_by_key(|e| e.0);
                e
            })
            .collect()
    }

    /// Rows from the transpose `ut` (sorted by column, any values) in the physical order `rowl` gives (test builds of the old structures).
    #[cfg(test)]
    pub fn from_nested(rowl: &[Vec<usize>], ut: &[Vec<(usize, R)>]) -> Self {
        let mut r = Self::new();
        r.reset(rowl.len());
        for (y, cols) in rowl.iter().enumerate() {
            r.reserve_row(y, cols.len() + 2);
            for &c in cols {
                let v = ut[y].iter().find(|e| e.0 == c).map_or(R::ZERO, |e| e.1);
                r.push(y, c, v);
            }
        }
        r
    }

    /// The audited invariants against the numeric columns `uc`: every row holds exactly the entries of the columns that have an entry at
    /// that row (the tail columns), no column twice, and the value bits of the `uc` entry.
    #[cfg(any(test, debug_assertions, feature = "factor-shadow"))]
    pub fn check(&self, uc: &[Vec<(usize, R)>], kind: &[i64]) {
        assert_eq!(self.span.len(), uc.len());
        let mut want: Vec<Vec<(usize, R)>> = vec![Vec::new(); uc.len()];
        for (c, col) in uc.iter().enumerate() {
            for &(y, v) in col {
                assert_eq!(kind[c], -2, "U entries of non-tail column {c}");
                want[y].push((c, v));
            }
        }
        for (y, w) in want.iter().enumerate() {
            let (cols, vals) = self.row(y);
            assert_eq!(cols.len(), w.len(), "length of U row {y}");
            let mut seen = std::collections::BTreeSet::new();
            for (&c, &v) in cols.iter().zip(vals) {
                assert!(seen.insert(c), "column {c} twice in U row {y}");
                let u = w.iter().find(|e| e.0 == c as usize).unwrap_or_else(|| panic!("column {c} not in U row {y}")).1;
                assert!(super::ab::same_bits(u, v), "value of column {c} in U row {y}");
            }
        }
    }
}

#[cfg(test)]
mod reserve_tests {
    use super::*;

    /// `reserve_all` gives every row the span (start, cap) and the arrays the length that `reserve_row` row by row gives.
    #[test]
    fn reserve_all_matches_reserve_row() {
        let mut seed = 0x2b_5e_17u64;
        for round in 0..200 {
            let m = 1 + round % 37;
            let caps: Vec<usize> = (0..m)
                .map(|_| {
                    seed ^= seed << 13;
                    seed ^= seed >> 7;
                    seed ^= seed << 17;
                    if seed % 3 == 0 { 0 } else { (seed % 9) as usize }
                })
                .collect();
            let mut a = URows::<f64>::new();
            a.reset(m);
            for (y, &n) in caps.iter().enumerate() {
                a.reserve_row(y, n);
            }
            let mut b = URows::<f64>::new();
            b.reset(m);
            b.reserve_all(caps.iter().copied());
            assert_eq!(a.col.len(), b.col.len());
            assert_eq!(a.val.len(), b.val.len());
            for y in 0..m {
                let (s, t) = (a.span[y], b.span[y]);
                assert_eq!((s.start, s.len, s.cap), (t.start, t.len, t.cap), "row {y} of {caps:?}");
            }
        }
    }
}
