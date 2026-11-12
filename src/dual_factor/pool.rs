//! Flat compressed-row pools for the immutable L columns and their transpose: one offset array, one index array and one value array,
//! the shape of the reference's L offsets/destinations/values (`H+0x38..0x60`), instead of a vector of small vectors.

use super::TOL;
use crate::real::Real;

/// Entries `idx[start[v]..start[v + 1]]` / `val[..]` of line `v`, in insertion order.
#[derive(Clone)]
pub(super) struct Csr<R> {
    pub start: Vec<u32>,
    pub idx: Vec<u32>,
    pub val: Vec<R>,
    /// One more than the largest entry index (0 when empty); `check_vector` relies on it.
    idx_bound: u32,
}

/// Forward L columns: line `v` lists the destinations `y` and multipliers of step `v`.
pub(super) type LPool<R> = Csr<R>;

impl<R: Copy> Csr<R> {
    /// No entries in the lines `[n1, loff)`.
    pub(super) fn gap_is_empty(&self, n1: usize, loff: usize) -> bool {
        self.start[n1] == self.start[loff]
    }

    #[cfg(test)]
    pub fn with_capacity(lines: usize, entries: usize) -> Self {
        let mut start = Vec::with_capacity(lines + 1);
        start.push(0);
        Csr { start, idx: Vec::with_capacity(entries), val: Vec::with_capacity(entries), idx_bound: 0 }
    }

    /// Append the next line.
    #[cfg(test)]
    pub fn push_line(&mut self, entries: impl Iterator<Item = (usize, R)>) {
        for (y, a) in entries {
            assert!(y < u32::MAX as usize);
            self.idx_bound = self.idx_bound.max(y as u32 + 1);
            self.idx.push(y as u32);
            self.val.push(a);
        }
        assert!(self.idx.len() <= u32::MAX as usize);
        self.start.push(self.idx.len() as u32);
    }

    /// L and its transpose from the structure's flat arrays (`l_start`: steps + 1 offsets, `l_idx`: row indices, `l_val`: multipliers),
    /// destinations relabelled through `step`. Same result as the reference `build_with_transpose` over the per-step slices: one pass
    /// relabels and counts, one `extend_from_slice` copies the values, the scatter fills the transpose in ascending source step
    /// (`start` is the prefix of the counts). Only the non-empty steps are walked. `forward` is true when every destination exceeds
    /// its source step.
    pub fn from_flat(lines: usize, l_start: &[usize], l_idx: &[usize], l_val: &[R], step: &[usize]) -> (Csr<R>, Csr<R>, bool) {
        let n = l_idx.len();
        let steps = l_start.len() - 1;
        assert!(l_start[0] == 0 && l_start[steps] == n && l_val.len() == n);
        assert!(n <= u32::MAX as usize && lines < u32::MAX as usize);
        let mut idx = vec![0u32; n];
        let mut tstart = vec![0u32; lines + 1];
        let (mut forward, mut top) = (true, 0usize);
        for s in 0..steps {
            let (a, b) = (l_start[s], l_start[s + 1]);
            if a == b {
                continue;
            }
            for (o, &r) in idx[a..b].iter_mut().zip(&l_idx[a..b]) {
                let y = step[r];
                top = top.max(y);
                forward &= y > s;
                tstart[y.wrapping_add(1)] += 1;
                *o = y as u32;
            }
        }
        // every destination below `lines` (a wrapped `usize::MAX` or an out-of-range label fails here)
        assert!(n == 0 || top < lines);
        let idx_bound = if n == 0 { 0 } else { top as u32 + 1 };
        let mut l = Csr { start: l_start.iter().map(|&e| e as u32).collect(), idx, val: Vec::with_capacity(n), idx_bound };
        l.val.extend_from_slice(l_val);
        for y in 0..lines {
            tstart[y + 1] += tstart[y];
        }
        let mut next: Vec<u32> = tstart[..lines].to_vec();
        let mut tidx = vec![0u32; n];
        let mut tval: Vec<R> = Vec::with_capacity(n);
        {
            let spare = tval.spare_capacity_mut();
            assert!(spare.len() >= n && tidx.len() == n && l.idx.len() == n && l.val.len() == n);
            for v in 0..steps {
                let (a, b) = (l_start[v], l_start[v + 1]);
                for k in a..b {
                    // SAFETY: `k < n` (offsets checked against `n` above and non-decreasing); `y < lines` (checked on `top`);
                    // `next[y]` runs through the `tstart[y]..tstart[y + 1]` slots counted for `y`, all below `n`, each used once.
                    unsafe {
                        let y = *l.idx.get_unchecked(k) as usize;
                        let slot = next.get_unchecked_mut(y);
                        let p = *slot as usize;
                        *slot += 1;
                        debug_assert!(p < n);
                        *tidx.get_unchecked_mut(p) = v as u32;
                        spare.get_unchecked_mut(p).write(*l.val.get_unchecked(k));
                    }
                }
            }
        }
        // SAFETY: the counts sum to `n` and every slot `0..n` was written exactly once above.
        unsafe { tval.set_len(n) };
        let idx_bound = if n == 0 { 0 } else { steps as u32 };
        let lt = Csr { start: tstart, idx: tidx, val: tval, idx_bound };
        (l, lt, forward)
    }

    #[cfg(test)]
    /// L and its transpose from the structure's steps in one counted pass each: the copy writes the pool directly (lines of known
    /// length, destinations relabelled through `step`) and counts the destinations on the way, the scatter then fills the transpose
    /// (`start` is the prefix of the counts). Returns `(l, lt, forward)`; `forward` is true when every destination exceeds its source.
    /// Same result as `push_line` per step, then `transpose(lines)`, then the `y > v` check.
    pub fn build_with_transpose<'s>(
        lines: usize,
        steps: usize,
        entries: usize,
        step: &[usize],
        line: impl Fn(usize) -> (&'s [usize], &'s [R]),
    ) -> (Csr<R>, Csr<R>, bool)
    where
        R: 's,
    {
        let mut l = Csr::with_capacity(lines, entries);
        let mut tstart = vec![0u32; lines + 1];
        let (mut forward, mut bound) = (true, 0u32);
        for s in 0..steps {
            let (rows, vals) = line(s);
            debug_assert_eq!(rows.len(), vals.len());
            for &r in rows {
                let y = step[r];
                assert!(y < u32::MAX as usize);
                bound = bound.max(y as u32 + 1);
                forward &= y > s;
                tstart[y + 1] += 1;
                l.idx.push(y as u32);
            }
            l.val.extend_from_slice(&vals[..rows.len()]);
            assert!(l.idx.len() <= u32::MAX as usize);
            l.start.push(l.idx.len() as u32);
        }
        l.idx_bound = bound;
        for y in 0..lines {
            tstart[y + 1] += tstart[y];
        }
        let mut next: Vec<u32> = tstart[..lines].to_vec();
        let mut idx = vec![0u32; l.idx.len()];
        let mut val: Vec<R> = Vec::with_capacity(l.idx.len());
        if let Some(&first) = l.val.first() {
            val.resize(l.idx.len(), first);
        }
        for v in 0..l.lines() {
            let (s, e) = (l.start[v] as usize, l.start[v + 1] as usize);
            for k in s..e {
                let y = l.idx[k] as usize;
                let p = next[y] as usize;
                next[y] += 1;
                idx[p] = v as u32;
                val[p] = l.val[k];
            }
        }
        let idx_bound = if idx.is_empty() { 0 } else { l.lines() as u32 };
        let lt = Csr { start: tstart, idx, val, idx_bound };
        (l, lt, forward)
    }

    #[inline(always)]
    pub fn lines(&self) -> usize {
        self.start.len() - 1
    }

    #[inline(always)]
    pub fn len(&self) -> usize {
        self.idx.len()
    }

    #[inline(always)]
    pub fn col(&self, v: usize) -> (&[u32], &[R]) {
        let (s, e) = (self.start[v] as usize, self.start[v + 1] as usize);
        (&self.idx[s..e], &self.val[s..e])
    }

    #[inline(always)]
    pub fn col_len(&self, v: usize) -> usize {
        (self.start[v + 1] - self.start[v]) as usize
    }

    /// Start of the entries of line `v`.
    #[inline(always)]
    pub fn begin(&self, v: usize) -> usize {
        self.start[v] as usize
    }

    /// The entries of line `v` as `(index, value)` pairs.
    #[inline(always)]
    pub fn iter(&self, v: usize) -> impl DoubleEndedIterator<Item = (usize, R)> + '_ {
        let (i, a) = self.col(v);
        i.iter().zip(a).map(|(&y, &a)| (y as usize, a))
    }

    /// Transpose with `lines` target lines by counting sort over the source lines in ascending order: each target line lists its
    /// sources `v` ascending (a repeated destination of one source keeps its entry order).
    #[cfg(test)]
    pub fn transpose(&self, lines: usize) -> Csr<R> {
        let mut start = vec![0u32; lines + 1];
        for &y in &self.idx {
            start[y as usize + 1] += 1;
        }
        for y in 0..lines {
            start[y + 1] += start[y];
        }
        let mut next: Vec<u32> = start[..lines].to_vec();
        let mut idx = vec![0u32; self.idx.len()];
        let mut val: Vec<R> = Vec::with_capacity(self.idx.len());
        // Placeholder fill with the first value (entries exist only when there is one); every slot is overwritten below.
        if let Some(&first) = self.val.first() {
            val.resize(self.idx.len(), first);
        }
        for v in 0..self.lines() {
            let (s, e) = (self.start[v] as usize, self.start[v + 1] as usize);
            for k in s..e {
                let y = self.idx[k] as usize;
                let p = next[y] as usize;
                next[y] += 1;
                idx[p] = v as u32;
                val[p] = self.val[k];
            }
        }
        let idx_bound = if idx.is_empty() { 0 } else { self.lines() as u32 };
        Csr { start, idx, val, idx_bound }
    }
}

/// `x[y] -= a * xv` for the entries `(idx[k], val[k])` in order, four per iteration and a 0-3 entry tail (the updates stay
/// sequential, so repeated destinations and the rounding are those of the plain loop).
///
/// # Safety
/// Every `idx[k] as usize` is `< x.len()`.
#[inline(always)]
unsafe fn sub_column<R: Real>(x: &mut [R], idx: &[u32], val: &[R], xv: R) {
    let n = idx.len().min(val.len());
    let mut k = 0;
    // SAFETY: `k + 4 <= n <= idx.len(), val.len()` in the loop, `k + r < n` in the tail; destinations are in range (caller).
    unsafe {
        while k + 4 <= n {
            *x.get_unchecked_mut(*idx.get_unchecked(k) as usize) -= *val.get_unchecked(k) * xv;
            *x.get_unchecked_mut(*idx.get_unchecked(k + 1) as usize) -= *val.get_unchecked(k + 1) * xv;
            *x.get_unchecked_mut(*idx.get_unchecked(k + 2) as usize) -= *val.get_unchecked(k + 2) * xv;
            *x.get_unchecked_mut(*idx.get_unchecked(k + 3) as usize) -= *val.get_unchecked(k + 3) * xv;
            k += 4;
        }
        match n - k {
            3 => {
                *x.get_unchecked_mut(*idx.get_unchecked(k) as usize) -= *val.get_unchecked(k) * xv;
                *x.get_unchecked_mut(*idx.get_unchecked(k + 1) as usize) -= *val.get_unchecked(k + 1) * xv;
                *x.get_unchecked_mut(*idx.get_unchecked(k + 2) as usize) -= *val.get_unchecked(k + 2) * xv;
            }
            2 => {
                *x.get_unchecked_mut(*idx.get_unchecked(k) as usize) -= *val.get_unchecked(k) * xv;
                *x.get_unchecked_mut(*idx.get_unchecked(k + 1) as usize) -= *val.get_unchecked(k + 1) * xv;
            }
            1 => *x.get_unchecked_mut(*idx.get_unchecked(k) as usize) -= *val.get_unchecked(k) * xv,
            _ => {}
        }
    }
}

impl<R: Copy> Csr<R> {
    /// Per-call O(1) check that `x` can be indexed by every line and every entry index of the pool; `start` is monotone and ends at
    /// `idx.len() == val.len()` by construction (`push_line` / `transpose` are the only writers).
    #[inline(always)]
    fn check_vector(&self, n: usize) {
        assert!(n >= self.lines() && n >= self.idx_bound as usize);
        debug_assert!(self.start.windows(2).all(|w| w[0] <= w[1]));
        debug_assert!(self.start.last().map(|&e| e as usize) == Some(self.idx.len()) && self.idx.len() == self.val.len());
        debug_assert!(self.idx.iter().all(|&y| y < self.idx_bound));
    }
}

/// Dense forward L sweep in ascending order over the steps `[0,n1)` and `[loff,lines)`, for a pool whose destinations all exceed
/// their source (`y > v`): `x[v]` is final when visited, so the `count_l` filter (`abs > TOL`) runs in the same pass. The steps
/// `[n1,loff)` are the singleton pivots of the factor, whose L columns are empty (`start[n1] == start[loff]`): they process nothing,
/// so skipping them leaves `x` and the length sum as the full sweep leaves them. Steps with `!(abs <= TOL)` (NaN included) are
/// processed. Returns `(count, sum of the processed column lengths)`; `count` is the number of `abs > TOL` entries of `x[..lines]`
/// when `count_gap` (the gap's entries are final after the sweep and are scanned), and `None` otherwise (the caller does not read it).
/// A pool whose gap holds entries (a structure that breaks the factor's phase invariants, e.g. a multi-entry column flagged as a slack
/// in a synthetic test) takes the full sweep, so the result is the full sweep's for every pool; a real factor has an empty gap.
#[inline(never)]
pub(super) fn fused_forward<R: Real>(
    l: &LPool<R>,
    x: &mut [R],
    n1: usize,
    loff: usize,
    count_gap: bool,
) -> (Option<usize>, usize) {
    // The per-slot count is read only with `count_gap`: the other calls run the variant without it.
    if count_gap {
        let (c, s) = fused_forward_impl::<R, true>(l, x, n1, loff);
        (Some(c), s)
    } else {
        (None, fused_forward_impl::<R, false>(l, x, n1, loff).1)
    }
}

#[inline(always)]
fn fused_forward_impl<R: Real, const COUNT: bool>(l: &LPool<R>, x: &mut [R], n1: usize, loff: usize) -> (usize, usize) {
    l.check_vector(x.len());
    assert!(n1 <= loff && loff <= l.lines());
    let (n1, loff) = if l.gap_is_empty(n1, loff) { (n1, loff) } else { (l.lines(), l.lines()) };
    let (idx, val) = (&l.idx[..], &l.val[..]);
    let tol = R::c(TOL);
    let (mut count, mut len_sum) = (0usize, 0usize);
    let start = &l.start[..];
    // One step: `x[v]` is final when visited (destinations exceed sources); `!(abs <= tol)` (NaN included) processes the column.
    macro_rules! visit {
        ($v:expr) => {{
            let v: usize = $v;
            // SAFETY: `v < lines <= x.len()` (checked above).
            let xv = unsafe { *x.get_unchecked(v) };
            let a = xv.abs();
            if COUNT {
                count += usize::from(a > tol);
            }
            if !(a <= tol) {
                // SAFETY: `v + 1 <= lines < start.len()`.
                let (s, e) = unsafe { (*start.get_unchecked(v) as usize, *start.get_unchecked(v + 1) as usize) };
                len_sum += e - s;
                // SAFETY: `s <= e <= idx.len() == val.len()` (monotone `start` ending at the entry count); destinations
                // `< idx_bound <= x.len()`.
                unsafe { sub_column(x, idx.get_unchecked(s..e), val.get_unchecked(s..e), xv) };
            }
        }};
    }
    // Four slots per guard: when all four are `abs <= tol` at the read, none of them processes a column, so none of the four
    // values changes before its turn (only a processed column writes) and the group is skipped (it adds nothing to the count);
    // otherwise the four are visited one by one in ascending order.
    macro_rules! range {
        ($lo:expr, $hi:expr) => {{
            let (mut v, hi): (usize, usize) = ($lo, $hi);
            while v + 4 <= hi {
                // SAFETY: `v + 3 < hi <= lines <= x.len()`.
                let (b0, b1, b2, b3) = unsafe {
                    (
                        x.get_unchecked(v).abs() <= tol,
                        x.get_unchecked(v + 1).abs() <= tol,
                        x.get_unchecked(v + 2).abs() <= tol,
                        x.get_unchecked(v + 3).abs() <= tol,
                    )
                };
                // All four `abs <= tol`: each contributes 0 to the count too.
                if !(b0 & b1 & b2 & b3) {
                    visit!(v);
                    visit!(v + 1);
                    visit!(v + 2);
                    visit!(v + 3);
                }
                v += 4;
            }
            while v < hi {
                visit!(v);
                v += 1;
            }
        }};
    }
    range!(0, n1);
    range!(loff, l.lines());
    if COUNT {
        // No step writes into the gap after its own turn (destinations exceed sources), and the steps after it have sources `>= loff`.
        count += x[n1..loff].iter().filter(|v| v.abs() > tol).count();
    }
    (count, len_sum)
}

/// Dense forward L sweep without a count (`!(abs <= TOL)` steps); returns the sum of the processed column lengths.
#[inline(never)]
pub(super) fn forward_sweep<R: Real>(l: &LPool<R>, x: &mut [R]) -> usize {
    l.check_vector(x.len());
    let (idx, val) = (&l.idx[..], &l.val[..]);
    let tol = R::c(TOL);
    let mut len_sum = 0usize;
    let start = &l.start[..];
    for v in 0..l.lines() {
        // SAFETY: `v < lines <= x.len()` (checked above).
        let xv = unsafe { *x.get_unchecked(v) };
        if xv.abs() <= tol {
            continue;
        }
        // SAFETY: `v + 1 <= lines < start.len()`.
        let (s, e) = unsafe { (*start.get_unchecked(v) as usize, *start.get_unchecked(v + 1) as usize) };
        len_sum += e - s;
        // SAFETY: as in `fused_forward`.
        unsafe { sub_column(x, idx.get_unchecked(s..e), val.get_unchecked(s..e), xv) };
    }
    len_sum
}

/// Dense L^T sweep over all rows, descending, processing `abs > TOL` only (NaN is skipped); returns the sum of the processed
/// row lengths.
#[inline(never)]
pub(super) fn reverse_sweep<R: Real>(lt: &Csr<R>, x: &mut [R]) -> usize {
    lt.check_vector(x.len());
    let (idx, val) = (&lt.idx[..], &lt.val[..]);
    let tol = R::c(TOL);
    let mut len_sum = 0usize;
    let start = &lt.start[..];
    let mut i = lt.lines();
    while i != 0 {
        i -= 1;
        // SAFETY: `i < lines <= x.len()` (checked above).
        let v = unsafe { *x.get_unchecked(i) };
        if !(v.abs() > tol) {
            continue;
        }
        // SAFETY: `i + 1 <= lines < start.len()`.
        let (s, e) = unsafe { (*start.get_unchecked(i) as usize, *start.get_unchecked(i + 1) as usize) };
        len_sum += e - s;
        // SAFETY: as in `fused_forward`.
        unsafe { sub_column(x, idx.get_unchecked(s..e), val.get_unchecked(s..e), v) };
    }
    len_sum
}

/// The plain checked kernels, kept as references for the differential tests.
#[cfg(test)]
pub(super) mod plain {
    use super::*;

    /// Bitwise equality (F80 by its representation).
    pub fn same_bits<R: Real>(a: R, b: R) -> bool {
        use std::any::Any;
        let (a, b) = (&a as &dyn Any, &b as &dyn Any);
        if let Some(a) = a.downcast_ref::<f64>() {
            a.to_bits() == b.downcast_ref::<f64>().unwrap().to_bits()
        } else {
            let (a, b) = (a.downcast_ref::<crate::f80::F80>().unwrap(), b.downcast_ref::<crate::f80::F80>().unwrap());
            a.repr() == b.repr()
        }
    }

    pub fn fused_forward<R: Real>(l: &LPool<R>, x: &mut [R]) -> (usize, usize) {
        let (mut count, mut len_sum) = (0usize, 0usize);
        for v in 0..l.lines() {
            let xv = x[v];
            let a = xv.abs();
            count += usize::from(a > R::c(TOL));
            if a <= R::c(TOL) {
                continue;
            }
            let (dest, mult) = l.col(v);
            len_sum += dest.len();
            for (&y, &a) in dest.iter().zip(mult) {
                x[y as usize] -= a * xv;
            }
        }
        (count, len_sum)
    }

    pub fn forward_sweep<R: Real>(l: &LPool<R>, x: &mut [R]) -> usize {
        let mut len_sum = 0usize;
        for v in 0..l.lines() {
            let xv = x[v];
            if xv.abs() <= R::c(TOL) {
                continue;
            }
            let (dest, mult) = l.col(v);
            len_sum += dest.len();
            for (&y, &a) in dest.iter().zip(mult) {
                x[y as usize] -= a * xv;
            }
        }
        len_sum
    }

    pub fn reverse_sweep<R: Real>(lt: &Csr<R>, x: &mut [R]) -> usize {
        let mut len_sum = 0usize;
        let mut i = lt.lines();
        while i != 0 {
            i -= 1;
            let v = x[i];
            if v.abs() > R::c(TOL) {
                let (src, mult) = lt.col(i);
                len_sum += src.len();
                for (&j, &lv) in src.iter().zip(mult) {
                    x[j as usize] -= lv * v;
                }
            }
        }
        len_sum
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::f80::F80;

    struct Rng(u64);
    impl Rng {
        fn next(&mut self) -> u64 {
            self.0 ^= self.0 << 13;
            self.0 ^= self.0 >> 7;
            self.0 ^= self.0 << 17;
            self.0
        }
        fn below(&mut self, n: usize) -> usize {
            (self.next() % n as u64) as usize
        }
    }

    fn bits_eq<R: Real>(a: R, b: R) -> bool {
        use std::any::Any;
        let (a, b) = (&a as &dyn Any, &b as &dyn Any);
        if let Some(a) = a.downcast_ref::<f64>() {
            a.to_bits() == b.downcast_ref::<f64>().unwrap().to_bits()
        } else {
            let (a, b) = (a.downcast_ref::<F80>().unwrap(), b.downcast_ref::<F80>().unwrap());
            a.repr() == b.repr()
        }
    }

    type Cols<R> = Vec<Vec<(usize, R)>>;

    fn value<R: Real>(rng: &mut Rng) -> R {
        let k = rng.below(12);
        let v = match k {
            0 => 0.0,
            1 => -0.0,
            2 => 1e-14,
            3 => -1e-13,
            4 => f64::NAN,
            _ => (rng.below(2000) as f64 - 1000.0) / 64.0,
        };
        R::c(v)
    }

    /// Random L: `forward` keeps every destination above its source; otherwise any destination (repeats allowed).
    fn random_cols<R: Real>(rng: &mut Rng, m: usize, forward: bool) -> Cols<R> {
        (0..m)
            .map(|v| {
                let n = if rng.below(4) == 0 { 0 } else { rng.below(6) };
                (0..n)
                    .filter_map(|_| {
                        let y = if forward {
                            if v + 1 >= m {
                                return None;
                            }
                            v + 1 + rng.below(m - v - 1)
                        } else {
                            rng.below(m)
                        };
                        Some((y, value::<R>(rng)))
                    })
                    .collect()
            })
            .collect()
    }

    fn pool_of<R: Real>(cols: &Cols<R>) -> LPool<R> {
        let mut p = LPool::with_capacity(cols.len(), 0);
        for c in cols {
            p.push_line(c.iter().copied());
        }
        p
    }

    fn same<R: Real>(a: &[R], b: &[R]) -> bool {
        a.len() == b.len() && a.iter().zip(b).all(|(&p, &q)| bits_eq(p, q))
    }

    fn reference_transpose<R: Real>(cols: &Cols<R>) -> Cols<R> {
        let mut lt = vec![Vec::new(); cols.len()];
        for (v, c) in cols.iter().enumerate() {
            for &(y, a) in c {
                lt[y].push((v, a));
            }
        }
        lt
    }

    fn reference_dense<R: Real>(cols: &Cols<R>, x: &mut [R], e: &mut f64) -> usize {
        for v in 0..cols.len() {
            let xv = x[v];
            if xv.abs() <= R::c(TOL) {
                continue;
            }
            *e = cols[v].len() as f64 * 2.0 + *e;
            for &(y, a) in &cols[v] {
                x[y] -= a * xv;
            }
        }
        (0..cols.len()).filter(|&v| x[v].abs() > R::c(TOL)).count()
    }

    fn reference_reverse<R: Real>(lt: &Cols<R>, x: &mut [R], e: &mut f64) {
        for i in (0..lt.len()).rev() {
            let v = x[i];
            if v.abs() > R::c(TOL) {
                *e = lt[i].len() as f64 * 2.0 + *e;
                for &(j, lv) in &lt[i] {
                    x[j] -= lv * v;
                }
            }
        }
    }

    /// `from_flat` against the per-step reference over random flat structures: empty steps, fewer steps than lines, repeated
    /// destinations (`step` is any map into the lines, a permutation in a third of the rounds), NaN and signed-zero values.
    fn check_flat<R: Real>(seed: u64) {
        let mut rng = Rng(seed);
        for round in 0..600 {
            let lines = 1 + rng.below(40);
            let steps = if round % 4 == 0 { rng.below(lines + 1) } else { lines };
            let mut step: Vec<usize> = (0..lines).collect();
            if round % 3 != 0 {
                for r in 0..lines {
                    step[r] = rng.below(lines);
                }
            } else {
                for i in (1..lines).rev() {
                    step.swap(i, rng.below(i + 1));
                }
            }
            let forward_only = round % 5 == 0;
            let (mut l_start, mut l_idx, mut l_val) = (vec![0usize], Vec::new(), Vec::new());
            for s in 0..steps {
                let n = if rng.below(3) == 0 { 0 } else { rng.below(7) };
                for _ in 0..n {
                    let r = if forward_only {
                        // a row whose step is above `s` when one exists
                        match (0..lines).find(|&r| step[r] > s && rng.below(2) == 0) {
                            Some(r) => r,
                            None => continue,
                        }
                    } else {
                        rng.below(lines)
                    };
                    l_idx.push(r);
                    l_val.push(value::<R>(&mut rng));
                }
                l_start.push(l_idx.len());
            }
            let (l, lt, fw) = LPool::from_flat(lines, &l_start, &l_idx, &l_val, &step);
            let (l0, lt0, fw0) =
                LPool::build_with_transpose(lines, steps, l_idx.len(), &step, |s| (&l_idx[l_start[s]..l_start[s + 1]], &l_val[l_start[s]..l_start[s + 1]]));
            assert_eq!(fw, fw0);
            for (a, b) in [(&l, &l0), (&lt, &lt0)] {
                assert_eq!(a.start, b.start);
                assert_eq!(a.idx, b.idx);
                assert_eq!(a.idx_bound, b.idx_bound);
                assert!(same(&a.val, &b.val));
            }
        }
    }

    #[test]
    fn from_flat_matches_reference() {
        check_flat::<f64>(0x9e3779b97f4a7c15);
        check_flat::<F80>(0x1234_5678_9abc_def1);
    }

    fn check<R: Real>(seed: u64) {
        let mut rng = Rng(seed);
        for round in 0..400 {
            let m = 1 + rng.below(40);
            let forward = round % 3 != 0;
            let cols = random_cols::<R>(&mut rng, m, forward);
            let pool = pool_of(&cols);
            assert_eq!(pool.lines(), m);
            for v in 0..m {
                let (i, a) = pool.col(v);
                assert_eq!(i.len(), cols[v].len());
                assert!(i.iter().zip(&cols[v]).all(|(&y, &(z, _))| y as usize == z));
                assert!(same(a, &cols[v].iter().map(|e| e.1).collect::<Vec<_>>()));
            }
            // transpose: sources ascending, repeated destinations in entry order, values untouched
            let lt = pool.transpose(m);
            let want = reference_transpose(&cols);
            for y in 0..m {
                let (i, a) = lt.col(y);
                assert_eq!(i.len(), want[y].len());
                assert!(i.iter().zip(&want[y]).all(|(&v, &(w, _))| v as usize == w));
                assert!(same(a, &want[y].iter().map(|e| e.1).collect::<Vec<_>>()));
            }
            let is_forward = cols.iter().enumerate().all(|(v, c)| c.iter().all(|&(y, _)| y > v));
            assert_eq!(is_forward, forward || is_forward);
            // repeated calls on fresh inputs, with the same pool
            for _ in 0..4 {
                let x0: Vec<R> = (0..m).map(|_| value::<R>(&mut rng)).collect();
                let e0 = (rng.below(1000)) as f64;
                // forward sweep without count
                let (mut xr, mut er) = (x0.clone(), e0);
                reference_dense(&cols, &mut xr, &mut er);
                let mut xn = x0.clone();
                let sum = forward_sweep(&pool, &mut xn);
                assert!(same(&xr, &xn));
                assert_eq!((sum as f64 * 2.0 + e0).to_bits(), er.to_bits());
                // fused count, valid only for forward pools
                if is_forward {
                    let (mut xr, mut er) = (x0.clone(), e0);
                    let count = reference_dense(&cols, &mut xr, &mut er);
                    let mut xn = x0.clone();
                    let (c, sum) = fused_forward(&pool, &mut xn, 0, 0, true);
                    assert_eq!(c, Some(count));
                    assert!(same(&xr, &xn));
                    assert_eq!((sum as f64 * 2.0 + e0).to_bits(), er.to_bits());
                }
                // plain checked kernels: bitwise identical vector, count and length sum
                {
                    let (mut xp, mut xn) = (x0.clone(), x0.clone());
                    assert_eq!(plain::forward_sweep(&pool, &mut xp), forward_sweep(&pool, &mut xn));
                    assert!(same(&xp, &xn));
                    if is_forward {
                        let (mut xp, mut xn) = (x0.clone(), x0.clone());
                        let (cp, sp) = plain::fused_forward(&pool, &mut xp);
                        assert_eq!((Some(cp), sp), fused_forward(&pool, &mut xn, 0, 0, true));
                        assert!(same(&xp, &xn));
                    }
                    let (mut xp, mut xn) = (x0.clone(), x0.clone());
                    assert_eq!(plain::reverse_sweep(&lt, &mut xp), reverse_sweep(&lt, &mut xn));
                    assert!(same(&xp, &xn));
                    // a longer vector is allowed; the tail is untouched
                    let mut xl = x0.clone();
                    xl.push(R::c(-0.0));
                    forward_sweep(&pool, &mut xl);
                    assert!(bits_eq(xl[m], R::c(-0.0)));
                }
                // reverse sweep on the transpose
                let (mut xr, mut er) = (x0.clone(), e0);
                reference_reverse(&want, &mut xr, &mut er);
                let mut xn = x0.clone();
                let sum = reverse_sweep(&lt, &mut xn);
                assert!(same(&xr, &xn));
                assert_eq!((sum as f64 * 2.0 + e0).to_bits(), er.to_bits());
            }
        }
    }

    #[test]
    fn pools_and_dense_sweeps_match_reference_f64() {
        check::<f64>(0x9e3779b97f4a7c15);
    }

    #[test]
    fn pools_and_dense_sweeps_match_reference_f80() {
        check::<F80>(0xd1b54a32d192ed03);
    }

    #[test]
    #[should_panic]
    fn short_vector_is_rejected() {
        let cols: Cols<f64> = vec![vec![(2, 1.0)], vec![], vec![]];
        let pool = pool_of(&cols);
        let mut x = vec![1.0, 0.0];
        forward_sweep(&pool, &mut x);
    }

    #[test]
    #[should_panic]
    fn index_beyond_vector_is_rejected() {
        // Two lines (`n >= lines` holds for n = 2) but an entry index 2 >= n: only the `idx_bound` check can reject it.
        let cols: Cols<f64> = vec![vec![(2, 1.0)], vec![]];
        let pool = pool_of(&cols);
        let mut x = vec![1.0, 0.0];
        assert!(x.len() >= pool.lines() && (x.len() as u32) < pool.idx_bound);
        forward_sweep(&pool, &mut x);
    }

    #[test]
    fn fused_count_skips_nan_but_processes_it() {
        let cols: Cols<f64> = vec![vec![(1, 2.0)], vec![]];
        let pool = pool_of(&cols);
        let mut x = vec![f64::NAN, 0.0];
        let (count, sum) = fused_forward(&pool, &mut x, 0, 0, true);
        assert_eq!((count, sum), (Some(0), 1));
        assert!(x[1].is_nan());
    }

    /// The restricted sweep (`[0,n1)` and `[loff,lines)`) against the full one on forward pools whose columns in `[n1,loff)` are
    /// emptied: vector bits, count (when asked for) and length sum agree; NaN and sub-tolerance values sit in the gap as well.
    fn gap_check<R: Real>(seed: u64) {
        let mut rng = Rng(seed);
        let (mut gaps, mut counted) = (0, 0);
        for round in 0..400 {
            let m = 1 + rng.below(40);
            let mut cols = random_cols::<R>(&mut rng, m, true);
            let n1 = rng.below(m + 1);
            let loff = n1 + rng.below(m - n1 + 1);
            for c in &mut cols[n1..loff] {
                c.clear();
            }
            gaps += usize::from(loff > n1);
            let pool = pool_of(&cols);
            assert!(pool.gap_is_empty(n1, loff));
            for count_gap in [false, true] {
                let x0: Vec<R> = (0..m).map(|_| value::<R>(&mut rng)).collect();
                let (mut xp, mut xn) = (x0.clone(), x0.clone());
                let (cp, sp) = plain::fused_forward(&pool, &mut xp);
                let (cn, sn) = fused_forward(&pool, &mut xn, n1, loff, count_gap);
                assert!(same(&xp, &xn), "round {round}");
                assert_eq!(sp, sn);
                assert_eq!(cn, count_gap.then_some(cp));
                counted += usize::from(count_gap && loff > n1);
            }
        }
        assert!(gaps > 100 && counted > 100);
    }

    #[test]
    fn restricted_sweep_matches_full_sweep_f64() {
        gap_check::<f64>(0x2545f4914f6cdd1d);
    }

    #[test]
    fn restricted_sweep_matches_full_sweep_f80() {
        gap_check::<F80>(0x9e3779b97f4a7c15);
    }

    /// Columns of length 0-9 with few distinct destinations (repeats inside a column), values with NaN, -0.0 and sub-tolerance
    /// entries, runs of small slots of every length and phase: the sweeps agree with the plain kernels bit for bit.
    fn short_column_check<R: Real>(seed: u64) {
        let mut rng = Rng(seed);
        for round in 0..600 {
            let m = 1 + rng.below(30);
            let cols: Cols<R> = (0..m)
                .map(|v| {
                    let n = rng.below(10);
                    let span = (m - v - 1).min(1 + rng.below(3));
                    (0..n).filter_map(|_| (span > 0).then(|| (v + 1 + rng.below(span), value::<R>(&mut rng)))).collect()
                })
                .collect();
            let pool = pool_of(&cols);
            // mostly small slots, so that groups of four are skipped and partly processed
            let x0: Vec<R> = (0..m).map(|_| if rng.below(4) == 0 { value::<R>(&mut rng) } else { R::c(if rng.below(2) == 0 { 0.0 } else { -1e-15 }) }).collect();
            let (mut xp, mut xn) = (x0.clone(), x0.clone());
            let (cp, sp) = plain::fused_forward(&pool, &mut xp);
            for count_gap in [false, true] {
                xn.copy_from_slice(&x0);
                let (cn, sn) = fused_forward(&pool, &mut xn, 0, 0, count_gap);
                assert!(same(&xp, &xn), "round {round}");
                assert_eq!((cn, sn), (count_gap.then_some(cp), sp));
            }
            let mut xs = x0.clone();
            assert_eq!(forward_sweep(&pool, &mut xs), sp);
            assert!(same(&xp, &xs));
            // the reverse sweep (any destinations) on the transpose, whose columns are the same short lengths
            let lt = pool.transpose(m);
            let (mut xp, mut xn) = (x0.clone(), x0.clone());
            assert_eq!(plain::reverse_sweep(&lt, &mut xp), reverse_sweep(&lt, &mut xn));
            assert!(same(&xp, &xn));
        }
    }

    #[test]
    fn short_columns_and_small_slot_runs_match_plain_f64() {
        short_column_check::<f64>(0x1234_5678_9abc_def1);
    }

    #[test]
    fn short_columns_and_small_slot_runs_match_plain_f80() {
        short_column_check::<F80>(0x0fed_cba9_8765_4321);
    }
}
