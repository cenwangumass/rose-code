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
    pub fn with_capacity(lines: usize, entries: usize) -> Self {
        let mut start = Vec::with_capacity(lines + 1);
        start.push(0);
        Csr { start, idx: Vec::with_capacity(entries), val: Vec::with_capacity(entries), idx_bound: 0 }
    }

    /// Append the next line.
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

/// `x[y] -= a * xv` for the entries `(idx[k], val[k])` in order, four per iteration (the updates stay sequential, so repeated
/// destinations and the rounding are those of the plain loop).
///
/// # Safety
/// Every `idx[k] as usize` is `< x.len()`.
#[inline(always)]
unsafe fn sub_column<R: Real>(x: &mut [R], idx: &[u32], val: &[R], xv: R) {
    let n = idx.len().min(val.len());
    let (idx, val) = (&idx[..n], &val[..n]);
    let (di, ri) = (idx.chunks_exact(4), idx.chunks_exact(4).remainder());
    let (dv, rv) = (val.chunks_exact(4), val.chunks_exact(4).remainder());
    for (y, a) in di.zip(dv) {
        // SAFETY: the caller guarantees every destination is in range.
        unsafe {
            *x.get_unchecked_mut(y[0] as usize) -= a[0] * xv;
            *x.get_unchecked_mut(y[1] as usize) -= a[1] * xv;
            *x.get_unchecked_mut(y[2] as usize) -= a[2] * xv;
            *x.get_unchecked_mut(y[3] as usize) -= a[3] * xv;
        }
    }
    for (&y, &a) in ri.iter().zip(rv) {
        // SAFETY: as above.
        unsafe { *x.get_unchecked_mut(y as usize) -= a * xv };
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

/// Dense forward L sweep over all steps, in ascending order, for a pool whose destinations all exceed their source (`y > v`):
/// `x[v]` is final when visited, so the `count_l` filter (`abs > TOL`) runs in the same pass. Steps with `!(abs <= TOL)` (NaN
/// included) are processed. Returns `(count, sum of the processed column lengths)`.
#[inline(never)]
pub(super) fn fused_forward<R: Real>(l: &LPool<R>, x: &mut [R]) -> (usize, usize) {
    l.check_vector(x.len());
    let (idx, val) = (&l.idx[..], &l.val[..]);
    let tol = R::c(TOL);
    let (mut count, mut len_sum) = (0usize, 0usize);
    let start = &l.start[..];
    for v in 0..l.lines() {
        // SAFETY: `v < lines <= x.len()` (checked above).
        let xv = unsafe { *x.get_unchecked(v) };
        let a = xv.abs();
        count += usize::from(a > tol);
        if a <= tol {
            continue;
        }
        // SAFETY: `v + 1 <= lines < start.len()`.
        let (s, e) = unsafe { (*start.get_unchecked(v) as usize, *start.get_unchecked(v + 1) as usize) };
        len_sum += e - s;
        // SAFETY: `s <= e <= idx.len() == val.len()` (monotone `start` ending at the entry count); destinations `< idx_bound <= x.len()`.
        unsafe { sub_column(x, idx.get_unchecked(s..e), val.get_unchecked(s..e), xv) };
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
mod plain {
    use super::*;

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
                    let (c, sum) = fused_forward(&pool, &mut xn);
                    assert_eq!(c, count);
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
                        assert_eq!(plain::fused_forward(&pool, &mut xp), fused_forward(&pool, &mut xn));
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
        let (count, sum) = fused_forward(&pool, &mut x);
        assert_eq!((count, sum), (0, 1));
        assert!(x[1].is_nan());
    }
}
