//! Structural mirror of the reference sparse FTRAN: the order in which it lists the nonzeros of `B^-1 a` (`reverse-engineer/specs/ftran-order.md`).
//!
//! Positions are the basis positions after the LU renumbering (the pivot step). A fresh factor has singleton positions (kind -1),
//! a nucleus list of diagonal-only positions (phase 3 pivots, kind = index in the list) and a triangular tail (phase 4 and the
//! Markowitz nucleus, slots `i21..m`) whose columns carry U entries split into tail positions (`a`) and nucleus-list positions (`b`).
//! Column replacements follow the update rules recovered from the GDB dumps. The factor is mirrored numerically as well (L multipliers,
//! U columns, diagonal, row etas, including the reference's lane-wise eta dot products), so the solve returns the reference's values
//! bit for bit and drops the entries that cancel exactly as it does.

use crate::lu_order::Structure;
use crate::real::Real;
use std::collections::BinaryHeap;

/// Sparse U solves visit only reached rows, in ascending basis position. Dense
/// solves retain their contiguous marker scan. Duplicate touches are idempotent.
enum FtranTouched {
    Dense(Vec<bool>),
    Sparse {
        words: Vec<u64>,
        occupied: Vec<usize>,
    },
}

impl FtranTouched {
    fn new(m: usize, dense: bool, occupied: Vec<usize>) -> Self {
        if dense {
            Self::Dense(vec![false; m])
        } else {
            Self::Sparse {
                words: vec![0; m.div_ceil(64)],
                occupied,
            }
        }
    }

    /// Returns true on the first visit, as needed by the tail dependency walk.
    #[inline]
    fn insert(&mut self, row: usize) -> bool {
        match self {
            Self::Dense(flags) => !std::mem::replace(&mut flags[row], true),
            Self::Sparse { words, occupied } => {
                let block = row / 64;
                let word = &mut words[block];
                if *word == 0 {
                    occupied.push(block);
                }
                let bit = 1u64 << (row % 64);
                let fresh = *word & bit == 0;
                *word |= bit;
                fresh
            }
        }
    }

    /// Consume a row in U's requested output order. Every numeric write has been
    /// marked by this point, so duplicate and zero outputs can consume the same bit.
    fn take_sparse(&mut self, row: usize) -> bool {
        let Self::Sparse { words, .. } = self else {
            unreachable!()
        };
        let bit = 1u64 << (row % 64);
        let word = &mut words[row / 64];
        let was_set = *word & bit != 0;
        *word &= !bit;
        was_set
    }

    fn finish(&mut self) {
        if let Self::Sparse { occupied, .. } = self {
            occupied.sort_unstable();
        }
    }

    #[inline]
    fn for_each(&self, mut visit: impl FnMut(usize)) {
        match self {
            Self::Dense(flags) => {
                for (row, &flag) in flags.iter().enumerate() {
                    if flag {
                        visit(row);
                    }
                }
            }
            Self::Sparse { words, occupied } => {
                for &block in occupied {
                    let mut bits = words[block];
                    while bits != 0 {
                        visit(block * 64 + bits.trailing_zeros() as usize);
                        bits &= bits - 1;
                    }
                }
            }
        }
    }
}

/// Per-factor scratch and an index of the append-only row etas. A cloned
/// factor builds its own scratch on demand; no caller can observe the cache.
#[derive(Default)]
#[cfg_attr(test, derive(Debug))]
struct LevelCache<R>(std::cell::RefCell<Option<Box<LevelWorkspace<R>>>>);

impl<R: Default> Clone for LevelCache<R> {
    fn clone(&self) -> Self {
        Self::default()
    }
}

#[derive(Default)]
#[cfg_attr(test, derive(Debug))]
struct LevelWorkspace<R> {
    values: Vec<R>,
    // 0 = untouched, 1 = queued for L, 2 = member of the sparse result.
    state: Vec<u8>,
    touched: Vec<usize>,
    buckets: Vec<Vec<usize>>,
    order: Vec<usize>,
    // For a newly active input, schedule only later etas that read it.
    uses: Vec<Vec<usize>>,
    term_ends: Vec<usize>,
    pending: Vec<u64>,
    gather: Vec<R>,
    incompatible: bool,
}

impl<R: Real> LevelWorkspace<R> {
    fn prepare(&mut self, f: &RefStruct<R>, nlevels: usize) -> bool {
        if self.incompatible {
            return false;
        }
        if self.values.is_empty() {
            self.values.resize(f.m, R::ZERO);
            self.state.resize(f.m, 0);
            self.buckets.resize_with(nlevels, Vec::new);
        }
        if self.term_ends.len() < f.etas.len() {
            self.uses.resize_with(f.m, Vec::new);
            for e in self.term_ends.len()..f.etas.len() {
                let eta = &f.etas[e];
                // The omitted eta would receive only positive zeros. Its
                // result rounds to canonical zero even with signed coefficients.
                if eta.ind.first() != Some(&eta.idx)
                    || !eta.val.iter().all(|v| v.to_f64().is_finite())
                {
                    self.incompatible = true;
                    return false;
                }
                let end = self.term_ends.last().copied().unwrap_or(0) + eta.ind.len();
                // The reference adds nonnegative integer WORK terms to a double.
                // This bound makes both its full sum and every partial sum exact.
                if end as u64 > (1u64 << 50) {
                    self.incompatible = true;
                    return false;
                }
                for &row in &eta.ind {
                    let uses = &mut self.uses[row];
                    if uses.is_empty() && f.first_eta[row] != e {
                        self.incompatible = true;
                        return false;
                    }
                    if uses.last() != Some(&e) {
                        uses.push(e);
                    }
                }
                self.term_ends.push(end);
            }
            self.pending.resize(f.etas.len().div_ceil(64), 0);
        }

        for row in self.touched.drain(..) {
            self.values[row] = R::ZERO;
            self.state[row] = 0;
        }
        self.order.clear();
        true
    }

    fn queue(&mut self, row: usize, level: usize) {
        if self.state[row] == 0 {
            self.state[row] = 1;
            self.touched.push(row);
            self.buckets[level].push(row);
        }
    }

    fn schedule(&mut self, row: usize, after: usize) {
        let uses = &self.uses[row];
        let start = uses.partition_point(|&index| index < after);
        for &e in &uses[start..] {
            self.pending[e / 64] |= 1u64 << (e % 64);
        }
    }

    fn solve(
        &mut self,
        f: &RefStruct<R>,
        levels: &[usize],
        maxlev: usize,
        rows: &[usize],
        vals: &[R],
    ) -> Option<(Vec<(usize, R)>, usize, f64)> {
        if !self.prepare(f, maxlev + 1) {
            return None;
        }
        let mut low = maxlev + 1;
        let mut high = 0;
        for (&r, &v) in rows.iter().zip(vals) {
            let i = f.perm[r];
            self.values[i] += v;
            self.queue(i, levels[i]);
            low = low.min(levels[i]);
            high = high.max(levels[i]);
        }
        let mut at_tolerance = false;
        while low <= high {
            let mut j = 0;
            while j < self.buckets[low].len() {
                let row = self.buckets[low][j];
                j += 1;
                let x = self.values[row];
                if x.abs() >= R::c(TOL) {
                    self.order.push(row);
                    self.state[row] = 2;
                    at_tolerance |= x.abs() == R::c(TOL);
                    for &(to, a) in f.lcolumn(row) {
                        self.queue(to, levels[to]);
                        high = high.max(levels[to]);
                        self.values[to] -= a * x;
                    }
                } else {
                    self.values[row] = R::ZERO;
                }
            }
            self.buckets[low].clear();
            low += 1;
        }
        let count = self.order.len();
        // The accounting suffix and numerical prefix may start with different
        // supports on a tolerance tie. Dense handoffs also retain raw residues.
        if !f.etas.is_empty() && (at_tolerance || count * f.dthr > f.m) {
            return None;
        }
        let mut work =
            ((f.nnz_l + f.m + 2 * f.m) as f64 * 5.0 * count as f64) / (f.m + 1) as f64 + 0.0;
        if crate::dump::on("lstage") {
            f.lsnap.borrow_mut().0.clone_from(&self.values);
        }
        if !f.etas.is_empty() {
            let mut first = f.etas.len();
            for j in 0..count {
                let row = self.order[j];
                first = first.min(f.first_eta[row]);
                self.schedule(row, 0);
            }
            for word in first / 64..self.pending.len() {
                while self.pending[word] != 0 {
                    let bit = self.pending[word].trailing_zeros() as usize;
                    self.pending[word] &= self.pending[word] - 1;
                    let e = word * 64 + bit;
                    let eta = &f.etas[e];
                    self.gather.clear();
                    self.gather
                        .extend(eta.ind.iter().map(|&row| self.values[row]));
                    let v = dot(&self.gather, &eta.val, eta.peel);
                    let row = eta.idx;
                    if v.abs() <= R::c(TOL) {
                        self.values[row] = R::ZERO;
                    } else {
                        self.values[row] = v;
                        if self.state[row] != 2 {
                            if self.state[row] == 0 {
                                self.touched.push(row);
                            }
                            self.state[row] = 2;
                            self.order.push(row);
                            self.schedule(row, e + 1);
                        }
                    }
                }
            }
            let prefix = first.checked_sub(1).map_or(0, |i| self.term_ends[i]);
            let terms = self.term_ends.last().unwrap() - prefix;
            let ew = 2.0 * terms as f64 + (f.etas.len() - first) as f64 * 3.0;
            work = ew + 3.0 * self.order.len() as f64 + work;
        }
        if crate::dump::on("lstage") {
            f.lsnap.borrow_mut().1.clone_from(&self.values);
        }
        let out = self
            .order
            .iter()
            .copied()
            .filter(|&i| self.values[i].abs() > R::c(TOL))
            .map(|i| (i, self.values[i]))
            .collect();
        Some((out, count, work))
    }
}

/// Scratch for independent BTRAN calls. An empty clone avoids copying dead work vectors.
#[derive(Default)]
#[cfg_attr(test, derive(Debug))]
struct BtranWork<R> {
    x: Vec<R>,
    tail: SparseBtranTail<R>,
    input_order: Vec<usize>,
    output_order: Vec<usize>,
    marks: Vec<u32>,
    generation: u32,
    eta_heads: Vec<usize>,
    eta_prev: Vec<usize>,
    pending: BinaryHeap<usize>,
}

impl<R: Default> Clone for BtranWork<R> {
    fn clone(&self) -> Self {
        Self::default()
    }
}

impl<R> BtranWork<R> {
    fn next_generation(&mut self, m: usize) -> u32 {
        self.marks.resize(m, 0);
        if self.generation == u32::MAX {
            self.marks.fill(0);
            self.generation = 0;
        }
        self.generation += 1;
        self.generation
    }

    fn index_etas(&mut self, etas: &[Eta<R>], m: usize) {
        self.eta_heads.resize(m, usize::MAX);
        for e in self.eta_prev.len()..etas.len() {
            self.eta_prev
                .push(std::mem::replace(&mut self.eta_heads[etas[e].idx], e));
        }
    }

    fn schedule_row(&mut self, row: usize, before: usize) {
        let mut e = self.eta_heads[row];
        while e != usize::MAX && e >= before {
            e = self.eta_prev[e];
        }
        if e != usize::MAX {
            self.pending.push(e);
        }
    }
}

#[cfg(test)]
#[path = "btran_oracle.rs"]
mod btran_oracle;

/// Working values are indexed by permutation slot. A clear occupancy bit represents +0,
/// so writes from previous right-hand sides need no numeric-array clear.
#[derive(Default)]
#[cfg_attr(test, derive(Debug))]
struct SparseBtranTail<R> {
    values: Vec<R>,
    occupied: Vec<u64>,
}

trait BtranTail<R: Real> {
    fn entry(&mut self, p: usize, slots: &[usize]) -> &mut R;

    fn next_slot(&self, from: usize, end: usize) -> Option<usize>;

    fn take(&mut self, p: usize, slot: usize) -> R;
}

impl<R: Real> BtranTail<R> for Vec<R> {
    fn entry(&mut self, p: usize, _: &[usize]) -> &mut R {
        &mut self[p]
    }

    fn next_slot(&self, from: usize, end: usize) -> Option<usize> {
        (from < end).then_some(from)
    }

    fn take(&mut self, p: usize, _: usize) -> R {
        std::mem::replace(&mut self[p], R::ZERO)
    }
}

impl<R: Real> BtranTail<R> for SparseBtranTail<R> {
    fn entry(&mut self, p: usize, slots: &[usize]) -> &mut R {
        let slot = slots[p];
        let word = &mut self.occupied[slot / 64];
        let bit = 1u64 << (slot % 64);
        if *word & bit == 0 {
            self.values[slot] = R::ZERO;
            *word |= bit;
        }
        &mut self.values[slot]
    }

    fn next_slot(&self, from: usize, _: usize) -> Option<usize> {
        let (word, skip) = (from / 64, from % 64);
        let words = self.occupied.get(word..)?;
        let mut mask = u64::MAX << skip;
        for (i, &bits) in words.iter().enumerate() {
            let pending = bits & mask;
            if pending != 0 {
                return Some((word + i) * 64 + pending.trailing_zeros() as usize);
            }
            mask = u64::MAX;
        }
        None
    }

    fn take(&mut self, _: usize, slot: usize) -> R {
        self.occupied[slot / 64] &= !(1 << (slot % 64));
        self.values[slot]
    }
}

impl<R: Real> SparseBtranTail<R> {
    fn reset(&mut self, m: usize) {
        self.values.resize(m, R::ZERO);
        self.occupied.resize(m.div_ceil(64), 0);
        self.occupied.fill(0);
    }
}

#[cfg(test)]
mod sparse_btran_tail_tests {
    use super::*;

    fn identity(m: usize) -> RefStruct {
        let cols: Vec<_> = (0..m).map(|i| (true, vec![i], vec![1.0])).collect();
        let st = crate::lu_order::factor_structure(m, &cols);
        RefStruct::fresh(m, m, &st, &cols)
    }

    fn check(f: &RefStruct, w: &mut SparseBtranTail<f64>, inp: &[(usize, f64)]) {
        let mut old_x = vec![0.0; f.m];
        let mut new_x = vec![0.0; f.m];
        let (old_order, old_e) = f.sparse_btran_u_scan(inp, &mut old_x);
        w.reset(f.m);
        let (new_order, new_e) = f.sparse_btran_u(inp, &mut new_x, w, Vec::new());
        let mut dense_x = vec![0.0; f.m];
        let (dense_order, dense_e) =
            f.sparse_btran_u(inp, &mut dense_x, &mut vec![0.0; f.m], Vec::new());
        assert_eq!(dense_order, old_order);
        assert_eq!(dense_e.to_bits(), old_e.to_bits());
        assert!(dense_x
            .iter()
            .zip(&old_x)
            .all(|(a, b)| a.to_bits() == b.to_bits()));
        assert_eq!(new_order, old_order);
        assert_eq!(new_e.to_bits(), old_e.to_bits());
        assert!(new_x
            .iter()
            .zip(old_x)
            .all(|(a, b)| a.to_bits() == b.to_bits()));
    }

    #[test]
    fn occupied_slots_match_scan_across_permutations_and_reused_right_hand_sides() {
        let mut seed = 0x3a6d4315_u64;
        let mut next = || {
            seed ^= seed << 13;
            seed ^= seed >> 7;
            seed ^= seed << 17;
            seed as usize
        };
        let vals = [
            1.0,
            -1.0,
            0.5,
            0.0,
            -0.0,
            TOL,
            -TOL,
            f64::from_bits(1),
            f64::MAX,
            f64::INFINITY,
            f64::NEG_INFINITY,
            f64::from_bits(0x7ff8000000000123),
            f64::from_bits(0x7ff8000000000456),
        ];
        for m in [1, 63, 64, 65, 127, 128, 129, 257] {
            let mut f = identity(m);
            let mut w = SparseBtranTail::default();
            for pass in 0..8 {
                for i in (1..m).rev() {
                    f.at.swap(i, next() % (i + 1));
                }
                for (s, &p) in f.at.iter().enumerate() {
                    f.slot[p] = s;
                }
                f.i21 = next() % m;
                f.h33 = if pass % 2 == 0 { 5.0 } else { 5.25 };
                let nval = if pass < 4 { 7 } else { vals.len() };
                for p in 0..m {
                    f.kind[p] = if f.slot[p] < f.i21 {
                        (p % 2) as i64 - 1
                    } else {
                        -2
                    };
                    f.diag[p] = [1.0, -2.0, f64::MIN_POSITIVE, 0.5][next() % 4];
                    let r = next() % m;
                    // Self, backward, repeated and forward writes occur in the same scatter.
                    f.ut[p] = [r, p, (p + m - 1) % m, r, (p + 1) % m]
                        .into_iter()
                        .map(|r| (r, vals[next() % nval]))
                        .collect();
                }
                for n in [0, 1, 2, 3, m / 2, m + 1] {
                    let inp: Vec<_> = (0..n).map(|_| (next() % m, vals[next() % nval])).collect();
                    check(&f, &mut w, &inp);
                }
            }
        }
    }

    #[test]
    fn repeated_large_btran_preserves_results_and_observable_state() {
        let mut f = identity(4096);
        assert!(f.bt_occupied);
        f.i21 = 8;
        for p in f.i21..f.m {
            f.kind[p] = -2;
        }
        f.ut[2] = vec![(63, 0.25), (127, -1.0), (63, 0.5)];
        f.ut[63] = vec![(64, 1.0), (4051, -0.25)];
        f.ut[64] = vec![(63, -1.0), (4095, 0.5), (64, 0.125)];
        let mut scan = f.clone();
        scan.bt_occupied = false;
        let rhs = [
            vec![(2, 1.0)],
            vec![(64, -2.0)],
            vec![],
            vec![(2, 1.0), (2, -1.0)],
        ];
        for dense_in in [false, true] {
            for input in &rhs {
                let bits = |v: Vec<(usize, f64)>| {
                    v.into_iter()
                        .map(|(p, x)| (p, x.to_bits()))
                        .collect::<Vec<_>>()
                };
                let expected = bits(scan.btran(input, true, dense_in));
                assert_eq!(bits(f.btran(input, true, dense_in)), expected);
                assert_eq!(f.work.to_bits(), scan.work.to_bits());
                assert_eq!(f.bt_dense, scan.bt_dense);
                assert_eq!(f.row_cached, scan.row_cached);
            }
        }
        let clone = f.clone();
        assert!(clone.bt_work.as_ref().unwrap().tail.values.is_empty());
    }

    #[test]
    fn captured_sc50a_transpose_matches_original_tail() {
        let records: serde_json::Value =
            serde_json::from_str(include_str!("../tests/fixtures/ftran_dense_tail.json")).unwrap();
        for rec in records.as_array().unwrap() {
            let diag: Vec<f64> = serde_json::from_value(rec["diag"].clone()).unwrap();
            let mut f = identity(diag.len());
            f.diag = diag;
            f.kind = serde_json::from_value(rec["kind"].clone()).unwrap();
            f.slot = serde_json::from_value(rec["slot"].clone()).unwrap();
            f.i21 = f.kind.iter().filter(|&&kind| kind != -2).count();
            for (p, &s) in f.slot.iter().enumerate() {
                f.at[s] = p;
            }
            let cols: Vec<Vec<(usize, f64)>> = serde_json::from_value(rec["uc"].clone()).unwrap();
            for (c, col) in cols.iter().enumerate() {
                for &(p, v) in col {
                    f.ut[p].push((c, v));
                }
            }
            let mut w = SparseBtranTail::default();
            for h in [2.0, 7.0] {
                f.h33 = h;
                for p in 0..f.m {
                    check(&f, &mut w, &[(p, 1.0), ((p + 1) % f.m, -0.25), (p, -0.5)]);
                }
            }
        }
    }
}

/// Direct numeric row lookup for tau dots. Stamps retain the first duplicate without clearing per row.
#[derive(Clone, Default)]
#[cfg_attr(test, derive(Debug))]
struct TauRow<R> {
    values: Vec<R>,
    stamps: Vec<usize>,
    epoch: usize,
}

impl<R: Real> TauRow<R> {
    fn prepare(&mut self, m: usize, row: &[(usize, R)]) {
        if self.values.len() != m {
            self.values.resize(m, R::ZERO);
            self.stamps.resize(m, 0);
        }
        self.epoch = self.epoch.wrapping_add(1);
        if self.epoch == 0 {
            self.stamps.fill(0);
            self.epoch = 1;
        }
        for &(c, v) in row {
            if self.stamps[c] != self.epoch {
                self.stamps[c] = self.epoch;
                self.values[c] = v;
            }
        }
    }

    fn get(&self, c: usize) -> Option<R> {
        (self.stamps[c] == self.epoch).then(|| self.values[c])
    }

    /// Append the row's surviving terms in virtual rowl order, including its duplicates.
    #[inline(never)]
    fn append_terms(&self, cols: &[usize], x: &[R], out: &mut Vec<(R, R)>) {
        let mut count = out.len();
        out.reserve(cols.len());
        let dst = out.as_mut_ptr();
        for &c in cols {
            if let Some(v) = self.get(c) {
                // Before each write, count <= original_len + visited < original_len + cols.len() <= capacity.
                // Reads stay checked, and each write initializes the next appended entry.
                unsafe {
                    dst.add(count).write((v, x[c]));
                }
                count += 1;
            }
        }
        // Appended entries are initialized in the same filtered order as the old extend.
        unsafe {
            out.set_len(count);
        }
    }
}

#[cfg(test)]
mod tau_gather_tests {
    use super::*;

    fn check<R: Real>(bits: impl Fn(R) -> (u64, u16)) {
        let values = [
            0.0,
            -0.0,
            f64::from_bits(1),
            -f64::MIN_POSITIVE,
            f64::MAX,
            -f64::MAX,
            f64::INFINITY,
            f64::NEG_INFINITY,
            f64::from_bits(0x7ff8_0000_0000_0001),
            f64::from_bits(0xfff8_0000_0000_0007),
        ];
        for m in [0, 1, 7, 31, 128] {
            let x: Vec<_> = (0..m)
                .map(|i| R::c(values[(i + 3) % values.len()]))
                .collect();
            let mut entries = Vec::new();
            for i in 0..m {
                if i % 3 != 2 {
                    entries.push((i, R::c(values[i % values.len()])));
                    entries.push((i, R::c(19.0)));
                }
            }
            let order: Vec<_> = if m == 0 {
                Vec::new()
            } else {
                (0..m * 2 + 7).rev().map(|i| (i * 7) % m).collect()
            };
            let mut row = TauRow::<R>::default();
            row.prepare(m, &entries);
            let prefix = (R::c(8.0), R::c(-13.0));
            let mut expected = vec![prefix];
            for &c in &order {
                if let Some(&(_, v)) = entries.iter().find(|e| e.0 == c) {
                    expected.push((v, x[c]));
                }
            }
            let mut actual = vec![prefix];
            row.append_terms(&order, &x, &mut actual);
            assert_eq!(
                actual
                    .iter()
                    .map(|&(v, x)| (bits(v), bits(x)))
                    .collect::<Vec<_>>(),
                expected
                    .iter()
                    .map(|&(v, x)| (bits(v), bits(x)))
                    .collect::<Vec<_>>()
            );
        }
        let mut row = TauRow::<R>::default();
        row.prepare(5, &[]);
        let mut out = Vec::new();
        row.append_terms(&[4, 3], &[], &mut out);
        assert!(out.is_empty());
    }

    #[test]
    fn ordered_gather_keeps_first_duplicates_holes_and_bits() {
        check::<f64>(|v| (v.to_bits(), 0));
        check::<crate::f80::F80>(|v| v.repr());
    }
}

#[derive(Clone)]
#[cfg_attr(test, derive(Debug))]
pub struct RefStruct<R: Real = f64> {
    m: usize,
    /// The pivot-row BTRAN of this position ran since the last update: the reference then finds the row-eta solve cached (`F+0x1d8`, `0xf193e`).
    row_cached: Option<usize>,
    kind: Vec<i64>,
    nuc: Vec<usize>,
    i21: usize,
    /// The reference's upper bound `[L+0xe8]` on the position of a tail column without entries, scanned downward at updates.
    mark: usize,
    a: Vec<Vec<usize>>,
    b: Vec<Vec<usize>>,
    /// Immutable forward L columns, preserving coefficient and destination order.
    lstart: Vec<usize>,
    lentries: Vec<(usize, R)>,
    perm: Vec<usize>,
    /// Immutable inverse row permutation, initialized at the first BTRAN output.
    row_of: std::cell::OnceCell<Vec<usize>>,
    dens: f64,
    dense_est: std::cell::Cell<f64>,
    /// Density threshold fixed at factorization time (`f104` in the notes): 2, or `5 (nnzL/m + 1)^2` capped at `m`.
    dthr: usize,
    /// Columns eligible for demotion, keyed by identity so slot walks need no index repairs.
    empty_tail: std::collections::BTreeSet<usize>,
    slot: Vec<usize>,
    at: Vec<usize>,
    /// Reused counts for ordering dense lists of permutation slots.
    slot_counts: std::cell::RefCell<Vec<usize>>,
    etas: Vec<Eta<R>>,
    bt_occupied: bool,
    bt_work: Option<Box<BtranWork<R>>>,
    /// Immutable L transpose, in ascending source-column order, shared by all BTRAN calls of this factor.
    lt: Vec<Vec<(usize, R)>>,
    /// Immutable forward levels and their maximum for level-bucket L solves (`0x12a990`). Used only
    /// for `m >= 250`, a small level count (`< m/10`) and a short Markowitz L part.
    levels: Option<(Vec<usize>, usize)>,
    level_cache: LevelCache<R>,
    /// Immutable transpose-graph levels and maximum, initialized when a BTRAN first uses level buckets.
    lt_levels: std::cell::OnceCell<(Vec<usize>, usize)>,
    lev_built: bool,
    diag: Vec<R>,
    /// Numeric U entries of each tail column `(position, value)`.
    uc: Vec<Vec<(usize, R)>>,
    /// Numeric U transpose, ordered by column, including duplicate entries in column order.
    ut: Vec<Vec<(usize, R)>>,
    /// Spike of the most recent solve: the L stage and eta result `(position, value)` in list order.
    spike: Vec<(usize, R)>,
    /// Whether the column replacement has used the spike since the last FTRAN.
    spike_used: bool,
    /// Reused by the two disjoint numeric phases of a column replacement.
    update_values: Vec<R>,
    update_seen: Vec<bool>,
    /// Dense x of the most recent solve after the L stage and after the row etas, for the `lstage` dump.
    lsnap: std::cell::RefCell<(Vec<R>, Vec<R>)>,
    /// Refactorization trigger inputs (`0x142350`): updates since the factorization, accumulated solve result counts, accumulated
    /// eta lengths, and the factor size `nnz(L) + nnz(U) + 2 m` fixed at factorization time.
    upd: usize,
    sum_alpha: f64,
    /// `upd` and `sum_alpha` before the most recent update: the scan window of `pricing0` counts the updates made before the one
    /// whose pivot triggered the refactorization.
    before_last: (usize, f64),
    sum_eta: f64,
    last_alpha: usize,
    znz: f64,
    /// Work-counter inputs (`crate::work`): phase-1 steps (`L+0xb8`), singleton-phase steps (`L+0xc4`), `nnz(L)`, and for each
    /// position the first row eta that lists it (`F+0x1a0`).
    n1: usize,
    loff: usize,
    nnz_l: usize,
    /// Nonzeros in the immutable L tail starting at loff.
    nnz_l_tail: usize,
    first_eta: Vec<usize>,
    /// Estimate of the most recent `ftran` (`0x120740`) or `btran` (`0x12de60`), before the factor object's scale 1.3.
    pub work: f64,
    /// The L-stage part and the U-stage part of `work` of the last `ftran`, and the L-stage estimate of the last `ftran_needed`.
    pub work_parts: (f64, f64),
    /// The work estimate of the last `ftran_keep`, which restores `work` itself.
    pub keep_work: f64,
    /// Whether the result of the last `ftran_keep` was dense (the reference's count -1).
    pub keep_dense: bool,
    pub tau_l: std::cell::Cell<f64>,
    /// The U-stage estimate of the last `ftran_needed` (`0x123690`).
    pub tau_u: std::cell::Cell<f64>,
    /// What the tau solve's U stage returns (`0x123690`, `0x1244df`): the number of nonzero tail columns, or the tail rows `m - i21` when it runs dense; the
    /// two-vector FTRAN adds it to the first result's count (`0x12963d`).
    pub tau_fin: std::cell::Cell<usize>,
    /// What the last `ftran` and `ftran_needed` leave for the joint U stages of the two-vector FTRAN (`joint_u`).
    f2: std::cell::RefCell<F2<R>>,
    /// Scratch terms of the second-vector nucleus row dot; keep capacity across rows and solves.
    tau_terms: std::cell::RefCell<Vec<(R, R)>>,
    /// Reusable direct view of one numeric U row; its values are read in virtual rowl order.
    tau_row: std::cell::RefCell<TauRow<R>>,
    /// Fill of the reference's pool for etas and the row-wise copy of U (`F+0x160`), and the capacity of each row of that copy.
    /// Only the parity of `fill` matters: it decides the alignment of the next eta's value array.
    fill: usize,
    row_cap: Vec<usize>,
    /// Parity of the start of each row in the row-wise copy (the value array of a row starts on a 16-byte boundary when it is 0): the rows of a
    /// fresh factor start on multiples of 4, a row moved to the end of the pool starts at the fill of that moment.
    row_par: Vec<u8>,
    /// The reference's row-wise copy: for each row, the tail columns whose A or B list holds it, in insertion order.
    rowl: Vec<Vec<usize>>,
    /// Alignment of the scratch eta behind each new diagonal (`ctx+0x1d0`): its value array sits `4 m` bytes after the index array in
    /// an arena whose offsets add `24 maxupd + 8 (m + ncols) + 8`, so the parity is `(maxupd + ncols + 1) % 2` whatever `m` is.
    scratch_peel: usize,
    /// Work of the last column replacement (`0xf15c0`), as `(site, estimate, scale)`.
    pub upd_work: Vec<(u32, f64, f64)>,
    /// Tail and nucleus U densities fixed at factorization time (`F+0x108`, `F+0x110`), and whether the level structure was built
    /// but has too many levels (`F+0x140 >= m/10`), which makes the L^T solve dense.
    h33: f64,
    h34: f64,
    lev_dense: bool,
    /// The factor's longest L path is at least `0.1 m` (`m >= 250`), whatever the strategy gate: the reference's L stage then
    /// hands on a dense vector (count -1) and the dense U routine runs.
    lev_dense_ft: bool,
    last_dense: bool,
    /// Whether the last result came out in ascending position order (the dense scan or the dense U routine): the reference's
    /// pivot column is then a dense vector for the Devex `wq` sum.
    last_ordered: bool,
    /// Whether the most recent `btran` returned a dense result (count `-1`): its L^T solve ran dense.
    pub bt_dense: bool,
    /// The factor object's update limit when the pricing lowered it (`F+0x188`, the dual path); the default is `max(20, m^0.7)`.
    pub upd_limit: Option<usize>,
}

/// The inputs of the joint U stages of the two-vector FTRAN: for the pivot column and for tau, the count of the L-stage list (`None`: dense, count -1)
/// and the solved values (the tail entries are the values after the tail solve); for the column also the positions of the L-stage list.
#[derive(Clone)]
#[cfg_attr(test, derive(Debug))]
struct F2<R> {
    col_n: Option<usize>,
    col_in: Vec<usize>,
    x_col: Vec<R>,
    /// All 64-row blocks written by the preceding sparse U solve. This list
    /// travels with its snapshot so a clone or replaced factor resets its own data.
    dirty_blocks: Option<Vec<usize>>,
    tau_n: Option<usize>,
    x_tau: Vec<R>,
    /// Whether the raw value of tau at a tail column exceeded the tolerance before the division (the joint routines test it first).
    tau_raw: Vec<bool>,
}

impl<R> Default for F2<R> {
    fn default() -> Self {
        F2 {
            col_n: None,
            col_in: Vec::new(),
            x_col: Vec::new(),
            dirty_blocks: None,
            tau_n: None,
            x_tau: Vec::new(),
            tau_raw: Vec::new(),
        }
    }
}

#[derive(Clone)]
#[cfg_attr(test, derive(Debug))]
struct Eta<R> {
    idx: usize,
    ind: Vec<usize>,
    val: Vec<R>,
    /// Elements before the first 16-byte boundary of the value array in the reference's pool (the parity of the pool fill at which
    /// the eta was stored).
    peel: usize,
}

const TOL: f64 = 1e-13;

/// Dump helpers: the values as doubles.
fn jv<R: Real>(v: &[R]) -> Vec<f64> {
    v.iter().map(|x| x.to_f64()).collect()
}

fn ju<R: Real>(v: &[Vec<(usize, R)>]) -> Vec<Vec<(usize, f64)>> {
    v.iter()
        .map(|x| x.iter().map(|e| (e.0, e.1.to_f64())).collect())
        .collect()
}

/// The row-wise dot product of the second-vector U stage (`0x123690`, loop `0x1241ed`): `s` minus the products `(value, x)`. Below 8 terms, or
/// when fewer than `peel + 8` are left, it is sequential; otherwise `peel` (0 or 1: the elements before the first 16-byte boundary of the
/// value array) terms are subtracted one by one, then blocks of eight go to four two-lane accumulators (the first one starts with `s` in lane
/// 0), combined as `(a0 + a1) + (a2 + a3)` lane-wise and then lane 0 + lane 1, and the last `(n - peel) % 8` terms are subtracted one by one.
fn row_dot<R: Real>(s: R, t: &[(R, R)], peel: usize) -> R {
    let n = t.len();
    let prod = |i: usize| t[i].0 * t[i].1;
    if !R::LANES || n < 8 || n < peel + 8 {
        return (0..n).fold(s, |a, i| a - prod(i));
    }
    let mut s = s;
    for i in 0..peel {
        s -= prod(i);
    }
    let end = n - ((n - peel) & 7);
    let mut acc = [
        [s, R::ZERO],
        [R::ZERO, R::ZERO],
        [R::ZERO, R::ZERO],
        [R::ZERO, R::ZERO],
    ];
    let mut i = peel;
    while i < end {
        for (q, a) in acc.iter_mut().enumerate() {
            a[0] -= prod(i + 2 * q);
            a[1] -= prod(i + 2 * q + 1);
        }
        i += 8;
    }
    let l = [acc[0][0] + acc[1][0], acc[0][1] + acc[1][1]];
    let h = [acc[2][0] + acc[3][0], acc[2][1] + acc[3][1]];
    let mut r = (l[0] + h[0]) + (l[1] + h[1]);
    for i in end..n {
        r -= prod(i);
    }
    r
}

/// The row dot of the extended-precision two-vector FTRAN output stage (`0x2642b0`, loop `0x2678c4`): for `n >= 8` terms, eight accumulators by term
/// position mod 8 (the first starting with `s`, the others at zero, each subtracting its terms), added in order `((((((a0 + a1) + a2) + ...) + a7`; the
/// remaining `n % 8` terms (all of them for `n < 8`) are then subtracted from the last to the first (the Duff cascade at `0x267b0d..0x267bbd`).
fn row_dot_q<R: Real>(s: R, t: &[(R, R)]) -> R {
    let n = t.len();
    let nb = n / 8 * 8;
    let prod = |i: usize| t[i].0 * t[i].1;
    let mut r = s;
    if nb > 0 {
        let mut a = [R::ZERO; 8];
        a[0] = s;
        for i in 0..nb {
            a[i & 7] -= prod(i);
        }
        r = a[0] + a[1] + a[2] + a[3] + a[4] + a[5] + a[6] + a[7];
    }
    for i in (nb..n).rev() {
        r -= prod(i);
    }
    r
}

/// The extended-precision eta dot product (`0x25e0e0`): blocks of eight terms go to eight accumulators (position mod 8), which are added
/// in order `((((((a0 + a1) + a2) + a3) + a4) + a5) + a6) + a7` (the decompiled order is misleading; see the stack code at `0x265587`); the remaining terms are then added one by one from the last to the first.
fn dot_q8<R: Real>(x: &[R], c: &[R]) -> R {
    let n = c.len();
    let nb = n / 8 * 8;
    let mut tot = R::ZERO;
    if nb > 0 {
        let mut a = [R::ZERO; 8];
        for i in 0..nb {
            a[i & 7] += x[i] * c[i];
        }
        tot = a[0] + a[1] + a[2] + a[3] + a[4] + a[5] + a[6] + a[7];
    }
    for i in (nb..n).rev() {
        tot += x[i] * c[i];
    }
    tot
}

/// The extended-precision dot of the update's new diagonal: two accumulators by position parity over the even part, added, then the odd last term (seed 2077: 6 terms match only this
/// order; seed 2147: 5 terms match only with the last term added after the accumulators).
fn dot_diag<R: Real>(x: &[R], c: &[R], peel: usize) -> R {
    if R::LANES {
        return dot(x, c, peel);
    }
    let n = c.len();
    let even = n & !1;
    let mut l = [R::ZERO; 2];
    for i in 0..even {
        l[i & 1] += x[i] * c[i];
    }
    let mut tot = l[0] + l[1];
    if even < n {
        // The odd last term is added after the two accumulators (seed 2147: 5 terms).
        tot += x[even] * c[even];
    }
    tot
}

/// Eta dot product as the reference evaluates it: the vectorised loop peels the elements before the first 16-byte boundary of
/// the value array (`peel` is 0 or 1) and keeps eight lane-wise partial sums; below `peel + 8` terms it is sequential. The extended-precision
/// code (`R::LANES` false) uses `dot_q8`.
fn dot<R: Real>(x: &[R], c: &[R], peel: usize) -> R {
    let n = c.len();
    let t = |i: usize| x[i] * c[i];
    if !R::LANES {
        return dot_q8(x, c);
    }
    if n < peel + 8 {
        let mut acc = R::ZERO;
        for i in 0..n {
            acc += t(i);
        }
        return acc;
    }
    let n8 = n - ((n - peel) & 7);
    let mut lanes = [R::ZERO; 8];
    for i in 0..peel {
        lanes[0] += t(i);
    }
    for i in peel..n8 {
        lanes[(i - peel) & 7] += t(i);
    }
    // The compiled reduction adds the four packed accumulators pairwise: (l0+l2 | l1+l3) + (l4+l6 | l5+l7), then the two halves.
    let mut tot = ((lanes[0] + lanes[2]) + (lanes[4] + lanes[6]))
        + ((lanes[1] + lanes[3]) + (lanes[5] + lanes[7]));
    for i in n8..n {
        tot += t(i);
    }
    tot
}

impl<R: Real> RefStruct<R> {
    /// Install all fields of a cached factor into an exclusively owned slot.
    /// Nested Vec::clone_from dispatch retains the old column and row buffers.
    pub(crate) fn restore(&mut self, fresh: &Self) {
        macro_rules! restore_fields {
            ($($name:ident),* $(,)?) => {{
                // Exhaustive destructuring keeps new fields from escaping restoration.
                let Self { $($name),* } = self;
                $( $name.clone_from(&fresh.$name); )*
            }};
        }
        restore_fields!(
            m,
            row_cached,
            kind,
            nuc,
            i21,
            mark,
            a,
            b,
            lstart,
            lentries,
            perm,
            row_of,
            dens,
            dense_est,
            dthr,
            empty_tail,
            slot,
            at,
            slot_counts,
            etas,
            bt_occupied,
            bt_work,
            lt,
            levels,
            level_cache,
            lt_levels,
            lev_built,
            diag,
            uc,
            ut,
            spike,
            spike_used,
            update_values,
            update_seen,
            lsnap,
            upd,
            sum_alpha,
            before_last,
            sum_eta,
            last_alpha,
            znz,
            n1,
            loff,
            nnz_l,
            nnz_l_tail,
            first_eta,
            work,
            work_parts,
            keep_work,
            keep_dense,
            tau_l,
            tau_u,
            tau_fin,
            f2,
            tau_terms,
            tau_row,
            fill,
            row_cap,
            row_par,
            rowl,
            scratch_peel,
            upd_work,
            h33,
            h34,
            lev_dense,
            lev_dense_ft,
            last_dense,
            last_ordered,
            bt_dense,
            upd_limit,
        );
    }

    /// Fresh factor from the pivot sequence; `cols[pos]` is the basis column at incoming position `pos` as `(is_slack, rows, values)`.
    pub fn fresh(
        m: usize,
        ncols: usize,
        st: &Structure<R>,
        cols: &[(bool, Vec<usize>, Vec<f64>)],
    ) -> RefStruct<R> {
        let mut step = vec![usize::MAX; m];
        for (s, &(r, _)) in st.pivots.iter().enumerate() {
            step[r] = s;
        }
        let n1 = st.phase.iter().filter(|&&p| p == 1 || p == 2).count();
        let n3 = st.phase.iter().filter(|&&p| p == 3).count();
        let nseq = st.loff;
        let i21 = n1 + n3;
        let mut kind = vec![-1i64; m];
        let mut nuc = Vec::new();
        for s in n1..i21 {
            kind[s] = nuc.len() as i64;
            nuc.push(s);
        }
        for s in i21..m {
            kind[s] = -2;
        }
        let mut lstart = Vec::with_capacity(m + 1);
        let mut lentries = Vec::with_capacity(st.l_vals.iter().map(Vec::len).sum());
        lstart.push(0);
        for (rows, vals) in st.l_rows.iter().zip(&st.l_vals) {
            lentries.extend(rows.iter().zip(vals).map(|(&r, &v)| (step[r], v)));
            lstart.push(lentries.len());
        }
        let lcolumn = |v: usize| &lentries[lstart[v]..lstart[v + 1]];
        let mut lt = vec![Vec::new(); m];
        for v in 0..m {
            for &(y, a) in lcolumn(v) {
                lt[y].push((v, a));
            }
        }
        // pivot-row values of the Markowitz steps, by column
        let mut mk: Vec<Vec<(usize, R)>> = vec![Vec::new(); m];
        for s in nseq..m {
            for &(c, v) in &st.u_rows[s] {
                mk[c].push((s, v));
            }
        }
        let mut a = vec![Vec::new(); m];
        let mut b = vec![Vec::new(); m];
        let mut uc: Vec<Vec<(usize, R)>> = vec![Vec::new(); m];
        let (mut sa, mut sb) = (0usize, 0usize);
        let mut xv = vec![R::ZERO; m];
        let mut reached = vec![false; m];
        let mut epos = vec![usize::MAX; m];
        let mut changed = vec![false; m];
        for t in i21..m {
            let c = st.pivots[t].1;
            let (_, idx, val) = &cols[c];
            let mut start: Vec<usize> = Vec::new();
            for (&i, &v) in idx.iter().zip(val) {
                let s = step[i];
                if s <= t {
                    xv[s] += R::c(v);
                    start.push(s);
                }
            }
            // Only the reached set is needed: the numeric solve below visits it in ascending order.
            let mut lst = Vec::new();
            for &v in &start {
                if !reached[v] {
                    reached[v] = true;
                    lst.push(v);
                }
            }
            let mut k = 0;
            while k < lst.len() {
                let v = lst[k];
                k += 1;
                for &(y, _) in lcolumn(v) {
                    if y <= t && !reached[y] {
                        reached[y] = true;
                        lst.push(y);
                    }
                }
            }
            // The numeric solve visits the steps in ascending order, not graph-discovery order: sums that cancel to rounding
            // noise then fall on the reference's side of the tolerance (PILOT.JA factor 45, column 925).
            let mut q = lst.clone();
            q.sort_unstable();
            for &v in &q {
                if v >= t {
                    continue;
                }
                let x = xv[v];
                // Noise at or below the tolerance does not propagate (a cancellation residue of 1.2e-13 made PEROLD's tail
                // column 600 carry a spurious entry).
                if x == R::ZERO || x.abs() <= R::c(TOL) {
                    continue;
                }
                for (e, &(y, multiplier)) in lcolumn(v).iter().enumerate() {
                    if y > t {
                        continue;
                    }
                    // A singleton-phase column subtracts `(x * a) / d` with the raw entry `a`, not `x * (a / d)`; the
                    // multiplier is the only form known for pair overrides and nucleus columns.
                    let a = st.l_raw[v][e];
                    let upd = if a.is_nan() {
                        multiplier * x
                    } else {
                        (x * a) / st.diag[v]
                    };
                    // A product at or below the tolerance does not create an entry (as in the Markowitz elimination's fill): PILOT.JA
                    // factor 48, column 922: 9.06e-14 + 5.95e-14 would make a spurious entry of -1.5e-13 at 901.
                    if xv[y] == R::ZERO && upd.abs() <= R::c(TOL) {
                        continue;
                    }
                    xv[y] -= upd;
                    changed[y] = true;
                    // A sum that cancels to the tolerance leaves the pattern, as an eliminated entry does (PILOT.JA factor 52,
                    // column 931: -3.2e-14 left by a cancellation plus a later 8.4e-14 must not become an entry of -1.2e-13).
                    if xv[y].abs() <= R::c(TOL) {
                        xv[y] = R::ZERO;
                    }
                }
            }
            let mut ents: Vec<(usize, R)> = Vec::new();
            let touched = q;
            for &v in &touched {
                if v < t
                    && (xv[v].abs() > R::c(TOL)
                        || (!changed[v] && xv[v] != R::ZERO && (kind[v] >= 0 || st.phase[v] == 4)))
                {
                    ents.push((v, xv[v]));
                }
            }
            for (i, e) in ents.iter().enumerate() {
                epos[e.0] = i;
            }
            for &(sidx, vv) in &mk[c] {
                if sidx < t {
                    if epos[sidx] != usize::MAX {
                        ents[epos[sidx]].1 = vv;
                    } else {
                        epos[sidx] = ents.len();
                        ents.push((sidx, vv));
                    }
                }
            }
            for e in &ents {
                epos[e.0] = usize::MAX;
            }
            // A singleton-row pivot copies active entries without a magnitude test (`0xe0872`, `0xe2a14`, `0xdc546`), including tail positions.
            // Unchanged nucleus entries are also retained (fx_e2_r25).
            ents.retain(|e| {
                e.1.abs() > R::c(TOL) || (!changed[e.0] && (kind[e.0] >= 0 || st.phase[e.0] == 4))
            });
            for &v in &lst {
                xv[v] = R::ZERO;
                reached[v] = false;
            }
            for &v in &touched {
                changed[v] = false;
            }
            for &v in &touched {
                xv[v] = R::ZERO;
            }
            ents.sort_unstable_by_key(|e| e.0);
            let av: Vec<usize> = ents
                .iter()
                .map(|e| e.0)
                .filter(|&x| kind[x] == -2)
                .collect();
            let mut bv: Vec<usize> = ents.iter().map(|e| e.0).filter(|&x| kind[x] >= 0).collect();
            bv.sort_unstable_by_key(|&x| std::cmp::Reverse(kind[x]));
            sa += av.len();
            sb += bv.len();
            // The reference lists the fresh tail column's A entries in descending order.
            a[t] = av.into_iter().rev().collect();
            b[t] = bv;
            uc[t] = ents;
        }
        let nt = m - i21;
        let nn = nuc.len();
        let h33 = if nt > 0 {
            (sa + nt) as f64 / nt as f64
        } else {
            1.5
        };
        let h34 = if nn > 0 {
            (sb + nn) as f64 / nn as f64
        } else {
            1.5
        };
        let nnz: usize = lentries.len();
        let dthr = if nnz <= m {
            2
        } else {
            let t = {
                let q = (nnz as f64 / m as f64) + 1.0;
                5.0 * (q * q)
            };
            if (m as f64) > t {
                t as usize
            } else {
                m
            }
        };
        let mut lev_dense = false;
        let mut lev_built = false;
        let mut lev_dense_ft = false;
        let nnz_tail: usize = lentries.len() - lstart[nseq];
        let levels = {
            let d = m - nseq;
            let outer = d < 6 || 10 * d * d <= m || nnz_tail <= 10 * d;
            if m >= 250 {
                let mut lev = vec![0usize; m];
                for v in 0..m {
                    for &(y, _) in lcolumn(v) {
                        lev[y] = lev[y].max(lev[v] + 1);
                    }
                }
                let maxlev = lev.iter().copied().max().unwrap_or(0);
                lev_dense_ft = (maxlev as f64) >= 0.1 * m as f64;
                lev_dense = lev_dense_ft;
                lev_built = true;
                if !lev_dense_ft && outer {
                    Some((lev, maxlev))
                } else {
                    None
                }
            } else {
                None
            }
        };
        let nnz_l_tot: usize = lentries.len();
        let perm = step;
        let znz = (nnz + sa + sb + 2 * m) as f64;
        let mark0 = (i21..m)
            .rev()
            .find(|&t| kind[t] == -2 && a[t].is_empty() && b[t].is_empty())
            .unwrap_or(i21.min(m - 1));
        let mut ut = vec![Vec::new(); m];
        for (c, col) in uc.iter().enumerate() {
            for &(row, v) in col {
                ut[row].push((c, v));
            }
        }
        RefStruct {
            upd: 0,
            sum_alpha: 0.0,
            before_last: (0, 0.0),
            sum_eta: 0.0,
            last_alpha: 0,
            znz,
            n1: st.phase.iter().filter(|&&p| p == 1).count(),
            loff: nseq,
            nnz_l: nnz_l_tot,
            nnz_l_tail: nnz_tail,
            first_eta: vec![usize::MAX; m],
            work: 0.0,
            work_parts: (0.0, 0.0),
            keep_work: 0.0,
            keep_dense: false,
            tau_l: std::cell::Cell::new(0.0),
            tau_u: std::cell::Cell::new(0.0),
            tau_fin: std::cell::Cell::new(0),
            f2: std::cell::RefCell::new(F2::default()),
            tau_terms: Default::default(),
            tau_row: Default::default(),
            scratch_peel: ((m as f64).powf(0.7) as usize)
                .max(20)
                .wrapping_add(ncols + 1)
                % 2,
            fill: 2 * (uc.iter().map(Vec::len).sum::<usize>() + nnz_l_tot),
            // A row lists its tail columns in the order of the factorization's U row (Markowitz steps), or in the order of the
            // basis columns for a row of the singleton phases.
            rowl: {
                let mut member = vec![Vec::new(); m];
                for t in i21..m {
                    if kind[t] == -2 {
                        for &y in a[t].iter().chain(b[t].iter()) {
                            member[y].push(t);
                        }
                    }
                }
                let mut pos_of_col = vec![usize::MAX; m];
                for (t, &(_, c)) in st.pivots.iter().enumerate() {
                    pos_of_col[c] = t;
                }
                let mut rl = vec![Vec::new(); m];
                let mut stamp = vec![usize::MAX; m];
                for y in 0..m {
                    if member[y].is_empty() {
                        continue;
                    }
                    let order: Vec<usize> = if let Some(u) = st.pair_urow.get(&y) {
                        u.iter().map(|&c| pos_of_col[c]).collect()
                    } else if st.u_rows[y].is_empty() {
                        // Filtering all basis columns in index order is equivalent to sorting just this row's members by
                        // column. Keep only the positions represented by pos_of_col, and emit each position once.
                        let mut row: Vec<usize> = member[y]
                            .iter()
                            .copied()
                            .filter(|&t| pos_of_col[st.pivots[t].1] == t)
                            .collect();
                        row.sort_unstable_by_key(|&t| st.pivots[t].1);
                        row.dedup();
                        rl[y] = row;
                        continue;
                    } else {
                        st.u_rows[y].iter().map(|&(c, _)| pos_of_col[c]).collect()
                    };
                    for &t in &member[y] {
                        stamp[t] = y;
                    }
                    rl[y] = order.into_iter().filter(|&t| stamp[t] == y).collect();
                }
                rl
            },
            row_cap: {
                let mut c = vec![0usize; m];
                for col in &uc {
                    for &(row, _) in col {
                        c[row] += 1;
                    }
                }
                c
            },
            // A fresh factor lays the rows of the row-wise copy out contiguously: the nonempty rows of the nucleus in the order of their list index
            // (`kind`), then the others by position, each row as long as its U entries (`ucol_len`); the parity of the start decides the alignment of the row's
            // value array in the dot product of the second FTRAN vector.
            row_par: {
                let mut len = vec![0usize; m];
                for col in &uc {
                    for &(row, _) in col {
                        len[row] += 1;
                    }
                }
                let mut order: Vec<usize> =
                    (0..m).filter(|&p| kind[p] >= 0 && len[p] > 0).collect();
                order.sort_by_key(|&p| kind[p]);
                order.extend((0..m).filter(|&p| kind[p] < 0 && len[p] > 0));
                let mut par = vec![0u8; m];
                let mut acc = 0usize;
                for p in order {
                    par[p] = (acc % 2) as u8;
                    acc += len[p];
                }
                par
            },
            upd_work: Vec::new(),
            h33,
            h34,
            lev_dense,
            lev_dense_ft,
            last_dense: false,
            last_ordered: false,
            upd_limit: None,
            bt_dense: false,
            row_cached: None,
            m,
            empty_tail: (0..m)
                .filter(|&x| kind[x] == -2 && a[x].is_empty() && b[x].is_empty())
                .collect(),
            kind,
            nuc,
            i21,
            mark: mark0,
            a,
            b,
            lstart,
            lentries,
            perm,
            row_of: Default::default(),
            dens: (h33 + h34) * (h33 + h34),
            dense_est: std::cell::Cell::new(0.0),
            dthr,
            slot: (0..m).collect(),
            at: (0..m).collect(),
            slot_counts: Default::default(),
            etas: Vec::new(),
            lt,
            levels,
            level_cache: LevelCache::default(),
            lt_levels: Default::default(),
            // Four or more occupancy words, and skipped zero divisions must always be exact.
            bt_occupied: !R::QUAD && m >= 256 && st.diag.iter().all(|d| d.to_f64().is_normal()),
            bt_work: None,
            lev_built,
            diag: st.diag.clone(),
            uc,
            ut,
            spike: Vec::new(),
            update_values: Vec::new(),
            update_seen: Vec::new(),
            spike_used: false,
            lsnap: Default::default(),
        }
    }

    #[inline]
    fn lcolumn(&self, v: usize) -> &[(usize, R)] {
        &self.lentries[self.lstart[v]..self.lstart[v + 1]]
    }

    fn lcolumns(&self) -> impl Iterator<Item = &[(usize, R)]> {
        self.lstart.windows(2).map(|s| &self.lentries[s[0]..s[1]])
    }

    /// Compact dense output without a conditional push for each surviving value.
    #[inline(never)]
    fn dense_spike(x: &[R]) -> Vec<(usize, R)> {
        let mut out = Vec::with_capacity(x.len());
        let dst: *mut (usize, R) = out.as_mut_ptr();
        let mut count = 0;
        for (i, &v) in x.iter().enumerate() {
            // Before iteration i, count <= i < x.len() <= capacity. Rejected
            // entries overwrite the next unused slot; each accepted entry
            // advances count, retaining exactly the original order and bits.
            unsafe {
                dst.add(count).write((i, v));
            }
            count += usize::from(v != R::ZERO);
        }
        // Every slot below count holds an accepted entry. R is Copy, so
        // overwriting or leaving the rejected entries outside len needs no drop.
        unsafe {
            out.set_len(count);
        }
        out
    }

    /// Pack a dense BTRAN result in row order without conditional vector growth.
    #[inline(never)]
    fn dense_btran_output(x: &[R], perm: &[usize]) -> Vec<(usize, R)> {
        let mut out = Vec::with_capacity(perm.len());
        let dst: *mut (usize, R) = out.as_mut_ptr();
        let mut count = 0;
        for (r, &i) in perm.iter().enumerate() {
            let v = x[i];
            // count <= r < perm.len() <= capacity before this write. Rejected entries
            // occupy the next unused slot until a later accepted value overwrites it.
            unsafe {
                dst.add(count).write((r, v));
            }
            count += usize::from(v != R::ZERO);
        }
        // All slots below count are initialized accepted entries in row order. R is Copy,
        // so overwritten and rejected entries outside len require no destruction.
        unsafe {
            out.set_len(count);
        }
        out
    }

    fn kahn(&self, start: &[usize]) -> Vec<usize> {
        let mut cnt = vec![0i32; self.m];
        let mut lst: Vec<usize> = start.to_vec();
        let n0 = lst.len();
        for &v in &lst {
            cnt[v] += 1;
        }
        let mut i = 0;
        while i < lst.len() {
            let v = lst[i];
            i += 1;
            for &(x, _) in self.lcolumn(v) {
                if cnt[x] == 0 {
                    lst.push(x);
                }
                cnt[x] += 1;
            }
        }
        let mut q: Vec<usize> = Vec::new();
        for &v in &lst[..n0] {
            cnt[v] -= 1;
            if cnt[v] == 0 {
                q.push(v);
            }
        }
        let mut i = 0;
        while i < q.len() {
            let v = q[i];
            i += 1;
            for &(x, _) in self.lcolumn(v) {
                cnt[x] -= 1;
                if cnt[x] == 0 {
                    q.push(x);
                }
            }
        }
        q
    }

    /// L solve (Kahn order) and row etas on the entering column; returns the nonzero list the U solve receives, with values.
    /// Also returns the work estimate of the L stage and the row etas (`0x12bd30` / `0x12a990` and the eta loop of `0x120740`).
    fn lstage_num(
        &self,
        rows: &[usize],
        vals: &[R],
        force_dense: bool,
    ) -> (Vec<(usize, R)>, usize, f64) {
        let m = self.m;
        let mut e = 0.0f64;
        let dense_l = rows.len() * self.dthr > self.m || force_dense || self.lev_dense_ft;
        if R::LANES && !dense_l {
            if let Some((levels, maxlev)) = &self.levels {
                let mut cache = self.level_cache.0.borrow_mut();
                if let Some(result) = cache
                    .get_or_insert_with(Box::default)
                    .solve(self, levels, *maxlev, rows, vals)
                {
                    return result;
                }
            }
        }
        let mut x = vec![R::ZERO; self.m];
        let mut start: Vec<usize> = if dense_l {
            Vec::new()
        } else {
            Vec::with_capacity(rows.len())
        };
        for (&r, &v) in rows.iter().zip(vals) {
            let s = self.perm[r];
            x[s] += v;
            if !dense_l {
                start.push(s);
            }
        }
        let (mut order, count_l) = if dense_l {
            // Dense L solve (`0x12bd30`): every node in ascending order; the list handed on is the ascending scan.
            e = m as f64 + rows.len() as f64 * 2.0 + e;
            e = (m - self.loff) as f64 * 2.0 + self.n1 as f64 + self.n1 as f64 + e;
            for v in 0..self.m {
                let xv = x[v];
                if xv.abs() <= R::c(TOL) {
                    continue;
                }
                e = self.lcolumn(v).len() as f64 * 2.0 + e;
                for &(y, a) in self.lcolumn(v) {
                    x[y] -= a * xv;
                }
            }
            let count_l = (0..self.m).filter(|&v| x[v].abs() > R::c(TOL)).count();
            (Vec::new(), count_l)
        } else if let Some((lev, maxlev)) = &self.levels {
            // Level-bucket L solve: nodes are taken level by level, each level in insertion order (inputs first, then the
            // nodes reached from lower levels); values below the tolerance are dropped without propagating.
            let maxlev = *maxlev;
            let mut buckets: Vec<Vec<usize>> = vec![Vec::new(); maxlev + 1];
            let mut active = vec![false; self.m];
            let mut l = usize::MAX;
            for &v in &start {
                if !active[v] {
                    active[v] = true;
                    buckets[lev[v]].push(v);
                    l = l.min(lev[v]);
                }
            }
            let mut remaining: usize = buckets.iter().map(Vec::len).sum();
            let mut out: Vec<usize> = Vec::new();
            while remaining > 0 && l <= maxlev {
                let mut i = 0;
                while i < buckets[l].len() {
                    let v = buckets[l][i];
                    i += 1;
                    remaining -= 1;
                    let xv = x[v];
                    if xv.abs() >= R::c(TOL) {
                        out.push(v);
                        for &(y, a) in self.lcolumn(v) {
                            if !active[y] {
                                active[y] = true;
                                buckets[lev[y]].push(y);
                                remaining += 1;
                            }
                            x[y] -= a * xv;
                        }
                    } else {
                        x[v] = R::ZERO;
                    }
                }
                l += 1;
            }
            let n = out.len();
            e = ((self.nnz_l + m + 2 * m) as f64 * 5.0 * n as f64) / (m + 1) as f64 + e;
            (out, n)
        } else {
            let q = self.kahn(&start);
            let x0 = x.clone();
            for &v in &q {
                let xv = x[v];
                // The sparse L solve drops a node at or below the tolerance without keeping or propagating it.
                if xv.abs() <= R::c(TOL) {
                    x[v] = R::ZERO;
                    continue;
                }
                for &(y, a) in self.lcolumn(v) {
                    x[y] -= a * xv;
                }
            }
            let count_l = q.iter().filter(|&&v| x[v].abs() > R::c(TOL)).count();
            let d = m - self.loff;
            let nnz_tail = self.nnz_l_tail;
            if !(d < 6 || 10 * d * d <= m || nnz_tail <= 10 * d) {
                // Third strategy of `0x12a990`: the head columns (below `loff`) in Kahn order, then the tail by an ascending
                // scan from the lowest tail index reached by an input or by a kept head node.
                let loff = self.loff;
                let head_cnt = q
                    .iter()
                    .filter(|&&v| v < loff && x[v].abs() > R::c(TOL))
                    .count();
                let mut lo = m;
                for &v in &start {
                    if v >= loff {
                        lo = lo.min(v);
                    }
                }
                for &v in q.iter().filter(|&&v| v < loff && x[v].abs() > R::c(TOL)) {
                    for &(y, _) in self.lcolumn(v) {
                        if y >= loff {
                            lo = lo.min(y);
                        }
                    }
                }
                let nnz_head = self.nnz_l - nnz_tail;
                let e1 = (((nnz_head + 2 * d) as f64) * 7.0 * head_cnt as f64) / (loff + 1) as f64;
                let mut e2 = (m - lo) as f64 * 2.0 + e1;
                for v in lo..m {
                    if x[v].abs() > R::c(TOL) {
                        e2 += self.lcolumn(v).len() as f64 * 2.0;
                    }
                }
                e = e2 + e;
            } else {
                e = (((self.nnz_l + m) as f64 + m as f64 + m as f64) * 7.0 * count_l as f64)
                    / (m + 1) as f64
                    + e;
            }
            if m >= 250 && count_l * self.dthr > m {
                // The reference's result is then dense and equals an ascending scan (PILOT4 call 219); the switch condition is
                // inferred from the data, not read from `0x12a990`.
                x = x0;
                for v in 0..self.m {
                    let xv = x[v];
                    if xv.abs() <= R::c(TOL) {
                        continue;
                    }
                    for &(y, a) in self.lcolumn(v) {
                        x[y] -= a * xv;
                    }
                }
                // With etas the final result is dense; without them this fallback
                // still returns the filtered ascending list, including its tolerance rule.
                (
                    if self.etas.is_empty() {
                        (0..self.m).collect()
                    } else {
                        Vec::new()
                    },
                    count_l,
                )
            } else {
                (q, count_l)
            }
        };
        if crate::dump::on("lstage") {
            self.lsnap.borrow_mut().0 = x.clone();
        }
        // The sparse L stage drops a node whose value is at or below the tolerance from its list; an eta that makes it nonzero again appends it
        // at the end (SCAGR7 iteration 173: `0x129aa0` start `[82, 125, 83]` leaves `[82, 49, 125]`).
        let output_dense = dense_l || (!self.etas.is_empty() && count_l * self.dthr > m);
        let mut have = Vec::new();
        if !output_dense {
            order.retain(|&v| x[v].abs() > R::c(TOL));
            have.resize(m, false);
            for &v in &order {
                have[v] = true;
            }
        }
        // Gather each eta into one buffer; its logical peel is stored in eta.peel.
        let mut eta_x: Vec<R> = Vec::new();
        if !self.etas.is_empty() {
            let netas = self.etas.len();
            if dense_l || count_l * self.dthr > m {
                e = self.sum_eta * 2.0 + e;
            } else {
                // Sparse eta pass: from the first eta that lists an input position; the list grows by the eta positions that
                // become nonzero.
                let mut inlist = vec![false; m];
                let mut c0 = 0usize;
                let mut first = netas;
                for &v in &order {
                    if x[v].abs() > R::c(TOL) && !inlist[v] {
                        inlist[v] = true;
                        c0 += 1;
                        first = first.min(self.first_eta[v]);
                    }
                }
                let mut xe = x.clone();
                let mut ee = (netas - first) as f64 * 3.0;
                for eta in &self.etas[first..] {
                    ee = eta.ind.len() as f64 * 2.0 + ee;
                    eta_x.clear();
                    eta_x.extend(eta.ind.iter().map(|&j| xe[j]));
                    let mut sv = dot(&eta_x, &eta.val, eta.peel);
                    if sv.abs() <= R::c(TOL) {
                        sv = R::ZERO;
                    } else if !inlist[eta.idx] {
                        inlist[eta.idx] = true;
                        c0 += 1;
                    }
                    xe[eta.idx] = sv;
                }
                e = ee + 3.0 * c0 as f64 + e;
            }
        }
        for eta in &self.etas {
            eta_x.clear();
            eta_x.extend(eta.ind.iter().map(|&j| x[j]));
            let mut s = dot(&eta_x, &eta.val, eta.peel);
            if s.abs() <= R::c(TOL) {
                s = R::ZERO;
            }
            if !output_dense && !have[eta.idx] && s != R::ZERO {
                order.push(eta.idx);
                have[eta.idx] = true;
            }
            x[eta.idx] = s;
        }
        if crate::dump::on("lstage") {
            self.lsnap.borrow_mut().1 = x.clone();
        }
        // A dense vector goes on as the raw array: the reference keeps the sub-tolerance residues of the L and eta stages.
        if output_dense {
            return (
                if R::LANES {
                    Self::dense_spike(&x)
                } else {
                    (0..m)
                        .filter(|&v| x[v] != R::ZERO)
                        .map(|v| (v, x[v]))
                        .collect()
                },
                count_l,
                e,
            );
        }
        (
            order
                .into_iter()
                .filter(|&v| x[v].abs() > R::c(TOL))
                .map(|v| (v, x[v]))
                .collect(),
            count_l,
            e,
        )
    }

    #[cfg(test)]
    fn frozen_lstage_num(
        &self,
        rows: &[usize],
        vals: &[R],
        force_dense: bool,
    ) -> (Vec<(usize, R)>, usize, f64) {
        let lcols: Vec<Vec<usize>> = self
            .lcolumns()
            .map(|c| c.iter().map(|e| e.0).collect())
            .collect();
        let lvals: Vec<Vec<R>> = self
            .lcolumns()
            .map(|c| c.iter().map(|e| e.1).collect())
            .collect();
        let m = self.m;
        let mut e = 0.0f64;
        let dense_l = rows.len() * self.dthr > self.m || force_dense || self.lev_dense_ft;
        let mut x = vec![R::ZERO; self.m];
        let mut start: Vec<usize> = if dense_l {
            Vec::new()
        } else {
            Vec::with_capacity(rows.len())
        };
        for (&r, &v) in rows.iter().zip(vals) {
            let s = self.perm[r];
            x[s] += v;
            if !dense_l {
                start.push(s);
            }
        }
        let (mut order, count_l) = if dense_l {
            // Dense L solve (`0x12bd30`): every node in ascending order; the list handed on is the ascending scan.
            e = m as f64 + rows.len() as f64 * 2.0 + e;
            e = (m - self.loff) as f64 * 2.0 + self.n1 as f64 + self.n1 as f64 + e;
            for v in 0..self.m {
                let xv = x[v];
                if xv.abs() <= R::c(TOL) {
                    continue;
                }
                e = lcols[v].len() as f64 * 2.0 + e;
                for (e, &y) in lcols[v].iter().enumerate() {
                    x[y] -= lvals[v][e] * xv;
                }
            }
            let count_l = (0..self.m).filter(|&v| x[v].abs() > R::c(TOL)).count();
            (Vec::new(), count_l)
        } else if let Some((lev, _)) = &self.levels {
            // Level-bucket L solve: nodes are taken level by level, each level in insertion order (inputs first, then the
            // nodes reached from lower levels); values below the tolerance are dropped without propagating.
            let maxlev = lev.iter().copied().max().unwrap_or(0);
            let mut buckets: Vec<Vec<usize>> = vec![Vec::new(); maxlev + 1];
            let mut active = vec![false; self.m];
            let mut l = usize::MAX;
            for &v in &start {
                if !active[v] {
                    active[v] = true;
                    buckets[lev[v]].push(v);
                    l = l.min(lev[v]);
                }
            }
            let mut remaining: usize = buckets.iter().map(Vec::len).sum();
            let mut out: Vec<usize> = Vec::new();
            while remaining > 0 && l <= maxlev {
                let mut i = 0;
                while i < buckets[l].len() {
                    let v = buckets[l][i];
                    i += 1;
                    remaining -= 1;
                    let xv = x[v];
                    if xv.abs() >= R::c(TOL) {
                        out.push(v);
                        for (e, &y) in lcols[v].iter().enumerate() {
                            if !active[y] {
                                active[y] = true;
                                buckets[lev[y]].push(y);
                                remaining += 1;
                            }
                            x[y] -= lvals[v][e] * xv;
                        }
                    } else {
                        x[v] = R::ZERO;
                    }
                }
                l += 1;
            }
            let n = out.len();
            e = ((self.nnz_l + m + 2 * m) as f64 * 5.0 * n as f64) / (m + 1) as f64 + e;
            (out, n)
        } else {
            let q = self.kahn(&start);
            let x0 = x.clone();
            for &v in &q {
                let xv = x[v];
                // The sparse L solve drops a node at or below the tolerance without keeping or propagating it.
                if xv.abs() <= R::c(TOL) {
                    x[v] = R::ZERO;
                    continue;
                }
                for (e, &y) in lcols[v].iter().enumerate() {
                    x[y] -= lvals[v][e] * xv;
                }
            }
            let count_l = q.iter().filter(|&&v| x[v].abs() > R::c(TOL)).count();
            let d = m - self.loff;
            let nnz_tail = self.nnz_l_tail;
            if !(d < 6 || 10 * d * d <= m || nnz_tail <= 10 * d) {
                // Third strategy of `0x12a990`: the head columns (below `loff`) in Kahn order, then the tail by an ascending
                // scan from the lowest tail index reached by an input or by a kept head node.
                let loff = self.loff;
                let head_cnt = q
                    .iter()
                    .filter(|&&v| v < loff && x[v].abs() > R::c(TOL))
                    .count();
                let mut lo = m;
                for &v in &start {
                    if v >= loff {
                        lo = lo.min(v);
                    }
                }
                for &v in q.iter().filter(|&&v| v < loff && x[v].abs() > R::c(TOL)) {
                    for &y in &lcols[v] {
                        if y >= loff {
                            lo = lo.min(y);
                        }
                    }
                }
                let nnz_head = self.nnz_l - nnz_tail;
                let e1 = (((nnz_head + 2 * d) as f64) * 7.0 * head_cnt as f64) / (loff + 1) as f64;
                let mut e2 = (m - lo) as f64 * 2.0 + e1;
                for v in lo..m {
                    if x[v].abs() > R::c(TOL) {
                        e2 += lcols[v].len() as f64 * 2.0;
                    }
                }
                e = e2 + e;
            } else {
                e = (((self.nnz_l + m) as f64 + m as f64 + m as f64) * 7.0 * count_l as f64)
                    / (m + 1) as f64
                    + e;
            }
            if m >= 250 && count_l * self.dthr > m {
                // The reference's result is then dense and equals an ascending scan (PILOT4 call 219); the switch condition is
                // inferred from the data, not read from `0x12a990`.
                x = x0;
                for v in 0..self.m {
                    let xv = x[v];
                    if xv.abs() <= R::c(TOL) {
                        continue;
                    }
                    for (e, &y) in lcols[v].iter().enumerate() {
                        x[y] -= lvals[v][e] * xv;
                    }
                }
                // With etas the final result is dense; without them this fallback
                // still returns the filtered ascending list, including its tolerance rule.
                (
                    if self.etas.is_empty() {
                        (0..self.m).collect()
                    } else {
                        Vec::new()
                    },
                    count_l,
                )
            } else {
                (q, count_l)
            }
        };
        if crate::dump::on("lstage") {
            self.lsnap.borrow_mut().0 = x.clone();
        }
        // The sparse L stage drops a node whose value is at or below the tolerance from its list; an eta that makes it nonzero again appends it
        // at the end (SCAGR7 iteration 173: `0x129aa0` start `[82, 125, 83]` leaves `[82, 49, 125]`).
        let output_dense = dense_l || (!self.etas.is_empty() && count_l * self.dthr > m);
        let mut have = Vec::new();
        if !output_dense {
            order.retain(|&v| x[v].abs() > R::c(TOL));
            have.resize(m, false);
            for &v in &order {
                have[v] = true;
            }
        }
        // Gather each eta into one buffer; its logical peel is stored in eta.peel.
        let mut eta_x: Vec<R> = Vec::new();
        if !self.etas.is_empty() {
            let netas = self.etas.len();
            if dense_l || count_l * self.dthr > m {
                e = self.sum_eta * 2.0 + e;
            } else {
                // Sparse eta pass: from the first eta that lists an input position; the list grows by the eta positions that
                // become nonzero.
                let mut inlist = vec![false; m];
                let mut c0 = 0usize;
                let mut first = netas;
                for &v in &order {
                    if x[v].abs() > R::c(TOL) && !inlist[v] {
                        inlist[v] = true;
                        c0 += 1;
                        first = first.min(self.first_eta[v]);
                    }
                }
                let mut xe = x.clone();
                let mut ee = (netas - first) as f64 * 3.0;
                for eta in &self.etas[first..] {
                    ee = eta.ind.len() as f64 * 2.0 + ee;
                    eta_x.clear();
                    eta_x.extend(eta.ind.iter().map(|&j| xe[j]));
                    let mut sv = dot(&eta_x, &eta.val, eta.peel);
                    if sv.abs() <= R::c(TOL) {
                        sv = R::ZERO;
                    } else if !inlist[eta.idx] {
                        inlist[eta.idx] = true;
                        c0 += 1;
                    }
                    xe[eta.idx] = sv;
                }
                e = ee + 3.0 * c0 as f64 + e;
            }
        }
        for eta in &self.etas {
            eta_x.clear();
            eta_x.extend(eta.ind.iter().map(|&j| x[j]));
            let mut s = dot(&eta_x, &eta.val, eta.peel);
            if s.abs() <= R::c(TOL) {
                s = R::ZERO;
            }
            if !output_dense && !have[eta.idx] && s != R::ZERO {
                order.push(eta.idx);
                have[eta.idx] = true;
            }
            x[eta.idx] = s;
        }
        if crate::dump::on("lstage") {
            self.lsnap.borrow_mut().1 = x.clone();
        }
        // A dense vector goes on as the raw array: the reference keeps the sub-tolerance residues of the L and eta stages.
        if output_dense {
            return (
                (0..m)
                    .filter(|&v| x[v] != R::ZERO)
                    .map(|v| (v, x[v]))
                    .collect(),
                count_l,
                e,
            );
        }
        (
            order
                .into_iter()
                .filter(|&v| x[v].abs() > R::c(TOL))
                .map(|v| (v, x[v]))
                .collect(),
            count_l,
            e,
        )
    }

    /// FTRAN of a sparse column: the nonzeros `(position, value)` of `B^-1 a` in the reference's output order. The spike is kept
    /// for the following `update`.
    pub fn ftran(&mut self, rows: &[usize], vals: &[R]) -> Vec<(usize, R)> {
        self.ftran_with(rows, vals, false)
    }

    /// `dense_in`: the second vector is handed over as a dense vector (count -1), so its L stage runs dense (`0x12bd30`) whatever its length.
    /// The second vector of the reference's two-vector FTRAN (`0x125f50`, its U stage `0x123690`): `B^-1 a` at the positions of
    /// the first vector's result only (`needed`), in that order. The L stage and the tail columns are solved as in `ftran`, the
    /// tail columns scattering into tail rows only. A needed singleton position is `x / diag`, a needed tail position its solved
    /// value, a needed nucleus position the row-wise dot product over the row copy (`rowl`), subtracted term by term from the L
    /// stage value and divided by the diagonal, without a tolerance test. Does not touch the spike of the last `ftran`.
    pub fn ftran_needed(
        &self,
        rows: &[usize],
        vals: &[R],
        needed: &[usize],
        dense_in: bool,
    ) -> Vec<(usize, R)> {
        let m = self.m;
        let (spike, count_l, e_l) = self.lstage_num(rows, vals, dense_in);
        self.tau_l.set(e_l);
        let tau_dense = rows.len() * self.dthr > self.m
            || dense_in
            || self.lev_dense_ft
            || (!self.etas.is_empty() && count_l * self.dthr > self.m);
        let tau_n = if tau_dense { None } else { Some(spike.len()) };
        // The inline dense joint loop of two dense vectors (`0x128c56..0x129152`) tests the raw value before the quotient: a tail column whose
        // raw value is at most the tolerance is zeroed even when the quotient exceeds it.
        // The same holds for the joint routine `0x12c300` (column sparse, tau dense; PILOT.JA call 249): whenever tau does not take the single-vector U stage
        // (`ok(n2)` of `joint_u` fails), the tail columns test the raw value as well.
        let joint_dense = !tau_n.is_some_and(|n| {
            self.nuc.len() + n < self.i21 || self.dens * (n as f64) < 2.0 * m as f64
        });
        let mut x = vec![R::ZERO; m];
        let mut touched = vec![false; m];
        let mut tail: Vec<usize> = Vec::new();
        for &(p, v) in &spike {
            x[p] = v;
            touched[p] = true;
            if !joint_dense && self.kind[p] == -2 {
                tail.push(p);
            }
        }
        let keep_tiny = !tau_dense && !joint_dense && self.tau_u_dense(&spike);
        if joint_dense || (!R::QUAD && keep_tiny) {
            // The ordinary dense tail loop (12403b..12405f) divides and stores
            // every slot before testing the quotient, including unreachable zero
            // slots. Their signed zeros are consumed by requested rows and F2.
            tail.clear();
            // Extra zero columns are zeroed by the joint dense kernel. Reached
            // columns keep the same descending slot order as the symbolic walk.
            tail.extend((self.i21..m).rev().map(|s| self.at[s]));
        } else {
            let mut k = 0;
            while k < tail.len() {
                let c = tail[k];
                k += 1;
                for &j in &self.a[c] {
                    if !touched[j] {
                        touched[j] = true;
                        tail.push(j);
                    }
                }
            }
            self.sort_slot_desc(&mut tail);
        }
        let mut tau_raw = vec![false; m];
        // The U stage `0x123690` solves the tail in one of three loops: the heap walk and the sparse-row loop zero a quotient of at most the tolerance, the dense
        // slot loop (`0x12403b`: `m <= h33^2 c` or `5 h33^2 t > m - i21`, see `tau_u_work`) stores the quotient and only skips its scatter (FIT2P call 14
        // keeps a tail value of 4e-16 that enters the row dots of the nucleus positions). A dense input (`tau_dense`) zeroes it (FIT2P call 49), and so does the
        // joint dense loop (`joint_dense`, fx_e2_r25 call 6: a tail quotient of -1.8e-14 is 0 in the reference).
        for &c in &tail {
            tau_raw[c] = x[c].abs() > R::c(TOL);
            let xc = x[c] / self.diag[c];
            if joint_dense && !tau_raw[c] {
                x[c] = R::ZERO;
                continue;
            }
            if xc.abs() <= R::c(TOL) {
                x[c] = if keep_tiny { xc } else { R::ZERO };
                continue;
            }
            x[c] = xc;
            for &(j, v) in &self.uc[c] {
                if self.kind[j] == -2 {
                    x[j] -= v * xc;
                }
            }
        }
        self.tau_u.set(self.tau_u_work(&spike, &x, needed.len()));
        {
            let mut f = self.f2.borrow_mut();
            f.tau_n = tau_n;
            f.x_tau = x.clone();
            f.tau_raw = tau_raw;
        }
        let mut terms = self.tau_terms.borrow_mut();
        let mut row = self.tau_row.borrow_mut();
        needed
            .iter()
            .map(|&p| {
                let v = match self.kind[p] {
                    -2 => x[p],
                    -1 => x[p] / self.diag[p],
                    _ => {
                        terms.clear();
                        row.prepare(m, &self.ut[p]);
                        row.append_terms(&self.rowl[p], &x, &mut terms);
                        let dot = if R::QUAD {
                            row_dot_q(x[p], &terms)
                        } else {
                            row_dot(x[p], &terms, self.row_par[p] as usize)
                        };
                        dot / self.diag[p]
                    }
                };
                (p, v)
            })
            .collect()
    }

    #[cfg(test)]
    fn frozen_ftran_needed(
        &self,
        rows: &[usize],
        vals: &[R],
        needed: &[usize],
        dense_in: bool,
    ) -> Vec<(usize, R)> {
        let m = self.m;
        let (spike, count_l, e_l) = self.lstage_num(rows, vals, dense_in);
        self.tau_l.set(e_l);
        let tau_dense = rows.len() * self.dthr > self.m
            || dense_in
            || self.lev_dense_ft
            || (!self.etas.is_empty() && count_l * self.dthr > self.m);
        let tau_n = if tau_dense { None } else { Some(spike.len()) };
        let mut x = vec![R::ZERO; m];
        let mut touched = vec![false; m];
        let mut tail: Vec<usize> = Vec::new();
        for &(p, v) in &spike {
            x[p] = v;
            touched[p] = true;
            if self.kind[p] == -2 {
                tail.push(p);
            }
        }
        let mut k = 0;
        while k < tail.len() {
            let c = tail[k];
            k += 1;
            for &j in &self.a[c] {
                if !touched[j] {
                    touched[j] = true;
                    tail.push(j);
                }
            }
        }
        self.sort_slot_desc(&mut tail);
        let mut tau_raw = vec![false; m];
        // The inline dense joint loop of two dense vectors (`0x128c56..0x129152`) tests the raw value before the quotient: a tail column whose
        // raw value is at most the tolerance is zeroed even when the quotient exceeds it.
        // The same holds for the joint routine `0x12c300` (column sparse, tau dense; PILOT.JA call 249): whenever tau does not take the single-vector U stage
        // (`ok(n2)` of `joint_u` fails), the tail columns test the raw value as well.
        let joint_dense = !tau_n.is_some_and(|n| {
            self.nuc.len() + n < self.i21 || self.dens * (n as f64) < 2.0 * m as f64
        });
        // The U stage `0x123690` solves the tail in one of three loops: the heap walk and the sparse-row loop zero a quotient of at most the tolerance, the dense
        // slot loop (`0x12403b`: `m <= h33^2 c` or `5 h33^2 t > m - i21`, see `tau_u_work`) stores the quotient and only skips its scatter (FIT2P call 14
        // keeps a tail value of 4e-16 that enters the row dots of the nucleus positions). A dense input (`tau_dense`) zeroes it (FIT2P call 49), and so does the
        // joint dense loop (`joint_dense`, fx_e2_r25 call 6: a tail quotient of -1.8e-14 is 0 in the reference).
        let keep_tiny = !tau_dense && !joint_dense && self.tau_u_dense(&spike);
        for &c in &tail {
            tau_raw[c] = x[c].abs() > R::c(TOL);
            let xc = x[c] / self.diag[c];
            if joint_dense && !tau_raw[c] {
                x[c] = R::ZERO;
                continue;
            }
            if xc.abs() <= R::c(TOL) {
                x[c] = if keep_tiny { xc } else { R::ZERO };
                continue;
            }
            x[c] = xc;
            for &(j, v) in &self.uc[c] {
                if self.kind[j] == -2 {
                    x[j] -= v * xc;
                }
            }
        }
        self.tau_u.set(self.tau_u_work(&spike, &x, needed.len()));
        {
            let mut f = self.f2.borrow_mut();
            f.tau_n = tau_n;
            f.x_tau = x.clone();
            f.tau_raw = tau_raw;
        }
        let mut terms = self.tau_terms.borrow_mut();
        let mut row = self.tau_row.borrow_mut();
        needed
            .iter()
            .map(|&p| {
                let v = match self.kind[p] {
                    -2 => x[p],
                    -1 => x[p] / self.diag[p],
                    _ => {
                        terms.clear();
                        row.prepare(m, &self.ut[p]);
                        row.append_terms(&self.rowl[p], &x, &mut terms);
                        let dot = if R::QUAD {
                            row_dot_q(x[p], &terms)
                        } else {
                            row_dot(x[p], &terms, self.row_par[p] as usize)
                        };
                        dot / self.diag[p]
                    }
                };
                (p, v)
            })
            .collect()
    }

    /// The U-stage work of the two-vector FTRAN when the reference does not run the two single-vector U stages: the joint routine `0x12c300` (the
    /// column is a sparse list and tau is not) or the inline dense one (`0x128742`, neither is). `None` for the usual pair. `out` holds the positions of the
    /// column's result. The selection follows `0x1286c3..0x12873c`.
    pub fn joint_u(&self, out: &[usize]) -> Option<f64> {
        let f = self.f2.borrow();
        let (m, i21) = (self.m, self.i21);
        let ok = |n: Option<usize>| {
            n.is_some_and(|n| self.nuc.len() + n < i21 || self.dens * (n as f64) < 2.0 * m as f64)
        };
        if ok(f.tau_n) {
            return None;
        }
        Some(if ok(f.col_n) {
            self.third_work(&f, out)
        } else {
            self.dense_pair_work(&f, out)
        })
    }

    /// Work of the output stage of both joint routines: `(row length + 2) 2` for each result position of the nucleus list.
    fn row_out_work(&self, out: &[usize]) -> f64 {
        let mut e = 0.0f64;
        for &p in out {
            if self.kind[p] >= 0 {
                e += (self.rowl[p].len() as f64 + 2.0) * 2.0;
            }
        }
        e
    }

    /// The joint routine `0x12c300` (column: sparse list, tau: dense).
    fn third_work(&self, f: &F2<R>, out: &[usize]) -> f64 {
        let (m, i21) = (self.m, self.i21);
        let mf = m as f64;
        let tail_rows = (m - i21) as f64;
        let cin = f.col_in.len();
        if cin == 0 {
            return 0.0;
        }
        let tau_pass = |c: usize| f.x_tau[c] != R::ZERO && f.tau_raw[c];
        let col_pass = |c: usize| f.x_col[c] != R::ZERO;
        let mut x2 = 0.0f64;
        if let Some(c1) = f.tau_n.filter(|&c| c > 0) {
            x2 = 2.0 * c1 as f64 + mf;
        }
        let (mut nuc_cnt, mut tail_cnt) = (0usize, 0usize);
        let sparse_in = mf > cin as f64 * self.dens;
        let mut nuc_in: Vec<usize> = Vec::new();
        if sparse_in {
            for &p in &f.col_in {
                match self.kind[p] {
                    -1 => {}
                    -2 => tail_cnt += 1,
                    _ => {
                        nuc_cnt += 1;
                        nuc_in.push(p);
                    }
                }
            }
        }
        let nuc_len = self.nuc.len();
        let use_dense_scan =
            !sparse_in || self.dens * (nuc_cnt + tail_cnt) as f64 > tail_rows + nuc_len as f64;
        let mut x5 = 0.0f64;
        let x1;
        if use_dense_scan {
            x2 += 5.0 * tail_rows;
            if sparse_in {
                x2 += 2.0 * (tail_cnt + nuc_cnt) as f64;
            }
            let mut x3 = tail_rows * 3.0 + x2;
            for t in (i21..m).rev() {
                let c = self.at[t];
                if tau_pass(c) {
                    x5 += (self.a[c].len() as f64 + 2.0) * 2.0;
                }
                if col_pass(c) {
                    x3 += 2.0 * (self.a[c].len() + self.b[c].len()) as f64;
                }
            }
            let x0 = 2.0 * nuc_len as f64 + x3;
            x1 = 2.0 * out.len() as f64 + x0;
        } else {
            let (mut r9, mut r15) = (nuc_cnt, nuc_cnt);
            let mut surv: Vec<usize> = Vec::new();
            for t in (i21..m).rev() {
                let c = self.at[t];
                let (tp, cp) = (tau_pass(c), col_pass(c));
                if tp || cp {
                    r9 += self.a[c].len();
                }
                if tp {
                    x5 += (self.a[c].len() as f64 + 2.0) * 2.0;
                }
                if cp {
                    r15 += self.b[c].len();
                    surv.push(c);
                }
            }
            let mut x0 = tail_rows * 5.0;
            x0 += x2;
            x0 += r9 as f64 * 3.0;
            if r15 as f64 * 2.5 > nuc_len as f64 {
                x0 += nuc_cnt as f64;
                for &c in &surv {
                    x0 += self.b[c].len() as f64 * 2.0 + 4.0;
                }
                x0 += nuc_len as f64 * 2.4;
            } else {
                let mut seen = vec![false; m];
                for &p in &nuc_in {
                    seen[p] = true;
                }
                let mut touched = nuc_cnt;
                for &c in &surv {
                    x0 += self.b[c].len() as f64 * 2.0 + 4.0;
                    for &j in &self.b[c] {
                        if !seen[j] {
                            seen[j] = true;
                            touched += 1;
                        }
                    }
                }
                x0 += 3.0 * touched as f64;
            }
            x1 = 2.0 * out.len() as f64 + x0;
        }
        x5 + (x1 + self.row_out_work(out))
    }

    /// The inline dense routine of `0x125f50` (`0x128742`, both vectors dense).
    fn dense_pair_work(&self, f: &F2<R>, out: &[usize]) -> f64 {
        let (m, i21) = (self.m, self.i21);
        let mf = m as f64;
        let mut x9 = 0.0f64;
        if let Some(n) = f.col_n.filter(|&n| n > 0) {
            x9 = mf + 2.0 * n as f64;
        }
        if let Some(n) = f.tau_n.filter(|&n| n > 0) {
            x9 += mf + 2.0 * n as f64;
        }
        let mut x5 = 5.0 * (m - i21) as f64;
        let mut x4 = 0.0f64;
        for t in (i21..m).rev() {
            let c = self.at[t];
            if f.x_col[c] != R::ZERO {
                x5 += ((self.a[c].len() + self.b[c].len()) as f64 + 2.0) * 2.0;
            }
            if f.x_tau[c] != R::ZERO && f.tau_raw[c] {
                x4 += (self.a[c].len() as f64 + 2.0) * 2.0;
            }
        }
        let last = if i21 == m {
            m as i64 - 1
        } else {
            i21 as i64 - 1
        };
        let mut x1 = (m - i21) as f64;
        x1 += last.max(0) as f64;
        x1 *= 5.0;
        x1 += x9;
        x4 += self.row_out_work(out);
        x1 += x5;
        x1 += x4;
        x1 += out.len() as f64;
        x1 += out.len() as f64;
        x1
    }

    /// Whether the U stage of the second vector runs its dense slot loop (`0x123690`: `m <= h33^2 c` for a spike of `c` entries, or `5 h33^2 t > m - i21` for `t` tail entries).
    fn tau_u_dense(&self, spike: &[(usize, R)]) -> bool {
        let (m, d) = (self.m as f64, self.h33);
        if m <= (d * d) * spike.len() as f64 {
            return true;
        }
        let t = spike.iter().filter(|e| self.kind[e.0] == -2).count() as f64;
        ((5.0 * d) * d) * t > (self.m - self.i21) as f64
    }

    /// Work of the tau solve's U stage (`0x123690`, which adds to `F+0x200`): `10 + 2 n + 6 c` for `n` needed positions and a spike of
    /// `c` entries; a dense spike (`m <= h33^2 c`) or one with many tail entries (`5 h33^2 t > m - i21`) books `3 (m - i21)` and `2 len` (its tail-row entries) per
    /// tail column with quotient above the tolerance (plus `2 t` in the second case) and `4 (m - i21)` at the end; otherwise the heap walk books `4 h33` per nonzero
    /// tail column and `4` per such column at the end. `(2 n) h34` is added last. `x` holds the solved tail values.
    fn tau_u_work(&self, spike: &[(usize, R)], x: &[R], n: usize) -> f64 {
        let m = self.m as f64;
        let (cnt, nf) = (spike.len() as f64, n as f64);
        let d = self.h33;
        let tail_rows = (self.m - self.i21) as f64;
        let mut e = ((10.0 + nf) + nf) + 6.0 * cnt;
        let mut dense = m <= (d * d) * cnt;
        if !dense {
            let t = spike.iter().filter(|e| self.kind[e.0] == -2).count() as f64;
            if ((5.0 * d) * d) * t > tail_rows {
                dense = true;
                e = 2.0 * t + e;
            }
        }
        debug_assert_eq!(dense, self.tau_u_dense(spike));
        let fin;
        if dense {
            e = 3.0 * tail_rows + e;
            for t in (self.i21..self.m).rev() {
                let c = self.at[t];
                if x[c].abs() > R::c(TOL) {
                    e += self.a[c].len() as f64 * 2.0;
                }
            }
            fin = tail_rows;
        } else {
            let nz = (self.i21..self.m)
                .filter(|&t| x[self.at[t]] != R::ZERO)
                .count() as f64;
            e = (d * 4.0) * nz + e;
            fin = nz;
        }
        self.tau_fin.set(fin as usize);
        ((nf + nf) * self.h34 + e) + 4.0 * fin
    }

    /// A plain FTRAN (`0x1205d0`) between the pivot column's solve and the column replacement, e.g. the bound-flip update: the
    /// result as `ftran` (`dense`: the input is a dense vector, as `ftran_basic`), but the spike of the pivot column and the
    /// per-solve statistics stay as they were.
    pub fn ftran_keep(&mut self, rows: &[usize], vals: &[R], dense: bool) -> Vec<(usize, R)> {
        let saved = (
            self.spike.clone(),
            self.spike_used,
            self.last_alpha,
            self.last_dense,
            self.last_ordered,
            self.work,
        );
        let out = self.ftran_with(rows, vals, dense);
        self.keep_work = self.work;
        self.keep_dense = self.last_dense;
        (
            self.spike,
            self.spike_used,
            self.last_alpha,
            self.last_dense,
            self.last_ordered,
            self.work,
        ) = saved;
        out
    }

    /// FTRAN of the right-hand side after a factorization: the reference always runs the dense L solve there, whatever the
    /// density, which fixes the order of the floating-point updates.
    pub fn ftran_basic(&mut self, rows: &[usize], vals: &[R]) -> Vec<(usize, R)> {
        // The solve does not leave a spike for the update routine: the discard of the next refactorization (`0xf461f`) sees the last pivot column's
        // spike, none on a fresh factor (`fx_cleanup_two`, two refactorizations without a pivot in between).
        let saved = (self.spike.clone(), self.spike_used);
        let out = self.ftran_with(rows, vals, true);
        (self.spike, self.spike_used) = saved;
        out
    }

    fn ftran_with(&mut self, rows: &[usize], vals: &[R], force_dense: bool) -> Vec<(usize, R)> {
        let m = self.m;
        let (spike, count_l, e_l) = self.lstage_num(rows, vals, force_dense);
        // The previous column snapshot has served every consumer. Reuse it as
        // numeric work and reset only blocks that could contain a prior write.
        let old = self.f2.get_mut();
        let mut x = std::mem::take(&mut old.x_col);
        let mut dirty = old.dirty_blocks.take();
        if x.len() != m {
            x = vec![R::ZERO; m];
            dirty = None;
        } else if let Some(blocks) = &mut dirty {
            for block in blocks.drain(..) {
                let start = block * 64;
                x[start..m.min(start + 64)].fill(R::ZERO);
            }
        } else {
            x.fill(R::ZERO);
        }
        for &(p, v) in &spike {
            x[p] = v;
        }
        // The dense U routine (`0x121ef0`) tests the raw value first; the sparse one (`0x1223e0`, at `0x12297b`) only the quotient.
        let u_dense = !R::LANES
            || rows.len() * self.dthr > m
            || (!self.etas.is_empty() && count_l * self.dthr > m)
            || !(self.nuc.len() + spike.len() < self.i21
                || self.dens * (spike.len() as f64) < 2.0 * m as f64)
            || self.lev_dense_ft;
        let bs = force_dense && crate::dump::on("bsolve");
        let mut touched = FtranTouched::new(m, u_dense, dirty.unwrap_or_default());
        let mut tail: Vec<usize> = Vec::new();
        for &(p, _) in &spike {
            touched.insert(p);
        }
        if u_dense && !bs {
            // Unreached columns contain zero and the dense kernel skips them before division.
            // Slot order preserves the arithmetic order of every reached column.
            tail.extend((self.i21..m).rev().map(|s| self.at[s]));
        } else {
            for &(p, _) in &spike {
                if self.kind[p] == -2 {
                    tail.push(p);
                }
            }
            let mut k = 0;
            while k < tail.len() {
                let c = tail[k];
                k += 1;
                for &j in &self.a[c] {
                    if touched.insert(j) {
                        tail.push(j);
                    }
                }
            }
            self.sort_slot_desc(&mut tail);
        }
        let x_l = if bs { x.clone() } else { Vec::new() };
        for &c in &tail {
            let mut xc = x[c];
            if u_dense && xc.abs() <= R::c(TOL) {
                x[c] = R::ZERO;
                continue;
            }
            xc /= self.diag[c];
            if xc.abs() <= R::c(TOL) {
                x[c] = R::ZERO;
                continue;
            }
            x[c] = xc;
            for &(j, v) in &self.uc[c] {
                x[j] -= v * xc;
                touched.insert(j);
            }
        }
        let x_t = if bs { x.clone() } else { Vec::new() };
        touched.finish();
        touched.for_each(|p| {
            if self.kind[p] != -2 {
                let v = x[p] / self.diag[p];
                x[p] = if (u_dense && x[p].abs() <= R::c(TOL)) || v.abs() <= R::c(TOL) {
                    R::ZERO
                } else {
                    v
                };
            }
        });
        if bs {
            crate::dump::emit(
                "bsolve",
                serde_json::json!({"l": jv(&x_l), "tail_done": jv(&x_t), "x": jv(&x), "tail": tail, "kind": self.kind, "slot": self.slot, "diag": jv(&self.diag), "u_dense": u_dense, "uc": ju(&self.uc), "lcols": self.lcolumns().map(|c| c.iter().map(|e| e.0).collect::<Vec<_>>()).collect::<Vec<_>>(), "lvals": self.lcolumns().map(|c| c.iter().map(|e| e.1.to_f64()).collect::<Vec<_>>()).collect::<Vec<_>>(), "perm": self.perm, "rows": rows, "vals": jv(&vals)}),
            );
        }
        // A dense result (the reference's `-1` count) is scanned by position index; a sparse one comes in the U solve's order.
        let dense = rows.len() * self.dthr > m
            || self.lev_dense_ft
            || (!self.etas.is_empty() && count_l * self.dthr > m)
            || !(self.nuc.len() + spike.len() < self.i21
                || self.dens * (spike.len() as f64) < 2.0 * m as f64);
        let inp: Vec<usize> = spike.iter().map(|e| e.0).collect();
        let vec_dense = rows.len() * self.dthr > m
            || self.lev_dense_ft
            || (!self.etas.is_empty() && count_l * self.dthr > m);
        let u_est = self.u_work(
            &inp,
            &x,
            !dense,
            if vec_dense { None } else { Some(spike.len()) },
        );
        self.work = u_est + e_l;
        self.work_parts = (e_l, u_est);
        self.last_dense = dense;
        self.last_ordered = dense || u_dense;
        let mut out = Vec::new();
        if dense || u_dense {
            out.extend(x.iter().copied().enumerate().filter(|&(_, v)| v != R::ZERO));
        } else {
            for p in self.uorder(&inp, &x) {
                if touched.take_sparse(p) && x[p] != R::ZERO {
                    out.push((p, x[p]));
                }
            }
            touched.for_each(|p| {
                if x[p] != R::ZERO {
                    out.push((p, x[p]));
                }
            });
        }
        let dirty_blocks = match touched {
            FtranTouched::Sparse { occupied, .. } => Some(occupied),
            FtranTouched::Dense(_) => None,
        };
        *self.f2.get_mut() = F2 {
            col_n: if vec_dense { None } else { Some(spike.len()) },
            col_in: inp,
            x_col: x,
            dirty_blocks,
            ..F2::<R>::default()
        };
        // The column replacement receives the spike in ascending order after a dense solve (the reference scans the dense vector),
        // and in L-stage order after a sparse one.
        let mut spike = spike;
        // Only a dense L-stage result is scanned by position; a sparse one keeps the L-stage order.
        if vec_dense {
            spike.sort_unstable_by_key(|e| e.0);
        }
        // The column replacement does not receive the residues.
        spike.retain(|e| e.1.abs() > R::c(TOL));
        self.spike = spike;
        self.spike_used = false;
        self.last_alpha = out.len();
        out
    }

    #[cfg(test)]
    fn frozen_ftran_with(
        &mut self,
        rows: &[usize],
        vals: &[R],
        force_dense: bool,
    ) -> Vec<(usize, R)> {
        let m = self.m;
        let (spike, count_l, e_l) = self.lstage_num(rows, vals, force_dense);
        let mut x = vec![R::ZERO; m];
        for &(p, v) in &spike {
            x[p] = v;
        }
        let mut touched = vec![false; m];
        let mut tail: Vec<usize> = Vec::new();
        for &(p, _) in &spike {
            touched[p] = true;
            if self.kind[p] == -2 {
                tail.push(p);
            }
        }
        let mut k = 0;
        while k < tail.len() {
            let c = tail[k];
            k += 1;
            for &j in &self.a[c] {
                if !touched[j] {
                    touched[j] = true;
                    tail.push(j);
                }
            }
        }
        self.sort_slot_desc(&mut tail);
        // The dense U routine (`0x121ef0`) tests the raw value first; the sparse one (`0x1223e0`, at `0x12297b`) only the quotient.
        let u_dense = !R::LANES
            || rows.len() * self.dthr > m
            || (!self.etas.is_empty() && count_l * self.dthr > m)
            || !(self.nuc.len() + spike.len() < self.i21
                || self.dens * (spike.len() as f64) < 2.0 * m as f64)
            || self.lev_dense_ft;
        let bs = force_dense && crate::dump::on("bsolve");
        let x_l = if bs { x.clone() } else { Vec::new() };
        for &c in &tail {
            let mut xc = x[c];
            if u_dense && xc.abs() <= R::c(TOL) {
                x[c] = R::ZERO;
                continue;
            }
            xc /= self.diag[c];
            if xc.abs() <= R::c(TOL) {
                x[c] = R::ZERO;
                continue;
            }
            x[c] = xc;
            for &(j, v) in &self.uc[c] {
                x[j] -= v * xc;
                touched[j] = true;
            }
        }
        let x_t = if bs { x.clone() } else { Vec::new() };
        for p in 0..m {
            if touched[p] && self.kind[p] != -2 {
                let v = x[p] / self.diag[p];
                x[p] = if (u_dense && x[p].abs() <= R::c(TOL)) || v.abs() <= R::c(TOL) {
                    R::ZERO
                } else {
                    v
                };
            }
        }
        if bs {
            crate::dump::emit(
                "bsolve",
                serde_json::json!({"l": jv(&x_l), "tail_done": jv(&x_t), "x": jv(&x), "tail": tail, "kind": self.kind, "slot": self.slot, "diag": jv(&self.diag), "u_dense": u_dense, "uc": ju(&self.uc), "lcols": self.lcolumns().map(|c| c.iter().map(|e| e.0).collect::<Vec<_>>()).collect::<Vec<_>>(), "lvals": self.lcolumns().map(|c| c.iter().map(|e| e.1.to_f64()).collect::<Vec<_>>()).collect::<Vec<_>>(), "perm": self.perm, "rows": rows, "vals": jv(&vals)}),
            );
        }
        // A dense result (the reference's `-1` count) is scanned by position index; a sparse one comes in the U solve's order.
        let dense = rows.len() * self.dthr > m
            || self.lev_dense_ft
            || (!self.etas.is_empty() && count_l * self.dthr > m)
            || !(self.nuc.len() + spike.len() < self.i21
                || self.dens * (spike.len() as f64) < 2.0 * m as f64);
        let inp: Vec<usize> = spike.iter().map(|e| e.0).collect();
        let vec_dense = rows.len() * self.dthr > m
            || self.lev_dense_ft
            || (!self.etas.is_empty() && count_l * self.dthr > m);
        let u_est = self.u_work(
            &inp,
            &x,
            !dense,
            if vec_dense { None } else { Some(spike.len()) },
        );
        *self.f2.borrow_mut() = F2 {
            col_n: if vec_dense { None } else { Some(spike.len()) },
            col_in: inp.clone(),
            x_col: x.clone(),
            ..F2::<R>::default()
        };
        self.work = u_est + e_l;
        self.work_parts = (e_l, u_est);
        self.last_dense = dense;
        self.last_ordered = dense || u_dense;
        let ord = if dense || u_dense {
            (0..m).collect()
        } else {
            self.uorder(&inp, &x)
        };
        let mut out: Vec<(usize, R)> = Vec::new();
        let mut done = vec![false; m];
        for &p in &ord {
            if !done[p] {
                done[p] = true;
                if x[p] != R::ZERO {
                    out.push((p, x[p]));
                }
            }
        }
        for p in 0..m {
            if touched[p] && !done[p] && x[p] != R::ZERO {
                out.push((p, x[p]));
            }
        }
        // The column replacement receives the spike in ascending order after a dense solve (the reference scans the dense vector),
        // and in L-stage order after a sparse one.
        let mut spike = spike;
        // Only a dense L-stage result is scanned by position; a sparse one keeps the L-stage order.
        if vec_dense {
            spike.sort_unstable_by_key(|e| e.0);
        }
        // The column replacement does not receive the residues.
        spike.retain(|e| e.1.abs() > R::c(TOL));
        self.spike = spike;
        self.spike_used = false;
        self.last_alpha = out.len();
        out
    }

    /// The reference's refactorization decision, evaluated before an iteration's own update: at the update limit
    /// `max(20, m^0.7)`, or from `U >= m^0.45` when the accumulated work `(1.5 S + 2 c) m` exceeds `min((Z - 1.5 m)^2, 4 m^2 U0)`.
    /// `(nucleus list length, first tail slot, dens)`; for diagnostics against the reference's `i3`, `i21` and `(hdr33+hdr34)^2`.
    /// The fresh factor's numerics for the `luout` dump: the diagonal and each position's L column `(position, value)`.
    pub fn numerics(&self) -> (&[R], Vec<Vec<(usize, R)>>) {
        let l = self.lcolumns().map(|c| c.to_vec()).collect();
        (&self.diag, l)
    }

    /// Fresh U columns in the stored solve order, indexed by pivot position.
    pub fn u_numerics(&self) -> &[Vec<(usize, R)>] {
        &self.uc
    }

    pub fn blists(&self) -> &[Vec<usize>] {
        &self.b
    }

    pub fn kinds(&self) -> &[i64] {
        &self.kind
    }

    pub fn dims(&self) -> (usize, usize, f64) {
        (self.nuc.len(), self.i21, self.dens)
    }

    pub fn updates(&self) -> usize {
        self.upd
    }

    /// Sum of the FTRAN result counts of the updates since the factorization (`F+0x190`).
    pub fn window_stats(&self) -> (usize, f64) {
        self.before_last
    }

    /// A refactorization decided without an update of its own iteration (the empty scan): every update counts in the statistics.
    pub fn count_all_updates(&mut self) {
        self.before_last = (self.upd, self.sum_alpha);
    }

    pub fn sum_alpha(&self) -> f64 {
        self.sum_alpha
    }

    /// The update limit `max(20, floor(m^0.7))`.
    pub fn max_updates(&self) -> usize {
        ((self.m as f64).powf(0.7) as usize).max(20)
    }

    pub fn refactor_due(&self) -> bool {
        let m = self.m as f64;
        let maxupd = self.upd_limit.unwrap_or((m.powf(0.7) as usize).max(20));
        if self.upd >= maxupd {
            return true;
        }
        let c46c = m.powf(0.45) as usize;
        if self.upd < c46c {
            return false;
        }
        let d6 = {
            let q = self.znz - 1.5 * m;
            q * q
        };
        let d7 = (4.0 * m * m * c46c as f64).min(d6);
        (1.5 * self.sum_alpha + 2.0 * self.sum_eta) * m - d7 > 0.0
    }

    /// Both sparse U-transpose schedules visit reached slots in the original order.
    fn sparse_btran_u<W: BtranTail<R>>(
        &self,
        inp: &[(usize, R)],
        x: &mut [R],
        w: &mut W,
        mut out: Vec<usize>,
    ) -> (Vec<usize>, f64) {
        let m = self.m;
        let cnt = inp.len();
        let ut = &self.ut;
        out.clear();
        let mut eu = 0.0f64;
        let push = |w: &mut W, row: &[(usize, R)], d: R| {
            for &(c, v) in row {
                let value = w.entry(c, &self.slot);
                *value = -(v * d) + *value;
            }
        };
        let mut lowest = m;
        if self.h33 > 5.0 {
            eu = cnt as f64 * 3.0;
        }
        for &(p, v) in inp {
            let d = v / self.diag[p];
            let tail = if self.h33 > 5.0 {
                self.kind[p] == -2
            } else {
                self.slot[p] >= self.i21
            };
            if tail {
                if self.h33 <= 5.0 || d.abs() > R::c(TOL) {
                    lowest = lowest.min(self.slot[p]);
                    *w.entry(p, &self.slot) += v;
                }
            } else if d.abs() > R::c(TOL) {
                out.push(p);
                x[p] = d;
                if self.h33 > 5.0 && self.kind[p] >= 0 {
                    eu = ut[p].len() as f64 * 2.0 + eu;
                }
                push(w, &ut[p], d);
                for &(c, _) in &ut[p] {
                    lowest = lowest.min(self.slot[c]);
                }
            }
        }
        if self.h33 > 5.0 {
            eu = (m - lowest) as f64 * 2.0 + eu;
        }
        // An untouched +0 divided by a normal diagonal is exact signed zero. With the
        // occupancy schedule those slots cannot pass the threshold and need no division.
        let mut next = lowest;
        while let Some(t) = w.next_slot(next, m) {
            next = t + 1;
            let p = self.at[t];
            let d = w.take(p, t) / self.diag[p];
            if d.abs() > R::c(TOL) {
                eu = ut[p].len() as f64 * 2.0 + eu;
                out.push(p);
                x[p] = d;
                push(w, &ut[p], d);
            }
        }
        (out, eu)
    }

    /// Original sparse U-transpose tail sweep, also used as the differential oracle.
    #[cfg(test)]
    fn sparse_btran_u_scan(&self, inp: &[(usize, R)], x: &mut [R]) -> (Vec<usize>, f64) {
        let m = self.m;
        let cnt = inp.len();
        let ut = &self.ut;
        let mut out: Vec<usize> = Vec::new();
        let mut w = vec![R::ZERO; m];
        let mut eu = 0.0f64;
        let push = |w: &mut Vec<R>, row: &[(usize, R)], d: R| {
            for &(c, v) in row {
                w[c] = -(v * d) + w[c];
            }
        };
        let mut lowest = m;
        if self.h33 > 5.0 {
            eu = cnt as f64 * 3.0;
        }
        for &(p, v) in inp {
            let d = v / self.diag[p];
            let tail = if self.h33 > 5.0 {
                self.kind[p] == -2
            } else {
                self.slot[p] >= self.i21
            };
            if tail {
                if self.h33 <= 5.0 || d.abs() > R::c(TOL) {
                    lowest = lowest.min(self.slot[p]);
                    w[p] += v;
                }
            } else if d.abs() > R::c(TOL) {
                out.push(p);
                x[p] = d;
                if self.h33 > 5.0 && self.kind[p] >= 0 {
                    eu = ut[p].len() as f64 * 2.0 + eu;
                }
                push(&mut w, &ut[p], d);
                for &(c, _) in &ut[p] {
                    lowest = lowest.min(self.slot[c]);
                }
            }
        }
        if self.h33 > 5.0 {
            eu = (m - lowest) as f64 * 2.0 + eu;
        }
        for t in lowest..m {
            let p = self.at[t];
            let d = w[p] / self.diag[p];
            w[p] = R::ZERO;
            if d.abs() > R::c(TOL) {
                eu = ut[p].len() as f64 * 2.0 + eu;
                out.push(p);
                x[p] = d;
                push(&mut w, &ut[p], d);
            }
        }
        (out, eu)
    }

    /// The reference BTRAN (`0x12de60`) of a sparse input given by basis position: U^T solve (`0x130b70` or the dense loop), the row
    /// etas in reverse, the L^T solve, then the rows. Returns `(row, value)` with the reference's drops, and leaves its work
    /// estimate in `work`. `direct` is the update-saving entry (`0x131710`), which goes to the sparse U^T solve without the density
    /// test unless the first position's U^T row is empty. `dense_in` passes the input as a dense vector (count `-1`).
    pub fn btran(&mut self, inp: &[(usize, R)], direct: bool, dense_in: bool) -> Vec<(usize, R)> {
        let mut ws = self.bt_work.take().unwrap_or_default();
        let result = self.btran_with(inp, direct, dense_in, &mut ws);
        self.bt_work = Some(ws);
        result
    }

    fn btran_with(
        &mut self,
        inp: &[(usize, R)],
        direct: bool,
        dense_in: bool,
        ws: &mut BtranWork<R>,
    ) -> Vec<(usize, R)> {
        let m = self.m;
        let ut = &self.ut;
        let mut x = std::mem::take(&mut ws.x);
        x.resize(m, R::ZERO);
        let mut e = 0.0f64;
        let cnt = inp.len();
        self.bt_dense = false;
        if direct && inp.len() == 1 {
            self.row_cached = Some(inp[0].0);
        }
        let hsum = (self.h33 + self.h34) as usize;
        let direct = direct && !inp.is_empty() && !ut[inp[0].0].is_empty();
        // The current vector: a list (`Some`) or dense (`None`).
        let mut list: Option<Vec<usize>>;
        if !dense_in && (direct || hsum * cnt <= m) {
            let order = std::mem::take(&mut ws.input_order);
            let (out, eu) = if self.bt_occupied && cnt <= m / 64 {
                ws.tail.reset(m);
                self.sparse_btran_u(inp, &mut x, &mut ws.tail, order)
            } else {
                self.sparse_btran_u(inp, &mut x, &mut vec![R::ZERO; m], order)
            };
            e = eu + e;
            if out.is_empty() {
                ws.x = x;
                ws.input_order = out;
                self.work = e;
                return Vec::new();
            }
            list = Some(out);
        } else {
            for &(p, v) in inp {
                x[p] = v;
            }
            if !dense_in && cnt >= 1 {
                e = cnt as f64 * 2.0 + m as f64 + e;
            }
            e = m as f64 + m as f64 + e;
            for t in 0..m {
                let p = self.at[t];
                let xv = x[p];
                let d = xv / self.diag[p];
                if xv.abs() <= R::c(TOL) || d.abs() <= R::c(TOL) {
                    x[p] = R::ZERO;
                    continue;
                }
                x[p] = d;
                e = ut[p].len() as f64 * 2.0 + e;
                for &(c, v) in &ut[p] {
                    x[c] -= v * d;
                }
            }
            list = None;
        }
        // Row etas, last first.
        if !self.etas.is_empty() {
            let netas = self.etas.len();
            match &mut list {
                Some(l) if l.len() * self.dthr <= m => {
                    let c0 = l.len();
                    let stamp = ws.next_generation(m);
                    // A sparse active list can schedule reverse eta visits by pivot row.
                    // The heap is ordered by eta number, exactly as the original reverse scan.
                    let indexed = self.bt_occupied && netas / 4 > c0;
                    if indexed {
                        ws.index_etas(&self.etas, m);
                    }
                    for &p in l.iter() {
                        if indexed && ws.marks[p] != stamp {
                            ws.schedule_row(p, netas);
                        }
                        ws.marks[p] = stamp;
                    }
                    let mut es = netas as f64 + netas as f64 + c0 as f64 * 5.0;
                    let mut cursor = netas;
                    loop {
                        let i = if indexed {
                            let Some(e) = ws.pending.pop() else { break };
                            let prev = ws.eta_prev[e];
                            if prev != usize::MAX {
                                ws.pending.push(prev);
                            }
                            e
                        } else {
                            if cursor == 0 {
                                break;
                            }
                            cursor -= 1;
                            cursor
                        };
                        let eta = &self.etas[i];
                        let xr = x[eta.idx];
                        if xr == R::ZERO {
                            continue;
                        }
                        es = eta.ind.len() as f64 * 3.0 + es;
                        for (k, &j) in eta.ind.iter().enumerate() {
                            if j == eta.idx {
                                x[j] = eta.val[k] * xr;
                            } else {
                                x[j] = eta.val[k] * xr + x[j];
                                if ws.marks[j] != stamp {
                                    ws.marks[j] = stamp;
                                    l.push(j);
                                    if indexed {
                                        ws.schedule_row(j, i);
                                    }
                                }
                            }
                        }
                    }
                    let n = l.len();
                    e = es + (n - c0) as f64 * 2.0 + 3.0 * n as f64 + e;
                    l.retain(|&p| {
                        if x[p].abs() >= R::c(TOL) {
                            true
                        } else {
                            x[p] = R::ZERO;
                            false
                        }
                    });
                    if l.is_empty() {
                        ws.input_order = std::mem::take(l);
                        ws.x = x;
                        self.work = e;
                        return Vec::new();
                    }
                }
                _ => {
                    for eta in self.etas.iter().rev() {
                        let xr = x[eta.idx];
                        if xr.abs() <= R::c(TOL) {
                            x[eta.idx] = R::ZERO;
                            continue;
                        }
                        for (k, &j) in eta.ind.iter().enumerate() {
                            if j == eta.idx {
                                x[j] = eta.val[k] * xr;
                            } else {
                                x[j] = eta.val[k] * xr + x[j];
                            }
                        }
                    }
                    e = self.sum_eta * 2.0 + e;
                    if let Some(order) = list.take() {
                        ws.input_order = order;
                    }
                }
            }
        }
        // L^T: row i of L pushes into the earlier nodes.
        let lt = &self.lt;
        let mut res = std::mem::take(&mut ws.output_order);
        res.clear();
        let mut partial = false;
        let sparse_l = match &list {
            Some(l) => l.len() * self.dthr <= m && !self.lev_dense,
            None => false,
        };
        if sparse_l {
            let l = list.as_ref().unwrap();
            let cnt = l.len();
            let (n1, nseq) = (self.n1, self.loff);
            let d = m - nseq;
            let nnz_tail = self.nnz_l_tail;
            let mut hi = 0usize;
            if d >= 6 && m < 10 * d * d && n1 < 2 * d && 10 * d < nnz_tail && cnt > 0 {
                hi = l.iter().filter(|&&i| i >= nseq).count();
            }
            let dn = if d > 0 {
                nnz_tail as f64 * 0.5 / d as f64
            } else {
                0.0
            };
            if hi as f64 * dn * dn > d as f64 {
                // Partial: inputs in the middle range directly, then the top range and the first-phase range, both descending.
                let mut el = cnt as f64 * 3.0;
                let mut w = vec![R::ZERO; m];
                let mut maxi = 0usize;
                for &i in l {
                    let v = x[i];
                    if i < nseq {
                        if i < n1 {
                            w[i] += v;
                        } else {
                            res.push(i);
                            el = lt[i].len() as f64 * 2.0 + el;
                            for &(j, lv) in &lt[i] {
                                w[j] -= lv * v;
                            }
                        }
                    } else {
                        maxi = maxi.max(i);
                        w[i] = v;
                    }
                }
                for &i in &res {
                    w[i] = x[i];
                }
                el = (n1 as f64 - nseq as f64 + maxi as f64) * 2.0 + el;
                let ranges: Vec<usize> = if nseq <= maxi {
                    (nseq..=maxi).rev().collect()
                } else {
                    Vec::new()
                };
                for i in ranges.into_iter().chain((0..n1).rev()) {
                    let v = w[i];
                    if v.abs() > R::c(TOL) {
                        el = lt[i].len() as f64 * 2.0 + el;
                        res.push(i);
                        for &(j, lv) in &lt[i] {
                            w[j] -= lv * v;
                        }
                    } else {
                        w[i] = R::ZERO;
                    }
                }
                x = w;
                partial = true;
                e = el + e;
            } else {
                // Level buckets of the transposed graph when the factor has levels (`0x12bd30`), else a max-heap on the node index; every
                // reached node is kept when at least the tolerance.
                let stamp = ws.next_generation(m);
                if self.levels.is_some() {
                    // Level buckets of the transposed graph: a node comes after the nodes whose rows it appears in.
                    let (lev, maxlev) = self.lt_levels.get_or_init(|| {
                        let mut lev = vec![0usize; m];
                        for v in (0..m).rev() {
                            for &(y, _) in self.lcolumn(v) {
                                lev[v] = lev[v].max(lev[y] + 1);
                            }
                        }
                        let maxlev = lev.iter().copied().max().unwrap_or(0);
                        (lev, maxlev)
                    });
                    let maxlev = *maxlev;
                    let mut buckets: Vec<Vec<usize>> = vec![Vec::new(); maxlev + 1];
                    let mut lo = usize::MAX;
                    for &i in l {
                        ws.marks[i] = stamp;
                        buckets[lev[i]].push(i);
                        lo = lo.min(lev[i]);
                    }
                    let mut lv_i = lo;
                    while lv_i <= maxlev {
                        let mut k = 0;
                        while k < buckets[lv_i].len() {
                            let i = buckets[lv_i][k];
                            k += 1;
                            let v = x[i];
                            if v.abs() >= R::c(TOL) {
                                res.push(i);
                                for &(j, lv) in &lt[i] {
                                    x[j] -= lv * v;
                                }
                                // The reference walks the row from its last entry.
                                for &(j, _) in lt[i].iter().rev() {
                                    if ws.marks[j] != stamp {
                                        ws.marks[j] = stamp;
                                        buckets[lev[j]].push(j);
                                    }
                                }
                            } else {
                                x[i] = R::ZERO;
                            }
                        }
                        lv_i += 1;
                    }
                } else {
                    let mut heap: BinaryHeap<usize> = BinaryHeap::new();
                    for &i in l {
                        ws.marks[i] = stamp;
                        heap.push(i);
                    }
                    while let Some(i) = heap.pop() {
                        let v = x[i];
                        if v.abs() >= R::c(TOL) {
                            res.push(i);
                            for &(j, lv) in &lt[i] {
                                x[j] -= lv * v;
                                if ws.marks[j] != stamp {
                                    ws.marks[j] = stamp;
                                    heap.push(j);
                                }
                            }
                        } else {
                            x[i] = R::ZERO;
                        }
                    }
                }
                let f1e0 = (self.nnz_l + m) as f64;
                let k = if self.levels.is_some() { 5.0 } else { 7.0 };
                e = (res.len() as f64 * k * (m as f64 + m as f64 + f1e0)) / (m + 1) as f64 + e;
            }
        } else {
            self.bt_dense = true;
            if let Some(l) = &list {
                if !l.is_empty() {
                    e = l.len() as f64 * 2.0 + m as f64 + e;
                }
            }
            e = m as f64 + m as f64 + e;
            for i in (0..m).rev() {
                let v = x[i];
                if v.abs() > R::c(TOL) {
                    e = lt[i].len() as f64 * 2.0 + e;
                    for &(j, lv) in &lt[i] {
                        x[j] -= lv * v;
                    }
                }
            }
        }
        self.work = e;
        let row_of = self.row_of.get_or_init(|| {
            let mut row_of = vec![0usize; m];
            for (r, &v) in self.perm.iter().enumerate() {
                row_of[v] = r;
            }
            row_of
        });
        if let Some(order) = list {
            ws.input_order = order;
        }
        let result = if self.bt_dense {
            if R::LANES {
                Self::dense_btran_output(&x, &self.perm[..m])
            } else {
                (0..m)
                    .filter_map(|r| {
                        let v = x[self.perm[r]];
                        (v != R::ZERO).then_some((r, v))
                    })
                    .collect()
            }
        } else {
            res.iter()
                .copied()
                .filter(|&i| x[i] != R::ZERO)
                .map(|i| (row_of[i], x[i]))
                .collect()
        };
        if self.bt_dense || partial {
            x.fill(R::ZERO);
        } else {
            // Every reached L node was either filtered to +0 or emitted, including duplicate seeds.
            for &i in &res {
                x[i] = R::ZERO;
            }
        }
        ws.x = x;
        ws.output_order = res;
        result
    }

    /// Density threshold of the solves (`F+0x104`).
    pub fn dthr(&self) -> usize {
        self.dthr
    }

    /// Whether the most recent `ftran` returned a dense result (the reference's count `-1`).
    /// Work of the factor update routine (`0xf15c0`) when it installs this fresh factor, as `(site, estimate)`. The routine
    /// keeps a running sum that it stores at several points; a store it skips leaves its amount to the next one.
    pub fn fresh_work(&self) -> Vec<(u32, f64)> {
        let m = self.m;
        let (dc, x100) = (self.i21 - self.nuc.len(), self.i21);
        let mut out = Vec::new();
        let mut acc = 9.0 * m as f64;
        if dc > 0 {
            out.push((0xf4f9b, acc));
            acc = 0.0;
        }
        acc += dc as f64;
        let r11 = if dc < x100 {
            out.push((0xf4ffd, acc));
            acc = 0.0;
            x100
        } else {
            dc
        };
        acc += 2.0 * r11 as f64;
        if r11 < m {
            out.push((0xf507a, acc));
            acc = 0.0;
        }
        // The highest tail slot whose U column is empty, scanning down to the tail start.
        let mut r13 = m - 1;
        while r13 > x100 && !self.uc[self.at[r13]].is_empty() {
            r13 -= 1;
        }
        acc += (m - 1 - r13) as f64;
        out.push((0xf51af, acc));
        out
    }

    /// Work of the update routine clearing the pending FTRAN spike when the factor is replaced (`0xf15c0`, site `0xf461f`).
    pub fn discard_work(&self) -> f64 {
        2.0 * self.spike.len() as f64
    }

    /// Work of discarding the pending FTRAN spike when the clean-up pass factorizes after the last iteration: the spike counts only if no
    /// column replacement has used it (a bound flip was the last FTRAN).
    pub fn discard_work_unused(&self) -> f64 {
        if self.spike_used {
            0.0
        } else {
            2.0 * self.spike.len() as f64
        }
    }

    /// `nnz(U)` of the fresh factor (the reference's `F+0xb0` pointer array total): the factor size less `nnz(L)` and the two diagonals.
    pub fn u_nnz(&self) -> usize {
        self.znz as usize - self.nnz_l - 2 * self.m
    }

    /// `nnz(L)` of the fresh factor (the reference's `F+0x58` pointer array total).
    pub fn l_nnz(&self) -> usize {
        self.nnz_l
    }

    /// Whether the fresh factor's L solve goes level by level (`0x12a990`).
    pub fn has_levels(&self) -> bool {
        self.lev_built
    }

    /// Dense x after the L stage and after the row etas, for the `lstage` dump.
    pub fn lstage_snap(&self) -> (Vec<R>, Vec<R>) {
        self.lsnap.borrow().clone()
    }

    /// U-stage state of the most recent solve for the `ustage` dump: the spike (the L-stage result), the diagonal, the tail
    /// columns' U entries and the slots.
    pub fn ustage_state(
        &self,
    ) -> (
        &[(usize, R)],
        &[R],
        &[Vec<(usize, R)>],
        &[usize],
        &[i64],
        usize,
    ) {
        (
            &self.spike,
            &self.diag,
            &self.uc,
            &self.slot,
            &self.kind,
            self.i21,
        )
    }

    pub fn last_dense(&self) -> bool {
        self.last_dense
    }

    pub fn last_ordered(&self) -> bool {
        self.last_ordered
    }

    /// Work estimate of an FTRAN with a dense input (the reference's count `-1`, e.g. the basic solution after a refactorization):
    /// dense L solve, all row etas, dense U solve. `b` is indexed by row.
    pub fn ftran_dense_work(&self, b: &[R]) -> f64 {
        let m = self.m;
        let mut x = vec![R::ZERO; m];
        for r in 0..m {
            x[self.perm[r]] += b[r];
        }
        let mut e = (m - self.loff) as f64 * 2.0 + self.n1 as f64 + self.n1 as f64;
        for v in 0..m {
            let xv = x[v];
            if xv.abs() <= R::c(TOL) {
                continue;
            }
            e = self.lcolumn(v).len() as f64 * 2.0 + e;
            for &(y, a) in self.lcolumn(v) {
                x[y] -= a * xv;
            }
        }
        if !self.etas.is_empty() {
            e = self.sum_eta * 2.0 + e;
            for eta in &self.etas {
                let xs: Vec<R> = eta.ind.iter().map(|&j| x[j]).collect();
                let sv = dot(&xs, &eta.val, eta.peel);
                x[eta.idx] = if sv.abs() <= R::c(TOL) { R::ZERO } else { sv };
            }
        }
        for t in (self.i21..m).rev() {
            let c = self.at[t];
            let xc = x[c] / self.diag[c];
            if x[c].abs() <= R::c(TOL) || xc.abs() <= R::c(TOL) {
                x[c] = R::ZERO;
                continue;
            }
            x[c] = xc;
            for &(j, v) in &self.uc[c] {
                x[j] -= v * xc;
            }
        }
        self.u_work(&[], &x, false, None) + e
    }

    /// Work estimate of the U stage: the sparse solve `0x1223e0` (`sparse`) or the dense one `0x121ef0`. `x` holds the solved
    /// values, so a tail column survived the tolerance test when its entry is nonzero. `cnt` is the input count when the input
    /// is a sparse list.
    fn u_work(&self, inp: &[usize], x: &[R], sparse: bool, cnt: Option<usize>) -> f64 {
        let m = self.m;
        let i21 = self.i21;
        let i3 = self.nuc.len();
        let len = |c: usize| self.uc[c].len();
        // Every tail slot from the top down, then the nucleus list (shared by the dense paths).
        let full_tail = |mut e: f64| {
            e = (m - i21) as f64 * 3.0 + e;
            for t in i21..m {
                let c = self.at[t];
                if x[c] != R::ZERO {
                    e += len(c) as f64 * 2.0;
                }
            }
            e
        };
        if !sparse {
            let mut e = 0.0f64;
            if let Some(n) = cnt {
                e = m as f64 + n as f64 * 2.0 + e;
            }
            e = (m - i21) as f64 * 3.0 + e;
            for t in (i21..m).rev() {
                let c = self.at[t];
                if x[c] != R::ZERO {
                    e = (len(c) as f64 + 2.0) * 2.0 + e;
                }
            }
            return 3.0 * i21.saturating_sub(1) as f64 + e;
        }
        let n = inp.len();
        let d2 = self.dens;
        if (m as f64) <= n as f64 * d2 {
            let e = full_tail(0.0);
            return i3 as f64 + i3 as f64 + e;
        }
        let tail: Vec<usize> = inp
            .iter()
            .copied()
            .filter(|&p| self.kind[p] == -2)
            .collect();
        let nuc_in = inp.iter().filter(|&&p| self.kind[p] >= 0).count();
        if d2 * ((nuc_in + tail.len()) as f64) > (i3 + (m - i21)) as f64 {
            let e = (2 * (nuc_in + tail.len())) as f64;
            let e = full_tail(e);
            return i3 as f64 + i3 as f64 + e;
        }
        let mut e = (tail.len() * 2) as f64;
        let mut mark = vec![false; m];
        let mut heap: BinaryHeap<usize> = BinaryHeap::new();
        for &p in &tail {
            mark[p] = true;
            heap.push(self.slot[p]);
        }
        let (mut lb, mut ld) = (nuc_in, nuc_in);
        let mut surv: Vec<usize> = Vec::new();
        let visit = |c: usize, lb: &mut usize, ld: &mut usize, surv: &mut Vec<usize>| -> bool {
            if x[c] != R::ZERO {
                *lb += self.b[c].len();
                *ld += self.a[c].len();
                surv.push(c);
                true
            } else {
                false
            }
        };
        while let Some(&top) = heap.peek() {
            let size = heap.len();
            if size > 1 && ((top - i21) as f64) < 6.0 * size as f64 {
                // Remaining slots are close together: scan the range from the top slot down.
                e = (top - i21) as f64 * 5.0 + size as f64 + size as f64 + e;
                for t in (i21..=top).rev() {
                    visit(self.at[t], &mut lb, &mut ld, &mut surv);
                }
                break;
            }
            heap.pop();
            let c = self.at[top];
            if visit(c, &mut lb, &mut ld, &mut surv) {
                for &j in &self.a[c] {
                    if !mark[j] {
                        mark[j] = true;
                        heap.push(self.slot[j]);
                    }
                }
            }
            e += 8.0;
        }
        e = ld as f64 * 3.0 + e;
        if lb as f64 * 2.5 <= i3 as f64 {
            let mut nm = vec![false; m];
            let mut cntn = 0usize;
            for &p in inp.iter().filter(|&&p| self.kind[p] >= 0) {
                if !nm[p] {
                    nm[p] = true;
                    cntn += 1;
                }
            }
            for &c in &surv {
                e = e + self.b[c].len() as f64 * 2.0 + 4.0;
                for &y in &self.b[c] {
                    if !nm[y] {
                        nm[y] = true;
                        cntn += 1;
                    }
                }
            }
            e + 3.0 * cntn as f64
        } else {
            e += nuc_in as f64;
            for &c in &surv {
                e = e + self.b[c].len() as f64 * 2.0 + 4.0;
            }
            e + 3.0 * i3 as f64
        }
    }

    /// Descending slot order without repeatedly reading indirect keys during comparisons. `slot` and `at` are inverse
    /// permutations, so sorting slots and mapping back also preserves repeated positions exactly.
    fn sort_slot_desc(&self, positions: &mut [usize]) {
        if positions.len() >= 64 && positions.len() > self.m / 8 {
            let mut counts = self.slot_counts.borrow_mut();
            counts.resize(self.m, 0);
            counts.fill(0);
            for &p in positions.iter() {
                counts[self.slot[p]] += 1;
            }
            let mut i = 0;
            for slot in (0..self.m).rev() {
                for _ in 0..counts[slot] {
                    positions[i] = self.at[slot];
                    i += 1;
                }
            }
            return;
        }
        for p in positions.iter_mut() {
            *p = self.slot[*p];
        }
        positions.sort_unstable_by(|a, b| b.cmp(a));
        for p in positions {
            *p = self.at[*p];
        }
    }

    /// Order of the nonzero positions of the solution, given the L-stage list.
    pub fn uorder(&self, inp: &[usize], x: &[R]) -> Vec<usize> {
        let m = self.m;
        let mut seen = vec![false; m];
        let mut inp_d: Vec<usize> = Vec::with_capacity(inp.len());
        for &x in inp {
            if !seen[x] {
                seen[x] = true;
                inp_d.push(x);
            }
        }
        let n = inp_d.len();
        let sing: Vec<usize> = inp_d
            .iter()
            .copied()
            .filter(|&x| self.kind[x] == -1)
            .collect();
        let tri: Vec<usize> = inp_d
            .iter()
            .copied()
            .filter(|&x| self.kind[x] == -2)
            .collect();
        let nuc: Vec<usize> = inp_d
            .iter()
            .copied()
            .filter(|&x| self.kind[x] >= 0)
            .collect();
        let i3 = self.nuc.len();
        let d2 = self.dens;
        let by_slot_desc = |v: &mut Vec<usize>| self.sort_slot_desc(v);
        let nuc_sorted = |set: Vec<usize>| {
            let mut v = set;
            v.sort_unstable_by_key(|&x| self.kind[x]);
            v
        };
        let mut have = vec![false; m];
        let dense = (m as f64) <= (n as f64) * d2
            || d2 * ((nuc.len() + tri.len()) as f64) > (i3 + (m - self.i21)) as f64;
        if dense {
            let mut t: Vec<usize> = tri.clone();
            for &x in &tri {
                have[x] = true;
            }
            let mut k = 0;
            while k < t.len() {
                let c = t[k];
                k += 1;
                for &j in &self.a[c] {
                    if !have[j] {
                        have[j] = true;
                        t.push(j);
                    }
                }
            }
            let mut nm = vec![false; m];
            let mut ns: Vec<usize> = Vec::new();
            for &x in nuc.iter().chain(t.iter().flat_map(|&c| self.b[c].iter())) {
                if !nm[x] {
                    nm[x] = true;
                    ns.push(x);
                }
            }
            by_slot_desc(&mut t);
            let mut out = sing;
            out.extend(t);
            out.extend(nuc_sorted(ns));
            return out;
        }
        let mut heap: BinaryHeap<(usize, usize)> = BinaryHeap::new();
        for &x in &tri {
            have[x] = true;
            heap.push((self.slot[x], x));
        }
        let mut order: Vec<usize> = Vec::new();
        let mut lb = nuc.len();
        while let Some((_, k)) = heap.pop() {
            order.push(k);
            // The reference's counter and nucleus-list discovery only see tail columns that survive the tolerance test.
            if x[k] != R::ZERO {
                lb += self.b[k].len();
            }
            for &j in &self.a[k] {
                if !have[j] {
                    have[j] = true;
                    heap.push((self.slot[j], j));
                }
            }
        }
        let mut nm = vec![false; m];
        let mut lst: Vec<usize> = Vec::new();
        for &y in nuc.iter().chain(
            order
                .iter()
                .filter(|&&k| x[k] != R::ZERO)
                .flat_map(|&k| self.b[k].iter()),
        ) {
            if !nm[y] {
                nm[y] = true;
                lst.push(y);
            }
        }
        let mut out = sing;
        out.extend(order);
        if (lb as f64) * 2.5 <= i3 as f64 {
            out.extend(lst);
        } else {
            out.extend(nuc_sorted(lst));
        }
        out
    }

    fn drop_nuc(&mut self, x: usize) {
        let idx = self.kind[x] as usize;
        let last = self.nuc.pop().unwrap();
        if last != x {
            self.nuc[idx] = last;
            self.kind[last] = idx as i64;
        }
    }

    /// Swap `r` (sitting in slot `hole`) up the dependency walk starting at slot `cur` (the loop at module offset 0xf2f7d of the
    /// reference update); when a column's list is as long as the distance to the hole the rest of the range shifts down instead.
    fn walk(&mut self, r: usize, hole: usize, mut cur: usize) {
        // Column lengths leave out row r, which the reference has already removed from U.
        let len = |s: &Self, pos: usize| {
            s.a[pos]
                .iter()
                .chain(&s.b[pos])
                .filter(|&&x| x != r)
                .count()
        };
        // The reference's work: each jump scans the column (2 per entry and 5) on the factor's scale; the final shift
        // moves the slots in between.
        let mut jumps = 0.0;
        while cur != hole {
            let pos = self.at[cur];
            if len(self, pos) >= cur - hole {
                self.upd_work
                    .push((0xf3d83, 2.0 * (cur - hole) as f64, crate::work::CTX));
                let occ = self.at[hole];
                for t in hole..cur {
                    self.at[t] = self.at[t + 1];
                    let p = self.at[t];
                    self.slot[p] = t;
                }
                self.at[cur] = occ;
                self.slot[occ] = cur;
                break;
            }
            jumps += 2.0 * len(self, pos) as f64 + 5.0;
            let mut nxt = hole;
            for &x in &self.a[pos] {
                let sx = self.slot[x];
                if sx > hole && sx < cur && sx > nxt {
                    nxt = sx;
                }
            }
            let occ = self.at[hole];
            self.at[hole] = pos;
            self.at[cur] = occ;
            self.slot[pos] = hole;
            self.slot[occ] = cur;
            cur = nxt;
        }
        if jumps > 0.0 {
            // Stored before the shift.
            let at = self.upd_work.len()
                - usize::from(self.upd_work.last().is_some_and(|w| w.0 == 0xf3d83));
            self.upd_work.insert(at, (0xf3ca4, jumps, crate::work::FAC));
        }
    }

    /// Row elimination of a column replacement at `r`: the eta multipliers `(position, multiplier)` in elimination order.
    fn eliminate(&self, r: usize, w: &mut [R], seen: &mut [bool]) -> (Vec<(usize, R)>, usize) {
        let m = self.m;
        let d0 = self.diag[r];
        // Tail entries of the maintained transpose, in the same ascending column order as the CSR rebuild.
        let rows = |i: usize| self.ut[i].iter().filter(|&&(k, _)| self.kind[k] == -2);
        let mut touched = Vec::new();
        let mut heap: std::collections::BinaryHeap<std::cmp::Reverse<(usize, usize)>> =
            std::collections::BinaryHeap::new();
        for &(k, v) in rows(r) {
            w[k] = v * (R::ONE / d0);
            touched.push(k);
            seen[k] = true;
            heap.push(std::cmp::Reverse((self.slot[k], k)));
        }
        let mut eta = Vec::new();
        // The reference's solve (`0x131760`) books 2 per entry of the row-wise copy of every accepted position, `r` first.
        // Count only the rows reached by the solve. The transpose is sorted by column,
        // so adjacent duplicates count once, just as in the full column/stamp rebuild.
        let rlen = |y: usize| {
            let mut prev = usize::MAX;
            let mut count = 0;
            for &(c, _) in &self.ut[y] {
                if c != prev && self.kind[c] == -2 {
                    count += 1;
                }
                prev = c;
            }
            count
        };
        // Only positions in the tail (slot >= i21) are visited in the heap phase of the solve that books; the earlier ones are
        // handled in the list phase.
        let mut visited = if self.slot[r] >= self.i21 { rlen(r) } else { 0 };
        // Dense branch (`F+0x108 > 5`): `3` for the one input entry, `2 len` for `r` (a tail column is handled by the slot walk
        // like the rest), `2 (m - minslot)` for the slot walk over every position from the lowest slot reached, and `2 len` for
        // every accepted position of the walk.
        let mut dense = 3.0 + 2.0 * rlen(r) as f64;
        let mut minslot = if self.kind[r] == -2 { self.slot[r] } else { m };
        for &(k, _) in rows(r) {
            minslot = minslot.min(self.slot[k]);
        }
        while let Some(std::cmp::Reverse((_, k))) = heap.pop() {
            let c = -w[k] / self.diag[k];
            if c.abs() <= R::c(TOL) {
                continue;
            }
            eta.push((k, c));
            let len = rlen(k);
            dense += 2.0 * len as f64;
            if self.slot[k] >= self.i21 {
                visited += len;
            }
            for &(l, v) in rows(k) {
                w[l] += c * v;
                if !seen[l] {
                    touched.push(l);
                    seen[l] = true;
                    heap.push(std::cmp::Reverse((self.slot[l], l)));
                }
            }
        }
        dense += 2.0 * (m - minslot) as f64;
        // The solve returns at once, booking nothing, when row `r` of the row-wise copy is empty.
        self.dense_est
            .set(if rows(r).next().is_none() { 0.0 } else { dense });
        for i in touched {
            w[i] = R::ZERO;
            seen[i] = false;
        }
        (eta, visited)
    }

    /// Replace a numeric U column and update its ordered row copy without doing any arithmetic.
    fn replace_u_column(&mut self, c: usize, col: Vec<(usize, R)>) {
        let old = std::mem::replace(&mut self.uc[c], col);
        for &(row, _) in &old {
            self.ut[row].retain(|e| e.0 != c);
        }
        for &(row, v) in &self.uc[c] {
            let pos = self.ut[row].partition_point(|e| e.0 <= c);
            self.ut[row].insert(pos, (c, v));
        }
    }

    /// Replace the column at position `r` by the entering column of the last `ftran`.
    pub fn update(&mut self, r: usize) {
        let m = self.m;
        if crate::dump::on("fstruct") {
            static CALLS: std::sync::atomic::AtomicUsize = std::sync::atomic::AtomicUsize::new(0);
            let n = CALLS.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
            crate::dump::emit(
                "fstruct",
                serde_json::json!({"r": r, "rowl": (0..m).filter(|&y| !self.rowl[y].is_empty()).map(|y| (y.to_string(), self.rowl[y].clone())).collect::<std::collections::BTreeMap<_, _>>(), "mark": self.mark, "n": n, "m": m, "i3": self.nuc.len(), "i21": self.i21, "kind": self.kind, "nuc": self.nuc, "seqA": self.slot, "fill": self.fill, "diag": jv(&self.diag), "uvals": (0..m).filter(|&k| self.kind[k] == -2).map(|k| (k.to_string(), (self.uc[k].iter().filter(|e| self.kind[e.0] == -2).map(|e| e.0).collect::<Vec<_>>(), self.uc[k].iter().filter(|e| self.kind[e.0] == -2).map(|e| e.1.to_f64()).collect::<Vec<_>>()))).collect::<std::collections::BTreeMap<_, _>>(), "cols": (0..m).filter(|&k| self.kind[k] == -2).map(|k| (k.to_string(), [&self.a[k], &self.b[k]])).collect::<std::collections::BTreeMap<_, _>>()}),
            );
        }
        self.before_last = (self.upd, self.sum_alpha);
        self.upd += 1;
        self.upd_work.clear();
        // The eta is stored first, at the current pool fill; its value array sits on a 16-byte boundary when the fill is even.
        let peel = self.fill % 2;
        // The old column's entries leave the row-wise copy after the new column's entries have been added.
        let old_rows: Vec<usize> = if self.kind[r] == -2 {
            self.uc[r].iter().map(|e| e.0).collect()
        } else {
            Vec::new()
        };
        // The new U column's entries (`0xf15c0`, sites `0xf2521` and `0xf2e3c`).
        let k = self.spike.iter().filter(|e| e.0 != r).count();
        if k > 0 {
            self.upd_work
                .push((0xf2521, 3.0 * k as f64, crate::work::CTX));
        }
        // A tail column leaves the row-wise copy of U entry by entry.
        if self.kind[r] == -2 && !self.uc[r].is_empty() {
            let len = self.uc[r].len();
            self.upd_work.push((
                if k == 0 {
                    0xf2521
                } else if len < k {
                    0xf2044
                } else {
                    0xf207a
                },
                5.0 * len as f64,
                crate::work::CTX,
            ));
        }
        self.sum_alpha += self.last_alpha as f64;
        let mut values = std::mem::take(&mut self.update_values);
        let mut seen = std::mem::take(&mut self.update_seen);
        values.resize(m, R::ZERO);
        seen.resize(m, false);
        let (eta_cand, solve_len) = self.eliminate(r, &mut values, &mut seen);
        // The solve runs only when a spike entry lies below `max(slot[r], i21 - 1)` (`0xf1859..0xf191f`).
        let solve_floor = self.slot[r].max(self.i21.saturating_sub(1));
        let cached = self.row_cached.take() == Some(r);
        if !cached
            && (solve_len > 0 || (self.h33 > 5.0 && self.dense_est.get() > 0.0))
            && self
                .spike
                .iter()
                .any(|e| e.0 != r && self.slot[e.0] > solve_floor)
        {
            let est = if self.h33 > 5.0 {
                self.dense_est.get()
            } else {
                2.0 * solve_len as f64
            };
            self.upd_work.push((0x1322ed, est, crate::work::FAC));
        }
        let mut sp = Vec::with_capacity(self.spike.len());
        for &(p, v) in &self.spike {
            values[p] = v;
            if !seen[p] {
                seen[p] = true;
                sp.push(p);
            }
        }
        let d0 = self.diag[r];
        if let Some(i) = sp.iter().position(|&x| x == r) {
            // The reference keeps the spike as a tail group and a nucleus group; the leaving column is deleted from the nucleus
            // group by moving that group's last entry into its place.
            let last = sp.pop().unwrap();
            if last != r {
                sp[i] = last;
            }
        }
        let old = self.kind[r];
        let mut single_rows: Vec<usize> = Vec::new();
        if old == -2 {
            // The old column leaves the rows it appears in (`0xf1e7a`): the entries before it in each row list are passed.
            let rows: Vec<usize> = self.a[r].iter().chain(self.b[r].iter()).copied().collect();
            for y in rows {
                if let Some(p) = self.rowl[y].iter().position(|&t| t == r) {
                    if self.rowl[y].len() > 1 && p > 0 {
                        self.upd_work.push((0xf1f82, p as f64, crate::work::CTX));
                    }
                    if self.rowl[y].len() == 1 {
                        single_rows.push(y);
                    }
                    self.rowl[y].swap_remove(p);
                }
            }
        }
        // A new column without entries besides its diagonal is not a tail column: it stays among the non-tail columns (at the
        // last non-tail slot) and keeps its place in the nucleus lists that refer to it.
        let bare = k == 0 && old != -2;
        let mut demoted_front = Vec::new();
        let mut convslot: Option<usize> = None;
        let mut conv: Option<usize> = None;
        // The reference demotes a bare tail column only when the new column depends on the tail; with no tail entries in the spike
        // the tail simply grows by one slot.
        if old != -2
            && self.slot[r] + 1 != self.i21
            && sp.iter().any(|&y| y != r && self.kind[y] == -2)
        {
            conv = self
                .empty_tail
                .iter()
                .copied()
                .max_by_key(|&x| self.slot[x]);
            if let Some(x) = conv {
                self.empty_tail.remove(&x);
                convslot = Some(self.slot[x]);
                // Only a column some tail column refers to joins the nucleus list; the others become singletons.
                let mut refs: Vec<usize> = self.rowl[x]
                    .iter()
                    .copied()
                    .filter(|&k| self.kind[k] == -2 && self.a[k].contains(&x))
                    .collect();
                refs.sort_unstable();
                refs.dedup();
                let held = !refs.is_empty();
                if held {
                    self.kind[x] = self.nuc.len() as i64;
                    self.nuc.push(x);
                } else {
                    self.kind[x] = -1;
                }
                for k in refs {
                    if self.kind[k] == -2 && self.a[k].contains(&x) {
                        // The last entry of A fills the hole, and the demoted column takes the slot A gives up (the front of B).
                        let p = self.a[k].iter().position(|&y| y == x).unwrap();
                        self.a[k].swap_remove(p);
                        if p > 0 {
                            self.upd_work.push((0xf3eac, p as f64, crate::work::CTX));
                        }
                        self.b[k].insert(0, x);
                        demoted_front.push(k);
                    }
                }
            }
        }
        let mut xlim: i64;
        {
            // The reference scans the slots from the mark down to the first tail column without entries (`0xf1c60`), 3 units per
            // slot passed, when a leaving nucleus column's spike reaches above the last nucleus slot.
            let oldslot0 = self.slot[r] as i64;
            let bt = self.i21 as i64 - 1;
            let floor = bt.max(oldslot0);
            let top = self
                .spike
                .iter()
                .filter(|e| e.0 != r)
                .map(|e| self.slot[e.0] as i64)
                .max()
                .unwrap_or(floor)
                .max(floor);
            xlim = bt;
            if top == floor {
                if bt >= oldslot0 && k > 0 {
                    xlim = bt - 1;
                }
            } else if floor <= oldslot0 {
                if oldslot0 == bt {
                    xlim = bt - 1 + i64::from(k == 0);
                }
            } else {
                let stop = convslot.map_or(floor, |c| c as i64);
                let passed = (self.mark as i64 - stop).max(0);
                if passed > 0 {
                    self.upd_work
                        .push((0xf1d02, 3.0 * passed as f64, crate::work::CTX));
                }
                self.mark = stop as usize;
                if convslot.is_none() {
                    xlim = bt - 1;
                }
            }
        }
        if old >= 0 {
            let held = self.rowl[r]
                .iter()
                .any(|&x| x != r && self.kind[x] == -2 && self.b[x].contains(&r));
            if bare && crate::dump::on("bare") {
                let tail = (0..m).filter(|&x| self.kind[x] == -2).count();
                crate::dump::emit(
                    "bare",
                    serde_json::json!({"r": r, "old": old, "held": held, "tail": tail}),
                );
            }
            if !bare || !held {
                // GUROBI: a bare `old >= 0` replacement whose `r` no tail column's B list holds (`!held`) is never reached; `0xf15c0` has no code that
                // removes a bare `r` from the nucleus list (count `[rbx+0xec]` decremented only at `0xf1eb7`, `0xf2a7f`), so the reference likely keeps it;
                // no witness (`specs/lu-factor.md`, "Bare replacement"; dump kind `bare`, all held).
                self.drop_nuc(r);
                if bare {
                    self.kind[r] = -1;
                }
            }
        } else if old == -2 {
            let old_b = std::mem::take(&mut self.b[r]);
            self.a[r].clear();
            self.replace_u_column(r, Vec::new());
            self.kind[r] = -1;
            for x in old_b {
                if self.kind[x] >= 0
                    && !self.rowl[x]
                        .iter()
                        .any(|&k| self.kind[k] == -2 && self.b[k].contains(&x))
                {
                    self.drop_nuc(x);
                    self.kind[x] = -1;
                }
            }
        }
        for &x in &sp {
            if self.kind[x] == -1 {
                self.kind[x] = self.nuc.len() as i64;
                self.nuc.push(x);
            }
        }
        let av: Vec<usize> = sp.iter().copied().filter(|&x| self.kind[x] == -2).collect();
        let bv: Vec<usize> = sp.iter().copied().filter(|&x| self.kind[x] >= 0).collect();
        let mut refcols = self.rowl[r].clone();
        refcols.sort_unstable();
        refcols.dedup();
        refcols.retain(|&x| {
            x != r && self.kind[x] == -2 && (self.a[x].contains(&r) || self.b[x].contains(&r))
        });
        // Only a surviving row-r reference can raise the mark after the slot walk.
        let refcols: Vec<_> = refcols.into_iter().map(|k| (k, self.slot[k])).collect();
        let oldslot = self.slot[r];
        let top: i64 = av.iter().map(|&x| self.slot[x] as i64).max().unwrap_or(-1);
        let mut moved = true;
        let mut tail_rebuild = false;
        if old == -2 {
            if top >= oldslot as i64 {
                self.walk(r, oldslot, top as usize);
            } else {
                moved = false;
            }
        } else if let Some(cs) = convslot {
            // the demoted column fills r's old slot; r takes the vacated tail slot and walks up from the highest spike entry
            let x = conv.unwrap();
            self.slot[x] = oldslot;
            self.at[oldslot] = x;
            self.slot[r] = cs;
            self.at[cs] = r;
            if top > cs as i64 {
                self.walk(r, cs, top as usize);
            }
        } else {
            // the tail grows by one slot: the occupant of the last non-tail slot fills r's old slot, r takes the new bottom slot
            let bt = self.i21 - 1;
            let q = self.at[bt];
            if q != r {
                self.slot[q] = oldslot;
                self.at[oldslot] = q;
            }
            self.slot[r] = bt;
            self.at[bt] = r;
            if !bare {
                self.i21 = bt;
            }
            if !av.is_empty() {
                self.walk(r, bt, top as usize);
            }
        }
        let _ = &mut tail_rebuild;
        // A row whose only entry was the old column leaves the singleton list when its position is within the reference's limit.
        // GUROBI: the reference removes the row-single nucleus column from the old tail column's lists only if `slot[y] <= [rsp+0x168]` (`0xf1ea5`);
        // rose removes it unconditionally above; no witness of the gated-off case (1693 corpus drops, `specs/lu-factor.md`, "latent gate").
        for y in single_rows {
            if xlim >= self.slot[y] as i64 {
                self.upd_work.push((0xf1f37, 3.0, crate::work::CTX));
            }
        }
        if !bare {
            self.kind[r] = -2;
        }
        // A demotion (`conv`) happens only when the spike references the tail, so the new column is not kept in place even when the demoted column was
        // the spike's only tail entry (FIT2P call 2263: column 175, the reference creates the row eta of the leaving column).
        let keep = (old != -2 && av.is_empty() && conv.is_none()) || (old == -2 && !moved);
        // Row r of U: the tail columns that reference r (`0xf2a51`), unless the column stays in its slot or the empty spike skips the row solve.
        // The maintained row list avoids scanning unrelated columns. Restore the old scan's ascending order and uniqueness.
        let refs = refcols.len();
        if refs > 0 && !sp.is_empty() && !(keep && !moved) {
            self.upd_work
                .push((0xf2a51, 3.0 * refs as f64, crate::work::CTX));
        }
        let mut elim: Vec<usize> = Vec::new();
        for (k, old_ref_slot) in refcols {
            if k == r || self.kind[k] != -2 {
                continue;
            }
            let in_a = self.a[k].contains(&r);
            let in_b = self.b[k].contains(&r);
            if bare {
                continue;
            }
            if in_a || in_b {
                // The reference scans the lists from their start for `r`; the booking is the number of entries passed.
                let pa = self.a[k].iter().position(|&y| y == r).unwrap_or(0);
                let pb = self.b[k]
                    .iter()
                    .position(|&y| y == r)
                    .unwrap_or(0)
                    .saturating_sub(usize::from(demoted_front.binary_search(&k).is_ok()));
                let (site, units) = if keep {
                    (0xf4073, if old != -2 && in_b { pb } else { 0 })
                } else if in_a {
                    (0xf2969, pa)
                } else {
                    (0xf3f71, pb)
                };
                if units > 0 {
                    self.upd_work.push((site, units as f64, crate::work::CTX));
                }
                if keep {
                    if !in_a {
                        self.a[k].push(r);
                    }
                    remove_front_fill(&mut self.b[k], r);
                } else {
                    if in_a {
                        // A and B share one array in the reference (A first): A's last entry fills the hole and B's last entry
                        // takes the slot A gives up, which is the front of B.
                        let p = self.a[k].iter().position(|&y| y == r).unwrap();
                        self.a[k].swap_remove(p);
                        if let Some(last) = self.b[k].pop() {
                            self.b[k].insert(0, last);
                        }
                    }
                    if in_b {
                        remove_last_fill(&mut self.b[k], r);
                    }
                    self.uc[k].retain(|e| e.0 != r);
                    self.ut[r].retain(|e| e.0 != k);
                    elim.push(k);
                    // A column left without entries raises the mark to its slot (`0xf29d5`).
                    if self.a[k].is_empty() && self.b[k].is_empty() {
                        self.empty_tail.insert(k);
                        self.mark = self.mark.max(old_ref_slot);
                    }
                }
            }
        }
        if !keep {
            self.rowl[r].clear();
        }
        let mut newdiag = values[r];
        if !elim.is_empty() {
            let mut ind = vec![r];
            let mut val = vec![R::ONE / d0];
            for &(k, c) in &eta_cand {
                ind.push(k);
                val.push(c);
            }
            let xs: Vec<R> = ind.iter().map(|&j| values[j]).collect();
            // The diagonal is summed from the scratch eta, whose alignment is fixed per factor object.
            newdiag = dot_diag(&xs, &val, self.scratch_peel);
            self.sum_eta += ind.len() as f64;
            self.upd_work
                .push((0xf1beb, 3.0 * ind.len() as f64, crate::work::CTX));
            let e = self.etas.len();
            for &j in &ind {
                self.first_eta[j] = self.first_eta[j].min(e);
            }
            let esz = ind.len();
            self.etas.push(Eta {
                idx: r,
                ind,
                val,
                peel,
            });
            self.fill += esz;
        }
        for &(p, _) in &self.spike {
            values[p] = R::ZERO;
            seen[p] = false;
        }
        self.update_values = values;
        self.update_seen = seen;
        if k == 0 {
            // A new column without entries is a candidate for the scan (`0xf2e58`).
            self.mark = self.mark.max(oldslot);
        }
        if k > 0 {
            self.upd_work
                .push((0xf2e3c, 5.0 * k as f64, crate::work::CTX));
        }
        // The new U column goes to the pool unless it fits the old column's slot, and each of its entries to a row of the row-wise copy, which moves to the end of the
        // pool with room to grow when it is full.
        let fill_in = self.fill;
        let mut fp_rows: Vec<(usize, usize, usize, usize)> = Vec::new();
        if k > old_rows.len() {
            self.fill += k;
        }
        if k > 0 {
            // A row that already holds an entry of the old column takes the new entry in place.
            // The rows are visited in the spike's order with the leaving column replaced by the last entry (`sp`).
            for &row in sp.iter().filter(|&&y| y != r && !old_rows.contains(&y)) {
                // The transpose includes duplicates in the same multiplicity as the column scan.
                let count = self.ut[row]
                    .iter()
                    .filter(|&&(c, _)| c != r && self.kind[c] == -2)
                    .count();
                if count >= self.row_cap[row] {
                    // The growth is capped by the updates left before the refactorization (`0xf2cd1`: `F+0x188 - F+0x18c`).
                    let left = (self
                        .upd_limit
                        .unwrap_or(((m as f64).powf(0.7) as usize).max(20))
                        + 1)
                    .saturating_sub(self.upd);
                    let extra = (count / 2).max(4).min(left);
                    fp_rows.push((row, count, self.row_cap[row], extra));
                    self.row_par[row] = (self.fill % 2) as u8;
                    self.fill += count + extra;
                    self.row_cap[row] = count + extra;
                }
            }
        }
        if crate::dump::on("fillparts") {
            static FP: std::sync::atomic::AtomicUsize = std::sync::atomic::AtomicUsize::new(0);
            let n = FP.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
            crate::dump::emit(
                "fillparts",
                serde_json::json!({"n": n, "k": k, "old": old_rows.len(), "col": if k > old_rows.len() { k } else { 0 }, "fill_after_eta": fill_in, "rows": fp_rows}),
            );
        }
        if self.bt_occupied && !newdiag.to_f64().is_normal() {
            self.bt_occupied = false;
        }
        self.diag[r] = newdiag;
        self.replace_u_column(r, self.spike.iter().copied().filter(|e| e.0 != r).collect());
        self.a[r] = av;
        let mut bb = bv;
        bb.reverse();
        for &y in self.a[r].iter().chain(bb.iter()) {
            self.rowl[y].push(r);
        }
        self.b[r] = bb;
        if self.kind[r] == -2 && self.a[r].is_empty() && self.b[r].is_empty() {
            self.empty_tail.insert(r);
        } else {
            self.empty_tail.remove(&r);
        }
        self.spike_used = true;
    }
}

/// Delete `r` from a column's nucleus list the way the reference does: the entry at the front of this (reversed) list fills the hole.
fn remove_front_fill(v: &mut Vec<usize>, r: usize) {
    if let Some(p) = v.iter().position(|&y| y == r) {
        if p != 0 {
            v[p] = v[0];
        }
        v.remove(0);
    }
}

fn remove_last_fill(v: &mut Vec<usize>, r: usize) {
    if let Some(p) = v.iter().position(|&y| y == r) {
        v.swap_remove(p);
    }
}

#[cfg(test)]
mod transpose_tests {
    use super::*;

    fn replacements<R: Real>(bits: impl Fn(R) -> (u64, u16)) {
        let mut seed = 0x519de38a7u64;
        let mut next = || {
            seed ^= seed << 13;
            seed ^= seed >> 7;
            seed ^= seed << 17;
            seed
        };
        for m in [1, 7, 31] {
            let cols: Vec<_> = (0..m).map(|i| (true, vec![i], vec![1.0])).collect();
            let st = crate::lu_order::factor_structure_r::<R>(m, &cols);
            let mut f = RefStruct::fresh(m, m, &st, &cols);
            for round in 0..128 {
                let c = next() as usize % m;
                let mut col = Vec::new();
                for i in 0..(round % (2 * m + 1)) {
                    let row = next() as usize % m;
                    let value = match i % 7 {
                        0 => 0.0,
                        1 => -0.0,
                        2 => f64::from_bits(1),
                        3 => f64::MAX,
                        _ => f64::from_bits(next() & 0xffef_ffff_ffff_ffff),
                    };
                    col.push((row, R::c(value)));
                }
                f.replace_u_column(c, col);
                // Rebuild from columns independently, checking duplicate ordering and every stored bit.
                let mut expected = vec![Vec::new(); m];
                for (k, col) in f.uc.iter().enumerate() {
                    for &(row, v) in col {
                        expected[row].push((k, bits(v)));
                    }
                }
                let actual: Vec<Vec<_>> =
                    f.ut.iter()
                        .map(|row| row.iter().map(|&(k, v)| (k, bits(v))).collect())
                        .collect();
                assert_eq!(actual, expected);
                let mut lookup = TauRow::<R>::default();
                for row in 0..m {
                    // Exercise epoch wrap as well as ordinary reuse across different and empty rows.
                    if row == m / 2 {
                        lookup.epoch = usize::MAX;
                    }
                    lookup.prepare(m, &f.ut[row]);
                    for c in 0..m {
                        let expected = f.uc[c].iter().find(|e| e.0 == row).map(|e| bits(e.1));
                        assert_eq!(lookup.get(c).map(&bits), expected);
                    }
                }
            }
        }
    }

    #[test]
    fn ordered_u_transpose_replacements() {
        replacements::<f64>(|v| (v.to_bits(), 0));
        replacements::<crate::f80::F80>(|v| v.repr());
    }
}

#[cfg(test)]
mod dense_tail_tests {
    use super::*;

    fn check_needed<R: Real>(f: &RefStruct<R>, x: &[R], bits: impl Fn(R) -> (u64, u16)) {
        let mut prepared = f.clone();
        prepared.ut = vec![Vec::new(); f.m];
        prepared.rowl = vec![Vec::new(); f.m];
        for c in 0..f.m {
            for &(j, v) in &f.uc[c] {
                prepared.ut[j].push((c, v));
            }
            for &j in f.a[c].iter().chain(&f.b[c]) {
                prepared.rowl[j].push(c);
            }
        }
        let rows: Vec<_> = (0..x.len()).filter(|&i| x[i] != R::ZERO).collect();
        let vals: Vec<_> = rows.iter().map(|&i| x[i]).collect();
        let needed: Vec<_> = (0..f.m).rev().collect();
        for dense in [false, true] {
            let mut old = prepared.clone();
            old.lev_dense_ft = dense;
            let new = old.clone();
            let expected = old.frozen_ftran_needed(&rows, &vals, &needed, dense);
            let actual = new.ftran_needed(&rows, &vals, &needed, dense);
            let encode =
                |v: Vec<(usize, R)>| v.into_iter().map(|(i, v)| (i, bits(v))).collect::<Vec<_>>();
            assert_eq!(encode(actual), encode(expected));
            assert_eq!(new.tau_l.get().to_bits(), old.tau_l.get().to_bits());
            assert_eq!(new.tau_u.get().to_bits(), old.tau_u.get().to_bits());
            let a = new.f2.borrow();
            let b = old.f2.borrow();
            assert_eq!(a.tau_n, b.tau_n);
            assert_eq!(a.tau_raw, b.tau_raw);
            assert_eq!(
                a.x_tau.iter().copied().map(&bits).collect::<Vec<_>>(),
                b.x_tau.iter().copied().map(&bits).collect::<Vec<_>>()
            );
        }
    }

    fn check<R: Real>(f: &RefStruct<R>, x: &[R], bits: impl Fn(R) -> (u64, u16)) {
        check_needed(f, x, &bits);
        let rows: Vec<_> = (0..x.len()).filter(|&i| x[i] != R::ZERO).collect();
        let vals: Vec<_> = rows.iter().map(|&i| x[i]).collect();
        let mut old = f.clone();
        let mut new = f.clone();
        let expected = old.frozen_ftran_with(&rows, &vals, false);
        let actual = new.ftran_with(&rows, &vals, false);
        let encode =
            |v: Vec<(usize, R)>| v.into_iter().map(|(i, v)| (i, bits(v))).collect::<Vec<_>>();
        assert_eq!(encode(actual), encode(expected));
        assert_eq!(new.work.to_bits(), old.work.to_bits());
        assert_eq!(encode(new.spike), encode(old.spike));
        assert_eq!(new.last_alpha, old.last_alpha);
        assert_eq!(new.last_dense, old.last_dense);
        assert_eq!(new.last_ordered, old.last_ordered);
        let a = new.f2.borrow();
        let b = old.f2.borrow();
        assert_eq!(
            a.x_col.iter().copied().map(&bits).collect::<Vec<_>>(),
            b.x_col.iter().copied().map(&bits).collect::<Vec<_>>()
        );
    }

    fn identity<R: Real>(m: usize) -> RefStruct<R> {
        let cols: Vec<_> = (0..m).map(|i| (true, vec![i], vec![1.0])).collect();
        let st = crate::lu_order::factor_structure_r::<R>(m, &cols);
        let mut f = RefStruct::fresh(m, m, &st, &cols);
        f.lev_dense_ft = true;
        f
    }

    fn cases<R: Real>(bits: impl Fn(R) -> (u64, u16) + Copy) {
        let mut seed = 0xdef8514302u64;
        let mut next = || {
            seed ^= seed << 13;
            seed ^= seed >> 7;
            seed ^= seed << 17;
            seed
        };
        let values = [
            0.0,
            -0.0,
            f64::from_bits(1),
            -f64::from_bits(1),
            1e-200,
            -1e-200,
            1e200,
            -1e200,
            1.0,
            -1.0,
            TOL,
            -TOL,
        ];
        // Empty tails, singleton factors, arbitrary slot permutations and cancellation-heavy scatters.
        for m in [1, 2, 7, 33] {
            for round in 0..128 {
                let mut f = identity::<R>(m);
                let i21 = round % (m + 1);
                f.i21 = i21;
                for i in (1..m).rev() {
                    let j = next() as usize % (i + 1);
                    f.at.swap(i, j);
                }
                for t in 0..m {
                    let p = f.at[t];
                    f.slot[p] = t;
                    f.kind[p] = if t >= i21 { -2 } else { -1 };
                }
                for t in i21..m {
                    let c = f.at[t];
                    f.diag[c] = R::c(values[4 + next() as usize % 6]);
                    for u in 0..t {
                        if next() % 3 == 0 {
                            let j = f.at[u];
                            f.uc[c].push((j, R::c(values[next() as usize % values.len()])));
                            if u >= i21 {
                                f.a[c].push(j);
                            }
                        }
                    }
                }
                let x: Vec<_> = (0..m)
                    .map(|_| R::c(values[next() as usize % values.len()]))
                    .collect();
                check(&f, &x, bits);
                check(&f, &vec![R::ZERO; m], bits);
            }
        }
        // SC50A dual-release dense U inputs captured from the parent binary's BSOLVE dump.
        let records: serde_json::Value =
            serde_json::from_str(include_str!("../tests/fixtures/ftran_dense_tail.json")).unwrap();
        for v in records.as_array().unwrap() {
            let numbers = |key: &str| {
                v[key]
                    .as_array()
                    .unwrap()
                    .iter()
                    .map(|x| R::c(x.as_f64().unwrap()))
                    .collect::<Vec<_>>()
            };
            let x = numbers("l");
            let m = x.len();
            let mut f = identity::<R>(m);
            f.diag = numbers("diag");
            f.kind = serde_json::from_value(v["kind"].clone()).unwrap();
            f.slot = serde_json::from_value(v["slot"].clone()).unwrap();
            f.nuc = (0..m).filter(|&p| f.kind[p] >= 0).collect();
            f.nuc.sort_unstable_by_key(|&p| f.kind[p]);
            f.i21 = (0..m)
                .filter(|&p| f.kind[p] == -2)
                .map(|p| f.slot[p])
                .min()
                .unwrap_or(m);
            for p in 0..m {
                f.at[f.slot[p]] = p;
            }
            for c in 0..m {
                f.uc[c] = v["uc"][c]
                    .as_array()
                    .unwrap()
                    .iter()
                    .map(|e| {
                        (
                            e[0].as_u64().unwrap() as usize,
                            R::c(e[1].as_f64().unwrap()),
                        )
                    })
                    .collect();
                f.a[c] = f.uc[c]
                    .iter()
                    .filter(|e| f.kind[e.0] == -2)
                    .map(|e| e.0)
                    .collect();
                f.b[c] = f.uc[c]
                    .iter()
                    .filter(|e| f.kind[e.0] >= 0)
                    .map(|e| e.0)
                    .collect();
            }
            check(&f, &x, bits);
        }
    }

    #[test]
    fn sparse_touches_match_frozen_ftran_across_block_boundaries() {
        let m = 1025;
        let mut f = identity::<f64>(m);
        f.lev_dense_ft = false;
        f.dthr = 1;
        f.i21 = m - 2;
        for p in f.i21..m {
            f.kind[p] = -2;
        }
        // Unsorted touches, repeated indices, cancellations and two tail columns.
        f.a[m - 1] = vec![m - 2];
        f.uc[m - 1] = vec![(512, 2.0), (m - 2, 1.0), (64, 1.0), (512, -2.0)];
        f.uc[m - 2] = vec![(63, 1.0), (0, -1.0)];
        let mut x = vec![0.0; m];
        for i in [1, 63, 64, 511, 960, m - 1] {
            x[i] = 1.0;
        }
        check(&f, &x, |v| (v.to_bits(), 0));
        check(&f, &vec![0.0; m], |v| (v.to_bits(), 0));
    }

    #[test]
    fn column_snapshot_recycles_only_written_blocks() {
        let m = 513;
        let cols: Vec<_> = (0..m).map(|r| (true, vec![r], vec![1.0])).collect();
        let st = crate::lu_order::factor_structure_r::<f64>(m, &cols);
        let mut new = RefStruct::fresh(m, m, &st, &cols);
        let mut old = new.clone();
        let encode =
            |v: &[(usize, f64)]| v.iter().map(|&(p, x)| (p, x.to_bits())).collect::<Vec<_>>();
        for round in 0..12 {
            if round == 5 {
                new = new.clone();
                old = old.clone();
            }
            let rows = [round * 39, 64, 512, round * 39, 0];
            let values = [3.5, -0.0, -(round as f64), -3.5, f64::from_bits(1)];
            let dense = round % 4 == 2;
            new.lev_dense_ft = dense;
            old.lev_dense_ft = dense;
            let a = new.ftran_with(&rows, &values, dense);
            let b = old.frozen_ftran_with(&rows, &values, dense);
            assert_eq!(encode(&a), encode(&b));
            assert_eq!(new.work.to_bits(), old.work.to_bits());
            assert_eq!(encode(&new.spike), encode(&old.spike));
            let nf = new.f2.borrow();
            let of = old.f2.borrow();
            assert_eq!(nf.dirty_blocks.is_none(), dense);
            assert_eq!(nf.col_in, of.col_in);
            assert_eq!(nf.col_n, of.col_n);
            let bits = |v: &[f64]| v.iter().map(|x| x.to_bits()).collect::<Vec<_>>();
            assert_eq!(bits(&nf.x_col), bits(&of.x_col));
        }
    }

    #[test]
    fn needed_dense_tail_divides_unreachable_zero_slots_before_drop() {
        // Minimal kernel witness of native 12403b..12405f: no path from the
        // sparse RHS to the negative diagonal, but the dense slot loop visits it.
        let mut f = identity::<f64>(8);
        f.lev_dense_ft = false;
        f.i21 = 4;
        f.kind[4..].fill(-2);
        f.h33 = 1.0;
        f.h34 = 0.0;
        f.dens = 1.0;
        f.diag[5] = -1.0;
        assert!(f.tau_u_dense(&[(4, 1.0)]));
        let legacy = f.frozen_ftran_needed(&[4], &[1.0], &[5, 6], false);
        assert_eq!(legacy[0].1.to_bits(), 0.0f64.to_bits());
        let actual = f.ftran_needed(&[4], &[1.0], &[5, 6], false);
        assert_eq!(actual[0].1.to_bits(), (-0.0f64).to_bits());
        assert_eq!(actual[1].1.to_bits(), 0.0f64.to_bits());
        let pair = f.f2.borrow();
        assert_eq!(pair.x_tau[5].to_bits(), (-0.0f64).to_bits());
        assert!(!pair.tau_raw[5]);
        assert_eq!(f.tau_fin.get(), 4);
    }

    #[test]
    fn dense_tail_matches_frozen_ftran() {
        cases::<f64>(|v| (v.to_bits(), 0));
        cases::<crate::f80::F80>(|v| v.repr());
    }
}

#[cfg(test)]
mod forward_l_tests {
    use super::*;

    fn factor<R: Real>(m: usize) -> RefStruct<R> {
        let size = m.max(1);
        let cols: Vec<_> = (0..size).map(|i| (true, vec![i], vec![1.0])).collect();
        let st = crate::lu_order::factor_structure_r::<R>(size, &cols);
        let mut f = RefStruct::fresh(size, size, &st, &cols);
        f.m = m;
        f.loff = 0;
        f.n1 = 0;
        f.lev_dense_ft = false;
        f.lstart.truncate(m + 1);
        f
    }

    fn store<R: Real>(f: &mut RefStruct<R>, columns: Vec<Vec<(usize, R)>>) {
        f.lstart.clear();
        f.lentries.clear();
        f.lstart.push(0);
        for col in columns {
            f.lentries.extend(col);
            f.lstart.push(f.lentries.len());
        }
    }

    fn check<R: Real>(
        f: &RefStruct<R>,
        r0: &[usize],
        v0: &[R],
        r1: &[usize],
        v1: &[R],
        bits: impl Fn(R) -> (u64, u16) + Copy,
    ) {
        let encode =
            |v: Vec<(usize, R)>| v.into_iter().map(|(i, v)| (i, bits(v))).collect::<Vec<_>>();
        for dense in [false, true] {
            for (rows, vals) in [(r0, v0), (r1, v1)] {
                let old = f.frozen_lstage_num(rows, vals, dense);
                let new = f.lstage_num(rows, vals, dense);
                assert_eq!(encode(old.0), encode(new.0));
                assert_eq!(old.1, new.1);
                assert_eq!(old.2.to_bits(), new.2.to_bits());
            }
        }
    }

    #[test]
    fn restored_quad_factors_match_independent_installs_after_updates() {
        type Q = crate::f80::F80;
        let make = |m, shift| {
            let cols: Vec<_> = (0..m)
                .map(|i| {
                    (
                        false,
                        vec![i, (i + 1) % m],
                        vec![1.0 + (shift + i) as f64, 0.125],
                    )
                })
                .collect();
            let st = crate::lu_order::factor_structure_r::<Q>(m, &cols);
            RefStruct::fresh(m, m, &st, &cols)
        };
        let bases = [make(7, 1), make(9, 3), make(5, 0)];
        let mut live = bases[0].clone();
        let encode = |v: Vec<(usize, Q)>| {
            v.into_iter()
                .map(|(j, x)| (j, x.repr()))
                .collect::<Vec<_>>()
        };
        for &n in &[0, 0, 1, 2, 1, 0, 2, 2] {
            let base = &bases[n];
            live.restore(base);
            let mut reference = base.clone();
            for pass in 0..3 {
                let r = (2 * pass + 1) % base.m;
                for (rows, v) in [
                    (vec![r], vec![Q::c(3.)]),
                    (
                        vec![r, (r + 2) % base.m, r],
                        vec![Q::c(-0.0), Q::c(-0.125), Q::c(0.5)],
                    ),
                ] {
                    let a = live.ftran(&rows, &v);
                    let b = reference.ftran(&rows, &v);
                    assert_eq!(encode(a), encode(b));
                    assert_eq!(live.work.to_bits(), reference.work.to_bits());
                    let a = live.btran(&[(r, Q::c(-0.75))], true, false);
                    let b = reference.btran(&[(r, Q::c(-0.75))], true, false);
                    assert_eq!(encode(a), encode(b));
                    assert_eq!(live.work.to_bits(), reference.work.to_bits());
                }
                live.update(r);
                reference.update(r);
                assert_eq!(live.upd_work, reference.upd_work);
            }
        }
    }

    #[test]
    fn level_bucket_range_survives_factor_clone() {
        let m = 250;
        let cols: Vec<_> = (0..m)
            .map(|i| {
                if i < 4 {
                    (false, vec![i, i + 1], vec![1.0, 0.5])
                } else {
                    (true, vec![i], vec![1.0])
                }
            })
            .collect();
        let st = crate::lu_order::factor_structure_r::<f64>(m, &cols);
        let f = RefStruct::fresh(m, m, &st, &cols);
        let (levels, max_level) = f.levels.as_ref().unwrap();
        assert_eq!(*max_level, 4);
        assert_eq!(levels.iter().copied().max(), Some(*max_level));
        for factor in [&f, &f.clone()] {
            check(factor, &[0], &[1.0], &[4, 0, 0], &[0.5, -1.0, 1.0], |v| {
                (v.to_bits(), 0)
            });
        }
    }

    #[test]
    fn indexed_level_etas_preserve_order_values_and_work_across_updates() {
        let m = 321;
        let mut f = factor::<f64>(m);
        let mut cols = vec![Vec::new(); m];
        for r in 0..m - 1 {
            if r % 8 < 6 {
                cols[r].push((r + 1, if r % 2 == 0 { 0.5 } else { -0.25 }));
            }
        }
        store(&mut f, cols);
        f.levels = Some(((0..m).map(|i| i % 8).collect(), 7));
        f.dthr = 2;
        f.nnz_l = f.lentries.len();
        let bits = |v: f64| (v.to_bits(), 0);
        for i in 0..137 {
            let r = (37 * i) % m;
            let a = (r + 1) % m;
            let b = (r + 65) % m;
            let ind = vec![r, a, a, b];
            for &j in &ind {
                f.first_eta[j] = f.first_eta[j].min(i);
            }
            f.etas.push(Eta {
                idx: r,
                ind,
                val: vec![1.0, 0.5, -0.5, 0.25],
                peel: i % 2,
            });
            f.sum_eta += 4.0;
            if [0, 1, 63, 64, 65, 136].contains(&i) {
                check(
                    &f,
                    &[0, 64, 0, r],
                    &[1.0, -3.25, 0.75, -0.0],
                    &[r, b, b],
                    &[2.0, 3.0, -3.0],
                    bits,
                );
                check(&f, &[], &[], &[8, 11], &[TOL, f64::from_bits(1)], bits);
                check(&f, &[b], &[f64::INFINITY], &[r], &[f64::NEG_INFINITY], bits);
            }
        }
        let rows: Vec<usize> = (0..161).step_by(2).collect();
        check(
            &f,
            &rows,
            &vec![1.0; rows.len()],
            &[60, 79],
            &[1.0, -1.0],
            bits,
        );
        check(&f.clone(), &[0, 1], &[0.5, -0.125], &[218], &[5.0], bits);
        let i = f.etas.len();
        let ind = vec![32, 100];
        for &j in &ind {
            f.first_eta[j] = f.first_eta[j].min(i);
        }
        f.etas.push(Eta {
            idx: 32,
            ind,
            val: vec![1.0, f64::INFINITY],
            peel: 0,
        });
        f.sum_eta += 2.0;
        check(&f, &[40], &[1.0], &[64], &[3.0], bits);
    }

    #[test]
    fn level_workspace_matches_captured_factors_and_random_update_streams() {
        let captures: serde_json::Value =
            serde_json::from_str(include_str!("../tests/fixtures/ftran_l_captured.json")).unwrap();
        let mut seed = 0x9f2324abcd_u64;
        let mut next = || {
            seed ^= seed << 13;
            seed ^= seed >> 7;
            seed ^= seed << 17;
            seed as usize
        };
        for rec in captures.as_array().unwrap() {
            let m = 321;
            let mut f = factor::<f64>(m);
            let perm: Vec<usize> = serde_json::from_value(rec["perm"].clone()).unwrap();
            f.perm[..perm.len()].copy_from_slice(&perm);
            let destinations: Vec<Vec<usize>> =
                serde_json::from_value(rec["lcols"].clone()).unwrap();
            let coeff: Vec<Vec<f64>> = serde_json::from_value(rec["lvals"].clone()).unwrap();
            let mut cols = vec![Vec::new(); m];
            for (column, (indices, values)) in destinations.into_iter().zip(coeff).enumerate() {
                cols[column] = indices.into_iter().zip(values).collect();
            }
            let mut levels = vec![0; m];
            for (i, column) in cols.iter().enumerate() {
                for &(j, _) in column {
                    levels[j] = levels[j].max(levels[i] + 1);
                }
            }
            let maxlev = levels.iter().copied().max().unwrap();
            f.levels = Some((levels, maxlev));
            store(&mut f, cols);
            f.nnz_l = f.lentries.len();
            f.dthr = 2;
            let captured_rows: Vec<usize> = serde_json::from_value(rec["rows"].clone()).unwrap();
            let captured_vals: Vec<f64> = serde_json::from_value(rec["vals"].clone()).unwrap();
            for e in 0..90 {
                let row = next() % m;
                let mut ind = vec![row];
                for _ in 0..next() % 12 {
                    ind.push(next() % m);
                }
                for &j in &ind {
                    f.first_eta[j] = f.first_eta[j].min(e);
                }
                let val = (0..ind.len())
                    .map(|_| (next() % 17) as f64 * 0.125 - 1.0)
                    .collect();
                f.sum_eta += ind.len() as f64;
                f.etas.push(Eta {
                    idx: row,
                    ind,
                    val,
                    peel: e % 2,
                });
                let rows: Vec<_> = (0..1 + next() % 13).map(|_| next() % m).collect();
                let vals: Vec<_> = rows
                    .iter()
                    .map(|_| (next() % 127) as f64 * 0.25 - 8.0)
                    .collect();
                check(&f, &captured_rows, &captured_vals, &rows, &vals, |v| {
                    (v.to_bits(), 0)
                });
            }
        }
    }

    fn cases<R: Real>(bits: impl Fn(R) -> (u64, u16) + Copy) {
        let mut seed = 0x564af2783u64;
        let mut next = || {
            seed ^= seed << 13;
            seed ^= seed >> 7;
            seed ^= seed << 17;
            seed
        };
        let values = [
            0.0,
            -0.0,
            1.0,
            -1.0,
            0.125,
            -0.25,
            TOL,
            -TOL,
            f64::from_bits(1),
            -f64::MIN_POSITIVE,
            1e-200,
            1e200,
            f64::MAX,
            -f64::MAX,
            f64::INFINITY,
            f64::NEG_INFINITY,
            f64::from_bits(0x7ff8_0000_0000_0123),
        ];
        for m in [0, 1, 2, 7, 31] {
            for round in 0..32 {
                let sample_len = if round < 16 { 8 } else { values.len() };
                let mut f = factor::<R>(m);
                for i in (1..m).rev() {
                    f.perm.swap(i, next() as usize % (i + 1));
                }
                let mut columns = vec![Vec::new(); m];
                for v in 0..m {
                    for y in v + 1..m {
                        if next() % 3 == 0 {
                            columns[v].push((y, R::c(values[next() as usize % sample_len])));
                        }
                    }
                }
                store(&mut f, columns);
                if m > 0 {
                    for k in 0..round % 4 {
                        let ind: Vec<_> = (0..m).filter(|_| next() % 2 == 0).collect();
                        let val: Vec<_> = ind
                            .iter()
                            .map(|_| R::c(values[next() as usize % sample_len]))
                            .collect();
                        f.sum_eta += ind.len() as f64;
                        f.etas.push(Eta {
                            idx: next() as usize % m,
                            ind,
                            val,
                            peel: k % 2,
                        });
                    }
                }
                let r0: Vec<_> = (0..m + round % (m + 1))
                    .map(|_| next() as usize % m)
                    .collect();
                let r1: Vec<_> = (0..round % (2 * m + 1))
                    .map(|_| next() as usize % m)
                    .collect();
                let v0: Vec<_> = r0
                    .iter()
                    .map(|_| R::c(values[next() as usize % sample_len]))
                    .collect();
                let v1: Vec<_> = r1
                    .iter()
                    .map(|_| R::c(values[next() as usize % sample_len]))
                    .collect();
                check(&f, &r0, &v0, &r1, &v1, bits);
            }
        }
        let captured: serde_json::Value =
            serde_json::from_str(include_str!("../tests/fixtures/ftran_l_captured.json")).unwrap();
        for rec in captured.as_array().unwrap() {
            let mut f = factor::<R>(rec["perm"].as_array().unwrap().len());
            f.perm = serde_json::from_value(rec["perm"].clone()).unwrap();
            let rows_l: Vec<Vec<usize>> = serde_json::from_value(rec["lcols"].clone()).unwrap();
            let columns: Vec<Vec<(usize, R)>> = rows_l
                .into_iter()
                .zip(rec["lvals"].as_array().unwrap())
                .map(|(rows, vals)| {
                    rows.into_iter()
                        .zip(vals.as_array().unwrap())
                        .map(|(r, v)| (r, R::c(v.as_f64().unwrap())))
                        .collect()
                })
                .collect();
            store(&mut f, columns);
            let rows: Vec<usize> = serde_json::from_value(rec["rows"].clone()).unwrap();
            let vals: Vec<_> = rec["vals"]
                .as_array()
                .unwrap()
                .iter()
                .map(|v| R::c(v.as_f64().unwrap()))
                .collect();
            // SC50A dual-release BSOLVE arrays from the validated parent, capped at two iterations.
            if !R::QUAD {
                let mut raw = vec![0.0; f.m];
                for (i, v) in f.frozen_lstage_num(&rows, &vals, true).0 {
                    raw[i] = v.to_f64();
                }
                let captured: Vec<_> = rec["l"]
                    .as_array()
                    .unwrap()
                    .iter()
                    .map(|v| v.as_f64().unwrap().to_bits())
                    .collect();
                assert_eq!(
                    raw.iter().map(|v| v.to_bits()).collect::<Vec<_>>(),
                    captured
                );
            }
            check(&f, &rows, &vals, &rows, &vals, bits);
        }
    }

    #[test]
    fn dense_compaction_preserves_bits() {
        let mut seed = 0xa7c9e21356u64;
        for m in [0, 1, 2, 7, 31, 128, 511] {
            for pattern in 0..5 {
                let x: Vec<f64> = (0..m)
                    .map(|i| {
                        seed ^= seed << 13;
                        seed ^= seed >> 7;
                        seed ^= seed << 17;
                        match pattern {
                            0 => 0.0,
                            1 => -0.0,
                            2 if i % 2 == 0 => 0.0,
                            3 => f64::from_bits(0x7ff8_0000_0000_0001 + (i as u64 & 1)),
                            _ => f64::from_bits(seed),
                        }
                    })
                    .collect();
                let expected: Vec<_> = x
                    .iter()
                    .enumerate()
                    .filter(|(_, v)| **v != 0.0)
                    .map(|(i, v)| (i, v.to_bits()))
                    .collect();
                let actual: Vec<_> = RefStruct::<f64>::dense_spike(&x)
                    .into_iter()
                    .map(|(i, v)| (i, v.to_bits()))
                    .collect();
                assert_eq!(actual, expected);
            }
        }
    }

    #[test]
    fn dense_btran_compaction_preserves_row_order_and_bits() {
        let mut seed = 0x29ace71f3bu64;
        for m in [0, 1, 2, 7, 31, 128, 511] {
            for pattern in 0..5 {
                let x: Vec<f64> = (0..m)
                    .map(|i| {
                        seed ^= seed << 13;
                        seed ^= seed >> 7;
                        seed ^= seed << 17;
                        match pattern {
                            0 => 0.0,
                            1 => -0.0,
                            2 if i % 2 == 0 => 0.0,
                            3 => f64::from_bits(0x7ff8_0000_0000_0001 + (i as u64 & 1)),
                            _ => f64::from_bits(seed),
                        }
                    })
                    .collect();
                let mut perm: Vec<_> = (0..m).collect();
                for i in (1..m).rev() {
                    seed ^= seed << 13;
                    seed ^= seed >> 7;
                    seed ^= seed << 17;
                    perm.swap(i, seed as usize % (i + 1));
                }
                let mut expected = Vec::with_capacity(m);
                for r in 0..m {
                    let i = perm[r];
                    if x[i] != 0.0 {
                        expected.push((r, x[i]));
                    }
                }
                let actual = RefStruct::<f64>::dense_btran_output(&x, &perm);
                assert_eq!(actual.capacity(), expected.capacity());
                assert_eq!(
                    actual
                        .iter()
                        .map(|&(i, v)| (i, v.to_bits()))
                        .collect::<Vec<_>>(),
                    expected
                        .iter()
                        .map(|&(i, v)| (i, v.to_bits()))
                        .collect::<Vec<_>>(),
                );
            }
        }
        let x = [0.0, f64::from_bits(0xfff8_0000_0000_0007), -0.0, 3.0];
        let out = RefStruct::<f64>::dense_btran_output(&x, &[3, 1, 3, 2, 0]);
        assert_eq!(
            out.iter()
                .map(|&(r, v)| (r, v.to_bits()))
                .collect::<Vec<_>>(),
            [
                (0, 3.0f64.to_bits()),
                (1, x[1].to_bits()),
                (2, 3.0f64.to_bits())
            ]
        );
    }

    #[test]
    fn forward_columns_match_frozen_l() {
        cases::<f64>(|v| (v.to_bits(), 0));
        cases::<crate::f80::F80>(|v| v.repr());
    }
}
