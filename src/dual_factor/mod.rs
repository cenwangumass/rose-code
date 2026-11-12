//! Dual-only factor backend, initially a faithful copy of the reference structural mirror: the order in which it lists the nonzeros of `B^-1 a` (`reverse-engineer/specs/ftran-order.md`).
//!
//! Positions are the basis positions after the LU renumbering (the pivot step). A fresh factor has singleton positions (kind -1),
//! a nucleus list of diagonal-only positions (phase 3 pivots, kind = index in the list) and a triangular tail (phase 4 and the
//! Markowitz nucleus, slots `i21..m`) whose columns carry U entries split into tail positions (`a`) and nucleus-list positions (`b`).
//! Column replacements follow the update rules recovered from the GDB dumps. The factor is mirrored numerically as well (L multipliers,
//! U columns, diagonal, row etas, including the reference's lane-wise eta dot products), so the solve returns the reference's values
//! bit for bit and drops the entries that cancel exactly as it does.
//!
//! Scaffold copied from `ftran_order.rs` at ef273dec. Numeric kernels, operation
//! order, work estimates and update rules are unchanged; primal-only helpers and
//! frozen test oracles remain in that module, not in this dual backend.

use crate::lu_order::Structure;
use crate::real::Real;
use ab::AbLists;
use pool::{Csr, LPool};
use urows::URows;
use std::collections::BinaryHeap;
mod ab;
mod btran_ws;
mod pool;
mod urows;
#[cfg(test)]
mod urows_tests;
#[cfg(feature = "factor-shadow")]
mod adapter;
#[cfg(test)]
mod ab_tests;
#[cfg(test)]
mod btran_tests;
#[cfg(test)]
mod handoff_tests;
#[cfg(test)]
mod needed_old;
#[cfg(test)]
mod needed_tests;
#[cfg(test)]
thread_local! {
    /// Dense L solves of this test thread that carried a count (the `force_dense` readers' case).
    static LIVE_COUNTS: std::cell::Cell<usize> = const { std::cell::Cell::new(0) };
    /// Routes of `ftran_with_sparse` on this test thread: calls, full sweeps, heap walks that stayed a heap, heap walks that fell back to a
    /// sweep, nucleus lists sorted, nucleus lists in discovery order.
    static SPARSE_ROUTES: std::cell::Cell<[usize; 6]> = const { std::cell::Cell::new([0; 6]) };
}
#[cfg(test)]
mod fresh_tests;
#[cfg(test)]
mod update_old;
#[cfg(test)]
mod full_old;
#[cfg(test)]
mod full_tests;
#[cfg(test)]
mod ordering_old;
#[cfg(test)]
mod ordering_tests;

#[cfg(feature = "factor-shadow")]
pub(crate) mod projection;
#[cfg(feature = "factor-shadow")]
pub mod shadow;

/// The sums `ftran_needed` books while it finishes the tail columns (all integer-valued).
#[derive(Default)]
struct TailCounts {
    /// Survivors (heap-walk booking).
    nz: usize,
    /// Tail-row entries of the survivors above the tolerance (dense booking).
    abs_a: usize,
    /// `2 len + 4` over the raw survivors of the joint dense loops, and the `len` of those whose snapshot entry is zero.
    tp_a2: usize,
    tp_only_a: usize,
    /// Columns with a nonzero final value (checked against the vector in debug builds).
    #[cfg(any(test, debug_assertions, feature = "factor-shadow"))]
    surv_nz: usize,
}

impl TailCounts {
    /// One finished column with the nonzero final value `fin`, `na` entries in its A list and, for the joint loops, its snapshot entry.
    #[inline(always)]
    fn book<const JD: bool, const TD: bool, R: Real>(&mut self, na: usize, fin: R, raw: bool, col: Option<&R>) {
        #[cfg(any(test, debug_assertions, feature = "factor-shadow"))]
        {
            self.surv_nz += 1;
        }
        if TD {
            if fin.abs() > R::c(TOL) {
                self.abs_a += na;
            }
        } else {
            self.nz += 1;
        }
        if JD && raw {
            self.tp_a2 += 2 * na + 4;
            if col.map_or(true, |v| *v == R::ZERO) {
                self.tp_only_a += na;
            }
        }
    }
}

impl<R: Real> DualFactor<R> {
    pub fn work(&self) -> f64 { self.work }
    pub fn work_parts(&self) -> (f64, f64) { self.work_parts }
    pub fn keep_work(&self) -> f64 { self.keep_work }
    pub fn keep_dense(&self) -> bool { self.keep_dense }
    pub fn tau_l(&self) -> f64 { self.tau_l.get() }
    pub fn tau_u(&self) -> f64 { self.tau_u.get() }
    pub fn tau_fin(&self) -> usize { self.tau_fin.get() }
    pub fn bt_dense(&self) -> bool { self.bt_dense }
    pub fn upd_limit(&self) -> Option<usize> { self.upd_limit }
    pub fn set_upd_limit(&mut self, limit: Option<usize>) { self.upd_limit = limit; }
    pub fn clear_update_work(&mut self) { self.upd_work.clear(); }
    pub fn take_update_work(&mut self) -> Vec<(u32, f64, f64)> { std::mem::take(&mut self.upd_work) }
    /// Hands the list taken by `take_update_work` back (emptied) so the next `update` does not allocate one.
    pub fn give_update_work(&mut self, mut list: Vec<(u32, f64, f64)>) {
        list.clear();
        if self.upd_work.capacity() == 0 {
            self.upd_work = list;
        }
    }
}

/// The sparse U heap hands over to a sweep of the remaining slots once `top - loff < SWEEP_RATIO * heap size` (the reference's constant is 6.0).
const SWEEP_RATIO: usize = 6;

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

    /// Marks every row of the list (the result of each `insert` is not needed), the variant tested once.
    fn insert_all(&mut self, rows: &[usize]) {
        match self {
            Self::Dense(flags) => {
                for &row in rows {
                    flags[row] = true;
                }
            }
            Self::Sparse { words, occupied } => {
                for &row in rows {
                    let block = row / 64;
                    let word = &mut words[block];
                    if *word == 0 {
                        occupied.push(block);
                    }
                    *word |= 1u64 << (row % 64);
                }
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

/// Retained result lists `(position, value)`: a solve takes one (empty, capacity kept) instead of allocating, and the consumer hands
/// it back with `DualFactor::give_list`. Only capacities circulate; every list is cleared when taken.
#[derive(Default)]
struct ListPool<R>(std::cell::RefCell<Vec<Vec<(usize, R)>>>);

impl<R> Clone for ListPool<R> {
    fn clone(&self) -> Self {
        Self(Default::default())
    }
}

impl<R> ListPool<R> {
    /// Lists kept at most (a few solves are alive at once).
    const KEEP: usize = 8;

    fn take(&self) -> Vec<(usize, R)> {
        let mut v = self.0.borrow_mut().pop().unwrap_or_default();
        v.clear();
        v
    }

    fn give(&self, mut v: Vec<(usize, R)>) {
        if v.capacity() == 0 {
            return;
        }
        v.clear();
        let mut pool = self.0.borrow_mut();
        if pool.len() < Self::KEEP {
            pool.push(v);
        }
    }
}

/// Retained work buffers of the L stage (`lstage_num`) and of its row etas. Between calls `active` is all false and the buckets are
/// empty; `order` and `eta_x` are cleared where they are taken.
struct LScratch<R> {
    buckets: Vec<Vec<usize>>,
    active: Vec<bool>,
    order: Vec<usize>,
    eta_x: Vec<R>,
}

impl<R> Default for LScratch<R> {
    fn default() -> Self {
        Self { buckets: Vec::new(), active: Vec::new(), order: Vec::new(), eta_x: Vec::new() }
    }
}

impl<R> Clone for LScratch<R> {
    fn clone(&self) -> Self {
        Self::default()
    }
}

/// Per-factor scratch and an index of the append-only row etas. A cloned
/// factor builds its own scratch on demand; no caller can observe the cache.
#[derive(Default)]

struct LevelCache<R>(std::cell::RefCell<Option<Box<LevelWorkspace<R>>>>);

impl<R: Default> Clone for LevelCache<R> {
    fn clone(&self) -> Self {
        Self::default()
    }
}

#[derive(Default)]

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
    /// Forget the state of a retired factor (its etas, rows and scratch values), keeping every allocation.
    fn reset(&mut self) {
        self.values.clear();
        self.state.clear();
        self.touched.clear();
        self.buckets.iter_mut().for_each(Vec::clear);
        self.order.clear();
        self.uses.iter_mut().for_each(Vec::clear);
        self.term_ends.clear();
        self.pending.clear();
        self.gather.clear();
        self.incompatible = false;
    }

    fn prepare(&mut self, f: &DualFactor<R>, nlevels: usize) -> bool {
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
        f: &DualFactor<R>,
        levels: &[usize],
        maxlev: usize,
        rows: &[usize],
        vals: &[R],
        pooled: bool,
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
                    let (dest, mult) = f.lcolumn(row);
                    for (&to, &a) in dest.iter().zip(mult) {
                        let to = to as usize;
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
        let mut out = if pooled { f.lists.take() } else { Vec::new() };
        out.extend(
            self.order
                .iter()
                .copied()
                .filter(|&i| self.values[i].abs() > R::c(TOL))
                .map(|i| (i, self.values[i])),
        );
        Some((out, count, work))
    }
}

/// Scratch for independent BTRAN calls. An empty clone avoids copying dead work vectors.
#[derive(Default)]

struct BtranWork<R> {
    x: Vec<R>,
    dense_tail: Vec<R>,
    partial: Vec<R>,
    buckets: Vec<Vec<usize>>,
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

/// What the L stage hands to the U stage: the sparse list `(position, value)`, or the dense work vector itself. The consumer's
/// first pass over `x` turns every `v == 0` slot (+0.0 or -0.0) into +0.0, the bits the old list round trip produced.
enum LOut<R> {
    List(Vec<(usize, R)>),
    Dense(Vec<R>),
}

// Test hook: hand the dense L-stage vector on through the old list round trip (the `dense_spike` list), to compare the two paths.
#[cfg(test)]
thread_local! {
    static FORCE_LIST_HANDOFF: std::cell::Cell<bool> = const { std::cell::Cell::new(false) };
    /// Dense L-stage vectors this test thread handed on as the array itself.
    static DENSE_HANDOFFS: std::cell::Cell<usize> = const { std::cell::Cell::new(0) };
}

/// Working values are indexed by permutation slot. A clear occupancy bit represents +0,
/// so writes from previous right-hand sides need no numeric-array clear.
#[derive(Default)]

struct SparseBtranTail<R> {
    values: Vec<R>,
    occupied: Vec<u64>,
    touched_words: Vec<usize>,
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
            if *word == 0 {
                self.touched_words.push(slot / 64);
            }
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
        for word in self.touched_words.drain(..) {
            self.occupied[word] = 0;
        }
        self.values.resize(m, R::ZERO);
        self.occupied.resize(m.div_ceil(64), 0);
    }
}

/// The set of tail columns without nucleus entries. Membership flags plus an unordered list: the only order-dependent query is the
/// demotion's maximum by slot, whose key `(slot, column)` is unique, and `iter` (dumps) yields ascending columns as the former `BTreeSet`.
#[derive(Clone, Default)]
struct EmptyTail {
    list: Vec<usize>,
    has: Vec<bool>,
}

impl EmptyTail {
    fn insert(&mut self, x: usize) {
        if x >= self.has.len() {
            self.has.resize(x + 1, false);
        }
        if !self.has[x] {
            self.has[x] = true;
            self.list.push(x);
        }
    }

    fn remove(&mut self, x: &usize) {
        if self.has.get(*x).copied().unwrap_or(false) {
            self.has[*x] = false;
            let p = self.list.iter().position(|y| y == x).unwrap();
            self.list.swap_remove(p);
        }
    }

    /// The member with the highest slot (ties, impossible for a slot permutation, broken by the larger column).
    fn max_by_slot(&self, slot: &[usize]) -> Option<usize> {
        self.list.iter().copied().max_by_key(|&x| (slot[x], x))
    }

    #[cfg(any(test, feature = "factor-shadow"))]
    fn iter(&self) -> std::vec::IntoIter<&usize> {
        let mut v: Vec<&usize> = self.list.iter().collect();
        v.sort_unstable();
        v.into_iter()
    }
}

impl FromIterator<usize> for EmptyTail {
    fn from_iter<I: IntoIterator<Item = usize>>(it: I) -> Self {
        let mut s = EmptyTail::default();
        for x in it {
            s.insert(x);
        }
        s
    }
}

#[derive(Clone)]

pub struct DualFactor<R: Real = f64> {
    m: usize,
    btran_ws: btran_ws::Workspace<R>,
    /// The pivot-row BTRAN of this position ran since the last update: the reference then finds the row-eta solve cached (`F+0x1d8`, `0xf193e`).
    row_cached: Option<usize>,
    kind: Vec<i64>,
    nuc: Vec<usize>,
    i21: usize,
    /// The reference's upper bound `[L+0xe8]` on the position of a tail column without entries, scanned downward at updates.
    mark: usize,
    /// The A and B lists of the tail columns with their coefficients (the bits of the matching `uc` entries).
    ab: AbLists<R>,
    /// Immutable forward L columns, preserving coefficient and destination order.
    lpool: LPool<R>,
    /// Every L destination is above its source step (`y > v`), so a dense forward scan sees final values.
    l_forward: bool,
    perm: Vec<usize>,
    /// Immutable inverse row permutation, initialized at the first BTRAN output.
    row_of: std::cell::OnceCell<Vec<usize>>,
    dens: f64,
    dense_est: std::cell::Cell<f64>,
    /// Density threshold fixed at factorization time (`f104` in the notes): 2, or `5 (nnzL/m + 1)^2` capped at `m`.
    dthr: usize,
    /// Columns eligible for demotion, keyed by identity so slot walks need no index repairs.
    empty_tail: EmptyTail,
    slot: Vec<usize>,
    at: Vec<usize>,
    etas: Vec<Eta<R>>,
    bt_occupied: bool,
    bt_work: Option<Box<BtranWork<R>>>,
    /// Immutable L transpose, in ascending source-column order, shared by all BTRAN calls of this factor.
    lt: Csr<R>,
    /// Immutable forward levels and their maximum for level-bucket L solves (`0x12a990`). Used only
    /// for `m >= 250`, a small level count (`< m/10`) and a short Markowitz L part.
    levels: Option<(Vec<usize>, usize)>,
    level_cache: LevelCache<R>,
    /// A reset workspace of a retired factor, adopted by `level_cache` on its first use (kept apart so the cache holds only live state).
    level_spare: LevelCache<R>,
    /// Retained result lists and L-stage buffers (see `ListPool`, `LScratch`); they survive `recycle`.
    lists: ListPool<R>,
    lscr: std::cell::RefCell<LScratch<R>>,
    /// Immutable transpose-graph levels and maximum, initialized when a BTRAN first uses level buckets.
    lt_levels: std::cell::OnceCell<(Vec<usize>, usize)>,
    lev_built: bool,
    diag: Vec<R>,
    /// Whether some diagonal entry is zero or NaN: a zero slot of the tail then divides to NaN instead of vanishing, so the sparse U stage
    /// must visit it (`ftran_with` keeps the closure walk for such a factor, which only a singular basis gives).
    bad_diag: bool,
    /// Numeric U entries of each tail column `(position, value)`.
    uc: Vec<Vec<(usize, R)>>,
    /// The row-wise copy of the numeric U entries of the tail columns, rows in the reference's physical order (one pool).
    urows: URows<R>,
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
    work: f64,
    /// The L-stage part and the U-stage part of `work` of the last `ftran`, and the L-stage estimate of the last `ftran_needed`.
    work_parts: (f64, f64),
    /// The work estimate of the last `ftran_keep`, which restores `work` itself.
    keep_work: f64,
    /// Whether the result of the last `ftran_keep` was dense (the reference's count -1).
    keep_dense: bool,
    tau_l: std::cell::Cell<f64>,
    /// The U-stage estimate of the last `ftran_needed` (`0x123690`).
    tau_u: std::cell::Cell<f64>,
    /// What the tau solve's U stage returns (`0x123690`, `0x1244df`): the number of nonzero tail columns, or the tail rows `m - i21` when it runs dense; the
    /// two-vector FTRAN adds it to the first result's count (`0x12963d`).
    tau_fin: std::cell::Cell<usize>,
    /// What the last `ftran` and `ftran_needed` leave for the joint U stages of the two-vector FTRAN (`joint_u`).
    f2: std::cell::RefCell<F2<R>>,
    /// Fill of the reference's pool for etas and the row-wise copy of U (`F+0x160`), and the capacity of each row of that copy.
    /// Only the parity of `fill` matters: it decides the alignment of the next eta's value array.
    fill: usize,
    row_cap: Vec<usize>,
    /// Parity of the start of each row in the row-wise copy (the value array of a row starts on a 16-byte boundary when it is 0): the rows of a
    /// fresh factor start on multiples of 4, a row moved to the end of the pool starts at the fill of that moment.
    row_par: Vec<u8>,
    /// Alignment of the scratch eta behind each new diagonal (`ctx+0x1d0`): its value array sits `4 m` bytes after the index array in
    /// an arena whose offsets add `24 maxupd + 8 (m + ncols) + 8`, so the parity is `(maxupd + ncols + 1) % 2` whatever `m` is.
    scratch_peel: usize,
    /// Work of the last column replacement (`0xf15c0`), as `(site, estimate, scale)`.
    upd_work: Vec<(u32, f64, f64)>,
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
    bt_dense: bool,
    /// The factor object's update limit when the pricing lowered it (`F+0x188`, the dual path); the default is `max(20, m^0.7)`.
    upd_limit: Option<usize>,
}

/// Sums over the tail columns that kept a nonzero value in a U solve (see `u_work`).
struct Surv<'a> {
    list: &'a [usize],
    n: usize,
    sa: usize,
    sb: usize,
    su: usize,
}

/// The inputs of the joint U stages of the two-vector FTRAN: for the pivot column and for tau, the count of the L-stage list (`None`: dense, count -1)
/// and the solved values (the tail entries are the values after the tail solve); for the column also the positions of the L-stage list.
#[derive(Clone)]

struct F2<R> {
    col_n: Option<usize>,
    col_in: Vec<usize>,
    x_col: Vec<R>,
    /// All 64-row blocks written by the preceding sparse U solve. This list
    /// travels with its snapshot so a clone or replaced factor resets its own data.
    dirty_blocks: Option<Vec<usize>>,
    tau_n: Option<usize>,
    #[cfg(any(test, feature = "factor-shadow"))]
    x_tau: Vec<R>,
    /// Whether the raw value of tau at a tail column exceeded the tolerance before the division (the joint routines test it first).
    /// Only the test and shadow builds keep the snapshot; the work estimates are booked inside the U loops.
    #[cfg(any(test, feature = "factor-shadow"))]
    tau_raw: Vec<bool>,
    /// The tail columns of the column's U solve with a nonzero final value, in the solve's (descending slot) order, and the sums of
    /// their `a` and `b` lengths. Integer-valued work terms are sums over exactly these columns.
    col_surv: Vec<usize>,
    col_sa: usize,
    col_sb: usize,
    /// Tau's tail columns with a nonzero final value whose raw value passed the tolerance (`tau_pass`): the sum of `2 len(a) + 4`, and
    /// the sum of `len(a)` over those whose column value is zero. Only kept when the joint routines run.
    tp_a2: usize,
    tp_only_a: usize,
    /// Retained buffer of the full U solves' output list (empty between calls; only its capacity travels).
    out_scr: Vec<(usize, R)>,
}

impl<R> Default for F2<R> {
    fn default() -> Self {
        F2 {
            col_n: None,
            col_in: Vec::new(),
            x_col: Vec::new(),
            dirty_blocks: None,
            tau_n: None,
            #[cfg(any(test, feature = "factor-shadow"))]
            x_tau: Vec::new(),
            #[cfg(any(test, feature = "factor-shadow"))]
            tau_raw: Vec::new(),
            col_surv: Vec::new(),
            col_sa: 0,
            col_sb: 0,
            tp_a2: 0,
            tp_only_a: 0,
            out_scr: Vec::new(),
        }
    }
}

#[derive(Clone)]

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
#[cfg(test)]
fn row_dot<R: Real>(s: R, t: &[(R, R)], peel: usize) -> R {
    row_dot_by(s, t.len(), peel, |i| t[i].0 * t[i].1)
}

/// `row_dot` over `n` terms given by their products `prod(i)`.
#[inline(always)]
fn row_dot_by<R: Real>(s: R, n: usize, peel: usize, prod: impl Fn(usize) -> R) -> R {
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
#[cfg(test)]
fn row_dot_q<R: Real>(s: R, t: &[(R, R)]) -> R {
    row_dot_q_by(s, t.len(), |i| t[i].0 * t[i].1)
}

/// `row_dot_q` over `n` terms given by their products `prod(i)`.
#[inline(always)]
fn row_dot_q_by<R: Real>(s: R, n: usize, prod: impl Fn(usize) -> R) -> R {
    let nb = n / 8 * 8;
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

/// Temporaries of `DualFactor::fresh` that survive from one factorization to the next. Between two columns `reached` and `changed` are all
/// false and `epos` is all `usize::MAX`.
#[derive(Default)]
struct FreshScratch {
    reached: Vec<bool>,
    changed: Vec<bool>,
    epos: Vec<usize>,
    lst: Vec<usize>,
}

/// The lists of a retired factor (`ab`, `uc`, `urows`), kept with their allocations for the next `fresh`
/// (`DualFactor::recycle`): the refactorization then neither frees the old rows nor allocates the new ones.
struct Recycled<R> {
    ab: AbLists<R>,
    uc: Vec<Vec<(usize, R)>>,
    urows: URows<R>,
    /// Scratch that carries only capacities: the level workspace (reset), the result lists and the L-stage buffers.
    level: Option<Box<LevelWorkspace<R>>>,
    lists: ListPool<R>,
    lscr: LScratch<R>,
    /// The level buckets of the BTRAN workspace (empty rows).
    bt_buckets: Vec<Vec<usize>>,
}

impl<R: Real> Recycled<R> {
    /// The retained lists, or new empty ones, as `m` empty rows each.
    fn take(m: usize) -> Self {
        fn rows<T>(v: &mut Vec<Vec<T>>, m: usize) {
            v.resize_with(m, Vec::new);
            v.iter_mut().for_each(Vec::clear);
        }
        let kept = RECYCLED.with(|r| r.take()).and_then(|b| b.downcast::<Self>().ok());
        let mut r = kept.map_or_else(
            || Self {
                ab: AbLists::new(),
                uc: Vec::new(),
                urows: URows::new(),
                level: None,
                lists: ListPool::default(),
                lscr: LScratch::default(),
                bt_buckets: Vec::new(),
            },
            |b| *b,
        );
        r.ab.reset(m);
        rows(&mut r.uc, m);
        r.urows.reset(m);
        r
    }
}

thread_local! {
    static FRESH_SCRATCH: std::cell::Cell<FreshScratch> = std::cell::Cell::new(FreshScratch::default());
    static RECYCLED: std::cell::Cell<Option<Box<dyn std::any::Any>>> = const { std::cell::Cell::new(None) };
}

impl<R: Real> DualFactor<R> {
    /// Retire the factor: its lists (`ab`, `uc`, `urows`) are kept, allocations included, for the next `fresh` of this
    /// scalar type on this thread (the rest is dropped). Only capacities carry over; `fresh` clears every row.
    pub fn recycle(self) {
        let Self { ab, uc, urows, level_cache, lists, lscr, bt_work, .. } = self;
        let mut bt_buckets = bt_work.map(|w| w.buckets).unwrap_or_default();
        bt_buckets.iter_mut().for_each(Vec::clear);
        let mut level = level_cache.0.into_inner();
        if let Some(l) = level.as_mut() {
            l.reset();
        }
        let lscr = lscr.into_inner();
        RECYCLED.with(|r| r.set(Some(Box::new(Recycled { ab, uc, urows, level, lists, lscr, bt_buckets }))));
    }

    /// Fresh factor from the pivot sequence; `cols[pos]` is the basis column at incoming position `pos` as `(is_slack, rows, values)`.
    ///
    /// Every structure is built to its final size from counting passes (the nested lists with their exact capacities), and the
    /// per-column temporaries live in a scratch that survives from one factorization to the next.
    pub fn fresh<'a>(
        m: usize,
        ncols: usize,
        st: &Structure<R>,
        cols: impl Into<crate::lu_order::Cols<'a>>,
    ) -> DualFactor<R> {
        let cols: crate::lu_order::Cols<'a> = cols.into();
        let mut step = vec![usize::MAX; m];
        for (s, &(r, _)) in st.pivots.iter().enumerate() {
            step[r] = s;
        }
        let (mut n1, mut n12, mut n3) = (0usize, 0usize, 0usize);
        for &p in &st.phase {
            n1 += usize::from(p == 1);
            n12 += usize::from(p == 1 || p == 2);
            n3 += usize::from(p == 3);
        }
        let nseq = st.loff;
        let nsing = st.phase.iter().filter(|&&p| p != 0).count();
        let i21 = n12 + n3;
        let mut kind = vec![-1i64; m];
        let mut nuc = Vec::with_capacity(n3.min(m));
        for s in n12..i21 {
            kind[s] = nuc.len() as i64;
            nuc.push(s);
        }
        for s in i21..m {
            kind[s] = -2;
        }
        let mut lpool = LPool::with_capacity(m, st.l_nnz());
        for s in 0..st.steps() {
            lpool.push_line(st.l_rows(s).iter().zip(st.l_vals(s)).map(|(&r, &v)| (step[r], v)));
        }
        let lcolumn = |v: usize| lpool.iter(v);
        let lt = lpool.transpose(m);
        let l_forward = (0..lpool.lines()).all(|v| lpool.col(v).0.iter().all(|&y| y as usize > v));
        // Pivot-row values of the Markowitz steps by column, `(step, value)` in ascending step order: counting sort.
        let mut mk_start = vec![0usize; m + 1];
        for s in nseq..m {
            for &(c, _) in st.u_row(s) {
                mk_start[c + 1] += 1;
            }
        }
        for c in 0..m {
            mk_start[c + 1] += mk_start[c];
        }
        let mut mk: Vec<(usize, R)> = Vec::new();
        if let Some(&(_, first)) = (nseq..m).find_map(|s| st.u_row(s).first()) {
            mk.resize(mk_start[m], (0, first));
            let mut next: Vec<usize> = mk_start[..m].to_vec();
            for s in nseq..m {
                for &(c, v) in st.u_row(s) {
                    mk[next[c]] = (s, v);
                    next[c] += 1;
                }
            }
        }
        let Recycled { mut ab, mut uc, mut urows, level, lists, lscr, bt_buckets } = Recycled::<R>::take(m);
        let (mut sa, mut sb) = (0usize, 0usize);
        let mut xv = vec![R::ZERO; m];
        let mut ents: Vec<(usize, R)> = Vec::new();
        let (mut av, mut bv): (Vec<(u32, R)>, Vec<(u32, R)>) = (Vec::new(), Vec::new());
        let mut sc = FRESH_SCRATCH.with(|s| s.take());
        // Invariants between columns (restored below): `reached` and `changed` all false, `epos` all MAX.
        sc.reached.resize(m, false);
        sc.changed.resize(m, false);
        sc.epos.resize(m, usize::MAX);
        let FreshScratch { reached, changed, epos, lst } = &mut sc;
        let mut dir: Vec<(usize, R)> = Vec::new();
        for t in i21..m {
            let c = st.pivots[t].1;
            let (_, idx, val) = &cols[c];
            // A tail column without a raw entry in a pair row (`nseq <= s < nsing`) is built directly: the raw copies of the
            // singleton rows, the elimination's own values (`mk`) above the tolerance. The L-solve replay is needed only for
            // the pair rows' values (and for the x87 scalar, where the rule is not established).
            // Justification and corpus check: plans/dual-factor/K0c-result.md section 1. The release guards below fall back to the
            // replay wherever its argument does not hold structurally: a raw entry of a singleton row whose L column reaches a
            // step <= t (it would propagate), a raw entry in a row the replay keeps but the direct rule drops (not kind >= 0 or
            // phase 4), and a repeated row in the column (unmerged MPS input; the replay sums the copies).
            let mut direct = !R::QUAD && !idx.iter().any(|&i| step[i] < t && step[i] >= nseq && step[i] < nsing);
            dir.clear();
            if direct {
                for (&i, &v) in idx.iter().zip(val) {
                    let s = step[i];
                    if s < t && s < nseq && R::c(v) != R::ZERO {
                        if !(kind[s] >= 0 || st.phase[s] == 4) || lcolumn(s).any(|(y, _)| y <= t) {
                            direct = false;
                            break;
                        }
                        dir.push((s, R::c(v)));
                    }
                }
            }
            if direct {
                for &(sidx, vv) in &mk[mk_start[c]..mk_start[c + 1]] {
                    if sidx < t && vv.abs() > R::c(TOL) {
                        dir.push((sidx, vv));
                    }
                }
                dir.sort_unstable_by_key(|e| e.0);
                if dir.windows(2).any(|w| w[0].0 == w[1].0) {
                    direct = false;
                }
            }
            if !direct || cfg!(any(debug_assertions, feature = "factor-shadow")) {
                lst.clear();
                for (&i, &v) in idx.iter().zip(val) {
                    let s = step[i];
                    if s <= t {
                        xv[s] += R::c(v);
                        if !reached[s] {
                            reached[s] = true;
                            lst.push(s);
                        }
                    }
                }
                // Only the reached set is needed: the numeric solve below visits it in ascending order.
                let mut k = 0;
                while k < lst.len() {
                    let v = lst[k];
                    k += 1;
                    for (y, _) in lcolumn(v) {
                        if y <= t && !reached[y] {
                            reached[y] = true;
                            lst.push(y);
                        }
                    }
                }
                // The numeric solve visits the steps in ascending order, not graph-discovery order: sums that cancel to rounding
                // noise then fall on the reference's side of the tolerance (PILOT.JA factor 45, column 925).
                lst.sort_unstable();
                for &v in lst.iter() {
                    if v >= t {
                        continue;
                    }
                    let x = xv[v];
                    // Noise at or below the tolerance does not propagate (a cancellation residue of 1.2e-13 made PEROLD's tail
                    // column 600 carry a spurious entry).
                    if x == R::ZERO || x.abs() <= R::c(TOL) {
                        continue;
                    }
                    for (e, (y, multiplier)) in lcolumn(v).enumerate() {
                        if y > t {
                            continue;
                        }
                        // A singleton-phase column subtracts `(x * a) / d` with the raw entry `a`, not `x * (a / d)`; the
                        // multiplier is the only form known for pair overrides and nucleus columns.
                        let a = st.l_raw(v)[e];
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
                ents.clear();
                for &v in lst.iter() {
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
                for &(sidx, vv) in &mk[mk_start[c]..mk_start[c + 1]] {
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
                for &v in lst.iter() {
                    xv[v] = R::ZERO;
                    reached[v] = false;
                    changed[v] = false;
                }
                ents.sort_unstable_by_key(|e| e.0);
            }
            if direct {
                #[cfg(any(debug_assertions, feature = "factor-shadow"))]
                assert!(
                    dir.len() == ents.len() && dir.iter().zip(&ents).all(|(d, e)| d.0 == e.0 && d.1.to_f64().to_bits() == e.1.to_f64().to_bits()),
                    "direct tail column {t} differs from the replay"
                );
                ents.clear();
                ents.extend_from_slice(&dir);
            }
            let na = ents.iter().filter(|e| kind[e.0] == -2).count();
            let nb = ents.iter().filter(|e| kind[e.0] >= 0).count();
            sa += na;
            sb += nb;
            // The reference lists the fresh tail column's A entries in descending order.
            av.clear();
            av.reserve(na);
            av.extend(ents.iter().rev().filter(|e| kind[e.0] == -2).map(|e| (e.0 as u32, e.1)));
            bv.clear();
            bv.reserve(nb);
            bv.extend(ents.iter().filter(|e| kind[e.0] >= 0).map(|e| (e.0 as u32, e.1)));
            bv.sort_unstable_by_key(|e| std::cmp::Reverse(kind[e.0 as usize]));
            ab.set(t, &av, &bv);
            uc[t].extend_from_slice(&ents);
        }
        FRESH_SCRATCH.with(|s| s.set(sc));
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
        let nnz: usize = lpool.len();
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
        let nnz_tail: usize = lpool.len() - lpool.begin(nseq);
        let levels = {
            let d = m - nseq;
            let outer = d < 6 || 10 * d * d <= m || nnz_tail <= 10 * d;
            if m >= 250 {
                let mut lev = vec![0usize; m];
                for v in 0..m {
                    for (y, _) in lcolumn(v) {
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
        let nnz_l_tot: usize = lpool.len();
        let perm = step;
        let znz = (nnz + sa + sb + 2 * m) as f64;
        let empty_list: Vec<usize> =
            (0..m).filter(|&x| kind[x] == -2 && ab.is_empty(x)).collect();
        // Only tail positions (`kind == -2`) are in the set: its largest element is the last one the downward scan from the end finds.
        let mark0 = empty_list.last().copied().unwrap_or(i21.min(m - 1));
        let empty_tail: EmptyTail = empty_list.into_iter().collect();
        // Entries per row of U: the capacity of each row of the row-wise copy, of `urows` and the layout of the rows.
        let mut row_cap = vec![0usize; m];
        for col in &uc {
            for &(row, _) in col {
                row_cap[row] += 1;
            }
        }
        let u_total: usize = row_cap.iter().sum();
        // `update` appends to the rows of `urows`: nonempty rows get room for a few entries (capacity only, never visible).
        let spare = |n: usize| if n == 0 { 0 } else { n + (n >> 2) + 2 };
        for (row, &n) in row_cap.iter().enumerate() {
            if n > 0 {
                urows.reserve_row(row, spare(n));
            }
        }
        // Filled in ascending column order; the rows are put in the reference's physical order below.
        for (c, col) in uc.iter().enumerate() {
            for &(row, v) in col {
                urows.push(row, c, v);
            }
        }
        // A row lists its tail columns in the order of the factorization's U row (Markowitz steps), or in the order of the
        // basis columns for a row of the singleton phases.
        {
            // member lists (the tail columns whose A or B list holds the row, ascending) as one counting-sorted array
            let mut mem_start = vec![0usize; m + 1];
            for t in i21..m {
                for &y in ab.both(t) {
                    mem_start[y as usize + 1] += 1;
                }
            }
            for y in 0..m {
                mem_start[y + 1] += mem_start[y];
            }
            let mut member = vec![0usize; mem_start[m]];
            let mut next: Vec<usize> = mem_start[..m].to_vec();
            for t in i21..m {
                for &y in ab.both(t) {
                    member[next[y as usize]] = t;
                    next[y as usize] += 1;
                }
            }
            let mut pos_of_col = vec![usize::MAX; m];
            for (t, &(_, c)) in st.pivots.iter().enumerate() {
                pos_of_col[c] = t;
            }
            // Marker words of the rows with a pair order: the hash is read only for those.
            let mut has_urow = vec![0u64; m.div_ceil(64)];
            for &y in st.pair_urow.keys() {
                if y < m {
                    has_urow[y / 64] |= 1 << (y % 64);
                }
            }
            let mut stamp = vec![usize::MAX; m];
            let mut tmp: Vec<usize> = Vec::new();
            for y in 0..m {
                let mem = &member[mem_start[y]..mem_start[y + 1]];
                if mem.is_empty() {
                    continue;
                }
                tmp.clear();
                if let Some(u) = if has_urow[y / 64] >> (y % 64) & 1 == 0 { None } else { st.pair_urow.get(&y) } {
                    for &t in mem {
                        stamp[t] = y;
                    }
                    tmp.extend(u.iter().map(|&c| pos_of_col[c]).filter(|&t| stamp[t] == y));
                } else if st.u_row(y).is_empty() {
                    // Filtering all basis columns in index order is equivalent to sorting just this row's members by
                    // column. Keep only the positions represented by pos_of_col, and emit each position once.
                    tmp.extend(mem.iter().copied().filter(|&t| pos_of_col[st.pivots[t].1] == t));
                    tmp.sort_unstable_by_key(|&t| st.pivots[t].1);
                    tmp.dedup();
                } else {
                    for &t in mem {
                        stamp[t] = y;
                    }
                    tmp.extend(st.u_row(y).iter().map(|&(c, _)| pos_of_col[c]).filter(|&t| stamp[t] == y));
                }
                // The row holds the same columns (audited): `xv` carries the values over to the new order.
                let (cs, vs) = urows.row_mut(y);
                assert_eq!(cs.len(), tmp.len(), "fresh: rowl members differ from the row");
                for (&c, &v) in cs.iter().zip(vs.iter()) {
                    xv[c as usize] = v;
                }
                for ((c, v), &t) in cs.iter_mut().zip(vs.iter_mut()).zip(&tmp) {
                    *c = t as u32;
                    *v = xv[t];
                }
            }
        }
        // A fresh factor lays the rows of the row-wise copy out contiguously: the nonempty rows of the nucleus in the order of their list index
        // (`kind`), then the others by position, each row as long as its U entries (`ucol_len`); the parity of the start decides the alignment of the row's
        // value array in the dot product of the second FTRAN vector.
        let row_par = {
            let mut par = vec![0u8; m];
            let mut acc = 0usize;
            // `nuc` lists the steps with `kind >= 0` in `kind` order.
            for p in nuc.iter().copied().chain((0..m).filter(|&p| kind[p] < 0)) {
                if row_cap[p] > 0 {
                    par[p] = (acc % 2) as u8;
                    acc += row_cap[p];
                }
            }
            par
        };
        let fresh = DualFactor {
            upd: 0,
            sum_alpha: 0.0,
            before_last: (0, 0.0),
            sum_eta: 0.0,
            last_alpha: 0,
            znz,
            n1,
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
            scratch_peel: ((m as f64).powf(0.7) as usize)
                .max(20)
                .wrapping_add(ncols + 1)
                % 2,
            fill: 2 * (u_total + nnz_l_tot),
            urows,
            row_cap,
            row_par,
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
            empty_tail,
            kind,
            nuc,
            i21,
            mark: mark0,
            ab,
            lpool,
            l_forward,
            perm,
            row_of: Default::default(),
            dens: (h33 + h34) * (h33 + h34),
            dense_est: std::cell::Cell::new(0.0),
            dthr,
            slot: (0..m).collect(),
            at: (0..m).collect(),
            btran_ws: Default::default(),
            etas: Vec::new(),
            lt,
            levels,
            level_cache: LevelCache::default(),
            level_spare: LevelCache(std::cell::RefCell::new(level)),
            lists,
            lscr: std::cell::RefCell::new(lscr),
            lt_levels: Default::default(),
            // Four or more occupancy words, and skipped zero divisions must always be exact.
            bt_occupied: !R::QUAD && m >= 256 && st.diag.iter().all(|d| d.to_f64().is_normal()),
            bt_work: Some(Box::new(BtranWork { buckets: bt_buckets, ..Default::default() })),
            lev_built,
            diag: st.diag.clone(),
            bad_diag: st.diag.iter().any(|d| !(d.abs() > R::ZERO)),
            uc,
            spike: Vec::new(),
            update_values: Vec::new(),
            update_seen: Vec::new(),
            spike_used: false,
            lsnap: Default::default(),
        };
        // Every factor any test builds is also compared with the construction it replaced.
        #[cfg(test)]
        fresh_tests::assert_same(&fresh, &Self::fresh_old(m, ncols, st, cols));
        #[cfg(any(test, debug_assertions, feature = "factor-shadow"))]
        {
            fresh.ab.check(&fresh.uc, &fresh.kind);
            fresh.urows.check(&fresh.uc, &fresh.kind);
        }
        fresh
    }

    /// Destinations and multipliers of L column `v`.
    #[inline(always)]
    fn lcolumn(&self, v: usize) -> (&[u32], &[R]) {
        self.lpool.col(v)
    }

    /// L columns as `(destination, multiplier)` lists (projection and dump views).
    pub(crate) fn lcolumns(&self) -> Vec<Vec<(usize, R)>> {
        (0..self.m).map(|v| self.lpool.iter(v).collect()).collect()
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

    /// Counts the entries with `v != 0`. The dense L-stage vector holds no -0.0: it starts as +0.0, the L stage and the first-touch
    /// scatter only add or subtract (a sum is -0.0 only when both addends are), and an eta result at or below the tolerance is stored as +0.0.
    /// The old pass stored +0.0 into every zero slot (a -0.0 would have read as +0.0 through the list round trip); the claim is checked
    /// in the test, debug and shadow builds.
    #[inline(never)]
    fn normalize_count(x: &mut [R]) -> usize {
        #[cfg(any(test, debug_assertions, feature = "factor-shadow"))]
        assert!(x.iter().all(|&v| v != R::ZERO || ab::same_bits(v, R::ZERO)), "dense L-stage vector holds a -0.0");
        x.iter().map(|&v| usize::from(v != R::ZERO)).sum()
    }

    /// One pass over a dense L-stage vector (no -0.0, see `normalize_count`): `inp` receives the positions of
    /// the entries with `v != 0` and `sp` the `(position, value)` entries above the tolerance, both in ascending position, without a
    /// conditional push per entry. Returns the number of entries with `v != 0`.
    #[inline(never)]
    fn dense_lists(x: &[R], inp: &mut Vec<usize>, sp: &mut Vec<(usize, R)>) -> usize {
        #[cfg(any(test, debug_assertions, feature = "factor-shadow"))]
        assert!(x.iter().all(|&v| v != R::ZERO || ab::same_bits(v, R::ZERO)), "dense L-stage vector holds a -0.0");
        inp.clear();
        sp.clear();
        inp.reserve(x.len() + 1);
        let pi = inp.as_mut_ptr();
        let mut ci = 0usize;
        for (i, &v) in x.iter().enumerate() {
            // `ci <= i < x.len()` before the write and `inp` holds `x.len() + 1` entries; a rejected position is overwritten by the next
            // one. The type is Copy, so entries beyond the length need no drop.
            unsafe { pi.add(ci).write(i) };
            ci += usize::from(v != R::ZERO);
        }
        unsafe { inp.set_len(ci) };
        // The spike is the part of the list above the tolerance, with the values (`cs < ci + 1 <= capacity` at every write).
        sp.reserve(ci + 1);
        let ps = sp.as_mut_ptr();
        let mut cs = 0usize;
        for &i in inp.iter() {
            let v = x[i];
            unsafe { ps.add(cs).write((i, v)) };
            cs += usize::from(v.abs() > R::c(TOL));
        }
        unsafe { sp.set_len(cs) };
        ci
    }

    /// Pack a dense BTRAN result in row order without conditional vector growth.
    #[inline(never)]
    fn dense_btran_output(x: &[R], perm: &[usize], mut out: Vec<(usize, R)>) -> Vec<(usize, R)> {
        out.clear();
        out.reserve(perm.len());
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
            for &x in self.lcolumn(v).0 {
                let x = x as usize;
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
            for &x in self.lcolumn(v).0 {
                let x = x as usize;
                cnt[x] -= 1;
                if cnt[x] == 0 {
                    q.push(x);
                }
            }
        }
        q
    }

    /// `count_l * dthr > m` for a count the caller reads: `lstage_num` returns `None` only when `rows_dense`, `lev_dense_ft` or (for
    /// the `force_dense` readers) the dense request already decides the test, or there are no etas, and the readers test those first.
    #[inline(always)]
    fn count_dense(count_l: Option<usize>, dthr: usize, m: usize) -> bool {
        count_l.expect("count_l is read only where lstage_num counted it") * dthr > m
    }

    /// Test oracle for the dense L sweep: the full sweep over every step on a copy of the input (`before`) against the restricted
    /// one (`after`), bit for bit, with the count where it was asked for and the length sum.
    #[cfg(test)]
    fn assert_full_sweep(&self, before: &[R], after: &[R], count: Option<usize>, len_sum: usize) {
        let mut full = before.to_vec();
        let (c, s) = pool::plain::fused_forward(&self.lpool, &mut full);
        assert!(full.iter().zip(after).all(|(&a, &b)| pool::plain::same_bits(a, b)), "x");
        assert_eq!(s, len_sum);
        if let Some(n) = count {
            assert_eq!(n, c);
        }
        if count.is_some() && self.lpool.gap_is_empty(self.n1, self.loff) && self.n1 < self.loff {
            LIVE_COUNTS.with(|k| k.set(k.get() + 1));
        }
    }

    /// L solve (Kahn order) and row etas on the entering column; returns the nonzero list the U solve receives, with values.
    /// Also returns the work estimate of the L stage and the row etas (`0x12bd30` / `0x12a990` and the eta loop of `0x120740`).
    fn lstage_num(&self, rows: &[usize], vals: &[R], force_dense: bool) -> (LOut<R>, Option<usize>, f64) {
        self.lstage_list(rows, vals, force_dense, false)
    }

    /// `lstage_num`; with `pooled` a list result is taken from the retained pool, which only a caller that gives it back
    /// (`ftran_needed`) may ask for.
    fn lstage_list(&self, rows: &[usize], vals: &[R], force_dense: bool, pooled: bool) -> (LOut<R>, Option<usize>, f64) {
        let m = self.m;
        let mut e = 0.0f64;
        let dense_l = rows.len() * self.dthr > self.m || force_dense || self.lev_dense_ft;
        if R::LANES && !dense_l {
            if let Some((levels, maxlev)) = &self.levels {
                let mut cache = self.level_cache.0.borrow_mut();
                if let Some(result) = cache
                    .get_or_insert_with(|| self.level_spare.0.take().unwrap_or_default())
                    .solve(self, levels, *maxlev, rows, vals, pooled)
                {
                    return (LOut::List(result.0), Some(result.1), result.2);
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
        // `count_l` is `None` where no reader can see it: the readers (here and in `ftran_needed` / `ftran_with`) test it only with
        // etas present and after `rows_dense`, `lev_dense_ft` (and, in `ftran_needed`, `force_dense`) have been ruled out.
        let (mut order, count_l) = if dense_l {
            // Dense L solve (`0x12bd30`): every node in ascending order; the list handed on is the ascending scan.
            let need_count = force_dense && rows.len() * self.dthr <= self.m && !self.lev_dense_ft && !self.etas.is_empty();
            e = m as f64 + rows.len() as f64 * 2.0 + e;
            e = (m - self.loff) as f64 * 2.0 + self.n1 as f64 + self.n1 as f64 + e;
            let count_l = if self.l_forward {
                // Every destination is above its source, so `x[v]` is final when visited and `count_l` is counted in the same pass.
                // `e` is integer-valued here (sums of `m`, counts and doubled counts), so the doubled column lengths are added at once.
                debug_assert!(e.fract() == 0.0);
                // The singleton steps `[n1, loff)` have no L entries and are not visited; their nonzeros are scanned only if read.
                #[cfg(test)]
                let before = x.clone();
                let (count_l, len_sum) = pool::fused_forward(&self.lpool, &mut x, self.n1, self.loff, need_count);
                #[cfg(test)]
                self.assert_full_sweep(&before, &x, count_l, len_sum);
                e = len_sum as f64 * 2.0 + e;
                count_l
            } else {
                let (n1, loff) = if self.lpool.gap_is_empty(self.n1, self.loff) { (self.n1, self.loff) } else { (self.m, self.m) };
                for v in (0..n1).chain(loff..self.m) {
                    let xv = x[v];
                    if xv.abs() <= R::c(TOL) {
                        continue;
                    }
                    let (dest, mult) = self.lcolumn(v);
                    e = dest.len() as f64 * 2.0 + e;
                    for (&y, &a) in dest.iter().zip(mult) {
                        x[y as usize] -= a * xv;
                    }
                }
                need_count.then(|| (0..self.m).filter(|&v| x[v].abs() > R::c(TOL)).count())
            };
            (Vec::new(), count_l)
        } else if let Some((lev, maxlev)) = &self.levels {
            // Level-bucket L solve: nodes are taken level by level, each level in insertion order (inputs first, then the
            // nodes reached from lower levels); values below the tolerance are dropped without propagating.
            let maxlev = *maxlev;
            let (mut buckets, mut active, mut out) = {
                let mut scr = self.lscr.borrow_mut();
                (std::mem::take(&mut scr.buckets), std::mem::take(&mut scr.active), std::mem::take(&mut scr.order))
            };
            // Retained buffers: `active` all false and the buckets empty between calls (restored below), `out` cleared here.
            #[cfg(test)]
            assert!(active.iter().all(|&a| !a) && buckets.iter().all(Vec::is_empty), "stale L-stage scratch");
            buckets.resize_with(maxlev + 1, Vec::new);
            active.resize(self.m, false);
            out.clear();
            let mut l = usize::MAX;
            for &v in &start {
                if !active[v] {
                    active[v] = true;
                    buckets[lev[v]].push(v);
                    l = l.min(lev[v]);
                }
            }
            let mut remaining: usize = buckets.iter().map(Vec::len).sum();
            while remaining > 0 && l <= maxlev {
                let mut i = 0;
                while i < buckets[l].len() {
                    let v = buckets[l][i];
                    i += 1;
                    remaining -= 1;
                    let xv = x[v];
                    if xv.abs() >= R::c(TOL) {
                        out.push(v);
                        let (dest, mult) = self.lcolumn(v);
                        for (&y, &a) in dest.iter().zip(mult) {
                            let y = y as usize;
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
            // Every queued node sits in its bucket: clear the markers through them, then the buckets.
            for b in buckets.iter_mut() {
                for &v in b.iter() {
                    active[v] = false;
                }
                b.clear();
            }
            {
                let mut scr = self.lscr.borrow_mut();
                scr.buckets = buckets;
                scr.active = active;
            }
            let n = out.len();
            e = ((self.nnz_l + m + 2 * m) as f64 * 5.0 * n as f64) / (m + 1) as f64 + e;
            (out, Some(n))
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
                let (dest, mult) = self.lcolumn(v);
                for (&y, &a) in dest.iter().zip(mult) {
                    x[y as usize] -= a * xv;
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
                    for &y in self.lcolumn(v).0 {
                        let y = y as usize;
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
                        e2 += self.lpool.col_len(v) as f64 * 2.0;
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
                    let (dest, mult) = self.lcolumn(v);
                    for (&y, &a) in dest.iter().zip(mult) {
                        x[y as usize] -= a * xv;
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
                    Some(count_l),
                )
            } else {
                (q, Some(count_l))
            }
        };
        if crate::dump::on("lstage") {
            self.lsnap.borrow_mut().0 = x.clone();
        }
        // The sparse L stage drops a node whose value is at or below the tolerance from its list; an eta that makes it nonzero again appends it
        // at the end (SCAGR7 iteration 173: `0x129aa0` start `[82, 125, 83]` leaves `[82, 49, 125]`).
        let output_dense = dense_l || (!self.etas.is_empty() && Self::count_dense(count_l, self.dthr, m));
        let mut have = Vec::new();
        if !output_dense {
            order.retain(|&v| x[v].abs() > R::c(TOL));
            have.resize(m, false);
            for &v in &order {
                have[v] = true;
            }
        }
        // Gather each eta into one buffer; its logical peel is stored in eta.peel.
        let mut eta_x: Vec<R> = std::mem::take(&mut self.lscr.borrow_mut().eta_x);
        if !self.etas.is_empty() {
            let netas = self.etas.len();
            if dense_l || Self::count_dense(count_l, self.dthr, m) {
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
        self.lscr.borrow_mut().eta_x = eta_x;
        if crate::dump::on("lstage") {
            self.lsnap.borrow_mut().1 = x.clone();
        }
        // A dense vector goes on as the raw array: the reference keeps the sub-tolerance residues of the L and eta stages.
        if output_dense {
            order.clear();
            self.lscr.borrow_mut().order = order;
            if !R::LANES {
                return (LOut::List((0..m).filter(|&v| x[v] != R::ZERO).map(|v| (v, x[v])).collect()), count_l, e);
            }
            // The U stage reads the array itself (the reference keeps a dense result dense).
            #[cfg(test)]
            if FORCE_LIST_HANDOFF.with(std::cell::Cell::get) {
                return (LOut::List(Self::dense_spike(&x)), count_l, e);
            }
            #[cfg(test)]
            DENSE_HANDOFFS.with(|k| k.set(k.get() + 1));
            return (LOut::Dense(x), count_l, e);
        }
        let mut list = if pooled { self.lists.take() } else { Vec::new() };
        list.extend(order.iter().copied().filter(|&v| x[v].abs() > R::c(TOL)).map(|v| (v, x[v])));
        order.clear();
        self.lscr.borrow_mut().order = order;
        (LOut::List(list), count_l, e)
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
    /// value, a needed nucleus position the row-wise dot product over the row pool (`urows`), subtracted term by term from the L
    /// stage value and divided by the diagonal, without a tolerance test. Does not touch the spike of the last `ftran`.
    pub fn ftran_needed(
        &self,
        rows: &[usize],
        vals: &[R],
        needed: &[usize],
        dense_in: bool,
    ) -> Vec<(usize, R)> {
        let (m, i21) = (self.m, self.i21);
        let (lout, count_l, e_l) = self.lstage_list(rows, vals, dense_in, true);
        // A dense L-stage vector is the U stage's work vector itself (it implies `tau_dense`); only its counts are needed.
        let (spike, dense_x) = match lout {
            LOut::List(list) => (list, None),
            LOut::Dense(mut x) => {
                let nz = Self::normalize_count(&mut x);
                (Vec::new(), Some((x, nz)))
            }
        };
        self.tau_l.set(e_l);
        let tau_dense = rows.len() * self.dthr > self.m
            || dense_in
            || self.lev_dense_ft
            || (!self.etas.is_empty() && Self::count_dense(count_l, self.dthr, self.m));
        debug_assert!(dense_x.is_none() || tau_dense);
        let tau_n = if tau_dense { None } else { Some(spike.len()) };
        // The inline dense joint loop of two dense vectors (`0x128c56..0x129152`) tests the raw value before the quotient: a tail column whose
        // raw value is at most the tolerance is zeroed even when the quotient exceeds it.
        // The same holds for the joint routine `0x12c300` (column sparse, tau dense; PILOT.JA call 249): whenever tau does not take the single-vector U stage
        // (`ok(n2)` of `joint_u` fails), the tail columns test the raw value as well.
        let joint_dense = !tau_n.is_some_and(|n| {
            self.nuc.len() + n < self.i21 || self.dens * (n as f64) < 2.0 * m as f64
        });
        // The spike's entry count and its tail-entry count (read only where `m > h33^2 n`), before the tail loop changes `x`.
        let (spike_n, spike_t) = match &dense_x {
            Some((x, nz)) => (*nz, if self.tail_count_needed(*nz) { self.tail_nonzero(x) } else { 0 }),
            None => (
                spike.len(),
                if self.tail_count_needed(spike.len()) { spike.iter().filter(|e| self.kind[e.0] == -2).count() } else { 0 },
            ),
        };
        let tdense = self.tau_u_dense_nt(spike_n, spike_t);
        let keep_tiny = !tau_dense && !joint_dense && tdense;
        // The dense slot loop (`12403b..12405f`, and the joint loops) divides and stores every tail slot, descending, before testing the
        // quotient, including unreachable zero slots: their signed zeros are consumed by requested rows and F2. Otherwise the reached
        // columns are the closure of the spike's tail positions under the `a` lists, in the same descending slot order.
        let dense_order = joint_dense || (!R::QUAD && keep_tiny);
        let mut tail: Vec<usize> = Vec::new();
        let mut x = match dense_x {
            Some((x, _)) => x,
            None => {
                let mut x = vec![R::ZERO; m];
                for &(p, v) in &spike {
                    x[p] = v;
                    if !dense_order && self.kind[p] == -2 {
                        tail.push(p);
                    }
                }
                x
            }
        };
        if !dense_order {
            let mut touched = vec![false; m];
            for &(p, _) in &spike {
                touched[p] = true;
            }
            let mut k = 0;
            while k < tail.len() {
                let c = tail[k];
                k += 1;
                for &j in self.ab.a(c) {
                    let j = j as usize;
                    if !touched[j] {
                        touched[j] = true;
                        tail.push(j);
                    }
                }
            }
            self.sort_tail_desc(&mut tail);
        }
        #[cfg(any(test, feature = "factor-shadow"))]
        let mut tau_raw = vec![false; m];
        #[cfg(not(any(test, feature = "factor-shadow")))]
        let mut tau_raw: [bool; 0] = [];
        // The U stage `0x123690` solves the tail in one of three loops: the heap walk and the sparse-row loop zero a quotient of at most the tolerance, the dense
        // slot loop (`0x12403b`: `m <= h33^2 c` or `5 h33^2 t > m - i21`, see `tau_u_work`) stores the quotient and only skips its scatter (FIT2P call 14
        // keeps a tail value of 4e-16 that enters the row dots of the nucleus positions). A dense input (`tau_dense`) zeroes it (FIT2P call 49), and so does the
        // joint dense loop (`joint_dense`, fx_e2_r25 call 6: a tail quotient of -1.8e-14 is 0 in the reference).
        //
        // The work estimates of the U stage are booked in the loop, where a tail column is finished (a scatter only hits lower slots, so its value
        // is final): every term is an integer-valued count (`len`, `2 len + 4`), hence the sums are exact and order free.
        let counts = {
            let f2 = self.f2.borrow();
            let x_col = &f2.x_col[..];
            let dense_slots = || self.at[i21..m].iter().rev().copied();
            macro_rules! pass {
                ($jd:literal, $kt:literal, $td:literal, $order:expr) => {
                    self.needed_tail::<$jd, $kt, $td>(&mut x, $order, x_col, &mut tau_raw)
                };
            }
            match (joint_dense, keep_tiny, tdense, dense_order) {
                (true, _, false, _) => pass!(true, false, false, dense_slots()),
                (true, _, true, _) => pass!(true, false, true, dense_slots()),
                (false, true, false, true) => pass!(false, true, false, dense_slots()),
                (false, true, true, true) => pass!(false, true, true, dense_slots()),
                (false, true, false, false) => pass!(false, true, false, tail.iter().copied()),
                (false, true, true, false) => pass!(false, true, true, tail.iter().copied()),
                (false, false, false, _) => pass!(false, false, false, tail.iter().copied()),
                (false, false, true, _) => pass!(false, false, true, tail.iter().copied()),
            }
        };
        // Invariant (P4-2): every nonzero tail slot was reached through the a-lists and processed, so the nonzero tail values left are exactly the
        // survivors counted in the loop.
        #[cfg(any(test, debug_assertions, feature = "factor-shadow"))]
        assert_eq!(
            (self.i21..m).filter(|&s| x[self.at[s]] != R::ZERO).count(),
            counts.surv_nz,
            "ftran_needed: nonzero tail slots differ from the processed survivors"
        );
        let (u_est, fin) = self.tau_u_finish(spike_n, spike_t, needed.len(), counts.nz, counts.abs_a);
        self.tau_u.set(u_est);
        self.tau_fin.set(fin);
        #[cfg(test)]
        {
            let (v, f) = self.tau_u_work_ref(spike_n, spike_t, &x, needed.len());
            assert_eq!(v.to_bits(), u_est.to_bits(), "tau_u");
            assert_eq!(f, fin, "tau_fin");
        }
        {
            let mut f = self.f2.borrow_mut();
            f.tau_n = tau_n;
            #[cfg(any(test, feature = "factor-shadow"))]
            {
                f.x_tau = x.clone();
                f.tau_raw = tau_raw;
            }
            f.tp_a2 = counts.tp_a2;
            f.tp_only_a = counts.tp_only_a;
        }
        let mut result = self.lists.take();
        result.extend(needed.iter().map(|&p| {
            let v = match self.kind[p] {
                -2 => x[p],
                -1 => x[p] / self.diag[p],
                _ => self.nucleus_dot(p, &x) / self.diag[p],
            };
            (p, v)
        }));
        self.lists.give(spike);
        result
    }

    /// Hands a result list of `btran` or `ftran_needed` back for reuse once its contents are no longer read.
    pub fn give_list(&self, list: Vec<(usize, R)>) {
        self.lists.give(list);
    }

    /// Number of tail slots whose position holds a nonzero value in the dense vector `x` (the tail positions are `at[i21..m]`).
    fn tail_nonzero(&self, x: &[R]) -> usize {
        debug_assert!((self.i21..self.m).all(|s| self.kind[self.at[s]] == -2));
        self.at[self.i21..self.m].iter().filter(|&&p| x[p] != R::ZERO).count()
    }

    /// The tail loop of `ftran_needed` over the columns `order` (descending slots): divide, store, test, scatter over the column's A list,
    /// with the work-estimate counts of the U stage booked as each column is finished. `JD`: the joint dense loops (raw test first, a
    /// zeroed column is neither divided nor scattered); `KT`: a quotient of at most the tolerance is kept (and not scattered); `TD`: the
    /// dense booking (`abs_a`) instead of the heap-walk booking (`nz`).
    #[inline(always)]
    fn needed_tail<const JD: bool, const KT: bool, const TD: bool>(
        &self,
        x: &mut [R],
        order: impl Iterator<Item = usize>,
        x_col: &[R],
        tau_raw: &mut [bool],
    ) -> TailCounts {
        let mut t = TailCounts::default();
        let diag = &self.diag[..];
        let ab = &self.ab;
        let tol = R::c(TOL);
        for c in order {
            let xr = x[c];
            let raw = xr.abs() > tol;
            #[cfg(any(test, feature = "factor-shadow"))]
            {
                tau_raw[c] = raw;
            }
            #[cfg(not(any(test, feature = "factor-shadow")))]
            let _ = &tau_raw;
            if JD && !raw {
                x[c] = R::ZERO;
                continue;
            }
            let xc = xr / diag[c];
            let tiny = xc.abs() <= tol;
            let fin = if tiny && !KT { R::ZERO } else { xc };
            x[c] = fin;
            if tiny {
                // A kept tiny quotient is a survivor of the dense slot loop (never scattered); a dropped one books nothing.
                if KT && fin != R::ZERO {
                    t.book::<JD, TD, R>(ab.na(c), fin, raw, x_col.get(c));
                }
                continue;
            }
            let (ai, av) = ab.a_pair(c);
            t.book::<JD, TD, R>(ai.len(), fin, raw, x_col.get(c));
            for (&j, &v) in ai.iter().zip(av) {
                x[j as usize] -= v * xc;
            }
        }
        t
    }

    /// The row-wise dot product of nucleus position `p` (`0x1241ed`): `x[p]` minus the products of the row's values with `x` at its columns,
    /// in the row's physical order.
    #[inline(always)]
    fn nucleus_dot(&self, p: usize, x: &[R]) -> R {
        let (cols, vals) = self.urows.row(p);
        let vals = &vals[..cols.len()];
        let prod = |i: usize| vals[i] * x[cols[i] as usize];
        if R::QUAD {
            row_dot_q_by(x[p], cols.len(), prod)
        } else {
            row_dot_by(x[p], cols.len(), self.row_par[p] as usize, prod)
        }
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
        let v = if ok(f.col_n) {
            self.third_work(&f, out)
        } else {
            self.dense_pair_work(&f, out)
        };
        #[cfg(test)]
        {
            let want = if ok(f.col_n) {
                self.third_work_ref(&f, out)
            } else {
                self.dense_pair_work_ref(&f, out)
            };
            assert_eq!(v.to_bits(), want.to_bits(), "joint_u");
        }
        Some(v)
    }

    /// Work of the output stage of both joint routines: `(row length + 2) 2` for each result position of the nucleus list.
    fn row_out_work(&self, out: &[usize]) -> f64 {
        let mut e = 0.0f64;
        for &p in out {
            if self.kind[p] >= 0 {
                e += (self.urows.len(p) as f64 + 2.0) * 2.0;
            }
        }
        e
    }

    /// The joint routine `0x12c300` (column: sparse list, tau: dense). The terms that the reference accumulates while it solves the tail
    /// columns are integer-valued counts, summed in the two U loops (`col_*` by `ftran`, `tp_*` by `ftran_needed`); the non-integral
    /// factors (`2.4`) multiply or follow the finished counts exactly as in `third_work_ref`.
    fn third_work(&self, f: &F2<R>, out: &[usize]) -> f64 {
        let (m, i21) = (self.m, self.i21);
        let mf = m as f64;
        let tail_rows = (m - i21) as f64;
        let cin = f.col_in.len();
        if cin == 0 {
            return 0.0;
        }
        let x5 = f.tp_a2 as f64;
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
        let x1;
        if use_dense_scan {
            x2 += 5.0 * tail_rows;
            if sparse_in {
                x2 += 2.0 * (tail_cnt + nuc_cnt) as f64;
            }
            let x3 = (tail_rows * 3.0 + x2) + (2 * (f.col_sa + f.col_sb)) as f64;
            let x0 = 2.0 * nuc_len as f64 + x3;
            x1 = 2.0 * out.len() as f64 + x0;
        } else {
            let r9 = nuc_cnt + f.col_sa + f.tp_only_a;
            let r15 = nuc_cnt + f.col_sb;
            let surv = &f.col_surv;
            let mut x0 = tail_rows * 5.0;
            x0 += x2;
            x0 += r9 as f64 * 3.0;
            let surv_b = (2 * f.col_sb + 4 * surv.len()) as f64;
            if r15 as f64 * 2.5 > nuc_len as f64 {
                x0 += nuc_cnt as f64;
                x0 += surv_b;
                x0 += nuc_len as f64 * 2.4;
            } else {
                let mut seen = vec![false; m];
                for &p in &nuc_in {
                    seen[p] = true;
                }
                let mut touched = nuc_cnt;
                x0 += surv_b;
                for &c in surv {
                    for &j in self.ab.b(c) {
                        let j = j as usize;
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

    /// The inline dense routine of `0x125f50` (`0x128742`, both vectors dense), from the same integer sums as `third_work`.
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
        let x5 = 5.0 * (m - i21) as f64 + (2 * (f.col_sa + f.col_sb) + 4 * f.col_surv.len()) as f64;
        let mut x4 = f.tp_a2 as f64;
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

    /// Work of the tau solve's U stage from the counts of the tail loop (see `tau_u_work_ref`): `abs_a` is the sum of `len(a)` over the tail
    /// columns whose solved value exceeds the tolerance (dense loop), `nz` the number of nonzero tail columns (heap walk). Returns the
    /// estimate and the finished count.
    fn tau_u_finish(&self, spike_n: usize, t: usize, n: usize, nz: usize, abs_a: usize) -> (f64, usize) {
        let m = self.m as f64;
        let (cnt, nf) = (spike_n as f64, n as f64);
        let d = self.h33;
        let tail_rows = (self.m - self.i21) as f64;
        let mut e = ((10.0 + nf) + nf) + 6.0 * cnt;
        let mut dense = m <= (d * d) * cnt;
        if !dense {
            let t = t as f64;
            if ((5.0 * d) * d) * t > tail_rows {
                dense = true;
                e = 2.0 * t + e;
            }
        }
        debug_assert_eq!(dense, self.tau_u_dense_nt(spike_n, t));
        let fin;
        if dense {
            e = 3.0 * tail_rows + e;
            e += (abs_a * 2) as f64;
            fin = tail_rows;
        } else {
            let nz = nz as f64;
            e = (d * 4.0) * nz + e;
            fin = nz;
        }
        ((((nf + nf) * self.h34 + e) + 4.0 * fin), fin as usize)
    }

    /// The joint routine `0x12c300` (reference copy).
    #[cfg(test)]
    fn third_work_ref(&self, f: &F2<R>, out: &[usize]) -> f64 {
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
                    x5 += (self.ab.na(c) as f64 + 2.0) * 2.0;
                }
                if col_pass(c) {
                    x3 += 2.0 * (self.ab.na(c) + self.ab.nb(c)) as f64;
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
                    r9 += self.ab.na(c);
                }
                if tp {
                    x5 += (self.ab.na(c) as f64 + 2.0) * 2.0;
                }
                if cp {
                    r15 += self.ab.nb(c);
                    surv.push(c);
                }
            }
            let mut x0 = tail_rows * 5.0;
            x0 += x2;
            x0 += r9 as f64 * 3.0;
            if r15 as f64 * 2.5 > nuc_len as f64 {
                x0 += nuc_cnt as f64;
                for &c in &surv {
                    x0 += self.ab.nb(c) as f64 * 2.0 + 4.0;
                }
                x0 += nuc_len as f64 * 2.4;
            } else {
                let mut seen = vec![false; m];
                for &p in &nuc_in {
                    seen[p] = true;
                }
                let mut touched = nuc_cnt;
                for &c in &surv {
                    x0 += self.ab.nb(c) as f64 * 2.0 + 4.0;
                    for &j in self.ab.b(c) {
                        let j = j as usize;
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
    #[cfg(test)]
    fn dense_pair_work_ref(&self, f: &F2<R>, out: &[usize]) -> f64 {
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
                x5 += ((self.ab.na(c) + self.ab.nb(c)) as f64 + 2.0) * 2.0;
            }
            if f.x_tau[c] != R::ZERO && f.tau_raw[c] {
                x4 += (self.ab.na(c) as f64 + 2.0) * 2.0;
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

    /// Whether the U stage of the second vector runs its dense slot loop (`0x123690`: `m <= h33^2 c` for a spike of `c` entries, or `5 h33^2 t > m - i21` for `t` tail entries),
    /// from the spike's entry count `n` and tail-entry count `t` (`t` is read only where `tail_count_needed(n)`).
    fn tau_u_dense_nt(&self, n: usize, t: usize) -> bool {
        let (m, d) = (self.m as f64, self.h33);
        if m <= (d * d) * n as f64 {
            return true;
        }
        ((5.0 * d) * d) * t as f64 > (self.m - self.i21) as f64
    }

    /// Whether the tail-entry count of a spike of `n` entries decides the dense slot loop (`m > h33^2 n`).
    fn tail_count_needed(&self, n: usize) -> bool {
        let d = self.h33;
        self.m as f64 > (d * d) * n as f64
    }

    /// Work of the tau solve's U stage (`0x123690`, which adds to `F+0x200`): `10 + 2 n + 6 c` for `n` needed positions and a spike of
    /// `c` entries; a dense spike (`m <= h33^2 c`) or one with many tail entries (`5 h33^2 t > m - i21`) books `3 (m - i21)` and `2 len` (its tail-row entries) per
    /// tail column with quotient above the tolerance (plus `2 t` in the second case) and `4 (m - i21)` at the end; otherwise the heap walk books `4 h33` per nonzero
    /// tail column and `4` per such column at the end. `(2 n) h34` is added last. `x` holds the solved tail values.
    #[cfg(test)]
    fn tau_u_work_ref(&self, spike_n: usize, t: usize, x: &[R], n: usize) -> (f64, usize) {
        let m = self.m as f64;
        let (cnt, nf) = (spike_n as f64, n as f64);
        let d = self.h33;
        let tail_rows = (self.m - self.i21) as f64;
        let mut e = ((10.0 + nf) + nf) + 6.0 * cnt;
        let mut dense = m <= (d * d) * cnt;
        if !dense {
            let t = t as f64;
            if ((5.0 * d) * d) * t > tail_rows {
                dense = true;
                e = 2.0 * t + e;
            }
        }
        debug_assert_eq!(dense, self.tau_u_dense_nt(spike_n, t));
        let fin;
        if dense {
            e = 3.0 * tail_rows + e;
            for t in (self.i21..self.m).rev() {
                let c = self.at[t];
                if x[c].abs() > R::c(TOL) {
                    e += self.ab.na(c) as f64 * 2.0;
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
        (((nf + nf) * self.h34 + e) + 4.0 * fin, fin as usize)
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
        let (lout, count_l, e_l) = self.lstage_num(rows, vals, force_dense);
        let spike = match lout {
            LOut::List(list) => list,
            LOut::Dense(x) => {
                let mut inp = std::mem::take(&mut self.f2.get_mut().col_in);
                let mut sp = std::mem::take(&mut self.spike);
                let nz = Self::dense_lists(&x, &mut inp, &mut sp);
                let u_dense = rows.len() * self.dthr > m
                    || self.lev_dense_ft
                    || (!self.etas.is_empty() && Self::count_dense(count_l, self.dthr, m))
                    || !(self.nuc.len() + nz < self.i21 || self.dens * (nz as f64) < 2.0 * m as f64);
                if u_dense && !(force_dense && crate::dump::on("bsolve")) {
                    return self.ftran_with_dense(rows, x, inp, sp, nz, count_l, e_l);
                }
                // The sparse U routine and the bsolve dump take the list.
                Self::dense_spike(&x)
            }
        };
        // The previous column snapshot has served every consumer. Reuse it as
        // numeric work and reset only blocks that could contain a prior write.
        let old = self.f2.get_mut();
        let mut x = std::mem::take(&mut old.x_col);
        let mut dirty = old.dirty_blocks.take();
        let mut col_surv = std::mem::take(&mut old.col_surv);
        col_surv.clear();
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
            || self.lev_dense_ft
            || (!self.etas.is_empty() && Self::count_dense(count_l, self.dthr, m))
            || !(self.nuc.len() + spike.len() < self.i21
                || self.dens * (spike.len() as f64) < 2.0 * m as f64);
        let bs = force_dense && crate::dump::on("bsolve");
        if !u_dense && !bs && !self.bad_diag {
            return self.ftran_with_sparse(rows, spike, x, dirty, col_surv, count_l, e_l);
        }
        let mut touched = FtranTouched::new(m, u_dense, dirty.unwrap_or_default());
        let mut inp = std::mem::take(&mut old.col_in);
        inp.clear();
        inp.extend(spike.iter().map(|e| e.0));
        touched.insert_all(&inp);
        // The sparse solve runs the closure of the list's tail positions under the `a` lists in descending slot order, which is also the
        // tail part of the output order: `uorder_open` produces it once and `uorder_close` finishes the order after the tail loop.
        let opened = !(u_dense && !bs);
        let tail: Vec<usize> = if opened {
            let tail = self.uorder_open(&inp);
            touched.insert_all(&tail);
            tail
        } else {
            // Unreached columns contain zero and the dense kernel skips them before division.
            // Slot order preserves the arithmetic order of every reached column.
            (self.i21..m).rev().map(|s| self.at[s]).collect()
        };
        let x_l = if bs { x.clone() } else { Vec::new() };
        // The U-stage work estimate sums integer-valued lengths over the tail columns that keep a nonzero value; they are collected where
        // the column is finished (a scatter only hits lower slots) instead of by a later scan of the tail (`u_work_ref`).
        let (mut sa, mut sb, mut su) = (0usize, 0usize, 0usize);
        let tol = R::c(TOL);
        if let FtranTouched::Sparse { words, occupied } = &mut touched {
            // The sparse solve (not `u_dense`): the quotient alone is tested, and a scatter marks its target block without the first-visit test
            // (the tail dependency walk is done; the set of marks is what the readers see).
            let (diag, uc) = (&self.diag[..], &self.uc[..]);
            for &c in &tail {
                let xc = x[c] / diag[c];
                if xc.abs() <= tol {
                    x[c] = R::ZERO;
                    continue;
                }
                x[c] = xc;
                col_surv.push(c);
                let (na, nb) = self.ab.counts(c);
                sa += na;
                sb += nb;
                let col = &uc[c];
                su += col.len();
                for &(j, v) in col {
                    x[j] -= v * xc;
                    let block = j / 64;
                    let word = &mut words[block];
                    if *word == 0 {
                        occupied.push(block);
                    }
                    *word |= 1u64 << (j % 64);
                }
            }
        } else {
            for &c in &tail {
                let mut xc = x[c];
                if xc.abs() <= tol {
                    x[c] = R::ZERO;
                    continue;
                }
                xc /= self.diag[c];
                if xc.abs() <= tol {
                    x[c] = R::ZERO;
                    continue;
                }
                x[c] = xc;
                col_surv.push(c);
                sa += self.ab.na(c);
                sb += self.ab.nb(c);
                su += self.uc[c].len();
                for &(j, v) in &self.uc[c] {
                    x[j] -= v * xc;
                    touched.insert(j);
                }
            }
        }
        // Invariant (P4-2): the nonzero tail slots are exactly the processed survivors (`col_surv`).
        #[cfg(any(test, debug_assertions, feature = "factor-shadow"))]
        assert_eq!(
            (self.i21..m).filter(|&s| x[self.at[s]] != R::ZERO).count(),
            col_surv.len(),
            "ftran_with: nonzero tail slots differ from the processed survivors"
        );
        let x_t = if bs { x.clone() } else { Vec::new() };
        touched.finish();
        {
            let (kind, diag) = (&self.kind[..], &self.diag[..]);
            touched.for_each(|p| {
                if kind[p] != -2 {
                    let v = x[p] / diag[p];
                    x[p] = if (u_dense && x[p].abs() <= tol) || v.abs() <= tol { R::ZERO } else { v };
                }
            });
        }
        if bs {
            crate::dump::emit(
                "bsolve",
                serde_json::json!({"l": jv(&x_l), "tail_done": jv(&x_t), "x": jv(&x), "tail": tail, "kind": self.kind, "slot": self.slot, "diag": jv(&self.diag), "u_dense": u_dense, "uc": ju(&self.uc), "lcols": (0..self.m).map(|v| self.lcolumn(v).0.iter().map(|&y| y as usize).collect::<Vec<_>>()).collect::<Vec<_>>(), "lvals": (0..self.m).map(|v| self.lcolumn(v).1.iter().map(|e| e.to_f64()).collect::<Vec<_>>()).collect::<Vec<_>>(), "perm": self.perm, "rows": rows, "vals": jv(&vals)}),
            );
        }
        // A dense result (the reference's `-1` count) is scanned by position index; a sparse one comes in the U solve's order.
        let dense = rows.len() * self.dthr > m
            || self.lev_dense_ft
            || (!self.etas.is_empty() && Self::count_dense(count_l, self.dthr, m))
            || !(self.nuc.len() + spike.len() < self.i21
                || self.dens * (spike.len() as f64) < 2.0 * m as f64);
        let vec_dense = rows.len() * self.dthr > m
            || self.lev_dense_ft
            || (!self.etas.is_empty() && Self::count_dense(count_l, self.dthr, m));
        let u_cnt = if vec_dense { None } else { Some(spike.len()) };
        let u_est = self.u_work(&inp, &x, !dense, u_cnt, &Surv { list: &col_surv, n: col_surv.len(), sa, sb, su });
        #[cfg(test)]
        assert_eq!(u_est.to_bits(), self.u_work_ref(&inp, &x, !dense, u_cnt).to_bits(), "u_work");
        self.work = u_est + e_l;
        self.work_parts = (e_l, u_est);
        self.last_dense = dense;
        self.last_ordered = dense || u_dense;
        let mut scratch = std::mem::take(&mut self.f2.get_mut().out_scr);
        scratch.clear();
        // A dense result is scanned by position, so only the sparse output reads the order `uorder_close` completes.
        let ordered = if opened && !(dense || u_dense) {
            Some(self.uorder_close(tail, &x))
        } else {
            if opened {
                self.uorder_discard(tail);
            }
            None
        };
        if dense || u_dense {
            scratch.extend(x.iter().copied().enumerate().filter(|&(_, v)| v != R::ZERO));
            if let Some(o) = ordered {
                self.uorder_done(o);
            }
        } else {
            let ordered = ordered.unwrap();
            for &p in &ordered {
                if touched.take_sparse(p) && x[p] != R::ZERO {
                    scratch.push((p, x[p]));
                }
            }
            self.uorder_done(ordered);
            touched.for_each(|p| {
                if x[p] != R::ZERO {
                    scratch.push((p, x[p]));
                }
            });
        }
        // The result is handed out with its exact length; the scratch keeps its capacity for the next call.
        let out = scratch.clone();
        scratch.clear();
        let dirty_blocks = match touched {
            FtranTouched::Sparse { occupied, .. } => Some(occupied),
            FtranTouched::Dense(_) => None,
        };
        *self.f2.get_mut() = F2 {
            col_n: if vec_dense { None } else { Some(spike.len()) },
            col_in: inp,
            x_col: x,
            dirty_blocks,
            col_surv,
            col_sa: sa,
            col_sb: sb,
            out_scr: scratch,
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

    /// The sparse-list U stage of `ftran_with` (`u_dense` false, no `bsolve` dump): the reference's `0x1223e0` and the tail of its joint routine.
    /// `x` holds the L-stage list `spike` scattered over a cleared vector (`dirty` lists the blocks the previous column left nonzero).
    ///
    /// The tail columns are visited in descending slot order, the order of the output and of every update. Where the dense gate of
    /// `uorder_open` holds (the reference sweeps, `12c6fa`/`12c804`/`1225df`) all tail slots are swept, a slot holding zero costing a test;
    /// otherwise the reached slots come from a max-heap of slots filled by the scatters themselves, which abandons itself for a sweep of the
    /// remaining slots once the heap is dense in its range (the reference's `6.0 * size > top - loff` at `1228c1`). A column holding zero
    /// scatters nothing, so which of the two walks visits it is not observable: same updates in the same order as the closure walk of
    /// `uorder_open` (`ftran_with_prev` in `ordering_old`), and no closure, key sort or tail list.
    ///
    /// The output is the singletons in input order, the surviving tail columns, the nucleus positions (in nucleus-list order when the
    /// reference lists the whole nucleus, in discovery order otherwise) and the other touched positions in ascending position.
    #[allow(clippy::too_many_arguments)]
    fn ftran_with_sparse(
        &mut self,
        rows: &[usize],
        spike: Vec<(usize, R)>,
        mut x: Vec<R>,
        dirty: Option<Vec<usize>>,
        mut col_surv: Vec<usize>,
        count_l: Option<usize>,
        e_l: f64,
    ) -> Vec<(usize, R)> {
        let (m, i21, i3) = (self.m, self.i21, self.nuc.len());
        let tol = R::c(TOL);
        let mut inp = std::mem::take(&mut self.f2.get_mut().col_in);
        inp.clear();
        inp.extend(spike.iter().map(|e| e.0));
        let mut touched = FtranTouched::new(m, false, dirty.unwrap_or_default());
        let mut order_scratch = self.btran_ws.order.borrow_mut();
        let btran_ws::Order { heap, bits, tri, ns, .. } = &mut *order_scratch;
        let (mut sa, mut sb, mut su) = (0usize, 0usize, 0usize);
        let sweep_all;
        let nuc_in;
        {
            let FtranTouched::Sparse { words, occupied } = &mut touched else { unreachable!() };
            let (diag, uc, kind, slot, at, ab) = (&self.diag[..], &self.uc[..], &self.kind[..], &self.slot[..], &self.at[..], &self.ab);
            // Mark the L-stage positions (a repeated position counts once) and count them by class; the tail ones seed the heap.
            let (mut n, mut nuc_n) = (0usize, 0usize);
            tri.clear();
            for &p in &inp {
                let word = &mut words[p / 64];
                if *word == 0 {
                    occupied.push(p / 64);
                }
                let bit = 1u64 << (p % 64);
                if *word & bit != 0 {
                    continue;
                }
                *word |= bit;
                n += 1;
                match kind[p] {
                    -1 => {}
                    -2 => tri.push(slot[p]),
                    _ => nuc_n += 1,
                }
            }
            nuc_in = nuc_n;
            let d2 = self.dens;
            sweep_all = (m as f64) <= (n as f64) * d2 || d2 * ((nuc_n + tri.len()) as f64) > (i3 + (m - i21)) as f64;
            // One tail column: the quotient alone is tested (the sparse routine), a zero slot is stored as +0.0; `$mark` records a scatter target.
            macro_rules! solve_column {
                ($c:expr, |$j:ident| $mark:block) => {{
                    let c = $c;
                    let xr = x[c];
                    if xr == R::ZERO {
                        x[c] = R::ZERO;
                    } else {
                        let xc = xr / diag[c];
                        if xc.abs() <= tol {
                            x[c] = R::ZERO;
                        } else {
                            x[c] = xc;
                            col_surv.push(c);
                            let (na, nb) = ab.counts(c);
                            sa += na;
                            sb += nb;
                            let col = &uc[c];
                            su += col.len();
                            for &($j, v) in col {
                                x[$j] -= v * xc;
                                $mark
                            }
                        }
                    }
                }};
            }
            macro_rules! mark {
                ($j:ident) => {{
                    let word = &mut words[$j / 64];
                    if *word == 0 {
                        occupied.push($j / 64);
                    }
                    *word |= 1u64 << ($j % 64);
                }};
            }
            col_surv.clear();
            #[cfg(test)]
            SPARSE_ROUTES.with(|r| {
                let mut v = r.get();
                v[0] += 1;
                v[1] += usize::from(sweep_all);
                r.set(v);
            });
            let mut top = i21;
            if sweep_all {
                for t in (i21..m).rev() {
                    solve_column!(at[t], |j| { mark!(j) });
                }
            } else {
                heap.clear();
                heap.extend(tri.iter().copied());
                while let Some(&t) = heap.peek() {
                    if heap.len() > 1 && t - i21 < SWEEP_RATIO * heap.len() {
                        top = t + 1;
                        #[cfg(test)]
                        SPARSE_ROUTES.with(|r| {
                            let mut v = r.get();
                            v[3] += 1;
                            r.set(v);
                        });
                        break;
                    }
                    heap.pop();
                    solve_column!(at[t], |j| {
                        let word = &mut words[j / 64];
                        if *word == 0 {
                            occupied.push(j / 64);
                        }
                        let bit = 1u64 << (j % 64);
                        if *word & bit == 0 {
                            *word |= bit;
                            if kind[j] == -2 {
                                heap.push(slot[j]);
                            }
                        }
                    });
                }
                #[cfg(test)]
                if top == i21 {
                    SPARSE_ROUTES.with(|r| {
                        let mut v = r.get();
                        v[2] += 1;
                        r.set(v);
                    });
                }
                heap.clear();
                for t in (i21..top).rev() {
                    solve_column!(at[t], |j| { mark!(j) });
                }
            }
        }
        // Invariant (P4-2): the nonzero tail slots are exactly the processed survivors (`col_surv`).
        #[cfg(any(test, debug_assertions, feature = "factor-shadow"))]
        assert_eq!(
            (self.i21..m).filter(|&s| x[self.at[s]] != R::ZERO).count(),
            col_surv.len(),
            "ftran_with: nonzero tail slots differ from the processed survivors"
        );
        touched.finish();
        // The nucleus positions come out in nucleus-list order when the reference lists the whole nucleus (the dense gate, or the survivors'
        // B lists reach more than 0.4 of it), otherwise in discovery order.
        let sorted_nuc = sweep_all || (nuc_in + sb) as f64 * 2.5 > i3 as f64;
        #[cfg(test)]
        SPARSE_ROUTES.with(|r| {
            let mut v = r.get();
            v[if sorted_nuc { 4 } else { 5 }] += 1;
            r.set(v);
        });
        ns.clear();
        {
            let (kind, diag) = (&self.kind[..], &self.diag[..]);
            touched.for_each(|p| {
                if kind[p] != -2 {
                    let v = x[p] / diag[p];
                    let w = if v.abs() <= tol { R::ZERO } else { v };
                    x[p] = w;
                    if sorted_nuc && kind[p] >= 0 && w != R::ZERO {
                        ns.push(p);
                    }
                }
            });
        }
        let vec_dense = rows.len() * self.dthr > m
            || self.lev_dense_ft
            || (!self.etas.is_empty() && Self::count_dense(count_l, self.dthr, m));
        debug_assert!(
            !(vec_dense || !(self.nuc.len() + spike.len() < self.i21 || self.dens * (spike.len() as f64) < 2.0 * m as f64)),
            "the sparse U route has a sparse result"
        );
        let u_cnt = if vec_dense { None } else { Some(spike.len()) };
        let u_est = self.u_work(&inp, &x, true, u_cnt, &Surv { list: &col_surv, n: col_surv.len(), sa, sb, su });
        #[cfg(test)]
        assert_eq!(u_est.to_bits(), self.u_work_ref(&inp, &x, true, u_cnt).to_bits(), "u_work");
        self.work = u_est + e_l;
        self.work_parts = (e_l, u_est);
        self.last_dense = false;
        self.last_ordered = false;
        let mut scratch = std::mem::take(&mut self.f2.get_mut().out_scr);
        scratch.clear();
        {
            let kind = &self.kind[..];
            for &p in &inp {
                if kind[p] == -1 && touched.take_sparse(p) && x[p] != R::ZERO {
                    scratch.push((p, x[p]));
                }
            }
            for &c in &col_surv {
                let marked = touched.take_sparse(c);
                debug_assert!(marked, "a surviving tail column was never marked");
                scratch.push((c, x[c]));
            }
            if sorted_nuc {
                Self::sort_distinct_keys(bits, ns, i3, false, |p| kind[p] as usize, |k| self.nuc[k]);
                for &p in ns.iter() {
                    touched.take_sparse(p);
                    scratch.push((p, x[p]));
                }
            } else {
                for &p in &inp {
                    if kind[p] >= 0 && touched.take_sparse(p) && x[p] != R::ZERO {
                        scratch.push((p, x[p]));
                    }
                }
                for &c in &col_surv {
                    for &y in self.ab.b(c) {
                        let y = y as usize;
                        if touched.take_sparse(y) && x[y] != R::ZERO {
                            scratch.push((y, x[y]));
                        }
                    }
                }
            }
            touched.for_each(|p| {
                if x[p] != R::ZERO {
                    scratch.push((p, x[p]));
                }
            });
        }
        drop(order_scratch);
        // The result is handed out with its exact length; the scratch keeps its capacity for the next call.
        let out = scratch.clone();
        scratch.clear();
        let FtranTouched::Sparse { occupied, .. } = touched else { unreachable!() };
        *self.f2.get_mut() = F2 {
            col_n: u_cnt,
            col_in: inp,
            x_col: x,
            dirty_blocks: Some(occupied),
            col_surv,
            col_sa: sa,
            col_sb: sb,
            out_scr: scratch,
            ..F2::<R>::default()
        };
        // The column replacement receives the spike in L-stage order after a sparse solve (ascending after a dense L-stage vector).
        let mut spike = spike;
        if vec_dense {
            spike.sort_unstable_by_key(|e| e.0);
        }
        spike.retain(|e| e.1.abs() > R::c(TOL));
        self.spike = spike;
        self.spike_used = false;
        self.last_alpha = out.len();
        out
    }

    /// `ftran_with` for a dense L-stage vector `x` (`nz` entries with `v != 0`, every other slot +0.0) that takes the dense U routine
    /// (`0x121ef0`): `x` is the U work vector and then the column snapshot, with no spike list, re-scatter or second zero-filled vector.
    /// The tail loop and the final pass are those of `ftran_with` for `u_dense`; a position the old marker set did not contain holds
    /// +0.0 and comes out as +0.0 from the final pass.
    #[allow(clippy::too_many_arguments)]
    fn ftran_with_dense(
        &mut self,
        rows: &[usize],
        mut x: Vec<R>,
        inp: Vec<usize>,
        spike: Vec<(usize, R)>,
        nz: usize,
        count_l: Option<usize>,
        e_l: f64,
    ) -> Vec<(usize, R)> {
        let m = self.m;
        // `inp` (the L-stage entries' positions) and `spike` (those above the tolerance, in ascending position as the column replacement
        // receives them after a dense solve) come from `dense_lists`.
        let old = self.f2.get_mut();
        let mut col_surv = std::mem::take(&mut old.col_surv);
        let mut scratch = std::mem::take(&mut old.out_scr);
        col_surv.clear();
        let (mut sa, mut sb, mut su) = (0usize, 0usize, 0usize);
        let (i21, tol) = (self.i21, R::c(TOL));
        // The tail columns, one lookup of each column's lists, walked in descending slot order: test the raw value, divide, test the
        // quotient, store, count, scatter over the whole U column.
        {
            let (diag, uc) = (&self.diag[..], &self.uc[..]);
            col_surv.reserve(m - i21);
            for &c in self.at[i21..m].iter().rev() {
                let mut xc = x[c];
                if xc.abs() <= tol {
                    x[c] = R::ZERO;
                    continue;
                }
                xc /= diag[c];
                if xc.abs() <= tol {
                    x[c] = R::ZERO;
                    continue;
                }
                x[c] = xc;
                col_surv.push(c);
                let (na, nb) = self.ab.counts(c);
                sa += na;
                sb += nb;
                let col = &uc[c];
                su += col.len();
                for &(j, v) in col {
                    x[j] -= v * xc;
                }
            }
        }
        #[cfg(any(test, debug_assertions, feature = "factor-shadow"))]
        assert_eq!(
            (self.i21..m).filter(|&s| x[self.at[s]] != R::ZERO).count(),
            col_surv.len(),
            "ftran_with: nonzero tail slots differ from the processed survivors"
        );
        // The tail slots are final; every other slot is zeroed at or below the tolerance and divided otherwise (the old marker set held
        // the L entries and the scatter targets; every other slot is +0.0 and stays so). The pass also publishes the nonzero slots in
        // ascending position into the retained list, which is then copied out at its exact length.
        #[cfg(any(test, debug_assertions, feature = "factor-shadow"))]
        assert!((0..m).all(|p| (self.kind[p] == -2) == (self.slot[p] >= self.i21)), "tail kind and slot differ");
        scratch.clear();
        scratch.reserve(m);
        let mut count = 0usize;
        {
            let dst = scratch.as_mut_ptr();
            let (kind, diag) = (&self.kind[..], &self.diag[..]);
            for (p, xp) in x.iter_mut().enumerate() {
                let v = *xp;
                let w = if v.abs() <= tol {
                    R::ZERO
                } else if kind[p] != -2 {
                    let q = v / diag[p];
                    if q.abs() <= tol { R::ZERO } else { q }
                } else {
                    v
                };
                *xp = w;
                if w != R::ZERO {
                    // `count <= p < m <= capacity`; the entry type is Copy.
                    unsafe { dst.add(count).write((p, w)) };
                    count += 1;
                }
            }
        }
        unsafe { scratch.set_len(count) };
        let vec_dense = rows.len() * self.dthr > m
            || self.lev_dense_ft
            || (!self.etas.is_empty() && Self::count_dense(count_l, self.dthr, m));
        let u_cnt = if vec_dense { None } else { Some(nz) };
        let u_est = self.u_work(&inp, &x, false, u_cnt, &Surv { list: &col_surv, n: col_surv.len(), sa, sb, su });
        #[cfg(test)]
        assert_eq!(u_est.to_bits(), self.u_work_ref(&inp, &x, false, u_cnt).to_bits(), "u_work");
        self.work = u_est + e_l;
        self.work_parts = (e_l, u_est);
        self.last_dense = true;
        self.last_ordered = true;
        let out = scratch.clone();
        scratch.clear();
        *self.f2.get_mut() = F2 {
            col_n: u_cnt,
            col_in: inp,
            x_col: x,
            dirty_blocks: None,
            col_surv,
            col_sa: sa,
            col_sb: sb,
            out_scr: scratch,
            ..F2::<R>::default()
        };
        self.spike = spike;
        self.spike_used = false;
        self.last_alpha = out.len();
        out
    }

    /// The diagonal of the fresh factor (the first part of `numerics`, without building L).
    pub fn diag(&self) -> &[R] {
        &self.diag
    }

    /// The fresh factor's numerics for the `luout` dump: the diagonal and each position's L column `(position, value)`.
    pub fn numerics(&self) -> (&[R], Vec<Vec<(usize, R)>>) {
        let l = self.lcolumns();
        (&self.diag, l)
    }

    /// Fresh U columns in the stored solve order, indexed by pivot position.
    pub fn u_numerics(&self) -> &[Vec<(usize, R)>] {
        &self.uc
    }

    /// The reference's refactorization decision, evaluated before an iteration's own update: at the update limit
    /// `max(20, m^0.7)`, or from `U >= m^0.45` when the accumulated work `(1.5 S + 2 c) m` exceeds `min((Z - 1.5 m)^2, 4 m^2 U0)`.
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
        let urows = &self.urows;
        out.clear();
        let mut eu = 0.0f64;
        let push = |w: &mut W, p: usize, d: R| {
            let (cols, vals) = urows.row(p);
            for (&c, &v) in cols.iter().zip(vals) {
                let value = w.entry(c as usize, &self.slot);
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
                    eu = urows.len(p) as f64 * 2.0 + eu;
                }
                push(w, p, d);
                for &c in urows.cols(p) {
                    lowest = lowest.min(self.slot[c as usize]);
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
                eu = urows.len(p) as f64 * 2.0 + eu;
                out.push(p);
                x[p] = d;
                push(w, p, d);
            }
        }
        (out, eu)
    }

    /// The dense U^T loop of the BTRAN over the slots in order; `e` is the work estimate so far, returned with the rows' share added.
    fn dense_btran_u(&self, x: &mut [R], mut e: f64) -> f64 {
        for t in 0..self.m {
            let p = self.at[t];
            let xv = x[p];
            let d = xv / self.diag[p];
            if xv.abs() <= R::c(TOL) || d.abs() <= R::c(TOL) {
                x[p] = R::ZERO;
                continue;
            }
            x[p] = d;
            e = self.urows.len(p) as f64 * 2.0 + e;
            let (cols, vals) = self.urows.row(p);
            for (&c, &v) in cols.iter().zip(vals) {
                x[c as usize] -= v * d;
            }
        }
        e
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
        let urows = &self.urows;
        let mut x = std::mem::take(&mut ws.x);
        x.resize(m, R::ZERO);
        let mut e = 0.0f64;
        let cnt = inp.len();
        self.bt_dense = false;
        if direct && inp.len() == 1 {
            self.row_cached = Some(inp[0].0);
        }
        let hsum = (self.h33 + self.h34) as usize;
        let direct = direct && !inp.is_empty() && !urows.is_empty(inp[0].0);
        // The current vector: a list (`Some`) or dense (`None`).
        let mut list: Option<Vec<usize>>;
        if !dense_in && (direct || hsum * cnt <= m) {
            let order = std::mem::take(&mut ws.input_order);
            let (out, eu) = if self.bt_occupied && cnt <= m / 64 {
                ws.tail.reset(m);
                self.sparse_btran_u(inp, &mut x, &mut ws.tail, order)
            } else {
                ws.dense_tail.resize(m, R::ZERO);
                self.sparse_btran_u(inp, &mut x, &mut ws.dense_tail, order)
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
            e = self.dense_btran_u(&mut x, e);
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
                let mut w = std::mem::take(&mut ws.partial);
                w.resize(m, R::ZERO);
                let mut maxi = 0usize;
                for &i in l {
                    let v = x[i];
                    if i < nseq {
                        if i < n1 {
                            w[i] += v;
                        } else {
                            res.push(i);
                            let (src, mult) = lt.col(i);
                            el = src.len() as f64 * 2.0 + el;
                            for (&j, &lv) in src.iter().zip(mult) {
                                w[j as usize] -= lv * v;
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
                for i in (nseq..maxi + 1).rev().chain((0..n1).rev()) {
                    let v = w[i];
                    if v.abs() > R::c(TOL) {
                        let (src, mult) = lt.col(i);
                        el = src.len() as f64 * 2.0 + el;
                        res.push(i);
                        for (&j, &lv) in src.iter().zip(mult) {
                            w[j as usize] -= lv * v;
                        }
                    } else {
                        w[i] = R::ZERO;
                    }
                }
                x.fill(R::ZERO);
                ws.partial = x;
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
                            for &y in self.lcolumn(v).0 {
                                lev[v] = lev[v].max(lev[y as usize] + 1);
                            }
                        }
                        let maxlev = lev.iter().copied().max().unwrap_or(0);
                        (lev, maxlev)
                    });
                    let maxlev = *maxlev;
                    ws.buckets.resize_with(maxlev + 1, Vec::new);
                    let buckets = &mut ws.buckets;
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
                                let (src, mult) = lt.col(i);
                                for (&j, &lv) in src.iter().zip(mult) {
                                    x[j as usize] -= lv * v;
                                }
                                // The reference walks the row from its last entry.
                                for &j in src.iter().rev() {
                                    let j = j as usize;
                                    if ws.marks[j] != stamp {
                                        ws.marks[j] = stamp;
                                        buckets[lev[j]].push(j);
                                    }
                                }
                            } else {
                                x[i] = R::ZERO;
                            }
                        }
                        buckets[lv_i].clear();
                        lv_i += 1;
                    }
                } else {
                    let heap = &mut ws.pending;
                    for &i in l {
                        ws.marks[i] = stamp;
                        heap.push(i);
                    }
                    while let Some(i) = heap.pop() {
                        let v = x[i];
                        if v.abs() >= R::c(TOL) {
                            res.push(i);
                            let (src, mult) = lt.col(i);
                            for (&j, &lv) in src.iter().zip(mult) {
                                let j = j as usize;
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
            if e.fract() == 0.0 && e < 4503599627370496.0 {
                // `e` is integer-valued (below 2^52): adding the doubled row lengths at once is exact.
                let len_sum = pool::reverse_sweep(lt, &mut x);
                e = len_sum as f64 * 2.0 + e;
            } else {
                for i in (0..m).rev() {
                    let v = x[i];
                    if v.abs() > R::c(TOL) {
                        let (src, mult) = lt.col(i);
                        e = src.len() as f64 * 2.0 + e;
                        for (&j, &lv) in src.iter().zip(mult) {
                            x[j as usize] -= lv * v;
                        }
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
                Self::dense_btran_output(&x, &self.perm[..m], self.lists.take())
            } else {
                let mut out = self.lists.take();
                out.extend((0..m).filter_map(|r| {
                    let v = x[self.perm[r]];
                    (v != R::ZERO).then_some((r, v))
                }));
                out
            }
        } else {
            let mut out = self.lists.take();
            out.extend(res.iter().copied().filter(|&i| x[i] != R::ZERO).map(|i| (row_of[i], x[i])));
            out
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

    /// The sources of step `v` in the forward reset list of the pricing weights (the reference's `F+0x70/0x80`): the steps whose L column holds
    /// the row of step `v`, ascending (the reference walks the list descending).
    #[inline(always)]
    pub fn reset_forward(&self, v: usize) -> &[u32] {
        self.lt.col(v).0
    }

    /// The tail columns of row `y` of the numeric U rows, in physical order (the reference's backward reset list `F+0xb0/0xd0`). Valid for the
    /// pricing reset only on a factor without updates.
    #[inline(always)]
    pub fn reset_backward(&self, y: usize) -> &[u32] {
        self.urows.cols(y)
    }

    /// `nnz(U)` of the fresh factor (the reference's `F+0xb0` pointer array total): the factor size less `nnz(L)` and the two diagonals.
    pub fn u_nnz(&self) -> usize {
        self.znz as usize - self.nnz_l - 2 * self.m
    }

    /// Whether the fresh factor's L solve goes level by level (`0x12a990`).
    pub fn has_levels(&self) -> bool {
        self.lev_built
    }

    /// Whether the most recent `ftran` returned a dense result (the reference's count `-1`).
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
        let mut scratch = self.btran_ws.dense.borrow_mut();
        let btran_ws::Dense { x, xs } = &mut *scratch;
        x.resize(m, R::ZERO);
        for r in 0..m {
            x[self.perm[r]] += b[r];
        }
        let mut e = (m - self.loff) as f64 * 2.0 + self.n1 as f64 + self.n1 as f64;
        // `e` is integer-valued here, so the doubled column lengths are added at once.
        debug_assert!(e.fract() == 0.0);
        let len_sum = pool::forward_sweep(&self.lpool, x);
        e = len_sum as f64 * 2.0 + e;
        if !self.etas.is_empty() {
            e = self.sum_eta * 2.0 + e;
            for eta in &self.etas {
                xs.clear();
                xs.extend(eta.ind.iter().map(|&j| x[j]));
                let sv = dot(xs, &eta.val, eta.peel);
                x[eta.idx] = if sv.abs() <= R::c(TOL) { R::ZERO } else { sv };
            }
        }
        let (mut n_surv, mut su) = (0usize, 0usize);
        for t in (self.i21..m).rev() {
            let c = self.at[t];
            let xc = x[c] / self.diag[c];
            if x[c].abs() <= R::c(TOL) || xc.abs() <= R::c(TOL) {
                x[c] = R::ZERO;
                continue;
            }
            x[c] = xc;
            n_surv += 1;
            su += self.uc[c].len();
            for &(j, v) in &self.uc[c] {
                x[j] -= v * xc;
            }
        }
        let result = self.u_work(&[], x, false, None, &Surv { list: &[], n: n_surv, sa: 0, sb: 0, su }) + e;
        #[cfg(test)]
        assert_eq!(result.to_bits(), (self.u_work_ref(&[], x, false, None) + e).to_bits(), "dense u_work");
        x.fill(R::ZERO);
        result
    }

    /// `u_work_ref` from the sums of the U loop: `s` describes the tail columns that kept a nonzero value (count, their `a`, `b` and
    /// `uc` lengths, and for the sparse solve the columns in the solve's descending slot order). Every summed term is an integer, so
    /// the sums equal the reference's scans bit for bit; the heap walk still decides the walk's own terms from the solved values.
    fn u_work(&self, inp: &[usize], x: &[R], sparse: bool, cnt: Option<usize>, s: &Surv<'_>) -> f64 {
        let m = self.m;
        let i21 = self.i21;
        let i3 = self.nuc.len();
        if !sparse {
            let mut e = 0.0f64;
            if let Some(n) = cnt {
                e = m as f64 + n as f64 * 2.0 + e;
            }
            e = (m - i21) as f64 * 3.0 + e;
            e += (2 * s.su + 4 * s.n) as f64;
            return 3.0 * i21.saturating_sub(1) as f64 + e;
        }
        let n = inp.len();
        let d2 = self.dens;
        if (m as f64) <= n as f64 * d2 {
            let e = (m - i21) as f64 * 3.0 + 0.0 + (2 * s.su) as f64;
            return i3 as f64 + i3 as f64 + e;
        }
        let mut scratch = self.btran_ws.work.borrow_mut();
        let btran_ws::Work { mark, nm, tail, heap, .. } = &mut *scratch;
        tail.clear();
        tail.extend(inp.iter().copied().filter(|&p| self.kind[p] == -2));
        let nuc_in = inp.iter().filter(|&&p| self.kind[p] >= 0).count();
        if d2 * ((nuc_in + tail.len()) as f64) > (i3 + (m - i21)) as f64 {
            let e = (2 * (nuc_in + tail.len())) as f64;
            let e = (m - i21) as f64 * 3.0 + e + (2 * s.su) as f64;
            tail.clear();
            return i3 as f64 + i3 as f64 + e;
        }
        let mut e = (tail.len() * 2) as f64;
        mark.prepare(m);
        heap.clear();
        for &p in tail.iter() {
            mark.insert(p);
            heap.push(self.slot[p]);
        }
        while let Some(&top) = heap.peek() {
            let size = heap.len();
            if size > 1 && ((top - i21) as f64) < 6.0 * size as f64 {
                // Remaining slots are close together: the reference scans the range from the top slot down.
                e = (top - i21) as f64 * 5.0 + size as f64 + size as f64 + e;
                break;
            }
            heap.pop();
            let c = self.at[top];
            if x[c] != R::ZERO {
                for &j in self.ab.a(c) {
                    let j = j as usize;
                    if mark.insert(j) {
                        heap.push(self.slot[j]);
                    }
                }
            }
            e += 8.0;
        }
        let (lb, ld) = (nuc_in + s.sb, nuc_in + s.sa);
        e = ld as f64 * 3.0 + e;
        let result = if lb as f64 * 2.5 <= i3 as f64 {
            nm.prepare(m);
            let mut cntn = 0usize;
            for &p in inp.iter().filter(|&&p| self.kind[p] >= 0) {
                if nm.insert(p) {
                    cntn += 1;
                }
            }
            e = e + (2 * s.sb + 4 * s.list.len()) as f64;
            for &c in s.list {
                for &y in self.ab.b(c) {
                    let y = y as usize;
                    if nm.insert(y) {
                        cntn += 1;
                    }
                }
            }
            e + 3.0 * cntn as f64
        } else {
            e += nuc_in as f64;
            e = e + (2 * s.sb + 4 * s.list.len()) as f64;
            e + 3.0 * i3 as f64
        };
        scratch.clear();
        result
    }

    /// Work estimate of the U stage: the sparse solve `0x1223e0` (`sparse`) or the dense one `0x121ef0`. `x` holds the solved
    /// values, so a tail column survived the tolerance test when its entry is nonzero. `cnt` is the input count when the input
    /// is a sparse list.
    #[cfg(test)]
    fn u_work_ref(&self, inp: &[usize], x: &[R], sparse: bool, cnt: Option<usize>) -> f64 {
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
        let mut scratch = self.btran_ws.work.borrow_mut();
        let btran_ws::Work { mark, nm, tail, surv, heap } = &mut *scratch;
        tail.clear();
        tail.extend(inp.iter().copied().filter(|&p| self.kind[p] == -2));
        let nuc_in = inp.iter().filter(|&&p| self.kind[p] >= 0).count();
        if d2 * ((nuc_in + tail.len()) as f64) > (i3 + (m - i21)) as f64 {
            let e = (2 * (nuc_in + tail.len())) as f64;
            let e = full_tail(e);
            tail.clear();
            return i3 as f64 + i3 as f64 + e;
        }
        let mut e = (tail.len() * 2) as f64;
        mark.prepare(m);
        heap.clear();
        for &p in tail.iter() {
            mark.insert(p);
            heap.push(self.slot[p]);
        }
        let (mut lb, mut ld) = (nuc_in, nuc_in);
        surv.clear();
        let visit = |c: usize, lb: &mut usize, ld: &mut usize, surv: &mut Vec<usize>| -> bool {
            if x[c] != R::ZERO {
                *lb += self.ab.nb(c);
                *ld += self.ab.na(c);
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
                    visit(self.at[t], &mut lb, &mut ld, surv);
                }
                break;
            }
            heap.pop();
            let c = self.at[top];
            if visit(c, &mut lb, &mut ld, surv) {
                for &j in self.ab.a(c) {
                    let j = j as usize;
                    if mark.insert(j) {
                        heap.push(self.slot[j]);
                    }
                }
            }
            e += 8.0;
        }
        e = ld as f64 * 3.0 + e;
        let result = if lb as f64 * 2.5 <= i3 as f64 {
            nm.prepare(m);
            let mut cntn = 0usize;
            for &p in inp.iter().filter(|&&p| self.kind[p] >= 0) {
                if nm.insert(p) {
                    cntn += 1;
                }
            }
            for &c in surv.iter() {
                e = e + self.ab.nb(c) as f64 * 2.0 + 4.0;
                for &y in self.ab.b(c) {
                    let y = y as usize;
                    if nm.insert(y) {
                        cntn += 1;
                    }
                }
            }
            e + 3.0 * cntn as f64
        } else {
            e += nuc_in as f64;
            for &c in surv.iter() {
                e = e + self.ab.nb(c) as f64 * 2.0 + 4.0;
            }
            e + 3.0 * i3 as f64
        };
        scratch.clear();
        result
    }

    /// Sorts distinct positions into descending slot order (a tail list built with a visited marker has no repeats), by one bit per slot
    /// for a long list and by the sorted slots for a short one: no counting pass over all `m` slots.
    fn sort_tail_desc(&self, positions: &mut [usize]) {
        let mut scratch = self.btran_ws.order.borrow_mut();
        Self::sort_distinct_keys(&mut scratch.bits, positions, self.m, true, |p| self.slot[p], |s| self.at[s]);
    }

    /// The previous tail sort (oracle of the tests; its counting-sort branch for long lists gave the same order): descending slot order.
    /// `slot` and `at` are inverse permutations, so sorting slots and mapping back also preserves repeated positions exactly.
    #[cfg(test)]
    fn sort_slot_desc(&self, positions: &mut [usize]) {
        for p in positions.iter_mut() {
            *p = self.slot[*p];
        }
        positions.sort_unstable_by(|a, b| b.cmp(a));
        for p in positions {
            *p = self.at[*p];
        }
    }

    /// Sorts the distinct keys `key(v)` of `v` ascending (descending with `desc`) and maps each back through `val`, which must
    /// invert `key`. The keys must be distinct and below `nkeys`, so the order does not depend on the sort algorithm. A long list is
    /// ordered by setting one bit per key and reading the words back (the words are left zero).
    fn sort_distinct_keys(
        bits: &mut Vec<u64>,
        v: &mut [usize],
        nkeys: usize,
        desc: bool,
        key: impl Fn(usize) -> usize,
        val: impl Fn(usize) -> usize,
    ) {
        if v.len() < 24 {
            for p in v.iter_mut() {
                *p = key(*p);
            }
            if desc {
                v.sort_unstable_by(|a, b| b.cmp(a));
            } else {
                v.sort_unstable();
            }
            for p in v.iter_mut() {
                *p = val(*p);
            }
            return;
        }
        if bits.len() < nkeys.div_ceil(64) {
            bits.resize(nkeys.div_ceil(64), 0);
        }
        let (mut lo, mut hi) = (usize::MAX, 0);
        for &p in v.iter() {
            let k = key(p);
            let w = k >> 6;
            #[cfg(any(test, debug_assertions, feature = "factor-shadow"))]
            assert!(bits[w] >> (k & 63) & 1 == 0, "duplicate key");
            bits[w] |= 1 << (k & 63);
            lo = lo.min(w);
            hi = hi.max(w);
        }
        let mut i = 0;
        if desc {
            for w in (lo..=hi).rev() {
                let mut word = std::mem::take(&mut bits[w]);
                while word != 0 {
                    let bit = 63 - word.leading_zeros() as usize;
                    word &= !(1 << bit);
                    v[i] = val(w << 6 | bit);
                    i += 1;
                }
            }
        } else {
            for w in lo..=hi {
                let mut word = std::mem::take(&mut bits[w]);
                while word != 0 {
                    let bit = word.trailing_zeros() as usize;
                    word &= word - 1;
                    v[i] = val(w << 6 | bit);
                    i += 1;
                }
            }
        }
        debug_assert_eq!(i, v.len());
    }

    /// Order of the nonzero positions of the solution, given the L-stage list. The result is the retained `Order::out`, taken from the
    /// scratch; hand it back with `uorder_done`. Kept as the composition of `uorder_open` and `uorder_close` (the solve itself calls the
    /// two halves around its tail loop).
    #[cfg(test)]
    fn uorder(&self, inp: &[usize], x: &[R]) -> Vec<usize> {
        let tail = self.uorder_open(inp);
        self.uorder_close(tail, x)
    }

    /// First half of `uorder`, which needs no solution values: classifies the L-stage list and returns the tail columns the U solve
    /// processes, the closure of the list's tail positions under the `a` lists, in descending slot order. This is the order the solve's
    /// tail loop runs in and also the tail part of the output order (slots are distinct, so the heap and the key sort agree). Hand the
    /// vector back to `uorder_close`.
    fn uorder_open(&self, inp: &[usize]) -> Vec<usize> {
        let m = self.m;
        let mut scratch = self.btran_ws.order.borrow_mut();
        let btran_ws::Order { seen, have, nm, heap, bits, inp_d, sing, tri, nuc, order, dense, .. } = &mut *scratch;
        seen.prepare(m);
        have.prepare(m);
        nm.prepare(m);
        inp_d.clear();
        sing.clear();
        tri.clear();
        nuc.clear();
        for &x in inp {
            if seen.insert(x) {
                inp_d.push(x);
                match self.kind[x] {
                    -1 => sing.push(x),
                    -2 => tri.push(x),
                    _ => nuc.push(x),
                }
            }
        }
        let n = inp_d.len();
        let i3 = self.nuc.len();
        let d2 = self.dens;
        *dense = (m as f64) <= (n as f64) * d2
            || d2 * ((nuc.len() + tri.len()) as f64) > (i3 + (m - self.i21)) as f64;
        if *dense {
            let mut t = std::mem::take(order);
            t.clear();
            t.extend_from_slice(tri);
            for &x in tri.iter() {
                have.insert(x);
            }
            let mut k = 0;
            while k < t.len() {
                let c = t[k];
                k += 1;
                for &j in self.ab.a(c) {
                    let j = j as usize;
                    if have.insert(j) {
                        t.push(j);
                    }
                }
            }
            // The positions of `t` are distinct, hence so are their slots (`slot` and `at` are inverse permutations).
            Self::sort_distinct_keys(bits, &mut t, m, true, |p| self.slot[p], |s| self.at[s]);
            t
        } else {
            // Slots are distinct, so the heap of slots pops in the same order as the (slot, position) pairs.
            heap.clear();
            for &x in tri.iter() {
                have.insert(x);
                heap.push(self.slot[x]);
            }
            let mut order = std::mem::take(order);
            order.clear();
            while let Some(s) = heap.pop() {
                let k = self.at[s];
                order.push(k);
                for &j in self.ab.a(k) {
                    let j = j as usize;
                    if have.insert(j) {
                        heap.push(self.slot[j]);
                    }
                }
            }
            order
        }
    }

    /// Second half of `uorder`: the singletons, the tail columns (`tail`, from `uorder_open`) and the nucleus positions they reach. `x`
    /// holds the solution after the tail loop (the reference's counter and nucleus-list discovery only see the tail columns that
    /// survive the tolerance test).
    fn uorder_close(&self, tail: Vec<usize>, x: &[R]) -> Vec<usize> {
        let mut scratch = self.btran_ws.order.borrow_mut();
        let btran_ws::Order { nm, bits, sing, nuc, ns, out, dense, order, .. } = &mut *scratch;
        let i3 = self.nuc.len();
        // The nucleus positions are in nucleus-list order (`kind` is the index into `self.nuc`, so the keys are distinct).
        let nuc_sorted = |bits: &mut Vec<u64>, v: &mut Vec<usize>| {
            Self::sort_distinct_keys(
                bits,
                v,
                i3,
                false,
                |x| {
                    debug_assert!(self.kind[x] >= 0, "nucleus key of a non-nucleus position");
                    self.kind[x] as usize
                },
                |k| self.nuc[k],
            );
        };
        out.clear();
        out.extend_from_slice(sing);
        ns.clear();
        if *dense {
            // Every tail column counts, and the sorted nucleus list does not depend on the discovery order.
            for y in nuc.iter().copied().chain(tail.iter().flat_map(|&c| self.ab.b(c).iter().map(|&y| y as usize))) {
                if nm.insert(y) {
                    ns.push(y);
                }
            }
            out.extend_from_slice(&tail);
            nuc_sorted(bits, ns);
            out.extend_from_slice(ns);
            *order = tail;
        } else {
            let mut lb = nuc.len();
            for &k in tail.iter().filter(|&&k| x[k] != R::ZERO) {
                lb += self.ab.nb(k);
            }
            for y in nuc.iter().copied().chain(
                tail.iter()
                    .filter(|&&k| x[k] != R::ZERO)
                    .flat_map(|&k| self.ab.b(k).iter().map(|&y| y as usize)),
            ) {
                if nm.insert(y) {
                    ns.push(y);
                }
            }
            out.extend_from_slice(&tail);
            if (lb as f64) * 2.5 > i3 as f64 {
                nuc_sorted(bits, ns);
            }
            out.extend_from_slice(ns);
            *order = tail;
        }
        scratch.clear();
        std::mem::take(&mut scratch.out)
    }

    /// Hands the tail buffer of `uorder_open` back to the scratch when the output order is not read (the `uorder_close` bookkeeping
    /// without the singleton, tail and nucleus lists).
    fn uorder_discard(&self, tail: Vec<usize>) {
        let mut scratch = self.btran_ws.order.borrow_mut();
        scratch.order = tail;
        scratch.clear();
    }

    /// Returns the buffer of `uorder` to the scratch.
    fn uorder_done(&self, mut out: Vec<usize>) {
        out.clear();
        self.btran_ws.order.borrow_mut().out = out;
    }

    /// The previous `uorder`, the oracle of the differential test.
    #[cfg(test)]
    fn uorder_ref(&self, inp: &[usize], x: &[R]) -> Vec<usize> {
        let m = self.m;
        let mut scratch = self.btran_ws.order.borrow_mut();
        let btran_ws::Order { seen, have, nm, .. } = &mut *scratch;
        seen.prepare(m);
        have.prepare(m);
        nm.prepare(m);
        let mut inp_d: Vec<usize> = Vec::with_capacity(inp.len());
        for &x in inp {
            if seen.insert(x) {
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
        let dense = (m as f64) <= (n as f64) * d2
            || d2 * ((nuc.len() + tri.len()) as f64) > (i3 + (m - self.i21)) as f64;
        if dense {
            let mut t: Vec<usize> = tri.clone();
            for &x in &tri {
                have.insert(x);
            }
            let mut k = 0;
            while k < t.len() {
                let c = t[k];
                k += 1;
                for &j in self.ab.a(c) {
                    let j = j as usize;
                    if have.insert(j) {
                        t.push(j);
                    }
                }
            }
            let mut ns: Vec<usize> = Vec::new();
            for x in nuc.iter().copied().chain(t.iter().flat_map(|&c| self.ab.b(c).iter().map(|&y| y as usize))) {
                if nm.insert(x) {
                    ns.push(x);
                }
            }
            by_slot_desc(&mut t);
            let mut out = sing;
            out.extend(t);
            out.extend(nuc_sorted(ns));
            scratch.clear();
            return out;
        }
        let mut heap: BinaryHeap<(usize, usize)> = BinaryHeap::new();
        for &x in &tri {
            have.insert(x);
            heap.push((self.slot[x], x));
        }
        let mut order: Vec<usize> = Vec::new();
        let mut lb = nuc.len();
        while let Some((_, k)) = heap.pop() {
            order.push(k);
            // The reference's counter and nucleus-list discovery only see tail columns that survive the tolerance test.
            if x[k] != R::ZERO {
                lb += self.ab.nb(k);
            }
            for &j in self.ab.a(k) {
                let j = j as usize;
                if have.insert(j) {
                    heap.push((self.slot[j], j));
                }
            }
        }
        let mut lst: Vec<usize> = Vec::new();
        for y in nuc.iter().copied().chain(
            order
                .iter()
                .filter(|&&k| x[k] != R::ZERO)
                .flat_map(|&k| self.ab.b(k).iter().map(|&y| y as usize)),
        ) {
            if nm.insert(y) {
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
        scratch.clear();
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
    fn walk(&mut self, r: usize, hole: usize, mut cur: usize, refcols: &[(usize, usize)]) {
        // Column lengths leave out row r, which the reference has already removed from U. `r` occurs only in the lists of the tail
        // columns of `rowl[r]` (`refcols`, sorted by column), so any other tail column has its full length; the rest is counted.
        let len = |s: &Self, pos: usize| {
            let count = || s.ab.both(pos).iter().filter(|&&x| x as usize != r).count();
            let fast = if s.kind[pos] == -2 && refcols.binary_search_by_key(&pos, |o| o.0).is_err() {
                s.ab.na(pos) + s.ab.nb(pos)
            } else {
                count()
            };
            debug_assert_eq!(fast, count());
            #[cfg(feature = "factor-shadow")]
            assert_eq!(fast, count(), "walk column length");
            fast
        };
        // The reference's work: each jump scans the column (2 per entry and 5) on the factor's scale; the final shift
        // moves the slots in between.
        let mut jumps = 0.0;
        while cur != hole {
            let pos = self.at[cur];
            let plen = len(self, pos);
            if plen >= cur - hole {
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
            jumps += 2.0 * plen as f64 + 5.0;
            let mut nxt = hole;
            for &x in self.ab.a(pos) {
                let x = x as usize;
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
    fn eliminate(&self, r: usize, w: &mut [R], seen: &mut [bool], scratch: &mut btran_ws::Update<R>) -> usize {
        let m = self.m;
        let d0 = self.diag[r];
        // Every entry of a row belongs to a tail column (audited), so the rows need no `kind` filter. The order inside a row is irrelevant here:
        // the entries of a row are distinct columns, and the heap is keyed by the unique `(slot, column)`.
        let urows = &self.urows;
        let btran_ws::Update { touched, heap, eta, .. } = scratch;
        touched.clear();
        heap.clear();
        eta.clear();
        let inv = R::ONE / d0;
        let (cols, vals) = urows.row(r);
        for (&k, &v) in cols.iter().zip(vals) {
            let k = k as usize;
            w[k] = v * inv;
            touched.push(k);
            seen[k] = true;
            heap.push(std::cmp::Reverse((self.slot[k], k)));
        }
        // The reference's solve (`0x131760`) books 2 per entry of the row-wise copy of every accepted position, `r` first.
        // Count only the rows reached by the solve.
        let rlen = |y: usize| urows.len(y);
        // Only positions in the tail (slot >= i21) are visited in the heap phase of the solve that books; the earlier ones are
        // handled in the list phase.
        let mut visited = if self.slot[r] >= self.i21 { rlen(r) } else { 0 };
        // Dense branch (`F+0x108 > 5`): `3` for the one input entry, `2 len` for `r` (a tail column is handled by the slot walk
        // like the rest), `2 (m - minslot)` for the slot walk over every position from the lowest slot reached, and `2 len` for
        // every accepted position of the walk.
        let mut dense = 3.0 + 2.0 * rlen(r) as f64;
        let mut minslot = if self.kind[r] == -2 { self.slot[r] } else { m };
        for &k in cols {
            minslot = minslot.min(self.slot[k as usize]);
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
            let (lcols, lvals) = urows.row(k);
            for (&l, &v) in lcols.iter().zip(lvals) {
                let l = l as usize;
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
            .set(if urows.is_empty(r) { 0.0 } else { dense });
        for &i in touched.iter() {
            w[i] = R::ZERO;
            seen[i] = false;
        }
        visited
    }

    /// Replace a numeric U column. Its entries left the row-wise copy and the new ones join it at the places `update` handles them (the
    /// rows of the old column at `0xf1e7a`, those of the new one at the end). Returns the replaced column so its buffer can be recycled.
    fn replace_u_column(&mut self, c: usize, col: Vec<(usize, R)>) -> Vec<(usize, R)> {
        std::mem::replace(&mut self.uc[c], col)
    }

    /// Replace the column at position `r` by the entering column of the last `ftran`.
    pub fn update(&mut self, r: usize) {
        let m = self.m;
        let mut scratch = self.btran_ws.update.take();
        if crate::dump::on("fstruct") {
            let n = crate::dump::FACTOR_STRUCT_CALLS.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
            crate::dump::emit(
                "fstruct",
                serde_json::json!({"r": r, "rowl": (0..m).filter(|&y| !self.urows.is_empty(y)).map(|y| (y.to_string(), self.urows.cols(y).to_vec())).collect::<std::collections::BTreeMap<_, _>>(), "mark": self.mark, "n": n, "m": m, "i3": self.nuc.len(), "i21": self.i21, "kind": self.kind, "nuc": self.nuc, "seqA": self.slot, "fill": self.fill, "diag": jv(&self.diag), "uvals": (0..m).filter(|&k| self.kind[k] == -2).map(|k| (k.to_string(), (self.uc[k].iter().filter(|e| self.kind[e.0] == -2).map(|e| e.0).collect::<Vec<_>>(), self.uc[k].iter().filter(|e| self.kind[e.0] == -2).map(|e| e.1.to_f64()).collect::<Vec<_>>()))).collect::<std::collections::BTreeMap<_, _>>(), "cols": (0..m).filter(|&k| self.kind[k] == -2).map(|k| (k.to_string(), [self.ab.a(k), self.ab.b(k)])).collect::<std::collections::BTreeMap<_, _>>()}),
            );
        }
        self.before_last = (self.upd, self.sum_alpha);
        self.upd += 1;
        self.upd_work.clear();
        if self.upd_work.capacity() == 0 {
            // The caller takes the list each call; one allocation instead of the doubling from empty.
            self.upd_work.reserve(32);
        }
        // The eta is stored first, at the current pool fill; its value array sits on a 16-byte boundary when the fill is even.
        let peel = self.fill % 2;
        // The old column's entries leave the row-wise copy after the new column's entries have been added.
        scratch.old_rows.clear();
        if self.kind[r] == -2 {
            scratch.old_rows.extend(self.uc[r].iter().map(|e| e.0));
        }
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
        let solve_len = self.eliminate(r, &mut values, &mut seen, &mut scratch);
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
        let mut sp = std::mem::take(&mut scratch.sp);
        sp.clear();
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
        let single_rows = &mut scratch.single_rows;
        single_rows.clear();
        if old == -2 {
            // The old column leaves the rows it appears in (`0xf1e7a`): the entries before it in each row list are passed.
            for &y in self.ab.both(r) {
                let y = y as usize;
                if let Some(p) = self.urows.position(y, r) {
                    let len = self.urows.len(y);
                    if len > 1 && p > 0 {
                        self.upd_work.push((0xf1f82, p as f64, crate::work::CTX));
                    }
                    if len == 1 {
                        single_rows.push(y);
                    }
                    self.urows.swap_remove(y, p);
                }
            }
        }
        // A new column without entries besides its diagonal is not a tail column: it stays among the non-tail columns (at the
        // last non-tail slot) and keeps its place in the nucleus lists that refer to it.
        let bare = k == 0 && old != -2;
        let demoted_front = &mut scratch.demoted_front;
        demoted_front.clear();
        let mut convslot: Option<usize> = None;
        let mut conv: Option<usize> = None;
        // The reference demotes a bare tail column only when the new column depends on the tail; with no tail entries in the spike
        // the tail simply grows by one slot.
        if old != -2
            && self.slot[r] + 1 != self.i21
            && sp.iter().any(|&y| y != r && self.kind[y] == -2)
        {
            conv = self.empty_tail.max_by_slot(&self.slot);
            if let Some(x) = conv {
                self.empty_tail.remove(&x);
                // Invariant behind the row lengths: every entry of `urows` belongs to a tail column (`kind == -2`). Only
                // columns of `empty_tail` leave the tail (here), so they must already have left the row-wise copy.
                debug_assert!(self.uc[x].is_empty(), "demoted empty_tail column {x} still has U entries");
                #[cfg(feature = "factor-shadow")]
                assert!(self.uc[x].is_empty(), "demoted empty_tail column {x} still has U entries");
                convslot = Some(self.slot[x]);
                // Only a column some tail column refers to joins the nucleus list; the others become singletons.
                let refs = &mut scratch.conversion_refs;
                refs.clear();
                refs.extend(self.urows.cols(x).iter().map(|&k| k as usize).filter(|&k| self.kind[k] == -2 && self.ab.a_pos(k, x).is_some()));
                refs.sort_unstable();
                refs.dedup();
                let held = !refs.is_empty();
                if held {
                    self.kind[x] = self.nuc.len() as i64;
                    self.nuc.push(x);
                } else {
                    self.kind[x] = -1;
                }
                for &k in refs.iter() {
                    let pa = if self.kind[k] == -2 { self.ab.a_pos(k, x) } else { None };
                    if let Some(p) = pa {
                        // The last entry of A fills the hole, and the demoted column takes the slot A gives up (the front of B).
                        self.ab.demote(k, p);
                        if p > 0 {
                            self.upd_work.push((0xf3eac, p as f64, crate::work::CTX));
                        }
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
            let held = self
                .urows
                .cols(r)
                .iter()
                .any(|&x| x as usize != r && self.kind[x as usize] == -2 && self.ab.b_pos(x as usize, r).is_some());
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
            let mut old_b = std::mem::take(&mut scratch.old_b);
            old_b.clear();
            old_b.extend_from_slice(self.ab.b(r));
            self.ab.clear(r);
            let old_col = self.replace_u_column(r, Vec::new());
            if old_col.capacity() > scratch.spare_col.capacity() {
                scratch.spare_col = old_col;
            }
            self.kind[r] = -1;
            for &x in &old_b {
                let x = x as usize;
                if self.kind[x] >= 0
                    && !self
                        .urows
                        .cols(x)
                        .iter()
                        .any(|&k| self.kind[k as usize] == -2 && self.ab.b_pos(k as usize, x).is_some())
                {
                    self.drop_nuc(x);
                    self.kind[x] = -1;
                }
            }
            scratch.old_b = old_b;
        }
        for &x in &sp {
            if self.kind[x] == -1 {
                self.kind[x] = self.nuc.len() as i64;
                self.nuc.push(x);
            }
        }
        // Built in recycled buffers (the lists replaced at the end of this call), with the same content and order as a fresh collect.
        // The lists carry the coefficients of the new U column's entries (`values` holds the spike at the positions of `sp`).
        let mut av = std::mem::take(&mut scratch.spare_a);
        av.clear();
        av.reserve(sp.len());
        av.extend(sp.iter().copied().filter(|&x| self.kind[x] == -2).map(|x| (x as u32, values[x])));
        let mut bv = std::mem::take(&mut scratch.spare_b);
        bv.clear();
        bv.reserve(sp.len());
        bv.extend(sp.iter().copied().filter(|&x| self.kind[x] >= 0).map(|x| (x as u32, values[x])));
        let refcols = &mut scratch.refcols;
        refcols.clear();
        refcols.extend(self.urows.cols(r).iter().map(|&k| (k as usize, self.slot[k as usize])));
        refcols.sort_unstable();
        refcols.dedup();
        refcols.retain(|&(x, _)| x != r && self.kind[x] == -2 && (self.ab.a_pos(x, r).is_some() || self.ab.b_pos(x, r).is_some()));
        // Only a surviving row-r reference can raise the mark after the slot walk.
        let oldslot = self.slot[r];
        let top: i64 = av.iter().map(|&(x, _)| self.slot[x as usize] as i64).max().unwrap_or(-1);
        let mut moved = true;
        let mut tail_rebuild = false;
        if old == -2 {
            if top >= oldslot as i64 {
                self.walk(r, oldslot, top as usize, &refcols[..]);
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
                self.walk(r, cs, top as usize, &refcols[..]);
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
                self.walk(r, bt, top as usize, &refcols[..]);
            }
        }
        let _ = &mut tail_rebuild;
        // A row whose only entry was the old column leaves the singleton list when its position is within the reference's limit.
        // GUROBI: the reference removes the row-single nucleus column from the old tail column's lists only if `slot[y] <= [rsp+0x168]` (`0xf1ea5`);
        // rose removes it unconditionally above; no witness of the gated-off case (1693 corpus drops, `specs/lu-factor.md`, "latent gate").
        for &y in single_rows.iter() {
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
        let mut eliminated = false;
        for &(k, old_ref_slot) in refcols.iter() {
            if k == r || self.kind[k] != -2 {
                continue;
            }
            let (pos_a, pos_b) = (self.ab.a_pos(k, r), self.ab.b_pos(k, r));
            let (in_a, in_b) = (pos_a.is_some(), pos_b.is_some());
            debug_assert!(!(in_a && in_b), "row {r} in both lists of column {k}");
            if bare {
                continue;
            }
            if in_a || in_b {
                // The reference scans the lists from their start for `r`; the booking is the number of entries passed.
                let pa = pos_a.unwrap_or(0);
                let pb = pos_b
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
                    // `r` moves from B to the end of A (the front slot of B, which A takes), B's front entry filling its hole.
                    if let Some(p) = pos_b {
                        self.ab.b_to_a(k, p);
                    }
                } else {
                    if let Some(p) = pos_a {
                        // A and B share one array in the reference (A first): A's last entry fills the hole and B's last entry
                        // takes the slot A gives up, which is the front of B.
                        self.ab.drop_a(k, p);
                    }
                    if let Some(p) = pos_b {
                        self.ab.drop_b(k, p);
                    }
                    // Row `r` of the row-wise copy is emptied below (`!keep`), which takes `k`'s entry with it.
                    self.uc[k].retain(|e| e.0 != r);
                    eliminated = true;
                    // A column left without entries raises the mark to its slot (`0xf29d5`).
                    if self.ab.is_empty(k) {
                        self.empty_tail.insert(k);
                        self.mark = self.mark.max(old_ref_slot);
                    }
                }
            }
        }
        if !keep {
            self.urows.clear_row(r);
        }
        let mut newdiag = values[r];
        if eliminated {
            let mut ind = Vec::with_capacity(1 + scratch.eta.len());
            ind.push(r);
            let mut val = Vec::with_capacity(1 + scratch.eta.len());
            val.push(R::ONE / d0);
            for &(k, c) in &scratch.eta {
                ind.push(k);
                val.push(c);
            }
            let xs = &mut scratch.xs;
            xs.clear();
            xs.extend(ind.iter().map(|&j| values[j]));
            // The diagonal is summed from the scratch eta, whose alignment is fixed per factor object.
            newdiag = dot_diag(xs, &val, self.scratch_peel);
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
        if k > scratch.old_rows.len() {
            self.fill += k;
        }
        if k > 0 {
            // A row that already holds an entry of the old column takes the new entry in place.
            // The rows are visited in the spike's order with the leaving column replaced by the last entry (`sp`).
            for &row in sp.iter().filter(|&&y| y != r && !scratch.old_rows.contains(&y)) {
                // Every entry of `urows` belongs to a tail column (audited) and `r`'s own entries sit only in `old_rows`, skipped above.
                let count = self.urows.len(row);
                let counted = || self.urows.cols(row).iter().filter(|&&c| c as usize != r && self.kind[c as usize] == -2).count();
                debug_assert_eq!(count, counted());
                #[cfg(feature = "factor-shadow")]
                assert_eq!(count, counted(), "update count scan");
                if count >= self.row_cap[row] {
                    // The growth is capped by the updates left before the refactorization (`0xf2cd1`: `F+0x188 - F+0x18c`).
                    let left = (self
                        .upd_limit
                        .unwrap_or(((m as f64).powf(0.7) as usize).max(20))
                        + 1)
                    .saturating_sub(self.upd);
                    let extra = (count / 2).max(4).min(left);
                    if crate::dump::on("fillparts") {
                        fp_rows.push((row, count, self.row_cap[row], extra));
                    }
                    self.row_par[row] = (self.fill % 2) as u8;
                    self.fill += count + extra;
                    self.row_cap[row] = count + extra;
                }
            }
        }
        if crate::dump::on("fillparts") {
            let n = crate::dump::FACTOR_FILL_CALLS.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
            crate::dump::emit(
                "fillparts",
                serde_json::json!({"n": n, "k": k, "old": scratch.old_rows.len(), "col": if k > scratch.old_rows.len() { k } else { 0 }, "fill_after_eta": fill_in, "rows": fp_rows}),
            );
        }
        if self.bt_occupied && !newdiag.to_f64().is_normal() {
            self.bt_occupied = false;
        }
        self.bad_diag |= !(newdiag.abs() > R::ZERO);
        self.diag[r] = newdiag;
        let mut col = std::mem::take(&mut scratch.spare_col);
        col.clear();
        col.reserve(self.spike.len());
        col.extend(self.spike.iter().copied().filter(|e| e.0 != r));
        let old_col = self.replace_u_column(r, col);
        if old_col.capacity() > scratch.spare_col.capacity() {
            scratch.spare_col = old_col;
        }
        bv.reverse();
        self.ab.set(r, &av, &bv);
        // The new column joins the end of each of its rows, with the coefficient bits of its A and B entries.
        for (&y, &v) in self.ab.both(r).iter().zip(self.ab.both_val(r)) {
            self.urows.push(y as usize, r, v);
        }
        scratch.spare_a = av;
        scratch.spare_b = bv;
        if self.kind[r] == -2 && self.ab.is_empty(r) {
            self.empty_tail.insert(r);
        } else {
            self.empty_tail.remove(&r);
        }
        self.spike_used = true;
        #[cfg(any(test, debug_assertions, feature = "factor-shadow"))]
        {
            self.ab.check(&self.uc, &self.kind);
            self.urows.check(&self.uc, &self.kind);
        }
        scratch.sp = sp;
        scratch.clear();
        self.btran_ws.update.replace(scratch);
    }
}

#[cfg(test)]
/// Delete `r` from a column's nucleus list the way the reference does: the entry at the front of this (reversed) list fills the hole.
fn remove_front_fill(v: &mut Vec<usize>, r: usize) {
    if let Some(p) = v.iter().position(|&y| y == r) {
        if p != 0 {
            v[p] = v[0];
        }
        v.remove(0);
    }
}

#[cfg(test)]
fn remove_last_fill(v: &mut Vec<usize>, r: usize) {
    if let Some(p) = v.iter().position(|&y| y == r) {
        v.swap_remove(p);
    }
}
