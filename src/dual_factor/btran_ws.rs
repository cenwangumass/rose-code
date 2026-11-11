//! Retained scratch for BTRAN, ordering/work estimates and column replacement.
use crate::real::Real;
use std::cell::RefCell;
use std::collections::BinaryHeap;

/// Generation-stamped markers; an epoch change logically clears every position.
#[derive(Default)]
pub(super) struct Marks {
    values: Vec<u32>,
    generation: u32,
}

impl Marks {
    pub fn prepare(&mut self, m: usize) {
        self.clear();
        self.values.resize(m, 0);
    }

    fn clear(&mut self) {
        if self.generation == u32::MAX {
            self.values.fill(0);
            self.generation = 0;
        }
        self.generation += 1;
    }

    /// Returns true only on the first visit in this epoch.
    pub fn insert(&mut self, i: usize) -> bool {
        let first = self.values[i] != self.generation;
        self.values[i] = self.generation;
        first
    }
}

#[derive(Default)]
pub(super) struct Work {
    pub mark: Marks,
    pub nm: Marks,
    pub tail: Vec<usize>,
    pub surv: Vec<usize>,
    pub heap: BinaryHeap<usize>,
}

impl Work {
    pub fn clear(&mut self) {
        self.mark.clear();
        self.nm.clear();
        self.tail.clear();
        self.surv.clear();
        self.heap.clear();
    }
}

#[derive(Default)]
pub(super) struct Order {
    pub seen: Marks,
    pub have: Marks,
    pub nm: Marks,
    pub heap: BinaryHeap<usize>,
    /// One bit per key for ordering long lists; all zero between calls.
    pub bits: Vec<u64>,
    pub inp_d: Vec<usize>,
    pub sing: Vec<usize>,
    pub tri: Vec<usize>,
    pub nuc: Vec<usize>,
    pub t: Vec<usize>,
    pub ns: Vec<usize>,
    pub order: Vec<usize>,
    pub out: Vec<usize>,
}

impl Order {
    pub fn clear(&mut self) {
        self.seen.clear();
        self.have.clear();
        self.nm.clear();
    }
}

#[derive(Default)]
pub(super) struct Dense<R> {
    pub x: Vec<R>,
    pub xs: Vec<R>,
}

#[derive(Default)]
pub(super) struct Update<R> {
    pub touched: Vec<usize>,
    pub heap: BinaryHeap<std::cmp::Reverse<(usize, usize)>>,
    pub eta: Vec<(usize, R)>,
    pub sp: Vec<usize>,
    pub old_rows: Vec<usize>,
    pub single_rows: Vec<usize>,
    pub demoted_front: Vec<usize>,
    pub refcols: Vec<(usize, usize)>,
    pub conversion_refs: Vec<usize>,
    pub xs: Vec<R>,
    /// Buffers recycled between updates: the replaced column lists and U column of the previous call (contents irrelevant, capacity kept).
    pub spare_a: Vec<(u32, R)>,
    pub spare_b: Vec<(u32, R)>,
    pub old_b: Vec<u32>,
    pub spare_col: Vec<(usize, R)>,
}

impl<R> Update<R> {
    pub fn clear(&mut self) {
        self.touched.clear();
        self.heap.clear();
        self.eta.clear();
        self.sp.clear();
        self.old_rows.clear();
        self.single_rows.clear();
        self.demoted_front.clear();
        self.refcols.clear();
        self.conversion_refs.clear();
        self.xs.clear();
        self.spare_a.clear();
        self.spare_b.clear();
        self.old_b.clear();
        self.spare_col.clear();
    }
}

#[derive(Default)]
pub(super) struct Workspace<R> {
    pub work: RefCell<Work>,
    pub order: RefCell<Order>,
    pub dense: RefCell<Dense<R>>,
    pub update: RefCell<Update<R>>,
}

impl<R: Real> Clone for Workspace<R> {
    fn clone(&self) -> Self { Self::default() }
}

#[cfg(test)]
mod tests {
    use super::Marks;

    #[test]
    fn stamps_reset_on_epoch_wrap_and_repeated_calls() {
        let mut marks = Marks::default();
        for _ in 0..8 {
            marks.prepare(129);
            assert!(marks.insert(64));
            assert!(!marks.insert(64));
            marks.clear();
        }
        marks.generation = u32::MAX;
        marks.values[64] = u32::MAX;
        marks.prepare(129);
        assert_eq!(marks.generation, 1);
        assert!(marks.insert(64));
        assert!(!marks.insert(64));
    }
}

#[cfg(feature = "factor-shadow")]
impl<R: Real> super::DualFactor<R> {
    pub(crate) fn check_btran_workspace(&self) {
        let same_bits = |a, b| super::projection::bits(a) == super::projection::bits(b);
        let clean = |marks: &Marks| marks.values.iter().all(|&v| v == 0 || v < marks.generation);
        let work = self.btran_ws.work.borrow();
        assert!(clean(&work.mark) && clean(&work.nm), "factor-shadow BTRAN WORK markers not clean");
        assert!(work.tail.is_empty() && work.surv.is_empty() && work.heap.is_empty(), "factor-shadow BTRAN WORK lists not clean");
        let order = self.btran_ws.order.borrow();
        assert!(clean(&order.seen) && clean(&order.have) && clean(&order.nm), "factor-shadow U-order markers not clean");
        let dense = self.btran_ws.dense.borrow();
        assert!(dense.x.iter().all(|&v| same_bits(v, R::ZERO)), "factor-shadow dense WORK values not clean");
        if let Some(ws) = &self.bt_work {
            assert!(ws.dense_tail.iter().chain(&ws.partial).all(|&v| same_bits(v, R::ZERO)),
                "factor-shadow BTRAN retained dense values not clean");
            assert!(ws.buckets.iter().all(Vec::is_empty), "factor-shadow BTRAN buckets not clean");
        }
        let update = self.btran_ws.update.borrow();
        assert!(update.touched.is_empty() && update.heap.is_empty() && update.eta.is_empty()
            && update.sp.is_empty() && update.old_rows.is_empty() && update.single_rows.is_empty()
            && update.demoted_front.is_empty() && update.refcols.is_empty() && update.conversion_refs.is_empty()
            && update.xs.is_empty(), "factor-shadow update scratch lists not clean");
    }
}
