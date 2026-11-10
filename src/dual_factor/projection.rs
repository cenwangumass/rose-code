//! Representation-neutral, bit-exact shadow boundary. Pool offsets and scratch
//! allocation addresses are deliberately not part of this logical projection.
use crate::real::Real;
use std::any::Any;

pub(crate) fn bits<R: Real>(value: R) -> String {
    let value = &value as &dyn Any;
    if let Some(x) = value.downcast_ref::<f64>() {
        format!("0x{:016x}", x.to_bits())
    } else if let Some(x) = value.downcast_ref::<crate::f80::F80>() {
        let (negative, exponent, significand) = x.factor_bits();
        format!("0x{:01x}:{:08x}:{:016x}", negative as u8, exponent as u32, significand)
    } else {
        panic!("factor-shadow: unsupported Real type")
    }
}

pub(crate) fn numbers<R: Real>(values: &[R]) -> Vec<String> {
    values.iter().map(|&x| bits(x)).collect()
}

// The hot boundary must not allocate a hex string for every entry on success.
pub(crate) fn same_bits<R: Real>(a: R, b: R) -> bool {
    let a = &a as &dyn Any;
    let b = &b as &dyn Any;
    if let Some(a) = a.downcast_ref::<f64>() {
        a.to_bits() == b.downcast_ref::<f64>().unwrap().to_bits()
    } else if let Some(a) = a.downcast_ref::<crate::f80::F80>() {
        a.factor_bits() == b.downcast_ref::<crate::f80::F80>().unwrap().factor_bits()
    } else {
        panic!("factor-shadow: unsupported Real type")
    }
}

pub(crate) fn check_numbers<R: Real>(method: &str, call: usize, field: &str, new: &[R], old: &[R]) {
    check(method, call, &format!("{field}.len"), &new.len(), &old.len());
    for (i, (&a, &b)) in new.iter().zip(old).enumerate() {
        if !same_bits(a, b) {
            check(method, call, &format!("{field}[{i}]"), &bits(a), &bits(b));
        }
    }
}

pub(crate) fn sparse<R: Real>(values: &[(usize, R)]) -> Vec<(usize, String)> {
    values.iter().map(|&(i, x)| (i, bits(x))).collect()
}

pub(crate) fn columns<R: Real>(values: &[Vec<(usize, R)>]) -> Vec<Vec<(usize, String)>> {
    values.iter().map(|v| sparse(v)).collect()
}

/// Returns the first difference, with a path that includes field and vector index.
pub(crate) trait Exact {
    fn difference(&self, old: &Self, path: &str) -> Option<String>;
}

macro_rules! atom {
    ($($t:ty),*) => { $(impl Exact for $t {
        fn difference(&self, old: &Self, path: &str) -> Option<String> {
            (self != old).then(|| format!("{path}: new={self:?}, old={old:?}"))
        }
    })* };
}
atom!(bool, String);

macro_rules! integer {
    ($($t:ty),*) => { $(impl Exact for $t {
        fn difference(&self, old: &Self, path: &str) -> Option<String> {
            (self != old).then(|| format!("{path}: new={self:#x}, old={old:#x}"))
        }
    })* };
}
integer!(usize, u32, u8, i64);

impl Exact for () {
    fn difference(&self, _: &Self, _: &str) -> Option<String> { None }
}

impl<T: Exact> Exact for Vec<T> {
    fn difference(&self, old: &Self, path: &str) -> Option<String> {
        if self.len() != old.len() {
            return Some(format!("{path}.len: new={}, old={}", self.len(), old.len()));
        }
        for (i, (a, b)) in self.iter().zip(old).enumerate() {
            if let Some(d) = a.difference(b, &format!("{path}[{i}]")) { return Some(d); }
        }
        None
    }
}

impl<T: Exact> Exact for Option<T> {
    fn difference(&self, old: &Self, path: &str) -> Option<String> {
        match (self, old) {
            (Some(a), Some(b)) => a.difference(b, path),
            (None, None) => None,
            _ => Some(format!("{path}.is_some: new={}, old={}", self.is_some(), old.is_some())),
        }
    }
}

impl<A: Exact, B: Exact> Exact for (A, B) {
    fn difference(&self, old: &Self, path: &str) -> Option<String> {
        self.0.difference(&old.0, &format!("{path}.0"))
            .or_else(|| self.1.difference(&old.1, &format!("{path}.1")))
    }
}

impl<A: Exact, B: Exact, C: Exact> Exact for (A, B, C) {
    fn difference(&self, old: &Self, path: &str) -> Option<String> {
        self.0.difference(&old.0, &format!("{path}.0"))
            .or_else(|| self.1.difference(&old.1, &format!("{path}.1")))
            .or_else(|| self.2.difference(&old.2, &format!("{path}.2")))
    }
}

pub(crate) fn check<T: Exact>(method: &str, call: usize, field: &str, new: &T, old: &T) {
    if let Some(d) = new.difference(old, field) {
        panic!("factor-shadow {method} call {call}: {d}");
    }
}

macro_rules! projection {
    ($name:ident { $($field:ident: $ty:ty),* $(,)? }) => {
        #[derive(Debug)]
        pub(crate) struct $name { $(pub(crate) $field: $ty),* }
        impl Exact for $name {
            fn difference(&self, old: &Self, path: &str) -> Option<String> {
                $(if let Some(d) = self.$field.difference(&old.$field, &format!("{path}.{}", stringify!($field))) {
                    return Some(d);
                })*
                None
            }
        }
    };
}

projection!(EtaState {
    target: usize,
    indices: Vec<usize>,
    values: Vec<String>,
    length: usize,
    alignment_peel: usize,
});

/// Dense snapshots contain every value the old consumer could read, including
/// signed zeros outside sparse supports. Reset bookkeeping is validated locally,
/// not compared: a dense clear and a covered sparse clear have the same meaning.
pub(crate) struct F2State<R> {
    pub(crate) column_count: Option<usize>,
    pub(crate) column_input: Vec<usize>,
    pub(crate) column: Vec<R>,
    pub(crate) tau_count: Option<usize>,
    pub(crate) tau: Vec<R>,
    pub(crate) tau_raw: Vec<bool>,
}

projection!(Canonical {
    m: usize,
    i21: usize,
    mark: usize,
    kind: Vec<i64>,
    nucleus: Vec<usize>,
    slot: Vec<usize>,
    at: Vec<usize>,
    diagonal: Vec<String>,
    lower_columns: Vec<Vec<(usize, String)>>,
    upper_columns: Vec<Vec<(usize, String)>>,
    upper_rows: Vec<Vec<(usize, String)>>,
    tail_destinations: Vec<Vec<usize>>,
    nucleus_destinations: Vec<Vec<usize>>,
    row_columns: Vec<Vec<usize>>,
    row_capacity: Vec<usize>,
    row_alignment: Vec<u8>,
    pool_fill: usize,
    etas: Vec<EtaState>,
    first_eta: Vec<usize>,
    updates: usize,
    alpha_sum: String,
    eta_sum: String,
    previous_window: (usize, String),
    last_alpha: usize,
    factor_size: String,
    tail_density: String,
    nucleus_density: String,
    lower_density_threshold: usize,
    cached_pivot: Option<usize>,
    spike: Vec<(usize, String)>,
    spike_consumed: bool,
    permutation: Vec<usize>,
    lower_transpose: Vec<Vec<(usize, String)>>,
    lower_levels: Option<(Vec<usize>, usize)>,
    empty_tail_positions: Vec<usize>,
    scratch_alignment: usize,
});

projection!(Observables {
    work: String,
    work_parts: (String, String),
    keep_work: String,
    keep_dense: bool,
    tau_lower_work: String,
    tau_upper_work: String,
    tau_count: usize,
    btran_dense: bool,
    last_dense: bool,
    last_ordered: bool,
    update_work: Vec<(u32, String, String)>,
    update_limit: Option<usize>,
    refactor_due: bool,
    upper_nnz: usize,
    has_levels: bool,
    updates: usize,
    alpha_sum: String,
    eta_sum: String,
    previous_window: (usize, String),
    last_alpha: usize,
    factor_size: String,
    density: String,
    dense_estimate: String,
    tail_density: String,
    nucleus_density: String,
    lower_density_threshold: usize,
    cached_pivot: Option<usize>,
    spike_consumed: bool,
    btran_occupied: bool,
    transpose_levels_dense: bool,
    forward_levels_dense: bool,
    phase_one_steps: usize,
    singleton_steps: usize,
    lower_nnz: usize,
    lower_tail_nnz: usize,
});

// Only the initial representation adapters know the old field names. A later
// backend can replace its adapter without changing either projection's schema.
macro_rules! impl_factor_projection {
    ($factor:ident) => {
        impl<R: crate::real::Real> $factor<R> {
            /// Validate caches against their source, rather than requiring equal
            /// allocation/lazy initialization histories on the two backends.
            pub(crate) fn check_retained(&self, method: &str, call: usize) {
                use crate::dual_factor::projection::same_bits;
                let require = |ok: bool, field: &str| {
                    assert!(ok, "factor-shadow {method} call {call}: {}.{field} invariant", stringify!($factor));
                };
                if let Some(inverse) = self.row_of.get() {
                    require(inverse.len() == self.m, "row_of.len");
                    for (r, &p) in self.perm.iter().enumerate() {
                        require(inverse[p] == r, "row_of");
                    }
                }
                if let Some((levels, maximum)) = self.lt_levels.get() {
                    require(levels.len() == self.m, "lt_levels.len");
                    for p in 0..self.m {
                        let expected = self.lcolumn(p).iter().map(|&(to, _)| levels[to] + 1).max().unwrap_or(0);
                        require(levels[p] == expected, "lt_levels");
                    }
                    require(*maximum == levels.iter().copied().max().unwrap_or(0), "lt_levels.maximum");
                }
                let f2 = self.f2.borrow();
                // A dirty-block list means all unlisted snapshot entries will
                // survive the next clear, and must already be positive zero.
                if let Some(blocks) = &f2.dirty_blocks {
                    let mut dirty = vec![false; self.m.div_ceil(64)];
                    for &block in blocks {
                        require(block < dirty.len(), "f2.dirty_blocks.range");
                        dirty[block] = true;
                    }
                    for (i, &value) in f2.x_col.iter().enumerate() {
                        require(dirty[i / 64] || same_bits(value, R::ZERO), "f2.dirty_blocks.coverage");
                    }
                }
                require(self.update_values.iter().all(|&v| same_bits(v, R::ZERO)), "update_values.clean");
                require(self.update_seen.iter().all(|&v| !v), "update_seen.clean");
                if let Some(ws) = &self.bt_work {
                    require(ws.x.is_empty() || ws.x.len() == self.m, "bt_work.x.len");
                    require(ws.x.iter().all(|&v| same_bits(v, R::ZERO)), "bt_work.x.clean");
                    // Tail values behind clear occupancy bits are logically +0;
                    // entry() initializes them before a subsequent read.
                    require(ws.tail.occupied.iter().all(|&v| v == 0), "bt_work.tail.clean");
                    require(ws.pending.is_empty(), "bt_work.pending.clean");
                    require(ws.marks.is_empty() || ws.marks.len() == self.m, "bt_work.marks.len");
                    require(ws.marks.iter().all(|&v| v <= ws.generation), "bt_work.marks.epoch");
                    require(ws.eta_prev.len() <= self.etas.len(), "bt_work.eta_index.len");
                    let mut heads = vec![usize::MAX; self.m];
                    for (e, &previous) in ws.eta_prev.iter().enumerate() {
                        let row = self.etas[e].idx;
                        require(previous == heads[row], "bt_work.eta_index.previous");
                        heads[row] = e;
                    }
                    require(ws.eta_heads.is_empty() || ws.eta_heads == heads, "bt_work.eta_index.heads");
                }
                if let Some(ws) = self.level_cache.0.borrow().as_ref() {
                    require(ws.term_ends.len() <= self.etas.len(), "level_cache.eta_index.len");
                    // prepare() may stop partway through an incompatible eta.
                    // That cache is permanently bypassed, including its scratch.
                    if ws.incompatible {
                        let mut rebuilt = LevelWorkspace::<R>::default();
                        let nlevels = self.levels.as_ref().map_or(0, |(_, maximum)| maximum + 1);
                        require(!rebuilt.prepare(self, nlevels), "level_cache.incompatible");
                    } else {
                        require(ws.values.len() == self.m && ws.state.len() == self.m, "level_cache.values.len");
                        require(ws.buckets.len() == self.levels.as_ref().map_or(0, |(_, maximum)| maximum + 1),
                            "level_cache.buckets.len");
                        let mut uses = vec![Vec::new(); self.m];
                        let mut end = 0;
                        for (e, &term_end) in ws.term_ends.iter().enumerate() {
                            let eta = &self.etas[e];
                            end += eta.ind.len();
                            require(term_end == end, "level_cache.term_ends");
                            require(eta.ind.first() == Some(&eta.idx)
                                && eta.val.iter().all(|v| v.to_f64().is_finite())
                                && end as u64 <= (1u64 << 50), "level_cache.compatible");
                            for &row in &eta.ind {
                                if uses[row].is_empty() {
                                    require(self.first_eta[row] == e, "level_cache.first_eta");
                                }
                                if uses[row].last() != Some(&e) { uses[row].push(e); }
                            }
                        }
                        require((ws.term_ends.is_empty() && ws.uses.is_empty()) || ws.uses == uses, "level_cache.uses");
                        require(ws.pending.len() == ws.term_ends.len().div_ceil(64), "level_cache.pending.len");
                        require(ws.pending.iter().all(|&v| v == 0), "level_cache.pending.clean");
                        require(ws.buckets.iter().all(Vec::is_empty), "level_cache.buckets.clean");
                        // Values/state of touched rows are deliberately retained;
                        // prepare clears precisely these rows before the next solve.
                        let mut touched = vec![false; self.m];
                        for &row in &ws.touched {
                            require(row < self.m, "level_cache.touched.range");
                            touched[row] = true;
                        }
                        require(ws.values.len() == ws.state.len(), "level_cache.state.len");
                        for (i, (&v, &state)) in ws.values.iter().zip(&ws.state).enumerate() {
                            require(state <= 2 && (touched[i] || (state == 0 && same_bits(v, R::ZERO))),
                                "level_cache.reset_coverage");
                        }
                    }
                }
                let row = self.tau_row.borrow();
                // prepare bumps the epoch (clearing stamps on wrap), and terms
                // are cleared before append_terms. Neither retained numeric
                // buffer is logically live between calls.
                require(row.values.len() == row.stamps.len(), "tau_row.stamps.len");
                require(row.stamps.iter().all(|&v| v <= row.epoch), "tau_row.stamps.epoch");
                require(row.epoch != 0 || row.stamps.iter().all(|&v| v == 0), "tau_row.epoch.zero");
            }

            pub(crate) fn f2_projection(&self) -> crate::dual_factor::projection::F2State<R> {
                let f = self.f2.borrow();
                crate::dual_factor::projection::F2State {
                    column_count: f.col_n, column_input: f.col_in.clone(), column: f.x_col.clone(),
                    tau_count: f.tau_n, tau: f.x_tau.clone(), tau_raw: f.tau_raw.clone(),
                }
            }

            pub(crate) fn canonical(&self) -> crate::dual_factor::projection::Canonical {
                use crate::dual_factor::projection::{bits, columns, numbers, sparse, Canonical, EtaState};
                Canonical {
                    m: self.m, i21: self.i21, mark: self.mark,
                    kind: self.kind.clone(), nucleus: self.nuc.clone(),
                    slot: self.slot.clone(), at: self.at.clone(), diagonal: numbers(&self.diag),
                    lower_columns: self.lcolumns().map(sparse).collect(),
                    upper_columns: columns(&self.uc), upper_rows: columns(&self.ut),
                    tail_destinations: self.a.clone(), nucleus_destinations: self.b.clone(),
                    row_columns: self.rowl.clone(), row_capacity: self.row_cap.clone(), row_alignment: self.row_par.clone(),
                    pool_fill: self.fill,
                    etas: self.etas.iter().map(|e| EtaState {
                        target: e.idx, indices: e.ind.clone(), values: numbers(&e.val),
                        length: e.ind.len(), alignment_peel: e.peel,
                    }).collect(),
                    first_eta: self.first_eta.clone(), updates: self.upd,
                    alpha_sum: bits(self.sum_alpha), eta_sum: bits(self.sum_eta),
                    previous_window: (self.before_last.0, bits(self.before_last.1)),
                    last_alpha: self.last_alpha, factor_size: bits(self.znz),
                    tail_density: bits(self.h33), nucleus_density: bits(self.h34), lower_density_threshold: self.dthr,
                    cached_pivot: self.row_cached, spike: sparse(&self.spike), spike_consumed: self.spike_used,
                    permutation: self.perm.clone(), lower_transpose: columns(&self.lt), lower_levels: self.levels.clone(),
                    empty_tail_positions: self.empty_tail.iter().copied().collect(), scratch_alignment: self.scratch_peel,
                }
            }

            pub(crate) fn observables(&self) -> crate::dual_factor::projection::Observables {
                use crate::dual_factor::projection::{bits, Observables};
                Observables {
                    work: bits(self.work), work_parts: (bits(self.work_parts.0), bits(self.work_parts.1)),
                    keep_work: bits(self.keep_work), keep_dense: self.keep_dense,
                    tau_lower_work: bits(self.tau_l.get()), tau_upper_work: bits(self.tau_u.get()), tau_count: self.tau_fin.get(),
                    btran_dense: self.bt_dense, last_dense: self.last_dense, last_ordered: self.last_ordered,
                    update_work: self.upd_work.iter().map(|&(s, e, w)| (s, bits(e), bits(w))).collect(),
                    update_limit: self.upd_limit, refactor_due: self.refactor_due(), upper_nnz: self.u_nnz(),
                    has_levels: self.has_levels(), updates: self.upd,
                    alpha_sum: bits(self.sum_alpha), eta_sum: bits(self.sum_eta),
                    previous_window: (self.before_last.0, bits(self.before_last.1)),
                    last_alpha: self.last_alpha, factor_size: bits(self.znz), density: bits(self.dens),
                    dense_estimate: bits(self.dense_est.get()), tail_density: bits(self.h33), nucleus_density: bits(self.h34),
                    lower_density_threshold: self.dthr, cached_pivot: self.row_cached, spike_consumed: self.spike_used,
                    btran_occupied: self.bt_occupied, transpose_levels_dense: self.lev_dense,
                    forward_levels_dense: self.lev_dense_ft, phase_one_steps: self.n1, singleton_steps: self.loff,
                    lower_nnz: self.nnz_l, lower_tail_nnz: self.nnz_l_tail,
                }
            }
        }
    };
}
pub(crate) use impl_factor_projection;

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn retained_numbers_are_bit_exact() {
        assert!(!same_bits(0.0f64, -0.0));
        let nan = f64::from_bits(0x7ff8_0000_0000_0123);
        assert!(same_bits(nan, nan));
        assert!(!same_bits(nan, f64::from_bits(0x7ff8_0000_0000_0124)));
        assert!(!same_bits(crate::f80::F80::from_f64(0.0), crate::f80::F80::from_f64(-0.0)));
    }

    #[test]
    #[should_panic(expected = "f2.x_col[1]")]
    fn retained_snapshot_reports_the_index() {
        check_numbers("ftran", 2, "f2.x_col", &[1.0f64, -0.0], &[1.0, 0.0]);
    }
}
