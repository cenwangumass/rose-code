//! Projection adapter of `DualFactor`: the field-specific part of `canonical()`, `observables()`, `f2_projection()` and
//! `check_retained()` (a copy of `impl_factor_projection!`, which `RefStruct` keeps using). Representation changes of
//! `DualFactor` edit only this file; the projected schema stays the oracle's.

use super::*;
use crate::real::Real;

impl<R: Real> DualFactor<R> {
    /// L transpose rows as `(source, multiplier)` lists.
    fn lt_columns(&self) -> Vec<Vec<(usize, R)>> {
        (0..self.m).map(|v| self.lt.iter(v).collect()).collect()
    }

    /// Validate caches against their source, rather than requiring equal
    /// allocation/lazy initialization histories on the two backends.
    pub(crate) fn check_retained(&self, method: &str, call: usize) {
        use crate::dual_factor::projection::same_bits;
        let require = |ok: bool, field: &str| {
            assert!(ok, "factor-shadow {method} call {call}: {}.{field} invariant", "DualFactor");
        };
        // A/B lists: the audited sets, no duplicates, coefficients carrying the bits of the matching `uc` entries.
        self.ab.check(&self.uc, &self.kind);
        if let Some(inverse) = self.row_of.get() {
            require(inverse.len() == self.m, "row_of.len");
            for (r, &p) in self.perm.iter().enumerate() {
                require(inverse[p] == r, "row_of");
            }
        }
        if let Some((levels, maximum)) = self.lt_levels.get() {
            require(levels.len() == self.m, "lt_levels.len");
            for p in 0..self.m {
                let expected = self.lcolumn(p).0.iter().map(|&to| levels[to as usize] + 1).max().unwrap_or(0);
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
        self.urows.check(&self.uc, &self.kind);
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
        let (tail_destinations, nucleus_destinations) = self.ab.to_vecs();
        Canonical {
            m: self.m, i21: self.i21, mark: self.mark,
            kind: self.kind.clone(), nucleus: self.nuc.clone(),
            slot: self.slot.clone(), at: self.at.clone(), diagonal: numbers(&self.diag),
            lower_columns: columns(&self.lcolumns()),
            upper_columns: columns(&self.uc), upper_rows: columns(&self.urows.to_ut()),
            tail_destinations, nucleus_destinations,
            row_columns: self.urows.to_rowl(), row_capacity: self.row_cap.clone(), row_alignment: self.row_par.clone(),
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
            permutation: self.perm.clone(), lower_transpose: columns(&self.lt_columns()), lower_levels: self.levels.clone(),
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
