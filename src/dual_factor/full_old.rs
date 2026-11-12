//! The `ftran_with` / `ftran_with_dense` of the stack before F3b (per-column length lookups, `FtranTouched::insert` on every scatter,
//! a fresh output `Vec`, a separate final pass and output scan), frozen as the oracle of the differential tests in `full_tests`.
use super::*;

impl<R: Real> DualFactor<R> {
    pub(super) fn ftran_with_old(&mut self, rows: &[usize], vals: &[R], force_dense: bool) -> Vec<(usize, R)> {
        let m = self.m;
        let (lout, count_l, e_l) = self.lstage_num(rows, vals, force_dense);
        let spike = match lout {
            LOut::List(list) => list,
            LOut::Dense(mut x) => {
                let mut inp = std::mem::take(&mut self.f2.get_mut().col_in);
                let mut sp = std::mem::take(&mut self.spike);
                let nz = Self::dense_lists(&x, &mut inp, &mut sp);
                let u_dense = rows.len() * self.dthr > m
                    || self.lev_dense_ft
                    || (!self.etas.is_empty() && Self::count_dense(count_l, self.dthr, m))
                    || !(self.nuc.len() + nz < self.i21 || self.dens * (nz as f64) < 2.0 * m as f64);
                if u_dense && !(force_dense && crate::dump::on("bsolve")) {
                    return self.ftran_with_dense_old(rows, x, inp, sp, nz, count_l, e_l);
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
        let mut touched = FtranTouched::new(m, u_dense, dirty.unwrap_or_default());
        let inp: Vec<usize> = spike.iter().map(|e| e.0).collect();
        for &p in &inp {
            touched.insert(p);
        }
        // The sparse solve runs the closure of the list's tail positions under the `a` lists in descending slot order, which is also the
        // tail part of the output order: `uorder_open` produces it once and `uorder_close` finishes the order after the tail loop.
        let opened = !(u_dense && !bs);
        let tail: Vec<usize> = if opened {
            let tail = self.uorder_open(&inp);
            for &c in &tail {
                touched.insert(c);
            }
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
            col_surv.push(c);
            sa += self.ab.na(c);
            sb += self.ab.nb(c);
            su += self.uc[c].len();
            for &(j, v) in &self.uc[c] {
                x[j] -= v * xc;
                touched.insert(j);
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
        let mut out = Vec::new();
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
            out.extend(x.iter().copied().enumerate().filter(|&(_, v)| v != R::ZERO));
            if let Some(o) = ordered {
                self.uorder_done(o);
            }
        } else {
            let ordered = ordered.unwrap();
            for &p in &ordered {
                if touched.take_sparse(p) && x[p] != R::ZERO {
                    out.push((p, x[p]));
                }
            }
            self.uorder_done(ordered);
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
            col_surv,
            col_sa: sa,
            col_sb: sb,
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

    /// `ftran_with` for a dense L-stage vector `x` (`nz` entries with `v != 0`, every other slot +0.0) that takes the dense U routine
    /// (`0x121ef0`): `x` is the U work vector and then the column snapshot, with no spike list, re-scatter or second zero-filled vector.
    /// The tail loop and the final pass are those of `ftran_with` for `u_dense`; a position the old marker set did not contain holds
    /// +0.0 and comes out as +0.0 from the final pass.
    #[allow(clippy::too_many_arguments)]
    fn ftran_with_dense_old(
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
        let mut col_surv = std::mem::take(&mut self.f2.get_mut().col_surv);
        col_surv.clear();
        let (mut sa, mut sb, mut su) = (0usize, 0usize, 0usize);
        for s in (self.i21..m).rev() {
            let c = self.at[s];
            let mut xc = x[c];
            if xc.abs() <= R::c(TOL) {
                x[c] = R::ZERO;
                continue;
            }
            xc /= self.diag[c];
            if xc.abs() <= R::c(TOL) {
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
            }
        }
        #[cfg(any(test, debug_assertions, feature = "factor-shadow"))]
        assert_eq!(
            (self.i21..m).filter(|&s| x[self.at[s]] != R::ZERO).count(),
            col_surv.len(),
            "ftran_with: nonzero tail slots differ from the processed survivors"
        );
        // The tail slots are final; every other slot is zeroed at or below the tolerance and divided otherwise (the old marker set held
        // the L entries and the scatter targets; every other slot is +0.0 and stays so).
        #[cfg(any(test, debug_assertions, feature = "factor-shadow"))]
        assert!((0..m).all(|p| (self.kind[p] == -2) == (self.slot[p] >= self.i21)), "tail kind and slot differ");
        for (p, xp) in x.iter_mut().enumerate() {
            let v = *xp;
            if v.abs() <= R::c(TOL) {
                *xp = R::ZERO;
            } else if self.kind[p] != -2 {
                let q = v / self.diag[p];
                *xp = if q.abs() <= R::c(TOL) { R::ZERO } else { q };
            }
        }
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
        let out: Vec<(usize, R)> = x.iter().copied().enumerate().filter(|&(_, v)| v != R::ZERO).collect();
        *self.f2.get_mut() = F2 {
            col_n: u_cnt,
            col_in: inp,
            x_col: x,
            dirty_blocks: None,
            col_surv,
            col_sa: sa,
            col_sb: sb,
            ..F2::<R>::default()
        };
        self.spike = spike;
        self.spike_used = false;
        self.last_alpha = out.len();
        out
    }
}
