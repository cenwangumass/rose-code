//! The pre-S3-6 `update` (per-call temporaries), frozen verbatim as the oracle of `update_without_temporaries_matches_original`.
use super::*;

impl<R: Real> DualFactor<R> {
    /// Replace a numeric U column and update its ordered row copy without doing any arithmetic.
    fn replace_u_column_old(&mut self, c: usize, col: Vec<(usize, R)>) {
        let old = std::mem::replace(&mut self.uc[c], col);
        for &(row, _) in &old {
            let entries = &mut self.ut[row];
            let first = entries.partition_point(|e| e.0 < c);
            let end = first + entries[first..].partition_point(|e| e.0 == c);
            entries.drain(first..end);
        }
        for &(row, v) in &self.uc[c] {
            let pos = self.ut[row].partition_point(|e| e.0 <= c);
            self.ut[row].insert(pos, (c, v));
        }
    }

    /// Replace the column at position `r` by the entering column of the last `ftran`.
    pub(super) fn update_old(&mut self, r: usize) {
        let m = self.m;
        let mut scratch = self.btran_ws.update.take();
        if crate::dump::on("fstruct") {
            let n = crate::dump::FACTOR_STRUCT_CALLS.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
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
            for &y in self.a[r].iter().chain(self.b[r].iter()) {
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
            conv = self
                .empty_tail
                .iter()
                .copied()
                .max_by_key(|&x| self.slot[x]);
            if let Some(x) = conv {
                self.empty_tail.remove(&x);
                convslot = Some(self.slot[x]);
                // Only a column some tail column refers to joins the nucleus list; the others become singletons.
                let refs = &mut scratch.conversion_refs;
                refs.clear();
                refs.extend(self.rowl[x].iter().copied().filter(|&k| self.kind[k] == -2 && self.a[k].contains(&x)));
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
            self.replace_u_column_old(r, Vec::new());
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
        let refcols = &mut scratch.refcols;
        refcols.clear();
        refcols.extend(self.rowl[r].iter().copied().map(|k| (k, self.slot[k])));
        refcols.sort_unstable();
        refcols.dedup();
        refcols.retain(|&(x, _)| x != r && self.kind[x] == -2 && (self.a[x].contains(&r) || self.b[x].contains(&r)));
        // Only a surviving row-r reference can raise the mark after the slot walk.
        let oldslot = self.slot[r];
        let top: i64 = av.iter().map(|&x| self.slot[x] as i64).max().unwrap_or(-1);
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
                    let row = &mut self.ut[r];
                    let first = row.partition_point(|e| e.0 < k);
                    let end = first + row[first..].partition_point(|e| e.0 == k);
                    row.drain(first..end);
                    eliminated = true;
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
        if eliminated {
            let mut ind = vec![r];
            let mut val = vec![R::ONE / d0];
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
                // The transpose includes duplicates in the same multiplicity as the column scan.
                // Every entry of `ut` belongs to a tail column (audited) and `r`'s own entries sit only in `old_rows`, skipped above.
                let count = self.ut[row].len();
                let counted = || self.ut[row].iter().filter(|&&(c, _)| c != r && self.kind[c] == -2).count();
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
        self.diag[r] = newdiag;
        self.replace_u_column_old(r, self.spike.iter().copied().filter(|e| e.0 != r).collect());
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
        scratch.sp = sp;
        scratch.clear();
        self.btran_ws.update.replace(scratch);
    }
}
