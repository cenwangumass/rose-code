//! Structural mirror of the reference sparse FTRAN: the order in which it lists the nonzeros of `B^-1 a` (`reverse-engineer/specs/ftran-order.md`).
//!
//! Positions are the basis positions after the LU renumbering (the pivot step). A fresh factor has singleton positions (kind -1),
//! a nucleus list of diagonal-only positions (phase 3 pivots, kind = index in the list) and a triangular tail (phase 4 and the
//! Markowitz nucleus, slots `i21..m`) whose columns carry U entries split into tail positions (`a`) and nucleus-list positions (`b`).
//! Column replacements follow the update rules recovered from the GDB dumps. The factor is mirrored numerically as well (L multipliers,
//! U columns, diagonal, row etas, including the reference's lane-wise eta dot products), so the solve returns the reference's values
//! bit for bit and drops the entries that cancel exactly as it does.

use crate::lu_order::Structure;
use std::collections::BinaryHeap;

pub struct RefStruct {
    m: usize,
    kind: Vec<i64>,
    nuc: Vec<usize>,
    i21: usize,
    a: Vec<Vec<usize>>,
    b: Vec<Vec<usize>>,
    lcols: Vec<Vec<usize>>,
    perm: Vec<usize>,
    dens: f64,
    /// Density threshold fixed at factorization time (`f104` in the notes): 2, or `5 (nnzL/m + 1)^2` capped at `m`.
    dthr: usize,
    slot: Vec<usize>,
    at: Vec<usize>,
    etas: Vec<Eta>,
    lvals: Vec<Vec<f64>>,
    /// Longest-path level of each L node, when the reference solves L level by level (`0x12a990`): only for `m >= 250`, a
    /// small level count (`< m/10`) and a short Markowitz L part. Otherwise the L solve is the Kahn FIFO.
    levels: Option<Vec<usize>>,
    diag: Vec<f64>,
    /// Numeric U entries of each tail column `(position, value)`.
    uc: Vec<Vec<(usize, f64)>>,
    /// Spike of the most recent solve: the L stage and eta result `(position, value)` in list order.
    spike: Vec<(usize, f64)>,
    /// Refactorization trigger inputs (`0x142350`): updates since the factorization, accumulated solve result counts, accumulated
    /// eta lengths, and the factor size `nnz(L) + nnz(U) + 2 m` fixed at factorization time.
    upd: usize,
    sum_alpha: f64,
    sum_eta: f64,
    last_alpha: usize,
    znz: f64,
    /// Work-counter inputs (`crate::work`): phase-1 steps (`L+0xb8`), singleton-phase steps (`L+0xc4`), `nnz(L)`, and for each
    /// position the first row eta that lists it (`F+0x1a0`).
    n1: usize,
    loff: usize,
    nnz_l: usize,
    first_eta: Vec<usize>,
    /// Estimate of the most recent `ftran` (`0x120740`) or `btran` (`0x12de60`), before the factor object's scale 1.3.
    pub work: f64,
    /// Fill of the reference's pool for etas and the row-wise copy of U (`F+0x160`), and the capacity of each row of that copy.
    /// Only the parity of `fill` matters: it decides the alignment of the next eta's value array.
    fill: usize,
    row_cap: Vec<usize>,
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
    last_dense: bool,
    /// Whether the most recent `btran` returned a dense result (count `-1`): its L^T solve ran dense.
    pub bt_dense: bool,
}

struct Eta {
    idx: usize,
    ind: Vec<usize>,
    val: Vec<f64>,
    /// Elements before the first 16-byte boundary of the value array in the reference's pool (the parity of the pool fill at which
    /// the eta was stored).
    peel: usize,
}

const TOL: f64 = 1e-13;

/// Eta dot product as the reference evaluates it: the vectorised loop peels the elements before the first 16-byte boundary of
/// the value array (`peel` is 0 or 1) and keeps eight lane-wise partial sums; below `peel + 8` terms it is sequential.
fn dot(x: &[f64], c: &[f64], peel: usize) -> f64 {
    let n = c.len();
    let t = |i: usize| x[i] * c[i];
    if n < peel + 8 {
        let mut acc = 0.0;
        for i in 0..n {
            acc += t(i);
        }
        return acc;
    }
    let n8 = n - ((n - peel) & 7);
    let mut lanes = [0.0f64; 8];
    for i in 0..peel {
        lanes[0] += t(i);
    }
    for i in peel..n8 {
        lanes[(i - peel) & 7] += t(i);
    }
    // The compiled reduction adds the four packed accumulators pairwise: (l0+l2 | l1+l3) + (l4+l6 | l5+l7), then the two halves.
    let mut tot = ((lanes[0] + lanes[2]) + (lanes[4] + lanes[6])) + ((lanes[1] + lanes[3]) + (lanes[5] + lanes[7]));
    for i in n8..n {
        tot += t(i);
    }
    tot
}

impl RefStruct {
    /// Fresh factor from the pivot sequence; `cols[pos]` is the basis column at incoming position `pos` as `(is_slack, rows, values)`.
    pub fn fresh(m: usize, ncols: usize, st: &Structure, cols: &[(bool, Vec<usize>, Vec<f64>)]) -> RefStruct {
        let mut step = vec![usize::MAX; m];
        for (s, &(r, _)) in st.pivots.iter().enumerate() {
            step[r] = s;
        }
        let n1 = st.phase.iter().filter(|&&p| p == 1 || p == 2).count();
        let n3 = st.phase.iter().filter(|&&p| p == 3).count();
        let nseq = st.phase.iter().filter(|&&p| p != 0).count();
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
        let lcols: Vec<Vec<usize>> = st.l_rows.iter().map(|rows| rows.iter().map(|&r| step[r]).collect()).collect();
        let lvals = st.l_vals.clone();
        // pivot-row values of the Markowitz steps, by column
        let mut mk: Vec<Vec<(usize, f64)>> = vec![Vec::new(); m];
        for s in nseq..m {
            for &(c, v) in &st.u_rows[s] {
                mk[c].push((s, v));
            }
        }
        let mut a = vec![Vec::new(); m];
        let mut b = vec![Vec::new(); m];
        let mut uc: Vec<Vec<(usize, f64)>> = vec![Vec::new(); m];
        let (mut sa, mut sb) = (0usize, 0usize);
        let mut xv = vec![0.0f64; m];
        let mut cnt = vec![0i32; m];
        for t in i21..m {
            let c = st.pivots[t].1;
            let (_, idx, val) = &cols[c];
            let mut start: Vec<usize> = Vec::new();
            for (&i, &v) in idx.iter().zip(val) {
                let s = step[i];
                if s <= t {
                    xv[s] += v;
                    start.push(s);
                }
            }
            let mut lst = start.clone();
            for &v in &lst {
                cnt[v] += 1;
            }
            let mut k = 0;
            while k < lst.len() {
                let v = lst[k];
                k += 1;
                for &y in &lcols[v] {
                    if y > t {
                        continue;
                    }
                    if cnt[y] == 0 {
                        lst.push(y);
                    }
                    cnt[y] += 1;
                }
            }
            let mut q: Vec<usize> = Vec::new();
            for &v in &lst[..start.len()] {
                cnt[v] -= 1;
                if cnt[v] == 0 {
                    q.push(v);
                }
            }
            let mut k = 0;
            while k < q.len() {
                let v = q[k];
                k += 1;
                for &y in &lcols[v] {
                    if y > t {
                        continue;
                    }
                    cnt[y] -= 1;
                    if cnt[y] == 0 {
                        q.push(y);
                    }
                }
            }
            for &v in &q {
                if v >= t {
                    continue;
                }
                let x = xv[v];
                if x == 0.0 {
                    continue;
                }
                for (e, &y) in lcols[v].iter().enumerate() {
                    if y > t {
                        continue;
                    }
                    xv[y] -= lvals[v][e] * x;
                }
            }
            let mut ents: Vec<(usize, f64)> = Vec::new();
            let mut touched: Vec<usize> = lst.clone();
            touched.sort_unstable();
            touched.dedup();
            for &v in &touched {
                if v < t && xv[v].abs() > TOL {
                    ents.push((v, xv[v]));
                }
            }
            for &(sidx, vv) in &mk[c] {
                if sidx < t {
                    if let Some(e) = ents.iter_mut().find(|e| e.0 == sidx) {
                        e.1 = vv;
                    } else {
                        ents.push((sidx, vv));
                    }
                }
            }
            for &v in &lst {
                xv[v] = 0.0;
                cnt[v] = 0;
            }
            for &v in &touched {
                xv[v] = 0.0;
            }
            ents.sort_unstable_by_key(|e| e.0);
            let av: Vec<usize> = ents.iter().map(|e| e.0).filter(|&x| kind[x] == -2).collect();
            let mut bv: Vec<usize> = ents.iter().map(|e| e.0).filter(|&x| kind[x] >= 0).collect();
            bv.sort_unstable_by_key(|&x| std::cmp::Reverse(kind[x]));
            sa += av.len();
            sb += bv.len();
            a[t] = av;
            b[t] = bv;
            uc[t] = ents;
        }
        let nt = m - i21;
        let nn = nuc.len();
        let h33 = if nt > 0 { (sa + nt) as f64 / nt as f64 } else { 1.5 };
        let h34 = if nn > 0 { (sb + nn) as f64 / nn as f64 } else { 1.5 };
        let nnz: usize = lcols.iter().map(|c| c.len()).sum();
        let dthr = if nnz <= m {
            2
        } else {
            let t = 5.0 * ((nnz as f64 / m as f64) + 1.0).powi(2);
            if (m as f64) > t { t as usize } else { m }
        };
        let mut lev_dense = false;
        let levels = {
            let d = m - nseq;
            let nnz_tail: usize = lcols[nseq..].iter().map(|c| c.len()).sum();
            let outer = d < 6 || 10 * d * d <= m || nnz_tail <= 10 * d;
            if m >= 250 && outer {
                let mut lev = vec![0usize; m];
                for v in 0..m {
                    for &y in &lcols[v] {
                        lev[y] = lev[y].max(lev[v] + 1);
                    }
                }
                let maxlev = lev.iter().copied().max().unwrap_or(0);
                lev_dense = (maxlev as f64) >= 0.1 * m as f64;
                if !lev_dense { Some(lev) } else { None }
            } else {
                None
            }
        };
        let nnz_l_tot: usize = lcols.iter().map(|c| c.len()).sum();
        let perm = step;
        let znz = (nnz + sa + sb + 2 * m) as f64;
        RefStruct {
            upd: 0,
            sum_alpha: 0.0,
            sum_eta: 0.0,
            last_alpha: 0,
            znz,
            n1: st.phase.iter().filter(|&&p| p == 1).count(),
            loff: nseq,
            nnz_l: nnz_l_tot,
            first_eta: vec![usize::MAX; m],
            work: 0.0,
            scratch_peel: ((m as f64).powf(0.7) as usize).max(20).wrapping_add(ncols + 1) % 2,
            fill: 2 * (uc.iter().map(Vec::len).sum::<usize>() + nnz_l_tot),
            row_cap: {
                let mut c = vec![0usize; m];
                for col in &uc {
                    for &(row, _) in col {
                        c[row] += 1;
                    }
                }
                c
            },
            upd_work: Vec::new(),
            h33,
            h34,
            lev_dense,
            last_dense: false,
            bt_dense: false,
            m,
            kind,
            nuc,
            i21,
            a,
            b,
            lcols,
            perm,
            dens: (h33 + h34) * (h33 + h34),
            dthr,
            slot: (0..m).collect(),
            at: (0..m).collect(),
            etas: Vec::new(),
            lvals,
            levels,
            diag: st.diag.clone(),
            uc,
            spike: Vec::new(),
        }
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
            for &x in &self.lcols[v] {
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
            for &x in &self.lcols[v] {
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
    fn lstage_num(&self, rows: &[usize], vals: &[f64], force_dense: bool) -> (Vec<(usize, f64)>, usize, f64) {
        let m = self.m;
        let mut e = 0.0f64;
        let mut x = vec![0.0f64; self.m];
        let mut start: Vec<usize> = Vec::with_capacity(rows.len());
        for (&r, &v) in rows.iter().zip(vals) {
            let s = self.perm[r];
            x[s] += v;
            start.push(s);
        }
        let dense_l = rows.len() * self.dthr > self.m || force_dense;
        let (mut order, count_l) = if dense_l {
            // Dense L solve (`0x12bd30`): every node in ascending order; the list handed on is the ascending scan.
            e = m as f64 + rows.len() as f64 * 2.0 + e;
            e = (m - self.loff) as f64 * 2.0 + self.n1 as f64 + self.n1 as f64 + e;
            for v in 0..self.m {
                let xv = x[v];
                if xv.abs() <= TOL {
                    continue;
                }
                e = self.lcols[v].len() as f64 * 2.0 + e;
                for (e, &y) in self.lcols[v].iter().enumerate() {
                    x[y] -= self.lvals[v][e] * xv;
                }
            }
            let count_l = (0..self.m).filter(|&v| x[v].abs() > TOL).count();
            ((0..self.m).collect(), count_l)
        } else if let Some(lev) = &self.levels {
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
                    if xv.abs() >= TOL {
                        out.push(v);
                        for (e, &y) in self.lcols[v].iter().enumerate() {
                            if !active[y] {
                                active[y] = true;
                                buckets[lev[y]].push(y);
                                remaining += 1;
                            }
                            x[y] -= self.lvals[v][e] * xv;
                        }
                    } else {
                        x[v] = 0.0;
                    }
                }
                l += 1;
            }
            let n = out.len();
            e = ((self.nnz_l + m + 2 * m) as f64 * 5.0 * n as f64) / (m + 1) as f64 + e;
            (out, n)
        } else {
            let q = self.kahn(&start);
            for &v in &q {
                let xv = x[v];
                if xv == 0.0 {
                    continue;
                }
                for (e, &y) in self.lcols[v].iter().enumerate() {
                    x[y] -= self.lvals[v][e] * xv;
                }
            }
            let count_l = q.iter().filter(|&&v| x[v].abs() > TOL).count();
            let d = m - self.loff;
            let nnz_tail: usize = self.lcols[self.loff..].iter().map(|c| c.len()).sum();
            if !(d < 6 || 10 * d * d <= m || nnz_tail <= 10 * d) && std::env::var_os("ROSE_WORK").is_some() {
                eprintln!("WARN partial L solve (0x12a990 third strategy) not modelled");
            }
            e = (((self.nnz_l + m) as f64 + m as f64 + m as f64) * 7.0 * count_l as f64) / (m + 1) as f64 + e;
            (q, count_l)
        };
        let mut have = vec![false; self.m];
        for &v in &order {
            have[v] = true;
        }
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
                    if x[v].abs() > TOL && !inlist[v] {
                        inlist[v] = true;
                        c0 += 1;
                        first = first.min(self.first_eta[v]);
                    }
                }
                let mut xe = x.clone();
                let mut ee = (netas - first) as f64 * 3.0;
                for eta in &self.etas[first..] {
                    ee = eta.ind.len() as f64 * 2.0 + ee;
                    let xs: Vec<f64> = eta.ind.iter().map(|&j| xe[j]).collect();
                    let mut sv = dot(&xs, &eta.val, eta.peel);
                    if sv.abs() <= TOL {
                        sv = 0.0;
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
            let xs: Vec<f64> = eta.ind.iter().map(|&j| x[j]).collect();
            let mut s = dot(&xs, &eta.val, eta.peel);
            if s.abs() <= TOL {
                s = 0.0;
            }
            if !have[eta.idx] && s != 0.0 {
                order.push(eta.idx);
                have[eta.idx] = true;
            }
            x[eta.idx] = s;
        }
        (order.into_iter().filter(|&v| x[v].abs() > TOL).map(|v| (v, x[v])).collect(), count_l, e)
    }

    /// FTRAN of a sparse column: the nonzeros `(position, value)` of `B^-1 a` in the reference's output order. The spike is kept
    /// for the following `update`.
    pub fn ftran(&mut self, rows: &[usize], vals: &[f64]) -> Vec<(usize, f64)> {
        self.ftran_with(rows, vals, false)
    }

    /// FTRAN of the right-hand side after a factorization: the reference always runs the dense L solve there, whatever the
    /// density, which fixes the order of the floating-point updates.
    pub fn ftran_basic(&mut self, rows: &[usize], vals: &[f64]) -> Vec<(usize, f64)> {
        self.ftran_with(rows, vals, true)
    }

    fn ftran_with(&mut self, rows: &[usize], vals: &[f64], force_dense: bool) -> Vec<(usize, f64)> {
        let m = self.m;
        let (spike, count_l, e_l) = self.lstage_num(rows, vals, force_dense);
        let mut x = vec![0.0f64; m];
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
        tail.sort_unstable_by_key(|&p| std::cmp::Reverse(self.slot[p]));
        for &c in &tail {
            let xc = x[c];
            if xc.abs() <= TOL {
                x[c] = 0.0;
                continue;
            }
            let xc = xc / self.diag[c];
            if xc.abs() <= TOL {
                x[c] = 0.0;
                continue;
            }
            x[c] = xc;
            for &(j, v) in &self.uc[c] {
                x[j] -= v * xc;
                touched[j] = true;
            }
        }
        for p in 0..m {
            if touched[p] && self.kind[p] != -2 {
                let v = x[p] / self.diag[p];
                x[p] = if x[p].abs() <= TOL || v.abs() <= TOL { 0.0 } else { v };
            }
        }
        // A dense result (the reference's `-1` count) is scanned by position index; a sparse one comes in the U solve's order.
        let dense = rows.len() * self.dthr > m
            || (!self.etas.is_empty() && count_l * self.dthr > m)
            || !(self.nuc.len() + spike.len() < self.i21 || self.dens * (spike.len() as f64) < 2.0 * m as f64);
        let inp: Vec<usize> = spike.iter().map(|e| e.0).collect();
        let vec_dense = rows.len() * self.dthr > m || (!self.etas.is_empty() && count_l * self.dthr > m);
        self.work = self.u_work(&inp, &x, !dense, if vec_dense { None } else { Some(spike.len()) }) + e_l;
        self.last_dense = dense;
        let ord = if dense { (0..m).collect() } else { self.uorder(&inp, &x) };
        let mut out: Vec<(usize, f64)> = Vec::new();
        let mut done = vec![false; m];
        for &p in &ord {
            if !done[p] {
                done[p] = true;
                if x[p] != 0.0 {
                    out.push((p, x[p]));
                }
            }
        }
        for p in 0..m {
            if touched[p] && !done[p] && x[p] != 0.0 {
                out.push((p, x[p]));
            }
        }
        // The column replacement receives the spike in ascending order after a dense solve (the reference scans the dense vector),
        // and in L-stage order after a sparse one.
        let mut spike = spike;
        if dense {
            spike.sort_unstable_by_key(|e| e.0);
        }
        self.spike = spike;
        self.last_alpha = out.len();
        out
    }

    /// The reference's refactorization decision, evaluated before an iteration's own update: at the update limit
    /// `max(20, m^0.7)`, or from `U >= m^0.45` when the accumulated work `(1.5 S + 2 c) m` exceeds `min((Z - 1.5 m)^2, 4 m^2 U0)`.
    /// `(nucleus list length, first tail slot, dens)`; for diagnostics against the reference's `i3`, `i21` and `(hdr33+hdr34)^2`.
    /// The fresh factor's numerics for the `luout` dump: the diagonal and each position's L column `(position, value)`.
    pub fn numerics(&self) -> (&[f64], Vec<Vec<(usize, f64)>>) {
        let l = (0..self.m).map(|v| self.lcols[v].iter().copied().zip(self.lvals[v].iter().copied()).collect()).collect();
        (&self.diag, l)
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

    pub fn refactor_due(&self) -> bool {
        let m = self.m as f64;
        let maxupd = (m.powf(0.7) as usize).max(20);
        if self.upd == maxupd {
            return true;
        }
        let c46c = m.powf(0.45) as usize;
        if self.upd < c46c {
            return false;
        }
        let d6 = (self.znz - 1.5 * m).powi(2);
        let d7 = (4.0 * m * m * c46c as f64).min(d6);
        (1.5 * self.sum_alpha + 2.0 * self.sum_eta) * m - d7 > 0.0
    }

    /// The reference BTRAN (`0x12de60`) of a sparse input given by basis position: U^T solve (`0x130b70` or the dense loop), the row
    /// etas in reverse, the L^T solve, then the rows. Returns `(row, value)` with the reference's drops, and leaves its work
    /// estimate in `work`. `direct` is the update-saving entry (`0x131710`), which goes to the sparse U^T solve without the density
    /// test unless the first position's U^T row is empty. `dense_in` passes the input as a dense vector (count `-1`).
    pub fn btran(&mut self, inp: &[(usize, f64)], direct: bool, dense_in: bool) -> Vec<(usize, f64)> {
        let m = self.m;
        let mut ut: Vec<Vec<(usize, f64)>> = vec![Vec::new(); m];
        for c in 0..m {
            for &(j, v) in &self.uc[c] {
                ut[j].push((c, v));
            }
        }
        let mut x = vec![0.0f64; m];
        let mut e = 0.0f64;
        let cnt = inp.len();
        self.bt_dense = false;
        let hsum = (self.h33 + self.h34) as usize;
        let direct = direct && !inp.is_empty() && !ut[inp[0].0].is_empty();
        // The current vector: a list (`Some`) or dense (`None`).
        let mut list: Option<Vec<usize>>;
        if !dense_in && (direct || hsum * cnt <= m) {
            let mut out: Vec<usize> = Vec::new();
            let mut w = vec![0.0f64; m];
            let mut eu = 0.0f64;
            let push = |w: &mut Vec<f64>, row: &[(usize, f64)], d: f64| {
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
                let tail = if self.h33 > 5.0 { self.kind[p] == -2 } else { self.slot[p] >= self.i21 };
                if tail {
                    if self.h33 <= 5.0 || d.abs() > TOL {
                        lowest = lowest.min(self.slot[p]);
                        w[p] += v;
                    }
                } else if d.abs() > TOL {
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
                w[p] = 0.0;
                if d.abs() > TOL {
                    eu = ut[p].len() as f64 * 2.0 + eu;
                    out.push(p);
                    x[p] = d;
                    push(&mut w, &ut[p], d);
                }
            }
            e = eu + e;
            if out.is_empty() {
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
                if xv.abs() <= TOL || d.abs() <= TOL {
                    x[p] = 0.0;
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
                    let mut inl = vec![false; m];
                    for &p in l.iter() {
                        inl[p] = true;
                    }
                    let mut es = netas as f64 + netas as f64 + c0 as f64 * 5.0;
                    for eta in self.etas.iter().rev() {
                        let xr = x[eta.idx];
                        if xr == 0.0 {
                            continue;
                        }
                        es = eta.ind.len() as f64 * 3.0 + es;
                        for (k, &j) in eta.ind.iter().enumerate() {
                            if j == eta.idx {
                                x[j] = eta.val[k] * xr;
                            } else {
                                x[j] = eta.val[k] * xr + x[j];
                                if !inl[j] {
                                    inl[j] = true;
                                    l.push(j);
                                }
                            }
                        }
                    }
                    let n = l.len();
                    e = es + (n - c0) as f64 * 2.0 + 3.0 * n as f64 + e;
                    l.retain(|&p| {
                        if x[p].abs() >= TOL {
                            true
                        } else {
                            x[p] = 0.0;
                            false
                        }
                    });
                    if l.is_empty() {
                        self.work = e;
                        return Vec::new();
                    }
                }
                _ => {
                    for eta in self.etas.iter().rev() {
                        let xr = x[eta.idx];
                        if xr.abs() <= TOL {
                            x[eta.idx] = 0.0;
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
                    list = None;
                }
            }
        }
        // L^T: row i of L pushes into the earlier nodes.
        let mut lt: Vec<Vec<(usize, f64)>> = vec![Vec::new(); m];
        for v in 0..m {
            for (k, &y) in self.lcols[v].iter().enumerate() {
                lt[y].push((v, self.lvals[v][k]));
            }
        }
        let mut res: Vec<usize> = Vec::new();
        let sparse_l = match &list {
            Some(l) => l.len() * self.dthr <= m && !self.lev_dense,
            None => false,
        };
        if sparse_l {
            let l = list.as_ref().unwrap();
            let cnt = l.len();
            let (n1, nseq) = (self.n1, self.loff);
            let d = m - nseq;
            let nnz_tail: usize = self.lcols[nseq..].iter().map(|c| c.len()).sum();
            let mut hi = 0usize;
            if d >= 6 && m < 10 * d * d && n1 < 2 * d && 10 * d < nnz_tail && cnt > 0 {
                hi = l.iter().filter(|&&i| i >= nseq).count();
            }
            let dn = if d > 0 { nnz_tail as f64 * 0.5 / d as f64 } else { 0.0 };
            if hi as f64 * dn * dn > d as f64 {
                // Partial: inputs in the middle range directly, then the top range and the first-phase range, both descending.
                let mut el = cnt as f64 * 3.0;
                let mut w = vec![0.0f64; m];
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
                let ranges: Vec<usize> = if nseq <= maxi { (nseq..=maxi).rev().collect() } else { Vec::new() };
                for i in ranges.into_iter().chain((0..n1).rev()) {
                    let v = w[i];
                    if v.abs() > TOL {
                        el = lt[i].len() as f64 * 2.0 + el;
                        res.push(i);
                        for &(j, lv) in &lt[i] {
                            w[j] -= lv * v;
                        }
                    } else {
                        w[i] = 0.0;
                    }
                }
                x = w;
                e = el + e;
            } else {
                // Level buckets of the transposed graph when the factor has levels (`0x12bd30`), else a max-heap on the node index; every
                // reached node is kept when at least the tolerance.
                let mut seen = vec![false; m];
                if self.levels.is_some() {
                    // Level buckets of the transposed graph: a node comes after the nodes whose rows it appears in.
                    let mut lev = vec![0usize; m];
                    for v in (0..m).rev() {
                        for &y in &self.lcols[v] {
                            lev[v] = lev[v].max(lev[y] + 1);
                        }
                    }
                    let maxlev = lev.iter().copied().max().unwrap_or(0);
                    let mut buckets: Vec<Vec<usize>> = vec![Vec::new(); maxlev + 1];
                    let mut lo = usize::MAX;
                    for &i in l {
                        seen[i] = true;
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
                            if v.abs() >= TOL {
                                res.push(i);
                                for &(j, lv) in &lt[i] {
                                    x[j] -= lv * v;
                                }
                                // The reference walks the row from its last entry.
                                for &(j, _) in lt[i].iter().rev() {
                                    if !seen[j] {
                                        seen[j] = true;
                                        buckets[lev[j]].push(j);
                                    }
                                }
                            } else {
                                x[i] = 0.0;
                            }
                        }
                        lv_i += 1;
                    }
                } else {
                    let mut heap: BinaryHeap<usize> = BinaryHeap::new();
                    for &i in l {
                        seen[i] = true;
                        heap.push(i);
                    }
                    while let Some(i) = heap.pop() {
                        let v = x[i];
                        if v.abs() >= TOL {
                            res.push(i);
                            for &(j, lv) in &lt[i] {
                                x[j] -= lv * v;
                                if !seen[j] {
                                    seen[j] = true;
                                    heap.push(j);
                                }
                            }
                        } else {
                            x[i] = 0.0;
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
                if v.abs() > TOL {
                    e = lt[i].len() as f64 * 2.0 + e;
                    res.push(i);
                    for &(j, lv) in &lt[i] {
                        x[j] -= lv * v;
                    }
                }
            }
        }
        self.work = e;
        let mut row_of = vec![0usize; m];
        for (r, &v) in self.perm.iter().enumerate() {
            row_of[v] = r;
        }
        if self.bt_dense {
            // The dense result keeps the sub-tolerance entries the L^T pass left in place.
            // and lists them by ascending row.
            return (0..m).filter_map(|r| { let i = self.perm[r]; (x[i] != 0.0).then(|| (r, x[i])) }).collect();
        }
        res.into_iter().filter(|&i| x[i] != 0.0).map(|i| (row_of[i], x[i])).collect()
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

    /// Whether the fresh factor's L solve goes level by level (`0x12a990`).
    pub fn has_levels(&self) -> bool {
        self.levels.is_some()
    }

    pub fn last_dense(&self) -> bool {
        self.last_dense
    }

    /// Work estimate of an FTRAN with a dense input (the reference's count `-1`, e.g. the basic solution after a refactorization):
    /// dense L solve, all row etas, dense U solve. `b` is indexed by row.
    pub fn ftran_dense_work(&self, b: &[f64]) -> f64 {
        let m = self.m;
        let mut x = vec![0.0f64; m];
        for r in 0..m {
            x[self.perm[r]] += b[r];
        }
        let mut e = (m - self.loff) as f64 * 2.0 + self.n1 as f64 + self.n1 as f64;
        for v in 0..m {
            let xv = x[v];
            if xv.abs() <= TOL {
                continue;
            }
            e = self.lcols[v].len() as f64 * 2.0 + e;
            for (k, &y) in self.lcols[v].iter().enumerate() {
                x[y] -= self.lvals[v][k] * xv;
            }
        }
        if !self.etas.is_empty() {
            e = self.sum_eta * 2.0 + e;
            for eta in &self.etas {
                let xs: Vec<f64> = eta.ind.iter().map(|&j| x[j]).collect();
                let sv = dot(&xs, &eta.val, eta.peel);
                x[eta.idx] = if sv.abs() <= TOL { 0.0 } else { sv };
            }
        }
        for t in (self.i21..m).rev() {
            let c = self.at[t];
            let xc = x[c] / self.diag[c];
            if x[c].abs() <= TOL || xc.abs() <= TOL {
                x[c] = 0.0;
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
    fn u_work(&self, inp: &[usize], x: &[f64], sparse: bool, cnt: Option<usize>) -> f64 {
        let m = self.m;
        let i21 = self.i21;
        let i3 = self.nuc.len();
        let len = |c: usize| self.uc[c].len();
        // Every tail slot from the top down, then the nucleus list (shared by the dense paths).
        let full_tail = |mut e: f64| {
            e = (m - i21) as f64 * 3.0 + e;
            for t in i21..m {
                let c = self.at[t];
                if x[c] != 0.0 {
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
                if x[c] != 0.0 {
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
        let tail: Vec<usize> = inp.iter().copied().filter(|&p| self.kind[p] == -2).collect();
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
            if x[c] != 0.0 {
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

    /// Order of the nonzero positions of the solution, given the L-stage list.
    pub fn uorder(&self, inp: &[usize], x: &[f64]) -> Vec<usize> {
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
        let sing: Vec<usize> = inp_d.iter().copied().filter(|&x| self.kind[x] == -1).collect();
        let tri: Vec<usize> = inp_d.iter().copied().filter(|&x| self.kind[x] == -2).collect();
        let nuc: Vec<usize> = inp_d.iter().copied().filter(|&x| self.kind[x] >= 0).collect();
        let i3 = self.nuc.len();
        let d2 = self.dens;
        let by_slot_desc = |v: &mut Vec<usize>| v.sort_unstable_by_key(|&x| std::cmp::Reverse(self.slot[x]));
        let nuc_sorted = |set: Vec<usize>| {
            let mut v = set;
            v.sort_unstable_by_key(|&x| self.kind[x]);
            v
        };
        let mut have = vec![false; m];
        let dense = (m as f64) <= (n as f64) * d2 || d2 * ((nuc.len() + tri.len()) as f64) > (i3 + (m - self.i21)) as f64;
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
            if x[k] != 0.0 {
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
        for &y in nuc.iter().chain(order.iter().filter(|&&k| x[k] != 0.0).flat_map(|&k| self.b[k].iter())) {
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
        let len = |s: &Self, pos: usize| s.a[pos].iter().chain(&s.b[pos]).filter(|&&x| x != r).count();
        // The reference's work: each jump scans the column (2 per entry and 5) on the factor's scale; the final shift
        // moves the slots in between.
        let mut jumps = 0.0;
        while cur != hole {
            let pos = self.at[cur];
            if len(self, pos) >= cur - hole {
                self.upd_work.push((0xf3d83, 2.0 * (cur - hole) as f64, crate::work::CTX));
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
            let at = self.upd_work.len() - usize::from(self.upd_work.last().is_some_and(|w| w.0 == 0xf3d83));
            self.upd_work.insert(at, (0xf3ca4, jumps, crate::work::FAC));
        }
    }

    /// Row elimination of a column replacement at `r`: the eta multipliers `(position, multiplier)` in elimination order.
    fn eliminate(&self, r: usize) -> Vec<(usize, f64)> {
        let m = self.m;
        let d0 = self.diag[r];
        let mut rows: Vec<Vec<(usize, f64)>> = vec![Vec::new(); m];
        for k in 0..m {
            if self.kind[k] == -2 {
                for &(row, v) in &self.uc[k] {
                    rows[row].push((k, v));
                }
            }
        }
        let mut w = vec![0.0f64; m];
        let mut seen = vec![false; m];
        let mut heap: std::collections::BinaryHeap<std::cmp::Reverse<(usize, usize)>> = std::collections::BinaryHeap::new();
        for &(k, v) in &rows[r] {
            w[k] = v * (1.0 / d0);
            seen[k] = true;
            heap.push(std::cmp::Reverse((self.slot[k], k)));
        }
        let mut eta = Vec::new();
        while let Some(std::cmp::Reverse((_, k))) = heap.pop() {
            let c = -w[k] / self.diag[k];
            if c.abs() <= TOL {
                continue;
            }
            eta.push((k, c));
            for &(l, v) in &rows[k] {
                w[l] += c * v;
                if !seen[l] {
                    seen[l] = true;
                    heap.push(std::cmp::Reverse((self.slot[l], l)));
                }
            }
        }
        eta
    }

    /// Replace the column at position `r` by the entering column of the last `ftran`.
    pub fn update(&mut self, r: usize) {
        let m = self.m;
        if crate::dump::on("fstruct") {
            static CALLS: std::sync::atomic::AtomicUsize = std::sync::atomic::AtomicUsize::new(0);
            let n = CALLS.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
            crate::dump::emit("fstruct", serde_json::json!({"n": n, "m": m, "i3": self.nuc.len(), "i21": self.i21, "kind": self.kind, "nuc": self.nuc, "seqA": self.slot, "fill": self.fill, "diag": self.diag, "uvals": (0..m).filter(|&k| self.kind[k] == -2).map(|k| (k.to_string(), (self.uc[k].iter().filter(|e| self.kind[e.0] == -2).map(|e| e.0).collect::<Vec<_>>(), self.uc[k].iter().filter(|e| self.kind[e.0] == -2).map(|e| e.1).collect::<Vec<_>>()))).collect::<std::collections::BTreeMap<_, _>>(), "cols": (0..m).filter(|&k| self.kind[k] == -2).map(|k| (k.to_string(), [&self.a[k], &self.b[k]])).collect::<std::collections::BTreeMap<_, _>>()}));
        }
        self.upd += 1;
        self.upd_work.clear();
        // The eta is stored first, at the current pool fill; its value array sits on a 16-byte boundary when the fill is even.
        let peel = self.fill % 2;
        // The old column's entries leave the row-wise copy after the new column's entries have been added.
        let old_rows: Vec<usize> = if self.kind[r] == -2 { self.uc[r].iter().map(|e| e.0).collect() } else { Vec::new() };
        // The new U column's entries (`0xf15c0`, sites `0xf2521` and `0xf2e3c`).
        let k = self.spike.iter().filter(|e| e.0 != r).count();
        if k > 0 {
            self.upd_work.push((0xf2521, 3.0 * k as f64, crate::work::CTX));
        }
        // A tail column leaves the row-wise copy of U entry by entry.
        if self.kind[r] == -2 && !self.uc[r].is_empty() {
            let len = self.uc[r].len();
            self.upd_work.push((if len < k { 0xf2044 } else { 0xf207a }, 5.0 * len as f64, crate::work::CTX));
        }
        self.sum_alpha += self.last_alpha as f64;
        let spike: Vec<usize> = self.spike.iter().map(|e| e.0).collect();
        let spike = &spike[..];
        let eta_cand = self.eliminate(r);
        let mut spk = vec![0.0f64; m];
        for &(p, v) in &self.spike {
            spk[p] = v;
        }
        let d0 = self.diag[r];
        let mut sp: Vec<usize> = Vec::with_capacity(spike.len());
        for &x in spike {
            if !sp.contains(&x) {
                sp.push(x);
            }
        }
        if let Some(i) = sp.iter().position(|&x| x == r) {
            // The reference keeps the spike as a tail group and a nucleus group; the leaving column is deleted from the nucleus
            // group by moving that group's last entry into its place.
            let last = sp.pop().unwrap();
            if last != r {
                sp[i] = last;
            }
        }
        let old = self.kind[r];
        let mut convslot: Option<usize> = None;
        let mut conv: Option<usize> = None;
        // The reference demotes a bare tail column only when the new column depends on the tail; with no tail entries in the spike
        // the tail simply grows by one slot.
        if old != -2 && self.slot[r] + 1 != self.i21 && sp.iter().any(|&y| y != r && self.kind[y] == -2) {
            conv = (0..m).filter(|&x| self.kind[x] == -2 && self.a[x].is_empty() && self.b[x].is_empty()).max_by_key(|&x| self.slot[x]);
            if let Some(x) = conv {
                convslot = Some(self.slot[x]);
                // Only a column some tail column refers to joins the nucleus list; the others become singletons.
                let held = (0..m).any(|k| self.kind[k] == -2 && self.a[k].contains(&x));
                if held {
                    self.kind[x] = self.nuc.len() as i64;
                    self.nuc.push(x);
                } else {
                    self.kind[x] = -1;
                }
                for k in 0..m {
                    if self.kind[k] == -2 && self.a[k].contains(&x) {
                        self.a[k].retain(|&y| y != x);
                        self.b[k].insert(0, x);
                    }
                }
            }
        }
        if old >= 0 {
            self.drop_nuc(r);
        } else if old == -2 {
            let old_b = std::mem::take(&mut self.b[r]);
            self.a[r].clear();
            self.uc[r].clear();
            self.kind[r] = -1;
            let mut refd = vec![false; m];
            for k in 0..m {
                if self.kind[k] == -2 {
                    for &y in &self.b[k] {
                        refd[y] = true;
                    }
                }
            }
            for x in old_b {
                if !refd[x] && self.kind[x] >= 0 {
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
            self.i21 = bt;
            if !av.is_empty() {
                self.walk(r, bt, top as usize);
            }
        }
        let _ = &mut tail_rebuild;
        self.kind[r] = -2;
        let keep = (old != -2 && av.is_empty()) || (old == -2 && !moved);
        // Row r of U: the tail columns that reference r (`0xf2a51`), unless the column stays in its slot.
        let refs = (0..m).filter(|&x| x != r && self.kind[x] == -2 && (self.a[x].contains(&r) || self.b[x].contains(&r))).count();
        if refs > 0 && !(keep && !moved) {
            self.upd_work.push((0xf2a51, 3.0 * refs as f64, crate::work::CTX));
        }
        let mut elim: Vec<usize> = Vec::new();
        for k in 0..m {
            if k == r || self.kind[k] != -2 {
                continue;
            }
            let in_a = self.a[k].contains(&r);
            let in_b = self.b[k].contains(&r);
            if in_a || in_b {
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
                    elim.push(k);
                }
            }
        }
        let mut newdiag = spk[r];
        if !elim.is_empty() && !eta_cand.is_empty() {
            let mut ind = vec![r];
            let mut val = vec![1.0 / d0];
            for &(k, c) in &eta_cand {
                ind.push(k);
                val.push(c);
            }
            let xs: Vec<f64> = ind.iter().map(|&j| spk[j]).collect();
            // The diagonal is summed from the scratch eta, whose alignment is fixed per factor object.
            newdiag = dot(&xs, &val, self.scratch_peel);
            self.sum_eta += ind.len() as f64;
            self.upd_work.push((0xf1beb, 3.0 * ind.len() as f64, crate::work::CTX));
            let e = self.etas.len();
            for &j in &ind {
                self.first_eta[j] = self.first_eta[j].min(e);
            }
            let esz = ind.len();
            self.etas.push(Eta { idx: r, ind, val, peel });
            self.fill += esz;
        }
        if k > 0 {
            self.upd_work.push((0xf2e3c, 5.0 * k as f64, crate::work::CTX));
        }
        // The new U column goes to the pool unless it fits the old column's slot, and each of its entries to a row of the row-wise copy, which moves to the end of the
        // pool with room to grow when it is full.
        if k > old_rows.len() {
            self.fill += k;
        }
        if k > 0 {
            let mut cnt = vec![0usize; m];
            for c in (0..m).filter(|&c| c != r && self.kind[c] == -2) {
                for &(row, _) in &self.uc[c] {
                    cnt[row] += 1;
                }
            }
            // A row that already holds an entry of the old column takes the new entry in place.
            for &(row, _) in self.spike.iter().filter(|e| e.0 != r && !old_rows.contains(&e.0)) {
                if cnt[row] >= self.row_cap[row] {
                    let extra = (cnt[row] / 2).max(4);
                    self.fill += cnt[row] + extra;
                    self.row_cap[row] = cnt[row] + extra;
                }
            }
        }
        self.diag[r] = newdiag;
        self.uc[r] = self.spike.iter().copied().filter(|e| e.0 != r).collect();
        self.a[r] = av;
        let mut bb = bv;
        bb.reverse();
        self.b[r] = bb;
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

