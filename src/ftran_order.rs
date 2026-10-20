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
}

struct Eta {
    idx: usize,
    ind: Vec<usize>,
    val: Vec<f64>,
}

const TOL: f64 = 1e-13;

/// Eta dot product as the reference evaluates it: sequential below eight terms, otherwise one element peeled (its value arrays sit
/// 8 bytes off a 16-byte boundary) and eight lane-wise partial sums.
fn dot(x: &[f64], c: &[f64]) -> f64 {
    let n = c.len();
    let t = |i: usize| x[i] * c[i];
    if n < 9 {
        let mut acc = 0.0;
        for i in 0..n {
            acc += t(i);
        }
        return acc;
    }
    let n8 = n - ((n - 1) & 7);
    let mut lanes = [0.0f64; 8];
    lanes[0] = t(0);
    for i in 1..n8 {
        lanes[(i - 1) & 7] += t(i);
    }
    let mut tot = (((lanes[0] + lanes[2]) + lanes[4]) + lanes[6]) + (((lanes[1] + lanes[3]) + lanes[5]) + lanes[7]);
    for i in n8..n {
        tot += t(i);
    }
    tot
}

impl RefStruct {
    /// Fresh factor from the pivot sequence; `cols[pos]` is the basis column at incoming position `pos` as `(is_slack, rows, values)`.
    pub fn fresh(m: usize, st: &Structure, cols: &[(bool, Vec<usize>, Vec<f64>)]) -> RefStruct {
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
        let perm = step;
        let znz = (nnz + sa + sb + 2 * m) as f64;
        RefStruct {
            upd: 0,
            sum_alpha: 0.0,
            sum_eta: 0.0,
            last_alpha: 0,
            znz,
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
    fn lstage_num(&self, rows: &[usize], vals: &[f64]) -> (Vec<(usize, f64)>, usize) {
        let mut x = vec![0.0f64; self.m];
        let mut start: Vec<usize> = Vec::with_capacity(rows.len());
        for (&r, &v) in rows.iter().zip(vals) {
            let s = self.perm[r];
            x[s] += v;
            start.push(s);
        }
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
        let mut order = q;
        let mut have = vec![false; self.m];
        for &v in &order {
            have[v] = true;
        }
        for eta in &self.etas {
            let xs: Vec<f64> = eta.ind.iter().map(|&j| x[j]).collect();
            let mut s = dot(&xs, &eta.val);
            if s.abs() <= TOL {
                s = 0.0;
            }
            if !have[eta.idx] && s != 0.0 {
                order.push(eta.idx);
                have[eta.idx] = true;
            }
            x[eta.idx] = s;
        }
        (order.into_iter().filter(|&v| x[v].abs() > TOL).map(|v| (v, x[v])).collect(), count_l)
    }

    /// FTRAN of a sparse column: the nonzeros `(position, value)` of `B^-1 a` in the reference's output order. The spike is kept
    /// for the following `update`.
    pub fn ftran(&mut self, rows: &[usize], vals: &[f64]) -> Vec<(usize, f64)> {
        let m = self.m;
        let (spike, count_l) = self.lstage_num(rows, vals);
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
        let ord = if dense { (0..m).collect() } else { self.uorder(&inp) };
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
        self.spike = spike;
        self.last_alpha = out.len();
        out
    }

    /// The reference's refactorization decision, evaluated before an iteration's own update: at the update limit
    /// `max(20, m^0.7)`, or from `U >= m^0.45` when the accumulated work `(1.5 S + 2 c) m` exceeds `min((Z - 1.5 m)^2, 4 m^2 U0)`.
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

    /// Order of the nonzero positions of the solution, given the L-stage list.
    pub fn uorder(&self, inp: &[usize]) -> Vec<usize> {
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
            lb += self.b[k].len();
            for &j in &self.a[k] {
                if !have[j] {
                    have[j] = true;
                    heap.push((self.slot[j], j));
                }
            }
        }
        let mut nm = vec![false; m];
        let mut lst: Vec<usize> = Vec::new();
        for &x in nuc.iter().chain(order.iter().flat_map(|&k| self.b[k].iter())) {
            if !nm[x] {
                nm[x] = true;
                lst.push(x);
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
        let _ = r;
        while cur != hole {
            let pos = self.at[cur];
            if self.a[pos].len() + self.b[pos].len() >= cur - hole {
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
        self.upd += 1;
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
            let last = sp.pop().unwrap();
            if i < sp.len() {
                sp[i] = last;
            }
        }
        let old = self.kind[r];
        let mut convslot: Option<usize> = None;
        let mut conv: Option<usize> = None;
        if old != -2 {
            conv = (0..m).filter(|&x| self.kind[x] == -2 && self.a[x].is_empty() && self.b[x].is_empty()).max_by_key(|&x| self.slot[x]);
            if let Some(x) = conv {
                convslot = Some(self.slot[x]);
                self.kind[x] = self.nuc.len() as i64;
                self.nuc.push(x);
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
                    self.b[k].retain(|&y| y != r);
                } else {
                    self.a[k].retain(|&y| y != r);
                    self.b[k].retain(|&y| y != r);
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
            newdiag = dot(&xs, &val);
            self.sum_eta += ind.len() as f64;
            self.etas.push(Eta { idx: r, ind, val });
        }
        self.diag[r] = newdiag;
        self.uc[r] = self.spike.iter().copied().filter(|e| e.0 != r).collect();
        self.a[r] = av;
        let mut bb = bv;
        bb.reverse();
        self.b[r] = bb;
    }
}
