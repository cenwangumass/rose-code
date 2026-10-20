//! Structural mirror of the reference sparse FTRAN: the order in which it lists the nonzeros of `B^-1 a` (`reverse-engineer/specs/ftran-order.md`).
//!
//! Positions are the basis positions after the LU renumbering (the pivot step). A fresh factor has singleton positions (kind -1),
//! a nucleus list of diagonal-only positions (phase 3 pivots, kind = index in the list) and a triangular tail (phase 4 and the
//! Markowitz nucleus, slots `i21..m`) whose columns carry U entries split into tail positions (`a`) and nucleus-list positions (`b`).
//! Column replacements follow the update rules recovered from the GDB dumps. Only the structure is mirrored; values come from the
//! solver's own factorization, entries that cancel numerically are dropped by the caller.

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
    slot: Vec<usize>,
    at: Vec<usize>,
    etas: Vec<(usize, Vec<usize>)>,
}

impl RefStruct {
    /// Fresh structure from the pivot sequence; `cols[pos]` are the rows of the basis column at incoming position `pos`.
    pub fn fresh(m: usize, st: &Structure, cols: &[Vec<usize>]) -> RefStruct {
        let mut step = vec![usize::MAX; m];
        for (s, &(r, _)) in st.pivots.iter().enumerate() {
            step[r] = s;
        }
        let n1 = st.phase.iter().filter(|&&p| p == 1 || p == 2).count();
        let n3 = st.phase.iter().filter(|&&p| p == 3).count();
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
        let mut a = vec![Vec::new(); m];
        let mut b = vec![Vec::new(); m];
        let (mut sa, mut sb) = (0usize, 0usize);
        let mut mark = vec![usize::MAX; m];
        for t in i21..m {
            let c = st.pivots[t].1;
            let mut seen: Vec<usize> = Vec::new();
            for &r in &cols[c] {
                let x = step[r];
                if x < t && mark[x] != t {
                    mark[x] = t;
                    seen.push(x);
                }
            }
            let mut k = 0;
            while k < seen.len() {
                let v = seen[k];
                k += 1;
                for &x in &lcols[v] {
                    if x < t && mark[x] != t {
                        mark[x] = t;
                        seen.push(x);
                    }
                }
            }
            let mut av: Vec<usize> = seen.iter().copied().filter(|&x| kind[x] == -2).collect();
            av.sort_unstable();
            let mut bv: Vec<usize> = seen.iter().copied().filter(|&x| kind[x] >= 0).collect();
            bv.sort_unstable_by_key(|&x| std::cmp::Reverse(kind[x]));
            sa += av.len();
            sb += bv.len();
            a[t] = av;
            b[t] = bv;
        }
        let nt = m - i21;
        let nn = nuc.len();
        let h33 = if nt > 0 { (sa + nt) as f64 / nt as f64 } else { 1.5 };
        let h34 = if nn > 0 { (sb + nn) as f64 / nn as f64 } else { 1.5 };
        let perm = step;
        RefStruct {
            m,
            kind,
            nuc,
            i21,
            a,
            b,
            lcols,
            perm,
            dens: (h33 + h34) * (h33 + h34),
            slot: (0..m).collect(),
            at: (0..m).collect(),
            etas: Vec::new(),
        }
    }

    /// L solve in the reference's order (Kahn FIFO over the L columns) followed by the row etas: the list handed to the U solve.
    pub fn lstage(&self, rows: &[usize]) -> Vec<usize> {
        let mut cnt = vec![0i32; self.m];
        let mut lst: Vec<usize> = rows.iter().map(|&r| self.perm[r]).collect();
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
        let mut have = vec![false; self.m];
        for &v in &q {
            have[v] = true;
        }
        for (idx, ind) in &self.etas {
            if !have[*idx] && ind.iter().any(|&j| have[j]) {
                q.push(*idx);
                have[*idx] = true;
            }
        }
        q
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

    pub fn ftran_order(&self, rows: &[usize]) -> Vec<usize> {
        self.uorder(&self.lstage(rows))
    }

    fn drop_nuc(&mut self, x: usize) {
        let idx = self.kind[x] as usize;
        let last = self.nuc.pop().unwrap();
        if last != x {
            self.nuc[idx] = last;
            self.kind[last] = idx as i64;
        }
    }

    /// Replace the column at position `r` by the entering column whose L-stage list is `spike`.
    pub fn update(&mut self, r: usize, spike: &[usize]) {
        let m = self.m;
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
        if old != -2 {
            let conv = (0..m).filter(|&x| self.kind[x] == -2 && self.a[x].is_empty() && self.b[x].is_empty()).max_by_key(|&x| self.slot[x]);
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
        let mut tail: Vec<usize> = (0..m).filter(|&p| self.kind[p] == -2 && p != r).collect();
        tail.sort_unstable_by_key(|&p| self.slot[p]);
        let top: i64 = av.iter().map(|&x| self.slot[x] as i64).max().unwrap_or(-1);
        let mut moved = true;
        let idx = if old == -2 {
            if top < oldslot as i64 {
                moved = false;
                tail.iter().filter(|&&p| self.slot[p] < oldslot).count()
            } else {
                tail.iter().filter(|&&p| (self.slot[p] as i64) <= top).count()
            }
        } else if matches!(convslot, Some(cs) if top < cs as i64) {
            let cs = convslot.unwrap();
            tail.iter().filter(|&&p| self.slot[p] < cs).count()
        } else if !av.is_empty() {
            tail.iter().filter(|&&p| (self.slot[p] as i64) <= top).count()
        } else {
            0
        };
        tail.insert(idx, r);
        let i21n = m - tail.len();
        let mut in_tail = vec![false; m];
        for &p in &tail {
            in_tail[p] = true;
        }
        let mut taken = vec![false; i21n];
        let mut misplaced: Vec<usize> = Vec::new();
        for p in 0..m {
            if in_tail[p] {
                continue;
            }
            if self.slot[p] < i21n {
                taken[self.slot[p]] = true;
            } else {
                misplaced.push(p);
            }
        }
        misplaced.sort_unstable_by_key(|&p| self.slot[p]);
        let holes: Vec<usize> = (0..i21n).filter(|&t| !taken[t]).collect();
        for (&p, &t) in misplaced.iter().zip(holes.iter()) {
            self.slot[p] = t;
            self.at[t] = p;
        }
        self.i21 = i21n;
        for (s, &p) in tail.iter().enumerate() {
            self.slot[p] = i21n + s;
            self.at[i21n + s] = p;
        }
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
                    elim.push(k);
                }
            }
        }
        if !elim.is_empty() {
            let mut ind = vec![r];
            ind.extend(elim);
            self.etas.push((r, ind));
        }
        self.a[r] = av;
        let mut bb = bv;
        bb.reverse();
        self.b[r] = bb;
    }
}
