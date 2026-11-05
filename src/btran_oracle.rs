use super::*;
#[rustfmt::skip]
impl<R: Real> RefStruct<R> {
    fn frozen_btran(&mut self, inp: &[(usize, R)], direct: bool, dense_in: bool) -> Vec<(usize, R)> {
        let m = self.m;
        let ut = &self.ut;
        let mut x = vec![R::ZERO; m];
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
                    let mut inl = vec![false; m];
                    for &p in l.iter() {
                        inl[p] = true;
                    }
                    let mut es = netas as f64 + netas as f64 + c0 as f64 * 5.0;
                    for eta in self.etas.iter().rev() {
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
                        if x[p].abs() >= R::c(TOL) {
                            true
                        } else {
                            x[p] = R::ZERO;
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
                    list = None;
                }
            }
        }
        // L^T: row i of L pushes into the earlier nodes.
        let lt = &self.lt;
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
                e = el + e;
            } else {
                // Level buckets of the transposed graph when the factor has levels (`0x12bd30`), else a max-heap on the node index; every
                // reached node is kept when at least the tolerance.
                let mut seen = vec![false; m];
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
                            if v.abs() >= R::c(TOL) {
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
                                x[i] = R::ZERO;
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
                        if v.abs() >= R::c(TOL) {
                            res.push(i);
                            for &(j, lv) in &lt[i] {
                                x[j] -= lv * v;
                                if !seen[j] {
                                    seen[j] = true;
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
        if self.bt_dense {
            // The dense result keeps the sub-tolerance entries the L^T pass left in place.
            // and lists them by ascending row.
            if R::LANES {
                return Self::dense_btran_output(&x, &self.perm[..m]);
            }
            let mut out = Vec::with_capacity(m);
            for r in 0..m {
                let i = self.perm[r];
                if x[i] != R::ZERO {
                    out.push((r, x[i]));
                }
            }
            return out;
        }
        res.into_iter()
            .filter(|&i| x[i] != R::ZERO)
            .map(|i| (row_of[i], x[i]))
            .collect()
    }

}

#[test]
fn whole_workspace_matches_frozen_btran_through_route_changes() {
    fn check<R: Real>(bits: impl Fn(R) -> (u64, u16)) {
        let m = 512;
        let cols: Vec<_> = (0..m).map(|i| (true, vec![i], vec![1.0])).collect();
        let st = crate::lu_order::factor_structure_r::<R>(m, &cols);
        let mut f = RefStruct::fresh(m, m, &st, &cols);
        f.h33 = 0.0;
        f.h34 = 0.0;
        f.i21 = m - 2;
        f.kind[m - 2] = -2;
        f.kind[m - 1] = -2;
        f.ut[509] = vec![(510, R::c(0.5))];
        f.ut[499] = vec![(511, R::c(-0.25)), (510, R::ONE), (511, R::c(0.25))];
        f.ut[510] = vec![(511, R::c(0.5))];
        f.lstart.clear();
        f.lentries.clear();
        for i in 0..m {
            f.lstart.push(f.lentries.len());
            if i % 8 != 7 || i == 487 {
                for &v in &[R::c(-0.5), R::c(0.125)] {
                    f.lentries.push((i + 1, v));
                    f.lt[i + 1].push((i, v));
                }
            }
        }
        f.lstart.push(f.lentries.len());
        f.nnz_l = f.lentries.len();
        f.dthr = 2;
        f.lev_dense = false;
        for e in 0..64 {
            let r = [509, 507, 505, 503, 501, 499, 497, 495][e % 8];
            f.etas.push(Eta {
                idx: r,
                ind: vec![r, r - 2, r - 2],
                val: vec![R::ONE, R::c(0.25), R::c(-0.125)],
                peel: 0,
            });
            f.sum_eta += 3.0;
        }
        let mut old = f.clone();
        let cases = [
            vec![(509, R::ONE)],
            vec![(11, R::c(-0.0))],
            vec![],
            vec![(507, R::ONE), (507, R::c(-1.0)), (499, R::c(f64::MAX))],
            vec![
                (495, R::c(f64::from_bits(0x7ff8000000000123))),
                (510, R::ONE),
            ],
            vec![(508, R::c(TOL)), (11, R::c(-TOL))],
            vec![(488, R::ONE), (499, R::c(f64::INFINITY))],
        ];
        for pass in 0..30 {
            if pass == 12 {
                f = f.clone();
                old = old.clone();
            }
            if pass == 20 {
                f.bt_work.as_mut().unwrap().generation = u32::MAX;
            }
            for factor in [&mut f, &mut old] {
                factor.loff = if pass % 3 == 0 { m - 24 } else { m };
                factor.nnz_l_tail = if pass % 3 == 0 { 384 } else { 0 };
                factor.n1 = 0;
                factor.lev_dense = pass % 7 == 1;
                factor.levels = if pass % 2 == 1 {
                    Some((vec![0; m], 0))
                } else {
                    None
                };
                if pass == 13 {
                    factor.etas.push(Eta {
                        idx: 509,
                        ind: vec![509, 497],
                        val: vec![R::c(0.5), R::ONE],
                        peel: 0,
                    });
                    factor.sum_eta += 2.0;
                }
            }
            for (i, rhs) in cases.iter().enumerate() {
                let dense = (pass + i) % 11 == 1;
                let direct = i % 2 == 0;
                let aa = old.frozen_btran(rhs, direct, dense);
                let bb = f.btran(rhs, direct, dense);
                assert_eq!(
                    aa.iter().map(|&(i, v)| (i, bits(v))).collect::<Vec<_>>(),
                    bb.iter().map(|&(i, v)| (i, bits(v))).collect::<Vec<_>>(),
                    "pass {pass} case {i}"
                );
                assert_eq!(f.bt_dense, old.bt_dense);
                assert_eq!(f.work.to_bits(), old.work.to_bits());
                assert_eq!(f.row_cached, old.row_cached);
                assert_eq!(f.row_of.get(), old.row_of.get());
                assert!(f
                    .bt_work
                    .as_ref()
                    .unwrap()
                    .x
                    .iter()
                    .all(|&v| bits(v) == bits(R::ZERO)));
            }
        }
    }
    check::<f64>(|v| (v.to_bits(), 0));
    check::<crate::f80::F80>(|v| v.repr());
}
