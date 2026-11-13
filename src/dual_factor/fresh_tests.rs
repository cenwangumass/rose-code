//! Differential test: the flat `DualFactor::fresh` against the construction it replaced (kept verbatim here), every field bitwise.
use super::*;

impl<R: Real> DualFactor<R> {
    /// The construction as it was before the flat build (differential oracle for `fresh`); `cols[pos]` is the basis column at incoming position `pos` as `(is_slack, rows, values)`.
    pub(super) fn fresh_old<'a>(
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
        let n1 = st.phase.iter().filter(|&&p| p == 1 || p == 2).count();
        let n3 = st.phase.iter().filter(|&&p| p == 3).count();
        let nseq = st.loff;
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
        let mut lpool = LPool::with_capacity(m, st.l_nnz());
        for s in 0..st.steps() {
            lpool.push_line(st.l_rows(s).iter().zip(st.l_vals(s)).map(|(&r, &v)| (step[r], v)));
        }
        let lcolumn = |v: usize| lpool.iter(v);
        let lt = lpool.transpose(m);
        let l_forward = (0..lpool.lines()).all(|v| lpool.col(v).0.iter().all(|&y| y as usize > v));
        // pivot-row values of the Markowitz steps, by column
        let mut mk: Vec<Vec<(usize, R)>> = vec![Vec::new(); m];
        for s in nseq..m {
            for &(c, v) in st.u_row(s) {
                mk[c].push((s, v));
            }
        }
        let mut a = vec![Vec::new(); m];
        let mut b = vec![Vec::new(); m];
        let mut uc: Vec<Vec<(usize, R)>> = vec![Vec::new(); m];
        let (mut sa, mut sb) = (0usize, 0usize);
        let mut xv = vec![R::ZERO; m];
        let mut reached = vec![false; m];
        let mut epos = vec![usize::MAX; m];
        let mut changed = vec![false; m];
        for t in i21..m {
            let c = st.pivots[t].1;
            let (_, idx, val) = &cols[c];
            let mut start: Vec<usize> = Vec::new();
            for (&i, &v) in idx.iter().zip(val) {
                let s = step[i];
                if s <= t {
                    xv[s] += R::c(v);
                    start.push(s);
                }
            }
            // Only the reached set is needed: the numeric solve below visits it in ascending order.
            let mut lst = Vec::new();
            for &v in &start {
                if !reached[v] {
                    reached[v] = true;
                    lst.push(v);
                }
            }
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
            let mut q = lst.clone();
            q.sort_unstable();
            for &v in &q {
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
            let mut ents: Vec<(usize, R)> = Vec::new();
            let touched = q;
            for &v in &touched {
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
            for &(sidx, vv) in &mk[c] {
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
            for &v in &lst {
                xv[v] = R::ZERO;
                reached[v] = false;
            }
            for &v in &touched {
                changed[v] = false;
            }
            for &v in &touched {
                xv[v] = R::ZERO;
            }
            ents.sort_unstable_by_key(|e| e.0);
            let av: Vec<usize> = ents
                .iter()
                .map(|e| e.0)
                .filter(|&x| kind[x] == -2)
                .collect();
            let mut bv: Vec<usize> = ents.iter().map(|e| e.0).filter(|&x| kind[x] >= 0).collect();
            bv.sort_unstable_by_key(|&x| std::cmp::Reverse(kind[x]));
            sa += av.len();
            sb += bv.len();
            // The reference lists the fresh tail column's A entries in descending order.
            a[t] = av.into_iter().rev().collect();
            b[t] = bv;
            uc[t] = ents;
        }
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
        // The lists as the old construction built them, with the coefficients looked up in `uc` (an independent derivation).
        let mut old_ab = AbLists::new();
        old_ab.reset(m);
        for t in 0..m {
            let with_coeff = |l: &[usize]| -> Vec<(u32, R)> {
                l.iter().map(|&x| (x as u32, uc[t].iter().find(|e| e.0 == x).map_or(R::ZERO, |e| e.1))).collect()
            };
            old_ab.set(t, &with_coeff(&a[t]), &with_coeff(&b[t]));
        }
        let mark0 = (i21..m)
            .rev()
            .find(|&t| kind[t] == -2 && a[t].is_empty() && b[t].is_empty())
            .unwrap_or(i21.min(m - 1));
        let mut ut = vec![Vec::new(); m];
        for (c, col) in uc.iter().enumerate() {
            for &(row, v) in col {
                ut[row].push((c, v));
            }
        }
        DualFactor {
            upd: 0,
            sum_alpha: 0.0,
            before_last: (0, 0.0),
            sum_eta: 0.0,
            last_alpha: 0,
            znz,
            n1: st.phase.iter().filter(|&&p| p == 1).count(),
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
            fill: 2 * (uc.iter().map(Vec::len).sum::<usize>() + nnz_l_tot),
            // A row lists its tail columns in the order of the factorization's U row (Markowitz steps), or in the order of the
            // basis columns for a row of the singleton phases.
            urows: URows::from_nested(&{
                let mut member = vec![Vec::new(); m];
                for t in i21..m {
                    if kind[t] == -2 {
                        for &y in a[t].iter().chain(b[t].iter()) {
                            member[y].push(t);
                        }
                    }
                }
                let mut pos_of_col = vec![usize::MAX; m];
                for (t, &(_, c)) in st.pivots.iter().enumerate() {
                    pos_of_col[c] = t;
                }
                let mut rl = vec![Vec::new(); m];
                let mut stamp = vec![usize::MAX; m];
                for y in 0..m {
                    if member[y].is_empty() {
                        continue;
                    }
                    let order: Vec<usize> = if let Some(u) = st.pair_urow.get(&y) {
                        u.iter().map(|&c| pos_of_col[c]).collect()
                    } else if st.u_row(y).is_empty() {
                        // Filtering all basis columns in index order is equivalent to sorting just this row's members by
                        // column. Keep only the positions represented by pos_of_col, and emit each position once.
                        let mut row: Vec<usize> = member[y]
                            .iter()
                            .copied()
                            .filter(|&t| pos_of_col[st.pivots[t].1] == t)
                            .collect();
                        row.sort_unstable_by_key(|&t| st.pivots[t].1);
                        row.dedup();
                        rl[y] = row;
                        continue;
                    } else {
                        st.u_row(y).iter().map(|&(c, _)| pos_of_col[c]).collect()
                    };
                    for &t in &member[y] {
                        stamp[t] = y;
                    }
                    rl[y] = order.into_iter().filter(|&t| stamp[t] == y).collect();
                }
                rl
            }, &ut),
            row_cap: {
                let mut c = vec![0usize; m];
                for col in &uc {
                    for &(row, _) in col {
                        c[row] += 1;
                    }
                }
                c
            },
            // A fresh factor lays the rows of the row-wise copy out contiguously: the nonempty rows of the nucleus in the order of their list index
            // (`kind`), then the others by position, each row as long as its U entries (`ucol_len`); the parity of the start decides the alignment of the row's
            // value array in the dot product of the second FTRAN vector.
            row_par: {
                let mut len = vec![0usize; m];
                for col in &uc {
                    for &(row, _) in col {
                        len[row] += 1;
                    }
                }
                let mut order: Vec<usize> =
                    (0..m).filter(|&p| kind[p] >= 0 && len[p] > 0).collect();
                order.sort_by_key(|&p| kind[p]);
                order.extend((0..m).filter(|&p| kind[p] < 0 && len[p] > 0));
                let mut par = vec![0u8; m];
                let mut acc = 0usize;
                for p in order {
                    par[p] = (acc % 2) as u8;
                    acc += len[p];
                }
                par
            },
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
            empty_tail: (0..m)
                .filter(|&x| kind[x] == -2 && a[x].is_empty() && b[x].is_empty())
                .collect(),
            kind,
            nuc,
            i21,
            mark: mark0,
            ab: old_ab,
            lpool,
            l_forward,
            perm_bound: perm.iter().map(|&p| p + 1).max().unwrap_or(0),
            perm,
            row_of: Default::default(),
            dens: (h33 + h34) * (h33 + h34),
            dense_est: std::cell::Cell::new(0.0),
            dthr,
            slot: (0..m).collect(),
            at: (0..m).collect(),
            btran_ws: Default::default(),
            etas: Vec::new(),
            lt_distinct: lt_lines_distinct(&lt),
            lt,
            levels,
            level_cache: LevelCache::default(),
            level_spare: LevelCache::default(),
            lists: Default::default(),
            lscr: Default::default(),
            lt_levels: Default::default(),
            // Four or more occupancy words, and skipped zero divisions must always be exact.
            bt_occupied: !R::QUAD && m >= 256 && st.diag.iter().all(|d| d.to_f64().is_normal()),
            bt_work: None,
            lev_built,
            diag: st.diag.clone(),
            bad_diag: st.diag.iter().any(|d| !(d.abs() > R::ZERO)),
            uc,
            spike: Vec::new(),
            update_values: Vec::new(),
            update_seen: Vec::new(),
            spike_used: false,
            lsnap: Default::default(),
        }
    }
}

/// The bytes of a number; the 16-byte x87 extended value has 6 padding bytes that are not part of it.
fn raw<T: Copy>(v: &T) -> &[u8] {
    let n = std::mem::size_of::<T>();
    unsafe { std::slice::from_raw_parts(v as *const T as *const u8, if n == 16 { 10 } else { n }) }
}

fn same_vals<R: Real>(x: &[R], y: &[R], what: &str) {
    assert_eq!(x.len(), y.len(), "{what} length");
    for (i, (p, q)) in x.iter().zip(y).enumerate() {
        assert_eq!(raw(p), raw(q), "{what}[{i}]");
    }
}

fn same_pairs<R: Real>(x: &[Vec<(usize, R)>], y: &[Vec<(usize, R)>], what: &str) {
    assert_eq!(x.len(), y.len(), "{what} length");
    for (i, (p, q)) in x.iter().zip(y).enumerate() {
        assert_eq!(p.len(), q.len(), "{what}[{i}] length");
        for (k, (e, f)) in p.iter().zip(q).enumerate() {
            assert_eq!(e.0, f.0, "{what}[{i}][{k}] index");
            assert_eq!(raw(&e.1), raw(&f.1), "{what}[{i}][{k}] value");
        }
    }
}

fn same_csr<R: Real>(x: &Csr<R>, y: &Csr<R>, what: &str) {
    assert_eq!(x.start, y.start, "{what} start");
    assert_eq!(x.idx, y.idx, "{what} idx");
    same_vals(&x.val, &y.val, what);
}

/// Lists, order and coefficient bits of every column.
pub(super) fn assert_same_ab<R: Real>(n: &AbLists<R>, o: &AbLists<R>) {
    assert_eq!(n.to_vecs(), o.to_vecs(), "A/B lists");
    for c in 0..n.to_vecs().0.len() {
        same_vals(n.a_val(c), o.a_val(c), "A coefficients");
        same_vals(n.b_val(c), o.b_val(c), "B coefficients");
    }
}

/// Every field a fresh factor sets, bitwise (the rest is the default state, equal by construction).
pub(super) fn assert_same<R: Real>(n: &DualFactor<R>, o: &DualFactor<R>) {
    assert_eq!((n.m, n.i21, n.mark, n.loff, n.n1), (o.m, o.i21, o.mark, o.loff, o.n1));
    assert_eq!((n.kind.clone(), n.nuc.clone()), (o.kind.clone(), o.nuc.clone()));
    assert_eq!((n.urows.to_rowl(), &n.row_cap, &n.row_par), (o.urows.to_rowl(), &o.row_cap, &o.row_par));
    assert_same_ab(&n.ab, &o.ab);
    assert_eq!((&n.perm, &n.first_eta, &n.slot, &n.at), (&o.perm, &o.first_eta, &o.slot, &o.at));
    assert_eq!((n.empty_tail.iter().copied().collect::<Vec<usize>>(), n.l_forward, n.dthr, n.lev_dense, n.lev_dense_ft, n.lev_built), (
        o.empty_tail.iter().copied().collect::<Vec<usize>>(),
        o.l_forward,
        o.dthr,
        o.lev_dense,
        o.lev_dense_ft,
        o.lev_built
    ));
    assert_eq!((n.nnz_l, n.nnz_l_tail, n.fill, n.scratch_peel, n.bt_occupied), (o.nnz_l, o.nnz_l_tail, o.fill, o.scratch_peel, o.bt_occupied));
    for (x, y, w) in [(n.znz, o.znz, "znz"), (n.h33, o.h33, "h33"), (n.h34, o.h34, "h34"), (n.dens, o.dens, "dens")] {
        assert_eq!(x.to_bits(), y.to_bits(), "{w}");
    }
    assert_eq!(n.levels, o.levels);
    same_csr(&n.lpool, &o.lpool, "lpool");
    same_csr(&n.lt, &o.lt, "lt");
    same_vals(&n.diag, &o.diag, "diag");
    same_pairs(&n.uc, &o.uc, "uc");
    same_pairs(&n.urows.to_ut(), &o.urows.to_ut(), "ut");
    // the physical order carries the same values entry for entry
    for y in 0..n.m {
        same_vals(n.urows.row(y).1, o.urows.row(y).1, "urows values");
    }
    n.urows.check(&n.uc, &n.kind);
}

/// Leave junk in every row, as the edits of a long update sequence would.
fn dirty(f: &mut DualFactor<f64>) {
    for i in 0..f.m {
        f.ab.set(i, &[(7, -1.0), (9, 2.0), (11, 0.5)], &[(5, 3.0)]);
        f.uc[i].push((3, -1.5));
        for c in 0..40 {
            f.urows.push(i, c, 2.5);
        }
    }
}

fn check(m: usize, cols: &[(bool, Vec<usize>, Vec<f64>)]) {
    let st = crate::lu_order::factor_structure(m, cols);
    // The builds after the first run on the retained scratch and on the lists of the previous factor, which has a different size
    // and content (`recycle`); the last one runs on lists that a column replacement has edited.
    for round in 0..3 {
        let mut n = DualFactor::<f64>::fresh(m, m, &st, cols);
        let o = DualFactor::<f64>::fresh_old(m, m, &st, cols);
        assert_same(&n, &o);
        if round == 1 {
            dirty(&mut n);
        }
        n.recycle();
    }
    let st = crate::lu_order::factor_structure_r::<crate::f80::F80>(m, cols);
    let n = DualFactor::<crate::f80::F80>::fresh(m, m, &st, cols);
    let o = DualFactor::<crate::f80::F80>::fresh_old(m, m, &st, cols);
    assert_same(&n, &o);
    n.recycle();
}

#[test]
fn flat_fresh_matches_the_construction_it_replaced() {
    let mut seed = 0x9e3779b97f4a7c15u64;
    let mut next = || {
        seed ^= seed << 13;
        seed ^= seed >> 7;
        seed ^= seed << 17;
        seed as usize
    };
    // sizes alternate so the retained scratch grows, shrinks and is reused; densities from singleton-only to heavy nucleus
    for (round, &m) in [3usize, 40, 17, 321, 65, 300, 1, 120, 300, 512].iter().enumerate() {
        for dens in [0usize, 1, 2, 4, 7] {
            let cols: Vec<_> = (0..m)
                .map(|i| {
                    let mut rows = vec![i];
                    let mut vals = vec![2.0 + (next() % 3) as f64];
                    for _ in 0..dens {
                        let r = next() % m;
                        if !rows.contains(&r) {
                            rows.push(r);
                            vals.push([0.125, -0.25, 1.0, -3.0, 1e-14][next() % 5]);
                        }
                    }
                    if (next() + round) % 4 == 0 {
                        // shuffled positions: the diagonal is not the pivot
                        rows.reverse();
                        vals.reverse();
                    }
                    (dens == 0 || vals.len() == 1, rows, vals)
                })
                .collect();
            check(m, &cols);
        }
    }
}

/// Columns in a scrambled order (the basis column of a tail position is neither ascending in the position nor in the row): the
/// singleton-phase rows list their members in ascending basis column, which the fill order must give.
#[test]
fn scrambled_basis_columns_match_the_construction_it_replaced() {
    let mut seed = 0x2545f4914f6cdd1du64;
    let mut next = || {
        seed ^= seed << 13;
        seed ^= seed >> 7;
        seed ^= seed << 17;
        seed as usize
    };
    for &m in &[5usize, 23, 64, 150, 301] {
        for dens in [1usize, 2, 3, 6] {
            let mut cols: Vec<(bool, Vec<usize>, Vec<f64>)> = (0..m)
                .map(|i| {
                    let mut rows = vec![i];
                    let mut vals = vec![2.0 + (next() % 3) as f64];
                    for _ in 0..dens {
                        let r = next() % m;
                        if !rows.contains(&r) {
                            rows.push(r);
                            vals.push([0.125, -0.25, 1.0, -3.0, 1e-14][next() % 5]);
                        }
                    }
                    (vals.len() == 1, rows, vals)
                })
                .collect();
            for i in (1..m).rev() {
                cols.swap(i, next() % (i + 1));
            }
            check(m, &cols);
        }
    }
}
