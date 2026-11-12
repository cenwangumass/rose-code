//! The construction before the single-pass front end (K2a): `row_columns`, the phases with their own `rc` build, the separate singleton
//! L loop, the `nuc_ent` copy and the scans of `Structure::work`. Kept as the oracle of the tests and of the `factor-shadow` build.

use super::*;

pub(super) fn factor_structure_r<'a, R: Real>(m: usize, cols: impl Into<Cols<'a>>) -> Structure<R> {
    let cols: Cols<'a> = cols.into();
    let (mut seq, mut phase, mut ov, loff, pair_gate_record) = phases(m, cols);
    let nsing = seq.len();
    let mut done_r = vec![false; m];
    let mut done_c = vec![false; cols.len()];
    for &(r, c) in &seq {
        done_r[r] = true;
        done_c[c] = true;
    }
    // The nucleus columns, flat: each column's kept entries in reverse input order.
    let mut nuc_pos: Vec<usize> = Vec::new();
    let mut nuc_start: Vec<usize> = vec![0];
    let mut nuc_ent: Vec<(usize, R)> = Vec::new();
    for (pos, (_, idx, val)) in cols.iter().enumerate() {
        if done_c[pos] {
            continue;
        }
        let first = nuc_ent.len();
        nuc_ent.extend(
            idx.iter()
                .copied()
                .zip(val.iter().copied())
                // The reference's matrix leaves out entries below the absolute tolerance (fx_e2_r25: a scaled entry 2^-44 in a row of the nucleus changed its count).
                .filter(|(i, v)| !done_r[*i] && v.abs() >= ABS_TOL)
                .map(|(i, v)| (i, R::c(v))),
        );
        nuc_ent[first..].reverse();
        nuc_pos.push(pos);
        nuc_start.push(nuc_ent.len());
    }
    // The work estimate counts every input entry in active rows and columns
    // before the numeric nucleus applies its magnitude filter (0xd4e4b,
    // 0xdf521..0xdf5ed); scaled entries below ABS_TOL still contribute.
    let nuc_nnz = cols
        .iter()
        .enumerate()
        .filter(|(pos, _)| !done_c[*pos])
        .map(|(_, (slack, rows, _))| {
            if *slack {
                1
            } else {
                rows.iter().filter(|r| !done_r[**r]).count()
            }
        })
        .sum();
    let nnz_cols: usize = cols.iter().map(|c| c.1.len()).sum();
    let nuc = nucleus(m, nnz_cols, &nuc_pos, &nuc_start, &nuc_ent);
    let nuc_work = nuc.work;
    seq.extend(nuc.pivots.iter().copied());
    phase.resize(seq.len(), 0);
    let mut step = vec![usize::MAX; m];
    for (s, &(r, _)) in seq.iter().enumerate() {
        step[r] = s;
    }
    let nsteps = seq.len();
    let mut l_start: Vec<usize> = Vec::with_capacity(nsteps + 1);
    let mut l_idx: Vec<usize> = Vec::with_capacity(nnz_cols);
    let mut l_val: Vec<R> = Vec::with_capacity(nnz_cols);
    let mut l_raw: Vec<R> = Vec::with_capacity(nnz_cols);
    let mut diag: Vec<R> = Vec::with_capacity(nsteps);
    l_start.push(0);
    // Marker words of the singleton steps with an override: the hash is read only for those.
    let mut has_ov = vec![0u64; nsing.div_ceil(64)];
    for &s in ov.l.keys() {
        if s < nsing {
            has_ov[s / 64] |= 1 << (s % 64);
        }
    }
    for (s, &(prow, c)) in seq.iter().enumerate().take(nsing) {
        let (_, idx, val) = &cols[c];
        let pv = R::c(val[idx.iter().position(|&i| i == prow).unwrap()]);
        // A pair override replaces the column (its multipliers are all that is known: no raw entries).
        if let Some((r, v)) = if has_ov[s / 64] >> (s % 64) & 1 == 0 { None } else { ov.l.remove(&s) } {
            l_idx.extend_from_slice(&r);
            l_raw.extend(std::iter::repeat_n(R::c(f64::NAN), v.len()));
            l_val.extend(v);
        } else {
            for (&i, &x) in idx.iter().zip(val) {
                if step[i] > s {
                    l_idx.push(i);
                    l_val.push(R::c(x) / pv);
                    l_raw.push(R::c(x));
                }
            }
        }
        l_start.push(l_idx.len());
        diag.push(pv);
    }
    let pair_unz = ov.unz.clone();
    let pair_urow = ov.urow.clone();
    for (&s, &d) in &ov.diag {
        diag[s] = d;
    }
    let l_base = l_idx.len();
    l_idx.extend_from_slice(&nuc.l_idx);
    l_raw.extend(std::iter::repeat_n(R::c(f64::NAN), nuc.l_val.len()));
    l_val.extend_from_slice(&nuc.l_val);
    l_start.extend(nuc.l_start[1..].iter().map(|&e| l_base + e));
    diag.extend_from_slice(&nuc.diag);
    // Singleton steps have no U row.
    let mut u_start: Vec<usize> = Vec::with_capacity(nsteps + 1);
    u_start.resize(nsing + 1, 0);
    u_start.extend(nuc.u_start[1..].iter().copied());
    let u_ent = nuc.u_ent;
    let mut st = Structure {
        nnz_b: 0,
        u_sing: 0,
        nsing,
        pivots: seq,
        phase,
        l_start,
        l_idx,
        l_val,
        l_raw,
        diag,
        u_start,
        u_ent,
        nuc_nnz,
        nuc_work,
        loff,
        pair_unz,
        pair_urow,
        dense_n: DENSE_N.with(|c| c.get()),
        m,
        pair_gate_record,
    };
    st.nnz_b = cols.iter().map(|c| if c.0 { 1 } else { c.1.len() }).sum();
    st.u_sing = u_sing_scan(&st, cols);
    st
}


/// Row incidence in ascending column order, retaining every occurrence of a duplicate row.
pub(super) fn row_columns(m: usize, cols: Cols<'_>) -> (Vec<usize>, Vec<usize>) {
    let mut offsets = vec![0; m + 1];
    for (_, rows, _) in cols.iter() {
        for &r in rows {
            offsets[r + 1] += 1;
        }
    }
    for r in 0..m {
        offsets[r + 1] += offsets[r];
    }
    let mut next = offsets[..m].to_vec();
    let mut columns = vec![0; offsets[m]];
    for (c, (_, rows, _)) in cols.iter().enumerate() {
        for &r in rows {
            columns[next[r]] = c;
            next[r] += 1;
        }
    }
    (offsets, columns)
}

pub(super) fn phases<R: Real>(
    m: usize,
    cols: Cols<'_>,
) -> (
    Vec<(usize, usize)>,
    Vec<u8>,
    crate::lu_pairs::Overrides<R>,
    usize,
    Option<PairGateRecord>,
) {
    let n = cols.len();
    let mut rc = vec![0usize; m];
    for (is_slack, idx, _) in cols.iter() {
        if *is_slack {
            rc[idx[0]] += 1 + m;
        } else {
            for &i in idx {
                rc[i] += 1;
            }
        }
    }
    let (row_offsets, rowcols) = row_columns(m, cols);
    let mut row_done = vec![false; m];
    let mut col_done = vec![false; n];
    let mut seq = Vec::with_capacity(m);
    let mut phase: Vec<u8> = Vec::with_capacity(m);
    let mut ph = 1u8;
    macro_rules! pivot {
        ($i:expr, $p:expr) => {{
            row_done[$i] = true;
            col_done[$p] = true;
            seq.push(($i, $p));
            phase.push(ph);
        }};
    }
    // Phase 1: rows with one structural column and no slack, LIFO stack.
    let mut st: Vec<usize> = (0..m).filter(|&i| rc[i] == 1).collect();
    while let Some(i) = st.pop() {
        if rc[i] == 0 || row_done[i] {
            continue;
        }
        rc[i] = 0;
        let Some(pos) = rowcols[row_offsets[i]..row_offsets[i + 1]]
            .iter()
            .copied()
            .find(|&p| !col_done[p] && !cols[p].0)
        else {
            continue;
        };
        pivot!(i, pos);
        for &r2 in &cols[pos].1 {
            if !row_done[r2] {
                rc[r2] -= 1;
                if rc[r2] == 1 {
                    st.push(r2);
                }
            }
        }
    }
    // Phase 2: slack columns whose row has only the slack left.
    ph = 2;
    for pos in 0..n {
        if col_done[pos] || !cols[pos].0 {
            continue;
        }
        let i = cols[pos].1[0];
        if rc[i] == m + 1 {
            rc[i] = 1;
            pivot!(i, pos);
        }
    }
    // Phase 3: one ascending pass.
    ph = 3;
    for pos in 0..n {
        if col_done[pos] {
            continue;
        }
        let (is_slack, idx, _) = &cols[pos];
        if !*is_slack {
            if idx.len() == 1 && rc[idx[0]] <= m && !row_done[idx[0]] {
                let i = idx[0];
                rc[i] = 1;
                pivot!(i, pos);
            }
        } else {
            let i = idx[0];
            if !row_done[i] {
                rc[i] = rc[i].wrapping_sub(m);
                pivot!(i, pos);
            }
        }
    }
    let (row_done3, col_done3) = (row_done.clone(), col_done.clone());
    // Phase 4: column singletons over the unassigned rows.
    ph = 4;
    let mut cnt = vec![0i64; n];
    for (pos, (_, idx, _)) in cols.iter().enumerate() {
        if !col_done[pos] {
            cnt[pos] = idx.iter().filter(|&&i| !row_done[i]).count() as i64;
        }
    }
    let mut st: Vec<usize> = (0..n).filter(|&p| !col_done[p] && cnt[p] == 1).collect();
    while let Some(pos) = st.pop() {
        if cnt[pos] <= 0 || col_done[pos] {
            continue;
        }
        let i = *cols[pos].1.iter().find(|&&i| !row_done[i]).unwrap();
        pivot!(i, pos);
        cnt[pos] = 0;
        for &c2 in &rowcols[row_offsets[i]..row_offsets[i + 1]] {
            if !col_done[c2] && c2 != pos {
                cnt[c2] -= 1;
                if cnt[c2] == 1 {
                    st.push(c2);
                }
            }
        }
    }
    let mut ov = crate::lu_pairs::Overrides::<R>::default();
    let loff = seq.len();
    let pair_gate_record = pair_gate(m - seq.len());
    if pair_gate_record.is_some() && std::env::var_os("ROSE_NO_PAIRS").is_none() {
        let before = seq.len();
        let run = crate::lu_pairs::reference::run::<R>;
        ov = run(m, cols, &mut seq, &mut phase, &row_done3, &col_done3, &row_done, &col_done);
        pair_gate_end(seq.len() - before);
    }
    (seq, phase, ov, loff, pair_gate_record)
}


/// The entries of the singleton steps' columns that go to U plus the U entries of the singleton steps' rows (the scan of the former
/// `Structure::work`).
fn u_sing_scan<R: Real>(st: &Structure<R>, cols: Cols<'_>) -> usize {
    let nsing = st.phase.iter().filter(|&&p| p != 0).count();
    let mut step_r = vec![usize::MAX; st.m];
    let mut step_c = vec![usize::MAX; cols.len()];
    for (k, &(r, c)) in st.pivots.iter().enumerate() {
        step_r[r] = k;
        step_c[c] = k;
    }
    let mut per = vec![0usize; nsing];
    for (pos, c) in cols.iter().enumerate() {
        for &i in &c.1 {
            if step_r[i] < nsing && step_c[pos] > step_r[i] {
                per[step_r[i]] += 1;
            }
        }
    }
    // A column pair's second pivot row is updated by the first pivot's row: the reference counts the entries that
    // result (fill included, cancelled ones left out).
    for (&s, &n) in &st.pair_unz {
        per[s] = n;
    }
    per.iter().sum::<usize>()
}

/// `Structure::work` with the scans over the basis columns.
pub(super) fn work<R: Real>(st: &Structure<R>, cols: Cols<'_>, levels: bool) -> Vec<(u32, f64)> {
    let m = st.m as f64;
    let nnz_b: usize = cols.iter().map(|c| if c.0 { 1 } else { c.1.len() }).sum();
    let nsing = st.phase.iter().filter(|&&p| p != 0).count();
    let nnz_l: usize = st.l_idx.len();
    let nnz_l_sing: usize = st.l_start[nsing];
    let u_sing = u_sing_scan(st, cols);
    let mut out = vec![
        (0xd4ded, 5.0 * m + 100.0 + 3.0 * (nnz_l_sing + u_sing) as f64),
        (0xd5300, 10.0 * m + 3.0 * nnz_b as f64 + 7.0 * st.nuc_nnz.max(1) as f64),
    ];
    if st.nuc_work > 0.0 {
        out.push((0xd85e6, st.nuc_work));
    }
    if st.dense_n > 0 {
        let n = st.dense_n as f64;
        out.push((0xd8973, 0.15 * n * (n + 1.0)));
    }
    let u_rows_nnz: usize = st.u_ent.iter().filter(|e| e.1.abs() > R::c(1e-13)).count();
    let nnz_lu = (nnz_l + u_sing + u_rows_nnz) as f64;
    let a = 24.0 * nnz_lu;
    let g = if a <= 2e6 {
        1.0
    } else if a < 4e6 {
        1.0 + (a - 2e6) * 9.0 / 2e6
    } else {
        10.0
    };
    out.push((0xdc7b2, (2.0 * m + nnz_lu) * 2.0 * g));
    if levels {
        out.push((0xdcb7b, 2.0 * (nnz_l as f64 + m)));
    }
    out
}

/// Whether two structures are the same bit for bit (the shadow oracle and the differential tests).
pub(super) fn assert_same<R: Real>(new: &Structure<R>, old: &Structure<R>, cols: Cols<'_>) {
    fn same_bits<R: Real>(a: R, b: R) -> bool {
        use std::any::Any;
        let (a, b) = (&a as &dyn Any, &b as &dyn Any);
        if let Some(a) = a.downcast_ref::<f64>() {
            a.to_bits() == b.downcast_ref::<f64>().unwrap().to_bits()
        } else {
            let a = a.downcast_ref::<crate::f80::F80>().expect("unsupported Real type");
            a.repr() == b.downcast_ref::<crate::f80::F80>().unwrap().repr()
        }
    }
    fn same<R: Real>(a: &[R], b: &[R], what: &str) {
        assert!(a.len() == b.len() && a.iter().zip(b).all(|(&x, &y)| same_bits(x, y)), "lu_order shadow: {what}");
    }
    assert_eq!(new.pivots, old.pivots, "lu_order shadow: pivots");
    assert_eq!(new.phase, old.phase, "lu_order shadow: phase");
    assert_eq!(new.l_start, old.l_start, "lu_order shadow: l_start");
    assert_eq!(new.l_idx, old.l_idx, "lu_order shadow: l_idx");
    same(&new.l_val, &old.l_val, "l_val");
    same(&new.l_raw, &old.l_raw, "l_raw");
    same(&new.diag, &old.diag, "diag");
    assert_eq!(new.u_start, old.u_start, "lu_order shadow: u_start");
    assert!(
        new.u_ent.len() == old.u_ent.len()
            && new.u_ent.iter().zip(&old.u_ent).all(|(x, y)| x.0 == y.0 && same_bits(x.1, y.1)),
        "lu_order shadow: u_ent"
    );
    assert_eq!(new.nuc_nnz, old.nuc_nnz, "lu_order shadow: nuc_nnz");
    assert_eq!(new.nuc_work.to_bits(), old.nuc_work.to_bits(), "lu_order shadow: nuc_work");
    assert_eq!(new.loff, old.loff, "lu_order shadow: loff");
    assert_eq!(new.pair_unz, old.pair_unz, "lu_order shadow: pair_unz");
    assert_eq!(new.pair_urow, old.pair_urow, "lu_order shadow: pair_urow");
    assert_eq!(new.dense_n, old.dense_n, "lu_order shadow: dense_n");
    assert_eq!(new.m, old.m);
    assert_eq!((new.nsing, new.nnz_b, new.u_sing), (old.nsing, old.nnz_b, old.u_sing), "lu_order shadow: work inputs");
    assert_eq!(
        new.pair_gate_record.map(|r| (r.remaining, r.att, r.gate, r.fac_no)),
        old.pair_gate_record.map(|r| (r.remaining, r.att, r.gate, r.fac_no)),
        "lu_order shadow: pair gate record"
    );
    for levels in [false, true] {
        let a = new.work(cols, levels);
        let b = work(old, cols, levels);
        assert!(
            a.len() == b.len() && a.iter().zip(&b).all(|(x, y)| x.0 == y.0 && x.1.to_bits() == y.1.to_bits()),
            "lu_order shadow: work (levels {levels})"
        );
    }
}

/// The shadow check of a factorization: the old construction on the same input and gate state, compared with `new`.
#[cfg(feature = "factor-shadow")]
pub(super) fn check<R: Real>(new: &Structure<R>, m: usize, cols: Cols<'_>, gate_before: (i64, i64)) {
    let gate_after = gate_get();
    gate_set(gate_before);
    let old = crate::dump::muted(|| factor_structure_r::<R>(m, cols));
    assert_eq!(gate_get(), gate_after, "lu_order shadow: pair gate state");
    assert_same(new, &old, cols);
}
