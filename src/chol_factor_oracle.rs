// Original complete assembly/child/block/factor pipeline, frozen before execution ownership changes.
use super::*;
pub fn assemble(
    f: &StdForm,
    rows: &[Vec<(usize, usize)>],
    d: &[f64],
    mp: &[usize],
    col: usize,
    w: usize,
    c: usize,
) -> (Vec<f64>, Vec<f64>) {
    let (idx, val, st) = (&f.row_idx, &f.val, &f.col_start);
    let (n1, n2) = (f.m, f.n2);
    let mut mm = vec![0.0; w * w];
    let mut r = vec![0.0; (c - w) * w];
    for i in col..col + w {
        for &(cj, pos) in &rows[i] {
            let s = st[cj];
            let e = st[cj + 1];
            if !(pos == s || idx[pos - 1] < col) {
                continue;
            }
            let mut l24 = pos;
            loop {
                l24 += 1;
                if e <= l24 || mp[idx[l24]] >= w {
                    break;
                }
            }
            let dd = d[cj];
            let mut j = pos;
            let sj = f.scale.as_ref().map_or(0.0, |v| v[cj]);
            if sj != 0.0 {
                // Scaled column (`Q+0x648[cj] != 0`): every contribution is the number (s s) d (triples) or s (d s) (the rest), added without a product with A.
                let (sv3, sv1) = ((sj * sj) * dd, sj * (dd * sj));
                while j + 2 < l24 {
                    let (a, b, cc) = (mp[idx[j]], mp[idx[j + 1]], mp[idx[j + 2]]);
                    mm[a * w + a] = rf(sv3 + mm[a * w + a]);
                    mm[b * w + a] = rf(sv3 + mm[b * w + a]);
                    mm[b * w + b] = rf(sv3 + mm[b * w + b]);
                    for q in j + 2..l24 {
                        let mq = mp[idx[q]];
                        add_m(&mut mm, &mut r, w, a, mq, sv3);
                        add_m(&mut mm, &mut r, w, b, mq, sv3);
                        add_m(&mut mm, &mut r, w, cc, mq, sv3);
                    }
                    for q in l24..e {
                        let mq = mp[idx[q]];
                        if mq < n1 {
                            add_m(&mut mm, &mut r, w, a, mq, sv3);
                            add_m(&mut mm, &mut r, w, b, mq, sv3);
                            add_m(&mut mm, &mut r, w, cc, mq, sv3);
                        }
                    }
                    j += 3;
                }
                while j < l24 {
                    let a = mp[idx[j]];
                    mm[a * w + a] = rf(sv1 + mm[a * w + a]);
                    for q in j + 1..l24 {
                        add_m(&mut mm, &mut r, w, a, mp[idx[q]], sv1);
                    }
                    for q in l24..e {
                        let mq = mp[idx[q]];
                        if mq < n1 {
                            add_m(&mut mm, &mut r, w, a, mq, sv1);
                        }
                    }
                    j += 1;
                }
                continue;
            }
            while j + 2 < l24 {
                let (a, b, cc) = (mp[idx[j]], mp[idx[j + 1]], mp[idx[j + 2]]);
                let (a0, a1, a2) = (val[j], val[j + 1], val[j + 2]);
                let (t0, t1, t2) = (a0 * dd, a1 * dd, a2 * dd);
                mm[a * w + a] = rf(a0 * t0 + mm[a * w + a]);
                mm[b * w + a] = rf(a1 * t0 + mm[b * w + a]);
                mm[b * w + b] = rf(a1 * t1 + mm[b * w + b]);
                for q in j + 2..l24 {
                    let (mq, v) = (mp[idx[q]], val[q]);
                    add_m(&mut mm, &mut r, w, a, mq, t0 * v);
                    add_m(&mut mm, &mut r, w, b, mq, t1 * v);
                    add_m(&mut mm, &mut r, w, cc, mq, t2 * v);
                }
                for q in l24..e {
                    let (mq, v) = (mp[idx[q]], val[q]);
                    if mq < n1 {
                        add_m(&mut mm, &mut r, w, a, mq, t0 * v);
                        add_m(&mut mm, &mut r, w, b, mq, t1 * v);
                        add_m(&mut mm, &mut r, w, cc, mq, t2 * v);
                    }
                }
                j += 3;
            }
            while j < l24 {
                let (a, v0) = (mp[idx[j]], val[j]);
                let t = dd * v0;
                mm[a * w + a] = rf(v0 * t + mm[a * w + a]);
                for q in j + 1..l24 {
                    add_m(&mut mm, &mut r, w, a, mp[idx[q]], val[q] * t);
                }
                for q in l24..e {
                    let mq = mp[idx[q]];
                    if mq < n1 {
                        add_m(&mut mm, &mut r, w, a, mq, val[q] * t);
                    }
                }
                j += 1;
            }
        }
        if f.sense[i] != b'=' {
            let k = (i - col) * w + (i - col);
            mm[k] = rf(mm[k] + d[n2 + i]);
        }
    }
    (mm, r)
}

/// `0x2f77f0` (used for `w < 3`, `c * w > 30000`): for every supernode column `i` and every column `cj` of row `i` (ascending), `t = D[cj] A[i][cj]`,
/// then `acc[row_q] = A[q] t + acc[row_q]` over the entries of `cj` from `i`'s position on; column `i` of the front is `acc` (the diagonal gets the slack `D`).
pub fn assemble_generic(
    f: &StdForm,
    rows: &[Vec<(usize, usize)>],
    d: &[f64],
    mp: &[usize],
    col: usize,
    w: usize,
    c: usize,
) -> (Vec<f64>, Vec<f64>) {
    let (idx, val, st) = (&f.row_idx, &f.val, &f.col_start);
    let n2 = f.n2;
    let mut mm = vec![0.0; w * w];
    let mut r = vec![0.0; (c - w) * w];
    let mut acc = vec![0.0; c];
    let mut seen = vec![false; c];
    let mut touched: Vec<usize> = Vec::new();
    for i in col..col + w {
        let a = i - col;
        for &(cj, pos) in &rows[i] {
            let e = st[cj + 1];
            let sj = f.scale.as_ref().map_or(0.0, |v| v[cj]);
            // Scaled column (`0x2f81fd`): every contribution is the number (s s) D[cj].
            let t = if sj != 0.0 {
                (sj * sj) * d[cj]
            } else {
                d[cj] * val[pos]
            };
            for q in pos..e {
                let rr = idx[q];
                let m = mp[rr];
                if m < c && rr >= i {
                    if !seen[m] {
                        seen[m] = true;
                        touched.push(m);
                    }
                    acc[m] = if sj != 0.0 {
                        t + acc[m]
                    } else {
                        val[q] * t + acc[m]
                    };
                }
            }
        }
        // The slack of an inequality row is added to the double sum before the one rounding to float.
        let slack = f.sense[i] != b'=';
        if slack && !seen[a] {
            seen[a] = true;
            touched.push(a);
        }
        for &m in &touched {
            if m < w {
                mm[m * w + a] = rf(if slack && m == a {
                    acc[m] + d[n2 + i]
                } else {
                    acc[m]
                });
            } else {
                r[(m - w) * w + a] = rf(acc[m]);
            }
            acc[m] = 0.0;
            seen[m] = false;
        }
        touched.clear();
    }
    (mm, r)
}

fn factor_block(
    mm: &mut [f64],
    r: &mut [f64],
    w: usize,
    c: usize,
    orig: &[f64],
    reg: Option<(usize, usize)>,
) -> (Vec<f64>, Vec<f64>) {
    let dense = (1..w).any(|i| (0..i).any(|j| mm[i * w + j] != 0.0));
    let (d, rcp) = if dense {
        dense_ldl(mm, w, orig, reg)
    } else {
        diag_only(mm, w, orig, reg)
    };
    let rows = c - w;
    if rows != 0 {
        if dense {
            let big32 = SINGLE.with(|c| c.get()) && w >= 2;
            let l32: Vec<f32> = if big32 {
                mm.iter().map(|&v| v as f32).collect()
            } else {
                Vec::new()
            };
            for rr in 0..rows {
                let x = if big32 {
                    let mut x: Vec<f32> =
                        r[rr * w..(rr + 1) * w].iter().map(|&v| v as f32).collect();
                    crate::sgemm32::strsm_llu_col(&l32, w, &mut x, rr, rows);
                    x.iter().map(|&v| v as f64).collect()
                } else if SINGLE.with(|c| c.get()) {
                    strsm_seq(mm, w, &r[rr * w..(rr + 1) * w])
                } else {
                    dtrsm_llu_col(mm, w, &r[rr * w..(rr + 1) * w], rr, rows)
                };
                r[rr * w..(rr + 1) * w].copy_from_slice(&x);
            }
        }
        for rr in 0..rows {
            for k in 0..w {
                r[rr * w + k] = rf(r[rr * w + k] * if rcp[k] != 0.0 { rf(rcp[k]) } else { 1.0 });
            }
        }
    }
    (d, rcp)
}

fn child_update_w1(
    ch: &Child,
    wp: usize,
    mp: &[usize],
    mpm: &mut [f64],
    rp: &mut [f64],
    cnt: usize,
    e: usize,
) {
    let (cc, l, dd) = (ch.c, ch.packed, ch.d[0]);
    let mpr: Vec<usize> = ch.rows.iter().map(|&r| mp[r]).collect();
    if wp < 5 {
        for t in cnt..e {
            let f = (-l[t]) * dd;
            for u in t..cc {
                if u < e {
                    let k = mpr[u] * wp + mpr[t];
                    mpm[k] = rf(l[u] * f + mpm[k]);
                } else {
                    let k = (mpr[u] - wp) * wp + mpr[t];
                    rp[k] = rf(f * l[u] + rp[k]);
                }
            }
        }
        return;
    }
    for t in cnt..cc {
        let f = (-l[t]) * dd;
        if t < e {
            for u in cnt..=t {
                let k = mpr[t] * wp + mpr[u];
                mpm[k] = rf(l[u] * f + mpm[k]);
            }
        } else {
            for u in cnt..e {
                let k = (mpr[t] - wp) * wp + mpr[u];
                rp[k] = rf(f * l[u] + rp[k]);
            }
        }
    }
}

/// `0x2fcbd0`, child `ch` into the front `(Mp, Rp)` of the parent (width `wp`, first column `col`). Rows of the child's list from `cnt` (first row inside
/// the parent's columns) to `e` (end of the parent's columns) are the `m` update rows, the rest the rect rows. `temp[t][k] = (-L[t][k]) d[k]`;
/// DGEMM('T','N', k = child width) blocks of 24 columns for the diagonal part (contiguous map: directly into `Mp`, else into a zero scratch and
/// scatter-add), then blocks of rect rows.
fn child_update(
    ch: &Child,
    wp: usize,
    col: usize,
    sym80: usize,
    mp: &[usize],
    mpm: &mut [f64],
    rp: &mut [f64],
) {
    let (wc, cc, rows_c, lp, d) = (ch.w, ch.c, ch.rows, ch.packed, ch.d);
    let cnt = (wc..cc)
        .find(|&t| rows_c[t] >= col)
        .expect("child without a row in the parent");
    let mut e = cnt;
    while e < cc && rows_c[e] < col + wp {
        e += 1;
    }
    if wc == 1 {
        child_update_w1(ch, wp, mp, mpm, rp, cnt, e);
        return;
    }
    let m = e - cnt;
    let tri = wc * (wc + 1) / 2;
    let lrow = |t: usize| &lp[tri + (t - wc) * wc..tri + (t - wc + 1) * wc];
    let mut temp = Vec::with_capacity(m * wc);
    for t in cnt..e {
        let lr = lrow(t);
        for k in 0..wc {
            temp.push(rf((-lr[k]) * d[k]));
        }
    }
    let first = mp[rows_c[cnt]];
    let contiguous = mp[rows_c[e - 1]] == first + m - 1;
    let mut j0 = 0;
    while j0 < m {
        let nb = (m - j0).min(24);
        let mut b = Vec::with_capacity(nb * wc);
        for t in cnt + j0..cnt + j0 + nb {
            b.extend_from_slice(lrow(t));
        }
        if contiguous {
            let mmr = j0 + nb;
            let mut cm = Vec::with_capacity(mmr * nb);
            for jj in 0..nb {
                for ii in 0..mmr {
                    cm.push(mpm[(first + j0 + jj) * wp + first + ii]);
                }
            }
            dgemm_t_n(mmr, nb, wc, &temp, wc, &b, wc, &mut cm, mmr);
            for jj in 0..nb {
                for ii in 0..mmr {
                    mpm[(first + j0 + jj) * wp + first + ii] = rf(cm[jj * mmr + ii]);
                }
            }
        } else {
            let mut cn = vec![0.0; m * nb];
            // Single precision: the call covers the rows `0..j0 + nb` only (the float code does not compute the part above the diagonal block).
            let mr = if SINGLE.with(|c| c.get()) { j0 + nb } else { m };
            dgemm_t_n(mr, nb, wc, &temp, wc, &b, wc, &mut cn, m);
            for jj in 0..nb {
                let t = j0 + jj;
                for u in 0..=t {
                    let (r, cl) = (mp[rows_c[cnt + t]], mp[rows_c[cnt + u]]);
                    mpm[r * wp + cl] = rf(mpm[r * wp + cl] + cn[jj * m + u]);
                }
            }
        }
        j0 += nb;
    }
    let mut t0 = e;
    let blk = if m != 0 {
        (0xc0usize).min(sym80 / m).max(1)
    } else {
        1
    };
    while t0 < cc {
        let mut t1 = (t0 + blk).min(cc);
        let r_first = mp[rows_c[t0]];
        if contiguous && mp[rows_c[t1 - 1]] == r_first + (t1 - t0) - 1 {
            while t1 < cc && mp[rows_c[t1]] == r_first + (t1 - t0) {
                t1 += 1;
            }
            let nr = t1 - t0;
            let mut b = Vec::with_capacity(nr * wc);
            for t in t0..t1 {
                b.extend_from_slice(lrow(t));
            }
            let mut cm = Vec::with_capacity(m * nr);
            for jj in 0..nr {
                for ii in 0..m {
                    cm.push(rp[(r_first - wp + jj) * wp + first + ii]);
                }
            }
            dgemm_t_n(m, nr, wc, &temp, wc, &b, wc, &mut cm, m);
            for jj in 0..nr {
                for ii in 0..m {
                    rp[(r_first - wp + jj) * wp + first + ii] = rf(cm[jj * m + ii]);
                }
            }
        } else {
            let nr = t1 - t0;
            let mut b = Vec::with_capacity(nr * wc);
            for t in t0..t1 {
                b.extend_from_slice(lrow(t));
            }
            let mut cn = vec![0.0; m * nr];
            dgemm_t_n(m, nr, wc, &temp, wc, &b, wc, &mut cn, m);
            for jj in 0..nr {
                let r = mp[rows_c[t0 + jj]];
                for u in 0..m {
                    let cl = mp[rows_c[cnt + u]];
                    rp[(r - wp) * wp + cl] = rf(rp[(r - wp) * wp + cl] + cn[jj * m + u]);
                }
            }
        }
        t0 = t1;
    }
}

pub fn numeric_factor(
    f: &StdForm,
    rows: &[Vec<(usize, usize)>],
    plan: &Plan,
    d: &[f64],
    zx: Option<&[f64]>,
    mut on_sn: Option<&mut dyn FnMut(usize, &Supernode, &[f64], &[f64])>,
) -> Factor {
    let m = f.m;
    let _guard = SingleGuard::new(plan.single);
    // Dense columns: the factor has `nt = m + nd` nodes (the dense nodes last), pivots regularise (`SYM+0x110`, `Q+0x218 > 0`), and the sparse supernodes are
    // assembled as without dense columns (`0x2f83e0` for `w >= 3`, else the generic routine; `fx_bar_densecol`); the dense supernode starts from the diagonal `-1/D[j]` of its columns, and a row's dense-node entries are the raw
    // coefficients `A[i][j]` (no `D`). The diagonal of a dense node is `-zx_j` (`z/x + v/u` of the iterate; the starting point has no `zx`: `-1/D_j`). The iterate uses `-zx_j` even when `D` is clamped to `1e12`.
    let nt = plan.n;
    let reg = !f.dense.is_empty();
    let mut dnode_col = vec![NONE; nt];
    // Per dense column its coefficient in every row (0 where absent), indexed by the node `m + k`.
    let mut dval: Vec<Vec<f64>> = vec![Vec::new(); nt];
    for (r, &j) in f.dense.iter().enumerate() {
        dnode_col[m + f.dense_k[r]] = j;
        let mut v = vec![0.0; m];
        for q in f.col_start[j]..f.col_start[j + 1] {
            v[f.row_idx[q]] = f.val[q];
        }
        dval[m + f.dense_k[r]] = v;
    }
    let mut vals: Vec<Vec<f64>> = Vec::with_capacity(plan.sns.len());
    let mut dvec: Vec<Vec<f64>> = Vec::with_capacity(plan.sns.len());
    let mut diag = vec![0.0; nt];
    let mut rcp = vec![0.0; nt];
    let mut mp = vec![NONE; nt];
    let mut minratio = 1e300f64;
    for (k, s) in plan.sns.iter().enumerate() {
        let (col, w, c) = (s.col, s.w, s.c);
        for (t, &r) in s.rows.iter().enumerate() {
            mp[r] = t;
        }
        let (mut mm, mut r) = if reg && col >= m {
            let mut mm = vec![0.0; w * w];
            for t in 0..w {
                let jd = dnode_col[col + t];
                mm[t * w + t] = match zx {
                    Some(z) if z[jd] > 0.0 => -z[jd],
                    _ => -1.0 / d[jd],
                };
            }
            (mm, vec![0.0; (c - w) * w])
        } else if w >= 3 && c * w <= 30000 {
            assemble(f, rows, d, &mp, col, w, c)
        } else {
            assemble_generic(f, rows, d, &mp, col, w, c)
        };
        if reg && col < m {
            for (t, &rr) in s.rows.iter().enumerate() {
                if rr >= m && t >= w {
                    for a in 0..w {
                        r[(t - w) * w + a] = rf(dval[rr][col + a]);
                    }
                }
            }
        }
        let orig: Vec<f64> = (0..w).map(|i| mm[i * w + i]).collect();
        let front_col = if dump::on("barrier_front") {
            std::env::var("FRONT_COL")
                .ok()
                .and_then(|v| v.parse::<usize>().ok())
        } else {
            None
        };
        let front_ent: Vec<(usize, usize)> = std::env::var("FRONT_ENT")
            .unwrap_or_else(|_| "0:0,1:0,1:1".into())
            .split(',')
            .filter_map(|t| t.split_once(':'))
            .filter_map(|(a, b)| Some((a.parse().ok()?, b.parse().ok()?)))
            .collect();
        let front_rec = |tag: &str, child: i64, mm: &[f64]| {
            let v: Vec<String> = front_ent
                .iter()
                .map(|&(i, j)| format!("{:08x}", (mm[i * w + j] as f32).to_bits()))
                .collect();
            dump::emit_tag(
                "barrier_front",
                tag,
                json!({"col": col, "child": child, "m": v}),
            );
        };
        if front_col == Some(col) {
            front_rec("ASM", -1, &mm);
        }
        for &ci in &s.upd {
            let cs = &plan.sns[ci];
            let ch = Child {
                w: cs.w,
                c: cs.c,
                rows: &cs.rows,
                packed: &vals[ci],
                d: &dvec[ci],
            };
            child_update(&ch, w, col, plan.sym80, &mp, &mut mm, &mut r);
            if front_col == Some(col) {
                front_rec("UPD", ci as i64, &mm);
            }
        }
        let (dd, rc) = factor_block(
            &mut mm,
            &mut r,
            w,
            c,
            &orig,
            if reg { Some((col, m)) } else { None },
        );
        let pk = pack(&mm, &r, w);
        for t in 0..w {
            diag[col + t] = dd[t];
            rcp[col + t] = rc[t];
            if orig[t] != 0.0 {
                minratio = minratio.min(dd[t].abs() / orig[t].abs());
            }
        }
        if let Some(cb) = on_sn.as_deref_mut() {
            cb(k, s, &pk, &dd);
        }
        vals.push(pk);
        dvec.push(dd);
        for &r_ in &s.rows {
            mp[r_] = NONE;
        }
    }
    let mut par = Vec::with_capacity(vals.len());
    let mut off = 0usize;
    for v in &vals {
        par.push(off & 1);
        off += v.len();
    }
    Factor {
        n: nt,
        vals,
        diag,
        rcp,
        par,
        minratio,
    }
}

