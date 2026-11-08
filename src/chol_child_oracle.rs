// Complete original child kernels, frozen before execution planning.
use super::*;
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
pub(super) fn child_update(
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

