//! Differential tests of the numeric U row pool (`URows`) against the three structures it replaced: the transpose `ut` (sorted by column),
//! the member lists `rowl` (physical order) and the stamped gather `TauRow`. The old kernels are retained here as nested-list copies:
//! the needed nucleus dot (gather, then `row_dot`/`row_dot_q`), the sparse and dense U^T solves of the BTRAN, and the row elimination of
//! `update` (`update_old::eliminate_old`). Random factors and update sequences, `f64` and extended precision, signed zeros in the rows
//! and in the vectors.
use super::ab::same_bits;
use super::*;
use crate::f80::F80;

/// The old row dot: the row's `(value, x)` terms in `rowl` order, the value of the first duplicate of a column (none occur), then the
/// dot product kernels on the gathered terms.
fn dot_old<R: Real>(f: &DualFactor<R>, ut: &[Vec<(usize, R)>], rowl: &[Vec<usize>], p: usize, x: &[R]) -> R {
    let terms: Vec<(R, R)> = rowl[p].iter().filter_map(|&c| ut[p].iter().find(|e| e.0 == c).map(|e| (e.1, x[c]))).collect();
    if R::QUAD {
        row_dot_q(x[p], &terms)
    } else {
        row_dot(x[p], &terms, f.row_par[p] as usize)
    }
}

/// The sparse U^T solve as it was on the nested transpose (`sparse_btran_u` before the row pool).
fn sparse_btran_u_old<R: Real, W: BtranTail<R>>(
    f: &DualFactor<R>,
    ut: &[Vec<(usize, R)>],
    inp: &[(usize, R)],
    x: &mut [R],
    w: &mut W,
) -> (Vec<usize>, f64) {
    let m = f.m;
    let cnt = inp.len();
    let mut out: Vec<usize> = Vec::new();
    let mut eu = 0.0f64;
    let push = |w: &mut W, row: &[(usize, R)], d: R| {
        for &(c, v) in row {
            let value = w.entry(c, &f.slot);
            *value = -(v * d) + *value;
        }
    };
    let mut lowest = m;
    if f.h33 > 5.0 {
        eu = cnt as f64 * 3.0;
    }
    for &(p, v) in inp {
        let d = v / f.diag[p];
        let tail = if f.h33 > 5.0 { f.kind[p] == -2 } else { f.slot[p] >= f.i21 };
        if tail {
            if f.h33 <= 5.0 || d.abs() > R::c(TOL) {
                lowest = lowest.min(f.slot[p]);
                *w.entry(p, &f.slot) += v;
            }
        } else if d.abs() > R::c(TOL) {
            out.push(p);
            x[p] = d;
            if f.h33 > 5.0 && f.kind[p] >= 0 {
                eu = ut[p].len() as f64 * 2.0 + eu;
            }
            push(w, &ut[p], d);
            for &(c, _) in &ut[p] {
                lowest = lowest.min(f.slot[c]);
            }
        }
    }
    if f.h33 > 5.0 {
        eu = (m - lowest) as f64 * 2.0 + eu;
    }
    let mut next = lowest;
    while let Some(t) = w.next_slot(next, m) {
        next = t + 1;
        let p = f.at[t];
        let d = w.take(p, t) / f.diag[p];
        if d.abs() > R::c(TOL) {
            eu = ut[p].len() as f64 * 2.0 + eu;
            out.push(p);
            x[p] = d;
            push(w, &ut[p], d);
        }
    }
    (out, eu)
}

/// The dense U^T loop of the BTRAN on the nested transpose.
fn dense_btran_u_old<R: Real>(f: &DualFactor<R>, ut: &[Vec<(usize, R)>], x: &mut [R], mut e: f64) -> f64 {
    for t in 0..f.m {
        let p = f.at[t];
        let xv = x[p];
        let d = xv / f.diag[p];
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
    e
}

fn same_vec<R: Real>(a: &[R], b: &[R], what: &str) {
    assert_eq!(a.len(), b.len(), "{what} length");
    for (i, (&p, &q)) in a.iter().zip(b).enumerate() {
        assert!(same_bits(p, q), "{what}[{i}]");
    }
}

struct Counts {
    dots: usize,
    long_dots: usize,
    sparse: usize,
    dense: usize,
    elim: usize,
    elim_eta: usize,
}

/// Every retained old kernel against the pool kernel on the current state of `f` (a clone with some values replaced by signed zeros and
/// large or tiny numbers, which no kernel treats specially but a different order or reduction tree would show).
fn compare_kernels<R: Real>(f: &DualFactor<R>, next: &mut impl FnMut() -> usize, n: &mut Counts) {
    let m = f.m;
    let mut g = f.clone();
    let pool = [R::ZERO, -R::ZERO, R::c(0.5), R::c(-0.75), R::c(1e-14), R::c(3e7), R::c(-1.0)];
    for y in 0..m {
        for v in g.urows.row_mut(y).1 {
            if next() % 5 == 0 {
                *v = pool[next() % pool.len()];
            }
        }
    }
    let (ut, rowl) = (g.urows.to_ut(), g.urows.to_rowl());
    // Needed nucleus dot.
    for _ in 0..3 {
        let x: Vec<R> = (0..m).map(|_| if next() % 4 == 0 { pool[next() % pool.len()] } else { R::c((next() % 17) as f64 * 0.25 - 2.0) }).collect();
        for p in 0..m {
            let (new, old) = (g.nucleus_dot(p, &x), dot_old(&g, &ut, &rowl, p, &x));
            assert!(same_bits(new, old), "nucleus dot of row {p} ({} terms)", rowl[p].len());
            n.dots += 1;
            n.long_dots += usize::from(rowl[p].len() >= 8);
        }
    }
    // Sparse and dense U^T solves.
    for h in [2.0, 7.0] {
        g.h33 = h;
        for _ in 0..8 {
            let inp: Vec<(usize, R)> = (0..1 + next() % 6).map(|_| (next() % m, pool[next() % pool.len()])).collect();
            let (mut xo, mut xs, mut xd) = (vec![R::ZERO; m], vec![R::ZERO; m], vec![R::ZERO; m]);
            let old = sparse_btran_u_old(&g, &ut, &inp, &mut xo, &mut vec![R::ZERO; m]);
            let dense = g.sparse_btran_u(&inp, &mut xd, &mut vec![R::ZERO; m], Vec::new());
            let mut tail = SparseBtranTail::default();
            tail.reset(m);
            let sparse = g.sparse_btran_u(&inp, &mut xs, &mut tail, Vec::new());
            for new in [&dense, &sparse] {
                assert_eq!(new.0, old.0, "sparse U^T order");
                assert_eq!(new.1.to_bits(), old.1.to_bits(), "sparse U^T work");
            }
            same_vec(&xd, &xo, "sparse U^T (dense tail) output");
            same_vec(&xs, &xo, "sparse U^T (retained tail) output");
            n.sparse += 1;
        }
        let mut xo: Vec<R> = (0..m).map(|_| if next() % 3 == 0 { pool[next() % pool.len()] } else { R::ZERO }).collect();
        let mut xn = xo.clone();
        let (eo, en) = (dense_btran_u_old(&g, &ut, &mut xo, 5.0), g.dense_btran_u(&mut xn, 5.0));
        assert_eq!(eo.to_bits(), en.to_bits(), "dense U^T work");
        same_vec(&xn, &xo, "dense U^T output");
        n.dense += 1;
    }
    // Row elimination.
    for _ in 0..6 {
        let r = next() % m;
        let (mut wn, mut wo) = (vec![R::ZERO; m], vec![R::ZERO; m]);
        let (mut sn, mut so) = (vec![false; m], vec![false; m]);
        let (mut un, mut uo) = (btran_ws::Update::<R>::default(), btran_ws::Update::<R>::default());
        let vn = g.eliminate(r, &mut wn, &mut sn, &mut un);
        let en = g.dense_est.get();
        let vo = g.eliminate_old(&ut, r, &mut wo, &mut so, &mut uo);
        let eo = g.dense_est.get();
        assert_eq!(vn, vo, "eliminate: solve length");
        assert_eq!(en.to_bits(), eo.to_bits(), "eliminate: dense estimate");
        assert_eq!(un.eta.len(), uo.eta.len(), "eliminate: eta length");
        for (a, b) in un.eta.iter().zip(&uo.eta) {
            assert!(a.0 == b.0 && same_bits(a.1, b.1), "eliminate: eta entry");
        }
        same_vec(&wn, &wo, "eliminate: work vector");
        assert_eq!(sn, so, "eliminate: markers");
        n.elim += 1;
        n.elim_eta += un.eta.len();
    }
}

fn run<R: Real>(quad: bool) {
    let mut seed = 0x2545f4914f6cdd1du64 ^ u64::from(quad);
    let mut next = || {
        seed ^= seed << 13;
        seed ^= seed >> 7;
        seed ^= seed << 17;
        seed as usize
    };
    let mut n = Counts { dots: 0, long_dots: 0, sparse: 0, dense: 0, elim: 0, elim_eta: 0 };
    let (mut tail_seen, mut grown) = (0, 0);
    for m in [24, 70, 200] {
        for round in 0..8 {
            let mut cols: Vec<_> = (0..m)
                .map(|i| {
                    let mut rows = vec![i];
                    let mut vals = vec![4.0];
                    for _ in 0..(1 + round % 6) * (1 + usize::from(round % 3 == 2) * 3) {
                        let r = if round % 2 == 1 {
                            if i == 0 {
                                continue;
                            }
                            next() % i
                        } else {
                            next() % m
                        };
                        if !rows.contains(&r) {
                            rows.push(r);
                            vals.push(if next() % 2 == 0 { 0.5 } else { -0.75 });
                        }
                    }
                    (true, rows, vals)
                })
                .collect();
            let st = crate::lu_order::factor_structure_r::<R>(m, &cols);
            let mut f = DualFactor::<R>::fresh(m, m, &st, &cols);
            compare_kernels(&f, &mut next, &mut n);
            for step in 0..m.min(60) {
                let r = next() % m;
                let (rows, vals): (Vec<usize>, Vec<f64>) = if next() % 6 == 0 {
                    (cols[r].1.clone(), cols[r].2.iter().map(|v| 2.0 * v).collect())
                } else {
                    let mut rows = vec![r];
                    let mut vals = vec![5.0 + step as f64];
                    for _ in 0..(1 + next() % 8) {
                        let q = next() % m;
                        if !rows.contains(&q) {
                            rows.push(q);
                            vals.push([0.25, -0.5, 1.5][next() % 3]);
                        }
                    }
                    (rows, vals)
                };
                let rv: Vec<R> = vals.iter().map(|&v| R::c(v)).collect();
                f.ftran(&rows, &rv);
                tail_seen += usize::from(f.kind[r] == -2);
                let before = (0..m).map(|y| f.urows.len(y)).sum::<usize>();
                f.update(r);
                grown += usize::from((0..m).map(|y| f.urows.len(y)).sum::<usize>() > before);
                f.urows.check(&f.uc, &f.kind);
                cols[r] = (true, rows, vals);
                if step % 7 == 0 {
                    compare_kernels(&f, &mut next, &mut n);
                }
            }
            compare_kernels(&f, &mut next, &mut n);
            f.recycle();
        }
    }
    assert!(n.dots > 1000 && n.long_dots > 100, "dots {} (>= 8 terms {})", n.dots, n.long_dots);
    assert!(n.sparse > 100 && n.dense > 20 && n.elim > 100 && n.elim_eta > 50, "sparse {} dense {} elim {} eta {}", n.sparse, n.dense, n.elim, n.elim_eta);
    assert!(tail_seen > 0 && grown > 0, "tail replacements {tail_seen}, growing updates {grown}");
}

#[test]
fn row_pool_kernels_match_the_nested_lists_f64() {
    run::<f64>(false);
}

#[test]
fn row_pool_kernels_match_the_nested_lists_extended() {
    run::<F80>(true);
}

/// The pool itself: pushes past the row room move the row, `swap_remove` and `clear_row` act as on `Vec`.
#[test]
fn row_pool_edits_follow_vec_semantics() {
    let mut seed = 0x9e3779b97f4a7c15u64;
    let mut next = || {
        seed ^= seed << 13;
        seed ^= seed >> 7;
        seed ^= seed << 17;
        seed as usize
    };
    let m = 12;
    let mut pool = URows::<f64>::new();
    pool.reset(m);
    let mut want: Vec<Vec<(usize, f64)>> = vec![Vec::new(); m];
    for step in 0..4000 {
        let y = next() % m;
        match next() % 7 {
            0 | 1 | 2 | 3 => {
                let c = next() % 64;
                if !want[y].iter().any(|e| e.0 == c) {
                    let v = step as f64 * if next() % 2 == 0 { 1.0 } else { -0.0 };
                    if want[y].is_empty() && next() % 2 == 0 {
                        pool.reserve_row(y, 1 + next() % 5);
                    }
                    pool.push(y, c, v);
                    want[y].push((c, v));
                }
            }
            4 | 5 => {
                if !want[y].is_empty() {
                    let p = next() % want[y].len();
                    assert_eq!(pool.position(y, want[y][p].0), Some(p));
                    pool.swap_remove(y, p);
                    want[y].swap_remove(p);
                }
            }
            _ => {
                if next() % 8 == 0 {
                    pool.clear_row(y);
                    want[y].clear();
                }
            }
        }
        for y in 0..m {
            let (c, v) = pool.row(y);
            assert_eq!(c.len(), want[y].len());
            for ((&c, &v), w) in c.iter().zip(v).zip(&want[y]) {
                assert!(c as usize == w.0 && v.to_bits() == w.1.to_bits());
            }
        }
    }
}
