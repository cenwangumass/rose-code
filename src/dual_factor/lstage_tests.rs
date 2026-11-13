//! L-stage kernels against the code they replaced: the eta gather, the permuted scatter, and the level solve (unchecked accesses)
//! against the checked one, bit for bit on random factors with etas.
use super::handoff_tests::{all_bits, bits};
use super::*;

struct Rng(u64);

impl Rng {
    fn next(&mut self) -> usize {
        self.0 ^= self.0 << 13;
        self.0 ^= self.0 >> 7;
        self.0 ^= self.0 << 17;
        self.0 as usize
    }
    /// Values whose products and sums are NaN of either sign (the sign of a NaN that meets another is part of the bits).
    fn special(&mut self) -> f64 {
        [0.0, -0.0, 1.0, -1.5, f64::INFINITY, f64::NEG_INFINITY, f64::from_bits(0x7ff8_0000_0000_0000), f64::from_bits(0xfff8_0000_0000_0000), 2.5, 0.0, 0.0][self.next() % 11]
    }
    fn value(&mut self) -> f64 {
        [0.25, -1.5, 1e-15, 3.0, -0.0, 0.0, TOL, -TOL, 1e-14, 7.1e3, -2.9e-3, f64::MIN_POSITIVE][self.next() % 12]
    }
}

fn gathered_dots<R: Real>(special: bool) {
    let mut g = Rng(0x2545f4914f6cdd1d);
    let mut buf: Vec<R> = Vec::new();
    for n in 0..80usize {
        for peel in 0..2usize {
            for _ in 0..40 {
                let m = n + 5;
                let mut v = |g: &mut Rng| if special { g.special() } else { g.value() };
                let x: Vec<R> = (0..m).map(|_| R::c(v(&mut g))).collect();
                let c: Vec<R> = (0..n).map(|_| R::c(v(&mut g))).collect();
                let ind: Vec<usize> = (0..n).map(|_| g.next() % m).collect();
                let gathered: Vec<R> = ind.iter().map(|&j| x[j]).collect();
                let want = dot(&gathered, &c, peel).to_f64().to_bits();
                assert_eq!(dot_gather(&x, &ind, &c, peel, &mut buf).to_f64().to_bits(), want, "dot_gather n {n} peel {peel}");
            }
        }
    }
}

#[test]
fn gathered_dot_matches_dot_f64() {
    gathered_dots::<f64>(false);
}

/// Products and sums that are NaN of either sign: the sign of a NaN meeting another is part of the bits.
#[test]
fn gathered_dot_keeps_the_nan_bits_of_dot_f64() {
    gathered_dots::<f64>(true);
}

#[test]
fn gathered_dot_matches_dot_f80() {
    gathered_dots::<crate::f80::F80>(false);
}

#[test]
fn gather_into_matches_the_collected_map() {
    let mut g = Rng(0x2545f4914f6cdd1d);
    let mut buf: Vec<f64> = vec![9.0; 3];
    for n in 0..80usize {
        let m = n + 5;
        let x: Vec<f64> = (0..m).map(|_| if g.next() % 7 == 0 { g.special() } else { g.value() }).collect();
        let ind: Vec<usize> = (0..n).map(|_| g.next() % m).collect();
        gather_into(&x, &ind, &mut buf);
        let want: Vec<f64> = ind.iter().map(|&j| x[j]).collect();
        assert_eq!(all_bits(&buf), all_bits(&want), "n {n}");
    }
}

#[test]
fn permuted_scatter_matches_the_plain_loop() {
    let mut g = Rng(0x9e3779b97f4a7c15);
    for m in [1usize, 5, 40, 300] {
        let mut perm: Vec<usize> = (0..m).collect();
        for i in (1..m).rev() {
            perm.swap(i, g.next() % (i + 1));
        }
        for n in [0, 1, 3, m, 2 * m] {
            // Repeated rows are summed in order.
            let rows: Vec<usize> = (0..n).map(|_| g.next() % m).collect();
            let vals: Vec<f64> = (0..n).map(|_| g.value()).collect();
            let mut want = vec![0.0f64; m];
            let mut want_start = Vec::new();
            for (&r, &v) in rows.iter().zip(&vals) {
                want[perm[r]] += v;
                want_start.push(perm[r]);
            }
            let mut x = vec![0.0f64; m];
            scatter_perm(&mut x, &perm, &rows, &vals);
            assert_eq!(all_bits(&x), all_bits(&want));
            let mut x = vec![0.0f64; m];
            let mut start = Vec::new();
            scatter_perm_list(&mut x, &perm, &rows, &vals, &mut start);
            assert_eq!(all_bits(&x), all_bits(&want));
            assert_eq!(start, want_start);
        }
    }
}

/// `LevelWorkspace::solve` with the checked accesses and the gathered eta copy (the code before the unchecked rewrite).
fn solve_old<R: Real>(
    ws: &mut LevelWorkspace<R>,
    f: &DualFactor<R>,
    levels: &[usize],
    maxlev: usize,
    rows: &[usize],
    vals: &[R],
) -> Option<(Vec<(usize, R)>, usize, f64)> {
    if !ws.prepare(f, maxlev + 1) {
        return None;
    }
    let mut low = maxlev + 1;
    let mut high = 0;
    for (&r, &v) in rows.iter().zip(vals) {
        let i = f.perm[r];
        ws.values[i] += v;
        ws.queue(i, levels[i]);
        low = low.min(levels[i]);
        high = high.max(levels[i]);
    }
    let mut at_tolerance = false;
    while low <= high {
        let mut j = 0;
        while j < ws.buckets[low].len() {
            let row = ws.buckets[low][j];
            j += 1;
            let x = ws.values[row];
            if x.abs() >= R::c(TOL) {
                ws.order.push(row);
                ws.state[row] = 2;
                at_tolerance |= x.abs() == R::c(TOL);
                let (dest, mult) = f.lcolumn(row);
                for (&to, &a) in dest.iter().zip(mult) {
                    let to = to as usize;
                    ws.queue(to, levels[to]);
                    high = high.max(levels[to]);
                    ws.values[to] -= a * x;
                }
            } else {
                ws.values[row] = R::ZERO;
            }
        }
        ws.buckets[low].clear();
        low += 1;
    }
    let count = ws.order.len();
    if !f.etas.is_empty() && (at_tolerance || count * f.dthr > f.m) {
        return None;
    }
    let mut work = ((f.nnz_l + f.m + 2 * f.m) as f64 * 5.0 * count as f64) / (f.m + 1) as f64 + 0.0;
    if !f.etas.is_empty() {
        let mut first = f.etas.len();
        for j in 0..count {
            let row = ws.order[j];
            first = first.min(f.first_eta[row]);
            ws.schedule(row, 0);
        }
        for word in first / 64..ws.pending.len() {
            while ws.pending[word] != 0 {
                let bit = ws.pending[word].trailing_zeros() as usize;
                ws.pending[word] &= ws.pending[word] - 1;
                let e = word * 64 + bit;
                let eta = &f.etas[e];
                let gather: Vec<R> = eta.ind.iter().map(|&row| ws.values[row]).collect();
                let v = dot(&gather, &eta.val, eta.peel);
                let row = eta.idx;
                if v.abs() <= R::c(TOL) {
                    ws.values[row] = R::ZERO;
                } else {
                    ws.values[row] = v;
                    if ws.state[row] != 2 {
                        if ws.state[row] == 0 {
                            ws.touched.push(row);
                        }
                        ws.state[row] = 2;
                        ws.order.push(row);
                        ws.schedule(row, e + 1);
                    }
                }
            }
        }
        let prefix = first.checked_sub(1).map_or(0, |i| ws.term_ends[i]);
        let terms = ws.term_ends.last().unwrap() - prefix;
        let ew = 2.0 * terms as f64 + (f.etas.len() - first) as f64 * 3.0;
        work = ew + 3.0 * ws.order.len() as f64 + work;
    }
    let out: Vec<(usize, R)> =
        ws.order.iter().copied().filter(|&i| ws.values[i].abs() > R::c(TOL)).map(|i| (i, ws.values[i])).collect();
    Some((out, count, work))
}

fn level_solves<R: Real>(seeds: &[u64]) -> (usize, usize) {
    let (mut solved, mut with_etas) = (0, 0);
    for &seed in seeds {
        let mut g = Rng(seed);
        for m in [9, 40, 130, 300] {
            let cols: Vec<_> = (0..m)
                .map(|i| match g.next() % 5 {
                    0 => (true, vec![i], vec![1.0]),
                    1 => (false, vec![i], vec![2.0]),
                    _ => {
                        let mut rows = vec![i];
                        let mut vals = vec![4.0];
                        for _ in 0..(1 + g.next() % 4) {
                            let r = g.next() % m;
                            if !rows.contains(&r) {
                                rows.push(r);
                                vals.push([0.5, -0.75, 1.25, 2.0, -1.0][g.next() % 5]);
                            }
                        }
                        (false, rows, vals)
                    }
                })
                .collect();
            let st = crate::lu_order::factor_structure_r::<R>(m, &cols);
            let mut f = DualFactor::<R>::fresh(m, m, &st, &cols);
            let (mut new_ws, mut old_ws) = (LevelWorkspace::<R>::default(), LevelWorkspace::<R>::default());
            for step in 0..m.min(48) {
                for _ in 0..4 {
                    let n = [1, 2, 3, 5, m / 8 + 1].map(|k: usize| k.min(m))[g.next() % 5];
                    let mut rows: Vec<usize> = Vec::new();
                    for _ in 0..n {
                        let r = g.next() % m;
                        if !rows.contains(&r) {
                            rows.push(r);
                        }
                    }
                    let vals: Vec<R> = rows.iter().map(|_| R::c(g.value())).collect();
                    let Some((levels, maxlev)) = f.levels.clone() else { continue };
                    let a = new_ws.solve(&f, &levels, maxlev, &rows, &vals, false);
                    let b = solve_old(&mut old_ws, &f, &levels, maxlev, &rows, &vals);
                    assert_eq!(a.is_some(), b.is_some(), "m {m} step {step}: route");
                    if let (Some(a), Some(b)) = (a, b) {
                        assert_eq!(bits(&a.0), bits(&b.0), "m {m} step {step}: list");
                        assert_eq!((a.1, a.2.to_bits()), (b.1, b.2.to_bits()), "m {m} step {step}: count, work");
                        solved += 1;
                        with_etas += usize::from(!f.etas.is_empty());
                    }
                    assert_eq!(all_bits(&new_ws.values), all_bits(&old_ws.values), "values");
                    assert_eq!((&new_ws.state, &new_ws.touched, &new_ws.order), (&old_ws.state, &old_ws.touched, &old_ws.order));
                    assert_eq!(new_ws.pending, old_ws.pending);
                }
                if step % 3 == 2 {
                    let r = g.next() % m;
                    let mut rows = vec![r];
                    let mut vals = vec![R::c(5.0 + step as f64)];
                    for _ in 0..(1 + g.next() % 6) {
                        let q = g.next() % m;
                        if !rows.contains(&q) {
                            rows.push(q);
                            vals.push(R::c([0.25, -0.5, 1.5][g.next() % 3]));
                        }
                    }
                    f.ftran(&rows, &vals);
                    f.update(r);
                }
            }
        }
    }
    (solved, with_etas)
}

#[test]
fn level_solve_matches_the_checked_solve_f64() {
    let (solved, with_etas) = level_solves::<f64>(&[0x2545f4914f6cdd1d, 0x9e3779b97f4a7c15, 0x43e941132abcu64, 0xdeadbeefcafe1234, 0x1234567, 0xabcdef12345, 0x77777777]);
    assert!(solved > 200 && with_etas > 40, "solved {solved}, with etas {with_etas}");
}
