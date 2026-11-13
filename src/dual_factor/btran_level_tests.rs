//! Differential test of the fused `level_sweep` against the previous two-loop level machinery.
use super::*;

struct Rng(u64);
impl Rng {
    fn next(&mut self) -> u64 {
        self.0 ^= self.0 << 13;
        self.0 ^= self.0 >> 7;
        self.0 ^= self.0 << 17;
        self.0
    }
    fn below(&mut self, n: usize) -> usize {
        (self.next() % n as u64) as usize
    }
}

/// The previous loop: scatter forward, then mark and push from the last entry.
#[allow(clippy::too_many_arguments)]
fn level_sweep_old(
    lt: &Csr<f64>,
    lev: &[usize],
    maxlev: usize,
    lo: usize,
    x: &mut [f64],
    marks: &mut [u32],
    stamp: u32,
    buckets: &mut [Vec<usize>],
    res: &mut Vec<usize>,
) {
    let mut lv_i = lo;
    while lv_i <= maxlev {
        let mut k = 0;
        while k < buckets[lv_i].len() {
            let i = buckets[lv_i][k];
            k += 1;
            let v = x[i];
            if v.abs() >= TOL {
                res.push(i);
                let (src, mult) = lt.col(i);
                for (&j, &lv) in src.iter().zip(mult) {
                    x[j as usize] -= lv * v;
                }
                for &j in src.iter().rev() {
                    let j = j as usize;
                    if marks[j] != stamp {
                        marks[j] = stamp;
                        buckets[lev[j]].push(j);
                    }
                }
            } else {
                x[i] = 0.0;
            }
        }
        buckets[lv_i].clear();
        lv_i += 1;
    }
}

#[test]
fn level_sweep_matches_two_loops() {
    let mut rng = Rng(0x1234_5678_9abc_def1);
    for m in [1usize, 2, 5, 40, 200] {
        for round in 0..60 {
            let mut lt = Csr::<f64>::with_capacity(m, 0);
            let dens = 1 + rng.below(6);
            for v in 0..m {
                let mut ys: Vec<usize> = Vec::new();
                if v > 0 {
                    for _ in 0..rng.below(dens + 1) {
                        let y = rng.below(v);
                        if !ys.contains(&y) {
                            ys.push(y);
                        }
                    }
                }
                lt.push_line(ys.into_iter().map(|y| {
                    let a = match rng.below(6) {
                        0 => 0.0,
                        1 => -0.0,
                        _ => (rng.next() % 1000) as f64 / 13.0 - 38.0,
                    };
                    (y, a)
                }));
            }
            let mut lev = vec![0usize; m];
            for u in (0..m).rev() {
                for &y in lt.col(u).0 {
                    lev[y as usize] = lev[y as usize].max(lev[u] + 1);
                }
            }
            let maxlev = lev.iter().copied().max().unwrap();
            let mut x0 = vec![0.0f64; m];
            let mut seeds = Vec::new();
            for i in (0..m).rev() {
                if rng.below(4) == 0 {
                    seeds.push(i);
                    x0[i] = match rng.below(6) {
                        0 => 1e-12,
                        1 => f64::NAN,
                        2 => -0.0,
                        _ => (rng.next() % 100) as f64 / 3.0 - 16.0,
                    };
                }
            }
            if seeds.is_empty() {
                continue;
            }
            let stamp = 1 + round as u32;
            let run = |new: bool| {
                let mut x = x0.clone();
                let mut marks = vec![0u32; m];
                let mut buckets: Vec<Vec<usize>> = vec![Vec::new(); maxlev + 1];
                let mut res = Vec::new();
                let mut lo = usize::MAX;
                for &i in &seeds {
                    marks[i] = stamp;
                    buckets[lev[i]].push(i);
                    lo = lo.min(lev[i]);
                }
                if new {
                    level_sweep(&lt, &lev, maxlev, lo, &mut x, &mut marks, stamp, &mut buckets, &mut res);
                } else {
                    level_sweep_old(&lt, &lev, maxlev, lo, &mut x, &mut marks, stamp, &mut buckets, &mut res);
                }
                assert!(buckets.iter().all(|b| b.is_empty()));
                (x.iter().map(|v| v.to_bits()).collect::<Vec<_>>(), marks, res)
            };
            assert!(run(true) == run(false), "m {m} round {round}");
        }
    }
}
