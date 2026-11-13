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
    level_sweep_two_loops(lt, lev, maxlev, lo, x, marks, stamp, buckets, res);
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

/// A repeated destination in an `lt` line: the guard must route to the two-loop code (bit for bit), where the fused loop differs.
#[test]
fn repeated_destination_uses_two_loops() {
    // node 2 subtracts a then b into node 0: (2 - a) - b = 1 - 2^-53 but (2 - b) - a = 1 (2 - 2^-53 rounds to 2).
    let (a, b) = (1.0f64, 2f64.powi(-53));
    let mut lt = Csr::<f64>::with_capacity(3, 0);
    lt.push_line(std::iter::empty());
    lt.push_line(std::iter::empty());
    lt.push_line([(0usize, a), (0usize, b)].into_iter());
    assert!(!lt_lines_distinct(&lt));
    let lev = vec![1usize, 0, 0];
    let maxlev = 1;
    let run = |mode: u8| {
        let mut x = vec![2.0f64, 0.0, 1.0];
        let mut marks = vec![0u32; 3];
        let mut buckets: Vec<Vec<usize>> = vec![Vec::new(); 2];
        let mut res = Vec::new();
        marks[2] = 1;
        buckets[0].push(2);
        match mode {
            0 => level_sweep_guarded(&lt, lt_lines_distinct(&lt), &lev, maxlev, 0, &mut x, &mut marks, 1, &mut buckets, &mut res),
            1 => level_sweep_old(&lt, &lev, maxlev, 0, &mut x, &mut marks, 1, &mut buckets, &mut res),
            _ => level_sweep(&lt, &lev, maxlev, 0, &mut x, &mut marks, 1, &mut buckets, &mut res),
        }
        (x.iter().map(|v| v.to_bits()).collect::<Vec<_>>(), marks, res)
    };
    assert!(run(0) == run(1));
    assert!(run(0).0 != run(2).0, "fused loop expected to differ on a repeated destination");
}
