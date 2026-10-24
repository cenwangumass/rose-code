//! The reference's doubleton pairing pass (`0xdf800`, after the singleton phases when 50 or more rows remain): pairs of rows
//! with the same two active columns are eliminated as 2x2 blocks, then pairs of columns with the same two active rows, with
//! the singletons that each elimination creates drained right away (`reverse-engineer/specs/lu-factor.md`).

use std::collections::HashMap;

const SING_TOL: f64 = 1e-10;
const DROP_TOL: f64 = 1e-13;
const UNZ_TOL: f64 = 1e-13;

/// What the pass computes that the generic singleton formulas (`factor_structure`) do not: the Schur pivot of each pair, and
/// for a row pair the second column's L multipliers.
#[derive(Default)]
pub struct Overrides {
    pub diag: HashMap<usize, f64>,
    pub l: HashMap<usize, (Vec<usize>, Vec<f64>)>,
    /// Number of nonzero entries in the U row of a column pair's second pivot (the pivot row `hi` after the elimination,
    /// fill from the first pivot's row included).
    pub unz: HashMap<usize, usize>,
    /// Columns of the U row of a column pair's two pivots in the reference's storage order (it fixes the order of the
    /// factor update's row lists): the first row lists its other columns in row order and the partner last; the second row
    /// lists its own columns, then the fill that the first row adds.
    pub urow: HashMap<usize, Vec<usize>>,
}

struct Pass<'a> {
    m: usize,
    cols: &'a [(bool, Vec<usize>, Vec<f64>)],
    /// Active counts (> 0 while active, -1 or 0 once assigned).
    rowc: Vec<i64>,
    colc: Vec<i64>,
    /// Row-wise copies (column positions and values) of the rows left after phase 3, in the reference's array order.
    ridx: Vec<Vec<usize>>,
    rval: Vec<Vec<f64>>,
    temp: Vec<f64>,
    /// Hash entries: owner id and key; the key's pending entry waits for a partner.
    ids: Vec<usize>,
    keys: Vec<usize>,
    pending: HashMap<usize, usize>,
    pairs: Vec<(usize, usize)>,
    seq: Vec<(usize, usize)>,
    phase: Vec<u8>,
    ov: Overrides,
}

impl Pass<'_> {
    fn register(&mut self, id: usize, key: usize) {
        let e = self.ids.len();
        self.ids.push(id);
        self.keys.push(key);
        if let Some(prev) = self.pending.remove(&key) {
            self.pairs.push((prev, e));
        } else {
            self.pending.insert(key, e);
        }
    }

    /// A row that reached two active columns: compacts its active entries to the front, orders the last two by column count
    /// (smaller first, ties by index) and files the pair.
    fn register_row(&mut self, r: usize) {
        let mut k = 0;
        for j in 0..self.ridx[r].len() {
            if self.colc[self.ridx[r][j]] > 0 {
                if k < j {
                    self.ridx[r].swap(k, j);
                    self.rval[r][k] = self.rval[r][j];
                }
                k += 1;
            }
        }
        if k < 2 {
            return;
        }
        let (last, prev) = (self.ridx[r][k - 1], self.ridx[r][k - 2]);
        let (cl, cp) = (self.colc[last], self.colc[prev]);
        let key = if cl < cp || (cl == cp && last < prev) {
            self.ridx[r].swap(k - 2, k - 1);
            self.rval[r].swap(k - 2, k - 1);
            last + self.m * prev
        } else {
            prev + self.m * last
        };
        self.register(r, key);
    }

    /// A column that reached two active rows: the first two active rows in column order, smaller count first.
    fn register_col(&mut self, c: usize) {
        let mut first = usize::MAX;
        let mut second = usize::MAX;
        for &r in &self.cols[c].1 {
            if self.rowc[r] > 0 {
                if first != usize::MAX {
                    second = r;
                    break;
                }
                first = r;
            }
        }
        if second == usize::MAX {
            return;
        }
        let (cs, cf) = (self.rowc[second], self.rowc[first]);
        let key = if cs < cf || (cs == cf && second < first) {
            second + self.m * first
        } else {
            first + self.m * second
        };
        self.register(c, key);
    }

    fn push_step(&mut self, row: usize, col: usize) -> usize {
        self.seq.push((row, col));
        self.phase.push(4);
        self.seq.len() - 1
    }

    /// A row with one active column left: pivot there; the column's other active rows lose one entry.
    fn row_singleton(&mut self, r: usize, rsing: &mut Vec<usize>) {
        let Some(j) = (0..self.ridx[r].len()).find(|&j| self.colc[self.ridx[r][j]] >= 1) else {
            return;
        };
        let c = self.ridx[r][j];
        self.rowc[r] = 0;
        self.colc[c] = 0;
        self.push_step(r, c);
        let col = self.cols[c].1.clone();
        for &rr in &col {
            if self.rowc[rr] > 0 {
                self.rowc[rr] -= 1;
                if self.rowc[rr] == 1 {
                    rsing.push(rr);
                } else if self.rowc[rr] == 2 {
                    self.register_row(rr);
                }
            }
        }
    }

    fn col_singleton(&mut self, c: usize, csing: &mut Vec<usize>) {
        let Some(&r) = self.cols[c].1.iter().find(|&&r| self.rowc[r] >= 1) else {
            return;
        };
        self.colc[c] = 0;
        self.rowc[r] = 0;
        self.push_step(r, c);
        for j in 0..self.ridx[r].len() {
            let c2 = self.ridx[r][j];
            if self.colc[c2] > 0 {
                self.colc[c2] -= 1;
                if self.colc[c2] == 1 {
                    csing.push(c2);
                } else if self.colc[c2] == 2 {
                    self.register_col(c2);
                }
            }
        }
    }

    /// Two rows with the same two active columns: pivot on the larger entry of the first column, then on the Schur element.
    fn row_pair(&mut self, r1: usize, r2: usize) {
        let (c0, c1) = (self.ridx[r1][0], self.ridx[r1][1]);
        let (v0, v1) = (self.rval[r1][0], self.rval[r1][1]);
        let (w0, w1) = (self.rval[r2][0], self.rval[r2][1]);
        let (prow, p, u12, orow, o0, o1) = if v0.abs() < w0.abs() {
            (r2, w0, w1, r1, v0, v1)
        } else {
            (r1, v0, v1, r2, w0, w1)
        };
        let schur = o1 - (o0 * u12) / p;
        if schur.abs() < SING_TOL {
            return;
        }
        let mut rsing: Vec<usize> = Vec::new();
        self.rowc[prow] = 0;
        self.colc[c0] = 0;
        self.push_step(prow, c0);
        let mut l1: Vec<usize> = Vec::new();
        let col0 = self.cols[c0].clone();
        for (&r, &v) in col0.1.iter().zip(&col0.2) {
            if self.rowc[r] > 0 {
                self.temp[r] = ((-v) * u12) / p;
                l1.push(r);
                self.rowc[r] -= 1;
                if self.rowc[r] == 1 {
                    rsing.push(r);
                } else if self.rowc[r] == 2 {
                    self.register_row(r);
                }
            }
        }
        self.rowc[orow] = 0;
        self.colc[c1] = 0;
        let s2 = self.push_step(orow, c1);
        self.temp[orow] = 0.0;
        let mut lr: Vec<usize> = Vec::new();
        let mut lv: Vec<f64> = Vec::new();
        let col1 = self.cols[c1].clone();
        for (&r, &v) in col1.1.iter().zip(&col1.2) {
            if self.rowc[r] > 0 {
                let x = (v + self.temp[r]) / schur;
                self.temp[r] = 0.0;
                if x.abs() > DROP_TOL {
                    lr.push(r);
                    lv.push(x);
                }
                self.rowc[r] -= 1;
                if self.rowc[r] == 1 {
                    rsing.push(r);
                } else if self.rowc[r] == 2 {
                    self.register_row(r);
                }
            }
        }
        for &r in &l1 {
            let x = self.temp[r] / schur;
            if x.abs() > DROP_TOL {
                lr.push(r);
                lv.push(self.temp[r] / schur);
            }
            self.temp[r] = 0.0;
        }
        self.ov.diag.insert(s2, schur);
        self.ov.l.insert(s2, (lr, lv));
        while let Some(r) = rsing.pop() {
            if self.rowc[r] > 0 {
                self.row_singleton(r, &mut rsing);
            }
        }
    }

    /// Two columns with the same two active rows: pivot the row with the smaller count on its larger entry, then the other row on
    /// the Schur element.
    fn col_pair(&mut self, d1: usize, d2: usize, key: usize) {
        let (lo, hi) = (key % self.m, key / self.m);
        let mut v = [0.0f64; 2];
        let mut others_lo: Vec<(usize, f64)> = Vec::new();
        for j in 0..self.ridx[lo].len() {
            let c = self.ridx[lo][j];
            if self.colc[c] > 0 {
                if c == d1 {
                    v[0] = self.rval[lo][j];
                } else if c == d2 {
                    v[1] = self.rval[lo][j];
                } else {
                    others_lo.push((c, self.rval[lo][j]));
                }
            }
        }
        let (pc, p, nc, u12) = if v[0].abs() < v[1].abs() {
            (d2, v[1], d1, v[0])
        } else {
            (d1, v[0], d2, v[1])
        };
        let mut hv = [0.0f64; 2];
        let mut others_hi: Vec<(usize, f64)> = Vec::new();
        for j in 0..self.ridx[hi].len() {
            let c = self.ridx[hi][j];
            if self.colc[c] > 0 {
                if c == pc {
                    hv[0] = self.rval[hi][j];
                } else if c == nc {
                    hv[1] = self.rval[hi][j];
                } else {
                    others_hi.push((c, self.rval[hi][j]));
                }
            }
        }
        let (x_p, x_n) = (hv[0], hv[1]);
        let schur = x_n - (u12 * x_p) / p;
        if schur.abs() < SING_TOL {
            return;
        }
        let mut csing: Vec<usize> = Vec::new();
        self.rowc[lo] = 0;
        self.colc[pc] = 0;
        let s1 = self.push_step(lo, pc);
        self.ov.urow.insert(
            s1,
            others_lo
                .iter()
                .map(|e| e.0)
                .chain(std::iter::once(nc))
                .collect(),
        );
        for &(c, u) in &others_lo {
            self.temp[c] = ((-x_p) * u) / p;
            self.colc[c] -= 1;
            if self.colc[c] == 1 {
                csing.push(c);
            } else if self.colc[c] == 2 {
                self.register_col(c);
            }
        }
        self.rowc[hi] = 0;
        self.colc[nc] = 0;
        let s2 = self.push_step(hi, nc);
        self.ov.diag.insert(s2, schur);
        self.ov.l.insert(s2, (Vec::new(), Vec::new()));
        {
            let mut vals: HashMap<usize, f64> = HashMap::new();
            for &(c, v) in &others_hi {
                vals.insert(c, v);
            }
            for &(c, _) in &others_lo {
                let t = self.temp[c];
                let e = vals.entry(c).or_insert(0.0);
                *e += t;
            }
            self.ov
                .unz
                .insert(s2, vals.values().filter(|v| v.abs() > UNZ_TOL).count());
            let own: Vec<usize> = others_hi
                .iter()
                .map(|e| e.0)
                .filter(|c| vals[c].abs() > UNZ_TOL)
                .collect();
            let fill = others_lo
                .iter()
                .map(|e| e.0)
                .filter(|c| !others_hi.iter().any(|h| h.0 == *c) && vals[c].abs() > UNZ_TOL);
            self.ov
                .urow
                .insert(s2, own.into_iter().chain(fill).collect());
        }
        for &(c, _) in &others_hi {
            self.colc[c] -= 1;
            if self.colc[c] == 1 {
                csing.push(c);
            } else if self.colc[c] == 2 {
                self.register_col(c);
            }
            self.temp[c] = 0.0;
        }
        for &(c, _) in &others_lo {
            self.temp[c] = 0.0;
        }
        while let Some(c) = csing.pop() {
            if self.colc[c] > 0 {
                self.col_singleton(c, &mut csing);
            }
        }
    }
}

/// Runs the pass after phase 4. `done3` are the rows and columns assigned by phases 1-3 (the row lists are built from what was
/// left then), `row_done`/`col_done` the state after phase 4. Appends the pass's pivots to `seq` and `phase`.
pub fn run(
    m: usize,
    cols: &[(bool, Vec<usize>, Vec<f64>)],
    seq: &mut Vec<(usize, usize)>,
    phase: &mut Vec<u8>,
    row_done3: &[bool],
    col_done3: &[bool],
    row_done: &[bool],
    col_done: &[bool],
) -> Overrides {
    let n = cols.len();
    let mut ridx: Vec<Vec<usize>> = vec![Vec::new(); m];
    let mut rval: Vec<Vec<f64>> = vec![Vec::new(); m];
    for pos in 0..n {
        if col_done3[pos] {
            continue;
        }
        for (&r, &v) in cols[pos].1.iter().zip(&cols[pos].2) {
            if !row_done3[r] {
                ridx[r].push(pos);
                rval[r].push(v);
            }
        }
    }
    let mut rowc = vec![-1i64; m];
    let mut colc = vec![-1i64; n];
    for pos in 0..n {
        if !col_done[pos] {
            colc[pos] = cols[pos].1.iter().filter(|&&r| !row_done[r]).count() as i64;
            for &r in &cols[pos].1 {
                if !row_done[r] {
                    if rowc[r] < 0 {
                        rowc[r] = 0;
                    }
                    rowc[r] += 1;
                }
            }
        }
    }
    let mut st = Pass {
        m,
        cols,
        rowc,
        colc,
        ridx,
        rval,
        temp: vec![0.0; m.max(n)],
        ids: Vec::new(),
        keys: Vec::new(),
        pending: HashMap::new(),
        pairs: Vec::new(),
        seq: std::mem::take(seq),
        phase: std::mem::take(phase),
        ov: Overrides::default(),
    };
    // Row pairs.
    for r in 0..m {
        if !row_done[r] && st.rowc[r] == 2 {
            st.register_row(r);
        }
    }
    while let Some((e1, e2)) = st.pairs.pop() {
        let (r1, r2) = (st.ids[e1], st.ids[e2]);
        if st.rowc[r1] > 1 && st.rowc[r2] > 1 {
            st.row_pair(r1, r2);
        }
    }
    // Column pairs.
    st.pending.clear();
    for c in 0..n {
        if st.colc[c] == 2 {
            st.register_col(c);
        }
    }
    while let Some((e1, e2)) = st.pairs.pop() {
        let (d1, d2) = (st.ids[e1], st.ids[e2]);
        if st.colc[d1] > 1 && st.colc[d2] > 1 {
            let key = st.keys[e1];
            st.col_pair(d1, d2, key);
        }
    }
    *seq = st.seq;
    *phase = st.phase;
    st.ov
}
