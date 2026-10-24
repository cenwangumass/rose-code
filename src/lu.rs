//! Sparse LU factorization of a basis matrix (left-looking, threshold partial pivoting) with product-form updates.

const DROP: f64 = 1e-14;
const NONE: usize = usize::MAX;

/// A sparse column given as `(row, value)` pairs.
pub type Column = Vec<(usize, f64)>;

pub struct Factor {
    m: usize,
    /// Position of the basis column eliminated at each step.
    col_order: Vec<usize>,
    /// Row chosen as pivot at each step, and its inverse.
    pivot_row: Vec<usize>,
    step_of_row: Vec<usize>,
    /// Multipliers of step `k`: `(row, l)`.
    l: Vec<Vec<(usize, f64)>>,
    /// Entries of `U` column `k` above the diagonal: `(step, value)`.
    u: Vec<Vec<(usize, f64)>>,
    diag: Vec<f64>,
    /// Product-form updates: pivot position and the full alpha column (by position).
    etas: Vec<(usize, Vec<f64>)>,
    pub nnz_lu: usize,
    pub nnz_eta: usize,
}

impl Factor {
    /// Factor the matrix whose column at position `p` is `cols[p]`. Returns the singular position on failure.
    pub fn new<C: std::borrow::Borrow<Column>>(m: usize, cols: &[C]) -> Result<Factor, usize> {
        let cols: Vec<&Column> = cols.iter().map(|c| c.borrow()).collect();
        // Eliminate sparser columns first (slack singletons first).
        let mut col_order: Vec<usize> = (0..m).collect();
        col_order.sort_by_key(|&p| (cols[p].len(), p));
        let mut step_of_row = vec![NONE; m];
        let mut pivot_row = Vec::with_capacity(m);
        let mut l: Vec<Vec<(usize, f64)>> = Vec::with_capacity(m);
        let mut u: Vec<Vec<(usize, f64)>> = Vec::with_capacity(m);
        let mut diag = Vec::with_capacity(m);
        let mut row_count = vec![0usize; m];
        for c in &cols {
            for &(r, _) in c.iter() {
                row_count[r] += 1;
            }
        }
        let mut w = vec![0.0f64; m];
        // Rows whose `w` may be nonzero (marker plus list), visited in ascending order: the scans below see the same rows in the
        // same order as full scans over `0..m` would, minus the ones that cannot contribute.
        let mut touched_mark = vec![false; m];
        let mut touched: Vec<usize> = Vec::new();
        // Pivot steps whose pivot row is touched and not yet processed (min-heap): ascending step order as in the dense loop.
        let mut step_mark = vec![false; m];
        let mut heap: std::collections::BinaryHeap<std::cmp::Reverse<usize>> =
            std::collections::BinaryHeap::new();
        for k in 0..m {
            let p = col_order[k];
            macro_rules! touch {
                ($r:expr) => {{
                    let r = $r;
                    if !touched_mark[r] {
                        touched_mark[r] = true;
                        touched.push(r);
                        let s = step_of_row[r];
                        if s != NONE && !step_mark[s] {
                            step_mark[s] = true;
                            heap.push(std::cmp::Reverse(s));
                        }
                    }
                }};
            }
            for &(r, a) in cols[p] {
                w[r] += a;
                touch!(r);
            }
            let mut ucol = Vec::new();
            while let Some(std::cmp::Reverse(j)) = heap.pop() {
                step_mark[j] = false;
                let t = w[pivot_row[j]];
                if t != 0.0 {
                    ucol.push((j, t));
                    for &(r, lv) in &l[j] {
                        w[r] -= lv * t;
                        touch!(r);
                    }
                }
            }
            touched.sort_unstable();
            // Threshold pivoting: among candidates within 0.1 of the maximum prefer the sparsest row.
            let mut maxabs = 0.0f64;
            for &r in &touched {
                if step_of_row[r] == NONE {
                    maxabs = maxabs.max(w[r].abs());
                }
            }
            if maxabs < 1e-11 {
                return Err(p);
            }
            let mut best = NONE;
            for &r in &touched {
                if step_of_row[r] == NONE && w[r].abs() >= 0.1 * maxabs {
                    if best == NONE
                        || row_count[r] < row_count[best]
                        || (row_count[r] == row_count[best] && w[r].abs() > w[best].abs())
                    {
                        best = r;
                    }
                }
            }
            let d = w[best];
            let mut lcol = Vec::new();
            for &r in &touched {
                if step_of_row[r] == NONE && r != best && w[r].abs() > DROP {
                    lcol.push((r, w[r] / d));
                }
            }
            for &r in &touched {
                w[r] = 0.0;
                touched_mark[r] = false;
            }
            touched.clear();
            step_of_row[best] = k;
            pivot_row.push(best);
            l.push(lcol);
            u.push(ucol);
            diag.push(d);
        }
        let nnz_lu =
            l.iter().map(Vec::len).sum::<usize>() + u.iter().map(Vec::len).sum::<usize>() + m;
        Ok(Factor {
            m,
            col_order,
            pivot_row,
            step_of_row,
            l,
            u,
            diag,
            etas: Vec::new(),
            nnz_lu,
            nnz_eta: 0,
        })
    }

    pub fn num_updates(&self) -> usize {
        self.etas.len()
    }

    /// Solve `B x = b` in place (`b` indexed by row on entry, by basis position on exit).
    pub fn ftran(&self, b: &mut Vec<f64>) {
        let m = self.m;
        for j in 0..m {
            let t = b[self.pivot_row[j]];
            if t != 0.0 {
                for &(r, lv) in &self.l[j] {
                    b[r] -= lv * t;
                }
            }
        }
        let mut x = vec![0.0; m];
        for k in (0..m).rev() {
            let xk = b[self.pivot_row[k]] / self.diag[k];
            x[self.col_order[k]] = xk;
            if xk != 0.0 {
                for &(j, uv) in &self.u[k] {
                    b[self.pivot_row[j]] -= uv * xk;
                }
            }
        }
        for (r, alpha) in &self.etas {
            let t = x[*r] / alpha[*r];
            if t != 0.0 {
                for (i, a) in alpha.iter().enumerate() {
                    if i != *r && *a != 0.0 {
                        x[i] -= a * t;
                    }
                }
            }
            x[*r] = t;
        }
        *b = x;
    }

    /// Solve `B' z = c` in place (`c` indexed by basis position on entry, by row on exit).
    pub fn btran(&self, c: &mut Vec<f64>) {
        let m = self.m;
        for (r, alpha) in self.etas.iter().rev() {
            let mut s = c[*r];
            for (i, a) in alpha.iter().enumerate() {
                if i != *r && *a != 0.0 {
                    s -= a * c[i];
                }
            }
            c[*r] = s / alpha[*r];
        }
        let mut v = vec![0.0; m];
        for k in 0..m {
            let mut s = c[self.col_order[k]];
            for &(j, uv) in &self.u[k] {
                s -= uv * v[j];
            }
            v[k] = s / self.diag[k];
        }
        for j in (0..m).rev() {
            let mut s = v[j];
            for &(r, lv) in &self.l[j] {
                s -= lv * v[self.step_of_row[r]];
            }
            v[j] = s;
        }
        let mut z = vec![0.0; m];
        for k in 0..m {
            z[self.pivot_row[k]] = v[k];
        }
        *c = z;
    }

    /// Product-form update: the column at position `r` is replaced; `alpha = B^-1 a_q` (by position).
    pub fn update(&mut self, r: usize, alpha: Vec<f64>) {
        self.nnz_eta += alpha.iter().filter(|a| **a != 0.0).count();
        self.etas.push((r, alpha));
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn dense_cols(a: &[&[f64]]) -> Vec<Column> {
        let m = a.len();
        (0..m)
            .map(|c| {
                (0..m)
                    .filter(|&r| a[r][c] != 0.0)
                    .map(|r| (r, a[r][c]))
                    .collect()
            })
            .collect()
    }

    #[test]
    fn solves_and_updates_against_dense_residuals() {
        let a: [&[f64]; 4] = [
            &[2.0, 0.0, 1.0, 0.0],
            &[1.0, 3.0, 0.0, 0.5],
            &[0.0, 1.0, 4.0, 0.0],
            &[0.0, 0.0, 1.0, 5.0],
        ];
        let cols = dense_cols(&a);
        let mut f = Factor::new(4, &cols).unwrap();
        let b = vec![1.0, 2.0, 3.0, 4.0];
        let mut x = b.clone();
        f.ftran(&mut x);
        for r in 0..4 {
            let lhs: f64 = (0..4).map(|c| a[r][c] * x[c]).sum();
            assert!((lhs - b[r]).abs() < 1e-12);
        }
        let mut z = b.clone();
        f.btran(&mut z);
        for c in 0..4 {
            let lhs: f64 = (0..4).map(|r| a[r][c] * z[r]).sum();
            assert!((lhs - b[c]).abs() < 1e-12);
        }
        // Replace position 1 by a new column and compare with a fresh factorization.
        let newcol = vec![1.0, 1.0, 2.0, 1.0];
        let mut alpha = newcol.clone();
        f.ftran(&mut alpha);
        f.update(1, alpha);
        let mut cols2 = cols.clone();
        cols2[1] = newcol.iter().enumerate().map(|(r, v)| (r, *v)).collect();
        let g = Factor::new(4, &cols2).unwrap();
        let (mut x1, mut x2) = (b.clone(), b.clone());
        f.ftran(&mut x1);
        g.ftran(&mut x2);
        let (mut z1, mut z2) = (b.clone(), b.clone());
        f.btran(&mut z1);
        g.btran(&mut z2);
        for i in 0..4 {
            assert!((x1[i] - x2[i]).abs() < 1e-12, "ftran {i}");
            assert!((z1[i] - z2[i]).abs() < 1e-12, "btran {i}");
        }
    }

    #[test]
    fn reports_singular_position() {
        let cols: Vec<Column> = vec![vec![(0, 1.0)], vec![(0, 2.0)]];
        assert!(Factor::new(2, &cols).is_err());
    }
}
