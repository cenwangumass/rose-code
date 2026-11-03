//! Sparse LU factorization of a basis matrix (left-looking, threshold partial pivoting) with product-form updates.

const DROP: f64 = 1e-14;
const NONE: usize = usize::MAX;

/// A sparse column given as `(row, value)` pairs.
pub type Column = Vec<(usize, f64)>;

enum Eta {
    Dense(usize, Vec<f64>),
    Sparse {
        r: usize,
        diagonal: f64,
        entries: Vec<(usize, f64)>,
    },
}

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
    /// Product-form updates, with sparse entries in ascending basis position.
    etas: Vec<Eta>,
    pub nnz_lu: usize,
    pub nnz_eta: usize,
}

impl Factor {
    /// Factor the matrix whose column at position `p` is `cols[p]`. Returns the singular position on failure.
    pub fn new<C: std::borrow::Borrow<Column>>(m: usize, cols: &[C]) -> Result<Factor, usize> {
        Self::with_floor(m, cols, 1e-11)
    }

    /// As `new`, with the absolute floor below which a column without an acceptable pivot is singular (the dual decides singularity with the reference's
    /// structure mirror and factors an ill-conditioned basis anyway: `fx_midsing_1390` has a diagonal entry of 5e-13).
    pub fn with_floor<C: std::borrow::Borrow<Column>>(
        m: usize,
        cols: &[C],
        floor: f64,
    ) -> Result<Factor, usize> {
        Self::factor(m, cols, floor, None)
    }

    /// Factor columns already in an accepted pivot order using the corresponding
    /// pivot rows. Singularity has been decided by the reference factorization;
    /// only an exactly zero auxiliary pivot is rejected here.
    pub fn with_pivot_rows<C: std::borrow::Borrow<Column>>(
        m: usize,
        cols: &[C],
        rows: &[usize],
    ) -> Result<Factor, usize> {
        Self::factor(m, cols, 0.0, Some(rows))
    }

    fn factor<C: std::borrow::Borrow<Column>>(
        m: usize,
        cols: &[C],
        floor: f64,
        rows: Option<&[usize]>,
    ) -> Result<Factor, usize> {
        let cols: Vec<&Column> = cols.iter().map(|c| c.borrow()).collect();
        // Eliminate sparser columns first (slack singletons first).
        let mut col_order: Vec<usize> = (0..m).collect();
        if rows.is_none() {
            col_order.sort_by_key(|&p| (cols[p].len(), p));
        }
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
            if maxabs < floor || maxabs == 0.0 {
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
            if let Some(rows) = rows {
                best = rows[k];
                if w[best] == 0.0 {
                    return Err(p);
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
        for eta in &self.etas {
            match eta {
                Eta::Dense(r, alpha) => {
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
                Eta::Sparse {
                    r,
                    diagonal,
                    entries,
                } => {
                    let t = x[*r] / diagonal;
                    if t != 0.0 {
                        for (i, a) in entries {
                            x[*i] -= a * t;
                        }
                    }
                    x[*r] = t;
                }
            }
        }
        *b = x;
    }

    /// Solve `B' z = c` in place (`c` indexed by basis position on entry, by row on exit).
    pub fn btran(&self, c: &mut Vec<f64>) {
        let m = self.m;
        for eta in self.etas.iter().rev() {
            match eta {
                Eta::Dense(r, alpha) => {
                    let mut s = c[*r];
                    for (i, a) in alpha.iter().enumerate() {
                        if i != *r && *a != 0.0 {
                            s -= a * c[i];
                        }
                    }
                    c[*r] = s / alpha[*r];
                }
                Eta::Sparse {
                    r,
                    diagonal,
                    entries,
                } => {
                    let mut s = c[*r];
                    for (i, a) in entries {
                        s -= a * c[*i];
                    }
                    c[*r] = s / diagonal;
                }
            }
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

    /// Accept the raw FTRAN column before pivot removal or tolerance filtering.
    pub(crate) fn update_column<R: crate::real::Real>(
        &mut self,
        r: usize,
        column: Vec<(usize, R)>,
    ) {
        let mut column: Vec<(usize, f64)> =
            column.into_iter().map(|(i, v)| (i, v.to_f64())).collect();
        // The dense kernel fixes NaN-payload selection through its multiply operand order.
        if column.len() >= self.m.div_ceil(2) || column.iter().any(|e| !e.1.is_finite()) {
            let mut alpha = vec![0.0; self.m];
            for (i, v) in column {
                alpha[i] = v;
            }
            self.update(r, alpha);
            return;
        }
        if !column.windows(2).all(|w| w[0].0 < w[1].0) {
            column.sort_by_key(|e| e.0);
            // Dense staging keeps the last occurrence of a repeated position.
            column.dedup_by(|later, earlier| {
                if later.0 == earlier.0 {
                    *earlier = *later;
                    true
                } else {
                    false
                }
            });
        }
        let mut diagonal = 0.0;
        column.retain(|&(i, a)| {
            if i == r {
                diagonal = a;
                false
            } else {
                a != 0.0
            }
        });
        self.nnz_eta += column.len() + usize::from(diagonal != 0.0);
        self.etas.push(Eta::Sparse {
            r,
            diagonal,
            entries: column,
        });
    }

    /// Product-form update: the column at position `r` is replaced; `alpha = B^-1 a_q` (by position).
    pub fn update(&mut self, r: usize, alpha: Vec<f64>) {
        self.nnz_eta += alpha.iter().filter(|a| **a != 0.0).count();
        self.etas.push(Eta::Dense(r, alpha));
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

#[cfg(test)]
impl Factor {
    fn frozen_ftran(&self, b: &mut Vec<f64>, etas: &[(usize, Vec<f64>)]) {
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
        for (r, alpha) in etas {
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

    fn frozen_btran(&self, c: &mut Vec<f64>, etas: &[(usize, Vec<f64>)]) {
        let m = self.m;
        for (r, alpha) in etas.iter().rev() {
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
}

#[cfg(test)]
mod eta_storage_tests {
    use super::*;

    fn check(m: usize, etas: &[(usize, Vec<f64>)], rhs: &[f64]) {
        let cols: Vec<Column> = (0..m).map(|i| vec![(i, 1.0)]).collect();
        let mut f = Factor::new(m, &cols).unwrap();
        for (r, alpha) in etas {
            f.update(*r, alpha.clone());
        }
        assert_eq!(
            f.nnz_eta,
            etas.iter()
                .map(|(_, a)| a.iter().filter(|a| **a != 0.0).count())
                .sum::<usize>()
        );
        let bits = |v: Vec<f64>| v.into_iter().map(f64::to_bits).collect::<Vec<_>>();
        let (mut old, mut new) = (rhs.to_vec(), rhs.to_vec());
        f.frozen_ftran(&mut old, etas);
        f.ftran(&mut new);
        assert_eq!(bits(old), bits(new), "forward");
        let (mut old, mut new) = (rhs.to_vec(), rhs.to_vec());
        f.frozen_btran(&mut old, etas);
        f.btran(&mut new);
        assert_eq!(bits(old), bits(new), "backward");
        let mut raw = Factor::new(m, &cols).unwrap();
        for (r, alpha) in etas {
            let mut column: Vec<_> = alpha
                .iter()
                .copied()
                .enumerate()
                .filter(|(i, a)| *i == *r || *a != 0.0)
                .collect();
            column.reverse();
            raw.update_column(*r, column);
        }
        assert_eq!(raw.nnz_eta, f.nnz_eta);
        let (mut old, mut new) = (rhs.to_vec(), rhs.to_vec());
        f.frozen_ftran(&mut old, etas);
        raw.ftran(&mut new);
        assert_eq!(bits(old), bits(new), "sparse forward");
        let (mut old, mut new) = (rhs.to_vec(), rhs.to_vec());
        f.frozen_btran(&mut old, etas);
        raw.btran(&mut new);
        assert_eq!(
            bits(old),
            bits(new),
            "sparse backward: etas={etas:?}, rhs={rhs:?}"
        );
    }

    fn check_columns<R: crate::real::Real>(
        m: usize,
        columns: &[(usize, Vec<(usize, R)>)],
        rhs: &[f64],
    ) {
        let basis: Vec<Column> = (0..m).map(|i| vec![(i, 1.0)]).collect();
        let mut f = Factor::new(m, &basis).unwrap();
        let mut etas = Vec::new();
        for (r, column) in columns {
            let mut alpha = vec![R::ZERO; m];
            for &(i, v) in column {
                alpha[i] = v;
            }
            etas.push((*r, alpha.into_iter().map(R::to_f64).collect::<Vec<_>>()));
            f.update_column(*r, column.clone());
            let diagonal = match f.etas.last().unwrap() {
                Eta::Dense(r, alpha) => alpha[*r],
                Eta::Sparse { diagonal, .. } => *diagonal,
            };
            assert_eq!(diagonal.to_bits(), etas.last().unwrap().1[*r].to_bits());
        }
        assert_eq!(
            f.nnz_eta,
            etas.iter()
                .map(|(_, a)| a.iter().filter(|a| **a != 0.0).count())
                .sum::<usize>()
        );
        let bits = |v: Vec<f64>| v.into_iter().map(f64::to_bits).collect::<Vec<_>>();
        let (mut old, mut new) = (rhs.to_vec(), rhs.to_vec());
        f.frozen_ftran(&mut old, &etas);
        f.ftran(&mut new);
        assert_eq!(bits(old), bits(new), "raw column forward");
        let (mut old, mut new) = (rhs.to_vec(), rhs.to_vec());
        f.frozen_btran(&mut old, &etas);
        f.btran(&mut new);
        assert_eq!(bits(old), bits(new), "raw column backward");
    }

    #[test]
    fn raw_columns_preserve_conversion_duplicates_and_diagonal_bits() {
        use crate::real::Real;

        fn cases<R: Real>(extra: &[R]) {
            let mut values: Vec<_> = [
                0.0,
                -0.0,
                1.0,
                -1.0,
                f64::from_bits(1),
                f64::MAX,
                f64::INFINITY,
                f64::from_bits(0x7ff8000000000123),
                f64::from_bits(0xfff8000000000456),
            ]
            .into_iter()
            .map(R::c)
            .collect();
            values.extend_from_slice(extra);
            for m in [1, 2, 7, 33] {
                for round in 0..values.len() {
                    let r = round % m;
                    let j = (r + 1) % m;
                    let mut columns = vec![
                        (r, Vec::new()),
                        (r, vec![(r, values[round])]),
                        (r, vec![(j, R::ONE), (r, values[round]), (j, R::ZERO)]),
                    ];
                    let dense: Vec<_> = (0..m)
                        .rev()
                        .map(|i| (i, values[(round + i) % values.len()]))
                        .collect();
                    columns.push((r, dense));
                    let rhs: Vec<_> = (0..m)
                        .map(|i| values[(i + round) % values.len()].to_f64())
                        .collect();
                    for n in 1..=columns.len() {
                        check_columns(m, &columns[..n], &rhs);
                    }
                }
            }
        }
        cases::<f64>(&[]);
        cases::<crate::f80::F80>(&[
            crate::f80::F80::from_repr(1, 0),
            crate::f80::F80::from_repr(1, 0x8000),
            crate::f80::F80::from_repr(u64::MAX, 0x7ffe),
        ]);
    }

    #[test]
    fn eta_solves_match_frozen_dense_kernels() {
        let values = [
            0.0,
            -0.0,
            1.0,
            -1.0,
            0.125,
            -0.25,
            f64::from_bits(1),
            -f64::MIN_POSITIVE,
            1e-200,
            -1e200,
            f64::MAX,
            f64::INFINITY,
            f64::NEG_INFINITY,
            f64::from_bits(0x7ff8000000000123),
            f64::from_bits(0xfff8000000000456),
        ];
        let mut seed = 0x97fa412357u64;
        let mut next = || {
            seed ^= seed << 13;
            seed ^= seed >> 7;
            seed ^= seed << 17;
            seed
        };
        for m in [0, 1, 2, 7, 33] {
            for round in 0..64 {
                let rhs: Vec<_> = (0..m)
                    .map(|_| values[next() as usize % values.len()])
                    .collect();
                let mut etas = Vec::new();
                if m > 0 {
                    for _ in 0..round % 5 {
                        let r = next() as usize % m;
                        let alpha: Vec<_> = (0..m)
                            .map(|_| {
                                if next() % 4 == 0 {
                                    values[next() as usize % values.len()]
                                } else {
                                    0.0
                                }
                            })
                            .collect();
                        etas.push((r, alpha));
                    }
                }
                check(m, &etas, &rhs);
            }
        }
        let captured: serde_json::Value =
            serde_json::from_str(include_str!("../tests/fixtures/aux_eta_ftran.json")).unwrap();
        let m = captured["m"].as_u64().unwrap() as usize;
        let mut etas = Vec::new();
        for col in captured["columns"].as_array().unwrap() {
            let mut alpha = vec![0.0; m];
            for entry in col.as_array().unwrap() {
                alpha[entry[0].as_u64().unwrap() as usize] = entry[1].as_f64().unwrap();
            }
            let r = col[0][0].as_u64().unwrap() as usize;
            etas.push((r, alpha.clone()));
            check(m, &etas, &alpha);
        }
    }
}
