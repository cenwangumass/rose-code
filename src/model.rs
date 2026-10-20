//! Continuous LP model stored by sparse columns.

pub const INF: f64 = 1e30;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Sense {
    Minimize,
    Maximize,
}

/// `min/max c'x + offset` subject to `row_lo <= A x <= row_up`, `col_lo <= x <= col_up`.
/// Bounds with magnitude `>= INF` are treated as infinite.
#[derive(Clone, Debug)]
pub struct LpModel {
    pub name: String,
    pub sense: Sense,
    pub offset: f64,
    pub obj: Vec<f64>,
    pub col_lo: Vec<f64>,
    pub col_up: Vec<f64>,
    pub row_lo: Vec<f64>,
    pub row_up: Vec<f64>,
    /// Column `j` holds `row_idx[col_start[j]..col_start[j + 1]]` with matching `val`.
    pub col_start: Vec<usize>,
    pub row_idx: Vec<usize>,
    pub val: Vec<f64>,
    pub col_names: Vec<String>,
    pub row_names: Vec<String>,
}

impl LpModel {
    pub fn num_cols(&self) -> usize {
        self.obj.len()
    }

    pub fn num_rows(&self) -> usize {
        self.row_lo.len()
    }

    pub fn col(&self, j: usize) -> impl Iterator<Item = (usize, f64)> + '_ {
        let (s, e) = (self.col_start[j], self.col_start[j + 1]);
        self.row_idx[s..e].iter().copied().zip(self.val[s..e].iter().copied())
    }

    /// Row activities `A x`.
    pub fn activities(&self, x: &[f64]) -> Vec<f64> {
        let mut r = vec![0.0; self.num_rows()];
        for j in 0..self.num_cols() {
            for (i, a) in self.col(j) {
                r[i] += a * x[j];
            }
        }
        r
    }

    pub fn objective(&self, x: &[f64]) -> f64 {
        self.obj.iter().zip(x).map(|(c, v)| c * v).sum::<f64>() + self.offset
    }
}
