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
    /// The range `|R|` of a ranged row as the file gives it (0 for none, or when absent); the solver's range column is `[0, row_range]`,
    /// which differs in the last bits from `row_up - row_lo`. Used only when it agrees with the bounds.
    pub row_range: Vec<f64>,
    /// The row type letters the file gives (`b'E'`, `b'L'`, `b'G'`); empty when unknown. A ranged G row is told from the other ranged rows only here.
    pub row_kind: Vec<u8>,
    /// The ranged rows in the order of the RANGES section (the order of the reference's range columns); empty when unknown (row order is used).
    pub range_order: Vec<usize>,
    /// Every RANGES entry on a constraint row, `(row, R)` in file order, zeros and signs as given (the barrier's range columns follow the entries).
    pub range_entries: Vec<(usize, f64)>,
    /// A row has two RANGES entries (the later one wins in the bounds; `crate::barrier::std_form` refuses such a file).
    pub range_dup: bool,
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

    /// Width of the range column of row `i` (`row_up - row_lo`, or the file's exact `|R|` when it agrees with the bounds).
    pub fn range_width(&self, i: usize) -> f64 {
        let w = self.row_up[i] - self.row_lo[i];
        match self.row_range.get(i) {
            Some(&r) if r > 0.0 && (w - r).abs() <= 1e-9 * r => r,
            _ => w,
        }
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
