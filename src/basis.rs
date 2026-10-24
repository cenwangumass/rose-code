//! Simplex basis over the internal variable space: columns `0..n`, then one slack per row `n..n+m`.

/// Basic variable per position plus the bound side of every nonbasic variable.
#[derive(Clone, Debug, PartialEq)]
pub struct Basis {
    /// `basic[p]` is the variable occupying basis position `p` (length `m`).
    pub basic: Vec<usize>,
    /// For nonbasic variables: `true` if at the upper bound. Length `n + m`; ignored for basic variables.
    pub at_upper: Vec<bool>,
}

impl Basis {
    /// All-slack basis with every column at its lower bound.
    pub fn slack(n: usize, m: usize) -> Basis {
        Basis {
            basic: (n..n + m).collect(),
            at_upper: vec![false; n + m],
        }
    }

    pub fn validate(&self, n: usize, m: usize) -> Result<(), String> {
        if self.basic.len() != m || self.at_upper.len() != n + m {
            return Err("basis dimensions do not match the model".into());
        }
        let mut seen = vec![false; n + m];
        for &v in &self.basic {
            if v >= n + m || seen[v] {
                return Err("basis has an invalid or repeated variable".into());
            }
            seen[v] = true;
        }
        Ok(())
    }
}
