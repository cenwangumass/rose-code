//! Representation-neutral, bit-exact shadow boundary. Pool offsets and scratch
//! allocation addresses are deliberately not part of this logical projection.
use crate::real::Real;
use std::any::Any;

pub(crate) fn bits<R: Real>(value: R) -> String {
    let value = &value as &dyn Any;
    if let Some(x) = value.downcast_ref::<f64>() {
        format!("0x{:016x}", x.to_bits())
    } else if let Some(x) = value.downcast_ref::<crate::f80::F80>() {
        let (negative, exponent, significand) = x.factor_bits();
        format!("0x{:01x}:{:08x}:{:016x}", negative as u8, exponent as u32, significand)
    } else {
        panic!("factor-shadow: unsupported Real type")
    }
}

pub(crate) fn numbers<R: Real>(values: &[R]) -> Vec<String> {
    values.iter().map(|&x| bits(x)).collect()
}

pub(crate) fn sparse<R: Real>(values: &[(usize, R)]) -> Vec<(usize, String)> {
    values.iter().map(|&(i, x)| (i, bits(x))).collect()
}

pub(crate) fn columns<R: Real>(values: &[Vec<(usize, R)>]) -> Vec<Vec<(usize, String)>> {
    values.iter().map(|v| sparse(v)).collect()
}

/// Returns the first difference, with a path that includes field and vector index.
pub(crate) trait Exact {
    fn difference(&self, old: &Self, path: &str) -> Option<String>;
}

macro_rules! atom {
    ($($t:ty),*) => { $(impl Exact for $t {
        fn difference(&self, old: &Self, path: &str) -> Option<String> {
            (self != old).then(|| format!("{path}: new={self:?}, old={old:?}"))
        }
    })* };
}
atom!(bool, String);

macro_rules! integer {
    ($($t:ty),*) => { $(impl Exact for $t {
        fn difference(&self, old: &Self, path: &str) -> Option<String> {
            (self != old).then(|| format!("{path}: new={self:#x}, old={old:#x}"))
        }
    })* };
}
integer!(usize, u32, u8, i64);

impl Exact for () {
    fn difference(&self, _: &Self, _: &str) -> Option<String> { None }
}

impl<T: Exact> Exact for Vec<T> {
    fn difference(&self, old: &Self, path: &str) -> Option<String> {
        if self.len() != old.len() {
            return Some(format!("{path}.len: new={}, old={}", self.len(), old.len()));
        }
        for (i, (a, b)) in self.iter().zip(old).enumerate() {
            if let Some(d) = a.difference(b, &format!("{path}[{i}]")) { return Some(d); }
        }
        None
    }
}

impl<T: Exact> Exact for Option<T> {
    fn difference(&self, old: &Self, path: &str) -> Option<String> {
        match (self, old) {
            (Some(a), Some(b)) => a.difference(b, path),
            (None, None) => None,
            _ => Some(format!("{path}.is_some: new={}, old={}", self.is_some(), old.is_some())),
        }
    }
}

impl<A: Exact, B: Exact> Exact for (A, B) {
    fn difference(&self, old: &Self, path: &str) -> Option<String> {
        self.0.difference(&old.0, &format!("{path}.0"))
            .or_else(|| self.1.difference(&old.1, &format!("{path}.1")))
    }
}

impl<A: Exact, B: Exact, C: Exact> Exact for (A, B, C) {
    fn difference(&self, old: &Self, path: &str) -> Option<String> {
        self.0.difference(&old.0, &format!("{path}.0"))
            .or_else(|| self.1.difference(&old.1, &format!("{path}.1")))
            .or_else(|| self.2.difference(&old.2, &format!("{path}.2")))
    }
}

pub(crate) fn check<T: Exact>(method: &str, call: usize, field: &str, new: &T, old: &T) {
    if let Some(d) = new.difference(old, field) {
        panic!("factor-shadow {method} call {call}: {d}");
    }
}

macro_rules! projection {
    ($name:ident { $($field:ident: $ty:ty),* $(,)? }) => {
        #[derive(Debug)]
        pub(crate) struct $name { $(pub(crate) $field: $ty),* }
        impl Exact for $name {
            fn difference(&self, old: &Self, path: &str) -> Option<String> {
                $(if let Some(d) = self.$field.difference(&old.$field, &format!("{path}.{}", stringify!($field))) {
                    return Some(d);
                })*
                None
            }
        }
    };
}

projection!(EtaState {
    target: usize,
    indices: Vec<usize>,
    values: Vec<String>,
    length: usize,
    alignment_peel: usize,
});

projection!(Canonical {
    m: usize,
    i21: usize,
    mark: usize,
    kind: Vec<i64>,
    nucleus: Vec<usize>,
    slot: Vec<usize>,
    at: Vec<usize>,
    diagonal: Vec<String>,
    lower_columns: Vec<Vec<(usize, String)>>,
    upper_columns: Vec<Vec<(usize, String)>>,
    upper_rows: Vec<Vec<(usize, String)>>,
    tail_destinations: Vec<Vec<usize>>,
    nucleus_destinations: Vec<Vec<usize>>,
    row_columns: Vec<Vec<usize>>,
    row_capacity: Vec<usize>,
    row_alignment: Vec<u8>,
    pool_fill: usize,
    etas: Vec<EtaState>,
    first_eta: Vec<usize>,
    updates: usize,
    alpha_sum: String,
    eta_sum: String,
    previous_window: (usize, String),
    last_alpha: usize,
    factor_size: String,
    tail_density: String,
    nucleus_density: String,
    lower_density_threshold: usize,
    cached_pivot: Option<usize>,
    spike: Vec<(usize, String)>,
    spike_consumed: bool,
    permutation: Vec<usize>,
    lower_transpose: Vec<Vec<(usize, String)>>,
    lower_levels: Option<(Vec<usize>, usize)>,
    empty_tail_positions: Vec<usize>,
    scratch_alignment: usize,
});

projection!(Observables {
    work: String,
    work_parts: (String, String),
    keep_work: String,
    keep_dense: bool,
    tau_lower_work: String,
    tau_upper_work: String,
    tau_count: usize,
    btran_dense: bool,
    last_dense: bool,
    last_ordered: bool,
    update_work: Vec<(u32, String, String)>,
    update_limit: Option<usize>,
    refactor_due: bool,
    upper_nnz: usize,
    has_levels: bool,
    updates: usize,
    alpha_sum: String,
    eta_sum: String,
    previous_window: (usize, String),
    last_alpha: usize,
    factor_size: String,
    density: String,
    dense_estimate: String,
    tail_density: String,
    nucleus_density: String,
    lower_density_threshold: usize,
    cached_pivot: Option<usize>,
    spike_consumed: bool,
    btran_occupied: bool,
    transpose_levels_dense: bool,
    forward_levels_dense: bool,
    phase_one_steps: usize,
    singleton_steps: usize,
    lower_nnz: usize,
    lower_tail_nnz: usize,
});

// Only the initial representation adapters know the old field names. A later
// backend can replace its adapter without changing either projection's schema.
macro_rules! impl_factor_projection {
    ($factor:ident) => {
        impl<R: crate::real::Real> $factor<R> {
            pub(crate) fn canonical(&self) -> crate::dual_factor::projection::Canonical {
                use crate::dual_factor::projection::{bits, columns, numbers, sparse, Canonical, EtaState};
                Canonical {
                    m: self.m, i21: self.i21, mark: self.mark,
                    kind: self.kind.clone(), nucleus: self.nuc.clone(),
                    slot: self.slot.clone(), at: self.at.clone(), diagonal: numbers(&self.diag),
                    lower_columns: self.lcolumns().map(sparse).collect(),
                    upper_columns: columns(&self.uc), upper_rows: columns(&self.ut),
                    tail_destinations: self.a.clone(), nucleus_destinations: self.b.clone(),
                    row_columns: self.rowl.clone(), row_capacity: self.row_cap.clone(), row_alignment: self.row_par.clone(),
                    pool_fill: self.fill,
                    etas: self.etas.iter().map(|e| EtaState {
                        target: e.idx, indices: e.ind.clone(), values: numbers(&e.val),
                        length: e.ind.len(), alignment_peel: e.peel,
                    }).collect(),
                    first_eta: self.first_eta.clone(), updates: self.upd,
                    alpha_sum: bits(self.sum_alpha), eta_sum: bits(self.sum_eta),
                    previous_window: (self.before_last.0, bits(self.before_last.1)),
                    last_alpha: self.last_alpha, factor_size: bits(self.znz),
                    tail_density: bits(self.h33), nucleus_density: bits(self.h34), lower_density_threshold: self.dthr,
                    cached_pivot: self.row_cached, spike: sparse(&self.spike), spike_consumed: self.spike_used,
                    permutation: self.perm.clone(), lower_transpose: columns(&self.lt), lower_levels: self.levels.clone(),
                    empty_tail_positions: self.empty_tail.iter().copied().collect(), scratch_alignment: self.scratch_peel,
                }
            }

            pub(crate) fn observables(&self) -> crate::dual_factor::projection::Observables {
                use crate::dual_factor::projection::{bits, Observables};
                Observables {
                    work: bits(self.work), work_parts: (bits(self.work_parts.0), bits(self.work_parts.1)),
                    keep_work: bits(self.keep_work), keep_dense: self.keep_dense,
                    tau_lower_work: bits(self.tau_l.get()), tau_upper_work: bits(self.tau_u.get()), tau_count: self.tau_fin.get(),
                    btran_dense: self.bt_dense, last_dense: self.last_dense, last_ordered: self.last_ordered,
                    update_work: self.upd_work.iter().map(|&(s, e, w)| (s, bits(e), bits(w))).collect(),
                    update_limit: self.upd_limit, refactor_due: self.refactor_due(), upper_nnz: self.u_nnz(),
                    has_levels: self.has_levels(), updates: self.upd,
                    alpha_sum: bits(self.sum_alpha), eta_sum: bits(self.sum_eta),
                    previous_window: (self.before_last.0, bits(self.before_last.1)),
                    last_alpha: self.last_alpha, factor_size: bits(self.znz), density: bits(self.dens),
                    dense_estimate: bits(self.dense_est.get()), tail_density: bits(self.h33), nucleus_density: bits(self.h34),
                    lower_density_threshold: self.dthr, cached_pivot: self.row_cached, spike_consumed: self.spike_used,
                    btran_occupied: self.bt_occupied, transpose_levels_dense: self.lev_dense,
                    forward_levels_dense: self.lev_dense_ft, phase_one_steps: self.n1, singleton_steps: self.loff,
                    lower_nnz: self.nnz_l, lower_tail_nnz: self.nnz_l_tail,
                }
            }
        }
    };
}
pub(crate) use impl_factor_projection;
