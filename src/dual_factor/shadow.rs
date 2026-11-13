//! Lockstep oracle. Only the new backend can publish diagnostics or advance
//! diagnostic counters. Comparisons never round an extended value to binary64.
use super::projection::{bits, check, columns, numbers, sparse, Exact};
use super::DualFactor;
use crate::dump::muted;
use crate::ftran_order::RefStruct;
use crate::lu_order::Structure;
use crate::real::Real;
use std::cell::Cell;

pub struct Shadow<R: Real = f64> {
    new: DualFactor<R>,
    old: RefStruct<R>,
    calls: Cell<usize>,
    /// The old backend's unfiltered result of the last `btran_rho`.
    old_full: Vec<(usize, R)>,
}

// The encoder is supplied explicitly, so generic Real values cannot accidentally
// use numeric equality (which equates signed zeros and rejects identical NaNs).
macro_rules! forward {
    ($name:ident ($($arg:ident: $ty:ty),*) -> $ret:ty, $encode:expr) => {
        pub fn $name(&self, $($arg: $ty),*) -> $ret {
            let new = self.new.$name($($arg),*);
            let old = muted(|| self.old.$name($($arg),*));
            let encode = $encode;
            self.finish(stringify!($name), &encode(&new), &encode(&old), false);
            new
        }
    };
    (mut $name:ident ($($arg:ident: $ty:ty),*) -> $ret:ty, $encode:expr) => {
        pub fn $name(&mut self, $($arg: $ty),*) -> $ret {
            let new = self.new.$name($($arg),*);
            let old = muted(|| self.old.$name($($arg),*));
            let encode = $encode;
            self.finish(stringify!($name), &encode(&new), &encode(&old), false);
            new
        }
    };
}

impl<R: Real> Shadow<R> {
    /// See `DualFactor::give_list`.
    pub fn give_list(&self, list: Vec<(usize, R)>) {
        self.new.give_list(list);
    }

    /// See `DualFactor::recycle`; the oracle is dropped.
    pub fn recycle(self) {
        self.new.recycle();
    }

    pub fn fresh<'a>(m: usize, ncols: usize, st: &Structure<R>, cols: impl Into<crate::lu_order::Cols<'a>>) -> Self {
        let cols: crate::lu_order::Cols<'a> = cols.into();
        let new = DualFactor::fresh(m, ncols, st, cols);
        let old = muted(|| RefStruct::fresh(m, ncols, st, cols));
        let result = Self { new, old, calls: Cell::new(0), old_full: Vec::new() };
        result.finish("fresh", &(), &(), true);
        result
    }

    fn finish<T: Exact>(&self, method: &str, new: &T, old: &T, canonical: bool) {
        let call = self.calls.get() + 1;
        self.calls.set(call);
        self.new.check_btran_workspace();
        check(method, call, "return", new, old);
        check(method, call, "observable", &self.new.observables(), &muted(|| self.old.observables()));
        // Every solver transition, not just updates, can leave state consumed
        // by a later call. Accessors only change already-compared observables.
        if canonical || matches!(method, "ftran" | "ftran_needed" | "ftran_basic" | "ftran_keep" | "btran") {
            self.new.check_retained(method, call);
            muted(|| self.old.check_retained(method, call));
            let a = self.new.f2_projection();
            let b = self.old.f2_projection();
            check(method, call, "f2.col_n", &a.column_count, &b.column_count);
            if a.column_input != b.column_input { check(method, call, "f2.col_in", &a.column_input, &b.column_input); }
            super::projection::check_numbers(method, call, "f2.x_col", &a.column, &b.column);
            check(method, call, "f2.tau_n", &a.tau_count, &b.tau_count);
            super::projection::check_numbers(method, call, "f2.x_tau", &a.tau, &b.tau);
            if a.tau_raw != b.tau_raw { check(method, call, "f2.tau_raw", &a.tau_raw, &b.tau_raw); }
        }
        if canonical {
            check(method, call, "canonical", &self.new.canonical(), &muted(|| self.old.canonical()));
        }
    }

    forward!(mut ftran(rows: &[usize], vals: &[R]) -> Vec<(usize, R)>, |v: &Vec<(usize, R)>| sparse(v));
    forward!(ftran_needed(rows: &[usize], vals: &[R], needed: &[usize], dense_in: bool) -> Vec<(usize, R)>,
        |v: &Vec<(usize, R)>| sparse(v));
    forward!(mut ftran_basic(rows: &[usize], vals: &[R]) -> Vec<(usize, R)>, |v: &Vec<(usize, R)>| sparse(v));
    forward!(mut ftran_keep(rows: &[usize], vals: &[R], dense: bool) -> Vec<(usize, R)>, |v: &Vec<(usize, R)>| sparse(v));
    forward!(mut btran(inp: &[(usize, R)], direct: bool, dense_in: bool) -> Vec<(usize, R)>,
        |v: &Vec<(usize, R)>| sparse(v));
    /// See `DualFactor::btran_rho`; the old backend's `btran` result, filtered, must equal it.
    pub fn btran_rho(&mut self, r: usize, tol: R, rows: &mut Vec<usize>, vals: &mut Vec<R>) {
        self.new.btran_rho(r, tol, rows, vals);
        let old = muted(|| self.old.btran(&[(r, R::ONE)], true, false));
        let kept: Vec<(usize, R)> = old.iter().copied().filter(|e| e.1.abs() > tol).collect();
        let new: Vec<(usize, R)> = rows.iter().copied().zip(vals.iter().copied()).collect();
        self.finish("btran", &sparse(&new), &sparse(&kept), false);
        self.old_full = old;
    }

    /// See `DualFactor::rho_full`; compared with the old backend's unfiltered list.
    pub fn rho_full(&mut self, rows: &[usize], vals: &[R]) -> Vec<(usize, R)> {
        let new = self.new.rho_full(rows, vals);
        let old = std::mem::take(&mut self.old_full);
        self.finish("rho_full", &sparse(&new), &sparse(&old), false);
        new
    }

    forward!(joint_u(out: &[usize]) -> Option<f64>, |v: &Option<f64>| v.map(bits));
    forward!(ftran_dense_work(b: &[R]) -> f64, |v: &f64| bits(*v));
    forward!(discard_work() -> f64, |v: &f64| bits(*v));
    forward!(fresh_work() -> Vec<(u32, f64)>, |v: &Vec<(u32, f64)>| {
        v.iter().map(|&(site, work)| (site, bits(work))).collect::<Vec<_>>()
    });
    forward!(refactor_due() -> bool, |v: &bool| *v);
    forward!(has_levels() -> bool, |v: &bool| *v);
    forward!(u_nnz() -> usize, |v: &usize| *v);
    forward!(last_dense() -> bool, |v: &bool| *v);
    forward!(last_ordered() -> bool, |v: &bool| *v);

    pub fn reset_forward(&self, v: usize) -> &[u32] {
        self.new.reset_forward(v)
    }

    pub fn reset_backward(&self, y: usize) -> &[u32] {
        self.new.reset_backward(y)
    }

    pub fn diag(&self) -> &[R] {
        let new = self.new.diag();
        let old = muted(|| self.old.numerics().0);
        self.finish("diag", &numbers(new), &numbers(old), false);
        new
    }

    pub fn numerics(&self) -> (&[R], Vec<Vec<(usize, R)>>) {
        let new = self.new.numerics();
        let old = muted(|| self.old.numerics());
        self.finish("numerics", &(numbers(new.0), columns(&new.1)), &(numbers(old.0), columns(&old.1)), false);
        new
    }

    pub fn u_numerics(&self) -> Vec<Vec<(usize, R)>> {
        let new = self.new.u_numerics();
        let old = muted(|| self.old.u_numerics());
        self.finish("u_numerics", &columns(&new), &columns(old), false);
        new
    }

    pub fn update(&mut self, row: usize) {
        self.new.update(row);
        muted(|| self.old.update(row));
        self.finish("update", &(), &(), true);
    }

    pub fn work(&self) -> f64 {
        let new = self.new.work();
        let old = muted(|| self.old.work);
        self.finish("work", &bits(new), &bits(old), false);
        new
    }
    pub fn work_parts(&self) -> (f64, f64) {
        let new = self.new.work_parts();
        let old = muted(|| self.old.work_parts);
        self.finish("work_parts", &(bits(new.0), bits(new.1)), &(bits(old.0), bits(old.1)), false);
        new
    }
    /// Pure functions of the new backend's route and work; the old backend has no counterpart.
    pub fn dense_route(&self, nrows: usize) -> bool {
        self.new.dense_route(nrows)
    }
    pub fn dense_work_of_solve(&self, solved: f64, nrows: usize) -> f64 {
        self.new.dense_work_of_solve(solved, nrows)
    }
    pub fn keep_work(&self) -> f64 {
        let new = self.new.keep_work();
        let old = muted(|| self.old.keep_work);
        self.finish("keep_work", &bits(new), &bits(old), false);
        new
    }
    pub fn keep_dense(&self) -> bool {
        let new = self.new.keep_dense();
        let old = muted(|| self.old.keep_dense);
        self.finish("keep_dense", &new, &old, false);
        new
    }
    pub fn tau_l(&self) -> f64 {
        let new = self.new.tau_l();
        let old = muted(|| self.old.tau_l.get());
        self.finish("tau_l", &bits(new), &bits(old), false);
        new
    }
    pub fn tau_u(&self) -> f64 {
        let new = self.new.tau_u();
        let old = muted(|| self.old.tau_u.get());
        self.finish("tau_u", &bits(new), &bits(old), false);
        new
    }
    pub fn tau_fin(&self) -> usize {
        let new = self.new.tau_fin();
        let old = muted(|| self.old.tau_fin.get());
        self.finish("tau_fin", &new, &old, false);
        new
    }
    pub fn bt_dense(&self) -> bool {
        let new = self.new.bt_dense();
        let old = muted(|| self.old.bt_dense);
        self.finish("bt_dense", &new, &old, false);
        new
    }
    pub fn upd_limit(&self) -> Option<usize> {
        let new = self.new.upd_limit();
        let old = muted(|| self.old.upd_limit);
        self.finish("upd_limit", &new, &old, false);
        new
    }
    pub fn set_upd_limit(&mut self, limit: Option<usize>) {
        self.new.set_upd_limit(limit);
        muted(|| self.old.upd_limit = limit);
        self.finish("set_upd_limit", &(), &(), false);
    }
    pub fn clear_update_work(&mut self) {
        self.new.clear_update_work();
        muted(|| self.old.upd_work.clear());
        self.finish("clear_update_work", &(), &(), false);
    }
    /// See `DualFactor::give_update_work`.
    pub fn give_update_work(&mut self, list: Vec<(u32, f64, f64)>) {
        self.new.give_update_work(list);
    }
    pub fn take_update_work(&mut self) -> Vec<(u32, f64, f64)> {
        let new = self.new.take_update_work();
        let old = muted(|| std::mem::take(&mut self.old.upd_work));
        let encode = |v: &Vec<(u32, f64, f64)>| {
            v.iter().map(|&(site, work, scale)| (site, bits(work), bits(scale))).collect::<Vec<_>>()
        };
        self.finish("take_update_work", &encode(&new), &encode(&old), false);
        new
    }
}

impl<R: Real> Clone for Shadow<R> {
    fn clone(&self) -> Self {
        let new = self.new.clone();
        let old = muted(|| self.old.clone());
        self.finish("clone", &(), &(), true);
        let result = Self { new, old, calls: Cell::new(self.calls.get()), old_full: Vec::new() };
        result.finish("clone result", &(), &(), true);
        result
    }
}
