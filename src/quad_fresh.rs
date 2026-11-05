//! Numerical reconstructions belonging to the existing bounded compact factor owners.
use super::{quad_main, QuadRhs, RhsMemo, Solver, NONE};
use crate::f80::F80;
use std::sync::{Arc, Mutex};

pub(super) struct DualSnapshot {
    basic: Vec<usize>,
    signs: Vec<i8>,
    costfree: bool,
    cleanup: bool,
    flags: u32,
    y: Vec<f64>,
    yq: Vec<F80>,
    rank: Vec<usize>,
}
impl DualSnapshot {
    pub(super) fn matches(
        &self,
        basic: &[usize],
        signs: &[i8],
        cf: bool,
        cleanup: bool,
        flags: u32,
    ) -> bool {
        self.costfree == cf
            && self.cleanup == cleanup
            && self.flags & !flags == 0
            && self.basic == basic
            && self.signs == signs
    }
    pub(super) fn capture(s: &Solver, signs: Vec<i8>, costfree: bool, flags: u32) -> Self {
        Self {
            basic: s.basic.clone(),
            signs,
            costfree,
            cleanup: s.cleanup,
            flags,
            y: s.y_ref.as_ref().unwrap().clone(),
            yq: s.y_q.as_ref().unwrap().clone(),
            rank: s.y_rank.clone(),
        }
    }
    pub(super) fn restore(&self, s: &mut Solver) {
        s.y_ref.get_or_insert_with(Vec::new).clone_from(&self.y);
        s.y_q.get_or_insert_with(Vec::new).clone_from(&self.yq);
        s.y_rank.clone_from(&self.rank);
    }
}

pub(super) fn primal_provider(s: &Solver) -> Option<Arc<Mutex<RhsMemo>>> {
    if s.p.m > 24
        || s.p.n + s.p.m > 64
        || !s
            .mirror
            .as_ref()
            .is_some_and(|m| m.is_q() && m.updates() == 0)
        || quad_main::status().is_none()
    {
        return None;
    }
    match s.mirror_q.as_ref()? {
        QuadRhs::Shared(p) => Some(p.clone()),
        _ => None,
    }
}

#[derive(Default)]
pub(super) struct PrimalMemo {
    pub(super) generation: u64,
    data: [Option<PrimalSnapshot>; 2],
}
struct PrimalSnapshot {
    basic: Vec<usize>,
    nonbasic: u64,
    upper: u64,
    shift: Vec<f64>,
    shiftq: Vec<F80>,
    cleanup: bool,
    x: Vec<f64>,
    xq: Vec<F80>,
    dst: Vec<f64>,
    dstq: Vec<F80>,
    drift: Option<F80>,
    work: f64,
    flags: u32,
}
fn masks(s: &Solver) -> (u64, u64) {
    let mut nonbasic = 0u64;
    let mut upper = 0u64;
    for (j, (&pos, &up)) in s.pos_of.iter().zip(&s.at_upper).enumerate() {
        if pos == NONE {
            nonbasic |= 1 << j;
            if up {
                upper |= 1 << j;
            }
        }
    }
    (nonbasic, upper)
}
impl PrimalMemo {
    // Both original RHS memo lookups must remain hits. A changed provider
    // invalidates all reconstruction snapshots before it mutates either cache.
    pub(super) fn invalidate(&mut self) {
        self.generation = self.generation.wrapping_add(1);
        self.data = Default::default();
    }
    pub(super) fn restore(&self, s: &mut Solver) -> bool {
        let Some(p) = &self.data[usize::from(s.want_drift)] else {
            return false;
        };
        let Some(flags) = quad_main::status() else {
            return false;
        };
        if p.flags & !flags != 0
            || p.cleanup != s.cleanup
            || p.basic != s.basic
            || p.shift.len() != s.rhs_shift.len()
            || p.shiftq.len() != s.rhs_shift_q.len()
            || !p
                .shift
                .iter()
                .zip(&s.rhs_shift)
                .all(|(a, b)| a.to_bits() == b.to_bits())
            || !p
                .shiftq
                .iter()
                .zip(&s.rhs_shift_q)
                .all(|(a, b)| a.same_bits(*b))
            || masks(s) != (p.nonbasic, p.upper)
        {
            return false;
        }
        s.work
            .add_from(0x1219d6, 0x120673, p.work, crate::work::FAC);
        s.x.clone_from(&p.x);
        s.x_q.get_or_insert_with(Vec::new).clone_from(&p.xq);
        let d = s.dst_q.get_or_insert_with(Vec::new);
        d.clear();
        d.extend(s.dst.iter().map(|&v| F80::from_f64(v)));
        for (k, &j) in s.basic.iter().enumerate() {
            s.dst[j] = p.dst[k];
            d[j] = p.dstq[k];
        }
        if s.want_drift {
            s.drift_err = p.drift;
        }
        true
    }
    pub(super) fn record(&mut self, s: &Solver, work: f64) {
        let slot = &mut self.data[usize::from(s.want_drift)];
        if slot.is_some() {
            return;
        }
        let Some(flags) = quad_main::status() else {
            return;
        };
        let (nonbasic, upper) = masks(s);
        *slot = Some(PrimalSnapshot {
            basic: s.basic.clone(),
            nonbasic,
            upper,
            cleanup: s.cleanup,
            shift: s.rhs_shift.clone(),
            shiftq: s.rhs_shift_q.clone(),
            x: s.x.clone(),
            xq: s.x_q.as_ref().unwrap().clone(),
            dst: s.basic.iter().map(|&j| s.dst[j]).collect(),
            dstq: s
                .basic
                .iter()
                .map(|&j| s.dst_q.as_ref().unwrap()[j])
                .collect(),
            drift: s.drift_err,
            work,
            flags,
        });
    }
}

#[cfg(test)]
pub(super) mod audit {
    use super::super::{Mir, SolveOptions, Status};
    use super::*;
    use std::cell::Cell;
    thread_local! {static ENABLED:Cell<bool>=const{Cell::new(false)};static EVENTS:Cell<[usize;6]>=const{Cell::new([0;6])};}
    fn count(k: usize) -> usize {
        EVENTS.with(|b| {
            let mut v = b.get();
            v[k] += 1;
            b.set(v);
            v[k]
        })
    }
    fn sample(k: usize) -> bool {
        if !ENABLED.get() {
            return false;
        }
        let n = count(k);
        n <= 48 || n % 89 == 0
    }
    fn bits(v: &[f64]) -> Vec<u64> {
        v.iter().map(|n| n.to_bits()).collect()
    }
    fn longs(v: &[F80]) -> Vec<(u64, u16, bool)> {
        v.iter()
            .map(|x| {
                let (m, e) = x.repr();
                (m, e, x.is_sign_negative())
            })
            .collect()
    }
    pub(in super::super) struct Before {
        x: Vec<f64>,
        q: Option<Vec<F80>>,
        d: Vec<f64>,
        dq: Option<Vec<F80>>,
        drift: Option<F80>,
        w: f64,
    }
    pub(in super::super) fn start_primal(s: &Solver) -> Option<Before> {
        sample(0).then(|| Before {
            x: s.x.clone(),
            q: s.x_q.clone(),
            d: s.dst.clone(),
            dq: s.dst_q.clone(),
            drift: s.drift_err,
            w: s.work.w,
        })
    }
    pub(in super::super) fn primal(s: &mut Solver, before: Option<Before>) {
        if !ENABLED.get() {
            return;
        }
        count(2);
        let Some(b) = before else {
            return;
        };
        count(4);
        let obs = |s: &Solver| {
            (
                bits(&s.x),
                longs(s.x_q.as_ref().unwrap()),
                bits(&s.dst),
                longs(s.dst_q.as_ref().unwrap()),
                s.drift_err.map(|f| longs(&[f])),
                s.work.w.to_bits(),
                quad_main::status(),
            )
        };
        let fast = obs(s);
        s.x = b.x;
        s.x_q = b.q;
        s.dst = b.d;
        s.dst_q = b.dq;
        s.drift_err = b.drift;
        s.work.w = b.w;
        s.recompute_primal_direct();
        assert_eq!(obs(s), fast, "all partial primal writes and ordered WORK");
    }
    pub(in super::super) struct DualBefore {
        factor: quad_main::QuadMain,
        w: f64,
        y: Option<Vec<f64>>,
        yq: Option<Vec<F80>>,
        rank: Vec<usize>,
    }
    pub(in super::super) fn start_dual(s: &Solver) -> Option<DualBefore> {
        if !sample(1) {
            return None;
        }
        let Some(Mir::Q(b)) = s.mirror.as_ref() else {
            return None;
        };
        Some(DualBefore {
            factor: quad_main::QuadMain::Direct((**b).clone()),
            w: s.work.w,
            y: s.y_ref.clone(),
            yq: s.y_q.clone(),
            rank: s.y_rank.clone(),
        })
    }
    pub(in super::super) fn dual(s: &mut Solver, site: &str, before: Option<DualBefore>) {
        if !ENABLED.get() {
            return;
        }
        count(3);
        let Some(b) = before else {
            return;
        };
        count(5);
        let obs = |s: &Solver| {
            let Some(Mir::Q(f)) = s.mirror.as_ref() else {
                panic!()
            };
            (
                format!("{:?}", (**f).clone()),
                bits(s.y_ref.as_ref().unwrap()),
                longs(s.y_q.as_ref().unwrap()),
                s.y_rank.clone(),
                s.y_costfree,
                s.work.w.to_bits(),
                quad_main::status(),
            )
        };
        let fast = obs(s);
        let after = s.mirror.replace(Mir::Q(b.factor));
        s.work.w = b.w;
        s.y_ref = b.y;
        s.y_q = b.yq;
        s.y_rank = b.rank;
        s.dual_btran_work(site);
        assert_eq!(
            obs(s),
            fast,
            "complete numeric BTRAN factor, outputs and ordered WORK"
        );
        s.mirror = after;
    }
    #[test]
    #[cfg(all(target_arch = "x86_64", target_os = "linux"))]
    fn fresh_reconstructions_match_the_live_unmodified_recurrences() {
        for (model, status, iters, refs) in [
            ("fx_quad_obj_5371", Status::Infeasible, 5319, 5),
            ("primal_5959", Status::NumericalFailure, 15688, 7842),
            ("quad_2717", Status::Optimal, 16010, 3022),
        ] {
            let model = crate::mps::parse(
                &std::fs::read_to_string(format!("tests/fixtures/{model}.mps")).unwrap(),
            )
            .unwrap();
            ENABLED.set(true);
            EVENTS.set([0; 6]);
            let result = crate::solve(&model, None, &SolveOptions::dual_release()).unwrap();
            ENABLED.set(false);
            assert_eq!(
                (result.status, result.iterations, result.refactorizations),
                (status, iters, refs)
            );
            let [_, _, ph, dh, checked_primal, checked_dual] = EVENTS.get();
            assert!(ph > 500 && dh > 500 && checked_primal >= 10 && checked_dual >= 10);
            eprintln!(
                "FRESH_ORACLE {} {} {} {}",
                ph, dh, checked_primal, checked_dual
            );
        }
    }
}
