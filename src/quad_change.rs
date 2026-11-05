//! The general phase-one cost update of a compact, replayable main factor.
use super::{quad_main, Mir, QuadRhs, RhsMemo, Solver, NONE};
use crate::f80::F80;
use std::{
    rc::Rc,
    sync::{Arc, Mutex, Weak},
};

pub(super) fn eligible(
    s: &Solver,
    dc: &[(usize, f64)],
    q: usize,
    q_only: bool,
) -> Option<(Rc<quad_main::State>, u32)> {
    if q_only || s.p.m > 24 || s.p.n + s.p.m > 64 || dc.len() > 64 || s.list.len() > 128 {
        return None;
    }
    // The separate entering-only update reads a prior row and may recurse.
    if q != NONE && s.pos_of[q] != NONE && dc.len() > 1 && dc.iter().any(|&(r, _)| r == s.pos_of[q])
    {
        return None;
    }
    let Some(Mir::Q(main)) = s.mirror.as_ref() else {
        return None;
    };
    let Some(QuadRhs::Shared(_)) = s.mirror_q.as_ref() else {
        return None;
    };
    let state = main.cursor()?;
    if s.rrows.iter().map(Vec::len).sum::<usize>() > 1024 {
        return None;
    }
    Some((state, quad_main::status()?))
}

pub(super) struct Key {
    dc: Vec<(usize, u64)>,
    q: usize,
    leave: usize,
    alpha: u64,
    basic: Vec<usize>,
    pos: Vec<usize>,
    upper: u64,
    excluded: Option<usize>,
    excl6: u64,
    d: Vec<f64>,
    df: Vec<F80>,
    w: Vec<f64>,
    rows: Vec<Vec<usize>>,
    list: Vec<usize>,
    words: Vec<u64>,
    level: usize,
    p1: bool,
    p3: u64,
    p0: bool,
    valid: bool,
    provider: Weak<Mutex<RhsMemo>>,
}

fn bits(v: &[bool]) -> u64 {
    v.iter()
        .enumerate()
        .fold(0, |m, (i, &b)| m | (u64::from(b) << i))
}

fn doubles(a: &[f64], b: &[f64]) -> bool {
    a.len() == b.len() && a.iter().zip(b).all(|(x, y)| x.to_bits() == y.to_bits())
}

impl Key {
    pub(super) fn capture(
        s: &Solver,
        dc: &[(usize, f64)],
        q: usize,
        leave: usize,
        alpha: f64,
    ) -> Self {
        let Some(QuadRhs::Shared(provider)) = s.mirror_q.as_ref() else {
            unreachable!()
        };
        Self {
            dc: dc.iter().map(|&(r, v)| (r, v.to_bits())).collect(),
            q,
            leave,
            alpha: alpha.to_bits(),
            basic: s.basic.clone(),
            pos: s.pos_of.clone(),
            upper: bits(&s.at_upper),
            excluded: s.excluded,
            excl6: bits(&s.excl6),
            d: s.d.clone(),
            df: s.df.clone(),
            w: s.w.clone(),
            rows: s.rrows.clone(),
            list: s.list.clone(),
            words: s.pivot_order.words.clone(),
            level: s.level,
            p1: s.p1_tol_on,
            p3: s.p3c0.to_bits(),
            p0: s.p0.is_some(),
            valid: s.row_work_valid,
            provider: Arc::downgrade(provider),
        }
    }

    fn matches(&self, s: &Solver, dc: &[(usize, f64)], q: usize, leave: usize, alpha: f64) -> bool {
        let Some(QuadRhs::Shared(provider)) = s.mirror_q.as_ref() else {
            return false;
        };
        self.q == q
            && self.leave == leave
            && self.alpha == alpha.to_bits()
            && self.level == s.level
            && self.excluded == s.excluded
            && self.p1 == s.p1_tol_on
            && self.p3 == s.p3c0.to_bits()
            && self.p0 == s.p0.is_some()
            && self.valid == s.row_work_valid
            && self.provider.as_ptr() == Arc::as_ptr(provider)
            && self.dc.len() == dc.len()
            && self
                .dc
                .iter()
                .zip(dc)
                .all(|(&(i, a), &(j, b))| i == j && a == b.to_bits())
            && self.basic == s.basic
            && self.pos == s.pos_of
            && self.list == s.list
            && self.words == s.pivot_order.words
            && self.upper == bits(&s.at_upper)
            && self.excl6 == bits(&s.excl6)
            && doubles(&self.d, &s.d)
            && doubles(&self.w, &s.w)
            && self.df.len() == s.df.len()
            && self.df.iter().zip(&s.df).all(|(a, b)| a.same_bits(*b))
            && self.rows == s.rrows
    }
}

pub(super) struct Effect {
    key: Key,
    rank: Weak<()>,
    d: Vec<f64>,
    df: Vec<F80>,
    list: Vec<usize>,
    words: Vec<u64>,
    level: usize,
    work: [f64; 3],
    flags: u32,
}

impl Effect {
    pub(super) fn capture(
        s: &Solver,
        key: Key,
        rank: Weak<()>,
        work: [f64; 3],
        flags: u32,
    ) -> Self {
        Self {
            key,
            rank,
            d: s.d.clone(),
            df: s.df.clone(),
            list: s.list.clone(),
            words: s.pivot_order.words.clone(),
            level: s.level,
            work,
            flags,
        }
    }

    pub(super) fn replay(
        &self,
        s: &mut Solver,
        state: &quad_main::State,
        dc: &[(usize, f64)],
        q: usize,
        leave: usize,
        alpha: f64,
        flags: u32,
    ) -> bool {
        // The rank memo's key/value is immutable. It alone holds a strong
        // reference to the witness. Eviction forbids replay of the skipped hit.
        if self.flags & !flags != 0
            || self.rank.strong_count() == 0
            || !self.key.matches(s, dc, q, leave, alpha)
        {
            return false;
        }
        let Some(Mir::Q(main)) = &mut s.mirror else {
            return false;
        };
        if !state.advance_change(main, flags) {
            return false;
        }
        s.work
            .add_from(0x12de50, 0x10804b, self.work[0], crate::work::FAC);
        s.work
            .add_from(0x135e92, 0x10820b, self.work[1], crate::work::CTX);
        s.d.clone_from(&self.d);
        s.df.clone_from(&self.df);
        s.pivot_order.words.clone_from(&self.words);
        s.level = self.level;
        s.list.clone_from(&self.list);
        s.work
            .add_from(0x1087b1, 0x108280, self.work[2], crate::work::CTX);
        true
    }
}

#[cfg(test)]
pub(super) mod audit {
    use super::super::{SolveOptions, Status};
    use super::*;
    use std::cell::Cell;

    thread_local! {static RUN:Cell<bool>=const{Cell::new(false)};static COUNTS:Cell<[usize;3]>=const{Cell::new([0;3])};}

    fn note(k: usize) -> usize {
        COUNTS.with(|b| {
            let mut a = b.get();
            a[k] += 1;
            b.set(a);
            a[k]
        })
    }

    pub(in super::super) struct Before {
        d: Vec<f64>,
        df: Vec<F80>,
        list: Vec<usize>,
        level: usize,
        words: Vec<u64>,
        w: f64,
        factor: quad_main::QuadMain,
    }

    pub(in super::super) fn start(s: &Solver) -> Option<Before> {
        if !RUN.get() {
            return None;
        }
        let n = note(0);
        if n > 32 && n % 71 != 0 {
            return None;
        }
        let Some(Mir::Q(f)) = &s.mirror else {
            return None;
        };
        Some(Before {
            d: s.d.clone(),
            df: s.df.clone(),
            list: s.list.clone(),
            level: s.level,
            words: s.pivot_order.words.clone(),
            w: s.work.w,
            factor: quad_main::QuadMain::Direct(Box::new((**f).clone())),
        })
    }

    pub(in super::super) fn check(
        s: &mut Solver,
        b: Option<Before>,
        dc: &[(usize, f64)],
        q: usize,
        leave: usize,
        alpha: f64,
        effect: &Effect,
    ) {
        if !RUN.get() {
            return;
        }
        note(1);
        let Some(b) = b else {
            return;
        };
        note(2);
        let obs = |s: &Solver| {
            let Some(Mir::Q(f)) = &s.mirror else { panic!() };
            (
                s.d.iter().map(|v| v.to_bits()).collect::<Vec<_>>(),
                s.df.iter()
                    .map(|v| {
                        let (m, e) = v.repr();
                        (m, e, v.is_sign_negative())
                    })
                    .collect::<Vec<_>>(),
                s.level,
                s.list.clone(),
                s.pivot_order.words.clone(),
                format!("{:?}", (**f).clone()),
                s.last_upd.clone(),
                s.last_row
                    .iter()
                    .map(|(i, n)| (*i, n.to_bits()))
                    .collect::<Vec<_>>(),
                s.work.w.to_bits(),
                quad_main::status(),
            )
        };
        let fast = obs(s);
        let owner = s.mirror.replace(Mir::Q(b.factor));
        s.d = b.d;
        s.df = b.df;
        s.list = b.list;
        s.level = b.level;
        s.pivot_order.words = b.words;
        s.work.w = b.w;
        // Signed-zero or a one-bit input change must not match a different operation.
        let original = s.d[0];
        s.d[0] = f64::from_bits(original.to_bits() ^ 1);
        assert!(!effect.key.matches(s, dc, q, leave, alpha));
        s.d[0] = original;
        assert!(effect.key.matches(s, dc, q, leave, alpha));
        assert!(!effect
            .key
            .matches(s, dc, q, leave, f64::from_bits(alpha.to_bits() ^ 1)));
        let mut witness = None;
        let work = s
            .cost_change_original(dc, q, leave, alpha, false, Some(&mut witness))
            .unwrap();
        assert_eq!(work.map(f64::to_bits), effect.work.map(f64::to_bits));
        assert!(Weak::ptr_eq(&witness.unwrap(), &effect.rank));
        assert_eq!(
            obs(s),
            fast,
            "all partial writes, ordered WORK and complete factor state"
        );
        s.mirror = owner;
    }

    #[test]
    #[cfg(all(target_arch = "x86_64", target_os = "linux"))]
    fn complete_phase_cost_updates_match_unmodified_controller() {
        for (name, status, iters, refs) in [
            ("quad_2717", Status::Optimal, 16010, 3022),
            ("fx_quad_obj_5371", Status::Infeasible, 5319, 5),
        ] {
            let model = crate::mps::parse(
                &std::fs::read_to_string(format!("tests/fixtures/{name}.mps")).unwrap(),
            )
            .unwrap();
            RUN.set(true);
            COUNTS.set([0; 3]);
            let sol = crate::solve(&model, None, &SolveOptions::dual_release()).unwrap();
            RUN.set(false);
            assert_eq!(
                (sol.status, sol.iterations, sol.refactorizations),
                (status, iters, refs)
            );
            let [_, hits, compared] = COUNTS.get();
            eprintln!("CHANGE_ORACLE {name} {hits} {compared}");
            assert!(hits > 500 && compared > 10);
        }
    }
}
