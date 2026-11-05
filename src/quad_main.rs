//! Complete compact main-factor states, rooted at one exact installation.
use crate::{f80::F80, ftran_order::RefStruct};
use std::{borrow::Cow, cell::OnceCell, rc::Rc};
type Factor = RefStruct<F80>;
type Values = Vec<(usize, F80)>;
type Work = Vec<(u32, f64, f64)>;

#[derive(Clone)]
pub(crate) enum QuadMain {
    Direct(Factor),
    Shared(Rc<State>),
}

pub(crate) struct State {
    factor: Factor,
    edge: OnceCell<Edge>,
    depth: usize,
    change: OnceCell<super::quad_change::Effect>,
    pricing_rows: super::quad_row::PricingRows,
    dual: OnceCell<super::quad_fresh::DualSnapshot>,
    #[cfg(all(target_arch = "x86_64", target_os = "linux"))]
    costs: super::quad_cost::CompactCosts,
}

struct Edge {
    input: Op<'static>,
    next: Rc<State>,
    out: Output,
    flags: u32,
}

#[derive(Clone)]
enum Output {
    Values(Values),
    Work(Work),
    Count,
}

enum Op<'a> {
    Column(Cow<'a, [usize]>, Cow<'a, [F80]>, bool),
    Transpose(Cow<'a, [(usize, F80)]>, bool, bool),
    Update(usize),
    Count,
}

impl Op<'_> {
    fn matches(&self, other: &Op<'_>) -> bool {
        match (self, other) {
            (Self::Column(r, v, b), Op::Column(s, w, c)) => {
                b == c
                    && r == s
                    && v.len() == w.len()
                    && v.iter().zip(w.iter()).all(|(a, b)| a.same_bits(*b))
            }
            (Self::Transpose(v, a, b), Op::Transpose(w, c, d)) => {
                a == c
                    && b == d
                    && v.len() == w.len()
                    && v.iter()
                        .zip(w.iter())
                        .all(|((i, a), (j, b))| i == j && a.same_bits(*b))
            }
            (Self::Update(a), Op::Update(b)) => a == b,
            (Self::Count, Op::Count) => true,
            _ => false,
        }
    }

    fn bounded(&self) -> bool {
        match self {
            Self::Column(r, v, _) => r.len() <= 100 && v.len() <= 100,
            Self::Transpose(v, _, _) => v.len() <= 100,
            _ => true,
        }
    }

    fn owned(self) -> Op<'static> {
        match self {
            Self::Column(r, v, b) => {
                Op::Column(Cow::Owned(r.into_owned()), Cow::Owned(v.into_owned()), b)
            }
            Self::Transpose(v, a, b) => Op::Transpose(Cow::Owned(v.into_owned()), a, b),
            Self::Update(r) => Op::Update(r),
            Self::Count => Op::Count,
        }
    }

    fn run(&self, m: &mut Factor) -> Output {
        match self {
            Self::Column(r, v, false) => Output::Values(m.ftran(r, v)),
            Self::Column(r, v, true) => Output::Values(m.ftran_basic(r, v)),
            Self::Transpose(v, a, b) => Output::Values(m.btran(v, *a, *b)),
            Self::Update(r) => {
                m.update(*r);
                Output::Work(std::mem::take(&mut m.upd_work))
            }
            Self::Count => {
                m.count_all_updates();
                Output::Count
            }
        }
    }
}

// Both exception units must retain all status flags observed after the original operation.
// Nondefault controls or an active factor dump always run the original numerical kernels.
pub(super) fn status() -> Option<u32> {
    #[cfg(all(target_arch = "x86_64", target_os = "linux"))]
    {
        let mut cw = 0u16;
        let mut sw = 0u16;
        let mut csr = 0u32;
        unsafe {
            std::arch::asm!("fnstcw [{cw}]", "fnstsw [{sw}]", "stmxcsr [{csr}]", cw=in(reg)&mut cw,sw=in(reg)&mut sw,csr=in(reg)&mut csr, options(nostack,preserves_flags));
        }
        if cw == 0x037f && csr & !63 == 0x1f80 && !crate::dump::any() {
            return Some((u32::from(sw & 63) << 6) | (csr & 63));
        }
    }
    None
}

impl std::ops::Deref for QuadMain {
    type Target = Factor;
    fn deref(&self) -> &Factor {
        match self {
            Self::Direct(f) => f,
            Self::Shared(s) => &s.factor,
        }
    }
}

impl QuadMain {
    pub(super) fn pricing_rows(&self) -> Option<&super::quad_row::PricingRows> {
        match self {
            Self::Shared(state) => Some(&state.pricing_rows),
            Self::Direct(_) => None,
        }
    }

    pub(super) fn new(f: Factor, eligible: bool) -> Self {
        if eligible && status().is_some() {
            Self::Shared(Rc::new(State {
                factor: f,
                edge: OnceCell::new(),
                depth: 0,
                change: OnceCell::new(),
                pricing_rows: Default::default(),
                dual: OnceCell::new(),
                #[cfg(all(target_arch = "x86_64", target_os = "linux"))]
                costs: Default::default(),
            }))
        } else {
            Self::Direct(f)
        }
    }

    #[cfg(all(target_arch = "x86_64", target_os = "linux"))]
    pub(super) fn costs(&self) -> Option<&super::quad_cost::CompactCosts> {
        match self {
            Self::Shared(state) => Some(&state.costs),
            Self::Direct(_) => None,
        }
    }

    pub(super) fn cursor(&self) -> Option<Rc<State>> {
        match self {
            Self::Shared(s) => Some(s.clone()),
            _ => None,
        }
    }

    pub(super) fn restart(&self) -> Self {
        match self {
            Self::Shared(s) if status().is_none() => Self::Direct(s.factor.clone()),
            _ => self.clone(),
        }
    }

    pub(super) fn reinstall(self, base: &Self) -> Self {
        match (self, base) {
            (Self::Direct(mut f), Self::Direct(b)) => {
                f.restore(b);
                Self::Direct(f)
            }
            _ => base.restart(),
        }
    }

    fn apply(&mut self, op: Op<'_>) -> Output {
        let state = match self {
            Self::Direct(f) => return op.run(f),
            Self::Shared(s) => s,
        };
        if let Some(flags) = status().filter(|_| op.bounded()) {
            if let Some(edge) = state.edge.get() {
                if edge.flags & !flags == 0 && edge.input.matches(&op) {
                    let out = edge.out.clone();
                    *state = edge.next.clone();
                    return out;
                }
            } else if state.depth < 48 {
                let mut next = state.factor.clone();
                let out = op.run(&mut next);
                let after = status().unwrap();
                let next = Rc::new(State {
                    factor: next,
                    edge: OnceCell::new(),
                    depth: state.depth + 1,
                    change: OnceCell::new(),
                    pricing_rows: Default::default(),
                    dual: OnceCell::new(),
                    #[cfg(all(target_arch = "x86_64", target_os = "linux"))]
                    costs: Default::default(),
                });
                state
                    .edge
                    .set(Edge {
                        input: op.owned(),
                        next: next.clone(),
                        out: out.clone(),
                        flags: after,
                    })
                    .ok()
                    .unwrap();
                *state = next;
                return out;
            }
        }
        let mut f = state.factor.clone();
        let out = op.run(&mut f);
        *self = Self::Direct(f);
        out
    }

    pub(super) fn btran(&mut self, v: &[(usize, F80)], direct: bool, dense: bool) -> Values {
        let Output::Values(out) = self.apply(Op::Transpose(Cow::Borrowed(v), direct, dense)) else {
            unreachable!()
        };
        out
    }
    pub(super) fn ftran(&mut self, r: &[usize], v: &[F80]) -> Values {
        self.column(r, v, false)
    }
    pub(super) fn ftran_basic(&mut self, r: &[usize], v: &[F80]) -> Values {
        self.column(r, v, true)
    }
    fn column(&mut self, r: &[usize], v: &[F80], basic: bool) -> Values {
        let Output::Values(out) = self.apply(Op::Column(Cow::Borrowed(r), Cow::Borrowed(v), basic))
        else {
            unreachable!()
        };
        out
    }
    pub(super) fn update(&mut self, r: usize) -> Work {
        let Output::Work(w) = self.apply(Op::Update(r)) else {
            unreachable!()
        };
        w
    }
    pub(super) fn count_all_updates(&mut self) {
        self.apply(Op::Count);
    }
}

#[cfg(all(test, target_arch = "x86_64", target_os = "linux"))]
mod tests {
    use super::*;
    fn factor(m: usize) -> Factor {
        let cols: Vec<_> = (0..m)
            .map(|i| (false, vec![i, (i + 1) % m], vec![2. + i as f64, 0.125]))
            .collect();
        let st = crate::lu_order::factor_structure_r::<F80>(m, &cols);
        RefStruct::fresh(m, m, &st, &cols)
    }
    fn fresh_controls() -> u32 {
        let mut csr = 0;
        unsafe {
            std::arch::asm!("stmxcsr [{p}]",p=in(reg)&mut csr,options(nostack,preserves_flags));
        }
        csr
    }
    fn set_controls(csr: u32) {
        unsafe {
            std::arch::asm!("ldmxcsr [{p}]",p=in(reg)&csr,options(nostack,preserves_flags));
        }
    }
    fn same(a: Output, b: Output) {
        match (a, b) {
            (Output::Values(a), Output::Values(b)) => {
                assert_eq!(a.len(), b.len());
                for ((i, x), (j, y)) in a.into_iter().zip(b) {
                    assert_eq!(i, j);
                    assert!(x.same_bits(y));
                }
            }
            (Output::Work(a), Output::Work(b)) => {
                assert_eq!(a.len(), b.len());
                for ((i, x, z), (j, y, w)) in a.into_iter().zip(b) {
                    assert_eq!((i, x.to_bits(), z.to_bits()), (j, y.to_bits(), w.to_bits()));
                }
            }
            (Output::Count, Output::Count) => (),
            _ => panic!("operation output changed"),
        }
    }
    fn compare(q: &mut QuadMain, f: &mut Factor, op: Op<'_>) {
        let expected = op.run(f);
        let actual = q.apply(op);
        same(expected, actual);
        // The original Clone intentionally clears independent BTRAN and level work areas.
        assert_eq!(format!("{:?}", (**q).clone()), format!("{:?}", f.clone()));
    }
    #[cfg(all(target_arch = "x86_64", target_os = "linux"))]
    #[test]
    fn compact_main_complete_transition_states() {
        let csr = fresh_controls();
        set_controls(0x1fbf);
        for m in [3, 7, 11] {
            let base = factor(m);
            let root = QuadMain::new(base.clone(), true);
            let QuadMain::Shared(top) = &root else {
                panic!("expected recorded factor")
            };
            let mut terminal = None;
            for _ in 0..3 {
                let mut f = base.clone();
                let mut q = root.restart();
                for n in 0..3 {
                    let r = (2 * n + 1) % m;
                    let rows = [r, (r + 2) % m, r];
                    let vals = [
                        F80::from_f64(-0.),
                        F80::from_f64(-0.125),
                        F80::from_f64(0.5),
                    ];
                    let rhs = [(r, F80::from_f64(-0.75))];
                    compare(
                        &mut q,
                        &mut f,
                        Op::Column((&rows[..]).into(), (&vals[..]).into(), false),
                    );
                    compare(
                        &mut q,
                        &mut f,
                        Op::Transpose((&rhs[..]).into(), true, false),
                    );
                    compare(&mut q, &mut f, Op::Update(r));
                    compare(&mut q, &mut f, Op::Count);
                    compare(
                        &mut q,
                        &mut f,
                        Op::Transpose((&rhs[..]).into(), false, true),
                    );
                    compare(
                        &mut q,
                        &mut f,
                        Op::Column((&rows[..]).into(), (&vals[..]).into(), true),
                    );
                }
                let QuadMain::Shared(end) = q else {
                    panic!("unexpected fallback")
                };
                assert_eq!(end.depth, 18);
                if let Some(want) = &terminal {
                    assert!(Rc::ptr_eq(want, &end));
                }
                terminal = Some(end);
            }
            assert!(top.edge.get().is_some());
            let mut branch = root.restart();
            let mut f = base.clone();
            let rows = [0; 101];
            let vals = [F80::ONE; 101];
            compare(
                &mut branch,
                &mut f,
                Op::Column((&rows[..]).into(), (&vals[..]).into(), false),
            );
            assert!(matches!(branch, QuadMain::Direct(_)));
            let mut change = root.restart();
            let mut f = base.clone();
            compare(
                &mut change,
                &mut f,
                Op::Transpose(Cow::Borrowed(&[]), false, false),
            );
            assert!(matches!(change, QuadMain::Direct(_)));
            let independent = QuadMain::new(base.clone(), true);
            let QuadMain::Shared(other) = independent else {
                panic!()
            };
            assert!(!Rc::ptr_eq(top, &other));
            assert!(other.edge.get().is_none());
            let p = F80::from_f64(0.);
            let z = F80::from_f64(-0.);
            let vals1 = [(0, p)];
            let vals2 = [(0, z)];
            assert!(
                !Op::Transpose((&vals1[..]).into(), true, true).matches(&Op::Transpose(
                    (&vals2[..]).into(),
                    true,
                    true
                ))
            );
            assert!(!F80::NAN.same_bits(-F80::NAN));
        }
        let root = QuadMain::new(factor(3), true);
        let mut q = root.restart();
        for _ in 0..48 {
            q.count_all_updates();
        }
        assert!(matches!(q, QuadMain::Shared(_)));
        q.count_all_updates();
        assert!(matches!(q, QuadMain::Direct(_)));
        let mut q = root.restart();
        set_controls(0x1f80);
        q.count_all_updates();
        assert!(matches!(q, QuadMain::Direct(_)));
        set_controls(0x3f80);
        assert!(matches!(root.restart(), QuadMain::Direct(_)));
        set_controls(csr);
    }
}

impl State {
    pub(super) fn try_dual<'a>(
        &'a self,
        live: &mut QuadMain,
        basic: &[usize],
        classes: &[i8],
        cf: bool,
        cleanup: bool,
        flags: u32,
    ) -> Option<&'a super::quad_fresh::DualSnapshot> {
        let data = self.dual.get()?;
        if !data.matches(basic, classes, cf, cleanup, flags) {
            return None;
        }
        let edge = self.edge.get()?;
        if edge.flags & !flags != 0 {
            return None;
        }
        *live = QuadMain::Shared(edge.next.clone());
        Some(data)
    }
    pub(super) fn may_record(&self, post: &QuadMain) -> bool {
        self.dual.get().is_none()
            && self.edge.get().is_some_and(|edge| {
                matches!(edge.input, Op::Transpose(_, false, true))
                    && matches!(post,QuadMain::Shared(after) if Rc::ptr_eq(after,&edge.next))
            })
    }
    pub(super) fn record_dual(&self, d: super::quad_fresh::DualSnapshot) {
        let _ = self.dual.set(d);
    }
}

impl State {
    pub(super) fn change(&self) -> Option<&super::quad_change::Effect> {
        self.change.get()
    }
    pub(super) fn advance_change(&self, live: &mut QuadMain, flags: u32) -> bool {
        let Some(edge) = self.edge.get() else {
            return false;
        };
        if edge.flags & !flags != 0 || !matches!(edge.input, Op::Transpose(_, false, false)) {
            return false;
        }
        *live = QuadMain::Shared(edge.next.clone());
        true
    }
    pub(super) fn record_change(&self, live: &QuadMain, effect: super::quad_change::Effect) {
        if let Some(edge) = self.edge.get() {
            if matches!(edge.input, Op::Transpose(_, false, false))
                && matches!(live,QuadMain::Shared(s) if Rc::ptr_eq(s,&edge.next))
            {
                let _ = self.change.set(effect);
            }
        }
    }
}
