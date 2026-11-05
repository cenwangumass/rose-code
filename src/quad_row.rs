//! Sparse pricing-row transitions at a retained compact main-factor state.
use super::{quad_main, PrimalRowInput, PrimalRowState};
use crate::f80::F80;
use std::{cell::RefCell, sync::Arc};

pub(super) type ResultRow = (bool, Vec<usize>, Vec<(usize, f64)>);

#[derive(Default)]
pub(super) struct PricingRows(RefCell<Option<Plan>>);

fn doubles(a: &[f64], b: &[f64]) -> bool {
    a.len() == b.len() && a.iter().zip(b).all(|(a, b)| a.to_bits() == b.to_bits())
}

fn extended(a: &[F80], b: &[F80]) -> bool {
    a.len() == b.len() && a.iter().zip(b).all(|(a, b)| a.same_bits(*b))
}

// Only returned row positions can read or modify the four numerical price arrays.
#[derive(Clone, Copy)]
struct Price {
    weight: f64,
    cost: f64,
    extended_weight: F80,
    extended_cost: F80,
}

impl Price {
    fn at(s: &PrimalRowState<'_>, j: usize) -> Self {
        Self {
            weight: s.w[j],
            cost: s.d[j],
            extended_weight: s.wf[j],
            extended_cost: s.df[j],
        }
    }

    fn matches(&self, s: &PrimalRowState<'_>, j: usize) -> bool {
        self.weight.to_bits() == s.w[j].to_bits()
            && self.cost.to_bits() == s.d[j].to_bits()
            && self.extended_weight.same_bits(s.wf[j])
            && self.extended_cost.same_bits(s.df[j])
    }

    fn restore(self, s: &mut PrimalRowState<'_>, j: usize) {
        s.w[j] = self.weight;
        s.d[j] = self.cost;
        s.wf[j] = self.extended_weight;
        s.df[j] = self.extended_cost;
    }
}

fn prices(s: &PrimalRowState<'_>, positions: &[usize]) -> Vec<Price> {
    positions.iter().map(|&j| Price::at(s, j)).collect()
}

struct Inputs {
    pos: Vec<usize>,
    fixed: Vec<bool>,
    flags: [u64; 10],
    products: Vec<f64>,
    extended_products: Vec<F80>,
    rho: Vec<f64>,
    u: Vec<f64>,
    tau: Vec<F80>,
    scalars: [F80; 3],
}

impl Inputs {
    fn flags(s: &PrimalRowState<'_>, a: &PrimalRowInput<'_>) -> [u64; 10] {
        [
            a.q as u64,
            a.leave as u64,
            s.excluded.unwrap_or(usize::MAX) as u64,
            u64::from(a.m2),
            u64::from(a.has_mirror),
            u64::from(a.faithful_weights),
            u64::from(a.cross.is_some()),
            a.theta_d.to_bits(),
            a.alpha_rq.to_bits(),
            a.wq.to_bits(),
        ]
    }

    fn scalars(a: &PrimalRowInput<'_>) -> [F80; 3] {
        let (alpha, weight) = a.qw.unwrap();
        [alpha, weight, a.theta_dq.unwrap()]
    }

    fn new(s: &PrimalRowState<'_>, a: &PrimalRowInput<'_>) -> Self {
        Self {
            pos: s.pos_of.to_vec(),
            fixed: s.fixed.to_vec(),
            flags: Self::flags(s, a),
            products: a.prods.unwrap().to_vec(),
            extended_products: a.prods_q.unwrap().to_vec(),
            rho: a.rho.to_vec(),
            u: a.u.to_vec(),
            tau: a.cross.map_or(&[][..], |x| x.tau).to_vec(),
            scalars: Self::scalars(a),
        }
    }

    fn matches(&self, s: &PrimalRowState<'_>, a: &PrimalRowInput<'_>) -> bool {
        self.flags == Self::flags(s, a)
            && extended(&self.scalars, &Self::scalars(a))
            && self.pos == s.pos_of
            && self.fixed == s.fixed
            && doubles(&self.products, a.prods.unwrap())
            && extended(&self.extended_products, a.prods_q.unwrap())
            && doubles(&self.rho, a.rho)
            && doubles(&self.u, a.u)
            && extended(&self.tau, a.cross.map_or(&[][..], |x| x.tau))
    }
}

const MAX_OVERLAYS: usize = 8;

struct Overlay {
    before: Vec<Price>,
    after: Vec<Price>,
    reset: bool,
    flags: u32,
}

struct Plan {
    inputs: Inputs,
    delta: Option<Arc<()>>,
    changed: Vec<usize>,
    row_values: Vec<(usize, f64)>,
    overlays: Vec<Overlay>,
}

impl PricingRows {
    pub(super) fn apply(
        &self,
        s: PrimalRowState<'_>,
        a: &PrimalRowInput<'_>,
        delta: Option<&Arc<()>>,
        run: impl FnOnce(PrimalRowState<'_>, bool) -> (ResultRow, bool),
    ) -> ResultRow {
        // Both admitted paths use the already-supplied products. Normal Devex
        // additionally requires its already-warmed, append-only delta entry.
        #[cfg(test)]
        if audit::original_only() {
            return run(s, false).0;
        }
        let n = s.w.len();
        let supported = n <= 100
            && s.wf.len() == n
            && s.d.len() == n
            && s.df.len() == n
            && a.rho.len() <= 24
            && a.u.len() <= 24
            && a.cross.is_none_or(|x| x.tau.len() <= 24)
            && s.pos_of.len() == n
            && s.fixed.len() == n
            && a.qw.is_some()
            && a.theta_dq.is_some()
            && a.prods.is_some_and(|v| v.len() == n)
            && a.prods_q.is_some_and(|v| v.len() == n)
            && (a.m2 || delta.is_some());
        let Some(flags) = quad_main::status().filter(|_| supported) else {
            return run(s, false).0;
        };
        let mut plan = self.0.borrow_mut();
        let mut before = if let Some(p) = plan.as_mut() {
            let same_delta = match (&p.delta, delta) {
                (None, None) => true,
                (Some(x), Some(y)) => Arc::ptr_eq(x, y),
                _ => false,
            };
            if !same_delta || !p.inputs.matches(&s, a) {
                return run(s, false).0;
            }
            if let Some(v) = p.overlays.iter().find(|v| {
                v.flags & !flags == 0
                    && v.before
                        .iter()
                        .zip(&p.changed)
                        .all(|(price, &j)| price.matches(&s, j))
            }) {
                #[cfg(test)]
                let expected = audit::comparing().then(|| audit::Expected::compute(&s, run));
                let mut s = s;
                for (price, &j) in v.after.iter().zip(&p.changed) {
                    price.restore(&mut s, j);
                }
                #[cfg(test)]
                if let Some(expected) = expected {
                    expected.check(&s, v.reset, p);
                }
                return (v.reset, p.changed.clone(), p.row_values.clone());
            }
            if p.overlays.len() == MAX_OVERLAYS {
                return run(s, false).0;
            }
            prices(&s, &p.changed)
        } else {
            (0..n).map(|j| Price::at(&s, j)).collect()
        };
        let state = PrimalRowState {
            w: &mut *s.w,
            d: &mut *s.d,
            wf: &mut *s.wf,
            df: &mut *s.df,
            pos_of: s.pos_of,
            fixed: s.fixed,
            excluded: s.excluded,
        };
        let (result, warm) = run(state, true);
        if warm {
            if let Some(flags) = quad_main::status() {
                let p = plan.get_or_insert_with(|| {
                    before = result.1.iter().map(|&j| before[j]).collect();
                    Plan {
                        inputs: Inputs::new(&s, a),
                        delta: delta.cloned(),
                        changed: result.1.clone(),
                        row_values: result.2.clone(),
                        overlays: Vec::new(),
                    }
                });
                // The ordered row support and coefficients only depend on the plan's
                // checked structural inputs. Keep the scalar path if those ever differ.
                if p.changed == result.1
                    && p.row_values.len() == result.2.len()
                    && p.row_values
                        .iter()
                        .zip(&result.2)
                        .all(|((j, x), (k, y))| j == k && x.to_bits() == y.to_bits())
                {
                    p.overlays.push(Overlay {
                        before,
                        after: prices(&s, &p.changed),
                        reset: result.0,
                        flags,
                    });
                }
            }
        }
        result
    }
}

#[cfg(all(test, target_arch = "x86_64", target_os = "linux"))]
mod tests {
    use super::*;
    use crate::simplex::{build, frozen_primal_row_update, primal_row_update, NONE};
    use std::cell::Cell;

    #[test]
    fn overlays_preserve_untouched_prices_and_bound_storage() {
        let model =
            crate::mps::parse(include_str!("../tests/fixtures/fx_quad_obj_5371.mps")).unwrap();
        let p = build(&model);
        let n = p.n + p.m;
        let mut products = vec![0.0; n];
        products[0] = 2.0;
        let scaled_products: Vec<_> = products.iter().copied().map(F80::from_f64).collect();
        let direction = vec![0.; p.m];
        let input = PrimalRowInput {
            p: &p,
            q: NONE,
            leave: NONE,
            m2: true,
            has_mirror: true,
            faithful_weights: true,
            prods: Some(&products),
            prods_q: Some(&scaled_products),
            rho: &direction,
            u: &direction,
            qw: Some((F80::ONE, F80::ONE)),
            cross: None,
            theta_d: 1.,
            alpha_rq: 1.,
            wq: 1.,
            theta_dq: Some(F80::ONE),
        };
        let rows = PricingRows::default();
        let executed = Cell::new(0);
        let positions = vec![NONE; n];
        let fixed = vec![false; n];
        for (round, value) in (1..=8).chain([1, 9, 9, 2]).enumerate() {
            let mut w = vec![1.; n];
            w[0] = value as f64;
            w[1] = f64::from_bits(0x7ff8000000000030 + round as u64);
            let mut d = vec![0.; n];
            d[0] = value as f64;
            d[1] = if round % 2 == 0 { -0. } else { 0. };
            let mut wf = w.iter().copied().map(F80::from_f64).collect::<Vec<_>>();
            let mut df = d.iter().copied().map(F80::from_f64).collect::<Vec<_>>();
            let (mut w0, mut d0, mut wf0, mut df0) = (w.clone(), d.clone(), wf.clone(), df.clone());
            macro_rules! state {
                ($w:ident,$d:ident,$wf:ident,$df:ident) => {
                    PrimalRowState {
                        w: &mut $w,
                        d: &mut $d,
                        wf: &mut $wf,
                        df: &mut $df,
                        pos_of: &positions,
                        fixed: &fixed,
                        excluded: None,
                    }
                };
            }
            let reference = frozen_primal_row_update(state!(w0, d0, wf0, df0), &input);
            let result = rows.apply(state!(w, d, wf, df), &input, None, |s, _| {
                executed.set(executed.get() + 1);
                (primal_row_update::<true>(s, &input), true)
            });
            assert_eq!(reference, result);
            assert!(
                doubles(&w, &w0) && doubles(&d, &d0) && extended(&wf, &wf0) && extended(&df, &df0)
            );
        }
        assert_eq!(executed.get(), 10);
        assert_eq!(rows.0.borrow().as_ref().unwrap().overlays.len(), 8);
    }
}

#[cfg(test)]
#[path = "quad_row_audit.rs"]
pub(super) mod audit;
