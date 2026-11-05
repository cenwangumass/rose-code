//! Test the lower-level Devex kernel and its higher-level numerical result independently.
use super::*;
use std::cell::Cell;

#[derive(Clone, Copy, PartialEq)]
enum Mode {
    Normal,
    Original,
    Compare,
}
thread_local! {
    static MODE: Cell<Mode> = const { Cell::new(Mode::Normal) };
    static VERIFIED: Cell<usize> = const { Cell::new(0) };
}

pub(in crate::simplex) struct Guard(Mode);

impl Drop for Guard {
    fn drop(&mut self) {
        MODE.set(self.0);
    }
}

pub(in crate::simplex) fn exercise_original() -> Guard {
    Guard(MODE.replace(Mode::Original))
}

pub(super) fn original_only() -> bool {
    MODE.get() == Mode::Original
}

pub(super) fn comparing() -> bool {
    MODE.get() == Mode::Compare
}

pub(super) struct Expected {
    weights: Vec<f64>,
    costs: Vec<f64>,
    extended_weights: Vec<F80>,
    extended_costs: Vec<F80>,
    row: ResultRow,
}

impl Expected {
    pub(super) fn compute(
        state: &PrimalRowState<'_>,
        original: impl FnOnce(PrimalRowState<'_>, bool) -> (ResultRow, bool),
    ) -> Self {
        let mut weights = state.w.to_vec();
        let mut costs = state.d.to_vec();
        let mut extended_weights = state.wf.clone();
        let mut extended_costs = state.df.to_vec();
        let inputs = PrimalRowState {
            w: &mut weights,
            d: &mut costs,
            wf: &mut extended_weights,
            df: &mut extended_costs,
            pos_of: state.pos_of,
            fixed: state.fixed,
            excluded: state.excluded,
        };
        // The activity counter below describes the solver under test. Keep
        // the independent oracle's genuine Devex term accesses out of it.
        let activity = crate::simplex::DEVEX_DELTA_ACTIVITY.get();
        let before = quad_main::status();
        let row = original(inputs, false).0;
        assert_eq!(quad_main::status(), before);
        crate::simplex::DEVEX_DELTA_ACTIVITY.set(activity);
        Self {
            weights,
            costs,
            extended_weights,
            extended_costs,
            row,
        }
    }
    pub(super) fn check(&self, state: &PrimalRowState<'_>, reset: bool, plan: &Plan) {
        assert!(doubles(&self.weights, state.w) && doubles(&self.costs, state.d));
        assert!(
            extended(&self.extended_weights, state.wf) && extended(&self.extended_costs, state.df)
        );
        assert_eq!(self.row.0, reset);
        assert_eq!(self.row.1, plan.changed);
        assert_eq!(self.row.2.len(), plan.row_values.len());
        for ((a, x), (b, y)) in self.row.2.iter().zip(&plan.row_values) {
            assert_eq!((*a, x.to_bits()), (*b, y.to_bits()));
        }
        VERIFIED.set(VERIFIED.get() + 1);
    }
}

#[cfg(all(target_arch = "x86_64", target_os = "linux"))]
#[test]
fn complete_compact_rows_match_with_pricing_overlays_active() {
    let _guard = Guard(MODE.replace(Mode::Compare));
    for (name, status, iterations) in [
        ("fx_quad_obj_5371", crate::Status::Infeasible, 5319),
        ("quad_2717", crate::Status::Optimal, 16010),
        ("primal_5959", crate::Status::NumericalFailure, 15688),
    ] {
        let model = crate::mps::parse(
            &std::fs::read_to_string(format!("tests/fixtures/{name}.mps")).unwrap(),
        )
        .unwrap();
        let mut opts = crate::simplex::SolveOptions::dual_release();
        opts.max_iterations = 20000;
        VERIFIED.set(0);
        let result = crate::solve(&model, None, &opts).unwrap();
        assert_eq!((result.status, result.iterations), (status, iterations));
        assert!(
            VERIFIED.get() > 1000,
            "insufficient real higher-level pricing hits on {name}"
        );
    }
}
