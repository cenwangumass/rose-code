//! Free-variable cost shift witnesses (`reverse-engineer/evidence/witness-dual-round4.md`, `witness-dual-round20.md`): the public start basis x0,x2,x3 under the
//! dual-release profile. The dual pivots once and ends with the shift flag raised; the scaled primal finds no blocking row for the free x4 (alpha floor `ctx+0x398`
//! = 3e-8, status -6), the end test sees an unscaled dual infeasibility of 1.9e-6 and the unscaled second pass, limited to `T + 0.2 max(T, cols/3 + rows)`, makes one
//! pivot that is undone but counted. Reference: Optimal in 2 iterations, objective, x and rc bit for bit.
use rose::{mps, solve, Basis, SolveOptions, Status};

fn run(src: &str) -> rose::SolveResult {
    let model = mps::parse(src).unwrap();
    let (n, m) = (model.num_cols(), model.num_rows());
    let basis = Basis {
        basic: vec![0, 2, 3],
        at_upper: vec![false; n + m],
    };
    solve(&model, Some(&basis), &SolveOptions::dual_release()).unwrap()
}

#[test]
fn free_shift_r4_second_pass() {
    let r = run(include_str!("fixtures/witness_dual_r4_free_shift.mps"));
    assert_eq!((r.status, r.iterations), (Status::Optimal, 2));
    assert_eq!(r.objective.to_bits(), 0.7499999390401881f64.to_bits());
    assert_eq!(r.x[1].to_bits(), 0.9999999187202508f64.to_bits());
    assert_eq!(
        r.reduced_cost[4].to_bits(),
        (-1.9073486328125e-6f64).to_bits()
    );
}

#[test]
fn free_shift_r20_second_pass() {
    let r = run(include_str!("fixtures/witness_dual_r20_free_second.mps"));
    assert_eq!((r.status, r.iterations), (Status::Optimal, 2));
    assert_eq!(r.objective.to_bits(), (-1.0000000379295306f64).to_bits());
    assert_eq!(r.x[1].to_bits(), 0.9999999187202508f64.to_bits());
    assert_eq!(
        r.reduced_cost[4].to_bits(),
        (-5.7220458984375e-6f64).to_bits()
    );
}
