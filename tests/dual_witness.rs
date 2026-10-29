use rose::{mps, solve, SolveOptions, Status};

#[test]
fn fresh_factor_cannot_absorb_at_a_nonzero_level_ends_unbounded() {
    // `fx_witness_dual_cc1ca` (branch `0xcc1ca`): the reference ends Unbounded (status 5) after 2 iterations, through the status -1004 exit and the primal.
    let model = mps::parse(include_str!("fixtures/fx_witness_dual_cc1ca.mps")).unwrap();
    let r = solve(&model, None, &SolveOptions::dual_release()).unwrap();
    assert_eq!((r.status, r.iterations), (Status::Unbounded, 2));
}

#[test]
fn quad_third_restart_starts_from_the_slack_basis() {
    // `fx_witness_dual_c6007` (branch `0xc6007`, rejected first row): the reference ends Infeasible after 14 iterations. The third restart of the quad pass has no crash
    // (one crash in the whole run), so it starts from the slack basis and makes three pivots.
    let model = mps::parse(include_str!("fixtures/fx_witness_dual_c6007.mps")).unwrap();
    let r = solve(&model, None, &SolveOptions::dual_release()).unwrap();
    assert_eq!((r.status, r.iterations), (Status::Infeasible, 14));
    let c = solve(&model, None, &SolveOptions::dual()).unwrap();
    assert_eq!((c.status, c.iterations), (Status::Infeasible, 3));
}

#[test]
fn controlled_dual_primal_end_refactors_a_fresh_factor_once() {
    // `fx_witness_dual_cc1ca` under the controlled profile: the reference primal started by the dual finds no entering variable on a fresh factor, requests a
    // refactorization (`0x1452ff`), and then prices the variable the dual left at an infinite upper bound: 3 iterations, Optimal.
    let model = mps::parse(include_str!("fixtures/fx_witness_dual_cc1ca.mps")).unwrap();
    let c = solve(&model, None, &SolveOptions::dual()).unwrap();
    assert_eq!((c.status, c.iterations), (Status::Optimal, 3));
}
