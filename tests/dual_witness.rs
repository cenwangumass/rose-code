use rose::{mps, solve, SolveOptions, Status};

#[test]
fn fresh_factor_cannot_absorb_at_a_nonzero_level_ends_unbounded() {
    // `fx_witness_dual_cc1ca` (branch `0xcc1ca`): the reference ends Unbounded (status 5) after 2 iterations, through the status -1004 exit and the primal.
    let model = mps::parse(include_str!("fixtures/fx_witness_dual_cc1ca.mps")).unwrap();
    let r = solve(&model, None, &SolveOptions::dual_release()).unwrap();
    assert_eq!((r.status, r.iterations), (Status::Unbounded, 2));
}
