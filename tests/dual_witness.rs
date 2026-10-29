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

#[test]
fn quad_dual_reads_out_the_long_double_duals() {
    // `fx_witness_dual_persistent_shift_71000` (dual-release, ill-conditioned basis, entries up to 4e6): the quad dual's `pi` is the long-double `y` of the reduced-cost
    // recompute; the double BTRAN gave duals 0.035% off the reference (and a reduced-cost error of 2e5 in the independent check).
    let model = mps::parse(include_str!("fixtures/fx_witness_dual_persistent_shift_71000.mps")).unwrap();
    let r = solve(&model, None, &SolveOptions::dual_release()).unwrap();
    assert_eq!((r.status, r.iterations), (Status::Optimal, 75));
    let pi = [727987659993.3921, -435.33547769125465, -18199685.544660836, -243054.8584814422, 1213314820.2932222, 3253.7126029372016];
    assert_eq!(r.dual, pi);
}

#[test]
#[ignore = "15000 rows, about 20 s; reads the fixture from ../reverse-engineer"]
fn mode_two_pricing_switch_reaches_the_reference_objective() {
    // `witness_dual_r2_mode_switch` (300 SC50A blocks): the reference switches the Devex mode 1 -> 2 at iteration 6011 (`0xf7bb2`, P50/P58 = 31999/97798) and ends
    // Optimal after 11393 iterations at -19372.523117569304; decisions agree on all 11395 pricing calls.
    let text = std::fs::read_to_string("../reverse-engineer/corpus/fixtures/witness_dual_r2_mode_switch.mps").unwrap();
    let r = solve(&mps::parse(&text).unwrap(), None, &SolveOptions::dual_release()).unwrap();
    assert_eq!((r.status, r.iterations, r.objective), (Status::Optimal, 11393, -19372.523117569304));
}
