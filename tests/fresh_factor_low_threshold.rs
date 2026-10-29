//! The public reoptimization sequence of `witness_primal_r2_fresh_low` (`reverse-engineer/evidence/witness-primal-round2.md`): a structural start basis
//! (the sixteen first columns), then the objective `c`, `-c`, `c`, each solve from the previous final basis. The reference ends the third solve after the
//! fresh-factor low-threshold rejection (`0x1453c7 -> 0x145403`: zero updates, threshold 1/128, residual 1.52587890625e-5, then threshold 0.5 and request 6).
//! Reference: `compare.py seq` (iterations 11, 14, 14; every pivot, objective, x, pi and rc bit for bit).
use rose::{mps, solve, Basis, SolveOptions, Status};

#[test]
fn fresh_factor_low_threshold_rejection_in_the_third_solve() {
    let mut model = mps::parse(include_str!("fixtures/witness_primal_r2_fresh_low.mps")).unwrap();
    let (n, m) = (model.num_cols(), model.num_rows());
    assert_eq!((n, m), (48, 16));
    let mut basis = Basis {
        basic: (0..m).collect(),
        at_upper: vec![false; n + m],
    };
    let want = [
        (11, 0xc00a662baa81be2c_u64),
        (14, 0xc03305bbd5fe13c8),
        (14, 0xc00a662c8d1899a8),
    ];
    let mut x = Vec::new();
    for (k, (iters, obj)) in want.into_iter().enumerate() {
        if k > 0 {
            model.obj.iter_mut().for_each(|c| *c = -*c);
        }
        let r = solve(&model, Some(&basis), &SolveOptions::release()).unwrap();
        assert_eq!((r.status, r.iterations, r.objective.to_bits()), (Status::Optimal, iters, obj), "solve {k}");
        basis = r.basis;
        x = r.x;
    }
    // The third solve's x (the rejection refactors at threshold 0.5; without it the objective stays -3.2998879738659266 of the first solve).
    let x_bits: Vec<u64> = x.iter().map(|v| v.to_bits()).collect();
    let want_x: [u64; 48] = [
        0x3ff9e7410807e2ce, 0x0000000000000000, 0x3fefffee74832a20,
        0x3ff0000a3f4f64e4, 0x3ffae9a4949d0ba4, 0x0000000000000000,
        0x3ff38be163395b71, 0x0000000000000000, 0x3ff9adb9ce568c89,
        0x0000000000000000, 0x3ff7444e85c9e37f, 0x0000000000000000,
        0x3ffef710a3115bc4, 0x0000000000000000, 0x3ff41a36b67b3ba2,
        0x0000000000000000, 0x0000000000000000, 0x3fe3ac80d815e453,
        0x0000000000000000, 0x0000000000000000, 0x0000000000000000,
        0x0000000000000000, 0x0000000000000000, 0x0000000000000000,
        0x0000000000000000, 0x3fdf002b7f07ecf3, 0x0000000000000000,
        0x0000000000000000, 0x0000000000000000, 0x3fd14b6e96da1030,
        0x0000000000000000, 0x0000000000000000, 0x3ff486d58e52039c,
        0x0000000000000000, 0x0000000000000000, 0x0000000000000000,
        0x0000000000000000, 0x3fe906af83794dda, 0x0000000000000000,
        0x0000000000000000, 0x0000000000000000, 0x0000000000000000,
        0x0000000000000000, 0x3fdf6371cbb049c1, 0x3fd20249254dcbc9,
        0x0000000000000000, 0x0000000000000000, 0x0000000000000000,
    ];
    assert_eq!(x_bits, want_x);
}
