//! The natural dual residual-scale threshold transition (`reverse-engineer/evidence/witness_dual_r3_residual_proof.md`): the structural start basis x2,x3,x8,x9,x10,x11
//! (the six first columns), then the objective `c`, `-c`, `c`, each solve from the previous final basis, under the dual-release profile (Method=1). In the third solve at
//! iteration 6 the first basic variable x3 = 1.0000097725718222 lies above its upper bound 0.999999; the reference takes the low arm of the residual-scale gate
//! (`0xc2025..0xc2044`: request 4, level 0, nonempty infeasible list, threshold below 0.5) three times and refactors at the thresholds 1/128, 1/32, 1/8 and 1/2 before
//! the fifth pivot of that solve (row 5, variable x3 leaves for x8). The low arm needs only the nonempty infeasible list, not the tiny-infeasibility test of the acceptance arm.
//! Reference: iterations 4, 3, 7; every pivot, objective, x, pi and rc bit for bit (capture `replay/gdb_dual_trace_multi.py` with `tools/witness_dual_r3_seq.c`).
use rose::{mps, solve, Basis, SolveOptions, Status};

#[test]
fn residual_low_threshold_ladder_in_the_third_solve() {
    let mut model = mps::parse(include_str!(
        "fixtures/witness_dual_r3_residual_low_small.mps"
    ))
    .unwrap();
    let (n, m) = (model.num_cols(), model.num_rows());
    assert_eq!((n, m), (16, 6));
    let mut basis = Basis {
        basic: (0..m).collect(),
        at_upper: vec![false; n + m],
    };
    let want: [(usize, u64, [u64; 16], [u64; 6], [u64; 16]); 3] = [
        (
            4,
            0xc014e4cf283c14c4,
            [
                0x3ff0000158437f14,
                0x3feffffde7210be9,
                0x4000000000000000,
                0x3f9263bd6a24f661,
                0x3ff7444e4c4307d6,
                0x0000000000000000,
                0x0000000000000000,
                0xbe9590d6b0c0db14,
                0x0000000000000000,
                0x0000000000000000,
                0x0000000000000000,
                0x3fcc484b78742668,
                0x0000000000000000,
                0x3fe906b04359d6c9,
                0x0000000000000000,
                0x0000000000000000,
            ],
            [
                0x41f19f9e9a8939b6,
                0xc1efb29d8c60de16,
                0x41da7bb2a8bcaef4,
                0xc1df2d2c8efa83b5,
                0x41fe5fe8b2d2d770,
                0xc201d883e7b7fc48,
            ],
            [
                0x0000000000000000,
                0xbfe6ce3200000000,
                0xbffd5a2400000000,
                0x0000000000000000,
                0x0000000000000000,
                0x4008e24100000000,
                0x400ea11e80000000,
                0x0000000000000000,
                0x40053a0700000000,
                0x4011065140000000,
                0x4004524f80000000,
                0x0000000000000000,
                0x400cd52e00000000,
                0x0000000000000000,
                0x40122dac80000000,
                0x4012489880000000,
            ],
        ),
        (
            3,
            0xbfe89220031390ea,
            [
                0x3ffbfaeeab523fbd,
                0x0000000000000000,
                0x3ff0000144b085df,
                0x3feffffe0d94d4ba,
                0x3ff00000b5529264,
                0x3feffffe77271c2f,
                0x3fb64eecc7adef55,
                0x0000000000000000,
                0x0000000000000000,
                0x0000000000000000,
                0x0000000000000000,
                0x0000000000000000,
                0x0000000000000000,
                0x0000000000000000,
                0x0000000000000000,
                0x0000000000000000,
            ],
            [
                0x4216f1e6588f0ca9,
                0xc214a268a8c3b58a,
                0x41f900b73fc0c6e7,
                0xc1fd6f07ae7c3902,
                0x41c0fb0d9e6d49cb,
                0xc1c3f3f9e66d7151,
            ],
            [
                0x0000000000000000,
                0x3ff0bdac00000000,
                0x0000000000000000,
                0x0000000000000000,
                0x0000000000000000,
                0x0000000000000000,
                0x0000000000000000,
                0x401afcc500000000,
                0x4026ce5e00000000,
                0x402547f100000000,
                0x400f079900000000,
                0x40209aeac0000000,
                0x3fd5082800000000,
                0x400fd16820000000,
                0x4006834870000000,
                0x3ff1a99180000000,
            ],
        ),
        (
            7,
            0xc014e4cf283c14c4,
            [
                0x3ff0000158437f14,
                0x3feffffde7210be9,
                0x4000000000000000,
                0x3f9263bd6a24f661,
                0x3ff7444e4c4307d6,
                0x0000000000000000,
                0x0000000000000000,
                0xbe9590d6b0c0db14,
                0x0000000000000000,
                0x0000000000000000,
                0x0000000000000000,
                0x3fcc484b78742668,
                0x0000000000000000,
                0x3fe906b04359d6c9,
                0x0000000000000000,
                0x0000000000000000,
            ],
            [
                0x41f19f9e9a8939b6,
                0xc1efb29d8c60de16,
                0x41da7bb2a8bcaef4,
                0xc1df2d2c8efa83b5,
                0x41fe5fe8b2d2d770,
                0xc201d883e7b7fc48,
            ],
            [
                0x0000000000000000,
                0xbfe6ce3200000000,
                0xbffd5a2400000000,
                0x0000000000000000,
                0x0000000000000000,
                0x4008e24100000000,
                0x400ea11e80000000,
                0x0000000000000000,
                0x40053a0700000000,
                0x4011065140000000,
                0x4004524f80000000,
                0x0000000000000000,
                0x400cd52e00000000,
                0x0000000000000000,
                0x40122dac80000000,
                0x4012489880000000,
            ],
        ),
    ];
    for (k, (iters, obj, x, pi, rc)) in want.into_iter().enumerate() {
        if k > 0 {
            model.obj.iter_mut().for_each(|c| *c = -*c);
        }
        let r = solve(&model, Some(&basis), &SolveOptions::dual_release()).unwrap();
        assert_eq!(
            (r.status, r.iterations, r.objective.to_bits()),
            (Status::Optimal, iters, obj),
            "solve {k}"
        );
        let bits = |v: &[f64]| v.iter().map(|a| a.to_bits()).collect::<Vec<u64>>();
        assert_eq!(bits(&r.x), x, "x, solve {k}");
        assert_eq!(bits(&r.dual), pi, "pi, solve {k}");
        assert_eq!(bits(&r.reduced_cost), rc, "rc, solve {k}");
        basis = r.basis;
    }
}
