# rose

A standalone LP solver (continuous LPs, no external solver): a primal and a dual simplex whose decisions were recovered from a reference
implementation (analysis and specifications: `reverse-engineer/` of the parent repository).

## Command line

    cargo build --release
    target/release/rose solve model.mps [--method primal|dual] [--output result.json]
                              [--basis-in basis.json] [--basis-out basis.json] [--max-iterations N]

Fixed and free format MPS (ranged rows, `RANGES`, bounds `UP LO FX FR MI PL`; integer markers and `BV LI UI SC` bounds are rejected). Progress goes to
stderr; the result is JSON in the `--output` file (stdout if absent):

    { "version": 1, "status": "optimal|infeasible|unbounded|iteration_limit|numerical_failure",
      "objective": ..., "x": [...], "dual": [...], "reduced_cost": [...],
      "basis": {"version": 1, "basic": [...], "at_upper": [...]},
      "farkas": [...], "ray": [...],
      "diagnostics": {"iterations": ..., "refactorizations": ..., "seconds": ..., "check": null} }

`objective` is null unless the status is optimal. `farkas` (one multiplier per row) accompanies `infeasible`, `ray` (one entry per column)
`unbounded`. `check` is the in-process `verify::check` of the result against the model (null when it passes). Exit code 0 for optimal, infeasible and
unbounded, 3 for iteration limit and numerical failure, 1 for invalid input, 2 for usage errors. The basis files are the format of `--basis-out`
(`basic[i]` is the variable in basis position i; variables are the columns, then the slacks); a supplied basis may be singular (it is repaired).

Library: `rose::{mps::parse, solve, SolveOptions, LpModel, Basis, Status}`; `SolveOptions::default()` (primal, controlled profile), `SolveOptions::dual()`,
`SolveOptions::dual_release()` and `SolveOptions::release()` select the reference profiles recovered so far.

## What the dual method supports

Cold start (a crash basis that is dual feasible by choice of bounds, artificial bounds of 1e30 for variables without the bound they need: the
reference's phase 1), warm start from any basis, equality, inequality, ranged and free rows, boxed, one-sided, fixed and free variables, both senses,
objective offsets, bound flipping in the ratio test, dual steepest edge pricing, cost shifts, refactorization by the reference's work-based trigger.
Statuses: optimal; infeasible with a Farkas certificate (the BTRAN row of the failed ratio test); unbounded with a recession ray and a feasible
point (the dual is infeasible; a zero-objective dual run proves the primal feasible, otherwise the status is infeasible); iteration limit (finite
point, not feasible); numerical failure is reported by the primal path only (a singular factorization during the dual raises an error).

Verified: the ten Netlib models AFIRO, SC50A, BLEND, ADLITTLE, BORE3D, BOEING2, DEGEN2, DEGEN3, PILOT4, PEROLD pass the independent exact
checker (`reverse-engineer/tools/check_rose.py`) and reach the reference's objective; the dual decisions equal the reference's on all of them
(`tests/dual_traces.rs` pins AFIRO, SC50A, BLEND, ADLITTLE from versioned fixtures); random-model agreement with the primal
(`tests/dual_random.rs`), the metamorphic checks and the general LP tests run under both methods.

## Recovery gaps (dual)

- 80BAU3B and FIT2P are not run in the dual loop (branches that no smaller model exercises are open: see `reverse-engineer/PLAN-dual.md`).
- Free-running schedules reach the reference's iteration counts on the ten models; the work counter is not exact at a few sites, and full floating
  point state agrees bit for bit only on AFIRO and SC50A (ulp-level differences elsewhere).
- The release profile (scaling, default pricing, perturbation) is partly recovered (`SolveOptions::dual_release()`); perturbation is not.
- Warm starts from a basis that is not dual feasible take a few more or fewer iterations than the reference (AFIRO 14/13, SC50A 33/34, ADLITTLE 65/59); the
  reference reports 0 to 3 iterations on the tiny unbounded fixtures where rose needs 0 to 2.
- A singular factorization in the middle of a dual run is an error, not a recovery (no refactorization with replaced columns).
- Mixed-integer models, presolve and crossover are out of scope.

## Tests

    cargo test --release

Fast (about 10 s on a quiet machine); the large reference traces are `#[ignore]`d (`cargo test --release -- --ignored`).
