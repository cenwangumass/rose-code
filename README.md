# rose

A standalone LP solver (continuous LPs, no external solver): a primal and a dual simplex whose decisions were recovered from a reference
implementation (analysis and specifications: `reverse-engineer/` of the parent repository).

## Command line

    cargo build --release
    target/release/rose solve model.mps [--method primal|dual] [--profile controlled|release] [--output result.json]
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

`--method` (default `primal`) picks the simplex method, `--profile` (default `controlled`) the reference configuration it follows:

| `--method` | `--profile controlled` (default) | `--profile release` |
|---|---|---|
| `primal` | `SolveOptions::default()`: no scaling, `SimplexPricing=1`, no perturbation | `SolveOptions::release()`: scaling, initial `SimplexPricing=0` pricing |
| `dual` | `SolveOptions::dual()`: no scaling, no perturbation | `SolveOptions::dual_release()`: scaling and cost perturbation `0.0002` (the library defaults) |

The controlled profiles match the reference's decisions on all twelve corpus models; the release profiles match on the ten smaller ones (primal: all
twelve) and are correct but not decision-identical on the dual's 80BAU3B and FIT2P. Pick `release` for the reference's default behaviour, `controlled`
(the default) for the best-verified path. The option order on the command line does not matter.

Library: `rose::{mps::parse, solve, SolveOptions, LpModel, Basis, Status}`; `SolveOptions::default()` (primal, controlled profile), `SolveOptions::dual()`,
`SolveOptions::dual_release()` and `SolveOptions::release()` select the reference profiles recovered so far.

## What the dual method supports

Cold start (a crash basis that is dual feasible by choice of bounds, artificial bounds of 1e30 for variables without the bound they need: the
reference's phase 1), warm start from any basis, equality, inequality, ranged and free rows, boxed, one-sided, fixed and free variables, both senses,
objective offsets, bound flipping in the ratio test, dual steepest edge pricing, cost shifts, refactorization by the reference's work-based trigger.
Statuses: optimal; infeasible with a Farkas certificate (the BTRAN row of the failed ratio test); unbounded with a recession ray and a feasible
point (the dual is infeasible; a zero-objective dual run proves the primal feasible, otherwise the status is infeasible); iteration limit (finite
point, not feasible); numerical failure is reported by the primal path only (a singular factorization during the dual raises an error).

Verified: all twelve Netlib models (AFIRO, SC50A, BLEND, ADLITTLE, BORE3D, BOEING2, DEGEN2, DEGEN3, PILOT4, PEROLD, 80BAU3B, FIT2P) pass the independent
exact checker (`reverse-engineer/tools/check_rose.py`) under both dual profiles and reach the reference's objective of the same profile within
`1e-7 + 1e-8 max(|z|)` (`reverse-engineer/benchmarks/README.md`); the controlled dual decisions equal the reference's on all twelve, the release
dual decisions on the ten smaller models (`tests/dual_traces.rs` pins AFIRO, SC50A, BLEND, ADLITTLE from versioned fixtures); random-model agreement with
the primal (`tests/dual_random.rs`), the metamorphic checks, the general LP tests and the fixture corpus (`tests/fixture_corpus.rs`: infeasible,
unbounded, near-dependent, badly scaled, tied, warm-started and singular-basis fixtures) run under the primal and both dual profiles.

## Recovery gaps (dual)

Details and disassembly addresses: `reverse-engineer/specs/dual-simplex.md`, "Unresolved recovery gaps of the dual path".

- Dual-release on 80BAU3B and FIT2P: correct results, but pivot sequences and iteration counts differ from the reference's (5154 against 5185 and
  8203 against 7183 iterations).
- Branches no corpus model exercises (free variables in the cost-shift code, the single-candidate and relaxed-fallback paths, the pivot-floor failure,
  the forced-row stack, the type-0 pricing fallback, the deferred second pass of the free-column pivots, the perturbation branches with
  `ctx+0x338 < 2`, status -1004): implemented from the disassembly or not at all.
- Full floating-point state agrees bit for bit only on AFIRO and SC50A (ulp-level differences elsewhere); the work counter is not exact at a few
  sites, so free-running refactorization schedules reach the reference's iteration counts on the smaller models without being guaranteed to elsewhere.
- Variables with an infinite lower bound, the reference's `INF_OR_UNBD` status and the zero-objective feasibility step that separates `unbounded`
  from `infeasible` are rose's own design, checked by the independent checker and the random tests, not by reference traces.
- A singular factorization in the middle of a dual run is an error, not a recovery (no refactorization with replaced columns).
- Four singular supplied-basis fixtures and the warm starts from non-optimal bases (AFIRO 14/13, SC50A 33/34, ADLITTLE 65/59 controlled; rose/reference)
  differ from the reference in iteration count only; results agree. The reference reports 0 to 3 iterations on the tiny unbounded fixtures, rose 0 to 2.
- Mixed-integer models, presolve and crossover are out of scope.

## Tests

    cargo test --release

Fast (about 10 s on a quiet machine); the large reference traces are `#[ignore]`d (`cargo test --release -- --ignored`).
