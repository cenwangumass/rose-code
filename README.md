# rose

A standalone LP solver (continuous LPs, no external solver): a primal and a dual simplex whose decisions were recovered from a reference
implementation (analysis and specifications: `reverse-engineer/` of the parent repository).

## Command line

    cargo build --release
    target/release/rose solve model.mps [--method primal|dual] [--profile controlled|release] [--output result.json]
                              [--basis-in basis.json] [--basis-out basis.json] [--max-iterations N]

On Linux, an optional CLI build can use the system TCMalloc allocator from gperftools:

    cargo build --release --bin rose --features system-tcmalloc
    # PGO build with the same allocator during training and final compilation:
    python3 -B perf/build_pgo.py --system-tcmalloc --check --test

This requires the gperftools development library (`libtcmalloc_minimal` for linking, and its shared library at runtime). The feature affects the CLI; library consumers keep their allocator. Ordinary builds have no TCMalloc dependency.

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

The original twelve models match the reference's iteration counts in all four simplex configurations and pass the independent checker.
Decision traces also agree on the expanded 26-model corpus in both profiles. Generated ill-conditioned witnesses still expose arithmetic and
restart differences; see `reverse-engineer/PLAN.md` and `STATE.md` for current coverage. Pick `release` for the reference's default behaviour,
`controlled` (the default) to disable scaling and perturbation. The option order on the command line does not matter.

Library: `rose::{mps::parse, solve, SolveOptions, LpModel, Basis, Status}`; `SolveOptions::default()` (primal, controlled profile), `SolveOptions::dual()`,
`SolveOptions::dual_release()` and `SolveOptions::release()` select the reference profiles recovered so far.

## What the dual method supports

Cold start (a crash basis that is dual feasible by choice of bounds, artificial bounds of 1e30 for variables without the bound they need: the
reference's phase 1), warm start from any basis, equality, inequality, ranged and free rows, boxed, one-sided, fixed and free variables, both senses,
objective offsets, bound flipping in the ratio test, dual steepest edge pricing, cost shifts, refactorization by the reference's work-based trigger.
Statuses follow the reference's exits: optimal; infeasible with a Farkas certificate; unbounded with a recession ray; iteration limit;
and numerical failure. An unbounded exit establishes an improving direction, but does not establish primal feasibility. A contradictory model
can therefore return `unbounded`, and `x` is not a feasible-point certificate for this status. Rose does not run an extra feasibility solve to
reclassify the reference's result. `verify::check` remains strict: it requires a feasible point as well as a valid ray when checking unboundedness,
so the CLI's `diagnostics.check` can report a failure for a faithfully reproduced reference exit.

Verified: all twelve original Netlib models (AFIRO, SC50A, BLEND, ADLITTLE, BORE3D, BOEING2, DEGEN2, DEGEN3, PILOT4, PEROLD, 80BAU3B, FIT2P)
pass the independent exact checker (`reverse-engineer/tools/check_rose.py`) under both dual profiles, with the reference's iteration counts
and objectives within `1e-7 + 1e-8 max(|z|)`. The expanded 26-model decision traces agree in both profiles. Certificate and ray fixtures,
random models, metamorphic checks and supplied-basis tests run in the release test suite.

## Recovery gaps (dual)

The current recovery plan is `reverse-engineer/PLAN.md`; rules and remaining branch details are in
`reverse-engineer/specs/dual-simplex.md`. Generated ill-conditioned models still differ in extended-precision arithmetic, primal handoffs
and restart control. Full-state and work-counter comparisons of the largest models cover bounded intervals; equal final iteration counts
do not establish every internal transition. Several forced-row and rare numerical branches remain without a witness.
Presolve, crossover, mixed-integer, quadratic and conic models are out of scope.

## Tests

    cargo test --release

Fast (about 10 s on a quiet machine); the large reference traces are `#[ignore]`d (`cargo test --release -- --ignored`).

    target/release/rose-parity check [--tier full]

Self-parity gate: every dump record and result of the runs in `parity/corpus.txt` must equal the frozen baseline's
(`parity/baseline/`). Run it after any change meant to keep the numbers; see `AGENTS.md`.
