# Performance state

Current state of the work in [PLAN.md](PLAN.md). Rewrite in place; no log (history is `jj log`).

## Numbers

Baseline solver: `8b645007` (same solver as `vlsokskq`). Current solver: `kyyvxyql` (dense LU membership flags and reusable seen generations).
Parity: full tier 987/987; release tests pass. Three timed runs after one warm-up, hyperfine 2.0.0; instruction counts are millions from perf stat.

`baseline wall` is the fixed initial process-wall median. `paired base wall` is the fresh baseline-binary median from alternating baseline/last-landed/candidate checks on one CPU (CPU 0 unless stated otherwise); it resolves host drift without changing the fixed baseline. `current wall` and `solve` use the new harness; **current/ref uses solve-only time**, `diagnostics.seconds`, against original Gurobi Runtime.

Reference provenance: original-library cached records `../reverse-engineer/runs/{degen2,degen3,perold,pilot4}-{profile}.json` for the ten previously valid cells. New three-run medians use `/tmp/rose-perf-reference/{primal_5959.primal-r,primal_5959.dual-r,quad_2717.dual-r,fx_quad_obj_5371.dual-r}.{1,2,3}.json`; summary `/tmp/rose-perf-reference/summary.json`. All original-library records identify Ryzen 9 5900X and library SHA256 `8e87fca221fdffc922331c8619426031b3e72dfd419446422f88af0cb0e98f9c`.
`-` refs are unreachable: original Gurobi rejects 80BAU3B, GREENBEA and mode-switch with license error 10010 (scratch `ID.error.json`, including direct 80BAU3B primal/dual attempts). The old 80BAU3B primal/dual refs used the patched library and have been removed; the old prose mode-switch time is unusable.

Regression checks: all 23 runs paired against the preceding solver on CPU 4. Apparent mode-switch and 80BAU3B barrier slowdowns resolve on CPU 2: mode-switch old/new wall 11.781/12.060 s (spreads 4.056/4.844 s), instructions 79720.2/79564.6M; barrier 0.387/0.386 s (spreads 0.012/0.014 s), unchanged instructions. Frozen-baseline audit has no new unresolved differences. Raw data `target/perf-members-full-{before,after}.txt`, `target/perf-members-mode-{before,after}.txt`, `target/perf-members-barrier-{before,after}.txt`.
Acceptance: primal_5959 primal full-sweep solve 0.790302 -> 0.724951 s (8.3%; spreads 0.024432/0.010984 s), instructions 4683.5 -> 4128.9M (-11.8%). Five-run alternating pair also clears spread: 0.784695 -> 0.721603 s (`target/perf-members-paired.txt`). Quad-objective solve 0.295456 -> 0.260088 s (12.0%), instructions 2741.2 -> 2344.8M (-14.5%). Only membership storage/lookup changes: numeric traversal and work estimates are unchanged.

| run | ref solve | baseline wall | paired base wall | current wall | solve | current/ref | instr |
|---|---|---|---|---|---|---|---|
| degen2.primal-c | 0.012711 | 0.182 | - | 0.094 | 0.089679 | 7.06 | 487.4 |
| degen2.dual-r | 0.009653 | 0.125 | - | 0.079 | 0.074861 | 7.76 | 364.2 |
| degen3.primal-c | 0.163137 | 3.237 | - | 1.504 | 1.482746 | 9.09 | 7779.5 |
| degen3.primal-r | 0.193685 | 3.872 | - | 1.089 | 1.078679 | 5.57 | 9628.4 |
| perold.primal-c | 0.160994 | 2.430 | - | 0.875 | 0.869745 | 5.40 | 7606.2 |
| pilot4.primal-r | 0.015787 | 0.253 | - | 0.090 | 0.086257 | 5.46 | 830.3 |
| 80bau3b.primal-c | - | 1.618 | - | 0.964 | 0.949674 | - | 10155.3 |
| primal_5959.primal-r | 0.067884 | 0.563 | 0.938 | 0.728 | 0.724951 | 10.68 | 4128.9 |
| degen3.dual-c | 0.082843 | 1.002 | - | 0.769 | 0.753034 | 9.09 | 3856.6 |
| degen3.dual-r | 0.082947 | 1.008 | 1.007 | 0.772 | 0.754825 | 9.10 | 3866.7 |
| perold.dual-r | 0.038532 | 0.355 | 0.360 | 0.345 | 0.337860 | 8.77 | 1845.8 |
| 80bau3b.dual-r | - | 1.103 | 1.688 | 1.163 | 1.130513 | - | 5307.2 |
| greenbea.dual-c | - | 2.870 | 2.363 | 2.043 | 2.022257 | - | 11007.5 |
| greenbea.dual-r | - | 3.598 | 3.592 | 1.337 | 1.322079 | - | 10915.3 |
| mode_switch.dual-r | - | 25.815 | - | 12.060 | 11.851636 | - | 79564.6 |
| primal_5959.dual-r | 0.098192 | 1.300 | - | 0.646 | 0.643861 | 6.56 | 5874.1 |
| quad_2717.dual-r | 0.092814 | 1.130 | - | 0.589 | 0.586502 | 6.32 | 4991.2 |
| fx_quad_obj_5371.dual-r | 0.038077 | 0.334 | 0.531 | 0.262 | 0.260088 | 6.83 | 2344.8 |
| degen3.barrier-r | 0.075214 | 0.608 | 0.856 | 0.620 | 0.578966 | 7.70 | 5500.7 |
| 80bau3b.barrier-r | - | 0.379 | 0.549 | 0.386 | 0.362489 | - | 4057.3 |
| greenbea.barrier-r | - | 1.569 | 2.493 | 2.219 | 2.168605 | - | 16052.5 |
| 80bau3b.pdhg | - | 4.552 | 6.526 | 5.349 | 5.326934 | - | 92389.4 |
| greenbea.pdhg | - | 2.801 | - | 1.912 | 1.900185 | - | 37740.6 |

## Hotspots

Per profiled bench run: the commit it was taken at and the top self-time functions (percent), at most
five lines each. Re-profile a run when a landed change touched its hot path.

- primal_5959.primal-r, `kyyvxyql` (profiling build, three bounded solves at 1999 Hz, `target/perf-members-primal5959.data`): F80 add_signed 6.76%, primal control 5.87%, F80 nucleus 5.65%, fresh F80 factor 5.52%, free 4.49% self cycles; malloc 4.07%. Current log `target/perf-members-primal5959.log` verifies numerical_failure, 10780 iterations, 5381 refactorizations. Original reference has status 12 and 10780 iterations, no refactor count. DWARF caller attribution is incomplete despite unwind FDEs; use a validated frame-pointer profile before relying on inclusive shares.
- fx_quad_obj_5371.dual-r, `kyyvxyql` (profiling build, three solves at 1999 Hz, `target/perf-members-fx.data`): F80 add_signed 8.06%, F80 nucleus 5.81%, multiply 5.74%, fresh F80 factor 5.17%, pack 3.83% self cycles. Solve remains infeasible, 5319 iterations, 5 refactorizations.
- degen3.dual-r, `yzqsoomv` (stale after LU membership changes; reprofile before choosing it; profiling build, three solves at 1999 Hz): iteration 10.59%, update 10.17%, L stage 8.00%, BTRAN 6.51%, ftran_with 6.51% self cycles. ftran_needed 5.64% plus term-buffer filling/column lookup 5.10%; the old per-row allocation is gone. Full solve: 2212 iterations, 24 refactorizations.
- mode_switch.dual-r, `xzuktkmx` (stale after F80/LU changes; reprofile before choosing it): update 22.21%, LU nucleus 15.90%, BTRAN 13.33%, fresh-factor row-list filtering 7.67%, factor_structure 3.36% self cycles. Full solve: 11393 iterations, 17 refactorizations.

## In progress

(none)

## Candidates

Ordered by expected payoff; one line each, with the bench run it targets. Start from PLAN.md "Where to
look first" and replace these with measured items after the first profiles.

- primal_5959.primal-r (10.68x, highest reachable ratio): establish full caller attribution (try frame-pointer profiling and parity-check that profiling binary) before choosing between F80 add, primal control and repeated factor/nucleus setup. Counts still match the reference's 10780 iterations; 5381 current refactorizations make allocation/recomputation relevant leads.
- DEGEN3 dual-r (9.10x), dual-c (9.09x), primal-c (9.09x), PEROLD dual-r (8.77x): follow in decreasing ratio; reprofile changed LU paths first.
- fx_quad_obj_5371.dual-r (6.83x): F80 arithmetic and factor construction remain hot; numeric rewrites require a frozen-kernel differential test first.
- Other reachable runs follow their table ratios. DEGEN3 row-dot filling/lookup remains an untried lead after buffer reuse; preserve accumulation order and work estimates.
- Runs without reachable references, after the reachable ratios: mode-switch update (22.21%) still scans A/B lists for other memberships and row-capacity accounting; BTRAN (13.33%) allocates scratch/recomputes L levels; LU nucleus (15.90%); fresh-factor row-list filtering (7.67%). Annotate before choosing a change.
- all: dump::on/record construction without dumps and allocation in FTRAN/BTRAN, pricing and ratio tests remain unmeasured leads.

## Tried, no gain

- Read eliminate work counts directly from maintained rowl lengths instead of rebuilding counts: parity full 987/987 and release tests pass, but mode_switch.dual-r 22.849 -> 21.851 s is smaller than the old 1.849 s spread (instructions 146886.1 -> 140945.4M, -4.0%). Alternating DEGEN2 primal 0.092/0.090 then 0.091/0.091 s; dual 0.068/0.067 s, marginal against spread. Reverted.

## Tried, breaks parity

(idea, first differing record from `rose-parity check`, why)
