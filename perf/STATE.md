# Performance state

Current state of the work in [PLAN.md](PLAN.md). Rewrite in place; no log (history is `jj log`).

## Numbers

Baseline solver: `8b645007` (same solver as `vlsokskq`). Current solver: `uprrtwlk` (build identical quad factor twins once, then deep-clone).
Parity: full tier 987/987; release tests pass. Three timed runs after one warm-up, hyperfine 2.0.0; instruction counts are millions from perf stat.

`baseline wall` is the fixed initial process-wall median. `paired base wall` is the fresh baseline-binary median from alternating baseline/last-landed/candidate checks on one CPU (CPU 0 unless stated otherwise); it resolves host drift without changing the fixed baseline. `current wall` and `solve` use the new harness; **current/ref uses solve-only time**, `diagnostics.seconds`, against original Gurobi Runtime.

Reference provenance: original-library cached records `../reverse-engineer/runs/{degen2,degen3,perold,pilot4}-{profile}.json` for the ten previously valid cells. New three-run medians use `/tmp/rose-perf-reference/{primal_5959.primal-r,primal_5959.dual-r,quad_2717.dual-r,fx_quad_obj_5371.dual-r}.{1,2,3}.json`; summary `/tmp/rose-perf-reference/summary.json`. All original-library records identify Ryzen 9 5900X and library SHA256 `8e87fca221fdffc922331c8619426031b3e72dfd419446422f88af0cb0e98f9c`.
`-` refs are unreachable: original Gurobi rejects 80BAU3B, GREENBEA and mode-switch with license error 10010 (scratch `ID.error.json`, including direct 80BAU3B primal/dual attempts). The old 80BAU3B primal/dual refs used the patched library and have been removed; the old prose mode-switch time is unusable.

Regression checks: all 23 runs paired against the preceding solver on CPU 4. PEROLD dual's small increase and quad_2717's host-speed shift resolve in five-run pairs on CPU 2: PEROLD old/new solve 0.216386/0.217990 s (spreads 0.003053/0.001343 s), instructions 1846.2/1846.1M; quad_2717 0.583504/0.546952 s, instructions 4993.8/4613.7M. Other increases lie within spread; no new unresolved frozen-baseline regression. Raw data `target/perf-twins-full-{before,after}.txt`, `target/perf-twins-recheck-{before,after}.txt`. Host speed varies substantially across the sweep; compare each immediate pair rather than medians from different sweeps.
Acceptance: primal_5959 primal solve 0.451306 -> 0.388362 s (13.9%; spreads 0.003942/0.007561 s), instructions 4128.8 -> 3521.0M (-14.7%). Primal_5959 dual 0.997306 -> 0.852086 s (-14.6%), instructions 5874.3 -> 4984.1M (-15.2%). CPU-2 quad-objective 0.259571 -> 0.229380 s (-11.6%), instructions 2344.8 -> 1997.5M (-14.8%). Identical fresh quad twins now share construction work through a deep copy; independent mutable factor state, PGATE records and all work estimates are preserved.

| run | ref solve | baseline wall | paired base wall | current wall | solve | current/ref | instr |
|---|---|---|---|---|---|---|---|
| degen2.primal-c | 0.012711 | 0.182 | - | 0.092 | 0.086527 | 6.81 | 487.8 |
| degen2.dual-r | 0.009653 | 0.125 | - | 0.079 | 0.074578 | 7.73 | 364.2 |
| degen3.primal-c | 0.163137 | 3.237 | - | 1.439 | 1.422194 | 8.72 | 7786.6 |
| degen3.primal-r | 0.193685 | 3.872 | - | 1.714 | 1.698625 | 8.77 | 9632.3 |
| perold.primal-c | 0.160994 | 2.430 | - | 1.383 | 1.375523 | 8.54 | 7613.5 |
| pilot4.primal-r | 0.015787 | 0.253 | - | 0.145 | 0.139482 | 8.84 | 831.2 |
| 80bau3b.primal-c | - | 1.618 | - | 0.958 | 0.942793 | - | 10200.0 |
| primal_5959.primal-r | 0.067884 | 0.563 | 0.938 | 0.390 | 0.388362 | 5.72 | 3521.0 |
| degen3.dual-c | 0.082843 | 1.002 | - | 0.500 | 0.489545 | 5.91 | 3856.7 |
| degen3.dual-r | 0.082947 | 1.008 | 1.007 | 0.499 | 0.488976 | 5.90 | 3866.4 |
| perold.dual-r | 0.038532 | 0.355 | 0.360 | 0.223 | 0.217990 | 5.66 | 1846.1 |
| 80bau3b.dual-r | - | 1.103 | 1.688 | 0.583 | 0.567675 | - | 5308.2 |
| greenbea.dual-c | - | 2.870 | 2.363 | 2.336 | 2.306902 | - | 11007.3 |
| greenbea.dual-r | - | 3.598 | 3.592 | 2.015 | 1.995987 | - | 10914.8 |
| mode_switch.dual-r | - | 25.815 | - | 12.112 | 11.887807 | - | 79567.5 |
| primal_5959.dual-r | 0.098192 | 1.300 | - | 0.854 | 0.852086 | 8.68 | 4984.1 |
| quad_2717.dual-r | 0.092814 | 1.130 | - | 0.549 | 0.546952 | 5.89 | 4613.7 |
| fx_quad_obj_5371.dual-r | 0.038077 | 0.334 | 0.531 | 0.231 | 0.229380 | 6.02 | 1997.5 |
| degen3.barrier-r | 0.075214 | 0.608 | 0.856 | 0.645 | 0.604927 | 8.04 | 5499.3 |
| 80bau3b.barrier-r | - | 0.379 | 0.549 | 0.390 | 0.365919 | - | 4057.5 |
| greenbea.barrier-r | - | 1.569 | 2.493 | 1.529 | 1.495408 | - | 16052.5 |
| 80bau3b.pdhg | - | 4.552 | 6.526 | 6.349 | 6.328385 | - | 92390.6 |
| greenbea.pdhg | - | 2.801 | - | 2.797 | 2.779878 | - | 37740.6 |

## Hotspots

Per profiled bench run: the commit it was taken at and the top self-time functions (percent), at most
five lines each. Re-profile a run when a landed change touched its hot path.

- primal_5959.primal-r, `kyyvxyql` (stale after twin reuse; frame-pointer build passes full 987/987 parity; three bounded solves at 1999 Hz, `target/perf-fp-primal5959.data`, report with `--no-inline`): primal control 6.32%, F80 nucleus 5.46%, fresh F80 factor 5.09%, F80 add_signed 6.67% self. Refactor 53.53%, quad mirrors 28.88% inclusive before twin reuse. Log verifies numerical_failure, 10780 iterations, 5381 refactorizations; original reference status 12, 10780 iterations, refactor count unavailable. Default inline expansion duplicates frames; use PLAN.md's validated frame-pointer workflow for caller attribution.
- fx_quad_obj_5371.dual-r, `kyyvxyql` (stale after twin reuse; profiling build, three solves at 1999 Hz, `target/perf-members-fx.data`): F80 add_signed 8.06%, F80 nucleus 5.81%, multiply 5.74%, fresh F80 factor 5.17%, pack 3.83% self cycles. Solve remains infeasible, 5319 iterations, 5 refactorizations.
- degen3.dual-r, `yzqsoomv` (stale after LU membership changes; reprofile before choosing it; profiling build, three solves at 1999 Hz): iteration 10.59%, update 10.17%, L stage 8.00%, BTRAN 6.51%, ftran_with 6.51% self cycles. ftran_needed 5.64% plus term-buffer filling/column lookup 5.10%; the old per-row allocation is gone. Full solve: 2212 iterations, 24 refactorizations.
- mode_switch.dual-r, `xzuktkmx` (stale after F80/LU changes; reprofile before choosing it): update 22.21%, LU nucleus 15.90%, BTRAN 13.33%, fresh-factor row-list filtering 7.67%, factor_structure 3.36% self cycles. Full solve: 11393 iterations, 17 refactorizations.

## In progress

(none)

## Candidates

Ordered by current/ref; reprofile any changed path before choosing a change. Ratios use the current
sweep's solve medians; paired old/new times determine whether a change gains or regresses.

- PILOT4 primal-r (8.84x), DEGEN3 primal-r (8.77x), primal-c (8.72x), primal_5959 dual-r (8.68x), PEROLD primal-c (8.54x): profile in decreasing ratio, starting with the small PILOT4 model. Fresh quad twins are already copied; avoid proposing the duplicate construction again.
- DEGEN3 barrier-r (8.04x), DEGEN2 dual-r (7.73x), primal-c (6.81x), fx_quad_obj_5371 dual-r (6.02x), DEGEN3 dual-c (5.91x), quad_2717 dual-r (5.89x), DEGEN3 dual-r (5.89x), primal_5959 primal-r (5.72x), PEROLD dual-r (5.66x): follow their measured ratios. DEGEN3 row-dot filling/lookup remains an untried lead after buffer reuse; preserve accumulation order and work estimates.
- Runs without reachable references, after the reachable ratios: mode-switch update (22.21%) still scans A/B lists for other memberships and row-capacity accounting; BTRAN (13.33%) allocates scratch/recomputes L levels; LU nucleus (15.90%); fresh-factor row-list filtering (7.67%). These old profiles are stale; annotate again before choosing a change.
- all: dump::on/record construction without dumps and allocation in FTRAN/BTRAN, pricing and ratio tests remain unmeasured leads.

## Tried, no gain

- Read eliminate work counts directly from maintained rowl lengths instead of rebuilding counts: parity full 987/987 and release tests pass, but mode_switch.dual-r 22.849 -> 21.851 s is smaller than the old 1.849 s spread (instructions 146886.1 -> 140945.4M, -4.0%). Alternating DEGEN2 primal 0.092/0.090 then 0.091/0.091 s; dual 0.068/0.067 s, marginal against spread. Reverted.

## Tried, breaks parity

(idea, first differing record from `rose-parity check`, why)
