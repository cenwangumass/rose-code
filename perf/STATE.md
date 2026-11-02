# Performance state

Current state of the work in [PLAN.md](PLAN.md). Rewrite in place; no log (history is `jj log`).

## Numbers

Baseline solver: `8b645007` (same solver as `vlsokskq`). Current solver: `munmovyt` (reuse fresh compact quad factors for repeated exact ordered inputs).
Parity: full tier 987/987; release tests pass. Three timed runs after one warm-up, hyperfine 2.0.0; instruction counts are millions from perf stat.

`baseline wall` is the fixed initial process-wall median. `paired base wall` is the fresh baseline-binary median from alternating baseline/last-landed/candidate checks on one CPU (CPU 0 unless stated otherwise); it resolves host drift without changing the fixed baseline. `current wall` and `solve` use the new harness; **current/ref uses solve-only time**, `diagnostics.seconds`, against original Gurobi Runtime.

Reference provenance: original-library cached records `../reverse-engineer/runs/{degen2,degen3,perold,pilot4}-{profile}.json` for the ten previously valid cells. New three-run medians use `/tmp/rose-perf-reference/{primal_5959.primal-r,primal_5959.dual-r,quad_2717.dual-r,fx_quad_obj_5371.dual-r}.{1,2,3}.json`; summary `/tmp/rose-perf-reference/summary.json`. All original-library records identify Ryzen 9 5900X and library SHA256 `8e87fca221fdffc922331c8619426031b3e72dfd419446422f88af0cb0e98f9c`.
`-` refs are unreachable: original Gurobi rejects 80BAU3B, GREENBEA and mode-switch with license error 10010 (scratch `ID.error.json`, including direct 80BAU3B primal/dual attempts). The old 80BAU3B primal/dual refs used the patched library and have been removed; the old prose mode-switch time is unusable.

Regression checks: all 23 runs paired against the preceding solver on CPU 4. Quad-objective's apparent slowdown resolves in CPU-2 frozen-baseline/prior/current pairs: solve 0.510253/0.361690/0.262031 s. Quad_2717 has 0.8% more instructions (4610.5 -> 4649.0M); solve 0.798910 -> 0.856417 s remains within the old 0.058426 s spread. PDHG's timing shift resolves against the frozen baseline in nine interleaved individual solves on CPU 0: frozen/current solve medians 6.552595/6.582182 s, spreads 0.124273/0.117833 s, instructions 92391.1/92390.6M. The main table retains PDHG's full-harness median; the interleaved perf-stat invocations are a separate baseline audit. No unresolved frozen-baseline regression. Raw data `target/perf-qcache-full-{before,after}.txt`, `target/perf-qcache-recheck-{base,before,after}.txt`, `target/perf-qcache-pdhg-interleaved.{json,log}`. Host speed varies substantially across the sweep; use immediate pairs and spreads.
Acceptance: primal_5959 dual solve 0.888601 -> 0.648892 s (-27.0%; spreads 0.016741/0.019458 s), instructions 4978.1 -> 4012.6M (-19.4%). Primal solve 0.399437 -> 0.285415 s (-28.5%; spreads 0.016081/0.004421 s), instructions 3516.9 -> 2858.1M (-18.7%). CPU-2 quad-objective solve 0.361690 -> 0.262031 s (-27.6%), instructions 1994.8 -> 1619.8M (-18.8%). The two-entry cache retains only fresh factors of <50-row matrices without a forced pair gate or numeric-build diagnostics. Keys compare exact ordered indices, value bits and threshold; independent clones preserve factor state, pool-layout records and work estimates.

| run | ref solve | baseline wall | paired base wall | current wall | solve | current/ref | instr |
|---|---|---|---|---|---|---|---|
| degen2.primal-c | 0.012711 | 0.182 | - | 0.097 | 0.091757 | 7.22 | 481.4 |
| degen2.dual-r | 0.009653 | 0.125 | - | 0.077 | 0.072103 | 7.47 | 358.3 |
| degen3.primal-c | 0.163137 | 3.237 | - | 1.364 | 1.349310 | 8.27 | 7652.8 |
| degen3.primal-r | 0.193685 | 3.872 | - | 1.659 | 1.643226 | 8.48 | 9401.1 |
| perold.primal-c | 0.160994 | 2.430 | - | 1.331 | 1.325291 | 8.23 | 7534.9 |
| pilot4.primal-r | 0.015787 | 0.253 | - | 0.135 | 0.130124 | 8.24 | 820.5 |
| 80bau3b.primal-c | - | 1.618 | - | 1.481 | 1.459521 | - | 10089.5 |
| primal_5959.primal-r | 0.067884 | 0.563 | 0.938 | 0.288 | 0.285415 | 4.20 | 2858.1 |
| degen3.dual-c | 0.082843 | 1.002 | - | 0.480 | 0.470084 | 5.67 | 3758.2 |
| degen3.dual-r | 0.082947 | 1.008 | 1.007 | 0.481 | 0.471135 | 5.68 | 3769.1 |
| perold.dual-r | 0.038532 | 0.355 | 0.360 | 0.215 | 0.209904 | 5.45 | 1819.7 |
| 80bau3b.dual-r | - | 1.103 | 1.688 | 0.561 | 0.546418 | - | 5154.6 |
| greenbea.dual-c | - | 2.870 | 2.363 | 1.570 | 1.554566 | - | 10815.8 |
| greenbea.dual-r | - | 3.598 | 3.592 | 2.005 | 1.985986 | - | 10717.1 |
| mode_switch.dual-r | - | 25.815 | - | 7.344 | 7.197301 | - | 76028.3 |
| primal_5959.dual-r | 0.098192 | 1.300 | - | 0.652 | 0.648892 | 6.61 | 4012.6 |
| quad_2717.dual-r | 0.092814 | 1.130 | - | 0.860 | 0.856417 | 9.23 | 4649.0 |
| fx_quad_obj_5371.dual-r | 0.038077 | 0.334 | 0.514 | 0.265 | 0.262031 | 6.88 | 1619.8 |
| degen3.barrier-r | 0.075214 | 0.608 | 0.856 | 0.907 | 0.845537 | 11.24 | 5500.8 |
| 80bau3b.barrier-r | - | 0.379 | 0.549 | 0.567 | 0.531044 | - | 4058.5 |
| greenbea.barrier-r | - | 1.569 | 2.493 | 2.241 | 2.192543 | - | 16052.7 |
| 80bau3b.pdhg | - | 4.552 | 6.526 | 6.409 | 6.389652 | - | 92390.6 |
| greenbea.pdhg | - | 2.801 | - | 2.691 | 2.674784 | - | 37741.2 |

## Hotspots

Per profiled bench run: the commit it was taken at and the top self-time functions (percent), at most
five lines each. Re-profile a run when a landed change touched its hot path.

- primal_5959.dual-r, `tplnupps` (stale after factor reuse; frame-pointer build passes full parity; three solves capped at 16000 iterations, `target/perf-dual5959-full.data`): F80 add_signed 7.50%, primal control 5.70%, F80 multiply 4.54%, F80 nucleus 3.88%, double nucleus 3.74% self. All three reach numerical_failure, 15688 iterations, 7842 refactorizations. Refactorization is 51.10% inclusive, quad mirrors 21.69% before factor reuse. LUIN capture `target/perf-dual5959-luin.log` has 7822 records and 40 input-head/threshold combinations; two-entry LRU estimate 99.3% hits. Wide-shift arithmetic rewrites are rejected below. The earlier 11000-iteration capture stopped early and is not the full-run profile.
- pilot4.primal-r, `uprrtwlk` (stale after row-capacity count reuse; frame-pointer build passes full parity; twenty solves capped at 700 iterations, `target/perf-pilot4-fp.data`): update 10.42%, BTRAN 7.84%, fresh factor 7.05%, FTRAN 6.89%, primal control 6.83% self. Refactor 28.11%, fresh factor 10.90% inclusive. Every solve is optimal, 690 iterations, 16 refactorizations. Update annotation: numeric U-entry recount about 18% of update time, structural rcount scan also hot (the rowl shortcut was already rejected below); do not conflate the two counts.
- primal_5959.primal-r, `kyyvxyql` (stale after twin reuse; frame-pointer build passes full 987/987 parity; three bounded solves at 1999 Hz, `target/perf-fp-primal5959.data`, report with `--no-inline`): primal control 6.32%, F80 nucleus 5.46%, fresh F80 factor 5.09%, F80 add_signed 6.67% self. Refactor 53.53%, quad mirrors 28.88% inclusive before twin reuse. Log verifies numerical_failure, 10780 iterations, 5381 refactorizations; original reference status 12, 10780 iterations, refactor count unavailable. Default inline expansion duplicates frames; use PLAN.md's validated frame-pointer workflow for caller attribution.
- fx_quad_obj_5371.dual-r, `kyyvxyql` (stale after twin reuse; profiling build, three solves at 1999 Hz, `target/perf-members-fx.data`): F80 add_signed 8.06%, F80 nucleus 5.81%, multiply 5.74%, fresh F80 factor 5.17%, pack 3.83% self cycles. Solve remains infeasible, 5319 iterations, 5 refactorizations.
- degen3.dual-r, `yzqsoomv` (stale after LU membership changes; reprofile before choosing it; profiling build, three solves at 1999 Hz): iteration 10.59%, update 10.17%, L stage 8.00%, BTRAN 6.51%, ftran_with 6.51% self cycles. ftran_needed 5.64% plus term-buffer filling/column lookup 5.10%; the old per-row allocation is gone. Full solve: 2212 iterations, 24 refactorizations.
- mode_switch.dual-r, `xzuktkmx` (stale after F80/LU changes; reprofile before choosing it): update 22.21%, LU nucleus 15.90%, BTRAN 13.33%, fresh-factor row-list filtering 7.67%, factor_structure 3.36% self cycles. Full solve: 11393 iterations, 17 refactorizations.

## In progress

(none)

## Candidates

Ordered by current/ref; reprofile any changed path before choosing a change. Ratios use the current
sweep's solve medians; paired old/new times determine whether a change gains or regresses.

- degen3.barrier-r (11.24x), quad_2717.dual-r (9.23x), degen3.primal-r (8.48x), degen3.primal-c (8.27x): reprofile in this order. Profiles of changed paths are stale. Quad_2717's +0.8% instruction cost is a lead for avoiding fresh-cache copies on non-repeated inputs; its timing remains within spread.
- Other reachable runs follow their table ratios. DEGEN3 row-dot filling/lookup remains an untried lead after buffer reuse; preserve accumulation order and work estimates. Quad twins share construction and compact repeated factors are cached; row-capacity growth already uses the numeric U transpose.
- Runs without reachable references, after the reachable ratios: mode-switch update (22.21%) still scans A/B lists for other memberships; BTRAN (13.33%) allocates scratch/recomputes L levels; LU nucleus (15.90%); fresh-factor row-list filtering (7.67%). These old profiles are stale; annotate again before choosing a change.
- all: dump::on/record construction without dumps and allocation in FTRAN/BTRAN, pricing and ratio tests remain unmeasured leads.

## Tried, no gain

- F80 add aligned with three guard bits, plus direct normal-result rounding: frozen old-add/pack differential test and fast parity 940/940 pass, instructions fall 4977.9 -> 4836.8M (-2.8%), but primal_5959 dual CPU-4 solve 0.555599 -> 0.557742 s (spreads 0.006716/0.011524 s). Alignment alone also has no timing gain on CPUs 2/4 despite -1.7% instructions. Reverted both (`target/perf-add-{target,recheck}.log`, `target/perf-add-direct-{target,recheck}.log`).
- Read eliminate work counts directly from maintained rowl lengths instead of rebuilding counts: parity full 987/987 and release tests pass, but mode_switch.dual-r 22.849 -> 21.851 s is smaller than the old 1.849 s spread (instructions 146886.1 -> 140945.4M, -4.0%). Alternating DEGEN2 primal 0.092/0.090 then 0.091/0.091 s; dual 0.068/0.067 s, marginal against spread. Reverted.

## Tried, breaks parity

(idea, first differing record from `rose-parity check`, why)
