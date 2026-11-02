# Performance state

Current state of the work in [PLAN.md](PLAN.md). Rewrite in place; no log (history is `jj log`).

## Numbers

Baseline solver: `8b645007` (same solver as `vlsokskq`). Current solver: `vnqmrmnq` (immutable coefficients stored alongside maintained dual row entries).
Parity: full tier 987/987; release tests pass. Three timed runs after one warm-up, hyperfine 2.0.0; instruction counts are millions from perf stat.

`baseline wall` is the fixed initial process-wall median. `paired base wall` is the fresh baseline-binary median from alternating baseline/last-landed/candidate checks on one CPU (CPU 0 unless stated otherwise); it resolves host drift without changing the fixed baseline. `current wall` and `solve` use the new harness; **current/ref uses solve-only time**, `diagnostics.seconds`, against original Gurobi Runtime.

Reference provenance: original-library cached records `../reverse-engineer/runs/{degen2,degen3,perold,pilot4}-{profile}.json` for the ten previously valid cells. New three-run medians use `/tmp/rose-perf-reference/{primal_5959.primal-r,primal_5959.dual-r,quad_2717.dual-r,fx_quad_obj_5371.dual-r}.{1,2,3}.json`; summary `/tmp/rose-perf-reference/summary.json`. All original-library records identify Ryzen 9 5900X and library SHA256 `8e87fca221fdffc922331c8619426031b3e72dfd419446422f88af0cb0e98f9c`.
`-` refs are unreachable: original Gurobi rejects 80BAU3B, GREENBEA and mode-switch with license error 10010 (scratch `ID.error.json`, including direct 80BAU3B primal/dual attempts). The old 80BAU3B primal/dual refs used the patched library and have been removed; the old prose mode-switch time is unusable.

Regression checks: all 23 runs paired against the preceding solver on CPU 4. Two apparent increases resolve in three alternating CPU-2 pairs: PEROLD primal-controlled solve 1.243203/1.213620 s (old/new spreads 0.209717/0.438725 s), with matching instructions; DEGEN3 dual-controlled 0.420798/0.398702 s (spreads 0.003287/0.002310 s), instructions 3348.6/3240.8M (-3.2%). No new fixed-baseline regression: PDHG is below the fixed baseline; 80BAU3B barrier remains below its fresh paired baseline wall. Raw data `target/perf-dualrows-full-{before,after}.txt`, `target/perf-dualrows-recheck-{before,after}.txt`. Host speed varies substantially; use immediate pairs and spreads.
Acceptance: CPU-2 PEROLD dual-release solve 0.322409 -> 0.312073 s (-3.2%; spreads 0.005671/0.007907 s), instructions 1768.7 -> 1717.9M (-2.9%) (`target/perf-dualrows-target.log`). The full sweep confirms 0.320436 -> 0.303641 s (-5.2%), beyond either spread. DEGEN3 dual-release improves 0.664241 -> 0.624230 s with -3.2% instructions in the target check. Maintained nonbasic row entries now carry the immutable coefficient previously found by searching their column on every product. Building, replacing, deleting, appending and importing rows preserve index order and first-match coefficient bits (including the missing-entry zero in imported rows). Arithmetic operations, accumulation order and WORK estimates remain unchanged.

| run | ref solve | baseline wall | paired base wall | current wall | solve | current/ref | instr |
|---|---|---|---|---|---|---|---|
| degen2.primal-c | 0.012711 | 0.182 | - | 0.089 | 0.084641 | 6.66 | 439.2 |
| degen2.dual-r | 0.009653 | 0.125 | - | 0.074 | 0.068767 | 7.12 | 329.3 |
| degen3.primal-c | 0.163137 | 3.237 | - | 1.222 | 1.206435 | 7.40 | 6569.1 |
| degen3.primal-r | 0.193685 | 3.872 | - | 0.926 | 0.915395 | 4.73 | 7824.1 |
| perold.primal-c | 0.160994 | 2.430 | - | 1.220 | 1.213620 | 7.54 | 6926.3 |
| pilot4.primal-r | 0.015787 | 0.253 | - | 0.087 | 0.083004 | 5.26 | 798.2 |
| 80bau3b.primal-c | - | 1.618 | - | 0.843 | 0.827655 | - | 8725.1 |
| primal_5959.primal-r | 0.067884 | 0.563 | 0.938 | 0.279 | 0.276695 | 4.08 | 2824.5 |
| degen3.dual-c | 0.082843 | 1.002 | - | 0.409 | 0.398702 | 4.81 | 3240.8 |
| degen3.dual-r | 0.082947 | 1.008 | 1.007 | 0.617 | 0.603048 | 7.27 | 3250.3 |
| perold.dual-r | 0.038532 | 0.355 | 0.360 | 0.310 | 0.303641 | 7.88 | 1718.1 |
| 80bau3b.dual-r | - | 1.103 | 1.688 | 0.736 | 0.713060 | - | 4260.4 |
| greenbea.dual-c | - | 2.870 | 2.363 | 1.718 | 1.698250 | - | 9197.0 |
| greenbea.dual-r | - | 3.598 | 3.592 | 1.698 | 1.678163 | - | 9099.7 |
| mode_switch.dual-r | - | 25.815 | - | 4.480 | 4.444583 | - | 50572.3 |
| primal_5959.dual-r | 0.098192 | 1.300 | - | 0.613 | 0.610127 | 6.21 | 3959.8 |
| quad_2717.dual-r | 0.092814 | 1.130 | - | 0.691 | 0.687945 | 7.41 | 4084.8 |
| fx_quad_obj_5371.dual-r | 0.038077 | 0.334 | 0.514 | 0.256 | 0.253052 | 6.65 | 1592.2 |
| degen3.barrier-r | 0.075214 | 0.608 | 0.856 | 0.471 | 0.408261 | 5.43 | 3779.9 |
| 80bau3b.barrier-r | - | 0.379 | 0.546 | 0.415 | 0.379118 | - | 3427.8 |
| greenbea.barrier-r | - | 1.569 | 2.493 | 1.501 | 1.450647 | - | 13029.7 |
| 80bau3b.pdhg | - | 4.552 | 6.526 | 4.475 | 4.461152 | - | 92388.9 |
| greenbea.pdhg | - | 2.801 | - | 2.739 | 2.721552 | - | 37741.1 |

## Hotspots

Per profiled bench run: the commit it was taken at and the top self-time functions (percent), at most
five lines each. Re-profile a run when a landed change touched its hot path.

- perold.dual-r, `szutlnuy` (stale after coefficient caching; frame-pointer build passes full 987/987 parity; ten solves capped at 1300 iterations, `target/perf-perold-dual.data`): iteration 13.15%, FTRAN 8.05%, update 7.93%, ftran_needed 6.03%, L stage 5.51% self. Iteration 54.50%, refactor 20.05%, FTRAN 16.42%, ftran_needed 16.31% inclusive. All ten optimal in 1273 iterations, 27 refactorizations. The annotated column searches (`target/perf-perold-dual-iteration-annotate.txt`) are now removed from the row product; reprofile before choosing further work. The work estimator also scans touched rows for a distinct-column count used only by sparse estimates.
- perold.primal-c, `yryslslw` (stale after dense work-count shortcut; frame-pointer build passes full 987/987 parity; three solves capped at 3500 iterations, `target/perf-perold.data`): BTRAN 8.71%, primal control 7.01%, row_product_work 6.39%, update 6.12%, fresh factor 6.10% self. Refactor 26.63%, FTRAN 11.82%, BTRAN 11.03% inclusive. All three optimal in 3399 iterations, 90 refactorizations; original reference also has 3399 iterations. The work-estimator annotation (`target/perf-perold-rowprod-annotate.txt`) attributes substantial time to checks in the touched-row scan; dense counts now avoid that scan. Reprofile before choosing the remaining sparse work.
- degen3.primal-r, `ltkyqoqr` (stale after member-based row ordering; frame-pointer build passes full 987/987 parity; three solves capped at 5300 iterations, `target/perf-primal3-next.data`): update 10.65%, FTRAN 9.69%, BTRAN 7.93%, fresh factor 4.59%, fresh row-list filtering 4.32% self. Signature classes remain 2.99%. All three optimal in 5165 iterations, 66 refactorizations. The full basis-order scan for fresh singleton-row lists is now replaced by sorting the row members; reprofile before choosing more work. Update annotation (`target/perf-primal3-next-update-annotate.txt`) still highlights eliminate row-count reconstruction, whose rowl-length replacement was already rejected below.
- quad_2717.dual-r, `zzvmznmp` (stale after cache sizing; frame-pointer build passes full parity; three solves capped at 20000 iterations, `target/perf-quad2717.data`): F80 add 13.58%, primal control 8.94%, F80 multiply 8.17%, pack 3.80%, rank_products_q 3.61% self. Refactorization 25.94%, quad mirrors 11.50%, rank_products_q 11.26%, cost_change_work 11.18% inclusive. All solves optimal in 16010 iterations with 3022 refactorizations. `target/perf-quad2717-luin.log`: 3015 LUIN records, 14 distinct input-head/threshold keys; FIFO hits 110/2998/3000 at capacities 2/4/8. The cache now holds four entries; larger capacities add only marginal hits. Cost-change work also updates reduced costs and cannot simply be removed.
- degen3.barrier-r, `munmovyt` (stale after kernel dispatch; frame-pointer build passes full parity; three solves capped at 100 iterations, `target/perf-barrier-degen3-valid.data`): FMA dispatch 26.46%, DGEMM 24.62%, existing hardware FMA implementation 14.77%, numeric_factor 7.40%, DTRSM 4.36% self. DGEMM 49.04%, numeric factor 81.45% inclusive; all three solves optimal in 14 iterations. The per-product FMA wrapper is an indirect call plus frame/return overhead (annotated `target/perf-barrier-fma-annotate.txt`); per-call dispatch now removes that per-product overhead without adding fused operations. BORE3D is a cheap witness (16 barrier iterations, 12 ms).
- primal_5959.dual-r, `tplnupps` (stale after factor reuse; frame-pointer build passes full parity; three solves capped at 16000 iterations, `target/perf-dual5959-full.data`): F80 add_signed 7.50%, primal control 5.70%, F80 multiply 4.54%, F80 nucleus 3.88%, double nucleus 3.74% self. All three reach numerical_failure, 15688 iterations, 7842 refactorizations. Refactorization is 51.10% inclusive, quad mirrors 21.69% before factor reuse. LUIN capture `target/perf-dual5959-luin.log` has 7822 records and 40 input-head/threshold combinations; two-entry LRU estimate 99.3% hits. Wide-shift arithmetic rewrites are rejected below. The earlier 11000-iteration capture stopped early and is not the full-run profile.
- pilot4.primal-r, `uprrtwlk` (stale after row-capacity count reuse; frame-pointer build passes full parity; twenty solves capped at 700 iterations, `target/perf-pilot4-fp.data`): update 10.42%, BTRAN 7.84%, fresh factor 7.05%, FTRAN 6.89%, primal control 6.83% self. Refactor 28.11%, fresh factor 10.90% inclusive. Every solve is optimal, 690 iterations, 16 refactorizations. Update annotation: numeric U-entry recount about 18% of update time, structural rcount scan also hot (the rowl shortcut was already rejected below); do not conflate the two counts.
- primal_5959.primal-r, `kyyvxyql` (stale after twin reuse; frame-pointer build passes full 987/987 parity; three bounded solves at 1999 Hz, `target/perf-fp-primal5959.data`, report with `--no-inline`): primal control 6.32%, F80 nucleus 5.46%, fresh F80 factor 5.09%, F80 add_signed 6.67% self. Refactor 53.53%, quad mirrors 28.88% inclusive before twin reuse. Log verifies numerical_failure, 10780 iterations, 5381 refactorizations; original reference status 12, 10780 iterations, refactor count unavailable. Default inline expansion duplicates frames; use PLAN.md's validated frame-pointer workflow for caller attribution.
- fx_quad_obj_5371.dual-r, `kyyvxyql` (stale after twin reuse; profiling build, three solves at 1999 Hz, `target/perf-members-fx.data`): F80 add_signed 8.06%, F80 nucleus 5.81%, multiply 5.74%, fresh F80 factor 5.17%, pack 3.83% self cycles. Solve remains infeasible, 5319 iterations, 5 refactorizations.
- degen3.dual-r, `sormxnyt` (stale after metadata caching; frame-pointer build passes full parity; three solves capped at 2300 iterations, `target/perf-btplan-degen3.data`): iteration 13.36%, update 8.89%, L stage 7.68%, BTRAN 7.21%, FTRAN 6.93% self; ftran_needed 5.34%, with the old column-scan filler gone. BTRAN annotation (`target/perf-btplan-btran-annotate.txt`) attributes over a quarter of its self samples to rebuilding transpose levels, and more to rebuilding the inverse row map. L columns, initial row permutation and loff are immutable between refactorizations. All solves optimal in 2212 iterations, 24 refactorizations; original reference also has 2212 iterations.
- mode_switch.dual-r, `xzuktkmx` (stale after F80/LU changes; reprofile before choosing it): update 22.21%, LU nucleus 15.90%, BTRAN 13.33%, fresh-factor row-list filtering 7.67%, factor_structure 3.36% self cycles. Full solve: 11393 iterations, 17 refactorizations.

## In progress

(none)

## Candidates

Ordered by current/ref; reprofile any changed path before choosing a change. Ratios use the current
sweep's solve medians; paired old/new times determine whether a change gains or regresses.

- perold.dual-r (7.88x), perold.primal-c (7.54x), quad_2717.dual-r (7.41x), degen3.primal-c (7.40x): reprofile in this order. Profiles of changed paths are stale.
- DEGEN3 dual: tau row-value lookup and immutable BTRAN metadata are settled. Reprofile before inspecting remaining iteration control, updates and scratch allocation.
- Other reachable runs follow their table ratios. Quad twins share construction and four compact repeated factors are cached; row-capacity growth already uses the numeric U transpose.
- Quad_2717: after profiling again, inspect F80 arithmetic and the ranked row products. Cost-change work also performs required reduced-cost updates; some extended row-product outputs may be unused, but masks must preserve every consumed value and all work bookings. The too-small factor-cache working set is settled.
- Barrier: profile again before choosing more work. DGEMM now dispatches existing fused operations once per call; other BLAS-like kernels still use their original scalar dispatch and accumulation order. A smaller remaining DGEMM gap may involve independent outputs and data reuse, which must preserve every per-output reduction order.
- Runs without reachable references, after the reachable ratios: mode-switch update (22.21%) still scans A/B lists for other memberships; BTRAN (13.33% before metadata caching) allocates scratch; LU nucleus (15.90%); fresh-factor row-list filtering (7.67%). These old profiles are stale; annotate again before choosing a change.
- all: dump::on/record construction without dumps and allocation in FTRAN/BTRAN, pricing and ratio tests remain unmeasured leads.

## Tried, no gain

- Precompute eligibility for row_product_work scans longer than n: fast parity 940/940 passes, PEROLD primal instructions -0.4%, solve 1.342457 -> 1.332925 s is within 0.249063/0.547923 s spreads; DEGEN3 primal-controlled has no timing gain and DEGEN2 instructions rise slightly (`target/perf-rowelig-long-target.log`). Reverted.

- Lazy per-column eligibility memoization in row_product_work: fast parity 940/940 passes, PEROLD primal instructions -1.1%, but CPU-2 solve 1.322285 -> 1.320622 s is within 0.013795 s spread; CPU-0 alternating recheck crosses speed bands with no clear aggregate gain (`target/perf-rowelig-{target,recheck}.log`). Reverted.

- Binary lookup in the U transpose for tau terms: fast parity 940/940 passes; DEGEN3 dual instructions 3769.4 -> 3753.4M (-0.4%), solve 0.739033 -> 0.716611 s is within old 0.025890 s spread. DEGEN2/PEROLD instructions rise slightly. Replaced by direct row-map trial (`target/perf-taulookup-target.log`).
- F80 add aligned with three guard bits, plus direct normal-result rounding: frozen old-add/pack differential test and fast parity 940/940 pass, instructions fall 4977.9 -> 4836.8M (-2.8%), but primal_5959 dual CPU-4 solve 0.555599 -> 0.557742 s (spreads 0.006716/0.011524 s). Alignment alone also has no timing gain on CPUs 2/4 despite -1.7% instructions. Reverted both (`target/perf-add-{target,recheck}.log`, `target/perf-add-direct-{target,recheck}.log`).
- Read eliminate work counts directly from maintained rowl lengths instead of rebuilding counts: parity full 987/987 and release tests pass, but mode_switch.dual-r 22.849 -> 21.851 s is smaller than the old 1.849 s spread (instructions 146886.1 -> 140945.4M, -4.0%). Alternating DEGEN2 primal 0.092/0.090 then 0.091/0.091 s; dual 0.068/0.067 s, marginal against spread. Reverted.

## Tried, breaks parity

(idea, first differing record from `rose-parity check`, why)
