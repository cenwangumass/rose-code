# Performance state

Current state of the work in [PLAN.md](PLAN.md). Rewrite in place; no log (history is `jj log`).

## Numbers

Baseline solver: `8b645007` (same solver as `vlsokskq`). Current solver: `rtnmszvk` (scan joint-dense tau tail slots without symbolic reachability).
Parity: full tier 987/987; release tests pass. Three timed runs after one warm-up, hyperfine 2.0.0; instruction counts are millions from perf stat.

`baseline wall` is the fixed initial process-wall median. `paired base wall` is the fresh baseline-binary median from alternating baseline/last-landed/candidate checks on one CPU (CPU 0 unless stated otherwise); it resolves host drift without changing the fixed baseline. `current wall` and `solve` use the new harness; **current/ref uses solve-only time**, `diagnostics.seconds`, against original Gurobi Runtime.

Reference provenance: original-library cached records `../reverse-engineer/runs/{degen2,degen3,perold,pilot4}-{profile}.json` for the ten previously valid cells. New three-run medians use `/tmp/rose-perf-reference/{primal_5959.primal-r,primal_5959.dual-r,quad_2717.dual-r,fx_quad_obj_5371.dual-r}.{1,2,3}.json`; summary `/tmp/rose-perf-reference/summary.json`. All original-library records identify Ryzen 9 5900X and library SHA256 `8e87fca221fdffc922331c8619426031b3e72dfd419446422f88af0cb0e98f9c`.
`-` refs are unreachable: original Gurobi rejects 80BAU3B, GREENBEA and mode-switch with license error 10010 (scratch `ID.error.json`, including direct 80BAU3B primal/dual attempts). The old 80BAU3B primal/dual refs used the patched library and have been removed; the old prose mode-switch time is unusable.

Regression checks: all 23 runs paired against the preceding solver on CPU 4, with pair order alternating by run. PEROLD dual's timing increase resolves in three alternating CPU-2 follow-up pairs: 0.241011 -> 0.230672 s (spreads 0.001261/0.004904 s), instructions -3.3%. No unresolved paired or fixed-baseline regressions. Raw data `target/perf-needtail-full-{before,after}.txt` and `target/perf-needtail-followup-{before,after}.txt`. Host speed varies substantially; use immediate pairs and spreads.
Acceptance: the full CPU-4 sweep gives DEGEN3 dual-controlled solve 0.321696 -> 0.313540 s (-2.5%; spreads 0.003753/0.002036 s), instructions 2636.4 -> 2574.2M (-2.4%); dual-release 0.321662 -> 0.313025 s (-2.7%; spreads 0.000790/0.001757 s), instructions -2.4%. DEGEN2 dual improves 0.059400 -> 0.056212 s (-5.4%; spreads 0.002064/0.001961 s), instructions -2.1%; the earlier target pair confirms 0.037215 -> 0.036790 s (-1.1%; spreads 0.000179/0.000156 s, `target/perf-needtail-target-{before,after}.txt`). Alternating PEROLD pairs improve 4.3% beyond both spreads. Joint-dense tau solves scan all tail slots in descending order; extra zero columns are zeroed by the unchanged kernel, preserving values and raw flags. Other tau paths retain their original graph walk. Frozen old/new comparisons cover outputs, tau state and WORK estimates on random/extreme f64/F80 and captured SC50A inputs; full parity 987/987 and release tests pass.

| run | ref solve | baseline wall | paired base wall | current wall | solve | current/ref | instr |
|---|---|---|---|---|---|---|---|
| degen2.primal-c | 0.012711 | 0.182 | - | 0.072 | 0.067726 | 5.33 | 355.7 |
| degen2.dual-r | 0.009653 | 0.125 | - | 0.061 | 0.056212 | 5.82 | 282.9 |
| degen3.primal-c | 0.163137 | 3.237 | - | 0.903 | 0.887450 | 5.44 | 5004.3 |
| degen3.primal-r | 0.193685 | 3.872 | - | 0.689 | 0.679683 | 3.51 | 5704.1 |
| perold.primal-c | 0.160994 | 2.430 | - | 0.615 | 0.610233 | 3.79 | 5525.1 |
| pilot4.primal-r | 0.015787 | 0.253 | - | 0.066 | 0.062193 | 3.94 | 622.8 |
| 80bau3b.primal-c | - | 1.618 | - | 0.783 | 0.768171 | - | 8064.7 |
| primal_5959.primal-r | 0.067884 | 0.563 | 0.938 | 0.223 | 0.221245 | 3.26 | 2280.4 |
| degen3.dual-c | 0.082843 | 1.002 | - | 0.323 | 0.313540 | 3.78 | 2574.2 |
| degen3.dual-r | 0.082947 | 1.008 | 1.007 | 0.323 | 0.313025 | 3.77 | 2583.8 |
| perold.dual-r | 0.038532 | 0.355 | 0.360 | 0.237 | 0.230672 | 5.99 | 1358.1 |
| 80bau3b.dual-r | - | 1.103 | 1.688 | 0.602 | 0.579686 | - | 3823.0 |
| greenbea.dual-c | - | 2.870 | 2.363 | 1.365 | 1.345210 | - | 7180.8 |
| greenbea.dual-r | - | 3.598 | 3.592 | 1.343 | 1.322534 | - | 7137.0 |
| mode_switch.dual-r | - | 25.815 | - | 4.001 | 3.948089 | - | 25106.6 |
| primal_5959.dual-r | 0.098192 | 1.300 | - | 0.308 | 0.306076 | 3.12 | 3167.1 |
| quad_2717.dual-r | 0.092814 | 1.130 | - | 0.359 | 0.356445 | 3.84 | 3176.4 |
| fx_quad_obj_5371.dual-r | 0.038077 | 0.334 | 0.514 | 0.118 | 0.116070 | 3.05 | 1214.7 |
| degen3.barrier-r | 0.075214 | 0.608 | 0.856 | 0.317 | 0.276323 | 3.67 | 3778.7 |
| 80bau3b.barrier-r | - | 0.379 | 0.546 | 0.276 | 0.251960 | - | 3425.8 |
| greenbea.barrier-r | - | 1.569 | 2.493 | 1.009 | 0.975538 | - | 13028.9 |
| 80bau3b.pdhg | - | 4.552 | 6.526 | 6.276 | 6.256346 | - | 92390.0 |
| greenbea.pdhg | - | 2.801 | - | 1.894 | 1.882034 | - | 37740.7 |

## Hotspots

Per profiled bench run: the commit it was taken at and the top self-time functions (percent), at most
five lines each. Re-profile a run when a landed change touched its hot path.

- degen2.primal-c, `rtyoxnwx` (stale after phase-one class reuse; frame-pointer build passes full 987/987 parity; forty solves capped at 2000 iterations, `target/perf-degen2c-next.data`): BTRAN 9.80%, primal control 7.95%, FTRAN 7.22%, update 4.98%, LU nucleus 4.26% self. Refactorization is 24.70%, FTRAN 13.22%, BTRAN 10.22% inclusive. All forty optimal in 531 iterations, 16 refactorizations. Signature construction is 1.26% and sig_class 1.11% self; repeated phase-one classification within the stable update window is now shared. Bounds-check removal and direct growth-diagonal access were measured and rejected below. Reprofile before choosing remaining BTRAN or refactor work.
- perold.dual-r, `tvqwvukl` (stale after dense FTRAN tail scan; frame-pointer build passes full 987/987 parity; ten solves capped at 1300 iterations, `target/perf-perold-next3.data`): FTRAN 9.67%, dual iteration 8.44%, LU nucleus 7.30%, ftran_needed 7.16%, L stage 6.70% self. FTRAN is 16.74%, ftran_needed 15.02%, refactorization 24.84% inclusive. All ten optimal in 1273 iterations, 27 refactorizations. FTRAN annotation (`target/perf-perold-next3-ftran-annotate.txt`) puts about 20% of its self samples in symbolic tail reachability and membership checks; dense FTRAN now skips that walk. Reprofile before choosing remaining FTRAN or refactorization work.
- perold.primal-c, `yryslslw` (stale after growth-factor reuse and the dense work-count shortcut; frame-pointer build passes full 987/987 parity; three solves capped at 3500 iterations, `target/perf-perold.data`): BTRAN 8.71%, primal control 7.01%, row_product_work 6.39%, update 6.12%, fresh factor 6.10% self. Refactor 26.63%, FTRAN 11.82%, BTRAN 11.03% inclusive. All three optimal in 3399 iterations, 90 refactorizations; original reference also has 3399 iterations. The work-estimator annotation (`target/perf-perold-rowprod-annotate.txt`) attributes substantial time to checks in the touched-row scan; dense counts now avoid that scan. Reprofile before choosing the remaining sparse work.
- degen3.primal-c, `oytzqrpy` (stale after growth-factor reuse; frame-pointer build passes full 987/987 parity; three solves capped at 10000 iterations, `target/perf-degen3c-next.data`): primal control 8.91%, BTRAN 8.09%, FTRAN 6.51%, LU nucleus 5.82%, fresh factor 5.25% self. Refactorization is 31.03%, reference_order_for 19.77%, BTRAN 11.48%, FTRAN 10.88% inclusive; dense BTRAN output compaction is another 2.58% self. All three optimal in 2492 iterations, 49 refactorizations. The accepted growth-check factor now carries through to installation; reprofile before choosing remaining reference-order or BTRAN work.
- degen3.primal-r, `xutnmxks` (stale after dense BTRAN output reservation; frame-pointer build passes full 987/987 parity; three solves capped at 5300 iterations, `target/perf-degen3r-rest.data`): FTRAN 12.37%, BTRAN 7.48%, update 7.08%, LU nucleus 5.90%, primal control 4.73% self. Refactorization is 24.83%, FTRAN 21.75%, BTRAN 11.26% inclusive; dense BTRAN output compaction is 3.45% inclusive/2.88% self. All three optimal in 5165 iterations, 66 refactorizations. FTRAN annotation (`target/perf-degen3r-rest-ftran-annotate.txt`) highlights tail graph traversal, scatter bounds checks and output done/fallback scans; prior dense-FTRAN output simplification was rejected below. Dense BTRAN output now reserves its maximum row count before scanning; reprofile before choosing remaining FTRAN graph, scatter or output work.
- quad_2717.dual-r, `ussxlpst` (stale after native ranked row products; frame-pointer build passes full 987/987 parity; eight solves capped at 17000 iterations, `target/perf-quad-next.data`): F80 add 15.13%, primal control 10.38%, F80 multiply 8.89%, pack 5.06%, rank_products_q 4.22% self. Row products are 13.20%, cost_change_work 12.77%, refactorization 17.26% inclusive. All eight optimal in 16010 iterations, 3022 refactorizations. Ranked row products now use native extended arithmetic under the required control mode; the software kernel remains the fallback. Reprofile before choosing remaining arithmetic or set_d_from_y work. The rejected masks, coefficient cache, distant-exponent shortcut and forced add inlining are listed below.
- degen3.barrier-r, `munmovyt` (stale after kernel dispatch; frame-pointer build passes full parity; three solves capped at 100 iterations, `target/perf-barrier-degen3-valid.data`): FMA dispatch 26.46%, DGEMM 24.62%, existing hardware FMA implementation 14.77%, numeric_factor 7.40%, DTRSM 4.36% self. DGEMM 49.04%, numeric factor 81.45% inclusive; all three solves optimal in 14 iterations. The per-product FMA wrapper is an indirect call plus frame/return overhead (annotated `target/perf-barrier-fma-annotate.txt`); per-call dispatch now removes that per-product overhead without adding fused operations. BORE3D is a cheap witness (16 barrier iterations, 12 ms).
- primal_5959.dual-r, `tplnupps` (stale after factor reuse; frame-pointer build passes full parity; three solves capped at 16000 iterations, `target/perf-dual5959-full.data`): F80 add_signed 7.50%, primal control 5.70%, F80 multiply 4.54%, F80 nucleus 3.88%, double nucleus 3.74% self. All three reach numerical_failure, 15688 iterations, 7842 refactorizations. Refactorization is 51.10% inclusive, quad mirrors 21.69% before factor reuse. LUIN capture `target/perf-dual5959-luin.log` has 7822 records and 40 input-head/threshold combinations; two-entry LRU estimate 99.3% hits. Wide-shift arithmetic rewrites are rejected below. The earlier 11000-iteration capture stopped early and is not the full-run profile.
- pilot4.primal-r, `lyvknlkq` (stale after LU bucket bounds; frame-pointer build passes full 987/987 parity; twenty solves capped at 700 iterations, `target/perf-pilot-next3.data`): primal control 8.90%, BTRAN 7.91%, update 7.32%, FTRAN 6.21%, LU nucleus 5.49% self. Refactorization is 29.17%, FTRAN 11.51%, update 10.62%, BTRAN 8.18% inclusive. All twenty optimal in 690 iterations, 16 refactorizations. LU annotation (`target/perf-pilot-next3-nucleus.txt`) puts about 20% of nucleus self samples in reverse scans over empty bucket suffixes; those scans now use occupied-prefix upper bounds. The two update membership shortcuts failed timing below. Reprofile before choosing remaining BTRAN, FTRAN, update or refactorization work.
- primal_5959.primal-r, `ozmqzstm` (stale after compact double-structure caching; frame-pointer build passes full 987/987 parity; five solves capped at 11000 iterations, `target/perf-primal5959-next.data`): primal control 9.81%, F80 add 8.85%, F80 multiply 4.39%, double LU nucleus 4.31%, set_d_from_y 3.92% self; malloc/free total 7.16%. Refactorization is 38.43%, reference_order_for 13.60%, double factor_structure 12.13% inclusive. All five numerical_failure in 10780 iterations, 5381 refactorizations, matching the original status and iteration count. Both extended factors and double structures are now cached for compact repeated bases. Reprofile before choosing remaining refactorization, F80 arithmetic or set_d_from_y work.
- fx_quad_obj_5371.dual-r, `kyyvxyql` (stale after twin reuse; profiling build, three solves at 1999 Hz, `target/perf-members-fx.data`): F80 add_signed 8.06%, F80 nucleus 5.81%, multiply 5.74%, fresh F80 factor 5.17%, pack 3.83% self cycles. Solve remains infeasible, 5319 iterations, 5 refactorizations.
- degen3.dual-c, `pkxouoxq` (stale after joint-dense tau tail scan; frame-pointer build passes full 987/987 parity; five solves capped at 2300 iterations, `target/perf-dual3c-next.data`): L stage 10.57%, dual iteration 7.73%, ftran_needed 7.66%, FTRAN 7.32%, BTRAN 5.98% self. Refactorization is 22.12%, ftran_needed 21.32%, L stage 14.88%, FTRAN 13.87% inclusive. All five solves are optimal in 2212 iterations, 24 refactorizations. Dense L output conversion is another 2.74% self. The L-stage annotation (`target/perf-dual3c-next-lstage.txt`) shows scatter coefficient indexing and repeated vector metadata loads; the zipped-row/value trial failed timing below. The needed-vector annotation (`target/perf-dual3c-next-needed.txt`) highlights symbolic tail traversal, now skipped for the joint-dense path, and the remaining tail membership test in numeric scatter. Reprofile before choosing more L-stage, scatter or iteration work.
- degen3.dual-r, `rxqrkvwm` (stale after eta gather-buffer reuse; frame-pointer build passes full 987/987 parity; five solves capped at 2300 iterations, `target/perf-dual3-rest.data`): L stage 9.24%, FTRAN 8.62%, ftran_needed 7.61%, iteration 7.02%, BTRAN 5.83% self. ftran_needed is 21.69%, L stage 13.85%, FTRAN 15.30%, refactorization 23.59% inclusive. All five optimal in 2212 iterations, 24 refactorizations. Dense L-stage output compaction is another 2.52% self; reserving it failed the timing gate below. Eta gathers now reuse one buffer per call. Reprofile before choosing more L-stage numeric scatter or output conversion work (`target/perf-dual3-rest-lstage-annotate.txt`).
- mode_switch.dual-r, `xzuktkmx` (stale after F80/LU changes; reprofile before choosing it): update 22.21%, LU nucleus 15.90%, BTRAN 13.33%, fresh-factor row-list filtering 7.67%, factor_structure 3.36% self cycles. Full solve: 11393 iterations, 17 refactorizations.

## In progress

None.

## Candidates

Ordered by current/ref; reprofile any changed path before choosing a change. Ratios use the current
sweep's solve medians; paired old/new times determine whether a change gains or regresses.

- perold.dual-r (5.99x), degen2.dual-r (5.82x), degen3.primal-c (5.44x), degen2.primal-c (5.33x), pilot4.primal-r (3.94x): work in this order. Reprofile changed paths before choosing a change; dual profiles are stale after the joint-dense tau shortcut. Host speed bands can reorder these snapshot ratios.
- DEGEN3 dual: tau row-value lookup and immutable BTRAN metadata are settled. Reprofile before inspecting remaining iteration control, updates and scratch allocation.
- Other reachable runs follow their table ratios. Quad twins share construction and four compact repeated factors and double structures are cached; row-capacity growth already uses the numeric U transpose.
- Quad_2717: ranked extended row products now use native x87 arithmetic with a software fallback. Reprofile before inspecting remaining F80 arithmetic and set_d_from_y; the old call shares no longer describe the current solver. Cost-change work performs required reduced-cost updates; the factor-cache working set is settled.
- Barrier: profile again before choosing more work. DGEMM now dispatches existing fused operations once per call; other BLAS-like kernels still use their original scalar dispatch and accumulation order. A smaller remaining DGEMM gap may involve independent outputs and data reuse, which must preserve every per-output reduction order.
- Runs without reachable references, after the reachable ratios: mode-switch update (22.21%) still scans A/B lists for other memberships; BTRAN (13.33% before metadata caching) allocates scratch; LU nucleus (15.90%); fresh-factor row-list filtering (7.67%). These old profiles are stale; annotate again before choosing a change.
- all: dump::on/record construction without dumps and allocation in FTRAN/BTRAN, pricing and ratio tests remain unmeasured leads.

## Tried, no gain

- Hold immutable kind/diag/U slices across the needed-vector scatter: fast parity 940/940 passes and DEGEN3 dual-controlled instructions -0.8%, but solve 0.461177 -> 0.488074 s has no gain (spreads 0.041088/0.008475 s); DEGEN2/PEROLD stay within spread. Reverted (`target/perf-tauviews-target-{before,after}.txt`).

- Zip L-stage row/value slices to remove coefficient index checks: fast parity 940/940 passes and DEGEN3 dual-controlled instructions -1.7%, but alternating CPU-2 solve 0.612644 -> 0.520240 s is within 0.144281/0.146594 s spreads; DEGEN2 and the first CPU-4 batch also stay within spread. Reverted (`target/perf-lzip-{target,recheck}-{before,after}.txt`).

- Query row lists for former B-column references instead of rebuilding all flags: fast parity 940/940 passes, PILOT4 instructions -0.3%, but alternating CPU-2 solve 0.095259 -> 0.097082 s has no gain (spreads 0.000879/0.004516 s). Reverted (`target/perf-rowrefs-{target,recheck}-{before,after}.txt`).

- Skip capacity counts when the unfiltered U-row length is below capacity: fast parity 940/940 passes, PILOT4 instructions -0.5%, but alternating CPU-2 solve 0.094382 -> 0.099878 s has no gain (spreads 0.002361/0.010112 s). Reverted (`target/perf-rowcapguard-{target,recheck}-{before,after}.txt`).

- Reserve the dense L-stage output capacity once: fast parity 940/940 passes and DEGEN3 dual instructions -0.8%, but CPU-4 solve 0.345439 -> 0.340777 s is within the old 0.007061 s spread; CPU-2 has no timing gain. Reverted (`target/perf-lout-{target,recheck}-{before,after}.txt`).

- Read growth-test diagonal without constructing unused diagnostic L columns: fast parity 940/940 passes and DEGEN2 instructions -0.4%, but CPU-4 solve medians are identical at 0.045557 s (spreads 0.000079/0.000667 s); CPU-2 timing also stays within spread. Reverted (`target/perf-diagread-{target,recheck}-{before,after}.txt`).

- Remove dense BTRAN L-transpose bounds checks using immutable-index invariants: fast parity 940/940 passes and DEGEN2 instructions -0.6%, but alternating CPU-4 solve 0.070389 -> 0.070534 s has no gain (spreads 0.001005/0.000560 s). Reverted (`target/perf-ltindex-recheck-{before,after}.txt`).

- Force-inline existing F80 add/sub: fast parity 940/940 passes, but Quad_2717 instructions +2.1% and solve 0.443852 -> 0.448784 s has no gain; fx_quad instructions +1.1%, solve +3.2%. Reverted (`target/perf-inlineadd-target-{before,after}.txt`).

- Return the larger finite F80 addend when exponent gap >= 66: frozen boundary/random/captured differential tests and fast parity 940/940 pass, but Quad_2717 instructions +1.4% and solve 0.454318 -> 0.475658 s (+4.7%, beyond both spreads). Reverted (`target/perf-distantadd-target-{before,after}.txt`).

- Fast F80 multiplication by finite powers of two for normal results: frozen/random/captured differential tests and fast parity 940/940 pass, but Quad_2717 instructions +0.5% and solve 0.709394 -> 0.712114 s has no gain. Reverted (`target/perf-powmul-target-{before,after}.txt`).

- Cache immutable extended row-wise coefficients on first use: frozen differential and fast parity 940/940 pass, Quad_2717 instructions -1.3%, but solve 0.708071 -> 0.702991 s is within 0.009744/0.035976 s spreads. Reverted (`target/perf-qcoef-target-{before,after}.txt`).

- Mask unused extended cost-change row products: frozen differential and fast parity 940/940 pass, Quad_2717 instructions -1.4%, but alternating CPU-4 solve 0.701696 -> 0.671770 s is within 0.079075/0.140497 s spreads. Reverted (`target/perf-qmask-recheck-{before,after}.txt`).

- Rebuild all elimination WORK row counts by sequential U-transpose scans: fast parity 940/940 passes and PILOT4 instructions -1.2%, but alternating CPU-4 solve 0.111296 -> 0.116587 s has no gain (spreads 0.009052/0.003302 s). Reverted in favor of on-demand counts (`target/perf-rcountrows-recheck-{before,after}.txt`).

- Emit dense FTRAN output directly without ord/done arrays: fast parity 940/940 passes, but PILOT4 solve 0.073238 -> 0.074108 s and instructions 737.9 -> 738.6M; DEGEN2 timing is within spread. Reverted (`target/perf-denseout-target-{before,after}.txt`).

- Move completed FTRAN x/input buffers into F2 after output construction instead of cloning: fast parity 940/940 passes, PILOT4 instructions 737.8 -> 737.3M (-0.06%), solve 0.073656 -> 0.073330 s is within old 0.000944 s spread; DEGEN2 is within spread and PEROLD changes host bands. Reverted (`target/perf-f2move-target.log`).

- Reuse BTRAN dense x/w buffers: fast parity 940/940 passes, PILOT4 instructions 747.0 -> 745.9M (-0.14%), but solve 0.113494 -> 0.116470 s has no gain (spreads 0.003324/0.002445 s); DEGEN2 is within spread and PEROLD changes host bands with little instruction change. Reverted (`target/perf-btbuf-target.log`).

- Precompute eligibility for row_product_work scans longer than n: fast parity 940/940 passes, PEROLD primal instructions -0.4%, solve 1.342457 -> 1.332925 s is within 0.249063/0.547923 s spreads; DEGEN3 primal-controlled has no timing gain and DEGEN2 instructions rise slightly (`target/perf-rowelig-long-target.log`). Reverted.

- Lazy per-column eligibility memoization in row_product_work: fast parity 940/940 passes, PEROLD primal instructions -1.1%, but CPU-2 solve 1.322285 -> 1.320622 s is within 0.013795 s spread; CPU-0 alternating recheck crosses speed bands with no clear aggregate gain (`target/perf-rowelig-{target,recheck}.log`). Reverted.

- Binary lookup in the U transpose for tau terms: fast parity 940/940 passes; DEGEN3 dual instructions 3769.4 -> 3753.4M (-0.4%), solve 0.739033 -> 0.716611 s is within old 0.025890 s spread. DEGEN2/PEROLD instructions rise slightly. Replaced by direct row-map trial (`target/perf-taulookup-target.log`).
- F80 add aligned with three guard bits, plus direct normal-result rounding: frozen old-add/pack differential test and fast parity 940/940 pass, instructions fall 4977.9 -> 4836.8M (-2.8%), but primal_5959 dual CPU-4 solve 0.555599 -> 0.557742 s (spreads 0.006716/0.011524 s). Alignment alone also has no timing gain on CPUs 2/4 despite -1.7% instructions. Reverted both (`target/perf-add-{target,recheck}.log`, `target/perf-add-direct-{target,recheck}.log`).
- Read eliminate work counts directly from maintained rowl lengths instead of rebuilding counts: parity full 987/987 and release tests pass, but mode_switch.dual-r 22.849 -> 21.851 s is smaller than the old 1.849 s spread (instructions 146886.1 -> 140945.4M, -4.0%). Alternating DEGEN2 primal 0.092/0.090 then 0.091/0.091 s; dual 0.068/0.067 s, marginal against spread. Reverted.

## Tried, breaks parity

(idea, first differing record from `rose-parity check`, why)
