# Performance state

Current state of the work in [PLAN.md](PLAN.md). Rewrite in place; no log (history is `jj log`).

## Numbers

Baseline solver: `8b645007` (same solver as `vlsokskq`). Current solver: `vsxqqkwy` (descending-slot sorts map through inverse slot/at permutations).
Parity: full tier 987/987; release tests pass. Three timed runs after one warm-up, hyperfine 2.0.0; instruction counts are millions from perf stat.

`baseline wall` is the fixed initial process-wall median. `paired base wall` is the fresh baseline-binary median from alternating baseline/last-landed/candidate checks on one CPU (CPU 0 unless stated otherwise); it resolves host drift without changing the fixed baseline. `current wall` and `solve` use the new harness; **current/ref uses solve-only time**, `diagnostics.seconds`, against original Gurobi Runtime.

Reference provenance: original-library cached records `../reverse-engineer/runs/{degen2,degen3,perold,pilot4}-{profile}.json` for the ten previously valid cells. New three-run medians use `/tmp/rose-perf-reference/{primal_5959.primal-r,primal_5959.dual-r,quad_2717.dual-r,fx_quad_obj_5371.dual-r}.{1,2,3}.json`; summary `/tmp/rose-perf-reference/summary.json`. All original-library records identify Ryzen 9 5900X and library SHA256 `8e87fca221fdffc922331c8619426031b3e72dfd419446422f88af0cb0e98f9c`.
`-` refs are unreachable: original Gurobi rejects 80BAU3B, GREENBEA and mode-switch with license error 10010 (scratch `ID.error.json`, including direct 80BAU3B primal/dual attempts). The old 80BAU3B primal/dual refs used the patched library and have been removed; the old prose mode-switch time is unusable.

Regression checks: all 23 runs paired against the preceding solver on CPU 4. Five apparent increases resolve in three alternating CPU-2 pairs: DEGEN3 primal-release solve 1.368385/1.319370 s; 80BAU3B primal-controlled 0.818883/1.262016 s (new spread 0.655586 s); mode-switch 6.672657/4.496711 s; DEGEN3 barrier 0.277706/0.276401 s; 80BAU3B PDHG 6.330111/6.375440 s (old spread 2.114518 s). All improve or remain within their paired spreads. Mode-switch/PDHG cross host-speed bands with nearly unchanged instructions; do not attribute their large timing shifts to the optimization. No new fixed-baseline regression: PDHG remains below its fresh paired baseline wall; 80BAU3B barrier is below the fixed baseline. Raw data `target/perf-slotsort-full-{before,after}.txt`, `target/perf-slotsort-followup-{before,after}.txt`. Host speed varies substantially; use immediate pairs and spreads.
Acceptance: CPU-2 PILOT4 primal-release solve 0.077738 -> 0.075848 s (-2.4%; spreads 0.000375/0.000302 s), instructions 759.2 -> 754.5M (-0.6%) (`target/perf-slotsort-target.log`). The full sweep confirms 0.120803 -> 0.117953 s, beyond either spread. PEROLD dual improves 0.191253 -> 0.181738 s (-5.0%) with -1.8% instructions; DEGEN2 primal improves 0.053022 -> 0.052033 s. FTRAN/FTRAN-needed/U-order descending-slot sorts map positions to their slots in place, sort the integers, then map back through at. Because slot/at are inverse permutations, the exact position order and repeated positions are preserved without repeated indirect comparator-key reads or extra allocations. Floating-point operations, traversal order, dumps and WORK remain unchanged.

| run | ref solve | baseline wall | paired base wall | current wall | solve | current/ref | instr |
|---|---|---|---|---|---|---|---|
| degen2.primal-c | 0.012711 | 0.182 | - | 0.055 | 0.052070 | 4.10 | 432.9 |
| degen2.dual-r | 0.009653 | 0.125 | - | 0.044 | 0.040758 | 4.22 | 321.4 |
| degen3.primal-c | 0.163137 | 3.237 | - | 0.739 | 0.729153 | 4.47 | 6471.0 |
| degen3.primal-r | 0.193685 | 3.872 | - | 1.333 | 1.319370 | 6.81 | 7650.9 |
| perold.primal-c | 0.160994 | 2.430 | - | 1.176 | 1.169645 | 7.27 | 6793.9 |
| pilot4.primal-r | 0.015787 | 0.253 | - | 0.123 | 0.117953 | 7.47 | 754.5 |
| 80bau3b.primal-c | - | 1.618 | - | 1.284 | 1.262016 | - | 8710.6 |
| primal_5959.primal-r | 0.067884 | 0.563 | 0.938 | 0.454 | 0.450967 | 6.64 | 2822.4 |
| degen3.dual-c | 0.082843 | 1.002 | - | 0.590 | 0.572583 | 6.91 | 3150.5 |
| degen3.dual-r | 0.082947 | 1.008 | 1.007 | 0.598 | 0.583173 | 7.03 | 3161.5 |
| perold.dual-r | 0.038532 | 0.355 | 0.360 | 0.295 | 0.288614 | 7.49 | 1647.8 |
| 80bau3b.dual-r | - | 1.103 | 1.688 | 0.718 | 0.695368 | - | 4179.2 |
| greenbea.dual-c | - | 2.870 | 2.363 | 1.024 | 1.011133 | - | 8858.0 |
| greenbea.dual-r | - | 3.598 | 3.592 | 1.014 | 1.000465 | - | 8778.3 |
| mode_switch.dual-r | - | 25.815 | - | 4.531 | 4.496711 | - | 50474.9 |
| primal_5959.dual-r | 0.098192 | 1.300 | - | 0.640 | 0.636384 | 6.48 | 3959.5 |
| quad_2717.dual-r | 0.092814 | 1.130 | - | 0.692 | 0.689395 | 7.43 | 4085.4 |
| fx_quad_obj_5371.dual-r | 0.038077 | 0.334 | 0.514 | 0.264 | 0.260667 | 6.85 | 1591.7 |
| degen3.barrier-r | 0.075214 | 0.608 | 0.856 | 0.317 | 0.276401 | 3.67 | 3779.3 |
| 80bau3b.barrier-r | - | 0.379 | 0.546 | 0.281 | 0.252521 | - | 3425.3 |
| greenbea.barrier-r | - | 1.569 | 2.493 | 1.006 | 0.972808 | - | 13029.9 |
| 80bau3b.pdhg | - | 4.552 | 6.526 | 6.394 | 6.375440 | - | 92390.0 |
| greenbea.pdhg | - | 2.801 | - | 2.766 | 2.748431 | - | 37741.1 |

## Hotspots

Per profiled bench run: the commit it was taken at and the top self-time functions (percent), at most
five lines each. Re-profile a run when a landed change touched its hot path.

- perold.dual-r, `vnqmrmnq` (stale after count reuse; frame-pointer build passes full 987/987 parity; ten solves capped at 1300 iterations, `target/perf-dualwork.data`): iteration 9.49%, FTRAN 8.55%, update 7.80%, ftran_needed 6.68%, L stage 5.73% self. All ten optimal in 1273 iterations, 27 refactorizations. The duplicate distinct-column count in the annotated iteration (`target/perf-dualwork-iteration-annotate.txt`) is now removed, using the count before pivot-row filtering. Reprofile before choosing further work.
- perold.primal-c, `yryslslw` (stale after dense work-count shortcut; frame-pointer build passes full 987/987 parity; three solves capped at 3500 iterations, `target/perf-perold.data`): BTRAN 8.71%, primal control 7.01%, row_product_work 6.39%, update 6.12%, fresh factor 6.10% self. Refactor 26.63%, FTRAN 11.82%, BTRAN 11.03% inclusive. All three optimal in 3399 iterations, 90 refactorizations; original reference also has 3399 iterations. The work-estimator annotation (`target/perf-perold-rowprod-annotate.txt`) attributes substantial time to checks in the touched-row scan; dense counts now avoid that scan. Reprofile before choosing the remaining sparse work.
- degen3.primal-r, `ltkyqoqr` (stale after member-based row ordering; frame-pointer build passes full 987/987 parity; three solves capped at 5300 iterations, `target/perf-primal3-next.data`): update 10.65%, FTRAN 9.69%, BTRAN 7.93%, fresh factor 4.59%, fresh row-list filtering 4.32% self. Signature classes remain 2.99%. All three optimal in 5165 iterations, 66 refactorizations. The full basis-order scan for fresh singleton-row lists is now replaced by sorting the row members; reprofile before choosing more work. Update annotation (`target/perf-primal3-next-update-annotate.txt`) still highlights eliminate row-count reconstruction, whose rowl-length replacement was already rejected below.
- quad_2717.dual-r, `zzvmznmp` (stale after cache sizing; frame-pointer build passes full parity; three solves capped at 20000 iterations, `target/perf-quad2717.data`): F80 add 13.58%, primal control 8.94%, F80 multiply 8.17%, pack 3.80%, rank_products_q 3.61% self. Refactorization 25.94%, quad mirrors 11.50%, rank_products_q 11.26%, cost_change_work 11.18% inclusive. All solves optimal in 16010 iterations with 3022 refactorizations. `target/perf-quad2717-luin.log`: 3015 LUIN records, 14 distinct input-head/threshold keys; FIFO hits 110/2998/3000 at capacities 2/4/8. The cache now holds four entries; larger capacities add only marginal hits. Cost-change work also updates reduced costs and cannot simply be removed.
- degen3.barrier-r, `munmovyt` (stale after kernel dispatch; frame-pointer build passes full parity; three solves capped at 100 iterations, `target/perf-barrier-degen3-valid.data`): FMA dispatch 26.46%, DGEMM 24.62%, existing hardware FMA implementation 14.77%, numeric_factor 7.40%, DTRSM 4.36% self. DGEMM 49.04%, numeric factor 81.45% inclusive; all three solves optimal in 14 iterations. The per-product FMA wrapper is an indirect call plus frame/return overhead (annotated `target/perf-barrier-fma-annotate.txt`); per-call dispatch now removes that per-product overhead without adding fused operations. BORE3D is a cheap witness (16 barrier iterations, 12 ms).
- primal_5959.dual-r, `tplnupps` (stale after factor reuse; frame-pointer build passes full parity; three solves capped at 16000 iterations, `target/perf-dual5959-full.data`): F80 add_signed 7.50%, primal control 5.70%, F80 multiply 4.54%, F80 nucleus 3.88%, double nucleus 3.74% self. All three reach numerical_failure, 15688 iterations, 7842 refactorizations. Refactorization is 51.10% inclusive, quad mirrors 21.69% before factor reuse. LUIN capture `target/perf-dual5959-luin.log` has 7822 records and 40 input-head/threshold combinations; two-entry LRU estimate 99.3% hits. Wide-shift arithmetic rewrites are rejected below. The earlier 11000-iteration capture stopped early and is not the full-run profile.
- pilot4.primal-r, `mzxvtkxy` (stale after slot sorting; frame-pointer build passes full 987/987 parity; twenty solves capped at 700 iterations, `target/perf-pilot-rest.data`): update 9.37%, FTRAN 8.05%, primal control 7.65%, BTRAN 6.79%, fresh factor 5.17% self. FTRAN slot sorting adds 2.01% self; row-count reconstruction no longer appears among the top samples after revalidation. All twenty optimal in 690 iterations, 16 refactorizations. Repeated comparator-key lookups in descending-slot sorts are now replaced by in-place slot/at mapping; reprofile before choosing further work.
- primal_5959.primal-r, `kyyvxyql` (stale after twin reuse; frame-pointer build passes full 987/987 parity; three bounded solves at 1999 Hz, `target/perf-fp-primal5959.data`, report with `--no-inline`): primal control 6.32%, F80 nucleus 5.46%, fresh F80 factor 5.09%, F80 add_signed 6.67% self. Refactor 53.53%, quad mirrors 28.88% inclusive before twin reuse. Log verifies numerical_failure, 10780 iterations, 5381 refactorizations; original reference status 12, 10780 iterations, refactor count unavailable. Default inline expansion duplicates frames; use PLAN.md's validated frame-pointer workflow for caller attribution.
- fx_quad_obj_5371.dual-r, `kyyvxyql` (stale after twin reuse; profiling build, three solves at 1999 Hz, `target/perf-members-fx.data`): F80 add_signed 8.06%, F80 nucleus 5.81%, multiply 5.74%, fresh F80 factor 5.17%, pack 3.83% self cycles. Solve remains infeasible, 5319 iterations, 5 refactorizations.
- degen3.dual-r, `sormxnyt` (stale after metadata caching; frame-pointer build passes full parity; three solves capped at 2300 iterations, `target/perf-btplan-degen3.data`): iteration 13.36%, update 8.89%, L stage 7.68%, BTRAN 7.21%, FTRAN 6.93% self; ftran_needed 5.34%, with the old column-scan filler gone. BTRAN annotation (`target/perf-btplan-btran-annotate.txt`) attributes over a quarter of its self samples to rebuilding transpose levels, and more to rebuilding the inverse row map. L columns, initial row permutation and loff are immutable between refactorizations. All solves optimal in 2212 iterations, 24 refactorizations; original reference also has 2212 iterations.
- mode_switch.dual-r, `xzuktkmx` (stale after F80/LU changes; reprofile before choosing it): update 22.21%, LU nucleus 15.90%, BTRAN 13.33%, fresh-factor row-list filtering 7.67%, factor_structure 3.36% self cycles. Full solve: 11393 iterations, 17 refactorizations.

## In progress

(none)

## Candidates

Ordered by current/ref; reprofile any changed path before choosing a change. Ratios use the current
sweep's solve medians; paired old/new times determine whether a change gains or regresses.

- perold.dual-r (7.49x), pilot4.primal-r (7.47x), quad_2717.dual-r (7.43x), perold.primal-c (7.27x): reprofile in this order. Profiles of changed paths are stale.
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
