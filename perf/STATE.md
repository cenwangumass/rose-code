# Performance state

Current state of the work in [PLAN.md](PLAN.md). Rewrite in place; no log (history is `jj log`).

## Numbers

Baseline solver: `8b645007` (same solver as `vlsokskq`). Current solver: `szutlnuy` (fresh singleton-row lists sorted from their members rather than full basis-order scans).
Parity: full tier 987/987; release tests pass. Three timed runs after one warm-up, hyperfine 2.0.0; instruction counts are millions from perf stat.

`baseline wall` is the fixed initial process-wall median. `paired base wall` is the fresh baseline-binary median from alternating baseline/last-landed/candidate checks on one CPU (CPU 0 unless stated otherwise); it resolves host drift without changing the fixed baseline. `current wall` and `solve` use the new harness; **current/ref uses solve-only time**, `diagnostics.seconds`, against original Gurobi Runtime.

Reference provenance: original-library cached records `../reverse-engineer/runs/{degen2,degen3,perold,pilot4}-{profile}.json` for the ten previously valid cells. New three-run medians use `/tmp/rose-perf-reference/{primal_5959.primal-r,primal_5959.dual-r,quad_2717.dual-r,fx_quad_obj_5371.dual-r}.{1,2,3}.json`; summary `/tmp/rose-perf-reference/summary.json`. All original-library records identify Ryzen 9 5900X and library SHA256 `8e87fca221fdffc922331c8619426031b3e72dfd419446422f88af0cb0e98f9c`.
`-` refs are unreachable: original Gurobi rejects 80BAU3B, GREENBEA and mode-switch with license error 10010 (scratch `ID.error.json`, including direct 80BAU3B primal/dual attempts). The old 80BAU3B primal/dual refs used the patched library and have been removed; the old prose mode-switch time is unusable.

Regression checks: all 23 runs paired against the preceding solver on CPU 4. Two apparent increases resolve within the spreads of three alternating CPU-2 pairs: mode-switch solve 6.488801/6.760096 s (spreads 3.727506/1.741015 s), instructions 61052.6/50575.1M (-17.2%); 80BAU3B PDHG 4.469287/6.413085 s (old spread 1.979182 s), with unchanged instructions. The first PDHG pair is 6.410000/6.413085 s; later pairs cross host-speed bands. Do not attribute mode-switch's band changes to this optimization. The initial PEROLD primal increase resolves in a separate alternating check: 1.217551/1.212041 s (old spread 0.071061 s), -1.4% instructions. No new fixed-baseline regression: PDHG and 80BAU3B barrier remain below their fresh paired baseline walls. Raw data `target/perf-roworder-full-{before,after}.txt`, `target/perf-roworder-recheck-{before,after}.txt`, `target/perf-roworder-target-recheck-{before,after}.txt`. Host speed varies substantially; use immediate pairs and spreads.
Acceptance: CPU-2 DEGEN3 primal-release solve 1.536586 -> 1.429188 s (-7.0%; spreads 0.048893/0.011173 s), instructions 8738.6 -> 7824.2M (-10.5%) (`target/perf-roworder-target.log`). The full sweep confirms 0.995624 -> 0.914121 s (-8.2%), beyond either spread. For a fresh row without a stored U row, filtering the complete basis-column order is equivalent to sorting only its members by basis column and emitting each position once. The position map filters unrepresented positions as before; stored-U-row and pair-row paths remain unchanged. Every row list, number and WORK record stays exact.

| run | ref solve | baseline wall | paired base wall | current wall | solve | current/ref | instr |
|---|---|---|---|---|---|---|---|
| degen2.primal-c | 0.012711 | 0.182 | - | 0.090 | 0.085479 | 6.72 | 439.2 |
| degen2.dual-r | 0.009653 | 0.125 | - | 0.076 | 0.071426 | 7.40 | 337.6 |
| degen3.primal-c | 0.163137 | 3.237 | - | 0.765 | 0.754935 | 4.63 | 6567.7 |
| degen3.primal-r | 0.193685 | 3.872 | - | 0.924 | 0.914121 | 4.72 | 7823.9 |
| perold.primal-c | 0.160994 | 2.430 | - | 0.785 | 0.780529 | 4.85 | 6926.1 |
| pilot4.primal-r | 0.015787 | 0.253 | - | 0.086 | 0.082616 | 5.23 | 798.0 |
| 80bau3b.primal-c | - | 1.618 | - | 0.832 | 0.817150 | - | 8725.5 |
| primal_5959.primal-r | 0.067884 | 0.563 | 0.938 | 0.435 | 0.432304 | 6.37 | 2823.6 |
| degen3.dual-c | 0.082843 | 1.002 | - | 0.661 | 0.646835 | 7.81 | 3348.1 |
| degen3.dual-r | 0.082947 | 1.008 | 1.007 | 0.660 | 0.645836 | 7.79 | 3358.1 |
| perold.dual-r | 0.038532 | 0.355 | 0.360 | 0.331 | 0.323941 | 8.41 | 1768.0 |
| 80bau3b.dual-r | - | 1.103 | 1.688 | 0.787 | 0.764637 | - | 4361.9 |
| greenbea.dual-c | - | 2.870 | 2.363 | 1.833 | 1.814546 | - | 9737.2 |
| greenbea.dual-r | - | 3.598 | 3.592 | 1.880 | 1.860153 | - | 9607.9 |
| mode_switch.dual-r | - | 25.815 | - | 6.809 | 6.760096 | - | 50575.1 |
| primal_5959.dual-r | 0.098192 | 1.300 | - | 0.629 | 0.625231 | 6.37 | 3960.3 |
| quad_2717.dual-r | 0.092814 | 1.130 | - | 0.710 | 0.707115 | 7.62 | 4086.2 |
| fx_quad_obj_5371.dual-r | 0.038077 | 0.334 | 0.514 | 0.266 | 0.263757 | 6.93 | 1592.3 |
| degen3.barrier-r | 0.075214 | 0.608 | 0.856 | 0.470 | 0.410933 | 5.46 | 3778.3 |
| 80bau3b.barrier-r | - | 0.379 | 0.546 | 0.412 | 0.377377 | - | 3426.7 |
| greenbea.barrier-r | - | 1.569 | 2.493 | 1.005 | 0.971828 | - | 13029.4 |
| 80bau3b.pdhg | - | 4.552 | 6.526 | 6.432 | 6.413085 | - | 92390.3 |
| greenbea.pdhg | - | 2.801 | - | 2.720 | 2.703512 | - | 37741.2 |

## Hotspots

Per profiled bench run: the commit it was taken at and the top self-time functions (percent), at most
five lines each. Re-profile a run when a landed change touched its hot path.

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

- perold.dual-r (8.41x), degen3.dual-c (7.81x), degen3.dual-r (7.79x), quad_2717.dual-r (7.62x): reprofile in this order. Profiles of changed paths are stale.
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
