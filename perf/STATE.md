# Performance state

Current state of the work in [PLAN.md](PLAN.md). Rewrite in place; no log (history is `jj log`).

## Numbers

Baseline solver: `8b645007` (same solver as `vlsokskq`). Current solver: `ussxlpst` (elimination WORK counts only the rows reached by the solve).
Parity: full tier 987/987; release tests pass. Three timed runs after one warm-up, hyperfine 2.0.0; instruction counts are millions from perf stat.

`baseline wall` is the fixed initial process-wall median. `paired base wall` is the fresh baseline-binary median from alternating baseline/last-landed/candidate checks on one CPU (CPU 0 unless stated otherwise); it resolves host drift without changing the fixed baseline. `current wall` and `solve` use the new harness; **current/ref uses solve-only time**, `diagnostics.seconds`, against original Gurobi Runtime.

Reference provenance: original-library cached records `../reverse-engineer/runs/{degen2,degen3,perold,pilot4}-{profile}.json` for the ten previously valid cells. New three-run medians use `/tmp/rose-perf-reference/{primal_5959.primal-r,primal_5959.dual-r,quad_2717.dual-r,fx_quad_obj_5371.dual-r}.{1,2,3}.json`; summary `/tmp/rose-perf-reference/summary.json`. All original-library records identify Ryzen 9 5900X and library SHA256 `8e87fca221fdffc922331c8619426031b3e72dfd419446422f88af0cb0e98f9c`.
`-` refs are unreachable: original Gurobi rejects 80BAU3B, GREENBEA and mode-switch with license error 10010 (scratch `ID.error.json`, including direct 80BAU3B primal/dual attempts). The old 80BAU3B primal/dual refs used the patched library and have been removed; the old prose mode-switch time is unusable.

Regression checks: all 23 runs paired against the preceding solver on CPU 4. GREENBEA dual-controlled and PDHG needed alternating CPU-2 follow-ups; both increases are within the combined paired min-max spreads. Dual instructions fall 4.0%, and PDHG instructions stay flat. No unresolved fixed-baseline regression: 80BAU3B barrier remains below its fresh paired baseline wall, and both PDHG runs remain within the fixed or fresh paired baseline bands. Raw paired data `target/perf-rlen-full-{before,after}.txt` and `target/perf-rlen-followup-{before,after}.txt`. Table refreshed on CPU 2 in `target/perf-quad-checkpoint-bench.txt`; the restored release binary is byte-identical to the saved `ussxlpst` solver (SHA256 prefix `5878e9618cce`), with full 987/987 parity and release tests checked again. Every fixed/fresh-baseline comparison stays within its spread; 80BAU3B PDHG spans 6.582..8.306 s wall against fresh baseline 6.526 s. Host speed varies substantially; use immediate pairs and spreads.
Acceptance: CPU-2 PILOT4 primal-release solve 0.073956 -> 0.071406 s (-3.4%; spreads 0.000906/0.000706 s), instructions 737.9 -> 706.8M (-4.2%). DEGEN2 primal improves 0.050818 -> 0.049744 s (-2.1%; spreads 0.000225/0.000275 s), instructions -3.4% (`target/perf-rlen-target-{before,after}.txt`). Elimination uses the ordered numeric U transpose to count distinct tail columns only for rows accepted by the solve, avoiding a full column/stamp rebuild and two scratch allocations. Numeric operations, WORK formulas and bookings, output order and dumps stay unchanged. The full sweep confirms instruction reductions up to 8.0% on DEGEN3 dual; timing changes across host bands are not attributed to this change.

| run | ref solve | baseline wall | paired base wall | current wall | solve | current/ref | instr |
|---|---|---|---|---|---|---|---|
| degen2.primal-c | 0.012711 | 0.182 | - | 0.082 | 0.077265 | 6.08 | 408.8 |
| degen2.dual-r | 0.009653 | 0.125 | - | 0.066 | 0.061630 | 6.38 | 302.4 |
| degen3.primal-c | 0.163137 | 3.237 | - | 1.097 | 1.082140 | 6.63 | 5998.9 |
| degen3.primal-r | 0.193685 | 3.872 | - | 1.233 | 1.220107 | 6.30 | 6946.0 |
| perold.primal-c | 0.160994 | 2.430 | - | 1.129 | 1.122242 | 6.97 | 6388.4 |
| pilot4.primal-r | 0.015787 | 0.253 | - | 0.119 | 0.113744 | 7.20 | 706.7 |
| 80bau3b.primal-c | - | 1.618 | - | 0.811 | 0.795797 | - | 8516.7 |
| primal_5959.primal-r | 0.067884 | 0.563 | 0.938 | 0.280 | 0.277589 | 4.09 | 2791.8 |
| degen3.dual-c | 0.082843 | 1.002 | - | 0.350 | 0.340987 | 4.12 | 2842.7 |
| degen3.dual-r | 0.082947 | 1.008 | 1.007 | 0.351 | 0.341007 | 4.11 | 2853.1 |
| perold.dual-r | 0.038532 | 0.355 | 0.360 | 0.173 | 0.168878 | 4.38 | 1521.7 |
| 80bau3b.dual-r | - | 1.103 | 1.688 | 0.430 | 0.415081 | - | 3960.5 |
| greenbea.dual-c | - | 2.870 | 2.363 | 0.958 | 0.944783 | - | 8348.9 |
| greenbea.dual-r | - | 3.598 | 3.592 | 0.947 | 0.933619 | - | 8272.6 |
| mode_switch.dual-r | - | 25.815 | - | 4.185 | 4.143810 | - | 44569.3 |
| primal_5959.dual-r | 0.098192 | 1.300 | - | 0.622 | 0.619040 | 6.30 | 3913.8 |
| quad_2717.dual-r | 0.092814 | 1.130 | - | 0.695 | 0.691440 | 7.45 | 4035.2 |
| fx_quad_obj_5371.dual-r | 0.038077 | 0.334 | 0.514 | 0.258 | 0.255327 | 6.71 | 1575.9 |
| degen3.barrier-r | 0.075214 | 0.608 | 0.856 | 0.465 | 0.405322 | 5.39 | 3779.9 |
| 80bau3b.barrier-r | - | 0.379 | 0.546 | 0.401 | 0.365763 | - | 3425.7 |
| greenbea.barrier-r | - | 1.569 | 2.493 | 1.537 | 1.486657 | - | 13028.4 |
| 80bau3b.pdhg | - | 4.552 | 6.526 | 6.640 | 6.620413 | - | 92390.8 |
| greenbea.pdhg | - | 2.801 | - | 1.917 | 1.905439 | - | 37740.6 |

## Hotspots

Per profiled bench run: the commit it was taken at and the top self-time functions (percent), at most
five lines each. Re-profile a run when a landed change touched its hot path.

- perold.dual-r, `vsxqqkwy` (stale after dense slot counting; frame-pointer build passes full 987/987 parity; ten solves capped at 1300 iterations, `target/perf-perold-rest.data`): FTRAN 8.75%, update 8.05%, iteration 7.63%, ftran_needed 6.42%, LU nucleus 5.89% self; L stage 5.72%. Integer slot sorting remains 1.80% in small-sort networks plus 1.25% in quicksort. All ten optimal in 1273 iterations, 27 refactorizations. Large dense slot lists now use reusable counts and a descending slot scan; smaller/sparser lists retain comparison sorting. Reprofile before choosing further work.
- perold.primal-c, `yryslslw` (stale after dense work-count shortcut; frame-pointer build passes full 987/987 parity; three solves capped at 3500 iterations, `target/perf-perold.data`): BTRAN 8.71%, primal control 7.01%, row_product_work 6.39%, update 6.12%, fresh factor 6.10% self. Refactor 26.63%, FTRAN 11.82%, BTRAN 11.03% inclusive. All three optimal in 3399 iterations, 90 refactorizations; original reference also has 3399 iterations. The work-estimator annotation (`target/perf-perold-rowprod-annotate.txt`) attributes substantial time to checks in the touched-row scan; dense counts now avoid that scan. Reprofile before choosing the remaining sparse work.
- degen3.primal-r, `ltkyqoqr` (stale after member-based row ordering; frame-pointer build passes full 987/987 parity; three solves capped at 5300 iterations, `target/perf-primal3-next.data`): update 10.65%, FTRAN 9.69%, BTRAN 7.93%, fresh factor 4.59%, fresh row-list filtering 4.32% self. Signature classes remain 2.99%. All three optimal in 5165 iterations, 66 refactorizations. The full basis-order scan for fresh singleton-row lists is now replaced by sorting the row members; reprofile before choosing more work. Update annotation (`target/perf-primal3-next-update-annotate.txt`) still highlights eliminate row-count reconstruction, whose rowl-length replacement was already rejected below.
- quad_2717.dual-r, `ussxlpst` (frame-pointer build passes full 987/987 parity; eight solves capped at 17000 iterations, `target/perf-quad-next.data`): F80 add 15.13%, primal control 10.38%, F80 multiply 8.89%, pack 5.06%, rank_products_q 4.22% self. Row products are 13.20%, cost_change_work 12.77%, refactorization 17.26% inclusive. All eight optimal in 16010 iterations, 3022 refactorizations. Cost-change extended products are consumed only for nonbasic variables whose double product exceeds ROW_DROP_TOL; masking these outputs and caching coefficient conversions both failed the timing gate below. F80 add annotation (`target/perf-quad-next-add-annotate.txt`) still shows normalization and shift work; the tested distant-exponent shortcut and forced add inlining both regress.
- degen3.barrier-r, `munmovyt` (stale after kernel dispatch; frame-pointer build passes full parity; three solves capped at 100 iterations, `target/perf-barrier-degen3-valid.data`): FMA dispatch 26.46%, DGEMM 24.62%, existing hardware FMA implementation 14.77%, numeric_factor 7.40%, DTRSM 4.36% self. DGEMM 49.04%, numeric factor 81.45% inclusive; all three solves optimal in 14 iterations. The per-product FMA wrapper is an indirect call plus frame/return overhead (annotated `target/perf-barrier-fma-annotate.txt`); per-call dispatch now removes that per-product overhead without adding fused operations. BORE3D is a cheap witness (16 barrier iterations, 12 ms).
- primal_5959.dual-r, `tplnupps` (stale after factor reuse; frame-pointer build passes full parity; three solves capped at 16000 iterations, `target/perf-dual5959-full.data`): F80 add_signed 7.50%, primal control 5.70%, F80 multiply 4.54%, F80 nucleus 3.88%, double nucleus 3.74% self. All three reach numerical_failure, 15688 iterations, 7842 refactorizations. Refactorization is 51.10% inclusive, quad mirrors 21.69% before factor reuse. LUIN capture `target/perf-dual5959-luin.log` has 7822 records and 40 input-head/threshold combinations; two-entry LRU estimate 99.3% hits. Wide-shift arithmetic rewrites are rejected below. The earlier 11000-iteration capture stopped early and is not the full-run profile.
- pilot4.primal-r, `oqoyqxpq` (stale after on-demand elimination WORK counts; frame-pointer build passes full 987/987 parity; twenty solves capped at 700 iterations, `target/perf-pilot-next2.data`): update 9.97%, FTRAN 7.80%, primal control 7.77%, fresh factor 6.78%, BTRAN 5.55% self. malloc/free total 3.30% self. All twenty optimal in 690 iterations, 16 refactorizations. Update annotation (`target/perf-pilot-next2-update.txt`) puts about 30% of update self samples in the full elimination WORK row-count rebuild, now replaced by on-demand transpose scans. Reprofile before choosing the remaining update work.
- primal_5959.primal-r, `kyyvxyql` (stale after twin reuse; frame-pointer build passes full 987/987 parity; three bounded solves at 1999 Hz, `target/perf-fp-primal5959.data`, report with `--no-inline`): primal control 6.32%, F80 nucleus 5.46%, fresh F80 factor 5.09%, F80 add_signed 6.67% self. Refactor 53.53%, quad mirrors 28.88% inclusive before twin reuse. Log verifies numerical_failure, 10780 iterations, 5381 refactorizations; original reference status 12, 10780 iterations, refactor count unavailable. Default inline expansion duplicates frames; use PLAN.md's validated frame-pointer workflow for caller attribution.
- fx_quad_obj_5371.dual-r, `kyyvxyql` (stale after twin reuse; profiling build, three solves at 1999 Hz, `target/perf-members-fx.data`): F80 add_signed 8.06%, F80 nucleus 5.81%, multiply 5.74%, fresh F80 factor 5.17%, pack 3.83% self cycles. Solve remains infeasible, 5319 iterations, 5 refactorizations.
- degen3.dual-r, `sormxnyt` (stale after metadata caching; frame-pointer build passes full parity; three solves capped at 2300 iterations, `target/perf-btplan-degen3.data`): iteration 13.36%, update 8.89%, L stage 7.68%, BTRAN 7.21%, FTRAN 6.93% self; ftran_needed 5.34%, with the old column-scan filler gone. BTRAN annotation (`target/perf-btplan-btran-annotate.txt`) attributes over a quarter of its self samples to rebuilding transpose levels, and more to rebuilding the inverse row map. L columns, initial row permutation and loff are immutable between refactorizations. All solves optimal in 2212 iterations, 24 refactorizations; original reference also has 2212 iterations.
- mode_switch.dual-r, `xzuktkmx` (stale after F80/LU changes; reprofile before choosing it): update 22.21%, LU nucleus 15.90%, BTRAN 13.33%, fresh-factor row-list filtering 7.67%, factor_structure 3.36% self cycles. Full solve: 11393 iterations, 17 refactorizations.

## In progress

None.

## Candidates

Ordered by current/ref; reprofile any changed path before choosing a change. Ratios use the current
sweep's solve medians; paired old/new times determine whether a change gains or regresses.

- quad_2717.dual-r (7.45x), pilot4.primal-r (7.20x), perold.primal-c (6.97x), fx_quad_obj_5371.dual-r (6.71x), degen3.primal-c (6.63x): work in this order. Quad_2717 is freshly profiled above; reprofile other changed paths before choosing a change. Host speed bands can reorder these snapshot ratios.
- DEGEN3 dual: tau row-value lookup and immutable BTRAN metadata are settled. Reprofile before inspecting remaining iteration control, updates and scratch allocation.
- Other reachable runs follow their table ratios. Quad twins share construction and four compact repeated factors are cached; row-capacity growth already uses the numeric U transpose.
- Quad_2717: the current profile is refreshed above. F80 arithmetic and set_d_from_y remain targets; inspect callers and assembly for repeated conversion or state reconstruction. The measured/rejected trials are listed below. Cost-change work performs required reduced-cost updates; the factor-cache working set is settled.
- Barrier: profile again before choosing more work. DGEMM now dispatches existing fused operations once per call; other BLAS-like kernels still use their original scalar dispatch and accumulation order. A smaller remaining DGEMM gap may involve independent outputs and data reuse, which must preserve every per-output reduction order.
- Runs without reachable references, after the reachable ratios: mode-switch update (22.21%) still scans A/B lists for other memberships; BTRAN (13.33% before metadata caching) allocates scratch; LU nucleus (15.90%); fresh-factor row-list filtering (7.67%). These old profiles are stale; annotate again before choosing a change.
- all: dump::on/record construction without dumps and allocation in FTRAN/BTRAN, pricing and ratio tests remain unmeasured leads.

## Tried, no gain

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
