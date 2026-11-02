# Performance state

Current state of the work in [PLAN.md](PLAN.md). Rewrite in place; no log (history is `jj log`).

## Numbers

Baseline solver: `8b645007` (same solver as `vlsokskq`). Current solver: `twspzytl` (row-indexed references for A-list demotion).
Parity: full tier 987/987; release tests pass. Three timed runs after one warm-up, hyperfine 2.0.0; instruction counts are millions from perf stat.

`baseline wall` is the fixed initial process-wall median. `paired base wall` is the fresh baseline-binary median from alternating baseline/last-landed/candidate checks on one CPU (CPU 0 unless stated otherwise); it resolves host drift without changing the fixed baseline. `current wall` and `solve` use the new harness; **current/ref uses solve-only time**, `diagnostics.seconds`, against original Gurobi Runtime.

Reference provenance: original-library cached records `../reverse-engineer/runs/{degen2,degen3,perold,pilot4}-{profile}.json` for the ten previously valid cells. New three-run medians use `/tmp/rose-perf-reference/{primal_5959.primal-r,primal_5959.dual-r,quad_2717.dual-r,fx_quad_obj_5371.dual-r}.{1,2,3}.json`; summary `/tmp/rose-perf-reference/summary.json`. All original-library records identify Ryzen 9 5900X and library SHA256 `8e87fca221fdffc922331c8619426031b3e72dfd419446422f88af0cb0e98f9c`.
`-` refs are unreachable: original Gurobi rejects 80BAU3B, GREENBEA and mode-switch with license error 10010 (scratch `ID.error.json`, including direct 80BAU3B primal/dual attempts). The old 80BAU3B primal/dual refs used the patched library and have been removed; the old prose mode-switch time is unusable.

Regression checks: all 23 runs paired against the preceding solver on CPU 2, alternating pair order by run. Mode-switch dual and GREENBEA PDHG timing flags were followed up on CPU 4. PDHG resolves within spread (1.879243 -> 1.885962 s, effectively identical instructions). Mode-switch follow-ups cross speed bands: stable second pair 2.276170 -> 1.963831 s (-13.72%; spreads 0.022391/0.022752 s), other pairs mixed; aggregate timing is inconclusive, with no consistent regression and an 11% instruction reduction. Table uses final follow-up times for these two rows and full-sweep times elsewhere. No repeatable regression established. Raw data `target/perf-demoterefs-full-{before,after}.txt`, `target/perf-demoterefs-followup-{before,after}.txt`, `target/perf-demoterefs-pdhg-{before,after}.txt`.
Acceptance: PEROLD dual instructions fall about 0.65%. Stable confirmation pairs improve beyond both spreads: CPU 4 0.230011 -> 0.223584 s (spreads 0.002119/0.002179 s); CPU 2 0.231271 -> 0.227321 s (0.003152/0.002335 s) and 0.231521 -> 0.226497 s (0.004112/0.002962 s). Other confirmation pairs have wide spreads. Full-sweep PEROLD 0.147430 -> 0.146186 s is within its wide baseline spread. GREENBEA dual-release instructions fall about 1.8%; mode-switch about 11%. Demotion now filters the existing row index, sorts and deduplicates references, preserving ascending column order instead of scanning every column twice. Full parity 987/987 and release tests pass. Raw target/confirmation data `target/perf-demoterefs-{target,confirm,recheck}-{before,after}.txt`.

| run | ref solve | baseline wall | paired base wall | current wall | solve | current/ref | instr |
|---|---|---|---|---|---|---|---|
| degen2.primal-c | 0.012711 | 0.182 | - | 0.071 | 0.066407 | 5.22 | 346.4 |
| degen2.dual-r | 0.009653 | 0.125 | - | 0.060 | 0.055066 | 5.70 | 275.8 |
| degen3.primal-c | 0.163137 | 3.237 | - | 0.877 | 0.862224 | 5.29 | 4851.5 |
| degen3.primal-r | 0.193685 | 3.872 | - | 1.037 | 1.022661 | 5.28 | 5560.2 |
| perold.primal-c | 0.160994 | 2.430 | - | 0.916 | 0.908863 | 5.65 | 5344.0 |
| pilot4.primal-r | 0.015787 | 0.253 | - | 0.098 | 0.091441 | 5.79 | 604.7 |
| 80bau3b.primal-c | - | 1.618 | - | 1.166 | 1.144355 | - | 7874.3 |
| primal_5959.primal-r | 0.067884 | 0.563 | 0.938 | 0.352 | 0.348951 | 5.14 | 2264.8 |
| degen3.dual-c | 0.082843 | 1.002 | - | 0.490 | 0.476030 | 5.75 | 2489.1 |
| degen3.dual-r | 0.082947 | 1.008 | 1.007 | 0.494 | 0.479651 | 5.78 | 2499.0 |
| perold.dual-r | 0.038532 | 0.355 | 0.360 | 0.151 | 0.146186 | 3.79 | 1322.1 |
| 80bau3b.dual-r | - | 1.103 | 1.688 | 0.415 | 0.399764 | - | 3737.4 |
| greenbea.dual-c | - | 2.870 | 2.363 | 0.838 | 0.825135 | - | 6885.9 |
| greenbea.dual-r | - | 3.598 | 3.592 | 0.814 | 0.800805 | - | 6826.6 |
| mode_switch.dual-r | - | 25.815 | - | 3.009 | 2.962669 | - | 19478.0 |
| primal_5959.dual-r | 0.098192 | 1.300 | - | 0.464 | 0.460424 | 4.69 | 3141.7 |
| quad_2717.dual-r | 0.092814 | 1.130 | - | 0.510 | 0.507528 | 5.47 | 3150.2 |
| fx_quad_obj_5371.dual-r | 0.038077 | 0.334 | 0.514 | 0.182 | 0.178987 | 4.70 | 1205.5 |
| degen3.barrier-r | 0.075214 | 0.608 | 0.856 | 0.461 | 0.401268 | 5.34 | 3780.7 |
| 80bau3b.barrier-r | - | 0.379 | 0.546 | 0.401 | 0.366603 | - | 3427.5 |
| greenbea.barrier-r | - | 1.569 | 2.493 | 1.497 | 1.447997 | - | 13029.7 |
| 80bau3b.pdhg | - | 4.552 | 6.526 | 6.058 | 6.035930 | - | 92388.9 |
| greenbea.pdhg | - | 2.801 | - | 1.898 | 1.885962 | - | 37740.6 |

## Hotspots

All simplex profiles except PILOT4 primal below predate row-indexed demotion; reprofile before selecting another simplex change.

Per profiled bench run: the commit it was taken at and the top self-time functions (percent), at most
five lines each. Re-profile a run when a landed change touched its hot path.

- degen2.primal-c, `rtyoxnwx` (stale after phase-one class reuse; frame-pointer build passes full 987/987 parity; forty solves capped at 2000 iterations, `target/perf-degen2c-next.data`): BTRAN 9.80%, primal control 7.95%, FTRAN 7.22%, update 4.98%, LU nucleus 4.26% self. Refactorization is 24.70%, FTRAN 13.22%, BTRAN 10.22% inclusive. All forty optimal in 531 iterations, 16 refactorizations. Signature construction is 1.26% and sig_class 1.11% self; repeated phase-one classification within the stable update window is now shared. Bounds-check removal and direct growth-diagonal access were measured and rejected below. Reprofile before choosing remaining BTRAN or refactor work.
- perold.dual-r, `muxlvwxk` (stale after row-indexed demotion; frame-pointer build passes full 987/987 parity; ten solves capped at 1300 iterations, `target/perf-perold-next6.data`): iteration 8.91%, FTRAN 7.39%, L stage 6.52%, update 6.08%, BTRAN 5.57% self. Refactorization is 24.54%, FTRAN 13.33%, ftran_needed 12.68%, L stage 10.31%, update 8.48% inclusive. All ten optimal in 1273 iterations, 27 refactorizations. Update annotation `target/perf-perold-next6-update.txt`; spike deduplication uses repeated linear membership searches. Prior rejected trials remain below.
- perold.primal-c, `yryslslw` (stale after growth-factor reuse and the dense work-count shortcut; frame-pointer build passes full 987/987 parity; three solves capped at 3500 iterations, `target/perf-perold.data`): BTRAN 8.71%, primal control 7.01%, row_product_work 6.39%, update 6.12%, fresh factor 6.10% self. Refactor 26.63%, FTRAN 11.82%, BTRAN 11.03% inclusive. All three optimal in 3399 iterations, 90 refactorizations; original reference also has 3399 iterations. The work-estimator annotation (`target/perf-perold-rowprod-annotate.txt`) attributes substantial time to checks in the touched-row scan; dense counts now avoid that scan. Reprofile before choosing the remaining sparse work.
- degen3.primal-c, `oytzqrpy` (stale after growth-factor reuse; frame-pointer build passes full 987/987 parity; three solves capped at 10000 iterations, `target/perf-degen3c-next.data`): primal control 8.91%, BTRAN 8.09%, FTRAN 6.51%, LU nucleus 5.82%, fresh factor 5.25% self. Refactorization is 31.03%, reference_order_for 19.77%, BTRAN 11.48%, FTRAN 10.88% inclusive; dense BTRAN output compaction is another 2.58% self. All three optimal in 2492 iterations, 49 refactorizations. The accepted growth-check factor now carries through to installation; reprofile before choosing remaining reference-order or BTRAN work.
- degen3.primal-r, `xutnmxks` (stale after dense BTRAN output reservation; frame-pointer build passes full 987/987 parity; three solves capped at 5300 iterations, `target/perf-degen3r-rest.data`): FTRAN 12.37%, BTRAN 7.48%, update 7.08%, LU nucleus 5.90%, primal control 4.73% self. Refactorization is 24.83%, FTRAN 21.75%, BTRAN 11.26% inclusive; dense BTRAN output compaction is 3.45% inclusive/2.88% self. All three optimal in 5165 iterations, 66 refactorizations. FTRAN annotation (`target/perf-degen3r-rest-ftran-annotate.txt`) highlights tail graph traversal, scatter bounds checks and output done/fallback scans; prior dense-FTRAN output simplification was rejected below. Dense BTRAN output now reserves its maximum row count before scanning; reprofile before choosing remaining FTRAN graph, scatter or output work.
- quad_2717.dual-r, `ussxlpst` (stale after native ranked row products; frame-pointer build passes full 987/987 parity; eight solves capped at 17000 iterations, `target/perf-quad-next.data`): F80 add 15.13%, primal control 10.38%, F80 multiply 8.89%, pack 5.06%, rank_products_q 4.22% self. Row products are 13.20%, cost_change_work 12.77%, refactorization 17.26% inclusive. All eight optimal in 16010 iterations, 3022 refactorizations. Ranked row products now use native extended arithmetic under the required control mode; the software kernel remains the fallback. Reprofile before choosing remaining arithmetic or set_d_from_y work. The rejected masks, coefficient cache, distant-exponent shortcut and forced add inlining are listed below.
- degen3.barrier-r, `munmovyt` (stale after kernel dispatch; frame-pointer build passes full parity; three solves capped at 100 iterations, `target/perf-barrier-degen3-valid.data`): FMA dispatch 26.46%, DGEMM 24.62%, existing hardware FMA implementation 14.77%, numeric_factor 7.40%, DTRSM 4.36% self. DGEMM 49.04%, numeric factor 81.45% inclusive; all three solves optimal in 14 iterations. The per-product FMA wrapper is an indirect call plus frame/return overhead (annotated `target/perf-barrier-fma-annotate.txt`); per-call dispatch now removes that per-product overhead without adding fused operations. BORE3D is a cheap witness (16 barrier iterations, 12 ms).
- primal_5959.dual-r, `tplnupps` (stale after factor reuse; frame-pointer build passes full parity; three solves capped at 16000 iterations, `target/perf-dual5959-full.data`): F80 add_signed 7.50%, primal control 5.70%, F80 multiply 4.54%, F80 nucleus 3.88%, double nucleus 3.74% self. All three reach numerical_failure, 15688 iterations, 7842 refactorizations. Refactorization is 51.10% inclusive, quad mirrors 21.69% before factor reuse. LUIN capture `target/perf-dual5959-luin.log` has 7822 records and 40 input-head/threshold combinations; two-entry LRU estimate 99.3% hits. Wide-shift arithmetic rewrites are rejected below. The earlier 11000-iteration capture stopped early and is not the full-run profile.
- pilot4.primal-r, `twspzytl` (frame-pointer build passes full 987/987 parity; twenty solves capped at 700 iterations, `target/perf-pilot-next5.data`): primal control 8.76%, BTRAN 7.88%, update 7.06%, FTRAN 5.53%, LU nucleus 3.93% self. Refactorization is 28.48%, numeric factor 3.57% inclusive. All twenty optimal in 690 iterations, 16 refactorizations. Numeric factor still performs an auxiliary pivot search before overwriting the result with supplied pivot rows. Prior BTRAN storage and buffer trials remain below.
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

- pilot4.primal-r (5.79x), degen3.dual-r (5.78x), degen3.dual-c (5.75x), degen2.dual-r (5.70x), perold.primal-c (5.65x): work in this order; reprofile changed paths first. Host speed bands can reorder these snapshot ratios.
- DEGEN3 dual: tau row-value lookup and immutable BTRAN metadata are settled. Reprofile before inspecting remaining iteration control, updates and scratch allocation.
- Other reachable runs follow their table ratios. Quad twins share construction and four compact repeated factors and double structures are cached; row-capacity growth already uses the numeric U transpose.
- Quad_2717: ranked extended row products now use native x87 arithmetic with a software fallback. Reprofile before inspecting remaining F80 arithmetic and set_d_from_y; the old call shares no longer describe the current solver. Cost-change work performs required reduced-cost updates; the factor-cache working set is settled.
- Barrier: profile again before choosing more work. DGEMM now dispatches existing fused operations once per call; other BLAS-like kernels still use their original scalar dispatch and accumulation order. A smaller remaining DGEMM gap may involve independent outputs and data reuse, which must preserve every per-output reduction order.
- Runs without reachable references, after the reachable ratios: mode-switch update (22.21%) still scans A/B lists for other memberships; BTRAN (13.33% before metadata caching) allocates scratch; LU nucleus (15.90%); fresh-factor row-list filtering (7.67%). These old profiles are stale; annotate again before choosing a change.
- all: dump::on/record construction without dumps and allocation in FTRAN/BTRAN, pricing and ratio tests remain unmeasured leads.

## Tried, no gain

- Skip auxiliary pivot selection and row-count setup when numeric LU receives prescribed rows (analysis idea 7), retaining max-magnitude/floor and supplied-pivot zero checks: fast parity 940/940 passes, PILOT4 instructions -0.11%, DEGEN2 primal -0.22%. Three alternating CPU-2 PILOT4 pairs give median solve medians 0.089099 -> 0.090044 s with no gain (aggregate spreads 0.007874/0.006202 s). Reverted (`target/perf-forcedpivot-{target,confirm}-{before,after}.txt`).

- Maintain B-list reference counts to replace factor-update full-tail membership scans (analysis idea 4): fast parity 940/940 passes, PEROLD instructions -0.59%, DEGEN2 dual -0.67%. Three alternating CPU-4 pairs give median solve medians 0.146994 -> 0.146334 s within aggregate 0.001341/0.003619 s spreads; CPU-2 rechecks cross speed bands without a gain beyond spread. Reverted (`target/perf-brefs-{target,confirm,recheck}-{before,after}.txt`). This tests B counts only; maintained reverse A memberships remain untested.

- Index nonbasic row-list membership (analysis idea 5), retaining linear first-match lookup for duplicate/handoff rows: fast parity 940/940 passes for both variants. Per-row hash maps increase PEROLD instructions 0.37%, DEGEN2 dual 1.55%; immutable matrix-entry IDs with array positions still increase PEROLD 0.05% and DEGEN2 0.29%. Array target PEROLD solve 0.232249 -> 0.229859 s is within 0.013778/0.010542 s spreads; DEGEN2 is slower. Both reverted (`target/perf-rowindex{,2}-target-{before,after}.txt`). The coefficient searches remain unchanged in these trials.

- Shared dense two-vector U-tail traversal with independent column/tau tolerance and scatter rules (analysis idea 1): frozen f64/F80 random/extreme/captured SC50A differential tests and fast parity 940/940 pass. Hoisting per-vector scatter tests out of the coefficient loop reduces PEROLD instructions 0.55%, but three alternating CPU-2 pairs give median solve medians 0.146587 -> 0.147595 s (aggregate spreads 0.001225/0.001968 s), with no timing gain; DEGEN2 target solve also increases. Reverted (`target/perf-pairtail{,3}-target-{before,after}.txt`, `target/perf-pairtail3-confirm-{before,after}.txt`). This trial shares only the dense U tail; shared L traversal and other paired paths remain untested.

- Cache the immutable L-level maximum for FTRAN (analysis idea 8): fast parity 940/940 passes and instructions fall 0.64% on PEROLD, 1.15% on DEGEN2 dual. PEROLD initial solve 0.148098 -> 0.146727 s is within the 0.001604 s baseline spread; three CPU-4 confirmation pairs give median solve medians 0.148094 -> 0.148166 s. CPU-2 rechecks cross speed bands and show no gain beyond aggregate spreads. Reverted (`target/perf-maxlev-{target,confirm,recheck}-{before,after}.txt`). Timing evidence is inconclusive; the instruction reduction alone does not meet acceptance.

- Marker array for factor-update spike deduplication instead of repeated `contains`: fast parity 940/940 passes, PEROLD instructions -0.23%, but three alternating CPU-2 pairs give median solve medians 0.230458 -> 0.231593 s with no gain (aggregate spreads 0.087231/0.091101 s); final pair also has no improvement. Reverted (`target/perf-spikemark-{target,confirm}-{before,after}.txt`).

- Share sparse eta accounting with the numeric pass only when `first == 0` (analysis idea 2): frozen f64/F80 random/extreme/captured-input tests and fast parity 940/940 pass, but PEROLD instructions fall only 0.11% and solve 0.232151 -> 0.234504 s has no gain (spreads 0.017993/0.004013 s); DEGEN2 instructions are flat. Reverted (`target/perf-etawork-target-{before,after}.txt`). General prefix sharing remains untested: the accounting input omits prior eta updates and may differ on sub-tolerance residues.

- Reuse dual pivot-row accumulator and marker buffers with full clearing (analysis idea 6): fast parity 940/940 passes; initial PEROLD instructions -0.47%, but three alternating CPU-2 pairs have effectively flat instructions and no timing gain (median of solve medians 0.236570 -> 0.241332 s; aggregate spreads 0.012539/0.082480 s). Reverted (`target/perf-rowworkspace-{target,confirm}-{before,after}.txt`).

- Column-local dual row-product sums for dense, strictly ascending rho and unique matrix rows: frozen random/extreme/captured ADLITTLE f64/F80 tests and fast parity 940/940 pass, but PEROLD instructions +1.8% with separate presence/value arrays, +2.2% with an Option buffer. The latter regresses solve 0.236082 -> 0.240998 s beyond 0.003144/0.002701 s spreads. DEGEN3 instructions also increase; both variants reverted (`target/perf-colrow{,2}-target-{before,after}.txt`).

- Reserve FTRAN result capacity from the output-order length: fast parity 940/940 passes and instructions fall 0.4–0.6%, but the initial PEROLD timing gain crosses host speed bands; alternating confirmation 0.151656 -> 0.149696 s is within 0.093843/0.090559 s spreads, and DEGEN2 has no repeatable gain. Reverted (`target/perf-ftreserve-{target,recheck}-{before,after}.txt`).

- Borrow pre-sorted columns in dual true_d instead of cloning and sorting: fast parity 940/940 passes and instructions fall 0.4–0.6%, but alternating CPU-2 PEROLD solve 0.231300 -> 0.231688 s has no gain (spreads 0.017962/0.008958 s); DEGEN2 and DEGEN3 stay within wide spreads. Reverted (`target/perf-truecols-{target,recheck}-{before,after}.txt`).

- Split dense dual pivot-row marking from sparse first-touch ordering: fast parity 940/940 passes, PEROLD instructions -0.32%, but alternating CPU-2 solve 0.232569 -> 0.227484 s is within 0.009078/0.006368 s spreads; DEGEN2 has no repeatable timing gain. Reverted (`target/perf-rowmark-{target,recheck}-{before,after}.txt`).

- Separate immutable L-transpose index and coefficient arrays with contiguous row offsets: fast parity 940/940 passes, but PILOT4 instructions +2.7%, DEGEN2 +3.5%, PEROLD primal +4.5%, and all timing changes stay within spreads. Reverted (`target/perf-ltsoa-target-{before,after}.txt`).

- Four-entry ordered scalar dense BTRAN scatter (reference U-transpose 0x12e25e, L-transpose 0x130197): frozen random/extreme/duplicate/captured SC50A f64/F80 tests and fast parity 940/940 pass, but PILOT4 instructions +0.3%, DEGEN2 +0.8%, PEROLD primal +0.6%; no timing gain beyond spreads. Reverted (`target/perf-btunroll-target-{before,after}.txt`).

- Guard the sparse L rollback clone by m >= 250 and Kahn reached-count upper bound: fast parity 940/940 passes, but instructions are effectively flat (PILOT4 -0.01%, Quad_2717 +0.01%); DEGEN2 0.068212 -> 0.064283 s is within the candidate 0.004463 s spread. Reverted (`target/perf-lclone-target-{before,after}.txt`).

- Transfer dense L-stage x directly into U with zero canonicalization: frozen f64/F80 U tests and fast parity 940/940 pass. Combined collection/normalization lowers DEGEN3 dual instructions 2.1%, but CPU-2 solve 0.486473 -> 0.493595 s has no gain (spreads 0.012765/0.036416 s). DEGEN2 initial -2.6% does not repeat beyond spread; PILOT4 instructions rise 0.2%. A separate zero-normalization pass was worse. Reverted (`target/perf-lbuffer{,2}-target-{before,after}.txt`, `target/perf-lbuffer2-recheck-{before,after}.txt`).

- Pack immutable L-transpose rows into contiguous storage: fast parity 940/940 passes, but PILOT4 solve 0.095785 -> 0.096071 s is within 0.002578/0.006817 s spreads and instructions rise 0.3%; DEGEN2 instructions rise 0.15% with no timing gain. Reverted (`target/perf-ltcsr-target-{before,after}.txt`).

- Borrow Markowitz candidate chains and rotate buckets in place: fast parity 940/940 passes and PEROLD instructions -0.5%, but alternating CPU-2 solve 0.237249 -> 0.236897 s is within 0.006498/0.011017 s spreads; DEGEN3 primal also stays within wide spreads. Reverted (`target/perf-searchborrow-{target,recheck}-{before,after}.txt`).

- Cache first-match pivot-row/column coefficients in a stamped lookup: semantic tests and fast parity 940/940 pass, PEROLD instructions -0.6%, but solve 0.151853 -> 0.151628 s is within 0.004297/0.000431 s spreads; DEGEN2 has a small timing increase with flat instructions and DEGEN3 primal has no gain. Reverted (`target/perf-pivotvalues-target-{before,after}.txt`).

- Pair eight-element eta chunks to eliminate per-lane checks: frozen f64/F80 tests and fast parity 940/940 pass, dot code shrinks 1037 -> 879 bytes and PEROLD instructions -1.5%, but CPU-4 solve medians are 0.150879/0.150880 s; alternating CPU-2 changes stay inside wide spreads. Reverted (`target/perf-dotchunks-{target,recheck}-{before,after}.txt`).

- Explicit eight-lane eta blocks: frozen random/extreme/captured tests, full parity 987/987 and release tests pass; PEROLD instructions -0.9%, but the full CPU-4 solve 0.151302 -> 0.151144 s stays within 0.001056/0.000552 s spreads, and confirmation has no gain. PILOT4 confirmation regresses 0.094533 -> 0.097417 s (+3.1%; spreads 0.002414/0.001277 s). Reverted (`target/perf-dotblocks-{target,recheck,full,confirm}-{before,after}.txt`).

- Read eta inputs directly by index instead of gathering: frozen random/extreme/captured f64/F80 tests and fast parity 940/940 pass, PEROLD instructions -1.0%, but alternating CPU-2 solve 0.230032 -> 0.232767 s has no gain (spreads 0.006889/0.003601 s); DEGEN3 also stays within spread. Reverted (`target/perf-indexdot-{target,recheck}-{before,after}.txt`).

- Elide pivot-row scatter bounds checks using private row-list indices: fast parity 940/940 passes and PEROLD instructions -0.9%, but alternating CPU-2 solve 0.234652 -> 0.232350 s is within 0.002445/0.004874 s spreads; DEGEN3 has no repeatable timing gain. Reverted (`target/perf-rowindex-{target,recheck}-{before,after}.txt`).

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
