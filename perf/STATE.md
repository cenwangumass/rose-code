# Performance state

Current state of the work in [PLAN.md](PLAN.md). Rewrite in place; no log (history is `jj log`).

## Numbers

Baseline solver: `8b645007` (same solver as `vlsokskq`). Current solver: `ppxxsslu` (PGO builder with rustc-default online profile merging; solver source unchanged from `wvqvkoyp`).
Parity: full tier 987/987; release tests pass. Three timed runs after one warm-up, hyperfine 2.0.0; instruction counts are millions from perf stat.

`baseline wall` is the fixed initial process-wall median. `paired base wall` is the fresh baseline-binary median from alternating baseline/last-landed/candidate checks on one CPU (CPU 0 unless stated otherwise); it resolves host drift without changing the fixed baseline. `current wall` and `solve` use the new harness; **current/ref uses solve-only time**, `diagnostics.seconds`, against original Gurobi Runtime.

Reference provenance: original-library cached records `../reverse-engineer/runs/{degen2,degen3,perold,pilot4}-{profile}.json` for the ten previously valid cells. New three-run medians use `/tmp/rose-perf-reference/{primal_5959.primal-r,primal_5959.dual-r,quad_2717.dual-r,fx_quad_obj_5371.dual-r}.{1,2,3}.json`; summary `/tmp/rose-perf-reference/summary.json`. All original-library records identify Ryzen 9 5900X and library SHA256 `8e87fca221fdffc922331c8619426031b3e72dfd419446422f88af0cb0e98f9c`.
`-` refs are unreachable: original Gurobi rejects 80BAU3B, GREENBEA and mode-switch with license error 10010 (scratch `ID.error.json`, including direct 80BAU3B primal/dual attempts). The old 80BAU3B primal/dual refs used the patched library and have been removed; the old prose mode-switch time is unusable.

Regression checks: all 23 runs paired against the validated PGO CLI on CPU 2, alternating pair order. Primal_5959 primal and mode-switch dual timing flags have mixed CPU-4 follow-ups across host speed bands, with no repeatable regression; comparable pairs stay within spread and instructions are effectively unchanged. Table uses median follow-up medians for those two rows and full-sweep candidate times elsewhere. Do not credit host-band timing shifts to this build fix. Raw data `target/perf-pgo-online-full-{before,after}.txt`, `target/perf-pgo-online-followup-{before,after}.txt`.
Acceptance: the PGO builder now clears `LLVM_PROFILE_FILE` rather than forcing per-process raw files, retaining rustc's validated default shared online merge. Same-source controlled builds, same 39 training cases: separate files increase 80BAU3B barrier instructions 3.045B -> 3.111B (+2.17%) and solve 0.329723 -> 0.337777 s (+2.44%; spreads 0.003364/0.003766 s); GREENBEA instructions 11.489B -> 11.683B (+1.69%). Restoring default collection returns barrier instruction counts to the original validated baseline. The exact downstream optimization site is not yet isolated. Raw control data `target/perf-pgo-collection-{online,separate}.txt`; original-baseline comparison `target/perf-pgo-online-{before,after}.txt`. This fixes build reproducibility rather than changing a numeric kernel.
Build: `python3 perf/build_pgo.py --check --test`; compiler/flags/training manifest in `target/perf-pgo/build.json`. Tested profile SHA256 `84fab3f7183cf79ad66ef39342ddb172a5ad72775d4a42307a8768682e43b60c`; build `/tmp/rose-pgo-onlinemerge-build.txt`, full parity and tests `/tmp/rose-pgo-online-checks.txt`. Full parity 987/987 and release tests pass. Solver source and frozen parity baseline are unchanged. Preserve the default profile collection mode when retraining; profile-use metadata is needed for profiling this optimized CLI.

| run | ref solve | baseline wall | paired base wall | current wall | solve | current/ref | instr |
|---|---|---|---|---|---|---|---|
| degen2.primal-c | 0.012711 | 0.182 | - | 0.068 | 0.063352 | 4.98 | 319.4 |
| degen2.dual-r | 0.009653 | 0.125 | - | 0.054 | 0.049241 | 5.10 | 252.9 |
| degen3.primal-c | 0.163137 | 3.237 | - | 0.849 | 0.837103 | 5.13 | 4576.9 |
| degen3.primal-r | 0.193685 | 3.872 | - | 1.007 | 0.994089 | 5.13 | 5273.5 |
| perold.primal-c | 0.160994 | 2.430 | - | 0.886 | 0.879085 | 5.46 | 5048.5 |
| pilot4.primal-r | 0.015787 | 0.253 | - | 0.094 | 0.088366 | 5.60 | 571.6 |
| 80bau3b.primal-c | - | 1.618 | - | 1.126 | 1.106779 | - | 7855.8 |
| primal_5959.primal-r | 0.067884 | 0.563 | 0.938 | 0.280 | 0.277395 | 4.09 | 2087.4 |
| degen3.dual-c | 0.082843 | 1.002 | - | 0.474 | 0.461698 | 5.57 | 2338.2 |
| degen3.dual-r | 0.082947 | 1.008 | 1.007 | 0.473 | 0.459940 | 5.54 | 2347.2 |
| perold.dual-r | 0.038532 | 0.355 | 0.360 | 0.259 | 0.252321 | 6.55 | 1244.8 |
| 80bau3b.dual-r | - | 1.103 | 1.688 | 0.399 | 0.385130 | - | 3632.4 |
| greenbea.dual-c | - | 2.870 | 2.363 | 0.781 | 0.769562 | - | 6544.4 |
| greenbea.dual-r | - | 3.598 | 3.592 | 0.776 | 0.764257 | - | 6484.2 |
| mode_switch.dual-r | - | 25.815 | - | 3.147 | 3.100158 | - | 20225.2 |
| primal_5959.dual-r | 0.098192 | 1.300 | - | 0.393 | 0.390545 | 3.98 | 2877.2 |
| quad_2717.dual-r | 0.092814 | 1.130 | - | 0.454 | 0.451586 | 4.87 | 2912.0 |
| fx_quad_obj_5371.dual-r | 0.038077 | 0.334 | 0.514 | 0.158 | 0.154882 | 4.07 | 1086.1 |
| degen3.barrier-r | 0.075214 | 0.608 | 0.856 | 0.645 | 0.580414 | 7.72 | 3648.8 |
| 80bau3b.barrier-r | - | 0.379 | 0.546 | 0.367 | 0.337943 | - | 3046.0 |
| greenbea.barrier-r | - | 1.569 | 2.493 | 1.322 | 1.277033 | - | 11488.1 |
| 80bau3b.pdhg | - | 4.552 | 6.526 | 4.962 | 4.944283 | - | 74948.7 |
| greenbea.pdhg | - | 2.801 | - | 1.494 | 1.484601 | - | 31955.5 |

## Hotspots

Except for DEGEN3 barrier below, profiles predate the corrected PGO builder; use the current profile-use flag and reprofile before selecting another change.

Per profiled bench run: the commit it was taken at and the top self-time functions (percent), at most
five lines each. Re-profile a run when a landed change touched its hot path.

- degen2.primal-c, `rtyoxnwx` (stale after phase-one class reuse; frame-pointer build passes full 987/987 parity; forty solves capped at 2000 iterations, `target/perf-degen2c-next.data`): BTRAN 9.80%, primal control 7.95%, FTRAN 7.22%, update 4.98%, LU nucleus 4.26% self. Refactorization is 24.70%, FTRAN 13.22%, BTRAN 10.22% inclusive. All forty optimal in 531 iterations, 16 refactorizations. Signature construction is 1.26% and sig_class 1.11% self; repeated phase-one classification within the stable update window is now shared. Bounds-check removal and direct growth-diagonal access were measured and rejected below. Reprofile before choosing remaining BTRAN or refactor work.
- perold.dual-r, `muxlvwxk` (stale after row-indexed demotion; frame-pointer build passes full 987/987 parity; ten solves capped at 1300 iterations, `target/perf-perold-next6.data`): iteration 8.91%, FTRAN 7.39%, L stage 6.52%, update 6.08%, BTRAN 5.57% self. Refactorization is 24.54%, FTRAN 13.33%, ftran_needed 12.68%, L stage 10.31%, update 8.48% inclusive. All ten optimal in 1273 iterations, 27 refactorizations. Update annotation `target/perf-perold-next6-update.txt`; spike deduplication uses repeated linear membership searches. Prior rejected trials remain below.
- perold.primal-c, `wvqvkoyp` (stale after PGO collection-mode correction; PGO frame-pointer build passes full 987/987 parity; three solves capped at 3500 iterations, `target/perf-peroldc-pgo.data`): primal control 16.89%, BTRAN 9.75%, LU nucleus 4.81%, FTRAN 4.78%, score 4.48% self. Refactorization is 23.52%, BTRAN 9.99%, FTRAN 9.12%, row-product work 1.91% inclusive. All three optimal in 3399 iterations, 90 refactorizations. Score repeatedly checks immutable free-variable bounds; prior row-product and solve trials remain below.
- degen3.primal-c, `oytzqrpy` (stale after growth-factor reuse; frame-pointer build passes full 987/987 parity; three solves capped at 10000 iterations, `target/perf-degen3c-next.data`): primal control 8.91%, BTRAN 8.09%, FTRAN 6.51%, LU nucleus 5.82%, fresh factor 5.25% self. Refactorization is 31.03%, reference_order_for 19.77%, BTRAN 11.48%, FTRAN 10.88% inclusive; dense BTRAN output compaction is another 2.58% self. All three optimal in 2492 iterations, 49 refactorizations. The accepted growth-check factor now carries through to installation; reprofile before choosing remaining reference-order or BTRAN work.
- degen3.primal-r, `votqkmtm` (stale after PGO; ThinLTO frame-pointer build passes full 987/987 parity; three solves capped at 5300 iterations, `target/perf-degen3r-thin.data`): BTRAN 10.62%, FTRAN 10.16%, update 6.43%, primal control 5.46%, factor_structure 4.37% self. Refactorization is 24.75%, FTRAN 18.32%, BTRAN 11.12% inclusive. All three optimal in 5165 iterations, 66 refactorizations. Signature construction is 2.60% self and repeatedly recomputes immutable bound flags and widths. Prior output/solve micro-optimization trials remain below.
- quad_2717.dual-r, `ussxlpst` (stale after native ranked row products; frame-pointer build passes full 987/987 parity; eight solves capped at 17000 iterations, `target/perf-quad-next.data`): F80 add 15.13%, primal control 10.38%, F80 multiply 8.89%, pack 5.06%, rank_products_q 4.22% self. Row products are 13.20%, cost_change_work 12.77%, refactorization 17.26% inclusive. All eight optimal in 16010 iterations, 3022 refactorizations. Ranked row products now use native extended arithmetic under the required control mode; the software kernel remains the fallback. Reprofile before choosing remaining arithmetic or set_d_from_y work. The rejected masks, coefficient cache, distant-exponent shortcut and forced add inlining are listed below.
- degen3.barrier-r, `ppxxsslu` (PGO profile-use/frame-pointer build, full parity 987/987; three solves capped at 100 iterations, all optimal in 14 iterations; `target/perf-barrier-pgo.data`): DGEMM existing FMA 23.15%, numeric_factor 13.11%, DTRSM 8.99%, scalar FMA 6.36%, ordering 5.46% self. BORE3D remains a cheap witness (16 barrier iterations, about 12 ms).
- primal_5959.dual-r, `tplnupps` (stale after factor reuse; frame-pointer build passes full parity; three solves capped at 16000 iterations, `target/perf-dual5959-full.data`): F80 add_signed 7.50%, primal control 5.70%, F80 multiply 4.54%, F80 nucleus 3.88%, double nucleus 3.74% self. All three reach numerical_failure, 15688 iterations, 7842 refactorizations. Refactorization is 51.10% inclusive, quad mirrors 21.69% before factor reuse. LUIN capture `target/perf-dual5959-luin.log` has 7822 records and 40 input-head/threshold combinations; two-entry LRU estimate 99.3% hits. Wide-shift arithmetic rewrites are rejected below. The earlier 11000-iteration capture stopped early and is not the full-run profile.
- pilot4.primal-r, `twspzytl` (stale after ThinLTO; frame-pointer build passes full 987/987 parity; twenty solves capped at 700 iterations, `target/perf-pilot-next5.data`): primal control 8.76%, BTRAN 7.88%, update 7.06%, FTRAN 5.53%, LU nucleus 3.93% self. Refactorization is 28.48%, numeric factor 3.57% inclusive. All twenty optimal in 690 iterations, 16 refactorizations. Numeric factor still performs an auxiliary pivot search before overwriting the result with supplied pivot rows. Prior BTRAN storage and buffer trials remain below.
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

- degen3.barrier-r (7.72x), perold.dual-r (6.55x), pilot4.primal-r (5.60x), degen3.dual-c (5.57x), degen3.dual-r (5.54x): work in this order; reprofile changed paths first. Host speed bands can reorder these snapshot ratios.
- Unapplied primal free-variable flag cache: full parity/tests pass and local PEROLD pairs improve, but its barrier comparison was confounded by per-process PGO collection versus the baseline's online merge. Source performance remains unestablished; retest with consistent collection (`target/perf-scorefree-*`).
- Unapplied primal row-product prefilter: skip dropped precomputed products before eligibility checks, preserving order and NaNs. Full parity/tests pass; local PEROLD instructions -1.4% with stable pair improvements. Its full-sweep comparison also used inconsistent PGO collection, so retain it as an unverified source lead (`target/perf-rowprefilter-*`, candidate copy `/tmp/rose-rowprefilter-candidate.rs`).
- DEGEN3 dual: tau row-value lookup and immutable BTRAN metadata are settled. Reprofile before inspecting remaining iteration control, updates and scratch allocation.
- Other reachable runs follow their table ratios. Quad twins share construction and four compact repeated factors and double structures are cached; row-capacity growth already uses the numeric U transpose.
- Quad_2717: ranked extended row products now use native x87 arithmetic with a software fallback. Reprofile before inspecting remaining F80 arithmetic and set_d_from_y; the old call shares no longer describe the current solver. Cost-change work performs required reduced-cost updates; the factor-cache working set is settled.
- Barrier: current PGO profile above identifies numeric factorization, DTRSM and remaining scalar FMA dispatch costs. DGEMM annotation (`target/perf-barrier-pgo-dgemm-annotate.txt`) concentrates samples in the nocopy `Kind::S` serial FMA chain. Test independent-output tiling (analysis idea 10), retaining each chain and block boundary; frozen DGEMM random/extreme/captured-input tests already cover this path. Diagnostic preparation and child-buffer trials are rejected below. Other BLAS-like kernels still use scalar dispatch; DTRSM dispatch failed the extreme-input invariant.
- Runs without reachable references, after the reachable ratios: mode-switch update (22.21%) still scans A/B lists for other memberships; BTRAN (13.33% before metadata caching) allocates scratch; LU nucleus (15.90%); fresh-factor row-list filtering (7.67%). These old profiles are stale; annotate again before choosing a change.
- all: dump::on/record construction without dumps and allocation in FTRAN/BTRAN, pricing and ratio tests remain unmeasured leads.

## Tried, no gain

- Borrow contiguous child panels from packed factor storage instead of allocating/copying them before DGEMM (analysis idea 13): fast parity 940/940 passes and DEGEN3 barrier instructions fall about 0.9%. Three alternating CPU-4 pairs cross speed bands; CPU-2 confirmation median solve medians 0.405972 -> 0.396088 s remain within aggregate 0.096096/0.022184 s spreads. The final stable baseline/candidate pair is 0.397590 -> 0.396078 s, within 0.003082/0.006880 s spreads. BORE3D also has no gain beyond spread. Reverted (`target/perf-childpanel-*`); instruction savings alone do not meet acceptance.

- Reuse child-update temporary/panel/product buffers within each factorization (analysis idea 13): fast parity 940/940 passes; DEGEN3 barrier solve 0.393035 -> 0.387227 s is within 0.010109/0.006277 s spreads and instructions rise about 0.16%. BORE3D timing is also within spread with effectively unchanged instructions. Reverted (`target/perf-childscratch-*`).

- Skip unused front-dump configuration parsing (analysis idea 16): fast parity 940/940 and bounded BORE3D front-dump comparisons pass. Hoisting selection/entries to factorization entry raises DEGEN3 barrier instructions 0.60%; guarding only entry parsing for the selected supernode raises them about 0.58%. Three alternating CPU-2 pairs for the latter give median solve medians 0.398125 -> 0.405277 s across host speed bands; neither meets acceptance. Both reverted (`target/perf-frontcfg{,2}-*`).

- Combine FTRAN touched/emitted markers into one three-state byte array: frozen f64/F80 random/extreme/captured-input tests and fast parity 940/940 pass. DEGEN3 primal-release solve 0.651511 -> 0.646340 s improves 0.79% beyond 0.004325/0.004931 s spreads, but instructions rise 0.17% (DEGEN2 primal +0.08%), failing the below-3% instruction-reduction acceptance rule. Reverted (`target/perf-ftmarks-target-{before,after}.txt`). FTRAN annotation `target/perf-degen3r-thin-ftran.txt` identifies output marker checks, scatter indexing and divisions as remaining costs.

- Avoid materializing the dense FTRAN tail-index vector: direct slice traversal through a shared solve closure passes frozen f64/F80 random/extreme/captured-input tests and fast parity 940/940, but DEGEN3 primal-release instructions rise 1.29% and solve 0.648145 -> 0.656844 s exceeds both 0.001714/0.002821 s spreads. A single chained iterator passes fast parity but still adds 0.25% instructions and has no timing gain (1.019540 -> 1.024090 s within spread). Both reverted (`target/perf-tailview{,2}-target-{before,after}.txt`).

- Cache primal signature bound flags/widths (analysis idea 8): fast parity 940/940 passes. Interleaved `(kind,width)` cache reduces DEGEN3 primal-release instructions 1.76%, but three alternating CPU-2 pairs give median solve medians 1.005068 -> 1.001979 s within 0.042525/0.061157 s aggregate spreads. Split byte flags/width arrays reduce instructions 2.54%, but confirmation medians 1.003885 -> 1.313467 s show no gain and cross speed bands. Both reverted (`target/perf-sigbounds{,2}-{target,confirm}-{before,after}.txt`).

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

(idea, first differing record or frozen-kernel differential result, why)

- Dispatch existing DTRSM fused operations once per kernel/panel (analysis idea 11): rejected before corpus parity. Frozen-kernel random/extreme differential test first differs at width 12, one RHS, column 0: the final output is negative versus positive quiet NaN (`0xfff8000000000000` versus `0x7ff8000000000000`). Finite input guards are insufficient because intermediate overflow can produce non-finite values. Source and test scaffolding reverted; no performance claim. Failure evidence `/tmp/rose-trsmfma-test.txt`.
