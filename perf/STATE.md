# Performance state

Current state of the work in [PLAN.md](PLAN.md). Rewrite in place; no log (history is `jj log`).

## Numbers

Baseline solver: `8b645007` (same solver as `vlsokskq`). Current solver: `sormxnyt` (generation-stamped numeric U row lookup for tau dots).
Parity: full tier 987/987; release tests pass. Three timed runs after one warm-up, hyperfine 2.0.0; instruction counts are millions from perf stat.

`baseline wall` is the fixed initial process-wall median. `paired base wall` is the fresh baseline-binary median from alternating baseline/last-landed/candidate checks on one CPU (CPU 0 unless stated otherwise); it resolves host drift without changing the fixed baseline. `current wall` and `solve` use the new harness; **current/ref uses solve-only time**, `diagnostics.seconds`, against original Gurobi Runtime.

Reference provenance: original-library cached records `../reverse-engineer/runs/{degen2,degen3,perold,pilot4}-{profile}.json` for the ten previously valid cells. New three-run medians use `/tmp/rose-perf-reference/{primal_5959.primal-r,primal_5959.dual-r,quad_2717.dual-r,fx_quad_obj_5371.dual-r}.{1,2,3}.json`; summary `/tmp/rose-perf-reference/summary.json`. All original-library records identify Ryzen 9 5900X and library SHA256 `8e87fca221fdffc922331c8619426031b3e72dfd419446422f88af0cb0e98f9c`.
`-` refs are unreachable: original Gurobi rejects 80BAU3B, GREENBEA and mode-switch with license error 10010 (scratch `ID.error.json`, including direct 80BAU3B primal/dual attempts). The old 80BAU3B primal/dual refs used the patched library and have been removed; the old prose mode-switch time is unusable.

Regression checks: all 23 runs paired against the preceding solver on CPU 4. Six apparent increases resolve in CPU-2 follow-ups: DEGEN3 primal-controlled 0.880814/0.877701 s; primal_5959 primal 0.419684/0.424432 s (old spread 0.143651 s); DEGEN3 dual-controlled 0.725833/0.696667 s; primal_5959 dual 0.618439/0.619766 s (new spread 0.016885 s); quad_2717 0.700469/0.714816 s (new spread 0.017975 s). Mode-switch changes host-speed bands in the opposite direction: 11.102679/7.429264 s, with matching 76027.7/76029.1M instructions; do not attribute that timing shift to this change. No new fixed-baseline regression; 80BAU3B PDHG remains within its previously paired baseline range, with unchanged instructions. Raw data `target/perf-taurow-full-{before,after}.txt`, `target/perf-taurow-recheck-{before,after}.txt`, `target/perf-taurow-mode-{before,after}.txt`. Host speed varies substantially; use immediate pairs and spreads.
Acceptance: independent CPU-2 DEGEN3 dual-release solve 0.734087 -> 0.698014 s (-4.9%; spreads 0.007071/0.002878 s), instructions 3769.0 -> 3690.3M (-2.1%). The full sweep also falls 0.722820 -> 0.703693 s, though that delta is inside its wider new spread. DEGEN2 dual instructions fall slightly; PEROLD rises 0.15% with timing within spread. Tau dots now copy a numeric U row into reusable generation-stamped scratch, then gather terms in the unchanged virtual rowl order. First-duplicate selection, dot-product operation order and WORK estimates remain exact; tests compare the direct lookup with independent old-column scans for f64/F80, including generation wrap.

| run | ref solve | baseline wall | paired base wall | current wall | solve | current/ref | instr |
|---|---|---|---|---|---|---|---|
| degen2.primal-c | 0.012711 | 0.182 | - | 0.094 | 0.089170 | 7.02 | 481.2 |
| degen2.dual-r | 0.009653 | 0.125 | - | 0.077 | 0.072043 | 7.46 | 357.9 |
| degen3.primal-c | 0.163137 | 3.237 | - | 0.888 | 0.877701 | 5.38 | 7652.4 |
| degen3.primal-r | 0.193685 | 3.872 | - | 1.048 | 1.036732 | 5.35 | 9400.6 |
| perold.primal-c | 0.160994 | 2.430 | - | 0.858 | 0.853453 | 5.30 | 7535.2 |
| pilot4.primal-r | 0.015787 | 0.253 | - | 0.088 | 0.084476 | 5.35 | 820.5 |
| 80bau3b.primal-c | - | 1.618 | - | 0.942 | 0.927054 | - | 10091.0 |
| primal_5959.primal-r | 0.067884 | 0.563 | 0.938 | 0.427 | 0.424432 | 6.25 | 2863.5 |
| degen3.dual-c | 0.082843 | 1.002 | - | 0.712 | 0.696667 | 8.41 | 3679.7 |
| degen3.dual-r | 0.082947 | 1.008 | 1.007 | 0.718 | 0.703693 | 8.48 | 3689.6 |
| perold.dual-r | 0.038532 | 0.355 | 0.360 | 0.328 | 0.320998 | 8.33 | 1822.1 |
| 80bau3b.dual-r | - | 1.103 | 1.688 | 0.854 | 0.831567 | - | 5152.8 |
| greenbea.dual-c | - | 2.870 | 2.363 | 2.020 | 1.999138 | - | 10822.0 |
| greenbea.dual-r | - | 3.598 | 3.592 | 2.033 | 2.011256 | - | 10722.5 |
| mode_switch.dual-r | - | 25.815 | - | 7.576 | 7.429264 | - | 76029.1 |
| primal_5959.dual-r | 0.098192 | 1.300 | - | 0.623 | 0.619766 | 6.31 | 4011.6 |
| quad_2717.dual-r | 0.092814 | 1.130 | - | 0.718 | 0.714816 | 7.70 | 4202.1 |
| fx_quad_obj_5371.dual-r | 0.038077 | 0.334 | 0.514 | 0.165 | 0.163290 | 4.29 | 1621.0 |
| degen3.barrier-r | 0.075214 | 0.608 | 0.856 | 0.317 | 0.276776 | 3.68 | 3780.7 |
| 80bau3b.barrier-r | - | 0.379 | 0.546 | 0.272 | 0.248529 | - | 3428.0 |
| greenbea.barrier-r | - | 1.569 | 2.493 | 1.013 | 0.979056 | - | 13026.9 |
| 80bau3b.pdhg | - | 4.552 | 6.526 | 6.515 | 6.494430 | - | 92390.6 |
| greenbea.pdhg | - | 2.801 | - | 1.910 | 1.898216 | - | 37740.7 |

## Hotspots

Per profiled bench run: the commit it was taken at and the top self-time functions (percent), at most
five lines each. Re-profile a run when a landed change touched its hot path.

- quad_2717.dual-r, `zzvmznmp` (stale after cache sizing; frame-pointer build passes full parity; three solves capped at 20000 iterations, `target/perf-quad2717.data`): F80 add 13.58%, primal control 8.94%, F80 multiply 8.17%, pack 3.80%, rank_products_q 3.61% self. Refactorization 25.94%, quad mirrors 11.50%, rank_products_q 11.26%, cost_change_work 11.18% inclusive. All solves optimal in 16010 iterations with 3022 refactorizations. `target/perf-quad2717-luin.log`: 3015 LUIN records, 14 distinct input-head/threshold keys; FIFO hits 110/2998/3000 at capacities 2/4/8. The cache now holds four entries; larger capacities add only marginal hits. Cost-change work also updates reduced costs and cannot simply be removed.
- degen3.barrier-r, `munmovyt` (stale after kernel dispatch; frame-pointer build passes full parity; three solves capped at 100 iterations, `target/perf-barrier-degen3-valid.data`): FMA dispatch 26.46%, DGEMM 24.62%, existing hardware FMA implementation 14.77%, numeric_factor 7.40%, DTRSM 4.36% self. DGEMM 49.04%, numeric factor 81.45% inclusive; all three solves optimal in 14 iterations. The per-product FMA wrapper is an indirect call plus frame/return overhead (annotated `target/perf-barrier-fma-annotate.txt`); per-call dispatch now removes that per-product overhead without adding fused operations. BORE3D is a cheap witness (16 barrier iterations, 12 ms).
- primal_5959.dual-r, `tplnupps` (stale after factor reuse; frame-pointer build passes full parity; three solves capped at 16000 iterations, `target/perf-dual5959-full.data`): F80 add_signed 7.50%, primal control 5.70%, F80 multiply 4.54%, F80 nucleus 3.88%, double nucleus 3.74% self. All three reach numerical_failure, 15688 iterations, 7842 refactorizations. Refactorization is 51.10% inclusive, quad mirrors 21.69% before factor reuse. LUIN capture `target/perf-dual5959-luin.log` has 7822 records and 40 input-head/threshold combinations; two-entry LRU estimate 99.3% hits. Wide-shift arithmetic rewrites are rejected below. The earlier 11000-iteration capture stopped early and is not the full-run profile.
- pilot4.primal-r, `uprrtwlk` (stale after row-capacity count reuse; frame-pointer build passes full parity; twenty solves capped at 700 iterations, `target/perf-pilot4-fp.data`): update 10.42%, BTRAN 7.84%, fresh factor 7.05%, FTRAN 6.89%, primal control 6.83% self. Refactor 28.11%, fresh factor 10.90% inclusive. Every solve is optimal, 690 iterations, 16 refactorizations. Update annotation: numeric U-entry recount about 18% of update time, structural rcount scan also hot (the rowl shortcut was already rejected below); do not conflate the two counts.
- primal_5959.primal-r, `kyyvxyql` (stale after twin reuse; frame-pointer build passes full 987/987 parity; three bounded solves at 1999 Hz, `target/perf-fp-primal5959.data`, report with `--no-inline`): primal control 6.32%, F80 nucleus 5.46%, fresh F80 factor 5.09%, F80 add_signed 6.67% self. Refactor 53.53%, quad mirrors 28.88% inclusive before twin reuse. Log verifies numerical_failure, 10780 iterations, 5381 refactorizations; original reference status 12, 10780 iterations, refactor count unavailable. Default inline expansion duplicates frames; use PLAN.md's validated frame-pointer workflow for caller attribution.
- fx_quad_obj_5371.dual-r, `kyyvxyql` (stale after twin reuse; profiling build, three solves at 1999 Hz, `target/perf-members-fx.data`): F80 add_signed 8.06%, F80 nucleus 5.81%, multiply 5.74%, fresh F80 factor 5.17%, pack 3.83% self cycles. Solve remains infeasible, 5319 iterations, 5 refactorizations.
- degen3.dual-r, `tnpksmkx` (stale after direct lookup; frame-pointer build passes full parity; three solves capped at 2300 iterations, `target/perf-degen3dual.data`): iteration 12.06%, update 8.49%, BTRAN 7.16%, L stage 7.02%, FTRAN 6.46% self; ftran_needed 5.57% plus term filling 5.19%. FTRAN-needed is 21.14% inclusive, refactor 19.63%. All solves optimal in 2212 iterations, 24 refactorizations. The original cached Gurobi dual-controlled/release records also have 2212 iterations (no refactorization-count field); extra pivots do not explain the remaining gap. Annotated term filling (`target/perf-degen3dual-terms-annotate.txt`) attributes substantial time to repeated U-column scans, now replaced by the direct row map built from the numeric U transpose.

- mode_switch.dual-r, `xzuktkmx` (stale after F80/LU changes; reprofile before choosing it): update 22.21%, LU nucleus 15.90%, BTRAN 13.33%, fresh-factor row-list filtering 7.67%, factor_structure 3.36% self cycles. Full solve: 11393 iterations, 17 refactorizations.

## In progress

(none)

## Candidates

Ordered by current/ref; reprofile any changed path before choosing a change. Ratios use the current
sweep's solve medians; paired old/new times determine whether a change gains or regresses.

- degen3.dual-r (8.48x), degen3.dual-c (8.41x), perold.dual-r (8.33x), quad_2717.dual-r (7.70x): reprofile in this order. Profiles of changed paths are stale.
- DEGEN3 dual: tau row-value lookup is settled. Inspect remaining iteration control and immutable BTRAN metadata (transpose levels, inverse row permutation, tail nonzero counts are rebuilt/scanned in calls); preserve bucket visitation, arithmetic order and work estimates.
- Other reachable runs follow their table ratios. Quad twins share construction and four compact repeated factors are cached; row-capacity growth already uses the numeric U transpose.
- Quad_2717: after profiling again, inspect F80 arithmetic and the ranked row products. Cost-change work also performs required reduced-cost updates; some extended row-product outputs may be unused, but masks must preserve every consumed value and all work bookings. The too-small factor-cache working set is settled.
- Barrier: profile again before choosing more work. DGEMM now dispatches existing fused operations once per call; other BLAS-like kernels still use their original scalar dispatch and accumulation order. A smaller remaining DGEMM gap may involve independent outputs and data reuse, which must preserve every per-output reduction order.
- Runs without reachable references, after the reachable ratios: mode-switch update (22.21%) still scans A/B lists for other memberships; BTRAN (13.33%) allocates scratch/recomputes L levels; LU nucleus (15.90%); fresh-factor row-list filtering (7.67%). These old profiles are stale; annotate again before choosing a change.
- all: dump::on/record construction without dumps and allocation in FTRAN/BTRAN, pricing and ratio tests remain unmeasured leads.

## Tried, no gain

- Binary lookup in the U transpose for tau terms: fast parity 940/940 passes; DEGEN3 dual instructions 3769.4 -> 3753.4M (-0.4%), solve 0.739033 -> 0.716611 s is within old 0.025890 s spread. DEGEN2/PEROLD instructions rise slightly. Replaced by direct row-map trial (`target/perf-taulookup-target.log`).
- F80 add aligned with three guard bits, plus direct normal-result rounding: frozen old-add/pack differential test and fast parity 940/940 pass, instructions fall 4977.9 -> 4836.8M (-2.8%), but primal_5959 dual CPU-4 solve 0.555599 -> 0.557742 s (spreads 0.006716/0.011524 s). Alignment alone also has no timing gain on CPUs 2/4 despite -1.7% instructions. Reverted both (`target/perf-add-{target,recheck}.log`, `target/perf-add-direct-{target,recheck}.log`).
- Read eliminate work counts directly from maintained rowl lengths instead of rebuilding counts: parity full 987/987 and release tests pass, but mode_switch.dual-r 22.849 -> 21.851 s is smaller than the old 1.849 s spread (instructions 146886.1 -> 140945.4M, -4.0%). Alternating DEGEN2 primal 0.092/0.090 then 0.091/0.091 s; dual 0.068/0.067 s, marginal against spread. Reverted.

## Tried, breaks parity

(idea, first differing record from `rose-parity check`, why)
