# Performance state

Current state of the work in [PLAN.md](PLAN.md). Rewrite in place; no log (history is `jj log`).

## Numbers

Baseline solver: `8b645007` (same solver as `vlsokskq`). Current solver: `zzvmznmp` (dispatch support for existing fused DGEMM operations once per kernel call).
Parity: full tier 987/987; release tests pass. Three timed runs after one warm-up, hyperfine 2.0.0; instruction counts are millions from perf stat.

`baseline wall` is the fixed initial process-wall median. `paired base wall` is the fresh baseline-binary median from alternating baseline/last-landed/candidate checks on one CPU (CPU 0 unless stated otherwise); it resolves host drift without changing the fixed baseline. `current wall` and `solve` use the new harness; **current/ref uses solve-only time**, `diagnostics.seconds`, against original Gurobi Runtime.

Reference provenance: original-library cached records `../reverse-engineer/runs/{degen2,degen3,perold,pilot4}-{profile}.json` for the ten previously valid cells. New three-run medians use `/tmp/rose-perf-reference/{primal_5959.primal-r,primal_5959.dual-r,quad_2717.dual-r,fx_quad_obj_5371.dual-r}.{1,2,3}.json`; summary `/tmp/rose-perf-reference/summary.json`. All original-library records identify Ryzen 9 5900X and library SHA256 `8e87fca221fdffc922331c8619426031b3e72dfd419446422f88af0cb0e98f9c`.
`-` refs are unreachable: original Gurobi rejects 80BAU3B, GREENBEA and mode-switch with license error 10010 (scratch `ID.error.json`, including direct 80BAU3B primal/dual attempts). The old 80BAU3B primal/dual refs used the patched library and have been removed; the old prose mode-switch time is unusable.

Regression checks: all 23 runs paired against the preceding solver on CPU 4. PEROLD dual's small shift resolves in five-run CPU-2 pairs: old/new solve 0.323283/0.328716 s, spreads 0.010717/0.010160 s, unchanged instructions. GREENBEA dual-controlled shifts between host-speed bands; three interleaved single-run harness samples per binary on CPU 0 give old/new solve medians 1.316744/1.363929 s, spreads 0.706242/0.695403 s, matching instructions. The table uses that interleaved median (each sample warmed separately). Fresh CPU-2 80BAU3B barrier frozen-baseline/prior/current solves are 0.510930/0.521363/0.375277 s; no fixed-baseline regression. GREENBEA PDHG's 1.10 s wall spread covers its fixed-baseline difference, with unchanged instructions. No unresolved regression. Raw data `target/perf-dgemm-full-{before,after}.txt`, `target/perf-dgemm-recheck-{base,before,after}.txt`, `target/perf-dgemm-greenbea-interleaved.{json,log}`. Host speed varies substantially; compare immediate pairs and spreads.
Acceptance: DEGEN3 barrier solve 0.832103 -> 0.406835 s (-51.1%; spreads 0.068020/0.007927 s), instructions 5499.7 -> 3779.9M (-31.3%). CPU-2 80BAU3B barrier solve 0.521363 -> 0.375277 s (-28.0%), instructions 4058.3 -> 3425.3M (-15.6%). GREENBEA barrier solve 2.178043 -> 1.479103 s (-32.1%), instructions 16053.3 -> 13029.7M (-18.8%). BORE3D process wall 15.825 -> 13.483 ms, instructions 94.83 -> 84.35M. Existing fused and unfused operations retain their order; runtime CPU dispatch happens once per DGEMM. Non-finite input arrays use the portable sequence to preserve NaN payload priority. A frozen old kernel/helper differential test covers both paths, random magnitudes/shapes/padding, signed zeros, overflow/underflow, infinities, NaNs and cached real AFIRO inputs.

| run | ref solve | baseline wall | paired base wall | current wall | solve | current/ref | instr |
|---|---|---|---|---|---|---|---|
| degen2.primal-c | 0.012711 | 0.182 | - | 0.062 | 0.058835 | 4.63 | 481.5 |
| degen2.dual-r | 0.009653 | 0.125 | - | 0.050 | 0.046867 | 4.86 | 358.3 |
| degen3.primal-c | 0.163137 | 3.237 | - | 1.369 | 1.353218 | 8.29 | 7652.9 |
| degen3.primal-r | 0.193685 | 3.872 | - | 1.618 | 1.603129 | 8.28 | 9400.4 |
| perold.primal-c | 0.160994 | 2.430 | - | 1.296 | 1.289131 | 8.01 | 7535.7 |
| pilot4.primal-r | 0.015787 | 0.253 | - | 0.134 | 0.128468 | 8.14 | 820.6 |
| 80bau3b.primal-c | - | 1.618 | - | 1.468 | 1.444111 | - | 10090.2 |
| primal_5959.primal-r | 0.067884 | 0.563 | 0.938 | 0.464 | 0.460603 | 6.79 | 2858.3 |
| degen3.dual-c | 0.082843 | 1.002 | - | 0.483 | 0.472536 | 5.70 | 3758.8 |
| degen3.dual-r | 0.082947 | 1.008 | 1.007 | 0.485 | 0.474712 | 5.72 | 3769.1 |
| perold.dual-r | 0.038532 | 0.355 | 0.360 | 0.336 | 0.328716 | 8.53 | 1819.2 |
| 80bau3b.dual-r | - | 1.103 | 1.688 | 0.566 | 0.550954 | - | 5155.0 |
| greenbea.dual-c | - | 2.870 | 2.363 | 1.379 | 1.363929 | - | 10816.7 |
| greenbea.dual-r | - | 3.598 | 3.592 | 1.986 | 1.966605 | - | 10715.4 |
| mode_switch.dual-r | - | 25.815 | - | 7.299 | 7.151443 | - | 76028.5 |
| primal_5959.dual-r | 0.098192 | 1.300 | - | 0.631 | 0.627833 | 6.39 | 4012.8 |
| quad_2717.dual-r | 0.092814 | 1.130 | - | 0.862 | 0.859014 | 9.26 | 4649.2 |
| fx_quad_obj_5371.dual-r | 0.038077 | 0.334 | 0.514 | 0.293 | 0.289284 | 7.60 | 1619.8 |
| degen3.barrier-r | 0.075214 | 0.608 | 0.856 | 0.468 | 0.406835 | 5.41 | 3779.9 |
| 80bau3b.barrier-r | - | 0.379 | 0.546 | 0.411 | 0.375277 | - | 3425.3 |
| greenbea.barrier-r | - | 1.569 | 2.493 | 1.530 | 1.479103 | - | 13029.7 |
| 80bau3b.pdhg | - | 4.552 | 6.526 | 4.513 | 4.500036 | - | 92390.6 |
| greenbea.pdhg | - | 2.801 | - | 3.303 | 3.286155 | - | 37741.1 |

## Hotspots

Per profiled bench run: the commit it was taken at and the top self-time functions (percent), at most
five lines each. Re-profile a run when a landed change touched its hot path.

- degen3.barrier-r, `munmovyt` (stale after kernel dispatch; frame-pointer build passes full parity; three solves capped at 100 iterations, `target/perf-barrier-degen3-valid.data`): FMA dispatch 26.46%, DGEMM 24.62%, existing hardware FMA implementation 14.77%, numeric_factor 7.40%, DTRSM 4.36% self. DGEMM 49.04%, numeric factor 81.45% inclusive; all three solves optimal in 14 iterations. The per-product FMA wrapper is an indirect call plus frame/return overhead (annotated `target/perf-barrier-fma-annotate.txt`); per-call dispatch now removes that per-product overhead without adding fused operations. BORE3D is a cheap witness (16 barrier iterations, 12 ms).
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

- quad_2717.dual-r (9.26x), perold.dual-r (8.53x), degen3.primal-c (8.29x), degen3.primal-r (8.28x): reprofile in this order. Profiles of changed paths are stale. Quad_2717's prior +0.8% instruction cost is a lead for avoiding fresh-cache copies on non-repeated inputs; its timing remains within spread.
- Other reachable runs follow their table ratios. DEGEN3 row-dot filling/lookup remains an untried lead after buffer reuse; preserve accumulation order and work estimates. Quad twins share construction and compact repeated factors are cached; row-capacity growth already uses the numeric U transpose.
- Barrier: profile again before choosing more work. DGEMM now dispatches existing fused operations once per call; other BLAS-like kernels still use their original scalar dispatch and accumulation order. A smaller remaining DGEMM gap may involve independent outputs and data reuse, which must preserve every per-output reduction order.
- Runs without reachable references, after the reachable ratios: mode-switch update (22.21%) still scans A/B lists for other memberships; BTRAN (13.33%) allocates scratch/recomputes L levels; LU nucleus (15.90%); fresh-factor row-list filtering (7.67%). These old profiles are stale; annotate again before choosing a change.
- all: dump::on/record construction without dumps and allocation in FTRAN/BTRAN, pricing and ratio tests remain unmeasured leads.

## Tried, no gain

- F80 add aligned with three guard bits, plus direct normal-result rounding: frozen old-add/pack differential test and fast parity 940/940 pass, instructions fall 4977.9 -> 4836.8M (-2.8%), but primal_5959 dual CPU-4 solve 0.555599 -> 0.557742 s (spreads 0.006716/0.011524 s). Alignment alone also has no timing gain on CPUs 2/4 despite -1.7% instructions. Reverted both (`target/perf-add-{target,recheck}.log`, `target/perf-add-direct-{target,recheck}.log`).
- Read eliminate work counts directly from maintained rowl lengths instead of rebuilding counts: parity full 987/987 and release tests pass, but mode_switch.dual-r 22.849 -> 21.851 s is smaller than the old 1.849 s spread (instructions 146886.1 -> 140945.4M, -4.0%). Alternating DEGEN2 primal 0.092/0.090 then 0.091/0.091 s; dual 0.068/0.067 s, marginal against spread. Reverted.

## Tried, breaks parity

(idea, first differing record from `rose-parity check`, why)
