# Performance state

Current state of the work in [PLAN.md](PLAN.md). Rewrite in place; no log (history is `jj log`).

## Numbers

Baseline solver: `8b645007` (same solver as `vlsokskq`). Current solver: `tplnupps` (read row-capacity counts from the maintained numeric U transpose).
Parity: full tier 987/987; release tests pass. Three timed runs after one warm-up, hyperfine 2.0.0; instruction counts are millions from perf stat.

`baseline wall` is the fixed initial process-wall median. `paired base wall` is the fresh baseline-binary median from alternating baseline/last-landed/candidate checks on one CPU (CPU 0 unless stated otherwise); it resolves host drift without changing the fixed baseline. `current wall` and `solve` use the new harness; **current/ref uses solve-only time**, `diagnostics.seconds`, against original Gurobi Runtime.

Reference provenance: original-library cached records `../reverse-engineer/runs/{degen2,degen3,perold,pilot4}-{profile}.json` for the ten previously valid cells. New three-run medians use `/tmp/rose-perf-reference/{primal_5959.primal-r,primal_5959.dual-r,quad_2717.dual-r,fx_quad_obj_5371.dual-r}.{1,2,3}.json`; summary `/tmp/rose-perf-reference/summary.json`. All original-library records identify Ryzen 9 5900X and library SHA256 `8e87fca221fdffc922331c8619426031b3e72dfd419446422f88af0cb0e98f9c`.
`-` refs are unreachable: original Gurobi rejects 80BAU3B, GREENBEA and mode-switch with license error 10010 (scratch `ID.error.json`, including direct 80BAU3B primal/dual attempts). The old 80BAU3B primal/dual refs used the patched library and have been removed; the old prose mode-switch time is unusable.

Regression checks: all 23 runs paired against the preceding solver on CPU 4. DEGEN3 dual-release and mode-switch timing shifts resolve on CPU 2: DEGEN3 old/new solve 0.738030/0.722804 s, instructions 3867.2/3769.0M; mode-switch 11.771787/8.007936 s (spreads 1.613853/4.368599 s), instructions 79565.7/76028.8M. The mixed-speed DEGEN3 primal-controlled median is replaced by its CPU-2 pair (1.361215/1.319368 s). Other increases lie within spread; no new unresolved frozen-baseline regression. Raw data `target/perf-rowcap-full-{before,after}.txt`, `target/perf-rowcap-recheck-{before,after}.txt`. Host speed varies substantially across the sweep; compare immediate pairs rather than medians from different sweeps.
Acceptance: PILOT4 primal solve 0.085617 -> 0.083791 s (-2.1%; spreads 0.000726/0.000653 s), instructions 831.0 -> 820.6M (-1.3%). Five-run target pairs on CPUs 0, 2 and 4 all reduce instructions (`target/perf-rowcap-target.log`, `target/perf-rowcap-pilot-recheck.log`); CPU 4 also clears timing spread. DEGEN3 dual-controlled solve 0.486865 -> 0.472392 s (-3.0%), instructions 3856.4 -> 3758.8M (-2.5%). Row-capacity growth reads numeric transpose counts only for the rows receiving an entry, preserving duplicate multiplicities, virtual pool layout and work estimates.

| run | ref solve | baseline wall | paired base wall | current wall | solve | current/ref | instr |
|---|---|---|---|---|---|---|---|
| degen2.primal-c | 0.012711 | 0.182 | - | 0.093 | 0.088805 | 6.99 | 481.4 |
| degen2.dual-r | 0.009653 | 0.125 | - | 0.078 | 0.072954 | 7.56 | 358.3 |
| degen3.primal-c | 0.163137 | 3.237 | - | 1.335 | 1.319368 | 8.09 | 7652.8 |
| degen3.primal-r | 0.193685 | 3.872 | - | 1.626 | 1.599373 | 8.26 | 9400.1 |
| perold.primal-c | 0.160994 | 2.430 | - | 0.856 | 0.850545 | 5.28 | 7534.7 |
| pilot4.primal-r | 0.015787 | 0.253 | - | 0.088 | 0.083791 | 5.31 | 820.6 |
| 80bau3b.primal-c | - | 1.618 | - | 0.939 | 0.924169 | - | 10089.0 |
| primal_5959.primal-r | 0.067884 | 0.563 | 0.938 | 0.388 | 0.385742 | 5.68 | 3516.9 |
| degen3.dual-c | 0.082843 | 1.002 | - | 0.482 | 0.472392 | 5.70 | 3758.8 |
| degen3.dual-r | 0.082947 | 1.008 | 1.007 | 0.738 | 0.722804 | 8.71 | 3769.0 |
| perold.dual-r | 0.038532 | 0.355 | 0.360 | 0.298 | 0.291495 | 7.57 | 1819.8 |
| 80bau3b.dual-r | - | 1.103 | 1.688 | 0.863 | 0.839693 | - | 5154.8 |
| greenbea.dual-c | - | 2.870 | 2.363 | 2.009 | 1.989300 | - | 10817.0 |
| greenbea.dual-r | - | 3.598 | 3.592 | 1.984 | 1.963932 | - | 10715.9 |
| mode_switch.dual-r | - | 25.815 | - | 8.202 | 8.007936 | - | 76028.8 |
| primal_5959.dual-r | 0.098192 | 1.300 | - | 0.934 | 0.930954 | 9.48 | 4978.1 |
| quad_2717.dual-r | 0.092814 | 1.130 | - | 0.544 | 0.541217 | 5.83 | 4610.4 |
| fx_quad_obj_5371.dual-r | 0.038077 | 0.334 | 0.531 | 0.228 | 0.225452 | 5.92 | 1994.8 |
| degen3.barrier-r | 0.075214 | 0.608 | 0.856 | 0.610 | 0.569710 | 7.57 | 5501.3 |
| 80bau3b.barrier-r | - | 0.379 | 0.549 | 0.380 | 0.356297 | - | 4057.0 |
| greenbea.barrier-r | - | 1.569 | 2.493 | 1.549 | 1.515951 | - | 16049.9 |
| 80bau3b.pdhg | - | 4.552 | 6.526 | 6.414 | 6.393266 | - | 92390.8 |
| greenbea.pdhg | - | 2.801 | - | 1.901 | 1.889849 | - | 37740.7 |

## Hotspots

Per profiled bench run: the commit it was taken at and the top self-time functions (percent), at most
five lines each. Re-profile a run when a landed change touched its hot path.

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

- primal_5959.dual-r (9.48x), degen3.dual-r (8.71x), degen3.primal-r (8.26x), degen3.primal-c (8.09x): reprofile in this order. All existing profiles of these paths are stale after the factor and row-capacity changes.
- Other reachable runs follow their table ratios. DEGEN3 row-dot filling/lookup remains an untried lead after buffer reuse; preserve accumulation order and work estimates. Fresh quad twins already share construction, and row-capacity growth already uses the numeric U transpose.
- Runs without reachable references, after the reachable ratios: mode-switch update (22.21%) still scans A/B lists for other memberships; BTRAN (13.33%) allocates scratch/recomputes L levels; LU nucleus (15.90%); fresh-factor row-list filtering (7.67%). These old profiles are stale; annotate again before choosing a change.
- all: dump::on/record construction without dumps and allocation in FTRAN/BTRAN, pricing and ratio tests remain unmeasured leads.

## Tried, no gain

- Read eliminate work counts directly from maintained rowl lengths instead of rebuilding counts: parity full 987/987 and release tests pass, but mode_switch.dual-r 22.849 -> 21.851 s is smaller than the old 1.849 s spread (instructions 146886.1 -> 140945.4M, -4.0%). Alternating DEGEN2 primal 0.092/0.090 then 0.091/0.091 s; dual 0.068/0.067 s, marginal against spread. Reverted.

## Tried, breaks parity

(idea, first differing record from `rose-parity check`, why)
