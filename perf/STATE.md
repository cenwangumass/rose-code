# Performance state

Current state of the work in [PLAN.md](PLAN.md). Rewrite in place; no log (history is `jj log`).

## Numbers

Baseline solver: `8b645007` (same solver as `vlsokskq`). Current solver: `lsvzktvr` (specialized normal-result F80 multiplication).
Parity: full tier 987/987; release tests pass. Three timed runs after one warm-up, hyperfine 2.0.0; instruction counts are millions from perf stat.

`baseline wall` is the fixed initial process-wall median. `paired base wall` is the fresh baseline-binary median from alternating baseline/last-landed/candidate checks on one CPU (CPU 0 unless stated otherwise); it resolves host drift without changing the fixed baseline. `current wall` and `solve` use the new harness; **current/ref uses solve-only time**, `diagnostics.seconds`, against original Gurobi Runtime.

Reference provenance: original-library cached records `../reverse-engineer/runs/{degen2,degen3,perold,pilot4}-{profile}.json` for the ten previously valid cells. New three-run medians use `/tmp/rose-perf-reference/{primal_5959.primal-r,primal_5959.dual-r,quad_2717.dual-r,fx_quad_obj_5371.dual-r}.{1,2,3}.json`; summary `/tmp/rose-perf-reference/summary.json`. All original-library records identify Ryzen 9 5900X and library SHA256 `8e87fca221fdffc922331c8619426031b3e72dfd419446422f88af0cb0e98f9c`.
`-` refs are unreachable: original Gurobi rejects 80BAU3B, GREENBEA and mode-switch with license error 10010 (scratch `ID.error.json`, including direct 80BAU3B primal/dual attempts). The old 80BAU3B primal/dual refs used the patched library and have been removed; the old prose mode-switch time is unusable.

Regression checks: all 23 runs compared back to back against the preceding solver on CPU 0. Apparent mode-switch and GREENBEA barrier slowdowns resolve on CPU 2: mode-switch old/new wall 12.168/8.201 s with unchanged instructions; GREENBEA barrier preceding/new 2.215/2.243 s (spreads 0.081/0.090 s), frozen baseline 2.493 s. Prior frozen-baseline checks remain valid for the other rows; raw data `/tmp/rose-mul-full-{before,after}.txt`, `/tmp/rose-mul-mode-{before,after}.txt`, `/tmp/rose-mul-green-{baseline,before,after}.txt`.
Acceptance: quad-objective full-sweep solve 0.326327 -> 0.309840 s (spreads 0.001772/0.002270 s), instructions 3029.2 -> 2881.9M (-4.9%). Quad_2717 solve 0.705755 -> 0.657029 s (spreads 0.005557/0.006383 s), instructions 6143.0 -> 5746.7M (-6.4%). Five-run alternating target checks also agree (`/tmp/rose-mul-paired.txt`). New scalar kernel is protected by the frozen multiply, 100000 random pairs, signed zeros/extremes, a real QTAU column and existing x87 hardware tests.

| run | ref solve | baseline wall | paired base wall | current wall | solve | current/ref | instr |
|---|---|---|---|---|---|---|---|
| degen2.primal-c | 0.012711 | 0.182 | - | 0.099 | 0.093364 | 7.35 | 502.9 |
| degen2.dual-r | 0.009653 | 0.125 | - | 0.081 | 0.076921 | 7.97 | 380.8 |
| degen3.primal-c | 0.163137 | 3.237 | - | 1.478 | 1.462675 | 8.97 | 7990.7 |
| degen3.primal-r | 0.193685 | 3.872 | - | 1.109 | 1.098454 | 5.67 | 9842.2 |
| perold.primal-c | 0.160994 | 2.430 | - | 0.906 | 0.901287 | 5.60 | 7982.1 |
| pilot4.primal-r | 0.015787 | 0.253 | - | 0.094 | 0.090747 | 5.75 | 889.2 |
| 80bau3b.primal-c | - | 1.618 | - | 0.975 | 0.960620 | - | 10183.7 |
| primal_5959.primal-r | 0.067884 | 0.563 | 0.938 | 0.527 | 0.524859 | 7.73 | 4864.3 |
| degen3.dual-c | 0.082843 | 1.002 | - | 0.782 | 0.765662 | 9.24 | 4033.1 |
| degen3.dual-r | 0.082947 | 1.008 | 1.007 | 0.805 | 0.789787 | 9.52 | 4042.8 |
| perold.dual-r | 0.038532 | 0.355 | 0.360 | 0.357 | 0.349777 | 9.08 | 1950.0 |
| 80bau3b.dual-r | - | 1.103 | 1.688 | 0.930 | 0.907127 | - | 5371.3 |
| greenbea.dual-c | - | 2.870 | 2.363 | 2.076 | 2.055692 | - | 11223.3 |
| greenbea.dual-r | - | 3.598 | 3.592 | 2.121 | 2.099393 | - | 11122.0 |
| mode_switch.dual-r | - | 25.815 | - | 8.201 | 8.001374 | - | 79723.8 |
| primal_5959.dual-r | 0.098192 | 1.300 | - | 0.754 | 0.751796 | 7.66 | 6945.2 |
| quad_2717.dual-r | 0.092814 | 1.130 | - | 0.659 | 0.657029 | 7.08 | 5746.7 |
| fx_quad_obj_5371.dual-r | 0.038077 | 0.334 | 0.531 | 0.312 | 0.309840 | 8.14 | 2881.9 |
| degen3.barrier-r | 0.075214 | 0.608 | 0.856 | 0.623 | 0.581122 | 7.73 | 5500.1 |
| 80bau3b.barrier-r | - | 0.379 | 0.385 | 0.386 | 0.362451 | - | 4058.3 |
| greenbea.barrier-r | - | 1.569 | 2.493 | 2.243 | 2.195781 | - | 16053.2 |
| 80bau3b.pdhg | - | 4.552 | 6.526 | 6.576 | 6.556314 | - | 92390.1 |
| greenbea.pdhg | - | 2.801 | - | 1.924 | 1.912006 | - | 37741.2 |

## Hotspots

Per profiled bench run: the commit it was taken at and the top self-time functions (percent), at most
five lines each. Re-profile a run when a landed change touched its hot path.

- degen3.dual-r, `lsvzktvr` (profiling build, 999 Hz, 555 samples): update 10.10%, Dual::iteration 10.05%, L stage 7.04%, BTRAN 6.75%, ftran_needed 5.57% self cycles. Term-vector construction inside ftran_needed is a further 5.46%, ftran_with 5.09%. Full solve: 2212 iterations, 24 refactorizations; inspect inclusive callers/annotation before choosing the next change.
- fx_quad_obj_5371.dual-r, `lsvzktvr` (profiling build, three solves at 1999 Hz): F80 add_signed 9.66%, F80 LU nucleus 6.31%, fresh F80 factor 5.03%, multiply 4.98%, pack 4.41% self cycles. Multiply was 9.85% in the preceding profile; solve remains infeasible, 5319 iterations, 5 refactorizations.
- mode_switch.dual-r, `xzuktkmx` (stale after F80 arithmetic change; reprofile before choosing it): update 22.21%, LU nucleus 15.90%, BTRAN 13.33%, fresh-factor row-list filtering 7.67%, factor_structure 3.36% self cycles. Full solve: 11393 iterations, 17 refactorizations.

## In progress

(none)

## Candidates

Ordered by expected payoff; one line each, with the bench run it targets. Start from PLAN.md "Where to
look first" and replace these with measured items after the first profiles.

- degen3.dual-r (9.52x, highest reachable ratio): inspect iteration/update and L-stage loops; ftran_needed builds temporary term vectors (5.46% self time) before each row-dot product. Preserve accumulation order and booked work.
- degen3.dual-c (9.24x), PEROLD dual-r (9.08x), DEGEN3 primal-c (8.97x): follow in decreasing ratio if the top run has no actionable lead.
- fx_quad_obj_5371.dual-r (8.14x): remaining F80 add/packing and factor construction; any further numeric rewrite needs a frozen-kernel differential test first.
- Other reachable runs follow in decreasing table ratio, including primal_5959.primal-r (7.73x) and DEGEN3 barrier-r (7.73x).
- Runs without reachable references, after the reachable ratios: mode-switch update (22.21%) still scans A/B lists for other memberships and row-capacity accounting; BTRAN (13.33%) allocates scratch/recomputes L levels; LU nucleus (15.90%); fresh-factor row-list filtering (7.67%). Annotate before choosing a change.
- all: dump::on/record construction without dumps and allocation in FTRAN/BTRAN, pricing and ratio tests remain unmeasured leads.

## Tried, no gain

- Read eliminate work counts directly from maintained rowl lengths instead of rebuilding counts: parity full 987/987 and release tests pass, but mode_switch.dual-r 22.849 -> 21.851 s is smaller than the old 1.849 s spread (instructions 146886.1 -> 140945.4M, -4.0%). Alternating DEGEN2 primal 0.092/0.090 then 0.091/0.091 s; dual 0.068/0.067 s, marginal against spread. Reverted.

## Tried, breaks parity

(idea, first differing record from `rose-parity check`, why)
