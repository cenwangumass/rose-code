# Performance state

Current state of the work in [PLAN.md](PLAN.md). Rewrite in place; no log (history is `jj log`).

## Numbers

Baseline solver: `8b645007` (same solver as `vlsokskq`). Current solver: `xzuktkmx` (ordered leaving-row reference enumeration).
Parity: full tier 987/987; release tests pass. Three timed runs after one warm-up, hyperfine 2.0.0; instruction counts are millions from perf stat.

`baseline wall` is the fixed initial process-wall median. `paired base wall` is the fresh baseline-binary median from CPU-0 alternating baseline/last-landed/candidate checks; it resolves host drift without changing the fixed baseline. `current wall` and `solve` use the new harness; **current/ref uses solve-only time**, `diagnostics.seconds`, against original Gurobi Runtime.

Reference provenance: original-library cached records `../reverse-engineer/runs/{degen2,degen3,perold,pilot4}-{profile}.json` for the ten previously valid cells. New three-run medians use `/tmp/rose-perf-reference/{primal_5959.primal-r,primal_5959.dual-r,quad_2717.dual-r,fx_quad_obj_5371.dual-r}.{1,2,3}.json`; summary `/tmp/rose-perf-reference/summary.json`. All original-library records identify Ryzen 9 5900X and library SHA256 `8e87fca221fdffc922331c8619426031b3e72dfd419446422f88af0cb0e98f9c`.
`-` refs are unreachable: original Gurobi rejects 80BAU3B, GREENBEA and mode-switch with license error 10010 (scratch `ID.error.json`, including direct 80BAU3B primal/dual attempts). The old 80BAU3B primal/dual refs used the patched library and have been removed; the old prose mode-switch time is unusable.

Regressions resolved in two CPU-0 alternating passes: baseline/current wall medians quad-objective 0.531/0.508 s, primal_5959 primal 0.938/0.880 s, DEGEN3 barrier 0.856/0.887 s (difference below the baseline 0.047 s spread). Remaining apparent regressions versus the last landed solver also resolve within spread; 80BAU3B barrier's first 3 ms difference resolves to 0.385/0.386 s, below the 0.006 s baseline spread. 80BAU3B PDHG's refreshed 6.712 s also resolves in a fresh baseline/current pair, 6.526/4.548 s, with unchanged instructions (92391.0/92390.1M). Raw pairs: `/tmp/rose-regress-{baseline,landed,candidate}-{1,2}.txt` and `/tmp/rose-pdhg-regress-{before,after}.txt`.
Acceptance versus preceding solver: pinned mode-switch wall median 14.128 -> 7.839 s (spreads 0.091/3.007 s), instructions 89073.5 -> 79723.6M (-10.5%). Full bench coverage: `/tmp/rose-refcols-rest-{before,after}.txt`, `/tmp/rose-refcols-mode-{before,after}.txt`; solve-only refresh `/tmp/rose-solve-rest.txt` plus the alternating candidate rows.

| run | ref solve | baseline wall | paired base wall | current wall | solve | current/ref | instr |
|---|---|---|---|---|---|---|---|
| degen2.primal-c | 0.012711 | 0.182 | - | 0.065 | 0.061235 | 4.82 | 502.8 |
| degen2.dual-r | 0.009653 | 0.125 | - | 0.052 | 0.049222 | 5.10 | 380.8 |
| degen3.primal-c | 0.163137 | 3.237 | - | 0.933 | 0.921702 | 5.65 | 7991.0 |
| degen3.primal-r | 0.193685 | 3.872 | - | 1.733 | 1.709889 | 8.83 | 9843.2 |
| perold.primal-c | 0.160994 | 2.430 | - | 1.427 | 1.420118 | 8.82 | 7980.8 |
| pilot4.primal-r | 0.015787 | 0.253 | - | 0.150 | 0.144509 | 9.15 | 889.2 |
| 80bau3b.primal-c | - | 1.618 | - | 1.504 | 1.480420 | - | 10185.5 |
| primal_5959.primal-r | 0.067884 | 0.563 | 0.938 | 0.880 | 0.876730 | 12.92 | 5052.0 |
| degen3.dual-c | 0.082843 | 1.002 | - | 0.810 | 0.793303 | 9.58 | 4033.1 |
| degen3.dual-r | 0.082947 | 1.008 | 1.007 | 0.524 | 0.513408 | 6.19 | 4042.7 |
| perold.dual-r | 0.038532 | 0.355 | 0.360 | 0.231 | 0.225523 | 5.85 | 1949.9 |
| 80bau3b.dual-r | - | 1.103 | 1.688 | 0.930 | 0.907305 | - | 5372.5 |
| greenbea.dual-c | - | 2.870 | 2.363 | 1.379 | 1.364218 | - | 11224.2 |
| greenbea.dual-r | - | 3.598 | 3.592 | 2.112 | 2.090778 | - | 11120.9 |
| mode_switch.dual-r | - | 25.815 | - | 11.936 | 11.721287 | - | 79721.4 |
| primal_5959.dual-r | 0.098192 | 1.300 | - | 0.784 | 0.781992 | 7.96 | 7213.9 |
| quad_2717.dual-r | 0.092814 | 1.130 | - | 0.693 | 0.691301 | 7.45 | 6143.9 |
| fx_quad_obj_5371.dual-r | 0.038077 | 0.334 | 0.531 | 0.508 | 0.505127 | 13.27 | 3030.7 |
| degen3.barrier-r | 0.075214 | 0.608 | 0.856 | 0.887 | 0.823458 | 10.95 | 5501.2 |
| 80bau3b.barrier-r | - | 0.379 | 0.385 | 0.386 | 0.361727 | - | 4057.7 |
| greenbea.barrier-r | - | 1.569 | - | 1.561 | 1.527513 | - | 16051.9 |
| 80bau3b.pdhg | - | 4.552 | 6.526 | 4.548 | 4.534188 | - | 92390.1 |
| greenbea.pdhg | - | 2.801 | - | 3.130 | 3.098889 | - | 37741.7 |

## Hotspots

Per profiled bench run: the commit it was taken at and the top self-time functions (percent), at most
five lines each. Re-profile a run when a landed change touched its hot path.

- fx_quad_obj_5371.dual-r, `xzuktkmx` (profiling build, 999 Hz, 624 samples): F80 multiply 9.85%, F80 add_signed 8.00%, F80 LU nucleus 5.27%, F80 fresh factor 4.99%, free 3.61% self cycles. Full solve: infeasible, 5319 iterations, 5 refactorizations. Inspect inclusive callers and allocation before choosing the next change.
- mode_switch.dual-r, `xzuktkmx` (profiling build, 199 Hz): factor update 22.21%, LU nucleus 15.90%, BTRAN 13.33%, fresh-factor row-list filtering 7.67%, factor_structure 3.36% self cycles. Full solve: 11393 iterations, 17 refactorizations.

## In progress

(none)

## Candidates

Ordered by expected payoff; one line each, with the bench run it targets. Start from PLAN.md "Where to
look first" and replace these with measured items after the first profiles.

- fx_quad_obj_5371.dual-r (13.27x, highest reachable ratio): F80 multiply/add and F80 factor construction; inspect inclusive profiles and allocation, then protect any numeric rewrite with differential tests before changing it.
- primal_5959.primal-r (12.92x): profile next if the higher-ratio run has no actionable lead.
- degen3.barrier-r (10.95x): profile dense/sparse factor kernels; other reachable runs follow in decreasing table ratio.
- Runs without reachable references, after the reachable ratios: mode-switch update (22.21%) still scans A/B lists for other memberships and row-capacity accounting; BTRAN (13.33%) allocates scratch/recomputes L levels; LU nucleus (15.90%); fresh-factor row-list filtering (7.67%). Annotate before choosing a change.
- all: dump::on/record construction without dumps and allocation in FTRAN/BTRAN, pricing and ratio tests remain unmeasured leads.

## Tried, no gain

- Read eliminate work counts directly from maintained rowl lengths instead of rebuilding counts: parity full 987/987 and release tests pass, but mode_switch.dual-r 22.849 -> 21.851 s is smaller than the old 1.849 s spread (instructions 146886.1 -> 140945.4M, -4.0%). Alternating DEGEN2 primal 0.092/0.090 then 0.091/0.091 s; dual 0.068/0.067 s, marginal against spread. Reverted.

## Tried, breaks parity

(idea, first differing record from `rose-parity check`, why)
