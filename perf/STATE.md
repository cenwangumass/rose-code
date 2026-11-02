# Performance state

Current state of the work in [PLAN.md](PLAN.md). Rewrite in place; no log (history is `jj log`).

## Numbers

Baseline solver: `8b645007` (same solver as `vlsokskq`). Current solver: `yzqsoomv` (factor-owned scratch terms for second-vector row dots).
Parity: full tier 987/987; release tests pass. Three timed runs after one warm-up, hyperfine 2.0.0; instruction counts are millions from perf stat.

`baseline wall` is the fixed initial process-wall median. `paired base wall` is the fresh baseline-binary median from alternating baseline/last-landed/candidate checks on one CPU (CPU 0 unless stated otherwise); it resolves host drift without changing the fixed baseline. `current wall` and `solve` use the new harness; **current/ref uses solve-only time**, `diagnostics.seconds`, against original Gurobi Runtime.

Reference provenance: original-library cached records `../reverse-engineer/runs/{degen2,degen3,perold,pilot4}-{profile}.json` for the ten previously valid cells. New three-run medians use `/tmp/rose-perf-reference/{primal_5959.primal-r,primal_5959.dual-r,quad_2717.dual-r,fx_quad_obj_5371.dual-r}.{1,2,3}.json`; summary `/tmp/rose-perf-reference/summary.json`. All original-library records identify Ryzen 9 5900X and library SHA256 `8e87fca221fdffc922331c8619426031b3e72dfd419446422f88af0cb0e98f9c`.
`-` refs are unreachable: original Gurobi rejects 80BAU3B, GREENBEA and mode-switch with license error 10010 (scratch `ID.error.json`, including direct 80BAU3B primal/dual attempts). The old 80BAU3B primal/dual refs used the patched library and have been removed; the old prose mode-switch time is unusable.

Regression checks: all 23 runs paired against the preceding solver on CPU 2. Apparent primal_5959 primal and DEGEN3 controlled-dual slowdowns resolve in CPU-0 five-run pairs: 0.849/0.852 s (spreads 0.035/0.021 s) and 0.804/0.770 s. A fresh CPU-0 frozen-baseline 80BAU3B barrier pair is 0.549/0.562 s, below the baseline 0.014 s spread. Prior frozen-baseline checks remain valid for other rows. Raw data: `/tmp/rose-terms-full-{before,after}.txt`, `/tmp/rose-terms-rerun-{before,after}.txt`, `/tmp/rose-terms-barrier-{base,after}.txt`.
Acceptance: DEGEN3 dual-r five-run alternating solve median 0.505061 -> 0.496431 s (1.7%; spreads 0.002573/0.003129 s), instructions 4042.7 -> 3951.0M (-2.3%). The full sweep also reduces instructions on both DEGEN3 dual profiles by 2.3%; its release timing is noisier. Raw target checks: `/tmp/rose-terms-paired.txt`. Only allocation lifetime changes; row-dot kernels, term order, virtual peel and work estimates stay identical.

| run | ref solve | baseline wall | paired base wall | current wall | solve | current/ref | instr |
|---|---|---|---|---|---|---|---|
| degen2.primal-c | 0.012711 | 0.182 | - | 0.064 | 0.061189 | 4.81 | 502.7 |
| degen2.dual-r | 0.009653 | 0.125 | - | 0.052 | 0.048873 | 5.06 | 374.4 |
| degen3.primal-c | 0.163137 | 3.237 | - | 0.930 | 0.919499 | 5.64 | 7990.6 |
| degen3.primal-r | 0.193685 | 3.872 | - | 1.698 | 1.681676 | 8.68 | 9842.5 |
| perold.primal-c | 0.160994 | 2.430 | - | 1.398 | 1.390921 | 8.64 | 7981.4 |
| pilot4.primal-r | 0.015787 | 0.253 | - | 0.149 | 0.143329 | 9.08 | 889.1 |
| 80bau3b.primal-c | - | 1.618 | - | 1.479 | 1.457079 | - | 10185.9 |
| primal_5959.primal-r | 0.067884 | 0.563 | 0.938 | 0.852 | 0.849147 | 12.51 | 4863.9 |
| degen3.dual-c | 0.082843 | 1.002 | - | 0.770 | 0.753972 | 9.10 | 3941.3 |
| degen3.dual-r | 0.082947 | 1.008 | 1.007 | 0.794 | 0.777990 | 9.38 | 3951.2 |
| perold.dual-r | 0.038532 | 0.355 | 0.360 | 0.229 | 0.223798 | 5.81 | 1941.5 |
| 80bau3b.dual-r | - | 1.103 | 1.688 | 0.595 | 0.579966 | - | 5358.7 |
| greenbea.dual-c | - | 2.870 | 2.363 | 1.372 | 1.357692 | - | 11169.3 |
| greenbea.dual-r | - | 3.598 | 3.592 | 1.526 | 1.510660 | - | 11072.9 |
| mode_switch.dual-r | - | 25.815 | - | 7.946 | 7.791706 | - | 79718.1 |
| primal_5959.dual-r | 0.098192 | 1.300 | - | 1.196 | 1.193249 | 12.15 | 6945.4 |
| quad_2717.dual-r | 0.092814 | 1.130 | - | 1.005 | 1.001623 | 10.79 | 5750.6 |
| fx_quad_obj_5371.dual-r | 0.038077 | 0.334 | 0.531 | 0.489 | 0.486393 | 12.77 | 2882.1 |
| degen3.barrier-r | 0.075214 | 0.608 | 0.856 | 0.865 | 0.804553 | 10.70 | 5496.6 |
| 80bau3b.barrier-r | - | 0.379 | 0.549 | 0.562 | 0.526811 | - | 4058.4 |
| greenbea.barrier-r | - | 1.569 | 2.493 | 2.278 | 2.226387 | - | 16052.7 |
| 80bau3b.pdhg | - | 4.552 | 6.526 | 5.680 | 5.665724 | - | 92391.1 |
| greenbea.pdhg | - | 2.801 | - | 2.733 | 2.716738 | - | 37741.7 |

## Hotspots

Per profiled bench run: the commit it was taken at and the top self-time functions (percent), at most
five lines each. Re-profile a run when a landed change touched its hot path.

- fx_quad_obj_5371.dual-r, `yzqsoomv` (profiling build, three solves at 1999 Hz): F80 add_signed 10.03%, F80 nucleus 5.58%, F80 multiply 5.29%, fresh F80 factor 4.91%, pack 4.72% self cycles; infeasible, 5319 iterations, 5 refactorizations.
- degen3.dual-r, `yzqsoomv` (profiling build, three solves at 1999 Hz): iteration 10.59%, update 10.17%, L stage 8.00%, BTRAN 6.51%, ftran_with 6.51% self cycles. ftran_needed 5.64% plus term-buffer filling/column lookup 5.10%; the old per-row allocation is gone. Full solve: 2212 iterations, 24 refactorizations.
- mode_switch.dual-r, `xzuktkmx` (stale after F80 arithmetic change; reprofile before choosing it): update 22.21%, LU nucleus 15.90%, BTRAN 13.33%, fresh-factor row-list filtering 7.67%, factor_structure 3.36% self cycles. Full solve: 11393 iterations, 17 refactorizations.

## In progress

(none)

## Candidates

Ordered by expected payoff; one line each, with the bench run it targets. Start from PLAN.md "Where to
look first" and replace these with measured items after the first profiles.

- fx_quad_obj_5371.dual-r (12.77x, highest reachable ratio): inspect F80 add_signed and pack (10.03% and 4.72% self) and F80 nucleus/factor construction. Protect any numeric rewrite with a frozen-kernel differential test first.
- primal_5959.primal-r (12.51x) and primal_5959.dual-r (12.15x): profile next if the highest-ratio run has no actionable lead.
- quad_2717.dual-r (10.79x) and DEGEN3 barrier-r (10.70x) follow in decreasing ratio.
- degen3.dual-r (9.38x): row-dot term lookup/filling still costs 5.10% even with buffer reuse; inspect iteration/update/L-stage loops as well. Other reachable runs follow their table ratios.
- Runs without reachable references, after the reachable ratios: mode-switch update (22.21%) still scans A/B lists for other memberships and row-capacity accounting; BTRAN (13.33%) allocates scratch/recomputes L levels; LU nucleus (15.90%); fresh-factor row-list filtering (7.67%). Annotate before choosing a change.
- all: dump::on/record construction without dumps and allocation in FTRAN/BTRAN, pricing and ratio tests remain unmeasured leads.

## Tried, no gain

- Read eliminate work counts directly from maintained rowl lengths instead of rebuilding counts: parity full 987/987 and release tests pass, but mode_switch.dual-r 22.849 -> 21.851 s is smaller than the old 1.849 s spread (instructions 146886.1 -> 140945.4M, -4.0%). Alternating DEGEN2 primal 0.092/0.090 then 0.091/0.091 s; dual 0.068/0.067 s, marginal against spread. Reverted.

## Tried, breaks parity

(idea, first differing record from `rose-parity check`, why)
