# Performance state

Current state of the work in [PLAN.md](PLAN.md). Rewrite in place; no log (history is `jj log`).

## Numbers

Baseline solver: `8b645007` (same solver as `vlsokskq`). Current solver: `ukkruzyl` (fixed-width rounding for normal F80 pack results).
Parity: full tier 987/987; release tests pass. Three timed runs after one warm-up, hyperfine 2.0.0; instruction counts are millions from perf stat.

`baseline wall` is the fixed initial process-wall median. `paired base wall` is the fresh baseline-binary median from alternating baseline/last-landed/candidate checks on one CPU (CPU 0 unless stated otherwise); it resolves host drift without changing the fixed baseline. `current wall` and `solve` use the new harness; **current/ref uses solve-only time**, `diagnostics.seconds`, against original Gurobi Runtime.

Reference provenance: original-library cached records `../reverse-engineer/runs/{degen2,degen3,perold,pilot4}-{profile}.json` for the ten previously valid cells. New three-run medians use `/tmp/rose-perf-reference/{primal_5959.primal-r,primal_5959.dual-r,quad_2717.dual-r,fx_quad_obj_5371.dual-r}.{1,2,3}.json`; summary `/tmp/rose-perf-reference/summary.json`. All original-library records identify Ryzen 9 5900X and library SHA256 `8e87fca221fdffc922331c8619426031b3e72dfd419446422f88af0cb0e98f9c`.
`-` refs are unreachable: original Gurobi rejects 80BAU3B, GREENBEA and mode-switch with license error 10010 (scratch `ID.error.json`, including direct 80BAU3B primal/dual attempts). The old 80BAU3B primal/dual refs used the patched library and have been removed; the old prose mode-switch time is unusable.

Regression checks: all 23 runs paired against the preceding solver on CPU 4. Only mode-switch exceeded spread in that sweep; the CPU-2 repeat is old/new wall 12.925/13.924 s with 5.123/7.484 s spreads and unchanged instructions (79716.9/79720.6M), so the difference is within noise. Earlier frozen-baseline comparisons remain valid for the other rows. Raw data `/tmp/rose-pack-full-{before,after}.txt`, `/tmp/rose-pack-mode-{before,after}.txt`.
Acceptance: quad-objective full-sweep solve 0.313341 -> 0.298175 s (4.8%; spreads 0.009187/0.004654 s), instructions 2881.6 -> 2741.2M (-4.9%). Both CPU-4 short pairs also exceed their spreads: 0.472655 -> 0.452944 and 0.491316 -> 0.463440 s (`/tmp/rose-pack-cpu4.txt`). Quad_2717 full-sweep solve 0.653225 -> 0.634022 s, instructions 5747.9 -> 5397.2M. Normal packing uses fixed 64-bit extraction; subnormal precision, sticky/tie rounding, carry and overflow handling remain identical. Frozen packer tests cover 100000 random inputs, boundaries and dump-derived values, with existing x87 tests as an independent check.

| run | ref solve | baseline wall | paired base wall | current wall | solve | current/ref | instr |
|---|---|---|---|---|---|---|---|
| degen2.primal-c | 0.012711 | 0.182 | - | 0.097 | 0.091914 | 7.23 | 502.5 |
| degen2.dual-r | 0.009653 | 0.125 | - | 0.081 | 0.075775 | 7.85 | 374.4 |
| degen3.primal-c | 0.163137 | 3.237 | - | 0.933 | 0.922049 | 5.65 | 7988.3 |
| degen3.primal-r | 0.193685 | 3.872 | - | 1.103 | 1.092719 | 5.64 | 9838.1 |
| perold.primal-c | 0.160994 | 2.430 | - | 0.908 | 0.903945 | 5.61 | 7979.2 |
| pilot4.primal-r | 0.015787 | 0.253 | - | 0.095 | 0.090878 | 5.76 | 889.0 |
| 80bau3b.primal-c | - | 1.618 | - | 1.181 | 1.164780 | - | 10184.2 |
| primal_5959.primal-r | 0.067884 | 0.563 | 0.938 | 0.836 | 0.832888 | 12.27 | 4683.5 |
| degen3.dual-c | 0.082843 | 1.002 | - | 0.781 | 0.764445 | 9.23 | 3941.7 |
| degen3.dual-r | 0.082947 | 1.008 | 1.007 | 0.771 | 0.753315 | 9.08 | 3951.4 |
| perold.dual-r | 0.038532 | 0.355 | 0.360 | 0.359 | 0.351819 | 9.13 | 1941.4 |
| 80bau3b.dual-r | - | 1.103 | 1.688 | 0.908 | 0.884330 | - | 5359.1 |
| greenbea.dual-c | - | 2.870 | 2.363 | 2.095 | 2.073904 | - | 11169.4 |
| greenbea.dual-r | - | 3.598 | 3.592 | 1.350 | 1.334948 | - | 11072.6 |
| mode_switch.dual-r | - | 25.815 | - | 13.924 | 13.645143 | - | 79720.6 |
| primal_5959.dual-r | 0.098192 | 1.300 | - | 0.713 | 0.711132 | 7.24 | 6686.8 |
| quad_2717.dual-r | 0.092814 | 1.130 | - | 0.636 | 0.634022 | 6.83 | 5397.2 |
| fx_quad_obj_5371.dual-r | 0.038077 | 0.334 | 0.531 | 0.300 | 0.298175 | 7.83 | 2741.2 |
| degen3.barrier-r | 0.075214 | 0.608 | 0.856 | 0.623 | 0.581781 | 7.74 | 5498.0 |
| 80bau3b.barrier-r | - | 0.379 | 0.549 | 0.543 | 0.507582 | - | 4058.5 |
| greenbea.barrier-r | - | 1.569 | 2.493 | 2.242 | 2.192248 | - | 16053.5 |
| 80bau3b.pdhg | - | 4.552 | 6.526 | 4.447 | 4.433581 | - | 92389.5 |
| greenbea.pdhg | - | 2.801 | - | 1.899 | 1.887373 | - | 37741.5 |

## Hotspots

Per profiled bench run: the commit it was taken at and the top self-time functions (percent), at most
five lines each. Re-profile a run when a landed change touched its hot path.

- primal_5959.primal-r, `ukkruzyl` (profiling build, three bounded solves at 1999 Hz, `target/perf-pack-primal5959.data`): F80 add_signed 6.32%, F80 nucleus 5.86%, primal control 5.28%, free 4.88%, fresh F80 factor 4.42% self cycles; malloc is a further 3.66%. Current profile log `target/perf-pack-primal5959.log` verifies numerical_failure, 10780 iterations, 5381 refactorizations. Original reference's three records all have status 12 and 10780 iterations; no reference refactorization count is available.
- fx_quad_obj_5371.dual-r, `ukkruzyl` (profiling build, three solves at 1999 Hz): F80 add_signed 8.97%, F80 nucleus 7.17%, multiply 5.38%, fresh F80 factor 4.60%, pack 2.95% self cycles; F64 nucleus 4.15%. Pack was 4.72% in the preceding profile; solve remains infeasible, 5319 iterations, 5 refactorizations.
- degen3.dual-r, `yzqsoomv` (profiling build, three solves at 1999 Hz): iteration 10.59%, update 10.17%, L stage 8.00%, BTRAN 6.51%, ftran_with 6.51% self cycles. ftran_needed 5.64% plus term-buffer filling/column lookup 5.10%; the old per-row allocation is gone. Full solve: 2212 iterations, 24 refactorizations.
- mode_switch.dual-r, `xzuktkmx` (stale after F80 arithmetic change; reprofile before choosing it): update 22.21%, LU nucleus 15.90%, BTRAN 13.33%, fresh-factor row-list filtering 7.67%, factor_structure 3.36% self cycles. Full solve: 11393 iterations, 17 refactorizations.

## In progress

(none)

## Candidates

Ordered by expected payoff; one line each, with the bench run it targets. Start from PLAN.md "Where to
look first" and replace these with measured items after the first profiles.

- primal_5959.primal-r (12.27x, highest reachable ratio): inspect inclusive profile/annotation of F80 factor/nucleus construction, allocation and primal control. Both original reference and current rose take 10780 iterations; 5381 current refactorizations make repeated setup a lead to measure. Numeric changes require frozen-kernel differential tests.
- DEGEN3 dual-c (9.23x), PEROLD dual-r (9.13x), DEGEN3 dual-r (9.08x): follow in decreasing ratio.
- fx_quad_obj_5371.dual-r (7.83x): remaining F80 add/packing and factor construction; protect any further numeric rewrite with a frozen-kernel differential test first.
- Other reachable runs follow their table ratios. DEGEN3 term lookup/filling remains hot despite buffer reuse; any optimization must preserve row-dot order and work estimates.
- Runs without reachable references, after the reachable ratios: mode-switch update (22.21%) still scans A/B lists for other memberships and row-capacity accounting; BTRAN (13.33%) allocates scratch/recomputes L levels; LU nucleus (15.90%); fresh-factor row-list filtering (7.67%). Annotate before choosing a change.
- all: dump::on/record construction without dumps and allocation in FTRAN/BTRAN, pricing and ratio tests remain unmeasured leads.

## Tried, no gain

- Read eliminate work counts directly from maintained rowl lengths instead of rebuilding counts: parity full 987/987 and release tests pass, but mode_switch.dual-r 22.849 -> 21.851 s is smaller than the old 1.849 s spread (instructions 146886.1 -> 140945.4M, -4.0%). Alternating DEGEN2 primal 0.092/0.090 then 0.091/0.091 s; dual 0.068/0.067 s, marginal against spread. Reverted.

## Tried, breaks parity

(idea, first differing record from `rose-parity check`, why)
