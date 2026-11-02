# Performance state

Current state of the work in [PLAN.md](PLAN.md). Rewrite in place; no log (history is `jj log`).

## Numbers

Baseline solver: `8b645007` (same code as `vlsokskq`, rose-parity gate). Current solver: `ttklzlkq` (L transpose cached; numeric U transpose maintained across updates).
Parity: full tier 987/987; release tests pass. Baseline measurements: three timed runs, one warm-up, hyperfine 2.0.0; instruction counts from perf stat. Wall time is variable on this host; use paired reruns before accepting gains.

Seconds, `perf/bench.sh` medians (release, no dumps); `ref` is the reference's `result.runtime` from
`../reverse-engineer/runs/` (`-`: no reference run with a runtime). The baseline stays fixed. Current timings/instructions use the full sweep, superseded by the latest paired rerun where available; instruction counts are millions (`--instr`).

Acceptance pairs versus the preceding solver (three timed runs): DEGEN3 dual-r 1.438 -> 0.874 s (old/new spreads 0.386/0.010 s), 6002.6 -> 4218.8M instructions; GREENBEA dual-r 2.850 -> 1.433 s (spreads 0.079/0.048 s), 14461.6 -> 11573.5M. Mode-switch instructions fall 39.4%, 146885.8 -> 89076.5M; median 18.534 -> 12.040 s is below the old 8.779 s spread, so its wall-time gain remains unproven.
CPU-0 alternating reruns resolve the full sweep's apparent regressions: PEROLD dual-r 0.484/0.381 s, 2577.3/2010.8M instructions; 80BAU3B barrier 0.553/0.563 s (spreads 0.022/0.024 s), GREENBEA barrier 2.271/2.271 s. Barrier/PDHG instruction counts remain unchanged within noise.

| run | ref | baseline | current | current/ref | instr |
|---|---|---|---|---|---|
| degen2.primal-c | 0.013 | 0.182 | 0.108 | 8.3 | 516.7 |
| degen2.dual-r | 0.010 | 0.125 | 0.086 | 8.6 | 393.7 |
| degen3.primal-c | 0.163 | 3.237 | 1.591 | 9.8 | 8239.9 |
| degen3.primal-r | 0.194 | 3.872 | 1.887 | 9.7 | 10286.1 |
| perold.primal-c | 0.161 | 2.430 | 1.488 | 9.2 | 8167.8 |
| pilot4.primal-r | 0.016 | 0.253 | 0.160 | 10.0 | 911.2 |
| 80bau3b.primal-c | 0.164 | 1.618 | 1.506 | 9.2 | 10448.3 |
| primal_5959.primal-r | - | 0.563 | 0.841 | - | 5054.2 |
| degen3.dual-c | 0.083 | 1.002 | 0.850 | 10.2 | 4208.2 |
| degen3.dual-r | 0.083 | 1.008 | 0.874 | 10.5 | 4218.8 |
| perold.dual-r | 0.039 | 0.355 | 0.381 | 9.8 | 2010.8 |
| 80bau3b.dual-r | 0.134 | 1.103 | 1.006 | 7.5 | 5712.7 |
| greenbea.dual-c | - | 2.870 | 1.703 | - | 11651.8 |
| greenbea.dual-r | - | 3.598 | 1.433 | - | 11573.5 |
| mode_switch.dual-r | 0.4 (reverse-engineer STATE.md) | 25.815 | 12.040 | 30.1 | 89076.5 |
| primal_5959.dual-r | - | 1.300 | 1.175 | - | 7213.1 |
| quad_2717.dual-r | - | 1.130 | 1.096 | - | 6150.7 |
| fx_quad_obj_5371.dual-r | - | 0.334 | 0.726 | - | 3034.5 |
| degen3.barrier-r | 0.075 | 0.608 | 0.857 | 11.4 | 5499.1 |
| 80bau3b.barrier-r | - | 0.379 | 0.563 | - | 4058.9 |
| greenbea.barrier-r | - | 1.569 | 2.271 | - | 16050.9 |
| 80bau3b.pdhg | - | 4.552 | 4.494 | - | 92389.5 |
| greenbea.pdhg | - | 2.801 | 1.903 | - | 37741.2 |

## Hotspots

Per profiled bench run: the commit it was taken at and the top self-time functions (percent), at most
five lines each. Re-profile a run when a landed change touched its hot path.

- mode_switch.dual-r, `ttklzlkq` (profiling build, 199 Hz): factor update 30.10%, LU nucleus 15.98%, BTRAN 12.49%, fresh-factor list filtering 6.31%, factor_structure 2.93% self cycles. Full solve: 11393 iterations, 17 refactorizations. Other initial profiles still needed: degen3.dual-r, degen3.primal-c, quad_2717.dual-r.

## In progress

- mode_switch.dual-r: enumerate leaving-row references from sorted, deduplicated rowl[r] instead of scanning every column twice; keep all per-column work estimates and their order.

## Candidates

Ordered by expected payoff; one line each, with the bench run it targets. Start from PLAN.md "Where to
look first" and replace these with measured items after the first profiles.

- Unresolved regressions against the baseline (PLAN.md "Measuring"; resolve before the next change
  lands): fx_quad_obj_5371.dual-r 0.334 -> 0.726 s, primal_5959.primal-r 0.563 -> 0.841 s,
  degen3.barrier-r 0.608 -> 0.857 s. Pair-rerun with instruction counts; bisect the landed commits if real.
- Goal setup (PLAN.md "Goal"): solve-time column in `perf/bench.sh`, missing `ref` cells via
  `ref_run.py`, then `current/ref` from solve times.
- mode_switch.dual-r: factor update (30.10%) still scans full A/B lists for membership, references and row-capacity accounting; annotate the loops before choosing a change.
- mode_switch.dual-r: BTRAN (12.49%) still allocates dense scratch vectors and recomputes L levels; annotate before choosing a change.
- mode_switch.dual-r: LU nucleus (15.98%); profile its loops before choosing a change.
- mode_switch.dual-r: fresh-factor row-list filtering (6.31%); initialization filters a full column-order list for each nonempty row, including an all-position list for empty original U rows.
- quad runs: cost of `crate::f80` arithmetic in the quad kernels.
- all: `dump::on` and record building on hot paths with no dump selected.
- all: per-iteration allocation in FTRAN/BTRAN, pricing and ratio tests.

## Tried, no gain

- Read eliminate work counts directly from maintained rowl lengths instead of rebuilding counts: parity full 987/987 and release tests pass, but mode_switch.dual-r 22.849 -> 21.851 s is smaller than the old 1.849 s spread (instructions 146886.1 -> 140945.4M, -4.0%). Alternating DEGEN2 primal 0.092/0.090 then 0.091/0.091 s; dual 0.068/0.067 s, marginal against spread. Reverted.

## Tried, breaks parity

(idea, first differing record from `rose-parity check`, why)
