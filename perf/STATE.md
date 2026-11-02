# Performance state

Current state of the work in [PLAN.md](PLAN.md). Rewrite in place; no log (history is `jj log`).

## Numbers

Baseline solver: `8b645007` (same code as `vlsokskq`, rose-parity gate). Current solver: `yyqrnssy` (L transpose cached per factorization).
Parity: full tier 987/987; release tests pass. Baseline measurements: three timed runs, one warm-up, hyperfine 2.0.0; instruction counts from perf stat. Wall time is variable on this host; use paired reruns before accepting gains.

Seconds, `perf/bench.sh` medians (release, no dumps); `ref` is the reference's `result.runtime` from
`../reverse-engineer/runs/` (`-`: no reference run with a runtime). The baseline stays fixed. Current timings/instructions use the full sweep, superseded by the latest paired rerun where available; instruction counts are millions (`--instr`).

Acceptance pairs (old/new, three timed runs): DEGEN2 dual 0.125/0.106 s, 621.5/506.2M instructions; PEROLD dual 0.357/0.307 s, 3027.5/2577.5M; 80BAU3B dual 1.706/1.228 s, 9995.3/7290.1M. Mode-switch 25.838/17.926 s, 227581.9/146867.1M: instruction reduction 35.5%, but wall-time gain is smaller than the baseline 10.069 s spread.
CPU-0 alternating old/new reruns resolved apparent slowdowns in the full sweep: DEGEN3 barrier 0.621/0.631 s (within spread), 80BAU3B barrier 0.386/0.383 s, GREENBEA barrier 2.317/2.292 s, quad-objective 0.639/0.337 s (wide old spread). Barrier/PDHG instructions remain unchanged within measurement noise.

| run | ref | baseline | current | current/ref | instr |
|---|---|---|---|---|---|
| degen2.primal-c | 0.013 | 0.182 | 0.140 | 10.8 | 727.0 |
| degen2.dual-r | 0.010 | 0.125 | 0.105 | 10.5 | 506.2 |
| degen3.primal-c | 0.163 | 3.237 | 2.430 | 14.9 | 12436.7 |
| degen3.primal-r | 0.194 | 3.872 | 1.816 | 9.4 | 15501.0 |
| perold.primal-c | 0.161 | 2.430 | 1.294 | 8.0 | 11234.8 |
| pilot4.primal-r | 0.016 | 0.253 | 0.132 | 8.2 | 1245.4 |
| 80bau3b.primal-c | 0.164 | 1.618 | 1.148 | 7.0 | 12035.6 |
| primal_5959.primal-r | - | 0.563 | 0.561 | - | 5141.0 |
| degen3.dual-c | 0.083 | 1.002 | 0.753 | 9.1 | 5992.1 |
| degen3.dual-r | 0.083 | 1.008 | 0.763 | 9.2 | 6002.3 |
| perold.dual-r | 0.039 | 0.355 | 0.307 | 7.9 | 2577.5 |
| 80bau3b.dual-r | 0.134 | 1.103 | 1.228 | 9.2 | 7290.1 |
| greenbea.dual-c | - | 2.870 | 2.781 | - | 14436.1 |
| greenbea.dual-r | - | 3.598 | 2.755 | - | 14465.9 |
| mode_switch.dual-r | 0.4 (reverse-engineer STATE.md) | 25.815 | 17.926 | 44.8 | 146867.1 |
| primal_5959.dual-r | - | 1.300 | 1.336 | - | 7324.6 |
| quad_2717.dual-r | - | 1.130 | 1.099 | - | 6265.5 |
| fx_quad_obj_5371.dual-r | - | 0.334 | 0.337 | - | 3078.5 |
| degen3.barrier-r | 0.075 | 0.608 | 0.631 | 8.4 | 5498.1 |
| 80bau3b.barrier-r | - | 0.379 | 0.383 | - | 4057.4 |
| greenbea.barrier-r | - | 1.569 | 2.292 | - | 16053.3 |
| 80bau3b.pdhg | - | 4.552 | 4.492 | - | 92390.1 |
| greenbea.pdhg | - | 2.801 | 1.900 | - | 37741.2 |

## Hotspots

Per profiled bench run: the commit it was taken at and the top self-time functions (percent), at most
five lines each. Re-profile a run when a landed change touched its hot path.

- mode_switch.dual-r, `yyqrnssy` (profiling build, 199 Hz): factor update 24.28%, BTRAN 17.31%, LU nucleus 9.94%, fresh-factor list filtering 3.64%, free 3.27% self cycles. Full solve: 11393 iterations, 17 refactorizations. Other initial profiles still needed: degen3.dual-r, degen3.primal-c, quad_2717.dual-r.

## In progress

(none)

## Candidates

Ordered by expected payoff; one line each, with the bench run it targets. Start from PLAN.md "Where to
look first" and replace these with measured items after the first profiles.

- mode_switch.dual-r: factor update (24.28%); eliminate rebuilds numeric U rows in CSR on every update (two full U scans). Investigate a reusable row index with exact column order and invalidation; the row-count-only shortcut is rejected below.
- mode_switch.dual-r: BTRAN (17.31%) still rebuilds the U transpose each call (U changes at each pivot); investigate reuse with invalidation.
- mode_switch.dual-r: LU nucleus (9.94%); profile its loops before choosing a change.
- quad runs: cost of `crate::f80` arithmetic in the quad kernels.
- all: `dump::on` and record building on hot paths with no dump selected.
- all: per-iteration allocation in FTRAN/BTRAN, pricing and ratio tests.

## Tried, no gain

- Read eliminate work counts directly from maintained rowl lengths instead of rebuilding counts: parity full 987/987 and release tests pass, but mode_switch.dual-r 22.849 -> 21.851 s is smaller than the old 1.849 s spread (instructions 146886.1 -> 140945.4M, -4.0%). Alternating DEGEN2 primal 0.092/0.090 then 0.091/0.091 s; dual 0.068/0.067 s, marginal against spread. Reverted.

## Tried, breaks parity

(idea, first differing record from `rose-parity check`, why)
