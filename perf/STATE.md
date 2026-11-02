# Performance state

Current state of the work in [PLAN.md](PLAN.md). Rewrite in place; no log (history is `jj log`).

## Numbers

Baseline commit: `vlsokskq` (rose-parity gate). Current commit: same (no optimization landed yet).
Parity: full tier 987/987 at the baseline.

Seconds, `perf/bench.sh` medians (release, no dumps); `ref` is the reference's `result.runtime` from
`../reverse-engineer/runs/` (`-`: no reference run with a runtime). Fill `baseline` once, in the first
session, then keep it fixed; update `current` and `instr` (millions, `--instr`) after each landed change.

| run | ref | baseline | current | current/ref | instr |
|---|---|---|---|---|---|
| degen2.primal-c | 0.013 | | | | |
| degen2.dual-r | 0.010 | | | | |
| degen3.primal-c | 0.163 | | | | |
| degen3.primal-r | 0.194 | | | | |
| perold.primal-c | 0.161 | | | | |
| pilot4.primal-r | 0.016 | | | | |
| 80bau3b.primal-c | 0.164 | | | | |
| primal_5959.primal-r | - | | | | |
| degen3.dual-c | 0.083 | | | | |
| degen3.dual-r | 0.083 | | | | |
| perold.dual-r | 0.039 | | | | |
| 80bau3b.dual-r | 0.134 | | | | |
| greenbea.dual-c | - | | | | |
| greenbea.dual-r | - | | | | |
| mode_switch.dual-r | 0.4 (reverse-engineer STATE.md) | | | | |
| primal_5959.dual-r | - | | | | |
| quad_2717.dual-r | - | | | | |
| fx_quad_obj_5371.dual-r | - | | | | |
| degen3.barrier-r | 0.075 | | | | |
| 80bau3b.barrier-r | - | | | | |
| greenbea.barrier-r | - | | | | |
| 80bau3b.pdhg | - | | | | |
| greenbea.pdhg | - | | | | |

## Hotspots

Per profiled bench run: the commit it was taken at and the top self-time functions (percent), at most
five lines each. Re-profile a run when a landed change touched its hot path.

(none yet: profile mode_switch.dual-r, degen3.dual-r, degen3.primal-c, quad_2717.dual-r first)

## In progress

(none)

## Candidates

Ordered by expected payoff; one line each, with the bench run it targets. Start from PLAN.md "Where to
look first" and replace these with measured items after the first profiles.

- mode_switch.dual-r: find the per-iteration O(m)/O(n) work behind the ~50x ratio.
- quad runs: cost of `crate::f80` arithmetic in the quad kernels.
- all: `dump::on` and record building on hot paths with no dump selected.
- all: per-iteration allocation in FTRAN/BTRAN, pricing and ratio tests.

## Tried, no gain

(idea, bench run, measured result)

## Tried, breaks parity

(idea, first differing record from `rose-parity check`, why)
