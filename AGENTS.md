# Project instructions

This directory is a separate Jujutsu repository. Use `jj` for all version
control operations here, running commands from this directory.

## Testing

Always use the release build when testing: run `cargo test --release` (and
`cargo build --release` for the CLI/examples). Debug builds are too slow for the
Netlib-based tests and can appear to hang on cycling or large models.

## Benchmarks and debugging

Prefer quick benchmarks and small problems when benchmarking or debugging. Use the
small Netlib models (AFIRO, SC50A, BLEND, ADLITTLE, BOEING2, BORE3D) and cap
iterations (`--max-iterations`) instead of running slow ones; FIT2P takes minutes
and is opt-in (`tools/bench_rose.py ... --with-fit2p` in `reverse-engineer/`).
Reproduce a problem on the smallest model that shows it before using a large one.

## Bit-for-bit parity during performance work

Rose agrees with the reference bit for bit on the runs of `parity/corpus.txt`. A change that is meant to
leave the numbers alone (performance work, refactors) must keep every dump record and every result
identical to the frozen baseline:

    cargo build --release && target/release/rose-parity check            # fast tier, 940 runs, ~2 s
    target/release/rose-parity check --tier full                         # at commit points, 987 runs, ~40 s

`rose-parity` runs each corpus line with `ROSE_DUMP=all`, hashes every record per tag, and compares with
`parity/baseline/`. On a difference it reruns the frozen baseline binary (`target/parity-baseline/rose`)
next to the new one and prints the first differing record and field. `--only ID,..` selects runs.
`target/` is not versioned: on a fresh checkout, build the commit named in `parity/baseline/META.json`
and copy its `rose` there (`check` itself needs only the digests).

- Re-record the baseline (`rose-parity record --tier full --source "<commit>"`) only from a build that
  passes `check` and the reference comparisons, i.e. after a change that deliberately moves rose closer
  to the reference. A performance change never re-records. A new dump kind adds records: check with
  `--ignore NEWTAG`, then re-record with the new build.
- The work counter (`crate::work`) decides when to refactorize, so its increments are part of the
  numbers. A faster kernel must book the same estimates as the reference's loops, not as its own: when
  a kernel is rewritten, compute its estimate in a separate function of the inputs (or keep the old
  counting loop) and leave the WORK records unchanged.
- Rewriting a numeric kernel: first copy the current version into a `#[cfg(test)]` reference module
  next to it and add a differential test (old against new, bit for bit) on random inputs and on inputs
  taken from a dump, then change the copy that runs. The parity gate is the backstop, not the first check.
- Do not add `target-cpu`, fast-math style flags or `mul_add`, and do not use `f64::powi` (use
  `real::pow2` or explicit products). The toolchain is pinned in `rust-toolchain.toml`; move it only
  with a full parity check. `tests/libm_pins.rs` pins the system libm's results: when it fails, run the
  full parity check before trusting the build.
- Profile with `cargo build --profile profiling` (release code with line tables).

## Performance work

The plan and working loop are in `perf/PLAN.md`, the current numbers and tried ideas in `perf/STATE.md`
(read it first; rewrite it in place, no log). Benchmarks: `perf/bench.sh` over `perf/bench.txt`.
