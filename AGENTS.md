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
