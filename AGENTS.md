# Project instructions

This directory is a separate Jujutsu repository. Use `jj` for all version
control operations here, running commands from this directory.

## Testing

Always use the release build when testing: run `cargo test --release` (and
`cargo build --release` for the CLI/examples). Debug builds are too slow for the
Netlib-based tests and can appear to hang on cycling or large models.
