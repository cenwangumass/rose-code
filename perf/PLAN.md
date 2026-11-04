# Performance plan

Make rose faster without changing any number it computes. Rose agrees with the reference (Gurobi 13.0.3)
bit for bit on the parity corpus; every speedup must keep that agreement exactly. This file is the
standing plan and the working loop. The current numbers, hotspots and the ideas already tried are in
[STATE.md](STATE.md); read it first in every session, and keep it current (rules at the end).

## Goal

Rose at least as fast as the reference: `current/ref` ≤ 1.0 on every bench run, under the invariant
below. Prioritize measured removable work across ordinary dual, ordinary primal, extended-precision and barrier paths. Historical ratios identify gaps but do not rank opportunities reliably across host timing bands.

- **Like with like.** `ref` is Gurobi's `Runtime` (solve only, single thread, same host: AMD Ryzen 9
  5900X). Compare it with rose's solve time, `diagnostics.seconds` in the result JSON (median over the
  timed runs), not with the process wall time of `bench.sh`, which adds start-up and MPS parsing. The
  harness reports both wall and solve times; STATE.md uses solve time for `current/ref`.
- **Missing `ref` cells:** STATE.md records which original-library runtimes are reachable and which
  fail the available license. If more become reachable, run them with
  `../reverse-engineer/tools/ref_run.py` (see its docstring; the profile names match bench.txt's: `release`, `dual-release`, `barrier-release`, `pdhg-release`, ...)
  through gurobipy with the unpatched library, three times, and take the median `Runtime`. Records made by
  `build/native_run` or the patched library have no usable runtime. Write the medians into STATE.md and
  note the run files used; do not overwrite existing files in `../reverse-engineer/runs/` (write the new
  ones to a scratch directory).
- **Done** when every bench run with a `ref` has `current/ref` ≤ 1.0, or the remaining gap of a run is
  documented in STATE.md as needing something the invariant forbids (name the rule and the measured
  share of the time it accounts for). A run without a reachable reference time does not block the goal.

## The invariant

A change is acceptable only if `target/release/rose-parity check --tier full` passes: every dump
record, every result vector, every status and exit code of all 987 corpus runs equals the frozen
baseline. No exceptions, no "within a ULP", no re-recording.

- Never edit `parity/baseline/` or re-record it; never pass `--ignore` to make a check pass. If you add a
  dump kind or a corpus run, record only the new runs, and only with the frozen binary:
  `rose-parity record --only NEWID --bin target/parity-baseline/rose` (for a new dump kind, check with
  `--ignore NEWTAG` instead; do not add dump kinds as part of an optimization).
- The dump records are part of the behaviour. A data-structure change that alters a record that only
  describes rose's own layout (pool fill, capacities, slots) still fails the gate: those records mirror
  the reference's layout. Keep the records, compute them some other way if needed.
- The work counter (`crate::work`, WORK records) decides when to refactorize. A faster loop must book
  exactly the same estimates (same sites, same values, same order). When the estimate is counted inside
  a loop you want to shorten or skip, first move the count into a separate function of the inputs that
  reproduces it, check, then speed up the loop.
- Floating point: keep every operation and its order. No reassociation, no `mul_add`/FMA, no
  `target-cpu` or fast-math flags, no `f64::powi` (use `real::pow2` or products), no change of a sum's
  order (including SIMD lane order: `Real` documents that the double kernels add in vector lanes and the
  x87 ones one term after the other; keep each kernel's existing order). Do not swap libm calls for
  approximations. Do not move the toolchain pin.
- No threads, no nondeterminism (hash map iteration, time, addresses) in anything that reaches a
  decision or a record.

## Tools

- `cargo build --release`: the binaries (`target/release/rose`, `target/release/rose-parity`).
- `python3 perf/build_pgo.py --check --test`: opt-in PGO release CLI, trained on 39 bounded small-model
  runs across all methods and recovery exits. Uses the pinned compiler and matching `llvm-profdata`;
  override the latter with `LLVM_PROFDATA`. Generated profiles, training log and exact build flags stay
  in `target/perf-pgo/`; no training data is committed. The current performance table uses this build.
  Ordinary Cargo builds remain the ThinLTO baseline and can overwrite the optimized CLI; rebuild with
  the script before comparing PGO candidates. Save the preceding binary before rebuilding as usual.
  The script uses controlled Rust flags for both stages; it does not inherit `RUSTFLAGS`.
  Both stages request 64-byte loop alignment (`-Cllvm-args=--align-loops=64`) to keep hot loop heads
  stable as nearby code changes.
  `--system-tcmalloc` additionally enables the Linux CLI's optional system allocator feature
  in both stages and release tests; gperftools development/runtime libraries must be installed.
  The selected Cargo features are recorded in `build.json`; reproduce them alongside saved flags.
  Training also clears `LLVM_PROFILE_FILE` to retain rustc's default shared online-merge mode.
  Per-process raw profiles regressed barrier benchmarks in a same-source controlled comparison.

- `cargo build --profile profiling`: the same code with line tables (`target/profiling/rose`), for perf.
  It produces the same bits as release.
- `target/release/rose-parity check [--tier full] [--only ID,..]`: the gate. Fast tier (940 runs) takes
  seconds, full (987) under a minute. On a difference it prints the first differing record and field.
- `perf/bench.sh [--only ID,..] [--runs N] [--instr] [--bin PATH]`: times the runs of
  [bench.txt](bench.txt) with hyperfine (median of N after a warm-up); `--instr` adds user-space
  instruction counts (`perf stat`), which have far less noise than wall time.
- For profiling a PGO candidate, reproduce its `features` and `rustflags` recorded in
  `target/perf-pgo/build.json` in addition to the frame-pointer flag below, then parity-check that
  profiling binary. For example, pass `--features system-tcmalloc` when selected. Without the saved
  profile-use flag the profiling build describes the ordinary ThinLTO baseline.
- For reliable inclusive profiles, build with `RUSTFLAGS='-C force-frame-pointers=yes' cargo build
  --profile profiling`, then parity-check that binary with `--tier full --bin target/profiling/rose`.
  Capture with `perf record -F 1999 -g --call-graph fp -o target/X.data -- target/profiling/rose solve
  MODEL ARGS --output /dev/null`. Use `perf report -i target/X.data --no-inline --children
  --call-graph none --stdio` for inclusive shares, `--no-children` for self time. The default inline
  expansion duplicates frames here; DWARF unwinding also loses callers. Use `perf annotate` for hot
  functions. Bound the run with a bench model or `--max-iterations N`. `perf stat -e
  instructions:u,cycles:u,cache-misses:u,branch-misses:u` gives a quick counter picture.
- For instruction-family self profiles, use the accepted release CLI directly: self attribution does
  not need frame pointers, which can change instruction counts. Keep the frame-pointer build for
  caller ancestry. Audit lost/throttled samples and represented instructions against direct counters.
  Keep wrapper result parsing/validation outside the measured process, and separate solver samples
  from wrapper samples.
- `python3 perf/paired.py --before BASE --after CANDIDATE --run ID --pairs 32 --seed N
  --cpu C --max-iterations CAP --output target/NAME.json`: bounded randomized AB/BA
  solve-time screening with interleaved same-binary AA controls. The schedule, binary
  hashes and decision rule are saved before warm-up. Complete result objects and exit
  codes must agree; raw outputs and every timing sample are retained. Existing output
  paths are refused. A positive screen requires sign p <= .01, unbiased AA controls,
  at least 0.5% median improvement, a conservative 95% paired-median bound beyond
  the AA uncertainty, and improvement in both execution orders. It is a screening
  tool, not a replacement for parity, independent confirmation or the full benchmark.
- `cargo test --release` (about 35 s): the reference-derived tests. Run before each commit.
- Format with `RUSTFMT=~/.local/bin/rustfmt-custom cargo fmt -- --config blank_line_between_items=true`
  (never plain `cargo fmt`). Version control is `jj`, run from `rose/`.

## Measuring

- A probe built with a preceding source's PGO profile is preliminary. For a changed
  representation or constructor, verify important codegen conclusions with the matching
  fixed 39-case profile; preliminary counts can change even outside the changed runtime
  branch. Preserve both builds and measurements, and freeze final source/tests before
  official timing screens.
- Times are wall-clock medians from `perf/bench.sh` on the release build, no dumps. Measure the
  baseline and the candidate in the same session, back to back (`--bin` lets you time a saved copy of the
  previous binary, e.g. `cp target/release/rose /tmp/rose-before` before rebuilding).
- For grouped hyperfine comparisons, require a gain larger than the run's min-max
  spread and supporting instruction counts. When host shifts make this inconclusive,
  use the randomized paired screen above with predeclared seeds, bounds and sample
  counts. Require two independent positive screens before treating a small gain as
  established; retain every failed/inconclusive screen too. Stop an inconclusive
  session instead of extending it until it passes. An AA control showing bias or
  uncertainty comparable to the proposed gain makes the session inconclusive.
  For a change below 3% on its target, require the instruction count to drop.
  Full parity, release tests, complete no-dump result validation and the full paired
  benchmark remain required; resolve any unrelated regression before accepting. A
  reversed-role slowdown passing significance, AA health, conservative noise
  separation and both-order checks is a regression even below 0.5%; the minimum
  gain threshold applies to accepting speedups.
- A change must not slow any other bench run by more than the noise. Run the whole bench set before
  committing; time only the affected runs while iterating.
- Regressions against the baseline column: a run whose `current` is slower than `baseline` by more than
  its spread must be resolved before the next change lands. Rerun it paired (baseline binary and current
  binary alternately, same CPU, `taskset -c N`), with instruction counts. If it is real, find the landed
  change that caused it (time the intermediate commits) and fix or revert that change; if it is noise,
  replace the table entry with the paired numbers. Keep unresolved regressions at the top of STATE.md's
  candidates.
- Keep each measurement bounded: the slowest bench runs take tens of seconds; use `--only` and
  `--runs 3` while iterating. Run anything expected to take more than a few minutes in the background.

## The loop

Repeat until the goal above is met. Running out of candidates does not end the work: re-profile the run
with the largest ratio and look one level deeper (inclusive time of the callers, `perf annotate` of the
hot loops, cache and branch misses) before declaring a gap forbidden by the invariant.

1. **Pick a target.** Resolve verified regressions first. Then investigate factor production and
   installation ownership, followed by FTRAN representation boundaries, barrier execution planning,
   remaining whole extended-precision kernels and primal kernel specialization. Start with
   [factor-audit.md](factor-audit.md). Use the smallest witness; refresh only the attribution needed
   to estimate removable cost before implementing. Keep one ordinary dual, ordinary primal, quad
   and barrier representative. A historical large ratio alone does not justify another small trial.
2. **Understand it.** Read the hot function and its callers. Find out why it is slow: repeated full
   scans where the reference keeps lists, allocation per iteration, dense loops over sparse data,
   recomputation of invariant values, F80 software arithmetic in a loop that only needs a few values,
   string/JSON work on the hot path when no dump is selected, bounds checks in tight loops, cache-hostile
   layout. For a result cache, measure exact ordered input working sets before choosing its bound;
   a small cache can miss a highly repetitive cycle. Write the hypothesis into STATE.md "In progress" (one line).
3. **Protect it.** If the change rewrites a numeric kernel, first copy the current version into a
   `#[cfg(test)]` module next to it and add a differential test: old against new, bit for bit
   (`to_bits()`), on random inputs (include zeros, signed zeros, tiny and huge magnitudes, empty and
   one-element cases, the F80 type where the kernel is generic) and on inputs taken from a dump of a
   real run. If the change only restructures control flow or memory, the parity gate is the test.
4. **Change one thing.** Keep the diff to one idea. No drive-by refactors in the same change.
5. **Check.** `cargo build --release`, `rose-parity check` (fast); fix or revert on any difference.
   The first differing record tells you where behaviour moved; do not guess around it.
6. **Measure** as above. If the gain is not real, revert (`jj restore` / `jj abandon`) and record the idea
   under "Tried, no gain" with the measured number.
7. **Commit.** `rose-parity check --tier full`, `cargo test --release`, format, then one jj commit per
   optimization. Message: what changed and why it is faster, the bench numbers before/after for the
   affected runs, and "parity: full tier 987/987". Start the next change on a new empty commit (`jj new`).
8. **Update STATE.md** (rewrite the affected lines, do not append a log).

If a change fails the gate and the cause is not found within a few attempts, revert it, record the idea
under "Tried, breaks parity" with the first differing record, and move on. Never commit a failing gate,
never commit with debug prints, never leave the tree with uncommitted half-work at the end of a session
(either commit a passing change or restore).

## Where to look first

These are leads to verify with a profile, not facts. Remove each from STATE.md's candidates once it is
measured.

- **The mode-switch witness** (`mode_switch.dual-r`, 15000 rows): a large absolute-time run whose
  original-library reference currently fails the available license. Work on it after reachable ratios.
  Check per-iteration scans, repeated candidate/weight rebuilds and dense work on sparse data; use the
  current profile and instruction count in STATE.md rather than the old, unusable reference timing.
- **Large dual and primal runs** (DEGEN3, PEROLD, 80BAU3B, GREENBEA): pricing, ratio tests,
  FTRAN/BTRAN (`ftran_order.rs`), the LU update and refactorization (`lu*.rs`), the bound-flip and
  update loops in `dual.rs` / `simplex.rs`.
- **The quad passes** (`quad_2717`, `fx_quad_obj_5371`, `primal_5959.dual-r`): `crate::f80` is a
  software x87 extended-precision type; its add/mul/div are in every quad kernel. Bounded-normal
  to-double bit assembly is now qualified with the complete original conversion
  and all MXCSR modes as oracles; preserve its original boundary fallback.
  Eight-entry pure normal-Devex increment reuse is qualified per compact owner;
  full extended input keys and lazy column fills preserve live mutable state.
  Reset/owner-aware measurements justify that bound on small problems; larger
  problems and mode 2 retain the original path. The lazy mode-2 coefficient
  trial remains unapplied after an inconclusive full-session control.
  Refresh remaining BTRAN, row/control and software arithmetic after these changes.
  Long-only register eta dots pass exact differential checks, including all MXCSR
  modes, but miss conservative timing-noise separation; do not extend that screen.
  Current post-Devex-cache attribution still puts extended BTRAN first. Bounded
  stack scratch is verified but has no useful counter gain. Before another
  division fast path, isolate its actual call sites: dense U already inlines
  zero-numerator handling, and global zero counts do not imply removable ABI work.
  Global and per-column quotient-result tables lose at counters. Exact reciprocal
  parameters pass arithmetic and complete-BTRAN oracles, but shared four-entry
  tables do not reduce target instructions. Preparing parameters on diagonal
  writes removes lookup/locking and lowers instructions/cycles; its controlled
  target screen remains inconclusive. Do not extend or repeat that screen.
  Current accepted ordinary-primal PILOT4 self attribution puts solve_primal,
  BTRAN, FTRAN and primal row updates first. Packed pair/run products lose at
  counters. Exact ordered numeric streams have negligible full-result reuse;
  repeated row orders cover about half the terms, but the tested eight-entry
  gather-plan stream costs3.09 times the original standalone computation once
  lookup/construction are included. Do not implement those caches from hit rates
  alone. The ordinary positive-floor rewrite lowers instructions but its fixed
  controlled screen is inconclusive; do not extend it or use its median as a gain.
  Packed32-bit immutable ordinary L entries preserve numeric/WORK behavior but
  lose at matching-PGO counters; do not treat reduced stream bytes as a gain.
  Two-word F80 software layout now has independent full-field equivalence and
  real register-ABI evidence. Direct sign-word operations improve target counters,
  but the fixed controlled screen is inconclusive; do not repeat or extend it.
  The plain layout grows Option<F80>16→24 and optional pairs 32→40. Tagged
  NonZeroU64 software metadata independently preserves original optional sizes,
  all full-field arithmetic and zero/None semantics, with a safe tag invariant.
  Matching-PGO fx/Quad instructions/cycles show no useful gain; stopped before
  timing seeds 29217/29289. Both layouts remain unapplied. Do not infer solver
  gains from ABI, optional sizes or isolated arithmetic/initialization probes.
  Exact auxiliary-constructor streams show no repeats on the ordinary small
  witnesses or bounded PILOT4/PEROLD, but eight entries capture 1760/3006 fx/Quad
  repeats. Full key/lookup/copy feasibility is favorable; matching-PGO integration
  reduces target instructions. Fixed fx screen 30217 remains inconclusive under
  conservative AA separation; do not repeat or extend it. Deep-copy cache remains
  unapplied. Fresh auxiliary numeric fields are read-only after construction;
  immutable-base ownership sharing would be a distinct, unimplemented experiment.
  Full FTRAN/BTRAN dependency and version-aware input measurement is complete
  across 18 bounded cases: 132,997 nonfresh calls have only two repeated full-FTRAN
  keys and no needed/BTRAN repeats. Do not implement a general mutable-version
  cache from these streams. Fresh BTRAN overlaps the prior unqualified cache;
  do not repeat its failed screen. Fresh full FTRAN is a distinct lead: two
  entries capture 1730/2883 fx/Quad and 5273/5627 primal5959 repeats, predominantly
  pivot-column solves. Full lookup/effect replay/copy component ratios are
  0.111/0.112/0.199/0.198 with exact captured outputs and written state. The
  matching-PGO integration preserves the complete original kernel, replaying
  only written fields and keeping wrappers/caller WORK live. All focused and
  per-hit full-field oracle checks pass; fx/Quad instructions fall 4.021/1.621%.
  Fixed fx screen 31217 is inconclusive, including opposite order medians; do
  not repeat or extend it. This cache remains unapplied. Next pursue larger
  factor-production/installation ownership savings, with exact numerical and
  allocation evidence. Auxiliary immutable-base sharing remains unimplemented;
  general producer/pre-overlay numerical identity is now disproven. The paired
  SC50A k1 supplied subset matches, as do 880 small/DEGEN2 supplied values, but a
  deterministic 51x51 dense handoff differs by one ULP before the mandatory
  overwrite. Do not treat sparse matches as a global proof or bypass the overlay.
  The strictly-lower-L contiguous-suffix dependency cut is now numerically
  verified with an independent complete constructor, dense mismatch, last-write
  duplicates/extremes/fallbacks and all SSE modes. Matching fixed39 PGO reduces
  PILOT4 instructions 1.068%, but screen32217 is inconclusive; stop without
  seed32289 or full candidate gates. It remains unapplied. Do not repeat either
  this screen or the earlier arithmetic-only guarded-replay screen.
  The extended BTRAN conversion/operation-boundary audit is complete. Prepacking
  alone does not make the small component faster than software; rebuilding its
  cache per call is worse. Combined multiply/store/subtract blocks retain both
  roundings and improve the component, but matching-PGO integration raises
  fx/Quad cycles5.75/4.84% despite fewer instructions; rejected before33217/33289.
  Register-product range certificates pass1M raw triples but add only0.44%
  component improvement, insufficient evidence for another integration.
  Next quantify dead native divisions: its dense U loop executes x87 division
  before dropping zero/tiny numerators, whereas accepted software inlines its
  zero return. Protect complete output/WORK/state against the original BTRAN before testing a predicate-first native path; use bounded small
  kernels and actual call-site counts first. Do not infer the native opportunity
  from whole-solver zero-division totals or repeat coefficient-cache/register
  variants without new cost evidence. Auxiliary immutable-base sharing remains
  an unimplemented lead; further fresh-factor microvariants need measured costs.
  Comparator alignment stays centralized; normal lunum's build/cache transport
  still needs an explicit no-build/schema-current source option for direct reuse.
  Look for a faster exact implementation of the same rounding (the F80 tests and `../reverse-engineer/evidence/f80-subnormal-reference.json`
  define it), and for F80 work that runs where its results are not used.
- **Dump overhead with no dump selected**: about 70 `dump::on("kind")` sites in `dual.rs` and
  `simplex.rs`. Each call scans the selected kinds by string; check whether any of them sits in a
  per-iteration or per-element loop, and whether any record is built before the `on` test.
- **Allocation**: current two-period frame-pointer ancestry puts roughly 45–47%
  of PEROLD allocator self samples under refactorization and 29–30% under iteration.
  Direct nucleus input construction and dense-handoff extraction are rejected at
  counters. Chain-maintenance counts do not yet justify an intrusive-list rewrite.
  Scalar slot-only U-order/elimination heap keys are verified equivalent but
  timing-inconclusive; do not extend that screen. Numeric-LU step queues were
  tested separately and rejected. Prioritize a whole remaining numeric/control
  pass over further small heap or constructor buffer variants.
  Inspect the allocating callee: an observed ancestor return address can skip it.
  Profiling-build shares do not bound accepted release cost.
- **Barrier** (`chol.rs`, `barrier.rs`, `sgemm32.rs`) and **PDHG** (`pdhg.rs`):
  the predecessor release refresh puts numeric factorization31.24%, DGEMM14.82%,
  DTRSM10.58%, and ordering/partition/integer-sort14.75% of Rose instruction
  self samples. Fixed child-rounding dispatch is rejected at matching-PGO counters.
  Unique AA-neighbour construction is qualified with complete ordered-pattern
  equality and two independent positive timing screens. Quotient-order adjacency
  buffer reuse and completed-hash store lifetime are also qualified; current
  native hash loops retain the accumulator in registers with unchanged additions.
  Return to the larger ordinary dual/primal gaps: use the current accepted
  PEROLD/PILOT4 instruction profiles above before choosing another whole pass. Remaining
  partition/numeric-factor annotations predate these changes and are diagnostic.

## Out of scope

- Any change of results, iteration counts, statuses or dump records; parameters or defaults.
- Algorithmic improvements beyond the reference (better pricing, presolve, other orderings): they
  change behaviour.
- Multi-threading, GPU, unsafe SIMD that changes the lane order.
- Debug-build speed; FIT2P (minutes; profile it only with `--max-iterations`).

## STATE.md rules

STATE.md holds only the current state: the commit the numbers belong to, the bench table (baseline,
current, reference), the current hotspots per bench run, the candidate list, the idea in progress, and
the ideas tried and rejected (one line each with the measured reason, so they are not retried). It is not
a log: rewrite lines in place, delete what is settled (a landed change leaves STATE.md through the
updated numbers and hotspots; its history is in `jj log`). Keep it short enough to read in one go.
