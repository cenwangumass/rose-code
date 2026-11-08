# Performance plan

Make rose faster without changing any number it computes. Rose agrees with the reference (Gurobi 13.0.3)
bit for bit on the parity corpus; every speedup must keep that agreement exactly. This file is the
standing plan and the working loop. The current numbers, hotspots and the ideas already tried are in
[STATE.md](STATE.md); read it first in every session, and keep it current (rules at the end).

## Goal

Make rose faster on every fixed row of `bench.txt`. Compare solve times on the same machine and logical CPU. The reference is the supplied reference reconstructed at `target/perf-reference/reconstructed/libgurobi130.so`, SHA-256 `26a505171314b3300b754b4eff5f5719366b47ccf6f3d59770df99e3f256ee30`. Every row has a measured reference, including 80BAU3B, GREENBEA and the mode-switch witness.

Use Gurobi's `Runtime` with the single-threaded profile parameters in `ref_run.py` and the iteration limit supplied by the bench row. Compare with rose's `diagnostics.seconds`. Take the median of three measurements after one warm-up per solver. Preserve each result, exit code, timing and the library hash. Keep measurements from different hosts, allocators and build configurations labeled separately. The current-host runner is `perf/compare_reference.py`; full results and provenance are linked from STATE.md.

The goal is met when rose's median solve time is strictly below Gurobi's on every row and full numerical parity passes. Keep any measured gaps and work constrained by the numerical invariant in STATE.md.

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
  in `target/perf-pgo/`; no training data is committed. Historical performance tables used this build; the current-host configuration is explicit in STATE.md.
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
- `python3 perf/compare_reference.py --library target/perf-reference/reconstructed/libgurobi130.so --output target/NAME`:
  full fixed 23-row current-host comparison, one warmup/three samples per solver,
  CPU 2 and external 70s process deadlines. It records Gurobi Runtime, effective
  parameters, complete raw outputs/exits and stable repeated results. Uses the
  qualified `../reverse-engineer/tools/ref_run.py` native harness; no gurobipy
  install is required. Descriptive medians do not qualify optimization commits.
- `cargo test --release` (about 35 s): the reference-derived tests. Run before each commit.
- Format with `RUSTFMT=~/.local/bin/rustfmt-custom cargo fmt -- --config blank_line_between_items=true`
  (never plain `cargo fmt`). Version control is `jj`, run from `rose/`.

## Measuring

- Put temporary diagnostic probes in an isolated source copy/workspace with a
  separate Cargo target directory. Preserve source hashes and fixtures, and
  validate complete results. This avoids invalidating the accepted build and
  relinking all integration tests after restoring diagnostic edits. The isolated
  baseline for score reuse is byte-identical to the accepted CLI; retain that
  check when using copied builds for comparisons.
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

- **The mode-switch witness** (`mode_switch.dual-r`, 15,000 rows) has the largest relative gap on the measured EPYC. Start from the current FTRAN/BTRAN sparse-workspace profile and large-buffer request counts in STATE.md.
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
  unapplied. Fresh auxiliary numeric fields are read-only after construction,
  but the current owner already removes most of those historical repeats;
  see the current constructor census before proposing shared ownership.
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
  allocation evidence. Current post-owner auxiliary-constructor measurements
  show only 31/122 eight-entry hits on fx/Quad (0.158/0.362 ms of observed work);
  do not implement shared-base caching from the stale pre-owner repeat counts.
  General producer/pre-overlay numerical identity is now disproven. The paired
  SC50A k1 supplied subset matches, as do 880 small/DEGEN2 supplied values, but a
  deterministic 51x51 dense handoff differs by one ULP before the mandatory
  overwrite. Do not treat sparse matches as a global proof or bypass the overlay.
  The strictly-lower-L contiguous-suffix dependency cut is now numerically
  verified with an independent complete constructor, dense mismatch, last-write
  duplicates/extremes/fallbacks and all SSE modes. Matching fixed 39 PGO reduces
  PILOT4 instructions 1.068%, but screen32217 is inconclusive; stop without
  seed32289 or full candidate gates. It remains unapplied. Do not repeat either
  this screen or the earlier arithmetic-only guarded-replay screen.
  A distinct fresh-constructor whole-pass fusion removes symbolic closure and
  sorting via a certified ascending numeric frontier. Complete constructor oracle
  and fast parity pass, but the fixed current-host PEROLD screen42117 is
  nonpositive (+3.925% median); no42189 or extension. It remains unapplied.
  Do not infer another replay optimization from graph-pass savings alone.
  The extended BTRAN conversion/operation-boundary audit is complete. Prepacking
  alone does not make the small component faster than software; rebuilding its
  cache per call is worse. Combined multiply/store/subtract blocks retain both
  roundings and improve the component, but matching-PGO integration raises
  fx/Quad cycles5.75/4.84% despite fewer instructions; rejected before33217/33289.
  Register-product range certificates pass1M raw triples but add only0.44%
  component improvement, insufficient evidence for another integration.
  Native division-routing measurement is complete:44.0/50.8% of dense-U
  positions on successful fx/Quad native calls have zero numerators. Dropping them
  before division preserves complete original outputs/WORK/state and matching-PGO
  fast940 parity, and lowers instructions, but fixed screen34217 establishes a
  Quad slowdown beyond AA uncertainty. Reject; no34289 or extension. Do not
  repeat the tested native dense U/eta/L variants from zero counts or isolated
  component improvements. Current post-owner auxiliary reuse does not justify
  sharing from historical repeat counts; use the current census in STATE.md.
  The ordinary-primal ranked-row consumer audit and maintained eligibility
  component are complete. Ordinary pivot streams need only nonbasics, all slacks,
  q and leave; the cost-change consumer must retain the full product. PILOT4 can
  omit59.55% of entries, but scanning basis positions/sets still loses in the
  component. Producer-notified maintenance improves PILOT4's component18.45%
  while regressing SC50A/DEGEN2; no unconditional solver integration follows.
  The fixed half-matrix-density/eight-pass warm-up gate is now tested across
  complete histories, including sparse fallback notifications. It remains slower
  on SC50A/DEGEN2 despite a PILOT4 component gain. Flat bitset membership likewise
  has no useful captured-stream gain. Do not integrate these tested mechanisms
  from arithmetic savings alone; actual refactor/repair binding remains absent.
  Original-scatter accumulator bounds elision passes the complete original
  oracle, memory-invariant proof, matching PGO and fast940. Its fixed screen35217
  is inconclusive with biased AA; do not repeat/extend it or run35289.
  Packed ordinary whole-row drop scanning is qualified: raw-bit bias/shift
  mapping makes positive-normal surrogate comparisons exact and exception-safe;
  ascending filtered indices retain the original body. Short (<128 columns),
  extended and nonprecomputed cases keep the range loop. Matching fixed 39 PGO,
  full987, all255 release tests and both DEGEN2-controlled screens36217/36289
  pass; all23 full result/exit pairs and the three extra regression controls
  qualify. No screen/control is extended or repeated.
  Accepted-build ordinary row/control attribution is refreshed on bounded
  SC50A/DEGEN2/PILOT4/PEROLD runs. DEGEN2 candidate-list update is6.88% instruction
  self; its annotation samples the linear bucket-threshold scan as well as
  fresh-membership paths. Before another change, freeze the current bucket loop
  as an independent oracle and test a direct exponent/mantissa classifier,
  preserving strict power-of-two boundaries, cap58, negatives/zeros/NaNs/infinity
  and floating-environment behavior. Keeping the original first comparison
  against2^-19 would reject exceptional/below-range cases before bit extraction;
  the proposed normal result is exponent_bits-1003-(fraction==0), capped58.
  This is an unimplemented hypothesis, not a verified rule or runtime gain.
  Measure bounded small-model counters before own-PGO qualification/timing.
  Current non-level sparse-L fallback census is complete: DEGEN2/PILOT4/PEROLD
  never enter Kahn replay; their L calls choose dense or level solvers. Do not
  implement numeric replay reuse on these targets from the unexecuted branch.
  Generated immutable dense L forward/transpose programs are tested separately:
  captured/live numeric, WORK and FP checks pass, but even a one-pass emitter
  loses when actual uses include construction/retirement. CPU2 PILOT4/PEROLD
  inclusive replay is about 63/29% slower. No per-LU-factor integration follows.
  A code program whose topology persists across all barrier Newton factors is
  a distinct unimplemented lifetime. CholSolver retains its Plan while Factor
  numerical arrays change: use a refreshed pointer table, not embedded values
  or allocation addresses. Forward needs both single/double peel schedules;
  backward shares geometry. Preserve dsub8/dsub2 reduction trees, four-product
  grouping and dynamic zero-diagonal rules. Current scope is measured below;
  triangular code generation remains unimplemented and lower priority.
  Pre-executor scoped DEGEN3 measurements put complete forward/backward near 7%,
  versus numeric factor 71% and child updates 35% (including DGEMM 19%). These
  nested release-test scopes are not accepted-release removable bounds.
  The complete solver-owned child/direct-GEMM-destination/rectangular-RHS
  executor is qualified: fixed screens 45631/45697 both pass, full 987 parity and
  release tests pass, all 23 benchmark outputs/exits agree. Independent output
  cells retain their original arithmetic; discarded products must still execute.
  Admitted A/B ranges guarantee finite normal-or-zero product intermediates.
  Gathered C retains its full logical finite-range check. Zero-C scatter accepts
  arbitrary parents because each final addition has at most one nonfinite
  operand. Keep the single-pass k gate, masked-exception check, distinct physical
  destinations, all logical A/B extents, signed-zero addition and original
  single/extreme/unmasked fallbacks. Passing strided parent slices to public
  GEMM remains invalid. Full original factor orchestration is frozen as oracle.
  Exact qualification and predecessor controls are in
  evidence/barrier-factor-executor.json; prior isolated 43117/44217 screens
  remain inconclusive and must not be repeated.
  Accepted-source normal/controller attribution is now refreshed on the six
  small models, 80BAU3B and GREENBEA at default 1000 limits. Broad numeric factor
  scopes are 44.82/50.77%, child 19.95/26.37%, normal solves 15.50/24.84%.
  These nested test scopes are not removable bounds. A complete trial executor
  targets the measured 8.76/6.58% selection scope, passes its full original
  selector oracle and fast 940, and lowers 80BAU instructions 5.869%. Its first
  fixed screen 48811 fails conservative AA separation; stop without 48883 or
  same-source extension/retry. Source/accepted CLI restored. Evidence:
  evidence/barrier-normal-scope.json and evidence/barrier-trial-executor.json.
  Preserve actual compiled FP semantics when replacing conditional kernels:
  lower ratios branch around arithmetic, while upper/tail loops speculatively
  divide before signaling comparisons. LLVM can erase neutralizing masks unless
  the instructions are protected; opaque reference calls and status tests matter.
  Packed-R assembly/updates/solve ownership and disjoint factor recycling remain
  unimplemented. Recycle only after the early return check and preserve par from
  logical lengths. The exact Plan census counts 13.0 MB R packing on DEGEN3,
  versus 83.2 MB original zero-product buffers now covered by admitted direct
  destinations. These are logical byte counts, not removable time. Prior arena
  and owner trials remain rejected; do not restart them without larger measured
  scope or assume that the independent historical gains add.
  The complete pre-U forward-column chain measures 2.72–3.64% on SC50A through
  DEGEN2, including required filtering and clock overhead. Tiny AFIRO measures
  7.25% with particularly material clock sensitivity. Do not implement the owner
  rewrite or escalate it to PEROLD.
  Preserve pending pre-U values through replacement, keep/basic and diagnostic
  consumers; retaining the overwritten U vector cannot substitute for them.
  Post-U output ordering is larger but natural traversal reuse is an existing
  timing-inconclusive trial and must not be repeated. Exact scopes and caveats:
  evidence/forward-owner-scope.json.
  Complete inner GEMM routing measurement is settled without implementation:
  construction-inclusive route-only proxies save about 0.086 ms on DEGEN2 and
  1.512 ms on DEGEN3, not solver speedups. Accepted-release DEGEN2 admission
  is about 1.67% instructions/1.64% cycles; S4x4 is 3.76%/3.92%, including
  required arithmetic/B loads. Healthy samples represent about 95.5% of direct
  batch counts. Geometry-only per-child packing reuse and repeated-validation
  counts do not establish a large removable cost; no certificate, descriptor
  executor or persistent packing rewrite follows. Keep exact accumulation
  order, K-boundary additions, full logical scans, unused outputs/flags and
  single/extreme fallbacks. No multi-K direct scatter. Raw evidence and excluded
  parser-only/throttled recordings: evidence/gemm-routing-scope.json.
  Complete predictor/HSD pairing is measured without implementation. GREENBEA
  has 72 pairs but zero overlapping correction rounds: predictors take none,
  HSD takes 97. Initial solves occupy 7.65% of the instrumented solve, or 8.70%
  with initial residuals; whole transactions include unmatched tails. Small
  models have no HSD calls. See evidence/hsd-pair-scope.json. No paired executor
  follows; do not repeat triangular-only repacking as a whole-pair experiment.
  Immutable-L prefix reuse has no entering-column lead on AFIRO through DEGEN2
  (cap 300): no same-column/support recurrences within LUOUT epochs. An optimistic
  support-only key gives one possible SC50A hit, two BLEND hits, and zero elsewhere.
  Other L calls and cross-epoch owner sharing remain unmeasured, without a high
  reuse invariant. No cache or large-model escalation. See evidence/l-prefix-reuse-scope.json.
  Direct outgoing-factor scopes reject a consuming installer: formation plus
  all Structure retirement is 3.109% of mode-switch, including required numeric
  work. Fresh metadata outside shared-L/transpose and Markowitz conversion is
  about 0.08%. Original lifetimes and all complete results are preserved.
  See evidence/factor-output-scope.json. Do not credit the sampled 9.58% union
  or 5.48% approximate residual as removable time.
  Healthy accepted-release primal-controlled profiles now cover DEGEN2 and
  80BAU3B. The latter's classification/score/row-update families are substantial,
  but old feasibility tracking and global buckets remain settled. Whole-score
  exact query reuse is separately rejected: 77.64% division repeats still yield
  +8.677% target instructions and +5.861% solve median in fixed pilot 53101.
  No formal 53131/53197 screen or pilot extension. Evidence:
  primal-controlled-scope.json and score-query-reuse.json under evidence/.
  Ordinary-primal paired BTRAN is also settled without implementation:
  controlled 80BAU offers only 4.281% shared coefficient loads under explicit
  monotone U-tail/eta/LT merges; even offline weighted row-LCS reaches 10.872%.
  Most LT pairs use level buckets, and 724 rho calls have no tau partner.
  Numeric divisions, controls, WORK and output remain. Small-model overlap
  does not justify an executor or another path-specific packet trial. See
  evidence/paired-btran-scope.json. Do not revive query-cache key/mask tuning,
  direct BTRAN rank/layout ownership or rho/Devex matrix-product fusion.
  A fresh substantial execution phase must be identified before another
  ordinary primal/dual implementation; aggregate hot-function shares alone
  do not qualify one. Complete independent-output LDL batching is settled:
  the full scalar oracle passes generated/actual fronts and fast940, but fixed
  DEGEN3 pilot54109 saves only 2.268% solve time against required5%, with noisy
  AA. No54131/54197, extension or packet-width/snapshot tuning. Costs include
  admission, matrix snapshot, validation and full matrix/FP-status replay.
  This is distinct from the earlier within-output chunk trial; neither gives
  a qualified runtime gain. No solve305 routing/codegen composition follows
  from settled triangular-plan scope. See evidence/ldl-output-batch-scope.json.
  The reconstructed reference and native timer/parameter harness are qualified;
  the full current-host fixed 23-row comparison is complete with 0/23 Rose wins.
  See evidence/current-reference-comparison.json and perf/compare_reference.py.
  Mode-switch dual is 8.21x slower, GREENBEA dual-release 3.97x, controlled 80BAU
  primal 4.88x. The GREENBEA solve-scope comparison is now qualified: 5.4666B
  versus 1.3461B instructions, healthy 2728/2688 samples and 99.8% direct coverage.
  At least 10.70% Rose self is pivot-row work already covered by settled trials.
  The allocator-owner lead is now closed: full 987-qualified cached FP source
  has +.586% GREENBEA instruction skew, healthy 2736 samples and 99.515% coverage;
  its 542 allocator leaves reach Rose caller context. Leading sites are outgoing
  mirror, Structure retirement and true_d; all overlap settled owners. Do not
  reopen buffer/cache/arena work from the 19.65% accepted-release allocator mass.
  See evidence/greenbea-allocation-owners.json for disjoint physical ancestry,
  explicit other contexts and immediate-callee/return-address proof.
  Next finish isolated qualification of the existing fixed 39 PGO/system-TCMalloc
  configuration with matching compiler tools, full parity and release tests,
  then measure it separately on this host. Preserve the ordinary accepted CLI,
  source/dependency hashes and baseline; no architecture or speedup claim from
  build/training alone. Continue searching for genuinely distinct whole-stage
  numerical work, not a different name for the settled owners.
  The scalar prefilter and prior rank/layout trials remain
  rejected; do not infer further gains from isolated predicate costs.
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
