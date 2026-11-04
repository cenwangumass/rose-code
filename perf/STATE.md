# Performance state

Current state of the work in [PLAN.md](PLAN.md). Rewrite in place; no log (history is `jj log`).

## Numbers on AMD EPYC 9V74

The starting solver is `934213eb`. The shipped fixed39 PGO executable runs with genuine Ubuntu gperftools TCMalloc 2.15. Gurobi uses the supplied patched 13.0.3 library, SHA-256 `26a505171314b3300b754b4eff5f5719366b47ccf6f3d59770df99e3f256ee30`. All 23 timed rows use a fresh process, a shared logical CPU, one warm-up per solver and three interleaved timed pairs. The table compares medians of Gurobi Runtime and rose diagnostics.seconds. The fixed GREENBEA PDHG cap is 20,000. This configuration supersedes the preliminary glibc-allocator survey.

| Run | Gurobi Runtime (s) | Rose solve (s) | Rose / Gurobi |
|---|---:|---:|---:|
| degen2.primal-c | 0.019153 | 0.043959 | 2.30 |
| degen2.dual-r | 0.013366 | 0.034534 | 2.58 |
| degen3.primal-c | 0.197820 | 0.531533 | 2.69 |
| degen3.primal-r | 0.239772 | 0.626245 | 2.61 |
| perold.primal-c | 0.198119 | 0.551410 | 2.78 |
| pilot4.primal-r | 0.020133 | 0.055149 | 2.74 |
| 80bau3b.primal-c | 0.132110 | 0.623391 | 4.72 |
| primal_5959.primal-r | 0.056607 | 0.119254 | 2.11 |
| degen3.dual-c | 0.096471 | 0.279956 | 2.90 |
| degen3.dual-r | 0.108019 | 0.280635 | 2.60 |
| perold.dual-r | 0.043718 | 0.130428 | 2.98 |
| 80bau3b.dual-r | 0.092070 | 0.355276 | 3.86 |
| greenbea.dual-c | 0.256883 | 0.771987 | 3.01 |
| greenbea.dual-r | 0.227844 | 0.772700 | 3.39 |
| mode_switch.dual-r | 0.117143 | 1.652152 | 14.10 |
| primal_5959.dual-r | 0.084074 | 0.205056 | 2.44 |
| quad_2717.dual-r | 0.105847 | 0.259393 | 2.45 |
| fx_quad_obj_5371.dual-r | 0.042942 | 0.105406 | 2.45 |
| degen3.barrier-r | 0.093120 | 0.203347 | 2.18 |
| 80bau3b.barrier-r | 0.148597 | 0.240573 | 1.62 |
| greenbea.barrier-r | 0.543017 | 0.996415 | 1.83 |
| 80bau3b.pdhg | 2.956445 | 5.214276 | 1.76 |
| greenbea.pdhg | 1.284483 | 2.117332 | 1.65 |

Starting result: rose wins 0/23. On mode-switch both solvers reach optimality after 11,393 iterations. The pinned Rust 1.98.1 release rebuild using TCMalloc and matching fixed39 training passes frozen full parity 987/987. All model parameters, solver hashes and samples from the starting sweep are in the release evidence.

### Mode-switch profile

Both solvers reach optimality in 11,393 iterations on the 15,000-row block witness. Three sampled rose runs in the preliminary glibc configuration yielded 1,531 whole-process CPU samples. All 158 samples in the maximum reduction returned to sparse forward-L level selection. Two full-dimension FTRAN marker scans accounted for another 130 samples and the sparse BTRAN range loop, including arithmetic on reached positions, for 74. Together these regions account for 23.6% of samples. The libc share was 619/1,531 for rose and 22/507 for the native reference runner across 15 runs. The profile includes each program's input and output work.

A separate run intercepted direct executable requests of at least 15,000 bytes. Rose requested 21.377 GB through calloc, 4.923 GB through memcpy and 0.049 GB through memset. BTRAN directly requested 6.856 GB through calloc. Repeated requests count separately. The probe measures requested sizes and cannot determine physical traffic or a speedup. All sampled and instrumented rose runs reproduced the complete uninstrumented result except elapsed seconds.

## Historical Ryzen 9 5900X measurements

The original host used TCMalloc and a different reference convention. The table and profiles in this historical section retain their original scope and provenance.

Baseline solver: `8b645007` (same solver as `vlsokskq`). Current solver: `sozuzvuv` (packed ordinary whole-row drop-predicate traversal, compact pure normal-Devex increment reuse, bounded-normal F80-to-double bit assembly, completed ordering hashes stored before allocation, reused adjacency buffers during quotient ordering, unique AA-neighbour construction, compact pure extended ranked-product-result reuse, compact native extended reduced-cost-result reuse, compact extended RHS-result/WORK-estimate reuse, four-row/four-column SIMD barrier DGEMM tiles with scalar two-column tails, allocation-free row-tile selection, whole-loop primal row dispatch, sparse auxiliary etas, optional Linux CLI `system-tcmalloc`, PGO with 64-byte loop alignment and online profile merging).
Parity: full tier 987/987; release tests pass. Three timed runs after one warm-up, with direct paired follow-ups noted below, hyperfine 2.0.0; instruction counts are millions from perf stat.

`baseline wall` is the fixed initial process-wall median. Other wall columns
are the current CPU4 three-run paired sweep. Solve uses that sweep except three
large positive medians within broad spreads, which use fixed three-pair follow-up
medians after additional reversed-role controls. Host speed bands affect snapshot
ratios against cached reference Runtime; only controlled within-session changes
are credited as gains. The first14 pairs used /tmp; after a baseline result-write
quota failure, mode-switch onward used the same workspace temporary directory for
both binaries. Solve diagnostics exclude output serialization; wall times across
those segments use different temporary storage.

Reference provenance: original-library cached records `../reverse-engineer/runs/{degen2,degen3,perold,pilot4}-{profile}.json` for the ten previously valid cells. New three-run medians use `/tmp/rose-perf-reference/{primal_5959.primal-r,primal_5959.dual-r,quad_2717.dual-r,fx_quad_obj_5371.dual-r}.{1,2,3}.json`; summary `/tmp/rose-perf-reference/summary.json`. All original-library records identify Ryzen 9 5900X and library SHA256 `8e87fca221fdffc922331c8619426031b3e72dfd419446422f88af0cb0e98f9c`.
`-` refs are unreachable: original Gurobi rejects 80BAU3B, GREENBEA and mode-switch with license error 10010 (scratch `ID.error.json`, including direct 80BAU3B primal/dual attempts). The old 80BAU3B primal/dual refs used the patched library and have been removed; the old prose mode-switch time is unusable.

Reference instruction measurement (`1a430e43` Rose CLI, original Gurobi 13.0.3, CPU 4, three alternating one/four-model pairs): median reference slope `(four - one) / 3` is 3.676M instructions/model on ADLITTLE and 323.665M on PEROLD; Rose whole-CLI medians are 12.108M and 1190.296M. The reference slope excludes fixed Python/uv import startup but includes fresh model reading, parameters, solve, logging, attributes and disposal, with one shared Env/process. Rose includes CLI startup/parsing/output, so these are not solve-only or identically scoped counts. All runs are optimal in 73/1273 iterations, caps 300/2000 do not bind; params `Method=1, Presolve=0, NetworkAlg=0, Threads=1, Seed=0`, other solver params default. Actual mapped original-library path/hash verified after import and model read; `LD_LIBRARY_PATH`/`LD_PRELOAD` unset. Sequential, unwarmed timing medians are reference 0.000629/0.055369 s and Rose solve 0.002760/0.204713 s; these context measurements do not replace the official Runtime targets above. Raw data, exact commands/environment, runnable measurement and mapping proof: `/tmp/rose-ref-measure-20261103/{summary.json,commands.json,measure.py,mapped-library.json,provenance.py}`. Use a solve-family instruction-profile comparison before attributing the remaining gap to cache locality alone.

Instruction-family self profiles (accepted `1a430e43`, two period-1M PEROLD samples): Rose versus original-reference read/solve-slope mass is roughly 283–287M versus 42–45M in factorization/fresh structure and update wrappers, 332–342M versus 114–126M in FTRAN/L, 213–223M versus 88–101M in pricing/ratio/control, and 54–60M versus 41–45M in BTRAN. These are self classifications, not inclusive caller ancestry, with tens-of-millions uncertainty per major family; Rose parsing is another 7–8M and native unmapped remainder about 10–11M. The large libc remainder resolves to allocator internals plus malloc/free/calloc/realloc (213–215M), copying (15–17M), memset (12–13M), and only 1–3M getenv, using mapped ELF offsets, disassembly and resolved IFUNC entries. Allocator instruction mass remains unassigned to callers. Period-50k samples were discarded for severe underrepresentation; period-1M represented totals match the counted ~1190M/~324M scopes. Artifacts and runnable analyzers: `/tmp/rose-family-profile-{1m,repeat}/`; libc mapping `/tmp/rose-family-profile-1m/libc-primitives.json`. ADLITTLE was profiled first but is too small for useful families at this period.

Native release instruction-family self profiles (`lpnyrpqs`, stale after sparse auxiliary eta storage; exact preceding accepted PGO/TCMalloc CLI SHA256 `74de62b7f700148104ac2ca6a2b14adf07ebef3bd3415f230165e71e4d1616df`): two CPU-4 batches of 16 PEROLD solves, cap 1300, all optimal in 1273 iterations/27 refactors. Fixed period 4M gives 4217/4218 total samples, zero throttle/loss records, and represented batch totals 99.75/99.78% of direct counts; filter `comm=rose` for solver self costs. Rose versus cached original-reference read/solve-slope self mass, millions/model: FTRAN/L 356.5–357.8 versus 113.7–126.3, factor construction plus update wrappers 242.5–246.3 versus 42.3–45.3, pricing/ratio/control 201.0–206.3 versus 88.3–101.0, BTRAN 55.5–59.0 versus 40.7–45.0. Native allocator is 93.3–95.0M; other libc 38.8–42.5M, whose callers remain unassigned. Parsing 7.8–9.8M, setup 6.0–6.8M, Rose helpers/unmapped 27.3–32.3M and unknown/external remain separate. Reference classifications/scopes match the earlier cached profiles; these are coarse sampled self estimates with tens-of-millions uncertainty in major families, not inclusive costs or confidence intervals. Release function self estimates put first-vector FTRAN at 109–111M, needed-vector FTRAN 92–95M, L stage 92–94M, numeric auxiliary LU 49–51M and nucleus 51–56M. Prioritize the FTRAN/L and factor/update excesses rather than interpreting allocation counts as time. Direct release count median is 1050.501M; frame-pointer profiling CLI 1117.384M (+6.37%), so its separate instruction recordings are excluded. Wrapper-heavy exploratory recordings are also excluded; final wrapper validates result JSON outside the measured process. ADLITTLE was measured first but yields too few samples for families. Commands, source/classifiers, exact build/hash, direct counts, audits, summaries and disassembly: `target/perf-native-families/{collect.py,run.py,analyze.py,manifest.json,analysis.json,summary.json,symbols.json,lstage-instructions.asm,perold-*.data,frame-pointer/direct-counts.json}`. Reference source `/tmp/rose-family-profile-{1m,repeat}/analysis.json`; official Runtime targets unchanged.

Allocation caller sampling of the `1a430e43` frame-pointer PGO binary (ADLITTLE first, then PEROLD caps 300/1300; all optimal at 73/1273 iterations): PEROLD makes 563,368 malloc, 32,023 calloc and 221,924 realloc calls, with 593,333 frees, over the whole CLI. Temporary interposition samples stacks every 61 or 64 allocation calls separately per primitive, excludes recursive profiler allocations and drops no stacks; sampled caller rankings agree across periods. Top individual sites include dual true_d column clone/sort, phase row adjacency, nucleus search chain clones, pricing reset row adjacency and sparse-L level queue growth. The first and third already have rejected trials below. Sparse-L level queue pushes account for roughly 18,000 calls; this estimates allocation frequency, not the instruction cost of a caller. Instrumented times are diagnostic only. Reproducible sampler, commands, exact binary hash, raw stacks and inline-aware analysis: `/tmp/rose-alloc-profile/{manifest.json,run.py,sample.c,sample61.c,analyze.py,perold.alloc,perold61.alloc,*.summary.json}`.

Current native release instruction profile (`vypkuxry`, CLI SHA256
`998efcda56a3d020475cd730cd64a80766bbe5cf8b73fc49122d862ad095d630`):
ADLITTLE/DEGEN2 bounded batches precede sixteen PEROLD dual solves, cap1300,
all optimal at1273 iterations. Period4M records4302 samples, no loss/throttle,
and represents99.87% of the directly counted17.230B batch instructions.
Whole-batch self shares at Rose symbols are L-stage10.18%, first FTRAN8.83%,
needed FTRAN7.23%, iteration12.95%, BTRAN6.04%, nucleus5.25%, update5.97%,
auxiliary factor4.09% and fresh3.44%. Scope includes launcher, parsing,
startup and validation; these are self observations, not inclusive solve-only
bounds. L-stage annotation still shows stack-indexed eight-lane eta partial
sums, dense numeric scatter and a spilled dense WORK accumulator. Separate immutable
charge caching and complete-loop extraction below reduce little target work and
have no useful counter gain; extraction does retain WORK in a register. The long-only
SIMD eta trial below removes the indexed sums but has no qualifying timing
gain. Exact accepted build/hash, bounded runner, direct counts, recording,
event audit and both kernel annotations: `target/perf-ordinary-devex/`.

Native allocation caller sampling (`lpnyrpqs`, saved PGO/loop-alignment/frame-pointer flags and `system-tcmalloc` feature; ADLITTLE cap 300 before PEROLD cap 1300, both normal optimal exits): PEROLD makes 546,494 malloc, 32,077 calloc and 214,625 realloc calls, with 576,513 frees. Independent periods 61/64 give identical totals, 13,004/12,395 sampled stacks and no drops. Refactorization ancestry accounts for 58.34/58.31% of sampled allocations, iteration outside refactorization 19.21/19.14%, other/setup 22.45/22.54%; these are allocation frequencies, not time shares. Both periods rank true_d column copies, nucleus search clones, pricing-reset adjacency, sparse-L queues, auxiliary numeric LU and fresh-factor lists among the leading sites. Several already have rejected trials below; prioritize produced factor storage/basis-column copies before another broad iteration workspace. The interposer calls exported tc_malloc/tc_calloc/tc_realloc/tc_free directly, with dynamic binding and mapping evidence, preserving the accepted allocator. Child-only instrumentation passes fast parity 940/940; never preload the non-native parity harness itself (mixed allocator alignment paths can abort). Commands, exact binary hash/build, sampled stacks, summaries and bindings: `target/perf-native-alloc/{manifest.json,sample61.c,sample64.c,bindings.txt,perold{61,64}.{alloc,summary.json,stacks.json},fast-parity.txt}`. Instrumented timings are diagnostic only.

PILOT4 primal instruction comparison (`5b4e350f`, cap 800, all optimal in 690 iterations and 16 Rose refactorizations): direct Rose whole-CLI count 528.698M versus original-reference fresh-read/solve slope 152.956M. Corrected period-8M profiles use two independent reference `(33 - 1) / 32` batches and Rose `(9 - 1) / 8` CLI batches; the latter include wrapper dispatch/result handling (counted slope 534.892M). All ten recordings have zero throttle/loss records; represented batch totals differ from counted totals by 0.003–1.56%. Period-1M data was throttled and is excluded. Observed self masses, Rose versus reference, are pricing/control/row products 158–166M versus 47–49M, factor construction/fresh/update wrappers 126–128M versus 25–31M, libc heap management 86–91M versus approximately zero, FTRAN/L 49–50M versus 24–26M, and BTRAN 35–40M versus 32–33M. Prioritize the first three excesses; allocator caller ancestry remains unknown. These ranges are repeat spread, not confidence intervals; coarse sampling and unmatched scopes leave tens-of-millions uncertainty in large buckets. Native unmapped 11–14M and libc unmapped 5–6M remain explicit. Library provenance/mapping verified; official Runtime targets unchanged. Commands, audits and analyzers: `/tmp/rose-pilot4-primal-families-8m/{summary.json,audit.json,analysis.json,commands.json,run.py,analyze.py}`.

PILOT4 primal allocation sample (`lotvwumw`, preceding frame-pointer PGO build, ADLITTLE first; caps 300/800, both optimal, PILOT4 690 iterations): 186,115 malloc, 11,291 calloc and 100,651 realloc calls, 195,939 frees across the whole CLI. Period-61 primitive-specific stack sampling has zero drops; 1999/4889 samples (about 41%) originate in ftran_order.rs, 1232 in lu_order.rs and 893 in simplex.rs. Prominent individual sites include temporary fresh-factor adjacency, LU pair-row growth, nucleus construction/search, and copying index/value arrays for basis columns. Counts measure allocation frequency, not caller instruction cost. Raw `target/perf-pilot4-trsm-{adlittle,pilot4}.alloc`, summaries/stacks beside them, manifest `target/perf-pilot4-trsm-alloc-manifest.json`, analyzer `target/perf-pilot4-trsm-alloc-analyze.py`; the cached sampler is `/tmp/rose-alloc-profile/sample61.so`.

Compact shared providers retain up to 32 pure extended ranked-product results
for their Solver's immutable problem, independently of main-factor updates. Keys
are the participating rows in sorted `(rank, row)` order and each value's exact
extended representation. Rank magnitudes and unselected values do not affect
the kernel and are omitted. Signed zeros and full low bits remain distinct.
Only compact shared owners, at most 50 rows/500 variables, native nearest/64-bit
masked x87 mode and valid dimensions are eligible; selected NaNs and all dumps
use the original path. The numerical kernels are source-identical, and both
consumers retain every WORK call in place. Results and scratch stay independent
of mutable main factors. Release tests cover rank ties/order changes, low bits,
signed zeros/extremes/infinity, ignored values, clones, mode fallback and eviction.

The test-only working-set measurement records 25,046 eligible calls with 22 exact
inputs on Quad_2717, and 6,073 calls with 23 inputs on fx_quad_obj_5371. A four-entry
cache misses the repeating cycle and regresses; that version is retained as
negative evidence. The final 32-entry bound fits the observed working sets. The
trace is removed before source/test freeze and all official timing screens.

Earlier ranked-product reuse passed two independent screens (-7.52/-8.51%),
full parity and release tests, plus all 23 complete benchmark pairs and regression
controls. Its ownership/key rules remain unchanged. Earlier qualification and
negative four-entry-cache evidence: `target/perf-rank-memo/`, qualified `v2/`;
RHS/WORK and reduced-cost evidence: `target/perf-rhs-memo/final/` and
`target/perf-cost-memo/`. The ranked-product-qualified CLI is
`13eba41bca5463721d08e6a634d066449b3d2e5de5d224200ecad88472ebfb9b`.

AA-pattern construction now builds non-dense row-to-column incidence and marks
unique off-diagonal neighbours before sorting each row. Only unique neighbours
are stored; dense-column extension/order/multiplicity and final CSR order are
unchanged. The downstream ordering, degree arithmetic, numeric kernels and WORK
source remain identical. The independent complete original constructor checks
512 generated matrices with duplicates/zero/NaN coefficients, empty/all/repeated
and reordered dense selections, plus AFIRO/SC50A/BLEND/ADLITTLE/BORE3D/DEGEN2.
Both AA-pattern tests pass with the matching fixed39 PGO build.

Complete AA-pattern qualification and original-source proof: `target/perf-aat-unique/`.

Quotient ordering now reuses each row's owned adjacency buffer. Both segments
are compacted stably before inserting the pivot at their original boundary,
which prevents overwriting unread tail entries. Prefix and tail hash additions
keep their original order and conditions; list contents, degree updates and
WORK remain unchanged. The complete original ordering is independently frozen.
Tests compare permutations, NZ and exact Ops/WORK bits on 512 generated graphs,
multiple levels including large labels that exercise hash rounding, and six
small models with dense nodes and generated levels.

Complete row-buffer qualification: `target/perf-order-reuse/`. The paired predecessor CLI is `ba7d7a7031f73e68fd126e9f2ef8fcc265d7beefcb705dedda671216ad6e5e1d`.

Each completed row hash is now stored before vector truncation/insertion,
ending the accumulator's lifetime before a possible allocation. No hash read
intervenes at the old store position. Every hash addition, degree operation and
WORK operation retains its original order; the existing independent original
ordering tests remain unchanged. Matching release assembly retains the hash
accumulator in a register in both loops, removing their per-term stack
load/store round trips. Ordering code is 17,697 bytes versus 17,761 before.

Complete hash-lifetime qualification remains in `target/perf-order-hash-store/`.
Its qualified CLI `90878f676735bad3d8fab1f9292db88e3813231b2acd558ec927e47f2229aa93`
is the current paired predecessor.

F80-to-double conversion now directly assembles the already rounded double
bits for normalized nonzero values with exponent -970 through 1022. This has
the same 53-bit significand rounding, ties-to-even condition and carry as the
original conversion. Within this domain the original integer conversion,
power-of-two scale, product and sign operation are all exact, normal and
finite under every MXCSR rounding/DAZ/FTZ combination. Zeros, nonfinites and
all boundary/out-of-range exponents retain the source-identical original
routine. No solver arithmetic, operation ordering, WORK or dump booking changes.

The conversion differential helper makes its input opaque with black_box before
both calls. This prevents constant folding of boundary expressions across the
MXCSR changes in its mode test. A later division prototype exposed a
negative-overflow infinity versus maximum-finite comparison under directed
rounding; opaque input preserves every existing case and makes both paths
execute at runtime. Both the prototype and accepted production source pass
these unchanged cases after the test fix. Production code, CLI and PGO profile
are unchanged. Failure, source snapshots and accepted-source checks:
`target/perf-f80-div-normal/`, `target/perf-f80-mode-testfix/`.

The complete original conversion is frozen independently before the rewrite.
Two matching-PGO release tests compare all 65,536 stored sign/exponent words
with ten significands, one million random representations, captured SC50A
factor values plus low-bit variations, and all 16 MXCSR rounding/flush modes.
Exact bits agree, including ties, carries, subnormals, overflow and NaNs.
Native assembly uses integer rounding/bit assembly on the admitted path.

Complete bounded-normal conversion qualification remains in
`target/perf-f80-double-bits/`. Its qualified CLI
`569e13b50e004005676d6eb1a2334a3cf8ef6d00a118fbc81d9e875c95711c6c` is this trial's paired predecessor;
the opaque-input test fix changes no production code.

The compact shared provider now retains up to eight FIFO normal-Devex input
keys. Each key preserves the full extended representations of `a_rq`, `wq_q`,
every `pq` value and every cross `tau` value, including signed zeros and low
bits. Its Solver owns immutable problem columns. Per-column values are filled
lazily with exactly `ratio * (wq_q * ratio - 2 * cross.dot(column))`, where
`ratio = pq[j] / a_rq`. Every call still reads and adds the current old weight,
and retains all floors, costs, twins, state flags, reset/output order and WORK.
Misses compute ratio/dot at their original point before reading the old weight;
the remaining expression retains its original operation order.

Only native x86_64 Linux, packed-cross, normal mode with at most 50 rows and
100 total variables, valid dimensions, finite key inputs, nonzero divisor and
standard masked nearest/64-bit x87 control is admitted. Dumps, mode 2 and all
unsupported cases run the unchanged original kernel. No F80 primitive changes.
The complete original production and independent frozen row bodies remain
source-identical. A reversal proof checks the cached specialization against
the original. Every eligible row in release unit-test builds also executes the
independent frozen kernel and compares full mutable weight/cost/twin bits, reset, changed
order and complete ordered output. New tests exercise actual warm hits on
bounded Quad/primal5959 solves, plus fx rows, low key bits, eviction and control
fallbacks.

Initial complete bounded ADLITTLE/fx/Quad pairs agree. Quad instructions fall
5.90% and median cycles 8.31%. Two independent Quad screens (32 AB + 32 AA,
CPU 4, cap 17000, seeds 14221/14293) pass every criterion: median solve
-11.104/-9.825%, intervals [-13.671%, -9.874%] / [-11.300%, -8.986%],
AA [-2.811%, +0.339%] / [-2.080%, +1.219%]. All 32 pairs improve in both;
sign p 2.33e-10 each, both execution orders improve and AA controls are healthy.

Full frozen parity passes 987/987. All 37 release groups pass (91 library, 254 total tests).
Formatting passes and all 23 complete no-dump benchmark result/exit pairs agree.
Source/tests and the screened CLI remain frozen. The 2 full-sweep timing
flags each receive fixed three-pair follow-ups and one predeclared reversed-role
control (A=candidate, B=predecessor), rejecting established regressions even
below 0.5%:

degen3.primal-c: follow-up median solve +0.11%, cycles -0.00%.
greenbea.pdhg: follow-up median solve +0.62%, cycles +0.57%.

degen3.primal-c (seed 14419): predecessor median -0.037%, interval [-0.757%, +0.832%], AA [-1.124%, +0.528%], sign p 0.57.
greenbea.pdhg (seed 15039): predecessor median +0.261%, interval [-0.910%, +0.751%], AA [-0.786%, +0.590%], sign p 0.811.

No control establishes a regression; every AA control passes its bias check.
No samples are extended and no control is repeated.

Devex-increment source/build, source/ownership proof, differential tests, both
screens and every gate/result/control: `target/perf-devex-delta/`. Fixed 39-case PGO uses
online merge, 64-byte loop alignment and system TCMalloc.
Profile SHA256 `79f8f3009d54e96bd4d3db1caa16684a4bf4477e58f31120852c6f02f304a84a`;
CLI SHA256 `998efcda56a3d020475cd730cd64a80766bbe5cf8b73fc49122d862ad095d630`.
Frozen baselines and cached original-reference Runtime targets are unchanged.


Ordinary precomputed whole-row updates of at least128 columns now scan eight
product bit patterns at once on x86_64 and visit only values for which the old
`abs(arj) <= ROW_DROP_TOL` drop test fails. NaNs remain in the traversal; their
basic-slack/nonbasic routing stays in the unchanged body. Integer bias/shift
mapping turns every magnitude, including exceptional inputs, into a positive
normal surrogate before packed comparison. It preserves the exact predicate
without adding invalid/denormal exceptions. Indices stay ascending. Short,
extended and nonprecomputed entries retain the range traversal; every arithmetic
expression, mutable-field effect, returned list and caller WORK estimate remains.

The complete current original body is independently frozen; helper reversal
restores it exactly, and full-source reversal is exact modulo formatting. Existing
whole-state captured-AFIRO/random/extreme/partial-extended families now include
128-column boundaries. Predicate tests compare2095552 values across16 MXCSR
modes, tolerance-adjacent bits, signs, subnormals/infinities/two NaN payloads and
unmasked invalid/denormal traps. All weights/costs/twins, norm/reset and list/value
order agree. The diagnostic full-row component uses captured used inputs for
SC50A8, DEGEN2181 and PILOT4367 calls, but its inlining/reset boundaries differ
from the solver and are not qualification evidence.

Matching fixed39 PGO/TCMalloc target counters lower DEGEN2 controlled instructions
3.59% and cycles3.17%. Two independent32AB+32AA screens (CPU4, cap800,
seeds36217/36289) qualify: median -3.242/-4.510%, CI[-5.244,-1.919]% /
[-5.850,-3.329]%,27/32 and30/32 improving, sign p5.65e-5/1.23e-7. AA is healthy
(p .215/.860), intervals[-1.840,+1.714]% /[-1.570,+1.607]%; both execution
orders improve in each screen.

Full frozen parity passes987/987, all37 release groups pass (92 library,
255 total,22 default ignored), formatting passes, and all23 complete no-dump
benchmark result/exit pairs agree. The sweep has no beyond-spread slowdown flags.
Three large positive medians within broad spreads get conservative extra checks:
fixed three-pair follow-ups and one already-predeclared32AB+32AA reversed-role
control each (A=candidate, B=predecessor), rejecting established regressions even
below0.5%:

degen2.primal-c (seed 36531): predecessor median +4.452%, interval [+2.841%, +6.420%], AA [-1.913%, +0.971%], AA p 0.377.
degen3.primal-r (seed 36624): predecessor median +2.804%, interval [+2.210%, +2.954%], AA [-0.253%, +0.777%], AA p 0.86.
greenbea.dual-c (seed 36903): predecessor median -0.118%, interval [-0.649%, +0.206%], AA [-0.547%, +1.736%], AA p 0.597.

No control establishes a regression; every AA control is healthy. No screen or
control is repeated or extended. The accepted baseline mode-switch solve first
hit /tmp quota while writing its result, before any valid benchmark measurement;
that failure is preserved. First14 completed pairs are retained, remaining pairs
resume using a workspace TMPDIR, and all complete output checks pass. No numeric
baseline, benchmark definition or original-reference Runtime target changes.

Qualified source/build, exact body/mapping/ownership proof, predicate/full-state
checks, both target screens, complete benchmark/follow-ups/controls and quota
recovery: `target/perf-row-mask/`, subdirectory `solver/`.
Profile SHA256 `b549084be3fb4e0937cab94702f2951818b2c8fb393dbc4f827377a0575eabc0`;
CLI SHA256 `db2bff5ae64aab9068f47ac2ee5b0395a6103916262df1760c18ffcfc36dcb89`.

| run | ref solve | baseline wall | paired base wall | current wall | solve | current/ref | instr |
|---|---|---|---|---|---|---|---|
| degen2.primal-c | 0.012711 | 0.182 | 0.064 | 0.070 | 0.049555 | 3.90 | 255.0 |
| degen2.dual-r | 0.009653 | 0.125 | 0.060 | 0.063 | 0.052485 | 5.44 | 203.2 |
| degen3.primal-c | 0.163137 | 3.237 | 1.009 | 0.906 | 0.884517 | 5.42 | 3664.0 |
| degen3.primal-r | 0.193685 | 3.872 | 0.816 | 1.104 | 0.772399 | 3.99 | 4177.5 |
| perold.primal-c | 0.160994 | 2.430 | 0.819 | 0.488 | 0.480174 | 2.98 | 4273.3 |
| pilot4.primal-r | 0.015787 | 0.253 | 0.055 | 0.054 | 0.047892 | 3.03 | 469.2 |
| 80bau3b.primal-c | - | 1.618 | 0.612 | 0.523 | 0.508453 | - | 5023.4 |
| primal_5959.primal-r | 0.067884 | 0.563 | 0.122 | 0.124 | 0.120107 | 1.77 | 1247.7 |
| degen3.dual-c | 0.082843 | 1.002 | 0.246 | 0.246 | 0.235856 | 2.85 | 1828.4 |
| degen3.dual-r | 0.082947 | 1.008 | 0.248 | 0.249 | 0.238040 | 2.87 | 1836.5 |
| perold.dual-r | 0.038532 | 0.355 | 0.128 | 0.128 | 0.120915 | 3.14 | 1046.1 |
| 80bau3b.dual-r | - | 1.103 | 0.312 | 0.312 | 0.297011 | - | 2475.0 |
| greenbea.dual-c | - | 2.870 | 0.661 | 0.840 | 0.997072 | - | 4916.4 |
| greenbea.dual-r | - | 3.598 | 1.031 | 1.032 | 1.011234 | - | 4856.9 |
| mode_switch.dual-r | - | 25.815 | 1.290 | 1.279 | 1.254626 | - | 12067.9 |
| primal_5959.dual-r | 0.098192 | 1.300 | 0.174 | 0.176 | 0.170853 | 1.74 | 1769.7 |
| quad_2717.dual-r | 0.092814 | 1.130 | 0.183 | 0.186 | 0.181244 | 1.95 | 1801.6 |
| fx_quad_obj_5371.dual-r | 0.038077 | 0.334 | 0.076 | 0.076 | 0.072023 | 1.89 | 750.4 |
| degen3.barrier-r | 0.075214 | 0.608 | 0.206 | 0.206 | 0.171778 | 2.28 | 2447.5 |
| 80bau3b.barrier-r | - | 0.379 | 0.225 | 0.226 | 0.204538 | - | 2626.7 |
| greenbea.barrier-r | - | 1.569 | 1.033 | 1.129 | 1.088304 | - | 9885.1 |
| 80bau3b.pdhg | - | 4.552 | 5.153 | 5.155 | 5.133938 | - | 74939.5 |
| greenbea.pdhg | - | 2.801 | 1.593 | 1.585 | 1.572746 | - | 31947.1 |

## Hotspots

Each profile below identifies its saved solver/build. Entries naming an earlier
source predate current changes and are historical attribution, not current cost
bounds. Reproduce their saved flags/features when consulting annotations; refresh
attribution on a changed path before selecting another optimization.

Per profiled bench run: the commit it was taken at and the top self-time functions (percent), at most
five lines each. Re-profile a run when a landed change touched its hot path.

- degen2.dual-r, `uoskoqqp` (stale after dense L-output compaction; PGO plus loop-alignment and frame-pointer flags; full profiling parity 987/987; eighty solves capped at 2000 iterations, all optimal in 634 iterations and 16 refactorizations; `target/perf-degen2-current.data`): FTRAN 8.39%, iteration 8.25%, L stage 5.44%, BTRAN 5.41%, update 4.64% self. Iteration 44.68%, refactorization 21.11%, FTRAN 14.34%, needed-vector FTRAN 12.65%, L stage 9.82% inclusive. Dense L output collection is another 2.32% self; malloc/free total 3.97%, realloc 0.85%. FTRAN annotation (`target/perf-degen2-current-ftran-annotate.txt`) samples full marker scans, numeric scatter and division. L-stage annotation (`target/perf-degen2-current-lstage-annotate.txt`) samples dense traversal, scatter and ordered work accumulation; sparse eta-accounting setup/gather is cold, so general prefix sharing remains untested. Numeric LU is 4.08% self/4.98% inclusive; its annotation (`target/perf-degen2-current-factor-annotate.txt`) samples scatter indexing, threshold-pivot scans and L divisions. Bitset ordering of touched rows is rejected below; pivot buckets already use indexed circular lists. Update annotation (`target/perf-degen2-current-update-annotate.txt`) samples row-capacity counts, structural filters and the empty-tail maximum scan. A/B membership queries already start from maintained row lists; a separate reverse A index remains untested. Dense paired L/U interleaving and individual U-estimate deferral are rejected below; remaining numerical traversals are leads.
- degen2.primal-c, `rtyoxnwx` (stale after phase-one class reuse; frame-pointer build passes full 987/987 parity; forty solves capped at 2000 iterations, `target/perf-degen2c-next.data`): BTRAN 9.80%, primal control 7.95%, FTRAN 7.22%, update 4.98%, LU nucleus 4.26% self. Refactorization is 24.70%, FTRAN 13.22%, BTRAN 10.22% inclusive. All forty optimal in 531 iterations, 16 refactorizations. Signature construction is 1.26% and sig_class 1.11% self; repeated phase-one classification within the stable update window is now shared. Bounds-check removal and direct growth-diagonal access were measured and rejected below. Reprofile before choosing remaining BTRAN or refactor work.
- perold.dual-r, `kpttrrms` (accepted native-TCMalloc/PGO release CLI; ADLITTLE first, then 24 dual-release solves capped at 1300, all optimal in 1273 iterations with matching complete outputs): period-4M instruction self profile has 6418 samples, zero loss/throttle events and represented instructions 99.73% of direct 25.742B whole-batch count. Iteration is 11.55%, FTRAN 9.26%, L stage 9.15%, factor update 7.70%, needed-vector FTRAN 7.43%, BTRAN 5.17%, nucleus 5.03%, auxiliary factor 4.19%, post-inner control 3.63%, fresh factor 3.55% and joint-U WORK accounting 2.45% self. Scope includes Python launcher, parsing, startup and validation; these are self shares, not inclusive or removable bounds. F2 snapshot allocation/ownership reuse is rejected below. Current build/hash/commands, validated results, direct counts, raw-event audit and report: `target/perf-f2-workspace/`. Earlier FTRAN/L/needed annotations remain in `target/perf-forward-current/`; use the current profile for cost claims.
- perold.primal-c, `wvqvkoyp` (stale after PGO collection-mode correction; PGO frame-pointer build passes full 987/987 parity; three solves capped at 3500 iterations, `target/perf-peroldc-pgo.data`): primal control 16.89%, BTRAN 9.75%, LU nucleus 4.81%, FTRAN 4.78%, score 4.48% self. Refactorization is 23.52%, BTRAN 9.99%, FTRAN 9.12%, row-product work 1.91% inclusive. All three optimal in 3399 iterations, 90 refactorizations. Score repeatedly checks immutable free-variable bounds; prior row-product and solve trials remain below.
- degen3.primal-c, `oytzqrpy` (stale after growth-factor reuse; frame-pointer build passes full 987/987 parity; three solves capped at 10000 iterations, `target/perf-degen3c-next.data`): primal control 8.91%, BTRAN 8.09%, FTRAN 6.51%, LU nucleus 5.82%, fresh factor 5.25% self. Refactorization is 31.03%, reference_order_for 19.77%, BTRAN 11.48%, FTRAN 10.88% inclusive; dense BTRAN output compaction is another 2.58% self. All three optimal in 2492 iterations, 49 refactorizations. The accepted growth-check factor now carries through to installation; reprofile before choosing remaining reference-order or BTRAN work.
- degen3.primal-r, `votqkmtm` (stale after PGO; ThinLTO frame-pointer build passes full 987/987 parity; three solves capped at 5300 iterations, `target/perf-degen3r-thin.data`): BTRAN 10.62%, FTRAN 10.16%, update 6.43%, primal control 5.46%, factor_structure 4.37% self. Refactorization is 24.75%, FTRAN 18.32%, BTRAN 11.12% inclusive. All three optimal in 5165 iterations, 66 refactorizations. Signature construction is 2.60% self and repeatedly recomputes immutable bound flags and widths. Prior output/solve micro-optimization trials remain below.
- quad_2717.dual-r, `nzspmlyk` (stale after native Devex cross-term dots; all saved PGO flags plus frame pointers; full profiling parity 987/987; eight solves capped at 17000, all optimal in 16010 iterations and 3022 refactorizations; `target/perf-quad2717-taugather.data`): primal control 17.78%, F80 add_signed 13.78%, F80 multiply 6.00%, rank_products_q 5.45%, extended BTRAN 4.93%, set_d_from_y 3.91% self. Refactorization 12.13%, extended BTRAN 12.02% inclusive; libc malloc/free total 6.17% self. Current annotations `target/perf-quad2717-taugather-{primal,rank}.txt`. The native ranked-product loop stores/reloads each 80-bit product before adding it to its accumulator; removing that traffic did not establish a repeatable timing gain (see below). Caller-site attribution assigns 5.60% of total cycles to F80 addition directly from primal control, about 3.20% within the Devex dot/reduction now ported; analysis `target/perf-quad2717-{stacks-offset.txt,add-sites.json}`, script `target/perf-quad-caller-sites.py`. Reprofile before assigning the remaining software arithmetic.
- degen3.barrier-r, `ztrylnzw` (accepted native-TCMalloc/PGO release CLI, BORE3D first; eight barrier-release solves capped at 100, all optimal with matching complete outputs): period-8M instruction self profile has 2639 samples, zero loss/throttle events and represented instructions 99.93% of the direct 21.128B whole-batch count. Numeric factorization is 31.19%, DGEMM 14.85%, DTRSM 12.01%, ordering 8.11%, solve 6.33%, partition 3.71%, index sorting 3.49% and assembly 1.40% self. These are batch fractions including launcher, parsing, startup and validation, not inclusive/solve-only shares or removable bounds. Prior eight-row, two-column, copy-Chain, row-scatter and LDL-chunk rejections still apply. Exact build/hash, runner, validated results, raw samples, direct count and event audit: `target/perf-barrier-0436/`.
- primal_5959.dual-r, `tplnupps` (stale after factor reuse; frame-pointer build passes full parity; three solves capped at 16000 iterations, `target/perf-dual5959-full.data`): F80 add_signed 7.50%, primal control 5.70%, F80 multiply 4.54%, F80 nucleus 3.88%, double nucleus 3.74% self. All three reach numerical_failure, 15688 iterations, 7842 refactorizations. Refactorization is 51.10% inclusive, quad mirrors 21.69% before factor reuse. LUIN capture `target/perf-dual5959-luin.log` has 7822 records and 40 input-head/threshold combinations; two-entry LRU estimate 99.3% hits. Wide-shift arithmetic rewrites are rejected below. The earlier 11000-iteration capture stopped early and is not the full-run profile.
- pilot4.primal-r, `ztrylnzw` (accepted compact-reduced-cost/native-TCMalloc/PGO release CLI, ADLITTLE first; 32 primal-release solves capped at 800, all optimal in 690 iterations with matching complete outputs): period-4M instruction self profile has 3934 samples, zero loss/throttle events and represented instructions 99.70% of the direct 15.783B whole-batch count. Primal control is 13.29%, ordinary row update 11.06%, BTRAN 10.27%, factor update 8.97%, FTRAN 5.69%, nucleus 5.54%, fresh factor 5.08%, candidate-list update 2.41% and auxiliary factor 1.50% self; native realloc/malloc/delete[] are 2.59/2.03/1.37%. Scope includes launcher, parsing, startup and validation; these self fractions are not inclusive or removable bounds. Row annotation samples short scalar Devex dots, the recurrence/norm check, invariant mode/product branches and the basic-variable scan. Primal control still samples ranked-product scatter. Prior short-dot/ranked-output filtering rejections remain applicable. Exact build/hash, runner, commands, complete results, direct count/event audit, reports and both annotations: `target/perf-primal-0436/`.
- primal_5959.primal-r, `zxxynqtw` (stale after native extended reduced costs; PGO with all saved flags plus frame pointers; full profiling parity 987/987; twenty solves capped at 11000 iterations, all numerical_failure in 10780 iterations and 5381 refactorizations; `target/perf-primal5959-current.data`): primal control 14.30%, F80 add 12.74%, F80 multiply 5.64%, set_d_from_y 5.36%, malloc/free 7.72% self. Refactorization is 28.35%, set_d_from_y 9.62%, recompute_primal 8.18%, extended BTRAN 7.84% inclusive. Reprofile remaining extended BTRAN, eta products, primal control and refactorization before another change.
- fx_quad_obj_5371.dual-r, `prxvqvrr` (preceding compact-RHS/native-TCMalloc/PGO release CLI, stale after reduced-cost-result reuse; ADLITTLE first; 64 dual-release solves capped at 5500, all infeasible in 5319 iterations with matching complete outputs): period-8M instruction self profile has 6409 samples, zero loss/throttle events and represented instructions 99.61% of direct 51.475B whole-batch count. Primal control is 9.19%, extended BTRAN 8.16%, FTRAN 6.93%, set_d_from_y 6.07%, ranked products 4.88%, factor update 4.51%, unsigned 128-bit division helper 3.76%, L stage 2.96%, auxiliary factor 2.92%, F80 division 2.47%, to_f64 2.43%, multiply 2.08%, add 1.93% and sub 1.65% self. Native malloc/delete[]/realloc are 3.85/2.90/2.32%. Scope includes launcher/parsing/startup/validation; inlined arithmetic is charged to callers and these fractions are not removable or solve-only bounds. Raw data, validated results, exact build/hash/commands, direct counts and event audit: `target/perf-rhs-current/`. The preceding profile remains in `target/perf-packed-number/profile8m/`; its period-2M exploratory capture was throttled and is excluded.
- degen3.dual-c, `nzspmlyk` (all saved PGO flags plus frame pointers; full profiling parity 987/987; five solves capped at 2300, all optimal in 2212 iterations and 24 refactorizations; `target/perf-degen3c-taugather.data`): iteration 10.43%, FTRAN 8.64%, L stage 7.38%, needed FTRAN 6.57%, numeric LU 4.76%, BTRAN 4.36%, update 4.20% self. Iteration 49.92%, refactorization 21.27%, needed FTRAN 16.05%, FTRAN 14.05%, L stage 9.79% inclusive. The new term gather is 0.94% self. Current annotations `target/perf-degen3c-taugather-{needed,gather}.txt`; inspect remaining joint FTRAN/L work and iteration control while accounting for rejected traversal/estimate trials below.
- degen3.dual-r, `rxqrkvwm` (stale after eta gather-buffer reuse; frame-pointer build passes full 987/987 parity; five solves capped at 2300 iterations, `target/perf-dual3-rest.data`): L stage 9.24%, FTRAN 8.62%, ftran_needed 7.61%, iteration 7.02%, BTRAN 5.83% self. ftran_needed is 21.69%, L stage 13.85%, FTRAN 15.30%, refactorization 23.59% inclusive. All five optimal in 2212 iterations, 24 refactorizations. Dense L-stage output compaction is another 2.52% self; reserving it failed the timing gate below. Eta gathers now reuse one buffer per call. Reprofile before choosing more L-stage numeric scatter or output conversion work (`target/perf-dual3-rest-lstage-annotate.txt`).
- mode_switch.dual-r, `yxpknytr` (stale after loop alignment; PGO release CLI, self-only profile; three solves capped at 12000 iterations, all optimal in 11393 iterations and 17 refactorizations; `target/perf-densecols-mode-before.data`): factor_structure 12.81%, usize maximum reduction 9.20%, post_inner 7.99%, FTRAN 5.96%, BTRAN 5.85% self. The FTRAN annotation samples two full marker scans; both hot loop heads are on 64-byte boundaries in this binary. The landed loop-alignment setting keeps these marker heads at boundaries as nearby code changes. The earlier marker-block trial is recorded below. Annotate the maximum reduction's caller before assigning it to a field.

- quad_2717.dual-r, `vypkuxry` (pre-packed-predicate accepted compact Devex increment reuse; release PGO/TCMalloc, ADLITTLE/fx first): sixteen solves capped17000 are optimal in16010 iterations, with complete results equal within the batch. Period4M has7231 samples (7201 Rose,30 launcher), zero loss/throttle and represents99.87% of the direct28.961B batch instructions. Whole-batch self shares at Rose symbols are extended BTRAN13.40%, primal control10.52%, native reduced costs5.17%, factor update4.33%, FTRAN3.68%, cached row update3.62%, F80-to-double3.50%, cost-change WORK2.61%, software subtraction2.23%, addition2.16% and integer division helper2.07%. These are sampled locations, not inclusive or solve-only bounds. The old row-update share is superseded after caching. BTRAN stack storage is verified but has no useful runtime gain below. Current BTRAN and factor-update annotations, saved build/hash, bounded runner, full outputs, counters and event/command audit: `target/perf-quad-devex/`; current PEROLD update annotation: `target/perf-ordinary-devex/update.asm`.

- quad_2717.dual-r, `ntpwppnk` (preceding compact Devex increment reuse; bounded-normal conversion, native release PGO/TCMalloc, ADLITTLE/fx first): sixteen solves capped at 17000 are optimal in 16010 iterations with complete outputs equal. Period 4M records 7682 samples with zero loss/throttle and represents 99.87% of directly counted 30.768B batch instructions. Rose-only self shares are extended BTRAN 13.69%, primal control 9.78%, extended row update 8.80%, native reduced costs 4.30%, leading F80-to-double copy 2.88%, integer division helper 2.54% and F80 division 1.43%. The recording includes launcher/validation outside Rose shares; CLI parsing/startup remain. These are sampled locations, not inclusive/solve-only bounds. Release division assembly has two quotient-helper calls before packing; the distinct one-quotient normal-result trial is timing-inconclusive below and must not be extended. Exact build, bounded runner, complete results, direct counts, recording, audit and division annotation: `target/perf-quad-3fcc/`.

- quad_2717.dual-r, `vpqyvlsq` (preceding bounded-normal conversion; native release PGO/TCMalloc, ADLITTLE/fx first): sixteen capped17000 solves are optimal in16010 iterations with complete outputs equal. Period4M has7935 samples, zero loss/throttle and represents99.91% of direct31.767B batch instructions. Rose-only self shares: extended BTRAN13.41%, extended row update10.21%, primal control9.50%, native reduced costs4.78%, leading F80-to-double copy4.68%, software subtraction3.43%, factor update3.13%, FTRAN3.07%. The recording also includes launcher/validation, excluded from the Rose-only shares; CLI parsing/startup remain. These are sampled instruction locations, not solve-only or inclusive bounds. This profile motivates the qualified conversion above; refresh the remaining changed paths before another kernel choice. Build, runner, complete results, counters, recording, audit and conversion annotation: `target/perf-quad-260d/`.

- quad_2717.dual-r, `ztrylnzw` (accepted native-TCMalloc/PGO release CLI, ADLITTLE first; sixteen dual-release solves capped at 17000, all optimal in 16010 iterations and 3022 refactorizations with matching complete outputs): period-4M instruction self profile has 8270 samples, zero loss/throttle events and represented instructions 99.95% of direct 33.095B whole-batch count. Extended BTRAN is 12.47%, ranked products 8.81%, primal control 8.81%, extended row update 7.59%, native extended reduced costs 4.41%, factor update 3.60% and FTRAN 3.57% self. Software F80 to-f64/sub/mul/add/div remain separately attributed, with inlined operations charged to callers; native malloc/delete/realloc are 2.85/3.01/1.48%. Scope includes launcher, parsing, startup and validation; these fractions are not inclusive or removable bounds. Dual true_d is not a leading symbol here. Exact build/hash/commands, complete results, counters, raw event audit and self report: `target/perf-quad-0436/`.

Allocator ancestry screen on accepted `kpttrrms` source: the matching-PGO,
TCMalloc line-table/frame-pointer profiling build passes full parity 987/987
and complete DEGEN2/PEROLD result checks. DEGEN2 whole-batch instructions
rise 11.28% versus the accepted CLI, so these samples locate ancestry only, not
accepted-release cost bounds. Two sixteen-solve CPU-4 PEROLD batches, cap 1300,
periods 4,000,000 / 3,791,733, have 4492/4736 Rose samples and 341/400 allocator
self samples. Rose ancestors resolve for 341/386 allocator samples. Disjoint
routes among allocator samples: refactorization 152/186 (44.6/46.5%), iteration
104/115 (30.5/28.8%), other post-inner 51/44 (15.0/11.0%), setup/other 34/55.
Represented batch instructions are 99.83/99.94% of the direct 16-solve count;
both have zero loss/throttle events. Plain-release DWARF on DEGEN2 resolves only
Rose ancestors for only 40/92 allocator samples; increasing stack capture 8K to 32K does not help.
Frame-pointer stacks can omit the immediate allocator caller: disassembly shows
the dominant observed `factor_structure_r` return site calls `nucleus`, rather
than allocating directly. Other frequent observed ancestors include fresh-factor
installation, FTRAN/L and dual refactor/controller work. ELF file offsets are
translated through PT_LOAD before source mapping. No new solver speedup follows.
Exact builds, full profiling parity, all results, raw stacks/counters, two-period
audit, ancestry grouping and verified disassembly: `target/perf-alloc-ancestry/`.
The accepted release CLI and build metadata remain byte-identical.

Nucleus chain-maintenance diagnostic on accepted source (temporary standard
`nucleuswork` dump, removed): complete bounded dual-release results equal the
accepted CLI on ADLITTLE/DEGEN2/PEROLD, caps 300/640/1300. The 9/16/27 nucleus
calls shift 2,103/86,995/1,424,407 IDs through front insertions and copy
2,832/33,771/541,599 IDs into chain snapshots. PEROLD cleans 184,154 IDs,
removing 40,353 dead IDs; elimination scans 350,277 IDs, performs 88,400 existing
numeric updates and computes 26,807 fill terms. Eleven PEROLD nucleus calls exit
through dense handoff; the small models have none. These are operation/movement
counts, not cycle costs or a linked-list speedup proof. Original source, accepted
CLI and metadata are restored exactly. Raw dumps/results/source/build and runner:
`target/perf-nucleus-chains/`.

Accepted-release disassembly has a 31,086-byte ordinary nucleus function and
0x6c8-byte local stack frame, with dense LU itself already separate. Matching-PGO
frame-pointer annotation (238 self samples) places about 18.9% of local sampled
instructions at stack-addressing instructions; this neither attributes stalls nor
bounds removable accepted-release cost. The kernel contains the full dense
handoff setup/output path alongside sparse elimination. Exact entry disassembly,
annotation screen and ancestry evidence are retained beside the chain counts.

Predecessor barrier instruction refresh (`kpttrrms`, accepted CLI13e, CPU4;
eight DEGEN3 barrier-release solves capped at100, complete outputs agree):
2,638 period-8M samples represent 99.89% of the direct 21.128B batch count,
with zero loss/throttle events. Rose-only self attribution puts numeric_factor
31.24%, DGEMM14.82%, DTRSM10.58%, order7.77%, partition3.64% and integer
quicksort3.34%. CLI read/startup and launcher/validation are separate scopes;
these are instruction locations, not inclusive solve-time or removable bounds.
BORE3D and DEGEN2 bounded results are retained too. The already full-parity
qualified line-table/frame-pointer CLI12603 has matching complete outputs on
all three models. Its 633 numeric-factor self samples locate child-update
scatter, rounding tests and dense LDL loops; codegen differs from release.
The fixed child-rounding dispatch is tested and rejected below. Raw profiles,
counts, complete results, build provenance and annotation: `target/perf-barrier-refresh/`.

F80 add-shape diagnostic on accepted production source (temporary test-only
counters, removed before candidate freeze): ADLITTLE has no finite-nonzero
F80 additions; fx_quad_obj_5371 counts 1,044,229, including 156,186 equal-exponent
opposite-sign cases; Quad_2717 counts 2,785,649, including 651,694 such cases
(14.96%/23.39%). Nonzero cancellation results number 152,586/625,468; adjacent
exponent/opposite-sign cases number 43,258/247,975, without a ratio predicate.
These are direct-solve operation frequencies, not cycle/removable-time shares.
All bounded statuses/iterations/refactors match the CLI; optimal ADLITTLE and
Quad objective/primal/dual/reduced-cost fields match exactly. Infeasible fx CLI
suppresses those fields to null, so no numeric-vector assertion is made there.
Original source, installed CLI and PGO metadata are restored exactly. The
separately tested normal-cancellation implementation is rejected below; counts
alone do not justify retrying it. Source, validating test, raw values/counts,
summary and original Sub annotation: `target/perf-f80-add-shapes/`.

Extended Devex ratio diagnostic (temporary test-only trace, removed before
candidate freeze): fx_quad_obj_5371 has 5,270 mode-2 row calls and 100,414
extended mode-2 weight updates, versus 11 normal-mode calls and 221 ratios.
Quad_2717 has 15,988 normal-mode calls and 285,143 ratios; six mode-2 calls
have no extended weight updates. ADLITTLE has none. Original accepted assembly
retains both `(wq_q / a_rq) / a_rq` Div calls inside the mode-2 loop; GOT 364550
resolves the two sites to F80::Div. The lazy coefficient trial below is
unapplied after an inconclusive full-session regression control.

Reset/owner-aware test-only traces retain provider Arcs during measurement to
prevent identity reuse, identify Solver starts explicitly, and capture full
`(a_rq,wq_q,pq[],tau[])` keys plus actual used-column masks. Quad has eight
owners with five distinct keys each. Eight FIFO entries per owner save
284,403/285,143 terms (four entries save only 20,223); primal5959 primal/dual
save 168,720/170,311 and 167,264/168,455 terms. Fx has no normal-mode hits.
These measurements justify the current bound; the earlier global scalar-only
512-entry simulation did not model ownership and is superseded. All keys are
finite with nonzero divisor in these traces; the production guard checks them.
Optimal ADLITTLE/Quad full numeric fields and all status/iteration/refactor
counts match the accepted CLI. Nonoptimal numeric fields are CLI-suppressed.
Instrumentation is removed before qualification. Counts, exact ordered streams,
Solver/owner scope, source snapshots, validating tests and restoration proof:
`target/perf-devex-reuse/`; original mode-2/scalar ratio diagnostic and assembly:
`target/perf-devex-ratios/`.

F80 division-shape measurement uses direct bounded ADLITTLE/fx/Quad/primal5959
dual-release solves. Only test-only counters are added; the independent Devex
shadow is disabled during counting to exclude duplicate oracle arithmetic.
All status/iteration/refactor counts match accepted runs; full optimal objective,
primal/dual/reduced-cost bits equal the accepted CLI. Total divisions / zero
numerators are0/0,417,598/25,396,599,835/131,998 and817,754/172,580. Every zero
is eligible for an exact finite-zero shortcut; no denominator zeros or nonfinites
occur. These are operation counts, not costs or proof of removable calls.
Accepted release BTRAN dense U already inlines the numerator-zero test before
its two quotient helpers. Sparse input uses an out-of-line Div call; its zero
frequency is not isolated. Do not credit all22.0/21.1% Quad/primal zero counts as
ABI savings. Temporary counters and shadow changes are removed and source/CLI
restored exactly. Source snapshots, shape summary, full numeric validation,
bounded command/build, accepted assembly/relocations and restoration scope:
`target/perf-f80-div-shapes/`, `target/perf-quad-devex/btran.asm`.

Exact scalar-quotient working sets include full raw sign/exponent/mantissa of
both finite nonzero inputs and span every internal Solver/factor reset within
one CLI-equivalent solve. Test-only Devex shadow arithmetic is excluded; all
status/iteration/refactor counts and optimal full numeric fields match the CLI.
FX/Quad/primal5959 eligible counts are392,202/467,837/645,174, with2342/746/4470
distinct pairs. A predeclared four-way hash with512 entries hits389,578/407,488/
638,531; a full FIFO512 simulation hits389,759/466,068/638,567. Bounds16–512
and raw scope/source/build remain in `target/perf-div-reuse/`. These operation
streams are diagnostic, not accepted-release cycle costs. The corresponding
global cache implementation is rejected below despite high reuse.

Power-of-two divisors with normal-result range occur25,040/71,493/63,624 times.
Zero first remainder occurs25,185/73,498/63,993 times, including equal mantissas
3738/17,421/8303 times. These counts are too small to treat cheap integer
shortcuts as the leading unmeasured cost; no shortcut is implemented. Original
dense-U assembly still calls the second quotient helper after forming a zero
remainder. The remainder classifier's exponent-difference field is explicitly
not a boundary-complete normalized-result predicate. Bounds, all statuses and
optimal numeric-vector checks, exact source/commands and restored artifacts:
`target/perf-div-powers/`, including `remainders/`.

Owner/column-aware BTRAN division measurement uses monotonic test-only fresh
factor IDs retained across Clone; IDs are never pointer-based or reused. Only
the three diagonal-division sites are tagged, and both raw operands remain in
every key across updates. FX/Quad/primal5959 have59,284/231,284/222,780 eligible
divisions across22/14/66 owners and166/96/428 owner-columns, with455/492/1455
distinct pairs. Maximum per-column cardinalities are8/14/12. FIFO bounds1/2/4/
8/16 give Quad hits42,404/53,950/79,948/148,158/230,792; sixteen-entry FX/primal
hits are58,829/221,325. Thus16 entries nearly eliminate misses in the measured
column scope. Both per-column implementations are now verified but lose
runtime at counters, as recorded below. Exact records for all bounds and models,
full scope/owner rules, numeric validation, source snapshots and restoration:
`target/perf-diag-div-reuse/`. Temporary hooks are removed and accepted source,
CLI and metadata restored exactly.

Divisor-only BTRAN measurement retains the same owner/column lifetimes and
three sites but keys only the normalized denominator mantissa, for reusable
integer reciprocal parameters rather than quotient results. FX/Quad/primal5959
have202/127/562 distinct keys, at most3/4/4 per column. A four-entry per-column
FIFO hits59,082/231,157/222,218 of59,284/231,284/222,780 eligible queries.
One-entry hits are49,576/213,858/198,971 and two-entry hits56,490/220,587/206,724.
All statuses/counts and optimal full numeric fields agree with the accepted CLI;
nonoptimal solves expose only status/counts through that CLI. No production
arithmetic changes. Instrumented source, scope/validation, all bounds and restored
source proof: `target/perf-diag-den-reuse/`.

A standalone release probe reuses the unmodified original F80 module as oracle
and computes a normal quotient with a precomputed exact integer reciprocal.
For normalized d and B=2^64, R=floor((B²-1)/d)=B+v gives a multiplication-based
quotient estimate at most one below the exact quotient; the exact remainder
supplies one correction and nearest-even rounding. Zero, exceptional,
unnormalized and boundary outputs retain original division. All2,000,050 integer
quotient/remainder comparisons and1,037,249 F80 comparisons pass, including
random full significands/exponents, boundaries, nonfinites and captured QTAU
values. Eight fixed alternating CPU-4 pairs each perform262,144 operations with
parameters already computed; both variants have identical nonzero checksums.
The arithmetic-only median ratio is approximately0.45. This excludes cache,
parameter construction and solver costs and is not an official timing result
or an accepted solver speedup. The integrated designs below preserve the
arithmetic but do not qualify for adoption. Complete source/original hash, own integer
proof, compiler metadata, raw results and summaries: `target/perf-preinverse-probe/`.

Accepted `sozuzvuv` native-release instruction refresh keeps the qualified CLI
unchanged (SHA db2bff5a, fixed39 PGO/TCMalloc, profile b549084b). Eight SC50A
controlled solves capped300 precede sixteen DEGEN2 controlled/PILOT4 release/
PEROLD dual-release solves capped800/800/1300. Every solve is optimal and complete
results/exit codes agree with the accepted full-benchmark outputs where present.
DEGEN2/PILOT4/PEROLD finish531/690/1273 iterations with16/16/27 refactors.

At period4M, DEGEN2 has1075 samples and98.84% represented instructions, PILOT4
1971/99.61%, PEROLD4309/99.88%; all have zero lost/throttle/unthrottle events.
SC50A's47 samples/91.50% coverage are smoke only. Fixed-period child endpoints
can underrepresent; no sample totals here imply exact removable costs. The
whole-batch denominator includes launcher, parsing, startup and validation;
shares at Rose symbols are self locations, not inclusive or solve-only bounds.

DEGEN2 self shares: primal control14.51%, BTRAN9.95%, candidate-list update6.88%,
packed row update5.58%, nucleus5.12%, factor update4.74%, FTRAN3.91%, fresh3.26%.
PILOT4: primal control12.94%, BTRAN10.35%, packed row update8.63%, FTRAN7.51%,
factor update7.10%, nucleus3.91%, fresh3.65%. PEROLD: iteration13.41%,
FTRAN9.72%, L-stage8.10%, factor update7.91%, needed FTRAN5.64%, BTRAN5.29%,
nucleus5.22%, auxiliary factor4.11%, post-inner4.06%.

DEGEN2's candidate-list annotation samples threshold construction/comparison in
`bucket` and membership paths in `fresh.contains`. The current bucket source
scans strict thresholds1..58; no replacement or claimed cost saving is verified.
PILOT4 packed-row annotation and exact symbols are retained too. No solver source
or build changes are made by this refresh. Commands, bounded results, immutable
build/hash binding, direct counts, raw event audits and reports/annotations:
`target/perf-row-accepted-refresh/`. The prior ordinary-primal profile remains in
`target/perf-primal-refresh/` as pre-packed-predicate evidence.

Ordinary ranked-product working sets are now measured with test-only monotonic
immutable-model owner IDs, spanning each bounded solve. ADLITTLE/BORE3D primal
release at cap300 have no ordinary rank_products calls; DEGEN2 release's cap300
prefix also has none. This does not describe their controlled profile or later
DEGEN2 work. PILOT4 primal release at cap800 has367 calls under one model owner,
1,401,403 terms and359 distinct exact ordered numeric keys. FIFO1..16 hits none;
32/64/128/512 hit1/3/6/8 calls, saving only5/19/47/57 terms. Full-result reuse
has negligible measured scope. The canonical key is the sorted (rank,row) order
projected to (row,full vi bits); ignored values and literal rank numbers cannot
change this arithmetic. Complete raw coefficient bits are retained per owner.
All statuses/counts and optimal objective/x/dual/reduced-cost bits match the
accepted CLI; nonoptimal numeric vectors are CLI-suppressed. Diagnostic reversal
recovers the original source exactly apart from final whitespace.

Ignoring numeric bits yields203 row orders. FIFO1/2/4/8/16 hits117/135/145/152/
156 calls, covering553,542/635,322/684,823/704,634/704,690 terms. Additional space
has almost no term benefit. Among366 adjacent queries,117 orders are identical;
72 changed orders need at most two insertions/deletions, while the median changed
order needs12 edits and has a261-row changed span. These are structural reuse
measurements, not construction-cost savings. Exact streams, every bound,
transition analysis, source/build/scope and validation: `target/perf-rank-reuse/`.

A standalone release gather-plan probe preserves the complete original ranked
product as its independent oracle. Plans group terms by output while retaining
each output's exact original term order; nonfinite input/coefficient cases use
the original scalar traversal. All5984 complete-vector comparisons pass across
16 MXCSR rounding/flush modes, all367 captured PILOT4 queries and boundary/
nonfinite vectors. Eight fixed alternating CPU-4 pairs replay four cold
solve-equivalent streams per variant. With order/plans precomputed and common
output allocation, median gather/scatter arithmetic ratio is0.916. Including
original order collection/sort plus actual eight-entry FIFO lookup/construction,
median plan/original ratio is3.091. The rolling checksums agree and are nonzero.
The standalone ThinLTO probe has no PGO and excludes solver integration; neither
ratio is an official solver screen or a universal lower bound. The initial
unpinned smoke is retained separately and excluded from these ratios. No plan
cache is installed. Original source proof, complete captured inputs, code/build/
compiler, raw timings and summaries: `target/perf-rank-plan-probe/`.

Producer-U identity audit: controlled-dual SC50A/ADLITTLE k0 has no nucleus or
tail columns. The first nontrivial SC50A dual factor is k1 (5 nucleus rows,
10 tail columns; existing geometry marks 4 guarded columns/2 structural visits).
The older perf-fresh SC50A/ADLITTLE case0 fixtures are primal-controlled.

The current two-call SC50A dual-controlled refcap has 51 p68 endpoints, resolving
the old m-entry endpoint gap in the legacy LIM400 cache. An existing no-build
trace_dump (SHA 763e5f2d…e02, system-TCMalloc/align64/accepted-profile fingerprint;
separate from the qualified CLI) capped at 8 iterations agrees with k1: input
basis-column IDs, pivot rows, 50 diagonal bits and 29 stored L bits/index order;
13 final U values agree by pivot(column,row). All 92 values are finite/nonzero;
physical U storage order and zero/nonfinite behavior are untested. The scratch
adapter supplies cached records to unchanged compare.py cmd_lunum alignment;
normal lunum rebuilds and still prefers the legacy cache. No new alignment rule
or broad source/parity identity for that example binary is claimed.

Producer/pre-overlay observation now has a nontrivial paired witness and a
counterexample. SC50A dual-controlled factors 0/1 match the existing reference
input columns bit-for-bit and match pivot row/column orders. Factor 1 has 5
Markowitz phase0 rows, 4 phase3 special-row positions and 5 phase4 singleton rows:
i21=40, loff=45, tail 10. Earlier geometry's “nucleus 5” means the Markowitz block;
it is different from RefStruct.nuc's phase3 block. SC50A's four producer-supplied
entries match pre-overlay raw bits and remain the installed values, but nine
other replay entries are absent from that producer subset. Existing final
reference U agreement covers 13 entries, not a complete producer replay proof.

Across 43 observed/reconstructed small/DEGEN2 factors (30 nontrivial; includes
live/cached overlap), 880 supplied values match replay bits. The 26 cached cases
are SC50A/ADLITTLE/BLEND/BORE3D; DEGEN2 is bounded at 300 controlled-dual iterations
and 16 fresh calls (eight observed). Small live probes use eight iterations and
two fresh calls; AFIRO hits that factor limit, and some initial factors are
vacuous. No observed sparse-factor retained Markowitz entry is missing from
the producer subset; 2596 earlier-row replay entries are not supplied there.
These observations do not establish general numerical identity.

A deterministic 51x51 double matrix exercises the actual dense handoff at its
minimum permitted size (`remaining>50`). At the first difference, tail step 8,
incoming column 42, pivot-row step 6 (original row 44), producer bits are
`0xbf51473431609794` (-0.0010545739837737088); replay/pre-entry bits are
`0xbf51473431609793` (-0.0010545739837737086): one ULP. The original final
column installs the producer bits. Observation stops at that column, with no
later column executed. The dense producer uses lane-wise dot kernels while
replay is scalar; no exact intermediate-cause model is claimed. This is a
Rose-internal counterexample, not a new reference-solver comparison.
Unconditional numerical producer/replay identity is false; preserve the original
overlay and do not use the matching sparse cases as a global proof.

A separate data-flow lead does not require identity: with strictly lower L,
a contiguous producer-covered suffix [cut,t) cannot feed an earlier output. If
every last-effective supplied value in that suffix survives tolerance retention,
all its final outputs are overwritten. Skipping that region's input/graph/numeric
work could preserve the earlier replay and original merge/retention/order logic.
Observed sparse-prefix coverage marks 254/878 columns, 933/2210 structural graph
visits and 502/671 executed update products in covered suffixes. The nine-column
dense counterexample prefix marks 8 columns, 2080/2130 graph visits and 120/120
update products. These are unweighted observation/structural counts, not saved
time or a qualified speedup. The arithmetic-only shortcut and the complete
contiguous-suffix graph-pruning trial both remain unapplied; neither timing
screen establishes a speedup. Do not repeat or extend either screen.

The graph-pruning trial preserves original merge/retention/order and WORK, using
strictly lower L and a last-write, above-tolerance suffix certificate. Its frozen
complete constructor matches every field and raw f64/F80 value on 26 cached
small factors per precision, randomized/extreme/nonfinite/duplicate cases,
nontriangular fallbacks, the dense51 mismatch and all 16 SSE modes. All ten
FTRAN/L tests and fast940 parity pass under matching fixed39 PGO/TCMalloc.
Eight bounded ADLITTLE/BORE3D/DEGEN2 and three PILOT4/PEROLD result/exit pairs
agree. PILOT4 instructions/cycles change -1.068/-0.700%; PEROLD -0.494/-1.108%.
The predeclared PILOT4 primal-release screen (CPU4, cap800, seed32217,
32AB+32AA) gives median -0.189%, CI [-2.199%, +1.965%], 17/32 improving,
sign p .430; order medians +1.105/-0.318%. Healthy AA spans
[-1.380%, +0.641%], p .597. Gain/sign/order/noise separation fail: stop without
seed32289, extension or full candidate gates. Source, accepted CLI/metadata and
profile are restored exactly; restored formatting, fast940/full987 and all
254 release tests pass (22 default ignored). Candidate/frozen oracle, raw-field
proof, build/freezes, counters, complete results, timing and restoration:
`target/perf-suffix-replay/`.

Original arithmetic/control/WORK reverses byte-for-byte after removing three
test-only observer hooks and raw-bit accessors. All probes are removed; accepted
source/CLI/PGO remain exact. Formatting, fast 940/full 987 parity and 254 release
tests pass. Raw producer/pre-overlay/final data, deterministic matrix/seed,
exact input/pivot scope checker, geometry definitions, first-difference stop,
source/binary/build provenance, suffix-count analysis and restoration:
`target/perf-pre-overlay/`.

F80 ABI audit and software-layout feasibility are complete. Actual accepted
multiply/divide use operand pointers and an indirect result pointer; actual
two-word software operators use scalar registers. The same16-byte value size
holds in both. The new metadata word stores the original unbiased i32 exponent
in low32 bits and sign in bit32. Mechanical constructor/field-read edits reverse
exactly; all original arithmetic/rounding/packing expressions remain unchanged.
Debug retains the original fields/order. A complete untouched original type is
the independent oracle, including all internal sign/exponent/mantissa bits.

The standalone release probe passes1,581,056 arithmetic/comparison/conversion
checks, boundary/zero/nonfinite/captured-QTAU values and16 MXCSR modes. Eight
fixed alternating CPU-4 pairs use4096 input pairs×16 repetitions, with equal
nonzero rolling checksums. Ordinary operator-dispatch new/old time ratios for
add/mul/div/mul-add are1.193/0.715/0.921/0.898; forced opaque call boundaries give
1.147/0.426/0.969/0.824. These are mixed isolated results without PGO, not solver
speedup claims. The two matching-PGO integrations are recorded below; neither
qualifies. Original/transformed code, reversible edit records, full-field checks,
compiler/flags, raw timing and real operator disassembly:
`target/perf-f80-two-word/`. Earlier toy/accepted-ABI evidence remains in
`target/perf-f80-abi-audit/`; prior packed-native/x87 probes are a different
encoding/CW/TLS/conversion design.

The optional-layout cost is independently checked: Option<F80> grows 16→24 in
plain two-word storage, and optional pairs 32→40. DevexDeltaEntry.values and
primal inputs use these types. Tagged NonZeroU64 software metadata retains the
original 16/16/32 value/optional/pair sizes; bit 63 is always set, low 32 retain the
unbiased exponent, bit 32 retains sign. Safe constructors, sign toggle, Copy and
Clone preserve the invariant. None stays distinct from Some(logical zero).
All original numerical arithmetic remains unchanged; Debug retains field order.
The independent full original oracle passes 1,581,056 counted checks, plus 736
sqrt and optional/default/tag assertions, boundary/captured inputs and 16 MXCSR
modes. No raw-zero F80 initialization was found in the source audit.

Eight fixed CPU-4 standalone pairs remain mixed: direct add/mul/div/mul-add
ratios 1.199/0.598/0.936/1.030; opaque 1.154/0.665/0.941/0.990. Value-vector
initialization ratios for lengths 7/32/410 are 1.116/1.328/1.005, optional vectors
0.759/0.737/0.710 (allocation included, default allocator, no PGO). These are
isolated feasibility measurements. Matching fixed 39 PGO passes 10 F80 and 6
FTRAN/L tests, fast 940 and all bounded complete result/exit pairs. Actual
multiply/divide still use scalar register arguments. fx instructions/cycles
rise 0.323/3.084%; Quad rises 0.008/3.656%, with branch misses+19.812%.
Stopped before seeds 29217/29289; no full candidate gates or speedup claim.
Tagged encoding improves optional sizes but gives no useful target counter gain;
no attribution of the plain layout's timing deficit to Option expansion follows.
Exact source/invariant proofs, standalone results, assembly, matching build,
source freeze, tests, bounded results/counters and restoration:
`target/perf-f80-niche/`. Prior toy layout/ABI evidence:
`target/perf-f80-two-word/option-layout.*`, `target/perf-f80-two-word/optional-layout-audit.json`.

Auxiliary-factor exact-input reuse is measured across isolated bounded solves.
AFIRO, SC50A, BLEND, ADLITTLE, BORE3D and DEGEN2 primal/dual each have no repeated
keys (4–13 calls per case); PILOT4 primal has 17 distinct calls and PEROLD dual 27.
Compact fx_quad (cap 5500) has 1786 calls/25 unique keys; Quad (cap 17000) has 3024/18.
An eight-entry FIFO hits 1760/3006, versus four-entry 1758/3000; every hit succeeds.
Keys retain m, all ordered columns/row indices/value bits including duplicates
and signed zeros, floor bits, and None/Some prescribed rows. Complete CLI result
objects agree with the accepted binary for all 16 cases. Test-only capture hooks
reverse exactly to the accepted source; no new reference comparison is claimed.

Original construction versus deep-copy fixed 8-pair CPU-4 probes give ratios
0.142/0.198 on fx/Quad. Full eight-entry key lookup, miss insertion and factor
copy cost 0.261/0.262 times reconstruction, with every stream output checked against
the unchanged constructor. PGO flags match the accepted source; libtcmalloc is
explicitly retained like the CLI. These isolated component timings include the
inactive test-only capture check and are feasibility evidence, not solver gains.
Source mutation audit shows that fresh numeric fields remain read-only; only
etas and nnz_eta change through updates. Immutable-base sharing is unimplemented.
All input streams, costs, simulations, allocator linkage, exact source proof,
complete results, bounded timings and probe code: `target/perf-aux-inputs/`.

The matching fixed 39-PGO integration caches eight fresh auxiliary constructor
results per thread for m<=50, returning deep owned clones with independent eta
updates. Dumps, nonfinite inputs, large matrices and alternate/unmasked/flush
FP environments retain the original constructor. No arithmetic, caller failure
point or WORK changes. An independent full original LU oracle passes all eight
LU tests (four new), covering full captured streams, key/boundary/failure cases,
random constructors, owned dense/sparse updates, FTRAN/BTRAN continuations and
all 16 MXCSR modes. All six FTRAN/L tests, fast 940 parity and bounded complete
result/exit pairs pass. fx instructions/cycles fall 3.162/1.910%; Quad 1.464/5.062%.

Fixed 32-AB/32-AA CPU-4 fx screen 30217 is inconclusive: median -4.682%, conservative
CI [-6.989,-2.210], 26/32 improving pairs, sign p 0.000268, order medians -4.941/-4.682%.
AA is healthy (CI [-1.539,+2.338], p 1.0), but conservative separation fails.
No 30289 confirmation, sample extension or full candidate gates. The candidate
remains unapplied; no speedup is accepted. Exact candidate/frozen oracle, captured
fixture, original source, source/build freeze, tests, counters, screen and
restoration: `target/perf-aux-cache/`.

FTRAN/BTRAN exact-input reuse is measured through 194,701 calls across 18
isolated bounded cases: primal/dual AFIRO, SC50A, BLEND, ADLITTLE, BORE3D and
DEGEN2; compact fx/Quad; PILOT4 primal, PEROLD dual; and primal5959 primal/dual
at cap 11000. A test-only monotonic matrix token survives Clone and changes on
every fresh construction and at entry to every attempted update. No pointer
identity is used. Full ordered indices, f64 bits or internal F80 m/e/sign fields,
method flags, needed order and both FP controls form exact keys. All 18 complete
CLI result objects agree with the accepted binary; original arithmetic/control
flow/WORK reverses exactly after removing the test-only hooks.

Of 132,997 calls after at least one update, only two repeat an exact key, both
full FTRAN on fx. Needed FTRAN and BTRAN have no nonfresh repeats, even with an
unbounded per-version table. A general mutable-version cache is not justified
by these streams. Fresh BTRAN reuse overlaps the earlier unqualified compact
BTRAN-result cache and does not justify repeating its failed screen.

Full FTRAN replaces all F2 fields, work/parts, last flags, filtered spike and
last_alpha; needed FTRAN replaces only F2 tau fields and tau L/U/fin estimates;
BTRAN sets work/density and conditionally row_cached. Lazy row/transpose-level
initialization is separately tagged at its exact branch. Across repeated keys,
output order/bits and all projected written semantic fields/lazy actions have
zero mismatches. Untouched caller state must remain live; scratch is not payload.
The ftran_keep/basic wrappers must still execute their original restoration and
keep-statistic logic. This is a dependency audit and measurement, not proof of
a complete general result cache.

A separate fresh full-FTRAN lead is favorable: two entries per fresh F80 matrix
version capture 1730/2883 fx/Quad repeats and 5273/5627 primal5959 primal/dual.
Most are force_dense=false pivot-column solves, beyond the accepted pure RHS
provider's force-dense result reuse. Complete key lookup, miss insertion, exact
F2/spike/scalar effect replay and output copying cost 0.111/0.112/0.199/0.198
of original stream time on those four streams (fixed 8 alternating CPU-4 pairs,
16 repetitions, matching accepted-source PGO flags and TCMalloc). Original
component stream medians are 3.901/4.006/5.279/5.386 ms. Factor matrices are cloned
once per stream equally, outside individual queries. All four captures reproduce
exactly on the component build, excluding instrumented elapsed time, and every
component call's full written-field payload equals the captured original.
These isolated measurements are not solver speedups or qualification screens.

The two-entry fresh full-FTRAN effect cache is numerically verified but
unapplied. It shares exact ordered row/full F80-bit/force-dense keys only through
the existing compact fresh template owner (m<50), with independent mutable factor
state. It replays complete F2, spike, work/parts and last-field values; original
wrappers and caller WORK sites remain unchanged. Updates and direct U replacement
detach the child cache before mutation. Dumps, updated/large factors, malformed/
nonfinite inputs and unsupported SSE/x87 modes retain the original kernel.
The complete current kernel is independently frozen; token reversal verifies all
other original constructor/helper/update/WORK and Real/F80 arithmetic unchanged.

Matching fixed 39 PGO passes ten FTRAN/L and nine F80 tests, fast 940 parity and
all bounded complete result/exit pairs. Every library-test hit also executes the
whole frozen kernel on a deep clone and compares ordered full-field result bits
and every original factor field value, including scratch and untouched caller
state. Captured input values, cold/warm clone histories, owned updates/failures,
eviction, raw low bits/zeros/duplicates/nonfinites, all 16 MXCSR modes, alternate
x87 controls and a bounded fx solve with >1000 shadow-checked hits pass.

fx instructions/cycles fall 4.021/3.285%; Quad 1.621/2.369%. The fixed 32-AB/32-AA
CPU-4 fx screen 31217 is inconclusive: median -0.355%, CI [-2.043,+1.454], 17/32
improving pairs, sign p 0.430; order medians -1.338/+0.142%. AA is healthy
(CI [-3.360,+0.425], p 0.215). The gain, significance, order and conservative
noise-separation criteria fail. Stop without confirmation 31289, extension or
full candidate gates; no solver speedup is accepted. Source/CLI/metadata/profile
are restored exactly; formatting, fast 940/full 987 parity and 254 release tests
pass. Exact original/current/frozen sources, input fixture, source proofs,
matching build/freeze, tests, counters/results, complete screen and restoration:
`target/perf-fresh-ftran/`.

Exact source reversal, raw streams, per-version/FIFO bounds, unchanged CLI results,
projected effects, component source/binary/timings, dependency audit and restoration:
`target/perf-version-stream/`.

## In progress

The EPYC release builds now use Rust 1.98.1, LLVM 22.1.8 and genuine TCMalloc. The forward L level scan is the first measured target. The frozen digests remain at their original contents.

## Extended BTRAN feasibility

Feasibility-screen baseline CLI hash: `e28b651026b477f0e5741e881efce27d640e2bf3f330f85a4f76b0e53c8804fc`.
The saved CPU-4 instruction self profiles use eight fx_quad_obj_5371 solves capped
at 5500 and sixteen quad_2717 solves capped at 17000. Every result has expected
status, iterations and refactorizations (infeasible/5319/5 and optimal/16010/3022).
Extended BTRAN is 7.01%/11.39% of batch instruction samples respectively, distinct
from the older inclusive cycle shares. There are 3539/8686 samples at periods
2M/4M; both have zero throttle/lost records. Python launcher/startup remain in
the denominator; represented-instruction audits and dense-branch coverage are
not yet established. Profiles do not attribute out-of-line arithmetic to callers.
Runner/manifests/results, reports, audit and BTRAN annotation are in
`target/perf-btran-native/`.

`cargo run --release --example btran_native_probe` screens the complete dense
L-transpose arithmetic pass, preserving multiply/store/subtract rounding and
including input/output and per-coefficient representation conversions. Its 168
bitwise differential cases cover captured SC50A L coefficients, finite random
F80 low bits, empty/singleton vectors, signed zeros, subnormal/min-normal
coefficients, and multiple transpose shapes. Four alternating CPU-4 timing pairs
on this mixed batch give software 0.0932–0.0960 s versus native 0.0293–0.0309 s,
roughly 3.1x. These are isolated kernel times, not solver speedups or a prediction
for complete BTRAN. Raw `probe-cpu4.log`, `probe-expanded.log`, and binary hashes
in `audit.json`; ordinary release probe build leaves the accepted CLI untouched.

The probe remains outside solver dispatch. A dense-only integration was tested
and rejected on solver timing (see the rejected-trial entry). It guarded the x87
control word, nonfinite inputs/coefficients/products/results, left caller input
untouched on failure, and preserved each WORK addition. Native helper and frozen
complete-BTRAN differential tests, including f64/F80, sparse/dense/direct inputs,
row etas, subnormals and exceptional-value fallback, passed; fast parity passed
940/940. Eight CPU-4 AB/BA result pairs per fixture matched all output fields
except seconds. This does not prove full parity or profitability for a complete
native BTRAN. The dense isolated timing remains a feasibility result only.

### Native conversion and operation-boundary audit

Small isolated dense U/eta/L tests use 80 finite cases with synthetic U/etas,
SC50A captured L and sizes8/11/32. Every output retains full raw F80 fields.
Eight alternating CPU4 pairs show old on-demand packed code +11.09% versus
software, prepacked coefficients +7.55%, rebuilding the cache each call +24.32%.
Combining multiply/store/subtract in one assembly block preserves both roundings
and improves the on-demand component by15.36% versus software. These are
standalone pinned-release component times, without solver PGO or ownership/cache
maintenance costs; they do not predict solver speedup.

The matching fixed39-PGO/TCMalloc integration keeps the complete original BTRAN
as an independent oracle and unchanged fallback. It uses canonical normal/zero
eligibility, rejects exceptional/subnormal intermediates, preserves all WORK
additions, output ordering, row-cache/density fields and lazy inverse-row state.
All three new differential tests pass, including full raw f64/F80 outputs,
26 cached small factors per precision, duplicates/extremes/low bits, all16 SSE
modes and unsupported x87 mode dispatch. Fast940 parity passes. Eight bounded
ADLITTLE/fx and three Quad complete result/exit pairs agree. Instructions change
-0.052/-2.876/-3.844%; cycles +3.518/+5.755/+4.842%. Reject before seeds33217/
33289 or full candidate gates. Accepted source/CLI/metadata/profile are restored
exactly; restored formatting, fast940/full987 and all254 release tests pass.

Short CPU4 user-cycle profiles of eight fx solves each (cap5500) have no lost or
throttled events. The candidate helper has179 samples (6.89% batch self), with
remaining BTRAN2.54%; baseline BTRAN8.32%. Samples fall around multiply/divide,
ten-byte loads/stores and surrounding representation traffic. This is sampled
location evidence, not precise stall attribution or comparable absolute costs.
The initial candidate profile attempt failed before solver execution because
the saved executable lacked its execute bit; it is retained as invalid. The
`after-valid` recording includes all eight validated results.

A separate register-product component checks1,048,576 raw full-precision triples,
including exponent boundaries. Under guarded nearest/64-bit x87 mode, a normal
finite product exponent certificate permits keeping the rounded product in a
register; other products use the stored-product path. All80 complete component
cases also agree. Eight additional alternating pairs give stored-product
on-demand -16.18% versus software, register-product -16.54%; the extra register
variant is only0.44% faster than the stored component and has no integrated or
qualified solver gain. Neither coefficient caching nor register retention alone
justifies another solver trial from these measurements. The separate early-drop
measurement and matching-PGO integration below reject the remaining tested
native division-routing variant too.

Frozen sources, arithmetic/component variants, inputs, compiler commands and all
component timings: `target/perf-native-cache-audit/`. Complete candidate/original
source proof, tests/build/freezes, all solver result/counter pairs, assembly,
profiles/event audit, decision and restoration: `target/perf-native-scatter/`.

### Native dense-U division before dropping values

Test-only counters on the prior combined-block native kernel observe5231/30836
native attempts on bounded fx/Quad solves, all successful. Dense U visits57541/
246688 positions;25317/125250 numerators are zero and fail tolerance (44.00/
50.77%). ADLITTLE has no native attempts. These are native-call operation counts,
not whole-solver division totals or removable-time bounds. Original production
arithmetic/WORK are unchanged during counting. Status/iteration/refactor counts
match accepted CLI results; all optimal objective/primal/dual/reduced-cost bits
match. Source is restored after counting and before candidate installation.

The separate predicate-first candidate drops finite at/below-tolerance numerators
before computing their quotient. The original loop drops them regardless of the
quotient, including a NaN diagonal; WORK has no increment for evaluating the
division. All ten focused FTRAN/L tests pass, retaining the independent complete
original BTRAN and full-field/mode tests; a targeted NaN-diagonal case verifies
the newly successful native route. Matching fixed39 PGO/TCMalloc and fast940
parity pass. All80 component cases agree; eight alternating CPU4 pairs improve
the component6.80% versus the prior native ordering, without a solver gain claim.
Eight bounded ADLITTLE/fx and three Quad complete result/exit pairs agree.
Instructions change -0.047/-3.214/-4.480%; cycles +2.264/+0.916/+0.441%.

The predeclared Quad dual-release screen (CPU4, cap17000, seed34217,32AB+32AA)
gives median +3.859%, CI [+3.186%, +5.494%], only4/32 improving; the one-sided
slowdown sign p is9.65e-6. Both order medians slow (+3.524/+4.975%). AA is healthy
(p .860), CI [-0.669%, +0.882%]; the entire slowdown interval exceeds that
uncertainty. Reject, with no34289, extension or full candidate gates. Accepted
source/CLI/metadata/profile restored exactly; restored formatting, fast940/
full987 and all254 release tests pass (22 default ignored). No native BTRAN
implementation is accepted. Source/counting driver and reversal, component
variants/timings, original/complete oracle, matching build/freezes, all bounded
results/counters, screen and restoration: `target/perf-native-drop/`.

### BTRAN branch and stage counts

A temporary `btranpath` dump kind on the accepted source counts arithmetic terms
for calls reaching L-transpose. Captures are bounded at 5500/17000 iterations;
both complete output objects match the accepted CLI except seconds. Numeric
operations and ordering are unchanged. This instrumentation is removed, and the
accepted CLI hash restored. Instrumented times are not benchmarks; fixed-profile
PGO hash-mismatch warnings make this binary unsuitable for timing attribution.

| Fixture / L path | Calls reaching L | U scatter terms | Eta terms | L scatter terms | Rows |
|---|---:|---:|---:|---:|---:|
| fx_quad_obj_5371 dense | 11320 | 123735 | 56464 | 186969 | 11 |
| fx_quad_obj_5371 heap | 16 | 66 | 38 | 56 | 11 |
| quad_2717 dense | 53960 | 463569 | 329434 | 616740 | 8 |

Quad_2717 has no heap/level/partial L calls in this measurement. The dense
L branch has at most 14 scatter terms per call on Quad_2717, 20 on the smaller
quad fixture. Zero-term dense calls number 5780/5 respectively. These are term
counts, not cycle costs: U divisions, graph/list management, packing, allocations
and threshold tests are excluded. Early returns before L produce no record;
stage totals therefore describe only calls reaching L, not every BTRAN entry.
The evidence contradicts extending the native kernel into sparse L as the next
optimization on these targets. A whole U/eta/L packed workspace would amortize
conversions across more work, but must preserve the distinct sparse U expression
`-(v*d)+w`, dense U subtraction, eta assignment/addition, threshold tests, ordered
supports and each WORK accumulation.

Raw captures/results, validated runner, source snapshots and binary provenance:
`target/perf-btbranch/{run.py,summary.json,allstages-summary.json,*-allstages.dump,
*-allstages.json,*-accepted.json,instrumented-ftran.rs,instrumented-dump.rs}`.
Build logs `/tmp/rose-btbranch-build{,2}.log`. No frozen baselines were changed;
no reference comparison or new solver speedup is claimed.

### Maintained eligible ranked-row traversal

The ordinary pivot consumes nonbasic structural products (including fixed), all
slacks, q and leave. A separate cost-change caller retains its complete original
product. A maintained row stream can remove other basic structural entries once
per basis change, avoiding the rejected per-entry eligibility test. Strictly
ascending column IDs certify that each row has no duplicates; duplicate/unsorted
rows retain the original full traversal. Swapping distinct column entries changes
no output's row accumulation order, so all retained arithmetic stays ordered.

Test-only caller captures preserve all original arithmetic/WORK. ADLITTLE,
BORE3D and DEGEN2 capped300 have no ordinary-pivot calls; SC50A capped300 and
DEGEN2 capped800 provide smaller witnesses (8/181 calls, optimal32/690 iterations,
7/17 refactors). PILOT4 capped800 has367 calls, optimal690/16. Its keys match all
367 calls of the cached general ranked-product diagnostic exactly. All bounded
statuses/counts and optimal full numeric vectors match the accepted CLI.

PILOT4 can omit834595/1401403 entries (59.55%). Full basis-position comparison
causes3943 changed slots and56385 entry updates; comparing membership sets needs
only7503 updates. Both complete scan-based component plans still lose roughly
32–35% in preliminary timing. A separate producer-notified plan applies the known
pivot and resynchronizes only at inferred basis-reorder events. Its final frozen
standalone release probe includes construction once per stream, result allocation
and every membership/resync update. Sorted row-order preparation is common and
outside timing; this uses the system allocator, no solver PGO. Eight alternating
CPU4 pairs of64 complete streams give SC50A+408.82%, DEGEN2+40.17%,
PILOT4-18.45%. These are component costs, not solver speedups; generation events
are inferred from the captured predicted post-pivot basis, not integrated with
the solver's refactor/repair lifecycle. Unconditional integration is unprofitable.

All492530 retained captured f64 values and76723 edge-case values agree bit for
bit. Edge cases cover16 MXCSR modes, sizes1/2/8/31, signed zeros/subnormals/extremes,
infinities/two NaN payloads, row orders, basis replacement/reordering and stamp
overflow; duplicate/unsorted fallback eligibility is checked. All diagnostics are
removed; accepted source/CLI/metadata stay exact. Formatting, full987 parity and
all254 release tests pass (22 default ignored). No candidate solver build or
official timing screen is run. Captures/baseline scope, all three plan sources,
preliminary negative timings, final source/binary freeze, all final timings/raw
checks, generation-scope limits and restoration: `target/perf-active-rank/`.

### Ranked-traversal gate and original bounds check

The separate gate admits a query after cumulative dense work reaches eight
immutable-matrix passes; dense means at least half the matrix's entries are
visited. Every sparse fallback still notifies an existing plan of the committed
pivot, or invalidates its generation. Complete64-stream, eight alternating CPU4
component pairs include all predicate/counting/construction/maintenance costs.
SC50A activates0/8 calls, DEGEN242/181 and PILOT4281/367. Changes versus the
original component are +6.04/+73.85/-9.70%; the gate is not broadly profitable.
All492530 retained captured outputs and76723 edge-case outputs agree exactly.
A bounded independent PEROLD primal-release prefix (cap600, iteration limit,
14 refactors) has no calls to this path; it supplies no generation/performance
coverage. Its test-only actual-refactor-count capture is removed and the accepted
source/CLI stay exact. Source, fixed gate, freeze, all samples and validation:
`target/perf-gated-rank/`.

A different flat active-entry bitset retains immutable coefficients, replaces
swap removal/location arrays with per-column bit updates and walks set bits in
original entry order, including duplicates. Preliminary captured-stream costs
are SC50A+202.43%, DEGEN2+59.19%, PILOT4+4.06%; no solver integration follows.
The initial edge fixture incorrectly demanded an unused basic-column output;
it is corrected to the retained-output contract, and76723 edge checks now pass,
including duplicate/unsorted entry order. The earlier timings ran despite that
fixture failure and remain preliminary; no qualified gain is asserted. Original
and corrected fixture/binaries, timing freeze, all samples and limits:
`target/perf-bit-rank/`.

The original ordinary ranked scatter's accumulator access is a separate trial.
`build` enumerates exactly n+m columns into private rwise; production dimensions
and that matrix stay paired and unchanged. Removing only its per-entry result
bounds check retains all entries, numeric expressions/order, callers and WORK.
The complete independent original method is source-identical, and its new
captured-AFIRO plus random/duplicate/nonfinite/two-NaN/all16-MXCSR tests pass.
Matching fixed39 PGO/TCMalloc and fast940 parity pass. Eight bounded SC50A/
DEGEN2 and three PILOT4 complete result/exit pairs agree. Instructions change
-0.059/-0.243/-0.738%; cycles -1.840/+0.372/-0.314%.

The fixed PILOT4 screen35217 (32AB+32AA, CPU4, cap800) is inconclusive: median
+0.729%, CI [-2.617%, +3.146%],14/32 improving, sign p .811; order medians
+1.793/+0.592%. AA is biased (median+1.411%, CI[+0.288%, +2.244%], p .007).
Stop without35289, extension or full candidate gates. Source, accepted CLI/
metadata/profile restored exactly; restored formatting, fast940/full987 and all
254 release tests pass (22 default ignored). The bounds-elision trial remains
unapplied. Full original/source/memory proof, matching build/freeze/tests, all
results/counters, screen and restoration: `target/perf-rank-unchecked/`.

### Pre-packed-predicate primal source-line screen

A no-frame-pointer, line-table profiling build using the accepted PGO flags and
system allocator passes full parity 987/987. Its complete ADLITTLE/PILOT4 results
match the accepted CLI except seconds, but direct instructions rise 8.88/6.44%.
Do not use its source sample shares as release cost bounds. Eight ADLITTLE and
32 PILOT4 primal-release solves are capped at 300/800; PILOT4 is optimal in
690 iterations and 16 refactorizations. Sampled instruction locations identify
primal control, BTRAN/FTRAN, the specialized row-update kernel and factor
construction. No source line in the cost-change region appears in the PILOT4
controller annotation; this is a source-attribution observation, not proof of
zero inclusive cost. Current source inspection also identifies repeated complete
structural-column length sums; replacing these from existing prefix offsets
is tested and rejected as timing-inconclusive below.

Build/manifest, full profiling parity, annotated controller, self report, event
audit and direct instruction/result comparison: `target/perf-costchange-audit/`;
build log `/tmp/rose-costchange-profile-build.log`. Direct counters were taken
after recording, not simultaneously with the sample batches. Represented-event
coverage was not established, so sample totals do not support absolute costs.
The accepted release CLI hash is restored and unchanged.

### Randomized timing control

`perf/paired.py` implements the session review's bounded randomized AB/BA screen
with interleaved same-binary AA controls. It freezes the schedule/decision rule
before warm-up, retains every raw result and sample, checks complete result/exit
identity outside the solve timer, and refuses to overwrite evidence. A positive
screen is not a full acceptance gate; independent confirmation and all correctness
and regression checks remain necessary. Rules are in PLAN.md and the tool's help.

On PILOT4 primal-release (CPU 4, cap 800, 32 pairs of each kind, seed 137), the saved
prefix-count candidate versus accepted binary has median -1.42%, conservative 95%
paired-median interval [-2.41%, +1.23%], 20/32 improving pairs (sign p .108), and
AA interval [-1.37%, +2.27%]. This is inconclusive under the predeclared rule; it
reinforces the original decision not to land the candidate. No samples were removed.
The accepted binary compared with itself (seed 271) reports +0.11%, interval
[-0.07%, +0.38%], and correctly fails the gain screen. Four deterministic rule
checks cover zero effect, a stable gain, an execution-order confound and biased
AA controls. Neither timing session modified solver sources or the accepted CLI.

Exact commands, schedules, hashes, complete result objects and decisions:
`target/perf-paired-nnz-pilot4-137.json` and its `.raw/` directory;
`target/perf-paired-controls/{identical-271.json,identical-271.raw/,statistics.json}`.
The prefix-count candidate remains archived and unaccepted. These screens establish
no new solver speedup or reference timing ratio.

## Candidates

The first structural target on the EPYC is FTRAN/BTRAN on mode-switch. Reuse factor-level traversal metadata and workspace across calls, and maintain reached positions so sparse solves initialize, divide and emit only the entries they must visit. Preserve ascending slot order, bucket insertion order, signed-zero and exceptional cases, all F2/dump state, and the original full-range WORK contributions. A standalone L-level maximum cache, blockwise marker skips and several workspace-reuse variants have prior timing rejections on other models. Measure a combined sparse path against the actual starting allocator and run the full gate before accepting it. The profile identifies costs. It does not establish a gain for this proposal.

The remaining candidates below derive from the Ryzen measurements and their historical profiles.

- Ordinary primal/dual snapshot ratios are in Numbers. Host speed bands can reorder them; use paired old/new measurements to establish gains. The current PEROLD dual profile above refreshes the largest measured ratio; the long-eta SIMD trial below is inconclusive, so inspect another measured cost before another dot variant.
- primal_5959.primal-r: reprofile remaining extended BTRAN, eta products, primal control and refactorization with current PGO flags before another change.
- DEGEN2 joint FTRAN: dense interleaving, sparse-column/dense-tau traversal and individual/combined U-estimate deferral are rejected below. Combined deferral has also been tested with shared runtime-policy kernels, unchanged factor fields and the original column vector return; it reduces instructions without a repeatable PEROLD timing gain. The final mode-switch variants do not isolate the cause of the earlier const-generic regression. Reprofile remaining FTRAN/L traversals before another change. Reference control-flow evidence: `target/perf-pairpack-reference-audit.txt`.
- Forward-L layout evidence: cached SC50A has 10/29 entries in consecutive-row runs of at least two after sorting, none in runs of four. Bounded restored-release BSOLVE captures (ADLITTLE cap 300, BORE3D 500, DEGEN2 640; all reach their normal optimal exits) show no duplicate or non-lower L entries. In DEGEN2's largest captured L matrix, 584/1578 entries are in runs of two and 138/1578 in runs of four; the smaller models have little long-run coverage. A sorted contiguous-run cache remains untested and its constructor/storage cost needs justification before implementation. Geometry/data: `target/perf-lranges-{full-summary.json,adlittle-300.log,bore3d-500.log,degen2-640.log}`.
- DEGEN3 dual: tau row-value lookup and immutable BTRAN metadata are settled. Reprofile before inspecting remaining iteration control, updates and scratch allocation.
- PILOT4 primal-release: the pre-packed-predicate accepted-release self profile puts primal control (13.29%), ordinary row update (11.06%), BTRAN (10.27%) and factor update (8.97%) first. The cached original-reference comparison identifies pricing/control and factor construction/update as the largest excesses; its historical shares are not current solver costs. BTRAN has a much smaller instruction excess. Direct BTRAN-list traversal to avoid the row-product order scan/sort is rejected below after a quad primal regression. The remaining hot primal samples lie in row-product scatter and the Devex recurrence/norm test. Short Devex dots already use scalar length-unrolled code; the long eight-lane path is cold in the current annotation. Non-precise branch samples cannot assign misses to the norm check, and branch-stack sampling is unavailable. The paired short-dot, static row-policy and finite-guarded two-output weight trials are rejected below. The current BTRAN annotation samples its dense L-transpose scatter strongly; candidate-list membership has no samples in this recording, while scoring dominates the list-update samples. Those are self-sampling observations, not zero-cost proofs or removable bounds. A/B lookup reuse and batching U-transpose deletions are rejected below; their tiny instruction changes do not justify a larger index solely from redundant-looking source scans. Fresh-factor temporary adjacency, per-tail-column traversal scratch and invariant structural entry counts are rejected below. Nucleus annotation samples immutable coefficient searches and dense-handoff row hashing; the shared coordinate-workspace variant is timing-inconclusive below. Produced factor data and basis-column copies remain constructor leads.
- Factor dependency walk: current accepted-release annotation samples repeated A/B length counts and the slot-maximum reduction (`target/perf-update-current/update.asm`). Count caching and maximum algebra have no established runtime gain; reusing the existing reverse row copy saves instructions but misses independent timing confirmation. All three implementations are restored and recorded below. Reprofile before another local variant rather than treating these medians as accepted improvements.
- Ranked double products: both per-term filtering and maintained required-entry traversal are rejected below. The latter preserves coefficient/duplicate order across ordinary pivots and arbitrary basis differences, but its 32-pair screen has no established timing benefit despite fewer instructions. Do not assume that removing these unused outputs accounts for a comparable runtime share. Reprofile remaining primal control or factor-update costs before another representation variant; the independent cost-change caller still uses the full routine.
- Other reachable runs follow their table ratios. Four compact extended factors and double structures remain cached; providers reuse exact RHS solve/estimate and pure fresh reduced-cost results. Main factors still clone mutable state; row-capacity growth already uses the numeric U transpose.
- Quad_2717: the accepted `vypkuxry` refresh follows pure increment reuse and still puts extended BTRAN first. Its fixed numeric/marker scratch can use initialized stack prefixes with exact outputs/WORK, but the tested implementation increases cycles and remains unapplied. Division-shape counts below distinguish zero numerators from their call-site codegen; dense U already inlines its zero return before quotient helpers. Native dense/packed BTRAN and the compact BTRAN-result cache remain unqualified. Ranked-product register retention/outlining and global power-of-two multiplication remain rejected; use current attribution before another arithmetic/storage variant.
- Barrier: reprofile before choosing more matrix-product, assembly or triangular work. Eight-row/four-column broadcast sharing and four-row/eight-column load sharing have no established timing gains below; reprofile before another tile width or copy-mode change; the four-output copy Chain trial is timing-inconclusive below. Preserve every per-output reduction and block boundary. Four-row/four-column SIMD S tiles with scalar two-column tails and allocation-free row-tile selection are accepted above. The two-column-only vector trial has no established timing gain below. Local scatter row/slice hoisting is timing-inconclusive below; dense-LDL chunks have no established timing gain. Refresh attribution before another local variant. DTRSM dispatch is settled; reprofile remaining FMA callers before choosing another whole kernel.
- Allocation: the cached native caller sample attributes about 58% of allocation frequency to refactorization, 19% to iteration outside refactorization. Flat auxiliary LU storage and dual factor-input reuse reduce allocation calls without a timing gain and are rejected below. Current release instruction profiles put FTRAN/L and factor/update traversals ahead of more small buffer changes; allocator instruction ancestry is now grouped in the frame-pointer screen above, with immediate-site and release-cost limits. Native sampling dispatch and bindings are verified; the rejected CLI-local small-block cache remains below.
- Dense BTRAN transpose geometry: forming the immutable transpose in constructor order from cached ADLITTLE/BORE3D/DEGEN2 BSOLVE factors gives 376/763, 972/2663 and 4206/14468 entries in consecutive destination pairs; four-wide groups cover only 100/763, 244/2663 and 328/14468. These are structural counts across captured factors, not dense-branch frequency or time savings. The structural query keeps rows in their original order. Stable descending-column gather equivalence is tested separately below and is slower in the standalone probe; no gather solver integration is retained. Data/source locations: `target/perf-primal-0436/transpose-geometry.json`. Existing transpose layout, scalar unroll and bounds-check rejections still apply.
- Dual tau staging: tau_dense separately serves only the pivot diagnostic and gathers at the column output positions. Current accepted PEROLD annotation has two scatter-loop instruction samples out of4303 total; allocator cost is separately unassigned. Avoiding this expansion remains untested, but the visible local work does not justify prioritizing it over pivot-row scatter. Annotation: `target/perf-ordinary-260d/iteration.asm`; accepted-source refresh: `target/perf-ordinary-devex/iteration.asm`. Fixed first-touch storage in the larger scatter is timing-inconclusive below; do not extend its screen.
- Dumps disabled: source audit of the simplex, dual, factor and PDHG hot paths finds JSON emission/serialization guarded by dump selection or an existing recording option; WORK estimates outside those guards remain solver inputs. The cached accepted-release PEROLD self report has no samples under the selected serde/dump/OnceLock symbol names, but inlined checks and preparation remain charged to callers. This is not an inclusive/total overhead bound or proof that all preparation is free. No dump optimization is justified by this screen; source anchors and the selected-symbol query are in `target/perf-dump-audit/`.

## Tried, no gain

- Memoize pure F80 quotients in a const thread-local table of512 entries (128 four-way FIFO sets,32KiB table). Eligibility is finite nonzero inputs; raw full operand fields and exact-key verification protect fingerprint collisions. The complete original division remains independently frozen, and source reversal removes only lookup/store around its unchanged arithmetic. Original zero/NaN/infinity paths remain numerically unchanged; all stores use original packed results. Collision/eviction and cold/warm random/full-bit/boundary/captured comparisons check internal sign/exponent/mantissa. Matching fixed39 PGO, all eleven F80 tests and fast940 parity pass; eight ADLITTLE/fx and three Quad complete result/exit pairs agree. Despite high simulated reuse, fx instructions fall only0.473% with cycles+0.845%; Quad instructions+0.230%, cycles+5.671%. Stop at counters without seeds21217/21289 or full gates. Assembly/relocation audit confirms both baseline and candidate retain only one out-of-line Div call site in BTRAN; the original dense-U division remains inlined, so loss of that inlining does not explain the regression. Source, accepted CLI and metadata restored exactly; restored fast940 parity passes. Full source/frozen oracle, proof, matching build, tests, assembly, call audit and all paired result/counters: `target/perf-div-memo/`. This rejects the tested global table; the separately tested column scope also loses at counters, as recorded below.

- Borrow magnitude-ordered F80 operands instead of copying their fields. Full original add_signed is independently frozen; reversal reproduces it exactly after the selector substitution. All arithmetic, shift/round/pack expressions, exceptional cases and WORK remain unchanged. New complete-routine differential tests compare internal sign/exponent/mantissa on random/full representations, exponent-gap/zero/subnormal/overflow/nonfinite boundaries, both effective signs and captured SC50A values. Matching fixed39 PGO, all ten F80 tests and fast940 parity pass; eight ADLITTLE/fx and three Quad full result/exit pairs agree per variant. V1 instructions fall0.264% on Quad, but cycles rise5.31/5.87/5.69% and branch misses rise21.9–29.4%; native ordering introduces branches and indirect sign loads. V2 uses std::hint::select_unpredictable for the same two references (availability compiled on pinned1.98.1; local core implementation inspected). It removes most of that miss penalty, but Quad instructions rise0.247% and median cycles1.269%; fx instructions rise0.058% with cycles-0.306%. Neither variant has useful target gain. Stop at counters without either reserved timing pair or full gates. Source, accepted CLI and build metadata restored exactly; restored fast940 parity passes. Original/frozen/candidate source, exact reversal and mode-independent value proof, both matching builds, tests, native assembly and complete paired outputs/counters: `target/perf-f80-add-borrow/`, v2 under `v2/`.

- Use initialized stack prefixes for five fixed-size extended BTRAN numeric/marker buffers on factors of at most32 rows; other cases retain heap storage. Dynamic support/result lists stay unchanged. The partial L handoff becomes an exact slice Copy. Source reversal reconstructs the complete original BTRAN; its independent body remains byte-identical (only signature name/visibility/format differ). Boundary, untouched-suffix, F80 low-bit partial-path and complete f64/F80 output/WORK/cache/density/initialization comparisons pass across sizes1/4/13/31/32/33/65, duplicate/extreme/two-NaN-payload inputs and sparse/dense/direct/eta cases. Matching fixed39 PGO, all nine FTRAN/L tests and fast940 parity pass; eight ADLITTLE/fx and three Quad complete result/exit pairs agree. Instructions change -0.004/-0.142/-0.396%; fx/Quad median cycles rise1.46/1.50% and solve medians rise1.91/2.12%. Stop at counters without seeds18217/18289 or full gates. Source, accepted CLI and build metadata restored exactly; restored fast940 parity passes. Full source/oracle, initialized-prefix safety/reversal proof, test/build provenance, native assembly and every paired result/counter: `target/perf-btran-stack/`.

- Cache exact immutable L-column WORK charges at fresh construction. Counts up to 2^53 retain the original `len as f64 * 2.0` value in every MXCSR mode; larger counts use the original expression. All additions, traversal and numeric scatter stay in place. Independent original L-stage and dense-WORK routines verify complete f64/F80 values/counts/WORK on generated/extreme/captured inputs, and a synthetic-offset test covers the exact-conversion boundary without huge allocations. Matching fixed39 PGO, all seven FTRAN/L tests and fast parity940 pass. Eight ADLITTLE/DEGEN2 and three PEROLD complete result/exit pairs agree. Instructions change -0.001/-0.096/-0.073%; PEROLD paired cycles -3.392/+0.334/+0.365% and solve -10.954/+0.732/+0.673%. The first pair crosses a host speed band; the remaining pairs show no useful gain. Stop at counters without seeds16217/16289 or full gates. Source, accepted CLI and metadata restored exactly. Complete source/oracles, bounds/ownership proof, matching build, tests, assembly and counters: `target/perf-lwork-charges/`.

- Extract the complete dense L numeric/WORK traversal into a non-inlined helper, used at the original L-stage and dense-WORK call sites. No factor fields, arithmetic expressions, ordered additions, tolerances or scatter order change; reversal reproduces both complete original methods and the helper loop exactly. Independent full L-stage and dense-WORK comparisons pass for f64/F80 generated/extreme/captured inputs. Matching fixed39 PGO, all six FTRAN/L tests and fast parity940 pass; eight ADLITTLE/DEGEN2 and three PEROLD complete result/exit pairs agree. Native ordinary assembly retains WORK in XMM0 instead of a stack slot. Instructions fall0.097/0.410/0.190%; PEROLD paired cycles -0.033/+2.951/-0.224% and solve -0.928/+8.529/+0.548%, with no useful target gain. Stop at counters without seeds17137/17209 or full gates. Source, accepted CLI and metadata restored exactly. Full source/reversal proof, original oracles, build/test provenance, assembly and all counter/result pairs: `target/perf-dense-l-outline/`.

- Select child-update rounding mode once per call, specializing both complete width-one and blocked child kernels without changing arithmetic expressions, update order or DGEMM/block boundaries. Complete original child bodies are frozen independently; reconstruction matches after method/rounding-selector substitution only. Original-versus-new tests cover double/single modes, 320 synthetic shape/value families, signed zeros/subnormals/extremes/two NaN payloads/captured kernel values and real ADLITTLE/SC50A/BORE3D front updates. Matching fixed39 PGO, both tests and fast parity 940/940 pass. Eight BORE3D/DEGEN2 and three DEGEN3 barrier complete result/exit pairs agree. Instructions rise +0.57/+1.01/+1.28%; cycles +1.27/+1.56/+2.14%. Stop at counters without seeds7151/7213 or full gates. Source/qualified CLI/build restored exactly. Frozen bodies, complete candidate/build, tests and measurements: `target/perf-child-rounding/`.

- Store scalar permutation slots instead of (slot, position) tuples in immutable-slot U-order and row-elimination heaps, decoding positions with at[slot]. Slot/at bijection preserves priority and duplicate multiplicity. Independently frozen original uorder/eliminate bodies compare sparse/dense inputs, random permutations, duplicate initial entries, f64/F80 values/order and exact visited/dense WORK bits within the existing complete-FTRAN test. Source reconstruction confirms only heap types/keys/decode change. Matching fixed39 PGO, that test and fast parity 940/940 pass; eight ADLITTLE/DEGEN2 and three PEROLD/PILOT4 complete result/exit pairs agree. Instructions fall 0.029/0.113/0.249/0.144%; cycles +0.246/+0.392/+0.077/-0.337%. The predeclared PEROLD screen (seed7001, CPU4, cap1300, 32AB+32AA) is inconclusive: median -0.124%, interval [-0.414%, +0.164%], AA [-0.457%, +0.243%], 19/32 improving, sign p0.189, both orders improve. Stop without seed7069, extension or full gates. Qualified source/CLI/build restored exactly. Candidate, source/frozen-body proof, tests and all raw measurements: `target/perf-slot-heaps/`; preceding heap source annotation `target/perf-heap-audit/`.

- Extract the dense nucleus handoff into a non-inlined helper without changing any arithmetic, WORK, buffer contents/order or failure branch. Reconstructing the original kernel and handoff is source-identical after whitespace and the boolean return adaptation; dense-LU kernels/tests are unchanged. Matching fixed39 PGO, frozen f64/F80 dense-LU tests and fast parity 940/940 pass. Ordinary nucleus code shrinks 31,086 -> 26,157 bytes and local stack 1,736 -> 1,496 bytes; dense LU inlines into the new helper. Eight ADLITTLE/DEGEN2 and three PEROLD/PILOT4 complete result/exit pairs agree. Instructions change -0.0065/-0.018/-0.023/-0.019%; cycles +0.106/+0.154/+0.041/+0.400%. No useful target counter gain; stopped before seeds6883/6947 or full gates. Source/qualified CLI/build restored exactly. Candidate, normalized source proof, exact disassembly/size/frame, tests and all results/counters: `target/perf-nucleus-outline/`.

- Construct nucleus nodes directly from active incoming columns, retaining forward conversion order and reversing each completed node segment to reproduce exact original node IDs and column chains. Removes temporary filtered numeric column vectors; numeric kernel/WORK source is identical. Independent frozen input construction compares complete nodes, chains, buckets and keys on 1,024 families in both f64/F80, including duplicates, inactive rows/columns, threshold edges and nonfinites. Matching fixed39 PGO, that test and fast parity 940/940 pass. Eight ADLITTLE/DEGEN2 and three PEROLD/PILOT4 complete result/exit pairs agree. Instructions fall 0.17/0.12/0.16/0.21%; cycles rise 1.30/0.17/1.21/0.38%. No useful target counter gain; stopped before seeds6761/6827 or full gates. Source/qualified CLI/build restored exactly. Candidate, original producer/initializer, source-identical kernel proof, tests, build and all results/counters: `target/perf-nucleus-input/`.

- Lazily retain exact rounded ranked-product vectors alongside existing compact extended cache entries. Original numeric kernels/conversion and WORK sites are unchanged; tests compare original trait conversion bits for random/special inputs, lazy fill, clones and fallback. Diagnostic Quad_2717 calls reuse 461,310 rounded values versus 1,080 computed; fx reuses 171,072 versus 693. Matching fixed39 PGO, all six cache tests and fast parity 940/940 pass; eight ADLITTLE/fx and three quad complete result/exit pairs agree. Instructions fall 1.37% on quad and 1.22% on fx. The predeclared quad screen (seed6563, 32AB+32AA, CPU4, cap17000) has median -1.27%, interval [-2.79%, -0.235%], AA [-0.718%, +0.439%], 26/32 improving, p0.000268, both orders improve. Its conservative bound does not clear AA uncertainty: inconclusive. Stop without seed6629, extension or full gates. Source/qualified CLI/build restored exactly; candidate, frozen source/kernel proof, profile, diagnostic and all raw samples remain in `target/perf-rank-rounded/`.

- Sparse-L level queues: a per-node linked FIFO replaces per-level vector growth, then a flat FIFO uses cached exact level capacities and phase-local reuse of slot-sorting scratch. The original frozen L body remains source-identical. Both preserve level/first-insertion order, including duplicate input/destination entries, same-level appends, backward-level assignments and threshold/nonfinite behavior. Frozen f64/F80 tests exercise captured SC50A data through levels; final tests also prove original node-prefix/exact capacity metadata and alternate queue use with slot sorting. Both pass fast parity 940/940 and eight ADLITTLE/DEGEN2 plus three PEROLD/PILOT4 complete result/exit pairs. Linked queues with the preceding PGO profile add PEROLD instructions 0.40%, cycles are essentially flat; stopped before own PGO/timing. Flat queues with that preceding profile add ADLITTLE/DEGEN2/PEROLD/PILOT4 instructions 2.11/1.70/1.21/1.16%. Matching fixed39-case PGO changes those to +0.008/-2.00/-1.15/-0.36%, with PEROLD cycles -1.07%; preliminary old-profile counts do not represent final codegen. The frozen own-PGO PEROLD screen (32 AB + 32 AA, CPU4, cap1300, seed6323) is inconclusive: median -0.4954%, interval [-0.78%, -0.15%], AA [-0.36%, +0.41%], 24/32 improving, sign p .00350 and both orders improving. It fails the declared minimum gain and noise separation; stopped without confirmation6389, sample extension or full gates. Qualified source/CLI/build restored exactly. Linked implementation/source/binary/tests/counters: `target/perf-level-links/linked/`; flat final source/tests, source/ownership proofs, both profile builds/counters, complete controlled schedule/samples and restoration: `target/perf-level-flat/`. No qualifying queue speedup is credited.

- Direct first-match sequential ordinary tau dots when both virtual and numeric U rows have fewer than eight entries, skipping row-map preparation and temporary term construction. The original generic and extended row-dot kernels and frozen complete needed-FTRAN body remain source-identical; WORK sites are untouched. Test-only row counts on accepted `kpttrrms` show 418/418 ADLITTLE, 16,714/19,765 DEGEN2 and 33,243/36,518 PEROLD nucleus dots have short virtual rows. These counts do not measure eligibility of the additional numeric-row guard or removable cost. The final PGO build passes 100,000 bitwise original-pipeline cases with independently generated coefficients (including differing duplicate values, missing entries, signed zeros, full random doubles, extremes and NaN payloads), derived captured SC50A rows, existing frozen f64/F80 FTRAN tests and fast parity 940/940. Eight ADLITTLE/DEGEN2 and three PEROLD/PILOT4 complete result/exit pairs agree. Final-PGO instructions change -0.43/-0.93/-0.56/+0.002%; initial cycle changes -1.19/+11.26/+0.94/-1.60% include host-band shifts and are not credited as gains. The fixed PEROLD screen (32 AB + 32 AA, CPU4, cap1300, seed6007) is inconclusive: median +0.06%, interval [-2.46%, +1.51%], AA [-1.79%, +0.87%], 16/32 improving pairs, sign p .570, one execution order slower. Stopped without reserved confirmation6073, sample extension or full gates. Stronger duplicate-value tests are frozen before official timing; final PGO is regenerated from that source. Source, qualified CLI and PGO metadata restored exactly. Full source/binary/build versions, diagnostic distribution, unchanged-oracle proof, all tests/results/counters, source freeze and complete failed screen/schedule: `target/perf-short-tau/`. No qualifying speedup; the common short-row shape does not alone justify replacing the current lookup path.

- Replace F2's numeric column snapshot with the existing output bitmap's nonzero predicate, retaining the tau snapshot. Original full/needed numeric snapshots and joint-U/third/dense WORK bodies are independently frozen; five bodies match after identifier/visibility/whitespace normalization only. Tests compare repeated f64/F80 outputs, stored predicates and exact WORK on exceptional/random/captured SC50A cases plus 3,072 randomized original accounting states. V1 revisits zero output positions; the final three-state marker preserves zero visitation. Both pass fast parity 940/940 and eight ADLITTLE/DEGEN2 plus three PEROLD/PILOT4 complete result/exit pairs. Preceding-profile final instructions rise +0.20/+1.02/+1.02/+0.93%. The unchanged final source/tests are then checked once with their first matching fixed39-case PGO profile, prompted by the level-queue experiment's measured profile effect: instructions become -0.12/+0.10/+0.07/+0.076%; cycles +0.23/-0.63/-0.51/-0.43%. Matching-profile original-snapshot/WORK tests and fast parity pass. PEROLD fails the predeclared instruction-improvement gate, so neither controlled seed6449/6503 nor full gates runs. The predicate equivalence is verified, but no storage implementation qualifies. Source/qualified CLI/build restored exactly. Parent-profile source/binary/oracles/tests/counters: `target/perf-f2-colflags/`, V1 under `v1/`; matching-profile source freeze/build/tests/counters and stop/restoration proof: `target/perf-f2-colflags-pgo/`.

- Reuse existing F2 column/input/tau/raw snapshot allocations as the next solve's zeroed workspace, then move the exact new snapshots back before return. No factor fields, numeric operations or WORK sites change. Existing frozen full/needed FTRAN bodies remain source-identical; the release differential test is expanded across repeated full/needed histories, all snapshot values/flags, f64/F80 extremes and captured SC50A inputs. That test and fast parity 940/940 pass. Eight ADLITTLE/DEGEN2 and three PEROLD/PILOT4 complete result/exit pairs agree. Under the accepted parent PGO profile, instructions rise +0.17/+0.52/+1.09/+0.94%; cycle medians rise +1.79/+1.08/+1.11/+0.77%. The complete implementation adds measured cost despite eliminating allocations/copies; rejected before new PGO training, controlled screens or full gates. Accepted source, CLI and build are restored exactly. This tests persistent workspace plus delayed snapshot installation, beyond the earlier move-only variant; it does not establish that all data in the numeric snapshots is necessary. Frozen-oracle audit, source/binary/build context, tests, all bounded results/counters, accepted-profile refresh and restoration proof: `target/perf-f2-workspace/`.

- Larger compact RHS/WORK/reduced-cost or fresh dense-BTRAN result-cache bounds: test-only exact input traces on accepted `kpttrrms` measure FIFO hit counts independently per provider generation, excluding pointer reuse. AFIRO/ADLITTLE have no eligible calls. Fx has 1,740 RHS/WORK/fresh-dense-BTRAN calls with 1,729 four-entry hits, and 2,612 reduced-cost calls with 2,601 hits. Quad_2717 has 2,894 RHS/WORK/fresh-dense-BTRAN calls with 2,884 hits, and 5,786 reduced-cost calls with 5,776 hits. Primal_5959 primal has 5,331 RHS/WORK/fresh-dense-BTRAN calls with 5,278 hits and 5,345 reduced-cost calls with 5,294 hits; dual has 7,792/7,753 and 7,802/7,763, respectively. Every count is identical at capacities 4/8/16/32/64. Runs are bounded and retain expected statuses/iterations (fx infeasible/5319, quad optimal/16010, primal_5959 numerical_failure/10780 and /15688). These are diagnostic traces/simulations, not solver timing or evidence of profitability of the earlier dense-BTRAN cache. No resizing or optimization integration; diagnostic code is removed, and accepted source/CLI/build remain exact. Source, test log, per-owner/global distinct counts, all simulations and restoration proof: `target/perf-cache-working-set/`.

- Four-entry compact-owner extended ranked-product memo: exact sorted participating-row/value keys preserve the unchanged numerical kernel, but the bound thrashes the observed cycles. Release differential and fast parity 940/940 pass, and eight ADLITTLE/fx plus three Quad_2717 complete result/exit pairs agree. With the accepted parent PGO profile, instructions change +0.025/+0.558/+1.444% and cycle medians +0.045/+3.156/+5.365%, respectively. Stopped before own-PGO training or controlled screens. A test-only exact-key trace then measures 22/23 distinct inputs over 25,046/6,073 calls; this supports the separately qualified 32-entry implementation above. V1 source/binary/build context, tests, all results/counters and diagnostic trace: `target/perf-rank-memo/`; V2 qualification is isolated under `v2/`.

- Prepack four A rows into contiguous four-lane records once per sequential-S DGEMM grid, replacing repeated strided scalar loads/shuffles with vector loads while retaining every existing FMA and C-update operation. The unchanged current four-column kernel is the independent oracle. A standalone probe passes 512 exact comparisons on captured AFIRO values, random finite representations, signed zeros, extremes and boundary K sizes. Both paths include output allocation; the candidate additionally includes allocation and construction of packed A on every call. After an initial 12M-term calibration (0.7–3.5ms samples), a fixed 100x workload with six alternating pairs per shape gives median changes 8x8x7 +4.97%, 16x16x32 +7.02%, 24x48x128 +2.17%, 48x48x192 +3.35%, 64x96x192 -3.37%, 128x128x192 +14.86%. The improving 64x96 grid routes through gemm_path2 in production, so it is not eligible for this nocopy-only mechanism. Synthetic grids use sequential-S tiles throughout and do not model full mixed-mode DGEMM; these timings are feasibility evidence, not solver gains. No integration or full solver gates. Solver source, qualified CLI and build remain unchanged. Frozen/generated kernels, captured-input provenance, both workload versions, raw timings and summary: `target/perf-gemm-packed-a/`.

- Restrict fresh dense extended BTRAN result reuse to the compact RHS owner and the dual-recomputation consumer. The independent main factor runs the unchanged numerical BTRAN on misses; hits restore exact WORK, density and inverse-row initialization, retaining update-saving row metadata. Four-entry exact ordered input keys exclude NaNs, dumps, updates and row etas. A release test covers 64 input families with signed zeros, low bits, extremes/infinity, duplicates and clone histories, followed by FTRAN/update/BTRAN continuation; fast parity passes 940/940. Using the accepted parent PGO profile for this initial screen, eight complete ADLITTLE/fx result/exit pairs and three primal_5959 primal pairs agree. Instructions change +0.14/-0.77/-1.66%; cycle medians +0.76/+1.63/-0.11%. No useful initial target cycle gain, with fx cycles higher; rejected before new PGO training, controlled screens or full gates. Source, accepted CLI and build metadata restored. This is narrower ownership than the earlier general main-factor cache and adds no RefStruct fields; neither implementation establishes a solver timing gain. Source snapshots, original-kernel audit, passing test/parity logs, all bounded outputs/counters and restoration proof: `target/perf-owned-btran/`.

- Pack software F80's sign and unchanged normalized exponent into one metadata word beside the mantissa, retaining the integer arithmetic and rounding expressions. Production disassembly shows the original operators use pointer arguments and memory returns; the standalone packed prototype uses the two-word ABI. Numeric differential checks pass for 100,157 pairs (edge cases including negated NaNs, captured QTAU values and random representations), plus 100,000 from-f64 and 100,000 from-i64 conversions. Checks include representation, sign, classes, f64 conversion bits, arithmetic, square root, comparisons, min/max and abs/neg; Debug formatting and full solver parity are not tested. Four fixed alternating CPU-4 pairs of 640,000 recurrence terms are all slower: median +19.82% with inlined operations and +8.17% with explicitly outlined add/multiply calls. These are standalone non-PGO feasibility timings, not controlled solver screens. Rejected before integration; production source and accepted CLI/PGO metadata remain unchanged. Transform, original/generated source, validation and raw timing logs, compiler metadata, production/probe disassembly and summary: `target/perf-f80-layout/`. This rejects normalized metadata packing for the tested software arithmetic; it does not establish a bound for other representations or whole kernels.

- Four-row/eight-column sequential DGEMM S tiles share the four-row A pack across eight independent column chains; every existing per-output FMA term, initial zero, K-block boundary, mode and C update remains unchanged. Overlapping output columns retain the existing dispatch. The current four-column kernel is independently source-verified frozen. All eight final-PGO Cholesky/BLAS release tests pass, including direct wide-versus-frozen tile comparisons (random/full finite representations, signed zeros, extremes, captured AFIRO values, boundary K sizes) and expanded whole-driver shapes; fast parity passes 940/940. Eight BORE3D and three DEGEN3 complete result/exit pairs agree. Instructions BORE3D -0.20%, DEGEN3 -0.74%; target cycles -0.91/-1.13/-1.23%. First predeclared screen (32 AB + 32 AA, CPU 4, cap 100, seed 4871) passes: median -1.58%, interval [-1.91%, -1.31%], AA [-0.63%, +0.21%], 29/32 improving, p 1.28e-6, both orders better. Independent confirmation 4931 is inconclusive: median -1.23%, interval [-2.42%, +0.67%], AA [-0.36%, +0.56%], 22/32 improving, p .0251; both orders better but significance/noise separation fail. Stopped without sample extension/repeated confirmation or full gates; no qualifying gain is credited. Source, accepted CLI and PGO metadata restored. Frozen/candidate source, binaries/builds, all tests and bounded counter/result pairs, complete controlled schedules/raw outputs and restoration proof: `target/perf-gemm4x8/`. This tests wider-column load sharing with the existing per-step row pack; it does not test persistent packed A storage.

- Dense BTRAN SSE2 scatter feasibility: pair adjacent destinations at runtime, precompute consecutive suffixes, or pair arbitrary distinct destinations with sequential duplicate fallback. All standalone variants retain per-destination multiply/subtract order and WORK bits on captured ADLITTLE/BORE3D/DEGEN2/PILOT4 factors, random representations/nonfinites/extremes and duplicate rows: 1,156 checks for runtime adjacency, 1,284 for suffix plans (including additional dense NaN-coefficient rows), 1,156 for arbitrary pairs. Eight alternating kernel pairs per factor retain identical vector resets; suffix setup is outside timing. Runtime adjacency medians are +13.44/+25.48/+37.45/+12.84%; arbitrary pairs +5.35/+13.19/+11.77/+13.05%. The suffix plan is eligible only on the largest PILOT4 factor and gives -1.14% within wide pair variation, establishing no gain; other factors use the scalar fallback. A structural four-source/two-target intersection screen covers 0/0/24/440 of 106/464/1578/2016 entries, not runtime savings; no blocked kernel was built. None is integrated into the solver, and no solver-speedup claim is made. Source/build, exact checks, all timing samples and structural coverage: `target/perf-bt-simd/`, subdirectories `planned/` and `arbitrary/`; factor provenance remains `target/perf-bt-gather/`.

- Replace immutable per-pivot nucleus L/U coefficient linear searches with first-occurrence position arrays, clearing only touched coordinates after each pivot. Final variant reuses the cleared row workspace for terminal dense handoff, retaining the previous HashMap last-write row-index rule. Arithmetic, node visitation, matrix stores and WORK sites/order are unchanged; no persisted factor fields are added. Both variants pass fast parity 940/940, four LU release tests including frozen f64/F80 dense kernels, eight ADLITTLE/BORE3D canonical result/exit pairs and three PILOT4 pairs. Sparse-only variant instructions PILOT4 -0.99%, DEGEN2 primal -0.025%, PEROLD dual -0.61%, with no initial cycle gain; the complete shared-workspace variant saves PILOT4 instructions 1.46% and initial median cycles about 1.2%. Small final instructions are +0.09/-0.005%. The predeclared final 32-AB/32-AA screen (CPU 4, cap 800, seed 4663) gives median -0.63%, interval [-0.95%, -0.13%], AA [-0.46%, +0.42%], 23/32 improving pairs (sign p .0100308), both orders improving. It fails sign significance and noise separation; stopped without reserved confirmation/full gates. No qualifying timing gain is credited; source, accepted CLI and PGO metadata restored. Ownership/lookup semantics, variant source/build/binary snapshots, tests, all bounded result/counter pairs and complete controlled schedule/raw samples: `target/perf-nucleus-index/`, final variant under `final/`.

- Dense BTRAN L-transpose gather feasibility: stable descending-source columns accumulate each output in a register, preserving the original per-output term order, duplicate order, tolerance filtering and descending WORK additions. A standalone release probe compares the unchanged scatter with the gather on cached ADLITTLE/BORE3D/DEGEN2 factors and a bounded current PILOT4 factor capture (cap 800, optimal 690 iterations). All 1,156 checks preserve output and WORK bits, including real factor/vector data, random full double representations, signed zeros/nonfinites/extremes and duplicate/unsorted strictly-lower columns. Metadata construction is outside timing; both timed kernels include the same vector reset. Eight alternating within-process pairs per largest captured factor show gather median +20.39/+17.07/+39.27/+52.05%, all pairs slower. Stopped before solver integration. This proves a tested dataflow equivalence, not solver/reference parity or a universal lower bound on gather implementations. The solver/qualified CLI are unchanged. Source/build, factor provenance, cached/new inputs, all validation and timing samples: `target/perf-bt-gather/`.

- Batch two independent ordinary Devex weight recurrences with SSE2, preserving each column dot and every multiply/subtract/add boundary. The original scalar row body and complete frozen oracle remain source-identical. The first helper-level nonfinite test exposes a NaN payload difference in an inlined scalar fallback (round 9); the final entry admits finite immutable coefficients and finite state/vector/scalar inputs, otherwise executes the original complete row kernel. Immutable coefficient eligibility is computed once; per-call state/vector checks remain. Two final-PGO release tests pass: 50,000 finite recurrence pairs (including finite inputs producing overflow/cancellation), plus captured scaled-AFIRO and expanded random whole-row cases covering modes, original exceptional fallback, signed zeros/NaN payloads, empty dimensions and extended state. Fast parity passes 940/940. Eight ADLITTLE/BORE3D and three PILOT4 complete result/exit pairs agree. Small instructions rise 0.07/0.13%; PILOT4 instructions rise 1.19% and initial median cycles about 3.19%. The full machinery has no benefit; stopped before controlled screens/full gates and restored source, accepted CLI and PGO metadata. Disassembly verifies packed multiply/divide/subtract/add and no new FMA. First/final source, failing/passing tests, exact binaries/builds, all counter/result pairs, oracle proof, assembly and restoration proof: `target/perf-primal-weight2/`. This rejects the complete finite-guarded two-output weight implementation; it does not bound every possible batched kernel.

- Guard the ordinary primal row entry once for present mirror/products, mode 1, faithful weights and no excluded variable. A const-specialized loop retains the original arithmetic; the complete independent frozen row oracle remains byte-identical. Parent/final-PGO release differential tests explicitly exercise the guarded entry on captured scaled AFIRO data and random/empty/extreme/signed-zero/NaN/duplicate cases. Fast parity passes 940/940; eight ADLITTLE/BORE3D canonical result/exit pairs and three PILOT4 pairs agree. Small instructions are flat; PILOT4 instructions fall 1.68%, initial median cycles about 0.97%, and the hot row body shrinks 3057 -> 2394 bytes. The predeclared 32-AB/32-AA screen (CPU 4, cap 800, seed 4421) gives median -0.83%, interval [-4.11%, +0.89%], AA [-2.60%, +0.62%], 19/32 improving pairs, sign p .189; both orders improve, but noise separation fails. Stopped without reserved confirmation/full gates; source, accepted CLI and PGO metadata restored. No qualifying gain is credited. Exact source/binary/build, oracle proof, tests, small/direct counters/results, symbol sizes and complete controlled schedule/raw outputs: `target/perf-primal-prepared/`.

- Two-word F80 software layout with unchanged m/exponent/sign semantics and original integer arithmetic. Complete original type, reversible source edits and logical sign-word proof independently protect the representation. Both matching fixed39 PGO variants pass all10 F80 release tests (new full-field oracle plus existing conversion/x87 tests), all6 FTRAN/L tests and fast940; eight bounded ADLITTLE/fx plus three Quad complete result/exit pairs agree. Initial mechanical layout: fx instructions-0.013%, cycles-3.182%; Quad instructions+0.043%, cycles+0.220%, branch misses+16.676%. Direct sign-bit negation and order-independent sign comparison retain the same numerical branches while avoiding metadata reconstruction: fx instructions-0.526%, cycles-2.508%; Quad instructions-0.440%, cycles-1.759%, branch misses+9.920%. The predeclared32-AB/32-AA Quad screen (CPU4, cap17000, seed28217) is inconclusive despite27/32 improving pairs: median-2.203%, interval[-4.210%,-1.191%], sign p5.65e-5, both orders-2.881%/-1.709%, healthy AA[-2.943%,+2.425%]. Conservative noise separation fails. Stop without confirmation28289, extension or full candidate gates; no solver gain is accepted. Accepted source/CLI/metadata/profile restored exactly, temporary oracle file removed; restored formatting, fast940 and full987 pass. Full sources/original oracle/proof, builds/tests/freezes, all bounded outputs/counters and complete screen: `target/perf-f80-two-word/solver/` and `solver-v2/`. This establishes exact software-layout feasibility, not a qualified runtime benefit.

- Packed ordinary immutable L-entry indices: fitting ordinary factors store checked32-bit row indices with unchanged numeric values in12-byte packed entries instead of16-byte usize/f64 tuples; extended and out-of-range inputs retain wide storage. Logical entry counts, order, snapshots/dumps and WORK are unchanged. Column views expose exact Copy values; source reversal preserves fresh arithmetic and complete L-stage/FTRAN/needed/Kahn bodies apart from field installation and iteration binding. Existing complete L tests now read independent wide storage for the frozen oracle. All8 matching fixed39-PGO FTRAN/L tests pass, including exact payload/index encoding, clone/large-index/extended fallback, full/needed solves and captured SC50A L arrays. Fast940 and eight bounded ADLITTLE/DEGEN2 plus three PEROLD complete result/exit pairs agree. Instructions rise0.453/0.739/0.501%; PEROLD cycles+0.600%. Stop at counters without timing seeds27217/27289 or full candidate gates; no compression variant is accepted. Original source/accepted CLI/metadata/profile restored exactly; restored formatting, fast940 and full987 pass. Source/proof, codec and wide-storage oracles, matching builds/tests/freezes, all outputs/counters and native assembly: `target/perf-lpacked/`. This tests logical stream compression, not an allocator-footprint claim or direct numerical producer-U reuse.

- Ordinary positive weight-floor specialization: the floor is always finite and strictly positive (W_MIN or1), allowing register-only SSE2 maximum with original NaN selection; the extended path keeps its existing max. The independent complete row oracle is byte-identical; reversal changes only the one ordinary floor selection. A frozen runtime primitive passes2,097,696 full-bit/boundary/quiet-and-signaling-NaN comparisons across all16 MXCSR modes; the complete captured/random full-row test and fast940 pass with matching fixed39 PGO. Eight controlled ADLITTLE/BORE3D and three release PILOT4 complete result/exit pairs agree. Instructions fall0.143/0.050/0.381%; PILOT4 cycles+0.221% are effectively flat. Its predeclared32-AB/32-AA PILOT4 screen (CPU4, cap800, seed26217) is inconclusive: median-0.952%, interval[-3.249%,+2.688%],17/32 improving, sign p.430, healthy AA[-2.561%,+1.320%], execution orders-1.820%/+1.342%. Stop without confirmation26289 or full candidate gates. No qualifying runtime gain. Source/accepted CLI/metadata/profile restored exactly; restored formatting, fast940 and full987 pass. Native accepted-row annotation, complete source/oracle/proof, builds/tests/freezes, bounded results/counters and controlled screen: `target/perf-weight-floor/`; accepted annotation is also retained in `target/perf-rank-reuse/primal-row.asm`.

- Consecutive-destination SSE2 ordinary ranked-row products. Static fixed-pair coverage on scaled ADLITTLE/BORE3D/DEGEN2/PILOT4 is86/209,379/770,177/2122,1539/2647 pairs (no duplicate pairs in those models). Immutable model metadata admits only finite-coefficient rows with at least8 entries and two-thirds consecutive pairs; nonfinite input values retain the original scalar loop. Each packed lane updates a distinct output with the same multiply/add and per-output term order; gaps, duplicate destinations, rank ties and exceptional inputs retain their scalar behavior. The independently frozen complete original product, captured AFIRO rho fixture, random finite/extreme/duplicate/rank-tie inputs, existing nonfinite payloads and all16 MXCSR rounding/flush combinations pass; matching fixed39 PGO and fast940 pass for both variants. Eight bounded ADLITTLE/BORE3D pairs and three PILOT4 pairs preserve complete result/exit objects. V1 checks pairs and two-element slices dynamically: PILOT4 instructions+0.944%, cycles+0.495%; small instructions+0.057/+0.078%. V2 stores flat ordered consecutive-run ranges once per immutable model row, validates one destination slice per run and packs its independent pairs. Independent plan tests verify complete ordered partitioning, including duplicates and unsorted gaps. It still loses: PILOT4 instructions+0.814%, cycles+1.396%; small instructions+0.022/+0.119%. Native code confirms packed arithmetic, V1 pair/slice checks and V2 planned-run loops; no exclusive cause is assigned to whole-program counter changes. Neither enters official timing: seeds25217/25289 unused, no full candidate gates. Accepted source/CLI/metadata/profile restored exactly; formatting, fast940 and full987 pass. Static measurement, original/candidate source/proofs, tests/builds/freezes, complete result/counter pairs, native assembly and decisions: `target/perf-rank-pairs/` and `v2/`. This rejects the tested pair/run metadata traversal, not every independent-destination implementation.

- Exact F80 reciprocal parameters at BTRAN's three diagonal divisions. The complete original division and BTRAN are independently frozen; the public/default division stays unchanged. All three matching fixed39 PGO designs pass11 F80 tests (2,000,050 integer and1,037,249 complete-F80 oracle comparisons, boundary/nonfinite/captured inputs), all8 FTRAN/L tests, fast940 parity, and eight ADLITTLE/fx plus three Quad complete bounded result/exit pairs. A four-entry per-column FIFO shared across fresh-factor clones preserves live sign/exponent/numerator fields; Quad instructions+0.249%, cycles-0.537%; fx+0.059%/-0.747%. Earlier fallback checks, merged normalized-input checks and proved removal of unreachable rounding carry reduce overhead but still give Quad instructions+0.089%, cycles-2.688%; fx-0.035%/+2.349%. Neither lookup version enters timing qualification. A separate design prepares parameters at the only production diagonal write, copies them with factors and eliminates locking/lookup; clone/update/fallback tests verify exact parameter pairing. Quad instructions-0.280%, cycles-3.631%; fx instructions-0.216%, cycles-1.280%. Its predeclared32-AB/32-AA Quad screen (CPU4, cap17000, seed24217) is inconclusive: median-0.631%, interval[-4.598%,+2.406%],17/32 improving, sign p.430, healthy AA[-2.053%,+0.456%], execution-order medians+0.269%/-1.667%. Stop without confirmation24289 or full candidate gates. No accepted solver gain. Every numerical and WORK operation outside the three divisions remains original; source-reversal proofs cover complete original arithmetic/BTRAN. Accepted source/CLI/metadata/profile restored exactly; restored formatting, fast940 and full987 pass. Complete sources, tests/proofs, builds/freezes, counters/results, annotations, complete screen and decisions: `target/perf-diag-preinverse/`, lookup simplification `v2/`, prepared design `v3/`.

- Sixteen-entry pure per-column BTRAN quotient reuse, shared by clones of a compact fresh-factor owner and keyed by both raw operands across updates, preserves the complete original numerical body and every WORK site. Matching fixed39 PGO builds pass the complete original BTRAN differential oracle (double/extended, direct/dense/sparse, updates, extremes/NaNs/duplicates), collision/eviction/clone/key tests and fast940 parity. Eight ADLITTLE/fx and three Quad bounded pairs retain complete canonical results/exits. Linear fingerprint lookup: fx instructions+0.256%, cycles+6.482%; Quad instructions+0.710%, cycles+3.368%, branch misses+11.82%. A metadata-only integer SIMD mask variant passes its independent scalar-mask test and all10 FTRAN/L tests plus fast940; full keys still verify collisions. It reduces the branch penalty but still loses: fx instructions+0.387%, cycles+5.827%; Quad instructions+0.771%, cycles+3.000%, branch misses+0.976%. No official screen or full gates follow either variant; seeds22217/22289 and23217/23289 are unused. Original arithmetic/source, accepted CLI/metadata restored exactly; restored formatting and fast940 pass. Frozen/candidate sources, source proofs, matching builds, tests, disassembly, complete result/counter pairs and decisions: `target/perf-diag-div-memo/` and `v2/`. This settles the measured full-quotient-result table variants; divisor-parameter reuse has a smaller key/working set but its tested integrations also fail qualification above.

- Memoize complete BTRAN results across clones in the existing four-entry compact extended factor cache, only before any update. Exact input bits (including signed zeros and full F80 low bits), order/duplicates and both mode flags key a bounded 16-entry Arc/Mutex cache; NaNs/unsupported bit keys/oversized inputs fall back. Nonempty cached outputs carry exact WORK/density; row_cached follows the current call, and inverse-row initialization is preserved. The original numerical BTRAN body is source-verified unchanged. Release tests cover AFIRO-based mixed factors, modes, clone histories, bounds, updates, empty/extreme/duplicate inputs, exact output/WORK bits and Send compatibility; three tests in the final variant also force a fingerprint collision and verify full-key comparison. Both variants pass fast parity 940/940 and eight bounded ADLITTLE/fx_quad_obj_5371 canonical result/exit pairs per model. V1 allocates each query key: primal_5959 primal instructions -4.73%, direct cycles +2.41/-1.72/-0.36%; grouped +54% solve crosses host bands and is not reproduced by direct pairs. V2 uses bounded initialized stack key storage and a fingerprint plus exact comparison, allocating retained keys only on misses. ADLITTLE instructions +0.09%, fx -2.54%, target -4.80%; target direct cycles -1.74/-2.12/-1.14%, while small fx median cycles rise 10.4%. The predeclared target screen (32 AB + 32 AA pairs, CPU 4, cap 11000, seed 3863) gives median -2.47%, interval [-3.01%, -0.90%], 25/32 improving pairs (p .00105), and both orders improve. It fails the declared noise-separation rule because AA spans [-1.13%, +0.87%]. Stopped without reserved confirmation/full gates or repeat-until-pass; no qualifying solver speedup is credited. Source, accepted CLI and PGO restored. Original and final source/binary/build, tests, all small/direct counter/result pairs, complete screen schedule/raw results and decision: `target/perf-btran-memo/`, final variant under `v2/`. Fresh RHS-result/estimate reuse is accepted above; this BTRAN variant remains unaccepted.

- Specialize F80's two 128-by-64 integer quotient/remainder steps with x86-64 DIV, retaining portable arithmetic fallback, precondition assertions and the unchanged final F80 pack/rounding. Normalized divisor significands are at least 2^63, so high(m << 63) is below the divisor; the first remainder proves the second quotient also fits u64. The complete previous division is independently frozen (source-normalized copy verified before editing); 100,000 helper quotient/remainder comparisons, 100,000 random extended divisions, captured QTAU/boundary/nonfinite/signed-zero cases and all nine F80 release tests pass. Fast parity 940/940 and eight bounded ADLITTLE/fx_quad_obj_5371 canonical result/exit pairs per model agree. Candidate F80 division has two DIV instructions and shrinks 807 -> 739 bytes. Small instructions are effectively flat on ADLITTLE and -1.93% on fx_quad; median cycles on the latter are nearly flat. Grouped fx solve 0.131176 -> 0.135346 s stays within spread. The predeclared 32-AB/32-AA screen (CPU 4, cap 5500, seed 3751) gives median -0.73%, interval [-1.25%, +0.10%], 21/32 improving pairs (sign p .0551), and biased AA controls (median -0.61%, interval [-1.16%, -0.04%], p .0201). Stopped without reserved confirmation/full gates: no established solver runtime benefit. Source, accepted CLI and PGO metadata restored. Candidate/frozen source and build, disassemblies, tests, small counters, complete result pairs and all controlled schedules/results: `target/perf-divrem/`.

- Packed native extended-number feasibility: the unchanged software F80 independently checks 101,044 captured/random/extreme input pairs, including full low bits, zeros, subnormals, nonfinites and six x87 control modes. Guarded and scoped add/sub/mul/div agree bit for bit; scoped validation retains a per-operation TLS flag and is confined to a probe that controls all calls and CW. Four bounded AB/BA CPU-4 pairs of 640,000-term recurrences reject the final two-u64 layout with cold software fallback: per-operation CW native multiply/add is 64.7% slower and scoped validation 50.4% slower; native division with every input/output conversion is 33.5/37.6% slower respectively. Earlier explicit-padding and outlined-path variants are retained; the first unpinned run is not used for these CPU-4 conclusions. Disassembly identified memory-return/copy overhead in the initial large fallback implementation, but correcting layout and outlining fallback did not establish profitability. These are isolated recurrence results, not whole-solver performance bounds or proof against every packed implementation; no number representation or native dispatch was installed. Sources/binaries, versions, bitwise checks, raw timings and exact build metadata: `target/perf-packed-number/probe-{v1,v2,v3,v4}*` and `probe-manifest.json`.

- Replace ThinLTO with full LTO in both stages of the fixed 39-case PGO build, leaving numeric source, compiler, allocator, training arguments and 64-byte loop alignment unchanged. Explicit `CARGO_PROFILE_RELEASE_LTO=fat` is confirmed by logged rustc commands; the build takes 97.9 s. Candidate CLI shrinks 4,289,304 -> 3,673,904 bytes (-14.35%) and passes fast parity 940/940. Eight bounded CPU-4 ADLITTLE/capped-DEGEN2 canonical result/exit pairs per model agree; instructions rise 0.49/0.66%. PEROLD instructions rise 1.44%; grouped solve 0.189210 -> 0.190470 s stays within spread. The predeclared 32-AB/32-AA screen (CPU 4, cap 1300, seed 3631) preserves every complete result/exit but gives median +1.57%, interval [+1.32%, +1.77%], AA [-0.45%, +0.49%], with 31/32 pairs slower and both execution orders slower. Rejected without reserved confirmation/full gates. Accepted source/CLI/PGO restored. Exact build wrapper/environment and verbose rustc log, candidate and accepted binaries/profile metadata, controlled schedule/raw results, counters and small pairs: `target/perf-fatlto/`. The supplementary `build-settings.json` is required in addition to `candidate-build.json` to reproduce full LTO; no build option or Cargo profile change is retained.

- Skip fresh-factor numerical updates into the Markowitz nucleus only when L is strictly lower and every reached nucleus output is overwritten by a producer U value above tolerance. The cutoff is the first phase-0 pivot, not loff (which precedes phase-4 singletons). Independent frozen-constructor tests compare all fresh factor fields, exact f64/F80 bits and defaults on 26 cached small-model factors per precision plus randomized/extreme/duplicate/nonfinite/missing-output/nontriangular fallback cases; both release tests pass. Fast parity 940/940 and eight ADLITTLE/BORE3D full-result/exit pairs agree. The guarded structural screen covers 492/719 tail columns across the 26 double factors, with 572/655 potential update visits skippable; these are structural counts, not executed arithmetic or attainable time. Small instructions rise 0.10/0.14%; PILOT4 instructions fall 0.57%. Grouped PILOT4 solve 0.074889 -> 0.074054 s stays within spread. The predeclared 32-AB/32-AA screen (CPU 4, cap 800, seed 3511) gives median +0.26%, interval [-0.53%, +1.78%], AA [-2.03%, +1.59%], sign p .811 and mixed execution-order medians. Stopped without reserved confirmation/full gates; no runtime benefit established. Source, accepted CLI and PGO metadata restored. Candidate/frozen source, builds, tests, counters, small/full-result pairs and all timing data: `target/perf-fresh-replay/`; coverage source/reports: `target/perf-fresh-coverage/`. This tests guarded arithmetic elimination while retaining reachability, sorting and entry materialization, not their elimination or direct producer/consumer equality.

- Reuse the existing sorted reverse-row copy in dependency-walk length queries: a tail column absent from rowl[r] has no r entry in A/B, so raw lengths suffice; references retain the exact filter. Move the already-existing clone/sort before the walk (rowl[r] does not change in between), retain its later filtering, and cache each step's length once. Fast parity 940/940 and eight ADLITTLE/BORE3D full-result/exit pairs agree; small instructions -0.08/-0.13%. PILOT4 instructions -1.20%, grouped solve 0.048497 -> 0.047413 s improves beyond both spreads. Two predeclared independent 32-AB/32-AA screens (CPU 4, cap 800, seeds 3371/3433) retain all canonical result/exit pairs, medians -1.51/-1.76%, intervals [-1.80%, -0.92%] / [-3.00%, -0.08%], AA intervals [-0.46%, +0.12%] / [-1.73%, +1.38%]. The first passes. Confirmation has only 23/32 improving pairs (one-sided sign p .0100308, above .01) and fails separation from AA uncertainty. Stopped as timing-inconclusive without extending/repeating the screen or running full gates. Source, accepted CLI and PGO metadata restored. Source/binary/build, logs, counters, both schedules and every controlled raw result: `target/perf-walkrefs/`. This does not establish a qualifying speedup despite promising medians; do not repeat it until it passes.

- Separate dependency-walk candidate validity from its running maximum: replace slots >= cur by hole before taking max, using the invariant nxt >= hole to remove a redundant lower-bound test. Fast parity 940/940 and eight ADLITTLE/BORE3D result/exit pairs agree. PILOT4 instructions -0.27%; grouped solve 0.073487 -> 0.075760 s remains within candidate spread. The predeclared 32-AB/32-AA screen (CPU 4, cap 800, seed 3271) gives median -0.18%, interval [-0.56%, +0.37%], AA [-0.59%, +0.05%], sign p .298 and mixed order medians. Stopped without reserved confirmation/full gates; no established runtime gain. Source, accepted CLI and PGO metadata restored. Candidate source/binary/build, logs, counter pairs, target table, disassembly and exact controlled schedule/raw data: `target/perf-walkmax/`.

- Compute dependency-walk A/B length once per jump and reuse it for the branch and WORK, with no intervening list mutation. Fast parity 940/940 and eight ADLITTLE/BORE3D full-result/exit pairs pass. PILOT4 instructions -0.66%; grouped solve 0.095928 -> 0.098167 s stays within baseline spread. The predeclared 32-AB/32-AA screen (CPU 4, cap 800, seed 3167) gives median -0.30%, interval [-1.88%, +1.35%], AA [-1.44%, +2.74%], sign p .430 and opposite execution-order medians. Stopped without reserved confirmation/full gates: no established runtime gain. Source, accepted CLI and PGO metadata restored. Candidate source/binary/build, logs, small-model counters, target table and exact controlled schedule/raw data: `target/perf-walklen/`. Cached accepted-release update disassembly/annotation identifies separate counts and the slot-max walk; `target/perf-update-current/update.asm`, using `target/perf-primal-current/pilot4.data` with its existing event/count audit. This trial does not change the slot-maximum reduction.

- Maintain coefficient rows containing only required ordinary double-pivot entries, keeping q/leave, all slacks and fixed nonbasics. Stable removal and sorted insertion retain duplicate-entry order; a basis snapshot checks the expected one-column pivot and falls back to complete differences after repairs/reordering. The cost-change path keeps the full product. The old full product is independently frozen: random/extreme/nonfinite/duplicate/rank-tie tests compare both row contents and retained bits across pivots, rotations and multiple-column repairs, with captured scaled-AFIRO inputs. The test and fast parity 940/940 pass. Eight CPU-4 ADLITTLE/BORE3D result/exit pairs agree, instructions +0.05/+0.07%. PILOT4 instructions -1.61%; the grouped solve 0.069438 -> 0.086644 s crosses host bands and initially flags +24.8%. Three direct AB/BA counter pairs retain complete optimal 690-iteration results and do not reproduce that magnitude: cycles -0.23/-0.42/+0.29%, branches about -1.41%, branch misses +1.23/+1.06/+2.58%. The predeclared 32-AB/32-AA screen (CPU 4, cap 800, seed 3061) gives median -0.09%, interval [-0.25%, +0.37%], AA [-0.63%, +0.47%], sign p .298. Stopped without reserved confirmation/full gates: no established runtime gain. Source, accepted CLI and PGO metadata restored. Candidate source/binary/build, logs, initial/target counters, exact direct-pair commands/results and complete controlled schedule/raw data: `target/perf-activerows/`. This establishes numerical feasibility on the checked cases, not profitability; the maintained data and snapshot/update costs remain part of the measured solver.

- Omit unused basic structural outputs only in the double ordinary-pivot ranked product, retaining all nonbasics (including fixed), all slacks and q/leave; the independent cost-change caller keeps the original full product. A copied frozen full kernel checks retained bits on scaled-AFIRO captured rho/ranks plus random duplicate entries/rank ties/signed zeros/extremes/two NaN payloads. The test and fast parity 940/940 pass; bounded ADLITTLE/BORE3D complete-result/exit pairs agree. PILOT4 instructions rise 2.21%, grouped solve 0.048283 -> 0.050302 s (+4.18%) exceeds both 0.000545/0.000208 s spreads. Three CPU-4 AB/BA direct counter pairs preserve every complete PILOT4 result/exit (cap 800, optimal in 690 iterations), with cycles +3.84/+4.15/+4.84%, branches +6.47% each and branch misses +7.15/+8.69/+8.47%. These whole-CLI counters support increased predicate/control cost; they do not attribute all extra cycles to misprediction. Stopped before randomized confirmation/full gates; source, accepted CLI and PGO metadata restored. Candidate source/binary/build, build/tests, small pairs, target logs and all direct counter commands/results: `target/perf-pivotproducts/`. Retained outputs are numerically feasible on these tests, but direct per-entry masking is unprofitable; an eligible-entry traversal avoiding that predicate is a separate untested mechanism.

- Lazy compact double ranked-product rows: flatten immutable row entries into separate u32 index/f64 value arrays, validate destinations once, retain the original sorted rank/entry order and keep extended products on their existing path. The old double implementation is frozen independently. Captured scaled-AFIRO rho inputs and random duplicate entries/rank ties/signed zeros/extremes/two NaN payloads agree bit for bit, including repeated cached calls; both variants pass fast parity 940/940. The first variant retains checked row slices and saves 0.36% PILOT4 instructions; solve 0.048028 -> 0.047851 s stays within spread. The second removes redundant range checks on the private validated cache and guards its width: PILOT4 instructions -0.49%, solve 0.047990 -> 0.047528 s stays within spread. ADLITTLE/BORE3D bounded full-result/exit comparisons agree; v2 instructions +0.02/+0.04%. The predeclared 32-AB/32-AA PILOT4 screen (CPU 4, cap 800, seed 2843) gives median -0.99%, interval [-2.08%, +0.40%], AA [-0.65%, +1.85%], sign p .108, only 20/32 improving pairs. Stopped without reserved confirmation/full gates: no established timing gain. Source, accepted CLI and PGO metadata restored. Source/binary/build snapshots for both variants, build/tests, initial/target counters, exact screen schedule and all controlled raw outputs: `target/perf-rankrows/`. Initial small-model raw files were superseded by v2; v1 summaries remain in `small.log`, with v1 source/binary/build and target counter logs retained. This tests an additional compact validated cache for the double product, not replacing the original row storage or omitting unused output columns.

- Pair terms in the short double Devex dot to reduce scalar length tests: the independent frozen full dot_lanes oracle (random/extreme/duplicate indices, signed zeros, two signed NaN payloads, underflow/overflow, both peels and captured scaled-AFIRO tau inputs) first rejects plain Rust pairing: `[inf * -0, qNaN(0x456) * minsub]` changes the accumulated NaN from `0xfff8000000000000` to `0xfff8000000000456`. Disassembly shows the second add uses the product as destination, reversing the original accumulator-first NaN priority. Explicit register-only SSE2 multiply/add destinations pass the oracle with parent and final PGO flags and fast parity 940/940, but controlled ADLITTLE/BORE3D instructions rise about 0.24/0.10%; PILOT4 +0.60%, DEGEN2 primal +0.29%, PEROLD primal +0.70%. Initial solve changes +3.95/+0.82/+6.33% stay within spread and establish no gain. Source/oracle/fixture restored. Raw `target/perf-shortdot-{fast.txt,small*,initial-*}`, disassembly `target/perf-shortdot/{candidate,frozen}.asm`, tests `/tmp/rose-shortdot-{parent-tests,kernel-tests-parent-flags,sse-kernel-tests-parent-flags,kernel-tests}.txt`; candidate `/tmp/rose-shortdot-candidate{,.rs,-build.json}`, fixture `/tmp/rose-shortdot-fixture.json`. Restored full parity 987/987 and release suite pass (`/tmp/rose-shortdot-restored-{parity,tests}.txt`).

- Batch U-transpose row deletions after factor-update column elimination using the already sorted unique eliminated-column list, retaining a single-column fast path: fast parity 940/940 passes. Small controlled instructions are effectively flat (about 0.00–0.02% higher); PILOT4 -0.03%, DEGEN2 primal -0.01%, PEROLD dual flat. Initial solve changes -1.54/-1.33/+6.78% all stay within spread, with no measurable gain. Source restored; raw `target/perf-utbatch-{fast.txt,small*,initial-*}`; candidate `/tmp/rose-utbatch-candidate{,.rs,-build.json}`. The old repeated retains do not account for a material instruction excess on these bounded cases. Restored full parity 987/987 and release suite pass (`/tmp/rose-utbatch-restored-{parity,tests}.txt`).

- Reuse first A/B lookup positions within update membership, work-booking and A removal, retaining B re-search after mutations: fast parity 940/940 passes. AFIRO instructions are flat; ADLITTLE/BORE3D controlled instructions fall about 0.01–0.04%; three CPU-4 PILOT4 instruction changes are -0.07%, DEGEN2 primal -0.03% and PEROLD dual -0.04%. Every solve change stays within spread (PILOT4 +1.70/-5.15/+7.26%, DEGEN2 +1.93/-2.69/+0.97%, PEROLD +0.65/+2.59/-7.22%). No measurable gain; source restored. Raw `target/perf-abpos-{fast.txt,small*,target-*}`; candidate `/tmp/rose-abpos-candidate{,.rs,-build.json}`. This does not test a maintained reverse A index.

- Reuse the existing BTRAN list order for double primal row products instead of rebuilding/sorting ranks, retaining last-write duplicates and filtering: frozen parent random/extreme/duplicate/captured-AFIRO tests and fast parity 940/940 pass. Controlled small-model instructions fall 0.07–0.08% AFIRO, 0.42–0.44% ADLITTLE and 0.61% BORE3D. Three CPU-4 DEGEN2 primal pairs improve 2.09/2.18/2.28% beyond spread, instructions -0.99%; DEGEN3 primal has two gains beyond spread, instructions -0.95%. PILOT4 instructions fall 0.40%, but all target timings remain within spread. Reject after the 23-case sweep flags primal_5959 primal +2.50% and 80BAU3B PDHG +49.34%; three follow-up pairs stay within broad spreads, yet direct primal_5959 cycles increase 1.06/3.02/6.88% with instructions +0.024%. PDHG direct cycles are mixed; full result objects equal the accepted CLI bit for bit except seconds. No cause is assigned to the quad regression. Source, fixture and oracle restored; raw `target/perf-rankorder-{small-controlled*,target-*,full-*,followup-*,regression/summary.json}`, `/tmp/rose-rankorder-{parent-tests,kernel-tests}.txt`; candidate `/tmp/rose-rankorder-candidate{,.rs,-build.json}`, fixture `/tmp/rose-rankorder-fixture.json`. Restored full parity 987/987 and release suite pass (`/tmp/rose-rankorder-restored-{parity,tests}.txt`).

- Reuse the unchanged full-precision row infeasibility between touch_rows tolerance and priority calculations, preserving the original f64 tolerance conversion: release disassembly confirms duplicate bound tests/subtraction, and the shared value shrinks f64 touch_rows 1182 -> 1122 bytes (F80 5996 -> 5968). Fast parity 940/940 passes, but native PEROLD instructions fall only 0.18%, DEGEN2 dual 0.10–0.11%, DEGEN3 dual 0.07%; bounded small counts are flat to -0.09%. Three alternating PEROLD solve changes +22.83/-21.87/+0.64% cross host speed bands and establish no repeatable gain; DEGEN2 +19.81/+0.89/-1.46% and DEGEN3 +21.44/+4.05/-0.36% include initial regressions and no gain beyond spread. Source restored (`target/perf-touchvalue-{small-summary.json,target-*,fast.txt,kernels.txt}`, `target/perf-touchrows-before.asm`; candidate `/tmp/rose-touchvalue-candidate{,.rs,-build.json}`).

- Defer both unused individual U estimates through shared runtime-policy kernels, testing stored Option estimates, unchanged factor fields with tuple returns, then unchanged fields/vector return with an L-work output parameter: all three pass frozen f64/F80 random/extreme/captured FTRAN tests and fast parity 940/940; emitted kernels have one instance per numeric type. Final PEROLD instructions fall 1.61–1.62%, DEGEN2 0.97–0.98%, DEGEN3 1.19%. Final PEROLD CPU-4 solve changes -0.18/+8.44/-0.92% have one gain beyond spread; CPU-2 confirmations -1.21/-2.55/-46.22% do not establish repeatability (first two within spread, last crosses speed bands). DEGEN3 has two comparable gains beyond spread, but other pairs stay within spread or regress +28.63%. Final mode-switch instructions are +0.012%, cycles +0.34/+1.38/-1.53%, and paired timing -0.15/-32.36/+0.43% (large middle shift uncredited); all complete results equal the accepted binary apart from seconds, 11393 iterations/17 refactors. Stored-Option mode-switch cycles +5.29/+1.69/+2.25% and tuple-return -1.00/+4.51/+18.34% also fail to establish an acceptable change. No repeatable primary-target gain; source/tests restored. Sharing kernels and restoring field/return layouts does not by itself qualify deferral (`target/perf-{runtimework,returnwork,outwork}-{kernels.txt,fast.txt,target-*,confirm-*,mode-*}`, `/tmp/rose-{runtimework-option,returnwork-candidate,outwork-candidate}/`, `/tmp/rose-{runtimework,returnwork,outwork}-kernel-tests.txt`).

- Share sparse eta accounting with the f64 numeric pass after explicit skipped-prefix evaluation proves each target starts and remains positive zero (analysis idea 2): a fresh exact-parent L-stage oracle, explicit nonempty-prefix/sub-tolerance/NaN/signed-zero/F80-fallback cases, captured SC50A tests and fast parity 940/940 pass. Three alternating CPU-4 PEROLD solve changes +4.03/+6.66/+30.25% establish no gain (last beyond spread), despite instructions -0.41%. DEGEN2 instructions fall only 0.06–0.07%, with solve -0.18/-6.07/-0.25% (only the middle beyond spread); DEGEN3 instructions -0.04%, solve -1.79/-4.88/-7.95% (only the middle beyond spread). Source and temporary oracle/tests restored; no repeatable target timing gain (`target/perf-etaprefix-{small-summary.json,target-*,fast.txt}`, `/tmp/rose-etaprefix-{parent-tests.txt,kernel-tests.txt,candidate.rs,candidate,candidate-build.json}`). This tests a sufficient unchanged-prefix predicate, not sharing after prefixes that alter residues.

- Reuse dual factorization input-column row/value buffers across refactors, growth retries and repairs, retaining model-entry order and unchanged numeric kernels: fast parity 940/940 passes. PEROLD malloc calls fall by 30,024, offset by 2,024 additional reallocations (28,000 net, 3.53%); instructions fall 0.35% (DEGEN2 dual 0.49–0.50%, DEGEN3 dual 0.31%). Three alternating CPU-4 PEROLD solve changes +5.22/-0.68/+0.60% include a regression beyond spread; DEGEN2 +1.37/-0.33/+0.34% all stay within spread, and DEGEN3 +0.44/+2.12/+28.00% includes two regressions beyond spread. No target timing gain; source restored (`target/perf-inputreuse-{small-summary.json,target-*,fast.txt}`, `target/perf-native-alloc/inputreuse-*`; candidate `/tmp/rose-inputreuse-candidate{,.rs,-build.json}`).

- Store auxiliary numeric LU L/U in flat buffers with column offsets, eliminating per-column vectors while preserving every entry and numeric operation: release LU tests and fast parity 940/940 pass. Native PEROLD allocation/reallocation calls fall by 34,201 (4.31%), but instructions fall only 0.27% (DEGEN2 0.46%, DEGEN3 dual 0.42–0.43%). Initial PEROLD solve changes -19.36/-0.10/-0.59% include a host speed-band shift; CPU-2 confirmations -0.70/-0.04/-0.80% all remain within spread. DEGEN3 has one qualifying -0.82% confirmation, but its other two and all three initial pairs remain within spread. No repeatable target timing gain; source restored (`target/perf-flatlu-{small-summary.json,target-*,confirm-*,fast.txt}`, `target/perf-native-alloc/flatlu-*`; candidate `/tmp/rose-flatlu-candidate{,.rs,-build.json}`).

- Elide internal-index U scatter bounds checks in first/second-vector FTRAN: frozen checked-kernel tests cover f64/F80 random/extreme/NaN/duplicate/captured SC50A inputs; all seven FTRAN tests and fast parity 940/940 pass for inline and outlined variants. Inline PEROLD instructions fall 2.08%, but initial solve changes -1.10/+10.48/+0.55% and CPU-2 confirmations +2.93/-0.12/+30.78% establish no repeatable gain and include regressions beyond spread. Outlined helpers reduce PEROLD instructions 1.45% and DEGEN3 dual 1.69–1.70%, but all three paired timings remain within spread (PEROLD -0.04/+2.95/+0.66%; DEGEN3 +2.95/-36.74/-0.12%, with a host speed-band shift). Source and temporary oracles restored; accepted binary unchanged (`target/perf-uscatter-{target-*,confirm-*,outline-target-*,outline-fast.txt}`, `/tmp/rose-uscatter-{tests.txt,outline-tests.txt}`).

- Cache ordered tail-only numeric U entries for the second-vector scatter: a fresh frozen complete ftran_needed oracle and cache consistency tests cover f64/F80, captured SC50A, random/extreme inputs, duplicate values, signed zeros/NaN payloads, column replacements and transient/final row-class changes; fast parity 940/940 passes. Net class changes refresh through the maintained U transpose at update completion. PEROLD instructions fall about 1.11%, DEGEN2 0.95% and DEGEN3 dual 2.78%; small-model counts range from flat to +0.3%. Three CPU-4 PEROLD timing changes -1.68/+22.94/+0.61% and CPU-2 confirmations +3.03/-21.39/-1.64% all stay within spread. DEGEN3 has one initial -5.27% pair beyond spread, but its other initial pairs and two confirmations remain within spread; its -23.21% confirmation crosses speed bands and is not credited. No repeatable target timing gain; source/oracles/cache tests restored (`target/perf-tailcache-{small-summary.json,target-*,confirm-*,fast.txt}`, `/tmp/rose-tailcache-{tests.txt,candidate,candidate.rs}`).

- Reuse the existing structural-column offset prefix for primal pivot-row order, row-product work and Devex work (analysis idea 8), with the original zero-row scan fallback: fast parity 940/940 passes; PILOT4 instructions fall 0.40–0.42%, but three CPU-4 solve changes -0.31/-4.22/+3.58% all stay within spread. DEGEN3 primal instructions fall about 0.35%, with +6.57/-1.46/-1.13% timing changes also within spread. Bounded small-model counts are flat. No qualifying timing gain; source restored (`target/perf-nnz-{small-summary.json,target-*,fast.txt}`, `/tmp/rose-nnz-candidate{,.rs}`).

- Retain dense L storage through U, replacing U's zero allocation plus spike scatter with a transferred buffer; sparse filtering stays unchanged and omitted signed zeros are normalized to positive zero. Existing frozen L/complete-FTRAN oracles are untouched; extended f64/F80 random/extreme/NaN/signed-zero/captured-SC50A representation checks and all six FTRAN/L release tests pass, as does fast parity 940/940. Eight bounded CPU-4 AB/BA small pairs preserve complete no-dump results; instructions ADLITTLE +0.28%, capped DEGEN2 +0.82%. Three alternating CPU-4 five-run PEROLD pairs regress solve 1.69/2.46/1.56%, all beyond spreads, instructions +1.61% each. Rejected before full benchmarking; restored source, accepted CLI and PGO metadata. Keeping dense storage still materializes the ordered spike and adds a normalization pass; avoiding the allocation/scatter alone is insufficient. Candidate source/binary/build `/tmp/rose-lu-transfer-candidate{.rs,,-build.json}`; build/test logs `/tmp/rose-lu-transfer-{build,tests}.log`; pairs and counters `target/perf-lu-transfer/{transfer-pairs.json,perold-*-{before,after}.txt}`.

- Transfer solved FTRAN `x` and input support into F2 after output extraction instead of cloning both arrays: fast parity 940/940 passes. ADLITTLE/capped-DEGEN2 eight AB/BA pairs show instructions +0.036/+0.026%; their timings overlap the fast-parity job and are diagnostic only. Sequential CPU-unpinned PEROLD three-run screen has solve 0.187584 -> 0.190345 s within spreads 0.051039/0.007805 s, instructions +0.424%. No qualifying gain; stopped before full benchmark and restored source, accepted CLI and PGO metadata. This tests retained-accounting ownership only, not dense L-to-U transfer. Candidate source/binary/build `/tmp/rose-factor-transfer-candidate{.rs,,-build.json}`; evidence `target/perf-factor-audit/{transfer-fast.txt,transfer-pairs.json,transfer-perold-{before,after}.txt}`.

- Reuse fresh-factor per-tail-column start, graph-discovery and sorted-order vectors: fast parity 940/940 passes and PILOT4 allocation/reallocation calls fall by 15,723 (about 5.3%), but instructions only fall 0.51–0.55%. Three initial CPU-4 PILOT4 solve changes -3.02/-0.10/-3.08% and three CPU-2 confirmations -2.17/-3.77/-1.32% all stay within spread; DEGEN3 primal instructions fall about 0.69%, with all three timing changes within spread too. No qualifying timing gain; source restored. Raw `target/perf-factor-scratch-{small-summary.json,target-*,confirm-*,after.alloc}`, baseline allocation totals `target/perf-factor-buckets-before.alloc`; candidate `/tmp/rose-factor-scratch-candidate{,.rs}`.

- Stable compact temporary fresh-factor Markowitz-column and tail-row membership adjacency: fast parity 940/940 passes. PILOT4 instructions fall 0.51–0.55%, but three CPU-4 solve changes +0.55/-0.69/+2.99% establish no gain (last beyond spread); DEGEN3 primal instructions fall about 0.32% with all timing changes within spread. ADLITTLE/BORE3D instructions fall about 0.5%; AFIRO is flat. Source restored; original/candidate allocation totals `target/perf-factor-buckets-alloc.json`, pairs `target/perf-factor-buckets-{small-summary.json,target-*}`, candidate `/tmp/rose-factor-buckets-candidate{,.rs}`.

- Native extended ranked products with register-held products under the existing control-word guard: frozen native/software differential tests and fast parity 940/940 pass. Register-only CPU-4 confirmations reduce Quad_2717 instructions 0.22–0.26%, but solve changes -1.59/-1.12/+0.97% all remain within spread; FX instructions fall 0.45% with no gain beyond spread. Separately outlining the native pass reduces Quad instructions 0.19–0.28%; CPU-2 solve changes -1.81/-3.80/-3.37%, only the last beyond spread. No repeatable timing gain established; both variants restored. Evidence: `target/perf-rankreg-confirm-*` and `/tmp/rose-rankreg-current-kernel-tests.txt`, `target/perf-rankoutline-{fast.txt,target-*,confirm-*}`; register-only source/binary retained in `/tmp/rose-rankreg-register-only{.rs,-build.json}` and `/tmp/rose-rankreg-register-only`.

- CLI-local bounded small-allocation cache (16-byte classes through 1024 bytes, 32 blocks/class/thread, System fallback and thread-exit reclamation): full parity 987/987, complete PGO release tests and additional zeroing/alignment/reallocation/cross-thread/teardown tests pass. Three alternating CPU-4 target pairs improve PEROLD 1.77–2.27% with instructions about -3.4%, DEGEN2 5.26–5.58% (-7.26%) and DEGEN3 3.08–3.46% (-4.05%). Rejected by unresolved mode-switch regression: all-23 sweep flags +50.68%; first CPU-2 follow-up is within spread, but fixed three-pair confirmation gives +4.13% beyond spread, +1.72% within spread, then a -32.68% timing-band shift. Three cycle pairs add 0.66/1.21/1.81% despite instructions -0.63%; a no-regression result is not established. PDHG also crosses wide timing bands with slightly more instructions (+0.25/+0.39%). Source/tests restored; no speedup credited (`target/perf-smallalloc-{target-*,confirm-*,full-*,followup-*,long-confirm-*,counters.json}`, `/tmp/rose-smallalloc-{full-parity.txt,tests.txt,cleanup-test.txt,source/}`).

- Replace sparse-L per-level vectors with flat FIFO successor links, preserving insertion order and all arithmetic/work expressions: fast parity 940/940 passes. Bounded small-model instructions are effectively flat. PEROLD whole-CLI allocation/reallocation calls fall by 21,297, but instructions only fall about 0.4%; all three alternating CPU-4 solve changes (-1.16/-4.26/-2.68%) stay within their spreads. CPU-2 has one qualifying pair and two within spread; no repeatable target timing gain. Source restored (`target/perf-levelqueue-{fast.txt,small-summary.json,target-*,confirm-*,alternating-*}`, `/tmp/rose-alloc-profile/levelqueue-{before,after}.alloc`).

- Construct nucleus row chains by reversing column/ID traversal and appending, exactly reproducing front-insertion order without column-chain clones or shifts: fast parity 940/940 passes. Bounded ADLITTLE instructions fall 0.16–0.20%; AFIRO varies around flat and BORE3D rises about 0.09%. PEROLD instructions fall 0.16–0.23%, DEGEN2 dual about 0.15–0.16%, DEGEN3 dual about 0.14%, but three alternating CPU-4 PEROLD solve changes +0.49/+0.13/+0.56% all stay within spread. DEGEN2 has one gain pair, with its other pairs and all DEGEN3 pairs within spread. No repeatable timing gain on the target; source restored (`target/perf-rowbuild-{fast.txt,small-summary.json,target-*,confirm-*}`, `/tmp/rose-rowbuild-candidate.rs`).

- Combine pivot-row values and touched flags in one indexed tuple buffer, retaining every numeric operation, first-touch order and count: fast parity 940/940 passes, but PEROLD instructions rise 0.68%, DEGEN2 dual 0.59%, DEGEN3 dual 0.52%. Initial solve changes +2.21/-7.68/+2.94%; PEROLD's increase exceeds spread and no acceptable gain is established. A bounded auxiliary f64 constructor disassembly already uses zeroed allocation and aliases the bool and word-marker constructors, so a missing zero-initialization fast path is not supported as the cause. Source restored (`target/perf-rowcells-{fast.txt,target-*,iteration.asm,init.asm}`, `/tmp/rose-rowcells-candidate.rs`).

- Keep pivot_row separate from the large dual iteration function to improve scatter register allocation: fast parity 940/940 passes, but PEROLD instructions rise 0.17%, DEGEN2 dual 0.12%, DEGEN3 dual 0.13%. Initial solve changes +42.13/+5.33/-1.66% include host speed-band shifts and establish no acceptable gain. The current PEROLD profile samples repeated metadata loads at the inlined scatter, but removing that context does not reduce total instructions. Attribute restored; candidate source and standalone assembly retained (`target/perf-pivotoutline-{fast.txt,target-*,kernel.asm}`; `/tmp/rose-pivotoutline-candidate.rs`).

- Defer the column U estimate to the single-vector fallback and omit the unused tau estimate on joint FTRAN paths, with explicit unavailable estimates and separate const-generic column kernels: frozen f64/F80 random/extreme/captured tests compare deferred column estimates bit for bit and retain tau counts; full parity 987/987 and the full PGO release suite pass. DEGEN2 instructions fall about 0.73%; comparable CPU-2 solve pairs 0.032116 -> 0.031811 s and 0.032348 -> 0.031618 s exceed their spreads. The 23-case sweep also improves DEGEN3 dual-release 0.273587 -> 0.270648 s (-0.93% instructions) and PEROLD dual-release 0.131784 -> 0.130330 s (-1.25% instructions). Wide speed-band shifts in other pairs are not credited. Reject because final mode-switch counter pairs increase cycles 1.86/6.07/1.80% and instructions 0.50/0.51/0.51%, despite earlier timing pairs staying within spread. DEGEN3 primal's large sweep timing flag does not repeat beyond spread in three follow-ups; final DEGEN2 cycles change -0.46/-4.45/+0.23%, so its isolated timing increase is not repeated consistently. Source/test changes restored (`target/perf-pairwork-{fast.txt,target-*,confirm-*,recheck-*,full-*,followup-*,final-counters.json,perold-counters.json}`; `/tmp/rose-pairwork-{full-parity,tests,kernel-tests}.txt`; exact benchmark CLI/source `/tmp/rose-pairwork-candidate{,.rs}`). Shared runtime-policy, factor-field and return-layout variants are also rejected above; none establishes a repeatable primary-target gain.

- Skip the unused individual tau U estimate when the joint U estimate is booked, retaining the original selector and exact tail result count: frozen f64/F80 random/extreme/captured tests and fast parity 940/940 pass. Three alternating CPU-4 pairs reduce DEGEN2 instructions about 0.45%, DEGEN3 0.57%, PEROLD 0.80%; DEGEN2 solve changes -0.25/-0.20/-0.15%, all within spread. Only one DEGEN3 pair exceeds its spread; the other DEGEN3 and all PEROLD pairs do not. Median solve medians DEGEN2 0.032242 -> 0.032161 s, DEGEN3 0.272987 -> 0.271977 s, PEROLD 0.132004 -> 0.131314 s do not establish a repeatable gain. Source/test changes restored (`target/perf-tauwork-{fast.txt,target-*,confirm-*}`, `/tmp/rose-tauwork-kernel-tests.txt`; candidate `/tmp/rose-tauwork-candidate.rs`). The audit confirms no other production reader of tau's individual estimate; the column estimate still serves standalone FTRAN callers and was not deferred in this trial.

- Share a sparse column's reachability/output order with the dense-tau slot walk (analysis idea 1), retaining separate scalar arrays and all individual/joint work estimates: frozen parent numeric and output-order oracles, f64/F80 random/exceptional/captured tests, an explicitly active sparse case and late-dependency fallback pass; fast parity 940/940 and targeted full dumps for DEGEN2/DEGEN3/PEROLD dual-release agree. DEGEN2 instructions fall about 0.20%, PEROLD 1.08%; DEGEN3 rises about 0.05–0.09%. Fixed three-pair CPU-4 confirmation has no timing gain beyond spread; large apparent shifts cross timing bands and are not credited. A bounded forty-solve release profile samples the new kernel at 1.12% self, so it is active. Source/helper/oracle changes restored (`target/perf-sparsejoint-{small.txt,target-*,confirm-*,self.data,self.txt}`; candidate `/tmp/rose-sparsejoint-candidate.rs`). This settles the sparse-column/dense-tau joint-walk lead; it does not test deferring unused individual work estimates.

- Interleave eligible f64 FTRAN vectors for packed dense scatter/division (analysis idea 1), with scalar replay for NaNs and unchanged sparse/quad paths: parent kernels copied exactly into frozen oracles; f64/F80 random/extreme/two-NaN-payload/finite-overflow/native-F80/captured-SC50A tests compare outputs, snapshots used by joint work, estimates and counts exactly. Tail-only packing passes fast 940/940 parity and targeted DEGEN2 full-dump parity, but initially adds DEGEN2/DEGEN3 instructions 1.29/2.01%; vectorized unpacking and packed division reduce this to +0.05/+0.20%, with all final target timings within spread. Carrying the interleaved buffer through dense L, etas and U also passes frozen tests and fast parity, but stable DEGEN2 solve 0.032377 -> 0.032791 s (+1.28%) exceeds 0.000115/0.000281 s spreads, instructions +0.18%; DEGEN3 0.272676 -> 0.275899 s (+1.18%, spreads 0.001136/0.002931 s), instructions +0.36%. PEROLD instructions/timing are effectively flat. Source/oracle/test changes restored (`target/perf-pairpack-{target-*,unzip-target-*,div-target-*,fullstages-target-*,*-kernel.asm,self.data,annotate.txt}`; final candidate `/tmp/rose-pairpack-fullstages-candidate.rs`). This tests dense paired stages; sparse-column/dense-tau traversal is rejected above.

- Find the highest-slot empty tail column by reverse slot traversal and early termination instead of a full-column maximum scan: slot-permutation uniqueness preserves the chosen column, and fast parity 940/940 passes. DEGEN2 dual instructions fall about 0.20%, primal about 0.17–0.24%, PEROLD dual 0.31%. All initial and fixed three-pair CPU-4 confirmation timings stay within their spreads: median solve medians DEGEN2 dual 0.049636 -> 0.048572 s, primal 0.058258 -> 0.059844 s, PEROLD dual 0.266693 -> 0.267239 s. No qualifying timing gain; source restored (`target/perf-convscan-{small.txt,target-*,confirm-*}`).

- Emit numeric-LU touched rows from an indexed bitset instead of comparison sorting, preserving ascending numerical traversal and the existing pending-step heap: independent ordered-set membership/drain tests and LU release tests pass; both variants pass fast parity 940/940. Replacing byte markers entirely adds DEGEN2 dual instructions 0.14%, primal 0.19% and PEROLD dual 0.41%, with all target timing changes within spread. Retaining byte markers for repeated coefficient touches still adds DEGEN2 dual 0.17% and primal 0.31%; dual timing -1.22% stays within spread. PEROLD instructions fall 0.24% in that form but timing -0.69% stays within spread. Both source variants and helper-test scaffolding restored (`target/perf-factorrows-{small.txt,target-*,hybrid-small.txt,hybrid-target-*}`). This concerns touched-row ordering; analysis idea 3's linked pivot buckets are already implemented in `0cb3b6c3`.

- Convert coefficient-list lengths through signed i64 for FTRAN/BTRAN work estimates, using the non-ZST slice-size bound while retaining each floating-point addition and its order: the single dense-L cast passes frozen f64/F80 random/extreme/captured-input tests and fast parity 940/940. Disassembly replaces one unsigned conversion sequence with one signed conversion at unchanged L-stage address/size. DEGEN2 instructions fall about 0.16%, but three alternating CPU-4 pairs change solve -0.20/-0.16/+1.11%, all within spread. Broadening to 22 other work-list casts passes FTRAN differential tests and fast parity; target DEGEN2 0.051517 -> 0.050306 s stays within 0.002494/0.011824 s spreads with instructions -0.17%; PEROLD 0.270595 -> 0.269949 s is also within spread (-0.13% instructions). Both variants and frozen-test updates restored (`target/perf-lencast-{small.txt,target-*,confirm-*,layout.json,*-lstage.asm,wide-target-*}`).

- Skip the FTRAN output fallback scan after dense output ordering has set every done flag: generic shortcut passes full parity 987/987 and PGO release tests; DEGEN2 instructions -0.53%, DEGEN3 dual -0.82%. Two of three CPU-4 target pairs improve beyond spread (DEGEN3 0.428918 -> 0.421134 s and 0.422879 -> 0.415789 s), and DEGEN3 cycles fall 0.88–1.30% in three counter pairs. Rejected by the all-23 paired sweep, fixed three-pair follow-ups and counters: quad-heavy cases repeat roughly 6–11% more cycles despite nearly flat instructions. Restricting the shortcut to f64 passes fast 940/940 parity and reduces that collateral cost, but Quad_2717 still adds 2.04/3.08/4.00% cycles and fx_quad adds 1.79/3.17/4.91% in alternating pairs. Both source variants restored; the cause of the remaining code-placement/code-generation sensitivity is unproved (`target/perf-ftfallback-{target-*,confirm-*,counters.json,full-*,followup-*,regression-counters.json,f64-*}`).

- Enumerate FTRAN markers directly from safe packed eight-byte boolean words (at least 128 positions, current aligned PGO baseline): fast parity 940/940 passes, DEGEN2 instructions fall about 0.83%, but CPU-2 solve 0.048298 -> 0.048804 s has no gain. Three alternating CPU-4 pairs have median solve medians 0.049778 -> 0.049636 s within aggregate 0.003476/0.003237 s spreads. Restored source and helper-test scaffolding (`target/perf-markerbits-{small.txt,target-*,confirm-*}`); no speedup credited.

- Reuse scratch across the dual iteration (analysis idea 6): row accumulators/order, ratio/contribution/index arrays, framework construction, FTRAN inputs, tau expansion and bound-flip buffers. Both inline and boxed workspace variants pass fast parity 940/940. Inline DEGEN2 0.049531 -> 0.049507 s is within 0.000578/0.003967 s spreads; PEROLD slows 0.199162 -> 0.209494 s with instructions +0.70%. Boxing avoids moving the workspace on retries, but fixed three-pair CPU-4 follow-ups still do not establish a repeatable DEGEN2 gain beyond spread (solve medians 0.032353 -> 0.032106 s, instructions -0.79%); PEROLD medians 0.131737 -> 0.132612 s and instructions +0.65%. Owned pivot/output buffers and forced/shift-selection buffers were not recycled in this trial. Restored source; no acceptance claim (`target/perf-workspace-{small.txt,target-*,box-target-*,box-confirm-*}`).

- Force 64-byte function entries in both PGO stages (analysis idea 20), retaining the solver source and 39-case online training: fast parity 940/940 passes and FTRAN has the same code size, but instructions remain effectively flat and target timings do not establish a standalone gain. GREENBEA barrier 0.838744 -> 0.854044 s (+1.82%) exceeds 0.005861/0.001459 s spreads. Function alignment does not ensure aligned internal loop heads; this flag is not used; the landed build aligns loop heads directly (`target/perf-align64-{small.json,target-*}`).

- Skip zero blocks in the two FTRAN marker scans: the all-size boolean-reduction iterator passes fast parity/FTRAN differential tests and lowers mode-switch instructions about 9.9%, but adds 0.6–0.7% instructions on small models. A long-array callback guard preserves small-model counts but its reduction compiles scalar, with only about 3% fewer mode-switch instructions and no timing gain. The final safe eight-byte word checks, enabled at 1024 markers with 64-marker blocks, pass full parity 987/987, release tests and word/block/tail boundary tests. Mode-switch instructions fall about 9.2%; comparable pairs improve 1.853019 -> 1.795435 s (-3.11%, spreads 0.007134/0.015916 s) and 2.767822 -> 2.725229 s (-1.54%, spreads 0.040573/0.035595 s). Rejected by the all-23 paired sweep and follow-ups: GREENBEA barrier repeats 0.846091 -> 0.857920 s (+1.40%, spreads 0.002392/0.007563 s) and, on another CPU, 0.841683 -> 0.853643 s (+1.42%, spreads 0.003414/0.005546 s), with flat instructions; three counter pairs add about 0.3–0.8% cycles. DEGEN3 dual's timing flag does not repeat in its final comparable pair. Source and test scaffolding restored; no table speedup credited (`target/perf-markerblocks-{target,small.json}*`, `target/perf-markerblocks-long-target-*`, `target/perf-markerwords-*`; final source `/tmp/rose-markerwords-candidate.rs`).

- Primal precomputed-row-product prefilter, retested with consistent online PGO merging: fast 940/940, full parity 987/987 and PGO release tests pass for both forms, each followed by an all-23 paired sweep. The iterator form improves DEGEN3 primal-controlled in three alternating follow-ups (3.3–5.0%, about 3.35% fewer instructions), but mode-switch remains unresolved: two of three counter pairs add about 2.5–2.7% cycles with flat instructions. A direct early continue improves a stable DEGEN3 pair 0.514192 -> 0.496084 s (-3.52%, spreads 0.003268/0.001467 s; instructions about -2.24%) and has mostly flat mode-switch counter rechecks, but repeats other regressions: 80BAU3B dual 0.366054 -> 0.370932 s (+1.33%) and 0.560780 -> 0.577604 s (+3.00%), both beyond spread; GREENBEA barrier 0.841948 -> 0.853973 s (+1.43%, spreads 0.000916/0.003669 s), also beyond spread in its first follow-up. These regressions have effectively unchanged instructions. Existing `ord` dumps show no primal row updates in mode-switch, so its timing cannot be attributed to executing the filter. Both source forms restored; no speedup credited. This settles the earlier lead whose full sweep used inconsistent profile collection (`target/perf-rowprefilter-{retest,early}-*`).

- Cache immutable fixed/lower-infinite/upper-infinite flags in the existing one-byte fixed-variable array (analysis idea 8), reusing them in pricing, feasibility classification and primal control: fast 940/940, full parity 987/987 and PGO release tests pass with consistent online profile merging. Stable CPU-2 PILOT4 pairs improve 0.055503 -> 0.054886 s and 0.055530 -> 0.054632 s (1.1–1.6%, beyond both spreads); instructions fall about 1.26%. PEROLD primal also improves in stable pairs with about 1.11% fewer instructions. Rejected by the all-23 paired regression sweep and alternating CPU-4 follow-ups: mode-switch 1.800502 -> 1.901218 s (+5.59%, spreads 0.017340/0.014060 s), with effectively unchanged instructions; the first follow-up also regresses beyond spread. DEGEN3 dual flags stay within follow-up spreads; quad_2717 has mixed follow-up evidence. Source restored, no table speedup credited (`target/perf-boundflags-all-{target,confirm,full,followup}-{before,after}.txt`; small data is `target/perf-boundflags-all-small.txt`).

- Separate free-variable byte array used only by score (analysis idea 8), retested with consistent online PGO merging: fast parity 940/940 passes, but PILOT4 instructions fall only about 0.10% and median solve medians 0.085830 -> 0.084329 s stay within aggregate 0.007276/0.012024 s spreads. PEROLD instructions fall about 0.14% and 0.845668 -> 0.841948 s stays within 0.090099/0.066434 s spreads. Reverted (`target/perf-freecache-retest-{small,target-before,target-after}.txt`). This settles the earlier cache lead whose comparison used inconsistent profile collection.

- Numeric LU pending pivot steps in a bitset priority queue instead of the binary heap/markers: unique ascending ordering matches an independent ordered set over duplicate, backward and sparse-word insertions, and fast parity 940/940 passes. PILOT4 instructions fall only about 0.1%; three alternating CPU-2 pairs give median solve medians 0.086697 -> 0.088784 s inside 0.014614/0.031085 s aggregate spreads, with no gain. ADLITTLE is also within spread. Reverted, including queue-test scaffolding (`target/perf-stepqueue-*`).

- Interleaved forward L destinations/coefficients in separate per-column vectors: frozen f64/F80 random/extreme/captured tests and fast parity 940/940 pass. DEGEN2 dual instructions fall about 1.1%, or 1.7% after eliminating the temporary old-value clone, but three alternating CPU-4 pairs for the latter give median solve medians 0.049033 -> 0.049886 s, within 0.002410/0.003751 s aggregate spreads. Nested-vector variants are not accepted; the flat-buffer variant is measured separately (`target/perf-lforward{,2}-target-{before,after}.txt`).

- Share the dense L traversal of entering/pricing FTRAN vectors when both initially select dense L (analysis idea 1), retaining separate eta passes and U solves: frozen f64/F80 random/extreme/captured SC50A tests and fast parity 940/940 pass. Initial version adds about 0.25% DEGEN2 dual-release instructions, with no timing gain. Combining NaN detection with the required nonzero-count scans lowers instructions only about 0.04%; three alternating CPU-4 pairs give median solve medians 0.050441 -> 0.050018 s, within aggregate 0.004460/0.006422 s spreads. Initial interleaving changes NaN payload priority, requiring the original sequential stages whenever dense L produces a NaN, including finite-input overflow. Candidate self profile samples the shared path (2.27% self), so it is active but saves too little. Both variants and test scaffolding reverted (`target/perf-pairl{,2}-target-{before,after}.txt`, `target/perf-pairl-self.txt`). Sparse-to-dense fallback and other jointly traversable stages remain untested.


- Share each column broadcast between two independent row vectors in an eight-row/four-column S tile, retaining four-row and scalar two-column tails and every output's FMA order, block boundary and final C addition. The accepted four-row SIMD driver is frozen before rewriting; scalar/grouped random/extreme/captured-input comparisons and fast parity 940/940 pass. Eight BORE3D canonical result/exit pairs agree, instructions -0.30%, cycles -0.09%. DEGEN3 barrier instructions -1.46%; grouped solve 0.180750 -> 0.178481 s initially improves beyond spread. The predeclared 32-AB/32-AA screen (CPU 4, cap 100, seed 2749) gives median -1.15%, interval [-2.00%, +0.07%], AA [-1.89%, +0.91%], sign p .0251. Stopped as timing-inconclusive without the reserved seed or full gates; source, accepted CLI and PGO metadata restored. Candidate source/binary/build, logs, counter pairs, schedule/raw samples and disassembly: `target/perf-gemm8/`. Instruction reduction alone is insufficient; wider output-row tiling has not established another speedup.

- Pack the accepted four-row/two-column S tile into two vectors, one independent output row per lane, with unchanged FMA term order and final C updates. The accepted scalar tile and complete driver are independently frozen; scalar/grouped random/extreme/captured-AFIRO comparisons and fast parity 940/940 pass. Eight BORE3D canonical result/exit pairs agree, instructions -0.61%, cycles +0.78%. DEGEN3 barrier instructions -3.10%; grouped solve 0.279525 -> 0.275116 s initially improves beyond spread. The predeclared 32-AB/32-AA screen (CPU 4, cap 100, seed 2459) gives median -0.46%, interval [-1.70%, +1.37%], AA [-1.59%, +1.31%], sign p .189 and opposite order medians. Stopped without reserved confirmation/full gates: no established runtime gain. Source, accepted CLI and PGO metadata restored. Candidate source/binary/build, logs, counter pairs, exact controlled samples and disassembly: `target/perf-gemm-lanes/`. Accepted-release instruction sampling (BORE3D first, eight DEGEN3 solves capped at 100, period 8M) has 2806 samples, no lost/throttle events and represented instructions 99.90% of the 22.470B direct batch count; DGEMM 19.81%, numeric factor 30.40%, DTRSM 9.76% self. These are batch diagnostic fractions including launcher/parsing/startup, not solve-only/removable bounds. Disassembly confirms two dependent packed accumulator chains; instruction reduction alone does not establish a benefit.

- Bound double dense-LDL input rows once and traverse eight-element array chunks, preserving all lane chains, scalar tails and the reduction tree. The current reduction is frozen independently before rewriting; old/new random, signed-zero/extreme/nonfinite and captured AFIRO-wrapper array tests pass. Initial auto-vectorization changes NaN payload priority; final code replays the original scalar reduction on nonfinite results, after which frozen tests and fast parity 940/940 pass. Eight CPU-4 BORE3D canonical result/exit pairs agree, instructions -1.35%, cycles effectively flat. DEGEN3 barrier instructions -3.13%; grouped solve 0.292376 -> 0.294792 s has no gain beyond spread. The predeclared randomized screen (32 AB + 32 AA pairs, CPU 4, cap 100, seed 2131) gives median -0.25%, conservative interval [-1.00%, +0.20%], AA interval [-0.28%, +0.81%], sign p .0551. Stopped for no established timing benefit without running the reserved independent seed or full gates. Disassembly confirms new packed multiply/add instructions, but instruction reduction alone is insufficient. Source, accepted CLI and PGO metadata restored. Candidate source `/tmp/rose-ldlchunks-candidate.rs`; exact measured binary hashes are retained in the paired report (binary/manifest archival hit the scratch quota, so the partial files were removed). Rebuild with `python3 -B perf/build_pgo.py --system-tcmalloc`; logs `/tmp/rose-ldlchunks-{build,build2,tests,tests2}.log`; raw `target/perf-ldlchunks/{transfer-pairs.json,target-*,randomized-degen3-2131.*,numeric-*.asm}`. This tests bounded chunk traversal with nonfinite replay, not batching independent LDL outputs or changing reduction order.

- Hoist noncontiguous child-update destination rows and traverse bounded source slices in lockstep, preserving scatter order and every addition/rounding without adding cached maps. Fast parity 940/940 and eight CPU-4 BORE3D canonical full-result/exit pairs pass. BORE3D instructions -0.23%; DEGEN3 barrier -1.25%. The initial grouped target solve 0.292165 -> 0.283405 s stays within candidate spread. Two predeclared randomized DEGEN3 sessions (32 AB + 32 AA pairs each, CPU 4, cap 100, seeds 1993/2017) preserve canonical results and give medians -2.13/-2.04%, conservative intervals [-2.50%, -1.96%] / [-2.55%, -1.27%]. The first passes; the second fails noise separation because its AA interval is [-2.27%, +0.60%]. Stopped as timing-inconclusive without extending/repeating the screen or running full gates. Source, accepted CLI and PGO metadata restored. Candidate `/tmp/rose-scatter-candidate{,.rs,-build.json}`; build `/tmp/rose-scatter-build.log`; raw `target/perf-scatter/{fast.txt,transfer-pairs.json,target-*,randomized-degen3-*}`. Accepted-release BORE3D-first instruction profiles and numeric-factor disassembly: `target/perf-factor-current/` (CLI SHA256 `71c94ad078d3808c4b93195f1a576a11eb46533826ce426d88d5fdb514312820`). DEGEN3's 5709 samples contain one throttle/unthrottle pair and no lost records, so its self shares are diagnostic only; disassembly identifies repeated mapped scatter indexing but does not establish removable time. This tests local row/slice hoisting, distinct from the rejected cached-child-map trial below.

- Batch four copy-based DGEMM outputs when all four use sequential Chain accumulation, retaining alpha's `1.0 * A`, per-output FMA order and final C addition; other modes keep their kernels. The current grouped implementation is independently frozen before editing, and both frozen/current and earlier scalar DGEMM tests pass, as does fast parity 940/940. Current accepted release instruction sampling (BORE3D first, eight DEGEN3 solves capped at 100, 5836 samples at period 4M) locates 21.01% self samples in DGEMM and 0.79% in the allocating row-tile selector; these are batch self fractions, not inclusive/removable bounds. BORE3D instructions +0.05% / no benefit; DEGEN3 -1.22%, GREENBEA +0.01%. Three grouped timing pairs remain within spread or cross bands. Two predeclared DEGEN3 randomized sessions (32 pairs each AB/AA, CPU 4, cap 100, seeds 1379/1453) preserve canonical full results and give medians -1.73/-1.29%, conservative intervals [-3.06%, -1.10%] / [-2.54%, -0.71%], with sign p 9.65e-6 / .000268. The first passes; the second fails the noise-separation criterion because its AA interval is [-1.70%, +1.10%]. Stopped as timing-inconclusive under the declared rule, without repeating until pass or running full gates. Source, accepted CLI and PGO metadata restored. Candidate `/tmp/rose-copy4-candidate{,.rs,-build.json}`; logs `/tmp/rose-copy4-{build,tests}.log`; raw `target/perf-copy4/{transfer-pairs.json,target-*,randomized-degen3-*}`. Current release profile runner/manifests/raw data/report: `target/perf-rowtile/`. The small copy-path effect does not establish a profitable multi-mode or multi-column tile.

- Eliminate the ordinary quad RHS twin: source audit finds its only consumer is recompute_primal while the main factor has zero updates. Run the unchanged force-dense RHS routine on the main factor with spike, WORK/parts, result flags/count, F2, L snapshots and ordering scratch moved aside and restored afterward. Keep the original two-build/two-factor path for numeric-construction diagnostics and repeat PGATE records as before. An independent twin-versus-borrowed test checks complete f64/F80 result bits and retained state after prior FTRAN, needed-FTRAN and BTRAN history, empty/extreme/NaN inputs, repeated calls and seeded snapshots; it and fast parity 940/940 pass. Eight CPU-4 AB/BA pairs per quad fixture preserve complete result objects; instructions fx_quad_obj_5371 -2.15%, Quad_2717 -1.17%, median cycles -0.22/-0.30%. A predeclared randomized Quad screen (32 pairs each AB/AA, seed 631, cap 17000) verifies complete canonical result/exit identity: median +0.47%, interval [-0.77%, +1.73%], 14/32 improving pairs (p .811), AA interval [-1.04%, +0.94%]. No gain; source, accepted CLI and PGO metadata restored before full benchmarking. The prototype establishes numerical feasibility of isolated borrowing for tested cases, not its profitability or full parity. Candidate source/binary/build `/tmp/rose-rhstwin-candidate/`; logs `/tmp/rose-rhstwin-{build,tests}.log`; raw `target/perf-rhstwin/{transfer-pairs.json,randomized-quad-631.json,randomized-quad-631.raw/}`. The existing fresh factor cache still requires an independent main copy before mutation.

- Share immutable L columns, L transpose, row permutation and L levels between cloned RefStruct factors in one Rc-owned bundle. Source audit confirms these fields do not change during production updates; U, diagonals, etas, solve caches, counters and snapshots remain independently owned. Existing test-only factor edits use copy-on-write. All 76 library release tests and fast parity 940/940 pass. Eight CPU-4 AB/BA pairs per quad fixture preserve complete no-dump results: instructions fx_quad_obj_5371 -0.62%, Quad_2717 -0.31%; median cycles -1.63/+2.58%. A predeclared randomized Quad_2717 screen (32 pairs each AB/AA, seed 491, cap 17000) checks complete canonical result/exit identity and reports median -0.84%, conservative interval [-3.21%, +0.50%], 20/32 improving pairs (p .108), AA interval [-1.10%, +2.11%]. Inconclusive; no repeat-until-pass, full benchmark or speedup claim. Source, accepted CLI and PGO metadata restored. The prototype proves the immutable/mutable separation for these fields but does not remove mutable-U cloning or construction. Candidate `/tmp/rose-sharefactor-candidate{,.rs,-build.json}`; logs `/tmp/rose-sharefactor-{build,tests}.log`; raw `target/perf-sharefactor/{transfer-pairs.json,randomized-quad-491.json,randomized-quad-491.raw/}`. Rc Send-compatibility would also need addressing before retaining this prototype.

- Replace repeated primal structural-nonzero scans with the existing immutable `offs[n]` prefix total (with an empty/no-slack fallback); arithmetic and WORK formulas are unchanged. All 76 library release tests and fast parity 940/940 pass. Eight CPU-4 ADLITTLE/capped-DEGEN2 primal pairs preserve complete no-dump results; instructions -0.14/-0.38%. Three alternating CPU-4 five-run full DEGEN2/PEROLD-primal/PILOT4 pairs remove about 0.99/1.05/0.77% instructions; each has one timing gain beyond spread, others are within spread or cross bands. Three CPU-2 confirmations remain inconclusive: DEGEN2 -5.72/-3.79/-0.61%, PEROLD +23.73/-3.28/-1.90%, PILOT4 +26.76/+3.14/-2.95%; the first PEROLD/PILOT4 pairs cross bands and the rest stay within spread. No repeatable confirmed gain; stopped before full benchmarking and restored source, accepted CLI and PGO metadata. Candidate `/tmp/rose-nnzcache-candidate{,.rs,-build.json}`; logs `/tmp/rose-nnzcache-{build,tests}.log`; raw `target/perf-nnzcache/{transfer-pairs.json,target-*,confirm-*}`. This is a timing-inconclusive work reduction, not evidence that the scans are necessary.

- Share immutable split LU input columns across basis builds and factor-cache entries: lazily create row/value Rc slices once per Internal, retain incoming-entry order, and generalize LU/pair/fresh/WORK input views to accept either owned Vecs or shared slices. This also makes the pair-pass singleton's defensive row clone shallow. No numerical kernel changes. All 76 library release tests and fast parity 940/940 pass; a fresh matching PGO build passes fast parity again. Eight CPU-4 ADLITTLE/capped-DEGEN2 AB/BA pairs preserve complete no-dump results; instructions +0.35/-0.35%. Three alternating CPU-4 five-run PILOT4/PEROLD-dual/Quad_2717 pairs save only 0.07/0.32/0.25% instructions. PEROLD solve changes +6.74/-0.57/-36.91% (first regresses beyond spread, last crosses host bands); PILOT4 +10.03/-7.06/-37.06% and Quad +3.49/+2.73/-46.45% show no repeatable gain (third pairs cross bands). No material target instruction reduction or acceptable timing result; source, accepted CLI and PGO metadata restored before full benchmarking. This tests shared split-column ownership, not reuse of numerical factor production. Candidate source/binary/build `/tmp/rose-sharedcols-candidate/`; build/tests `/tmp/rose-sharedcols-{build,tests,pgo}.log`; raw `target/perf-sharedcols/{fast.txt,transfer-pairs.json,target-*-{before,after}.txt}`. Input cache allocation costs include all model columns once; public input signatures were broadened only in this rejected prototype.

- Combine native packed BTRAN multiply/update assembly blocks, retaining the ten-byte intermediate product store and exceptional-product/result fallback. A CPU-4 cycle profile of the preceding packed-threshold candidate (12 bounded Quad_2717 solves, validated statuses/iterations/refactors) places 8.52% self samples in the packed kernel; sampled stack moves after the product store motivate eliminating intermediate Option/packed-result transfers. This is sampled attribution, not proof of a specific hardware stall. Primitive differential tests compare the combined blocks with both separate native operations and software F80 on 4096 random/full-low-bit/extreme/subnormal cases; packed threshold and frozen complete-BTRAN tests pass, as does fast parity 940/940. Eight CPU-4 AB/BA pairs per fixture preserve complete no-dump results. Instructions fx_quad_obj_5371 -2.25%, Quad_2717 -2.92%, but median cycles +2.94/+1.88%, solve medians +3.24/+0.71%; cycles rise in 7/8 and 6/8 pairs respectively. No target gain; stopped before full benchmarking and restored source, accepted CLI and PGO metadata. The combined blocks still store and reload the product to preserve its rounding. Candidate `/tmp/rose-packedfuse{,-f80.rs,-ftran.rs,-real.rs,-build.json}`; logs `/tmp/rose-packedfuse-{build,tests}.log`; raw `target/perf-packedfuse/transfer-pairs.json`. Preceding-candidate cycle recording, annotation and validating runner: `target/perf-packedbt-audit/{quad.data,kernel.asm,profile.py,quad_2717-manifest.json}`. No full-parity or solver-speedup claim.

- Keep a packed x87 workspace across dense U-transpose, reverse etas and dense L-transpose, with per-operation 80-bit rounding and nonfinite fallback. The sparse U branch stays on software. An independent frozen complete BTRAN oracle checks outputs, WORK, row-cache and density flags on f64/F80 random/extreme/duplicate/sparse/dense/direct/eta cases; the test and fast parity 940/940 pass. Eight CPU-4 AB/BA pairs per fixture preserve complete no-dump results; first variant instructions fx_quad_obj_5371 -1.93%, quad_2717 -2.46%, median cycles +3.47/+2.68%. Two CPU-2 five-run confirmations show no repeatable gain: Quad changes -42.89/+5.61% (first crosses bands, second within spread), fx +3.09/+2.63% within spread. A second variant replaces threshold unpacking with direct canonical packed magnitude comparisons; threshold-adjacent low-bit/subnormal/extreme/sign tests, complete-BTRAN differential and fast parity pass. Its eight CPU-4 pairs preserve complete results, instructions -2.21/-2.86%, median cycles +6.71/+7.64% and solve medians +6.81/+7.87%. No qualifying speedup; both variants restored before full benchmarking. This tests dense U/eta/L together, not sparse U or a cache of packed coefficients. Candidates `/tmp/rose-packedbt-v{1,2}{,-ftran.rs,-f80.rs,-real.rs,-build.json}`; build/test logs `/tmp/rose-packedbt-{build,tests}.log`, `/tmp/rose-packedbt-v2-{build,tests}.log`; raw `target/perf-packedbt{,-v2}/{transfer-pairs.json,confirm-*-{before,after}.txt}` (confirmations only in the first directory). Accepted CLI and PGO metadata restored.

- Native x87 dense BTRAN L-transpose pass through Real dispatch, with untouched software fallback and guarded control word/nonfinite input, coefficient, product and result handling. The entire parent BTRAN is frozen as an independent oracle; two release differential tests pass (random full representations, f64/F80, empty/singleton/duplicate inputs, sparse/dense/direct modes, row etas, WORK, row-cache/density flags, subnormals and overflow fallback), as does fast parity 940/940. Eight CPU-4 AB/BA pairs per quad fixture preserve full no-dump results; instructions fx_quad_obj_5371 -1.11%, quad_2717 -1.41%, median cycles +1.18/+0.93%. Three alternating CPU-4 five-run quad timing changes +2.11/-10.33/-1.95% all stay within spread; fx changes +2.68/-28.16/+8.44%, with the middle crossing speed bands. No repeatable target gain; restored source, accepted CLI and PGO metadata before full benchmarking. This tests the dense L branch only, not sparse/partial L or U-transpose/etas. Candidate `/tmp/rose-nativebt-candidate{,-ftran.rs,-f80.rs,-real.rs,-build.json}`; logs `/tmp/rose-nativebt-{build,tests}.log`; raw `target/perf-nativebt/{transfer-pairs.json,target-*-{before,after}.txt}`. The standalone finite native probe remains available; its 3.1x isolated result does not predict solver benefit.

- Precompute barrier assembly column ranges, triple-loop split points and local row indices once in CholSolver's Plan. The existing public assembly kernel remains the independent oracle/fallback. Differential tests on AFIRO/BLEND/BORE3D fronts, both precision modes, scaled/unscaled columns and extreme/signed-zero/NaN inputs pass bit for bit; fast parity 940/940 passes. Eight CPU-4 BORE3D AB/BA pairs preserve complete no-dump results with instructions -0.11% and no timing gain. Three alternating CPU-4 five-run DEGEN3 pairs have instructions +0.19%, solve changes +0.18/-4.03/-0.04%, all within spread; the middle pair crosses speed bands. No qualifying gain; source, accepted CLI and PGO metadata restored. This tests cached assembly ranges/indices, not fully compiled numeric contribution destinations or block schedules. Candidate `/tmp/rose-asmplan-candidate{.rs,,-build.json}`, logs `/tmp/rose-asmplan-{build,tests}.log`; raw `target/perf-asmplan/{transfer-pairs.json,degen3-*-{before,after}.txt}`. No full benchmark or full parity claim.

- Precompute child-to-parent row maps and child column ranges in the barrier symbolic Plan: numeric loops, update order, single-mode rounding and arithmetic remain unchanged. Width-one updates stop allocating a mapping on each factorization; wider updates use direct mapped indices. Fast parity 940/940 passes. Eight CPU-4 BORE3D AB/BA pairs preserve complete no-dump results, instructions -0.26%; DEGEN3 instructions -0.77%. Three alternating CPU-4 five-run DEGEN3 changes +0.29/-3.85/-4.18% include two gains beyond spread. Two CPU-2 confirmations -0.22/-16.20% stay within spread; the latter crosses speed bands. No repeatable confirmed timing gain; session stopped as inconclusive, source/accepted CLI/PGO metadata restored. This tests only child ranges and row mappings, not assembly ranges or rectangular block schedules. Candidate `/tmp/rose-childplan-candidate{.rs,,-build.json}`, build `/tmp/rose-childplan-build.log`; raw `target/perf-childplan/{fast.txt,transfer-pairs.json,degen3-*,confirm-*}`. No full benchmark or commit-gate claim.

- Borrow contiguous child panels from packed factor storage instead of allocating/copying them before DGEMM (analysis idea 13): fast parity 940/940 passes and DEGEN3 barrier instructions fall about 0.9%. Three alternating CPU-4 pairs cross speed bands; CPU-2 confirmation median solve medians 0.405972 -> 0.396088 s remain within aggregate 0.096096/0.022184 s spreads. The final stable baseline/candidate pair is 0.397590 -> 0.396078 s, within 0.003082/0.006880 s spreads. BORE3D also has no gain beyond spread. Reverted (`target/perf-childpanel-*`); instruction savings alone do not meet acceptance.

- Reuse child-update temporary/panel/product buffers within each factorization (analysis idea 13): fast parity 940/940 passes; DEGEN3 barrier solve 0.393035 -> 0.387227 s is within 0.010109/0.006277 s spreads and instructions rise about 0.16%. BORE3D timing is also within spread with effectively unchanged instructions. Reverted (`target/perf-childscratch-*`).

- Skip unused front-dump configuration parsing (analysis idea 16): fast parity 940/940 and bounded BORE3D front-dump comparisons pass. Hoisting selection/entries to factorization entry raises DEGEN3 barrier instructions 0.60%; guarding only entry parsing for the selected supernode raises them about 0.58%. Three alternating CPU-2 pairs for the latter give median solve medians 0.398125 -> 0.405277 s across host speed bands; neither meets acceptance. Both reverted (`target/perf-frontcfg{,2}-*`).

- Combine FTRAN touched/emitted markers into one three-state byte array: frozen f64/F80 random/extreme/captured-input tests and fast parity 940/940 pass. DEGEN3 primal-release solve 0.651511 -> 0.646340 s improves 0.79% beyond 0.004325/0.004931 s spreads, but instructions rise 0.17% (DEGEN2 primal +0.08%), failing the below-3% instruction-reduction acceptance rule. Reverted (`target/perf-ftmarks-target-{before,after}.txt`). FTRAN annotation `target/perf-degen3r-thin-ftran.txt` identifies output marker checks, scatter indexing and divisions as remaining costs.

- Avoid materializing the dense FTRAN tail-index vector: direct slice traversal through a shared solve closure passes frozen f64/F80 random/extreme/captured-input tests and fast parity 940/940, but DEGEN3 primal-release instructions rise 1.29% and solve 0.648145 -> 0.656844 s exceeds both 0.001714/0.002821 s spreads. A single chained iterator passes fast parity but still adds 0.25% instructions and has no timing gain (1.019540 -> 1.024090 s within spread). Both reverted (`target/perf-tailview{,2}-target-{before,after}.txt`).

- Cache primal signature bound flags/widths (analysis idea 8): fast parity 940/940 passes. Interleaved `(kind,width)` cache reduces DEGEN3 primal-release instructions 1.76%, but three alternating CPU-2 pairs give median solve medians 1.005068 -> 1.001979 s within 0.042525/0.061157 s aggregate spreads. Split byte flags/width arrays reduce instructions 2.54%, but confirmation medians 1.003885 -> 1.313467 s show no gain and cross speed bands. Both reverted (`target/perf-sigbounds{,2}-{target,confirm}-{before,after}.txt`).

- Skip auxiliary pivot selection and row-count setup when numeric LU receives prescribed rows (analysis idea 7), retaining max-magnitude/floor and supplied-pivot zero checks: fast parity 940/940 passes, PILOT4 instructions -0.11%, DEGEN2 primal -0.22%. Three alternating CPU-2 PILOT4 pairs give median solve medians 0.089099 -> 0.090044 s with no gain (aggregate spreads 0.007874/0.006202 s). Reverted (`target/perf-forcedpivot-{target,confirm}-{before,after}.txt`).

- Maintain B-list reference counts to replace factor-update full-tail membership scans (analysis idea 4): fast parity 940/940 passes, PEROLD instructions -0.59%, DEGEN2 dual -0.67%. Three alternating CPU-4 pairs give median solve medians 0.146994 -> 0.146334 s within aggregate 0.001341/0.003619 s spreads; CPU-2 rechecks cross speed bands without a gain beyond spread. Reverted (`target/perf-brefs-{target,confirm,recheck}-{before,after}.txt`). This tests B counts only; maintained reverse A memberships remain untested.

- Index nonbasic row-list membership (analysis idea 5), retaining linear first-match lookup for duplicate/handoff rows: fast parity 940/940 passes for both variants. Per-row hash maps increase PEROLD instructions 0.37%, DEGEN2 dual 1.55%; immutable matrix-entry IDs with array positions still increase PEROLD 0.05% and DEGEN2 0.29%. Array target PEROLD solve 0.232249 -> 0.229859 s is within 0.013778/0.010542 s spreads; DEGEN2 is slower. Both reverted (`target/perf-rowindex{,2}-target-{before,after}.txt`). The coefficient searches remain unchanged in these trials.

- Shared dense two-vector U-tail traversal with independent column/tau tolerance and scatter rules (analysis idea 1): frozen f64/F80 random/extreme/captured SC50A differential tests and fast parity 940/940 pass. Hoisting per-vector scatter tests out of the coefficient loop reduces PEROLD instructions 0.55%, but three alternating CPU-2 pairs give median solve medians 0.146587 -> 0.147595 s (aggregate spreads 0.001225/0.001968 s), with no timing gain; DEGEN2 target solve also increases. Reverted (`target/perf-pairtail{,3}-target-{before,after}.txt`, `target/perf-pairtail3-confirm-{before,after}.txt`). This trial shares only the dense U tail; the paired-stage trials above cover other measured paths.

- Cache the immutable L-level maximum for FTRAN (analysis idea 8): fast parity 940/940 passes and instructions fall 0.64% on PEROLD, 1.15% on DEGEN2 dual. PEROLD initial solve 0.148098 -> 0.146727 s is within the 0.001604 s baseline spread; three CPU-4 confirmation pairs give median solve medians 0.148094 -> 0.148166 s. CPU-2 rechecks cross speed bands and show no gain beyond aggregate spreads. Reverted (`target/perf-maxlev-{target,confirm,recheck}-{before,after}.txt`). Timing evidence is inconclusive; the instruction reduction alone does not meet acceptance.

- Marker array for factor-update spike deduplication instead of repeated `contains`: fast parity 940/940 passes, PEROLD instructions -0.23%, but three alternating CPU-2 pairs give median solve medians 0.230458 -> 0.231593 s with no gain (aggregate spreads 0.087231/0.091101 s); final pair also has no improvement. Reverted (`target/perf-spikemark-{target,confirm}-{before,after}.txt`).

- Share sparse eta accounting with the numeric pass only when `first == 0` (analysis idea 2): frozen f64/F80 random/extreme/captured-input tests and fast parity 940/940 pass, but PEROLD instructions fall only 0.11% and solve 0.232151 -> 0.234504 s has no gain (spreads 0.017993/0.004013 s); DEGEN2 instructions are flat. Reverted (`target/perf-etawork-target-{before,after}.txt`). Guarded positive-zero-prefix sharing is also rejected above; the accounting input omits prior eta updates, so sharing after changed prefixes still needs a rule for sub-tolerance residues.

- Reuse dual pivot-row accumulator and marker buffers with full clearing (analysis idea 6): fast parity 940/940 passes; initial PEROLD instructions -0.47%, but three alternating CPU-2 pairs have effectively flat instructions and no timing gain (median of solve medians 0.236570 -> 0.241332 s; aggregate spreads 0.012539/0.082480 s). Reverted (`target/perf-rowworkspace-{target,confirm}-{before,after}.txt`).

- Column-local dual row-product sums for dense, strictly ascending rho and unique matrix rows: frozen random/extreme/captured ADLITTLE f64/F80 tests and fast parity 940/940 pass, but PEROLD instructions +1.8% with separate presence/value arrays, +2.2% with an Option buffer. The latter regresses solve 0.236082 -> 0.240998 s beyond 0.003144/0.002701 s spreads. DEGEN3 instructions also increase; both variants reverted (`target/perf-colrow{,2}-target-{before,after}.txt`).

- Reserve FTRAN result capacity from the output-order length: fast parity 940/940 passes and instructions fall 0.4–0.6%, but the initial PEROLD timing gain crosses host speed bands; alternating confirmation 0.151656 -> 0.149696 s is within 0.093843/0.090559 s spreads, and DEGEN2 has no repeatable gain. Reverted (`target/perf-ftreserve-{target,recheck}-{before,after}.txt`).

- Borrow pre-sorted columns in dual true_d instead of cloning and sorting: fast parity 940/940 passes and instructions fall 0.4–0.6%, but alternating CPU-2 PEROLD solve 0.231300 -> 0.231688 s has no gain (spreads 0.017962/0.008958 s); DEGEN2 and DEGEN3 stay within wide spreads. Reverted (`target/perf-truecols-{target,recheck}-{before,after}.txt`).

- Split dense dual pivot-row marking from sparse first-touch ordering: fast parity 940/940 passes, PEROLD instructions -0.32%, but alternating CPU-2 solve 0.232569 -> 0.227484 s is within 0.009078/0.006368 s spreads; DEGEN2 has no repeatable timing gain. Reverted (`target/perf-rowmark-{target,recheck}-{before,after}.txt`).

- Separate immutable L-transpose index and coefficient arrays with contiguous row offsets: fast parity 940/940 passes, but PILOT4 instructions +2.7%, DEGEN2 +3.5%, PEROLD primal +4.5%, and all timing changes stay within spreads. Reverted (`target/perf-ltsoa-target-{before,after}.txt`).

- Four-entry ordered scalar dense BTRAN scatter (reference U-transpose 0x12e25e, L-transpose 0x130197): frozen random/extreme/duplicate/captured SC50A f64/F80 tests and fast parity 940/940 pass, but PILOT4 instructions +0.3%, DEGEN2 +0.8%, PEROLD primal +0.6%; no timing gain beyond spreads. Reverted (`target/perf-btunroll-target-{before,after}.txt`).

- Guard the sparse L rollback clone by m >= 250 and Kahn reached-count upper bound: fast parity 940/940 passes, but instructions are effectively flat (PILOT4 -0.01%, Quad_2717 +0.01%); DEGEN2 0.068212 -> 0.064283 s is within the candidate 0.004463 s spread. Reverted (`target/perf-lclone-target-{before,after}.txt`).

- Transfer dense L-stage x directly into U with zero canonicalization: frozen f64/F80 U tests and fast parity 940/940 pass. Combined collection/normalization lowers DEGEN3 dual instructions 2.1%, but CPU-2 solve 0.486473 -> 0.493595 s has no gain (spreads 0.012765/0.036416 s). DEGEN2 initial -2.6% does not repeat beyond spread; PILOT4 instructions rise 0.2%. A separate zero-normalization pass was worse. Reverted (`target/perf-lbuffer{,2}-target-{before,after}.txt`, `target/perf-lbuffer2-recheck-{before,after}.txt`).

- Pack immutable L-transpose rows into contiguous storage: fast parity 940/940 passes, but PILOT4 solve 0.095785 -> 0.096071 s is within 0.002578/0.006817 s spreads and instructions rise 0.3%; DEGEN2 instructions rise 0.15% with no timing gain. Reverted (`target/perf-ltcsr-target-{before,after}.txt`).

- Borrow Markowitz candidate chains and rotate buckets in place: fast parity 940/940 passes and PEROLD instructions -0.5%, but alternating CPU-2 solve 0.237249 -> 0.236897 s is within 0.006498/0.011017 s spreads; DEGEN3 primal also stays within wide spreads. Reverted (`target/perf-searchborrow-{target,recheck}-{before,after}.txt`).

- Cache first-match pivot-row/column coefficients in a stamped lookup: semantic tests and fast parity 940/940 pass, PEROLD instructions -0.6%, but solve 0.151853 -> 0.151628 s is within 0.004297/0.000431 s spreads; DEGEN2 has a small timing increase with flat instructions and DEGEN3 primal has no gain. Reverted (`target/perf-pivotvalues-target-{before,after}.txt`).

- Lazy reuse of the extended mode-2 coefficient `(wq_q / a_rq) / a_rq` per primal row call. Both immutable operands, division order, first-use location and each column's multiplication/addition order are retained; primitive kernels and WORK sites remain unchanged. Complete frozen row body stays source-identical. Matching fixed39 PGO passes existing whole-row and new repeated-column differential tests, full parity 987/987, all 37 release groups (253 tests), formatting and all 23 complete benchmark result/exit pairs. Fx instructions fall 3.11%; target screens 13217/13289 (32 AB + 32 AA each, CPU 4, cap 5500) both qualify: median -3.620/-4.568%, intervals [-5.541%, -2.697%] / [-5.132%, -4.079%], AA [-2.193%, +1.256%] / [-0.108%, +1.543%], 30/32 and 32/32 improving, both orders improve. Full sweep flags DEGEN3 primal-r +0.76%, PEROLD dual-r +59.10% and GREENBEA dual-c +16.09%, all with unchanged instructions. Fixed three-pair follow-ups are mixed; the predeclared reversed-role controls establish no regression. However DEGEN3's AA control is biased: median +0.675%, interval [+0.077%, +2.751%], p 0.0201. Its predecessor AB median is -0.330%, interval [-0.801%, -0.138%], sign p 0.01003 and no conservative noise separation. PEROLD/GREENBEA controls have healthy AA and no established slowdown. The final control had already completed when stop was requested; all three complete controls are retained. Qualification is inconclusive under the healthy-AA rule; no screen/control is extended or repeated. Source, accepted CLI/build and opaque-input test fix are restored exactly; no optimization is accepted. Source/proof, frozen oracle/new test, diagnostic provenance, matching build/assembly, both positive target screens and every full gate/result/control: `target/perf-devex-mode2/`.

- Exact normal F80 cancellation for opposite effective signs and equal exponents, or adjacent exponents with magnitude ratio below two: a 64-bit difference and normalization replace the general 128-bit align/pack only when no rounding is required. The complete original add/sub primitive is frozen independently; reconstruction reverses only the new branch. Matching fixed39 PGO passes all eleven F80 tests (two million random/full-representation/close-magnitude comparisons, captured SC50A factors, signs/zeros/nonfinites/boundaries and existing x87/conversion tests), fast parity 940/940, and eight ADLITTLE/fx plus three Quad_2717 complete result/exit pairs. Quad instructions -2.11%, but cycles rise in all three pairs (+3.63/+2.05/+1.70%) and branch misses +13.69/+8.33/+5.11%. Fx instructions -1.44%, cycles rise in all eight pairs (+0.61 through +9.50%). These whole-CLI counters do not attribute all extra cost to branch prediction, but establish no useful target gain. Stopped before predeclared timing seeds 12217/12289 and full gates. Source, accepted CLI/build and opaque-input test fix are restored exactly. Candidate/frozen source, algebra/rounding proof, matching build, tests, all outputs/counters and restoration: `target/perf-f80-cancel/`.

- Normal F80 division normalizes the quotient before one exact 128-by-64 DIV, then rounds from its exact remainder; boundary/exceptional results retain the complete original two-step expansion/packer. Complete original division is frozen independently; source reconstruction removes only the normal branch. Final matching-PGO checks pass all eleven F80 tests (one million random divisions, captured SC50A values, boundary/nonfinite/signed-zero cases, 100,000 helper comparisons and existing x87/conversion tests), plus fast parity 940/940. Eight bounded ADLITTLE/fx and three Quad_2717 pairs preserve complete no-dump result/exit objects. Quad instructions -2.40%, cycles -4.35%; fx instructions -3.18%, cycles -5.88%. Native normal code has one DIV; fallback retains both helper calls. The first predeclared Quad screen (seed 11411, 32 AB + 32 AA, CPU 4, cap 17000) is inconclusive: median -4.060%, interval [-5.100%, -1.557%], AA [-0.251%, +1.760%], 28/32 improving, sign p9.65e-6, both orders improve, healthy AA. Conservative noise separation fails. The reserved confirmation 11483 was accidentally started after a shell assertion because the shell lacked stop-on-error; it was interrupted (exit 130), and partial schedule/raw data are retained without a confirmation decision. No extension or full gates. Source, qualified CLI/build and separately committed opaque-input test fix are restored exactly. Initial failed mode-test evidence, final source/proof/build/tests/assembly/counters, complete failed screen and partial confirmation: `target/perf-f80-div-normal/`.

- Fixed initialized slots for dual pivot-row first-touch order remove list growth from the coefficient scatter, preserving ordered indices, numeric expressions and WORK exactly. Reconstruction reverses only storage/count/truncation changes; producer/update/handoff bounds prove every structural index is below n and first-touch count is at most n. Matching fixed39 PGO and fast940 parity pass; eight ADLITTLE/DEGEN2 and three PEROLD pairs preserve complete no-dump result/exit objects. Instructions fall0.48/0.63/0.71%; PEROLD cycles fall1.20%. Native scatter keeps the row multiplier in XMM0. The predeclared PEROLD screen (seed10007, 32AB+32AA, CPU4, cap1300) gives median -1.056%, interval [-1.585%, -0.680%], AA [-0.761%, -0.050%], 30/32 improving, sign p1.23e-7, both orders improve. AA is biased (p0.0201) and conservative noise separation fails: inconclusive. Stop without seed10069, extension or full gates. Source/qualified CLI/build restored exactly; complete source/proof/profile/assembly/counters/schedule/results and restoration: `target/perf-pivot-fixedorder/`.

- Schedule all sparse-L successor queue insertions for each column before its numeric scatter, keeping identical queue insertion/active/remaining updates and the same numeric operations/order. The complete frozen L oracle passes all three focused release tests (f64/F80 random/extreme/captured inputs, output order/value bits, density and WORK estimate), and fast parity passes940/940. Eight bounded ADLITTLE/DEGEN2 pairs and three PEROLD pairs preserve complete no-dump result/exit objects. PEROLD instructions +0.108%, cycles +0.184%, median solve +0.800%; DEGEN2 instructions +0.120%, cycles +0.325%. No useful target counter gain; stopped before predeclared timing seeds9601/9667 or full gates. Source, accepted CLI and PGO metadata restored exactly. Candidate source/build, original source, frozen-source hash, oracle logs, assembly, complete results/counters and setup-build provenance: `target/perf-lschedule/`.

The first outlined four-register SSE2 eta-dot trial stopped before timing at
n=4, peel0 with two distinct NaN payloads in its unchanged scalar fallback.
The same nonopaque original-versus-frozen test also fails on unchanged accepted
source at exactly that case. That observation is not proof of a SIMD arithmetic
mismatch. Making both function pointers opaque preserves every test input and
comparison and makes the original/frozen test pass before any kernel rewrite.
The first candidate remains unqualified; its evidence is in
`target/perf-eta-sse2/`.

A separate long-only eta trial places native dispatch after the unchanged short
and extended branches. Four SSE2 registers retain the original eight scalar
chains and reduction tree; a nonfinite native result falls through the complete
original long body. Source reversal reproduces the original dot exactly, and
all WORK sites remain unchanged. Matching fixed39 PGO passes the validated
opaque original oracle on random/full-bit/extreme/two-NaN-payload/captured
inputs and all16 MXCSR rounding/flush modes. All eight FTRAN/L tests and fast
parity940 pass; eight bounded ADLITTLE/DEGEN2 and three PEROLD canonical
result/exit pairs agree. Native assembly has register accumulators and packed
multiply/add, without reassociation or FMA. Instructions fall0.11/0.49/1.83%;
PEROLD median cycles fall0.73%. The predeclared PEROLD screen (seed15431,
32AB+32AA, CPU4, cap1300) is inconclusive: median -1.506%, interval
[-2.550%, -0.759%], AA [-0.882%, +1.017%], 24/32 improving, sign p0.00350,
both orders improve and AA is healthy. Conservative noise separation fails.
Stop without seed15503, extension or full gates. Solver sources, accepted CLI
and build metadata are restored exactly; restored fast parity940 passes.
The new test/oracle snapshots are retained only as evidence, with no production
or committed test change. Baseline failure/pass, source proof, matching build,
all differential tests, assembly, counters, complete schedule/results and
restoration: `target/perf-eta-long/`.

- Pair eight-element eta chunks to eliminate per-lane checks: frozen f64/F80 tests and fast parity 940/940 pass, dot code shrinks 1037 -> 879 bytes and PEROLD instructions -1.5%, but CPU-4 solve medians are 0.150879/0.150880 s; alternating CPU-2 changes stay inside wide spreads. Reverted (`target/perf-dotchunks-{target,recheck}-{before,after}.txt`).

- Explicit eight-lane eta blocks: frozen random/extreme/captured tests, full parity 987/987 and release tests pass; PEROLD instructions -0.9%, but the full CPU-4 solve 0.151302 -> 0.151144 s stays within 0.001056/0.000552 s spreads, and confirmation has no gain. PILOT4 confirmation regresses 0.094533 -> 0.097417 s (+3.1%; spreads 0.002414/0.001277 s). Reverted (`target/perf-dotblocks-{target,recheck,full,confirm}-{before,after}.txt`).

- Read eta inputs directly by index instead of gathering: frozen random/extreme/captured f64/F80 tests and fast parity 940/940 pass, PEROLD instructions -1.0%, but alternating CPU-2 solve 0.230032 -> 0.232767 s has no gain (spreads 0.006889/0.003601 s); DEGEN3 also stays within spread. Reverted (`target/perf-indexdot-{target,recheck}-{before,after}.txt`).

- Elide pivot-row scatter bounds checks using private row-list indices: fast parity 940/940 passes and PEROLD instructions -0.9%, but alternating CPU-2 solve 0.234652 -> 0.232350 s is within 0.002445/0.004874 s spreads; DEGEN3 has no repeatable timing gain. Reverted (`target/perf-rowindex-{target,recheck}-{before,after}.txt`).

- Hold immutable kind/diag/U slices across the needed-vector scatter: fast parity 940/940 passes and DEGEN3 dual-controlled instructions -0.8%, but solve 0.461177 -> 0.488074 s has no gain (spreads 0.041088/0.008475 s); DEGEN2/PEROLD stay within spread. Reverted (`target/perf-tauviews-target-{before,after}.txt`).

- Zip L-stage row/value slices to remove coefficient index checks: fast parity 940/940 passes and DEGEN3 dual-controlled instructions -1.7%, but alternating CPU-2 solve 0.612644 -> 0.520240 s is within 0.144281/0.146594 s spreads; DEGEN2 and the first CPU-4 batch also stay within spread. Reverted (`target/perf-lzip-{target,recheck}-{before,after}.txt`).

- Query row lists for former B-column references instead of rebuilding all flags: fast parity 940/940 passes, PILOT4 instructions -0.3%, but alternating CPU-2 solve 0.095259 -> 0.097082 s has no gain (spreads 0.000879/0.004516 s). Reverted (`target/perf-rowrefs-{target,recheck}-{before,after}.txt`).

- Skip capacity counts when the unfiltered U-row length is below capacity: fast parity 940/940 passes, PILOT4 instructions -0.5%, but alternating CPU-2 solve 0.094382 -> 0.099878 s has no gain (spreads 0.002361/0.010112 s). Reverted (`target/perf-rowcapguard-{target,recheck}-{before,after}.txt`).

- Reserve the dense L-stage output capacity once: fast parity 940/940 passes and DEGEN3 dual instructions -0.8%, but CPU-4 solve 0.345439 -> 0.340777 s is within the old 0.007061 s spread; CPU-2 has no timing gain. Reverted (`target/perf-lout-{target,recheck}-{before,after}.txt`).

- Read growth-test diagonal without constructing unused diagnostic L columns: fast parity 940/940 passes and DEGEN2 instructions -0.4%, but CPU-4 solve medians are identical at 0.045557 s (spreads 0.000079/0.000667 s); CPU-2 timing also stays within spread. Reverted (`target/perf-diagread-{target,recheck}-{before,after}.txt`).

- Remove dense BTRAN L-transpose bounds checks using immutable-index invariants: fast parity 940/940 passes and DEGEN2 instructions -0.6%, but alternating CPU-4 solve 0.070389 -> 0.070534 s has no gain (spreads 0.001005/0.000560 s). Reverted (`target/perf-ltindex-recheck-{before,after}.txt`).

- Force-inline existing F80 add/sub: fast parity 940/940 passes, but Quad_2717 instructions +2.1% and solve 0.443852 -> 0.448784 s has no gain; fx_quad instructions +1.1%, solve +3.2%. Reverted (`target/perf-inlineadd-target-{before,after}.txt`).

- Return the larger finite F80 addend when exponent gap >= 66: frozen boundary/random/captured differential tests and fast parity 940/940 pass, but Quad_2717 instructions +1.4% and solve 0.454318 -> 0.475658 s (+4.7%, beyond both spreads). Reverted (`target/perf-distantadd-target-{before,after}.txt`).

- Fast F80 multiplication by finite powers of two for normal results: frozen/random/captured differential tests and fast parity 940/940 pass, but Quad_2717 instructions +0.5% and solve 0.709394 -> 0.712114 s has no gain. Reverted (`target/perf-powmul-target-{before,after}.txt`).

- Cache immutable extended row-wise coefficients on first use: frozen differential and fast parity 940/940 pass, Quad_2717 instructions -1.3%, but solve 0.708071 -> 0.702991 s is within 0.009744/0.035976 s spreads. Reverted (`target/perf-qcoef-target-{before,after}.txt`).

- Mask unused extended cost-change row products: frozen differential and fast parity 940/940 pass, Quad_2717 instructions -1.4%, but alternating CPU-4 solve 0.701696 -> 0.671770 s is within 0.079075/0.140497 s spreads. Reverted (`target/perf-qmask-recheck-{before,after}.txt`).

- Rebuild all elimination WORK row counts by sequential U-transpose scans: fast parity 940/940 passes and PILOT4 instructions -1.2%, but alternating CPU-4 solve 0.111296 -> 0.116587 s has no gain (spreads 0.009052/0.003302 s). Reverted in favor of on-demand counts (`target/perf-rcountrows-recheck-{before,after}.txt`).

- Emit dense FTRAN output directly without ord/done arrays: fast parity 940/940 passes, but PILOT4 solve 0.073238 -> 0.074108 s and instructions 737.9 -> 738.6M; DEGEN2 timing is within spread. Reverted (`target/perf-denseout-target-{before,after}.txt`).

- Move completed FTRAN x/input buffers into F2 after output construction instead of cloning: fast parity 940/940 passes, PILOT4 instructions 737.8 -> 737.3M (-0.06%), solve 0.073656 -> 0.073330 s is within old 0.000944 s spread; DEGEN2 is within spread and PEROLD changes host bands. Reverted (`target/perf-f2move-target.log`).

- Reuse BTRAN dense x/w buffers: fast parity 940/940 passes, PILOT4 instructions 747.0 -> 745.9M (-0.14%), but solve 0.113494 -> 0.116470 s has no gain (spreads 0.003324/0.002445 s); DEGEN2 is within spread and PEROLD changes host bands with little instruction change. Reverted (`target/perf-btbuf-target.log`).

- Precompute eligibility for row_product_work scans longer than n: fast parity 940/940 passes, PEROLD primal instructions -0.4%, solve 1.342457 -> 1.332925 s is within 0.249063/0.547923 s spreads; DEGEN3 primal-controlled has no timing gain and DEGEN2 instructions rise slightly (`target/perf-rowelig-long-target.log`). Reverted.

- Lazy per-column eligibility memoization in row_product_work: fast parity 940/940 passes, PEROLD primal instructions -1.1%, but CPU-2 solve 1.322285 -> 1.320622 s is within 0.013795 s spread; CPU-0 alternating recheck crosses speed bands with no clear aggregate gain (`target/perf-rowelig-{target,recheck}.log`). Reverted.

- Binary lookup in the U transpose for tau terms: fast parity 940/940 passes; DEGEN3 dual instructions 3769.4 -> 3753.4M (-0.4%), solve 0.739033 -> 0.716611 s is within old 0.025890 s spread. DEGEN2/PEROLD instructions rise slightly. Replaced by direct row-map trial (`target/perf-taulookup-target.log`).
- F80 add aligned with three guard bits, plus direct normal-result rounding: frozen old-add/pack differential test and fast parity 940/940 pass, instructions fall 4977.9 -> 4836.8M (-2.8%), but primal_5959 dual CPU-4 solve 0.555599 -> 0.557742 s (spreads 0.006716/0.011524 s). Alignment alone also has no timing gain on CPUs 2/4 despite -1.7% instructions. Reverted both (`target/perf-add-{target,recheck}.log`, `target/perf-add-direct-{target,recheck}.log`).
- Read eliminate work counts directly from maintained rowl lengths instead of rebuilding counts: parity full 987/987 and release tests pass, but mode_switch.dual-r 22.849 -> 21.851 s is smaller than the old 1.849 s spread (instructions 146886.1 -> 140945.4M, -4.0%). Alternating DEGEN2 primal 0.092/0.090 then 0.091/0.091 s; dual 0.068/0.067 s, marginal against spread. Reverted.

- Borrow nucleus elimination row/column chains instead of cloning them: fast parity 940/940 passes, but CPU-2 three-run PEROLD solve 0.130605 -> 0.129989 s (-0.47%) stays within the old 0.001630 s spread despite instructions -0.24%. DEGEN2/DEGEN3 dual-release timing also stays within spread, with instructions only -0.05/-0.03%. Reverted (`target/perf-elimborrow-target-{before,after}.txt`, `target/perf-elimborrow-target-summary.json`).

## Tried, breaks parity

None.
