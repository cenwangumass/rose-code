# Factor production and ownership audit

Ownership audit source: b421f8ab; guarded replay screen: d897aa91;
compact RHS reuse: prxvqvrr; pure reduced-cost reuse: ztrylnzw; pure ranked-product reuse: kpttrrms.
Cost claims below identify their measured scope. Reusing producer U values in
the consumer factor remains unproven;
guarded elimination of overwritten replay arithmetic is feasible on tested
cases but has no established runtime benefit (below).

| Producer | Consumers and lifetime | Consequence |
|---|---|---|
| `lu_order::Structure` | Pivot order, phases, L rows/multipliers/raw entries, diagonal and partial U feed `RefStruct::fresh`; structure and incoming columns remain live for WORK and pricing installation. | Consuming the structure at the current fresh call would invalidate later consumers. Narrow ownership transfers need an installation payload first. |
| `RefStruct::fresh` | Copies diagonal; remaps L destinations into flat forward storage and builds its transpose. Replays elimination for tail U, then overlays Markowitz entries. | Producer L row IDs and consumer step IDs differ. Moving nested producer vectors directly cannot replace the remap. Direct numerical U reuse requires intermediate equality evidence; guarded skipping of overwritten replay values is a separate case below. |
| Primal `reference_order_for` | Returns the accepted growth-check mirror for installation (`simplex.rs`, `growth_factor`). Retried attempts book their own WORK. | Duplicate accepted growth-check construction is already eliminated. Do not count it as remaining savings. |
| Primal `quad_mirrors` | Each compact cache entry retains a pure main-factor starting point and an isolated RHS provider. Main factors deep-clone before mutation; repeated RHS owners share the provider. Numeric-build diagnostics retain their two builds; uncacheable factors keep direct RHS ownership. | Duplicate ordinary construction and repeated compact RHS-factor cloning are eliminated. RHS results, pure WORK estimates and reduced costs have separate exact keys and four-entry bounds; pure ranked products have a 32-entry bound. Mutable main U/counters/scratch remain independent. |
| Auxiliary `Factor` | Constructed after mirror acceptance; used in primal fallback, cleanup, dual deferred-free pricing and stall-objective evaluation. Update counts also participate in controller decisions. Construction failure can restart primal or return singular basis in dual. | Deletion or laziness requires preserving both numerical consumers and the original failure point. Matching pivot rows is insufficient. |
| Incoming basis columns | `reference_order_for` splits immutable `p.fcols` into owned row/value arrays on each cache miss. Retained for structure construction, fresh replay and WORK. | Borrowed column views could remove copying without changing arithmetic; audit all factor interfaces and cache ownership before changing the input type. |

Source anchors: `src/lu_order.rs` Structure; `src/ftran_order.rs` RefStruct and
fresh; `src/simplex.rs` reference_order, quad_mirrors, refactor, reference_order_for;
`src/dual.rs` lu_attempt and refactor. In fresh, singleton arithmetic uses raw
`(x*a)/d` where available, not `x*(a/d)`; tolerance tests and Markowitz overlays
prevent treating the replay as redundant by inspection.

## Cost evidence and limits

Saved `target/perf-native-families/summary.json` profiles predate sparse auxiliary
etas and current primal dispatch. PEROLD factor/update self mass is about 243–246M
instructions of 1050.5M whole-CLI instructions (23%). The session review attributes
38–44M to fresh and 49–51M to auxiliary construction, about 8–9% combined. These
are historical screening bounds, not current solve-only shares. Even deleting both
would leave most of the recorded gap. Allocator costs are separate; do not count
all allocations as factor ownership savings.

The accepted release CLI hash was rechecked as
`e28b651026b477f0e5741e881efce27d640e2bf3f330f85a4f76b0e53c8804fc`.
CPU-4 same-binary runs (12 each, dual release, cap 300, dumps disabled):

| Model | Median solve | Min–max | Result |
|---|---:|---:|---|
| ADLITTLE | 2.201 ms | 1.975–2.361 ms | optimal, 73 iterations, 9 refactorizations, exit 0 |
| DEGEN2 prefix | 12.853 ms | 11.766–13.366 ms | iteration_limit, 300 iterations, 8 refactorizations, exit 3 |

All complete result objects match within each model after removing seconds.
Raw results and exact command/binary provenance: `target/perf-factor-audit/baseline.json`
and adjacent result files. These short runs take under a second together, but their
18%/12% min–max spans relative to median cannot resolve a small speedup. They do
not refresh factor self-cost attribution or establish reference equality.

## Current instruction screen

The accepted CLI was sampled directly (no frame pointers), CPU 4, 128 independent
CLI solves per batch, same dual/release/cap-300 arguments as above. Each result
was checked for expected status, iteration count, refactorization count and exit
code. ADLITTLE used instruction period 200,000 (3544 samples); DEGEN2 used
1,000,000 (8673 samples). Both reports show zero lost samples. These are batch
self-sample fractions, including Python launcher, parsing and process startup;
they are not solve-only percentages or inclusive costs. Throttle-event and
represented-instruction audits remain necessary for quantitative cost bounds.

| Self symbol | ADLITTLE batch | DEGEN2 batch |
|---|---:|---:|
| fresh | 1.27% | 1.40% |
| auxiliary Factor::factor | 1.13% | 2.76% |
| factor_structure_r | 0.68% | 2.06% |
| ftran_with | 2.06% | 7.83% |
| lstage_num | 1.69% | 5.40% |

Nucleus is separately 2.03% on ADLITTLE. Allocator ancestry is not resolved by
these self profiles. The evidence does not justify a broad factor-input rewrite
yet; FTRAN representation has a larger directly visible target on this prefix.
Raw data/reports and validating batch runner: `target/perf-factor-audit/{adlittle,degen2}.data`,
`*-report.txt`, `profile.py`. No reference comparison was performed.

A first FTRAN ownership prototype deferred the F2 installation until output
extraction finished, then moved the dense result and support instead of cloning.
All arithmetic, output ordering and WORK calls stayed in place; the existing
frozen FTRAN implementation was untouched. Fast parity passed 940/940. Small
instruction counts were essentially flat; a sequential PEROLD screen showed
+0.424% instructions and no timing gain beyond spread. The prototype was
rejected and the accepted source/binary/build metadata restored. Details are in
STATE.md's rejected-trial entry. This does not test elimination of the L-to-U
round trip itself.

A second prototype retained dense L storage across U, avoiding the fresh U
allocation and scatter while retaining the ordered spike for accounting and
updates. Sparse results kept the existing expansion. Dense omitted signed zeros
were normalized to positive zero, preserving the old expansion's representation.
Six frozen f64/F80 FTRAN/L tests and fast parity 940/940 passed. Eight CPU-4 small
AB/BA pairs matched complete no-dump results. Instructions rose 0.28% on ADLITTLE
and 0.82% on capped DEGEN2. On PEROLD, three alternating CPU-4 five-run pairs
regressed solve time 1.69/2.46/1.56%, all beyond spread, with instructions +1.61%.
The prototype was rejected and the accepted source/binary/PGO metadata restored.
This rules out transfer plus a separate signed-zero normalization pass; it does
not establish that removing the ordered spike or fusing compaction is safe or
profitable. Evidence and exact candidate artifacts are linked in STATE.md.

A shared-column ownership prototype lazily split immutable model columns once,
then selected Rc row/value slices for each basis and shared them with factor-cache
entries. Generic dereferenceable input views left LU/pair/fresh numerical kernels
unchanged and made one pair-pass defensive row clone shallow. All 76 library
release tests and fast parity 940/940 passed; the fresh PGO build passed fast
parity as well. Small no-dump results matched completely except seconds. On
PILOT4 primal, PEROLD dual and Quad_2717, instructions fell only 0.07/0.32/0.25%.
Three alternating five-run pairs established no repeatable timing gain; large
third-pair improvements crossed host speed bands. The prototype was rejected,
and the accepted source/CLI/build metadata restored. This tests immutable split
input ownership only; it does not establish the cost of numerical producer
reuse or immutable factor-data sharing. Exact evidence is in STATE.md.

A forward-factor sharing prototype grouped immutable L starts/entries, L
transpose, row permutation and L levels behind one shared owner, so fresh quad
twins and cached-factor clones retained these arrays without deep copies. U,
diagonals, etas, counters, snapshots and solve scratch remained independent.
Source audit found no production mutation of the shared fields. All 76 library
release tests and fast parity 940/940 passed. Target instruction savings were
only 0.62% on fx_quad_obj_5371 and 0.31% on Quad_2717. A predeclared randomized
Quad_2717 screen retained all canonical results and timing samples: median
-0.84%, conservative interval [-3.21%, +0.50%], AA interval [-1.10%, +2.11%].
The session stopped as inconclusive and restored the accepted implementation.
This isolates immutable forward storage; mutable-U cloning is still present.
The Rc prototype was not retained, including its change to Send compatibility.
Evidence is in STATE.md and `target/perf-sharefactor/`.

## Compact RHS reuse

The RHS factor has one solver consumer, `recompute_primal`, which reads returned
solve/estimate values and never updates it. Its provider also owns the separate
pure reduced-cost cache below; that cache does not inspect or mutate the factor. Compact factor-cache entries
now retain an isolated `QuadRhs` provider. Sharing it avoids the RHS deep clone
on repeated basis hits; the main factor remains independently owned. Each provider
retains a separate RHS factor plus at most four extended-input/result pairs and
four double-input/scalar WORK estimates. Retained storage increases by one RHS
factor per compact cache entry; the existing four-factor, fewer-than-50-row gate
bounds the scope.

The result key preserves ordered row indices, duplicates, signed zeros and every
extended significand bit. The estimate key preserves the original double RHS
bits, because extended perturbations can make the actual solve's RHS differ.
WORK booking stays at its original call site and order. Provider factors never
receive column updates; their numerical FTRAN/work-estimate routines are unchanged.
NaN, oversized or mismatched solve inputs use the old path. Dumps execute live
numerics so all internal stage records remain available. Larger and diagnostic
builds retain direct RHS ownership.

Independent comparisons use an AFIRO mixed basis, finite random full-low-bit
values, signed zeros, subnormals, nonfinites, repeated/duplicate/empty inputs,
clone histories, independent main updates, bounded eviction and input lengths
0, 1, 2, 27, 49, 99, 100 and 101. Exact output and estimate bits agree; Send
compatibility is retained. Full parity passes 987/987, all 37 release groups pass,
and all 23 complete no-dump benchmark result/exit pairs agree. Two independent
predeclared primal_5959 primal screens improve median solve by 7.31/4.79% and
meet every control rule. Instructions fall 10.90% there, 11.39% on its dual,
3.98% on Quad_2717 and 8.88% on fx_quad_obj_5371. These are within-session
measurements, not a claim that Rose meets all original-reference targets.

The earlier isolated RHS-twin borrowing experiment restored every touched main
factor field, passed its tested numerical cases and fast parity, but had no
established runtime benefit. It did not reuse returned solves or scalar estimates.
Its rejection remains in STATE.md; it is distinct from the accepted provider.
Qualified source/build, raw controls/results, full gates, regression checks and
implementation audit: `target/perf-rhs-memo/final/`.

## Compact reduced-cost reuse

On a fresh main factor, `set_d_from_y` can reuse pure native extended reduced-cost
outputs from the existing shared provider. Each provider belongs to one Solver;
its `p` columns/costs remain immutable. Ordered full extended dual bits, the full
nonbasic mask, exclusion and costfree flag key up to four FIFO entries. At most
500 variables are retained. Other control modes, dumps, NaNs and malformed sizes
fall back. Main factor metadata is neither read from cached outputs nor restored.
The numerical kernel and WORK call sites are unchanged.

The old kernel writes every double entry, but only nonbasic non-excluded extended
entries. Hits preserve current extended values at all other positions, rather
than replaying stale values there. Differential release tests compare the unchanged
kernel across clones, history, low-bit/signed-zero/extreme inputs, mask equivalence,
key changes and eviction. Full parity passes 987/987, all 245 release tests pass,
and all 23 complete no-dump result/exit pairs agree. Two independent controlled
fx_quad_obj_5371 screens improve solve by 6.30/4.10%; instructions fall 1.24%.
The full sweep's four timing flags and their fixed follow-ups/controls are resolved
in STATE.md. This establishes the local gain, not the overall reference goal.
Source/binary/build, kernel/ownership proof and all qualification data:
`target/perf-cost-memo/`.

## Compact ranked-product reuse

Each shared compact provider belongs to one Solver whose problem and row-wise
coefficient arrays are immutable. Pure extended ranked-product vectors can be
reused across main-factor updates: the kernel reads only that problem and its
ordered participating input rows/values. The key preserves sorted `(rank,row)`
row order and full extended value bits, including signed zero. Rank magnitudes
and unselected values do not reach arithmetic and are omitted. Selected NaNs,
dumps, unsupported native control modes, malformed sizes and larger problems
keep the original path. Both consumers keep all WORK calls in place; numerical
kernels are source-identical and mutable factor state is neither read nor replayed.

The bound is 32 FIFO entries per existing provider. A test-only working-set
measurement shows only 22/23 exact inputs over 25,046/6,073 eligible calls on
Quad_2717/fx; the first four-entry version thrashes and regresses. Diagnostic
tracing is removed before frozen source/tests and official timing screens.
The final release test covers ties/order changes, exact low bits/signed zeros,
extremes/infinity, ignored values, clones, control-mode fallback and eviction.

At `kpttrrms`, full parity is 987/987, all 246 release tests pass and all 23
complete no-dump result/exit pairs agree. Two predeclared controlled Quad_2717
screens improve solve by 7.52/8.51%, each with all 32 pairs improving and full
noise/significance criteria passing. Instructions fall 4.04%. The full sweep's
two timing flags have fixed follow-ups and one reversed-role control each;
neither establishes regression. STATE.md retains exact statistics and current
snapshot ratios. Source/binary/build, ownership/kernel audit, both variants,
working-set measurement and all qualification data: `target/perf-rank-memo/`,
qualified variant under `v2/`. This establishes a local gain, not the overall
reference performance goal or safety of numerical producer-U reuse.

## Compact normal-Devex increment reuse

The existing shared provider also owns an eight-entry FIFO of pure extended
normal-Devex increments. Its Solver owns immutable problem columns; full-bit
`a_rq`, `wq_q`, `pq[]` and `tau[]` keys identify every arithmetic input. Entries
fill only used columns. Each call still adds its live old weight and performs
all original cost/twin/reset/state/order updates and WORK. The bound is at most
50 rows/100 total variables; mode 2, dumps, malformed/nonfinite inputs and
unsupported native/x87 modes retain the original production kernel.

Reset/owner-aware measurements show five keys per Quad owner and >99% saved
terms with eight entries. Native cross has no address-dependent peel or FP
status input. Production/frozen row kernels remain source-identical; reversal
proof and release-test per-row full-state differential checks cover hits and
fallbacks. At `vypkuxry`, two fixed target screens gain 11.10/9.83%; instructions
fall 5.90%. Full 987-run parity, all 254 release tests, formatting and 23 complete
benchmark pairs pass; fixed regression controls are healthy and establish no
slowdown. Scope/proof and every qualification record:
`target/perf-devex-delta/`; ownership-aware diagnostic:
`target/perf-devex-reuse/`. The overall reference performance goal is unmet.

## Nucleus coordinate workspace screen

The constructor repeatedly searches immutable pivot L/U vectors for a coefficient.
Temporary position arrays can preserve the first occurrence of duplicate keys;
clearing each pivot's touched coordinates makes the row array reusable for terminal
dense assembly. That assembly keeps last-write row positions, matching its old
HashMap. Node visits, matrix stores, arithmetic and WORK order remain unchanged.
Both tested variants preserve fast parity and bounded complete outputs; the final
shared workspace reduces PILOT4 instructions 1.46%, but its fixed timing screen
fails the declared significance/noise separation rules. The source and qualified
CLI are restored; there is no accepted speedup from this workspace. Exact evidence
and the sparse-only/final variant distinction are recorded in STATE.md and
`target/perf-nucleus-index/`.

## Guarded fresh replay screen

A strictly lower L cannot propagate a nucleus destination back into a singleton
output. A prototype checked this property during transpose construction, marked
producer U rows for each tail column, and required every reached nucleus output
to have a replacement with magnitude above tolerance. Under this guard it skipped
all numerical updates into the nucleus; uncovered/tiny/NaN/nontriangular cases
kept the original replay. The boundary was the first phase-0 pivot, not `loff`,
which precedes phase-4 singleton steps. Singleton arithmetic and order stayed
unchanged. The guard did not assume equality of replay and producer values.

The frozen full constructor and candidate agree in every fresh state field with
exact f64/F80 bits on the 26 cached ADLITTLE/BLEND/BORE3D/SC50A factors per
precision, plus randomized structures with extreme/nonfinite/duplicate values,
missing producer outputs and nontriangular fallback. Fast parity passes 940/940.
All 26 double factors have strictly lower L; 492/719 tail columns meet coverage,
and 572/655 structural update visits could be skipped. These are potential
visits, not executed operations or a runtime bound. Eight ADLITTLE/BORE3D
canonical result/exit pairs agree; instructions rise 0.10/0.14% there. PILOT4
primal instructions fall 0.57%, but its predeclared randomized screen establishes
no timing gain: median +0.26%, conservative interval [-0.53%, +1.78%], AA
[-2.03%, +1.59%]. The session stopped without confirmation/full gates, and
accepted source/CLI/PGO metadata are restored. Reachability, sorting and entry
materialization remained; eliminating them needs a separate proof and cost
screen. Artifacts: `target/perf-fresh-replay/` and
`target/perf-fresh-coverage/`; precise bounds and build hashes are in STATE.md
and the saved manifest. This is not evidence that direct numerical reuse is safe.

## Sparse-L level FIFO storage screen

The original level solver uses per-level growing vectors. Two representations
preserve its order: a per-node linked FIFO, and a flat FIFO with exact immutable
level capacities. The flat metadata retains the original per-node level prefix;
its cached starts count all nodes assigned to each level. Since each node is
enqueued at most once, every bucket fits that capacity. Queue reads use increasing
per-level indices, including same-level appends, while late appends to earlier
levels remain unprocessed as in the original solver. The numeric loop and WORK
expressions/order are unchanged; the full frozen L body remains source-identical.

The flat queue reuses `slot_counts` only during the level branch. It resets every
per-level used count and writes entries before publishing them; the later slot
sort resizes and zeroes every logical count before reading. Tests alternate both
consumers, prove original node levels/exact capacities, and compare complete
f64/F80 outputs/counts/WORK against the original on random duplicate/exceptional,
same/backward-level and captured SC50A cases. Final code adds no factor fields or
unsafe operations; the level-field header size is compile-time checked.

Both implementations pass fast parity and bounded complete outputs. Linked queues
add 0.40% PEROLD instructions under the preceding PGO profile. The flat candidate's
preceding-profile ADLITTLE instruction regression (+2.11%, on a model without the
level path) disappears with matching fixed39-case PGO (+0.008%). Final instructions
fall 2.00% on DEGEN2 and 1.15% on PEROLD; PEROLD initial cycles fall 1.07%.
Its fixed 32AB/32AA screen remains inconclusive: median -0.4954%, interval
[-0.78%, -0.15%], AA [-0.36%, +0.41%], 24/32 improving. No reserved confirmation,
extension or full gates follow. Source/CLI/build are restored; neither storage
change is accepted. Source/oracle/ownership proofs, all test/build/binary versions,
counters and controlled evidence: `target/perf-level-links/linked/` and
`target/perf-level-flat/`. These measurements do not establish a broader bound
on sparse L or paired-stage optimizations.

## Short ordinary tau-row screen

The virtual nucleus U row is below eight entries for 91.0% of PEROLD's 36,518
row dots and 84.6% of DEGEN2's 19,765 in bounded complete dual-release solves.
A candidate also bounds the numeric U row below eight and uses its first matching
coefficient directly for each virtual entry, preserving duplicate selection,
filtering and sequential subtraction order. Longer and extended rows keep the
original path; every WORK call remains in place. Original row-dot and frozen
needed-FTRAN bodies are source-identical. The length distribution is not a
measurement of the combined guard's coverage or removable cost.

The frozen final source/tests pass exact pipeline checks on 100,000 random/full
representations with independently generated duplicate coefficients, signed zeros,
NaN payloads/extremes and captured SC50A data; frozen f64/F80 FTRAN, fast parity
and bounded complete results also pass. Matching fixed39-case PGO reduces PEROLD
instructions 0.56%, but its predeclared 32AB/32AA timing screen is inconclusive
(median +0.06%, interval [-2.46%, +1.51%], 16/32 improving). No independent
confirmation, sample extension or full gates follow. Qualified source/CLI/build
are restored; there is no accepted short-row specialization. Exact distributions,
source/oracle proof, all versions/tests/results/counters and controlled samples:
`target/perf-short-tau/`.

## F2 snapshot lifetime and consumer audit

The existing F2 payload retains column/tau solved numeric arrays, input positions,
counts and raw-tau flags for later joint-U WORK accounting. The previous snapshot
is not read during the next full/needed solve's arithmetic. A persistent-workspace
prototype takes those buffers, zeroes/reuses them and moves the exact new snapshots
back after output extraction, retaining all numeric operations and WORK calls.
Repeated frozen f64/F80 comparisons and fast parity pass, but the initial accepted-
profile screen increases PEROLD instructions 1.09% and cycles 1.11% (PILOT4
+0.94/+0.77%). It is rejected before new PGO or controlled screens; qualified
source/CLI/build are restored. This extends the earlier move-only experiment;
allocation removal alone has no demonstrated benefit. Evidence and restoration:
`target/perf-f2-workspace/`.

Production reads of `x_col` and `x_tau` occur only in `third_work` and
`dense_pair_work` as `x_col[c] != ZERO` and `x_tau[c] != ZERO && tau_raw[c]`.
Ordered `col_in`, counts and factor structure remain separate required inputs.
No production dump emits these numeric F2 arrays. The consumed predicates are
now independently checked against original numeric state and original joint-U,
third and dense WORK routines: five original numeric/WORK bodies remain identical
after oracle/state identifiers, visibility and whitespace normalization only.
Repeated full/needed f64/F80 tests cover captured SC50A inputs, exact predicate
bits and WORK; 3,072 random original states exercise branch/count/list conditions,
signed zeros and nonfinites.

Two column-only storage prototypes reuse the output bitmap instead of cloning
`x_col`. V1's nonzero-only marker revisits zeros in the fallback output scan.
The final marker distinguishes unvisited/visited-zero/visited-nonzero, preserving
the original visitation test while exposing the exact nonzero predicate. Both
preserve fast parity and complete bounded outputs. The final implementation's
preceding-profile instruction increase is about 1% on PEROLD/PILOT4. Its unchanged
source/tests are subsequently verified once with matching fixed39-case PGO,
following the level-queue experiment's demonstrated profile effect. Final own-PGO
PEROLD instructions still rise 0.07%, cycles fall 0.51%; PILOT4 instructions rise
0.076%. This fails the predeclared instruction gate, so no controlled timing or
full gates follow. Qualified source/CLI/build are restored, and no storage variant
is accepted. Predicate equivalence remains verified; profitable compaction does
not follow from it.
Numeric snapshot/copy costs are not isolated by the full FTRAN self shares.
Source proofs, independent oracles and both complete screens:
`target/perf-f2-colflags/` and `target/perf-f2-colflags-pgo/`; the source dependency inventory remains in
`target/perf-f2-workspace/f2-consumer-audit.json`.

## Rounded ranked-product conversion reuse

The existing compact provider's two consumers have distinct needs: cost-change
uses extended products only, while the primal iteration also rounds every output
to double. A bounded prototype lazily retains those exact converted bits per
existing 32-entry result. Kernel and conversion bodies and WORK calls remain
unchanged. Original-trait differential checks cover special values, lazy fill,
clones and fallback; all six provider tests and fast parity pass with matching
fixed39 PGO. Diagnostic counters find 461,310 reused versus 1,080 computed values
on Quad_2717, and 171,072 versus 693 on fx_quad_obj_5371. These counts describe
eligible rounded calls, not all conversions in the solve.

Own-PGO instructions fall 1.37% on quad and 1.22% on fx. The first predeclared
32AB/32AA screen gives median -1.27%, interval [-2.79%, -0.235%] with AA
[-0.718%, +0.439%]. The bound fails noise separation despite sign p0.000268
and both orders improving. The experiment stops without confirmation, extension
or full gates; no conversion cache is accepted. Qualified source/CLI/build are
restored. Evidence: `target/perf-rank-rounded/`.

## Allocator instruction ancestry

Current accepted-source profiling uses matching PGO/TCMalloc, line tables and
frame pointers, after full profiling parity 987/987. Plain release DWARF on
DEGEN2 resolves only 40/92 sampled allocator Rose ancestors, and larger captured
stacks do not improve it. The frame-pointer build raises DEGEN2 whole-batch
instructions 11.28%; its shares are ancestry evidence, not release cost bounds.

Two bounded sixteen-solve PEROLD profiles (CPU4, dual release, cap1300; periods
4M and 3,791,733) put 152/341 and 186/400 allocator self samples under the
refactorization route, versus 104/341 and 115/400 under iteration. This resolves
broad instruction ancestry beyond the old allocation-frequency counts. Both
profiles represent 99.83/99.94% of direct batch instructions with no loss or
throttle. Complete results agree, and accepted release artifacts are unchanged.

Immediate allocation sites are still not proven by frame-pointer ancestry.
Disassembly verifies that the leading observed factor_structure return address
follows a call to nucleus: the immediate caller can be skipped by unwinding.
Debug return-site lines alone would falsely blame the subsequent sequence extend.
ELF file offsets must first map through PT_LOAD virtual addresses. Frequent
observed ancestry includes nucleus/structure, fresh installation, FTRAN/L and
controller/refactor work; inspect the allocating callee before proposing reuse.
Evidence and runnable analyzers: `target/perf-alloc-ancestry/`.

## Direct nucleus input construction

A complete prototype eliminates filtered temporary numeric columns by building
nodes directly, converting in original forward order and reversing node segments
to preserve node IDs/chains exactly. The original numerical kernel and WORK body
are source-identical. Frozen producer/initializer checks match complete state
across 1,024 f64/F80 families, and matching-PGO tests/fast parity pass. Bounded
complete results also match. PEROLD instructions fall only 0.16%, while cycle
median rises 1.21%; small-model instruction reductions are also below 0.2%.
The trial stops at initial counters without controlled screens or full gates.
Qualified source/CLI/build are restored; no input-allocation change is accepted.
Evidence: `target/perf-nucleus-input/`.

## Nucleus chain-maintenance counts and code size

A temporary standard dump kind counts list work while preserving every numeric
operation and WORK call. Complete ADLITTLE/DEGEN2/PEROLD bounded results equal
the accepted CLI; source, dump registry, CLI and build are then restored exactly.

| Model | Nucleus calls / dense exits | Front-shift IDs | Snapshot IDs | Clean visits / removed dead IDs |
|---|---:|---:|---:|---:|
| ADLITTLE | 9 / 0 | 2,103 | 2,832 | 1,900 / 489 |
| DEGEN2 | 16 / 0 | 86,995 | 33,771 | 22,900 / 7,908 |
| PEROLD | 27 / 11 | 1,424,407 | 541,599 | 184,154 / 40,353 |

PEROLD elimination scans 350,277 IDs for 88,400 existing numeric updates and
26,807 candidate fill terms. The measured ID movement alone does not establish
that an intrusive-list rewrite pays off; it has no cycle attribution. Accepted
release ordinary nucleus code is 31,086 bytes with a 0x6c8 local stack frame;
dense LU is already a separate function. The full dense handoff construction
and output installation remain inside nucleus. Frame-pointer annotation has
238 nucleus self samples with 18.9% on stack-addressing instructions, but its
codegen differs from release. Splitting the handoff can be tested without
changing numeric operations, list layout or bookkeeping; code size and stack
changes need direct release verification before timing claims.
Evidence: `target/perf-nucleus-chains/`, with source annotation in
`target/perf-alloc-ancestry/perold-fp.nucleus.asm`.

## Dense nucleus handoff extraction

The dense handoff was extracted into a non-inlined helper with identical source
statements except the boolean success return. Reconstructing the original kernel
matches normalized source exactly; numeric kernels and WORK are untouched.
Matching-PGO frozen dense-LU tests and fast parity pass, and all bounded complete
results agree. Ordinary nucleus code shrinks from 31,086 to 26,157 bytes and its
stack allocation from 1,736 to 1,496 bytes; the new helper absorbs dense LU.
PEROLD instructions fall only 0.023%, cycles rise 0.041%. This fails the initial
useful-counter gate: no controlled timing/full gates follow, and all accepted
artifacts are restored. Code size/stack reduction alone does not establish a gain.
Evidence: `target/perf-nucleus-outline/`.

## Immutable-slot heap-key compaction

U-order and row-elimination borrow the factor immutably. Slots and positions
are inverse permutations, so slot priority alone preserves the original tuple
priority, including identical duplicate entries. A complete two-site prototype
halves heap entry size and decodes each popped position through at[slot]. Frozen
original bodies check both consumers over random slot permutations, sparse/dense
inputs, duplicate initial entries and f64/F80 results/order/WORK. Reconstruction
matches source except key representation/decode; matching-PGO tests and fast
parity pass. Bounded complete results agree.

PEROLD instructions fall 0.249%, cycles rise 0.077%. The first predeclared screen
(seed7001, 32AB/32AA) gives median -0.124%, interval [-0.414%, +0.164%], AA
[-0.457%, +0.243%], sign p0.189. It is inconclusive; no confirmation, extension
or full gates follow. All accepted artifacts are restored. Equivalence is
verified, profitability is not. This differs from the earlier numeric-LU pending
step bitset trial. Frame-pointer annotation locates about 0.39/0.37% Rose self
samples in tuple max/min heap pops and 4.8% of update-local self samples in
inlined heap pushes; these are profiling-build locations, not release cost bounds.
Evidence: `target/perf-slot-heaps/`, `target/perf-heap-audit/`.

## Unique AA-pattern neighbour production

Before sorting, row-incidence generation marks now retain each non-dense
neighbour once and omit the diagonal. Original dense-column extension and CSR
output remain untouched. Complete Pattern equality against the frozen original
constructor preserves all downstream symbolic choices, floating degree/work
arithmetic and numerical factorization inputs. Both generated-input and six-model
release comparisons pass. Matching fixed39 PGO gives target instructions -2.40%;
independent DEGEN3 barrier screens qualify median solve gains of 1.71/1.93%.
Full parity 987/987, all release tests, formatting and all 23 complete benchmark
pairs pass; predeclared reversed-role controls establish no unrelated regression.
No numeric kernel, dump or baseline changes. Current proof/build/measurements:
`target/perf-aat-unique/`; predecessor sort-caller attribution:
`target/perf-barrier-refresh/`.

## Quotient-order row ownership

Row rebuilding takes ownership of the existing adjacency buffer, stably compacts
its element and variable segments, then inserts the pivot at the saved segment
boundary. Deferring insertion prevents clobbering unread tail values. All ordered
hash additions, degree updates, resulting lists and WORK expressions are retained.
The independently frozen complete original ordering verifies permutations/NZ and
Ops/WORK bits on generated graphs (including large floating-hash seeds) and real
small-model patterns/levels. Matching-PGO tests, full parity 987/987, all 250
release tests and 23 complete benchmark pairs pass. Independent DEGEN3 barrier
screens establish median gains of 0.71/0.64%, with instructions -2.20%.
Predeclared reversed-role controls establish no unrelated regression. Native
assembly still spills the hash accumulator per term; this gain does not prove
register retention. Rule/proof/build and complete measurements:
`target/perf-order-reuse/`.

## Completed ordering-hash lifetime

The completed hsh[i] write precedes vector insertion, ending the hash local's
lifetime before a possible allocation. Its next read remains after row rebuilding;
all floating additions and work/degree operations keep their original order.
Existing full-original-order tests include rounded-hash seeds and model levels.
Matching release assembly now keeps both hash-loop accumulators in registers.
Independent target screens qualify median solve gains of 2.02/2.10%, with
instructions -0.40% and cycles -3.91%. Full parity, all release tests and complete
benchmark pairs pass; fixed controls establish no unrelated regression.
Source/build, exact store-relocation proof, assembly and all raw qualification:
`target/perf-order-hash-store/`.

## Current ordinary L-stage refresh

Accepted release `vypkuxry` is sampled directly after bounded ADLITTLE/DEGEN2
checks, then sixteen PEROLD dual solves, cap1300. All complete results agree.
The period4M recording has4302 samples, zero loss/throttle, and represents99.87%
of the counted17.230B batch instructions. Whole-batch sample shares include
L-stage10.18%, first/needed FTRAN8.83/7.23%, nucleus5.25%, auxiliary factor4.09%
and fresh3.44%; launcher/startup/parsing/validation remain in scope. These are
self observations, not inclusive or removable cost bounds. Annotation still
locates numeric scatter, dense WORK accumulation and stack-indexed eta lanes.

A long-only four-register eta prototype preserves every scalar accumulation
chain, reduction and WORK site. Opaque original/frozen calls validate before
the rewrite and pass with matching PGO, all16 MXCSR modes, nonfinites/extremes
and captured values. The earlier nonopaque scalar-fallback NaN test failure
also reproduces on unchanged accepted source, so it cannot establish a SIMD
bug. Eight FTRAN/L tests and fast940 parity pass. PEROLD instructions fall1.83%,
but its fixed32AB/32AA screen fails conservative noise separation: median
-1.506%, interval[-2.550%,-0.759%], AA[-0.882%,+1.017%]. No confirmation,
extension or full gates follow; all accepted artifacts are restored exactly.
Source/proof, unchanged-source oracle failure/pass, mode tests, assembly and
full timing evidence: `target/perf-eta-long/`. Accepted release profile and
counter/event audit: `target/perf-ordinary-devex/`.

## Dense L WORK lifetime screens

Exact immutable per-column charge caching retains every WORK addition and all
numeric traversal. A guard limits stored `len as f64 * 2.0` charges to counts
at most 2^53, whose conversion/scaling is exact in every MXCSR mode; larger
columns retain the original expression. Independent complete L-stage and dense
WORK bodies verify generated/extreme/captured f64/F80 cases. Matching fixed39
PGO, seven FTRAN/L tests and fast940 parity pass; all bounded no-dump result/exit
pairs agree. PEROLD instructions fall only 0.073%, with two of three paired
cycle changes positive. The cache is restored without timing screens.

A separate complete-loop extraction leaves factor layout and expressions
unchanged. Source reversal reconstructs both production methods and the helper
from the original loop. Its ordinary assembly retains WORK in XMM0; matching
PGO, six FTRAN/L tests and fast940 parity pass. PEROLD instructions fall0.190%
and cycles change -0.033/+2.951/-0.224% across the three direct pairs. The desired
register lifetime is established, but a useful runtime gain is not. No official
timing or full gates follow; all accepted artifacts are restored exactly.
Ownership/bound and source proofs, independent oracles, build/test provenance,
assembly and complete paired results: `target/perf-lwork-charges/` and
`target/perf-dense-l-outline/`. These local observations do not bound total
FTRAN bookkeeping cost or establish safety of producer-U reuse.

## Fixed BTRAN scratch storage

The current accepted Quad profile retains extended BTRAN as the largest self
symbol (13.40% of batch instruction samples, including launcher/startup in scope).
The period4M recording has7231 samples, zero loss/throttle and represents99.87%
of counted28.961B batch instructions;7201 samples belong to Rose. All bounded
complete results agree. Current data/annotations: `target/perf-quad-devex/`.

A trial moves five fixed-size numeric/marker buffers to initialized stack prefixes
only for extended factors with at most32 rows. Dynamic support/output lists keep
their original storage. Prefix initialization, exclusive arrays, Copy element
types and checked bounds justify slice creation from MaybeUninit. The partial
L-stage handoff copies the exact initialized values instead of moving its Vec.
Reversal reproduces the complete original BTRAN; independent whole-body tests
check result/WORK bits, row cache, density and lazy inverse/level state. Separate
boundary/suffix and low-bit partial-path cases pass, as do all nine FTRAN/L
tests and fast940 parity with matching fixed39 PGO. Quad/fx instructions fall
0.396/0.142%, but cycles rise1.50/1.46%. No official timing/full gates follow;
all accepted artifacts are restored. Source/proof, oracle/tests, build/assembly
and full result/counter pairs: `target/perf-btran-stack/`.

A separate test-only division measurement excludes independent Devex-shadow
arithmetic. Quad/primal5959 have131,998/172,580 zero numerators out of599,835/
817,754 divisions. All solver statuses/counts and optimal numeric fields match
the accepted CLI. Accepted dense-U assembly already inlines its zero return
before quotient helpers, so these counts do not establish removable ABI cost.
The sparse-input out-of-line call remains distinct; no production change is
made from these counts. Exact scope/counters/source/restoration and symbol/
relocation proof: `target/perf-f80-div-shapes/`.

## F80 operand-order lifetime screen

Accepted BTRAN instruction annotation includes conditional field selection in
inlined F80 addition/subtraction. These locations alone do not establish a
removable cycle bound. Replacing value tuples with ordered references preserves
all numeric expressions; an independent complete original routine and source
reversal protect every path. All ten F80 tests (including random/full-bit/
boundary/captured internal-field comparisons) and fast940 parity pass with
matching fixed39 PGO.

The first variant saves0.264% Quad instructions but increases cycles5.3–5.9%
in all three pairs, with sharply higher branch misses. A second uses explicit
unpredictable selection for the same references. Its branch-miss penalty is
largely removed, but target instructions rise0.247% and median cycles1.269%.
All bounded complete results agree. No official timing/full gates follow;
source, accepted CLI and metadata are restored exactly. This does not establish
field-copy removal as profitable or warrant assigning the sampled conditional
moves their apparent instruction-sample share as removable runtime. Both complete
variants/proofs/builds/tests, native assembly and paired results:
`target/perf-f80-add-borrow/`, including `v2/`.

## Pure quotient reuse scope

Global per-CLI-equivalent streams retain full raw F80 operand fields across all
Solver resets and exclude test-only Devex-oracle divisions. A fixed512-entry
four-way table would reuse99.3/87.1/99.0% of eligible fx/Quad/primal5959 divisions.
The corresponding32KiB const thread-local implementation preserves the complete
original arithmetic, exact keys and all WORK. Frozen original, collision,
eviction and cold/warm boundary/random/captured tests pass, along with all
eleven F80 tests and fast940 parity. It still loses cycles: fx+0.845%, Quad
+5.671%. Full bounded outputs agree; no official timing/full gates follow.
The accepted source/CLI/metadata are restored. Exact scope/simulations:
`target/perf-div-reuse/`; complete prototype/proof/build/tests/counters and
relocation/call audit: `target/perf-div-memo/`. The dense-U division remains
inlined in both binaries.

A separate diagnostic tags only BTRAN's three diagonal-division sites. Fresh
factor IDs are monotonic and survive Clone; no pointer identity is used. Full
raw numerator/divisor keys are retained across all updates. Quad has231,284
eligible queries,492 distinct (owner,column,operand-pair) keys across14 owners
and96 columns, and a maximum14 keys per column. FIFO1/2/4/8/16 hits are42,404/
53,950/79,948/148,158/230,792. Sixteen-entry fx/primal hits are58,829/221,325
of59,284/222,780, with per-column maxima8/12. This justified the tested16-entry
bound; both column-indexed implementations now lose at counters, as recorded
below. A cache shared by
fresh-factor clones must still key both numerical operands across updates;
its values are pure quotients, not mutable factor results. Original arithmetic
and every WORK site remain required. Data, lifetime rules, exact tagged source,
full optimal numeric/status/count validation and restored-source proof:
`target/perf-diag-div-reuse/`.

The per-column implementation shares only pure full-key quotient values through
one Arc/Mutex per compact fresh owner (one lock per BTRAN). The complete original
BTRAN and unchanged arithmetic are independently verified; collision, eviction,
updates, clone sharing, raw low bits, nonfinites and mode tests pass. Matching PGO
and fast940 pass. Linear fingerprints give Quad instructions+0.710%, cycles
+3.368%, branch misses+11.82%; integer SIMD fingerprint masks give+0.771%,
+3.000%,+0.976%. FX cycles rise6.482/5.827%. All bounded full result/exit pairs
agree. Neither variant reaches official timing/full gates. Both are restored;
source/proofs/build/tests/assembly/raw counters: `target/perf-diag-div-memo/`
and `v2/`.

A denominator-only diagnostic finds at most3/4/4 distinct normalized divisor
mantissas per FX/Quad/primal5959 owner-column. Four-entry FIFO parameters would
hit59,082/231,157/222,218 of59,284/231,284/222,780 eligible divisions. Mantissa
alone is a valid key for integer reciprocal parameters; it is insufficient for
quotient-result reuse. Diagnostic scope, all bounds, status/count and optimal
full-field checks: `target/perf-diag-den-reuse/`.

A standalone optimized probe keeps the complete original F80 module as oracle.
For normalized denominators, a precomputed65-bit reciprocal gives an exact
quotient/remainder after at most one correction; normal-result rounding follows
from that exact remainder and other inputs retain original division. Its own
integer proof,2,000,050 integer and1,037,249 F80 comparisons pass. Fixed eight
CPU-4 alternating pairs give about0.45 arithmetic time ratio with parameters
precomputed and equal nonzero checksums. This establishes arithmetic feasibility
only. Integrated storage, lookup and prepared-on-write designs are now tested
below; none qualifies for adoption.
Full sources/proof/compiler/results: `target/perf-preinverse-probe/`.

The integrated four-entry per-column parameter cache shares one owner lock per
BTRAN and keys the normalized divisor mantissa; live inputs still determine
sign, exponent and rounding. Complete original division/BTRAN oracles, all11
F80 and8 FTRAN/L tests, fast940 and complete bounded result/exit pairs pass for
three matching-PGO designs. Original arithmetic outside three division sites
and all WORK remain unchanged. Initial lookup raises Quad instructions0.249%;
earlier fallback checks and removal of proved-unreachable normal rounding carry
still raise them0.089%. The latter removes no original exceptional boundary
path. Neither reaches official timing.

The prepared variant pairs one reciprocal with each current diagonal entry.
Fresh factors prepare it, Clone copies it, and a single setter routes the only
production diagonal mutation and refreshes the parameter. This removes lookup
and shared synchronization; ordinary/large factors retain original division.
Clone/update/extreme/fallback checks validate pairing. Quad instructions fall
0.280%, cycles3.631%; fx0.216%/1.280%. The fixed32-AB/32-AA Quad screen24217 is
inconclusive (median-0.631%, interval[-4.598%,+2.406%], sign p.430; healthy AA,
execution orders disagree). No confirmation24289 or full candidate gates;
accepted source/CLI/metadata/profile restored, formatting and full987 pass.
The arithmetic feasibility remains valid but its tested integrations have no
qualified runtime benefit. Source/proofs/builds/tests, ownership/mutation
checks, complete counters/results and screen: `target/perf-diag-preinverse/`
with `v2/` and `v3/`.

Current accepted ordinary-primal PILOT4 profiling (ADLITTLE/DEGEN2 first) gives
13.92/10.49/8.37/8.17% whole-batch self shares to solve_primal/BTRAN/FTRAN/row
update, respectively. All16 solves are optimal690 iterations/16 refactors.
Period4M has1983 samples, no loss/throttle,99.69% instruction coverage;95
samples are Python dispatch. The annotated ranked-row-shaped scatter inside
solve_primal is an inferred source attribution. This is current native self
cost, not allocator ancestry or proof that a scatter rewrite pays off. Existing
rank-row storage/active-row/rank-order failures still apply. Exact commands,
results/hash/build/count/audit and annotations: `target/perf-primal-refresh/`.

## Remaining proof obligations

Before a factor ownership prototype, separate column-copy and immutable-factor
allocation costs, and inventory auxiliary
reads between refactorizations. Replay numerical producer/consumer values only
through the standard dump/comparison tools; stop at the first mismatch. Retain
existing arithmetic and WORK booking for the first ownership experiment. The accepted pure-result reuse above does not establish numerical producer-U
identity or safety of broader factor sharing.
