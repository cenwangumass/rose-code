# Factor production and ownership audit

Ownership audit source: b421f8ab; guarded replay screen: d897aa91;
compact RHS reuse: prxvqvrr; pure reduced-cost reuse: ztrylnzw.
Cost claims below identify their measured scope. Reusing producer U values in
the consumer factor remains unproven;
guarded elimination of overwritten replay arithmetic is feasible on tested
cases but has no established runtime benefit (below).

| Producer | Consumers and lifetime | Consequence |
|---|---|---|
| `lu_order::Structure` | Pivot order, phases, L rows/multipliers/raw entries, diagonal and partial U feed `RefStruct::fresh`; structure and incoming columns remain live for WORK and pricing installation. | Consuming the structure at the current fresh call would invalidate later consumers. Narrow ownership transfers need an installation payload first. |
| `RefStruct::fresh` | Copies diagonal; remaps L destinations into flat forward storage and builds its transpose. Replays elimination for tail U, then overlays Markowitz entries. | Producer L row IDs and consumer step IDs differ. Moving nested producer vectors directly cannot replace the remap. Direct numerical U reuse requires intermediate equality evidence; guarded skipping of overwritten replay values is a separate case below. |
| Primal `reference_order_for` | Returns the accepted growth-check mirror for installation (`simplex.rs`, `growth_factor`). Retried attempts book their own WORK. | Duplicate accepted growth-check construction is already eliminated. Do not count it as remaining savings. |
| Primal `quad_mirrors` | Each compact cache entry retains a pure main-factor starting point and an isolated RHS provider. Main factors deep-clone before mutation; repeated RHS owners share the provider. Numeric-build diagnostics retain their two builds; uncacheable factors keep direct RHS ownership. | Duplicate ordinary construction and repeated compact RHS-factor cloning are eliminated. RHS results and pure WORK estimates have separate exact keys and four-entry bounds. Mutable main U/counters/scratch remain independent. |
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

## Remaining proof obligations

Before a factor ownership prototype, separate column-copy and immutable-factor
allocation costs, and inventory auxiliary
reads between refactorizations. Replay numerical producer/consumer values only
through the standard dump/comparison tools; stop at the first mismatch. Retain
existing arithmetic and WORK booking for the first ownership experiment. No
solver change or speedup is claimed by this audit.
