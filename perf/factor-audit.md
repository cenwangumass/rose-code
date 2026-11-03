# Factor production and ownership audit

Ownership audit source: b421f8ab; guarded replay screen: accepted d897aa91.
This is a source/lifetime audit, not a measurement of attainable speedup. Direct numerical reuse remains unproven;
guarded elimination of overwritten replay arithmetic is feasible on tested
cases but has no established runtime benefit (below).

| Producer | Consumers and lifetime | Consequence |
|---|---|---|
| `lu_order::Structure` | Pivot order, phases, L rows/multipliers/raw entries, diagonal and partial U feed `RefStruct::fresh`; structure and incoming columns remain live for WORK and pricing installation. | Consuming the structure at the current fresh call would invalidate later consumers. Narrow ownership transfers need an installation payload first. |
| `RefStruct::fresh` | Copies diagonal; remaps L destinations into flat forward storage and builds its transpose. Replays elimination for tail U, then overlays Markowitz entries. | Producer L row IDs and consumer step IDs differ. Moving nested producer vectors directly cannot replace the remap. Direct numerical U reuse requires intermediate equality evidence; guarded skipping of overwritten replay values is a separate case below. |
| Primal `reference_order_for` | Returns the accepted growth-check mirror for installation (`simplex.rs`, `growth_factor`). Retried attempts book their own WORK. | Duplicate accepted growth-check construction is already eliminated. Do not count it as remaining savings. |
| Primal `quad_mirrors` | One fresh extended factor is deep-cloned for the second twin; diagnostic mode retains two builds. A four-entry small-model cache also clones factors. | Duplicate ordinary twin numerical construction is already eliminated. Sharing immutable L/permutation/transpose storage is still a candidate; mutable U, counters, scratch, spikes and cells must stay independent. |
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

The RHS twin has only one consumer, `recompute_primal`, guarded by zero main
factor updates. A prototype removed the ordinary second factor and ran the
unchanged force-dense RHS solve on the main factor, moving aside and restoring
all touched solve/accounting state. Diagnostic mode retained both original
builds and factors. The main numerical factor remained unchanged. Differential
checks against an independent twin passed for f64/F80 with prior solve history,
full result bits, retained state and repeated/extreme inputs; fast parity passed
940/940. Instructions fell 2.15% on fx_quad_obj_5371 and 1.17% on Quad_2717,
but a predeclared randomized Quad screen reported +0.47% median time, interval
[-0.77%, +1.73%], and no gain. The prototype was rejected and the accepted
source/CLI/PGO metadata restored. Thus isolated borrowing is numerically feasible
on tested cases but has not established a runtime benefit; the fresh-factor
cache's main copy and the RHS routine's temporary work remain present. Exact
artifacts and bounds are in STATE.md.

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
