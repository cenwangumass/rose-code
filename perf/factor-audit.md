# Factor production and ownership audit

Source inspected: accepted solver b421f8ab. This is a source/lifetime audit, not a
measurement of attainable speedup. Numerical reuse remains unproven.

| Producer | Consumers and lifetime | Consequence |
|---|---|---|
| `lu_order::Structure` | Pivot order, phases, L rows/multipliers/raw entries, diagonal and partial U feed `RefStruct::fresh`; structure and incoming columns remain live for WORK and pricing installation. | Consuming the structure at the current fresh call would invalidate later consumers. Narrow ownership transfers need an installation payload first. |
| `RefStruct::fresh` | Copies diagonal; remaps L destinations into flat forward storage and builds its transpose. Replays elimination for tail U, then overlays Markowitz entries. | Producer L row IDs and consumer step IDs differ. Moving nested producer vectors directly cannot replace the remap. Reusing numerical U requires intermediate equality evidence. |
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

## Remaining proof obligations

Before a factor ownership prototype, separate column-copy and immutable-factor
allocation costs, and inventory auxiliary
reads between refactorizations. Replay numerical producer/consumer values only
through the standard dump/comparison tools; stop at the first mismatch. Retain
existing arithmetic and WORK booking for the first ownership experiment. No
solver change or speedup is claimed by this audit.
