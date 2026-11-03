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

## Remaining proof obligations

Before prototyping, measure current factor self cost on these bounded witnesses,
separate column-copy and immutable-factor allocation costs, and inventory auxiliary
reads between refactorizations. Replay numerical producer/consumer values only
through the standard dump/comparison tools; stop at the first mismatch. Retain
existing arithmetic and WORK booking for the first ownership experiment. No
solver change or speedup is claimed by this audit.
