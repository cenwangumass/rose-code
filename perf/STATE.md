# Performance state

Current state of the work in [PLAN.md](PLAN.md). Rewrite in place; no log (history is `jj log`).

## Dual simplex instruction stack (local, unpushed)

Stack on `main` `a5050f4a` to `23039b0a` (14 commits: walk-list rebuild, BTRAN level sweep with its
repeated-entry guard, BTRAN output packers, slim FTRAN WORK accounting, dense-route WORK from the real
solve, pricing with cached basic bounds, MPS parse without per-line allocations, pivot-row scatter and
word-block skip, L-stage level solve and eta dots, needed-U and full-U scatters). Every commit is bit-identical:
release tests 372 passed with 22 ignored, full parity normal and `factor-shadow` 987/987, corpus diff 0 over
564 runs. `perf stat -e instructions:u -r 3`, dual release, current remote host, against the reference:

| Model | `main` | `23039b0a` | Change | Ratio |
|---|--:|--:|--:|--:|
| GREENBEA | 2573.4 M | 2237.9 M | -13.04% | 1.554 |
| DEGEN3 | 985.2 M | 892.6 M | -9.41% | 1.549 |
| DEGEN2 | 124.12 M | 112.39 M | -9.45% | 1.690 |
| 80BAU3B | 1157.7 M | 1026.7 M | -11.32% | 1.664 |
| PEROLD | 688.8 M | 613.5 M | -10.93% | 1.834 |

Instruction-ratio geomean 1.655 (`main` 1.856). Fixed 23-row wall-time table on `23039b0a`: 0/23 Rose wins,
dual rows geomean 1.614 (degen3.dual-r 1.307, greenbea.dual-r 1.545, degen2.dual-r 1.706). Remaining
GREENBEA gap by family: FTRAN ~215 M, factor construction ~143 M, BTRAN ~70 M, update 59 M, pricing and
glue up to 195 M; the exact micro-rewrites left are each below 0.5%.

## Replacement candidate

The retained needed-U numeric-row component and ordinary dense-tail slot fix are
implemented; the full backend replacement remains incomplete and unaccepted.
Strict state checks pass738 calls in ten complete small-model/profile runs;
raw output/RHS kernels pass924 pivots on AFIRO, SC50A, BLEND, BORE3D and DEGEN2
(release), including signed zeros. Numeric factor projections and update
structure agree within the comparator's documented coordinate/order mappings.
ADLITTLE raw rho row50 remains native+1/Rose-1 at pivot48 (G-row coordinates).
Native physical layout/visitation, shared pool lifetimes, constructor/auxiliary
duplication and cached update elimination remain unresolved.

The matching current-host fixed23-row reference benchmark completed with0/23
Rose wins. Dual ratios range1.161–7.533; DEGEN2 is65.925/17.934ms
(3.676x). These three-sample medians establish a large remaining gap, not a
candidate-versus-baseline speedup. The full measurement includes all methods;
parameters, binaries, raw outputs and exits are saved in
`target/factor-rows/reference-full`. The corpus frozen baseline is unchanged:
932/940 fast cases are identical; the other8 have exactly29 +0-to--0 tau changes
in18 of44019 paired records, all other44001 records identical and final results
bit-identical except elapsed time. Full release tests pass295 with22 ignored, including all131 library tests.
The full frozen gate has970 unchanged and17 changed runs out of987; all changed
streams are DTAU and their first differences are tau +0-to--0. Only the8 fast
cases were exhaustively classified. The additional9 full-tier changes have not
all been independently checked against native captures; baseline unchanged.
The numerical traversal fix is committed at0c00d2ce; the retained-row component
remains an uncommitted candidate without a qualified performance gain. Evidence:
[evidence/dual-factor-candidate.json](evidence/dual-factor-candidate.json).

The accepted baseline and its prior forensic measurements below remain separate
from this unaccepted candidate.

## Current host: AMD EPYC virtual CPU, family 23

Source Rose `main` = `789467e2` (dual factor capture-gaps C1X+C2X+C3X on `a416b9c8`), pinned Rust 1.98.1,
ordinary ThinLTO release (no `factor-shadow`), glibc, CPU 2. CLI SHA256
`7aa21f98f0c817f7e96ac27013201460f33c9fd921ca60f0405c2d1150a8224c` (previous table: `18061fe9`, CLI SHA256
`8483049968e1aa8be28e7003ec5a0acdeb80a77a07fea05b50c6142a54d38410`). Host line printed by the runner:
AMD EPYC Processor (with IBPB), kernel 7.0.0-38-generic (remote timing machine via `rtime`). The ordinary
baseline remains intact. The separately qualified fixed39 PGO/system-TCMalloc configuration below has its own
ELF and table. Historical EPYC 9V74/TCMalloc tables remain separate configurations.

Optional PGO build (measured 2026-11-13 on `04fbad99`): `perf/build_pgo.py` passes every gate (release tests, full parity normal and
`factor-shadow` 987/987, corpus diff 0) and is offered as an optional release build, not the default. Against plain release
(remote `perf stat -r 3`, dual release) it cuts instructions by about 1-3% and cycles on most rows (primal and barrier included), but
`witness_dual_r2_mode_switch` is 7-33% slower in cycles in every batch with unchanged instructions, also when the training set adds
capped mode_switch and GREENBEA dual runs. The profile is valid only for the toolchain it was trained with (rustc 1.98.1, LLVM
22.1.8) and must be retrained, and the gates rerun, after any source or toolchain change. `codegen-units = 1` and fat LTO are not
usable: both change NaN payloads or signs (release tests `chol` dgemm/trsm and `gathered_dot` fail).

The current fixed 23-row reference comparison (re-run 2026-11-11 on `789467e2`) is complete: Rose is faster on
**0/23**. Each solver has one warmup and three measured fresh-process runs per
row, rotating order on the same CPU. Timers are Rose `diagnostics.seconds` and
Gurobi `Runtime`; all 23 rows' complete results/exits are stable within each solver. Valid nonoptimal statuses
and Rose exit 3 are preserved (primal_5959 rows, greenbea.barrier-r, greenbea.pdhg). Barrier
uses its normal 1000 iteration limit; GREENBEA PDHG uses the specified 20000.
Reference logging is explicitly off. Three-sample medians are descriptive,
not an optimization acceptance screen or cross-solver numerical-alignment proof.

The dual rows are what the dual factor work since `18061fe9` changed (the whole dual factor backend effort, not only
the capture-gaps stack; the stack alone is -5.6% to -6.5% instructions): Rose ms fell on all dual rows with a dual-heavy
workload (for example degen2.dual-r 58.7 -> 47.4, degen3.dual 462-489 -> 388-414, 80bau3b.dual-r 566 -> 446,
greenbea.dual-c 1332 -> 1068, mode_switch.dual-r 1344 -> 1062). Primal/barrier/PDHG rows are unchanged code;
their Rose ms moved by -8% to +34% and the Gurobi column also rose 0-25% (for example 80bau3b.primal-c Rose
+34%, Gurobi +20%; mode_switch.dual-r Gurobi 164 -> 283), so this run's host was slower or noisier than the
previous one. Compare the ratio column, and treat single-row ms deltas as noise-limited.
"Rose ms before" is the previous table; the ratio is from this run only (previous ratios in the JSON).

| Row | Rose ms | Rose ms before | Gurobi ms | Rose / Gurobi |
|---|---:|---:|---:|---:|
| degen2.primal-c | 70.017 | 71.128 | 22.629 | 3.094 |
| degen2.dual-r | 47.361 | 58.745 | 19.122 | 2.477 |
| degen3.primal-c | 964.372 | 883.168 | 313.855 | 3.073 |
| degen3.primal-r | 1118.214 | 1216.082 | 398.545 | 2.806 |
| perold.primal-c | 1040.952 | 1024.754 | 308.519 | 3.374 |
| pilot4.primal-r | 111.932 | 93.306 | 35.008 | 3.197 |
| 80bau3b.primal-c | 1233.757 | 922.550 | 226.345 | 5.451 |
| primal_5959.primal-r | 108.012 | 99.369 | 74.215 | 1.455 |
| degen3.dual-c | 414.245 | 489.323 | 173.096 | 2.393 |
| degen3.dual-r | 387.979 | 462.506 | 162.951 | 2.381 |
| perold.dual-r | 220.981 | 259.579 | 76.904 | 2.873 |
| 80bau3b.dual-r | 446.100 | 565.991 | 153.476 | 2.907 |
| greenbea.dual-c | 1068.384 | 1331.731 | 411.515 | 2.596 |
| greenbea.dual-r | 1333.392 | 1394.194 | 458.471 | 2.908 |
| mode_switch.dual-r | 1062.425 | 1343.815 | 283.299 | 3.750 |
| primal_5959.dual-r | 160.409 | 127.878 | 114.422 | 1.402 |
| quad_2717.dual-r | 159.920 | 153.781 | 125.055 | 1.279 |
| fx_quad_obj_5371.dual-r | 65.010 | 47.779 | 47.433 | 1.371 |
| degen3.barrier-r | 317.058 | 271.010 | 150.438 | 2.108 |
| 80bau3b.barrier-r | 443.008 | 405.132 | 200.357 | 2.211 |
| greenbea.barrier-r | 1860.167 | 1651.639 | 787.428 | 2.362 |
| 80bau3b.pdhg | 4685.953 | 3744.162 | 4118.993 | 1.138 |
| greenbea.pdhg | 2261.199 | 1898.783 | 1809.566 | 1.250 |


The reference was reconstructed from the supplied archive and hash-checked
patches at `target/perf-reference/reconstructed/libgurobi130.so`, exact SHA256
`26a505171314b3300b754b4eff5f5719366b47ccf6f3d59770df99e3f256ee30`.
Environment initialization and AFIRO profiles are qualified. The native harness
now checks Runtime/status/iteration reads, honours parameter precedence, records
effective overrides and preserves failures/timeouts. Parent harness commits:
`7e089564`, `a86024e2`. The first sweep prefix stopped on an incorrect harness
exit check and is excluded; the corrected complete panel is used exclusively.
Evidence: [current-reference-comparison-789467e2.json](evidence/current-reference-comparison-789467e2.json) (this table, with all samples, previous ms and previous Gurobi ms),
[current-reference-comparison.json](evidence/current-reference-comparison.json) (previous table),
all manifests, raw output, source/runner/library/model hashes and qualified
small/error witnesses in `target/perf-reference/`.
[reference-availability.json](evidence/reference-availability.json) records reconstruction.
The earlier single-binary host baseline in evidence/current-host.json remains
separate from this paired solver comparison. Full frozen parity 987/987 and
all 292 runnable release tests pass; numerical source and accepted CLI unchanged.

### Dual factor FTRAN redesign stack (`4a5e881a`), two runs

Same method and host as the table above (`compare_reference.py`, reconstructed libgurobi130.so sha 26a50517, medians of
3 after 1 warmup). Commits since `93540e37`: A/B coefficient lists, dense L-to-U handoff, sparse U order without touched
sets, numeric U row pool, needed and full U kernels, U ordering by tail sweep/heap as the reference does, pivot row
residual and scatter, retained per-iteration scratch, L kernels. Numbers unchanged (parity 987/987 normal and shadow,
corpus diff 0). GREENBEA dual instructions 3560.9 M -> 2747.3 M (reference 1436 M).

| Row | Rose ms run 1 / run 2 | Gurobi ms run 1 / run 2 | Rose / Gurobi run 1 / run 2 | Ratio in the table above |
|---|---:|---:|---:|---:|
| degen2.dual-r | 36.9 / 35.7 | 18.5 / 17.2 | 1.995 / 2.079 | 2.477 |
| degen3.dual-c | 269.5 / 254.4 | 165.0 / 151.8 | 1.634 / 1.676 | 2.393 |
| degen3.dual-r | 266.7 / 255.6 | 163.3 / 144.2 | 1.633 / 1.773 | 2.381 |
| perold.dual-r | 149.4 / 164.1 | 72.5 / 71.0 | 2.062 / 2.312 | 2.873 |
| 80bau3b.dual-r | 317.9 / 319.4 | 145.8 / 140.2 | 2.181 / 2.278 | 2.907 |
| greenbea.dual-c | 715.8 / 701.4 | 413.6 / 374.6 | 1.731 / 1.872 | 2.596 |
| greenbea.dual-r | 696.6 / 687.3 | 358.9 / 370.2 | 1.941 / 1.857 | 2.908 |
| mode_switch.dual-r | 739.3 / 716.2 | 223.7 / 178.3 | 3.305 / 4.017 | 3.750 |
| primal_5959.dual-r | 148.5 / 131.9 | 105.0 / 105.5 | 1.414 / 1.251 | 1.402 |
| quad_2717.dual-r | 146.3 / 138.6 | 120.9 / 107.4 | 1.209 / 1.291 | 1.279 |
| fx_quad_obj_5371.dual-r | 61.3 / 55.5 | 43.6 / 45.9 | 1.406 / 1.211 | 1.371 |

Family mean ratios (run 1 / run 2 vs the table above): dual 1.865 / 1.965 vs 2.394; primal 3.220 / 3.197 vs 3.207;
barrier 2.420 / 2.264 vs 2.227; PDHG 1.202 / 1.144 vs 1.194. Primal, barrier and PDHG do not use the dual factor and are
unchanged within run-to-run noise (up to 28% on the smallest row). Rose is slower on 23/23 rows. Evidence:
[f4-stack7-run1.json](evidence/f4-stack7-run1.json), [f4-stack7-run2.json](evidence/f4-stack7-run2.json) (all 23 rows).

### Current optimized DEGEN2 dual stage accounting

Qualified fixed39 PGO/system-TCMalloc ELF55af6617 and reconstructed reference
26a50517, CPU2, dual-release. Central compare.py agrees on the covered637 pricing
state fields/decisions at zero ULP (inactive priority payload is unpaired), all15 post-setup factor boundaries/causes and final
vector bits. Both have634 pivots and16 factors including initial. The earlier
747 full FTRAN/634 needed/1094 BTRAN census counted corresponding/booked mirror
work, not every physical solve: current probe finds747 mirror full,634 needed,
1095 mirror BTRAN,5 auxiliary FTRAN and1 auxiliary BTRAN. Five auxiliary FTRANs
serve fresh_objective; both additional read-out BTRAN implementations belong to
termination. Both factor implementations update621 times. Reference update
wrapper637 and special helper88 are not successful-update counts.

Accepted solve medians45.784/16.696ms and solve-only181.193/53.198M instructions,
137.330/48.015M cycles remain separate from the diagnostic build. Healthy varied
accepted-ELF sampling supports48.026M instruction/19.616M cycle FTRAN/BTRAN self
excess and10.494M/7.857M factor-update self excess. The remaining pooled gap is
69.475M instructions/61.843M cycles. These are approximate self boundaries, not
inclusive per-stage/removable costs. Independent coverage93.40/94.25% instruction
and95.16/100.78% cycle limits precision; fixed-period/prime sensitivity and
throttled4M/frequency captures remain excluded. All512 retained profile and32
counter outputs agree. Evidence: degen2-stage-accounting.json under evidence/.

The exclusive-stage follow-up is complete without solver optimization. A plain
isolated release rebuild is byte-identical to55af6617. The probe keeps the same
fixed39 profile/TCMalloc, reconstructs all five modified original source files
literally after removing instrumentation, passes all17201 bounded DEGEN2 dump
records and full987 frozen parity with PMU instrumentation enabled. Delegated
central comparator cap20/640 and15 factor boundaries pass, with complete final
reference/accepted vectors bitwise equal including signed zero.

Eight shuffled CPU2 repeats per mode produce this diagnostic-only arithmetic
mean exclusive timer table (stage means sum49.415431ms):

| Stage | ms | Share | Allocation + growth calls |
|---|---:|---:|---:|
| Setup and initial factorization | 1.434 | 2.90% | 11032 |
| Pricing / ratio / pricing touch | 1.470 | 2.97% | 5310 |
| FTRAN / BTRAN | 17.257 | 34.92% | 30105 |
| Pivot products / state / controller | 11.216 | 22.70% | 27938 |
| Factor updates / subsequent factorization | 17.535 | 35.48% | 127777 |
| Termination / read-out / destruction | 0.224 | 0.45% | 557 |
| Unattributed loop/control/retirement | 0.279 | 0.57% | 8 |

Setup and termination override nested operation scopes; in-loop F/B/factor
scopes override pivot/controller parents. Initial reset belongs setup, later
pricing-weight reset belongs factor, selection/ratio/touch belongs pricing.
All allocation, copying, bookkeeping and publication/destruction stay in their
owning scopes; the residual includes pivot retirement outside guards. GlobalAlloc
counts146329 allocations/56398 growths/146321 frees,64.358MB summed requested
sizes (not live bytes, copied bytes or RSS). Factor owns63.03% allocation/growth
calls but its35.48% time includes required numerical/structural work, not a
removable allocation lifetime. Native-private allocation is not separately proven.

Plain/off/timer-noalloc/timer-allocation whole-process instruction skew is
0/+4.148/+6.203/+7.485%. Timer-allocation median50.162ms versus plain48.422ms;
run spread prevents an exact net wall-overhead claim. Empty timer transitions
cost52.8ns, about0.682ms for12915 transitions. PMU transitions cost4.422us,
roughly57ms; all PMU wall attribution is excluded (median116.350ms). Separate
user-counter stage means sum196.399M instructions/165.327M cycles with no
multiplexing, including instrumentation; do not transfer these to accepted ELF.
Stage-entry averages are descriptive mixed scopes, not homogeneous per-kernel
costs. Source/layout/allocator tracking and PGO-flow warnings remain explicit.

This closes diagnostic lifecycle attribution while leaving the accepted
cross-solver inclusive stage gap unmeasured. Factor/F/B/pivot lead; tiny setup,
termination or pricing alone cannot represent the broad residual in this build.
No arena/workspace/publication, eta, dense-L or row-scatter trial is reopened.
The14stage factor split below now separates update, construction, pricing-weight
reset and publication at unchanged lifetimes; existing rejected owner experiments
remain closed.
Evidence: [degen2-exclusive-stage-followup.json](evidence/degen2-exclusive-stage-followup.json),
raw isolated source/builds/commands/calibration/checks in target/perf-degen2-followup/.
Both accepted binary hashes and reference library are unchanged. The earlier
objcopy metadata rewrite was restored byte-for-byte before this follow-up.

### Current DEGEN2 triangular-kernel diagnosis

The current architecture and field/lifetime contract are in
[dual-factor-backend.md](../../reverse-engineer/specs/dual-factor-backend.md),
with provenance in [dual-factor-architecture.json](evidence/dual-factor-architecture.json).
A delegated bounded native SC50A capture now spans initial construction, all
six updates/seven pivot solves and next refactor's actual return (context7,
request10). It completes in2.1s/334events with no errors; elapsed is not performance
evidence. Central comparison verifies48 installed snapshots, L/LT immutability,
raw column/row coefficient transport, cache publication/consumption/retirement,
and seven source pivot projections. Initial50weight BTRANs,7rho and1true-d,
7FTRANpairs,3standaloneFTRANs and both constructors have paired returns.
Native keeps the same factor/pool/cache allocations; fresh views pack numeric
rows/L/U/LT into shared index/value pools. Fill58→142 resets/rebuilds to84,
capacity2778 unchanged. Cache pivots25/33/40 publish2/2/3 entries after U returns,
retain them until update and clear count/key atf377d/f3783. Slot/at alias existing
allocations at installation. Mark is F+e8=H+c8, not H+e8=F108 density; the spec
corrects the earlier base error. No monitored private allocation/free calls occur
inside this interval; outer allocation, exhaustion and ultimate teardown remain
unmeasured. This establishes a bounded native lifetime without implying complete
arithmetic identity or a speedup. No solver code/build/restart changed.

The separate pricing-reset internal comparison stops earlier at reset3/context82:
all516 term-role projections and444 weights per phase agree, but SUM differs
8ULP (native409ca37aef572848/Rose409ca37aef572850). Literal ADD/MUL order differs
from term0. Reset4 actual comparison was not continued; cached graph/final-weight
proofs remain separate. Whole15 Rose resets cost5.199606M instructions, including
3.509121M graph preparation and0.588160M cleanup; all later costs are source-only.
This is not an F/B kernel cost or a speedup bound. Evidence:
[degen2-reset-cost.json](evidence/degen2-reset-cost.json). No optimization.

The current bounded cached-evidence synthesis retains the accepted binaries and
cap20 complete result/exit/signed-zero control. Complete matched joint135 costs
81377/23001 instructions; sparse/dense rho BTRAN costs1226/504 and18179/14215.
Dense-L extra empty-node visits explain7018 of8956 walk instructions, with the
previous empty-gap trial closed. Valid complete inclusive cycles and exhaustive
per-purpose arithmetic/lifecycle equivalence remain missing. See the current
request synthesis in [degen2-fb-diagnostic.json](evidence/degen2-fb-diagnostic.json).

Strict central compare.py dual-kernels pairs162 pivot outputs bitwise, then
stops at needed pivot163/pricingk164/numericalit162: basis22 is native−0
(8000000000000000), Rose+0 (0000000000000000). Nativeimmediate123690 and final
compactslot0 agree on−0; fullcolumn values/support match. No comparisons after
this difference. Earlier637state comparator usednumeric equality for zeros;
that is not allinternalbit parity. Finalvector bit equality remains separate.
Do not fix or normalize the zero as part of performance work.

The native full/needed634pair census is192normal(162sparsefullU/30densefullU),
105sparsecolumn/densetau joint and337inlinejointdense. These whole-run counts
are structural, not strictnumeric equivalence after163. Native1094BTRANs split
634rho(sparseU;291sparseLT/343denseLT),444initialweights(sparse/sparse),16true-d
recomputes(dense/dense). Rose adds one mirror and one auxiliaryreadoutBTRAN,
plusfiveauxiliaryobjectiveFTRANs. Flip97/basic16 rawFTRANpairing and exhaustive
initial/recomputeBTRAN outputs remain unproved. Selectedrho1/20/115 bits pass.

Centralized nativepublication mappings explain firstposition-orderdifference115:
densefullU121ef0 returns−1 count, wrappercompactsascendingexcludingpivot then
appends pivotlast. Inlinejoint has its own provenpublicationblock. Needed
requestedposition order can differ; do not normalize numericalentryvisits away.
Selectedinput/Ustate/rowoperand proofs cover1/20/115 and pricing136/native135:
15rows29orderedterms and21rows25terms respectively. Dense neededU zero slots
and all L/eta/U visited-term histories remain unproved.

The concrete local excess is representation/traversal work. At jointpivot135,
Rose checks40numericUentries for8actualtail scatters, while native visits8A
coefficients;32nucleus-destination checks are discarded before multiplication.
Within-column destinationorder differs independently. Rose also rebuilds a
stamped numeric-row lookup and gathers virtual-row terms before needed dots;
native reads its storednumericrow directly. Both use scalar mul/sub onshort
rows and packedoperations onlongrows; no general arithmetic-speed claim follows.
Acceptedneededtail filter+scatter region has6.933M instruction/3.511Mcycle self,
including required arithmetic/skid, not a removable2.56% time share. Shorttau,
two-vector fusion, WORKdeferral, F2workspace/predicate, outputordering, denseL,
packedL and long-eta trials stay settled.

Existing acceptedprofiles reclassified by actualhelper owners give known
FTRAN/needed/sharedL/WORKpublication self66.597M instructions/39.052Mcycles
versus17.121M/18.356M native; BTRAN8.187M/6.810M versus5.801M/5.067M. Native
125f50 is the two-vector wrapper, not needed-only. Selfpool differences total
51.863M instructions/22.439Mcycles;76.132M/66.876M of whole-solve gap lies
outside these pools or remains unowned/sampling residual. They are not matched
inclusive per-kernel excesses. No allocation/helper/caller cleanup is assumedzero.

An isolated cost-only releaseprobe preserves originals modulo insertedguards,
passes17201original bounded dump records and full987 frozen parity. Samefixed39
profile/TCMalloc, plaincontrol byte-identical55af6617. SixshuffledCPU2 repeats
separate off/timer-noalloc/timer/PMU modes. Whole-process instructions skew
+3.229/+4.776/+5.113%; emptyinside/outer timers65.2/125.3ns, PMU4.557/9.100us.
AllPMUwalltimes excluded. Diagnosticinclusive rootmeans, not acceptedshares:
full7477.882ms;needed6346.896ms;bookedBTRAN10943.046ms;jointWORK6340.748ms.
Needed442joint-route calls take6.546ms,192sparse0.350ms. FiveauxFTRANs cost
0.064ms, extra mirror/auxreadoutBTRANs0.011/0.013ms individually. Do not infer
retained auxiliaryfactor construction/621updates/retirement cost from these.
Callee local allocation/copy/drop included; caller returned-vector retirement
and retainedfactor lifecycle remain outside and explicitly unassigned. TraceELF
is separate, selected-only, excluded from all performance claims.

Delegated central zero-origin comparison explains the first difference: native
visits a dense tail zero slot and divides +0 by -1; Rose's sparse traversal
omits basis22 and retains +0. The requested tau22 is consumed, but the observed
DSE subtraction erases its sign and both weight updates produce26.25 exactly.
Raw needed outputs and DSX/DTAU still differ; strict pairing still stops163.
This local consumer observation does not establish later-state equivalence.
Evidence: [degen2-zero-origin-and-search.json](evidence/degen2-zero-origin-and-search.json).

The completed joint-pivot135 transaction census establishes an additional
mechanism: cached forward L completes its numerical pass, rejects the dense
handoff with etas, then the general path rebuilds buckets and repeats that L.
Accepted-ELF runtime observations and delegated central comparison prove:

| Selected L solve | Native terms | Rose cached + fallback terms |
|---|---:|---:|
| Full | 46 | 46 + 46 |
| Needed | 132 | 132 + 132 |

Each ordered destination, coefficient, source and destination-before bit agrees;
all444 raw L-result slots agree with native logical materialization and both
Rose passes. This proves duplication in this pair, not general removal safety.
A bounded Rose census through162 pivots finds27 completed-prefix density
fallbacks:7 full/265 repeated terms and20 needed/1208 repeated terms. AFIRO3
and SC50A60 have none. No new capture crosses the first divergence163.

Accepted Rose full/needed roots make36/39 allocation or growth requests in the
selected pair, of which28/33 belong to L. Native's entire selected11f620 wrapper
makes zero observed exported allocator calls, reuses existing vector buffers and
clears two28-byte marker spans. These are dynamic counts at stated boundaries,
not complete instruction/cycle deltas or proof of zero private allocation.
TauRow and term scratch do not allocate in this call; they still scan25 numeric
entries, gather25 virtual entries and evaluate25 terms for21 nucleus rows.
The earlier needed-tail40checks/8scatters versus native8coefficients remains.

Construction and cleanup are kept with owners: needed's624-byte result retires
in its caller, while3552-byte x_tau and444-byte tau_raw remain inF2. Full's
1024-byte result transfers through pending_alpha into auxiliary sparse eta
storage, rather than being freed at return; its retirement belongs factor
lifecycle. L workspace/eta-index buffers and mirror spike also remain live.
No retained factor lifecycle is assigned by a mixed per-call average.

The completed rejection-cost pass uses target-pid PMU descriptors on unchanged
accepted55af; no rebuilt probe or codegen skew. A fixed-PC Rose-only whole-solve
census finds789 cached attempts and236 completed-prefix density rejections:
182 full and54 needed. Numerical pairing still stops163; these later observations
are instruction accounting only, with no numerical vectors or reference pairing.

| Rejected family | Calls | Cached attempt instructions | Fallback instructions before late cleanup |
|---|---:|---:|---:|
| Full FTRAN | 182 | 1.123M | 1.042M |
| Needed FTRAN | 54 | 0.456M | 0.452M |

Fallback includes x/start/bucket/active/output construction, repeated L, local
bucket/active destruction and WORK conversion. Original order/start/x destruction
after required eta/output adds0.0276M instructions; these frees are counted once
with their creating stage. Raw total1.522017M, approximately1.520129M after
calibrated probe instructions, is about0.84% of accepted181.193M solve instructions.
It is a loose opportunity ceiling: required WORK and handoff/control costs remain,
some cleanup remains necessary, and cached numerical work is needed to obtain the
result. Do not claim both L passes or their retained-cache construction removable.

At selected135 the full/needed fallback front costs6447/13530 raw instructions.
An existing3552-byte F2clone malloc+copy costs411 raw instructions; it is only a
handoff-copy proxy. Actual handoff support copying and cache behavior remain to
be tested. Correct late3552-byte x free is147707..14770d; earlier measured
1476af..1476b5 is32-byte eta scratch, not x. All236 late cleanup sets are observed.

Calibration measures7 instructions across a known5-instruction sequence. Detailed
root probes add16.125/19.875 instructions versus expected16/20, with no PMU
multiplexing. Added stops inflate full/needed cycles roughly48%/101%; all new
cycle/timer data are excluded. Whole3fixed instruction repeats differ by a few
instructions; all29 complete results and exits agree, retaining signed zero.
The dynamic-return whole-solve probes exceeded15s/45s deadlines and are excluded;
confirmed terminal before the simpler fixed-PC replacement, which finishes~1.1s.
No runtime optimization was implemented during this bounded investigation.

Existing F/B self pools already include this mechanism; do not add its cost again
to their51.863M instruction/22.439M cycle excess. The residual outside these pools
remains76.132M/66.876M. The earlier BTRAN prose6.743M was a transcription error;
raw6.810M is retained. No accepted time saving is established and the2.74x target
gap remains. This small mechanism cannot by itself close that gap.

The completed-prefix L handoff trial is now REJECTED. Isolated source freezes
complete original L and cached solver, with a reversible source proof. Original
and candidate19-test release runs pass, including captured135 full/needed,
SC50A, generated/extreme/tie/duplicate/nonfinite/incompatible/warm-state cases,
f64/F80 values/count/WORK/cache state and captured MXCSR results/stickyflags.
Final ownfixed39 PGO/systemTCMalloc build111c41c0 passes fast940 parity and all
17201 original bounded DEGEN2 dump records. All34 counter result/exit objects agree.

The predeclared CPU2 eight-pair DEGEN2 counter gate gives wholeprocess
instructions196.496M→195.394M(−0.5606%), cycles134.526M→139.909M(+4.0016%).
Individual cycle pairs have both signs; this is a failed descriptive counter gate,
not a causal microarchitectural proof or formal timing claim. AFIRO/SC50A/ADLITTLE
instructions+0.008/+0.070/+0.093%, cycles−6.078/+0.475/+4.895%. Stop without
reserved74117/74189, retry/extension, PGO tuning, enum/flag/layout or buffer-moving
rescue. No full candidate gates or production source promotion. Candidate
1329a7ff, source/frozen oracle/tests/capture/own profile and every result remain
isolated. The numerical feasibility finding is retained; no speedup is credited.
Evidence: [degen2-l-handoff.json](evidence/degen2-l-handoff.json).

Existing F/B/cost evidence: [degen2-l-rejection-cost.json](evidence/degen2-l-rejection-cost.json),
[degen2-joint-transaction.json](evidence/degen2-joint-transaction.json),
[degen2-triangular-kernel-diagnosis.json](evidence/degen2-triangular-kernel-diagnosis.json).
The bounded [F/B diagnostic synthesis](evidence/degen2-fb-diagnostic.json) separates
purpose/path counts, diagnostic inclusive roots and accepted sampled self pools.
Its current complete-boundary synthesis incorporates the transaction costs below:
the largest selected localization is Rose L+eta56868I versus a deliberately wider
native input/L/inline-eta/dispatch prefix12211I. This does not isolate removable
arithmetic or support construction. BTRAN still lacks matched complete per-purpose
root/helper/cleanup costs and exhaustive visited-entry histories. No new profiling
or solver instrumentation is needed for the synthesis; an accepted cap20 result
again matches cache including zero signs, with accepted binary hashes unchanged.
The L+eta instruction partition is now complete on unchanged accepted Rose.
Whole747full/634needed roots sum25.387227M instructions,210I above the earlier
root census. Dense L numeric walk9.344032M across592calls dominates; dense
support recount.870241M, dense-spike construction/publication4.130082M across828,
general eta gather1.767646M/dot2.065352M across14887etas, cached numeric/scheduling
1.718040M across789attempts and repeated fallback numeric.991367M across236.
Successful cached suffix(eta/output included)1.219813M and other input/prepare/
WORK/dispatch/cleanup3.280654M close the disjoint ledger exactly. Costs include
required operations and are not time/removal shares. No native pairing after163.

At selected135, cache preparation/reset5846I, cached traversal18763I, fallback
traversal15574I and dense spike helper9932I(two4966Icalls) explain most56868I.
Each general eta pass has one4term eta, gather99I/dot33I; eta arithmetic is tiny
in this pair. Queue calls cost6460I including helpers but whole accepted sampled
queue self is only.398M I/.533M cycles(.22%solveI), so no standalone queue trial.
Selected arithmetic/input/support pairing still uses central compare.py. Root/
detail residual0–2I, whole correction115350I for57675intervals; all debugger
cycles/wall attribution excluded. Sixteen complete result/exit/signedzero controls
pass, small AFIRO/SC50A first; both accepted binaries/reference remain unchanged.
See [degen2-l-stage-partition.json](evidence/degen2-l-stage-partition.json).
Delegated native135 final fixed-PC partitions pass7central boundary gates and
reconcile23001I wrapper/12211I pre-U exactly. Native full/needed L=1985/4806I:
setup/reset/input185/320, numeric/node-support/output1776/4462, WORK/return24each.
Post-L5132I comprises denseclear858, scatter/setup810, both eta dots+metadata/
WORK127, spike/admission3300 and dispatch37. One444slot spike scan/publication
costs3285I; the two Rose dense-spike helpers cost9932I total, with different
necessary publication responsibilities. Do not treat the difference as removable.
Initial native nested-entry probes had a2I deficit and remain excluded; final
probes omit those stops and reconcile without ad hoc correction.

The early sparse135 pair is unrepresentative of the largest whole exposure:
9.111M of9.344M dense numeric instructions occur afterpricing163, unpaired.
The dense witness source149/native148 is now measured on accepted binaries:
full14/needed389ordered coefficient/source/destination-before/product/result
projections agree centrally. Needed444rawslots before/input/L/aftereta agree, includingzerosigns; fullsparse
afterL uses an explicit logicalembedding, not native444outputmemory. Fullaftereta
19publishedslots/order/bits and sourceworkspace atthosepositions agree; native
unpublishedslots remainunpaired. All16source result/exit/zerosigncontrols and
4native centralcostcaptures pass, root/detail both14018I. Visitation
differs: native dense scans112positions(0..106,439..443), Rose444. The332extra
columns107..438 all haveempty Lspans;202inactive/130active visits cost7018–7019I
in fixed-PC brackets, matching9Iinactive/40Iactive staticpaths. These comprise
78.36% of the selected densewalk excess14435−5479=8956I. Remaining1938I mixes
index/bounds/unrolling/WORK/codegen, without causalcycle attribution.

Complete source149 F/B transaction72460I versusnative14814018I includes owned
transient cleanup and excludescallercheck: roots67988/13547 difference54441,
outside4472/471 difference4001, total58442I. Native fullL1131/denseL6205I;
denseL contains numeric/node/rowWORK5479I, setup/scatter/clear718I, return8I.
Native preU9480I and source L+eta31099I have deliberatelydifferent widths;
not an exact arithmetic-only delta. Retainedfactor lifetimes stay separate.

WholeRose-only census over592densecalls finds119415emptygapvisits/41338active,
16geometries with everygap empty and allLdestinations strictlyforward. No gapNaNs
observed; this does not remove the preservation requirement for nonfinites.
Static originalpath exposure2.356201M I is1.30038%solveI inside9.344M densewalk,
not a removable/timebound. New directdense9.344025M differsprior by7I.
See [degen2-dense-L-diagnosis.json](evidence/degen2-dense-L-diagnosis.json).

The isolated empty-middle specialization is now INCONCLUSIVE and closed. It
retains original gapthreshold tests and +0WORK additions in roworder, avoiding
knownempty Lspan lookup/check/conversion; invalid/nonempty gaps use originalwalk.
Complete original Lsource is independently frozen, source reversal exact. The
unchanged harness passes19focusedrelease tests; expanded candidate and ownprofile
pass20, including captured149/SC50A, generated f64/F80 extremes/duplicates/gap
boundaries, exact values/order/count/WORK/cache/snapshots and16MXCSR modes for
both formats. All661ordered immutable Lcoefficient/index bits from cached135
agree with149. Ownfixed39 PGO/systemTCMalloc profilebf1ee89a, ELF8ead3797, passes
fast940 and all17201originalDEGEN2dump records. Delegated central candidateprefix
162PASS/first163expected signedzeroFAIL unchanged; selected149outputs/rho/context
PASS and164extendedDSXrecords match acceptedprefix. Candidateinnertrace unclaimed.

Predeclared CPU2countergate passes: targetwholeprocess196.501M→195.095M I
(−.715457%),145.989M→144.137M cycles(−1.268227%). Small AFIRO/SC50A/ADLITTLE
instructions+.0031/−.0386/−.1071%, cycles−4.303/−2.144/+.921%; all34complete
result/exitpairs agree. These counters are descriptive, not solver-time proof.
Screen83117(32AB+32AA) is INCONCLUSIVE: median−2.92064%, CI[−6.73678,+4.99806]%,
18/32,p.2983; AA CI[−6.15317,+6.77304]%,p.3771, bothorders improve. Sign and
conservative AA-separation gates fail. Stop without83189, extension, tuning,
loop/helper/threshold/layout/PGO/buffer rescue orpromotion. Full987/allrelease/
fullpanel not run afterfailedscreen. Acceptedsource/binaries/reference unchanged.
This closes the tested empty-gap design despite exact numerical feasibility and
instruction improvement. Evidence: [degen2-empty-L-gap.json](evidence/degen2-empty-L-gap.json).

The bounded accepted-ELF U-cost diagnostic now separates full33556I and needed
34432–34433I at source149/native148. Full L+eta7819I leaves25714I after L;
needed L+eta23280I leaves11133–11134I. Source jointWORK1693I is already inside
outside-root4472I. The complete72460/14018I transaction difference58442I stays
reconciled; native jointU/publication4050I cannot be split arbitrarily by purpose.
Source full output-order construction costs16622I including malloc6/calloc3/
realloc11/free8, with5012requested allocation/growth bytes (not traffic). Actual
19inputpositions (11singleton/7tail/1nucleus) expand to92ordered candidates and
23nonzero outputs including pivot138; ordered filtering/buffer retirement2857I
and leftover touched scan81I are separate. Full F2publication153I/spike206I and
finalcleanup103I remain included. Needed post-L input/tail construction3914I,
tail3469I, individualWORK520I, F2clone/publication440I, requestedrows2612I,
publication9I andcleanup169–170I are disjoint; retained factor lifetimes remain
separate. These are complete scopes, not arithmetic-only or removable costs.

Two fine repeats each AFIROcap3/SC50Acap3/DEGEN2cap150 and the separate numerical
order observation preserve complete results/exits/signedzeros. Fixed-PC5I→7I
calibration subtracts2I per interval; full/needed root versus fine residual0–1I.
All debuggercycles/times excluded. One initial fine capture used an invalid
instruction-boundary PC and stopped the disposable inferior; it is retained and
excluded, active PCs corrected/validated. Accepted on-disk binaries unchanged.
The large concrete mechanisms are extra Lcolumn visits, support reconstruction,
and output ordering/publication; no broad slower-f64-codegen claim is supported.
Central compare.py dual-u-boundary passes cachedeta/nativeactualUentry and
source19positions/kinds plus published23orderedrawbits; wrapperremoval yields22.
Native153full/39needed actualdivisions include130/3zeroquotients. Actual source
Uoperation projections are now centrally checked below; global operation order
and exact gathered-row column identities remain separate from endpoint equality.
The selected BTRAN boundary follow-up is complete on unchanged55af/26a505.
Firstdense-LT native148/sourcepricing149 (it147,r138,upd16) is inside the
strict prefix; all147earlier rho routes use sparseLT. Complete matching caller
transactions, including row-length bookkeeping and Rose deferred raw/filtered
Vec cleanup, cost source1226/native504I sparse1 and18179/14215I dense149.
Roots alone are1022/457I and15958/10276I; native dense postfilter costs3923I,
so root-only publication excess4080I overstates the complete3964I difference.
Rose dense packing6297I includes malloc7104B and444x14I loop; separate filter
1701I/rowlength408I/WORK21I/twofrees80I remain owned. Whole360pack calls
cost2.266946M I(1.251%solveI):343rho2.159875M,16true-d.100774M,one mirror
readout.006297M. This includes required publication, not removable/time share.
Central compare.py dual-btran-boundary passes raw input/head through accepted
DSX, U10entry publication, dense444prepared/afterLT/output and444orderednode
value bits. Both walk444nodes/65terms, but ordered term comparison stops at
first row441 destination native98/Rose15. No65ordered-arithmetic parity claim.
The observer's bad RefStruct-relative head field is retained unpaired; cached
accepted DSX supplies actual Dual basis state. No solver difference normalized.
All17source full result/exit/zero-sign controls pass. Coarse/detail/trace roots
agree; per-node source+1I/LT+2I and native dense detail+2I are explicit uncertainty.
Calibration5I→7I; debuggercycles inflate materially and all new cycle/time claims
are excluded. Retained/private factor lifetime and other BTRAN families remain
unmeasured. The U/eta follow-up now accounts for the selected pre-L difference:
Rose2133I U/setup +1345I post-U/eta versus native875+981I, difference1622I.
Native981I also owns cachedUrow materialization/inputconversion/WORK/cleanup,
so this is a wider localization, not a pure corresponding eta-kernel delta.
Source18occupied-entry calls486I and three-entry input slot rescan43I include
required addressing/bitset semantics. Exclusive tail1684I partitions into
schedule415/divide-threshold45/output-WORK521/scatter621/advance72/setup10;
eta903I scan splits check102/WORK90/term345/support185/advance170/setup11.
Central compare.py dual-btran-u-eta stops at first Uscatter: after matched
unitDIV138, native targets419/Rose22. Later16nativeDIVvs10Rose (six extra
nativezeroDIV),18products/updates each and8eta/6active/23MUL17ADD are structural
counts only, no later ordered arithmetic pairing. First literal native update
uses MUL/signXOR/ADD, Rose MUL/SUB; no arithmetic-throughput cause follows.
Whole1095mirror roots reconcile11.924293M I: U2.196653M, eta1.868727M,
LT/publication/cleanup7.858913M. Rho634 U+eta3.658169M is2.019%solveI;
all-purpose U+eta4.065380M is2.244%, including required work, not removable/time.
All17source full result/exit/zero controls pass. Coarse/fine/helper/partition
U2133/eta1345 agree; numerical trace undercounts46/17I and has no cost claims.
Native initial observer with nine invalid after-PCs is excluded; corrected wide
instruction decode validates all PCs and repeats small qualifications. No solver
optimization/new build or speedup. See
[degen2-btran-u-eta.json](evidence/degen2-btran-u-eta.json). No further local
BTRAN prototype follows from this scope; return to complete non-empty constructor
accounting with source-only late costs and small prefrontier native witnesses.
Evidence: [degen2-btran-boundary.json](evidence/degen2-btran-boundary.json).
See [degen2-U-diagnosis.json](evidence/degen2-U-diagnosis.json).

The unchanged-ELF output-order/lifetime census is complete across747innerfull
calls (634pivot/97flip/16basic),364sparse/383dense. Both compiled branches have
probes; all observed no-dump calls use clone2. Disjoint instruction scopes:

| Inner full-FTRAN scope | Calls | Instructions M |
|---|---:|---:|
| Entry, L+eta and pre-L setup | 747 | 9.982448 |
| Post-L input/tail/numeric/individualWORK before output | 747 | 10.208631 |
| Sparse output-order construction and local helper cleanup | 364 | 4.315376 |
| Sparse filter/result growth/order-buffer retirement | 364 | 0.848464 |
| Sparse remaining-support scan | 364 | 0.023377 |
| Dense result scan/growth | 383 | 2.136981 |
| Output flags/F2 publication/prior-buffer retirement | 747 | 0.120048 |
| Spike publication/prior-spike retirement | 747 | 0.330844 |
| Return publication/local cleanup/epilogue | 747 | 0.069661 |
| Total | 747 | 28.035830 |

Order construction is2.381646% of accepted181.193M solveinstructions; even
free deletion would not bridge the gap. Required ordering and replacement costs
remain, so this is not a removable/time share. Actual4728inputs expand to21036
orderedcandidates and10013sparse outputs. Order scope observes1092calloc,
1766malloc,1987realloc,2494free;364returnedorderbuffers retire in filter.
Requested1.285928MB allocation/growth,543192B memcpy and609984B memset include
nested allocator requests; not physicaltraffic/livebytes. Costs include helpers.
Output/publication/localcleanup sum7.844751M I. CallerCALL/RET, keep/basicwrapper
spike preservation/restore, eventual returnedcolumn auxiliaryEta retirement and
final retainedfactor destruction are outside this explicitly inner boundary.

Root-only28.035779M versus detailed28.035830M residual51I,94calls differ and
maximum4I;13,408I stopped-interval correction applied. Small AFIRO3/SC50A3/
DEGEN2cap20 root/detail totals agree exactly. All13complete result/exit/zerosign
controls pass; longest detailed4.262s, separate memoryobservation8.749s, neither
is timing evidence. No cycle/codegen/solver/binary changes. Ordercost0.977023M
throughcoveredsourcepricing163 and3.338353M later; later costs are source-only,
not matched native excess. Existing central gate reused, no new numericalpairing.
See [degen2-output-order-census.json](evidence/degen2-output-order-census.json).
Actual Uoperation diagnosis at matched149/native148 is complete. Source fresh
full19logicalinputbits/needed444rawinputbits and kind arrays pass delegated
central compare.py dual-u-operations, reusing the existing nativecapture:

| Actual division family | Rose | Native | Native-only zero DIV |
|---|---:|---:|---:|
| Full singleton | 11 | 11 | 0 |
| Full tail | 11 | 20 | 9 |
| Full nucleus | 11 | 122 | 111 |
| Needed singleton | 11 | 11 | 0 |
| Needed tail | 19 | 19 | 0 |
| Needed nucleus | 9 | 9 | 0 |

All72corresponding numerator/diagonal/quotientbits pass. All20full+32needed
scatter from-column/destination coefficient/source/before/product/after tuples
pass; native ADD/signflip actualoperands remain recorded separately. All21ordered
rowdot coefficient/value/accumulator/product/result projections pass for9nucleus
rows. Gathered sourcepairs lack observedcolumnIDs, so exactrowdotentryidentity
is unproved. Globaldivision/scatter interleaving differs, retained literally;
this is not identical instructionstream parity. An initial comparator field-label
mistake for tuple(coefficient,value) was corrected from actualsource/disassembly,
not a numerical difference. No normalization or later-pivot comparison.

Needed tail visits237numericUCcandidates for32scatters;205nucleus candidates
are rejected before coefficientload/MUL. Native visits32relevant coefficient
entries in presegmented[start,split) storage (12d510/12d518,neededloop12d72f),
without per-entry kindfiltering; full-only suffix[split,end) is separate. Partition
construction/maintenance and other support/indexwork are not free or accounted
by this selected scatterloop. Accepted fixed-PC measurement gives205×7Irejectcore=1435I exactly,
205×3Iadvance/static615I versus616Imeasured:2050static/2051measured instructions.
This is59.095% of selected3469I tailroot. Detailed3470I versusroot3469I residual1I;
496marks/990I stopped-interval correction, allcycles excluded. Remaining1419I
includes actualnumeric/index/control/WORK responsibilities, not a reference gap.
New2050I is already within34432Ineeded/72460Icomplete transaction, not additive.
Rose executes120fewer full zeroDIVs, so extra Uarithmetic is not this gapmechanism.
Tailentries are noncontiguous inUC (column315 offsets20/22/33of34), so no free
contiguous-slice shortcut is proved; indexedmetadata lifecycle costs are required.

Separate rawtrace AFIRO3/SC50A3/DEGEN2cap150 and fourcostcontrols preserve all7
complete results/exits/zerosigns. Small controls do not exercise selectedUCloop.
Trace solverclock6.294s versusplain.005735s (~1097x) is distorted/excluded;
actual register/memory bits only, no modeled floating arithmetic. Acceptedsource,
binaries, WORK/dumps/frozenparity unchanged. No new referencecapture or optimizer.
See [degen2-U-operation-diagnosis.json](evidence/degen2-U-operation-diagnosis.json).
Whole neededU census now covers634calls/44322executednonzero-column visits:
339979UCcandidates,161569tail scatters,178410rejectednucleus entries. Read-only
bulk index/kind observation at14360d avoids per-entry traps; independent count/
geometry modes agree percall andselected149reproduces237/32/205. At10I/reject,
whole exposure1.784100M I is0.984641% of acceptedsolveinstructions, before any
replacementcost. Coveredprefix162calls gives31590I; later472calls1.752510M I
is source-only exposure, not matchednative excess. No debuggercycles/timeclaims.

Lifecycle capture includes637actualstructuralstates:16constructors,621updates,
and finalfactor16 with no neededUcalls. Targetedfreshreturncapture includes its
217tailcolumns/929UCentries ratherthan omitting the unusedconstructor. All11
complete result/exit/zerosigncontrols pass; longest count13.139s/geometry13.916s,
finalcapture.670s. No solver/numeric/source/build changes or newnativecapture.

Standalone integer-only indexlifecycle replay is separate non-PGO release/
Rust1.98.1/systemTCMalloc, not the qualifiedsolver. All637stateindexvectors verify,
includingdeletedcolumns and16factor retirements. Three balanced CPU2rounds of64
lifetimes each, controls subtract common fixtureparse/eventdispatch:

| Index lifecycle policy | Columns refreshed | UC entries scanned | Indices written | Instructions M | Cycles M |
|---|---:|---:|---:|---:|---:|
| Eager rebuild | 59043 | 463057 | 202766 | 7.936083 | 12.079160 |
| Changed columns, free invalidations | 3425 | 41297 | 22830 | 0.937037 | 1.299930 |

The replay includes header/buffer construction,refresh andcleanup, retains
capacity and supplies exactoutputcounts/dirtysets free. Actualdirty detection,
updatehooks andnumericlookupindirection over161569retainedscatters are missing;
fixtureheapstate/codegen differ fromsolver. No subtraction fromacceptedcycles/
time or verifiedsolvergain. Eager testeddesign is unpromising. Optimisticdirty
budget arithmetic leaves.847M I(.4675%solveI) before missingcosts; crossbuild
illustration only, not a measurednetbenefit or proofagainst all representations.
No solverprototype/rescue follows thisdiagnostic. See
[degen2-U-candidate-census.json](evidence/degen2-U-candidate-census.json).
The same-pivot U/eliminate followup now passes all621captured sourceupdates;
strict native pairing remains bounded at162. Frozen guard studies and complete
component lifecycle cost are recorded below. No unconditional reuse or speedup
claim; no failedscreenretry.
Natural-output reuse,slot heaps,WORKdeferral,F2workspace/owner/arenas andemptyLgap
trials stayclosed. No optimization; the accepted2.74xgap remains.
The [complete transaction cost evidence](evidence/degen2-fb-transaction-cost.json)
now measures unchanged accepted ELFs with target-pid instruction counters. Matched
full+needed transactions at source pricing1/20/115/136 cost Rose13020/8592/43495/81377
instructions versus native1785/1696/18138/23001. Source136 maps native135. Disjoint
Rose scopes exclude caller pivotcheck and intervening state updates, include input
construction, helpers, WORK, publication and deferred transient cleanup. Retained
Lcache/F2/spike/auxEta lifetimes stay outside, not assumed free. No new cycle/time
claim: debugger distortion excludes both. Native30 central comparisons and all24
Rose result/exit/signed-zero controls pass; strict pairing still stops163.

Selected135 root excess is54197 instructions, preparation/WORK/publication/cleanup
excess4179, totaling58376. Rose L+eta roots56868 versus native's deliberately wider
input/L/inline-eta/dispatch prefix12211 localizes most selected excess. This is not
an exact L-only delta or exhaustive eta arithmetic equivalence. Native21 complete
numeric-row scopes cost996I; Rose prepare/gather/dot portions3196I omit dispatch/
output control, so pure arithmetic instruction differences remain unmeasured.

Whole-Rose fixed-PC census finds747 full L+eta roots9.966M instructions and634
needed roots15.421M, combined25.387M. Needed root total32.788M separates numeric-row
prepare1.609M, gather2.350M, dot/divide1.114M and other12.293M after L+eta. Numeric
row preparation/gather3.959M is2.185% of whole-solve instructions including required
work, not time or removable cost. Of19765 numeric rows,15665 occur in value mode;
do not extrapolate early matched135 to this later unpaired phase. The whole-census
instruction correction is198488I; selected detailed/root controls agree within1I,
while extra-detailed native wrappers retain1–4I residual. Adaptive per-output
breakpoints understated the selected root by20I; those root costs are excluded,
and fixed output breakpoints reproduce43354I before the whole census. Allocator-owner traces
are excluded from all instruction costs; fixed-PC cleanup includes guard/free costs.
The repeated-L fallback ceiling remains0.84% of solve instructions and its handoff
trial is rejected. No rowmap/workspace/owner or handoff trial is reopened. New costs
are inside existing self/whole pools, never additive; whole solve gap stays2.74x.
All accepted binary/reference hashes remain unchanged. Whole target gap2.74x
remains; the active goal is not achieved.

### Current DEGEN2 factor-stage split

The completed14stage probe is a separate original-source/fixed39a018/systemTCMalloc
ELFaad3107e. Reversing its four scoped edits literally reproduces the prior
verified7stage source; all40 complete results/exits and17201 dump records agree.
The new14stage probe has not newly run full987; the prior7stage source passed it.
These are diagnostic measurements, not accepted-binary or candidate cost shares.
Initial factorization belongs setup; the table counts15 later constructions and
621 updates. Eight shuffled CPU2 repeats give exclusive means:

| Factor substage | Entries | ms | Allocation + growth |
|---|---:|---:|---:|
| Ordered/numeric factor construction | 15 | 4.424 | 37744 |
| Main mirror construction | 15 | 2.704 | 25762 |
| Auxiliary factor construction | 15 | 2.497 | 8227 |
| Pricing weight reset | 15 | 1.394 | 17070 |
| Main mirror update | 621 | 4.669 | 22595 |
| Auxiliary update | 621 | 0.412 | 72 |
| Outgoing main mirror retirement | 15 | 0.292 | 0 |
| Other parent bookkeeping/local and auxiliary retirement | 15 | 1.435 | 16307 |
| Total | mixed | 17.826 | 127777 |

Allocation/growth total exactly reconciles the prior broad factor stage.
Returned Structure/mirror/auxiliary objects retain their owners; callee-local
cleanup is included, outgoing mirror drop is separate, and outgoing auxiliary
factor/local Structure drop remain in other. No lifecycle cost is inferred from
five objective FTRAN averages. The9.625ms construction sum and4.669ms mirror
update are not removable bounds.

Whole-process instruction skew off/timers-noalloc/timers is+5.024/+7.164/+8.447%.
All14 timer stage means sum49.540500ms versus own solve-clock49.539508ms; the~1us
edge difference remains explicit. PMU stage totals include probe/read/codegen
costs; PMU wallmean119.280ms is excluded. No accepted per-substage cross-solver
gap follows from these values.

Native d4190(caller f34f2, singleton/pair helper df800) produces current numeric
CSR/diag/L, while update headers p70/p78/kind/nuc at constructor return mix old
state. Complete installation occurs later in parent f15c0, after pricing reset;
do not equate constructor return alone with the complete Rose fresh root.
Before strict divergence163, cached factor events are pricing0/16/37/82/132;
132mapsRoseiteration131/request4. Rose's separate left-looking auxiliaryfactor
has no separately matched native auxiliary constructor. Reset maps132460/caller
f48fd; full16count is structural only. No new comparison crosses divergence.

The first nontrivial factor16 census is complete. Delegated centralcompare
matches444 basis/pivot positions, all584 scaledinput coefficientbits,444diag
bits,140ordered L entries and emptyU. This factor is20phase1/424phase2 with no
Markowitz nucleus or fresh tail. Sourceprobe reverses literally to original
numeric sources, preserves1002original cap20dumprecords and complete result/exit.
Only registered isolatedFC tags are excluded; no baseline is re-recorded.

| Selected operation family | Rose | Native |
|---|---:|---:|
| Logical input coefficients | 584 | 584 |
| Main L coefficient divisions | 140 | 140 |
| Main fresh-U coefficient products | 0 | 0 |
| Auxiliary elimination products / L divisions | 2 / 6 | no matched auxiliary constructor |

Rose performs three584-entry singletonfilters plus510pivotlookup visits,
materializes584integer row-adjacency indices after584count visits, maps140L,
transposes140LT andvisits140Llevel edges. Native scans160structural indices twice,
copies28eligible numericrowadjacency entries, stores20explicit/424implicitdiagonals,
and divides140Lcoefficients. Integerrowlists and numericadjacency are distinct
families: their counts cannot be subtracted as duplicate numeric materialization.
Output equality and mapped operand bits do not prove division execution-order
identity. No substantial extra numerical elimination pass is established.

The independent Rose auxiliary factor orders the basis differently: L6/U136,
584inputadds,136heap pushes/pops and444 touched sorts over586total entries. It is
additional work with a different representation, not corresponding repeated
main-factor arithmetic. Diagnostic allocations/growth ordered139+emptyNucleus12,
fresh126,aux60 stay separate from accepted costs; native0exportedallocatorcalls
inside d4190 exclude persistent workspace/driver lifetime. Returned objects and
later retirement remain with their owning scopes.

Original accepted-ELF target-pid instruction measurements rule out empty-nucleus
setup as a substantial lead: threeempty calls total39247rawI(~0.0217%solveI),
including their construction and cleanup. All16nucleus calls total7.294644MI,
including required later arithmetic. No new empty-case or Markowitz/arena trial
follows; no debug cycle/time saving is inferred.

A distinct integer-ordering cost is now measured: auxiliary col_order stable
sort by(column_length,incoming_position), c852a..c853e, takes1.025200M rawI over16
builds(~0.566%accepted181.193MI). It includes sorthelper allocation/cleanup;
original identityvector creation and replacement histogram/scatter/cleanup are
outside the opportunity. Target-pid calibration implies~2extraI per bracket,
no multiplexing; all debuggercycles/timers are excluded. SC50A qualifies first;
all seven source/counter result checks match. WholeRose counter runs do not
extend numerical pairing past163.

The bounded length-bucket ordering trial is now INCONCLUSIVE. It freezes the
original integer keysort, preserves initialidentityvector and produces the exact
(length,incomingposition) order with a bounded histogram/prefix/scatter; small
inputs/widekeyranges keep originalsort. Numericfactor body literally reconstructs
unchanged. Sixfocused release tests pass, includingcapturedfactor16permutation,
random/tie/empty/boundary/fallback orders andexistingLU tests. Ownfixed39PGO/
TCMalloc passes fast940,17201originalbounded dumps andall34counterresults/exits.

CPU2 DEGEN2 counters improve−0.3951% instructions/−1.2937% cycles. The first
predeclared32AB+32AA screen75117 is inconclusive: median−2.1457%, interval
[−4.9610%,+2.3988%],17/32improving,p.4300; AB−4.6035%,BA+2.0725%. AA is healthy
(p.8601) but interval[−4.5697%,+4.3121%] prevents conservative separation.
All128measured+2warmup outputs/exits andbegin/end hashes agree. Stop without75189,
extension/retry/key/threshold/histogram/radix/allocator/codegen rescue orfull
candidategates. Candidateed4051b7/source/profile/tests/results stayisolated;
accepted source/binaries unchanged, no gain credited. This settlesauxiliary
columnordering; rejectedtrue_d coefficientpresorting remains separatelyclosed.
Evidence: [degen2-aux-length-order.json](evidence/degen2-aux-length-order.json).

The accepted-ELF complete constructor accounting now reconciles a42.509657M I
root union:15refactor parents plus initial reference-order/fresh/aux roots.
Initial outer setup and final retained-factor retirement remain outside this union.
Disjoint owners (initial construction included; later reset/retirement only):

| Owner | Calls | M instructions |
|---|---:|---:|
| Basis input preparation/order publication | 16 | 1.297683 |
| Ordered structure outside nucleus | 16 | 7.180906 |
| Numeric nucleus including local construction/cleanup | 16 | 7.312739 |
| Main mirror fresh | 16 | 8.688487 |
| Auxiliary LU | 16 | 8.480433 |
| Pricing reset | 15 | 5.200393 |
| Old mirror / old auxiliary / local Structure drop | 15 each | 1.877117 |
| Factor WORK | 15 | 0.767629 |
| Attempt tuple/caller publication | 15 | 0.007622 |
| Refactor other/quality/copy/local bookkeeping | 15 | 1.696648 |

Ordered-only partition14.493645M includes pre-nucleus structural/input4.572M,
nucleus7.313M, singleton L materialization2.376M and other.232M. These include
required numerics, indexing, allocation/copy/cleanup, not removable/time shares.
Current nested nucleus7.312739M is18,131I above minimal fixed-PC7.294608M;
method/observer/private-path sensitivity is unassigned (~.010%whole solveI).
Use the internally closed nested ledger; do not mix its parents with minimalchildren.
All21source complete result/exit/zero controls pass. Sixteen5I sequences read7;
coarse/detail parent totals vary by a few hundredI, native root/detail by0–1I.
All debugger cycles/timers are excluded. One incomplete leaf RET observer is
excluded; fixed CALL/return scopes include normal tailcalled cleanup paths.

Central compare.py extends the earlier empty-nucleus proof: source82/native82
has961scaledinputbits,L491/U25; source132/native131 has1138bits,L661/U32.
All444diag,basis/pivots/thresholds,orderedL/logicalU pass, with proven physical
U reversal transport, not internal arithmetic-order identity. At actual parent
return, installed i3/i21/kind/slot/diag match first source FSTRUCT of sameepoch;
full A/B/rowl/transpose/cache/pool identity remains unproved. Complete parent
instructions source1,259,082/native178,471 and1,412,585/198,156 show selected
excess1,080,611/1,214,429I. Local helpers and cleanup included, retained lifetime
outside return unmeasured; tiny3/5-column nuclei do not represent late whole costs.
Selected pricing-reset roots223,933/265,617I versus native17,272/19,291I warrant
actual geometry/recurrence accounting before any candidate. The ptype2 discarded
recurrence hypothesis does not explain this gap: native also executes recurrences
and sum before its type check. The unused growth-test L copy is the previously
rejected direct-diagonal-read trial and stays closed. No optimization/new build.
Evidence: [degen2-construction-accounting.json](evidence/degen2-construction-accounting.json).

Factorcensus evidence: [degen2-factor-census.json](evidence/degen2-factor-census.json),
[degen2-factor-split.json](evidence/degen2-factor-split.json). No new arithmetic,
WORK or dump/baseline change has been promoted; target2.74x gap remains.

### Current DEGEN2 pivot/state accounting

The isolated22stage source reconstructs the14stage numerical source literally.
Same originalfixed39 PGO/systemTCMalloc, separateELF9de5c7c4. ThirtyshuffledCPU2
runs(6permode) and17201originalDEGEN2dumprecords agree; releaseprobe with timers
passes all987frozenparityruns(351.4s). No solver/baseline/WORK change or promotion.

| Exclusive diagnostic purpose | Mean ms |
|---|---:|
| Row products and WORK | 3.857 |
| Input/norm/count/check | 1.498 |
| Basis/rrows | 0.420 |
| Dual costs | 0.460 |
| DSE weights | 0.143 |
| Basic values/flags/priorities | 1.550 |
| Recomputation outside numerical children | 1.421 |
| Object construction/driver/deferred cleanup | 0.672 |
| Remaining ROW control/local retirement | 1.878 |

Sum11.898ms includes object/driver work formerly outside14stageROW; do not subtract
old11.042ms as ifsameELF/scope. All22stages sum49.909047ms versus ownclock49.908257ms
(edge0.790us); phasebin sums reconcile exactly. Allocations146329,growth56398,
frees146321,requested64.358MB match oldcensus. Deferredoriginalinputlocals still
drop underROWresidual; objectretirement isexplicit, notsilentlyfree or charged
as completeinputlifecycle. Returnedfullcolumn remainsauxEta/factor-owned.

Off/timer-noalloc/timer-alloc wholeprocess instruction skew+5.793/+11.214/+13.773%.
Emptyinside54.5ns/outer113.5ns;13543entries imply~1.537ms emptytransitionbudget,
not automatic subtraction. DSE0.143ms has~.0346ms inside-equivalentclockcost and
isoverheadsensitive. PMUouter9.285us andmean186.933ms invalidateallPMUwallattribution.
ProbeI/C includeinstrumentation/codegen andstayseparatefromaccepted55afcosts.

Accepted55af fixedPCtarget-pid counters measure634rowproduct+WORK roots12.610M
instructions (10.857M in356valuephasecalls),634basis/rrows roots0.799M,621fast
DSE/basic/WORK roots6.538M including0.824M DSE. The13before-refactorweightupdates
are outsidefastcensus. Selectedpricing136/native135 is10157/4394I rowproduct,
319/190I basis/rrows and4423/1650I combinedfastDSE/basic. Native dualupdate841I;
itsRoseacceptedcounterpart isnotseparatelymeasured. Caller/clamp/return differences
remainexplicit; globalnative selfpools are notinclusivecostcounterparts.

Delegatedcompare.py dual-row-cost qualifies consumedinputs/rowoutput, basis/status,
dualoutputs, W/xb/flags againstcachedsource136→137 andnativePIN136→137. Rawpriority
comparison FAIL ispreserved: all444P+18slots differ, firstslot native1 versusRose0.
Disassembly provesP40=0selectedlistwalk neverreads thatpayload; nativeownpayload
staysunchanged. This explicitinactive/unpairedfield rule isnotnormalization or
fullstatePASS andmustnot beusedinvalue mode. Strictkernel−0/+0 frontier163 remains;
no newnativecapture crosses selected135. Native3root/3finaldetailcontrols agree.

Phase-at-entry counts remain132positivelevel/147level0listwalk/358valuestates;
rowproduct diagnostic3.205of3.857ms isvaluephase, unpaired withnativeafter163.
Framework_rhs isinactive throughoutDEGEN2pmode1. BroadROWresidual1.878ms andactual
coefficient-visitation versusconstruction costs arestillunresolved. No substantial
newremovablepass oracceptedtimegain isproved; whole targetgap stays2.74x.

The row-product operation census below closes this question withouta newoptimizer.
Workspace/touchedreset,wholepivotowner,prefilter/weightedtraversal,true_dpresorting,
rowmapcache andordinaryscatterboundselision trialsremainclosed. Evidence:
[degen2-row-split.json](evidence/degen2-row-split.json), rawsource/build/calibration/
parity/acceptedcounters in target/perf-degen2-row-split/; nativecentralreports in
 target/perf-degen2-native-row-cost/.

### Current DEGEN2 row-product operation diagnosis

The immutable accepted binaries execute 179 MUL and 179 ADD operations at Rose
pricing136/native pivot135. Central compare.py matches 14 ordered input rows and
179 coefficient/rho/product projections. Native reverses product operand roles
in seven cases; actual placements are retained. Both publish 105 row entries.

Strict ADD comparison stops at its first term (row295, column1): native really
adds its 1e-200 first-touch sentinel (bits16687e92154ef7ac), whereas Rose adds +0.
No normalization or later accumulator-state comparison. The published endpoint
proof is separate; this does not assert a full internal-state prefix. Native
uniform precomputes are zero, so the disassembly-only product-reuse lead is closed.
Selected complete instruction counts remain Rose10157/native4394; an equivalent
native scatter-only counter boundary is missing. No new cycle/time claim.

The whole-Rose 634-call instruction ledger reconciles to 12.610334M:

| Region | M instructions |
|---|---:|
| Coefficient arithmetic, first touch and sparse order growth | 4.011 |
| Dense support order and tolerance filtering | 3.410 |
| Slack arithmetic, indexing and order growth | 2.201 |
| Output gather, status, sign and multiplication | 1.927 |
| Length estimate, explicit allocation, cleanup and WORK | 1.060 |

There are 310 sparse-order and 324 dense-order calls, with 220946 coefficient
visits over 52238 rho slots. Explicit array/result allocation and local cleanup
cost 604149I (~.333% of solve instructions); growth stays inside the other loops.
Logical zeroing covers 5.580MB, not measured traffic or removable time. One product
per distinct row coefficient could avoid at most 156628 scalar MUL instructions
(~.0864% of solve instructions) before costs: no prototype follows. Seven complete
result/exit/signed-zero controls pass, root/detail counts agree, and accepted hashes
are unchanged. No optimization or parity/WORK/baseline change. See
[degen2-row-product.json](evidence/degen2-row-product.json) for full boundaries,
phase counts, raw first-difference evidence and centralized commands.

Mirror-update accounting below closes the diagnostic-only fp_rows question.
Rejected workspace, combined-cell, traversal, ownership and storage trials remain
closed; the whole-solve 2.74x runtime gap is unchanged.

### Current DEGEN2 mirror-update accounting

Immutable accepted55af/reference26a505, dual-release. Eleven new Rose result/exit/
signed-zero controls pass; no source/binary/WORK/dump/baseline change. Source global
update131 (pricing136/native pivot135) has row419, q140, update3→4, fill1599→1684,
i3 51→63 and i21 436→435. The next source update starts in the same factor epoch.

The complete selected boundary is Rose15337–15340 instructions versus native3182:
mirror update plus deferred WORK publication/vector retirement, pricing touch and
counter bookkeeping; auxiliaryFactor.update_column is separate. Native132a10
pricing touch contributes311I. Its support38 order, P60 argument14, m/r/q and all
444 P68 values before/after match centrally. Spike24 ordered bits, logical dense
spike444, eta indices/values, diagonal bits and captured factor projections also
match. Native nucleus numeric coefficients/physical row-pool layout remain
uncaptured; no full factor-state or arithmetic-visit equivalence claim. Native
observes four memcpy totaling60B and no exported allocator/free helpers inside
this boundary; retained pools have a separate construction/retirement lifetime.

Whole-Rose621-update calibrated instruction ledger:

| Region | M instructions |
|---|---:|
| Preparation before elimination | 0.206 |
| Row elimination, helpers and local cleanup | 2.470 |
| Spike/structural transitions/new diagonal and eta | 7.331 |
| Pool maintenance, including diagnostic metadata | 1.935 |
| Numeric/structural installation and local cleanup | 2.918 |
| WORK publication, retirement and pricing touch | 0.582 |

Total15.441646M I reconciles exactly with the complete roots and phase sums.
Source fp_rows has4055 tuples: initializer/push/cleanup223479I, already inside the
ledger, only.1233% of solve instructions. No fp_rows gating prototype follows.
Row elimination is2.469725M I (1.363% of solveI):188 empty-output calls26914I and
433 nonempty calls2.442811M I. All621 source row_cached tags match the pivot, but
eliminate runs before taking that marker; the marker suppresses WORK booking only.
Native131760 is not called at the selected update and cached eta data is consumed.
This is a dispatch difference, not yet a proof that Rose can reuse its BTRAN data.

Counters attach to the original inferior and exclude all debugger cycles/timers.
Root/detail selected counts differ by at most3I; coarse/fine whole roots differ141I.
No new codegen or wall-time claim. Raw ROW sentinel/P18 failures and triangular
output frontier163 remain preserved; new native captures stay within selected135.
Existing self-pool gaps are not additive with this ledger and cannot substitute
for a whole-run native inclusive update cost. Evidence:
[degen2-mirror-update.json](evidence/degen2-mirror-update.json), raw commands/counters/
snapshots in target/perf-degen2-mirror-update/ and centralized native reports in
 target/perf-degen2-native-update-cost/.

The original-stage comparison now passes159DEGEN2updatepairs through actual
sourcepricing163/nativepivot162. Boundcap162 produces this extra pricingcall;
source164/native163 remains firstknownreference signedzero divergence. AFIRO3/
SC50A3 and SC50A20/DEGEN2cap20 gates also pass, including nonemptyetas up to6.
Centralcompare.py dual-btran-eliminate pairs pricing/generation/selfpointer and
old updateepoch; updateincrements the counter beforeeliminate. CapturedUfactor
fields andorderednonpivotUvalues match independenteliminateeta values; pivotr
is preservedseparately. MXCSR/fctrl/fstat and scratch/seenrestoration match.
No normalization, wholefactor-state or fullarithmetic-history identity claim.

The prefix has162Ucaptures/159updates:43nonemptyetas/112entries,max9; three
unusedcaptures atpricing16/37/82 remain explicit lifecyclecosts. Logical sparse
valuecopy1792B is not traffic/time; completecapture/install/invalidation/retirement
costs are now measured by the actual-solver shadow study below. Selected136 givesU[419:+1,440:−1,441:−1,442:−1] versus
eta[440:−1,441:−1,442:−1]; Uwork12 is not stored dense_est27 (visited6).
Existing acceptedcost ledger shows covered159updates cost68746I, versus later
462updates2.400979M I, totaling2.469725M. The capturedprefix alone is2.78% of
eliminatecost andcannot justify a prefix-onlyoptimizer or whole-stage savings.
Sixcomplete result/exit/zerosigncontrols pass, no newreferencecapture.

An independent copy of completeoriginalkernels finds a real rounding restriction.
Only an appendeddiagnostic test changes the isolatedordinaryrelease testbuild;
removing it literally reproduces source57c86340. Three-node triangularU has
unitdiagonals, UTrow0[(1,1),(2,1)],row1[(2,0.1)],unitinput0. Actualoriginal
sparseBTRANU/eliminate kernels produce these eta-position2 bits:

| MXCSR rounding | BTRAN U | Eliminate |
|---|---|---|
| Nearest | bfeccccccccccccd | bfeccccccccccccd |
| Downward | bfeccccccccccccd | bfeccccccccccccc |
| Upward | bfeccccccccccccc | bfeccccccccccccd |
| Toward zero | bfeccccccccccccc | bfeccccccccccccc |

Eachcase startsfreshscratch/resetstatus andrestorescallerMXCSR. Exceptionstatus
matches but directed-rounding values differonebit. This test passes by asserting
the counterexample, not by provingreuse. BTRANpush andeliminate use different
sign/addition sequences; row_cached alone isinsufficient. Thebuild is separate
non-PGO release/glibc, functionalonly; accepted55af remainsunchanged. Allcaptured
LPpairs useddefaultMXCSR1fa0. No solveroptimization or numericrewrite.
See [degen2-btran-eliminate-reuse.json](evidence/degen2-btran-eliminate-reuse.json).
The complete frozen-kernel differential eligibility study now passes23,888
comparisons (f64/F80, both sparse-tail schedules, generated/extreme/nonfinite/
duplicate/threshold inputs and16 FP contexts). The diagnostic guard admits905
comparisons:581generated f64 and324captured-prefix schedule cases; all admitted
ordered values, visited/dense accounting, scratch restoration and FP-flag coverage
match. Quotient tracking exactly preserves original kernel output/WORK/flags in
all23,888 cases. F80 always falls back. Directed rounding, dropped pivots,
duplicate initial entries and non-unit input have concrete fallback witnesses.
These tests do not prove general reuse safety or provide a cheap runtime guard.

A read-only accepted-ELF source census qualifies all634 pivot-U captures and621
consuming updates under a stronger offline finite-factor/triangular/unique-row/
unit/retained-pivot/tail-only/FP-state predicate. Thirteen captures are discarded
at factor replacement. All4209 nonpivot entries include4088 consumed and121
unused entries;67344logical copy bytes are not traffic or timing. Complete
AFIRO3/SC50A20/DEGEN2cap20/640 result/exit/zero-sign controls match; diagnostic
plus plain-control wall times0.68/0.66/0.81/3.10s are excluded from cost claims.
The first script attempt failed syntax before solving and is retained excluded.
Whole-run eligibility remains source-only beyond the strict native frontier.

Evidence: [degen2-btran-eligibility.json](evidence/degen2-btran-eligibility.json).
The complete source-only U/eliminate comparison now passes all621updates:
433nonempty outputs/4088entries, with634captures and13unusedoutputs/121entries.
Delegated compare.py preserves original-order binary64 dense bookkeeping and
matches all orderedbits, factor epochs, FPstate and scratch; nativefrontier stays
162/sourcepricing163. Four complete capture/ordinary result/exit/zero-sign
controls pass, fullcapture12.97s; no debugger cost measurements.

An isolated complete-original-kernel release/TCmalloc component charges quotient
tracking, practical guards including source row_cached, capture/copy, faithful
visited/dense bookkeeping, installation and all cleanup, including unusedcaptures.
All four source-history replays and owner/epoch/row/FP negative guards pass.
Three fixedCPU2 rounds (64full histories/round) give instruction medians after
common originalU subtraction:

| Component scope | Instructions M |
|---|---:|
| Original eliminate,621calls | 3.077533 |
| TrackedU lowering difference | 0.150793 |
| Capture guards/copy/allocation/cleanup,634calls | 0.206598 |
| Consume/bookkeeping/install,621calls | 0.617185 |
| Complete lifecycle replacement component | 0.974576 |

The ledger sums; netcomponent instruction margin2.10293–2.10296M. This is separate
non-PGO codegen, not accepted costs: original component is24.61% above accepted
eliminate2.469725M. The tracking delta includes changed routine lowering, not just
finite-classification costs. Empty PMU/clock interval397I gives6.20I/history gross
probe overhead at64; no per-term instrumentation/multiplexing. Cycledifferences
havebothsigns andtimers varystrongly; noacceptedcycle/timegain. AFIRO3 andDEGEN2cap20
instruction margins are−186/−317I, SC50A20+4999.875I; no broad unconditional gain.

Evidence: [degen2-btran-lifecycle.json](evidence/degen2-btran-lifecycle.json).
The matching actual-solver shadow diagnostic is now complete. Its ownfixed39
PGO/systemTCMalloc build retains original eliminate and every original numeric
line/order. Full987 frozen parity passes inverify mode; all17201 original DEGEN2
dump records agree. All16 plain complete results/exits/signedzeros and72 counter
objects agree. Actual fullsolve has16constructors,634captures/4209entries,
621verifiedconsumes/4088entries,13retiredunusedcaptures andzero guardfailures;
no clones/restores occur onthis path. Separate release lifetime tests cover
clone/restore clearing, movedowner rejection, one-shotconsume, epoch/row/FP
failures, overwrite andDrop. Initial2x2 testfixture chose denseU; corrected4x4
fixture reaches sparseU, withinitialfailure source/binary/log retained. Measured
probe source/profile/costs didnot change.

Fixed-PC instruction scopes inthis diagnostic measure original eliminate621times:
2.492544M raw off and2.492545M on; after4markerI+2debuggerI/call,2.488818/2.488819M.
A known2-instruction calibration observes4I; debuggercycles/timers areexcluded.
Source/helper/localcleanup areincluded; outputeta retirement remainswithlater
update. Same-ELF helperinstruction deltas, includingallocation/free/telemetry:

| Actual shadow helper | Invocations | On−off instructions |
|---|---:|---:|
| BTRAN entry/invalidation | 1095 | 30259 |
| U capture/guards/copy | 1078 | 163482 |
| Guarded consume/faithful bookkeeping | 621 | 549435 |
| Shadow-result retirement after original solve | 621 | 30619 |
| Total | | 773795 |

Six balanced wholeCLI triplets give median recording−off774234.5I, range756400–
787011; allPMU enabled==running. The439.5I medianresidual isbelow~31k processI
variation, soUtracking/ctor/Drop/unusedcleanup andallocator effectscannotbe split
precisely. Offmedian194.582799M I versusaccepted194.409926M (+.172873M) isnotzero
instrumentation: loopclassification14190a..141949/141ddb..141e1c usesintegerchecks
andcmov evenwhenoff. Inputchecks140a03..140a25 brancharoundthem. Baselineskew also
containschangedoriginalcodegen/profile/layout; no causal exactsplit isproved.
Recordingmedian195.357524M I isdiagnostic, not an optimization. Cycle differences
havebothsigns; noacceptedtime/cyclegain or repeat/extension. No allocationobserver
or per-term counters; lifecycletelemetry cost remainsincluded.

Evidence: [degen2-btran-shadow.json](evidence/degen2-btran-shadow.json).
Accepted/source/baseline hashes unchanged, nativefrontier162/source163 unchanged;
no newnativecapture. Shadowcorrectness/costcovers testeddefaultmaskedrelease
FPstate. A futureoptimizer mustexplicitlyreject unsupportedunmaskedFPcontrols
andall priorunsafe cases, retaincompleteoriginalfallback andexactWORK/denseorder.
Onlyorder1.5–1.7M I (<1%solveI) remains asillustrativebudget, notremovabletime or
actualspeedup. The2.74xgap andother-stage/unowned remainder remain unchanged.
The one fixed same-live-pivot transfer trial is now REJECTED atitscountergate.
IndependentcompleteoriginalU/eliminate copies firstmatchunchangedruntime, then
candidatepasses24832releasecomparisons:1849admitted(581generatedf64+1268captured
634Utimes2schedules),11782F80fallbacks, generatedextreme/nonfinite/duplicate/
threshold/signedzero cases and16FPcontexts. All64exceptionmaskcombosreject
unsupportedcontrols; clone/restoreclearing,movedownedstate,one-shotconsume,
overwrite,epoch/row/FPguards pass. Frozenfallback andallWORKsites/orderunchanged.

Finalownfixed39PGO/systemTCMalloc ELFf8db0260 passesfast940 andall17201original
fullDEGEN2dumprecords (also1002cap20/4049cap149). All34counter completeobjects/
exits/signedzerosagree. Actualcompiledcaller branchcensus proves621reuse/0fallback
onfullsolve. Delegatedcentralcomparison stillpassesnative1..162/source163 and
preservesfirstnative163/source164basis22−0/+0; no normalization/newnativecapture.

The predeclaredCPU2 eight-paircountergate fails:

| Wholeprocess median | Accepted55af | Candidatef8db | Change |
|---|---:|---:|---:|
| Instructions M | 196.504071 | 194.880191 | −0.826385% |
| Cycles M | 145.405439 | 147.210458 | +1.241369% |

Everyinstructionpair improves,cyclepairs havebothsigns(4each); thisdescriptive
failedgate isnot a causalregressionproof. SmallAFIRO/SC50A/ADLITTLEinstruction
changes+.0383/−.0943/−.3134%,cycles+2.042/−6.396/−9.645%. Stopwithoutreserved
84117/84189, retry/extension,helper/flag/layout/PGO/WORKpolicyrescue. Fullcandidate
987/allrelease/full23panel notrunafterfailedgate; noproductionpromotion. Source,
completeoracles/tests,ownprofile/binary andallrawresults remainisolated. Accepted
source/binaries/baselines unchanged; no runtimegaincredited,2.74xgap remains.
Evidence: [degen2-row-transfer.json](evidence/degen2-row-transfer.json).

The unchanged accepted55af structural-update diagnostic is complete. Whole source
has510walks/7195steps (7006jumps,189shifts),127742first-length and124183repeated
predicate inputs,83272actualdependency-loop entries and2417shiftedpositions.
Low-overhead roots cost3.341818M I; finer disjoint parts sum3.342011M (193I/
.0058% residual retained). Beforestrictfrontier sourcepricing<=163 only65walks/
33749I; later445walks/3.308262M I are source-only, not pairednativeexcess.

| Accepted walk part | Instructions M | Solve instruction exposure |
|---|---:|---:|
| First A/B length filtering | .990797 | .547% |
| Repeated length for WORK | .862783 | .476% |
| Dependency slot-maximum scan | 1.200087 | .662% |
| Slot shifts and swaps, including jump WORK add | .157802 | .087% |
| Other setup/control/WORK formation | .130542 | .072% |

Nucleuspop/swap396calls cost7602I(.0042%solveI), excludingheld-decisionsearch.
SC50A20qualifies9walks4809I;AFIRO3hasnone. Sixcoarse/fourfinecomplete output/
exit/zero-signcontrols matchaccepted. WholefineGDB23.24s vscoarse2.34s isoverhead,
not timeattribution; allcycles/debuggertimers excluded, PMUno multiplexing.
Fiveoriginalinstructionsobserve7I; corrected2Iperstoppedinterval. Selectedcoarse/
finerootboth304I; earliersemantictrace326Ihas22Iextraobservercostandisexcluded.

Delegated actualsource136/native135 consumedprojection passes: drop row419,
index1,last117,size51→50; walkhole435/cur443/column428/length9>=distance8;
1shift/8orderedactualwrites,0dependency scans/jumps. Nativewidewalk131I includes
setup/inlineWORK/finalrestore, Rose304I excludesdeferredWORKpublication: not an
exactcompleteinclusivedelta. Pop21/15I bothomitcallerremovedpivotkindstore.
Raw14kinddifferences/13earlysourcepromotions remainexplicit; nofullintermediate
factoridentity ornewfloatinghistoryclaim. Nativecapture reproducescachedvalidated
nativeendpoints; strictnative162/source163frontierunchanged. No newsolverbuild.
Evidence: [degen2-structural-update.json](evidence/degen2-structural-update.json).

Dominantlengthfiltering1.854M I overlapsclosedwalklen/walkrefs; slotmaximum1.200M
walkmaxalso closed. Freshslot/pop mechanisms total<.1%solveI before requiredwork,
so no newbounds/copy/outlineprototype. These scopes arealreadyinside7.330873M
spike/transitionregion and15.441646Mupdateledger; ~3.981M transitionI remains
outsidewalk/pop, notunattributedsolvework oraremovablebound. Accepted2.74xgap
remains. PriorUtransfer/BTRANmemo/Lhandoff/fusion/rowmap/workspace/arena closed.
Next singleexperiment: costordinaryfactorB-listfrontinsert/frontfill-delete
moves/searches/growth/cleanup onacceptedELF, includeall16factor lifetimesbefore
anynewrepresentation. Benefitunknownuntil actualowner-inclusivecost isknown;
not rejectedB-referencecounts/spikemarkdedup orLUconstructorchainqueue trials.

### Current GREENBEA dual solve attribution

Matched bounded release-profile solve scopes count 5.4666B Rose instructions
versus 1.3461B reference (4.06x). Accepted Rose ELF is unchanged; a clock hook
asserts its two solve-boundary calls. The reference brackets GRBoptimize.
External read/setup/output is excluded; internal setup remains. API-call scope
is not proven identical to Runtime's internal timer boundary. FIFO boundary
work remains uncalibrated instrumentation cost. Every complete result/exit
matches its own original full-panel output after removing only the timer.

GREENBEA has 2,728 Rose and 2,688 reference samples at period 16M, zero loss/
throttle, and 99.806%/99.843% matching direct-event coverage. Per-process maps
and ELF PT_LOAD/symbol sizes recover samples missed by perf's disabled mmap
tracking; independent offset checks agree. Six Rose high-address samples stay
explicit unknown; no samples are dropped. Initial missing-mapping and ack
framing probes are retained and excluded. Small DEGEN2 qualification resolves
all 272 samples; its 90.23% reference coverage reflects short-process tails and is
not a quantitative phase bound.

Rose libc allocation has 536 samples (19.65%, approximately 1072M sampled
instructions/solve), versus no reference libc allocation samples. This does
not prove zero private/inlined reference allocation. At least 292 samples
(10.70%) belong to pivot-row work: 226 inlined mainline samples plus 66 separate
symbols. That accounts for 62.60% of iteration's 361 self samples, so the whole
iteration self share cannot be called generic control overhead. Arithmetic,
filtering/order and publication remain; skid prevents assigning individual
checks or stack reloads as causes. FTRAN WORK (joint_u) and output publication
(dense_spike) are separated from execution; broad support mass is 1138M versus
502.5M known reference FTRAN/L/U. These are self estimates with different
routine boundaries, not exact operation-equivalent excesses or removable time.

The allocation-owner lead is closed without implementation. The cached
frame-pointer ELF 68dc54d9 passes fast 940 and full 987 parity; nine ordinary/CLI
source files match current source exactly, while only older chol files differ.
Its own instruction skew is +0.789% on DEGEN2 and +0.586% on GREENBEA. This
qualifies owner discovery, not transferring its cost shares to accepted release.

The healthy GREENBEA FP profile has 2736 samples, zero loss/throttle and 99.515% direct
event coverage. All 542 libc allocator leaves reach a Rose-context caller;
14 remain Rose other rather than a specific buffer. Physical ancestry partitions
296 under refactor, 220 iteration/post and 26 driver-context leaves. The last group
contains 23 startup/end and 3 iteration/post return sites; it is not all lifecycle
work. Leading returns are e72b2 (61 leaves): outgoing RefStruct retirement at
mirror replacement; e76e3 (54): local Structure retirement; 10b0d1 (32):
recompute_d returning from true_d, before dual_residual_check. Disassembly
confirms callees because line attribution can name the following statement.
F80-looking aliases do not establish extended execution. Those 147 samples are
5.37% of FP self, not a removable-time bound.

All leading owners overlap settled cross-epoch banks/generation arenas,
Structure publication/retirement and reduced-cost rebuild/presorted-column
studies. Remaining construction/pricing, row/framework and triangular-storage
owners supply no distinct lifetime. Existing dense-row, workspace, row-cell,
outlining, FTRAN owner and whole pivot handoff trials remain settled. No new
buffer/cache/arena trial follows. Evidence:
[greenbea-allocation-owners.json](evidence/greenbea-allocation-owners.json),
all FP/provenance/parity/counters/frames/source-site proofs in
`target/perf-greenbea-owners/`. Release solve scopes remain separately recorded
in [greenbea-phase-profile.json](evidence/greenbea-phase-profile.json).

### Current host: existing PGO/system-TCMalloc configuration

The existing fixed39 build passes full frozen parity **987/987**, all **292**
runnable release tests (22 ignored), and complete cross-build comparison of
92 Rose plus 92 reference results/exits. Only timer fields are removed;
canonical JSON comparison preserves signed zero. Source is unchanged. Its CLI
SHA256 is `55af661754b952817fe9438dbeb157413f1063d105c27413a4c502acf541f6fb`.
Pinned Rust 1.98.1 and matching LLVM 22.1.8 profile tools use the existing
39 bounded training cases, 64-byte loop alignment and system-tcmalloc feature.
Profile SHA256 is `a018c80fcf31c7a35a87d9f7be872b8f6174a92969af2c801ed4dc8a2581e46b`.
Ubuntu gperftools 2.18.1-1 is hash-checked and extracted into the task sysroot;
no system package was installed. Loader bindings confirm genuine TCMalloc.

The fixed 23-row current-host comparison is complete: **2/23** descriptive
median wins, both PDHG. Same CPU 2, one warmup plus three measured fresh-process
runs per solver/row, rotating order, identical parameters and timer definitions.
GREENBEA PDHG remains iteration-limited in both solvers. These measurements are
not two independent 32AB+32AA optimization acceptance screens. The ordinary
accepted CLI is unchanged; this opt-in configuration is not promoted as a new
runtime fix. The build/training/test workflow took 737 seconds in isolation.

| Row | Rose ms | Gurobi ms | Rose / Gurobi |
|---|---:|---:|---:|
| degen2.primal-c | 53.430 | 22.085 | 2.419 |
| degen2.dual-r | 48.166 | 16.650 | 2.893 |
| degen3.primal-c | 836.768 | 323.327 | 2.588 |
| degen3.primal-r | 881.971 | 338.821 | 2.603 |
| perold.primal-c | 783.889 | 289.241 | 2.710 |
| pilot4.primal-r | 74.208 | 25.504 | 2.910 |
| 80bau3b.primal-c | 716.392 | 177.107 | 4.045 |
| primal_5959.primal-r | 83.765 | 69.112 | 1.212 |
| degen3.dual-c | 401.439 | 148.556 | 2.702 |
| degen3.dual-r | 398.572 | 144.588 | 2.757 |
| perold.dual-r | 201.585 | 65.795 | 3.064 |
| 80bau3b.dual-r | 434.474 | 135.702 | 3.202 |
| greenbea.dual-c | 1084.460 | 376.091 | 2.884 |
| greenbea.dual-r | 1067.992 | 355.073 | 3.008 |
| mode_switch.dual-r | 1047.013 | 160.296 | 6.532 |
| primal_5959.dual-r | 111.958 | 98.496 | 1.137 |
| quad_2717.dual-r | 110.624 | 98.811 | 1.120 |
| fx_quad_obj_5371.dual-r | 38.710 | 36.871 | 1.050 |
| degen3.barrier-r | 244.120 | 118.492 | 2.060 |
| 80bau3b.barrier-r | 308.504 | 183.202 | 1.684 |
| greenbea.barrier-r | 1341.209 | 763.236 | 1.757 |
| 80bau3b.pdhg | 3141.776 | 3381.691 | 0.929 |
| greenbea.pdhg | 1596.285 | 1667.711 | 0.957 |

Evidence: [optimized-build-qualification.json](evidence/optimized-build-qualification.json).
All source/dependency/profile hashes, logs, loader proof and raw outputs are in
`target/perf-optimized-build/`. Run the saved `rose-pgo-tcmalloc` with
`LD_LIBRARY_PATH=$PWD/target/perf-optimized-build/sysroot/usr/lib/x86_64-linux-gnu`.
The isolated build command and environment are in `build-command.json` there.
The reference remains glibc and uses the same reconstructed library. No
numerical architecture or individual optimization gain follows from this sweep.

Complete `solve305` ownership is rejected without implementation. The remaining
owner work is dense-mask metadata, temporary reuse and output publication;
original triangular/matrix arithmetic remains. Cached per-record diagnostic
residual scopes are 1.562% of whole 80BAU and 1.352% GREENBEA, including required
arithmetic and instrumentation. They are not removable bounds or current PGO
measurements. No composed triangular/route/packing/buffer prototype follows.
[solve305-owner-review.json](evidence/solve305-owner-review.json).

The same qualified optimized ELF now has solve-only mode-switch attribution.
After healthy DEGEN2 qualification, eight mode-switch instruction recordings,
eight cycle recordings and eight direct counter solves preserve all complete
outputs/exits. Period32M yields 856 instruction and 811 cycle samples, zero
loss/throttle; all instruction symbols resolve, while seven high cycle IPs
remain explicit unmapped. Represented instructions are 99.992% of an independent
counter batch; represented cycles are 90.545%, so cycle shares are coarse and
sensitive to cross-batch variability. No timing gain follows from these probes.

TCMalloc and libc memory respectively occupy 20.56% and 18.57% of instruction
self; these totals do not reveal removable lifetimes. The newly isolated
`empty_tail` maximum fold has 41/856 instruction but only 12/811 cycle samples
(4.79%/1.48%). A maintained slot bitmap is distinct from the earlier dependency
walk maximum and rejected reverse-slot scan, but must transport membership
through every slot rotation/swap, transient demotion and final replacement.
That maintenance has no demonstrated substantial whole-lifetime benefit; no
prototype/census follows. Do not infer a new arena, workspace or cache trial
from the memory shares. Evidence:
[optimized-mode-scope.json](evidence/optimized-mode-scope.json), raw profiles,
per-process maps, full validation and excluded throttled small probe in
`target/perf-optimized-scope/`. Numerical source and both qualified CLIs remain
unchanged; next evidence must identify a distinct substantial execution cost.

Earlier small dual probes and PEROLD ancestry remain separately labeled in
their evidence. The unthrottled PEROLD frame-pointer batch (24 full solves,
cap 1300, all optimal 1273/27 refactors) has 3736 samples, 99.57% direct coverage,
and refactor 27.03% inclusive; FTRAN/full-needed self families do not by themselves
justify another ownership trial. Old throttled exploratory recordings remain
excluded. The current GREENBEA solve-scope profile is the next attribution source.

### Current mode-switch and refactor attribution

The accepted release CLI now has a healthy full mode-switch profile (cap 12000,
optimal after 11,393 iterations, 17 refactorizations): 4,896 instruction samples
at period 8M, zero loss/throttle events, and 99.932% direct-event coverage.
Release self samples put memset at 7.41%, fresh at 4.92%, nucleus at 4.62%,
needed FTRAN at 3.96%, update at 3.27%, and ordinary FTRAN at 2.49%. Substantial
native allocator/copy work is separate. These whole-batch instruction fractions
include startup, serialization, and bash/timeout launchers. JSON validation is
outside PMU measurement. These are not solve-time or removable bounds.

The preceding frame-pointer ELF has unchanged ordinary source and matches all
mode-switch outputs. Its healthy period-32M recording has 2,464 samples, zero
loss/throttle events, and 99.502% direct-event coverage. Refactorization is
45.37% inclusive. Disjoint nearest physical stages include nucleus at 8.85%,
fresh at 8.12%, refactor remainder at 8.12%, pricing reset at 6.01%, structure
remainder at 5.19%, and auxiliary LU at 3.69%. These totals contain necessary
numerical work and allocation. The original period-8M stack recording has six
throttle/unthrottle pairs and is retained but excluded quantitatively.

Cached healthy PEROLD stacks reproduce 1,010/3,736 refactor samples (27.03% of
the batch): nucleus 7.47%, fresh 6.58%, auxiliary LU 3.80%, refactor remainder
3.75%, and pricing reset 1.85%. Several mode-switch memory samples resolve only
to ancestor calls into iteration, needed FTRAN, or framework RHS; those return
addresses do not identify particular buffers. Identical-code folding can give
ordinary drop/generic code F80 names. Source-line attribution and fixed-period
sampling are approximate. Do not combine different builds' fractions.

All 194 current profiling/counter output objects and exits agree with the
accepted sweep. The profiling ELF has fast parity 940/940; release full parity
987/987 and all 292 runnable release tests pass. Numerical source and accepted
CLI are unchanged. [mode-switch-scope.json](evidence/mode-switch-scope.json), all
recordings, excluded samples, source/offset maps, and raw outcomes are in
`target/perf-mode-scope/`; PEROLD reproduction is in `target/perf-refactor-stages/`.

### Markowitz summary scope

The distinct incremental clean/maximum certificate proposal is closed at complete
cost scope, before an eligibility census or implementation. Existing linked-chain
and borrowed-search trials do not test unchanged numerical maximum reuse. A
timer-only release-library probe retains every numerical/search/WORK statement;
full source reconstructs literally after removing test scopes. Six small models
cap300 precede PEROLD cap1300 and mode-switch cap12000. All18 complete paired
results agree, including explicit bits for every public floating-point value.

On mode-switch, complete search including helper descendants occupies3.76–4.13%
of instrumented solve time. All cleanup plus ordered maximum folds occupies
3.004–3.591%; maximum alone is0.439–0.458%. These include required first/changed
work, not merely potential unchanged repeats. PEROLD search is2.669–2.749% and
cleanup plus maximum is1.320–1.394%. Candidate costs depend on opposite-axis
counts even when the own chain/maximum is unchanged, so caching a winner cannot
follow from these certificates. A maximum reuse path also needs an independent
FP-status/control proof. Nested timers/inlining affect these ordinary ThinLTO/
glibc diagnostic scopes; they are not exact accepted-PGO removable bounds.
No sufficient substantial residual is established. No shadow certificate/key
census, cache/borrowed-search rescue, prototype or formal screen follows.
[Summary scope](evidence/markowitz-summary-scope.json), complete source/probe,
release-test ELF, counts and outputs: `target/perf-markowitz-scope/`.

### Current compact extended execution

The remaining numerical extended BTRAN is cold after compact state replay.
A CPU2 ordinary-release library census first checks all six small models
(cap300), then complete fx/Quad rows (caps5500/17000). All20 complete results
match plain solves and the numerical source reconstructs literally after
removing test scopes. Across four runs per extended row, all uncached F80
BTRAN occupies0.407–0.447% of fx and0.116–0.142% of Quad measured solve time;
dense-U intersections are0.190–0.207%/0.060–0.067%. Fresh-only intersections
are smaller still. These diagnostic scopes include original helpers/work/output
and are not exact optimized-build removable bounds. Provider caches are created
per Solver; this does not measure a retained cross-solve cache.

The distinct whole-kernel native transaction is closed without a prototype.
Canonical admission, specified binary80 operation boundaries, cumulative status,
staged state and software fallback could replace per-primitive validation, but
current actual call scope cannot justify it. Historical BTRAN shares/attempt
counts below do not transfer to this source. Conditional subnormal/state proof
and its limitations are retained, without claiming a tested native executor.
[Transaction scope](evidence/btran-transaction-scope.json), all probes/logs and
release-test ELF in `target/perf-btran-transaction/`.

The accepted optimized Quad CLI now has717 solve-only cycle samples, zero
loss/throttle,716 resolved samples and96.89% independent-counter coverage.
All64 profiling/counter outputs and exits match the complete optimized panel.
The numerical F80 BTRAN range has zero self samples. Current self attribution
is solve_primal23.57%, QuadRhs::update_row4.60%, rank_products_observed3.07%,
QuadMain::apply2.09%, and memcmp2.09%; these include required work, not removable
owners. Allocation cannot be assigned to query keys: hits already borrow inputs,
use stack keys or retained scratch and allocate retained keys only on misses.
Returned output clones are separate. Persistent component storage also has no
new qualified residual beyond the settled pivot/basis/arena trials.

An isolated source-matched PGO/system-TCMalloc profiling build with line tables
and frame pointers passes full frozen parity987/987. Eight DEGEN2 solves precede
32 Quad profile and32 counter solves; all72 complete objects/exits match the
optimized panel. Quad has765 cycle samples, zero loss/throttle and764 resolved
leaves. Direct instructions are1.228% above the accepted optimized build; sampled
events/direct cycles are103.06% across the separate batches. Keep this ELF's
fractions separate from the accepted build. Of174 solve_primal self samples,
30 have outermost source at the Devex update and25 at breakpoint ratio; inline
chains and unknown locations are retained. Physical external callers remain
separate from leaf/source attribution, with library unwinding limitations.
No removable whole-owner bound or runtime optimization follows from these
fractions. Read-only architecture review closes the controller lead without a
prototype: ratio/point/classification/pricing fusion crosses required numerical
and state boundaries; storage replacements overlap settled ownership trials.
In particular41 external leaves reach only the return from solve_primal and
cannot identify particular buffers or retirement costs.
[Controller profile](evidence/quad-current-controller-profile.json),
commands, separate target, binary and logs are in `target/perf-quad-profile/`;
normalized stacks/source maps are in `target/perf-optimized-scope/`.
Both accepted CLIs and numerical source are unchanged.

### Factor geometry recurrence

A bounded ordinary construction census covers AFIRO, SC50A, BLEND, ADLITTLE,
BORE3D and DEGEN2 under both controlled simplex methods, cap300. All12 complete
library results match their plain same-build solves. None of82 constructors
repeats its ordered input geometry within a solve. The literal numerical
constructor source reconstructs after removing only test observation hooks.

An offline order-preserving component census finds3331 components and3180
strict geometry repetitions, but zero repeated sparse Cartesian weight. All63
positive-weight components are unique; even an optimistic best-prior matching
prefix covers zero of4717 proxy units. Repetitions include3117 one-by-one and63
nontrivial zero-sparse components. Singleton/pair work is unweighted; there are
no dense cases. These are necessary geometry comparisons, not execution traces
or cost bounds. They do not test arbitrary row isomorphisms/internal subtraces,
full/large histories or numerical independence from global search and WORK.
No large capture or executable program follows this small-history result.

[Factor geometry evidence](evidence/factor-geometry-recurrence.json) records keys,
source restoration, independent review and limitations. All raw records,
probe/analyzers and release-test ELF remain in `target/perf-factor-geometry/`.
Production source and both qualified CLIs are unchanged.

### Barrier factor executor

The whole-register x87 Kahan transaction is timing-inconclusive and unapplied.
It retains all five rounded intermediates, software zero identities, complete
numeric state restoration and original software fallback. This differs from
the historical guarded-product trial and the forbidden interval shortcut.
An independent complete current kernel precedes edits. Generated64896 terms
and855947 live terms across224 actual dots match every intermediate sign,
mantissa and exponent; all seven library histories match complete results/float
bits. The corrected FP harness passes896 environment cases with actual native
execution asserted on admitted inputs, including all eight caller registers.
Its initial standalone restore lacked XMM clobber declarations; the separate
corrected proof adds them without a production timing variant or retry.

Admission/state-inclusive component ratios on80BAU are0.306–0.357. Matching39-case
PGO/system-TCMalloc passes fast940. The fixed4mixed/2AA cost pilot preserves
all28 complete results/exits and gives80BAU solve-9.288%, instructions-6.945%,
cycles-8.162%. The first fixed32AB/32AA screen58131 gives median-5.213%, interval
[-7.359%,-0.0473%],23/32 improving, sign p0.0100308; both orders improve. AA has
p0.05010 and interval[-2.848%,+7.562%]. Sign significance and conservative noise
separation fail. Stop without58197, extension, zero/schedule/admission rescue,
full candidate gates or a new benchmark win claim. Root source and both accepted
CLIs remain unchanged. [Kahan transaction](evidence/kahan-register-transaction.json),
all original/candidate sources, complete proofs, corrected tests, binaries,
profiles, schedules and raw outputs: `target/perf-kahan-transaction/`.

The separate generic assembly program lead is closed without a prototype.
The earlier compiled assembly experiment covers only the standard branch;
its all-front test calls that standard implementation even on generic-dispatch
shapes. New complete entry scopes split the two actual production branches.
All18 instrumented library results match their disabled plain solves across
the six small models first, then80BAU3B/GREENBEA (cap1000). The numerical source
reconstructs literally after removing test-only scopes. Entire generic assembly
occupies2.20–2.28%/3.88–3.99% of the two large instrumented solves, including
required arithmetic/output work. All generic calls are narrow fronts; no
oversized-front dispatch occurs. This does not support a substantial removable
owner. Ordinary ThinLTO/glibc library scopes include instrumentation effects and
are not exact accepted-PGO removable bounds. No executor, matching PGO or formal
timing screen follows. [Generic assembly scope](evidence/generic-assembly-scope.json),
original/probe/source reconstruction, release-test ELF, all shapes and runs:
`target/perf-generic-assembly-scope/`.

The solver compiles child-update geometry once and retains its scratch across
Newton factors. Logical GEMM destinations update the parent directly without
materializing gathered or zero-product buffers when admitted by the exact
range/context checks. Discarded upper outputs still execute. Four independent
double RHS lanes share each live rectangular-triangular coefficient. All original
per-output reductions, precision paths, fallback inputs and WORK calls remain.

Two fixed CPU 2 DEGEN3 screens (32 AB + 32 AA each, seeds 45631/45697, cap 100)
pass: median solve reductions 13.286% and 14.379%, confidence intervals
[-15.253%,-9.897%] and [-16.585%,-12.314%], 31/32 and 32/32 improving. Both
orders improve and healthy AA uncertainty is clear of both intervals.
DEGEN3 instructions fall 14.485%; DEGEN2 barrier instructions fall 9.246%.
These counter reductions describe the complete executor, not additive gains
from its components. The materialized compiled-plan/packet control reduces
DEGEN3 instructions 14.638%, so no independent direct-C timing gain is claimed.

The complete 23-row sweep preserves all 184 complete result/exit objects. DEGEN3
barrier solve median changes 0.328336 to 0.271563 s (-17.291%); 80BAU3B and
GREENBEA barrier changes remain within their measured spreads. A separated
primal_5959.primal-r short-sweep slowdown receives a fixed reversed-role screen:
median accepted/candidate -2.127%, interval [-7.474%,+1.816%], sign p 0.0551;
AA is healthy (p 0.2153), with interval [-5.966%,+0.506%]. No regression is
established; this does not establish timing equivalence. Other median increases
remain within the sweep spreads. Full parity 987/987 and all 292 runnable
release tests pass (22 ignored).

Configuration is ordinary ThinLTO/glibc, with no global target-CPU flags or PGO.
This acceptance sweep precedes the reconstructed-reference comparison above;
it supplies no independent Gurobi win count or speed ratio. Source hashes, all 23 before/after rows,
controls, source variants and exact qualification are in
[barrier-factor-executor.json](evidence/barrier-factor-executor.json); raw data
and proof are in `target/perf-factor-engine/`. Prior isolated schedule/packet
screens remain inconclusive and were not repeated.

## Historical AMD EPYC 9V74 measurements

A recurring compact primal basis now restores its accepted variable order into the live buffer and installs its two factor handles directly. Fresh, singular and diagnostic factorizations retain the original work and pivot rules.

The full primal_5959.primal-r sweep took 0.094669 → 0.078778 seconds versus the preceding commit `e4df849b` (-16.8%). The prespecified randomized timing comparisons measured -11.52% and -10.51%. Frozen parity passed 987/987, 283 runnable release tests passed, and all full-length mode-switch WORK records matched. Full outputs and exit codes agreed on every benchmark. 4 reversed-role controls covered flagged regressions; none established a slowdown.

The table gives same-CPU medians after one warm-up and three fresh-process timed solves for each of the four executables, with rotating order. All rose binaries use genuine TCMalloc. The patched Gurobi library hash is `26a505171314b3300b754b4eff5f5719366b47ccf6f3d59770df99e3f256ee30`. Current rose wins **5/23**.

| Run | Gurobi (s) | Uploaded rose (s) | Preceding rose (s) | Current rose (s) | Current / Gurobi |
|---|---:|---:|---:|---:|---:|
| degen2.primal-c | 0.019700 | 0.070753 | 0.051559 | 0.086430 | 4.39 |
| degen2.dual-r | 0.018759 | 0.055531 | 0.043750 | 0.061261 | 3.27 |
| degen3.primal-c | 0.285311 | 0.863128 | 0.810041 | 0.782957 | 2.74 |
| degen3.primal-r | 0.324561 | 0.911660 | 1.049446 | 0.787360 | 2.43 |
| perold.primal-c | 0.271103 | 0.759709 | 1.040564 | 0.823650 | 3.04 |
| pilot4.primal-r | 0.042927 | 0.087877 | 0.094110 | 0.093569 | 2.18 |
| 80bau3b.primal-c | 0.175074 | 0.989292 | 0.795513 | 0.806673 | 4.61 |
| primal_5959.primal-r | 0.090778 | 0.207468 | 0.094669 | 0.078778 | 0.87 |
| degen3.dual-c | 0.146547 | 0.433089 | 0.374213 | 0.381341 | 2.60 |
| degen3.dual-r | 0.137805 | 0.445857 | 0.433228 | 0.487388 | 3.54 |
| perold.dual-r | 0.067197 | 0.211987 | 0.265044 | 0.225320 | 3.35 |
| 80bau3b.dual-r | 0.127239 | 0.603060 | 0.557311 | 0.475941 | 3.74 |
| greenbea.dual-c | 0.354868 | 1.187191 | 1.196558 | 1.176313 | 3.31 |
| greenbea.dual-r | 0.311835 | 1.201368 | 1.261721 | 1.301069 | 4.17 |
| mode_switch.dual-r | 0.237340 | 2.336878 | 1.144451 | 1.191617 | 5.02 |
| primal_5959.dual-r | 0.129710 | 0.269287 | 0.117423 | 0.096123 | 0.74 |
| quad_2717.dual-r | 0.116162 | 0.331867 | 0.124672 | 0.123069 | 1.06 |
| fx_quad_obj_5371.dual-r | 0.045304 | 0.110665 | 0.041111 | 0.040029 | 0.88 |
| degen3.barrier-r | 0.135277 | 0.276130 | 0.284263 | 0.297667 | 2.20 |
| 80bau3b.barrier-r | 0.168752 | 0.361854 | 0.320291 | 0.363318 | 2.15 |
| greenbea.barrier-r | 0.699639 | 1.461157 | 1.342396 | 1.373002 | 1.96 |
| 80bau3b.pdhg | 3.719941 | 6.932070 | 3.035841 | 2.795014 | 0.75 |
| greenbea.pdhg | 1.524285 | 2.748687 | 1.449413 | 1.465858 | 0.96 |

### Compact primal basis installation

A complete immutable installation retains the variable indices in the accepted pivot order. A recurring basis copies that order into its existing storage. A fresh installation moves the order into place. Direct mutable factors live behind owned pointers, so passing a shared installation between the builder and solver transfers a small handle. All pivot choices, factor lifetimes, cache matches and eviction rules are unchanged. The standalone CLI reads its constant factor diagnostics inside the solve timer. Library callers retain live environment checks.

The two disjoint installation routines saved a median 8.91% of the full parent primal_5959.primal-r solver time across eight paired instrumented runs. That figure includes routine entry and return, handling cache misses and reinstalling the factor. All four original compact solves and 2,696,537 chronological WORK records agreed across the two modes. An independent check changed diagnostic settings within an actual library solve and compared factor counts and diagnostic streams.

The full ordinary primal cases added 457,580 matching WORK records on 6 rows.

The first whole-solve comparison measured -11.52% with interval [-12.99%, -10.85%]. Rose won 60/64 pairs. The confirm whole-solve comparison measured -10.51% with interval [-11.65%, -9.57%]. Rose won 107/128 pairs. 

Fresh CPU6 comparisons against the supplied patched Gurobi used 48 mixed and 48 identical-Gurobi pairs per row.

| Row | Rose vs Gurobi | 95% interval | Faster pairs | Paired test |
|---|---:|---:|---:|---|
| primal_5959.dual-r | -6.65% | [-12.62%, +0.60%] | 31/48 | Inconclusive |
| fx_quad_obj_5371.dual-r | -16.11% | [-18.39%, -6.79%] | 38/48 | Passed |
| quad_2717.dual-r | -7.29% | [-11.09%, -1.44%] | 33/48 | Inconclusive |
| primal_5959.primal-r | -6.76% | [-10.53%, -0.46%] | 32/48 | Inconclusive |

An independent CPU3 comparison of the same executable used 48 mixed and 48 identical-Gurobi pairs on each original PDHG solve.

- greenbea.pdhg measured -6.25%, interval [-8.81%, -1.79%], 37/48 faster pairs. Inconclusive.
- 80bau3b.pdhg measured -19.30%, interval [-22.23%, -17.23%], 48/48 faster pairs. Passed.

Rose won 14/69 individual timed comparisons against Gurobi in the 23-row sweep.

### Preceding release e4df849b

The compact controller can reuse the original general cost-change BTRAN when its ordered inputs, pricing state, rank owner and factor cursor agree. It advances to the exact original successor, applies the original partial writes and issues the same chronological WORK calls.

The full quad_2717.dual-r sweep took 0.167139 → 0.164750 seconds versus the preceding commit `2d278c18` (-1.4%). The prespecified randomized timing comparisons measured -13.60% and -10.87%. Frozen parity passed 987/987, 283 runnable release tests passed, and all full-length mode-switch WORK records matched. Full outputs and exit codes agreed on every benchmark. 8 reversed-role controls covered flagged regressions; none established a slowdown.

The table gives same-CPU medians after one warm-up and three fresh-process timed solves for each of the four executables, with rotating order. All rose binaries use genuine TCMalloc. The patched Gurobi library hash is `26a505171314b3300b754b4eff5f5719366b47ccf6f3d59770df99e3f256ee30`. Current rose wins **5/23**.

### Phase-one compact cost updates

A repeated general cost BTRAN can restore the original sparse writes and pricing list and advance to its exact factor successor. Each original WORK addition still happens in order. On the independent pre-integration numerical branch redoing 211 real hits verified the full factor, numerical outputs, list and chronological WORK operands. The phase occupied 25.65 ms before the change and 5.97 ms afterward on full quad_2717. The final price and phase-one releases both exercise the independent full-model Devex and higher-level pricing checks.

The first full quad_2717.dual-r comparison measured -13.60%, 95% interval [-14.41%, -12.63%], with 169/192 pairs faster. The second full quad_2717.dual-r comparison measured -10.87%, 95% interval [-13.20%, -6.42%], with 135/192 pairs faster. 

In the 23-row sweep, Rose won 13/69 timed comparisons. It finished faster in all three rounds on 80bau3b.pdhg, greenbea.pdhg.

A separate CPU2 comparison used 48 mixed and 48 identical-Gurobi pairs on each compact run. Every original output and exit code was checked.

| Run | Rose vs Gurobi | 95% interval | Faster pairs | Outcome |
|---|---:|---:|---:|---|
| primal_5959.dual-r | -0.53% | [-5.77%, +7.98%] | 24/48 | Inconclusive |
| fx_quad_obj_5371.dual-r | -5.16% | [-9.09%, +5.84%] | 28/48 | Inconclusive |
| quad_2717.dual-r | +0.75% | [-4.47%, +7.22%] | 23/48 | Inconclusive |
| primal_5959.primal-r | +4.69% | [+0.12%, +11.91%] | 15/48 | Gurobi faster |

A separate CPU8 follow-up checked every additional ordinary or barrier row that had beaten Gurobi in the sweep.

- degen2.primal-c measured +130.40% against patched Gurobi, with 95% interval [+121.44%, +144.98%] and 0/48 pairs faster. The prespecified test did not pass.

An independent CPU6 direct check of this committed solver used 48 mixed and 48 identical-Gurobi pairs per PDHG row. On 80bau3b.pdhg Rose was 19.00% faster, interval [-21.63%, -17.01%], and won 46/48 pairs. On greenbea.pdhg it was 6.31% faster, interval [-8.83%, +1.30%], and won 29/48. Only 80BAU3B passed the prespecified direct comparison. All 388 complete original outputs and exit statuses were reverified. Evidence in next/committed-pdhg-reference.

### Preceding release 2d278c18

The compact primal solver retains the ordered pricing row independently of mutable prices. A recurring transition restores only the costs, weights and extended twins on its original support. A miss runs the original row update.

The full fx_quad_obj_5371.dual-r sweep took 0.066168 → 0.046607 seconds versus the preceding commit `caf46745` (-29.6%). The prespecified randomized timing comparisons measured -20.79% and -20.77%. Frozen parity passed 987/987, 282 runnable release tests passed, and all full-length mode-switch WORK records matched. Full outputs and exit codes agreed on every benchmark. 1 reversed-role controls covered flagged regressions; none established a slowdown.

The table gives same-CPU medians after one warm-up and three fresh-process timed solves for each of the four executables, with rotating order. All rose binaries use genuine TCMalloc. The patched Gurobi library hash is `26a505171314b3300b754b4eff5f5719366b47ccf6f3d59770df99e3f256ee30`. Current rose wins **3/23**.

### Integrated compact pricing rows

A compact pivot row retains its checked geometry and up to eight sparse numerical price transitions. A hit copies the four price values only at participating positions. The earlier diagnostic reran 45,888 original operations and checked all four arrays, including untouched positions. The existing Devex test exercises the original lower cache directly and still requires more than 1,000 actual lower-kernel cache hits. Separately, the integrated full-model test keeps pricing reuse enabled, checks the entire numerical state and ordered row against the original kernel on every verified hit and requires more than 1,000 genuine hits per model. The oracle excludes its own work from the live cache-hit counter. A new independent PGO build and new timings were required after the earlier prototype failed the original lower-cache coverage assertion.

The first full fx_quad_obj_5371.dual-r comparison measured -20.79%, 95% interval [-21.96%, -19.34%], with 172/192 pairs faster. The second full fx_quad_obj_5371.dual-r comparison measured -20.77%, 95% interval [-23.43%, -18.44%], with 164/192 pairs faster. 

In the 23-row sweep, Rose won 10/69 timed comparisons. It finished faster in all three rounds on 80bau3b.pdhg, greenbea.pdhg.

A separate CPU2 comparison used 48 mixed and 48 identical-Gurobi pairs on each compact run. Every original output and exit code was checked.

| Run | Rose vs Gurobi | 95% interval | Faster pairs | Outcome |
|---|---:|---:|---:|---|
| primal_5959.dual-r | -3.61% | [-11.88%, +3.67%] | 27/48 | Inconclusive |
| fx_quad_obj_5371.dual-r | +3.26% | [-6.31%, +12.70%] | 21/48 | Inconclusive |
| quad_2717.dual-r | +17.82% | [+13.26%, +21.23%] | 5/48 | Gurobi faster |
| primal_5959.primal-r | +4.47% | [-9.66%, +22.88%] | 19/48 | Inconclusive |

### Preceding release caf46745

A fresh compact main factor can own the numerical effects of reconstructing its primal and dual values. Matching basis, ordered inputs and floating-point state allow the solver to restore the original writes and advance to the exact existing transpose successor. It still evaluates the drift and pricing decisions and books each original WORK site.

The full primal_5959.dual-r sweep took 0.138375 → 0.122820 seconds versus the preceding commit `f4fca9d5` (-11.2%). The prespecified randomized timing comparisons measured -11.93% and -9.16%. Frozen parity passed 987/987, 280 runnable release tests passed, and all full-length mode-switch WORK records matched. Full outputs and exit codes agreed on every benchmark. 8 reversed-role controls covered flagged regressions; none established a slowdown.

The table gives same-CPU medians after one warm-up and three fresh-process timed solves for each of the four executables, with rotating order. All rose binaries use genuine TCMalloc. The patched Gurobi library hash is `26a505171314b3300b754b4eff5f5719366b47ccf6f3d59770df99e3f256ee30`. Current rose wins **1/23**.

### Compact primal and dual reconstruction

On the independent compact LU branch the original primal and dual reconstruction routines consumed a combined 29% of full 5959 dual time. Restoring their verified numerical writes and advancing the exact BTRAN successor saved 10.60% of solve time in the instrumented comparison, including lookup and restore costs. The two independent whole-solve comparisons on that source measured 10.71% and 10.77% gains. The production comparisons below use the actual reduced-cost parent.

The first whole 5959 dual comparison changed time by -11.93%, 95% interval [-14.51%, -9.37%]. Rose improved in 145/192 pairs. The confirm whole 5959 dual comparison changed time by -9.16%, 95% interval [-11.77%, -6.22%]. Rose improved in 143/192 pairs. 

The CPU7 columns use one warm-up and three full runs per role. The last three columns use 48 mixed Rose/Gurobi pairs and 48 identical-Gurobi pairs on CPU2, with all original results and exit codes checked.

| Run | CPU7 Rose (s) | CPU7 Gurobi (s) | CPU2 paired change | 95% interval | Test |
|---|---:|---:|---:|---:|---|
| primal_5959.dual-r | 0.200750 | 0.247125 | +13.45% | [+9.72%, +19.98%] | Gurobi faster |
| fx_quad_obj_5371.dual-r | 0.056303 | 0.040817 | +20.47% | [+17.08%, +28.04%] | Gurobi faster |
| quad_2717.dual-r | 0.157666 | 0.140523 | +28.31% | [+22.46%, +35.95%] | Gurobi faster |
| primal_5959.primal-r | 0.211758 | 0.160609 | +13.93% | [+7.78%, +19.93%] | Gurobi faster |

### Preceding release f4fca9d5

At a compact LU state, an exact repeat of the extended dual vector, active columns and cost-free mode can reuse the earlier native reduced costs. The solver restores only the slots written by the original pass and executes the final active column again to preserve x87 effects. Any incompatible floating-point or basis state runs the original pass.

The full primal_5959.dual-r sweep took 0.146243 → 0.108262 seconds versus the preceding commit `2731b202` (-26.0%). The prespecified randomized timing comparisons measured -8.77% and -8.93%. Frozen parity passed 987/987, 279 runnable release tests passed, and all full-length mode-switch WORK records matched. Full outputs and exit codes agreed on every benchmark. 3 reversed-role controls covered flagged regressions; none established a slowdown.

The table gives same-CPU medians after one warm-up and three fresh-process timed solves for each of the four executables, with rotating order. All rose binaries use genuine TCMalloc. The patched Gurobi library hash is `26a505171314b3300b754b4eff5f5719366b47ccf6f3d59770df99e3f256ee30`. Current rose wins **2/23**.

### Native reduced costs on the preceding solver

More than 99% of the original slow-path calls on the four compact models used an already observed extended dual vector and active column set. In a complete 5959 dual run on the independent LU branch, the native reduced-cost pass fell from 13.03 ms to 4.40 ms. That comparison includes checks on both misses and hits, result copies and re-execution of the final column.

The first full 5959 dual comparison on the direct parent measured -8.77% with interval [-11.10%, -6.69%] and 145/192 pairs faster. The confirm full 5959 dual comparison on the direct parent measured -8.93% with interval [-10.16%, -8.08%] and 159/192 pairs faster. 

Independent CPU5 times for the four affected original rows used one warm-up and three measured solves per executable.

| Row | Parent (s) | Current (s) | Patched Gurobi (s) |
|---|---:|---:|---:|
| primal_5959.primal-r | 0.115778 | 0.105213 | 0.082683 |
| primal_5959.dual-r | 0.159208 | 0.138820 | 0.114746 |
| quad_2717.dual-r | 0.142661 | 0.110191 | 0.099937 |
| fx_quad_obj_5371.dual-r | 0.052750 | 0.062587 | 0.035242 |

### Preceding release 2731b202

Each compact LU installation owns a bounded sequence of complete main-factor states. FTRAN, BTRAN, factor updates and counter changes advance to an existing successor when the entire ordered input and numerical environment agree. A different input resumes with a private mutable factor. Each caller receives the same ordered output and books the same WORK increments.

The full quad_2717.dual-r sweep took 0.184991 → 0.141771 seconds versus the preceding commit `ba40e67f` (-23.4%). The prespecified randomized timing comparisons measured -31.66% and -31.82%. Frozen parity passed 987/987, 278 runnable release tests passed, and all full-length mode-switch WORK records matched. Full outputs and exit codes agreed on every benchmark. 5 reversed-role controls covered flagged regressions; none established a slowdown.

The table gives same-CPU medians after one warm-up and three fresh-process timed solves for each of the four executables, with rotating order. All rose binaries use genuine TCMalloc. The patched Gurobi library hash is `26a505171314b3300b754b4eff5f5719366b47ccf6f3d59770df99e3f256ee30`. Current rose wins **1/23**.

### Compact simplex factor reuse on the preceding solver

The solver observed 259,079 recurring factor transitions in the four original compact solves. In the diagnostic run, redoing 2,269 of them checked the complete successor factor, return values, ordered WORK and floating-point flags. The formatted production source passed a separate check of 1,437 calls for floating-point environment side effects. The preceding LU release and its direct parent agreed on 11 alternative floating-point control settings. All four complete solves reproduced the original 1,319,417 records with replay enabled and 1,377,120 records with the extended diagnostic enabled.

| Original full solve | CPU8 paired change | 95% interval | Faster pairs |
|---|---:|---:|---:|
| quad_2717.dual-r, first | -31.66% | [-33.25%, -30.56%] | 191/192 |
| quad_2717.dual-r, repeat | -31.82% | [-33.16%, -30.56%] | 188/192 |
| fx_quad_obj_5371.dual-r | -31.88% | [-35.86%, -23.89%] | 30/32 |
| primal_5959.dual-r | -29.53% | [-32.40%, -25.24%] | 32/32 |

On CPU5 the independent panel used a warm-up and three fresh solves per binary for each affected row. All 48 results and exit codes matched the corresponding original.

| Run | Parent (s) | Current (s) | Patched Gurobi (s) |
|---|---:|---:|---:|
| primal_5959.primal-r | 0.121862 | 0.122861 | 0.080842 |
| primal_5959.dual-r | 0.220272 | 0.158302 | 0.120692 |
| quad_2717.dual-r | 0.253557 | 0.162094 | 0.106417 |
| fx_quad_obj_5371.dual-r | 0.112011 | 0.070354 | 0.049463 |

The separate CPU5 reference comparisons used 48 mixed pairs and 48 identical-Gurobi pairs per run.

- greenbea.pdhg measured -5.48% against patched Gurobi, 95% interval [-14.39%, +1.99%], with 30/48 pairs faster. The prespecified test did not pass.
- 80bau3b.pdhg measured -21.78% against patched Gurobi, 95% interval [-25.53%, -15.47%], with 48/48 pairs faster. The prespecified test passed.

### Preceding release ba40e67f

The barrier solver still computes every rounded extended product, error and sum of its original Kahan recurrence. For same-sign operands at the same exponent, it can carry their exact signed significand difference directly into the following compensation subtraction, normalizing or falling back at the original precision boundaries. Exceptional operands retain the original operations. The two binary64 conversions and final homogeneous expression stay in their original order.

The full 80bau3b.barrier-r sweep took 0.266418 → 0.241987 seconds versus the preceding commit `db13d89c` (-9.2%). The prespecified randomized timing comparisons measured -5.13% and -5.62%. Frozen parity passed 987/987, 277 runnable release tests passed, and all full-length mode-switch WORK records matched. Full outputs and exit codes agreed on every benchmark. 5 reversed-role controls covered flagged regressions; none established a slowdown.

The table gives same-CPU medians after one warm-up and three fresh-process timed solves for each of the four executables, with rotating order. All rose binaries use genuine TCMalloc. The patched Gurobi library hash is `26a505171314b3300b754b4eff5f5719366b47ccf6f3d59770df99e3f256ee30`. Current rose wins **2/23**.

### Preceding release db13d89c

The compact primal solver retains each ordered basis and its unscaled factor charges. When it revisits that basis, it restores the main numerical factor into its owned buffers and reuses a compatible auxiliary factor after clearing its eta updates. Each installation applies its current work scale in the original order.

The full primal_5959.dual-r sweep took 0.257346 → 0.162680 seconds versus the preceding commit `c337c24d` (-36.8%). The prespecified randomized timing comparisons measured -12.36% and -11.26%. Frozen parity passed 987/987, 275 runnable release tests passed, and all full-length mode-switch WORK records matched. Full outputs and exit codes agreed on every benchmark. 2 reversed-role controls covered flagged regressions; none established a slowdown.

The table gives same-CPU medians after one warm-up and three fresh-process timed solves for each of the four executables, with rotating order. All rose binaries use genuine TCMalloc. The patched Gurobi library hash is `26a505171314b3300b754b4eff5f5719366b47ccf6f3d59770df99e3f256ee30`. Current rose wins **2/23**.

### Preceding release c337c24d

The forward operator compiles contiguous input blocks in the retained reflected-primal coordinates. It keeps the original scalar prefix, eight independent arithmetical lanes, horizontal reduction tree and tail, then feeds each small group of complete row products directly to the dual projection. On x86-64, gather dispatch bounds both dimensions before using signed vector offsets.

The full greenbea.pdhg sweep took 1.975066 → 1.384803 seconds versus the preceding commit `52a90742` (-29.9%). The prespecified randomized timing comparisons measured -18.63% and -20.44%. Frozen parity passed 987/987, 271 runnable release tests passed, and all full-length mode-switch WORK records matched. Full outputs and exit codes agreed on every benchmark. 5 reversed-role controls covered flagged regressions; none established a slowdown.

The table gives same-CPU medians after one warm-up and three fresh-process timed solves for each of the four executables, with rotating order. All rose binaries use genuine TCMalloc. The patched Gurobi library hash is `26a505171314b3300b754b4eff5f5719366b47ccf6f3d59770df99e3f256ee30`. Current rose wins **3/23**.

### Independent barrier comparisons on ba40e67f

On CPU4 one warm-up and three complete independent solves per role measured 0.242183 seconds for the ba40e67f 36-step 80BAU3B barrier solver, 0.273153 seconds for its direct parent, and 0.117189 seconds for the patched reference. The original 38, 105 and 16 barrier records and the numerical outputs on 80BAU3B, GREENBEA and DEGEN3 all matched. The parent and candidate also matched on 185,880 public result bits and floating-point environment records, including exceptional operands, with identical unmasked traps.

The original CPU2 four-way sweep put GREENBEA barrier 9.45% above its parent. A separate 96-pair reversed-role trial on the complete 103-step solve put the parent at +1.10% against the new solver, 95% interval [-1.45%, +3.30%]. Its identical-binary control sign p was 1.0000. The comparison did not establish a regression.

### Direct PDHG measurements on the preceding db13d89c solver

The preceding compact LU executable ran 48 mixed rose/Gurobi pairs and 48 Gurobi/Gurobi controls per PDHG row on CPU5. Each sample used a fresh process, the original benchmark settings and each solver’s own timer. Every rose output and exit matched the original rose. Every Gurobi result matched the original patched reference. The evidence contains the 388 raw outputs, run orders and numerical audits.

| Run | Iterations | Rose vs Gurobi, paired median | 95% interval | Rose faster | Identical-Gurobi interval | Prespecified criteria |
|---|---:|---:|---:|---:|---:|---|
| 80bau3b.pdhg | 42,154 | -22.30% | [-25.78%, -17.93%] | 48/48 | [-5.85%, +0.16%] | Passed |
| greenbea.pdhg | 20,000 | -9.29% | [-11.28%, -3.39%] | 37/48 | [-3.29%, +3.60%] | Inconclusive |

The preceding c82803ab executable recorded −23.48% and −9.81% on these two runs in its CPU8 randomized comparison. Its CPU2 80BAU3B barrier win did not repeat on CPU3. Its raw observations remain in reference-confirmation.tar.xz.

### Preceding release 52a90742

Short columns with the same number of terms share an instruction sequence for four independent transpose dots. Their results feed the box projection directly. Ordinary steps retain the primal vector, anchor and reflected vector in this order, while the forward matrix uses their remapped indices. Observation and restart boundaries restore the original coordinate order.

The full 80bau3b.pdhg sweep took 4.750526 → 3.324115 seconds versus the preceding commit `2969bd40` (-30.0%). The prespecified randomized timing comparisons measured -28.27% and -26.50%. Frozen parity passed 987/987, 270 runnable release tests passed, and all full-length mode-switch WORK records matched. Full outputs and exit codes agreed on every benchmark. 9 reversed-role controls covered flagged regressions; none established a slowdown.

The table gives same-CPU medians after one warm-up and three fresh-process timed solves for each of the four executables, with rotating order. All rose binaries use genuine TCMalloc. The patched Gurobi library hash is `26a505171314b3300b754b4eff5f5719366b47ccf6f3d59770df99e3f256ee30`. Current rose wins **0/23**.

### Preceding release 2969bd40

The primal solver retains the structural nonzero count and uses the changed variables to assemble its ordered candidate history. The full 5,022-pivot 80BAU3B trajectory matched all 72,071 WORK, 3,949 ORD and 249 CORD records and every output byte after timing normalization.

The full 80bau3b.primal-c sweep took 1.355297 → 1.757710 seconds versus the preceding commit `8eb01cef` (+29.7%). The prespecified randomized timing comparisons measured -15.87% and -14.49%. Frozen parity passed 987/987, 268 runnable release tests passed, and all full-length mode-switch WORK records matched. Full outputs and exit codes agreed on every benchmark. 9 reversed-role controls covered flagged regressions; none established a slowdown.

The three-round CPU2 sweep put 80BAU3B controlled primal 29.69% above its predecessor. The two independent CPU8 comparisons each found a gain. In the prespecified follow-up on CPU6, 48 mixed drift-balanced quartets and 48 identical-binary quartets put candidate solve time 13.33% lower, 95% interval [−17.61%, −8.24%], and child CPU time 12.28% lower. On primal_5959 dual, the corresponding +2.64% solve-time and +1.43% process-CPU differences were within the identical-binary intervals. The table retains the original warmup-plus-three measurements in full, including these two cases.

The table gives same-CPU medians after one warm-up and three fresh-process timed solves for each of the four executables, with rotating order. All rose binaries use genuine TCMalloc. The patched Gurobi library hash is `26a505171314b3300b754b4eff5f5719366b47ccf6f3d59770df99e3f256ee30`. Current rose wins **1/23**.

### Preceding release 8eb01cef

During ordinary PDHG steps, dense primal projection and Halpern averaging use matching slices, while the dual pass consumes the resulting matrix product in independent coordinates. This lets the compiler vectorize the dense loops. Checks, restarts, residual calculation and termination use the original materialized states and the original sparse products.

The full 80bau3b.pdhg sweep took 4.794536 → 3.246851 seconds versus the preceding commit `f23a2eb0` (-32.3%). The prespecified randomized timing comparisons measured -31.66% and -31.88%. Frozen parity passed 987/987, 266 runnable release tests passed, and all full-length mode-switch WORK records matched. Full outputs and exit codes agreed on every benchmark. 5 reversed-role controls covered flagged regressions; none established a slowdown.

The table gives same-CPU medians after one warm-up and three fresh-process timed solves for each of the four executables, with rotating order. All rose binaries use genuine TCMalloc. The patched Gurobi library hash is `26a505171314b3300b754b4eff5f5719366b47ccf6f3d59770df99e3f256ee30`. Current rose wins **0/23**.

### Preceding release f23a2eb0

BTRAN reuses its sparse numeric workspace through the upper transpose, reverse row updates and lower transpose. On eligible ordinary factors the upper solve visits occupied slots in slot order, and the reverse eta index visits updates reached from active pivot rows in descending order. Dense and partial paths preserve their clearing rules. The unchanged full BTRAN kernel is the test oracle. A separate mode-switch instrumentation run measured 6.856 GB fewer large zero-allocation requests than the preceding solver.

The full mode_switch.dual-r sweep took 1.097172 → 0.965360 seconds versus the preceding commit `8da37a69` (-12.0%). The prespecified randomized timing comparisons measured -19.06% and -19.66%. Frozen parity passed 987/987, 264 runnable release tests passed, and all full-length mode-switch WORK records matched. Full outputs and exit codes agreed on every benchmark. 3 reversed-role controls covered flagged regressions; none established a slowdown.

The table gives same-CPU medians after one warm-up and three fresh-process timed solves for each of the four executables, with rotating order. All rose binaries use genuine TCMalloc. The patched Gurobi library hash is `26a505171314b3300b754b4eff5f5719366b47ccf6f3d59770df99e3f256ee30`. Current rose wins **0/23**.

### Preceding factor-update change

The factor update now retains its numerical workspace, snapshots just the referenced slots and maintains the set of empty tail columns that qualify for demotion. Selection preserves the reference slot order. On mode-switch, an allocation trace of the complete solve counted 45,520 fewer large zero-allocation requests (3.073 GB) and 11,380 fewer large copies (1.366 GB). The replacement issued 15 additional large fills totalling 1.8 MB.

Against `6a5f7ade`, full 11,393-pivot mode-switch took 1.398917 → 1.119417 s, a -20.0% change. Two independently seeded 64-pair comparisons measured -13.34% and -12.29%, with 59 and 50 improving pairs and healthy controls. All 202,772 WORK records agreed with the uploaded solver. The fixed39 PGO executable passed 987/987 frozen runs and 260 runnable release tests. All 23 full output and exit comparisons agreed. The sweep flagged 4 additional reversed-role checks. None established a slowdown under the predeclared criteria.

All four executables ran in fresh processes on logical CPU2, with one warm-up, three timed solves and rotating order. All rose builds used genuine TCMalloc. The Gurobi reference is the supplied patched library with SHA-256 `26a505171314b3300b754b4eff5f5719366b47ccf6f3d59770df99e3f256ee30`. The new solver wins **0/23** against that reference.

### Preceding forward-L change

The preceding forward-L change retains its sparse forward-L workspace and level queues across calls. An index of eta input rows schedules the updates reached by a right-hand side in their original order. An exact bounded suffix sum books the original sparse WORK estimate, including skipped terms. Nonfinite coefficients, exact tolerance ties with etas and dense handoffs retain the general path. The preceding solver `60298604` already caches the level bound and reuses the forward U-column workspace.

On the full 11,393-pivot mode-switch solve, the new stage took 1.388262 → 1.191710 seconds versus `60298604` (-14.2%). Two randomized 48-pair CPU2 screens measured -11.38% and -12.89%, each with 45/48 improving and healthy same-binary controls. All 202,772 WORK records agreed byte for byte. The captured-factor and incremental-eta differential tests, 987/987 frozen parity and all 260 runnable release tests passed. Eleven reversed-role controls covered rows flagged by either sweep or by the raw-file audit. None established a slowdown. On GREENBEA dual the identical-binary control was biased (p=.021), and the mixed pairs split 8–8. A separate fixed, drift-balanced 48-quartet CPU5 comparison plus 48 same-binary quartets gave -1.95% solve time and -0.83% whole-process CPU time. Both excluded the flagged increase, with healthy identical-binary controls. See the retained raw quartets and intervals.

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
3.59% and cycles3.17%. Two independent32AB+32AA screens (CPU4, cap 800,
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
DEGEN2 work. PILOT4 primal release at cap 800 has367 calls under one model owner,
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
The predeclared PILOT4 primal-release screen (CPU4, cap 800, seed32217,
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
etas and nnz_eta change through updates. These historical repeat counts predate
the accepted AuxiliaryOwner and do not justify another constructor cache.
Current-source exact-key measurements after that owner find zero repeats in
twelve small primal/dual solves (AFIRO through DEGEN2, cap 800). Full bounded
fx/Quad have only 57/140 remaining constructor calls, 25/18 unique keys and 31/122
eight-entry FIFO hits. Their hit constructors total 0.158/0.362 ms in a single
release-test diagnostic; this excludes lookup/storage/ownership overhead and is
not a solver speedup or accepted-release cost bound. All 14 complete SolveResult
objects match uninstrumented solves. No shared-base cache is implemented from
these stale historical counts. Evidence: [auxiliary-current-reuse.json](evidence/auxiliary-current-reuse.json).
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

## Remaining performance gaps

The numerical source is `18061fe9`, including the accepted barrier factor
executor. Current-host fixed23 medians show ordinary ThinLTO/glibc **0/23**
wins and separately qualified fixed39 PGO/system-TCMalloc **2/23**, both PDHG.
The latter's largest gaps are mode-switch dual 6.532x, controlled 80BAU primal
4.045x, PEROLD dual 3.064x and GREENBEA dual-release 3.008x. These three-sample
medians are descriptive, not optimization acceptance screens. Current exact
configuration tables and raw evidence are at the top; older host/allocator
results are historical. Full frozen parity and all runnable release tests pass
for both unchanged qualified binaries. No current large structural candidate
has qualified; closed experiments below must not be retried under new names.

## Extended BTRAN feasibility

Current compact execution is measured above: actual remaining numerical BTRAN
is cold. The following baseline/native experiments are historical feasibility
and rejected-trial records, not current cost shares or priority evidence.

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

The fixed PILOT4 screen35217 (32AB+32AA, CPU4, cap 800) is inconclusive: median
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


The complete pre-U forward-column publication/expansion/replacement chain is
2.72–3.64% on SC50A through DEGEN2 in bounded release-test solves, including
required filtering and clock overhead. Tiny AFIRO measures 7.25%; its short
call intervals are particularly sensitive to instrumentation overhead.
DEGEN2 cap300 has 405 full FTRAN, 300 needed FTRAN
and 294 updates; actual producer tuples total 390,992 logical bytes. The clock
control's empty scoped interval is about 39 ns (94 ns total scope overhead),
material at this scale. No owner rewrite or PEROLD escalation is justified.
Post-U publication is 7.24% on DEGEN2, with output ordering alone 4.71%; natural
traversal reuse there was already timing-inconclusive below and is not retried.
All 36 final pinned enabled/disabled solve results agree. These diagnostic
nested test scopes are not accepted-release removable bounds. Numerical source
and accepted CLI are restored byte for byte.
[forward-owner-scope.json](evidence/forward-owner-scope.json), raw source,
calibration and both exploratory/final measurements in `target/perf-forward-owner/`.

Accepted-source complete normal/controller scopes cover all six small models,
then default-limit 80BAU3B (optimal/36) and GREENBEA (suboptimal/103). Numeric
factorization is about 44.82%/50.77%, with child updates 19.95%/26.37% and
assembly 9.96%/10.41%; normal solves are 15.50%/24.84%, including complete
forward/backward 7.39%/13.55%. Refinement is 2.34%/9.38%, including correction
solves 0%/5.58%. 80BAU3B has no HSD raw or correction solves; GREENBEA has
72 raw and 131 correction solves. Trial selection occupies 8.76%/6.58%; the
complete packet executor tested against this scope is inconclusive below.
These are nested release-test wall scopes, not removable accepted-release
bounds. All 48 final paired complete solve objects agree. Earlier GREENBEA
cap 150 observations are excluded: a nondefault limit bypasses later termination
branches and reaches iteration_limit/150 rather than the benchmark's exit.
[barrier-normal-scope.json](evidence/barrier-normal-scope.json), raw probes,
source, outputs and earlier observations in `target/perf-barrier-normal/`.

The complete inner GEMM action census does not justify a routing rewrite.
Cached AFIRO through DEGEN3 geometry agrees with all prior call counts. DEGEN2
has 12,834 ordered actions/factor; DEGEN3 has 109,213, requiring 0.87 MB even
with optimistic eight-byte descriptors. The arithmetic-free release proxy,
including one construction per unique shape, saves only about 0.086 ms over
DEGEN2's factor history and 1.512 ms over DEGEN3's. This proxy is not a solver
speedup or removable bound: arithmetic, inlining and descriptor loads matter.
Healthy accepted-release DEGEN2 self samples place scatter admission at
1.67% instructions/1.64% cycles and S4x4 at 3.76%/3.92%; whole native scatter
is 5.85%/7.67%. These include launcher/parsing/startup. Both recordings have
zero lost/throttle events but represent only about 95.5% of matching direct
batch counts; short-process sampling tails limit attribution.

The optimistic four-row packing census distinguishes per-GEMM and per-child
lifetimes: DEGEN2 needs 5,491 versus 2,287 four-value records/factor; DEGEN3
102,680 versus 23,175. BLEND/ADLITTLE have no cross-batch reuse and BORE3D has
714 versus 510 records. Geometry reuse alone does not establish packing cost
or eligibility; the current small-model S4x4 scope includes required arithmetic
and B loads. Repeated A/B scans exceed producer R/D elements by 4.18x on
DEGEN2 and 6.48x on DEGEN3, but admission's measured scope is small. No
producer certificate, persistent packing or descriptor executor is implemented.
Keep original temp arithmetic, full logical validation, flags and fallbacks.
All final profiling/direct-counter batches agree on complete results/exits;
full parity987/987 and release tests pass. Source and accepted CLI are unchanged.
Parser-only and throttled exploratory
recordings are retained and excluded. [gemm-routing-scope.json](evidence/gemm-routing-scope.json),
source/proxy/census and all raw evidence in `target/perf-gemm-routing/`.

Complete predictor/HSD pairing has no overlapping correction rounds on
GREENBEA. Its 72 pairs occur at iterations 31–102. Every predictor takes zero
corrections; HSD takes 97 attempted corrections across the pairs. Of the pairs,
52 refine and 20 do not; all are double, with early CG off. Other calls account
for the remaining 34 of 131 corrections. Pairable initial solves occupy 7.65%
of the instrumented whole solve, or 8.70% with initial residuals. Whole paired
transactions occupy 14.90% but contain unmatched HSD tails. No paired executor
follows from this scope. All six small models have no HSD calls. All 21
enabled/disabled diagnostic pairs agree on complete SolveResult objects.
These nested release-test scopes are not removable accepted-release bounds.
Source and CLI are restored byte for byte. [hsd-pair-scope.json](evidence/hsd-pair-scope.json),
per-round counts, scopes, source, and outputs are in `target/perf-hsd-pair/`.

An entering-only immutable-L prefix census gives no cache implementation lead.
AFIRO, SC50A, BLEND, ADLITTLE, BORE3D, and DEGEN2 primal release (cap 300) have
14/32/64/104/101/300 entering FTRAN calls. None repeats the same entering column
and ordered input support within a LUOUT epoch. Omitting values and entering
identity gives one possible support hit on SC50A, two on BLEND, and zero
elsewhere. Across epochs, the same column/support repeats 0/2/8/22/15/15 times;
these do not establish reuse within one immutable L lifetime. All 12 dumped/plain
complete output objects and exits agree. LUOUT epochs approximate owners; other
L calls and cross-epoch owner sharing are unmeasured. Counts are not cost shares
or exact cache hit rates. The inclusive L profile also contains mutable etas.
No cache or larger-model escalation follows. [l-prefix-reuse-scope.json](evidence/l-prefix-reuse-scope.json),
raw evidence and scripts are in `target/perf-l-prefix/`.

Direct output-lifetime scopes rule out the consuming factor producer/installer.
Formation plus all Structure retirement is 3.109% of the full mode-switch solve,
including required singleton arithmetic. Fresh-prefix metadata outside the
previously shared L/transpose and Markowitz conversion scopes is about 0.08%.
SC50A through DEGEN2 formation/retirement scopes are 3.2–4.0%; tiny AFIRO is
5.85% with more material instrumentation overhead. All 21 enabled/disabled
pairs agree on complete SolveResult objects. Formation/retirement counts balance;
mode-switch has 17 formations, 15 successful-refactor retirements, two startup/
other retirements, and no retry/error retirement. Test-only zero-sized first/last
fields bracket existing field destruction without moving retirement earlier.
These are nested release-test wall scopes, not accepted-release removable bounds;
marker overhead is included. No consuming installer follows. Source and accepted
CLI are restored byte for byte. [factor-output-scope.json](evidence/factor-output-scope.json),
source, calibration, marker limits, and raw data are in `target/perf-factor-output/`.

The current accepted-release controlled-primal profiles cover DEGEN2 first, then
80BAU3B. Healthy 80BAU period-32M self samples put feasibility classification at
13.76%, score at 11.87%, ordinary row update at 7.21%, BTRAN at 6.21%, dot_lanes
at 4.82%, and candidate-list update at 3.78%. DEGEN2 has 2,434 healthy samples
and 97.64% direct-event coverage; 80BAU has 2,384 samples and 99.70% coverage.
The original 80BAU period-8M recording throttled and is retained but excluded.
All 178 profiling/counter complete outputs and exits match the accepted sweep.
These are whole-batch instruction fractions, not solve-time or removable bounds.
[primal-controlled-scope.json](evidence/primal-controlled-scope.json), all raw
recordings, outcomes, counters, and scorer annotation are in `target/perf-primal-scope/`.
The prior feasibility-tracking and global pricing-bucket trials remain settled.

Query-time exact whole-score reuse was tested independently of global buckets.
On 80BAU, 2,647,720 of 3,410,347 division-producing queries repeat complete keys
and satisfy conservative FP-status admission (77.64%). The isolated prototype
retains original list traversal, scalar fallback, and WORK. Cache-protocol tests
cover all 16 rounding/flush modes and special/random values; every live small/
80BAU query is checked against the frozen scalar, with real cache hits asserted.
Fast parity passes 940/940. The fixed cost pilot (seed 53101, four mixed and two
AA pairs per model) rejects it: 80BAU instructions +8.677%, cycles +5.821%, solve
median +5.861%, both orders slower. DEGEN2 instructions rise 2.982%.
No formal 53131/53197 screen, extension, full candidate gates, or sweep follows.
The accepted source/CLI are unchanged. The initial FP oracle's optimized-away
call and its opaque-call correction are retained with the other evidence.
[score-query-reuse.json](evidence/score-query-reuse.json), complete count/prototype/
frozen-oracle/pilot artifacts are in `target/perf-score-{reuse,candidate}/`.

The distinct complete batched pricing executor is also rejected at the cheap
small-model component gate. Four independent AVX lanes compute eligible scores;
original ordered full/list/changed/old-list transitions and WORK remain. Frozen
complete score/price/best/scanned/update methods precede edits and reconstruct
literally. All12small-model histories across AFIRO through DEGEN2 match each
original complete pass's result, list/level/failure state, WORKw bits and final
MXCSR, plus the complete plain SolveResult. No exceptional/generated full oracle
or full987 candidate qualification follows because no useful cost is shown.
The two DEGEN2 complete-pass replay ratios are1.172/1.267; CPU was unpinned in
this initial cheap discriminator, so this is not a controlled regression claim.
No80BAU prototype replay, matchingPGO build, formal screen, packet-width or mask
rescue follows. Runtime source and qualified binaries remain unchanged.

Current optimized solve-only80BAU cycle attribution has569samples, zero loss/
throttle,568resolved symbols and98.785%of the independent direct count. All24
DEGEN2/80BAU profiling/counter complete outputs match the optimized full panel.
Scoring is inlined into price/update_list:23+41samples(11.25%) include necessary
list/bucket work; the old separate11.87%instruction share cannot be transferred.
Evidence: [primal-score-batch.json](evidence/primal-score-batch.json), independent
originals/candidate/release-test ELF and all component data in
`target/perf-score-batch/`; current raw maps/offset annotations/counters are in
`target/perf-optimized-scope/`. Require a distinct substantial execution mechanism
before another primal implementation; neither arithmetic-only packed division
nor old scalar shares establish its economics.

Ordinary-primal rho/tau BTRAN traversal pairing is rejected at scope. Ordered
coefficient-load overlap is 40.65% on SC50A, 34.13% on BLEND and 26.97% on
DEGEN2 release, but only 4.281% on controlled 80BAU. The latter has 3,225 pairs
and 724 rho-only transactions; its pairs perform 4,252,758 coefficient visits,
with 182,068 shared under monotone U-tail, reverse-eta and descending LT merges.
Even the offline weighted row-LCS upper bound shares only 462,348 (10.872%).
LT has 1,836 level/level pairs, 1,330 level/dense and 59 dense/dense; common
physical support does not supply a cheap order-preserving bucket schedule.
Divisions, predicates, sparse scheduling, WORK, output/clearing and one-lane
arithmetic remain. These are active-row counts, not removable time bounds.
All four enabled/disabled complete SolveResult objects agree; root source and
accepted CLI are unchanged. No paired executor, path-specific packed trial,
timing screen or larger-model escalation follows. Dynamic-rank rho products
and rho-eligible column-order Devex dots remain separate. Evidence:
[paired-btran-scope.json](evidence/paired-btran-scope.json); original row sequences,
probe, analysis and accepted gates are in `target/perf-btran-pair/`.

Complete independent-output dense LDL batching is rejected at its fixed cost
pilot. Two-row geometry covers 95.56% of DEGEN2 and 96.66% of DEGEN3 dot terms.
The isolated whole-stage SSE2 executor retains scalar diagonal pivots, each
row's original eight chains/fold/tail, separate reciprocal scaling, and the next
temp boundary. All-masked/finite admission, full matrix snapshot, per-result
checks and matrix/MXCSR restoration before scalar replay are included in costs.
The frozen complete stage passes 2,944 generated cases across all 16 rounding/
flush modes, 1,272 captured small fronts, and 2,130 actual DEGEN3 fronts with
full matrix/d/rcp/FP-status equality. Real admission/pair/replay counts are
asserted. Fast parity passes 940/940; the isolated baseline is byte-identical
to the accepted CLI. The fixed small counter discriminator admits DEGEN3 after
DEGEN2 instructions fall 9.690%, cycles 0.483%; its solve timing is noisy.
DEGEN3 cost pilot 54109 (four mixed/two AA pairs) measures instructions -5.756%,
cycles -2.513%, solve -2.268%, below the declared 5% minimum; AA median +9.070%.
All 28 pilot complete results/exits match accepted outputs. Stop without
54131/54197, extension, width/snapshot tuning, full candidate gates or sweep.
No runtime change accepted. Root source/CLI remain unchanged; full parity987
and all292 runnable release tests pass. Whole solve305 routing/codegen remains
unqualified beyond the settled triangular plan. Evidence:
[ldl-output-batch-scope.json](evidence/ldl-output-batch-scope.json); complete
prototype/oracle/input/tests/assembly/pilot artifacts in `target/perf-ldl-output/`.

The exact Plan census counts 1,008 wide GEMM calls/factor on DEGEN3, 799 zero-C;
fifteen factors describe 83.2 MB of logical cn buffers, 37.4% unused cells whose
arithmetic is retained. Original R-packing copies cover 13.0 MB; contiguous-C
copies 5.46 MB. These are logical volumes, not time bounds. Direct destinations
now cover the admitted double-precision boundary; packed-R ownership and
factor-generation recycling remain unimplemented, and prior arena-only trials
provide no qualifying gain. [factor-materialization.json](evidence/factor-materialization.json).

Current non-level sparse-L replay census rules out that duplication on the
measured targets. DEGEN2, PILOT4 and PEROLD have no Kahn calls and no dense
fallback replay. PEROLD's 2,687 L entries split into 1,237 initially dense and
1,450 level-based calls; PILOT4's 707 split into 431 and 276. All twelve small
primal/dual runs and both larger bounded runs match complete uninstrumented
SolveResult objects. These are branch/product counts, not cost shares. Do not
implement a replay-removal certificate for these models from the code path
alone. [l-replay-census.json](evidence/l-replay-census.json), raw source/probe
and outputs in `target/perf-l-replay-census/`.

The full forward U solve carries its numeric column and dirty blocks across calls. The current-host PEROLD profile above refreshes whole FTRAN/L and refactorization attribution. Inspect a whole representation or ownership boundary there before another local kernel trial; reproduce PGO/TCMalloc qualification separately when those builds are available. The occupied-slot BTRAN, touched-only pivot-row and PDHG packet/fused-step experiments below have completed their prescribed screens. Initial unit-weight setup is about 3.3% of full mode-switch time in the recorded baseline sample.

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

### EPYC structural trials

- Auxiliary column order by bounded stable length histogram/prefix/scatter,
  exactly preserving(length,incomingposition) and originalsmall/widekeyfallback.
  Numericfactor body literally unchanged;6release tests,fast940 and17201dumps pass.
  Ownfixed39PGO DEGEN2 counterpairs give−0.3951%I/−1.2937%C. Firstscreen75117:
  median−2.1457%,CI[−4.9610,+2.3988]%,17/32,p.4300; AB−4.6035/BA+2.0725%.
  HealthyAA,p.8601,CI[−4.5697,+4.3121] preventsseparation. INCONCLUSIVE; stopwithout
  75189,extension/retry/key/threshold/histogram/radix/allocator/codegenrescue orfull
  candidategates. Acceptedsource/binaries unchanged. Source/proof/oracle/tests,
  profile,counters,all130timingoutputs andcompletefailedscreen in
  target/perf-aux-length-order/; [evidence](evidence/degen2-aux-length-order.json).


- Completed level-L density handoff reuses the exact already-computed prefix and
  retains the original general eta/output suffix and WORK, leaving ties and
  incompatible cases on fallback. Independent complete-original frozen solver,
  captured135 full/needed+SC50A and generated/extreme/cache/MXCSR tests pass19/19;
  matching ownfixed39 PGO passes fast940 and17201 original dump records. The
  fixed eight-pair DEGEN2 counter gate gives−0.5606% wholeprocess instructions
  but+4.0016% cycles. REJECT before74117/74189; no screen, retuning, layout/buffer
  rescue or full acceptance gates. Main source/accepted binaries unchanged.
  Evidence: [degen2-l-handoff.json](evidence/degen2-l-handoff.json), complete
  implementation/oracle/proof/capture/tests/profile/results in target/perf-l-handoff/.


- Ordinary phase-one dead-extended slope pipeline: accepted optimized machine
  code executes the software-F80 subtraction at232ec8 on ordinary DEGEN2,
  confirmed by bounded hash-checked GDB witness (parenttool009da8db) and exact
  complete result/exit comparison. Seven guarded qa/qb updates and the unused
  partition argument preserve every hardware FP operation; frozen ratio,
  partition and both lane sums reconstruct literally. All12468live full-kernel
  calls,32completeSolveResults and3072generated both-mode/16MXCSR cases agree
  on result bits, passes, callbacks and FPstatus. Whole original breakpoint
  kernel is only1.361ms per80BAUsolve in CPU2 isolated ordinary-release replay;
  role ratios cross, so no controlled gain is established. Broadsolve_primal
  self16.34%cannot be assigned to this cold kernel. Reject without production
  integration,matchingPGO,formal screens/full candidate gates or const-mode
  rescue. Both qualified CLIs and runtime source unchanged. Evidence:
  [phase1-dead-extended.json](evidence/phase1-dead-extended.json), all source,
  frozenoracles, releaseELFs, rawcaptures/disassembly and tests in
  `target/perf-phase1-dispatch/`.


- Complete corrector trial executor: AVX256 holds primal/dual recurrences for two
  weights; an odd weight uses AVX128 without padded lanes. Original coordinate,
  upper and homogeneous order, exact weight/winner expressions and speculative
  FP-status effects are preserved. The complete original selector is frozen.
  Four release tests pass: 9,936 full-representation/mode cases, 192 captured
  checks from 12 original small-model inputs, opaque integrated dispatch and
  exact inactive/upper-hole/decline regressions. Fast parity passes 940/940;
  every bounded counter and fixed-screen full result/exit object agrees.
  80BAU3B instructions fall 5.869%. Fixed CPU 2 screen 48811 (32 AB + 32 AA,
  default cap 1000) is inconclusive: median -3.344%, CI[-8.732%,-1.222%],
  24/32 improving, p 0.00350 and both orders favorable, but healthy AA
  [-2.157%,+4.462%] defeats conservative separation. Stop without 48883,
  extension/retry or candidate full gates. Accepted source/CLI are restored
  exactly. [barrier-trial-executor.json](evidence/barrier-trial-executor.json),
  all source, fixtures, failed/corrected tests, instruction proofs and raw timing
  in `target/perf-trial-executor/`.


- Transfer the dense forward-L buffer through U with positive-zero normalization fused into spike compaction (current virtual EPYC host, ordinary ThinLTO/glibc). All 16 independent FTRAN/L/needed-FTRAN release tests and fast parity 940/940 pass; transfer arrays match expansion bit for bit on random and captured f64/F80 inputs. Four AB/BA pairs per model preserve all 32 complete results and exits. Instructions ADLITTLE +0.015%, DEGEN2 cap300 +0.270%, PEROLD dual +0.092%, PILOT4 primal +0.132%: no initial reduction. Counter collection overlapped fast parity, so cycle/timing changes are diagnostic. Rejected before matching-PGO, formal timing or full candidate gates; source and accepted CLI restored exactly. This closes the fused-normalization variant left by the earlier separate-pass transfer trial. Do not repeat either variant. [Evidence](evidence/dense-l-transfer.json); candidate, independent test extension and raw counters/results in `target/perf-dense-transfer/`.

- A shared incoming-basis representation supplied ordinary dual structural pivoting, pricing resets, singleton WORK counts and the independently pivoted auxiliary factor. The candidate passed the full 987 frozen runs, all eight ordinary dual outputs, 690,060 exact WORK records, 46 full LU input/output pairs and a 6,144-case auxiliary factor oracle. On 80BAU3B dual, time inside the complete refactor improved in all 16 matched pairs, saving a median 3.70% of parent full-solve time. The whole-solve median improved 8.56%. The required refactor contribution was 5.00%, so the candidate stopped at this pilot. Evidence in next/dual-basis-layout and next/dual-basis-layout-review.


- A complete allocation-lifetime census of the released compact primal solver observed 572,651 same-pivot nonrefactor allocations on Quad and 481,739 on 5959 dual. The genuine TCMalloc checkout calibration used the actual size/site mix. Reusing the identifiable pivot arrays and their helpers projected 3.19% on Quad, 3.22% on 5959 dual, 2.33% on FX and 1.94% on 5959 primal. Even crediting every identified short-lived helper across the solve projected only 3.08–4.80% before filling the new buffers. No production arena was admitted. Evidence in next/compact-pivot-lifetimes.

- A single barrier factor owner sent its original row-major coefficients to all outgoing Schur consumers and stored its solution coefficients in four-destination packets for both triangular directions. All 156 complete factors and 1,352 actual forward/backward vectors matched on the three original barrier rows. On the full 103-step GREENBEA solve, 5/9 pairs improved for a median child CPU saving of 2.50%. The independent 10% and seven-pair requirements failed. Evidence in next/barrier-factor-solve-owner.

- Adding the current real rho, ordered numerical pivot row, threshold counts and each preceding local pivot's column, weight and tau to the installed component factor history still left 21 conflicting results among 6,007 exact cross-block matches. Adding the absolute live factor update count gave 10,554 unique keys in 10,554 post-install pivots. The shorter key would change answers; the longer key provided no reuse. Evidence in next/mode-history-inputs.

- Controlled ordinary primal row updates and their Devex products ran in ordinary precision on all 7,271 accepted 80BAU3B and PEROLD pivots studied. Across 5.22 million active columns, 1,280 repeated the necessary numerical inputs for both products. Even treating all measured ranked products and Devex dots as free gives 8.34% and 7.07% gross bounds before validation. Evidence in next/ordinary-primal-owner.

- On full DEGEN3 primal-r the entering FTRAN took 21.76% of solve time across 5,232 calls. Although 2,988 raw right-hand sides repeated, only one had a previously seen identical ordered output. PEROLD primal-c had nine such recurrences in 3,490 calls, with an optimistic 0.00953% whole-solve reuse bound. For DEGEN3 pricing0 the entire method took 12.50%, including its 5.55% numerical BTRAN. None of 3,694 pricing BTRANs repeated even an unordered input support within its live factor version. Evidence in next/primal-full-scope.

- Across 4,500 physical component installations on mode-switch, the main and auxiliary factors never shared an entire local pivot order. Their common prefixes covered 59/145,484 auxiliary products. Three disjoint local Markowitz, fresh-U and auxiliary-column loops accounted for 174/1,637 full-solve samples, or 10.63%. Reaching a net 10% improvement from those three loops would require removing about 94% of their combined cost before paying for validation and materializing the ordered outputs. Evidence in next/mode-local-elimination.

- On three complete connected dual models the original U scatter and transpose solve loops occupied 7.38% on GREENBEA, 9.26% on DEGEN3 and 9.96% on PEROLD. Grouping four nearby mutable U slots retained the products and their order. On DEGEN3 it could merge 544,473 destination updates but would have to maintain panels containing 2,548,460 affected coefficient entries over the complete pivot sequence. No production panel was admitted. Evidence in next/connected-u-panels.

- A barrier triangular-solve plan stored four independent rectangle rows coefficient-major and updated the diagonal block in place. The original forward and backward Rust kernels agreed bit for bit on four captured numerical factors per model and 16 generated RHS vectors per factor. Excluding setup, complete triangular time on GREENBEA was 23.29% lower in seven paired batches; 80BAU3B measured 24.88% lower. The GREENBEA threshold was 40% before packet construction. The probe stopped there. These are kernel replay measurements. Evidence in next/barrier-triangular-program.


- A transposed matrix program evaluated independent connected ordinary dual pivot-row outputs together while preserving row arithmetic and updating nonbasic membership. It passed 987/987, 448 independent FP-control cases, 5,919,040 original row-output comparisons and 508,011 original WORK records. Seven paired GREENBEA dual-controlled runs measured at most 0.715% of the full solve saved by the native interval (four pairs favorable). The predeclared 5% minimum and five-pair requirement failed. Evidence in next/connected-pivot-program.

- A right-looking barrier factor delivered each completed front directly to its consumers in the original child order. Every packed numerical factor, pivot, reciprocal, parity offset and minimum ratio matched across 104 GREENBEA, 15 DEGEN3 and 37 80BAU3B factor installations. On all 103 natural GREENBEA steps the median CPU saving was 5.32%, eight of nine faster, below its prespecified 10% cutoff. The diagnostic stopped before production qualification. Evidence in next/barrier-newton-graph.


- A component-coordinate ordinary dual pivot prototype reduced full mode-switch requested zero/copy/fill ranges from 7,799,440,529 to 191,905,238 bytes (97.54%) at the actual pivot entry. Both f64/F80 crossed-handoff marker and frozen-factor/spill oracles agreed. All 11,393 pivot payloads, 202,772 chronological WORK records and frozen full 987 matched. The 12-pair full CPU pilot measured −12.18% with both orders favorable, but 8/12 improved against the prespecified 9 required. No formal comparison, confirmation or landing. The independent prescribed reference panel measured parent 0.901255, candidate 0.748702, and genuine patched Gurobi 0.138421 seconds (one warm-up and three measured solves per role). Every result and exit matched its original for both solvers. Complete observations and hashes are in mode_switch.dual-r/report.json. Evidence in next/ordinary-pivot-view.

- The solve-owned range certificate for generated ordinary PDHG kernels replaced the prior per-step full-coordinate scan. It passed 987/987, 480 exceptional-state comparisons and 3,960 consecutive expiry/restart tests. On full GREENBEA the first prespecified screen gained 11.27% (56/64 faster). Confirmation gained 6.72% (48/64), interval [-12.15%, -4.57%] against the original parent, missing both the required 10% and noise separation. All 516 full outputs/exits agreed. The independent 23-row engineering sweep measured 80BAU3B PDHG 15.6% slower than its parent, although it beat the patched reference in all three rounds there. The trial stopped after this confirmation. Evidence in next/pdhg-envelope-final.


On the corrected f2e71888 ELF, 42 native scoped full solves and 14 source-observed full solves all matched their original endpoints/exits. All 68 auxiliary numeric reads were fresh: 54 actual stall-objective FTRANs and 14 final BTRAN vectors demonstrably overwritten by active mirror results. Across all 14 ordinary rows, even granting all dense/raw updates, complete eta retirement and all auxiliary triangular solves for free accounted for only 0.19–1.12% of the whole solve. The separately observed 572 constructors still make their immediate singular decisions. A new signed-incidence full-basis certificate had only 2/572 complete bases even eligible for its liberal column-support pretest, together 0.014% of the one affected full row. No separate solver patch.

The existing whole-factor-generation arena covered every small allocation in all 15 complete mode-switch refactors and measured 5.61% faster solver time, 4.89% faster child CPU, versus its required 8%; connected PEROLD measured 11.38% slower. A fresh 40-solve lifetime clock panel on f2e71888 bracketed all residual larger-than-256-byte TCMalloc requests, their later reallocation/copy and freeing, with all endpoints matching and zero outstanding allocations. Their entire measured underlying cost was 1.24% of mode-switch, 2.15% of PEROLD, 1.03% of GREENBEA and 0.80% of 80BAU3B. All-family Rust-side native drop-glue was independently 0.73% of mode-switch. Even treating both extra costs as entirely removable gives an indicative 7.58% mode-switch projection before any new backing, clearing, copying or ownership overhead. No new arena experiment.

The complete dual-rebuild owner was stopped on the corrected 7dba613b release. On full 5959 dual, 7,742 fresh and 15,474 nonfresh keys repeated exactly. The proof validates the immediately preceding BTRAN and actual resident reduced-cost owner, restores all d and only the written df, executes the native last column on updated factors, and retains the live phase, tolerance/counter and WORK decisions. The new fixed39 TCMalloc executable 5d61e75f passed 987/987, both four-row WORK modes (1,319,417 replay-active plus 1,377,120 dump events) and eleven alternate FP settings. The prespecified seven full paired native-scope comparisons counted 23,391 BTRANs and 23,385 post-recomputes per solver per run. Post-recompute medians were 9.435 ms before and 3.780 ms after; BTRAN medians were 5.451 ms before and 5.604 ms after. Median paired aggregate scope savings divided by the full parent solve were **4.75035%**, below the frozen **5.00%** threshold. The individual pair estimates ranged from −12.91% to +10.06%. Only this mechanism comparison ran. The branch did not advance to formal timing or landing.

The independent release-lib shadow also passed on the three complete 5959, quad and fx trajectories. It exercised 40,793 actual restores, reran and checked 418 original numerical post-recomputes (including full d/df, tolerance/counter, sum, x87 CW/SW and MXCSR) and independently checked 429 original WORK recomputations.

- Across 56 original ordinary solves covering all 14 rows, the combined dense/sparse base-L operations accounted for 5.82–11.50% of solve time. Building a replacement two-direction layout cost 0.95–2.48% per full solve. Even an optimistic per-install selector that removes all indicated triangular arithmetic, scheduling and WORK had a timer-corrected 9.51% ceiling on the most promising row, 80BAU3B dual. The dense-only share never exceeded 9.81%. The lower-only proposal stopped before production.

- Full compact post-pivot row geometry and candidate-list reuse passed complete output and WORK checks. Its primary nine-pair diagnostic put whole 5959 dual CPU 5.34% lower with seven favorable pairs. Independent inclusive instrumentation over the original row producer, factor update, complete list transition and next price call attributed only 3.17 ms of 155.31 ms, or 2.04% of whole-solve CPU, to the change. It stopped before PGO and formal comparison.


- The compact pivot direction prototype retained exact FTRAN successor state, ordered/signed binary64 and extended outputs and both dense representations. The point update reused its checked numerical result for repeated inputs. Across four original model solves, re-executing the original verified 46,594 directions, 46,591 point steps and 331 complete successor-factor snapshots. Both WORK modes and 987/987 frozen cases agreed. The primary full fx viability test changed by -3.33%, interval [-14.40%, +4.31%]. It stopped below the 5% minimum.

- A matrix-specific PDHG program emitted ordinary sparse products in their original arithmetic order at solve start. A same-build feasibility comparison was 34.23% faster for GREENBEA and 26.43% for 80BAU3B. Production inserted mandatory finite-domain checks before any state mutation and preserved the earlier operator on exceptional inputs and protection failures. Full parity passed 987/987, the exact production/parent exceptional-state audit matched 552,480 numerical values and 480 floating environments, and all 59,679 eligible ordinary benchmark steps took the compiled path. The fixed39 production test on original GREENBEA measured -1.97%, interval [-6.42%, +6.56%], with 18/32 pairs faster. It stopped without confirmation or a release commit.


- Ordinary dual pivot transaction, measured on 2731b202. Solve-only native stacks, exact ELF-scoped interposition and 24 complete phase profiles separated earlier pivot-frame, secondary FTRAN, output-order/WORK and sparse-level costs from the remaining pivot handoffs. The earlier regions cover 97.8% of 7.799 GB of mode-switch transfer requests and 92.4% of 1.047 GB on 80BAU3B. PEROLD has only seven iteration libc leaves in 990 solve samples. Remaining packing, alignment, flip preparation, basis changes and row canonicalization cost 1.65%/2.67%/2.77% on mode-switch/PEROLD/80BAU3B. The full ratio test, basic/Devex update and reduced-cost update still perform required arithmetic. Even granting all three for free in addition to free handoffs gives 5.31%/9.92%/10.78%, before owner bookkeeping. The ordered ratio-coefficient vector repeated 0/11,393, 1/1,185 and 9/5,185 times respectively, disregarding additional live cost dependencies. No dense/joint FTRAN outputs occurred in mode-switch. No qualifying whole-pivot rewrite was implemented. Complete raw controls, outputs, source lines, scripts, and disjoint estimates: next/pivot-transaction/review.md.

- On 80BAU3B barrier, 12,022 to 12,061 of 12,061 input weights changed between adjacent factors. In 20,052 front transitions, no front retained the full set of its own required diagonal inputs. Across different rows in one factor, all 1,503 DEGEN3 diagonal supports and all 2,262 80BAU3B diagonal supports were distinct. The additional ordering call in the CLI runs after the end of the solve timer. These checks did not support reuse of complete numerical fronts or credit for removing the post-solve ordering call.

- A PDHG program classified homogeneous operand-address sequences and common projection bounds. It matched both original full outcomes and focused exceptional cases. Nine complete GREENBEA paired trials gave 0.18% less CPU time, with five of nine improving. The exact-ELF primal phase, 22.49% of solve time, was 1.40% slower in the separate phase test. The production build and formal comparisons were held.


- Ordinary dual factor storage across epochs, 2731b202. Solve-only native stacks put dual refactor at 42.6% of mode-switch and 30.5% of PEROLD dual. A six-model exact row/capacity census and two actual allocator-pointer campaigns separated typed storage from numerical work and unrelated allocations. The 17 studied families include split basis inputs, structural output, main mirror and transpose, auxiliary L/U and pricing reset. The allocations, reallocations and frees of buffers fitting an alternating same-slot bank, including every previously untagged retiring buffer, accounted for an estimated 7.02% of the entire mode-switch solve (95% solve-bootstrap 5.78–8.35%) and 2.33% of PEROLD (1.55–3.24%) before owner bookkeeping. The complete local numeric factor never repeated in 4,200 between-refactor block comparisons under physical coordinate translation. Identical complete physical Markowitz rows covered 9.0% of its local cross products. All diagnostic results and exits were canonical. No production patch or new speedup claim. Evidence and scripts in next/refactor-slot-lifetime and next/refactor-retired-slots.

- A joint PDHG layout retained both primal and dual coordinates between ordinary iterations and evaluated four long outputs together. It passed the original endpoints and focused checkpoint/restart and exceptional-value checks. On nine matched complete GREENBEA pairs the median was 1.77% faster and six pairs improved. It stopped below the prespecified 5% minimum and seven-pair requirement. On separate captured late ordinary steps, the complete inputs to 8.97% of forward terms and 0.83% of transpose terms were bit-identical to the preceding step. The corresponding 80BAU3B fractions were 7.80% and 1.65%.


- Letting each compact factor state own one complete pricing-row result preserved full 987/987 parity, all four complete outputs and both WORK modes. Re-executing every hit verified all changed arrays and flags. Only 864 of 5,205 updates hit on the prespecified fx target. Its 24 mixed and 24 same-binary pairs measured +0.36%, 95% interval [-2.55%, +5.41%]. The trial stopped there.


- A basis-ordered feasibility program grouped four simple bounds while preserving the original loss-of-feasibility decisions. On 80BAU3B 87.81% of visited positions fit the four-wide groups. The complete captured classifier was 25.81% faster across nine paired replays, with eight pairs favorable. In ten full solves of the current compact LU candidate, 166 of 1,622 in-solve samples belonged to those checks. The projected whole-solve gain was 2.64%, below the 5% threshold, so the experiment ended before a production build. All three captured model results and original classifications agreed. The exceptional scalar fallback passed 32,768 numerical and MXCSR cases plus six unmasked-trap or early-exit cases.


- Reusing the natural numerical traversal as the sparse FTRAN U output order reduced GREENBEA dual allocations by 62,141 and preserved every output in 987 parity runs. Its 160 mixed and 160 identical-binary pairs gave −3.19% [−5.21%, −1.18%] and a same-binary interval of [−4.62%, +1.53%]. The intervals overlapped. No second comparison ran.
- Matching disconnected main-factor blocks by normalized pivot rank, exact input coefficient bits, ordered lower factors and Markowitz entries covered 8,670 of 69,863 numerical U products on mode-switch. All matching complete numerical and ordered virtual U outputs agreed. The observed matching production time was about 1.0% of the instrumented solve before paying any matching or remapping cost. Five ordinary controls were all below 0.003%.
- Storing extended values in their native ten-byte encoding across whole factor lifetimes and using guarded native arithmetic in FTRAN/BTRAN passed full 987/987. Distinct software and native implementations matched on 300,060 operand pairs and on all 576 floating-environment checks including the full occupied x87 register stack. The predeclared complete fx_quad_obj_5371 pilot was 35.53% slower. Both presentation orders slowed. No formal timing followed.
- Sharing the ordinary lower factor, its transpose and pricing-reset predecessor geometry saved 211,138 allocations on the full 80BAU3B dual solve. Full 987/987, the public mutable-pivot fallback and all recorded factor/solve checks passed. Its first 96-pair comparison was −1.78% [−5.92%, +0.50%]. No confirmation followed.
- A single generation-marked frame for the sparse dual pivot row, Devex framework RHS and tau column passed full 987/987 and all 202,772 mode-switch work records. On the complete mode-switch witness its fixed 16-pair viability trial measured −3.09% [−14.52%, +9.35%], with 9 of 16 faster. It missed the 5% target and the second order did not improve.


- A main-factor buffer restore with separate shared row/work installation measured −6.59% and −7.71% on the released c82803ab solver. The second identical-binary control was biased (p=.0146). Qualification stopped. A separate parked auxiliary-factor cache missed the noise criterion on its own parent. The active compact installer tests a different ownership lifetime.
- A fresh U producer merged the ordered live numerical frontier with the Markowitz overlays. It preserved 987/987 and all 642 original full mode-switch results and exits. The 160-pair comparison measured −3.04%, interval [−6.09%, −0.47%], versus an identical-binary interval [−3.66%, +2.78%]. Its sign p=.01079 missed the fixed threshold, and the intervals overlapped. No confirmation ran.
- A certified integer interval for the two barrier complementarity roundings matched all 987 frozen runs, all 159 original barrier records and an independent raw-bit and floating-environment audit. Full 80BAU3B barrier was 11.33% and 11.28% faster than this rose parent in independent 384-pair comparisons. The direct CPU0 medians were 0.288851 seconds for the experiment and 0.154557 for patched Gurobi. The certificate skips the intermediate rounded Kahan stages and falls outside the operation-order rule in PLAN.md. It remains an experiment. A separate exact software implementation reproduced the original intermediate values but its integrated component projection of 3.48% missed the predeclared 4% threshold.


- Retaining the previous auxiliary LU base by exclusive ownership and discarding its eta history on reinstall passed 987/987. The two full fx_quad_obj_5371 comparisons measured −5.36% and −3.72%. The latter interval [−5.01%, −1.37%] did not clear the identical-binary error band, so no commit followed.
- Maintaining global primal pricing score buckets avoided 2,612,079 of 4,514,341 observed full/best/old-list score evaluations on 80BAU3B. Full parity and the 76,269-record target trace matched, but its 96-pair full medians were +1.13%, interval [−2.40%, +3.05%]. Both presentation orders were slower. No confirmation.
- A direct-coordinate Markowitz producer passed 987/987 and its independent sparse/dense/duplicate oracle. On PEROLD primal the full 12-pair viability trial regressed 3.81% (3/12 wins). The four-RHS auxiliary panel also preserved its independent complete factors, but cost 9.22% more on the PEROLD dual constructor replay. Neither continued to formal timing.
- An end-to-end barrier factor arena with direct child and solve views reduced allocation calls 68.9% on GREENBEA and 58.2% on 80BAU3B. It reproduced all recorded child updates and factors. The median whole-process CPU projection was 4.658% on GREENBEA and 2.045% on 80BAU3B, below the frozen 5% gate. No production integration.


- Direct dense dual-row accumulation removed the membership scan and passed its f64/F80 oracle, full 987 and PEROLD WORK trace. Its 192-pair target comparison measured −1.63% [−4.19%, +1.28%]; the 107 improving pairs did not establish a gain.
- Tracking ordinary basic feasibility classes through sparse pivots passed 987/987 and a per-observation complete-rescan oracle on two models. At 80BAU3B the full 96-pair median was −10.26%, interval [−17.19%, −5.08%], versus an identical-binary band reaching +6.30%. No independent second screen ran.
- Owning the sparse primal pivot column through eta installation preserved the full 5,022-pivot trace and removed 92,008,555 requested bytes on 80BAU3B. The fixed 24-pair viability comparison measured −1.66% [−12.98%, +8.61%], below its required gain and significance. No formal trial or landing.


- Width-one barrier leaves removed 492,128 dynamic allocations over the genuine 103-iteration GREENBEA benchmark. Full parity, independent arithmetic checks and all 1,026 exact-output samples passed. The fixed comparison measured −3.048% with interval [−4.749%, −1.210%] and identical-binary interval [−2.580%, +2.074%]. The intervals did not separate, so no confirmation ran. The two-column factor/triangular owner is a separate trial.
- Reusing constant patterns across adjacent batches of the ordered transpose passed 987/987. The bounded GREENBEA pilot improved by 5.38%, but the full comparison had 1.576651 s for the new code, 1.458054 s for its ordered parent and 1.314887 s for patched Gurobi. It did not qualify.


- Build a prescribed-pivot auxiliary validator and defer materializing its full factor. PEROLD primal. The full corpus and differential checks passed. The comparison on its designated parent measured -2.01%, 95% [-5.42%, +3.19%], identical-binary interval [-2.99%, +4.73%]. The prespecified test was inconclusive.

- Share immutable compact-basis installation rows and pivot structure across installations. primal_5959 dual. The full corpus and differential checks passed. The comparison on its designated parent measured -4.82%, 95% [-7.41%, -0.38%], identical-binary interval [-0.74%, +2.88%]. The prespecified test was inconclusive. The previous source had passed independent comparisons at −5.03% and −4.90% against the indexed predecessor. The result on the BTRAN parent governs landing.

- Sparse secondary FTRAN retained its numerical lane and visited the touched slots for WORK. Full frozen parity passed. Its 96-pair full mode-switch comparison measured −3.29% [−10.77%, +0.88%] with identical-binary interval [−4.94%, +5.46%] and sign p=0.131. The independent confirmation was not run.
- Component-local dual row, framework RHS and pivot-column cross-term storage together passed parity. The first full comparison measured −8.56%. The independent median of −9.86% did not separate from its identical-binary interval, so the combination was set aside.
- Computing the sparse forward product directly during the stationary-column PDHG scatter preserved the captured trajectory, all four focused tests and 987/987. The 8,000-iteration pilot was 29.97% slower than its stationary parent and stopped at its preregistered rule.


- Compiling sparse forward PDHG rows from exactly zero reflected coordinates kept their original lanes, reduction tree and scalar tail. At runtime a nonzero input invalidated only dependent rows. The captured-matrix differential and full frozen 987 passed. The fixed 8,000-iteration 80BAU3B pilot regressed 17.92% against the dense-state source, triggering the prespecified stop. No formal comparison or commit followed.

- Component-local dual row and framework RHS alone passed 987/987. The two full mode-switch comparisons measured −7.37% and −7.11%. The latter interval [−10.18%, −3.84%] did not separate from the identical-binary interval [−3.64%, +5.25%]. The combined basis-slot experiment above has a distinct second-column producer and consumer.

- Interning exact ordered normal-equation assembly expressions on DEGEN3 preserved frozen parity and complete results. Its 1,024-pair comparison measured −1.71%, interval [−2.43%, −1.0878%], identical-binary interval [−1.0911%, +0.26%]. It missed the predeclared noise separation and did not qualify.



- Removing the second dense FTRAN result expansion and aligning the two handed-off columns. Full frozen parity passed. The first fixed comparison measured -2.58%, interval [-8.65%, +0.85%], AA [-0.63%, +5.31%]. The predeclared timing test was inconclusive.

- Keeping sparse dual pivot rows in packed first-touch order. Full frozen parity passed. The first fixed comparison measured -5.06%, interval [-7.98%, -0.67%], AA [-3.50%, +4.86%]. The predeclared timing test was inconclusive.

- Replacing gathers with vector loads on contiguous PDHG matrix intervals. Full frozen parity passed. The first fixed comparison measured -4.27%, interval [-7.51%, -1.61%], AA [-4.77%, +0.87%]. The predeclared timing test was inconclusive.

- Avoiding barrier iteration history in final-result-only solves. Full frozen parity passed. The first fixed comparison measured -2.88%, interval [-4.23%, -1.74%], AA [-3.13%, +0.23%]. The predeclared timing test was inconclusive.

- Link live Markowitz coefficients in both row and column order, eliminating front shifts and stale-chain cleanup. PEROLD dual. The full frozen corpus passed. The first randomized timing comparison measured -2.54% with a 95% interval of [-5.83%, +0.09%]. Same-binary interval [-0.83%, +6.91%]. The predeclared screen was inconclusive.

- Reuse the primal row-product workspaces and visit reachable columns for sparse BTRAN support. DEGEN2 primal. The full frozen corpus passed. The first randomized timing comparison measured -3.65% with a 95% interval of [-5.33%, -1.58%]. Same-binary interval [-2.70%, +1.46%]. The predeclared screen was inconclusive.

- Guarded native x87 Kahan products for the barrier. This GREENBEA trial explicitly set 500 and ended at the iteration limit, whereas the unmodified benchmark ends suboptimal at 103. The full frozen corpus passed. The first randomized timing comparison measured -4.77% with a 95% interval of [-10.20%, +0.71%]. Same-binary interval [-3.31%, +5.81%]. The predeclared screen was inconclusive.

- Compute only consumed entries in eligible child frontal updates. 80BAU3B barrier. The full frozen corpus passed. The first randomized timing comparison measured -7.31% with a 95% interval of [-13.59%, -1.53%]. Same-binary interval [-5.04%, +15.78%]. The predeclared screen was inconclusive.

- Touched-only factor update storage. The completed screen gave -25.95% median, 95% interval [-33.24%, +6.28%], AA [-12.50%, +3.43%]. The candidate remains unqualified.

- In-place barrier substitution. This GREENBEA trial likewise used the 500-iteration endpoint. The completed screen gave -0.24% median, 95% interval [-1.78%, +2.44%], AA [-2.45%, +5.61%]. The candidate remains unqualified.

- Cross-row repeated product and prefix reuse on full GREENBEA PDHG. The completed screen gave -12.61% median, 95% interval [-18.09%, -6.96%], AA [-3.80%, +12.07%]. The candidate remains unqualified.

- Equal-length AVX2 sparse operator groups with long-row peel scheduling. The completed screen gave -12.93% median, 95% interval [-20.00%, -2.21%], AA [-12.06%, +5.35%]. The candidate remains unqualified.

- Sparse FTRAN occupied-block traversal alone preserved 987/987 and the complete solver result. The first capped mode-switch screen passed. The second gave -17.88% paired median, 95% interval [-20.61%, -10.07%], AA [-10.18%, +7.82%], 27/32 pairs improving. The conservative noise-separation criterion failed. The subsequent column-workspace experiment also changes numeric buffer ownership, reset and output bookkeeping.

- Sparse BTRAN U-transpose workspace preserved 987/987 and the exact full mode-switch WORK trace. The fixed mode-switch screen gave -8.05% paired median, 95% interval [-10.84%, -4.78%], AA [-5.87%, +3.37%], 27/32 pairs improving. The conservative noise bound failed.

- Dual pivot-row accumulation reused coefficient and membership arrays and reset only touched coordinates. Full parity passed. The mode-switch screen gave -5.31% paired median, 95% interval [-12.54%, +4.46%], AA [-3.30%, +5.32%], 18/32 pairs improving. The interval includes zero.

- Packed short-row PDHG grouped independent equal-length rows in four lanes in both sparse multiplications. Full parity passed. The GREENBEA screen gave -14.32% paired median, 95% interval [-22.17%, -7.28%], AA [-4.02%, +23.02%], 28/32 pairs improving. The same-binary interval was too wide to qualify the change.

- Initial dual weights from structurally isolated unit right-hand sides reduced a prior-profile mode-switch zero-pivot pilot by 41.6%. The four complete-run ratios crossed in both directions and the full median increased 5.5%. The original entire weight-initialization loop was about 3.3% of this solve. The preliminary source was left out of the solver.

- Auxiliary LU ownership audit found no identical complete constructor inputs on PILOT4 primal (17 constructors) or PEROLD dual (27). Canonical identical independent components could reuse zero of 223,269 elimination products on PILOT4 and 81 of 1,081,134 on PEROLD dual. Deferring construction still has to preserve the exact singular-factor decision.

- PDHG ordinary-step fusion combined each sparse matrix multiply with the projection and Halpern update. Checks, restarts and limit-boundary steps kept their original paths. Full parity passed. The fixed GREENBEA screen gave -9.42% paired median, 95% interval [-15.06%, -2.20%], AA [-9.83%, +11.96%], 24/32 pairs improving. Its interval overlapped the same-binary uncertainty. No confirmation was run.

<!-- end EPYC structural trials -->

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

- Ordinary positive weight-floor specialization: the floor is always finite and strictly positive (W_MIN or1), allowing register-only SSE2 maximum with original NaN selection; the extended path keeps its existing max. The independent complete row oracle is byte-identical; reversal changes only the one ordinary floor selection. A frozen runtime primitive passes2,097,696 full-bit/boundary/quiet-and-signaling-NaN comparisons across all16 MXCSR modes; the complete captured/random full-row test and fast940 pass with matching fixed39 PGO. Eight controlled ADLITTLE/BORE3D and three release PILOT4 complete result/exit pairs agree. Instructions fall0.143/0.050/0.381%; PILOT4 cycles+0.221% are effectively flat. Its predeclared32-AB/32-AA PILOT4 screen (CPU4, cap 800, seed26217) is inconclusive: median-0.952%, interval[-3.249%,+2.688%],17/32 improving, sign p.430, healthy AA[-2.561%,+1.320%], execution orders-1.820%/+1.342%. Stop without confirmation26289 or full candidate gates. No qualifying runtime gain. Source/accepted CLI/metadata/profile restored exactly; restored formatting, fast940 and full987 pass. Native accepted-row annotation, complete source/oracle/proof, builds/tests/freezes, bounded results/counters and controlled screen: `target/perf-weight-floor/`; accepted annotation is also retained in `target/perf-rank-reuse/primal-row.asm`.

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

- Compiled standard barrier assembly into ordered scalar contribution destinations preserves the complete frozen assembly oracle on all fronts of AFIRO/SC50A/BLEND/ADLITTLE/BORE3D, including live coefficients/scaling and single/double rounding. Construction costs 0.653 ms across the small models versus 0.650 ms saved by 100 repeated assemblies. Whole-solver instructions increase on all five measured models (+0.235% to +1.491%); DEGEN2 cycles increase 9.63%. All 40 complete result/exit objects agree. Rejected before formal timing/PGO integration. Production source restored; full parity passes 987/987. Evidence: [compiled-assembly.json](evidence/compiled-assembly.json), raw candidate/oracle/results in `target/perf-compiled-assembly/`.

- Fuse fresh-factor structural L closure, cloning/sorting and numeric replay into an ascending write-discovered bitset frontier, certified strictly lower over every L entry; nonlower graphs keep the original closure. The complete original constructor oracle passes all three release tests (generated f64/F80, small-model factors, all 16 SSE modes, directed cleanup/drop/cancellation/duplicate/overlay/nonfinite cases, full fields/U bits/WORK/flags). Fast parity passes 940/940. All 48 complete bounded counter results/exits agree. PILOT4/PEROLD instructions fall only 0.387/0.254%; cycles move +3.84/-7.55% in three diagnostic pairs. Fixed ordinary-release PEROLD screen42117 (CPU2, cap1300, 32AB+32AA) is nonpositive: median +3.925%, interval [-1.954%, +5.371%], 12/32 improving, both run orders slower; healthy but broad AA interval [-1.713%, +3.661%]. Stop without42189, extension or full candidate gates. Source and accepted CLI restored exactly; restored fast parity passes940/940. No optimization accepted. Evidence: [fresh-numeric-frontier.json](evidence/fresh-numeric-frontier.json), raw source/oracle/proof/build/tests/counters/schedule/results in `target/perf-fresh-frontier/`.

- Generated complete immutable dense forward-L and L-transpose SSE2 programs with fixed workspace displacements and exact coefficient/WORK literals. All 8,192 captured SC50A/mode/input cases and 4,888 actual dense-pass inputs across eight bounded models preserve numeric/WORK bits and MXCSR; all complete solver results match uninstrumented runs. Execution-only tiny-factor probes improve, but construction-inclusive lifetime replay loses. Histories optimistically share all identical geometry/coefficients within one solve and exclude capture/cache costs. A one-pass mapped-image emitter removes staging/relocation and produces byte-identical code on all 182 observed geometry-direction groups, but final CPU2 replay still loses about 63% on PILOT4 and 29% on PEROLD; each aggregate direction also loses. Do not integrate this per-LU-factor program or infer solver gains from repeated-kernel timing. Source restored exactly; no production patch/full candidate gates. Evidence: [l-generated-program.json](evidence/l-generated-program.json), complete probes/source/checks/replays in `target/perf-l-program/`. Programs with topology retained across all barrier Newton factors would have a different, untested lifetime.

- The isolated four-double-RHS rectangle prototype is numerically verified and
preserved on `experiments/barrier-rhs-packets`, but remains unqualified for
acceptance. Independent checks cover5,624rectangles and1,488complete blocks;
full frozen parity passes987/987 and285release tests pass (22ignored). Original
G4 schedules, live coefficients, separate nonfused operations, scalar tails,
nonfinite replay and unmasked-environment fallback are retained. DEGEN3
instructions/cycles fall9.224/14.029%. First fixed screen44217 (CPU2,cap100,
32AB+32AA) measures-10.103%, interval[-13.000,-4.208],27/32improving,p5.65e-5,
bothordersimprove; healthy AA[-4.300,+2.931] defeats conservative separation.
No44289,extension or same-source retry. This isolated packet experiment was not accepted.
[rectangular-rhs-packet.json](evidence/rectangular-rhs-packet.json), exact
implementation/proofs/full checks and every timing result in
`target/perf-rhs-packet/`. The combined factor executor is qualified above;
these isolated measurements remain inconclusive and are not added together.

- Retain complete barrier child execution row schedules and temp/product scratch across Newton factors. Geometry includes both width-one orders and wide GEMM calls, exact logical slice extents, contiguous gather/store versus zero/add, 24-column groups and rectangle extension. Live coefficients, D and precision stay original; public numeric_factor retains its fallback. All 7,200 independent child comparisons and 30 complete small-model factor/callback cases pass, including flags and single-to-double poisoned-scratch reuse. Fast parity940 passes and all42final bounded result/exit pairs agree. DEGEN2/DEGEN3 instructions fall4.803/5.246%, cycles10.404/3.831%. The fixed ordinary-release DEGEN3 screen43117 (CPU2,cap100,32AB+32AA) is inconclusive: median-6.561%, interval[-9.205,-2.991],26/32improving,p0.0002675,bothordersimprove; healthy AA interval[-1.949,+4.326] defeats conservative separation. Stop without43189,extension or full candidate gates. Source/acceptedCLI restored; no optimization accepted. Evidence: [barrier-child-program.json](evidence/barrier-child-program.json), complete implementation/oracle/tests/builds/counters/screen in `target/perf-child-program/`. The combined factor executor is qualified above; this isolated schedule screen remains inconclusive.

## Tried, breaks parity

None.
