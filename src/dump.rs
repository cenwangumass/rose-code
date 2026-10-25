//! Diagnostic dumps for comparison with the reference captures, selected with `ROSE_DUMP=kind,kind` (or `all`).
//!
//! Each event is one line on stderr, `TAG {json}`, with floats written to round-trip exactly. The tags and record fields
//! follow the reference captures of `reverse-engineer/tools/refcap.py`; `reverse-engineer/tools/compare.py` runs rose and
//! compares the two. Add a kind here rather than an ad hoc `eprintln!` behind its own variable.
use std::sync::OnceLock;

/// Dump kinds, their line tags and what they record.
pub const KINDS: &[(&str, &str, &str)] = &[
    ("status", "STATUS", "once per solve, written by `examples/trace_dump` whatever ROSE_DUMP selects: final status, iterations, objective"),
    ("dualacc", "DACC", "dual: each x_B recompute (`0xc1861`): factorization count, request, R = max(1e-9, max|rhs|), max|x_B|, residual e2 (0 when not read), the threshold candidate and `ctx+0x390` before it"),
    ("header", "PHDR", "once per solve: n, m and the fixed variables (lo == up)"),
    ("pricing", "PLIST", "each pricing call: call index, candidate list, d, w, basic set and x"),
    ("scale", "SCALE", "once per solve: the scale factors, columns (structurals, then range columns) and rows"),
    ("pricing0", "P0", "`pricing0` scan: list, position, duals and stored d at entry (tag in), and the chosen variable with its d (tag out)"),
    ("ftran", "FTR", "pivot-column FTRAN input rows and result in output order"),
    ("ustage", "USTG", "U-stage state of each pivot-column FTRAN: the L-stage result (spike, ascending), the diagonal, the tail columns' U entries and slots"),
    ("lstage", "LST", "dense x (by position) of each pivot-column FTRAN after the L stage (l) and after the row etas (e)"),
    ("rho", "BTR", "pivot-row BTRAN result in output order"),
    ("tau", "BTT", "Devex cross-term BTRAN result in output order"),
    ("dual", "BTY", "dual BTRAN after a factorization, in output order"),
    ("ord", "ORD", "pivot-row input list and changed columns per pivot (CORD for the cost-change pass)"),
    ("xb", "XB", "right-hand side and x_B of the basic-solution recompute"),
    ("fstruct", "FSTRUCT", "U-factor structure (kinds, nucleus list, tail start, slots) before each column replacement"),
    ("luin", "LUIN", "basis column order handed to each factorization (the reference's lu_in `cols`), and the order it returns"),
    ("luout", "LUOUT", "numerics of each fresh factor: the diagonal and the L columns (the reference's lu_in `diag`, `Lcp`/`Lidx`/`Lval`)"),
    ("ent", "ENT", "phase 1 ratio-test inputs: the signed FTRAN list, and distance, range and type per basis position"),
    ("ent2", "ENT2", "phase 2 ratio-test inputs: the signed FTRAN list, and distance to the lower bound, range and free flag per basis position"),
    ("bsolve", "BSOLVE", "basic-solution FTRAN (after a factorization): L-stage result, x after the tail columns, final x, kinds, slots, diag"),
    ("switch", "SWITCH", "method decision after each refactorization (`0x142500`): the counters it reads and whether it switched to Devex pricing"),
    ("dcr", "DCR", "the dual start basis: the (column, row) assignments of the crash in order"),
    ("dsx", "DSX", "each call of the dual iteration function: state at entry (head, x_B, d in the reference frame, W, status, flags, priorities, candidate list, threshold, mode, cursor) and the decision (row, leaving variable, entering variable, flips, pivot row, rho)"),
    ("dsxl", "DSX", "lean variant of dsx for whole solves: per call only k, the leaving row r and variable, the entering variable q and the bound flips (`compare.py trace --lean`)"),
    ("dtau", "DTAU", "each dual pivot at the weight update: leaving position, alpha, w (rho norm), the pivot column without the pivot entry (position, value) and tau at the same positions"),
    ("p68", "P68", "each reset of the pricing weights (`0x132460`): the type, the U and L column structure (steps) and the weights P+0x68 by position"),
    ("work", "WORK", "the work counter W (`crate::work`): each increment (k inc: d, site, caller), each refactorization check (trig: W), the time clause (clause), and the smoothed factorization cost (w158); `ROSE_WORK=1` also selects it"),
    ("pdhg", "PDHG", "PDHG: the state at each residual test (every 100th iteration, and the iteration-limit exit): counters, eta, omega, merit R, the residual scalars and the status"),
    ("barrier_struct", "PAT", "barrier setup structure (`examples/barrier_struct`): records PAT (AA' pattern), SYMIN (n1, dense count), ORDIN (level flags), ORD, ORDFIX, SYM, and BSTRUCT (status: `unsupported` names an undecoded routine)"),
    ("barrier", "BARRIER", "barrier iteration (`examples/barrier_replay`): per iteration the state at its start (x, z, y, u, v, mu, objectives, max residuals), the predictor steps ap/ad, mu_aff, target, the selected steps tp/td with their weights, gamma and the termination code"),
    ("barrier_adapt", "AD", "release barrier: per iteration the work-window totals (`w138` factor, `w130` predictor, `w148` first corrector, `w140` later correctors, `w128` update), the ratio, the mean step product, the corrector limit before and after, and the trial weights of the next iteration (shaped like `replay/gdb_barrier_adapt.py`)"),
    ("barrier_factor", "FSN", "barrier numeric factor (`crate::chol`), per factor f (0 = starting point, then one per iteration) and supernode in column order: FSN (f, col, w, c, row list) and FOUT (f, col, packed factor values, pivots), shaped like the captures of `replay/gdb_barrier_chol.py`; `compare.py chol --rose`"),
    ("barrier_best", "BEST", "barrier best-iterate tracker, one record per call (flag 0 per iteration, 1 at the end): it, q88 (saved iterate), qa8, the running minima, Q+0x220 / 0x230 (norms of the updated iterate), the history arrays h48..h80 at it and at q88; `compare.py barrier_best`"),
    ("barrier_solve", "RSOLVE", "barrier solve of the normal equations (`0x305f10`), per call s (counted from 0, the starting point first): rx, ry, D and the results dx, dy as hex bit patterns, shaped like the SOLVE / SOLVEOUT records of `replay/gdb_barrier_solve.py` (rose orders the solves of an iteration differently, so match them by input)"),
    ("fillparts", "FILLP", "per factor update: pool fill on entry, eta size, new column entries added, and the row growths (row, cnt, cap, extra)"),
];

fn kinds() -> &'static [String] {
    static K: OnceLock<Vec<String>> = OnceLock::new();
    K.get_or_init(|| {
        std::env::var("ROSE_DUMP")
            .map(|v| {
                v.split(',')
                    .map(|s| s.trim().to_string())
                    .filter(|s| !s.is_empty())
                    .collect()
            })
            .unwrap_or_default()
    })
}

/// Whether the dump `kind` is selected.
pub fn on(kind: &str) -> bool {
    debug_assert!(
        KINDS.iter().any(|k| k.0 == kind),
        "unknown dump kind {kind}"
    );
    kinds().iter().any(|k| k == kind || k == "all")
}

/// Whether any dump is selected.
pub fn any() -> bool {
    !kinds().is_empty()
}

/// Writes one record of `kind` under its tag.
pub fn emit(kind: &str, rec: serde_json::Value) {
    let tag = KINDS.iter().find(|k| k.0 == kind).map_or("DUMP", |k| k.1);
    // One write per record: stderr is unbuffered and `Value`'s Display writes piece by piece.
    let line = format!("{tag} {}\n", serde_json::to_string(&rec).unwrap());
    std::io::Write::write_all(&mut std::io::stderr().lock(), line.as_bytes()).unwrap();
}

/// A sparse vector as `[[index, value], ...]`.
pub fn sparse(v: &[(usize, f64)]) -> serde_json::Value {
    serde_json::Value::Array(v.iter().map(|&(i, x)| serde_json::json!([i, x])).collect())
}

/// Whether the per-call dump `kind` is wanted for call `k`: `ROSE_FROM=k0` leaves out the calls before `k0` (a window of a long run).
pub fn on_from(kind: &str, k: usize) -> bool {
    static F: OnceLock<usize> = OnceLock::new();
    on(kind)
        && k >= *F.get_or_init(|| {
            std::env::var("ROSE_FROM")
                .ok()
                .and_then(|v| v.parse().ok())
                .unwrap_or(0)
        })
}

/// Writes one record of `kind` under an explicit tag (kinds that emit several record types).
pub fn emit_tag(kind: &str, tag: &str, rec: serde_json::Value) {
    debug_assert!(
        KINDS.iter().any(|k| k.0 == kind),
        "unknown dump kind {kind}"
    );
    let line = format!("{tag} {}\n", serde_json::to_string(&rec).unwrap());
    std::io::Write::write_all(&mut std::io::stderr().lock(), line.as_bytes()).unwrap();
}
