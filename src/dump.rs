//! Diagnostic dumps for comparison with the reference captures, selected with `ROSE_DUMP=kind,kind` (or `all`).
//!
//! Each event is one line on stderr, `TAG {json}`, with floats written to round-trip exactly. The tags and record fields
//! follow the reference captures of `reverse-engineer/tools/refcap.py`; `reverse-engineer/tools/compare.py` runs rose and
//! compares the two. Add a kind here rather than an ad hoc `eprintln!` behind its own variable.
use std::sync::OnceLock;

/// Dump kinds, their line tags and what they record.
pub const KINDS: &[(&str, &str, &str)] = &[
    ("header", "PHDR", "once per solve: n, m and the fixed variables (lo == up)"),
    ("pricing", "PLIST", "each pricing call: call index, candidate list, d, w, basic set and x"),
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
    ("fillparts", "FILLP", "per factor update: pool fill on entry, eta size, new column entries added, and the row growths (row, cnt, cap, extra)"),
];

fn kinds() -> &'static [String] {
    static K: OnceLock<Vec<String>> = OnceLock::new();
    K.get_or_init(|| {
        std::env::var("ROSE_DUMP")
            .map(|v| v.split(',').map(|s| s.trim().to_string()).filter(|s| !s.is_empty()).collect())
            .unwrap_or_default()
    })
}

/// Whether the dump `kind` is selected.
pub fn on(kind: &str) -> bool {
    debug_assert!(KINDS.iter().any(|k| k.0 == kind), "unknown dump kind {kind}");
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
