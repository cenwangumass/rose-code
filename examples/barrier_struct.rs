//! `ROSE_DUMP=barrier_struct barrier_struct model.mps [dense_min]`: the barrier setup structure (AA' pattern, levels, ordering,
//! symbolic factorization) as PAT / SYMIN / ORDIN / ORD / ORDFIX / SYM records (stderr), shaped like the reference captures.
//! A model that needs an undecoded routine (the partition `0x325850`) gets PAT, SYMIN and a `BSTRUCT {"unsupported": ...}` record only.
use rose::dump::{emit_tag, on};
use rose::{mps, order};
use serde_json::json;

fn main() {
    let a: Vec<String> = std::env::args().collect();
    let model = mps::parse(&std::fs::read_to_string(&a[1]).unwrap()).unwrap();
    let dense_min = a
        .get(2)
        .map_or(order::DENSE_COL_MIN, |s| s.parse().unwrap());
    let k = "barrier_struct";
    let dense = order::dense_columns(&model, dense_min);
    let (pat, n1) = order::aat_pattern(&model, &dense);
    if on(k) {
        emit_tag(
            k,
            "PAT",
            json!({"n": pat.n, "nnz": pat.nnz(), "colptr": pat.colptr, "rowind": pat.rowind}),
        );
        emit_tag(
            k,
            "SYMIN",
            json!({"n": pat.n, "n1": n1, "rows": n1, "dense": dense.len()}),
        );
    }
    match order::structure(&model, dense_min) {
        Ok(s) => {
            if on(k) {
                emit_tag(k, "ORDIN", json!({"n": s.pattern.n, "perm": s.levels}));
                emit_tag(
                    k,
                    "ORD",
                    json!({"n": s.pattern.n, "perm": s.order.perm, "nz": s.order.nz, "ops": s.order.ops}),
                );
                emit_tag(k, "ORDFIX", json!({"nz": s.sym.exact_nz}));
                let y = &s.sym;
                emit_tag(
                    k,
                    "SYM",
                    json!({"n": s.pattern.n, "perm": y.perm, "nz": y.nz, "nidx": y.idx.len(), "ops": y.ops,
                    "width": y.width, "count": y.count, "rowoff": y.rowoff, "idx": y.idx}),
                );
            }
            println!(
                "Barrier statistics: AA' NZ {} Factor NZ {} Factor Ops {:e}",
                pat.nnz() / 2,
                s.sym.nz,
                s.sym.ops
            );
        }
        Err((_, _, _, e)) => {
            if on(k) {
                emit_tag(k, "BSTRUCT", json!({"unsupported": e.to_string()}));
            }
            println!("unsupported: {e}");
        }
    }
}
