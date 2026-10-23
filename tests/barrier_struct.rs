//! Barrier setup structure (`rose::order`) against the reference captures (`tests/fixtures/barrier-struct-controlled.json`:
//! AFIRO, BLEND, BORE3D; exact permutations, counts and row lists of the ordering and the symbolic factorization).
use rose::order::{self, Pattern};
use serde_json::Value;

fn ints(v: &Value) -> Vec<i64> {
    v.as_array().unwrap().iter().map(|x| x.as_i64().unwrap()).collect()
}

#[test]
fn structure_matches_reference_captures() {
    let fx: Value = serde_json::from_str(include_str!("fixtures/barrier-struct-controlled.json")).unwrap();
    for (name, mps_text) in [
        ("afiro", include_str!("fixtures/afiro.mps")),
        ("blend", include_str!("fixtures/blend.mps")),
        ("bore3d", include_str!("fixtures/bore3d.mps")),
    ] {
        let want = &fx[name];
        let model = rose::mps::parse(mps_text).unwrap();
        let s = order::structure(&model, order::DENSE_COL_MIN).unwrap_or_else(|e| panic!("{name}: {}", e.3));
        assert_eq!(s.pattern.nnz() as i64, want["aat_nnz"].as_i64().unwrap(), "{name} AA' nnz");
        assert_eq!(s.order.perm, ints(&want["ord_perm"]), "{name} ORD perm");
        assert_eq!(s.order.nz, want["ord_nz"].as_i64().unwrap(), "{name} ORD nz");
        assert_eq!(s.order.ops, want["ord_ops"].as_f64().unwrap(), "{name} ORD ops");
        assert_eq!(s.sym.exact_nz, want["exact_nz"].as_i64().unwrap(), "{name} exact nz");
        assert_eq!(s.sym.perm, ints(&want["sym_perm"]), "{name} SYM perm");
        assert_eq!(s.sym.nz, want["sym_nz"].as_i64().unwrap(), "{name} Factor NZ");
        assert_eq!(s.sym.ops, want["sym_ops"].as_f64().unwrap(), "{name} Factor Ops");
        assert_eq!(s.sym.width, ints(&want["width"]), "{name} widths");
        assert_eq!(s.sym.count, ints(&want["count"]), "{name} counts");
        assert_eq!(s.sym.idx, ints(&want["idx"]), "{name} row lists");
    }
}

#[test]
fn afiro_log_numbers() {
    let model = rose::mps::parse(include_str!("fixtures/afiro.mps")).unwrap();
    let s = order::structure(&model, order::DENSE_COL_MIN).unwrap();
    // the reference log: AA' NZ 63, Factor NZ 218, Factor Ops 2.146e+03
    assert_eq!((s.pattern.nnz() / 2, s.sym.nz, s.sym.ops), (63, 218, 2146.0));
}

/// A connected component larger than max(600, n/100) goes through the partition 0x325850 (a path of 700 nodes: the separator it finds is a set of
/// level -2 nodes; labels supplied by the caller replace it).
#[test]
fn large_component_is_partitioned() {
    let n = 700;
    let mut colptr = vec![0];
    let mut rowind = Vec::new();
    for i in 0..n {
        if i > 0 {
            rowind.push(i - 1);
        }
        if i + 1 < n {
            rowind.push(i + 1);
        }
        colptr.push(rowind.len());
    }
    let p = Pattern { n, colptr, rowind };
    let lv = order::levels(&p, n, None).unwrap();
    assert!(lv.iter().any(|&x| x == -2) && lv.iter().any(|&x| x == -1));
    // with labels supplied by the caller, the separator nodes (label 2) move to level -2
    let mut b = |sz: usize, _: &[usize], _: &[usize]| (0..sz).map(|i| if i == 350 { 2 } else { 0 }).collect::<Vec<i64>>();
    let lv = order::levels(&p, n, Some(&mut b)).unwrap();
    assert_eq!(lv.iter().filter(|&&x| x == -2).count(), 1);
    assert_eq!(lv[350], -2);
}

/// The partition 0x325850 of PEROLD's 625-node component (captured by `replay/gdb_barrier_part.py`, generator state 0): the labels agree.
#[test]
fn partition_matches_reference_labels() {
    let fx: Value = serde_json::from_str(include_str!("fixtures/partition-perold.json")).unwrap();
    let colptr: Vec<usize> = ints(&fx["colptr"]).iter().map(|&x| x as usize).collect();
    let rowind: Vec<usize> = ints(&fx["rowind"]).iter().map(|&x| x as usize).collect();
    let mut rng = 0u64;
    let lab = order::partition(fx["size"].as_u64().unwrap() as usize, &colptr, &rowind, &mut rng);
    assert_eq!(lab, ints(&fx["labels"]));
}
