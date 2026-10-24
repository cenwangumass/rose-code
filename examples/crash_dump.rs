//! `crash_dump model.mps`: prints the crash basis as the reference's status array (basis position, -1 lower, -2 upper) on one line.
use rose::{crash_basis, mps};
fn main() {
    let a: Vec<String> = std::env::args().collect();
    let model = mps::parse(&std::fs::read_to_string(&a[1]).unwrap()).unwrap();
    let b = crash_basis(&model);
    let mut st: Vec<i64> = b
        .at_upper
        .iter()
        .map(|&u| if u { -2 } else { -1 })
        .collect();
    for (p, &j) in b.basic.iter().enumerate() {
        st[j] = p as i64;
    }
    println!("{}", serde_json::to_string(&st).unwrap());
}
