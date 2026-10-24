//! `lu_cols luin.out k`: prints the mirror's pivot (row, position) sequence and phases for factorization k.
use rose::lu_order::factor_structure;

fn main() {
    let a: Vec<String> = std::env::args().collect();
    let k: u64 = a[2].parse().unwrap();
    for line in std::fs::read_to_string(&a[1]).unwrap().lines() {
        let Some(js) = line.strip_prefix("LUIN ") else {
            continue;
        };
        let r: serde_json::Value = serde_json::from_str(js).unwrap();
        if r["k"].as_u64().unwrap() != k {
            continue;
        }
        let m = r["m"].as_u64().unwrap() as usize;
        let n = r["n"].as_u64().unwrap();
        let mut cols = Vec::new();
        for c in r["mat"].as_array().unwrap() {
            let v = c[0].as_u64().unwrap();
            let rows: Vec<usize> = c[1]
                .as_array()
                .unwrap()
                .iter()
                .map(|x| x.as_u64().unwrap() as usize)
                .collect();
            let vals: Vec<f64> = c[2]
                .as_array()
                .unwrap()
                .iter()
                .map(|x| x.as_f64().unwrap())
                .collect();
            cols.push((v >= n, rows, vals));
        }
        let st = factor_structure(m, &cols);
        println!("{}", serde_json::to_string(&(st.pivots, st.phase)).unwrap());
    }
}
