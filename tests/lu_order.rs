//! The pivot order of the reference LU factorization (singleton phases and Markowitz nucleus) on 26 factorizations captured from
//! the reference (ADLITTLE, BLEND, SC50A, BORE3D, controlled profile): every pivot-row sequence must be identical.

use rose::lu_order::pivot_order;

#[test]
fn pivot_rows_match_the_reference() {
    let text = std::fs::read_to_string(format!("{}/tests/fixtures/lu-pivot-orders.json", env!("CARGO_MANIFEST_DIR"))).unwrap();
    let data: serde_json::Value = serde_json::from_str(&text).unwrap();
    let mut total = 0;
    for (model, cases) in data.as_object().unwrap() {
        for (k, case) in cases.as_array().unwrap().iter().enumerate() {
            let m = case["m"].as_u64().unwrap() as usize;
            let n = case["n"].as_u64().unwrap() as usize;
            let cols: Vec<(bool, Vec<usize>, Vec<f64>)> = case["mat"]
                .as_array()
                .unwrap()
                .iter()
                .map(|c| {
                    let c = c.as_array().unwrap();
                    let v = c[0].as_u64().unwrap() as usize;
                    let rows = c[1].as_array().unwrap().iter().map(|x| x.as_u64().unwrap() as usize).collect();
                    let vals = c[2].as_array().unwrap().iter().map(|x| x.as_f64().unwrap()).collect();
                    (v >= n, rows, vals)
                })
                .collect();
            let want: Vec<usize> = case["p20"].as_array().unwrap().iter().map(|x| x.as_u64().unwrap() as usize).collect();
            let got: Vec<usize> = pivot_order(m, &cols).into_iter().map(|(r, _)| r).collect();
            assert_eq!(got, want, "{model} factorization {k}");
            total += 1;
        }
    }
    assert_eq!(total, 26);
}
