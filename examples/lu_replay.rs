//! `lu_replay luin.out`: factors the basis matrices captured from the reference (`reverse-engineer/replay/gdb_lu_in.py` lines `LUIN {json}`)
//! with the pivot-order mirror and prints, per factorization, the first step whose pivot row, L column (as a set of absolute values) or
//! diagonal differs from the reference.
use rose::lu_order::factor_structure;

fn main() {
    let a: Vec<String> = std::env::args().collect();
    for line in std::fs::read_to_string(&a[1]).unwrap().lines() {
        let Some(js) = line.strip_prefix("LUIN ") else {
            continue;
        };
        let r: serde_json::Value = serde_json::from_str(js).unwrap();
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
        let p20: Vec<usize> = r["p20"]
            .as_array()
            .unwrap()
            .iter()
            .map(|x| x.as_u64().unwrap() as usize)
            .collect();
        let lcp: Vec<usize> = r["Lcp"]
            .as_array()
            .unwrap()
            .iter()
            .map(|x| x.as_u64().unwrap() as usize)
            .collect();
        let lval: Vec<f64> = r["Lval"]
            .as_array()
            .unwrap()
            .iter()
            .map(|x| x.as_f64().unwrap())
            .collect();
        let diag: Vec<f64> = r["diag"]
            .as_array()
            .unwrap()
            .iter()
            .map(|x| x.as_f64().unwrap())
            .collect();
        let mut first_row = None;
        let mut first_l = None;
        let mut first_d = None;
        for s in 0..m {
            if first_row.is_none() && st.pivots[s].0 != p20[s] {
                first_row = Some(s);
            }
            let mut a: Vec<u64> = lval[lcp[s]..lcp[s + 1]]
                .iter()
                .map(|v| v.abs().to_bits())
                .collect();
            let mut b: Vec<u64> = st.l_vals(s).iter().map(|v| v.abs().to_bits()).collect();
            a.sort();
            b.sort();
            if first_l.is_none() && a != b {
                first_l = Some(s);
            }
            if first_d.is_none() && diag[s].abs() != st.diag[s].abs() {
                first_d = Some(s);
            }
        }
        if let Some(s) = first_l {
            println!(
                "  step {s}: ref {:?} rose {:?}",
                &lval[lcp[s]..lcp[s + 1]],
                st.l_vals(s)
            );
        }
        println!(
            "k={} m={} nucleus={} first row {:?} first L {:?} first diag {:?}",
            r["k"], m, st.nuc_nnz, first_row, first_l, first_d
        );
    }
}
