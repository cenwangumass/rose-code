//! `quad_lu_replay quad_lu.jsonl`: factors the basis matrices captured from the reference's quad LU (`replay/gdb_quad_lu.py`, lines `QLUIN {json}`)
//! with the pivot-order mirror in long double and prints, per factorization, the first step whose pivot row, L column or diagonal differs
//! from the reference (bit for bit; the L column as a multiset).
use rose::f80::F80;
use rose::lu_order::{factor_structure_r, set_threshold};

fn ld(v: &serde_json::Value) -> F80 {
    let h = |s: &str| u64::from_str_radix(s.trim_start_matches("0x"), 16).unwrap();
    F80::from_repr(h(v[0].as_str().unwrap()), h(v[1].as_str().unwrap()) as u16)
}

fn main() {
    let a: Vec<String> = std::env::args().collect();
    for line in std::fs::read_to_string(&a[1]).unwrap().lines() {
        let Some(js) = line.strip_prefix("QLUIN ") else {
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
        set_threshold(r["thr"].as_f64().unwrap());
        let st = factor_structure_r::<F80>(m, &cols);
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
        let lval: Vec<F80> = r["Lval"].as_array().unwrap().iter().map(ld).collect();
        let diag: Vec<F80> = r["diag"].as_array().unwrap().iter().map(ld).collect();
        let (mut fr, mut fl, mut fd) = (None, None, None);
        let steps = st.pivots.len();
        if r["ret"].as_i64() != Some(0) {
            println!("k={} singular (ret {}), {} pivots", r["k"], r["ret"], steps);
        }
        for s in 0..steps {
            if fr.is_none() && st.pivots[s].0 != p20[s] {
                fr = Some(s);
            }
            let mut x: Vec<(u64, u16)> = lval[lcp[s]..lcp[s + 1]]
                .iter()
                .map(|v| v.abs().repr())
                .collect();
            let mut y: Vec<(u64, u16)> = st.l_vals(s).iter().map(|v| v.abs().repr()).collect();
            x.sort();
            y.sort();
            if fl.is_none() && x != y {
                fl = Some(s);
                println!("  step {s}: ref {x:x?} rose {y:x?}");
            }
            if fd.is_none() && diag[s].abs().repr() != st.diag[s].abs().repr() {
                fd = Some(s);
                println!(
                    "  diag {s}: ref {:x?} rose {:x?}",
                    diag[s].repr(),
                    st.diag[s].repr()
                );
            }
        }
        println!(
            "k={} m={} first row {:?} first L {:?} first diag {:?}",
            r["k"], m, fr, fl, fd
        );
    }
}
