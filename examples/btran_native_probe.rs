//! Isolated dense extended-BTRAN L-transpose feasibility screen.
//! Run with `cargo run --release --example btran_native_probe`.
//! This is not a solver implementation: inputs are finite, the control word must
//! have masked exceptions/nearest rounding/64-bit precision, and exceptional
//! intermediate results are not supported. Integration requires a software fallback.
use rose::f80::F80;
use std::{arch::asm, hint::black_box, time::Instant};

type Packed = [u8; 16];

fn pack(v: F80) -> Packed {
    let (m, s) = v.repr();
    let mut b = [0; 16];
    b[..8].copy_from_slice(&m.to_le_bytes());
    b[8..10].copy_from_slice(&s.to_le_bytes());
    b
}

fn unpack(b: Packed) -> F80 {
    F80::from_repr(
        u64::from_le_bytes(b[..8].try_into().unwrap()),
        u16::from_le_bytes(b[8..10].try_into().unwrap()),
    )
}

fn software(mut x: Vec<F80>, lt: &[Vec<(usize, F80)>]) -> Vec<F80> {
    for i in (0..x.len()).rev() {
        let v = x[i];
        if v.abs() > F80::from(1e-13) {
            for &(j, a) in &lt[i] {
                x[j] -= a * v;
            }
        }
    }
    x
}

fn native(x: Vec<F80>, lt: &[Vec<(usize, F80)>]) -> Vec<F80> {
    let mut x: Vec<_> = x.into_iter().map(pack).collect();
    let tol = pack(F80::from(1e-13));
    for i in (0..x.len()).rev() {
        // The prototype keeps the software threshold predicate to isolate arithmetic
        // and includes every representation conversion in its timed work.
        if !(unpack(x[i]).abs() > unpack(tol)) {
            continue;
        }
        let v = x[i];
        for &(j, a) in &lt[i] {
            let a = pack(a);
            let mut product = [0u8; 16];
            // SAFETY: each input/output owns 16 bytes for a ten-byte load/store.
            // The product store keeps multiply and subtract roundings separate;
            // the x87 stack is balanced and both used registers are clobbered.
            unsafe {
                asm!("fld tbyte ptr [{a}]", "fld tbyte ptr [{v}]", "fmulp st(1), st",
                    "fstp tbyte ptr [{p}]", "fld tbyte ptr [{x}]", "fld tbyte ptr [{p}]",
                    "fsubp st(1), st", "fstp tbyte ptr [{x}]",
                    a=in(reg) a.as_ptr(),v=in(reg) v.as_ptr(),p=in(reg) product.as_mut_ptr(),
                    x=in(reg) x[j].as_mut_ptr(),out("st(0)") _,out("st(1)") _,options(nostack,preserves_flags));
            }
        }
    }
    x.into_iter().map(unpack).collect()
}

fn main() {
    let mut cw = 0u16;
    unsafe {
        asm!("fnstcw [{cw}]",cw=in(reg)&mut cw,options(nostack,preserves_flags));
    }
    assert_eq!(cw & 0x0f3f, 0x033f);
    let captured: serde_json::Value =
        serde_json::from_str(include_str!("../tests/fixtures/ftran_l_captured.json")).unwrap();
    let mut cases = Vec::new();
    let mut seed = 0x985abc125u64;
    let mut next = || {
        seed ^= seed << 13;
        seed ^= seed >> 7;
        seed ^= seed << 17;
        seed
    };
    for rec in captured.as_array().unwrap() {
        let rows: Vec<Vec<usize>> = serde_json::from_value(rec["lcols"].clone()).unwrap();
        let mut lt = vec![Vec::new(); rows.len()];
        for (i, rs) in rows.iter().enumerate() {
            for (e, &j) in rs.iter().enumerate() {
                lt[j].push((i, F80::from(rec["lvals"][i][e].as_f64().unwrap())));
            }
        }
        let x = vec![F80::from(1.25); rows.len()];
        cases.push((x, lt));
    }
    for n in [0, 1, 8, 32, 64] {
        for round in 0..32 {
            let mut val = || {
                let k = next();
                if round == 0 {
                    F80::from(if k & 1 == 0 { 0.0 } else { -0.0 })
                } else {
                    F80::from_repr(
                        (1 << 63) | (k >> 1),
                        (16370 + (k % 26)) as u16 | ((k as u16 & 1) << 15),
                    )
                }
            };
            let x = (0..n).map(|_| val()).collect();
            let lt = (0..n)
                .map(|i| {
                    (0..i)
                        .filter(|j| (i + j) % 3 == 0)
                        .map(|j| (j, val()))
                        .collect()
                })
                .collect();
            cases.push((x, lt));
        }
    }
    for (mant, exp) in [(1, 0), (u64::MAX >> 1, 0), (1 << 63, 1)] {
        for sign in [0, 0x8000] {
            cases.push((
                vec![F80::from(-0.0), F80::from(1.0)],
                vec![vec![], vec![(0, F80::from_repr(mant, exp | sign))]],
            ));
        }
    }
    for (x, lt) in &cases {
        let old = software(x.clone(), lt);
        let new = native(x.clone(), lt);
        assert_eq!(
            old.iter().map(|v| v.repr()).collect::<Vec<_>>(),
            new.iter().map(|v| v.repr()).collect::<Vec<_>>()
        );
    }
    println!("differential cases {}", cases.len());
    for pair in 0..4 {
        for name in if pair % 2 == 0 {
            ["software", "native"]
        } else {
            ["native", "software"]
        } {
            let start = Instant::now();
            for _ in 0..100 {
                for (x, lt) in &cases {
                    let x = black_box(x.clone());
                    black_box(if name == "native" {
                        native(x, black_box(lt))
                    } else {
                        software(x, black_box(lt))
                    });
                }
            }
            println!("pair {pair} {name} {}", start.elapsed().as_secs_f64());
        }
    }
}
