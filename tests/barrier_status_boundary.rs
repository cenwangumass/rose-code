//! Reference 0x2e2fe0 inputs from fx_bar_both3rhs, controlled iteration 17.
//! GDB changed only Q+0x148 to -Q+0x138 or its adjacent binary64 values;
//! each recorded status is the reference return before the solve resumes.
use rose::barrier::{status, TermIn};
use serde_json::Value;

#[test]
fn combined_sign_boundary_from_reference_state() {
    let cases: Vec<Value> =
        serde_json::from_str(include_str!("fixtures/barrier-status-boundary.json")).unwrap();
    for c in cases {
        let q = &c["q"];
        let p = &c["P"];
        let word = |map: &Value, off: &str| map[off].as_u64().unwrap();
        let f = |off: &str| f64::from_bits(word(q, off));
        let i = |off: &str| word(q, off) as u32 as i32;
        let h = |off: &str| {
            c["hist"][off]
                .as_array()
                .unwrap()
                .iter()
                .map(|v| v.as_f64().unwrap())
                .collect()
        };
        let t = TermIn {
            it: c["it"].as_i64().unwrap() as i32,
            q130: f("0x130"),
            q138: f("0x138"),
            q140: f("0x140"),
            q148: f("0x148"),
            q128: f("0x128"),
            q1e8: f("0x1e8"),
            off: f64::from_bits(c["m60"].as_u64().unwrap()),
            q1f0: f("0x1f0"),
            q1f8: f("0x1f8"),
            tol: f("0x258"),
            ncomp: word(q, "0x210") as f64,
            tau: f("0x690"),
            kappa: f("0x698"),
            q198: f("0x198"),
            q1a8: f("0x1a8"),
            q1d8: f("0x1d8"),
            q1e0: f("0x1e0"),
            q1a0: f("0x1a0"),
            q1b0: f("0x1b0"),
            q1c0: f("0x1c0"),
            q1d0: f("0x1d0"),
            q220: f("0x220"),
            q228: f("0x228"),
            q160: f("0x160"),
            q170: f("0x170"),
            q260: f("0x260"),
            q268: f("0x268"),
            q270: f("0x270"),
            q278: f("0x278"),
            q2a0: i("0x2a0"),
            q5e0: f("0x5e0"),
            q5e8: f("0x5e8"),
            i88: i("0x88"),
            if8: i("0xf8"),
            ia8: i("0xa8"),
            if4: (word(q, "0xf0") >> 32) as u32 as i32,
            p2dd4: word(p, "0x2dd4") as u32 != 0,
            p2460: f64::from_bits(word(p, "0x2460")),
            p2bbc: word(p, "0x2bbc") as u32 as i32,
            p2c10: word(p, "0x2c10") as u32 as i32,
            time_up: false,
            h50: h("0x50"),
            h68: h("0x68"),
            h78: h("0x78"),
            h80: h("0x80"),
        };
        let expected = c["ret"].as_i64().unwrap() as i32;
        assert_eq!(status(&t), (expected, 6), "delta {}", c["delta"]);
    }
}
