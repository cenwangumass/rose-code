//! Pins the system libm results of the transcendental functions rose calls (`powf`, `ln`, `exp`, `log2`, `log10`).
//!
//! These go to the platform's libm (glibc here), whose results are not specified to the last bit and change between libm versions; rose's
//! agreement with the reference was established with the values pinned below (glibc 2.44, x86-64). A failure means this libm gives other bits
//! somewhere on these inputs: rerun the parity gate (`rose-parity check --tier full`) before trusting results built on it, and re-pin only
//! when it passes.
//!
//! Each case hashes the result bits over a fixed deterministic input set covering the ranges the call sites see (integer sizes for the
//! `m^p` limits, wide positive ranges for the numeric ones).

fn mix(h: u64, x: f64) -> u64 {
    (h ^ x.to_bits())
        .wrapping_mul(0x9e37_79b9_7f4a_7c15)
        .rotate_left(29)
}

/// Deterministic doubles spread over `[2^lo, 2^hi)`: random mantissas and exponents (xorshift).
fn spread(lo: i32, hi: i32, n: usize) -> Vec<f64> {
    let mut s = 0x2545_f491_4f6c_dd1du64;
    (0..n)
        .map(|_| {
            s ^= s << 13;
            s ^= s >> 7;
            s ^= s << 17;
            let e = lo + (s % (hi - lo) as u64) as i32;
            let m = f64::from_bits(0x3ff0_0000_0000_0000 | (s >> 12));
            m * rose::real::pow2(e)
        })
        .collect()
}

fn ints(n: u32) -> Vec<f64> {
    (1..=n).map(f64::from).collect()
}

fn pin(name: &str, xs: &[f64], f: impl Fn(f64) -> f64) -> u64 {
    let h = xs.iter().fold(0u64, |h, &x| mix(h, f(x)));
    println!("{name}: {h:#018x}");
    h
}

/// (name, hash) of every case; `cargo test --release --test libm_pins -- --nocapture` prints the current values.
fn cases() -> Vec<(&'static str, u64)> {
    let sizes = ints(200_000);
    let wide = spread(-60, 60, 200_000);
    let small = spread(-40, 7, 200_000);
    vec![
        (
            "powf 0.7 sizes",
            pin("powf 0.7 sizes", &sizes, |x| x.powf(0.7)),
        ),
        (
            "powf 0.65 sizes",
            pin("powf 0.65 sizes", &sizes, |x| x.powf(0.65)),
        ),
        (
            "powf 0.45 sizes",
            pin("powf 0.45 sizes", &sizes, |x| x.powf(0.45)),
        ),
        (
            "powf 0.35 wide",
            pin("powf 0.35 wide", &wide, |x| x.powf(0.35)),
        ),
        (
            "powf 0.25 wide",
            pin("powf 0.25 wide", &wide, |x| x.powf(0.25)),
        ),
        ("log2 wide", pin("log2 wide", &wide, f64::log2)),
        (
            "log2 sizes",
            pin("log2 sizes", &sizes, |x| (x + 1.0).log2()),
        ),
        ("log10 sizes", pin("log10 sizes", &sizes, f64::log10)),
        ("ln wide", pin("ln wide", &wide, f64::ln)),
        ("exp small", pin("exp small", &small, |x| (x - 64.0).exp())),
        ("exp signed", pin("exp signed", &small, |x| (-x).exp())),
    ]
}

#[test]
fn libm_results_are_the_pinned_ones() {
    // glibc 2.44, x86-64
    let pinned: &[(&str, u64)] = &[
        ("powf 0.7 sizes", 0xc456dd4979d933a4),
        ("powf 0.65 sizes", 0x52649c496e363e10),
        ("powf 0.45 sizes", 0xdbc32911048b2687),
        ("powf 0.35 wide", 0xe4b2d4dfa601bbd9),
        ("powf 0.25 wide", 0x5b6cb710fb59fda2),
        ("log2 wide", 0x5d8d22b86001222a),
        ("log2 sizes", 0x047443fecc5d664b),
        ("log10 sizes", 0x21d026b9f2b97ddb),
        ("ln wide", 0x32eff004a9f53607),
        ("exp small", 0x959a2e381cfae9ea),
        ("exp signed", 0xb4623b7c58619769),
    ];
    let now = cases();
    let mut bad = Vec::new();
    for (name, h) in &now {
        match pinned.iter().find(|p| p.0 == *name) {
            Some(p) if p.1 == *h => {}
            Some(p) => bad.push(format!("{name}: pinned {:#018x}, this libm {h:#018x}", p.1)),
            None => bad.push(format!("{name}: not pinned (this libm {h:#018x})")),
        }
    }
    assert!(
        bad.is_empty(),
        "libm results differ from the pinned ones:\n{}",
        bad.join("\n")
    );
}
