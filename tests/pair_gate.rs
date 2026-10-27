use rose::{mps, solve, Basis, SolveOptions};
use std::process::Command;

#[test]
fn singular_pairgate_child() {
    let Ok(profile) = std::env::var("ROSE_PAIRGATE_TEST") else {
        return;
    };
    let model = mps::parse(include_str!("fixtures/fx_sing_pairgate.mps")).unwrap();
    let v: serde_json::Value =
        serde_json::from_str(include_str!("fixtures/warm_fx_sing_pairgate.basis.json")).unwrap();
    let basis = Basis {
        basic: v["basic"]
            .as_array()
            .unwrap()
            .iter()
            .map(|x| x.as_u64().unwrap() as usize)
            .collect(),
        at_upper: v["at_upper"]
            .as_array()
            .unwrap()
            .iter()
            .map(|x| x.as_bool().unwrap())
            .collect(),
    };
    let mut opts = if profile == "release" {
        SolveOptions::release()
    } else {
        SolveOptions::default()
    };
    opts.max_iterations = 100;
    solve(&model, Some(&basis), &opts).unwrap();
}

#[test]
fn singular_start_keeps_pair_attempts_and_factor_numbers() {
    for profile in ["controlled", "release"] {
        let output = Command::new(std::env::current_exe().unwrap())
            .args(["--exact", "singular_pairgate_child", "--nocapture"])
            .env("ROSE_PAIRGATE_TEST", profile)
            .env("ROSE_DUMP", "pairgate")
            .output()
            .unwrap();
        assert!(output.status.success(), "{profile}: child solve failed");
        let stderr = String::from_utf8(output.stderr).unwrap();
        let states: Vec<_> = stderr
            .lines()
            .filter_map(|l| l.strip_prefix("PGATE "))
            .map(|l| serde_json::from_str::<serde_json::Value>(l).unwrap())
            .map(|r| {
                (
                    r["rem"].as_u64().unwrap(),
                    r["att"].as_u64().unwrap(),
                    r["gate"].as_i64().unwrap(),
                    r["c2c0"].as_u64().unwrap(),
                )
            })
            .collect();
        assert_eq!(
            states,
            [
                (59, 0, 0, 1),
                (if profile == "release" { 56 } else { 58 }, 1, 0, 1),
                (56, 2, 0, 2),
                (56, 3, 0, 3)
            ],
            "{profile}"
        );
    }
}
