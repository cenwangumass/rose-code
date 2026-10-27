use rose::{mps, solve, SolveOptions};
use std::process::Command;

#[test]
fn restart_threshold_child() {
    if std::env::var_os("ROSE_RESTART_TEST").is_none() {
        return;
    }
    let model = mps::parse(include_str!("fixtures/randlp50154.mps")).unwrap();
    let mut opts = SolveOptions::dual();
    opts.max_iterations = 25;
    opts.refactor_schedule = Some(vec![
        4,
        11,
        usize::MAX,
        4,
        11,
        18,
        usize::MAX,
        4,
        11,
        18,
        usize::MAX,
        4,
        11,
        usize::MAX,
        4,
        11,
        18,
        usize::MAX,
        4,
        11,
        18,
    ]);
    let _ = solve(&model, None, &opts);
}

#[test]
fn singular_restart_resets_large_value_pricing_gate() {
    let output = Command::new(std::env::current_exe().unwrap())
        .args(["--exact", "restart_threshold_child", "--nocapture"])
        .env("ROSE_RESTART_TEST", "1")
        .env("ROSE_DUMP", "dsx")
        .output()
        .unwrap();
    assert!(output.status.success(), "bounded child solve failed");
    let stderr = String::from_utf8(output.stderr).unwrap();
    let calls: Vec<serde_json::Value> = stderr
        .lines()
        .filter_map(|l| l.strip_prefix("DSX "))
        .map(|l| serde_json::from_str(l).unwrap())
        .collect();
    let event = &calls[46];
    assert_eq!(event["thr"].as_f64(), Some(10.0));
    assert_eq!(event["list"], serde_json::json!([1, 3, 5, 8, 9, 11]));
}

#[test]
fn repaired_confirmation_child() {
    if std::env::var_os("ROSE_REPAIR_CONFIRM_TEST").is_none() {
        return;
    }
    let model = mps::parse(include_str!("fixtures/randlp50154.mps")).unwrap();
    let mut opts = SolveOptions::dual();
    opts.max_iterations = 150;
    opts.refactor_schedule = Some(vec![
        4,
        11,
        usize::MAX,
        4,
        11,
        18,
        usize::MAX,
        4,
        11,
        18,
        usize::MAX,
        4,
        11,
        usize::MAX,
        4,
        11,
        18,
        usize::MAX,
        4,
        11,
        18,
    ]);
    let _ = solve(&model, None, &opts);
}

#[test]
fn repaired_confirmation_hands_off_without_extra_dual_pricing() {
    let output = Command::new(std::env::current_exe().unwrap())
        .args(["--exact", "repaired_confirmation_child", "--nocapture"])
        .env("ROSE_REPAIR_CONFIRM_TEST", "1")
        .env("ROSE_DUMP", "dsx,lustart")
        .output()
        .unwrap();
    assert!(output.status.success());
    let stderr = String::from_utf8(output.stderr).unwrap();
    assert_eq!(
        stderr.lines().filter(|l| l.starts_with("DSX ")).count(),
        142
    );
    assert!(stderr
        .lines()
        .filter_map(|l| l.strip_prefix("LUSTART "))
        .map(|l| serde_json::from_str::<serde_json::Value>(l).unwrap())
        .any(|v| v["primal"] == true));
}
