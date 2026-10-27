//! Reference counts for the primal scans made during a dual handoff.
use rose::{mps, solve, SolveOptions};
use std::process::Command;

#[test]
fn primal_handoff_child() {
    if std::env::var_os("ROSE_HANDOFF_TEST").is_none() {
        return;
    }
    let model = mps::parse(include_str!("fixtures/fx_midsing_1031.mps")).unwrap();
    let mut opts = SolveOptions::dual();
    opts.max_iterations = 100;
    opts.refactor_schedule = Some(vec![
        3,
        usize::MAX,
        3,
        usize::MAX,
        3,
        usize::MAX,
        3,
        usize::MAX,
        3,
        usize::MAX,
        3,
        3,
    ]);
    solve(&model, None, &opts).unwrap();
}

#[test]
fn failed_scan_refactors_before_primal_ratio() {
    let output = Command::new(std::env::current_exe().unwrap())
        .args(["--exact", "primal_handoff_child", "--nocapture"])
        .env("ROSE_HANDOFF_TEST", "1")
        .env("ROSE_DUMP", "ent")
        .output()
        .unwrap();
    assert!(output.status.success());
    let stderr = String::from_utf8(output.stderr).unwrap();
    let calls = stderr
        .lines()
        .filter(|line| line.starts_with("ENT "))
        .count();
    assert_eq!(
        calls, 3,
        "expected the three reference phase-one ratio calls"
    );
}
