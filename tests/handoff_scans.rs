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

#[test]
fn optimal_handoff_child() {
    let Ok(case) = std::env::var("ROSE_E2_TEST") else {
        return;
    };
    let (text, schedule) = if case == "r03125" {
        (
            include_str!("fixtures/fx_e2_r03125.mps"),
            vec![5, 14, 22, 30],
        )
    } else {
        (
            include_str!("fixtures/fx_e2_r0625.mps"),
            vec![5, 14, 22, 29],
        )
    };
    let model = mps::parse(text).unwrap();
    let mut opts = SolveOptions::dual();
    opts.max_iterations = 100;
    opts.refactor_schedule = Some(schedule);
    let result = solve(&model, None, &opts).unwrap();
    assert_eq!(result.status, rose::Status::Optimal);
    assert_eq!(result.iterations, if case == "r03125" { 36 } else { 35 });
}

#[test]
fn optimal_primal_handoff_ends_without_an_extra_dual_pass() {
    for case in ["r03125", "r0625"] {
        let output = Command::new(std::env::current_exe().unwrap())
            .args(["--exact", "optimal_handoff_child", "--nocapture"])
            .env("ROSE_E2_TEST", case)
            .env("ROSE_DUMP", "dsx")
            .output()
            .unwrap();
        assert!(output.status.success(), "{case}: child solve failed");
        let stderr = String::from_utf8(output.stderr).unwrap();
        let calls = stderr
            .lines()
            .filter(|line| line.starts_with("DSX "))
            .count();
        assert_eq!(
            calls, 36,
            "{case}: expected the reference's dual call count"
        );
    }
}
