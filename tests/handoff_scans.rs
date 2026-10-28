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

#[test]
fn terminal_handoff_child() {
    if std::env::var_os("ROSE_TERMINAL_TEST").is_none() {
        return;
    }
    let model = mps::parse(include_str!("fixtures/fx_midsing_1496.mps")).unwrap();
    let mut opts = SolveOptions::dual();
    opts.max_iterations = 200;
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
    ]);
    solve(&model, None, &opts).unwrap();
}

#[test]
fn terminal_fresh_scan_is_not_a_pricing_call() {
    let output = Command::new(std::env::current_exe().unwrap())
        .args(["--exact", "terminal_handoff_child", "--nocapture"])
        .env("ROSE_TERMINAL_TEST", "1")
        .env("ROSE_DUMP", "pricing")
        .output()
        .unwrap();
    assert!(output.status.success());
    let stderr = String::from_utf8(output.stderr).unwrap();
    let mut calls = 0;
    for line in stderr
        .lines()
        .filter_map(|line| line.strip_prefix("PLIST "))
    {
        let record: serde_json::Value = serde_json::from_str(line).unwrap();
        if record["retract"] == true {
            calls -= 1;
        } else {
            calls += 1;
        }
    }
    assert_eq!(calls, 25, "expected the reference's pricing call count");
}

#[test]
fn primal_residual_child() {
    if std::env::var_os("ROSE_RESID_TEST").is_none() {
        return;
    }
    let model = mps::parse(include_str!("fixtures/fx_midsing_1106.mps")).unwrap();
    let mut opts = SolveOptions::dual();
    opts.max_iterations = 200;
    opts.refactor_schedule = Some(vec![
        4,
        10,
        16,
        16,
        16,
        usize::MAX,
        4,
        10,
        usize::MAX,
        4,
        10,
        usize::MAX,
        4,
        10,
        usize::MAX,
        4,
        11,
        18,
        usize::MAX,
        7,
        9,
        9,
    ]);
    solve(&model, None, &opts).unwrap();
}

#[test]
fn primal_handoff_residuals_use_reference_accumulation_order() {
    let output = Command::new(std::env::current_exe().unwrap())
        .args(["--exact", "primal_residual_child", "--nocapture"])
        .env("ROSE_RESID_TEST", "1")
        .env("ROSE_DUMP", "dualtol")
        .output()
        .unwrap();
    assert!(output.status.success());
    let stderr = String::from_utf8(output.stderr).unwrap();
    let vals: Vec<f64> = stderr
        .lines()
        .filter_map(|line| line.strip_prefix("DTOL "))
        .map(|line| serde_json::from_str::<serde_json::Value>(line).unwrap())
        .filter(|r| r["src"] == "primal-cmp")
        .map(|r| r["meas"].as_f64().unwrap())
        .collect();
    assert_eq!(
        vals,
        [7.75, 9.375, 7.75, 9.375, 7.75, 9.375, 13.1875, 7.75, 9.375, 13.1875]
    );
}

#[test]
fn dual_handoff_returns_before_controller_dispatches_primal() {
    let output = Command::new(std::env::current_exe().unwrap())
        .args(["--exact", "optimal_handoff_child", "--nocapture"])
        .env("ROSE_E2_TEST", "r03125")
        .env("ROSE_DUMP", "callback")
        .output().unwrap();
    assert!(output.status.success());
    let stderr = String::from_utf8(output.stderr).unwrap();
    let events: Vec<serde_json::Value> = stderr.lines()
        .filter_map(|line| line.strip_prefix("CALLBACK "))
        .map(|line| serde_json::from_str(line).unwrap()).collect();
    let labels: Vec<_> = events.iter().map(|event| event["event"].as_str().unwrap()).collect();
    assert_eq!(labels, ["enter", "exit", "transfer", "enter", "exit", "primal_outcome"]);
    assert_eq!(events[2]["status"], -1004);
    assert_eq!(events[2]["selector"], 0);
    assert_eq!(events[3]["selector"], 0);
    assert_eq!(events[5]["error"], false);
}

#[test]
fn artificial_handoff_child() {
    if std::env::var_os("ROSE_ARTIFICIAL_TEST").is_none() {
        return;
    }
    let model = mps::parse(include_str!("fixtures/primal_5959.mps")).unwrap();
    let mut opts = SolveOptions::dual_release();
    opts.max_iterations = 1000;
    let _ = solve(&model, None, &opts);
}

#[test]
fn artificial_handoff_returns_before_primal_entry() {
    let output = Command::new(std::env::current_exe().unwrap())
        .args(["--exact", "artificial_handoff_child", "--nocapture"])
        .env("ROSE_ARTIFICIAL_TEST", "1")
        .env("ROSE_DUMP", "callback")
        .output().unwrap();
    assert!(output.status.success());
    let stderr = String::from_utf8(output.stderr).unwrap();
    let events: Vec<serde_json::Value> = stderr.lines()
        .filter_map(|line| line.strip_prefix("CALLBACK "))
        .map(|line| serde_json::from_str(line).unwrap()).collect();
    let i = events.iter().position(|event| event["route"] == "artificial").unwrap();
    assert_eq!(events[i - 1]["event"], "exit");
    assert_eq!(events[i - 1]["selector"], 1);
    assert_eq!(events[i]["status"], -1004);
    assert_eq!(events[i + 1]["event"], "enter");
    assert_eq!(events[i + 1]["selector"], 0);
    assert_eq!(events[i - 1]["carry"], events[i + 1]["carry"]);
    assert_eq!(events[i + 2]["event"], "exit");
    assert_eq!(events[i + 3]["event"], "primal_outcome");
    assert_eq!(events[i + 3]["error"], false);
}
