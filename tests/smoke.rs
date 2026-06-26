//! End-to-end smoke tests driving the built binary on the golden fixture.

use std::path::PathBuf;
use std::process::Command;

fn bin() -> PathBuf {
    PathBuf::from(env!("CARGO_BIN_EXE_rsomics-vcf-trio-stats"))
}

fn golden_vcf() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/golden/trio.vcf")
}

fn golden_ped() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/golden/trio.ped")
}

#[test]
fn ped_run_emits_two_trio_rows() {
    let out = Command::new(bin())
        .arg(golden_vcf())
        .arg("--ped")
        .arg(golden_ped())
        .output()
        .unwrap();
    assert!(out.status.success());
    let text = String::from_utf8(out.stdout).unwrap();
    let rows: Vec<&str> = text.lines().filter(|l| l.starts_with("FLT0")).collect();
    assert_eq!(rows.len(), 2);
    assert!(text.contains("DEF\tFLT0\tall"));
}

#[test]
fn json_envelope() {
    let out = Command::new(bin())
        .arg(golden_vcf())
        .arg("--ped")
        .arg(golden_ped())
        .arg("--json")
        .output()
        .unwrap();
    assert!(out.status.success());
    let v: serde_json::Value = serde_json::from_slice(&out.stdout).unwrap();
    assert_eq!(v["status"], "ok");
    assert_eq!(v["result"]["rows"].as_array().unwrap().len(), 2);
    assert_eq!(v["result"]["rows"][0]["child"], "CH1");
    assert_eq!(v["result"]["rows"][0]["npass"], 19);
}

#[test]
fn pfm_run() {
    let out = Command::new(bin())
        .arg(golden_vcf())
        .arg("--pfm")
        .arg("CH1,FA1,MO1")
        .output()
        .unwrap();
    assert!(out.status.success());
    let text = String::from_utf8(out.stdout).unwrap();
    assert_eq!(text.lines().filter(|l| l.starts_with("FLT0")).count(), 1);
}

#[test]
fn missing_trio_spec_fails() {
    let out = Command::new(bin()).arg(golden_vcf()).output().unwrap();
    assert!(!out.status.success());
    let err = String::from_utf8(out.stderr).unwrap();
    assert!(err.contains("--ped") || err.contains("--pfm"));
}

#[test]
fn unknown_sample_in_pfm_fails() {
    let out = Command::new(bin())
        .arg(golden_vcf())
        .arg("--pfm")
        .arg("NOPE,FA1,MO1")
        .output()
        .unwrap();
    assert!(!out.status.success());
}

#[test]
fn help_exits_zero() {
    let out = Command::new(bin()).arg("--help").output().unwrap();
    assert!(out.status.success());
}
