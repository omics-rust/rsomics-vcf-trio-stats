//! End-to-end smoke tests driving the built binary on the golden fixture.

use std::io::Write;
use std::path::PathBuf;
use std::process::Command;

use serde::Deserialize;

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
    // A single JSON document: Deserializer::end() errors on any trailing bytes,
    // so a double-printed envelope would fail here.
    let mut de = serde_json::Deserializer::from_slice(&out.stdout);
    let v = serde_json::Value::deserialize(&mut de).unwrap();
    de.end().unwrap();
    assert_eq!(v["status"], "ok");
    let rows = v["result"]["rows"].as_array().unwrap();
    assert_eq!(rows.len(), 2);
    assert_eq!(rows[0]["child"], "CH1");
    assert_eq!(rows[0]["npass"], 19);
    assert_eq!(rows[0]["ndnm_recurrent"], 1);
}

/// A record whose FORMAT lacks a GT subfield must fail loudly: non-zero exit and
/// a diagnostic on stderr, never a silent wrong result.
#[test]
fn vcf_without_gt_fails_loud() {
    let dir = tempfile::tempdir().unwrap();
    let vcf = dir.path().join("nogt.vcf");
    let mut f = std::fs::File::create(&vcf).unwrap();
    write!(
        f,
        "##fileformat=VCFv4.2\n\
         ##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Depth\">\n\
         #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tCH1\tFA1\tMO1\n\
         chr1\t100\t.\tA\tG\t30\tPASS\t.\tDP\t50\t50\t50\n"
    )
    .unwrap();
    let out = Command::new(bin())
        .arg(&vcf)
        .arg("--pfm")
        .arg("CH1,FA1,MO1")
        .output()
        .unwrap();
    assert!(!out.status.success());
    assert!(!out.stderr.is_empty(), "expected a diagnostic on stderr");
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
