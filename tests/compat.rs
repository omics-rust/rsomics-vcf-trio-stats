//! Value-exact compatibility with `bcftools +trio-stats`.
//!
//! The primary assertion compares our `FLT0` rows against frozen values
//! captured from `bcftools +trio-stats 1.23.1` on the committed golden, so it
//! runs without bcftools installed. A second, version-gated test re-runs the
//! live oracle when a matching bcftools is on PATH and diffs the `FLT0` rows
//! field for field.

use std::path::PathBuf;
use std::process::Command;

use rsomics_vcf_trio_stats::{TrioSpec, trio_stats};

fn golden_vcf() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/golden/trio.vcf")
}

fn golden_ped() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/golden/trio.ped")
}

/// The `FLT0` rows produced by `bcftools +trio-stats 1.23.1` on the golden, in
/// `child father mother npass nnon_ref nmendel nnovel nuntrans ntrans nts ntv
/// tstv nhom_dnm nrecur_dnm` order.
const EXPECTED_FLT0: &[&str] = &[
    "FLT0\tCH1\tFA1\tMO1\t19\t19\t3\t2\t1\t8\t8\t10\t0.80\t1\t1",
    "FLT0\tCH2\tFA2\tMO2\t20\t8\t1\t0\t0\t6\t4\t4\t1.00\t1\t1",
];

fn our_flt0() -> Vec<String> {
    let table = trio_stats(&golden_vcf(), &TrioSpec::Ped(&golden_ped())).unwrap();
    table
        .to_text()
        .lines()
        .filter(|l| l.starts_with("FLT0"))
        .map(str::to_string)
        .collect()
}

#[test]
fn flt0_rows_match_frozen_bcftools() {
    let ours = our_flt0();
    assert_eq!(ours.len(), EXPECTED_FLT0.len(), "trio row count mismatch");
    for (got, want) in ours.iter().zip(EXPECTED_FLT0) {
        assert_eq!(got, want, "per-trio stats row differs from bcftools");
    }
}

fn bcftools_version() -> Option<String> {
    let out = Command::new("bcftools").arg("--version").output().ok()?;
    if !out.status.success() {
        return None;
    }
    let s = String::from_utf8_lossy(&out.stdout);
    let first = s.lines().next()?;
    first.split_whitespace().nth(1).map(str::to_string)
}

/// Live differential against the oracle, gated on bcftools 1.23.x being present;
/// otherwise it loudly skips so a wrong version cannot silently pass.
#[test]
fn live_differential_against_bcftools() {
    let Some(version) = bcftools_version() else {
        eprintln!("SKIP live_differential: bcftools not on PATH");
        return;
    };
    if !version.starts_with("1.23") {
        eprintln!("SKIP live_differential: bcftools {version} != 1.23.x");
        return;
    }
    let out = Command::new("bcftools")
        .arg("+trio-stats")
        .arg("-p")
        .arg(golden_ped())
        .arg(golden_vcf())
        .output()
        .expect("run bcftools +trio-stats");
    assert!(
        out.status.success(),
        "bcftools +trio-stats failed: {}",
        String::from_utf8_lossy(&out.stderr)
    );
    let oracle: Vec<String> = String::from_utf8_lossy(&out.stdout)
        .lines()
        .filter(|l| l.starts_with("FLT0"))
        .map(str::to_string)
        .collect();
    assert_eq!(our_flt0(), oracle, "live bcftools differential mismatch");
}

/// The single-trio `--pfm` path must agree with bcftools `-P`.
#[test]
fn pfm_matches_bcftools() {
    let Some(version) = bcftools_version() else {
        eprintln!("SKIP pfm_matches: bcftools not on PATH");
        return;
    };
    if !version.starts_with("1.23") {
        eprintln!("SKIP pfm_matches: bcftools {version} != 1.23.x");
        return;
    }
    let ours: Vec<String> = trio_stats(&golden_vcf(), &TrioSpec::Pfm("CH1,FA1,MO1"))
        .unwrap()
        .to_text()
        .lines()
        .filter(|l| l.starts_with("FLT0"))
        .map(str::to_string)
        .collect();
    let out = Command::new("bcftools")
        .arg("+trio-stats")
        .arg("-P")
        .arg("CH1,FA1,MO1")
        .arg(golden_vcf())
        .output()
        .expect("run bcftools +trio-stats -P");
    let oracle: Vec<String> = String::from_utf8_lossy(&out.stdout)
        .lines()
        .filter(|l| l.starts_with("FLT0"))
        .map(str::to_string)
        .collect();
    assert_eq!(ours, oracle, "pfm differential mismatch");
}
