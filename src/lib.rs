//! Per-trio transmission and de-novo statistics from a VCF + PED, value-exact
//! with `bcftools +trio-stats`.
//!
//! The PED defines trios (a sample with both parents present in the VCF); the
//! VCF is streamed once, folding each record into every trio's counters
//! (`triostats`). The rendered table reproduces the plugin's `FLT0` output: the
//! `#`-prefixed column legend, the `DEF` line, and one `FLT0` row per trio.

mod ped;
mod triostats;
mod vcf;

use std::path::Path;

use rsomics_common::Result;

pub use ped::Trio;
pub use triostats::{Gt, Record, TrioStats, acgt2int};

/// How the trios were specified.
pub enum TrioSpec<'a> {
    /// Path to a 6-column PED file.
    Ped(&'a Path),
    /// A `proband,father,mother` sample-name triple.
    Pfm(&'a str),
}

/// One rendered trio row: the three sample names and the counters.
#[derive(Debug, Clone, serde::Serialize)]
pub struct TrioRow {
    pub child: String,
    pub father: String,
    pub mother: String,
    #[serde(flatten)]
    pub stats: TrioStats,
    pub tstv: f64,
}

/// The full per-trio table.
#[derive(Debug, Clone, serde::Serialize)]
pub struct TrioTable {
    pub rows: Vec<TrioRow>,
}

impl TrioTable {
    /// Render the table in `bcftools +trio-stats` text form: the `#` legend, the
    /// `DEF` line, then one `FLT0` row per trio. The Ts/Tv column uses the
    /// plugin's `%.2f`, printing `inf` when there are no transversions.
    #[must_use]
    pub fn to_text(&self) -> String {
        let mut out = String::new();
        out.push_str("# DEF line defines the single (unfiltered) expression\n");
        out.push_str("# FLT* lines report numbers for every trio:\n");
        out.push_str("#   1) filter id\n");
        out.push_str("#   2) child\n");
        out.push_str("#   3) father\n");
        out.push_str("#   4) mother\n");
        out.push_str("#   5) number of valid trio genotypes\n");
        out.push_str("#   6) number of non-reference trio genotypes\n");
        out.push_str("#   7) number of DNMs/Mendelian errors\n");
        out.push_str("#   8) number of novel singleton alleles in the child\n");
        out.push_str("#   9) number of untransmitted trio singletons\n");
        out.push_str("#   10) number of transmitted trio singletons\n");
        out.push_str("#   11) number of transitions\n");
        out.push_str("#   12) number of transversions\n");
        out.push_str("#   13) overall ts/tv\n");
        out.push_str("#   14) number of homozygous DNMs/Mendelian errors\n");
        out.push_str("#   15) number of recurrent DNMs/Mendelian errors\n");
        out.push_str("DEF\tFLT0\tall\n");
        for r in &self.rows {
            let s = &r.stats;
            let tstv = if s.ntv == 0 {
                "inf".to_string()
            } else {
                format!("{:.2}", r.tstv)
            };
            out.push_str(&format!(
                "FLT0\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{tstv}\t{}\t{}\n",
                r.child,
                r.father,
                r.mother,
                s.npass,
                s.nnon_ref,
                s.nmendel_err,
                s.nnovel,
                s.nsingleton,
                s.ndoubleton,
                s.nts,
                s.ntv,
                s.ndnm_hom,
                s.ndnm_recurrent,
            ));
        }
        out
    }
}

/// Compute the per-trio statistics for `vcf` given the trio specification.
pub fn trio_stats(vcf: &Path, spec: &TrioSpec) -> Result<TrioTable> {
    let names = vcf::read_sample_names(vcf)?;
    let trios = match spec {
        TrioSpec::Ped(p) => ped::parse_ped(p, &names)?,
        TrioSpec::Pfm(s) => ped::pfm_trio(s, &names)?,
    };
    let stats = vcf::read_stats(vcf, &trios)?;
    let rows = trios
        .into_iter()
        .zip(stats)
        .map(|(t, s)| {
            let tstv = s.tstv();
            TrioRow {
                child: t.child_name,
                father: t.father_name,
                mother: t.mother_name,
                stats: s,
                tstv,
            }
        })
        .collect();
    Ok(TrioTable { rows })
}
