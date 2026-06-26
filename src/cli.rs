use std::path::PathBuf;

use clap::Parser;
use rsomics_common::{CommonFlags, Result, RsomicsError, ToolMeta};
use rsomics_vcf_trio_stats::{TrioSpec, trio_stats};

pub const META: ToolMeta = ToolMeta {
    name: env!("CARGO_PKG_NAME"),
    version: env!("CARGO_PKG_VERSION"),
};

#[derive(Parser, Debug)]
#[command(
    name = "rsomics-vcf-trio-stats",
    version,
    about = "Per-trio VCF transmission and de-novo statistics — valid/non-ref GTs, Mendelian errors, novel/untransmitted/transmitted singletons, Ts/Tv, hom and recurrent DNMs (bcftools +trio-stats)"
)]
pub struct Cli {
    /// Input VCF or VCF.gz file.
    #[arg(value_name = "INPUT")]
    pub vcf: PathBuf,

    /// PED file defining the trios (6-column: family sample paternal maternal
    /// sex phenotype). A trio is a sample with both parents present in the VCF.
    #[arg(short = 'p', long, value_name = "FILE")]
    pub ped: Option<PathBuf>,

    /// A single trio given as `proband,father,mother` sample names, as an
    /// alternative to `--ped`.
    #[arg(short = 'P', long, value_name = "P,F,M", conflicts_with = "ped")]
    pub pfm: Option<String>,

    #[command(flatten)]
    pub common: CommonFlags,
}

impl Cli {
    pub fn execute(self) -> Result<()> {
        self.common.install_rayon_pool()?;
        let spec = match (self.ped.as_deref(), self.pfm.as_deref()) {
            (Some(p), None) => TrioSpec::Ped(p),
            (None, Some(s)) => TrioSpec::Pfm(s),
            (None, None) => {
                return Err(RsomicsError::InvalidInput(
                    "missing the --ped (-p) or --pfm (-P) option".into(),
                ));
            }
            (Some(_), Some(_)) => unreachable!("clap conflicts_with prevents this"),
        };
        let table = trio_stats(&self.vcf, &spec)?;

        if self.common.json {
            let env = serde_json::json!({
                "schema_version": rsomics_common::SCHEMA_VERSION,
                "tool": META.name,
                "tool_version": META.version,
                "status": "ok",
                "result": table,
            });
            println!("{}", serde_json::to_string(&env).unwrap_or_default());
        } else {
            print!("{}", table.to_text());
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use clap::CommandFactory;

    #[test]
    fn cli_debug_assert() {
        Cli::command().debug_assert();
    }
}
