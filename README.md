# rsomics-vcf-trio-stats

Per-trio transmission and de-novo statistics from a VCF and a PED, value-exact
with `bcftools +trio-stats`.

A PED file (or an explicit `proband,father,mother` triple) defines the trios; a
trio is a sample with both its parents present in the VCF. The VCF is streamed
once and each record is folded into every trio's counters, producing one row per
trio with the same columns `bcftools +trio-stats` emits.

## Usage

```sh
rsomics-vcf-trio-stats in.vcf[.gz] --ped family.ped
rsomics-vcf-trio-stats in.vcf[.gz] --pfm child,father,mother
```

`-t/--threads`, `-q/--quiet` and `--json` are accepted via the shared flags.

## Output columns

Each `FLT0` row carries the child, father and mother sample names followed by:

| # | Column | Meaning |
|---|--------|---------|
| 5 | nGood | valid trio genotypes (all three members non-missing) |
| 6 | nNonRef | trio genotypes with at least one alternate allele |
| 7 | nMendelErr | DNMs / Mendelian errors (a child allele derivable from neither parent) |
| 8 | nNovel | novel singleton allele present only in the child |
| 9 | nUntrans | untransmitted trio singletons (one alt allele in one parent only) |
| 10 | nTrans | transmitted trio singletons (one alt allele in a parent and the child) |
| 11 | nTs | transitions, over the distinct ALT alleles in the trio |
| 12 | nTv | transversions, over the distinct ALT alleles in the trio |
| 13 | Ts/Tv | nTs / nTv (`inf` when nTv is 0) |
| 14 | nHomDNM | homozygous DNMs |
| 15 | nRecurDNM | recurrent DNMs (the non-inherited allele also seen in other samples) |

### Classification

- **Mendelian error / DNM** — for the child's ordered alleles `(c0, c1)`, the
  trio is consistent if `c0` is one of the father's alleles and `c1` is one of
  the mother's, or vice versa. When neither orientation holds, the site is a
  Mendelian error. A homozygous child error is a hom-DNM. The non-inherited
  ("culprit") allele drives the recurrent test: it is recurrent when its
  population allele count exceeds 1 (het DNM) or 2 (hom DNM).
- **Transmission singletons** — counted per allele from its in-trio occurrence
  count. An allele seen once and carried by the child is *novel*; seen once and
  not in the child is an *untransmitted* singleton; seen exactly twice as a
  het-parent → het-child pair is a *transmitted* singleton (doubleton).
- **Ts/Tv** — per site (per trio), set if any distinct single-base ALT allele in
  the trio is a transition (`|acgt2int(ref) - acgt2int(alt)| == 2`, i.e. A↔G or
  C↔T) or a transversion; indel and `*` alleles are excluded.

Population allele counts use `INFO/AC`+`INFO/AN` when both are present, else are
derived from the genotypes of every sample, matching htslib
`bcf_calc_ac(BCF_UN_INFO|BCF_UN_FMT)`. Ploidy is respected: a hemizygous call
(e.g. a chrX/Y `GT` of `1`) contributes a single allele copy to this count, so
it does not inflate the culprit allele's count in the recurrent-DNM test.

## Origin

This crate is an independent Rust reimplementation of the `bcftools +trio-stats`
plugin. The plugin is MIT-licensed, so its source was read and cited directly:
the per-trio counting loop, the Mendelian/de-novo classification, the
`bcf_acgt2int`-based Ts/Tv test, the PED parsing and trio detection, and the
output column set all follow `plugins/trio-stats.c`.

Test fixtures are independently authored VCF + PED files; the expected counts in
the compatibility test were captured from `bcftools +trio-stats 1.23.1`.

License: MIT OR Apache-2.0.
Upstream credit: bcftools <https://github.com/samtools/bcftools> (MIT, plugin by
Petr Danecek, Genome Research Ltd).
