//! Per-trio transmission and de-novo counters, value-exact with bcftools
//! `+trio-stats`.
//!
//! For each record every complete trio's three genotypes (child, father,
//! mother) are folded into one of the trio's counters. The arithmetic mirrors
//! `plugins/trio-stats.c::process_record`: the per-trio allele tally over the
//! six GT slots, the non-ref test, the `bcf_acgt2int` transition test
//! `abs(ref-alt)==2` over the distinct ALT alleles in the trio, the Mendelian
//! cross-check, the homozygous/recurrent DNM split keyed on the population
//! allele count, and the singleton/doubleton transmission test against the
//! per-trio allele count.

/// One trio's accumulated counters. Field order matches the plugin's `FLT` row
/// after the three sample names.
#[derive(Debug, Clone, Default, PartialEq, Eq, serde::Serialize)]
pub struct TrioStats {
    pub npass: u32,
    pub nnon_ref: u32,
    pub nmendel_err: u32,
    pub nnovel: u32,
    pub nsingleton: u32,
    pub ndoubleton: u32,
    pub nts: u32,
    pub ntv: u32,
    pub ndnm_hom: u32,
    pub ndnm_recurrent: u32,
}

impl TrioStats {
    /// The trio's Ts/Tv ratio, infinite when there are no transversions,
    /// matching the plugin's `ntv ? (float)nts/ntv : INFINITY`.
    #[must_use]
    pub fn tstv(&self) -> f64 {
        if self.ntv == 0 {
            f64::INFINITY
        } else {
            f64::from(self.nts) / f64::from(self.ntv)
        }
    }
}

/// htslib `bcf_acgt2int`: A/a→0, C/c→1, G/g→2, T/t→3, anything else→-1. The
/// plugin's `abs(ref-alt)==2` test then selects exactly the A↔G and C↔T
/// transitions; an `N` or other ambiguous base maps to -1, so every comparison
/// against it lands in the transversion branch.
#[must_use]
pub fn acgt2int(b: u8) -> i32 {
    match b {
        b'A' | b'a' => 0,
        b'C' | b'c' => 1,
        b'G' | b'g' => 2,
        b'T' | b't' => 3,
        _ => -1,
    }
}

/// The allele indices of one sample's GT after `parse_genotype`: `None` is a
/// missing genotype; a haploid call is treated as homozygous diploid for the
/// per-trio tally and Mendelian test (the plugin sets `als[1]=als[0]`), while
/// `ploidy` records the true call length so the population allele count over
/// `bcf_calc_ac` charges a hemizygous call one copy, not two.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Gt {
    pub als: [i32; 2],
    pub ploidy: u8,
}

/// A record's parsed alleles, the per-allele population count `ac` (indexed by
/// allele number, computed over every sample), the numeric REF code for the
/// Ts/Tv test (`-1` when REF is not a single base), and which allele index, if
/// any, is the `*` overlapping-deletion symbol. All borrowed from streaming
/// scratch so a record costs no allocation.
pub struct Record<'a> {
    pub alleles: &'a [&'a [u8]],
    pub ac: &'a [i64],
    pub ref_code: i32,
    pub star_allele: Option<i32>,
}

impl<'a> Record<'a> {
    #[must_use]
    pub fn borrowed(alleles: &'a [&'a [u8]], ac: &'a [i64]) -> Self {
        let ref_allele = alleles[0];
        let ref_code = if ref_allele.len() == 1 {
            acgt2int(ref_allele[0])
        } else {
            -1
        };
        let star_allele = alleles
            .iter()
            .position(|a| a.len() == 1 && a[0] == b'*')
            .map(|i| i as i32);
        Record {
            alleles,
            ac,
            ref_code,
            star_allele,
        }
    }

    /// Is `allele` a single-base SNV allele (not the `*` symbol, not an indel)?
    fn is_snv_base(&self, allele: i32) -> bool {
        let a = &self.alleles[allele as usize];
        a.len() == 1 && a[0] != b'*'
    }
}

/// Fold one record's trio (child, father, mother genotypes) into the trio's
/// accumulating stats, exactly as `process_record`'s per-trio body does. A
/// missing genotype in any member skips the trio at this site. `ac_trio` is a
/// reusable scratch buffer of length `rec.alleles.len()`.
pub fn accumulate(
    stats: &mut TrioStats,
    rec: &Record,
    child: Option<Gt>,
    father: Option<Gt>,
    mother: Option<Gt>,
    ac_trio: &mut [i64],
) {
    let (child, father, mother) = match (child, father, mother) {
        (Some(c), Some(f), Some(m)) => (c, f, m),
        _ => return,
    };

    stats.npass += 1;

    let als = [
        child.als[0],
        child.als[1],
        father.als[0],
        father.als[1],
        mother.als[0],
        mother.als[1],
    ];
    let star = rec.star_allele;

    ac_trio.iter_mut().for_each(|c| *c = 0);
    let mut has_star_allele = false;
    let mut has_nonref = false;
    for &a in &als {
        if Some(a) == star {
            has_star_allele = true;
            continue;
        }
        if a != 0 {
            has_nonref = true;
        }
        ac_trio[a as usize] += 1;
    }
    if !has_nonref {
        return;
    }

    stats.nnon_ref += 1;

    if rec.ref_code != -1 {
        let mut has_ts = false;
        let mut has_tv = false;
        for &a in &als {
            if a == 0 || Some(a) == star {
                continue;
            }
            if !rec.is_snv_base(a) {
                continue;
            }
            let alt = acgt2int(rec.alleles[a as usize][0]);
            if (rec.ref_code - alt).abs() == 2 {
                has_ts = true;
            } else {
                has_tv = true;
            }
        }
        if has_ts {
            stats.nts += 1;
        }
        if has_tv {
            stats.ntv += 1;
        }
    }

    // The star allele was already accounted for at the deletion's primary
    // record; skip the remaining per-site stats to avoid double-counting.
    if has_star_allele {
        return;
    }

    let cf0 = child.als[0] == father.als[0] || child.als[0] == father.als[1];
    let cm1 = child.als[1] == mother.als[0] || child.als[1] == mother.als[1];
    if !cf0 || !cm1 {
        let cm0 = child.als[0] == mother.als[0] || child.als[0] == mother.als[1];
        let cf1 = child.als[1] == father.als[0] || child.als[1] == father.als[1];
        if !cm0 || !cf1 {
            stats.nmendel_err += 1;

            let mut dnm_hom = false;
            if child.als[0] == child.als[1] {
                stats.ndnm_hom += 1;
                dnm_hom = true;
            }

            let culprit = if !cf0 && !cm0 {
                child.als[0]
            } else if !cf1 && !cm1 {
                child.als[1]
            } else if rec.ac[child.als[0] as usize] < rec.ac[child.als[1] as usize] {
                child.als[0]
            } else {
                child.als[1]
            };
            let ac_culprit = rec.ac[culprit as usize];
            if (!dnm_hom && ac_culprit > 1) || (dnm_hom && ac_culprit > 2) {
                stats.ndnm_recurrent += 1;
            }
        }
    }

    for (j, &cnt) in ac_trio.iter().enumerate() {
        if cnt == 0 {
            continue;
        }
        let j = j as i32;
        if cnt == 1 {
            if child.als[0] == j || child.als[1] == j {
                stats.nnovel += 1;
            } else {
                stats.nsingleton += 1;
            }
        } else if cnt == 2 {
            let in_child = child.als[0] == j || child.als[1] == j;
            let child_hom = child.als[0] == j && child.als[1] == j;
            if !in_child || child_hom {
                continue;
            }
            let parent_hom = (father.als[0] == j && father.als[1] == j)
                || (mother.als[0] == j && mother.als[1] == j);
            if parent_hom {
                continue;
            }
            stats.ndoubleton += 1;
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn rec<'a>(alleles: &'a [&'a [u8]], ac: &'a [i64]) -> Record<'a> {
        Record::borrowed(alleles, ac)
    }

    fn gt(a: i32, b: i32) -> Option<Gt> {
        Some(Gt {
            als: [a, b],
            ploidy: 2,
        })
    }

    fn hap(a: i32) -> Option<Gt> {
        Some(Gt {
            als: [a, a],
            ploidy: 1,
        })
    }

    fn run(rec: &Record, c: Option<Gt>, f: Option<Gt>, m: Option<Gt>) -> TrioStats {
        let mut s = TrioStats::default();
        let mut ac_trio = vec![0i64; rec.alleles.len()];
        accumulate(&mut s, rec, c, f, m, &mut ac_trio);
        s
    }

    #[test]
    fn acgt2int_table() {
        assert_eq!(acgt2int(b'A'), 0);
        assert_eq!(acgt2int(b'C'), 1);
        assert_eq!(acgt2int(b'G'), 2);
        assert_eq!(acgt2int(b'T'), 3);
        assert_eq!(acgt2int(b'a'), 0);
        assert_eq!(acgt2int(b'N'), -1);
        assert_eq!((acgt2int(b'A') - acgt2int(b'G')).abs(), 2);
        assert_eq!((acgt2int(b'C') - acgt2int(b'T')).abs(), 2);
        assert_ne!((acgt2int(b'A') - acgt2int(b'C')).abs(), 2);
    }

    #[test]
    fn missing_member_skips_trio() {
        let r = rec(&[b"A", b"G"], &[3, 3]);
        let s = run(&r, gt(0, 1), None, gt(0, 1));
        assert_eq!(s, TrioStats::default());
    }

    #[test]
    fn ref_only_trio_passes_but_not_nonref() {
        let r = rec(&[b"A", b"G"], &[6, 0]);
        let s = run(&r, gt(0, 0), gt(0, 0), gt(0, 0));
        assert_eq!(s.npass, 1);
        assert_eq!(s.nnon_ref, 0);
    }

    #[test]
    fn clean_transmission_het_father() {
        // Father 0/1, mother 0/0, child 0/1. A>G transition. allele1 ac in
        // trio = 2 (father + child) -> transmitted doubleton.
        let r = rec(&[b"A", b"G"], &[5, 3]);
        let s = run(&r, gt(0, 1), gt(0, 1), gt(0, 0));
        assert_eq!(s.npass, 1);
        assert_eq!(s.nnon_ref, 1);
        assert_eq!(s.nmendel_err, 0);
        assert_eq!(s.nts, 1);
        assert_eq!(s.ntv, 0);
        assert_eq!(s.ndoubleton, 1);
        assert_eq!(s.nnovel, 0);
        assert_eq!(s.nsingleton, 0);
    }

    #[test]
    fn untransmitted_singleton() {
        // Father carries the only alt; child does not inherit it. allele1 ac in
        // trio = 1 (father only), not in child -> untransmitted singleton.
        let r = rec(&[b"A", b"C"], &[5, 1]);
        let s = run(&r, gt(0, 0), gt(0, 1), gt(0, 0));
        assert_eq!(s.nsingleton, 1);
        assert_eq!(s.ndoubleton, 0);
        assert_eq!(s.nnovel, 0);
        assert_eq!(s.nmendel_err, 0);
        assert_eq!(s.ntv, 1);
    }

    #[test]
    fn het_denovo_recurrent() {
        // Both parents 0/0, child 0/1: child allele 1 cannot come from a
        // parent -> Mendelian error, het (not hom). Population ac[1]=4 > 1 so
        // recurrent. allele1 ac in trio = 1, present in child -> novel.
        let r = rec(&[b"A", b"G"], &[8, 4]);
        let s = run(&r, gt(0, 1), gt(0, 0), gt(0, 0));
        assert_eq!(s.nmendel_err, 1);
        assert_eq!(s.ndnm_hom, 0);
        assert_eq!(s.ndnm_recurrent, 1);
        assert_eq!(s.nnovel, 1);
        assert_eq!(s.nsingleton, 0);
    }

    #[test]
    fn hom_denovo_non_recurrent() {
        // Parents 0/0, child 1/1: hom DNM. Population ac[1]=2; hom needs ac>2
        // to be recurrent, so not recurrent here.
        let r = rec(&[b"A", b"T"], &[10, 2]);
        let s = run(&r, gt(1, 1), gt(0, 0), gt(0, 0));
        assert_eq!(s.nmendel_err, 1);
        assert_eq!(s.ndnm_hom, 1);
        assert_eq!(s.ndnm_recurrent, 0);
        // allele1 ac in trio = 2 but child is hom for it -> neither novel nor
        // doubleton.
        assert_eq!(s.nnovel, 0);
        assert_eq!(s.ndoubleton, 0);
    }

    #[test]
    fn hom_denovo_recurrent_when_ac_high() {
        let r = rec(&[b"A", b"T"], &[10, 3]);
        let s = run(&r, gt(1, 1), gt(0, 0), gt(0, 0));
        assert_eq!(s.ndnm_hom, 1);
        assert_eq!(s.ndnm_recurrent, 1);
    }

    #[test]
    fn transition_and_transversion_hetaa() {
        // G>A,T het child 1/2, parents transmit. G->A transition, G->T
        // transversion: both flags set, ts and tv each +1.
        let r = rec(&[b"G", b"A", b"T"], &[2, 2, 2]);
        let s = run(&r, gt(1, 2), gt(0, 1), gt(0, 2));
        assert_eq!(s.nts, 1);
        assert_eq!(s.ntv, 1);
    }

    #[test]
    fn star_allele_skips_post_stats_but_counts_nonref_and_tstv_ref_only() {
        // Trio carries '*' (allele 2). Mendel/singleton stats are skipped, but
        // nnon_ref is incremented if a real alt is present.
        let r = rec(&[b"A", b"G", b"*"], &[5, 1, 2]);
        let s = run(&r, gt(0, 2), gt(0, 1), gt(0, 0));
        assert_eq!(s.npass, 1);
        assert_eq!(s.nnon_ref, 1);
        // Star present -> no mendel/singleton/doubleton/novel accounting.
        assert_eq!(s.nmendel_err, 0);
        assert_eq!(s.nsingleton, 0);
        assert_eq!(s.ndoubleton, 0);
        assert_eq!(s.nnovel, 0);
        // G is an alt SNV in the trio -> transition counted.
        assert_eq!(s.nts, 1);
    }

    #[test]
    fn indel_not_counted_in_tstv() {
        // AT>A indel: not a single-base alt, excluded from Ts/Tv.
        let r = rec(&[b"AT", b"A"], &[5, 1]);
        let s = run(&r, gt(0, 1), gt(0, 1), gt(0, 0));
        assert_eq!(s.nts, 0);
        // REF length 2 -> ref_code -1 -> Ts/Tv block skipped entirely.
        assert_eq!(s.ntv, 0);
        assert_eq!(s.ndoubleton, 1);
    }

    #[test]
    fn mendel_culprit_ambiguous_picks_lower_af_allele() {
        // Child 1/2, father 0/0, mother 0/0: both alleles are de novo. Neither
        // a0F/a0M nor a1F/a1M, so culprit falls to the AF tiebreak: ac[1]=1 <
        // ac[2]=5 -> culprit allele 1, ac 1, not recurrent (het).
        let r = rec(&[b"A", b"G", b"T"], &[6, 1, 5]);
        let s = run(&r, gt(1, 2), gt(0, 0), gt(0, 0));
        assert_eq!(s.nmendel_err, 1);
        assert_eq!(s.ndnm_hom, 0);
        assert_eq!(s.ndnm_recurrent, 0);
    }

    #[test]
    fn haploid_child_hom_dnm_not_recurrent_at_ac_two() {
        // chrX: child hemizygous alt (1), father hemizygous ref (0), mother het
        // (0/1). The child's alt is non-inherited -> hom DNM. With the population
        // count charging the haploid child one copy, ac[1]=2 (child+mother), and
        // a hom DNM needs ac>2 to be recurrent -> not recurrent.
        let r = rec(&[b"A", b"G"], &[2, 2]);
        let s = run(&r, hap(1), hap(0), gt(0, 1));
        assert_eq!(s.nmendel_err, 1);
        assert_eq!(s.ndnm_hom, 1);
        assert_eq!(s.ndnm_recurrent, 0);
    }

    #[test]
    fn tstv_ratio() {
        let s = TrioStats {
            nts: 6,
            ntv: 0,
            ..Default::default()
        };
        assert!(s.tstv().is_infinite());
        let s2 = TrioStats {
            nts: 6,
            ntv: 3,
            ..Default::default()
        };
        assert!((s2.tstv() - 2.0).abs() < 1e-12);
    }
}
