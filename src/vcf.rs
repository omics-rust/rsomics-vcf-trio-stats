//! Streaming VCF reader driving the per-trio stats.
//!
//! The `#CHROM` line fixes the sample column order, against which the PED trios
//! are resolved. Every record is walked once over its tab fields: REF/ALT build
//! the allele list, INFO yields AC/AN (preferred when both are present, matching
//! htslib `bcf_calc_ac`'s `BCF_UN_INFO|BCF_UN_FMT` priority), and each sample's
//! GT token is decoded into allele indices. The population allele count drives
//! the recurrent-DNM test and is computed over every sample in the file. Only
//! the genotype columns belonging to a trio member are decoded per record.

use std::collections::BTreeSet;
use std::fs::File;
use std::io::{BufRead, BufReader, Read};
use std::path::Path;

use flate2::read::MultiGzDecoder;
use rsomics_common::{Result, RsomicsError};

use crate::ped::Trio;
use crate::triostats::{Gt, Record, TrioStats, accumulate};

fn open(path: &Path) -> Result<Box<dyn BufRead>> {
    let file = File::open(path).map_err(|e| {
        RsomicsError::Io(std::io::Error::new(
            e.kind(),
            format!("cannot open {}: {e}", path.display()),
        ))
    })?;
    let is_gz = path
        .extension()
        .is_some_and(|e| e.eq_ignore_ascii_case("gz"));
    let reader: Box<dyn Read> = if is_gz {
        Box::new(MultiGzDecoder::new(file))
    } else {
        Box::new(file)
    };
    Ok(Box::new(BufReader::new(reader)))
}

/// All sample names in `#CHROM` column order.
pub fn read_sample_names(path: &Path) -> Result<Vec<String>> {
    let mut reader = open(path)?;
    let mut buf = String::new();
    loop {
        buf.clear();
        let n = reader
            .read_line(&mut buf)
            .map_err(|e| RsomicsError::Io(std::io::Error::new(e.kind(), "reading VCF header")))?;
        if n == 0 {
            break;
        }
        if buf.starts_with("#CHROM") {
            let fields: Vec<&str> = buf.trim_end().split('\t').collect();
            if fields.len() < 10 {
                return Err(RsomicsError::InvalidInput(
                    "#CHROM header has no sample columns; trio-stats needs genotypes".into(),
                ));
            }
            return Ok(fields[9..].iter().map(|s| (*s).to_string()).collect());
        }
        if !buf.starts_with('#') {
            break;
        }
    }
    Err(RsomicsError::InvalidInput(
        "VCF has no #CHROM header line".into(),
    ))
}

/// Decode a GT token (e.g. `0/1`, `1`, `./.`, `1|2`) into allele indices,
/// reproducing `parse_genotype`: a missing first or second allele yields `None`;
/// a haploid call (single allele, or one terminated early) is homozygous with
/// both slots equal; otherwise the diploid pair is returned. Phasing is ignored.
fn parse_gt(gt: &[u8]) -> Option<Gt> {
    let mut als = [0i32; 2];
    let mut slot = 0usize;
    let mut cur: i64 = -1;
    let mut cur_missing = false;
    let mut seen = false;

    let mut flush = |slot: &mut usize, cur: i64, missing: bool, seen: bool| -> Option<()> {
        if !seen || *slot >= 2 {
            return Some(());
        }
        if missing {
            return None;
        }
        als[*slot] = cur as i32;
        *slot += 1;
        Some(())
    };

    for &b in gt {
        match b {
            b'/' | b'|' => {
                flush(&mut slot, cur, cur_missing, seen)?;
                cur = -1;
                cur_missing = false;
                seen = false;
            }
            b'0'..=b'9' => {
                if cur < 0 {
                    cur = 0;
                }
                cur = cur * 10 + i64::from(b - b'0');
                seen = true;
            }
            b'.' => {
                cur_missing = true;
                seen = true;
            }
            _ => {}
        }
    }
    flush(&mut slot, cur, cur_missing, seen)?;

    match slot {
        0 => None,
        1 => Some(Gt {
            als: [als[0], als[0]],
        }),
        _ => Some(Gt { als }),
    }
}

/// Locate the `GT` subfield index in a FORMAT field.
fn gt_index(format: &[u8]) -> Option<usize> {
    format.split(|&b| b == b':').position(|f| f == b"GT")
}

/// Extract the GT token from a sample field given the GT subfield index.
fn gt_token(sample: &[u8], gt_idx: usize) -> &[u8] {
    sample.split(|&b| b == b':').nth(gt_idx).unwrap_or(b".")
}

/// Pull a comma-list integer INFO tag out of the INFO column. `None` if absent.
fn info_ints(info: &[u8], key: &[u8]) -> Option<Vec<i64>> {
    if info == b"." {
        return None;
    }
    for kv in info.split(|&b| b == b';') {
        let (k, v) = match kv.iter().position(|&b| b == b'=') {
            Some(p) => (&kv[..p], &kv[p + 1..]),
            None => (kv, &b""[..]),
        };
        if k == key {
            let mut out = Vec::new();
            for tok in v.split(|&b| b == b',') {
                let n: i64 = std::str::from_utf8(tok).ok()?.trim().parse().ok()?;
                out.push(n);
            }
            return Some(out);
        }
    }
    None
}

/// Fill `ac` with the per-allele population count, reusing the buffer. Prefers
/// INFO/AC+AN when both parse (`ac[0]=AN-sum(AC)`, `ac[i+1]=AC[i]`); otherwise
/// counts every non-missing allele occurrence across all samples, mirroring
/// `bcf_calc_ac(BCF_UN_INFO|BCF_UN_FMT)`.
fn fill_allele_counts(ac: &mut Vec<i64>, n_allele: usize, info: &[u8], all_gts: &[Option<Gt>]) {
    ac.clear();
    ac.resize(n_allele, 0);

    if let (Some(an), Some(ac_alt)) = (info_ints(info, b"AN"), info_ints(info, b"AC"))
        && an.len() == 1
        && ac_alt.len() == n_allele - 1
    {
        let sum_alt: i64 = ac_alt.iter().sum();
        ac[0] = an[0] - sum_alt;
        for (i, &c) in ac_alt.iter().enumerate() {
            ac[i + 1] = c;
        }
        return;
    }
    for gt in all_gts.iter().flatten() {
        for &a in &gt.als {
            if a >= 0 && (a as usize) < n_allele {
                ac[a as usize] += 1;
            }
        }
    }
}

/// Accumulate per-trio stats over the whole VCF. The scratch buffers
/// (`all_gts`, `alleles`, `ac`, `ac_trio`) are reused across records so the
/// streaming loop performs no per-line heap allocation.
pub fn read_stats(path: &Path, trios: &[Trio]) -> Result<Vec<TrioStats>> {
    // The set of columns any trio touches; only these are decoded per record.
    let needed: BTreeSet<usize> = trios
        .iter()
        .flat_map(|t| [t.child, t.father, t.mother])
        .collect();
    let max_col = needed.iter().copied().max().unwrap_or(0);

    let mut stats = vec![TrioStats::default(); trios.len()];
    let mut reader = open(path)?;
    let mut all_gts: Vec<Option<Gt>> = Vec::new();
    let mut ac: Vec<i64> = Vec::new();
    let mut ac_trio: Vec<i64> = Vec::new();

    let mut block = vec![0u8; 1 << 20];
    let mut carry: Vec<u8> = Vec::new();
    loop {
        let n = reader
            .read(&mut block)
            .map_err(|e| RsomicsError::Io(std::io::Error::new(e.kind(), "reading VCF")))?;
        if n == 0 {
            break;
        }
        let mut data = &block[..n];
        if !carry.is_empty() {
            if let Some(nl) = data.iter().position(|&b| b == b'\n') {
                carry.extend_from_slice(&data[..nl]);
                process_line(
                    &carry,
                    trios,
                    max_col,
                    &mut stats,
                    &mut all_gts,
                    &mut ac,
                    &mut ac_trio,
                )?;
                carry.clear();
                data = &data[nl + 1..];
            } else {
                carry.extend_from_slice(data);
                continue;
            }
        }
        let mut start = 0usize;
        while let Some(off) = data[start..].iter().position(|&b| b == b'\n') {
            let line = &data[start..start + off];
            process_line(
                line,
                trios,
                max_col,
                &mut stats,
                &mut all_gts,
                &mut ac,
                &mut ac_trio,
            )?;
            start += off + 1;
        }
        carry.extend_from_slice(&data[start..]);
    }
    if !carry.is_empty() {
        process_line(
            &carry,
            trios,
            max_col,
            &mut stats,
            &mut all_gts,
            &mut ac,
            &mut ac_trio,
        )?;
    }
    Ok(stats)
}

#[allow(clippy::too_many_arguments)]
fn process_line(
    line: &[u8],
    trios: &[Trio],
    max_col: usize,
    stats: &mut [TrioStats],
    all_gts: &mut Vec<Option<Gt>>,
    ac: &mut Vec<i64>,
    ac_trio: &mut Vec<i64>,
) -> Result<()> {
    if line.first() == Some(&b'#') {
        return Ok(());
    }
    let mut bytes = line;
    while bytes.last() == Some(&b'\n') || bytes.last() == Some(&b'\r') {
        bytes = &bytes[..bytes.len() - 1];
    }
    if bytes.is_empty() {
        return Ok(());
    }

    let mut reff: &[u8] = b"";
    let mut alt: &[u8] = b"";
    let mut info: &[u8] = b".";
    all_gts.clear();

    let mut gt_idx = usize::MAX;
    let mut col = 0usize;
    let mut start = 0usize;
    for i in 0..=bytes.len() {
        if i == bytes.len() || bytes[i] == b'\t' {
            let field = &bytes[start..i];
            match col {
                3 => reff = field,
                4 => alt = field,
                7 => info = field,
                8 => {
                    gt_idx = gt_index(field).ok_or_else(|| {
                        RsomicsError::InvalidInput("FORMAT has no GT field".into())
                    })?;
                }
                c if c >= 9 => all_gts.push(parse_gt(gt_token(field, gt_idx))),
                _ => {}
            }
            start = i + 1;
            col += 1;
        }
    }
    if col < 10 {
        return Err(RsomicsError::InvalidInput(
            "trio-stats needs per-sample genotypes; record has fewer than 10 columns".into(),
        ));
    }
    if all_gts.len() <= max_col {
        return Err(RsomicsError::InvalidInput(format!(
            "sample column {max_col} out of range in record"
        )));
    }

    let two: [&[u8]; 2] = [reff, alt];
    let mut multi: Vec<&[u8]> = Vec::new();
    let allele_slice: &[&[u8]] = if alt == b"." {
        &two[..1]
    } else if !alt.contains(&b',') {
        &two
    } else {
        multi.push(reff);
        for a in alt.split(|&b| b == b',') {
            multi.push(a);
        }
        &multi
    };
    let n_allele = allele_slice.len();

    fill_allele_counts(ac, n_allele, info, all_gts);
    if ac_trio.len() != n_allele {
        ac_trio.resize(n_allele, 0);
    }
    let rec = Record::borrowed(allele_slice, ac);

    for (out_i, trio) in trios.iter().enumerate() {
        accumulate(
            &mut stats[out_i],
            &rec,
            all_gts[trio.child],
            all_gts[trio.father],
            all_gts[trio.mother],
            ac_trio,
        );
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parse_gt_diploid() {
        assert_eq!(parse_gt(b"0/1"), Some(Gt { als: [0, 1] }));
        assert_eq!(parse_gt(b"1|2"), Some(Gt { als: [1, 2] }));
    }

    #[test]
    fn parse_gt_haploid_is_hom() {
        assert_eq!(parse_gt(b"1"), Some(Gt { als: [1, 1] }));
        assert_eq!(parse_gt(b"0"), Some(Gt { als: [0, 0] }));
    }

    #[test]
    fn parse_gt_missing() {
        assert_eq!(parse_gt(b"./."), None);
        assert_eq!(parse_gt(b"."), None);
        assert_eq!(parse_gt(b"0/."), None);
        assert_eq!(parse_gt(b"./1"), None);
    }

    #[test]
    fn info_ac_an() {
        let info = b"DP=30;AC=2,1;AN=6";
        assert_eq!(info_ints(info, b"AC"), Some(vec![2, 1]));
        assert_eq!(info_ints(info, b"AN"), Some(vec![6]));
        assert_eq!(info_ints(info, b"XX"), None);
    }

    #[test]
    fn ac_from_genotypes_when_no_info() {
        let gts = vec![Some(Gt { als: [0, 1] }), Some(Gt { als: [1, 1] })];
        let mut ac = Vec::new();
        fill_allele_counts(&mut ac, 2, b".", &gts);
        assert_eq!(ac, vec![1, 3]);
    }
}
