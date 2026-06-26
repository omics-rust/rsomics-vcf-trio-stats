//! PED parsing and trio detection, matching `trio-stats.c::parse_ped`.
//!
//! A PED line is whitespace-split into at least four columns: family, sample,
//! paternal, maternal (sex/phenotype and beyond are ignored). A complete trio
//! requires the sample, paternal and maternal IDs to all be present in the VCF
//! header. Duplicate trios (identical child/father/mother triples) are skipped
//! with a one-time warning; a child listed in two distinct trios is a fatal
//! error. Trios are sorted by their smallest VCF sample column so the streaming
//! scan touches genotype columns roughly in order, reproducing the plugin's
//! `qsort(cmp_trios)`.

use std::collections::HashMap;
use std::path::Path;

use rsomics_common::{Result, RsomicsError};

/// One detected trio: the VCF sample column indices of child, father and mother
/// plus their names for the report rows.
#[derive(Debug, Clone)]
pub struct Trio {
    pub child: usize,
    pub father: usize,
    pub mother: usize,
    pub child_name: String,
    pub father_name: String,
    pub mother_name: String,
}

fn index_of<'a>(names: &'a HashMap<&'a str, usize>, id: &str) -> Option<usize> {
    names.get(id).copied()
}

/// Parse the PED file and return the complete trios found in `sample_names`
/// (the VCF `#CHROM` sample order), in the plugin's sorted order.
pub fn parse_ped(path: &Path, sample_names: &[String]) -> Result<Vec<Trio>> {
    let text = std::fs::read_to_string(path).map_err(|e| {
        RsomicsError::Io(std::io::Error::new(
            e.kind(),
            format!("cannot read PED {}: {e}", path.display()),
        ))
    })?;

    let index: HashMap<&str, usize> = sample_names
        .iter()
        .enumerate()
        .map(|(i, n)| (n.as_str(), i))
        .collect();

    let mut trios: Vec<Trio> = Vec::new();
    let mut seen_child: HashMap<String, ()> = HashMap::new();
    let mut seen_trio: HashMap<String, ()> = HashMap::new();
    let mut dup_warned = false;

    for line in text.lines() {
        let cols: Vec<&str> = line.split_whitespace().collect();
        if cols.is_empty() {
            continue;
        }
        if cols.len() < 4 {
            return Err(RsomicsError::InvalidInput(format!(
                "could not parse the PED line: {line}"
            )));
        }
        let (sample, pat, mat) = (cols[1], cols[2], cols[3]);

        let father = match index_of(&index, pat) {
            Some(i) => i,
            None => continue,
        };
        let mother = match index_of(&index, mat) {
            Some(i) => i,
            None => continue,
        };
        let child = match index_of(&index, sample) {
            Some(i) => i,
            None => continue,
        };

        let trio_key = format!("{sample} {pat} {mat}");
        if seen_trio.contains_key(&trio_key) {
            if !dup_warned {
                eprintln!(
                    "Warning: the trio \"{trio_key}\" is listed multiple times, skipping. \
                     (This message is printed only once.)"
                );
                dup_warned = true;
            }
            continue;
        }
        if seen_child.contains_key(sample) {
            return Err(RsomicsError::InvalidInput(format!(
                "the child \"{sample}\" is listed in two trios"
            )));
        }
        seen_child.insert(sample.to_string(), ());
        seen_trio.insert(trio_key, ());

        trios.push(Trio {
            child,
            father,
            mother,
            child_name: sample.to_string(),
            father_name: pat.to_string(),
            mother_name: mat.to_string(),
        });
    }

    if trios.is_empty() {
        return Err(RsomicsError::InvalidInput(
            "no complete trio identified in the VCF for this PED".into(),
        ));
    }

    trios.sort_by_key(|t| t.child.min(t.father).min(t.mother));
    Ok(trios)
}

/// Build a single trio from explicit proband,father,mother sample names (`-P`).
pub fn pfm_trio(pfm: &str, sample_names: &[String]) -> Result<Vec<Trio>> {
    let parts: Vec<&str> = pfm.split(',').collect();
    if parts.len() != 3 {
        return Err(RsomicsError::InvalidInput(format!(
            "could not parse --pfm {pfm:?}; expected proband,father,mother"
        )));
    }
    let index: HashMap<&str, usize> = sample_names
        .iter()
        .enumerate()
        .map(|(i, n)| (n.as_str(), i))
        .collect();
    let lookup = |id: &str| {
        index_of(&index, id)
            .ok_or_else(|| RsomicsError::InvalidInput(format!("no such sample: {id:?}")))
    };
    let child = lookup(parts[0])?;
    let father = lookup(parts[1])?;
    let mother = lookup(parts[2])?;
    Ok(vec![Trio {
        child,
        father,
        mother,
        child_name: parts[0].to_string(),
        father_name: parts[1].to_string(),
        mother_name: parts[2].to_string(),
    }])
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;

    fn names(v: &[&str]) -> Vec<String> {
        v.iter().map(|s| (*s).to_string()).collect()
    }

    fn write_ped(content: &str) -> tempfile::NamedTempFile {
        let mut f = tempfile::NamedTempFile::new_in(std::env::temp_dir()).unwrap();
        f.write_all(content.as_bytes()).unwrap();
        f.flush().unwrap();
        f
    }

    #[test]
    fn detects_complete_trio() {
        let ped = write_ped("FAM C F M 1 0\n");
        let trios = parse_ped(ped.path(), &names(&["F", "M", "C"])).unwrap();
        assert_eq!(trios.len(), 1);
        assert_eq!(trios[0].child_name, "C");
        assert_eq!(trios[0].child, 2);
        assert_eq!(trios[0].father, 0);
        assert_eq!(trios[0].mother, 1);
    }

    #[test]
    fn skips_incomplete_trio() {
        // Mother not in VCF -> not a complete trio -> no trios -> error.
        let ped = write_ped("FAM C F MISSING 1 0\n");
        let err = parse_ped(ped.path(), &names(&["F", "C"])).unwrap_err();
        assert!(format!("{err}").contains("no complete trio"));
    }

    #[test]
    fn duplicate_trio_skipped_not_fatal() {
        let ped = write_ped("FAM C F M\nFAM C F M\n");
        let trios = parse_ped(ped.path(), &names(&["F", "M", "C"])).unwrap();
        assert_eq!(trios.len(), 1);
    }

    #[test]
    fn child_in_two_trios_is_fatal() {
        let ped = write_ped("FAM C F M\nFAM C F2 M2\n");
        let err = parse_ped(ped.path(), &names(&["F", "M", "C", "F2", "M2"])).unwrap_err();
        assert!(format!("{err}").contains("two trios"));
    }

    #[test]
    fn trios_sorted_by_min_column() {
        // Trio A children at columns {5,6,7}, trio B at {0,1,2}: B sorts first.
        let ped = write_ped("FAM CA FA MA\nFAM CB FB MB\n");
        let trios = parse_ped(
            ped.path(),
            &names(&["FB", "MB", "CB", "x", "y", "FA", "MA", "CA"]),
        )
        .unwrap();
        assert_eq!(trios[0].child_name, "CB");
        assert_eq!(trios[1].child_name, "CA");
    }

    #[test]
    fn pfm_explicit() {
        let trios = pfm_trio("C,F,M", &names(&["F", "M", "C"])).unwrap();
        assert_eq!(trios.len(), 1);
        assert_eq!(trios[0].child, 2);
    }
}
