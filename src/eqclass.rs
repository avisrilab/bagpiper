//! Equivalence-class molecules: the packed key built from aligned reads.

use std::fs::File;
use std::io::{self, Write};
use std::path::Path;

use noodles::bam;

use crate::seq::{CellId, Umi};

/// One molecule: cell, UMI, and its transcript-id set (sorted, reference order). Exact equality of
/// all three fields defines a PCR duplicate.
#[derive(Clone, PartialEq, Eq, PartialOrd, Ord, Hash, Debug)]
pub struct Molecule {
    pub cell: CellId,
    pub umi: Umi,
    pub txps: Vec<u32>,
}

/// Parse a barcoded read name into its cell and molecular UMI. Two layouts, distinguished by field
/// count (the original id is assumed underscore-free):
/// - v4 `origid_CB_UMI` (3 fields): the UMI is the 3' UMI.
/// - v5 `origid_CB_algoUMI_umi1_umi2` (5 fields): the molecular UMI is the two 5' TSO 6-mers
///   concatenated, dropping the barcode-step algo-UMI. With `v5_binid`, the 3 bp binning index (the
///   last 3 nt of the 12 nt algo-UMI slot) is appended, giving the 15 nt key the `collapse` stage
///   slices back.
///
/// None on a malformed barcode/UMI, an unexpected field count, or (under `v5_binid`) an algo-UMI slot
/// that is not 12 nt.
pub fn parse_key(read_name: &[u8], v5_binid: bool) -> Option<(CellId, Umi)> {
    let f: Vec<&[u8]> = read_name.split(|&b| b == b'_').collect();
    match f.len() {
        3 => Some((CellId::from_ascii(f[1])?, Umi::from_ascii(f[2])?)),
        5 => {
            let cell = CellId::from_ascii(f[1])?;
            let mut key = f[3].to_vec();
            key.extend_from_slice(f[4]);
            if v5_binid {
                if f[2].len() != 12 {
                    return None;
                }
                key.extend_from_slice(&f[2][9..12]);
            }
            Some((cell, Umi::from_ascii(&key)?))
        }
        _ => None,
    }
}

/// An eqclass built from a BAM: every reference transcript (name, length) in reference order (index
/// = transcript id), and one molecule per kept read group.
pub struct EqClass {
    pub transcripts: Vec<(String, u32)>,
    pub molecules: Vec<Molecule>,
}

impl EqClass {
    /// Write the text form: transcript count, one `name\tlength` line per transcript, then one
    /// `CB\tUMI\ttxp...` line per molecule with transcripts named in reference order.
    pub fn write_text<W: Write>(&self, out: &mut W) -> io::Result<()> {
        writeln!(out, "{}", self.transcripts.len())?;
        for (name, len) in &self.transcripts {
            writeln!(out, "{}\t{}", name, len)?;
        }
        for m in &self.molecules {
            write!(out, "{}\t{}", m.cell.render(), m.umi.render())?;
            for &tid in &m.txps {
                write!(out, "\t{}", self.transcripts[tid as usize].0)?;
            }
            writeln!(out)?;
        }
        Ok(())
    }
}

/// Build the eqclass from a name-grouped BAM. Reads are grouped by name; supplementary alignments
/// are ignored; a group whose first non-supplementary alignment is unmapped is dropped. The
/// transcript-id set of the remaining alignments is sorted into the molecule. The read name (v4 or
/// v5) gives the cell and molecular UMI; `v5_binid` packs the V5 binning index (see [`parse_key`]).
pub fn read_bam<P: AsRef<Path>>(path: P, v5_binid: bool) -> io::Result<EqClass> {
    let mut reader = File::open(path).map(bam::io::Reader::new)?;
    let header = reader.read_header()?;

    let transcripts: Vec<(String, u32)> = header
        .reference_sequences()
        .iter()
        .map(|(name, rs)| (name.to_string(), rs.length().get() as u32))
        .collect();

    let mut molecules = Vec::new();
    let mut cur_name: Vec<u8> = Vec::new();
    let mut cur_tids: Vec<u32> = Vec::new();
    let mut first_unmapped: Option<bool> = None;
    let mut have_group = false;

    for result in reader.records() {
        let record = result?;
        let name: &[u8] = match record.name() {
            Some(n) => n.as_ref(),
            None => continue,
        };

        if !have_group || name != cur_name.as_slice() {
            if have_group {
                flush_group(&cur_name, &mut cur_tids, first_unmapped, v5_binid, &mut molecules);
            }
            cur_name.clear();
            cur_name.extend_from_slice(name);
            cur_tids.clear();
            first_unmapped = None;
            have_group = true;
        }

        let flags = record.flags();
        if flags.is_supplementary() {
            continue;
        }
        if first_unmapped.is_none() {
            first_unmapped = Some(flags.is_unmapped());
        }
        if let Some(id) = record.reference_sequence_id() {
            cur_tids.push(id? as u32);
        }
    }
    if have_group {
        flush_group(&cur_name, &mut cur_tids, first_unmapped, v5_binid, &mut molecules);
    }

    Ok(EqClass {
        transcripts,
        molecules,
    })
}

/// Emit one molecule for a completed read group, unless it was dropped: no non-supplementary
/// alignment, first non-supplementary unmapped, or a malformed key.
fn flush_group(
    name: &[u8],
    tids: &mut Vec<u32>,
    first_unmapped: Option<bool>,
    v5_binid: bool,
    out: &mut Vec<Molecule>,
) {
    if first_unmapped != Some(false) {
        return;
    }
    let (cell, umi) = match parse_key(name, v5_binid) {
        Some(key) => key,
        None => return,
    };
    let mut txps = std::mem::take(tids);
    txps.sort_unstable();
    out.push(Molecule { cell, umi, txps });
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parse_v4_key() {
        let (cell, umi) = parse_key(b"readid_AAACGTTGCAGAACAC_ACGTACGTACGT", false).unwrap();
        assert_eq!(cell.render(), "AAACGTTGCAGAACAC");
        assert_eq!(umi.render(), "ACGTACGTACGT");
    }

    #[test]
    fn parse_v5_key_concatenates_tso_umis() {
        // origid_CB_algoUMI_umi1_umi2 -> CB, umi1+umi2 (algo-UMI dropped)
        let name = b"rid_AAACGTTGCAGAACAC_AAAAAACGTACG_GTAGTT_AGGGGG";
        let (cell, umi) = parse_key(name, false).unwrap();
        assert_eq!(cell.render(), "AAACGTTGCAGAACAC");
        assert_eq!(umi.render(), "GTAGTTAGGGGG");
    }

    #[test]
    fn parse_v5_binid_appends_bin_index() {
        // with v5_binid, append the last 3 nt of the 12 nt algo-UMI slot (AAAAAACGTACG -> ACG)
        let name = b"rid_AAACGTTGCAGAACAC_AAAAAACGTACG_GTAGTT_AGGGGG";
        assert_eq!(parse_key(name, true).unwrap().1.render(), "GTAGTTAGGGGGACG");
        // a wrong-length algo-UMI slot is dropped, not packed
        assert!(parse_key(b"rid_AAACGTTGCAGAACAC_ACG_GTAGTT_AGGGGG", true).is_none());
    }

    #[test]
    fn parse_rejects_unexpected_field_count() {
        assert!(parse_key(b"readid_CB_UMI_extra", false).is_none());
        assert!(parse_key(b"nounderscore", false).is_none());
    }

    #[test]
    fn parse_rejects_non_acgt() {
        assert!(parse_key(b"readid_AAACGTTGCAGAACAC_ACGTNCGTACGT", false).is_none());
    }
}
