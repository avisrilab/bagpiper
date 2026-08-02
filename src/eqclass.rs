//! Equivalence-class molecules: the packed key built from aligned reads.

use std::fs::File;
use std::io::{self, Write};
use std::path::Path;

use noodles::bam;

use crate::seq::{CellId, Umi};

/// One molecule: cell, UMI, and its aligned transcript ids in BAM (alignment) order. Exact equality
/// of all three fields, transcript order included, defines a PCR duplicate: the reference dedup is a
/// whole-line string match, so two reads of one molecule whose shared transcript set is listed in a
/// different order are kept as distinct molecules.
#[derive(Clone, PartialEq, Eq, PartialOrd, Ord, Hash, Debug)]
pub struct Molecule {
    pub cell: CellId,
    pub umi: Umi,
    pub txps: Vec<u32>,
}

/// Parse a v4 read name `origid_CB_UMI` into its cell and UMI. None on a malformed barcode/UMI or
/// an unexpected field count; the original id is assumed underscore-free.
pub fn parse_key(read_name: &[u8]) -> Option<(CellId, Umi)> {
    let fields: Vec<&[u8]> = read_name.split(|&b| b == b'_').collect();
    if fields.len() != 3 {
        return None;
    }
    let cell = CellId::from_ascii(fields[1])?;
    let umi = Umi::from_ascii(fields[2])?;
    Some((cell, umi))
}

/// An eqclass built from a BAM: every reference transcript (name, length) in reference order (index
/// = transcript id), and one molecule per kept read group.
pub struct EqClass {
    pub transcripts: Vec<(String, u32)>,
    pub molecules: Vec<Molecule>,
}

impl EqClass {
    /// Write the text form: transcript count, one `name\tlength` line per transcript, then one
    /// `CB\tUMI\ttxp...` line per molecule with transcripts named in alignment order.
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
/// transcript ids of the remaining alignments are recorded in BAM (alignment) order.
pub fn read_bam<P: AsRef<Path>>(path: P) -> io::Result<EqClass> {
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
                flush_group(&cur_name, &mut cur_tids, first_unmapped, &mut molecules);
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
        flush_group(&cur_name, &mut cur_tids, first_unmapped, &mut molecules);
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
    out: &mut Vec<Molecule>,
) {
    if first_unmapped != Some(false) {
        return;
    }
    let (cell, umi) = match parse_key(name) {
        Some(key) => key,
        None => return,
    };
    // Alignment order, unsorted: dedup keys on the whole molecule (order included), matching the
    // reference. The transcript set is sorted only later, when count builds the eqclass key.
    let txps = std::mem::take(tids);
    out.push(Molecule { cell, umi, txps });
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parse_v4_key() {
        let (cell, umi) = parse_key(b"readid_AAACGTTGCAGAACAC_ACGTACGTACGT").unwrap();
        assert_eq!(cell.render(), "AAACGTTGCAGAACAC");
        assert_eq!(umi.render(), "ACGTACGTACGT");
    }

    #[test]
    fn parse_rejects_unexpected_field_count() {
        assert!(parse_key(b"readid_CB_UMI_extra").is_none());
        assert!(parse_key(b"nounderscore").is_none());
    }

    #[test]
    fn parse_rejects_non_acgt() {
        assert!(parse_key(b"readid_AAACGTTGCAGAACAC_ACGTNCGTACGT").is_none());
    }
}
