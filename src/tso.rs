//! V5 TSO seal: the 5' dual-UMI carried on the template-switch oligo. The read is fit against the
//! anchor `CGCAGAGT V6 TACTG V6 GAAT` (V = any base) with the shared [`crate::sw`] aligner; the two
//! 6 nt UMIs are read off the anchor, the cDNA is trimmed to after the seal, and a score floor
//! rejects spurious low-score seals. The extracted read id gains `_umi1_umi2` so the V5 eqclass key
//! is `CB + umi1umi2`.

use std::fs::File;
use std::io::{self, BufWriter, Write};
use std::path::Path;

use flate2::write::GzEncoder;
use flate2::Compression;
use log::info;

use crate::fastq;
use crate::parallel;
use crate::sw::{advance, fit};

/// TSO anchor: CGCAGAGT, a 6 nt UMI, TACTG, a 6 nt UMI, GAAT.
const TSO: &[u8] = b"CGCAGAGTVVVVVVTACTGVVVVVVGAAT";
/// Minimum seal score. Perfect anchors score 22 (17 anchor bases * +2, minus the 12 UMI positions
/// that never match the V wildcards, * -1); 16 allows up to two anchor mismatches (each costs 3),
/// covering ONT error while cutting the spurious low-score tail.
const MIN_SEAL_SCORE: i32 = 16;
/// A seal starting past this many nt is a likely chimera (a second molecule's TSO): still emitted,
/// but counted so the residual is visible rather than silent.
const FAR_SEAL_START: usize = 100;
const MIN_READ_LEN: usize = 40;
const UMI_LEN: usize = 6;

/// Extract the two 6 nt 5' UMIs and the cDNA (after the seal) from a read, with the seal start (for
/// the chimera flag). None if the read is short, the seal scores below the floor, or the walk drifts
/// to a wrong-length UMI (an indel in or around it), whose key can't be trusted for dedup.
pub fn seal_extraction(read: &[u8]) -> Option<(Vec<u8>, Vec<u8>, Vec<u8>, usize)> {
    if read.len() < MIN_READ_LEN {
        return None;
    }
    let (xstart, xend, ops, score) = fit(read, TSO)?;
    if score < MIN_SEAL_SCORE {
        return None;
    }
    let txp = read[xend..].to_vec();

    let mut it = ops.iter();
    let mut offset = xstart;
    advance(&mut offset, 8, &mut it)?; // CGCAGAGT
    let s1 = offset;
    advance(&mut offset, 6, &mut it)?; // UMI1
    let umi1 = read[s1..offset].to_vec();
    advance(&mut offset, 5, &mut it)?; // TACTG
    let s2 = offset;
    advance(&mut offset, 6, &mut it)?; // UMI2
    let umi2 = read[s2..offset].to_vec();

    if umi1.len() != UMI_LEN || umi2.len() != UMI_LEN {
        return None;
    }
    Some((umi1, umi2, txp, xstart))
}

/// Per-run TSO counts.
#[derive(Default, Clone, Copy, Debug)]
pub struct Stats {
    pub total: u64,
    pub matched: u64,
    pub small: u64,
    pub no_seal: u64,
    pub far_seal: u64, // subset of matched: seal started past FAR_SEAL_START (possible chimera)
}

fn gz_writer(path: std::path::PathBuf) -> io::Result<GzEncoder<BufWriter<File>>> {
    Ok(GzEncoder::new(BufWriter::new(File::create(path)?), Compression::default()))
}

fn write_fasta<W: Write>(w: &mut W, id: &[u8], seq: &[u8]) -> io::Result<()> {
    w.write_all(b">")?;
    w.write_all(id)?;
    w.write_all(b"\n")?;
    w.write_all(seq)?;
    w.write_all(b"\n")
}

/// Extract the 5' UMI from each barcoded read, writing `>origid_CB_UMI_umi1_umi2` + trimmed cDNA to
/// the passed sink and unassigned reads unchanged to the failed sink.
pub fn run_tso(r1: &Path, out_dir: &Path) -> io::Result<Stats> {
    let passed = gz_writer(out_dir.join("passed.tso.nanopore.fa.gz"))?;
    let failed = gz_writer(out_dir.join("failed.tso.nanopore.fa.gz"))?;
    let mut reader = fastq::open_gz(r1)?;

    parallel::run(
        || {
            reader.next().map(|rec| {
                rec.map(|r| (fastq::read_name(r.id()).to_vec(), r.seq().to_vec()))
                    .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
            })
        },
        parallel::default_workers(),
        || (),
        |_: &mut (), (name, seq): (Vec<u8>, Vec<u8>)| {
            let e = seal_extraction(&seq);
            (name, seq, e)
        },
        move |rx| -> io::Result<Stats> {
            let mut passed = passed;
            let mut failed = failed;
            let mut stats = Stats::default();
            for (name, seq, e) in rx {
                stats.total += 1;
                if stats.total % 1_000_000 == 0 {
                    info!("tso {}M reads", stats.total / 1_000_000);
                }
                match e {
                    Some((umi1, umi2, txp, seal_start)) if seq.len() >= MIN_READ_LEN => {
                        let mut id = name;
                        id.push(b'_');
                        id.extend_from_slice(&umi1);
                        id.push(b'_');
                        id.extend_from_slice(&umi2);
                        write_fasta(&mut passed, &id, &txp)?;
                        stats.matched += 1;
                        stats.far_seal += (seal_start > FAR_SEAL_START) as u64;
                    }
                    _ if seq.len() < MIN_READ_LEN => {
                        write_fasta(&mut failed, &name, &seq)?;
                        stats.small += 1;
                    }
                    _ => {
                        write_fasta(&mut failed, &name, &seq)?;
                        stats.no_seal += 1;
                    }
                }
            }
            passed.finish()?;
            failed.finish()?;
            Ok(stats)
        },
    )
}

#[cfg(test)]
mod tests {
    use super::*;

    fn tso_read(leftover: &str, umi1: &str, umi2: &str, cdna: &str) -> String {
        format!("{leftover}CGCAGAGT{umi1}TACTG{umi2}GAAT{cdna}")
    }

    #[test]
    fn extracts_both_umis_and_trims_cdna() {
        let (leftover, umi1, umi2) = ("ACGTACGTAC", "AACCGG", "TTGGCA");
        let cdna = "GAGAGAGAGAGAGAGAGAGAGAGAGAGAGA"; // 30 nt, no anchor substrings
        let read = tso_read(leftover, umi1, umi2, cdna);
        let (g1, g2, txp, start) = seal_extraction(read.as_bytes()).expect("seal");
        assert_eq!(g1, umi1.as_bytes());
        assert_eq!(g2, umi2.as_bytes());
        assert_eq!(txp, cdna.as_bytes(), "cDNA trimmed to after the seal");
        assert_eq!(start, leftover.len(), "seal start = leftover length");
    }

    #[test]
    fn rejects_read_without_a_seal() {
        // No TSO anchor anywhere: the best fit scores below the floor.
        let read = "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT";
        assert!(seal_extraction(read.as_bytes()).is_none());
    }

    #[test]
    fn rejects_short_read() {
        assert!(seal_extraction(b"CGCAGAGTAACCGGTACTG").is_none());
    }
}
