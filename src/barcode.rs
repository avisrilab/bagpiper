//! Barcode assignment: per-read extraction of the cell barcode, UMI, and cDNA.
//!
//! Nanopore is a cascade: reverse-strand regex, forward-strand regex, then the SW seal on the read
//! and on its reverse complement. Illumina matches the reverse-strand regex on revcomp(R1) and takes
//! the cDNA from R2. This module returns the assignment; the I/O driver writes it and normalizes the
//! strand.

use std::io;
use std::path::{Path, PathBuf};

use log::info;
use regex::Regex;

use crate::chemistry::Chemistry;
use crate::fastq;
use crate::parallel;
use crate::seal;
use crate::seq::{revcomp, CellId};
use crate::whitelist::{is_ambiguous, Whitelist};

/// Outcome of assigning one read (pair). On the non-polyA strand (`is_polya` false) the driver
/// reverse-complements the UMI and cDNA so all output lands on one strand.
#[derive(Clone, PartialEq, Eq, Debug)]
pub enum Assign {
    Matched {
        cell: CellId,
        umi: Vec<u8>,
        cdna: Vec<u8>,
        is_polya: bool,
    },
    Small,
    NoRegex,
    MultiRegex,
    Mismatch,
    Ambiguous,
}

const MIN_READ_LEN: usize = 40;
const MIN_CDNA_LEN: usize = 50;

/// Per-stage result inside the nanopore cascade: assigned (stop), ambiguous (remember, keep trying),
/// or a miss (keep trying).
enum Stage {
    Matched(Assign),
    Ambiguous,
    Miss,
}

/// Resolve the four captured segments to a CellId, reverse-complementing the observed segments on
/// the forward strand (matching the reference). Err(true) = ambiguous, Err(false) = not found.
fn resolve(caps: &regex::Captures, wl: &Whitelist, is_forward: bool) -> Result<CellId, bool> {
    let orient = |b: &[u8]| if is_forward { revcomp(b) } else { b.to_vec() };
    let seg = |n: &str| orient(caps.name(n).unwrap().as_str().as_bytes());
    let owned = [seg("bc1"), seg("bc2"), seg("bc3"), seg("bc4")];
    let rows = wl.resolve([&owned[0], &owned[1], &owned[2], &owned[3]]);
    CellId::from_rows(rows).ok_or_else(|| is_ambiguous(&rows))
}

/// One regex stage: single match, cDNA length floor, whitelist resolution. `is_forward` sets both
/// the barcode orientation and the returned strand (reverse strand carries the polyA).
fn regex_stage(re: &Regex, wl: &Whitelist, read: &str, is_forward: bool) -> Stage {
    let mut it = re.captures_iter(read);
    let caps = match it.next() {
        Some(c) => c,
        None => return Stage::Miss,
    };
    if it.next().is_some() {
        return Stage::Miss; // more than one placement
    }
    let cdna = caps.name("seq").unwrap().as_str().as_bytes();
    if cdna.len() < MIN_CDNA_LEN {
        return Stage::Miss;
    }
    match resolve(&caps, wl, is_forward) {
        Ok(cell) => Stage::Matched(Assign::Matched {
            cell,
            umi: caps.name("umi").unwrap().as_str().as_bytes().to_vec(),
            cdna: cdna.to_vec(),
            is_polya: !is_forward,
        }),
        Err(true) => Stage::Ambiguous,
        Err(false) => Stage::Miss,
    }
}

/// One SW-seal stage on `read` (already in the orientation to search). The seal barcodes are used in
/// stored orientation (no revcomp), matching the reference `get_long_cb_from_string`.
fn seal_stage(read: &[u8], wl: &Whitelist) -> Stage {
    let (bc, umi, cdna) = match seal::seal_extraction(read) {
        Some(x) => x,
        None => return Stage::Miss,
    };
    let rows = wl.resolve([&bc[0], &bc[1], &bc[2], &bc[3]]);
    match CellId::from_rows(rows) {
        Some(cell) => Stage::Matched(Assign::Matched {
            cell,
            umi,
            cdna,
            is_polya: true,
        }),
        None => {
            if is_ambiguous(&rows) {
                Stage::Ambiguous
            } else {
                Stage::Miss
            }
        }
    }
}

/// Assign one nanopore read: reverse regex, forward regex, seal(read), seal(revcomp).
pub fn assign_nanopore(read: &[u8], rev_re: &Regex, fwd_re: &Regex, wl: &Whitelist) -> Assign {
    if read.len() < MIN_READ_LEN {
        return Assign::Small;
    }
    let s = std::str::from_utf8(read).expect("read is ACGT");
    let mut ambiguous = false;

    for (re, is_forward) in [(rev_re, false), (fwd_re, true)] {
        match regex_stage(re, wl, s, is_forward) {
            Stage::Matched(a) => return a,
            Stage::Ambiguous => ambiguous = true,
            Stage::Miss => {}
        }
    }

    let rc = revcomp(read);
    for cand in [read, &rc[..]] {
        match seal_stage(cand, wl) {
            Stage::Matched(a) => return a,
            Stage::Ambiguous => ambiguous = true,
            Stage::Miss => {}
        }
    }

    if ambiguous {
        Assign::Ambiguous
    } else {
        Assign::Mismatch
    }
}

/// Assign one illumina pair: match the reverse-strand regex on revcomp(R1), correct the barcode, and
/// carry the cDNA from R2. Illumina output is written as-is (`is_polya` true).
pub fn assign_illumina(r1: &[u8], r2: &[u8], rev_re: &Regex, wl: &Whitelist) -> Assign {
    let r1rc = revcomp(r1);
    let s = std::str::from_utf8(&r1rc).expect("R1 is ACGT");
    let mut it = rev_re.captures_iter(s);
    let caps = match it.next() {
        Some(c) => c,
        None => return Assign::NoRegex,
    };
    if it.next().is_some() {
        return Assign::MultiRegex;
    }
    match resolve(&caps, wl, true) {
        Ok(cell) => Assign::Matched {
            cell,
            umi: caps.name("umi").unwrap().as_str().as_bytes().to_vec(),
            cdna: r2.to_vec(),
            is_polya: true,
        },
        Err(true) => Assign::Ambiguous,
        Err(false) => Assign::Mismatch,
    }
}

/// Per-run assignment counts.
#[derive(Default, Clone, Copy, Debug)]
pub struct Stats {
    pub total: u64,
    pub matched: u64,
    pub small: u64,
    pub no_regex: u64,
    pub multi: u64,
    pub ambiguous: u64,
    pub mismatch: u64,
}

fn record_id(name: &[u8], cell: CellId, umi: &[u8]) -> Vec<u8> {
    let mut id = name.to_vec();
    id.push(b'_');
    id.extend_from_slice(cell.render().as_bytes());
    id.push(b'_');
    id.extend_from_slice(umi);
    id
}

/// Nanopore: assign each read via the cascade, writing `>origid_CB_UMI` + cDNA to the passed sink
/// (UMI and cDNA reverse-complemented on the non-polyA strand) and unassigned reads unchanged to the
/// failed sink.
pub fn run_nanopore(
    r1: &[PathBuf],
    wl_path: &Path,
    out_dir: &Path,
    workers: usize,
) -> io::Result<Stats> {
    let wl = Whitelist::from_csv(wl_path)?;
    let rev = Chemistry::PipV4.barcode_regex(true);
    let fwd = Chemistry::PipV4.barcode_regex(false);
    let passed = fastq::gz_writer(out_dir.join("passed.bcd.nanopore.fa.gz"))?;
    let failed = fastq::gz_writer(out_dir.join("failed.bcd.nanopore.fa.gz"))?;
    let mut reader = fastq::MultiReader::open(r1);

    parallel::run(
        || {
            reader
                .next_seq()
                .map(|res| res.map(|(id, seq)| (fastq::read_name(&id).to_vec(), seq)))
        },
        workers,
        || (),
        |_: &mut (), (name, seq): (Vec<u8>, Vec<u8>)| {
            let a = assign_nanopore(&seq, &rev, &fwd, &wl);
            (name, seq, a)
        },
        move |rx| -> io::Result<Stats> {
            let mut passed = passed;
            let mut failed = failed;
            let mut stats = Stats::default();
            for (name, seq, a) in rx {
                stats.total += 1;
                if stats.total % 1_000_000 == 0 {
                    info!("barcoded {}M reads", stats.total / 1_000_000);
                }
                match a {
                    Assign::Matched {
                        cell,
                        umi,
                        cdna,
                        is_polya,
                    } => {
                        if is_polya {
                            fastq::write_fasta(&mut passed, &record_id(&name, cell, &umi), &cdna)?;
                        } else {
                            let (umi, cdna) = (revcomp(&umi), revcomp(&cdna));
                            fastq::write_fasta(&mut passed, &record_id(&name, cell, &umi), &cdna)?;
                        }
                        stats.matched += 1;
                    }
                    Assign::Small => {
                        fastq::write_fasta(&mut failed, &name, &seq)?;
                        stats.small += 1;
                    }
                    Assign::Ambiguous => {
                        fastq::write_fasta(&mut failed, &name, &seq)?;
                        stats.ambiguous += 1;
                    }
                    _ => {
                        fastq::write_fasta(&mut failed, &name, &seq)?;
                        stats.mismatch += 1;
                    }
                }
            }
            passed.finish()?;
            failed.finish()?;
            Ok(stats)
        },
    )
}

/// Illumina: match the reverse regex on revcomp(R1), pair the barcode with R2 as cDNA, and write
/// `>origid_CB_UMI` + R2 for assigned pairs. Unassigned pairs are only counted.
pub fn run_illumina(
    r1: &Path,
    r2: &Path,
    wl_path: &Path,
    out_dir: &Path,
    workers: usize,
) -> io::Result<Stats> {
    let wl = Whitelist::from_csv(wl_path)?;
    let rev = Chemistry::PipV4.barcode_regex(true);
    let out = fastq::gz_writer(out_dir.join("read.bcd.illumina.fa.gz"))?;
    let mut r1r = fastq::open_gz(r1)?;
    let mut r2r = fastq::open_gz(r2)?;

    parallel::run(
        || match (r1r.next(), r2r.next()) {
            (Some(a), Some(b)) => Some(
                a.and_then(|ra| b.map(|rb| (ra, rb)))
                    .map(|(ra, rb)| {
                        (
                            fastq::read_name(ra.id()).to_vec(),
                            ra.seq().to_vec(),
                            rb.seq().to_vec(),
                        )
                    })
                    .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e)),
            ),
            (None, None) => None,
            _ => Some(Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "R1/R2 record count mismatch",
            ))),
        },
        workers,
        || (),
        |_: &mut (), (name, r1seq, r2seq): (Vec<u8>, Vec<u8>, Vec<u8>)| {
            let a = assign_illumina(&r1seq, &r2seq, &rev, &wl);
            (name, a)
        },
        move |rx| -> io::Result<Stats> {
            let mut out = out;
            let mut stats = Stats::default();
            for (name, a) in rx {
                stats.total += 1;
                if stats.total % 1_000_000 == 0 {
                    info!("barcoded {}M reads", stats.total / 1_000_000);
                }
                match a {
                    Assign::Matched {
                        cell, umi, cdna, ..
                    } => {
                        fastq::write_fasta(&mut out, &record_id(&name, cell, &umi), &cdna)?;
                        stats.matched += 1;
                    }
                    Assign::NoRegex => stats.no_regex += 1,
                    Assign::MultiRegex => stats.multi += 1,
                    Assign::Ambiguous => stats.ambiguous += 1,
                    _ => stats.mismatch += 1,
                }
            }
            out.finish()?;
            Ok(stats)
        },
    )
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::chemistry::Chemistry;
    use std::io::Write;

    fn rc(s: &str) -> String {
        String::from_utf8(revcomp(s.as_bytes())).unwrap()
    }

    fn write_wl(rows: &[[&str; 4]]) -> std::path::PathBuf {
        let p = std::env::temp_dir().join(format!("bp_bc_wl_{}.csv", std::process::id()));
        let mut f = std::fs::File::create(&p).unwrap();
        writeln!(f, "barcode_4,barcode_3,barcode_2,barcode_1").unwrap();
        for r in rows {
            writeln!(f, "{},{},{},{}", r[0], r[1], r[2], r[3]).unwrap();
        }
        p
    }

    /// A reverse-strand nanopore read for row-0 of the test whitelist: the observed segments equal
    /// the stored (revcomp'd) whitelist keys, so no forward revcomp is needed to match.
    fn reverse_read(cdna: &str, umi: &str, bc4: &str, bc3: &str, bc2: &str, bc1: &str) -> String {
        format!("{cdna}{umi}{bc4}CTCGA{bc3}CTC{bc2}CAT{bc1}AA")
    }

    #[test]
    fn nanopore_reverse_regex_assigns_row0() {
        let p = write_wl(&[
            ["AAAAAAAA", "AAAAAA", "CCCCCC", "CCCCCCCC"], // row 0
            ["GGGGGGGG", "GGGGGG", "TTTTTT", "TTTTTTTT"], // row 1
        ]);
        let wl = Whitelist::from_csv(&p).unwrap();
        let rev = Chemistry::PipV4.barcode_regex(true);
        let fwd = Chemistry::PipV4.barcode_regex(false);

        // Observed reverse-strand segments = revcomp of the whitelist fields (stored orientation).
        let umi = "ACACACGTGTGT";
        let read = reverse_read(
            &"A".repeat(55),
            umi,
            &rc("AAAAAAAA"), // bc4 -> TTTTTTTT
            &rc("AAAAAA"),   // bc3 -> TTTTTT
            &rc("CCCCCC"),   // bc2 -> GGGGGG
            &rc("CCCCCCCC"), // bc1 -> GGGGGGGG
        );

        match assign_nanopore(read.as_bytes(), &rev, &fwd, &wl) {
            Assign::Matched {
                cell,
                umi: u,
                cdna,
                is_polya,
            } => {
                assert_eq!(
                    cell.render(),
                    "A".repeat(16),
                    "row 0 in every segment -> all-A barcode"
                );
                assert_eq!(u, umi.as_bytes());
                assert_eq!(cdna, "A".repeat(55).as_bytes());
                assert!(is_polya, "reverse strand carries the polyA");
            }
            other => panic!("expected a match, got {other:?}"),
        }
        let _ = std::fs::remove_file(&p);
    }

    #[test]
    fn nanopore_off_whitelist_is_mismatch() {
        let p = write_wl(&[["AAAAAAAA", "AAAAAA", "CCCCCC", "CCCCCCCC"]]);
        let wl = Whitelist::from_csv(&p).unwrap();
        let rev = Chemistry::PipV4.barcode_regex(true);
        let fwd = Chemistry::PipV4.barcode_regex(false);
        // bc1 = ACGTACGT is neither exact nor edit-1 of the single row.
        let read = reverse_read(
            &"A".repeat(55),
            "ACACACGTGTGT",
            "TTTTTTTT",
            "TTTTTT",
            "GGGGGG",
            "ACGTACGT",
        );
        assert_eq!(
            assign_nanopore(read.as_bytes(), &rev, &fwd, &wl),
            Assign::Mismatch
        );
        let _ = std::fs::remove_file(&p);
    }

    #[test]
    fn nanopore_short_read_is_small() {
        let p = write_wl(&[["AAAAAAAA", "AAAAAA", "CCCCCC", "CCCCCCCC"]]);
        let wl = Whitelist::from_csv(&p).unwrap();
        let rev = Chemistry::PipV4.barcode_regex(true);
        let fwd = Chemistry::PipV4.barcode_regex(false);
        assert_eq!(
            assign_nanopore(b"ACGTACGTACGT", &rev, &fwd, &wl),
            Assign::Small
        );
        let _ = std::fs::remove_file(&p);
    }
}
