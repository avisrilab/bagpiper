//! Barcode assignment: per-read extraction of the cell barcode, UMI, and cDNA.
//!
//! Nanopore is a cascade: reverse-strand regex, forward-strand regex, then the SW seal on the read
//! and on its reverse complement. Illumina matches the reverse-strand regex on revcomp(R1) and takes
//! the cDNA from R2. This module returns the assignment; the I/O driver writes it and normalizes the
//! strand.

use regex::Regex;

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
            Assign::Matched { cell, umi: u, cdna, is_polya } => {
                assert_eq!(cell.render(), "A".repeat(16), "row 0 in every segment -> all-A barcode");
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
        let read = reverse_read(&"A".repeat(55), "ACACACGTGTGT", "TTTTTTTT", "TTTTTT", "GGGGGG", "ACGTACGT");
        assert_eq!(assign_nanopore(read.as_bytes(), &rev, &fwd, &wl), Assign::Mismatch);
        let _ = std::fs::remove_file(&p);
    }

    #[test]
    fn nanopore_short_read_is_small() {
        let p = write_wl(&[["AAAAAAAA", "AAAAAA", "CCCCCC", "CCCCCCCC"]]);
        let wl = Whitelist::from_csv(&p).unwrap();
        let rev = Chemistry::PipV4.barcode_regex(true);
        let fwd = Chemistry::PipV4.barcode_regex(false);
        assert_eq!(assign_nanopore(b"ACGTACGTACGT", &rev, &fwd, &wl), Assign::Small);
        let _ = std::fs::remove_file(&p);
    }
}
