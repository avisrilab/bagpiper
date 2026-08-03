//! Barcode-side seal: a Smith-Waterman rescue for reads the linker regex misses.
//!
//! The read is fit against the anchor pattern `CTCGA V6 CTC V6 CAT` (V = any base, so those columns
//! score as mismatches) with the shared [`crate::sw`] aligner, then the four barcode segments and the
//! UMI are read off fixed offsets around the anchor.

use crate::sw::{advance, fit};

/// Barcode-side anchor: CTCGA, a 6 nt segment, CTC, a 6 nt segment, CAT.
const SEAL: &[u8] = b"CTCGAVVVVVVCTCVVVVVVCAT";

/// Rescue barcode segments, UMI, and cDNA from a read via the barcode-side seal. Returns the four
/// segments in read order [bc1, bc2, bc3, bc4], the 12 nt UMI, and the cDNA prefix; None if the seal
/// aligns too near an end to carry the flanking UMI and terminal barcodes.
pub fn seal_extraction(read: &[u8]) -> Option<(Vec<Vec<u8>>, Vec<u8>, Vec<u8>)> {
    let (xstart, xend, ops, _score) = fit(read, SEAL)?;
    if xstart < 40 || xend + 8 > read.len() {
        return None;
    }

    let txp = read[0..xstart - 20].to_vec();
    let umi = read[xstart - 20..xstart - 8].to_vec();
    let mut barcodes = vec![read[xstart - 8..xstart].to_vec()]; // bc adjacent to CTCGA

    let mut it = ops.iter();
    let mut offset = xstart;
    advance(&mut offset, 5, &mut it)?; // CTCGA
    let seg = offset;
    advance(&mut offset, 6, &mut it)?; // segment
    barcodes.push(read[seg..offset].to_vec());
    advance(&mut offset, 3, &mut it)?; // CTC
    let seg = offset;
    advance(&mut offset, 6, &mut it)?; // segment
    barcodes.push(read[seg..offset].to_vec());
    advance(&mut offset, 3, &mut it)?; // CAT
    if offset + 8 > read.len() {
        return None;
    }
    barcodes.push(read[offset..offset + 8].to_vec()); // terminal barcode
    barcodes.reverse();

    Some((barcodes, umi, txp))
}

#[cfg(test)]
mod tests {
    use super::*;

    // Read layout the seal expects: [txp][umi:12][bcA:8] CTCGA [bcB:6] CTC [bcC:6] CAT [bcD:8][tail]
    fn seal_read(txp: &str, umi: &str, a: &str, b: &str, c: &str, d: &str, tail: &str) -> String {
        format!("{txp}{umi}{a}CTCGA{b}CTC{c}CAT{d}{tail}")
    }

    #[test]
    fn extracts_barcodes_umi_and_txp() {
        let txp = "GAGAGAGAGAGAGAGAGAGAGAGAGAGAGA"; // 30 nt, no anchor substrings
        let umi = "ACACACGTGTGT";
        let (a, b, c, d) = ("AAGGTTCC", "GGAACC", "TTGGAA", "CCTTGGAA");
        let read = seal_read(txp, umi, a, b, c, d, "TATATATATA");
        let (barcodes, got_umi, got_txp) = seal_extraction(read.as_bytes()).expect("extract");
        assert_eq!(got_txp, txp.as_bytes());
        assert_eq!(got_umi, umi.as_bytes());
        assert_eq!(
            barcodes,
            vec![d.as_bytes(), c.as_bytes(), b.as_bytes(), a.as_bytes()]
        );
    }

    #[test]
    fn rejects_seal_too_early() {
        let read = seal_read(
            "GAGAG",
            "ACACACGTGTGT",
            "AAGGTTCC",
            "GGAACC",
            "TTGGAA",
            "CCTTGG",
            "TA",
        );
        assert!(seal_extraction(read.as_bytes()).is_none());
    }

    #[test]
    fn rejects_no_room_for_last_barcode() {
        let read = format!(
            "{}{}{}CTCGA{}CTC{}CAT{}",
            "GAGAGAGAGAGAGAGAGAGAGAGAGAGAGA", "ACACACGTGTGT", "AAGGTTCC", "GGAACC", "TTGGAA", "CCT"
        );
        assert!(seal_extraction(read.as_bytes()).is_none());
    }
}
