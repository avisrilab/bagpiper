//! Barcode-side seal: a Smith-Waterman rescue for reads the linker regex misses.
//!
//! The read is fit against the anchor pattern `CTCGA V6 CTC V6 CAT` (V = any base, so those columns
//! score as mismatches), placing the whole pattern into a substring of the read (pattern global,
//! read free-clipped both ends) with affine gaps. The alignment walk then reads the four barcode
//! segments and the UMI off fixed offsets around the anchor. Scoring, gap model, and tie-break order
//! match rust-bio's `Aligner::custom` for these clip parameters (verified by `tests/seal_diff.rs`).

const MATCH: i32 = 2;
const MISMATCH: i32 = -1;
const GAP_OPEN: i32 = -1;
const GAP_EXTEND: i32 = -2;
const NEG_INF: i32 = i32::MIN / 4;

/// Barcode-side anchor: CTCGA, a 6 nt segment, CTC, a 6 nt segment, CAT.
const SEAL: &[u8] = b"CTCGAVVVVVVCTCVVVVVVCAT";

#[derive(Clone, Copy, PartialEq, Eq, Debug)]
enum Op {
    Match,
    Subst,
    Ins, // gap in the pattern: consume a read base
    Del, // gap in the read: consume a pattern base
}

/// Fit `pat` entirely into a substring of `read` (pattern global, read free-clipped both ends),
/// affine gaps. Returns (xstart, xend, ops from xstart to xend). None if the pattern cannot fit.
fn fit(read: &[u8], pat: &[u8]) -> Option<(usize, usize, Vec<Op>)> {
    let (m, n) = (read.len(), pat.len());
    if n == 0 || m < n {
        return None;
    }
    // S: best ending aligned/clipped; ix: gap in pattern (consume read); dy: gap in read (consume pattern).
    let mut s = vec![vec![NEG_INF; n + 1]; m + 1];
    let mut ix = vec![vec![NEG_INF; n + 1]; m + 1];
    let mut dy = vec![vec![NEG_INF; n + 1]; m + 1];
    for row in s.iter_mut() {
        row[0] = 0; // free read-prefix clip
    }
    for j in 1..=n {
        let g = GAP_OPEN + GAP_EXTEND * j as i32;
        s[0][j] = g;
        dy[0][j] = g;
    }
    for i in 1..=m {
        for j in 1..=n {
            ix[i][j] = (ix[i - 1][j] + GAP_EXTEND).max(s[i - 1][j] + GAP_OPEN + GAP_EXTEND);
            dy[i][j] = (dy[i][j - 1] + GAP_EXTEND).max(s[i][j - 1] + GAP_OPEN + GAP_EXTEND);
            let diag = s[i - 1][j - 1] + if read[i - 1] == pat[j - 1] { MATCH } else { MISMATCH };
            s[i][j] = diag.max(ix[i][j]).max(dy[i][j]);
        }
    }

    // Best full-pattern alignment ends at the smallest row i maximizing s[i][n] (free read-suffix clip).
    let mut xend = 0;
    let mut best = NEG_INF;
    for i in 0..=m {
        if s[i][n] > best {
            best = s[i][n];
            xend = i;
        }
    }

    // Traceback, preferring match/subst > ins > del, and gap-open over gap-extend on ties (bio order).
    let mut ops = Vec::new();
    let (mut i, mut j) = (xend, n);
    #[derive(PartialEq)]
    enum Layer {
        S,
        I,
        D,
    }
    let mut layer = Layer::S;
    while j > 0 {
        match layer {
            Layer::S => {
                let diag = if i > 0 {
                    s[i - 1][j - 1] + if read[i - 1] == pat[j - 1] { MATCH } else { MISMATCH }
                } else {
                    NEG_INF
                };
                if i > 0 && s[i][j] == diag {
                    ops.push(if read[i - 1] == pat[j - 1] { Op::Match } else { Op::Subst });
                    i -= 1;
                    j -= 1;
                } else if s[i][j] == ix[i][j] {
                    layer = Layer::I;
                } else {
                    layer = Layer::D;
                }
            }
            Layer::I => {
                ops.push(Op::Ins);
                if ix[i][j] == s[i - 1][j] + GAP_OPEN + GAP_EXTEND {
                    layer = Layer::S; // opened
                }
                i -= 1;
            }
            Layer::D => {
                ops.push(Op::Del);
                if dy[i][j] == s[i][j - 1] + GAP_OPEN + GAP_EXTEND {
                    layer = Layer::S; // opened
                }
                j -= 1;
            }
        }
    }
    let xstart = i;
    ops.reverse();
    Some((xstart, xend, ops))
}

/// Advance `offset` (a read index) past `xlim` pattern columns, following the alignment ops. Mirrors
/// the reference walk: a deletion consumes a pattern column without a read base, an insertion a read
/// base without a column.
fn advance(offset: &mut usize, xlim: usize, ops: &mut std::slice::Iter<Op>) -> Option<()> {
    let mut check = 0;
    loop {
        match ops.next()? {
            Op::Match | Op::Subst => check += 1,
            Op::Del => {
                *offset = offset.checked_sub(1)?;
                check += 1;
            }
            Op::Ins => *offset += 1,
        }
        if check == xlim {
            *offset += xlim;
            return Some(());
        }
    }
}

/// Rescue barcode segments, UMI, and cDNA from a read via the barcode-side seal. Returns the four
/// segments in read order [bc1, bc2, bc3, bc4], the 12 nt UMI, and the cDNA prefix; None if the seal
/// aligns too near an end to carry the flanking UMI and terminal barcodes.
pub fn seal_extraction(read: &[u8]) -> Option<(Vec<Vec<u8>>, Vec<u8>, Vec<u8>)> {
    let (xstart, xend, ops) = fit(read, SEAL)?;
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
        let read = seal_read("GAGAG", "ACACACGTGTGT", "AAGGTTCC", "GGAACC", "TTGGAA", "CCTTGG", "TA");
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
