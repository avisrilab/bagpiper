//! Barcode chemistry: each method's read layout and how a barcode read is parsed, behind a seam
//! that new methods join by adding a [`Chemistry`] variant.
//!
//! Only PIP-seq V4 exists today. V5, 10x, and split-pool are expected to add their own variant with
//! their own extraction (10x slices by fixed position rather than matching linkers), so the barcode
//! stage dispatches on the chemistry instead of hard-coding one layout.

use regex::Regex;

use crate::seq::revcomp;

/// A supported barcode chemistry. Add a variant per method; keep speculative methods out until
/// they are implemented.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum Chemistry {
    /// PIP-seq V4: four cell-barcode segments (8-6-6-8 nt) and a 12 nt UMI joined by fixed linkers.
    PipV4,
}

impl Chemistry {
    /// The barcode-read regex with named captures. Forward strand yields
    /// `NL bc1 bc2 bc3 bc4 umi seq`; reverse yields `seq umi bc4 bc3 bc2 bc1 NR`. The cDNA flank is
    /// always `seq`; linkers are matched but not captured.
    pub fn barcode_regex(&self, reverse: bool) -> Regex {
        match self {
            Chemistry::PipV4 => pipseq_v4_regex(reverse),
        }
    }
}

const V4_BC_LENS: [usize; 4] = [8, 6, 6, 8];
const V4_UMI_LEN: usize = 12;
const V4_LINKERS: [&str; 3] = ["ATG", "GAG", "TCGAG"]; // bc1-bc2, bc2-bc3, bc3-bc4

fn pipseq_v4_regex(reverse: bool) -> Regex {
    let nt = |n: usize| format!("[ACGT]{{{}}}", n);
    let (b1, b2, b3, b4, umi) = (
        nt(V4_BC_LENS[0]),
        nt(V4_BC_LENS[1]),
        nt(V4_BC_LENS[2]),
        nt(V4_BC_LENS[3]),
        nt(V4_UMI_LEN),
    );
    let rc = |s: &str| String::from_utf8(revcomp(s.as_bytes())).unwrap();
    let pattern = if reverse {
        format!(
            "(?<seq>[ACGT]*)(?<umi>{umi})(?<bc4>{b4})(?:{l3})(?<bc3>{b3})(?:{l2})(?<bc2>{b2})(?:{l1})(?<bc1>{b1})(?<NR>[ACGT]*)",
            l1 = rc(V4_LINKERS[0]),
            l2 = rc(V4_LINKERS[1]),
            l3 = rc(V4_LINKERS[2]),
        )
    } else {
        format!(
            "(?<NL>[ACGT]*)(?<bc1>{b1})(?:{l1})(?<bc2>{b2})(?:{l2})(?<bc3>{b3})(?:{l3})(?<bc4>{b4})(?<umi>{umi})(?<seq>[ACGT]*)",
            l1 = V4_LINKERS[0],
            l2 = V4_LINKERS[1],
            l3 = V4_LINKERS[2],
        )
    };
    Regex::new(&pattern).unwrap()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn v4_forward_captures_barcodes_umi_and_cdna() {
        let re = Chemistry::PipV4.barcode_regex(false);
        // NL, then bc1 ATG bc2 GAG bc3 TCGAG bc4 umi, then cDNA. Segments avoid the linker motifs.
        let read = "GG\
            CACACACA\
            ATG\
            TCTCTC\
            GAG\
            AGAGAG\
            TCGAG\
            CTCTCTCT\
            ACACACGTGTGT\
            TTTTCCCCAAAAGGGGTTTTCCCC";
        let caps = re.captures(read).expect("forward read should match");
        assert_eq!(&caps["bc1"], "CACACACA");
        assert_eq!(&caps["bc2"], "TCTCTC");
        assert_eq!(&caps["bc3"], "AGAGAG");
        assert_eq!(&caps["bc4"], "CTCTCTCT");
        assert_eq!(&caps["umi"], "ACACACGTGTGT");
        assert_eq!(&caps["seq"], "TTTTCCCCAAAAGGGGTTTTCCCC");
    }

    #[test]
    fn v4_reverse_captures_with_revcomp_linkers() {
        let re = Chemistry::PipV4.barcode_regex(true);
        // seq, umi, bc4, CTCGA, bc3, CTC, bc2, CAT, bc1, NR (reverse-strand linkers = revcomp).
        let read = "TTTTCCCCAAAAGGGGTTTTCCCC\
            ACACACGTGTGT\
            CTCTCTCT\
            CTCGA\
            AGAGAG\
            CTC\
            TCTCTC\
            CAT\
            CACACACA\
            GG";
        let caps = re.captures(read).expect("reverse read should match");
        assert_eq!(&caps["bc4"], "CTCTCTCT");
        assert_eq!(&caps["bc3"], "AGAGAG");
        assert_eq!(&caps["bc2"], "TCTCTC");
        assert_eq!(&caps["bc1"], "CACACACA");
        assert_eq!(&caps["umi"], "ACACACGTGTGT");
    }
}
