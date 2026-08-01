//! Packed sequence primitives and molecular-key types.

/// 2-bit base code: A=0, C=1, G=2, T=3. None for any other byte.
#[inline]
pub fn encode_base(b: u8) -> Option<u8> {
    match b {
        b'A' => Some(0),
        b'C' => Some(1),
        b'G' => Some(2),
        b'T' => Some(3),
        _ => None,
    }
}

/// Inverse of `encode_base` on the low 2 bits.
#[inline]
pub fn decode_base(code: u8) -> u8 {
    [b'A', b'C', b'G', b'T'][(code & 3) as usize]
}

/// Render the low `2*len` bits of `value` as `len` bases, most-significant base first.
pub fn render_2bit(value: u64, len: usize) -> String {
    let bytes: Vec<u8> = (0..len)
        .map(|i| decode_base((value >> ((len - 1 - i) * 2)) as u8))
        .collect();
    String::from_utf8(bytes).unwrap()
}

/// Cell barcode identity: the base-96 pack of four whitelist row indices (segment 0 least
/// significant), rendered to a 16 nt string for `barcodes.tsv`.
#[derive(Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash, Debug)]
pub struct CellId(pub u64);

impl CellId {
    const RENDER_LEN: usize = 16;
    const N_ROWS: u64 = 96;

    /// Pack four whitelist row indices. None if any index is not a valid row (0..96), which covers
    /// the not-found and ambiguous sentinels.
    pub fn from_rows(rows: [u64; 4]) -> Option<CellId> {
        let mut value = 0u64;
        let mut place = 1u64;
        for row in rows {
            if row >= Self::N_ROWS {
                return None;
            }
            value += row * place;
            place *= Self::N_ROWS;
        }
        Some(CellId(value))
    }

    /// The `barcodes.tsv` string.
    pub fn render(&self) -> String {
        render_2bit(self.0, Self::RENDER_LEN)
    }
}

/// Packed UMI, 2 bits per base, up to 32 nt. Upstream extraction emits ACGT only; a non-ACGT base
/// yields None and the caller drops the read. Length is part of identity.
#[derive(Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash, Debug)]
pub struct Umi {
    bits: u64,
    len: u8,
}

impl Umi {
    pub fn from_ascii(seq: &[u8]) -> Option<Umi> {
        if seq.len() > 32 {
            return None;
        }
        let mut bits = 0u64;
        for &b in seq {
            bits = (bits << 2) | encode_base(b)? as u64;
        }
        Some(Umi {
            bits,
            len: seq.len() as u8,
        })
    }

    pub fn len(&self) -> usize {
        self.len as usize
    }

    pub fn is_empty(&self) -> bool {
        self.len == 0
    }

    pub fn render(&self) -> String {
        render_2bit(self.bits, self.len as usize)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn render_2bit_is_msb_first() {
        assert_eq!(render_2bit(0, 16), "A".repeat(16));
        assert_eq!(render_2bit(1, 16), format!("{}C", "A".repeat(15)));
        assert_eq!(render_2bit(2, 16), format!("{}G", "A".repeat(15)));
        assert_eq!(render_2bit(3, 16), format!("{}T", "A".repeat(15)));
        assert_eq!(render_2bit(3, 4), "AAAT");
    }

    #[test]
    fn cellid_matches_base96_pack() {
        // segment 0 least significant; equals the reference short-CB encoding.
        assert_eq!(CellId::from_rows([0, 0, 0, 0]).unwrap().render(), "A".repeat(16));
        assert_eq!(
            CellId::from_rows([1, 0, 0, 0]).unwrap().render(),
            format!("{}C", "A".repeat(15))
        );
        assert_eq!(CellId::from_rows([0, 1, 0, 0]).unwrap(), CellId(96));
    }

    #[test]
    fn cellid_rejects_invalid_row() {
        assert!(CellId::from_rows([96, 0, 0, 0]).is_none());
        assert!(CellId::from_rows([u64::MAX, 0, 0, 0]).is_none());
        assert!(CellId::from_rows([0, 0, 0, u64::MAX - 1]).is_none());
    }

    #[test]
    fn umi_round_trips_and_preserves_length() {
        let u = Umi::from_ascii(b"ACGTACGTACGT").unwrap();
        assert_eq!(u.len(), 12);
        assert_eq!(u.render(), "ACGTACGTACGT");
        assert_ne!(
            Umi::from_ascii(b"AC").unwrap(),
            Umi::from_ascii(b"AAC").unwrap()
        );
    }

    #[test]
    fn umi_rejects_non_acgt() {
        assert!(Umi::from_ascii(b"ACGTN").is_none());
    }
}
