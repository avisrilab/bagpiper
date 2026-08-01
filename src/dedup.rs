//! Exact PCR deduplication.

use crate::eqclass::Molecule;

/// Collapse exact PCR duplicates in place: sort by packed (cell, UMI, transcript set) and drop
/// adjacent equals. Molecules equal in all three fields are one PCR-duplicate group. The surviving
/// set equals the reference exact-dedup; order is by packed key.
pub fn exact(molecules: &mut Vec<Molecule>) {
    molecules.sort_unstable();
    molecules.dedup();
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::seq::{CellId, Umi};

    fn mol(cb: &str, umi: &str, txps: &[u32]) -> Molecule {
        Molecule {
            cell: CellId::from_ascii(cb.as_bytes()).unwrap(),
            umi: Umi::from_ascii(umi.as_bytes()).unwrap(),
            txps: txps.to_vec(),
        }
    }

    #[test]
    fn collapses_identical_keeps_distinct() {
        // Three identical molecules collapse to one; a distinct cell and a distinct UMI survive.
        let mut m = vec![
            mol("AAAAAAAAAAAAAAAA", "ACGTACGTACGT", &[0]),
            mol("AAAAAAAAAAAAAAAA", "ACGTACGTACGT", &[0]),
            mol("AAAAAAAAAAAAAAAA", "ACGTACGTACGT", &[0]),
            mol("CAAAAAAAAAAAAAAA", "ACGTACGTACGT", &[0]),
            mol("AAAAAAAAAAAAAAAA", "TTTTTTTTTTTT", &[1]),
        ];
        exact(&mut m);
        assert_eq!(m.len(), 3);
    }

    #[test]
    fn same_cell_umi_different_txps_survive() {
        // Identical cell and UMI but a different transcript set is a distinct molecule.
        let mut m = vec![
            mol("AAAAAAAAAAAAAAAA", "ACGTACGTACGT", &[0]),
            mol("AAAAAAAAAAAAAAAA", "ACGTACGTACGT", &[1]),
        ];
        exact(&mut m);
        assert_eq!(m.len(), 2);
    }
}
