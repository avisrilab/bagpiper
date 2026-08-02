//! Cell-barcode whitelist: exact and edit-1 correction of the four barcode segments.
//!
//! Each whitelist segment is stored reverse-complemented, with its full edit-1 neighborhood in a
//! secondary table. A neighbor shared by two rows is marked ambiguous rather than assigned. Lookup
//! returns the four whitelist rows, which [`crate::seq::CellId::from_rows`] packs into the barcode.

use std::collections::HashMap;
use std::fs::File;
use std::io::{self, BufRead, BufReader};
use std::path::Path;

use crate::seq::revcomp;

/// A segment observed nowhere in the whitelist (not exact, not edit-1).
pub const NOT_FOUND: u64 = u64::MAX;
/// A segment edit-1 from two different rows, so not correctable unambiguously.
pub const AMBIGUOUS: u64 = u64::MAX - 1;

/// True if any resolved segment was left ambiguous (distinguishes ambiguity from a plain miss).
pub fn is_ambiguous(rows: &[u64]) -> bool {
    rows.iter().any(|&r| r == AMBIGUOUS)
}

/// Four-segment PIP-seq whitelist. `primary[i]` maps the reverse-complemented segment i to its row;
/// `secondary[i]` maps each edit-1 neighbor to its row (or [`AMBIGUOUS`]).
pub struct Whitelist {
    primary: [HashMap<Vec<u8>, u64>; 4],
    secondary: [HashMap<Vec<u8>, u64>; 4],
}

impl Whitelist {
    /// Load from a CSV whose columns are `barcode_4,barcode_3,barcode_2,barcode_1` (one header line,
    /// then one row per whitelist entry). Column i (0=bc4) is stored under segment index `3 - i`.
    pub fn from_csv<P: AsRef<Path>>(path: P) -> io::Result<Whitelist> {
        let mut primary: Vec<HashMap<Vec<u8>, u64>> = vec![HashMap::new(); 4];
        let mut secondary: Vec<HashMap<Vec<u8>, u64>> = vec![HashMap::new(); 4];

        let mut lines = BufReader::new(File::open(path)?).lines();
        lines.next(); // header

        let mut row: u64 = 0;
        for line in lines {
            let line = line?;
            for (i, field) in line.split(',').enumerate() {
                let seg = 3 - i;
                let key = revcomp(field.as_bytes());
                primary[seg].insert(key.clone(), row);
                for neighbor in edit1_combinations(&key) {
                    secondary[seg]
                        .entry(neighbor)
                        .and_modify(|existing| {
                            if *existing != row {
                                *existing = AMBIGUOUS;
                            }
                        })
                        .or_insert(row);
                }
            }
            row += 1;
        }

        Ok(Whitelist {
            primary: primary.try_into().unwrap(),
            secondary: secondary.try_into().unwrap(),
        })
    }

    /// Resolve four observed segments (in stored, i.e. reverse-complemented, orientation) to their
    /// whitelist rows. Exact match wins; otherwise the edit-1 table; otherwise [`NOT_FOUND`].
    pub fn resolve(&self, segments: [&[u8]; 4]) -> [u64; 4] {
        let mut rows = [NOT_FOUND; 4];
        for i in 0..4 {
            rows[i] = match self.primary[i].get(segments[i]) {
                Some(&r) => r,
                None => *self.secondary[i].get(segments[i]).unwrap_or(&NOT_FOUND),
            };
        }
        rows
    }
}

/// All strings within edit distance 1 of `s` (deletions, insertions, substitutions), plus the
/// length-restored variants the reference generates so an indel'd read still matches. Duplicates are
/// harmless: the secondary table folds them under the ambiguity guard.
fn edit1_combinations(s: &[u8]) -> Vec<Vec<u8>> {
    const NTS: [u8; 4] = [b'A', b'T', b'C', b'G'];
    let mut out = Vec::new();
    let n = s.len();

    // Deletion, and the deletion re-padded at either end.
    for i in 0..n {
        let mut d = s.to_vec();
        d.remove(i);
        out.push(d.clone());
        for &c in &NTS {
            let mut x = d.clone();
            x.push(c);
            out.push(x);
        }
        for &c in &NTS {
            let mut x = d.clone();
            x.insert(0, c);
            out.push(x);
        }
    }

    // Insertion, and the insertion trimmed at either end.
    for i in 0..=n {
        for &c in &NTS {
            let mut ins = s.to_vec();
            ins.insert(i, c);
            out.push(ins.clone());
            let mut front = ins.clone();
            front.remove(0);
            out.push(front);
            let mut back = ins;
            back.remove(back.len() - 1);
            out.push(back);
        }
    }

    // Substitution.
    for i in 0..n {
        for &c in &NTS {
            let mut x = s.to_vec();
            x[i] = c;
            out.push(x);
        }
    }

    out
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::seq::{revcomp, CellId};
    use std::io::Write;

    fn write_wl(tag: &str, rows: &[[&str; 4]]) -> std::path::PathBuf {
        let p = std::env::temp_dir().join(format!("bp_wl_{}_{}.csv", std::process::id(), tag));
        let mut f = File::create(&p).unwrap();
        writeln!(f, "barcode_4,barcode_3,barcode_2,barcode_1").unwrap();
        for r in rows {
            writeln!(f, "{},{},{},{}", r[0], r[1], r[2], r[3]).unwrap();
        }
        p
    }

    /// Segments in stored orientation: resolve expects revcomp of the observed bc1..bc4.
    fn stored(bc4: &str, bc3: &str, bc2: &str, bc1: &str) -> [Vec<u8>; 4] {
        [
            revcomp(bc1.as_bytes()),
            revcomp(bc2.as_bytes()),
            revcomp(bc3.as_bytes()),
            revcomp(bc4.as_bytes()),
        ]
    }

    fn refs(v: &[Vec<u8>; 4]) -> [&[u8]; 4] {
        [&v[0], &v[1], &v[2], &v[3]]
    }

    #[test]
    fn exact_edit1_and_unknown() {
        let p = write_wl(
            "exact",
            &[
                ["AAAAAAAA", "AAAAAA", "CCCCCC", "CCCCCCCC"], // row 0
                ["GGGGGGGG", "GGGGGG", "TTTTTT", "TTTTTTTT"], // row 1
            ],
        );
        let wl = Whitelist::from_csv(&p).unwrap();

        let r0 = stored("AAAAAAAA", "AAAAAA", "CCCCCC", "CCCCCCCC");
        assert_eq!(wl.resolve(refs(&r0)), [0, 0, 0, 0]);
        // exact rows pack to the all-A barcode (row 0 in every segment).
        assert_eq!(CellId::from_rows(wl.resolve(refs(&r0))).unwrap().render(), "A".repeat(16));

        let r1 = stored("GGGGGGGG", "GGGGGG", "TTTTTT", "TTTTTTTT");
        assert_eq!(wl.resolve(refs(&r1)), [1, 1, 1, 1]);

        // one substitution in bc1 still corrects to row 0.
        let mut e = stored("AAAAAAAA", "AAAAAA", "CCCCCC", "CCCCCCCC");
        e[0][7] = b'C'; // revcomp(CCCCCCCC)=GGGGGGGG -> GGGGGGGC, edit-1 from row 0
        assert_eq!(wl.resolve(refs(&e))[0], 0);

        // a segment absent and not edit-1 -> NOT_FOUND.
        let mut junk = stored("AAAAAAAA", "AAAAAA", "CCCCCC", "CCCCCCCC");
        junk[0] = b"ACGTACGT".to_vec();
        assert_eq!(wl.resolve(refs(&junk))[0], NOT_FOUND);

        let _ = std::fs::remove_file(&p);
    }

    #[test]
    fn edit1_ambiguity_is_rejected() {
        // row 0 bc1 revcomp = GGGGGGGG, row 1 bc1 revcomp = TGGGGGGG; "AGGGGGGG" is edit-1 from both
        // and exact-matches neither -> ambiguous.
        let p = write_wl(
            "ambig",
            &[
                ["AAAAAAAA", "AAAAAA", "CCCCCC", "CCCCCCCC"], // row 0
                ["GGGGGGGG", "GGGGGG", "TTTTTT", "CCCCCCCA"], // row 1
            ],
        );
        let wl = Whitelist::from_csv(&p).unwrap();

        let mut q = stored("AAAAAAAA", "AAAAAA", "CCCCCC", "CCCCCCCC");
        q[0] = b"AGGGGGGG".to_vec();
        let got = wl.resolve(refs(&q));
        assert_eq!(got[0], AMBIGUOUS);
        assert_eq!(&got[1..], &[0, 0, 0]);

        let _ = std::fs::remove_file(&p);
    }
}
