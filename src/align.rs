//! Internal alignment: map barcoded reads to the transcriptome with an embedded minimap2 and build
//! the eqclass directly, replacing the external `minimap2 -ax map-ont --for-only -N 200 -p 0.9` plus
//! BAM step. The pinned recipe is baked in so the wrong-flags mistake is impossible: `with_cigar()`
//! matches the CLI's base-level alignment (without it `-p 0.9` filters on chaining scores and the
//! retained secondary set diverges), `best_n = 200` is `-N`, `pri_ratio = 0.9` is `-p`, and the
//! FOR_ONLY flag is `--for-only`.

use std::ffi::CStr;
use std::io;
use std::path::Path;

use minimap2::{Aligner, Built};

use crate::eqclass::{parse_key, EqClass, Molecule};
use crate::fastq;
use crate::parallel;

const MM_F_FOR_ONLY: i64 = 0x100000; // minimap2 --for-only

/// Prune a read's competing `(transcript-id, internal-gap)` hits. `internal-gap` is the
/// inserted+deleted bases in that hit's alignment CIGAR (structural disagreement, excluding
/// mismatches and terminal clips). Keeps transcripts within `margin` bases of the smallest gap and
/// drops the rest, but only when the best fit clears `floor` (else nothing fits cleanly and all hits
/// are kept). Output is sorted and deduped. Empty in, empty out.
fn prune(gaps: &[(u32, u32)], margin: u32, floor: u32) -> Vec<u32> {
    let Some(best) = gaps.iter().map(|&(_, g)| g).min() else {
        return Vec::new();
    };
    let mut kept: Vec<u32> = if best > floor {
        gaps.iter().map(|&(t, _)| t).collect()
    } else {
        gaps.iter()
            .filter(|&&(_, g)| g - best < margin)
            .map(|&(t, _)| t)
            .collect()
    };
    kept.sort_unstable();
    kept.dedup();
    kept
}

/// Structural-resolution parameters: `margin` and `floor` in internal-indel bases (see [`prune`]).
#[derive(Clone, Copy)]
pub struct Resolve {
    pub margin: u32,
    pub floor: u32,
}

/// Build the eqclass by mapping the barcoded reads in `reads` to `reference`. Transcripts are the
/// reference sequences in index order (the BAM `@SQ` order), so `count` sees identical matrix
/// dimensions and transcript ids either way. Each mapped read becomes one molecule carrying its
/// non-supplementary target ids, sorted (matching the BAM path's tid order). With `resolve` set,
/// each read's hits are pruned to the structurally-best-fitting transcripts (see [`prune`]).
pub fn align_to_eqclass<P: AsRef<Path>>(
    reads: P,
    reference: P,
    v5_binid: bool,
    resolve: Option<Resolve>,
    workers: usize,
) -> io::Result<EqClass> {
    let mut aligner = Aligner::builder()
        .map_ont()
        .with_cigar()
        .with_index(reference.as_ref(), None)
        .map_err(|e| io::Error::new(io::ErrorKind::Other, format!("minimap2 index: {e}")))?;
    aligner.mapopt.best_n = 200;
    aligner.mapopt.pri_ratio = 0.9;
    aligner.mapopt.flag |= MM_F_FOR_ONLY;

    let n_seq = aligner.n_seq() as usize;
    let mut transcripts: Vec<(String, u32)> = Vec::with_capacity(n_seq);
    for i in 0..n_seq {
        let s = aligner
            .get_seq(i)
            .ok_or_else(|| io::Error::new(io::ErrorKind::Other, format!("ref seq {i} missing")))?;
        let name = unsafe { CStr::from_ptr(s.name) }
            .to_str()
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?
            .to_string();
        transcripts.push((name, s.len));
    }

    let mut reader = fastq::open_gz(reads)?;
    let molecules = parallel::run(
        || {
            reader.next().map(|rec| {
                rec.map(|r| (fastq::read_name(r.id()).to_vec(), r.seq().to_vec()))
                    .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
            })
        },
        workers,
        || (),
        |_: &mut (), (name, seq): (Vec<u8>, Vec<u8>)| {
            map_one(&aligner, &name, &seq, v5_binid, resolve)
        },
        |rx| -> io::Result<Vec<Molecule>> { Ok(rx.into_iter().flatten().collect()) },
    )?;

    Ok(EqClass {
        transcripts,
        molecules,
    })
}

/// Internal-indel gap of a hit's alignment: total inserted + deleted bases in the CIGAR (op 1 = I,
/// op 2 = D). `None` if the hit carries no CIGAR. Terminal soft-clips and mismatches are excluded,
/// so the gap measures structural disagreement, not read length or base-call noise.
fn internal_gap(m: &minimap2::Mapping) -> Option<u32> {
    let cigar = m.alignment.as_ref()?.cigar.as_ref()?;
    Some(
        cigar
            .iter()
            .filter(|&&(_, op)| op == 1 || op == 2)
            .map(|&(len, _)| len)
            .sum(),
    )
}

/// Map one read: `None` if its key is malformed or it is unmapped (skipped exactly as the BAM path
/// skips unmapped reads). Target ids come back sorted, so the molecule matches the BAM-derived one.
/// With `resolve` set, the hits are pruned to the structurally-best-fitting transcripts by
/// [`prune`]; a hit with no CIGAR cannot be assessed, so such a read is left unpruned.
fn map_one(
    aligner: &Aligner<Built>,
    name: &[u8],
    seq: &[u8],
    v5_binid: bool,
    resolve: Option<Resolve>,
) -> Option<Molecule> {
    let (cell, umi) = parse_key(name, v5_binid)?;
    let hits = aligner
        .map(seq, false, false, None, None, Some(name))
        .ok()?;
    let mapped = hits
        .iter()
        .filter(|m| !m.is_supplementary)
        .filter(|m| m.target_id >= 0);
    let mut txps: Vec<u32> = match resolve {
        Some(r) => {
            // One pass: (tid, internal-gap) per hit. `assessable` goes false if any hit lacks a
            // CIGAR, in which case the read is left unpruned (never drop on absent evidence).
            let mut gaps: Vec<(u32, u32)> = Vec::new();
            let mut assessable = true;
            for m in mapped {
                let tid = m.target_id as u32;
                match internal_gap(m) {
                    Some(g) => gaps.push((tid, g)),
                    None => {
                        assessable = false;
                        gaps.push((tid, 0)); // gap unused; the read is left unpruned below
                    }
                }
            }
            if assessable {
                prune(&gaps, r.margin, r.floor)
            } else {
                let mut t: Vec<u32> = gaps.into_iter().map(|(t, _)| t).collect();
                t.sort_unstable();
                t.dedup();
                t
            }
        }
        None => mapped.map(|m| m.target_id as u32).collect(),
    };
    if txps.is_empty() {
        return None;
    }
    txps.sort_unstable();
    Some(Molecule { cell, umi, txps })
}

#[cfg(test)]
mod tests {
    use super::*;
    use flate2::write::GzEncoder;
    use flate2::Compression;
    use std::io::Write;

    /// Deterministic high-complexity ACGT so minimap2 seeds cleanly (no external randomness).
    fn synth(seed: u64, n: usize) -> String {
        const B: [u8; 4] = [b'A', b'C', b'G', b'T'];
        let mut x = seed;
        (0..n)
            .map(|_| {
                x = x
                    .wrapping_mul(6364136223846793005)
                    .wrapping_add(1442695040888963407);
                B[((x >> 33) & 3) as usize] as char
            })
            .collect()
    }

    #[test]
    fn maps_read_to_its_transcript() {
        let dir = std::env::temp_dir().join(format!("bp_align_{}", std::process::id()));
        let _ = std::fs::remove_dir_all(&dir);
        std::fs::create_dir_all(&dir).unwrap();

        // Two unrelated transcripts; a barcoded read equal to TXP0 should map to tid 0 only.
        let (txp0, txp1) = (synth(1, 600), synth(2, 600));
        let refp = dir.join("ref.fa");
        std::fs::File::create(&refp)
            .unwrap()
            .write_all(format!(">TXP0\n{txp0}\n>TXP1\n{txp1}\n").as_bytes())
            .unwrap();

        let readsp = dir.join("reads.fa.gz");
        let mut gz = GzEncoder::new(
            std::fs::File::create(&readsp).unwrap(),
            Compression::default(),
        );
        gz.write_all(format!(">r1_AAACGTTGCAGAACAC_ACGTACGTACGT\n{txp0}\n").as_bytes())
            .unwrap();
        gz.finish().unwrap();

        let eq = align_to_eqclass(
            &readsp,
            &refp,
            false,
            None,
            crate::parallel::default_workers(),
        )
        .unwrap();
        assert_eq!(eq.transcripts.len(), 2);
        assert_eq!(eq.transcripts[0].0, "TXP0");
        assert_eq!(eq.molecules.len(), 1, "one mapped molecule");
        let m = &eq.molecules[0];
        assert_eq!(m.cell.render(), "AAACGTTGCAGAACAC");
        assert_eq!(m.umi.render(), "ACGTACGTACGT");
        assert_eq!(m.txps, vec![0], "maps to TXP0 (tid 0) only");

        let _ = std::fs::remove_dir_all(&dir);
    }

    #[test]
    fn prune_clean_cassette_drops_the_gapped_isoform() {
        assert_eq!(prune(&[(0, 0), (1, 150)], 20, 20), vec![0]);
    }

    #[test]
    fn prune_shared_region_keeps_both() {
        assert_eq!(prune(&[(0, 0), (1, 0)], 20, 20), vec![0, 1]);
    }

    #[test]
    fn prune_sub_margin_noise_keeps_both() {
        assert_eq!(prune(&[(0, 0), (1, 5)], 20, 20), vec![0, 1]);
    }

    #[test]
    fn prune_partial_drops_only_the_far_one() {
        assert_eq!(prune(&[(0, 0), (1, 2), (2, 150)], 20, 20), vec![0, 1]);
    }

    #[test]
    fn prune_margin_boundary_is_exclusive() {
        // gap - best == margin is dropped; strictly within margin is kept.
        assert_eq!(prune(&[(0, 0), (1, 20)], 20, 20), vec![0]);
    }

    #[test]
    fn prune_no_clean_fit_keeps_all() {
        // best gap 100 exceeds the floor, so nothing is trusted and all hits stay.
        assert_eq!(prune(&[(0, 100), (1, 300)], 20, 20), vec![0, 1]);
    }

    #[test]
    fn prune_dedups_repeated_transcript() {
        // A transcript can appear twice (two alignments); the output is a set.
        assert_eq!(prune(&[(0, 0), (0, 5), (1, 150)], 20, 20), vec![0]);
    }

    #[test]
    fn prune_single_hit_passes_through() {
        assert_eq!(prune(&[(0, 0)], 20, 20), vec![0]);
    }

    #[test]
    fn prune_empty_is_empty() {
        assert!(prune(&[], 20, 20).is_empty());
    }

    #[test]
    fn structural_resolve_drops_the_spanned_cassette_only() {
        // TXP0 carries a 120 bp cassette between long shared flanks; TXP1 skips it. A read that
        // spans the cassette aligns to TXP1 with a 120 bp insertion, so resolution drops TXP1. A
        // read living in the shared flank fits both cleanly and stays ambiguous.
        let dir = std::env::temp_dir().join(format!("bp_struct_{}", std::process::id()));
        let _ = std::fs::remove_dir_all(&dir);
        std::fs::create_dir_all(&dir).unwrap();

        let (fl, cassette, fr) = (synth(20, 2000), synth(21, 120), synth(22, 2000));
        let ra = format!("{fl}{cassette}{fr}");
        let ro = format!("{fl}{fr}");
        let refp = dir.join("ref.fa");
        std::fs::write(&refp, format!(">TXP0\n{ra}\n>TXP1\n{ro}\n")).unwrap();

        let txps = |seq: &str, resolve: Option<Resolve>| -> Vec<u32> {
            let readsp = dir.join("r.fa.gz");
            let mut gz = GzEncoder::new(
                std::fs::File::create(&readsp).unwrap(),
                Compression::default(),
            );
            gz.write_all(format!(">r_AAACGTTGCAGAACAC_ACGTACGTACGT\n{seq}\n").as_bytes())
                .unwrap();
            gz.finish().unwrap();
            let w = crate::parallel::default_workers();
            align_to_eqclass(&readsp, &refp, false, resolve, w)
                .unwrap()
                .molecules[0]
                .txps
                .clone()
        };

        let r = Some(Resolve {
            margin: 20,
            floor: 20,
        });
        assert_eq!(
            txps(&ra, r),
            vec![0],
            "cassette read resolves to its isoform"
        );
        assert_eq!(txps(&ra, None), vec![0, 1], "default keeps both");
        assert_eq!(
            txps(&fr, r),
            vec![0, 1],
            "shared-region read stays ambiguous"
        );

        let _ = std::fs::remove_dir_all(&dir);
    }
}
