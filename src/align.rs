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

/// Build the eqclass by mapping the barcoded reads in `reads` to `reference`. Transcripts are the
/// reference sequences in index order (the BAM `@SQ` order), so `count` sees identical matrix
/// dimensions and transcript ids either way. Each mapped read becomes one molecule carrying its
/// non-supplementary target ids, sorted (matching the BAM path's tid order).
pub fn align_to_eqclass<P: AsRef<Path>>(
    reads: P,
    reference: P,
    v5_binid: bool,
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
        |_: &mut (), (name, seq): (Vec<u8>, Vec<u8>)| map_one(&aligner, &name, &seq, v5_binid),
        |rx| -> io::Result<Vec<Molecule>> { Ok(rx.into_iter().flatten().collect()) },
    )?;

    Ok(EqClass {
        transcripts,
        molecules,
    })
}

/// Map one read: `None` if its key is malformed or it is unmapped (skipped exactly as the BAM path
/// skips unmapped reads). Target ids come back sorted, so the molecule matches the BAM-derived one.
fn map_one(aligner: &Aligner<Built>, name: &[u8], seq: &[u8], v5_binid: bool) -> Option<Molecule> {
    let (cell, umi) = parse_key(name, v5_binid)?;
    let hits = aligner
        .map(seq, false, false, None, None, Some(name))
        .ok()?;
    let mut txps: Vec<u32> = hits
        .iter()
        .filter(|m| !m.is_supplementary)
        .filter_map(|m| (m.target_id >= 0).then_some(m.target_id as u32))
        .collect();
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

        let eq =
            align_to_eqclass(&readsp, &refp, false, crate::parallel::default_workers()).unwrap();
        assert_eq!(eq.transcripts.len(), 2);
        assert_eq!(eq.transcripts[0].0, "TXP0");
        assert_eq!(eq.molecules.len(), 1, "one mapped molecule");
        let m = &eq.molecules[0];
        assert_eq!(m.cell.render(), "AAACGTTGCAGAACAC");
        assert_eq!(m.umi.render(), "ACGTACGTACGT");
        assert_eq!(m.txps, vec![0], "maps to TXP0 (tid 0) only");

        let _ = std::fs::remove_dir_all(&dir);
    }
}
