//! Differential test: the V5 TSO seal must reproduce a rust-bio `Aligner::custom` extraction on a
//! battery of synthetic reads (clean, substituted, indel'd anchors). bio is a dev-dependency only.
//! The reference is a transcription of the pre-rewrite tso extraction on bio; any divergence fails.

use bio::alignment::pairwise::{Aligner, Scoring, MIN_SCORE};
use bio::alignment::AlignmentOperation;

/// Reference TSO seal on bio: score-gate at 16, extract the two 6 nt UMIs off the anchor, trim the
/// cDNA to after the seal, reject a wrong-length UMI. Returns (umi1, umi2, txp, seal_start).
fn reference_tso(seq1: &str) -> Option<(String, String, String, usize)> {
    if seq1.len() < 40 {
        return None;
    }
    let x = seq1.as_bytes();
    let y = b"CGCAGAGTVVVVVVTACTGVVVVVVGAAT";
    let scoring = Scoring::from_scores(-1, -2, 2, -1)
        .xclip(0)
        .yclip(MIN_SCORE);
    let mut aligner = Aligner::with_scoring(scoring);
    let alignment = aligner.custom(x, y);
    if alignment.score < 16 {
        return None;
    }
    let seal_start = alignment.xstart;
    let txp_seq = seq1[alignment.xend..seq1.len()].to_string();

    let mut step = alignment.operations.iter();
    let mut offset = match step.next().unwrap() {
        AlignmentOperation::Xclip(ct) => *ct,
        _ => 0,
    };
    let mut update = |offset: &mut usize, xlim: usize| {
        let mut check = 0;
        loop {
            match step.next().unwrap() {
                AlignmentOperation::Match | AlignmentOperation::Subst => check += 1,
                AlignmentOperation::Del => {
                    *offset -= 1;
                    check += 1;
                }
                AlignmentOperation::Ins => *offset += 1,
                _ => unreachable!("{:?}", alignment.operations),
            }
            if check == xlim {
                *offset += xlim;
                break;
            }
        }
    };
    update(&mut offset, 8); // CGCAGAGT
    let s = offset;
    update(&mut offset, 6); // UMI1
    let umi1 = seq1[s..offset].to_string();
    update(&mut offset, 5); // TACTG
    let s = offset;
    update(&mut offset, 6); // UMI2
    let umi2 = seq1[s..offset].to_string();

    if umi1.len() != 6 || umi2.len() != 6 {
        return None;
    }
    Some((umi1, umi2, txp_seq, seal_start))
}

fn mine(seq: &str) -> Option<(String, String, String, usize)> {
    let (u1, u2, txp, start) = bagpiper::tso::seal_extraction(seq.as_bytes())?;
    Some((
        String::from_utf8(u1).unwrap(),
        String::from_utf8(u2).unwrap(),
        String::from_utf8(txp).unwrap(),
        start,
    ))
}

struct Rng(u64);
impl Rng {
    fn next(&mut self) -> u64 {
        self.0 = self
            .0
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        self.0 >> 33
    }
    fn under(&mut self, n: usize) -> usize {
        (self.next() as usize) % n
    }
    fn dna(&mut self, n: usize) -> String {
        (0..n).map(|_| b"ACGT"[self.under(4)] as char).collect()
    }
}

#[test]
fn tso_matches_bio_on_synthetic_reads() {
    let mut rng = Rng(0x243F6A8885A308D3);
    let mut checked = 0;
    let mut disagree: Vec<String> = Vec::new();

    for case in 0..1800 {
        let leftover_len = 5 + rng.under(30); // >=5 so bio always emits a leading Xclip
        let leftover = rng.dna(leftover_len);
        let (umi1, umi2) = (rng.dna(6), rng.dna(6));
        let cdna_len = 15 + rng.under(40);
        let cdna = rng.dna(cdna_len);
        // Anchors, perturbed a third of the time each way; the aligner is what could diverge.
        let (mut a, mut b, mut c) = (
            "CGCAGAGT".to_string(),
            "TACTG".to_string(),
            "GAAT".to_string(),
        );
        match case % 3 {
            1 => {
                // 1-2 substitutions across the anchors (the aligner is what could diverge)
                for _ in 0..1 + rng.under(2) {
                    let seg = match rng.under(3) {
                        0 => &mut a,
                        1 => &mut b,
                        _ => &mut c,
                    };
                    let pos = rng.under(seg.len());
                    let mut by = seg.as_bytes().to_vec();
                    by[pos] = b"ACGT"[rng.under(4)];
                    *seg = String::from_utf8(by).unwrap();
                }
            }
            2 => {
                // one indel in an anchor
                let seg = match rng.under(3) {
                    0 => &mut a,
                    1 => &mut b,
                    _ => &mut c,
                };
                let pos = rng.under(seg.len());
                let mut by = seg.clone().into_bytes();
                if rng.under(2) == 0 {
                    by.insert(pos, b"ACGT"[rng.under(4)]);
                } else {
                    by.remove(pos);
                }
                *seg = String::from_utf8(by).unwrap();
            }
            _ => {}
        }

        let read = format!("{leftover}{a}{umi1}{b}{umi2}{c}{cdna}");
        checked += 1;
        if mine(&read) != reference_tso(&read) && disagree.len() < 5 {
            disagree.push(format!(
                "read={read}\n  mine={:?}\n  ref ={:?}",
                mine(&read),
                reference_tso(&read)
            ));
        }
    }

    eprintln!(
        "tso differential: checked {checked}, disagreements {}",
        disagree.len()
    );
    assert!(
        disagree.is_empty(),
        "tso diverged from bio on {} reads:\n{}",
        disagree.len(),
        disagree.join("\n")
    );
}
