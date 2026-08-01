//! Differential test: the ported seal aligner must reproduce rust-bio's `Aligner::custom` extraction
//! on a battery of synthetic reads (clean, substituted, indel'd). bio is a dev-dependency only, so
//! the shipped tool carries no aligner dependency. Any divergence fails the test with the read.

use bio::alignment::pairwise::{Aligner, Scoring, MIN_SCORE};
use bio::alignment::AlignmentOperation;

/// Reference seal, a verbatim transcription of the pre-rewrite `algo::seal_extraction` on bio.
fn reference_seal(seq1: &str) -> Option<(Vec<String>, String, String)> {
    let x = seq1.as_bytes();
    let y = b"CTCGAVVVVVVCTCVVVVVVCAT";
    let scoring = Scoring::from_scores(-1, -2, 2, -1).xclip(0).yclip(MIN_SCORE);
    let mut aligner = Aligner::with_scoring(scoring);
    let alignment = aligner.custom(x, y);

    let xstart = alignment.xstart;
    if xstart < 40 || (alignment.xend + 8) > seq1.len() {
        return None;
    }
    let txp_seq = seq1[0..(xstart - 20)].to_string();
    let umi = seq1[(xstart - 20)..(xstart - 8)].to_string();
    let mut barcodes = vec![seq1[(xstart - 8)..xstart].to_string()];

    let mut step = alignment.operations.iter();
    let mut offset = match step.next().unwrap() {
        AlignmentOperation::Xclip(ct) => *ct,
        _ => unreachable!("{:?}", alignment.operations),
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
    update(&mut offset, 5);
    let s = offset;
    update(&mut offset, 6);
    barcodes.push(seq1[s..offset].to_string());
    update(&mut offset, 3);
    let s = offset;
    update(&mut offset, 6);
    barcodes.push(seq1[s..offset].to_string());
    update(&mut offset, 3);
    barcodes.push(seq1[offset..offset + 8].to_string());
    barcodes.reverse();
    Some((barcodes, umi, txp_seq))
}

fn mine(seq: &str) -> Option<(Vec<String>, String, String)> {
    let (bcs, umi, txp) = bagpiper::seal::seal_extraction(seq.as_bytes())?;
    Some((
        bcs.into_iter().map(|b| String::from_utf8(b).unwrap()).collect(),
        String::from_utf8(umi).unwrap(),
        String::from_utf8(txp).unwrap(),
    ))
}

struct Rng(u64);
impl Rng {
    fn next(&mut self) -> u64 {
        self.0 = self.0.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
        self.0 >> 33
    }
    fn under(&mut self, n: usize) -> usize {
        (self.next() as usize) % n
    }
    fn dna(&mut self, n: usize) -> String {
        (0..n).map(|_| b"ACGT"[self.under(4)] as char).collect()
    }
}

fn assemble(txp: &str, umi: &str, a: &str, b: &str, c: &str, d: &str, tail: &str) -> String {
    format!("{txp}{umi}{a}CTCGA{b}CTC{c}CAT{d}{tail}")
}

#[test]
fn seal_matches_bio_on_synthetic_reads() {
    let mut rng = Rng(0x9E3779B97F4A7C15);
    let mut checked = 0;
    let mut disagree: Vec<String> = Vec::new();

    for case in 0..1800 {
        // txp >= 22 so xstart = txp+12+8 >= 42 (clears the 40 nt floor).
        let txp_len = 22 + rng.under(40);
        let txp = rng.dna(txp_len);
        let umi = rng.dna(12);
        let (mut a, mut b, mut c, mut d) =
            (rng.dna(8), rng.dna(6), rng.dna(6), rng.dna(8));
        let tail_len = 8 + rng.under(12);
        let tail = rng.dna(tail_len);

        // A third clean, a third substituted, a third with a single indel in a barcode segment.
        let kind = case % 3;
        if kind == 1 {
            let n_sub = 1 + rng.under(3);
            for _ in 0..n_sub {
                let which = rng.under(4);
                let seg = match which {
                    0 => &mut a,
                    1 => &mut b,
                    2 => &mut c,
                    _ => &mut d,
                };
                let pos = rng.under(seg.len());
                let mut bytes = seg.as_bytes().to_vec();
                bytes[pos] = b"ACGT"[rng.under(4)];
                *seg = String::from_utf8(bytes).unwrap();
            }
        } else if kind == 2 {
            // insert or delete one base in b or c (leave anchors intact)
            let which = if rng.under(2) == 0 { &mut b } else { &mut c };
            if rng.under(2) == 0 {
                let pos = rng.under(which.len());
                let mut bytes = which.clone().into_bytes();
                bytes.insert(pos, b"ACGT"[rng.under(4)]);
                *which = String::from_utf8(bytes).unwrap();
            } else {
                let pos = rng.under(which.len());
                let mut bytes = which.clone().into_bytes();
                bytes.remove(pos);
                *which = String::from_utf8(bytes).unwrap();
            }
        }

        let read = assemble(&txp, &umi, &a, &b, &c, &d, &tail);
        checked += 1;
        if mine(&read) != reference_seal(&read) {
            if disagree.len() < 5 {
                disagree.push(format!(
                    "read={read}\n  mine={:?}\n  ref ={:?}",
                    mine(&read),
                    reference_seal(&read)
                ));
            }
        }
    }

    eprintln!("seal differential: checked {checked}, disagreements {}", disagree.len());
    assert!(
        disagree.is_empty(),
        "seal diverged from bio on {} reads:\n{}",
        disagree.len(),
        disagree.join("\n")
    );
}
