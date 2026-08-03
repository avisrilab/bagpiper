//! Golden test for the barcode drivers: run them on the synthetic fixtures and compare the
//! decompressed output to the committed golden record sets (sorted, since order is not contractual).
//! The fixtures are all regex-path (no seal); the seal is covered by `seal_diff.rs`.

use std::io::Read;
use std::path::{Path, PathBuf};

fn md() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
}

/// Decompress a `.fa.gz` and fold each 2-line FASTA record into one sorted `>id\tseq` line.
fn records_sorted(fa_gz: &Path) -> Vec<String> {
    let mut s = String::new();
    flate2::read::GzDecoder::new(std::fs::File::open(fa_gz).unwrap())
        .read_to_string(&mut s)
        .unwrap();
    let lines: Vec<&str> = s.lines().collect();
    let mut recs: Vec<String> = lines
        .chunks(2)
        .map(|c| format!("{}\t{}", c[0], c.get(1).copied().unwrap_or("")))
        .collect();
    recs.sort();
    recs
}

fn golden(name: &str) -> Vec<String> {
    let p = md().join("tests/fixtures/golden").join(name);
    let mut v: Vec<String> = std::fs::read_to_string(&p)
        .unwrap_or_else(|e| panic!("read {:?}: {}", p, e))
        .lines()
        .map(str::to_string)
        .collect();
    v.sort();
    v
}

fn assert_records_eq(got: &[String], want: &[String], label: &str) {
    if got == want {
        return;
    }
    let only_got: Vec<&String> = got.iter().filter(|r| !want.contains(r)).collect();
    let only_want: Vec<&String> = want.iter().filter(|r| !got.contains(r)).collect();
    panic!(
        "{label} drifted from golden\n  in output not golden ({}): {:#?}\n  in golden not output ({}): {:#?}",
        only_got.len(),
        only_got,
        only_want.len(),
        only_want
    );
}

#[test]
fn nanopore_matches_golden() {
    let out = std::env::temp_dir().join(format!("bp_gold_np_{}", std::process::id()));
    let _ = std::fs::remove_dir_all(&out);
    std::fs::create_dir_all(&out).unwrap();
    bagpiper::barcode::run_nanopore(
        &[md().join("tests/fixtures/barcode_nanopore.fq.gz")],
        &md().join("tests/whitelist/synthetic_barcodes.csv"),
        &out,
    )
    .unwrap();
    assert_records_eq(
        &records_sorted(&out.join("passed.bcd.nanopore.fa.gz")),
        &golden("passed.tsv"),
        "passed.bcd.nanopore",
    );
    assert_records_eq(
        &records_sorted(&out.join("failed.bcd.nanopore.fa.gz")),
        &golden("failed.tsv"),
        "failed.bcd.nanopore",
    );
    let _ = std::fs::remove_dir_all(&out);
}

#[test]
fn illumina_matches_golden() {
    let out = std::env::temp_dir().join(format!("bp_gold_ill_{}", std::process::id()));
    let _ = std::fs::remove_dir_all(&out);
    std::fs::create_dir_all(&out).unwrap();
    bagpiper::barcode::run_illumina(
        &md().join("tests/fixtures/barcode_illumina_r1.fq.gz"),
        &md().join("tests/fixtures/barcode_illumina_r2.fq.gz"),
        &md().join("tests/whitelist/synthetic_barcodes.csv"),
        &out,
    )
    .unwrap();
    assert_records_eq(
        &records_sorted(&out.join("read.bcd.illumina.fa.gz")),
        &golden("illumina.tsv"),
        "read.bcd.illumina",
    );
    let _ = std::fs::remove_dir_all(&out);
}
