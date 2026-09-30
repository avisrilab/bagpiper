//! End-to-end test: run the built binary on the simulated fixture in `tests/e2e/fixture` and compare
//! against the recorded output of the same release in `tests/e2e/expected/<version>/`.
//!
//! Legs: (A) `barcode` on the fixture reads; (B) `count --b1` on the frozen BAM; (C) `count --r1
//! --reference` on leg A's output, only when the binary has the internal aligner. Legs B and C also
//! check the matrix against the simulated truth, so a recorded output cannot be wrong unnoticed.
//!
//! BAGPIPER_BIN=<path> tests another build (CI uses it to test older releases with these tests);
//! the release is read from `--version`. BAGPIPER_E2E_BLESS=1 writes the expected files instead of
//! comparing; run it only for a new release, from that release's own build.

use std::collections::{BTreeMap, HashMap};
use std::fs;
use std::io::{BufRead, BufReader, Read};
use std::path::{Path, PathBuf};
use std::process::Command;

/// Relative tolerance for matrix values: the EM runs in f32.
const REL_TOL: f64 = 1e-6;
/// Dominant-isoform accuracy required on genes with >= 10% unique sequence.
const MIN_IDENTIFIABLE_ACC: f64 = 0.95;

fn md() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
}

fn fixture(name: &str) -> PathBuf {
    md().join("tests/e2e/fixture").join(name)
}

fn bin() -> PathBuf {
    std::env::var("BAGPIPER_BIN")
        .map(PathBuf::from)
        .unwrap_or_else(|_| PathBuf::from(env!("CARGO_BIN_EXE_bagpiper")))
}

fn bless() -> bool {
    std::env::var("BAGPIPER_E2E_BLESS").is_ok_and(|v| v == "1")
}

/// Run the binary; panic on failure; return stdout + stderr.
fn run(args: &[&str]) -> String {
    let out = Command::new(bin())
        .args(args)
        .output()
        .expect("run bagpiper");
    let text = format!(
        "{}{}",
        String::from_utf8_lossy(&out.stdout),
        String::from_utf8_lossy(&out.stderr)
    );
    assert!(out.status.success(), "bagpiper {args:?} failed:\n{text}");
    text
}

fn version() -> String {
    run(&["--version"])
        .split_whitespace()
        .nth(1)
        .expect("version string")
        .to_string()
}

fn expected_dir() -> PathBuf {
    md().join("tests/e2e/expected").join(version())
}

fn tmp(leg: &str) -> PathBuf {
    let d = std::env::temp_dir().join(format!("bp_e2e_{leg}_{}", std::process::id()));
    let _ = fs::remove_dir_all(&d);
    fs::create_dir_all(&d).unwrap();
    d
}

fn read_gz(p: &Path) -> String {
    let mut s = String::new();
    flate2::read::MultiGzDecoder::new(fs::File::open(p).unwrap())
        .read_to_string(&mut s)
        .unwrap();
    s
}

/// FNV-1a 64: a stable digest for recorded sequences (std's hasher is not stable across releases).
fn fnv(s: &str) -> u64 {
    s.bytes().fold(0xcbf29ce484222325, |h, b| {
        (h ^ b as u64).wrapping_mul(0x100000001b3)
    })
}

/// Compare `got` lines against the recorded file, or record them in bless mode.
fn check_lines(name: &str, got: &[String]) {
    let p = expected_dir().join(name);
    if bless() {
        fs::create_dir_all(expected_dir()).unwrap();
        fs::write(&p, got.join("\n") + "\n").unwrap();
        return;
    }
    let want: Vec<String> = fs::read_to_string(&p)
        .unwrap_or_else(|_| panic!("no recorded output {} for this release", p.display()))
        .lines()
        .map(String::from)
        .collect();
    let diff: Vec<_> = got
        .iter()
        .zip(&want)
        .filter(|(g, w)| g != w)
        .take(5)
        .collect();
    assert!(
        got.len() == want.len() && diff.is_empty(),
        "{name}: {} lines vs {} recorded; first differences {diff:?}",
        got.len(),
        want.len()
    );
}

/// Passed records as `read  CB  UMI  cDNA-length  cDNA-digest`, failed records as the read name.
fn barcode_records(out: &Path) -> (Vec<String>, Vec<String>) {
    let fa = |name: &str| -> Vec<(String, String)> {
        let text = read_gz(&out.join(name));
        let mut it = text.lines();
        let mut v = Vec::new();
        while let (Some(h), Some(s)) = (it.next(), it.next()) {
            v.push((
                h[1..].split_whitespace().next().unwrap().to_string(),
                s.to_string(),
            ));
        }
        v
    };
    let mut passed: Vec<String> = fa("passed.bcd.nanopore.fa.gz")
        .into_iter()
        .map(|(h, s)| {
            let mut p = h.rsplitn(3, '_');
            let (umi, cb, read) = (p.next().unwrap(), p.next().unwrap(), p.next().unwrap());
            format!("{read}\t{cb}\t{umi}\t{}\t{:016x}", s.len(), fnv(&s))
        })
        .collect();
    let mut failed: Vec<String> = fa("failed.bcd.nanopore.fa.gz")
        .into_iter()
        .map(|(h, _)| h)
        .collect();
    passed.sort();
    failed.sort();
    (passed, failed)
}

fn run_barcode(out: &Path) -> String {
    run(&[
        "barcode",
        "--r1",
        fixture("reads.fq.gz").to_str().unwrap(),
        "--whitelist",
        md().join("tests/whitelist/synthetic_barcodes.csv")
            .to_str()
            .unwrap(),
        "--nanopore",
        "-o",
        out.to_str().unwrap(),
    ])
}

/// `(barcode, transcript) -> value` from a count output directory.
fn matrix(dir: &Path) -> BTreeMap<(String, String), f64> {
    let bcs: Vec<String> = read_gz(&dir.join("barcodes.tsv.gz"))
        .lines()
        .map(String::from)
        .collect();
    let fts: Vec<String> = read_gz(&dir.join("features.tsv.gz"))
        .lines()
        .map(|l| l.split('\t').next().unwrap().to_string())
        .collect();
    let mtx = read_gz(&dir.join("matrix.mtx.gz"));
    let mut lines = mtx.lines().skip_while(|l| l.starts_with('%'));
    lines.next(); // dimensions
    lines
        .map(|l| {
            let f: Vec<&str> = l.split_whitespace().collect();
            let (r, c): (usize, usize) = (f[0].parse().unwrap(), f[1].parse().unwrap());
            (
                (bcs[r - 1].clone(), fts[c - 1].clone()),
                f[2].parse().unwrap(),
            )
        })
        .collect()
}

fn check_matrix(name: &str, got: &BTreeMap<(String, String), f64>) {
    let lines: Vec<String> = got
        .iter()
        .map(|((b, t), v)| format!("{b}\t{t}\t{v}"))
        .collect();
    let p = expected_dir().join(name);
    if bless() {
        return check_lines(name, &lines);
    }
    let want: BTreeMap<(String, String), f64> = fs::read_to_string(&p)
        .unwrap_or_else(|_| panic!("no recorded output {} for this release", p.display()))
        .lines()
        .map(|l| {
            let f: Vec<&str> = l.split('\t').collect();
            ((f[0].to_string(), f[1].to_string()), f[2].parse().unwrap())
        })
        .collect();
    assert_eq!(
        got.keys().collect::<Vec<_>>(),
        want.keys().collect::<Vec<_>>(),
        "{name}: nonzero entries differ"
    );
    for (k, w) in &want {
        let g = got[k];
        assert!(
            (g - w).abs() <= REL_TOL * g.abs().max(w.abs()),
            "{name}: {k:?} = {g}, recorded {w}"
        );
    }
}

fn tsv(name: &str) -> Vec<Vec<String>> {
    BufReader::new(fs::File::open(fixture(name)).unwrap())
        .lines()
        .skip(1)
        .map(|l| l.unwrap().split('\t').map(String::from).collect())
        .collect()
}

/// Truth check: every simulated cell is called (a barcode with >= 10% of the largest barcode's
/// counts) and mapped to it, and on genes with >= 10% unique sequence the dominant isoform is right
/// in >= MIN_IDENTIFIABLE_ACC of cell-genes (an empty call counts as wrong).
fn check_truth(m: &BTreeMap<(String, String), f64>, cb_cell: &HashMap<String, String>) {
    let mut per_bc: HashMap<&str, f64> = HashMap::new();
    for ((b, _), v) in m {
        *per_bc.entry(b).or_default() += v;
    }
    let max = per_bc.values().cloned().fold(0.0, f64::max);
    let mut est: HashMap<(String, String), f64> = HashMap::new();
    let mut called = std::collections::BTreeSet::new();
    for ((b, t), v) in m {
        if per_bc[b.as_str()] >= 0.1 * max {
            let cell = cb_cell
                .get(b)
                .unwrap_or_else(|| panic!("called barcode {b} maps to no cell"));
            called.insert(cell.clone());
            *est.entry((cell.clone(), t.clone())).or_default() += v;
        }
    }
    let cells: Vec<String> = tsv("cells.tsv").into_iter().map(|r| r[0].clone()).collect();
    assert_eq!(
        called.into_iter().collect::<Vec<_>>(),
        cells,
        "called cells"
    );

    let mut truth: HashMap<(String, String), f64> = HashMap::new();
    for r in tsv("truth.tsv") {
        truth.insert((r[0].clone(), r[1].clone()), r[2].parse().unwrap());
    }
    let mut genes: BTreeMap<String, (f64, Vec<String>)> = BTreeMap::new();
    for r in tsv("genes.tsv") {
        let e = genes
            .entry(r[1].clone())
            .or_insert((r[4].parse().unwrap(), Vec::new()));
        e.1.push(r[0].clone());
    }
    let (mut ok, mut n, mut ok_all, mut n_all) = (0, 0, 0, 0);
    for cell in &cells {
        for (uniq, txs) in genes.values() {
            let get = |h: &HashMap<(String, String), f64>, t: &String| {
                h.get(&(cell.clone(), t.clone())).cloned().unwrap_or(0.0)
            };
            let argmax = |h: &HashMap<(String, String), f64>| {
                txs.iter()
                    .filter(|t| get(h, t) > 0.0)
                    .max_by(|a, b| get(h, a).partial_cmp(&get(h, b)).unwrap())
            };
            let right = argmax(&est).is_some() && argmax(&est) == argmax(&truth);
            n_all += 1;
            ok_all += right as usize;
            if *uniq >= 0.10 {
                n += 1;
                ok += right as usize;
            }
        }
    }
    eprintln!("dominant isoform: identifiable {ok}/{n}, all genes {ok_all}/{n_all}");
    assert!(
        ok as f64 >= MIN_IDENTIFIABLE_ACC * n as f64,
        "identifiable-gene accuracy {ok}/{n} below {MIN_IDENTIFIABLE_ACC}"
    );
}

/// The per-run assignment-path counters, when the release prints them, must sum to `matched`.
fn check_counters(log: &str) {
    let Some(line) = log.lines().find(|l| l.starts_with("Assignment path")) else {
        return;
    };
    let nums: Vec<u64> = line
        .split(|c: char| !c.is_ascii_digit() && c != '.')
        .filter(|s| !s.is_empty() && !s.contains('.'))
        .map(|s| s.parse().unwrap())
        .collect();
    // "(of N matched): ... : a (...) ... : b (...) ... : c (...)"
    assert_eq!(
        nums[1] + nums[2] + nums[3],
        nums[0],
        "counters do not sum: {line}"
    );
}

#[test]
fn leg_a_barcode() {
    let out = tmp("a");
    let log = run_barcode(&out);
    check_counters(&log);
    // The run totals, normalized across releases: 0.1.x prints "Total: N  Matched: N ...", 0.2.0
    // logs "[INFO] total N  matched N ...".
    let totals: Vec<String> = log
        .lines()
        .filter_map(|l| {
            let l = l.to_lowercase();
            let i = l.find("total:").or_else(|| l.find("] total "))?;
            let tail = l[i..].trim_start_matches("] ").replace(':', "");
            Some(tail.split_whitespace().collect::<Vec<_>>().join(" "))
        })
        .collect();
    check_lines("totals.txt", &totals);
    let (passed, failed) = barcode_records(&out);
    check_lines("passed.tsv", &passed);
    check_lines("failed.tsv", &failed);
    let _ = fs::remove_dir_all(&out);
}

#[test]
fn leg_b_count_frozen_bam() {
    let out = tmp("b");
    run(&[
        "count",
        "--b1",
        fixture("aligned.bam").to_str().unwrap(),
        "-o",
        out.to_str().unwrap(),
    ]);
    let m = matrix(&out);
    check_matrix("matrix.tsv", &m);
    let cb_cell: HashMap<String, String> = tsv("aligned.cb_cell.tsv")
        .into_iter()
        .map(|r| (r[0].clone(), r[1].clone()))
        .collect();
    check_truth(&m, &cb_cell);
    let _ = fs::remove_dir_all(&out);
}

#[test]
fn leg_c_count_internal_aligner() {
    if !run(&["count", "--help"]).contains("--r1") {
        eprintln!(
            "release {} has no internal aligner; leg C skipped",
            version()
        );
        return;
    }
    let bc = tmp("c_bc");
    run_barcode(&bc);
    let out = tmp("c");
    run(&[
        "count",
        "--r1",
        bc.join("passed.bcd.nanopore.fa.gz").to_str().unwrap(),
        "--reference",
        fixture("reference.fa").to_str().unwrap(),
        "-o",
        out.to_str().unwrap(),
    ]);
    let m = matrix(&out);
    check_matrix("matrix_internal.tsv", &m);
    // barcode -> cell by majority of the simulated read names behind it
    let mut votes: HashMap<String, HashMap<String, usize>> = HashMap::new();
    for rec in barcode_records(&bc).0 {
        let f: Vec<&str> = rec.split('\t').collect();
        let cell = f[0].split('.').next().unwrap().to_string();
        *votes
            .entry(f[1].to_string())
            .or_default()
            .entry(cell)
            .or_default() += 1;
    }
    let cb_cell = votes
        .into_iter()
        .map(|(b, v)| (b, v.into_iter().max_by_key(|(_, n)| *n).unwrap().0))
        .collect();
    check_truth(&m, &cb_cell);
    let _ = fs::remove_dir_all(&bc);
    let _ = fs::remove_dir_all(&out);
}
