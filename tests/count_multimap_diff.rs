//! Synthetic multi-mapping differential test for the count/EM path.
//!
//! k562_real.bam has no multi-transcript molecules, so the count oracle never exercises the
//! length-weighted EM split. This builds eqclasses with unique anchors on both transcripts of a
//! shared class, which forces a fractional split instead of a collapse to the shortest transcript.
//!
//! Always-on arm: the multi-mapping path must produce fractional abundances and conserve mass.
//! Gated arm (BP_TEST_REF_BIN = path to the reference binary; its runtime libs must be on the
//! inherited DYLD path): my matrix must match the reference `count` within EM tolerance.

use std::collections::HashMap;
use std::io::{Read, Write};
use std::path::Path;

use bagpiper::eqclass::{EqClass, Molecule};
use bagpiper::seq::{render_2bit, CellId, Umi};

fn transcripts() -> Vec<(String, u32)> {
    vec![
        ("T0".to_string(), 100),
        ("T1".to_string(), 200),
        ("T2".to_string(), 150),
    ]
}

/// (barcode, [(sorted transcript set, molecule multiplicity)]) per cell. Each multiplicity becomes
/// that many molecules with distinct UMIs. Shared classes carry unique anchors on both members.
fn spec() -> Vec<(&'static str, Vec<(Vec<u32>, u32)>)> {
    vec![
        (
            "AAAAAAAAAAAAAAAA",
            vec![(vec![0], 6), (vec![1], 4), (vec![0, 1], 5)],
        ),
        (
            "AAAAAAAAAAAAAAAC",
            vec![(vec![1], 3), (vec![2], 2), (vec![1, 2], 4), (vec![0], 1)],
        ),
    ]
}

fn build_eqclass() -> EqClass {
    let mut molecules = Vec::new();
    let mut umi_ix = 0u64;
    for (cb, groups) in spec() {
        let cell = CellId::from_ascii(cb.as_bytes()).unwrap();
        for (txps, mult) in groups {
            for _ in 0..mult {
                let umi = Umi::from_ascii(render_2bit(umi_ix, 12).as_bytes()).unwrap();
                umi_ix += 1;
                molecules.push(Molecule {
                    cell,
                    umi,
                    txps: txps.clone(),
                });
            }
        }
    }
    let mut eq = EqClass {
        transcripts: transcripts(),
        molecules,
    };
    bagpiper::dedup::exact(&mut eq.molecules); // cell-sort, as in the real pipeline
    eq
}

fn read_to_string(path: &Path, gz: bool) -> String {
    let f = std::fs::File::open(path).unwrap();
    let mut s = String::new();
    if gz {
        flate2::read::GzDecoder::new(f).read_to_string(&mut s).unwrap();
    } else {
        std::io::BufReader::new(f).read_to_string(&mut s).unwrap();
    }
    s
}

/// Read barcodes.tsv + matrix.mtx from a count output dir into a (barcode, tid) -> value map.
fn read_values(dir: &Path, gz: bool) -> HashMap<(String, u32), f32> {
    let ext = if gz { ".gz" } else { "" };
    let barcodes: Vec<String> = read_to_string(&dir.join(format!("barcodes.tsv{ext}")), gz)
        .lines()
        .map(str::to_string)
        .collect();
    let mtx = read_to_string(&dir.join(format!("matrix.mtx{ext}")), gz);
    let mut out = HashMap::new();
    for line in mtx.lines().skip(3) {
        let f: Vec<&str> = line.split('\t').collect();
        let row: usize = f[0].parse().unwrap();
        let tid: u32 = f[1].parse().unwrap();
        let val: f32 = f[2].parse().unwrap();
        out.insert((barcodes[row - 1].clone(), tid), val);
    }
    out
}

#[test]
fn multimapping_em_produces_fractional_abundances() {
    let eq = build_eqclass();
    let dir = std::env::temp_dir().join(format!("bp_mm_{}", std::process::id()));
    let _ = std::fs::remove_dir_all(&dir);
    std::fs::create_dir_all(&dir).unwrap();
    bagpiper::count::write_matrix(&eq, &dir).unwrap();

    let vals = read_values(&dir, true);
    assert!(
        vals.values().any(|v| v.fract().abs() > 1e-6),
        "multi-mapping EM must yield at least one fractional abundance: {vals:?}"
    );
    // Mass is conserved: total assigned ~= number of molecules (25), up to floored tiny alphas.
    let total: f32 = vals.values().sum();
    assert!((total - 25.0).abs() < 0.5, "total mass {total} should be ~25");

    let _ = std::fs::remove_dir_all(&dir);
}

#[test]
fn matches_reference_count_within_tolerance() {
    let ref_bin = match std::env::var("BP_TEST_REF_BIN") {
        Ok(b) => b,
        _ => return,
    };
    let eq = build_eqclass();
    let dir = std::env::temp_dir().join(format!("bp_mmdiff_{}", std::process::id()));
    let _ = std::fs::remove_dir_all(&dir);
    let mine = dir.join("mine");
    let reff = dir.join("ref");
    std::fs::create_dir_all(&mine).unwrap();
    std::fs::create_dir_all(&reff).unwrap();

    // my output (gzipped)
    bagpiper::count::write_matrix(&eq, &mine).unwrap();

    // the same molecules as a reference sorted_eqclass.txt (cell-sorted by dedup already)
    let eqtxt = dir.join("sorted_eqclass.txt");
    {
        let mut w = std::io::BufWriter::new(std::fs::File::create(&eqtxt).unwrap());
        writeln!(w, "{}", eq.transcripts.len()).unwrap();
        for (name, len) in &eq.transcripts {
            writeln!(w, "{name}\t{len}").unwrap();
        }
        for m in &eq.molecules {
            write!(w, "{}\t{}", m.cell.render(), m.umi.render()).unwrap();
            for &tid in &m.txps {
                write!(w, "\t{}", eq.transcripts[tid as usize].0).unwrap();
            }
            writeln!(w).unwrap();
        }
    }

    let status = std::process::Command::new(&ref_bin)
        .args(["count", "--eq"])
        .arg(&eqtxt)
        .arg("-o")
        .arg(&reff)
        .status()
        .expect("run reference count");
    assert!(status.success(), "reference count failed");

    let mine_v = read_values(&mine, true);
    let ref_v = read_values(&reff, false);

    let mut keys: std::collections::BTreeSet<(String, u32)> = std::collections::BTreeSet::new();
    keys.extend(mine_v.keys().cloned());
    keys.extend(ref_v.keys().cloned());
    let mut max_abs = 0.0f32;
    for k in &keys {
        let a = mine_v.get(k).copied().unwrap_or(0.0);
        let b = ref_v.get(k).copied().unwrap_or(0.0);
        max_abs = max_abs.max((a - b).abs());
    }
    eprintln!("multimap differential: {} positions, max_abs_diff = {max_abs}", keys.len());
    assert!(max_abs < 1e-2, "matrix must match reference within EM tolerance, got {max_abs}");

    let _ = std::fs::remove_dir_all(&dir);
}
