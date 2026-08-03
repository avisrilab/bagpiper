//! Opt-in guarded binning-index molecule collapse for the V5 TSO-UMI chemistry.
//!
//! The 5' TSO UMI over-splits a true molecule when the 5' seal mis-anchors on a noisy long-read
//! start, inflating the molecule count. The 3 bp binning index (packed onto the key by `eqclass
//! --v5-binid`) is an independent weak per-molecule tag from a cleaner part of the read; this stage
//! uses it to merge over-split 5' UMIs, GUARDED so it fires only where it is safe:
//!
//!   per (cell barcode, gene): let n5 = number of distinct 5' UMIs
//!     n5 <= T  ->  molecules = groups of reads sharing a BIN-id  (rescue the over-split)
//!     n5 >  T  ->  molecules = distinct 5' UMIs                   (fall back; the 3 bp BIN-id
//!                                                                  collides at high expression)
//!
//! The guard is per-GENE because collision risk on the small BIN-id space is set by how many
//! molecules of that gene sit in the cell. Merges never cross a gene, so within-gene isoform
//! proportions are preserved; only absolute molecule counts change, which is why it is opt-in.
//!
//! Operates in memory on the cell-sorted molecules (see [`crate::dedup::exact`]): the 15 nt
//! `--v5-binid` UMI is sliced into `umi5[..12]` + `binid[12..15]`, and a molecule's gene comes from
//! the transcript->gene map applied to its transcript set. Each merged group becomes one molecule
//! carrying the group's dominant transcript set, which `count` then keys on.

use std::collections::HashMap;
use std::fs::File;
use std::io::{self, BufRead, BufReader};
use std::path::Path;

use crate::eqclass::{EqClass, Molecule};
use crate::seq::CellId;

const UMI5_LEN: usize = 12;
const PACKED_KEY_LEN: usize = 15;

/// Load a transcript->gene map (TSV `transcript_id <tab> gene_id ...`; a header row whose first
/// field is `transcript_id` and any extra columns are ignored). Transcripts absent from the map fall
/// back to their own name as the gene, so a partial map degrades gracefully.
pub fn load_t2g<P: AsRef<Path>>(path: P) -> io::Result<HashMap<String, String>> {
    let mut map = HashMap::new();
    for line in BufReader::new(File::open(path)?).lines() {
        let line = line?;
        let mut it = line.split_whitespace();
        if let (Some(tx), Some(gene)) = (it.next(), it.next()) {
            if tx == "transcript_id" {
                continue;
            }
            map.insert(tx.to_string(), gene.to_string());
        }
    }
    Ok(map)
}

/// Resolve the gene of a transcript set: the shared gene if every transcript maps to the same one,
/// else `None` (multi-gene, never BIN-merged). Unknown transcripts map to their own name.
fn resolve_gene(
    txps: &[u32],
    transcripts: &[(String, u32)],
    t2g: &HashMap<String, String>,
) -> Option<String> {
    let mut g: Option<&str> = None;
    for &t in txps {
        let name = transcripts[t as usize].0.as_str();
        let gene = t2g.get(name).map(String::as_str).unwrap_or(name);
        match g {
            None => g = Some(gene),
            Some(prev) if prev == gene => {}
            Some(_) => return None,
        }
    }
    g.map(str::to_string)
}

/// The dominant transcript set among a group (indices into `mols`): most frequent, ties broken by
/// the lexicographically smallest set (deterministic), matching the reference.
fn dominant_txps(group: &[usize], mols: &[Molecule]) -> Vec<u32> {
    let mut counts: HashMap<&[u32], usize> = HashMap::new();
    for &i in group {
        *counts.entry(mols[i].txps.as_slice()).or_insert(0) += 1;
    }
    counts
        .into_iter()
        .max_by(|a, b| a.1.cmp(&b.1).then_with(|| b.0.cmp(a.0)))
        .map(|(txps, _)| txps.to_vec())
        .unwrap()
}

/// Emit one molecule per group (groups keyed by umi5 or binid), in sorted key order for determinism.
/// The molecule carries the group's dominant transcript set; its UMI is the first member's (count
/// keys on the transcript set, not the UMI, so the representative only needs to be valid).
fn emit_groups(
    groups: &HashMap<&str, Vec<usize>>,
    cell_mols: &[Molecule],
    out: &mut Vec<Molecule>,
) {
    let mut keys: Vec<&str> = groups.keys().copied().collect();
    keys.sort_unstable();
    for k in keys {
        let group = &groups[k];
        out.push(Molecule {
            cell: cell_mols[group[0]].cell,
            umi: cell_mols[group[0]].umi,
            txps: dominant_txps(group, cell_mols),
        });
    }
}

/// Collapse the molecules of ONE cell under the guarded-BIN-id rule, appending merged molecules.
fn collapse_cell(
    cell_mols: &[Molecule],
    transcripts: &[(String, u32)],
    t2g: &HashMap<String, String>,
    threshold: usize,
    out: &mut Vec<Molecule>,
) {
    let rendered: Vec<String> = cell_mols.iter().map(|m| m.umi.render()).collect();
    let umi5 = |i: usize| &rendered[i][..UMI5_LEN.min(rendered[i].len())];
    let binid = |i: usize| {
        let r = &rendered[i];
        if r.len() >= PACKED_KEY_LEN {
            &r[UMI5_LEN..PACKED_KEY_LEN]
        } else {
            ""
        }
    };

    // Partition by gene; multi-gene molecules never BIN-merge.
    let mut by_gene: HashMap<String, Vec<usize>> = HashMap::new();
    let mut multi_gene: Vec<usize> = Vec::new();
    for (i, m) in cell_mols.iter().enumerate() {
        match resolve_gene(&m.txps, transcripts, t2g) {
            Some(g) => by_gene.entry(g).or_default().push(i),
            None => multi_gene.push(i),
        }
    }

    let mut genes: Vec<&String> = by_gene.keys().collect();
    genes.sort_unstable();
    for gene in genes {
        let idxs = &by_gene[gene];
        let mut umi5_groups: HashMap<&str, Vec<usize>> = HashMap::new();
        for &i in idxs {
            umi5_groups.entry(umi5(i)).or_default().push(i);
        }
        if umi5_groups.len() <= threshold {
            let mut bin_groups: HashMap<&str, Vec<usize>> = HashMap::new();
            for &i in idxs {
                bin_groups.entry(binid(i)).or_default().push(i);
            }
            emit_groups(&bin_groups, cell_mols, out);
        } else {
            emit_groups(&umi5_groups, cell_mols, out);
        }
    }

    if !multi_gene.is_empty() {
        let mut umi5_groups: HashMap<&str, Vec<usize>> = HashMap::new();
        for &i in &multi_gene {
            umi5_groups.entry(umi5(i)).or_default().push(i);
        }
        emit_groups(&umi5_groups, cell_mols, out);
    }
}

/// Guarded-BIN-id collapse over the whole eqclass, in place. `eq.molecules` must be cell-sorted
/// (dedup leaves them so); the molecules must carry the 15 nt `--v5-binid` key. Streams one cell-run
/// at a time, so peak extra memory is the collapsed output (smaller than the input) plus the map.
pub fn collapse(eq: &mut EqClass, t2g: &HashMap<String, String>, threshold: usize) {
    let mols = std::mem::take(&mut eq.molecules);
    let mut out: Vec<Molecule> = Vec::new();
    let mut i = 0;
    while i < mols.len() {
        let cell: CellId = mols[i].cell;
        let start = i;
        while i < mols.len() && mols[i].cell == cell {
            i += 1;
        }
        collapse_cell(&mols[start..i], &eq.transcripts, t2g, threshold, &mut out);
    }
    eq.molecules = out;
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::seq::{CellId, Umi};

    fn t2g(pairs: &[(&str, &str)]) -> HashMap<String, String> {
        pairs
            .iter()
            .map(|(t, g)| (t.to_string(), g.to_string()))
            .collect()
    }

    // A molecule with a 15 nt (umi5 + binid) key and a single transcript id.
    fn mol(umi5: &str, binid: &str, tid: u32) -> Molecule {
        Molecule {
            cell: CellId(0),
            umi: Umi::from_ascii(format!("{umi5}{binid}").as_bytes()).unwrap(),
            txps: vec![tid],
        }
    }

    fn eqc(mols: Vec<Molecule>) -> EqClass {
        EqClass {
            // tid 0 -> TXP0 (gene GENE_A), tid 1 -> TXP1 (gene GENE_B)
            transcripts: vec![("TXP0".into(), 100), ("TXP1".into(), 100)],
            molecules: mols,
        }
    }

    fn run(mols: Vec<Molecule>, threshold: usize) -> usize {
        let mut eq = eqc(mols);
        collapse(
            &mut eq,
            &t2g(&[("TXP0", "GENE_A"), ("TXP1", "GENE_B")]),
            threshold,
        );
        eq.molecules.len()
    }

    #[test]
    fn low_expression_merges_over_split_umis_sharing_a_binid() {
        // One gene, three reads: two distinct 5' UMIs but the same BIN-id -> one molecule at T>=2.
        let mols = vec![
            mol("AAAAAACCCCCC", "ACG", 0),
            mol("GGGGGGTTTTTT", "ACG", 0), // over-split 5' UMI, same BIN-id
            mol("AAAAAACCCCCC", "TGC", 0), // different BIN-id -> its own molecule
        ];
        // n5 = 2 <= T: merge by BIN-id -> {ACG, TGC} = 2 molecules.
        assert_eq!(run(mols, 50), 2);
    }

    #[test]
    fn high_expression_falls_back_to_the_5prime_umi() {
        // One gene, two reads with distinct 5' UMIs but the same BIN-id; T=1 forces fallback.
        let mols = vec![mol("AAAAAACCCCCC", "ACG", 0), mol("GGGGGGTTTTTT", "ACG", 0)];
        // n5 = 2 > T=1: keep the 5' UMI -> 2 molecules (no BIN merge).
        assert_eq!(run(mols, 1), 2);
    }

    #[test]
    fn merges_never_cross_a_gene() {
        // Same 5' UMI and BIN-id, but different genes: two molecules, never merged.
        let mols = vec![mol("AAAAAACCCCCC", "ACG", 0), mol("AAAAAACCCCCC", "ACG", 1)];
        assert_eq!(run(mols, 50), 2);
    }
}
