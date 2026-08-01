//! Length-weighted EM abundance estimation.
//!
//! Ported from alevin-fry (https://github.com/COMBINE-lab/alevin-fry, src/em.rs), BSD-3-Clause.
//! Eqclasses are folded in sorted key order (BTreeMap) so the f32 sum is deterministic.

use std::collections::BTreeMap;

const MIN_OUTPUT_ALPHA: f32 = 0.01;
const ALPHA_CHECK_CUTOFF: f32 = 1e-2;
const MIN_ITER: u32 = 2;
const MAX_ITER: u32 = 100;
const REL_DIFF_TOLERANCE: f32 = 1e-2;

/// One EM update: split each multi-transcript eqclass's count across its transcripts by
/// length-weighted current abundance; a single-transcript eqclass adds its whole count.
fn update(
    alphas_in: &[f32],
    alphas_out: &mut [f32],
    txp_lengths: &[u32],
    eqclasses: &BTreeMap<Vec<u32>, u32>,
) {
    for (labels, &count) in eqclasses {
        if labels.len() > 1 {
            let min_len = labels
                .iter()
                .map(|&l| txp_lengths[l as usize])
                .min()
                .unwrap() as f32;
            let mut denominator = 0.0f32;
            let mut conditionals = Vec::with_capacity(labels.len());
            for &label in labels {
                let prob = alphas_in[label as usize] * min_len / txp_lengths[label as usize] as f32;
                denominator += prob;
                conditionals.push(prob);
            }
            if denominator > 0.0 {
                let inv_denominator = count as f32 / denominator;
                for (i, &label) in labels.iter().enumerate() {
                    alphas_out[label as usize] += inv_denominator * conditionals[i];
                }
            }
        } else {
            alphas_out[labels[0] as usize] += count as f32;
        }
    }
}

/// Estimate abundances for one cell's eqclasses, leaving the result in `alphas`. Both `alphas` and
/// `scratch` must be zeroed and sized to the transcriptome; `scratch` is reused working space.
/// Abundances below the output floor are set to zero.
pub fn optimize(
    eqclasses: &BTreeMap<Vec<u32>, u32>,
    txp_lengths: &[u32],
    alphas: &mut [f32],
    scratch: &mut [f32],
) {
    for (labels, &count) in eqclasses {
        if labels.len() == 1 {
            alphas[labels[0] as usize] += count as f32;
        }
    }
    alphas.iter_mut().for_each(|a| *a = (*a + 0.5) * 1e-3);

    let mut it_num: u32 = 0;
    let mut converged = true;
    while it_num < MIN_ITER || (it_num < MAX_ITER && !converged) {
        update(alphas, scratch, txp_lengths, eqclasses);

        converged = true;
        for i in 0..alphas.len() {
            if scratch[i] > ALPHA_CHECK_CUTOFF && (alphas[i] - scratch[i]).abs() > REL_DIFF_TOLERANCE {
                converged = false;
            }
            alphas[i] = scratch[i];
            scratch[i] = 0.0;
        }
        it_num += 1;
    }

    alphas.iter_mut().for_each(|a| {
        if *a < MIN_OUTPUT_ALPHA {
            *a = 0.0;
        }
    });
}

#[cfg(test)]
mod tests {
    use super::*;

    fn run(eq: &[(Vec<u32>, u32)], lengths: &[u32]) -> Vec<f32> {
        let map: BTreeMap<Vec<u32>, u32> = eq.iter().cloned().collect();
        let mut alphas = vec![0.0f32; lengths.len()];
        let mut scratch = vec![0.0f32; lengths.len()];
        optimize(&map, lengths, &mut alphas, &mut scratch);
        alphas
    }

    #[test]
    fn unique_counts_pass_through() {
        // Single-transcript eqclasses carry their whole count; EM does not move them.
        let a = run(&[(vec![0], 3), (vec![1], 5)], &[100, 100]);
        assert!((a[0] - 3.0).abs() < 1e-4, "{:?}", a);
        assert!((a[1] - 5.0).abs() < 1e-4, "{:?}", a);
    }

    #[test]
    fn equal_length_shared_eqclass_splits_evenly() {
        // One eqclass shared by two equal-length transcripts splits its count in half.
        let a = run(&[(vec![0, 1], 4)], &[100, 100]);
        assert!((a[0] - 2.0).abs() < 1e-4, "{:?}", a);
        assert!((a[1] - 2.0).abs() < 1e-4, "{:?}", a);
    }

    #[test]
    fn deterministic_across_runs() {
        // The sorted-key fold makes repeated runs bit-identical.
        let eq = [(vec![0, 1], 7), (vec![0], 2), (vec![1, 2], 5), (vec![2], 3)];
        let lengths = [120, 300, 210];
        assert_eq!(run(&eq, &lengths), run(&eq, &lengths));
    }
}
