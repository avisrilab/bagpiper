//! Smith-Waterman fitting alignment: place a pattern entirely into a substring of a read (pattern
//! global, read free-clipped at both ends) with affine gaps. Shared by the barcode seal
//! ([`crate::seal`]) and the V5 TSO seal ([`crate::tso`]). Scoring, gap model, and tie-break order
//! match rust-bio's `Aligner::custom` for these clip parameters (verified by `tests/seal_diff.rs`).

const MATCH: i32 = 2;
const MISMATCH: i32 = -1;
const GAP_OPEN: i32 = -1;
const GAP_EXTEND: i32 = -2;
const NEG_INF: i32 = i32::MIN / 4;

#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub(crate) enum Op {
    Match,
    Subst,
    Ins, // gap in the pattern: consume a read base
    Del, // gap in the read: consume a pattern base
}

/// Fit `pat` entirely into a substring of `read` (pattern global, read free-clipped both ends),
/// affine gaps. Returns (xstart, xend, ops from xstart to xend, score). None if the pattern cannot
/// fit. `score` is the alignment score under the shared scoring (V columns in the pattern never match
/// an ACGT read base, so they score as mismatches).
pub(crate) fn fit(read: &[u8], pat: &[u8]) -> Option<(usize, usize, Vec<Op>, i32)> {
    let (m, n) = (read.len(), pat.len());
    if n == 0 || m < n {
        return None;
    }
    // S: best ending aligned/clipped; ix: gap in pattern (consume read); dy: gap in read (consume pattern).
    let mut s = vec![vec![NEG_INF; n + 1]; m + 1];
    let mut ix = vec![vec![NEG_INF; n + 1]; m + 1];
    let mut dy = vec![vec![NEG_INF; n + 1]; m + 1];
    for row in s.iter_mut() {
        row[0] = 0; // free read-prefix clip
    }
    for j in 1..=n {
        let g = GAP_OPEN + GAP_EXTEND * j as i32;
        s[0][j] = g;
        dy[0][j] = g;
    }
    for i in 1..=m {
        for j in 1..=n {
            ix[i][j] = (ix[i - 1][j] + GAP_EXTEND).max(s[i - 1][j] + GAP_OPEN + GAP_EXTEND);
            dy[i][j] = (dy[i][j - 1] + GAP_EXTEND).max(s[i][j - 1] + GAP_OPEN + GAP_EXTEND);
            let diag = s[i - 1][j - 1]
                + if read[i - 1] == pat[j - 1] {
                    MATCH
                } else {
                    MISMATCH
                };
            s[i][j] = diag.max(ix[i][j]).max(dy[i][j]);
        }
    }

    // Best full-pattern alignment ends at the smallest row i maximizing s[i][n] (free read-suffix clip).
    let mut xend = 0;
    let mut best = NEG_INF;
    for i in 0..=m {
        if s[i][n] > best {
            best = s[i][n];
            xend = i;
        }
    }

    // Traceback, preferring match/subst > ins > del, and gap-open over gap-extend on ties (bio order).
    let mut ops = Vec::new();
    let (mut i, mut j) = (xend, n);
    #[derive(PartialEq)]
    enum Layer {
        S,
        I,
        D,
    }
    let mut layer = Layer::S;
    while j > 0 {
        match layer {
            Layer::S => {
                let diag = if i > 0 {
                    s[i - 1][j - 1]
                        + if read[i - 1] == pat[j - 1] {
                            MATCH
                        } else {
                            MISMATCH
                        }
                } else {
                    NEG_INF
                };
                if i > 0 && s[i][j] == diag {
                    ops.push(if read[i - 1] == pat[j - 1] {
                        Op::Match
                    } else {
                        Op::Subst
                    });
                    i -= 1;
                    j -= 1;
                } else if s[i][j] == ix[i][j] {
                    layer = Layer::I;
                } else {
                    layer = Layer::D;
                }
            }
            Layer::I => {
                ops.push(Op::Ins);
                if ix[i][j] == s[i - 1][j] + GAP_OPEN + GAP_EXTEND {
                    layer = Layer::S; // opened
                }
                i -= 1;
            }
            Layer::D => {
                ops.push(Op::Del);
                if dy[i][j] == s[i][j - 1] + GAP_OPEN + GAP_EXTEND {
                    layer = Layer::S; // opened
                }
                j -= 1;
            }
        }
    }
    let xstart = i;
    ops.reverse();
    Some((xstart, xend, ops, best))
}

/// Advance `offset` (a read index) past `xlim` pattern columns, following the alignment ops. A
/// deletion consumes a pattern column without a read base, an insertion a read base without a column.
pub(crate) fn advance(
    offset: &mut usize,
    xlim: usize,
    ops: &mut std::slice::Iter<Op>,
) -> Option<()> {
    let mut check = 0;
    loop {
        match ops.next()? {
            Op::Match | Op::Subst => check += 1,
            Op::Del => {
                *offset = offset.checked_sub(1)?;
                check += 1;
            }
            Op::Ins => *offset += 1,
        }
        if check == xlim {
            *offset += xlim;
            return Some(());
        }
    }
}
