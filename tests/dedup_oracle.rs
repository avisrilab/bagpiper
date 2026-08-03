//! Real-data oracle: `read_bam` + exact dedup must reproduce the reference `sorted_eqclass.txt` as
//! the same set of unique molecule lines. Line order may differ (the reference sorts lines as
//! strings, the rewrite by packed key); both are valid intermediates under the count-matrix
//! contract. Parameterized by env so it is skipped when the data is not present:
//!   BP_TEST_BAM         path to the name-grouped BAM
//!   BP_TEST_REF_SORTED  path to the reference sorted_eqclass.txt (header + deduped body)
//!   BP_TEST_OUT         optional: write the rendered output here for an out-of-band diff

#[test]
fn dedup_matches_reference_sorted() {
    let (bam, ref_sorted) = match (
        std::env::var("BP_TEST_BAM"),
        std::env::var("BP_TEST_REF_SORTED"),
    ) {
        (Ok(a), Ok(b)) => (a, b),
        _ => return,
    };

    let mut eq = bagpiper::eqclass::read_bam(&bam, false).expect("read_bam");
    bagpiper::dedup::exact(&mut eq.molecules);
    let mut buf = Vec::new();
    eq.write_text(&mut buf).unwrap();
    let got = String::from_utf8(buf).unwrap();
    if let Ok(out) = std::env::var("BP_TEST_OUT") {
        std::fs::write(&out, &got).unwrap();
    }
    let want = std::fs::read_to_string(&ref_sorted).unwrap();

    let mut got_lines: Vec<&str> = got.lines().collect();
    let mut want_lines: Vec<&str> = want.lines().collect();
    assert_eq!(
        got_lines.len(),
        want_lines.len(),
        "line count differs (got {}, want {})",
        got_lines.len(),
        want_lines.len()
    );
    got_lines.sort_unstable();
    want_lines.sort_unstable();
    assert_eq!(
        got_lines, want_lines,
        "deduped molecule set must match the reference"
    );
}
