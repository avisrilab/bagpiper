//! Real-data oracle: `read_bam` on a name-grouped BAM must reproduce the reference eqclass.
//! Parameterized by env so it is skipped when the data is not present:
//!   BP_TEST_BAM     path to the BAM
//!   BP_TEST_REF_EQ  path to the reference eqclass text (header + `CB\tUMI\ttxp...` body)

#[test]
fn bam_eqclass_matches_reference() {
    let (bam, ref_eq) = match (
        std::env::var("BP_TEST_BAM"),
        std::env::var("BP_TEST_REF_EQ"),
    ) {
        (Ok(a), Ok(b)) => (a, b),
        _ => return,
    };

    let eq = bagpiper::eqclass::read_bam(&bam, false).expect("read_bam");
    let mut buf = Vec::new();
    eq.write_text(&mut buf).unwrap();
    let got = String::from_utf8(buf).unwrap();
    if let Ok(out) = std::env::var("BP_TEST_OUT") {
        std::fs::write(&out, &got).unwrap();
    }
    let want = std::fs::read_to_string(&ref_eq).unwrap();

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
        "eqclass content must match the reference"
    );
}
