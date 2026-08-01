//! Real-data oracle: read_bam + exact dedup + count must produce a well-formed gzipped matrix
//! whose declared nnz is the true triplet count. Values are checked against the reference matrix
//! out-of-band (the EM's f32 fold differs by summation order run to run). Parameterized by env so
//! it is skipped when the data is not present:
//!   BP_TEST_BAM       path to the name-grouped BAM
//!   BP_TEST_OUT_DIR   directory to write barcodes.tsv.gz / features.tsv.gz / matrix.mtx.gz into

use std::io::Read;

fn read_gz(path: String) -> String {
    let mut s = String::new();
    flate2::read::GzDecoder::new(std::fs::File::open(path).unwrap())
        .read_to_string(&mut s)
        .unwrap();
    s
}

#[test]
fn count_matrix_is_well_formed() {
    let (bam, out_dir) = match (
        std::env::var("BP_TEST_BAM"),
        std::env::var("BP_TEST_OUT_DIR"),
    ) {
        (Ok(a), Ok(b)) => (a, b),
        _ => return,
    };
    std::fs::create_dir_all(&out_dir).unwrap();

    let mut eq = bagpiper::eqclass::read_bam(&bam).expect("read_bam");
    let n_txp = eq.transcripts.len();
    bagpiper::dedup::exact(&mut eq.molecules);
    bagpiper::count::write_matrix(&eq, &out_dir).expect("write_matrix");

    let barcodes = read_gz(format!("{out_dir}/barcodes.tsv.gz"));
    let n_cells = barcodes.lines().count();
    let features = read_gz(format!("{out_dir}/features.tsv.gz"));
    assert_eq!(features.lines().count(), n_txp, "one feature per transcript");

    let mtx = read_gz(format!("{out_dir}/matrix.mtx.gz"));
    let lines: Vec<&str> = mtx.lines().collect();
    assert_eq!(lines[0], "%%MatrixMarket matrix coordinate real general");
    let dims: Vec<&str> = lines[2].split('\t').collect();
    assert_eq!(dims[0].parse::<usize>().unwrap(), n_cells, "rows == cells");
    assert_eq!(dims[1].parse::<usize>().unwrap(), n_txp, "cols == transcripts");
    let declared_nnz: usize = dims[2].parse().unwrap();
    assert_eq!(declared_nnz, lines.len() - 3, "declared nnz == actual triplet lines (true nnz)");

    // Every triplet is in range and carries a positive value.
    for line in &lines[3..] {
        let f: Vec<&str> = line.split('\t').collect();
        let r: usize = f[0].parse().unwrap();
        let c: usize = f[1].parse().unwrap();
        let v: f32 = f[2].parse().unwrap();
        assert!(r >= 1 && r <= n_cells, "row in range");
        assert!(c >= 1 && c <= n_txp, "col in range");
        assert!(v > 0.0, "value positive");
    }
}
