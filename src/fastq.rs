//! Gzipped FASTQ/FASTA I/O: needletail's parser fed a flate2 gzip stream for input, and a gzip FASTA
//! writer for the barcode/tso stages. needletail's own compression backends (bzip2/xz/zstd, which
//! compile C) stay disabled; PIP-seq data is gzip, decoded pure-Rust by flate2.

use std::fs::File;
use std::io::{self, BufReader, BufWriter, Write};
use std::path::Path;

use flate2::read::GzDecoder;
use flate2::write::GzEncoder;
use flate2::Compression;
use needletail::parser::{parse_fastx_reader, FastxReader};

/// Open a gzip-compressed FASTQ file for streaming records.
pub fn open_gz<P: AsRef<Path>>(path: P) -> io::Result<Box<dyn FastxReader>> {
    let gz = GzDecoder::new(BufReader::new(File::open(path)?));
    parse_fastx_reader(gz).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
}

/// The read id up to the first space or tab (the original read name, dropping any description).
pub fn read_name(id: &[u8]) -> &[u8] {
    id.split(|&b| b == b' ' || b == b'\t').next().unwrap_or(id)
}

/// A gzip writer over a new file for FASTA output (deterministic header: flate2 defaults mtime to 0).
pub fn gz_writer(path: std::path::PathBuf) -> io::Result<GzEncoder<BufWriter<File>>> {
    Ok(GzEncoder::new(
        BufWriter::new(File::create(path)?),
        Compression::default(),
    ))
}

/// Write one 2-line FASTA record (`>id\nseq\n`).
pub fn write_fasta<W: Write>(w: &mut W, id: &[u8], seq: &[u8]) -> io::Result<()> {
    w.write_all(b">")?;
    w.write_all(id)?;
    w.write_all(b"\n")?;
    w.write_all(seq)?;
    w.write_all(b"\n")
}

#[cfg(test)]
mod tests {
    use super::*;
    use flate2::read::GzDecoder;
    use std::io::Read;

    #[test]
    fn write_fasta_round_trips_through_gz() {
        let dir = std::env::temp_dir().join(format!("bp_fastq_{}", std::process::id()));
        let _ = std::fs::remove_dir_all(&dir);
        std::fs::create_dir_all(&dir).unwrap();
        let p = dir.join("out.fa.gz");

        let mut w = gz_writer(p.clone()).unwrap();
        write_fasta(&mut w, b"read1_CB_UMI", b"ACGTACGT").unwrap();
        write_fasta(&mut w, b"read2_CB_UMI", b"TTTT").unwrap();
        w.finish().unwrap();

        let mut s = String::new();
        GzDecoder::new(File::open(&p).unwrap())
            .read_to_string(&mut s)
            .unwrap();
        assert_eq!(s, ">read1_CB_UMI\nACGTACGT\n>read2_CB_UMI\nTTTT\n");
        let _ = std::fs::remove_dir_all(&dir);
    }
}
