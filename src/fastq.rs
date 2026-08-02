//! Gzipped FASTQ input: needletail's parser fed a flate2 gzip stream. needletail's own compression
//! backends (bzip2/xz/zstd, which compile C) stay disabled; PIP-seq data is gzip, decoded pure-Rust
//! by flate2.

use std::fs::File;
use std::io::{self, BufReader};
use std::path::Path;

use flate2::read::GzDecoder;
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
