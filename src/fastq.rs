//! Gzipped FASTQ/FASTA I/O: needletail's parser fed a flate2 gzip stream for input, and a gzip FASTA
//! writer for the barcode/tso stages. needletail's own compression backends (bzip2/xz/zstd, which
//! compile C) stay disabled; PIP-seq data is gzip, decoded pure-Rust by flate2.

use std::fs::File;
use std::io::{self, BufReader, BufWriter, Write};
use std::path::{Path, PathBuf};

use flate2::read::MultiGzDecoder;
use flate2::write::GzEncoder;
use flate2::Compression;
use needletail::parser::{parse_fastx_reader, FastxReader};

/// Open a gzip-compressed FASTQ file for streaming records. Uses `MultiGzDecoder`, so a
/// multi-member gzip (e.g. two flow-cell fastqs concatenated, or dorado's multi-member output) is
/// read in full rather than truncated at the first member.
pub fn open_gz<P: AsRef<Path>>(path: P) -> io::Result<Box<dyn FastxReader>> {
    let gz = MultiGzDecoder::new(BufReader::new(File::open(path)?));
    parse_fastx_reader(gz).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
}

/// Streams records from one or more gzip files in sequence, so a library split across several fastqs
/// (e.g. two flow-cell files) is read as one stream. Yields owned `(id, seq)` so the per-file reader
/// can be swapped at a record boundary without returning a borrow across it.
pub struct MultiReader {
    paths: std::vec::IntoIter<PathBuf>,
    cur: Option<Box<dyn FastxReader>>,
}

impl MultiReader {
    pub fn open(paths: &[PathBuf]) -> MultiReader {
        MultiReader {
            paths: paths.to_vec().into_iter(),
            cur: None,
        }
    }

    /// The next record as owned `(id, seq)`, advancing to the next file at end-of-file. None once
    /// every file is exhausted; an open or parse error is surfaced as the item.
    pub fn next_seq(&mut self) -> Option<io::Result<(Vec<u8>, Vec<u8>)>> {
        loop {
            if self.cur.is_none() {
                match open_gz(self.paths.next()?) {
                    Ok(r) => self.cur = Some(r),
                    Err(e) => return Some(Err(e)),
                }
            }
            let exhausted = match self.cur.as_mut().unwrap().next() {
                Some(Ok(r)) => return Some(Ok((r.id().to_vec(), r.seq().to_vec()))),
                Some(Err(e)) => return Some(Err(io::Error::new(io::ErrorKind::InvalidData, e))),
                None => true,
            };
            if exhausted {
                self.cur = None;
            }
        }
    }
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

    #[test]
    fn open_gz_reads_across_concatenated_gzip_members() {
        // A flow cell can split one library into two fastqs; concatenated they are a multi-member
        // gzip. open_gz must read across the boundary (MultiGzDecoder), not stop at the first member
        // as a plain GzDecoder would (which would silently drop the second file's reads).
        let dir = std::env::temp_dir().join(format!("bp_multigz_{}", std::process::id()));
        let _ = std::fs::remove_dir_all(&dir);
        std::fs::create_dir_all(&dir).unwrap();

        let write_member = |name: &str, ids: &[&[u8]]| {
            let p = dir.join(name);
            let mut w = gz_writer(p.clone()).unwrap();
            for id in ids {
                write_fasta(&mut w, id, b"ACGT").unwrap();
            }
            w.finish().unwrap();
            std::fs::read(p).unwrap()
        };
        let mut bytes = write_member("a.fa.gz", &[b"r1", b"r2"]);
        bytes.extend(write_member("b.fa.gz", &[b"r3", b"r4"]));
        let cat = dir.join("cat.fa.gz");
        std::fs::write(&cat, &bytes).unwrap();

        let mut reader = open_gz(&cat).unwrap();
        let mut n = 0;
        while let Some(rec) = reader.next() {
            rec.unwrap();
            n += 1;
        }
        assert_eq!(n, 4, "must read all 4 records across the two gzip members");
        let _ = std::fs::remove_dir_all(&dir);
    }

    #[test]
    fn multireader_chains_files_in_order() {
        // A split library (two separate fastqs) must read as one stream, in file then record order.
        let dir = std::env::temp_dir().join(format!("bp_multi_{}", std::process::id()));
        let _ = std::fs::remove_dir_all(&dir);
        std::fs::create_dir_all(&dir).unwrap();
        let write = |name: &str, ids: &[&[u8]]| {
            let p = dir.join(name);
            let mut w = gz_writer(p.clone()).unwrap();
            for id in ids {
                write_fasta(&mut w, id, b"ACGT").unwrap();
            }
            w.finish().unwrap();
            p
        };
        let a = write("a.fa.gz", &[b"r1", b"r2"]);
        let b = write("b.fa.gz", &[b"r3", b"r4"]);

        let mut mr = MultiReader::open(&[a.clone(), b]);
        let mut ids = Vec::new();
        while let Some(res) = mr.next_seq() {
            ids.push(String::from_utf8(res.unwrap().0).unwrap());
        }
        assert_eq!(
            ids,
            vec!["r1", "r2", "r3", "r4"],
            "chained in file then record order"
        );

        // one file: same as reading that file alone
        let mut one = MultiReader::open(std::slice::from_ref(&a));
        let mut n = 0;
        while let Some(res) = one.next_seq() {
            res.unwrap();
            n += 1;
        }
        assert_eq!(n, 2);

        // no files: nothing to read
        assert!(MultiReader::open(&[]).next_seq().is_none());
        let _ = std::fs::remove_dir_all(&dir);
    }
}
