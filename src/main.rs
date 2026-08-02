use std::path::PathBuf;

use clap::{Parser, Subcommand};

#[derive(Parser)]
#[command(name = "bagpiper", version, about = "Process PIP-seq data.")]
struct Cli {
    #[command(subcommand)]
    command: Cmd,
}

#[derive(Subcommand)]
enum Cmd {
    /// Assign cell barcodes and UMIs from FASTQ, writing barcoded reads as gzipped FASTA.
    Barcode {
        /// input read 1 FASTQ (.fq.gz)
        #[arg(long)]
        r1: PathBuf,
        /// input read 2 FASTQ (.fq.gz); illumina only
        #[arg(long)]
        r2: Option<PathBuf>,
        /// barcode whitelist CSV
        #[arg(long)]
        whitelist: PathBuf,
        /// output directory
        #[arg(short, long)]
        output: PathBuf,
        /// nanopore (long-read) mode: the reverse/forward regex + Smith-Waterman seal cascade
        #[arg(long)]
        nanopore: bool,
    },
}

fn main() -> std::io::Result<()> {
    match Cli::parse().command {
        Cmd::Barcode {
            r1,
            r2,
            whitelist,
            output,
            nanopore,
        } => {
            std::fs::create_dir_all(&output)?;
            let stats = if nanopore {
                bagpiper::barcode::run_nanopore(&r1, &whitelist, &output)?
            } else {
                let r2 = r2.ok_or_else(|| {
                    std::io::Error::new(std::io::ErrorKind::InvalidInput, "illumina mode needs --r2")
                })?;
                bagpiper::barcode::run_illumina(&r1, &r2, &whitelist, &output)?
            };
            eprintln!(
                "Total: {}  Matched: {}  Small: {}  Ambiguous: {}  Mismatch: {}",
                stats.total, stats.matched, stats.small, stats.ambiguous, stats.mismatch
            );
        }
    }
    Ok(())
}
