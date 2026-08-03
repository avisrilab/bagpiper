use std::path::PathBuf;

use clap::{Parser, Subcommand};
use log::info;

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
    /// Count matrix from a name-grouped BAM: exact-dedup molecules and run the per-cell EM, writing
    /// gzipped MatrixMarket + barcodes/features.
    Count {
        /// pre-aligned, name-grouped BAM
        #[arg(long)]
        b1: PathBuf,
        /// output directory
        #[arg(short, long)]
        output: PathBuf,
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
            bagpiper::logging::init(&output, "barcode")?;
            info!(
                "barcode {} r1={}",
                if nanopore { "nanopore" } else { "illumina" },
                r1.display()
            );
            let stats = if nanopore {
                bagpiper::barcode::run_nanopore(&r1, &whitelist, &output)?
            } else {
                let r2 = r2.ok_or_else(|| {
                    std::io::Error::new(std::io::ErrorKind::InvalidInput, "illumina mode needs --r2")
                })?;
                bagpiper::barcode::run_illumina(&r1, &r2, &whitelist, &output)?
            };
            info!(
                "total {}  matched {}  small {}  ambiguous {}  mismatch {}",
                stats.total, stats.matched, stats.small, stats.ambiguous, stats.mismatch
            );
        }
        Cmd::Count { b1, output } => {
            std::fs::create_dir_all(&output)?;
            bagpiper::logging::init(&output, "count")?;
            info!("count b1={}", b1.display());
            let mut eq = bagpiper::eqclass::read_bam(&b1)?;
            let raw = eq.molecules.len();
            bagpiper::dedup::exact(&mut eq.molecules);
            bagpiper::count::write_matrix(&eq, &output)?;
            info!(
                "transcripts {}  molecules {} -> deduped {}  output {}",
                eq.transcripts.len(),
                raw,
                eq.molecules.len(),
                output.display()
            );
        }
    }
    Ok(())
}
