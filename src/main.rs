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
    /// Count matrix from aligned reads: either a pre-aligned name-grouped BAM (`--b1`) or barcoded
    /// reads aligned internally against a transcriptome (`--r1` + `--reference`). Exact-dedups
    /// molecules, runs the per-cell EM, and writes gzipped MatrixMarket + barcodes/features.
    Count {
        /// pre-aligned, name-grouped BAM (alternative to --r1 + --reference)
        #[arg(long)]
        b1: Option<PathBuf>,
        /// barcoded reads (gzipped FASTA) to align internally; requires --reference
        #[arg(long)]
        r1: Option<PathBuf>,
        /// transcriptome FASTA for internal alignment; requires --r1
        #[arg(long)]
        reference: Option<PathBuf>,
        /// output directory
        #[arg(short, long)]
        output: PathBuf,
        /// V5: pack the 3 bp binning index onto the molecular key (for the opt-in collapse stage)
        #[arg(long)]
        v5_binid: bool,
    },
    /// Extract the V5 5' dual-UMI (TSO seal) from barcoded reads, writing
    /// `>origid_CB_UMI_umi1_umi2` + trimmed cDNA. Run after `barcode`, before alignment, for V5.
    Tso {
        /// barcoded reads (gzipped FASTA) from `barcode`
        #[arg(long)]
        r1: PathBuf,
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
        Cmd::Count {
            b1,
            r1,
            reference,
            output,
            v5_binid,
        } => {
            std::fs::create_dir_all(&output)?;
            bagpiper::logging::init(&output, "count")?;
            let mut eq = match (b1, r1, reference) {
                (Some(b1), None, None) => {
                    info!("count b1={} v5_binid={}", b1.display(), v5_binid);
                    bagpiper::eqclass::read_bam(&b1, v5_binid)?
                }
                (None, Some(r1), Some(reference)) => {
                    info!("count align r1={} ref={}", r1.display(), reference.display());
                    bagpiper::align::align_to_eqclass(&r1, &reference, v5_binid)?
                }
                _ => {
                    return Err(std::io::Error::new(
                        std::io::ErrorKind::InvalidInput,
                        "count needs either --b1, or --r1 with --reference",
                    ))
                }
            };
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
        Cmd::Tso { r1, output } => {
            std::fs::create_dir_all(&output)?;
            bagpiper::logging::init(&output, "tso")?;
            info!("tso r1={}", r1.display());
            let stats = bagpiper::tso::run_tso(&r1, &output)?;
            info!(
                "total {}  matched {}  small {}  no_seal {}  far_seal {}",
                stats.total, stats.matched, stats.small, stats.no_seal, stats.far_seal
            );
        }
    }
    Ok(())
}
