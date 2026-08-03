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
        /// V5: run the guarded-BIN-id collapse (requires --t2g); implies packing the BIN-id key
        #[arg(long)]
        collapse: bool,
        /// transcript->gene TSV for --collapse
        #[arg(long)]
        t2g: Option<PathBuf>,
        /// --collapse guard: BIN-merge only where a cell-gene has <= T distinct 5' UMIs
        #[arg(long, default_value_t = 50)]
        collapse_t: usize,
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
                    std::io::Error::new(
                        std::io::ErrorKind::InvalidInput,
                        "illumina mode needs --r2",
                    )
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
            collapse,
            t2g,
            collapse_t,
        } => {
            std::fs::create_dir_all(&output)?;
            bagpiper::logging::init(&output, "count")?;
            let source = match (b1, r1, reference) {
                (Some(b1), None, None) => bagpiper::count::Source::Bam(b1),
                (None, Some(reads), Some(reference)) => {
                    bagpiper::count::Source::Align { reads, reference }
                }
                _ => {
                    return Err(std::io::Error::new(
                        std::io::ErrorKind::InvalidInput,
                        "count needs either --b1, or --r1 with --reference",
                    ))
                }
            };
            let collapse_cfg = if collapse {
                let t2g = t2g.ok_or_else(|| {
                    std::io::Error::new(
                        std::io::ErrorKind::InvalidInput,
                        "--collapse requires --t2g",
                    )
                })?;
                Some((t2g, collapse_t))
            } else {
                None
            };
            info!("count collapse={}", collapse);
            let c = bagpiper::count::run(source, collapse_cfg, &output)?;
            info!(
                "transcripts {}  molecules {} -> deduped {} -> final {}  output {}",
                c.transcripts,
                c.raw,
                c.deduped,
                c.final_molecules,
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
