//! bagpiper: PIP-seq barcode assignment and equivalence-class quantification.

// Core primitives
pub mod parallel; // reusable feeder-consumer
pub mod seq; // 2-bit packed CellId / Umi
pub mod whitelist; // exact + edit-1 barcode correction

// I/O
pub mod fastq;
pub mod logging;

// Chemistry (read layouts behind a seam)
pub mod chemistry;

// Pipeline stages
pub mod barcode; // FASTQ -> assigned CB/UMI + cDNA
pub mod seal; // Smith-Waterman barcode rescue
pub mod eqclass; // BAM -> packed molecules
pub mod dedup; // exact PCR dedup
pub mod em; // length-weighted EM
pub mod count; // molecules -> count matrix
