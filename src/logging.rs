//! Unified run log. Terminal output is ephemeral, so every subcommand tees its log records to both
//! stderr and a per-run file `<output_dir>/bagpiper.<subcommand>.log` in the run's own output folder.
//! One logger, set up once per run from `main`.

use std::fs::File;
use std::io::{self, BufWriter, Write};
use std::path::Path;
use std::sync::Mutex;
use std::time::SystemTime;

use log::{LevelFilter, Log, Metadata, Record};

struct TeeLogger {
    file: Mutex<BufWriter<File>>,
    level: LevelFilter,
}

impl Log for TeeLogger {
    fn enabled(&self, metadata: &Metadata) -> bool {
        metadata.level() <= self.level
    }

    fn log(&self, record: &Record) {
        if !self.enabled(record.metadata()) {
            return;
        }
        let ts = humantime::format_rfc3339_seconds(SystemTime::now());
        let line = format!("{} [{}] {}", ts, record.level(), record.args());
        eprintln!("{}", line);
        if let Ok(mut f) = self.file.lock() {
            let _ = writeln!(f, "{}", line);
            let _ = f.flush();
        }
    }

    fn flush(&self) {
        if let Ok(mut f) = self.file.lock() {
            let _ = f.flush();
        }
    }
}

/// Initialize logging for a run: records at INFO and above go to stderr and to
/// `<output_dir>/bagpiper.<subcommand>.log`. Call once, after the output folder exists; the file
/// mtime records when the run happened. Data outputs are untouched, so the log is additive only.
pub fn init(output_dir: &Path, subcommand: &str) -> io::Result<()> {
    let level = LevelFilter::Info;
    let path = output_dir.join(format!("bagpiper.{}.log", subcommand));
    let file = File::create(&path)?;
    let logger = TeeLogger {
        file: Mutex::new(BufWriter::new(file)),
        level,
    };
    log::set_boxed_logger(Box::new(logger))
        .map_err(|e| io::Error::new(io::ErrorKind::Other, e.to_string()))?;
    log::set_max_level(level);
    Ok(())
}
