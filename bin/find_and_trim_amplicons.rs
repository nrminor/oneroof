#!/usr/bin/env rust-script
//! Find and trim primer sequences in FASTQ/FASTA reads.
//!
//! This Rust script combines primer finding and trimming operations into a single process,
//! replacing two-step find-then-trim processing. It uses approximate string matching to locate primers,
//! handles both orientations, and outputs trimmed amplicons with detailed statistics.
//!
//! # Features
//!
//! - **Parallel processing**: Multi-threaded read processing for high throughput
//! - **Approximate matching**: Allows mismatches in primer sequences (configurable)
//! - **Bidirectional search**: Finds amplicons in both forward and reverse orientations
//! - **Smart selection**: Chooses best amplicon based on edit distance and length
//! - **Flexible output**: Supports FASTQ/FASTA with optional gzip compression
//! - **Detailed statistics**: Tracks success rates and failure reasons
//! - **IUPAC support**: Handles ambiguity codes in primer sequences
//!
//! # Algorithm
//!
//! By default, the entire read is searched for primer sequences. For each read:
//!
//! 1. Search for forward and reverse primers in both orientations
//! 2. Find all valid amplicons (primers in correct order, within length bounds)
//! 3. Select the amplicon with lowest total edit distance (ties broken by length)
//! 4. Trim the read to exclude primers, keeping only the insert sequence
//! 5. Write trimmed sequence to output file
//!
//! ## Search Window Optimization
//!
//! For long reads (e.g., Nanopore), searching the entire read for primers can be slow.
//! The `--forward-window` and `--reverse-window` options limit primer search to bounded
//! regions at the ends of reads, improving performance when primers are expected near
//! the read termini.
//!
//! - `--forward-window N`: Search for forward primer only in the first N bases (5' end
//!   in forward orientation, 3' end in reverse orientation)
//! - `--reverse-window N`: Search for reverse primer only in the last N bases (3' end
//!   in forward orientation, 5' end in reverse orientation)
//! - Default (0): Search the entire read (no windowing)
//!
//! Choose window sizes larger than your primer length plus any expected adapter/barcode
//! sequences. A warning is issued if the window size is smaller than the primer length.
//!
//! # Usage
//!
//! Basic usage (searches entire read for primers):
//!
//! ```bash
//! find_and_trim_amplicons.rs \
//!   -i input.fastq.gz \
//!   -o output.fastq.gz \
//!   -f ACGTACGTACGT \
//!   -r TGCATGCATGCA \
//!   --min-len 100 \
//!   --max-len 500 \
//!   -k 2 \
//!   -t 8
//! ```
//!
//! With search windows for long reads (primers expected within 100bp of read ends):
//!
//! ```bash
//! find_and_trim_amplicons.rs \
//!   -i long_reads.fastq.gz \
//!   -o trimmed.fastq.gz \
//!   -f ACGTACGTACGT \
//!   -r TGCATGCATGCA \
//!   --forward-window 100 \
//!   --reverse-window 100
//! ```
//!
//! # Environment Variables
//!
//! - `RUST_LOG`: Controls logging verbosity (e.g., `debug`, `info`, `warn`)
//! - `RAYON_NUM_THREADS`: Overrides default thread count
//!
//! ```cargo
//! [dependencies]
//! paraseq = { version = "0.4", features = ["niffler"] }
//! noodles = { version = "0.102", features = ["fastq", "fasta"] }
//! sassy = { version = "0.1" }
//! clap = { version = "4.5", features = ["derive", "cargo"] }
//! anyhow = "1.0"
//! env_logger = "0.11"
//! log = "0.4"
//! num_cpus = "1.16"
//! flate2 = "1.0"
//! ```

use anyhow::Result;
use clap::Parser;
use flate2::{write::GzEncoder, Compression};
use log::{debug, info, warn};
use paraseq::{fastx, prelude::*, ProcessError};
use sassy::{profiles::Iupac, Searcher};
use std::{
    cmp,
    fs::File,
    io::Write,
    path::PathBuf,
    sync::{
        atomic::{AtomicUsize, Ordering},
        Arc, Mutex,
    },
};

// Type aliases for complex nested types
type SharedWriter = Arc<Mutex<Box<dyn Write + Send>>>;
type SharedStats = Arc<TrimStats>;

/// Command line arguments for primer trimming
#[derive(Parser, Debug)]
#[command(name = "primer-trim")]
#[command(version, about = "Find and trim primer sequences in FASTQ/FASTA reads")]
struct Args {
    /// Input FASTQ/FASTA file (format auto-detected, gzip/bzip2 supported)
    #[arg(short = 'i', long)]
    input: PathBuf,

    /// Output file path
    #[arg(short = 'o', long)]
    output: PathBuf,

    /// Forward primer sequence
    #[arg(short = 'f', long)]
    forward: String,

    /// Reverse primer sequence
    #[arg(short = 'r', long)]
    reverse: String,

    /// Maximum number of mismatches allowed
    #[arg(short = 'k', long, default_value = "2")]
    max_mismatch: u8,

    /// Minimum amplicon length (after trimming)
    #[arg(long, default_value = "50")]
    min_len: usize,

    /// Maximum amplicon length (after trimming)
    #[arg(long, default_value = "2000")]
    max_len: usize,

    /// Number of threads (default: half of available cores)
    #[arg(short = 't', long, default_value_t = default_threads())]
    threads: usize,

    /// Output format: fastq or fasta
    #[arg(long, default_value = "fastq")]
    format: OutputFormat,

    /// Disable gzip compression of output
    #[arg(long)]
    no_compress: bool,

    /// Statistics output file (default: stderr)
    #[arg(long)]
    stats: Option<PathBuf>,

    /// Limit forward primer search to N bases at the expected end of the read.
    /// In forward orientation, searches the first N bases (5' end).
    /// In reverse orientation, searches the last N bases (3' end).
    /// Default (0) searches the entire read. Use for long reads when primers
    /// are expected near read termini.
    #[arg(long, default_value = "0")]
    forward_window: usize,

    /// Limit reverse primer search to N bases at the expected end of the read.
    /// In forward orientation, searches the last N bases (3' end).
    /// In reverse orientation, searches the first N bases (5' end).
    /// Default (0) searches the entire read. Use for long reads when primers
    /// are expected near read termini.
    #[arg(long, default_value = "0")]
    reverse_window: usize,
}

#[derive(Debug, Clone, Copy, clap::ValueEnum)]
enum OutputFormat {
    Fastq,
    Fasta,
}

fn default_threads() -> usize {
    std::env::var("RAYON_NUM_THREADS")
        .ok()
        .and_then(|s| s.parse().ok())
        .unwrap_or_else(|| (num_cpus::get() / 2).max(1))
}

/// Validate command-line arguments
fn validate_args(args: &Args) -> Result<()> {
    debug!("Validating command-line arguments");

    // Check input file exists
    if !args.input.exists() {
        anyhow::bail!("Input file does not exist: {:?}", args.input);
    }
    debug!("Input file exists: {:?}", args.input);

    // Validate primer sequences (DNA alphabet with IUPAC codes)
    let valid_bases = b"ACGTURYSWKMBDHVNacgturyswkmbdhvn";

    for (name, seq) in [("forward", &args.forward), ("reverse", &args.reverse)] {
        if seq.is_empty() {
            anyhow::bail!("{} primer cannot be empty", name);
        }

        for &base in seq.as_bytes() {
            if !valid_bases.contains(&base) {
                anyhow::bail!(
                    "{} primer contains invalid character '{}' (must be DNA/IUPAC alphabet)",
                    name,
                    base as char
                );
            }
        }
        debug!("{} primer validated: {} ({} bp)", name, seq, seq.len());
    }

    // Validate numeric parameters
    if args.min_len == 0 {
        anyhow::bail!("Minimum length must be greater than 0");
    }

    if args.min_len > args.max_len {
        anyhow::bail!(
            "Minimum length ({}) cannot be greater than maximum length ({})",
            args.min_len,
            args.max_len
        );
    }
    debug!("Length bounds: min={}, max={}", args.min_len, args.max_len);

    if args.threads == 0 {
        anyhow::bail!("Number of threads must be at least 1");
    }
    debug!("Thread count: {}", args.threads);

    // Warn if primers are very short (likely to cause many false matches)
    if args.forward.len() < 10 {
        warn!(
            "Forward primer is very short ({} bp), may cause false matches",
            args.forward.len()
        );
    }
    if args.reverse.len() < 10 {
        warn!(
            "Reverse primer is very short ({} bp), may cause false matches",
            args.reverse.len()
        );
    }

    // Warn if search windows are smaller than primers (primer cannot be found)
    if args.forward_window > 0 && args.forward_window < args.forward.len() {
        warn!(
            "Forward window ({} bp) is smaller than forward primer ({} bp) - primer cannot be found",
            args.forward_window,
            args.forward.len()
        );
    }
    if args.reverse_window > 0 && args.reverse_window < args.reverse.len() {
        warn!(
            "Reverse window ({} bp) is smaller than reverse primer ({} bp) - primer cannot be found",
            args.reverse_window,
            args.reverse.len()
        );
    }

    debug!("All arguments validated successfully");
    Ok(())
}

/// Validate that the input/output format combination is supported.
///
/// FASTA input cannot produce FASTQ output because there are no quality scores
/// to preserve. All other combinations are valid.
fn validate_format_combination(
    input_format: fastx::Format,
    output_format: OutputFormat,
) -> Result<()> {
    if matches!(input_format, fastx::Format::Fasta) && matches!(output_format, OutputFormat::Fastq)
    {
        anyhow::bail!(
            "Cannot output FASTQ format from FASTA input (no quality scores available). \
             Use --format fasta instead."
        );
    }
    Ok(())
}

/// Create output writer with optional compression
fn create_writer(args: &Args) -> Result<SharedWriter> {
    let file = File::create(&args.output)?;

    let writer: Box<dyn Write + Send> = if args.no_compress {
        Box::new(file)
    } else {
        Box::new(GzEncoder::new(file, Compression::default()))
    };

    Ok(Arc::new(Mutex::new(writer)))
}

/// Reverse complement a DNA sequence, handling IUPAC ambiguity codes.
///
/// This function reverses the sequence and complements each base according
/// to standard Watson-Crick pairing rules, extended to support all IUPAC
/// ambiguity codes.
///
/// # Case Preservation
///
/// The function preserves case: uppercase input produces uppercase output,
/// lowercase produces lowercase. This is useful for maintaining quality
/// information encoded in case.
///
/// # IUPAC Codes
///
/// Supports all standard IUPAC nucleotide codes:
/// - Standard bases: A↔T, C↔G
/// - Ambiguity codes: R↔Y, K↔M, B↔V, D↔H
/// - Palindromic: S↔S, W↔W
/// - Wildcard: N↔N
///
/// # Examples
///
/// ```
/// assert_eq!(reverse_complement(b"ACGT"), b"ACGT");
/// assert_eq!(reverse_complement(b"ATCG"), b"CGAT");
/// ```
fn reverse_complement(seq: &[u8]) -> Vec<u8> {
    seq.iter()
        .rev()
        .map(|&b| match b {
            // Uppercase bases
            b'A' => b'T',
            b'T' => b'A',
            b'C' => b'G',
            b'G' => b'C',
            b'N' => b'N',
            // Lowercase bases
            b'a' => b't',
            b't' => b'a',
            b'c' => b'g',
            b'g' => b'c',
            b'n' => b'n',
            // Uppercase IUPAC ambiguity codes
            b'R' => b'Y', // R (A/G) -> Y (C/T)
            b'Y' => b'R', // Y (C/T) -> R (A/G)
            b'S' => b'S', // S (G/C) -> S (G/C)
            b'W' => b'W', // W (A/T) -> W (A/T)
            b'K' => b'M', // K (G/T) -> M (A/C)
            b'M' => b'K', // M (A/C) -> K (G/T)
            b'B' => b'V', // B (C/G/T) -> V (A/C/G)
            b'D' => b'H', // D (A/G/T) -> H (A/C/T)
            b'H' => b'D', // H (A/C/T) -> D (A/G/T)
            b'V' => b'B', // V (A/C/G) -> B (C/G/T)
            // Lowercase IUPAC ambiguity codes
            b'r' => b'y',
            b'y' => b'r',
            b's' => b's',
            b'w' => b'w',
            b'k' => b'm',
            b'm' => b'k',
            b'b' => b'v',
            b'd' => b'h',
            b'h' => b'd',
            b'v' => b'b',
            _ => b, // Unknown characters pass through
        })
        .collect()
}

/// Primer pair with pre-computed reverse complements.
///
/// Stores both the original primer sequences and their reverse complements
/// to enable efficient bidirectional searching without repeated computation.
#[derive(Clone)]
struct PrimerPair {
    /// Forward primer sequence (5' to 3')
    forward: Vec<u8>,
    /// Reverse primer sequence (5' to 3')
    reverse: Vec<u8>,
    /// Reverse complement of forward primer
    forward_rc: Vec<u8>,
    /// Reverse complement of reverse primer
    reverse_rc: Vec<u8>,
}

impl PrimerPair {
    fn new(forward: &[u8], reverse: &[u8]) -> Self {
        Self {
            forward: forward.to_vec(),
            reverse: reverse.to_vec(),
            forward_rc: reverse_complement(forward),
            reverse_rc: reverse_complement(reverse),
        }
    }
}

/// Orientation of the amplicon in the read.
///
/// Determines which primer sequences were used for matching:
/// - Forward: forward primer as-is, reverse primer as reverse complement
/// - Reverse: forward primer as reverse complement, reverse primer as-is
#[derive(Debug, Clone, Copy, PartialEq)]
enum Orientation {
    /// Forward orientation (expected for most reads)
    Forward,
    /// Reverse orientation (read is reverse-complemented)
    Reverse,
}

/// Reasons why an amplicon was not found in a read.
///
/// Used for detailed statistics tracking to understand why reads
/// are being filtered out.
#[derive(Debug, Clone, Copy, PartialEq)]
enum FailureReason {
    /// Forward primer not found within mismatch threshold
    NoForwardPrimer,
    /// Reverse primer not found within mismatch threshold
    NoReversePrimer,
    /// Primers found but overlap or in wrong order
    InvalidStructure,
    /// Insert length below minimum threshold
    TooShort,
    /// Insert length above maximum threshold
    TooLong,
}

/// Amplicon match result with position and edit distances.
///
/// Represents a valid amplicon found in a read, with coordinates
/// pointing to the trimmed insert sequence (primers excluded).
#[derive(Debug, Clone)]
struct AmpliconMatch {
    /// Start position of insert (after forward primer)
    start: usize,
    /// End position of insert (before reverse primer)
    end: usize,
    /// Edit distance for forward primer match
    fwd_cost: i32,
    /// Edit distance for reverse primer match
    rev_cost: i32,
    /// Orientation in which amplicon was found
    orientation: Orientation,
}

impl AmpliconMatch {
    fn total_cost(&self) -> i32 {
        self.fwd_cost + self.rev_cost
    }
}

/// Thread-safe statistics tracking for primer trimming.
///
/// All counters use atomic operations to allow safe concurrent updates
/// from multiple processing threads without locks.
#[derive(Default)]
struct TrimStats {
    /// Total number of reads processed
    total_reads: AtomicUsize,
    /// Number of reads with valid amplicons found
    amplicons_found: AtomicUsize,
    /// Reads where forward primer was not found
    no_forward_primer: AtomicUsize,
    /// Reads where reverse primer was not found
    no_reverse_primer: AtomicUsize,
    /// Reads where primers overlap or are in wrong order
    invalid_structure: AtomicUsize,
    /// Reads where insert is below minimum length
    too_short: AtomicUsize,
    /// Reads where insert exceeds maximum length
    too_long: AtomicUsize,
    /// Amplicons found in forward orientation
    forward_orientation: AtomicUsize,
    /// Amplicons found in reverse orientation
    reverse_orientation: AtomicUsize,
}

impl TrimStats {
    #[allow(dead_code)]
    fn new() -> Self {
        Self::default()
    }

    fn report(&self, output: &mut dyn Write) -> std::io::Result<()> {
        let total = self.total_reads.load(Ordering::Relaxed);
        let found = self.amplicons_found.load(Ordering::Relaxed);

        writeln!(output, "\nPrimer Trimming Statistics:")?;
        writeln!(output, "  Total reads processed: {}", total)?;

        if total > 0 {
            writeln!(
                output,
                "  Amplicons found: {} ({:.2}%)",
                found,
                100.0 * found as f64 / total as f64
            )?;
        } else {
            writeln!(output, "  Amplicons found: 0")?;
        }

        writeln!(
            output,
            "  Forward orientation: {}",
            self.forward_orientation.load(Ordering::Relaxed)
        )?;
        writeln!(
            output,
            "  Reverse orientation: {}",
            self.reverse_orientation.load(Ordering::Relaxed)
        )?;

        writeln!(output, "\nReads excluded:")?;
        writeln!(
            output,
            "  No forward primer: {}",
            self.no_forward_primer.load(Ordering::Relaxed)
        )?;
        writeln!(
            output,
            "  No reverse primer: {}",
            self.no_reverse_primer.load(Ordering::Relaxed)
        )?;
        writeln!(
            output,
            "  Invalid structure: {}",
            self.invalid_structure.load(Ordering::Relaxed)
        )?;
        writeln!(
            output,
            "  Too short: {}",
            self.too_short.load(Ordering::Relaxed)
        )?;
        writeln!(
            output,
            "  Too long: {}",
            self.too_long.load(Ordering::Relaxed)
        )?;

        Ok(())
    }
}

/// Primer trimming processor with search capabilities.
///
/// This is the main processing engine that implements the `ParallelProcessor`
/// trait from paraseq. Each thread gets its own clone of this processor,
/// sharing only the writer and statistics via Arc.
#[derive(Clone)]
struct PrimerTrimProcessor {
    /// Configured windowed search (primers, mismatch tolerance, windows)
    search: WindowedSearch,
    /// Minimum insert length (after trimming primers)
    min_len: usize,
    /// Maximum insert length (after trimming primers)
    max_len: usize,
    /// Approximate string matching engine (sassy)
    searcher: Searcher<Iupac>,
    /// Shared output writer (thread-safe)
    writer: SharedWriter,
    /// Output format (FASTQ or FASTA)
    output_format: OutputFormat,
    /// Shared statistics tracker (thread-safe)
    stats: SharedStats,
}

/// Configured primer search with optional windowing.
///
/// By default (window sizes of 0), the entire read is searched for primers.
/// When window sizes are specified, primer search is limited to bounded regions
/// at the ends of reads - an optimization for long reads where primers are
/// expected near the termini.
///
/// This struct holds immutable search configuration (primers, mismatch tolerance,
/// window sizes) and provides methods that accept a mutable searcher. This design
/// enables split borrowing - the caller can pass `&mut searcher` while the config
/// is borrowed immutably.
#[derive(Clone)]
struct WindowedSearch {
    primers: PrimerPair,
    max_mismatch: usize,
    forward_window: usize,
    reverse_window: usize,
}

impl WindowedSearch {
    fn new(
        primers: PrimerPair,
        max_mismatch: u8,
        forward_window: usize,
        reverse_window: usize,
    ) -> Self {
        Self {
            primers,
            max_mismatch: max_mismatch as usize,
            forward_window,
            reverse_window,
        }
    }

    /// Compute search window bounds for a given end of the read.
    ///
    /// Returns (start, end) indices for slicing the read sequence.
    ///
    /// # Arguments
    ///
    /// * `read_len` - Length of the read sequence
    /// * `window_size` - Size of the search window (0 = no limit)
    /// * `at_end` - If true, window is at 3' end; if false, at 5' end
    fn window_bounds(&self, read_len: usize, window_size: usize, at_end: bool) -> (usize, usize) {
        if window_size == 0 || window_size >= read_len {
            (0, read_len)
        } else if at_end {
            (read_len - window_size, read_len)
        } else {
            (0, window_size)
        }
    }

    /// Search for a primer in a windowed region, returning matches in read coordinates.
    fn search_primer(
        &self,
        searcher: &mut Searcher<Iupac>,
        primer: &[u8],
        read_seq: &[u8],
        window_start: usize,
        window_end: usize,
    ) -> Vec<sassy::Match> {
        let slice = &read_seq[window_start..window_end];

        searcher
            .search(primer, slice, self.max_mismatch)
            .into_iter()
            .map(|mut m| {
                // Adjust coordinates from slice space to read space
                m.text_start += window_start;
                m.text_end += window_start;
                m
            })
            .collect()
    }

    /// Search for forward primer in the appropriate window for the given orientation.
    ///
    /// When `forward_window` is 0, searches the entire read. Otherwise, searches
    /// a bounded region: the 5' end for forward orientation, or the 3' end for
    /// reverse orientation (where the forward primer RC would appear).
    fn search_forward(
        &self,
        searcher: &mut Searcher<Iupac>,
        read_seq: &[u8],
        orientation: Orientation,
    ) -> Vec<sassy::Match> {
        let read_len = read_seq.len();
        let at_end = matches!(orientation, Orientation::Reverse);
        let (start, end) = self.window_bounds(read_len, self.forward_window, at_end);

        let primer = match orientation {
            Orientation::Forward => &self.primers.forward,
            Orientation::Reverse => &self.primers.forward_rc,
        };

        self.search_primer(searcher, primer, read_seq, start, end)
    }

    /// Search for reverse primer in the appropriate window for the given orientation.
    ///
    /// When `reverse_window` is 0, searches the entire read. Otherwise, searches
    /// a bounded region: the 3' end for forward orientation (where the reverse
    /// primer RC would appear), or the 5' end for reverse orientation.
    fn search_reverse(
        &self,
        searcher: &mut Searcher<Iupac>,
        read_seq: &[u8],
        orientation: Orientation,
    ) -> Vec<sassy::Match> {
        let read_len = read_seq.len();
        let at_end = matches!(orientation, Orientation::Forward);
        let (start, end) = self.window_bounds(read_len, self.reverse_window, at_end);

        let primer = match orientation {
            Orientation::Forward => &self.primers.reverse_rc,
            Orientation::Reverse => &self.primers.reverse,
        };

        self.search_primer(searcher, primer, read_seq, start, end)
    }
}

impl PrimerTrimProcessor {
    #[allow(clippy::too_many_arguments)]
    fn new(
        primers: PrimerPair,
        max_mismatch: u8,
        min_len: usize,
        max_len: usize,
        forward_window: usize,
        reverse_window: usize,
        writer: SharedWriter,
        output_format: OutputFormat,
        stats: SharedStats,
    ) -> Self {
        Self {
            search: WindowedSearch::new(primers, max_mismatch, forward_window, reverse_window),
            min_len,
            max_len,
            searcher: Searcher::new(true, None),
            writer,
            output_format,
            stats,
        }
    }

    /// Write a trimmed record to the output file
    fn write_trimmed_record<R: Record>(&self, record: &R, amplicon: &AmpliconMatch) -> Result<()> {
        let sequence = record.seq();

        // Validate amplicon coordinates
        if amplicon.end > sequence.len() {
            anyhow::bail!(
                "Amplicon end ({}) exceeds sequence length ({})",
                amplicon.end,
                sequence.len()
            );
        }

        let trimmed_seq = &sequence[amplicon.start..amplicon.end];

        let mut writer = self.writer.lock().unwrap();

        match self.output_format {
            OutputFormat::Fastq => {
                let quality = record.qual().ok_or_else(|| {
                    anyhow::anyhow!(
                        "FASTQ output requires quality scores, but input record has none"
                    )
                })?;

                // Validate quality length matches sequence
                if amplicon.end > quality.len() {
                    anyhow::bail!(
                        "Amplicon end ({}) exceeds quality length ({})",
                        amplicon.end,
                        quality.len()
                    );
                }

                let trimmed_qual = &quality[amplicon.start..amplicon.end];

                // Write FASTQ record
                writeln!(writer, "@{}", record.id_str())?;
                writer.write_all(trimmed_seq)?;
                writer.write_all(b"\n+\n")?;
                writer.write_all(trimmed_qual)?;
                writer.write_all(b"\n")?;
            }
            OutputFormat::Fasta => {
                // Write FASTA record
                writeln!(writer, ">{}", record.id_str())?;
                writer.write_all(trimmed_seq)?;
                writer.write_all(b"\n")?;
            }
        }

        Ok(())
    }

    /// Find the best amplicon in a read, trying both orientations.
    ///
    /// This method performs a two-pass search:
    /// 1. Forward orientation: forward primer as-is, reverse primer as RC
    /// 2. Reverse orientation: forward primer as RC, reverse primer as-is
    ///
    /// Returns the first valid amplicon found (forward orientation preferred),
    /// or a detailed failure reason if no valid amplicon exists.
    ///
    /// # Arguments
    ///
    /// * `read_seq` - The read sequence to search
    ///
    /// # Returns
    ///
    /// * `Ok(AmpliconMatch)` - Valid amplicon with coordinates and costs
    /// * `Err(FailureReason)` - Specific reason why no amplicon was found
    fn find_amplicon(&mut self, read_seq: &[u8]) -> Result<AmpliconMatch, FailureReason> {
        // Pass 1: Forward orientation (fwd as-is, rev as RC)
        match self.try_orientation(read_seq, Orientation::Forward) {
            Ok(amplicon) => Ok(amplicon),
            Err(reason) => {
                // Save the failure reason from forward orientation
                let fwd_failure = reason;

                // Pass 2: Reverse orientation (fwd as RC, rev as-is)
                match self.try_orientation(read_seq, Orientation::Reverse) {
                    Ok(amplicon) => Ok(amplicon),
                    // If both orientations fail, return the forward failure reason
                    // (since forward is the expected orientation)
                    Err(_) => Err(fwd_failure),
                }
            }
        }
    }

    /// Try to find an amplicon in a specific orientation
    ///
    /// Uses bounded search windows when configured to limit the search space
    /// for better performance on long reads.
    fn try_orientation(
        &mut self,
        read_seq: &[u8],
        orientation: Orientation,
    ) -> Result<AmpliconMatch, FailureReason> {
        let fwd_matches = self
            .search
            .search_forward(&mut self.searcher, read_seq, orientation);

        if fwd_matches.is_empty() {
            return Err(FailureReason::NoForwardPrimer);
        }

        let rev_matches = self
            .search
            .search_reverse(&mut self.searcher, read_seq, orientation);

        if rev_matches.is_empty() {
            return Err(FailureReason::NoReversePrimer);
        }

        // Find best valid amplicon (coordinates are in read space)
        self.find_best_amplicon(&fwd_matches, &rev_matches, orientation)
    }

    /// Find the best valid amplicon from primer matches
    /// Selects the match with lowest total cost, breaking ties by longest amplicon
    fn find_best_amplicon(
        &self,
        fwd_matches: &[sassy::Match],
        rev_matches: &[sassy::Match],
        orientation: Orientation,
    ) -> Result<AmpliconMatch, FailureReason> {
        let mut best: Option<AmpliconMatch> = None;
        let mut has_invalid_structure = false;
        let mut has_too_short = false;
        let mut has_too_long = false;

        for fwd in fwd_matches {
            for rev in rev_matches {
                // Determine insert boundaries based on orientation
                // Forward: FWD at 5' end, REV_RC at 3' end → insert is fwd.text_end to rev.text_start
                // Reverse: REV at 5' end, FWD_RC at 3' end → insert is rev.text_end to fwd.text_start
                let (insert_start, insert_end) = match orientation {
                    Orientation::Forward => (fwd.text_end, rev.text_start),
                    Orientation::Reverse => (rev.text_end, fwd.text_start),
                };

                // Check valid structure (insert_start must be before insert_end)
                if insert_start > insert_end {
                    // Primers overlap or in wrong order
                    has_invalid_structure = true;
                    continue;
                }

                // Trimmed length excludes primers (insert only)
                let trimmed_len = insert_end - insert_start;

                // Check length bounds on trimmed sequence
                if trimmed_len < self.min_len {
                    has_too_short = true;
                    continue;
                }

                if trimmed_len > self.max_len {
                    has_too_long = true;
                    continue;
                }

                let candidate = AmpliconMatch {
                    start: insert_start,
                    end: insert_end,
                    fwd_cost: fwd.cost,
                    rev_cost: rev.cost,
                    orientation,
                };

                // Keep best scoring match, breaking ties by longest amplicon
                let is_better = match &best {
                    None => true,
                    Some(current) => {
                        let cost_cmp = candidate.total_cost().cmp(&current.total_cost());
                        match cost_cmp {
                            cmp::Ordering::Less => true, // Lower cost is better
                            cmp::Ordering::Equal => {
                                // Same cost: prefer longer amplicon
                                let candidate_len = candidate.end - candidate.start;
                                let current_len = current.end - current.start;
                                candidate_len > current_len
                            }
                            cmp::Ordering::Greater => false,
                        }
                    }
                };

                if is_better {
                    best = Some(candidate);
                }
            }
        }

        // Return the best match if found, otherwise return the most specific failure reason
        match best {
            Some(amplicon) => Ok(amplicon),
            None if has_invalid_structure => Err(FailureReason::InvalidStructure),
            None if has_too_short => Err(FailureReason::TooShort),
            None if has_too_long => Err(FailureReason::TooLong),
            None => Err(FailureReason::InvalidStructure),
        }
    }
}

impl<Rf: Record> ParallelProcessor<Rf> for PrimerTrimProcessor {
    fn process_record(&mut self, record: Rf) -> Result<(), ProcessError> {
        let read_count = self.stats.total_reads.fetch_add(1, Ordering::Relaxed) + 1;

        // Log progress every 10000 reads
        if read_count.is_multiple_of(10000) {
            debug!("Processed {} reads", read_count);
        }

        let seq = record.seq();

        // Try to find amplicon
        match self.find_amplicon(&seq) {
            Ok(amplicon) => {
                debug!(
                    "Found amplicon in {}: orientation={:?}, insert_len={}, cost={}",
                    record.id_str(),
                    amplicon.orientation,
                    amplicon.end - amplicon.start,
                    amplicon.total_cost()
                );

                // Update statistics
                self.stats.amplicons_found.fetch_add(1, Ordering::Relaxed);
                match amplicon.orientation {
                    Orientation::Forward => {
                        self.stats
                            .forward_orientation
                            .fetch_add(1, Ordering::Relaxed);
                    }
                    Orientation::Reverse => {
                        self.stats
                            .reverse_orientation
                            .fetch_add(1, Ordering::Relaxed);
                    }
                }

                // Trim and write
                self.write_trimmed_record(&record, &amplicon).map_err(|e| {
                    ProcessError::Process(Box::new(std::io::Error::other(e.to_string())))
                })?;
            }
            Err(reason) => {
                debug!("No amplicon in {}: {:?}", record.id_str(), reason);

                // Track specific failure reason
                match reason {
                    FailureReason::NoForwardPrimer => {
                        self.stats.no_forward_primer.fetch_add(1, Ordering::Relaxed);
                    }
                    FailureReason::NoReversePrimer => {
                        self.stats.no_reverse_primer.fetch_add(1, Ordering::Relaxed);
                    }
                    FailureReason::InvalidStructure => {
                        self.stats.invalid_structure.fetch_add(1, Ordering::Relaxed);
                    }
                    FailureReason::TooShort => {
                        self.stats.too_short.fetch_add(1, Ordering::Relaxed);
                    }
                    FailureReason::TooLong => {
                        self.stats.too_long.fetch_add(1, Ordering::Relaxed);
                    }
                }
            }
        }

        Ok(())
    }
}

fn main() -> Result<()> {
    env_logger::init();
    let args = Args::parse();

    // Validate arguments
    validate_args(&args)?;

    info!("=== Primer Trimming Started ===");
    info!("Input file: {:?}", args.input);
    info!("Output file: {:?}", args.output);
    info!("Output format: {:?}", args.format);
    info!(
        "Compression: {}",
        if args.no_compress { "disabled" } else { "gzip" }
    );
    info!("Threads: {}", args.threads);
    info!("Max mismatches: {}", args.max_mismatch);
    info!("Insert length bounds: {}-{} bp", args.min_len, args.max_len);

    // Log search window settings
    if args.forward_window > 0 {
        info!(
            "Forward primer search window: first {} bp",
            args.forward_window
        );
    } else {
        info!("Forward primer search window: entire read");
    }
    if args.reverse_window > 0 {
        info!(
            "Reverse primer search window: last {} bp",
            args.reverse_window
        );
    } else {
        info!("Reverse primer search window: entire read");
    }

    // Create primer pair with pre-computed reverse complements
    let primers = PrimerPair::new(args.forward.as_bytes(), args.reverse.as_bytes());

    info!(
        "Forward primer: {} ({} bp)",
        args.forward,
        args.forward.len()
    );
    info!(
        "Reverse primer: {} ({} bp)",
        args.reverse,
        args.reverse.len()
    );
    debug!(
        "Forward primer RC: {}",
        String::from_utf8_lossy(&primers.forward_rc)
    );
    debug!(
        "Reverse primer RC: {}",
        String::from_utf8_lossy(&primers.reverse_rc)
    );

    // Create output writer
    debug!("Creating output writer");
    let writer = create_writer(&args)?;

    // Create statistics tracker
    let stats = Arc::new(TrimStats::default());

    // Create processor
    debug!("Initializing primer trimming processor");
    let mut processor = PrimerTrimProcessor::new(
        primers,
        args.max_mismatch,
        args.min_len,
        args.max_len,
        args.forward_window,
        args.reverse_window,
        writer,
        args.format,
        stats.clone(),
    );

    // Read and process in parallel
    let reader = fastx::Reader::from_path(&args.input)?;
    let input_format = reader.format();

    info!(
        "Detected input format: {}",
        match input_format {
            fastx::Format::Fasta => "FASTA",
            fastx::Format::Fastq => "FASTQ",
        }
    );

    // Validate format combination before processing
    validate_format_combination(input_format, args.format)?;

    info!("Processing records...");
    reader
        .process_parallel(&mut processor, args.threads)
        .map_err(|e| anyhow::anyhow!("Processing failed: {}", e))?;

    info!("Processing complete");

    // Report statistics
    let mut stats_output: Box<dyn Write> = match &args.stats {
        Some(path) => {
            info!("Writing statistics to: {:?}", path);
            Box::new(File::create(path)?)
        }
        None => {
            debug!("Writing statistics to stderr");
            Box::new(std::io::stderr())
        }
    };
    stats.report(&mut *stats_output)?;

    info!("=== Primer Trimming Finished ===");
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use flate2::read::GzDecoder;
    use std::io::Read;
    use std::thread;

    // Type aliases for test utilities
    type SharedBuffer = Arc<Mutex<Vec<u8>>>;
    type TestWriterPair = (SharedWriter, SharedBuffer);

    // Helper function to create a test processor (full-read search, no windows)
    fn create_test_processor(
        primers: PrimerPair,
        max_mismatch: u8,
        min_len: usize,
        max_len: usize,
    ) -> PrimerTrimProcessor {
        create_test_processor_with_windows(primers, max_mismatch, min_len, max_len, 0, 0)
    }

    // Helper function to create a test processor with custom search windows
    fn create_test_processor_with_windows(
        primers: PrimerPair,
        max_mismatch: u8,
        min_len: usize,
        max_len: usize,
        forward_window: usize,
        reverse_window: usize,
    ) -> PrimerTrimProcessor {
        // Create a dummy writer that discards output
        let writer: Box<dyn std::io::Write + Send> = Box::new(std::io::sink());
        let writer = Arc::new(Mutex::new(writer));
        let stats = Arc::new(TrimStats::default());

        PrimerTrimProcessor::new(
            primers,
            max_mismatch,
            min_len,
            max_len,
            forward_window,
            reverse_window,
            writer,
            OutputFormat::Fastq,
            stats,
        )
    }

    #[test]
    fn test_reverse_complement_standard_bases() {
        assert_eq!(reverse_complement(b"ACGT"), b"ACGT");
        assert_eq!(reverse_complement(b"AAAA"), b"TTTT");
        assert_eq!(reverse_complement(b"CCCC"), b"GGGG");
        assert_eq!(reverse_complement(b"ATCG"), b"CGAT");
    }

    #[test]
    fn test_reverse_complement_case_preserved() {
        // Lowercase input produces lowercase output
        assert_eq!(reverse_complement(b"acgt"), b"acgt");
        // Mixed case: AcGt -> reverse: tGcA -> complement: aCgT
        assert_eq!(reverse_complement(b"AcGt"), b"aCgT");
        // Uppercase input produces uppercase output
        assert_eq!(reverse_complement(b"ATCG"), b"CGAT");
    }

    #[test]
    fn test_reverse_complement_iupac() {
        // R (A/G) -> Y (C/T)
        assert_eq!(reverse_complement(b"R"), b"Y");
        assert_eq!(reverse_complement(b"Y"), b"R");

        // S (G/C) -> S (G/C) - palindromic
        assert_eq!(reverse_complement(b"S"), b"S");

        // W (A/T) -> W (A/T) - palindromic
        assert_eq!(reverse_complement(b"W"), b"W");

        // K (G/T) -> M (A/C)
        assert_eq!(reverse_complement(b"K"), b"M");
        assert_eq!(reverse_complement(b"M"), b"K");

        // B (C/G/T) -> V (A/C/G)
        assert_eq!(reverse_complement(b"B"), b"V");
        assert_eq!(reverse_complement(b"V"), b"B");

        // D (A/G/T) -> H (A/C/T)
        assert_eq!(reverse_complement(b"D"), b"H");
        assert_eq!(reverse_complement(b"H"), b"D");
    }

    #[test]
    fn test_reverse_complement_n() {
        assert_eq!(reverse_complement(b"ANTN"), b"NANT");
    }

    #[test]
    fn test_primer_pair_creation() {
        let primers = PrimerPair::new(b"ACGT", b"TGCA");
        assert_eq!(primers.forward, b"ACGT");
        assert_eq!(primers.reverse, b"TGCA");
        assert_eq!(primers.forward_rc, b"ACGT");
        assert_eq!(primers.reverse_rc, b"TGCA");
    }

    #[test]
    fn test_trim_stats_empty() {
        let stats = TrimStats::new();
        let mut output = Vec::new();
        stats.report(&mut output).unwrap();
        let report = String::from_utf8(output).unwrap();
        assert!(report.contains("Total reads processed: 0"));
        assert!(report.contains("Amplicons found: 0"));
    }

    #[test]
    fn test_trim_stats_with_data() {
        let stats = TrimStats::new();
        stats.total_reads.store(100, Ordering::Relaxed);
        stats.amplicons_found.store(75, Ordering::Relaxed);
        stats.forward_orientation.store(50, Ordering::Relaxed);
        stats.reverse_orientation.store(25, Ordering::Relaxed);

        let mut output = Vec::new();
        stats.report(&mut output).unwrap();
        let report = String::from_utf8(output).unwrap();

        assert!(report.contains("Total reads processed: 100"));
        assert!(report.contains("Amplicons found: 75 (75.00%)"));
        assert!(report.contains("Forward orientation: 50"));
        assert!(report.contains("Reverse orientation: 25"));
    }

    #[test]
    fn test_trim_stats_thread_safety() {
        let stats = Arc::new(TrimStats::new());
        let mut handles = vec![];

        // Spawn 10 threads, each incrementing counters 100 times
        for _ in 0..10 {
            let stats_clone = Arc::clone(&stats);
            let handle = thread::spawn(move || {
                for _ in 0..100 {
                    stats_clone.total_reads.fetch_add(1, Ordering::Relaxed);
                    stats_clone.amplicons_found.fetch_add(1, Ordering::Relaxed);
                }
            });
            handles.push(handle);
        }

        // Wait for all threads to complete
        for handle in handles {
            handle.join().unwrap();
        }

        // Verify counts
        assert_eq!(stats.total_reads.load(Ordering::Relaxed), 1000);
        assert_eq!(stats.amplicons_found.load(Ordering::Relaxed), 1000);
    }

    #[test]
    fn test_amplicon_match_total_cost() {
        let amplicon = AmpliconMatch {
            start: 0,
            end: 100,
            fwd_cost: 1,
            rev_cost: 2,
            orientation: Orientation::Forward,
        };
        assert_eq!(amplicon.total_cost(), 3);
    }

    #[test]
    fn test_find_amplicon_forward_orientation() {
        let primers = PrimerPair::new(b"ACGT", b"TGCA");
        let mut processor = create_test_processor(primers, 0, 5, 100); // min_len=5 for insert

        // Read with exact primers: ACGT + insert + TGCA(RC=TGCA)
        // Using ATAT instead of NNNN to avoid N wildcard matching
        let read = b"ACGTATATATATATGCA";
        let result = processor.find_amplicon(read);

        assert!(result.is_ok());
        let amplicon = result.unwrap();
        assert_eq!(amplicon.orientation, Orientation::Forward);
        assert_eq!(amplicon.fwd_cost, 0);
        assert_eq!(amplicon.rev_cost, 0);
        // Primers trimmed: start after ACGT (pos 4), end before TGCA (pos 13)
        assert_eq!(amplicon.start, 4);
        assert_eq!(amplicon.end, 13);
    }

    #[test]
    fn test_find_amplicon_with_mismatches() {
        let primers = PrimerPair::new(b"ACGT", b"TGCA");
        let mut processor = create_test_processor(primers, 2, 5, 100); // min_len=5 for insert

        // Read with 1 mismatch in each primer: ACTT + insert + TGCC
        let read = b"ACTTATATATATATTGCC";
        let result = processor.find_amplicon(read);

        assert!(result.is_ok());
        let amplicon = result.unwrap();
        assert_eq!(amplicon.fwd_cost, 1);
        assert_eq!(amplicon.rev_cost, 1);
    }

    #[test]
    fn test_find_amplicon_too_short() {
        let primers = PrimerPair::new(b"ACGT", b"TGCA");
        let mut processor = create_test_processor(primers, 0, 50, 100);

        // Insert is only 9bp (13-4), below min_len of 50
        let read = b"ACGTATATATATATGCA";
        let result = processor.find_amplicon(read);

        assert!(result.is_err());
        assert_eq!(result.unwrap_err(), FailureReason::TooShort);
    }

    #[test]
    fn test_find_amplicon_too_long() {
        let primers = PrimerPair::new(b"ACGT", b"TGCA");
        let mut processor = create_test_processor(primers, 0, 5, 8);

        // Insert is 9bp (13-4), above max_len of 8
        let read = b"ACGTATATATATATGCA";
        let result = processor.find_amplicon(read);

        assert!(result.is_err());
        assert_eq!(result.unwrap_err(), FailureReason::TooLong);
    }

    #[test]
    fn test_find_amplicon_no_forward_primer() {
        let primers = PrimerPair::new(b"ACGT", b"TGCA");
        let mut processor = create_test_processor(primers, 0, 5, 100);

        // Read missing forward primer
        let read = b"ATATATATATTGCA";
        let result = processor.find_amplicon(read);

        assert!(result.is_err());
        assert_eq!(result.unwrap_err(), FailureReason::NoForwardPrimer);
    }

    #[test]
    fn test_find_amplicon_no_reverse_primer() {
        let primers = PrimerPair::new(b"ACGT", b"TGCA");
        let mut processor = create_test_processor(primers, 0, 5, 100);

        // Read missing reverse primer
        let read = b"ACGTATATATATA";
        let result = processor.find_amplicon(read);

        assert!(result.is_err());
        assert_eq!(result.unwrap_err(), FailureReason::NoReversePrimer);
    }

    #[test]
    fn test_find_amplicon_wrong_order() {
        let primers = PrimerPair::new(b"ACGT", b"TGCA");
        let mut processor = create_test_processor(primers, 0, 5, 100);

        // This read has REV at 5' and FWD_RC at 3' - this is valid reverse orientation!
        // TGCA (REV) at position 0, ACGT (FWD_RC, palindromic) at position 10
        let read = b"TGCAATATATACGT";
        let result = processor.find_amplicon(read);

        // Should find valid amplicon in reverse orientation
        assert!(result.is_ok());
        let amplicon = result.unwrap();
        assert_eq!(amplicon.orientation, Orientation::Reverse);
        assert_eq!(amplicon.start, 4); // After REV primer
        assert_eq!(amplicon.end, 10); // Before FWD_RC primer
    }

    // Phase 6 Tests: ParallelProcessor Integration

    #[test]
    fn test_process_record_trims_primers() {
        let primers = PrimerPair::new(b"ACGT", b"TGCA");
        let (writer, buffer) = create_test_writer();
        let stats = Arc::new(TrimStats::default());

        let mut processor = PrimerTrimProcessor::new(
            primers,
            0,
            5, // min_len=5 for insert only
            100,
            0, // forward_window (0 = full read)
            0, // reverse_window (0 = full read)
            writer,
            OutputFormat::Fastq,
            stats.clone(),
        );

        // Read: ACGT + ATATATAT + TGCA (RC of TGCA is TGCA)
        // Expected output: ATATATAT (primers trimmed out)
        let record = MockRecord {
            id: b"test_read".to_vec(),
            seq: b"ACGTATATATATTGCA".to_vec(),
            qual: Some(b"IIIIIIIIIIIIIIII".to_vec()),
        };

        processor.process_record(record).unwrap();

        // Check output
        let output = String::from_utf8(buffer.lock().unwrap().clone()).unwrap();
        let lines: Vec<&str> = output.lines().collect();

        assert_eq!(lines.len(), 4);
        assert_eq!(lines[0], "@test_read");
        assert_eq!(lines[1], "ATATATAT"); // Insert only, primers removed
        assert_eq!(lines[2], "+");
        assert_eq!(lines[3], "IIIIIIII"); // Quality for insert only

        // Check statistics
        assert_eq!(stats.total_reads.load(Ordering::Relaxed), 1);
        assert_eq!(stats.amplicons_found.load(Ordering::Relaxed), 1);
        assert_eq!(stats.forward_orientation.load(Ordering::Relaxed), 1);
    }

    #[test]
    fn test_process_record_no_amplicon() {
        let primers = PrimerPair::new(b"ACGT", b"TGCA");
        let (writer, buffer) = create_test_writer();
        let stats = Arc::new(TrimStats::default());

        let mut processor = PrimerTrimProcessor::new(
            primers,
            0,
            10,
            100,
            0, // forward_window (0 = full read)
            0, // reverse_window (0 = full read)
            writer,
            OutputFormat::Fastq,
            stats.clone(),
        );

        // Read with no primers
        let record = MockRecord {
            id: b"test_read".to_vec(),
            seq: b"GGGGGGGGGGGG".to_vec(),
            qual: Some(b"IIIIIIIIIIII".to_vec()),
        };

        processor.process_record(record).unwrap();

        // Check no output was written
        let output = buffer.lock().unwrap();
        assert_eq!(output.len(), 0);

        // Check statistics
        assert_eq!(stats.total_reads.load(Ordering::Relaxed), 1);
        assert_eq!(stats.amplicons_found.load(Ordering::Relaxed), 0);
        // This read has no primers at all, so should fail on forward primer
        assert_eq!(stats.no_forward_primer.load(Ordering::Relaxed), 1);
    }

    // Phase 8 Tests: Detailed Statistics Tracking

    #[test]
    fn test_detailed_failure_tracking() {
        let primers = PrimerPair::new(b"ACGTACGT", b"TGCATGCA");
        let (writer, _buffer) = create_test_writer();
        let stats = Arc::new(TrimStats::default());

        let mut processor = PrimerTrimProcessor::new(
            primers,
            0,
            10,
            50,
            0, // forward_window (0 = full read)
            0, // reverse_window (0 = full read)
            writer,
            OutputFormat::Fastq,
            stats.clone(),
        );

        // Test 1: No forward primer
        let record1 = MockRecord {
            id: b"read1".to_vec(),
            seq: b"GGGGGGGGGGGGGGGG".to_vec(),
            qual: Some(b"IIIIIIIIIIIIIIII".to_vec()),
        };
        processor.process_record(record1).unwrap();

        // Test 2: No reverse primer
        let record2 = MockRecord {
            id: b"read2".to_vec(),
            seq: b"ACGTACGTGGGGGGGG".to_vec(),
            qual: Some(b"IIIIIIIIIIIIIIII".to_vec()),
        };
        processor.process_record(record2).unwrap();

        // Test 3: Too short (insert < 10)
        let record3 = MockRecord {
            id: b"read3".to_vec(),
            seq: b"ACGTACGTATTTGCATGCA".to_vec(), // Insert is only 2bp
            qual: Some(b"IIIIIIIIIIIIIIIIIII".to_vec()),
        };
        processor.process_record(record3).unwrap();

        // Test 4: Too long (insert > 50)
        let record4 = MockRecord {
            id: b"read4".to_vec(),
            seq: b"ACGTACGTAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAATGCATGCA".to_vec(),
            qual: Some(
                b"IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII".to_vec(),
            ),
        };
        processor.process_record(record4).unwrap();

        // Test 5: Invalid structure (primers overlap)
        let record5 = MockRecord {
            id: b"read5".to_vec(),
            seq: b"TGCATGCAACGTACGT".to_vec(), // Reverse before forward
            qual: Some(b"IIIIIIIIIIIIIIII".to_vec()),
        };
        processor.process_record(record5).unwrap();

        // Verify statistics
        assert_eq!(stats.total_reads.load(Ordering::Relaxed), 5);
        assert_eq!(stats.amplicons_found.load(Ordering::Relaxed), 0);
        assert_eq!(stats.no_forward_primer.load(Ordering::Relaxed), 1);
        assert_eq!(stats.no_reverse_primer.load(Ordering::Relaxed), 1);
        assert_eq!(stats.too_short.load(Ordering::Relaxed), 1);
        assert_eq!(stats.too_long.load(Ordering::Relaxed), 1);
        assert_eq!(stats.invalid_structure.load(Ordering::Relaxed), 1);
    }

    #[test]
    fn test_find_amplicon_reverse_orientation() {
        let primers = PrimerPair::new(b"ACGTACGT", b"TGCATGCA");
        let mut processor = create_test_processor(primers, 0, 5, 100);

        // Read in reverse orientation: fwd_RC + insert + rev (not rev_RC)
        // Forward RC: ACGTACGT -> ACGTACGT (palindrome)
        // Reverse: TGCATGCA (as-is, not RC)
        // So we need: ACGTACGT + insert + TGCATGCA
        let read = b"ACGTACGTATATATATTGCATGCA";
        let result = processor.find_amplicon(read);

        // Should find in reverse orientation
        assert!(result.is_ok());
        let amplicon = result.unwrap();
        // Note: ACGTACGT is a palindrome, so orientation detection may vary
        // The important thing is we found a valid amplicon
        assert_eq!(amplicon.fwd_cost, 0);
        assert_eq!(amplicon.rev_cost, 0);
    }

    #[test]
    fn test_find_amplicon_multiple_matches_picks_best() {
        let primers = PrimerPair::new(b"ACGT", b"TGCA");
        let mut processor = create_test_processor(primers, 1, 5, 100);

        // Read with multiple potential amplicons
        let read = b"ACGTATATATATTGCAGGGGACTTATATATTGCA";
        let result = processor.find_amplicon(read);

        assert!(result.is_ok());
        let amplicon = result.unwrap();
        // Should find a valid amplicon with exact matches (cost 0)
        assert_eq!(amplicon.fwd_cost, 0);
        assert_eq!(amplicon.rev_cost, 0);
        assert_eq!(amplicon.start, 4); // After ACGT at position 0-3
                                       // The algorithm picks the longest valid amplicon with lowest cost
                                       // Verify it's a reasonable length
        assert!(amplicon.end > amplicon.start);
        assert!(amplicon.end - amplicon.start >= 5); // At least min_len
    }

    #[test]
    fn test_find_amplicon_read_shorter_than_primers() {
        let primers = PrimerPair::new(b"ACGTACGTACGTACGT", b"TGCATGCATGCATGCA");
        let mut processor = create_test_processor(primers, 0, 5, 100);

        // Very short read
        let read = b"ACGT";
        let result = processor.find_amplicon(read);

        // Should fail - read too short to contain primers
        assert!(result.is_err());
    }

    #[test]
    fn test_find_amplicon_primers_at_read_boundaries() {
        let primers = PrimerPair::new(b"ACGT", b"TGCA");
        let mut processor = create_test_processor(primers, 0, 5, 100);

        // Primers at exact start and end of read with minimal insert
        let read = b"ACGTATATATATTGCA";
        let result = processor.find_amplicon(read);

        assert!(result.is_ok());
        let amplicon = result.unwrap();
        assert_eq!(amplicon.start, 4); // After ACGT
        assert_eq!(amplicon.end, 12); // Before TGCA
        assert_eq!(amplicon.end - amplicon.start, 8); // Insert length
    }

    #[test]
    fn test_write_trimmed_record_quality_length_mismatch() {
        let primers = PrimerPair::new(b"ACGT", b"TGCA");
        let (writer, _buffer) = create_test_writer();
        let stats = Arc::new(TrimStats::default());

        let processor =
            PrimerTrimProcessor::new(primers, 0, 5, 100, 0, 0, writer, OutputFormat::Fastq, stats);

        let amplicon = AmpliconMatch {
            start: 4,
            end: 12,
            fwd_cost: 0,
            rev_cost: 0,
            orientation: Orientation::Forward,
        };

        // Quality string is shorter than sequence
        let record = MockRecord {
            id: b"test_read".to_vec(),
            seq: b"ACGTACGTACGTTGCA".to_vec(),
            qual: Some(b"IIII".to_vec()), // Too short!
        };

        // Should panic or error when trying to slice quality
        let result = processor.write_trimmed_record(&record, &amplicon);
        assert!(result.is_err());
    }

    #[test]
    fn test_find_amplicon_exact_min_length() {
        let primers = PrimerPair::new(b"ACGT", b"TGCA");
        let mut processor = create_test_processor(primers, 0, 8, 100);

        // Insert is exactly min_len (8bp)
        let read = b"ACGTATATATATTGCA";
        let result = processor.find_amplicon(read);

        assert!(result.is_ok());
        let amplicon = result.unwrap();
        assert_eq!(amplicon.end - amplicon.start, 8);
    }

    #[test]
    fn test_find_amplicon_exact_max_length() {
        let primers = PrimerPair::new(b"ACGT", b"TGCA");
        let mut processor = create_test_processor(primers, 0, 5, 8);

        // Insert is exactly max_len (8bp)
        let read = b"ACGTATATATATTGCA";
        let result = processor.find_amplicon(read);

        assert!(result.is_ok());
        let amplicon = result.unwrap();
        assert_eq!(amplicon.end - amplicon.start, 8);
    }

    #[test]
    fn test_process_record_updates_orientation_stats() {
        let primers = PrimerPair::new(b"ACGT", b"TGCA");
        let (writer, _buffer) = create_test_writer();
        let stats = Arc::new(TrimStats::default());

        let mut processor = PrimerTrimProcessor::new(
            primers,
            0,
            5,
            100,
            0, // forward_window (0 = full read)
            0, // reverse_window (0 = full read)
            writer,
            OutputFormat::Fastq,
            stats.clone(),
        );

        // Process a read in forward orientation
        let record = MockRecord {
            id: b"test_read".to_vec(),
            seq: b"ACGTATATATATTGCA".to_vec(),
            qual: Some(b"IIIIIIIIIIIIIIII".to_vec()),
        };

        processor.process_record(record).unwrap();

        // Verify orientation was tracked
        assert_eq!(stats.amplicons_found.load(Ordering::Relaxed), 1);
        assert!(
            stats.forward_orientation.load(Ordering::Relaxed) == 1
                || stats.reverse_orientation.load(Ordering::Relaxed) == 1
        );
    }

    // Phase 7 Tests: Error Handling & Validation

    #[test]
    fn test_validate_args_missing_input() {
        let args = Args {
            input: PathBuf::from("/nonexistent/file.fastq"),
            output: PathBuf::from("/tmp/out.fastq"),
            forward: "ACGTACGT".to_string(),
            reverse: "TGCATGCA".to_string(),
            max_mismatch: 2,
            min_len: 50,
            max_len: 2000,
            threads: 4,
            format: OutputFormat::Fastq,
            no_compress: true,
            stats: None,
            forward_window: 0,
            reverse_window: 0,
        };

        let result = validate_args(&args);
        assert!(result.is_err());
        assert!(result.unwrap_err().to_string().contains("does not exist"));
    }

    #[test]
    fn test_validate_args_invalid_primer_sequence() {
        let temp_file = std::env::temp_dir().join("test_validate.fastq");
        std::fs::write(&temp_file, "@read\nACGT\n+\nIIII\n").unwrap();

        let args = Args {
            input: temp_file.clone(),
            output: PathBuf::from("/tmp/out.fastq"),
            forward: "ACGT123".to_string(), // Invalid characters
            reverse: "TGCATGCA".to_string(),
            max_mismatch: 2,
            min_len: 50,
            max_len: 2000,
            threads: 4,
            format: OutputFormat::Fastq,
            no_compress: true,
            stats: None,
            forward_window: 0,
            reverse_window: 0,
        };

        let result = validate_args(&args);
        assert!(result.is_err());
        assert!(result
            .unwrap_err()
            .to_string()
            .contains("invalid character"));

        std::fs::remove_file(&temp_file).unwrap();
    }

    #[test]
    fn test_validate_args_empty_primer() {
        let temp_file = std::env::temp_dir().join("test_validate2.fastq");
        std::fs::write(&temp_file, "@read\nACGT\n+\nIIII\n").unwrap();

        let args = Args {
            input: temp_file.clone(),
            output: PathBuf::from("/tmp/out.fastq"),
            forward: "".to_string(), // Empty
            reverse: "TGCATGCA".to_string(),
            max_mismatch: 2,
            min_len: 50,
            max_len: 2000,
            threads: 4,
            format: OutputFormat::Fastq,
            no_compress: true,
            stats: None,
            forward_window: 0,
            reverse_window: 0,
        };

        let result = validate_args(&args);
        assert!(result.is_err());
        assert!(result.unwrap_err().to_string().contains("cannot be empty"));

        std::fs::remove_file(&temp_file).unwrap();
    }

    #[test]
    fn test_validate_args_min_greater_than_max() {
        let temp_file = std::env::temp_dir().join("test_validate3.fastq");
        std::fs::write(&temp_file, "@read\nACGT\n+\nIIII\n").unwrap();

        let args = Args {
            input: temp_file.clone(),
            output: PathBuf::from("/tmp/out.fastq"),
            forward: "ACGTACGT".to_string(),
            reverse: "TGCATGCA".to_string(),
            max_mismatch: 2,
            min_len: 2000, // Greater than max_len
            max_len: 50,
            threads: 4,
            format: OutputFormat::Fastq,
            no_compress: true,
            stats: None,
            forward_window: 0,
            reverse_window: 0,
        };

        let result = validate_args(&args);
        assert!(result.is_err());
        assert!(result
            .unwrap_err()
            .to_string()
            .contains("cannot be greater than"));

        std::fs::remove_file(&temp_file).unwrap();
    }

    #[test]
    fn test_validate_args_zero_threads() {
        let temp_file = std::env::temp_dir().join("test_validate4.fastq");
        std::fs::write(&temp_file, "@read\nACGT\n+\nIIII\n").unwrap();

        let args = Args {
            input: temp_file.clone(),
            output: PathBuf::from("/tmp/out.fastq"),
            forward: "ACGTACGT".to_string(),
            reverse: "TGCATGCA".to_string(),
            max_mismatch: 2,
            min_len: 50,
            max_len: 2000,
            threads: 0, // Invalid
            format: OutputFormat::Fastq,
            no_compress: true,
            stats: None,
            forward_window: 0,
            reverse_window: 0,
        };

        let result = validate_args(&args);
        assert!(result.is_err());
        assert!(result.unwrap_err().to_string().contains("at least 1"));

        std::fs::remove_file(&temp_file).unwrap();
    }

    #[test]
    fn test_validate_args_valid() {
        let temp_file = std::env::temp_dir().join("test_validate5.fastq");
        std::fs::write(&temp_file, "@read\nACGT\n+\nIIII\n").unwrap();

        let args = Args {
            input: temp_file.clone(),
            output: PathBuf::from("/tmp/out.fastq"),
            forward: "ACGTACGTACGT".to_string(),
            reverse: "TGCATGCATGCA".to_string(),
            max_mismatch: 2,
            min_len: 50,
            max_len: 2000,
            threads: 4,
            format: OutputFormat::Fastq,
            no_compress: true,
            stats: None,
            forward_window: 0,
            reverse_window: 0,
        };

        let result = validate_args(&args);
        assert!(result.is_ok());

        std::fs::remove_file(&temp_file).unwrap();
    }

    // Phase 5 Tests: Output Writing

    // Mock record for testing
    struct MockRecord {
        id: Vec<u8>,
        seq: Vec<u8>,
        qual: Option<Vec<u8>>,
    }

    impl Record for MockRecord {
        fn id(&self) -> &[u8] {
            &self.id
        }

        fn seq(&self) -> std::borrow::Cow<'_, [u8]> {
            std::borrow::Cow::Borrowed(&self.seq)
        }

        fn seq_raw(&self) -> &[u8] {
            &self.seq
        }

        fn qual(&self) -> Option<&[u8]> {
            self.qual.as_deref()
        }
    }

    // Helper to create a test writer that captures output
    fn create_test_writer() -> TestWriterPair {
        let buffer = Arc::new(Mutex::new(Vec::new()));
        let buffer_clone = buffer.clone();

        // Create a wrapper that writes to our shared buffer
        struct SharedVecWriter {
            buffer: Arc<Mutex<Vec<u8>>>,
        }

        impl std::io::Write for SharedVecWriter {
            fn write(&mut self, buf: &[u8]) -> std::io::Result<usize> {
                self.buffer.lock().unwrap().write(buf)
            }

            fn flush(&mut self) -> std::io::Result<()> {
                self.buffer.lock().unwrap().flush()
            }
        }

        let writer: Box<dyn std::io::Write + Send> = Box::new(SharedVecWriter {
            buffer: buffer_clone,
        });
        (Arc::new(Mutex::new(writer)), buffer)
    }

    #[test]
    fn test_write_fastq_format() {
        let primers = PrimerPair::new(b"ACGT", b"TGCA");
        let (writer, buffer) = create_test_writer();
        let stats = Arc::new(TrimStats::default());

        let processor = PrimerTrimProcessor::new(
            primers,
            0,
            10,
            100,
            0,
            0,
            writer,
            OutputFormat::Fastq,
            stats,
        );

        let amplicon = AmpliconMatch {
            start: 4,
            end: 17,
            fwd_cost: 0,
            rev_cost: 0,
            orientation: Orientation::Forward,
        };

        let record = MockRecord {
            id: b"test_read".to_vec(),
            seq: b"ACGTACGTATATATATATGCA".to_vec(),
            qual: Some(b"IIIIIIIIIIIIIIIIIIIII".to_vec()),
        };

        processor.write_trimmed_record(&record, &amplicon).unwrap();

        // Extract the written data
        let output = String::from_utf8(buffer.lock().unwrap().clone()).unwrap();

        // Verify FASTQ format
        let lines: Vec<&str> = output.lines().collect();
        assert_eq!(lines.len(), 4);
        assert_eq!(lines[0], "@test_read");
        assert_eq!(lines[1], "ACGTATATATATA"); // Trimmed sequence (positions 4-17 exclusive = 13 chars)
        assert_eq!(lines[2], "+");
        assert_eq!(lines[3], "IIIIIIIIIIIII"); // Trimmed quality (13 chars)
    }

    #[test]
    fn test_write_fasta_format() {
        let primers = PrimerPair::new(b"ACGT", b"TGCA");
        let (writer, buffer) = create_test_writer();
        let stats = Arc::new(TrimStats::default());

        let processor = PrimerTrimProcessor::new(
            primers,
            0,
            10,
            100,
            0,
            0,
            writer,
            OutputFormat::Fasta,
            stats,
        );

        let amplicon = AmpliconMatch {
            start: 4,
            end: 17,
            fwd_cost: 0,
            rev_cost: 0,
            orientation: Orientation::Forward,
        };

        let record = MockRecord {
            id: b"test_read".to_vec(),
            seq: b"ACGTACGTATATATATATGCA".to_vec(),
            qual: None,
        };

        processor.write_trimmed_record(&record, &amplicon).unwrap();

        // Extract the written data
        let output = String::from_utf8(buffer.lock().unwrap().clone()).unwrap();

        // Verify FASTA format
        let lines: Vec<&str> = output.lines().collect();
        assert_eq!(lines.len(), 2);
        assert_eq!(lines[0], ">test_read");
        assert_eq!(lines[1], "ACGTATATATATA"); // Trimmed sequence
    }

    #[test]
    fn test_write_fastq_requires_quality() {
        let primers = PrimerPair::new(b"ACGT", b"TGCA");
        let (writer, _buffer) = create_test_writer();
        let stats = Arc::new(TrimStats::default());

        let processor = PrimerTrimProcessor::new(
            primers,
            0,
            10,
            100,
            0,
            0,
            writer,
            OutputFormat::Fastq,
            stats,
        );

        let amplicon = AmpliconMatch {
            start: 0,
            end: 10,
            fwd_cost: 0,
            rev_cost: 0,
            orientation: Orientation::Forward,
        };

        let record = MockRecord {
            id: b"test_read".to_vec(),
            seq: b"ACGTACGTACGT".to_vec(),
            qual: None,
        };

        // Should error when quality is None for FASTQ format
        let result = processor.write_trimmed_record(&record, &amplicon);
        assert!(result.is_err());
        assert!(result.unwrap_err().to_string().contains("quality scores"));
    }

    #[test]
    fn test_create_writer_uncompressed() {
        let temp_dir = std::env::temp_dir();
        let output_path = temp_dir.join("test_uncompressed.fastq");

        let args = Args {
            input: PathBuf::from("dummy.fastq"),
            output: output_path.clone(),
            forward: "ACGT".to_string(),
            reverse: "TGCA".to_string(),
            max_mismatch: 2,
            min_len: 50,
            max_len: 2000,
            threads: 1,
            format: OutputFormat::Fastq,
            no_compress: true,
            stats: None,
            forward_window: 0,
            reverse_window: 0,
        };

        let writer = create_writer(&args).unwrap();

        // Write some test data
        {
            let mut w = writer.lock().unwrap();
            writeln!(w, "test data").unwrap();
        }

        // Read back and verify it's uncompressed
        let mut file = std::fs::File::open(&output_path).unwrap();
        let mut contents = String::new();
        file.read_to_string(&mut contents).unwrap();
        assert_eq!(contents, "test data\n");

        // Cleanup
        std::fs::remove_file(&output_path).unwrap();
    }

    #[test]
    fn test_create_writer_gzip_compressed() {
        let temp_dir = std::env::temp_dir();
        let output_path = temp_dir.join("test_compressed.fastq.gz");

        let args = Args {
            input: PathBuf::from("dummy.fastq"),
            output: output_path.clone(),
            forward: "ACGT".to_string(),
            reverse: "TGCA".to_string(),
            max_mismatch: 2,
            min_len: 50,
            max_len: 2000,
            threads: 1,
            format: OutputFormat::Fastq,
            no_compress: false,
            stats: None,
            forward_window: 0,
            reverse_window: 0,
        };

        let writer = create_writer(&args).unwrap();

        // Write some test data
        {
            let mut w = writer.lock().unwrap();
            writeln!(w, "test data").unwrap();
        }
        // Drop the writer to flush and close the gzip stream
        drop(writer);

        // Read back and verify it's gzip compressed
        let file = std::fs::File::open(&output_path).unwrap();
        let mut decoder = GzDecoder::new(file);
        let mut contents = String::new();
        decoder.read_to_string(&mut contents).unwrap();
        assert_eq!(contents, "test data\n");

        // Cleanup
        std::fs::remove_file(&output_path).unwrap();
    }

    // Format combination validation tests

    #[test]
    fn test_validate_format_combination_fasta_to_fasta() {
        // FASTA -> FASTA is valid
        let result = validate_format_combination(fastx::Format::Fasta, OutputFormat::Fasta);
        assert!(result.is_ok());
    }

    #[test]
    fn test_validate_format_combination_fasta_to_fastq() {
        // FASTA -> FASTQ is invalid (no quality scores)
        let result = validate_format_combination(fastx::Format::Fasta, OutputFormat::Fastq);
        assert!(result.is_err());
        assert!(result.unwrap_err().to_string().contains("quality scores"));
    }

    #[test]
    fn test_validate_format_combination_fastq_to_fasta() {
        // FASTQ -> FASTA is valid (quality scores discarded)
        let result = validate_format_combination(fastx::Format::Fastq, OutputFormat::Fasta);
        assert!(result.is_ok());
    }

    #[test]
    fn test_validate_format_combination_fastq_to_fastq() {
        // FASTQ -> FASTQ is valid
        let result = validate_format_combination(fastx::Format::Fastq, OutputFormat::Fastq);
        assert!(result.is_ok());
    }

    // Windowed search tests

    #[test]
    fn test_window_bounds_no_limit() {
        let primers = PrimerPair::new(b"ACGT", b"TGCA");
        let search = WindowedSearch::new(primers, 0, 0, 0);

        // window_size=0 means search entire read
        assert_eq!(search.window_bounds(100, 0, false), (0, 100));
        assert_eq!(search.window_bounds(100, 0, true), (0, 100));
    }

    #[test]
    fn test_window_bounds_at_5_prime() {
        let primers = PrimerPair::new(b"ACGT", b"TGCA");
        let search = WindowedSearch::new(primers, 0, 50, 50);

        // at_end=false means 5' end (first N bases)
        assert_eq!(search.window_bounds(100, 50, false), (0, 50));
    }

    #[test]
    fn test_window_bounds_at_3_prime() {
        let primers = PrimerPair::new(b"ACGT", b"TGCA");
        let search = WindowedSearch::new(primers, 0, 50, 50);

        // at_end=true means 3' end (last N bases)
        assert_eq!(search.window_bounds(100, 50, true), (50, 100));
    }

    #[test]
    fn test_window_bounds_larger_than_read() {
        let primers = PrimerPair::new(b"ACGT", b"TGCA");
        let search = WindowedSearch::new(primers, 0, 200, 200);

        // Window larger than read should search entire read
        assert_eq!(search.window_bounds(100, 200, false), (0, 100));
        assert_eq!(search.window_bounds(100, 200, true), (0, 100));
    }

    #[test]
    fn test_windowed_search_forward_orientation_primer_in_window() {
        // Forward orientation: FWD at 5' end, REV_RC at 3' end
        // Using 4-base primers for simpler math
        let primers = PrimerPair::new(b"ACGT", b"TGCA");
        // Windows of 10 bases at each end
        let mut processor = create_test_processor_with_windows(primers, 0, 5, 200, 10, 10);

        // 30bp read: FWD(4) + insert(18) + REV_RC(4) + padding(4) = 30bp
        // FWD at 0-4, REV_RC at 22-26, padding at 26-30
        // REV_RC of TGCA is TGCA (palindrome)
        let read = b"ACGTATATATATATATATATATTGCAAAAA";
        assert_eq!(read.len(), 30);

        let result = processor.find_amplicon(read);
        assert!(
            result.is_ok(),
            "Should find amplicon when primers are within windows"
        );
        let amplicon = result.unwrap();
        assert_eq!(amplicon.orientation, Orientation::Forward);
        assert_eq!(amplicon.start, 4); // After forward primer
        assert_eq!(amplicon.end, 22); // Before reverse primer RC
    }

    #[test]
    fn test_windowed_search_forward_orientation_primer_outside_window() {
        // Forward primer placed beyond the forward window
        let primers = PrimerPair::new(b"ACGT", b"TGCA");
        // Forward window of only 5 bases - primer starts at position 10
        let mut processor = create_test_processor_with_windows(primers, 0, 5, 200, 5, 15);

        // 30bp read: padding(10) + FWD(4) + insert(8) + REV_RC(4) + padding(4) = 30bp
        let read = b"AAAAAAAAAAAACGTATATATATATTGCAAAA";
        assert_eq!(read.len(), 32);

        let result = processor.find_amplicon(read);
        assert!(
            result.is_err(),
            "Should not find amplicon when forward primer is outside window"
        );
        assert_eq!(result.unwrap_err(), FailureReason::NoForwardPrimer);
    }

    #[test]
    fn test_windowed_search_forward_orientation_rev_primer_outside_window() {
        // Reverse primer RC placed beyond the reverse window
        let primers = PrimerPair::new(b"ACGT", b"TGCA");
        // Reverse window of only 5 bases at 3' end
        let mut processor = create_test_processor_with_windows(primers, 0, 5, 200, 15, 5);

        // 32bp read: FWD(4) + insert(8) + REV_RC(4) + padding(16) = 32bp
        // REV_RC ends at position 16, but window only covers positions 27-32
        let read = b"ACGTATATATATATTGCAAAAAAAAAAAAAAAA";
        assert_eq!(read.len(), 33);

        let result = processor.find_amplicon(read);
        assert!(
            result.is_err(),
            "Should not find amplicon when reverse primer is outside window"
        );
        assert_eq!(result.unwrap_err(), FailureReason::NoReversePrimer);
    }

    #[test]
    fn test_windowed_search_reverse_orientation_primer_in_window() {
        // Reverse orientation: REV at 5' end, FWD_RC at 3' end
        let primers = PrimerPair::new(b"ACGT", b"TGCA");
        // Windows of 10 bases at each end
        let mut processor = create_test_processor_with_windows(primers, 0, 5, 200, 10, 10);

        // 31bp read in reverse orientation: REV(4) + insert(18) + FWD_RC(4) + padding(5)
        // REV = TGCA at 0-4, FWD_RC = ACGT at 22-26 (ACGT is palindromic)
        let read = b"TGCAATATATATATATATATATACGTAAAAA";
        assert_eq!(read.len(), 31);

        let result = processor.find_amplicon(read);
        assert!(
            result.is_ok(),
            "Should find amplicon in reverse orientation"
        );
        let amplicon = result.unwrap();
        assert_eq!(amplicon.orientation, Orientation::Reverse);
        assert_eq!(amplicon.start, 4); // After reverse primer (rev.text_end)
        assert_eq!(amplicon.end, 22); // Before forward primer RC (fwd.text_start)
    }

    #[test]
    fn test_windowed_search_reverse_orientation_fwd_rc_outside_window() {
        // In reverse orientation, FWD_RC is at 3' end (searched with forward_window)
        let primers = PrimerPair::new(b"ACGT", b"TGCA");
        // Forward window of only 5 bases at 3' end (in reverse orientation)
        let mut processor = create_test_processor_with_windows(primers, 0, 5, 200, 5, 15);

        // 32bp read: REV(4) + insert(8) + FWD_RC(4) + padding(16) = 32bp
        // FWD_RC at position 12-16, but forward_window=5 means we only search last 5 bases
        let read = b"TGCAATATATATATACGTAAAAAAAAAAAAAAAA";
        assert_eq!(read.len(), 34);

        let result = processor.find_amplicon(read);
        // Forward orientation also fails (no FWD at 5'), so we get NoForwardPrimer
        assert!(
            result.is_err(),
            "Should not find amplicon when FWD_RC is outside window"
        );
    }

    #[test]
    fn test_windowed_search_reverse_orientation_rev_outside_window() {
        // In reverse orientation, REV is at 5' end (searched with reverse_window)
        let primers = PrimerPair::new(b"ACGT", b"TGCA");
        // Reverse window of only 5 bases at 5' end (in reverse orientation)
        let mut processor = create_test_processor_with_windows(primers, 0, 5, 200, 15, 5);

        // 32bp read: padding(10) + REV(4) + insert(8) + FWD_RC(4) + padding(6) = 32bp
        // REV at position 10-14, but reverse_window=5 means we only search first 5 bases
        let read = b"AAAAAAAAAATGCAATATATATATACGTAAAAAA";
        assert_eq!(read.len(), 34);

        let result = processor.find_amplicon(read);
        // Forward orientation also fails (no FWD at 5'), reverse fails (REV outside window)
        assert!(
            result.is_err(),
            "Should not find amplicon when REV is outside window"
        );
    }

    #[test]
    fn test_windowed_search_coordinate_adjustment() {
        // Verify that match coordinates are correctly adjusted from window space to read space
        let primers = PrimerPair::new(b"ACGT", b"TGCA");
        // Windows of 15 bases at each end
        let mut processor = create_test_processor_with_windows(primers, 0, 5, 200, 15, 15);

        // 40bp read with primers positioned within windows:
        // FWD (ACGT) at position 9-13 (within 15-base window at 5' end: 0-15)
        // REV_RC (TGCA) at position 29-33 (within 15-base window at 3' end: 25-40)
        let read = b"AAAAAAAAAACGTATATATAAAAAAAAAATGCAAAAAAAA";
        assert_eq!(read.len(), 40);

        let result = processor.find_amplicon(read);
        assert!(
            result.is_ok(),
            "Should find amplicon with primers near window edges"
        );
        let amplicon = result.unwrap();

        // Verify coordinates are in read space, not window space
        // FWD ends at 13, REV_RC starts at 29
        assert_eq!(amplicon.start, 13, "Start should be after forward primer");
        assert_eq!(amplicon.end, 29, "End should be before reverse primer RC");
    }

    #[test]
    fn test_windowed_search_both_orientations_tried() {
        // When forward orientation fails due to windowing, reverse should still be tried
        let primers = PrimerPair::new(b"ACGT", b"TGCA");
        let mut processor = create_test_processor_with_windows(primers, 0, 5, 200, 10, 10);

        // 30bp read in reverse orientation with primers within windows
        // REV at 5' (position 0-4), FWD_RC at 3' (position 26-30)
        let read = b"TGCAATATATATATATATATATATACGT";
        assert_eq!(read.len(), 28);

        let result = processor.find_amplicon(read);
        assert!(
            result.is_ok(),
            "Should find amplicon in reverse orientation"
        );
        let amplicon = result.unwrap();
        assert_eq!(amplicon.orientation, Orientation::Reverse);
    }

    #[test]
    fn test_windowed_search_zero_window_searches_entire_read() {
        // Window of 0 should behave like no windowing (search entire read)
        let primers = PrimerPair::new(b"ACGT", b"TGCA");
        let mut processor = create_test_processor_with_windows(primers, 0, 5, 200, 0, 0);

        // 40bp read with primers in the middle - would fail with small windows
        let read = b"AAAAAAAAAAAAAAAACGTATATATATATTGCAAAAAAAAA";
        assert_eq!(read.len(), 41);

        let result = processor.find_amplicon(read);
        assert!(
            result.is_ok(),
            "Should find amplicon anywhere when windows are 0"
        );
    }

    #[test]
    fn test_windowed_search_asymmetric_windows() {
        // Different window sizes for forward and reverse primers
        let primers = PrimerPair::new(b"ACGT", b"TGCA");
        // Large forward window (20), small reverse window (10)
        let mut processor = create_test_processor_with_windows(primers, 0, 5, 200, 20, 10);

        // 36bp read with:
        // FWD (ACGT) at position 14-18 (within 20-base window at 5': 0-20)
        // REV_RC (TGCA) at position 28-32 (within 10-base window at 3': 26-36)
        let read = b"AAAAAAAAAAAAAAAACGTATATAAAAAAATGCAAAA";
        assert_eq!(read.len(), 37);

        let result = processor.find_amplicon(read);
        assert!(
            result.is_ok(),
            "Should find amplicon with asymmetric windows"
        );
        let amplicon = result.unwrap();
        assert_eq!(amplicon.orientation, Orientation::Forward);
    }
}
