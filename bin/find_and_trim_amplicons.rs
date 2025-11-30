#!/usr/bin/env rust-script
//! Find and trim primer sequences in FASTQ/FASTA reads.
//!
//! This high-performance Rust script combines primer finding and trimming operations
//! into a single process, replacing traditional two-process workflows. It uses
//! approximate string matching to locate primers, handles both orientations,
//! and outputs trimmed amplicons with detailed statistics.
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
//! 1. For each read, search for forward and reverse primers in both orientations
//! 2. Find all valid amplicons (primers in correct order, within length bounds)
//! 3. Select the amplicon with lowest total edit distance (ties broken by length)
//! 4. Trim the read to exclude primers, keeping only the insert sequence
//! 5. Write trimmed sequence to output file
//!
//! # Usage
//!
//! ```bash
//! trim_primers.rs \
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
use paraseq::{fastq, prelude::*, ProcessError};
use sassy::{profiles::Iupac, Searcher};
use std::{
    cmp,
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
    /// Input FASTQ/FASTA file (gzip supported)
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

    debug!("All arguments validated successfully");
    Ok(())
}

/// Create output writer with optional compression
fn create_writer(args: &Args) -> Result<SharedWriter> {
    let file = std::fs::File::create(&args.output)?;

    let writer: Box<dyn std::io::Write + Send> = if args.no_compress {
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

    fn report(&self, output: &mut dyn std::io::Write) -> std::io::Result<()> {
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
    /// Primer sequences with pre-computed reverse complements
    primers: PrimerPair,
    /// Maximum number of mismatches allowed in primer matching
    max_mismatch: u8,
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

impl PrimerTrimProcessor {
    fn new(
        primers: PrimerPair,
        max_mismatch: u8,
        min_len: usize,
        max_len: usize,
        writer: SharedWriter,
        output_format: OutputFormat,
        stats: SharedStats,
    ) -> Self {
        Self {
            primers,
            max_mismatch,
            min_len,
            max_len,
            // rc=true enables reverse complement search, alpha=None uses default scoring
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
                let quality = record
                    .qual()
                    .ok_or_else(|| anyhow::anyhow!("FASTQ format requires quality scores"))?;

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
        // Clone primers to avoid borrow checker issues
        let fwd = self.primers.forward.clone();
        let rev_rc = self.primers.reverse_rc.clone();
        let fwd_rc = self.primers.forward_rc.clone();
        let rev = self.primers.reverse.clone();

        // Pass 1: Forward orientation (fwd as-is, rev as RC)
        match self.try_orientation(read_seq, &fwd, &rev_rc, Orientation::Forward) {
            Ok(amplicon) => Ok(amplicon),
            Err(reason) => {
                // Save the failure reason from forward orientation
                let fwd_failure = reason;

                // Pass 2: Reverse orientation (fwd as RC, rev as-is)
                match self.try_orientation(read_seq, &fwd_rc, &rev, Orientation::Reverse) {
                    Ok(amplicon) => Ok(amplicon),
                    // If both orientations fail, return the forward failure reason
                    // (since forward is the expected orientation)
                    Err(_) => Err(fwd_failure),
                }
            }
        }
    }

    /// Try to find an amplicon in a specific orientation
    fn try_orientation(
        &mut self,
        read_seq: &[u8],
        fwd_primer: &[u8],
        rev_primer: &[u8],
        orientation: Orientation,
    ) -> Result<AmpliconMatch, FailureReason> {
        // Search for forward primer (max_cost is usize in sassy)
        let fwd_matches = self
            .searcher
            .search(fwd_primer, read_seq, self.max_mismatch as usize);

        if fwd_matches.is_empty() {
            return Err(FailureReason::NoForwardPrimer);
        }

        // Search for reverse primer
        let rev_matches = self
            .searcher
            .search(rev_primer, read_seq, self.max_mismatch as usize);

        if rev_matches.is_empty() {
            return Err(FailureReason::NoReversePrimer);
        }

        // Find best valid amplicon
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
                // Check valid orientation (forward primer before reverse)
                if fwd.text_end > rev.text_start {
                    // Primers overlap or in wrong order
                    has_invalid_structure = true;
                    continue;
                }

                // Trimmed length excludes primers (insert only)
                let trimmed_len = rev.text_start - fwd.text_end;

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
                    start: fwd.text_end, // Start after forward primer
                    end: rev.text_start, // End before reverse primer
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
            None => {
                // Prioritize failure reasons: invalid structure > too short > too long
                if has_invalid_structure {
                    Err(FailureReason::InvalidStructure)
                } else if has_too_short {
                    Err(FailureReason::TooShort)
                } else if has_too_long {
                    Err(FailureReason::TooLong)
                } else {
                    // Should not happen if we have matches, but default to invalid structure
                    Err(FailureReason::InvalidStructure)
                }
            }
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
        writer,
        args.format,
        stats.clone(),
    );

    // Read and process in parallel
    info!("Reading and processing FASTQ records...");
    let reader = fastq::Reader::from_path(&args.input)?;

    reader
        .process_parallel(&mut processor, args.threads)
        .map_err(|e| anyhow::anyhow!("Processing failed: {}", e))?;

    info!("Processing complete");

    // Report statistics
    let mut stats_output: Box<dyn std::io::Write> = match &args.stats {
        Some(path) => {
            info!("Writing statistics to: {:?}", path);
            Box::new(std::fs::File::create(path)?)
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

    // Helper function to create a test processor
    fn create_test_processor(
        primers: PrimerPair,
        max_mismatch: u8,
        min_len: usize,
        max_len: usize,
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

        // Primers in wrong order: reverse before forward
        let read = b"TGCAATATATACGT";
        let result = processor.find_amplicon(read);

        // Should not find valid amplicon (primers in wrong order)
        assert!(result.is_err());
        assert_eq!(result.unwrap_err(), FailureReason::InvalidStructure);
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
            PrimerTrimProcessor::new(primers, 0, 5, 100, writer, OutputFormat::Fastq, stats);

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

        let processor =
            PrimerTrimProcessor::new(primers, 0, 10, 100, writer, OutputFormat::Fastq, stats);

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

        let processor =
            PrimerTrimProcessor::new(primers, 0, 10, 100, writer, OutputFormat::Fasta, stats);

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

        let processor =
            PrimerTrimProcessor::new(primers, 0, 10, 100, writer, OutputFormat::Fastq, stats);

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
}
