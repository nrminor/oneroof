#!/usr/bin/env -S cargo +nightly -Zscript
---
[dependencies]
polars = { version = "0.44", features = ["lazy", "csv", "strings", "regex", "fmt"] }
clap = { version = "4.5", features = ["derive"] }
env_logger = "0.11"
log = "0.4"
itertools = "0.13"
regex = "1.11"
anyhow = "1.0"
---

//! Resplice primers - finds all possible combinations of spike-in primers in an amplicon scheme.
//!
//! This script reads a BED file containing PCR primer coordinates and generates all valid
//! primer pair combinations when spike-in primers are present. It handles forward and reverse
//! primers, assigns unique indices, and outputs a new BED file with all possible amplicon
//! combinations.
//!

use anyhow::Result;
use clap::Parser;
use itertools::Itertools;
use log::{debug, error, info, warn};
use polars::prelude::*;
use regex::Regex;
use std::collections::HashMap;
use std::path::{Path, PathBuf};

/// Command line arguments for primer resplicing
#[derive(Parser, Debug)]
#[command(
    name = "resplice_primers",
    about = "Finds all possible combinations of spike-in primers in an amplicon scheme"
)]
struct Args {
    /// BED file with one-off spike-in primers to be respliced into possible amplicons
    #[arg(short = 'i', long = "input_bed")]
    input_bed: PathBuf,

    /// Output prefix for final respliced amplicon BED file
    #[arg(short = 'o', long = "output_prefix", default_value = "respliced")]
    output_prefix: String,

    /// The suffix to be expected in the names for forward primers
    #[arg(short = 'f', long = "fwd_suffix", default_value = "_LEFT")]
    fwd_suffix: String,

    /// The suffix to be expected in the names for reverse primers
    #[arg(short = 'r', long = "rev_suffix", default_value = "_RIGHT")]
    rev_suffix: String,

    /// The symbol used to delimit the index of a spike-in primer
    #[arg(short = 'd', long = "idx_delim", default_value = "-")]
    idx_delim: String,

    /// The position where the primer spike-in index should be expected
    #[arg(short = 'p', long = "idx_position", default_value_t = -1)]
    idx_position: i32,

    /// Increase verbosity level
    #[arg(short = 'v', long = "verbose", action = clap::ArgAction::Count)]
    verbose: u8,
}

/// Check if a BED file exists at the given path
fn check_bed_existence(path: &Path) -> Result<()> {
    info!("Parsing input BED file at {:?}...", path);

    if !path.is_file() {
        let parent = path.parent().unwrap_or(Path::new("."));
        let bed_files: Vec<_> = std::fs::read_dir(parent)?
            .filter_map(|entry| entry.ok())
            .filter(|entry| {
                entry
                    .path()
                    .extension()
                    .and_then(|ext| ext.to_str())
                    .map(|ext| ext == "bed")
                    .unwrap_or(false)
            })
            .map(|entry| entry.path())
            .collect();

        if bed_files.is_empty() {
            error!(
                "The provided input BED file, {:?}, does not exist or is not a file. Aborting.",
                path
            );
        } else {
            error!(
                "The provided input BED file, {:?}, does not exist or is not a file. \
                Perhaps you meant one of the following files?\n{:#?}",
                path, bed_files
            );
        }
        anyhow::bail!("Input BED file not found");
    }

    Ok(())
}

/// Check primer names for multiple occurrences of the index delimiter symbol
fn check_idx_delims(df: &DataFrame, idx_delim: &str) -> Result<()> {
    let names = df.column("NAME")?.str()?;
    let mut warning_list = Vec::new();

    for opt_name in names.into_iter() {
        if let Some(name) = opt_name {
            let delim_count = name.matches(idx_delim).count();
            if delim_count > 1 {
                warning_list.push(name.to_string());
            }
        }
    }

    if !warning_list.is_empty() {
        warn!(
            "{} primer names contained more than one of the symbol used to delimit spike-in primer index, '{}', \
            which is a special value in this program. Unexpected behavior, such as merging multiple amplicons \
            into one or incorrect resplicing, is likely to occur.\n\n\
            To fix it, either use a different delimiter symbol than '{}' or modify the primer names so that \
            the symbol '{}' is only used to denote a spike-in index.\n\n\
            Below is the list of {} primers that generated this warning:\n{:#?}",
            warning_list.len(), idx_delim, idx_delim, idx_delim, warning_list.len(), warning_list
        );
    }

    Ok(())
}

/// Partition a BED DataFrame by amplicon name
fn partition_by_amplicon(
    df: DataFrame,
    fwd_suffix: &str,
    rev_suffix: &str,
) -> Result<Vec<DataFrame>> {
    // Create regex to remove spike-in indices
    let idx_regex = Regex::new(r"-\d+")?;

    let df = df
        .lazy()
        .with_column(col("NAME").alias("ORIG_NAME"))
        .with_column(
            col("NAME")
                .str()
                .replace_all(lit(idx_regex.to_string()), lit(""), false)
                .alias("NAME"),
        )
        .with_column(
            col("NAME")
                .str()
                .replace_all(lit(fwd_suffix), lit(""), false)
                .str()
                .replace_all(lit(rev_suffix), lit(""), false)
                .alias("Amplicon"),
        )
        .select([
            col("Ref"),
            col("Start Position"),
            col("Stop Position"),
            col("ORIG_NAME"),
            col("NAME"),
            col("INDEX"),
            col("SENSE"),
            col("Amplicon"),
        ])
        .collect()?;

    // Partition by amplicon
    let amplicon_col = df.column("Amplicon")?.str()?;
    let unique_amplicons = amplicon_col.unique()?;
    let amplicons = unique_amplicons.sort(false);

    let mut partitions = Vec::new();
    for amplicon_opt in amplicons.into_iter() {
        if let Some(amplicon) = amplicon_opt {
            let mask = df.column("Amplicon")?.str()?.equal(amplicon);
            let partition = df.filter(&mask)?;
            partitions.push(partition);
        }
    }

    Ok(partitions)
}

/// Assign new indices to primers within an amplicon
fn assign_new_indices(
    mut df: DataFrame,
    fwd_suffix: &str,
    rev_suffix: &str,
    idx_delim: &str,
) -> Result<DataFrame> {
    let height = df.height();
    let indices: Vec<String> = (1..=height).map(|i| i.to_string()).collect();

    let df = df.with_column(Series::new("index".into(), indices))?;

    // Build new NAME column

    let new_names = df
        .clone()
        .lazy()
        .select([when(col("NAME").str().contains(lit(fwd_suffix), false))
            .then(col("Amplicon") + lit(fwd_suffix) + lit(idx_delim) + col("index"))
            .otherwise(col("Amplicon") + lit(rev_suffix) + lit(idx_delim) + col("index"))
            .alias("NEW_NAME")])
        .collect()?
        .column("NEW_NAME")?
        .clone();

    Ok(df
        .drop("NAME")?
        .with_column(new_names.with_name("NAME".into()))?
        .select([
            "Ref",
            "Start Position",
            "Stop Position",
            "ORIG_NAME",
            "NAME",
            "INDEX",
            "SENSE",
            "Amplicon",
        ])?)
}

/// Normalize indices for all amplicons
fn normalize_indices(
    partitioned: Vec<DataFrame>,
    fwd_suffix: &str,
    rev_suffix: &str,
) -> Result<HashMap<String, DataFrame>> {
    let mut normalized_dfs = Vec::new();

    for primer_df in partitioned {
        // Partition by NAME to group forward/reverse primers
        let name_col = primer_df.column("NAME")?.str()?;
        let unique_names = name_col.unique()?;
        let names = unique_names.sort(false);

        let mut name_dfs = Vec::new();
        for name_opt in names.into_iter() {
            if let Some(name) = name_opt {
                let mask = primer_df.column("NAME")?.str()?.equal(name);
                let name_df = primer_df.filter(&mask)?;
                let corrected = assign_new_indices(name_df, fwd_suffix, rev_suffix, &"-")?;
                name_dfs.push(corrected);
            }
        }

        let concatenated = if name_dfs.len() > 1 {
            concat(
                name_dfs
                    .iter()
                    .map(|df| df.clone().lazy())
                    .collect::<Vec<_>>(),
                UnionArgs::default(),
            )?
            .collect()?
        } else {
            name_dfs[0].clone()
        };

        normalized_dfs.push(concatenated);
    }

    // Concatenate all and partition by amplicon into a map
    let all_normalized = if normalized_dfs.len() > 1 {
        concat(
            normalized_dfs
                .iter()
                .map(|df| df.clone().lazy())
                .collect::<Vec<_>>(),
            UnionArgs::default(),
        )?
        .collect()?
    } else {
        normalized_dfs[0].clone()
    };

    let amplicon_col = all_normalized.column("Amplicon")?.str()?;
    let amplicons = amplicon_col.unique()?;

    let mut amplicon_map = HashMap::new();
    for amplicon_opt in amplicons.into_iter() {
        if let Some(amplicon) = amplicon_opt {
            let mask = all_normalized.column("Amplicon")?.str()?.equal(amplicon);
            let partition = all_normalized.filter(&mask)?;
            amplicon_map.insert(amplicon.to_string(), partition);
        }
    }

    Ok(amplicon_map)
}

/// Check if a string can be converted to an integer
fn convertable_to_int(s: &str) -> bool {
    s.parse::<i32>().is_ok()
}

/// Resolve primer names to generate new labels for all combinations
fn resolve_primer_names(
    old_fwd_primers: Vec<&str>,
    old_rev_primers: Vec<&str>,
    idx_delim: &str,
    idx_position: i32,
) -> Result<(Vec<String>, Vec<String>)> {
    let all_pairs: Vec<_> = old_fwd_primers
        .iter()
        .cartesian_product(old_rev_primers.iter())
        .collect();

    let mut new_primer_pairs = Vec::new();
    let mut old_primer_pairs = Vec::new();

    for (i, (&fwd, &rev)) in all_pairs.iter().enumerate() {
        let fwd_parts: Vec<&str> = fwd.split(idx_delim).collect();
        let rev_parts: Vec<&str> = rev.split(idx_delim).collect();

        let idx_pos = if idx_position < 0 {
            fwd_parts.len() as i32 + idx_position
        } else {
            idx_position
        } as usize;

        let fwd_final = fwd_parts
            .get(idx_pos)
            .ok_or_else(|| anyhow::anyhow!("Invalid index position for forward primer"))?;
        let rev_final = rev_parts
            .get(idx_pos)
            .ok_or_else(|| anyhow::anyhow!("Invalid index position for reverse primer"))?;

        if !convertable_to_int(fwd_final) {
            error!(
                "The primer {}, which has been paired with {} does not end with a {}-delimited integer, \
                e.g. '-1', which is required for properly handling different possible primer combinations. \
                Printing the set of combinations being handled:\n{:#?}",
                fwd, rev, idx_delim, all_pairs
            );
            anyhow::bail!("Invalid primer format");
        }

        if !convertable_to_int(rev_final) {
            error!(
                "The primer {}, which has been paired with {} does not end with a {}-delimited integer, \
                e.g. '-1', which is required for properly handling different possible primer combinations. \
                Printing the set of combinations being handled:\n{:#?}",
                rev, fwd, idx_delim, all_pairs
            );
            anyhow::bail!("Invalid primer format");
        }

        let new_fwd =
            fwd.replace(&format!("{}{}", idx_delim, fwd_final), "") + &format!("_splice{}", i + 1);
        let new_rev =
            rev.replace(&format!("{}{}", idx_delim, rev_final), "") + &format!("_splice{}", i + 1);

        old_primer_pairs.push((fwd.to_string(), rev.to_string()));
        new_primer_pairs.push((new_fwd, new_rev));
    }

    // Flatten the pairs
    let old_names: Vec<String> = old_primer_pairs
        .iter()
        .flat_map(|(f, r)| vec![f.clone(), r.clone()])
        .collect();

    let new_names: Vec<String> = new_primer_pairs
        .iter()
        .flat_map(|(f, r)| vec![f.clone(), r.clone()])
        .collect();

    Ok((old_names, new_names))
}

/// Resplice primers to generate all possible combinations
fn resplice_primers(
    amplicon_dfs: HashMap<String, DataFrame>,
    fwd_suffix: &str,
    rev_suffix: &str,
    idx_delim: &str,
    idx_position: i32,
) -> Result<Vec<DataFrame>> {
    let mut mutated_frames = Vec::new();

    for (amplicon, primer_df) in amplicon_dfs {
        if primer_df.height() == 2 {
            // Standard two-primer amplicon
            let names_series = primer_df.column("NAME")?.str()?;
            let pair_labels: Vec<_> = names_series
                .into_iter()
                .filter_map(|opt| opt.map(|s| s.to_string()))
                .collect();

            debug!(
                "Pair of primers within the amplicon {} detected: {:?}. \
                No resplicing will be needed here, though double check that the necessary \
                forward and reverse suffixes, {} and {}, are present.",
                amplicon, pair_labels, fwd_suffix, rev_suffix
            );

            assert!(
                pair_labels.iter().any(|p| p.contains(fwd_suffix)),
                "The forward suffix {} is missing in the provided primer pairs: {:?}. Aborting.",
                fwd_suffix,
                pair_labels
            );

            assert!(
                pair_labels.iter().any(|p| p.contains(rev_suffix)),
                "The reverse suffix {} is missing in the provided primer pairs: {:?}. Aborting.",
                rev_suffix,
                pair_labels
            );

            mutated_frames.push(primer_df);
            continue;
        }

        if primer_df.height() == 1 {
            error!(
                "There is a single primer without an amplicon running around! \
                Here's the parsed BED file for this amplicon:\n{:?}",
                primer_df
            );
            anyhow::bail!("Single primer without pair found");
        }

        // Extract forward and reverse primers
        let names_series = primer_df.column("NAME")?.str()?;
        let primers: Vec<String> = names_series
            .into_iter()
            .filter_map(|opt| opt.map(|s| s.to_string()))
            .collect();

        let fwd_primers: Vec<&str> = primers
            .iter()
            .filter(|p| p.contains(fwd_suffix))
            .map(|s| s.as_str())
            .collect();

        let rev_primers: Vec<&str> = primers
            .iter()
            .filter(|p| p.contains(rev_suffix))
            .map(|s| s.as_str())
            .collect();

        debug!(
            "Resplicing will be performed for the following primers:\n{:#?}",
            primers
        );

        let (old_primer_names, new_primer_names) =
            resolve_primer_names(fwd_primers, rev_primers, idx_delim, idx_position)?;

        assert_eq!(
            old_primer_names.len(),
            new_primer_names.len(),
            "Insufficient number of replacement names generated for amplicon {}: {:?}",
            amplicon,
            primer_df
        );

        // Create mapping dataframe
        let mapping_df = DataFrame::new(vec![Series::new("NAME".into(), old_primer_names).into()])?;

        // Join with original data
        let mut joined = mapping_df.join(
            &primer_df,
            ["NAME"],
            ["NAME"],
            JoinArgs::new(JoinType::Left),
        )?;

        // Replace NAME column with new names
        let new_names_series = Series::new("NAME".into(), new_primer_names);
        let new_df = joined.with_column(new_names_series)?.select([
            "Ref",
            "Start Position",
            "Stop Position",
            "ORIG_NAME",
            "NAME",
            "INDEX",
            "SENSE",
            "Amplicon",
        ])?;

        mutated_frames.push(new_df);
    }

    Ok(mutated_frames)
}

/// Finalize primer pairings and remove invalid combinations
fn finalize_primer_pairings(
    mutated_frames: Vec<DataFrame>,
    fwd_suffix: &str,
    rev_suffix: &str,
) -> Result<DataFrame> {
    let mut final_frames = Vec::new();

    for df in mutated_frames {
        let names_series = df.column("NAME")?.str()?;
        let names: Vec<String> = names_series
            .into_iter()
            .filter_map(|opt| opt.map(|s| s.to_string()))
            .collect();

        let fwd_keepers: Vec<_> = names.iter().filter(|n| n.contains(fwd_suffix)).collect();

        let rev_keepers: Vec<_> = names.iter().filter(|n| n.contains(rev_suffix)).collect();

        if fwd_keepers.is_empty() || rev_keepers.is_empty() {
            warn!(
                "Incorrect splicing occurred for the following sets of primers. \
                The amplicon they are derived from will be skipped:\n{:?}\n{:?}",
                fwd_keepers, rev_keepers
            );
            continue;
        }

        final_frames.push(df);
    }

    let final_df = if final_frames.len() > 1 {
        concat(
            final_frames
                .iter()
                .map(|df| df.clone().lazy())
                .collect::<Vec<_>>(),
            UnionArgs::default(),
        )?
        .collect()?
    } else {
        final_frames[0].clone()
    };

    Ok(final_df)
}

fn main() -> Result<()> {
    let args = Args::parse();

    // Set up logging
    let log_level = match args.verbose {
        0 => "warn",
        1 => "info",
        2 => "debug",
        _ => "trace",
    };

    env_logger::Builder::from_env(env_logger::Env::default().default_filter_or(log_level)).init();

    info!("Logging set to {}", log_level);

    // Check that input file exists
    check_bed_existence(&args.input_bed)?;

    // Read BED file
    let mut parse_options = CsvParseOptions::default();
    parse_options.separator = b'\t';

    let bed_df = CsvReadOptions::default()
        .with_has_header(false)
        .with_parse_options(parse_options)
        .try_into_reader_with_file_path(Some(args.input_bed.clone()))?
        .finish()?
        .lazy()
        .rename(
            [
                "column_1", "column_2", "column_3", "column_4", "column_5", "column_6",
            ],
            [
                "Ref",
                "Start Position",
                "Stop Position",
                "NAME",
                "INDEX",
                "SENSE",
            ],
            false,
        )
        .collect()?;

    // Check for multiple delimiter occurrences
    check_idx_delims(&bed_df, &args.idx_delim)?;

    // Partition by amplicon
    let partitioned_bed =
        partition_by_amplicon(bed_df.clone(), &args.fwd_suffix, &args.rev_suffix)?;
    info!(
        "{} input primers split into {} discrete amplicons.",
        bed_df.height(),
        partitioned_bed.len()
    );

    if bed_df.height() % partitioned_bed.len() != 0 {
        info!(
            "There appear to be {} spike-in primers.",
            bed_df.height() % partitioned_bed.len()
        );
    } else {
        info!(
            "There don't appear to be any spike-in primers. Assuming that amplicons can be parsed properly, \
            this means that no resplicing will occur, though primer indices will still be updated."
        );
    }

    // Normalize indices
    info!(
        "Re-assigning 1-based indices to account for current or future spike-ins amongst the {} input primers.",
        bed_df.height()
    );
    let indexed_primer_dfs =
        normalize_indices(partitioned_bed, &args.fwd_suffix, &args.rev_suffix)?;

    info!(
        "Re-indexing successful. Proceeding to resplicing any implicit combinations \
        that are not accounted for with a row in the input BED file."
    );

    // Resplice primers
    let respliced_dfs = resplice_primers(
        indexed_primer_dfs.clone(),
        &args.fwd_suffix,
        &args.rev_suffix,
        &args.idx_delim,
        args.idx_position,
    )?;

    if respliced_dfs.len() != indexed_primer_dfs.len() {
        warn!(
            "The number of amplicon primer sets that made it through resplicing, {}, \
            does not match the number of input amplicons, {}. Data loss may have occurred for this primer set.",
            respliced_dfs.len(), indexed_primer_dfs.len()
        );
    }

    // Finalize primer pairings
    info!(
        "Finalizing primer pairings amongst {} amplicons by verifying that \
        every forward primer comes with a reverse primer.",
        respliced_dfs.len()
    );
    let final_df = finalize_primer_pairings(respliced_dfs, &args.fwd_suffix, &args.rev_suffix)?;

    info!(
        "The input BED file contained {} primers. After searching for any spike-ins to resplice, \
        the final BED file contains {}",
        bed_df.height(),
        final_df.height()
    );

    // Drop unnecessary columns and sort
    let output = final_df
        .lazy()
        .drop(["Amplicon", "ORIG_NAME"])
        .sort_by_exprs(
            vec![col("Ref"), col("Start Position"), col("Stop Position")],
            SortMultipleOptions::default(),
        )
        .collect()?;

    // Write output
    let output_path = format!("{}.bed", args.output_prefix);
    let mut file = std::fs::File::create(&output_path)?;
    CsvWriter::new(&mut file)
        .include_header(false)
        .with_separator(b'\t')
        .finish(&mut output.clone())?;

    info!(
        "Reading and processing of primers in '{:?}' was successful. \
        The output BED file was written to '{}'. Goodbye!",
        args.input_bed, output_path
    );

    Ok(())
}
