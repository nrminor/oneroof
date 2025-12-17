/*
 * PRIMER_HANDLING: Amplicon-based read processing workflow
 *
 * This subworkflow handles primer preparation, amplicon finding/trimming,
 * and per-sample read merging for amplicon sequencing data.
 *
 * Supports two input modes:
 *   1. BED + Reference FASTA: Full validation, spike-in resplicing, sequence extraction
 *   2. Primer TSV: Direct use of pre-computed primer sequences
 *
 * The consolidated PREPARE_PRIMERS process replaces what was previously 8 separate
 * processes (VALIDATE_PRIMER_BED, RESPLICE_PRIMERS, SPLIT_PRIMER_COMBOS,
 * GET_PRIMER_SEQS, CREATE_PRIMER_TSV, COLLECT_PRIMER_TSV, GET_PRIMER_PATTERNS).
 */

include {
    PREPARE_PRIMERS ;
    PREPARE_PRIMERS_FROM_TSV
} from "../modules/prepare_primers"
include { REPORT_REFERENCES        } from "../modules/reporting"
include { FIND_AND_TRIM_AMPLICONS  } from "../modules/find_and_trim_amplicons"
include { FAIDX                    } from "../modules/samtools"
include {
    MERGE_BY_SAMPLE ;
    AMPLICON_STATS
} from "../modules/seqkit"
include { READ_DOWNSAMPLING        } from "../modules/rasusa"
include { CREATE_AMPLICON_TSV      } from "../modules/primer_patterns"


workflow PRIMER_HANDLING {
    take:
    ch_basecalls // tuple(sample_id, reads_path)
    ch_primer_bed // path to primer BED file (optional)
    ch_refseq // path to reference FASTA (required if using BED)
    ch_primer_tsv // path to primer TSV file (optional, alternative to BED)

    main:

    // Filter out samples with too few reads
    ch_filtered_basecalls = ch_basecalls
        .map { sample_id, reads -> tuple(sample_id, file(reads), file(reads).countFasta()) }
        .filter { item -> item[2] > 10 }
        .map { sample_id, reads, _count -> tuple(sample_id, file(reads)) }

    // Prepare primers based on input type
    if (params.primer_bed && params.primer_bed != "") {
        // BED + FASTA input: full validation and resplicing
        PREPARE_PRIMERS(ch_primer_bed, ch_refseq)

        ch_primer_pairs = PREPARE_PRIMERS.out.primer_pairs
        ch_respliced_bed = PREPARE_PRIMERS.out.respliced_bed

        // Copy reference assets for provenance
        REPORT_REFERENCES(ch_primer_bed.mix(ch_refseq))
    }
    else if (params.primer_tsv && params.primer_tsv != "") {
        // TSV input: validation and passthrough
        PREPARE_PRIMERS_FROM_TSV(ch_primer_tsv)

        ch_primer_pairs = PREPARE_PRIMERS_FROM_TSV.out.primer_pairs
        ch_respliced_bed = Channel.empty()
    }

    // Split primer pairs TSV into channel of amplicon tuples
    // Each row becomes: (amplicon_name, fwd_sequence, rev_sequence)
    ch_amplicons = ch_primer_pairs
        .splitCsv(header: true, sep: '\t', strip: true)
        .map { row ->
            assert row.amplicon_name : "Missing amplicon_name in primer pairs TSV"
            assert row.fwd_sequence : "Missing fwd_sequence for ${row.amplicon_name}"
            assert row.rev_sequence : "Missing rev_sequence for ${row.amplicon_name}"
            tuple(row.amplicon_name, row.fwd_sequence, row.rev_sequence)
        }

    // Combine samples × amplicons for parallel processing
    // Result: (sample_id, reads, amplicon_name, fwd_sequence, rev_sequence)
    ch_sample_amplicons = ch_filtered_basecalls
        .combine(ch_amplicons)
        .map { sample_id, reads, amplicon_name, fwd_seq, rev_seq ->
            tuple(sample_id, reads, amplicon_name, fwd_seq, rev_seq)
        }

    // Find and trim amplicons
    FIND_AND_TRIM_AMPLICONS(ch_sample_amplicons)

    // Filter empty outputs and index
    ch_trimmed = FIND_AND_TRIM_AMPLICONS.out
        .map { sample_id, fasta -> tuple(sample_id, fasta, file(fasta).countFasta()) }
        .filter { item -> item[2] > 0 }
        .map { sample_id, fasta, _count -> tuple(sample_id, file(fasta)) }

    FAIDX(ch_trimmed)

    // Downsample reads per amplicon
    READ_DOWNSAMPLING(FAIDX.out)

    // Compute per-amplicon statistics (grouped by sample)
    AMPLICON_STATS(READ_DOWNSAMPLING.out.groupTuple(by: 0))

    // Create amplicon summary TSV joining stats with positions
    CREATE_AMPLICON_TSV(
        AMPLICON_STATS.out.collect(),
        ch_primer_pairs,
    )

    // Merge all amplicons per sample
    MERGE_BY_SAMPLE(READ_DOWNSAMPLING.out.groupTuple(by: 0))

    emit:
    merged_reads     = MERGE_BY_SAMPLE.out // tuple(sample_id, merged_fasta)
    primer_pairs     = ch_primer_pairs // path to primer_pairs.tsv
    respliced_bed    = ch_respliced_bed // path to respliced.bed (empty if TSV input)
    amplicon_stats   = AMPLICON_STATS.out // per-sample amplicon statistics
    amplicon_summary = CREATE_AMPLICON_TSV.out.summary_tsv // amplicon_summary.tsv
}
