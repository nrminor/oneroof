process GENERATE_AMPLICON_SUMMARY {
    /*
    Summarize amplicon coverage by joining stats with primer position data.
    
    Uses primer_pairs.tsv (preferred) which already contains amplicon positions,
    avoiding the need to re-parse BED files.
    */

    publishDir params.primer_handling, mode: 'copy'

    errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }
    maxRetries 2

    input:
    path amplicon_stats_files, stageAs: "stats_*.tsv"
    path primer_pairs_tsv

    output:
    path "amplicon_summary.tsv", emit: summary_tsv

    script:
    """
    generate_amplicon_summary.py --primer-tsv ${primer_pairs_tsv} --pattern "stats_*.tsv"
    """
}
