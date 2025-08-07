process CREATE_AMPLICON_TSV {

    publishDir params.primer_handling, mode: 'copy'
    
    input:
    path amplicon_stats_files, stageAs: "stats_*.tsv"
    path bed_file
    
    output:
    path "amplicon_summary.tsv", emit: summary_tsv
    
    script:
    """
   amplicon_tsv.py --bed ${bed_file} --pattern "stats_*.tsv"
    """
}