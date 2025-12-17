process PREPARE_PRIMERS {
    /*
    Consolidated primer preparation process.
    
    Replaces: VALIDATE_PRIMER_BED, RESPLICE_PRIMERS, SPLIT_PRIMER_COMBOS,
              GET_PRIMER_SEQS, CREATE_PRIMER_TSV, COLLECT_PRIMER_TSV, GET_PRIMER_PATTERNS
    
    Handles spike-in primer resplicing and sequence extraction in a single process.
    */

    publishDir params.primer_handling, mode: 'copy', overwrite: true

    errorStrategy 'finish'

    input:
    path primer_bed
    path reference

    output:
    path "primer_pairs.tsv", emit: primer_pairs
    path "respliced.bed", emit: respliced_bed

    script:
    """
    prepare_primers.py \\
        --input-bed ${primer_bed} \\
        --reference ${reference} \\
        --output-tsv primer_pairs.tsv \\
        --output-bed respliced.bed \\
        --fwd-suffix ${params.fwd_suffix} \\
        --rev-suffix ${params.rev_suffix}
    """
}


process PREPARE_PRIMERS_FROM_TSV {
    /*
    Primer preparation from existing TSV (bypass BED processing).
    
    Used when users provide a pre-computed primer_pairs.tsv instead of BED + FASTA.
    */

    publishDir params.primer_handling, mode: 'copy', overwrite: true

    errorStrategy 'finish'

    input:
    path primer_tsv

    output:
    path "primer_pairs.tsv", emit: primer_pairs

    script:
    """
    prepare_primers.py \\
        --input-tsv ${primer_tsv} \\
        --output-tsv primer_pairs.tsv
    """
}
