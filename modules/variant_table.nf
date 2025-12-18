process COLLECT_FULL_VARIANT_TABLE {

    publishDir params.variants, mode: 'copy', overwrite: true

    errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }
    maxRetries 2

    input:
    path "input_dir/*"

    output:
    tuple path("full_variant_table.tsv"), path("full_variant_table.parquet")

    script:
    """
    collect_full_variant_table.py \\
    --input-dir input_dir \\
    --output full_variant_table \\
    --consensus-threshold ${params.min_consensus_freq}
    """
}
