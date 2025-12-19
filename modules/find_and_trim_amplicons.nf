process FIND_AND_TRIM_AMPLICONS {
    /*
    Find and trim amplicon sequences from reads using primer sequences.
    
    Takes primer sequences directly as values (from primer_pairs.tsv),
    eliminating the need for intermediate pattern files.
    */

    tag "${sample_id}, ${amplicon_name}"

    errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }
    maxRetries 2

    cpus 4

    input:
    tuple val(sample_id), path(reads), val(amplicon_name), val(fwd_sequence), val(rev_sequence)

    output:
    tuple val(sample_id), path("${sample_id}.${amplicon_name}.fasta.gz")

    script:
    """
    export RUST_LOG=find_and_trim_amplicons=info

    # Use pre-compiled binary in container, fall back to rust-script locally
    if command -v find_and_trim_amplicons &> /dev/null; then
        TRIM_CMD=find_and_trim_amplicons
    else
        TRIM_CMD=find_and_trim_amplicons.rs
    fi

    \$TRIM_CMD \\
        --input ${reads} \\
        --forward ${fwd_sequence} \\
        --reverse ${rev_sequence} \\
        --max-mismatch ${params.max_mismatch} \\
        --min-len ${params.min_len} \\
        --max-len ${params.max_len} \\
        --forward-window ${params.forward_window} \\
        --reverse-window ${params.reverse_window} \\
        ${params.no_primer_trim ? '--no-trim' : ''} \\
        --threads ${task.cpus} \\
        --format fasta \\
        --no-compress \\
        --stats stats.txt \\
        --output /dev/stdout \\
    | bgzip -c > ${sample_id}.${amplicon_name}.fasta.gz
    """
}
