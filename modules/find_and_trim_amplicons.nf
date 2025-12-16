process FIND_AND_TRIM_AMPLICONS {

	tag "${barcode}, ${amplicon}"

	errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }
	maxRetries 2

	cpus 4

    input:
	tuple path(reads), path(patterns)

    output:
    tuple val(barcode), path("${barcode}.${amplicon}.fasta.gz")

    script:
	barcode = file(reads).getSimpleName()
	amplicon = file(patterns).getSimpleName()
    """
	export RUST_LOG=find_and_trim_amplicons=info 

    FORWARD_PATTERN=\$(head -n 1 ${patterns})
    REVERSE_PATTERN=\$(tail -n 1 ${patterns})

    find_and_trim_amplicons.rs \\
    --input ${reads} \\
    --forward \${FORWARD_PATTERN} \\
    --reverse \${REVERSE_PATTERN} \\
    --max-mismatch ${params.max_mismatch} \\
    --min-len ${params.min_len} \\
    --max-len ${params.max_len} \\
    --forward-window ${params.forward_window} \\
    --reverse-window ${params.reverse_window} \\
    --threads ${task.cpus} \\
    --format fasta \\
	--no-compress \\
    --stats stats.txt \\
	--output /dev/stdout \\
	| bgzip -c > ${barcode}.${amplicon}.fasta.gz
    """

}
