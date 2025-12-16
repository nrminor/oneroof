process GET_PRIMER_PATTERNS {

    /* */

	errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }
	maxRetries 2

	input:
	path primer_fasta

	output:
	path "${primer_combo}.txt"

	script:
	primer_combo = file(primer_fasta.toString()).getSimpleName()
	"""
	make_primer_patterns.py \
	-i ${primer_combo}.fasta \
	-o ${primer_combo} \
	-f "" \
	-r ""
	"""

}

process WRITE_PRIMER_FASTA {

	errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }
	maxRetries 2

    input:
    tuple val(name), val(sequences)

    output:
    path "${name}.fasta"

    script:
    """
    echo -e "${sequences}" > ${name}.fasta
    """
}

process CREATE_AMPLICON_TSV {
    publishDir params.primer_handling, mode: 'copy'

	errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }
	maxRetries 2

    input:
    path amplicon_stats_files, stageAs: "stats_*.tsv"
    path bed_file

    output:
    path "amplicon_summary.tsv", emit: summary_tsv

    script:
    """
    create_amplicon_tsv.py --bed ${bed_file} --pattern "stats_*.tsv"
    """
}

process CREATE_PRIMER_TSV {

    errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }
    maxRetries 2

    input:
        path primer_seq_fasta

    output:
        path "*.tsv", emit: primer_pair

    script:
    def uuid = UUID.randomUUID().toString()
    """
    amplicon_name=\$(basename "$primer_seq_fasta" | sed 's/\\.[^.]*\$//')

    fwd=\$(awk 'NR==2' "$primer_seq_fasta")
    rev=\$(awk 'NR==4' "$primer_seq_fasta")

    echo -e "\$amplicon_name\t\$fwd\t\$rev" > "${uuid}.tsv"
    """
}

process COLLECT_PRIMER_TSV {
    publishDir params.primer_handling, mode: 'copy', overwrite: true

    errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }
    maxRetries 2

    input:
    path all_primers

    output:
    path "primer_pairs.tsv", emit: primer_pairs

    script:
    """
    awk -F"\t" 'BEGIN {print "amplicon_name\tfwd_sequence\treverse_sequence"} \$1 != "amplicon_name" {print \$0}' *.tsv > primer_pairs.tsv
    """
}
