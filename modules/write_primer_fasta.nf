process WRITE_PRIMER_FASTA {
    input:
    tuple val(name), val(sequences)

    output:
    path "${name}.fasta"

    script:
    """
    echo -e "${sequences}" > ${name}.fasta
    
    """
}
