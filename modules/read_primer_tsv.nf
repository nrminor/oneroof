process READ_PRIMER_TSV{
input:
    tuple val(name), val(sequences)

output:
    path "${name}.fasta"

script:
"""
    echo -e "${sequences}" > ${name}.fasta
    
"""

}