process RESPLICE_PRIMERS {

    /*
    */

    publishDir params.respliced, mode: 'copy', overwrite: true

    errorStrategy 'finish'

    input:
    path bed_file

    output:
    path "respliced.bed"

    script:
    """
    resplice_primers.py \
    --input_bed ${bed_file} \
    --fwd_suffix ${params.fwd_suffix} \
    --rev_suffix ${params.rev_suffix}
    """

}
