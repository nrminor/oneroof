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
    resplice_primers.py -i ${bed_file}
    """

}
