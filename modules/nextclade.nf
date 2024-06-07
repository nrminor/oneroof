process DOWNLOAD_DATASET {

    storeDir params.nextclade_cache

    output:
    path "dataset_label/"

    when:
    params.dev == false

    script:
    dataset_label = params.nextclade_dataset.split("/")[-1]
    """
    nextclade dataset get \
    --name ${params.nextclade_cache} \
    --output-dir ${dataset_label}
    """

}

process RUN_NEXTCLADE {

    publishDir params.nextclade, mode: 'copy'

    input:
    path sequences
    each path("${nextclade_dataset}") 

    output:
    path "${label}/"

    when:
    params.dev == false

    script:
    label = file(sequences.toString()).getSimpleName()
    """
    nextclade run \
    --input-dataset ${nextclade_dataset} \
    --output-all=${label}/ \
    ${sequences}
    """

}