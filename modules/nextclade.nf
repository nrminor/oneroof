process DOWNLOAD_DATASET {

    storeDir params.nextclade_cache

    output:
    path "dataset_label/"

    script:
    dataset_label = params.nextclade_dataset.split("/")[-1]
    """
    nextclade dataset get \
    --name ${params.nextclade_cache} \
    --output-dir ${dataset_label}
    """

}

process NEXTCLADE_RUN {

    publishDir params.nextclade, mode: 'copy'

    input:
    path sequences
    each path("${nextclade_dataset}") 

    output:
    path "${label}/"

    script:
    label = file(sequences.toString()).getSimpleName()
    """
    nextclade run \
    --input-dataset ${nextclade_dataset} \
    --output-all=${label}/ \
    ${sequences}
    """

}