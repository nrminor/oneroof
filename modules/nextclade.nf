process CHECK_DATASET {

    errorStrategy 'ignore'

    output:
    val "${params.nextclade_dataset}"

    when:
    params.nextclade_dataset

    shell:
    '''
    dataset_options=$(nextclade datasets list)

    if grep -Fxq !{params.nextclade_dataset} <<< "$dataset_options"; then
        echo "!{params.nextclade_dataset} found in supported Nextclade datasets. It is safe to proceed."
    else
        echo "!{params.nextclade_dataset} not found in the list of supported Nextclade datasets."
        exit 1
    fi
    '''

}

process DOWNLOAD_DATASET {

    storeDir params.nextclade_cache

    input:
    val valid_dataset

    output:
    path "dataset_label/"

    script:
    dataset_label = params.nextclade_dataset.split("/")[-1]
    """
    nextclade dataset get \
    --name ${valid_dataset} \
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

    script:
    label = file(sequences.toString()).getSimpleName()
    """
    nextclade run \
    --input-dataset ${nextclade_dataset} \
    --output-all=${label}/ \
    ${sequences}
    """

}
