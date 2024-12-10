process ESTIMATE_HB_PERC {

    tag "${meta.sample}"
    label "process_low"
    container "${params.containers.base}"

    input:
    tuple val(meta), path(gene_counts)
    path hb_genes

    output:
    path ("${meta.sample}_perc_mapping.json"), emit: tsv

    script:
    """
    calculate_perc_mapping.py \\
        --target_genes "${hb_genes}" \\
        --gene_counts "${gene_counts}" \\
        --strandedness "${meta.strandedness}" \\
        --out_json "${meta.sample}_perc_mapping.json"
    """

    stub:
    """
    touch "${meta.sample}_perc_mapping.json"
    """
}
