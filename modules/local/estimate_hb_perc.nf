process ESTIMATE_HB_PERC {
    tag "estimate_hb_perc"
    label "process_low"

    conda "conda-forge::python=3.8.3"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.8.3' :
        'biocontainers/python:3.8.3' }"

    input:
    tuple val(meta), path(gene_counts)
    path hb_genes

    output:
    tuple val(meta), path("${meta.sample}_perc_mapping.json"), emit: tsv
    path("versions.yml"), emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    calculate_perc_mapping.py \\
        --target_genes "${hb_genes}" \\
        --gene_counts "${gene_counts}" \\
        --strandedness "${meta.strandedness}" \\
        --out_json "${meta.sample}_perc_mapping.json"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        calculate_perc_mapping: \$(calculate_perc_mapping.py --version)
    END_VERSIONS
    """

    stub:
    """
    touch "${meta.sample}_perc_mapping.json"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        calculate_perc_mapping: \$(calculate_perc_mapping.py --version)
    END_VERSIONS
    """
}
