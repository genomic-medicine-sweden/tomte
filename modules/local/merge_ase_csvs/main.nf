process MERGE_ASE_CSVS {
    tag "$meta.id"
    label 'process_single'

    conda "conda-forge::gawk=5.3.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ubuntu:22.04' :
        'nf-core/ubuntu:22.04' }"

    input:
    tuple val(meta), path(csvs)

    output:
    tuple val(meta), path("*_ase.csv"), emit: csv
    path "versions.yml",                emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    # Merge CSVs: keep header from first file, skip header in subsequent files
    awk 'FNR==1 && NR>1 {next} {print}' ${csvs} > ${prefix}_ase.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        awk: \$(awk --version 2>&1 | head -1 | sed 's/[^0-9.]//g; s/^\\.//; s/\\.$//; s/\\.\\.*/./g')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_ase.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        awk: 5.3.0
    END_VERSIONS
    """
}
