process BOOTSTRAPANN {
    tag '$meta.id'
    label 'process_low'

    conda "conda-forge::python=2.7"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:2.7' :
        'quay.io/biocontainers/python:2.7' }"

    input:
    tuple val(meta), path(vcf), path(tbi)
    tuple val(meta), path(csv)

    output:
    tuple val(meta), path("*.vcf")

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    $baseDir/bin/BootstrapAnn.py \\
        --vcf ${vcf} \\
        --ase ${csv} \\
        > ${prefix}_ase.vcf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_ase.vcf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
