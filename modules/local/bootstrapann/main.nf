process BOOTSTRAPANN {
    tag "$meta.id"
    label 'process_low'

    conda "bioconda::pydamage=0.70"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pydamage:0.70--pyhdfd78af_0' :
        'biocontainers/pydamage:0.70--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(vcf), path(tbi), path(csv)

    output:
    tuple val(meta), path("*.vcf"), emit: vcf
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    BootstrapAnn.py \\
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
