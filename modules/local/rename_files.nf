process RENAME_FILES {
    tag "$meta.id"
    label 'process_single'

    conda "conda-forge::coreutils=8.31"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gnu-wget:1.18--h36e9172_9' :
        'biocontainers/gnu-wget:1.18--h36e9172_9' }"

    input:
    tuple val(meta), path(input)

    output:
    tuple val(meta), path("${output}"), emit: output
    path "versions.yml",                emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def simplename = input.getSimpleName()
    output = input.toString().replace("${simplename}","${prefix}")
    """
    ln -s $input ${output}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ln: \$(echo \$(ln --version 2>&1 | head -n 1 | cut -d ' ' -f4))
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def simplename = input.getSimpleName()
    output = input.toString().replace("${simplename}","${prefix}")
    """
    touch ${output}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ln: \$(echo \$(ln --version 2>&1 | head -n 1 | cut -d ' ' -f4))
    END_VERSIONS
    """
}
