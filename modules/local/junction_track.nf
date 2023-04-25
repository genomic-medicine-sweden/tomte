process JUNCTION_TRACK {
    tag 'junction_track'
    label 'process_low'

    conda "conda-forge::python=3.8.3"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.8.3' :
        'quay.io/biocontainers/python:3.8.3' }"

    input:
    path(sj)

    output:
    path('*.bed')    , emit: bed
    path "versions.yml" , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    $baseDir/bin/sj2bed.py --sj ${sj} --output ${prefix}_junction.bed

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sj2bed: v1.0
    END_VERSIONS
    """

    stub:
    """
    touch rrna.gtf
    touch rrna.bed

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sj2bed: v1.0
    END_VERSIONS
    """
}
