process WGET_DOWNLOAD {
    tag "wget_gunzip"
    label 'process_low'

    conda "bioconda::gnu-wget=1.18"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gnu-wget:1.18--h5bf99c6_5' :
        'quay.io/biocontainers/gnu-wget:1.18--h5bf99c6_5' }"

    input:
    tuple val(meta), val(path_wget)

    output:
    path(meta), emit: downloaded_file
    path "versions.yml"  , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def filename = path_wget.split('/').last()
    """

    wget $path_wget

    if [[ $filename == *.tar.gz ]]; then
        tar -xzf $filename
    else
        mv $filename $meta
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        wget: \$(echo wget -V 2>&1 | grep "GNU Wget" | cut -d" " -f3 > versions.yml)
    END_VERSIONS
    """

    stub:
    """
    touch $meta
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        wget: \$(echo wget -V 2>&1 | grep "GNU Wget" | cut -d" " -f3 > versions.yml)
    END_VERSIONS
    """

}
