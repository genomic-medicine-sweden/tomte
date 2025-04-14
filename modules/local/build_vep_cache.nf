process BUILD_VEP_CACHE {
    tag "build_vep_cache"
    label 'process_short'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bcftools:1.20--h8b25389_0':
        'biocontainers/bcftools:1.20--h8b25389_0' }"

    input:
    path(vep_reference_files)
    path(gnomad_files)

    output:
    path("vep_cache")           , emit: vep_cache
    path("vep_plugin_files.csv"), emit: vep_plugin_file
    path "versions.yml"         , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    # Create file listing all vep plugins to use. Note that file extension must be .csv

    echo "vep_files" > vep_plugin_files.csv
    mkdir vep_cache
    mkdir vep_cache/vep_plugins

    for file_path in ${vep_reference_files} ${gnomad_files}; do
        file_name=\$(basename "\$file_path")
        if [ -d "\${file_path}" ]; then
            # If it's a directory, copy it to vep_cache
            cp -r "\${file_path}" vep_cache/.
            echo "vep_cache/\${file_name}" >> vep_plugin_files.csv
        elif [ -f "\${file_path}" ]; then
            # If it's a file, copy it to vep_cache/vep_plugins
            cp "\${file_path}" vep_cache/vep_plugins/.
            echo "vep_cache/vep_plugins/\${file_name}" >> vep_plugin_files.csv
        fi
    done

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        wget: \$(echo wget -V 2>&1 | grep "GNU Wget" | cut -d" " -f3 )
    END_VERSIONS
    """

    stub:
    """
    mkdir vep_cache
    touch vep_plugin_files.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        wget: \$(echo wget -V 2>&1 | grep "GNU Wget" | cut -d" " -f3 )
    END_VERSIONS
    """

}
