process SPLIT_BED_BY_CHROM {
    tag "$meta.id"
    label 'process_single'

    conda "conda-forge::gawk=5.3.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ubuntu:22.04' :
        'nf-core/ubuntu:22.04' }"

    input:
    tuple val(meta), path(bed)

    output:
    path "*.bed",        emit: intervals
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    # Split BED file into per-chromosome files, keeping only standard chromosomes
    awk '\$1 ~ /^(chr[0-9]+|chr[XYM]|chrMT|[1-9][0-9]?|X|Y|MT?)$/ {print > \$1".bed"}' ${bed}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        awk: \$(awk --version 2>&1 | head -1 | sed 's/[^0-9.]//g; s/^\\.//; s/\\.$//; s/\\.\\.*/./g')
    END_VERSIONS
    """

    stub:
    """
    touch chr1.bed

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        awk: 5.3.0
    END_VERSIONS
    """
}
