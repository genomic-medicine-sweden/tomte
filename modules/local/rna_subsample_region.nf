process RNA_DOWNSAMPLE {
    tag "$meta.id"
    label 'process_low'

    conda "bioconda::samtools=1.16.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/samtools:1.16.1--h6899075_1' :
        'quay.io/biocontainers/samtools:1.16.1--h6899075_1' }"

    input:
    tuple val(meta), path(bam)
    path subsample_bed

    output:
    tuple val(meta), path("*_subsamp.bam")                                , emit: bam
    tuple val(meta), path("*_subsamp.bam.bai")                            , emit: bai
    tuple val(meta), path("*_subsamp.bam"), path("*_subsamp.bam.bai") , emit: bam_bai
    path  "versions.yml"                                                      , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def seed_frac       = (args.contains('-s')) ? '' : '-s 0.001'

    """
    samtools view -@ $task.cpus -b -U non_select.bam -L ${subsample_bed} ${bam} | samtools view $seed_frac -@ $task.cpus -b -o select.bam
    samtools merge -u non_select.bam select.bam -o ${prefix}_subsamp.bam
    samtools index ${prefix}_subsamp.bam
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix    ?: "${meta.id}"
    """
    touch ${prefix}.bam
    touch ${prefix}.bam.bai

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """
}
