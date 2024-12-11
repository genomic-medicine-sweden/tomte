process GATK4_ASEREADCOUNTER {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gatk4:4.5.0.0--py36hdfd78af_0':
        'biocontainers/gatk4:4.5.0.0--py36hdfd78af_0' }"

    input:
    tuple val(meta),  path(bam), path(bai), path(vcf), path(tbi)
    tuple val(meta2), path(fasta)
    tuple val(meta3), path(fai)
    tuple val(meta4), path(dict)
    path intervals

    output:
    tuple val(meta), file("*_ase.csv"), emit: csv
    path "versions.yml"               , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def reference_command = fasta ? "--reference $fasta" : ""
    def dictionary_command = fasta ? "--sequence-dictionary $dict" : ""
    def intervals_command = intervals ? "--intervals $intervals" : ""

    def avail_mem = 3072
    if (!task.memory) {
        log.info '[GATK ASEReadCounter] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.'
    } else {
        avail_mem = (task.memory.mega*0.8).intValue()
    }
    """
    gatk --java-options "-Xmx${avail_mem}M -XX:-UsePerfData" \\
        ASEReadCounter \\
        --output ${prefix}_ase.csv \\
        --input ${bam} \\
        --variant ${vcf} \\
        $reference_command \\
        $intervals_command \\
        $dictionary_command \\
        --tmp-dir . \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    touch ${prefix}_ase.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
    END_VERSIONS
    """
}
