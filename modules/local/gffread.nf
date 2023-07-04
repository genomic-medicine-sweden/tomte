process GFFREAD {
    tag "$gff"
    label 'process_low'

    conda "bioconda::gffread=0.12.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gffread:0.12.1--h8b12597_0' :
        'biocontainers/gffread:0.12.1--h8b12597_0' }"

    input:
    tuple val(meta), path(gff)
    tuple val(meta), path(fasta)

    output:
    path "*.gtf"        , emit: gtf,        optional: true
    path "*.fa"         , emit: tr_fasta,   optional: true
    path "versions.yml" , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args   ?: ''
    def genome = fasta ? "-g ${fasta}" : ''
    def prefix = task.ext.prefix ?: gff.getSimpleName()
    def outp   = args.contains("-w ") ? "" : "-o ${prefix}.gtf"
    """
    gffread \\
        $gff \\
        $genome \\
        $args \\
        $outp
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gffread: \$(gffread --version 2>&1)
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: gff.getSimpleName()
    """
    touch ${prefix}.gtf
    touch ${prefix}.fa

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gffread: \$(gffread --version 2>&1)
    END_VERSIONS
    """
}
