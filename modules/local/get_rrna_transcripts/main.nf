process GET_RRNA_TRANSCRIPTS {
    tag 'get_rrna_bed'
    label 'process_low'

    conda "bioconda::pirate=1.0.4 bioconda::perl-bioperl=1.7.2"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pirate:1.0.4--hdfd78af_2' :
        'biocontainers/pirate:1.0.4--hdfd78af_2' }"

    input:
    tuple val(meta), path(gtf)

    output:
    path('rrna.gtf')    , emit: rrnagtf
    path('rrna.bed')    , emit: bed
    path "versions.yml" , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    get_rrna_transcripts.py --gtf ${gtf} --output rrna.gtf

    gtf2bed rrna.gtf > rrna.bed

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        get_rrna_transcripts: v1.0
    END_VERSIONS
    """

    stub:
    """
    touch rrna.gtf
    touch rrna.bed

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        get_rrna_transcripts: v1.0
    END_VERSIONS
    """
}
