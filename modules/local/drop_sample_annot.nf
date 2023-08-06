process DROP_ANNOTATION {
    tag "DROP_annot"
    label 'process_low'

    conda "bioconda::drop=1.3.3"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/drop:1.3.3--pyhdfd78af_0' :
        'biocontainers/drop:1.3.3--pyhdfd78af_0' }"

    input:
    path(processed_gene_counts)
    path(gtf)

    output:
    path('sample_annotation.tsv'), emit: sample_annotation_drop
    path "versions.yml"              , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def gtf_name   = gtf ? gtf.getBaseName() : ""

    """
    drop_sample_annot.py \\
        --count_file $processed_gene_counts \\
        --gtf $gtf_name \\
        --output sample_annotation.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        drop_sample_annot: v1.0
    END_VERSIONS
    """

    stub:
    """
    touch sample_annotation.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        drop_sample_annot: v1.0
    END_VERSIONS
    """
}
