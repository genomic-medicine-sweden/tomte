process GENERATE_ANNOTATION_DROP {
    tag "DROP_annot"
    label 'process_low'

    conda "bioconda::pydamage=0.70"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pydamage:0.70--pyhdfd78af_0' :
        'biocontainers/pydamage:0.70--pyhdfd78af_0' }"

    input:
    path(processed_gene_counts)
    path(gtf)
    path(reference_count_file)

    output:
    path('sample_annotation.tsv'), emit: sample_annotation_drop
    path "versions.yml"              , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def ref_counts = reference_count_file ? "--ref_count_file $reference_count_file" : ""
    def gtf_name   = gtf ? gtf.getBaseName() : ""

    """
    generate_drop_sample_annot.py \\
        --count_file $processed_gene_counts \\
        --gtf $gtf_name \\
        $ref_counts \\
        --output sample_annotation.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """

    stub:
    """
    touch sample_annotation.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
