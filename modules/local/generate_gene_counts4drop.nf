process GENERATE_COUNTS_DROP {
    tag "DROP_counts"
    label 'process_low'

    conda "bioconda::pydamage=0.70"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pydamage:0.70--pyhdfd78af_0' :
        'biocontainers/pydamage:0.70--pyhdfd78af_0' }"

    input:
    path(counts)
    val(samples)
    path(gtf)
    path(reference_count_file)

    output:
    path('processed_geneCounts.tsv'), emit: processed_gene_counts
    path "versions.yml"             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def ref_counts = reference_count_file ? "--ref_count_file $reference_count_file" : ""
    def strandedness = samples ? "--strandedness ${samples.strandedness}" : ""
    def input_samples = samples ? "${samples.id}" : ""
    """
    $baseDir/bin/generate_gene_counts.py \\
        --star ${counts} \\
        --sample $input_samples \\
        $strandedness \\
        $ref_counts \\
        --output processed_geneCounts.tsv \\
        --gtf $gtf \\

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """

    stub:
    """
    touch processed_geneCounts.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
