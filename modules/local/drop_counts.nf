process DROP_COUNTS {
    tag "DROP_counts"
    label 'process_low'

    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        exit 1, "Local DROP module does not support Conda. Please use Docker / Singularity / Podman instead."
    }

    container "docker.io/clinicalgenomics/drop:1.3.3"

    input:
    path(counts)
    val(samples)
    path(gtf)
    path(reference_count_file)

    output:
    path('processed_geneCounts.tsv.gz'), emit: processed_gene_counts
    path "versions.yml"             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def ids = "${samples.id}".replace("[","").replace("]","").replace(",","")
    def strandedness = "${samples.strandedness}".replace("[","").replace("]","").replace(",","")
    """
    $baseDir/bin/drop_counts.py \\
        --star ${counts} \\
        --samples $ids \\
        --strandedness $strandedness \\
        --ref_count_file ${reference_count_file} \\
        --output processed_geneCounts.tsv.gz \\
        --gtf ${gtf} \\

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        drop_counts: \$(\$baseDir/bin/drop_counts.py --version )
    END_VERSIONS
    """

    stub:
    """
    touch processed_geneCounts.tsv.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        drop_counts: \$(\$baseDir/bin/drop_counts.py --version )
    END_VERSIONS
    """
}
