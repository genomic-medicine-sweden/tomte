process DROP_SAMPLE_ANNOT {
    tag "DROP_sample_annot"
    label 'process_low'

    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        exit 1, "Local DROP module does not support Conda. Please use Docker / Singularity / Podman instead."
    }

    container "docker.io/clinicalgenomics/drop:1.3.3"

    input:
    path(bam)
    val(samples)
    path(ref_gene_counts)
    path(ref_annot)
    val(drop_group_samples_ae)
    val(drop_group_samples_as)

    output:
    path('sample_annotation.tsv'), emit: drop_annot
    path "versions.yml"   , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def ids = "${samples.id}".replace("[","").replace("]","").replace(",","")
    def strandedness = "${samples.strandedness}".replace("[","").replace("]","").replace(",","")
    def single_end = "${samples.single_end}".replace("[","").replace("]","").replace(",","")
    def drop_group = "${drop_group_samples_ae},${drop_group_samples_as}".replace(" ","").replace("[","").replace("]","")
    """
    $baseDir/bin/drop_sample_annot.py \\
        --bam ${bam} \\
        --samples $ids \\
        --strandedness $strandedness \\
        --single_end $single_end \\
        --ref_count_file ${ref_gene_counts} \\
        --ref_annot ${ref_annot} \\
        --drop_group_sample $drop_group \\
        --output sample_annotation.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        drop_sample_annot: \$(\$baseDir/bin/drop_sample_annot --version )
    END_VERSIONS
    """

    stub:
    """
    touch sample_annotation.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        drop_sample_annot: \$(\$baseDir/bin/drop_sample_annot --version )
    END_VERSIONS
    """
}
