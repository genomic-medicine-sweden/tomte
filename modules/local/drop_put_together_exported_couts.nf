process DROP_PUT_TOGETHER_EXPORTED_COUNTS {
    tag "DROP_put_together_exported_couts"
    label 'process_low'

    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        exit 1, "Local DROP module does not support Conda. Please use Docker / Singularity / Podman instead."
    }

    container "docker.io/clinicalgenomics/drop:1.4.0"

    input:
    path(exported_counts_ae)
    path(exported_counts_as)
    tuple val(meta), path(gtf)

    output:
    path('exported_counts'), emit: drop_exported_counts
    path "versions.yml"    , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def ae_run = exported_counts_ae ? true : false
    def as_run = exported_counts_as ? true : false
    def gtf_no_extension = gtf.baseName

    """

    mkdir -p exported_counts
    if [[ "$ae_run" == "true" ]];then
        cp ${exported_counts_ae}/* exported_counts/.
    fi

    if [[ "$as_run" == "true" ]];then
        cp ${exported_counts_as}/* exported_counts/.
    fi

    mv exported_counts/sample_annotation.tsv .

    $baseDir/bin/drop_sample_annot_exported_counts.py \\
        --sample_annot "sample_annotation.tsv" \\
        --ae_run $ae_run \\
        --as_run $as_run \\
        --gtf $gtf_no_extension

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        drop_sample_annotation_exported_counts: \$(\$baseDir/bin/drop_sample_annotation_exported_counts.py --version )
    END_VERSIONS

    """

    stub:
    """
    mkdir -p exported_counts

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        drop_sample_annotation_exported_counts: \$(\$baseDir/bin/drop_sample_annotation_exported_counts.py --version )
    END_VERSIONS
    """
}
