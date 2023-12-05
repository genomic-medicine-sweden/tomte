process DROP_FILTER_RESULTS {
    tag "DROP_FILTER_RESULTS"
    label 'process_low'

    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        exit 1, "Local DROP module does not support Conda. Please use Docker / Singularity / Podman instead."
    }

    container "docker.io/clinicalgenomics/drop:1.3.3"

    input:
    val(samples)
    path gene_panel_clinical_filter
    path out_drop_ae_rds_in
    path out_drop_gene_name_in
    path out_drop_as_tsv_in

    output:
    path('OUTRIDER_provided_samples_top_hits.tsv')         , optional: true, emit: ae_out_unfiltered
    path('OUTRIDER_provided_samples_top_hits_filtered.tsv'), optional: true, emit: ae_out_filtered
    path('FRASER_provided_samples_top_hits.tsv')           , optional: true, emit: as_out_unfiltered
    path('FRASER_provided_samples_top_hits_filtered.tsv')  , optional: true, emit: as_out_filtered
    path "versions.yml"                                    , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def ids = "${samples.id}".replace("[","").replace("]","").replace(",","")
    def gene_panel_filter = gene_panel_clinical_filter ? "--gene_panel ${gene_panel_clinical_filter}" : ''
    def drop_ae_rds = out_drop_ae_rds_in ? "--drop_ae_rds ${out_drop_ae_rds_in}" : ''
    def out_drop_gene_name = out_drop_gene_name_in ? "--out_drop_gene_name ${out_drop_gene_name_in}" : ''
    def out_drop_as_tsv = out_drop_as_tsv_in ? "--out_drop_as_tsv ${out_drop_as_tsv_in}" : ''

    """
    $baseDir/bin/drop_filter_results.py \\
        --samples $ids \\
        $gene_panel_filter \\
        $drop_ae_rds \\
        $out_drop_gene_name \\
        $out_drop_as_tsv \\
        --output_file_subfix_as "provided_samples_top_hits_filtered" \\
        --output_file_subfix_ae "provided_samples_top_hits_filtered"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        drop_filter_results: \$(\$baseDir/bin/drop_filter_results --version )
    END_VERSIONS
    """

    stub:
    """
    touch OUTRIDER_provided_samples_top_hits.tsv
    touch OUTRIDER_provided_samples_top_hits_filtered.tsv
    touch FRASER_provided_samples_top_hits.tsv
    touch FRASER_provided_samples_top_hits_filtered.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        drop_filter_results: \$(\$baseDir/bin/drop_filter_results --version )
    END_VERSIONS
    """
}
