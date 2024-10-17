process DROP_FILTER_RESULTS {
    tag "DROP_FILTER_RESULTS"
    label 'process_low'

    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        exit 1, "Local DROP module does not support Conda. Please use Docker / Singularity / Podman instead."
    }

    container "docker.io/clinicalgenomics/drop:1.4.0"

    input:
    val(case_info)
    path gene_panel_clinical_filter
    path out_drop_ae_rds_in
    path out_drop_gene_name_in
    path out_drop_as_tsv_in

    output:
    path('*outrider_top_hits_research.tsv') , optional: true, emit: ae_out_research
    path('*outrider_top_hits_clinical.tsv') , optional: true, emit: ae_out_clinical
    path('*fraser_top_hits_research.tsv')   , optional: true, emit: as_out_research
    path('*fraser_top_hits_clinical.tsv')   , optional: true, emit: as_out_clinical
    path "versions.yml"                                      , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def ids = "${case_info.probands}".replace("[","").replace("]","").replace(",","")
    def case_id = "${case_info.id}".replace("[","").replace("]","").replace(",","")
    def gene_panel_filter = gene_panel_clinical_filter ? "--gene_panel ${gene_panel_clinical_filter}" : ''
    def drop_ae_rds = out_drop_ae_rds_in ? "--drop_ae_rds ${out_drop_ae_rds_in}" : ''
    def out_drop_gene_name = out_drop_gene_name_in ? "--out_drop_gene_name ${out_drop_gene_name_in}" : ''
    def out_drop_as_tsv = out_drop_as_tsv_in ? "--out_drop_as_tsv ${out_drop_as_tsv_in}" : ''

    """
    drop_filter_results.py \\
        --samples $ids \\
        $gene_panel_filter \\
        $drop_ae_rds \\
        $out_drop_gene_name \\
        $out_drop_as_tsv \\
        --case $case_id \\

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        drop_filter_results: \$(drop_filter_results --version )
    END_VERSIONS
    """

    stub:
    """
    touch outrider_top_hits_research.tsv
    touch outrider_top_hits_clinical.tsv
    touch fraser_top_hits_research.tsv
    touch fraser_top_hits_clinical.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        drop_filter_results: \$(drop_filter_results.py --version )
    END_VERSIONS
    """
}
