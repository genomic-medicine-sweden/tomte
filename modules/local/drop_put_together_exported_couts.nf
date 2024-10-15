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
    #!/bin/bash
    mkdir -p exported_counts
    if ($ae_run) ;then
        cp ${exported_counts_ae}/* exported_counts/.
    fi

    if ($as_run) ;then
        cp ${exported_counts_as}/* exported_counts/.
    fi

    cd exported_counts/

    awk -F'\t' 'NR==1 {for(i=1; i<=NF; i++) if(\$i=="RNA_BAM_FILE") rna_col=i} NR>1 {\$rna_col="NA"} 1' OFS='\t' sample_annotation.tsv > sampleAnnotation.tsv
    awk -F'\t' 'NR==1 {for(i=1; i<=NF; i++) if(\$i=="SPLICE_COUNTS_DIR") rna_col=i} NR>1 {\$rna_col="exported_counts"} 1' OFS='\t' sampleAnnotation.tsv > tmpfile && mv tmpfile sampleAnnotation.tsv
    awk -v gene_annot="${gtf_no_extension}" -F'\t' 'NR==1 {for (i=1; i<=NF; i++) if (\$i=="GENE_ANNOTATION") rna_col=i} NR>1 {if (\$rna_col == "NA") \$rna_col = gene_annot} 1' OFS='\t' sampleAnnotation.tsv > tmpfile && mv tmpfile sampleAnnotation.tsv
    awk -F'\t' 'NR==1 {for(i=1; i<=NF; i++) if(\$i=="GENE_COUNTS_FILE") rna_col=i} NR>1 {\$rna_col="exported_counts/geneCounts.tsv.gz"} 1' OFS='\t' sampleAnnotation.tsv > tmpfile && mv tmpfile sampleAnnotation.tsv

    rm sample_annotation.tsv

    cd ..

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gawk: \$(awk -Wversion | sed '1!d; s/.*Awk //; s/,.*//')
    END_VERSIONS
    
    """

    stub:
    """
    mkdir -p exported_counts

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gawk: \$(awk -Wversion | sed '1!d; s/.*Awk //; s/,.*//')
    END_VERSIONS
    """
}
