process DROP_CONFIG_RUN_AE {
    tag "DROP_CONFIG_RUN_AE"
    label 'process_drop'

    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        exit 1, "Local DROP module does not support Conda. Please use Docker / Singularity / Podman instead."
    }

    container "docker.io/clinicalgenomics/drop:1.4.0"

    input:
    tuple val(meta), path(fasta), path(fai)
    tuple val(meta2), path(gtf)
    path sample_annotation
    tuple path(bam), path(bai)
    path ref_drop_count_file
    path ref_splice_folder
    val(genome)
    val(drop_group_samples_ae)
    val(drop_group_samples_as)
    val(drop_padjcutoff_ae)
    val(drop_zScoreCutoff)
    val(skip_export_counts_drop)

    output:
    path('config.yaml')             , emit: config_drop
    path('output')                  , emit: drop_ae_out
    path('OUTRIDER_results_all.Rds'), emit: drop_ae_rds
    path('gene_name_mapping*')      , emit: drop_gene_name
    path('geneCounts.tsv.gz')       , emit: gene_counts_ae, optional: true
    path "versions.yml"             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def genome_assembly = "${genome}".contains("h37") ? "hg19" : "${genome}"
    def drop_group = "${drop_group_samples_ae}".replace(" ","")
    def drop_other_group_samples = "${drop_group_samples_as}".replace(" ","")
    def zscorecutoff = drop_zScoreCutoff ? "--zscorecutoff ${drop_zScoreCutoff}" : ''

    """
    TMPDIR=\$PWD
    HOME=\$PWD

    drop init

    $baseDir/bin/drop_config.py \\
        --genome_fasta ${fasta} \\
        --gtf ${gtf} \\
        --drop_module AE \\
        --genome_assembly $genome_assembly \\
        --drop_group_samples $drop_group \\
        --drop_other_group_samples $drop_other_group_samples \\
        --padjcutoff ${drop_padjcutoff_ae} \\
        $zscorecutoff \\
        --output config.yaml

    snakemake aberrantExpression --cores ${task.cpus} --rerun-triggers mtime $args

    if [[ !skip_export_counts_drop ]]; then
        snakemake exportCounts --cores 1
        mkdir -p exported_counts
        cp sample_annotation.tsv exported_counts/.
        cp output/processed_results/exported_counts/*/geneCounts.tsv.gz exported_counts/.
    fi

    cp output/processed_results/aberrant_expression/*/outrider/outrider/OUTRIDER_results_all.Rds .
    cp output/processed_data/preprocess/*/gene_name_mapping_*.tsv .

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        drop_config: \$(\$baseDir/bin/drop_config.py --version )
        drop: v\$(echo \$(drop --version) |  sed -n 's/drop, version //p')
    END_VERSIONS
    """

    stub:
    """
    touch config.yaml
    touch OUTRIDER_results_all.Rds
    touch gene_name_mapping_.tsv
    mkdir output
    mkdir exported_counts

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        drop_config: \$(\$baseDir/bin/drop_config.py --version )
        drop: v\$(echo \$(drop --version) |  sed -n 's/drop, version //p')
    END_VERSIONS
    """
}
