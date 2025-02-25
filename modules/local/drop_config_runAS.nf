process DROP_CONFIG_RUN_AS {
    tag "DROP_CONFIG_RUN_AS"
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
    val(drop_group_samples_as)
    val(drop_group_samples_ae)
    val(drop_padjcutoff_as)
    val(skip_export_counts_drop)

    output:
    path('config.yaml')             , emit: config_drop
    path('output')                  , emit: drop_as_out
    path('FRASER_results_fraser--*'), emit: drop_as_tsv
    path('gene_name_mapping*')      , emit: drop_gene_name
    path('exported_counts_as')      , emit: gene_counts_as, optional: true
    path "versions.yml"             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def genome_assembly = "${genome}".contains("h37") ? "hg19" : "${genome}"
    def drop_group = "${drop_group_samples_as}".replace(" ","")
    def drop_other_group_samples = "${drop_group_samples_ae}".replace(" ","")
    """
    TMPDIR=\$PWD
    HOME=\$PWD


    drop init

    drop_config.py \\
        --genome_fasta ${fasta} \\
        --gtf ${gtf}\\
        --drop_module AS \\
        --genome_assembly $genome_assembly \\
        --drop_group_samples $drop_group \\
        --drop_other_group_samples $drop_other_group_samples \\
        --padjcutoff ${drop_padjcutoff_as} \\
        --output config.yaml

    snakemake aberrantSplicing --cores ${task.cpus} --rerun-triggers mtime $args

    if [[ $skip_export_counts_drop == false ]]; then
        snakemake exportCounts --cores 1
        mkdir -p exported_counts_as
        cp sample_annotation.tsv exported_counts_as/.
        cp output/processed_results/exported_counts/*/*.gz exported_counts_as/.
    fi

    cp output/html/AberrantSplicing/FRASER_results_fraser--*.tsv .
    cp output/processed_data/preprocess/*/gene_name_mapping_*.tsv .

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        drop_config: \$(drop_config.py --version)
        drop: \$(echo \$(drop --version) |  sed -n 's/drop, version //p')
    END_VERSIONS
    """

    stub:
    """
    touch config.yaml
    touch FRASER_results_fraser--.tsv
    touch gene_name_mapping_.tsv
    mkdir output
    if [[ $skip_export_counts_drop == false ]]; then
        mkdir -p exported_counts_as
        touch exported_counts_as/k_j_counts.tsv.gz
        touch exported_counts_as/k_theta_counts.tsv.gz
        touch exported_counts_as/n_psi3_counts.tsv.gz
        touch exported_counts_as/n_psi5_counts.tsv.gz
        touch exported_counts_as/n_theta_counts.tsv.gz
        touch exported_counts_as/sample_annotation.tsv
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        drop_config: \$(drop_config.py --version)
        drop: \$(echo \$(drop --version) |  sed -n 's/drop, version //p')
    END_VERSIONS
    """
}
