process DROP_CONFIG_RUN_AS {
    tag "DROP_CONFIG_RUN_AS"
    label 'process_high'

    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        exit 1, "Local DROP module does not support Conda. Please use Docker / Singularity / Podman instead."
    }

    container "docker.io/clinicalgenomics/drop:1.3.3"

    input:
    tuple val(meta), path(fasta), path(fai)
    tuple val(meta2), path(gtf)
    path sample_annotation
    tuple path(bam), path(bai)
    path ref_drop_count_file
    path ref_splice_folder
    val(genome)
    val(drop_group_samples_as)
    val(drop_padjcutoff_as)

    output:
    path('config.yaml')             , emit: config_drop
    path('output')                  , emit: drop_as_out
    path('FRASER_results_fraser--*'), emit: drop_as_tsv
    path('gene_name_mapping*')      , emit: drop_gene_name
    path "versions.yml"             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def genome_assembly = "${genome}".contains("h37") ? "hg19" : "${genome}"
    def drop_group = "${drop_group_samples_as}".replace(" ","")
    """
    TMPDIR=\$PWD
    HOME=\$PWD


    drop init

    $baseDir/bin/drop_config.py \\
        --genome_fasta ${fasta} \\
        --gtf ${gtf}\\
        --drop_module AS \\
        --genome_assembly $genome_assembly \\
        --drop_group_samples $drop_group \\
        --padjcutoff ${drop_padjcutoff_as} \\
        --output config.yaml

    snakemake aberrantSplicing --cores ${task.cpus} --rerun-triggers mtime $args

    cp output/html/AberrantSplicing/FRASER_results_fraser--*.tsv .
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
    touch FRASER_results_fraser--.tsv
    touch gene_name_mapping_.tsv
    mkdir output

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        drop_config: \$(\$baseDir/bin/drop_config.py --version )
        drop: v\$(echo \$(drop --version) |  sed -n 's/drop, version //p')
    END_VERSIONS
    """
}
