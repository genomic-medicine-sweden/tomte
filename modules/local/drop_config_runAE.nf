process DROP_CONFIG_RUN_AE {
    tag "DROP_CONFIG_RUN_AE"
    label 'process_high'

    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        exit 1, "Local DROP module does not support Conda. Please use Docker / Singularity / Podman instead."
    }

    container "docker.io/clinicalgenomics/drop:1.3.3"

    input:
    tuple val(meta), path(fasta), path(fai)
    path gtf
    path sample_annotation
    path gene_counts
    val(genome)
    val(drop_padjcutoff_ae)
    val(drop_zScoreCutoff)

    output:
    path('config.yaml')             , emit: config_drop
    path('output')                  , emit: drop_ae_out
    path('OUTRIDER_results_all.Rds'), emit: drop_ae_rds
    path('gene_name_mapping*')      , emit: drop_gene_name
    path "versions.yml"             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def genome_assembly = "${genome}".contains("h37") ? "hg19" : "${genome}"
    def zscorecutoff = drop_zScoreCutoff ? "--zscorecutoff ${drop_zScoreCutoff}" : ''

    """
    TMPDIR=\$PWD

    drop init

    $baseDir/bin/drop_config.py \\
        --genome_fasta ${fasta} \\
        --gtf ${gtf} \\
        --drop_module AE \\
        --genome_assembly $genome_assembly \\
        --padjcutoff ${drop_padjcutoff_ae} \\
        $zscorecutoff \\
        --output config.yaml

    snakemake aberrantExpression --cores ${task.cpus} --rerun-triggers mtime

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

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        drop_config: \$(\$baseDir/bin/drop_config.py --version )
        drop: v\$(echo \$(drop --version) |  sed -n 's/drop, version //p')
    END_VERSIONS
    """
}
