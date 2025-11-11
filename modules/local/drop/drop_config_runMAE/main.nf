process DROP_CONFIG_RUN_MAE {
    tag "DROP_CONFIG_RUN_MAE"
    label 'process_drop'

    input:
    tuple val(meta), path(fasta), path(fai)
    tuple val(meta2), path(gtf)
    tuple val(meta), path(dict)
    path sample_annotation
    val(genome)
    tuple path(bam), path(bai)
    tuple path(vcf), path(vcf_tbi)
    tuple path(drop_mae_high_q_vcf), path(drop_mae_high_q_vcf_tbi)
    val(drop_add_af)

    output:
    path('config.yaml')       , emit: config_drop
    path('output')            , emit: drop_mae_out
    path('MAE_results_*')     , emit: drop_mae_tsv
    path('gene_name_mapping*'), emit: drop_gene_name
    path "versions.yml"       , emit: versions

    when:
    // Exit if running this module with -profile conda / -profile mamba
    workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() == 0 ||
        { log.error("Local DROP module does not support Conda. Please use Docker / Singularity / Podman instead."); return false }

    script:
    def args = task.ext.args ?: ''
    def genome_assembly = "${genome}".contains("h37") ? "hg19" : "${genome}"
    def gtf_basename = gtf.getName().replaceAll(/\.(gtf|gff3?|GTF|GFF3?)$/, '')
    def addAF = drop_add_af ? "--drop_add_af" : ''

    """
    TMPDIR=\$PWD
    HOME=\$PWD

    # Fill missing GENE_ANNOTATION entries with the GTF basename (required for MAE analysis)
    awk -v val="${gtf_basename}" 'BEGIN{FS=OFS="\t"}
    NR==1 {for(i=1;i<=NF;i++) if(\$i=="GENE_ANNOTATION") col=i}
    NR>1 && \$col=="NA" {\$col=val}
    {print}
    ' ${sample_annotation} > tmp && mv tmp ${sample_annotation}

    drop init

    drop_config.py \\
        --genome_fasta ${fasta} \\
        --gtf ${gtf} \\
        --drop_module MAE \\
        --genome_assembly $genome_assembly \\
        --drop_mae_high_q_vcf $drop_mae_high_q_vcf \\
        $addAF \\
        --output config.yaml

    snakemake mae --cores ${task.cpus} --rerun-triggers mtime $args

    cp output/processed_results/mae/mae/MAE_results_*.tsv .
    rm MAE_results_*_rare.tsv
    cp output/processed_data/*/gene_name_mapping_*.tsv .

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        drop_config: \$(drop_config.py --version)
        drop: \$(echo \$(drop --version) |  sed -n 's/drop, version //p')
    END_VERSIONS
    """

    stub:
    """
    touch config.yaml
    touch MAE_results_genome.tsv
    touch gene_name_mapping_.tsv
    mkdir output

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        drop_config: \$(drop_config.py --version)
        drop: \$(echo \$(drop --version) |  sed -n 's/drop, version //p')
    END_VERSIONS
    """
}
