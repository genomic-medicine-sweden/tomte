process DROP_CONFIG_RUN_MAE {
    tag "DROP_CONFIG_RUN_MAE"
    label 'process_drop'

    input:
    tuple val(meta), path(fasta), path(fai)
    tuple val(meta2), path(gtf)
    path sample_annotation
    val(genome)
    tuple path(bam), path(bai)
    tuple path(vcf), path(vcf_tbi)
    tuple path(ref_vcf), path(ref_vcf_tbi)

    output:
    path('config.yaml')             , emit: config_drop
    path('output')                  , emit: drop_ae_out
    path('exported_counts_ae')      , emit: gene_counts_ae, optional: true
    path "versions.yml"             , emit: versions

    when:
    // Exit if running this module with -profile conda / -profile mamba
    workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() == 0 ||
        { log.error("Local DROP module does not support Conda. Please use Docker / Singularity / Podman instead."); return false }

    script:
    def args = task.ext.args ?: ''
    def genome_assembly = "${genome}".contains("h37") ? "hg19" : "${genome}"

    """
    TMPDIR=\$PWD
    HOME=\$PWD

    drop init

    drop_config.py \\
        --genome_fasta ${fasta} \\
        --gtf ${gtf} \\
        --drop_module MAE \\
        --genome_assembly $genome_assembly \\
        --ref_vcf $ref_vcf \\
        --output config.yaml

    snakemake mae --cores ${task.cpus} --rerun-triggers mtime $args

    cp output/processed_results/mae/mae/MAE_results_* .

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
    touch MAE_results_genome.tsv.gz
    touch MAE_results_genome_rare.tsv
    mkdir output

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        drop_config: \$(drop_config.py --version)
        drop: \$(echo \$(drop --version) |  sed -n 's/drop, version //p')
    END_VERSIONS
    """
}
