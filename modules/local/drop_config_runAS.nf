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
    path gtf
    path sample_annotation
    tuple val(meta), path(bam), path(bai)
    path ref_splice_folder

    output:
    path('config.yaml'), emit: config_drop
    path('output')     , emit: drop_ae_out
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    TMPDIR=\$PWD

    drop init

    $baseDir/bin/drop_config.py \\
        --genome_fasta $fasta \\
        --gtf $gtf \\
        --drop_module AS \\
        --output config.yaml

    snakemake aberrantSplicing --cores ${task.cpus} --rerun-triggers mtime

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        drop_config: v1.0
        drop: v\$(echo \$(drop --version) |  sed -n 's/drop, version //p')
    END_VERSIONS
    """

    stub:
    """
    touch config.yaml
    mkdir output

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        drop_config: v1.0
        drop: v\$(echo \$(drop --version) |  sed -n 's/drop, version //p')
    END_VERSIONS
    """
}
