process GENERATE_CONFIG_RUN_AE_DROP {
    tag "DROP_AE"
    label 'process_high'

    //conda "bioconda::drop=1.3.3"
    //container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    //    'https://depot.galaxyproject.org/singularity/drop:1.3.3--pyhdfd78af_0' :
    //    'biocontainers/drop:1.3.3--pyhdfd78af_0' }"

    //containerOptions {
    //    (workflow.containerEngine == 'singularity') ?
    //        "--writable" :
    //        "--privileged"
    //    }

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

    $baseDir/bin/generate_drop_config.py \\
        --genome_fasta $fasta \\
        --gtf $gtf \\
        --output config.yaml

    snakemake aberrantExpression --cores ${task.cpus} --rerun-triggers mtime

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        generate_drop_config: v1.0
        drop: v\$(echo \$(drop --version) |  sed -n 's/drop, version //p')
    END_VERSIONS
    """

    stub:
    """
    touch config.yaml
    mkdir output

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        generate_drop_config: v1.0
        drop: v\$(echo \$(drop --version) |  sed -n 's/drop, version //p')
    END_VERSIONS
    """
}
