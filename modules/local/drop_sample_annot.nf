process DROP_SAMPLE_ANNOT {
    tag "DROP_sample_annot"
    label 'process_low'

    conda "bioconda::drop=1.3.3"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/drop:1.3.3--pyhdfd78af_0' :
        'biocontainers/drop:1.3.3--pyhdfd78af_0' }"

    input:
    path(bam)
    val(samples)
    path(processed_gene_counts)
    path(ref_annot)
    path(gtf)

    output:
    path('sample_annotation.tsv'), emit: drop_annot
    path "versions.yml"   , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def strandedness = samples ? "--strandedness ${samples.strandedness}" : ""
    def input_samples = samples ? "${samples.id}" : ""
    def single_end = samples ? "${samples.single_end}" : ""

    """
    $baseDir/bin/drop_sample_annot.py \\
        --bam ${bam} \\
        --sample $input_samples \\
        $strandedness \\
        --single_end $single_end \\
        --gtf $gtf \\
        --count_file ${processed_gene_counts} \\
        --ref_annot ${ref_annot} \\
        --output sample_annotation.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        drop_sample_annot: v1.0
    END_VERSIONS
    """

    stub:
    """
    touch sample_annotation.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        drop_sample_annot: v1.0
    END_VERSIONS
    """
}
