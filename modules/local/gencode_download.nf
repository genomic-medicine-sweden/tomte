process GENECODE_DOWNLOAD {
    tag "ensembl"
    label 'process_low'

    conda "bioconda::gnu-wget=1.18"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gnu-wget:1.18--h5bf99c6_5' :
        'quay.io/biocontainers/gnu-wget:1.18--h5bf99c6_5' }"

    input:
    val ensembl_version
    val gencode_version
    val meta

    output:
    tuple val(meta), path("*.fa")  , emit: fasta
    tuple val(meta), path("*.gtf") , emit: gtf
    path "versions.yml"            , emit: versions
    

    script:
    """
    wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_${gencode_version}/${genome}.primary_assembly.genome.fa.gz
    wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_${gencode_version}/gencode.v${gencode_version}.primary_assembly.annotation.gtf.gz

    gunzip ${genome}.primary_assembly.genome.fa.gz
    gunzip gencode.v${gencode_version}.primary_assembly.annotation.gtf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        wget: \$(echo wget -V 2>&1 | grep "GNU Wget" | cut -d" " -f3 > versions.yml)
    END_VERSIONS
    """

    stub:
    """
    touch ${genome}.primary_assembly.genome.fa
    touch gencode.v${gencode_version}.primary_assembly.annotation.gtf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        wget: \$(echo wget -V 2>&1 | grep "GNU Wget" | cut -d" " -f3 > versions.yml)
    END_VERSIONS
    """

}