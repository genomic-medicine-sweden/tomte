process VEP_GNOMAD_DOWNLOAD {
    tag "vep_gnomad_download"
    label 'process_very_long'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bcftools:1.20--h8b25389_0':
        'biocontainers/bcftools:1.20--h8b25389_0' }"

    input:
    val genome
    val vep_cache_version

    output:
    tuple path("vcf.gz"), path("vcf.gz.tbi"), emit: gnomad_vcf_tbi
    path "versions.yml"  , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def gnomad_version2download = "${genome}".contains("38") ? "4.0": "2.1.1"
    def base_gnomad_path="https://storage.googleapis.com/gcp-public-data--gnomad/release/${gnomad_version2download}/vcf/genomes/gnomad.genomes."

    """
    # Gnomad
    if [[ "${genome}" == *"38"* ]]; then
        for chr in {1..22} X Y; do
            wget ${base_gnomad_path}v${gnomad_version2download}.sites.chr\${chr}.vcf.bgz --timeout=0
            wget ${base_gnomad_path}v${gnomad_version2download}.sites.chr\${chr}.vcf.bgz.tbi
        done
    else
        wget ${base_gnomad_path}r${gnomad_version2download}.sites.vcf.bgz --timeout=0
        wget ${base_gnomad_path}r${gnomad_version2download}.sites.vcf.bgz.tbi
    fi
    if [[ "${genome}" == *"38"* ]]; then
        bcftools concat gnomad.genomes.*${gnomad_version2download}.sites.chr{1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,X,Y}.vcf.bgz | bcftools annotate --output-type z --output gnomad_v${gnomad_version2download}.vcf.gz --include 'FILTER="PASS"' --remove ^INFO/AF,INFO/AF_grpmax,INFO/AF_popmax --threads $task.cpus --write-index=tbi
    else
        bcftools annotate --output-type z --output gnomad_v${gnomad_version2download}.vcf.gz --include 'FILTER="PASS"' --remove ^INFO/AF,INFO/AF_grpmax,INFO/AF_popmax --threads $task.cpus --write-index=tbi gnomad.genomes.*${gnomad_version2download}.sites.vcf.bgz
    fi
    rm gnomad.genomes.*${gnomad_version2download}.*.vcf.bgz*

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        wget: \$(echo wget -V 2>&1 | grep "GNU Wget" | cut -d" " -f3 > versions.yml)
    END_VERSIONS
    """

    stub:
    """
    touch gnomad_v${gnomad_version2download}.vcf.gz
    touch gnomad_v${gnomad_version2download}.vcf.gz.tbi
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        wget: \$(echo wget -V 2>&1 | grep "GNU Wget" | cut -d" " -f3 > versions.yml)
        bcftools: \$(bcftools --version 2>&1 | head -n1 | sed 's/^.*bcftools //; s/ .*\$//')
        tabix: \$(echo \$(tabix -h 2>&1) | sed 's/^.*Version: //; s/ .*\$//')
    END_VERSIONS
    """

}
