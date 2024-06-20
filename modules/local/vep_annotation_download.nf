process GENECODE_DOWNLOAD {
    tag "vep"
    label 'process_low'

    conda "bioconda::gnu-wget=1.18"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gnu-wget:1.18--h5bf99c6_5' :
        'quay.io/biocontainers/gnu-wget:1.18--h5bf99c6_5' }"

    input:
    val vep_cache_version
    val genome
    val gnomad_version2download

    output:
    tuple val(meta), path("vep_cache")  , emit: vep_cache
    path "versions.yml"                 , emit: versions


    script:

    """
    vep_url="ftp://ftp.ensembl.org/pub/release-${vep_cache_version}/variation/indexed_vep_cache/homo_sapiens_merged_vep_${vep_cache_version}_${genome}.tar.gz"
    clinvar_vcf_url="ftp://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_${genome}/weekly/clinvar.vcf.gz"
    clinvar_tbi_url="ftp://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_${genome}/weekly/clinvar.vcf.gz.tbi"
    current_date=/$(date -I)
    
    # Vep cache
    mkdir vep_cache; cd vep_cache
    wget -O homo_sapiens_merged_vep.tar.gz $vep_url && tar xvf homo_sapiens_merged_vep.tar.gz && rm homo_sapiens_merged_vep.tar.gz

    # Vep plugins
    mkdir vep_plugins; cd vep_plugins

    # Clinvar
    wget -O clinvar_${current_date}.vcf.gz $clinvar_vcf_url
    wget -O clinvar_${current_date}.vcf.gz.tbi $clinvar_tbi_url

    # Create file listing all vep plugins to use. Note that file extension must be .csv
    gnomad_vcf="gnomad_v${gnomad_version2download}.vcf.gz"
    echo "vep_files" > vep_files.csv
    echo "vep_plugins/clinvar_${current_date}.vcf.gz" >> vep_files.csv
    echo "vep_plugins/clinvar_${current_date}.vcf.gz.tbi" >> vep_files.csv
    echo "vep_plugins/${gnomad_vcf}" >> vep_files.csv
    echo "vep_plugins/${gnomad_vcf}.tbi" >> vep_files.csv

    gsutil rsync -r gs://gcp-public-data--gnomad/release/${gnomad_version2download}/vcf/genomes/ .
    bcftools concat gnomad.genomes.*${gnomad_version2download}.sites.*{1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,X,Y}.vcf.bgz | bcftools annotate --output-type z --output gnomad_v${gnomad_version2download}.vcf.gz --include 'FILTER="PASS"' --remove ^INFO/AF,INFO/AF_grpmax,INFO/AF_popmax
    tabix -p vcf gnomad_v${gnomad_version2download}.vcf.gz
    rm gnomad.genomes.*${gnomad_version2download}.*.vcf.bgz*


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        wget: \$(echo wget -V 2>&1 | grep "GNU Wget" | cut -d" " -f3 > versions.yml)
    END_VERSIONS
    """

    stub:
    """
    touch vep_cache

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        wget: \$(echo wget -V 2>&1 | grep "GNU Wget" | cut -d" " -f3 > versions.yml)
    END_VERSIONS
    """

}