process VEP_DOWNLOAD {
    tag "vep"
    label 'process_long'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bcftools:1.20--h8b25389_0':
        'biocontainers/bcftools:1.20--h8b25389_0' }"

    input:
    val genome
    val vep_cache_version

    output:
    path("vep_cache")      , emit: vep_cache
    path("vep_files.csv")  , emit: plugin_file
    path "versions.yml"    , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def gnomad_version2download = "${genome}".contains("38") ? "4.0": "2.1.1"
    def current_date  = new java.util.Date().format( 'yyyy-MM-dd')
    def gnomad_vcf="gnomad_v${gnomad_version2download}.vcf.gz"
    def vep_url="ftp://ftp.ensembl.org/pub/release-${vep_cache_version}/variation/indexed_vep_cache/homo_sapiens_merged_vep_${vep_cache_version}_${genome}.tar.gz"
    def clinvar_vcf_url="ftp://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_${genome}/weekly/clinvar.vcf.gz"
    def clinvar_tbi_url="ftp://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_${genome}/weekly/clinvar.vcf.gz.tbi"
    def base_gnomad_path="https://storage.googleapis.com/gcp-public-data--gnomad/release/${gnomad_version2download}/vcf/genomes/gnomad.genomes."

    """
    # Create file listing all vep plugins to use. Note that file extension must be .csv
    echo "vep_files" > vep_files.csv
    echo "vep_cache/vep_plugins/clinvar_${current_date}.vcf.gz" >> vep_files.csv
    echo "vep_cache/vep_plugins/clinvar_${current_date}.vcf.gz.tbi" >> vep_files.csv
    echo "vep_cache/vep_plugins/${gnomad_vcf}" >> vep_files.csv
    echo "vep_cache/vep_plugins/${gnomad_vcf}.tbi" >> vep_files.csv

    # Vep cache
    mkdir vep_cache; cd vep_cache
    wget -O homo_sapiens_merged_vep.tar.gz $vep_url && tar xvf homo_sapiens_merged_vep.tar.gz && rm homo_sapiens_merged_vep.tar.gz

    # Vep plugins
    mkdir vep_plugins; cd vep_plugins

    # Clinvar
    wget -O clinvar_${current_date}.vcf.gz $clinvar_vcf_url
    wget -O clinvar_${current_date}.vcf.gz.tbi $clinvar_tbi_url

    # Gnomad
    for chr in {1..22} X Y; do
        if [[ "${genome}" == *"38"* ]]; then
            wget ${base_gnomad_path}v${gnomad_version2download}.sites.chr\${chr}.vcf.bgz
            wget ${base_gnomad_path}v${gnomad_version2download}.sites.chr\${chr}.vcf.bgz.tbi
        else
            if [[ "\$chr" != "Y" ]]; then
                wget ${base_gnomad_path}r${gnomad_version2download}.sites.\${chr}.vcf.bgz
                wget ${base_gnomad_path}r${gnomad_version2download}.sites.\${chr}.vcf.bgz.tbi
            fi
        fi
    done

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
    touch vep_plugins.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        wget: \$(echo wget -V 2>&1 | grep "GNU Wget" | cut -d" " -f3 > versions.yml)
    END_VERSIONS
    """

}
