//
// Annotating SNVs
//

include { ENSEMBLVEP               } from '../../modules/local/ensemblvep'
include { TABIX_TABIX as TABIX_VEP } from '../../modules/nf-core/tabix/tabix/main'


workflow ANNOTATE_SNV {
    take:
    vcf
    val_vep_genome
    val_vep_cache_version
    ch_vep_cache
    ch_fasta

    main:
        ch_versions = Channel.empty()

        // Annotate with VEP
        ENSEMBLVEP(
            vcf,
            val_vep_genome,
            "homo_sapiens",
            val_vep_cache_version,
            ch_vep_cache,
            ch_fasta,
            []
        )
        ch_versions = ch_versions.mix(ENSEMBLVEP.out.versions.first())

        TABIX_VEP(ENSEMBLVEP.out.vcf_gz)
        ch_versions = ch_versions.mix(TABIX_VEP.out.versions.first())

    emit:
        json     = ENSEMBLVEP.out.json
        vcf_gz   = ENSEMBLVEP.out.vcf_gz
        tbi_gz   = TABIX_VEP.out.tbi
        tab_gz   = ENSEMBLVEP.out.tab_gz
        json_gz  = ENSEMBLVEP.out.json_gz
        report   = ENSEMBLVEP.out.report
        versions = ch_versions // channel: [ path(versions.yml) ]
}

