//
// Annotating SNVs
//

include { ENSEMBLVEP               } from '../../modules/local/ensemblvep'
include { TABIX_TABIX as TABIX_VEP } from '../../modules/nf-core/tabix/tabix/main'


workflow ANNOTATE_SNV {
    take:
    vcf                   // channel (mandatory): [ val(meta), path(ann_vcf) ]
    val_vep_genome        // channel (mandatory): [ val(genome) ]
    val_vep_cache_version // channel (mandatory): [ val(cache_version) ]
    ch_vep_cache          // channel (mandatory): [ path(cache) ]
    ch_fasta              // channel (mandatory): [ path(fasta) ]

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

        TABIX_VEP(ENSEMBLVEP.out.vcf_gz)

        ch_versions = ch_versions.mix(ENSEMBLVEP.out.versions.first())
        ch_versions = ch_versions.mix(TABIX_VEP.out.versions.first())

    emit:
        json     = ENSEMBLVEP.out.json    // channel: [ val(meta), path(json) ]
        vcf_gz   = ENSEMBLVEP.out.vcf_gz  // channel: [ val(meta), path(vcf_gz) ]
        tbi_gz   = TABIX_VEP.out.tbi      // channel: [ val(meta), path(vcf_gz_tbi) ]
        tab_gz   = ENSEMBLVEP.out.tab_gz  // channel: [ val(meta), path(tab_gz) ]
        json_gz  = ENSEMBLVEP.out.json_gz // channel: [ val(meta), path(json_gz) ]
        report   = ENSEMBLVEP.out.report  // channel: [ path(summary_html) ]
        versions = ch_versions            // channel: [ path(versions.yml) ]
}

