//
// Annotating SNVs
//

include { ENSEMBLVEP_VEP           } from '../../modules/nf-core/ensemblvep/vep/main'
include { TABIX_TABIX as TABIX_VEP } from '../../modules/nf-core/tabix/tabix/main'


workflow ANNOTATE_SNV {
    take:
    vcf                   // channel:   [mandatory] [ val(meta), path(vcf), path(tbi) ]
    val_vep_genome        // parameter: [mandatory] 'GRCh37' or 'GRCh38'
    val_vep_cache_version // parameter: [mandatory] default: 110
    ch_vep_cache          // channel:   [mandatory] [ path(cache) ]
    ch_fasta              // channel:   [mandatory] [ val(meta), path(fasta) ]
    ch_vep_extra_files    // channel:   [mandatory] [ path(files) ]

    main:
        ch_versions = Channel.empty()

        // Annotate with VEP
        ENSEMBLVEP_VEP(
            vcf.map{ meta, vcf -> [ meta, vcf, [] ] },
            val_vep_genome,
            "homo_sapiens",
            val_vep_cache_version,
            ch_vep_cache,
            ch_fasta,
            ch_vep_extra_files
        )
        ch_versions = ch_versions.mix(ENSEMBLVEP_VEP.out.versions.first())

        TABIX_VEP(ENSEMBLVEP_VEP.out.vcf)
        ch_versions = ch_versions.mix(TABIX_VEP.out.versions.first())

    emit:
        vcf_gz   = ENSEMBLVEP_VEP.out.vcf    // channel: [ val(meta), path(vcf.gz) ]
        tbi_gz   = TABIX_VEP.out.tbi         // channel: [ val(meta), path(tbi) ]
        tab_gz   = ENSEMBLVEP_VEP.out.tab    // channel: [ val(meta), path(tab.gz) ]
        json_gz  = ENSEMBLVEP_VEP.out.json   // channel: [ val(meta), path(json.gz) ]
        report   = ENSEMBLVEP_VEP.out.report // channel: [ path(html) ]
        versions = ch_versions               // channel: [ path(versions.yml) ]
}

