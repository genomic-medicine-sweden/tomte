//
// Call variants
//

// Modules
include { BCFTOOLS_MPILEUP   } from '../../modules/nf-core/bcftools/mpileup/main'
include { BCFTOOLS_MERGE     } from '../../modules/nf-core/bcftools/merge/main'
include { TABIX_TABIX        } from '../../modules/nf-core/tabix/tabix/main'

// Subworkflows
include { CALL_VARIANTS_GATK } from './call_variants_gatk.nf'

workflow CALL_VARIANTS {
    take:
        ch_bam_bai     // channel (mandatory): [ val(meta), [ path(bam), path(bai) ] ]
        ch_fasta       // channel (mandatory): [ path(fasta) ]
        ch_fai         // channel (mandatory): [ path(fai) ]
        ch_dict        // channel (mandatory): [ path(dict) ]
        variant_caller // string (mandatory)
        ch_case_info   // string (mandatory)

    main:

        ch_versions = Channel.empty()
        ch_vcf = Channel.empty()
        ch_tbi = Channel.empty()
        ch_stats = Channel.empty()

        switch (variant_caller) {

            case 'gatk':

                CALL_VARIANTS_GATK(
                    ch_bam_bai,
                    ch_fasta,
                    ch_fai,
                    ch_dict,
                )

                CALL_VARIANTS_GATK.out.vcf
                    .collect{it[1]}
                    .ifEmpty([])
                    .toList()
                    .set { file_list_vcf }

                CALL_VARIANTS_GATK.out.tbi
                    .collect{it[1]}
                    .ifEmpty([])
                    .toList()
                    .set { file_list_tbi }

                ch_stats = ch_stats.mix(CALL_VARIANTS_GATK.out.stats)
                ch_versions = ch_versions.mix(CALL_VARIANTS_GATK.out.versions.first())

                break

            case 'bcftools':

                BCFTOOLS_MPILEUP(
                    ch_bam_bai.map{ meta, bam, bai -> [ meta, bam, [] ]},
                    ch_fasta.map { it -> [[:], it] },
                    false
                )

                BCFTOOLS_MPILEUP.out.vcf
                    .collect{it[1]}
                    .ifEmpty([])
                    .toList()
                    .set { file_list_vcf }

                BCFTOOLS_MPILEUP.out.tbi
                    .collect{it[1]}
                    .ifEmpty([])
                    .toList()
                    .set { file_list_tbi }

                ch_stats = ch_stats.mix(BCFTOOLS_MPILEUP.out.stats)
                ch_versions = ch_versions.mix(BCFTOOLS_MPILEUP.out.versions.first())

                break

            default:
                exit 1, "Unknown variantcaller: ${variant_caller}"
        }

        ch_case_info
            .combine(file_list_vcf)
            .combine(file_list_tbi)
            .set { ch_vcf_tbi }

        ch_vcf_tbi.branch {
            meta, vcf, tbi ->
                single: vcf.size() == 1
                    return [meta, vcf]
                multiple: vcf.size() > 1
                    return [meta, vcf, tbi]
            }.set { ch_case_vcf }

        BCFTOOLS_MERGE( ch_case_vcf.multiple,
            ch_fasta.map { it -> [[:], it] },
            ch_fai.map { it -> [[:], it] },
            []
        )

        ch_vcf =  BCFTOOLS_MERGE.out.merged_variants.mix(ch_case_vcf.single)
        TABIX_TABIX(ch_vcf)
        ch_tbi=TABIX_TABIX.out.tbi

    emit:
        vcf      = ch_vcf              // channel: [ val(meta), path(vcf) ]
        tbi      = ch_tbi              // channel: [ val(meta), path(tbi) ]
        vcf_tbi  = ch_vcf.join(ch_tbi) // channel: [ val(meta), path(vcf), path(tbi) ]
        stats    = ch_stats            // channel: [ val(meta), path(stats) ]
        versions = ch_versions         // channel: [ path(versions.yml) ]
}
