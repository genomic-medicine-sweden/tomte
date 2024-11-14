//
// Call variants
//

// Modules
include { BCFTOOLS_MPILEUP                     } from '../../modules/nf-core/bcftools/mpileup/main'
include { BCFTOOLS_NORM as SPLIT_MULTIALLELICS } from '../../modules/nf-core/bcftools/norm/main'
include { BCFTOOLS_NORM as REMOVE_DUPLICATES   } from '../../modules/nf-core/bcftools/norm/main'
include { ADD_FOUND_IN_TAG                     } from '../../modules/local/add_found_in_tag/main'

// Subworkflows
include { CALL_VARIANTS_GATK } from './call_variants_gatk.nf'

workflow CALL_VARIANTS {
    take:
        ch_bam_bai         // channel:   [mandatory] [ val(meta), [ path(bam), path(bai) ] ]
        ch_fasta           // channel:   [mandatory] [ val(meta), path(fasta) ]
        ch_fai             // channel:   [mandatory] [ val(meta),  path(fai) ]
        ch_dict            // channel:   [mandatory] [ val(meta), path(dict) ]
        variant_caller     // parameter: [mandatory] default: 'bcftools'

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

                ch_vcf = ch_vcf.mix(CALL_VARIANTS_GATK.out.vcf)
                ch_tbi = ch_tbi.mix(CALL_VARIANTS_GATK.out.tbi)
                ch_stats = ch_stats.mix(CALL_VARIANTS_GATK.out.stats)
                ch_versions = ch_versions.mix(CALL_VARIANTS_GATK.out.versions.first())

                break

            case 'bcftools':

                BCFTOOLS_MPILEUP(
                    ch_bam_bai.map{ meta, bam, bai -> [ meta, bam, [] ]},
                    ch_fasta,
                    false
                )

                ch_vcf = ch_vcf.mix(BCFTOOLS_MPILEUP.out.vcf)
                ch_tbi = ch_tbi.mix(BCFTOOLS_MPILEUP.out.tbi)
                ch_stats = ch_stats.mix(BCFTOOLS_MPILEUP.out.stats)

                ch_versions = ch_versions.mix(BCFTOOLS_MPILEUP.out.versions.first())

                break

            default:
                exit 1, "Unknown variantcaller: ${variant_caller}"
        }

        ch_in_split_multi = ch_vcf.join(ch_tbi)
        SPLIT_MULTIALLELICS(ch_in_split_multi, ch_fasta)

        ch_remove_dup_in = SPLIT_MULTIALLELICS.out.vcf.join(SPLIT_MULTIALLELICS.out.tbi)
        REMOVE_DUPLICATES(ch_remove_dup_in, ch_fasta)

        ADD_FOUND_IN_TAG(REMOVE_DUPLICATES.out.vcf.join(REMOVE_DUPLICATES.out.tbi), variant_caller)

        ch_vcf_tbi = ADD_FOUND_IN_TAG.out.vcf.join(ADD_FOUND_IN_TAG.out.tbi)

        ch_versions = ch_versions.mix( SPLIT_MULTIALLELICS.out.versions.first() )
        ch_versions = ch_versions.mix( REMOVE_DUPLICATES.out.versions.first() )
        ch_versions = ch_versions.mix( ADD_FOUND_IN_TAG.out.versions.first() )

    emit:
        vcf      = ADD_FOUND_IN_TAG.out.vcf // channel: [ val(meta), path(vcf) ]
        tbi      = ADD_FOUND_IN_TAG.out.tbi // channel: [ val(meta), path(tbi) ]
        vcf_tbi  = ch_vcf_tbi               // channel: [ val(meta), path(vcf), path(tbi) ]
        stats    = ch_stats                 // channel: [ val(meta), path(stats) ]
        versions = ch_versions              // channel: [ path(versions.yml) ]
}
