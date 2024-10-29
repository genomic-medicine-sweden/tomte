//
// Call variants
//

// Modules
include { BCFTOOLS_MPILEUP                     } from '../../modules/nf-core/bcftools/mpileup/main'
include { BCFTOOLS_NORM as SPLIT_MULTIALLELICS } from '../../modules/nf-core/bcftools/norm/main'
include { BCFTOOLS_NORM as REMOVE_DUPLICATES   } from '../../modules/nf-core/bcftools/norm/main'
include { ADD_VARCALLER_TO_BED                 } from '../../modules/local/add_varcallername_to_bed'
include { BCFTOOLS_ANNOTATE                    } from '../../modules/nf-core/bcftools/annotate/main'

// Subworkflows
include { CALL_VARIANTS_GATK } from './call_variants_gatk.nf'

workflow CALL_VARIANTS {
    take:
        ch_bam_bai         // channel:   [mandatory] [ val(meta), [ path(bam), path(bai) ] ]
        ch_fasta           // channel:   [mandatory] [ val(meta), path(fasta) ]
        ch_fai             // channel:   [mandatory] [ val(meta),  path(fai) ]
        ch_dict            // channel:   [mandatory] [ val(meta), path(dict) ]
        variant_caller     // parameter: [mandatory] default: 'bcftools'
        ch_foundin_header  // channel:   [mandatory] [ path(header) ]
        ch_genome_chrsizes // channel:   [mandatory] [ path(chrsizes) ]

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

        ch_genome_chrsizes.flatten().map { chromsizes ->
            return [[id:variant_caller], chromsizes]
            }
            .set { ch_varcallerinfo }

        ADD_VARCALLER_TO_BED(ch_varcallerinfo)
        ch_bed_in_annot = ADD_VARCALLER_TO_BED.out.gz.collect()
        ch_tbi_in_annot = ADD_VARCALLER_TO_BED.out.tbi.collect()

        bed_with_tbi = ch_bed_in_annot.join(ch_tbi_in_annot)
        in_bcftools_annot = REMOVE_DUPLICATES.out.vcf.join(REMOVE_DUPLICATES.out.tbi).join(bed_with_tbi)

        BCFTOOLS_ANNOTATE(
            in_bcftools_annot,
            ch_foundin_header
        )
        ch_vcf_tbi = BCFTOOLS_ANNOTATE.out.vcf.join(BCFTOOLS_ANNOTATE.out.tbi)

        ch_versions = ch_versions.mix( BCFTOOLS_MPILEUP.out.versions.first() )
        ch_versions = ch_versions.mix( SPLIT_MULTIALLELICS.out.versions.first() )
        ch_versions = ch_versions.mix( REMOVE_DUPLICATES.out.versions.first() )
        ch_versions = ch_versions.mix( ADD_VARCALLER_TO_BED.out.versions.first() )
        ch_versions = ch_versions.mix( BCFTOOLS_ANNOTATE.out.versions.first() )

    emit:
        vcf      = BCFTOOLS_ANNOTATE.out.vcf // channel: [ val(meta), path(vcf) ]
        tbi      = BCFTOOLS_ANNOTATE.out.tbi // channel: [ val(meta), path(tbi) ]
        vcf_tbi  = ch_vcf_tbi                // channel: [ val(meta), path(vcf), path(tbi) ]
        stats    = ch_stats                  // channel: [ val(meta), path(stats) ]
        versions = ch_versions               // channel: [ path(versions.yml) ]
}
