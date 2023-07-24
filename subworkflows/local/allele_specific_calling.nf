//
// Allele specific variant calling
//

include { BCFTOOLS_VIEW        } from '../../modules/nf-core/bcftools/view/main'
include { BCFTOOLS_INDEX       } from '../../modules/nf-core/bcftools/index/main'
include { GATK4_ASEREADCOUNTER } from '../../modules/nf-core/gatk4/asereadcounter/main'
include { BOOTSTRAPANN         } from '../../modules/local/bootstrapann'


workflow ALLELE_SPECIFIC_CALLING {
    take:
        vcf_tbi   // channel: [ val(meta), path(vcf), path(tbi) ]
        bam_bai   // channel: [ val(meta), path(bam), path(bai) ]
        fasta     // channel: [ path(fasta) ]
        fai       // channel: [ path(fai) ]
        dict      // channel: [ val(meta), [ path(dict) ] ]
        intervals // channel: [ path(interval_list) ]

    main:
        ch_versions = Channel.empty()

        // Keep only does variants in the vcf that are SNVs and are heterozygote
        BCFTOOLS_VIEW(
            vcf_tbi,
            [],
            [],
            []
        )

        BCFTOOLS_INDEX(
            BCFTOOLS_VIEW.out.vcf
        )

        dict_no_meta = dict.map{ meta, it -> [it] }.collect()
        GATK4_ASEREADCOUNTER(
            bam_bai,
            BCFTOOLS_VIEW.out.vcf.join(BCFTOOLS_INDEX.out.tbi),
            fasta,
            fai,
            dict_no_meta,
            intervals
        )

        BOOTSTRAPANN(
            vcf_tbi,
            GATK4_ASEREADCOUNTER.out.csv
        )

        ch_versions = ch_versions.mix(BCFTOOLS_VIEW.out.versions.first())
        ch_versions = ch_versions.mix(BCFTOOLS_INDEX.out.versions.first())
        ch_versions = ch_versions.mix(GATK4_ASEREADCOUNTER.out.versions.first())
        ch_versions = ch_versions.mix(BOOTSTRAPANN.out.versions.first())

    emit:
        vcf  = BOOTSTRAPANN.out.vcf // channel: [ val(meta), path(ann_vcf)]
        versions = ch_versions      // channel: [ path(versions.yml) ]
}
