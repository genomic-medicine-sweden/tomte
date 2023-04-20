//
// Allele specific variant calling
//

include { BCFTOOLS_VIEW  } from '../../modules/nf-core/bcftools/view/main' 
include { BCFTOOLS_INDEX } from '../../modules/nf-core/bcftools/index/main'


workflow ALLELE_SPECIFIC_CALLING {
    take:
        vcf_tbi

    main:
        ch_versions = Channel.empty()

        // Keep only does variants in the vcf that are SNVs and are heterozygote
        BCFTOOLS_VIEW(vcf_tbi, [], [], [])
        ch_versions = ch_versions.mix(BCFTOOLS_VIEW.out.versions.first())

        BCFTOOLS_INDEX(BCFTOOLS_VIEW.out.vcf)
        ch_versions = ch_versions.mix(BCFTOOLS_INDEX.out.versions.first())

    emit:
        vcf_tbi  = BCFTOOLS_VIEW.out.vcf.join(BCFTOOLS_INDEX.out.tbi)
        versions = ch_versions // channel: [ path(versions.yml) ]
}
