//
// Allele specific variant calling
//

include { BCFTOOLS_VIEW        } from '../../modules/nf-core/bcftools/view/main' 
include { BCFTOOLS_INDEX       } from '../../modules/nf-core/bcftools/index/main'
include { GATK4_ASEREADCOUNTER } from '../../modules/nf-core/gatk4/asereadcounter/main'
include { BOOTSTRAPANN         } from '../../modules/local/bootstrapann'


workflow ALLELE_SPECIFIC_CALLING {
    take:
        vcf_tbi
        bam_bai
        fasta
        fai
        dict
        intervals

    main:
        ch_versions = Channel.empty()

        // Keep only does variants in the vcf that are SNVs and are heterozygote
        BCFTOOLS_VIEW(
            vcf_tbi,
            [],
            [],
            []
        )
        ch_versions = ch_versions.mix(BCFTOOLS_VIEW.out.versions.first())

        BCFTOOLS_INDEX(
            BCFTOOLS_VIEW.out.vcf
        )
        ch_versions = ch_versions.mix(BCFTOOLS_INDEX.out.versions.first())
        
        dict_no_meta = dict.map{ meta, it -> [it] }.collect()
        GATK4_ASEREADCOUNTER(
            bam_bai,
            BCFTOOLS_VIEW.out.vcf.join(BCFTOOLS_INDEX.out.tbi), 
            fasta,
            fai,
            dict_no_meta,
            intervals
        )
        ch_versions = ch_versions.mix(GATK4_ASEREADCOUNTER.out.versions.first())

        BOOTSTRAPANN(
            vcf_tbi,
            GATK4_ASEREADCOUNTER.out.csv
        )
        ch_versions = ch_versions.mix(BOOTSTRAPANN.out.versions.first())


    emit:
        vcf  = BOOTSTRAPANN.out.vcf
        versions = ch_versions // channel: [ path(versions.yml) ]
}
