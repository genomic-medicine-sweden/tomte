//
// Allele specific variant calling
//

include { BCFTOOLS_VIEW        } from '../../modules/nf-core/bcftools/view/main'
include { BCFTOOLS_INDEX       } from '../../modules/nf-core/bcftools/index/main'
include { GATK4_ASEREADCOUNTER } from '../../modules/nf-core/gatk4/asereadcounter/main'
include { BOOTSTRAPANN         } from '../../modules/local/bootstrapann'
include { TABIX_BGZIPTABIX     } from '../../modules/nf-core/tabix/bgziptabix/main'
include { BCFTOOLS_MERGE       } from '../../modules/nf-core/bcftools/merge/main'
include { RENAME_FILES         } from '../../modules/local/rename_files'
include { TABIX_TABIX          } from '../../modules/nf-core/tabix/tabix/main'


workflow ALLELE_SPECIFIC_CALLING {
    take:
        ch_ind_vcf_tbi // channel: [mandatory] [ val(meta), [ path(vcf), path(tbi) ] ]
        ch_bam_bai     // channel: [mandatory] [ val(meta), [ path(bam), path(bai) ] ]
        ch_fasta       // channel: [mandatory] [ val(meta), path(fasta) ]
        ch_fai         // channel: [mandatory] [ val(meta), path(fai) ]
        ch_dict        // channel: [mandatory] [ val(meta), path(dict) ]
        ch_intervals   // channel: [mandatory] [ path(intervals) ]
        ch_case_info   // channel: [mandatory] [ val(case_info) ]

    main:
        ch_versions = Channel.empty()

        // Keep only does variants in the vcf that are SNVs and are heterozygote
        BCFTOOLS_VIEW(
            ch_ind_vcf_tbi,
            [],
            [],
            []
        )

        BCFTOOLS_INDEX(
            BCFTOOLS_VIEW.out.vcf
        )

        GATK4_ASEREADCOUNTER(
            ch_bam_bai,
            BCFTOOLS_VIEW.out.vcf.join(BCFTOOLS_INDEX.out.tbi),
            ch_fasta,
            ch_fai,
            ch_dict,
            ch_intervals
        )

        BOOTSTRAPANN(
            ch_ind_vcf_tbi,
            GATK4_ASEREADCOUNTER.out.csv
        )

        TABIX_BGZIPTABIX(BOOTSTRAPANN.out.vcf)

        TABIX_BGZIPTABIX.out.gz_tbi
                    .collect{it[1]}
                    .ifEmpty([])
                    .toList()
                    .set { file_list_vcf }

        TABIX_BGZIPTABIX.out.gz_tbi
                    .collect{it[2]}
                    .ifEmpty([])
                    .toList()
                    .set { file_list_tbi }

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
            ch_fasta,
            ch_fai,
            []
        )

        RENAME_FILES( ch_case_vcf.single)

        BCFTOOLS_MERGE.out.merged_variants
            .mix( RENAME_FILES.out.output )
            .set { ch_vcf }

        TABIX_TABIX( ch_vcf )
        ch_tbi = TABIX_TABIX.out.tbi

        ch_versions = ch_versions.mix(BCFTOOLS_VIEW.out.versions.first())
        ch_versions = ch_versions.mix(BCFTOOLS_INDEX.out.versions.first())
        ch_versions = ch_versions.mix(GATK4_ASEREADCOUNTER.out.versions.first())
        ch_versions = ch_versions.mix(BOOTSTRAPANN.out.versions.first())
        ch_versions = ch_versions.mix(TABIX_BGZIPTABIX.out.versions.first())
        ch_versions = ch_versions.mix( BCFTOOLS_MERGE.out.versions.first() )
        ch_versions = ch_versions.mix( RENAME_FILES.out.versions.first() )
        ch_versions = ch_versions.mix( TABIX_TABIX.out.versions.first() )

    emit:
        vcf      = ch_vcf      // channel: [ val(meta), [ path(vcf) ] ]
        tbi      = ch_tbi      // channel: [ val(meta), [ path(tbi) ] ]
        versions = ch_versions // channel: [ path(versions.yml) ]
}
