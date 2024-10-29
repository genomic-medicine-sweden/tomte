//
// Allele specific variant calling
//
include { BCFTOOLS_NORM                        } from '../../modules/nf-core/bcftools/norm/main'
include { BCFTOOLS_VIEW                        } from '../../modules/nf-core/bcftools/view/main'
include { GATK4_ASEREADCOUNTER                 } from '../../modules/nf-core/gatk4/asereadcounter/main'
include { BOOTSTRAPANN                         } from '../../modules/local/bootstrapann'
include { TABIX_BGZIPTABIX                     } from '../../modules/nf-core/tabix/bgziptabix/main'
include { BCFTOOLS_MERGE                       } from '../../modules/nf-core/bcftools/merge/main'
include { RENAME_FILES                         } from '../../modules/local/rename_files'
include { TABIX_TABIX                          } from '../../modules/nf-core/tabix/tabix/main'
include { BCFTOOLS_NORM as SPLIT_MULTIALLELICS } from '../../modules/nf-core/bcftools/norm/main'
include { BCFTOOLS_NORM as REMOVE_DUPLICATES   } from '../../modules/nf-core/bcftools/norm/main'


workflow ALLELE_SPECIFIC_CALLING {
    take:
        ch_ind_vcf_tbi     // channel: [mandatory] [ val(meta), [ path(vcf), path(tbi) ] ]
        ch_bam_bai         // channel: [mandatory] [ val(meta), [ path(bam), path(bai) ] ]
        ch_fasta           // channel: [mandatory] [ val(meta), path(fasta) ]
        ch_fai             // channel: [mandatory] [ val(meta), path(fai) ]
        ch_dict            // channel: [mandatory] [ val(meta), path(dict) ]
        ch_case_info       // channel: [mandatory] [ val(case_info) ]

    main:
        ch_versions = Channel.empty()

        // Keep only one variant per position in the vcf
        BCFTOOLS_NORM(
            ch_ind_vcf_tbi,
            ch_fasta
        )

        // Keep only does variants in the vcf that are SNVs and are heterozygote
        BCFTOOLS_VIEW(
            BCFTOOLS_NORM.out.vcf.join(BCFTOOLS_NORM.out.tbi),
            [],
            [],
            []
        )

        ch_vcf_tbi_sample = BCFTOOLS_VIEW.out.vcf.join(BCFTOOLS_VIEW.out.tbi)
        ch_bam_bai_vcf_tbi = ch_bam_bai.join(ch_vcf_tbi_sample)
        GATK4_ASEREADCOUNTER(
            ch_bam_bai_vcf_tbi,
            ch_fasta,
            ch_fai,
            ch_dict,
            []
        )

        BOOTSTRAPANN(
            ch_ind_vcf_tbi.join(GATK4_ASEREADCOUNTER.out.csv),
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

        BCFTOOLS_MERGE.out.vcf
            .mix( RENAME_FILES.out.output )
            .set { ch_vcf_merged }

        TABIX_TABIX( ch_vcf_merged )

        ch_in_split_multi = ch_vcf_merged.join(TABIX_TABIX.out.tbi)
        SPLIT_MULTIALLELICS(ch_in_split_multi, ch_fasta)

        ch_remove_dup_in = SPLIT_MULTIALLELICS.out.vcf.join(SPLIT_MULTIALLELICS.out.tbi)
        REMOVE_DUPLICATES(ch_remove_dup_in, ch_fasta)

        ch_versions = ch_versions.mix( BCFTOOLS_NORM.out.versions.first() )
        ch_versions = ch_versions.mix( BCFTOOLS_VIEW.out.versions.first() )
        ch_versions = ch_versions.mix( GATK4_ASEREADCOUNTER.out.versions.first() )
        ch_versions = ch_versions.mix( BOOTSTRAPANN.out.versions.first() )
        ch_versions = ch_versions.mix( TABIX_BGZIPTABIX.out.versions.first() )
        ch_versions = ch_versions.mix( BCFTOOLS_MERGE.out.versions.first() )
        ch_versions = ch_versions.mix( RENAME_FILES.out.versions.first() )
        ch_versions = ch_versions.mix( TABIX_TABIX.out.versions.first() )
        ch_versions = ch_versions.mix( SPLIT_MULTIALLELICS.out.versions.first() )
        ch_versions = ch_versions.mix( REMOVE_DUPLICATES.out.versions.first() )

    emit:
        vcf      = REMOVE_DUPLICATES.out.vcf // channel: [ val(meta), [ path(vcf) ] ]
        tbi      = REMOVE_DUPLICATES.out.tbi // channel: [ val(meta), [ path(tbi) ] ]
        versions = ch_versions               // channel: [ path(versions.yml) ]
}
