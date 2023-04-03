//
// Allignment
//

include { CAT_FASTQ            } from '../../modules/nf-core/cat/fastq/main'
include { FASTP                } from '../../modules/nf-core/fastp/main'
include { STAR_ALIGN           } from '../../modules/nf-core/star/align/main'
include { RNA_DOWNSAMPLE       } from '../../modules/local/rna_downsample'
include { RNA_SUBSAMPLE_REGION } from '../../modules/local/rna_downsample'

workflow ALIGNMENT {
    take:
        reads
        star_index
        gtf
        platform
        downsample_bed

    main:
        ch_versions = Channel.empty()

        CAT_FASTQ(reads)
        ch_versions = ch_versions.mix(CAT_FASTQ.out.versions.first())

        FASTP(CAT_FASTQ.out.reads,[],false,false)
        ch_versions = ch_versions.mix(FASTP.out.versions.first())

        STAR_ALIGN(FASTP.out.reads, star_index, gtf, false, 'illumina', false)
        ch_versions = ch_versions.mix(STAR_ALIGN.out.versions.first())

        RNA_DOWNSAMPLE( STAR_ALIGN.out.bam, downsample_bed )
        ch_versions = ch_versions.mix(RNA_DOWNSAMPLE.out.versions.first())

        RNA_SUBSAMPLE_REGION( RNA_DOWNSAMPLE.out.bam_bai )
        ch_versions = ch_versions.mix(RNA_SUBSAMPLE_REGION.out.versions.first())

    emit:
        merged_reads   = CAT_FASTQ.out.reads
        fastp_report   = FASTP.out.json
        bam_ds         = RNA_SUBSAMPLE_REGION.out.bam
        bam_ds_bai     = RNA_SUBSAMPLE_REGION.out.bam_bai
        gene_counts    = STAR_ALIGN.out.tab
        star_log_final = STAR_ALIGN.out.log_final
        versions       = ch_versions
}
