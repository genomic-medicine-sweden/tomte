//
// Allele specific variant calling
//

include { BCFTOOLS_VIEW  } from '../modules/nf-core/bcftools/view/main' 
include { BCFTOOLS_INDEX } from '../modules/nf-core/bcftools/index/main'


workflow ALLELE_SPECIFIC_CALLING {
    take:
        reads
        star_index
        gtf
        platform
        subsample_bed
        seed_frac
        num_reads
        subsample_region_switch
        downsample_switch

    main:
        ch_versions = Channel.empty()

        // Keep only does variants in the vcf that are SNVs and are heterozygote
        CAT_FASTQ(reads)
        ch_versions = ch_versions.mix(CAT_FASTQ.out.versions.first())

    emit:
        merged_reads   = CAT_FASTQ.out.reads
        fastp_report   = FASTP.out.json
        bam            = STAR_ALIGN.out.bam
        bam_bai        = STAR_ALIGN.out.bam.join(SAMTOOLS_INDEX.out.bai)
        bam_ds_bai     = ch_bam_bai_out
        gene_counts    = STAR_ALIGN.out.tab
        star_log_final = STAR_ALIGN.out.log_final
        versions       = ch_versions
}
