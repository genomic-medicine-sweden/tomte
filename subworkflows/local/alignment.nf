//
// Allignment
//

include { CAT_FASTQ            } from '../../modules/nf-core/cat/fastq/main'
include { FASTP                } from '../../modules/nf-core/fastp/main'
include { STAR_ALIGN           } from '../../modules/nf-core/star/align/main'
include { SAMTOOLS_INDEX       } from '../../modules/nf-core/samtools/index/main'
include { RNA_DOWNSAMPLE       } from '../../modules/local/rna_downsample'
include { RNA_SUBSAMPLE_REGION } from '../../modules/local/rna_subsample_region.nf'
include { SALMON_QUANT         } from '../../modules/nf-core/salmon/quant/main'
include { SAMTOOLS_VIEW        } from '../../modules/nf-core/samtools/view/main'

workflow ALIGNMENT {
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
        salmon_index
        ch_genome_fasta

    main:
        ch_versions = Channel.empty()

        ch_fastq = branchFastqToSingleAndMulti(reads)

        CAT_FASTQ(ch_fastq.multiple_fq)
        .reads.mix(ch_fastq.single_fq)
        .set { ch_cat_fastq }

        FASTP(ch_cat_fastq, [], false, false)

        STAR_ALIGN(FASTP.out.reads, star_index, gtf, false, 'illumina', false)

        SAMTOOLS_INDEX( STAR_ALIGN.out.bam )

        ch_bam_bai = Channel.empty()
        ch_bam_bai_out = Channel.empty()

        if (subsample_region_switch) {
            RNA_SUBSAMPLE_REGION( STAR_ALIGN.out.bam, subsample_bed, seed_frac)
            ch_bam_bai = ch_bam_bai.mix(RNA_SUBSAMPLE_REGION.out.bam_bai)
            if (!downsample_switch) {
                ch_bam_bai_out = RNA_SUBSAMPLE_REGION.out.bam_bai
            } else {
                RNA_DOWNSAMPLE( ch_bam_bai, num_reads)
                ch_bam_bai_out = RNA_DOWNSAMPLE.out.bam_bai
            }
        } else {
            ch_bam_bai = ch_bam_bai.mix(STAR_ALIGN.out.bam.join(SAMTOOLS_INDEX.out.bai))
             if (!downsample_switch) {
                ch_bam_bai_out = STAR_ALIGN.out.bam.join(SAMTOOLS_INDEX.out.bai)
            } else {
                RNA_DOWNSAMPLE( ch_bam_bai, num_reads)
                ch_bam_bai_out = RNA_DOWNSAMPLE.out.bam_bai
            }
        }

        SAMTOOLS_VIEW( STAR_ALIGN.out.bam.join(SAMTOOLS_INDEX.out.bai), ch_genome_fasta, [] )


        SALMON_QUANT( FASTP.out.reads, salmon_index, gtf, [], false, 'A')

        ch_versions = ch_versions.mix(CAT_FASTQ.out.versions.first())
        ch_versions = ch_versions.mix(FASTP.out.versions.first())
        ch_versions = ch_versions.mix(STAR_ALIGN.out.versions.first())
        ch_versions = ch_versions.mix(SAMTOOLS_INDEX.out.versions.first())
        ch_versions = ch_versions.mix(RNA_SUBSAMPLE_REGION.out.versions.first())
        ch_versions = ch_versions.mix(SAMTOOLS_VIEW.out.versions.first())
        ch_versions = ch_versions.mix(RNA_DOWNSAMPLE.out.versions.first())
        ch_versions = ch_versions.mix(SALMON_QUANT.out.versions.first())

    emit:
        merged_reads    = CAT_FASTQ.out.reads
        fastp_report    = FASTP.out.json
        bam             = STAR_ALIGN.out.bam
        bam_bai         = ch_bam_bai
        bam_ds_bai      = ch_bam_bai_out
        gene_counts     = STAR_ALIGN.out.read_per_gene_tab
        spl_junc        = STAR_ALIGN.out.spl_junc_tab
        star_log_final  = STAR_ALIGN.out.log_final
        star_wig        = STAR_ALIGN.out.wig
        salmon_result   = SALMON_QUANT.out.results
        salmon_info     = SALMON_QUANT.out.json_info
        versions        = ch_versions
}


// Custom functions

/**
* Branch the read channel into differnt channels,
* depending on whether the sample has multiple fastq files or not.
* The resulting channels gets the original sample id in meta.
*
* @param ch_reads Channel containing meta and fastq reads
* @return Channel containing meta with original id and branched on number of fastq files
*/
def branchFastqToSingleAndMulti(ch_reads) {

    return ch_reads
        .map {
            meta, fastq ->
                original_id = meta.id.split('_T')[0..-2].join('_')
                [ meta + [id: original_id], fastq ]
        }
        .groupTuple()
        .branch {
            meta, fastq ->
                single_fq: fastq.size() == 1
                    return [ meta, fastq.flatten() ]
                multiple_fq: fastq.size() > 1
                    return [ meta, fastq.flatten() ]
        }
}
