//
// Alignment
//

include { CAT_FASTQ            } from '../../../modules/nf-core/cat/fastq'
include { FASTP                } from '../../../modules/nf-core/fastp'
include { STAR_ALIGN           } from '../../../modules/nf-core/star/align'
include { SAMTOOLS_FASTQ       } from '../../../modules/nf-core/samtools/fastq'
include { SAMTOOLS_INDEX       } from '../../../modules/nf-core/samtools/index'
include { RNA_DOWNSAMPLE       } from '../../../modules/local/rna_downsample'
include { RNA_SUBSAMPLE_REGION } from '../../../modules/local/rna_subsample_region'
include { SALMON_QUANT         } from '../../../modules/nf-core/salmon/quant'
include { SAMTOOLS_VIEW        } from '../../../modules/nf-core/samtools/view'

workflow ALIGNMENT {
    take:
    ch_fastq_input_reads   //   channel: [optional]  [ val(meta), [path(reads)] ]
    ch_bam_bai_input_reads //   channel: [optional]  [ val(meta), [path(bam) path(bai)] ]
    star_index             //   channel: [mandatory] [ val(meta), path(star_index) ]
    ch_gtf                 //   channel: [mandatory] [ val(meta), path(gtf) ]
    ch_platform            //   channel: [mandatory] [ val(platform) ]
    subsample_bed          //   channel: [optional]  [ path(subsample_bed) ]
    seed_frac              // parameter: [optional]  default: 0.001
    num_reads              // parameter: [optional]  default: 120000000
    skip_subsample_region  // parameter: [mandatory] default: true
    skip_downsample        // parameter: [mandatory] default: true
    salmon_index           //   channel: [mandatory] [ path(salmon_index) ]
    ch_genome_fasta        //   channel: [mandatory] [ val(meta), path(fasta) ]
    save_mapped_as_cram    // parameter: [mandatory] default: true

    main:
    ch_versions = Channel.empty()

    // Branch inputs based on meta.is_fastq flag
    ch_all_inputs = ch_fastq_input_reads.mix(ch_bam_bai_input_reads)

    ch_inputs_branched = ch_all_inputs.branch { meta, files ->
        fastq: meta.is_fastq == true
            return [meta, files]
        bam: meta.is_fastq == false
            return [meta, files]
    }

    // Process FASTQ inputs
    ch_fastq = branchFastqToSingleAndMulti(ch_inputs_branched.fastq)

    CAT_FASTQ(ch_fastq.multiple_fq)
    ch_cat_fastq = CAT_FASTQ.out.reads.mix(ch_fastq.single_fq)

    FASTP(ch_cat_fastq, [], false, false, false)

    STAR_ALIGN(FASTP.out.reads, star_index, ch_gtf, false, ch_platform, false)

    // Process BAM inputs - extract BAM files and convert to FASTQ
    ch_bam_input_reads = ch_inputs_branched.bam.map { meta, bambai -> [ meta, bambai[0] ] }
    ch_bai_input_reads = ch_inputs_branched.bam.map { meta, bambai -> [ meta, bambai[1] ] }

    SAMTOOLS_FASTQ(ch_bam_input_reads, false)

    // Combine all FASTQ inputs for Salmon
    ch_salmon_input = FASTP.out.reads.mix(SAMTOOLS_FASTQ.out.fastq)

    // Run SALMON_QUANT
    SALMON_QUANT(ch_salmon_input, salmon_index, ch_gtf.map{ _meta, gtf -> gtf }, [], false, 'A')

    // Continue with BAM processing
    ch_bam_from_star = STAR_ALIGN.out.bam_sorted_aligned
    ch_bam_2_process = ch_bam_input_reads.mix(ch_bam_from_star)

    SAMTOOLS_INDEX(ch_bam_from_star)
    ch_bai = ch_bai_input_reads.mix(SAMTOOLS_INDEX.out.bai)

    // Subsampling and downsampling logic
    if (!skip_subsample_region) {
        RNA_SUBSAMPLE_REGION(ch_bam_2_process, subsample_bed, seed_frac)
        ch_bam_bai_not_downsamp = RNA_SUBSAMPLE_REGION.out.bam_bai
        ch_versions = ch_versions.mix(RNA_SUBSAMPLE_REGION.out.versions.first())
        if (skip_downsample) {
            ch_bam_bai_input_drop = RNA_SUBSAMPLE_REGION.out.bam_bai
        } else {
            RNA_DOWNSAMPLE(ch_bam_bai_not_downsamp, num_reads)
            ch_bam_bai_input_drop = RNA_DOWNSAMPLE.out.bam_bai
            ch_versions = ch_versions.mix(RNA_DOWNSAMPLE.out.versions.first())
        }
    } else {
        ch_bam_bai_not_downsamp = ch_bam_2_process.join(ch_bai)
        if (skip_downsample) {
            ch_bam_bai_input_drop = ch_bam_2_process.join(ch_bai)
        } else {
            RNA_DOWNSAMPLE(ch_bam_bai_not_downsamp, num_reads)
            ch_bam_bai_input_drop = RNA_DOWNSAMPLE.out.bam_bai
            ch_versions = ch_versions.mix(RNA_DOWNSAMPLE.out.versions.first())
        }
    }

    if (save_mapped_as_cram) {
        SAMTOOLS_VIEW(ch_bam_2_process.join(ch_bai), ch_genome_fasta, [], 'crai')
        ch_versions = ch_versions.mix(SAMTOOLS_VIEW.out.versions.first())
        ch_cram_crai = SAMTOOLS_VIEW.out.cram.join(SAMTOOLS_VIEW.out.crai)
    } else {
        ch_cram_crai = Channel.empty()
    }

    // Collect versions
    ch_versions = ch_versions.mix(CAT_FASTQ.out.versions.first())
    ch_versions = ch_versions.mix(FASTP.out.versions.first())
    ch_versions = ch_versions.mix(STAR_ALIGN.out.versions.first())
    ch_versions = ch_versions.mix(SAMTOOLS_FASTQ.out.versions.first())
    ch_versions = ch_versions.mix(SAMTOOLS_INDEX.out.versions.first())
    ch_versions = ch_versions.mix(SALMON_QUANT.out.versions.first())

    emit:
    merged_reads    = CAT_FASTQ.out.reads              // channel: [ val(meta), path(fastq) ]
    fastp_report    = FASTP.out.json                   // channel: [ val(meta), path(json) ]
    bam             = ch_bam_2_process                 // channel: [ val(meta), path(bam) ]
    bam_bai         = ch_bam_bai_not_downsamp          // channel: [ val(meta), path(bam), path(bai) ]
    bam_ds_bai      = ch_bam_bai_input_drop            // channel: [ val(meta), path(bam), path(bai) ]
    gene_counts     = STAR_ALIGN.out.read_per_gene_tab // channel: [ val(meta), path(tsv) ]
    spl_junc        = STAR_ALIGN.out.spl_junc_tab      // channel: [ val(meta), path(tsv) ]
    star_log_final  = STAR_ALIGN.out.log_final         // channel: [ val(meta), path(log) ]
    star_wig        = STAR_ALIGN.out.wig               // channel: [ val(meta), path(wig) ]
    salmon_result   = SALMON_QUANT.out.results         // channel: [ val(meta), path(results) ]
    salmon_info     = SALMON_QUANT.out.json_info       // channel: [ val(meta), path(json) ]
    cram_crai       = ch_cram_crai                     // channel: [ val(meta), path(cram), path(crai) ]
    versions        = ch_versions                      // channel: [ path(versions.yml) ]
}


// Custom functions

/**
* Branch the read channel into different channels,
* depending on whether the sample has multiple fastq files or not.
* The resulting channels gets the original sample id in meta.
*
* @param ch_reads Channel containing meta and fastq reads
* @return Channel containing meta with original id and branched on number of fastq files
*/
def branchFastqToSingleAndMulti(ch_reads) {

    return ch_reads
        .map { meta, fastqs ->
            return [ groupKey( meta + [id:meta.sample], meta.fq_pairs ), fastqs ]
        }
        .groupTuple()
        .branch {
            meta, fastqs ->
                single_fq: fastqs.size() == 1
                    return [ meta, fastqs.flatten() ]
                multiple_fq: fastqs.size() > 1
                    return [ meta, fastqs.flatten() ]
        }
}
