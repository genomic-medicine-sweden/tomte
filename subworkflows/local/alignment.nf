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
        reads                 // channel:   [mandatory] [ val(meta), [path(reads)]  ]
        star_index            // channel:   [mandatory] [ val(meta), path(star_index) ]
        ch_gtf                // channel:   [mandatory] [ val(meta), path(gtf) ]
        ch_platform           // channel:   [mandatory] [ val(platform) ]
        subsample_bed         // channel:   [optional]  [ path(subsample_bed) ]
        seed_frac             // parameter: [optional]  default: 0.001
        num_reads             // parameter: [optional]  default: 120000000
        skip_subsample_region // parameter: [mandatory] default: true
        skip_downsample       // parameter: [mandatory] default: true
        salmon_index          // channel:   [mandatory] [ path(salmon_index) ]
        ch_genome_fasta       // channel:   [mandatory] [ val(meta), path(fasta) ]

    main:
        ch_versions = Channel.empty()

        reads
            .branch {
                fastq: it[1].any { it.toString().endsWith('.fastq.gz') || it.toString().endsWith('.fq.gz') }
                bam:   it[1].any { it.toString().endsWith('.bam') }
            }
            .set { ch_input_branch }
        
        ch_bam_reads = ch_input_branch.bam
        ch_fastq_reads = ch_input_branch.fastq

        ch_input_branch.fastq.subscribe { println "Fastq branch: $it" }
        ch_input_branch.bam.subscribe { println "Bam branch: $it" }

        ch_fastq = branchFastqToSingleAndMulti(ch_fastq_reads)


        CAT_FASTQ(ch_fastq.multiple_fq)
            .reads.mix(ch_fastq.single_fq)
            .set { ch_cat_fastq }

        FASTP(ch_cat_fastq, [], false, false, false)

        STAR_ALIGN(FASTP.out.reads, star_index, ch_gtf, false, ch_platform, false)

        ch_bam_aligned=ch_bam_reads.mix(STAR_ALIGN.out.bam_sorted_aligned)

        SAMTOOLS_INDEX( ch_bam_aligned )

        ch_bam_bai = Channel.empty()
        ch_bam_bai_out = Channel.empty()

        if (!skip_subsample_region) {
            RNA_SUBSAMPLE_REGION( ch_bam_aligned, subsample_bed, seed_frac)
            ch_bam_bai = ch_bam_bai.mix(RNA_SUBSAMPLE_REGION.out.bam_bai)
            ch_versions = ch_versions.mix(RNA_SUBSAMPLE_REGION.out.versions.first())
            if (skip_downsample) {
                ch_bam_bai_out = RNA_SUBSAMPLE_REGION.out.bam_bai
            } else {
                RNA_DOWNSAMPLE( ch_bam_bai, num_reads)
                ch_bam_bai_out = RNA_DOWNSAMPLE.out.bam_bai
                ch_versions = ch_versions.mix(RNA_DOWNSAMPLE.out.versions.first())
            }
        } else {
            ch_bam_bai = ch_bam_bai.mix(ch_bam_aligned.join(SAMTOOLS_INDEX.out.bai))
            if (skip_downsample) {
                ch_bam_bai_out = ch_bam_aligned.join(SAMTOOLS_INDEX.out.bai)
            } else {
                RNA_DOWNSAMPLE( ch_bam_bai, num_reads)
                ch_bam_bai_out = RNA_DOWNSAMPLE.out.bam_bai
                ch_versions = ch_versions.mix(RNA_DOWNSAMPLE.out.versions.first())
            }
        }

        SAMTOOLS_VIEW( ch_bam_aligned.join(SAMTOOLS_INDEX.out.bai), ch_genome_fasta, [] )

        SALMON_QUANT( FASTP.out.reads, salmon_index, ch_gtf.map{ meta, gtf ->  gtf  }, [], false, 'A')

        ch_versions = ch_versions.mix(CAT_FASTQ.out.versions.first())
        ch_versions = ch_versions.mix(FASTP.out.versions.first())
        ch_versions = ch_versions.mix(STAR_ALIGN.out.versions.first())
        ch_versions = ch_versions.mix(SAMTOOLS_INDEX.out.versions.first())
        ch_versions = ch_versions.mix(SAMTOOLS_VIEW.out.versions.first())
        ch_versions = ch_versions.mix(SALMON_QUANT.out.versions.first())

    emit:
        merged_reads    = CAT_FASTQ.out.reads                // channel: [ val(meta), path(fastq) ]
        fastp_report    = FASTP.out.json                     // channel: [ val(meta), path(json) ]
        bam             = ch_bam_aligned  // channel: [ val(meta), path(bam) ]
        bam_bai         = ch_bam_bai                         // channel: [ val(meta), path(bam), path(bai) ]
        bam_ds_bai      = ch_bam_bai_out                     // channel: [ val(meta), path(bam), path(bai) ]
        gene_counts     = STAR_ALIGN.out.read_per_gene_tab   // channel: [ val(meta), path(tsv) ]
        spl_junc        = STAR_ALIGN.out.spl_junc_tab        // channel: [ val(meta), path(tsv) ]
        star_log_final  = STAR_ALIGN.out.log_final           // channel: [ val(meta), path(log) ]
        star_wig        = STAR_ALIGN.out.wig                 // channel: [ val(meta), path(wig) ]
        salmon_result   = SALMON_QUANT.out.results           // channel: [ val(meta), path(results) ]
        salmon_info     = SALMON_QUANT.out.json_info         // channel: [ val(meta), path(json) ]
        versions        = ch_versions                        // channel: [ path(versions.yml) ]
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
