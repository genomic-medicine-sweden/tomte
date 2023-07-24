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
        reads                   // channel (mandatory): [ val(meta), [ path(reads) ] ]
        star_index              // channel (mandatory): [ path(star_index) ]
        gtf                     // channel (mandatory): [ path(gtf) ]
        subsample_bed           // channel: [ path(subsample_bed) ] MODIFY TO INCLUDE META!!!!
        seed_frac               // channel: [ val(seed_frac) ]
        num_reads               // channel: [ val(num_reads) ]
        subsample_region_switch // channel: [ bol(subsample_region_switch) ]
        downsample_switch       // channel: [ bol(downsample_switch) ]
        salmon_index            // channel (mandatory): [ path(salmon_index) ]
        ch_genome_fasta         // channel (mandatory): [ val(meta), [ path(fasta) ] ]

    main:
        ch_versions = Channel.empty()

        CAT_FASTQ(reads)

        FASTP(
            CAT_FASTQ.out.reads,
            [],
            false,
            false
        )

        STAR_ALIGN(
            FASTP.out.reads,
            star_index,
            gtf,
            false,
            'illumina',
            false
        )

        SAMTOOLS_INDEX( STAR_ALIGN.out.bam )

        ch_bam_bai = Channel.empty()
        ch_bam_bai_out = Channel.empty()

        // SHOULD BE PRETIFIED!!!
        // If subsample_region_switch is true and subsample_bed (i.e. hemoglobin regions) is provided,
        // only a fraction of the reads in the area will be kept
        // If downsample_switch is true the number of reads remaining after subsampleing will be
        // downsampled
        if (subsample_region_switch) {
            RNA_SUBSAMPLE_REGION(
                STAR_ALIGN.out.bam,
                subsample_bed,
                seed_frac
            )
            ch_bam_bai = ch_bam_bai.mix(RNA_SUBSAMPLE_REGION.out.bam_bai)
            if (!downsample_switch) {
                ch_bam_bai_out = ch_bam_bai.mix(RNA_SUBSAMPLE_REGION.out.bam_bai)
            }
        } else {
            ch_bam_bai = ch_bam_bai.mix(STAR_ALIGN.out.bam.join(SAMTOOLS_INDEX.out.bai))
            if (!downsample_switch) {
                ch_bam_bai_out = STAR_ALIGN.out.bam.join(SAMTOOLS_INDEX.out.bai)
            }
        }

        // If save_mapped_as_cram is true compress bam to cram
        SAMTOOLS_VIEW(
            STAR_ALIGN.out.bam.join(SAMTOOLS_INDEX.out.bai),
            ch_genome_fasta,
            []
        )

        // If downsample_switch is true the number of reads remaining after subsampleing will be
        // downsampled
        if (downsample_switch) {
            RNA_DOWNSAMPLE(
                ch_bam_bai,
                num_reads
            )
            ch_bam_bai_out = ch_bam_bai.mix(RNA_DOWNSAMPLE.out.bam_bai)
        }

        SALMON_QUANT(
            FASTP.out.reads,
            salmon_index,
            gtf,
            [],
            false,
            'A'
        )

        ch_versions = ch_versions.mix(CAT_FASTQ.out.versions.first())
        ch_versions = ch_versions.mix(FASTP.out.versions.first())
        ch_versions = ch_versions.mix(STAR_ALIGN.out.versions.first())
        ch_versions = ch_versions.mix(SAMTOOLS_INDEX.out.versions.first())
        ch_versions = ch_versions.mix(RNA_SUBSAMPLE_REGION.out.versions.first())
        ch_versions = ch_versions.mix(SAMTOOLS_VIEW.out.versions.first())
        ch_versions = ch_versions.mix(RNA_DOWNSAMPLE.out.versions.first())
        ch_versions = ch_versions.mix(SALMON_QUANT.out.versions.first())

    emit:
        merged_reads    = CAT_FASTQ.out.reads              // channel: [ val(meta), [ path(reads) ] ]
        fastp_report    = FASTP.out.json                   // channel: [ val(meta), [ path(json) ] ]
        bam             = STAR_ALIGN.out.bam               // channel: [ path(bam) ]
        bam_bai         = ch_bam_bai                       // channel: [ val(meta), [ path(bam) ], [ path(bai) ] ]
        bam_ds_bai      = ch_bam_bai_out                   // channel: [ val(meta), [ path(bam) ], [ path(bai) ] ]
        gene_counts     = STAR_ALIGN.out.read_per_gene_tab // channel: [ val(meta), [ path(read_per_gene_tab)] ]
        spl_junc        = STAR_ALIGN.out.spl_junc_tab      // channel: [ val(meta), [ path(spl_junc_tab)] ]
        star_log_final  = STAR_ALIGN.out.log_final         // channel: [ val(meta), [ path(log_final)] ]
        star_wig        = STAR_ALIGN.out.wig
        salmon_result   = SALMON_QUANT.out.results         // channel: [ val(meta), [ path(quant)] ]
        salmon_info     = SALMON_QUANT.out.json_info       // channel: [ val(meta), [ path(json)] ]
        versions        = ch_versions
}
