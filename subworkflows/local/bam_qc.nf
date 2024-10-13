//
// BAM QC
//

include { PICARD_COLLECTRNASEQMETRICS     } from '../../modules/nf-core/picard/collectrnaseqmetrics/main'
include { PICARD_COLLECTINSERTSIZEMETRICS } from '../../modules/nf-core/picard/collectinsertsizemetrics/main'

workflow BAM_QC {
    take:
        ch_bam            // channel: [mandatory] [ val(meta), path(bam) ]
        ch_fasta          // channel: [mandatory] [ val(meta), path(fasta) ]
        ch_refflat        // channel: [mandatory] [ path(refflat) ]
        ch_rrna_intervals // channel: [mandatory] [ path(intervals) ]

    main:
        ch_versions = Channel.empty()

        PICARD_COLLECTRNASEQMETRICS(
            ch_bam,
            ch_refflat,
            ch_fasta.map{ meta, fasta -> fasta },
            ch_rrna_intervals
        )

        PICARD_COLLECTINSERTSIZEMETRICS(
            ch_bam
        )

        ch_versions = ch_versions.mix(PICARD_COLLECTRNASEQMETRICS.out.versions.first())
        ch_versions = ch_versions.mix(PICARD_COLLECTINSERTSIZEMETRICS.out.versions.first())

    emit:
        metrics_general_rna = PICARD_COLLECTRNASEQMETRICS.out.metrics     // channel: [  val(meta), path(txt) ]
        metrics_insert_size = PICARD_COLLECTINSERTSIZEMETRICS.out.metrics // channel: [  val(meta), path(txt) ]
        versions = ch_versions                                            // channel: [ path(versions.yml) ]
}
