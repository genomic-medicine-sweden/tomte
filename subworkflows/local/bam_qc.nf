//
// BAM QC
//

include { PICARD_COLLECTRNASEQMETRICS } from '../../modules/nf-core/picard/collectrnaseqmetrics/main'

workflow BAM_QC {
    take:
        bam
        fasta_no_meta
        refflat
        rrna_intervals

    main:
        ch_versions = Channel.empty()

        PICARD_COLLECTRNASEQMETRICS(
            bam,
            refflat,
            fasta_no_meta,
            rrna_intervals
        )
        ch_versions = ch_versions.mix(PICARD_COLLECTRNASEQMETRICS.out.versions.first())

    emit:
        metrics  = PICARD_COLLECTRNASEQMETRICS.out.metrics
        versions = ch_versions
}
