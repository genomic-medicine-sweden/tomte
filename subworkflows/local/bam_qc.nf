//
// BAM QC
//

include { PICARD_COLLECTRNASEQMETRICS } from '../../modules/nf-core/picard/collectrnaseqmetrics/main'

workflow BAM_QC {
    take:
        bam             // channel (mandatory): [ path(bam) ]
        fasta_no_meta   // channel (mandatory): [ path(fasta) ]
        refflat         // channel (mandatory): [ path(refflat) ]
        rrna_intervals  // channel (mandatory): [ path(interval_list) ]

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
        metrics  = PICARD_COLLECTRNASEQMETRICS.out.metrics // channel: [ val(meta), path(rna_metrics) ]
        versions = ch_versions                             // channel: [ path(versions.yml) ]
}
