//
// IGV TRACKS
//

include { JUNCTION_TRACK } from '../../modules/nf-core/local/junction_track'

workflow IGV_TRACKS {
    take:
        bam
        fasta_no_meta
        refflat
        rrna_intervals

    main:
        ch_versions = Channel.empty()
        

        JUNCTION_TRACK(
            bam,
            refflat,
            fasta_no_meta,
            rrna_intervals
        )
        ch_versions = ch_versions.mix(JUNCTION_TRACK.out.versions.first())

    emit:
        metrics  = JUNCTION_TRACK.out.bed
        versions = ch_versions
}
