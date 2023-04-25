//
// IGV TRACKS
//

include { JUNCTION_TRACK   } from '../../modules/nf-core/local/junction_track'
include { UCSC_WIGTOBIGWIG } from '../modules/nf-core/ucsc/wigtobigwig/main'                                                                          

workflow IGV_TRACKS {
    take:
        bam
        fasta_no_meta
        refflat
        rrna_intervals

    main:
        ch_versions = Channel.empty()

        UCSC_WIGTOBIGWIG(
            
        )
        ch_versions = ch_versions.mix(UCSC_WIGTOBIGWIG.out.versions.first())

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
