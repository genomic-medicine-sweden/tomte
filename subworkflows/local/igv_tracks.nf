//
// IGV TRACKS
//

include { JUNCTION_TRACK   } from '../../modules/local/junction_track'
include { UCSC_WIGTOBIGWIG } from '../../modules/nf-core/ucsc/wigtobigwig/main'                                                                          

workflow IGV_TRACKS {
    take:
        wig
        path_sizes
        gene_counts

    main:
        ch_versions = Channel.empty()

        UCSC_WIGTOBIGWIG(
            wig,
            path_sizes
        )
        ch_versions = ch_versions.mix(UCSC_WIGTOBIGWIG.out.versions.first())

        sj=gene_counts
            .map{meta, tab ->
            return[meta, tab[1]] }
        sj.view()
        JUNCTION_TRACK(
            sj
        )
        ch_versions = ch_versions.mix(JUNCTION_TRACK.out.versions.first())



    emit:
        bw       = UCSC_WIGTOBIGWIG.out.bw
        bed      = JUNCTION_TRACK.out.bed
        versions = ch_versions
}