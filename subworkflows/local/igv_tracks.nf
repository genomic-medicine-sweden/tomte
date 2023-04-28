//
// IGV TRACKS
//

include { JUNCTION_TRACK   } from '../../modules/local/junction_track'
include { UCSC_WIGTOBIGWIG } from '../../modules/nf-core/ucsc/wigtobigwig/main'                                                                          

workflow IGV_TRACKS {
    take:
        wig
        chrom_sizes
        spl_junc

    main:
        ch_versions = Channel.empty()

        // Selecting the Signal.UniqueMultiple.str1.out.wig file
        ch_wig = wig
            .map{ meta, wigs ->
            return[meta, wigs[1]] }

        UCSC_WIGTOBIGWIG(
            ch_wig,
            chrom_sizes
        )
        ch_versions = ch_versions.mix(UCSC_WIGTOBIGWIG.out.versions.first())

        JUNCTION_TRACK(
            spl_junc
        )
        ch_versions = ch_versions.mix(JUNCTION_TRACK.out.versions.first())



    emit:
        bw       = UCSC_WIGTOBIGWIG.out.bw
        bed      = JUNCTION_TRACK.out.bed
        versions = ch_versions
}
