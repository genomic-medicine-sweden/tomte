//
// IGV TRACKS
//

include { JUNCTION_TRACK   } from '../../modules/local/junction_track'
include { UCSC_WIGTOBIGWIG } from '../../modules/nf-core/ucsc/wigtobigwig/main'
include { TABIX_BGZIPTABIX } from '../../modules/nf-core/tabix/bgziptabix/main'

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

        // Making Bigwig file
        UCSC_WIGTOBIGWIG(
            ch_wig,
            chrom_sizes
        )
        ch_versions = ch_versions.mix(UCSC_WIGTOBIGWIG.out.versions.first())

        // Making junction track
        JUNCTION_TRACK(
            spl_junc
        )
        ch_versions = ch_versions.mix(JUNCTION_TRACK.out.versions.first())

        // Bgziping and creating tbi file for junction track
        TABIX_BGZIPTABIX(
            JUNCTION_TRACK.out.bed
        )
        ch_versions = ch_versions.mix(TABIX_BGZIPTABIX.out.versions.first())


    emit:
        bw       = UCSC_WIGTOBIGWIG.out.bw
        bed      = TABIX_BGZIPTABIX.out.gz_tbi
        versions = ch_versions
}
