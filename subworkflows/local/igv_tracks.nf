//
// IGV TRACKS
//

include { JUNCTION_TRACK   } from '../../modules/local/junction_track'
include { UCSC_WIGTOBIGWIG } from '../../modules/nf-core/ucsc/wigtobigwig/main'
include { TABIX_BGZIPTABIX } from '../../modules/nf-core/tabix/bgziptabix/main'

workflow IGV_TRACKS {
    take:
        wig         //
        chrom_sizes // channel (mandatory): [ path(sizes) ]
        spl_junc    // channel (mandatory): [ val(meta), [ path(spl_junc_tab)] ]

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

        // Making junction track
        JUNCTION_TRACK(
            spl_junc
        )

        // Bgziping and creating tbi file for junction track
        TABIX_BGZIPTABIX(
            JUNCTION_TRACK.out.bed
        )

        ch_versions = ch_versions.mix(UCSC_WIGTOBIGWIG.out.versions.first())
        ch_versions = ch_versions.mix(JUNCTION_TRACK.out.versions.first())
        ch_versions = ch_versions.mix(TABIX_BGZIPTABIX.out.versions.first())


    emit:
        bw       = UCSC_WIGTOBIGWIG.out.bw     // channel: [ val(meta), [ path(bigwig) ] ]
        bed      = TABIX_BGZIPTABIX.out.gz_tbi // channel: [ val(meta), [ path(bed_gz) ], [ path(bed_gz_tbi) ] ]
        versions = ch_versions                 // channel: [ path(versions.yml) ]
}
