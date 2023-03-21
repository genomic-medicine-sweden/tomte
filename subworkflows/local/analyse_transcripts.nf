//
// ANALYSE TRANSCRITPS
//

 include { STRINGTIE_STRINGTIE } from '../../modules/nf-core/stringtie/stringtie/main'
 include { GFFCOMPARE          } from '../../modules/nf-core/gffcompare/main'

workflow ANALYSE_TRANSCRIPTS {
    take:
        ch_bam            // channel (mandatory): [ val(meta), [ path(bam) ] ]
        ch_gtf            // channel (mandatory): [ path(gtf) ]
        ch_fasta_fai_meta // channel (mandatory): [ val(meta), [ path(fasta), path(fai) ]

    main:
        ch_versions = Channel.empty()

        STRINGTIE_STRINGTIE(
            ch_bam,
            ch_gtf
        )
        ch_versions = ch_versions.mix(STRINGTIE_STRINGTIE.out.versions.first())

        GFFCOMPARE(
            STRINGTIE_STRINGTIE.out.transcript_gtf,
            ch_fasta_fai_meta,
            ch_gtf.map{ gtf -> [ [id:gtf.simpleName], gtf ] }
        )
        ch_versions = ch_versions.mix(GFFCOMPARE.out.versions.first())

    emit:
        transcript_gtf = STRINGTIE_STRINGTIE.out.transcript_gtf // channel: [ val(meta), [ path(transctript_gtf)] ]
        abundance      = STRINGTIE_STRINGTIE.out.abundance      // channel: [ val(meta), [ path(abundance) ] ]
        coverage_gtf   = STRINGTIE_STRINGTIE.out.coverage_gtf   // channel: [ val(meta), [ path(coverage_gtf) ] ]
        annotated_gtf  = GFFCOMPARE.out.annotated_gtf           // channel: [ val(meta), [ path(annotated_gtf) ] ]
        stats_gtf      = GFFCOMPARE.out.stats                   // channel: [ val(meta), [ path(stats) ] ]
        versions       = ch_versions                            // channel: [ path(versions.yml) ]
}
