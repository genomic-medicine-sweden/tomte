//
// Prepare reference files
//

include { GUNZIP as GUNZIP_FASTA                       } from '../../../modules/nf-core/gunzip'
include { GATK4_BEDTOINTERVALLIST as BEDTOINTERVALLIST } from '../../../modules/nf-core/gatk4/bedtointervallist'
include { GATK4_CREATESEQUENCEDICTIONARY as BUILD_DICT } from '../../../modules/nf-core/gatk4/createsequencedictionary'
include { GUNZIP as GUNZIP_GTF                         } from '../../../modules/nf-core/gunzip'
include { GET_RRNA_TRANSCRIPTS                         } from '../../../modules/local/get_rrna_transcripts'
include { UCSC_GTFTOGENEPRED                           } from '../../../modules/nf-core/ucsc/gtftogenepred'
include { GET_CHROM_SIZES                              } from '../../../modules/local/get_chrom_sizes'
include { GUNZIP as GUNZIP_TRFASTA                     } from '../../../modules/nf-core/gunzip'
include { GFFREAD                                      } from '../../../modules/nf-core/gffread'
include { SAMTOOLS_FAIDX as SAMTOOLS_FAIDX_GENOME      } from '../../../modules/nf-core/samtools/faidx'
include { UNTAR as UNTAR_STAR_INDEX                    } from '../../../modules/nf-core/untar'
include { STAR_GENOMEGENERATE as BUILD_STAR_GENOME     } from '../../../modules/nf-core/star/genomegenerate'
include { UNTAR as UNTAR_VEP_CACHE                     } from '../../../modules/nf-core/untar'
include { UNTAR as UNTAR_SALMON_INDEX                  } from '../../../modules/nf-core/untar'
include { SALMON_INDEX as SALMON_INDEX                 } from '../../../modules/nf-core/salmon/index'

workflow PREPARE_REFERENCES {
    take:
    ch_fasta                  // channel: [mandatory] [ val(meta), path(fasta) ]
    ch_fai_input              // channel: [optional]  [ val(meta), path(fai) ]
    ch_star_index_input       // channel: [optional]  [ val(meta), path(star_index) ]
    ch_gtf                    // channel: [mandatory] [ val(meta), path(gtf) ]
    ch_vep_cache_input        // channel: [optional]  [ path(vep_cache) ]
    ch_transcript_fasta_input // channel: [optional]  [ path(transcript_fasta) ]
    ch_salmon_index_input     // channel: [optional]  [ path(salmon_index) ]
    ch_sequence_dict_input    // channel: [optional]  [ val(meta), path(dict) ]
    gunzip_fasta              // boolean: should we gunzip fasta
    gunzip_gtf                // boolean: should we gunzip gtf
    gunzip_transcript_fasta   // boolean: should we gunzip transcript_fasta
    untar_star_index          // boolean: should we untar star_index
    untar_salmon_index        // boolean: should we untar salmon_index
    untar_vep_cache           // boolean: should we untar vep_cache
    build_fai                 // boolean: should we build fai
    build_sequence_dict       // boolean: should we build sequence_dict
    build_transcript_fasta    // boolean: should we build transcript_fasta
    build_star_index          // boolean: should we build star_index
    build_salmon_index        // boolean: should we build salmon_index


    main:
    ch_versions = Channel.empty()

    // Gunzip fasta if necessary
    if ( gunzip_fasta ) {
        GUNZIP_FASTA( ch_fasta )
        ch_fasta_final = GUNZIP_FASTA.out.gunzip.collect()
        ch_versions = ch_versions.mix( GUNZIP_FASTA.out.versions )
    } else {
        ch_fasta_final = branchChannelToCompressedAndUncompressed(ch_fasta).uncompressed.collect()
    }
    // If no genome indices, create it
    if ( build_fai ) {
        SAMTOOLS_FAIDX_GENOME( ch_fasta_final,[[],[]] )
        ch_fai = SAMTOOLS_FAIDX_GENOME.out.fai.collect()
        ch_versions = ch_versions.mix( SAMTOOLS_FAIDX_GENOME.out.versions )
    } else {
        ch_fai = ch_fai_input.collect()
    }

    // If no dictionary, create it
    if ( build_sequence_dict ) {
        BUILD_DICT( ch_fasta_final )
        ch_dict_final = BUILD_DICT.out.dict.collect()
        ch_versions = ch_versions.mix( BUILD_DICT.out.versions )
    } else {
        ch_dict_final = ch_sequence_dict_input.collect()
    }

    // Get chrom sizes
    GET_CHROM_SIZES( ch_fai )

    // Gunzip gtf if necessary
    if ( gunzip_gtf ) {
        GUNZIP_GTF( ch_gtf )
        ch_gtf_final = GUNZIP_GTF.out.gunzip.collect()
        ch_versions = ch_versions.mix( GUNZIP_GTF.out.versions )
    } else {
        ch_gtf_final = branchChannelToCompressedAndUncompressed( ch_gtf ).uncompressed.collect()
    }

    // If no star index, create it
    if ( build_star_index ) {
        BUILD_STAR_GENOME( ch_fasta_final, ch_gtf_final )
        ch_star_final = BUILD_STAR_GENOME.out.index.collect()
        ch_versions = ch_versions.mix( BUILD_STAR_GENOME.out.versions )
    } else {
        // Untar star index if necessary
        if ( untar_star_index ) {
            UNTAR_STAR_INDEX( ch_star_index_input )
            ch_star_final = UNTAR_STAR_INDEX.out.untar.collect()
            ch_versions = ch_versions.mix( UNTAR_STAR_INDEX.out.versions )
        } else {
            ch_star_final = branchChannelToCompressedAndUncompressed( ch_star_index_input ).uncompressed.collect()
        }
    }

    // Convert gtf to refflat for picard
    UCSC_GTFTOGENEPRED(ch_gtf_final)

    // Get rRNA transcripts and convert to interval_list format
    GET_RRNA_TRANSCRIPTS(ch_gtf_final)
    BEDTOINTERVALLIST( GET_RRNA_TRANSCRIPTS.out.bed.map { it -> [ [id:it.name], it ] }, ch_dict_final )
    ch_interval = BEDTOINTERVALLIST.out.interval_list.map{ _meta, interv -> [interv] }.collect()

    // Prepare fasta fai
    ch_fasta_fai = ch_fasta_final.join( ch_fai ).collect()

    // Preparing transcript fasta
    if ( build_transcript_fasta ) {
        GFFREAD( ch_gtf_final, ch_fasta_final.map{ _meta, fasta -> fasta } )
        ch_transcript_fasta_final = GFFREAD.out.gffread_fasta.map{ _meta, gffread_fa -> gffread_fa }.collect()
        ch_versions = ch_versions.mix( GFFREAD.out.versions )
    } else {
        // Gunzip transcript fasta if necessary
        if ( gunzip_transcript_fasta ) {
            GUNZIP_TRFASTA ( ch_transcript_fasta_input.map { it -> [[:], it] } )
            ch_transcript_fasta_final = GUNZIP_TRFASTA.out.gunzip.map{ _meta, index -> index }.collect()
            ch_versions = ch_versions.mix( GUNZIP_TRFASTA.out.versions )
        } else {
            ch_transcript_fasta_final = branchChannelToCompressedAndUncompressed( ch_transcript_fasta_input ).uncompressed.collect()
        }
    }


    // If no salmon index, create it
    if ( build_salmon_index ) {
        SALMON_INDEX( ch_fasta_final.map{ _meta, fasta -> fasta }, ch_transcript_fasta_final )
        ch_salmon_final = SALMON_INDEX.out.index.collect()
        ch_versions = ch_versions.mix( SALMON_INDEX.out.versions )
    } else {
        // Untar salmon index if necessary
        if ( untar_salmon_index ) {
            UNTAR_SALMON_INDEX( ch_salmon_index_input.map { it -> [[:], it] } )
            ch_salmon_final = UNTAR_SALMON_INDEX.out.untar.map{ _meta, index -> index }.collect()
            ch_versions = ch_versions.mix( UNTAR_SALMON_INDEX.out.versions )
        } else {
            ch_salmon_final = branchChannelToCompressedAndUncompressed( ch_salmon_index_input ).uncompressed.collect()
        }
    }

    // Untar vep chache is necesary
    if ( untar_vep_cache ) {
        UNTAR_VEP_CACHE ( ch_vep_cache_input.map { vep_cache -> [ [ id:'vep_cache' ], vep_cache ] } )
        ch_final_vep = UNTAR_VEP_CACHE.out.untar.map{ _meta, vep_cache -> vep_cache }.collect()
        ch_versions = ch_versions.mix( UNTAR_VEP_CACHE.out.versions )
    } else {
        ch_final_vep = branchChannelToCompressedAndUncompressed( ch_vep_cache_input ).uncompressed.collect()
    }

    ch_versions = ch_versions.mix(GET_CHROM_SIZES.out.versions)
    ch_versions = ch_versions.mix(UCSC_GTFTOGENEPRED.out.versions)
    ch_versions = ch_versions.mix(GET_RRNA_TRANSCRIPTS.out.versions)
    ch_versions = ch_versions.mix(BEDTOINTERVALLIST.out.versions)

    emit:
    chrom_sizes   = GET_CHROM_SIZES.out.sizes.collect()    // channel: [ path(sizes) ]
    fasta         = ch_fasta_final                         // channel: [ val(meta), path(fasta) ]
    fai           = ch_fai                                 // channel: [ val(meta), path(fai) ]
    fasta_fai     = ch_fasta_fai                           // channel: [ val(meta), path(fasta), path(fai) ]
    gtf           = ch_gtf_final                           // channel: [ val(meta), path(gtf) ]
    sequence_dict = ch_dict_final                          // channel: [ val(meta), path(dict) ]
    star_index    = ch_star_final                          // channel: [ val(meta), path(star_index) ]
    salmon_index  = ch_salmon_final                        // channel: [ path(salmon_index) ]
    refflat       = UCSC_GTFTOGENEPRED.out.refflat
        .map{ _meta, refflat -> refflat }.collect()        // channel: [ path(refflat) ]
    rrna_bed      = GET_RRNA_TRANSCRIPTS.out.bed.collect() // channel: [ path(bed) ]
    interval_list = ch_interval                            // channel: [ path(interval) ]
    vep_cache     = ch_final_vep.collect()                 // channel: [ path(cache) ]
    versions      = ch_versions                            // channel: [ path(versions.yml) ]
}
// Custom functions
/**
* Branch a channel into different channels,
* depending on whether the path is compressed or not.
* The resulting channels get meta only if the original one had it.
*
* @param Channel that may contain meta
* @return Channel branched on whether the file is compressed or uncompressed
*/
def branchChannelToCompressedAndUncompressed(ch) {
    if (ch.flatten().count().map{it == (1)}) {
        return ch.branch { file ->
                compressed: file.join("").endsWith(".gz") // If the file ends with .gz
                    return file
                uncompressed: !(file.join("").endsWith(".gz")) // If the file doesn't end with .gz
                    return file
                }
    } else if (ch.flatten().count().map{it == (2)}) {
        return ch.branch{ meta, file ->
                compressed: file.join("").endsWith(".gz") // If the file ends with .gz
                    return [ meta, file ]
                uncompressed: !(file.join("").endsWith(".gz")) // If the file doesn't end with .gz
                    return [ meta, file ]
                }
    }
}
