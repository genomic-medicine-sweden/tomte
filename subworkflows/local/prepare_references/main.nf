//
// Prepare reference files
//

include { GUNZIP as GUNZIP_FASTA                       } from '../../../modules/nf-core/gunzip/main'
include { GATK4_BEDTOINTERVALLIST as BEDTOINTERVALLIST } from '../../../modules/nf-core/gatk4/bedtointervallist/main'
include { GATK4_CREATESEQUENCEDICTIONARY as BUILD_DICT } from '../../../modules/nf-core/gatk4/createsequencedictionary/main'
include { GUNZIP as GUNZIP_GTF                         } from '../../../modules/nf-core/gunzip/main'
include { GET_RRNA_TRANSCRIPTS                         } from '../../../modules/local/get_rrna_transcripts/main'
include { GTFTOGENEPRED_REFFLAT as GTF_TO_REFFLAT      } from '../../../modules/local/gtftorefflat'
include { GET_CHROM_SIZES                              } from '../../../modules/local/get_chrom_sizes/main'
include { GUNZIP as GUNZIP_TRFASTA                     } from '../../../modules/nf-core/gunzip/main'
include { GFFREAD                                      } from '../../../modules/local/gffread'
include { SAMTOOLS_FAIDX as SAMTOOLS_FAIDX_GENOME      } from '../../../modules/nf-core/samtools/faidx/main'
include { UNTAR as UNTAR_STAR_INDEX                    } from '../../../modules/nf-core/untar/main'
include { STAR_GENOMEGENERATE as BUILD_STAR_GENOME     } from '../../../modules/nf-core/star/genomegenerate/main'
include { UNTAR as UNTAR_VEP_CACHE                     } from '../../../modules/nf-core/untar/main'
include { UNTAR as UNTAR_SALMON_INDEX                  } from '../../../modules/nf-core/untar/main'
include { SALMON_INDEX as SALMON_INDEX                 } from '../../../modules/nf-core/salmon/index/main'

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

    main:
        ch_versions = Channel.empty()

        // Gunzip fasta if necessary
        GUNZIP_FASTA(ch_fasta)
        ch_fasta_mix = branchChannelToCompressedAndUncompressed(ch_fasta)

        ch_fasta_final = ch_fasta_mix.uncompressed.mix(GUNZIP_FASTA.out.gunzip.collect())

        // If no genome indices, create it
        SAMTOOLS_FAIDX_GENOME(ch_fasta_final,[[],[]])
        ch_fai = ch_fai_input.mix(SAMTOOLS_FAIDX_GENOME.out.fai).collect()

        // If no dictionary, create it
        BUILD_DICT(ch_fasta_final)
        ch_dict = ch_sequence_dict_input.mix(BUILD_DICT.out.dict).collect()

        // Get chrom sizes
        GET_CHROM_SIZES( ch_fai )

        // Gunzip gtf if necessary
        GUNZIP_GTF(ch_gtf)
        ch_gtf_mix = branchChannelToCompressedAndUncompressed(ch_gtf)
        ch_gtf_final = ch_gtf_mix.uncompressed.mix(GUNZIP_GTF.out.gunzip.collect())

        // If no star index, create it
        BUILD_STAR_GENOME( ch_fasta_final, ch_gtf_final )
        // Untar star index if necessary
        UNTAR_STAR_INDEX(ch_star_index_input)
        ch_star_mix = branchChannelToCompressedAndUncompressed(ch_star_index_input)
        ch_star_mixed = ch_star_mix.uncompressed.mix(UNTAR_STAR_INDEX.out.untar.collect())
        ch_star_final = ch_star_mixed.mix(BUILD_STAR_GENOME.out.index.collect())

        // Convert gtf to refflat for picard
        GTF_TO_REFFLAT(ch_gtf_final)

        // Get rRNA transcripts and convert to interval_list format
        GET_RRNA_TRANSCRIPTS(ch_gtf_final)
        BEDTOINTERVALLIST( GET_RRNA_TRANSCRIPTS.out.bed.map { it -> [ [id:it.name], it ] }, ch_dict )
        ch_interval = BEDTOINTERVALLIST.out.interval_list.map{ meta, interv -> [interv] }.collect()

        // Preparing transcript fasta
        ch_fasta_fai = ch_fasta_final.join(ch_fai).collect()
        GFFREAD(ch_gtf_final, ch_fasta_fai)

        // Gunzip transcript fasta if necessary
        GUNZIP_TRFASTA ( ch_transcript_fasta_input.map { it -> [[:], it] } )
        ch_transcript_fasta_mix = branchChannelToCompressedAndUncompressed(ch_transcript_fasta_input)
        ch_transcript_fasta_mixed = ch_transcript_fasta_mix.uncompressed.mix(GUNZIP_TRFASTA.out.gunzip.map{meta, index -> index}.collect())
        ch_transcript_fasta_final = ch_transcript_fasta_mixed.mix(GFFREAD.out.tr_fasta.collect())

        // If no salmon index, create it
        SALMON_INDEX(ch_fasta_final.map{ meta, fasta -> [ fasta ] }, ch_transcript_fasta_final)
        // Untar salmon index if necessary
        UNTAR_SALMON_INDEX( ch_salmon_index_input.map { it -> [[:], it] } )
        ch_salmon_mix = branchChannelToCompressedAndUncompressed(ch_salmon_index_input)
        ch_salmon_mixed = ch_salmon_mix.uncompressed.mix(UNTAR_SALMON_INDEX.out.untar.map{meta, index -> index}.collect())
        ch_salmon_final = ch_salmon_mixed.mix(SALMON_INDEX.out.index.collect())

        // Untar vep chache is necesary
        UNTAR_VEP_CACHE (ch_vep_cache_input.map { vep_cache -> [[id:'vep_cache'], vep_cache] })
        ch_vep_cache_mix = branchChannelToCompressedAndUncompressed(ch_vep_cache_input)
        ch_final_vep = ch_vep_cache_mix.uncompressed.mix(UNTAR_VEP_CACHE.out.untar.map{meta, vep_cache -> vep_cache}.collect())

        ch_versions = ch_versions.mix(GUNZIP_FASTA.out.versions)
        ch_versions = ch_versions.mix(SAMTOOLS_FAIDX_GENOME.out.versions)
        ch_versions = ch_versions.mix(BUILD_DICT.out.versions)
        ch_versions = ch_versions.mix(GET_CHROM_SIZES.out.versions)
        ch_versions = ch_versions.mix(GUNZIP_GTF.out.versions)
        ch_versions = ch_versions.mix(GUNZIP_TRFASTA.out.versions)
        ch_versions = ch_versions.mix(UNTAR_STAR_INDEX.out.versions)
        ch_versions = ch_versions.mix(BUILD_STAR_GENOME.out.versions)
        ch_versions = ch_versions.mix(GTF_TO_REFFLAT.out.versions)
        ch_versions = ch_versions.mix(GET_RRNA_TRANSCRIPTS.out.versions)
        ch_versions = ch_versions.mix(BEDTOINTERVALLIST.out.versions)
        ch_versions = ch_versions.mix(GFFREAD.out.versions)
        ch_versions = ch_versions.mix(UNTAR_SALMON_INDEX.out.versions)
        ch_versions = ch_versions.mix(SALMON_INDEX.out.versions)
        ch_versions = ch_versions.mix(UNTAR_VEP_CACHE.out.versions)

    emit:
        chrom_sizes   = GET_CHROM_SIZES.out.sizes.collect()    // channel: [ path(sizes) ]
        fasta         = ch_fasta_final.collect()               // channel: [ val(meta), path(fasta) ]
        fai           = ch_fai                                 // channel: [ val(meta), path(fai) ]
        fasta_fai     = ch_fasta_fai                           // channel: [ val(meta), path(fasta), path(fai) ]
        gtf           = ch_gtf_final.collect()                 // channel: [ val(meta), path(gtf) ]
        sequence_dict = ch_dict                                // channel: [ val(meta), path(dict) ]
        star_index    = ch_star_final.collect()                // channel: [ val(meta), path(star_index) ]
        salmon_index  = ch_salmon_final.collect()              // channel: [ path(salmon_index) ]
        refflat       = GTF_TO_REFFLAT.out.refflat.collect()   // channel: [ path(refflat) ]
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
