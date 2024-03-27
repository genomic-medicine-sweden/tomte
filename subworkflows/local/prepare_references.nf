//
// Prepare reference files
//

include { GATK4_BEDTOINTERVALLIST as BEDTOINTERVALLIST } from '../../modules/nf-core/gatk4/bedtointervallist/main'
include { GATK4_CREATESEQUENCEDICTIONARY as BUILD_DICT } from '../../modules/nf-core/gatk4/createsequencedictionary/main'
include { GET_RRNA_TRANSCRIPTS                         } from '../../modules/local/get_rrna_transcripts'
include { GTFTOGENEPRED_REFFLAT as GTF_TO_REFFLAT      } from '../../modules/local/gtftorefflat'
include { GET_CHROM_SIZES                              } from '../../modules/local/get_chrom_sizes'
include { GFFREAD                                      } from '../../modules/local/gffread'
include { SAMTOOLS_FAIDX as SAMTOOLS_FAIDX_GENOME      } from '../../modules/nf-core/samtools/faidx/main'
include { STAR_GENOMEGENERATE as BUILD_STAR_GENOME     } from '../../modules/nf-core/star/genomegenerate/main'
include { UNTAR as UNTAR_VEP_CACHE                     } from '../../modules/nf-core/untar/main'
include { SALMON_INDEX as SALMON_INDEX                 } from '../../modules/nf-core/salmon/index/main'

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

        // If no genome indices, create it
        SAMTOOLS_FAIDX_GENOME(ch_fasta,[[],[]])
        ch_fai = ch_fai_input.mix(SAMTOOLS_FAIDX_GENOME.out.fai).collect()

        // If no dictionary, create it
        BUILD_DICT(ch_fasta)
        ch_dict = ch_sequence_dict_input.mix(BUILD_DICT.out.dict).collect()

        // Get chrom sizes
        GET_CHROM_SIZES( ch_fai )

        // If no star index, create it
        BUILD_STAR_GENOME ( ch_fasta, ch_gtf )
        ch_star_index = ch_star_index_input.mix(BUILD_STAR_GENOME.out.index).collect()

        // Convert gtf to refflat for picard
        GTF_TO_REFFLAT(ch_gtf)

        // Get rRNA transcripts and convert to interval_list format
        GET_RRNA_TRANSCRIPTS(ch_gtf)
        BEDTOINTERVALLIST( GET_RRNA_TRANSCRIPTS.out.bed.map { it -> [ [id:it.name], it ] }, ch_dict )
        ch_interval = BEDTOINTERVALLIST.out.interval_list.map{ meta, interv -> [interv] }.collect()

        // Preparing transcript fasta
        ch_fasta_fai = ch_fasta.mix(ch_fai.map{meta, fai -> fai}).collect()
        GFFREAD(ch_gtf, ch_fasta_fai)
        ch_transcript_fasta = ch_transcript_fasta_input.mix(GFFREAD.out.tr_fasta.collect()).collect()

        // If no Salmon index, create it
        SALMON_INDEX(ch_fasta.map{ meta, fasta -> [ fasta ] }, ch_transcript_fasta)
        ch_salmon_index = ch_salmon_index_input.mix(SALMON_INDEX.out.index).collect()

        // Untar vep chache is necesary
        UNTAR_VEP_CACHE (ch_vep_cache_input.map { it -> [[id:'vep_cache'], it] })
        ch_untar_vep = UNTAR_VEP_CACHE.out.untar.map{ meta, files -> [files] }.collect()
 
        ch_versions = ch_versions.mix(SAMTOOLS_FAIDX_GENOME.out.versions)
        ch_versions = ch_versions.mix(BUILD_DICT.out.versions)
        ch_versions = ch_versions.mix(GET_CHROM_SIZES.out.versions)
        ch_versions = ch_versions.mix(BUILD_STAR_GENOME.out.versions)
        ch_versions = ch_versions.mix(GTF_TO_REFFLAT.out.versions)
        ch_versions = ch_versions.mix(GET_RRNA_TRANSCRIPTS.out.versions)
        ch_versions = ch_versions.mix(BEDTOINTERVALLIST.out.versions)
        ch_versions = ch_versions.mix(GFFREAD.out.versions)
        ch_versions = ch_versions.mix(SALMON_INDEX.out.versions)
        ch_versions = ch_versions.mix(UNTAR_VEP_CACHE.out.versions)

    emit:
        chrom_sizes   = GET_CHROM_SIZES.out.sizes.collect()    // channel: [ path(sizes) ]
        fasta         = ch_fasta                               // channel: [ val(meta), path(fasta) ]
        fai           = ch_fai                                 // channel: [ val(meta), path(fai) ]
        fasta_fai     = ch_fasta_fai                           // channel: [ val(meta), path(fasta), path(fai) ]
        sequence_dict = ch_dict                                // channel: [ val(meta), path(dict) ]
        star_index    = ch_star_index                          // channel: [ val(meta), path(star_index) ]
        salmon_index  = ch_salmon_index                        // channel: [ path(salmon_index) ]
        refflat       = GTF_TO_REFFLAT.out.refflat.collect()   // channel: [ path(refflat) ]
        rrna_bed      = GET_RRNA_TRANSCRIPTS.out.bed.collect() // channel: [ path(bed) ]
        interval_list = ch_interval                            // channel: [ path(interval) ]
        vep_resources = ch_untar_vep                           // channel: [ path(cache) ]
        versions      = ch_versions                            // channel: [ path(versions.yml) ]
}
