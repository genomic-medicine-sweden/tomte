//
// Prepare reference files
//

include { GATK4_BEDTOINTERVALLIST as BEDTOINTERVALLIST } from '../../modules/nf-core/gatk4/bedtointervallist/main'
include { GATK4_CREATESEQUENCEDICTIONARY as BUILD_DICT } from '../../modules/nf-core/gatk4/createsequencedictionary/main'
include { GET_RRNA_TRANSCRIPTS                         } from '../../modules/local/get_rrna_transcripts'
include { GTFTOGENEPRED_REFFLAT as GTF_TO_REFFLAT      } from '../../modules/local/gtftorefflat'
include { GET_CHROM_SIZES                              } from '../../modules/local/get_chrom_sizes'
include { GFFREAD                                      } from '../../modules/local/gffread'
include { GUNZIP as GUNZIP_FASTA                       } from '../../modules/nf-core/gunzip/main'
include { GUNZIP as GUNZIP_GTF                         } from '../../modules/nf-core/gunzip/main'
include { GUNZIP as GUNZIP_TRFASTA                     } from '../../modules/nf-core/gunzip/main'
include { SAMTOOLS_FAIDX as SAMTOOLS_FAIDX_GENOME      } from '../../modules/nf-core/samtools/faidx/main'
include { STAR_GENOMEGENERATE as BUILD_STAR_GENOME     } from '../../modules/nf-core/star/genomegenerate/main'
include { UNTAR as UNTAR_STAR_INDEX                    } from '../../modules/nf-core/untar/main'
include { UNTAR as UNTAR_SALMON_INDEX                  } from '../../modules/nf-core/untar/main'
include { UNTAR as UNTAR_VEP_CACHE                     } from '../../modules/nf-core/untar/main'
include { SALMON_INDEX as SALMON_INDEX                 } from '../../modules/nf-core/salmon/index/main'

workflow PREPARE_REFERENCES {
    take:
        fasta             // channel (mandatory): [ val(meta), [ path(fasta) ] ]
        fai               // channel (mandatory): [ val(meta), [ path(fai) ] ]
        star_index        // channel : [ path(star_index) ]
        gtf               // channel : [ path(gtf) ]
        ch_vep_cache      // channel : [ val(meta), [ path(vep_cache) ] ]
        transcript_fasta  // channel : [ path(transcript_fasta) ]
        salmon_index      // channel : [ path(salmon_index) ]

    main:
        ch_versions = Channel.empty()

        // Unzip fasta if it ends in .gz
        GUNZIP_FASTA(fasta)
        ch_fasta = GUNZIP_FASTA.out.gunzip ? GUNZIP_FASTA.out.gunzip.collect() : fasta

        // Create genome indices if not provided
        SAMTOOLS_FAIDX_GENOME(ch_fasta,[[],[]])
        ch_fai = Channel.empty().mix(fai, SAMTOOLS_FAIDX_GENOME.out.fai).collect()

        // Create fasta.fai if not provided
        BUILD_DICT(ch_fasta)
        ch_dict = BUILD_DICT.out.dict.collect()

        // Unzip gtf if it ends in .gz
        gtf_meta=channel.of(gtf).map{it -> [[id:it[0]], it]}.collect()
        GUNZIP_GTF(gtf_meta)
        ch_gtf_no_meta  = GUNZIP_GTF.out.gunzip ? GUNZIP_GTF.out.gunzip.map{ meta, gtf -> [gtf] }.collect() : Channel.fromPath(gtf)

        // Get chrom sizes
        GET_CHROM_SIZES( ch_fai )

        ch_fasta_no_meta =  ch_fasta.map{ meta, fasta -> [ fasta ] }

        ch_star = star_index ? Channel.fromPath(star_index).collect() : Channel.empty()
        // Create star if not provided
        BUILD_STAR_GENOME (ch_fasta_no_meta, ch_gtf_no_meta)
        // Untar star index if it ends in .gz
        UNTAR_STAR_INDEX( ch_star.map { it -> [[:], it] } )
        ch_star_index = (!star_index) ?  BUILD_STAR_GENOME.out.index.collect() :
                                        (star_index.endsWith(".gz") ? UNTAR_STAR_INDEX.out.untar.map { it[1] }.collect() : star_index)

        // Convert gtf to refflat for picard
        GTF_TO_REFFLAT(ch_gtf_no_meta)

        // Get rRNA transcripts and convert to interval_list format
        GET_RRNA_TRANSCRIPTS(ch_gtf_no_meta)

        BEDTOINTERVALLIST( GET_RRNA_TRANSCRIPTS.out.bed.map {it -> [ [id:it.name], it ]}, ch_dict )

        // Untar VEP cache if it ends in .gz
        UNTAR_VEP_CACHE (ch_vep_cache)

        // Preparing transcript fasta
        ch_fasta_fai = ch_fasta.join(ch_fai).collect()
        ch_tr_fasta = transcript_fasta ? Channel.fromPath(transcript_fasta).map {it -> [[id:it[0].simpleName], it]}.collect() : Channel.empty()
        // Create transcript_fasta if not provided
        GFFREAD(ch_gtf_no_meta.map{it -> [[id:it[0].simpleName], it]},ch_fasta_fai)
        // Unzip transcript_fasta if it ends in .gz
        GUNZIP_TRFASTA(ch_tr_fasta)
        transcript_fasta_no_meta = (!transcript_fasta) ? GFFREAD.out.tr_fasta :
                                   (transcript_fasta.endsWith(".gz") ? GUNZIP_TRFASTA.out.gunzip.map{ meta, fasta -> [ fasta ] } : ch_tr_fasta.map{ meta, fasta -> [ fasta ] })

        // Setting up Salmon index
        ch_salmon = salmon_index ? Channel.fromPath(salmon_index).collect() : Channel.empty()
        // Untar salmon_index if it ends in .gz
        UNTAR_SALMON_INDEX( ch_salmon.map { it -> [[:], it] } )
        // Create salmon_index if not provided
        SALMON_INDEX(ch_fasta_no_meta, transcript_fasta_no_meta)
        ch_salmon_index = (!salmon_index) ? SALMON_INDEX.out.index.collect() :
                                            (salmon_index.endsWith(".gz") ? UNTAR_SALMON_INDEX.out.untar.map { it[1] }.collect() : salmon_index)

        ch_versions = ch_versions.mix(GUNZIP_FASTA.out.versions)
        ch_versions = ch_versions.mix(SAMTOOLS_FAIDX_GENOME.out.versions)
        ch_versions = ch_versions.mix(BUILD_DICT.out.versions)
        ch_versions = ch_versions.mix(GET_CHROM_SIZES.out.versions)
        ch_versions = ch_versions.mix(BUILD_STAR_GENOME.out.versions)
        ch_versions = ch_versions.mix(GTF_TO_REFFLAT.out.versions)
        ch_versions = ch_versions.mix(GET_RRNA_TRANSCRIPTS.out.versions)
        ch_versions = ch_versions.mix(BEDTOINTERVALLIST.out.versions)
        ch_versions = ch_versions.mix(GFFREAD.out.versions)
        ch_versions = ch_versions.mix(UNTAR_VEP_CACHE.out.versions)
        ch_versions = ch_versions.mix(SALMON_INDEX.out.versions)

    emit:
        chrom_sizes      = GET_CHROM_SIZES.out.sizes.collect()                                           // channel: [ path(sizes) ]
        fasta_meta       = ch_fasta                                                                      // channel: [ val(meta), [ path(fasta) ] ]
        fasta_no_meta    = ch_fasta_no_meta                                                              // channel: [ path(fasta) ]
        fai              = ch_fai                                                                        // channel: [ val(meta), [ path(fai) ] ]
        fai_no_meta      = ch_fai.map{ meta, fai -> [fai] }.collect()                                    // channel: [ path(fai) ]
        fasta_fai_meta   = ch_fasta.join(ch_fai).collect()                                               // channel: [ val(meta), [ path(fasta) ], [ path(fai) ] ]
        sequence_dict    = BUILD_DICT.out.dict.collect()                                                 // channel: [ val(meta), [ path(dict) ] ]
        gtf              = ch_gtf_no_meta                                                                // channel: [ path(gtf) ]
        star_index       = ch_star_index                                                                 // channel: [ path(star_index) ]
        salmon_index     = ch_salmon_index                                                               // channel: [ path(salmon_index) ]
        refflat          = GTF_TO_REFFLAT.out.refflat.collect()                                          // channel: [ path(refflat) ]
        rrna_bed         = GET_RRNA_TRANSCRIPTS.out.bed.collect()                                        // channel: [ path(rrna_bed) ]
        interval_list    = BEDTOINTERVALLIST.out.interval_list.map{ meta, interv -> [interv] }.collect() // channel: [ path(interval_list) ]
        vep_resources    = UNTAR_VEP_CACHE.out.untar.map{meta, files -> [files]}.collect()               // channel: [ path(cache) ]
        versions         = ch_versions                                                                   // channel: [ path(versions.yml) ]
}
