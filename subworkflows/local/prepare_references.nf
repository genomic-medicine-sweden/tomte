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
        fasta            // parameter: [mandatory] path to fasta
        fai              // channel:   [optional]  [ val(meta), path(fai) ]
        star_index       // parameter: [optional]  path to star_index
        gtf              // parameter: [mandatory] path to gtf
        ch_vep_cache     // channel:   [optional]  [ path(vep_cache) ]
        transcript_fasta // parameter: [optional]  path to transcript_fasta
        salmon_index     // parameter: [optional]  path to salmon_index

    main:
        ch_versions = Channel.empty()

        fasta_meta = Channel.fromPath(fasta).map { it -> [ [id:it[0] ], it] }.collect()
        GUNZIP_FASTA(fasta_meta)
        ch_fasta = fasta.endsWith(".gz") ? GUNZIP_FASTA.out.gunzip.collect() : fasta_meta


        // If no genome indices, create it
        SAMTOOLS_FAIDX_GENOME(ch_fasta,[[],[]])
        ch_fai = fai.mix(SAMTOOLS_FAIDX_GENOME.out.fai).collect()

        BUILD_DICT(ch_fasta)
        ch_dict = BUILD_DICT.out.dict.collect()

        gtf_meta = Channel.fromPath(gtf).map{ it -> [ [id:it[0]], it ] }.collect()
        GUNZIP_GTF(gtf_meta)
        ch_gtf = gtf.endsWith(".gz") ? GUNZIP_GTF.out.gunzip.collect() : gtf_meta.collect()

        // Get chrom sizes
        GET_CHROM_SIZES( ch_fai )

        ch_fasta_no_meta = ch_fasta.map{ meta, fasta -> [ fasta ] }

        ch_star = star_index ?
            Channel.fromPath(star_index).collect().map { it -> [ [id:it[0].simpleName], it ] }
            : Channel.empty()

        BUILD_STAR_GENOME (ch_fasta, ch_gtf )
        UNTAR_STAR_INDEX( ch_star.map { it -> [[:], it[1]] } )
        ch_star_index = (!star_index) ?  BUILD_STAR_GENOME.out.index.collect() :
                                        (star_index.endsWith(".gz") ? UNTAR_STAR_INDEX.out.untar.map { it[1] } : ch_star)
        // Convert gtf to refflat for picard
        GTF_TO_REFFLAT(ch_gtf)

        // Get rRNA transcripts and convert to interval_list format
        GET_RRNA_TRANSCRIPTS(ch_gtf)

        BEDTOINTERVALLIST( GET_RRNA_TRANSCRIPTS.out.bed.map { it -> [ [id:it.name], it ] }, ch_dict )
        ch_interval = BEDTOINTERVALLIST.out.interval_list.map{ meta, interv -> [interv] }.collect()

        UNTAR_VEP_CACHE (ch_vep_cache)
        ch_untar_vep = UNTAR_VEP_CACHE.out.untar.map{ meta, files -> [files] }.collect()

        // Preparing transcript fasta
        ch_fasta_fai = ch_fasta.mix(ch_fai.map{meta, fai -> fai}).collect()

        GFFREAD(ch_gtf, ch_fasta_fai)
        transcript_fasta_no_meta = (!transcript_fasta) ? GFFREAD.out.tr_fasta.collect() :
                                    (transcript_fasta.endsWith(".gz") ? GUNZIP_TRFASTA.out.gunzip.collect().map{ meta, fasta -> [ fasta ] } : transcript_fasta)

        // Setting up Salmon index
        ch_salmon = salmon_index ? Channel.fromPath(salmon_index).collect() : Channel.empty()
        UNTAR_SALMON_INDEX( ch_salmon.map { it -> [[:], it] } )
        SALMON_INDEX(ch_fasta_no_meta, transcript_fasta_no_meta)

        ch_salmon_index = (!salmon_index) ? SALMON_INDEX.out.index.collect() :
                                            (salmon_index.endsWith(".gz") ? UNTAR_SALMON_INDEX.out.untar.map { it[1] } : salmon_index.collect())

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
        chrom_sizes   = GET_CHROM_SIZES.out.sizes.collect()        // channel: [ path(sizes) ]
        fasta         = ch_fasta                                   // channel: [ val(meta), path(fasta) ]
        fai           = ch_fai                                     // channel: [ val(meta), path(fai) ]
        fasta_fai     = ch_fasta_fai                               // channel: [ val(meta), path(fasta), path(fai) ]
        sequence_dict = BUILD_DICT.out.dict.collect()              // channel: [ val(meta), path(dict) ]
        gtf           = ch_gtf                                     // channel: [ val(meta), path(gtf) ]
        star_index    = ch_star_index                              // channel: [ val(meta), path(star_index) ]
        salmon_index  = ch_salmon_index                            // channel: [ path(salmon_index) ]
        refflat       = GTF_TO_REFFLAT.out.refflat.collect()       // channel: [ path(refflat) ]
        rrna_bed      = GET_RRNA_TRANSCRIPTS.out.bed.collect()     // channel: [ path(bed) ]
        interval_list = ch_interval                                // channel: [ path(interval) ]
        vep_resources = ch_untar_vep                               // channel: [ path(cache) ]
        versions      = ch_versions                                // channel: [ path(versions.yml) ]
}
