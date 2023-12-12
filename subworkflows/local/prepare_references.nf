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
        fasta
        fai
        star_index
        gtf
        ch_vep_cache
        transcript_fasta
        salmon_index

    main:
        ch_versions = Channel.empty()

        fasta_meta = Channel.fromPath(fasta).map { it -> [ [id:it[0] ], it] }.collect()
        GUNZIP_FASTA(fasta_meta)
        ch_fasta = fasta.endsWith(".gz") ? GUNZIP_FASTA.out.gunzip.collect() : fasta_meta


        // If no genome indices, create it
        SAMTOOLS_FAIDX_GENOME(ch_fasta,[[],[]])
        ch_fai = fai.mix(SAMTOOLS_FAIDX_GENOME.out.fai).collect()
        //ch_fai.view()

        BUILD_DICT(ch_fasta)
        ch_dict = BUILD_DICT.out.dict.collect()

        gtf_meta = Channel.fromPath(gtf).map{ it -> [ [id:it[0]], it ] }.collect()
        GUNZIP_GTF(gtf_meta)
        ch_gtf_no_meta = gtf.endsWith(".gz") ? GUNZIP_GTF.out.gunzip.map{ meta, gtf -> [gtf] }.collect() : Channel.fromPath(gtf)

        // Get chrom sizes
        GET_CHROM_SIZES( ch_fai )

        ch_fasta_no_meta = ch_fasta.map{ meta, fasta -> [ fasta ] }

        ch_gtf=ch_gtf_no_meta.map { it -> [[:], it] }
        ch_star = star_index ?
            Channel.fromPath(star_index).collect().map { it -> [ [id:it[0].simpleName], it ] }
            : Channel.empty()

        BUILD_STAR_GENOME (ch_fasta, ch_gtf )
        UNTAR_STAR_INDEX( ch_star.map { it -> [[:], it[1]] } )
        ch_star_index = (!star_index) ?  BUILD_STAR_GENOME.out.index.collect() :
                                        (star_index.endsWith(".gz") ? UNTAR_STAR_INDEX.out.untar.map { it[1] } : ch_star)
        // Convert gtf to refflat for picard
        GTF_TO_REFFLAT(ch_gtf_no_meta)

        // Get rRNA transcripts and convert to interval_list format
        GET_RRNA_TRANSCRIPTS(ch_gtf_no_meta)

        BEDTOINTERVALLIST( GET_RRNA_TRANSCRIPTS.out.bed.map { it -> [ [id:it.name], it ] }, ch_dict )

        UNTAR_VEP_CACHE (ch_vep_cache)

        // Preparing transcript fasta
        ch_fai.view()
        ch_fasta_fai = ch_fasta.mix(ch_fai.map{meta, fai -> fai}).collect()
        ch_fasta_fai.view()
        //ch_fasta.view()
        GFFREAD(ch_gtf_no_meta.map{ it -> [ [id:it[0].simpleName], it ] },ch_fasta_fai)
        transcript_fasta_no_meta = (!transcript_fasta) ? GFFREAD.out.tr_fasta.collect() :
                                    (transcript_fasta.endsWith(".gz") ? GUNZIP_TRFASTA.out.gunzip.collect().map{ meta, fasta -> [ fasta ] } : transcript_fasta)

        // Setting up Salmon index
        ch_salmon = salmon_index ? Channel.fromPath(salmon_index).collect() : Channel.empty()
        UNTAR_SALMON_INDEX( ch_salmon.map { it -> [[:], it] } )
        SALMON_INDEX(ch_fasta_no_meta, transcript_fasta_no_meta)

        ch_salmon_index = (!salmon_index) ? SALMON_INDEX.out.index.collect() :
                                            (salmon_index.endsWith(".gz") ? UNTAR_SALMON_INDEX.out.untar.map { it[1] } : salmon_index)

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
        chrom_sizes      = GET_CHROM_SIZES.out.sizes.collect()                                 // channel: [ path(sizes) ]
        fasta_meta       = ch_fasta
        fasta_no_meta    = ch_fasta_no_meta
        fai              = ch_fai
        fai_no_meta      = ch_fai.map{ meta, fai -> [fai] }.collect()
        fasta_fai_meta   = ch_fasta_fai
        sequence_dict    = BUILD_DICT.out.dict.collect()
        gtf              = ch_gtf_no_meta
        star_index       = ch_star_index
        salmon_index     = ch_salmon_index
        refflat          = GTF_TO_REFFLAT.out.refflat.collect()
        rrna_bed         = GET_RRNA_TRANSCRIPTS.out.bed.collect()
        interval_list    = BEDTOINTERVALLIST.out.interval_list.map{ meta, interv -> [interv] }.collect()
        vep_resources    = UNTAR_VEP_CACHE.out.untar.map{ meta, files -> [files] }.collect()     // channel: [ path(cache) ]
        versions         = ch_versions
}
