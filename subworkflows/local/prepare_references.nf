//
// Prepare reference files
//

include { GATK4_BEDTOINTERVALLIST as BEDTOINTERVALLIST } from '../../modules/nf-core/gatk4/bedtointervallist/main'
include { GATK4_CREATESEQUENCEDICTIONARY as BUILD_DICT } from '../../modules/nf-core/gatk4/createsequencedictionary/main'
include { GET_RRNA_TRANSCRIPTS                         } from '../../modules/local/get_rrna_transcripts'
include { GTFTOGENEPRED_REFFLAT as GTF_TO_REFFLAT      } from '../../modules/local/gtftorefflat'
include { GET_CHROM_SIZES                              } from '../../modules/local/get_chrom_sizes'
include { GUNZIP as GUNZIP_FASTA                       } from '../../modules/nf-core/gunzip/main'
include { GUNZIP as GUNZIP_GTF                         } from '../../modules/nf-core/gunzip/main'
include { SAMTOOLS_FAIDX as SAMTOOLS_FAIDX_GENOME      } from '../../modules/nf-core/samtools/faidx/main'
include { STAR_GENOMEGENERATE as BUILD_STAR_GENOME     } from '../../modules/nf-core/star/genomegenerate/main'
include { UNTAR as UNTAR_STAR_INDEX                    } from '../../modules/nf-core/untar/main'
include { UNTAR as UNTAR_VEP_CACHE                     } from '../../modules/nf-core/untar/main'
include { SALMON_INDEX as SALMON_INDEX                 } from '../../modules/nf-core/salmon/index/main'
include { GUNZIP as GUNZIP_TRFASTA                     } from '../../modules/nf-core/gunzip/main'

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

        GUNZIP_FASTA(fasta)
        ch_fasta = GUNZIP_FASTA.out.gunzip ? GUNZIP_FASTA.out.gunzip.collect() : fasta
        ch_versions = ch_versions.mix(GUNZIP_FASTA.out.versions)

        // If no genome indices, create it
        SAMTOOLS_FAIDX_GENOME(ch_fasta,[[],[]])
        ch_versions = ch_versions.mix(SAMTOOLS_FAIDX_GENOME.out.versions)
        ch_fai = Channel.empty().mix(fai, SAMTOOLS_FAIDX_GENOME.out.fai).collect()

        BUILD_DICT(ch_fasta)
        ch_dict = BUILD_DICT.out.dict.collect()
        ch_versions = ch_versions.mix(BUILD_DICT.out.versions)

        gtf_meta=channel.of(gtf).map{it -> [[id:it[0]], it]}.collect()
        GUNZIP_GTF(gtf_meta)
        ch_gtf_no_meta  = GUNZIP_GTF.out.gunzip ? GUNZIP_GTF.out.gunzip.map{ meta, gtf -> [gtf] }.collect() : Channel.fromPath(gtf)
        ch_versions = ch_versions.mix(GUNZIP_GTF.out.versions)

        // Get chrom sizes
        GET_CHROM_SIZES( ch_fai )
        ch_versions = ch_versions.mix(GET_CHROM_SIZES.out.versions)

        ch_fasta_no_meta =  ch_fasta.map{ meta, fasta -> [ fasta ] }
        BUILD_STAR_GENOME (ch_fasta_no_meta, ch_gtf_no_meta)
        UNTAR_STAR_INDEX( star_index.map { it -> [[:], it] } )
        ch_star_index = UNTAR_STAR_INDEX.out.untar ? UNTAR_STAR_INDEX.out.untar.map {meta, it -> [it] }.collect()
                                                    : ( star_index ? star_index : BUILD_STAR_GENOME.out.index)

        // Convert gtf to refflat for picard
        GTF_TO_REFFLAT(ch_gtf_no_meta)
        ch_versions = ch_versions.mix(GTF_TO_REFFLAT.out.versions)

        // Get rRNA transcripts and convert to interval_list format
        GET_RRNA_TRANSCRIPTS(ch_gtf_no_meta)
        ch_versions = ch_versions.mix(GET_RRNA_TRANSCRIPTS.out.versions)

        BEDTOINTERVALLIST( GET_RRNA_TRANSCRIPTS.out.bed.map {it -> [ [id:it.name], it ]}, ch_dict )
        ch_versions = ch_versions.mix(BEDTOINTERVALLIST.out.versions)

        UNTAR_VEP_CACHE (ch_vep_cache)
        ch_versions = ch_versions.mix(UNTAR_VEP_CACHE.out.versions)

        // Setting up Salmon index
        GUNZIP_TRFASTA(transcript_fasta)
        ch_transcript_fasta = GUNZIP_TRFASTA.out.gunzip ? GUNZIP_TRFASTA.out.gunzip : transcript_fasta
        ch_transcript_fasta_no_meta = ch_transcript_fasta.map{ meta, fasta -> [ fasta ] }

        ch_salmon_index = salmon_index ? Channel.fromPath( salmon_index ).collect() : Channel.empty()
        if (!transcript_fasta) {
                // We would need to add gffread here but the module needs to be changed a lot, it can be done later
            } 
        if ( !salmon_index ) {
            ch_salmon_index = SALMON_INDEX(ch_fasta_no_meta, ch_transcript_fasta_no_meta).index
            ch_versions = ch_versions.mix(SALMON_INDEX.out.versions)
        }  else if( salmon_index && salmon_index.endsWith(".gz") ) {
            ch_salmon_index = UNTAR_SALMON_INDEX( ch_salmon_index.map { it -> [[:], it] } ).untar.map { it[1] }
        }

    emit:
        chrom_sizes      = GET_CHROM_SIZES.out.sizes.collect()                                 // channel: [ path(sizes) ]
        fasta_meta       = ch_fasta
        fasta_no_meta    = ch_fasta_no_meta
        fai              = ch_fai
        fai_no_meta      = ch_fai.map{ meta, fai -> [fai] }.collect()
        fasta_fai_meta   = ch_fasta.join(ch_fai).collect()
        sequence_dict    = BUILD_DICT.out.dict.collect()
        gtf              = ch_gtf_no_meta
        star_index       = ch_star_index
        salmon_index     = ch_salmon_index
        transcript_fasta = ch_transcript_fasta_no_meta
        refflat          = GTF_TO_REFFLAT.out.refflat.collect()
        rrna_bed         = GET_RRNA_TRANSCRIPTS.out.bed.collect()
        interval_list    = BEDTOINTERVALLIST.out.interval_list.map{ meta, interv -> [interv] }.collect()
        vep_resources    = UNTAR_VEP_CACHE.out.untar.map{meta, files -> [files]}.collect()     // channel: [ path(cache) ]
        versions         = ch_versions
}
