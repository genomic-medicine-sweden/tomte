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


workflow PREPARE_REFERENCES {
    take:
        fasta
        star_index
        gtf
        ch_vep_cache

    main:
        ch_versions = Channel.empty()

        // Prepare fasta file
        ch_fasta_meta = Channel.fromPath(fasta).map{ it -> [ [id:it.simpleName], it ] }.collect()
        if ( fasta.endsWith(".gz") ) {
            ch_fasta_meta = GUNZIP_FASTA(ch_fasta_meta).gunzip
            ch_versions = ch_versions.mix(GUNZIP_FASTA.out.versions)
        }
        ch_fasta_no_meta =  ch_fasta_meta.map{ meta, fasta -> [ fasta ] }

        // Genome indices
        SAMTOOLS_FAIDX_GENOME(ch_fasta_meta)
        ch_versions = ch_versions.mix(SAMTOOLS_FAIDX_GENOME.out.versions)

        BUILD_DICT(ch_fasta_no_meta)
        ch_dict = BUILD_DICT.out.dict.collect()
        ch_versions = ch_versions.mix(BUILD_DICT.out.versions)

        gtf_meta=channel.of(gtf).map{it -> [[id:it[0]], it]}.collect()
        GUNZIP_GTF(gtf_meta)
        ch_gtf  = GUNZIP_GTF.out.gunzip ? GUNZIP_GTF.out.gunzip.map{ meta, gtf -> [gtf] }.collect() : Channel.fromPath(gtf)
        ch_versions = ch_versions.mix(GUNZIP_GTF.out.versions)

        // Get chrom sizes
        GET_CHROM_SIZES( SAMTOOLS_FAIDX_GENOME.out.fai )
        ch_versions = ch_versions.mix(GET_CHROM_SIZES.out.versions)

        // Setting up STAR index channel
        ch_star_index = star_index ? Channel.fromPath(star_index).collect() : Channel.empty()
        if ( !star_index ) {
            ch_star_index = BUILD_STAR_GENOME (ch_fasta_no_meta, ch_gtf).index
            ch_versions = ch_versions.mix(BUILD_STAR_GENOME.out.versions)
        }
        else if( star_index && star_index.endsWith(".gz") ) {
            ch_star_index = UNTAR_STAR_INDEX( ch_star_index.map { it -> [[:], it] } ).untar.map { it[1] }
            ch_versions = ch_versions.mix(UNTAR_STAR_INDEX.out.versions)
        }

        // Convert gtf to refflat for picard
        GTF_TO_REFFLAT(ch_gtf)
        ch_versions = ch_versions.mix(GTF_TO_REFFLAT.out.versions)

        // Get rRNA transcripts and convert to interval_list format
        GET_RRNA_TRANSCRIPTS(ch_gtf)
        ch_versions = ch_versions.mix(GET_RRNA_TRANSCRIPTS.out.versions)

        BEDTOINTERVALLIST( GET_RRNA_TRANSCRIPTS.out.bed.map {it -> [ [id:it.name], it ]}, ch_dict )
        ch_versions = ch_versions.mix(BEDTOINTERVALLIST.out.versions)

        UNTAR_VEP_CACHE (ch_vep_cache)
        ch_versions = ch_versions.mix(UNTAR_VEP_CACHE.out.versions)

    emit:
        chrom_sizes    = GET_CHROM_SIZES.out.sizes.collect()                                 // channel: [ path(sizes) ]
        fasta_meta     = ch_fasta_meta
        fasta_no_meta  = ch_fasta_no_meta
        fasta_fai      = SAMTOOLS_FAIDX_GENOME.out.fai.map{ meta, fai -> [fai] }.collect()
        fasta_fai_meta = ch_fasta_meta.join(SAMTOOLS_FAIDX_GENOME.out.fai).collect()
        sequence_dict  = BUILD_DICT.out.dict.collect()
        gtf            = ch_gtf
        star_index     = ch_star_index
        refflat        = GTF_TO_REFFLAT.out.refflat.collect()
        rrna_bed       = GET_RRNA_TRANSCRIPTS.out.bed.collect()
        interval_list  = BEDTOINTERVALLIST.out.interval_list.map{ meta, interv -> [interv] }.collect()
        vep_resources  = UNTAR_VEP_CACHE.out.untar.map{meta, files -> [files]}.collect()     // channel: [ path(cache) ]
        versions       = ch_versions
}
