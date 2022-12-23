//
// Prepare reference files
//

include { SAMTOOLS_FAIDX as SAMTOOLS_FAIDX_GENOME      } from '../../modules/nf-core/samtools/faidx/main'
include { GATK4_CREATESEQUENCEDICTIONARY as BUILD_DICT } from '../../modules/nf-core/gatk4/createsequencedictionary/main'
include { UNTAR as UNTAR_STAR_INDEX                    } from '../../modules/nf-core/untar/main'
include { GUNZIP as GUNZIP_GTF                         } from '../../modules/nf-core/gunzip/main'
include { GTFTOGENEPRED_REFFLAT as GTF_TO_REFFLAT      } from '../../modules/local/gtftorefflat'
include { GET_RRNA_TRANSCRIPTS                         } from '../../modules/local/get_rrna_transcripts'
include { GATK4_BEDTOINTERVALLIST as BEDTOINTERVALLIST } from '../../modules/nf-core/gatk4/bedtointervallist/main'

workflow PREPARE_REFERENCES {
    take:
        fasta_no_meta
        fasta_meta
        star_index
        gtf


    main:
        ch_versions = Channel.empty()
        
        // Genome indices
        SAMTOOLS_FAIDX_GENOME(fasta_meta)
        ch_versions = ch_versions.mix(SAMTOOLS_FAIDX_GENOME.out.versions)

        BUILD_DICT(fasta_no_meta)
        ch_dict = BUILD_DICT.out.dict.collect()
        ch_versions = ch_versions.mix(BUILD_DICT.out.versions)

        star_index_meta=channel.of(star_index).map{it -> [[id:it[0]], it]}.collect()
        UNTAR_STAR_INDEX(star_index_meta)
        ch_versions = ch_versions.mix(UNTAR_STAR_INDEX.out.versions)

        gtf_meta=channel.of(gtf).map{it -> [[id:it[0]], it]}.collect()
        GUNZIP_GTF(gtf_meta)
        ch_gtf  = GUNZIP_GTF.out.gunzip ? GUNZIP_GTF.out.gunzip.map{ meta, gtf -> [gtf] }.collect() : Channel.fromPath(gtf)
        ch_versions = ch_versions.mix(GUNZIP_GTF.out.versions)

        GTF_TO_REFFLAT(ch_gtf)
        ch_versions = ch_versions.mix(GTF_TO_REFFLAT.out.versions)

        GET_RRNA_TRANSCRIPTS(ch_gtf)
        ch_versions = ch_versions.mix(GET_RRNA_TRANSCRIPTS.out.versions)

        bed_meta = GET_RRNA_TRANSCRIPTS.out.bed.map {it -> [[id:it[0]], it]}.collect()
        BEDTOINTERVALLIST(bed_meta,ch_dict)
        ch_versions = ch_versions.mix(BEDTOINTERVALLIST.out.versions)

    emit:
        fasta_fai     = SAMTOOLS_FAIDX_GENOME.out.fai.map{ meta, fai -> [fai] }.collect()
        sequence_dict = BUILD_DICT.out.dict.collect()
        star_index    = UNTAR_STAR_INDEX.out.untar.map{ meta, star -> [star] }.collect()
        gtf           = ch_gtf
        refflat       = GTF_TO_REFFLAT.out.refflat.collect()
        rrna_bed      = GET_RRNA_TRANSCRIPTS.out.bed.collect()
        interval_list = BEDTOINTERVALLIST.out.interval_list.map{ meta, interv -> [interv] }.collect()
}
