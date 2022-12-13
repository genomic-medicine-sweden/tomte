//
// Allignment
//

include { CAT_FASTQ  } from '../../modules/nf-core/cat/fastq/main'
include { FASTP      } from '../../modules/nf-core/fastp/main'
include { STAR_ALIGN } from '../../modules/nf-core/star/align/main'

workflow ALLIGNMENT {
    take:
        reads
        star_index
        gtf

    main:
        ch_versions = Channel.empty()
        
        CAT_FASTQ(reads)
        ch_versions = ch_versions.mix(CAT_FASTQ.out.versions.first())

        FASTP(CAT_FASTQ.out.reads,[],false,false)
        ch_versions = ch_versions.mix(FASTP.out.versions.first())

        STAR_ALIGN(FASTP.out.reads, star_index, gtf, false, 'illumina', false)
        ch_versions = ch_versions.mix(STAR_ALIGN.out.versions.first())

    emit:
        ch_merged_reads = CAT_FASTQ.out.reads
        ch_bam = STAR_ALIGN.out.bam
}
