//
// Call variants GATK
//

// Modules
include { GATK4_HAPLOTYPECALLER   } from '../../../modules/nf-core/gatk4/haplotypecaller/main'
include { GATK4_SPLITNCIGARREADS  } from '../../../modules/nf-core/gatk4/splitncigarreads/main'
include { GATK4_VARIANTFILTRATION } from '../../../modules/nf-core/gatk4/variantfiltration/main'
include { SAMTOOLS_INDEX          } from '../../../modules/nf-core/samtools/index/main'
include { BCFTOOLS_STATS          } from '../../../modules/nf-core/bcftools/stats/main'

workflow CALL_VARIANTS_GATK {
    take:
        ch_bam_bai // channel: [mandatory] [ val(meta), [ path(bam), path(bai) ] ]
        ch_fasta   // channel: [mandatory] [ val(meta), path(fasta) ]
        ch_fai     // channel: [mandatory] [ val(meta), path(fai) ]
        ch_dict    // channel: [mandatory] [ val(meta), path(dict) ]

    main:

        ch_versions = Channel.empty()

        GATK4_SPLITNCIGARREADS(
            ch_bam_bai.map{ meta, bam, bai -> [meta, bam, bai, []] },
            ch_fasta,
            ch_fai,
            ch_dict
        )
        ch_versions = ch_versions.mix(GATK4_SPLITNCIGARREADS.out.versions.first())

        SAMTOOLS_INDEX(
            GATK4_SPLITNCIGARREADS.out.bam
        )
        ch_versions = ch_versions.mix(SAMTOOLS_INDEX.out.versions.first())

        ch_split_bam_bai = GATK4_SPLITNCIGARREADS.out.bam.join(SAMTOOLS_INDEX.out.bai)

        GATK4_HAPLOTYPECALLER(
            ch_split_bam_bai.map{ meta, bam, bai -> [meta, bam, bai, [], []] },
            ch_fasta,
            ch_fai,
            ch_dict,
            [[],[]],
            [[],[]],
        )
        ch_versions = ch_versions.mix(GATK4_HAPLOTYPECALLER.out.versions.first())

        GATK4_VARIANTFILTRATION(
            GATK4_HAPLOTYPECALLER.out.vcf.join(GATK4_HAPLOTYPECALLER.out.tbi),
            ch_fasta,
            ch_fai,
            ch_dict
        )
        ch_versions = ch_versions.mix(GATK4_VARIANTFILTRATION.out.versions.first())

        BCFTOOLS_STATS(
            GATK4_HAPLOTYPECALLER.out.vcf.join(GATK4_HAPLOTYPECALLER.out.tbi),
            [[],[]],
            [[],[]],
            [[],[]],
            [[],[]],
            [[],[]],
        )
        ch_versions = ch_versions.mix(BCFTOOLS_STATS.out.versions.first())

    emit:
        vcf      = GATK4_HAPLOTYPECALLER.out.vcf // channel: [ val(meta), path(vcf) ]
        tbi      = GATK4_HAPLOTYPECALLER.out.tbi // channel: [ val(meta), path(tbi) ]
        stats    = BCFTOOLS_STATS.out.stats      // channel: [ val(meta), path(stats) ]
        versions = ch_versions                   // channel: [ path(versions.yml) ]
}
