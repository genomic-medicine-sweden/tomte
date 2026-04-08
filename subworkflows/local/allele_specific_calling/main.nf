//
// Allele specific variant calling
//
include { BCFTOOLS_NORM                        } from '../../../modules/nf-core/bcftools/norm/main'
include { BCFTOOLS_VIEW                        } from '../../../modules/nf-core/bcftools/view/main'
include { GATK4_ASEREADCOUNTER                 } from '../../../modules/nf-core/gatk4/asereadcounter/main'
include { BOOTSTRAPANN                         } from '../../../modules/local/bootstrapann/main'
include { TABIX_BGZIPTABIX                     } from '../../../modules/nf-core/tabix/bgziptabix/main'
include { BCFTOOLS_MERGE                       } from '../../../modules/nf-core/bcftools/merge/main'
include { RENAME_FILES                         } from '../../../modules/local/rename_files/main'
include { TABIX_TABIX                          } from '../../../modules/nf-core/tabix/tabix/main'
include { BCFTOOLS_NORM as SPLIT_MULTIALLELICS } from '../../../modules/nf-core/bcftools/norm/main'
include { BCFTOOLS_NORM as REMOVE_DUPLICATES   } from '../../../modules/nf-core/bcftools/norm/main'
include { SPLIT_FAI_INTERVALS                  } from '../../../modules/local/split_fai_intervals/main'
include { SPLIT_BED_BY_CHROM                   } from '../../../modules/local/split_bed_by_chrom/main'
include { MERGE_ASE_CSVS                       } from '../../../modules/local/merge_ase_csvs/main'


workflow ALLELE_SPECIFIC_CALLING {
    take:
    ch_ind_vcf_tbi     // channel: [mandatory] [ val(meta), [ path(vcf), path(tbi) ] ]
    ch_bam_bai         // channel: [mandatory] [ val(meta), [ path(bam), path(bai) ] ]
    ch_fasta           // channel: [mandatory] [ val(meta), path(fasta) ]
    ch_fai             // channel: [mandatory] [ val(meta), path(fai) ]
    ch_dict            // channel: [mandatory] [ val(meta), path(dict) ]
    ch_case_info       // channel: [mandatory] [ val(case_info) ]
    ch_ase_intervals   // channel: [optional]  [ val(meta), path(bed) ] or empty channel
    use_intervals      // boolean: true if ase_intervals param is provided

    main:
    ch_versions = Channel.empty()

    // Keep only one variant per position in the vcf
    BCFTOOLS_NORM(
        ch_ind_vcf_tbi,
        ch_fasta
    )

    // Keep only does variants in the vcf that are SNVs and are heterozygote
    BCFTOOLS_VIEW(
        BCFTOOLS_NORM.out.vcf.join(BCFTOOLS_NORM.out.tbi),
        [],
        [],
        []
    )

    ch_vcf_tbi_sample = BCFTOOLS_VIEW.out.vcf.join(BCFTOOLS_VIEW.out.tbi)
    ch_bam_bai_vcf_tbi = ch_bam_bai.join(ch_vcf_tbi_sample)

    // Create per-chromosome interval BED files for scatter-gather parallelisation.
    // If ase_intervals is provided: split that BED by chromosome — each job sees only
    // exonic variants on one chromosome, reducing memory and runtime per job.
    // Otherwise: derive whole-chromosome BEDs from the FAI as fallback.
    if (use_intervals) {
        SPLIT_BED_BY_CHROM(ch_ase_intervals)
        ch_interval_files = SPLIT_BED_BY_CHROM.out.intervals.flatten()
        ch_versions = ch_versions.mix( SPLIT_BED_BY_CHROM.out.versions )
    } else {
        SPLIT_FAI_INTERVALS(ch_fai)
        ch_interval_files = SPLIT_FAI_INTERVALS.out.intervals.flatten()
        ch_versions = ch_versions.mix( SPLIT_FAI_INTERVALS.out.versions )
    }

    // Scatter: combine each sample with each interval → one job per sample per chromosome
    ch_scattered = ch_bam_bai_vcf_tbi
        .combine(ch_interval_files)
        .map { meta, bam, bai, vcf, tbi, interval ->
            def new_meta = meta + [interval: interval.baseName]
            [[new_meta, bam, bai, vcf, tbi], interval]
        }
        .multiMap { bam_vcf_input, interval ->
            bam_vcf:   bam_vcf_input
            intervals: interval
        }

    GATK4_ASEREADCOUNTER(
        ch_scattered.bam_vcf,
        ch_fasta,
        ch_fai,
        ch_dict,
        ch_scattered.intervals
    )

    // Gather: group per-interval CSVs back by original sample ID, then merge into one
    ch_ase_csvs = GATK4_ASEREADCOUNTER.out.csv
        .map { meta, csv -> [meta - meta.subMap('interval'), csv] }
        .groupTuple()

    MERGE_ASE_CSVS(ch_ase_csvs)

    BOOTSTRAPANN(
        ch_ind_vcf_tbi
            .map { meta, vcf, tbi -> [meta.id, meta, vcf, tbi] }
            .join(MERGE_ASE_CSVS.out.csv.map { meta, csv -> [meta.id, csv] })
            .map { id, meta, vcf, tbi, csv -> [meta, vcf, tbi, csv] }
    )

    TABIX_BGZIPTABIX(BOOTSTRAPANN.out.vcf)

    TABIX_BGZIPTABIX.out.gz_tbi
        .collect{it[1]}
        .ifEmpty([])
        .toList()
        .set { file_list_vcf }

    TABIX_BGZIPTABIX.out.gz_tbi
        .collect{it[2]}
        .ifEmpty([])
        .toList()
        .set { file_list_tbi }

    ch_case_info
        .combine(file_list_vcf)
        .combine(file_list_tbi)
        .set { ch_vcf_tbi }

    ch_vcf_tbi.branch {
        meta, vcf, tbi ->
            single: vcf.size() == 1
                return [meta, vcf]
            multiple: vcf.size() > 1
                return [meta, vcf, tbi]
        }.set { ch_case_vcf }

    BCFTOOLS_MERGE( ch_case_vcf.multiple,
        ch_fasta,
        ch_fai,
        [[],[]]
    )

    RENAME_FILES( ch_case_vcf.single)

    BCFTOOLS_MERGE.out.vcf
        .mix( RENAME_FILES.out.output )
        .set { ch_vcf_merged }

    TABIX_TABIX( ch_vcf_merged )

    ch_in_split_multi = ch_vcf_merged.join(TABIX_TABIX.out.tbi)
    SPLIT_MULTIALLELICS(ch_in_split_multi, ch_fasta)

    ch_remove_dup_in = SPLIT_MULTIALLELICS.out.vcf.join(SPLIT_MULTIALLELICS.out.tbi)
    REMOVE_DUPLICATES(ch_remove_dup_in, ch_fasta)

    ch_versions = ch_versions.mix( BCFTOOLS_NORM.out.versions.first() )
    ch_versions = ch_versions.mix( BCFTOOLS_VIEW.out.versions.first() )
    ch_versions = ch_versions.mix( GATK4_ASEREADCOUNTER.out.versions.first() )
    ch_versions = ch_versions.mix( MERGE_ASE_CSVS.out.versions.first() )
    ch_versions = ch_versions.mix( BOOTSTRAPANN.out.versions.first() )
    ch_versions = ch_versions.mix( TABIX_BGZIPTABIX.out.versions.first() )
    ch_versions = ch_versions.mix( BCFTOOLS_MERGE.out.versions.first() )
    ch_versions = ch_versions.mix( RENAME_FILES.out.versions.first() )
    ch_versions = ch_versions.mix( TABIX_TABIX.out.versions.first() )
    ch_versions = ch_versions.mix( SPLIT_MULTIALLELICS.out.versions.first() )
    ch_versions = ch_versions.mix( REMOVE_DUPLICATES.out.versions.first() )

    emit:
    vcf      = REMOVE_DUPLICATES.out.vcf // channel: [ val(meta), [ path(vcf) ] ]
    tbi      = REMOVE_DUPLICATES.out.tbi // channel: [ val(meta), [ path(tbi) ] ]
    versions = ch_versions               // channel: [ path(versions.yml) ]
}
