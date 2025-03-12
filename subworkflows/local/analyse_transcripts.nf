//
// ANALYSE TRANSCRITPS
//

include { STRINGTIE_STRINGTIE               } from '../../modules/nf-core/stringtie/stringtie/main'
include { GFFCOMPARE                        } from '../../modules/nf-core/gffcompare/main'
include { DROP_SAMPLE_ANNOT                 } from '../../modules/local/drop_sample_annot'
include { DROP_CONFIG_RUN_AE                } from '../../modules/local/drop_config_runAE'
include { DROP_CONFIG_RUN_AS                } from '../../modules/local/drop_config_runAS'
include { DROP_FILTER_RESULTS               } from '../../modules/local/drop_filter_results'
include { DROP_PUT_TOGETHER_EXPORTED_COUNTS } from '../../modules/local/drop_put_together_exported_couts.nf'

workflow ANALYSE_TRANSCRIPTS {
    take:
        ch_bam_bai                    // channel:   [mandatory] [ val(meta), [ path(bam) ],[ path(bai) ] ]
        ch_bam_ds_bai                 // channel:   [mandatory] [ val(meta), [ path(bam) ],[ path(bai) ] ]
        ch_gtf                        // channel:   [mandatory] [ val(meta), [ path(gtf) ] ]
        ch_fasta_fai                  // channel:   [mandatory] [ val(meta), [ path(fasta), path(fai) ]
        gene_counts                   // channel:   [optional]  [ val(meta), path(tsv) ]
        ch_ref_drop_count_file        // channel:   [optional]  [ path(tsv) ]
        ch_ref_drop_annot_file        // channel:   [optional]  [ path(tsv) ]
        ch_ref_drop_splice_folder     // channel:   [optional]  [ path(folder) ]
        genome                        // parameter: [optional]  'hg19', 'GRCh37', 'hg38' or 'GRCh38'
        drop_group_samples_ae         // parameter: [optional]  default: 'outrider'
        drop_group_samples_as         // parameter: [optional]  default: 'fraser'
        drop_padjcutoff_ae            // parameter: [optional]  default: 0.05
        drop_padjcutoff_as            // parameter: [optional]  default: 0.1
        drop_zscorecutoff             // parameter: [optional]  default: 0
        ch_gene_panel_clinical_filter // channel:   [optional]  [ path(tsv) ]
        case_info                     // channel:   [optional]  [ val(case_id) ]
        skip_drop_ae                  // parameter: [mandatory] default: 'false'
        skip_export_counts_drop       // parameter: [mandatory] default: 'true'

    main:
        ch_versions = Channel.empty()

        // DROP
        // Generates count files for samples and merges them with reference count file

        // Generates sample annotation
        ch_bam_files = ch_bam_ds_bai.collect{it[1]}

        ch_bam_ds_bai
            .map { meta, bam, bai ->
            [ meta.id, meta.single_end, meta.strandedness, meta.sex, bam, bai ]
            }
            .collect(flat:false)
            .map { it.sort { a, b -> a[0] <=> b[0] } } // Sort on ID
            .map { it.transpose() }
        .set { ch_bam_files_annot }

        DROP_SAMPLE_ANNOT(
            ch_bam_files_annot,
            ch_ref_drop_count_file.ifEmpty([]),
            ch_ref_drop_annot_file.ifEmpty([]),
            drop_group_samples_ae,
            drop_group_samples_as
        )


        // Generates config file and runs Aberrant expression module

        // Sort bam and bai files for stable order
        // ch_bam_files_sorted = ch_bam_files.toList().map { list ->
        //     list.sort { a, b -> a.getName() <=> b.getName() }
        // }
        // ch_bai_files_sorted = ch_bam_ds_bai.collect{ it[2] }.toList().map { list ->
        //     list.sort { a, b -> a.getName() <=> b.getName() }
        // }
        // ch_bam_bai_files = ch_bam_files_sorted.combine(ch_bai_files_sorted)
        ch_bam_files_annot
            .map { _id, _single_end, _strandedness, _sex, bam, bai ->
                [ bam, bai ] 
            }
        .set{ ch_bam_bai_files }

        DROP_CONFIG_RUN_AE(
            ch_fasta_fai,
            ch_gtf,
            DROP_SAMPLE_ANNOT.out.drop_annot,
            ch_bam_bai_files,
            ch_ref_drop_count_file.ifEmpty([]),
            ch_ref_drop_splice_folder.ifEmpty([]),
            genome,
            drop_group_samples_ae,
            drop_group_samples_as,
            drop_padjcutoff_ae,
            drop_zscorecutoff,
            skip_export_counts_drop
        )

        // Generates config file and runs Aberrant splicing module
        DROP_CONFIG_RUN_AS(
            ch_fasta_fai,
            ch_gtf,
            DROP_SAMPLE_ANNOT.out.drop_annot,
            ch_bam_bai_files,
            ch_ref_drop_count_file.ifEmpty([]),
            ch_ref_drop_splice_folder.ifEmpty([]),
            genome,
            drop_group_samples_as,
            drop_group_samples_ae,
            drop_padjcutoff_as,
            skip_export_counts_drop
        )

        // Generates a folder with exported_counts if required
        DROP_PUT_TOGETHER_EXPORTED_COUNTS(
            DROP_CONFIG_RUN_AE.out.gene_counts_ae.ifEmpty([]),
            DROP_CONFIG_RUN_AS.out.gene_counts_as.ifEmpty([]),
            ch_gtf
        )

        ch_out_drop_ae_rds       = DROP_CONFIG_RUN_AE.out.drop_ae_rds       ? DROP_CONFIG_RUN_AE.out.drop_ae_rds.collect()
                                                                            : Channel.empty()
        ch_out_drop_gene_name_ae = DROP_CONFIG_RUN_AE.out.drop_gene_name    ? DROP_CONFIG_RUN_AE.out.drop_gene_name.collect()
                                                                            : Channel.empty()
        ch_out_drop_gene_name_as = DROP_CONFIG_RUN_AS.out.drop_gene_name    ? DROP_CONFIG_RUN_AS.out.drop_gene_name.collect()
                                                                            : Channel.empty()
        ch_out_drop_as_tsv       = DROP_CONFIG_RUN_AS.out.drop_as_tsv       ? DROP_CONFIG_RUN_AS.out.drop_as_tsv.collect()
                                                                            : Channel.empty()
        ch_out_drop_gene_name    = (!skip_drop_ae) ? ch_out_drop_gene_name_ae : ch_out_drop_gene_name_as

        DROP_FILTER_RESULTS(
            case_info,
            ch_gene_panel_clinical_filter.ifEmpty([]),
            ch_out_drop_ae_rds.ifEmpty([]),
            ch_out_drop_gene_name,
            ch_out_drop_as_tsv.ifEmpty([])
        )

        // Stringtie
        ch_bam = ch_bam_bai.map{ meta, bam, bai -> [meta, [bam]] }
        STRINGTIE_STRINGTIE(
            ch_bam,
            ch_gtf.map{ meta, gtf -> gtf }
        )

        // Compare stringtie results to reference
        GFFCOMPARE(
            STRINGTIE_STRINGTIE.out.transcript_gtf,
            ch_fasta_fai,
            ch_gtf
        )

        ch_versions = ch_versions.mix(DROP_SAMPLE_ANNOT.out.versions)
        ch_versions = ch_versions.mix(DROP_CONFIG_RUN_AE.out.versions)
        ch_versions = ch_versions.mix(DROP_CONFIG_RUN_AS.out.versions)
        ch_versions = ch_versions.mix(DROP_FILTER_RESULTS.out.versions)
        ch_versions = ch_versions.mix(STRINGTIE_STRINGTIE.out.versions.first())
        ch_versions = ch_versions.mix(GFFCOMPARE.out.versions.first())

    emit:
        transcript_gtf        = STRINGTIE_STRINGTIE.out.transcript_gtf   // channel: [ val(meta), [ path(transctript_gtf)] ]
        abundance             = STRINGTIE_STRINGTIE.out.abundance        // channel: [ val(meta), [ path(abundance) ] ]
        coverage_gtf          = STRINGTIE_STRINGTIE.out.coverage_gtf     // channel: [ val(meta), [ path(coverage_gtf) ] ]
        annotated_gtf         = GFFCOMPARE.out.annotated_gtf             // channel: [ val(meta), [ path(annotated_gtf) ] ]
        stats_gtf             = GFFCOMPARE.out.stats                     // channel: [ val(meta), [ path(stats) ] ]
        annotation_drop       = DROP_SAMPLE_ANNOT.out.drop_annot         // channel: [ path(sample_annotation.tsv) ]
        config_drop_ae        = DROP_CONFIG_RUN_AE.out.config_drop       // channel: [ path(confg_file.yml) ]
        drop_ae_out           = DROP_CONFIG_RUN_AE.out.drop_ae_out       // channel: [ path(drop_output_AE) ]
        config_drop_as        = DROP_CONFIG_RUN_AS.out.config_drop       // channel: [ path(confg_file.yml) ]
        drop_as_out           = DROP_CONFIG_RUN_AS.out.drop_as_out       // channel: [ path(drop_output_AS) ]
        drop_ae_out_clinical  = DROP_FILTER_RESULTS.out.ae_out_clinical  // channel: [ path(drop_AE_clinical.tsv) ]
        drop_ae_out_research  = DROP_FILTER_RESULTS.out.ae_out_research  // channel: [ path(drop_AE_research.tsv) ]
        drop_as_out_clinical  = DROP_FILTER_RESULTS.out.as_out_clinical  // channel: [ path(drop_AS_clinical.tsv) ]
        drop_as_out_research  = DROP_FILTER_RESULTS.out.as_out_research  // channel: [ path(drop_AS_research.tsv) ]
        versions              = ch_versions                              // channel: [ path(versions.yml) ]
}
