//
// ANALYSE TRANSCRITPS
//

include { STRINGTIE_STRINGTIE               } from '../../../modules/nf-core/stringtie/stringtie'
include { GFFCOMPARE                        } from '../../../modules/nf-core/gffcompare'
include { DROP_SAMPLE_ANNOT                 } from '../../../modules/local/drop/drop_sample_annot'
include { DROP_CONFIG_RUN_AE                } from '../../../modules/local/drop/drop_config_runAE'
include { DROP_CONFIG_RUN_AS                } from '../../../modules/local/drop/drop_config_runAS'
include { DROP_CONFIG_RUN_MAE                } from '../../../modules/local/drop/drop_config_runMAE'
include { DROP_FILTER_RESULTS               } from '../../../modules/local/drop/drop_filter_results'
include { DROP_PUT_TOGETHER_EXPORTED_COUNTS } from '../../../modules/local/drop/drop_put_together_exported_counts'

workflow ANALYSE_TRANSCRIPTS {
    take:
    ch_bam_bai                    //   channel: [mandatory] [ val(meta), [ path(bam) , path(bai) ] ]
    ch_bam_ds_bai                 //   channel: [mandatory] [ val(meta), [ path(bam) , path(bai) ] ]
    ch_gtf                        //   channel: [mandatory] [ val(meta), [ path(gtf) ] ]
    ch_fasta_fai                  //   channel: [mandatory] [ val(meta), [ path(fasta), path(fai) ]
    ch_dict                       //   channel:   [mandatory] [ val(meta), path(dict) ]
    ch_ref_drop_count_file        //   channel: [optional]  [ path(tsv) ]
    ch_ref_drop_annot_file        //   channel: [optional]  [ path(tsv) ]
    ch_ref_drop_splice_folder     //   channel: [optional]  [ path(folder) ]
    genome                        // parameter: [optional]  'hg19', 'GRCh37', 'hg38' or 'GRCh38'
    drop_group_samples_ae         // parameter: [optional]  default: 'outrider'
    drop_group_samples_as         // parameter: [optional]  default: 'fraser'
    drop_padjcutoff_ae            // parameter: [optional]  default: 0.05
    drop_padjcutoff_as            // parameter: [optional]  default: 0.1
    drop_zscorecutoff             // parameter: [optional]  default: 0
    ch_gene_panel_clinical_filter //   channel: [optional]  [ path(tsv) ]
    case_info                     //   channel: [optional]  [ val(case_id) ]
    skip_drop_ae                  // parameter: [mandatory] default: false
    skip_drop_as                  // parameter: [mandatory] default: false
    skip_drop_mae                 // parameter: [mandatory] calculated from samplesheet
    skip_export_counts_drop       // parameter: [mandatory] default: true
    ch_drop_mae_high_q_vcf_tbi    //   channel: [ path(ch_high_q_vcf), path(ch_high_q_vcf_tbi) ]
    skip_stringtie                // parameter: [mandatory] default: false

    main:
    ch_versions = Channel.empty()

    log.info("[MAE-DEBUG] Hello world")
    log.info "[MAE-DEBUG] skip_drop_ae=${skip_drop_ae}  skip_drop_as=${skip_drop_as}  skip_drop_mae=${skip_drop_mae}"
    log.info "[MAE-DEBUG] genome=${genome}"

    ch_bam_ds_bai
        .map { meta, bam, bai ->
        [ meta.id, meta.single_end, meta.strandedness, meta.sex, meta.vcf, meta.vcf_tbi, meta.dna_id, bam, bai ]
        }
        .collect(flat:false)
        .map { it.sort { a, b -> a[0] <=> b[0] } } // Sort on ID
        .map { it.transpose() }
        .set { ch_bam_files_annot }

    ch_bam_files_annot
        .map { _id, _single_end, _strandedness, _sex, _vcf, _vcf_tbi, _dna_id, bam, bai ->
            [ bam, bai ]
        }
        .set{ ch_bam_bai_files }

    // Make channel containing only vcf and its index, removing NAs in case those files are present only for some samples
    ch_bam_ds_bai
        .filter { meta, bam, bai ->
            meta.vcf && meta.vcf != "NA" && meta.vcf.toString() != "" && meta.vcf.toString().trim() != ""
        }
        .map { meta, bam, bai ->
            [ meta.vcf, meta.vcf_tbi ]
        }
        .set { ch_vcf_tbi_files }

    // DROP
    if ( !skip_drop_ae | !skip_drop_as | !skip_drop_mae ) {
        // Generates count files for samples and merges them with reference count file
        DROP_SAMPLE_ANNOT(
            ch_bam_files_annot,
            ch_ref_drop_count_file.ifEmpty([]),
            ch_ref_drop_annot_file.ifEmpty([]),
            drop_group_samples_ae,
            drop_group_samples_as
        )

        ch_bam_files_annot.view { x -> "[MAE-DEBUG] ch_bam_files_annot item: ${x[0]} (ids)  -> bam files: ${x[7]?.size() ?: 'NA'}" }
        ch_bam_bai_files.view   { x -> "[MAE-DEBUG] ch_bam_bai_files item: ${x.collect{ it -> it.getName() }}" }
        ch_bam_ds_bai.view { meta, _bam, _bai ->
            "[MAE-DEBUG] sample=${meta.id}  vcf=${meta.vcf ?: 'null'}  tbi=${meta.vcf_tbi ?: 'null'}"
        }
        ch_vcf_tbi_files.view { vcf, tbi -> "[MAE-DEBUG] ch_vcf_tbi_files: vcf=${vcf}  tbi=${tbi}" }
        ch_fasta_fai.view    { meta, fasta, fai -> "[MAE-DEBUG] ch_fasta_fai: fasta=${fasta} fai=${fai}" }
        ch_gtf.view          { meta, gtf -> "[MAE-DEBUG] ch_gtf: gtf=${gtf}" }
        ch_dict.view         { meta, dict -> "[MAE-DEBUG] ch_dict: dict=${dict}" }
        ch_drop_mae_high_q_vcf_tbi.view { v, t ->
            "[MAE-DEBUG] ch_drop_mae_high_q_vcf_tbi: vcf=${v} tbi=${t}"
        }


        // Generates config file and runs Aberrant expression module
        if ( !skip_drop_ae ) {

            log.info "[MAE-DEBUG] Entering MAE block: invoking DROP_CONFIG_RUN_MAE"

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
            ch_versions = ch_versions.mix( DROP_CONFIG_RUN_AE.out.versions )
        }

        // Generates config file and runs Aberrant splicing module
        if ( !skip_drop_as ) {
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
            ch_versions = ch_versions.mix( DROP_CONFIG_RUN_AS.out.versions )
        }

        // Generates config file and runs monoallelic expression module
        if ( !skip_drop_mae ) {
            DROP_CONFIG_RUN_MAE(
                ch_fasta_fai,
                ch_gtf,
                ch_dict,
                DROP_SAMPLE_ANNOT.out.drop_annot,
                genome,
                ch_bam_bai_files,
                ch_vcf_tbi_files,
                ch_drop_mae_high_q_vcf_tbi
            )
            ch_versions = ch_versions.mix( DROP_CONFIG_RUN_MAE.out.versions )
        }

        ch_out_drop_gene_name = !skip_drop_ae ? DROP_CONFIG_RUN_AE.out.drop_gene_name.collect()
            : !skip_drop_as  ? DROP_CONFIG_RUN_AS.out.drop_gene_name.collect()
            : DROP_CONFIG_RUN_MAE.out.drop_gene_name.collect()

        ch_out_drop_ae_rds    = !skip_drop_ae ? DROP_CONFIG_RUN_AE.out.drop_ae_rds.collect()
            : Channel.empty()
        ch_out_drop_as_tsv    = !skip_drop_as ? DROP_CONFIG_RUN_AS.out.drop_as_tsv.collect()
            : Channel.empty()
        ch_out_drop_mae_tsv   = !skip_drop_mae ? DROP_CONFIG_RUN_MAE.out.drop_mae_tsv.collect()
            : Channel.empty()

        DROP_FILTER_RESULTS(
            case_info,
            ch_gene_panel_clinical_filter.ifEmpty([]),
            ch_out_drop_ae_rds.ifEmpty([]),
            ch_out_drop_gene_name,
            ch_out_drop_as_tsv.ifEmpty([]),
            ch_out_drop_mae_tsv.ifEmpty([])
        )

        ch_versions = ch_versions.mix( DROP_SAMPLE_ANNOT.out.versions )
        ch_versions = ch_versions.mix( DROP_FILTER_RESULTS.out.versions )

        // Generates a folder with exported_counts if required
        if ( !skip_export_counts_drop ) {
            DROP_PUT_TOGETHER_EXPORTED_COUNTS(
                DROP_CONFIG_RUN_AE.out.gene_counts_ae.ifEmpty([]),
                DROP_CONFIG_RUN_AS.out.gene_counts_as.ifEmpty([]),
                ch_gtf
            )
            ch_versions = ch_versions.mix(DROP_PUT_TOGETHER_EXPORTED_COUNTS.out.versions)
        }
    }

    // Stringtie
    ch_bam = ch_bam_bai.map{ meta, bam, _bai -> [meta, [bam]] }
    if ( !skip_stringtie ){
        STRINGTIE_STRINGTIE(
            ch_bam,
            ch_gtf.map{ _meta, gtf -> gtf }
        )

        // Compare stringtie results to reference
        GFFCOMPARE(
            STRINGTIE_STRINGTIE.out.transcript_gtf,
            ch_fasta_fai,
            ch_gtf
        )
        ch_versions = ch_versions.mix(STRINGTIE_STRINGTIE.out.versions.first())
        ch_versions = ch_versions.mix(GFFCOMPARE.out.versions.first())
    }

    emit:
    transcript_gtf        = !skip_stringtie                  // channel: [ val(meta), [ path(transctript_gtf)] ]
        ? STRINGTIE_STRINGTIE.out.transcript_gtf
        : Channel.empty()
    abundance            = !skip_stringtie                   // channel: [ val(meta), [ path(abundance) ] ]
        ? STRINGTIE_STRINGTIE.out.abundance
        : Channel.empty()
    coverage_gtf         = !skip_stringtie                   // channel: [ val(meta), [ path(coverage_gtf) ] ]
        ? STRINGTIE_STRINGTIE.out.coverage_gtf
        : Channel.empty()
    annotated_gtf        = !skip_stringtie                   // channel: [ val(meta), [ path(annotated_gtf) ] ]
        ? GFFCOMPARE.out.annotated_gtf
        : Channel.empty()
    stats_gtf            = !skip_stringtie                   // channel: [ val(meta), [ path(stats) ] ]
        ? GFFCOMPARE.out.stats
        : Channel.empty()
    annotation_drop      = ( !skip_drop_ae | !skip_drop_as ) // channel: [ path(sample_annotation.tsv) ]
        ? DROP_SAMPLE_ANNOT.out.drop_annot
        : Channel.empty()
    config_drop_ae       = !skip_drop_ae                     // channel: [ path(confg_file.yml) ]
        ? DROP_CONFIG_RUN_AE.out.config_drop
        : Channel.empty()
    drop_ae_out          = !skip_drop_ae                     // channel: [ path(drop_output_AE) ]
        ? DROP_CONFIG_RUN_AE.out.drop_ae_out
        : Channel.empty()
    config_drop_as       = !skip_drop_as                     // channel: [ path(confg_file.yml) ]
        ? DROP_CONFIG_RUN_AS.out.config_drop
        : Channel.empty()
    drop_as_out          = !skip_drop_as                     // channel: [ path(drop_output_AS) ]
        ? DROP_CONFIG_RUN_AS.out.drop_as_out
        : Channel.empty()
    drop_ae_out_clinical = !skip_drop_ae                     // channel: [ path(drop_AE_clinical.tsv) ]
        ? DROP_FILTER_RESULTS.out.ae_out_clinical
        : Channel.empty()
    drop_ae_out_research = !skip_drop_ae                     // channel: [ path(drop_AE_research.tsv) ]
        ? DROP_FILTER_RESULTS.out.ae_out_research
        : Channel.empty()
    drop_as_out_clinical = !skip_drop_as                     // channel: [ path(drop_AS_clinical.tsv) ]
        ? DROP_FILTER_RESULTS.out.as_out_clinical
        : Channel.empty()
    drop_as_out_research = !skip_drop_as                     // channel: [ path(drop_AS_research.tsv) ]
        ? DROP_FILTER_RESULTS.out.as_out_research
        : Channel.empty()
    versions             = ch_versions                       // channel: [ path(versions.yml) ]
}
