/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// SUBWORKFLOW: nf-core
//
include { paramsSummaryMultiqc   } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'

//
// SUBWORKFLOW: local
//
include { PREPARE_REFERENCES      } from '../subworkflows/local/prepare_references'
include { ALIGNMENT               } from '../subworkflows/local/alignment'
include { BAM_QC                  } from '../subworkflows/local/bam_qc'
include { ANALYSE_TRANSCRIPTS     } from '../subworkflows/local/analyse_transcripts'
include { CALL_VARIANTS           } from '../subworkflows/local/call_variants'
include { ALLELE_SPECIFIC_CALLING } from '../subworkflows/local/allele_specific_calling'
include { ANNOTATE_SNV            } from '../subworkflows/local/annotate_snv'
include { IGV_TRACKS              } from '../subworkflows/local/igv_tracks'
include { methodsDescriptionText  } from '../subworkflows/local/utils_nfcore_tomte_pipeline'

//
// MODULE: nf-core
//

//
// MODULE: local
//
include { FASTQC                 } from '../modules/nf-core/fastqc/main'
include { MULTIQC                } from '../modules/nf-core/multiqc/main'
include { paramsSummaryMap       } from 'plugin/nf-validation'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow TOMTE {

    take:
    ch_samplesheet // channel: samplesheet read in from --input

    main:

    ch_versions      = Channel.empty()
    ch_multiqc_files = Channel.empty()

    ch_samples   = ch_samplesheet.map { meta, fastqs -> meta }
    ch_case_info = ch_samples.toList().map { create_case_channel(it) }

    ch_vep_cache_unprocessed      = params.vep_cache                    ? Channel.fromPath(params.vep_cache).map { it -> [[id:'vep_cache'], it] }.collect()
                                                                        : Channel.value([[],[]])
    ch_vep_filters                = params.vep_filters                  ? Channel.fromPath(params.vep_filters).collect()
                                                                        : Channel.value([])
    fai                           = params.fai                          ? Channel.fromPath(params.fai).map {it -> [[id:it[0].simpleName], it]}.collect()
                                                                        : Channel.empty()
    ch_ref_drop_count_file        = params.reference_drop_count_file    ? Channel.fromPath(params.reference_drop_count_file).collect()
                                                                        : Channel.empty()
    ch_ref_drop_splice_folder     = params.reference_drop_splice_folder ? Channel.fromPath(params.reference_drop_splice_folder).collect()
                                                                        : Channel.empty()
    ch_ref_drop_annot_file        = params.reference_drop_annot_file    ? Channel.fromPath(params.reference_drop_annot_file).collect()
                                                                        : Channel.empty()
    ch_gene_panel_clinical_filter = params.gene_panel_clinical_filter   ? Channel.fromPath(params.gene_panel_clinical_filter).collect()
                                                                        : Channel.empty()
    ch_vep_extra_files_unsplit    = params.vep_plugin_files             ? Channel.fromPath(params.vep_plugin_files).collect()
                                                                        : Channel.value([])
    ch_platform                   = params.platform.collect()

    // Read and store paths in the vep_plugin_files file
    ch_vep_extra_files_unsplit.splitCsv ( header:true )
        .map { row ->
            f = file(row.vep_files[0])
            if(f.isFile() || f.isDirectory()){
                return [f]
            } else {
                error("\nVep database file ${f} does not exist.")
            }
        }
        .collect()
        .set {ch_vep_extra_files}

    PREPARE_REFERENCES(
        params.fasta,
        fai,
        params.star_index,
        params.gtf,
        ch_vep_cache_unprocessed,
        params.transcript_fasta,
        params.salmon_index
    ).set { ch_references }

    // Gather built indices or get them from the params
    ch_chrom_sizes      = ch_references.chrom_sizes
    ch_sequence_dict    = params.sequence_dict          ? Channel.fromPath(params.sequence_dict).collect()
                                                        : ( ch_references.sequence_dict            ?: Channel.empty() )
    ch_subsample_bed    = params.subsample_bed          ? Channel.fromPath(params.subsample_bed).collect()
                                                        : Channel.empty()
    ch_vep_cache        = ( params.vep_cache && params.vep_cache.endsWith("tar.gz") )  ? ch_references.vep_resources
                                                        : ( params.vep_cache  ? Channel.fromPath(params.vep_cache).collect() : Channel.value([]) )

    FASTQC (
        ch_samplesheet
    )
    ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect{it[1]})
    ch_versions = ch_versions.mix(FASTQC.out.versions.first())

    ALIGNMENT(
        ch_samplesheet,
        ch_references.star_index,
        ch_references.gtf,
        ch_platform,
        params.subsample_bed,
        params.seed_frac,
        params.num_reads,
        params.switch_subsample_region,
        params.switch_downsample,
        ch_references.salmon_index,
        ch_references.fasta_meta
    ).set { ch_alignment }
    ch_versions = ch_versions.mix(ALIGNMENT.out.versions)

    BAM_QC(
        ch_alignment.bam,
        ch_references.fasta_no_meta,
        ch_references.refflat,
        ch_references.interval_list
    )
    ch_versions = ch_versions.mix(BAM_QC.out.versions)

    ANALYSE_TRANSCRIPTS(
        ch_alignment.bam_bai,
        ch_alignment.bam_ds_bai,
        ch_references.gtf,
        ch_references.fasta_fai_meta,
        ch_alignment.gene_counts,
        ch_ref_drop_count_file,
        ch_ref_drop_annot_file,
        ch_ref_drop_splice_folder,
        params.genome,
        params.drop_group_samples_ae,
        params.drop_group_samples_as,
        params.drop_padjcutoff_ae,
        params.drop_padjcutoff_as,
        params.drop_zscorecutoff,
        ch_gene_panel_clinical_filter,
        ch_case_info
    )
    ch_versions = ch_versions.mix(ANALYSE_TRANSCRIPTS.out.versions)

    CALL_VARIANTS(
        ch_alignment.bam_bai,
        ch_references.fasta_no_meta,
        ch_references.fai_no_meta,
        ch_references.sequence_dict,
        params.variant_caller
    )
    ch_versions = ch_versions.mix(CALL_VARIANTS.out.versions)

    ALLELE_SPECIFIC_CALLING(
        CALL_VARIANTS.out.vcf_tbi,
        ch_alignment.bam_bai,
        ch_references.fasta_no_meta,
        ch_references.fai_no_meta,
        ch_references.sequence_dict,
        ch_references.interval_list,
        ch_case_info
    )
    ch_versions = ch_versions.mix(ALLELE_SPECIFIC_CALLING.out.versions)

    ANNOTATE_SNV (
        ALLELE_SPECIFIC_CALLING.out.vcf,
        params.genome,
        params.vep_cache_version,
        ch_vep_cache,
        ch_references.fasta_meta,
        ch_vep_extra_files,
    )
    ch_versions = ch_versions.mix(ANNOTATE_SNV.out.versions)

    IGV_TRACKS(
        ch_alignment.star_wig,
        ch_chrom_sizes,
        ch_alignment.spl_junc
    )
    ch_versions = ch_versions.mix(IGV_TRACKS.out.versions)

    //
    // Collate and save software versions
    //
    softwareVersionsToYAML(ch_versions)
        .collectFile(storeDir: "${params.outdir}/pipeline_info", name: 'tomte_software_mqc_versions.yml', sort: true, newLine: true)
        .set { ch_collated_versions }

    //
    // MODULE: MultiQC
    //

    ch_multiqc_config                     = Channel.fromPath("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
    ch_multiqc_custom_config              = params.multiqc_config ? Channel.fromPath(params.multiqc_config, checkIfExists: true) : Channel.empty()
    ch_multiqc_logo                       = params.multiqc_logo ? Channel.fromPath(params.multiqc_logo, checkIfExists: true) : Channel.empty()
    summary_params                        = paramsSummaryMap(workflow, parameters_schema: "nextflow_schema.json")
    ch_workflow_summary                   = Channel.value(paramsSummaryMultiqc(summary_params))
    ch_multiqc_custom_methods_description = params.multiqc_methods_description ? file(params.multiqc_methods_description, checkIfExists: true) : file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)
    ch_methods_description                = Channel.value(methodsDescriptionText(ch_multiqc_custom_methods_description))
    ch_multiqc_files                      = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_files                      = ch_multiqc_files.mix(ch_collated_versions)
    ch_multiqc_files                      = ch_multiqc_files.mix(ch_methods_description.collectFile(name: 'methods_description_mqc.yaml', sort: false))
    ch_multiqc_files                      = ch_multiqc_files.mix(FASTQC.out.zip.collect{it[1]}.ifEmpty([]))
    ch_multiqc_files                      = ch_multiqc_files.mix(ALIGNMENT.out.fastp_report.collect{it[1]}.ifEmpty([]))
    ch_multiqc_files                      = ch_multiqc_files.mix(ALIGNMENT.out.star_log_final.collect{it[1]}.ifEmpty([]))
    ch_multiqc_files                      = ch_multiqc_files.mix(ALIGNMENT.out.gene_counts.collect{it[1]}.ifEmpty([]))
    ch_multiqc_files                      = ch_multiqc_files.mix(ALIGNMENT.out.salmon_info.collect{it[1]}.ifEmpty([]))
    ch_multiqc_files                      = ch_multiqc_files.mix(BAM_QC.out.metrics.collect{it[1]}.ifEmpty([]))
    ch_multiqc_files                      = ch_multiqc_files.mix(ANALYSE_TRANSCRIPTS.out.stats_gtf.collect{it[1]}.ifEmpty([]))
    ch_multiqc_files                      = ch_multiqc_files.mix(CALL_VARIANTS.out.stats.collect{it[1]}.ifEmpty([]))
    ch_multiqc_files                      = ch_multiqc_files.mix(ANNOTATE_SNV.out.report.collect{it[1]}.ifEmpty([]))

    MULTIQC (
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList()
    )

    emit:
    multiqc_report = MULTIQC.out.report.toList() // channel: /path/to/multiqc_report.html
    multiqc_report = Channel.empty()
    versions       = ch_versions                 // channel: [ path(versions.yml) ]
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Function to get a list of metadata (e.g. case id) for the case [ meta ]
def create_case_channel(List rows) {
    def case_info = [:]
    def probands = []

    for (item in rows) {
        probands.add(item.sample)
    }

    case_info.probands = probands.unique()
    case_info.id = rows[0].case

    return case_info
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
