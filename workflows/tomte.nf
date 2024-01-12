/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    PRINT PARAMS SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { paramsSummaryLog; paramsSummaryMap; fromSamplesheet } from 'plugin/nf-validation'

def logo = NfcoreTemplate.gmsLogo(workflow, params.monochrome_logs)
def citation = '\n' + WorkflowMain.citation(workflow) + '\n'
def summary_params = paramsSummaryMap(workflow)

// Print parameter summary log to screen
log.info logo + paramsSummaryLog(workflow) + citation

WorkflowTomte.initialise(params, log)

// TODO nf-core: change this to use the plugin nf-validation
// Check input path parameters to see if they exist
def checkPathParamList = [
    params.input,
    params.multiqc_config,
    params.fasta,
    params.fai,
    params.sequence_dict,
    params.star_index,
    params.salmon_index,
    params.transcript_fasta,
    params.gtf,
    params.reference_drop_annot_file,
    params.reference_drop_splice_folder,
    params.reference_drop_count_file,
    params.gene_panel_clinical_filter,
    params.subsample_bed,
    params.vep_filters,
    params.vep_cache
]

for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CONFIG FILES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

ch_multiqc_config          = Channel.fromPath("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
ch_multiqc_custom_config   = params.multiqc_config ? Channel.fromPath( params.multiqc_config, checkIfExists: true ) : Channel.empty()
ch_multiqc_logo            = params.multiqc_logo   ? Channel.fromPath( params.multiqc_logo, checkIfExists: true ) : Channel.empty()
ch_multiqc_custom_methods_description = params.multiqc_methods_description ? file(params.multiqc_methods_description, checkIfExists: true) : file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

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

//
// MODULE: local
//

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// SUBWORKFLOW: nf-core/subworkflows
//

//
// MODULE: Installed directly from nf-core/modules
//
include { FASTQC                      } from '../modules/nf-core/fastqc/main'
include { MULTIQC                     } from '../modules/nf-core/multiqc/main'
include { CUSTOM_DUMPSOFTWAREVERSIONS } from '../modules/nf-core/custom/dumpsoftwareversions/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Info required for completion email and summary
def multiqc_report = []

workflow TOMTE {

    ch_versions = Channel.empty()

    // Initialize input channels
    ch_input = Channel.fromPath(params.input)
    Channel.fromSamplesheet("input")
        .tap { ch_original_input }
        .map { meta, fastq1, fastq2 -> meta.id }
        .reduce([:]) { counts, sample -> //get counts of each sample in the samplesheet - for groupTuple
            counts[sample] = (counts[sample] ?: 0) + 1
            counts
        }
        .combine( ch_original_input )
        .map { counts, meta, fastq1, fastq2 ->
            new_meta = meta + [num_lanes:counts[meta.id],
                        read_group:"\'@RG\\tID:"+ fastq1.toString().split('/')[-1] + "\\tPL:ILLUMINA\\tSM:"+meta.id+"\'"]
            if (!fastq2) {
                return [ new_meta + [ single_end:true ], [ fastq1 ] ]
            } else {
                return [ new_meta + [ single_end:false ], [ fastq1, fastq2 ] ]
            }
        }
        .tap{ ch_input_counts }
        .map { meta, fastqs -> fastqs }
        .reduce([:]) { counts, fastqs -> //get line number for each row to construct unique sample ids
            counts[fastqs] = counts.size() + 1
            return counts
        }
        .combine( ch_input_counts )
         .map { lineno, meta, fastqs -> //append line number to sampleid
             new_meta = meta + [id:meta.id+"_T"+lineno[fastqs]]
             return [ new_meta, fastqs ]
         }
         .set { ch_reads }

    ch_samples   = ch_reads.map { meta, fastqs -> meta}
    ch_case_info = ch_samples.toList().map { create_case_channel(it) }

    //if (params.input) {
    //    ch_input = Channel.fromPath(params.input)
    //    CHECK_INPUT (ch_input)
    //    ch_versions = ch_versions.mix(CHECK_INPUT.out.versions)
    //} else {
    //    exit 1, 'Input samplesheet not specified!'
    //}


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

    PREPARE_REFERENCES(
        params.fasta,
        fai,
        params.star_index,
        params.gtf,
        ch_vep_cache_unprocessed,
        params.transcript_fasta,
        params.salmon_index
    ).set { ch_references }
    ch_versions = ch_versions.mix(PREPARE_REFERENCES.out.versions)

    // Gather built indices or get them from the params
    ch_chrom_sizes           = ch_references.chrom_sizes
    ch_sequence_dict         = params.sequence_dict           ? Channel.fromPath(params.sequence_dict).collect()
                                                              : ( ch_references.sequence_dict            ?: Channel.empty() )
    ch_subsample_bed         = params.subsample_bed           ? Channel.fromPath(params.subsample_bed).collect()
                                                              : Channel.empty()
    ch_vep_cache             = ( params.vep_cache && params.vep_cache.endsWith("tar.gz") )  ? ch_references.vep_resources
                                                                : ( params.vep_cache  ? Channel.fromPath(params.vep_cache).collect() : Channel.value([]) )
    //
    // SUBWORKFLOW: Read in samplesheet, validate and stage input files
    //
    // INPUT_CHECK (
    //    file(params.input)
    // )
    // ch_versions = ch_versions.mix(INPUT_CHECK.out.versions)
    // TODO: OPTIONAL, you can use nf-validation plugin to create an input channel from the samplesheet with Channel.fromSamplesheet("input")
    // See the documentation https://nextflow-io.github.io/nf-validation/samplesheets/fromSamplesheet/
    // ! There is currently no tooling to help you write a sample sheet schema

    //
    // MODULE: Run FastQC
    //

    FASTQC (
        ch_reads
    )
    ch_versions = ch_versions.mix(FASTQC.out.versions.first())


    // Alignment
    ALIGNMENT(
        ch_reads,
        ch_references.star_index,
        ch_references.gtf,
        params.platform,
        params.subsample_bed,
        params.seed_frac,
        params.num_reads,
        params.switch_subsample_region,
        params.switch_downsample,
        ch_references.salmon_index,
        ch_references.fasta_meta
    ).set {ch_alignment}
    ch_versions = ch_versions.mix(ALIGNMENT.out.versions)

    // BAM QC
    BAM_QC(
        ch_alignment.bam,
        ch_references.fasta_no_meta,
        ch_references.refflat,
        ch_references.interval_list
    )
    ch_versions = ch_versions.mix(BAM_QC.out.versions)

    // Analyse transcripts
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

    // Call variants
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
        ch_references.interval_list
    )
    ch_versions = ch_versions.mix(ALLELE_SPECIFIC_CALLING.out.versions)

    ANNOTATE_SNV (
        ALLELE_SPECIFIC_CALLING.out.vcf,
        params.genome,
        params.vep_cache_version,
        ch_vep_cache,
        ch_references.fasta_no_meta,
    )
    ch_versions = ch_versions.mix(ANNOTATE_SNV.out.versions)

    IGV_TRACKS(
        ch_alignment.star_wig,
        ch_chrom_sizes,
        ch_alignment.spl_junc
    )

    CUSTOM_DUMPSOFTWAREVERSIONS (
        ch_versions.unique().collectFile(name: 'collated_versions.yml')
    )

    //
    // MODULE: MultiQC
    //
    workflow_summary    = WorkflowTomte.paramsSummaryMultiqc(workflow, summary_params)
    ch_workflow_summary = Channel.value(workflow_summary)

    methods_description    = WorkflowTomte.methodsDescriptionText(workflow, ch_multiqc_custom_methods_description, params)
    ch_methods_description = Channel.value(methods_description)

    ch_multiqc_files = Channel.empty()
    ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_files = ch_multiqc_files.mix(ch_methods_description.collectFile(name: 'methods_description_mqc.yaml'))
    ch_multiqc_files = ch_multiqc_files.mix(CUSTOM_DUMPSOFTWAREVERSIONS.out.mqc_yml.collect())
    ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect{it[1]}.ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(ALIGNMENT.out.fastp_report.collect{it[1]}.ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(ALIGNMENT.out.star_log_final.collect{it[1]}.ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(ALIGNMENT.out.gene_counts.collect{it[1]}.ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(ALIGNMENT.out.salmon_info.collect{it[1]}.ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(BAM_QC.out.metrics.collect{it[1]}.ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(ANALYSE_TRANSCRIPTS.out.stats_gtf.collect{it[1]}.ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(CALL_VARIANTS.out.stats.collect{it[1]}.ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(ANNOTATE_SNV.out.report.collect{it[1]}.ifEmpty([]))

    MULTIQC (
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList()
    )
    multiqc_report = MULTIQC.out.report.toList()
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
        probands.add(item.id.split("_T")[0])
    }

    case_info.probands = probands
    case_info.id = rows[0].case_id

    return case_info
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    COMPLETION EMAIL AND SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow.onComplete {
    if (params.email || params.email_on_fail) {
        NfcoreTemplate.email(workflow, params, summary_params, projectDir, log, multiqc_report)
    }
    NfcoreTemplate.dump_parameters(workflow, params)
    NfcoreTemplate.summary(workflow, params, log)
    if (params.hook_url) {
        NfcoreTemplate.IM_notification(workflow, params, summary_params, projectDir, log)
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
