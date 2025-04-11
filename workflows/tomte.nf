/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { FASTQC                          } from '../modules/nf-core/fastqc/main'
include { MULTIQC                         } from '../modules/nf-core/multiqc/main'
include { PEDDY                           } from '../modules/nf-core/peddy/main'
include { CREATE_PEDIGREE_FILE            } from '../modules/local/create_pedigree_file'
include { ESTIMATE_HB_PERC                } from '../modules/local/estimate_hb_perc'
include { paramsSummaryMap                } from 'plugin/nf-schema'
include { paramsSummaryMultiqc            } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML          } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText          } from '../subworkflows/local/utils_nfcore_tomte_pipeline'
include { SAMTOOLS_CONVERT as CRAM_TO_BAM } from '../modules/nf-core/samtools/convert/main'

//
// SUBWORKFLOW: local
//
include { DOWNLOAD_REFERENCES     } from '../subworkflows/local/download_references'
include { PREPARE_REFERENCES      } from '../subworkflows/local/prepare_references/main'
include { ALIGNMENT               } from '../subworkflows/local/alignment/main'
include { BAM_QC                  } from '../subworkflows/local/bam_qc'
include { ANALYSE_TRANSCRIPTS     } from '../subworkflows/local/analyse_transcripts/main'
include { CALL_VARIANTS           } from '../subworkflows/local/call_variants'
include { ALLELE_SPECIFIC_CALLING } from '../subworkflows/local/allele_specific_calling'
include { ANNOTATE_SNV            } from '../subworkflows/local/annotate_snv'
include { IGV_TRACKS              } from '../subworkflows/local/igv_tracks/main'

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

    // Mandatory
    ch_samples        = ch_samplesheet.map { meta, fastqs -> meta }
    ch_case_info      = ch_samples.toList().map { create_case_channel(it) }
    ch_platform       = Channel.from(params.platform).collect()

    // Optional
    ch_vep_refs_download_unprocessed = params.vep_refs_download         ? Channel.fromPath(params.vep_refs_download)
                                                                        : Channel.empty()

    DOWNLOAD_REFERENCES(
        params.genome,
        params.gencode_annotation_version,
        ch_vep_refs_download_unprocessed,
        params.vep_cache_version
    ).set { downloads }
    ch_versions = ch_versions.mix(DOWNLOAD_REFERENCES.out.versions)

    // Optional
    ch_fasta                      = params.fasta                        ? Channel.fromPath(params.fasta).map {it -> [[id:it.getSimpleName()], it]}.collect()
                                                                        : downloads.fasta.map {it -> [[id:it.getSimpleName()], it]}.collect()
    ch_gtf                        = params.gtf                          ? Channel.fromPath(params.gtf).map {it -> [[id:it.getSimpleName()], it]}.collect()
                                                                        : downloads.gtf.map {it -> [[id:it.getSimpleName()], it]}.collect()
    ch_vep_cache_unprocessed      = params.vep_cache                    ? Channel.fromPath(params.vep_cache)
                                                                        : Channel.empty().mix(downloads.vep_cache)
    ch_vep_extra_files_unsplit    = params.vep_plugin_files             ? Channel.fromPath(params.vep_plugin_files)
                                                                        : Channel.empty().mix(downloads.vep_plugin)
    ch_fai                        = params.fai                          ? Channel.fromPath(params.fai).map {it -> [[id:it.getSimpleName()], it]}.collect()
                                                                        : Channel.empty()
    ch_gene_panel_clinical_filter = params.gene_panel_clinical_filter   ? Channel.fromPath(params.gene_panel_clinical_filter).collect()
                                                                        : Channel.empty()
    ch_ref_drop_annot_file        = params.reference_drop_annot_file    ? Channel.fromPath(params.reference_drop_annot_file).collect()
                                                                        : Channel.empty()
    ch_ref_drop_count_file        = params.reference_drop_count_file    ? Channel.fromPath(params.reference_drop_count_file).collect()
                                                                        : Channel.empty()
    ch_ref_drop_splice_folder     = params.reference_drop_splice_folder ? Channel.fromPath(params.reference_drop_splice_folder).collect()
                                                                        : Channel.empty()
    ch_salmon_index               = params.salmon_index                 ? Channel.fromPath(params.salmon_index)
                                                                        : Channel.empty()
    ch_star_index                 = params.star_index                   ? Channel.fromPath(params.star_index).map {it -> [[id:it.getSimpleName()], it]}.collect()
                                                                        : Channel.empty()
    ch_transcript_fasta           = params.transcript_fasta             ? Channel.fromPath(params.transcript_fasta)
                                                                        : Channel.empty()
    ch_sequence_dict              = params.sequence_dict                ? Channel.fromPath(params.sequence_dict).map{ it -> [[id:it.getSimpleName()], it] }.collect()
                                                                        : Channel.empty()
    ch_subsample_bed              = params.subsample_bed                ? Channel.fromPath(params.subsample_bed).collect()
                                                                        : Channel.empty()
    ch_hb_genes                   = params.hb_genes                     ? Channel.fromPath(params.hb_genes).collect()
                                                                        : Channel.empty()

    // Read and store paths in the vep_plugin_files file
    ch_vep_extra_files_unsplit.splitCsv(header: true)
    .flatMap { row ->
        row.vep_files.split(',').collect { file(it.trim()) }
    }
    .map { path ->
        if (params.skip_download_vep) {
            if(path.isFile() || path.isDirectory()) {
                return path
            } else {
                error("\nVep database file ${path} does not exist.")
            }
        } else {
            return path
        }
    }
    .collect()
    .set { ch_vep_extra_files }

    PREPARE_REFERENCES(
        ch_fasta,
        ch_fai,
        ch_star_index,
        ch_gtf,
        ch_vep_cache_unprocessed,
        ch_transcript_fasta,
        ch_salmon_index,
        ch_sequence_dict
    ).set { ch_references }
    ch_versions = ch_versions.mix(PREPARE_REFERENCES.out.versions)

    // Prepare input
    ch_samplesheet
        .branch {
            fastq: it[1].any { it.toString().endsWith('.fastq.gz') || it.toString().endsWith('.fq.gz') }
            cram:  it[1].any { it.toString().endsWith('.cram') }
            bam:   it[1].any { it.toString().endsWith('.bam') }
        }
        .set { ch_input_branch }
    ch_cram_reads = ch_input_branch.cram.map{ meta, cram_crai -> [meta, cram_crai[0], cram_crai[1]] }

    // Convet cram to bam
    CRAM_TO_BAM(ch_cram_reads, ch_references.fasta, ch_references.fai)
    ch_bam_from_cram = CRAM_TO_BAM.out.bam.join(CRAM_TO_BAM.out.bai).map{ meta, bam, bai -> [meta, [bam, bai]]}

    ch_bam_reads = ch_input_branch.bam.mix(ch_bam_from_cram)
    ch_fastq_reads = ch_input_branch.fastq

    FASTQC (
        ch_fastq_reads
    )
    ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect{it[1]})
    ch_versions = ch_versions.mix(FASTQC.out.versions)

    ALIGNMENT(
        ch_fastq_reads,
        ch_bam_reads,
        ch_references.star_index,
        ch_references.gtf,
        ch_platform,
        ch_subsample_bed,
        params.seed_frac,
        params.num_reads,
        params.skip_subsample_region,
        params.skip_downsample,
        ch_references.salmon_index,
        ch_references.fasta
    ).set { ch_alignment }
    ch_versions = ch_versions.mix(ALIGNMENT.out.versions)

    BAM_QC(
        ch_alignment.bam,
        ch_references.fasta,
        ch_references.refflat,
        ch_references.interval_list
    )
    ch_versions = ch_versions.mix(BAM_QC.out.versions)


    ANALYSE_TRANSCRIPTS(
        ch_alignment.bam_bai,
        ch_alignment.bam_ds_bai,
        ch_references.gtf,
        ch_references.fasta_fai,
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
        ch_case_info,
        params.skip_drop_ae,
        params.skip_export_counts_drop
    )
    ch_versions = ch_versions.mix(ANALYSE_TRANSCRIPTS.out.versions)

    CALL_VARIANTS(
        ch_alignment.bam_bai,
        ch_references.fasta,
        ch_references.fai,
        ch_references.sequence_dict,
        params.variant_caller
    )
    ch_versions = ch_versions.mix(CALL_VARIANTS.out.versions)

    ALLELE_SPECIFIC_CALLING(
        CALL_VARIANTS.out.vcf_tbi,
        ch_alignment.bam_bai,
        ch_references.fasta,
        ch_references.fai,
        ch_references.sequence_dict,
        ch_case_info
    )
    ch_versions = ch_versions.mix(ALLELE_SPECIFIC_CALLING.out.versions)

    ANNOTATE_SNV (
        ALLELE_SPECIFIC_CALLING.out.vcf,
        params.genome,
        params.vep_cache_version,
        ch_references.vep_cache,
        ch_references.fasta,
        ch_vep_extra_files,
        ch_gene_panel_clinical_filter
    )
    ch_versions = ch_versions.mix(ANNOTATE_SNV.out.versions)

    IGV_TRACKS(
        ch_alignment.star_wig,
        ch_references.chrom_sizes,
        ch_alignment.spl_junc
    )
    ch_versions = ch_versions.mix(IGV_TRACKS.out.versions)

    ch_pedfile = CREATE_PEDIGREE_FILE(ch_samples.toList()).ped
    ch_versions = ch_versions.mix(CREATE_PEDIGREE_FILE.out.versions)
    PEDDY (
        CALL_VARIANTS.out.vcf_tbi,
        ch_pedfile
    )
    ch_versions = ch_versions.mix(PEDDY.out.versions)

    ESTIMATE_HB_PERC(ALIGNMENT.out.gene_counts, ch_hb_genes)
    ch_versions = ch_versions.mix(ESTIMATE_HB_PERC.out.versions)

    /*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    COLLECT SOFTWARE VERSIONS & MultiQC
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

    softwareVersionsToYAML(ch_versions)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name:  'tomte_'  + 'pipeline_software_' +  'mqc_'  + 'versions.yml',
            sort: true,
            newLine: true
        ).set { ch_collated_versions }


    //
    // MODULE: MultiQC
    //

    summary_params                        = paramsSummaryMap(workflow, parameters_schema: "nextflow_schema.json")
    ch_workflow_summary                   = Channel.value(paramsSummaryMultiqc(summary_params))
    ch_multiqc_config                     = Channel.fromPath("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
    ch_multiqc_custom_config              = params.multiqc_config ? Channel.fromPath(params.multiqc_config, checkIfExists: true) : Channel.empty()
    ch_multiqc_logo                       = params.multiqc_logo ? Channel.fromPath(params.multiqc_logo, checkIfExists: true) : Channel.empty()
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
    ch_multiqc_files                      = ch_multiqc_files.mix(BAM_QC.out.metrics_general_rna.collect{it[1]}.ifEmpty([]))
    ch_multiqc_files                      = ch_multiqc_files.mix(BAM_QC.out.metrics_insert_size.collect{it[1]}.ifEmpty([]))
    ch_multiqc_files                      = ch_multiqc_files.mix(ANALYSE_TRANSCRIPTS.out.stats_gtf.collect{it[1]}.ifEmpty([]))
    ch_multiqc_files                      = ch_multiqc_files.mix(CALL_VARIANTS.out.stats.collect{it[1]}.ifEmpty([]))
    ch_multiqc_files                      = ch_multiqc_files.mix(ANNOTATE_SNV.out.report.collect{it[1]}.ifEmpty([]))

    MULTIQC (
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList(),
        [],
        []
    )

    emit:
    multiqc_report         = MULTIQC.out.report.toList() // channel: /path/to/multiqc_report.html
    vcf_tbi                = CALL_VARIANTS.out.vcf_tbi   // channel: [ val(meta), path(vcf), path(tbi) ]
    junction_bed           = IGV_TRACKS.out.bed          // channel: [ val(meta), path(bed.gz), path(tbi) ]
    drop_ae_out_clinical   = ANALYSE_TRANSCRIPTS.out.drop_ae_out_clinical // channel: [ path(drop_AE_clinical.tsv) ]
    drop_ae_out_research   = ANALYSE_TRANSCRIPTS.out.drop_ae_out_research // channel: [ path(drop_AE_research.tsv) ]
    drop_as_out_clinical   = ANALYSE_TRANSCRIPTS.out.drop_as_out_clinical // channel: [ path(drop_AS_clinical.tsv) ]
    drop_as_out_research   = ANALYSE_TRANSCRIPTS.out.drop_as_out_research // channel: [ path(drop_AS_research.tsv) ]
    bigwig                 = IGV_TRACKS.out.bw           // channel: [ val(meta), path(bw) ]
    ped                    = ch_pedfile                  // channel: [ path(ped_file) ]
    multiqc_data           = MULTIQC.out.data            // channel: [ path(multiqc_data) ]
    hb_estimates           = ESTIMATE_HB_PERC.out.json   // channel: [ val(meta), path(json) ]
    bam_bai                = ALIGNMENT.out.bam_bai       // channel: [ val(meta), path(bam), path(bai) ]
    versions               = ch_versions                 // channel: [ path(versions.yml) ]
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
