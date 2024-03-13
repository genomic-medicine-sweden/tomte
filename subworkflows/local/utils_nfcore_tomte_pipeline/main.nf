//
// Subworkflow with functionality specific to the genomic-medicine-sweden/tomte pipeline
//

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { UTILS_NFVALIDATION_PLUGIN } from '../../nf-core/utils_nfvalidation_plugin'
include { paramsSummaryMap          } from 'plugin/nf-validation'
include { fromSamplesheet           } from 'plugin/nf-validation'
include { UTILS_NEXTFLOW_PIPELINE   } from '../../nf-core/utils_nextflow_pipeline'
include { completionEmail           } from '../../nf-core/utils_nfcore_pipeline'
include { completionSummary         } from '../../nf-core/utils_nfcore_pipeline'
include { dashedLine                } from '../../nf-core/utils_nfcore_pipeline'
include { logColours                } from '../../nf-core/utils_nfcore_pipeline'
include { getWorkflowVersion        } from '../../nf-core/utils_nfcore_pipeline'
include { imNotification            } from '../../nf-core/utils_nfcore_pipeline'
include { UTILS_NFCORE_PIPELINE     } from '../../nf-core/utils_nfcore_pipeline'
include { workflowCitation          } from '../../nf-core/utils_nfcore_pipeline'

/*
========================================================================================
    SUBWORKFLOW TO INITIALISE PIPELINE
========================================================================================
*/

workflow PIPELINE_INITIALISATION {

    take:
    version           // boolean: Display version and exit
    help              // boolean: Display help text
    validate_params   // boolean: Boolean whether to validate parameters against the schema at runtime
    monochrome_logs   // boolean: Do not use coloured log outputs
    nextflow_cli_args //   array: List of positional nextflow CLI args
    outdir            //  string: The output directory where the results will be saved
    input             //  string: Path to input samplesheet

    main:

    ch_versions = Channel.empty()

    //
    // Print version and exit if required and dump pipeline parameters to JSON file
    //
    UTILS_NEXTFLOW_PIPELINE (
        version,
        true,
        outdir,
        workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1
    )

    //
    // Validate parameters and generate parameter summary to stdout
    //
    pre_help_text = gmsTomteHatLogo2(monochrome_logs)
    post_help_text = '\n' + workflowCitation() + '\n' + dashedLine(monochrome_logs)
    def String workflow_command = "nextflow run ${workflow.manifest.name} -profile <docker/singularity/.../institute> --input samplesheet.csv --outdir <OUTDIR>"
    UTILS_NFVALIDATION_PLUGIN (
        help,
        workflow_command,
        pre_help_text,
        post_help_text,
        validate_params,
        "nextflow_schema.json"
    )

    //
    // Check config provided to the pipeline
    //
    UTILS_NFCORE_PIPELINE (
        nextflow_cli_args
    )
    //
    // Custom validation for pipeline parameters
    //
    validateInputParameters()

    //
    // Create channel from input file provided through params.input
    //

    Channel
    .fromSamplesheet("input")
    .tap { ch_original_input }
    .map { meta, fastq_1, fastq_2 -> meta.id }
    .reduce([:]) { counts, sample -> // get counts of each sample in the samplesheet - for groupTuple
        counts[sample] = (counts[sample] ?: 0) + 1
        counts
    }
    .combine ( ch_original_input )
    .map { counts, meta, fastq_1, fastq_2 ->
        if (!fastq_2) {
            return [ meta + [ single_end:true, fq_pairs:counts[meta.id] ], [ fastq_1 ] ]
        } else {
            return [ meta + [ single_end:false, fq_pairs:counts[meta.id] ], [ fastq_1, fastq_2 ] ]
        }
    }
    .tap { ch_input_counts }
    .map { meta, fastqs -> fastqs }
    .reduce([:]) { counts, fastqs -> // get number of fastq sets in the run - for creating unique ID:s
        counts[fastqs] = counts.size() + 1
        return counts
    }
    .combine( ch_input_counts )
    .map { lineno, meta, fastqs -> // append line number to sample id for unique set ids
        new_meta = meta + [id:meta.id+"_id"+lineno[fastqs]]
        return [ new_meta, fastqs ]
    }
    .tap { ch_samplesheet } // Output, the rest is just for validation
    .map { meta, fastqs ->
        return [ meta.sample, groupKey( meta + [id:meta.sample], meta.fq_pairs ), fastqs ]
    }
    .groupTuple()
    .map {
        validateInputSamplesheet(it)
    }

    emit:
    samplesheet = ch_samplesheet
    versions    = ch_versions
}

/*
========================================================================================
    SUBWORKFLOW FOR PIPELINE COMPLETION
========================================================================================
*/

workflow PIPELINE_COMPLETION {

    take:
    email           //  string: email address
    email_on_fail   //  string: email address sent on pipeline failure
    plaintext_email // boolean: Send plain-text email instead of HTML
    outdir          //    path: Path to output directory where results will be published
    monochrome_logs // boolean: Disable ANSI colour codes in log output
    hook_url        //  string: hook URL for notifications
    multiqc_report  //  string: Path to MultiQC report

    main:

    summary_params = paramsSummaryMap(workflow, parameters_schema: "nextflow_schema.json")

    //
    // Completion email and summary
    //
    workflow.onComplete {
        if (email || email_on_fail) {
            completionEmail(summary_params, email, email_on_fail, plaintext_email, outdir, monochrome_logs, multiqc_report.toList())
        }

        completionSummary(monochrome_logs)

        if (hook_url) {
            imNotification(summary_params, hook_url)
        }
    }
}

/*
========================================================================================
    FUNCTIONS
========================================================================================
*/
//
// Check and validate pipeline parameters
//
def validateInputParameters() {
    genomeExistsError()
}

//
// Validate channels from input samplesheet
//
def validateInputSamplesheet(input) {
    def (metas, fastqs) = input[1..2]

    // Check that multiple runs of the same sample are of the same datatype i.e. single-end / paired-end
    def endedness_ok = metas.collect{ it.single_end }.unique().size == 1
    if (!endedness_ok) {
        error("Please check input samplesheet -> Multiple runs of a sample must be of the same datatype i.e. single-end or paired-end: ${metas[0].id}")
    }

    return [ metas[0], fastqs ]
}
//
// Get attribute from genome config file e.g. fasta
//
def getGenomeAttribute(attribute) {
    if (params.genomes && params.genome && params.genomes.containsKey(params.genome)) {
        if (params.genomes[ params.genome ].containsKey(attribute)) {
            return params.genomes[ params.genome ][ attribute ]
        }
    }
    return null
}

//
// Exit pipeline if incorrect --genome key provided
//
def genomeExistsError() {
    if (params.genomes && params.genome && !params.genomes.containsKey(params.genome)) {
        def error_string = "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n" +
            "  Genome '${params.genome}' not found in any config files provided to the pipeline.\n" +
            "  Currently, the available genome keys are:\n" +
            "  ${params.genomes.keySet().join(", ")}\n" +
            "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
        error(error_string)
    }
}

//
// Generate methods description for MultiQC
//
def toolCitationText() {
    def citation_text = [
            "Tools used in the workflow included:",
            "BCFtools (Danecek et al. 2021),",
            "DROP (Yépez et al. 2021),",
            params.switch_vep ? "EnsemblVEP (McLaren et al. 2016)," : "",
            "fastp (Chen et al. 2018),",
            "FastQC (Andrews 2010),",
            params.switch_drop_as ? "FRASER (Mertes et al 2021)," : "",
            "GATK (McKenna et al. 2010),",
            params.switch_stringtie ? "GFFCompare (Pertea et al. 2020), StringTie (Pertea et al. 2015)," : "",
            "MultiQC (Ewels et al. 2016),",
            params.switch_drop_ae ? "OUTRIDER (Brechtmann et al. 2018)," : "",
            "SAMtools (Danecek et al. 2021),",
            "Salmon (Patro et al. 2017),",
            "STAR (Dobin et al. 2012),",
            params.switch_build_tracks ? "UCSC tools (Kent et al. 2010)" : "",
            "."
        ].join(' ').trim()

    return citation_text
}

def toolBibliographyText() {
    def reference_text = [
            "<li>Andrews S, (2010) FastQC, URL: https://www.bioinformatics.babraham.ac.uk/projects/fastqc/).</li>",
            params.switch_drop_ae ? "<li>Brechtmann F, Mertes C, Matusevičiūtė A, et al. OUTRIDER: A Statistical Method for Detecting Aberrantly Expressed Genes in RNA Sequencing Data. The American Journal of Human Genetics. 12 2018;103:907-917. doi:10.1016/J.AJHG.2018.10.025</li>" : "",
            "<li>Chen S, Zhou Y, Chen Y, Gu J. fastp: an ultra-fast all-in-one FASTQ preprocessor. Bioinformatics (Oxford, England). 9 2018;34:i884-i890. doi:10.1093/BIOINFORMATICS/BTY560</li>",
            "<li>Dale R, Grüning B, Sjödin A, et al. Bioconda: sustainable and comprehensive software distribution for the life sciences. Nature methods. 7 2018;15:475-476. doi:10.1038/S41592-018-0046-7</li>",
            "<li>Danecek P, Bonfield JK, Liddle J, et al. Twelve years of SAMtools and BCFtools. GigaScience. 1 2021;10:1-4. doi:10.1093/GIGASCIENCE/GIAB008</li>",
            "<li>Dobin A, Davis CA, Schlesinger F, et al. STAR: ultrafast universal RNA-seq aligner. Bioinformatics. 10 2012;29:15-21. doi:10.1093/bioinformatics/bts635</li>",
            "<li>Ewels P, Magnusson M, Lundin S, Käller M. MultiQC: summarize analysis results for multiple tools and samples in a single report. Bioinformatics (Oxford, England). 10 2016;32:3047-3048. doi:10.1093/BIOINFORMATICS/BTW354</li>",
            "<li>Ewels PA, Peltzer A, Fillinger S, et al. The nf-core framework for community-curated bioinformatics pipelines. Nature biotechnology. 3 2020;38:276-278. doi:10.1038/S41587-020-0439-X</li>",
            "<li>Kent WJ, Zweig AS, Barber G, Hinrichs AS, Karolchik D. BigWig and BigBed: enabling browsing of large distributed datasets. Bioinformatics. 9 2010;26:2204-2207. doi:10.1093/BIOINFORMATICS/BTQ351</li>",
            "<li>Kurtzer GM, Sochat V, Bauer MW. Singularity: Scientific containers for mobility of compute. PloS one. 5 2017;12. doi:10.1371/JOURNAL.PONE.0177459</li>",
            "<li>Leprevost FDV, Grüning BA, Aflitos SA, et al. BioContainers: an open-source and community-driven framework for software standardization. Bioinformatics (Oxford, England). 8 2017;33:2580-2582. doi:10.1093/BIOINFORMATICS/BTX192</li>",
            "<li>McKenna A, Hanna M, Banks E, et al. The Genome Analysis Toolkit: A MapReduce framework for analyzing next-generation DNA sequencing data. Genome Research. 9 2010;20:1297-1303. doi:10.1101/GR.107524.110</li>",
            params.switch_vep ? "<li>McLaren W, Gil L, Hunt SE, et al. The Ensembl Variant Effect Predictor. Genome biology. 6 2016;17. doi:10.1186/S13059-016-0974-4</li>" : "",
            "<li>MerkelDirk. Docker. Linux Journal. Published online 3 2014. doi:10.5555/2600239.2600241</li>",
            "<li>Mertes C, Scheller IF, Yépez VA, et al. Detection of aberrant splicing events in RNA-seq data using FRASER. Nature Communications 2021 12:1. 1 2021;12:1-13. doi:10.1038/s41467-020-20573-7</li>",
            "<li>Patro R, Duggal G, Love MI, Irizarry RA, Kingsford C. Salmon provides fast and bias-aware quantification of transcript expression. Nature methods. 2017;14:417-419. doi:10.1038/NMETH.4197</li>",
            params.switch_stringtie ? "<li>Pertea M, Pertea G. GFF Utilities: GffRead and GffCompare. F1000Research. 9 2020;9:304. doi:10.12688/F1000RESEARCH.23297.1</li>" : "",
            params.switch_stringtie ? "<li>Pertea M, Pertea GM, Antonescu CM, Chang TC, Mendell JT, Salzberg SL. StringTie enables improved reconstruction of a transcriptome from RNA-seq reads. Nature biotechnology. 2015;33:290-295. doi:10.1038/NBT.3122</li>": ""
        ].join(' ').trim()

    return reference_text
}

def methodsDescriptionText(mqc_methods_yaml) {
    // Convert  to a named map so can be used as with familar NXF ${workflow} variable syntax in the MultiQC YML file
    def meta = [:]
    meta.workflow = workflow.toMap()
    meta["manifest_map"] = workflow.manifest.toMap()

    // Pipeline DOI
    meta["doi_text"] = meta.manifest_map.doi ? "(doi: <a href=\'https://doi.org/${meta.manifest_map.doi}\'>${meta.manifest_map.doi}</a>)" : ""
    meta["nodoi_text"] = meta.manifest_map.doi ? "": "<li>If available, make sure to update the text to include the Zenodo DOI of version of the pipeline used. </li>"

    // Tool references
    meta["tool_citations"] = toolCitationText().replaceAll(", \\.", ".").replaceAll("\\. \\.", ".").replaceAll(", \\.", ".")
    meta["tool_bibliography"] = toolBibliographyText()


    def methods_text = mqc_methods_yaml.text

    def engine =  new groovy.text.SimpleTemplateEngine()
    def description_html = engine.createTemplate(methods_text).make(meta)

    return description_html.toString()
}

    //
    // GMS tomte logo (WIP)
    //
def gmsTomteLogo(monochrome_logs=true) {
    Map colors = logColours(monochrome_logs)
    String.format(
        """\n
        ${dashedLine(monochrome_logs)}
        ${colors.blue}   _____ __  __  _____            _                  _${colors.reset}
        ${colors.blue}  / ____|  \\/  |/ ____|          | |                | |${colors.reset}
        ${colors.blue} | |  __| \\  / | (___    ______  | |_ ___  _ __ ___ | |_ ___${colors.reset}
        ${colors.blue} | | |_ | |\\/| |\\___ \\  |______| | __/ _ \\| '_ ` _ \\| __/ _ \\${colors.reset}
        ${colors.blue} | |__| | |  | |____) |          | || (_) | | | | | | ||  __/${colors.reset}
        ${colors.blue}  \\_____|_|  |_|_____/            \\__\\___/|_| |_| |_|\\__\\___|${colors.reset}

        ${colors.purple}  ${workflow.manifest.name} ${getWorkflowVersion()}${colors.reset}
        ${dashedLine(monochrome_logs)}
        """.stripIndent()
    )
}

def gmsTomteHatLogo(monochrome_logs=true) {
    Map colors = logColours(monochrome_logs)
    String.format(
        """\n
        ${dashedLine(monochrome_logs)}
        ${colors.red}                                                                     @++@${colors.reset}
        ${colors.blue}   _____ __  __  _____            _                  _              ${colors.red}@    @${colors.reset}
        ${colors.blue}  / ____|  \\/  |/ ____|          | |                | |            ${colors.red}@      @@${colors.reset}
        ${colors.blue} | |  __| \\  / | (___    ______  | |_ ___  _ __ ___ | |_ ___      ${colors.red}@       @${colors.white}(+)${colors.reset}
        ${colors.blue} | | |_ | |\\/| |\\___ \\  |______| | __/ _ \\| '_ ` _ \\| __/ _ \\    ${colors.red}@        @${colors.reset}
        ${colors.blue} | |__| | |  | |____) |          | || (_) | | | | | | ||  __/    ${colors.red}@         @${colors.reset}
        ${colors.blue}  \\_____|_|  |_|_____/            \\__\\___/|_| |_| |_|\\__\\___|    ${colors.red}@         @${colors.reset}
        ${colors.white}                                                                 @===${colors.yellow}(@)${colors.white}===@${colors.reset}
        ${colors.purple}  ${workflow.manifest.name} ${getWorkflowVersion()}${colors.reset}
        ${dashedLine(monochrome_logs)}
        """.stripIndent()
    )
}

def gmsTomteHatLogo2(monochrome_logs=true) {
    Map colors = logColours(monochrome_logs)
    String.format(
        """\n
        ${dashedLine(monochrome_logs)}                 ${colors.red}_${colors.reset}
        ${colors.red}                                                                     /  \\${colors.reset}
        ${colors.blue}   _____ __  __  _____            _                  _              ${colors.red}/     \\${colors.reset}
        ${colors.blue}  / ____|  \\/  |/ ____|          | |                | |            ${colors.red}/      |\\${colors.reset}
        ${colors.blue} | |  __| \\  / | (___    ______  | |_ ___  _ __ ___ | |_ ___      ${colors.red}/       |${colors.white}(+)${colors.reset}
        ${colors.blue} | | |_ | |\\/| |\\___ \\  |______| | __/ _ \\| '_ ` _ \\| __/ _ \\    ${colors.red}|        \\${colors.reset}
        ${colors.blue} | |__| | |  | |____) |          | || (_) | | | | | | ||  __/    ${colors.red}|         |${colors.reset}
        ${colors.blue}  \\_____|_|  |_|_____/            \\__\\___/|_| |_| |_|\\__\\___|    ${colors.red}|_________|${colors.reset}
        ${colors.white}                                                                 {_________}${colors.reset}
        ${colors.purple}  ${workflow.manifest.name} ${getWorkflowVersion()}${colors.reset}
        ${dashedLine(monochrome_logs)}
        """.stripIndent()
    )
}
