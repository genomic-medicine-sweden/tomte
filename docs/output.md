# genomic-medicine-sweden/tomte: Output

## Introduction

This document describes the output produced by the pipeline. Most of the plots are taken from the MultiQC report, which summarises results at the end of the pipeline.

The directories listed below will be created in the results directory after the pipeline has finished. All paths are relative to the top-level results directory.

## Pipeline overview

The pipeline is built using [Nextflow](https://www.nextflow.io/) and processes data using the following steps:

- [`Trimming`](#trimming)
  - [`FASTP`](#fastp) trims reads
- [`Transcript quantification`](#transcript-quantification)
  - [`Salmon`](#salmon) quantifies transcripts
- [`Allignment`](#allignment)
  - [`STAR`](#star) aligns reads to the genome
- [`Tracks`](#tracks)
  - [`Tracks`](#tracks-1) outputs tracks
- [`Transcript analysis`](#transcript-analysis)
  - [`DROP`](#drop) aberrant expression and aberrant splicing discovery
  - [`StringTie`](#stringtie) guided transcript assembly
  - [`GffCompare`](#gffcompare) annnotation of guided transcript assembly
- [`Variant Calling`](#variant-calling)
  - [`BCFtools Mpileups`](#mpileups) single nucleotide variation calling
  - [`GATK best practices SNV Calling`](#gatk-best-practices-snv-calling)
- [`Allele specific variant Calling`](#allele-specific-variant-calling)
  - [`ASEReadCounter`](#asereadcounter) allele Specific Read Counter
  - [`BootstrapAnn`](#bootstrapann) assesses allelic imbalance
- [`Variant annotation`](#variant-annotation)
  - [`VEP`](#vep) annotation
- [`Pipeline information and QCs`](#pipeline-information-and-qcs)
  - [Pipeline information](#pipeline-information) - report metrics generated during the workflow execution
  - [`Picard CollectRnaSeqMetrics`](#picard-collectrnaseqmetrics) alignment QC
  - [`Picard CollectInsertSizeMetrics`](#picard-collectinsertsizemetrics) insert size
  - [`MultiQC`](#multiqc) presents QCs

### Trimming

#### FASTP

[FASTP](https://github.com/OpenGene/fastp) is a fastq preprocessing tool that gives general quality metrics about your sequenced reads and trims adapters from them. For further reading and documentation see the [FASTP documentation](https://github.com/OpenGene/fastp).

<details markdown="1">
<summary>Output files</summary>

- `trimming/`
  - `*.fastp.html`: a report consisting on a standalone HTML file that can be viewed in your web browser.
  - `*.fastp.log`: run log.
  - `*.fastp.json`: a report containing the same information as the html as a json file.
  - `*.fastp.fastq.gz`: gzip compressed trimmed reads.

</details>

### Transcript quantification

#### Salmon

[`Salmon`](https://salmon.readthedocs.io/en/latest/) quantifies reads.

<details markdown="1">
<summary>Output files</summary>

- `alignment/sample`
  - `quant.sf`: quantification file.
  - `quant.genes.sf`: quantification file per gene.
  - `logs/salmon_quant.log`: log file.
  - `cmd_info.json`: main command line parameters with which Salmon was run.

</details>

### Alignment

#### STAR

[`STAR`](https://github.com/alexdobin/STAR) aligns reads to the genome reference. For further reading and documentation see the [STAR manual](https://physiology.med.cornell.edu/faculty/skrabanek/lab/angsd/lecture_notes/STARmanual.pdf).

<details markdown="1">
<summary>Output files</summary>

- `alignment/`
  - `*.SJ.out.tab`: the high confidence collapssed junctions.
  - `*.ReadsPerGene.out.tab`: read count per gene.
  - `*.Log.progress.out`: run progress statistics report updated every minute.
  - `*.Log.out`: log file containing run details.
  - `*.Log.final.out`: a summary of the mapping statistics. It is calculated indivisually per read and then averaged.
  - `*.Aligned.out.bam`: Aligned reads.

</details>

### Tracks

#### Tracks

Outputs both junction tracks and bigwig files. For wigToBigWig [`UCSC wigToBigWig`](https://genome.ucsc.edu/goldenPath/help/bigWig.html) is used.

<details markdown="1">
<summary>Output files</summary>

- `ucsc/`
  - `*.bw`: track in bigwig format.
  - `*_junction.bed`: junction bed.
  - `*_bed.gz`: bed file with sample data.
  - `*_bed.gz.tbi`: index for bed file with sample data.

</details>

### Transcript analysis

#### DROP

[`DROP`](https://github.com/gagneurlab/drop/) is a pipleine that detects aberrant expression, aberrant spliceing, and monoallelic expression. For the time being, aberrant expression and aberrant splicing modules are run. Afterwards another script is run to filter results.

<details markdown="1">
<summary>Output files</summary>

- `analyse_transcripts/drop`
  - `OUTRIDER_provided_samples_top_hits.tsv`: provides at least the top 20 most significant events reported by OUTRIDER in each sample.
  - `OUTRIDER_provided_samples_top_hits_filtered.tsv`: filters OUTRIDER_provided_samples_top_hits according to genes provided by gene_panel_clinical_filter.
  - `FRASER_provided_samples_top_hits.tsv`: provides the aberrant spliced events reported by FRASER.
  - `FRASER_provided_samples_top_hits_filtered.tsv`: filters FRASER_provided_samples_top_hits according to genes provided by gene_panel_clinical_filter.

</details>

#### StringTie

[`StringTie`](https://ccb.jhu.edu/software/stringtie/) will perform guided transcript assembly.

<details markdown="1">
<summary>Output files</summary>

- `analyse_transcripts`
  - `*.coverage.gtf`: coverage on the sample.
  - `*.gene.abundance.txt`: gene abundance on the sample.
  - `*.transcripts.gtf`: transcripts assembled on the sample

</details>

#### GffCompare

[`GffCompare`](https://github.com/gpertea/gffcompare) annotates stringtie results with the reference transcripts, marking each assembled transcript as either normal or aberrant.

<details markdown="1">
<summary>Output files</summary>

- `analyse_transcripts`
  - `*.stats`: data summary and accuracy estimation.
  - `*.annotated.gtf`: annotated gtf file.
  - `*.tracking`: transcripts assembled on the sample
  - `*.transcripts.gtf.refmap`: list for each reference transcript what query transcript partially or fully matches it.
  - `*.transcripts.gtf.tmap`: list the most similar reference transcript to each query transcript.

</details>

### Variant Calling

#### Mpileups

[`BCFtools Mpileups`](https://samtools.github.io/bcftools/bcftools.html#mpileup) SNV calling. Default SNV caller.

<details markdown="1">
<summary>Output files</summary>

- `call variants`
  - `*.vcf.gz`: file in vcf format containing variants found in the patient.
  - `*.vcf.gz.tbi`: index for .vcf.gz file.
  - `*.bcftools_stats.txt`: stats on non-reference allele frequency, depth distribution, stats by quality and per-sample counts, singleton stats, etc.

</details>

#### GATK best practices SNV Calling

[`GATK best practices SNV Calling`](https://gatk.broadinstitute.org/hc/en-us/articles/360035531192-RNAseq-short-variant-discovery-SNPs-Indels-) SNV calling will only be activated by setting parameter variant_caller
to "gatk". Involves several steps: [`SplitN Cigar Reads`](https://gatk.broadinstitute.org/hc/en-us/articles/360036858811-SplitNCigarReads), [`Haplotype Caller`](https://gatk.broadinstitute.org/hc/en-us/articles/360037225632-HaplotypeCaller), [`Variant Filtration`](https://gatk.broadinstitute.org/hc/en-us/articles/360037434691-VariantFiltration) and [`BCFtools stats`](https://samtools.github.io/bcftools/bcftools.html#stats).

<details markdown="1">
<summary>Output files</summary>

- `call variants`
  - `*.vcf.gz`: file in vcf format containing variants found in the patient.
  - `*.vcf.gz.tbi`: index for .vcf.gz file.
  - `*.bcftools_stats.txt`: stats on non-reference allele frequency, depth distribution, stats by quality and per-sample counts, singleton stats, etc.

</details>

### Allele specific variant calling

#### ASEReadCounter

[`ASEReadCounter`](https://gatk.broadinstitute.org/hc/en-us/articles/360037428291-ASEReadCounter) allele Specific Read Counter.

#### BootstrapAnn

[`BootstrapAnn`](https://github.com/J35P312/BootstrapAnn#bootstrapann) detects expression imbalance between alleles.

<details markdown="1">
<summary>Output files</summary>

- `bootstrapann`
  - `*ase.vcf`: annotated vcf where allelic imbalance is marked

</details>

### Variant annotation

#### VEP

[`VEP`](https://github.com/Ensembl/ensembl-vep) annotates vcfs.

<details markdown="1">
<summary>Output files</summary>

- `annotate_vep`
  - `*ase_vep.vcf.gz`: annotated vcf
  - `*ase_vep.vcf.gz.tbi`: index for annotated vcf

</details>

### Pipeline information and QCs

#### Pipeline information

[Nextflow](https://www.nextflow.io/docs/latest/tracing.html) provides excellent functionality for generating various reports relevant to the running and execution of the pipeline. This will allow you to troubleshoot errors with the running of the pipeline, and also provide you with other information such as launch commands, run times and resource usage.

<details markdown="1">
<summary>Output files</summary>

- `pipeline_info/`
  - Reports generated by Nextflow: `execution_report.html`, `execution_timeline.html`, `execution_trace.txt` and `pipeline_dag.dot`/`pipeline_dag.svg`.
  - Reports generated by the pipeline: `pipeline_report.html`, `pipeline_report.txt` and `software_versions.yml`. The `pipeline_report*` files will only be present if the `--email` / `--email_on_fail` parameter's are used when running the pipeline.
  - Reformatted samplesheet files used as input to the pipeline: `samplesheet.valid.csv`.
  - Parameters used by the pipeline run: `params.json`.

</details>

#### Picard CollectRnaSeqMetrics

[`Picard CollectRnaSeqMetrics`](https://broadinstitute.github.io/picard/) alignment QC

<details markdown="1">
<summary>Output files</summary>

- `bam_qc/`
  - `*rna_metrics`: metrics describing the distribution of the bases within the transcripts.
  - `*_insert_size.txt`: metrics describing the insert size.

#### Picard CollectInsertSizeMetrics

[`Picard CollectInsertSizeMetrics`](https://broadinstitute.github.io/picard/)

<details markdown="1">
<summary>Output files</summary>

- `bam_qc/`
  - `*rna_metrics`: metrics describing the distribution of the bases within the transcripts.
  - `sample.txt`: metrics describing the insert size.

</details>

#### MultiQC

[MultiQC](http://multiqc.info) is a visualization tool that generates a single HTML report summarising all samples in your project. Most of the pipeline QC results are visualised in the report and further statistics are available in the report data directory.

Results generated by MultiQC collate pipeline QC from supported tools e.g. FastQC. The pipeline has special steps which also allow the software versions to be reported in the MultiQC output for future traceability. For more information about how to use MultiQC reports, see <http://multiqc.info>.

<details markdown="1">
<summary>Output files</summary>

- `multiqc/`
  - `multiqc_report.html`: a standalone HTML file that can be viewed in your web browser.
  - `multiqc_data/`: directory containing parsed statistics from the different tools used in the pipeline.
  - `multiqc_plots/`: directory containing static images from the report in various formats.

</details>
