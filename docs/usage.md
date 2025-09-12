# genomic-medicine-sweden/tomte: Usage

> _Documentation of pipeline parameters is generated automatically from the pipeline schema and can no longer be found in markdown files._

Table of contents:

- [genomic-medicine-sweden/tomte: Usage](#genomic-medicine-swedentomte-usage)
  - [Introduction](#introduction)
  - [Prerequisites](#prerequisites)
  - [Run genomic-medicine-sweden/tomte with test data](#run-genomic-medicine-swedentomte-with-test-data)
    - [Updating the pipeline](#updating-the-pipeline)
  - [Run genomic-medicine-sweden/tomte with your data](#run-genomic-medicine-swedentomte-with-your-data)
    - [Samplesheet](#samplesheet)
    - [Reference files and parameters](#reference-files-and-parameters)
      - [Alignment and pseudo quantification](#1-alignment)
      - [Junction track and bigwig generation](#2-junction-track-and-bigwig)
      - [Region subsampling](#3-subsample-region)
      - [Variant calling](#4-variant-calling---snv)
      - [SNV annotation](#5-snv-annotation-ensembl-vep)
      - [Stringtie & gffcompare](#6-stringtie-and-gffcompare)
      - [DROP](#7-drop)
        - [Preparing DROP input](#preparing-input-for-drop)
  - [Run the pipeline](#run-the-pipeline)
    - [Direct input in CLI](#direct-input-in-cli)
    - [Import from a config file (recommended)](#import-from-a-config-file-recommended)
- [Best practices](#best-practices)
- [Core Nextflow arguments](#core-nextflow-arguments)
  - [`-profile`](#-profile)
  - [`-resume`](#-resume)
  - [`-c`](#-c)
- [Custom configuration](#custom-configuration)
  - [Changing resources](#changing-resources)
  - [Custom Containers](#custom-containers)
  - [Custom Tool Arguments](#custom-tool-arguments)
    - [nf-core/configs](#nf-coreconfigs)
  - [Azure Resource Requests](#azure-resource-requests)
  - [Running in the background](#running-in-the-background)
  - [Nextflow memory requirements](#nextflow-memory-requirements)
  - [Running the pipeline without Internet access](#running-the-pipeline-without-internet-access)

## Introduction

**tomte** is a bioinformatics best-practice analysis pipeline for analysing RNAseq data from patients with rare diseases.

The pipeline is built using [Nextflow](https://www.nextflow.io), a workflow tool to run tasks across multiple compute infrastructures in a very portable manner. It uses Docker/Singularity containers making installation trivial and results highly reproducible. The [Nextflow DSL2](https://www.nextflow.io/docs/latest/dsl2.html) implementation of this pipeline uses one container per process which makes it much easier to maintain and update software dependencies. Where possible, these processes have been submitted to and installed from [nf-core/modules](https://github.com/nf-core/modules) in order to make them available to all nf-core pipelines, and to everyone within the Nextflow community!

## Prerequisites

1. Install Nextflow (>=22.10.1) using the instructions [here.](https://nextflow.io/docs/latest/getstarted.html#installation)
2. Install one of the following technologies for full pipeline reproducibility: Docker, Singularity, Podman, Shifter or Charliecloud.
   > Almost all nf-core pipelines give you the option to use conda as well. However, some tools used in the tomte pipeline do not have a conda package so we do not support conda at the moment.

### Updating the pipeline

When you run the above command, Nextflow automatically pulls the pipeline code from GitHub and stores it as a cached version. When running the pipeline after this, it will always use the cached version if available - even if the pipeline has been updated since. To make sure that you're running the latest version of the pipeline, make sure that you regularly update the cached version of the pipeline:

```bash
nextflow pull genomic-medicine-sweden/tomte
```

## Run genomic-medicine-sweden/tomte with test data

Before running the pipeline with your data, we recommend running it with the test dataset available in the test_data folder provided with the pipeline and [here](https://github.com/nf-core/test-datasets/tree/raredisease). You do not need to download any of the data as part of it came directly with the pipeline and the other part will be fetched automatically for you when you use the test profile.

Run the following command, where YOURPROFILE is the package manager you installed on your machine. For example, `-profile test,docker` or `-profile test,singularity`:

```
nextflow run genomic-medicine-sweden/tomte \
    -revision dev -profile test,<YOURPROFILE> \
    --outdir <OUTDIR>
```

> Check [nf-core/configs](https://github.com/nf-core/configs/tree/master/conf) to see if a custom config file to run nf-core pipelines already exists for your institute. If so, you can simply use `-profile test,<institute>` in your command. This enables the appropriate package manager and sets the appropriate execution settings for your machine.
> NB: The order of profiles is important! They are loaded in sequence, so later profiles can overwrite earlier profiles.

Running the command creates the following files in your working directory:

```
work                # Directory containing the Nextflow working files
<OUTDIR>            # Finished results in specified location (defined with --outdir)
.nextflow_log       # Log file from Nextflow
# Other Nextflow hidden files, like history of pipeline logs.
```

Test profile runs the pipeline with a case containing three samples, but if you would like to test the pipeline with one sample, use `-profile test_one_sample,<YOURPROFILE>`.

> Note that the default cpu and memory configurations used in tomte are written keeping the test profile (&dataset, which is tiny) in mind. You should override these values in configs to get it to work on larger datasets. Check the section `custom-configuration` below to know more about how to configure resources for your platform.

## Run genomic-medicine-sweden/tomte with your data

Running the pipeline involves three steps:

1. Prepare a samplesheet
2. Gather all required references
3. Supply samplesheet and references, and run the command

#### Samplesheet

A samplesheet is used to pass the information about the sample(s), such as the path to the FASTQ/BAM/CRAM files and other meta data (sex, phenotype, etc.,) to the pipeline in csv format.

genomic-medicine-sweden/tomte will requires the information given bellow.

| Fields         | Description                                                                                                                                                                            | Mandatory?                         |
| -------------- | -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- | ---------------------------------- |
| `case`         | Case ID, for the analysis used when generating a family VCF.                                                                                                                           | Mandatory                          |
| `sample`       | Custom sample name. This entry will be identical for multiple sequencing libraries/runs from the same sample. Spaces in sample names are automatically converted to underscores (`_`). | Mandatory                          |
| `strandedness` | Sample strandness                                                                                                                                                                      | Mandatory                          |
| `fastq_1`      | Absolute path to FASTQ file for Illumina short reads 1. File has to be gzipped and have the extension ".fastq.gz" or ".fq.gz".                                                         | Provide either fastq_1 or bam/cram |
| `fastq_2`      | Absolute path to FASTQ file for Illumina short reads 2. File has to be gzipped and have the extension ".fastq.gz" or ".fq.gz".                                                         | Only if paired fastqs are provided |
| `bam_cram`     | Full path to BAM/CRAM file.                                                                                                                                                            | Provide either fastq_1 or bam/cram |
| `bai_crai`     | Full path to BAM/CRAM index file.                                                                                                                                                      | Only when bam/cram is provided     |
| `Paternal`     | Father's custom sample name. If there is no paternal sample available, it can be left empty.                                                                                           | Optional                           |
| `Maternal`     | Mother's custom sample name. If there is no maternal sample available, it can be left empty.                                                                                           | Optional                           |
| `Sex`          | Sample sex. The valid input is M or 1 for male; F or 2 for female; NA, 0, or other if unknown                                                                                          | Optional                           |
| `dna_vcf`      | Full path to DNA vcf file to run DROP's Mono Allelic expression (MAE) module file.                                                                                                     | Only if you want to run MAE        |
| `dna_vcf_tbi`  | Full path to DNA vcf file's index file.                                                                                                                                                | Only if you want to run MAE        |

It is also possible to include multiple runs of the same sample in a samplesheet. For example, when you have re-sequenced the same sample more than once to increase sequencing depth. In that case, the `sample` identifiers in the samplesheet have to be the same. The pipeline will align the raw read/read-pairs independently before merging the alignments belonging to the same sample. Below is an example for a trio with the proband sequenced across two lanes:

| case  | sample       | strandedness | fastq_1                          | fastq_2                          | Paternal     | Maternal     | Sex | dna_vcf             | dna_vcf_tbi             |
| ----- | ------------ | ------------ | -------------------------------- | -------------------------------- | ------------ | ------------ | --- | ------------------- | ----------------------- |
| fam_1 | CONTROL_REP1 | reverse      | AEG588A1_S1_L002_R1_001.fastq.gz | AEG588A1_S1_L002_R2_001.fastq.gz |              |              | M   | AEG588A1_DNA.vcf.gz | AEG588A1_DNA.vcf.gz.tbi |
| fam_1 | CONTROL_REP2 | reverse      | AEG588A2_S1_L003_R1_001.fastq.gz | AEG588A2_S1_L003_R2_001.fastq.gz |              |              | F   |                     |                         |
| fam_1 | PATIENT_1    | reverse      | AEG588A3_S1_L001_R1_001.fastq.gz | AEG588A3_S1_L001_R2_001.fastq.gz | CONTROL_REP1 | CONTROL_REP2 | M   |                     |                         |
| fam_1 | PATIENT_1    | reverse      | AEG588A3_S1_L002_R1_001.fastq.gz | AEG588A3_S1_L002_R2_001.fastq.gz | CONTROL_REP1 | CONTROL_REP2 | M   |                     |                         |

Here is an example of a samplesheet where BAM files are provided:

| case  | sample       | strandedness | fastq_1 | fastq_2 | bam_cram     | bai_crai         |
| ----- | ------------ | ------------ | ------- | ------- | ------------ | ---------------- |
| fam_1 | CONTROL_REP1 | reverse      |         |         | AEG588A1.bam | AEG588A1.bam.bai |
| fam_1 | CONTROL_REP2 | reverse      |         |         | AEG588A2.bam | AEG588A2.bam.bai |

Here is an example of a samplesheet where CRAM files are provided:

| case  | sample       | strandedness | fastq_1 | fastq_2 | bam_cram      | bai_crai           |
| ----- | ------------ | ------------ | ------- | ------- | ------------- | ------------------ |
| fam_1 | CONTROL_REP1 | reverse      |         |         | AEG588A1.cram | AEG588A1.cram.crai |
| fam_1 | CONTROL_REP2 | reverse      |         |         | AEG588A2.cram | AEG588A2.cram.crai |

If you would like to see more examples of what a typical samplesheet looks like for a trio, follow this link, [sample_sheet](https://github.com/genomic-medicine-sweden/tomte/blob/master/test_data/samplesheet_chr21.csv)

#### Reference files and parameters

In genomic-medicine-sweden/tomte, references can be supplied using parameters. We have also introduced the possiblility of using the `--igenomes_base` parameter to point to a path where genome specific reference files are placed (fasta, fai, gtf, star_index, salmon_index, subsample_bed). To make sure that the names of the reference files match those in your directory, check [igenomes.config](https://github.com/genomic-medicine-sweden/tomte/blob/master/conf/igenomes.config).

If no references are provided by the user the pipeline will automatically download a fasta and a gtf file. The user can select the desired genome and gencode version using `--genome` and `--genome_annotation_version`. If the user also wants to download vep cache and vep plugins references they will have to set `--skip_download_vep false`. The user will have to provide a comma separated file containing the plugins they want to download `--vep_refs_download`, this file should NOT contain the path to gnomad database. If the user also wants to download the gnomad database they will have to set `--skip_download_gnomad false`, bare in mind that about ~900GB of data will be downloaded, so storage space and time are needed. The data will then be processed and its size significantly reduced to under 40GB.

Note that the pipeline is modular in architecture. It offers you the flexibility to choose between different tools. For example, you can call SNVs either with BCFtools or with GATK. You also have the option to turn off sections of the pipeline if you do not want to run them. For example, drop aberrant expression module can be turned off by setting `--skip_drop_ae true`. This flexibility means that in any given analysis run, a combination of tools included in the pipeline will not be executed. So the pipeline is written in a way that can account for these differences while working with reference parameters. If a tool is not going to be executed during the course of a run, parameters used only by that tool need not be provided. For example, if you are not running DROP aberrant splicing, you do not need to provide `--reference_drop_splice_folder`.

genomic-medicine-sweden/tomte consists of several tools used for various purposes. For convenience, we have grouped those tools under the following categories:

1. Alignment and pseudo quantification (STAR & Salmon)
2. Junction track and bigwig
3. Subsample region (Samtools)
4. Variant calling - SNV (BCFTools or GATK's HaplotypeCaller)
5. SNV annotation (ensembl VEP)
6. Stringtie & gffcompare
7. DROP

> We have only listed the groups that require at least one input from the user. For example, the pipeline also runs WigToBigWig, but it does not require any input other than the bam files passed by the pipeline. Hence, it is not mentioned in the list above. To know more about the tools used in the pipeline check the [README](../README.md).

The mandatory and optional parameters for each category are tabulated below.

> Alignment, QC stats, repeat expansions, SNV variant calling and ensembl VEP are run by default. Hence, the mandatory parameters used by those features will always have to be provided to the pipeline.

##### 1. Alignment

| Mandatory | Optional                       |
| --------- | ------------------------------ |
|           | fasta<sup>1</sup>              |
|           | gtf<sup>1</sup>                |
|           | fasta_fai<sup>2</sup>          |
|           | sequence_dict<sup>2</sup>      |
|           | salmon_index<sup>2</sup>       |
|           | star_index<sup>2</sup>         |
|           | transcript_fasta<sup>2</sup>   |
|           | genome<sup>3</sup>             |
|           | platform<sup>4</sup>           |
|           | min_trimmed_length<sup>5</sup> |
|           | star_two_pass_mode<sup>6</sup> |

<sup>1</sup> If the parameter is not provided by the user, it will be downloaded.<br />
<sup>2</sup> If the parameter is not provided by the user, it will be generated from the fasta and gtf files.<br />
<sup>3</sup> If it is not provided by the user, the default value is GRCh38.<br />
<sup>4</sup> If it is not provided by the user, the default value is illumina.<br />
<sup>5</sup> If it is not provided by the user, the default value is 40.<br />
<sup>6</sup> If it is not provided by the user, the default value is Basic.

##### 2. Junction track and bigwig

| Mandatory | Optional                       |
| --------- | ------------------------------ |
|           | skip_build_tracks <sup>1</sup> |

<sup>1</sup> If it is not provided by the user, the default value is false

##### 3. Subsample region

| Mandatory     | Optional                          |
| ------------- | --------------------------------- |
| subsample_bed | skip_subsample_region<sup>1</sup> |
|               | seed_frac<sup>2</sup>             |

<sup>1</sup> If it is not provided by the user, the default value is false
<sup>2</sup> If it is not provided by the user, the default value is 0.001

##### 4. Variant calling - SNV

| Mandatory | Optional                         |
| --------- | -------------------------------- |
|           | variant_caller<sup>1</sup>       |
|           | bcftools_caller_mode<sup>2</sup> |
|           | skip_variant_calling<sup>3</sup> |

<sup>1</sup> If it is not provided by the user, the default value is bcftools<br />
<sup>2</sup> If it is not provided by the user, the default value is multiallelic<br />
<sup>3</sup> If it is not provided by the user, the default value is false

#### 5. SNV annotation (ensembl VEP)

| Mandatory | Optional                         |
| --------- | -------------------------------- |
|           | skip_vep<sup>1</sup>             |
|           | vep_plugin_files<sup>2</sup>     |
|           | vep_cache<sup>2</sup>            |
|           | vep_cache_version<sup>3</sup>    |
|           | skip_download_vep<sup>4</sup>    |
|           | skip_download_gnomad<sup>4</sup> |
|           | vep_refs_download                |
|           | gene_panel_clinical_filter       |

<sup>1</sup> If it is not provided by the user, the default value is false<br />
<sup>2</sup> VEP cache and plugins can be automatically downloaded by the pipeline by setting `--skip_download_vep false`, `--skip_download_gnomad false` and providing a lcsv with a list of files to download `--vep_refs_download` as done [here](https://github.com/genomic-medicine-sweden/tomte/blob/dev/test_data/vep_to_download.csv). VEP caches can also be downloaded [here](https://www.ensembl.org/info/docs/tools/vep/script/vep_cache.html#cache). VEP plugins may also be installed in the cache directory, and the plugin pLI is mandatory to install. To supply files required by VEP plugins, use `vep_plugin_files` parameter. See example cache [here](https://raw.githubusercontent.com/nf-core/test-datasets/raredisease/reference/vep_cache_and_plugins.tar.gz).<br />
<sup>3</sup> If it is not provided by the user, the default value is 112, supported values are 107, 110, and 112
<sup>4</sup> If it is not provided by the user, the default value true

#### 6. Stringtie and gffcompare

| Mandatory | Optional                   |
| --------- | -------------------------- |
| fasta     | skip_stringtie<sup>1</sup> |
| gtf       |                            |

<sup>1</sup> If it is not provided by the user, the default value is false

#### 7. DROP

DROP - aberrant expression

| Mandatory                             | Optional                            |
| ------------------------------------- | ----------------------------------- |
| reference_drop_annot_file<sup>1</sup> | skip_drop_ae<sup>2</sup>            |
| reference_drop_count_file             | drop_group_samples_ae<sup>3</sup>   |
| fasta                                 | drop_padjcutoff_ae<sup>4</sup>      |
| gtf                                   | drop_zscorecutoff<sup>5</sup>       |
|                                       | gene_panel_clinical_filter          |
|                                       | skip_downsample<sup>6</sup>         |
|                                       | num_reads<sup>7</sup>               |
|                                       | genome<sup>8</sup>                  |
|                                       | skip_export_counts_drop<sup>9</sup> |

<sup>1</sup> To get more information on how to format it, see below<br />
<sup>2</sup> If it is not provided by the user, the default value is false<br />
<sup>3</sup> If it is not provided by the user, the default value is outrider<br />
<sup>4</sup> If it is not provided by the user, the default value is 0.05<br />
<sup>5</sup> If it is not provided by the user, the default value is 0<br />
<sup>6</sup> If it is not provided by the user, the default value is false<br />
<sup>7</sup> If it is not provided by the user, the default value is 120000000<br />
<sup>8</sup> If it is not provided by the user, the default value is GRCh38
<sup>9</sup> If it is not provided by the user, the default value is true<br />

DROP - aberrant splicing

| Mandatory                             | Optional                            |
| ------------------------------------- | ----------------------------------- |
| reference_drop_annot_file<sup>1</sup> | skip_drop_as<sup>2</sup>            |
| reference_drop_splice_folder          | drop_group_samples_as<sup>3</sup>   |
|                                       | drop_padjcutoff_as<sup>4</sup>      |
|                                       | gene_panel_clinical_filter          |
|                                       | skip_downsample<sup>5</sup>         |
|                                       | num_reads<sup>6</sup>               |
|                                       | genome<sup>7</sup>                  |
|                                       | skip_export_counts_drop<sup>8</sup> |

<sup>1</sup> To get more information on how to format it, see below<br />
<sup>2</sup> If it is not provided by the user, the default value is false<br />
<sup>3</sup> If it is not provided by the user, the default value is fraser<br />
<sup>4</sup> If it is not provided by the user, the default value is 0.1<br />
<sup>5</sup> If it is not provided by the user, the default value is false<br />
<sup>6</sup> If it is not provided by the user, the default value is 120000000<br />
<sup>7</sup> If it is not provided by the user, the default value is GRCh38
<sup>8</sup> If it is not provided by the user, the default value is true<br />

DROP - monoallelic expression

| Mandatory                                     | Optional                            |
| --------------------------------------------- | ----------------------------------- |
| reference_drop_annot_file<sup>1</sup>         | drop_mae_high_q_vcf<sup>2</sup>     |
| variant calling from WGS (vcf/vcf.gz)         | drop_mae_high_q_vcf_tbi<sup>2</sup> |
| variant calling from WGS (vcf.tbi/vcf.gz.tbi) | gene_panel_clinical_filter          |
|                                               | genome<sup>3</sup>                  |

<sup>1</sup> To get more information on how to format it, see below<br />
<sup>2</sup> If it is not provided by the user, the user can chose to download it by ` --skip_download_drop_mae_high_q_vcf false`<br />
<sup>7</sup> If it is not provided by the user, the default value is GRCh38

##### Preparing input for DROP

If you want to run [DROP](https://github.com/gagneurlab/drop) aberrant expression or aberrant splicing you have to provide reference counts, splice counts, and a sample sheet. The sample sheet should contain the columns as those in the [test sample annotation](../test_data/drop_data/sampleAnnotation.tsv), you can also add an optional sex column. You do not need to include the samples you are running through the pipeline in the sample sheet.

IIf you want to run the DROP Mono Allelic Expression (MAE) module, you do NOT need to provide reference counts, splice counts, or a traditional sample sheet. Instead, the sample sheet should include the VCF or VCF.GZ file obtained from the DNA sample corresponding to your RNA samples. Additionally, you can provide a file specifying regions used to confirm that both DNA and RNA are from the same individual. If this file is not provided, it will be downloaded automatically if `--skip_download_drop_mae_high_q_vcf false` is specified

###### Preparing your DROP control database

You have several options on how to create such a database. You can either build it or download it from one of the [available databases](https://github.com/gagneurlab/drop#datasets).

To build your own database you will need at least 50 samples for aberrant expression, if you only run aberrant splicing 30 samples will suffice but DROP authors recommend to have at least around 100 for both modules. You can use Tomte to build your own database, to do so we recommend to run with the following parameters:

- `--skip_export_counts_drop false` this switch will ensure that a folder in references called export_counts is created
- `--skip_drop_as false` if you want to get a database for aberrant splicing
- `--skip_drop_ae false` if you want to get a database for aberrant expression
- `--skip_subsample_region false` if you have sequenced any material with overrepresented regions (such as hemoglobin in whole blood) we recommend to remove it by setting this parameter to false and providing a bed with the overrepresented region with `--subsample_bed`
- `--skip_downsample false` if you have very deeply sequenced samples, we recommend to downsample, the default is 60M read pairs (or 120M reads).
- `--skip_build_tracks true`, `--skip_stringtie true`, `--skip_variant_calling`, `--skip_vep true` as most users will be interested in getting the database rather than other downstream results, which will require a considerable amount of resources given the number of samples run.

Running DROP with many samples requires a lot of time and a lot of memory, that is why we recommend to subsample overrepresented regions and downsample if you have deeply sequenced samples. If your run fails, we recommend increasing memory without increasing the number of cores. If it fails and has been running for a while, try to relaunch it from the work directory where DROP was run so that DROP continues from the point where it failed (if you restart the pipeline with `-resume` it will begin from the start).

To restart DROP, start by finding the work directory where it was run. You can do so by opening the execution trace file in the pipeline_info folder and looking at the hash of the processes with name `TOMTE:ANALYSE_TRANSCRIPTS:DROP_CONFIG_RUN_AE` and `TOMTE:ANALYSE_TRANSCRIPTS:DROP_CONFIG_RUN_AS`. The work directory used to run the pipeline followed by the hash should be enough information to find the folder where DROP was run. Tomte should have set everything up in that directory so go into it and restart the run by running from the container created by Tomte the script `.command.sh`. If you want to run it with slurm remember to add a header with number of cores, time...

If you want to add samples to an existing database, follow the same steps described above, making sure that you also provide the database you want to add samples to by using `--reference_drop_annot_file` and `--reference_drop_count_file` and/or `--reference_drop_splice_folder`. In this case scenerio, make sure that you have used the same references for the database as for the new set of samples.

If you prefer to run DROP locally outside from Tomte follow instructions given by the [authors of DROP](https://github.com/gagneurlab/drop)

###### Running MonoAllelic Expression module

To run DROP MAE, you must provide the following files for each sample, generated after variant calling on Whole Genome Sequencing (WGS) data:
`vcf/vcf.gz` and `vcf.tbi/vcf.gz.tbi`. Tomte will automatically attempt to run MAE for any sample in the samplesheet that includes these files.
Additionally, DROP MAE requires a VCF and index file containing positions where SNPs are usually called at high quality. This positions will be used to verify that the DNA VCF and RNA BAM file originate from the same individual. To obtain these files you can:

- Download them yourself from [TUM's publid repository](https://www.cmm.in.tum.de/public/paper/drop_analysis/resource/) and provide them via `--drop_mae_high_q_vcf` and `--drop_mae_high_q_vcf_tbi`
- Set `--skip_download_drop_mae_high_q_vcf false` and have the pipeline do it for you

## Running the pipeline

The typical command for running the pipeline is as follows:

```bash
nextflow run genomic-medicine-sweden/tomte --input ./samplesheet.csv --outdir ./results --genome GRCh37 -profile docker
```

This will launch the pipeline with the `docker` configuration profile. See below for more information about profiles.

Note that the pipeline will create the following files in your working directory:

```bash
work                # Directory containing the nextflow working files
<OUTDIR>            # Finished results in specified location (defined with --outdir)
.nextflow_log       # Log file from Nextflow
# Other nextflow hidden files, eg. history of pipeline runs and old logs.
```

If you wish to repeatedly use the same parameters for multiple runs, rather than specifying each flag in the command, you can specify these in a params file.

Pipeline settings can be provided in a `yaml` or `json` file via `-params-file <file>`.

> [!WARNING]
> Do not use `-c <file>` to specify parameters as this will result in errors. Custom config files specified with `-c` must only be used for [tuning process resource specifications](https://nf-co.re/docs/usage/configuration#tuning-workflow-resources), other infrastructural tweaks (such as output directories), or module arguments (args).

The above pipeline run specified with a params file in yaml format:

```bash
nextflow run genomic-medicine-sweden/tomte -profile docker -params-file params.yaml
```

with:

```yaml title="params.yaml"
input: './samplesheet.csv'
outdir: './results/'
genome: 'GRCh37'
<...>
```

You can also generate such `YAML`/`JSON` files via [nf-core/launch](https://nf-co.re/launch).

### Updating the pipeline

When you run the above command, Nextflow automatically pulls the pipeline code from GitHub and stores it as a cached version. When running the pipeline after this, it will always use the cached version if available - even if the pipeline has been updated since. To make sure that you're running the latest version of the pipeline, make sure that you regularly update the cached version of the pipeline:

```bash
nextflow pull genomic-medicine-sweden/tomte
```

### Reproducibility

It is a good idea to specify the pipeline version when running the pipeline on your data. This ensures that a specific version of the pipeline code and software are used when you run your pipeline. If you keep using the same tag, you'll be running the same version of the pipeline, even if there have been changes to the code since.

First, go to the [genomic-medicine-sweden/tomte releases page](https://github.com/genomic-medicine-sweden/tomte/releases) and find the latest pipeline version - numeric only (eg. `1.3.1`). Then specify this when running the pipeline with `-r` (one hyphen) - eg. `-r 1.3.1`. Of course, you can switch to another version by changing the number after the `-r` flag.

This version number will be logged in reports when you run the pipeline, so that you'll know what you used when you look back in the future. For example, at the bottom of the MultiQC reports.

To further assist in reproducibility, you can use share and reuse [parameter files](#running-the-pipeline) to repeat pipeline runs with the same settings without having to write out a command with every single parameter.

> [!TIP]
> If you wish to share such profile (such as upload as supplementary material for academic publications), make sure to NOT include cluster specific paths to files, nor institutional specific profiles.

## Core Nextflow arguments

> [!NOTE]
> These options are part of Nextflow and use a _single_ hyphen (pipeline parameters use a double-hyphen)

### `-profile`

Use this parameter to choose a configuration profile. Profiles can give configuration presets for different compute environments.

Several generic profiles are bundled with the pipeline which instruct the pipeline to use software packaged using different methods (Docker, Singularity, Podman, Shifter, Charliecloud, Apptainer, Conda) - see below.

> [!IMPORTANT]
> We highly recommend the use of Docker or Singularity containers for full pipeline reproducibility, however when this is not possible, Conda is also supported.

The pipeline also dynamically loads configurations from [https://github.com/nf-core/configs](https://github.com/nf-core/configs) when it runs, making multiple config profiles for various institutional clusters available at run time. For more information and to check if your system is supported, please see the [nf-core/configs documentation](https://github.com/nf-core/configs#documentation).

Note that multiple profiles can be loaded, for example: `-profile test,docker` - the order of arguments is important!
They are loaded in sequence, so later profiles can overwrite earlier profiles.

If `-profile` is not specified, the pipeline will run locally and expect all software to be installed and available on the `PATH`. This is _not_ recommended, since it can lead to different results on different machines dependent on the computer environment.

- `test`
  - A profile with a complete configuration for automated testing
  - Includes links to test data so needs no other parameters
- `docker`
  - A generic configuration profile to be used with [Docker](https://docker.com/)
- `singularity`
  - A generic configuration profile to be used with [Singularity](https://sylabs.io/docs/)
- `podman`
  - A generic configuration profile to be used with [Podman](https://podman.io/)
- `shifter`
  - A generic configuration profile to be used with [Shifter](https://nersc.gitlab.io/development/shifter/how-to-use/)
- `charliecloud`
  - A generic configuration profile to be used with [Charliecloud](https://hpc.github.io/charliecloud/)
- `apptainer`
  - A generic configuration profile to be used with [Apptainer](https://apptainer.org/)
- `wave`
  - A generic configuration profile to enable [Wave](https://seqera.io/wave/) containers. Use together with one of the above (requires Nextflow ` 24.03.0-edge` or later).
- `conda`
  - A generic configuration profile to be used with [Conda](https://conda.io/docs/). Please only use Conda as a last resort i.e. when it's not possible to run the pipeline with Docker, Singularity, Podman, Shifter, Charliecloud, or Apptainer.

### `-resume`

Specify this when restarting a pipeline. Nextflow will use cached results from any pipeline steps where the inputs are the same, continuing from where it got to previously. For input to be considered the same, not only the names must be identical but the files' contents as well. For more info about this parameter, see [this blog post](https://www.nextflow.io/blog/2019/demystifying-nextflow-resume.html).

You can also supply a run name to resume a specific run: `-resume [run-name]`. Use the `nextflow log` command to show previous run names.

### `-c`

Specify the path to a specific config file (this is a core Nextflow command). See the [nf-core website documentation](https://nf-co.re/usage/configuration) for more information.

## Custom configuration

### Resource requests

Whilst the default requirements set within the pipeline will hopefully work for most people and with most input data, you may find that you want to customise the compute resources that the pipeline requests. Each step in the pipeline has a default set of requirements for number of CPUs, memory and time. For most of the pipeline steps, if the job exits with any of the error codes specified [here](https://github.com/nf-core/rnaseq/blob/4c27ef5610c87db00c3c5a3eed10b1d161abf575/conf/base.config#L18) it will automatically be resubmitted with higher resources request (2 x original, then 3 x original). If it still fails after the third attempt then the pipeline execution is stopped.

To change the resource requests, please see the [max resources](https://nf-co.re/docs/usage/configuration#max-resources) and [tuning workflow resources](https://nf-co.re/docs/usage/configuration#tuning-workflow-resources) section of the nf-core website.

### Custom Containers

In some cases, you may wish to change the container or conda environment used by a pipeline steps for a particular tool. By default, nf-core pipelines use containers and software from the [biocontainers](https://biocontainers.pro/) or [bioconda](https://bioconda.github.io/) projects. However, in some cases the pipeline specified version maybe out of date.

To use a different container from the default container or conda environment specified in a pipeline, please see the [updating tool versions](https://nf-co.re/docs/usage/configuration#updating-tool-versions) section of the nf-core website.

### Custom Tool Arguments

A pipeline might not always support every possible argument or option of a particular tool used in pipeline. Fortunately, nf-core pipelines provide some freedom to users to insert additional parameters that the pipeline does not include by default.

To learn how to provide additional arguments to a particular tool of the pipeline, please see the [customising tool arguments](https://nf-co.re/docs/usage/configuration#customising-tool-arguments) section of the nf-core website.

### nf-core/configs

In most cases, you will only need to create a custom config as a one-off but if you and others within your organisation are likely to be running nf-core pipelines regularly and need to use the same settings regularly it may be a good idea to request that your custom config file is uploaded to the `nf-core/configs` git repository. Before you do this please can you test that the config file works with your pipeline of choice using the `-c` parameter. You can then create a pull request to the `nf-core/configs` repository with the addition of your config file, associated documentation file (see examples in [`nf-core/configs/docs`](https://github.com/nf-core/configs/tree/master/docs)), and amending [`nfcore_custom.config`](https://github.com/nf-core/configs/blob/master/nfcore_custom.config) to include your custom profile.

See the main [Nextflow documentation](https://www.nextflow.io/docs/latest/config.html) for more information about creating your own configuration files.

If you have any questions or issues please send us a message on [Slack](https://nf-co.re/join/slack) on the [`#configs` channel](https://nfcore.slack.com/channels/configs).

## Running in the background

Nextflow handles job submissions and supervises the running jobs. The Nextflow process must run until the pipeline is finished.

The Nextflow `-bg` flag launches Nextflow in the background, detached from your terminal so that the workflow does not stop if you log out of your session. The logs are saved to a file.

Alternatively, you can use `screen` / `tmux` or similar tool to create a detached session which you can log back into at a later time.
Some HPC setups also allow you to run nextflow within a cluster job submitted your job scheduler (from where it submits more jobs).

## Nextflow memory requirements

In some cases, the Nextflow Java virtual machines can start to request a large amount of memory.
We recommend adding the following line to your environment to limit this (typically in `~/.bashrc` or `~./bash_profile`):

```bash
NXF_OPTS='-Xms1g -Xmx4g'
```
