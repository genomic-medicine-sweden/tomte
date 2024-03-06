# genomic-medicine-sweden/tomte: Usage

## :warning: Please read this documentation on github website: [tomte usage](https://github.com/genomic-medicine-sweden/tomte)

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
      - [Region subsampling](#2-subsample-region)
      - [Variant calling](#3-variant-calling---snv)
      - [SNV annotation](#4-snv-annotation-ensembl-vep)
      - [DROP](#5-drop)
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

A samplesheet is used to pass the information about the sample(s), such as the path to the FASTQ files and other meta data (sex, phenotype, etc.,) to the pipeline in csv format.

genomic-medicine-sweden/tomte will requires the information given bellow.

| Fields         | Description                                                                                                                                                                            |
| -------------- | -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `case`         | Case ID, for the analysis used when generating a family VCF.                                                                                                                           |
| `sample`       | Custom sample name. This entry will be identical for multiple sequencing libraries/runs from the same sample. Spaces in sample names are automatically converted to underscores (`_`). |
| `fastq_1`      | Absolute path to FASTQ file for Illumina short reads 1. File has to be gzipped and have the extension ".fastq.gz" or ".fq.gz".                                                         |
| `fastq_2`      | Absolute path to FASTQ file for Illumina short reads 2. File has to be gzipped and have the extension ".fastq.gz" or ".fq.gz".                                                         |
| `strandedness` | Sample strandness                                                                                                                                                                      |

It is also possible to include multiple runs of the same sample in a samplesheet. For example, when you have re-sequenced the same sample more than once to increase sequencing depth. In that case, the `sample` identifiers in the samplesheet have to be the same. The pipeline will align the raw read/read-pairs independently before merging the alignments belonging to the same sample. Below is an example for a trio with the proband sequenced across two lanes:

| case  | sample       | fastq_1                          | fastq_2                          | strandedness |
| ----- | ------------ | -------------------------------- | -------------------------------- | ------------ |
| fam_1 | CONTROL_REP1 | AEG588A1_S1_L002_R1_001.fastq.gz | AEG588A1_S1_L002_R2_001.fastq.gz | reverse      |
| fam_1 | CONTROL_REP2 | AEG588A2_S1_L003_R1_001.fastq.gz | AEG588A2_S1_L003_R2_001.fastq.gz | reverse      |
| fam_1 | PATIENT_1    | AEG588A3_S1_L001_R1_001.fastq.gz | AEG588A3_S1_L001_R2_001.fastq.gz | reverse      |
| fam_1 | PATIENT_1    | AEG588A3_S1_L002_R1_001.fastq.gz | AEG588A3_S1_L002_R2_001.fastq.gz | reverse      |

If you would like to see more examples of what a typical samplesheet looks like for a duo, follow this links, [sample_sheet](https://github.com/genomic-medicine-sweden/tomte/blob/master/test_data/samplesheet_chr21.csv)

#### Reference files and parameters

In genomic-medicine-sweden/tomte, references can be supplied using parameters. We have also introduced the possiblility of using the ´--igenomes_base´ parameter to point to a path where genome specific refrence files are placed (fasta, fai, gtf, star_index, salmon_index, subsample_bed). To make sure that the names of the reference files match those in your directory, check [igenomes.config](https://github.com/genomic-medicine-sweden/tomte/blob/master/conf/igenomes.config).

Note that the pipeline is modular in architecture. It offers you the flexibility to choose between different tools. For example, you can call SNVs either with BCFtools or with GATK. You also have the option to turn off sections of the pipeline if you do not want to run them. For example, drop aberrant expression module can be turned off by setting `--switch_drop_ae FALSE`. This flexibility means that in any given analysis run, a combination of tools included in the pipeline will not be executed. So the pipeline is written in a way that can account for these differences while working with reference parameters. If a tool is not going to be executed during the course of a run, parameters used only by that tool need not be provided. For example, if you are not running DROP aberrant splicing, you do not need to provide `--reference_drop_splice_folder`.

genomic-medicine-sweden/tomte consists of several tools used for various purposes. For convenience, we have grouped those tools under the following categories:

1. Alignment and pseudo quantification (STAR & Salmon)
2. Junction track and bigwig
3. Subsample_region (Samtools)
4. Variant calling - SNV (BCFTools or GATK's GermlineCNVCaller)
5. SNV annotation (ensembl VEP)
6. Stringtie & gffcompare
7. DROP

> We have only listed the groups that require at least one input from the user. For example, the pipeline also runs WigToBigWig, but it does not require any input other than the bam files passed by the pipeline. Hence, it is not mentioned in the list above. To know more about the tools used in the pipeline check the [README](../README.md).

The mandatory and optional parameters for each category are tabulated below.

> Alignment, QC stats, repeat expansions, SNV variant calling and ensembl VEP are run by default. Hence, the mandatory parameters used by those features will always have to be provided to the pipeline.

##### 1. Alignment

| Mandatory | Optional                       |
| --------- | ------------------------------ |
| fasta     | fasta_fai<sup>1</sup>          |
| gtf       | sequence_dict<sup>1</sup>      |
|           | salmon_index<sup>1</sup>       |
|           | star_index<sup>1</sup>         |
|           | transcript_fasta<sup>1</sup>   |
|           | genome<sup>2</sup>             |
|           | platform<sup>3</sup>           |
|           | min_trimmed_length<sup>4</sup> |
|           | star_two_pass_mode<sup>4</sup> |

<sup>1</sup> If the parameter is not provided by the user, it will be generated from the fasta and gtf files.<br />
<sup>2</sup> If it is not provided by the user, the default value is GRCh38.<br />
<sup>3</sup> If it is not provided by the user, the default value is illumina.<br />
<sup>4</sup> If it is not provided by the user, the default value is 40.<br />
<sup>5</sup> If it is not provided by the user, the default value is Basic.

##### 2. Junction track and bigwig

| Mandatory | Optional                         |
| --------- | -------------------------------- |
|           | switch_build_tracks <sup>1</sup> |

<sup>1</sup> If it is not provided by the user, the default value is true

##### 3. Subsample region

| Mandatory     | Optional                            |
| ------------- | ----------------------------------- |
| subsample_bed | switch_subsample_region<sup>1</sup> |
|               | seed_frac<sup>2</sup>               |

<sup>1</sup> If it is not provided by the user, the default value is true
<sup>2</sup> If it is not provided by the user, the default value is 0.001

##### 4. Variant calling - SNV

| Mandatory | Optional                         |
| --------- | -------------------------------- |
|           | variant_caller<sup>1</sup>       |
|           | bcftools_caller_mode<sup>2</sup> |

<sup>1</sup> If it is not provided by the user, the default value is bcftools<br />
<sup>2</sup> If it is not provided by the user, the default value is multiallelic

#### 5. SNV annotation (ensembl VEP)

| Mandatory                    | Optional               |
| ---------------------------- | ---------------------- |
| vep_plugin_files<sup>1</sup> | switch_vep<sup>2</sup> |
|                              | vep_cache<sup>3</sup>  |
|                              | vep_cache_version      |
|                              | vep_filters            |

<sup>1</sup> VEP caches can be downloaded [here](https://www.ensembl.org/info/docs/tools/vep/script/vep_cache.html#cache). VEP plugins may be installed in the cache directory, and the plugin pLI is mandatory to install. To supply files required by VEP plugins, use `vep_plugin_files` parameter. See example cache [here](https://raw.githubusercontent.com/nf-core/test-datasets/raredisease/reference/vep_cache_and_plugins.tar.gz).<br />
<sup>2</sup> If it is not provided by the user, the default value is true<br />
<sup>3</sup> If it is not provided by the user, the default value is 110, supported values are 107 and 110 <br />

#### 6. Stringtie & gffcompare

| Mandatory | Optional                     |
| --------- | ---------------------------- |
| fasta     | switch_stringtie<sup>1</sup> |
| gtf       |                              |

<sup>1</sup> If it is not provided by the user, the default value is true

#### 7. DROP

DROP - aberrant expression

| Mandatory                             | Optional                          |
| ------------------------------------- | --------------------------------- |
| reference_drop_annot_file<sup>1</sup> | switch_drop_ae<sup>2</sup>        |
| reference_drop_count_file             | drop_group_samples_ae<sup>3</sup> |
| fasta                                 | drop_padjcutoff_ae<sup>4</sup>    |
| gtf                                   | drop_zscorecutoff<sup>5</sup>     |
|                                       | gene_panel_clinical_filter        |
|                                       | switch_downsample<sup>6</sup>     |
|                                       | num_reads<sup>7</sup>             |
|                                       | genome<sup>8</sup>                |

<sup>1</sup> To get more information on how to format it, see below<br />
<sup>2</sup> If it is not provided by the user, the default value is true<br />
<sup>3</sup> If it is not provided by the user, the default value is outrider<br />
<sup>4</sup> If it is not provided by the user, the default value is 0.05<br />
<sup>5</sup> If it is not provided by the user, the default value is 0<br />
<sup>6</sup> If it is not provided by the user, the default value is true<br />
<sup>7</sup> If it is not provided by the user, the default value is 120000000<br />
<sup>8</sup> If it is not provided by the user, the default value is GRCh38

DROP - aberrant splicing

| Mandatory                             | Optional                          |
| ------------------------------------- | --------------------------------- |
| reference_drop_annot_file<sup>1</sup> | switch_drop_as<sup>2</sup>        |
| reference_drop_splice_folder          | drop_group_samples_as<sup>3</sup> |
|                                       | drop_padjcutoff_as<sup>4</sup>    |
|                                       | gene_panel_clinical_filter        |
|                                       | switch_downsample<sup>5</sup>     |
|                                       | num_reads<sup>6</sup>             |
|                                       | genome<sup>7</sup>                |

<sup>1</sup> To get more information on how to format it, see below<br />
<sup>2</sup> If it is not provided by the user, the default value is true<br />
<sup>3</sup> If it is not provided by the user, the default value is fraser<br />
<sup>4</sup> If it is not provided by the user, the default value is 0.1<br />
<sup>5</sup> If it is not provided by the user, the default value is true<br />
<sup>6</sup> If it is not provided by the user, the default value is 120000000<br />
<sup>7</sup> If it is not provided by the user, the default value is GRCh38

##### Preparing input for DROP

If you want to run [DROP](https://github.com/gagneurlab/drop) aberrant expression or aberrant splicing you have to provide reference counts, splice counts and a sample sheet. The sample sheet should contain the columns as those in the [test sample annotation](../test_data/drop_data/sampleAnnotation.tsv), you do not need to include the samples you are running through the pipeline in the sample sheet.

To obtain the gene counts and splice counts you will have to download the counts from one of the [available databases](https://github.com/gagneurlab/drop#datasets) or run drop locally with your own samples. If you choose the second option, you should start by runnig the module(s) you want to export counts for. Afterwards, you need to run the exportCounts module. Make sure that your config has only the modules you want to export and have already run as <run: true> , that only existing groups are mentioned in the config, and that exportCounts excludGroups is null or contains a group of samples you want to exclude. Finally, run:

```console
snakemake exportCounts --cores 1
```

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

:::warning
Do not use `-c <file>` to specify parameters as this will result in errors. Custom config files specified with `-c` must only be used for [tuning process resource specifications](https://nf-co.re/docs/usage/configuration#tuning-workflow-resources), other infrastructural tweaks (such as output directories), or module arguments (args).
:::

The above pipeline run specified with a params file in yaml format:

```bash
nextflow run genomic-medicine-sweden/tomte -profile docker -params-file params.yaml
```

with `params.yaml` containing:

```yaml
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

It is a good idea to specify a pipeline version when running the pipeline on your data. This ensures that a specific version of the pipeline code and software are used when you run your pipeline. If you keep using the same tag, you'll be running the same version of the pipeline, even if there have been changes to the code since.

First, go to the [genomic-medicine-sweden/tomte releases page](https://github.com/genomic-medicine-sweden/tomte/releases) and find the latest pipeline version - numeric only (eg. `1.3.1`). Then specify this when running the pipeline with `-r` (one hyphen) - eg. `-r 1.3.1`. Of course, you can switch to another version by changing the number after the `-r` flag.

This version number will be logged in reports when you run the pipeline, so that you'll know what you used when you look back in the future. For example, at the bottom of the MultiQC reports.

To further assist in reproducbility, you can use share and re-use [parameter files](#running-the-pipeline) to repeat pipeline runs with the same settings without having to write out a command with every single parameter.

:::tip
If you wish to share such profile (such as upload as supplementary material for academic publications), make sure to NOT include cluster specific paths to files, nor institutional specific profiles.
:::

## Core Nextflow arguments

:::note
These options are part of Nextflow and use a _single_ hyphen (pipeline parameters use a double-hyphen).
:::

### `-profile`

Use this parameter to choose a configuration profile. Profiles can give configuration presets for different compute environments.

Several generic profiles are bundled with the pipeline which instruct the pipeline to use software packaged using different methods (Docker, Singularity, Podman, Shifter, Charliecloud, Apptainer, Conda) - see below.

:::info
We highly recommend the use of Docker or Singularity containers for full pipeline reproducibility, however when this is not possible, Conda is also supported.
:::

The pipeline also dynamically loads configurations from [https://github.com/nf-core/configs](https://github.com/nf-core/configs) when it runs, making multiple config profiles for various institutional clusters available at run time. For more information and to see if your system is available in these configs please see the [nf-core/configs documentation](https://github.com/nf-core/configs#documentation).

Note that multiple profiles can be loaded, for example: `-profile test,docker` - the order of arguments is important!
They are loaded in sequence, so later profiles can overwrite earlier profiles.

If `-profile` is not specified, the pipeline will run locally and expect all software to be installed and available on the `PATH`. This is _not_ recommended, since it can lead to different results on different machines dependent on the computer enviroment.

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
- `conda`
  - A generic configuration profile to be used with [Conda](https://conda.io/docs/). Please only use Conda as a last resort i.e. when it's not possible to run the pipeline with Docker, Singularity, Podman, Shifter, Charliecloud, or Apptainer.

### `-resume`

Specify this when restarting a pipeline. Nextflow will use cached results from any pipeline steps where the inputs are the same, continuing from where it got to previously. For input to be considered the same, not only the names must be identical but the files' contents as well. For more info about this parameter, see [this blog post](https://www.nextflow.io/blog/2019/demystifying-nextflow-resume.html).

You can also supply a run name to resume a specific run: `-resume [run-name]`. Use the `nextflow log` command to show previous run names.

### `-c`

Specify the path to a specific config file (this is a core Nextflow command). See the [nf-core website documentation](https://nf-co.re/usage/configuration) for more information.

## Custom configuration

### Resource requests

Whilst the default requirements set within the pipeline will hopefully work for most people and with most input data, you may find that you want to customise the compute resources that the pipeline requests. Each step in the pipeline has a default set of requirements for number of CPUs, memory and time. For most of the steps in the pipeline, if the job exits with any of the error codes specified [here](https://github.com/nf-core/rnaseq/blob/4c27ef5610c87db00c3c5a3eed10b1d161abf575/conf/base.config#L18) it will automatically be resubmitted with higher requests (2 x original, then 3 x original). If it still fails after the third attempt then the pipeline execution is stopped.

To change the resource requests, please see the [max resources](https://nf-co.re/docs/usage/configuration#max-resources) and [tuning workflow resources](https://nf-co.re/docs/usage/configuration#tuning-workflow-resources) section of the nf-core website.

### Custom Containers

In some cases you may wish to change which container or conda environment a step of the pipeline uses for a particular tool. By default nf-core pipelines use containers and software from the [biocontainers](https://biocontainers.pro/) or [bioconda](https://bioconda.github.io/) projects. However in some cases the pipeline specified version maybe out of date.

To use a different container from the default container or conda environment specified in a pipeline, please see the [updating tool versions](https://nf-co.re/docs/usage/configuration#updating-tool-versions) section of the nf-core website.

### Custom Tool Arguments

A pipeline might not always support every possible argument or option of a particular tool used in pipeline. Fortunately, nf-core pipelines provide some freedom to users to insert additional parameters that the pipeline does not include by default.

To learn how to provide additional arguments to a particular tool of the pipeline, please see the [customising tool arguments](https://nf-co.re/docs/usage/configuration#customising-tool-arguments) section of the nf-core website.

### nf-core/configs

In most cases, you will only need to create a custom config as a one-off but if you and others within your organisation are likely to be running nf-core pipelines regularly and need to use the same settings regularly it may be a good idea to request that your custom config file is uploaded to the `nf-core/configs` git repository. Before you do this please can you test that the config file works with your pipeline of choice using the `-c` parameter. You can then create a pull request to the `nf-core/configs` repository with the addition of your config file, associated documentation file (see examples in [`nf-core/configs/docs`](https://github.com/nf-core/configs/tree/master/docs)), and amending [`nfcore_custom.config`](https://github.com/nf-core/configs/blob/master/nfcore_custom.config) to include your custom profile.

See the main [Nextflow documentation](https://www.nextflow.io/docs/latest/config.html) for more information about creating your own configuration files.

If you have any questions or issues please send us a message on [Slack](https://nf-co.re/join/slack) on the [`#configs` channel](https://nfcore.slack.com/channels/configs).

## Azure Resource Requests

To be used with the `azurebatch` profile by specifying the `-profile azurebatch`.
We recommend providing a compute `params.vm_type` of `Standard_D16_v3` VMs by default but these options can be changed if required.

Note that the choice of VM size depends on your quota and the overall workload during the analysis.
For a thorough list, please refer the [Azure Sizes for virtual machines in Azure](https://docs.microsoft.com/en-us/azure/virtual-machines/sizes).

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
