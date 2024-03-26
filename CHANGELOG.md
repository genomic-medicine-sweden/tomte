# genomic-medicine-sweden/tomte: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

### `Added`

- Added automatic tests to test the pipeline with all switches set to false [#100](https://github.com/genomic-medicine-sweden/tomte/pull/100)

### `Fixed`

- Subsample and downsample switches [#97](https://github.com/genomic-medicine-sweden/tomte/pull/97)
- Now all reference files come with meta to avoid confusion [#101](https://github.com/genomic-medicine-sweden/tomte/pull/101)
- GATK4_ASEREADCOUNTER and GATK4_SPLITNCIGARREADS have been updated [#101](https://github.com/genomic-medicine-sweden/tomte/pull/101)

## 1.1.0 - Rudolph [2024-03-11]

Release of genomic-medicine-sweden/tomte, created with the [nf-core](https://nf-co.re/) template.

### `Added`

- switch_vep, switch_build_tracks and switch_stringtie to make the pipeline more versatile [#61](https://github.com/genomic-medicine-sweden/tomte/pull/61)
- Updated template and nf-tools [#65](https://github.com/genomic-medicine-sweden/tomte/pull/65)
- Use `nf-validation` plugin for parameter and samplesheet validation [#66](https://github.com/genomic-medicine-sweden/tomte/pull/66)
- Installed the nf-core version of ensemblvep/vep module [#67](https://github.com/genomic-medicine-sweden/tomte/pull/67)
- A new parameter `vep_plugin_files` to supply files required by vep plugins [#67](https://github.com/genomic-medicine-sweden/tomte/pull/67)
- The possibility of using `igenomes_base` to point to a path where genome specific reference files are placed (fasta, fai, gtf, star_index, salmon_index, subsample_bed) [#76](https://github.com/genomic-medicine-sweden/tomte/pull/76)
- Merging of case's vcf files [#80](https://github.com/genomic-medicine-sweden/tomte/pull/80)
- Reference list to MultiQC report [#88](https://github.com/genomic-medicine-sweden/tomte/pull/88)
- Added module to calculate insert size and added results to MultiQC report [#90](https://github.com/genomic-medicine-sweden/tomte/pull/90)

### `Fixed`

- Renamed the other switches (subsample_region_switch, downsample_switch, run_drop_ae_switch and run_drop_as_switch) so that they all start with switch\* (switch_subsample_region, switch_downsample, switch_drop_ae and switch_drop_as) [#61](https://github.com/genomic-medicine-sweden/tomte/pull/61)
- Separated modules.config into smaller configs [#61](https://github.com/genomic-medicine-sweden/tomte/pull/61)
- Missing fasta_fai channel when fai file is given [#63](https://github.com/genomic-medicine-sweden/tomte/pull/63)
- DROP output file columns, removing duplicate column and adding same ids to both AE and AS [#68](https://github.com/genomic-medicine-sweden/tomte/pull/68)
- Patch tools update and case ID parsing [#71](https://github.com/genomic-medicine-sweden/tomte/pull/71)
- Naming of DROP output files [#72](https://github.com/genomic-medicine-sweden/tomte/pull/72)
- VEP plugin schema to allow for directories [#74](https://github.com/genomic-medicine-sweden/tomte/pull/74)
- Made params.platform into a channel [#75](https://github.com/genomic-medicine-sweden/tomte/pull/75)
- Changed name of salmon's quant.nf to include sample id [#78](https://github.com/genomic-medicine-sweden/tomte/pull/78)
- Shortened name of DROP output files [#79](https://github.com/genomic-medicine-sweden/tomte/pull/79)
- Merging of vcfs has been moved to after bootstrapAnn [#81](https://github.com/genomic-medicine-sweden/tomte/pull/81)
- Substituted bgzip and tabix modules by bgzip_tabix module [#85](https://github.com/genomic-medicine-sweden/tomte/pull/85)
- Updated module input channels in the GATK variant calling subworkflow [#89](https://github.com/genomic-medicine-sweden/tomte/pull/89)

### `Dependencies`

### `Deprecated`

## 1.0.0 - Nisse [2023-11-06]

### `Added`

- Trim reads with FASTP
- Read mapping with STAR
- Transcript quantification with Salmon
- Output junction tracks
- Output bigwig
- Choice to subsample overrepresented regions with Samtools
- Choice to downsample number of reads with Samtools
- Detection of aberrant expression with DROP
- Detection of aberrant splicing with DROP
- Filter aberrant expression and aberrant splicing results
- Guided transcript assembly with StringTie
- Filtering results of guided transcript assembly with GffCompare
- SNVs calling with GATK or BCFtools Mpileups
- Allele Specific Read Counter with ASEReadCounter
- Assess allelic imbalance with BootstrapAnn
- Annotation with VEP
- Alignment QC with Picard CollectRnaSeqMetrics
- Present QCs with MultiQC
