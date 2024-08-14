# genomic-medicine-sweden/tomte: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

# X.X.X - [XXXX-XX-XX]

### `Added`

- Fasta, gtf, vep cache and plugins can now be downloaded automatically by the pipeline if they are not provided by the user [#149](https://github.com/genomic-medicine-sweden/tomte/pull/149)
- Added `--gencode_annotation_version`, the version of the gencode reference version to download if fasta or gtf is not provided [#149](https://github.com/genomic-medicine-sweden/tomte/pull/149)
- Added the possibility to provide `--vep_refs_download`, a comma separated csv determining the vep references that should be downloaded (excluding gnomad ones) alongside with a switch `--skip_download_vep` for the vep reference download in general and `--skip_download_gnomad` for gnomad in particular [#149](https://github.com/genomic-medicine-sweden/tomte/pull/149)

### `Fixed`

### `Parameters`

| Old parameter | New parameter                  |
| ------------- | ------------------------------ |
|               | `--gencode_annotation_version` |
|               | `--vep_refs_download`          |
|               | `--skip_download_vep`          |
|               | `--skip_download_gnomad`       |

:::note Parameter has been updated if both old and new parameter information is present. Parameter has been added if just thenew parameter information is present. Parameter has been removed if new parameter information isn't present. :::

### `Changed`

## 2.1.0 - Elf [2024-06-26]

### `Added`

- Installed bcftools/norm [#127](https://github.com/genomic-medicine-sweden/tomte/pull/127)
- Installed bcftools/annotate [#127](https://github.com/genomic-medicine-sweden/tomte/pull/127)

### `Fixed`

- One line per call on vcf to make vcf suitable for Scout [#127](https://github.com/genomic-medicine-sweden/tomte/pull/127)
- Added variant caller to vcf to make vcf suitable for Scout [#127](https://github.com/genomic-medicine-sweden/tomte/pull/127)
- Normalised calls in vcf [#127](https://github.com/genomic-medicine-sweden/tomte/pull/127)
- Removed regions parameter from ASEReadCounter to obtain ASE for all regions [#129](https://github.com/genomic-medicine-sweden/tomte/pull/129)
- Updated container for GET_CHROM_SIZES [#132](https://github.com/genomic-medicine-sweden/tomte/pull/132)
- Updated container for RENAME_FILES [#132](https://github.com/genomic-medicine-sweden/tomte/pull/132)
- Fixed GATK's ASEReadCounter by adding bcftools norm to avoid having duplicated positions in vcf [#137](https://github.com/genomic-medicine-sweden/tomte/pull/137)

### `Parameters`

### `Changed`

- Updated template to v2.14.1 [#123](https://github.com/genomic-medicine-sweden/tomte/pull/123)
- Changed DROP output column names to camel case [#132](https://github.com/genomic-medicine-sweden/tomte/pull/132)
- Updated bcftools view and norm to be able to create index within the actual module [#137](https://github.com/genomic-medicine-sweden/tomte/pull/137)
- Updated bcftools annot to be able to create index within the actual module [#140](https://github.com/genomic-medicine-sweden/tomte/pull/140)
- Updated modules bcftools/merge, bcftools/mpileup, bcftools/stats, ensemblvep/vep, fastp, gawk, multiqc, salmon/quant, samtools/faidx, samtools/index and samtools/view [#141](https://github.com/genomic-medicine-sweden/tomte/pull/141)
- Removed unused modules bcftools/index and tabix/bgzip [#141](https://github.com/genomic-medicine-sweden/tomte/pull/141)

| Tool             | Old version | New version |
| ---------------- | ----------- | ----------- |
| bcftools/index   | 1.18        |             |
| bcftools/merge   | 1.18        | 1.20        |
| bcftools/mpileup | 1.18        | 1.20        |
| bcftools/norm    | 1.18        | 1.20        |
| bcftools/stats   | 1.18        | 1.20        |
| bcftools/view    | 1.18        | 1.20        |
| fastp            | 0.23.4      | 0.23.4      |
| gawk             | 5.1.0       | 5.3.0       |
| multiqc          | 1.21        | 1.22.3      |
| salmon/quant     | 1.10.1      | 1.10.1      |
| samtools/faidx   | 1.19.2      | 1.20        |
| samtools/index   | 1.19.2      | 1.20        |
| samtools/view    | 1.19.2      | 1.20        |
| tabix/bgzip      | 1.19.1      |             |

:::note Version has been updated if both old and new version information is present. Version has been added if just the new version information is present. Version has been removed if new version information isn't present. :::

## 2.0.1 - Grinch [2024-04-25]

### `Added`

### `Fixed`

- Vep annotated research results will be published [#115](https://github.com/genomic-medicine-sweden/tomte/pull/115)

### `Parameters`

## 2.0.0 - Santa [2024-04-19]

### `Added`

- Added automatic tests to test the pipeline with all switches set to false [#100](https://github.com/genomic-medicine-sweden/tomte/pull/100)
- Added better documentation on subworkflow input [#101](https://github.com/genomic-medicine-sweden/tomte/pull/101)
- Added option to add extra arguments to DROP aberrant expression and aberrant splicing [#104](https://github.com/genomic-medicine-sweden/tomte/pull/104)
- Added a function to branch references into compressed/uncompressed [#107](https://github.com/genomic-medicine-sweden/tomte/pull/107)
- Added nf-core modules gawk and filter vep to create a clinical vcf [#109](https://github.com/genomic-medicine-sweden/tomte/pull/109)

### `Fixed`

- Subsample and downsample switches [#97](https://github.com/genomic-medicine-sweden/tomte/pull/97)
- Now all reference files come with meta to avoid confusion [#101](https://github.com/genomic-medicine-sweden/tomte/pull/101)
- GATK4_ASEREADCOUNTER and GATK4_SPLITNCIGARREADS have been updated [#101](https://github.com/genomic-medicine-sweden/tomte/pull/101)
- Updated GATK4_ASEREADCOUNTER, now bam and vcf will be given as one channel [#103](https://github.com/genomic-medicine-sweden/tomte/pull/103)
- Prepare reference subworkflow has been reformated and simplified [#105](https://github.com/genomic-medicine-sweden/tomte/pull/105)
- FastQC have been updated to correctly allocate memory [#106](https://github.com/genomic-medicine-sweden/tomte/pull/106)
- vep_filters is now extracted from gene_panel_clinical_filter [#109](https://github.com/genomic-medicine-sweden/tomte/pull/109)
- Updated modules bcftools/stats, ensemblvep/vep, fastp, gatk4/bedtointervallist, samtools/faidx [#110](https://github.com/genomic-medicine-sweden/tomte/pull/110)

### `Parameters`

- Removed `--vep_filters`, it will now be automatically extracted from the `--gene_panel_clinical_filter`[#109](https://github.com/genomic-medicine-sweden/tomte/pull/109)
  | Old parameter | New parameter |
  | --------------- | ------------- |
  | `--vep_filters` | |

- Updated parameter names to make their use easier and more clear, changing the names from `switch` to `skip` and their default value from `true` to `false` [#108](https://github.com/genomic-medicine-sweden/tomte/pull/108)

| Old parameter               | New parameter             |
| --------------------------- | ------------------------- |
| `--switch_subsample_region` | `--skip_subsample_region` |
| `--switch_downsample`       | `--skip_downsample`       |
| `--switch_build_tracks`     | `--skip_build_tracks`     |
| `--switch_stringtie`        | `--skip_stringtie`        |
| `--switch_vep`              | `--skip_vep`              |
| `--switch_drop_ae`          | `--skip_drop_ae`          |
| `--switch_drop_as`          | `--skip_drop_as`          |

:::note Parameter has been updated if both old and new parameter information is present. Parameter has been added if just the new parameter information is present. Parameter has been removed if new parameter information isn't present. :::

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
