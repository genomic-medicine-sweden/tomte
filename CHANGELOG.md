# genomic-medicine-sweden/tomte: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## X.X.X - TBD [XXXX-XX-XX]

### `Added`

- Added the option of starting from a bam or cram file [#193](https://github.com/genomic-medicine-sweden/tomte/pull/193)
- Optionally run Peddy for per-sample sex- and heterozygosity checks [#190](https://github.com/genomic-medicine-sweden/tomte/pull/190)
- Optionally calculate percentage mapping to hemoglobin genes (or any other set of genes provided) [#190](https://github.com/genomic-medicine-sweden/tomte/pull/190)
- Added the option of providing sex as 0, 1, or 2 as in the raredisease pipeline [#192](https://github.com/genomic-medicine-sweden/tomte/pull/192)
- Nf-test for pipeline [#207](https://github.com/genomic-medicine-sweden/tomte/pull/207)
- Output channels to Tomte, allowing it to run as part of a larger pipeline [#206](https://github.com/genomic-medicine-sweden/tomte/pull/206)
- Nf-test for alignment subworkflow [#209](https://github.com/genomic-medicine-sweden/tomte/pull/209)
- Nf-test for prepare_references subworkflow [#214](https://github.com/genomic-medicine-sweden/tomte/pull/214)
- Nf-test for analyse_transcripts subworkflow [#215](https://github.com/genomic-medicine-sweden/tomte/pull/215)
- Nf-test for call_variants subworkflow [#218](https://github.com/genomic-medicine-sweden/tomte/pull/218)
- Nf-test for bam_qc subworkflow [#219](https://github.com/genomic-medicine-sweden/tomte/pull/219)
- Nf-test for analyse_transcripts subworkflow allele_specific_calling [#220](https://github.com/genomic-medicine-sweden/tomte/pull/220)
- Nf-test for annotate_snv subworkflow [#221](https://github.com/genomic-medicine-sweden/tomte/pull/221)
- Nf-test for igv_tracks subworkflow [#222](https://github.com/genomic-medicine-sweden/tomte/pull/222)
- Nf-test for download_references subworkflow [#223](https://github.com/genomic-medicine-sweden/tomte/pull/223)
- Nf-test for bootstrapann, create_pedigree_file, estimate_hb_perc, get_chrom_sizes, get_rrna_transcripts,junction_track, rename_files, rna_downsample, and rna_subsample_region modules [#224](https://github.com/genomic-medicine-sweden/tomte/pull/224)

### `Parameters`

| Old parameter | New parameter              |
| ------------- | -------------------------- |
|               | `--skip_peddy`             |
|               | `--skip_calculate_hb_frac` |
|               | `--hb_genes`               |
|               | `--trace_report_suffix`    |

### `Changed`

- Updated nf-core/tools template to v3.1.2 [#204](https://github.com/genomic-medicine-sweden/tomte/pull/204)
- Updated nf-core/tools template to v3.2.0 [#205](https://github.com/genomic-medicine-sweden/tomte/pull/205)
- Updated multiqc module [#205](https://github.com/genomic-medicine-sweden/tomte/pull/205)
- Updated gatk4/splitncigarreads module [#206](https://github.com/genomic-medicine-sweden/tomte/pull/206)
- Updated local gffread module to nf-core one [#227](https://github.com/genomic-medicine-sweden/tomte/pull/227)

| Tool                   | Old version | New version |
| ---------------------- | ----------- | ----------- |
| multiqc                | 1.25.1      | 1.27        |
| gatk4/splitncigarreads | 4.5.0.0     | 4.6.1.0     |
| gffread                | 0.12.1      | 0.12.7      |

### `Fixed`

- DROP expression outliers previously returned max 20 lowest p-value entries. Now it returns those classified as having aberrant expression first, and then fills up to 20 from those with lowest p-values. [#206](https://github.com/genomic-medicine-sweden/tomte/pull/206)
- GENCODE_DOWNLOAD stub crashes by removing double versions.yml write. [#206](https://github.com/genomic-medicine-sweden/tomte/pull/206)
- Default peddy parent ID to "0" if not present. [#206](https://github.com/genomic-medicine-sweden/tomte/pull/206)
- BUILD_VEP_CACHE, VEP_GNOMAD_DOWNLOAD and WGET_DOWNLOAD stubs crash, fixed by removing double versions.yml write. [#223](https://github.com/genomic-medicine-sweden/tomte/pull/223)

## 3.0.1 - Snowman [2025-03-26]

### `Fixed`

- Making a release to fix to the bug solved in [#206](https://github.com/genomic-medicine-sweden/tomte/pull/206). DROP expression outliers previously returned max 20 lowest p-value entries. Now it returns those classified as having aberrant expression first, and then fills up to 20 from those with lowest p-values. [#210](https://github.com/genomic-medicine-sweden/tomte/pull/210)

## 3.0.0 - Three Kings [2024-11-18]

### `Added`

- Functionality to create DROP databases and to add samples to existing ones [#147](https://github.com/genomic-medicine-sweden/tomte/pull/147)
- A switch `--skip_variant_calling` for variant calling [#169](https://github.com/genomic-medicine-sweden/tomte/pull/169)
- Functionality to output DROP databases in references folder with a working sample annotation sheet [#172](https://github.com/genomic-medicine-sweden/tomte/pull/172)
- Added optional sex info col to samplesheet, used in DROP[#168](https://github.com/genomic-medicine-sweden/tomte/pull/168)
- Added more documentation regarding DROP [#178](https://github.com/genomic-medicine-sweden/tomte/pull/178)

### `Fixed`

- Versions for all modules involving drop will now be outputed in version.yml and multiqc file [#174](https://github.com/genomic-medicine-sweden/tomte/pull/174)
- Fixed bug when running variant calling with gatk [#182](https://github.com/genomic-medicine-sweden/tomte/pull/182)

### `Parameters`

| Old parameter                    | New parameter            |
| -------------------------------- | ------------------------ |
| `--max_cpus`                     |                          |
| `--max_memory`                   |                          |
| `--max_time`                     |                          |
| `--validationShowHiddenParams`   |                          |
| `--validationSkipDuplicateCheck` |                          |
| `--validationS3PathCheck`        |                          |
| `--monochromeLogs`               | `--monochrome_logs`      |
|                                  | `--skip_variant_calling` |

### `Changed`

- Updated modules ensemblvep/filtervep, ensemblvep/vep [#159](https://github.com/genomic-medicine-sweden/tomte/pull/159)
- Updated gencode version from 37 to 46 [#159](https://github.com/genomic-medicine-sweden/tomte/pull/159)
- Updated modules using drop drop_config_runAE, drop_config_runAS, drop_sample_annot, and drop_filter_results [#147](https://github.com/genomic-medicine-sweden/tomte/pull/147)
- Updated nf-core/tools template to v3.0.2 [#167](https://github.com/genomic-medicine-sweden/tomte/pull/167)
- Updated multiqc version to 1.25.1 [#167](https://github.com/genomic-medicine-sweden/tomte/pull/167)
- Updated modules bcftools/annotate, bcftools/merge, bcftools/mpileup, bcftools/norm, bcftools/stats, bcftools/view, cat/fastq, ensemblvep/filtervep, ensemblvep/vep, fastp, gatk4/asereadcounter, gatk4/bedtointervallist, gatk4/createsequencedictionary, gatk4/haplotypecaller, gatk4/splitncigarreads, gatk4/variantfiltration, gawk, gffcompare, gunzip, picard/collectinsertsizemetrics, picard/collectrnaseqmetrics, salmon/index, salmon/quant, samtools/faidx, samtools/index, samtools/view, star/align, star/genomegenerate, stringtie/stringtie, tabix/bgziptabix, tabix/tabix, ucsc/wigtobigwig, untar [#177](https://github.com/genomic-medicine-sweden/tomte/pull/177)
- Changed how variant caller is added to the vcf, it is now done using the local module add_found_in_tag [#184](https://github.com/genomic-medicine-sweden/tomte/pull/184)

| Tool                           | Old version | New version |
| ------------------------------ | ----------- | ----------- |
| ensemblvep/filtervep           | 110         | 113         |
| ensemblvep/vep                 | 110         | 113         |
| DROP                           | 1.3.3       | 1.4.0       |
| multiqc                        | 1.21        | 1.25.1      |
| cat/fastq                      | 8.30        | 9.5         |
| picard/collectinsertsizemtrics | 3.2.0       | 3.3.0       |
| salmon/index                   | 1.10.1      | 1.10.3      |
| salmon/quant                   | 1.10.1      | 1.10.3      |
| samtools/faidx                 | 1.20        | 1.21        |
| samtools/index                 | 1.20        | 1.21        |
| samtools/view                  | 1.20        | 1.21        |
| star/align                     | 2.7.10a     | 2.7.11b     |
| star/genomegenerate            | 2.7.10a     | 2.7.11b     |
| stringtie/stringtie            | 2.2.1       | 2.2.3       |

## 2.2.1 - Scrooge [2024-08-28]

### `Fixed`

- After an update, MultiQC was not outputing data for RnaSeqMetrics so an earlier version will be used [#156](https://github.com/genomic-medicine-sweden/tomte/pull/156)

### `Changed`

- Downgraded multiqc version [#156](https://github.com/genomic-medicine-sweden/tomte/pull/156)

| Tool    | Old version | New version |
| ------- | ----------- | ----------- |
| multiqc | 1.24.1      | 1.21        |

## 2.2.0 - TioDeNadal [2024-08-27]

### `Added`

- Fasta, gtf, vep cache and plugins can now be downloaded automatically by the pipeline if they are not provided by the user [#149](https://github.com/genomic-medicine-sweden/tomte/pull/149)
- Added `--gencode_annotation_version`, the version of the gencode reference version to download if fasta or gtf is not provided [#149](https://github.com/genomic-medicine-sweden/tomte/pull/149)
- Added the possibility to provide `--vep_refs_download`, a comma separated csv determining the vep references that should be downloaded (excluding gnomad ones) alongside with a switch `--skip_download_vep` for the vep reference download in general and `--skip_download_gnomad` for gnomad in particular [#149](https://github.com/genomic-medicine-sweden/tomte/pull/149)

### `Fixed`

- Input to BootstrapAnn is now supplied in a single channel. Previously they were supplied in separate channels, which could cause mix-ups if more than one sample was supplied [#151](https://github.com/genomic-medicine-sweden/tomte/pull/151)

### `Parameters`

| Old parameter | New parameter                  |
| ------------- | ------------------------------ |
|               | `--gencode_annotation_version` |
|               | `--vep_refs_download`          |
|               | `--skip_download_vep`          |
|               | `--skip_download_gnomad`       |

> [!NOTE]
> Parameter has been updated if both old and new parameter information is present.
> Parameter has been added if just thenew parameter information is present.
> Parameter has been removed if new parameter information isn't present.

### `Changed`

- Updated modules bcftools/annotate, bcftools/mpileup, bcftools/view, cat/fastq, ensemblvep/filtervep, fastp, fastqc, gatk4/haplotypecaller, gatk4/splitncigarreads, gunzip, multiqc, picard/collectrnaseqmetrics, samtools/index, star/align, star/genomegenerate, stringtie/stringtie, tabix/bgziptabix, tabix/tabix and untar [#153](https://github.com/genomic-medicine-sweden/tomte/pull/153)

| Tool                            | Old version | New version |
| ------------------------------- | ----------- | ----------- |
| gunzip                          | 20.04       | 22.04       |
| multiqc                         | 1.22.3      | 1.24.1      |
| picard/collectinsertsizemetrics | 3.1.1       | 3.2.0       |
| tabix/bgziptabix                | 1.19.1      | 1.20        |
| tabix/tabix                     | 1.19.1      | 1.20        |
| untar                           | 20.04       | 22.04       |

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
