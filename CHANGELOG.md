# genomic-medicine-sweden/tomte: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## 1.1.0 - Rudolph [xxxx-xx-xx]

Initial release of genomic-medicine-sweden/tomte, created with the [nf-core](https://nf-co.re/) template.

### `Added`

- switch_vep, switch_build_tracks and switch_stringtie to make the pipeline more versatile [#61](https://github.com/genomic-medicine-sweden/tomte/pull/61)
- Updated template and nf-tools [#65](https://github.com/genomic-medicine-sweden/tomte/pull/65)
- Use `nf-validation` plugin for parameter and samplesheet validation [#66](https://github.com/genomic-medicine-sweden/tomte/pull/66)
- Installed the nf-core version of ensemblvep/vep module [#67](https://github.com/genomic-medicine-sweden/tomte/pull/67)
- A new parameter `vep_plugin_files` to supply files required by vep plugins [#67](https://github.com/genomic-medicine-sweden/tomte/pull/67)

### `Fixed`

- Renamed the other switches (subsample_region_switch, downsample_switch, run_drop_ae_switch and run_drop_as_switch) so that they all start with switch\* (switch_subsample_region, switch_downsample, switch_drop_ae and switch_drop_as)
- Separated modules.config into smaller configs
- Patch tools update [#71](https://github.com/genomic-medicine-sweden/tomte/pull/71)

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
