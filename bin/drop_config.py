#!/usr/bin/env python3

import argparse
from pathlib import Path
import yaml


def get_config() -> dict:
    """Return default config dictionary used as input for DROP"""
    return yaml.safe_load(
        """
        projectTitle: "DROP: Detection of RNA Outliers Pipeline"
        root:             # root directory of all intermediate output
        htmlOutputPath:   # path for HTML rendered reports
        indexWithFolderName: true # whether the root base name should be part of the index name

        hpoFile: null  # if null, downloads it from webserver
        sampleAnnotation: # path to sample annotation (see documenation on how to create it)
        
        geneAnnotation:
            gtf: null
            # multiple annotations with custom names are possible
            # <annotation_name> : <path to gencode v29 annotation>
            # v37:  /path/to/gencode29.gtf.gz # example

        genomeAssembly: genome  # either hg19/hs37d5 or hg38/GRCh38
        exportCounts:
            # specify which gene annotations to include and which
            # groups to exclude when exporting counts
            geneAnnotations:
                - gtf
            excludeGroups:
                - null
        genome: # path to genome sequence in fasta format.
            # You can define multiple reference genomes in yaml format, ncbi: path_to_ncbi, ucsc: path_to_ucsc
            # the keywords that define the path should be in the GENOME column of the SA table

        exportCounts:
            # specify which gene annotations to include and which
            # groups to exclude when exporting counts
            geneAnnotations:
                - gtf
            excludeGroups:
                - null
        
        aberrantExpression:
            run: false
            groups:
                - outrider
            fpkmCutoff: 1
            implementation: autoencoder
            padjCutoff: 0.05
            zScoreCutoff: 0
            genesToTest: null
            maxTestedDimensionProportion: 3
            yieldSize: 2000000
        
        aberrantSplicing:
            run: false
            groups:
                - fraser
            recount: false
            longRead: false
            keepNonStandardChrs: true
            filter: true
            minExpressionInOneSample: 20
            quantileMinExpression: 10
            minDeltaPsi: 0.05
            implementation: PCA
            padjCutoff: 0.1
            maxTestedDimensionProportion: 6
            genesToTest: null
            FRASER_version: "FRASER2"
            deltaPsiCutoff : 0.1
            quantileForFiltering: 0.75

        mae:
            run: false
            groups:
                - mae
            gatkIgnoreHeaderCheck: true
            padjCutoff: .05
            allelicRatioCutoff: 0.8
            addAF: true
            maxAF: .001
            maxVarFreqCohort: .05
            # VCF-BAM matching
            qcVcf: null
            qcGroups: mae
            dnaRnaMatchCutoff: 0.85
            
        rnaVariantCalling:
            run: false
            groups:
                - batch_0
            highQualityVCFs:
                - Data/Mills_and_1000G_gold_standard.indels.hg19.sites.chrPrefix.vcf.gz
                - Data/1000G_phase1.snps.high_confidence.hg19.sites.chrPrefix.vcf.gz
            dbSNP: Data/00-All.vcf.gz
            repeat_mask: Data/hg19_repeatMasker_sorted.chrPrefix.bed
            createSingleVCF: true
            addAF: true
            maxAF: 0.001
            maxVarFreqCohort: 0.05
            hcArgs: ""
            minAlt: 3
            yieldSize: 100000

        tools:
            gatkCmd: gatk
            bcftoolsCmd: bcftools
            samtoolsCmd: samtools
        """
    )


def update_config(
    yaml_object: dict,
    genome: str,
    gtf: str,
    genome_assembly: str,
    padjcutoff: float,
    zscorecutoff: float,
    drop_module: str,
) -> dict:
    """
    Updates config file to add correct genome, gtf,
    adjusted p-value, and Z-score for the module to be run.
    """
    gtf_name = Path(gtf).name
    gtf_without_ext = Path(gtf).stem
    genome_name = Path(genome).name

    yaml_object["genome"] = genome_name
    yaml_object["root"] = "output"
    yaml_object["htmlOutputPath"] = "output/html"
    yaml_object["sampleAnnotation"] = "sample_annotation.tsv"
    yaml_object["geneAnnotation"][gtf_without_ext] = gtf_name
    yaml_object["geneAnnotation"].pop("gtf", None)
    yaml_object["exportCounts"]["geneAnnotations"] = [gtf_without_ext]
    yaml_object["genomeAssembly"] = genome_assembly

    ## Export counts
    if drop_module == "AE":
        yaml_object["aberrantExpression"]["run"] = ["true"]
        yaml_object["aberrantExpression"]["padjCutoff"] = padjcutoff
        yaml_object["aberrantExpression"]["zScoreCutoff"] = zscorecutoff
    elif drop_module == "AS":
        yaml_object["aberrantSplicing"]["run"] = ["true"]
        yaml_object["aberrantSplicing"]["padjCutoff"] = padjcutoff
    elif drop_module == "MAE":
        yaml_object["mae"]["run"] = ["true"]
    return yaml_object


def write_yaml(out_path: str, yaml_object: dict):
    with open(out_path, "w") as outfile:
        yaml.dump(yaml_object, outfile, sort_keys=False)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        formatter_class=argparse.MetavarTypeHelpFormatter,
        description="""Generate config file for DROP.""",
    )

    parser.add_argument("--genome_fasta", type=str, help="Specify genome fasta base name", required=True)

    parser.add_argument("--gtf", type=str, help="Specify gtf file name", required=True)

    parser.add_argument("--output", type=str, help="Specify output file", required=True)

    parser.add_argument(
        "--genome_assembly",
        type=str,
        help="Specify genome for drop can be either hg19/hs37d5 or hg38/GRCh38",
        required=True,
    )

    parser.add_argument("--padjcutoff", type=float, help="Specify adjusted p-value cut-off", required=True)

    parser.add_argument(
        "--zscorecutoff",
        default=0,
        type=float,
        help="Specify z-score cut-off, this is an optional value",
        required=False,
    )

    parser.add_argument("--drop_module", type=str, help="Specify module to run: AE, AS or MAE", required=True)

    args = parser.parse_args()
    master_config = update_config(
        yaml_object=get_config(),
        genome=args.genome_fasta,
        gtf=args.gtf,
        genome_assembly=args.genome_assembly,
        padjcutoff=args.padjcutoff,
        zscorecutoff=args.zscorecutoff,
        drop_module=args.drop_module,
    )
    write_yaml(out_path=args.output, yaml_object=master_config)
