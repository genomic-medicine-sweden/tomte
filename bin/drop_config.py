#!/usr/bin/env python3

import argparse
from pathlib import Path
import yaml
from typing import Dict, Any
from copy import deepcopy

SCRIPT_VERSION = "v2.0"
CONFIG_YAML = {
    "projectTitle": "DROP: Detection of RNA Outliers Pipeline",
    "root": None,
    "htmlOutputPath": None,
    "indexWithFolderName": True,
    "hpoFile": None,
    "sampleAnnotation": None,
    "geneAnnotation": {"gtf": None},
    "genomeAssembly": "genome",
    "exportCounts": {"geneAnnotations": ["gtf"], "excludeGroups": [None]},
    "genome": None,
    "aberrantExpression": {
        "run": False,
        "groups": ["outrider"],
        "fpkmCutoff": 1,
        "implementation": "autoencoder",
        "padjCutoff": 0.05,
        "zScoreCutoff": 0,
        "genesToTest": None,
        "maxTestedDimensionProportion": 3,
        "yieldSize": 2000000,
    },
    "aberrantSplicing": {
        "run": False,
        "groups": ["fraser"],
        "recount": False,
        "longRead": False,
        "keepNonStandardChrs": True,
        "filter": True,
        "minExpressionInOneSample": 20,
        "quantileMinExpression": 10,
        "minDeltaPsi": 0.05,
        "implementation": "PCA",
        "padjCutoff": 0.1,
        "maxTestedDimensionProportion": 6,
        "genesToTest": None,
        "FRASER_version": "FRASER2",
        "deltaPsiCutoff": 0.1,
        "quantileForFiltering": 0.75,
    },
    "mae": {
        "run": False,
        "groups": ["mae"],
        "gatkIgnoreHeaderCheck": True,
        "padjCutoff": 0.05,
        "allelicRatioCutoff": 0.8,
        "addAF": True,
        "maxAF": 0.001,
        "maxVarFreqCohort": 0.05,
        "qcVcf": None,
        "qcGroups": "mae",
        "dnaRnaMatchCutoff": 0.85,
    },
    "rnaVariantCalling": {
        "run": False,
        "groups": ["batch_0"],
        "highQualityVCFs": [
            "Data/Mills_and_1000G_gold_standard.indels.hg19.sites.chrPrefix.vcf.gz",
            "Data/1000G_phase1.snps.high_confidence.hg19.sites.chrPrefix.vcf.gz",
        ],
        "dbSNP": "Data/00-All.vcf.gz",
        "repeat_mask": "Data/hg19_repeatMasker_sorted.chrPrefix.bed",
        "createSingleVCF": True,
        "addAF": True,
        "maxAF": 0.001,
        "maxVarFreqCohort": 0.05,
        "hcArgs": "",
        "minAlt": 3,
        "yieldSize": 100000,
    },
    "tools": {"gatkCmd": "gatk", "bcftoolsCmd": "bcftools", "samtoolsCmd": "samtools"},
}


def update_config(
    genome: Path,
    gtf: Path,
    genome_assembly: str,
    drop_group_samples: str,
    drop_other_group_samples: str,
    padjcutoff: float,
    zscorecutoff: float,
    drop_module: str,
) -> Dict:
    """
    Updates config file to add correct genome, gtf,
    adjusted p-value, and Z-score for the module to be run.
    """
    gtf_without_ext = gtf.stem
    genome_name = genome.name
    config_copy: Dict[str, Any] = deepcopy(CONFIG_YAML)

    config_copy["genome"] = genome_name
    config_copy["root"] = "output"
    config_copy["htmlOutputPath"] = "output/html"
    config_copy["sampleAnnotation"] = "sample_annotation.tsv"
    config_copy["geneAnnotation"][gtf_without_ext] = str(gtf)
    config_copy["geneAnnotation"].pop("gtf", None)
    config_copy["exportCounts"]["geneAnnotations"] = [gtf_without_ext]
    config_copy["genomeAssembly"] = genome_assembly
    config_copy["exportCounts"]["excludeGroups"] = [drop_other_group_samples]

    # Export counts
    if drop_module == "AE":
        config_copy["aberrantExpression"]["run"] = ["true"]
        config_copy["aberrantExpression"]["groups"] = [drop_group_samples]
        config_copy["aberrantExpression"]["padjCutoff"] = padjcutoff
        config_copy["aberrantExpression"]["zScoreCutoff"] = zscorecutoff
    elif drop_module == "AS":
        config_copy["aberrantSplicing"]["run"] = ["true"]
        config_copy["aberrantSplicing"]["groups"] = [drop_group_samples]
        config_copy["aberrantSplicing"]["padjCutoff"] = padjcutoff
    elif drop_module == "MAE":
        config_copy["mae"]["run"] = ["true"]
    return config_copy


def write_yaml(out_path: str, yaml_object: dict):
    with open(out_path, "w") as outfile:
        yaml.dump(yaml_object, outfile, sort_keys=False)


def parse_args(argv=None):
    """Define and immediately parse command line arguments."""
    parser = argparse.ArgumentParser(
        formatter_class=argparse.MetavarTypeHelpFormatter,
        description="""Generate config file for DROP.""",
    )

    parser.add_argument(
        "--genome_fasta",
        type=Path,
        help="Specify genome fasta base name",
        required=True,
    )
    parser.add_argument(
        "--gtf",
        type=Path,
        help="Specify gtf file name",
        required=True,
    )
    parser.add_argument(
        "--output",
        type=str,
        help="Specify output file",
        required=True,
    )
    parser.add_argument(
        "--genome_assembly",
        type=str,
        help="Specify genome for drop can be either hg19/hs37d5 or hg38/GRCh38.",
        choices=["hg19", "hs37d5", "hg38", "GRCh38"],
        required=True,
    )
    parser.add_argument(
        "--drop_group_samples",
        type=str,
        help="Specify drop group to analyse",
        required=True,
    )
    parser.add_argument(
        "--drop_other_group_samples",
        type=str,
        help="Specify the drop group to exclude in exportCounts",
        required=True,
    )
    parser.add_argument(
        "--padjcutoff",
        type=float,
        help="Specify adjusted p-value cut-off",
        required=True,
    )
    parser.add_argument(
        "--zscorecutoff",
        default=0,
        type=float,
        help="Specify z-score cut-off, this is an optional value",
        required=False,
    )
    parser.add_argument(
        "--drop_module",
        type=str,
        help="Specify module to run: AE, AS or MAE",
        required=True,
    )
    parser.add_argument(
        "--version",
        action="version",
        version=SCRIPT_VERSION,
    )
    return parser.parse_args(argv)


def main():
    """Coordinate argument parsing and program execution."""
    args = parse_args()
    master_config: Dict[str, Any] = update_config(
        genome=args.genome_fasta,
        gtf=args.gtf,
        genome_assembly=args.genome_assembly,
        drop_group_samples=args.drop_group_samples,
        drop_other_group_samples=args.drop_other_group_samples,
        padjcutoff=args.padjcutoff,
        zscorecutoff=args.zscorecutoff,
        drop_module=args.drop_module,
    )
    write_yaml(out_path=args.output, yaml_object=master_config)


if __name__ == "__main__":
    main()
