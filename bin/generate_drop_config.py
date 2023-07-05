#!/usr/bin/env python3

import argparse
from pathlib import Path

import yaml


def read_config():
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

        genomeAssembly: hg38  # either hg19/hs37d5 or hg38/GRCh38
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

        aberrantExpression:
            run: true
            groups:
                - outrider
            fpkmCutoff: 1
            implementation: autoencoder
            padjCutoff: 0.05
            zScoreCutoff: 0
            maxTestedDimensionProportion: 3
            dassie:
                tssWindow: 500
                pasWindow: 1000

        aberrantSplicing:
            run: false
            groups:
                - outrider
            recount: false
            longRead: false
            keepNonStandardChrs: true
            filter: true
            minExpressionInOneSample: 20
            minDeltaPsi: 0.05
            implementation: PCA-BB-Decoder
            padjCutoff: 0.05
            zScoreCutoff: 0
            deltaPsiCutoff: 0.3
            maxTestedDimensionProportion: 6

        mae:
            run: false
            groups:
                - outrider
            gatkIgnoreHeaderCheck: true
            padjCutoff: .05
            allelicRatioCutoff: 0.8
            addAF: false
            maxAF: .001
            maxVarFreqCohort: .05
            # VCF-BAM matching
            qcVcf: null
            qcGroups: mae

        tools:
            gatkCmd: gatk
            bcftoolsCmd: bcftools
            samtoolsCmd: samtools
        """
    )


def update_config(yaml_object, genome, gtf):

    gtf_name = Path(gtf).name
    genome_name = Path(genome).name

    yaml_object["genome"] = genome_name
    yaml_object["root"] = "output"
    yaml_object["htmlOutputPath"] = "output/html"
    yaml_object["sampleAnnotation"] = "sample_annotation.tsv"
    yaml_object["geneAnnotation"][gtf_name] = gtf_name
    yaml_object["geneAnnotation"].pop("gtf", None)

    ## Export counts
    yaml_object["exportCounts"]["geneAnnotations"] = [gtf_name]

    ## Expression module
    yaml_object["aberrantExpression"]["groups"] = ["outrider"]

    ## Splicing module

    ## MAE module
    return yaml_object


def write_yaml(out_path, yaml_object):
    with open(out_path, "w") as outfile:
        yaml.dump(yaml_object, outfile, sort_keys=False)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        formatter_class=argparse.MetavarTypeHelpFormatter,
        description="""Generate config file for DROP.""",
    )

    parser.add_argument(
        "--genome_fasta",
        type=str,
        help="Specify genome fasta base name",
    )

    parser.add_argument(
        "--gtf",
        type=str,
        help="Specify gtf file name",
    )

    parser.add_argument(
        "--output",
        type=str,
        help="Specify output file",
    )

    args = parser.parse_args()

    yaml_object = read_config()
    master_config = update_config(yaml_object=yaml_object, genome=args.genome_fasta, gtf=args.gtf)
    write_yaml(out_path=args.output, yaml_object=master_config)
