#!/usr/bin/env python3

import argparse
import re
from collections import OrderedDict
from pathlib import Path
import csv
import os
import pandas as pd


class SampleAnnotation:
    """SampleAnnotation class"""

    SAMPLE_ANNOTATION_COLUMNS = [
        "RNA_ID",
        "RNA_BAM_FILE",
        "DNA_VCF_FILE",
        "DNA_ID",
        "DROP_GROUP",
        "PAIRED_END",
        "COUNT_MODE",
        "COUNT_OVERLAPS",
        "SPLICE_COUNTS_DIR",
        "STRAND",
        "HPO_TERMS",
        "GENE_COUNTS_FILE",
        "GENE_ANNOTATION",
        "GENOME",
    ]

    def __init__(self, bam, sample, strandedness, single_end, gtf, count_file, out_file):
        """Write the Sample Annotation tsv file"""
        with open(out_file, "w") as tsv_file:
            fieldnames = self.SAMPLE_ANNOTATION_COLUMNS
            writer = csv.DictWriter(tsv_file, fieldnames=fieldnames, delimiter="\t")

            writer.writeheader()

            for index, id in enumerate(sample):
                sa_dict = {}.fromkeys(fieldnames, "NA")
                sa_dict["RNA_ID"] = id.strip("[],")
                sa_dict["DROP_GROUP"] = "outrider,fraser"
                sa_dict["GENE_COUNTS_FILE"] = count_file
                sa_dict["GENE_ANNOTATION"] = Path(gtf).stem
                sa_dict["STRAND"] = strandedness[index].strip("[],")
                paired_end_func = lambda x: True if x.strip("[],").lower() == "false" else False
                sa_dict["PAIRED_END"] = paired_end_func(single_end[index])
                sa_dict["RNA_BAM_FILE"] = bam[index]
                writer.writerow(sa_dict)


def final_annot(count_file: Path, ref_annot: Path, out_file: Path):
    """Concatinates the Sample Annotation produced by SampleAnnotation with the one
    provided for the reference samples, checking for duplicate sample IDs"""
    df_samples = pd.read_csv("drop_annotation_given_samples.tsv", sep="\t")
    df_reference = pd.read_csv(ref_annot, sep="\t")
    df_reference["GENE_COUNTS_FILE"] = count_file
    df_samples["COUNT_OVERLAPS"] = df_reference["COUNT_OVERLAPS"].iloc[0]
    df_samples["COUNT_MODE"] = df_reference["COUNT_MODE"].iloc[0]
    df_samples["HPO_TERMS"] = df_reference["HPO_TERMS"].iloc[0]
    for id in df_samples["RNA_ID"]:
        df_reference = df_reference[df_reference["RNA_ID"].str.contains(id) == False]
    df = pd.concat([df_samples, df_reference]).reset_index(drop=True)
    df.fillna("NA", inplace=True)
    df.to_csv(out_file, index=False, sep="\t")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        formatter_class=argparse.MetavarTypeHelpFormatter,
        description="""Generate DROP sample annotation for patients.""",
    )

    parser.add_argument(
        "--bam",
        type=str,
        nargs="+",
        help="bam files for the patient",
        required=True,
    )

    parser.add_argument(
        "--sample",
        type=str,
        nargs="+",
        help="corresponding sample name",
        required=True,
    )

    parser.add_argument(
        "--strandedness",
        type=str,
        nargs="+",
        help="strandedness of RNA",
        required=True,
    )

    parser.add_argument(
        "--single_end",
        type=str,
        nargs="+",
        help="is the sample paired end?",
        required=True,
    )

    parser.add_argument(
        "--gtf",
        type=str,
        help="Transcript annotation file in gtf format",
        required=True,
    )

    parser.add_argument(
        "--count_file", type=str, help="A tsv file of gene counts for all processed samples.", required=True
    )

    parser.add_argument(
        "--ref_annot",
        type=str,
        help="Path to reference annotation tsv",
        required=True,
    )

    parser.add_argument(
        "--output",
        type=str,
        help="Path to save to",
        required=True,
    )

    args = parser.parse_args()
    SampleAnnotation(
        args.bam, args.sample, args.strandedness, args.single_end, args.gtf, args.count_file, "drop_annotation_given_samples.tsv"
    )
    final_annot(args.count_file, args.ref_annot, args.output)
