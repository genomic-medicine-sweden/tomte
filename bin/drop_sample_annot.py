#!/usr/bin/env python3

import argparse
from pathlib import Path
import csv
from pandas import read_csv, DataFrame, concat
import sys

SCRIPT_VERSION = "v1.0"
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


def write_sample_annotation_to_tsv(
    bam: str, samples: str, strandedness: str, single_end: str, gtf: str, count_file: str, out_file: str
):
    """Write the Sample Annotation tsv file."""
    with open(out_file, "w") as tsv_file:
        fieldnames = SAMPLE_ANNOTATION_COLUMNS
        writer = csv.DictWriter(tsv_file, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for index, id in enumerate(samples):
            sa_dict: dict = {}.fromkeys(fieldnames, "NA")
            sa_dict["RNA_ID"] = id
            sa_dict["DROP_GROUP"] = "outrider,fraser"
            sa_dict["GENE_COUNTS_FILE"] = count_file
            sa_dict["GENE_ANNOTATION"] = Path(gtf).stem
            sa_dict["STRAND"] = strandedness[index]
            sa_dict["PAIRED_END"] = is_paired_end(single_end[index])
            sa_dict["RNA_BAM_FILE"] = bam[index]
            writer.writerow(sa_dict)


def is_paired_end(single_end: str) -> bool:
    """Logical funciton to determine if a sample is paired end"""
    if single_end.lower() == "false":
        return True
    return False


def write_final_annot_to_tsv(count_file: str, ref_annot: str, out_file: str):
    """
    Concatenates the Sample Annotation produced by SampleAnnotation with the one
    provided for the reference samples, checking for duplicate sample IDs
    """
    df_samples: DataFrame = read_csv("drop_annotation_given_samples.tsv", sep="\t")
    df_reference: DataFrame = read_csv(ref_annot, sep="\t")
    df_reference["GENE_COUNTS_FILE"] = count_file
    df_samples["COUNT_OVERLAPS"] = df_reference["COUNT_OVERLAPS"].iloc[0]
    df_samples["COUNT_MODE"] = df_reference["COUNT_MODE"].iloc[0]
    df_samples["HPO_TERMS"] = df_reference["HPO_TERMS"].iloc[0]
    for id in df_samples["RNA_ID"]:
        df_reference = df_reference[df_reference["RNA_ID"].str.contains(id) == False]
    df: DataFrame = concat([df_samples, df_reference]).reset_index(drop=True)
    df.fillna("NA", inplace=True)
    df.to_csv(out_file, index=False, sep="\t")


def parse_args(argv=None):
    """Define and immediately parse command line arguments."""
    parser = argparse.ArgumentParser(
        formatter_class=argparse.MetavarTypeHelpFormatter,
        description="""Generate DROP sample annotation for patients.""",
    )
    parser.add_argument("--bam", type=str, nargs="+", help="bam files for the patient", required=True)
    parser.add_argument("--samples", type=str, nargs="+", help="corresponding sample name", required=True)
    parser.add_argument("--strandedness", type=str, nargs="+", help="strandedness of RNA", required=True)
    parser.add_argument("--single_end", type=str, nargs="+", help="is the sample paired end?", required=True)
    parser.add_argument("--gtf", type=str, help="Transcript annotation file in gtf format", required=True)
    parser.add_argument(
        "--count_file", type=str, help="A tsv file of gene counts for all processed samples.", required=True
    )
    parser.add_argument("--ref_annot", type=str, help="Path to reference annotation tsv", required=True)
    parser.add_argument("--output", type=str, help="Path to save to", required=True)
    parser.add_argument("--version", action="version", version=SCRIPT_VERSION)
    return parser.parse_args(argv)


def main(argv=None):
    """Coordinate argument parsing and program execution."""
    args = parse_args(argv)
    write_sample_annotation_to_tsv(
        bam=args.bam,
        samples=args.samples,
        strandedness=args.strandedness,
        single_end=args.single_end,
        gtf=args.gtf,
        count_file=args.count_file,
        out_file="drop_annotation_given_samples.tsv",
    )
    write_final_annot_to_tsv(count_file=args.count_file, ref_annot=args.ref_annot, out_file=args.output)


if __name__ == "__main__":
    main()
