#!/usr/bin/env python3

import argparse
import csv
from pandas import read_csv, DataFrame, concat, isna
import os
import re

SCRIPT_VERSION = "1.3"
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
    "SEX",
]


def write_sample_annotation_to_tsv(
    bam: str,
    dna_vcf: str,
    samples: str,
    strandedness: str,
    single_end: str,
    sex: str,
    drop_group_sample: str,
    out_file: str,
    gtf: str,
):
    """Write the Sample Annotation tsv file."""
    with open(out_file, "w") as tsv_file:
        writer = csv.DictWriter(tsv_file, fieldnames=SAMPLE_ANNOTATION_COLUMNS, delimiter="\t")
        writer.writeheader()
        for index, id in enumerate(samples):
            sa_dict: dict = {}.fromkeys(SAMPLE_ANNOTATION_COLUMNS, "NA")
            sa_dict["RNA_ID"] = id
            sa_dict["DNA_ID"] = id
            sa_dict["STRAND"] = is_stranded(strandedness[index])
            sa_dict["SEX"] = sex[index]
            sa_dict["PAIRED_END"] = is_paired_end(single_end[index])
            sa_dict["RNA_BAM_FILE"] = bam[index]
            value = dna_vcf[index].strip()
            sa_dict["DNA_VCF_FILE"] = os.path.basename(value) if value not in ("", "NA") else "NA"
            sa_dict["DROP_GROUP"] = drop_group_sample + ",mae" if sa_dict["DNA_VCF_FILE"] != "NA" else drop_group_sample
            writer.writerow(sa_dict)


def is_paired_end(single_end: str) -> bool:
    """Logical funciton to determine if a sample is paired end"""
    if single_end.lower() == "false":
        return True
    return False


def is_stranded(strandedness: str) -> str:
    """Logical funciton to determine sample strandness"""
    if strandedness.lower() == "reverse":
        return "reverse"
    elif strandedness.lower() == "forward":
        return "yes"
    else:
        return "no"


def count_mode(sample_count_mode: str) -> str:
    """Logical function to determine if count mode is given or default "IntersectionStrict" should be used"""
    if isna(sample_count_mode) or sample_count_mode == "" or sample_count_mode == "NA":
        return "IntersectionStrict"
    else:
        return sample_count_mode


def count_overlaps(sample_count_overlap: str) -> str:
    """Logical funciton to determine if count overlap is given or default "TRUE" should be used"""
    if isna(sample_count_overlap) or sample_count_overlap == "" or sample_count_overlap == "NA":
        return True
    else:
        return sample_count_overlap


def write_final_annot_to_tsv(ref_count_file: str, ref_annot: str, out_file: str):
    """
    Concatenates the Sample Annotation produced by SampleAnnotation with the one
    provided for the reference samples, if one is provided, checking for duplicate sample IDs
    """
    df_samples: DataFrame = read_csv("drop_annotation_given_samples.tsv", sep="\t")

    # Remove sex column if no non-NA values are provided
    if df_samples["SEX"].count() == 0:
        df_samples.drop(columns=["SEX"], inplace=True)

    if ref_annot == "None" or ref_count_file == "None":
        print(
            "No reference samples were provided by the user see usage of --ref_count_file and --ref_annot if you want to provide reference samples"
        )
        if df_samples.shape[0] < 50:
            print("At least 30 samples are required for Aberrant Splicing and 50 for Aberrant expression")
            print(f"Only {df_samples.shape[0]} samples were provided by the user")
        df_samples.fillna("NA", inplace=True)
        df_samples["COUNT_MODE"] = "IntersectionStrict"
        df_samples["COUNT_OVERLAPS"] = True
        df_samples.to_csv(out_file, index=False, sep="\t")
    else:
        df_reference: DataFrame = read_csv(ref_annot, sep="\t")
        df_reference["GENE_COUNTS_FILE"] = ref_count_file
        df_reference["SPLICE_COUNTS_DIR"] = (
            df_reference["SPLICE_COUNTS_DIR"]
            .fillna("NA")
            .str.rstrip("/")
            .apply(lambda x: os.path.basename(x) if x != "NA" else x)
        )
        df_reference["DROP_GROUP"] = df_reference["DROP_GROUP"].str.replace(" ", "")
        df_samples["COUNT_OVERLAPS"] = count_overlaps(df_reference["COUNT_OVERLAPS"].iloc[0])
        df_samples["COUNT_MODE"] = count_mode(df_reference["COUNT_MODE"].iloc[0])
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
    parser.add_argument(
        "--bam",
        type=str,
        nargs="+",
        help="bam files for the analyzed samples",
        required=True,
    )
    parser.add_argument(
        "--dna_vcf",
        type=str,
        nargs="+",
        help="DNA VCF files to perform MAE on the analyzed samples",
        required=True,
    )
    parser.add_argument(
        "--samples",
        type=str,
        nargs="+",
        help="corresponding sample name",
        required=True,
    )
    parser.add_argument("--strandedness", type=str, nargs="+", help="strandedness of RNA", required=True)
    parser.add_argument("--sex", type=str, nargs="+", help="Sex of samples", required=True)
    parser.add_argument(
        "--single_end",
        type=str,
        nargs="+",
        help="is the sample paired end?",
        required=True,
    )
    parser.add_argument(
        "--ref_count_file",
        type=str,
        default="None",
        help="A tsv file of gene counts for reference samples.",
        required=False,
    )
    parser.add_argument(
        "--ref_annot",
        type=str,
        default="None",
        help="Path to reference annotation tsv",
        required=False,
    )
    parser.add_argument(
        "--drop_group_sample",
        type=str,
        default="None",
        help="Drop group of analyzed samples",
        required=False,
    )
    parser.add_argument(
        "--gtf",
        type=str,
        help="Specify gtf file name used to run",
        required=True,
    )
    parser.add_argument("--output", type=str, help="Path to save to", required=True)
    parser.add_argument("--version", action="version", version=SCRIPT_VERSION)
    return parser.parse_args(argv)


def main():
    """Coordinate argument parsing and program execution."""
    args = parse_args()
    write_sample_annotation_to_tsv(
        bam=args.bam,
        dna_vcf=args.dna_vcf,
        samples=args.samples,
        strandedness=args.strandedness,
        single_end=args.single_end,
        sex=args.sex,
        gtf=args.gtf,
        drop_group_sample=args.drop_group_sample,
        out_file="drop_annotation_given_samples.tsv",
    )

    write_final_annot_to_tsv(
        ref_count_file=args.ref_count_file,
        ref_annot=args.ref_annot,
        out_file=args.output,
    )


if __name__ == "__main__":
    main()
