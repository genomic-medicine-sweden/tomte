#!/usr/bin/env python3

import argparse
import csv
from pandas import read_csv, DataFrame
from pathlib import Path

SCRIPT_VERSION = "v1.0"
GENE_COUNTS_PATH = Path("exported_counts/geneCounts.tsv.gz")
SPLICE_COUNTS_DIR = Path("exported_counts")
OUTPUT_FILE = Path("exported_counts/sampleAnnotation.tsv")


def modify_gene_counts_df(df: DataFrame, col_name: str, run: str, value_in: str):
    """Modifies column col_name in df if run is true to make all
    rows equal value_in. If run is false it will make all rows NA"""
    if run == "true":
        df[col_name] = value_in
    else:
        df[col_name] = "NA"
    return df


def modify_and_write_sample_annotation(
    sample_annot: str, ae_run: str, as_run: str, gtf: str
):
    """
    Modifies and writes Sample Annotation produced by DROP to make one
    that can be used as input for Tomte
    """
    df_samples: DataFrame = read_csv(sample_annot, sep="\t")
    df_samples["RNA_BAM_FILE"] = "NA"
    df_samples["GENE_ANNOTATION"] = df_samples["GENE_ANNOTATION"].fillna(gtf)
    df_samples = modify_gene_counts_df(
        df=df_samples,
        col_name="GENE_COUNTS_FILE",
        run=ae_run,
        value_in=str(GENE_COUNTS_PATH),
    )
    df_samples = modify_gene_counts_df(
        df=df_samples,
        col_name="SPLICE_COUNTS_DIR",
        run=as_run,
        value_in=str(SPLICE_COUNTS_DIR),
    )
    df_samples.to_csv(OUTPUT_FILE, index=False, sep="\t")


def parse_args(argv=None):
    """Define and immediately parse command line arguments."""
    parser = argparse.ArgumentParser(
        formatter_class=argparse.MetavarTypeHelpFormatter,
        description="""Generate DROP sample annotation for exported db.""",
    )
    parser.add_argument(
        "--sample_annot",
        type=str,
        help="original sample annotation in export_counts folder",
        required=True,
    )
    parser.add_argument(
        "--ae_run",
        type=str,
        help="Was aberrant expression run?",
        required=True,
    )
    parser.add_argument(
        "--as_run",
        type=str,
        help="Was aberrant splicing run?",
        required=True,
    )
    parser.add_argument(
        "--gtf",
        type=str,
        help="Specify gtf file name used to run",
        required=True,
    )
    parser.add_argument("--version", action="version", version=SCRIPT_VERSION)
    return parser.parse_args(argv)


def main():
    """Coordinate argument parsing and program execution."""
    args = parse_args()
    modify_and_write_sample_annotation(
        sample_annot=args.sample_annot,
        ae_run=args.ae_run,
        as_run=args.as_run,
        gtf=args.gtf,
    )


if __name__ == "__main__":
    main()
