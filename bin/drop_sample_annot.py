#!/usr/bin/env python3

import argparse
import csv
import os
from pathlib import Path


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
        "STRAND",
        "HPO_TERMS",
        "GENE_COUNTS_FILE",
        "GENE_ANNOTATION",
        "GENOME",
    ]

    def __init__(self, cnts_file, out_file, gtf_file):
        """Create SampleAnnotation given the parameters"""
        self.cnts_file = cnts_file
        self.out_file = out_file
        self.gtf_file = gtf_file

    def parse_header(self):
        """Parse the first line of gene counts file"""
        samples = []
        with open(self.cnts_file) as file_object:
            samples = file_object.readline().split()

        del samples[0]  # remove GeneID field

        sample_cnt_file = {sample: os.path.basename(self.cnts_file) for sample in samples}

        return samples, sample_cnt_file

    def write_table(self):
        """Write the Sample Annotation tsv file"""
        with open(self.out_file, "w") as tsv_file:
            fieldnames = self.SAMPLE_ANNOTATION_COLUMNS
            sample_ids, sample_cnt_file = self.parse_header()
            writer = csv.DictWriter(tsv_file, fieldnames=fieldnames, delimiter="\t")

            writer.writeheader()
            for id in sample_ids:
                sa_dict = {}.fromkeys(fieldnames, "NA")
                sa_dict["RNA_ID"] = id
                sa_dict["DROP_GROUP"] = "outrider"
                sa_dict["GENE_COUNTS_FILE"] = sample_cnt_file[id]
                sa_dict["GENE_ANNOTATION"] = Path(self.gtf_file).name
                writer.writerow(sa_dict)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        formatter_class=argparse.MetavarTypeHelpFormatter,
        description="""Generate sample annotation file for DROP.""",
    )

    parser.add_argument(
        "--count_file",
        type=str,
        help="A tsv file of gene counts for all processed samples.",
        required=True,
    )

    parser.add_argument(
        "--output",
        type=str,
        help="Path to save to",
        required=True,
    )

    parser.add_argument(
        "--gtf",
        type=str,
        help="Specify gtf file name",
        required=True,
    )

    args = parser.parse_args()

    SampleAnnotation(args.count_file, args.output, args.gtf).write_table()
