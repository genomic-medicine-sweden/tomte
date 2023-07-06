#!/usr/bin/env python3

import argparse
import re
from collections import OrderedDict
from pathlib import Path

import pandas as pd

translator = {
    "unstranded": 1,
    "forward": 2,
    "reverse": 3,
}


def get_non_std_genes(gtf: str) -> set[str]:
    """Create list of genes not belonging to chr1-21 or chrM"""

    std_chrs = [
        "chr1",
        "chr2",
        "chr3",
        "chr4",
        "chr5",
        "chr6",
        "chr7",
        "chr8",
        "chr9",
        "chr10",
        "chr11",
        "chr12",
        "chr13",
        "chr14",
        "chr15",
        "chr16",
        "chr17",
        "chr18",
        "chr19",
        "chr20",
        "chr21",
        "chr22",
        "chrX",
        "chrY",
        "chrM",
    ]

    gene_id_regex = re.compile('gene_id "(.+?)"')
    genes_to_exclude = []

    with open(gtf, "r") as gtf_file:
        for line in gtf_file:
            if line.startswith("#"):
                continue

            if line.split()[0] in std_chrs:
                continue

            gene_id = re.search(gene_id_regex, line)
            genes_to_exclude.append(gene_id.group(1))

    return set(genes_to_exclude)


def read_star_gene_cnts(sample: str, star: Path, strandedness: str) -> dict:
    """Read gene count file(s) from STAR output."""
    sample_ids = {}
    gene_ids = {}
    with open(star) as in_tab:
        for line in in_tab:
            if not line.startswith("N_"):
                gene_id = line.split()[0]
                strand = translator[strandedness.lower()]
                counts = line.split()[strand]
                gene_ids[gene_id] = int(counts)
    gene_ids = OrderedDict(sorted(gene_ids.items()))
    sample_ids[sample] = gene_ids
    return sample_ids


def transform_to_table(gene_ids_dict, outfile, genes_to_exlude, ref_count_file: Path):
    """Transform in dictionary into tsv friendly."""

    one_sample = next(iter(gene_ids_dict))
    gene_list = list(gene_ids_dict[one_sample].keys())
    genes = {}

    for gene in gene_list:
        genes[gene] = []
        for sample in gene_ids_dict:
            genes[gene].append(gene_ids_dict[sample][gene])

    count_table = pd.DataFrame.from_dict(genes, orient="index", columns=gene_ids_dict.keys())
    count_table.index.name = "geneID"

    final_table = None

    if ref_count_file:
        ref_table = pd.read_csv(
            ref_count_file,
            sep="\t",
            header=0,
            index_col=0,
        )
        final_table = count_table.combine_first(ref_table)
    else:
        final_table = count_table

    final_table.drop(genes_to_exclude, inplace=True)
    final_table.to_csv(outfile, sep="\t", header=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        formatter_class=argparse.MetavarTypeHelpFormatter,
        description="""Generate collated gene counts from each STAR output.""",
    )
    parser.add_argument("--star", type=str, nargs="+", help="*ReadsPerGene.out.tab from STAR", required=True)
    parser.add_argument("--sample", type=str, nargs="+", help="corresponding sample name", required=True)
    parser.add_argument("--strandedness", type=str, help="strandedness of RNA")
    parser.add_argument("--output", type=str, help="output tsv file name", required=True)
    parser.add_argument("--gtf", type=str, help="Transcript annotation file in gtf format", required=True)
    parser.add_argument("--ref_count_file", type=str, help="Optional reference count set")

    args = parser.parse_args()
    master_dict = {}
    for index, sample_id in enumerate(args.sample):
        sample_id = re.sub(r"[\[\],]", "", sample_id)
        master_dict.update(read_star_gene_cnts(sample=sample_id, star=args.star[index], strandedness=args.strandedness))

    genes_to_exclude = get_non_std_genes(args.gtf)
    transform_to_table(master_dict, args.output, genes_to_exclude, args.ref_count_file)
