#!/usr/bin/env python3

import argparse
from pathlib import Path
from typing import Optional, Set, Tuple, Dict
import csv
from collections import defaultdict
import json

VERSION = "1.0.0"


def main(
    hbs_path: Path,
    gene_counts_path: Path,
    strandedness: Optional[str],
    out_json: Optional[Path],
    verbose: bool,
):

    target_genes = get_target_genes(hbs_path)

    if strandedness not in ["forward", "reverse"] and strandedness is not None:
        raise ValueError(
            f"strandedness option should be forward, reverse or not assinged. Found: {strandedness}"
        )

    (nbr_genes, total_count, target_genes_counts) = summarize_gene_counts(
        gene_counts_path, target_genes, strandedness, verbose
    )

    metrics = dict()
    metrics["nbr_genes"] = nbr_genes
    metrics["total_count"] = total_count
    metrics["average_count"] = total_count / nbr_genes
    metrics["target_genes_count"] = sum(list(target_genes_counts.values()))
    metrics["target_genes_frac"] = metrics["target_genes_count"] / metrics["total_count"]
    metrics["per_target_gene_count"] = target_genes_counts

    if verbose:
        print(json.dumps(metrics))

    if out_json is not None:
        with out_json.open("w") as out_fh:
            out_fh.write(json.dumps(metrics))


def summarize_gene_counts(
    gene_counts_path: Path, target_genes: Set[str], strandedness: Optional[str], verbose: bool
) -> Tuple[int, int, Dict[str, int]]:
    nbr_genes = 0
    total_count = 0
    target_genes_counts = defaultdict(int)
    with gene_counts_path.open("r") as in_fh:
        for line in in_fh:
            (gene_id_raw, both, fw, rv) = line.rstrip().split("\t")

            # Remove ensembl gene version part
            gene_id = gene_id_raw.split(".")[0]

            if gene_id.startswith("N_"):
                if verbose:
                    print(f"Skipping N_ prefixed row: {line.rstrip()}")
                continue

            if strandedness == "reverse":
                value = int(rv)
            elif strandedness == "forward":
                value = int(fw)
            else:
                value = int(both)

            total_count += value
            if gene_id in target_genes:
                target_genes_counts[gene_id] += value
            nbr_genes += 1
    return (nbr_genes, total_count, target_genes_counts)


def get_target_genes(hb_map_path: Path) -> Set[str]:
    hb_genes = set()
    with hb_map_path.open("r") as hb_in:
        hb_genes_reader = csv.DictReader(hb_in, delimiter="\t")
        for row in hb_genes_reader:
            ensembl_val = row.get("ensembl")
            if ensembl_val is None:
                raise ValueError(f"Expected ensembl value, found None for row: {row}")
            if ensembl_val in hb_genes:
                raise ValueError(f"Value {ensembl_val} has already been added to {hb_genes}")
            hb_genes.add(ensembl_val)
    return hb_genes


def parse_arguments():
    parser = argparse.ArgumentParser(
        description="Calculate percentage reads mapping to a group of genes"
    )
    parser.add_argument(
        "--version",
        action="version",
        version=VERSION,
        help="Show program's version number and exit.",
    )
    parser.add_argument("--target_genes", required=True, help="TSV file with one header 'ensembl'")
    parser.add_argument("--gene_counts", required=True, help="htseq-lib like output")
    parser.add_argument(
        "--strandedness", required=True, help="forward, reverse or none", default=None
    )
    parser.add_argument("--out_json", help="Output QC JSON path")
    parser.add_argument("--verbose", action="store_true", help="Print additional information")
    args = parser.parse_args()
    return args


if __name__ == "__main__":
    args = parse_arguments()
    out_json = Path(args.out_json) if args.out_json is not None else None
    main(Path(args.target_genes), Path(args.gene_counts), args.strandedness, out_json, args.verbose)
