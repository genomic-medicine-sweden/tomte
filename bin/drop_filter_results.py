#!/usr/bin/env python3

import argparse
import re
from collections import OrderedDict
from pathlib import Path
import csv
import os
import pandas as pd
import pyreadr


def FilteringResults(sample, gene_panel, out_drop_ae_rds, out_drop_gene_name, out_drop_as_tsv):
    # Read gene panel if it has been provided
    if gene_panel != "None":
        df_panel = pd.read_csv(gene_panel, sep="\t", skiprows=26)

    if out_drop_ae_rds != "None":
        rds_AE = pyreadr.read_r(out_drop_ae_rds)
        df_results_AE = rds_AE[None]
        # Keep only samples provided to tomte
        df_results_family_AE = df_results_AE.loc[df_results_AE["sampleID"].isin(sample)]
        # Count how many events are significant per provided sample
        df_family_AE_top20 = pd.DataFrame()

        for id in sample:
            df_id = df_results_family_AE[df_results_family_AE["sampleID"] == id]
            if sum(df_id["padjust"] < 0.05) < 20:
                df_id = df_id.sort_values(by=["pValue"]).reset_index()
                df_family_AE_top20 = pd.concat([df_family_AE_top20, df_id[0:20]], ignore_index=True, sort=False)
            else:
                df_family_AE_top20 = pd.concat([df_family_AE_top20, df_id], ignore_index=True, sort=False)
        df_family_AE_top20 = df_family_AE_top20.drop(columns=["index"])

        # Annotate with hgnc
        df_genes = pd.read_csv(out_drop_gene_name)
        df_genes.rename(columns={"gene_name": "hgncSymbol"}, inplace=True)
        df_family_annotated_AE_top20 = df_genes.merge(df_family_AE_top20, left_on="gene_id", right_on="geneID")
        df_family_annotated_AE_top20.to_csv("OUTRIDER_provided_sample_top20.tsv", sep="\t", index=False, header=True)

        # Only genes from panel
        if gene_panel != "None":
            df_clinical_AE = df_family_annotated_AE_top20.loc[
                df_family_annotated_AE_top20["hgncSymbol"].isin(df_panel["hgnc_symbol"])
            ]
            df_clinical_AE = df_panel[["hgnc_symbol", "hgnc_id"]].merge(
                df_family_annotated_AE_top20, left_on="hgnc_symbol", right_on="hgncSymbol"
            )
            df_clinical_AE = df_clinical_AE.drop(columns=["hgnc_symbol"])
            df_clinical_AE.to_csv("OUTRIDER_provided_sample_top20_filtered.tsv", sep="\t", index=False, header=True)

    if out_drop_as_tsv != "None":
        # Patient family data only
        df_results_AS = pd.read_csv(out_drop_as_tsv, sep="\t")
        df_results_family_AS = df_results_AS.loc[df_results_AS["sampleID"].isin(sample)]
        df_results_family_AS.to_csv("FRASER_provided_sample.tsv", sep="\t", index=False, header=True)

        if gene_panel != "None":
            # Only genes from panel
            df_panel = pd.read_csv(gene_panel, sep="\t", skiprows=26)
            df_clinical_AS = df_results_family_AS.loc[df_results_family_AS["hgncSymbol"].isin(df_panel["hgnc_symbol"])]
            df_clinical_AS.to_csv("FRASER_provided_sample_filtered.tsv", sep="\t", index=False, header=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        formatter_class=argparse.MetavarTypeHelpFormatter,
        description="""Filter DROP results to keep only the patient samples and the genes within the panel if provided""",
    )

    parser.add_argument(
        "--sample",
        type=str,
        nargs="+",
        help="corresponding sample name",
        required=True,
    )

    parser.add_argument(
        "--gene_panel",
        type=str,
        default="None",
        help="Path to gene panel for filtering",
        required=False,
    )

    parser.add_argument(
        "--drop_ae_rds",
        type=str,
        default="None",
        help="Path to RDS output from DROP AE",
        required=False,
    )

    parser.add_argument(
        "--out_drop_gene_name",
        type=str,
        default="None",
        help="Path to gene name annotion, output from DROP AE",
        required=False,
    )

    parser.add_argument(
        "--out_drop_as_tsv",
        type=str,
        default="None",
        help="Path to tsv output from DROP AE",
        required=False,
    )

    args = parser.parse_args()
    FilteringResults(
        args.sample,
        args.gene_panel,
        args.drop_ae_rds,
        args.out_drop_gene_name,
        args.out_drop_as_tsv,
    )
