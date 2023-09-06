#!/usr/bin/env python3

import argparse
import pandas as pd
import pyreadr

GENE_PANEL_HEADER = ["chromosome", "gene_start", "gene_stop", "hgnc_id", "hgnc_symbol"]
GENE_PANEL_COLUMNS_TO_KEEP = ["hgnc_symbol", "hgnc_id"]


def get_top_hits(id: str, df_results_family_AE: pd) -> pd:
    """
    Filter results to get only those with the id provided.
    If there are >= 20 hits it will output all of them.
    If there are less than 20 hits, it will output the 20 with the lowest p-value.
    """
    df_id: pd = df_results_family_AE[df_results_family_AE["sampleID"] == id]
    if sum(df_id["aberrant"] == "True") >= 20:
        return df_id
    df_id = df_id.sort_values(by=["pValue"]).reset_index()
    return df_id[:20]


def annotate_with_hgnc(df_family_AE_top_hits: pd, out_drop_gene_name: str) -> pd:
    """
    Annotate results from DROP AE with hgnc ids.
    """
    df_genes = pd.read_csv(out_drop_gene_name)
    df_genes.rename(columns={"gene_name": "hgncSymbol"}, inplace=True)
    return df_genes.merge(df_family_AE_top_hits, left_on="gene_id", right_on="geneID")


def keep_only_hits_in_gene_panel(df_family_top_hits: pd, gene_panel: str, module_name: str) -> pd:
    """
    Filter out from results any gene that is not present in the provided gene panel.
    """
    if gene_panel != "None":
        gene_panel_header = GENE_PANEL_HEADER
        df_panel = pd.read_csv(gene_panel, sep="\t", names=gene_panel_header, header=None, comment="#", index_col=False)
        df_clinical = df_family_top_hits.loc[df_family_top_hits["hgncSymbol"].isin(df_panel["hgnc_symbol"])]
        df_clinical = df_panel[GENE_PANEL_COLUMNS_TO_KEEP].merge(
            df_family_top_hits, left_on="hgnc_symbol", right_on="hgncSymbol"
        )
        df_clinical = df_clinical.drop(columns=["hgnc_symbol"])
        file_name = f"{module_name}_provided_samples_top_hits_filtered.tsv"
        df_clinical.to_csv(file_name, sep="\t", index=False, header=True)


def outrider_results_filter(samples: list, gene_panel: str, out_drop_ae_rds: str, out_drop_gene_name: str):
    """
    Filter results to get only those from the sample(s) provided.
    If there are >= 20 hits per sample it will output all of them.
    If there are less than 20 hits per sample, it will output the 20 with the lowest p-value.
    The results will be annotated with hgnc id.
    Two tsvs will be outputed:
        - One filtered to keep gense in the gene panel (if provided).
        - Another that is unfilterd.
    """
    rds_AE = pyreadr.read_r(out_drop_ae_rds)
    df_results_AE = rds_AE[None]
    # Keep only samples provided to tomte
    df_results_family_AE = df_results_AE.loc[df_results_AE["sampleID"].isin(samples)]
    df_family_AE_top_hits = pd.DataFrame()
    for id in samples:
        df_sample_AE_top_hits = get_top_hits(id, df_results_family_AE)
        df_family_AE_top_hits = pd.concat([df_family_AE_top_hits, df_sample_AE_top_hits], ignore_index=True, sort=False)
        df_family_AE_top_hits = df_family_AE_top_hits.drop(columns=["index"])
    df_family_annotated_AE_top_hits = annotate_with_hgnc(df_family_AE_top_hits, out_drop_gene_name)
    keep_only_hits_in_gene_panel(df_family_annotated_AE_top_hits, gene_panel, "OUTRIDER")


def fraser_results_filter(samples: list, gene_panel: str, out_drop_as_tsv: str):
    """
    Filter results to get only those from the sample(s) provided.
    Two tsvs will be outputed:
        - One filtered to keep gense in the gene panel (if provided).
        - Another that is unfilterd.
    """
    df_results_AS = pd.read_csv(out_drop_as_tsv, sep="\t")
    # Keep only samples provided to tomte
    df_results_family_AS = df_results_AS.loc[df_results_AS["sampleID"].isin(samples)]
    df_results_family_AS.to_csv("FRASER_provided_samples_top_hits.tsv", sep="\t", index=False, header=True)
    keep_only_hits_in_gene_panel(df_results_family_AS, gene_panel, "FRASER")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        formatter_class=argparse.MetavarTypeHelpFormatter,
        description="""Filter DROP results to keep only the patient samples and the genes within the panel if provided""",
    )

    parser.add_argument(
        "--samples",
        type=str,
        nargs="+",
        help="corresponding samples name",
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

    outrider_results_filter(
        args.samples,
        args.gene_panel,
        args.drop_ae_rds,
        args.out_drop_gene_name,
    )

    fraser_results_filter(
        args.samples,
        args.gene_panel,
        args.out_drop_as_tsv,
    )
