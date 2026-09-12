"""Subsetting cells by threshold/phenotype/clonotype, and differential
expression between "specific" vs "bystander" groups.

Port of ``Expression.R``.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import anndata as ad
import scanpy as sc

__all__ = ["get_specific_cells", "get_gene_signature"]


def get_specific_cells(
    adata: ad.AnnData,
    expression: str,
    phenotype: str | None,
    cell_types: list[str],
    goi: str,
    cell_type_col: str = "cdr3_na",
) -> ad.AnnData:
    """Subset ``adata`` to cells with a given clonotype/cell-type
    membership, expression threshold, and (optionally) CD4/CD8 phenotype.

    Port of ``GetSpecificCells`` (``Expression.R``).

    Parameters
    ----------
    expression:
        Value to match in ``obs[f'Threshold_{goi}']`` (typically ``"high"``
        or ``"unassigned"``).
    phenotype:
        ``"CD4"``, ``"CD8"``, or ``None`` to skip phenotype filtering.
    cell_types:
        Clonotype/cell-type labels to include.
    goi:
        Gene of interest; looks up ``obs[f'Threshold_{goi}']``.

    Returns
    -------
    anndata.AnnData
        Subset of ``adata`` (a view is *not* returned; a ``.copy()`` is).
    """
    threshold_col = f"Threshold_{goi}"
    if threshold_col not in adata.obs:
        raise KeyError(f"{threshold_col!r} not found in obs; run threshold_marker_genes first")

    mask = adata.obs[cell_type_col].isin(cell_types) & (adata.obs[threshold_col] == expression)
    if phenotype == "CD4":
        mask &= adata.obs.get("CD4cells", 0) == 1
    elif phenotype == "CD8":
        mask &= adata.obs.get("CD8cells", 0) == 1

    return adata[mask].copy()


def get_gene_signature(
    adata: ad.AnnData,
    goi: str,
    specific_clonotypes: dict[str, list[str]],
    bystander_clonotypes: dict[str, list[str]],
    cell_type_col: str = "cdr3_na",
    fc_lim: float = 0.0,
    min_pct: float = 0.25,
    pval_cutoff: float = 0.05,
) -> dict[str, pd.DataFrame]:
    """Differential expression between "specific" and "bystander" CD4/CD8
    subsets, for a gene-of-interest signature.

    Simplified port of ``GetGeneSignature`` (``Expression.R``): rather than
    the original's bespoke parsing of a multi-condition annotation
    spreadsheet to decide which clonotypes are "specific" vs "bystander",
    the caller passes those two clonotype lists directly (e.g. taken from
    :func:`paocseq.annotate.annotate_cell_types`), for each of ``"CD4"`` and
    ``"CD8"``. Uses :func:`scanpy.tl.rank_genes_groups` (Wilcoxon, matching
    the spirit of the original's ``FindMarkers(test.use="bimod")`` call) in
    place of Seurat's ``FindMarkers``.

    Parameters
    ----------
    specific_clonotypes, bystander_clonotypes:
        ``{"CD4": [...], "CD8": [...]}``-style dicts of clonotype labels.

    Returns
    -------
    dict[str, pandas.DataFrame]
        ``{"CD4": df, "CD8": df}`` (only for phenotypes present in both
        input dicts), each with columns ``gene``, ``log2fc``, ``pval_adj``.
    """
    results: dict[str, pd.DataFrame] = {}
    for pheno in ("CD4", "CD8"):
        specific = specific_clonotypes.get(pheno, [])
        bystander = bystander_clonotypes.get(pheno, [])
        if not specific or not bystander:
            continue

        pheno_col = f"{pheno}cells"
        group_col = "__paocseq_specificity_group__"
        groups = pd.Series("other", index=adata.obs_names)
        in_pheno = adata.obs.get(pheno_col, 0) == 1
        groups[in_pheno & adata.obs[cell_type_col].isin(specific)] = "specific"
        groups[in_pheno & adata.obs[cell_type_col].isin(bystander)] = "bystander"
        adata.obs[group_col] = groups.values

        if (groups == "specific").sum() < 2 or (groups == "bystander").sum() < 2:
            del adata.obs[group_col]
            continue

        sc.tl.rank_genes_groups(
            adata, groupby=group_col, groups=["specific"], reference="bystander", method="wilcoxon"
        )
        names = np.asarray(adata.uns["rank_genes_groups"]["names"]["specific"])
        lfc = np.asarray(adata.uns["rank_genes_groups"]["logfoldchanges"]["specific"])
        padj = np.asarray(adata.uns["rank_genes_groups"]["pvals_adj"]["specific"])
        del adata.obs[group_col]

        df = pd.DataFrame({"gene": names, "log2fc": lfc, "pval_adj": padj})
        df = df[(df["log2fc"].abs() > fc_lim) & (df["pval_adj"] < pval_cutoff)]
        df = df.sort_values("log2fc", ascending=False).reset_index(drop=True)
        results[pheno] = df

    return results
