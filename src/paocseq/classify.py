"""Marker-gene thresholding and clonotype/cell-type classification.

Two pieces of the original R package live here:

1. Marker-gene "high"/"unassigned" thresholding -- the section embedded in
   ``CombineData`` (``Import.R``, ``preset``/``threshold.cutoff`` logic) that
   produces ``Threshold_<gene>`` metadata columns, factored out here as
   :func:`threshold_marker_genes` so it can be re-run independently (this is
   also what the R vignette calls ``ClassifyClonotypes``).
2. :func:`classify_cell_types` -- a simplified, dependency-light stand-in
   for ``ClassifyCellTypes`` (``Classification.R``). The original supports a
   ``method="ZINB"`` option built on a vendored copy of the Bioconductor
   ``DEsingle`` package (a zero-inflated negative-binomial DE test); rather
   than re-implement that whole statistical model, ``method="wilcoxon"``
   here uses :func:`scanpy.tl.rank_genes_groups` (Wilcoxon rank-sum test) to
   get the same "is this clonotype specifically expressing the gene of
   interest" answer with a well-tested, off-the-shelf test. The
   ``"mahalanobis"``/``"taxicab"``/``"zscore"``/``"isoforest"`` methods are
   unchanged and call straight through to :mod:`paocseq.distances`.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import anndata as ad
import scanpy as sc

from . import distances as _distances

__all__ = ["threshold_marker_genes", "classify_cell_types"]


def threshold_marker_genes(
    adata: ad.AnnData,
    marker_genes: list[str],
    layer: str | None = None,
    threshold_cutoff: float = 0.975,
    preset: int = 1,
    threshold_entry: float | list[float] | None = None,
    mask: list[str] | None = None,
    high_expr_genes: list[str] | None = None,
) -> ad.AnnData:
    """Add a ``Threshold_<gene>`` "high"/"unassigned" metadata column per
    marker gene.

    Port of the per-sample thresholding block embedded in ``CombineData``
    (``Import.R``): normalized expression of each marker gene is compared to
    a per-gene cutoff, and every cell above the cutoff is labelled
    ``"high"``.

    Parameters
    ----------
    adata:
        Query AnnData. Uses ``adata.X`` (or ``adata.layers[layer]`` if
        given) as the expression matrix to threshold against directly
        (equivalent to the R code's ``SCT`` slot -- pass in an already
        log-normalized layer, e.g. via ``scanpy.pp.normalize_total`` +
        ``scanpy.pp.log1p``, for comparable behaviour).
    marker_genes:
        Genes to threshold.
    threshold_cutoff:
        Quantile (0-1) above which a cell is called "high" for a gene, when
        ``preset in (1, 2)``.
    preset:
        ``1``: per-gene quantile cutoff computed on all cells in ``adata``
        (the "control is part of the dataset" case).
        ``0``: a single fixed cutoff (``threshold_entry``, a scalar or a
        length-1 list) applied to every marker gene.
        ``2``: fixed per-gene cutoffs, one entry of ``threshold_entry`` per
        gene.
    mask:
        Subset of ``marker_genes`` that should instead be thresholded on raw
        counts (equivalent to the R code's UMI/RNA-slot branch, used e.g.
        for surface-protein markers that shouldn't be SCT-normalized).
    high_expr_genes:
        Deprecated alias for ``mask`` kept for readability; if both are
        given, ``mask`` wins.

    Returns
    -------
    anndata.AnnData
        ``adata`` with new ``obs['Threshold_<gene>']`` columns.
    """
    mask = mask if mask is not None else (high_expr_genes or [])
    x = adata.X if layer is None else adata.layers[layer]
    x = x.toarray() if hasattr(x, "toarray") else np.asarray(x)
    raw = adata.X
    raw = raw.toarray() if hasattr(raw, "toarray") else np.asarray(raw)

    for i, gene in enumerate(marker_genes):
        if gene not in adata.var_names:
            continue
        gi = adata.var_names.get_loc(gene)
        use_raw = gene in mask
        values = raw[:, gi] if use_raw else x[:, gi]

        if preset == 1:
            cutoff = np.quantile(values, threshold_cutoff)
        elif preset == 0:
            if threshold_entry is None:
                raise ValueError("preset=0 requires threshold_entry")
            cutoff = threshold_entry if np.isscalar(threshold_entry) else threshold_entry[0]
        elif preset == 2:
            if threshold_entry is None:
                raise ValueError("preset=2 requires threshold_entry (one value per gene)")
            cutoff = threshold_entry[i]
        else:
            raise ValueError("preset must be 0, 1, or 2")

        labels = np.where(values > cutoff, "high", "unassigned")
        adata.obs[f"Threshold_{gene}"] = labels

    return adata


def _cd4_cd8_majority(adata: ad.AnnData, cell_type_col: str, clone: str, min_cells: int = 5) -> str:
    sub = adata.obs[adata.obs[cell_type_col] == clone]
    if len(sub) <= min_cells:
        return "-"
    cd4 = int((sub.get("CD4cells", pd.Series(dtype=int)) == 1).sum())
    cd8 = int((sub.get("CD8cells", pd.Series(dtype=int)) == 1).sum())
    if cd4 > cd8:
        return "CD4"
    if cd8 > cd4:
        return "CD8"
    return "-"


def classify_cell_types(
    adata: ad.AnnData,
    clonotypes: list[str],
    cell_type_col: str = "cdr3_na",
    method: str = "wilcoxon",
    goi: str = "IFNG",
    percentile: float = 0.01,
    reference_data: pd.DataFrame | None = None,
    gene_list: list[str] | None = None,
) -> pd.DataFrame:
    """Classify each clonotype/cell type in ``clonotypes`` as CD4/CD8 and
    (optionally) as specifically responding to ``goi``.

    Simplified port of ``ClassifyCellTypes`` (``Classification.R``); see the
    module docstring for how ``method`` maps onto the original R
    ``method="ZINB"|"Mahalanobis"|"Taxi cab"|"z-score"|"IsoForest"`` options.

    Returns
    -------
    pandas.DataFrame
        One row per clonotype, with columns ``cell_type`` (the clonotype
        label), ``phenotype`` (``CD4``/``CD8``/``-``), and, when
        ``method="wilcoxon"``, ``pvalue``/``log2fc`` for ``goi`` in that
        clonotype vs. all other cells.
    """
    rows = []
    for clone in clonotypes:
        phenotype = _cd4_cd8_majority(adata, cell_type_col, clone)
        row = {"cell_type": clone, "phenotype": phenotype}

        if method == "wilcoxon":
            if goi not in adata.var_names:
                row["pvalue"] = np.nan
                row["log2fc"] = np.nan
            else:
                group_col = "__paocseq_group__"
                adata.obs[group_col] = np.where(adata.obs[cell_type_col] == clone, "clone", "rest")
                if (adata.obs[group_col] == "clone").sum() >= 2 and (adata.obs[group_col] == "rest").sum() >= 2:
                    sc.tl.rank_genes_groups(
                        adata, groupby=group_col, groups=["clone"], reference="rest", method="wilcoxon"
                    )
                    names = adata.uns["rank_genes_groups"]["names"]["clone"]
                    idx = np.where(np.asarray(names) == goi)[0]
                    if idx.size:
                        i = idx[0]
                        row["pvalue"] = float(adata.uns["rank_genes_groups"]["pvals_adj"]["clone"][i])
                        row["log2fc"] = float(adata.uns["rank_genes_groups"]["logfoldchanges"]["clone"][i])
                    else:
                        row["pvalue"] = np.nan
                        row["log2fc"] = np.nan
                else:
                    row["pvalue"] = np.nan
                    row["log2fc"] = np.nan
                del adata.obs[group_col]

        elif method in ("mahalanobis", "taxicab", "zscore", "isoforest"):
            if reference_data is None or gene_list is None:
                raise ValueError(f"method={method!r} requires reference_data and gene_list")
            dist_code = {"mahalanobis": 0, "taxicab": 1, "zscore": 2, "isoforest": 3}[method]
            sub = adata[adata.obs[cell_type_col] == clone].copy()
            if sub.n_obs == 0:
                row["distance"] = np.nan
            else:
                sub = _distances.classify_cells(sub, reference_data, gene_list, distance=dist_code)
                col = {"mahalanobis": "Mdist", "taxicab": "Mdist", "zscore": "Mdist", "isoforest": "Isoforest"}[
                    method
                ]
                row["distance"] = float(np.nanmean(sub.obs[col]))
        else:
            raise ValueError(f"Unknown method: {method!r}")

        rows.append(row)

    return pd.DataFrame(rows)
