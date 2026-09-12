"""Plotting helpers, port of ``Plotting.R``.

``SegmentPlot`` in R drew a circos/chord diagram (via the ``circlize``
package) comparing top clonotype frequencies between two samples; here
:func:`segment_plot` draws the same comparison as a paired horizontal bar
chart with matplotlib, to avoid adding a circos-plotting dependency.
``UMAPReduce`` (Seurat integration + UMAP) becomes :func:`umap_reduce`
using the equivalent Scanpy pipeline. ``QCPlot`` becomes :func:`qc_plot`.
``SaveHeatmap`` (which in R hard-codes indices/column counts for one
specific dataset) becomes a general-purpose :func:`expression_heatmap`.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import anndata as ad
import matplotlib.pyplot as plt
import scanpy as sc

__all__ = ["qc_plot", "umap_reduce", "segment_plot", "expression_heatmap"]


def qc_plot(
    adata: ad.AnnData,
    n_feature_lower: float | None = None,
    n_feature_upper: float | None = None,
    percent_mt_upper: float | None = None,
    ax=None,
):
    """Three-panel QC violin plot: n_features, total counts, percent
    mitochondrial, with optional cutoff lines. Port of ``QCPlot``
    (``Plotting.R``)."""
    x = adata.X.toarray() if hasattr(adata.X, "toarray") else np.asarray(adata.X)
    n_features = (x > 0).sum(axis=1)
    n_counts = x.sum(axis=1)
    pct_mt = adata.obs["percent_mt"] if "percent_mt" in adata.obs else np.zeros(adata.n_obs)

    fig, axes = plt.subplots(1, 3, figsize=(10, 4)) if ax is None else (ax[0].figure, ax)
    for a, values, title, cutoff in zip(
        axes,
        [n_features, n_counts, pct_mt],
        ["Features per cell", "Counts per cell", "Percent mitochondria"],
        [(n_feature_lower, n_feature_upper), (None, None), (None, percent_mt_upper)],
    ):
        a.violinplot(values, showmeans=False, showmedians=True)
        for c in cutoff:
            if c is not None:
                a.axhline(c, color="black", linewidth=1)
        a.set_title(title, fontsize=9)
        a.set_xticks([])
    fig.tight_layout()
    return fig


def umap_reduce(
    adata_list: list[ad.AnnData],
    color_by: list[str] | None = None,
    n_top_genes: int = 2000,
    n_pcs: int = 30,
    random_state: int = 356,
    save_path: str | None = None,
):
    """Normalize, integrate (batch-corrected concatenation if more than one
    sample), and UMAP-embed a list of ``AnnData`` samples, plotting
    ``color_by`` metadata columns.

    Port of ``UMAPReduce`` (``Plotting.R``): the original used Seurat's
    CCA-based ``FindIntegrationAnchors``/``IntegrateData``; here, samples are
    concatenated and Scanpy's Harmony-free PCA + neighbors + UMAP pipeline is
    used directly (call :func:`scanpy.external.pp.harmony_integrate` on the
    result yourself first if you need explicit batch correction).

    Returns
    -------
    anndata.AnnData
        The combined, embedded object (``.obsm['X_umap']`` populated).
    """
    color_by = color_by or ["sampleref"]
    if len(adata_list) > 1:
        combined = ad.concat(adata_list, join="inner", label="batch", keys=[str(i) for i in range(len(adata_list))])
    else:
        combined = adata_list[0].copy()

    sc.pp.normalize_total(combined, target_sum=1e4)
    sc.pp.log1p(combined)
    sc.pp.highly_variable_genes(combined, n_top_genes=n_top_genes)
    combined_hvg = combined[:, combined.var["highly_variable"]].copy()
    sc.pp.scale(combined_hvg, max_value=10)
    sc.tl.pca(combined_hvg, n_comps=min(n_pcs, combined_hvg.n_obs - 1, combined_hvg.n_vars - 1), random_state=random_state)
    sc.pp.neighbors(combined_hvg, random_state=random_state)
    sc.tl.umap(combined_hvg, random_state=random_state)
    combined.obsm["X_umap"] = combined_hvg.obsm["X_umap"]
    combined.obsm["X_pca"] = combined_hvg.obsm["X_pca"]

    present = [c for c in color_by if c in combined.obs]
    if present:
        fig = sc.pl.umap(combined, color=present, show=False, return_fig=True)
        if save_path:
            fig.savefig(save_path, dpi=150, bbox_inches="tight")

    return combined


def segment_plot(
    summary_table: pd.DataFrame,
    sample_columns: tuple[str, str],
    segment_names: tuple[str, str] = ("S1", "S2"),
    n_segments: int = 20,
    save_path: str | None = None,
):
    """Paired top-N clonotype-frequency bar chart comparing two samples.

    Matplotlib re-implementation of ``SegmentPlot`` (``Plotting.R``), which
    drew a circos/chord diagram via the ``circlize`` R package; the
    frequency-based content (top-N clonotypes per sample, ranked, normalized
    to proportions) is unchanged, only the chart type differs so this
    package doesn't need to add a circos-plotting dependency.

    Parameters
    ----------
    summary_table:
        A table like the output of :func:`paocseq.annotate.annotate_cell_types`,
        with a ``clonotype`` column and one abundance/frequency column per
        sample.
    sample_columns:
        The two column names in ``summary_table`` to compare.
    """
    col1, col2 = sample_columns
    df = summary_table[(summary_table[col1] > 0) | (summary_table[col2] > 0)].copy()
    df[col1] = df[col1] / df[col1].sum()
    df[col2] = df[col2] / df[col2].sum()

    top1 = df.nlargest(n_segments, col1)[["clonotype", col1]]
    top2 = df.nlargest(n_segments, col2)[["clonotype", col2]]

    fig, axes = plt.subplots(1, 2, figsize=(10, 6), sharey=False)
    axes[0].barh(top1["clonotype"].astype(str), top1[col1], color="green")
    axes[0].invert_yaxis()
    axes[0].set_title(segment_names[0])
    axes[0].set_xlabel("Proportion of cells")

    axes[1].barh(top2["clonotype"].astype(str), top2[col2], color="blue")
    axes[1].invert_yaxis()
    axes[1].set_title(segment_names[1])
    axes[1].set_xlabel("Proportion of cells")

    fig.tight_layout()
    if save_path:
        fig.savefig(save_path, dpi=150, bbox_inches="tight")
    return fig


def expression_heatmap(
    adata: ad.AnnData,
    gene_list: list[str],
    group_col: str = "cdr3_na",
    n_groups: int = 20,
    save_path: str | None = None,
):
    """Gene x clonotype mean-expression heatmap.

    General-purpose replacement for ``SaveHeatmap`` (``Plotting.R``), which
    hard-coded dataset-specific column ranges and gene-panel sizes; here the
    top ``n_groups`` most abundant groups (by cell count) in
    ``obs[group_col]`` are shown against ``gene_list``.
    """
    present_genes = [g for g in gene_list if g in adata.var_names]
    x = adata[:, present_genes].X
    x = x.toarray() if hasattr(x, "toarray") else np.asarray(x)
    df = pd.DataFrame(x, columns=present_genes)
    df[group_col] = adata.obs[group_col].to_numpy()

    top_groups = df[group_col].value_counts().nlargest(n_groups).index
    means = df[df[group_col].isin(top_groups)].groupby(group_col)[present_genes].mean()
    means = means.loc[top_groups]

    fig, ax = plt.subplots(figsize=(max(6, len(present_genes) * 0.5), max(4, n_groups * 0.3)))
    im = ax.imshow(means.to_numpy(), aspect="auto", cmap="viridis")
    ax.set_xticks(range(len(present_genes)))
    ax.set_xticklabels(present_genes, rotation=90, fontsize=8)
    ax.set_yticks(range(len(means)))
    ax.set_yticklabels(means.index, fontsize=8)
    fig.colorbar(im, ax=ax, label="Mean expression")
    fig.tight_layout()
    if save_path:
        fig.savefig(save_path, dpi=150, bbox_inches="tight")
    return fig
