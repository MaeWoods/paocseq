"""Distance-based cell classification against a reference signature.

Port of ``Distances.R``. ``classify_cells`` mirrors ``ClassifyCells``: it
normalizes a query ``AnnData``'s raw counts (proportional-fitting + log1p,
i.e. ``PFlog1pPF``, matching the R code exactly), restricts to a gene
signature, and scores each cell against a reference matrix using one of:

* ``distance=0`` -- Mahalanobis distance (:func:`paocseq.._solvers.mahalanobis_distance`)
* ``distance=1`` -- mean taxicab (L1) distance to every reference cell
* ``distance=2`` -- z-score distance (single-gene signature only)
* ``distance=3`` -- isolation-forest outlier score, aggregated per clonotype
  (:func:`paocseq.._solvers.isolation_forest_score`)
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import anndata as ad

from . import _solvers

__all__ = ["pf_log1p_pf", "classify_cells", "add_distances", "percent_outlier"]


def pf_log1p_pf(counts: np.ndarray) -> np.ndarray:
    """Proportional-fitting + log1p + proportional-fitting normalization.

    Direct port of the repeated::

        PF = log(1 + t(t(counts) / mean(colSums(counts))))
        PFlog1pPF = t(t(PF) / mean(colSums(PF)))

    block used throughout ``Distances.R``. ``counts`` is genes x cells.
    """
    counts = np.asarray(counts, dtype=np.float64)
    col_sums = counts.sum(axis=0)
    mean_col_sum = col_sums.mean() if col_sums.mean() != 0 else 1.0
    pf = np.log1p(counts / mean_col_sum)
    pf_col_sums = pf.sum(axis=0)
    mean_pf_col_sum = pf_col_sums.mean() if pf_col_sums.mean() != 0 else 1.0
    return pf / mean_pf_col_sum


def _gene_panel(adata: ad.AnnData, gene_list: list[str]) -> tuple[np.ndarray, list[str]]:
    """Raw counts (genes x cells) restricted to ``gene_list``, NaN-safe."""
    present = [g for g in gene_list if g in adata.var_names]
    x = adata[:, present].X
    x = x.toarray() if hasattr(x, "toarray") else np.asarray(x)
    return x.T, present


def classify_cells(
    adata: ad.AnnData,
    reference_data: pd.DataFrame,
    gene_list: list[str],
    cell_types: pd.Series | None = None,
    distance: int = 0,
    scramble: bool = False,
    n_trees: int = 10,
    max_height: int = 20,
    random_state: int | None = None,
) -> ad.AnnData:
    """Score each cell in ``adata`` against ``reference_data`` for a gene
    signature, adding the result as ``obs`` metadata. Port of
    ``ClassifyCells`` (``Distances.R``).

    Parameters
    ----------
    adata:
        Query AnnData; raw counts are used (``adata.X``), matching the R
        function's use of ``@assays$RNA@counts``.
    reference_data:
        ``(genes, reference_cells)`` DataFrame, gene-indexed the same as
        ``gene_list``.
    gene_list:
        Gene signature to score on.
    cell_types:
        Per-cell clonotype/cell-type labels (required for ``distance=3``);
        if not given, ``adata.obs['cdr3_na']`` is used when present.
    distance:
        ``0`` Mahalanobis, ``1`` taxicab, ``2`` z-score (first gene of
        ``gene_list`` only, matching the R behaviour), ``3`` isolation
        forest (mean outlier score per clonotype).
    scramble:
        If True, score against a random set of genes of the same size
        instead of ``gene_list`` (used in the R package to sanity-check the
        signature's specificity).

    Returns
    -------
    anndata.AnnData
        ``adata`` with an added ``obs`` column (``'Mdist'``, ``'Mdistscr'``,
        or ``'Isoforest'`` depending on ``distance``/``scramble``).
    """
    rng = np.random.default_rng(random_state)
    counts, present = _gene_panel(adata, gene_list)
    normalized = pf_log1p_pf(counts)
    normalized = np.nan_to_num(normalized, nan=0.0)

    ref_genes = [g for g in gene_list if g in reference_data.index]
    ref_panel = reference_data.loc[ref_genes].to_numpy(dtype=np.float64)

    if scramble:
        n_pick = len(present)
        idx = rng.integers(0, normalized.shape[0], size=n_pick)
        test_panel = normalized[idx, :]
    else:
        test_panel = normalized

    if distance == 0:
        dist = _solvers.mahalanobis_distance(ref_panel, test_panel)
        col = "Mdistscr" if scramble else "Mdist"
        adata.obs[col] = dist
    elif distance == 1:
        # mean absolute difference to every reference cell, per gene, averaged
        dist = np.mean(
            np.abs(ref_panel[:, :, None] - test_panel[:, None, :]).mean(axis=1), axis=0
        )
        col = "Mdistscr" if scramble else "Mdist"
        adata.obs[col] = dist
    elif distance == 2:
        sig = ref_panel[0, :]
        mean_sig, std_sig = sig.mean(), sig.std(ddof=1)
        dist = np.abs(test_panel[0, :] - mean_sig) / std_sig
        col = "Mdistscr" if scramble else "Mdist"
        adata.obs[col] = dist
    elif distance == 3:
        if cell_types is None:
            if "cdr3_na" not in adata.obs:
                raise ValueError("distance=3 requires cell_types or adata.obs['cdr3_na']")
            cell_types = adata.obs["cdr3_na"]
        scores = np.full(adata.n_obs, np.nan)
        for ctype in pd.unique(cell_types):
            mask = (cell_types == ctype).to_numpy()
            if not mask.any():
                continue
            query = test_panel[:, mask]
            per_cell = _solvers.isolation_forest_score(
                ref_panel, query, n_trees=n_trees, max_height=max_height, random_state=random_state
            )
            scores[mask] = per_cell.mean()
        adata.obs["Isoforest"] = scores
    else:
        raise ValueError("distance must be 0, 1, 2, or 3")

    return adata


def add_distances(
    adata: ad.AnnData,
    signature_ref: pd.DataFrame,
    gene_list: list[str],
    cell_types: pd.Series | None = None,
    distance: int = 0,
    scramble: bool = False,
) -> ad.AnnData:
    """Convenience wrapper matching ``AddDistances`` (``Distances.R``):
    normalizes ``adata`` and calls :func:`classify_cells`."""
    return classify_cells(
        adata,
        reference_data=signature_ref,
        gene_list=gene_list,
        cell_types=cell_types,
        distance=distance,
        scramble=scramble,
    )


def percent_outlier(
    query_cells: np.ndarray,
    reference_data: np.ndarray,
    n_trees: int = 10,
    max_height: int = 20,
    threshold: float = 0.75,
    random_state: int | None = None,
) -> dict:
    """Fraction of ``query_cells`` flagged as outliers relative to
    ``reference_data`` by the isolation-forest score. Port of
    ``PercentOutlier`` (``Distances.R``), using the vectorized
    :func:`paocseq._solvers.isolation_forest_score` in place of the R
    function's per-cell for-loop.

    Parameters
    ----------
    query_cells, reference_data:
        ``(genes, cells)`` arrays.
    threshold:
        Normalized isolation score above which a cell is called an outlier.

    Returns
    -------
    dict
        ``{'clone_height': ..., 'clone_AS': ..., 'outlier_fraction': ...}``
    """
    scores = _solvers.isolation_forest_score(
        reference_data, query_cells, n_trees=n_trees, max_height=max_height, random_state=random_state
    )
    return {
        "clone_AS": float(np.mean(scores)),
        "outlier_fraction": float(np.mean(scores > threshold)),
    }
