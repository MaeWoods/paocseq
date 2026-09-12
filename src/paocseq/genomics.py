"""Clonotype ranking by marker-gene UMI counts or by arbitrary metadata.

Port of ``RankTCRs`` (``Genomics.R``). The original ``GetTCRs`` (pulling
full-length TCR sequences out of a Cell Ranger ``all_contig_annotations.json``
file) and ``MakeReceptor`` (stitching full receptor sequences using live
UniProt/BLAST-style lookups over the network) are not ported: they are
one-off, network-dependent utilities orthogonal to the array/matrix-based
analysis this package focuses on, and are not needed for the feature-barcode
matrix + contig-annotation-CSV workflow described in the vignettes.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import anndata as ad

__all__ = ["rank_tcrs"]


def rank_tcrs(
    adata: ad.AnnData,
    gene: str = "IFNG",
    meta_data_col: str | None = None,
    cell_type_col: str = "cdr3_na",
) -> pd.DataFrame:
    """Rank clonotypes by the max / mean value of a gene's raw counts (or an
    arbitrary ``obs`` column) across the cells belonging to that clonotype.

    Port of ``RankTCRs`` (``Genomics.R``).

    Parameters
    ----------
    adata:
        Query AnnData with clonotype labels in ``obs[cell_type_col]``.
    gene:
        Gene to rank by (raw counts), used when ``meta_data_col`` is None.
    meta_data_col:
        If given, rank by this ``obs`` column instead of a gene's counts.
    cell_type_col:
        ``obs`` column holding the clonotype/cell-type label.

    Returns
    -------
    pandas.DataFrame
        Columns: ``clonotype``, ``clone_size``, ``max_value``, ``mean_value``,
        sorted by ``max_value`` descending. Call ``.sort_values('mean_value',
        ascending=False)`` for the R function's second (avg-ranked) table.
    """
    if meta_data_col is not None:
        if meta_data_col not in adata.obs:
            raise KeyError(meta_data_col)
        values = adata.obs[meta_data_col].to_numpy(dtype=np.float64)
    else:
        if gene not in adata.var_names:
            raise KeyError(gene)
        x = adata[:, gene].X
        values = (x.toarray() if hasattr(x, "toarray") else np.asarray(x)).ravel()

    clones = adata.obs[cell_type_col].to_numpy()
    df = pd.DataFrame({"clonotype": clones, "value": values})
    summary = (
        df.groupby("clonotype")["value"]
        .agg(clone_size="size", max_value="max", mean_value="mean")
        .reset_index()
    )
    return summary.sort_values("max_value", ascending=False).reset_index(drop=True)
