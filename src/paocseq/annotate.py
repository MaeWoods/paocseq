"""Clonotype / cell-type summary tables.

Port of ``AnnotateCellTypes`` (``Annotation.R``) -- also referred to as
``AnnotateClonotypes`` in the R package's vignette. Produces one row per
clonotype (or other cell-type label) across a set of samples, with the
percentage of CD4 and CD8 cells calling ``Threshold_<goi>`` == "high" in
each sample, plus overall clone frequency and whether the clonotype is
shared across samples.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import anndata as ad

__all__ = ["annotate_cell_types"]


def annotate_cell_types(
    adata_list: list[ad.AnnData],
    goi: str,
    cell_type_col: str = "cdr3_na",
    sample_names: list[str] | None = None,
    path: str | None = None,
) -> pd.DataFrame:
    """Build a per-clonotype summary table of gene-of-interest positivity
    across samples.

    Port of ``AnnotateCellTypes`` (``Annotation.R``). Requires
    ``adata.obs[cell_type_col]`` (clonotype id, e.g. ``'cdr3_na'``) and
    ``adata.obs[f'Threshold_{goi}']`` (``"high"``/``"unassigned"``, from
    :func:`paocseq.classify.threshold_marker_genes`) on each sample, and
    ``adata.obs['CD4cells']`` / ``adata.obs['CD8cells']`` (from
    :func:`paocseq.io.combine_data`).

    Parameters
    ----------
    adata_list:
        One ``AnnData`` per sample.
    goi:
        Gene of interest; looks for the column ``Threshold_<goi>``.
    cell_type_col:
        ``obs`` column holding the clonotype/cell-type label.
    sample_names:
        Optional display names for each entry of ``adata_list``; falls back
        to ``adata.obs['sampleref'].iloc[0]`` or ``f"sample{k}"``.
    path:
        If given, the resulting table is also written to this CSV path.

    Returns
    -------
    pandas.DataFrame
        Columns: ``clonotype``, ``shared``, ``frequency``, ``avg_clone_size``,
        and, per sample, ``<sample>_CD4_pct``, ``<sample>_CD8_pct``,
        ``<sample>_all_pct``, ``<sample>_cells``.
    """
    threshold_col = f"Threshold_{goi}"
    n_batch = len(adata_list)
    if sample_names is None:
        sample_names = []
        for k, a in enumerate(adata_list):
            if "sampleref" in a.obs and len(a.obs) > 0:
                sample_names.append(str(a.obs["sampleref"].iloc[0]))
            else:
                sample_names.append(f"sample{k + 1}")

    per_sample_clonotypes: list[set] = []
    for a in adata_list:
        vals = set(a.obs[cell_type_col].astype(str).unique()) - {"unassigned"}
        per_sample_clonotypes.append(vals)

    all_clonotypes: set = set()
    intersect_clonotypes: set | None = None
    for vals in per_sample_clonotypes:
        all_clonotypes |= vals
        intersect_clonotypes = vals if intersect_clonotypes is None else (intersect_clonotypes & vals)
    intersect_clonotypes = intersect_clonotypes or set()

    # Clone frequency: first non-zero countcln value found for the clonotype
    # across samples (matches the R `Clonefreq` assignment loop).
    freq_lookup: dict[str, int] = {}
    for a in adata_list:
        if "countcln" not in a.obs:
            continue
        for ct, size in zip(a.obs[cell_type_col], a.obs["countcln"]):
            if ct not in freq_lookup or freq_lookup[ct] == 0:
                freq_lookup[ct] = int(size)

    clonotypes_sorted = sorted(all_clonotypes, key=lambda c: freq_lookup.get(c, 0), reverse=True)

    rows = []
    for clone in clonotypes_sorted:
        row: dict = {
            "clonotype": clone,
            "shared": "yes" if clone in intersect_clonotypes else "no",
            "frequency": freq_lookup.get(clone, 0),
        }
        clone_sizes = []
        for k, a in enumerate(adata_list):
            name = sample_names[k]
            obs = a.obs
            in_clone = obs[cell_type_col] == clone
            n_total = int(in_clone.sum())
            clone_sizes.append(n_total)

            has_thresh = threshold_col in obs
            if n_total == 0:
                cd4_pct = cd8_pct = all_pct = 0.0
            else:
                if has_thresh:
                    high_mask = obs[threshold_col] == "high"
                else:
                    high_mask = pd.Series(False, index=obs.index)

                cd4_mask = obs.get("CD4cells", pd.Series(0, index=obs.index)) == 1
                cd8_mask = obs.get("CD8cells", pd.Series(0, index=obs.index)) == 1

                n_cd4 = int((in_clone & cd4_mask).sum())
                n_cd8 = int((in_clone & cd8_mask).sum())
                n_cd4_high = int((in_clone & cd4_mask & high_mask).sum())
                n_cd8_high = int((in_clone & cd8_mask & high_mask).sum())
                n_all_high = int((in_clone & high_mask).sum())

                cd4_pct = 100.0 * n_cd4_high / n_cd4 if n_cd4 > 0 else 0.0
                cd8_pct = 100.0 * n_cd8_high / n_cd8 if n_cd8 > 0 else 0.0
                all_pct = 100.0 * n_all_high / n_total if n_total > 0 else 0.0

            row[f"{name}_CD4_pct"] = round(cd4_pct, 2)
            row[f"{name}_CD8_pct"] = round(cd8_pct, 2)
            row[f"{name}_all_pct"] = round(all_pct, 2)
            row[f"{name}_cells"] = n_total

        row["avg_clone_size"] = float(np.ceil(np.mean(clone_sizes))) if clone_sizes else 0.0
        rows.append(row)

    df = pd.DataFrame(rows)
    if len(df):
        df = df.sort_values("avg_clone_size", ascending=False).reset_index(drop=True)

    if path is not None:
        df.to_csv(path, index=False)

    return df
