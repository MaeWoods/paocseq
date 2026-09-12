"""Data import: 10x Genomics feature-barcode matrices + VDJ contig
annotations, combined and QC'd into ``AnnData`` objects.

Port of ``CombineData`` in ``Import.R``. Where the R version returned a list
of Seurat objects, :func:`combine_data` returns a list of ``AnnData``
objects (one per sample, or one per hashtag if ``demultiplex=True``) -- the
AnnData/Scanpy stack is the direct analogue of the Seurat object used
throughout the rest of aocseq.

Each returned ``AnnData`` has:

* ``.X`` -- raw counts (genes as ``.var``, cells as ``.obs``), same as
  Seurat's ``RNA`` counts assay.
* ``.obs['CD4cells']`` / ``.obs['CD8cells']`` -- 0/1 indicator from raw CD4 /
  CD8A counts (mutually exclusive positivity), same logic as the R
  ``CD4cells`` / ``CD8cells`` metadata columns.
* ``.obs['sampleref']`` -- the sample name.
* ``.obs['clonotype']`` / ``.obs['countcln']`` / ``.obs['cdr3_na']`` /
  ``.obs['cdr3']`` -- clonotype id, clone size, nucleotide CDR3 and amino
  acid CDR3 of the dominant (highest-read) productive TRB chain for that
  cell, when ``vdj_path`` is supplied.
"""

from __future__ import annotations

from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc

from . import demux as _demux

__all__ = ["read_10x_sample", "combine_data"]


def read_10x_sample(path: str | Path, sample_name: str | None = None) -> ad.AnnData:
    """Read a single 10x Genomics ``filtered_feature_bc_matrix`` directory
    (or ``.h5``) into an ``AnnData`` object. Thin wrapper around
    :func:`scanpy.read_10x_mtx` / :func:`scanpy.read_10x_h5`, analogous to
    ``Seurat::Read10X`` + ``CreateSeuratObject`` in the R package.
    """
    path = Path(path)
    if path.is_file() and path.suffix == ".h5":
        adata = sc.read_10x_h5(path)
    else:
        adata = sc.read_10x_mtx(path, var_names="gene_symbols", cache=False)
    adata.var_names_make_unique()
    adata.obs_names_make_unique()
    if sample_name is not None:
        adata.obs["sampleref"] = sample_name
    return adata


def _infer_sample_name(path: str | Path) -> str:
    """Port of the trailing-path-component parser in ``CombineData`` --
    equivalent to ``os.path.basename(path.rstrip('/'))``."""
    return Path(str(path).rstrip("/\\")).name


def _percent_mito(adata: ad.AnnData) -> np.ndarray:
    mito = adata.var_names.str.upper().str.startswith("MT-")
    if not mito.any():
        return np.zeros(adata.n_obs)
    x = adata[:, mito].X
    x = x.toarray() if hasattr(x, "toarray") else np.asarray(x)
    total = adata.X.toarray() if hasattr(adata.X, "toarray") else np.asarray(adata.X)
    denom = total.sum(axis=1)
    with np.errstate(divide="ignore", invalid="ignore"):
        pct = np.where(denom > 0, 100.0 * x.sum(axis=1) / denom, 0.0)
    return np.asarray(pct).ravel()


def _qc_filter(
    adata: ad.AnnData,
    n_feature_lower_q: float,
    n_feature_upper_q: float,
    percent_mt_upper_q: float,
) -> ad.AnnData:
    """Quantile-based QC filter, matching the ``upperQ``/``lowerQ`` branch of
    ``CombineData`` (default behaviour when ``QC_plots`` is not used to
    supply fixed cutoffs)."""
    #x = adata.X.toarray() if hasattr(adata.X, "toarray") else np.asarray(adata.X)
    #print("got here 1")
    #print(adata)
    #commenting out this section of the Sonnet port due to it failing for
    #large values of n_features. Replacing with a method that works for typical 
    #single cell data that was not provided
    #n_features = (x > 0).sum(axis=1)
    
    chunksize=10000
    chunks=[]
    
    for i in range(0,adata.shape[0],chunksize):
        chunk_view=adata[i:i+chunksize]
        chunk=chunk_view.copy()
        x = chunk.X.toarray() if hasattr(chunk.X, "toarray") else np.asarray(chunk.X)
        n_features = (x > 0).sum(axis=1)
        pct_mt = _percent_mito(chunk)
        chunk.obs["n_features"] = n_features
        chunk.obs["percent_mt"] = pct_mt

        lo = np.quantile(n_features, n_feature_lower_q)
        hi = np.quantile(n_features, n_feature_upper_q)
        mt_hi = np.quantile(pct_mt, percent_mt_upper_q)

        keep = (n_features > lo) & (n_features < hi) & (pct_mt < mt_hi)
        chunks.append(chunk[keep].copy())
        
    adata_patched = ad.concat(chunks, axis=0, join="outer")
    
    #lo = np.percentile(adata.X.toarray(), 5, axis=1, keepdims=True)
    #hi = np.percentile(adata.X.toarray(), 95, axis=1, keepdims=True)
    
    #mask=(adata.X >= lo) & (adata.X <= hi)
    #if isinstance(adata.X,np.array):
    #    adata.X[~mask]=0
    #else:
    #    adata.X=adata.X.multiply(mask).tocsr()
    #print("got here 1a")
    #pct_mt = _percent_mito(adata)
    #print("got here 1b")
    #adata.obs["n_features"] = n_features
    #adata.obs["percent_mt"] = pct_mt
    #print("got here 2")

    #lo = np.quantile(n_features, n_feature_lower_q)
    #hi = np.quantile(n_features, n_feature_upper_q)
    #mt_hi = np.quantile(pct_mt, percent_mt_upper_q)
    #print("got here 3")
    #keep = (n_features > lo) & (n_features < hi) & (pct_mt < mt_hi)
    return adata_patched.copy()


def _add_cd4_cd8_indicator(adata: ad.AnnData) -> None:
    """Mutually-exclusive CD4/CD8A raw-count positivity indicator, matching
    the ``CD4cells``/``CD8cells`` loop in ``CombineData``."""
    x = adata.X.toarray() if hasattr(adata.X, "toarray") else np.asarray(adata.X)
    genes = adata.var_names
    cd4_idx = genes.get_loc("CD4") if "CD4" in genes else None
    cd8_idx = genes.get_loc("CD8A") if "CD8A" in genes else None

    cd4cells = np.zeros(adata.n_obs, dtype=int)
    cd8cells = np.zeros(adata.n_obs, dtype=int)
    if cd4_idx is not None and cd8_idx is not None:
        cd4_counts = x[:, cd4_idx]
        cd8_counts = x[:, cd8_idx]
        cd4cells[(cd8_counts == 0) & (cd4_counts > 0)] = 1
        cd8cells[(cd8_counts > 0) & (cd4_counts == 0)] = 1
    adata.obs["CD4cells"] = cd4cells
    adata.obs["CD8cells"] = cd8cells


def _annotate_trb_clonotypes(
    adata_list: list[ad.AnnData],
    vdj_tables: list[pd.DataFrame],
    productive_value: str = "true",
) -> None:
    """Assign clonotype id / clone size / CDR3 nucleotide & amino-acid
    sequence to each cell from its dominant productive TRB chain, and
    consolidate clonotype names/sizes across all samples.

    Direct port of the "Annotate cells with clonotypes (TCRB)" section at
    the end of ``CombineData`` in ``Import.R``. Mutates each ``AnnData`` in
    ``adata_list`` in place, adding ``obs['cdr3_na']``, ``obs['cdr3']``,
    ``obs['clonotype']``, ``obs['countcln']``.
    """
    n_samples = len(adata_list)

    mono_trb_na: list[np.ndarray] = []
    mono_trb_aa: list[np.ndarray] = []
    for k in range(n_samples):
        barcodes = adata_list[k].obs_names
        trb = vdj_tables[k]
        trb = trb[(trb["productive"].astype(str).str.lower() == productive_value.lower()) & (trb["chain"] == "TRB")]
        na = np.full(len(barcodes), "unassigned", dtype=object)
        aa = np.full(len(barcodes), "unassigned", dtype=object)
        grouped = {bc: sub for bc, sub in trb.groupby("barcode")}
        for j, bc in enumerate(barcodes):
            sub = grouped.get(bc)
            if sub is None or len(sub) == 0:
                continue
            if len(sub) == 1:
                na[j] = sub["cdr3_nt"].iloc[0]
                aa[j] = sub["cdr3"].iloc[0]
            else:
                pos = sub["reads"].to_numpy().argmax()
                na[j] = sub["cdr3_nt"].iloc[pos]
                aa[j] = sub["cdr3"].iloc[pos]
        mono_trb_na.append(na)
        mono_trb_aa.append(aa)

    intersect_tcrs = set(mono_trb_na[0]) - {"unassigned"}
    for h in range(1, n_samples):
        intersect_tcrs &= set(mono_trb_na[h]) - {"unassigned"}
    intersect_tcrs_list = sorted(intersect_tcrs)

    imtcrs: list[list[str]] = []
    sizes: list[np.ndarray] = []
    for h in range(n_samples):
        tail = sorted(set(mono_trb_na[h]) - {"unassigned"} - set(intersect_tcrs_list))
        combined = intersect_tcrs_list + tail
        imtcrs.append(combined)
        trb = vdj_tables[h]
        trb = trb[(trb["productive"].astype(str).str.lower() == productive_value.lower()) & (trb["chain"] == "TRB")]
        counts = trb.groupby("cdr3_nt")["barcode"].nunique()
        sizes.append(np.array([counts.get(t, 0) for t in combined]))

    nclono: list[dict[str, str]] = []
    for q in range(n_samples):
        order = np.argsort(-sizes[q], kind="stable")
        names = [""] * len(imtcrs[q])
        for rank, f in enumerate(order):
            if f >= len(intersect_tcrs_list):
                names[f] = f"clonotype{rank + 1}_{q + 1}"
            else:
                names[f] = f"clonotype{rank + 1}"
        nclono.append(dict(zip(imtcrs[q], names)))

    for k in range(n_samples):
        trb = vdj_tables[k]
        trb = trb[(trb["productive"].astype(str).str.lower() == productive_value.lower()) & (trb["chain"] == "TRB")]
        freq_by_nt = trb.groupby("cdr3_nt")["barcode"].nunique().to_dict()

        clonotype = np.full(adata_list[k].n_obs, "unassigned", dtype=object)
        countcln = np.zeros(adata_list[k].n_obs, dtype=int)
        for j, nt in enumerate(mono_trb_na[k]):
            if nt == "unassigned":
                continue
            countcln[j] = freq_by_nt.get(nt, 0)
            clonotype[j] = nclono[k].get(nt, "unassigned")

        adata_list[k].obs["cdr3_na"] = mono_trb_na[k]
        adata_list[k].obs["cdr3"] = mono_trb_aa[k]
        adata_list[k].obs["clonotype"] = clonotype
        adata_list[k].obs["countcln"] = countcln


def combine_data(
    gex_path: list[str | Path],
    marker_gene: list[str],
    vdj_path: list[str | Path] | None = None,
    sample_name: list[str] | None = None,
    n_feature_lower_q: float = 0.05,
    n_feature_upper_q: float = 0.95,
    percent_mt_upper_q: float = 0.95,
    demultiplex: bool = False,
    demultiplex_index: list[int] | None = None,
    nameshashtags: list[str] | None = None,
    n_ht_per_sample: int = 1,
    file_saved: str | Path | None = None,
    verbose: bool = True,
) -> list[ad.AnnData]:
    """Read 10x Genomics feature-barcode matrices (+ optional VDJ contig
    annotation CSVs), QC filter, and annotate with marker-gene / CD4-CD8 /
    clonotype metadata.

    Port of ``CombineData`` (``Import.R``). Returns a list of ``AnnData``
    objects, the AnnData analogue of the list-of-Seurat-objects returned by
    the R function.

    Parameters
    ----------
    gex_path:
        One 10x ``filtered_feature_bc_matrix`` directory (or ``.h5`` file)
        path per gene-expression sample.
    marker_gene:
        Gene(s) of interest to check are present in the panel (kept as
        metadata for downstream thresholding -- see
        :mod:`paocseq.classify`).
    vdj_path:
        One ``all_contig_annotations.csv`` path per sample (same order as
        ``gex_path``), or ``None`` to skip VDJ/clonotype annotation.
    sample_name:
        Optional explicit sample names; inferred from the trailing path
        component of each ``gex_path`` entry otherwise.
    n_feature_lower_q, n_feature_upper_q, percent_mt_upper_q:
        Quantile QC cutoffs, matching the ``lowerQ``/``upperQ`` defaults in
        the R function.
    demultiplex:
        If True, treat the antibody-capture panel of each ``gex_path`` entry
        as containing hashtag oligos and demultiplex with
        :func:`paocseq.demux.gmm_demux` before splitting each sample into
        one ``AnnData`` per hashtag.
    demultiplex_index:
        0-based feature indices of the hashtags within the antibody-capture
        panel (only used if ``demultiplex=True``).
    nameshashtags:
        Sample names for each hashtag (only used if ``demultiplex=True``).
    n_ht_per_sample:
        Number of hashtags multiplexed per ``gex_path`` entry.
    file_saved:
        If given, the returned list is also pickled to this path (the
        AnnData analogue of R's ``saveRDS(Dataset, file.saved)``).

    Returns
    -------
    list[anndata.AnnData]
    """
    n_gex = len(gex_path)
    if sample_name is None:
        sample_name = [_infer_sample_name(p) for p in gex_path]

    raw = [read_10x_sample(p, sample_name=sample_name[k]) for k, p in enumerate(gex_path)]

    if demultiplex:
        if demultiplex_index is None or nameshashtags is None:
            raise ValueError("demultiplex=True requires demultiplex_index and nameshashtags")
        adata_list: list[ad.AnnData] = []
        for q in range(n_gex):
            adata = raw[q]
            if "Antibody Capture" not in adata.var.get("feature_types", pd.Series(dtype=str)).unique().tolist():
                raise ValueError(
                    f"gex_path[{q}] does not have an 'Antibody Capture' feature type; "
                    "cannot demultiplex hashtags."
                )
            ab_mask = adata.var["feature_types"] == "Antibody Capture"
            ab = adata[:, ab_mask]
            ab_counts = pd.DataFrame(
                ab.X.toarray() if hasattr(ab.X, "toarray") else np.asarray(ab.X),
                index=ab.obs_names,
                columns=ab.var_names,
            )
            start = q * n_ht_per_sample
            these_idx = demultiplex_index[start : start + n_ht_per_sample]
            these_names = nameshashtags[start : start + n_ht_per_sample]
            labels = _demux.gmm_demux(ab_counts, these_idx, these_names)

            gex_only = adata[:, ~ab_mask].copy()
            gex_only.obs["hashtag_label"] = labels.reindex(gex_only.obs_names).values
            for name in these_names:
                sub = gex_only[gex_only.obs["hashtag_label"] == name].copy()
                sub.obs["sampleref"] = name
                adata_list.append(sub)
        n_samples = len(adata_list)
    else:
        adata_list = raw
        n_samples = n_gex

    filtered: list[ad.AnnData] = []
    for k in range(n_samples):
        adata = adata_list[k]
        if verbose:
            print(f"QC filtering sample {k + 1}/{n_samples} ({adata.n_obs} cells)...")
        adata = _qc_filter(adata, n_feature_lower_q, n_feature_upper_q, percent_mt_upper_q)
        _add_cd4_cd8_indicator(adata)
        missing = [g for g in marker_gene if g not in adata.var_names]
        if missing and verbose:
            print(f"  Warning: marker genes not found in panel: {missing}")
        adata.uns["marker_gene"] = list(marker_gene)
        filtered.append(adata)

    if vdj_path:
        if demultiplex:
            expanded_vdj: list[str | Path] = []
            for q in range(n_gex):
                expanded_vdj.extend([vdj_path[q]] * n_ht_per_sample)
            vdj_path_use = expanded_vdj
        else:
            vdj_path_use = list(vdj_path)

        vdj_tables = []
        for k in range(n_samples):
            df = pd.read_csv(vdj_path_use[k])
            cell_mask = df["is_cell"].astype(str).str.lower() == "true"
            productive_mask = df["productive"].astype(str).str.lower() == "true"
            df = df[cell_mask & productive_mask]
            barcode_set = set(filtered[k].obs_names)
            testintersection = list(set(df["barcode"].tolist()) & set(list(barcode_set)))
            df = df[df["barcode"].isin(barcode_set)]
            keep_barcodes = sorted(barcode_set & set(df["barcode"]))
            filtered[k] = filtered[k][filtered[k].obs_names.isin(keep_barcodes)].copy()
            vdj_tables.append(df)

        # Re-subset vdj tables now that filtered[k] has been trimmed to the
        # shared barcode set.
        for k in range(n_samples):
            vdj_tables[k] = vdj_tables[k][vdj_tables[k]["barcode"].isin(set(filtered[k].obs_names))]

        _annotate_trb_clonotypes(filtered, vdj_tables)

    if file_saved is not None:
        import pickle

        with open(file_saved, "wb") as fh:
            pickle.dump(filtered, fh)

    return filtered
