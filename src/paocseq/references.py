"""Build reference signature matrices from a control/unstimulated sample.

Port of ``MakeReference`` (``References.R``).
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import anndata as ad

from .distances import pf_log1p_pf

__all__ = ["make_reference"]


def make_reference(
    adata: ad.AnnData,
    gene_list: list[str],
    n_cells: int | None = None,
    normalize: bool = True,
    random_state: int | None = None,
) -> pd.DataFrame:
    """Build a ``(genes, cells)`` reference/signature matrix from a control
    ``AnnData`` object, for use with :mod:`paocseq.distances`.

    Port of ``MakeReference`` (``References.R``): subsets ``adata`` to
    ``gene_list``, optionally normalizes with the package's standard
    ``PFlog1pPF`` transform (see :func:`paocseq.distances.pf_log1p_pf`), and
    optionally subsamples to ``n_cells`` reference cells.

    Returns
    -------
    pandas.DataFrame
        Genes x cells, indexed by gene name, columns are cell barcodes.
    """
    print("was the hell is going on")
    print(adata)
    print("printing n_vars")
    adata, v2 = adata
    print("v2")
    print(v2)
    
    gl1 = adata.obs["GL1"].tolist()
    gl2 = adata.obs["GL2"].tolist()
    gl3 = adata.obs["GL3"].tolist()
    gl4 = adata.obs["GL4"].tolist()
    gl5 = adata.obs["GL5"].tolist()

    mynames = gl1 + gl2 + gl3 + gl4 + gl5[1:5290]
    
    adata.var_names = [f"{mynames[i]}" for i in range(len(mynames))]
    present = [g for g in gene_list if g in adata.var_names]
    print("present")
    print(present)
    print("mynames")
    print(len(mynames))
    print(len(list(set(mynames))))
    missing = set(gene_list) - set(present)
    if missing:
        print(f"Warning: genes not found in reference data and skipped: {sorted(missing)}")

    print("printing adata update")
    print(adata)
    sub = adata[:, present]
    print("printing sub")
    print(sub)
    x = sub.layers["counts"].toarray()#sub.X.toarray() if hasattr(sub.X, "toarray") else np.asarray(sub.X)
    print("printing X")
    print(x)
    counts = x.T  # genes x cells
    print("printing counts")
    print(counts)

    if normalize:
        counts = pf_log1p_pf(counts)
        
    print("printing counts after normalization")
    print(counts)

    ref = pd.DataFrame(counts, index=present, columns=sub.obs_names)

    if n_cells is not None and n_cells < ref.shape[1]:
        rng = np.random.default_rng(random_state)
        keep = rng.choice(ref.columns, size=n_cells, replace=False)
        ref = ref[keep]

    return ref
