"""Pure-Python/numpy fallback for the numeric solvers ported from the
original R package ``aocseq``'s Rcpp extensions.

This module is used automatically by :mod:`paocseq._solvers` when the
compiled C++ extension (``paocseq._solvers_ext``, built from
``csrc/mahalanobis.cpp`` + ``csrc/isolation_forest.cpp`` -- direct,
dependency-free translations of the originals) is not available, e.g.
because no C++ compiler was present at install time. It is kept as a
readable, dependency-light reference implementation and as a safety net;
prefer the compiled extension when possible, since it matches the original
algorithms' control flow (and performance characteristics) most closely.

Two numerical routines are covered here:

1. ``GetMahalanobis`` (``csrc/mahalanobis.cpp``, ported from the corrected
   ``GetMahalanobis.cpp``) -- the mean vector and gene-by-gene covariance
   matrix are computed from the reference/signature panel only (the query
   cell being scored is *not* folded into those reference statistics), the
   covariance matrix is inverted (here with ``np.linalg.inv`` in place of
   the hand-written Gauss-Jordan elimination), and the Mahalanobis distance
   of each query cell to that fixed reference distribution is returned.
2. ``isoForest`` / ``tree`` (``csrc/isolation_forest.cpp``, ported from the
   legacy ``solvers.cpp``) -- a bespoke isolation-forest variant. Instead of
   choosing a random feature to split on (as in a standard isolation
   forest), it always splits on the *gene with the largest excess kurtosis*
   in the reference+query panel, then recurses on randomly-thresholded
   halves of the (reference + one query cell) samples, returning a
   normalized outlier score for the query cell (1.0 == deeply isolated /
   likely outlier).

Both routines are re-implemented here as vectorized numpy array
operations, keeping the same statistical algorithm and inputs/outputs as
the compiled extension, but without the C++ dependency.
"""

from __future__ import annotations

import numpy as np

__all__ = ["mahalanobis_distance", "isolation_forest_score"]


def mahalanobis_distance(
    signature_cells: np.ndarray,
    test_cells: np.ndarray,
) -> np.ndarray:
    """Mahalanobis distance of each query cell to a fixed reference panel.

    Mirrors ``GetMahalanobis`` (``csrc/mahalanobis.cpp``): the mean vector
    and gene-by-gene covariance matrix are computed once from
    ``signature_cells`` only, then every query cell in ``test_cells`` is
    scored against that fixed reference distribution (the query cell does
    *not* contribute to the reference mean/covariance).

    Parameters
    ----------
    signature_cells:
        ``(n_genes, n_reference_cells)`` array. Reference / signature panel,
        genes on rows, cells on columns (matches the R ``reference.data``
        layout: genes x cells).
    test_cells:
        ``(n_genes, n_query_cells)`` array. Query panel to score, same gene
        order as ``signature_cells``.

    Returns
    -------
    np.ndarray
        ``(n_query_cells,)`` array of Mahalanobis distances.
    """
    signature_cells = np.asarray(signature_cells, dtype=np.float64)
    test_cells = np.asarray(test_cells, dtype=np.float64)
    if signature_cells.ndim == 1:
        signature_cells = signature_cells[None, :]
    if test_cells.ndim == 1:
        test_cells = test_cells[None, :]

    n_genes, n_ref = signature_cells.shape
    n_query = test_cells.shape[1]

    # Reference mean/covariance/inverse only depend on signature_cells, so
    # (unlike the legacy leave-one-out version) this is computed once and
    # reused for every query cell.
    means = signature_cells.mean(axis=1)
    centered = signature_cells - means[:, None]
    cov = (centered @ centered.T) / n_ref
    try:
        inv_cov = np.linalg.inv(cov)
    except np.linalg.LinAlgError:
        inv_cov = np.linalg.pinv(cov)

    distances = np.empty(n_query, dtype=np.float64)
    for z in range(n_query):
        diff = test_cells[:, z] - means
        distances[z] = np.sqrt(max(float(diff @ inv_cov @ diff), 0.0))

    return distances


def _excess_kurtosis(matrix: np.ndarray) -> np.ndarray:
    """Per-row excess-kurtosis-like statistic, matching the C++ formula in isolation_forest.cpp.

    The C++ code computes, per gene ``g``::

        fourth_moment = (mean(x - mean(x)))^4          # note: NOT E[(x-mean)^4]
        standard_dev  = (mean((x - mean(x))^2))^2
        kurtosis[g]   = fourth_moment / standard_dev

    This is reproduced verbatim (rather than substituted with the textbook
    kurtosis) so the gene chosen to split on matches the original algorithm's
    behaviour as closely as possible.
    """
    n = matrix.shape[1]
    mean = matrix.mean(axis=1, keepdims=True)
    centered = matrix - mean
    fourth_moment_fill = centered.mean(axis=1)
    standard_dev_fill = (centered**2).mean(axis=1)
    fourth_moment = fourth_moment_fill**4
    standard_dev = standard_dev_fill**2
    with np.errstate(divide="ignore", invalid="ignore"):
        kurt = np.where(standard_dev > 0, fourth_moment / standard_dev, 0.0)
    return np.nan_to_num(kurt, nan=0.0, posinf=0.0, neginf=0.0)


def _isolation_tree_heights(
    values: np.ndarray,
    max_height: int,
    rng: np.random.Generator,
) -> np.ndarray:
    """Recursively split ``values`` on random thresholds, recording the
    tree-depth ("height") at which each sample becomes isolated.

    This mirrors the recursive ``tree`` function in ``csrc/isolation_forest.cpp``: at each
    node, a split value is drawn uniformly between the min and max of the
    remaining values, samples are partitioned into left/right, and recursion
    continues until a node holds a single point, all points are identical,
    or ``max_height`` is reached.
    """
    heights = np.zeros(values.shape[0], dtype=np.float64)

    def _recurse(idx: np.ndarray, depth: int) -> None:
        if idx.size <= 1 or depth >= max_height:
            heights[idx] = depth
            return
        subset = values[idx]
        lo, hi = subset.min(), subset.max()
        if lo == hi:
            heights[idx] = depth
            return
        # NOTE: the original isolation_forest.cpp uses `mindata + runif(0,1)*maxdata`
        # (not the more usual `mindata + runif(0,1)*(maxdata-mindata)`).
        # Reproduced verbatim here to match the compiled extension / original
        # algorithm's behaviour exactly, quirk and all.
        split_value = lo + rng.uniform(0, 1) * hi
        left_mask = subset <= split_value
        left_idx, right_idx = idx[left_mask], idx[~left_mask]
        if left_idx.size == 0 or right_idx.size == 0:
            heights[idx] = depth
            return
        _recurse(left_idx, depth + 1)
        _recurse(right_idx, depth + 1)

    _recurse(np.arange(values.shape[0]), 0)
    return heights


def _path_length_normalizer(n: int) -> float:
    """c(n) normalization constant used by isolation forests (average path
    length of an unsuccessful search in a BST of ``n`` points), matching the
    formula used in both ``isoForest`` (C++) and ``NormalizationScore`` (R)."""
    if n <= 1:
        return 1.0
    return 2.0 * (np.log(n - 1) + np.euler_gamma) - (2.0 * (np.log(n - 1) / np.log(n)))


def isolation_forest_score(
    reference_cells: np.ndarray,
    query_cells: np.ndarray,
    n_trees: int = 10,
    max_height: int = 30,
    random_state: int | None = None,
) -> np.ndarray:
    """Outlier score for each query cell relative to a reference panel.

    Reimplements ``isoForest``/``tree`` from ``csrc/isolation_forest.cpp``: for every
    query cell, the cell is appended to the reference panel, the gene with
    the largest kurtosis (see :func:`_excess_kurtosis`) is selected, an
    isolation forest of ``n_trees`` trees is built on that one gene's
    values, and the average isolation depth of the query cell (last column)
    is converted to a normalized score in (0, 1], where values close to 1
    indicate the query cell is an outlier relative to the reference panel.

    Parameters
    ----------
    reference_cells:
        ``(n_genes, n_reference_cells)`` reference / signature panel.
    query_cells:
        ``(n_genes, n_query_cells)`` panel of cells to score.
    n_trees:
        Number of isolation trees to average over per query cell.
    max_height:
        Maximum recursion depth per tree.
    random_state:
        Optional seed for reproducibility.

    Returns
    -------
    np.ndarray
        ``(n_query_cells,)`` array of outlier scores in (0, 1].
    """
    reference_cells = np.asarray(reference_cells, dtype=np.float64)
    query_cells = np.asarray(query_cells, dtype=np.float64)
    if reference_cells.ndim == 1:
        reference_cells = reference_cells[None, :]
    if query_cells.ndim == 1:
        query_cells = query_cells[None, :]

    n_genes, n_ref = reference_cells.shape
    n_query = query_cells.shape[1]
    n_panel = n_ref + 1
    rng = np.random.default_rng(random_state)

    scores = np.empty(n_query, dtype=np.float64)
    c_norm = _path_length_normalizer(n_panel)

    for r in range(n_query):
        panel = np.concatenate([reference_cells, query_cells[:, r : r + 1]], axis=1)
        kurt = _excess_kurtosis(panel)
        gene_idx = int(np.argmax(kurt))
        values = panel[gene_idx, :]

        avg_height = 0.0
        for _ in range(n_trees):
            heights = _isolation_tree_heights(values, max_height, rng)
            avg_height += heights[-1] / n_trees

        scores[r] = 2.0 ** (-avg_height / c_norm)

    return scores
