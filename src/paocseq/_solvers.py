"""Numeric solvers ported from the original R package ``aocseq``'s Rcpp
extensions: a Mahalanobis distance (:func:`mahalanobis_distance`) and a
kurtosis-driven isolation forest (:func:`isolation_forest_score`). See
:mod:`paocseq.distances` and :mod:`paocseq.classify` for how these are used.

This module dispatches to whichever implementation is available:

1. **Preferred:** ``paocseq._solvers_ext``, a compiled CPython extension
   built from ``csrc/mahalanobis.cpp`` + ``csrc/isolation_forest.cpp`` --
   direct, dependency-free C++ ports of the original aocseq
   ``GetMahalanobis.cpp`` and (legacy) ``solvers.cpp``'s isolation forest
   (see each file's header comment for exactly what changed going from
   Rcpp to a plain CPython/NumPy C-API extension). Built automatically by
   ``pip install -e .`` / ``uv sync`` if a C++ compiler is available (see
   ``setup.py``).
2. **Fallback:** :mod:`paocseq._solvers_py`, a pure numpy re-implementation
   of the same two algorithms, used automatically if the compiled extension
   isn't available (e.g. no C++ compiler at install time).

Both implementations expose the same two functions with the same
signatures, so the rest of paocseq (and your own code) can just
``from paocseq import _solvers`` / ``import paocseq.solvers`` and call
``_solvers.mahalanobis_distance(...)`` / ``_solvers.isolation_forest_score(...)``
without caring which backend is active. Check :data:`USING_COMPILED_EXTENSION`
if you want to know which one you got.
"""

from __future__ import annotations

import warnings

import numpy as np

__all__ = ["mahalanobis_distance", "isolation_forest_score", "USING_COMPILED_EXTENSION"]

try:
    from . import _solvers_ext as _ext  # compiled C++ extension, see csrc/

    USING_COMPILED_EXTENSION = True
except ImportError:  # pragma: no cover - exercised whenever no C++ toolchain was available at install time
    _ext = None
    USING_COMPILED_EXTENSION = False
    warnings.warn(
        "paocseq's compiled C++ solver extension (paocseq._solvers_ext) is not "
        "available -- falling back to the pure-Python/numpy implementation in "
        "paocseq._solvers_py. Results are numerically equivalent but slower. "
        "To build the extension, make sure a C++ compiler is available and "
        "reinstall with `pip install -e . --force-reinstall --no-deps` (or "
        "`uv sync --reinstall-package paocseq`).",
        RuntimeWarning,
        stacklevel=2,
    )
    from . import _solvers_py as _fallback


def mahalanobis_distance(signature_cells: np.ndarray, test_cells: np.ndarray) -> np.ndarray:
    """Leave-one-out Mahalanobis distance of each query cell to a reference set.

    Direct port of ``GetMahalanobis`` (``csrc/mahalanobis.cpp``, the
    corrected version): the mean vector and gene-by-gene covariance matrix
    are computed from the reference/signature panel only, inverted by
    Gauss-Jordan elimination, and the Mahalanobis distance of each query
    cell to that fixed reference distribution is returned.

    Parameters
    ----------
    signature_cells:
        ``(n_genes, n_reference_cells)`` array. Reference / signature panel,
        genes on rows, cells on columns.
    test_cells:
        ``(n_genes, n_query_cells)`` array. Query panel to score, same gene
        order as ``signature_cells``.

    Returns
    -------
    np.ndarray
        ``(n_query_cells,)`` array of Mahalanobis distances (float64).
    """
    signature_cells = np.asarray(signature_cells, dtype=np.float32)
    test_cells = np.asarray(test_cells, dtype=np.float32)
    if signature_cells.ndim == 1:
        signature_cells = signature_cells[None, :]
    if test_cells.ndim == 1:
        test_cells = test_cells[None, :]

    if USING_COMPILED_EXTENSION:
        print("USES COMPILED VERSION!!!!! YEY!!!!")
        out = _ext.mahalanobis_distance(np.ascontiguousarray(signature_cells), np.ascontiguousarray(test_cells))
        return np.asarray(out, dtype=np.float64)
    return _fallback.mahalanobis_distance(signature_cells, test_cells)


def isolation_forest_score(
    reference_cells: np.ndarray,
    query_cells: np.ndarray,
    n_trees: int = 10,
    max_height: int = 30,
    random_state: int | None = None,
) -> np.ndarray:
    """Outlier score for each query cell relative to a reference panel.

    Direct port of ``isoForest``/``tree`` from ``csrc/isolation_forest.cpp``:
    for every query cell, the cell is appended to the reference panel, the
    gene with the largest kurtosis is selected, an isolation forest of
    ``n_trees`` trees is built on that one gene's values (re-selecting the
    max-kurtosis gene at every node), and the average isolation depth of the
    query cell is converted to a normalized score in (0, 1], where values
    close to 1 indicate the query cell is an outlier relative to the
    reference panel.

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
        ``(n_query_cells,)`` array of outlier scores in (0, 1] (float64).
    """
    reference_cells = np.asarray(reference_cells, dtype=np.float32)
    query_cells = np.asarray(query_cells, dtype=np.float32)
    if reference_cells.ndim == 1:
        reference_cells = reference_cells[None, :]
    if query_cells.ndim == 1:
        query_cells = query_cells[None, :]

    if USING_COMPILED_EXTENSION:
        out = _ext.isolation_forest_score(
            np.ascontiguousarray(reference_cells),
            np.ascontiguousarray(query_cells),
            n_trees=n_trees,
            max_height=max_height,
            random_state=random_state,
        )
        return np.asarray(out, dtype=np.float64)
    return _fallback.isolation_forest_score(
        reference_cells, query_cells, n_trees=n_trees, max_height=max_height, random_state=random_state
    )
