"""Two-component Gaussian-mixture EM algorithm and its two aocseq use cases:

* :func:`gmm_demux` -- hashtag-oligo (HTO) demultiplexing ("GMMDemux" in R).
* :func:`phenotype_markers` -- CD4 / CD8 phenotype calling from a 1-D PCA
  embedding of CD4/CD8A/CD8B counts ("PhenotypeMarkers" in R).

Both are ports of the corresponding functions in ``Import.R``. The
underlying two-component EM routine (``EStep`` / ``MStep`` in R) is exposed
here as :func:`e_step` / :func:`m_step` for reuse/testing.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd
from scipy.stats import norm

__all__ = [
    "EMResult",
    "e_step",
    "m_step",
    "fit_two_component_em",
    "gmm_demux",
    "phenotype_markers",
]


@dataclass
class EMResult:
    mu: np.ndarray  # (2,)
    var: np.ndarray  # (2,)
    alpha: np.ndarray  # (2,)
    loglik_trace: list[float]


def e_step(x: np.ndarray, mu: np.ndarray, sd: np.ndarray, alpha: np.ndarray):
    """Expectation step of a two-component 1-D Gaussian mixture EM.

    Direct port of ``EStep`` in ``Import.R``.

    Returns
    -------
    loglik: float
        Sum of log(mixture density) over non-NaN entries of ``x``.
    posterior: np.ndarray
        ``(n, 2)`` posterior responsibilities for each component.
    """
    x = np.asarray(x, dtype=np.float64)
    comp1 = norm.pdf(x, mu[0], sd[0]) * alpha[0]
    comp2 = norm.pdf(x, mu[1], sd[1]) * alpha[1]
    total = comp1 + comp2
    with np.errstate(divide="ignore", invalid="ignore"):
        post1 = comp1 / total
        post2 = comp2 / total
        loglik_terms = np.log(total)
    loglik = np.nansum(loglik_terms)
    posterior = np.column_stack([post1, post2])
    return loglik, posterior


def m_step(x: np.ndarray, posterior: np.ndarray):
    """Maximization step. Direct port of ``MStep`` in ``Import.R``."""
    x = np.asarray(x, dtype=np.float64)
    n1 = np.nansum(posterior[:, 0])
    n2 = np.nansum(posterior[:, 1])

    mu1 = np.nansum(posterior[:, 0] * x) / n1
    mu2 = np.nansum(posterior[:, 1] * x) / n2

    var1 = np.nansum(posterior[:, 0] * (x - mu1) ** 2) / n1
    var2 = np.nansum(posterior[:, 1] * (x - mu2) ** 2) / n2

    alpha1 = n1 / len(x)
    alpha2 = n2 / len(x)

    return {
        "mu": np.array([mu1, mu2]),
        "var": np.array([var1, var2]),
        "alpha": np.array([alpha1, alpha2]),
    }


def fit_two_component_em(
    x: np.ndarray,
    max_iter: int = 50,
    tol: float = 1e-6,
    random_state: int | None = None,
) -> EMResult:
    """Fit a two-component 1-D Gaussian mixture by EM, kmeans-initialized.

    Port of the fit loop embedded in ``GMMDemux`` / ``PhenotypeMarkers``:
    the two components are first seeded with a 1-D 2-means clustering (done
    here with a tiny closed-form Lloyd's algorithm, since scikit-learn is
    intentionally not a dependency), then refined with EM until the
    log-likelihood changes by less than ``tol`` or ``max_iter`` is reached.
    """
    x = np.asarray(x, dtype=np.float64)
    valid = x[~np.isnan(x)]
    rng = np.random.default_rng(random_state)

    # --- 1-D 2-means initialization (mirrors R's kmeans(na.omit(wait), 2)) ---
    centers = np.array(
        [valid.min() + 0.25 * (valid.max() - valid.min()), valid.min() + 0.75 * (valid.max() - valid.min())]
    )
    if centers[0] == centers[1]:
        centers = np.array([valid.mean() - 1.0, valid.mean() + 1.0])
    for _ in range(100):
        d0 = np.abs(valid - centers[0])
        d1 = np.abs(valid - centers[1])
        cluster = (d1 < d0).astype(int)
        new_centers = np.array(
            [
                valid[cluster == 0].mean() if np.any(cluster == 0) else centers[0],
                valid[cluster == 1].mean() if np.any(cluster == 1) else centers[1],
            ]
        )
        if np.allclose(new_centers, centers):
            centers = new_centers
            break
        centers = new_centers

    mu0 = np.array(
        [
            valid[cluster == 0].mean() if np.any(cluster == 0) else valid.mean(),
            valid[cluster == 1].mean() if np.any(cluster == 1) else valid.mean(),
        ]
    )
    var0 = np.array(
        [
            valid[cluster == 0].var(ddof=1) if np.sum(cluster == 0) > 1 else valid.var(ddof=1),
            valid[cluster == 1].var(ddof=1) if np.sum(cluster == 1) > 1 else valid.var(ddof=1),
        ]
    )
    sd0 = np.sqrt(var0)
    size = np.array([np.sum(cluster == 0), np.sum(cluster == 1)], dtype=np.float64)
    alpha0 = size / size.sum()

    loglik_trace: list[float] = []
    mu, var, alpha = mu0, var0, alpha0
    cur_loglik = None
    for i in range(max_iter):
        sd = np.sqrt(var)
        loglik, posterior = e_step(x, mu, sd, alpha)
        m = m_step(x, posterior)
        mu, var, alpha = m["mu"], m["var"], m["alpha"]
        loglik_trace.append(loglik)
        if cur_loglik is not None and abs(cur_loglik - loglik) < tol:
            break
        cur_loglik = loglik

    return EMResult(mu=mu, var=var, alpha=alpha, loglik_trace=loglik_trace)


def _clr_transform(counts_plus_one: np.ndarray) -> np.ndarray:
    """Centered log-ratio transform of a (cells x hashtags) count matrix
    (each entry already has +1 added), matching the ``clrfn`` closure in
    ``GMMDemux``: ``log(x / geometric_mean(x[x>0]))`` per hashtag column."""
    out = np.full_like(counts_plus_one, np.nan, dtype=np.float64)
    for j in range(counts_plus_one.shape[1]):
        col = counts_plus_one[:, j]
        positive = col[col > 0]
        if positive.size == 0:
            continue
        denom = np.exp(np.mean(np.log(positive)))
        with np.errstate(divide="ignore", invalid="ignore"):
            out[:, j] = np.log(col / denom)
    out[~np.isfinite(out)] = np.nan
    return out


def gmm_demux(
    hto_counts: pd.DataFrame,
    hashtag_index: list[int],
    names_hashtags: list[str],
    random_state: int | None = None,
) -> pd.Series:
    """Assign each cell to the hashtag (sample) whose antibody-capture
    counts most strongly identify it, via a per-hashtag two-component GMM
    fit on centered-log-ratio-transformed counts.

    Port of ``GMMDemux`` (``Import.R``), simplified to only handle the
    labelling step: given an already-loaded hashtag/antibody-capture count
    matrix (cells x hashtags), it returns a per-cell label. Use
    :func:`paocseq.io.combine_data` with ``demultiplex=True`` to run this
    directly against a 10x feature-barcode matrix.

    Parameters
    ----------
    hto_counts:
        ``(n_cells, n_hashtags)`` DataFrame of raw HTO/antibody counts,
        indexed by cell barcode.
    hashtag_index:
        0-based column positions (within ``hto_counts``) of the hashtags to
        demultiplex on.
    names_hashtags:
        Sample name for each entry of ``hashtag_index``, in the same order.
    random_state:
        Optional seed for the EM initialization.

    Returns
    -------
    pd.Series
        Per-cell label (one of ``names_hashtags`` or ``"none"``/``"unassigned"``),
        indexed like ``hto_counts``.
    """
    if len(hashtag_index) != len(names_hashtags):
        raise ValueError("hashtag_index and names_hashtags must be the same length")

    counts_plus_one = hto_counts.to_numpy(dtype=np.float64) + 1.0
    clr = _clr_transform(counts_plus_one)

    n_cells = clr.shape[0]
    fits: list[EMResult] = []
    for col in hashtag_index:
        fit = fit_two_component_em(clr[:, col], random_state=random_state)
        fits.append(fit)

    donor_label = np.full(n_cells, "unassigned", dtype=object)

    for j in range(n_cells):
        high_probs = []
        low_probs = []
        for k_pos, col in enumerate(hashtag_index):
            fit = fits[k_pos]
            mu, var, alpha = fit.mu, fit.var, fit.alpha
            x = clr[j, col]
            if np.isnan(x):
                high_probs.append(0.0)
                low_probs.append(1.0)
                continue
            high_i, low_i = (1, 0) if mu[1] > mu[0] else (0, 1)
            sd = np.sqrt(var)
            p_high = norm.pdf(x, mu[high_i], sd[high_i])
            p_low = norm.pdf(x, mu[low_i], sd[low_i])
            px = p_low * alpha[low_i] + p_high * alpha[high_i]
            if px == 0 or np.isnan(px):
                high_probs.append(0.0)
                low_probs.append(1.0)
                continue
            high_probs.append((p_high * alpha[high_i]) / px)
            low_probs.append((p_low * alpha[low_i]) / px)

        n_ht = len(hashtag_index)
        probs = np.zeros(n_ht + 1)
        for i in range(n_ht):
            others_low = [low_probs[m] for m in range(n_ht) if m != i]
            probs[i] = high_probs[i] * float(np.prod(others_low)) if others_low else high_probs[i]
        probs[n_ht] = 1.0 - probs[:n_ht].sum()

        best = int(np.argmax(probs))
        if best < n_ht:
            donor_label[j] = names_hashtags[best]
        else:
            donor_label[j] = "none"

    return pd.Series(donor_label, index=hto_counts.index, name="hashtag_label")


def phenotype_markers(
    pc1: np.ndarray,
    cd4: np.ndarray,
    cd8a: np.ndarray,
    cd8b: np.ndarray,
    random_state: int | None = None,
) -> np.ndarray:
    """Assign each cell a CD4 / CD8 phenotype from a 1-D PCA embedding of
    CD4/CD8A/CD8B counts.

    Port of ``PhenotypeMarkers`` (``Import.R``): fits a two-component EM on
    ``pc1``, uses the fitted mixture to assign each cell that expresses at
    least one of CD4/CD8A/CD8B to the higher- or lower-mean component, and
    then resolves which component corresponds to CD8 (the component with
    both the higher max and higher total CD8A count is called CD8; the
    other CD4).

    Parameters
    ----------
    pc1:
        ``(n_cells,)`` first principal component of [CD4, CD8A, CD8B] counts.
    cd4, cd8a, cd8b:
        ``(n_cells,)`` raw counts for the three genes.

    Returns
    -------
    np.ndarray
        ``(n_cells,)`` array of ``"CD4"``, ``"CD8"``, or ``"unassigned"``.
    """
    pc1 = np.asarray(pc1, dtype=np.float64)
    cd4 = np.asarray(cd4, dtype=np.float64)
    cd8a = np.asarray(cd8a, dtype=np.float64)
    cd8b = np.asarray(cd8b, dtype=np.float64)

    fit = fit_two_component_em(pc1, max_iter=25, random_state=random_state)
    mu, var, alpha = fit.mu, fit.var, fit.alpha
    sd = np.sqrt(var)
    high_i, low_i = (1, 0) if mu[1] > mu[0] else (0, 1)

    phenotype = np.full(pc1.shape[0], "unassigned", dtype=object)
    expresses_any = (cd4 > 0) | (cd8a > 0) | (cd8b > 0)

    p_high = norm.pdf(pc1, mu[high_i], sd[high_i])
    p_low = norm.pdf(pc1, mu[low_i], sd[low_i])
    px = p_low * alpha[low_i] + p_high * alpha[high_i]
    with np.errstate(divide="ignore", invalid="ignore"):
        post_high = np.where(px > 0, (p_high * alpha[high_i]) / px, np.nan)
        post_low = np.where(px > 0, (p_low * alpha[low_i]) / px, np.nan)

    class_a = expresses_any & (post_high > post_low)
    class_b = expresses_any & (post_high < post_low)
    phenotype[class_a] = "ClassA"
    phenotype[class_b] = "ClassB"

    max_a_cd8 = cd8a[class_a].max() if class_a.any() else -np.inf
    max_b_cd8 = cd8a[class_b].max() if class_b.any() else -np.inf
    sum_a_cd8 = cd8a[class_a].sum() if class_a.any() else -np.inf
    sum_b_cd8 = cd8a[class_b].sum() if class_b.any() else -np.inf

    out = phenotype.copy()
    if max_a_cd8 > max_b_cd8 and sum_a_cd8 > sum_b_cd8:
        out[class_a] = "CD8"
        out[class_b] = "CD4"
    elif max_a_cd8 < max_b_cd8 and sum_a_cd8 < sum_b_cd8:
        out[class_a] = "CD4"
        out[class_b] = "CD8"
    else:
        out[class_a] = "ClassA"
        out[class_b] = "ClassB"

    return out
