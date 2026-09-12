// -*- C++ -*-
//===-------------------------- mahalanobis.cpp -----------------------------===//
//
// Part of the aocseq / paocseq Project
//
// Direct C++ port of the aocseq `GetMahalanobis.cpp`. The
// algorithm here is unchanged from that file: build the reference panel's
// mean vector and gene-by-gene covariance matrix from the `n_ref` reference
// cells only (the query/test cell is *not* folded into the reference
// statistics, unlike the legacy version), invert the covariance matrix by
// Gauss-Jordan elimination (forward pass to upper-triangular-with-unit-
// diagonal, backward pass to eliminate above the diagonal -- both passes
// use a single scalar elimination factor captured once per row, and take the 
// quadratic form.
//
// Differences from the uploaded GetMahalanobis.cpp, and why:
//   * Rcpp's `std::vector<float>` in/out marshalling to R, and the large
//     number of pre-sized scratch-buffer arguments the R wrapper used to
//     allocate (`dist`, `meansinp`, `colsv`, `covinpt1..4`, `sigtest`,
//     etc.), are gone -- there's no R wrapper here, so this function
//     allocates its own scratch space, and the Python-facing wrapper in
//     py_solvers_module.cpp does the numpy marshalling instead.
//   * The reference mean/covariance/inverse only depend on the reference
//     panel, not on which query cell is being scored (this was already true
//     in the uploaded file -- notice its mean/covariance loops only ever
//     touch `specificCells`, i.e. columns `0..dimN-1` of `SignatureCellsTest`,
//     never the appended query column at index `dimN`). The uploaded R/Rcpp
//     version still recomputed that mean/covariance/inversion from scratch
//     on every iteration of its outer `z` (query-cell) loop, even though the
//     result never changes across `z`. Here that computation is hoisted
//     outside the query loop and done once, then reused for every query
//     cell -- purely a performance optimization, the arithmetic itself is
//     identical to what the original loop computed on each pass.
//   * One addition, offered and agreed on: the pivoting block from the
//     *legacy* solvers.cpp is gone in the corrected GetMahalanobis.cpp (it's
//     present only as a commented-out block), so an exactly-zero pivot
//     (`denom_divide == 0`, e.g. a reference gene with zero variance) would
//     silently divide by zero and propagate NaN/Inf through the rest of the
//     row. A minimal safety net is added here -- if the pivot magnitude is
//     smaller than a small epsilon, it's nudged away from zero by that same
//     epsilon (a tiny ridge term) before dividing, rather than dividing by
//     exactly 0 -- everything else about the elimination (which entries get
//     updated, in what order, with what factor) is unchanged from the
//     uploaded file.
//
//===----------------------------------------------------------------------===//

#include "mahalanobis.hpp"

#include <cmath>
#include <stdexcept>

namespace paocseq {

namespace {
constexpr double kPivotEpsilon = 1e-12;
}

std::vector<float> GetMahalanobis(int n_genes, int n_ref, int n_query,
                                   const std::vector<float> &reference,
                                   const std::vector<float> &query) {
  if (static_cast<int>(reference.size()) != n_genes * n_ref) {
    throw std::invalid_argument("reference size does not match n_genes*n_ref");
  }
  if (static_cast<int>(query.size()) != n_genes * n_query) {
    throw std::invalid_argument("query size does not match n_genes*n_query");
  }

  // --- MeansVec: mean of each gene over the n_ref reference cells only ---
  // (matches the uploaded file's `for(n=0; n<dimN; n++) MeansVec[m] +=
  // SignatureCellsTest[m*(dimN+1)+n]/dimN`, restricted to specificCells --
  // the query cell is not included, unlike the legacy version.)
  std::vector<double> means(n_genes, 0.0);
  for (int g = 0; g < n_genes; ++g) {
    double sum = 0.0;
    for (int c = 0; c < n_ref; ++c) sum += reference[static_cast<size_t>(g) * n_ref + c];
    means[g] = sum / n_ref;
  }

  // --- Covariance_matrix: gene-by-gene covariance over the reference cells,
  // divisor n_ref (matches `covsum/(dimN)` in the uploaded file). ---
  std::vector<double> cov(static_cast<size_t>(n_genes) * n_genes, 0.0);
  for (int i = 0; i < n_genes; ++i) {
    for (int j = 0; j < n_genes; ++j) {
      double covsum = 0.0;
      for (int k = 0; k < n_ref; ++k) {
        covsum += (reference[static_cast<size_t>(i) * n_ref + k] - means[i]) *
                  (reference[static_cast<size_t>(j) * n_ref + k] - means[j]);
      }
      cov[static_cast<size_t>(i) * n_genes + j] = covsum / n_ref;
    }
  }

  // --- Gauss-Jordan elimination, forward pass: make Covariance_matrix
  // upper-triangular with a unit diagonal, tracking the same operations on
  // InvCovariance_matrix (starts as the identity). Same structure as the
  // uploaded file's forward `while(h<dimgene && k<dimgene)` loop. ---
  std::vector<double> inv(static_cast<size_t>(n_genes) * n_genes, 0.0);
  for (int i = 0; i < n_genes; ++i) inv[static_cast<size_t>(i) * n_genes + i] = 1.0;

  int h = 0, k = 0;
  while (h < n_genes && k < n_genes) {
    double denom_divide = cov[static_cast<size_t>(h) * n_genes + k];
    if (std::fabs(denom_divide) < kPivotEpsilon) {
      // Safety net not present in the uploaded file (see header comment):
      // nudge an exactly/near-zero pivot away from zero instead of dividing
      // by it.
      denom_divide = std::copysign(kPivotEpsilon, denom_divide == 0.0 ? 1.0 : denom_divide);
    }

    for (int d = 0; d < n_genes; ++d) {
      inv[static_cast<size_t>(h) * n_genes + d] /= denom_divide;
      if (d == k) {
        cov[static_cast<size_t>(h) * n_genes + d] = 1.0;
      } else {
        cov[static_cast<size_t>(h) * n_genes + d] /= denom_divide;
      }
    }

    for (int i = h + 1; i < n_genes; ++i) {
      double valpre = cov[static_cast<size_t>(i) * n_genes + k];
      for (int d = 0; d < n_genes; ++d) {
        inv[static_cast<size_t>(i) * n_genes + d] -= valpre * inv[static_cast<size_t>(h) * n_genes + d];
      }
      for (int d = 0; d < n_genes; ++d) {
        cov[static_cast<size_t>(i) * n_genes + d] -= valpre * cov[static_cast<size_t>(h) * n_genes + d];
      }
    }

    h += 1;
    k += 1;
  }

  // --- Backward pass: eliminate above the diagonal, same structure as the
  // uploaded file's `while(p<dimgene)` reverse loop. ---
  int p = 0;
  while (p < n_genes) {
    int h1 = n_genes - (p + 1);
    if (h1 > 0) {
      for (int m = 1; m < h1 + 1; ++m) {
        int i = h1 - m;
        double valpretwo = cov[static_cast<size_t>(i) * n_genes + h1];
        for (int s = 0; s < n_genes; ++s) {
          int d = n_genes - (s + 1);
          inv[static_cast<size_t>(i) * n_genes + d] -= valpretwo * inv[static_cast<size_t>(h1) * n_genes + d];
        }
        for (int s = 0; s < n_genes; ++s) {
          int d = n_genes - (s + 1);
          cov[static_cast<size_t>(i) * n_genes + d] -= valpretwo * cov[static_cast<size_t>(h1) * n_genes + d];
        }
      }
    }
    p += 1;
  }

  // --- Per-query-cell quadratic form: distance = sqrt(diff^T * inv * diff),
  // same as the uploaded file's columnVec/distances loops. ---
  std::vector<float> distances(n_query, 0.0f);
  std::vector<double> diff(n_genes);
  for (int z = 0; z < n_query; ++z) {
    for (int g = 0; g < n_genes; ++g) diff[g] = query[static_cast<size_t>(g) * n_query + z] - means[g];

    double dist_sq = 0.0;
    for (int kk = 0; kk < n_genes; ++kk) {
      double column_kk = 0.0;
      for (int i = 0; i < n_genes; ++i) column_kk += diff[i] * inv[static_cast<size_t>(i) * n_genes + kk];
      dist_sq += column_kk * diff[kk];
    }
    distances[z] = static_cast<float>(std::sqrt(std::max(dist_sq, 0.0)));
  }

  return distances;
}

} // namespace paocseq
