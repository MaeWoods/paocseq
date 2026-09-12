// -*- C++ -*-
//===----------------------- isolation_forest.cpp ---------------------------===//
//
// Part of the aocseq / paocseq Project
//
// This is a direct C++ port of the original aocseq package's isolation
// forest (`tree()` / `isoForest()` in the legacy `solvers.cpp`, an Rcpp
// extension). Split into its own file/module (separate from the
// Mahalanobis-distance solver in mahalanobis.cpp) since combining both
// algorithms in one Rcpp source file was causing build issues on the R
// side -- the same split is mirrored here for consistency, even though a
// single CPython extension can still export both.
//
// Differences from the original, and why:
//   * Rcpp's `std::vector<float>` in/out marshalling to R is replaced with
//     plain C++ vectors/arrays; the Python-facing wrapper in
//     py_solvers_module.cpp does the equivalent marshalling to/from numpy.
//   * The original function took a large number of pre-sized scratch
//     buffers as arguments (because the R wrapper allocated them). Since
//     there is no R wrapper here, IsoForest() below allocates and
//     initializes its own scratch space internally; the per-iteration
//     algorithm inside is otherwise unchanged.
//   * `R::runif(0, 1)` is replaced with a seeded `std::mt19937` +
//     `std::uniform_real_distribution<float>`, matching R's default
//     Mersenne Twister generator family (not bit-for-bit identical to R's
//     stream, but the same distribution and generator family).
//
//===----------------------------------------------------------------------===//

#include "isolation_forest.hpp"

#include <algorithm>
#include <cmath>
#include <random>
#include <stdexcept>

namespace paocseq {

// ---------------------------------------------------------------------
// Isolation forest: kurtosis-driven recursive splitting
// ---------------------------------------------------------------------

namespace {

// Recompute, for every gene, the original's non-standard "kurtosis" over
// the currently-active columns only (active[q] == true), then return the
// index of the gene with the largest value. Formula reproduced verbatim
// from solvers.cpp:
//   fourth_moment_fill = mean(x - mean(x))         [note: not squared/etc. before the mean]
//   standard_dev_fill  = mean((x - mean(x))^2)
//   kurtosis           = fourth_moment_fill^4 / standard_dev_fill^2
int SelectSplitGene(const std::vector<float> &panel, const std::vector<char> &active, int n_genes, int panel_n) {
  int best_gene = 0;
  double best_kurtosis = -1.0;
  bool have_best = false;

  for (int g = 0; g < n_genes; ++g) {
    double mean = 0.0;
    int count = 0;
    for (int q = 0; q < panel_n; ++q) {
      if (active[q]) {
        mean += panel[static_cast<size_t>(g) * panel_n + q];
        ++count;
      }
    }
    if (count == 0) continue;
    mean /= count;

    double fourth_moment_fill = 0.0, standard_dev_fill = 0.0;
    for (int q = 0; q < panel_n; ++q) {
      if (active[q]) {
        double d = panel[static_cast<size_t>(g) * panel_n + q] - mean;
        fourth_moment_fill += d;
        standard_dev_fill += d * d;
      }
    }
    fourth_moment_fill /= count;
    standard_dev_fill /= count;

    double fourth_moment = fourth_moment_fill * fourth_moment_fill * fourth_moment_fill * fourth_moment_fill;
    double standard_dev = standard_dev_fill * standard_dev_fill;
    double kurtosis = (standard_dev > 0.0) ? (fourth_moment / standard_dev) : 0.0;

    if (!have_best || kurtosis > best_kurtosis) {
      best_kurtosis = kurtosis;
      best_gene = g;
      have_best = true;
    }
  }
  return best_gene;
}

// Recursive isolation-tree split, mirroring `tree()` in solvers.cpp:
// re-select the max-kurtosis gene among the still-active columns at every
// node, split on a random threshold, recurse, and record the depth at
// which each column stops being split (a leaf, or max_height reached).
void IsoTreeRecurse(const std::vector<float> &panel, std::vector<char> active, std::vector<float> &heights,
                     int panel_n, int n_genes, int dim_split, int cur_height, int max_height,
                     std::mt19937 &rng) {
  if (dim_split <= 1 || cur_height >= max_height) {
    for (int q = 0; q < panel_n; ++q) {
      if (active[q]) heights[q] = static_cast<float>(cur_height);
    }
    return;
  }

  const int gene = SelectSplitGene(panel, active, n_genes, panel_n);

  float mindata = 0.0f, maxdata = 0.0f;
  bool first = true;
  for (int q = 0; q < panel_n; ++q) {
    if (active[q]) {
      float v = panel[static_cast<size_t>(gene) * panel_n + q];
      if (first) {
        mindata = maxdata = v;
        first = false;
      } else {
        if (v < mindata) mindata = v;
        if (v > maxdata) maxdata = v;
      }
    }
  }

  std::vector<char> left(panel_n, 0), right(panel_n, 0);
  int counter_left = 0, counter_right = 0;

  if (mindata == maxdata) {
    // Degenerate split (all active values identical): mirror the original's
    // fallback of dividing the active set roughly in half by encounter order.
    int split_left_count = (dim_split + 1) / 2; // ceil(dim_split / 2)
    int counter = 0;
    for (int q = 0; q < panel_n; ++q) {
      if (active[q]) {
        if (counter < split_left_count) {
          left[q] = 1;
          ++counter_left;
        } else {
          right[q] = 1;
          ++counter_right;
        }
        ++counter;
      }
    }
  } else {
    std::uniform_real_distribution<float> unif(0.0f, 1.0f);
    // NOTE: preserved verbatim from the original solvers.cpp:
    // `split_value = mindata + runif(0,1) * maxdata` (not
    // `mindata + runif(0,1) * (maxdata - mindata)`).
    float split_value = mindata + unif(rng) * maxdata;
    for (int q = 0; q < panel_n; ++q) {
      if (active[q]) {
        float v = panel[static_cast<size_t>(gene) * panel_n + q];
        if (v <= split_value) {
          left[q] = 1;
          ++counter_left;
        } else {
          right[q] = 1;
          ++counter_right;
        }
      }
    }
  }

  if (counter_left > 0) {
    IsoTreeRecurse(panel, left, heights, panel_n, n_genes, counter_left, cur_height + 1, max_height, rng);
  }
  if (counter_right > 0) {
    IsoTreeRecurse(panel, right, heights, panel_n, n_genes, counter_right, cur_height + 1, max_height, rng);
  }
}

} // namespace

// reference: n_genes x n_ref, gene-major. query: n_genes x n_query, gene-major.
// For every query cell r: append it as the last column of the reference
// panel, build `num_trees` isolation trees (re-selecting the max-kurtosis
// gene at every node), average the isolation depth of the query cell across
// trees, and convert to a normalized outlier score with the standard
// isolation-forest path-length normalization constant c(n), matching the
// original's `H_O[r] = 2^(-avg_height / c)` formula.
std::vector<float> IsoForest(int num_trees, int n_genes, int n_ref, int n_query,
                              const std::vector<float> &reference, const std::vector<float> &query,
                              int max_height, unsigned int seed) {
  if (static_cast<int>(reference.size()) != n_genes * n_ref) {
    throw std::invalid_argument("reference size does not match n_genes*n_ref");
  }
  if (static_cast<int>(query.size()) != n_genes * n_query) {
    throw std::invalid_argument("query size does not match n_genes*n_query");
  }

  const int panel_n = n_ref + 1; // "dimN" in the original
  std::vector<float> H_O(n_query, 0.0f);
  std::mt19937 rng(seed);

  std::vector<float> panel(static_cast<size_t>(n_genes) * panel_n);
  for (int g = 0; g < n_genes; ++g) {
    for (int c = 0; c < n_ref; ++c) {
      panel[static_cast<size_t>(g) * panel_n + c] = reference[static_cast<size_t>(g) * n_ref + c];
    }
  }

  const double c_norm =
      (panel_n > 2)
          ? (2.0 * (std::log(panel_n - 1) + 0.5772156649015329) - (2.0 * (std::log(panel_n - 1) / std::log(static_cast<double>(panel_n)))))
          : 1.0;

  for (int r = 0; r < n_query; ++r) {
    for (int g = 0; g < n_genes; ++g) {
      panel[static_cast<size_t>(g) * panel_n + (panel_n - 1)] = query[static_cast<size_t>(g) * n_query + r];
    }

    std::vector<float> avg_height(panel_n, 0.0f);
    for (int t = 0; t < num_trees; ++t) {
      std::vector<char> active(panel_n, 1);
      std::vector<float> heights(panel_n, -1.0f);
      IsoTreeRecurse(panel, active, heights, panel_n, n_genes, panel_n, 0, max_height, rng);
      for (int q = 0; q < panel_n; ++q) avg_height[q] += heights[q] / static_cast<float>(num_trees);
    }

    H_O[r] = static_cast<float>(std::pow(2.0, -static_cast<double>(avg_height[panel_n - 1]) / c_norm));
  }

  return H_O;
}

} // namespace paocseq
