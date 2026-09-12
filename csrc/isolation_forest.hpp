// -*- C++ -*-
//===----------------------- isolation_forest.hpp ---------------------------===//
// Part of the aocseq / paocseq Project, under the MIT license.
//===----------------------------------------------------------------------===//
#pragma once

#include <vector>

namespace paocseq {

// Kurtosis-driven isolation-forest outlier score of each query cell
// relative to a reference panel. See isolation_forest.cpp for full
// documentation. reference is n_genes x n_ref, query is n_genes x n_query,
// both gene-major. Returns a vector of length n_query, each entry in (0, 1].
std::vector<float> IsoForest(int num_trees, int n_genes, int n_ref, int n_query,
                              const std::vector<float> &reference, const std::vector<float> &query,
                              int max_height, unsigned int seed);

} // namespace paocseq
