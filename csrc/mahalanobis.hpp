// -*- C++ -*-
//===------------------------- mahalanobis.hpp -----------------------------===//
// Part of the aocseq / paocseq Project, under the MIT license.
//===----------------------------------------------------------------------===//
#pragma once

#include <vector>

namespace paocseq {

// Mahalanobis distance of each query cell to a fixed reference panel's
// mean/covariance. See mahalanobis.cpp for full documentation. reference is
// n_genes x n_ref, query is n_genes x n_query, both gene-major (row-major
// with genes as the outer dimension). Returns a vector of length n_query.
std::vector<float> GetMahalanobis(int n_genes, int n_ref, int n_query,
                                   const std::vector<float> &reference,
                                   const std::vector<float> &query);

} // namespace paocseq
