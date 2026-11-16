// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_LINALG_LOW_RANK_HPP
#define PALACE_LINALG_LOW_RANK_HPP

#include <complex>
#include <span>
#include <vector>
#include <mfem.hpp>
#include "linalg/vector.hpp"
#include "utils/communication.hpp"

namespace palace
{

//
// Sum of complex rank-one updates y += a Σ_k g_k v_k (v_kᴴ x) (Conj = true) or
// y += a Σ_k g_k v_k (v_kᵀ x) (Conj = false). The term vectors are copied at construction
// and stored packed on the union of their nonzero local entries, so a product costs
// O(K * nnz) local work and a single MPI reduction for all K terms.
//
template <bool Conj>
class LowRankComplexUpdate
{
private:
  MPI_Comm comm;

  // Number of terms (identical on all ranks) and their scalar weights.
  int K;
  std::vector<std::complex<double>> g;

  // Local indices where any term is nonzero, and the term values there (K x nnz,
  // term-major).
  mfem::Array<int> idx;
  Vector vr, vi;

  // Workspace for the local and global partial dot products (4 per term) and the combined
  // update coefficients (2 per term).
  mutable Vector dots, coeff;

public:
  LowRankComplexUpdate(MPI_Comm comm, std::span<const ComplexVector *const> v,
                       std::span<const std::complex<double>> g);

  [[nodiscard]] int NumTerms() const { return K; }
  [[nodiscard]] int NumPackedEntries() const { return idx.Size(); }

  void AddMult(const ComplexVector &x, ComplexVector &y,
               const std::complex<double> a = 1.0) const;
};

}  // namespace palace

#endif  // PALACE_LINALG_LOW_RANK_HPP
