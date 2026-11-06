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
// and stored packed on their nonzero local entries, with the terms grouped by connected
// support (typically one group per port), so a product costs O(Σ_groups K_group *
// nnz_group) local work and a single MPI reduction for all K terms.
//
template <bool Conj>
class LowRankComplexUpdate
{
private:
  MPI_Comm comm;

  // Local size, number of terms (identical on all ranks) and their scalar weights.
  int n, K;
  std::vector<std::complex<double>> g;

  // The terms with local support, ordered by group. Group i holds the local terms
  // [group_terms[i], group_terms[i + 1]) (global indices perm[t]), the packed entries
  // [group_ents[i], group_ents[i + 1]) (local indices idx[j]), and their values
  // term-major from group_vals[i].
  mfem::Array<int> perm, group_terms, group_ents, group_vals, idx;
  Vector vr, vi;

  // Device lookup tables: the group of each local term and of each packed entry, and the
  // number of blocks per term of the dot products.
  mfem::Array<int> term_group, ent_group;
  int nb;

  // Workspace for the block partials of the device dot products (2 per block), the local
  // and global dot products (2 per term), and the update coefficients of the local terms (2
  // per term).
  mutable Vector partial, dots, coeff;

public:
  // The term vectors v have local size n.
  LowRankComplexUpdate(MPI_Comm comm, int n, std::span<const ComplexVector *const> v,
                       std::span<const std::complex<double>> g);

  [[nodiscard]] int NumTerms() const { return K; }
  [[nodiscard]] int NumGroups() const { return group_terms.Size() - 1; }
  [[nodiscard]] int NumPackedEntries() const { return idx.Size(); }

  void AddMult(const ComplexVector &x, ComplexVector &y,
               const std::complex<double> a = 1.0) const;
};

}  // namespace palace

#endif  // PALACE_LINALG_LOW_RANK_HPP
