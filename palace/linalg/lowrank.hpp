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
// and stored in blocks, typically one per port: each block packs its terms on the union of
// their nonzero local entries, so a product costs O(Σ_blocks K_block * nnz_block) local
// work and a single MPI reduction for all K terms. Blocks may share entries.
//
template <bool Conj>
class LowRankComplexUpdate
{
private:
  MPI_Comm comm;

  // Local size, number of terms (identical on all ranks) and their scalar weights.
  int n, K;
  std::vector<std::complex<double>> g;

  // A block with local entries: its terms (global indices), the local indices of its packed
  // entries, and the values of its terms there (term-major). The device dot products of its
  // terms use nb blocks of threads per term, with partials from offset, and its update
  // coefficients start at first.
  struct Block
  {
    std::vector<int> terms;
    mfem::Array<int> idx;
    Vector vr, vi;
    int nb, offset, first;
  };
  std::vector<Block> blocks;

  // Workspace for the block partials of the device dot products (2 per block of threads),
  // the local and global dot products (2 per term), and the update coefficients of the
  // terms of the local blocks (2 per term).
  mutable Vector partial, dots, coeff;

public:
  // The term vectors v have local size n, and term k belongs to the block with id block[k].
  LowRankComplexUpdate(MPI_Comm comm, int n, std::span<const ComplexVector *const> v,
                       std::span<const std::complex<double>> g, std::span<const int> block);

  [[nodiscard]] int NumTerms() const { return K; }
  [[nodiscard]] int NumBlocks() const { return static_cast<int>(blocks.size()); }
  [[nodiscard]] int NumPackedEntries() const;

  void AddMult(const ComplexVector &x, ComplexVector &y,
               const std::complex<double> a = 1.0) const;
};

}  // namespace palace

#endif  // PALACE_LINALG_LOW_RANK_HPP
