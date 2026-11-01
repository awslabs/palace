// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_LINALG_HODLR_HPP
#define PALACE_LINALG_HODLR_HPP

#include <vector>

namespace palace
{

// Hierarchical off-diagonal low-rank (HODLR) representation of the symmetric interface
// operator S_E. The interface DOFs are permuted (recursive coordinate-median clustering);
// diagonal leaf blocks are stored dense, and each internal node's off-diagonal coupling
// S[I1,I2] is stored as a low-rank factor pair U V^T (symmetry gives the transpose block V
// U^T). Storage and apply are O(nG * (leaf + rank * log nG)) instead of O(nG^2). Replicated
// across ranks (small once compressed); serializable for MPI broadcast and file I/O.
struct Hodlr
{
  struct Leaf
  {
    int s, m;               // start (permuted), size
    std::vector<double> D;  // m x m, row-major
  };
  struct Block
  {
    int s1, m1, s2, m2, r;  // row/col starts+sizes (permuted), rank
    std::vector<double> U;  // m1 x r, row-major
    std::vector<double> V;  // m2 x r, row-major (B[I1,I2] ~= U V^T)
  };
  int n = 0;
  std::vector<int> perm;  // permuted position -> interface global index
  std::vector<Leaf> leaves;
  std::vector<Block> blocks;

  long long Storage() const;

  // y = S x, both in the permuted ordering (length n).
  void Mult(const double *xp, double *yp) const;

  // Flatten to a double buffer (indices packed as doubles) for MPI_Bcast and file I/O.
  std::vector<double> Serialize() const;

  static Hodlr Deserialize(const std::vector<double> &b);
};

// Recursively build a Hodlr from a dense symmetric S (n x n, row-major): split the index
// set `idx` at the median of its widest coordinate axis, compress the off-diagonal block
// via a truncated SVD to relative tolerance `tol`, recurse on the diagonal blocks. Leaf
// blocks (size <= `leaf`) are stored dense. `base` is the permuted offset of this block.
// `coords` is n x 3 (global-index order).
void BuildHodlr(const std::vector<double> &S, int n, const std::vector<int> &idx, int base,
                const std::vector<double> &coords, double tol, int leaf, Hodlr &h);

}  // namespace palace

#endif  // PALACE_LINALG_HODLR_HPP
