// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_LINALG_LOWER_TRIANGLE_HPP
#define PALACE_LINALG_LOWER_TRIANGLE_HPP

#include <complex>
#include <vector>
#include <mfem.hpp>

namespace palace
{

// Lower-triangle COO pattern (1-based global indices, one entry per position, columns
// ascending in each row) of the local rows of the union of the given symmetric matrices
// (same row distribution; null ones skipped). The entries of each matrix in pattern order
// are returned in vals (empty for a null one), and the first pattern entry of each local
// row in row_ptr (n_loc + 1).
void LowerTrianglePattern(const std::vector<const mfem::HypreParMatrix *> &A,
                          std::vector<int> &irn, std::vector<int> &jcn,
                          std::vector<std::vector<double>> &vals,
                          std::vector<int> &row_ptr);

// The rows and columns of a symmetric operator of local true DOFs it keeps (the others
// pinned), numbered by rank (the system): the lower-triangle pattern of the union of its
// parts (same row distribution; null ones skipped) on the kept DOFs, as COO entries in
// system indices (1-based), with the entries of each part on it, so that the system's
// values are a combination of its parts. Collective.
class SideSystem
{
public:
  SideSystem(const std::vector<const mfem::HypreParMatrix *> &parts,
             const std::vector<char> &keep);

  // The system indices (0-based) of global true DOFs (replicated, in their order; -1 for
  // one not kept), each rank mapping those it owns. Collective.
  std::vector<HYPRE_BigInt> Indices(const std::vector<HYPRE_BigInt> &tdofs) const;

  // Σ_p coef[p] parts[p] on the pattern (null parts skipped).
  template <typename T>
  std::vector<T> Values(const std::vector<T> &coef) const;

  // The local true DOF of each local row of the system, and the global size of the system.
  std::vector<int> rows;
  HYPRE_BigInt n_glob = 0;
  // The pattern: its entries in system indices, and on the local true DOFs (as
  // LowerTrianglePattern, with global true-DOF columns), for AddToPattern.
  std::vector<int> irn, jcn, row_ptr, tdof_jcn;
  std::vector<std::vector<double>> parts;

private:
  MPI_Comm comm;
  HYPRE_BigInt tstart = 0;
  std::vector<HYPRE_BigInt> index;  // the system index of each local true DOF (-1: pinned)
};

// val += a X on a pattern from LowerTrianglePattern, for a matrix X with entries inside the
// pattern (or zero).
void AddToPattern(const mfem::HypreParMatrix &X, double a, const std::vector<int> &row_ptr,
                  const std::vector<int> &jcn, std::vector<double> &val);
void AddToPattern(const mfem::HypreParMatrix &X, std::complex<double> a,
                  const std::vector<int> &row_ptr, const std::vector<int> &jcn,
                  std::vector<std::complex<double>> &val);

}  // namespace palace

#endif  // PALACE_LINALG_LOWER_TRIANGLE_HPP
