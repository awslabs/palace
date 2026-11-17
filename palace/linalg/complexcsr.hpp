// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_LINALG_COMPLEX_CSR_HPP
#define PALACE_LINALG_COMPLEX_CSR_HPP

#include <complex>
#include <cstddef>
#include <vector>
#include <mfem.hpp>
#include "linalg/operator.hpp"

namespace palace
{

//
// Local rows of a complex-valued operator, whose real and imaginary parts are
// HypreParMatrix or ParOperator objects with the same row distribution, merged into a
// complex-valued CSR matrix with global column indices (sorted within each row). Used by
// the sparse direct solvers which factor complex-valued matrices in complex arithmetic.
//
// With drop_small_entries (Solver.Linear.DropSmallEntries), the entries of a row which are
// smaller in magnitude than drop_small_entries_tol times the 2-norm of the row are not
// stored, except on the diagonal. It is the rule which MfemWrapperSolver applies to its
// real-valued matrix, here for the complex-valued entries: when the imaginary part is zero,
// the entries which are dropped are those which are dropped from the real part.
//
struct ComplexLocalCSR
{
  // Global number of rows and global index of the first local row.
  HYPRE_BigInt n_glob = 0, first_row = 0;

  // Local number of rows, row offsets (size n_loc + 1), global column indices, and values.
  int n_loc = 0;
  std::vector<HYPRE_BigInt> I, J;
  std::vector<std::complex<double>> V;

  explicit ComplexLocalCSR(const ComplexOperator &op, bool drop_small_entries = false);

  // Replaces hash by a hash of the sparsity pattern, and returns whether it has not changed
  // on any process: the solvers keep it to reuse the reordering for the next matrix with
  // the same sparsity pattern. Collective.
  bool SamePattern(MPI_Comm comm, std::size_t &hash) const;
};

// Whether the operator has no imaginary part, or an imaginary part whose parallel assembled
// matrix has no nonzero entry on any process, as for a model without losses (the imaginary
// part is assembled, with explicit zeros). Collective.
bool HasZeroImag(const ComplexOperator &op);

// Release the parallel assembled matrices of the ParOperator real and imaginary parts of
// the operator, for the solvers which copy the matrix (as MfemWrapperSolver does for the
// sparse direct solvers). They are assembled again if needed.
void ReleaseParallelAssemble(const ComplexOperator &op);

}  // namespace palace

#endif  // PALACE_LINALG_COMPLEX_CSR_HPP
