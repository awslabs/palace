// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_HYPRE_HPP
#define PALACE_HYPRE_HPP

#include <mfem.hpp>

namespace mfem
{

// Convenience wrapper for casting away the const on the pointers and dispatching onto the
// original function that has the argument type: mfem::Array2D<mfem::HypreParMatrix*> &.
mfem::HypreParMatrix *
HypreParMatrixFromBlocks(mfem::Array2D<const mfem::HypreParMatrix *> &blocks,
                         mfem::Array2D<double> *blockCoeff = nullptr);

}  // namespace mfem

namespace palace::hypre
{

//
// Extensions to Hypre linear algebra routines.
//

// Version 1: Eliminates (including from the sparsity pattern) the given list of
//   rows/columns from the square matrix and sets the diagonal value according to
//   diag_policy. A combination of mfem::HypreParMatrix:: EliminateRowsCols and
//   hypre_ParCSRMatrixDropSmallEntriesHost. Specialized for host operation currently.
// Version 2: A mfem::HypreParMatrix::EliminateBC with option to specify a general scalar
//   for eliminated rows.
// The specified rows/columns should be in local numbering.
enum class DiagonalPolicy
{
  USER,
  ZERO,
  ONE,
  KEEP
};
void hypreParCSREliminateRowsCols(hypre_ParCSRMatrix *A, const mfem::Array<int> &rows_cols,
                                  hypre::DiagonalPolicy diag_policy,
                                  HYPRE_Complex diag = 0.0, bool ignore_rows = false);
void hypreParCSREliminateRowsColsv2(hypre_ParCSRMatrix *A,
                                    const mfem::Array<int> &rows_cols,
                                    hypre::DiagonalPolicy diag_policy,
                                    HYPRE_Complex diag = 0.0, bool ignore_rows = false);

// Eliminates (zeros) the given list of rows (columns), and also eliminates all columns
// (rows) which contain a nonzero in the specified rows (columns) to be eliminated. From
// Hypre's hypre_AMESetup. Returns as a new matrix (leaves the old matrix intact). The
// specified rows (columns) should be in local numbering.
hypre_ParCSRMatrix *hypreParCSREliminateRowsWithCols(hypre_ParCSRMatrix *A,
                                                     const mfem::Array<int> &rows);
hypre_ParCSRMatrix *hypreParCSREliminateColsWithRows(hypre_ParCSRMatrix *A,
                                                     const mfem::Array<int> &cols);

// Copy the entries from A into B, for sparsity(B) ⊆ sparsity(A).
void hypreParCSRCopy(hypre_ParCSRMatrix *A, hypre_ParCSRMatrix *B);

// Get the row sums (with absolute value) of the local rows of the matrix.
void hypreParCSRRowSums(hypre_ParCSRMatrix *A, mfem::Vector &rowsums);

// Compute the matrix infinity norm for a complex matrix stored with separate real and
// imaginary parts, for sparsity(Ai) ⊆ sparsity(Ar).
void hypreParCSRInfNorm(hypre_ParCSRMatrix *Ar, hypre_ParCSRMatrix *Ai, HYPRE_Real *norm);

}  // namespace palace::hypre

#endif  // PALACE_HYPRE_HPP
