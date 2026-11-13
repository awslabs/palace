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

// val += a X on a pattern from LowerTrianglePattern, for a matrix X with entries inside the
// pattern (or zero).
void AddToPattern(const mfem::HypreParMatrix &X, double a, const std::vector<int> &row_ptr,
                  const std::vector<int> &jcn, std::vector<double> &val);
void AddToPattern(const mfem::HypreParMatrix &X, std::complex<double> a,
                  const std::vector<int> &row_ptr, const std::vector<int> &jcn,
                  std::vector<std::complex<double>> &val);

}  // namespace palace

#endif  // PALACE_LINALG_LOWER_TRIANGLE_HPP
