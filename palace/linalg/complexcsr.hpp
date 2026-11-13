// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_LINALG_COMPLEX_CSR_HPP
#define PALACE_LINALG_COMPLEX_CSR_HPP

#include <complex>
#include <vector>
#include <mfem.hpp>
#include "linalg/operator.hpp"

namespace palace
{

//
// Local rows of a complex-valued operator, whose real and imaginary parts are
// HypreParMatrix or ParOperator objects with the same row distribution, merged into a
// complex-valued CSR matrix with global column indices (sorted within each row). Used by
// the sparse direct solvers which factor complex-valued matrices in complex arithmetic. The
// parallel assembled matrices of ParOperator parts are released (the solvers copy the
// matrix).
//
struct ComplexLocalCSR
{
  // Global number of rows and global index of the first local row.
  HYPRE_BigInt n_glob = 0, first_row = 0;

  // Local number of rows, row offsets (size n_loc + 1), global column indices, and values.
  int n_loc = 0;
  std::vector<int> I;
  std::vector<HYPRE_BigInt> J;
  std::vector<std::complex<double>> V;

  explicit ComplexLocalCSR(const ComplexOperator &op);
};

}  // namespace palace

#endif  // PALACE_LINALG_COMPLEX_CSR_HPP
