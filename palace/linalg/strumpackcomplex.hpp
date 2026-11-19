// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_LINALG_STRUMPACK_COMPLEX_HPP
#define PALACE_LINALG_STRUMPACK_COMPLEX_HPP

#include <mfem.hpp>

#if defined(MFEM_USE_STRUMPACK)

#include <complex>
#include <cstddef>
#include <memory>
#include <vector>
#include "linalg/operator.hpp"
#include "linalg/solver.hpp"
#include "utils/labels.hpp"

namespace strumpack
{

template <typename scalar_t, typename integer_t>
class SparseSolverMPIDist;

}  // namespace strumpack

namespace palace
{

//
// STRUMPACK sparse direct solver for complex-valued operators in complex arithmetic. Its
// factors need about 40% less memory than those of the equivalent real-valued system of
// twice the size used by MfemWrapperSolver. The compression options and the node-shared
// permutation are those of StrumpackSolver, and small entries are dropped as described for
// ComplexLocalCSR. The integer type is HYPRE_BigInt, as for StrumpackSolver
// (mfem::STRUMPACKSolver).
//
class StrumpackComplexSolver : public Solver<ComplexOperator>
{
private:
  MPI_Comm comm;
  SymbolicFactorization reorder;
  SparseCompression compression;
  double lr_tol;
  int butterfly_l, lossy_prec;
  bool drop_small_entries, reorder_reuse;
  int print;
  std::unique_ptr<strumpack::SparseSolverMPIDist<std::complex<double>, HYPRE_BigInt>>
      solver;
  mutable std::vector<std::complex<double>> b, x;

  // Hash of the sparsity pattern of the last factored matrix, for the reuse of the
  // reordering.
  std::size_t pattern = 0;

public:
  StrumpackComplexSolver(MPI_Comm comm, SymbolicFactorization reorder,
                         SparseCompression compression, double lr_tol, int butterfly_l,
                         int lossy_prec, bool drop_small_entries, bool reorder_reuse,
                         int print);
  ~StrumpackComplexSolver() override;

  // Access the underlying solver (after SetOperator).
  const auto &GetSolver() const { return *solver; }

  void SetOperator(const ComplexOperator &op) override;

  void Mult(const ComplexVector &x, ComplexVector &y) const override;
};

}  // namespace palace

#endif

#endif  // PALACE_LINALG_STRUMPACK_COMPLEX_HPP
