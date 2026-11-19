// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_LINALG_STRUMPACK_HPP
#define PALACE_LINALG_STRUMPACK_HPP

#include <mfem.hpp>

#if defined(MFEM_USE_STRUMPACK)

#include "linalg/operator.hpp"
#include "utils/iodata.hpp"

namespace palace
{

// Set the compression options of STRUMPACK, for the real-valued solvers and for the solver
// in complex arithmetic (StrumpackComplexSolver).
template <typename ScalarType>
void ConfigureStrumpackCompression(strumpack::SPOptions<ScalarType> &opts,
                                   SparseCompression compression, double lr_tol,
                                   int butterfly_l, int lossy_prec);

// Enable the STRUMPACK option which stores the fill-reducing permutation of the distributed
// solver and its inverse (two integers per matrix row) once per node, in an MPI shared
// memory window, instead of on every process (Solver.Linear.STRUMPACKSharedPermutation),
// for the real-valued solvers and for the solver in complex arithmetic. STRUMPACK does this
// after a sequential ordering (not ParMETIS or PT-Scotch) without compression, and falls
// back to a copy per process, with a warning, if the window cannot be allocated. The option
// is a patch in the STRUMPACK builds of Palace (superbuild and Spack package), and
// STRUMPACK defines no macro for it: this does nothing with a STRUMPACK which lacks it.
template <typename ScalarType>
inline void EnableStrumpackSharedPermutation(strumpack::SPOptions<ScalarType> &opts)
{
  if constexpr (requires { opts.enable_shared_permutation(); })
  {
    opts.enable_shared_permutation();
  }
}

//
// A wrapper for the STRUMPACK direct solver package.
//
template <typename StrumpackSolverType>
class StrumpackSolverBase : public StrumpackSolverType
{
private:
  MPI_Comm comm;

public:
  StrumpackSolverBase(MPI_Comm comm, SymbolicFactorization reorder,
                      SparseCompression compression, double lr_tol, int butterfly_l,
                      int lossy_prec, bool reorder_reuse, bool shared_permutation,
                      int print);

  StrumpackSolverBase(const IoData &iodata, MPI_Comm comm, int print)
    : StrumpackSolverBase(comm, iodata.solver.linear.sym_factorization,
                          iodata.solver.linear.strumpack_compression_type,
                          iodata.solver.linear.strumpack_lr_tol,
                          iodata.solver.linear.strumpack_butterfly_l,
                          iodata.solver.linear.strumpack_lossy_precision,
                          iodata.solver.linear.reorder_reuse,
                          iodata.solver.linear.strumpack_shared_permutation, print)
  {
  }

  void SetOperator(const Operator &op) override;

  void SetReorderReuse(bool reorder_reuse);
};

using StrumpackSolver = StrumpackSolverBase<mfem::STRUMPACKSolver>;

using StrumpackMixedPrecisionSolver =
    StrumpackSolverBase<mfem::STRUMPACKMixedPrecisionSolver>;

}  // namespace palace

#endif

#endif  // PALACE_LINALG_STRUMPACK_HPP
