// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_LINALG_STRUMPACK_HPP
#define PALACE_LINALG_STRUMPACK_HPP

#include <mfem.hpp>

#if defined(MFEM_USE_STRUMPACK)

#include "utils/iodata.hpp"

namespace palace
{

//
// A wrapper for the STRUMPACK direct solver package.
//
template <typename StrumpackSolverType>
class StrumpackSolverBase : public StrumpackSolverType
{
private:
  MPI_Comm comm;

public:
  StrumpackSolverBase(MPI_Comm comm, config::LinearSolverData::SymFactType reorder,
                      config::LinearSolverData::CompressionType compression, double lr_tol,
                      int butterfly_l, int lossy_prec, int print);

  StrumpackSolverBase(MPI_Comm comm, const IoData &iodata, int print)
    : StrumpackSolverBase(comm, iodata.solver.linear.sym_fact_type,
                          iodata.solver.linear.strumpack_compression_type,
                          iodata.solver.linear.strumpack_lr_tol,
                          iodata.solver.linear.strumpack_butterfly_l,
                          iodata.solver.linear.strumpack_lossy_precision, print_lvl)
  {
  }

  // Sets matrix associated with the STRUMPACK solver.
  void SetOperator(const mfem::Operator &op) override;
};

using StrumpackSolver = StrumpackSolverBase<mfem::STRUMPACKSolver>;

#if STRUMPACK_VERSION_MAJOR >= 6 && STRUMPACK_VERSION_MINOR >= 3 && \
    STRUMPACK_VERSION_PATCH > 1
using StrumpackMixedPrecisionSolver =
    StrumpackSolverBase<mfem::STRUMPACKMixedPrecisionSolver>;
#endif

}  // namespace palace

#endif

#endif  // PALACE_LINALG_STRUMPACK_HPP
