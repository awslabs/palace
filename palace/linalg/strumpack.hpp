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

  StrumpackSolverBase(const IoData &iodata, MPI_Comm comm, int print)
    : StrumpackSolverBase(comm, iodata.solver.linear.sym_fact_type,
                          iodata.solver.linear.strumpack_compression_type,
                          iodata.solver.linear.strumpack_lr_tol,
                          iodata.solver.linear.strumpack_butterfly_l,
                          iodata.solver.linear.strumpack_lossy_precision, print)
  {
  }

  void SetOperator(const Operator &op) override;
};

using StrumpackSolver = StrumpackSolverBase<mfem::STRUMPACKSolver>;

using StrumpackMixedPrecisionSolver =
    StrumpackSolverBase<mfem::STRUMPACKMixedPrecisionSolver>;

}  // namespace palace

#endif

#endif  // PALACE_LINALG_STRUMPACK_HPP
