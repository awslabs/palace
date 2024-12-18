// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_LINALG_MUMPS_HPP
#define PALACE_LINALG_MUMPS_HPP

#include <mfem.hpp>

#if defined(MFEM_USE_MUMPS)

#include "utils/iodata.hpp"

namespace palace
{

//
// A wrapper for the MUMPS direct solver package.
//
class MumpsSolver : public mfem::MUMPSSolver
{
public:
  MumpsSolver(MPI_Comm comm, mfem::MUMPSSolver::MatType sym,
              config::LinearSolverData::SymFactType reorder, double blr_tol, int print);
  MumpsSolver(const IoData &iodata, MPI_Comm comm, int print)
    : MumpsSolver(comm,
                  (iodata.solver.linear.pc_mat_shifted ||
                   iodata.problem.type == config::ProblemData::Type::TRANSIENT ||
                   iodata.problem.type == config::ProblemData::Type::ELECTROSTATIC ||
                   iodata.problem.type == config::ProblemData::Type::MAGNETOSTATIC)
                      ? mfem::MUMPSSolver::SYMMETRIC_POSITIVE_DEFINITE
                      : mfem::MUMPSSolver::SYMMETRIC_INDEFINITE,
                  iodata.solver.linear.sym_fact_type,
                  (iodata.solver.linear.strumpack_compression_type ==
                   config::LinearSolverData::CompressionType::BLR)
                      ? iodata.solver.linear.strumpack_lr_tol
                      : 0.0,
                  print)
  {
  }
};

}  // namespace palace

#endif

#endif  // PALACE_LINALG_MUMPS_HPP
