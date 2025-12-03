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
  MumpsSolver(MPI_Comm comm, mfem::MUMPSSolver::MatType sym, SymbolicFactorization reorder,
              double blr_tol, bool reorder_reuse, int print);
  MumpsSolver(const IoData &iodata, MPI_Comm comm, int print)
    : MumpsSolver(
          comm,
          (iodata.solver.linear.pc_mat_shifted ||
           iodata.problem.type == ProblemType::TRANSIENT ||
           iodata.problem.type == ProblemType::ELECTROSTATIC ||
           iodata.problem.type == ProblemType::MAGNETOSTATIC)
              ? mfem::MUMPSSolver::SYMMETRIC_POSITIVE_DEFINITE
          : iodata.boundaries.periodic.wave_vector == std::array<double, 3>{0.0, 0.0, 0.0}
              ? mfem::MUMPSSolver::SYMMETRIC_INDEFINITE
              : mfem::MUMPSSolver::UNSYMMETRIC,
          iodata.solver.linear.sym_factorization,
          (iodata.solver.linear.strumpack_compression_type == SparseCompression::BLR)
              ? iodata.solver.linear.strumpack_lr_tol
              : 0.0,
          iodata.solver.linear.reorder_reuse, print)
  {
  }

  void SetReorderReuse(bool reorder_reuse);
};

}  // namespace palace

#endif

#endif  // PALACE_LINALG_MUMPS_HPP
