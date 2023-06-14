// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "mumps.hpp"

#if defined(MFEM_USE_MUMPS)

#include "linalg/rap.hpp"

namespace palace
{

MumpsSolver::MumpsSolver(MPI_Comm comm, mfem::MUMPSSolver::MatType sym,
                         config::LinearSolverData::SymFactType reorder, double blr_tol,
                         int print)
  : mfem::MUMPSSolver(comm)
{
  // Configure the solver (must be called before SetOperator).
  SetPrintLevel(print);
  SetMatrixSymType(sym);
  if (reorder == config::LinearSolverData::SymFactType::METIS)
  {
    SetReorderingStrategy(mfem::MUMPSSolver::METIS);
  }
  else if (reorder == config::LinearSolverData::SymFactType::PARMETIS)
  {
    SetReorderingStrategy(mfem::MUMPSSolver::PARMETIS);
  }
  else if (reorder == config::LinearSolverData::SymFactType::SCOTCH)
  {
    SetReorderingStrategy(mfem::MUMPSSolver::SCOTCH);
  }
  else if (reorder == config::LinearSolverData::SymFactType::PTSCOTCH)
  {
    SetReorderingStrategy(mfem::MUMPSSolver::PTSCOTCH);
  }
  else
  {
    SetReorderingStrategy(mfem::MUMPSSolver::AUTOMATIC);  // MUMPS should have good defaults
  }
  SetReorderingReuse(true);  // Repeated calls use same sparsity pattern
  if (blr_tol > 0.0)
  {
    SetBLRTol(blr_tol);
  }
}

void MumpsSolver::SetOperator(const Operator &op)
{
  const auto *PtAP = dynamic_cast<const ParOperator *>(&op);
  if (PtAP)
  {
    mfem::MUMPSSolver::SetOperator(const_cast<ParOperator *>(PtAP)->ParallelAssemble());
  }
  else
  {
    mfem::MUMPSSolver::SetOperator(op);
  }
}

}  // namespace palace

#endif
