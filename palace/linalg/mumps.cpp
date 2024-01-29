// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "mumps.hpp"

#if defined(MFEM_USE_MUMPS)

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
  switch (reorder)
  {
    case config::LinearSolverData::SymFactType::METIS:
      SetReorderingStrategy(mfem::MUMPSSolver::METIS);
      break;
    case config::LinearSolverData::SymFactType::PARMETIS:
      SetReorderingStrategy(mfem::MUMPSSolver::PARMETIS);
      break;
    case config::LinearSolverData::SymFactType::SCOTCH:
      SetReorderingStrategy(mfem::MUMPSSolver::SCOTCH);
      break;
    case config::LinearSolverData::SymFactType::PTSCOTCH:
      SetReorderingStrategy(mfem::MUMPSSolver::PTSCOTCH);
      break;
    case config::LinearSolverData::SymFactType::PORD:
      SetReorderingStrategy(mfem::MUMPSSolver::PORD);
      break;
    case config::LinearSolverData::SymFactType::AMD:
    case config::LinearSolverData::SymFactType::RCM:
      SetReorderingStrategy(mfem::MUMPSSolver::AMD);
      break;
    case config::LinearSolverData::SymFactType::DEFAULT:
      SetReorderingStrategy(mfem::MUMPSSolver::AUTOMATIC);  // Should have good default
      break;
  }
  SetReorderingReuse(true);  // Repeated calls use same sparsity pattern
  if (blr_tol > 0.0)
  {
    SetBLRTol(blr_tol);
  }
}

}  // namespace palace

#endif
