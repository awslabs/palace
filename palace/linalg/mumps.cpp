// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "mumps.hpp"

#if defined(MFEM_USE_MUMPS)

namespace palace
{

MumpsSolver::MumpsSolver(MPI_Comm comm, mfem::MUMPSSolver::MatType sym,
                         SymbolicFactorization reorder, double blr_tol, bool reorder_reuse,
                         int print)
  : mfem::MUMPSSolver(comm)
{
  // Configure the solver (must be called before SetOperator).
  SetPrintLevel(print);
  SetMatrixSymType(sym);
  switch (reorder)
  {
    case SymbolicFactorization::METIS:
      SetReorderingStrategy(mfem::MUMPSSolver::METIS);
      break;
    case SymbolicFactorization::PARMETIS:
      SetReorderingStrategy(mfem::MUMPSSolver::PARMETIS);
      break;
    case SymbolicFactorization::SCOTCH:
      SetReorderingStrategy(mfem::MUMPSSolver::SCOTCH);
      break;
    case SymbolicFactorization::PTSCOTCH:
      SetReorderingStrategy(mfem::MUMPSSolver::PTSCOTCH);
      break;
    case SymbolicFactorization::PORD:
      SetReorderingStrategy(mfem::MUMPSSolver::PORD);
      break;
    case SymbolicFactorization::AMD:
    case SymbolicFactorization::RCM:
      SetReorderingStrategy(mfem::MUMPSSolver::AMD);
      break;
    case SymbolicFactorization::DEFAULT:
      SetReorderingStrategy(mfem::MUMPSSolver::AUTOMATIC);  // Should have good default
      break;
  }
  SetReorderingReuse(reorder_reuse);  // If true repeated calls use same sparsity pattern
  if (blr_tol > 0.0)
  {
    SetBLRTol(blr_tol);
  }
}

}  // namespace palace

#endif
