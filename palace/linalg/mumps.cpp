// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "mumps.hpp"

#if defined(MFEM_USE_MUMPS)

namespace palace
{

MumpsSolver::MumpsSolver(MPI_Comm comm, mfem::MUMPSSolver::MatType sym, int sym_fact_type,
                         double blr_tol, int print)
  : mfem::MUMPSSolver(comm)
{
  // Configure the solver (must be called before SetOperator).
  SetPrintLevel(print);
  SetMatrixSymType(sym);
  if (sym_fact_type == 2)
  {
    SetReorderingStrategy(mfem::MUMPSSolver::PARMETIS);
  }
  else if (sym_fact_type == 1)
  {
    SetReorderingStrategy(mfem::MUMPSSolver::METIS);
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

void MumpsSolver::SetOperator(const ParOperator &op)
{
  mfem::MUMPSSolver::SetOperator(const_cast<ParOperator *>(&op)->ParallelAssemble());
}

}  // namespace palace

#endif
