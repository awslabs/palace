// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "mumps.hpp"

#if defined(MFEM_USE_MUMPS)

#include "utils/iodata.hpp"

namespace palace

{

namespace
{

mfem::MUMPSSolver::MatType GetMumpsMatType(MatrixSymmetry sym)
{
  switch (sym)
  {
    case MatrixSymmetry::SPD:
      return mfem::MUMPSSolver::SYMMETRIC_POSITIVE_DEFINITE;
    case MatrixSymmetry::SYMMETRIC:
      return mfem::MUMPSSolver::SYMMETRIC_INDEFINITE;
    case MatrixSymmetry::UNSYMMETRIC:
      return mfem::MUMPSSolver::UNSYMMETRIC;
  }
  return mfem::MUMPSSolver::UNSYMMETRIC;
}

}  // namespace

MumpsSolver::MumpsSolver(MPI_Comm comm, MatrixSymmetry sym, SymbolicFactorization reorder,
                         double blr_tol, bool reorder_reuse, int print)
  : mfem::MUMPSSolver(comm)
{
  // Configure the solver (must be called before SetOperator).
  SetPrintLevel(print);
  SetMatrixSymType(GetMumpsMatType(sym));
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

MumpsSolver::MumpsSolver(const IoData &iodata, MPI_Comm comm, int print)
  : MumpsSolver(comm, iodata.solver.linear.pc_mat_sym,
                iodata.solver.linear.sym_factorization,
                (iodata.solver.linear.strumpack_compression_type == SparseCompression::BLR)
                    ? iodata.solver.linear.strumpack_lr_tol
                    : 0.0,
                iodata.solver.linear.reorder_reuse, print)
{
}

void MumpsSolver::SetReorderReuse(bool reorder_reuse)
{
  SetReorderingReuse(reorder_reuse);  // If true repeated calls use same sparsity pattern
}

}  // namespace palace

#endif
