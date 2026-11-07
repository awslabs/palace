// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "amg.hpp"

#include <algorithm>

#include <_hypre_parcsr_mv.h>

namespace palace
{

BoomerAmgSolver::BoomerAmgSolver(int cycle_it, int smooth_it, bool agg_coarsen, int print)
  : mfem::HypreBoomerAMG(), cycle_it(cycle_it)
{
  HYPRE_BoomerAMGSetPrintLevel(*this, (print > 1) ? print - 1 : 0);
  HYPRE_BoomerAMGSetMaxIter(*this, cycle_it);
  HYPRE_BoomerAMGSetTol(*this, 0.0);

  // Set additional BoomerAMG options.
  int agg_levels = agg_coarsen ? 1 : 0;  // Number of aggressive coarsening levels
  double theta = 0.5;  // AMG strength parameter = 0.25 is 2D optimal (0.5-0.8 for 3D)
  int relax_type = 8;  // 8 = l1-symm. GS, 13 = l1-GS, 18 = l1-Jacobi, 16 = Chebyshev
  if (mfem::Device::Allows(mfem::Backend::DEVICE_MASK))
  {
    // Modify options for GPU-supported features.
    agg_levels = 0;
    relax_type = 18;
  }

  HYPRE_BoomerAMGSetAggNumLevels(*this, agg_levels);
  HYPRE_BoomerAMGSetStrongThreshold(*this, theta);
  HYPRE_BoomerAMGSetRelaxType(*this, relax_type);
  HYPRE_BoomerAMGSetNumSweeps(*this, smooth_it);

  // int coarse_relax_type = 8;  // l1-symm. GS (inexact coarse solve)
  // HYPRE_BoomerAMGSetCycleRelaxType(*this, coarse_relax_type, 3);
}

void BoomerAmgSolver::Mult(const mfem::HypreParVector &b, mfem::HypreParVector &x) const
{
  if (!iterative_mode)
  {
    x.HypreWrite();
    hypre_ParVectorSetZeros(x);
  }
  mfem::HypreBoomerAMG::Mult(b, x);
}

void BoomerAmgSolver::SetOperator(const mfem::Operator &op)
{
  // When the solver already has an operator, mfem::HypreBoomerAMG::SetOperator recreates
  // the hypre solver with its coarsening and relaxation options, but with one V-cycle.
  mfem::HypreBoomerAMG::SetOperator(op);
  HYPRE_BoomerAMGSetMaxIter(*this, cycle_it);
  B_mv.reset();
  X_mv.reset();
}

bool BoomerAmgSolver::MultiComponentSupported()
{
#if defined(HYPRE_BOOMERAMG_MULTICOMPONENT) && !defined(HYPRE_USING_PERSISTENT_COMM)
  return hypre_NumThreads() == 1 && !mfem::Device::Allows(mfem::Backend::DEVICE_MASK);
#else
  return false;
#endif
}

void BoomerAmgSolver::ArrayMult(const mfem::Array<const mfem::Vector *> &X,
                                mfem::Array<mfem::Vector *> &Y) const
{
  MFEM_VERIFY(X.Size() == Y.Size(),
              "Mismatch in the number of vectors in BoomerAmgSolver::ArrayMult!");
  const int nv = X.Size();
  if (nv < 2 || iterative_mode || !MultiComponentSupported())
  {
    mfem::HypreBoomerAMG::ArrayMult(X, Y);
    return;
  }
  MFEM_VERIFY(A, "HypreParMatrix A is missing for BoomerAmgSolver::ArrayMult!");
  const int n = A->Height();
  for (int k = 0; k < nv; k++)
  {
    MFEM_ASSERT(X[k] && X[k]->Size() == n && Y[k] && Y[k]->Size() == n,
                "Vector size mismatch in BoomerAmgSolver::ArrayMult!");
  }
  if (!setup_called)
  {
    Setup(*X[0], *Y[0]);
  }

  // The multicomponent vectors are reset with the operator.
  hypre_ParCSRMatrix *hA = *A;
  if (!B_mv || hypre_ParVectorNumVectors(B_mv.get()) != nv)
  {
    for (auto *v : {&B_mv, &X_mv})
    {
      v->reset(hypre_ParMultiVectorCreate(hypre_ParCSRMatrixComm(hA),
                                          hypre_ParCSRMatrixGlobalNumRows(hA),
                                          hypre_ParCSRMatrixRowStarts(hA), nv));
      hypre_ParVectorInitialize(v->get());
    }
  }
  HYPRE_Real *b = hypre_VectorData(hypre_ParVectorLocalVector(B_mv.get()));
  for (int k = 0; k < nv; k++)
  {
    const auto *XK = X[k]->HostRead();
    std::copy(XK, XK + n, b + k * n);
  }
  hypre_ParVectorSetZeros(X_mv.get());

  // Same error handling as mfem::HypreSolver::Mult.
  const HYPRE_Int err = HYPRE_BoomerAMGSolve(*this, *A, (HYPRE_ParVector)B_mv.get(),
                                             (HYPRE_ParVector)X_mv.get());
  if (error_mode == WARN_HYPRE_ERRORS)
  {
    if (err)
    {
      MFEM_WARNING("Error during solve! Error code: " << err);
    }
  }
  else if (error_mode == ABORT_HYPRE_ERRORS)
  {
    MFEM_VERIFY(!err, "Error during solve! Error code: " << err);
  }
  hypre_error_flag = 0;

  const HYPRE_Real *x = hypre_VectorData(hypre_ParVectorLocalVector(X_mv.get()));
  for (int k = 0; k < nv; k++)
  {
    auto *YK = Y[k]->HostWrite();
    std::copy(x + k * n, x + (k + 1) * n, YK);
  }
}

}  // namespace palace
