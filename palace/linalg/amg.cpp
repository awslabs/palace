// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "amg.hpp"

#include <algorithm>

#include <_hypre_parcsr_mv.h>

namespace palace
{

BoomerAmgSolver::BoomerAmgSolver(int cycle_it, int smooth_it, bool agg_coarsen, int print)
  : mfem::HypreBoomerAMG()
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

BoomerAmgSolver::~BoomerAmgSolver()
{
  for (auto *v : {B_mv, X_mv})
  {
    if (v)
    {
      hypre_ParVectorDestroy((hypre_ParVector *)v);
    }
  }
}

void BoomerAmgSolver::ArrayMult(const mfem::Array<const mfem::Vector *> &X,
                                mfem::Array<mfem::Vector *> &Y) const
{
  // Multicomponent vectors are only supported by the host relaxation and coarse solve.
  const int nv = X.Size();
  if (nv < 2 || iterative_mode || mv_failed ||
      mfem::Device::Allows(mfem::Backend::DEVICE_MASK))
  {
    mfem::HypreBoomerAMG::ArrayMult(X, Y);
    return;
  }
  MFEM_VERIFY(A, "HypreParMatrix A is missing for BoomerAmgSolver::ArrayMult!");
  if (!setup_called)
  {
    Setup(*X[0], *Y[0]);
  }
  hypre_ParCSRMatrix *hA = *A;
  const int n = A->Height();
  auto *hB = (hypre_ParVector *)B_mv;
  if (!hB || hypre_ParVectorNumVectors(hB) != nv ||
      hypre_VectorSize(hypre_ParVectorLocalVector(hB)) != n)
  {
    for (auto *v : {&B_mv, &X_mv})
    {
      if (*v)
      {
        hypre_ParVectorDestroy((hypre_ParVector *)*v);
      }
      auto *hv = hypre_ParMultiVectorCreate(hypre_ParCSRMatrixComm(hA),
                                            hypre_ParCSRMatrixGlobalNumRows(hA),
                                            hypre_ParCSRMatrixRowStarts(hA), nv);
      hypre_ParVectorInitialize(hv);
      *v = (HYPRE_ParVector)hv;
    }
    hB = (hypre_ParVector *)B_mv;
  }
  auto *hX = (hypre_ParVector *)X_mv;
  HYPRE_Real *b = hypre_VectorData(hypre_ParVectorLocalVector(hB));
  for (int k = 0; k < nv; k++)
  {
    const auto *XK = X[k]->HostRead();
    std::copy(XK, XK + n, b + k * n);
  }
  hypre_ParVectorSetZeros(hX);
  if (HYPRE_BoomerAMGSolve(*this, *A, B_mv, X_mv))
  {
    // Not supported for this configuration (relaxation type): use separate solves.
    hypre_error_flag = 0;
    mv_failed = true;
    mfem::HypreBoomerAMG::ArrayMult(X, Y);
    return;
  }
  const HYPRE_Real *x = hypre_VectorData(hypre_ParVectorLocalVector(hX));
  for (int k = 0; k < nv; k++)
  {
    auto *YK = Y[k]->HostWrite();
    std::copy(x + k * n, x + (k + 1) * n, YK);
  }
}

}  // namespace palace
