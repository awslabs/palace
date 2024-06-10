// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "amg.hpp"

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

}  // namespace palace
