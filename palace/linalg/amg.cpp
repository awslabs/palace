// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "amg.hpp"

namespace palace
{

BoomerAmgSolver::BoomerAmgSolver(int cycle_it, int smooth_it, int print)
  : mfem::HypreBoomerAMG()
{
  SetPrintLevel((print > 1) ? print - 1 : 0);
  SetMaxIter(cycle_it);
  SetTol(0.0);

  // Set additional BoomerAMG options.
  double theta = 0.5;  // AMG strength parameter = 0.25 is 2D optimal (0.5-0.8 for 3D)
  int agg_levels = 1;  // Number of aggressive coarsening levels

  SetStrengthThresh(theta);
  SetAggressiveCoarsening(agg_levels);
  HYPRE_BoomerAMGSetNumSweeps(*this, smooth_it);

  // int coarse_relax_type = 8;  // l1-symm. GS (inexact coarse solve)
  // HYPRE_BoomerAMGSetCycleRelaxType(*this, coarse_relax_type, 3);
}

void BoomerAmgSolver::SetOperator(const ParOperator &op)
{
  mfem::HypreBoomerAMG::SetOperator(const_cast<ParOperator *>(&op)->ParallelAssemble());
}

}  // namespace palace
