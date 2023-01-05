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
  SetNumSweeps(smooth_it);
  Init();
}

void BoomerAmgSolver::Init()
{
  double theta = 0.5;  // AMG strength parameter = 0.25 is 2D optimal (0.5-0.8 for 3D)
  int agg_levels = 1;  // Number of aggressive coarsening levels
  SetStrengthThresh(theta);
  SetAggressiveCoarsening(agg_levels);
}

void BoomerAmgSolver::SetNumSweeps(int relax_sweeps)
{
  HYPRE_BoomerAMGSetNumSweeps(*this, relax_sweeps);
}

void BoomerAmgSolver::SetCoarseRelaxType(int relax_type)
{
  HYPRE_BoomerAMGSetCycleRelaxType(*this, relax_type, 3);
}

}  // namespace palace
