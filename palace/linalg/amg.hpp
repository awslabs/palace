// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_LINALG_AMG_HPP
#define PALACE_LINALG_AMG_HPP

#include <mfem.hpp>
#include "utils/iodata.hpp"

namespace palace
{

//
// A wrapper for Hypre's BoomerAMG solver.
//
class BoomerAmgSolver : public mfem::HypreBoomerAMG
{
public:
  BoomerAmgSolver(int cycle_it = 1, int smooth_it = 1, bool agg_coarsen = true,
                  int print = 0);
  BoomerAmgSolver(const IoData &iodata, bool coarse_solver, int print)
    : BoomerAmgSolver(coarse_solver ? 1 : iodata.solver.linear.mg_cycle_it,
                      iodata.solver.linear.mg_smooth_it,
                      iodata.solver.linear.amg_agg_coarsen, print)
  {
  }
};

}  // namespace palace

#endif  // PALACE_LINALG_AMG_HPP
