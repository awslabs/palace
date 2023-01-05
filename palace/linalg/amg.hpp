// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_AMG_HPP
#define PALACE_AMG_HPP

#include <mfem.hpp>
#include "utils/iodata.hpp"

namespace palace
{

//
// A wrapper for Hypre's BoomerAMG solver.
//
class BoomerAmgSolver : public mfem::HypreBoomerAMG
{
private:
  // Helper function for setting common settings.
  void Init();

public:
  BoomerAmgSolver(int cycle_it = 1, int smooth_it = 1, int print = 0);
  BoomerAmgSolver(const IoData &iodata, int print)
    : BoomerAmgSolver((iodata.solver.linear.mat_gmg) ? 1 : iodata.solver.linear.mg_cycle_it,
                      iodata.solver.linear.mg_smooth_it, print)
  {
  }

  // Set the number of smoothing iterations to be performed at each level.
  void SetNumSweeps(int relax_sweeps);

  // Set the relaxation type on the coarsest level. Useful for specifying to not use a
  // direct solve when the coarse matrix may be singular(relax_type = 8 is the AMS
  // default).
  void SetCoarseRelaxType(int relax_type);
};

}  // namespace palace

#endif  // PALACE_AMG_HPP
