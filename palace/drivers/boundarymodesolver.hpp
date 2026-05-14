// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_DRIVERS_BOUNDARY_MODE_SOLVER_HPP
#define PALACE_DRIVERS_BOUNDARY_MODE_SOLVER_HPP

#include <memory>
#include <mfem.hpp>
#include "drivers/basesolver.hpp"

namespace palace
{

// Driver class for 2D waveguide mode analysis.
class BoundaryModeSolver : public BaseSolver
{
public:
  BoundaryModeSolver(const IoData &iodata, bool root, int size = 0, int num_thread = 0,
                     const char *git_tag = nullptr);

  // If original mesh is 3D, extract 2D submesh from boundary attributes and changes
  // frame of IoData to match 2D coordinate system.
  void Preprocess(IoData &iodata, std::unique_ptr<mfem::Mesh> &smesh,
                  MPI_Comm comm) const override;

  std::pair<ErrorIndicator, long long int>
  Solve(const std::vector<std::unique_ptr<Mesh>> &mesh) const override;
};

}  // namespace palace

#endif  // PALACE_DRIVERS_BOUNDARY_MODE_SOLVER_HPP
