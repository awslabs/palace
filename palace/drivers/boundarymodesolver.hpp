// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_DRIVERS_BOUNDARY_MODE_SOLVER_HPP
#define PALACE_DRIVERS_BOUNDARY_MODE_SOLVER_HPP

#include "drivers/basesolver.hpp"
#include "models/boundarymodeoperator.hpp"

namespace palace
{

//
// Driver class for 2D waveguide mode analysis. Constructs a BoundaryModeOperator (which
// owns the FE spaces, material operator, and eigenvalue solver), then orchestrates the
// solve and postprocessing. When configured with a 3D parent mesh and
// Solver.BoundaryMode.Attributes, PreprocessMesh extracts the 2D submesh during the
// read-mesh stage so the parallel pipeline (and AMR) treats the cross-section as the
// primary solve mesh.
//
class BoundaryModeSolver : public BaseSolver
{
public:
  using BaseSolver::BaseSolver;

  void PreprocessMesh(mesh::SerialMesh &smesh, MPI_Comm comm) const override;

  std::pair<ErrorIndicator, long long int>
  Solve(const std::vector<std::unique_ptr<Mesh>> &mesh) const override;

private:
  // Projection frame for the extracted submesh. Populated by PreprocessMesh when a 3D
  // cross-section is being analysed; empty otherwise. Mutable so the const
  // PreprocessMesh hook can write to it; read by Solve.
  mutable SubmeshFrame frame_;
  mutable bool from_submesh_ = false;
};

}  // namespace palace

#endif  // PALACE_DRIVERS_BOUNDARY_MODE_SOLVER_HPP
