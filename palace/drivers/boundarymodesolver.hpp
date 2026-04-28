// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_DRIVERS_BOUNDARY_MODE_SOLVER_HPP
#define PALACE_DRIVERS_BOUNDARY_MODE_SOLVER_HPP

#include <memory>
#include <mfem.hpp>
#include "drivers/basesolver.hpp"
#include "utils/geodata.hpp"

namespace palace
{

//
// Driver class for 2D waveguide mode analysis. Constructs a BoundaryModeOperator and a
// ModeEigenSolver on top of it, then orchestrates the solve and postprocessing. When
// configured with a 3D parent mesh and Solver.BoundaryMode.Attributes, PreprocessMesh
// extracts the 2D submesh during the read-mesh stage so the parallel pipeline (and AMR)
// treats the cross-section as the primary solve mesh.
//
class BoundaryModeSolver : public BaseSolver
{
public:
  BoundaryModeSolver(const IoData &iodata, bool root, int size = 0, int num_thread = 0,
                     const char *git_tag = nullptr);
  ~BoundaryModeSolver() override = default;

  double PreprocessMesh(std::unique_ptr<mfem::Mesh> &smesh, MPI_Comm comm) const override;

  std::pair<ErrorIndicator, long long int>
  Solve(const std::vector<std::unique_ptr<Mesh>> &mesh) const override;

private:
  // Non-null when the solve mesh was extracted from a 3D parent. Populated by
  // PreprocessMesh in nondim units (centroid pre-scaled by Lc), consumed by Solve for
  // material-tensor rotation and iodata-path projection. Mutable so the const
  // PreprocessMesh hook can populate it; write-once per solver instance.
  mutable std::unique_ptr<mesh::SubmeshFrame> frame_;
};

}  // namespace palace

#endif  // PALACE_DRIVERS_BOUNDARY_MODE_SOLVER_HPP
