// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_DRIVERS_BOUNDARY_MODE_SOLVER_HPP
#define PALACE_DRIVERS_BOUNDARY_MODE_SOLVER_HPP

#include <memory>
#include <mfem.hpp>
#include "drivers/basesolver.hpp"

namespace palace
{

// Tangent frame of a 2D cross-section in the parent 3D coordinate system. Populated by
// BoundaryModeSolver::PreprocessMesh during submesh extraction and held by the driver
// for (a) rotating material tensors into the tangent plane when constructing the
// MaterialOperator in Solve, and (b) projecting 3D iodata path coordinates onto the 2D
// local frame. A direct-2D problem has no frame (frame_ stays null).
struct SubmeshFrame
{
  mfem::Vector centroid, e1, e2, normal;
};

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
  BoundaryModeSolver(const IoData &iodata, bool root, int size = 0, int num_thread = 0,
                     const char *git_tag = nullptr);
  ~BoundaryModeSolver() override = default;

  double PreprocessMesh(std::unique_ptr<mfem::Mesh> &smesh, MPI_Comm comm) const override;

  std::pair<ErrorIndicator, long long int>
  Solve(const std::vector<std::unique_ptr<Mesh>> &mesh) const override;

private:
  // Non-null when the solve mesh was extracted from a 3D parent. Populated by
  // PreprocessMesh; consumed by Solve for material-tensor rotation and iodata-path
  // projection. Mutable so the const PreprocessMesh hook can populate it; write-once per
  // solver instance.
  mutable std::unique_ptr<SubmeshFrame> frame_;
};

}  // namespace palace

#endif  // PALACE_DRIVERS_BOUNDARY_MODE_SOLVER_HPP
