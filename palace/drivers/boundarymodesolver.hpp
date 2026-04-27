// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_DRIVERS_BOUNDARY_MODE_SOLVER_HPP
#define PALACE_DRIVERS_BOUNDARY_MODE_SOLVER_HPP

#include <memory>
#include "drivers/basesolver.hpp"

namespace palace
{

struct SubmeshFrame;

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
  // Out-of-line ctor/dtor so that instantiations in main.cpp — where only the forward-
  // declared SubmeshFrame is visible — don't need the full type for default_delete.
  BoundaryModeSolver(const IoData &iodata, bool root, int size = 0, int num_thread = 0,
                     const char *git_tag = nullptr);
  ~BoundaryModeSolver() override;

  void PreprocessMesh(std::unique_ptr<mfem::Mesh> &smesh, MPI_Comm comm) const override;

  std::pair<ErrorIndicator, long long int>
  Solve(const std::vector<std::unique_ptr<Mesh>> &mesh) const override;

private:
  // Non-null when the solve mesh was extracted from a 3D parent. Populated by
  // PreprocessMesh; consumed by Solve. Mutable so the const PreprocessMesh hook can
  // populate it; the frame is itself write-once per solver instance.
  mutable std::unique_ptr<SubmeshFrame> frame_;
};

}  // namespace palace

#endif  // PALACE_DRIVERS_BOUNDARY_MODE_SOLVER_HPP
