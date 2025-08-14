// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_DRIVERS_SURFACE_CURL_SOLVER_HPP
#define PALACE_DRIVERS_SURFACE_CURL_SOLVER_HPP

#include <memory>
#include "linalg/vector.hpp"

namespace palace
{

class IoData;
class Mesh;
class FiniteElementSpace;

// Solve 2D surface curl problem for flux loop initial condition
std::unique_ptr<Vector> SolveSurfaceCurlProblem(const IoData &iodata,
                                                const Mesh &mesh,
                                                const FiniteElementSpace &nd_fespace);

// Verify flux through holes using computed magnetic field B
void VerifyFluxThroughHoles(const mfem::ParGridFunction &B_gf,
                           const std::vector<int> &hole_attributes,
                           const std::vector<double> &target_fluxes,
                           const Mesh &mesh,
                           MPI_Comm comm);

}  // namespace palace

#endif  // PALACE_DRIVERS_SURFACE_CURL_SOLVER_HPP