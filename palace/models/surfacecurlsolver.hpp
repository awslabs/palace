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

// Forward declaration
class SurfaceFluxData;

// Solve 2D surface curl problem for flux loop initial condition
Vector SolveSurfaceCurlProblem(const SurfaceFluxData &data,
                               const IoData &iodata, const Mesh &mesh,
                               const FiniteElementSpace &nd_fespace,
                               int flux_loop_idx);

// Overload with pre-allocated result vector
void SolveSurfaceCurlProblem(const SurfaceFluxData &data, const IoData &iodata,
                             const Mesh &mesh, const FiniteElementSpace &nd_fespace,
                             int flux_loop_idx, Vector &result);

// Verify flux through holes using computed magnetic field B
void VerifyFluxThroughHoles(const mfem::ParGridFunction &B_gf,
                            const std::vector<int> &hole_attributes,
                            const std::vector<double> &target_fluxes, const Mesh &mesh,
                            MPI_Comm comm);

// Verify flux through all holes in a multi flux setting
void VerifyFluxThroughAllHoles(const mfem::ParGridFunction &B_gf, const IoData &iodata,
                               int current_flux_loop_idx, const Mesh &mesh, MPI_Comm comm);

}  // namespace palace

#endif  // PALACE_DRIVERS_SURFACE_CURL_SOLVER_HPP