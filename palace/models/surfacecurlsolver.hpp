// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_DRIVERS_SURFACE_CURL_SOLVER_HPP
#define PALACE_DRIVERS_SURFACE_CURL_SOLVER_HPP

#include <memory>
#include "linalg/vector.hpp"
#include "utils/labels.hpp"

namespace palace
{

class IoData;
class Mesh;
class FiniteElementSpace;
class MaterialOperator;

// Forward declarations
class SurfaceFluxData;
class CurlCurlOperator;
template <ProblemType T>
class PostOperator;

// Solve 2D surface curl problem for flux loop initial condition
Vector SolveSurfaceCurlProblem(const SurfaceFluxData &flux_data, const IoData &iodata,
                               const Mesh &mesh, const FiniteElementSpace &nd_fespace,
                               int flux_loop_idx,
                               PostOperator<ProblemType::MAGNETOSTATIC> &post_op);

void SolveSurfaceCurlProblem(const SurfaceFluxData &flux_data, const IoData &iodata,
                             const Mesh &mesh, const FiniteElementSpace &nd_fespace,
                             int flux_loop_idx,
                             PostOperator<ProblemType::MAGNETOSTATIC> &post_op,
                             Vector &result);

void VerifyFluxThroughHoles(const mfem::ParGridFunction &B_gf,
                            const std::vector<int> &hole_attributes,
                            const std::vector<double> &target_fluxes, const Mesh &mesh,
                            const MaterialOperator &mat_op,
                            const mfem::Vector &flux_direction, MPI_Comm comm);

// Verify flux through all holes in a multi flux setting
void VerifyFluxThroughAllHoles(const mfem::ParGridFunction &B_gf, const IoData &iodata,
                               int current_flux_loop_idx, const Mesh &mesh,
                               const MaterialOperator &mat_op, MPI_Comm comm);

}  // namespace palace

#endif  // PALACE_DRIVERS_SURFACE_CURL_SOLVER_HPP