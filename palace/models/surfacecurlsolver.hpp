// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_DRIVERS_SURFACE_CURL_SOLVER_HPP
#define PALACE_DRIVERS_SURFACE_CURL_SOLVER_HPP

#include <vector>
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

// Build the curl-free cut cohomology generator (±Φ across a cut through the hole centroid)
// used as the London shifted-penalty drive a_h, on the 3D ND space (see
// BuildCutCohomologyGenerator).
Vector SolveSurfaceCurlProblem(const SurfaceFluxData &flux_data, const Mesh &mesh,
                               const FiniteElementSpace &nd_fespace,
                               PostOperator<ProblemType::MAGNETOSTATIC> &post_op);

void SolveSurfaceCurlProblem(const SurfaceFluxData &flux_data, const Mesh &mesh,
                             const FiniteElementSpace &nd_fespace,
                             PostOperator<ProblemType::MAGNETOSTATIC> &post_op,
                             Vector &result);

// Integrate the magnetic flux density B ⋅ n over the given boundary attributes, with the
// surface normal oriented according to flux_direction. Returns the global (summed) flux.
double ComputeFluxThroughSurface(const mfem::ParGridFunction &B_gf,
                                 const std::vector<int> &attributes, const Mesh &mesh,
                                 const MaterialOperator &mat_op,
                                 const mfem::Vector &flux_direction, MPI_Comm comm);

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
