// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_MODELS_SURFACE_CURL_SOLVER_HPP
#define PALACE_MODELS_SURFACE_CURL_SOLVER_HPP

#include <vector>
#include "linalg/vector.hpp"
#include "utils/labels.hpp"

namespace palace
{

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

// The flux of ComputeFluxThroughSurface as a linear functional of B: f, on the true DOFs of
// the B space, with f^T B the flux (B evaluated in the first neighboring element, which the
// normal continuity of B makes equal to the two-sided average).
Vector FluxThroughSurfaceFunctional(const FiniteElementSpace &rt_fespace,
                                    const std::vector<int> &attributes,
                                    const mfem::Vector &flux_direction);

}  // namespace palace

#endif  // PALACE_MODELS_SURFACE_CURL_SOLVER_HPP
