// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_FEM_METRIC_FIELD_HPP
#define PALACE_FEM_METRIC_FIELD_HPP

#include <string>
#include <mfem.hpp>
#include "fem/errorindicator.hpp"

namespace palace
{

//
// Utilities for converting per-element error indicators into per-vertex metric tensors
// and writing the result in libMeshb .solb/.meshb format. The metric tensor encodes
// desired element sizes for external mesh adaptation tools (Metris, refine).
//
namespace metric
{

// Convert per-element error indicators to per-vertex isotropic metric tensors.
//
// Algorithm (equidistribution principle):
//   1. Average element errors to vertices, weighted by element volume.
//   2. Compute desired element size via h_k = h_avg * (eta_k / eta_avg)^{-2/(2p+d)}
//      where p = polynomial order, d = spatial dimension.
//   3. Apply h_min/h_max bounds to prevent degenerate elements.
//   4. Build symmetric metric tensor M = (1/h_k)^2 * I per vertex.
//
// Returns a vector of metric tensor components laid out per vertex:
//   - 2D: 3 components per vertex (M11, M12, M22)
//   - 3D: 6 components per vertex (M11, M12, M13, M22, M23, M33)
mfem::Vector ComputeIsotropicMetric(mfem::ParMesh &mesh, const ErrorIndicator &indicator,
                                    int order, double h_min = 0.0, double h_max = 0.0);

// Write metric tensor field as a .solb ASCII file (libMeshb SolAtVertices format).
// The metric vector should contain ncomp entries per vertex in upper-triangle order
// (3 for 2D, 6 for 3D).
void WriteSolb(const std::string &filename, const mfem::Vector &metric, int num_vertices,
               int dim);

// Write the mesh as an ASCII .meshb file (libMeshb format) to accompany the .solb.
// Only the root MPI process writes; others return immediately.
void WriteMeshMeshb(const std::string &filename, mfem::ParMesh &mesh);

}  // namespace metric

}  // namespace palace

#endif  // PALACE_FEM_METRIC_FIELD_HPP
