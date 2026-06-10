// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_MODELS_MESH_SENSITIVITY_HPP
#define PALACE_MODELS_MESH_SENSITIVITY_HPP

#include <complex>
#include <string>
#include <vector>
#include <mfem.hpp>

namespace palace
{

class ComplexVector;
class IoData;
class SpaceOperator;

//
// Computes mesh sensitivities df/dx for mesh nodes.
//
// For the linear eigenvalue problem (lossless), analytic differentiation via Enzyme AD or
// forward-mode dual numbers produces exact mesh sensitivities in a single CeedOperator
// pass — same cost as one K*u evaluation. The differentiation mechanism (Enzyme vs dual
// numbers) is selected at compile time and encapsulated in the QFunction layer.
//
// For the quadratic EVP (lossy) and driven problems, a finite difference fallback is
// provided until the analytic path supports the damping term C.
//
class MeshSensitivity
{
public:
  // Result structure containing per-node sensitivities.
  struct SensitivityResult
  {
    // Number of boundary nodes (each has 3 spatial components).
    int num_bdr_nodes = 0;

    // The global vertex indices of the boundary nodes.
    std::vector<int> bdr_node_indices;

    // Sensitivity dJ/dx as a flat vector of size 3 * num_bdr_nodes.
    // Layout: [dx_0, dy_0, dz_0, dx_1, dy_1, dz_1, ...].
    mfem::Vector dJ_dx;
  };

  // Compute mesh sensitivities via finite differences for the eigenvalue problem.
  //
  // The objective is the eigenvalue sensitivity bilinear form:
  //   J(x) = u^H * (K(x) - omega^2 * M(x)) * u   (lossless linear EVP)
  //   J(x) = u^H * (K(x) + lambda*C(x) + lambda^2*M(x)) * u   (quadratic EVP)
  //
  // For the linear EVP with M-normalized eigenvectors, the eigenvalue sensitivity is:
  //   dmu/dx_i = u^H * (dK/dx_i - mu * dM/dx_i) * u
  //
  // Parameters:
  //   iodata      - configuration data
  //   space_op    - spatial discretization operator
  //   u           - eigenvector (forward solution)
  //   omega       - eigenfrequency (complex, from eigenvalue solver)
  //   has_damping - whether the problem is quadratic (C present)
  //   perturbation - FD step size; if <= 0, a default is computed from the mesh
  static SensitivityResult ComputeEigen(const IoData &iodata, SpaceOperator &space_op,
                                        const ComplexVector &u, std::complex<double> omega,
                                        bool has_damping, double perturbation);

  // Compute mesh sensitivities via finite differences for the driven (S-parameter) problem.
  //
  // The objective is the S-parameter adjoint bilinear form:
  //   J(x) = v^H * A(x) * u  where A = K + i*omega*C - omega^2*M
  //
  // Parameters:
  //   iodata      - configuration data
  //   space_op    - spatial discretization operator
  //   u           - forward solution (E-field)
  //   v           - adjoint solution
  //   omega       - angular frequency (real for driven problems)
  //   perturbation - FD step size; if <= 0, a default is computed from the mesh
  static SensitivityResult ComputeDriven(const IoData &iodata, SpaceOperator &space_op,
                                         const ComplexVector &u, const ComplexVector &v,
                                         double omega, double perturbation);

  // Collect the unique set of vertex indices on the mesh boundary. Works in parallel:
  // each process returns only its locally-owned boundary vertices.
  static std::vector<int> GetBoundaryVertexIndices(const mfem::ParMesh &mesh);

  // Compute the sesquilinear form v^H * A(x) * u for the current mesh configuration.
  // Resets geometry cache and reassembles K, C, M from scratch using SpaceOperator.
  static std::complex<double> EvalSystemForm(SpaceOperator &space_op,
                                             const ComplexVector &u, const ComplexVector &v,
                                             std::complex<double> omega, bool has_damping,
                                             MPI_Comm comm);

  // Compute a default FD step size from the mesh bounding box diagonal.
  static double DefaultPerturbation(const mfem::ParMesh &mesh, MPI_Comm comm);

  // Write the mesh sensitivity dJ/dx as a .solb vector field file.
  //
  // The SensitivityResult contains dJ/dx only for boundary nodes. This method expands
  // the result to ALL mesh vertices (interior vertices get zero sensitivity), producing
  // a .solb file compatible with the .meshb file written by metric::WriteMeshMeshb.
  //
  // The expansion is necessary for the sensitivity chain: the .solb vertex count must
  // match the mesh vertex count so that dJ/dx can be contracted with dx/dp from the
  // geometry (EGADS) side, which is also defined on all mesh vertices.
  //
  // Only MPI rank 0 writes; the method gathers the parallel mesh to rank 0 and maps
  // local boundary indices to serial vertex indices.
  //
  // Parameters:
  //   filename - output .solb file path
  //   result   - sensitivity result from ComputeEigen or ComputeDriven
  //   mesh     - the parallel mesh (used to determine total vertex count and gather)
  //   comm     - MPI communicator
  static void WriteSensitivitySolb(const std::string &filename,
                                   const SensitivityResult &result, mfem::ParMesh &mesh,
                                   MPI_Comm comm);
};

}  // namespace palace

#endif  // PALACE_MODELS_MESH_SENSITIVITY_HPP
