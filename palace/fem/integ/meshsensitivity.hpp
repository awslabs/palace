// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_FEM_INTEG_MESHSENSITIVITY_HPP
#define PALACE_FEM_INTEG_MESHSENSITIVITY_HPP

#include <memory>
#include <vector>
#include <mfem.hpp>
#include "fem/libceed/ceed.hpp"

namespace palace
{

class FiniteElementSpace;
class MaterialPropertyCoefficient;
class Mesh;

namespace ceed
{
class Operator;
}

// Assembles a CeedOperator that computes the mesh sensitivity df/dx_nodes for the combined
// CurlCurlMass eigenvalue integrand (f_K - lambda * f_M) in a single pass.
//
// The operator maps mesh node coordinates to per-node sensitivities:
//   sensitivity = B_mesh^T * QFunction(B_mesh * x_nodes, B_hcurl * u, attr, w)
//
// where the QFunction differentiates the integrand w.r.t. the Jacobian J at each
// quadrature point via Enzyme AD or forward-mode dual numbers. The chain rule contraction
// with mesh shape function gradients is handled by libCEED via CEED_EVAL_GRAD.
//
// This is NOT a BilinearFormIntegrator — it is a custom nonlinear functional with a
// different field layout (active: mesh nodes; passive: eigenvector, attributes).
class MeshSensitivityIntegrator
{
public:
  // Assemble and apply the sensitivity CeedOperator.
  //
  // Parameters:
  //   mesh       - Palace mesh (provides geometry data, attribute mappings)
  //   nd_fespace - ND (H(curl)) finite element space for the eigenvector
  //   Q_mass     - Mass coefficient (first in the coefficient pair)
  //   Q_stiff    - Stiffness coefficient (second in the coefficient pair, mu_inv)
  //   lambda     - Eigenvalue (omega^2 for linear EVP)
  //   eigenvector - Solution vector in the ND FE space (L-vector / true-dof vector)
  //   sensitivity - Output: per-node sensitivity in the mesh nodal FE space (same size as
  //                 mesh nodes vector). The caller must handle parallel reduction (P^T).
  static void Apply(Mesh &mesh, const FiniteElementSpace &nd_fespace,
                    const MaterialPropertyCoefficient *Q_mass,
                    const MaterialPropertyCoefficient *Q_stiff, double lambda,
                    const mfem::Vector &eigenvector, mfem::Vector &sensitivity);
};

}  // namespace palace

#endif  // PALACE_FEM_INTEG_MESHSENSITIVITY_HPP
