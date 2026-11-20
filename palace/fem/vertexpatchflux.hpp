// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_FEM_VERTEX_PATCH_FLUX_HPP
#define PALACE_FEM_VERTEX_PATCH_FLUX_HPP

#include <vector>
#include "linalg/vector.hpp"

namespace palace
{

class FiniteElementSpace;
class MaterialOperator;

namespace fem
{

// Lowest-order equilibrated flux of div(ε ∇u) = 0 from vertex patches (Braess and
// Schöberl, Ern and Vohralík). For each potential v (true DOFs of primal) with flux
// D_v = -ε ∇v, and each hat ψₐ of the lowest-order H1 space hats (on a nonconforming mesh
// its support also contains the cells of the hanging vertices it constrains), σₐ minimizes
// ‖ε^{-1/2} (σ - ψₐ D_v)‖ over the fields of the lowest-order Raviart-Thomas space rt0
// (broken across the interior Dirichlet boundaries, or not) which are supported in the
// support of ψₐ and have a zero normal trace on the faces on which ψₐ vanishes and on the
// Neumann faces (exterior faces without a Dirichlet attribute), subject to ∫_K div σ = ∫_K
// ∇ψₐ ⋅ D_v on each cell K. The hats sum to one, so D = Σₐ σₐ is divergence free with a
// zero normal trace on the Neumann faces. When a patch has no face with a free normal trace
// (a free hat), its cell data sum to zero by Galerkin orthogonality, up to the algebraic
// error of v: the integrals use the quadrature rule of the operator assembly (fem::
// DefaultIntegrationOrder), in which v is the Galerkin solution, and the patch problem is
// solved in the least-squares sense, which subtracts the mean of the data; this residual
// remains as the divergence of D (which the estimator accounts for). The patch problems are
// solved on the process owning the hat, from the cell data of its support sent there.
// Returns the true-DOF vectors of rt0. Collective.
std::vector<Vector> VertexPatchFlux(const FiniteElementSpace &primal,
                                    const FiniteElementSpace &hats,
                                    const FiniteElementSpace &rt0,
                                    const MaterialOperator &mat_op,
                                    const std::vector<int> &dirichlet_attr,
                                    const std::vector<Vector> &potentials);

}  // namespace fem

}  // namespace palace

#endif  // PALACE_FEM_VERTEX_PATCH_FLUX_HPP
