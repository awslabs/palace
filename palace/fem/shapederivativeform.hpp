// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_FEM_SHAPE_DERIVATIVE_FORM_HPP
#define PALACE_FEM_SHAPE_DERIVATIVE_FORM_HPP

#include <memory>
#include <vector>
#include <mfem.hpp>
#include "fem/integrator.hpp"
#include "fem/libceed/operator.hpp"
#include "linalg/vector.hpp"

namespace palace
{

class FiniteElementSpace;

//
// The derivative of a form a(u, w) with respect to the positions of the mesh nodes: for
// fields u and w, the gradient g = ∂a(u, w)/∂X on the true DOFs of the mesh nodal space, so
// that a displacement V of the nodes changes a(u, w) by g ⋅ V to first order.
//
class ShapeDerivativeOperator
{
  friend class ShapeDerivativeForm;

private:
  // Space of the fields u and w, and space of the mesh nodes.
  const FiniteElementSpace &fespace;
  const mfem::ParFiniteElementSpace &nodal_fespace;

  // Composite operator of each thread, with the mesh nodal L-vector as output. The fields
  // are passive inputs: one pair of vectors per Ceed context, shared by its sub-operators.
  ceed::Operator op;
  mutable std::vector<CeedVector> u_vec, w_vec;
  mutable Vector u_l, w_l, g_l;

  ShapeDerivativeOperator(const FiniteElementSpace &fespace,
                          const mfem::ParFiniteElementSpace &nodal_fespace);

public:
  ~ShapeDerivativeOperator();

  const auto &GetNodalSpace() const { return nodal_fespace; }

  // Compute g = ∂a(u, w)/∂X for the true DOF vectors u and w.
  void Mult(const Vector &u, const Vector &w, Vector &g) const;
};

//
// Integrators of the derivative of a form with respect to the positions of the mesh nodes,
// one for each integrator of the form (see fem/integrator.hpp). The form must be on a
// scalar H1 or H(curl) space of a mesh with continuous nodes.
//
class ShapeDerivativeForm
{
private:
  const FiniteElementSpace &fespace;
  std::vector<std::unique_ptr<ShapeDerivativeIntegrator>> domain_integs, boundary_integs;

public:
  ShapeDerivativeForm(const FiniteElementSpace &fespace) : fespace(fespace) {}

  template <typename T, typename... U>
  void AddDomainIntegrator(U &&...args)
  {
    domain_integs.push_back(std::make_unique<T>(std::forward<U>(args)...));
  }

  template <typename T, typename... U>
  void AddBoundaryIntegrator(U &&...args)
  {
    boundary_integs.push_back(std::make_unique<T>(std::forward<U>(args)...));
  }

  std::unique_ptr<ShapeDerivativeOperator> Assemble() const;
};

}  // namespace palace

#endif  // PALACE_FEM_SHAPE_DERIVATIVE_FORM_HPP
