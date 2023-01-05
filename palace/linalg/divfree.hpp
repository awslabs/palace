// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_DIV_FREE_HPP
#define PALACE_DIV_FREE_HPP

#include <memory>
#include <mfem.hpp>
#include "linalg/petsc.hpp"

namespace palace
{

class MaterialOperator;

//
// This solver implements a projection onto a divergence-free space satisfying Gᵀ M x = 0,
// where G represents the discrete gradient matrix with columns spanning the nullspace of
// the curl-curl operator.
//
class DivFreeSolver : public mfem::Solver
{
private:
  // Operators for the divergence-free projection.
  std::unique_ptr<mfem::Operator> WeakDiv, Grad;
  std::vector<std::unique_ptr<mfem::Operator>> M;

  // Linear solver and preconditioner for the projected linear system (Gᵀ M G) y = x.
  std::unique_ptr<mfem::IterativeSolver> ksp;
  std::unique_ptr<mfem::Solver> pc;

  // Workspace objects for solver application.
  mutable mfem::Vector psi, rhs, xr, xi;

  // Boundary condition dofs for essential BCs.
  mfem::Array<int> h1_bdr_tdof_list;

public:
  DivFreeSolver(const MaterialOperator &mat_op, const mfem::Array<int> &bdr_marker,
                mfem::ParFiniteElementSpace &nd_fespace,
                mfem::ParFiniteElementSpaceHierarchy &h1_fespaces, double tol, int max_it,
                int print);

  // Operator is set in constructor.
  void SetOperator(const mfem::Operator &op) override {}

  // Given a vector of Nedelec dofs for an arbitrary vector field, compute the Nedelec dofs
  // of the irrotational portion of this vector field. The resulting vector will satisfy
  // ∇ x x = 0.
  void Mult(mfem::Vector &x) const
  {
    // Compute the divergence of x.
    WeakDiv->Mult(x, rhs);

    // Apply essential BC and solve the linear system.
    psi = 0.0;
    rhs.SetSubVector(h1_bdr_tdof_list, 0.0);
    ksp->Mult(rhs, psi);

    // Compute the irrotational portion of x and subtract.
    Grad->AddMult(psi, x, 1.0);
  }
  void Mult(const mfem::Vector &x, mfem::Vector &y) const override
  {
    y = x;
    Mult(y);
  }
  void Mult(petsc::PetscParVector &x) const
  {
    x.GetToVectors(xr, xi);
    Mult(xr);
    Mult(xi);
    x.SetFromVectors(xr, xi);
  }
  void Mult(const petsc::PetscParVector &x, petsc::PetscParVector &y) const
  {
    y.Copy(x);
    Mult(y);
  }
  using mfem::Operator::Mult;
};

}  // namespace palace

#endif  // PALACE_DIV_FREE_HPP
