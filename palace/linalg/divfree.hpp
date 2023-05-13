// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_LINALG_DIV_FREE_HPP
#define PALACE_LINALG_DIV_FREE_HPP

#include <memory>
#include <vector>
#include <mfem.hpp>
#include "linalg/complex.hpp"
#include "linalg/ksp.hpp"
#include "linalg/operator.hpp"
#include "linalg/vector.hpp"

namespace palace
{

class MaterialOperator;
class KspSolver;

//
// This solver implements a projection onto a divergence-free space satisfying Gᵀ M x = 0,
// where G represents the discrete gradient matrix with columns spanning the nullspace of
// the curl-curl operator.
//
class DivFreeSolver : public mfem::Solver
{
private:
  // Operators for the divergence-free projection.
  std::unique_ptr<ParOperator> WeakDiv, Grad;
  std::vector<std::unique_ptr<ParOperator>> M;

  // Linear solver for the projected linear system (Gᵀ M G) y = x.
  std::unique_ptr<KspSolver> ksp;

  // Workspace objects for solver application.
  mutable Vector psi, rhs;

public:
  DivFreeSolver(const MaterialOperator &mat_op, mfem::ParFiniteElementSpace &nd_fespace,
                mfem::ParFiniteElementSpaceHierarchy &h1_fespaces,
                const std::vector<mfem::Array<int>> &h1_bdr_tdof_lists, double tol,
                int max_it, int print);

  void SetOperator(const Operator &op) override {}

  // Given a vector of Nedelec dofs for an arbitrary vector field, compute the Nedelec dofs
  // of the irrotational portion of this vector field. The resulting vector will satisfy
  // ∇ x y = 0.
  void Mult(Vector &y) const
  {
    // Compute the divergence of y.
    WeakDiv->Mult(y, rhs);

    // Apply essential BC and solve the linear system.
    if (M.back()->GetEssentialTrueDofs())
    {
      rhs.SetSubVector(*M.back()->GetEssentialTrueDofs(), 0.0);
    }
    ksp->Mult(rhs, psi);

    // Compute the irrotational portion of y and subtract.
    Grad->AddMult(psi, y, 1.0);
  }
  void Mult(const Vector &x, Vector &y) const override
  {
    y = x;
    Mult(y);
  }
  void Mult(ComplexVector &y) const
  {
    Mult(y.Real());
    Mult(y.Imag());
  }
  void Mult(const ComplexVector &x, ComplexVector &y) const
  {
    y = x;
    Mult(y);
  }
};

}  // namespace palace

#endif  // PALACE_LINALG_DIV_FREE_HPP
