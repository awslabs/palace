// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_LINALG_FLUX_PROJECTOR_HPP
#define PALACE_LINALG_FLUX_PROJECTOR_HPP

#include <memory>
#include <mfem.hpp>
#include "linalg/ksp.hpp"
#include "linalg/operator.hpp"
#include "linalg/vector.hpp"

namespace palace
{

class MaterialOperator;

//
// This solver implements a solver to compute a smooth reconstruction of a
// discontinuous flux. The difference between this resulting smooth flux and the
// original non-smooth flux provides a localizable error estimate. An instance
// of FluxProjector can be reused across solutions, thus the construction of the
// operator is separated from the construction of the flux RHS.
class FluxProjector
{
private:
  // Operator for the mass matrix inversion
  std::unique_ptr<Operator> M;

  // Linear solver and preconditioner for the projected linear system M σ = σ̂
  std::unique_ptr<KspSolver> ksp;

  mutable Vector tmp;

public:
  FluxProjector(mfem::ParFiniteElementSpaceHierarchy &smooth_flux_fes, double tol = 1e-12,
                int max_it = 200, int print_level = 1, int pa_order_threshold = 1);

  // Given a vector of dof defining the flux, compute the smooth flux IN PLACE.
  void Mult(Vector &x) const
  {
    tmp = x;
    Mult(tmp, x);
  }
  void Mult(const Vector &x, Vector &y) const { ksp->Mult(x, y); }
  void Mult(ComplexVector &x) const
  {
    Mult(x.Real());
    Mult(x.Imag());
  }
  void Mult(const ComplexVector &x, ComplexVector &y) const
  {
    y = x;
    Mult(y);
  }
};

}  // namespace palace

#endif  // PALACE_LINALG_FLUX_PROJECTOR_HPP
