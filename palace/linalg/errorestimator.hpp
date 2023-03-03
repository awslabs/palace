// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_ERROR_ESTIMATOR_HPP
#define PALACE_ERROR_ESTIMATOR_HPP

#include <memory>
#include <mfem.hpp>
#include "linalg/petsc.hpp"

namespace palace
{

class MaterialOperator;

//
// This solver implements a solver to compute a smooth reconstruction of a
// discontinuous flux. The difference between this resulting smooth flux and the
// original non-smooth flux provides a localizable error estimate. An instance
// of FluxProjector can be reused across solutions, thus the construction of the
// operator is separated from the construction of the flux RHS.
class FluxProjector : public mfem::Solver
{
private:
  // Operator for the mass matrix inversion
  std::vector<std::unique_ptr<mfem::Operator>> M;

  // Linear solver and preconditioner for the projected linear system M σ = σ̂
  std::unique_ptr<mfem::IterativeSolver> ksp;
  std::unique_ptr<mfem::Solver> pc;

  mutable mfem::Vector tmp;

public:
  FluxProjector(mfem::ParFiniteElementSpaceHierarchy &smooth_flux_fes, double tol,
                int max_it, int print);

  // Operator is set in constructor.
  void SetOperator(const mfem::Operator &op) override {}

  // Given a vector of dof defining the flux
  void Mult(mfem::Vector &x) const
  {
    tmp = x;
    Mult(tmp, x);
  }
  void Mult(const mfem::Vector &x, mfem::Vector &y) const override { ksp->Mult(x, y); }
  void Mult(petsc::PetscParVector &x) const
  {
    auto cv = x.GetToVectors();
    tmp = cv.real;
    Mult(tmp, cv.real);
    tmp = cv.imag;
    Mult(tmp, cv.imag);
    x.SetFromVectors(cv);
  }
  void Mult(const petsc::PetscParVector &x, petsc::PetscParVector &y) const
  {
    auto cv = x.GetToVectors();
    tmp = cv.real;
    Mult(tmp, cv.real);
    tmp = cv.imag;
    Mult(tmp, cv.imag);
    y.SetFromVectors(cv);
  }

  // TODO: Add in an ArrayMult optimization

  petsc::PetscParVector GetFluxRHS(const mfem::ParGridFunction &phi) const;
};

}  // namespace palace

#endif  // PALACE_ERROR_ESTIMATOR_HPP
