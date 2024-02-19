// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_LINALG_ERROR_ESTIMATOR_HPP
#define PALACE_LINALG_ERROR_ESTIMATOR_HPP

#include <memory>
#include <mfem.hpp>
#include "fem/errorindicator.hpp"
#include "fem/fespace.hpp"
#include "linalg/ksp.hpp"
#include "linalg/operator.hpp"
#include "linalg/vector.hpp"

namespace palace
{

class MaterialOperator;

//
// Classes used in the estimation of element-wise solution errors via a global L2 projection
// of a discontinuous flux onto a smooth space.
//

// This solver computes a smooth reconstruction of a discontinuous flux. The difference
// between this resulting smooth flux and the original non-smooth flux provides a
// localizable error estimate.
template <typename VecType>
class FluxProjector
{
  using OperType = typename std::conditional<std::is_same<VecType, ComplexVector>::value,
                                             ComplexOperator, Operator>::type;

private:
  // Operator for the mass matrix inversion.
  std::unique_ptr<OperType> Flux, M;

  // Linear solver and preconditioner for the projected linear system.
  std::unique_ptr<BaseKspSolver<OperType>> ksp;

  // Workspace object for solver application.
  mutable VecType rhs;

public:
  FluxProjector(const MaterialOperator &mat_op, const FiniteElementSpace &nd_fespace,
                double tol, int max_it, int print);
  FluxProjector(const MaterialOperator &mat_op, const FiniteElementSpace &h1_fespace,
                const FiniteElementSpace &h1d_fespace, double tol, int max_it, int print);

  void Mult(const VecType &x, VecType &y) const;
};

// Class used for computing curl flux error estimate, i.e. || μ⁻¹ ∇ × Uₕ - F ||_K where F
// denotes a smooth reconstruction of μ⁻¹ ∇ × Uₕ.
template <typename VecType>
class CurlFluxErrorEstimator
{
  // Reference to material property data (not owned).
  const MaterialOperator &mat_op;

  // Finite element space used to represent U and F.
  FiniteElementSpace &nd_fespace;

  // Global L2 projection solver.
  FluxProjector<VecType> projector;

  // Temporary vectors for error estimation.
  mutable VecType F, F_gf, U_gf;

public:
  CurlFluxErrorEstimator(const MaterialOperator &mat_op, FiniteElementSpace &nd_fespace,
                         double tol, int max_it, int print);

  // Compute elemental error indicators given a vector of true DOF and fold into an existing
  // indicator.
  void AddErrorIndicator(const VecType &U, ErrorIndicator &indicator) const;
};

// Class used for computing gradient flux error estimate, i.e. || ε ∇Uₕ - F ||_K, where F
// denotes a smooth reconstruction of ε ∇Uₕ.
class GradFluxErrorEstimator
{
  // Reference to material property data (not owned).
  const MaterialOperator &mat_op;

  // Finite element space used to represent U.
  FiniteElementSpace &h1_fespace;

  // Vector H1 space used to represent the components of F, ordered by component.
  std::unique_ptr<FiniteElementSpace> h1d_fespace;

  // Global L2 projection solver.
  FluxProjector<Vector> projector;

  // Temporary vectors for error estimation.
  mutable Vector F, F_gf, U_gf;

public:
  GradFluxErrorEstimator(const MaterialOperator &mat_op, FiniteElementSpace &h1_fespace,
                         double tol, int max_it, int print);

  // Compute elemental error indicators given a vector of true DOF and fold into an existing
  // indicator.
  void AddErrorIndicator(const Vector &U, ErrorIndicator &indicator) const;
};

}  // namespace palace

#endif  // PALACE_LINALG_ERROR_ESTIMATOR_HPP
