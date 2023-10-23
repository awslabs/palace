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
class FluxProjector
{
private:
  // Operator for the mass matrix inversion.
  std::unique_ptr<Operator> Flux, M;

  // Linear solver and preconditioner for the projected linear system.
  std::unique_ptr<KspSolver> ksp;

  // Workspace object for solver application.
  mutable Vector rhs;

public:
  FluxProjector(const MaterialOperator &mat_op,
                const FiniteElementSpaceHierarchy &nd_fespaces, double tol, int max_it,
                int print, int pa_order_threshold);
  FluxProjector(const MaterialOperator &mat_op,
                const FiniteElementSpaceHierarchy &h1_fespaces,
                const FiniteElementSpace &h1d_fespace, double tol, int max_it, int print,
                int pa_order_threshold);

  template <typename VecType>
  void Mult(const VecType &x, VecType &y) const;
};

// Class used for computing curl flux error estimate, i.e. || μ⁻¹ ∇ × Uₕ - F ||_K where F
// denotes a smooth reconstruction of μ⁻¹ ∇ × Uₕ.
template <typename VecType>
class CurlFluxErrorEstimator
{
  using GridFunctionType =
      typename std::conditional<std::is_same<VecType, ComplexVector>::value,
                                mfem::ParComplexGridFunction, mfem::ParGridFunction>::type;

  // Reference to input data (not owned).
  const MaterialOperator &mat_op;

  // Finite element space used to represent U and F.
  const FiniteElementSpace &nd_fespace;

  // Global L2 projection solver.
  FluxProjector projector;

  // Temporary vectors for error estimation.
  mutable VecType F;
  mutable GridFunctionType F_gf, U_gf;

public:
  CurlFluxErrorEstimator(const MaterialOperator &mat_op,
                         const FiniteElementSpaceHierarchy &nd_fespaces, double tol,
                         int max_it, int print, int pa_order_threshold);

  // Compute elemental error indicators given a vector of true DOF.
  ErrorIndicator ComputeIndicators(const VecType &U) const;

  // Compute elemental error indicators given a vector of true DOF and fold into an existing
  // indicator.
  void AddErrorIndicator(const VecType &U, ErrorIndicator &indicator) const
  {
    indicator.AddIndicator(ComputeIndicators(U));
  }
};

// Class used for computing gradient flux error estimate, i.e. || ε ∇Uₕ - F ||_K, where F
// denotes a smooth reconstruction of ε ∇Uₕ.
class GradFluxErrorEstimator
{
  // Reference to input data (not owned).
  const MaterialOperator &mat_op;

  // Finite element space used to represent U.
  const FiniteElementSpace &h1_fespace;

  // Vector H1 space used to represent the components of F, ordered by component.
  std::unique_ptr<FiniteElementSpace> h1d_fespace;

  // Global L2 projection solver.
  FluxProjector projector;

  // Temporary vectors for error estimation.
  mutable Vector F;
  mutable mfem::ParGridFunction F_gf, U_gf;

public:
  GradFluxErrorEstimator(const MaterialOperator &mat_op,
                         const FiniteElementSpaceHierarchy &h1_fespaces, double tol,
                         int max_it, int print, int pa_order_threshold);

  // Compute elemental error indicators given a vector of true DOF.
  ErrorIndicator ComputeIndicators(const Vector &U) const;

  // Compute elemental error indicators given a vector of true DOF and fold into an existing
  // indicator.
  void AddErrorIndicator(const Vector &U, ErrorIndicator &indicator) const
  {
    indicator.AddIndicator(ComputeIndicators(U));
  }
};

}  // namespace palace

#endif  // PALACE_LINALG_ERROR_ESTIMATOR_HPP
