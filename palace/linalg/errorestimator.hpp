// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_LINALG_ERROR_ESTIMATOR_HPP
#define PALACE_LINALG_ERROR_ESTIMATOR_HPP

#include <memory>
#include <mfem.hpp>
#include "fem/errorindicator.hpp"
#include "fem/fespace.hpp"
#include "fem/libceed/operator.hpp"
#include "linalg/ksp.hpp"
#include "linalg/operator.hpp"
#include "linalg/vector.hpp"

namespace palace
{

class MaterialPropertyCoefficient;
class MaterialOperator;

//
// Classes used in the estimation of element-wise solution errors via a global L2 projection
// of a discontinuous flux onto a smooth space (flux recovery).
//

template <typename VecType>
class TimeDependentFluxErrorEstimator;

// This solver computes a smooth recovery of a discontinuous flux. The difference between
// this resulting smooth flux and the original non-smooth flux provides a localizable error
// estimate.
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
  FluxProjector(const MaterialPropertyCoefficient &coeff,
                const FiniteElementSpaceHierarchy &smooth_fespaces,
                const FiniteElementSpace &rhs_fespace, double tol, int max_it, int print,
                bool use_mg);

  void Mult(const VecType &x, VecType &y) const;
};

// Class used for computing gradient flux error estimate, η_K = || ε Eₕ - D ||_K, where D
// denotes a smooth reconstruction of ε Eₕ = ε ∇Vₕ with continuous normal component.
template <typename VecType = Vector>
class GradFluxErrorEstimator
{
  friend class TimeDependentFluxErrorEstimator<VecType>;

private:
  // Finite element spaces used to represent E and the recovered D.
  const FiniteElementSpace &nd_fespace, &rt_fespace;

  // Global L2 projection solver.
  FluxProjector<VecType> projector;

  // Operator which performs the integration of the flux error on each element.
  ceed::Operator integ_op;

  // Temporary vectors for error estimation.
  mutable VecType E_gf, D, D_gf;

public:
  GradFluxErrorEstimator(const MaterialOperator &mat_op, FiniteElementSpace &nd_fespace,
                         FiniteElementSpaceHierarchy &rt_fespaces, double tol, int max_it,
                         int print, bool use_mg);

  // Compute elemental error indicators given the electric field as a vector of true dofs,
  // and fold into an existing indicator. The indicators are nondimensionalized using the
  // total field energy.
  void AddErrorIndicator(const VecType &E, double Et, ErrorIndicator &indicator) const;
};

// Class used for computing curl flux error estimate, η_K = || μ⁻¹ Bₕ - H ||_K where H
// denotes a smooth reconstruction of μ⁻¹ Bₕ = μ⁻¹ ∇ × Eₕ with continuous tangential
// component.
template <typename VecType = Vector>
class CurlFluxErrorEstimator
{
  friend class TimeDependentFluxErrorEstimator<VecType>;

private:
  // Finite element space used to represent B and the recovered H.
  const FiniteElementSpace &rt_fespace, &nd_fespace;

  // Global L2 projection solver.
  FluxProjector<VecType> projector;

  // Operator which performs the integration of the flux error on each element.
  ceed::Operator integ_op;

  // Temporary vectors for error estimation.
  mutable VecType B_gf, H, H_gf;

public:
  CurlFluxErrorEstimator(const MaterialOperator &mat_op, FiniteElementSpace &rt_fespace,
                         FiniteElementSpaceHierarchy &nd_fespaces, double tol, int max_it,
                         int print, bool use_mg);

  // Compute elemental error indicators given the magnetic flux density as a vector of true
  // dofs, and fold into an existing indicator. The indicators are nondimensionalized using
  // the total field energy.
  void AddErrorIndicator(const VecType &B, double Et, ErrorIndicator &indicator) const;
};

// Class used for computing sum of the gradient flux and curl flux error estimates,
// η²_K = || ε Eₕ - D ||²_K + || μ⁻¹ Bₕ - H ||²_K, where D and H denote a smooth
// reconstructions of ε Eₕ = ε ∇Vₕ with continuous normal component and μ⁻¹ Bₕ = μ⁻¹ ∇ × Eₕ
// with continuous tangential component.
template <typename VecType>
class TimeDependentFluxErrorEstimator
{
private:
  GradFluxErrorEstimator<VecType> grad_estimator;
  CurlFluxErrorEstimator<VecType> curl_estimator;

public:
  TimeDependentFluxErrorEstimator(const MaterialOperator &mat_op,
                                  FiniteElementSpaceHierarchy &nd_fespaces,
                                  FiniteElementSpaceHierarchy &rt_fespaces, double tol,
                                  int max_it, int print, bool use_mg);

  // Compute elemental error indicators given the electric field and magnetic flux density
  // as a vectors of true dofs, and fold into an existing indicator. The indicators are
  // nondimensionalized using the total field energy.
  void AddErrorIndicator(const VecType &E, const VecType &B, double Et,
                         ErrorIndicator &indicator) const;
};

}  // namespace palace

#endif  // PALACE_LINALG_ERROR_ESTIMATOR_HPP
