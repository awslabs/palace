// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_LINALG_ERROR_ESTIMATOR_HPP
#define PALACE_LINALG_ERROR_ESTIMATOR_HPP

#include <memory>
#include <vector>
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

// The flux recovery space of the error estimators: a view of the given space which is
// broken (discontinuous) across the interior boundaries formed by the boundary attributes
// in crack_attr_list (see FiniteElementSpace's broken-view constructor), or nullptr when
// the mesh has no such interior boundaries (then the recovery uses the given space itself).
// The estimators and the postprocessing of the recovered flux (PostOperator) build it with
// this one rule from the same space and list, so that they hold D on the same space
// (SD-403, SD-3202). Collective.
std::unique_ptr<FiniteElementSpace>
BuildBrokenFluxSpace(FiniteElementSpace &fespace, const std::vector<int> &crack_attr_list);

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
  using OperType =
      std::conditional_t<std::is_same_v<VecType, ComplexVector>, ComplexOperator, Operator>;

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

  // Recovery CG uses the preconditioned residual norm, not the Euclidean norm.
  bool GetConverged() const { return ksp->GetConverged(); }
  int GetTotalIterations() const { return ksp->NumTotalMultIterations(); }
  double GetFinalRelativeResidual() const { return ksp->GetFinalRelativeResidual(); }
  double GetInitialResidual() const { return ksp->GetInitialResidual(); }
  double GetFinalResidual() const { return ksp->GetFinalResidual(); }
};

// Forward declarations for friend access.
template <typename VecType>
class BoundaryModeFluxErrorEstimator;

// Class used for computing gradient flux error estimate, η_K = || ε Eₕ - D ||_K, where D
// denotes a smooth reconstruction of ε Eₕ = ε ∇Vₕ with continuous normal component.
template <typename VecType = Vector>
class GradFluxErrorEstimator
{
  friend class TimeDependentFluxErrorEstimator<VecType>;
  friend class BoundaryModeFluxErrorEstimator<VecType>;

private:
  // Recovery space which is discontinuous across interior boundaries, when there are any.
  std::unique_ptr<FiniteElementSpaceHierarchy> broken_rt_fespaces;

  // Finite element spaces used to represent E and the recovered D.
  const FiniteElementSpace &nd_fespace, &rt_fespace;

  // Global L2 projection solver.
  FluxProjector<VecType> projector;

  // Operator which performs the integration of the flux error on each element.
  ceed::Operator integ_op;

  // Temporary vectors for error estimation.
  mutable VecType E_gf, D, D_gf;

public:
  // The recovered D is allowed to be discontinuous across the interior boundaries formed
  // by the boundary attributes in crack_attr_list (which are treated as physically
  // separating the two sides, e.g. thin conducting sheets), using a broken finite element
  // space as for a mesh cut along these interior boundaries.
  GradFluxErrorEstimator(const MaterialOperator &mat_op, FiniteElementSpace &nd_fespace,
                         FiniteElementSpaceHierarchy &rt_fespaces, double tol, int max_it,
                         int print, bool use_mg,
                         const std::vector<int> &crack_attr_list = {});

  // Compute elemental error indicators given the electric field as a vector of true dofs,
  // and fold into an existing indicator. The indicators are nondimensionalized using the
  // total field energy.
  void AddErrorIndicator(const VecType &E, double Et, ErrorIndicator &indicator) const;
  void AddErrorIndicator(const VecType &E, const VecType &D, double Et,
                         ErrorIndicator &indicator) const;

  // Recover D = epsilon E in the H(div)-conforming flux space.
  void RecoverFlux(const VecType &E, VecType &D) const;

  // The space of the recovered D (true dofs): broken across the interior boundaries of
  // crack_attr_list when the mesh has any, else the finest space of rt_fespaces.
  const FiniteElementSpace &GetFluxSpace() const { return rt_fespace; }

  const FluxProjector<VecType> &GetFluxProjector() const { return projector; }
};

// Class used for computing curl flux error estimate, η_K = || μ⁻¹ Bₕ - H ||_K where H
// denotes a smooth reconstruction of μ⁻¹ Bₕ = μ⁻¹ ∇ × Eₕ with continuous tangential
// component.
template <typename VecType = Vector>
class CurlFluxErrorEstimator
{
  friend class TimeDependentFluxErrorEstimator<VecType>;
  friend class BoundaryModeFluxErrorEstimator<VecType>;

private:
  // Recovery space which is discontinuous across interior boundaries, when there are any.
  std::unique_ptr<FiniteElementSpaceHierarchy> broken_nd_fespaces;

  // Finite element space used to represent B and the recovered H.
  const FiniteElementSpace &rt_fespace, &nd_fespace;

  // Global L2 projection solver.
  FluxProjector<VecType> projector;

  // Operator which performs the integration of the flux error on each element.
  ceed::Operator integ_op;

  // Temporary vectors for error estimation.
  mutable VecType B_gf, H, H_gf;

public:
  // The recovered H is allowed to be discontinuous across the interior boundaries formed
  // by the boundary attributes in crack_attr_list, see GradFluxErrorEstimator.
  CurlFluxErrorEstimator(const MaterialOperator &mat_op, FiniteElementSpace &rt_fespace,
                         FiniteElementSpaceHierarchy &nd_fespaces, double tol, int max_it,
                         int print, bool use_mg,
                         const std::vector<int> &crack_attr_list = {});

  // Compute elemental error indicators given the magnetic flux density as a vector of true
  // dofs, and fold into an existing indicator. The indicators are nondimensionalized using
  // the total field energy.
  void AddErrorIndicator(const VecType &B, double Et, ErrorIndicator &indicator) const;
};

// Class used for computing sum of the gradient flux and curl flux error estimates for 3D,
// η²_K = || ε Eₕ - D ||²_K + || μ⁻¹ Bₕ - H ||²_K, where D and H denote smooth
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
                                  int max_it, int print, bool use_mg,
                                  const std::vector<int> &crack_attr_list = {});

  void AddErrorIndicator(const VecType &E, const VecType &B, double Et,
                         ErrorIndicator &indicator) const;
  void AddErrorIndicator(const VecType &E, const VecType &B, const VecType &D, double Et,
                         ErrorIndicator &indicator) const;
  void RecoverElectricFlux(const VecType &E, VecType &D) const;
  const FiniteElementSpace &GetElectricFluxSpace() const
  {
    return grad_estimator.GetFluxSpace();
  }
};

// 2D boundary mode error estimator. Owns its own single-level FE spaces for flux
// recovery. Uses L2 curl space (scalar Bz) + H1 hierarchy (smooth H projection)
// instead of RT/ND.
template <typename VecType>
class BoundaryModeFluxErrorEstimator
{
private:
  GradFluxErrorEstimator<VecType> grad_estimator;
  CurlFluxErrorEstimator<VecType> curl_estimator;

public:
  BoundaryModeFluxErrorEstimator(const MaterialOperator &mat_op,
                                 FiniteElementSpaceHierarchy &nd_fespaces,
                                 FiniteElementSpaceHierarchy &rt_fespaces,
                                 FiniteElementSpace &curl_fespace,
                                 FiniteElementSpaceHierarchy &h1_fespaces, double tol,
                                 int max_it, int print, bool use_mg,
                                 const std::vector<int> &crack_attr_list = {});

  void AddErrorIndicator(const VecType &E, const VecType &B, double Et,
                         ErrorIndicator &indicator) const;
  void AddErrorIndicator(const VecType &E, const VecType &B, const VecType &D, double Et,
                         ErrorIndicator &indicator) const;
  void RecoverElectricFlux(const VecType &E, VecType &D) const;
  const FiniteElementSpace &GetElectricFluxSpace() const
  {
    return grad_estimator.GetFluxSpace();
  }
};

}  // namespace palace

#endif  // PALACE_LINALG_ERROR_ESTIMATOR_HPP
