// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_LINALG_ERROR_ESTIMATOR_HPP
#define PALACE_LINALG_ERROR_ESTIMATOR_HPP

#include <complex>
#include <memory>
#include <utility>
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

//
// Classes used in the estimation of element-wise solution errors via a global L2 projection
// of a discontinuous flux onto a smooth space (flux recovery).
//

template <typename VecType>
class TimeDependentFluxErrorEstimator;

// This solver computes a smooth recovery of a discontinuous flux. The difference between
// this resulting smooth flux and the original non-smooth flux provides a localizable error
// estimate.
//
// The projection of a frequency-dependent flux, F(ω) = F₀ + Σₖ wₖ(ω) Fₖ, uses the
// frequency-independent mass matrix, solver, and preconditioner of the smooth space, with
// the constant flux F₀ and material-restricted flux operators Fₖ assembled once and only
// their complex scalar weights wₖ given at each frequency.
template <typename VecType>
class FluxProjector
{
  using OperType =
      std::conditional_t<std::is_same_v<VecType, ComplexVector>, ComplexOperator, Operator>;

private:
  // Operator for the mass matrix inversion.
  std::unique_ptr<OperType> Flux, M;

  // Flux operators Fₖ for the weighted terms.
  std::vector<std::unique_ptr<OperType>> Flux_k;

  // Linear solver and preconditioner for the projected linear system.
  std::unique_ptr<BaseKspSolver<OperType>> ksp;

  // Workspace object for solver application.
  mutable VecType rhs;

public:
  // The constant flux F₀ only, and with the weighted terms Fₖ. Not a default argument: the
  // coefficient type is incomplete here.
  FluxProjector(const MaterialPropertyCoefficient &coeff,
                const FiniteElementSpaceHierarchy &smooth_fespaces,
                const FiniteElementSpace &rhs_fespace, double tol, int max_it, int print,
                bool use_mg);
  FluxProjector(const MaterialPropertyCoefficient &coeff,
                const FiniteElementSpaceHierarchy &smooth_fespaces,
                const FiniteElementSpace &rhs_fespace, double tol, int max_it, int print,
                bool use_mg, const std::vector<MaterialPropertyCoefficient> &coeff_k);

  // y = M⁻¹ F₀ x.
  void Mult(const VecType &x, VecType &y) const;

  // y = M⁻¹ (F₀ + Σₖ wₖ Fₖ) x, with one weight per term.
  void Mult(const VecType &x, VecType &y, const std::vector<std::complex<double>> &w) const
    requires std::is_same_v<VecType, ComplexVector>;
};

// Retained libCEED QFunction contexts (one per element geometry sub-operator and libCEED
// context), owned for in-place updates of their values.
class CeedContextList
{
private:
  std::vector<std::pair<Ceed, CeedQFunctionContext>> ctx;

public:
  CeedContextList() = default;
  CeedContextList(CeedContextList &&other) noexcept : ctx(std::move(other.ctx))
  {
    other.ctx.clear();
  }
  CeedContextList &operator=(CeedContextList &&) = delete;
  CeedContextList(const CeedContextList &) = delete;
  CeedContextList &operator=(const CeedContextList &) = delete;
  ~CeedContextList();

  void Append(Ceed ceed, CeedQFunctionContext c);
  bool Empty() const { return ctx.empty(); }

  // Overwrite the values of every context, all of the given size.
  void Update(const void *data, std::size_t size) const;
};

// Forward declarations for friend access.
template <typename VecType>
class BoundaryModeFluxErrorEstimator;

// Class used for computing gradient flux error estimate, η_K = || ε Eₕ - D ||_K, where D
// denotes a smooth reconstruction of ε Eₕ = ε ∇Vₕ with continuous normal component, in
// the energy norm || F ||²_K = ∫_K Fᴴ ε⁻¹ F dx.
//
// With a frequency-dependent material (complex fields only), the discontinuous flux is
// the complete constitutive Dₕ = ε(iω) Eₕ, with the complex permittivity of every
// material, including the loss tangent of a constant one, and the recovered D its
// projection. The indicator is η²_K = ∫_K (Dₕ - D)ᴴ W⁻¹ (Dₕ - D) dx, with the positive
// stored-energy permittivity W(ω) of each model (Re{ε} for a constant material): a
// positive, energy-scaled recovery indicator, not a stored-energy error bound. W E is not a
// displacement (εE, not WE, has a continuous normal component), so it is not recovered.
// For a general PoleResidue model with W(ω) ≤ 0, the indicator falls back to ε∞ with a
// warning. Without a frequency-dependent material, the established indicator above with
// Re{ε} is unchanged.
template <typename VecType = Vector>
class GradFluxErrorEstimator
{
  friend class TimeDependentFluxErrorEstimator<VecType>;
  friend class BoundaryModeFluxErrorEstimator<VecType>;

private:
  // Finite element spaces used to represent E and the recovered D.
  const FiniteElementSpace &nd_fespace, &rt_fespace;

  // Material properties (not owned), and the globally supported frequency-dependent
  // materials (configured models), which select the complex constitutive route. With a
  // loss tangent, that route also projects Im{ε} of the constant materials.
  const MaterialOperator &mat_op;
  std::vector<std::size_t> fd_materials;
  bool fd_loss_tangent = false;

  // Global L2 projection solver.
  FluxProjector<VecType> projector;

  // Operator which performs the integration of the flux error on each element, and its
  // retained contexts (ε', ε'', W⁻¹ tables) on the complex constitutive route, at the
  // real angular frequency fd_omega.
  ceed::Operator integ_op;
  CeedContextList integ_ctx;
  mutable double fd_omega = 0.0;
  mutable std::vector<bool> fd_warned;

  // Temporary vectors for error estimation.
  mutable VecType E_gf, D, D_gf;

  // Squared elemental indicators, and the squared elemental energy norms of the recovered
  // flux of the last estimate (for a fallback normalization).
  Vector ComputeEstimates(const VecType &E) const;
  Vector ComputeEstimates(const VecType &E, double omega) const
    requires std::is_same_v<VecType, ComplexVector>;
  Vector RecoveredFluxNorms() const;
  void SetFrequency(double omega) const
    requires std::is_same_v<VecType, ComplexVector>;

public:
  GradFluxErrorEstimator(const MaterialOperator &mat_op, FiniteElementSpace &nd_fespace,
                         FiniteElementSpaceHierarchy &rt_fespaces, double tol, int max_it,
                         int print, bool use_mg);

  // Compute elemental error indicators given the electric field as a vector of true dofs,
  // and fold into an existing indicator. The indicators are nondimensionalized using the
  // total field energy. Not for a frequency-dependent material.
  void AddErrorIndicator(const VecType &E, double Et, ErrorIndicator &indicator) const;

  // The same at the real angular frequency ω of the field, for a frequency-dependent
  // material. A nonpositive or nonfinite Et (a general model with a negative stored
  // energy) uses the energy norm of the recovered flux instead, with a warning.
  void AddErrorIndicator(const VecType &E, double Et, double omega,
                         ErrorIndicator &indicator) const
    requires std::is_same_v<VecType, ComplexVector>;
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
  // Finite element space used to represent B and the recovered H.
  const FiniteElementSpace &rt_fespace, &nd_fespace;

  // Global L2 projection solver.
  FluxProjector<VecType> projector;

  // Operator which performs the integration of the flux error on each element.
  ceed::Operator integ_op;

  // Temporary vectors for error estimation.
  mutable VecType B_gf, H, H_gf;

  Vector ComputeEstimates(const VecType &B) const;
  Vector RecoveredFluxNorms() const;

public:
  CurlFluxErrorEstimator(const MaterialOperator &mat_op, FiniteElementSpace &rt_fespace,
                         FiniteElementSpaceHierarchy &nd_fespaces, double tol, int max_it,
                         int print, bool use_mg);

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
                                  int max_it, int print, bool use_mg);

  void AddErrorIndicator(const VecType &E, const VecType &B, double Et,
                         ErrorIndicator &indicator) const;

  // The same at the real angular frequency ω of the fields, required with a
  // frequency-dependent material (see GradFluxErrorEstimator).
  void AddErrorIndicator(const VecType &E, const VecType &B, double Et, double omega,
                         ErrorIndicator &indicator) const
    requires std::is_same_v<VecType, ComplexVector>;
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
                                 int max_it, int print, bool use_mg);

  void AddErrorIndicator(const VecType &E, const VecType &B, double Et,
                         ErrorIndicator &indicator) const;

  // The same at the real angular frequency ω of the fields, required with a
  // frequency-dependent material (see GradFluxErrorEstimator).
  void AddErrorIndicator(const VecType &E, const VecType &B, double Et, double omega,
                         ErrorIndicator &indicator) const
    requires std::is_same_v<VecType, ComplexVector>;
};

}  // namespace palace

#endif  // PALACE_LINALG_ERROR_ESTIMATOR_HPP
