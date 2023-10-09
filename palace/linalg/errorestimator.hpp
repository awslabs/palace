// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_LINALG_ERROR_ESTIMATOR_HPP
#define PALACE_LINALG_ERROR_ESTIMATOR_HPP

#include <memory>
#include <mfem.hpp>
#include "fem/errorindicator.hpp"
#include "linalg/ksp.hpp"
#include "linalg/operator.hpp"
#include "linalg/vector.hpp"

namespace palace
{

class IoData;
class MaterialOperator;
class PostOperator;

//
// This solver computes a smooth reconstruction of a discontinuous flux. The difference
// between this resulting smooth flux and the original non-smooth flux provides a
// localizable error estimate. An instance  of FluxProjector can be reused across solutions,
// thus the construction of the operator is separated from the construction of the flux RHS.

template <typename SmoothFluxFiniteElementCollection>
class FluxProjector
{
private:
  // Operator for the mass matrix inversion.
  std::unique_ptr<Operator> M;

  // Linear solver and preconditioner for the projected linear system M σ = σ̂.
  std::unique_ptr<KspSolver> ksp;

public:
  FluxProjector(mfem::ParFiniteElementSpaceHierarchy &smooth_flux_fes, double tol,
                int max_it, int print_level, int pa_order_threshold);

  inline void Mult(const Vector &x, Vector &y) const { ksp->Mult(x, y); }
  inline void Mult(const ComplexVector &x, ComplexVector &y) const
  {
    Mult(x.Real(), y.Real());
    Mult(x.Imag(), y.Imag());
  }
};

// Class used for computing curl flux error estimate, i.e. || μ⁻¹∇ × Vₕ - F ||_K where F
// denotes a smooth reconstruction of μ⁻¹∇ × Vₕ.
class CurlFluxErrorEstimator
{
  const MaterialOperator &mat_op;
  // The finite element space used to represent V, and F.
  mfem::ParFiniteElementSpaceHierarchy &nd_fespaces;
  FluxProjector<mfem::ND_FECollection> smooth_projector;
  mutable Vector smooth_flux, flux_rhs;
  mutable mfem::ParGridFunction field_gf, smooth_flux_gf;

public:
  // Constructor for using geometric and p multigrid.
  CurlFluxErrorEstimator(const IoData &iodata, const MaterialOperator &mat_op,
                         mfem::ParFiniteElementSpaceHierarchy &nd_fespaces);

  // Compute elemental error indicators given a vector of true DOF.
  template <typename VectorType>
  ErrorIndicator ComputeIndicators(const VectorType &v) const;

  // Compute elemental error indicators given a vector of true DOF, v, and fold into an
  // existing indicator. Optionally set the error indicator field within a PostOperator.
  template <typename VectorType>
  void AddErrorIndicator(ErrorIndicator &indicator, PostOperator &postop,
                         const VectorType &v) const;
  template <typename VectorType>
  void AddErrorIndicator(ErrorIndicator &indicator, const VectorType &v) const;
};

// Class used for computing grad flux error estimate, i.e. || ϵ ∇ ϕₕ - F ||_K where F
// denotes a smooth reconstruction of ϵ ∇ ϕₕ.
class GradFluxErrorEstimator
{
  const MaterialOperator &mat_op;
  // The finite element space used to represent ϕ, and components of F
  mfem::ParFiniteElementSpaceHierarchy &h1_fespaces;

  FluxProjector<mfem::H1_FECollection> smooth_projector;
  mutable Vector smooth_flux, flux_rhs;
  mutable mfem::ParGridFunction field_gf, smooth_flux_gf;

public:
  // Constructor for using geometric and p multigrid.
  GradFluxErrorEstimator(const IoData &iodata, const MaterialOperator &mat_op,
                         mfem::ParFiniteElementSpaceHierarchy &h1_fespaces);

  // Compute elemental error indicators given a vector of true DOF.
  ErrorIndicator ComputeIndicators(const Vector &v) const;

  // Compute elemental error indicators given a vector of true DOF, v, and fold into an
  // existing indicator.
  void AddErrorIndicator(ErrorIndicator &indicator, const Vector &v) const
  {
    indicator.AddIndicator(ComputeIndicators(v));
  }
};

}  // namespace palace

#endif  // PALACE_LINALG_ERROR_ESTIMATOR_HPP
