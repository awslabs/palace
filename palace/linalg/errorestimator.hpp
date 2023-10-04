// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_LINALG_ERROR_ESTIMATOR_HPP
#define PALACE_LINALG_ERROR_ESTIMATOR_HPP

#include <memory>
#include <mfem.hpp>
#include "linalg/ksp.hpp"
#include "linalg/operator.hpp"
#include "linalg/vector.hpp"
#include "utils/errorindicators.hpp"

namespace palace
{

class IoData;
class MaterialOperator;

//
// This solver computes a smooth reconstruction of a discontinuous flux. The difference
// between this resulting smooth flux and the original non-smooth flux provides a
// localizable error estimate. An instance  of FluxProjector can be reused across solutions,
// thus the construction of the operator is separated from the construction of the flux RHS.
class FluxProjector
{
private:
  // Operator for the mass matrix inversion.
  std::unique_ptr<Operator> M;

  // Linear solver and preconditioner for the projected linear system M σ = σ̂.
  std::unique_ptr<KspSolver> ksp;

  // Intermediate storage vector used in Mult.
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

// Class used for computing curl flux error estimate,
// i.e. || μ⁻¹∇ × Vₕ - F ||_K
// where F denotes a smooth reconstruction of μ⁻¹∇ × Vₕ.
class CurlFluxErrorEstimator
{
  const MaterialOperator &mat_op;
  // The finite element space used to represent V.
  mfem::ParFiniteElementSpace &fespace;

  std::vector<std::unique_ptr<mfem::ND_FECollection>> smooth_flux_fecs;
  mutable mfem::ParFiniteElementSpaceHierarchy smooth_flux_fespace;
  mutable FluxProjector smooth_projector;

  mfem::L2_FECollection coarse_flux_fec;
  mutable mfem::ParFiniteElementSpace coarse_flux_fespace;

  std::vector<mfem::DenseMatrix> scalar_mass_matrices;
  std::vector<mfem::DenseMatrix> smooth_to_coarse_embed;

  mutable ComplexVector complex_flux;
  mutable Vector real_flux;

public:
  // Constructor for using geometric and p multigrid.
  CurlFluxErrorEstimator(const IoData &iodata, const MaterialOperator &mat_op,
                         const std::vector<std::unique_ptr<mfem::ParMesh>> &mesh,
                         mfem::ParFiniteElementSpace &fespace);

  // Compute elemental error indicators given a complex vector of true DOF.
  IndicatorsAndNormalization ComputeIndicators(const ComplexVector &v,
                                        bool normalize) const;

  // Compute elemental error indicators given a vector of true DOF.
  IndicatorsAndNormalization ComputeIndicators(const Vector &v, bool normalize ) const;
};

// Class used for computing grad flux error estimate,
// i.e. || ϵ ∇ ϕₕ - F ||_K
// where F denotes a smooth reconstruction of ϵ ∇ ϕₕ.
class GradFluxErrorEstimator
{
  const MaterialOperator &mat_op;
  // The finite element space used to represent ϕ.
  mfem::ParFiniteElementSpace &fespace;

  // Collections and spaces for the smooth flux. Note the hierarchy uses the
  // SCALAR finite element space, whilst the true flux is in the VECTOR finite
  // element space. This allows for a component wise inversion.
  std::vector<std::unique_ptr<mfem::H1_FECollection>> smooth_flux_fecs;
  mutable mfem::ParFiniteElementSpaceHierarchy smooth_flux_component_fespace;
  mutable mfem::ParFiniteElementSpace smooth_flux_fespace;
  mutable FluxProjector smooth_projector;

  mfem::L2_FECollection coarse_flux_fec;
  mutable mfem::ParFiniteElementSpace coarse_flux_fespace;

  std::vector<mfem::DenseMatrix> scalar_mass_matrices;
  std::vector<mfem::DenseMatrix> smooth_to_coarse_embed;

public:
  // Constructor for using geometric and p multigrid.
  GradFluxErrorEstimator(const IoData &iodata, const MaterialOperator &mat_op,
                         const std::vector<std::unique_ptr<mfem::ParMesh>> &mesh,
                         mfem::ParFiniteElementSpace &fespace);

  // Compute elemental error indicators given a vector of true DOF.
  IndicatorsAndNormalization ComputeIndicators(const Vector &v, bool normalize) const;
};

}  // namespace palace

#endif  // PALACE_LINALG_ERROR_ESTIMATOR_HPP
