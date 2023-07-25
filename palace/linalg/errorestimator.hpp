// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_LINALG_ERROR_ESTIMATOR_HPP
#define PALACE_LINALG_ERROR_ESTIMATOR_HPP

#include <mfem.hpp>

#include "linalg/fluxprojector.hpp"

namespace palace
{

class IoData;
class MaterialOperator;

namespace petsc
{
class PetscParVector;
}  // namespace petsc

// Class used for computing curl flux error estimate,
// i.e. || μ⁻¹∇ × V - F ||_K
// where F denotes a smooth reconstruction of μ⁻¹∇ × V
class CurlFluxErrorEstimator
{
  const MaterialOperator &mat_op;
  // The finite element space used to represent V
  mfem::ParFiniteElementSpace &fes;

  std::vector<std::unique_ptr<mfem::ND_FECollection>> smooth_flux_fecs;
  mutable mfem::ParFiniteElementSpaceHierarchy smooth_flux_fes;
  mutable FluxProjector smooth_projector;

  mfem::L2_FECollection coarse_flux_fec;
  mutable mfem::ParFiniteElementSpace coarse_flux_fes;

  std::vector<mfem::DenseMatrix> scalar_mass_matrices;
  std::vector<mfem::DenseMatrix> smooth_to_coarse_embed;

public:
  // Constructor for using geometric and p multigrid.
  CurlFluxErrorEstimator(const IoData &iodata, const MaterialOperator &mat_op,
                         const std::vector<std::unique_ptr<mfem::ParMesh>> &mesh,
                         mfem::ParFiniteElementSpace &fes);

  // Compute elemental error indicators given a complex vector of true DOF.
  Vector operator()(const ComplexVector &v) const;

  // Compute elemental error indicators given a vector of true DOF.
  Vector operator()(const Vector &v) const;
};

// Class used for computing grad flux error estimate,
// i.e. || ϵ ∇ ϕ - F ||_K
// where F denotes a smooth reconstruction of ϵ ∇ ϕ
class GradFluxErrorEstimator
{
  const MaterialOperator &mat_op;
  // The finite element space used to represent ϕ
  mfem::ParFiniteElementSpace &fes;

  // Collections and spaces for the smooth flux. Note the hierarchy uses the
  // SCALAR finite element space, whilst the true flux is in the VECTOR finite
  // element space. This allows for a component wise inversion.
  std::vector<std::unique_ptr<mfem::H1_FECollection>> smooth_flux_fecs;
  mutable mfem::ParFiniteElementSpaceHierarchy smooth_flux_component_fes;
  mutable mfem::ParFiniteElementSpace smooth_flux_fes;
  mutable FluxProjector smooth_projector;

  mfem::L2_FECollection coarse_flux_fec;
  mutable mfem::ParFiniteElementSpace coarse_flux_fes;

  std::vector<mfem::DenseMatrix> scalar_mass_matrices;
  std::vector<mfem::DenseMatrix> smooth_to_coarse_embed;

public:
  // Constructor for using geometric and p multigrid.
  GradFluxErrorEstimator(const IoData &iodata, const MaterialOperator &mat_op,
                         const std::vector<std::unique_ptr<mfem::ParMesh>> &mesh,
                         mfem::ParFiniteElementSpace &fes);

  // Compute elemental error indicators given a vector of true DOF.
  Vector operator()(const Vector &v) const;
};

}  // namespace palace

#endif  // PALACE_LINALG_ERROR_ESTIMATOR_HPP
