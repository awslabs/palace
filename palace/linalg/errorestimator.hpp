// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_ERROR_ESTIMATOR_HPP
#define PALACE_ERROR_ESTIMATOR_HPP

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

  mfem::RT_FECollection flux_fec;
  mutable mfem::ParFiniteElementSpace flux_fes;
  std::vector<std::unique_ptr<mfem::ND_FECollection>> smooth_flux_fecs;
  mutable mfem::ParFiniteElementSpaceHierarchy smooth_flux_fes;
  mutable FluxProjector projector;

public:
  // Constructor for using geometric and p multigrid
  CurlFluxErrorEstimator(const IoData &iodata, const MaterialOperator &mat_op,
                         const std::vector<std::unique_ptr<mfem::ParMesh>> &mesh,
                         mfem::ParFiniteElementSpace &fes);
  // Constructor for using only p multigrid
  CurlFluxErrorEstimator(const IoData &iodata, const MaterialOperator &mat_op,
                         std::unique_ptr<mfem::ParMesh> &mesh,
                         mfem::ParFiniteElementSpace &fes);

  // Compute elemental error indicators given a vector of true DOF, and the
  // finite element space they are associated with
  mfem::Vector operator()(const petsc::PetscParVector &v, bool use_mfem = true) const;
};

// Class used for computing grad flux error estimate,
// i.e. || ϵ ∇ ϕ - F ||_K
// where F denotes a smooth reconstruction of ϵ ∇ u
class GradFluxErrorEstimator
{
  const MaterialOperator &mat_op;
  // The finite element space used to represent ϕ
  mfem::ParFiniteElementSpace &fes;

  mfem::L2_FECollection flux_fec;
  mutable mfem::ParFiniteElementSpace flux_fes;
  std::vector<std::unique_ptr<mfem::H1_FECollection>> smooth_flux_fecs;
  mutable mfem::ParFiniteElementSpaceHierarchy smooth_flux_fes;
  mutable FluxProjector projector;

public:
  // Constructor for using geometric and p multigrid
  GradFluxErrorEstimator(const IoData &iodata, const MaterialOperator &mat_op,
                         const std::vector<std::unique_ptr<mfem::ParMesh>> &mesh,
                         mfem::ParFiniteElementSpace &fes);
  // Constructor for using only p multigrid
  GradFluxErrorEstimator(const IoData &iodata, const MaterialOperator &mat_op,
                         std::unique_ptr<mfem::ParMesh> &mesh,
                         mfem::ParFiniteElementSpace &fes);

  // Compute elemental error indicators given a vector of true DOF, and the
  // finite element space they are associated with
  mfem::Vector operator()(const mfem::Vector &v, bool use_mfem = true) const;
};

}  // namespace palace

#endif  // PALACE_ERROR_ESTIMATOR_HPP
