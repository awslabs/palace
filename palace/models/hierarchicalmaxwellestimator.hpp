// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_MODELS_HIERARCHICALMAXWELLESTIMATOR_HPP
#define PALACE_MODELS_HIERARCHICALMAXWELLESTIMATOR_HPP

#include <complex>
#include <memory>
#include <vector>
#include <mfem.hpp>

#include "fem/hierarchicalerrorestimator.hpp"
#include "fem/singularassembly.hpp"
#include "fem/singulardofs.hpp"

namespace palace
{

class FiniteElementSpace;
class SpaceOperator;

// Estimator-only p+1 domain data for a singular Maxwell SpaceOperator. This workspace
// assembles uneliminated element matrices only: it never constructs a global fine operator
// or performs a global fine solve. The first production slice is serial; the local algebra
// and stable singular numbering are shared with the future parallel adapter.
class HierarchicalMaxwellDomainData
{
private:
  std::unique_ptr<mfem::ND_FECollection> fine_nd_collection;
  std::unique_ptr<mfem::H1_FECollection> fine_h1_collection;
  std::unique_ptr<FiniteElementSpace> fine_nd_space;
  std::unique_ptr<FiniteElementSpace> fine_h1_space;

  fem::hierarchical::SparseInjection injection;
  std::vector<bool> coarse_essential, fine_essential;
  std::vector<std::vector<int>> element_enrichment_guests;

  struct ElementMatrices
  {
    int support_element = -1;
    std::vector<int> dofs;
    mfem::DenseMatrix curl_curl;
    mfem::DenseMatrix mass_real;
    mfem::DenseMatrix mass_imag;
    mfem::DenseMatrix mass_abs;
  };
  std::vector<ElementMatrices> elements;
  std::vector<fem::singular::LocalNDBoundaryPatchMatrices> boundary_stiffness;
  std::vector<fem::singular::LocalNDBoundaryPatchMatrices> boundary_damping;
  std::vector<fem::singular::LocalNDBoundaryPatchMatrices> boundary_mass;
  std::vector<fem::singular::LocalNDBoundaryPatchMatrices> boundary_stiffness_abs;
  std::vector<fem::singular::LocalNDBoundaryPatchMatrices> boundary_damping_abs;
  std::vector<fem::singular::LocalNDBoundaryPatchMatrices> boundary_mass_abs;
  bool unsupported_polynomial_boundary = false;

public:
  explicit HierarchicalMaxwellDomainData(SpaceOperator &space_op);
  ~HierarchicalMaxwellDomainData();

  const FiniteElementSpace &GetFineNDSpace() const { return *fine_nd_space; }
  FiniteElementSpace &GetFineNDSpace() { return *fine_nd_space; }
  const auto &GetInjection() const { return injection; }
  const auto &GetCoarseEssentialMask() const { return coarse_essential; }
  const auto &GetFineEssentialMask() const { return fine_essential; }
  const auto &GetElementEnrichmentGuests() const { return element_enrichment_guests; }
  int GetFineStandardSize() const;
  int GetEnrichmentSize() const;

  // Exact domain part of
  //   K - omega^2 (Mr + i Mi)
  // for genuinely complex omega. Boundary damping, lumped ports, and A2(omega) are appended
  // by the driver-level facet adapter before residual assembly.
  std::vector<fem::hierarchical::ComplexLocalOperatorContribution>
  BuildComplexDomainContributions(std::complex<double> omega) const;

  // Positive graph metric K + |omega|^2 M_abs used only to lift the exact complex residual.
  std::vector<fem::hierarchical::LocalOperatorContribution>
  BuildDomainMetricContributions(std::complex<double> omega) const;

  // Complete polynomial domain-plus-boundary operator
  //   K + i omega C - omega^2 (Mr + i Mi)
  // including singular lumped-port and supported surface-impedance facets. General
  // frequency-dependent A2(omega) terms are appended by the next driver adapter.
  std::vector<fem::hierarchical::ComplexLocalOperatorContribution>
  BuildComplexPolynomialContributions(std::complex<double> omega) const;

  // Positive domain-plus-boundary graph metric with absolute physical coefficients:
  //   K_abs + |omega| C_abs + |omega|^2 M_abs.
  std::vector<fem::hierarchical::LocalOperatorContribution>
  BuildPolynomialMetricContributions(std::complex<double> omega) const;
};

}  // namespace palace

#endif  // PALACE_MODELS_HIERARCHICALMAXWELLESTIMATOR_HPP
