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
#include "linalg/vector.hpp"

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
  SpaceOperator *space_op;
  std::unique_ptr<mfem::HypreParMatrix> enrichment_prolongation;
  int enrichment_size = 0;

  fem::hierarchical::SparseInjection injection;
  std::vector<bool> coarse_essential, fine_essential;
  std::vector<std::vector<int>> element_enrichment_guests;

  // Isotropic material scaling makes the imaginary and absolute mass blocks exact scalar
  // multiples of one unweighted mass block per element, so only two dense matrices are
  // stored: curl_curl (inverse-permeability weighted) and mass (unweighted), plus the
  // three electric scalars.
  struct ElementMatrices
  {
    int support_element = -1;
    std::vector<int> dofs;
    mfem::Array<int> standard_dofs;
    mfem::Array<int> enrichment_dofs;
    mfem::DenseMatrix curl_curl;
    mfem::DenseMatrix mass;
    double electric_real = 1.0;
    double electric_imag = 0.0;
    double electric_abs = 1.0;
  };
  std::vector<ElementMatrices> elements;
  std::vector<fem::singular::LocalNDBoundaryPatchMatrices> boundary_stiffness;
  std::vector<fem::singular::LocalNDBoundaryPatchMatrices> boundary_damping;
  std::vector<fem::singular::LocalNDBoundaryPatchMatrices> boundary_mass;
  std::vector<fem::singular::LocalNDBoundaryPatchMatrices> boundary_stiffness_abs;
  std::vector<fem::singular::LocalNDBoundaryPatchMatrices> boundary_damping_abs;
  std::vector<fem::singular::LocalNDBoundaryPatchMatrices> boundary_mass_abs;
  bool unsupported_polynomial_boundary = false;
  bool entity_patches_available = false;
  mfem::Array<int> fine_standard_essential_true_dofs;
  mfem::Array<int> enrichment_essential_true_dofs;

  std::vector<fem::singular::LocalNDElementPatchMatrices>
  BuildPolynomialMetricElementPatches(std::complex<double> omega) const;

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
  int GetEnrichmentSize() const { return enrichment_size; }

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

  struct ParallelEstimate
  {
    Vector indicator_energy;
    double total_energy = 0.0;
  };

  // Patch shape used to lift the exact complex residual. ELEMENT is the MPI-capable
  // production smoother shape; ENTITY is the certified edge/face/interior engine, which
  // measures far stronger marking but is currently available only on one rank with
  // identity local-to-true layouts.
  enum class PatchShape
  {
    ELEMENT,
    ENTITY
  };
  bool EntityPatchesAvailable() const { return entity_patches_available; }

  // One uneliminated element contribution of the positive metric, expressed on global
  // combined true DOFs with local orientation signs folded into the matrix. Records are
  // rank-portable: any rank can scatter them without knowing the source rank's local
  // numbering, which is the transport primitive for cross-rank entity-patch support
  // unions. Combined global numbering places all standard true DOFs first and appends the
  // global enrichment numbering.
  struct TrueElementRecord
  {
    HYPRE_BigInt element_id = -1;  // Owner-rank global element identity.
    std::vector<HYPRE_BigInt> fine_dofs;
    std::vector<HYPRE_BigInt> coarse_dofs;
    std::vector<char> fine_essential;
    std::vector<char> coarse_essential;
    mfem::DenseMatrix metric;
    // Element block of the coarse-to-fine injection in the same sign-folded true
    // convention (rows: coarse_dofs, columns: fine_dofs). Receiving ranks cannot evaluate
    // transfer matrices for ghost elements, so guest columns assemble from these blocks.
    mfem::DenseMatrix injection;
  };

  // Extract this rank's local metric contributions as global-true records. Facet
  // contributions are merged into their support element's record so each element appears
  // exactly once.
  std::vector<TrueElementRecord>
  BuildTrueMetricElementRecords(std::complex<double> omega) const;

  // Exchange a three-layer halo of records with every mesh-neighbor rank, so shared-entity
  // patch owners hold the complete support union (owned complement modes reach one element
  // ring, coarse guest columns one more, and their support elements a third). Received
  // records are deduplicated by global element identity; the returned set contains the
  // local records first, then ghosts. ghost_count reports how many ghosts were appended.
  std::vector<TrueElementRecord> ExchangeHaloRecords(std::vector<TrueElementRecord> records,
                                                     int *ghost_count = nullptr) const;

  // Rank-count-agnostic certified entity-patch lifting over augmented records: sorted
  // true-DOF complement bases per owned edge/face/interior entity, dense Cholesky patch
  // solves over complete halo support unions, globally averaged partition-of-unity
  // guests, and the same bounded undamped sweeps as the serial engine, applied
  // identically to both complex residual components. Returns per-owned-element metric
  // energies of the lifted correction.
  ParallelEstimate
  LiftTrueComplexResidualByEntityPatches(std::complex<double> omega,
                                         const Vector &residual_true_real,
                                         const Vector &residual_true_imag) const;

private:
  // Shared lifting tail: consumes the uneliminated local complex residual and applies one
  // fixed lifting map of the selected patch shape identically to both components.
  ParallelEstimate LiftLocalComplexResidual(std::complex<double> omega,
                                            fem::hierarchical::ComplexResidual residual,
                                            PatchShape patch_shape) const;

  // Streaming form of the exact complex polynomial residual b - A(omega) x on the
  // uneliminated local combined layout: element and facet operator blocks are formed on
  // the stack one at a time instead of materializing two dense matrices per element.
  fem::hierarchical::ComplexResidual
  AssembleStreamingComplexResidual(std::complex<double> omega, const Vector &injected_real,
                                   const Vector &injected_imag) const;

  // Reduce one uneliminated local combined component to owned true DOFs, zeroing
  // essential equations.
  void LocalToTrue(const Vector &local, Vector &combined_true) const;

public:
  // Port-excitation load b(omega) = i omega RHS1 evaluated against the p+1 test space, in
  // the uneliminated local combined layout. Standard rows assemble the same lumped-port
  // and surface-current boundary coefficients as the coarse excitation; enrichment rows
  // are order-independent and splice exactly from the coarse combined right-hand side.
  // Public so tests can check Galerkin consistency against the coarse excitation.
  fem::hierarchical::ComplexResidual
  AssembleFineExcitation(int excitation_idx, double omega,
                         const ComplexVector &coarse_rhs) const;

  // Parallel exact polynomial eigen-residual estimate (zero RHS). The coarse combined true
  // field is injected to the local p+1 standard/singular layout, all local element/facet
  // residuals are summed to parallel true DOFs, and one fixed lifting map (of the selected
  // patch shape) is applied identically to both complex components. No global p+1 matrix
  // or solve is constructed.
  ParallelEstimate
  EstimatePolynomialEigenResidual(std::complex<double> omega,
                                  const ComplexVector &coarse_field,
                                  PatchShape patch_shape = PatchShape::ELEMENT) const;

  // Driven counterpart: the residual is b(omega) - A(omega) x for the coarse solution x
  // and the exact excitation functional of the selected port drive, so a converged coarse
  // solve has exactly the hierarchical p-enrichment residual left over.
  ParallelEstimate
  EstimatePolynomialDrivenResidual(double omega, const ComplexVector &coarse_field,
                                   int excitation_idx, const ComplexVector &coarse_rhs,
                                   PatchShape patch_shape = PatchShape::ELEMENT) const;
};

}  // namespace palace

#endif  // PALACE_MODELS_HIERARCHICALMAXWELLESTIMATOR_HPP
