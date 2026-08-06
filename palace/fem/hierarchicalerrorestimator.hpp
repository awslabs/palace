// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_FEM_HIERARCHICALERRORESTIMATOR_HPP
#define PALACE_FEM_HIERARCHICALERRORESTIMATOR_HPP

#include <set>
#include <vector>
#include <mfem.hpp>

namespace palace::fem::hierarchical
{
// Sparse coefficients of one complete conforming patch direction in a combined standard
// plus singular true-DOF numbering. The owning estimator is responsible for orientation and
// constraint transformations before constructing this column.
struct SparseColumn
{
  std::vector<int> dofs;
  std::vector<double> values;
};

// One uneliminated element or boundary-facet contribution. Several contributions may share
// the same support element (for example domain curl/mass plus impedance facets).
struct LocalOperatorContribution
{
  int support_element = -1;
  std::vector<int> dofs;
  mfem::DenseMatrix matrix;
  mfem::Vector rhs;
};

// Complex counterpart used for the complete frequency-domain residual. Splitting local
// matrices and loads into real and imaginary slots preserves all cross-coupling from
//
//   A(omega) = K + i omega C - omega^2 (Mr + i Mi) + A2(omega).
//
// The residual assembler below is agnostic to how those slots were formed, so genuinely
// complex frequencies and non-polynomial A2 terms use the exact same path.
struct ComplexLocalOperatorContribution
{
  int support_element = -1;
  std::vector<int> dofs;
  mfem::DenseMatrix matrix_real;
  mfem::DenseMatrix matrix_imag;
  mfem::Vector rhs_real;
  mfem::Vector rhs_imag;
};

struct ComplexResidual
{
  mfem::Vector real;
  mfem::Vector imag;
};

// Assemble b-Ax without constructing A. Essential equations are set to zero after all local
// contributions are scattered.
mfem::Vector AssembleResidual(int combined_size,
                              const std::vector<LocalOperatorContribution> &contributions,
                              const mfem::Vector &injected,
                              const std::vector<bool> &essential);

// Assemble the exact coupled complex residual. For A = Ar + i Ai and x = xr + i xi,
//
//   rr = br - Ar xr + Ai xi,
//   ri = bi - Ai xr - Ar xi.
//
// Essential equations are zeroed only after every element and facet contribution has been
// scattered.
ComplexResidual AssembleComplexResidual(
    int combined_size, const std::vector<ComplexLocalOperatorContribution> &contributions,
    const mfem::Vector &injected_real, const mfem::Vector &injected_imag,
    const std::vector<bool> &essential);

// Assemble B^T A B and B^T r over a complete support union. Returns the number of local
// matrix entries touched, a deterministic work proxy.
long long
AssembleRestrictedOperator(const std::vector<LocalOperatorContribution> &contributions,
                           const std::set<int> &support_elements,
                           const std::vector<SparseColumn> &basis,
                           const mfem::Vector &residual, mfem::DenseMatrix &restricted,
                           mfem::Vector &restricted_rhs);

// Evaluate x^T A x by local contributions. This is also used for support-local scalar line
// searches and element indicators.
double Energy(const std::vector<LocalOperatorContribution> &contributions,
              const mfem::Vector &vector);

// Sparse signed p -> p+1 injection between two spaces on the same mesh, assembled
// element-by-element with MFEM DOF transformations applied on both sides. Columns are
// indexed by unsigned coarse VDof; entries address unsigned fine VDofs. The consistency
// error is the largest disagreement between elements that share a coarse basis function
// and is an orientation-covariance diagnostic: it vanishes when signed VDofs and DOF
// transformations are handled correctly.
struct SparseInjection
{
  std::vector<SparseColumn> columns;
  double consistency_error = 0.0;
  int signed_coarse_dofs = 0;
  int signed_fine_dofs = 0;
  int nonidentity_transformations = 0;
};

SparseInjection BuildSparsePInjection(mfem::Mesh &mesh,
                                      mfem::FiniteElementSpace &coarse_space,
                                      mfem::FiniteElementSpace &fine_space);

struct PatchLiftingOptions
{
  // Total additive-Schwarz defect-correction passes, including the first patch solve.
  int sweeps = 4;
  // Sparse entries below this magnitude are dropped from complement basis columns.
  double drop_tolerance = 1.0e-13;
  // Residual-dependent scalar line searches improve a single real lift but make separate
  // real/imaginary lifts depend on the arbitrary global complex phase. Complex lifting
  // therefore forces this false and applies the same fixed linear patch inverse to both
  // components.
  bool residual_line_search = true;
};

struct PatchLiftingReport
{
  std::vector<double> indicator;
  double energy = 0.0;
  double work = 0.0;
  double maximum_patch_residual = 0.0;
  double maximum_patch_condition = 0.0;
  int edge_patches = 0;
  int face_patches = 0;
  int interior_patches = 0;
  int owned_modes = 0;
  int maximum_support_elements = 0;
  int maximum_patch_dimension = 0;
  int maximum_element_overlap = 0;
  // Number of patches hosting each coarse standard DOF (by unsigned coarse VDof) or
  // enrichment DOF as an averaged partition-of-unity guest. A zero count for a free DOF
  // means that direction is not covered by any patch.
  std::vector<int> coarse_guest_counts;
  std::vector<int> enrichment_guest_counts;
};

// Element-local hierarchical residual lifting on the combined
// [fine standard, enrichment] layout. Every conforming p -> p+1 complement direction is
// owned by exactly one edge, face (3D), or element-interior patch; incident coarse
// standard and enrichment DOFs join their owner patches as partition-of-unity guests
// averaged by stable DOF identity. Patch corrections solve the restriction of the
// coercive metric contributions by dense Cholesky over measured complete support unions;
// the (possibly indefinite) residual contributions only define the lifted right-hand
// side. No global fine matrix is assembled and no global fine solve is performed.
//
// The coarse combined vector uses the [coarse standard, enrichment] layout, and both
// essential masks cover their full combined ranges. element_enrichment_guests lists, for
// every mesh element, the enrichment DOFs a patch owning that element must host; family
// filtering (for example excluding rotational guests) is the caller's responsibility.
PatchLiftingReport EstimateByPatchLifting(
    mfem::Mesh &mesh, mfem::FiniteElementSpace &coarse_space,
    mfem::FiniteElementSpace &fine_space, const SparseInjection &injection,
    const std::vector<LocalOperatorContribution> &residual_contributions,
    const std::vector<LocalOperatorContribution> &metric_contributions,
    const std::vector<bool> &fine_essential, const std::vector<bool> &coarse_essential,
    const mfem::Vector &coarse_combined,
    const std::vector<std::vector<int>> &element_enrichment_guests,
    const PatchLiftingOptions &options = {});

// Lower-level lifting for one real residual. Complex callers should use
// LiftComplexResidualByPatches rather than calling this independently for each component:
// the latter disables residual-dependent scalar line searches so the summed indicator is
// invariant under an arbitrary global complex phase.
PatchLiftingReport LiftResidualByPatches(
    mfem::Mesh &mesh, mfem::FiniteElementSpace &coarse_space,
    mfem::FiniteElementSpace &fine_space, const SparseInjection &injection,
    const std::vector<LocalOperatorContribution> &metric_contributions,
    const std::vector<bool> &fine_essential, const std::vector<bool> &coarse_essential,
    const mfem::Vector &residual,
    const std::vector<std::vector<int>> &element_enrichment_guests,
    const PatchLiftingOptions &options = {});

struct ComplexPatchLiftingReport
{
  PatchLiftingReport real;
  PatchLiftingReport imag;
  std::vector<double> indicator;
  double energy = 0.0;
};

// Apply one fixed linear approximate B^-1 to both components of the exact complex
// residual, then sum their B-energy indicators. Because the same real linear map acts on
// both components, the result is invariant under r -> exp(i phi) r.
ComplexPatchLiftingReport LiftComplexResidualByPatches(
    mfem::Mesh &mesh, mfem::FiniteElementSpace &coarse_space,
    mfem::FiniteElementSpace &fine_space, const SparseInjection &injection,
    const std::vector<LocalOperatorContribution> &metric_contributions,
    const std::vector<bool> &fine_essential, const std::vector<bool> &coarse_essential,
    const ComplexResidual &residual,
    const std::vector<std::vector<int>> &element_enrichment_guests,
    const PatchLiftingOptions &options = {});

}  // namespace palace::fem::hierarchical

#endif  // PALACE_FEM_HIERARCHICALERRORESTIMATOR_HPP
