// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_FEM_SINGULARASSEMBLY_HPP
#define PALACE_FEM_SINGULARASSEMBLY_HPP

#include <cstddef>
#include <map>
#include <memory>
#include <set>
#include <string>
#include <vector>
#include <mfem.hpp>

#include "fem/singulardofs.hpp"
#include "fem/singularelements.hpp"
#include "fem/singulargeometry.hpp"

namespace palace
{

namespace fem
{

namespace singular
{

inline constexpr int H1DuffyReferenceOrder = 47;
inline constexpr int H1DuffyComparisonOrder = 39;
inline constexpr double H1DuffyRadialPower = 6.0;
inline constexpr double H1DuffyErrorSafetyFactor = 8.0;
// This is divisible by the denominators of both supported canonical exponents,
// nu = 1/2 and 2/3, so the transformed radial products are polynomial.
inline constexpr double TriangleDuffyRadialPower = 6.0;

struct AdaptiveAssemblyOptions
{
  int quadrature_order;
  double absolute_tolerance;
  double relative_tolerance;
  int maximum_subdivisions;
  bool fixed_subdivision = false;
  int subdivisions = 6;
  std::string reference_cache;
};

// Real, piecewise-constant isotropic operator coefficients for one element.
// The electric coefficient is shared by H1 diffusion and ND mass so material
// weighting preserves the enriched H1-to-ND exact sequence. It may be signed or
// zero to assemble one Cartesian component of a complex permittivity. The
// inverse-magnetic coefficient remains positive.
struct IsotropicMaterialCoefficients
{
  double electric;
  double inverse_magnetic;
};

struct ElementEnrichmentMatrices
{
  // For scalar singular potentials, diffusion is exactly the submatrix of the
  // ND mass matrix corresponding to their analytic gradient bases.
  mfem::DenseMatrix h1_diffusion;
  mfem::DenseMatrix h1_diffusion_estimated_absolute_error;
  mfem::DenseMatrix h1_mass;
  mfem::DenseMatrix h1_mass_estimated_absolute_error;

  mfem::DenseMatrix nd_mass;
  mfem::DenseMatrix nd_mass_estimated_absolute_error;
  mfem::DenseMatrix nd_curl_curl;
  mfem::DenseMatrix nd_curl_curl_estimated_absolute_error;

  std::size_t total_quadrature_leaf_count = 0;
  int maximum_subdivision_depth = 0;
};

struct ElementStandardEnrichmentMatrices
{
  // Rows of the standard-enrichment matrices use the MFEM element-local basis
  // ordering. The reverse blocks are stored as exact transposes.
  mfem::DenseMatrix h1_standard_enrichment;
  mfem::DenseMatrix h1_enrichment_standard;
  mfem::DenseMatrix h1_estimated_absolute_error;
  mfem::DenseMatrix h1_mass_standard_enrichment;
  mfem::DenseMatrix h1_mass_enrichment_standard;
  mfem::DenseMatrix h1_mass_estimated_absolute_error;

  mfem::DenseMatrix nd_mass_standard_enrichment;
  mfem::DenseMatrix nd_mass_enrichment_standard;
  mfem::DenseMatrix nd_mass_estimated_absolute_error;

  mfem::DenseMatrix nd_curl_curl_standard_enrichment;
  mfem::DenseMatrix nd_curl_curl_enrichment_standard;
  mfem::DenseMatrix nd_curl_curl_estimated_absolute_error;

  std::size_t total_quadrature_leaf_count = 0;
  int maximum_subdivision_depth = 0;
  std::size_t affine_nd_mass_contraction_count = 0;
  std::size_t affine_nd_mass_reintegration_count = 0;
  std::size_t affine_nd_mass_reintegration_batch_count = 0;
  std::size_t affine_nd_curl_contraction_count = 0;
  std::size_t affine_nd_curl_reintegration_count = 0;
  std::size_t affine_nd_curl_reintegration_batch_count = 0;
  double setup_time = 0.0;
  double nd_coupling_time = 0.0;
  double h1_gradient_coupling_time = 0.0;
  double h1_mass_coupling_time = 0.0;
};

struct ElementH1EnrichmentMatrices
{
  mfem::DenseMatrix enrichment_enrichment;
  mfem::DenseMatrix enrichment_enrichment_estimated_absolute_error;
  mfem::DenseMatrix standard_enrichment;
  mfem::DenseMatrix enrichment_standard;
  mfem::DenseMatrix standard_enrichment_estimated_absolute_error;

  std::size_t total_quadrature_leaf_count = 0;
  int maximum_subdivision_depth = 0;
};

struct LocalSparseOperatorBlocks
{
  std::unique_ptr<mfem::SparseMatrix> enrichment_enrichment;
  std::unique_ptr<mfem::SparseMatrix> standard_enrichment;
  std::unique_ptr<mfem::SparseMatrix> enrichment_standard;

  std::unique_ptr<mfem::SparseMatrix> enrichment_enrichment_estimated_absolute_error;
  std::unique_ptr<mfem::SparseMatrix> standard_enrichment_estimated_absolute_error;

  // Exact diagonal after subtracting the conforming standard-space interpolant
  // from every enrichment basis. Entries contain assembled element
  // contributions in the rank-local canonical enrichment ordering.
  std::unique_ptr<mfem::Vector> transformed_enrichment_diagonal;
};

// Material-weighted element matrices retained only for coupled Maxwell
// smoothing. Standard indices use MFEM's signed local-DOF convention;
// enrichment indices address the rank-local canonical enrichment vector.
struct LocalNDElementPatchMatrices
{
  int element = -1;
  mfem::Array<int> standard_dofs;
  mfem::Array<int> enrichment_dofs;
  mfem::DenseMatrix mass;
  mfem::DenseMatrix mass_estimated_absolute_error;
  mfem::DenseMatrix curl_curl;
  mfem::DenseMatrix curl_curl_estimated_absolute_error;
};

// One material-weighted boundary mass contribution in final unsigned local combined-DOF
// coordinates. Standard trace orientation and tetrahedral face DofTransformation have
// already been applied. Boundary terms sharing an adjacent support element remain separate
// so the hierarchical residual can include every lumped/impedance facet exactly once.
struct LocalNDBoundaryPatchMatrices
{
  int boundary = -1;
  int element = -1;
  std::vector<int> dofs;
  mfem::DenseMatrix mass;
};

struct LocalSparseEnrichmentMatrices
{
  LocalSparseOperatorBlocks h1_diffusion;
  LocalSparseOperatorBlocks h1_mass;
  LocalSparseOperatorBlocks nd_mass;
  LocalSparseOperatorBlocks nd_curl_curl;
  std::vector<LocalNDElementPatchMatrices> nd_element_patches;

  std::size_t total_quadrature_leaf_count = 0;
  int maximum_subdivision_depth = 0;
  std::size_t affine_reference_table_entries = 0;
  std::size_t affine_reference_cache_hits = 0;
  std::size_t affine_nd_mass_contraction_count = 0;
  std::size_t affine_nd_mass_reintegration_count = 0;
  std::size_t affine_nd_mass_reintegration_batch_count = 0;
  std::size_t affine_nd_curl_contraction_count = 0;
  std::size_t affine_nd_curl_reintegration_count = 0;
  std::size_t affine_nd_curl_reintegration_batch_count = 0;
  double enrichment_evaluation_time = 0.0;
  double standard_enrichment_evaluation_time = 0.0;
  double standard_reference_generation_time = 0.0;
  double standard_enrichment_setup_time = 0.0;
  double nd_coupling_time = 0.0;
  double h1_gradient_coupling_time = 0.0;
  double h1_mass_coupling_time = 0.0;
  double material_transformation_time = 0.0;
  double sparse_insertion_time = 0.0;
  double sparse_finalization_time = 0.0;
};

struct LocalSparseH1EnrichmentMatrices
{
  LocalSparseOperatorBlocks diffusion;

  std::size_t total_quadrature_leaf_count = 0;
  int maximum_subdivision_depth = 0;
  std::size_t duffy_reference_table_entries = 0;
  std::size_t duffy_reference_cache_hits = 0;
  std::size_t affine_reference_table_entries = 0;
  std::size_t affine_reference_cache_hits = 0;
  std::size_t affine_reference_pattern_count = 0;
  std::size_t affine_reference_pattern_hits = 0;
  std::size_t affine_reference_persistent_hits = 0;
  std::size_t affine_reference_persistent_writes = 0;
  std::size_t affine_reference_generated_leaf_count = 0;
  double affine_reference_generation_time = 0.0;
};

struct ParallelSparseOperatorBlocks
{
  std::unique_ptr<mfem::HypreParMatrix> enrichment_enrichment;
  std::unique_ptr<mfem::HypreParMatrix> standard_enrichment;
  std::unique_ptr<mfem::HypreParMatrix> enrichment_standard;

  std::unique_ptr<mfem::HypreParMatrix> enrichment_enrichment_estimated_absolute_error;
  std::unique_ptr<mfem::HypreParMatrix> standard_enrichment_estimated_absolute_error;

  // True-DOF assembly of LocalSparseOperatorBlocks::transformed_enrichment_diagonal.
  std::unique_ptr<mfem::Vector> transformed_enrichment_diagonal;
};

struct LocalInterpolatedNDDiagonals
{
  mfem::Vector mass;
  mfem::Vector curl_curl;
};

struct ParallelSparseEnrichmentMatrices
{
  ParallelSparseOperatorBlocks h1_diffusion;
  ParallelSparseOperatorBlocks h1_mass;
  ParallelSparseOperatorBlocks nd_mass;
  ParallelSparseOperatorBlocks nd_curl_curl;
};

// Return the constant physical barycentric gradients and positive Jacobian
// determinant for an actually affine three-dimensional tetrahedron. A
// high-order nodal representation of an affine map is accepted.
BarycentricGradients
GetAffineBarycentricGradients(mfem::ElementTransformation &transformation,
                              double &jacobian_determinant);

// Return constant physical barycentric gradients and the positive area
// Jacobian for an actually affine two-dimensional triangle. A high-order nodal
// representation of an affine map is accepted.
TriangleBarycentricGradients
GetAffineTriangleBarycentricGradients(mfem::ElementTransformation &transformation,
                                      double &jacobian_determinant);

// Construct the parallel map from uniquely owned enrichment true DOFs to the
// rank-local canonical enrichment DOFs used by element maps.
std::unique_ptr<mfem::HypreParMatrix>
BuildParallelEnrichmentProlongation(MPI_Comm comm, const TrueDofMap &dofs);

std::unique_ptr<mfem::Vector> AssembleParallelEnrichmentVector(MPI_Comm comm,
                                                               const TrueDofMap &dofs,
                                                               const mfem::Vector &local);

// Pull a true-DOF standard-space interpolant back to the two rank-local
// vectors used during element assembly. The result maps canonical local
// enrichment coefficients to standard local finite-element coefficients.
std::unique_ptr<mfem::SparseMatrix>
BuildLocalEnrichmentInterpolant(const mfem::HypreParMatrix &true_interpolant,
                                const mfem::ParFiniteElementSpace &standard_fespace,
                                const TrueDofMap &enrichment_numbering);

// Assemble d_e^T A_ss d_e elementwise for every column d_e of a local
// enrichment interpolant. One mass and curl-curl vector is returned for each
// material batch.
std::vector<LocalInterpolatedNDDiagonals> AssembleLocalInterpolatedNDDomainDiagonals(
    const DofTopology &topology, mfem::FiniteElementSpace &nd_fespace,
    const mfem::SparseMatrix &local_interpolant,
    const std::vector<std::vector<IsotropicMaterialCoefficients>> &material_batches);
std::vector<LocalInterpolatedNDDiagonals> AssembleLocalInterpolatedNDDomainDiagonals(
    const TriangleDofTopology &topology, mfem::FiniteElementSpace &nd_fespace,
    const mfem::SparseMatrix &local_interpolant,
    const std::vector<std::vector<IsotropicMaterialCoefficients>> &material_batches);

mfem::Vector AssembleLocalInterpolatedNDBoundaryDiagonal(
    const DofTopology &topology, mfem::FiniteElementSpace &nd_fespace,
    const mfem::SparseMatrix &local_interpolant,
    const std::map<int, double> &boundary_coefficients);
mfem::Vector AssembleLocalInterpolatedNDBoundaryDiagonal(
    const TriangleDofTopology &topology, mfem::FiniteElementSpace &nd_fespace,
    const mfem::SparseMatrix &local_interpolant,
    const std::map<int, double> &boundary_coefficients);

// Complete a local transformed diagonal,
//
//   diag(A_ee - D^T A_se - A_es D + D^T A_ss D),
//
// from an already assembled singular block and the elementwise standard-space
// energy above.
void SetLocalTransformedEnrichmentDiagonal(
    LocalSparseOperatorBlocks &blocks, const mfem::SparseMatrix &local_interpolant,
    const mfem::Vector &interpolated_standard_diagonal);

// Assemble enrichment-enrichment element matrices on one tetrahedron. The
// explicit-gradient overload is the affine reference-tensor path. The
// transformation overload detects high-order nodal maps which are actually
// affine and otherwise evaluates the physical Jacobian pointwise under
// singularity-aligned reference quadrature.
//
// Scalar material coefficients are deliberately excluded and must be applied
// by the caller. Every upper-triangular basis pair is integrated once and
// mirrored exactly.
//
// The routine throws if any adaptive reference-tensor entry fails its requested
// tolerance. Returned error matrices conservatively contract the entrywise
// reference-tensor estimates with the physical geometry factors.
ElementEnrichmentMatrices AssembleElementEnrichmentMatrices(
    const ElementDofMap &element_dofs, const BarycentricGradients &grad_lambda,
    double jacobian_determinant, const AdaptiveAssemblyOptions &options);

ElementEnrichmentMatrices
AssembleElementEnrichmentMatrices(const ElementDofMap &element_dofs,
                                  mfem::ElementTransformation &transformation,
                                  const AdaptiveAssemblyOptions &options);

ElementEnrichmentMatrices
AssembleTriangleElementEnrichmentMatrices(const TriangleElementDofMap &element_dofs,
                                          const TriangleBarycentricGradients &grad_lambda,
                                          double jacobian_determinant,
                                          const AdaptiveAssemblyOptions &options);

ElementEnrichmentMatrices
AssembleTriangleElementEnrichmentMatrices(const TriangleElementDofMap &element_dofs,
                                          mfem::ElementTransformation &transformation,
                                          const AdaptiveAssemblyOptions &options);

// Assemble the standard-enrichment coupling blocks using the actual MFEM H1
// and ND simplex basis conventions. The finite elements must have matching
// positive order. Scalar material coefficients and element DOF orientation
// transformations are deliberately excluded and must be applied by the
// caller.
//
// ND entries are integrated adaptively in physical coordinates. H1 coupling
// is formed as G^T M using MFEM's element-local discrete gradient G and the
// matching singular gradient columns of the ND mass coupling. This preserves
// the enriched H1-to-ND exact sequence algebraically.
ElementStandardEnrichmentMatrices AssembleElementStandardEnrichmentMatrices(
    const ElementDofMap &element_dofs, const mfem::FiniteElement &h1_fe,
    const mfem::FiniteElement &nd_fe, mfem::ElementTransformation &transformation,
    const AdaptiveAssemblyOptions &options);

ElementStandardEnrichmentMatrices AssembleTriangleElementStandardEnrichmentMatrices(
    const TriangleElementDofMap &element_dofs, const mfem::FiniteElement &h1_fe,
    const mfem::FiniteElement &nd_fe, mfem::ElementTransformation &transformation,
    const AdaptiveAssemblyOptions &options);

// Assemble only the H1 diffusion blocks needed by electrostatics. For a
// first-order standard H1 element, standard-enrichment and same-feature
// enrichment-enrichment entries use geometry-independent reference tensors
// from independently compared feature-aligned Duffy rules. Distinct-feature
// enrichment pairs use a two-chart partition-of-unity Duffy rule. Higher-order
// standard H1 coupling uses direct adaptive physical gradient inner products.
// The paired ND enrichment still defines the exact combined gradient map, but
// unused rotational and curl-curl blocks are deliberately excluded from this
// path.
ElementH1EnrichmentMatrices AssembleElementH1EnrichmentMatrices(
    const ElementDofMap &element_dofs, const mfem::FiniteElement &h1_fe,
    mfem::ElementTransformation &transformation, const AdaptiveAssemblyOptions &options);

// Apply the MFEM element-to-global dual transformations to the standard rows
// of all coupling blocks. Quadrature bounds are propagated with the absolute
// transformation matrices, and reverse blocks are reconstructed as exact
// transposes. Signed VDof indices are not part of DofTransformation and must
// still be honored when the returned matrices are scattered.
void ApplyStandardDofTransformations(const mfem::DofTransformation &h1_dof_transformation,
                                     const mfem::DofTransformation &nd_dof_transformation,
                                     ElementStandardEnrichmentMatrices &matrices);

// Apply one element's isotropic material coefficients to the geometry-only
// local blocks. Estimated absolute errors are multiplied by the corresponding
// positive coefficient. Invalid coefficients, malformed matrices, nonfinite
// entries, and overflow are rejected before any matrix is modified.
void ApplyIsotropicMaterialCoefficients(const IsotropicMaterialCoefficients &coefficients,
                                        ElementEnrichmentMatrices &matrices);
void ApplyIsotropicMaterialCoefficients(const IsotropicMaterialCoefficients &coefficients,
                                        ElementStandardEnrichmentMatrices &matrices);

// Assemble local sparse L-vector blocks over a conforming simplex mesh. The
// standard spaces and topology must refer to the same mesh, and one positive
// isotropic material pair is supplied per element. Standard signed VDofs and
// MFEM DOF transformations are applied to values; error bounds use the
// absolute transformations and unsigned VDof indices. Reverse coupling blocks
// are generated as exact sparse transposes.
//
// MFEM finite-element evaluation mutates transformation scratch state, so this
// correctness path is serial until an explicitly thread-safe evaluator exists.
LocalSparseEnrichmentMatrices AssembleLocalSparseEnrichmentMatrices(
    const DofTopology &topology, mfem::FiniteElementSpace &h1_fespace,
    mfem::FiniteElementSpace &nd_fespace,
    const std::vector<IsotropicMaterialCoefficients> &materials,
    const AdaptiveAssemblyOptions &options);

LocalSparseEnrichmentMatrices AssembleLocalSparseEnrichmentMatrices(
    const TriangleDofTopology &topology, mfem::FiniteElementSpace &h1_fespace,
    mfem::FiniteElementSpace &nd_fespace,
    const std::vector<IsotropicMaterialCoefficients> &materials,
    const AdaptiveAssemblyOptions &options);

// Assemble several material-weighted operators while evaluating each element's
// nonpolynomial basis tensors and adaptive quadrature only once. The outer result
// ordering matches material_batches. If retained_patch_batch is nonnegative,
// combined standard-plus-enrichment element matrices are retained only in that
// result entry for use by an overlapping Maxwell smoother. Pass
// RetainAllNDElementPatchBatches to retain them in every result. Pass
// RetainNDElementPatchStripsOnly for the estimator-only extraction mode: retain only the
// standard-enrichment and enrichment-enrichment columns in every result, omit
// reconstructible standard blocks and quadrature-error matrices, and leave the returned
// sparse operators and aggregate quadrature counts empty. Timing/cache diagnostics remain
// populated. This mode avoids duplicating data already owned by the production
// SpaceOperator.
inline constexpr int RetainNDElementPatchStripsOnly = -3;
inline constexpr int RetainAllNDElementPatchBatches = -2;
std::vector<LocalSparseEnrichmentMatrices> AssembleLocalSparseEnrichmentMatricesBatch(
    const DofTopology &topology, mfem::FiniteElementSpace &h1_fespace,
    mfem::FiniteElementSpace &nd_fespace,
    const std::vector<std::vector<IsotropicMaterialCoefficients>> &material_batches,
    const AdaptiveAssemblyOptions &options, int retained_patch_batch = -1);

std::vector<LocalSparseEnrichmentMatrices> AssembleLocalSparseEnrichmentMatricesBatch(
    const TriangleDofTopology &topology, mfem::FiniteElementSpace &h1_fespace,
    mfem::FiniteElementSpace &nd_fespace,
    const std::vector<std::vector<IsotropicMaterialCoefficients>> &material_batches,
    const AdaptiveAssemblyOptions &options, int retained_patch_batch = -1);

// Electrostatic specialization of AssembleLocalSparseEnrichmentMatrices. It
// assembles only the H1 diffusion block and therefore does not evaluate any
// rotational ND basis or curl-curl integral.
LocalSparseH1EnrichmentMatrices AssembleLocalSparseH1EnrichmentMatrices(
    const DofTopology &topology, mfem::FiniteElementSpace &h1_fespace,
    const std::vector<IsotropicMaterialCoefficients> &materials,
    const AdaptiveAssemblyOptions &options);

LocalSparseH1EnrichmentMatrices AssembleLocalSparseH1EnrichmentMatrices(
    const TriangleDofTopology &topology, mfem::FiniteElementSpace &h1_fespace,
    const std::vector<IsotropicMaterialCoefficients> &materials,
    const AdaptiveAssemblyOptions &options);

// Assemble the singular portions of a tangential H(curl) boundary mass
// operator
//
//   integral_Gamma coefficient (n x u) dot (n x v) dS.
//
// Coefficients are indexed by the original positive mesh boundary attribute.
// The standard-standard block remains assembled by Palace's ordinary boundary
// integrator; these routines return only the standard-enrichment and
// enrichment-enrichment blocks needed to complete the combined space.
//
// Boundary quadrature is feature aligned and rejects a selected boundary
// segment or face containing a nonintegrable singular edge. Curved boundary
// and volume transformations are evaluated pointwise.
LocalSparseOperatorBlocks AssembleLocalSparseNDBoundaryMassMatrices(
    const DofTopology &topology, mfem::FiniteElementSpace &nd_fespace,
    const std::map<int, double> &boundary_coefficients,
    const AdaptiveAssemblyOptions &options,
    std::vector<LocalNDBoundaryPatchMatrices> *retained_patches = nullptr,
    const std::set<int> &excluded_singular_attributes = {});

LocalSparseOperatorBlocks AssembleLocalSparseNDBoundaryMassMatrices(
    const TriangleDofTopology &topology, mfem::FiniteElementSpace &nd_fespace,
    const std::map<int, double> &boundary_coefficients,
    const AdaptiveAssemblyOptions &options);

// Assemble the singular portions of a scalar H1 boundary mass operator
//
//   integral_Gamma coefficient u v dS
//
// on a two-dimensional triangular mesh. This is the longitudinal Robin block
// required by BoundaryMode impedance boundaries.
LocalSparseOperatorBlocks AssembleLocalSparseH1BoundaryMassMatrices(
    const TriangleDofTopology &topology, mfem::FiniteElementSpace &h1_fespace,
    const std::map<int, double> &boundary_coefficients,
    const AdaptiveAssemblyOptions &options);

// Form true-DOF parallel blocks as P_test^T A_local P_trial. The custom
// enrichment prolongations are built from the canonical parallel numbering;
// the standard prolongations come from MFEM. Error bounds use entrywise
// absolute standard prolongations. Reverse value blocks are exact parallel
// transposes of the assembled forward blocks.
ParallelSparseEnrichmentMatrices
AssembleParallelSparseEnrichmentMatrices(const LocalSparseEnrichmentMatrices &local,
                                         const ParallelDofNumbering &parallel_numbering,
                                         const mfem::ParFiniteElementSpace &h1_fespace,
                                         const mfem::ParFiniteElementSpace &nd_fespace);

// Interpolate the finite scalar singular potentials into the standard nodal H1
// space. The returned true-DOF map has dimensions
//
//   standard H1 true DOFs x enrichment H1 true DOFs.
//
// Values on the singular feature are evaluated by their continuous zero limit,
// rather than by evaluating the singular gradient formula at rho = 0.
std::unique_ptr<mfem::HypreParMatrix>
BuildParallelH1EnrichmentInterpolant(const DofTopology &topology,
                                     const ParallelDofNumbering &parallel_numbering,
                                     const mfem::ParFiniteElementSpace &h1_fespace);

std::unique_ptr<mfem::HypreParMatrix>
BuildParallelH1EnrichmentInterpolant(const TriangleDofTopology &topology,
                                     const ParallelDofNumbering &parallel_numbering,
                                     const mfem::ParFiniteElementSpace &h1_fespace);

// Interpolate only the finite rotational singular vector bases into the
// standard Nedelec space. Gradient-family columns are exactly zero because
// their interpolants are built from the H1 potential above to preserve the
// discrete gradient relation. The returned true-DOF map has dimensions
//
//   standard ND true DOFs x enrichment ND true DOFs.
std::unique_ptr<mfem::HypreParMatrix> BuildParallelNDRotationalEnrichmentInterpolant(
    const DofTopology &topology, const ParallelDofNumbering &parallel_numbering,
    const mfem::ParFiniteElementSpace &nd_fespace);

std::unique_ptr<mfem::HypreParMatrix> BuildParallelNDRotationalEnrichmentInterpolant(
    const TriangleDofTopology &topology, const ParallelDofNumbering &parallel_numbering,
    const mfem::ParFiniteElementSpace &nd_fespace);

ParallelSparseOperatorBlocks
AssembleParallelSparseH1EnrichmentMatrices(const LocalSparseH1EnrichmentMatrices &local,
                                           const ParallelDofNumbering &parallel_numbering,
                                           const mfem::ParFiniteElementSpace &h1_fespace);

ParallelSparseOperatorBlocks
AssembleParallelSparseNDBoundaryMassMatrices(const LocalSparseOperatorBlocks &local,
                                             const ParallelDofNumbering &parallel_numbering,
                                             const mfem::ParFiniteElementSpace &nd_fespace);

ParallelSparseOperatorBlocks
AssembleParallelSparseH1BoundaryMassMatrices(const LocalSparseOperatorBlocks &local,
                                             const ParallelDofNumbering &parallel_numbering,
                                             const mfem::ParFiniteElementSpace &h1_fespace);

// Restrict a fine-level enrichment block through the standard-space
// prolongation P. The enrichment basis is identical on every polynomial
// level, so the exact combined transfer is diag(P, I):
//
//   A_se,c = P^T A_se,f,  A_es,c = A_se,c^T,  A_ee,c = A_ee,f.
//
// Error bounds use |P|^T. An empty input block produces an empty output block.
ParallelSparseOperatorBlocks
RestrictParallelSparseOperatorBlocks(const ParallelSparseOperatorBlocks &fine,
                                     const mfem::HypreParMatrix &standard_prolongation);

ParallelSparseEnrichmentMatrices RestrictParallelSparseEnrichmentMatrices(
    const ParallelSparseEnrichmentMatrices &fine,
    const mfem::HypreParMatrix &h1_standard_prolongation,
    const mfem::HypreParMatrix &nd_standard_prolongation);

// Project the enrichment-containing blocks of a symmetric ND operator through
// diag(G_standard, G_enrichment). The standard-standard block remains owned by
// Palace's ordinary auxiliary operator.
ParallelSparseOperatorBlocks
ProjectParallelSparseOperatorBlocksToH1(const ParallelSparseOperatorBlocks &nd,
                                        const mfem::HypreParMatrix &standard_gradient,
                                        const mfem::HypreParMatrix &enrichment_gradient);

// Assemble the enrichment entries of
//
//   integral_Gamma coefficient dot v dS
//
// into the process-local owned enrichment true-DOF vector. When supplied, the
// boundary marker avoids evaluating the coefficient and singular traces on
// attributes outside its support. Standard entries continue to use MFEM's
// boundary linear form.
void AssembleParallelNDBoundaryLinearForm(
    const DofTopology &topology, const ParallelDofNumbering &parallel_numbering,
    mfem::ParFiniteElementSpace &nd_fespace, mfem::VectorCoefficient &coefficient,
    const AdaptiveAssemblyOptions &options, mfem::Vector &enrichment_true_dofs,
    const mfem::Array<int> *attribute_marker = nullptr);

void AssembleParallelNDBoundaryLinearForm(
    const TriangleDofTopology &topology, const ParallelDofNumbering &parallel_numbering,
    mfem::ParFiniteElementSpace &nd_fespace, mfem::VectorCoefficient &coefficient,
    const AdaptiveAssemblyOptions &options, mfem::Vector &enrichment_true_dofs,
    const mfem::Array<int> *attribute_marker = nullptr);

// Construct the exact true-DOF map from enriched scalar potentials to their
// matching enriched ND gradient functions. Every H1 column has one coefficient
// +1 in the ND row identified by ParallelDofNumbering::h1_to_nd_true;
// rotational ND rows are exactly zero.
std::unique_ptr<mfem::HypreParMatrix>
BuildParallelEnrichmentGradient(MPI_Comm comm,
                                const ParallelDofNumbering &parallel_numbering);

}  // namespace singular

}  // namespace fem

}  // namespace palace

#endif  // PALACE_FEM_SINGULARASSEMBLY_HPP
