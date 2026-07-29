// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_FEM_SINGULARASSEMBLY_HPP
#define PALACE_FEM_SINGULARASSEMBLY_HPP

#include <cstddef>
#include <map>
#include <memory>
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
};

struct LocalSparseEnrichmentMatrices
{
  LocalSparseOperatorBlocks h1_diffusion;
  LocalSparseOperatorBlocks h1_mass;
  LocalSparseOperatorBlocks nd_mass;
  LocalSparseOperatorBlocks nd_curl_curl;

  std::size_t total_quadrature_leaf_count = 0;
  int maximum_subdivision_depth = 0;
};

struct LocalSparseH1EnrichmentMatrices
{
  LocalSparseOperatorBlocks diffusion;

  std::size_t total_quadrature_leaf_count = 0;
  int maximum_subdivision_depth = 0;
  std::size_t duffy_reference_table_entries = 0;
  std::size_t duffy_reference_cache_hits = 0;
};

struct ParallelSparseOperatorBlocks
{
  std::unique_ptr<mfem::HypreParMatrix> enrichment_enrichment;
  std::unique_ptr<mfem::HypreParMatrix> standard_enrichment;
  std::unique_ptr<mfem::HypreParMatrix> enrichment_standard;

  std::unique_ptr<mfem::HypreParMatrix> enrichment_enrichment_estimated_absolute_error;
  std::unique_ptr<mfem::HypreParMatrix> standard_enrichment_estimated_absolute_error;
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
    const AdaptiveAssemblyOptions &options);

LocalSparseOperatorBlocks AssembleLocalSparseNDBoundaryMassMatrices(
    const TriangleDofTopology &topology, mfem::FiniteElementSpace &nd_fespace,
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

ParallelSparseOperatorBlocks
AssembleParallelSparseH1EnrichmentMatrices(const LocalSparseH1EnrichmentMatrices &local,
                                           const ParallelDofNumbering &parallel_numbering,
                                           const mfem::ParFiniteElementSpace &h1_fespace);

ParallelSparseOperatorBlocks
AssembleParallelSparseNDBoundaryMassMatrices(const LocalSparseOperatorBlocks &local,
                                             const ParallelDofNumbering &parallel_numbering,
                                             const mfem::ParFiniteElementSpace &nd_fespace);

// Assemble the enrichment entries of
//
//   integral_Gamma coefficient dot v dS
//
// into the process-local owned enrichment true-DOF vector. The coefficient
// itself is responsible for restricting its support to the desired boundary
// attributes. Standard entries continue to use MFEM's boundary linear form.
void AssembleParallelNDBoundaryLinearForm(const DofTopology &topology,
                                          const ParallelDofNumbering &parallel_numbering,
                                          mfem::ParFiniteElementSpace &nd_fespace,
                                          mfem::VectorCoefficient &coefficient,
                                          const AdaptiveAssemblyOptions &options,
                                          mfem::Vector &enrichment_true_dofs);

void AssembleParallelNDBoundaryLinearForm(const TriangleDofTopology &topology,
                                          const ParallelDofNumbering &parallel_numbering,
                                          mfem::ParFiniteElementSpace &nd_fespace,
                                          mfem::VectorCoefficient &coefficient,
                                          const AdaptiveAssemblyOptions &options,
                                          mfem::Vector &enrichment_true_dofs);

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
