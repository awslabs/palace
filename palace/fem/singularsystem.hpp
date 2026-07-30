// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_FEM_SINGULARSYSTEM_HPP
#define PALACE_FEM_SINGULARSYSTEM_HPP

#include <memory>
#include <mfem.hpp>

#include "fem/singularassembly.hpp"
#include "linalg/vector.hpp"

namespace palace
{

namespace fem
{

namespace singular
{

// Form the monolithic true-DOF operator
//
//   [ A_ss  A_se ]
//   [ A_es  A_ee ]
//
// from an existing standard-standard matrix and the assembled enrichment
// blocks. On each rank, vectors use the process-local ordering
// [standard true DOFs, owned enrichment true DOFs].
std::unique_ptr<mfem::HypreParMatrix>
BuildParallelEnrichedOperator(const mfem::HypreParMatrix &standard_standard,
                              const ParallelSparseOperatorBlocks &enrichment);

// Compact exact off-diagonal zeros left behind by essential-DOF elimination.
// Diagonal storage is always retained, and no tolerance is applied.
HYPRE_BigInt RemoveExplicitZeros(mfem::HypreParMatrix &matrix);

// Form diag(G_standard, G_enrichment), the exact true-DOF gradient from the
// combined standard-plus-enriched H1 space to the corresponding ND space.
std::unique_ptr<mfem::HypreParMatrix>
BuildParallelEnrichedGradient(const mfem::HypreParMatrix &standard_gradient,
                              const mfem::HypreParMatrix &enrichment_gradient);

// Form diag(P_standard, I_enrichment), the true-DOF prolongation between two
// polynomial levels of a combined standard-plus-enriched space. The singular
// basis and its canonical true DOFs are identical on every p-level.
std::unique_ptr<mfem::HypreParMatrix>
BuildParallelEnrichedProlongation(const mfem::HypreParMatrix &standard_prolongation,
                                  const TrueDofMap &enrichment_numbering);

// Independently constrain the four blocks of an enriched SPD operator. The
// resulting blocks have the same action as symmetric row/column elimination
// applied after monolithic block assembly.
struct ParallelConstrainedOperatorBlocks
{
  std::unique_ptr<mfem::HypreParMatrix> standard_standard;
  std::unique_ptr<mfem::HypreParMatrix> standard_enrichment;
  std::unique_ptr<mfem::HypreParMatrix> enrichment_standard;
  std::unique_ptr<mfem::HypreParMatrix> enrichment_enrichment;
};

ParallelConstrainedOperatorBlocks BuildParallelConstrainedOperatorBlocks(
    const mfem::HypreParMatrix &standard_standard,
    const ParallelSparseOperatorBlocks &enrichment,
    const mfem::Array<int> &standard_essential_true_dofs,
    const mfem::Array<int> &enrichment_essential_true_dofs);

// Principal submatrix containing every enrichment true DOF and every standard
// true DOF with nonzero standard-enrichment coupling. The process-local index
// order is retained in true_dofs and is the restriction/prolongation order of
// matrix.
struct ParallelCoupledPatch
{
  std::unique_ptr<mfem::HypreParMatrix> matrix;
  mfem::Array<int> true_dofs;
  HYPRE_BigInt global_standard_dofs;
  HYPRE_BigInt global_enrichment_dofs;
};

ParallelCoupledPatch
BuildParallelCoupledPatch(const mfem::HypreParMatrix &constrained,
                          const mfem::HypreParMatrix &standard_enrichment,
                          int standard_local_size);

// Overlapping principal patches associated with straight singular features.
// Each patch contains every enrichment true DOF whose singular node or edge
// belongs to the feature and every standard true DOF directly coupled to one
// of those enrichment DOFs. Node functions at polygon corners intentionally
// occur in each incident feature patch.
struct ParallelFeaturePatch
{
  std::size_t feature;
  std::unique_ptr<mfem::HypreParMatrix> matrix;
  mfem::Array<int> true_dofs;
  HYPRE_BigInt global_standard_dofs;
  HYPRE_BigInt global_enrichment_dofs;
};

struct ParallelFeaturePatches
{
  std::vector<ParallelFeaturePatch> patches;
  HYPRE_BigInt sum_global_standard_dofs;
  HYPRE_BigInt sum_global_enrichment_dofs;
  HYPRE_BigInt maximum_global_standard_dofs;
  HYPRE_BigInt maximum_global_enrichment_dofs;
  int minimum_enrichment_multiplicity;
  int maximum_enrichment_multiplicity;
};

ParallelFeaturePatches BuildParallelFeaturePatches(
    const mfem::HypreParMatrix &constrained,
    const mfem::HypreParMatrix &standard_enrichment, int standard_local_size,
    const std::vector<std::vector<std::size_t>> &local_enrichment_features,
    const TrueDofMap &enrichment_numbering, std::size_t number_features);

// An invertible symmetric diagonal change of coordinates for an SPD system:
//
//   A_s = S A S,  b_s = S b,  x = S y,
//
// where S_ii = A_ii^(-1/2). The transformed matrix has a unit diagonal. This
// changes only the algebraic coordinates, not the represented finite-element
// function.
class SymmetricDiagonalScaling
{
private:
  std::unique_ptr<mfem::HypreParMatrix> matrix;
  Vector coordinate_scaling;
  double scaled_diagonal_minimum;
  double scaled_diagonal_maximum;

  void Apply(const Vector &input, Vector &output, bool inverse) const;

public:
  explicit SymmetricDiagonalScaling(const mfem::HypreParMatrix &unscaled);

  const mfem::HypreParMatrix &GetMatrix() const { return *matrix; }
  const Vector &GetCoordinateScaling() const { return coordinate_scaling; }
  double GetScaledDiagonalMinimum() const { return scaled_diagonal_minimum; }
  double GetScaledDiagonalMaximum() const { return scaled_diagonal_maximum; }

  // Transform a physical right-hand side b into b_s = S b.
  void ScaleRHS(const Vector &unscaled, Vector &scaled) const;

  // Transform a physical initial guess x into y = S^(-1) x.
  void ScaleInitialGuess(const Vector &unscaled, Vector &scaled) const;

  // Recover physical coefficients x = S y from transformed coordinates.
  void RecoverSolution(const Vector &scaled, Vector &unscaled) const;
};

struct ParallelDirichletSystem
{
  // The constrained matrix and the entries removed during symmetric
  // row/column elimination. Their sum is the original unconstrained matrix.
  std::unique_ptr<mfem::HypreParMatrix> constrained;
  std::unique_ptr<mfem::HypreParMatrix> eliminated;

  // Process-local indices in the combined true-DOF vector.
  mfem::Array<int> essential_true_dofs;

  // Apply nonzero essential values in x to an existing unconstrained right
  // hand side b.
  void EliminateRHS(const mfem::Vector &x, mfem::Vector &b) const;
};

// Symmetrically eliminate standard and enrichment essential true DOFs from a
// combined operator. Enrichment indices are local to the owned enrichment
// block and are shifted by standard_local_size.
ParallelDirichletSystem
BuildParallelDirichletSystem(std::unique_ptr<mfem::HypreParMatrix> &&matrix,
                             int standard_local_size,
                             const mfem::Array<int> &standard_essential_true_dofs,
                             const mfem::Array<int> &enrichment_essential_true_dofs);

}  // namespace singular

}  // namespace fem

}  // namespace palace

#endif  // PALACE_FEM_SINGULARSYSTEM_HPP
