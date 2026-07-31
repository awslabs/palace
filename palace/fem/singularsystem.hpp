// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_FEM_SINGULARSYSTEM_HPP
#define PALACE_FEM_SINGULARSYSTEM_HPP

#include <memory>
#include <mfem.hpp>

#include "fem/singularassembly.hpp"
#include "linalg/rap.hpp"
#include "linalg/vector.hpp"

namespace palace
{

namespace fem
{

namespace singular
{

// Symmetric restricted-additive-Schwarz inverse on enriched elements. Each
// dense patch couples all standard and singular ND basis functions on one
// element. True-to-local prolongations provide MPI exchange and the standard
// Nedelec orientation; diagonal overlap weights are applied on both sides.
class ParallelElementPatchInverse : public Operator
{
private:
  struct Patch
  {
    mfem::Array<int> local_dofs;
    mfem::Vector signs;
    mfem::DenseMatrix inverse;
    mutable Vector rhs, correction;
  };

  const mfem::HypreParMatrix *standard_prolongation;
  std::unique_ptr<mfem::HypreParMatrix> enrichment_prolongation;
  int standard_true_size;
  int standard_local_size;
  int enrichment_local_size;
  std::vector<Patch> patches;
  Vector true_weight;
  mutable Vector scaled_input, standard_local_rhs, enrichment_local_rhs,
      standard_local_correction, enrichment_local_correction;

  static int DecodeSignedDof(int dof);
  static double DecodeSignedDofSign(int dof);
  void ApplyTrueWeight(const Vector &input, Vector &output) const;

public:
  ParallelElementPatchInverse(
      const mfem::ParFiniteElementSpace &standard_fespace,
      const TrueDofMap &enrichment_numbering,
      const std::vector<LocalNDElementPatchMatrices> &element_matrices,
      double stiffness_coefficient, double mass_coefficient,
      const mfem::Array<int> &standard_essential_true_dofs,
      const mfem::Array<int> &enrichment_essential_true_dofs);

  std::size_t GetNumPatches() const { return patches.size(); }
  void Mult(const Vector &input, Vector &output) const override;
  void MultTranspose(const Vector &input, Vector &output) const override
  {
    Mult(input, output);
  }
};

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

// Matrix-free standard-standard block plus sparse enrichment couplings:
//
//   [ A_ss  A_se ]
//   [ A_es  A_ee ].
//
// The standard block can retain libCEED partial assembly. Essential rows and
// columns are eliminated independently in each block, exactly matching
// monolithic symmetric elimination.
class ParallelHybridEnrichedOperator : public Operator
{
private:
  std::unique_ptr<Operator> standard_standard;
  ParallelSparseOperatorBlocks enrichment;
  std::shared_ptr<const Operator> coupled_patch_inverse;
  int standard_size;

  void MakeInputBlocks(const Vector &input, Vector &standard, Vector &enriched) const;
  void MakeOutputBlocks(Vector &output, Vector &standard, Vector &enriched) const;

public:
  ParallelHybridEnrichedOperator(
      std::unique_ptr<Operator> &&standard_standard,
      const ParallelSparseOperatorBlocks &enrichment,
      const mfem::Array<int> &standard_essential_true_dofs,
      const mfem::Array<int> &enrichment_essential_true_dofs,
      Operator::DiagonalPolicy diagonal_policy,
      std::shared_ptr<const Operator> coupled_patch_inverse = {});
  ParallelHybridEnrichedOperator(
      std::unique_ptr<Operator> &&standard_standard,
      ParallelSparseOperatorBlocks &&enrichment,
      const mfem::Array<int> &standard_essential_true_dofs,
      const mfem::Array<int> &enrichment_essential_true_dofs,
      Operator::DiagonalPolicy diagonal_policy,
      std::shared_ptr<const Operator> coupled_patch_inverse = {});

  const Operator &GetStandardStandard() const { return *standard_standard; }
  const ParallelSparseOperatorBlocks &GetEnrichmentBlocks() const { return enrichment; }
  const Operator *GetCoupledPatchInverse() const { return coupled_patch_inverse.get(); }
  int GetStandardSize() const { return standard_size; }

  void AssembleDiagonal(Vector &diagonal) const override;
  void Mult(const Vector &input, Vector &output) const override;
  void MultTranspose(const Vector &input, Vector &output) const override;
  void AddMult(const Vector &input, Vector &output,
               double coefficient = 1.0) const override;
  void AddMultTranspose(const Vector &input, Vector &output,
                        double coefficient = 1.0) const override;
};

// Congruence-transformed hybrid operator
//
//   A_s = T^T A T,  T = [ I  -D ]
//                          [ 0   I ],
//
// where D is a conforming standard-space interpolant of the high-order part of
// the gradient enrichment. This is an invertible change of basis: it preserves
// the represented H(curl) space while reducing the angle between the standard
// and enrichment subspaces seen by p-multigrid.
class ParallelTransformedHybridEnrichedOperator : public Operator
{
private:
  std::unique_ptr<ParallelHybridEnrichedOperator> untransformed;
  const mfem::HypreParMatrix *coordinate_shift;
  int standard_size;
  mutable Vector transformed_input, transformed_action, standard_work, enrichment_work;

  void ApplyCoordinateTransform(const Vector &input, Vector &output) const;
  void ApplyCoordinateTransformTranspose(const Vector &input, Vector &output) const;

public:
  ParallelTransformedHybridEnrichedOperator(
      std::unique_ptr<ParallelHybridEnrichedOperator> &&untransformed,
      const mfem::HypreParMatrix &coordinate_shift);

  const ParallelHybridEnrichedOperator &GetUntransformedOperator() const
  {
    return *untransformed;
  }
  const mfem::HypreParMatrix &GetCoordinateShift() const { return *coordinate_shift; }
  int GetStandardSize() const { return standard_size; }

  void Mult(const Vector &input, Vector &output) const override;
  void MultTranspose(const Vector &input, Vector &output) const override;
  void AddMult(const Vector &input, Vector &output,
               double coefficient = 1.0) const override;
  void AddMultTranspose(const Vector &input, Vector &output,
                        double coefficient = 1.0) const override;
};

// Enrichment-enrichment principal block of a transformed hybrid operator:
//
//   A'_ee = A_ee - D^T A_se - A_es D + D^T A_ss D.
//
// Its action is extracted matrix-free from the complete transformed operator,
// preserving partial assembly in A_ss. AssembleDiagonal uses an exact
// elementwise transformed diagonal when all active components provide one,
// with diag(A_ee) retained as an explicit fallback.
class ParallelTransformedEnrichmentOperator : public Operator
{
private:
  const ParallelTransformedHybridEnrichedOperator *transformed;
  int standard_size;
  mutable Vector combined_input, combined_action;

public:
  explicit ParallelTransformedEnrichmentOperator(
      const ParallelTransformedHybridEnrichedOperator &transformed);

  void AssembleDiagonal(Vector &diagonal) const override;
  void Mult(const Vector &input, Vector &output) const override;
  void MultTranspose(const Vector &input, Vector &output) const override;
};

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

// Form the stabilized-GFEM transfer
//
//   [ P_standard  D_fine - P_standard D_coarse ]
//   [      0                    I                ].
//
// The supplied correction is the already-formed top-right block. Together
// with the matching level coordinate transformations, this transfer represents
// exactly the same physical coarse functions and preserves Galerkin nesting.
std::unique_ptr<mfem::HypreParMatrix> BuildParallelEnrichedProlongation(
    const mfem::HypreParMatrix &standard_prolongation,
    const mfem::HypreParMatrix &standard_enrichment_correction,
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
