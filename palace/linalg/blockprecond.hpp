// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_LINALG_BLOCK_PRECOND_HPP
#define PALACE_LINALG_BLOCK_PRECOND_HPP

#include <memory>
#include <vector>
#include <mfem.hpp>
#include "linalg/operator.hpp"
#include "linalg/solver.hpp"
#include "linalg/vector.hpp"

namespace palace
{

//
// Block lower-triangular preconditioner for a 2-block system:
//
//   P = [P0      0 ]
//       [L10    P1 ]
//
// where P0 and P1 are sub-solvers for the diagonal blocks, and L10 is the lower
// off-diagonal operator (not owned). When L10 is null, this reduces to a block-diagonal
// preconditioner.
//
// Application (forward solve of P z = r):
//   z0 = P0^{-1} r0
//   z1 = P1^{-1} (r1 - L10 z0)
//
// The input/output vectors are monolithic (size = block0_size + block1_size).
//
template <typename OperType>
class BlockDiagonalPreconditioner : public Solver<OperType>
{
  using VecType = typename Solver<OperType>::VecType;

private:
  int block0_size;
  std::unique_ptr<Solver<OperType>> pc0, pc1;
  const OperType *L10 = nullptr;  // Lower off-diagonal (not owned), may be null.
  mutable VecType x0, y0, x1, y1, t1;

  // Copy sub-vector from src[offset..offset+size) to dst[0..size) for each component.
  static void ExtractSubVector(const VecType &src, VecType &dst, int offset, int size);
  static void InsertSubVector(const VecType &src, VecType &dst, int offset, int size);

public:
  BlockDiagonalPreconditioner(int block0_size, std::unique_ptr<Solver<OperType>> &&pc0,
                              std::unique_ptr<Solver<OperType>> &&pc1);

  // Set operators on each block's sub-solver independently.
  void SetBlockOperators(const OperType &op0, const OperType &op1);

  // Set the lower off-diagonal operator L10 for block-triangular preconditioning.
  // When set, the preconditioner applies: z0 = P0^{-1} r0, z1 = P1^{-1}(r1 - L10*z0).
  void SetOffDiagonalOperator(const OperType *op10) { L10 = op10; }

  // SetOperator for the monolithic system (no-op; use SetBlockOperators instead).
  void SetOperator(const OperType &op) override;

  void Mult(const VecType &x, VecType &y) const override;
};

using BlockDiagonalPreconditionerReal = BlockDiagonalPreconditioner<Operator>;
using BlockDiagonalPreconditionerComplex = BlockDiagonalPreconditioner<ComplexOperator>;

//
// Additive correction on overlapping principal patches:
//
//   B_p = sum_i R_i^T B_i R_i.
//
// Every B_i must be symmetric positive definite for B_p to be symmetric
// positive semidefinite. Positive definiteness on a target subspace requires
// the caller to verify that the patches cover that subspace.
//
class AdditivePatchSolver : public mfem::Solver
{
private:
  struct Patch
  {
    mfem::Array<int> dofs;
    std::unique_ptr<mfem::Solver> solver;
    mutable Vector rhs, correction;
  };

  std::vector<Patch> patches;

  static void Restrict(const Vector &source, Vector &destination,
                       const mfem::Array<int> &indices);
  static void AddProlongation(const Vector &source, Vector &destination,
                              const mfem::Array<int> &indices);

public:
  AdditivePatchSolver(int size, std::vector<mfem::Array<int>> &&patch_dofs,
                      std::vector<std::unique_ptr<mfem::Solver>> &&patch_solvers);

  void SetPatchOperators(const std::vector<const Operator *> &patch_operators);
  void SetOperator(const Operator &full_operator) override;
  void Mult(const Vector &x, Vector &y) const override;
};

//
// Symmetric multiplicative subspace correction for an enriched SPD system.
// The first subspace is the complete standard block, while the second is a
// symmetric patch correction which is positive definite on the enrichment
// complement. One application is
//
//   y = B_s r + (I - B_s A) B_p (I - A B_s) r.
//
// If B_s is SPD on the standard subspace and B_p is symmetric positive
// semidefinite with positive action on every enrichment coordinate, this
// composition is SPD on the full space. The left factor is the exact adjoint
// of the right residual correction; no projection identity is assumed for the
// approximate standard solver.
//
class SymmetricPatchSubspacePreconditioner : public Solver<Operator>
{
private:
  int standard_size;
  std::unique_ptr<mfem::Solver> standard_pc, patch_pc;
  const Operator *op = nullptr;
  mutable Vector residual, action, standard_rhs, standard_correction, patch_correction;

  void AddStandardCorrection(const Vector &source, Vector &destination,
                             double coefficient = 1.0) const;
  void UpdateResidual(const Vector &rhs, const Vector &solution) const;

public:
  SymmetricPatchSubspacePreconditioner(int standard_size,
                                       std::unique_ptr<mfem::Solver> &&standard_pc,
                                       std::unique_ptr<mfem::Solver> &&patch_pc);

  void SetSubspaceOperators(const Operator &full_operator,
                            const Operator &standard_operator);
  void SetOperator(const Operator &full_operator) override;
  void Mult(const Vector &x, Vector &y) const override;
};

}  // namespace palace

#endif  // PALACE_LINALG_BLOCK_PRECOND_HPP
