// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_LINALG_BLOCK_PRECOND_HPP
#define PALACE_LINALG_BLOCK_PRECOND_HPP

#include <memory>
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

}  // namespace palace

#endif  // PALACE_LINALG_BLOCK_PRECOND_HPP
