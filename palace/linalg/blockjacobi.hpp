// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_LINALG_BLOCK_JACOBI_SMOOTHER_HPP
#define PALACE_LINALG_BLOCK_JACOBI_SMOOTHER_HPP

#include <complex>
#include <memory>
#include <string>
#include "linalg/aggregateblocks.hpp"
#include "linalg/jacobi.hpp"
#include "linalg/operator.hpp"
#include "linalg/vector.hpp"

namespace palace
{

//
// Block Jacobi smoother for a ParOperator or ComplexParOperator: the point Jacobi method
// where the diagonal of the operator is replaced by dense diagonal blocks of at most
// block_size DOFs, grouped along the directions in which the mesh elements are flat (see
// linalg/aggregateblocks.hpp), which makes the smoother robust in the aspect ratio of the
// elements. The blocks are constructed from the operator itself, which provides its
// element-local operator, its finite element space, the prolongation from the true DOFs,
// and its definite counterpart if it is not definite. True DOFs in no block keep the point
// Jacobi scaling, which is set to zero on the DOFs of the blocks: there the block product
// replaces it. Without blocks (on a 2D mesh, or on a mesh without flat elements) the
// smoother is the point Jacobi smoother.
//
// The blocks are real. For a complex symmetric operator they are those of its real part
// (with the signs of the real part, see linalg/aggregateblocks.hpp), applied between the
// diagonal matrix Ψ of the phases which the diagonal D of the operator has on top of the
// sign of its real part, Ψ² = sign(Re D) D̅ / |D|: M⁻¹ = Ψ M_r⁻¹ Ψ, which for a block of
// one DOF is the phase of the point Jacobi method.
//
// The groups of a nodal space which are too large for a block are solved approximately
// as a whole, once per smoothing sweep (the sweep correction of the relaxation), with the
// phase and sign of the diagonal of the operator on their DOFs.
//
// The blocks are applied on the host. The smoother is not a JacobiSmoother for the fused
// update of CgSolver (which checks for the exact type): the scaling it inherits is the
// one of the true DOFs in no block, and zero on the others.
//
template <OperatorType OperType>
class BlockJacobiSmoother : public JacobiSmoother<OperType>
{
  using VecType = typename Solver<OperType>::VecType;
  using JacobiSmoother<OperType>::comm;
  using JacobiSmoother<OperType>::dinv;

private:
  // The blocks for the operator (empty if the operator has none), the largest size of a
  // block, and the reduction asked of the approximate solves of the large nodal groups.
  std::unique_ptr<AggregateBlocks> blocks;
  const int block_size;
  const double sheet_reduction;

  // Phases (complex operators) and sign (real operators) for the blocks and the sheets,
  // and work vectors.
  std::vector<std::complex<double>> psi, psi_sheet;
  double sign_sheet = 1.0;
  mutable Vector wr, wi, zr, zi;

  // Whether to print what became of the blocks when the operator is set, and what was
  // printed last: an operator with the same result (the next frequency) prints nothing.
  const bool print;
  std::string summary;

  // d += s M⁻¹ r on the DOFs of the blocks, where the point scaling is zero.
  void AddBlockMult(double s, const VecType &r, VecType &d) const;

public:
  // The approximate solves of the large nodal groups reduce the error by 1 / (2k + 1) for
  // the Chebyshev order k of the smoother which uses this relaxation (see
  // linalg/aggregateblocks.cpp).
  BlockJacobiSmoother(MPI_Comm comm, int block_size, int cheby_order, bool print = false);
  ~BlockJacobiSmoother() override;

  void SetOperator(const OperType &op) override;

  void Mult(const VecType &x, VecType &y) const override;

  void MultTranspose(const VecType &x, VecType &y) const override { Mult(x, y); }

  void ScaledMult(double s, const VecType &r, VecType &d) const override;

  void ScaledAddMult(double sd, double sr, const VecType &r, VecType &d) const override;

  double EstimateMaxEigenvalue(const OperType &A, VecType &x0) const override;

  bool HasSweepCorrection() const override;

  void AddSweepCorrection(const VecType &r, VecType &y) const override;

  // Whether the operator has blocks or sheets: without, the smoother is point Jacobi.
  bool HasBlocks() const { return blocks != nullptr; }

  const AggregateBlocks *GetBlocks() const { return blocks.get(); }
};

}  // namespace palace

#endif  // PALACE_LINALG_BLOCK_JACOBI_SMOOTHER_HPP
