// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_LINALG_ENTITY_JACOBI_HPP
#define PALACE_LINALG_ENTITY_JACOBI_HPP

#include <array>
#include <vector>
#include "linalg/operator.hpp"
#include "linalg/solver.hpp"
#include "linalg/vector.hpp"

namespace palace
{

class FiniteElementSpace;

namespace ceed
{

class Operator;

}  // namespace ceed

//
// Block Jacobi preconditioner for a real SPD operator A = Pᵀ A_loc P on an H(div) or
// H(curl) space, with one block per mesh entity: the true DOFs which the process owns on
// each edge, face, and element interior (separately for each copy of a broken space) form a
// dense block of A. For high-order nodal bases, the basis functions of an entity are close
// to linearly dependent on anisotropic elements, which the Jacobi diagonal cannot correct.
// The blocks are those of the assembled A, computed from the entries of the libCEED
// element matrices of A_loc which can contribute to them (including through the
// constraints of a nonconforming mesh and from the elements of other processes). A true
// DOF without an unconstrained local L-DOF forms a block of its own.
//
template <OperatorType OperType>
class EntityBlockJacobiSolver : public Solver<OperType>
{
  using VecType = typename Solver<OperType>::VecType;

private:
  // True DOFs of each block (offsets block_ptr) and the inverse blocks (row-major, offsets
  // inv_ptr).
  std::vector<int> block_ptr, block_dofs;
  std::vector<std::size_t> inv_ptr;
  std::vector<double> inv;

  // Storage of the inverse blocks used by the applications: full blocks (inv), or their
  // upper triangles (row i holds the columns j >= i, offsets sym_ptr).
  bool upper = false;
  std::vector<std::size_t> sym_ptr;

  template <int NV>
  void Apply(const std::array<const double *, NV> &X,
             const std::array<double *, NV> &Y) const;

  // Inverse diagonal of the operator, as for Jacobi.
  Vector dinv;

public:
  // Construct from the process-local libCEED operator A_loc on the given space. Collective.
  EntityBlockJacobiSolver(const ceed::Operator &op, const FiniteElementSpace &fespace);

  // Whether the preconditioner is used for the given space: host execution without OpenMP,
  // and an H(div) or H(curl) space of order two or three on a mesh of tetrahedra and prisms
  // (in any combination) or of triangles, or of order four on a mesh of only tetrahedra or
  // only triangles. The result is the same on all processes.
  static bool Supported(const FiniteElementSpace &fespace);

  // The true DOFs of block b are block_dofs[block_ptr[b], block_ptr[b + 1]).
  const auto &GetBlockOffsets() const { return block_ptr; }
  const auto &GetBlockDofs() const { return block_dofs; }

  // The inverse diagonal D⁻¹ of the operator (its real part, for complex operators), as
  // used by JacobiSmoother, for the stopping norm of CgSolver (set by SetOperator).
  const Vector &GetInverseDiagonal() const { return dinv; }

  void SetOperator(const OperType &op) override;

  void Mult(const VecType &x, VecType &y) const override;

  void MultTranspose(const VecType &x, VecType &y) const override { Mult(x, y); }
};

}  // namespace palace

#endif  // PALACE_LINALG_ENTITY_JACOBI_HPP
