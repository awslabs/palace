// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_LINALG_AGGREGATE_BLOCKS_HPP
#define PALACE_LINALG_AGGREGATE_BLOCKS_HPP

#include <memory>
#include <vector>
#include <mfem.hpp>
#include "linalg/operator.hpp"
#include "linalg/vector.hpp"

namespace palace
{

class ParOperator;

//
// Blocks of the block Jacobi smoother which follow the vertex aggregates of the mesh (see
// linalg/vertexaggregates.hpp), for a real symmetric operator P ᵀ A P on a finite element
// space of a 3D mesh with a partially assembled A (or with an assembled local matrix A).
//
//  - The local DOFs of the mesh entities with a vertex in an aggregate are grouped by the
//    set of aggregates which the vertices of their entity belong to: over a stack of flat
//    elements, the DOFs of the entities above each other form a "column". The true DOFs of
//    a column form a block, and the blocks are merged where the conforming prolongation of
//    a hanging entity makes the functions which the strong part of the operator does not
//    see in two blocks depend on each other. All other DOFs are in no block.
//  - The matrix of a block is the exact diagonal block of P ᵀ A P: the entries of the
//    element matrices of A between DOFs of a common block are summed. Element matrices are
//    obtained by element-local applications of the libCEED operator, and nothing outside
//    of the blocks is assembled.
//  - A block has at most a given number of DOFs, a power of two up to 64. The true DOFs
//    of a larger group (an aggregate which is a sheet of flat elements without structure)
//    which a process owns are split there into connected pieces of at most that size,
//    cutting first the pairs of DOFs which are not tied across the flat direction and then
//    the weakest ties. A group which fits and has DOFs on more than one process is
//    gathered on one of them.
//  - The blocks are dense and are stored in size classes of 1, 2, 4, ..., 64 DOFs: a block
//    between two classes is padded into the upper one. What is stored is the inverse of
//    each block, from a symmetric eigenvalue decomposition: the pseudo-inverse where the
//    block is singular (the blocks of a curl-curl operator without a mass term have the
//    gradients supported in the block as kernel). The inverses of a class are applied by
//    one kernel, a product with a dense symmetric matrix for every block.
//  - For an operator which is not definite, the blocks are those of its definite
//    counterpart with the signs of the operator itself: with M the block of the
//    counterpart, A the block of the operator, and A u = μ M u with Uᵀ M U = I, the
//    inverse is U sign(μ) Uᵀ instead of M⁻¹ = U Uᵀ. The relaxed operator then has the
//    eigenvalues |μ| ≤ 1 inside a block: a gradient inside a block of a curl-curl
//    operator with a negative mass term is relaxed with its true sign, and a direction in
//    which the block of the operator is close to singular is left alone instead of
//    being amplified. This is the block form of the sign of a negative diagonal entry in
//    the point Jacobi method.
//  - On a nodal space, a group which is too large for a block is not split: its DOFs stay
//    on their processes and are solved approximately as a whole, with a fixed Chebyshev
//    polynomial of the matrix of the group preconditioned by one AMG V-cycle, as a
//    separate step of a smoothing sweep. The degree targets a reduction on an estimated
//    interval and is bounded (see SheetDegree and GetSheetLimits). The block size bounds
//    the dense blocks, not the size of such a group.
//
// The blocks are constructed and applied on the host.
//
class AggregateBlocks
{
private:
  struct Impl;
  std::unique_ptr<Impl> impl;

public:
  // Number of size classes: blocks of 1, 2, 4, ..., 64 DOFs.
  static constexpr int NUM_CLASSES = 7;
  static constexpr int MAX_BLOCK_SIZE = 1 << (NUM_CLASSES - 1);

  // What became of the groups of at least two DOFs: global numbers.
  struct Census
  {
    HYPRE_BigInt blocks[NUM_CLASSES] = {};    // Blocks of every size class,
    HYPRE_BigInt block_dofs = 0;              // their DOFs,
    HYPRE_BigInt entries = 0;                 // and stored entries (padded, one triangle),
    HYPRE_BigInt max_entries = 0;             // in all and on the process with the most
    HYPRE_BigInt split = 0;                   // Groups which were split into blocks
    HYPRE_BigInt gathered = 0;                // Blocks with DOFs of more than one process
    HYPRE_BigInt dropped = 0;                 // Singular directions of the blocks
    HYPRE_BigInt negative = 0;                // Directions with the negative sign
    HYPRE_BigInt sheets = 0, sheet_dofs = 0;  // Groups solved approximately as a whole,
    HYPRE_BigInt sheet_entries = 0;           // the entries of their AMG hierarchy,
    int sheet_degree = 0;                     // the degree of the polynomial,
    bool sheet_limited = false;               // and whether that is its limit
    HYPRE_BigInt skipped = 0, skipped_dofs = 0;  // Groups left to the point Jacobi method
  };

  // Blocks of at most block_size DOFs of the real symmetric operator A. If A is not
  // definite, A_def is its definite counterpart (otherwise nullptr). The approximate
  // solves of the large groups of a nodal space reduce the error by at least
  // sheet_reduction.
  AggregateBlocks(const ParOperator &A, const ParOperator *A_def, int block_size,
                  double sheet_reduction);
  ~AggregateBlocks();

  // Whether the blocks can be constructed for the operator (on a 3D mesh).
  static bool IsSupported(const ParOperator &A);

  // The smallest degree d for which the Chebyshev iteration on an interval [a, b] of
  // positive eigenvalues reduces the error by at least the given factor:
  // 2 ρᵈ / (1 + ρ²ᵈ) ≤ reduction, with ρ = (√κ - 1) / (√κ + 1) and κ = b / a.
  static int SheetDegree(double a, double b, double reduction);

  // The limits of the work of the approximate solves for a required reduction r: the
  // largest condition number κ = 1 / r² for which the degree is the one of the estimated
  // interval, the degree for it, which is the largest one used, and the largest number of
  // steps of the Lanczos run which estimates the interval, after which the error bound
  // of the conjugate gradient method for that κ is at the rounding level.
  struct SheetLimits
  {
    double condition;
    int degree, lanczos_steps;
  };
  static SheetLimits GetSheetLimits(double reduction);

  const Census &GetCensus() const;

  // The local true DOFs of the blocks.
  const std::vector<int> &GetBlockDofs() const;

  // z = M⁻¹ x on the DOFs of the blocks, with M⁻¹ the (signed, pseudo-) inverses of the
  // blocks. Other entries of z are left untouched.
  void MultInverse(const Vector &x, Vector &z) const;

  // The groups which are solved approximately are not part of M⁻¹ (their DOFs keep the
  // point Jacobi scaling in it). Their solves are a separate step y += S r on their DOFs,
  // where S is the fixed polynomial approximation of the inverse of their matrix, a
  // symmetric positive definite operator.
  bool HasSheets() const;
  const std::vector<int> &GetSheetDofs() const;
  void AddSheetSolve(const Vector &r, Vector &y) const;
};

}  // namespace palace

#endif  // PALACE_LINALG_AGGREGATE_BLOCKS_HPP
