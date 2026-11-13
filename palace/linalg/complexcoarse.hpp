// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_LINALG_COMPLEX_COARSE_HPP
#define PALACE_LINALG_COMPLEX_COARSE_HPP

#include <memory>
#include <vector>
#include <mpi.h>
#include "linalg/operator.hpp"
#include "linalg/solver.hpp"
#include "linalg/vector.hpp"

namespace palace
{

class FiniteElementSpace;

// A vector of entries ±1 for the true DOFs of a finite element space which, when the space
// is first order, is the same function on the mesh for every partition of it (at higher
// orders, the basis of the face degrees of freedom of tetrahedra depends on the face
// orientations of the partition): the signs of the degrees of freedom of a pseudo-random
// field, a hash of the position on a grid of 2⁻²⁶ (mesh coordinates are
// nondimensionalized), far below the mesh size and above the round-off of the coordinates,
// projected on each element and averaged over the elements of each degree of freedom. The
// signs follow the orientations of the degrees of freedom, as an operator does. For a
// broken space, every version of a true DOF takes the value of the true DOF of the
// underlying space.
Vector MeshNoise(const FiniteElementSpace &fespace);

//
// Sparse direct coarse solver for a complex-valued operator which chooses at setup between
// the real-valued approximation B = (Ar + Ai)⁻¹ and the exact complex-valued matrix
// Ar + i Ai factored in complex arithmetic (Solver.Linear.ComplexCoarseSolve omitted).
//
// The real-valued solve is exact where Ai vanishes; elsewhere its error E = I - B A =
// (1 - i) B Ai can exceed 1 (where Ai dominates, or Ar + Ai is nearly singular), and the
// solve then amplifies components which the coarse correction is meant to remove. The
// measure is ‖Eᵏ w‖ / ‖w‖ for k = 1, 2, 3 successive real-valued solves of a vector w which
// is the same function on the mesh for every partition when the space is first order
// (MeshNoise), so that the choice, made from global norms, does not depend on the number of
// processes. The measurement stops at the first value above 1 or not a number: then the
// real-valued solver is released and the complex-valued matrix is factored instead.
// Otherwise the cost is three solves per operator.
//
// The choice only goes one way and belongs to the solver: while it uses the real-valued
// approximation, it measures every operator it is given, and once it has chosen the
// complex-valued matrix, it factors the complex-valued matrix of all its later operators
// (the later frequencies of a sweep) without a measurement. A new solver (on the next mesh
// of an adaptive refinement) measures again.
//
class MeasuredComplexCoarseSolver : public Solver<ComplexOperator>
{
private:
  MPI_Comm comm;

  // Solver for the real-valued approximation, released when the complex-valued matrix is
  // chosen, and solver for the complex-valued matrix, which has no operator until then.
  std::unique_ptr<Solver<ComplexOperator>> real_pc, complex_pc;

  // Measurement vector, and errors of the most recent measurement.
  Vector w;
  std::vector<double> real_solve_error;

  // Whether a measurement of this solver has been printed.
  bool printed = false;

  std::vector<double> MeasureRealSolveError(const ComplexOperator &op) const;

public:
  MeasuredComplexCoarseSolver(MPI_Comm comm,
                              std::unique_ptr<Solver<ComplexOperator>> &&real_pc,
                              std::unique_ptr<Solver<ComplexOperator>> &&complex_pc,
                              Vector &&w);

  void SetInitialGuess(bool guess) override;

  void SetOperator(const ComplexOperator &op) override;

  void Mult(const ComplexVector &x, ComplexVector &y) const override;

  // Whether the solver uses the complex-valued matrix.
  bool UsesComplexMatrix() const { return !real_pc; }

  // The errors ‖Eᵏ w‖ / ‖w‖ of the most recent measurement (empty if none).
  const std::vector<double> &GetRealSolveError() const { return real_solve_error; }
};

}  // namespace palace

#endif  // PALACE_LINALG_COMPLEX_COARSE_HPP
