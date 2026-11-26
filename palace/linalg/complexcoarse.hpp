// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_LINALG_COMPLEX_COARSE_HPP
#define PALACE_LINALG_COMPLEX_COARSE_HPP

#include <functional>
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
// Sparse direct coarse solver for a complex-valued operator Ar + i Ai, which factors one
// of three matrices:
// - Ar in real arithmetic, when the imaginary part is zero (see HasZeroImag), as for a
//   model without losses: this is exact, needs half of the memory of the complex-valued
//   matrix, and the same factors solve the real and imaginary parts of the right-hand side;
// - Ar + i Ai in complex arithmetic (Solver.Linear.ComplexCoarseSolve);
// - the real-valued approximation Ar + Ai in real arithmetic, when ComplexCoarseSolve is
//   omitted (the solver is then given a measurement vector) and the measured error of this
//   approximation stays below 1.
//
// The real-valued solve B = (Ar + Ai)⁻¹ is exact where Ai vanishes; elsewhere its error
// E = I - B A = (1 - i) B Ai can exceed 1 (where Ai dominates, or Ar + Ai is nearly
// singular), and the solve then amplifies components which the coarse correction is meant
// to remove. The measure is ‖Eᵏ w‖ / ‖w‖ for k = 1, 2, 3 successive real-valued solves of
// a vector w which is the same function on the mesh for every partition when the space is
// first order (MeshNoise), so that the choice, made from global norms, does not depend on
// the number of processes. The measurement stops at the first value above 1 or not a
// number (a failed factorization): then the real-valued solver is released and the
// complex-valued matrix is factored instead. Otherwise the cost is three solves per
// operator.
//
// The choice between the real part and the other two is made for each operator, so it
// follows the operators of a frequency sweep. The measured choice only goes from the
// real-valued approximation to the complex-valued matrix: while the solver uses the
// approximation, it measures every operator with an imaginary part which it is given, and
// once a measurement has chosen the complex-valued matrix, the later operators with an
// imaginary part (the later frequencies of a sweep) are factored in complex arithmetic
// without a measurement. A new solver (on the next mesh of an adaptive refinement) measures
// again. Only the solver in use exists: it is built when it is first needed, and the other
// one is destroyed, with its factors, before the matrix is factored.
//
// Where the coarse-level operator is assembled differently for the real-valued
// approximation and for the complex-valued solve (the second-order absorbing boundary term
// and the PML terms, see SpaceOperator::AssemblePreconditioner), the measured choice is
// given the operator for the real-valued approximation and its difference from the
// complex-valued operator, and tells the assembler when it has chosen the latter (see
// SetChoiceContext), so that it measures and factors the same matrices as with
// ComplexCoarseSolve false and true.
//
class ComplexCoarseSolver : public Solver<ComplexOperator>
{
public:
  using SolverFactory = std::function<std::unique_ptr<Solver<ComplexOperator>>()>;

private:
  MPI_Comm comm;

  // Solver for a real-valued matrix (an MfemWrapperSolver, which factors the real part of
  // an operator without imaginary part and, for the measured choice, the sum of the real
  // and imaginary parts otherwise) and solver in complex arithmetic.
  SolverFactory make_real_pc, make_complex_pc;
  std::unique_ptr<Solver<ComplexOperator>> real_pc, complex_pc;

  // Measurement vector of the measured choice (empty without it), errors of its most recent
  // measurement, and whether a measurement has chosen the complex-valued matrix.
  Vector w;
  std::vector<double> real_solve_error;
  bool measured_complex = false;

  // Context of the next call of SetOperator (see SetChoiceContext, not owned).
  bool *complex_chosen = nullptr;
  const ComplexOperator *correction_op = nullptr;

  // Whether to print which matrix is factored, when it changes, and the matrix of the
  // previous operator.
  bool print;
  enum class Matrix : char
  {
    NONE,
    REAL_PART,
    REAL_APPROXIMATION,
    COMPLEX
  };
  Matrix matrix = Matrix::NONE;

  // The errors ‖Eᵏ w‖ / ‖w‖, k = 1, 2, 3, of the real-valued solver for the operator (plus
  // the correction, if given), up to the first one which is above 1 or not a number.
  std::vector<double> MeasureRealSolveError(const ComplexOperator &op,
                                            const ComplexOperator *correction) const;

public:
  // With an empty measurement vector, an operator with an imaginary part is always factored
  // in complex arithmetic.
  ComplexCoarseSolver(MPI_Comm comm, SolverFactory make_real_pc,
                      SolverFactory make_complex_pc, bool print, Vector &&w = Vector());

  // The solver in use: at most one of the two is not null.
  Solver<ComplexOperator> *GetRealSolver() { return real_pc.get(); }
  Solver<ComplexOperator> *GetComplexSolver() { return complex_pc.get(); }

  // The errors of the most recent measurement (empty if none).
  const std::vector<double> &GetRealSolveError() const { return real_solve_error; }

  void SetInitialGuess(bool guess) override;

  // For the measured choice with operators whose coarse level is assembled differently for
  // the real-valued approximation and for the complex-valued solve (see
  // SpaceOperator::AssemblePreconditioner), set before each call of SetOperator (the
  // pointers are only used during that call):
  // - a flag which is shared with the assembler of the operators. The solver sets it when
  //   it chooses the complex-valued matrix, and the later operators are then assembled
  //   for the complex-valued solve. If it is already set (by another solver with operators
  //   of the same assembler), the solver uses the complex-valued matrix without a
  //   measurement;
  // - a correction (null if none): the operator given to SetOperator is assembled for the
  //   real-valued approximation, whose solver factors the sum of its real and imaginary
  //   parts, and the complex-valued operator is that operator plus the correction. The
  //   error of the real-valued solve is measured for the complex-valued operator, which is
  //   the matrix factored in complex arithmetic.
  void SetChoiceContext(bool *chosen, const ComplexOperator *correction)
  {
    complex_chosen = chosen;
    correction_op = correction;
  }

  void SetOperator(const ComplexOperator &op) override;

  void Mult(const ComplexVector &x, ComplexVector &y) const override;
};

}  // namespace palace

#endif  // PALACE_LINALG_COMPLEX_COARSE_HPP
