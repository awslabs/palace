// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_LINALG_SOLVER_HPP
#define PALACE_LINALG_SOLVER_HPP

#include <type_traits>
#include <mfem.hpp>
#include "linalg/operator.hpp"
#include "linalg/vector.hpp"

namespace palace
{

//
// The base Solver<OperType> class is a templated version of mfem::Solver for operation with
// real- or complex-valued operators.
//

// Abstract base class for real-valued or complex-valued solvers.
template <typename OperType>
class Solver : public OperType
{
  static_assert(std::is_same<OperType, Operator>::value ||
                    std::is_same<OperType, ComplexOperator>::value,
                "Solver can only be defined for OperType = Operator or ComplexOperator!");

protected:
  using VecType = typename std::conditional<std::is_same<OperType, ComplexOperator>::value,
                                            ComplexVector, Vector>::type;

  // Whether or not to use the second argument of Mult() as an initial guess.
  bool initial_guess;

public:
  Solver(bool initial_guess = false) : OperType(), initial_guess(initial_guess) {}
  virtual ~Solver() = default;

  // Configure whether or not to use an initial guess when applying the solver.
  virtual void SetInitialGuess(bool guess) { initial_guess = guess; }

  // Set the operator associated with the solver, or update it if called repeatedly.
  virtual void SetOperator(const OperType &op) = 0;

  // Apply the solver for the transpose problem.
  void MultTranspose(const VecType &x, VecType &y) const override
  {
    MFEM_ABORT("MultTranspose() is not implemented for base class Solver<OperType>!");
  }

  // Apply the solver with a preallocated temporary storage vector.
  virtual void Mult2(const VecType &x, VecType &y, VecType &r) const
  {
    MFEM_ABORT("Mult2() with temporary storage vector is not implemented for base class "
               "Solver<OperType>!");
  }

  // Apply the solver for the transpose problem with a preallocated temporary storage
  //  vector.
  virtual void MultTranspose2(const VecType &x, VecType &y, VecType &r) const
  {
    MFEM_ABORT("MultTranspose2() with temporary storage vector is not implemented for base "
               "class Solver<OperType>!");
  }
};

// This solver wraps a real-valued mfem::Solver for application to complex-valued problems
// as a preconditioner inside of a Solver<OperType> or for assembling the matrix-free
// preconditioner operator as an mfem::HypreParMatrix.
template <typename OperType>
class MfemWrapperSolver : public Solver<OperType>
{
  using VecType = typename Solver<OperType>::VecType;

private:
  // The actual mfem::Solver.
  std::unique_ptr<mfem::Solver> pc;

  // System matrix A in parallel assembled form.
  std::unique_ptr<mfem::HypreParMatrix> A;

  // Whether or not to save the parallel assembled matrix after calling
  // mfem::Solver::SetOperator (some solvers copy their input).
  bool save_assembled;

  // Whether to use the exact complex-valued system matrix or the real-valued
  // approximation A = Ar + Ai.
  bool complex_matrix = true;

public:
  MfemWrapperSolver(std::unique_ptr<mfem::Solver> &&pc, bool save_assembled = true,
                    bool complex_matrix = true)
    : Solver<OperType>(pc->iterative_mode), pc(std::move(pc)),
      save_assembled(save_assembled), complex_matrix(complex_matrix)
  {
  }

  // Access the underlying solver.
  const mfem::Solver &GetSolver() { return *pc; }

  // Configure whether or not to save the assembled operator.
  void SetSaveAssembled(bool save) { save_assembled = save; }

  void SetInitialGuess(bool guess) override
  {
    Solver<OperType>::SetInitialGuess(guess);
    pc->iterative_mode = guess;
  }

  void SetOperator(const OperType &op) override;

  void Mult(const VecType &x, VecType &y) const override;
};

}  // namespace palace

#endif  // PALACE_LINALG_SOLVER_HPP
