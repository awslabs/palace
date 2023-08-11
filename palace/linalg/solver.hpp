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
class Solver
{
protected:

  static constexpr bool is_complex = std::is_same<OperType, ComplexOperator>::value;
  static constexpr bool is_real = std::is_same<OperType, Operator>::value;

  static_assert(is_complex || is_real,
                "Solver can only be defined for OperType = Operator or ComplexOperator!");

  typedef typename std::conditional<std::is_same<OperType, ComplexOperator>::value,
                                    ComplexVector, Vector>::type VecType;

  // Whether or not to use the second argument of Mult() as an initial guess.
  bool initial_guess;

public:
  Solver(bool initial_guess = false) : initial_guess(initial_guess) {}
  virtual ~Solver() = default;

  // Configure whether or not to use an initial guess when applying the solver.
  virtual void SetInitialGuess(bool guess) { initial_guess = guess; }

  // Set the operator associated with the solver, or update it if called repeatedly.
  virtual void SetOperator(const OperType &op) = 0;

  // Apply the solver.
  virtual void Mult(const VecType &x, VecType &y) const = 0;

  // Apply the solver for the transpose problem.
  virtual void MultTranspose(const VecType &x, VecType &y) const
  {
    MFEM_ABORT("MultTranspose() is not implemented for base class Solver<OperType>!");
  }
};

// This solver wraps a real-valued mfem::Solver for application to complex-valued problems
// as a preconditioner inside of a Solver<OperType>
template <typename OperType>
class WrapperSolver : public Solver<OperType>
{
protected:
  using Solver<OperType>::is_real;
  using Solver<OperType>::is_complex;

  typedef typename Solver<OperType>::VecType VecType;

  std::unique_ptr<mfem::Solver> pc;

public:
  WrapperSolver(std::unique_ptr<mfem::Solver> &&pc)
    : Solver<OperType>(pc->iterative_mode), pc(std::move(pc))
  {
  }

  void SetInitialGuess(bool guess) override
  {
    Solver<OperType>::SetInitialGuess(guess);
    pc->iterative_mode = guess;
  }

  void SetOperator(const OperType &op) override
  {
    if constexpr (is_real)
    {
      pc->SetOperator(op);
    }
    else if constexpr (is_complex)
    {
      MFEM_VERIFY(op.IsReal() && op.HasReal(),
                  "WrapperSolver::SetOperator requires an operator which is purely real for "
                  "mfem::Solver!");
      pc->SetOperator(*op.Real());
    }
  }

  void Mult(const VecType &x, VecType &y) const override
  {
    if constexpr (is_real)
    {
      pc->Mult(x, y);
    }
    else if constexpr (is_complex)
    {
      mfem::Array<const Vector *> X(2);
      mfem::Array<Vector *> Y(2);
      X[0] = &x.Real();
      X[1] = &x.Imag();
      Y[0] = &y.Real();
      Y[1] = &y.Imag();
      pc->ArrayMult(X, Y);
    }
  }
};

}  // namespace palace

#endif  // PALACE_LINALG_SOLVER_HPP
