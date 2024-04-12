// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_LINALG_KSP_HPP
#define PALACE_LINALG_KSP_HPP

#include <memory>
#include <type_traits>
#include "linalg/iterative.hpp"
#include "linalg/operator.hpp"
#include "linalg/solver.hpp"

namespace palace
{

class FiniteElementSpaceHierarchy;
class IoData;

//
// Linear solver class composing an iterative solver and preconditioner object.
//
template <typename OperType>
class BaseKspSolver
{
  static_assert(std::is_same<OperType, Operator>::value ||
                    std::is_same<OperType, ComplexOperator>::value,
                "Solver can only be defined for OperType = Operator or ComplexOperator!");

  using VecType = typename std::conditional<std::is_same<OperType, ComplexOperator>::value,
                                            ComplexVector, Vector>::type;

protected:
  // The actual solver and preconditioner objects.
  std::unique_ptr<IterativeSolver<OperType>> ksp;
  std::unique_ptr<Solver<OperType>> pc;

  // Counters for number of calls to Mult method for linear solves, and cumulative number
  // of iterations.
  mutable int ksp_mult, ksp_mult_it;

  // Enable timer contribution for Timer::KSP_PRECONDITIONER.
  bool use_timer;

public:
  BaseKspSolver(const IoData &iodata, FiniteElementSpaceHierarchy &fespaces,
                FiniteElementSpaceHierarchy *aux_fespaces = nullptr);
  BaseKspSolver(std::unique_ptr<IterativeSolver<OperType>> &&ksp,
                std::unique_ptr<Solver<OperType>> &&pc);

  int NumTotalMult() const { return ksp_mult; }
  int NumTotalMultIterations() const { return ksp_mult_it; }

  void SetOperators(const OperType &op, const OperType &pc_op);

  void Mult(const VecType &x, VecType &y) const;
};

using KspSolver = BaseKspSolver<Operator>;
using ComplexKspSolver = BaseKspSolver<ComplexOperator>;

}  // namespace palace

#endif  // PALACE_LINALG_KSP_HPP
