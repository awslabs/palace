// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_LINALG_JACOBI_SMOOTHER_HPP
#define PALACE_LINALG_JACOBI_SMOOTHER_HPP

#include "linalg/operator.hpp"
#include "linalg/solver.hpp"
#include "linalg/solverstats.hpp"
#include "linalg/vector.hpp"

namespace palace
{

//
// Simple Jacobi smoother using the diagonal vector from OperType::AssembleDiagonal(),
// which allows for (approximate) diagonal construction for matrix-free operators.
//
template <OperatorType OperType>
class JacobiSmoother : public Solver<OperType>
{
  using VecType = typename Solver<OperType>::VecType;

private:
  // MPI communicator associated with the solver operator and vectors.
  MPI_Comm comm;

  // Inverse diagonal scaling of the operator (real-valued for now).
  VecType dinv;

  // Damping factor and scaling factor for maximum eigenvalue.
  double omega, sf_max;

  // Performance tally id for the eigenvalue estimate (keyed by the role at construction).
  int perf_lambda;

public:
  JacobiSmoother(MPI_Comm comm, double omega = 1.0, double sf_max = 1.0)
    : Solver<OperType>(), comm(comm), omega(omega), sf_max(sf_max),
      perf_lambda((omega == 0.0) ? perf::LambdaMaxTally() : -1)
  {
  }

  void SetOperator(const OperType &op) override;

  void Mult(const VecType &x, VecType &y) const override;

  void MultTranspose(const VecType &x, VecType &y) const override { Mult(x, y); }
};

}  // namespace palace

#endif  // PALACE_LINALG_JACOBI_SMOOTHER_HPP
