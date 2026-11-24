// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_LINALG_RELAXATION_HPP
#define PALACE_LINALG_RELAXATION_HPP

#include "linalg/operator.hpp"
#include "linalg/solver.hpp"
#include "linalg/vector.hpp"

namespace palace
{

//
// Abstract base class for a stationary relaxation M⁻¹ ≈ A⁻¹ constructed from an operator A:
// the point Jacobi method, where M is the diagonal of A, or the block Jacobi method. It is
// what a polynomial smoother accelerates, for which it provides the scaled applications
// below and an estimate of the largest eigenvalue of M⁻¹A (Mult applies M⁻¹).
//
template <OperatorType OperType>
class Relaxation : public Solver<OperType>
{
protected:
  using VecType = typename Solver<OperType>::VecType;

public:
  // d = s M⁻¹ r.
  virtual void ScaledMult(double s, const VecType &r, VecType &d) const = 0;

  // d = sd d + sr M⁻¹ r.
  virtual void ScaledAddMult(double sd, double sr, const VecType &r, VecType &d) const = 0;

  // Estimate the largest eigenvalue of M⁻¹A, for the operator A given to SetOperator. A
  // vector x0 of the size of the operator is the starting vector of the estimate, and on
  // return x0 is the dominant vector it found (see linalg::SpectralNorm).
  virtual double EstimateMaxEigenvalue(const OperType &A, VecType &x0) const = 0;

  // A relaxation can leave a part of M⁻¹ out of the applications above, and apply it as a
  // correction y += S r of a smoothing sweep instead, where r = x - A y is the residual
  // after the polynomial of the sweep (before it, for the transposed sweep). The block
  // Jacobi method does this for the blocks it solves approximately, to solve them once per
  // sweep and not once per term of the polynomial.
  virtual bool HasSweepCorrection() const { return false; }
  virtual void AddSweepCorrection(const VecType &r, VecType &y) const {}
};

}  // namespace palace

#endif  // PALACE_LINALG_RELAXATION_HPP
