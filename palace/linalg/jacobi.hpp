// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_LINALG_JACOBI_SMOOTHER_HPP
#define PALACE_LINALG_JACOBI_SMOOTHER_HPP

#include "linalg/operator.hpp"
#include "linalg/relaxation.hpp"
#include "linalg/vector.hpp"

namespace palace
{

//
// Simple Jacobi smoother using the diagonal vector from OperType::AssembleDiagonal(),
// which allows for (approximate) diagonal construction for matrix-free operators.
//
template <OperatorType OperType>
class JacobiSmoother : public Relaxation<OperType>
{
protected:
  using VecType = typename Solver<OperType>::VecType;

  // MPI communicator associated with the solver operator and vectors.
  MPI_Comm comm;

  // Inverse diagonal scaling of the operator (real-valued for now).
  VecType dinv;

private:
  // Damping factor and scaling factor for maximum eigenvalue.
  double omega, sf_max;

  // Whether, on all processes, the scaling is real-valued and callers may apply it in their
  // own host loops.
  bool host_real;

public:
  JacobiSmoother(MPI_Comm comm, double omega = 1.0, double sf_max = 1.0)
    : Relaxation<OperType>(), comm(comm), omega(omega), sf_max(sf_max), host_real(false)
  {
  }

  void SetOperator(const OperType &op) override;

  void Mult(const VecType &x, VecType &y) const override;

  void MultTranspose(const VecType &x, VecType &y) const override { Mult(x, y); }

  void ScaledMult(double s, const VecType &r, VecType &d) const override;

  void ScaledAddMult(double sd, double sr, const VecType &r, VecType &d) const override;

  double EstimateMaxEigenvalue(const OperType &A, VecType &x0) const override;

  // Return the scaling ω D⁻¹ applied by Mult, if it is real-valued and callers may apply it
  // in their own single-threaded host loops (the same on all processes), or nullptr.
  const Vector *GetHostRealScaling() const;
};

}  // namespace palace

#endif  // PALACE_LINALG_JACOBI_SMOOTHER_HPP
