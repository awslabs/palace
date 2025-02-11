// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_LINALG_FLOQUET_CORR_HPP
#define PALACE_LINALG_FLOQUET_CORR_HPP

#include <memory>
#include <vector>
#include "linalg/ksp.hpp"
#include "linalg/operator.hpp"
#include "linalg/vector.hpp"

namespace mfem
{

template <typename T>
class Array;

}  // namespace mfem

namespace palace
{

class FiniteElementSpace;
class MaterialOperator;

//
// This solver calculates a correction for the magnetic flux density field
// when Floquet periodicity is imposed. The correction is the cross product
// of the Floquet wave vector with the electric field.
//
template <typename VecType>
class FloquetCorrSolver
{
  using OperType = typename std::conditional<std::is_same<VecType, ComplexVector>::value,
                                             ComplexOperator, Operator>::type;
  using ScalarType =
      typename std::conditional<std::is_same<OperType, ComplexOperator>::value,
                                std::complex<double>, double>::type;

private:
  // Operators for the floquet correction.
  std::unique_ptr<OperType> M, Cross;

  // Linear solver for the linear system M y = x.
  std::unique_ptr<BaseKspSolver<OperType>> ksp;

  // Workspace objects for solver application.
  mutable VecType rhs;

public:
  FloquetCorrSolver(const MaterialOperator &mat_op, FiniteElementSpace &nd_fespace,
                    FiniteElementSpace &rt_fespace, double tol, int max_it, int print);

  // Given a vector of Nedelec dofs for an arbitrary vector field, compute
  // the Raviart-Thomas space field y = [kp x] x, where [kp x] is a matrix
  // representing the action of the cross product with the Floquet wave vector.
  void Mult(const VecType &x, VecType &y) const;
  void AddMult(const VecType &x, VecType &y, ScalarType a = 1.0) const;
};

}  // namespace palace

#endif  // PALACE_LINALG_FLOQUET_CORR_HPP
