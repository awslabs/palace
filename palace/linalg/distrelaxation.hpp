// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_LINALG_DIST_RELAXATION_SMOOTHER_HPP
#define PALACE_LINALG_DIST_RELAXATION_SMOOTHER_HPP

#include <memory>
#include "linalg/operator.hpp"
#include "linalg/solver.hpp"
#include "linalg/vector.hpp"

namespace mfem
{

template <typename T>
class Array;
class ParFiniteElementSpace;

}  // namespace mfem

namespace palace
{

//
// Hiptmair distributive relaxation smoother applying smoothers to both the operator in the
// primary space as well as its projection into an auxiliary space.
// Reference: Hiptmair, Multigrid method for Maxwell's equations, SIAM J. Numer. Anal.
//            (1998).
//
template <typename OperType>
class DistRelaxationSmoother : public Solver<OperType>
{
private:
  // Number of smoother iterations.
  const int pc_it;

  // System matrix and its projection GᵀAG (not owned).
  const OperType *A, *A_G;
  const mfem::Array<int> *dbc_tdof_list_G;

  // Discrete gradient matrix.
  std::unique_ptr<Operator> G;

  // Point smoother objects for each matrix.
  mutable std::unique_ptr<Solver<OperType>> B;
  std::unique_ptr<Solver<OperType>> B_G;

  // Temporary vectors for smoother application.
  mutable VecType r, x_G, y_G;

public:
  DistRelaxationSmoother(mfem::ParFiniteElementSpace &nd_fespace,
                         mfem::ParFiniteElementSpace &h1_fespace, int smooth_it,
                         int cheby_smooth_it, int cheby_order);

  void SetOperator(const OperType &op) override
  {
    MFEM_ABORT("SetOperator with a single operator is not implemented for "
               "DistRelaxationSmoother, use the two argument signature instead!");
  }
  void SetOperators(const OperType &op, const OperType &op_G);

  void Mult(const VecType &x, VecType &y) const override;

  void MultTranspose(const VecType &x, VecType &y) const override;
};

}  // namespace palace

#endif  // PALACE_LINALG_DIST_RELAXATION_SMOOTHER_HPP
