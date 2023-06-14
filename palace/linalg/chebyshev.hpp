// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_LINALG_CHEBYSHEV_SMOOTHER_HPP
#define PALACE_LINALG_CHEBYSHEV_SMOOTHER_HPP

#include "linalg/operator.hpp"
#include "linalg/solver.hpp"
#include "linalg/vector.hpp"

namespace palace
{

//
// Matrix-free diagonally-scaled Chebyshev smoothing. This is largely the same as
// mfem::OperatorChebyshevSmoother allows a nonzero initial guess and uses alternative
// methods to estimate the largest eigenvalue. See also Phillips and Fischer, Optimal
// Chebyshev smoothers and one-sided V-cycles, arXiv:2210.03179v1 (2022) for reference on
// the 4th-kind Chebyshev polynomial smoother.
//
template <typename OperType>
class ChebyshevSmoother : public Solver<OperType>
{
protected:
  typedef typename Solver<OperType>::VecType VecType;

private:
  // Number of smoother iterations and polynomial order.
  const int pc_it, order;

  // System matrix (not owned).
  const OperType *A;

  // Inverse diagonal scaling of the operator (real-valued for now).
  Vector dinv;

  // Maximum operator eigenvalue for Chebyshev polynomial smoothing.
  double lambda_max;

  // Temporary vectors for smoother application.
  mutable VecType r, d;

public:
  ChebyshevSmoother(int smooth_it, int poly_order);

  void SetOperator(const OperType &op) override;

  void Mult(const VecType &x, VecType &y) const override;

  void MultTranspose(const VecType &x, VecType &y) const override
  {
    Mult(x, y);  // Assumes operator symmetry
  }
};

}  // namespace palace

#endif  // PALACE_LINALG_CHEBYSHEV_SMOOTHER_HPP
