// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_LINALG_CHEBYSHEV_SMOOTHER_HPP
#define PALACE_LINALG_CHEBYSHEV_SMOOTHER_HPP

#include <mfem.hpp>
#include "linalg/operator.hpp"
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
class ChebyshevSmoother : public mfem::Solver
{
private:
  // Number of smoother iterations and polynomial order.
  const int pc_it, order;

  // System matrix (not owned).
  const Operator *A;

  // Inverse diagonal scaling of the operator.
  Vector dinv;

  // Maximum operator eigenvalue for Chebyshev polynomial smoothing.
  double lambda_max;

  // Temporary vectors for smoother application.
  mutable Vector r, d;

public:
  ChebyshevSmoother(int smooth_it, int poly_order);

  void SetOperator(const Operator &op) override;

  void Mult(const Vector &x, Vector &y) const override;

  void MultTranspose(const Vector &x, Vector &y) const override
  {
    Mult(x, y);  // Assumes operator symmetry
  }

  // XX TODO REMOVE...
  //  void Mult(const mfem::Vector &x, mfem::Vector &y) const override
  //  {
  //    mfem::Array<const mfem::Vector *> X(1);
  //    mfem::Array<mfem::Vector *> Y(1);
  //    X[0] = &x;
  //    Y[0] = &y;
  //    ArrayMult(X, Y);
  //  }

  // void ArrayMult(const mfem::Array<const mfem::Vector *> &X,
  //                mfem::Array<mfem::Vector *> &Y) const override;

  // void ArrayMultTranspose(const mfem::Array<const mfem::Vector *> &X,
  //                         mfem::Array<mfem::Vector *> &Y) const override
  // {
  //   ArrayMult(X, Y);  // Assumes operator symmetry
  // }
};

}  // namespace palace

#endif  // PALACE_LINALG_CHEBYSHEV_SMOOTHER_HPP
