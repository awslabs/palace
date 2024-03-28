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
// methods to estimate the largest eigenvalue. We use a smoother based on Chebyshev
// polynomials of the 4th-kind as proposed in recent work.
// Reference: Phillips and Fischer, Optimal Chebyshev smoothers and one-sided V-cycles,
//            arXiv:2210.03179v1 (2022).
//
template <typename OperType>
class ChebyshevSmoother : public Solver<OperType>
{
  using VecType = typename Solver<OperType>::VecType;

private:
  // MPI communicator associated with the solver operator and vectors.
  MPI_Comm comm;

  // Number of smoother iterations and polynomial order.
  const int pc_it, order;

  // System matrix (not owned).
  const OperType *A;

  // Inverse diagonal scaling of the operator (real-valued for now).
  VecType dinv;

  // Maximum operator eigenvalue for Chebyshev polynomial smoothing.
  double lambda_max, sf_max;

  // Temporary vector for smoother application.
  mutable VecType d, r;

public:
  ChebyshevSmoother(MPI_Comm comm, int smooth_it, int poly_order, double sf_max);

  void SetOperator(const OperType &op) override;

  void Mult(const VecType &x, VecType &y) const override
  {
    if (r.Size() != y.Size())
    {
      r.SetSize(y.Size());
      r.UseDevice(true);
    }
    Mult2(x, y, r);
  }

  void MultTranspose(const VecType &x, VecType &y) const override
  {
    if (r.Size() != y.Size())
    {
      r.SetSize(y.Size());
      r.UseDevice(true);
    }
    MultTranspose2(x, y, r);
  }

  void Mult2(const VecType &x, VecType &y, VecType &r) const override;

  void MultTranspose2(const VecType &x, VecType &y, VecType &r) const override
  {
    Mult2(x, y, r);  // Assumes operator symmetry
  }
};

//
// Matrix-free diagonally-scaled Chebyshev smoothing using standard 1st-kind Chebyshev
// polynomials.
// Reference: Adams et al., Parallel multigrid smoothing: polynomial versus Gauss–Seidel,
//            JCP (2003).
//
template <typename OperType>
class ChebyshevSmoother1stKind : public Solver<OperType>
{
  using VecType = typename Solver<OperType>::VecType;

private:
  // MPI communicator associated with the solver operator and vectors.
  MPI_Comm comm;

  // Number of smoother iterations and polynomial order.
  const int pc_it, order;

  // System matrix (not owned).
  const OperType *A;

  // Inverse diagonal scaling of the operator (real-valued for now).
  VecType dinv;

  // Parameters depending on maximum and minimum operator eigenvalue estimates for Chebyshev
  // polynomial smoothing.
  double theta, delta, sf_max, sf_min;

  // Temporary vector for smoother application.
  mutable VecType d, r;

public:
  ChebyshevSmoother1stKind(MPI_Comm comm, int smooth_it, int poly_order, double sf_max,
                           double sf_min);

  void SetOperator(const OperType &op) override;

  void Mult(const VecType &x, VecType &y) const override
  {
    if (r.Size() != y.Size())
    {
      r.SetSize(y.Size());
      r.UseDevice(true);
    }
    Mult2(x, y, r);
  }

  void MultTranspose(const VecType &x, VecType &y) const override
  {
    if (r.Size() != y.Size())
    {
      r.SetSize(y.Size());
      r.UseDevice(true);
    }
    MultTranspose2(x, y, r);
  }

  void Mult2(const VecType &x, VecType &y, VecType &r) const override;

  void MultTranspose2(const VecType &x, VecType &y, VecType &r) const override
  {
    Mult2(x, y, r);  // Assumes operator symmetry
  }
};

}  // namespace palace

#endif  // PALACE_LINALG_CHEBYSHEV_SMOOTHER_HPP
