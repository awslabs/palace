// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_LINALG_CHEBYSHEV_SMOOTHER_HPP
#define PALACE_LINALG_CHEBYSHEV_SMOOTHER_HPP

#include <memory>
#include "linalg/operator.hpp"
#include "linalg/relaxation.hpp"
#include "linalg/solver.hpp"
#include "linalg/vector.hpp"

namespace palace
{

//
// Matrix-free Chebyshev smoothing, scaled by a relaxation M⁻¹ (diagonally-scaled for the
// point Jacobi method, which is the default). This is largely the same as
// mfem::OperatorChebyshevSmoother allows a nonzero initial guess and uses alternative
// methods to estimate the largest eigenvalue. We use a smoother based on Chebyshev
// polynomials of the 4th-kind as proposed in recent work.
// Reference: Phillips and Fischer, Optimal Chebyshev smoothers and one-sided V-cycles,
//            arXiv:2210.03179v1 (2022).
//
template <OperatorType OperType>
class ChebyshevSmoother : public Solver<OperType>
{
  using VecType = typename Solver<OperType>::VecType;

private:
  // Number of smoother iterations and polynomial order.
  const int pc_it, order;

  // System matrix (not owned).
  const OperType *A;

  // Relaxation M⁻¹ which scales the operator in the polynomial.
  std::unique_ptr<Relaxation<OperType>> M;

  // Maximum eigenvalue of M⁻¹A for Chebyshev polynomial smoothing.
  double lambda_max, sf_max;

  // Dominant vector of the latest maximum eigenvalue estimate, from which the estimate for
  // the next operator starts.
  VecType lambda_vec;

  // Temporary vector for smoother application.
  mutable VecType d, r;

  // Apply the polynomial to the residual r of y: y = y + p(M⁻¹A) M⁻¹ r.
  void Polynomial(VecType &y, VecType &r) const;

public:
  // Diagonally-scaled smoother (point Jacobi relaxation).
  ChebyshevSmoother(MPI_Comm comm, int smooth_it, int poly_order, double sf_max);

  // Smoother scaled by the given relaxation.
  ChebyshevSmoother(std::unique_ptr<Relaxation<OperType>> &&M, int smooth_it,
                    int poly_order, double sf_max);

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

  // Assumes operator symmetry.
  void MultTranspose2(const VecType &x, VecType &y, VecType &r) const override;
};

//
// Matrix-free Chebyshev smoothing scaled by a relaxation (see above) using standard
// 1st-kind Chebyshev polynomials. Reference: Adams et al., Parallel multigrid smoothing:
// polynomial versus Gauss–Seidel,
//            JCP (2003).
//
template <OperatorType OperType>
class ChebyshevSmoother1stKind : public Solver<OperType>
{
  using VecType = typename Solver<OperType>::VecType;

private:
  // Number of smoother iterations and polynomial order.
  const int pc_it, order;

  // System matrix (not owned).
  const OperType *A;

  // Relaxation M⁻¹ which scales the operator in the polynomial.
  std::unique_ptr<Relaxation<OperType>> M;

  // Parameters depending on maximum and minimum operator eigenvalue estimates for Chebyshev
  // polynomial smoothing.
  double theta, delta, sf_max, sf_min;

  // Dominant vector of the latest maximum eigenvalue estimate, from which the estimate for
  // the next operator starts.
  VecType lambda_vec;

  // Temporary vector for smoother application.
  mutable VecType d, r;

  // Apply the polynomial to the residual r of y: y = y + p(M⁻¹A) M⁻¹ r.
  void Polynomial(VecType &y, VecType &r) const;

public:
  // Diagonally-scaled smoother (point Jacobi relaxation).
  ChebyshevSmoother1stKind(MPI_Comm comm, int smooth_it, int poly_order, double sf_max,
                           double sf_min);

  // Smoother scaled by the given relaxation.
  ChebyshevSmoother1stKind(std::unique_ptr<Relaxation<OperType>> &&M, int smooth_it,
                           int poly_order, double sf_max, double sf_min);

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

  // Assumes operator symmetry.
  void MultTranspose2(const VecType &x, VecType &y, VecType &r) const override;
};

}  // namespace palace

#endif  // PALACE_LINALG_CHEBYSHEV_SMOOTHER_HPP
