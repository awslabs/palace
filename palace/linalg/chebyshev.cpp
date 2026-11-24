// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "chebyshev.hpp"

#include <cmath>
#include "linalg/jacobi.hpp"

namespace palace
{

namespace
{

template <bool Transpose = false>
inline void ApplyOp(const Operator &A, const Vector &x, Vector &y)
{
  A.Mult(x, y);
}

template <bool Transpose = false>
inline void ApplyOp(const ComplexOperator &A, const ComplexVector &x, ComplexVector &y)
{
  if constexpr (!Transpose)
  {
    A.Mult(x, y);
  }
  else
  {
    A.MultHermitianTranspose(x, y);
  }
}

template <bool Transpose = false>
inline void ApplyOp(const Operator &A, const Vector &x, Vector &y, const double a)
{
  A.AddMult(x, y, a);
}

template <bool Transpose = false>
inline void ApplyOp(const ComplexOperator &A, const ComplexVector &x, ComplexVector &y,
                    const double a)
{
  if constexpr (!Transpose)
  {
    A.AddMult(x, y, a);
  }
  else
  {
    A.AddMultHermitianTranspose(x, y, a);
  }
}

}  // namespace

template <OperatorType OperType>
ChebyshevSmoother<OperType>::ChebyshevSmoother(MPI_Comm comm, int smooth_it, int poly_order,
                                               double sf_max)
  : ChebyshevSmoother(std::make_unique<JacobiSmoother<OperType>>(comm), smooth_it,
                      poly_order, sf_max)
{
}

template <OperatorType OperType>
ChebyshevSmoother<OperType>::ChebyshevSmoother(std::unique_ptr<Relaxation<OperType>> &&M,
                                               int smooth_it, int poly_order, double sf_max)
  : Solver<OperType>(), pc_it(smooth_it), order(poly_order), A(nullptr), M(std::move(M)),
    lambda_max(0.0), sf_max(sf_max)
{
  MFEM_VERIFY(this->M, "Chebyshev smoothing requires a relaxation method!");
  MFEM_VERIFY(order > 0, "Polynomial order for Chebyshev smoothing must be positive!");
}

template <OperatorType OperType>
void ChebyshevSmoother<OperType>::SetOperator(const OperType &op)
{
  A = &op;
  d.SetSize(op.Height());
  d.UseDevice(true);
  M->SetOperator(op);

  // Set up Chebyshev coefficients using the computed maximum eigenvalue estimate. See
  // mfem::OperatorChebyshevSmoother or Adams et al. (2003). The estimate starts from the
  // dominant vector of the previous operator, which is usually a nearby one (for example,
  // the next frequency of a driven sweep).
  lambda_max = sf_max * M->EstimateMaxEigenvalue(*A, lambda_vec);
  MFEM_VERIFY(std::isfinite(lambda_max) && lambda_max > 0.0,
              "Encountered invalid maximum eigenvalue in Chebyshev smoother!");

  this->height = op.Height();
  this->width = op.Width();
}

template <OperatorType OperType>
void ChebyshevSmoother<OperType>::Mult2(const VecType &x, VecType &y, VecType &r) const
{
  // Apply smoother: y = y + p(A) (x - A y) .
  for (int it = 0; it < pc_it; it++)
  {
    if (this->initial_guess || it > 0)
    {
      ApplyOp(*A, y, r);
      linalg::AXPBY(1.0, x, -1.0, r);
    }
    else
    {
      r = x;
      y = 0.0;
    }

    Polynomial(y, r);

    // Correction of the relaxation on the residual after the polynomial, if it has one.
    if (M->HasSweepCorrection())
    {
      ApplyOp(*A, y, r);
      linalg::AXPBY(1.0, x, -1.0, r);
      M->AddSweepCorrection(r, y);
    }
  }
}

template <OperatorType OperType>
void ChebyshevSmoother<OperType>::MultTranspose2(const VecType &x, VecType &y,
                                                 VecType &r) const
{
  if (!M->HasSweepCorrection())
  {
    Mult2(x, y, r);  // Assumes operator symmetry
    return;
  }

  // Apply the transposed smoother, which has the correction of the relaxation before the
  // polynomial (both are symmetric).
  for (int it = 0; it < pc_it; it++)
  {
    if (this->initial_guess || it > 0)
    {
      ApplyOp(*A, y, r);
      linalg::AXPBY(1.0, x, -1.0, r);
    }
    else
    {
      r = x;
      y = 0.0;
    }
    M->AddSweepCorrection(r, y);
    ApplyOp(*A, y, r);
    linalg::AXPBY(1.0, x, -1.0, r);
    Polynomial(y, r);
  }
}

template <OperatorType OperType>
void ChebyshevSmoother<OperType>::Polynomial(VecType &y, VecType &r) const
{
  // 4th-kind Chebyshev smoother, from Phillips and Fischer or Lottes (with k -> k + 1
  // shift due to 1-based indexing): y = y + p(M⁻¹A) M⁻¹ r for the residual r of y.
  M->ScaledMult(4.0 / (3.0 * lambda_max), r, d);
  for (int k = 1; k < order; k++)
  {
    y += d;
    ApplyOp(*A, d, r, -1.0);
    const double sd = (2.0 * k - 1.0) / (2.0 * k + 3.0);
    const double sr = (8.0 * k + 4.0) / ((2.0 * k + 3.0) * lambda_max);
    M->ScaledAddMult(sd, sr, r, d);
  }
  y += d;
}

template <OperatorType OperType>
ChebyshevSmoother1stKind<OperType>::ChebyshevSmoother1stKind(MPI_Comm comm, int smooth_it,
                                                             int poly_order, double sf_max,
                                                             double sf_min)
  : ChebyshevSmoother1stKind(std::make_unique<JacobiSmoother<OperType>>(comm), smooth_it,
                             poly_order, sf_max, sf_min)
{
}

template <OperatorType OperType>
ChebyshevSmoother1stKind<OperType>::ChebyshevSmoother1stKind(
    std::unique_ptr<Relaxation<OperType>> &&M, int smooth_it, int poly_order, double sf_max,
    double sf_min)
  : Solver<OperType>(), pc_it(smooth_it), order(poly_order), A(nullptr), M(std::move(M)),
    theta(0.0), sf_max(sf_max), sf_min(sf_min)
{
  MFEM_VERIFY(this->M, "Chebyshev smoothing requires a relaxation method!");
  MFEM_VERIFY(order > 0, "Polynomial order for Chebyshev smoothing must be positive!");
}

template <OperatorType OperType>
void ChebyshevSmoother1stKind<OperType>::SetOperator(const OperType &op)
{
  A = &op;
  d.SetSize(op.Height());
  d.UseDevice(true);
  M->SetOperator(op);

  // Set up Chebyshev coefficients using the computed maximum eigenvalue estimate. The
  // optimized estimate of lambda_min comes from (2.24) of Phillips and Fischer (2022). See
  // ChebyshevSmoother::SetOperator for the starting vector.
  if (sf_min <= 0.0)
  {
    sf_min = 1.69 / (std::pow(order, 1.68) + 2.11 * order + 1.98);
  }
  const double lambda_max = sf_max * M->EstimateMaxEigenvalue(*A, lambda_vec);
  MFEM_VERIFY(std::isfinite(lambda_max) && lambda_max > 0.0,
              "Encountered invalid maximum eigenvalue in Chebyshev smoother!");
  const double lambda_min = sf_min * lambda_max;
  theta = 0.5 * (lambda_max + lambda_min);
  delta = 0.5 * (lambda_max - lambda_min);

  this->height = op.Height();
  this->width = op.Width();
}

template <OperatorType OperType>
void ChebyshevSmoother1stKind<OperType>::Mult2(const VecType &x, VecType &y,
                                               VecType &r) const
{
  // Apply smoother: y = y + p(A) (x - A y) .
  for (int it = 0; it < pc_it; it++)
  {
    if (this->initial_guess || it > 0)
    {
      ApplyOp(*A, y, r);
      linalg::AXPBY(1.0, x, -1.0, r);
    }
    else
    {
      r = x;
      y = 0.0;
    }

    Polynomial(y, r);

    // Correction of the relaxation on the residual after the polynomial, if it has one.
    if (M->HasSweepCorrection())
    {
      ApplyOp(*A, y, r);
      linalg::AXPBY(1.0, x, -1.0, r);
      M->AddSweepCorrection(r, y);
    }
  }
}

template <OperatorType OperType>
void ChebyshevSmoother1stKind<OperType>::MultTranspose2(const VecType &x, VecType &y,
                                                        VecType &r) const
{
  if (!M->HasSweepCorrection())
  {
    Mult2(x, y, r);  // Assumes operator symmetry
    return;
  }

  // Apply the transposed smoother, which has the correction of the relaxation before the
  // polynomial (both are symmetric).
  for (int it = 0; it < pc_it; it++)
  {
    if (this->initial_guess || it > 0)
    {
      ApplyOp(*A, y, r);
      linalg::AXPBY(1.0, x, -1.0, r);
    }
    else
    {
      r = x;
      y = 0.0;
    }
    M->AddSweepCorrection(r, y);
    ApplyOp(*A, y, r);
    linalg::AXPBY(1.0, x, -1.0, r);
    Polynomial(y, r);
  }
}

template <OperatorType OperType>
void ChebyshevSmoother1stKind<OperType>::Polynomial(VecType &y, VecType &r) const
{
  // 1th-kind Chebyshev smoother, from Phillips and Fischer or Adams: y = y + p(M⁻¹A) M⁻¹ r
  // for the residual r of y.
  M->ScaledMult(1.0 / theta, r, d);
  double rhop = delta / theta;
  for (int k = 1; k < order; k++)
  {
    y += d;
    ApplyOp(*A, d, r, -1.0);
    const double rho = 1.0 / (2.0 * theta / delta - rhop);
    const double sd = rho * rhop;
    const double sr = 2.0 * rho / delta;
    M->ScaledAddMult(sd, sr, r, d);
    rhop = rho;
  }
  y += d;
}

template class ChebyshevSmoother<Operator>;
template class ChebyshevSmoother<ComplexOperator>;

template class ChebyshevSmoother1stKind<Operator>;
template class ChebyshevSmoother1stKind<ComplexOperator>;

}  // namespace palace
