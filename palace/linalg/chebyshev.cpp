// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "chebyshev.hpp"

#include <mfem/general/forall.hpp>
#include "linalg/rap.hpp"

namespace palace
{

namespace
{

void GetInverseDiagonal(const ParOperator &A, Vector &dinv)
{
  dinv.SetSize(A.Height());
  A.AssembleDiagonal(dinv);
  dinv.Reciprocal();
}

void GetInverseDiagonal(const ComplexParOperator &A, ComplexVector &dinv)
{
  MFEM_VERIFY(A.HasReal() || A.HasImag(),
              "Invalid zero ComplexOperator for ChebyshevSmoother!");
  dinv.SetSize(A.Height());
  dinv.SetSize(A.Height());
  dinv = 0.0;
  if (A.HasReal())
  {
    A.Real()->AssembleDiagonal(dinv.Real());
  }
  if (A.HasImag())
  {
    A.Imag()->AssembleDiagonal(dinv.Imag());
  }
  dinv.Reciprocal();
}

double GetLambdaMax(MPI_Comm comm, const Operator &A, const Vector &dinv)
{
  DiagonalOperator Dinv(dinv);
  ProductOperator DinvA(Dinv, A);
  return linalg::SpectralNorm(comm, DinvA, false);
}

double GetLambdaMax(MPI_Comm comm, const ComplexOperator &A, const ComplexVector &dinv)
{
  ComplexDiagonalOperator Dinv(dinv);
  ComplexProductOperator DinvA(Dinv, A);
  return linalg::SpectralNorm(comm, DinvA, false);
}

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

template <bool Transpose = false>
inline void ApplyOrder0(double sr, const Vector &dinv, const Vector &r, Vector &d)
{
  const int N = d.Size();
  const auto *DI = dinv.Read();
  const auto *R = r.Read();
  auto *D = d.ReadWrite();
  mfem::forall(N, [=] MFEM_HOST_DEVICE(int i) { D[i] = sr * DI[i] * R[i]; });
}

template <bool Transpose = false>
inline void ApplyOrder0(const double sr, const ComplexVector &dinv, const ComplexVector &r,
                        ComplexVector &d)
{
  const int N = dinv.Size();
  const auto *DIR = dinv.Real().Read();
  const auto *DII = dinv.Imag().Read();
  const auto *RR = r.Real().Read();
  const auto *RI = r.Imag().Read();
  auto *DR = d.Real().ReadWrite();
  auto *DI = d.Imag().ReadWrite();
  if constexpr (!Transpose)
  {
    mfem::forall(N,
                 [=] MFEM_HOST_DEVICE(int i)
                 {
                   DR[i] = sr * (DIR[i] * RR[i] - DII[i] * RI[i]);
                   DI[i] = sr * (DII[i] * RR[i] + DIR[i] * RI[i]);
                 });
  }
  else
  {
    mfem::forall(N,
                 [=] MFEM_HOST_DEVICE(int i)
                 {
                   DR[i] = sr * (DIR[i] * RR[i] + DII[i] * RI[i]);
                   DI[i] = sr * (-DII[i] * RR[i] + DIR[i] * RI[i]);
                 });
  }
}

template <bool Transpose = false>
inline void ApplyOrderK(const double sd, const double sr, const Vector &dinv,
                        const Vector &r, Vector &d)
{
  const int N = dinv.Size();
  const auto *DI = dinv.Read();
  const auto *R = r.Read();
  auto *D = d.ReadWrite();
  mfem::forall(N, [=] MFEM_HOST_DEVICE(int i) { D[i] = sd * D[i] + sr * DI[i] * R[i]; });
}

template <bool Transpose = false>
inline void ApplyOrderK(const double sd, const double sr, const ComplexVector &dinv,
                        const ComplexVector &r, ComplexVector &d)
{
  const int N = dinv.Size();
  const auto *DIR = dinv.Real().Read();
  const auto *DII = dinv.Imag().Read();
  const auto *RR = r.Real().Read();
  const auto *RI = r.Imag().Read();
  auto *DR = d.Real().ReadWrite();
  auto *DI = d.Imag().ReadWrite();
  if constexpr (!Transpose)
  {
    mfem::forall(N,
                 [=] MFEM_HOST_DEVICE(int i)
                 {
                   DR[i] = sd * DR[i] + sr * (DIR[i] * RR[i] - DII[i] * RI[i]);
                   DI[i] = sd * DI[i] + sr * (DII[i] * RR[i] + DIR[i] * RI[i]);
                 });
  }
  else
  {
    mfem::forall(N,
                 [=] MFEM_HOST_DEVICE(int i)
                 {
                   DR[i] = sd * DR[i] + sr * (DIR[i] * RR[i] + DII[i] * RI[i]);
                   DI[i] = sd * DI[i] + sr * (-DII[i] * RR[i] + DIR[i] * RI[i]);
                 });
  }
}

}  // namespace

template <typename OperType>
ChebyshevSmoother<OperType>::ChebyshevSmoother(int smooth_it, int poly_order, double sf_max)
  : Solver<OperType>(), pc_it(smooth_it), order(poly_order), A(nullptr), lambda_max(0.0),
    sf_max(sf_max)
{
  MFEM_VERIFY(order > 0, "Polynomial order for Chebyshev smoothing must be positive!");
}

template <typename OperType>
void ChebyshevSmoother<OperType>::SetOperator(const OperType &op)
{
  using ParOperType =
      typename std::conditional<std::is_same<OperType, ComplexOperator>::value,
                                ComplexParOperator, ParOperator>::type;

  A = &op;
  r.SetSize(op.Height());
  d.SetSize(op.Height());
  this->height = op.Height();
  this->width = op.Width();

  const auto *PtAP = dynamic_cast<const ParOperType *>(&op);
  MFEM_VERIFY(PtAP,
              "ChebyshevSmoother requires a ParOperator or ComplexParOperator operator!");
  GetInverseDiagonal(*PtAP, dinv);

  // Set up Chebyshev coefficients using the computed maximum eigenvalue estimate. See
  // mfem::OperatorChebyshevSmoother or Adams et al. (2003).
  lambda_max = sf_max * GetLambdaMax(PtAP->GetComm(), *A, dinv);
}

template <typename OperType>
void ChebyshevSmoother<OperType>::Mult(const VecType &x, VecType &y) const
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

    // 4th-kind Chebyshev smoother, from Phillips and Fischer or Lottes (with k -> k + 1
    // shift due to 1-based indexing).
    ApplyOrder0(4.0 / (3.0 * lambda_max), dinv, r, d);
    for (int k = 1; k < order; k++)
    {
      y += d;
      ApplyOp(*A, d, r, -1.0);
      const double sd = (2.0 * k - 1.0) / (2.0 * k + 3.0);
      const double sr = (8.0 * k + 4.0) / ((2.0 * k + 3.0) * lambda_max);
      ApplyOrderK(sd, sr, dinv, r, d);
    }
    y += d;
  }
}

template <typename OperType>
ChebyshevSmoother1stKind<OperType>::ChebyshevSmoother1stKind(int smooth_it, int poly_order,
                                                             double sf_max, double sf_min)
  : Solver<OperType>(), pc_it(smooth_it), order(poly_order), A(nullptr), theta(0.0),
    sf_max(sf_max), sf_min(sf_min)
{
  MFEM_VERIFY(order > 0, "Polynomial order for Chebyshev smoothing must be positive!");
}

template <typename OperType>
void ChebyshevSmoother1stKind<OperType>::SetOperator(const OperType &op)
{
  using ParOperType =
      typename std::conditional<std::is_same<OperType, ComplexOperator>::value,
                                ComplexParOperator, ParOperator>::type;

  A = &op;
  r.SetSize(op.Height());
  d.SetSize(op.Height());

  const auto *PtAP = dynamic_cast<const ParOperType *>(&op);
  MFEM_VERIFY(
      PtAP,
      "ChebyshevSmoother1stKind requires a ParOperator or ComplexParOperator operator!");
  GetInverseDiagonal(*PtAP, dinv);

  // Set up Chebyshev coefficients using the computed maximum eigenvalue estimate. The
  // optimized estimate of lambda_min comes from (2.24) of Phillips and Fischer (2022).
  if (sf_min <= 0.0)
  {
    sf_min = 1.69 / (std::pow(order, 1.68) + 2.11 * order + 1.98);
  }
  const double lambda_max = sf_max * GetLambdaMax(PtAP->GetComm(), *A, dinv);
  const double lambda_min = sf_min * lambda_max;
  theta = 0.5 * (lambda_max + lambda_min);
  delta = 0.5 * (lambda_max - lambda_min);
}

template <typename OperType>
void ChebyshevSmoother1stKind<OperType>::Mult(const VecType &x, VecType &y) const
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

    // 1th-kind Chebyshev smoother, from Phillips and Fischer or Adams.
    ApplyOrder0(1.0 / theta, dinv, r, d);
    double rhop = delta / theta;
    for (int k = 1; k < order; k++)
    {
      y += d;
      ApplyOp(*A, d, r, -1.0);
      const double rho = 1.0 / (2.0 * theta / delta - rhop);
      const double sd = rho * rhop;
      const double sr = 2.0 * rho / delta;
      ApplyOrderK(sd, sr, dinv, r, d);
      rhop = rho;
    }
    y += d;
  }
}

template class ChebyshevSmoother<Operator>;
template class ChebyshevSmoother<ComplexOperator>;

template class ChebyshevSmoother1stKind<Operator>;
template class ChebyshevSmoother1stKind<ComplexOperator>;

}  // namespace palace
