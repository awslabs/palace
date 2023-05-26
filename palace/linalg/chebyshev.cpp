// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "chebyshev.hpp"

#include <vector>
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

}  // namespace

template <typename OperType>
ChebyshevSmoother<OperType>::ChebyshevSmoother(int smooth_it, int poly_order)
  : Solver<OperType>(), pc_it(smooth_it), order(poly_order), A(nullptr)
{
}

template <typename OperType>
void ChebyshevSmoother<OperType>::SetOperator(const OperType &op)
{
  typedef typename std::conditional<std::is_same<OperType, ComplexOperator>::value,
                                    ComplexParOperator, ParOperator>::type ParOperType;

  const auto *PtAP = dynamic_cast<const ParOperType *>(&op);
  MFEM_VERIFY(PtAP,
              "ChebyshevSmoother requires a ParOperator or ComplexParOperator operator!");
  GetInverseDiagonal(*PtAP, dinv);

  A = &op;
  r.SetSize(A->Height());
  d.SetSize(A->Height());

  // Set up Chebyshev coefficients using the computed maximum eigenvalue estimate. See
  // mfem::OperatorChebyshevSmoother or Adams et al., Parallel multigrid smoothing:
  // polynomial versus Gauss-Seidel, JCP (2003).
  BaseDiagonalOperator<OperType> Dinv(dinv);
  BaseProductOperator<OperType> DinvA(Dinv, *A);
  lambda_max = 1.01 * linalg::SpectralNorm(PtAP->GetComm(), DinvA, false);
}

namespace
{

template <bool Transpose = false>
inline void ApplyOp(const Operator &A, const Vector &x, Vector &y, const double a)
{
  if (a == 0.0)
  {
    A.Mult(x, y);
  }
  else
  {
    A.AddMult(x, y, a);
  }
}

template <bool Transpose = false>
inline void ApplyOp(const ComplexOperator &A, const ComplexVector &x, ComplexVector &y,
                    const double a)
{
  if constexpr (!Transpose)
  {
    if (a == 0.0)
    {
      A.Mult(x, y);
    }
    else
    {
      A.AddMult(x, y, a);
    }
  }
  else
  {
    if (a == 0.0)
    {
      A.MultHermitianTranspose(x, y);
    }
    else
    {
      A.AddMultHermitianTranspose(x, y, a);
    }
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
                   const double t = DII[i] * RR[i] + DIR[i] * RI[i];
                   DR[i] = sr * (DIR[i] * RR[i] - DII[i] * RI[i]);
                   DI[i] = sr * t;
                 });
  }
  else
  {
    mfem::forall(N,
                 [=] MFEM_HOST_DEVICE(int i)
                 {
                   const double t = -DII[i] * RR[i] + DIR[i] * RI[i];
                   DR[i] = sr * (DIR[i] * RR[i] + DII[i] * RI[i]);
                   DI[i] = sr * t;
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
                   const double t = DII[i] * RR[i] + DIR[i] * RI[i];
                   DR[i] = sd * DR[i] + sr * (DIR[i] * RR[i] - DII[i] * RI[i]);
                   DI[i] = sd * DI[i] + sr * t;
                 });
  }
  else
  {
    mfem::forall(N,
                 [=] MFEM_HOST_DEVICE(int i)
                 {
                   const double t = -DII[i] * RR[i] + DIR[i] * RI[i];
                   DR[i] = sd * DR[i] + sr * (DIR[i] * RR[i] + DII[i] * RI[i]);
                   DI[i] = sd * DI[i] + sr * t;
                 });
  }
}

}  // namespace

template <typename OperType>
void ChebyshevSmoother<OperType>::Mult(const VecType &x, VecType &y) const
{
  // Apply smoother: y = y + p(A) (x - A y) .
  for (int it = 0; it < pc_it; it++)
  {
    if (this->initial_guess || it > 0)
    {
      ApplyOp(*A, y, r, 0.0);
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

template class ChebyshevSmoother<Operator>;
template class ChebyshevSmoother<ComplexOperator>;

}  // namespace palace
