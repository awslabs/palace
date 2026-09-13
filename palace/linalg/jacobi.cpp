// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "jacobi.hpp"

#include <cmath>
#include <mfem/general/forall.hpp>

namespace palace
{

namespace
{

bool HasPositiveFiniteDiagonal(MPI_Comm comm, const Vector &real,
                               const Vector *imag = nullptr)
{
  int valid = 1;
  if (real.Size() > 0)
  {
    const auto *R = real.HostRead();
    const auto *I = imag ? imag->HostRead() : nullptr;
    for (int i = 0; i < real.Size(); i++)
    {
      if (!(std::isfinite(R[i]) && R[i] > 0.0) ||
          (I && !(std::isfinite(I[i]) && I[i] == 0.0)))
      {
        valid = 0;
        break;
      }
    }
  }
  Mpi::GlobalMin(1, &valid, comm);
  return valid;
}

double GetLambdaMax(MPI_Comm comm, const Operator &A, const Vector &dinv)
{
  // D⁻¹A is generally not Hermitian, but is similar to the Hermitian operator
  // D⁻¹ᐟ²AD⁻¹ᐟ² when A has a finite, strictly positive diagonal.
  MFEM_VERIFY(HasPositiveFiniteDiagonal(comm, dinv),
              "Jacobi smoother spectral estimation requires a finite, strictly positive "
              "operator diagonal!");
  Vector dinv_sqrt(dinv);
  linalg::Sqrt(dinv_sqrt);
  DiagonalOperator DinvSqrt(dinv_sqrt);
  ProductOperator ADinvSqrt(A, DinvSqrt);
  ProductOperator S(DinvSqrt, ADinvSqrt);
  return linalg::SpectralNorm(comm, S, true);
}

double GetLambdaMax(MPI_Comm comm, const ComplexOperator &A, const ComplexVector &dinv)
{
  if (A.IsReal())
  {
    MFEM_VERIFY(HasPositiveFiniteDiagonal(comm, dinv.Real(), &dinv.Imag()),
                "Jacobi smoother spectral estimation requires a finite, strictly positive "
                "real operator diagonal!");
    ComplexVector dinv_sqrt(dinv);
    linalg::Sqrt(dinv_sqrt.Real());
    ComplexDiagonalOperator DinvSqrt(dinv_sqrt);
    ComplexProductOperator ADinvSqrt(A, DinvSqrt);
    ComplexProductOperator S(DinvSqrt, ADinvSqrt);
    return linalg::SpectralNorm(comm, S, true);
  }
  ComplexDiagonalOperator Dinv(dinv);
  ComplexProductOperator DinvA(Dinv, A);
  return linalg::SpectralNorm(comm, DinvA, false);
}

template <bool Transpose = false>
inline void Apply(const Vector &dinv, const Vector &x, Vector &y)
{
  const bool use_dev = dinv.UseDevice() || x.UseDevice() || y.UseDevice();
  const int N = dinv.Size();
  const auto *DI = dinv.Read(use_dev);
  const auto *X = x.Read(use_dev);
  auto *Y = y.Write(use_dev);
  mfem::forall_switch(use_dev, N, [=] MFEM_HOST_DEVICE(int i) { Y[i] = DI[i] * X[i]; });
}

template <bool Transpose = false>
inline void Apply(const ComplexVector &dinv, const ComplexVector &x, ComplexVector &y)
{
  const bool use_dev = dinv.UseDevice() || x.UseDevice() || y.UseDevice();
  const int N = dinv.Size();
  const auto *DIR = dinv.Real().Read(use_dev);
  const auto *DII = dinv.Imag().Read(use_dev);
  const auto *XR = x.Real().Read(use_dev);
  const auto *XI = x.Imag().Read(use_dev);
  auto *YR = y.Real().Write(use_dev);
  auto *YI = y.Imag().Write(use_dev);
  if constexpr (!Transpose)
  {
    mfem::forall_switch(use_dev, N,
                        [=] MFEM_HOST_DEVICE(int i)
                        {
                          YR[i] = DIR[i] * XR[i] - DII[i] * XI[i];
                          YI[i] = DII[i] * XR[i] + DIR[i] * XI[i];
                        });
  }
  else
  {
    mfem::forall_switch(use_dev, N,
                        [=] MFEM_HOST_DEVICE(int i)
                        {
                          YR[i] = DIR[i] * XR[i] + DII[i] * XI[i];
                          YI[i] = -DII[i] * XR[i] + DIR[i] * XI[i];
                        });
  }
}

}  // namespace

template <typename OperType>
void JacobiSmoother<OperType>::SetOperator(const OperType &op)
{
  dinv.SetSize(op.Height());
  dinv.UseDevice(true);
  op.AssembleDiagonal(dinv);
  dinv.Reciprocal();

  // Damping factor. If the given damping is zero, estimate the spectral radius-minimizing
  // damping factor.
  if (omega == 0.0)
  {
    const auto lambda_max = GetLambdaMax(comm, op, dinv);
    MFEM_VERIFY(std::isfinite(lambda_max) && lambda_max > 0.0,
                "Encountered invalid maximum eigenvalue in Jacobi smoother!");
    const auto lambda_min = (sf_max - 1.0) * lambda_max;
    const auto denominator = lambda_min + lambda_max;
    MFEM_VERIFY(std::isfinite(denominator) && denominator > 0.0,
                "Automatic Jacobi damping requires a finite, strictly positive damping "
                "denominator!");
    omega = 2.0 / denominator;
  }
  if (omega != 1.0)
  {
    dinv *= omega;
  }

  this->height = op.Height();
  this->width = op.Width();
}

template <typename OperType>
void JacobiSmoother<OperType>::Mult(const VecType &x, VecType &y) const
{
  MFEM_ASSERT(!this->initial_guess, "JacobiSmoother does not use initial guess!");
  Apply(dinv, x, y);
}

template class JacobiSmoother<Operator>;
template class JacobiSmoother<ComplexOperator>;

}  // namespace palace
