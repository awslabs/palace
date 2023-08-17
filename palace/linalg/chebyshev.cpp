// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "chebyshev.hpp"

#include <vector>
#include <numeric>
#include <limits>
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

  // auto lower_bound = 1e-6 * linalg::Norml2(A.GetComm(), dinv);
  auto size = dinv.Size();
  Mpi::GlobalSum(1, &size, A.GetComm());
  auto min = dinv.Min();
  Mpi::GlobalMin(1, &min, A.GetComm());
  auto max = dinv.Max();
  Mpi::GlobalMax(1, &max, A.GetComm());

  auto min_mag = std::transform_reduce(dinv.begin(), dinv.end(),
    std::numeric_limits<double>::max(), [](auto x, auto y){return std::min(x,y); }, [](auto v){return std::abs(v);});
  Mpi::GlobalMin(1, &min_mag, A.GetComm());

  auto lower_bound = 1e-2 * linalg::Norml2(A.GetComm(), dinv) / size;

  Mpi::Print("norm {:.3e}, lower_bound {:.3e}, min {:.3e}, max {:.3e}, min_mag {:.3e}\n",
  linalg::Norml2(A.GetComm(), dinv), lower_bound, min, max, min_mag);

  constexpr double invp = 1.0/3;
  for (auto &x : dinv)
  {
    if (std::abs(x) < lower_bound)
    {
      // x = std::copysign(lower_bound, x);
      x = std::copysign(lower_bound * std::pow(std::abs(x)/lower_bound, invp), x);
    }
  }

  dinv.Reciprocal();
}

void GetInverseDiagonal(const ComplexParOperator &A, Vector &dinv)
{
  MFEM_VERIFY(A.HasReal() && !A.HasImag(),
              "ComplexOperator for ChebyshevSmoother must be real-valued for now!");
  dinv.SetSize(A.Height());
  A.Real()->AssembleDiagonal(dinv);

  // auto lower_bound = 1e-6 * linalg::Norml2(A.GetComm(), dinv);
  auto size = dinv.Size();
  Mpi::GlobalSum(1, &size, A.GetComm());
  auto min = dinv.Min();
  Mpi::GlobalMin(1, &min, A.GetComm());
  auto max = dinv.Max();
  Mpi::GlobalMax(1, &max, A.GetComm());

  auto min_mag = std::transform_reduce(dinv.begin(), dinv.end(),
    std::numeric_limits<double>::max(), [](auto x, auto y){return std::min(x,y); }, [](auto v){return std::abs(v);});
  Mpi::GlobalMin(1, &min_mag, A.GetComm());

  auto lower_bound = linalg::Norml2(A.GetComm(), dinv) / size;

  Mpi::Print("norm {:.3e}, lower_bound {:.3e}, min {:.3e}, max {:.3e}, min_mag {:.3e}\n",
  linalg::Norml2(A.GetComm(), dinv), lower_bound, min, max, min_mag);

  constexpr double invp = 1.0/3;
  for (auto &x : dinv)
  {
    if (std::abs(x) < lower_bound)
    {
      // x = std::copysign(lower_bound, x);
      x = std::copysign(lower_bound * std::pow(std::abs(x)/lower_bound, invp), x);
    }
  }

  dinv.Reciprocal();
  // dinv *= 0.5;
  // MFEM_VERIFY(A.HasReal() || A.HasImag(),
  //             "Invalid zero ComplexOperator for ChebyshevSmoother!");
  // dinv.SetSize(A.Height());
  // dinv.SetSize(A.Height());
  // dinv = 0.0;
  // if (A.HasReal())
  // {
  //   A.Real()->AssembleDiagonal(dinv.Real());
  // }
  // if (A.HasImag())
  // {
  //   A.Imag()->AssembleDiagonal(dinv.Imag());
  // }
  // dinv.Reciprocal();
}

double GetLambdaMax(MPI_Comm comm, const Operator &A, const Vector &dinv)
{
  DiagonalOperator Dinv(dinv);
  ProductOperator DinvA(Dinv, A);
  return linalg::SpectralNorm(comm, DinvA, false);
}

double GetLambdaMax(MPI_Comm comm, const ComplexOperator &A, const Vector &dinv)
{
  MFEM_VERIFY(A.HasReal() && !A.HasImag(),
              "ComplexOperator for ChebyshevSmoother must be real-valued for now!");
  DiagonalOperator Dinv(dinv);
  ProductOperator DinvA(Dinv, *A.Real());
  return linalg::SpectralNorm(comm, DinvA, false);
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
  using ParOperType =
      typename std::conditional<std::is_same<OperType, ComplexOperator>::value,
                                ComplexParOperator, ParOperator>::type;

  A = &op;
  r.SetSize(op.Height());
  d.SetSize(op.Height());

  const auto *PtAP = dynamic_cast<const ParOperType *>(&op);
  MFEM_VERIFY(PtAP,
              "ChebyshevSmoother requires a ParOperator or ComplexParOperator operator!");
  GetInverseDiagonal(*PtAP, dinv);

  // Set up Chebyshev coefficients using the computed maximum eigenvalue estimate. See
  // mfem::OperatorChebyshevSmoother or Adams et al., Parallel multigrid smoothing:
  // polynomial versus Gauss-Seidel, JCP (2003).
  lambda_max = 1.01 * GetLambdaMax(PtAP->GetComm(), *A, dinv);
}

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
inline void ApplyOrder0(const double sr, const Vector &dinv, const ComplexVector &r,
                        ComplexVector &d)
{
  const int N = dinv.Size();
  // const auto *DIR = dinv.Real().Read();
  // const auto *DII = dinv.Imag().Read();
  const auto *DIR = dinv.Read();
  const auto *RR = r.Real().Read();
  const auto *RI = r.Imag().Read();
  auto *DR = d.Real().ReadWrite();
  auto *DI = d.Imag().ReadWrite();
  if constexpr (!Transpose)
  {
    mfem::forall(N,
                 [=] MFEM_HOST_DEVICE(int i)
                 {
                   // DR[i] = sr * (DIR[i] * RR[i] - DII[i] * RI[i]);
                   // DI[i] = sr * (DII[i] * RR[i] + DIR[i] * RI[i]);
                   DR[i] = sr * DIR[i] * RR[i];
                   DI[i] = sr * DIR[i] * RI[i];
                 });
  }
  else
  {
    mfem::forall(N,
                 [=] MFEM_HOST_DEVICE(int i)
                 {
                   // DR[i] = sr * (DIR[i] * RR[i] + DII[i] * RI[i]);
                   // DI[i] = sr * (-DII[i] * RR[i] + DIR[i] * RI[i]);
                   DR[i] = sr * DIR[i] * RR[i];
                   DI[i] = sr * DIR[i] * RI[i];
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
inline void ApplyOrderK(const double sd, const double sr, const Vector &dinv,
                        const ComplexVector &r, ComplexVector &d)
{
  const int N = dinv.Size();
  // const auto *DIR = dinv.Real().Read();
  // const auto *DII = dinv.Imag().Read();
  const auto *DIR = dinv.Read();
  const auto *RR = r.Real().Read();
  const auto *RI = r.Imag().Read();
  auto *DR = d.Real().ReadWrite();
  auto *DI = d.Imag().ReadWrite();
  if constexpr (!Transpose)
  {
    mfem::forall(N,
                 [=] MFEM_HOST_DEVICE(int i)
                 {
                   // DR[i] = sd * DR[i] + sr * (DIR[i] * RR[i] - DII[i] * RI[i]);
                   // DI[i] = sd * DI[i] + sr * (DII[i] * RR[i] + DIR[i] * RI[i]);
                   DR[i] = sd * DR[i] + sr * DIR[i] * RR[i];
                   DI[i] = sd * DI[i] + sr * DIR[i] * RI[i];
                 });
  }
  else
  {
    mfem::forall(N,
                 [=] MFEM_HOST_DEVICE(int i)
                 {
                   // DR[i] = sd * DR[i] + sr * (DIR[i] * RR[i] + DII[i] * RI[i]);
                   // DI[i] = sd * DI[i] + sr * (-DII[i] * RR[i] + DIR[i] * RI[i]);
                   DR[i] = sd * DR[i] + sr * DIR[i] * RR[i];
                   DI[i] = sd * DI[i] + sr * DIR[i] * RI[i];
                 });
  }
}

}  // namespace

template <typename OperType>
void ChebyshevSmoother<OperType>::Mult(const VecType &x, VecType &y) const
{
  Mpi::Print("Cheby ||x|| = {:.4e}, ||y|| = {:.4e}\n",
    linalg::Norml2(MPI_COMM_WORLD, x),
    linalg::Norml2(MPI_COMM_WORLD, y));

  // Apply smoother: y = y + p(A) (x - A y) .
  for (int it = 0; it < pc_it; it++)
  {
    if (this->initial_guess || it > 0)
    {
      ApplyOp(*A, y, r);
      Mpi::Print("Cheby ||Ay|| = {:.4e}\n",
        linalg::Norml2(MPI_COMM_WORLD, r));

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

    Mpi::Print("Cheby ApplyOrder0 lambda_max = {:.4e}, ||dinv|| = {:.4e}, ||r|| = {:.4e}, ||d|| = {:.4e}\n", lambda_max,
    linalg::Norml2(MPI_COMM_WORLD, dinv), linalg::Norml2(MPI_COMM_WORLD, r), linalg::Norml2(MPI_COMM_WORLD, d));

    for (int k = 1; k < order; k++)
    {
      y += d;
      ApplyOp(*A, d, r, -1.0);

      Mpi::Print("Cheby ApplyOp, ||y|| = {:.4e}, ||r|| = {:.4e}, \n",
      linalg::Norml2(MPI_COMM_WORLD, y), linalg::Norml2(MPI_COMM_WORLD, r));

      const double sd = (2.0 * k - 1.0) / (2.0 * k + 3.0);
      const double sr = (8.0 * k + 4.0) / ((2.0 * k + 3.0) * lambda_max);
      ApplyOrderK(sd, sr, dinv, r, d);

      Mpi::Print("Cheby ApplyOrder{} sd {:.4e} sr {:.4e}, ||dinv|| = {:.4e}, ||r|| = {:.4e}, ||d|| = {:.4e}\n", k, sd, sr,
      linalg::Norml2(MPI_COMM_WORLD, dinv), linalg::Norml2(MPI_COMM_WORLD, r), linalg::Norml2(MPI_COMM_WORLD, d));

    }
    y += d;
    Mpi::Print("Cheby ||y|| = {:.4e}\n",
      linalg::Norml2(MPI_COMM_WORLD, y));
  }
}

template class ChebyshevSmoother<Operator>;
template class ChebyshevSmoother<ComplexOperator>;

}  // namespace palace
