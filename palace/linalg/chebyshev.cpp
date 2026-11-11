// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "chebyshev.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <utility>
#include <vector>
#include <mfem/general/forall.hpp>
#include "linalg/rap.hpp"
#include "linalg/solverstats.hpp"
#include "utils/perfstats.hpp"

namespace palace
{

namespace
{

bool HasFiniteNonzeroDiagonal(MPI_Comm comm, const Vector &real,
                              const Vector *imag = nullptr)
{
  int valid = 1;
  if (real.Size() > 0)
  {
    const auto *R = real.HostRead();
    const auto *I = imag ? imag->HostRead() : nullptr;
    for (int i = 0; i < real.Size(); i++)
    {
      if (!(std::isfinite(R[i]) && R[i] != 0.0) ||
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

double LambdaMax(MPI_Comm comm, const Operator &A, const Vector &dinv, Vector &x0)
{
  // D⁻¹A is generally not Hermitian. With Σ = sgn(D), it is similar to ΣS for the
  // Hermitian S = |D|⁻¹ᐟ²A|D|⁻¹ᐟ², so ρ(D⁻¹A) ≤ ‖S‖₂, with equality when Σ commutes with S:
  // for a diagonal of uniform sign (A definite of either sign), and for a negative definite
  // A with essential DOFs eliminated using a unit diagonal.
  Vector dinv_sqrt(dinv);
  dinv_sqrt.Abs();
  linalg::Sqrt(dinv_sqrt);
  DiagonalOperator DinvSqrt(dinv_sqrt);
  ProductOperator ADinvSqrt(A, DinvSqrt);
  ProductOperator S(DinvSqrt, ADinvSqrt);
  return linalg::SpectralNorm(comm, S, true, &x0);
}

double GetLambdaMax(MPI_Comm comm, const Operator &A, const Vector &dinv, Vector &x0)
{
  perf::ScopedTally t(perf::LambdaMaxTally());
  MFEM_VERIFY(HasFiniteNonzeroDiagonal(comm, dinv),
              "Chebyshev smoother spectral estimation requires a finite, nonzero operator "
              "diagonal!");
  return LambdaMax(comm, A, dinv, x0);
}

double GetLambdaMax(MPI_Comm comm, const ComplexOperator &A, const ComplexVector &dinv,
                    ComplexVector &x0)
{
  perf::ScopedTally t(perf::LambdaMaxTally());
  if (A.IsReal())
  {
    MFEM_VERIFY(HasFiniteNonzeroDiagonal(comm, dinv.Real(), &dinv.Imag()),
                "Chebyshev smoother spectral estimation requires a finite, nonzero real "
                "operator diagonal!");

    // Estimate in real arithmetic, which halves the cost, starting from and returning the
    // real part of x0.
    const int n = A.Height();
    Vector x0r;
    if (x0.Size() == n)
    {
      x0r.SetSize(n);
      x0r.UseDevice(true);
      x0r = x0.Real();
    }
    const double lambda_max = LambdaMax(comm, *A.Real(), dinv.Real(), x0r);
    if (x0r.Size() == n)
    {
      x0.SetSize(n);
      x0.UseDevice(true);
      x0.Real() = x0r;
      x0.Imag() = 0.0;
    }
    return lambda_max;
  }
  ComplexDiagonalOperator Dinv(dinv);
  ComplexProductOperator DinvA(Dinv, A);
  return linalg::SpectralNorm(comm, DinvA, false, &x0);
}

// Smallest and largest eigenvalues of the symmetric tridiagonal matrix with diagonal a and
// off-diagonal b, by bisection on Sturm sequence counts.
std::pair<double, double> TridiagonalExtremeEigenvalues(const std::vector<double> &a,
                                                        const std::vector<double> &b)
{
  const int n = static_cast<int>(a.size());
  double lo = a[0], hi = a[0], bmax = 0.0;
  for (int i = 0; i < n; i++)
  {
    const double bl = (i > 0) ? std::abs(b[i - 1]) : 0.0;
    const double br = (i < n - 1) ? std::abs(b[i]) : 0.0;
    lo = std::min(lo, a[i] - bl - br);
    hi = std::max(hi, a[i] + bl + br);
    bmax = std::max(bmax, br);
  }
  const double pivmin = std::numeric_limits<double>::min() * std::max(1.0, bmax * bmax);
  const double tol = 1.0e-8 * std::max(std::abs(lo), std::abs(hi));
  auto Bisect = [&](int k)
  {
    // Smallest x with at least k eigenvalues less than x (the k-th smallest eigenvalue).
    double l = lo, h = hi;
    for (int it = 0; it < 200 && h - l > tol; it++)
    {
      const double x = 0.5 * (l + h);
      int count = 0;
      double q = 1.0;
      for (int i = 0; i < n; i++)
      {
        q = a[i] - x - ((i > 0) ? b[i - 1] * b[i - 1] / q : 0.0);
        if (std::abs(q) < pivmin)
        {
          q = -pivmin;
        }
        count += (q < 0.0);
      }
      (count >= k ? h : l) = x;
    }
    return 0.5 * (l + h);
  };
  return {Bisect(1), Bisect(n)};
}

// Estimate the indefiniteness ratio ρ = -λₘᵢₙ(S) / λₘₐₓ(S) of S = |D|⁻¹ᐟ² A |D|⁻¹ᐟ² for a
// symmetric operator A, given |D|⁻¹, with at most max_it Lanczos iterations. The Ritz
// values lie inside the spectrum of S, so the estimate increases with the number of
// iterations; it stops as soon as it exceeds rho_max.
double EstimateIndefiniteness(MPI_Comm comm, const Operator &A, const Vector &dinv_abs,
                              double rho_max, int max_it)
{
  Vector dinv_sqrt(dinv_abs);
  linalg::Sqrt(dinv_sqrt);
  DiagonalOperator DinvSqrt(dinv_sqrt);
  ProductOperator ADinvSqrt(A, DinvSqrt);
  ProductOperator S(DinvSqrt, ADinvSqrt);
  Vector v(A.Height()), v0(A.Height()), w(A.Height());
  v.UseDevice(true);
  v0.UseDevice(true);
  w.UseDevice(true);
  linalg::SetRandom(comm, v);
  linalg::Normalize(comm, v);
  std::vector<double> alpha, beta;
  double rho = 0.0;
  for (int it = 0; it < max_it; it++)
  {
    S.Mult(v, w);
    alpha.push_back(linalg::Dot(comm, w, v));
    linalg::AXPBYPCZ(-alpha.back(), v, (it > 0) ? -beta.back() : 0.0, v0, 1.0, w);
    const auto [lambda_min, lambda_max] = TridiagonalExtremeEigenvalues(alpha, beta);
    rho = (lambda_max > 0.0) ? -lambda_min / lambda_max
                             : std::numeric_limits<double>::infinity();
    const double norm = linalg::Norml2(comm, w);
    if (!(rho <= rho_max) || norm <= 1.0e-12 * std::abs(lambda_max))
    {
      break;
    }
    beta.push_back(norm);
    v0 = v;
    v.Set(1.0 / norm, w);
  }
  return rho;
}

double GetIndefiniteness(MPI_Comm comm, const Operator &A, const Vector &dinv,
                         double rho_max, int max_it)
{
  Vector dinv_abs(dinv);
  dinv_abs.Abs();
  return EstimateIndefiniteness(comm, A, dinv_abs, rho_max, max_it);
}

double GetIndefiniteness(MPI_Comm comm, const ComplexOperator &A, const ComplexVector &dinv,
                         double rho_max, int max_it)
{
  // Use the Hermitian part of A (the real part for complex symmetric A), scaled with the
  // magnitude of the complex diagonal.
  if (!A.Real())
  {
    return std::numeric_limits<double>::infinity();
  }
  ComplexVector dinv_abs(dinv);
  dinv_abs.Abs();
  return EstimateIndefiniteness(comm, *A.Real(), dinv_abs.Real(), rho_max, max_it);
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
  const bool use_dev = dinv.UseDevice() || r.UseDevice() || d.UseDevice();
  const int N = d.Size();
  const auto *DI = dinv.Read(use_dev);
  const auto *R = r.Read(use_dev);
  auto *D = d.Write(use_dev);
  mfem::forall_switch(use_dev, N,
                      [=] MFEM_HOST_DEVICE(int i) { D[i] = sr * DI[i] * R[i]; });
}

template <bool Transpose = false>
inline void ApplyOrder0(const double sr, const ComplexVector &dinv, const ComplexVector &r,
                        ComplexVector &d)
{
  const bool use_dev = dinv.UseDevice() || r.UseDevice() || d.UseDevice();
  const int N = dinv.Size();
  const auto *DIR = dinv.Real().Read(use_dev);
  const auto *DII = dinv.Imag().Read(use_dev);
  const auto *RR = r.Real().Read(use_dev);
  const auto *RI = r.Imag().Read(use_dev);
  auto *DR = d.Real().Write(use_dev);
  auto *DI = d.Imag().Write(use_dev);
  if constexpr (!Transpose)
  {
    mfem::forall_switch(use_dev, N,
                        [=] MFEM_HOST_DEVICE(int i)
                        {
                          DR[i] = sr * (DIR[i] * RR[i] - DII[i] * RI[i]);
                          DI[i] = sr * (DII[i] * RR[i] + DIR[i] * RI[i]);
                        });
  }
  else
  {
    mfem::forall_switch(use_dev, N,
                        [=] MFEM_HOST_DEVICE(int i)
                        {
                          DR[i] = sr * (DIR[i] * RR[i] + DII[i] * RI[i]);
                          DI[i] = sr * (-DII[i] * RR[i] + DIR[i] * RI[i]);
                        });
  }
}

// Update the iterate with the previous search direction and compute the next one in a
// single pass: y += d, then d = sd d + sr D⁻¹ r.
template <bool Transpose = false>
inline void ApplyOrderK(const double sd, const double sr, const Vector &dinv,
                        const Vector &r, Vector &d, Vector &y)
{
  const bool use_dev = dinv.UseDevice() || r.UseDevice() || d.UseDevice() || y.UseDevice();
  const int N = dinv.Size();
  const auto *DI = dinv.Read(use_dev);
  const auto *R = r.Read(use_dev);
  auto *D = d.ReadWrite(use_dev);
  auto *Y = y.ReadWrite(use_dev);
  mfem::forall_switch(use_dev, N,
                      [=] MFEM_HOST_DEVICE(int i)
                      {
                        Y[i] += D[i];
                        D[i] = sd * D[i] + sr * DI[i] * R[i];
                      });
}

template <bool Transpose = false>
inline void ApplyOrderK(const double sd, const double sr, const ComplexVector &dinv,
                        const ComplexVector &r, ComplexVector &d, ComplexVector &y)
{
  const bool use_dev = dinv.UseDevice() || r.UseDevice() || d.UseDevice() || y.UseDevice();
  const int N = dinv.Size();
  const auto *DIR = dinv.Real().Read(use_dev);
  const auto *DII = dinv.Imag().Read(use_dev);
  const auto *RR = r.Real().Read(use_dev);
  const auto *RI = r.Imag().Read(use_dev);
  auto *DR = d.Real().ReadWrite(use_dev);
  auto *DI = d.Imag().ReadWrite(use_dev);
  auto *YR = y.Real().ReadWrite(use_dev);
  auto *YI = y.Imag().ReadWrite(use_dev);
  if constexpr (!Transpose)
  {
    mfem::forall_switch(use_dev, N,
                        [=] MFEM_HOST_DEVICE(int i)
                        {
                          YR[i] += DR[i];
                          YI[i] += DI[i];
                          DR[i] = sd * DR[i] + sr * (DIR[i] * RR[i] - DII[i] * RI[i]);
                          DI[i] = sd * DI[i] + sr * (DII[i] * RR[i] + DIR[i] * RI[i]);
                        });
  }
  else
  {
    mfem::forall_switch(use_dev, N,
                        [=] MFEM_HOST_DEVICE(int i)
                        {
                          YR[i] += DR[i];
                          YI[i] += DI[i];
                          DR[i] = sd * DR[i] + sr * (DIR[i] * RR[i] + DII[i] * RI[i]);
                          DI[i] = sd * DI[i] + sr * (-DII[i] * RR[i] + DIR[i] * RI[i]);
                        });
  }
}

// Apply the operator to the search direction d for the residual update r -= A d. A
// correction z of the iterate pending in the residual is applied together with it.
template <typename Op, typename VecType>
inline void ApplyDirection(const Op &A_dir, const VecType &d, const VecType *z, VecType &dz,
                           VecType &r)
{
  if (z)
  {
    if (dz.Size() != d.Size())
    {
      dz.SetSize(d.Size());
      dz.UseDevice(true);
    }
    linalg::AXPBYPCZ(1.0, d, 1.0, *z, 0.0, dz);
    ApplyOp(A_dir, dz, r, -1.0);
  }
  else
  {
    ApplyOp(A_dir, d, r, -1.0);
  }
}

}  // namespace

const Operator &GetDirectionOperator(const Operator &op)
{
  const auto *par = dynamic_cast<const ParOperator *>(&op);
  const auto *single = par ? par->GetSinglePrecisionOperator() : nullptr;
  return single ? *single : op;
}

const ComplexOperator &GetDirectionOperator(const ComplexOperator &op)
{
  const auto *par = dynamic_cast<const ComplexParOperator *>(&op);
  const auto *single = par ? par->GetSinglePrecisionOperator() : nullptr;
  return single ? *single : op;
}

template <OperatorType OperType>
ChebyshevSmoother<OperType>::ChebyshevSmoother(MPI_Comm comm, int smooth_it, int poly_order,
                                               double sf_max, int max_poly_order)
  : Solver<OperType>(), comm(comm), pc_it(smooth_it), min_order(poly_order),
    max_order(std::max(poly_order, max_poly_order)), order(poly_order), A(nullptr),
    A_dir(nullptr), lambda_max(0.0), sf_max(sf_max)
{
  MFEM_VERIFY(order > 0, "Polynomial order for Chebyshev smoothing must be positive!");
}

template <OperatorType OperType>
void ChebyshevSmoother<OperType>::SetOperator(const OperType &op)
{
  A = &op;
  A_dir = &GetDirectionOperator(op);
  d.SetSize(op.Height());
  dinv.SetSize(op.Height());
  d.UseDevice(true);
  dinv.UseDevice(true);
  op.AssembleDiagonal(dinv);
  dinv.Reciprocal();

  // Set up Chebyshev coefficients using the computed maximum eigenvalue estimate. See
  // mfem::OperatorChebyshevSmoother or Adams et al. (2003). The estimate starts from the
  // dominant vector of the previous operator, which is usually a nearby one (for example,
  // the next frequency of a driven sweep).
  lambda_max = sf_max * GetLambdaMax(comm, *A, dinv, lambda_vec);
  MFEM_VERIFY(std::isfinite(lambda_max) && lambda_max > 0.0,
              "Encountered invalid maximum eigenvalue in Chebyshev smoother!");

  // Use the maximum order only when the operator is nearly positive semi-definite: a
  // smoothing step of order k amplifies the error at the most negative eigenvalue -ρ λₘₐₓ
  // of D⁻¹ A by about cosh(k √(2ρ)).
  if (max_order > min_order)
  {
    constexpr double max_rho = 1.0e-2;
    constexpr int max_it = 30;
    order = (GetIndefiniteness(comm, op, dinv, max_rho, max_it) <= max_rho) ? max_order
                                                                            : min_order;
  }

  this->height = op.Height();
  this->width = op.Width();
}

template <OperatorType OperType>
void ChebyshevSmoother<OperType>::SetOrder(int poly_order)
{
  MFEM_VERIFY(poly_order > 0, "Polynomial order for Chebyshev smoothing must be positive!");
  order = poly_order;
}

template <OperatorType OperType>
void ChebyshevSmoother<OperType>::Mult2(const VecType &x, VecType &y, VecType &r) const
{
  // Apply smoother: y = y + p(A) (x - A y). With a residual input or output, the residual
  // is updated through the applications to the search directions, also from one iteration
  // to the next.
  const bool update = this->residual_in || this->residual_out;
  for (int it = 0; it < pc_it; it++)
  {
    if (update && (it > 0 || this->residual_in))
    {
      // The residual of the iterate is given, or updated by the previous iteration.
    }
    else if (this->initial_guess || it > 0)
    {
      {
        perf::ScopedTally t(perf_op);
        ApplyOp(*A, y, r);
      }
      linalg::AXPBY(1.0, x, -1.0, r);
    }
    else
    {
      r = x;
      y = 0.0;
    }
    const VecType *z = (it == 0) ? this->residual_pending : nullptr;
    if (z)
    {
      y += *z;
    }

    // 4th-kind Chebyshev smoother, from Phillips and Fischer or Lottes (with k -> k + 1
    // shift due to 1-based indexing).
    ApplyOrder0(4.0 / (3.0 * lambda_max), dinv, r, d);
    for (int k = 1; k < order; k++)
    {
      {
        perf::ScopedTally t(perf_op);
        ApplyDirection(*A_dir, d, (k == 1) ? z : nullptr, dz, r);
      }
      const double sd = (2.0 * k - 1.0) / (2.0 * k + 3.0);
      const double sr = (8.0 * k + 4.0) / ((2.0 * k + 3.0) * lambda_max);
      ApplyOrderK(sd, sr, dinv, r, d, y);
    }
    if (update && (it < pc_it - 1 || this->residual_out))
    {
      perf::ScopedTally t(perf_op);
      ApplyDirection(*A_dir, d, (order == 1) ? z : nullptr, dz, r);
    }
    y += d;
  }
}

template <OperatorType OperType>
ChebyshevSmoother1stKind<OperType>::ChebyshevSmoother1stKind(MPI_Comm comm, int smooth_it,
                                                             int poly_order, double sf_max,
                                                             double sf_min)
  : Solver<OperType>(), comm(comm), pc_it(smooth_it), order(poly_order), A(nullptr),
    A_dir(nullptr), theta(0.0), sf_max(sf_max), sf_min(sf_min)
{
  MFEM_VERIFY(order > 0, "Polynomial order for Chebyshev smoothing must be positive!");
}

template <OperatorType OperType>
void ChebyshevSmoother1stKind<OperType>::SetOperator(const OperType &op)
{
  A = &op;
  A_dir = &GetDirectionOperator(op);
  d.SetSize(op.Height());
  dinv.SetSize(op.Height());
  d.UseDevice(true);
  dinv.UseDevice(true);
  op.AssembleDiagonal(dinv);
  dinv.Reciprocal();

  // Set up Chebyshev coefficients using the computed maximum eigenvalue estimate. The
  // optimized estimate of lambda_min comes from (2.24) of Phillips and Fischer (2022). See
  // ChebyshevSmoother::SetOperator for the starting vector.
  if (sf_min <= 0.0)
  {
    sf_min = 1.69 / (std::pow(order, 1.68) + 2.11 * order + 1.98);
  }
  const double lambda_max = sf_max * GetLambdaMax(comm, *A, dinv, lambda_vec);
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
  // Apply smoother: y = y + p(A) (x - A y). See ChebyshevSmoother::Mult2 for the residual
  // input and output.
  const bool update = this->residual_in || this->residual_out;
  for (int it = 0; it < pc_it; it++)
  {
    if (update && (it > 0 || this->residual_in))
    {
      // The residual of the iterate is given, or updated by the previous iteration.
    }
    else if (this->initial_guess || it > 0)
    {
      {
        perf::ScopedTally t(perf_op);
        ApplyOp(*A, y, r);
      }
      linalg::AXPBY(1.0, x, -1.0, r);
    }
    else
    {
      r = x;
      y = 0.0;
    }
    const VecType *z = (it == 0) ? this->residual_pending : nullptr;
    if (z)
    {
      y += *z;
    }

    // 1th-kind Chebyshev smoother, from Phillips and Fischer or Adams.
    ApplyOrder0(1.0 / theta, dinv, r, d);
    double rhop = delta / theta;
    for (int k = 1; k < order; k++)
    {
      {
        perf::ScopedTally t(perf_op);
        ApplyDirection(*A_dir, d, (k == 1) ? z : nullptr, dz, r);
      }
      const double rho = 1.0 / (2.0 * theta / delta - rhop);
      const double sd = rho * rhop;
      const double sr = 2.0 * rho / delta;
      ApplyOrderK(sd, sr, dinv, r, d, y);
      rhop = rho;
    }
    if (update && (it < pc_it - 1 || this->residual_out))
    {
      perf::ScopedTally t(perf_op);
      ApplyDirection(*A_dir, d, (order == 1) ? z : nullptr, dz, r);
    }
    y += d;
  }
}

template class ChebyshevSmoother<Operator>;
template class ChebyshevSmoother<ComplexOperator>;

template class ChebyshevSmoother1stKind<Operator>;
template class ChebyshevSmoother1stKind<ComplexOperator>;

}  // namespace palace
