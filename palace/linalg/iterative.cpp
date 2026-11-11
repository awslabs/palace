// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "iterative.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <limits>
#include <string>
#include <typeinfo>
#include "linalg/jacobi.hpp"
#include "linalg/orthog.hpp"
#include "utils/communication.hpp"
#include "utils/perfstats.hpp"
#include "utils/timer.hpp"

namespace palace
{

namespace
{

// Apply update(i) for all i in [0, N) and return the sum of the terms f(i) after the
// updates, in one pass over memory: each chunk is updated by a simple loop, which compilers
// vectorize, and its terms are then summed from cache with several partial sums.
template <typename UpdateFunction, typename TermFunction>
double UpdateAndSum(int N, UpdateFunction &&update, TermFunction &&f)
{
  constexpr int chunk = 256;
  double sum = 0.0;
  for (int begin = 0; begin < N; begin += chunk)
  {
    const int end = std::min(begin + chunk, N);
    for (int i = begin; i < end; i++)
    {
      update(i);
    }
    sum += linalg::HostSums<1>(begin, end, f)[0];
  }
  return sum;
}

// r -= a z, z = D r, and return the local inner product (z, r), in one pass.
double CgUpdateJacobi(double a, Vector &r, Vector &z, const Vector &dinv)
{
  const int N = r.Size();
  MFEM_VERIFY(z.Size() == N && dinv.Size() == N, "Size mismatch for fused CG update!");
  const auto *D = dinv.HostRead();
  auto *R = r.HostReadWrite(), *Z = z.HostReadWrite();
  return UpdateAndSum(
      N,
      [=](int i)
      {
        R[i] -= a * Z[i];
        Z[i] = D[i] * R[i];
      },
      [=](int i) { return std::array{R[i] * Z[i]}; });
}

// The same for complex vectors and a real D, for which (z, r) = rᴴ D r is real. The step
// size a is complex unless the operator is Hermitian (and then only up to rounding).
std::complex<double> CgUpdateJacobi(std::complex<double> a, ComplexVector &r,
                                    ComplexVector &z, const Vector &dinv)
{
  const int N = r.Size();
  MFEM_VERIFY(z.Size() == N && dinv.Size() == N, "Size mismatch for fused CG update!");
  const double ar = a.real(), ai = a.imag();
  const auto *D = dinv.HostRead();
  auto *RR = r.Real().HostReadWrite(), *RI = r.Imag().HostReadWrite();
  auto *ZR = z.Real().HostReadWrite(), *ZI = z.Imag().HostReadWrite();
  return UpdateAndSum(
      N,
      [=](int i)
      {
        RR[i] -= ar * ZR[i] - ai * ZI[i];
        RI[i] -= ai * ZR[i] + ar * ZI[i];
        ZR[i] = D[i] * RR[i];
        ZI[i] = D[i] * RI[i];
      },
      [=](int i) { return std::array{RR[i] * ZR[i] + RI[i] * ZI[i]}; });
}

inline double Conj(double x)
{
  return x;
}

inline std::complex<double> Conj(std::complex<double> x)
{
  return std::conj(x);
}

template <typename T>
inline void CheckDot(T dot, const char *msg)
{
  MFEM_ASSERT(std::isfinite(dot) && dot >= 0.0, msg << dot << "!");
}

template <typename T>
inline void CheckDot(std::complex<T> dot, const char *msg)
{
  MFEM_ASSERT(std::isfinite(dot.real()) && std::isfinite(dot.imag()) && dot.real() >= 0.0,
              msg << dot << "!");
}

template <typename T>
inline void GeneratePlaneRotation(const T dx, const T dy, T &cs, T &sn)
{
  // See LAPACK's s/dlartg.
  const T safmin = std::numeric_limits<T>::min();
  const T safmax = 1 / safmin;
  const T root_min = std::sqrt(safmin);
  const T root_max = std::sqrt(safmax / 2);
  if (dy == 0.0)
  {
    cs = 1.0;
    sn = 0.0;
    return;
  }
  if (dx == 0.0)
  {
    cs = 0.0;
    sn = std::copysign(1.0, dy);
    return;
  }
  T dx1 = std::abs(dx);
  T dy1 = std::abs(dy);
  if (dx1 > root_min && dx1 < root_max && dy1 > root_min && dy1 < root_max)
  {
    T d = std::sqrt(dx * dx + dy * dy);
    cs = dx1 / d;
    sn = dy / std::copysign(d, dx);
  }
  else
  {
    T u = std::clamp(std::max(dx1, dy1), safmin, safmax);
    T dxs = dx / u;
    T dys = dy / u;
    T d = std::sqrt(dxs * dxs + dys * dys);
    cs = std::abs(dxs) / d;
    sn = dys / std::copysign(d, dx);
  }
}

template <typename T>
inline void GeneratePlaneRotation(const std::complex<T> dx, const std::complex<T> dy, T &cs,
                                  std::complex<T> &sn)
{
  // Generates a plane rotation so that:
  //   [  cs        sn ] [ dx ]  =  [ r ]
  //   [ -conj(sn)  cs ] [ dy ]     [ 0 ]
  // where cs is real and cs² + |sn|² = 1. See LAPACK's c/zlartg.
  const T safmin = std::numeric_limits<T>::min();
  const T safmax = 1 / safmin;
  if (dy == 0.0)
  {
    cs = 1.0;
    sn = 0.0;
    return;
  }
  if (dx == 0.0)
  {
    cs = 0.0;
    if (dy.real() == 0.0)
    {
      sn = std::conj(dy) / std::abs(dy.imag());
    }
    else if (dy.imag() == 0.0)
    {
      sn = std::conj(dy) / std::abs(dy.real());
    }
    else
    {
      const T root_min = std::sqrt(safmin);
      const T root_max = std::sqrt(safmax / 2);
      T dy1 = std::max(std::abs(dy.real()), std::abs(dy.imag()));
      if (dy1 > root_min && dy1 < root_max)
      {
        sn = std::conj(dy) / std::sqrt(dy.real() * dy.real() + dy.imag() * dy.imag());
      }
      else
      {
        T u = std::clamp(dy1, safmin, safmax);
        std::complex<T> dys = dy / u;
        sn = std::conj(dys) / std::sqrt(dys.real() * dys.real() + dys.imag() * dys.imag());
      }
    }
    return;
  }
  const T root_min = std::sqrt(safmin);
  const T root_max = std::sqrt(safmax / 4);
  T dx1 = std::max(std::abs(dx.real()), std::abs(dx.imag()));
  T dy1 = std::max(std::abs(dy.real()), std::abs(dy.imag()));
  if (dx1 > root_min && dx1 < root_max && dy1 > root_min && dy1 < root_max)
  {
    T dx2 = dx.real() * dx.real() + dx.imag() * dx.imag();
    T dy2 = dy.real() * dy.real() + dy.imag() * dy.imag();
    T dz2 = dx2 + dy2;
    if (dx2 >= dz2 * safmin)
    {
      cs = std::sqrt(dx2 / dz2);
      if (dx2 > root_min && dz2 < root_max * 2)
      {
        sn = std::conj(dy) * (dx / std::sqrt(dx2 * dz2));
      }
      else
      {
        sn = std::conj(dy) * ((dx / cs) / dz2);
      }
    }
    else
    {
      T d = std::sqrt(dx2 * dz2);
      cs = dx2 / d;
      sn = std::conj(dy) * (dx / d);
    }
  }
  else
  {
    T u = std::clamp(std::max(dx1, dy1), safmin, safmax), w;
    std::complex<T> dys = dy / u, dxs;
    T dy2 = dys.real() * dys.real() + dys.imag() * dys.imag(), dx2, dz2;
    if (dx1 / u < root_min)
    {
      T v = std::clamp(dx1, safmin, safmax);
      w = v / u;
      dxs = dx / v;
      dx2 = dxs.real() * dxs.real() + dxs.imag() * dxs.imag();
      dz2 = dx2 * w * w + dy2;
    }
    else
    {
      w = 1.0;
      dxs = dx / u;
      dx2 = dxs.real() * dxs.real() + dxs.imag() * dxs.imag();
      dz2 = dx2 + dy2;
    }
    if (dx2 >= dz2 * safmin)
    {
      cs = std::sqrt(dx2 / dz2);
      if (dx2 > root_min && dz2 < root_max * 2)
      {
        sn = std::conj(dys) * (dxs / std::sqrt(dx2 * dz2));
      }
      else
      {
        sn = std::conj(dys) * ((dxs / cs) / dz2);
      }
    }
    else
    {
      T d = std::sqrt(dx2 * dz2);
      cs = dx2 / d;
      sn = std::conj(dys) * (dxs / d);
    }
    cs *= w;
  }
}

template <typename T>
inline void ApplyPlaneRotation(T &dx, T &dy, const T cs, const T sn)
{
  T t = cs * dx + sn * dy;
  dy = -sn * dx + cs * dy;
  dx = t;
}

template <typename T>
inline void ApplyPlaneRotation(std::complex<T> &dx, std::complex<T> &dy, const T cs,
                               const std::complex<T> sn)
{
  std::complex<T> t = cs * dx + sn * dy;
  dy = -std::conj(sn) * dx + cs * dy;
  dx = t;
}

template <typename OperType, typename VecType>
inline void ApplyB(const Solver<OperType> *B, const VecType &x, VecType &y,
                   bool use_timer = true)
{
  BlockTimer bt(Timer::KSP_PRECONDITIONER, use_timer);
  MFEM_ASSERT(B, "Missing preconditioner in ApplyB!");
  B->Mult(x, y);
}

// Local inner products (z, r) and (D r, r) for a real diagonal D, in one pass on the host.
inline std::array<double, 2> LocalDotNorm(const Vector &z, const Vector &r, const Vector &d)
{
  const auto *Z = z.HostRead(), *R = r.HostRead(), *D = d.HostRead();
  return linalg::HostSums<2>(0, r.Size(), [=](int i)
                             { return std::array{R[i] * Z[i], D[i] * R[i] * R[i]}; });
}

inline std::array<std::complex<double>, 2>
LocalDotNorm(const ComplexVector &z, const ComplexVector &r, const Vector &d)
{
  const auto *ZR = z.Real().HostRead(), *ZI = z.Imag().HostRead();
  const auto *RR = r.Real().HostRead(), *RI = r.Imag().HostRead(), *D = d.HostRead();
  const auto s = linalg::HostSums<3>(
      0, r.Size(),
      [=](int i)
      {
        return std::array{RR[i] * ZR[i] + RI[i] * ZI[i], RR[i] * ZI[i] - RI[i] * ZR[i],
                          D[i] * (RR[i] * RR[i] + RI[i] * RI[i])};
      });
  return {std::complex<double>(s[0], s[1]), s[2]};
}

template <typename OperType, typename VecType>
inline void InitialResidual(PreconditionerSide side, const OperType *A,
                            const Solver<OperType> *B, const VecType &b, VecType &x,
                            VecType &r, VecType &z, bool initial_guess,
                            bool use_timer = true)
{
  if (B && side == PreconditionerSide::LEFT)
  {
    if (initial_guess)
    {
      A->Mult(x, z);
      linalg::AXPBY(1.0, b, -1.0, z);
      ApplyB(B, z, r, use_timer);
    }
    else
    {
      ApplyB(B, b, r, use_timer);
      x = 0.0;
    }
  }
  else  // !B || side == PreconditionerSide::RIGHT
  {
    if (initial_guess)
    {
      A->Mult(x, r);
      linalg::AXPBY(1.0, b, -1.0, r);
    }
    else
    {
      r = b;
      x = 0.0;
    }
  }
}

template <typename OperType, typename VecType>
inline void ApplyBA(PreconditionerSide side, const OperType *A, const Solver<OperType> *B,
                    const VecType &x, VecType &y, VecType &z, bool use_timer = true)
{
  if (B && side == PreconditionerSide::LEFT)
  {
    A->Mult(x, z);
    ApplyB(B, z, y, use_timer);
  }
  else if (B && side == PreconditionerSide::RIGHT)
  {
    ApplyB(B, x, z, use_timer);
    A->Mult(z, y);
  }
  else
  {
    A->Mult(x, y);
  }
}

}  // namespace

template <OperatorType OperType>
IterativeSolver<OperType>::IterativeSolver(MPI_Comm comm, int print)
  : Solver<OperType>(), comm(comm), A(nullptr), B(nullptr), final_restarts(0)
{
  print_opts.Warnings();
  if (print > 0)
  {
    print_opts.Summary();
    if (print > 1)
    {
      print_opts.Iterations();
      if (print > 2)
      {
        print_opts.All();
      }
    }
  }
  int_width = 3;
  tab_width = 0;

  rel_tol = abs_tol = 0.0;
  max_it = 100;

  converged = false;
  initial_res = 1.0;
  final_res = 0.0;
  final_it = 0;

  use_timer = false;
}

template <OperatorType OperType>
void IterativeSolver<OperType>::SetPerfTallies(const std::string &prefix,
                                               Timer::Index parent)
{
  perf_apply = perf::RegisterTally(prefix + "/ApplyBA", parent);
  perf_orth = perf::RegisterTally(prefix + "/Orthogonalize", parent);
}

template <OperatorType OperType>
void CgSolver<OperType>::Mult(const VecType &b, VecType &x) const
{
  // Set up workspace.
  ScalarType beta, beta_prev = 0.0, alpha, denom;
  RealType res, eps;
  MFEM_VERIFY(A, "Operator must be set for CgSolver::Mult!");
  MFEM_ASSERT(A->Width() == x.Size() && A->Height() == b.Size(),
              "Size mismatch for CgSolver::Mult!");
  r.SetSize(A->Height());
  z.SetSize(A->Height());
  p.SetSize(A->Height());
  r.UseDevice(true);
  z.UseDevice(true);
  p.UseDevice(true);

  // A Jacobi preconditioner may be applied in the same pass as the residual update and the
  // inner product (not a subclass, which may change the application), unless the stopping
  // test uses the Jacobi norm, which the unfused pass computes.
  const Vector *dinv = nullptr;
  if (B && !dinv_norm && typeid(*B) == typeid(JacobiSmoother<OperType>))
  {
    dinv = static_cast<const JacobiSmoother<OperType> *>(B)->GetHostRealScaling();
  }

  // Initialize.
  if (this->initial_guess)
  {
    A->Mult(x, r);
    linalg::AXPBY(1.0, b, -1.0, r);
  }
  else
  {
    r = b;
    x = 0.0;
  }
  if (B)
  {
    ApplyB(B, r, z, this->use_timer);
  }
  else
  {
    z = r;
  }
  // For a residual v and u = B v, the inner product (u, v) and the residual norm of the
  // stopping test, √|(u, v)| or the Jacobi norm √(D⁻¹ v, v), with one reduction.
  auto DotNorm = [this](const VecType &u, const VecType &v, RealType &norm)
  {
    if (!dinv_norm)
    {
      const auto dot = linalg::Dot(comm, u, v);
      norm = std::sqrt(std::abs(dot));
      return dot;
    }
    auto dot = LocalDotNorm(u, v, *dinv_norm);
    Mpi::GlobalSum(2, dot.data(), comm);
    norm = std::sqrt(std::abs(dot[1]));
    return dot[0];
  };
  beta = DotNorm(z, r, res);
  CheckDot(beta, "PCG preconditioner is not positive definite: (Br, r) = ");
  if (this->initial_guess && dinv_norm)
  {
    DotNorm(b, b, initial_res);
  }
  else if (this->initial_guess)
  {
    ScalarType beta_rhs;
    if (B)
    {
      ApplyB(B, b, p, this->use_timer);
      beta_rhs = linalg::Dot(comm, p, b);
    }
    else
    {
      beta_rhs = linalg::Dot(comm, b, b);
    }
    CheckDot(beta_rhs, "PCG preconditioner is not positive definite: (Bb, b) = ");
    initial_res = std::sqrt(std::abs(beta_rhs));
  }
  else
  {
    initial_res = res;
  }
  eps = std::max(rel_tol * initial_res, abs_tol);
  converged = (res <= eps);

  // Begin iterations. x is updated with each search direction in the same pass as the next
  // one is computed, and with the last one after the iterations.
  int it = 0;
  if (print_opts.iterations)
  {
    Mpi::Print(comm, "{}Residual norms for PCG solve\n",
               std::string(tab_width + int_width - 1, ' '));
  }
  for (; it < max_it && !converged; it++)
  {
    if (print_opts.iterations)
    {
      Mpi::Print(comm, "{}{:{}d} KSP residual norm ||r||_B = {:.6e}\n",
                 std::string(tab_width, ' '), it, int_width, res);
    }
    if (!it)
    {
      p = z;
    }
    else
    {
      linalg::AXPYAndAYPX(alpha, p, x, beta / beta_prev, z);
    }

    {
      perf::ScopedTally t(this->perf_apply);
      A->Mult(p, z);
    }
    denom = linalg::Dot(comm, z, p);
    CheckDot(denom, "PCG operator is not positive definite: (Ap, p) = ");
    MFEM_VERIFY(denom != 0.0,
                "PCG operator is not positive definite: (Ap, p) = " << denom << "!");
    alpha = beta / denom;

    beta_prev = beta;
    if (dinv)
    {
      // The fused loop costs no more than the residual update and inner product alone, so
      // it is not timed as a preconditioner application.
      beta = CgUpdateJacobi(alpha, r, z, *dinv);
      Mpi::GlobalSum(1, &beta, comm);
      res = std::sqrt(std::abs(beta));
    }
    else
    {
      r.Add(-alpha, z);
      if (B)
      {
        perf::ScopedTally t(this->perf_apply);
        ApplyB(B, r, z, this->use_timer);
      }
      else
      {
        z = r;
      }
      beta = DotNorm(z, r, res);
    }
    CheckDot(beta, "PCG preconditioner is not positive definite: (Br, r) = ");
    converged = (res <= eps);
  }
  if (it > 0)
  {
    x.Add(alpha, p);
  }
  if (print_opts.iterations)
  {
    Mpi::Print(comm, "{}{:{}d} KSP residual norm ||r||_B = {:.6e}\n",
               std::string(tab_width, ' '), it, int_width, res);
  }
  if (print_opts.summary || (print_opts.warnings && eps > 0.0 && !converged))
  {
    Mpi::Print(comm, "{}PCG solver {} in {:d} iteration{}", std::string(tab_width, ' '),
               converged ? "converged" : "did NOT converge", it, (it == 1) ? "" : "s");
    if (it > 0)
    {
      Mpi::Print(comm, " (avg. reduction factor: {:.3e})\n",
                 std::pow(res / initial_res, 1.0 / it));
    }
    else
    {
      Mpi::Print(comm, "\n");
    }
  }
  final_res = res;
  final_it = it;
}

template <OperatorType OperType>
void GmresSolver<OperType>::Initialize() const
{
  if (!V.empty())
  {
    MFEM_ASSERT(V.size() == static_cast<std::size_t>(max_dim + 1) &&
                    V[0].Size() == A->Height(),
                "Repeated solves with GmresSolver should not modify the operator size or "
                "restart dimension!");
    return;
  }
  if (max_dim < 0)
  {
    max_dim = max_it;
  }
  constexpr int init_size = 5;
  V.resize(max_dim + 1);
  for (int j = 0; j < std::min(init_size, max_dim + 1); j++)
  {
    V[j].SetSize(A->Height());
    V[j].UseDevice(true);
  }
  // H, s, cs, sn are allocated incrementally in Update() to avoid O(max_dim^2) upfront
  // cost when max_dim is large (e.g., defaulting to max_it).
  s.resize(std::min(init_size + 1, max_dim + 1));
  cs.resize(std::min(init_size + 1, max_dim + 1));
  sn.resize(std::min(init_size + 1, max_dim + 1));
  H.resize(static_cast<std::size_t>(max_dim + 1) * std::min(init_size, max_dim));
  Hraw.resize(H.size());
}

template <OperatorType OperType>
void GmresSolver<OperType>::Update(int j) const
{
  // Add storage for basis vectors, Hessenberg columns, and rotations in increments.
  constexpr int add_size = 10;
  for (int k = j + 1; k < std::min(j + 1 + add_size, max_dim + 1); k++)
  {
    V[k].SetSize(A->Height());
    V[k].UseDevice(true);
  }
  int needed_cols = std::min(j + 1 + add_size, max_dim);
  int current_cols = static_cast<int>(H.size()) / (max_dim + 1);
  if (needed_cols > current_cols)
  {
    H.resize(static_cast<std::size_t>(max_dim + 1) * needed_cols);
    Hraw.resize(H.size());
  }
  auto needed_size = static_cast<std::size_t>(std::min(j + 2 + add_size, max_dim + 1));
  if (needed_size > s.size())
  {
    s.resize(needed_size);
    cs.resize(needed_size);
    sn.resize(needed_size);
  }
}

template <OperatorType OperType>
void GmresSolver<OperType>::RotateColumn(int i) const
{
  ScalarType *Hi = H.data() + static_cast<std::size_t>(i) * (max_dim + 1);
  for (int k = 0; k < i; k++)
  {
    ApplyPlaneRotation(Hi[k], Hi[k + 1], cs[k], sn[k]);
  }
  GeneratePlaneRotation(Hi[i], Hi[i + 1], cs[i], sn[i]);
  ApplyPlaneRotation(Hi[i], Hi[i + 1], cs[i], sn[i]);
  ApplyPlaneRotation(s[i], s[i + 1], cs[i], sn[i]);
}

template <OperatorType OperType>
void GmresSolver<OperType>::DelayedCGS2Step(int j) const
{
  // Below this relative size of the Pythagorean difference (about half of the significant
  // digits lost), compute the norm of the new vector explicitly (one more reduction).
  constexpr RealType pythag_tol = 1.0e-6;
  const std::size_t ld = max_dim + 1;
  ScalarType *Hc = Hraw.data() + j * ld;
  VecType &y = V[j + 1];
  if (j == 0)
  {
    // V[0] is normalized: one Gram-Schmidt pass for y, with the norm of y.
    ScalarType e[2];
    linalg::LocalMultiDot(V, 1, y, nullptr, e, nullptr, e + 1);
    Mpi::GlobalSum(2, e, comm);
    Hc[0] = e[0];
    const RealType d = std::real(e[1]);
    RealType p2 = d - std::norm(e[0]);
    RealType nrm2 = linalg::MultiSubtract(V, 1, y, e, 1.0);
    if (!(p2 > pythag_tol * d))
    {
      Mpi::GlobalSum(1, &nrm2, comm);
      p2 = nrm2;
    }
    Hc[1] = std::sqrt(std::max(p2, RealType(0.0)));
    nu = 1.0;
    return;
  }

  // One pass over the basis and one reduction: a = Qᴴ u, b = Qᴴ y, uᴴ u, uᴴ y, yᴴ y.
  VecType &u = V[j];
  dots.resize(2 * j + 3);
  ScalarType *a = dots.data(), *b = a + j, *e = a + 2 * j;
  linalg::LocalMultiDot(V, j, u, &y, a, b, e);
  Mpi::GlobalSum(2 * j + 3, dots.data(), comm);
  const RealType alpha = std::real(e[0]), d = std::real(e[2]);
  RealType aa = 0.0;
  for (int l = 0; l < j; l++)
  {
    aa += std::norm(a[l]);
  }
  RealType nu2 = alpha - aa;
  const bool extra_pass = !(nu2 >= 0.5 * alpha);
  if (extra_pass)
  {
    // The reorthogonalization removed most of u: orthogonalize it once more, explicitly.
    linalg::MultiSubtract(V, j, u, a, 1.0);
    std::vector<ScalarType> a2(j);
    linalg::LocalMultiDot(V, j, u, nullptr, a2.data(), nullptr);
    Mpi::GlobalSum(j, a2.data(), comm);
    for (int l = 0; l < j; l++)
    {
      a[l] += a2[l];
    }
    nu2 = linalg::MultiSubtract(V, j, u, a2.data(), 1.0);
    Mpi::GlobalSum(1, &nu2, comm);
  }
  nu = std::sqrt(nu2);

  // With q = (u - Q a) / nu, M V[j - 1] = Q (h + a) + nu q completes column j - 1.
  ScalarType *Hp = Hraw.data() + (j - 1) * ld;
  for (int l = 0; l < j; l++)
  {
    Hp[l] += a[l];
  }
  Hp[j] = nu;

  // M q = (y - M Q a) / nu = (y - [Q q] Hbar a) / nu, where Hbar holds the completed
  // columns 0, ..., j - 1, so [Q q]ᴴ M q = (h1 - Hbar a) / nu with h1 = [Qᴴ y; qᴴ y] and
  // qᴴ y = (uᴴ y - aᴴ b) / nu.
  ScalarType ab = 0.0;
  for (int l = 0; l < j; l++)
  {
    ab += Conj(a[l]) * b[l];
  }
  const ScalarType g = (e[1] - ab) / nu;
  RealType hh = std::norm(g);
  for (int r = 0; r <= j; r++)
  {
    ScalarType t = g;
    if (r < j)
    {
      t = b[r];
      hh += std::norm(b[r]);
    }
    for (int l = std::max(r - 1, 0); l < j; l++)
    {
      t -= Hraw[l * ld + r] * a[l];
    }
    Hc[r] = t / nu;
  }

  // One more pass over the basis: V[j] = q and V[j + 1] = (y - Q b - g q) / nu, which is
  // M q - [Q q] Hc after one Gram-Schmidt pass. Its norm follows from the Pythagorean
  // theorem: ‖y - [Q q] h1‖² = yᴴ y - ‖h1‖².
  RealType nrm2 = linalg::MultiSubtract(V, j, u, extra_pass ? nullptr : a, 1.0 / nu, &y, b,
                                        g, 1.0 / nu);
  RealType p2 = (d - hh) / (nu * nu);
  if (!(d - hh > pythag_tol * d))
  {
    Mpi::GlobalSum(1, &nrm2, comm);
    p2 = nrm2;
  }
  Hc[j + 1] = std::sqrt(std::max(p2, RealType(0.0)));
}

template <OperatorType OperType>
void GmresSolver<OperType>::Arnoldi(int j) const
{
  const std::size_t ld = max_dim + 1;
  ScalarType *Hj = H.data() + j * ld;
  if (gs_orthog == Orthogonalization::CGS2)
  {
    if (Hraw.size() < H.size())
    {
      Hraw.resize(H.size());
    }
    {
      perf::ScopedTally t(this->perf_orth);
      DelayedCGS2Step(j);
    }
    if (j > 0)
    {
      // Column j - 1 is complete: rotate it again, from the entries of s before its
      // provisional rotation.
      std::copy_n(Hraw.data() + (j - 1) * ld, j + 1, H.data() + (j - 1) * ld);
      s[j - 1] = s_prev;
      s[j] = 0.0;
      RotateColumn(j - 1);
    }
    std::copy_n(Hraw.data() + j * ld, j + 2, Hj);
    s_prev = s[j];
    RotateColumn(j);
    return;
  }
  {
    perf::ScopedTally t(this->perf_orth);
    VecType &w = V[j + 1];
    if (gs_orthog == Orthogonalization::MGS)
    {
      linalg::OrthogonalizeColumnMGS(comm, V, w, Hj, j + 1);
      Hj[j + 1] = linalg::Norml2(comm, w);
    }
    else
    {
      // One pass over the basis for the projections and one for the update and the norm.
      linalg::LocalMultiDot(V, j + 1, w, nullptr, Hj, nullptr);
      Mpi::GlobalSum(j + 1, Hj, comm);
      RealType nrm2 = linalg::MultiSubtract(V, j + 1, w, Hj, 1.0);
      Mpi::GlobalSum(1, &nrm2, comm);
      Hj[j + 1] = std::sqrt(nrm2);
    }
    w *= 1.0 / Hj[j + 1];
  }
  if (Hraw.size() < H.size())
  {
    Hraw.resize(H.size());
  }
  std::copy_n(Hj, j + 2, Hraw.data() + j * ld);
  RotateColumn(j);
}

template <OperatorType OperType>
void GmresSolver<OperType>::Mult(const VecType &b, VecType &x) const
{
  // Set up workspace.
  RealType beta = 0.0, true_beta, eps = 0.0;
  MFEM_VERIFY(A, "Operator must be set for GmresSolver::Mult!");
  MFEM_ASSERT(A->Width() == x.Size() && A->Height() == b.Size(),
              "Size mismatch for GmresSolver::Mult!");
  r.SetSize(A->Height());
  r.UseDevice(true);
  Initialize();

  // Begin iterations. A left-preconditioned solve whose true residual check fails (see
  // SetTrueResidualCheck) is solved again from zero if it started from a nonzero initial
  // guess (phase 1). If the check still fails, the iterate is continued with right
  // preconditioning (phase 2), on the true residual, to test whether a nearby vector
  // satisfies the check; the iterate is kept if one does, else the solve is falsely
  // converged.
  converged = false;
  this->false_converged = false;
  stagnated = false;
  res_hist.clear();
  int it = 0, restart = 0, phase = 0;
  PreconditionerSide side = pc_side;
  bool nonzero_guess = false;
  RealType res_bound = 0.0, norm_b = 0.0, res_left = 0.0, pres_left = 0.0;
  VecType x_left;
  if (print_opts.iterations)
  {
    Mpi::Print(comm, "{}Residual norms for GMRES solve\n",
               std::string(tab_width + int_width - 1, ' '));
  }
  arnoldi_dim = 0;
  while (true)
  {
    const int restart0 = restart, it_max = (phase == 0) ? max_it : it + max_it;
    const bool from_x = (phase == 0) ? this->initial_guess : (phase == 2);
    for (; it < it_max; restart++)
    {
      // Initialize.
      arnoldi_dim = 0;
      InitialResidual(side, A, B, b, x, r, V[0], (from_x || restart > restart0),
                      this->use_timer);
      true_beta = linalg::Norml2(comm, r);
      CheckDot(true_beta, "GMRES residual norm is not valid: beta = ");
      if (it == 0)
      {
        if (this->initial_guess)
        {
          // True residual of the initial guess (see SetInitialGuessCheck), before V[0],
          // which holds it with left preconditioning, is overwritten.
          const bool check_guess = (guess_tol > 0.0);
          RealType res0 = 0.0, norm_b0 = 0.0;
          if (check_guess)
          {
            res0 = linalg::Norml2(comm, (B && side == PreconditionerSide::LEFT) ? V[0] : r);
            norm_b0 = linalg::Norml2(comm, b);
          }
          RealType beta_rhs;
          if (B && side == PreconditionerSide::LEFT)
          {
            ApplyB(B, b, V[0], this->use_timer);
            beta_rhs = linalg::Norml2(comm, V[0]);
          }
          else  // !B || side == PreconditionerSide::RIGHT
          {
            beta_rhs = linalg::Norml2(comm, b);
          }
          CheckDot(beta_rhs, "GMRES residual norm is not valid: beta_rhs = ");
          initial_res = beta_rhs;

          // The residual of a zero initial guess is the (preconditioned) right-hand side.
          nonzero_guess = (true_beta != beta_rhs);
          if (check_guess && nonzero_guess && norm_b0 > 0.0)
          {
            if (res0 >= guess_tol * norm_b0)
            {
              // Start from zero, whose (preconditioned) residual is already computed.
              x = 0.0;
              if (B && side == PreconditionerSide::LEFT)
              {
                r = V[0];
              }
              else
              {
                r = b;
              }
              true_beta = beta_rhs;
              nonzero_guess = false;
              if (print_opts.summary || print_opts.warnings)
              {
                Mpi::Print(comm,
                           "{}GMRES initial guess discarded: its true residual is {:.3e} "
                           "norm(b), solving from zero\n",
                           std::string(tab_width, ' '), res0 / norm_b0);
              }
            }
          }
        }
        else
        {
          initial_res = true_beta;
        }
        eps = std::max(rel_tol * initial_res, abs_tol);
      }
      else if (restart > restart0 && beta > 0.0 &&
               std::abs(beta - true_beta) > 0.1 * true_beta && print_opts.warnings)
      {
        Mpi::Print(
            comm,
            "{}GMRES residual at restart ({:.6e}) is far from the residual norm estimate "
            "from the recursion formula ({:.6e}) (initial residual = {:.6e})\n",
            std::string(tab_width, ' '), true_beta, beta, initial_res);
      }
      beta = true_beta;
      if (beta <= eps)
      {
        converged = true;
        break;
      }
      if (it == 0 || restart == restart0)
      {
        res_hist.push_back(beta);
      }

      V[0] = 0.0;
      V[0].Add(1.0 / beta, r);
      std::ranges::fill(s, 0.0);
      s[0] = beta;

      int j = 0;
      for (;; j++, it++)
      {
        if (print_opts.iterations)
        {
          Mpi::Print(comm, "{}{:{}d} (restart {:d}) KSP residual norm {:.6e}\n",
                     std::string(tab_width, ' '), it, int_width, restart, beta);
        }
        VecType &w = V[j + 1];
        if (w.Size() == 0)
        {
          Update(j);
        }
        {
          perf::ScopedTally t(this->perf_apply);
          ApplyBA(side, A, B, V[j], w, r, this->use_timer);
        }

        Arnoldi(j);

        beta = std::abs(s[j + 1]);
        CheckDot(beta, "GMRES residual norm is not valid: beta = ");
        converged = (beta <= eps);
        if (converged || Stagnated(beta) || j + 1 == max_dim || it + 1 == it_max ||
            (phase == 2 && !std::isfinite(beta)))
        {
          it++;
          break;
        }
      }

      // Columns 0, ..., j - 1 of the Hessenberg matrix and the basis vectors V[0], ...,
      // V[j] are final (with CGS2, column j and V[j + 1] miss their reorthogonalization).
      // The right-preconditioned continuation of a left-preconditioned solve has no Arnoldi
      // relation for B A (see GetArnoldiRelation).
      arnoldi_dim = (side == pc_side) ? j : 0;

      // Reconstruct the solution (for restart or due to convergence or maximum iterations).
      for (int i = j; i >= 0; i--)
      {
        ScalarType *Hi = H.data() + i * (max_dim + 1);
        s[i] /= Hi[i];
        for (int k = i - 1; k >= 0; k--)
        {
          s[k] -= Hi[k] * s[i];
        }
      }
      if (!B || side == PreconditionerSide::LEFT)
      {
        for (int k = 0; k <= j; k++)
        {
          x.Add(s[k], V[k]);
        }
      }
      else  // B && side == PreconditionerSide::RIGHT
      {
        r = 0.0;
        for (int k = 0; k <= j; k++)
        {
          r.Add(s[k], V[k]);
        }
        ApplyB(B, r, V[0], this->use_timer);
        x += V[0];
      }
      if (converged || stagnated || (phase == 2 && !std::isfinite(beta)))
      {
        break;
      }
    }
    if (phase == 0 && !(converged && it > 0 && B && side == PreconditionerSide::LEFT &&
                        true_res_tol > 0.0 && rel_tol * initial_res >= abs_tol))
    {
      break;
    }

    // Check the true residual (see SetTrueResidualCheck). A solve from zero (phase 1) which
    // does not converge leaves the converged iterate from the initial guess to continue.
    RealType res = res_left;
    if (phase == 1 && !converged)
    {
      x = x_left;
    }
    else
    {
      A->Mult(x, r);
      linalg::AXPBY(1.0, b, -1.0, r);
      res = linalg::Norml2(comm, r);
    }
    if (phase == 0)
    {
      norm_b = linalg::Norml2(comm, b);
      res_bound = std::max(true_res_tol, 1.0e7 * rel_tol) * norm_b;
    }
    if (phase < 2)
    {
      if (res <= res_bound)
      {
        break;
      }
      if (phase == 0 || converged)
      {
        // Keep the left-preconditioned iterate and its explicit preconditioned residual.
        if (x_left.Size() != x.Size())
        {
          x_left.SetSize(x.Size());
          x_left.UseDevice(true);
        }
        x_left = x;
        res_left = res;
        ApplyB(B, r, V[0], this->use_timer);
        pres_left = linalg::Norml2(comm, V[0]);
      }
      if (std::isfinite(res))
      {
        // A nonzero initial guess can leave residual components which the preconditioned
        // residual does not see: solve again from zero first.
        const bool resolve = (phase == 0 && nonzero_guess);
        Mpi::Warning(comm,
                     "GMRES converged on the preconditioned residual ({:.3e} relative) but "
                     "the true residual norm(b-Ax)/norm(b) is {:.3e}: {}!\n",
                     pres_left / initial_res, res / norm_b,
                     resolve
                         ? "solving again from zero"
                         : "continuing on the true residual with right preconditioning");
        phase = resolve ? 1 : 2;
        if (phase == 2)
        {
          side = PreconditionerSide::RIGHT;
          eps = 0.5 * res_bound;
        }
        beta = 0.0;
        converged = stagnated = false;
        res_hist.clear();
        restart++;
        continue;
      }
    }
    else
    {
      // The continuation must satisfy the check close to the iterate, which is kept (the
      // continuation only minimizes the true residual and can drift along the resonant
      // components which the preconditioned residual resolves).
      V[0] = x;
      V[0] -= x_left;
      const RealType dx = linalg::Norml2(comm, V[0]), norm_x = linalg::Norml2(comm, x_left);
      x = x_left;
      if (res <= res_bound && dx <= 1.0e-4 * norm_x)
      {
        Mpi::Warning(comm,
                     "GMRES with right preconditioning reduced the true residual "
                     "norm(b-Ax)/norm(b) to {:.3e} within {:.1e} of the solution, which is "
                     "kept (true residual {:.3e}, preconditioned residual {:.3e} relative, "
                     "{:d} iterations in total)!\n",
                     res / norm_b, dx / norm_x, res_left / norm_b, pres_left / initial_res,
                     it);
        converged = true;
        beta = pres_left;
        break;
      }
      Mpi::Warning(comm,
                   "GMRES with right preconditioning reached a true residual "
                   "norm(b-Ax)/norm(b) of {:.3e} at {:.1e} of the solution!\n",
                   res / norm_b, dx / norm_x);
    }
    converged = false;
    this->false_converged = true;
    Mpi::Warning(comm,
                 "GMRES converged on the preconditioned residual ({:.3e} relative) but the "
                 "true residual norm(b-Ax)/norm(b) is {:.3e}: the preconditioner is "
                 "unstable (false convergence)!\n",
                 pres_left / initial_res, res_left / norm_b);
    beta = pres_left;
    break;
  }
  if (print_opts.iterations)
  {
    Mpi::Print(comm, "{}{:{}d} (restart {:d}) KSP residual norm {:.6e}\n",
               std::string(tab_width, ' '), it, int_width, restart, beta);
  }
  if (print_opts.summary || (print_opts.warnings && eps > 0.0 && !converged))
  {
    Mpi::Print(comm, "{}GMRES solver {} in {:d} iteration{}", std::string(tab_width, ' '),
               converged ? "converged" : "did NOT converge", it, (it == 1) ? "" : "s");
    if (it > 0)
    {
      Mpi::Print(comm, " (avg. reduction factor: {:.3e}){}\n",
                 std::pow(beta / initial_res, 1.0 / it),
                 stagnated ? fmt::format(", stagnated: residual reduced by less than "
                                         "{:.1e} over the last {:d} iterations",
                                         stag_tol, stag_it)
                           : "");
    }
    else
    {
      Mpi::Print(comm, "\n");
    }
  }
  final_res = beta;
  final_it = it;
  final_restarts = converged ? restart : std::max(restart - 1, 0);
}

template <OperatorType OperType>
void FgmresSolver<OperType>::Initialize() const
{
  GmresSolver<OperType>::Initialize();
  if (gs_orthog == Orthogonalization::CGS2)
  {
    T.resize(H.size());
  }
  constexpr int init_size = 5;
  Z.resize(max_dim + 1);
  for (int j = 0; j < std::min(init_size, max_dim + 1); j++)
  {
    Z[j].SetSize(A->Height());
    Z[j].UseDevice(true);
  }
}

template <OperatorType OperType>
void FgmresSolver<OperType>::Update(int j) const
{
  // Add storage for basis vectors in increments.
  GmresSolver<OperType>::Update(j);
  if (gs_orthog == Orthogonalization::CGS2)
  {
    T.resize(H.size());
  }
  constexpr int add_size = 10;
  for (int k = j + 1; k < std::min(j + 1 + add_size, max_dim + 1); k++)
  {
    Z[k].SetSize(A->Height());
    Z[k].UseDevice(true);
  }
}

template <OperatorType OperType>
void FgmresSolver<OperType>::Mult(const VecType &b, VecType &x) const
{
  // Set up workspace.
  RealType beta = 0.0, true_beta, eps = 0.0;
  MFEM_VERIFY(A && B, "Operator and preconditioner must be set for FgmresSolver::Mult!");
  MFEM_ASSERT(A->Width() == x.Size() && A->Height() == b.Size(),
              "Size mismatch for FgmresSolver::Mult!");
  Initialize();

  // Begin iterations.
  converged = false;
  this->stagnated = false;
  this->res_hist.clear();
  int it = 0, restart = 0;
  if (print_opts.iterations)
  {
    Mpi::Print(comm, "{}Residual norms for FGMRES solve\n",
               std::string(tab_width + int_width - 1, ' '));
  }
  for (; it < max_it; restart++)
  {
    // Initialize.
    InitialResidual(PreconditionerSide::RIGHT, A, B, b, x, Z[0], V[0],
                    (this->initial_guess || restart > 0), this->use_timer);
    true_beta = linalg::Norml2(comm, Z[0]);
    CheckDot(true_beta, "FGMRES residual norm is not valid: beta = ");
    if (it == 0)
    {
      if (this->initial_guess)
      {
        auto beta_rhs = linalg::Norml2(comm, b);
        CheckDot(beta_rhs, "FGMRES residual norm is not valid: beta_rhs = ");
        initial_res = beta_rhs;
      }
      else
      {
        initial_res = true_beta;
      }
      eps = std::max(rel_tol * initial_res, abs_tol);
    }
    else if (beta > 0.0 && std::abs(beta - true_beta) > 0.1 * true_beta &&
             print_opts.warnings)
    {
      Mpi::Print(
          comm,
          "{}FGMRES residual at restart ({:.6e}) is far from the residual norm estimate "
          "from the recursion formula ({:.6e}) (initial residual = {:.6e})\n",
          std::string(tab_width, ' '), true_beta, beta, initial_res);
    }
    beta = true_beta;
    if (beta <= eps)
    {
      converged = true;
      break;
    }
    if (it == 0)
    {
      this->res_hist.push_back(beta);
    }

    V[0] = 0.0;
    V[0].Add(1.0 / beta, Z[0]);
    std::ranges::fill(s, 0.0);
    s[0] = beta;

    int j = 0;
    for (;; j++, it++)
    {
      if (print_opts.iterations)
      {
        Mpi::Print(comm, "{}{:{}d} (restart {:d}) KSP residual norm {:.6e}\n",
                   std::string(tab_width, ' '), it, int_width, restart, beta);
      }
      VecType &w = V[j + 1];
      if (w.Size() == 0)
      {
        Update(j);
      }
      {
        perf::ScopedTally t(this->perf_apply);
        ApplyBA(PreconditionerSide::RIGHT, A, B, V[j], w, Z[j], this->use_timer);
      }

      this->Arnoldi(j);
      if (gs_orthog == Orthogonalization::CGS2)
      {
        // Z[j] = B u for V[j] = u before its reorthogonalization, so the preconditioned
        // vector of the final V[j] = (u - V a) / nu is (Z[j] - Z T a) / nu (column j of Z
        // T).
        const std::size_t ld = max_dim + 1;
        if (T.size() < H.size())
        {
          T.resize(H.size());
        }
        ScalarType *Tj = T.data() + j * ld;
        std::fill_n(Tj, j + 1, ScalarType(0.0));
        Tj[j] = 1.0;
        for (int l = 0; l < j; l++)
        {
          const ScalarType *Tl = T.data() + l * ld;
          for (int i = 0; i <= l; i++)
          {
            Tj[i] -= Tl[i] * dots[l];
          }
        }
        for (int i = 0; i <= j; i++)
        {
          Tj[i] /= nu;
        }
      }

      beta = std::abs(s[j + 1]);
      CheckDot(beta, "FGMRES residual norm is not valid: beta = ");
      converged = (beta <= eps);
      if (converged || this->Stagnated(beta) || j + 1 == max_dim || it + 1 == max_it)
      {
        it++;
        break;
      }
    }

    // Reconstruct the solution (for restart or due to convergence or maximum iterations).
    for (int i = j; i >= 0; i--)
    {
      ScalarType *Hi = H.data() + i * (max_dim + 1);
      s[i] /= Hi[i];
      for (int k = i - 1; k >= 0; k--)
      {
        s[k] -= Hi[k] * s[i];
      }
    }
    if (gs_orthog == Orthogonalization::CGS2)
    {
      const std::size_t ld = max_dim + 1;
      for (int i = 0; i <= j; i++)
      {
        ScalarType c = 0.0;
        for (int l = i; l <= j; l++)
        {
          c += T[l * ld + i] * s[l];
        }
        x.Add(c, Z[i]);
      }
    }
    else
    {
      for (int k = 0; k <= j; k++)
      {
        x.Add(s[k], Z[k]);
      }
    }
    if (converged || this->stagnated)
    {
      break;
    }
  }
  if (print_opts.iterations)
  {
    Mpi::Print(comm, "{}{:{}d} (restart {:d}) KSP residual norm {:.6e}\n",
               std::string(tab_width, ' '), it, int_width, restart, beta);
  }
  if (print_opts.summary || (print_opts.warnings && eps > 0.0 && !converged))
  {
    Mpi::Print(comm, "{}FGMRES solver {} in {:d} iteration{}", std::string(tab_width, ' '),
               converged ? "converged" : "did NOT converge", it, (it == 1) ? "" : "s");
    if (it > 0)
    {
      Mpi::Print(comm, " (avg. reduction factor: {:.3e}){}\n",
                 std::pow(beta / initial_res, 1.0 / it),
                 this->stagnated ? fmt::format(", stagnated: residual reduced by less than "
                                               "{:.1e} over the last {:d} iterations",
                                               this->stag_tol, this->stag_it)
                                 : "");
    }
    else
    {
      Mpi::Print(comm, "\n");
    }
  }
  final_res = beta;
  final_it = it;
  final_restarts = converged ? restart : std::max(restart - 1, 0);
}

template class IterativeSolver<Operator>;
template class IterativeSolver<ComplexOperator>;
template class CgSolver<Operator>;
template class CgSolver<ComplexOperator>;
template class GmresSolver<Operator>;
template class GmresSolver<ComplexOperator>;
template class FgmresSolver<Operator>;
template class FgmresSolver<ComplexOperator>;

}  // namespace palace
