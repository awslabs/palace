// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "iterative.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <limits>
#include <span>
#include <string>
#include <typeinfo>
#include "linalg/jacobi.hpp"
#include "linalg/orthog.hpp"
#include "utils/communication.hpp"
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

// Replaces the nonzero initial guess x by α x, with α = (u, w) / (u, u) minimizing the
// initial residual r = w - α u in the norm of the stopping test (w = B b and u = B A x for
// left preconditioning, else w = b and u = A x), sets r, and returns ||w||, the residual
// norm of a zero initial guess. Then ||r||² = ||w||² - |(u, w)|² / ||u||² ≤ ||w||², and an
// exact guess stays exact (α = 1); products which are not finite give α = 0. With left
// preconditioning, u = w - r₁ reuses the residual r₁ = B (b - A x), so B is applied as for
// the unscaled guess, and r = (1 - α) w + α r₁ carries the rounding error of B relative to
// about |1 - α| ||w|| + |α| ||r₁|| (relative to ||w|| with u = B A x). That is amplified
// for |α| > 1 with an inexact preconditioner, so for |α| > 2, the residual B (b - α A x) is
// computed instead (one more application of A and B).
template <typename OperType, typename VecType>
inline double ScaledInitialResidual(MPI_Comm comm, PreconditionerSide side,
                                    const OperType *A, const Solver<OperType> *B,
                                    const VecType &b, VecType &x, VecType &r, VecType &z,
                                    bool use_timer = true)
{
  using ScalarType = decltype(linalg::LocalDot(r, r));
  const bool left = (B && side == PreconditionerSide::LEFT);
  if (left)
  {
    InitialResidual(side, A, B, b, x, r, z, true, use_timer);  // r = r₁
    ApplyB(B, b, z, use_timer);
    linalg::AXPBY(1.0, z, -1.0, r);  // r = u
  }
  else
  {
    A->Mult(x, r);  // r = u
  }
  const VecType &w = left ? z : b, &u = r;

  // The inner products (w, w), (u, w) = uᴴ w, and (u, u) in one reduction.
  ScalarType dot[3] = {linalg::LocalDot(w, w), linalg::LocalDot(w, u),
                       linalg::LocalDot(u, u)};
  Mpi::GlobalSum(3, dot, comm);
  const double uu = std::real(dot[2]);
  ScalarType alpha = (uu > 0.0) ? ScalarType(dot[1] / uu) : ScalarType(0.0);
  if (!std::isfinite(std::abs(alpha)))
  {
    alpha = 0.0;
  }
  if (alpha == 0.0)
  {
    x = 0.0;
    r = w;
  }
  else if (left && std::abs(alpha) > 2.0)
  {
    x *= alpha;
    InitialResidual(side, A, B, b, x, r, z, true, use_timer);
  }
  else
  {
    x *= alpha;
    linalg::AXPBY(ScalarType(1.0), w, -alpha, r);
  }
  return std::sqrt(std::abs(dot[0]));
}

// Sets the residual r of x at the start of a GMRES cycle and returns its norm. The first
// cycle starts from zero without an initial guess, else from the initial guess unless it is
// zero, scaled as above with scale_guess, and sets initial_res, the reference of the
// stopping test (the residual norm of a zero initial guess).
template <typename OperType, typename VecType>
inline double CycleResidual(MPI_Comm comm, PreconditionerSide side, const OperType *A,
                            const Solver<OperType> *B, const VecType &b, VecType &x,
                            VecType &r, VecType &z, bool first_cycle, bool initial_guess,
                            bool scale_guess, double &initial_res, bool use_timer = true)
{
  if (first_cycle && initial_guess && linalg::Norml2(comm, x) != 0.0)
  {
    if (scale_guess)
    {
      initial_res = ScaledInitialResidual(comm, side, A, B, b, x, r, z, use_timer);
    }
    else
    {
      InitialResidual(side, A, B, b, x, r, z, true, use_timer);
      const bool left = (B && side == PreconditionerSide::LEFT);
      if (left)
      {
        ApplyB(B, b, z, use_timer);
      }
      initial_res = linalg::Norml2(comm, left ? z : b);
    }
    return linalg::Norml2(comm, r);
  }
  InitialResidual(side, A, B, b, x, r, z, !first_cycle, use_timer);
  const double beta = linalg::Norml2(comm, r);
  if (first_cycle)
  {
    initial_res = beta;
  }
  return beta;
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
  : Solver<OperType>(), comm(comm), A(nullptr), B(nullptr)
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
  // inner product (not a subclass, which may change the application).
  const Vector *dinv = nullptr;
  if (B && typeid(*B) == typeid(JacobiSmoother<OperType>))
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
  beta = linalg::Dot(comm, z, r);
  CheckDot(beta, "PCG preconditioner is not positive definite: (Br, r) = ");
  res = std::sqrt(std::abs(beta));
  if (this->initial_guess)
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

  // Cumulative decrease of the least-squares objective, φ(x₀) - φ(xₖ), for the
  // least-squares convergence test (with x₀ = 0, φ(x₀) = ‖d‖²).
  const bool ls_test = (ls_tol >= 0.0);
  if (ls_test)
  {
    MFEM_VERIFY(!this->initial_guess,
                "The least-squares convergence test of PCG requires a zero initial guess!");
    ls_decrease.assign(1, 0.0);
  }

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

    A->Mult(p, z);
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
    }
    else
    {
      r.Add(-alpha, z);
      if (B)
      {
        ApplyB(B, r, z, this->use_timer);
      }
      else
      {
        z = r;
      }
      beta = linalg::Dot(comm, z, r);
    }
    CheckDot(beta, "PCG preconditioner is not positive definite: (Br, r) = ");
    res = std::sqrt(std::abs(beta));
    converged = (res <= eps);
    if (ls_test)
    {
      // The objective decreases by αₖ (rₖ, zₖ) = |(rₖ, zₖ)|² / (pₖ, A pₖ) in iteration k.
      // After m iterations, the squared error of xₘ is estimated by the decrease over the
      // iterations ⌈m/2⌉, ..., m - 1, extrapolated when it is more than half of that over
      // ⌈m/4⌉, ..., ⌈m/2⌉ - 1 (see SetLeastSquaresTol), and compared to tol² φ(xₘ). The
      // iterations also stop once the decrease is below the rounding error of φ.
      const int m = it + 1;
      ls_decrease.push_back(ls_decrease.back() +
                            std::abs(beta_prev) * std::abs(beta_prev) / std::abs(denom));
      const RealType phi = ls_norm2 - ls_decrease[m];
      const RealType last = ls_decrease[m] - ls_decrease[(m + 1) / 2];
      const RealType prev = ls_decrease[(m + 1) / 2] - ls_decrease[(m + 3) / 4];
      RealType err2 = std::numeric_limits<RealType>::infinity();
      if (last <= 0.5 * prev)
      {
        err2 = last;
      }
      else if (last < prev)
      {
        err2 = last * last / (prev - last);
      }
      converged = converged ||
                  (m >= 2 && (err2 <= ls_tol * ls_tol * phi ||
                              last <= std::numeric_limits<RealType>::epsilon() * ls_norm2));
    }
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
typename GmresSolver<OperType>::RealType GmresSolver<OperType>::CorrectionNorm(int j) const
{
  // Back substitution with the rotated Hessenberg matrix (upper triangular), on a copy of
  // the rotated right-hand side: O(j²) scalar operations, no vector operations.
  std::vector<ScalarType> y(s.begin(), s.begin() + j + 1);
  for (int i = j; i >= 0; i--)
  {
    const ScalarType *Hi = H.data() + i * (max_dim + 1);
    y[i] /= Hi[i];
    for (int k = i - 1; k >= 0; k--)
    {
      y[k] -= Hi[k] * y[i];
    }
  }
  RealType norm2 = 0.0;
  for (const auto &yi : y)
  {
    norm2 += std::norm(yi);
  }
  return std::sqrt(norm2);
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

  // Monitor the correction to a nonzero initial guess, if requested: x - x₀ = V y with
  // orthonormal V for left or no preconditioning (first cycle only). With a correction
  // limit, the guess is not scaled, so that the correction is to x₀ as given.
  correction_exceeded = false;
  RealType x0_norm = 0.0;
  const bool monitor_correction = max_rel_correction > 0.0 && this->initial_guess &&
                                  (!B || pc_side == PreconditionerSide::LEFT) &&
                                  (x0_norm = linalg::Norml2(comm, x)) > 0.0;

  // Begin iterations.
  converged = false;
  int it = 0, restart = 0;
  if (print_opts.iterations)
  {
    Mpi::Print(comm, "{}Residual norms for GMRES solve\n",
               std::string(tab_width + int_width - 1, ' '));
  }
  for (; it < max_it; restart++)
  {
    // Initialize.
    true_beta =
        CycleResidual(comm, pc_side, A, B, b, x, r, V[0], it == 0, this->initial_guess,
                      max_rel_correction <= 0.0, initial_res, this->use_timer);
    CheckDot(true_beta, "GMRES residual norm is not valid: beta = ");
    if (it == 0)
    {
      CheckDot(initial_res, "GMRES residual norm is not valid: initial_res = ");
      eps = std::max(rel_tol * initial_res, abs_tol);
    }
    else if (beta > 0.0 && std::abs(beta - true_beta) > 0.1 * true_beta &&
             print_opts.warnings)
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
      ApplyBA(pc_side, A, B, V[j], w, r, this->use_timer);

      ScalarType *Hj = H.data() + j * (max_dim + 1);
      Hj[j + 1] =
          linalg::OrthogonalizeColumn(gs_orthog, comm, std::span(V).first(j + 2), Hj);
      w *= 1.0 / Hj[j + 1];

      for (int k = 0; k < j; k++)
      {
        ApplyPlaneRotation(Hj[k], Hj[k + 1], cs[k], sn[k]);
      }
      GeneratePlaneRotation(Hj[j], Hj[j + 1], cs[j], sn[j]);
      ApplyPlaneRotation(Hj[j], Hj[j + 1], cs[j], sn[j]);
      ApplyPlaneRotation(s[j], s[j + 1], cs[j], sn[j]);

      beta = std::abs(s[j + 1]);
      CheckDot(beta, "GMRES residual norm is not valid: beta = ");
      converged = (beta <= eps);
      if (monitor_correction && restart == 0)
      {
        const RealType corr = CorrectionNorm(j);
        correction_exceeded = (corr > max_rel_correction * (x0_norm + corr));
      }
      if (converged || correction_exceeded || j + 1 == max_dim || it + 1 == max_it)
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
    if (!B || pc_side == PreconditionerSide::LEFT)
    {
      for (int k = 0; k <= j; k++)
      {
        x.Add(s[k], V[k]);
      }
    }
    else  // B && pc_side == PreconditionerSide::RIGHT
    {
      r = 0.0;
      for (int k = 0; k <= j; k++)
      {
        r.Add(s[k], V[k]);
      }
      ApplyB(B, r, V[0], this->use_timer);
      x += V[0];
    }
    if (converged || correction_exceeded)
    {
      break;
    }
  }
  if (print_opts.iterations)
  {
    Mpi::Print(comm, "{}{:{}d} (restart {:d}) KSP residual norm {:.6e}\n",
               std::string(tab_width, ' '), it, int_width, restart, beta);
  }
  if (correction_exceeded)
  {
    if (print_opts.summary || print_opts.warnings)
    {
      Mpi::Print(
          comm,
          "{}GMRES solver stopped in {:d} iteration{}: the correction to the initial "
          "guess exceeds {:.3e} of the solution\n",
          std::string(tab_width, ' '), it, (it == 1) ? "" : "s", max_rel_correction);
    }
  }
  else if (print_opts.summary || (print_opts.warnings && eps > 0.0 && !converged))
  {
    Mpi::Print(comm, "{}GMRES solver {} in {:d} iteration{}", std::string(tab_width, ' '),
               converged ? "converged" : "did NOT converge", it, (it == 1) ? "" : "s");
    if (it > 0)
    {
      Mpi::Print(comm, " (avg. reduction factor: {:.3e})\n",
                 std::pow(beta / initial_res, 1.0 / it));
    }
    else
    {
      Mpi::Print(comm, "\n");
    }
  }
  final_res = beta;
  final_it = it;
}

template <OperatorType OperType>
void FgmresSolver<OperType>::Initialize() const
{
  GmresSolver<OperType>::Initialize();
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
  int it = 0, restart = 0;
  if (print_opts.iterations)
  {
    Mpi::Print(comm, "{}Residual norms for FGMRES solve\n",
               std::string(tab_width + int_width - 1, ' '));
  }
  for (; it < max_it; restart++)
  {
    // Initialize.
    true_beta = CycleResidual(comm, PreconditionerSide::RIGHT, A, B, b, x, Z[0], V[0],
                              it == 0, this->initial_guess, max_rel_correction <= 0.0,
                              initial_res, this->use_timer);
    CheckDot(true_beta, "FGMRES residual norm is not valid: beta = ");
    if (it == 0)
    {
      CheckDot(initial_res, "FGMRES residual norm is not valid: initial_res = ");
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
      ApplyBA(PreconditionerSide::RIGHT, A, B, V[j], w, Z[j], this->use_timer);

      ScalarType *Hj = H.data() + j * (max_dim + 1);
      Hj[j + 1] =
          linalg::OrthogonalizeColumn(gs_orthog, comm, std::span(V).first(j + 2), Hj);
      w *= 1.0 / Hj[j + 1];

      for (int k = 0; k < j; k++)
      {
        ApplyPlaneRotation(Hj[k], Hj[k + 1], cs[k], sn[k]);
      }
      GeneratePlaneRotation(Hj[j], Hj[j + 1], cs[j], sn[j]);
      ApplyPlaneRotation(Hj[j], Hj[j + 1], cs[j], sn[j]);
      ApplyPlaneRotation(s[j], s[j + 1], cs[j], sn[j]);

      beta = std::abs(s[j + 1]);
      CheckDot(beta, "FGMRES residual norm is not valid: beta = ");
      converged = (beta <= eps);
      if (converged || j + 1 == max_dim || it + 1 == max_it)
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
    for (int k = 0; k <= j; k++)
    {
      x.Add(s[k], Z[k]);
    }
    if (converged)
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
      Mpi::Print(comm, " (avg. reduction factor: {:.3e})\n",
                 std::pow(beta / initial_res, 1.0 / it));
    }
    else
    {
      Mpi::Print(comm, "\n");
    }
  }
  final_res = beta;
  final_it = it;
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
