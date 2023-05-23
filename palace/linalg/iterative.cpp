// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "iterative.hpp"

#include <cmath>
#include <algorithm>
#include <limits>
#include <string>
#include "linalg/orthog.hpp"
#include "utils/communication.hpp"

namespace palace
{

namespace
{

template <typename T>
inline void CheckDot(T dot, std::string msg)
{
  MFEM_ASSERT(std::isfinite(dot) && dot >= 0.0, msg);
}

template <typename T>
inline void CheckDot(std::complex<T> dot, std::string msg)
{
  MFEM_ASSERT(std::isfinite(dot.real()) && std::is_finite(dot.imag()) && dot.real() >= 0.0,
              msg);
}

template <typename T>
inline constexpr T SafeMin()
{
  // Originally part of <T>LAPACK.
  // <T>LAPACK is free software: you can redistribute it and/or modify it under
  // the terms of the BSD 3-Clause license.
  //
  // Copyright (c) 2021-2023, University of Colorado Denver. All rights reserved.
  // Copyright (c) 2017-2021, University of Tennessee. All rights reserved.
  //
  // Original author: Weslley S Pereira, University of Colorado Denver, USA
  constexpr int fradix = std::numeric_limits<T>::radix;
  constexpr int expm = std::numeric_limits<T>::min_exponent;
  constexpr int expM = std::numeric_limits<T>::max_exponent;
  return std::max(std::pow(fradix, T(expm - 1)), std::pow(fradix, T(1 - expM)));
}

template <typename T>
inline constexpr T SafeMax()
{
  // Originally part of <T>LAPACK.
  // <T>LAPACK is free software: you can redistribute it and/or modify it under
  // the terms of the BSD 3-Clause license.
  //
  // Copyright (c) 2021-2023, University of Colorado Denver. All rights reserved.
  // Copyright (c) 2017-2021, University of Tennessee. All rights reserved.
  //
  // Original author: Weslley S Pereira, University of Colorado Denver, USA
  constexpr int fradix = std::numeric_limits<T>::radix;
  constexpr int expm = std::numeric_limits<T>::min_exponent;
  constexpr int expM = std::numeric_limits<T>::max_exponent;
  return std::min(std::pow(fradix, T(1 - expm)), std::pow(fradix, T(expM - 1)));
}

template <typename T>
inline void GeneratePlaneRotation(const T dx, const T dy, T &cs, T &sn)
{
  // See LAPACK's s/dlartg.
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
  const T root_min = std::sqrt(SafeMin());
  const T root_max = std::sqrt(SafeMax() / 2);
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
    T u = std::min(SafeMax(), std::max(SafeMin(), std::max(dx1, dy1)));
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
  //   [  cs        sn ] . [ dx ]  =  [ r ]
  //   [ -conj(sn)  cs ]   [ dy ]     [ 0 ]
  // where cs is real and cs² + |sn|² = 1. See LAPACK's c/zlartg.
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
      const T root_min = std::sqrt(SafeMin());
      const T root_max = std::sqrt(SafeMax() / 2);
      T dy1 = std::max(std::abs(dy.real()), std::abs(dy.imag()));
      if (dy1 > root_min && dy1 < root_max)
      {
        sn = std::conj(dy) / std::sqrt(dy.real() * dy.real() + dy.imag() * dy.imag());
      }
      else
      {
        T u = std::min(SafeMax(), std::max(SafeMin(), dy1));
        std::complex<T> dys = dy / u;
        sn = std::conj(dys) / std::sqrt(dys.real() * dys.real() + dys.imag() * dys.imag());
      }
    }
    return;
  }
  const T root_min = std::sqrt(SafeMin());
  const T root_max = std::sqrt(SafeMax() / 4);
  T dx1 = std::max(std::abs(dx.real()), std::abs(dx.imag()));
  T dy1 = std::max(std::abs(dy.real()), std::abs(dy.imag()));
  if (dx1 > root_min && dx1 < root_max && dy1 > root_min && dy1 < root_max)
  {
    T dx2 = dx.real() * dx.real() + dx.imag() * dx.imag();
    T dy2 = dy.real() * dy.real() + dy.imag() * dy.imag();
    T dz2 = dx2 + dy2;
    if (dx2 >= dz2 * SafeMin())
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
    T u = std::min(SafeMax(), std::max(SafeMin(), std::max(dx1, dy1))), w;
    std::complex<T> dys = dy / u, dxs;
    T dy2 = dys.real() * dys.real() + dys.imag() * dys.imag(), dx2, dz2;
    if (dx1 / u < root_min)
    {
      T v = std::min(SafeMax(), std::max(SafeMin(), dx1));
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
    if (dx2 >= dz2 * SafeMin())
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

}  // namespace

template <typename OperType>
IterativeSolver<OperType>::IterativeSolver(MPI_comm comm, int print)
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
  initial_res = final_res = 0.0;
  final_it = 0;
}

template <typename OperType>
void CgSolver<OperType>::Mult(const VecType &b, VecType &x) const
{
  // Set up workspace.
  ScalarType beta, beta_prev, alpha, denom;
  RealType res, eps;
  MFEM_VERIFY(A, "Operator must be set for CgSolver::Mult!");
  MFEM_ASSERT(A->Width() == x.Size() && A->Height() == b.Size(),
              "Size mismatch for CgSolver::Mult!");
  r.SetSize(A->Height());
  z.SetSize(A->Height());
  p.SetSize(A->Height());

  // Initialize.
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
  if (B)
  {
    B->Mult(r, z);
  }
  else
  {
    z = r;
  }
  beta = linalg::Dot(comm, z, r);
  CheckDot(beta, "PCG preconditioner is not positive definite: (Br, r) = " << beta << "!");
  res = initial_res = std::sqrt(std::abs(beta));
  eps = std::max(rel_tol * res, abs_tol);
  converged = (res < eps);

  // Begin iterations.
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
      Mpi::Print(comm, "{}{:{}d} iteration, residual (B r, r) = {:.6e}\n",
                 std::string(tab_width, ' '), it, int_width, beta);
    }
    if (!it)
    {
      p = z;
    }
    else
    {
      linalg::AXPBY(1.0, z, beta / beta_prev, p);
    }

    A->Mult(p, z);
    denom = linalg::Dot(comm, z, p);
    CheckDot(denom, "PCG operator is not positive definite: (Ap, p) = " << denom << "!");
    alpha = beta / denom;

    x.Add(alpha, p);
    r.Add(-alpha, z);

    beta_prev = beta;
    if (B)
    {
      B->Mult(r, z);
    }
    else
    {
      z = r;
    }
    beta = linalg::Dot(comm, z, r);
    CheckDot(beta,
             "PCG preconditioner is not positive definite: (Br, r) = " << beta << "!");
    res = std::sqrt(std::abs(beta));
    converged = (res < eps);
  }
  if (print_opts.iterations)
  {
    Mpi::Print(comm, "{}{:{}d} iteration, residual (B r, r) = {:.6e}\n",
               std::string(tab_width, ' '), it, int_width, beta);
  }
  if (print_opts.summary || (print_opts.warnings && !converged))
  {
    Mpi::Print(comm, "{}PCG solver {} with {:d} iteration{}", std::string(tab_width, ' '),
               converged ? "converged" : "did NOT converge", it, (it == 1) ? "" : "s");
    if (it > 0)
    {
      Mpi::Print(comm, " (avg. reduction factor: {:.6e})\n",
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

template <typename OperType>
void GmresSolver<OperType>::Initialize() const
{
  if (!V.empty())
  {
    MFEM_ASSERT(V.Size() == max_dim + 1 && V[0].Size() == A->Height(),
                "Repeated solves with GmresSolver should not modify the operator size or "
                "restart dimension!");
    return;
  }
  if (max_dim < 0)
  {
    max_dim = max_it;
  }
  V.resize(max_dim + 1);
  for (int j = 0; j < std::min(5, max_dim + 1); j++)
  {
    V[j].SetSize(A->Height());
  }
  if (flexible)
  {
    Z.resize(max_dim + 1);
    for (int j = 0; j < std::min(5, max_dim + 1); j++)
    {
      Z[j].SetSize(A->Height());
    }
  }
  else
  {
    r.SetSize(A->Height());
  }
  H.resize((max_dim + 1) * max_dim);
  s.resize(max_dim + 1);
  cs.resize(max_dim + 1);
  sn.resize(max_dim + 1);
}

template <typename OperType>
void GmresSolver<OperType>::Mult(const VecType &x, VecType &y) const
{
  // Set up workspace.
  RealType beta = 0.0, true_beta, eps;
  MFEM_VERIFY(A, "Operator must be set for GmresSolver::Mult!");
  MFEM_ASSERT(A->Width() == x.Size() && A->Height() == b.Size(),
              "Size mismatch for GmresSolver::Mult!");
  Initialize();

  // Begin iterations.
  converged = false;
  int it = 0, restart = 0;
  if (print_opts.iterations)
  {
    Mpi::Print(comm, "{}Residual norms for {}GMRES solve\n", flexible ? "F" : "",
               std::string(tab_width + int_width - 1, ' '));
  }
  for (; it < max_it && !converged; restart++)
  {
    // Initialize.
    if (B && pc_side == PrecSide::LEFT)
    {
      if (initial_guess || restart > 0)
      {
        A->Mult(x, V[0]);
        linalg::AXPBY(1.0, b, -1.0, V[0]);
        B->Mult(V[0], r);
      }
      else
      {
        B->Mult(b, r);
        x = 0.0;
      }
    }
    else  // !B || pc_side == PrecSide::RIGHT
    {
      if (initial_guess || restart > 0)
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
    true_beta = linalg::Norml2(comm, r);
    CheckDot(true_beta, "GMRES residual norm is not valid: ||Br|| = " << true_beta << "!");
    if (it == 0)
    {
      initial_res = true_beta;
      eps = std::max(rel_tol * true_beta, abs_tol);
    }
    else if (beta > 0.0 && std::abs(beta - true_beta) > 0.1 * initial_res &&
             print_opts.warnings)
    {
      Mpi::Print(
          comm,
          "{}{}GMRES residual at restart ({:.6e}) is far from the residual norm estimate "
          "from the recursion formula ({.6e}) (initial residual = {:.6e})\n",
          std::string(tab_width, ' '), flexible ? "F" : "", true_beta, beta, initial_res);
    }
    beta = true_beta;
    if (beta < eps)
    {
      converged = true;
      break;
    }

    V[0] = 0.0;
    V[0].Add(1.0 / beta, r);
    std::fill(s.begin(), s.end(), 0.0);
    s[0] = beta;

    int j = 0;
    for (; j < max_dim && it < max_it; j++, it++)
    {
      if (print_opts.iterations)
      {
        Mpi::Print(comm, "{}{:{}d} iteration ({:d} restarts), residual {:.6e}\n", it,
                   std::string(tab_width, ' '), int_width, restart, beta);
      }
      VecType &w = V[j + 1];
      if (w.Size() == 0)
      {
        // Add storage for basis vectors in increments.
        for (int k = j + 1; k < std::min(j + 11, max_dim + 1); k++)
        {
          V[k].SetSize(A->Height());
        }
        if (flexible)
        {
          for (int k = j + 1; k < std::min(j + 11, max_dim + 1); k++)
          {
            Z[k].SetSize(A->Height());
          }
        }
      }
      if (B && pc_side == PrecSide::LEFT)
      {
        A->Mult(V[j], r);
        B->Mult(r, w);
      }
      else if (B && pc_side == PrecSide::RIGHT)
      {
        if (!flexible)
        {
          B->Mult(V[j], r);
          A->Mult(r, w);
        }
        else
        {
          B->Mult(V[j], Z[j]);
          A->Mult(Z[j], w);
        }
      }
      else
      {
        A->Mult(V[j], w);
      }

      ScalarType *Hj = H.data() + j * (max_dim + 1);
      switch (orthog_type)
      {
        case OrthogType::MGS:
          linalg::OrthogonalizeColumnMGS(comm, V, w, Hj, j + 1);
          break;
        case OrthogType::CGS:
          linalg::OrthogonalizeColumnCGS(comm, V, w, Hj, j + 1);
          break;
        case OrthogType::CGS2:
          linalg::OrthogonalizeColumnCGS(comm, V, w, Hj, j + 1, true);
          break;
      }
      Hj[j + 1] = linalg::Norml2(comm, w);
      w *= 1.0 / Hj[j + 1];

      for (int k = 0; k < j; k++)
      {
        ApplyPlaneRotation(Hj[k], Hj[k + 1], cs[k], sn[k]);
      }
      GeneratePlaneRotation(Hj[j], Hj[j + 1], cs[j], sn[j]);
      ApplyPlaneRotation(Hj[j], Hj[j + 1], cs[j], sn[j]);
      ApplyPlaneRotation(s[j], s[j + 1], cs[j], sn[j]);

      beta = std::abs(s[j + 1]);
      CheckDot(beta, "GMRES residual norm is not valid: ||Br|| = " << beta << "!");
      if (beta < eps)
      {
        converged = true;
        break;
      }
    }

    // Reconstruct the solution (for restart or due to convergence or maximum iterations).
    for (int i = j; i >= 0; i--)
    {
      ScalarType *Hi = H.data() + i * (max_dim + 1);
      s[i] /= Hi[i];
      for (int k = 0; k < i; k++)
      {
        s[k] -= Hi[k] * s[i];
      }
    }
    if (!B || pc_side == PrecSide::LEFT)
    {
      for (int k = 0; k <= j; k++)
      {
        x.Add(s[k], V[k]);
      }
    }
    else  // B && pc_side == PrecSide::RIGHT
    {
      if (!flexible)
      {
        r = 0.0;
        for (int k = 0; k <= j; k++)
        {
          r.Add(s[k], V[k]);
        }
        B->Mult(r, V[0]);
        x += V[0];
      }
      else
      {
        for (int k = 0; k <= j; k++)
        {
          x.Add(s[k], Z[k]);
        }
      }
    }
  }
  if (print_opts.iterations)
  {
    Mpi::Print(comm, "{}{:{}d} iteration ({:d} restarts), residual {:.6e}\n", it, int_width,
               std::string(tab_width, ' '), restart, beta);
  }
  if (print_opts.summary || (print_opts.warnings && !converged))
  {
    Mpi::Print(comm, "{}{}GMRES solver {} with {:d} iteration{}", flexible ? "F" : "",
               std::string(tab_width, ' '), converged ? "converged" : "did NOT converge",
               it, (it == 1) ? "" : "s");
    if (it > 0)
    {
      Mpi::Print(comm, " (avg. reduction factor: {:.6e})\n",
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
