// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include <cmath>
#include <complex>
#include <limits>
#include <type_traits>
#include <vector>
#include <Eigen/Dense>
#include "linalg/operator.hpp"
#include "linalg/vector.hpp"
#include "utils/communication.hpp"

namespace palace::linalg
{

namespace
{

template <typename VecType, typename OperType>
void MultAdjoint(const OperType &A, const VecType &x, VecType &y)
{
  if constexpr (std::is_same_v<OperType, ComplexOperator>)
  {
    A.MultHermitianTranspose(x, y);
  }
  else
  {
    A.MultTranspose(x, y);
  }
}

template <typename VecType, typename OperType>
double LanczosSpectralNorm(MPI_Comm comm, const OperType &A, bool herm, VecType *x0,
                           double tol, int max_it)
{
  // Lanczos iteration for the eigenvalue of largest magnitude of B = A, or B = Aᴴ A when A
  // is not Hermitian. The extreme Ritz values converge from inside the spectrum, and much
  // faster than the power iteration. Without reorthogonalization, the basis loses
  // orthogonality only along converged Ritz vectors, which does not change the extreme
  // Ritz values. The basis restarts from the current Ritz vector when it is full.
  MFEM_VERIFY(A.Height() == A.Width(), "Spectral norm requires a square matrix!");
  constexpr int max_basis = 8;

  // A starting vector gets a small random component, so that an eigenvector which becomes
  // dominant but is orthogonal to it is still found after a few estimates.
  constexpr double start_random = 1.0e-2;

  const int n = A.Height();
  std::vector<VecType> V;
  V.reserve(max_basis);
  auto AddVector = [&]() -> VecType &
  {
    V.emplace_back(n);
    V.back().UseDevice(true);
    return V.back();
  };
  VecType w(n), t;
  w.UseDevice(true);
  if (!herm)
  {
    t.SetSize(n);
    t.UseDevice(true);
  }
  {
    // The choice is collective, since a process with no rows trivially matches any size.
    int warm = (x0 && x0->Size() == n);
    Mpi::GlobalMin(1, &warm, comm);
    VecType &v = AddVector();
    SetRandom(comm, v);
    if (warm)
    {
      v *= start_random / Norml2(comm, v);
      AXPY(1.0 / Norml2(comm, *x0), *x0, v);
    }
    Normalize(comm, v);
  }

  // Iterate until the Ritz value of largest magnitude has converged. Both ends of the Ritz
  // spectrum move outward, so for an indefinite operator the other end also has to converge
  // while it could still become the largest in magnitude. After a restart, the first Ritz
  // value is the Rayleigh quotient of the restart vector, so convergence is only checked
  // from the second iteration.
  std::vector<double> alpha, beta;
  Eigen::VectorXd s;
  double l = 0.0, l_min0 = 0.0, l_max0 = 0.0;
  int it = 0, k = 0;
  bool converged = false;
  auto RitzVector = [&](VecType &y)
  {
    if (k == 0)
    {
      y = V[0];
      return;
    }
    y = 0.0;
    for (int i = 0; i < k; i++)
    {
      AXPY(s(i), V[i], y);
    }
    Normalize(comm, y);
  };
  while (it < max_it)
  {
    const VecType &v = V[k];
    if (herm)
    {
      A.Mult(v, w);
    }
    else
    {
      A.Mult(v, t);
      MultAdjoint(A, t, w);
    }
    alpha.push_back(std::real(Dot(comm, w, v)));
    AXPY(-alpha.back(), v, w);
    if (k > 0)
    {
      AXPY(-beta.back(), V[k - 1], w);
    }
    const double b = Norml2(comm, w);
    it++;
    k++;
    if (!std::isfinite(alpha.back()) || !std::isfinite(b))
    {
      l = std::numeric_limits<double>::quiet_NaN();
      break;
    }

    // Extreme Ritz values, and the coefficients in the basis of the Ritz vector of the
    // largest in magnitude, plus the one at the other end when that end is tested (so that
    // a restart keeps both).
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eig;
    eig.computeFromTridiagonal(Eigen::Map<const Eigen::VectorXd>(alpha.data(), k),
                               Eigen::Map<const Eigen::VectorXd>(beta.data(), k - 1),
                               Eigen::ComputeEigenvectors);
    const double l_min = eig.eigenvalues()(0), l_max = eig.eigenvalues()(k - 1);
    l = std::max(std::abs(l_min), std::abs(l_max));
    const bool min_major = (-l_min > l_max);
    const bool test_min = min_major || (l_max > 0.0 && l_min < -0.5 * l);
    const bool test_max = !min_major || (l_min < 0.0 && l_max > 0.5 * l);
    s = eig.eigenvectors().col(min_major ? 0 : k - 1);
    if (k > 1 && test_min && test_max)
    {
      s += eig.eigenvectors().col(min_major ? k - 1 : 0);
    }
    if ((k > 1 && (!test_min || std::abs(l_min - l_min0) <= tol * l) &&
         (!test_max || std::abs(l_max - l_max0) <= tol * l)) ||
        b <= 1.0e-12 * l)
    {
      converged = true;
      break;
    }
    l_min0 = l_min;
    l_max0 = l_max;
    if (k < max_basis)
    {
      beta.push_back(b);
      VecType &u = (k < static_cast<int>(V.size())) ? V[k] : AddVector();
      u = w;
      u *= 1.0 / b;
    }
    else
    {
      RitzVector(w);
      V[0] = w;
      alpha.clear();
      beta.clear();
      k = 0;
    }
  }
  const double norm = herm ? l : std::sqrt(l);
  if (!converged)
  {
    Mpi::Warning(comm,
                 "Lanczos iteration did not converge in {:d} iterations, norm = {:.3e}!\n",
                 it, norm);
  }
  if (x0 && std::isfinite(l))
  {
    x0->SetSize(n);
    x0->UseDevice(true);
    RitzVector(*x0);
  }
  return norm;
}

}  // namespace

double SpectralNorm(MPI_Comm comm, const Operator &A, bool sym, Vector *x0, double tol,
                    int max_it)
{
  return LanczosSpectralNorm(comm, A, sym, x0, tol, max_it);
}

double SpectralNorm(MPI_Comm comm, const ComplexOperator &A, bool herm, ComplexVector *x0,
                    double tol, int max_it)
{
  return LanczosSpectralNorm(comm, A, herm, x0, tol, max_it);
}

}  // namespace palace::linalg
