// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_LINALG_ORTHOG_HPP
#define PALACE_LINALG_ORTHOG_HPP

#include <algorithm>
#include <array>
#include <cmath>
#include <complex>
#include <cstddef>
#include <initializer_list>
#include <span>
#include <type_traits>
#include <vector>
#include "linalg/operator.hpp"
#include "linalg/vector.hpp"
#include "utils/communication.hpp"
#include "utils/labels.hpp"

namespace palace::linalg
{

//
// Gram-Schmidt orthogonalization of the last vector w of V against the others, which are
// orthonormal. The projections are stored in H, and w is not normalized.
//

// Modified Gram-Schmidt, with one reduction per basis vector. Returns the norm of w after.
template <typename VecType, typename ScalarType>
inline double OrthogonalizeColumnMGS(MPI_Comm comm, std::span<VecType> V, ScalarType *H)
{
  MFEM_ASSERT(!V.empty(), "Missing vector for MGS orthogonalization!");
  auto &w = V.back();
  for (std::size_t j = 0; j + 1 < V.size(); j++)
  {
    // Note the order is important for complex vectors.
    H[j] = Dot(comm, w, V[j]);
    w.Add(-H[j], V[j]);
  }
  return Norml2(comm, w);
}

namespace internal
{

// Views of the real and imaginary parts of complex vectors.
inline std::span<const Vector> RealParts(std::span<Vector> V, std::vector<Vector> &)
{
  return V;
}

inline std::span<const Vector> RealParts(std::span<ComplexVector> V,
                                         std::vector<Vector> &parts)
{
  parts.resize(2 * V.size());
  for (std::size_t j = 0; j < V.size(); j++)
  {
    parts[2 * j].MakeRef(V[j].Real(), 0, V[j].Size());
    parts[2 * j + 1].MakeRef(V[j].Imag(), 0, V[j].Size());
  }
  return parts;
}

// The inner products G[i + n k] = Xᵢᴴ Yₖ of the n vectors Xᵢ with the vectors Yₖ, given by
// their parts X and Y (see RealParts), with one LocalDots and one reduction.
template <typename VecType>
inline auto InnerProducts(MPI_Comm comm, std::span<const Vector> X,
                          std::span<const Vector *const> Y)
{
  std::vector<double> dots(X.size() * Y.size());
  LocalDots(X, Y, dots.data());
  Mpi::GlobalSum(static_cast<int>(dots.size()), dots.data(), comm);
  if constexpr (std::is_same_v<VecType, Vector>)
  {
    return dots;
  }
  else
  {
    // (xᵣ + i xᵢ)ᴴ (yᵣ + i yᵢ) = xᵣ · yᵣ + xᵢ · yᵢ + i (xᵣ · yᵢ - xᵢ · yᵣ).
    const std::size_t n = X.size() / 2;
    std::vector<std::complex<double>> G(n * (Y.size() / 2));
    for (std::size_t q = 0; q < G.size(); q++)
    {
      const auto d = [&](std::size_t a, std::size_t b)
      { return dots[2 * (q % n) + a + (2 * (q / n) + b) * 2 * n]; };
      G[q] = {d(0, 0) + d(1, 1), d(0, 1) - d(1, 0)};
    }
    return G;
  }
}

// w[k] = s[k] (w[k] + Σⱼ c[j + n k] vⱼ), with s[k] = 1 without s, for the n = c.size() /
// w.size() vectors vⱼ given by their parts X (see RealParts), in one pass over all of them.
// A vⱼ may also be one of the w[k], for its value before the update.
template <typename VecType, typename ScalarType>
inline void AddCombination(std::span<const Vector> X, const std::vector<ScalarType> &c,
                           std::initializer_list<VecType *> w,
                           std::initializer_list<double> s = {})
{
  constexpr std::size_t p = std::is_same_v<VecType, ComplexVector> ? 2 : 1;
  const std::size_t K = p * w.size(), n = c.size() / w.size();
  std::vector<Vector *> Y(K);
  std::vector<double> C(K * p * n), S(K, 1.0);
  for (std::size_t k = 0; k < w.size(); k++)
  {
    VecType &wk = *w.begin()[k];
    const ScalarType *ck = c.data() + n * k;
    if (s.size() > 0)
    {
      std::fill_n(S.begin() + p * k, p, s.begin()[k]);
    }
    if constexpr (p == 1)
    {
      Y[k] = &wk;
      for (std::size_t j = 0; j < n; j++)
      {
        C[k + K * j] = ck[j];
      }
    }
    else
    {
      Y[2 * k] = &wk.Real();
      Y[2 * k + 1] = &wk.Imag();
      for (std::size_t j = 0; j < n; j++)
      {
        // c (xᵣ + i xᵢ) = Re(c) xᵣ - Im(c) xᵢ + i (Im(c) xᵣ + Re(c) xᵢ).
        double *Cr = C.data() + 2 * k + K * 2 * j, *Ci = Cr + K;  // Terms of xᵣ, xᵢ
        Cr[0] = std::real(ck[j]);
        Cr[1] = std::imag(ck[j]);
        Ci[0] = -std::imag(ck[j]);
        Ci[1] = std::real(ck[j]);
      }
    }
  }
  LinearCombination(X.first(p * n), C.data(), Y, S.data(), true);
}

}  // namespace internal

// Classical Gram-Schmidt, reorthogonalized if refine, in the inner product of the real
// weight operator W (Euclidean without W). Returns the norms of w before and after. Each
// pass applies W once and computes all inner products, including the squared norm of w, in
// one pass over the real and imaginary parts of V with one reduction, then updates w in one
// pass. With refine, V is orthonormal to round-off, so the squared norm after the second
// pass is the one before it minus those of the projections, which is accurate unless the
// pass removes most of w: only then, or without refine, does the norm after take one more
// reduction (and W application). With W, the workspace is required and must not alias any
// vector of V.
template <typename VecType, typename ScalarType, typename WeightOperator = Operator>
inline std::array<double, 2>
OrthogonalizeColumnCGS(MPI_Comm comm, std::span<VecType> V, ScalarType *H, bool refine,
                       const WeightOperator *W = nullptr, ComplexVector *work = nullptr)
{
  MFEM_ASSERT(!V.empty(), "Missing vector for CGS orthogonalization!");
  MFEM_ASSERT(!W || work, "Missing workspace for weighted CGS orthogonalization!");
  constexpr std::size_t p = std::is_same_v<VecType, ComplexVector> ? 2 : 1;
  const std::size_t m = V.size() - 1;
  std::vector<Vector> parts;
  const auto X = internal::RealParts(V, parts);  // vⱼ = Σₐ iᵃ X[p j + a]
  auto &w = V.back();
  if (W)
  {
    work->SetSize(w.Size());
    work->UseDevice(w.UseDevice());
  }
  auto apply = [&]()  // Parts of W w
  {
    std::array<const Vector *, p> Ww;
    for (std::size_t a = 0; a < p; a++)
    {
      Ww[a] = &X[p * m + a];
      if (W)
      {
        W->Mult(*Ww[a], a ? work->Imag() : work->Real());
        Ww[a] = a ? &work->Imag() : &work->Real();
      }
    }
    return Ww;
  };
  std::vector<std::conditional_t<p == 1, double, std::complex<double>>> c(m);
  std::fill_n(H, m, ScalarType(0.0));
  double pre_norm_sq = 0.0, norm_sq = 0.0;
  bool accurate = false;
  for (int pass = 0; pass < ((refine && m > 0) ? 2 : 1); pass++)
  {
    const auto h = internal::InnerProducts<VecType>(comm, X, apply());  // vⱼᴴ W w
    if (pass == 0)
    {
      pre_norm_sq = std::real(h[m]);
    }
    if (m == 0)
    {
      const double norm = std::sqrt(std::abs(pre_norm_sq));
      return {norm, norm};
    }
    double proj_sq = 0.0;
    for (std::size_t j = 0; j < m; j++)
    {
      H[j] += h[j];
      c[j] = -h[j];
      proj_sq += std::norm(h[j]);
    }
    norm_sq = std::real(h[m]) - proj_sq;
    accurate = refine && (norm_sq >= 0.5 * std::real(h[m]));
    internal::AddCombination(X, c, {&w});
  }
  if (!accurate)
  {
    norm_sq = 0.0;
    const auto Ww = apply();
    for (std::size_t a = 0; a < p; a++)
    {
      norm_sq += LocalDot(X[p * m + a], *Ww[a]);
    }
    Mpi::GlobalSum(1, &norm_sq, comm);
  }
  return {std::sqrt(std::abs(pre_norm_sq)), std::sqrt(std::abs(norm_sq))};
}

// Gram-Schmidt of the given type, Euclidean. Returns the norm of w after.
template <typename VecType, typename ScalarType>
inline double OrthogonalizeColumn(Orthogonalization type, MPI_Comm comm,
                                  std::span<VecType> V, ScalarType *H)
{
  MFEM_VERIFY(type != Orthogonalization::DEFAULT, "Unresolved default orthogonalization!");
  return (type == Orthogonalization::MGS)
             ? OrthogonalizeColumnMGS(comm, V, H)
             : OrthogonalizeColumnCGS(comm, V, H, type == Orthogonalization::CGS2)[1];
}

}  // namespace palace::linalg

#endif  // PALACE_LINALG_ORTHOG_HPP
