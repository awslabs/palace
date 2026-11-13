// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_LINALG_ORTHOG_HPP
#define PALACE_LINALG_ORTHOG_HPP

#include <algorithm>
#include <array>
#include <cmath>
#include <complex>
#include <cstddef>
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

}  // namespace internal

// Classical Gram-Schmidt, reorthogonalized if refine, in the inner product of the real
// weight operator W (Euclidean without W). Returns the norms of w before and after. Each
// pass applies W once and computes all inner products, including the squared norm of w, in
// one pass over the real and imaginary parts of V with one reduction. The workspace must
// not alias any vector of V.
template <typename VecType, typename ScalarType, typename WeightOperator = Operator>
inline std::array<double, 2>
OrthogonalizeColumnCGS(MPI_Comm comm, std::span<VecType> V, ScalarType *H, bool refine,
                       ComplexVector &work, const WeightOperator *W = nullptr)
{
  MFEM_ASSERT(!V.empty(), "Missing vector for CGS orthogonalization!");
  constexpr std::size_t p = std::is_same_v<VecType, ComplexVector> ? 2 : 1;
  constexpr std::complex<double> unit[2] = {1.0, {0.0, 1.0}};
  const std::size_t m = V.size() - 1;
  std::vector<Vector> parts;
  const auto X = internal::RealParts(V, parts);  // vⱼ = Σₐ iᵃ X[p j + a]
  auto &w = V.back();
  work.SetSize(w.Size());
  work.UseDevice(w.UseDevice());
  auto apply = [&]()  // Parts of W w
  {
    std::array<const Vector *, p> Ww;
    for (std::size_t a = 0; a < p; a++)
    {
      Ww[a] = &X[p * m + a];
      if (W)
      {
        W->Mult(*Ww[a], a ? work.Imag() : work.Real());
        Ww[a] = a ? &work.Imag() : &work.Real();
      }
    }
    return Ww;
  };
  std::vector<double> dots(p * p * (m + 1));
  std::vector<std::complex<double>> h(m + 1), c(p * m);
  std::fill_n(H, m, ScalarType(0.0));
  double pre_norm_sq = 0.0;
  for (int pass = 0; pass < ((refine && m > 0) ? 2 : 1); pass++)
  {
    LocalDots(X, apply(), dots.data());
    Mpi::GlobalSum(static_cast<int>(dots.size()), dots.data(), comm);
    for (std::size_t j = 0; j <= m; j++)
    {
      // vⱼᴴ W w = Σₐ Σᵦ conj(iᵃ) iᵇ X[p j + a]ᵀ W X[p m + b].
      h[j] = 0.0;
      for (std::size_t a = 0; a < p; a++)
      {
        for (std::size_t b = 0; b < p; b++)
        {
          h[j] += std::conj(unit[a]) * unit[b] * dots[p * j + a + b * p * (m + 1)];
        }
      }
    }
    if (pass == 0)
    {
      pre_norm_sq = h[m].real();
    }
    if (m > 0)
    {
      for (std::size_t j = 0; j < m; j++)
      {
        if constexpr (p == 1)
        {
          H[j] += h[j].real();
        }
        else
        {
          H[j] += h[j];
        }
        for (std::size_t a = 0; a < p; a++)
        {
          c[p * j + a] = -h[j] * unit[a];
        }
      }
      LinearCombination(X.first(p * m), c, work);
      if constexpr (p == 1)
      {
        AXPY(1.0, work.Real(), w);
      }
      else
      {
        w.AXPY(1.0, work);
      }
    }
  }
  double norm_sq = 0.0;
  const auto Ww = apply();
  for (std::size_t a = 0; a < p; a++)
  {
    norm_sq += LocalDot(X[p * m + a], *Ww[a]);
  }
  Mpi::GlobalSum(1, &norm_sq, comm);
  return {std::sqrt(std::abs(pre_norm_sq)), std::sqrt(std::abs(norm_sq))};
}

// Gram-Schmidt of the given type, Euclidean. Returns the norm of w after.
template <typename VecType, typename ScalarType>
inline double OrthogonalizeColumn(Orthogonalization type, MPI_Comm comm,
                                  std::span<VecType> V, ScalarType *H, ComplexVector &work)
{
  return (type == Orthogonalization::MGS)
             ? OrthogonalizeColumnMGS(comm, V, H)
             : OrthogonalizeColumnCGS(comm, V, H, type == Orthogonalization::CGS2, work)[1];
}

}  // namespace palace::linalg

#endif  // PALACE_LINALG_ORTHOG_HPP
