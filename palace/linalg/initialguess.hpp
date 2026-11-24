// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_LINALG_INITIAL_GUESS_HPP
#define PALACE_LINALG_INITIAL_GUESS_HPP

#include <algorithm>
#include <array>
#include <cstddef>
#include <span>
#include <type_traits>
#include <vector>
#include <mfem.hpp>
#include "linalg/operator.hpp"
#include "linalg/orthog.hpp"
#include "linalg/vector.hpp"
#include "utils/communication.hpp"

namespace palace::linalg
{

// Copies of the latest solutions of a sequence of solves, up to a capacity, and their
// minimal-residual combination as the initial guess of the next solve.
template <OperatorType OperType>
class SolutionHistory
{
  using VecType = typename SumOperatorTraits<OperType>::VecType;
  using ScalarType = typename SumOperatorTraits<OperType>::ScalarType;

private:
  std::vector<VecType> storage;
  std::size_t newest = 0, count = 0;

  // Sets y = Σₖ c[k] V[k] in one pass over the real (and imaginary) parts of the vectors.
  static void Combine(std::span<VecType *const> V, std::span<const ScalarType> c,
                      VecType &y)
  {
    constexpr std::size_t p = std::is_same_v<VecType, ComplexVector> ? 2 : 1;
    std::vector<Vector> parts(p * V.size());
    std::vector<ScalarType> d(p * V.size());
    for (std::size_t k = 0; k < V.size(); k++)
    {
      if constexpr (p == 2)
      {
        parts[2 * k].MakeRef(V[k]->Real(), 0, V[k]->Size());
        parts[2 * k + 1].MakeRef(V[k]->Imag(), 0, V[k]->Size());
        d[2 * k] = c[k];
        d[2 * k + 1] = c[k] * ScalarType(0.0, 1.0);
      }
      else
      {
        parts[k].MakeRef(*V[k], 0, V[k]->Size());
        d[k] = c[k];
      }
    }
    LinearCombination(parts, d, y);
  }

public:
  SolutionHistory(std::size_t capacity, int n) : storage(capacity)
  {
    for (auto &v : storage)
    {
      v.SetSize(n);
      v.UseDevice(true);
    }
  }

  // Add a copy of x, replacing the oldest solution when full.
  void Push(const VecType &x)
  {
    if (!storage.empty())
    {
      newest = (newest + storage.size() - 1) % storage.size();
      storage[newest] = x;
      count = std::min(count + 1, storage.size());
    }
  }

  // Forget the stored solutions, or also release their memory.
  void Clear() { count = 0; }
  void Release()
  {
    count = 0;
    storage.clear();
  }

  std::size_t Size() const { return count; }

  // The stored solutions, the most recent first.
  std::vector<const VecType *> Vectors() const
  {
    std::vector<const VecType *> X(count);
    for (std::size_t j = 0; j < count; j++)
    {
      X[j] = &storage[(newest + j) % storage.size()];
    }
    return X;
  }

  // Sets x = Σⱼ cⱼ Xⱼ, the combination of the stored solutions X (the most recent first)
  // which minimizes the residual ||b - A x|| in the Euclidean norm (not that of a
  // preconditioned stopping test), with a QR factorization of the images A Xⱼ by classical
  // Gram-Schmidt with reorthogonalization (about three reductions per solution, and one for
  // Qᴴ b). A solution whose image is not finite, or numerically dependent on those of the
  // previous ones (its orthogonal part is at most 10⁻¹⁰ times its norm), gets cⱼ = 0. The
  // coefficients of nearly dependent images can be large, and so the rounding errors with
  // them: the combination is only kept if its residual from the images, b - Σⱼ cⱼ A Xⱼ, is
  // smaller than that of the first kept solution alone, scaled to minimize it, which is the
  // result otherwise (a few more vector updates and one more reduction). x has the size of
  // b, and is zero if all images are zero or not finite. The images and their
  // orthonormalization are only allocated during the call.
  void MinimalResidualCombination(MPI_Comm comm, const OperType &A, const VecType &b,
                                  VecType &x)
  {
    constexpr double dependence_tol = 1.0e-10;
    const std::size_t m = count;

    // The images W of the kept solutions, and Q and R (column-major, m x m) of them.
    std::vector<VecType> W(m), Q(m);
    std::vector<ScalarType> R(m * m);
    std::vector<VecType *> kept, images;
    for (std::size_t j = 0; j < m; j++)
    {
      const std::size_t n = kept.size();
      VecType &X = storage[(newest + j) % storage.size()], &q = Q[n];
      for (auto *v : {&W[n], &q})
      {
        v->SetSize(b.Size());
        v->UseDevice(true);
      }
      A.Mult(X, W[n]);
      q = W[n];
      ScalarType *Rn = R.data() + n * m;
      const auto [norm_image, norm] =
          OrthogonalizeColumnCGS(comm, std::span(Q).first(n + 1), Rn, true);
      if (!(norm > dependence_tol * norm_image))
      {
        continue;
      }
      Rn[n] = norm;
      q *= 1.0 / norm;
      kept.push_back(&X);
      images.push_back(&W[n]);
    }
    const std::size_t n = kept.size();
    if (n == 0)
    {
      x = 0.0;
      return;
    }

    // The coefficients: Qᴴ b from the inner products of the real (and imaginary) parts in
    // one pass and one reduction, then the back substitution R c = Qᴴ b. alpha scales the
    // first kept solution alone.
    std::vector<ScalarType> c(n);
    {
      constexpr std::size_t p = std::is_same_v<VecType, ComplexVector> ? 2 : 1;
      std::vector<Vector> parts;
      std::array<const Vector *, p> bp;
      if constexpr (p == 2)
      {
        bp = {&b.Real(), &b.Imag()};
      }
      else
      {
        bp = {&b};
      }
      std::vector<double> d(p * p * n);  // d[k + l p n] = parts[k] · bp[l]
      LocalDots(internal::RealParts(std::span(Q).first(n), parts), bp, d.data());
      Mpi::GlobalSum(static_cast<int>(d.size()), d.data(), comm);
      for (std::size_t k = 0; k < n; k++)
      {
        if constexpr (p == 2)
        {
          c[k] = {d[2 * k] + d[2 * k + 1 + 2 * n], d[2 * k + 2 * n] - d[2 * k + 1]};
        }
        else
        {
          c[k] = d[k];
        }
      }
    }
    const ScalarType qb = c[0], alpha = c[0] / R[0];
    for (std::size_t p = n; p-- > 0;)
    {
      for (std::size_t q = p + 1; q < n; q++)
      {
        c[p] -= R[q * m + p] * c[q];
      }
      c[p] /= R[p * m + p];
    }
    if (n > 1)
    {
      // The residuals b - Σⱼ cⱼ A Xⱼ and b - alpha A X₀ = b - (Q₀ᴴ b) Q₀, in Q.
      Combine(images, c, Q[1]);
      linalg::AXPBY(1.0, b, -1.0, Q[1]);
      linalg::AXPBY(ScalarType(1.0), b, -qb, Q[0]);
      double res[2] = {std::real(LocalDot(Q[1], Q[1])), std::real(LocalDot(Q[0], Q[0]))};
      Mpi::GlobalSum(2, res, comm);
      if (res[0] < res[1])
      {
        Combine(kept, c, x);
        return;
      }
    }
    Combine(std::span(kept).first(1), std::span(&alpha, 1), x);
  }
};

}  // namespace palace::linalg

#endif  // PALACE_LINALG_INITIAL_GUESS_HPP
