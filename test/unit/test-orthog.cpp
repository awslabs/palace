// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include <algorithm>
#include <array>
#include <cmath>
#include <complex>
#include <numbers>
#include <span>
#include <type_traits>
#include <utility>
#include <vector>
#include <catch2/catch_template_test_macros.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <catch2/matchers/catch_matchers_range_equals.hpp>
#include "linalg/operator.hpp"
#include "linalg/orthog.hpp"
#include "linalg/vector.hpp"
#include "utils/communication.hpp"
#include "utils/labels.hpp"

using namespace palace;
using namespace Catch::Matchers;
using namespace Catch;

namespace
{

// A dense weight operator, applied on the device, which counts its applications.
struct CountingWeightOperator
{
  int n;
  Vector matrix;  // Column-major
  mutable int applications = 0;

  explicit CountingWeightOperator(int n) : n(n), matrix(n * n)
  {
    matrix = 0.0;
    for (int i = 0; i < n; i++)
    {
      matrix[i + n * i] = 2.0 + i;
      if (i > 0)
      {
        matrix[i + n * (i - 1)] = matrix[i - 1 + n * i] = 0.25;
      }
    }
    matrix.UseDevice(true);
  }

  explicit CountingWeightOperator(const mfem::DenseMatrix &A)
    : n(A.Height()), matrix(A.Height() * A.Width())
  {
    std::copy_n(A.GetData(), matrix.Size(), matrix.HostWrite());
    matrix.UseDevice(true);
  }

  void Mult(const Vector &x, Vector &y) const
  {
    applications++;
    const int n = this->n;
    const auto *A = matrix.Read();
    const auto *X = x.Read();
    auto *Y = y.Write();
    mfem::forall(n,
                 [=] MFEM_HOST_DEVICE(int i)
                 {
                   double sum = 0.0;
                   for (int j = 0; j < n; j++)
                   {
                     sum += A[i + n * j] * X[j];
                   }
                   Y[i] = sum;
                 });
  }

  void Mult(const ComplexVector &x, ComplexVector &y) const
  {
    Mult(x.Real(), y.Real());
    Mult(x.Imag(), y.Imag());
  }
};

// Random entries in [-1, 1] (and imaginary parts) with a fixed seed.
void Randomize(Vector &v, int seed)
{
  v.Randomize(seed);
}

void Randomize(ComplexVector &v, int seed)
{
  v.Real().Randomize(seed);
  v.Imag().Randomize(seed + 1);
}

// Per-vector CGS of w against the first m vectors of V, with hⱼ = dot(w, vⱼ) = vⱼᴴ W w.
template <typename VecType, typename ScalarType, typename Dot>
void ReferenceCGS(MPI_Comm comm, const std::vector<VecType> &V, VecType &w, ScalarType *H,
                  std::size_t m, bool refine, const Dot &dot)
{
  std::fill_n(H, m, ScalarType(0.0));
  for (int pass = 0; pass < (refine ? 2 : 1); pass++)
  {
    std::vector<ScalarType> h(m + 1);
    for (std::size_t j = 0; j < m; j++)
    {
      h[j] = dot(w, V[j]);
    }
    Mpi::GlobalSum(static_cast<int>(m), h.data(), comm);
    for (std::size_t j = 0; j < m; j++)
    {
      H[j] += h[j];
      w.Add(-h[j], V[j]);
    }
  }
}

}  // namespace

TEMPLATE_TEST_CASE("OrthogonalizeColumnCGS", "[orthog][Serial][Parallel][GPU]", Vector,
                   ComplexVector)
{
  // CGS of the last vector of V, with one weight application per pass and one for the
  // final norm (W applied to the real and imaginary parts for complex vectors), compared
  // with the per-vector CGS.
  using VecType = TestType;
  using ScalarType = decltype(linalg::LocalDot(std::declval<const VecType &>(),
                                               std::declval<const VecType &>()));
  constexpr int p = std::is_same_v<VecType, ComplexVector> ? 2 : 1;
  const auto m = GENERATE(std::size_t{0}, std::size_t{1}, std::size_t{4});
  const bool refine = GENERATE(false, true);
  const bool weighted = GENERATE(false, true);
  const bool dependent = GENERATE(false, true);
  const auto comm = Mpi::World();
  constexpr int n = 6;
  CountingWeightOperator W(n);
  VecType reference_work(n);
  reference_work.UseDevice(true);
  ComplexVector work;
  auto dot = [&](const VecType &x, const VecType &y)
  {
    if (!weighted)
    {
      return linalg::LocalDot(x, y);
    }
    W.Mult(x, reference_work);
    return linalg::LocalDot(reference_work, y);
  };
  auto norm = [&](const VecType &x)
  {
    auto norm_sq = dot(x, x);
    Mpi::GlobalSum(1, &norm_sq, comm);
    return std::sqrt(std::abs(norm_sq));
  };

  std::vector<VecType> V;
  V.reserve(m + 1);
  std::vector<ScalarType> H(m), H_ref(m);
  for (std::size_t j = 0; j <= m; j++)
  {
    auto &v = V.emplace_back(n);
    v.UseDevice(true);
    Randomize(v, 314159 + 17 * Mpi::Rank(comm) + 2 * static_cast<int>(j));
    if (j < m)
    {
      ReferenceCGS(comm, V, v, H.data(), j, true, dot);
      v *= 1.0 / norm(v);
    }
  }
  auto &w = V.back();
  if (dependent && m > 0)
  {
    // Nearly in the span of the others.
    w *= 1.0e-10;
    w.Add(1.0, V[0]);
  }
  VecType expected(w);
  expected.UseDevice(true);
  const double pre_norm_ref = norm(w);
  ReferenceCGS(comm, V, expected, H_ref.data(), m, refine, dot);
  W.applications = 0;
  const auto [pre_norm, post_norm] = linalg::OrthogonalizeColumnCGS(
      comm, std::span(V), H.data(), refine, weighted ? &W : nullptr, &work);
  if (weighted)
  {
    CHECK(W.applications == p * ((m == 0) ? 1 : (refine ? 3 : 2)));
  }
  CHECK_THAT(pre_norm, WithinRel(pre_norm_ref, 1.0e-14));
  CHECK_THAT(post_norm, WithinRel(norm(w), 1.0e-14));
  for (std::size_t j = 0; j < m; j++)
  {
    CHECK_THAT(std::abs(H[j] - H_ref[j]), WithinAbs(0.0, 1.0e-13));
  }
  expected.Add(-1.0, w);
  CHECK_THAT(linalg::Norml2(comm, expected), WithinAbs(0.0, 1.0e-13));
}

TEST_CASE("OrthogonalizeColumn - Real Empty", "[orthog][Serial][Parallel][GPU]")
{
  auto orthog_type =
      GENERATE(Orthogonalization::MGS, Orthogonalization::CGS, Orthogonalization::CGS2);

  int mpi_rank = Mpi::Rank(Mpi::World());

  std::vector<Vector> V(1);
  auto &w = V[0];
  w.SetSize(4);
  w.UseDevice(true);
  auto d_v = w.Write();
  mfem::forall(w.Size(), [=] MFEM_HOST_DEVICE(int i) { d_v[i] = mpi_rank + i; });

  Vector w_orig = w;  // Copy for comparison

  double H[1];  // Not used without other vectors
  linalg::OrthogonalizeColumn(orthog_type, Mpi::World(), std::span(V), H);

  // Vector should remain unchanged. RangeEquals reads the host arrays directly, so copy
  // the device data back to the host first.
  w.HostRead();
  w_orig.HostRead();
  CHECK_THAT(w, RangeEquals(w_orig));
}

TEST_CASE("OrthogonalizeColumn Parameterized - Real 1",
          "[orthog][Serial][Parallel][GPU]")
{
  auto orthog_type =
      GENERATE(Orthogonalization::MGS, Orthogonalization::CGS, Orthogonalization::CGS2);

  int mpi_rank = Mpi::Rank(Mpi::World());
  int mpi_size = Mpi::Size(Mpi::World());

  std::vector<double> H(mpi_size, 0.0);  // Orthogonalization coefficients
  std::vector<Vector> V(mpi_size + 1);
  for (auto &v : V)
  {
    v.UseDevice(true);
    v.SetSize(mpi_size + 1);
    auto d_v = v.Write();
    mfem::forall(v.Size(), [=] MFEM_HOST_DEVICE(int i) { d_v[i] = 0.0; });
  }

  // Rank-wise orthogonal basis vectors
  V.at(mpi_rank).HostReadWrite()[mpi_rank] = 1.0;

  // Pick random vector
  auto &w = V.at(mpi_size);
  w.Randomize();

  linalg::OrthogonalizeColumn(orthog_type, Mpi::World(), std::span(V), H.data());

  // We should have zeroed out the value at "mpi_rank"
  CHECK_THAT(w.HostRead()[mpi_rank], WithinAbs(0.0, 1e-12));

  // Check full orthogonalization
  for (int i = 0; i < mpi_size; i++)
  {
    auto dot = linalg::Dot(Mpi::World(), w, V.at(i));
    CHECK_THAT(dot, WithinAbs(0.0, 1e-12));
  }
}

TEST_CASE("OrthogonalizeColumn Parameterized - Real 2",
          "[orthog][Serial][Parallel][GPU]")
{
  auto orthog_type =
      GENERATE(Orthogonalization::MGS, Orthogonalization::CGS, Orthogonalization::CGS2);

  int mpi_rank = Mpi::Rank(Mpi::World());
  int mpi_size = Mpi::Size(Mpi::World());

  // Orthogonalization coefficients
  std::array<double, 2> H{};

  // Create basis vectors
  std::vector<Vector> V(3);
  for (auto &v : V)
  {
    v.UseDevice(true);
    v.SetSize(4);
    auto d_v = v.Write();
    mfem::forall(v.Size(), [=] MFEM_HOST_DEVICE(int i) { d_v[i] = 0.0; });
  }

  // Vector 1:

  // [1, 0, 0, 0] on all ranks so globally ~[1,0,0,0,1,0,0,0,....] if multiple mpi ranks
  V[0].HostReadWrite()[0] = 1.0;

  auto v0_norm = linalg::Norml2(Mpi::World(), V[0]);
  CHECK_THAT(v0_norm, WithinRel(std::sqrt(mpi_size)) || WithinAbs(0.0, 1e-14));
  V[0] *= 1 / v0_norm;

  // Vector 2:

  // [0, 1, 0, 0] on all ranks
  V[1].HostReadWrite()[1] = 1.0;
  linalg::OrthogonalizeColumn(orthog_type, Mpi::World(), std::span(V).first(2), H.data());

  // Should be exact in double as multiply by zero. OrthogonalizeColumn does not
  // normalize. RangeEquals reads the host arrays directly.
  V[1].HostRead();
  CHECK_THAT(V[1], RangeEquals(std::vector<double>{0.0, 1.0, 0.0, 0.0}));

  // Normalize by hand.
  auto v1_norm = linalg::Norml2(Mpi::World(), V[1]);
  CHECK_THAT(v1_norm, WithinRel(std::sqrt(mpi_size)) || WithinAbs(0.0, 1e-14));
  V[1] *= 1 / v1_norm;

  // Vector 3: Non-trivial vector to orthogonalize
  auto &w = V[2];
  auto d_v = w.Write();
  mfem::forall(w.Size(), [=] MFEM_HOST_DEVICE(int i) { d_v[i] = mpi_rank + i; });

  linalg::OrthogonalizeColumn(orthog_type, Mpi::World(), std::span(V), H.data());

  // Check orthogonality
  auto dot0 = linalg::Dot(Mpi::World(), w, V[0]);
  auto dot1 = linalg::Dot(Mpi::World(), w, V[1]);

  CHECK_THAT(dot0, WithinAbs(0.0, 1e-12));
  CHECK_THAT(dot1, WithinAbs(0.0, 1e-12));

  // Componenet 2 & 3 are untouched by orthogonalization so should be the same as before.
  const auto *h_w = w.HostRead();
  CHECK_THAT(h_w[2], WithinRel(mpi_rank + 2.));
  CHECK_THAT(h_w[3], WithinRel(mpi_rank + 3.));

  // H should be orthogonalization coefficients
  CHECK_THAT(H[0], WithinRel((mpi_size * (mpi_size - 1.0) / (2 * v0_norm))));
  CHECK_THAT(H[1], WithinRel((mpi_size * (mpi_size + 1.0) / (2 * v1_norm))));
}

TEST_CASE("OrthogonalizeColumn Parameterized - Complex 1",
          "[orthog][Serial][Parallel][GPU]")
{
  auto orthog_type =
      GENERATE(Orthogonalization::MGS, Orthogonalization::CGS, Orthogonalization::CGS2);

  int mpi_rank = Mpi::Rank(Mpi::World());
  int mpi_size = Mpi::Size(Mpi::World());

  std::vector<std::complex<double>> H(mpi_size, 0.0);  // Orthogonalization coefficients
  std::vector<ComplexVector> V(mpi_size + 1);
  for (auto &v : V)
  {
    v.UseDevice(true);
    v.SetSize(mpi_size + 1);
    v = 0.0;
  }

  // Rank-wise orthogonal basis vectors
  auto theta = (2. * std::numbers::pi * mpi_rank) / mpi_size;
  V.at(mpi_rank).Real().HostReadWrite()[mpi_rank] = std::cos(theta);
  V.at(mpi_rank).Imag().HostReadWrite()[mpi_rank] = std::sin(theta);

  // Pick random vector to orthogonalize
  auto &w = V.at(mpi_size);
  w.Real().Randomize();
  w.Imag().Randomize();
  const ComplexVector w0(w);

  const double norm =
      linalg::OrthogonalizeColumn(orthog_type, Mpi::World(), std::span(V), H.data());
  CHECK_THAT(norm, WithinRel(linalg::Norml2(Mpi::World(), w), 1e-14));

  // We should have zeroed out the value at "mpi_rank"
  CHECK_THAT(w.Real().HostRead()[mpi_rank], WithinAbs(0.0, 1e-12));
  CHECK_THAT(w.Imag().HostRead()[mpi_rank], WithinAbs(0.0, 1e-12));

  // The coefficients are the projections of the original vector.
  for (int i = 0; i < mpi_size; i++)
  {
    CHECK_THAT(std::abs(H[i] - linalg::Dot(Mpi::World(), w0, V.at(i))),
               WithinAbs(0.0, 1e-12));
  }

  // Check full orthogonalization
  for (int i = 0; i < mpi_size; i++)
  {
    auto dot = linalg::Dot(Mpi::World(), w, V.at(i));
    CHECK_THAT(std::abs(dot), WithinAbs(0.0, 1e-12));
  }
}

TEST_CASE("OrthogonalizeColumn Weighted - Real 1", "[orthog][Serial]")
{
  const bool refine = GENERATE(false, true);

  // Symmetric positive definite, so a valid inner product.
  mfem::DenseMatrix A(3, 3);
  A = 0.0;
  A(0, 0) = 2.0;
  A(1, 0) = 1.0;
  A(0, 1) = 1.0;
  A(1, 1) = 2.0;
  A(2, 2) = 2.0;
  const CountingWeightOperator W(A);

  // Create basis vectors
  std::vector<Vector> V(3);
  for (auto &v : V)
  {
    v.UseDevice(true);
    v.SetSize(3);
    v = 0.0;
  }
  V[0].HostReadWrite()[0] = 1.0 / std::numbers::sqrt2;  // Normalized w.r.t W
  V[1].HostReadWrite()[2] = 1.0 / std::numbers::sqrt2;

  auto &w = V[2];
  // Fixed seed for determinism. A seedless Randomize() uses time(0) (see
  // mfem::Vector::Randomize), making the test input differ every run.
  w.Randomize(314159);
  Vector Ww0(3);
  W.Mult(w, Ww0);
  const double pre_norm_ref = std::sqrt(linalg::Dot(Mpi::World(), w, Ww0));
  const std::array<double, 2> H_ref = {linalg::Dot(Mpi::World(), Ww0, V[0]),
                                       linalg::Dot(Mpi::World(), Ww0, V[1])};

  std::vector<double> H(2, 0.0);

  ComplexVector work;
  const auto [pre_norm, post_norm] = linalg::OrthogonalizeColumnCGS(
      Mpi::World(), std::span(V), H.data(), refine, &W, &work);

  // Check orthogonality with respect to weight matrix
  Vector WVj(3);  // Temporary workspace
  W.Mult(V[0], WVj);
  double dot0 = linalg::Dot(Mpi::World(), w, WVj);
  W.Mult(V[1], WVj);
  double dot1 = linalg::Dot(Mpi::World(), w, WVj);

  CHECK_THAT(dot0, WithinAbs(0.0, 1e-12));
  CHECK_THAT(dot1, WithinAbs(0.0, 1e-12));

  // The coefficients and norms are in the W inner product.
  W.Mult(w, WVj);
  CHECK_THAT(pre_norm, WithinRel(pre_norm_ref, 1e-14));
  CHECK_THAT(post_norm, WithinRel(std::sqrt(linalg::Dot(Mpi::World(), w, WVj)), 1e-14));
  CHECK_THAT(H[0], WithinRel(H_ref[0], 1e-14));
  CHECK_THAT(H[1], WithinRel(H_ref[1], 1e-14));
}

TEST_CASE("OrthogonalizeColumn Weighted - Complex 1", "[orthog][Serial]")
{
  const bool refine = GENERATE(false, true);

  // Symmetric positive definite, so a valid inner product.
  mfem::DenseMatrix A(3, 3);
  A = 0.0;
  A(0, 0) = 2.0;
  A(1, 0) = 1.0;
  A(0, 1) = 1.0;
  A(1, 1) = 2.0;
  A(2, 2) = 2.0;
  const CountingWeightOperator W(A);

  // Create basis vectors
  std::vector<ComplexVector> V(3);
  for (auto &v : V)
  {
    v.UseDevice(true);
    v.SetSize(3);
    v = 0.0;
  }
  V[0].Real().HostReadWrite()[0] = 1.0 / std::numbers::sqrt2;
  V[1].Imag().HostReadWrite()[2] = 1.0 / std::numbers::sqrt2;

  auto &w = V[2];
  // Fixed seeds for determinism (see note in the real-valued test above). Real and
  // imaginary parts use distinct seeds so they are not identical sequences.
  w.Real().Randomize(314159);
  w.Imag().Randomize(314160);
  ComplexVector Ww0(3);
  W.Mult(w, Ww0);
  const double pre_norm_ref = std::sqrt(linalg::Dot(Mpi::World(), w, Ww0).real());
  const std::array<std::complex<double>, 2> H_ref = {linalg::Dot(Mpi::World(), Ww0, V[0]),
                                                     linalg::Dot(Mpi::World(), Ww0, V[1])};

  std::vector<std::complex<double>> H(2, 0.0);

  ComplexVector work;
  const auto [pre_norm, post_norm] = linalg::OrthogonalizeColumnCGS(
      Mpi::World(), std::span(V), H.data(), refine, &W, &work);

  // Check orthogonality with respect to weight matrix
  ComplexVector WVj(3);  // Temporary workspace
  W.Mult(V[0], WVj);
  auto dot0 = linalg::Dot(Mpi::World(), w, WVj);

  W.Mult(V[1], WVj);
  auto dot1 = linalg::Dot(Mpi::World(), w, WVj);

  CHECK_THAT(dot0.real(), WithinAbs(0.0, 1e-12));
  CHECK_THAT(dot0.imag(), WithinAbs(0.0, 1e-12));
  CHECK_THAT(dot1.real(), WithinAbs(0.0, 1e-12));
  CHECK_THAT(dot1.imag(), WithinAbs(0.0, 1e-12));

  // The coefficients and norms are in the W inner product.
  W.Mult(w, WVj);
  CHECK_THAT(pre_norm, WithinRel(pre_norm_ref, 1e-14));
  CHECK_THAT(post_norm,
             WithinRel(std::sqrt(linalg::Dot(Mpi::World(), w, WVj).real()), 1e-14));
  CHECK_THAT(std::abs(H[0] - H_ref[0]), WithinAbs(0.0, 1e-14 * std::abs(H_ref[0])));
  CHECK_THAT(std::abs(H[1] - H_ref[1]), WithinAbs(0.0, 1e-14 * std::abs(H_ref[1])));
}
