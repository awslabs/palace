// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include <algorithm>
#include <cmath>
#include <complex>
#include <numbers>
#include <span>
#include <vector>
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

struct CountingWeightOperator
{
  mfem::DenseMatrix matrix;
  mutable int applications = 0;

  explicit CountingWeightOperator(int n) : matrix(n)
  {
    matrix = 0.0;
    for (int i = 0; i < n; i++)
    {
      matrix(i, i) = 2.0 + i;
      if (i > 0)
      {
        matrix(i, i - 1) = matrix(i - 1, i) = 0.25;
      }
    }
  }

  void Mult(const Vector &x, Vector &y) const
  {
    applications++;
    matrix.Mult(x, y);
  }
};

// Per-vector CGS of w against the first m vectors of V.
template <typename Dot>
void ReferenceCGS(MPI_Comm comm, const std::vector<Vector> &V, Vector &w, double *H,
                  std::size_t m, bool refine, const Dot &dot)
{
  std::fill_n(H, m, 0.0);
  for (int pass = 0; pass < (refine ? 2 : 1); pass++)
  {
    std::vector<double> h(m + 1);
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

TEST_CASE("OrthogonalizeColumnCGS", "[orthog][Serial][Parallel]")
{
  // CGS of the last vector of V, with one weight application per pass and one for the
  // final norm, compared with the per-vector CGS.
  const auto m = GENERATE(std::size_t{0}, std::size_t{1}, std::size_t{4});
  const bool refine = GENERATE(false, true);
  const bool weighted = GENERATE(false, true);
  const bool dependent = GENERATE(false, true);
  const auto comm = Mpi::World();
  constexpr int n = 6;
  CountingWeightOperator W(n);
  Vector reference_work(n);
  ComplexVector work;
  auto dot = [&](const Vector &x, const Vector &y)
  {
    if (!weighted)
    {
      return linalg::LocalDot(x, y);
    }
    W.Mult(x, reference_work);
    return linalg::LocalDot(reference_work, y);
  };
  auto norm = [&](const Vector &x)
  {
    auto norm_sq = dot(x, x);
    Mpi::GlobalSum(1, &norm_sq, comm);
    return std::sqrt(norm_sq);
  };

  std::vector<Vector> V;
  V.reserve(m + 1);
  std::vector<double> H(m), H_ref(m);
  for (std::size_t j = 0; j <= m; j++)
  {
    auto &v = V.emplace_back(n);
    v.Randomize(314159 + 17 * Mpi::Rank(comm) + 2 * j);
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
  Vector expected(w);
  const double pre_norm_ref = norm(w);
  ReferenceCGS(comm, V, expected, H_ref.data(), m, refine, dot);
  W.applications = 0;
  const auto [pre_norm, post_norm] = linalg::OrthogonalizeColumnCGS(
      comm, std::span(V), H.data(), refine, work, weighted ? &W : nullptr);
  if (weighted)
  {
    CHECK(W.applications == ((m > 0 && refine) ? 3 : 2));
  }
  CHECK_THAT(pre_norm, WithinRel(pre_norm_ref, 1.0e-14));
  CHECK_THAT(post_norm, WithinRel(norm(w), 1.0e-14));
  for (std::size_t j = 0; j < m; j++)
  {
    CHECK_THAT(H[j] - H_ref[j], WithinAbs(0.0, 1.0e-13));
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
  ComplexVector work;
  linalg::OrthogonalizeColumn(orthog_type, Mpi::World(), std::span(V), H, work);

  // Vector should remain unchanged. RangeEquals reads the host arrays directly, so copy
  // the device data back to the host first.
  w.HostRead();
  w_orig.HostRead();
  CHECK_THAT(w, RangeEquals(w_orig));
}

TEST_CASE("OrthogonalizeColumn Parameterized - Real 1", "[orthog][Serial][Parallel]")
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
  V.at(mpi_rank)[mpi_rank] = 1.0;

  // Pick random vector
  auto &w = V.at(mpi_size);
  w.Randomize();

  ComplexVector work;
  linalg::OrthogonalizeColumn(orthog_type, Mpi::World(), std::span(V), H.data(), work);

  // We should have zeroed out the value at "mpi_rank"
  CHECK_THAT(w[mpi_rank], WithinAbs(0.0, 1e-12));

  // Check full orthogonalization
  for (int i = 0; i < mpi_size; i++)
  {
    auto dot = linalg::Dot(Mpi::World(), w, V.at(i));
    CHECK_THAT(dot, WithinAbs(0.0, 1e-12));
  }
}

// TODO: This test could be GPU enabled, but the code below is not correct for that due to
// incorrect host syncing in set-up.
TEST_CASE("OrthogonalizeColumn Parameterized - Real 2", "[orthog][Serial][Parallel]")
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
  V[0][0] = 1.0;

  auto v0_norm = linalg::Norml2(Mpi::World(), V[0]);
  CHECK_THAT(v0_norm, WithinRel(std::sqrt(mpi_size)) || WithinAbs(0.0, 1e-14));
  V[0] *= 1 / v0_norm;

  // Vector 2:

  // [0, 1, 0, 0] on all ranks
  V[1][1] = 1.0;
  ComplexVector work;
  linalg::OrthogonalizeColumn(orthog_type, Mpi::World(), std::span(V).first(2), H.data(),
                              work);

  // Should be exact in double as multiply by zero. OrthogonalizeColumn does not
  // normalize.
  CHECK_THAT(V[1], RangeEquals(std::vector<double>{0.0, 1.0, 0.0, 0.0}));

  // Normalize by hand.
  auto v1_norm = linalg::Norml2(Mpi::World(), V[1]);
  CHECK_THAT(v1_norm, WithinRel(std::sqrt(mpi_size)) || WithinAbs(0.0, 1e-14));
  V[1] *= 1 / v1_norm;

  // Vector 3: Non-trivial vector to orthogonalize
  auto &w = V[2];
  auto d_v = w.Write();
  mfem::forall(w.Size(), [=] MFEM_HOST_DEVICE(int i) { d_v[i] = mpi_rank + i; });

  linalg::OrthogonalizeColumn(orthog_type, Mpi::World(), std::span(V), H.data(), work);

  // Check orthogonality
  auto dot0 = linalg::Dot(Mpi::World(), w, V[0]);
  auto dot1 = linalg::Dot(Mpi::World(), w, V[1]);

  CHECK_THAT(dot0, WithinAbs(0.0, 1e-12));
  CHECK_THAT(dot1, WithinAbs(0.0, 1e-12));

  // Componenet 2 & 3 are untouched by orthogonalization so should be the same as before.
  CHECK_THAT(w[2], WithinRel(mpi_rank + 2.));
  CHECK_THAT(w[3], WithinRel(mpi_rank + 3.));

  // H should be orthogonalization coefficients
  CHECK_THAT(H[0], WithinRel((mpi_size * (mpi_size - 1.0) / (2 * v0_norm))));
  CHECK_THAT(H[1], WithinRel((mpi_size * (mpi_size + 1.0) / (2 * v1_norm))));
}

TEST_CASE("OrthogonalizeColumn Parameterized - Complex 1", "[orthog][Serial][Parallel]")
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
  V.at(mpi_rank).Real()[mpi_rank] = std::cos(theta);
  V.at(mpi_rank).Imag()[mpi_rank] = std::sin(theta);

  // Pick random vector to orthogonalize
  auto &w = V.at(mpi_size);
  w.Real().Randomize();
  w.Imag().Randomize();
  const ComplexVector w0(w);

  ComplexVector work;
  const double norm =
      linalg::OrthogonalizeColumn(orthog_type, Mpi::World(), std::span(V), H.data(), work);
  CHECK_THAT(norm, WithinRel(linalg::Norml2(Mpi::World(), w), 1e-14));

  // We should have zeroed out the value at "mpi_rank"
  CHECK_THAT(w.Real()[mpi_rank], WithinAbs(0.0, 1e-12));
  CHECK_THAT(w.Imag()[mpi_rank], WithinAbs(0.0, 1e-12));

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

  mfem::DenseMatrix W(3, 3);
  W = 0.0;
  W(0, 0) = 2.0;
  W(1, 0) = 2.0;
  W(0, 1) = 2.0;
  W(1, 1) = 1.0;
  W(2, 2) = 2.0;

  // Create basis vectors
  std::vector<Vector> V(3);
  for (auto &v : V)
  {
    v.UseDevice(true);
    v.SetSize(3);
    v = 0.0;
  }
  V[0][0] = 1.0 / std::numbers::sqrt2;  // Normalized w.r.t W
  V[1][2] = 1.0 / std::numbers::sqrt2;

  auto &w = V[2];
  // Fixed seed for determinism. A seedless Randomize() uses time(0) (see
  // mfem::Vector::Randomize), making the test input differ every run.
  w.Randomize(314159);

  std::vector<double> H(2, 0.0);

  ComplexVector work;
  linalg::OrthogonalizeColumnCGS(Mpi::World(), std::span(V), H.data(), refine, work, &W);

  // Check orthogonality with respect to weight matrix
  Vector WVj(3);  // Temporary workspace
  W.Mult(V[0], WVj);
  double dot0 = linalg::Dot(Mpi::World(), w, WVj);
  W.Mult(V[1], WVj);
  double dot1 = linalg::Dot(Mpi::World(), w, WVj);

  CHECK_THAT(dot0, WithinAbs(0.0, 1e-12));
  CHECK_THAT(dot1, WithinAbs(0.0, 1e-12));
}

TEST_CASE("OrthogonalizeColumn Weighted - Complex 1", "[orthog][Serial]")
{
  const bool refine = GENERATE(false, true);

  mfem::DenseMatrix W(3, 3);
  W = 0.0;
  W(0, 0) = 2.0;
  W(1, 0) = 2.0;
  W(0, 1) = 2.0;
  W(1, 1) = 1.0;
  W(2, 2) = 2.0;

  // Create basis vectors
  std::vector<ComplexVector> V(3);
  for (auto &v : V)
  {
    v.UseDevice(true);
    v.SetSize(3);
    v = 0.0;
  }
  V[0].Real()[0] = 1.0 / std::numbers::sqrt2;
  V[1].Imag()[2] = 1.0 / std::numbers::sqrt2;

  auto &w = V[2];
  // Fixed seeds for determinism (see note in the real-valued test above). Real and
  // imaginary parts use distinct seeds so they are not identical sequences.
  w.Real().Randomize(314159);
  w.Imag().Randomize(314160);

  std::vector<std::complex<double>> H(2, 0.0);

  ComplexVector work;
  linalg::OrthogonalizeColumnCGS(Mpi::World(), std::span(V), H.data(), refine, work, &W);

  auto W_wrap = ComplexWrapperOperator(&W, nullptr);

  // Check orthogonality with respect to weight matrix
  ComplexVector WVj(3);  // Temporary workspace
  W_wrap.Mult(V[0], WVj);
  auto dot0 = linalg::Dot(Mpi::World(), w, WVj);

  W_wrap.Mult(V[1], WVj);
  auto dot1 = linalg::Dot(Mpi::World(), w, WVj);

  CHECK_THAT(dot0.real(), WithinAbs(0.0, 1e-12));
  CHECK_THAT(dot0.imag(), WithinAbs(0.0, 1e-12));
  CHECK_THAT(dot1.real(), WithinAbs(0.0, 1e-12));
  CHECK_THAT(dot1.imag(), WithinAbs(0.0, 1e-12));
}
