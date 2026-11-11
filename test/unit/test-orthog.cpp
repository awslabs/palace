// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include <complex>
#include <memory>
#include <numbers>
#include <type_traits>
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

  void Mult(const ComplexVector &x, ComplexVector &y) const
  {
    applications++;
    matrix.Mult(x.Real(), y.Real());
    matrix.Mult(x.Imag(), y.Imag());
  }
};

}  // namespace

TEMPLATE_TEST_CASE("Weighted CGS reuses the weight application within each pass",
                   "[orthog][Serial][Parallel]", Vector, ComplexVector)
{
  using VecType = TestType;
  using ScalarType =
      std::conditional_t<std::is_same_v<VecType, Vector>, double, std::complex<double>>;
  const auto m = GENERATE(std::size_t{0}, std::size_t{1}, std::size_t{4});
  const bool refine = GENERATE(false, true);
  const auto comm = Mpi::World();
  constexpr int n = 6;
  CountingWeightOperator W(n);
  VecType reference_work(n), work;
  auto dot = [&W, &reference_work](const VecType &x, const VecType &y)
  {
    W.Mult(x, reference_work);
    return linalg::LocalDot(reference_work, y);
  };

  std::vector<VecType> V;
  V.reserve(m + 1);
  std::vector<ScalarType> H(m), H_ref(m);
  for (std::size_t j = 0; j <= m; j++)
  {
    auto &v = V.emplace_back(n);
    const auto seed = 314159 + 17 * Mpi::Rank(comm) + 2 * j;
    if constexpr (std::is_same_v<VecType, Vector>)
    {
      v.Randomize(seed);
    }
    else
    {
      v.Real().Randomize(seed);
      v.Imag().Randomize(seed + 1);
    }
    if (j < m)
    {
      linalg::OrthogonalizeColumnMGS(comm, V, v, H.data(), j, dot);
      auto norm_sq = dot(v, v);
      Mpi::GlobalSum(1, &norm_sq, comm);
      v *= 1.0 / std::sqrt(std::abs(norm_sq));
    }
  }

  // Match the PROM use case, where the vector being orthogonalized is V.back().
  auto &w = V.back();
  VecType expected(w);
  linalg::OrthogonalizeColumnCGS(comm, V, expected, H_ref.data(), m, refine, dot);
  W.applications = 0;
  linalg::OrthogonalizeColumnWeightedCGS(comm, V, w, H.data(), m, W, work, refine);
  CHECK(W.applications == (m == 0 ? 0 : (refine ? 2 : 1)));

  for (std::size_t j = 0; j < m; j++)
  {
    CHECK_THAT(std::abs(H[j] - H_ref[j]), WithinAbs(0.0, 1.0e-13));
    auto projection = dot(w, V[j]);
    Mpi::GlobalSum(1, &projection, comm);
    CHECK_THAT(std::abs(projection), WithinAbs(0.0, 1.0e-12));
  }
  expected.Add(-1.0, w);
  CHECK_THAT(linalg::Norml2(comm, expected), WithinAbs(0.0, 1.0e-13));
}

class RealWeightedInnerProduct
{
  // Choose generic operator, although can improve by refining for specialized type.
  std::shared_ptr<Operator> weight_op;

  // Don't have inner product wrapper / implemented in Operator, so need to allocate a
  // vector as a workspace. (TODO: Optimize this away).
  mutable Vector v_workspace = {};

  void SetWorkspace(const Vector &blueprint) const
  {
    v_workspace.SetSize(blueprint.Size());
    v_workspace.UseDevice(blueprint.UseDevice());
  }

public:
  template <typename OpType>
  explicit RealWeightedInnerProduct(const std::shared_ptr<OpType> &weight_op_)
    : weight_op(weight_op_)
  {
    MFEM_VERIFY(weight_op->Height() == weight_op->Width(),
                "Real weight operator must be square! (" << weight_op->Height()
                                                         << " != " << weight_op->Width());
  }
  // Follow same conventions as Dot:  yᴴ x or yᵀ x (note y comes second in the arguments).
  double operator()(const Vector &x, const Vector &y) const
  {
    SetWorkspace(x);
    weight_op->Mult(x, v_workspace);
    return linalg::LocalDot(v_workspace, y);
  }

  std::complex<double> operator()(const ComplexVector &x, const ComplexVector &y) const
  {
    using namespace std::complex_literals;
    SetWorkspace(x.Real());

    // weight_op is real.
    weight_op->Mult(x.Real(), v_workspace);
    std::complex<double> dot = {+linalg::LocalDot(v_workspace, y.Real()),
                                -linalg::LocalDot(v_workspace, y.Imag())};

    weight_op->Mult(x.Imag(), v_workspace);
    dot += std::complex<double>{linalg::LocalDot(v_workspace, y.Imag()),
                                linalg::LocalDot(v_workspace, y.Real())};

    return dot;
  }
};

TEST_CASE("OrthogonalizeColumn - Real Empty", "[orthog][Serial][Parallel][GPU]")
{
  auto orthog_type =
      GENERATE(Orthogonalization::MGS, Orthogonalization::CGS, Orthogonalization::CGS2);

  int mpi_rank = Mpi::Rank(Mpi::World());

  std::vector<Vector> V{};

  Vector w(4);
  w.UseDevice(true);
  auto d_v = w.Write();
  mfem::forall(w.Size(), [=] MFEM_HOST_DEVICE(int i) { d_v[i] = mpi_rank + i; });

  Vector w_orig = w;  // Copy for comparison

  double H[1];  // Not used when m = 0
  linalg::OrthogonalizeColumn(orthog_type, Mpi::World(), V, w, H, 0);

  // Vector should remain unchanged
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

  linalg::OrthogonalizeColumn(orthog_type, Mpi::World(), V, w, H.data(), mpi_size);

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
  std::vector<Vector> V(2);
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
  linalg::OrthogonalizeColumn(orthog_type, Mpi::World(), V, V[1], H.data(), 1);

  // Should be exact in double as multiply by zero. OrthogonalizeColumn does not
  // normalize.
  CHECK_THAT(V[1], RangeEquals(std::vector<double>{0.0, 1.0, 0.0, 0.0}));

  // Normalize by hand.
  auto v1_norm = linalg::Norml2(Mpi::World(), V[1]);
  CHECK_THAT(v1_norm, WithinRel(std::sqrt(mpi_size)) || WithinAbs(0.0, 1e-14));
  V[1] *= 1 / v1_norm;

  // Vector 3: Non-trivial vector to orthogonalize
  Vector w(4);
  w.UseDevice(true);
  auto d_v = w.Write();
  mfem::forall(w.Size(), [=] MFEM_HOST_DEVICE(int i) { d_v[i] = mpi_rank + i; });

  linalg::OrthogonalizeColumn(orthog_type, Mpi::World(), V, w, H.data(), 2);

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

  linalg::OrthogonalizeColumn(orthog_type, Mpi::World(), V, w, H.data(), mpi_size);

  // We should have zeroed out the value at "mpi_rank"
  CHECK_THAT(w.Real()[mpi_rank], WithinAbs(0.0, 1e-12));
  CHECK_THAT(w.Imag()[mpi_rank], WithinAbs(0.0, 1e-12));

  // Check full orthogonalization
  for (int i = 0; i < mpi_size; i++)
  {
    auto dot = linalg::Dot(Mpi::World(), w, V.at(i));
    CHECK_THAT(std::abs(dot), WithinAbs(0.0, 1e-12));
  }
}

TEST_CASE("OrthogonalizeColumn Weighted - Real 1", "[orthog][Serial]")
{
  auto orthog_type =
      GENERATE(Orthogonalization::MGS, Orthogonalization::CGS, Orthogonalization::CGS2);

  mfem::DenseMatrix W(3, 3);
  W = 0.0;
  W(0, 0) = 2.0;
  W(1, 0) = 2.0;
  W(0, 1) = 2.0;
  W(1, 1) = 1.0;
  W(2, 2) = 2.0;

  // Create basis vectors
  std::vector<Vector> V(2);
  for (auto &v : V)
  {
    v.UseDevice(true);
    v.SetSize(3);
    v = 0.0;
  }
  V[0][0] = 1.0 / std::numbers::sqrt2;  // Normalized w.r.t W
  V[1][2] = 1.0 / std::numbers::sqrt2;

  Vector w(3);
  w.UseDevice(true);
  // Fixed seed for determinism. A seedless Randomize() uses time(0) (see
  // mfem::Vector::Randomize), making the test input differ every run.
  w.Randomize(314159);

  std::vector<double> H(2, 0.0);

  RealWeightedInnerProduct weight_op{std::make_shared<mfem::DenseMatrix>(W)};
  linalg::OrthogonalizeColumn(orthog_type, Mpi::World(), V, w, H.data(), 2, weight_op);

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
  auto orthog_type =
      GENERATE(Orthogonalization::MGS, Orthogonalization::CGS, Orthogonalization::CGS2);

  mfem::DenseMatrix W(3, 3);
  W = 0.0;
  W(0, 0) = 2.0;
  W(1, 0) = 2.0;
  W(0, 1) = 2.0;
  W(1, 1) = 1.0;
  W(2, 2) = 2.0;

  // Create basis vectors
  std::vector<ComplexVector> V(2);
  for (auto &v : V)
  {
    v.UseDevice(true);
    v.SetSize(3);
    v = 0.0;
  }
  V[0].Real()[0] = 1.0 / std::numbers::sqrt2;
  V[1].Imag()[2] = 1.0 / std::numbers::sqrt2;

  ComplexVector w(3);
  w.UseDevice(true);
  // Fixed seeds for determinism (see note in the real-valued test above). Real and
  // imaginary parts use distinct seeds so they are not identical sequences.
  w.Real().Randomize(314159);
  w.Imag().Randomize(314160);

  std::vector<std::complex<double>> H(2, 0.0);

  RealWeightedInnerProduct weight_op{std::make_shared<mfem::DenseMatrix>(W)};
  linalg::OrthogonalizeColumn(orthog_type, Mpi::World(), V, w, H.data(), 2, weight_op);

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

TEMPLATE_TEST_CASE("Fused Gram-Schmidt kernels", "[orthog][Serial][Parallel]", Vector,
                   ComplexVector)
{
  using VecType = TestType;
  using ScalarType = std::conditional_t<std::is_same_v<VecType, ComplexVector>,
                                        std::complex<double>, double>;
  const int n = GENERATE(1, 37, 2 * 8192 + 5);
  const int m = GENERATE(0, 1, 5);
  CAPTURE(n, m);
  MPI_Comm comm = Mpi::World();
  std::vector<VecType> V(m + 2);
  for (int j = 0; j < m + 2; j++)
  {
    V[j].SetSize(n);
    V[j].UseDevice(true);
    linalg::SetRandom(comm, V[j], j + 1);
  }
  VecType &x = V[m], &y = V[m + 1];
  auto close = [](ScalarType a, ScalarType b, double scale)
  { return std::abs(a - b) <= 1.0e-13 * scale; };
  const double nx = std::sqrt(std::abs(linalg::LocalDot(x, x)));
  const double ny = std::sqrt(std::abs(linalg::LocalDot(y, y)));

  // Inner products with and without y.
  std::vector<ScalarType> hx(m), hy(m), e(3), hx1(m), e1(3);
  linalg::LocalMultiDot(V, m, x, &y, hx.data(), hy.data(), e.data());
  linalg::LocalMultiDot(V, m, x, nullptr, hx1.data(), nullptr, e1.data());
  for (int j = 0; j < m; j++)
  {
    const double nv = std::sqrt(std::abs(linalg::LocalDot(V[j], V[j])));
    CHECK(close(hx[j], linalg::LocalDot(x, V[j]), nx * nv));
    CHECK(close(hy[j], linalg::LocalDot(y, V[j]), ny * nv));
    CHECK(close(hx1[j], hx[j], nx * nv));
  }
  CHECK(close(e[0], linalg::LocalDot(x, x), nx * nx));
  CHECK(close(e[1], linalg::LocalDot(y, x), nx * ny));
  CHECK(close(e[2], linalg::LocalDot(y, y), ny * ny));
  CHECK(close(e1[0], e[0], nx * nx));

  // Updates, against the vector operations.
  ScalarType g = 0.3;
  if constexpr (std::is_same_v<ScalarType, std::complex<double>>)
  {
    g = {0.3, -0.7};
  }
  // Size of the terms of the updates, for the rounding errors of different evaluation
  // orders.
  double sx = nx, sy = ny;
  for (int j = 0; j < m; j++)
  {
    const double nv = std::sqrt(std::abs(linalg::LocalDot(V[j], V[j])));
    sx += std::abs(hx[j]) * nv;
    sy += std::abs(hy[j]) * nv;
  }
  sy += std::abs(g) * sx;
  for (int variant = 0; variant < 3; variant++)
  {
    CAPTURE(variant);
    VecType x0 = x, y0 = y, x1 = x, y1 = y;
    double nrm2 = 0.0, ref = 0.0;
    if (variant == 0)
    {
      // x = sx (x - V hx).
      nrm2 = linalg::MultiSubtract(V, m, x0, hx.data(), 0.5);
      for (int j = 0; j < m; j++)
      {
        x1.Add(-hx[j], V[j]);
      }
      x1 *= 0.5;
      ref = std::real(linalg::LocalDot(x1, x1));
    }
    else
    {
      // x = sx (x - V hx) (or x = sx x) and y = sy (y - V hy - g x).
      nrm2 = linalg::MultiSubtract(V, m, x0, (variant == 1) ? hx.data() : nullptr, 0.5, &y0,
                                   hy.data(), g, 2.0);
      for (int j = 0; variant == 1 && j < m; j++)
      {
        x1.Add(-hx[j], V[j]);
      }
      x1 *= 0.5;
      for (int j = 0; j < m; j++)
      {
        y1.Add(-hy[j], V[j]);
      }
      y1.Add(-g, x1);
      y1 *= 2.0;
      ref = std::real(linalg::LocalDot(y1, y1));
    }
    x0 -= x1;
    y0 -= y1;
    CHECK(std::sqrt(std::abs(linalg::LocalDot(x0, x0))) <= 1.0e-14 * (1.0 + m) * sx);
    CHECK(std::sqrt(std::abs(linalg::LocalDot(y0, y0))) <= 1.0e-14 * (1.0 + m) * sy);
    CHECK_THAT(nrm2, WithinRel(ref, 1.0e-12) || WithinAbs(ref, 1.0e-300));
  }
}
