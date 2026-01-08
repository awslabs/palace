// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include <complex>
#include <memory>
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

// Wapper class to make iteration over orthogonalization methods easy.
class orthogonalize_wrapper
{
public:
  Orthogonalization orthgo_type;

  orthogonalize_wrapper(Orthogonalization orthgo_type_) : orthgo_type(orthgo_type_) {}

  template <typename VecType, typename ScalarType,
            typename InnerProductW = linalg::InnerProductStandard>
  void operator()(MPI_Comm comm, const std::vector<VecType> &V, VecType &w, ScalarType *H,
                  std::size_t m, const InnerProductW &dot_op = {}) const
  {
    switch (orthgo_type)
    {
      case Orthogonalization::MGS:
        linalg::OrthogonalizeColumnMGS(comm, V, w, H, m, dot_op);
        break;
      case Orthogonalization::CGS:
        linalg::OrthogonalizeColumnCGS(comm, V, w, H, m, false, dot_op);
        break;
      case Orthogonalization::CGS2:
        linalg::OrthogonalizeColumnCGS(comm, V, w, H, m, true, dot_op);
        break;
    }
  }
};

TEST_CASE("OrthogonalizeColumn - Real Empty", "[orthog][Serial][Parallel][GPU]")
{
  auto orthogonalize_fn = GENERATE(orthogonalize_wrapper(Orthogonalization::MGS),
                                   orthogonalize_wrapper(Orthogonalization::CGS),
                                   orthogonalize_wrapper(Orthogonalization::CGS2));

  int mpi_rank = Mpi::Rank(Mpi::World());

  std::vector<Vector> V{};

  Vector w(4);
  w.UseDevice(true);
  auto d_v = w.Write();
  mfem::forall(w.Size(), [=] MFEM_HOST_DEVICE(int i) { d_v[i] = mpi_rank + i; });

  Vector w_orig = w;  // Copy for comparison

  double H[1];  // Not used when m = 0
  orthogonalize_fn(Mpi::World(), V, w, H, 0);

  // Vector should remain unchanged
  CHECK_THAT(w, RangeEquals(w_orig));
}

TEST_CASE("OrthogonalizeColumn Parameterized - Real 1", "[orthog][Serial][Parallel]")
{
  auto orthogonalize_fn = GENERATE(orthogonalize_wrapper(Orthogonalization::MGS),
                                   orthogonalize_wrapper(Orthogonalization::CGS),
                                   orthogonalize_wrapper(Orthogonalization::CGS2));

  int mpi_rank = Mpi::Rank(Mpi::World());
  int mpi_size = Mpi::Size(Mpi::World());

  std::vector<double> H(mpi_size, 0.0);  // Orthogonalization coefficients
  std::vector<Vector> V(mpi_size + 1);
  for (auto &v : V)
  {
    v.UseDevice(true);
    v.SetSize(mpi_size + 1);
    v = 0.0;
  }

  // Rank-wise orthogonal basis vectors
  V.at(mpi_rank)[mpi_rank] = 1.0;

  // Pick random vector
  auto &w = V.at(mpi_size);
  w.Randomize();

  orthogonalize_fn(Mpi::World(), V, w, H.data(), mpi_size);

  // We should have zeroed out the value at "mpi_rank"
  CHECK_THAT(w[mpi_rank], WithinAbs(0.0, 1e-12));

  // Check full orthogonalization
  for (int i = 0; i < mpi_size; i++)
  {
    auto dot = linalg::Dot(Mpi::World(), w, V.at(i));
    CHECK_THAT(dot, WithinAbs(0.0, 1e-12));
  }
}

TEST_CASE("OrthogonalizeColumn Parameterized - Real 2", "[orthog][Serial][Parallel][GPU]")
{
  auto orthogonalize_fn = GENERATE(orthogonalize_wrapper(Orthogonalization::MGS),
                                   orthogonalize_wrapper(Orthogonalization::CGS),
                                   orthogonalize_wrapper(Orthogonalization::CGS2));

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
    v = 0.0;
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
  orthogonalize_fn(Mpi::World(), V, V[1], H.data(), 1);

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

  orthogonalize_fn(Mpi::World(), V, w, H.data(), 2);

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
  auto orthogonalize_fn = GENERATE(orthogonalize_wrapper(Orthogonalization::MGS),
                                   orthogonalize_wrapper(Orthogonalization::CGS),
                                   orthogonalize_wrapper(Orthogonalization::CGS2));

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
  auto theta = (2. * M_PI * mpi_rank) / mpi_size;
  V.at(mpi_rank).Real()[mpi_rank] = std::cos(theta);
  V.at(mpi_rank).Imag()[mpi_rank] = std::sin(theta);

  // Pick random vector to orthogonalize
  auto &w = V.at(mpi_size);
  w.Real().Randomize();
  w.Imag().Randomize();

  orthogonalize_fn(Mpi::World(), V, w, H.data(), mpi_size);

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
  auto orthogonalize_fn = GENERATE(orthogonalize_wrapper(Orthogonalization::MGS),
                                   orthogonalize_wrapper(Orthogonalization::CGS),
                                   orthogonalize_wrapper(Orthogonalization::CGS2));

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
  V[0][0] = 1.0 / std::sqrt(2);  // Normalized w.r.t W
  V[1][2] = 1.0 / std::sqrt(2);

  Vector w(3);
  w.UseDevice(true);
  w.Randomize();

  std::vector<double> H(2, 0.0);

  linalg::InnerProductRealWeight weight_op{std::make_shared<mfem::DenseMatrix>(W)};
  orthogonalize_fn(Mpi::World(), V, w, H.data(), 2, weight_op);

  // Check orthogonality with respect to weight matrix
  Vector WVj(3);  // Temporary workspace
  W.Mult(V[0], WVj);
  double dot0 = linalg::Dot(Mpi::World(), w, WVj);
  W.Mult(V[1], WVj);
  double dot1 = linalg::Dot(Mpi::World(), w, WVj);

  CHECK_THAT(dot0, WithinAbs(0.0, 1e-12));
  CHECK_THAT(dot1, WithinAbs(0.0, 1e-12));

  // Resulting vector is along (1, -1) direction
  CHECK_THAT(w[0], WithinRel(-w[1]));
}

TEST_CASE("OrthogonalizeColumn Weighted - Complex 1", "[orthog][Serial]")
{
  auto orthogonalize_fn = GENERATE(orthogonalize_wrapper(Orthogonalization::MGS),
                                   orthogonalize_wrapper(Orthogonalization::CGS),
                                   orthogonalize_wrapper(Orthogonalization::CGS2));

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
  V[0].Real()[0] = 1.0 / std::sqrt(2);
  V[1].Imag()[2] = 1.0 / std::sqrt(2);

  ComplexVector w(3);
  w.UseDevice(true);
  w.Real().Randomize();
  w.Imag().Randomize();

  std::vector<std::complex<double>> H(2, 0.0);

  linalg::InnerProductRealWeight weight_op{std::make_shared<mfem::DenseMatrix>(W)};
  orthogonalize_fn(Mpi::World(), V, w, H.data(), 2, weight_op);

  auto W_wrap = ComplexWrapperOperator(&W, nullptr);

  // Check orthogonality with respect to weight matrix
  ComplexVector WVj(3);  // Temporary workspace
  W_wrap.Mult(V[0], WVj);
  W_wrap.Mult(V[1], WVj);

  auto dot0 = linalg::Dot(Mpi::World(), w, WVj);
  auto dot1 = linalg::Dot(Mpi::World(), w, WVj);

  CHECK_THAT(dot0.real(), WithinAbs(0.0, 1e-12));
  CHECK_THAT(dot0.imag(), WithinAbs(0.0, 1e-12));
  CHECK_THAT(dot1.real(), WithinAbs(0.0, 1e-12));
  CHECK_THAT(dot1.imag(), WithinAbs(0.0, 1e-12));

  // Resulting vector is along (1, -1) direction
  CHECK_THAT(w.Real()[0], WithinRel(-w.Real()[1]));
  CHECK_THAT(w.Imag()[0], WithinRel(-w.Imag()[1]));
}
