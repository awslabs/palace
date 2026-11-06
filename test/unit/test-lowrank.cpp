// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include <algorithm>
#include <cmath>
#include <complex>
#include <memory>
#include <random>
#include <vector>
#include <catch2/catch_template_test_macros.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators_all.hpp>
#include "linalg/lowrank.hpp"
#include "linalg/vector.hpp"
#include "models/floquetportoperator.hpp"
#include "utils/communication.hpp"

namespace palace
{

namespace
{

void FillRandom(ComplexVector &v, std::mt19937 &gen)
{
  std::uniform_real_distribution<double> dist(-1.0, 1.0);
  auto *VR = v.Real().HostWrite();
  auto *VI = v.Imag().HostWrite();
  for (int i = 0; i < v.Size(); i++)
  {
    VR[i] = dist(gen);
    VI[i] = dist(gen);
  }
}

std::vector<std::complex<double>> Values(const ComplexVector &u)
{
  std::vector<std::complex<double>> w(u.Size());
  const auto *UR = u.Real().HostRead();
  const auto *UI = u.Imag().HostRead();
  for (int i = 0; i < u.Size(); i++)
  {
    w[i] = {UR[i], UI[i]};
  }
  return w;
}

// Largest entrywise error, relative to the largest reference entry (at least 1).
double MaxError(const ComplexVector &u, const std::vector<std::complex<double>> &ref)
{
  const auto w = Values(u);
  double err = 0.0, scale = 1.0;
  for (std::size_t i = 0; i < ref.size(); i++)
  {
    err = std::max({err, std::abs(w[i].real() - ref[i].real()),
                    std::abs(w[i].imag() - ref[i].imag())});
    scale = std::max({scale, std::abs(ref[i].real()), std::abs(ref[i].imag())});
  }
  return err / scale;
}

}  // namespace

// Compare the port-packed update against a reference in complex arithmetic (one global
// dot per term, then full-length updates).
TEMPLATE_TEST_CASE("LowRankComplexUpdate matches a complex reference",
                   "[lowrank][Serial][Parallel][GPU]", std::true_type, std::false_type)
{
  constexpr bool conj = TestType::value;
  MPI_Comm comm = Mpi::World();
  const int rank = Mpi::Rank(comm), size = Mpi::Size(comm);
  // One device block per term, uncapped and capped blocks per term, and the cap of one
  // block.
  const auto [n0, K] =
      GENERATE(table<int, int>({{23, 7}, {8200, 7}, {8200, 80}, {40, 1100}}));
  const int n = n0 + 3 * rank;

  // 0: sparse support on every rank, 1: support on rank 0 only, 2: full support, 3: three
  // disjoint ports of 1/6, 1/3 and 1/2 of the entries (term k on port k % 3, every other
  // entry for odd k), 4: as 3, with the last term also on the first entry of the next
  // port, which joins two of the ports, 5: term k only on rank k % size, on every third
  // entry.
  const int support = GENERATE(0, 1, 2, 3, 4, 5);
  std::mt19937 gen(1234 + 17 * rank + 101 * support + n0 + K);
  std::uniform_real_distribution<double> dist(-1.0, 1.0);

  auto port = [n](int i) { return (6 * i < n) ? 0 : (2 * i < n) ? 1 : 2; };
  auto first = [n](int p) { return (p == 0) ? 0 : (p == 1) ? (n + 5) / 6 : (n + 1) / 2; };
  std::vector<std::unique_ptr<ComplexVector>> v;
  std::vector<const ComplexVector *> v_ptr;
  std::vector<std::complex<double>> g;
  for (int k = 0; k < K; k++)
  {
    auto vk = std::make_unique<ComplexVector>(n);
    vk->UseDevice(true);
    auto *VR = vk->Real().HostWrite();
    auto *VI = vk->Imag().HostWrite();
    for (int i = 0; i < n; i++)
    {
      const int b = k % 3;
      const bool on =
          (support == 2) || (support == 0 && (i % 4 == 0 || (i % 4 == 1 && k % 2 == 0))) ||
          (support == 1 && rank == 0 && i < 4 + k) ||
          ((support == 3 || support == 4) && port(i) == b && (k % 2 == 0 || i % 2 == 0)) ||
          (support == 4 && k == K - 1 && i == first((b + 1) % 3)) ||
          (support == 5 && rank == k % size && i % 3 == 0);
      VR[i] = on ? dist(gen) : 0.0;
      VI[i] = on ? dist(gen) : 0.0;
    }
    g.emplace_back(0.3 + 0.1 * (k % 7), -0.2 + 0.05 * (k % 7));
    v_ptr.push_back(vk.get());
    v.push_back(std::move(vk));
  }

  ComplexVector x(n), y0(n);
  x.UseDevice(true);
  y0.UseDevice(true);
  FillRandom(x, gen);
  FillRandom(y0, gen);
  const std::complex<double> a(0.7, -0.4);

  // Reference results for y0 and, with x aliasing y, for x.
  const auto xc = Values(x);
  std::vector<std::vector<std::complex<double>>> vc(K);
  std::vector<double> d(2 * K, 0.0);
  for (int k = 0; k < K; k++)
  {
    vc[k] = Values(*v[k]);
    std::complex<double> dot = 0.0;
    for (int i = 0; i < n; i++)
    {
      dot += (conj ? std::conj(vc[k][i]) : vc[k][i]) * xc[i];
    }
    d[2 * k] = dot.real();
    d[2 * k + 1] = dot.imag();
  }
  Mpi::GlobalSum(2 * K, d.data(), comm);
  auto y_ref = Values(y0), x_ref = xc;
  for (int k = 0; k < K; k++)
  {
    const auto s = a * g[k] * std::complex<double>(d[2 * k], d[2 * k + 1]);
    for (int i = 0; i < n; i++)
    {
      y_ref[i] += s * vc[k][i];
      x_ref[i] += s * vc[k][i];
    }
  }

  LowRankComplexUpdate<conj> update(comm, n, v_ptr, g);
  const bool local = (support != 1 || rank == 0) && (support != 5 || rank < K);
  const int nnz = !local           ? 0
                  : (support == 0) ? (n + 3) / 4 + (n + 2) / 4
                  : (support == 1) ? std::min(n, K + 3)
                  : (support == 5) ? (n + 2) / 3
                                   : n;
  CHECK(update.NumTerms() == K);
  CHECK(update.NumGroups() == (!local ? 0 : (support == 3) ? 3 : (support == 4) ? 2 : 1));
  CHECK(update.NumPackedEntries() == nnz);

  // The copies must be independent of the source vectors.
  for (auto &vk : v)
  {
    *vk = 0.0;
  }

  ComplexVector y(y0);
  update.AddMult(x, y, a);
  CHECK(MaxError(y, y_ref) <= 1e-12);

  update.AddMult(x, x, a);
  CHECK(MaxError(x, x_ref) <= 1e-12);
}

TEST_CASE("LowRankComplexUpdate with no terms", "[lowrank][Serial][Parallel][GPU]")
{
  LowRankComplexUpdate<true> update(Mpi::World(), 5, {}, {});
  ComplexVector x(5), y(5);
  x = 1.0;
  y = 2.0;
  update.AddMult(x, y, 1.0);
  const auto *YR = y.Real().HostRead();
  for (int i = 0; i < 5; i++)
  {
    CHECK(YR[i] == 2.0);
  }
}

// The Floquet port operator uses the conjugated inner product: F v = g (vᴴ v) v = g ‖v‖² v.
TEST_CASE("Floquet low-rank operator uses the conjugated inner product",
          "[lowrank][Serial][Parallel][GPU]")
{
  MPI_Comm comm = Mpi::World();
  const int n = 11 + Mpi::Rank(comm);
  std::mt19937 gen(42 + Mpi::Rank(comm));
  ComplexVector v(n), y(n);
  v.UseDevice(true);
  y.UseDevice(true);
  FillRandom(v, gen);
  const std::complex<double> g(0.3, -0.2);
  LowRankComplexOperator op(comm, n, {&v}, {g});
  op.Mult(v, y);
  const double v_sq = std::pow(linalg::Norml2(comm, v), 2);
  auto ref = Values(v);
  for (auto &r : ref)
  {
    r *= g * v_sq;
  }
  CHECK(MaxError(y, ref) <= 1e-12);
}

}  // namespace palace
