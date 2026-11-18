// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include <complex>
#include <memory>
#include <random>
#include <vector>
#include <catch2/catch_template_test_macros.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators_all.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include "linalg/lowrank.hpp"
#include "linalg/vector.hpp"
#include "utils/communication.hpp"

namespace palace
{
using namespace Catch::Matchers;

namespace
{

void FillRandom(ComplexVector &v, std::mt19937 &gen, bool keep_zeros = false)
{
  std::uniform_real_distribution<double> dist(-1.0, 1.0);
  auto *VR = v.Real().HostReadWrite();
  auto *VI = v.Imag().HostReadWrite();
  for (int i = 0; i < v.Size(); i++)
  {
    if (keep_zeros && VR[i] == 0.0 && VI[i] == 0.0)
    {
      continue;
    }
    VR[i] = dist(gen);
    VI[i] = dist(gen);
  }
}

}  // namespace

// Compare the batched, port-packed update against the per-term reference it replaces (one
// global dot per real product, then full-length updates).
TEMPLATE_TEST_CASE("LowRankComplexUpdate matches per-term reference",
                   "[lowrank][Serial][Parallel][GPU]", std::true_type, std::false_type)
{
  constexpr bool conj = TestType::value;
  MPI_Comm comm = Mpi::World();
  const int rank = Mpi::Rank(comm);
  const int n = 23 + 3 * rank;
  const int K = 7;

  // 0: sparse support on every rank, 1: support on rank 0 only, 2: full support.
  const int support = GENERATE(0, 1, 2);
  std::mt19937 gen(1234 + 17 * rank + 101 * support);

  std::vector<std::unique_ptr<ComplexVector>> v;
  std::vector<const ComplexVector *> v_ptr;
  std::vector<std::complex<double>> g;
  for (int k = 0; k < K; k++)
  {
    auto vk = std::make_unique<ComplexVector>(n);
    vk->UseDevice(true);
    *vk = 0.0;
    {
      auto *VR = vk->Real().HostReadWrite();
      for (int i = 0; i < n; i++)
      {
        const bool on = (support == 2) ||
                        (support == 0 && (i % 4 == 0 || (i % 4 == 1 && k % 2 == 0))) ||
                        (support == 1 && rank == 0 && i < 4 + k);
        VR[i] = on ? 1.0 : 0.0;
      }
    }
    FillRandom(*vk, gen, /*keep_zeros=*/true);
    g.emplace_back(0.3 + 0.1 * k, -0.2 + 0.05 * k);
    v_ptr.push_back(vk.get());
    v.push_back(std::move(vk));
  }

  ComplexVector x(n), y0(n);
  x.UseDevice(true);
  y0.UseDevice(true);
  FillRandom(x, gen);
  FillRandom(y0, gen);
  const std::complex<double> a(0.7, -0.4);

  ComplexVector y_ref(y0);
  for (int k = 0; k < K; k++)
  {
    const auto &vk = *v[k];
    const double rr = linalg::Dot(comm, vk.Real(), x.Real());
    const double ii = linalg::Dot(comm, vk.Imag(), x.Imag());
    const double ri = linalg::Dot(comm, vk.Real(), x.Imag());
    const double ir = linalg::Dot(comm, vk.Imag(), x.Real());
    const auto dot = conj ? std::complex<double>(rr + ii, ri - ir)
                          : std::complex<double>(rr - ii, ri + ir);
    y_ref.Add(a * g[k] * dot, vk);
  }

  LowRankComplexUpdate<conj> update(comm, v_ptr, g);
  CHECK(update.NumTerms() == K);
  if (support == 2)
  {
    CHECK(update.NumPackedEntries() == n);
  }
  else if (support == 1 && rank > 0)
  {
    CHECK(update.NumPackedEntries() == 0);
  }
  else
  {
    CHECK(update.NumPackedEntries() < n);
  }

  // The copies must be independent of the source vectors.
  for (auto &vk : v)
  {
    *vk = 0.0;
  }

  ComplexVector y(y0);
  update.AddMult(x, y, a);

  const auto *YR = y.Real().HostRead();
  const auto *YI = y.Imag().HostRead();
  const auto *RR = y_ref.Real().HostRead();
  const auto *RI = y_ref.Imag().HostRead();
  for (int i = 0; i < n; i++)
  {
    CHECK_THAT(YR[i], WithinAbs(RR[i], 1e-12));
    CHECK_THAT(YI[i], WithinAbs(RI[i], 1e-12));
  }

  // Entries outside every term's support are left untouched.
  if (support == 1 && rank > 0)
  {
    const auto *Y0R = y0.Real().HostRead();
    const auto *Y0I = y0.Imag().HostRead();
    for (int i = 0; i < n; i++)
    {
      CHECK(YR[i] == Y0R[i]);
      CHECK(YI[i] == Y0I[i]);
    }
  }
}

TEST_CASE("LowRankComplexUpdate with no terms", "[lowrank][Serial][Parallel][GPU]")
{
  LowRankComplexUpdate<true> update(Mpi::World(), {}, {});
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

}  // namespace palace
