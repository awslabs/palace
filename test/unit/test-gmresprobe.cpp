// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include <cmath>
#include <functional>
#include <limits>
#include <utility>
#include <vector>
#include <catch2/catch_test_macros.hpp>
#include "linalg/gmresprobe.hpp"
#include "linalg/vector.hpp"
#include "utils/communication.hpp"

namespace palace
{

namespace
{

using Apply = std::function<void(const ComplexVector &, ComplexVector &)>;

constexpr int n_local = 64;

// Diagonal operator with entries from lo to hi over the global index.
Apply Diagonal(MPI_Comm comm, double lo, double hi)
{
  return [=](const ComplexVector &x, ComplexVector &y)
  {
    const int n = x.Size(), N = n * Mpi::Size(comm), offset = n * Mpi::Rank(comm);
    const double *xr = x.Real().HostRead(), *xi = x.Imag().HostRead();
    double *yr = y.Real().HostWrite(), *yi = y.Imag().HostWrite();
    for (int i = 0; i < n; i++)
    {
      const double d = lo + (hi - lo) * (offset + i) / (N - 1.0);
      yr[i] = d * xr[i];
      yi[i] = d * xi[i];
    }
  };
}

// Scaled identity, zero and non-finite operators, and the cyclic shift within the block of
// each process.
Apply Scaled(double a)
{
  return [=](const ComplexVector &x, ComplexVector &y)
  {
    y = x;
    y *= a;
  };
}

Apply Shift()
{
  return [](const ComplexVector &x, ComplexVector &y)
  {
    const int n = x.Size();
    const double *xr = x.Real().HostRead(), *xi = x.Imag().HostRead();
    double *yr = y.Real().HostWrite(), *yi = y.Imag().HostWrite();
    for (int i = 0; i < n; i++)
    {
      yr[i] = xr[(i + n - 1) % n];
      yi[i] = xi[(i + n - 1) % n];
    }
  };
}

ComplexVector RandomStart(MPI_Comm comm)
{
  ComplexVector v0(n_local);
  v0.UseDevice(true);
  linalg::SetRandom(comm, v0);
  linalg::Normalize(comm, v0);
  return v0;
}

// The probe with the candidate and reference operators, checking that the selection left
// is the decision and that all processes report the same result.
GmresProbeResult Probe(MPI_Comm comm, const ComplexVector &v0, const Apply &candidate,
                       const Apply &reference, const GmresProbeOptions &opts = {})
{
  bool use_reference = false;
  const auto result = CompareWithGmresProbe<ComplexVector>(
      comm, v0, [&](bool reference) { use_reference = reference; },
      [&](const ComplexVector &x, ComplexVector &y)
      { use_reference ? reference(x, y) : candidate(x, y); }, opts);
  CHECK(use_reference == result.use_reference);
  std::vector<double> v = {static_cast<double>(result.its),
                           result.res_candidate,
                           result.res_reference,
                           result.failed_candidate ? 1.0 : 0.0,
                           result.failed_reference ? 1.0 : 0.0,
                           result.use_reference ? 1.0 : 0.0};
  auto v_min = v, v_max = v;
  Mpi::GlobalMin(static_cast<int>(v.size()), v_min.data(), comm);
  Mpi::GlobalMax(static_cast<int>(v.size()), v_max.data(), comm);
  CHECK(v_min == v_max);
  return result;
}

}  // namespace

TEST_CASE("GMRES probe runs", "[gmresprobe][Serial][Parallel]")
{
  MPI_Comm comm = Mpi::World();
  const auto v0 = RandomStart(comm);
  auto Run = [&](const Apply &apply, int its)
  {
    ProbeGmres<ComplexVector> run(comm, v0, 20);
    std::vector<double> res;
    for (int k = 0; k < its && !run.Stopped(); k++)
    {
      run.Step(apply);
      res.push_back(run.Residual());
    }
    return std::make_pair(std::move(run), res);
  };
  SECTION("Invariant Krylov space")
  {
    auto [run, res] = Run(Scaled(2.0), 20);
    CHECK(run.Iterations() == 1);
    CHECK(run.Stopped());
    CHECK(!run.Failed());
    CHECK(run.Residual() < 1.0e-12);
  }
  SECTION("Start vector in the null space")
  {
    auto [run, res] = Run(Scaled(0.0), 20);
    CHECK(run.Iterations() == 1);
    CHECK(run.Stopped());
    CHECK(!run.Failed());
    CHECK(run.Residual() == 1.0);
  }
  SECTION("Non-finite values")
  {
    auto [run, res] = Run(Scaled(std::numeric_limits<double>::quiet_NaN()), 20);
    CHECK(run.Iterations() == 1);
    CHECK(run.Failed());
    CHECK(run.Stopped());
  }
  SECTION("Distinct eigenvalues")
  {
    auto [run, res] = Run(Diagonal(comm, 1.0, 100.0), 20);
    CHECK(run.Iterations() == 20);
    CHECK(!run.Stopped());
    for (std::size_t k = 1; k < res.size(); k++)
    {
      CHECK(res[k] <= res[k - 1]);
    }
    CHECK(res.back() < 0.1);
  }
}

TEST_CASE("GMRES probe decisions", "[gmresprobe][Serial][Parallel]")
{
  MPI_Comm comm = Mpi::World();
  const auto v0 = RandomStart(comm);
  const GmresProbeOptions opts;
  SECTION("Faster candidate")
  {
    // The candidate converges much faster: kept, decided after the minimum iterations.
    const auto result =
        Probe(comm, v0, Diagonal(comm, 1.0, 2.0), Diagonal(comm, 1.0, 100.0));
    CHECK(!result.use_reference);
    CHECK(result.its >= opts.min_keep_its);
    CHECK(result.its < opts.max_its);
    CHECK(result.ratio < 1.0);
  }
  SECTION("Slower candidate")
  {
    const auto result =
        Probe(comm, v0, Diagonal(comm, 1.0, 100.0), Diagonal(comm, 1.0, 2.0));
    CHECK(result.use_reference);
    CHECK(result.its < opts.max_its);
    CHECK(result.ratio > opts.max_ratio);
  }
  SECTION("Failed candidate")
  {
    const auto result = Probe(comm, v0, Scaled(std::numeric_limits<double>::quiet_NaN()),
                              Diagonal(comm, 1.0, 100.0));
    CHECK(result.failed_candidate);
    CHECK(!result.failed_reference);
    CHECK(result.use_reference);
    CHECK(result.its == 1);
  }
  SECTION("Failed reference")
  {
    const auto result = Probe(comm, v0, Diagonal(comm, 1.0, 100.0),
                              Scaled(std::numeric_limits<double>::quiet_NaN()));
    CHECK(!result.failed_candidate);
    CHECK(result.failed_reference);
    CHECK(!result.use_reference);
  }
  SECTION("Both failed")
  {
    const auto nan = Scaled(std::numeric_limits<double>::quiet_NaN());
    const auto result = Probe(comm, v0, nan, nan);
    CHECK(result.failed_candidate);
    CHECK(result.failed_reference);
    CHECK(result.use_reference);
  }
  SECTION("Start vector in the null space of both")
  {
    // No information: both runs stop at once, and the candidate is kept.
    const auto result = Probe(comm, v0, Scaled(0.0), Scaled(0.0));
    CHECK(result.its == 1);
    CHECK(result.res_candidate == 1.0);
    CHECK(result.res_reference == 1.0);
    CHECK(!result.use_reference);
  }
  SECTION("Total stagnation")
  {
    // GMRES makes no progress for the cyclic shift from the first unit vector of each
    // block: the probe stops after max_its iterations, not twice as many.
    ComplexVector e(n_local);
    e = 0.0;
    e.Real()[0] = 1.0;
    linalg::Normalize(comm, e);
    const auto result = Probe(comm, e, Shift(), Shift());
    CHECK(result.its == opts.max_its);
    CHECK(result.res_candidate > 1.0 - 1.0e-12);
    CHECK(!result.use_reference);
  }
  SECTION("Candidate differing on one process")
  {
    // As smoothing orders raised on some processes only: the candidate is faster on the
    // block of the last process only. All processes run the same iterations and take the
    // same decision (checked in Probe).
    const int last = Mpi::Size(comm) - 1;
    const auto fast = Diagonal(comm, 1.0, 2.0), slow = Diagonal(comm, 1.0, 100.0);
    Apply candidate = [&](const ComplexVector &x, ComplexVector &y)
    { (Mpi::Rank(comm) == last) ? fast(x, y) : slow(x, y); };
    Probe(comm, v0, candidate, slow);
  }
  SECTION("Candidate failing on one process")
  {
    const auto nan = Scaled(std::numeric_limits<double>::quiet_NaN());
    const auto slow = Diagonal(comm, 1.0, 100.0);
    Apply candidate = [&](const ComplexVector &x, ComplexVector &y)
    { (Mpi::Rank(comm) == 0) ? nan(x, y) : slow(x, y); };
    const auto result = Probe(comm, v0, candidate, slow);
    CHECK(result.failed_candidate);
    CHECK(result.use_reference);
  }
}

TEST_CASE("GMRES probe decision rule", "[gmresprobe][Serial]")
{
  // The reference is used when the candidate is predicted to need more than max_ratio times
  // its iterations, ln(r_r) / ln(r_c) > max_ratio, or when the candidate failed.
  constexpr double q = 1.25;
  CHECK(UseReferencePreconditioner(0.5, 0.4, false, false, q));
  CHECK(!UseReferencePreconditioner(0.5, 0.45, false, false, q));
  CHECK(!UseReferencePreconditioner(0.0, 0.0, false, false, q));
  CHECK(UseReferencePreconditioner(1.0, 0.99, false, false, q));
  CHECK(!UseReferencePreconditioner(1.0, 1.0, false, false, q));
  CHECK(UseReferencePreconditioner(0.1, 0.5, true, false, q));
  CHECK(!UseReferencePreconditioner(0.5, 0.1, false, true, q));
  CHECK(UseReferencePreconditioner(0.5, 0.1, true, true, q));
}

}  // namespace palace
