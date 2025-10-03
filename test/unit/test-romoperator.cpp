// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include <iterator>
#include <vector>
#include <fmt/format.h>
#include <fmt/ostream.h>
#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/benchmark/catch_benchmark_all.hpp>
#include <catch2/generators/catch_generators_all.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include "models/romoperator.hpp"

using namespace palace;

TEST_CASE("MinimalRationalInterpolation", "[romoperator][Serial][Parallel]")
{
  MPI_Comm comm = Mpi::World();

  auto fn_tan_shift = [](double z)
  { return std::tan(0.5 * M_PI * (z - std::complex<double>(1., 1.))); };

  // Test scalar case: 2 sample points for 4 x 2 vector
  MinimalRationalInterpolation mri_1(6);

  CHECK(mri_1.GetSamplePoints() == std::vector<double>{});
  CHECK_THROWS(mri_1.FindMaxError(1));

  for (double x_sample : {-2.0, -1.0, 1.0, 2.0})
  {
    auto tan_eval = fn_tan_shift(x_sample) / double(Mpi::Size(comm));
    std::vector<std::complex<double>> tmp = {x_sample * tan_eval, -tan_eval, 5. * tan_eval,
                                             10. * x_sample * tan_eval,
                                             2 * x_sample * x_sample * x_sample * tan_eval};
    ComplexVector c_vec(tmp.size());
    c_vec.Set(tmp.data(), tmp.size(), false);

    mri_1.AddSolutionSample(x_sample, c_vec, comm, Orthogonalization::MGS);
  }

  CHECK(mri_1.GetSamplePoints().size() == 4);
  CHECK(mri_1.GetSamplePoints() == std::vector<double>{-2.0, -1.0, 1.0, 2.0});

  // By symmetry of poles max error should be at zero.
  auto max_err_1 = mri_1.FindMaxError(5);
  REQUIRE(max_err_1.size() == 5);

  // By symmetry highest error should be at zero.
  CHECK_THAT(max_err_1[0], Catch::Matchers::WithinAbsMatcher(0.0, 1e-6));

  // Test that elements of max_error are unique.
  // TODO: get better test for multiple N.
  std::sort(max_err_1.begin(), max_err_1.end());
  CHECK(std::adjacent_find(max_err_1.begin(), max_err_1.end()) == max_err_1.end());

  // TODO: Add more stringent tests of MRI, including estimating poles.
}