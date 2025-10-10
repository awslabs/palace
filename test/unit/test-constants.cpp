// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include <sstream>
#include <string>
#include <fmt/format.h>
#include <catch2/catch_test_macros.hpp>
#include <nlohmann/json.hpp>
#include <catch2/benchmark/catch_benchmark_all.hpp>
#include <catch2/generators/catch_generators_all.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <catch2/matchers/catch_matchers_vector.hpp>
#include "utils/configfile.hpp"
#include "utils/iodata.hpp"

using namespace palace;
using namespace palace::electromagnetics;

TEST_CASE("EM Constant Check", "[units][Serial]")
{
  // This test Should be EXACT within double, so no matcher
  CHECK(epsilon0_ * mu0_ * c0_ * c0_ == 1.0);

  CHECK_THAT(epsilon0_, Catch::Matchers::WithinAbs(8.8541878188e-12, 14e-24));
  CHECK_THAT(Z0_, Catch::Matchers::WithinAbs(376.730313412, 59e-12));
  CHECK_THAT(mu0_, Catch::Matchers::WithinRel(4e-7 * M_PI, 1.6e-10));
}