// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include <catch2/catch_test_macros.hpp>
#include "drivers/eigensolver.hpp"

using palace::internal::MaximumWeightModeAssignment;

TEST_CASE("Eigenmode maximum-weight assignment", "[eigensolver][surface-response][Serial]")
{
  SECTION("Empty")
  {
    CHECK(MaximumWeightModeAssignment({}).empty());
  }

  SECTION("Global optimum differs from row-greedy matching")
  {
    const auto assignment = MaximumWeightModeAssignment({{0.90, 0.80}, {0.85, 0.10}});
    REQUIRE(assignment.size() == 2);
    CHECK(assignment[0] == 1);
    CHECK(assignment[1] == 0);
  }

  SECTION("Rectangular")
  {
    const auto assignment =
        MaximumWeightModeAssignment({{0.20, 0.95, 0.10}, {0.80, 0.10, 0.70}});
    REQUIRE(assignment.size() == 2);
    CHECK(assignment[0] == 1);
    CHECK(assignment[1] == 0);
  }

  SECTION("Invalid")
  {
    CHECK_THROWS(MaximumWeightModeAssignment({{1.0}, {0.5}}));
    CHECK_THROWS(MaximumWeightModeAssignment({{1.0, 0.0}, {0.5}}));
    CHECK_THROWS(MaximumWeightModeAssignment({{1.0, -0.1}}));
  }
}
