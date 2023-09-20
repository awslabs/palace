// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include <mfem.hpp>
#include <catch2/catch_test_macros.hpp>
#include "fem/bilinearform.hpp"
#include "fem/integrator.hpp"

namespace palace
{

namespace
{

void TestCeedOperator()
{

  // XX TODO WIP

  REQUIRE(true);
};

}  // namespace

TEST_CASE("libCEED Operators", "[libCEED]")
{

  // XX TODO MESH DIRECTORY?

  TestCeedOperator();

  // //XX TODO BENCHMARK
  // BENCHMARK("   ")
  // {

  //   return X;
  // }
};

}  // namespace palace
