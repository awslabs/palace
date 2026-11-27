// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "fixtures.hpp"
#include "surfaceresponse-fixtures.hpp"

#include <map>
#include <set>
#include <vector>
#include <catch2/catch_test_macros.hpp>
#include <nlohmann/json.hpp>
#include <catch2/generators/catch_generators.hpp>
#include "fem/mesh.hpp"
#include "models/surfaceresponseoperator.hpp"
#include "utils/iodata.hpp"

using namespace palace;
using json = nlohmann::json;

// The substrate check (tracer fixes T3, a warning only): the configured domains next to the
// target interfaces' boundary elements whose relative permittivity is neither vacuum nor
// the library's Fabrication.SubstratePermittivity. The island mesh's metal sheet (attribute
// 9) lies inside domain 1 (both sides).
TEST_CASE_METHOD(test::SurfaceResponseFiles,
                 "FindSubstratePermittivityMismatches names the domains next to the "
                 "target interfaces that are neither vacuum nor the library substrate",
                 "[surfaceresponseoperator][substratepermittivity][3d][Serial][Parallel]")
{
  // P4 S3-0 (review r2084 F-1): both crack forms (the legacy key; upstream's default = the
  // key unset), as P3-23.
  const auto crack_policy = GENERATE(test::CrackPolicy::Legacy, test::CrackPolicy::Default);
  INFO("crack policy: " << test::CrackPolicyName(crack_policy));
  auto Mismatches = [&](const json &permittivity, const std::set<int> &attributes)
  {
    json config = ConvexIslandConfig(crack_policy);
    if (!permittivity.is_null())
    {
      config["Domains"]["Materials"][0]["Permittivity"] = permittivity;
    }
    IoData iodata(config, false);
    Mesh mesh(MakeIslandMesh());
    return FindSubstratePermittivityMismatches(iodata, mesh.Get(), attributes, 11.45);
  };
  // Vacuum (the default 1) and the library substrate are not reported.
  CHECK(Mismatches(json(), {9}).empty());
  CHECK(Mismatches(11.45, {9}).empty());
  CHECK(Mismatches(1.0, {9}).empty());
  // The device's silicon at 11.47 next to the metal is reported, by material and attribute.
  const auto reported = Mismatches(11.47, {9});
  REQUIRE(reported.size() == 1);
  CHECK(reported.at(0) == std::vector<int>{1});
  // An anisotropic domain with one principal value off the substrate is reported.
  CHECK(Mismatches(json::array({11.45, 11.45, 11.6}), {9}).size() == 1);
  // No target boundary element: nothing is read.
  CHECK(Mismatches(11.47, {}).empty());
  // Vacuum is Palace's own test (IsVacuumMaterial, the vacuum side of MA and SA
  // interfaces): a domain of permittivity 1 with a permeability 2 is reported.
  json config = ConvexIslandConfig(crack_policy);
  config["Domains"]["Materials"][0]["Permeability"] = 2.0;
  IoData iodata(config, false);
  Mesh mesh(MakeIslandMesh());
  CHECK(FindSubstratePermittivityMismatches(iodata, mesh.Get(), {9}, 11.45).size() == 1);
}
