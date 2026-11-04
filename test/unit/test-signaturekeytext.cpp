// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "fixtures.hpp"
#include "surfaceresponse-fixtures.hpp"

#include <cmath>
#include <filesystem>
#include <fstream>
#include <string>
#include <vector>
#include <catch2/catch_test_macros.hpp>
#include <nlohmann/json.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>
#include "fem/mesh.hpp"
#include "models/surfaceresponseidentification.hpp"
#include "models/surfaceresponseoperator.hpp"
#include "utils/communication.hpp"
#include "utils/iodata.hpp"

namespace palace
{

namespace fs = std::filesystem;
using json = nlohmann::json;
using namespace Catch::Matchers;

namespace
{

// The KeyText contract of a Signature-bearing manifest entry (decision 605 (2) F-4 (c)):
// the text whose sha256 is the Hash, and whose parse is the Signature (Type included).
void CheckKeyText(const json &entry)
{
  REQUIRE(entry.contains("KeyText"));
  REQUIRE(entry["KeyText"].is_string());
  const std::string key_text = entry["KeyText"].get<std::string>();
  CHECK(Sha256Hex(key_text) == entry["Hash"].get<std::string>());
  const json parsed = json::parse(key_text);
  CHECK(parsed == entry["Signature"]);
  CHECK(parsed["Type"] == entry["Signature"]["Type"]);
}

}  // namespace

TEST_CASE("SignatureKeyAndHash key text carries the hashed lexemes",
          "[surfaceresponseidentification][keytext][Serial]")
{
  // The three doubles of the stored 7fd482e5adfc Signature whose nlohmann (Grisu2) lexemes
  // are not the shortest round-trip digits (the R2-failure diagnosis F-4, MEASURED):
  // Python's repr writes 7.378134999999999 / -2.953708 / -2.732016, so a float re-dump does
  // not hash to the key. The key TEXT does, by construction.
  const json signature = json::parse(R"({"EdgeCount":2,"Portions":[{"P":[7.378134999999999,
      -2.953708,-2.732016,0.5]}]})");
  const auto [key, hash] = SignatureKeyAndHash(signature, "SpatialEdgeCluster");
  CHECK_THAT(key, ContainsSubstring("7.3781349999999994"));
  CHECK_THAT(key, ContainsSubstring("-2.9537079999999998"));
  CHECK_THAT(key, ContainsSubstring("-2.7320159999999998"));
  CHECK(Sha256Hex(key) == hash);
  // The same fixture pins the Python side (test_signature_library.py KEYTEXT_FIXTURE): the
  // key text and its digest are the cross-language contract.
  CHECK(key == R"({"EdgeCount":2,"Portions":[{"P":[7.3781349999999994,-2.9537079999999998,)"
               R"(-2.7320159999999998,0.5]}],"Type":"SpatialEdgeCluster"})");
  CHECK(hash == "95423956626c93e2982f74f01c7a119a57a8cd44495a054243071caffa1b3779");
  json expected = signature;
  expected["Type"] = "SpatialEdgeCluster";
  CHECK(json::parse(key) == expected);
  // The shortest-digit text (what a float re-dump produces) is a different key.
  const std::string shortest =
      std::string(R"({"EdgeCount":2,"Portions":[{"P":[7.378134999999999,-2.953708,)") +
      R"(-2.732016,0.5]}],"Type":"SpatialEdgeCluster"})";
  CHECK(json::parse(shortest) == expected);
  CHECK(Sha256Hex(shortest) != hash);
  CHECK(Sha256Hex(shortest) ==
        "5e30744232102a28c1997259012c6465418308055fcdb63cb055f9f1fd99da83");
}

// R-9b (fix round P3; decision 1286 Q-5's CONDITION): the family keys the planner writes
// (prepare_surface_response_coupons.py --taper-families;
// test_prepare_surface_response_coupons TaperFamilyPlanner.FE6_FAMILY_KEYS) are what the
// C++ keys at those cross-sections: the FE6 f32f mean section (0 / 0.826394 / 2.203837
// / 3.030321 R) scaled by s and quantised to the 1e-6 R key grid (RoundTo's arithmetic)
// hashes to the same digests through SignatureKeyAndHash (the Grisu2
// lexemes 2.4301019999999998 / 0.7437549999999999 / 2.4242209999999997 / 3.6260459999999997
// included) - the cross-language contract of the family rows.
TEST_CASE("SignatureKeyAndHash R-9b family keys",
          "[surfaceresponseidentification][keytext][r9b][Serial]")
{
  const json mean = json::parse(R"({"Edges":[
      {"Conductor":1,"GapSide":1,"Interfaces":["MA","MS","SA"],"Law":"{\"Type\":\"PEC\"}","OffsetOverR":0.0},
      {"Conductor":1,"GapSide":-1,"Interfaces":["MA","MS","SA"],"Law":"{\"Type\":\"PEC\"}","OffsetOverR":0.826394},
      {"Conductor":1,"GapSide":1,"Interfaces":["MA","MS","SA"],"Law":"{\"Type\":\"PEC\"}","OffsetOverR":2.203837},
      {"Conductor":1,"GapSide":-1,"Interfaces":["MA","MS","SA"],"Law":"{\"Type\":\"PEC\"}","OffsetOverR":3.030321}]})");
  auto Quantize = [](double value)
  {
    const double r = std::round(value / 1.0e-6) * 1.0e-6;
    return r == 0.0 ? 0.0 : r;
  };
  const std::vector<std::pair<double, std::string>> family = {
      {0.801929, "6eb211ccce87ea1a"},
      {0.9, "a5ef2edf3ce909a6"},
      {1.1, "9915bcf3b40d8133"},
      {1.196588, "e339af1d9c33c94d"}};
  for (const auto &[scale, hash16] : family)
  {
    json scaled = mean;
    for (auto &edge : scaled["Edges"])
    {
      edge["OffsetOverR"] = Quantize(edge["OffsetOverR"].get<double>() * scale);
    }
    const auto [key, hash] = SignatureKeyAndHash(scaled, "ParallelEdgeCluster");
    INFO("scale " << scale << " key " << key);
    CHECK(hash.substr(0, 16) == hash16);
    CHECK(Sha256Hex(key) == hash);
    CHECK(json::parse(key)["Type"] == "ParallelEdgeCluster");
  }
  const auto [key, hash] = SignatureKeyAndHash(mean, "ParallelEdgeCluster");
  CHECK_THAT(key, ContainsSubstring("\"OffsetOverR\":0.826394}"));
  CHECK(hash == "f32f48576bfa3bf68f5b1066ebda2405f3b00b6a17bbffa7dec012d36a8aaf2d");
}

TEST_CASE_METHOD(test::SurfaceResponseFiles,
                 "Surface-response requirements manifest writes KeyText beside Hash",
                 "[surfaceresponseoperator][3d][features][keytext][Serial][Parallel]")
{
#if !defined(MFEM_USE_GSLIB)
  SKIP("SurfaceResponseOperator requires MFEM_USE_GSLIB");
#else
  // The features preflight of the island (every feature Missing against the legacy convex
  // library: 4 convex corners + 4 isolated edges, as the features patch construction case
  // asserts): both Signature-bearing writers of the manifest — the version-2 Requirements
  // records and the Identification.Features entries — carry KeyText.
  json island_config = IslandConfig();
  auto &correction = island_config["Solver"]["Electrostatic"]["ResponseCorrection"];
  correction.erase("PatchConstruction");
  correction["Library"] = convex_library_3d_path.string();
  IoData iodata(island_config, false);
  iodata.boundaries.cracked_attributes.insert(9);
  Mesh island_mesh(MakeIslandMesh());
  const auto manifest_path = temp.temp_dir / "surface-response-requirements-keytext.json";
  WriteSurfaceResponseRequirements(iodata, island_mesh, manifest_path.string());
  Mpi::Barrier(Mpi::World());
  std::ifstream manifest_input(manifest_path);
  REQUIRE(manifest_input);
  const json manifest = json::parse(manifest_input);
  REQUIRE(manifest["Version"] == 2);

  const auto &features = manifest["Identification"]["Features"];
  REQUIRE(features.size() == 8);
  for (const auto &feature : features)
  {
    INFO("feature " << feature["Id"] << " " << feature["Type"]);
    CheckKeyText(feature);
    CHECK(feature["Signature"]["Type"] == feature["Type"]);
  }

  const auto &requirements = manifest["Requirements"];
  REQUIRE(!requirements.empty());
  for (const auto &record : requirements)
  {
    INFO("requirement " << record["Topology"] << " " << record["Hash"]);
    REQUIRE(record.contains("Signature"));
    CheckKeyText(record);
    CHECK(record["Signature"]["Type"] == record["Topology"]);
  }
#endif
}

}  // namespace palace
