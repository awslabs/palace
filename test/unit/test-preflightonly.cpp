// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "fixtures.hpp"
#include "surfaceresponse-fixtures.hpp"

#include <filesystem>
#include <fstream>
#include <iterator>
#include <memory>
#include <string>
#include <vector>
#include <catch2/catch_test_macros.hpp>
#include <nlohmann/json.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>
#include "driver.hpp"
#include "fem/mesh.hpp"
#include "models/surfaceresponseoperator.hpp"
#include "utils/communication.hpp"
#include "utils/configfile.hpp"
#include "utils/iodata.hpp"
#include "utils/outputdir.hpp"

using namespace palace;
using namespace Catch::Matchers;
using json = nlohmann::json;
namespace fs = std::filesystem;

namespace
{

std::string ReadBytes(const fs::path &path)
{
  std::ifstream input(path, std::ios::binary);
  return {std::istreambuf_iterator<char>(input), std::istreambuf_iterator<char>()};
}

}  // namespace

// The configuration key PreflightOnly (tracer fixes T1): the response-correction block of
// the configured problem type selects the preflight exactly as --surface-response-preflight
// does; a key in the block of another problem type, or absent, selects nothing.
TEST_CASE_METHOD(test::SurfaceResponseFiles,
                 "PreflightOnly selects the surface-response preflight from the block of "
                 "the configured problem type",
                 "[surfaceresponseoperator][preflightonly][Serial]")
{
  json electrostatic = ConvexIslandConfig();
  json maxwell = ConvexMaxwellIslandConfig();
  CHECK(SurfaceResponsePreflightKey(IoData(electrostatic, false)).empty());
  CHECK(SurfaceResponsePreflightKey(IoData(maxwell, false)).empty());

  json electrostatic_key = electrostatic;
  electrostatic_key["Solver"]["Electrostatic"]["ResponseCorrection"]["PreflightOnly"] =
      true;
  CHECK(SurfaceResponsePreflightKey(IoData(electrostatic_key, false)) ==
        "Solver.Electrostatic.ResponseCorrection.PreflightOnly");
  json maxwell_key = maxwell;
  maxwell_key["Solver"]["SurfaceResponseCorrection"]["PreflightOnly"] = true;
  CHECK(SurfaceResponsePreflightKey(IoData(maxwell_key, false)) ==
        "Solver.SurfaceResponseCorrection.PreflightOnly");
  json maxwell_false = maxwell;
  maxwell_false["Solver"]["SurfaceResponseCorrection"]["PreflightOnly"] = false;
  CHECK(SurfaceResponsePreflightKey(IoData(maxwell_false, false)).empty());

  // The key in the block another problem type reads is ignored like every key of an unused
  // block (an electrostatic problem refuses the Maxwell block itself).
  json maxwell_other = maxwell;
  maxwell_other["Solver"]["Electrostatic"] = electrostatic_key["Solver"]["Electrostatic"];
  CHECK(SurfaceResponsePreflightKey(IoData(maxwell_other, false)).empty());

  json invalid = electrostatic;
  invalid["Solver"]["Electrostatic"]["ResponseCorrection"]["PreflightOnly"] = "Yes";
  CHECK_THROWS_WITH(IoData(invalid, false),
                    ContainsSubstring("\"PreflightOnly\" must be a boolean"));
}

// The key changes nothing inside the preflight: the same configuration with and without it
// writes byte-identical manifests and patch dry runs.
TEST_CASE_METHOD(test::SurfaceResponseFiles,
                 "PreflightOnly leaves the preflight's manifest and patch dry run "
                 "byte-identical",
                 "[surfaceresponseoperator][preflightonly][3d][Serial][Parallel]")
{
#if !defined(MFEM_USE_GSLIB)
  SKIP("SurfaceResponseOperator requires MFEM_USE_GSLIB");
#else
  MPI_Comm comm = Mpi::World();
  auto Preflight = [&](bool key, const fs::path &dir)
  {
    json config = ConvexIslandConfig();
    config["Problem"]["Output"] = dir.string();
    if (key)
    {
      config["Solver"]["Electrostatic"]["ResponseCorrection"]["PreflightOnly"] = true;
    }
    IoData iodata(config, false);
    iodata.boundaries.cracked_attributes.insert(9);
    MakeOutputFolder(iodata, comm);
    Mesh mesh(MakeIslandMesh());
    const auto outcome = WriteSurfaceResponseRequirements(
        iodata, mesh, (dir / "surface-response-requirements.json").string());
    Mpi::Barrier(comm);
    CHECK(outcome.fail_closed_summary.empty());
  };
  const fs::path without = temp.temp_dir / "preflight-without-key";
  const fs::path with = temp.temp_dir / "preflight-with-key";
  Preflight(false, without);
  Preflight(true, with);
  // The patch dry run byte for byte; the manifest equal but for its one run-metadata leaf
  // (the containment test's wall time).
  REQUIRE(fs::is_regular_file(without / "surface-response-patches.csv"));
  CHECK(ReadBytes(without / "surface-response-patches.csv") ==
        ReadBytes(with / "surface-response-patches.csv"));
  auto Manifest = [](const fs::path &dir)
  {
    json manifest = json::parse(ReadBytes(dir / "surface-response-requirements.json"));
    auto &diagnostics = manifest["Identification"]["Diagnostics"];
    if (diagnostics.contains("DomainBoundaryExclusions"))
    {
      diagnostics["DomainBoundaryExclusions"].erase("WallTime");
    }
    return manifest;
  };
  CHECK(Manifest(without) == Manifest(with));
#endif
}
