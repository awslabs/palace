// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "fixtures.hpp"
#include "surfaceresponse-fixtures.hpp"

#include <cstdio>
#include <filesystem>
#include <fstream>
#include <iterator>
#include <string>
#include <unistd.h>
#include <catch2/catch_test_macros.hpp>
#include <nlohmann/json.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>
#include "fem/mesh.hpp"
#include "models/surfaceresponseoperator.hpp"
#include "utils/communication.hpp"
#include "utils/iodata.hpp"
#include "utils/outputdir.hpp"

using namespace palace;
using namespace Catch::Matchers;
using json = nlohmann::json;
namespace fs = std::filesystem;

namespace
{

json ReadJson(const fs::path &path)
{
  std::ifstream input(path);
  REQUIRE(input);
  json data;
  input >> data;
  return data;
}

// The text a callable writes to the process's standard output (the root's log).
template <typename F>
std::string CaptureStdout(const fs::path &file, F &&run)
{
  std::fflush(stdout);
  const int saved = dup(fileno(stdout));
  std::FILE *target = std::fopen(file.c_str(), "w");
  REQUIRE(target);
  dup2(fileno(target), fileno(stdout));
  run();
  std::fflush(stdout);
  dup2(saved, fileno(stdout));
  close(saved);
  std::fclose(target);
  std::ifstream input(file);
  return {std::istreambuf_iterator<char>(input), std::istreambuf_iterator<char>()};
}

}  // namespace

// The Missing warning names the metal law when it is the cause (tracer fixes T4; visibility
// only): one extra line, only when unmatched features have a library model of the same
// geometry class and parameters under another metal boundary law.
TEST_CASE_METHOD(test::SurfaceResponseFiles,
                 "The Missing warning names the metal boundary law only when it is why "
                 "features have no model",
                 "[surfaceresponseoperator][lawcausedmissing][3d][Serial]")
{
#if !defined(MFEM_USE_GSLIB)
  SKIP("SurfaceResponseOperator requires MFEM_USE_GSLIB");
#else
  MPI_Comm comm = Mpi::World();
  // The convex library keyed on SA only like the island's features (all models perfect
  // conductors).
  const fs::path sa_library = temp.temp_dir / "fabrication-process-convex-3d-sa-law.json";
  {
    json library = ReadJson(convex_library_3d_path);
    library["Name"] = "unit-test-process-convex-3d-sa";
    for (auto &model : library["Models"])
    {
      model["Interfaces"] = {{{"Type", "SA"}, {"Coupon", 1}}};
    }
    std::ofstream output(sa_library);
    output << library.dump(2) << "\n";
  }
  // The Maxwell island: the outer box PEC, the metal sheet (attribute 9) PEC or impedance.
  auto Log = [&](bool impedance, const fs::path &library, const std::string &name)
  {
    json config = ConvexMaxwellIslandConfig();
    const fs::path dir = temp.temp_dir / name;
    config["Problem"]["Output"] = dir.string();
    auto &correction = config["Solver"]["SurfaceResponseCorrection"];
    correction["Library"] = library.string();
    correction["UnmatchedPolicy"] = "Warn";
    correction.erase("PatchConstruction");
    if (impedance)
    {
      config["Boundaries"]["Ground"]["Attributes"] = {1, 2, 3, 4, 5, 6};
      config["Boundaries"]["Impedance"] = {{{"Attributes", {9}}, {"Ls", 4.5e-14}}};
    }
    IoData iodata(config, false);
    iodata.boundaries.cracked_attributes.insert(9);
    MakeOutputFolder(iodata, comm);
    Mesh mesh(MakeIslandMesh());
    return CaptureStdout(temp.temp_dir / (name + ".log"),
                         [&]()
                         {
                           WriteSurfaceResponseRequirements(
                               iodata, mesh,
                               (dir / "surface-response-requirements.json").string());
                         });
  };
  const std::string law_line = "have library models of the same geometry class under "
                               "another metal boundary law";

  // All perfect conductors: no law line (the features match or are Missing for geometry).
  const std::string pec = Log(false, sa_library, "law-pec");
  CHECK_THAT(pec, !ContainsSubstring(law_line));

  // Impedance metal, perfect-conductor library: every unmatched feature the same library
  // serves under PEC is named with the device law and its attribute.
  const std::string impedance = Log(true, sa_library, "law-impedance");
  CHECK_THAT(impedance, ContainsSubstring("has no model for") &&
                            ContainsSubstring(law_line) &&
                            ContainsSubstring("device metal Impedance (Cs=0, Ls=4.5e-14, "
                                              "Rs=0) (boundary attribute(s) 9)") &&
                            ContainsSubstring("library models PEC."));

  // A geometry the library lacks under every law (the spatial-cluster library matches no
  // island feature): Missing, but no law line.
  const std::string missing = Log(true, spatial_cluster_library_3d_path, "law-missing");
  CHECK_THAT(missing,
             ContainsSubstring("has no model for") && !ContainsSubstring(law_line));
#endif
}
