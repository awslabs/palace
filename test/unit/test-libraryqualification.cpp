// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include <filesystem>
#include <fstream>
#include <functional>
#include <string>
#include <vector>
#include <catch2/catch_test_macros.hpp>
#include <nlohmann/json.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>
#include "fem/mesh.hpp"
#include "fixtures.hpp"
#include "models/surfaceresponseoperator.hpp"
#include "surfaceresponse-fixtures.hpp"
#include "utils/communication.hpp"
#include "utils/iodata.hpp"
#include "utils/outputdir.hpp"

using namespace palace;
using namespace Catch::Matchers;
using json = nlohmann::json;
namespace fs = std::filesystem;

// R-5 (fix round P2; OI-LIBRARY-INHERITED-PENDING-MODEL; USER decision 1230 Q8 REJECT,
// supervisor decision 1260; DESIGN-v3 s6): the loader enforces a model's explicit
// QualificationStatus. A status other than "Qualified" is refused at load unless the
// library's top-level AllowPendingModels lists the model by Name with a Decision; a listed
// model is applied exactly as before (bitwise continuity) and the admission is recorded on
// the manifest (Library.AllowPendingModels; QualificationStatus / AllowDecision on the
// SelectedModels entries and the Features[].Match of that model). A status-less model is
// the record form of the 2D planner: allowed, nothing recorded.

namespace
{

json ReadJson(const fs::path &path)
{
  std::ifstream input(path);
  REQUIRE(input);
  return json::parse(input);
}

// Every leaf path of a JSON document (the dry-run comparator's form).
void Leaves(const json &node, const std::string &prefix,
            std::vector<std::pair<std::string, json>> &out)
{
  if (node.is_object())
  {
    for (const auto &[key, value] : node.items())
    {
      Leaves(value, prefix + "/" + key, out);
    }
  }
  else if (node.is_array())
  {
    for (std::size_t i = 0; i < node.size(); i++)
    {
      Leaves(node[i], prefix + "[" + std::to_string(i) + "]", out);
    }
  }
  else
  {
    out.emplace_back(prefix, node);
  }
}

}  // namespace

TEST_CASE_METHOD(test::SurfaceResponseFiles,
                 "R-5 QualificationStatus at load: a non-Qualified model is refused unless "
                 "AllowPendingModels names it, and an admitted model is applied as before "
                 "and recorded (decisions 1230 Q8 / 1260)",
                 "[surfaceresponseoperator][libraryqualification][3d][Serial][Parallel]")
{
#if !defined(MFEM_USE_GSLIB)
  SKIP("SurfaceResponseOperator requires MFEM_USE_GSLIB");
#else
  MPI_Comm comm = Mpi::World();
  json config = ConvexIslandConfig();
  config["Solver"]["Electrostatic"]["ResponseCorrection"]["UnmatchedPolicy"] = "Warn";
  // The Features (identification) path with the convex library keyed on SA only like the
  // island's features (the sheet-library form of test-geometryfallback.cpp): Exact matches.
  config["Solver"]["Electrostatic"]["ResponseCorrection"].erase("PatchConstruction");
  auto WriteVariant =
      [&](const std::string &suffix, const std::function<void(json &)> &edit)
  {
    const fs::path path = temp.temp_dir / ("fabrication-process-r5-" + suffix + ".json");
    if (Mpi::Root(comm))
    {
      json library = ReadJson(convex_library_3d_path);
      library["Name"] = "unit-test-process-r5";
      for (auto &model : library["Models"])
      {
        model["Interfaces"] = {{{"Type", "SA"}, {"Coupon", 1}}};
      }
      edit(library);
      std::ofstream output(path);
      output << library.dump(2) << "\n";
    }
    Mpi::Barrier(comm);
    return path;
  };
  auto Preflight = [&](const fs::path &library, const std::string &tag)
  {
    const fs::path dir = temp.temp_dir / ("r5-" + tag);
    json run_config = config;
    run_config["Problem"]["Output"] = dir.string();
    run_config["Solver"]["Electrostatic"]["ResponseCorrection"]["Library"] =
        library.string();
    IoData iodata(run_config, false);
    iodata.boundaries.cracked_attributes.insert(9);
    MakeOutputFolder(iodata, comm);
    Mesh mesh(MakeIslandMesh());
    const fs::path manifest_path = dir / "surface-response-requirements.json";
    WriteSurfaceResponseRequirements(iodata, mesh, manifest_path.string());
    Mpi::Barrier(comm);
    return std::make_pair(ReadJson(manifest_path), dir / "surface-response-patches.csv");
  };
  auto ReadBytes = [](const fs::path &path)
  {
    std::ifstream input(path, std::ios::binary);
    REQUIRE(input);
    return std::string(std::istreambuf_iterator<char>(input), {});
  };
  // The base: every model status-less (the record form), no allow leaf.
  const auto [base, base_patches] = Preflight(WriteVariant("base", [](json &) {}), "base");
  std::string applied_model;
  for (const auto &requirement : base.at("Requirements"))
  {
    if (requirement.contains("SelectedModels") && !requirement["SelectedModels"].empty())
    {
      applied_model = requirement["SelectedModels"][0]["Name"].get<std::string>();
      break;
    }
  }
  REQUIRE_FALSE(applied_model.empty());
  const std::string decision =
      "decision 315 / 1260: PendingQualification admitted by the curator for the unit test";
  auto SetStatus = [&](json &library, const std::string &status)
  {
    bool found = false;
    for (auto &model : library["Models"])
    {
      if (model["Name"] == applied_model)
      {
        model["QualificationStatus"] = status;
        found = true;
      }
    }
    REQUIRE(found);
  };
  // The manifest minus the leaves R-5 adds and the library path (which names the variant).
  auto Strip = [](json manifest)
  {
    manifest["Library"].erase("Path");
    manifest["Library"].erase("AllowPendingModels");
    for (auto &requirement : manifest["Requirements"])
    {
      if (requirement.contains("SelectedModels"))
      {
        for (auto &model : requirement["SelectedModels"])
        {
          model.erase("QualificationStatus");
          model.erase("AllowDecision");
        }
      }
    }
    for (auto &feature : manifest["Identification"]["Features"])
    {
      if (feature.contains("Match"))
      {
        feature["Match"].erase("QualificationStatus");
        feature["Match"].erase("AllowDecision");
      }
    }
    // The run's own wall time.
    if (manifest["Identification"]["Diagnostics"].contains("DomainBoundaryExclusions"))
    {
      manifest["Identification"]["Diagnostics"]["DomainBoundaryExclusions"].erase(
          "WallTime");
    }
    return manifest;
  };

  SECTION("an explicit non-Qualified status without the allow leaf is refused at load")
  {
    const auto library =
        WriteVariant("pending", [&](json &l) { SetStatus(l, "PendingQualification"); });
    CHECK_THROWS_WITH(Preflight(library, "pending"),
                      ContainsSubstring("QualificationStatus \"PendingQualification\"") &&
                          ContainsSubstring("AllowPendingModels does not list it") &&
                          ContainsSubstring(applied_model));
  }
  SECTION("a malformed allow leaf is refused")
  {
    const auto library =
        WriteVariant("malformed",
                     [&](json &l)
                     {
                       SetStatus(l, "PendingQualification");
                       l["AllowPendingModels"] = {{{"Name", applied_model}}};
                     });
    CHECK_THROWS_WITH(Preflight(library, "malformed"),
                      ContainsSubstring("nonempty Name and a nonempty Decision"));
  }
  SECTION("an admitted model is applied as before and the admission is recorded")
  {
    const auto library = WriteVariant(
        "admitted",
        [&](json &l)
        {
          SetStatus(l, "PendingQualification");
          l["AllowPendingModels"] = {{{"Name", applied_model}, {"Decision", decision}}};
        });
    const auto [manifest, patches] = Preflight(library, "admitted");
    // The library block carries the leaf as read.
    REQUIRE(manifest["Library"].contains("AllowPendingModels"));
    CHECK(manifest["Library"]["AllowPendingModels"] ==
          json({{{"Name", applied_model}, {"Decision", decision}}}));
    // Every SelectedModels entry of that model carries the status and the decision, no
    // other entry does; the Features[].Match of that model likewise.
    int annotated_requirements = 0, annotated_features = 0;
    for (const auto &requirement : manifest["Requirements"])
    {
      for (const auto &model : requirement.value("SelectedModels", json::array()))
      {
        const bool pending = model["Name"] == applied_model;
        CHECK(model.contains("QualificationStatus") == pending);
        CHECK(model.contains("AllowDecision") == pending);
        if (pending)
        {
          CHECK(model["QualificationStatus"] == "PendingQualification");
          CHECK(model["AllowDecision"] == decision);
          CHECK(model["Weight"].get<double>() == 1.0);
          annotated_requirements++;
        }
      }
    }
    for (const auto &feature : manifest["Identification"]["Features"])
    {
      const auto &match = feature["Match"];
      const bool pending = match.contains("Model") && match["Model"] == applied_model;
      CHECK(match.contains("QualificationStatus") == pending);
      CHECK(match.contains("AllowDecision") == pending);
      annotated_features += pending ? 1 : 0;
    }
    CHECK(annotated_requirements > 0);
    CHECK(annotated_features > 0);
    // Bitwise continuity (decision 1260): the patches and every other manifest leaf are
    // the base's (the model is applied as before, weight 1, not weight 0).
    CHECK(ReadBytes(patches) == ReadBytes(base_patches));
    std::vector<std::pair<std::string, json>> a, b;
    Leaves(Strip(base), "", a);
    Leaves(Strip(manifest), "", b);
    CHECK(a == b);
    CHECK(manifest["Summary"]["Counts"] == base["Summary"]["Counts"]);
  }
  SECTION("an explicit Qualified status and a stale admission record nothing new")
  {
    const auto library = WriteVariant(
        "qualified",
        [&](json &l)
        {
          SetStatus(l, "Qualified");
          // A stale entry (the named model is Qualified): a warning, nothing admitted.
          l["AllowPendingModels"] = {{{"Name", applied_model}, {"Decision", decision}}};
        });
    const auto [manifest, patches] = Preflight(library, "qualified");
    CHECK(ReadBytes(patches) == ReadBytes(base_patches));
    for (const auto &requirement : manifest["Requirements"])
    {
      for (const auto &model : requirement.value("SelectedModels", json::array()))
      {
        CHECK_FALSE(model.contains("QualificationStatus"));
      }
    }
    for (const auto &feature : manifest["Identification"]["Features"])
    {
      CHECK_FALSE(feature["Match"].contains("QualificationStatus"));
    }
    // The leaf as read is still recorded (it is the library's data).
    CHECK(manifest["Library"].contains("AllowPendingModels"));
    std::vector<std::pair<std::string, json>> a, b;
    Leaves(Strip(base), "", a);
    Leaves(Strip(manifest), "", b);
    CHECK(a == b);
  }
  SECTION("a library without the key writes no AllowPendingModels leaf")
  {
    CHECK_FALSE(base["Library"].contains("AllowPendingModels"));
  }
#endif
}
