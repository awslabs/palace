// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "fixtures.hpp"
#include "surfaceresponse-fixtures.hpp"

#include <array>
#include <filesystem>
#include <fstream>
#include <memory>
#include <set>
#include <string>
#include <vector>
#include <catch2/catch_test_macros.hpp>
#include <nlohmann/json.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>
#include "fem/mesh.hpp"
#include "models/surfaceresponseidentification.hpp"
#include "models/surfaceresponseoperator.hpp"
#include "utils/communication.hpp"
#include "utils/iodata.hpp"
#include "utils/outputdir.hpp"

using namespace palace;
using namespace Catch::Matchers;
using json = nlohmann::json;
namespace fs = std::filesystem;

// OI-SCP-BOXGROW, the T5 FALLBACK through the identification, the matching pass and the
// requirements manifest (decisions 958 (a) / 960 R4 / 961 MINOR-2 / MINOR-3 / 969 Q3). The
// offset-corner-pair mesh with the straddling island (MakeOffsetCornerPairMesh(true), R =
// 0.2): the 6-edge cluster's T2 / T3 / T4-final box reaches z = 1.93 and x = -0.055; the
// third island's lower corners (0.125, 1.75) / (0.625, 1.75) sit 0.18 inside the z face,
// their +z arm windows exit it by 0.02 (inside 1.9 R, outside 0.1 R of each 2 R window:
// T5 fires on both), and the R-band of the left corner's +z arm reaches x = -0.075, 0.02
// beyond the x face: ONE step of 0.25 R on both faces (the cluster's frame x1 / y1 =
// device -x / +z). Three Signature-keyed libraries built from the mesh's own
// identification (one model per feature hash, the fixture model of the topology with the
// feature's Signature; the cluster model signature-only): the UNGROWN library carries the
// cluster's pre-T5 signature (SpatialSupport.CornerContainment.Ungrown.Signature: the key
// main's binary emits) -> the fallback places it (the feature's Hash / Signature read the
// ungrown key, Match.Model the ungrown model, CornerContainment.Fallback {GrownKey,
// PlacedKey, PlacedModel, Reason NoMatchedModel} + Grown; the manifest lists the grown key
// as a Missing-class requirement flagged FallbackPlaced with PlacedModel, excluded from
// Summary.Counts.Missing (0) and Complete (true), carried by Summary.FallbackPlaced); the
// GROWN library carries the grown signature -> an exact match, no fallback, no
// FallbackPlaced; a library WITHOUT the cluster -> the grown key Missing (counted), the
// Ungrown subtree kept, no Fallback. The dry-run patch tables of the UNGROWN and GROWN
// placements agree (one model geometry, one frame of record).
TEST_CASE_METHOD(test::SurfaceResponseFiles,
                 "T5 excluded-corner containment fallback through the matching pass and "
                 "the requirements manifest (decisions 960 R4 / 961 / 969)",
                 "[surfaceresponseoperator][cornercontainment][3d][Serial][Parallel]")
{
#if !defined(MFEM_USE_GSLIB)
  SKIP("SurfaceResponseOperator requires MFEM_USE_GSLIB");
#else
  constexpr double R = 0.2;
  MPI_Comm comm = Mpi::World();
  const fs::path mesh_path = temp.temp_dir / "t5-offset-corner-pair-straddling.mesh";
  const fs::path ungrown_library = temp.temp_dir / "surface-process-t5-ungrown.json";
  const fs::path grown_library = temp.temp_dir / "surface-process-t5-grown.json";
  const fs::path no_cluster_library = temp.temp_dir / "surface-process-t5-nocluster.json";
  std::vector<std::unique_ptr<Mesh>> meshes;
  meshes.push_back(std::make_unique<Mesh>(MakeOffsetCornerPairMesh(true)));
  auto Config = [&](const fs::path &library, const fs::path &output)
  {
    json config = ConvexMaxwellIslandConfig();
    config["Model"]["Mesh"] = mesh_path.string();
    config["Model"]["L0"] = 1.0;
    config["Problem"]["Output"] = output.string();
    config["Boundaries"]["Ground"]["Attributes"] = {1, 2, 3, 4, 5, 6, 9, 10, 11};
    config["Boundaries"]["Postprocessing"]["Dielectric"][0]["Attributes"] = {9, 10, 11};
    auto &correction = config["Solver"]["SurfaceResponseCorrection"];
    correction["Library"] = library.string();
    correction["CorrectionMode"] = "PostprocessOnly";
    correction["UnmatchedPolicy"] = "Warn";
    correction.erase("PatchConstruction");
    return config;
  };
  auto MakeIoData = [&](json config)
  {
    IoData iodata(std::move(config), false);
    for (const int attribute : {9, 10, 11})
    {
      iodata.boundaries.cracked_attributes.insert(attribute);
    }
    return iodata;
  };
  auto Preflight = [&](const fs::path &library, const char *name)
  {
    const fs::path dir = temp.temp_dir / name;
    IoData iodata = MakeIoData(Config(library, dir));
    MakeOutputFolder(iodata, comm);
    WriteSurfaceResponseRequirements(iodata, *meshes.back(),
                                     (dir / "surface-response-requirements.json").string());
    Mpi::Barrier(comm);
    std::ifstream input(dir / "surface-response-requirements.json");
    REQUIRE(input);
    return std::make_pair(dir, json::parse(input));
  };
  auto ClusterFeature = [](const json &manifest) -> const json &
  {
    const auto &features = manifest["Identification"]["Features"];
    const auto it = std::find_if(features.begin(), features.end(), [](const auto &feature)
                                 { return feature["Type"] == "SpatialEdgeCluster"; });
    REQUIRE(it != features.end());
    return *it;
  };

  // (0) The identification against the fixture's Edges-keyed library: every feature
  //     Missing; T5 acts on the one cluster (the straddling island's lower corners).
  if (Mpi::Root(comm))
  {
    mfem::Mesh serial = MakeSerialOffsetCornerPairMesh(true);
    std::ofstream output(mesh_path);
    serial.Print(output);
  }
  Mpi::Barrier(comm);
  const auto [identification_dir, identification] =
      Preflight(spatial_cluster_library_3d_path, "t5-identification");
  const json &cluster = ClusterFeature(identification);
  REQUIRE(cluster["SpatialSupport"].contains("CornerContainment"));
  const json containment = cluster["SpatialSupport"]["CornerContainment"];
  REQUIRE(containment.contains("Ungrown"));
  CHECK(!containment.contains("Fallback"));
  CHECK(containment["Steps"] == json({0, 0, 1, 1}));
  REQUIRE(containment["Vertices"].size() == 2);
  // Per vertex, its OWN band and steps (decision 972, review MAJOR-1): the left corner
  // (0.125, 1.75) is short on the x face too (its +z arm's band reaches x = -0.075), the
  // right corner (0.625, 1.75) on the z face only; the cluster's Steps is the max.
  std::multiset<std::array<int, 4>> vertex_steps;
  for (const auto &vertex : containment["Vertices"])
  {
    CHECK(vertex["Type"] == "ConvexCorner");
    CHECK_THAT(vertex["InsideOverR"].get<double>(), WithinAbs(1.9, 1.0e-3));
    CHECK_THAT(vertex["OutsideEveryOverR"].get<double>(), WithinAbs(0.1, 1.0e-3));
    vertex_steps.insert(vertex["Steps"].get<std::array<int, 4>>());
    // The band is the vertex's own: 2 R wide across its +z arm (x in [vx - R, vx + R]).
    const auto &band = vertex["BandOverR"];
    CHECK_THAT(band[2].get<double>() - band[0].get<double>(), WithinAbs(2.0, 1.0e-3));
  }
  CHECK(vertex_steps == std::multiset<std::array<int, 4>>{{0, 0, 1, 1}, {0, 0, 0, 1}});
  CHECK(cluster["SpatialSupport"]["Contract"] == 3);
  CHECK(cluster["SpatialSupport"]["Growth"]["Grown"] == true);
  CHECK(cluster["SpatialSupport"]["Growth"]["Steps"] == json({0, 0, 0, 0}));
  const std::string grown_hash = cluster["Hash"];
  const std::string ungrown_hash = containment["Ungrown"]["Key"];
  CHECK(containment["GrownKey"] == grown_hash);
  CHECK(grown_hash != ungrown_hash);
  CHECK(containment["Ungrown"]["Contract"] == 3);  // the islands are context
  CHECK(containment["Ungrown"]["SpatialSupport"]["Box"] ==
        cluster["SpatialSupport"]["ClaimsBox"]);
  // Every cluster requirement is the grown key.
  int cluster_requirements = 0;
  for (const auto &requirement : identification["Requirements"])
  {
    if (requirement["Topology"] == "SpatialEdgeCluster")
    {
      cluster_requirements++;
      CHECK(requirement["Hash"] == grown_hash);
      CHECK(requirement["Status"] == "Missing");
      CHECK(!requirement.contains("FallbackPlaced"));
    }
  }
  CHECK(cluster_requirements == 1);
  CHECK(!identification["Summary"].contains("FallbackPlaced"));

  // The three Signature-keyed libraries.
  if (Mpi::Root(comm))
  {
    std::ifstream library_input(spatial_cluster_library_3d_path);
    REQUIRE(library_input);
    const json fixture_library = json::parse(library_input);
    REQUIRE(fixture_library["Models"].size() == 3);
    REQUIRE(fixture_library["Models"][0]["Topology"] == "IsolatedEdge");
    REQUIRE(fixture_library["Models"][1]["Topology"] == "ConvexCorner");
    REQUIRE(fixture_library["Models"][2]["Name"] == "offset-corner-pair");
    auto WriteLibrary =
        [&](const json *cluster_signature, const fs::path &path, const char *suffix)
    {
      json library = fixture_library;
      library["Name"] = std::string("unit-test-process-spatial-cluster-t5-") + suffix;
      library["Models"] = json::array();
      std::set<std::string> hashes;
      for (const auto &feature : identification["Identification"]["Features"])
      {
        if (!hashes.insert(feature["Hash"].get<std::string>()).second)
        {
          continue;
        }
        const std::string type = feature["Type"];
        REQUIRE((type == "IsolatedEdge" || type == "ConvexCorner" ||
                 type == "SpatialEdgeCluster"));
        if (type == "SpatialEdgeCluster" && !cluster_signature)
        {
          continue;
        }
        json model = fixture_library["Models"][type == "IsolatedEdge"   ? 0
                                               : type == "ConvexCorner" ? 1
                                                                        : 2];
        model["Signature"] = feature["Signature"];
        if (type == "ConvexCorner")
        {
          model["Angle"] = feature["Signature"]["AngleDegrees"];
          model["Name"] = "convex-corner-90";
        }
        else if (type == "SpatialEdgeCluster")
        {
          model.erase("Edges");
          model["Name"] = std::string("offset-corner-pair-t5-") + suffix;
          model["Signature"] = *cluster_signature;
          json support = json::array();
          for (const double x : {-1.3, 1.3})
          {
            for (const double y : {-1.3, 1.3})
            {
              for (const double z : {-0.05, 0.05})
              {
                support.push_back({x, y, z});
              }
            }
          }
          model["SupportPoints"] = support;
        }
        library["Models"].push_back(std::move(model));
      }
      std::ofstream output(path);
      output << library.dump(2) << "\n";
    };
    const json ungrown_signature = containment["Ungrown"]["Signature"];
    const json grown_signature = cluster["Signature"];
    WriteLibrary(&ungrown_signature, ungrown_library, "ungrown");
    WriteLibrary(&grown_signature, grown_library, "grown");
    WriteLibrary(nullptr, no_cluster_library, "nocluster");
  }
  Mpi::Barrier(comm);

  // (1) The UNGROWN library: the fallback.
  const auto [ungrown_dir, ungrown] = Preflight(ungrown_library, "t5-ungrown");
  {
    const json &placed = ClusterFeature(ungrown);
    CHECK(placed["Hash"] == ungrown_hash);
    CHECK(placed["Signature"] == containment["Ungrown"]["Signature"]);
    REQUIRE(placed.contains("Match"));
    CHECK(placed["Match"]["Model"] == "offset-corner-pair-t5-ungrown");
    CHECK(placed["Match"]["Status"] != "Missing");
    CHECK_THAT(placed["Match"]["Note"].get<std::string>(),
               ContainsSubstring("T5 fallback (decision 960 R4)"));
    const json &support = placed["SpatialSupport"];
    CHECK(support["Contract"] == 3);
    CHECK(support["Box"] == support["ClaimsBox"]);  // the ungrown record is the record
    REQUIRE(support.contains("CornerContainment"));
    const json &fallback_containment = support["CornerContainment"];
    CHECK(!fallback_containment.contains("Ungrown"));
    REQUIRE(fallback_containment.contains("Fallback"));
    REQUIRE(fallback_containment.contains("Grown"));
    const json &fallback = fallback_containment["Fallback"];
    CHECK(fallback["GrownKey"] == grown_hash);
    CHECK(fallback["PlacedKey"] == ungrown_hash);
    CHECK(fallback["PlacedModel"] == "offset-corner-pair-t5-ungrown");
    CHECK(fallback["Reason"] == "NoMatchedModel");
    CHECK(fallback_containment["Grown"]["Key"] == grown_hash);
    CHECK(fallback_containment["Grown"]["Contract"] == 3);
    CHECK(fallback_containment["Grown"]["Signature"] == cluster["Signature"]);
    CHECK(fallback_containment["Grown"]["SpatialSupport"]["Box"] ==
          cluster["SpatialSupport"]["Box"]);
    CHECK(!fallback_containment["Grown"]["SpatialSupport"].contains("CornerContainment"));
    CHECK(fallback_containment["Steps"] == containment["Steps"]);
    CHECK(fallback_containment["Vertices"] == containment["Vertices"]);
    // The requirements: the grown key Missing-class, FallbackPlaced, excluded from the
    // counts; every other requirement Exact.
    int grown_entries = 0, missing_entries = 0;
    for (const auto &requirement : ungrown["Requirements"])
    {
      if (requirement["Hash"] == grown_hash)
      {
        grown_entries++;
        CHECK(requirement["Status"] == "Missing");
        CHECK(requirement["FallbackPlaced"] == true);
        CHECK(requirement["PlacedModel"] == "offset-corner-pair-t5-ungrown");
        CHECK(requirement["Count"] == 1);
        CHECK(requirement["Instances"] == 1);
        CHECK(requirement["Signature"] == cluster["Signature"]);
        CHECK(requirement["Geometry"]["Signature"] == cluster["Signature"]);
        CHECK(!requirement.contains("SelectedModels"));
      }
      else
      {
        // The ungrown key among them, Exact: today's entry (the placed model's).
        CHECK(requirement["Status"] == "Exact");
        CHECK(!requirement.contains("FallbackPlaced"));
      }
      missing_entries += requirement["Status"] == "Missing" ? 1 : 0;
    }
    CHECK(grown_entries == 1);
    CHECK(missing_entries == 1);
    CHECK(ungrown["Summary"]["Counts"]["Missing"] == 0);
    CHECK(ungrown["Summary"]["Counts"]["Exact"] ==
          identification["Summary"]["Counts"]["Missing"]);
    CHECK(ungrown["Complete"] == true);
    REQUIRE(ungrown["Summary"].contains("FallbackPlaced"));
    CHECK(ungrown["Summary"]["FallbackPlaced"]["Count"] == 1);
    REQUIRE(ungrown["Summary"]["FallbackPlaced"]["Keys"].size() == 1);
    CHECK(ungrown["Summary"]["FallbackPlaced"]["Keys"][0]["Hash"] == grown_hash);
    CHECK(ungrown["Summary"]["FallbackPlaced"]["Keys"][0]["PlacedModel"] ==
          "offset-corner-pair-t5-ungrown");
  }

  // (2) The GROWN library: an exact match, no fallback.
  const auto [grown_dir, grown] = Preflight(grown_library, "t5-grown");
  {
    const json &placed = ClusterFeature(grown);
    CHECK(placed["Hash"] == grown_hash);
    CHECK(placed["Match"]["Model"] == "offset-corner-pair-t5-grown");
    const json &support = placed["SpatialSupport"];
    CHECK(support["Box"] == cluster["SpatialSupport"]["Box"]);
    REQUIRE(support.contains("CornerContainment"));
    CHECK(support["CornerContainment"].contains("Ungrown"));
    CHECK(!support["CornerContainment"].contains("Fallback"));
    for (const auto &requirement : grown["Requirements"])
    {
      CHECK(requirement["Status"] == "Exact");
      CHECK(!requirement.contains("FallbackPlaced"));
    }
    CHECK(grown["Summary"]["Counts"]["Missing"] == 0);
    CHECK(grown["Summary"]["Counts"]["Exact"] == ungrown["Summary"]["Counts"]["Exact"]);
    CHECK(!grown["Summary"].contains("FallbackPlaced"));
    CHECK(grown["Complete"] == true);
  }
  // The dry-run patch tables of the two placements agree (one model geometry, one frame):
  // every field of every row but the cluster patches' Model name (column 3 of the header
  // Patch, Feature, Topology, Model, ...; ReadPatchRows skips the header).
  if (Mpi::Root(comm))
  {
    const auto ungrown_rows = ReadPatchRows(ungrown_dir / "surface-response-patches.csv");
    const auto grown_rows = ReadPatchRows(grown_dir / "surface-response-patches.csv");
    REQUIRE(ungrown_rows.size() == grown_rows.size());
    REQUIRE(ungrown_rows.size() > 1);
    constexpr std::size_t model_column = 3;
    int cluster_rows = 0;
    for (std::size_t i = 0; i < ungrown_rows.size(); i++)
    {
      REQUIRE(ungrown_rows[i].size() == grown_rows[i].size());
      for (std::size_t k = 0; k < ungrown_rows[i].size(); k++)
      {
        if (k == model_column)
        {
          if (ungrown_rows[i][k] == "offset-corner-pair-t5-ungrown")
          {
            cluster_rows++;
            CHECK(grown_rows[i][k] == "offset-corner-pair-t5-grown");
            continue;
          }
        }
        CHECK(ungrown_rows[i][k] == grown_rows[i][k]);
      }
    }
    CHECK(cluster_rows > 0);
  }

  // (3) No cluster model: the grown key Missing (counted), the Ungrown subtree kept.
  const auto [no_cluster_dir, no_cluster] = Preflight(no_cluster_library, "t5-nocluster");
  {
    const json &unplaced = ClusterFeature(no_cluster);
    CHECK(unplaced["Hash"] == grown_hash);
    CHECK(unplaced["SpatialSupport"]["CornerContainment"].contains("Ungrown"));
    CHECK(!unplaced["SpatialSupport"]["CornerContainment"].contains("Fallback"));
    int missing_entries = 0;
    for (const auto &requirement : no_cluster["Requirements"])
    {
      if (requirement["Status"] == "Missing")
      {
        missing_entries++;
        CHECK(requirement["Hash"] == grown_hash);
        CHECK(!requirement.contains("FallbackPlaced"));
      }
    }
    CHECK(missing_entries == 1);
    CHECK(no_cluster["Summary"]["Counts"]["Missing"] == 1);
    CHECK(no_cluster["Complete"] == false);
    CHECK(!no_cluster["Summary"].contains("FallbackPlaced"));
  }
#endif
}
