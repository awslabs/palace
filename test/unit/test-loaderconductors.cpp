// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "fixtures.hpp"
#include "surfaceresponse-fixtures.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <memory>
#include <set>
#include <string>
#include <string_view>
#include <vector>
#include <mfem.hpp>
#include <catch2/catch_test_macros.hpp>
#include <nlohmann/json.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>
#include "fem/gridfunction.hpp"
#include "fem/mesh.hpp"
#include "models/spaceoperator.hpp"
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

// A device whose spatial cluster has a CONTEXT-ONLY conductor (decision 282's class; the
// 448693d60a6f shape of decision 649 (1)): on the plane y = 0.5 of a 2 x 1 x 2 box (0.125
// hexahedra, R = 0.2 = 0.625 cells), a PEC island (attribute 9) x in [0.125, 1.125], z in
// [0.25, 0.75] with a one-cell step x in [0.125, 0.625], z in [0.75, 0.875]: the step's
// concave and convex corners, 0.625 R apart, form a three-edge spatial cluster whose
// claims all lie on conductor 1 (the island's other corners are 2.5 R away and stay
// isolated). A second PEC strip (attribute 10) x in [1.375, 1.875], z in [1.25, 1.5]
// has no vertex within 2 R of the cluster, so it claims nothing there, but its left corner
// lies inside the cluster's support box: two of its edges are foreign context pieces
// (Chain false, the next conductor label 2). The strip's own ends are further clusters of
// their own (two corners 1.25 R apart), left unmatched here. The strip sits 2.5 R above
// the island's right corner (1.125, 0.75) along z: under the oriented separating-axis
// vertex rule of the identification fix round (decision 1062 A.5 (a)) two axis-aligned
// +-R cubes 2.5 R apart along one axis are disjoint, so the island corner and the strip's
// left corner stay separate features (at the former z in [1.125, 1.375] the strip's corner
// lay 1.25 R x 1.875 R from it — inside its cube — and the rule clusters them).
std::unique_ptr<mfem::ParMesh> MakeStepWithForeignStripMesh()
{
  constexpr double extent = 2.0;
  mfem::Mesh serial = mfem::Mesh::MakeCartesian3D(16, 4, 16, mfem::Element::HEXAHEDRON,
                                                  extent, 1.0, extent);
  for (int face = 0; face < serial.GetNumFaces(); face++)
  {
    int element1, element2;
    serial.GetFaceElements(face, &element1, &element2);
    if (element1 < 0 || element2 < 0)
    {
      continue;
    }
    mfem::Array<int> vertices;
    serial.GetFaceVertices(face, vertices);
    bool on_plane = true;
    double xmin = extent, xmax = 0.0;
    double zmin = extent, zmax = 0.0;
    for (const int vertex : vertices)
    {
      const double *point = serial.GetVertex(vertex);
      on_plane = on_plane && std::abs(point[1] - 0.5) < 1.0e-12;
      xmin = std::min(xmin, point[0]);
      xmax = std::max(xmax, point[0]);
      zmin = std::min(zmin, point[2]);
      zmax = std::max(zmax, point[2]);
    }
    auto Inside = [&](double x0, double x1, double z0, double z1)
    {
      return xmin >= x0 - 1.0e-12 && xmax <= x1 + 1.0e-12 && zmin >= z0 - 1.0e-12 &&
             zmax <= z1 + 1.0e-12;
    };
    const bool island =
        Inside(0.125, 1.125, 0.25, 0.75) || Inside(0.125, 0.625, 0.75, 0.875);
    const bool strip = Inside(1.375, 1.875, 1.25, 1.5);
    if (on_plane && (island || strip))
    {
      serial.AddBdrElement(serial.GetFace(face)->Duplicate(&serial));
      serial.SetBdrAttribute(serial.GetNBE() - 1, island ? 9 : 10);
    }
  }
  serial.FinalizeTopology();
  serial.Finalize();
  return std::make_unique<mfem::ParMesh>(Mpi::World(), serial);
}

// The island step is inside conductor 1; the strip is conductor 2 (plan coordinates (x, z)
// of the mesh above).
bool InsideIslandStep(const std::array<double, 3> &p)
{
  return p[0] > 0.125 && p[0] < 0.625 && p[2] > 0.75 && p[2] < 0.875 &&
         std::abs(p[1] - 0.5) < 1.0e-9;
}
bool InsideStrip(const std::array<double, 3> &p)
{
  return p[0] > 1.375 && p[0] < 1.875 && p[2] > 1.25 && p[2] < 1.5 &&
         std::abs(p[1] - 0.5) < 1.0e-9;
}

// 0.1 R inside the metal of a straight signature piece {P, Gap} (the generator's reference
// rule: the gap direction points from the metal into the gap), library units (= R x the
// signature's units of R).
std::array<double, 3> ReferenceInsidePiece(const json &piece, double radius)
{
  const auto P = piece["P"].get<std::array<double, 4>>();
  const auto gap = piece["Gap"].get<std::array<double, 2>>();
  return {(0.5 * (P[0] + P[2]) - 0.1 * gap[0]) * radius,
          (0.5 * (P[1] + P[3]) - 0.1 * gap[1]) * radius, 0.0};
}

// Integer plan-view segments (the loader's PlanViewBoundary contract) of straight signature
// pieces on the 1e-6 R grid.
json PlanViewComponent(int conductor, const std::vector<json> &pieces, double radius)
{
  json segments = json::array();
  for (const auto &piece : pieces)
  {
    const auto P = piece["P"].get<std::array<double, 4>>();
    auto Q = [&](double value)
    { return static_cast<long long int>(std::llround(value * radius * 1.0e6)); };
    segments.push_back({{Q(P[0]), Q(P[1]), 0}, {Q(P[2]), Q(P[3]), 0}});
  }
  return {{"Conductor", conductor}, {"Segments", segments}};
}

void WriteJsonOnRoot(const fs::path &path, const json &value)
{
  if (Mpi::Root(Mpi::World()))
  {
    std::ofstream output(path);
    output << value.dump(2) << "\n";
  }
  Mpi::Barrier(Mpi::World());
}

json ReadJsonAfterBarrier(const fs::path &path)
{
  Mpi::Barrier(Mpi::World());
  std::ifstream input(path);
  REQUIRE(input);
  return json::parse(input);
}

}  // namespace

TEST_CASE_METHOD(test::SurfaceResponseFiles,
                 "SurfaceResponseOperator spatial cluster with a context-only conductor",
                 "[surfaceresponseoperator][3d][maxwell][spatial][placement]"
                 "[loaderconductors][Serial][Parallel]")
{
#if !defined(MFEM_USE_GSLIB)
  SKIP("SurfaceResponseOperator requires MFEM_USE_GSLIB");
#else
  // Decision 649 (1) (the rerun-3b diagnosis of 448693d60a6f, F-643-1): a device-plan
  // spatial cluster model may carry a context-only conductor, a plan conductor with no
  // claimed Edge. The generator emits one ConductorReference per PLAN conductor (the
  // response carries that conductor's state: ConductorStates 1 for a two-conductor
  // coupon, OpenContourPaths 1 -> 2), so the loader counts the canonical conductors over
  // the Edges AND the PlanViewBoundary components. Such a model LOADS, MATCHES its feature
  // by Signature and is PLACED with its second anchor inside the claim-less conductor; a
  // plan of one conductor with two references, or a single Reference with a plan of two
  // conductors, stays refused.
  json config = ConvexMaxwellIslandConfig();
  config["Boundaries"]["Ground"]["Attributes"] = {1, 2, 3, 4, 5, 6, 9, 10};
  config["Boundaries"]["Postprocessing"]["Dielectric"][0]["Attributes"] = {9, 10};
  auto &correction = config["Solver"]["SurfaceResponseCorrection"];
  correction["Library"] = spatial_cluster_library_3d_path.string();
  correction["UnmatchedPolicy"] = "Warn";
  correction.erase("PatchConstruction");  // Features (the default): the Signature path
  auto MakeIoData = [&](const fs::path &library_path)
  {
    auto library_config = config;
    library_config["Solver"]["SurfaceResponseCorrection"]["Library"] =
        library_path.string();
    IoData iodata(library_config, false);
    iodata.boundaries.cracked_attributes.insert(9);
    iodata.boundaries.cracked_attributes.insert(10);
    return iodata;
  };
  std::vector<std::unique_ptr<Mesh>> meshes;
  meshes.push_back(std::make_unique<Mesh>(MakeStepWithForeignStripMesh()));
  const Mesh &mesh = *meshes.back();

  // The identification of the device against the fixture library: the step cluster's
  // claims are all on conductor 1 and its context carries conductor 2 (Chain false).
  const auto missing_manifest_path =
      temp.temp_dir / "surface-response-requirements-context-only-missing.json";
  {
    IoData missing_iodata = MakeIoData(spatial_cluster_library_3d_path);
    WriteSurfaceResponseRequirements(missing_iodata, mesh, missing_manifest_path.string());
  }
  const json missing_manifest = ReadJsonAfterBarrier(missing_manifest_path);
  REQUIRE(missing_manifest["Version"] == 2);
  const double radius = missing_manifest["Library"]["MatchingRadius"].get<double>();
  CHECK_THAT(radius, WithinAbs(0.2, 1.0e-12));
  auto IsContextOnlyConductorCluster = [](const json &feature)
  {
    if (feature["Type"] != "SpatialEdgeCluster")
    {
      return false;
    }
    const auto &signature = feature["Signature"];
    const bool claims_on_one =
        std::all_of(signature["Portions"].begin(), signature["Portions"].end(),
                    [](const json &portion) { return portion["Conductor"] == 1; });
    const auto &context = signature.value("Context", json::array());
    const bool foreign_two =
        std::any_of(context.begin(), context.end(), [](const json &piece)
                    { return piece["Chain"] == false && piece["Conductor"] == 2; });
    return claims_on_one && foreign_two;
  };
  const auto &missing_features = missing_manifest["Identification"]["Features"];
  const auto feature = std::find_if(missing_features.begin(), missing_features.end(),
                                    IsContextOnlyConductorCluster);
  REQUIRE(feature != missing_features.end());
  CHECK((*feature)["Match"]["Status"] == "Missing");
  const json signature = (*feature)["Signature"];
  REQUIRE(signature["Portions"].size() == 3);
  REQUIRE(signature.contains("Box"));
  CHECK((*feature)["SpatialSupport"]["Contract"] == 3);
  CHECK((*feature)["SpatialSupport"]["Context"]["ForeignConductors"] == 1);
  std::vector<json> claimed_pieces, foreign_pieces;
  for (const auto &portion : signature["Portions"])
  {
    claimed_pieces.push_back(portion);
  }
  for (const auto &piece : signature["Context"])
  {
    if (piece["Chain"] == false)
    {
      CHECK(piece["Conductor"] == 2);
      foreign_pieces.push_back(piece);
    }
  }
  REQUIRE(foreign_pieces.size() == 2);

  // The context-only-conductor model: the fixture's two-conductor cluster coupon
  // (matrices with one conductor state, OpenContourPaths 1 -> 2, two references) keyed by
  // the feature's Signature, its Edges the three conductor-1 claims in the canonical
  // frame, its plan view naming conductors 1 and 2, its references 0.1 R inside the
  // claimed metal and inside the context-only conductor.
  std::ifstream fixture_input(spatial_cluster_library_3d_path);
  REQUIRE(fixture_input);
  auto library = json::parse(fixture_input);
  library["Name"] = "unit-test-process-spatial-cluster-context-only-conductor-3d";
  auto &model = library["Models"].back();
  REQUIRE(model["Name"] == "offset-corner-pair");
  REQUIRE(model["ConductorReferences"].size() == 2);
  REQUIRE(model["OpenContourPaths"].size() == 2);
  model["Name"] = "step-with-context-only-conductor";
  model["Signature"] = signature;
  model["Interfaces"] = SignatureInterfaces(signature);
  model["Edges"] = ChordedSignatureEdges(signature, radius);
  REQUIRE(model["Edges"].size() == 3);
  for (const auto &edge : model["Edges"])
  {
    CHECK(edge["Conductor"] == 1);
  }
  const auto reference_1 = ReferenceInsidePiece(claimed_pieces[0], radius);
  const auto reference_2 = ReferenceInsidePiece(foreign_pieces[0], radius);
  model["ConductorReferences"] = {reference_1, reference_2};
  const json two_conductor_plan_view = {PlanViewComponent(1, claimed_pieces, radius),
                                        PlanViewComponent(2, foreign_pieces, radius)};
  model["PlanViewBoundary"] = two_conductor_plan_view;
  model["MaskRegularization"] = {{"Version", 1},
                                 {"PhysicalBoundary", "TaperAndRound"},
                                 {"ContinuationBoundary", "Vertical"}};
  const auto library_path =
      temp.temp_dir / "surface-process-spatial-cluster-context-only-conductor-3d.json";
  WriteJsonOnRoot(library_path, library);
  constexpr std::string_view count_message =
      "one conductor reference for each canonical conductor label";

  // LOADS and MATCHES: the feature is Matched by key to the model (the rule before this
  // fix refused the library at the read: references 2 against the claimed label 1).
  const auto matched_manifest_path =
      temp.temp_dir / "surface-response-requirements-context-only-matched.json";
  {
    IoData matched_iodata = MakeIoData(library_path);
    WriteSurfaceResponseRequirements(matched_iodata, mesh, matched_manifest_path.string());
  }
  const json matched_manifest = ReadJsonAfterBarrier(matched_manifest_path);
  const auto &matched_identification = matched_manifest["Identification"];
  const auto &matched_features = matched_identification["Features"];
  const auto matched = std::find_if(matched_features.begin(), matched_features.end(),
                                    IsContextOnlyConductorCluster);
  REQUIRE(matched != matched_features.end());
  CHECK((*matched)["Hash"] == (*feature)["Hash"]);
  CHECK((*matched)["Match"]["Status"] == "Matched");
  CHECK((*matched)["Match"]["Model"] == "step-with-context-only-conductor");
  const auto &requirements = matched_manifest["Requirements"];
  const auto requirement =
      std::find_if(requirements.begin(), requirements.end(),
                   [&](const json &entry) { return entry["Hash"] == (*feature)["Hash"]; });
  REQUIRE(requirement != requirements.end());
  CHECK((*requirement)["Status"] == "Exact");
  CHECK((*requirement)["Topology"] == "SpatialEdgeCluster");

  // PLACES: the patch dry run places the model with the identity map on the feature's
  // frame (every model edge endpoint on a claimed portion), and its two conductor anchors
  // land inside the island step (conductor 1) and inside the strip (the claim-less
  // conductor 2).
  Mpi::Barrier(Mpi::World());
  const auto patch_rows = ReadPatchRows(temp.temp_dir / "surface-response-patches.csv");
  CheckPlacedModelEdges(*matched, matched_identification, patch_rows, model["Edges"],
                        1.0e-6 * radius);
  const auto patch =
      std::find_if(patch_rows.begin(), patch_rows.end(), [&](const auto &row)
                   { return std::stoi(row[1]) == (*matched)["Id"].get<int>(); });
  REQUIRE(patch != patch_rows.end());
  CHECK((*patch)[3] == "step-with-context-only-conductor");
  std::array<double, 3> origin{};
  std::array<std::array<double, 3>, 3> axes{};
  for (int d = 0; d < 3; d++)
  {
    origin[d] = std::stod((*patch)[13 + d]);
    for (int a = 0; a < 3; a++)
    {
      axes[a][d] = std::stod((*patch)[16 + 3 * a + d]);
    }
  }
  auto Anchor = [&](const std::array<double, 3> &reference)
  {
    std::array<double, 3> anchor = origin;
    for (int a = 0; a < 3; a++)
    {
      for (int d = 0; d < 3; d++)
      {
        anchor[d] += reference[a] * axes[a][d];
      }
    }
    return anchor;
  };
  const auto anchor_1 = Anchor(reference_1);
  const auto anchor_2 = Anchor(reference_2);
  INFO("anchor 1 " << anchor_1[0] << ' ' << anchor_1[1] << ' ' << anchor_1[2]);
  INFO("anchor 2 " << anchor_2[0] << ' ' << anchor_2[1] << ' ' << anchor_2[2]);
  CHECK(InsideIslandStep(anchor_1));
  CHECK(InsideStrip(anchor_2));
  CHECK_FALSE(InsideStrip(anchor_1));
  CHECK_FALSE(InsideIslandStep(anchor_2));

  // The operator: the model is placed once, its basis keeps the conductor-2 state (the
  // four trace knots + ConductorStates 1 of the fixture matrices), and the Maxwell response
  // of the placed coupon is finite.
  IoData iodata = MakeIoData(library_path);
  SpaceOperator space(iodata, meshes);
  {
    SurfaceResponseOperator response(iodata, space);
    CHECK(response.GetPatchCount() > 0);
    const auto statistics = response.GetStatistics();
    const auto &catalog = statistics["ModelCatalog"];
    const auto entry =
        std::find_if(catalog.begin(), catalog.end(), [](const json &record)
                     { return record["Name"] == "step-with-context-only-conductor"; });
    REQUIRE(entry != catalog.end());
    CHECK((*entry)["PatchCount"] == 1);
    CHECK((*entry)["BasisSize"] == 5);
    GridFunction field(space.GetNDSpace(), true);
    mfem::VectorConstantCoefficient field_coefficient = ConstantFieldCoefficient();
    field.Real().ProjectCoefficient(field_coefficient);
    field.Imag() = 0.0;
    const auto result = response.GetMaxwellResponse(field, 0.0);
    CHECK(std::isfinite(result.domain_correction));
    CHECK(result.matched_length_fraction > 0.0);
    CHECK(result.matched_length_fraction < 1.0);  // the strip's ends stay unmatched
  }

  // The geometric near-match PLACES its model (idfix DESIGN section 5 (b) / A.4; decisions
  // 1079 Q4 / 1088 / 1089): the same coupon keyed by a Signature whose CONTEXT differs from
  // the feature's inside the geometric context band (one context piece end moved 5e-2 R
  // along x; the Context entries REVERSED = permuted; its plan view and second reference
  // from the moved pieces) reads GeometricNearMatched AND is placed: the placement's
  // context checks H5a (the agreement) and H5b (the model's context piece ends on the
  // device edges) read the geometric band for that status, not the quantum near-match / the
  // signature tolerance (the exec-B twin b4bc28485f58, context 7.7e-2 R permuted, was
  // matched and then PlacementRejected at H5a before 1088). A context moved just OUTSIDE
  // the band (1.1e-1 R) is Missing: the match refuses it, so the placement never sees it
  // (both read ClusterGeometricNearDifference: the placement can accept nothing the match
  // refused).
  {
    auto NearContextLibrary =
        [&](double shift_over_R, bool permute, const std::string &name)
    {
      auto near_library = library;
      near_library["Name"] = "unit-test-process-spatial-cluster-" + name + "-3d";
      auto &near_model = near_library["Models"].back();
      near_model["Name"] = name;
      json near_signature = signature;
      auto &context = near_signature["Context"];
      REQUIRE(context.size() >= 2);
      bool shifted = false;
      std::vector<json> near_foreign;
      for (auto &piece : context)
      {
        if (!shifted && piece["Chain"] == false)
        {
          piece["P"][0] = piece["P"][0].get<double>() + shift_over_R;
          piece["P"][2] = piece["P"][2].get<double>() + shift_over_R;
          shifted = true;
        }
        if (piece["Chain"] == false)
        {
          near_foreign.push_back(piece);
        }
      }
      REQUIRE(shifted);
      if (permute)
      {
        std::reverse(context.begin(), context.end());
      }
      near_model["Signature"] = near_signature;
      near_model["ConductorReferences"] = {reference_1,
                                           ReferenceInsidePiece(near_foreign[0], radius)};
      near_model["PlanViewBoundary"] = {PlanViewComponent(1, claimed_pieces, radius),
                                        PlanViewComponent(2, near_foreign, radius)};
      const auto path =
          temp.temp_dir / ("surface-process-spatial-cluster-" + name + "-3d.json");
      WriteJsonOnRoot(path, near_library);
      return path;
    };
    auto Preflight = [&](const fs::path &path, const std::string &name)
    {
      const auto manifest_path =
          temp.temp_dir / ("surface-response-requirements-" + name + ".json");
      if (Mpi::Root(Mpi::World()))
      {
        fs::remove(temp.temp_dir / "surface-response-patches.csv");
      }
      Mpi::Barrier(Mpi::World());
      IoData near_iodata = MakeIoData(path);
      WriteSurfaceResponseRequirements(near_iodata, mesh, manifest_path.string());
      Mpi::Barrier(Mpi::World());
      return ReadJsonAfterBarrier(manifest_path);
    };
    auto ModelRows = [&](const std::string &name)
    {
      const auto rows = ReadPatchRows(temp.temp_dir / "surface-response-patches.csv");
      return std::count_if(rows.begin(), rows.end(),
                           [&](const auto &row) { return row[3] == name; });
    };
    const json near = Preflight(NearContextLibrary(5.0e-2, true, "near-context-permuted"),
                                "near-context-permuted");
    const auto &near_features = near["Identification"]["Features"];
    const auto near_feature = std::find_if(near_features.begin(), near_features.end(),
                                           IsContextOnlyConductorCluster);
    REQUIRE(near_feature != near_features.end());
    CHECK((*near_feature)["Hash"] == (*feature)["Hash"]);
    REQUIRE((*near_feature)["Match"]["Status"] == "GeometricNearMatched");
    CHECK((*near_feature)["Match"]["Model"] == "near-context-permuted");
    const auto &record = (*near_feature)["Match"]["GeometricNearMatch"];
    CHECK_THAT(record["ContextDeviationOverR"].get<double>(), WithinAbs(5.0e-2, 1.0e-6));
    CHECK(record["ContextPermuted"] == true);
    CHECK(record["ClaimsDeviationOverR"].get<double>() <= 1.0e-9);
    CHECK(near["Summary"]["GeometricNearMatched"]["Count"].get<int>() >= 1);
    CHECK(!near["Identification"]["Diagnostics"].contains("GeometryFallbacks"));
    CHECK(ModelRows("near-context-permuted") >= 1);
    const auto &near_requirements = near["Requirements"];
    const auto near_requirement = std::find_if(
        near_requirements.begin(), near_requirements.end(),
        [&](const json &entry) { return entry["Hash"] == (*feature)["Hash"]; });
    REQUIRE(near_requirement != near_requirements.end());
    CHECK((*near_requirement)["Status"] == "GeometricNearMatched");
    CHECK((*near_requirement)["SelectedModels"][0]["Name"] == "near-context-permuted");
    const json outside = Preflight(NearContextLibrary(1.1e-1, false, "far-context-outside"),
                                   "far-context-outside");
    const auto &outside_features = outside["Identification"]["Features"];
    const auto outside_feature = std::find_if(
        outside_features.begin(), outside_features.end(), IsContextOnlyConductorCluster);
    REQUIRE(outside_feature != outside_features.end());
    CHECK((*outside_feature)["Match"]["Status"] == "Missing");
    CHECK(!(*outside_feature)["Match"].contains("GeometricNearMatch"));
    CHECK(outside["Summary"]["GeometricNearMatched"]["Count"] == 0);
    CHECK(!outside["Identification"]["Diagnostics"].contains("GeometryFallbacks"));
    CHECK(ModelRows("far-context-outside") == 0);
  }

  // Fail closed, kept: the same model whose plan view names conductor 1 only (two
  // references for one canonical conductor).
  {
    auto one_conductor_plan_library = library;
    one_conductor_plan_library["Name"] =
        "unit-test-process-spatial-cluster-context-only-one-plan-3d";
    one_conductor_plan_library["Models"].back()["PlanViewBoundary"] = {
        PlanViewComponent(1, claimed_pieces, radius)};
    const auto path =
        temp.temp_dir / "surface-process-spatial-cluster-context-only-one-plan-3d.json";
    WriteJsonOnRoot(path, one_conductor_plan_library);
    IoData refused_iodata = MakeIoData(path);
    CHECK_THROWS_WITH(SurfaceResponseOperator(refused_iodata, space),
                      ContainsSubstring(std::string(count_message)));
  }

  // Fail closed, new: a single Reference with a plan of two conductors (one reference for
  // two canonical conductors; before this fix the count check passed and the model was
  // refused later, by a different message: the basis-size mismatch of the fixture matrices'
  // one conductor state against conductor_state_count 0 of the Reference branch, or
  // "ConductorReferences must contain at least two points" -- the count message now names
  // the defect at the library read).
  {
    auto single_reference_library = library;
    single_reference_library["Name"] =
        "unit-test-process-spatial-cluster-context-only-single-reference-3d";
    auto &single_reference_model = single_reference_library["Models"].back();
    single_reference_model["Reference"] = reference_1;
    single_reference_model.erase("ConductorReferences");
    single_reference_model.erase("OpenContourPaths");
    single_reference_model["ContourGroups"] = {4};
    const auto path =
        temp.temp_dir /
        "surface-process-spatial-cluster-context-only-single-reference-3d.json";
    WriteJsonOnRoot(path, single_reference_library);
    IoData refused_iodata = MakeIoData(path);
    CHECK_THROWS_WITH(SurfaceResponseOperator(refused_iodata, space),
                      ContainsSubstring(std::string(count_message)));
  }

  // Fail closed, kept: a plan-view label beyond the references (conductor 3 with no label
  // 2: three canonical conductors for two references).
  {
    auto skipped_label_library = library;
    skipped_label_library["Name"] =
        "unit-test-process-spatial-cluster-context-only-skipped-label-3d";
    skipped_label_library["Models"].back()["PlanViewBoundary"] = {
        PlanViewComponent(1, claimed_pieces, radius),
        PlanViewComponent(3, foreign_pieces, radius)};
    const auto path = temp.temp_dir /
                      "surface-process-spatial-cluster-context-only-skipped-label-3d.json";
    WriteJsonOnRoot(path, skipped_label_library);
    IoData refused_iodata = MakeIoData(path);
    CHECK_THROWS_WITH(SurfaceResponseOperator(refused_iodata, space),
                      ContainsSubstring(std::string(count_message)));
  }
#endif
}

}  // namespace palace
