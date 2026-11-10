// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "surfaceresponse-fixtures.hpp"

#include <array>
#include <cstddef>
#include <limits>
#include <set>
#include <string>
#include <vector>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>
#include "models/surfaceresponseoperator.hpp"
#include "utils/communication.hpp"

using namespace palace;
using namespace Catch::Matchers;

namespace
{

using ResponsePatchData = config::ElectrostaticSolverData::ResponseCorrectionPatchData;
using MirrorPlaneData =
    config::ElectrostaticSolverData::ResponseCorrectionData::MirrorPlaneData;
using GeometryFallbackData =
    config::ElectrostaticSolverData::ResponseCorrectionData::GeometryFallbackData;

void CheckSameExclusions(const DomainBoundaryExclusions &a,
                         const DomainBoundaryExclusions &b)
{
  CHECK(a.tested_patches == b.tested_patches);
  CHECK(a.tested_points == b.tested_points);
  CHECK(a.reflected_points == b.reflected_points);
  REQUIRE(a.patches.size() == b.patches.size());
  for (std::size_t i = 0; i < a.patches.size(); i++)
  {
    const auto &x = a.patches[i], &y = b.patches[i];
    CHECK(x.patch == y.patch);
    CHECK(x.tested_points == y.tested_points);
    CHECK(x.outside_points == y.outside_points);
    CHECK(x.nearest_outside_point == y.nearest_outside_point);
    CHECK(x.nearest_distance == y.nearest_distance);
    CHECK(x.reason == y.reason);
    CHECK(x.unmerged_topology == y.unmerged_topology);
    CHECK(x.unmerged_own == y.unmerged_own);
  }
  REQUIRE(a.mirrored.size() == b.mirrored.size());
  for (std::size_t i = 0; i < a.mirrored.size(); i++)
  {
    CHECK(a.mirrored[i].patch == b.mirrored[i].patch);
    CHECK(a.mirrored[i].tested_points == b.mirrored[i].tested_points);
    CHECK(a.mirrored[i].outside_points == b.mirrored[i].outside_points);
  }
  CHECK(a.unmerged_reach.size() == b.unmerged_reach.size());
}

void CheckSameLedger(const std::vector<GeometryFallbackData> &a,
                     const std::vector<GeometryFallbackData> &b)
{
  REQUIRE(a.size() == b.size());
  for (std::size_t i = 0; i < a.size(); i++)
  {
    CHECK(a[i].cls == b[i].cls);
    CHECK(a[i].site == b[i].site);
    CHECK(a[i].reason == b[i].reason);
    CHECK(a[i].what == b[i].what);
    CHECK(a[i].patches == b[i].patches);
    CHECK(a[i].models == b[i].models);
    CHECK(a[i].features == b[i].features);
    CHECK(a[i].length == b[i].length);
    CHECK(a[i].example_origin == b[i].example_origin);
  }
}

}  // namespace

// Lane fix-response-memory (decisions 1439 / 1446; DESIGN s2.5 / s3 I-7, commit C3): the
// containment test in candidate blocks of at most `block_points` tested points gives the
// SAME DomainBoundaryExclusions (patches, reasons, nearest points / distances, the
// Mirrored list, the reflected count), the SAME patch weights and the SAME
// geometry-fallback records as one block holding every point, for every block size down to
// ONE point (every candidate its own block, the strongest form), on a configuration with
// applied, cut (DomainBoundary MirrorOff), Mirrored (a Natural plane on the x = 1 face),
// misplaced (UnlocatedPatch) and all-points-outside patches of two models with different
// point counts.
TEST_CASE_METHOD(test::SurfaceResponseFiles,
                 "FindDomainBoundaryExclusions in candidate blocks equals the whole-list "
                 "form for every block size",
                 "[surfaceresponseoperator][domainboundary][memfix][3d][Serial][Parallel]")
{
#if !defined(MFEM_USE_GSLIB)
  SKIP("SurfaceResponseOperator requires MFEM_USE_GSLIB");
#else
  constexpr double R = 0.2;
  auto mesh = MakeIslandMesh();
  // Model 0: a 4-point cross-section (3 sections for a longitudinal cell: 3 x (4 + 1
  // reference) points); model 1: a spatial 8-point box (1 section: 8 + 1 points).
  const std::vector<std::array<double, 3>> edge_basis = {
      {-0.05, -0.05, 0.0}, {0.05, -0.05, 0.0}, {0.05, 0.05, 0.0}, {-0.05, 0.05, 0.0}};
  const std::vector<std::array<double, 3>> spatial_basis = {
      {-0.05, -0.05, -0.05}, {0.05, -0.05, -0.05}, {0.05, 0.05, -0.05},
      {-0.05, 0.05, -0.05},  {-0.05, -0.05, 0.05}, {0.05, -0.05, 0.05},
      {0.05, 0.05, 0.05},    {-0.05, 0.05, 0.05}};
  auto BasisPoints = [&](int model_idx) -> const std::vector<std::array<double, 3>> *
  { return model_idx == 1 ? &spatial_basis : &edge_basis; };
  auto Spatial = [](int model_idx) { return model_idx == 1; };
  auto Name = [](int model_idx) { return std::string(model_idx == 1 ? "box" : "edge"); };
  // A cell in the lead's frame (u from the metal into the gap = -x, w along the edge) with
  // its metal-edge reference at the origin.
  auto Patch = [&](const std::array<double, 3> &origin, int model, int feature,
                   bool with_reference = true)
  {
    ResponsePatchData patch;
    patch.model = model;
    patch.origin = origin;
    patch.axis_u = {-1.0, 0.0, 0.0};
    patch.axis_v = {0.0, 1.0, 0.0};
    patch.axis_w = {0.0, 0.0, 1.0};
    patch.longitudinal_cell =
        model == 1 ? std::array<double, 2>{0.0, 0.0} : std::array<double, 2>{-0.03, 0.03};
    patch.provenance.feature = feature;
    patch.provenance.coupon_depth = model == 1 ? 0.0 : R;
    // The default configuration carries the origin reference; a cell without one tests
    // the all-points-outside path.
    if (!with_reference)
    {
      patch.conductor_references.clear();
    }
    return patch;
  };
  // The Natural mirror plane on the x = 1 face (outward normal +x, offset 1) with a band
  // wide enough for the gap-side points at x = 1.05: reflected into the mesh -> Mirrored.
  MirrorPlaneData plane;
  plane.attribute = 7;
  plane.normal = {1.0, 0.0, 0.0};
  plane.offset = 1.0;
  plane.status = "Natural";
  plane.box_min = {1.0, 0.0, 0.0};
  plane.box_max = {1.0, 1.0, 1.0};
  const std::vector<MirrorPlaneData> planes = {plane};
  const double band = 0.1;

  // The configuration: applied cells of both models, the cut cells on the x = 1 face
  // (Mirrored with the plane, DomainBoundary MirrorOff without), a misplaced reference
  // (1e-2 R off the face), an all-points-outside cell without a reference, interleaved so
  // that the block boundaries fall between every class.
  const std::vector<ResponsePatchData> configuration = {
      Patch({0.25, 0.5, 0.30}, 0, 3), Patch({1.0, 0.5, 0.30}, 0, 7),
      Patch({0.25, 0.5, 0.50}, 1, 4), Patch({0.25, 0.5, 0.70}, 0, 3),
      Patch({1.0, 0.5, 0.60}, 1, 8),  Patch({1.0 + 1.0e-2 * R, 0.5, 0.45}, 0, 9),
      Patch({0.60, 0.5, 0.30}, 1, 5), Patch({2.0, 0.5, 0.30}, 0, 11, false),
      Patch({1.0, 0.5, 0.80}, 0, 7),  Patch({0.25, 0.5, 0.90}, 0, 3),
  };
  const std::set<std::size_t> skipped = {3};  // one skipped candidate (not tested)
  const std::size_t all_points = std::numeric_limits<std::size_t>::max();
  std::vector<GeometryFallbackData> unused_entries;
  auto Run = [&](std::size_t block_points, bool mirror,
                 std::vector<ResponsePatchData> &patches,
                 std::vector<GeometryFallbackData> *entries)
  {
    GeometryFallbackLedger ledger(GeometryFailurePolicy::RECORD_RAW, false,
                                  entries ? *entries : unused_entries);
    return FindDomainBoundaryExclusions(
        *mesh, patches, BasisPoints, Spatial, Name, 1.0, R, skipped,
        mirror ? planes : std::vector<MirrorPlaneData>{}, mirror ? band : 0.0, {},
        entries ? &ledger : nullptr, block_points);
  };
  for (const bool mirror : {false, true})
  {
    CAPTURE(mirror);
    // The reference: one block holding every point, with the ledger.
    std::vector<GeometryFallbackData> reference_entries;
    std::vector<ResponsePatchData> reference_patches = configuration;
    const auto reference = Run(all_points, mirror, reference_patches, &reference_entries);
    REQUIRE(reference.tested_patches == 9);
    // 5 edge cells with a reference (3 sections x 5 points), one without (3 x 4), 3
    // spatial boxes (8 + 1).
    REQUIRE(reference.tested_points == 5 * 15 + 12 + 3 * 9);
    REQUIRE(!reference.patches.empty());
    if (mirror)
    {
      REQUIRE(!reference.mirrored.empty());
      CHECK(reference.reflected_points > 0);
    }
    else
    {
      CHECK(reference.mirrored.empty());
    }
    // UnlocatedPatch records: the all-points-outside cell, and the misplaced reference
    // unless the plane reflects it into the mesh.
    REQUIRE(reference_entries.size() == (mirror ? 1u : 2u));
    for (const auto &entry : reference_entries)
    {
      CHECK(entry.cls == "UnlocatedPatch");
    }
    for (const std::size_t block_points :
         {std::size_t(1), std::size_t(9), std::size_t(15), std::size_t(16), std::size_t(40),
          std::size_t(100), all_points})
    {
      CAPTURE(block_points);
      std::vector<GeometryFallbackData> entries;
      std::vector<ResponsePatchData> patches = configuration;
      const auto blocked = Run(block_points, mirror, patches, &entries);
      CheckSameExclusions(blocked, reference);
      CheckSameLedger(entries, reference_entries);
      for (std::size_t p = 0; p < patches.size(); p++)
      {
        CHECK(patches[p].weight == reference_patches[p].weight);
      }
      // Every phase is present whatever the blocking (the keys of the timing print).
      for (const char *phase :
           {"points", "locator", "find", "reduce", "reflect", "distances"})
      {
        CHECK(blocked.phase_times.count(phase) == 1);
      }
    }
  }
  // Without a ledger the fail-closed abort fires with the same patch named, at any
  // blocking (the misplaced reference is candidate 4 of 9: its own block or the only one).
  for (const std::size_t block_points : {std::size_t(1), all_points})
  {
    CAPTURE(block_points);
    std::vector<ResponsePatchData> patches = configuration;
    CHECK_THROWS_WITH(Run(block_points, false, patches, nullptr),
                      ContainsSubstring("misplaced or mis-scaled coupon") &&
                          ContainsSubstring("patch 5 (0-based") &&
                          ContainsSubstring("model edge") &&
                          ContainsSubstring("metal-edge reference"));
  }
  // The default block size is the design's 2^23 unless the environment sets it.
  CHECK(DomainBoundaryContainmentBlockPoints(Mpi::World()) >= 1);
#endif
}
