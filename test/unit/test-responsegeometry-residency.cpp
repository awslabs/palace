// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "surfaceresponse-fixtures.hpp"

#include <cstddef>
#include <memory>
#include <utility>
#include <vector>
#include <catch2/catch_test_macros.hpp>
#include <nlohmann/json.hpp>
#include "fem/mesh.hpp"
#include "linalg/vector.hpp"
#include "models/laplaceoperator.hpp"
#include "models/surfaceresponseoperator.hpp"
#include "utils/communication.hpp"
#include "utils/iodata.hpp"

using namespace palace;
using json = nlohmann::json;

namespace
{

using ResponsePatchData = config::ElectrostaticSolverData::ResponseCorrectionPatchData;

// Every field of ResponseCorrectionPatchData set to a non-default value, every vector with
// two entries (configfile.hpp: a field added there and not to PackPatch / UnpackPatch shows
// here as a default on the unpacked side).
ResponsePatchData FullyPopulatedPatch()
{
  ResponsePatchData patch;
  patch.model = 7;
  patch.origin = {1.25, -2.5, 3.75};
  patch.axis_u = {0.6, 0.8, 0.0};
  patch.axis_v = {-0.8, 0.6, 0.0};
  patch.axis_w = {0.0, 0.0, -1.0};
  patch.conductor_references = {{0.1, 0.2, 0.3}, {-0.4, 0.5, -0.6}};
  patch.conductor_probe_edges = {{0.11, -1, 2, 3}, {-0.22, 1, 1, 2}};
  patch.conductor_probes = {{{1.0, 2.0, 3.0}, 1, 1}, {{-1.0, -2.0, -3.0}, 0, 0}};
  patch.weight = 0.375;
  patch.domain_weight = 0.625;
  patch.longitudinal_cell = {-0.03125, 0.046875};
  patch.interpolation_group = 5;
  patch.maxwell_conductor_anchors = {{9.0, 8.0, 7.0}, {6.0, 5.0, 4.0}};
  patch.maxwell_reference_is_pec = false;
  auto &provenance = patch.provenance;
  provenance.feature = 13;
  provenance.segment = 17;
  provenance.s0 = 1.0e-7;
  provenance.s1 = 2.5e-7;
  provenance.stretch = 3;
  provenance.edge_offset = -0.015625;
  provenance.own_cell = {{{1.1, 1.2, 1.3}, {2.1, 2.2, 2.3}}};
  provenance.has_own_cell = true;
  provenance.quadrature_weight = 0.9;
  provenance.model_weight = 0.45;
  provenance.side_factor = 0.5;
  provenance.coupon_depth = 1.9e-6;
  provenance.claims = {{4, {0.01, 0.02, 0.03}, {0.04, 0.05, 0.06}},
                       {5, {0.07, 0.08, 0.09}, {0.10, 0.11, 0.12}}};
  provenance.support_box = {-1.0, 2.0, -3.0, 4.0};
  provenance.has_support_box = true;
  provenance.chain = {{1.0, 2.0, 3.0, 4.0}, {5.0, 6.0, 7.0, 8.0}};
  provenance.foreign = {{-1.0, -2.0, -3.0, -4.0}, {-5.0, -6.0, -7.0, -8.0}};
  provenance.raw_claims = {{6, {0.13, 0.14, 0.15}, {0.16, 0.17, 0.18}},
                           {7, {0.19, 0.20, 0.21}, {0.22, 0.23, 0.24}}};
  provenance.mirror_formed_topology = 2;
  provenance.mirror_reach = 3.8e-6;
  provenance.mirror_edge_weights = {0.25, 0.75};
  provenance.probe_snaps = {{11, 12, 13, {1.0e-9, -2.0e-9, 3.0e-9}, 3.7e-9},
                            {14, 15, 16, {-4.0e-9, 5.0e-9, -6.0e-9}, 8.8e-9}};
  provenance.probe_snaps_computed = false;
  return patch;
}

void CheckSamePatch(const ResponsePatchData &a, const ResponsePatchData &b)
{
  CHECK(a.model == b.model);
  CHECK(a.origin == b.origin);
  CHECK(a.axis_u == b.axis_u);
  CHECK(a.axis_v == b.axis_v);
  CHECK(a.axis_w == b.axis_w);
  CHECK(a.conductor_references == b.conductor_references);
  REQUIRE(a.conductor_probe_edges.size() == b.conductor_probe_edges.size());
  for (std::size_t i = 0; i < a.conductor_probe_edges.size(); i++)
  {
    CHECK(a.conductor_probe_edges[i].offset == b.conductor_probe_edges[i].offset);
    CHECK(a.conductor_probe_edges[i].gap_direction ==
          b.conductor_probe_edges[i].gap_direction);
    CHECK(a.conductor_probe_edges[i].conductor == b.conductor_probe_edges[i].conductor);
    CHECK(a.conductor_probe_edges[i].piece == b.conductor_probe_edges[i].piece);
  }
  REQUIRE(a.conductor_probes.size() == b.conductor_probes.size());
  for (std::size_t i = 0; i < a.conductor_probes.size(); i++)
  {
    CHECK(a.conductor_probes[i].point == b.conductor_probes[i].point);
    CHECK(a.conductor_probes[i].edge == b.conductor_probes[i].edge);
    CHECK(a.conductor_probes[i].end_face == b.conductor_probes[i].end_face);
  }
  CHECK(a.weight == b.weight);
  CHECK(a.domain_weight == b.domain_weight);
  CHECK(a.longitudinal_cell == b.longitudinal_cell);
  CHECK(a.interpolation_group == b.interpolation_group);
  CHECK(a.maxwell_conductor_anchors == b.maxwell_conductor_anchors);
  CHECK(a.maxwell_reference_is_pec == b.maxwell_reference_is_pec);
  const auto &p = a.provenance, &q = b.provenance;
  CHECK(p.feature == q.feature);
  CHECK(p.segment == q.segment);
  CHECK(p.s0 == q.s0);
  CHECK(p.s1 == q.s1);
  CHECK(p.stretch == q.stretch);
  CHECK(p.edge_offset == q.edge_offset);
  CHECK(p.own_cell == q.own_cell);
  CHECK(p.has_own_cell == q.has_own_cell);
  CHECK(p.quadrature_weight == q.quadrature_weight);
  CHECK(p.model_weight == q.model_weight);
  CHECK(p.side_factor == q.side_factor);
  CHECK(p.coupon_depth == q.coupon_depth);
  for (const auto &pair :
       {std::make_pair(&p.claims, &q.claims), std::make_pair(&p.raw_claims, &q.raw_claims)})
  {
    REQUIRE(pair.first->size() == pair.second->size());
    for (std::size_t i = 0; i < pair.first->size(); i++)
    {
      CHECK((*pair.first)[i].segment == (*pair.second)[i].segment);
      CHECK((*pair.first)[i].p0 == (*pair.second)[i].p0);
      CHECK((*pair.first)[i].p1 == (*pair.second)[i].p1);
    }
  }
  CHECK(p.support_box == q.support_box);
  CHECK(p.has_support_box == q.has_support_box);
  CHECK(p.chain == q.chain);
  CHECK(p.foreign == q.foreign);
  CHECK(p.mirror_formed_topology == q.mirror_formed_topology);
  CHECK(p.mirror_reach == q.mirror_reach);
  CHECK(p.mirror_edge_weights == q.mirror_edge_weights);
  REQUIRE(p.probe_snaps.size() == q.probe_snaps.size());
  for (std::size_t i = 0; i < p.probe_snaps.size(); i++)
  {
    CHECK(p.probe_snaps[i].vertex == q.probe_snaps[i].vertex);
    CHECK(p.probe_snaps[i].segment == q.probe_snaps[i].segment);
    CHECK(p.probe_snaps[i].arc == q.probe_snaps[i].arc);
    CHECK(p.probe_snaps[i].displacement == q.probe_snaps[i].displacement);
    CHECK(p.probe_snaps[i].distance == q.probe_snaps[i].distance);
  }
  CHECK(p.probe_snaps_computed == q.probe_snaps_computed);
}

}  // namespace

// Lane fix-response-memory (decisions 1439 / 1446 (R-1: C6 approved); DESIGN A1.0 / A1.4,
// commit C6): the placed-patch list's residence. (1) PackPatch / UnpackPatch round-trip a
// fully populated patch field for field. (2) After a construction the shared geometry keeps
// this rank's LOCAL patches only (the round-robin count), the global count kept; a second
// construction on the same geometry (an AMR iteration's) re-assembles the full list and
// builds the SAME operator: patch count, root patch assignments, the uncovered / domain-
// boundary portions and a response evaluation bit-equal to the first's.
TEST_CASE_METHOD(test::SurfaceResponseFiles,
                 "SurfaceResponseGeometry keeps the local patches between constructions "
                 "and re-assembles the former list",
                 "[surfaceresponseoperator][memfix][3d][legacy][Serial][Parallel]")
{
#if !defined(MFEM_USE_GSLIB)
  SKIP("SurfaceResponseOperator requires MFEM_USE_GSLIB");
#else
  MPI_Comm comm = Mpi::World();
  const int rank = Mpi::Rank(comm), size = Mpi::Size(comm);

  SECTION("PackPatch / UnpackPatch round-trip every field")
  {
    const ResponsePatchData original = FullyPopulatedPatch();
    std::vector<double> packed;
    SurfaceResponseGeometry::PackPatch(packed, 42, original);
    SurfaceResponseGeometry::PackPatch(packed, 7, ResponsePatchData());
    std::size_t index = 0;
    ResponsePatchData unpacked;
    std::size_t offset = SurfaceResponseGeometry::UnpackPatch(packed, 0, index, unpacked);
    CHECK(index == 42);
    CheckSamePatch(unpacked, original);
    ResponsePatchData defaults;
    offset = SurfaceResponseGeometry::UnpackPatch(packed, offset, index, defaults);
    CHECK(index == 7);
    CheckSamePatch(defaults, ResponsePatchData());
    CHECK(offset == packed.size());
  }

  SECTION("the shared geometry across two constructions")
  {
    IoData iodata(IslandConfig(), false);
    iodata.boundaries.cracked_attributes.insert(9);
    std::vector<std::unique_ptr<Mesh>> meshes;
    meshes.push_back(std::make_unique<Mesh>(MakeIslandMesh()));
    LaplaceOperator laplace(iodata, meshes);
    std::shared_ptr<const SurfaceResponseGeometry> geometry;
    SurfaceResponseOperator first(iodata, laplace, &geometry);
    REQUIRE(geometry);
    REQUIRE(first.GetPatchCount() > 1);
    // The residence after the construction: local patches only, the global count kept.
    CHECK(geometry->PatchesLocal());
    CHECK(geometry->GlobalPatchCount() == static_cast<std::size_t>(first.GetPatchCount()));
    std::size_t local = 0;
    for (int p = 0; p < first.GetPatchCount(); p++)
    {
      local += (p % size == rank) ? 1 : 0;
    }
    CHECK(geometry->RetainedPatchCount() == local);
    if (size > 1)
    {
      CHECK(geometry->RetainedPatchCount() < geometry->GlobalPatchCount());
    }
    // The second construction (an AMR iteration's) re-assembles the list and releases it
    // again; the operator it builds is the first one.
    SurfaceResponseOperator second(iodata, laplace, &geometry);
    CHECK(geometry->PatchesLocal());
    CHECK(geometry->RetainedPatchCount() == local);
    REQUIRE(second.GetPatchCount() == first.GetPatchCount());
    CHECK(second.GetBasisSize() == first.GetBasisSize());
    CHECK(second.GetPatchWeight() == first.GetPatchWeight());
    const auto &a = first.GetPatchAssignments(), &b = second.GetPatchAssignments();
    REQUIRE(a.size() == b.size());  // root: every patch; the other ranks: none
    for (std::size_t p = 0; p < a.size(); p++)
    {
      CHECK(a[p].model == b[p].model);
      CHECK(a[p].global_index == b[p].global_index);
      CHECK(a[p].origin == b[p].origin);
      CHECK(a[p].axis_u == b[p].axis_u);
      CHECK(a[p].axis_v == b[p].axis_v);
      CHECK(a[p].axis_w == b[p].axis_w);
      CHECK(a[p].weight == b[p].weight);
    }
    CHECK(first.GetUncoveredPortions().size() == second.GetUncoveredPortions().size());
    CHECK(first.GetDomainBoundaryPortions().size() ==
          second.GetDomainBoundaryPortions().size());
    const int n = laplace.GetH1Space().GetTrueVSize();
    Vector x(n);
    for (int i = 0; i < n; i++)
    {
      x(i) = 0.21 + 0.017 * (i + 1) * (rank + 1);
    }
    const auto response_first = first.GetElectrostaticResponse(x, true, true);
    const auto response_second = second.GetElectrostaticResponse(x, true, true);
    CHECK(response_second.domain_correction == response_first.domain_correction);
    CHECK(response_second.fabricated_surface_energy ==
          response_first.fabricated_surface_energy);
    REQUIRE(response_second.patch_contributions.size() ==
            response_first.patch_contributions.size());
    for (std::size_t p = 0; p < response_first.patch_contributions.size(); p++)
    {
      CHECK(response_second.patch_contributions[p].patch ==
            response_first.patch_contributions[p].patch);
      CHECK(response_second.patch_contributions[p].weight ==
            response_first.patch_contributions[p].weight);
      CHECK(response_second.patch_contributions[p].domain_correction ==
            response_first.patch_contributions[p].domain_correction);
      CHECK(response_second.patch_contributions[p].fabricated_surface_energy ==
            response_first.patch_contributions[p].fabricated_surface_energy);
    }
  }
#endif
}
