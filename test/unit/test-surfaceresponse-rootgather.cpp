// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "surfaceresponse-fixtures.hpp"

#include <cstddef>
#include <memory>
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

void CheckSamePatchContributions(
    const std::vector<SurfaceResponseOperator::PatchContribution> &a,
    const std::vector<SurfaceResponseOperator::PatchContribution> &b)
{
  REQUIRE(a.size() == b.size());
  for (std::size_t i = 0; i < a.size(); i++)
  {
    CHECK(a[i].patch == b[i].patch);
    CHECK(a[i].model == b[i].model);
    CHECK(a[i].feature == b[i].feature);
    CHECK(a[i].weight == b[i].weight);
    CHECK(a[i].cell == b[i].cell);
    CHECK(a[i].domain_correction == b[i].domain_correction);
    CHECK(a[i].domain_correction_fixed_flux == b[i].domain_correction_fixed_flux);
    CHECK(a[i].fabricated_surface_energy == b[i].fabricated_surface_energy);
    CHECK(a[i].fabricated_surface_energy_fixed_flux ==
          b[i].fabricated_surface_energy_fixed_flux);
  }
}

void CheckSamePatchTraces(const std::vector<SurfaceResponseOperator::PatchTrace> &a,
                          const std::vector<SurfaceResponseOperator::PatchTrace> &b)
{
  REQUIRE(a.size() == b.size());
  for (std::size_t i = 0; i < a.size(); i++)
  {
    CHECK(a[i].patch == b[i].patch);
    CHECK(a[i].model == b[i].model);
    CHECK(a[i].contour_size == b[i].contour_size);
    CHECK(a[i].coefficients == b[i].coefficients);
  }
}

// The reduced scalars of a response (every rank, both residencies): bit-equal.
void CheckSameScalars(const SurfaceResponseOperator::ElectrostaticResponse &a,
                      const SurfaceResponseOperator::ElectrostaticResponse &b)
{
  CHECK(a.domain_correction == b.domain_correction);
  CHECK(a.domain_correction_fixed_flux == b.domain_correction_fixed_flux);
  CHECK(a.fabricated_surface_energy == b.fabricated_surface_energy);
  CHECK(a.fabricated_surface_energy_fixed_flux == b.fabricated_surface_energy_fixed_flux);
  CHECK(a.trace_closure_spread == b.trace_closure_spread);
  CHECK(a.maximum_trace_closure_spread == b.maximum_trace_closure_spread);
  REQUIRE(a.model_contributions.size() == b.model_contributions.size());
  for (std::size_t m = 0; m < a.model_contributions.size(); m++)
  {
    CHECK(a.model_contributions[m].model == b.model_contributions[m].model);
    CHECK(a.model_contributions[m].patch_count == b.model_contributions[m].patch_count);
    CHECK(a.model_contributions[m].patch_weight == b.model_contributions[m].patch_weight);
    CHECK(a.model_contributions[m].domain_correction ==
          b.model_contributions[m].domain_correction);
    CHECK(a.model_contributions[m].fabricated_surface_energy ==
          b.model_contributions[m].fabricated_surface_energy);
  }
}

}  // namespace

// Lane fix-response-memory (decisions 1439 / 1446 / 1447; DESIGN s2.2 / s3 (M1), commit
// C4): Mpi::GatherLarge gives the root the rank-order concatenation Mpi::AllgatherLarge
// gives every rank (ragged lengths, an empty rank); and the operator's ROOT_ONLY
// patch-record residency (the drivers' form) gives the root the REPLICATED gathers'
// per-patch contributions and traces to the bit, the other ranks empty lists, with every
// reduced scalar (the energies, the model contributions) identical on every rank in both
// modes.
TEST_CASE_METHOD(test::SurfaceResponseFiles,
                 "SurfaceResponseOperator ROOT_ONLY patch records equal the REPLICATED "
                 "gathers on the root and are empty elsewhere",
                 "[surfaceresponseoperator][memfix][2d][patch][Serial][Parallel]")
{
#if !defined(MFEM_USE_GSLIB)
  SKIP("SurfaceResponseOperator requires MFEM_USE_GSLIB");
#else
  MPI_Comm comm = Mpi::World();
  const int rank = Mpi::Rank(comm), size = Mpi::Size(comm);

  SECTION("GatherLarge is AllgatherLarge's rank-order concatenation on the root")
  {
    // Rank r holds 3 + 2 r values, rank 1 none (the empty rank; ragged lengths).
    std::vector<double> local;
    if (rank != 1)
    {
      for (int i = 0; i < 3 + 2 * rank; i++)
      {
        local.push_back(1000.0 * rank + i + 0.25);
      }
    }
    const auto everywhere = Mpi::AllgatherLarge(local, comm);
    const auto root_only = Mpi::GatherLarge(local, 0, comm);
    std::size_t expected = 0;
    for (int r = 0; r < size; r++)
    {
      expected += r == 1 ? 0 : static_cast<std::size_t>(3 + 2 * r);
    }
    REQUIRE(everywhere.size() == expected);
    if (rank == 0)
    {
      CHECK(root_only == everywhere);
    }
    else
    {
      CHECK(root_only.empty());
    }
  }

  SECTION("the operator's per-patch records by residency")
  {
    // Five explicit 2D patches of one model on the unit square (round-robin over the ranks:
    // 3 + 2 under two ranks), the explicit-patches fixture of test-surfaceresponseoperator.
    json patches = json::array();
    for (const auto &origin : {std::array<double, 2>{0.87, 0.35},
                               {0.62, 0.35},
                               {0.37, 0.35},
                               {0.12, 0.35},
                               {0.50, 0.60}})
    {
      patches.push_back({{"Model", 1},
                         {"Origin", {origin[0], origin[1], 0.0}},
                         {"AxisU", {1.0, 0.0, 0.0}},
                         {"AxisV", {0.0, 1.0, 0.0}},
                         {"Reference", {0.0, 0.0, 0.0}}});
    }
    json config = {
        {"Problem", {{"Type", "Electrostatic"}, {"Output", temp.temp_dir.string()}}},
        {"Model", {{"Mesh", "unused.msh"}}},
        {"Domains", {{"Materials", {{{"Attributes", {1}}}}}}},
        {"Boundaries",
         {{"Ground", {{"Attributes", {1, 3, 4}}}},
          {"Terminal", {{{"Index", 1}, {"Attributes", {2}}}}}}},
        {"Solver",
         {{"Order", 1},
          {"Electrostatic",
           {{"ResponseCorrection",
             {{"Models",
               {{{"Index", 1},
                 {"FabricatedMatrix", fabricated_path.string()},
                 {"ThinMatrix", thin_path.string()},
                 {"FabricatedSurfaceMatrix", fabricated_surface_path.string()},
                 {"ThinSurfaceMatrix", thin_surface_path.string()},
                 {"BasisPoints", points_path.string()},
                 {"Interfaces", {{{"Target", 4}, {"Coupon", 1}}}}}}},
              {"Patches", patches}}}}}}}};
    IoData iodata(config, false);
    mfem::Mesh serial_mesh =
        mfem::Mesh::MakeCartesian2D(8, 8, mfem::Element::TRIANGLE, false, 1.0, 1.0);
    while (serial_mesh.GetNE() < size)
    {
      serial_mesh.UniformRefinement();
    }
    std::vector<std::unique_ptr<Mesh>> meshes;
    meshes.push_back(
        std::make_unique<Mesh>(std::make_unique<mfem::ParMesh>(comm, serial_mesh)));
    LaplaceOperator laplace_op(iodata, meshes);
    SurfaceResponseOperator response(iodata, laplace_op);
    REQUIRE(response.GetPatchCount() == 5);
    CHECK(response.GetPatchRecordResidency() ==
          SurfaceResponseOperator::PatchRecordResidency::REPLICATED);
    const int n = laplace_op.GetH1Space().GetTrueVSize();
    Vector x(n);
    for (int i = 0; i < n; i++)
    {
      x(i) = 0.17 + 0.03 * (i + 1) * (rank + 1);
    }

    // The replicated reference (the former form, the default).
    const auto reference = response.GetElectrostaticResponse(x, true, true);
    const auto reference_traces = response.GetPatchTraces(x);
    REQUIRE(reference.patch_contributions.size() == 5);
    REQUIRE(reference_traces.size() == 5);
    for (std::size_t p = 1; p < 5; p++)
    {
      CHECK(reference.patch_contributions[p - 1].patch <
            reference.patch_contributions[p].patch);
    }

    // ROOT_ONLY (the drivers' form).
    response.SetPatchRecordResidency(
        SurfaceResponseOperator::PatchRecordResidency::ROOT_ONLY);
    const auto root_only = response.GetElectrostaticResponse(x, true, true);
    const auto root_only_traces = response.GetPatchTraces(x);
    CheckSameScalars(root_only, reference);
    if (rank == 0)
    {
      CheckSamePatchContributions(root_only.patch_contributions,
                                  reference.patch_contributions);
      CheckSamePatchTraces(root_only_traces, reference_traces);
    }
    else
    {
      CHECK(root_only.patch_contributions.empty());
      CHECK(root_only_traces.empty());
    }

    // Back to REPLICATED: the former lists again on every rank (the mode is a residency,
    // not a state).
    response.SetPatchRecordResidency(
        SurfaceResponseOperator::PatchRecordResidency::REPLICATED);
    const auto again = response.GetElectrostaticResponse(x, true, true);
    CheckSameScalars(again, reference);
    CheckSamePatchContributions(again.patch_contributions, reference.patch_contributions);
    CheckSamePatchTraces(response.GetPatchTraces(x), reference_traces);
  }
#endif
}
