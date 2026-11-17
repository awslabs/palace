// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "fixtures.hpp"
#include "surfaceresponse-fixtures.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <fstream>
#include <memory>
#include <string>
#include <vector>
#include <mfem.hpp>
#include <catch2/catch_test_macros.hpp>
#include <nlohmann/json.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include "fem/gridfunction.hpp"
#include "fem/mesh.hpp"
#include "linalg/vector.hpp"
#include "models/spaceoperator.hpp"
#include "models/surfaceresponseoperator.hpp"
#include "utils/communication.hpp"
#include "utils/geodata.hpp"
#include "utils/iodata.hpp"

namespace palace
{

namespace fs = std::filesystem;

using json = nlohmann::json;
using Catch::Matchers::WithinAbs;

namespace
{

// The tetrahedral island mesh (MakeSerialIslandMesh: the metal sheet 9 in the plane
// y = 0.5, the island x, z in [0.25, 0.75] or, for `aperture`, its complement) read
// through Palace's own mesh pipeline, so that under the legacy crack key the PEC sheet is
// CRACKED as in a device run of record (AddInterfaceBdrElements: duplicated faces, the two
// sides' vertices moved apart by CrackDisplacementFactor x h_min) and under upstream's
// default (the key unset) it is NOT (P3-23). Under MPI the partition is the sheet plane
// itself (every element below y = 0.5 on the lower half of the ranks, every element above
// on the upper half), so the sheet faces are shared by two ranks' locators.
std::unique_ptr<mfem::ParMesh> ReadTetIslandMesh(IoData &iodata, bool aperture,
                                                 const fs::path &directory,
                                                 const std::string &tag)
{
  const auto mesh_path = directory / ("island-tet-" + tag + ".mesh");
  const auto partition_path = directory / ("island-tet-" + tag + ".partition");
  const auto comm = Mpi::World();
  const int size = Mpi::Size(comm);
  if (Mpi::Root(comm))
  {
    mfem::Mesh serial =
        test::SurfaceResponseFiles::MakeSerialIslandMesh(false, true, aperture);
    serial.Save(mesh_path.string(), 17);
  }
  Mpi::Barrier(comm);
  iodata.model.mesh = mesh_path.string();
  auto serial_mesh = mesh::Load(iodata, comm);
  if (size > 1)
  {
    if (serial_mesh)
    {
      const int lower_ranks = size / 2;
      std::ofstream output(partition_path);
      output << "number_of_elements " << serial_mesh->GetNE() << "\n"
             << "number_of_processors " << size << "\n";
      mfem::Vector center(3);
      for (int element = 0; element < serial_mesh->GetNE(); element++)
      {
        serial_mesh->GetElementCenter(element, center);
        output << (center[1] < 0.5 ? element % lower_ranks
                                   : lower_ranks + element % (size - lower_ranks))
               << "\n";
      }
    }
    Mpi::Barrier(comm);
    iodata.model.partitioning = partition_path.string();
  }
  return mesh::Partition(iodata, std::move(serial_mesh), comm);
}

// One ring reading of the Maxwell in-plane line cases: the patches of `model_name` (the
// corner basis rings), the largest error of a knot-to-knot trace difference against
// -E0 . (q - p) on the constant field E0, and the largest in-sheet ring trace with the PEC
// dofs zeroed.
struct RingReading
{
  std::string name;
  int patches = 0;
  double max_difference_error = 0.0;
  double max_zeroed_in_plane_trace = 0.0;
};

RingReading ReadMaxwellRing(const json &config, bool aperture, const fs::path &directory,
                            const std::string &tag, const std::string &model_name,
                            test::CrackPolicy crack_policy)
{
  mfem::VectorConstantCoefficient field_coefficient =
      test::SurfaceResponseFiles::ConstantFieldCoefficient();
  const mfem::Vector &constant_field = field_coefficient.GetVec();
  IoData iodata(config, false);
  std::vector<std::unique_ptr<Mesh>> meshes;
  meshes.push_back(
      std::make_unique<Mesh>(ReadTetIslandMesh(iodata, aperture, directory, tag)));
  // The premise of the form: the PEC sheet 9 is cracked under the legacy key only.
  REQUIRE(iodata.boundaries.cracked_attributes.count(9) ==
          (crack_policy == test::CrackPolicy::Legacy ? 1u : 0u));
  SpaceOperator space(iodata, meshes);
  SurfaceResponseOperator response(iodata, space);
  CHECK(response.GetPatchCount() > 0);

  GridFunction field(space.GetNDSpace(), true);
  field.Real().ProjectCoefficient(field_coefficient);
  field.Imag() = 0.0;
  const auto result = response.GetMaxwellResponse(field, 0.0);
  CHECK(result.loop_residual < 1.0e-10);
  CHECK_THAT(result.matched_length_fraction, WithinAbs(1.0, 1.0e-12));

  Vector field_true, zeroed_true;
  field.Real().GetTrueDofs(field_true);
  zeroed_true = field_true;
  zeroed_true.SetSubVector(space.GetNDDbcTDofLists().back(), 0.0);
  const auto traces = response.GetPatchTraces(field_true);
  const auto zeroed_traces = response.GetPatchTraces(zeroed_true);
  const auto names = response.GetModelNames();
  RingReading reading{model_name};
  if (Mpi::Root(Mpi::World()))
  {
    const auto &assignments = response.GetPatchAssignments();
    REQUIRE(traces.size() == zeroed_traces.size());
    for (std::size_t t = 0; t < traces.size(); t++)
    {
      const auto &entry = traces[t];
      if (names.at(entry.model) != model_name)
      {
        continue;
      }
      const auto assignment =
          std::find_if(assignments.begin(), assignments.end(),
                       [&](const auto &a) { return a.global_index == entry.patch; });
      REQUIRE(assignment != assignments.end());
      REQUIRE(entry.coefficients.size() == 12);
      reading.patches++;
      // The corner basis points: three rings of four knots at local z = -0.02, 0, 0.02.
      auto Device = [&](int knot)
      {
        const std::array<double, 2> xy[4] = {
            {0.02, 0.02}, {0.08, 0.02}, {0.08, 0.08}, {0.02, 0.08}};
        const double z = (knot / 4 - 1) * 0.02;
        std::array<double, 3> p{};
        for (int d = 0; d < 3; d++)
        {
          p[d] = assignment->origin[d] + xy[knot % 4][0] * assignment->axis_u[d] +
                 xy[knot % 4][1] * assignment->axis_v[d] + z * assignment->axis_w[d];
        }
        return p;
      };
      for (int ring = 0; ring < 3; ring++)
      {
        for (int i = 0; i < 4; i++)
        {
          const int knot = 4 * ring + i;
          const int next = 4 * ring + (i + 1) % 4;
          const auto p = Device(knot), q = Device(next);
          double expected = 0.0;
          for (int d = 0; d < 3; d++)
          {
            expected -= constant_field[d] * (q[d] - p[d]);
          }
          const double difference = entry.coefficients[next] - entry.coefficients[knot];
          reading.max_difference_error =
              std::max(reading.max_difference_error, std::abs(difference - expected));
        }
      }
      for (int i = 4; i < 8; i++)
      {
        reading.max_zeroed_in_plane_trace = std::max(
            reading.max_zeroed_in_plane_trace, std::abs(zeroed_traces[t].coefficients[i]));
      }
    }
  }
  return reading;
}

}  // namespace

TEST_CASE_METHOD(test::SurfaceResponseFiles,
                 "SurfaceResponseOperator Maxwell lines in cracked and interior faces",
                 "[surfaceresponseoperator][3d][maxwell][corner][Serial][Parallel]")
{
#if !defined(MFEM_USE_GSLIB)
  SKIP("SurfaceResponseOperator requires MFEM_USE_GSLIB");
#else
  // Decisions 829 / 831 (lane inplane): the corner models' z = 0 ring lies in the sheet
  // plane. On the convex island that ring (and the anchor line to it) lies ON the PEC
  // sheet, i.e. inside the crack Palace opens between the sheet's two faces; on the
  // aperture the concave corners' ring lies in the gap, in uncracked interior faces of the
  // same plane. Both lines must be covered and integrated exactly once: on the constant
  // field E0 every knot trace is -E0 . (p - a) from the patch's anchor a, so the trace
  // DIFFERENCE between two knots of one ring is -E0 . (p_i - p_j) whatever the anchor, and
  // every closed ring has a zero loop residual. With the PEC dofs zeroed the in-sheet ring
  // reads zero (the tangential field vanishes on both sides of a PEC sheet): the stretch
  // contributes ~0.

  // (a) The PEC-sheet face: the convex island, convex-corner-90 (4 placements).
  const auto convex =
      ReadMaxwellRing(ConvexMaxwellIslandConfig(), false, temp.temp_dir, "convex",
                      "convex-corner-90", test::CrackPolicy::Legacy);
  // (b) The uncracked interior face: the aperture, concave-corner-90 (4 placements).
  auto concave_config = ConvexMaxwellIslandConfig();
  concave_config["Solver"]["SurfaceResponseCorrection"]["Library"] =
      concave_library_3d_path.string();
  const auto concave = ReadMaxwellRing(concave_config, true, temp.temp_dir, "aperture",
                                       "concave-corner-90", test::CrackPolicy::Legacy);
  if (Mpi::Root(Mpi::World()))
  {
    for (const auto &reading : {convex, concave})
    {
      INFO(reading.name);
      CHECK(reading.patches == 4);
      CHECK(reading.max_difference_error < 1.0e-10);
    }
    // Information (decision 831 MINOR-2): the in-sheet ring on the PEC dofs zeroed field.
    CHECK(convex.max_zeroed_in_plane_trace < 1.0e-10);
    WARN("in-sheet ring trace with the PEC dofs zeroed: convex "
         << convex.max_zeroed_in_plane_trace << ", aperture "
         << concave.max_zeroed_in_plane_trace);
  }
#endif
}

// P3-23 (the sibling of the case above for upstream's default crack form, SL8-RESOLUTION s4
// (a)): the same two rings with Model.CrackInternalBoundaryElements unset, so the PEC sheet
// 9 is NOT cracked (the premise is required in ReadMaxwellRing) and the convex corner's z =
// 0 ring lies in the sheet's single (shared) faces instead of inside a crack. The
// expectations are the physics of the constant field E0 (exact in the lowest-order Nedelec
// space), not the cracked case's values: every knot-to-knot trace difference is -E0 . (q -
// p), every closed ring has a zero loop residual, each corner model is placed on the
// island's four corners, and with the PEC dofs zeroed the in-sheet ring reads zero.
TEST_CASE_METHOD(test::SurfaceResponseFiles,
                 "SurfaceResponseOperator Maxwell lines in uncracked PEC sheet faces",
                 "[surfaceresponseoperator][3d][maxwell][corner][Serial][Parallel]")
{
#if !defined(MFEM_USE_GSLIB)
  SKIP("SurfaceResponseOperator requires MFEM_USE_GSLIB");
#else
  const auto convex = ReadMaxwellRing(ConvexMaxwellIslandConfig(test::CrackPolicy::Default),
                                      false, temp.temp_dir, "convex-uncracked",
                                      "convex-corner-90", test::CrackPolicy::Default);
  auto concave_config = ConvexMaxwellIslandConfig(test::CrackPolicy::Default);
  concave_config["Solver"]["SurfaceResponseCorrection"]["Library"] =
      concave_library_3d_path.string();
  const auto concave =
      ReadMaxwellRing(concave_config, true, temp.temp_dir, "aperture-uncracked",
                      "concave-corner-90", test::CrackPolicy::Default);
  if (Mpi::Root(Mpi::World()))
  {
    for (const auto &reading : {convex, concave})
    {
      INFO(reading.name);
      CHECK(reading.patches == 4);
      CHECK(reading.max_difference_error < 1.0e-10);
    }
    CHECK(convex.max_zeroed_in_plane_trace < 1.0e-10);
    WARN("uncracked: in-sheet ring trace with the PEC dofs zeroed: convex "
         << convex.max_zeroed_in_plane_trace << ", aperture "
         << concave.max_zeroed_in_plane_trace << "; max trace-difference error: convex "
         << convex.max_difference_error << ", aperture " << concave.max_difference_error);
  }
#endif
}

}  // namespace palace
