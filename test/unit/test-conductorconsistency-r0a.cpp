// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "fixtures.hpp"
#include "surfaceresponse-fixtures.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <map>
#include <memory>
#include <set>
#include <string>
#include <vector>
#include <fmt/format.h>
#include <mfem.hpp>
#include <catch2/catch_test_macros.hpp>
#include <nlohmann/json.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>
#include "fem/gridfunction.hpp"
#include "fem/mesh.hpp"
#include "linalg/ksp.hpp"
#include "linalg/vector.hpp"
#include "models/laplaceoperator.hpp"
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

// The SHORTED CPW STUB of the ITEM-0 test (idfix DESIGN B.2) in the unit box: a 2 x 1 x 2
// box (hexahedra, h = 0.125) whose plane y = 0.5 is PEC metal (attribute 9) EXCEPT two
// slots of bare dielectric, x in [0.5, 2] x z in [0.75, 0.875] and x in [0.5, 2] x z in
// [1.125, 1.25], running into the wall x = 2 (the open end); the strip z in [0.875, 1.125]
// between them is the ground plane itself continuing at x = 0.5 (the SHORT): ONE metal
// component, so every edge of the 4-edge cross-section (offsets 0 / 0.125 / 0.375 / 0.5,
// gaps 0.125, strip 0.25, all below 2R = 0.4) carries conductor 1 - the one-label
// ParallelEdgeCluster of the CHR-012 key a670c49a86d8 (gaps and strip 3 um = 1.58 R).
// The metal plane of the stub: attribute 9 on the interior faces of y = 0.5 outside the two
// slots; the mesh finalized.
void AddShortedCpwMetalPlane(mfem::Mesh &serial)
{
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
    double xmin = 2.0, xmax = 0.0, zmin = 2.0, zmax = 0.0;
    for (const int vertex : vertices)
    {
      const double *point = serial.GetVertex(vertex);
      on_plane = on_plane && std::abs(point[1] - 0.5) < 1.0e-12;
      xmin = std::min(xmin, point[0]);
      xmax = std::max(xmax, point[0]);
      zmin = std::min(zmin, point[2]);
      zmax = std::max(zmax, point[2]);
    }
    if (!on_plane)
    {
      continue;
    }
    const bool in_slot_x = xmin >= 0.5 - 1.0e-12;
    const bool in_slot_a = zmin >= 0.75 - 1.0e-12 && zmax <= 0.875 + 1.0e-12;
    const bool in_slot_b = zmin >= 1.125 - 1.0e-12 && zmax <= 1.25 + 1.0e-12;
    if (in_slot_x && (in_slot_a || in_slot_b))
    {
      continue;  // bare dielectric
    }
    serial.AddBdrElement(serial.GetFace(face)->Duplicate(&serial));
    serial.SetBdrAttribute(serial.GetNBE() - 1, 9);
  }
  serial.FinalizeTopology();
  serial.Finalize();
}

mfem::Mesh MakeSerialShortedCpwMesh(mfem::Element::Type type = mfem::Element::HEXAHEDRON)
{
  mfem::Mesh serial = mfem::Mesh::MakeCartesian3D(16, 4, 16, type, 2.0, 1.0, 2.0);
  AddShortedCpwMetalPlane(serial);
  return serial;
}

// The stub with a VOID (U-H14a), TETRAHEDRAL (the exact segment intersections the line
// coverage reads exist on simplices only; a hexahedral mesh samples the lines blindly and
// fails at the point location, the H9a site): the twelve tetrahedra of the two top-layer
// cells y in [0.75, 1] over the STRIP (z in [0.875, 1.125]) at the x-cell `x_cell` (x in
// [x_cell h, (x_cell + 1) h], h = 0.125) are removed, so the ring's top line (y = 0.8, z
// from 0.8125 to 1.1875) of a patch at that x crosses a region no element covers while its
// KNOTS (over the slots, z = 0.8125 / 1.1875) stay inside elements: the DomainBoundary pass
// (which locates the basis points) keeps the patch and the line coverage is what fails.
// The void's faces get the boundary attribute 7; the box walls keep the Cartesian
// attributes; the metal plane is untouched (the void does not reach y = 0.5).
mfem::Mesh MakeSerialShortedCpwMeshWithVoid(int x_cell)
{
  constexpr double h = 0.125;
  mfem::Mesh base =
      mfem::Mesh::MakeCartesian3D(16, 4, 16, mfem::Element::TETRAHEDRON, 2.0, 1.0, 2.0);
  auto InVoid = [&](const mfem::Vector &centre)
  {
    return centre[1] > 0.75 && centre[2] > 0.875 && centre[2] < 1.125 &&
           centre[0] > x_cell * h && centre[0] < (x_cell + 1) * h;
  };
  mfem::Mesh serial(3, base.GetNV(), base.GetNE() - 12, 0, 3);
  for (int v = 0; v < base.GetNV(); v++)
  {
    serial.AddVertex(base.GetVertex(v));
  }
  int removed = 0;
  for (int e = 0; e < base.GetNE(); e++)
  {
    mfem::Vector centre(3);
    base.GetElementCenter(e, centre);
    if (InVoid(centre))
    {
      removed++;
      continue;
    }
    serial.AddElement(base.GetElement(e)->Duplicate(&serial));
  }
  MFEM_VERIFY(removed == 12, "The stub void must remove exactly twelve tetrahedra!");
  serial.FinalizeTopology(true);  // the boundary elements of every exterior face
  // The wall attributes from the base mesh by the face's vertex set; a face no wall had
  // (the void's) gets 7.
  auto Key = [](const mfem::Array<int> &vertices)
  {
    std::vector<int> key(vertices.GetData(), vertices.GetData() + vertices.Size());
    std::sort(key.begin(), key.end());
    return key;
  };
  std::map<std::vector<int>, int> wall_attributes;
  for (int be = 0; be < base.GetNBE(); be++)
  {
    mfem::Array<int> vertices;
    base.GetBdrElementVertices(be, vertices);
    wall_attributes[Key(vertices)] = base.GetBdrAttribute(be);
  }
  for (int be = 0; be < serial.GetNBE(); be++)
  {
    mfem::Array<int> vertices;
    serial.GetBdrElementVertices(be, vertices);
    const auto wall = wall_attributes.find(Key(vertices));
    serial.SetBdrAttribute(be, wall != wall_attributes.end() ? wall->second : 7);
  }
  AddShortedCpwMetalPlane(serial);
  return serial;
}

// The synthetic quasi-TEM field of the stub (idfix DESIGN B.2 (c)): the strip at the
// voltage V(x) against both grounds, E_z = -dphi/dz across the two gaps (phi rises from 0
// on the left ground to V on the strip across slot A and falls back across slot B), 0
// elsewhere; V(x) = `amplitude` x the indicator of [x0, x1] (+ `tiny` on [x2, x3]). Real
// part only.
struct StubField
{
  double amplitude = 1.0, x0 = 1.5, x1 = 2.0;
  double tiny = 0.0, x2 = 0.5, x3 = 1.0;
  double VoltageAt(double x) const
  {
    if (x >= x0 - 1.0e-12 && x <= x1 + 1.0e-12)
    {
      return amplitude;
    }
    if (x >= x2 - 1.0e-12 && x <= x3 + 1.0e-12)
    {
      return tiny;
    }
    return 0.0;
  }
  void operator()(const mfem::Vector &x, mfem::Vector &E) const
  {
    E.SetSize(3);
    E = 0.0;
    const double voltage = VoltageAt(x[0]);
    if (x[2] > 0.75 && x[2] < 0.875)
    {
      E[2] = -voltage / 0.125;  // phi rises toward the strip: E points away from it
    }
    else if (x[2] > 1.125 && x[2] < 1.25)
    {
      E[2] = voltage / 0.125;
    }
  }
};

// The library of the stub tests: the convex library's header (R = 0.2, SA) with ONE
// one-label 4-edge ParallelEdgeCluster model (offsets 0 / 0.125 / 0.375 / 0.5, GapDirection
// + - + -, every Conductor 1, one Reference on the strip): a 4-knot ring contour around the
// cross-section (the top knots at v = `ring_top_v` above the plane, the bottom ones at
// v = -0.1 below it), diagonal response matrices. Written by the root; returns the path.
fs::path WriteStubLibrary(const fs::path &directory, const fs::path &template_library,
                          double ring_top_v, const std::string &suffix)
{
  const auto library_path =
      directory / fmt::format("fabrication-process-r0a-stub-3d-{}.json", suffix);
  const auto points_path = directory / fmt::format("r0a-stub-basis-points-{}.csv", suffix);
  const auto domain_path = directory / fmt::format("r0a-stub-domain-{}.csv", suffix);
  const auto surface_path = directory / fmt::format("r0a-stub-surface-{}.csv", suffix);
  if (Mpi::Root(Mpi::World()))
  {
    {
      std::ofstream output(points_path);
      output << "x,y,z\n"
             << "0.0625," << ring_top_v << ",0.0\n"
             << "0.4375," << ring_top_v << ",0.0\n"
             << "0.4375,-0.1,0.0\n"
             << "0.0625,-0.1,0.0\n";
    }
    {
      std::ofstream output(domain_path);
      output << "basis_i,basis_j,Q_ij (J)\n";
      for (int i = 1; i <= 4; i++)
      {
        for (int j = i; j <= 4; j++)
        {
          output << i << "," << j << "," << (i == j ? 1.0e-12 : 0.0) << "\n";
        }
      }
    }
    {
      std::ofstream output(surface_path);
      output << "interface,edge,R (m),basis_i,basis_j,Q_ij (J),Q_total_ij (J)\n";
      for (int edge = 1; edge <= 4; edge++)
      {
        for (int i = 1; i <= 4; i++)
        {
          for (int j = i; j <= 4; j++)
          {
            output << "1," << edge << ",0.2," << i << "," << j << ","
                   << (i == j ? 0.25e-12 : 0.0) << "," << (i == j ? 0.25e-12 : 0.0) << "\n";
          }
        }
      }
    }
    std::ifstream input(template_library);
    REQUIRE(input);
    json library = json::parse(input);
    library["Name"] = "unit-test-process-r0a-stub-3d";
    json model = library["Models"][0];
    model["Name"] = "cpw-one-label-four-edge";
    model["Topology"] = "ParallelEdgeCluster";
    model["Interfaces"] = {{{"Type", "SA"}, {"Coupon", 1}}};
    model["BasisPoints"] = points_path.string();
    model["FabricatedMatrix"] = domain_path.string();
    model["ThinMatrix"] = domain_path.string();
    model["FabricatedSurfaceMatrix"] = surface_path.string();
    model["ThinSurfaceMatrix"] = surface_path.string();
    model["Reference"] = {0.25, 0.0, 0.0};
    model["EdgeOffsetTolerance"] = 1.0e-6;
    model["Edges"] = {{{"Offset", 0.0}, {"GapDirection", 1}, {"Conductor", 1}},
                      {{"Offset", 0.125}, {"GapDirection", -1}, {"Conductor", 1}},
                      {{"Offset", 0.375}, {"GapDirection", 1}, {"Conductor", 1}},
                      {{"Offset", 0.5}, {"GapDirection", -1}, {"Conductor", 1}}};
    for (const char *key : {"Separation", "SeparationTolerance", "Angle", "AngleTolerance",
                            "TraceMesh", "TraceMeshTriangles", "ConductorReferences"})
    {
      model.erase(key);
    }
    library["Models"] = {model};
    std::ofstream output(library_path);
    output << library.dump(2) << "\n";
  }
  Mpi::Barrier(Mpi::World());
  return library_path;
}

}  // namespace

// U-R0a (the identification fix round, idfix DESIGN A.1 (c) / A.3 (b) / B.3 (g); decisions
// 1062 / 1064): the conductor-consistency gate on every model class. Part 1 = the LOCAL
// METAL PIECES of a 2D cross-section (one rule, shared by name with the library
// generators' LocalPiece, U-PIECES) and the ELECTROSTATIC path's translational probes:
// every translational patch carries the device metal point of every model Edge at both
// end faces of its cell, read against the reference its metal is held at; on the device's
// own solve they read the Dirichlet value (0 excluded, the energies unchanged), on a
// potential that varies along the edge they fire, and RecordOnly records without excluding.
TEST_CASE("LocalPiecesOfParallelEdges groups a cross-section's edges into metal pieces",
          "[surfaceresponseoperator][conductorconsistency][pieces][Serial]")
{
  // A CPW (ground-left / strip / ground-right): the record's 2 / 2 / 2 um stack and the
  // CHR-012 3 / 3 / 3 um key a670c49a86d8 (idfix DESIGN 0.4) have the same pattern.
  CHECK(LocalPiecesOfParallelEdges({0.0, 2.0, 4.0, 6.0}, {1, -1, 1, -1}) ==
        std::vector<int>{1, 2, 2, 3});
  CHECK(LocalPiecesOfParallelEdges({0.0, 3.0, 6.0, 9.0}, {1, -1, 1, -1}) ==
        std::vector<int>{1, 2, 2, 3});
  // A three-edge stack (ground / strip, the gap on one side only): two pieces.
  CHECK(LocalPiecesOfParallelEdges({0.0, 2.0, 4.0}, {1, -1, 1}) ==
        std::vector<int>{1, 2, 2});
  // A gap (metal outside): two pieces; a strip (metal between): one; an isolated edge: one.
  CHECK(LocalPiecesOfParallelEdges({-1.0, 1.0}, {1, -1}) == std::vector<int>{1, 2});
  CHECK(LocalPiecesOfParallelEdges({-1.0, 1.0}, {-1, 1}) == std::vector<int>{1, 1});
  CHECK(LocalPiecesOfParallelEdges({0.0}, {1}) == std::vector<int>{1});
  CHECK(LocalPiecesOfParallelEdges({0.0}, {-1}) == std::vector<int>{1});
  // The 9-edge v3-j stack fb9f050d1e24 (A.3 (c): 5 pieces vs 1 reference).
  CHECK(LocalPiecesOfParallelEdges({0, 1, 2, 3, 4, 5, 6, 7, 8},
                                   {1, -1, 1, -1, 1, -1, 1, -1, 1}) ==
        std::vector<int>{1, 2, 2, 3, 3, 4, 4, 5, 5});
  // The piece of a lateral offset: the strip's centre, the grounds, a gap (0).
  const std::vector<double> offsets = {0.0, 2.0, 4.0, 6.0};
  const std::vector<int> gaps = {1, -1, 1, -1};
  const auto pieces = LocalPiecesOfParallelEdges(offsets, gaps);
  CHECK(LocalPieceOfOffset(offsets, gaps, pieces, 3.0) == 2);
  CHECK(LocalPieceOfOffset(offsets, gaps, pieces, -1.9) == 1);
  CHECK(LocalPieceOfOffset(offsets, gaps, pieces, 0.0) == 1);
  CHECK(LocalPieceOfOffset(offsets, gaps, pieces, 2.0) == 2);
  CHECK(LocalPieceOfOffset(offsets, gaps, pieces, 7.5) == 3);
  CHECK(LocalPieceOfOffset(offsets, gaps, pieces, 1.0) == 0);
  CHECK(LocalPieceOfOffset(offsets, gaps, pieces, 5.0) == 0);
  // The one-reference gap model of record (same-conductor-gap-2um: Reference u = -1.95):
  // the left piece.
  CHECK(LocalPieceOfOffset({-1.0, 1.0}, {1, -1}, {1, 2}, -1.95) == 1);
  CHECK(LocalPieceOfOffset({-1.0, 1.0}, {1, -1}, {1, 2}, 1.95) == 2);
  // Not a cross-section: two consecutive gaps or two consecutive metals.
  CHECK_THROWS_WITH(LocalPiecesOfParallelEdges({0.0, 1.0}, {1, 1}),
                    ContainsSubstring("no metal between them"));
  CHECK_THROWS_WITH(LocalPiecesOfParallelEdges({0.0, 1.0}, {-1, -1}),
                    ContainsSubstring("no gap between them"));
  CHECK_THROWS_WITH(LocalPiecesOfParallelEdges({1.0, 0.0}, {1, -1}),
                    ContainsSubstring("strictly increasing"));
  CHECK_THROWS_WITH(LocalPiecesOfParallelEdges({0.0}, {2}), ContainsSubstring("-1 or 1"));
}

// The ELECTROSTATIC path on the PEC island (x, z in [0.25, 0.75] on y = 0.5, the terminal
// at 1 V; MakeIslandMesh) with the convex library (isolated edge + convex-corner-90, SA
// only) under the surface mortar: the 16 isolated-edge cells carry two probes each (their
// one edge at both end faces); the fixture's corner coupon has no conductor vertex on its
// trace mesh and stays recorded untestable (UnprobedModels), as today.
TEST_CASE_METHOD(
    test::SurfaceResponseFiles,
    "SurfaceResponseOperator conductor-consistency gate on translational "
    "patches (R0-a, electrostatic)",
    "[surfaceresponseoperator][conductorconsistency][r0a][3d][Serial][Parallel]")
{
#if !defined(MFEM_USE_GSLIB)
  SKIP("SurfaceResponseOperator requires MFEM_USE_GSLIB");
#else
  const auto library_path = temp.temp_dir / "fabrication-process-r0a-convex-3d.json";
  if (Mpi::Root(Mpi::World()))
  {
    std::ifstream input(convex_library_3d_path);
    REQUIRE(input);
    json library = json::parse(input);
    library["Name"] = "unit-test-process-r0a-convex-3d";
    for (auto &model : library["Models"])
    {
      model["Interfaces"] = {{{"Type", "SA"}, {"Coupon", 1}}};
    }
    std::ofstream output(library_path);
    output << library.dump(2) << "\n";
  }
  Mpi::Barrier(Mpi::World());
  json config = IslandConfig();
  auto &correction = config["Solver"]["Electrostatic"]["ResponseCorrection"];
  correction.erase("PatchConstruction");
  correction["Library"] = library_path.string();
  correction["TraceCoupling"] = "SurfaceMortar";
  correction["MortarOversampling"] = 2;
  auto Build = [&](const json &run_config, std::vector<std::unique_ptr<Mesh>> &meshes)
  {
    IoData iodata(run_config, false);
    iodata.boundaries.cracked_attributes.insert(9);
    meshes.push_back(std::make_unique<Mesh>(MakeIslandMesh()));
    return iodata;
  };

  SECTION("the device's own solve: every translational probe reads the Dirichlet value")
  {
    std::vector<std::unique_ptr<Mesh>> meshes;
    IoData iodata = Build(config, meshes);
    CHECK_FALSE(
        iodata.solver.electrostatic.response_correction->conductor_consistency_record_only);
    LaplaceOperator laplace(iodata, meshes);
    auto K = laplace.GetStiffnessMatrix();
    KspSolver ksp(iodata, laplace.GetH1Spaces());
    ksp.SetOperators(*K, *K);
    Vector V, RHS;
    laplace.GetExcitationVector(1, *K, V, RHS);
    ksp.Mult(RHS, V);
    REQUIRE(ksp.GetConverged());
    SurfaceResponseOperator response(iodata, laplace);
    REQUIRE(response.GetPatchCount() == 20);  // 16 isolated cells + 4 corners
    REQUIRE(response.HasConductorConsistencyProbes());
    const auto before = response.GetElectrostaticResponse(V);
    const auto records = response.ApplyConductorConsistencyGate(V, 1);
    REQUIRE(records.size() == 16);
    int translational = 0, spatial = 0;
    for (const auto &record : records)
    {
      CHECK(record.source == 1);
      CHECK_FALSE(record.excluded);
      CHECK_FALSE(record.record_only);
      CHECK_FALSE(record.maxwell);
      CHECK_THAT(record.excitation_potential, WithinRel(1.0, 1.0e-12));
      if (record.probe_kind == 1)
      {
        translational++;
        // One edge, two end faces: the isolated edge's single piece, conductor 1.
        CHECK(record.plane_knots == 2);
        CHECK(record.off_plane_knots == 0);
        CHECK(record.adjacent_knots == 0);
        CHECK(record.worst_edge == 0);
        CHECK(record.worst_piece == 1);
        CHECK((record.worst_end_face == 0 || record.worst_end_face == 1));
        CHECK(record.worst_conductor == 1);
        // The probe lies on the island perimeter (y = 0.5; x or z on the island edge) and
        // reads the 1-V Dirichlet sheet like the reference: deviation 0 to roundoff.
        CHECK_THAT(record.worst_point[1], WithinAbs(0.5, 1.0e-9));
        const bool on_x_edge = std::abs(record.worst_point[0] - 0.25) < 1.0e-9 ||
                               std::abs(record.worst_point[0] - 0.75) < 1.0e-9;
        const bool on_z_edge = std::abs(record.worst_point[2] - 0.25) < 1.0e-9 ||
                               std::abs(record.worst_point[2] - 0.75) < 1.0e-9;
        CHECK((on_x_edge || on_z_edge));
        CHECK(record.max_deviation < 1.0e-9);
        CHECK(record.max_ratio < 1.0e-9);
        CHECK(record.cell_length > 0.0);
      }
      else
      {
        spatial++;
        CHECK(record.worst_edge == -1);
        CHECK(record.worst_piece == 0);
        CHECK(record.worst_end_face == -1);
      }
    }
    CHECK(translational == 16);
    CHECK(spatial == 0);
    // The probes change no trace: the energies after the gate are the energies before.
    const auto after = response.GetElectrostaticResponse(V);
    CHECK(after.domain_correction == before.domain_correction);
    for (const auto &[interface, energy] : before.fabricated_surface_energy)
    {
      CHECK(after.fabricated_surface_energy.at(interface) == energy);
    }
    const auto statistics = response.GetStatistics();
    const auto &diagnostics = statistics["Diagnostics"]["ConductorConsistency"];
    CHECK(diagnostics["Count"].get<int>() == 0);
    CHECK(diagnostics["TestedPatches"].get<int>() == 16);
    CHECK(diagnostics["Records"].size() == 16);
    CHECK_FALSE(diagnostics.contains("RecordOnly"));
    // The corner coupon of this fixture: no conductor vertex on its trace mesh (recorded
    // untestable, never a silent zero); the isolated edge is probed, so it is not listed.
    REQUIRE(diagnostics["UnprobedModels"].size() == 1);
    CHECK(diagnostics["UnprobedModels"][0]["Model"] == "convex-corner-90");
    for (const char *leaf :
         {"ProbeRule", "MaxwellFloorRule", "DeviationNorm", "ProbePathRule"})
    {
      CHECK(diagnostics.contains(leaf));
    }
    int translational_rows = 0;
    for (const auto &entry : diagnostics["Records"])
    {
      if (entry["ProbeKind"] == "TranslationalEdge")
      {
        translational_rows++;
        CHECK(entry["WorstPiece"].get<int>() == 1);
        CHECK_FALSE(entry["Maxwell"].get<bool>());
        CHECK_FALSE(entry["RecordOnly"].get<bool>());
      }
    }
    CHECK(translational_rows == 16);
  }

  SECTION(
      "a potential varying along the edge fires the probes; RecordOnly keeps the weight")
  {
    for (const bool record_only : {false, true})
    {
      json run_config = config;
      if (record_only)
      {
        run_config["Solver"]["Electrostatic"]["ResponseCorrection"]
                  ["ConductorConsistency"] = {{"RecordOnly", true}};
      }
      std::vector<std::unique_ptr<Mesh>> meshes;
      IoData iodata = Build(run_config, meshes);
      CHECK(iodata.solver.electrostatic.response_correction
                ->conductor_consistency_record_only == record_only);
      LaplaceOperator laplace(iodata, meshes);
      SurfaceResponseOperator response(iodata, laplace);
      REQUIRE(response.GetPatchCount() == 20);  // 16 isolated cells + 4 corners
      // A potential growing along z, 0 on the island's z = 0.5 midline, +-0.25 at its
      // z = 0.25 / 0.75 edges: on the two edges x = const the cell-end probes read the
      // edge coordinate, which differs from the cell-centre reference by half the cell
      // (0.025) against a trace amplitude of order the mortar's span (R = 0.2) times the
      // slope: ratios of order 0.1 >> 0.02. Not a device solve: a coupon whose metal is
      // not at the potential its basis assumes, the gate's defining case.
      mfem::ParGridFunction potential(&laplace.GetH1Space().Get());
      mfem::FunctionCoefficient coefficient([](const mfem::Vector &x)
                                            { return x[2] - 0.5; });
      potential.ProjectCoefficient(coefficient);
      Vector V;
      potential.GetTrueDofs(V);
      const double weight_before = response.GetPatchWeight();
      const auto records = response.ApplyConductorConsistencyGate(V, 1);
      REQUIRE(records.size() == 16);
      int fired = 0;
      for (const auto &record : records)
      {
        if (record.probe_kind != 1)
        {
          continue;
        }
        CHECK(record.record_only == record_only);
        if (record.excluded)
        {
          fired++;
          CHECK(record.max_ratio > SurfaceResponseOperator::kConductorConsistencyTolerance);
          CHECK(record.worst_end_face >= 0);
        }
      }
      CHECK(fired > 0);
      const auto statistics = response.GetStatistics();
      const auto &diagnostics = statistics["Diagnostics"]["ConductorConsistency"];
      if (record_only)
      {
        // The verdicts stand in the record; nothing is excluded, nothing served raw.
        CHECK(diagnostics["RecordOnly"].get<bool>());
        CHECK(diagnostics["Count"].get<int>() == 0);
        CHECK(diagnostics["ExcludedPatches"].empty());
        CHECK(diagnostics["RecordOnlyVerdicts"].size() == static_cast<std::size_t>(fired));
        CHECK_FALSE(diagnostics.contains("RawPortions"));
        CHECK(response.GetPatchWeight() == weight_before);
        CHECK(response.GetConductorConsistencyPortions().empty());
      }
      else
      {
        CHECK(diagnostics["Count"].get<int>() == fired);
        CHECK(diagnostics["ExcludedPatches"].size() == static_cast<std::size_t>(fired));
        CHECK_FALSE(diagnostics.contains("RecordOnlyVerdicts"));
        CHECK(diagnostics.contains("RawPortions"));
        CHECK(response.GetPatchWeight() < weight_before);
        CHECK_FALSE(response.GetConductorConsistencyPortions().empty());
      }
    }
  }
#endif
}

// The MAXWELL path (idfix DESIGN B.3 (g) / SUPERVISOR AMENDMENT 3 C.1; decision 1064): the
// shorted CPW stub with the one-label 4-edge ParallelEdgeCluster model of its
// cross-section, collocated (the Maxwell form), under synthetic complex fields. (a) the
// strip at V, real: the strip probes read V against the ground anchor, A_p ~ V, ratio ~1 >>
// 0.02 -> EXCLUDED, FloorApplied false; (b) the strip at V e^{j pi/2}: the same D and
// verdict (the complex norm); (c) the strip at 0 (a field normal to the plane only): kept;
// (d) a ~0-field copy (A_p < F): FloorApplied true, ratio ~0, kept; (f) two sources gated
// in both orders give identical verdicts and FloorApplied flags (the floor's set is the
// PRE-GATE applied set).
TEST_CASE_METHOD(test::SurfaceResponseFiles,
                 "SurfaceResponseOperator conductor-consistency gate on the Maxwell path "
                 "(R0-a, the shorted CPW stub)",
                 "[surfaceresponseoperator][conductorconsistency][r0a][maxwell][3d][Serial]"
                 "[Parallel]")
{
#if !defined(MFEM_USE_GSLIB)
  SKIP("SurfaceResponseOperator requires MFEM_USE_GSLIB");
#else
  // The library (WriteStubLibrary): the one-label 4-edge cluster model, the ring's top
  // knots at v = 0.1.
  const auto library_path =
      WriteStubLibrary(temp.temp_dir, convex_library_3d_path, 0.1, "r0a");

  json config = ConvexMaxwellIslandConfig();
  config["Boundaries"]["Ground"]["Attributes"] = {1, 2, 3, 4, 5, 6, 9};
  auto &correction = config["Solver"]["SurfaceResponseCorrection"];
  correction.erase("PatchConstruction");
  correction["Library"] = library_path.string();
  correction["UnmatchedPolicy"] = "Warn";  // the short's corner cluster has no model: raw
  correction["TraceCoupling"] = "Collocated";
  auto Build = [&](const json &run_config, std::vector<std::unique_ptr<Mesh>> &meshes)
  {
    IoData iodata(run_config, false);
    iodata.boundaries.cracked_attributes.insert(9);
    mfem::Mesh serial = MakeSerialShortedCpwMesh();
    meshes.push_back(
        std::make_unique<Mesh>(std::make_unique<mfem::ParMesh>(Mpi::World(), serial)));
    return iodata;
  };
  auto Field = [](SpaceOperator &space, const StubField &stub, double phase)
  {
    GridFunction field(space.GetNDSpace(), true);
    mfem::VectorFunctionCoefficient coefficient(
        3, [&](const mfem::Vector &x, mfem::Vector &E) { stub(x, E); });
    field.Real().ProjectCoefficient(coefficient);
    field.Imag() = 0.0;
    ComplexVector E(field.Real().ParFESpace()->GetTrueVSize());
    field.Real().GetTrueDofs(E.Real());
    field.Imag().GetTrueDofs(E.Imag());
    if (phase != 0.0)
    {
      // E e^{j phase}: the imaginary part from the real one.
      E.Imag() = E.Real();
      E.Imag() *= std::sin(phase);
      E.Real() *= std::cos(phase);
    }
    return E;
  };
  // The CPW cells: the records with translational probes (the short's corner cluster is
  // unmatched, raw: no patch).
  auto StubRecords =
      [](const std::vector<SurfaceResponseOperator::ConductorConsistencyRecord> &records)
  {
    std::vector<SurfaceResponseOperator::ConductorConsistencyRecord> stub;
    for (const auto &record : records)
    {
      if (record.probe_kind == 1)
      {
        stub.push_back(record);
      }
    }
    return stub;
  };

  std::vector<std::unique_ptr<Mesh>> meshes;
  IoData iodata = Build(config, meshes);
  SpaceOperator space(iodata, meshes);
  SurfaceResponseOperator response(iodata, space);
  REQUIRE(response.GetPatchCount() > 0);
  REQUIRE(response.HasConductorConsistencyProbes());
  const auto tolerance = SurfaceResponseOperator::kConductorConsistencyTolerance;
  const auto floor = SurfaceResponseOperator::kConductorConsistencyAmplitudeFloor;
  const int patch_count = response.GetPatchCount();

  SECTION("(c) the strip at the ground potential (a field normal to the plane): kept")
  {
    GridFunction field(space.GetNDSpace(), true);
    mfem::VectorFunctionCoefficient coefficient(3,
                                                [](const mfem::Vector &, mfem::Vector &E)
                                                {
                                                  E.SetSize(3);
                                                  E = 0.0;
                                                  E[1] = 1.0;
                                                });
    field.Real().ProjectCoefficient(coefficient);
    field.Imag() = 0.0;
    ComplexVector E(field.Real().ParFESpace()->GetTrueVSize());
    field.Real().GetTrueDofs(E.Real());
    field.Imag().GetTrueDofs(E.Imag());
    const auto records =
        StubRecords(response.ApplyConductorConsistencyGate(E, 0, 5.0, false));
    REQUIRE(records.size() == static_cast<std::size_t>(patch_count));
    for (const auto &record : records)
    {
      CHECK(record.maxwell);
      CHECK_THAT(record.frequency, WithinAbs(5.0, 0.0));
      CHECK_FALSE(record.self_consistent_pre_gate);
      CHECK(record.plane_knots == 8);  // 4 edges x 2 end faces
      CHECK(record.amplitude > 0.0);   // the ring's vertical legs read E_y x 0.2
      CHECK(record.normalization == record.amplitude);
      // The w = 0 legs across the gaps run along z: E_y gives them nothing.
      CHECK(record.max_ratio < 1.0e-8);
      CHECK_FALSE(record.excluded);
    }
    CHECK(response.GetPatchWeight() > 0.0);
  }

  SECTION("(a) / (b) the strip at V, real and at e^{j pi/2}: excluded with the same D")
  {
    StubField stub;
    stub.amplitude = 1.0;
    stub.x0 = 0.5;
    stub.x1 = 2.0;
    std::vector<std::unique_ptr<Mesh>> phase_meshes;
    IoData phase_iodata = Build(config, phase_meshes);
    SpaceOperator phase_space(phase_iodata, phase_meshes);
    SurfaceResponseOperator phase_response(phase_iodata, phase_space);
    const auto real_records = StubRecords(
        response.ApplyConductorConsistencyGate(Field(space, stub, 0.0), 0, 5.0, false));
    const auto phase_records = StubRecords(phase_response.ApplyConductorConsistencyGate(
        Field(phase_space, stub, 0.5 * M_PI), 0, 5.0, false));
    REQUIRE(real_records.size() == static_cast<std::size_t>(patch_count));
    REQUIRE(phase_records.size() == real_records.size());
    int excluded = 0;
    for (std::size_t i = 0; i < real_records.size(); i++)
    {
      const auto &record = real_records[i];
      const auto &phased = phase_records[i];
      CHECK(record.patch == phased.patch);
      // The strip probes (edges 1 / 2, piece 2) read the gap voltage V against the ground
      // anchor (edge 0, piece 1): D ~ V, the trace amplitude ~ V (a ring knot above a gap
      // reads about half of it): ratio ~1 >> 0.02.
      CHECK(record.worst_piece == 2);
      CHECK((record.worst_edge == 1 || record.worst_edge == 2));
      CHECK_THAT(record.max_deviation, WithinRel(1.0, 0.15));
      CHECK(record.max_ratio > 10 * tolerance);
      CHECK(record.normalization == record.amplitude);  // FloorApplied false
      CHECK(record.excluded);
      // (b): |e^{j pi/2} (P - C)| = |P - C|.
      CHECK_THAT(phased.max_deviation, WithinRel(record.max_deviation, 1.0e-10));
      CHECK_THAT(phased.amplitude, WithinRel(record.amplitude, 1.0e-10));
      CHECK(phased.excluded == record.excluded);
      excluded++;
    }
    CHECK(excluded == patch_count);
    // RecordRaw: every stub patch excluded (weight 0, sticky), its raw claim served RAW.
    CHECK(response.GetPatchWeight() == 0.0);
    CHECK_FALSE(response.GetConductorConsistencyPortions().empty());
    const auto statistics = response.GetStatistics();
    const auto &diagnostics = statistics["Diagnostics"]["ConductorConsistency"];
    CHECK(diagnostics["Count"].get<int>() == patch_count);
    CHECK(diagnostics.contains("RawPortions"));
    for (const auto &entry : diagnostics["Records"])
    {
      CHECK(entry["Maxwell"].get<bool>());
      CHECK_FALSE(entry["FloorApplied"].get<bool>());
      CHECK(entry["ProbeKind"] == "TranslationalEdge");
    }
    // Sticky: a later source with the strip at 0 re-tests nothing (the weights stay 0).
    const auto later = StubRecords(response.ApplyConductorConsistencyGate(
        Field(space, StubField{0.0, 0.5, 2.0}, 0.0), 1, 6.0, false));
    CHECK(later.empty());
    CHECK(response.GetPatchWeight() == 0.0);
  }

  SECTION("(d) a ~0-field copy of the stub reads through the floor: kept, FloorApplied")
  {
    // The strip at V on x in [1.5, 2] (excluded) and at 5e-6 V on x in [0.5, 1] (A_p below
    // the floor 1e-3 x V: FloorApplied, ratio ~ 5e-3 < 0.02, kept); nothing on [1, 1.5]
    // (D = 0, kept).
    StubField stub;
    stub.amplitude = 1.0;
    stub.x0 = 1.5;
    stub.x1 = 2.0;
    stub.tiny = 5.0e-6;
    stub.x2 = 0.5;
    stub.x3 = 1.0;
    const auto records = StubRecords(
        response.ApplyConductorConsistencyGate(Field(space, stub, 0.0), 0, 5.0, false));
    REQUIRE(records.size() == static_cast<std::size_t>(patch_count));
    int excluded = 0, floored_kept = 0, zero_kept = 0;
    double largest = 0.0;
    for (const auto &record : records)
    {
      largest = std::max(largest, record.amplitude);
    }
    for (const auto &record : records)
    {
      CHECK_THAT(record.excitation_potential, WithinRel(largest, 1.0e-12));
      const bool floored = record.normalization > record.amplitude;
      if (record.excluded)
      {
        excluded++;
        CHECK_FALSE(floored);
        CHECK(record.worst_point[0] > 1.5 - 0.3);
      }
      else if (floored)
      {
        floored_kept++;
        CHECK_THAT(record.normalization, WithinRel(floor * largest, 1.0e-12));
        CHECK(record.max_ratio < tolerance);
        if (record.max_deviation > 1.0e-9)
        {
          CHECK(record.worst_point[0] < 1.0 + 0.3);
        }
        else
        {
          zero_kept++;
        }
      }
    }
    CHECK(excluded > 0);
    CHECK(floored_kept > 0);
    CHECK(zero_kept > 0);
    CHECK(excluded + floored_kept == patch_count);
  }

  SECTION("(f) two sources gated in either order give identical verdicts and floors")
  {
    // Source A: the strip at V on [1.5, 2] + 5e-6 V on [0.5, 1]; source B: V on [1.5, 2] +
    // 5e-6 V on [1, 1.5]. After A excludes the [1.5, 2] cells, B's floor still sees them
    // (C.1: the floor's set is the pre-gate applied set), so the [1, 1.5] cells pass
    // through the floor exactly as they do when B comes first.
    StubField a;
    a.amplitude = 1.0;
    a.x0 = 1.5;
    a.x1 = 2.0;
    a.tiny = 5.0e-6;
    a.x2 = 0.5;
    a.x3 = 1.0;
    StubField b = a;
    b.x2 = 1.0;
    b.x3 = 1.5;
    struct Verdict
    {
      bool excluded, floored;
    };
    auto Gate = [&](bool a_first)
    {
      std::vector<std::unique_ptr<Mesh>> order_meshes;
      IoData order_iodata = Build(config, order_meshes);
      SpaceOperator order_space(order_iodata, order_meshes);
      SurfaceResponseOperator order_response(order_iodata, order_space);
      std::map<std::pair<int, int>, Verdict> verdicts;  // (patch, source) -> verdict
      for (const int step : {0, 1})
      {
        const bool is_a = (step == 0) == a_first;
        const auto records = StubRecords(order_response.ApplyConductorConsistencyGate(
            Field(order_space, is_a ? a : b, 0.0), is_a ? 0 : 1, 5.0, false));
        for (const auto &record : records)
        {
          verdicts[{record.patch, record.source}] = {record.excluded, record.normalization >
                                                                          record.amplitude};
        }
      }
      return verdicts;
    };
    const auto a_first = Gate(true);
    const auto b_first = Gate(false);
    // Source A tests every patch in both orders; source B tests the ones A kept.
    int compared = 0;
    for (const auto &[key, verdict] : a_first)
    {
      const auto other = b_first.find(key);
      if (other == b_first.end())
      {
        continue;
      }
      CHECK(verdict.excluded == other->second.excluded);
      CHECK(verdict.floored == other->second.floored);
      compared++;
    }
    CHECK(compared > 0);
    int excluded_a = 0, excluded_b = 0;
    for (const auto &[key, verdict] : a_first)
    {
      excluded_a += verdict.excluded ? 1 : 0;
    }
    for (const auto &[key, verdict] : b_first)
    {
      excluded_b += verdict.excluded ? 1 : 0;
    }
    CHECK(excluded_a == excluded_b);
    CHECK(excluded_a > 0);
  }
#endif
}

// U-H14a (hardening-2 DESIGN row H14a; idfix DESIGN section 7 OI-MX-CONTOUR-COVERAGE-NC-
// REFINEMENT part (2); decisions 1072 Q2 / 1074): a Maxwell contour line that leaves the
// mesh beyond the crack tolerance no longer aborts under RecordRaw: its PATCH is excluded
// (weight 0, every line dropped from the quadrature) and recorded MaxwellPatchExcluded, its
// raw claim served RAW once as the DomainBoundary share; FailClosed keeps the abort (the
// ledger's deferred Enforce). The fixture = the shorted CPW stub with a VOID in the top
// layer over the strip and a ring whose top line (the knots at v = 0.3, y = 0.8, over the
// slots) runs through it: exactly the patches whose origin lies in the void's x-cell are
// excluded; the intact stub records nothing.
TEST_CASE_METHOD(test::SurfaceResponseFiles,
                 "SurfaceResponseOperator excludes a Maxwell patch whose contour line "
                 "leaves the mesh (H14a MaxwellPatchExcluded, the voided shorted CPW stub)",
                 "[surfaceresponseoperator][geometryfallback][maxwell][h14a][3d][Serial]"
                 "[Parallel]")
{
#if !defined(MFEM_USE_GSLIB)
  SKIP("SurfaceResponseOperator requires MFEM_USE_GSLIB");
#else
  constexpr double h = 0.125;
  const auto library_path =
      WriteStubLibrary(temp.temp_dir, convex_library_3d_path, 0.3, "h14a");
  json config = ConvexMaxwellIslandConfig();
  config["Boundaries"]["Ground"]["Attributes"] = {1, 2, 3, 4, 5, 6, 9};
  auto &correction = config["Solver"]["SurfaceResponseCorrection"];
  correction.erase("PatchConstruction");
  correction["Library"] = library_path.string();
  correction["UnmatchedPolicy"] = "Warn";
  correction["TraceCoupling"] = "Collocated";
  auto Build = [&](const json &run_config, mfem::Mesh serial,
                   std::vector<std::unique_ptr<Mesh>> &meshes)
  {
    IoData iodata(run_config, false);
    iodata.boundaries.cracked_attributes.insert(9);
    meshes.push_back(
        std::make_unique<Mesh>(std::make_unique<mfem::ParMesh>(Mpi::World(), serial)));
    return iodata;
  };
  auto StripAtGround = [](SpaceOperator &space)
  {
    GridFunction field(space.GetNDSpace(), true);
    mfem::VectorFunctionCoefficient coefficient(3,
                                                [](const mfem::Vector &, mfem::Vector &E)
                                                {
                                                  E.SetSize(3);
                                                  E = 0.0;
                                                  E[1] = 1.0;
                                                });
    field.Real().ProjectCoefficient(coefficient);
    field.Imag() = 0.0;
    ComplexVector E(field.Real().ParFESpace()->GetTrueVSize());
    field.Real().GetTrueDofs(E.Real());
    field.Imag().GetTrueDofs(E.Imag());
    return E;
  };
  auto InCell = [&](double x, int cell)
  { return x - cell * h > 1.0e-9 && (cell + 1) * h - x > 1.0e-9; };

  // The intact (tetrahedral) stub: the tall ring is covered everywhere, the ledger is
  // empty; the void's x-cell = the one holding a patch origin strictly inside it (the
  // root's patch table, broadcast), the expected exclusions = the patches whose origin
  // lies in that cell.
  int patch_count = 0;
  int x_cell = -1;
  int expected_excluded = 0;
  double intact_domain_boundary_length = 0.0;
  {
    std::vector<std::unique_ptr<Mesh>> meshes;
    IoData iodata =
        Build(config, MakeSerialShortedCpwMesh(mfem::Element::TETRAHEDRON), meshes);
    SpaceOperator space(iodata, meshes);
    SurfaceResponseOperator response(iodata, space);
    patch_count = response.GetPatchCount();
    REQUIRE(patch_count > 0);
    CHECK(response.GetGeometryFallbacks().empty());
    CHECK(response.GetMaxwellPatchExcludedPortionRecord().portions.empty());
    intact_domain_boundary_length = response.GetDomainBoundaryPortionRecord().length;
    if (Mpi::Root(Mpi::World()))
    {
      const auto &assignments = response.GetPatchAssignments();
      REQUIRE(static_cast<int>(assignments.size()) == patch_count);
      for (const auto &assignment : assignments)
      {
        const int cell = static_cast<int>(std::floor(assignment.origin[0] / h));
        if (cell >= 1 && cell <= 14 && InCell(assignment.origin[0], cell))
        {
          x_cell = cell;
          break;
        }
      }
      for (const auto &assignment : assignments)
      {
        expected_excluded += InCell(assignment.origin[0], x_cell) ? 1 : 0;
      }
    }
    Mpi::Broadcast(1, &x_cell, 0, Mpi::World());
    Mpi::Broadcast(1, &expected_excluded, 0, Mpi::World());
  }
  REQUIRE(x_cell >= 1);
  REQUIRE(expected_excluded >= 1);
  REQUIRE(expected_excluded < patch_count);

  SECTION(
      "RecordRaw: exactly the void-cell patches are excluded, recorded, served raw once")
  {
    std::vector<std::unique_ptr<Mesh>> meshes;
    IoData iodata = Build(config, MakeSerialShortedCpwMeshWithVoid(x_cell), meshes);
    SpaceOperator space(iodata, meshes);
    SurfaceResponseOperator response(iodata, space);
    // The identification and the placement are unchanged by the void (the same patches).
    CHECK(response.GetPatchCount() == patch_count);
    const auto &fallbacks = response.GetGeometryFallbacks();
    REQUIRE(static_cast<int>(fallbacks.size()) == expected_excluded);
    std::set<int> excluded_patches;
    double recorded_length = 0.0;
    for (const auto &entry : fallbacks)
    {
      CHECK(entry.cls == "MaxwellPatchExcluded");
      CHECK(entry.site == "H14a/ResolveMaxwellLineCoverage");
      CHECK(entry.what == "patch");
      REQUIRE(entry.patches.size() == 1);
      excluded_patches.insert(entry.patches.front());
      CHECK(entry.models == std::vector<std::string>{"cpw-one-label-four-edge"});
      CHECK_FALSE(entry.topology.empty());  // the model's topology label
      CHECK_THAT(entry.reason, ContainsSubstring("leaves the finite-element mesh") &&
                                   ContainsSubstring("excluded (weight 0)"));
      CHECK(entry.length > 0.0);  // the length served raw = the raw claim's
      recorded_length += entry.length;
    }
    CHECK(static_cast<int>(excluded_patches.size()) == expected_excluded);
    if (Mpi::Root(Mpi::World()))
    {
      // The excluded set = the void-cell patches, weight 0 in the patch table; the others
      // keep their weight.
      for (const auto &assignment : response.GetPatchAssignments())
      {
        const bool inside = InCell(assignment.origin[0], x_cell);
        CHECK(inside == (excluded_patches.count(assignment.global_index) > 0));
        CHECK((inside ? assignment.weight == 0.0 : assignment.weight > 0.0));
      }
    }
    // Counted ONCE (decision 1074, the P == S form): the class-named portions = the growth
    // of the DomainBoundary portions = the geometry-fallback portions; every portion names
    // the class.
    const auto &h14a = response.GetMaxwellPatchExcludedPortionRecord();
    REQUIRE_FALSE(h14a.portions.empty());
    CHECK(h14a.length > 0.0);
    CHECK_THAT(recorded_length, WithinRel(h14a.length, 1.0e-12));
    CHECK_THAT(response.GetDomainBoundaryPortionRecord().length,
               WithinRel(intact_domain_boundary_length + h14a.length, 1.0e-12));
    double fallback_length = 0.0;
    for (const auto &portion : response.GetGeometryFallbackPortions())
    {
      CHECK(portion.topology == "MaxwellPatchExcluded");
      double length2 = 0.0;
      for (int d = 0; d < 3; d++)
      {
        length2 += (portion.p1[d] - portion.p0[d]) * (portion.p1[d] - portion.p0[d]);
      }
      fallback_length += std::sqrt(length2);
    }
    CHECK_THAT(fallback_length, WithinRel(h14a.length, 1.0e-12));
    for (const auto &portion : h14a.portions)
    {
      CHECK(portion.topology == "MaxwellPatchExcluded");
    }
    // The record the banner reads: the class, the count, the length share per Type.
    const json record = response.GetGeometryFallbackRecord();
    REQUIRE_FALSE(record.is_null());
    CHECK(record["Count"].get<int>() == expected_excluded);
    CHECK(record["Classes"].contains("MaxwellPatchExcluded"));
    CHECK(record["Classes"]["MaxwellPatchExcluded"].get<int>() == expected_excluded);
    CHECK_FALSE(response.IsRawFallback());
    // The conductor-consistency gate skips the excluded patches (weight 0) and the field
    // evaluation runs on the kept ones.
    int stub_records = 0;
    for (const auto &gate_record :
         response.ApplyConductorConsistencyGate(StripAtGround(space), 0, 5.0, false))
    {
      if (gate_record.probe_kind == 1)
      {
        stub_records++;
        CHECK(excluded_patches.count(gate_record.patch) == 0);
        CHECK_FALSE(gate_record.excluded);
      }
    }
    CHECK(stub_records == patch_count - expected_excluded);
    CHECK(response.GetPatchWeight() > 0.0);
  }

  SECTION("FailClosed: the configuration aborts at its end, naming the class")
  {
    json fail_closed = config;
    fail_closed["Solver"]["SurfaceResponseCorrection"]["GeometryFailurePolicy"] =
        "FailClosed";
    std::vector<std::unique_ptr<Mesh>> meshes;
    IoData iodata = Build(fail_closed, MakeSerialShortedCpwMeshWithVoid(x_cell), meshes);
    SpaceOperator space(iodata, meshes);
    CHECK_THROWS_WITH(SurfaceResponseOperator(iodata, space),
                      ContainsSubstring("GeometryFailurePolicy FailClosed") &&
                          ContainsSubstring("MaxwellPatchExcluded"));
  }
#endif
}

// R0-b's reader (idfix DESIGN A.3 (b); decision 1064 (4)): a model whose conductor
// references are one per LOCAL METAL PIECE loads with the piece carrier (Edges[].LocalPiece
// on an edge cluster, the top-level LocalPieces leaf on a gap), the piece count must equal
// the reference count, a ParallelEdgeCluster's pieces must be the ones its GapDirections
// bound (the ONE rule), and a model without the carrier loads exactly as today.
TEST_CASE_METHOD(
    test::SurfaceResponseFiles,
    "Fabrication-process library local-piece conductor references (R0-b reader)",
    "[surfaceresponseoperator][conductorconsistency][r0a][r0b][3d][Serial]")
{
#if !defined(MFEM_USE_GSLIB)
  SKIP("SurfaceResponseOperator requires MFEM_USE_GSLIB");
#else
  json base_library;
  {
    std::ifstream input(coupled_library_3d_path);
    REQUIRE(input);
    base_library = json::parse(input);
  }
  // The template gap model: two references and two open contour paths (the version-2
  // DifferentConductorGap form of the fixture).
  json gap_template;
  json cluster_template;
  for (const auto &model : base_library["Models"])
  {
    if (model.value("Topology", "") == "DifferentConductorGap" &&
        model.contains("ConductorReferences"))
    {
      gap_template = model;
    }
  }
  REQUIRE(!gap_template.is_null());
  {
    std::ifstream input(parallel_cluster_library_3d_path);
    REQUIRE(input);
    const json cluster_library = json::parse(input);
    for (const auto &model : cluster_library["Models"])
    {
      if (model.value("Name", "") == "cpw-three-conductor-four-edge-cluster")
      {
        cluster_template = model;
      }
    }
  }
  REQUIRE(!cluster_template.is_null());
  json config = ConvexMaxwellIslandConfig();
  auto &correction = config["Solver"]["SurfaceResponseCorrection"];
  correction["UnmatchedPolicy"] = "Warn";
  correction.erase("PatchConstruction");
  int variant = 0;
  auto Load = [&](json model, const std::string &label)
  {
    json library = base_library;
    library["Name"] = "unit-test-r0b-" + label;
    library["MatchingRadius"] = 0.2;  // the island config's EdgeDistances
    library["Models"] = {model};
    const auto library_path =
        temp.temp_dir / ("fabrication-process-r0b-" + std::to_string(variant++) + ".json");
    if (Mpi::Root(Mpi::World()))
    {
      std::ofstream output(library_path);
      output << library.dump(2) << "\n";
    }
    Mpi::Barrier(Mpi::World());
    auto run_config = config;
    run_config["Solver"]["SurfaceResponseCorrection"]["Library"] = library_path.string();
    IoData iodata(run_config, false);
    iodata.boundaries.cracked_attributes.insert(9);
    Mesh mesh(MakeIslandMesh());
    const auto manifest_path =
        temp.temp_dir / ("surface-response-requirements-r0b-" + label + ".json");
    WriteSurfaceResponseRequirements(iodata, mesh, manifest_path.string());
  };

  SECTION("a SameConductorGap without the leaf loads as today; with LocalPieces 2 and two "
          "references it loads; a mismatch is refused by name")
  {
    json same_gap = gap_template;
    same_gap["Name"] = "same-gap-local-pieces";
    same_gap["Topology"] = "SameConductorGap";
    // Without the leaf: plural references on a SameConductorGap are refused (today's rule).
    CHECK_THROWS_WITH(Load(same_gap, "scg-plural-no-leaf"),
                      ContainsSubstring("ConductorReferences") &&
                          ContainsSubstring("SameConductorGap with the LocalPieces leaf"));
    // Today's form: one Reference, no paths: loads.
    json legacy_gap = same_gap;
    legacy_gap.erase("ConductorReferences");
    legacy_gap.erase("OpenContourPaths");
    legacy_gap["Reference"] = {-13.0, 0.0, 0.0};
    CHECK_NOTHROW(Load(legacy_gap, "scg-legacy"));
    // With LocalPieces = 2 + 2 references + the open paths: loads.
    same_gap["LocalPieces"] = 2;
    CHECK_NOTHROW(Load(same_gap, "scg-local-pieces"));
    // A mismatched count is refused by name.
    json three = same_gap;
    three["LocalPieces"] = 3;
    CHECK_THROWS_WITH(Load(three, "scg-three"), ContainsSubstring("LocalPieces must be 2"));
    json one_reference = same_gap;
    one_reference["ConductorReferences"] = {{-13.0, 0.0, 0.0}};
    CHECK_THROWS_WITH(Load(one_reference, "scg-one-reference"),
                      ContainsSubstring("at least two points"));
    // The leaf on another topology is refused.
    json on_cluster = cluster_template;
    on_cluster["LocalPieces"] = 3;
    CHECK_THROWS_WITH(Load(on_cluster, "pec-leaf"),
                      ContainsSubstring("LocalPieces is the SameConductorGap"));
  }

  SECTION("a ParallelEdgeCluster with Edges[].LocalPiece: one reference per piece, the "
          "pieces those of its GapDirections")
  {
    // The three-conductor four-edge template (labels 1 / 2 / 2 / 3, three references) read
    // as a one-label CPW with three LOCAL pieces (ground / strip / ground).
    json cluster = cluster_template;
    cluster["Name"] = "cpw-one-label-three-pieces";
    const std::vector<int> pieces = {1, 2, 2, 3};
    for (std::size_t k = 0; k < 4; k++)
    {
      cluster["Edges"][k]["Conductor"] = 1;
      cluster["Edges"][k]["LocalPiece"] = pieces[k];
    }
    cluster["OpenContourPaths"] = {
        {{"Indices", {1, 2}}, {"StartConductor", 1}, {"EndConductor", 2}},
        {{"Indices", {3, 4}}, {"StartConductor", 2}, {"EndConductor", 3}}};
    CHECK_NOTHROW(Load(cluster, "pec-local-pieces"));
    // Without the leaf a one-label cluster must carry ONE reference (today's rule).
    json no_leaf = cluster;
    for (std::size_t k = 0; k < 4; k++)
    {
      no_leaf["Edges"][k].erase("LocalPiece");
    }
    CHECK_THROWS_WITH(Load(no_leaf, "pec-no-leaf"),
                      ContainsSubstring("one conductor reference for each canonical "
                                        "conductor label"));
    // The leaf on some edges only is refused.
    json partial = cluster;
    partial["Edges"][2].erase("LocalPiece");
    CHECK_THROWS_WITH(Load(partial, "pec-partial"),
                      ContainsSubstring("every edge or on none"));
    // Pieces that are not the ones the GapDirections bound are refused (the ONE rule).
    json wrong = cluster;
    wrong["Edges"][2]["LocalPiece"] = 3;
    wrong["Edges"][3]["LocalPiece"] = 4;
    CHECK_THROWS_WITH(Load(wrong, "pec-wrong-pieces"),
                      ContainsSubstring("not the piece its GapDirections bound"));
    // Fewer references than pieces is refused by name.
    json two_references = cluster;
    two_references["ConductorReferences"] = {{0.0, 0.0, 0.0}, {12.0, 0.0, 0.0}};
    two_references["OpenContourPaths"] = {
        {{"Indices", {1, 2}}, {"StartConductor", 1}, {"EndConductor", 2}}};
    CHECK_THROWS_WITH(Load(two_references, "pec-two-references"),
                      ContainsSubstring("one reference per local piece is required"));
  }
#endif
}

TEST_CASE("Response-correction ConductorConsistency.RecordOnly is read and bounded",
          "[surfaceresponseoperator][conductorconsistency][r0a][config][Serial]")
{
  json config = {{"Problem", {{"Type", "Electrostatic"}, {"Output", "unused"}}},
                 {"Model", {{"Mesh", "unused.msh"}}},
                 {"Domains", {{"Materials", {{{"Attributes", {1}}}}}}},
                 {"Boundaries", {{"Ground", {{"Attributes", {1}}}}}},
                 {"Solver",
                  {{"Electrostatic",
                    {{"ResponseCorrection",
                      {{"Library", "unused.json"},
                       {"ConductorConsistency", {{"RecordOnly", true}}}}}}}}}};
  {
    IoData iodata(config, false);
    REQUIRE(iodata.solver.electrostatic.response_correction);
    CHECK(
        iodata.solver.electrostatic.response_correction->conductor_consistency_record_only);
  }
  config["Solver"]["Electrostatic"]["ResponseCorrection"].erase("ConductorConsistency");
  {
    IoData iodata(config, false);
    CHECK_FALSE(
        iodata.solver.electrostatic.response_correction->conductor_consistency_record_only);
  }
  config["Solver"]["Electrostatic"]["ResponseCorrection"]["ConductorConsistency"] = {
      {"Tolerance", 0.05}};
  CHECK_THROWS_WITH(IoData(config, false), ContainsSubstring("unknown key \"Tolerance\""));
}

}  // namespace palace
