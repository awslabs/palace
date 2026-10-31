// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "fixtures.hpp"
#include "surfaceresponse-fixtures.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iterator>
#include <map>
#include <memory>
#include <set>
#include <sstream>
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
using ResponseCorrectionData = config::ElectrostaticSolverData::ResponseCorrectionData;
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
// v = -0.1 below it), diagonal response matrices; `reference` = the model's Reference (the
// default on the strip; a reference in a slot is the ReferenceInGap fixture). Written by
// the root; returns the path.
fs::path WriteStubLibrary(const fs::path &directory, const fs::path &template_library,
                          double ring_top_v, const std::string &suffix,
                          std::array<double, 3> reference = {0.25, 0.0, 0.0})
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
    model["Reference"] = {reference[0], reference[1], reference[2]};
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

// The OPEN-ENDED CPW stub of the 1064 (2) spatial fixture (decision 1082 Q3): the shorted
// CPW stub's plane (AddShortedCpwMetalPlane's slots, x from 0.5 into the wall) on a 3 x 1 x
// 2 box whose STRIP ENDS at x = 2.25 — beyond it the two slots merge into one bare region
// z in [0.75, 1.25] — so the strip's two convex end corners with the ground edges 0.125
// away on both sides form a SpatialEdgeCluster at x = 2.25 holding THREE metal pieces
// (ground-left, strip, ground-right) in its box; the 4-edge stack x in (0.5, 2.25) stays
// the stub's translational PEC feature (1.75 long: the end windows leave a span) and the
// short's corner cluster at x = 0.5 is as in every stub case. `island` false = the strip is
// shorted to the ground at x = 0.5 (ONE metal component: the end cluster has ONE canonical
// label, the records' one-label multi-piece form); `island` true = the strip starts at
// x = 0.75 as its own component (the end cluster has TWO labels, the junction-lead form).
// The ground edges beyond the end face each other at 0.5 > 2R (no pair); the end is 0.75
// from the wall.
mfem::Mesh MakeSerialOpenCpwMesh(bool island)
{
  const double strip_begin = island ? 0.75 : 0.5;
  constexpr double strip_end = 2.25;
  mfem::Mesh serial =
      mfem::Mesh::MakeCartesian3D(24, 4, 16, mfem::Element::HEXAHEDRON, 3.0, 1.0, 2.0);
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
    double xmin = 3.0, xmax = 0.0, zmin = 2.0, zmax = 0.0;
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
    const bool in_strip_band = zmin >= 0.875 - 1.0e-12 && zmax <= 1.125 + 1.0e-12;
    const bool beyond_strip = in_strip_band && xmin >= strip_end - 1.0e-12;
    const bool before_strip =
        in_strip_band && xmin >= 0.5 - 1.0e-12 && xmax <= strip_begin + 1.0e-12;
    if ((in_slot_x && (in_slot_a || in_slot_b)) || beyond_strip || before_strip)
    {
      continue;  // bare dielectric
    }
    serial.AddBdrElement(serial.GetFace(face)->Duplicate(&serial));
    serial.SetBdrAttribute(serial.GetNBE() - 1, 9);
  }
  serial.FinalizeTopology();
  serial.Finalize();
  return serial;
}

// The strip-end cluster as the identification keys it: its Features[] entry of the
// requirements manifest (Signature / Frame), read by an electrostatic geometry-only
// preflight of the open-ended CPW mesh against a library without a spatial model (the
// decision-632 placeholder form of test-conductorconsistency's Signature-keyed case). The
// ONE SpatialEdgeCluster feature whose frame origin lies at x > 1.5 (the strip end at
// x = 2.25; the short's / the island-start cluster lies at x <= 0.75; the stack between
// them is translational).
json IdentifyStripEndCluster(const test::SurfaceResponseFiles &files,
                             const fs::path &library, const fs::path &directory,
                             bool island)
{
  json config = files.IslandConfig();
  auto &correction = config["Solver"]["Electrostatic"]["ResponseCorrection"];
  correction.erase("PatchConstruction");
  correction["Library"] = library.string();
  correction["UnmatchedPolicy"] = "Warn";
  IoData iodata(config, false);
  iodata.boundaries.cracked_attributes.insert(9);
  mfem::Mesh serial = MakeSerialOpenCpwMesh(island);
  Mesh mesh(std::make_unique<mfem::ParMesh>(Mpi::World(), serial));
  if (Mpi::Root(Mpi::World()))
  {
    fs::create_directories(directory);
  }
  Mpi::Barrier(Mpi::World());
  const auto manifest_path = directory / "surface-response-requirements.json";
  WriteSurfaceResponseRequirements(iodata, mesh, manifest_path.string());
  Mpi::Barrier(Mpi::World());
  std::ifstream input(manifest_path);
  REQUIRE(input);
  const json manifest = json::parse(input);
  json end_cluster;
  int end_clusters = 0;
  for (const auto &feature : manifest["Identification"]["Features"])
  {
    if (feature["Type"] == "SpatialEdgeCluster" &&
        feature["Frame"]["Origin"][0].get<double>() > 1.5)
    {
      end_cluster = feature;
      end_clusters++;
    }
  }
  REQUIRE(end_clusters == 1);
  REQUIRE(end_cluster["Match"]["Status"] == "Missing");
  // The key the identification gave the cluster (information for the log): its frame and
  // its support box (units of R; the face rules grow it along the slots).
  Mpi::Print("Strip-end cluster keyed: Frame {} Box {}\n", end_cluster["Frame"].dump(),
             end_cluster["Signature"].value("Box", json()).dump());
  return end_cluster;
}

// The frame of an identified feature (Features[].Frame: Origin + Axes u / v / n, mesh units
// here: the unit tests' coordinate scale is 1): device point -> model-frame (u, v, w).
struct FeatureFrame
{
  std::array<double, 3> origin{};
  std::array<std::array<double, 3>, 3> axes{};
  explicit FeatureFrame(const json &frame)
  {
    origin = frame["Origin"].get<std::array<double, 3>>();
    axes = frame["Axes"].get<std::array<std::array<double, 3>, 3>>();
  }
  std::array<double, 3> Local(const std::array<double, 3> &point) const
  {
    std::array<double, 3> local{};
    for (int a = 0; a < 3; a++)
    {
      for (int d = 0; d < 3; d++)
      {
        local[a] += (point[d] - origin[d]) * axes[a][d];
      }
    }
    return local;
  }
};

// The support box of an identified spatial cluster, [u0, v0, u1, v1] in units of R: the
// hashed Signature's Box (a contract-3 key) or, for a claims-only key, the SpatialSupport
// record's.
std::array<double, 4> FeatureBoxOverR(const json &feature)
{
  if (feature["Signature"].contains("Box"))
  {
    return feature["Signature"]["Box"].get<std::array<double, 4>>();
  }
  REQUIRE(feature.contains("SpatialSupport"));
  REQUIRE(feature["SpatialSupport"].contains("Box"));
  return feature["SpatialSupport"]["Box"].get<std::array<double, 4>>();
}

// The spatial-cluster coupon of the strip-end cluster, keyed by the identification's own
// Signature (the model's frame = the feature's canonical frame): BasisPoints = an 8-knot
// ring at w = +0.1 in a 1.5 R square about the frame origin (inside the mesh whatever the
// frame's in-plane orientation; the identification's support box is far larger here, its
// face rules grow it over the slots); an explicit TraceMesh whose conductor vertices are
// the given device points on the plane w = 0 (`ground_points` on the two grounds' metal,
// `strip_points` on the strip's); `two_references` false = ONE reference (`ground_anchor`)
// and every conductor vertex 'conductor 1' (the records' one-label multi-piece form); true
// = ConductorReferences {ground_anchor, strip_anchor}, the strip vertices 'conductor 2' and
// two OpenContourPaths (knots 1-4 and 8-5, conductor 1 -> 2). Diagonal response matrices.
// Written by the root; returns the library path.
fs::path WriteStripEndClusterLibrary(
    const fs::path &directory, const fs::path &template_library, const json &feature,
    const std::vector<std::array<double, 3>> &ground_points,
    const std::vector<std::array<double, 3>> &strip_points,
    const std::array<double, 3> &ground_anchor, const std::array<double, 3> &strip_anchor,
    bool two_references, const std::string &suffix)
{
  constexpr double R = 0.2, ring_w = 0.1, half = 1.5 * R;
  const FeatureFrame frame(feature["Frame"]);
  const std::vector<std::array<double, 2>> ring = {
      {-half, -half}, {0.0, -half}, {half, -half}, {half, 0.0},
      {half, half},   {0.0, half},  {-half, half}, {-half, 0.0}};
  const int basis_size = static_cast<int>(ring.size()) + (two_references ? 1 : 0);
  const auto library_path =
      directory / fmt::format("fabrication-process-r0a-stripend-3d-{}.json", suffix);
  const auto points_path =
      directory / fmt::format("r0a-stripend-basis-points-{}.csv", suffix);
  const auto vertices_path =
      directory / fmt::format("r0a-stripend-trace-vertices-{}.csv", suffix);
  const auto triangles_path =
      directory / fmt::format("r0a-stripend-trace-triangles-{}.csv", suffix);
  const auto domain_path = directory / fmt::format("r0a-stripend-domain-{}.csv", suffix);
  const auto surface_path = directory / fmt::format("r0a-stripend-surface-{}.csv", suffix);
  if (Mpi::Root(Mpi::World()))
  {
    {
      std::ofstream output(points_path);
      output << "x,y,z\n";
      for (const auto &knot : ring)
      {
        output << fmt::format("{:.16e},{:.16e},{:.16e}\n", knot[0], knot[1], ring_w);
      }
    }
    {
      std::ofstream output(vertices_path);
      output << "vertex,x,y,z,basis,conductor\n";
      int index = 0;
      for (const auto &knot : ring)
      {
        index++;
        output << fmt::format("{:d},{:.16e},{:.16e},{:.16e},{:d},0\n", index, knot[0],
                              knot[1], ring_w, index);
      }
      for (const auto &point : ground_points)
      {
        const auto local = frame.Local(point);
        output << fmt::format("{:d},{:.16e},{:.16e},{:.16e},0,1\n", ++index, local[0],
                              local[1], local[2]);
      }
      for (const auto &point : strip_points)
      {
        const auto local = frame.Local(point);
        output << fmt::format("{:d},{:.16e},{:.16e},{:.16e},0,{:d}\n", ++index, local[0],
                              local[1], local[2], two_references ? 2 : 1);
      }
    }
    {
      // A fan of the ring's edges onto the first conductor vertex (vertex 9).
      std::ofstream output(triangles_path);
      output << "triangle,vertex_i,vertex_j,vertex_k\n";
      for (std::size_t i = 0; i < ring.size(); i++)
      {
        output << fmt::format("{:d},{:d},{:d},9\n", i + 1, i + 1,
                              (i + 1) % ring.size() + 1);
      }
    }
    {
      std::ofstream output(domain_path);
      output << "basis_i,basis_j,Q_ij (J)\n";
      for (int i = 1; i <= basis_size; i++)
      {
        for (int j = i; j <= basis_size; j++)
        {
          output << i << "," << j << "," << (i == j ? 1.0e-12 : 0.0) << "\n";
        }
      }
    }
    {
      std::ofstream output(surface_path);
      output << "interface,edge,R (m),basis_i,basis_j,Q_ij (J),Q_total_ij (J)\n";
      for (int i = 1; i <= basis_size; i++)
      {
        for (int j = i; j <= basis_size; j++)
        {
          output << "1,1," << R << "," << i << "," << j << "," << (i == j ? 0.25e-12 : 0.0)
                 << "," << (i == j ? 0.25e-12 : 0.0) << "\n";
        }
      }
    }
    std::ifstream input(template_library);
    REQUIRE(input);
    json library = json::parse(input);
    library["Name"] = "unit-test-process-r0a-stripend-3d-" + suffix;
    json model = library["Models"][0];
    model["Name"] = "strip-end-cluster-" + suffix;
    model["Topology"] = "SpatialEdgeCluster";
    model["Interfaces"] = {{{"Type", "SA"}, {"Coupon", 1}}};
    model["Signature"] = feature["Signature"];
    model["BasisPoints"] = points_path.string();
    model["TraceMesh"] = {{"Vertices", vertices_path.string()},
                          {"Triangles", triangles_path.string()}};
    model["FabricatedMatrix"] = domain_path.string();
    model["ThinMatrix"] = domain_path.string();
    model["FabricatedSurfaceMatrix"] = surface_path.string();
    model["ThinSurfaceMatrix"] = surface_path.string();
    for (const char *key :
         {"Separation", "SeparationTolerance", "Angle", "AngleTolerance", "Edges",
          "Reference", "ConductorReferences", "OpenContourPaths", "ContourGroups"})
    {
      model.erase(key);
    }
    const auto ground = frame.Local(ground_anchor);
    if (two_references)
    {
      const auto strip = frame.Local(strip_anchor);
      model["ConductorReferences"] = {{ground[0], ground[1], 0.0},
                                      {strip[0], strip[1], 0.0}};
      model["OpenContourPaths"] = {
          {{"Indices", {1, 2, 3, 4}}, {"StartConductor", 1}, {"EndConductor", 2}},
          {{"Indices", {8, 7, 6, 5}}, {"StartConductor", 1}, {"EndConductor", 2}}};
    }
    else
    {
      model["Reference"] = {ground[0], ground[1], 0.0};
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

  SECTION("(g) the probe path is the end-face w = 0 leg, not a chord: a ramp V(x) = x")
  {
    // B.3 (g) / C.3 (decision 1082 Q4): with E_z = -x / 0.125 in slot A (+x / 0.125 in
    // slot B) at EVERY x (no jump, so the field is exactly representable in ND1 on the
    // hexes: E_z linear in x; curl E || y, zero flux through the ring's x = const plane),
    // the strip probe on the end face at x_e reads P = V(x_e) = x_e along the straight
    // w = 0 leg across slot A at that face, while a straight CHORD from the anchor at the
    // patch origin's cross-section (x_c) to the same probe point would read the mean
    // (x_c + x_e) / 2: the record's D equals x_e to roundoff and differs from the chord by
    // |x_e - x_c| / 2, so 'at the probe's v' is load-bearing (the rule discriminates).
    GridFunction field(space.GetNDSpace(), true);
    mfem::VectorFunctionCoefficient coefficient(3,
                                                [](const mfem::Vector &x, mfem::Vector &E)
                                                {
                                                  E.SetSize(3);
                                                  E = 0.0;
                                                  if (x[2] > 0.75 && x[2] < 0.875)
                                                  {
                                                    E[2] = -x[0] / 0.125;
                                                  }
                                                  else if (x[2] > 1.125 && x[2] < 1.25)
                                                  {
                                                    E[2] = x[0] / 0.125;
                                                  }
                                                });
    field.Real().ProjectCoefficient(coefficient);
    field.Imag() = 0.0;
    ComplexVector E(field.Real().ParFESpace()->GetTrueVSize());
    field.Real().GetTrueDofs(E.Real());
    field.Imag().GetTrueDofs(E.Imag());
    std::map<int, double> origin_x;  // the root's patch table
    if (Mpi::Root(Mpi::World()))
    {
      for (const auto &assignment : response.GetPatchAssignments())
      {
        origin_x[assignment.global_index] = assignment.origin[0];
      }
    }
    const auto records =
        StubRecords(response.ApplyConductorConsistencyGate(E, 0, 5.0, false));
    REQUIRE(records.size() == static_cast<std::size_t>(patch_count));
    int discriminated = 0;
    for (const auto &record : records)
    {
      // The worst probe: a strip edge (1 / 2, piece 2) at one end face, on the plane.
      CHECK(record.worst_piece == 2);
      CHECK((record.worst_edge == 1 || record.worst_edge == 2));
      CHECK((record.worst_end_face == 0 || record.worst_end_face == 1));
      CHECK_THAT(record.worst_point[1], WithinAbs(0.5, 1.0e-9));
      CHECK((std::abs(record.worst_point[2] - 0.875) < 1.0e-9 ||
             std::abs(record.worst_point[2] - 1.125) < 1.0e-9));
      const double x_e = record.worst_point[0];
      CHECK(x_e > 0.5);
      // D = V(x_e) along the end-face leg (exactly, the field being in the space).
      CHECK_THAT(record.max_deviation, WithinRel(x_e, 1.0e-9));
      CHECK(record.excluded);
      if (Mpi::Root(Mpi::World()))
      {
        REQUIRE(origin_x.count(record.patch) == 1);
        const double x_c = origin_x.at(record.patch);
        const double chord = 0.5 * (x_c + x_e);
        // The end face is off the origin's cross-section by a finite part of the cell.
        CHECK(std::abs(x_e - x_c) > 1.0e-3);
        CHECK(std::abs(x_e - x_c) <= record.cell_length + 1.0e-9);
        CHECK_THAT(std::abs(record.max_deviation - chord),
                   WithinRel(0.5 * std::abs(x_e - x_c), 1.0e-6));
        CHECK(std::abs(record.max_deviation - chord) > 1.0e-3);
        discriminated++;
      }
    }
    if (Mpi::Root(Mpi::World()))
    {
      CHECK(discriminated == patch_count);
    }
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
    // The FloorApplied COUNT leaf (decision 1082 Q2; B.3 (f) / C.4: read per record form):
    // the floored records of this pass (every one kept here).
    const auto statistics = response.GetStatistics();
    const auto &diagnostics = statistics["Diagnostics"]["ConductorConsistency"];
    CHECK(diagnostics["FloorAppliedPatches"].get<int>() == floored_kept);
    CHECK(diagnostics["TestedPatches"].get<int>() == patch_count);
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

// ReferenceInGap (idfix R0-a; decision 1082 Q1, the d1073 review's MINOR-3): a model
// whose conductor reference lies in a GAP of its own cross-section is a library defect
// recorded by name on BOTH paths at the placement (ONE rule, the GeometryFailurePolicy
// form): under RecordRaw the feature is not placed (its claims raw) and recorded; under
// FailClosed the configuration aborts at its end naming the class and the model. A
// reference ON an edge (inclusive) or on metal is placed as before.
TEST_CASE_METHOD(test::SurfaceResponseFiles,
                 "SurfaceResponseOperator records a model whose conductor reference lies "
                 "in a gap of its cross-section (ReferenceInGap, both paths)",
                 "[surfaceresponseoperator][geometryfallback][conductorconsistency][r0a]"
                 "[3d][Serial][Parallel]")
{
#if !defined(MFEM_USE_GSLIB)
  SKIP("SurfaceResponseOperator requires MFEM_USE_GSLIB");
#else
  // The Maxwell path on the shorted CPW stub: the one-label 4-edge model (edges at 0 /
  // 0.125 / 0.375 / 0.5) with its Reference in slot A (u = 0.0625).
  json config = ConvexMaxwellIslandConfig();
  config["Boundaries"]["Ground"]["Attributes"] = {1, 2, 3, 4, 5, 6, 9};
  auto &correction = config["Solver"]["SurfaceResponseCorrection"];
  correction.erase("PatchConstruction");
  correction["UnmatchedPolicy"] = "Warn";
  correction["TraceCoupling"] = "Collocated";
  auto Build = [&](json run_config, const fs::path &library,
                   std::vector<std::unique_ptr<Mesh>> &meshes)
  {
    run_config["Solver"]["SurfaceResponseCorrection"]["Library"] = library.string();
    IoData iodata(run_config, false);
    iodata.boundaries.cracked_attributes.insert(9);
    mfem::Mesh serial = MakeSerialShortedCpwMesh();
    meshes.push_back(
        std::make_unique<Mesh>(std::make_unique<mfem::ParMesh>(Mpi::World(), serial)));
    return iodata;
  };
  const auto in_gap = WriteStubLibrary(temp.temp_dir, convex_library_3d_path, 0.1,
                                       "reference-in-gap", {0.0625, 0.0, 0.0});
  const auto on_edge = WriteStubLibrary(temp.temp_dir, convex_library_3d_path, 0.1,
                                        "reference-on-edge", {0.125, 0.0, 0.0});

  SECTION("RecordRaw (a solve): the feature is not placed, its claims raw, recorded")
  {
    std::vector<std::unique_ptr<Mesh>> meshes;
    IoData iodata = Build(config, in_gap, meshes);
    SpaceOperator space(iodata, meshes);
    SurfaceResponseOperator response(iodata, space);
    // The only matched model is the defective one: nothing is placed (the short's corner
    // cluster is unmatched, raw, as in every stub case), so the run also carries the H1
    // whole-run NoMatchedFeature record; the ReferenceInGap record(s) name the feature.
    CHECK(response.GetPatchCount() == 0);
    std::vector<ResponseCorrectionData::GeometryFallbackData> fallbacks;
    bool no_matched_feature = false;
    for (const auto &entry : response.GetGeometryFallbacks())
    {
      if (entry.cls == "NoMatchedFeature")
      {
        no_matched_feature = true;
        continue;
      }
      fallbacks.push_back(entry);
    }
    CHECK(no_matched_feature);
    REQUIRE_FALSE(fallbacks.empty());
    for (const auto &entry : fallbacks)
    {
      CHECK(entry.cls == "ReferenceInGap");
      CHECK(entry.site == "R0a/BuildAutomaticResponseData3D:reference-in-gap");
      CHECK(entry.what == "feature");
      CHECK(entry.models == std::vector<std::string>{"cpw-one-label-four-edge"});
      CHECK(entry.topology == "ParallelEdgeCluster");
      CHECK(entry.length > 0.0);
      CHECK_THAT(entry.reason,
                 ContainsSubstring("conductor reference 1 at lateral offset 6.25") &&
                     ContainsSubstring("lies in a gap of its cross-section") &&
                     ContainsSubstring("cpw-one-label-four-edge"));
    }
    const auto record = response.GetStatistics()["Diagnostics"]["GeometryFallbacks"];
    CHECK(record["Classes"]["ReferenceInGap"].get<int>() ==
          static_cast<int>(fallbacks.size()));
  }

  SECTION("FailClosed: the configuration aborts at its end, naming the class and model")
  {
    json fail_closed = config;
    fail_closed["Solver"]["SurfaceResponseCorrection"]["GeometryFailurePolicy"] =
        "FailClosed";
    std::vector<std::unique_ptr<Mesh>> meshes;
    IoData iodata = Build(fail_closed, in_gap, meshes);
    SpaceOperator space(iodata, meshes);
    CHECK_THROWS_WITH(SurfaceResponseOperator(iodata, space),
                      ContainsSubstring("GeometryFailurePolicy FailClosed") &&
                          ContainsSubstring("ReferenceInGap"));
  }

  SECTION("a reference ON an edge lies on metal (inclusive): placed, no record")
  {
    std::vector<std::unique_ptr<Mesh>> meshes;
    IoData iodata = Build(config, on_edge, meshes);
    SpaceOperator space(iodata, meshes);
    SurfaceResponseOperator response(iodata, space);
    CHECK(response.GetPatchCount() > 0);
    CHECK(response.GetGeometryFallbacks().empty());
  }

  SECTION("the electrostatic path shares the rule: the island's isolated-edge model with "
          "its reference in the gap (u > 0) is recorded, not placed")
  {
    const auto library_path = temp.temp_dir / "fabrication-process-r0a-reference-gap.json";
    if (Mpi::Root(Mpi::World()))
    {
      std::ifstream input(convex_library_3d_path);
      REQUIRE(input);
      json library = json::parse(input);
      library["Name"] = "unit-test-process-r0a-reference-gap-3d";
      for (auto &model : library["Models"])
      {
        model["Interfaces"] = {{{"Type", "SA"}, {"Coupon", 1}}};
      }
      REQUIRE(library["Models"][0]["Topology"] == "IsolatedEdge");
      library["Models"][0]["Reference"] = {0.1, 0.0, 0.0};  // into the gap (+u)
      std::ofstream output(library_path);
      output << library.dump(2) << "\n";
    }
    Mpi::Barrier(Mpi::World());
    json island = IslandConfig();
    auto &island_correction = island["Solver"]["Electrostatic"]["ResponseCorrection"];
    island_correction.erase("PatchConstruction");
    island_correction["Library"] = library_path.string();
    island_correction["TraceCoupling"] = "SurfaceMortar";
    island_correction["MortarOversampling"] = 2;
    // A rejected feature is unmatched (raw): the run's UnmatchedPolicy applies to it as to
    // every Missing feature (Error would abort the run by that policy, not by this rule).
    island_correction["UnmatchedPolicy"] = "Warn";
    IoData iodata(island, false);
    iodata.boundaries.cracked_attributes.insert(9);
    std::vector<std::unique_ptr<Mesh>> meshes;
    meshes.push_back(std::make_unique<Mesh>(MakeIslandMesh()));
    LaplaceOperator laplace(iodata, meshes);
    SurfaceResponseOperator response(iodata, laplace);
    const auto &fallbacks = response.GetGeometryFallbacks();
    REQUIRE_FALSE(fallbacks.empty());
    for (const auto &entry : fallbacks)
    {
      CHECK(entry.cls == "ReferenceInGap");
      CHECK(entry.topology == "IsolatedEdge");
      CHECK(entry.models.front() == "isolated");
    }
  }
#endif
}

// U-R0a's SPATIAL case (decision 1064 (2); decision 1082 Q3; idfix DESIGN A.3 (c) / B.3 /
// C.3): the Maxwell probes of a spatial model = its TraceMesh's plane conductor vertices,
// each read along ONE straight plane leg from its OWN conductor's anchor. On the strip-end
// cluster of the open-ended CPW stub (three metal pieces in its box: ground-left, the
// strip, ground-right; the records' one-label multi-piece class, A.3 (c)) under the stub's
// synthetic field with the strip at V against both grounds (E_z across the two slots): G1 =
// the one-label cluster served by a ONE-reference model (the anchor on ground-left): the
// strip vertices are read from the ground anchor across slot A, D ~ V, A_p ~ V -> EXCLUDED,
// FloorApplied false (the ground-right vertices read V - V = 0 across both slots); G1b =
// the same model with the strip AT the ground potential (a field normal to the plane only):
// kept, D ~ 0; G2 = the strip as an island (two labels) served by a TWO-reference model
// {ground, strip}: every vertex is read from its own piece's anchor, D ~ 0 -> KEPT, the
// conductor state = V. The models are keyed by the identification's own Signature and
// placed in its Frame (the decision-632 placeholder form), so no hand key enters the test.
TEST_CASE_METHOD(test::SurfaceResponseFiles,
                 "SurfaceResponseOperator conductor-consistency gate on a spatial cluster "
                 "(R0-a, the one-label multi-piece strip end; decision 1064 (2))",
                 "[surfaceresponseoperator][conductorconsistency][r0a][maxwell][spatial]"
                 "[3d][Serial][Parallel]")
{
#if !defined(MFEM_USE_GSLIB)
  SKIP("SurfaceResponseOperator requires MFEM_USE_GSLIB");
#else
  constexpr double V = 1.0;
  // The identification library (the stub's one-label PEC model: no spatial model, so the
  // strip-end cluster is Missing and keyed; the 4-edge stack matches it as in every stub
  // case).
  const auto identification_library =
      WriteStubLibrary(temp.temp_dir, convex_library_3d_path, 0.1, "stripend-id");
  // The device points of the coupon's metal (mesh units; all on the plane y = 0.5): three
  // on ground-left (z < 0.75), three on ground-right (z > 1.25), three on the strip (z in
  // (0.875, 1.125), x < 2.25); the anchors on ground-left and on the strip.
  const std::vector<std::array<double, 3>> ground_points = {
      {2.125, 0.5, 0.6875}, {2.25, 0.5, 0.6875}, {2.375, 0.5, 0.6875},
      {2.125, 0.5, 1.3125}, {2.25, 0.5, 1.3125}, {2.375, 0.5, 1.3125}};
  const std::vector<std::array<double, 3>> strip_points = {
      {2.1875, 0.5, 1.0}, {2.125, 0.5, 0.9375}, {2.125, 0.5, 1.0625}};
  const std::array<double, 3> ground_anchor = {2.15625, 0.5, 0.6875};
  const std::array<double, 3> strip_anchor = {2.15625, 0.5, 1.0};
  auto OnPlane = [](const json &feature, const std::array<double, 3> &point)
  {
    const FeatureFrame frame(feature["Frame"]);
    return std::abs(frame.Local(point)[2]) < 1.0e-9;
  };
  json config = ConvexMaxwellIslandConfig();
  config["Boundaries"]["Ground"]["Attributes"] = {1, 2, 3, 4, 5, 6, 9};
  auto &correction = config["Solver"]["SurfaceResponseCorrection"];
  correction.erase("PatchConstruction");
  correction["UnmatchedPolicy"] = "Warn";
  correction["TraceCoupling"] = "Collocated";
  auto Build =
      [&](const fs::path &library, bool island, std::vector<std::unique_ptr<Mesh>> &meshes)
  {
    json run_config = config;
    run_config["Solver"]["SurfaceResponseCorrection"]["Library"] = library.string();
    IoData iodata(run_config, false);
    iodata.boundaries.cracked_attributes.insert(9);
    mfem::Mesh serial = MakeSerialOpenCpwMesh(island);
    meshes.push_back(
        std::make_unique<Mesh>(std::make_unique<mfem::ParMesh>(Mpi::World(), serial)));
    return iodata;
  };
  // The strip at V against both grounds on x in [1.75, 3] (the end cluster's ring and
  // beyond; StubField: E_z = -V / 0.125 in slot A, +V / 0.125 in slot B, zero elsewhere:
  // curl-free where the ring lies). `normal` = a field normal to the plane only (the strip
  // at the ground potential: nothing across the slots).
  auto Field = [](SpaceOperator &space, bool normal)
  {
    GridFunction field(space.GetNDSpace(), true);
    StubField stub;
    stub.amplitude = V;
    stub.x0 = 1.75;
    stub.x1 = 3.0;
    mfem::VectorFunctionCoefficient coefficient(
        3,
        [normal, &stub](const mfem::Vector &x, mfem::Vector &E)
        {
          if (normal)
          {
            E.SetSize(3);
            E = 0.0;
            E[1] = 1.0;
            return;
          }
          stub(x, E);
        });
    field.Real().ProjectCoefficient(coefficient);
    field.Imag() = 0.0;
    ComplexVector E(field.Real().ParFESpace()->GetTrueVSize());
    field.Real().GetTrueDofs(E.Real());
    field.Imag().GetTrueDofs(E.Imag());
    return E;
  };
  auto SpatialRecords =
      [](const std::vector<SurfaceResponseOperator::ConductorConsistencyRecord> &records)
  {
    std::vector<SurfaceResponseOperator::ConductorConsistencyRecord> spatial;
    for (const auto &record : records)
    {
      if (record.probe_kind == 0)
      {
        spatial.push_back(record);
      }
    }
    return spatial;
  };
  const auto tolerance = SurfaceResponseOperator::kConductorConsistencyTolerance;

  SECTION("G1 / G1b: one label, three pieces, one reference")
  {
    const json feature = IdentifyStripEndCluster(*this, identification_library,
                                                 temp.temp_dir / "stripend-g1", false);
    for (const auto &point : ground_points)
    {
      REQUIRE(OnPlane(feature, point));
    }
    for (const auto &point : strip_points)
    {
      REQUIRE(OnPlane(feature, point));
    }
    const auto library = WriteStripEndClusterLibrary(
        temp.temp_dir, convex_library_3d_path, feature, ground_points, strip_points,
        ground_anchor, strip_anchor, false, "g1");
    // G1: the strip at V -> the strip vertices read V from the ground anchor: EXCLUDED.
    {
      std::vector<std::unique_ptr<Mesh>> meshes;
      IoData iodata = Build(library, false, meshes);
      SpaceOperator space(iodata, meshes);
      SurfaceResponseOperator response(iodata, space);
      REQUIRE(response.HasConductorConsistencyProbes());
      const auto records = SpatialRecords(
          response.ApplyConductorConsistencyGate(Field(space, false), 0, 5.0, false));
      REQUIRE(records.size() == 1);
      const auto &record = records.front();
      CHECK(record.maxwell);
      CHECK(record.plane_knots == 9);
      CHECK(record.state == 0.0);  // a one-reference model has no conductor state
      CHECK(record.worst_conductor == 1);
      // The worst probe: a strip vertex, read across slot A from the ground anchor.
      CHECK(record.worst_point[2] > 0.875);
      CHECK(record.worst_point[2] < 1.125);
      CHECK_THAT(record.worst_point[1], WithinAbs(0.5, 1.0e-9));
      CHECK_THAT(record.max_deviation, WithinRel(V, 0.2));
      CHECK(record.max_ratio > 10 * tolerance);
      CHECK(record.normalization == record.amplitude);  // FloorApplied false
      CHECK(record.excluded);
      CHECK_FALSE(response.GetConductorConsistencyPortions().empty());
      const auto statistics = response.GetStatistics();
      const auto &diagnostics = statistics["Diagnostics"]["ConductorConsistency"];
      CHECK(diagnostics["Count"].get<int>() >= 1);
      CHECK(diagnostics["FloorAppliedPatches"].get<int>() == 0);
      for (const auto &entry : diagnostics["Records"])
      {
        if (entry["ProbeKind"] == "SpatialVertex")
        {
          CHECK(entry["Excluded"].get<bool>());
          CHECK_FALSE(entry["FloorApplied"].get<bool>());
        }
      }
    }
    // G1b: the strip at the ground potential (a field normal to the plane only): KEPT.
    {
      std::vector<std::unique_ptr<Mesh>> meshes;
      IoData iodata = Build(library, false, meshes);
      SpaceOperator space(iodata, meshes);
      SurfaceResponseOperator response(iodata, space);
      const auto records = SpatialRecords(
          response.ApplyConductorConsistencyGate(Field(space, true), 0, 5.0, false));
      REQUIRE(records.size() == 1);
      const auto &record = records.front();
      CHECK(record.plane_knots == 9);
      CHECK(record.amplitude > 0.0);  // the anchor line up to the ring reads E_y x 0.1
      CHECK(record.normalization == record.amplitude);
      CHECK(record.max_ratio < 1.0e-8);
      CHECK_FALSE(record.excluded);
    }
  }

  SECTION("G2: the strip as an island (two labels), two references: kept, the state = V")
  {
    const json feature = IdentifyStripEndCluster(*this, identification_library,
                                                 temp.temp_dir / "stripend-g2", true);
    for (const auto &point : ground_points)
    {
      REQUIRE(OnPlane(feature, point));
    }
    for (const auto &point : strip_points)
    {
      REQUIRE(OnPlane(feature, point));
    }
    const auto library = WriteStripEndClusterLibrary(
        temp.temp_dir, convex_library_3d_path, feature, ground_points, strip_points,
        ground_anchor, strip_anchor, true, "g2");
    std::vector<std::unique_ptr<Mesh>> meshes;
    IoData iodata = Build(library, true, meshes);
    SpaceOperator space(iodata, meshes);
    SurfaceResponseOperator response(iodata, space);
    REQUIRE(response.HasConductorConsistencyProbes());
    const auto records = SpatialRecords(
        response.ApplyConductorConsistencyGate(Field(space, false), 0, 5.0, false));
    REQUIRE(records.size() == 1);
    const auto &record = records.front();
    CHECK(record.plane_knots == 9);
    // The lift carries the ground-to-strip voltage as the conductor state; every vertex
    // reads ~0 against its own anchor.
    CHECK_THAT(record.state, WithinRel(V, 0.2));
    CHECK(record.amplitude >= record.state);
    CHECK(record.normalization == record.amplitude);
    CHECK(record.max_ratio < 1.0e-6);
    CHECK_FALSE(record.excluded);
    CHECK(response.GetConductorConsistencyPortions().empty());
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

// OI-COVERAGE-EXCLUSION-ENERGY (idfix DESIGN section 8; decisions 1064 (2) / 1072 Q3 /
// 1074): the config leaf Diagnostics.ExclusionEnergy; the identification's excluded /
// unclaimed portions by class close the perimeter partition (perimeter = assigned +
// excluded + unclaimed-other) and travel with the v21 cache as a required leaf; the
// electrostatic driver writes surface-response-exclusion-energy.csv when the leaf is on and
// nothing else changes (every other table bitwise ON vs OFF).
TEST_CASE("Response-correction Diagnostics.ExclusionEnergy is read and bounded",
          "[surfaceresponseoperator][exclusionenergy][config][Serial]")
{
  json config = {
      {"Problem", {{"Type", "Electrostatic"}, {"Output", "unused"}}},
      {"Model", {{"Mesh", "unused.msh"}}},
      {"Domains", {{"Materials", {{{"Attributes", {1}}}}}}},
      {"Boundaries", {{"Ground", {{"Attributes", {1}}}}}},
      {"Solver",
       {{"Electrostatic",
         {{"ResponseCorrection",
           {{"Library", "unused.json"}, {"Diagnostics", {{"ExclusionEnergy", true}}}}}}}}}};
  {
    IoData iodata(config, false);
    REQUIRE(iodata.solver.electrostatic.response_correction);
    CHECK(iodata.solver.electrostatic.response_correction->exclusion_energy);
  }
  config["Solver"]["Electrostatic"]["ResponseCorrection"]["Diagnostics"] = {
      {"ExclusionEnergy", false}};
  {
    IoData iodata(config, false);
    CHECK_FALSE(iodata.solver.electrostatic.response_correction->exclusion_energy);
  }
  config["Solver"]["Electrostatic"]["ResponseCorrection"].erase("Diagnostics");
  {
    // Off by default (decision 1076: the measured overhead exceeded the 10 % bar; the
    // stage E / F configurations turn it on).
    IoData iodata(config, false);
    CHECK_FALSE(iodata.solver.electrostatic.response_correction->exclusion_energy);
  }
  config["Solver"]["Electrostatic"]["ResponseCorrection"]["Diagnostics"] = {
      {"PatchEnergy", true}};
  CHECK_THROWS_WITH(IoData(config, false),
                    ContainsSubstring("unknown key \"PatchEnergy\""));
}

TEST_CASE_METHOD(test::SurfaceResponseFiles,
                 "SurfaceResponseOperator exclusion-energy portions close the perimeter "
                 "partition and travel with the geometry cache (section 8)",
                 "[surfaceresponseoperator][exclusionenergy][3d][Serial][Parallel]")
{
#if !defined(MFEM_USE_GSLIB)
  SKIP("SurfaceResponseOperator requires MFEM_USE_GSLIB");
#else
  // The electrostatic path on the shorted CPW stub (the one-label cluster library): the
  // ground plane's rim on the box walls is TruncationCut, the short's corners are unmatched
  // (uncovered, not unclaimed).
  const auto library_path =
      WriteStubLibrary(temp.temp_dir, convex_library_3d_path, 0.1, "s8");
  json config = ConvexIslandConfig();
  auto &correction = config["Solver"]["Electrostatic"]["ResponseCorrection"];
  correction.erase("PatchConstruction");
  correction["Library"] = library_path.string();
  correction["UnmatchedPolicy"] = "Warn";
  correction["Diagnostics"] = {{"ExclusionEnergy", true}};
  IoData iodata(config, false);
  iodata.boundaries.cracked_attributes.insert(9);
  const double scale = iodata.units.GetMeshLengthRelativeScale();
  // The preflight manifest's Totals (mesh-file units): the perimeter partition of record.
  const auto manifest_path = temp.temp_dir / "surface-response-requirements-s8.json";
  {
    mfem::Mesh serial = MakeSerialShortedCpwMesh();
    Mesh mesh(std::make_unique<mfem::ParMesh>(Mpi::World(), serial));
    WriteSurfaceResponseRequirements(iodata, mesh, manifest_path.string());
  }
  Mpi::Barrier(Mpi::World());
  std::ifstream manifest_input(manifest_path);
  REQUIRE(manifest_input);
  const json manifest = json::parse(manifest_input);
  const auto &totals = manifest.at("Identification").at("Totals");
  const double perimeter = totals.at("PerimeterLength").get<double>();
  const double assigned = totals.at("AssignedLength").get<double>();
  const double excluded = totals.at("ExcludedLength").get<double>();
  REQUIRE(perimeter > 0.0);
  REQUIRE(excluded > 0.0);  // the ground plane's rim and its untargeted perimeter

  // The operator writes the geometry cache (version 21) as it is built.
  const auto cache_path = temp.temp_dir / "s8-geometry-cache.json";
  std::vector<std::unique_ptr<Mesh>> meshes;
  {
    mfem::Mesh serial = MakeSerialShortedCpwMesh();
    meshes.push_back(
        std::make_unique<Mesh>(std::make_unique<mfem::ParMesh>(Mpi::World(), serial)));
  }
  LaplaceOperator laplace(iodata, meshes);
  std::vector<ResponseCorrectionData::UncoveredPortionData> identification_portions;
  {
    test::GeometryCacheEnvGuard cache_env(cache_path.string(), true);
    SurfaceResponseOperator response(iodata, laplace);
    REQUIRE(response.ExclusionEnergyRequested());
    identification_portions = response.GetIdentificationExclusionPortions();
    // The export set = the identification's portions (+ UnprobedModels /
    // MaxwellPatchExcluded, none here: every stub patch is probed), each named by class.
    const auto export_set = response.GetExclusionEnergyPortions();
    CHECK(export_set.size() == identification_portions.size());
    CHECK(response.GetUnprobedModelPortionRecord().portions.empty());
    for (const auto &portion : export_set)
    {
      CHECK_FALSE(portion.topology.empty());
    }
  }
  Mpi::Barrier(Mpi::World());
  auto Length = [](const auto &portion)
  {
    double length2 = 0.0;
    for (int d = 0; d < 3; d++)
    {
      length2 += (portion.p1[d] - portion.p0[d]) * (portion.p1[d] - portion.p0[d]);
    }
    return std::sqrt(length2);
  };
  std::map<std::string, double> length_by_class;
  for (const auto &portion : identification_portions)
  {
    length_by_class[portion.topology] += Length(portion) * scale;
  }
  // Every exclusion class of the manifest, by NAME and length (the portions ARE the
  // manifest's exclusions); no portion of a class the manifest lacks (the identification's
  // classes, CornerArmTrim and Unclaimed-other aside).
  std::map<std::string, double> manifest_by_class;
  for (const auto &exclusion : manifest.at("Identification").at("Exclusions"))
  {
    if (exclusion.at("Length").get<double>() > 0.0)
    {
      manifest_by_class[exclusion.at("Class").get<std::string>()] +=
          exclusion.at("Length").get<double>();
    }
  }
  REQUIRE_FALSE(manifest_by_class.empty());
  double excluded_classes = 0.0;
  for (const auto &[cls, length] : length_by_class)
  {
    if (cls == "Unclaimed-other" || cls == "CornerArmTrim")
    {
      continue;
    }
    excluded_classes += length;
    REQUIRE(manifest_by_class.count(cls) == 1);
    CHECK_THAT(length, WithinRel(manifest_by_class.at(cls), 1.0e-9));
  }
  CHECK(length_by_class.size() - length_by_class.count("Unclaimed-other") -
            length_by_class.count("CornerArmTrim") ==
        manifest_by_class.size());
  // The partition closes: the excluded classes = the manifest's ExcludedLength; the
  // unclaimed remainder = perimeter - assigned - excluded (no Unclaimed-other portion when
  // every non-excluded segment is claimed).
  CHECK_THAT(excluded_classes, WithinRel(excluded, 1.0e-9));
  const double unclaimed = length_by_class.count("Unclaimed-other")
                               ? length_by_class.at("Unclaimed-other")
                               : 0.0;
  CHECK_THAT(unclaimed, WithinAbs(perimeter - assigned - excluded, 1.0e-9 * perimeter));
  // The cache (version 21) carries the portions as a REQUIRED leaf: a round trip returns
  // them; a cache without the leaf is refused by name (decision 1072 Q3).
  REQUIRE(fs::exists(cache_path));
  const auto stale_path = temp.temp_dir / "s8-geometry-cache-stale.json";
  if (Mpi::Root(Mpi::World()))
  {
    std::ifstream input(cache_path);
    json cache = json::parse(input);
    REQUIRE(cache.at("Version").get<int>() == 21);
    REQUIRE(cache.at("ExclusionPortions").size() == identification_portions.size());
    cache.erase("ExclusionPortions");
    std::ofstream output(stale_path);
    output << cache.dump() << "\n";
  }
  Mpi::Barrier(Mpi::World());
  const auto &request = *iodata.solver.electrostatic.response_correction;
  const auto reread = ReadResponseGeometryCache(cache_path, request);
  REQUIRE(reread.exclusion_portions.size() == identification_portions.size());
  for (std::size_t i = 0; i < reread.exclusion_portions.size(); i++)
  {
    CHECK(reread.exclusion_portions[i].topology == identification_portions[i].topology);
    CHECK(reread.exclusion_portions[i].segment == identification_portions[i].segment);
  }
  CHECK_THROWS_WITH(ReadResponseGeometryCache(stale_path, request),
                    ContainsSubstring("without ExclusionPortions"));
#endif
}

TEST_CASE_METHOD(test::SurfaceResponseFiles,
                 "Electrostatic Diagnostics.ExclusionEnergy writes the exclusion-energy "
                 "table and changes nothing else (section 8)",
                 "[electrostaticsolver][surfaceresponseoperator][exclusionenergy][3d]"
                 "[Serial][Parallel]")
{
#if !defined(MFEM_USE_GSLIB)
  SKIP("SurfaceResponseOperator requires MFEM_USE_GSLIB");
#else
  const fs::path mesh_path = temp.temp_dir / "shorted-cpw-stub.mesh";
  if (Mpi::Root(Mpi::World()))
  {
    mfem::Mesh serial = MakeSerialShortedCpwMesh();
    std::ofstream output(mesh_path);
    serial.Print(output);
  }
  Mpi::Barrier(Mpi::World());
  const auto library_path =
      WriteStubLibrary(temp.temp_dir, convex_library_3d_path, 0.1, "es8");
  json config = ConvexIslandConfig();
  config["Model"]["Mesh"] = mesh_path.string();
  config["Model"]["L0"] = 1.0;
  config["Model"]["Lc"] = 1.0;
  // The plane (Terminal 9) inside the grounded box: its rim on the walls is NonManifold
  // (physical perimeter shared with the wall metal: in the interface's edge-distance tree,
  // so its within-R energy is attributable). A TruncationCut / Port cut is NOT physical
  // perimeter and lies in no edge tree: the partition attributes the field near a cut to
  // the nearest physical edge, so those classes read 0 here (the named limitation of the
  // export; the lane record carries it).
  auto &correction = config["Solver"]["Electrostatic"]["ResponseCorrection"];
  correction.erase("PatchConstruction");
  correction["Library"] = library_path.string();
  correction["UnmatchedPolicy"] = "Warn";
  config["Solver"]["Linear"] = {{"Tol", 1.0e-12}, {"MaxIts", 400}};
  const fs::path off_dir = temp.temp_dir / "exclusion-energy-off";
  const fs::path on_dir = temp.temp_dir / "exclusion-energy-on";
  config["Problem"]["Output"] = off_dir.string();  // the default: off (decision 1076)
  test::RunElectrostatic(config);
  correction["Diagnostics"] = {{"ExclusionEnergy", true}};
  config["Problem"]["Output"] = on_dir.string();
  test::RunElectrostatic(config);
  Mpi::Barrier(Mpi::World());
  if (!Mpi::Root(Mpi::World()))
  {
    return;
  }
  const fs::path table = "surface-response-exclusion-energy.csv";
  CHECK_FALSE(fs::exists(off_dir / table));
  REQUIRE(fs::exists(on_dir / table));
  // Every other CSV of the run is bitwise ON vs OFF (the export is a diagnostic beside).
  int compared = 0;
  for (const auto &entry : fs::directory_iterator(off_dir))
  {
    if (entry.path().extension() != ".csv")
    {
      continue;
    }
    const auto name = entry.path().filename();
    REQUIRE(fs::exists(on_dir / name));
    std::ifstream off(entry.path(), std::ios::binary), on(on_dir / name, std::ios::binary);
    const std::string off_bytes((std::istreambuf_iterator<char>(off)), {});
    const std::string on_bytes((std::istreambuf_iterator<char>(on)), {});
    CHECK(off_bytes == on_bytes);
    compared++;
  }
  CHECK(compared > 0);
  // The table: the header of record, one row per class x type for the single source; the
  // NonManifold class (the plane's rim on the grounded walls) reads > 0 against the
  // classified tree (decision 1076); every energy non-negative; P == S: the classes and
  // Claimed sum to the Total row (the within-R energy under the classified tree by one
  // integral) to 1e-12 relative, and no more than the type's interface raw energy.
  std::ifstream input(on_dir / table);
  std::string line;
  REQUIRE(std::getline(input, line));
  CHECK(line ==
        "source,class,type,class portions,class length (m),raw energy (J),interface "
        "raw energy (J),share");
  std::map<std::string, double> energy_by_class_type, total_by_type, interface_raw_by_type;
  std::map<std::string, double> energy_by_class;
  std::set<std::string> classes;
  int rows = 0;
  while (std::getline(input, line))
  {
    std::vector<std::string> fields;
    std::stringstream row(line);
    std::string field;
    while (std::getline(row, field, ','))
    {
      fields.push_back(field);
    }
    REQUIRE(fields.size() == 8);
    CHECK(fields[0] == "1");  // the one source
    const std::string &cls = fields[1], &type = fields[2];
    const double energy = std::stod(fields[5]);
    CHECK(energy >= 0.0);
    CHECK(std::stod(fields[7]) >= 0.0);
    interface_raw_by_type[type] = std::stod(fields[6]);
    energy_by_class[cls] += energy;
    if (cls == "Total")
    {
      total_by_type[type] = energy;
    }
    else
    {
      classes.insert(cls);
      energy_by_class_type[type] += energy;
      if (cls != "Claimed")
      {
        CHECK(std::stoi(fields[3]) > 0);
        CHECK(std::stod(fields[4]) > 0.0);
      }
    }
    rows++;
  }
  CHECK(rows > 0);
  CHECK(classes.count("NonManifold") == 1);
  CHECK(classes.count("Claimed") == 1);
  // The NonManifold rim's row carries energy (the classified tree attributes it there).
  CHECK(energy_by_class["NonManifold"] > 0.0);
  REQUIRE_FALSE(total_by_type.empty());
  for (const auto &[type, total] : total_by_type)
  {
    CHECK(total > 0.0);
    CHECK_THAT(energy_by_class_type[type], WithinRel(total, 1.0e-12));
    CHECK(total <= interface_raw_by_type[type] * (1.0 + 1.0e-9));
  }
#endif
}

}  // namespace palace
