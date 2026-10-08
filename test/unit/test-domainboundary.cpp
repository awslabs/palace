// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "fixtures.hpp"
#include "surfaceresponse-fixtures.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <map>
#include <memory>
#include <set>
#include <sstream>
#include <string>
#include <tuple>
#include <vector>
#include <fmt/format.h>
#include <mfem.hpp>
#include <catch2/catch_test_macros.hpp>
#include <nlohmann/json.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>
#include "fem/fespace.hpp"
#include "fem/gridfunction.hpp"
#include "fem/mesh.hpp"
#include "models/laplaceoperator.hpp"
#include "models/materialoperator.hpp"
#include "models/surfacepostoperator.hpp"
#include "models/surfaceresponseidentification.hpp"
#include "models/surfaceresponsemirror.hpp"
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

// A PEC lead (attribute 9) on the plane y = 0.5 of the unit cube meshed 8 x 8 x 8 (h =
// 0.125): x in [0.25, 0.75] from the bottom face z = 0 (where the lead is cut by the
// domain) to its end cap at z = 0.5 (two 90-degree convex corners). The mesh is then
// sheared out of the metal plane, z -> z + shear (y - 0.5), so that the metal stays
// planar and rectangular (the identification sees the unsheared lead) while the domain
// faces z = 0 and z = 1 tilt about x: the lead's long edges (along z) meet the cut
// obliquely, at the angle theta = atan(shear) to its normal (0, -shear, 1) / |.|, and a
// coupon cross-section at normal distance d from the cut sticks out of it where its
// vertical half-extent h_v exceeds d / sin(theta). With notch, the element (x in [0.125,
// 0.25], y in [0.25, 0.375], z in [0.25, 0.375]) is removed: a cavity between the knots
// of the left edge's cross-sections (x = 0.05 and 0.45 at y = 0.3) that its mortar
// samples along the ring segment between them enter.
mfem::Mesh MakeSerialObliqueCutLeadMesh(double shear, bool notch)
{
  constexpr int n = 8;
  constexpr double h = 1.0 / n;
  auto Vertex = [](int i, int j, int k) { return i + (n + 1) * (j + (n + 1) * k); };
  auto IsNotch = [&](int i, int j, int k) { return notch && i == 1 && j == 2 && k == 2; };
  int element_count = 0;
  for (int k = 0; k < n; k++)
  {
    for (int j = 0; j < n; j++)
    {
      for (int i = 0; i < n; i++)
      {
        element_count += IsNotch(i, j, k) ? 0 : 1;
      }
    }
  }
  mfem::Mesh serial(3, (n + 1) * (n + 1) * (n + 1), element_count, 0, 3);
  for (int k = 0; k <= n; k++)
  {
    for (int j = 0; j <= n; j++)
    {
      for (int i = 0; i <= n; i++)
      {
        serial.AddVertex(i * h, j * h, k * h);
      }
    }
  }
  for (int k = 0; k < n; k++)
  {
    for (int j = 0; j < n; j++)
    {
      for (int i = 0; i < n; i++)
      {
        if (IsNotch(i, j, k))
        {
          continue;
        }
        const int vertices[8] = {Vertex(i, j, k),
                                 Vertex(i + 1, j, k),
                                 Vertex(i + 1, j + 1, k),
                                 Vertex(i, j + 1, k),
                                 Vertex(i, j, k + 1),
                                 Vertex(i + 1, j, k + 1),
                                 Vertex(i + 1, j + 1, k + 1),
                                 Vertex(i, j + 1, k + 1)};
        serial.AddHex(vertices, 1);
      }
    }
  }
  serial.FinalizeTopology();  // generates the exterior boundary elements (attribute 1)
  // The Cartesian attributes (1 bottom z, 2 front y, 3 right x, 4 back y, 5 left x, 6 top
  // z); the notch's cavity faces keep attribute 1 (a grounded cavity).
  for (int be = 0; be < serial.GetNBE(); be++)
  {
    mfem::Array<int> vertices;
    serial.GetBdrElementVertices(be, vertices);
    std::array<double, 3> centroid{};
    for (const int vertex : vertices)
    {
      for (int d = 0; d < 3; d++)
      {
        centroid[d] += serial.GetVertex(vertex)[d] / vertices.Size();
      }
    }
    int attribute = 1;
    constexpr double tolerance = 1.0e-12;
    if (std::abs(centroid[2]) < tolerance)
    {
      attribute = 1;
    }
    else if (std::abs(centroid[1]) < tolerance)
    {
      attribute = 2;
    }
    else if (std::abs(centroid[0] - 1.0) < tolerance)
    {
      attribute = 3;
    }
    else if (std::abs(centroid[1] - 1.0) < tolerance)
    {
      attribute = 4;
    }
    else if (std::abs(centroid[0]) < tolerance)
    {
      attribute = 5;
    }
    else if (std::abs(centroid[2] - 1.0) < tolerance)
    {
      attribute = 6;
    }
    serial.SetBdrAttribute(be, attribute);
  }
  // The lead: the interior faces on y = 0.5 with x in [0.25, 0.75], z in [0, 0.5].
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
    double xmin = 1.0, xmax = 0.0, zmin = 1.0, zmax = 0.0;
    for (const int vertex : vertices)
    {
      const double *point = serial.GetVertex(vertex);
      on_plane = on_plane && std::abs(point[1] - 0.5) < 1.0e-12;
      xmin = std::min(xmin, point[0]);
      xmax = std::max(xmax, point[0]);
      zmin = std::min(zmin, point[2]);
      zmax = std::max(zmax, point[2]);
    }
    if (on_plane && xmin >= 0.25 - 1.0e-12 && xmax <= 0.75 + 1.0e-12 &&
        zmax <= 0.5 + 1.0e-12)
    {
      serial.AddBdrElement(serial.GetFace(face)->Duplicate(&serial));
      serial.SetBdrAttribute(serial.GetNBE() - 1, 9);
    }
  }
  serial.FinalizeTopology();
  serial.Finalize();
  for (int vertex = 0; vertex < serial.GetNV(); vertex++)
  {
    double *point = serial.GetVertex(vertex);
    point[2] += shear * (point[1] - 0.5);
  }
  return serial;
}

std::unique_ptr<mfem::ParMesh> MakeObliqueCutLeadMesh(double shear, bool notch)
{
  mfem::Mesh serial = MakeSerialObliqueCutLeadMesh(shear, notch);
  return std::make_unique<mfem::ParMesh>(Mpi::World(), serial);
}

// The M1 symmetry fixture (boundary-cut DESIGN 4.1): a pad on the metal plane y = 0.5 of
// the box [0, 1] x [0, 1] x [0, 2] meshed 8 x 8 x 8 (h = 0.125 in x, y; 0.25 in z) and
// sheared IN THE METAL PLANE, x -> x + s (z - 1), so that the wall x = 1 becomes the
// VERTICAL plane x = 1 + s (z - 1) (a natural cut perpendicular to the metal plane: the
// window configuration). The pad spans x in [0, 1] (from the wall x = 0 to the wall x = 1:
// a strip between two parallel cuts, no real corner), z in [0.75, 1.25] before the shear
// (0.75 > 3 R from the z walls: their images reach nothing): its long edges (along z,
// tilted to (s, 0, 1)) are parallel to the walls, its short edges (along x) meet each wall
// at theta = atan(1 / s) (s = 1: 45 degrees; the virtual corners with the images are 2
// theta = 90 degrees, convex at the top edge (metal in the acute wedge), concave at the
// bottom edge). The FULL mesh is the half plus its exact reflection across the wall
// (reflected vertices appended, the wall's vertices shared, hex vertex order mirrored,
// boundary attributes mapped, the pad faces reflected): the mirror-symmetric full problem
// whose restriction the half solves.
// A pad rectangle in unsheared coordinates: {x0, z0, x1, z1}; the metal is the union of the
// rectangles (an L shape, a slot).
using PadRectangle = std::array<double, 4>;
const std::vector<PadRectangle> kDefaultPads = {{0.0, 0.75, 1.0, 1.25}};

// A BEND of the pad's edges along x (decision 512, the chains lane): the mesh is shifted in
// z by t(x), piecewise linear with the given slope on every x cell (n = 8 entries, the
// unsheared x in [i h, (i + 1) h]); every edge along x becomes a polyline with a joint at
// each slope change (sub-noise joints: one chain with windowed curvature), the x walls stay
// the vertical planes x = 0 / x = 1 and the pad's edges meet them at atan(1 / slope). Only
// with s = 0 (the unshear is then x-independent).
using PadBend = std::vector<double>;

mfem::Mesh MakeSerialMirrorPadMesh(double s, bool full,
                                   const std::vector<PadRectangle> &pads = kDefaultPads,
                                   int nz = 8, const PadBend &bend = {})
{
  constexpr int n = 8;
  constexpr double h = 1.0 / n;
  const double hz = 2.0 / nz;
  constexpr double z_center = 1.0;
  constexpr double tolerance = 1.0e-12;
  REQUIRE((bend.empty() || (bend.size() == n && s == 0.0)));
  auto Bend = [&](double x)
  {
    double t = 0.0;
    for (std::size_t i = 0; i < bend.size(); i++)
    {
      t += bend[i] * std::clamp(x - i * h, 0.0, h);
    }
    return t;
  };
  // The wall x = 1 + s (z - 1): unit outward normal and offset.
  const double norm = std::sqrt(1.0 + s * s);
  const std::array<double, 3> normal = {1.0 / norm, 0.0, -s / norm};
  const double offset = (1.0 - z_center * s) / norm;
  auto Reflect = [&](const std::array<double, 3> &p)
  {
    const double d = normal[0] * p[0] + normal[1] * p[1] + normal[2] * p[2] - offset;
    return std::array<double, 3>{p[0] - 2.0 * d * normal[0], p[1] - 2.0 * d * normal[1],
                                 p[2] - 2.0 * d * normal[2]};
  };
  auto OnWall = [&](const std::array<double, 3> &p)
  {
    return std::abs(normal[0] * p[0] + normal[1] * p[1] + normal[2] * p[2] - offset) <
           tolerance;
  };
  // Half vertices (sheared) and hexes.
  std::vector<std::array<double, 3>> vertices;
  auto Vertex = [](int i, int j, int k) { return i + (n + 1) * (j + (n + 1) * k); };
  for (int k = 0; k <= nz; k++)
  {
    for (int j = 0; j <= n; j++)
    {
      for (int i = 0; i <= n; i++)
      {
        vertices.push_back({i * h + s * (k * hz - z_center), j * h, k * hz + Bend(i * h)});
      }
    }
  }
  std::vector<std::array<int, 8>> hexes;
  for (int k = 0; k < nz; k++)
  {
    for (int j = 0; j < n; j++)
    {
      for (int i = 0; i < n; i++)
      {
        hexes.push_back({Vertex(i, j, k), Vertex(i + 1, j, k), Vertex(i + 1, j + 1, k),
                         Vertex(i, j + 1, k), Vertex(i, j, k + 1), Vertex(i + 1, j, k + 1),
                         Vertex(i + 1, j + 1, k + 1), Vertex(i, j + 1, k + 1)});
      }
    }
  }
  const std::size_t half_vertices = vertices.size(), half_hexes = hexes.size();
  if (full)
  {
    // The reflected copy: the wall's vertices are shared, every other vertex reflected; a
    // reflected hex keeps a positive Jacobian with its first and second vertex pairs
    // swapped (the mirror of the reference cube's x axis).
    std::vector<int> image(half_vertices);
    for (std::size_t v = 0; v < half_vertices; v++)
    {
      if (OnWall(vertices[v]))
      {
        image[v] = static_cast<int>(v);
      }
      else
      {
        image[v] = static_cast<int>(vertices.size());
        vertices.push_back(Reflect(vertices[v]));
      }
    }
    for (std::size_t e = 0; e < half_hexes; e++)
    {
      const auto &hex = hexes[e];
      hexes.push_back({image[hex[1]], image[hex[0]], image[hex[3]], image[hex[2]],
                       image[hex[5]], image[hex[4]], image[hex[7]], image[hex[6]]});
    }
  }
  mfem::Mesh serial(3, static_cast<int>(vertices.size()), static_cast<int>(hexes.size()), 0,
                    3);
  for (const auto &vertex : vertices)
  {
    serial.AddVertex(vertex[0], vertex[1], vertex[2]);
  }
  for (const auto &hex : hexes)
  {
    serial.AddHex(hex.data(), 1);
  }
  serial.FinalizeTopology();  // the exterior boundary elements (attribute 1)
  // The Cartesian attributes of the unsheared half (1 bottom z, 2 front y, 3 right x = the
  // wall, 4 back y, 5 left x, 6 top z), the image's faces mapped to the same attributes
  // (the reflection preserves y: the ground y = 1 is one symmetric face); the wall is
  // interior in the full mesh.
  auto Unshear = [&](std::array<double, 3> p)
  {
    // A point of the image half is mapped back by the reflection first.
    const double d = normal[0] * p[0] + normal[1] * p[1] + normal[2] * p[2] - offset;
    if (d > tolerance)
    {
      p = Reflect(p);
    }
    p[0] -= s * (p[2] - z_center);
    p[2] -= Bend(p[0]);
    return p;
  };
  for (int be = 0; be < serial.GetNBE(); be++)
  {
    mfem::Array<int> face_vertices;
    serial.GetBdrElementVertices(be, face_vertices);
    std::array<double, 3> centroid{};
    for (const int vertex : face_vertices)
    {
      for (int d = 0; d < 3; d++)
      {
        centroid[d] += serial.GetVertex(vertex)[d] / face_vertices.Size();
      }
    }
    const auto u = Unshear(centroid);
    int attribute = 1;
    if (std::abs(u[2]) < tolerance)
    {
      attribute = 1;
    }
    else if (std::abs(u[1]) < tolerance)
    {
      attribute = 2;
    }
    else if (std::abs(u[0] - 1.0) < tolerance)
    {
      attribute = 3;
    }
    else if (std::abs(u[1] - 1.0) < tolerance)
    {
      attribute = 4;
    }
    else if (std::abs(u[0]) < tolerance)
    {
      attribute = 5;
    }
    else if (std::abs(u[2] - 2.0) < tolerance)
    {
      attribute = 6;
    }
    serial.SetBdrAttribute(be, attribute);
  }
  // The pad: the interior faces on y = 0.5 with (unsheared) x in [0, 1], z in [0.75,
  // 1.25], on both halves of the full mesh.
  for (int face = 0; face < serial.GetNumFaces(); face++)
  {
    int element1, element2;
    serial.GetFaceElements(face, &element1, &element2);
    if (element1 < 0 || element2 < 0)
    {
      continue;
    }
    mfem::Array<int> face_vertices;
    serial.GetFaceVertices(face, face_vertices);
    bool on_plane = true;
    std::array<double, 3> centroid{};
    for (const int vertex : face_vertices)
    {
      const double *point = serial.GetVertex(vertex);
      on_plane = on_plane && std::abs(point[1] - 0.5) < tolerance;
      for (int d = 0; d < 3; d++)
      {
        centroid[d] += point[d] / face_vertices.Size();
      }
    }
    // A face belongs to the pad when its (unsheared) centroid lies in one of the
    // rectangles.
    const auto u = Unshear(centroid);
    bool inside = false;
    for (const auto &pad : pads)
    {
      inside = inside || (u[0] >= pad[0] - tolerance && u[0] <= pad[2] + tolerance &&
                          u[2] >= pad[1] - tolerance && u[2] <= pad[3] + tolerance);
    }
    if (on_plane && inside)
    {
      serial.AddBdrElement(serial.GetFace(face)->Duplicate(&serial));
      serial.SetBdrAttribute(serial.GetNBE() - 1, 9);
    }
  }
  serial.FinalizeTopology();
  serial.Finalize();
  return serial;
}

// Whether a point (sheared coordinates) lies in the device of MakeObliqueCutLeadMesh.
bool InsideObliqueCutLead(const std::array<double, 3> &point, double shear, bool notch)
{
  const double z = point[2] - shear * (point[1] - 0.5);
  constexpr double tolerance = 1.0e-12;
  if (point[0] < -tolerance || point[0] > 1.0 + tolerance || point[1] < -tolerance ||
      point[1] > 1.0 + tolerance || z < -tolerance || z > 1.0 + tolerance)
  {
    return false;
  }
  if (notch && point[0] > 0.125 + tolerance && point[0] < 0.25 - tolerance &&
      point[1] > 0.25 + tolerance && point[1] < 0.375 - tolerance && z > 0.25 + tolerance &&
      z < 0.375 - tolerance)
  {
    return false;
  }
  return true;
}

struct DryRunRow
{
  std::size_t patch = 0;
  int feature = -1;
  std::string topology;
  double weight = 0.0;
  double s0 = 0.0, s1 = 0.0;
  std::array<double, 3> origin{}, axis_u{}, axis_v{}, axis_w{};
  std::array<double, 2> strip{};
};

std::vector<DryRunRow> ReadDryRun(const fs::path &path)
{
  std::ifstream input(path);
  REQUIRE(input);
  std::string line;
  std::getline(input, line);
  std::vector<std::string> header;
  {
    std::stringstream stream(line);
    std::string field;
    while (std::getline(stream, field, ','))
    {
      header.push_back(field);
    }
  }
  auto Column = [&](const std::string &name)
  {
    const auto it = std::find(header.begin(), header.end(), name);
    REQUIRE(it != header.end());
    return static_cast<std::size_t>(it - header.begin());
  };
  const std::size_t patch = Column("Patch"), feature = Column("Feature"),
                    topology = Column("Topology"), weight = Column("Weight"),
                    s0 = Column("S0"), s1 = Column("S1"), origin = Column("OriginX"),
                    axis_u = Column("AxisUX"), axis_v = Column("AxisVX"),
                    axis_w = Column("AxisWX"), strip = Column("StripBegin");
  std::vector<DryRunRow> rows;
  while (std::getline(input, line))
  {
    std::vector<std::string> fields;
    std::stringstream stream(line);
    std::string field;
    while (std::getline(stream, field, ','))
    {
      fields.push_back(field);
    }
    REQUIRE(fields.size() == header.size());
    DryRunRow row;
    row.patch = std::stoul(fields[patch]);
    row.feature = std::stoi(fields[feature]);
    row.topology = fields[topology];
    row.weight = std::stod(fields[weight]);
    row.s0 = std::stod(fields[s0]);
    row.s1 = std::stod(fields[s1]);
    for (int d = 0; d < 3; d++)
    {
      row.origin[d] = std::stod(fields[origin + d]);
      row.axis_u[d] = std::stod(fields[axis_u + d]);
      row.axis_v[d] = std::stod(fields[axis_v + d]);
      row.axis_w[d] = std::stod(fields[axis_w + d]);
    }
    row.strip = {std::stod(fields[strip]), std::stod(fields[strip + 1])};
    rows.push_back(row);
  }
  return rows;
}

}  // namespace

// Decision 258: a library-placed patch any of whose placed coupon points (the model's basis
// points and conductor references at the origin cross-section and at both longitudinal cell
// ends moved 1e-3 R inward) lies outside the device mesh is not applied and is recorded as
// a DomainBoundary exclusion — in the preflight manifest (Identification.Diagnostics, the
// inventory next to Missing, Weight 0 in the dry run) and by the operator, through one
// containment test whose decision is independent of the rank count ([Parallel]; the
// operator's point location routes every point to the ranks whose boxes contain it at any
// rank count, decision 346 (b)). Here the lead's long edges meet the
// tilted cut (shear 0.45: theta = 24.2 deg, R sin theta = 0.082) obliquely; each edge's
// 0.3 from the cut to the corner window is three mesh-segment portions of two Gauss cells
// (0.0625). The first cell [0, 0.0625] has its origin at 0.0264 along the edge (normal
// distance 0.024 < 0.082) and its begin section on the cut: excluded; the second cell
// [0.0625, 0.125] has its origin inside (0.0986 x cos theta = 0.090 > 0.082) but its begin
// section at normal distance 0.057 < 0.082 sticks out: excluded through the cell end
// alone (a longer cut portion would otherwise pass the preflight and abort on its strip
// slices); the third cell begins at 0.114 > 0.082: applied, as are the remaining cells,
// the cap portion and the corners. A point the operator cannot locate for any other
// reason (the notch: a cavity between the knots of a cross-section that only the mortar's
// ring samples enter) still fails closed, naming the patch.
TEST_CASE_METHOD(test::SurfaceResponseFiles,
                 "SurfaceResponseOperator domain-boundary exclusion",
                 "[surfaceresponseoperator][domainboundary][3d][Serial][Parallel]")
{
#if !defined(MFEM_USE_GSLIB)
  SKIP("SurfaceResponseOperator requires MFEM_USE_GSLIB");
#else
  constexpr double shear = 0.45;
  constexpr double R = 0.2;
  const double inset = kSignatureParameterToleranceOverRadius * R;
  // The isolated-edge coupon with basis points at the full box (+-R laterally and
  // vertically), keyed on SA only like the island's features.
  const auto basis_path = temp.temp_dir / "domain-boundary-basis-points.csv";
  const auto library_path = temp.temp_dir / "fabrication-process-domain-boundary-3d.json";
  if (Mpi::Root(Mpi::World()))
  {
    {
      std::ofstream output(basis_path);
      output << "x,y,z\n"
             << "-0.2,-0.2,0.0\n"
             << "0.2,-0.2,0.0\n"
             << "0.2,0.2,0.0\n"
             << "-0.2,0.2,0.0\n";
    }
    std::ifstream input(convex_library_3d_path);
    REQUIRE(input);
    json library = json::parse(input);
    library["Name"] = "unit-test-process-domain-boundary-3d";
    for (auto &model : library["Models"])
    {
      model["Interfaces"] = {{{"Type", "SA"}, {"Coupon", 1}}};
      if (model["Topology"] == "IsolatedEdge")
      {
        model["BasisPoints"] = basis_path.string();
      }
    }
    std::ofstream output(library_path);
    output << library.dump(2) << "\n";
  }
  Mpi::Barrier(Mpi::World());
  json config = IslandConfig();
  // The cut faces are artificial domain boundaries, not metal: like a window's walls they
  // carry no boundary condition (natural); the ground is the far face y = 1 only.
  config["Boundaries"]["Ground"]["Attributes"] = {4};
  auto &correction = config["Solver"]["Electrostatic"]["ResponseCorrection"];
  correction.erase("PatchConstruction");
  correction["Library"] = library_path.string();
  correction["TraceCoupling"] = "SurfaceMortar";
  correction["MortarOversampling"] = 2;
  IoData iodata(config, false);
  iodata.boundaries.cracked_attributes.insert(9);
  const auto manifest_path = temp.temp_dir / "surface-response-requirements-domain.json";
  const auto patches_path = temp.temp_dir / "surface-response-patches.csv";

  std::vector<std::array<double, 3>> basis_points = {
      {-0.2, -0.2, 0.0}, {0.2, -0.2, 0.0}, {0.2, 0.2, 0.0}, {-0.2, 0.2, 0.0}};
  // The rule's placed points of a dry-run row (the basis points and the reference at the
  // origin and at the inset cell ends), and whether all lie inside the device.
  auto Expected = [&](const DryRunRow &row, bool notch)
  {
    std::vector<double> sections = {0.0};
    if (row.strip[1] > row.strip[0])
    {
      sections.push_back(std::min(0.0, row.strip[0] + inset));
      sections.push_back(std::max(0.0, row.strip[1] - inset));
    }
    int outside = 0;
    for (const double section : sections)
    {
      std::vector<std::array<double, 3>> points;
      if (row.topology == "isolated edge")
      {
        for (const auto &local : basis_points)
        {
          std::array<double, 3> point{};
          for (int d = 0; d < 3; d++)
          {
            point[d] = row.origin[d] + local[0] * row.axis_u[d] + local[1] * row.axis_v[d] +
                       section * row.axis_w[d];
          }
          points.push_back(point);
        }
      }
      std::array<double, 3> reference{};
      for (int d = 0; d < 3; d++)
      {
        reference[d] = row.origin[d] + section * row.axis_w[d];
      }
      points.push_back(reference);
      for (const auto &point : points)
      {
        outside += InsideObliqueCutLead(point, shear, notch) ? 0 : 1;
      }
    }
    return outside;
  };

  auto Preflight = [&](bool notch)
  {
    Mesh dry_run_mesh(MakeObliqueCutLeadMesh(shear, notch));
    WriteSurfaceResponseRequirements(iodata, dry_run_mesh, manifest_path.string());
    Mpi::Barrier(Mpi::World());
    std::ifstream input(manifest_path);
    REQUIRE(input);
    return json::parse(input);
  };

  SECTION("preflight record, dry run and operator agree on 1 and 2 ranks")
  {
    const json manifest = Preflight(false);
    const auto rows = ReadDryRun(patches_path);
    // The lead's features: two long edges (cut end to corner window), the cap portion,
    // two convex corners; the long edges carry two cells each.
    std::set<int> isolated_features;
    int corner_rows = 0;
    for (const auto &row : rows)
    {
      if (row.topology == "isolated edge")
      {
        isolated_features.insert(row.feature);
      }
      else
      {
        REQUIRE(row.topology == "convex corner");
        corner_rows++;
      }
    }
    REQUIRE(corner_rows == 2);
    REQUIRE(isolated_features.size() == 3);
    REQUIRE(rows.size() == 18);

    // The record is exactly the rows whose placed sections leave the device: one per long
    // edge, the cell touching the cut; the rule's R |sin theta| in the sheared geometry.
    std::set<std::size_t> expected;
    for (const auto &row : rows)
    {
      if (Expected(row, false) > 0)
      {
        expected.insert(row.patch);
      }
    }
    REQUIRE(expected.size() == 4);
    const double sin_theta = shear / std::sqrt(1.0 + shear * shear);
    const double cos_theta = 1.0 / std::sqrt(1.0 + shear * shear);
    int long_edge_cells = 0;
    for (const auto &row : rows)
    {
      if (row.topology != "isolated edge" || std::abs(row.axis_w[2]) < 0.5)
      {
        continue;  // a corner, or the cap portion (along x)
      }
      long_edge_cells++;
      // The long edges run along z from the cut at z = 0 (on the metal plane the shear
      // vanishes); the cell's nearest end to the cut, at normal distance (z along the edge)
      // x cos theta, decides: excluded iff that distance is below the stick-out R sin
      // theta of the vertical half-extent R.
      CHECK_THAT(std::abs(row.axis_w[2]), WithinAbs(1.0, 1.0e-12));
      const double cell_begin = std::min(row.origin[2] + row.strip[0] * row.axis_w[2],
                                         row.origin[2] + row.strip[1] * row.axis_w[2]);
      CHECK((expected.count(row.patch) > 0) == (cell_begin * cos_theta < R * sin_theta));
      if (expected.count(row.patch))
      {
        CHECK_THAT(row.strip[1] - row.strip[0], WithinAbs(0.0625, 1.0e-9));
        // The first cell's origin section itself leaves the cut (the S1p pattern: 0.0264
        // cos theta = 0.024 < 0.082); the second cell's origin is inside (0.0986 cos theta
        // = 0.090 > 0.082) and the cell is excluded through its begin section alone.
        CHECK((row.origin[2] * cos_theta < R * sin_theta) == (cell_begin < 1.0e-9));
      }
    }
    REQUIRE(long_edge_cells == 12);
    const auto &diagnostics =
        manifest["Identification"]["Diagnostics"]["DomainBoundaryExclusions"];
    REQUIRE(diagnostics["Count"].get<int>() == 4);
    REQUIRE(diagnostics["Features"].get<int>() == 2);
    REQUIRE(diagnostics["TestedPatches"].get<int>() == 18);
    REQUIRE(diagnostics["WallTime"].get<double>() >= 0.0);
    std::set<std::size_t> recorded;
    for (const auto &entry : diagnostics["Patches"])
    {
      recorded.insert(entry["Patch"].get<std::size_t>());
      CHECK_THAT(entry["CellLength"].get<double>(), WithinAbs(0.0625, 1.0e-9));
      CHECK_THAT(entry["PortionLength"].get<double>(), WithinAbs(0.125, 1.0e-9));
      CHECK(entry["OutsidePoints"].get<int>() > 0);
      CHECK(entry["OutsidePoints"].get<int>() ==
            Expected(rows[entry["Patch"].get<std::size_t>()], false));
      CHECK(entry["TestedPoints"].get<int>() == 15);  // 3 sections x (4 knots + reference)
      CHECK(entry["Topology"] == "isolated edge");
      // The nearest outside point lies beyond the cut by less than the stick-out R sin
      // theta; its recorded distance (to the nearest element bounding box: a lower bound
      // under a tilted face) is below that too.
      const auto &point = entry["NearestOutsidePoint"];
      CHECK_FALSE(InsideObliqueCutLead(
          {point[0].get<double>(), point[1].get<double>(), point[2].get<double>()}, shear,
          false));
      CHECK(entry["NearestDistance"].get<double>() >= 0.0);
      CHECK(entry["NearestDistance"].get<double>() < R * sin_theta);
    }
    CHECK(recorded == expected);
    CHECK_THAT(diagnostics["CellLength"].get<double>(), WithinAbs(0.25, 1.0e-9));
    CHECK_THAT(diagnostics["PortionLength"].get<double>(), WithinAbs(0.5, 1.0e-9));
    CHECK_THAT(diagnostics["CellEndInsetOverRadius"].get<double>(),
               WithinAbs(kSignatureParameterToleranceOverRadius, 0.0));
    // The inventory next to Missing: the cell length, the portion sum alongside.
    const auto &summary = manifest["Summary"];
    CHECK(summary["Counts"]["DomainBoundary"].get<int>() == 4);
    CHECK_THAT(summary["TotalEdgeLengths"]["DomainBoundary"].get<double>(),
               WithinAbs(0.25, 1.0e-9));
    // The matched inventory is untouched (not a library gap): the long edges 0.3 + 0.3,
    // the cap 0.1, the two corners' arms 2 x 2R.
    CHECK_THAT(summary["TotalEdgeLengths"]["Exact"].get<double>(),
               WithinAbs(0.3 + 0.3 + 0.1 + 0.8, 1.0e-9));
    CHECK_THAT(summary["DomainBoundary"]["PortionLength"].get<double>(),
               WithinAbs(0.5, 1.0e-9));
    CHECK(summary["DomainBoundary"]["Patches"].get<int>() == 4);
    // The dry run: the excluded rows carry weight 0 (no new column), every other row a
    // positive weight.
    for (const auto &row : rows)
    {
      if (expected.count(row.patch))
      {
        CHECK(row.weight == 0.0);
      }
      else
      {
        CHECK(row.weight > 0.0);
      }
    }

    // The operator: the same decision (14 applied patches, the same record).
    {
      std::vector<std::unique_ptr<Mesh>> meshes;
      meshes.push_back(std::make_unique<Mesh>(MakeObliqueCutLeadMesh(shear, false)));
      LaplaceOperator laplace(iodata, meshes);
      SurfaceResponseOperator response(iodata, laplace);
      REQUIRE(response.GetPatchCount() == 14);
      const auto statistics = response.GetStatistics();
      const auto &operator_record = statistics["Diagnostics"]["DomainBoundaryExclusions"];
      REQUIRE(operator_record["Count"].get<int>() == 4);
      std::set<std::size_t> operator_patches;
      for (const auto &entry : operator_record["Patches"])
      {
        operator_patches.insert(entry["Patch"].get<std::size_t>());
      }
      CHECK(operator_patches == expected);
      CHECK_THAT(operator_record["CellLength"].get<double>(), WithinAbs(0.25, 1.0e-9));
      CHECK(operator_record["TestedPatches"].get<int>() == 18);
      for (const auto &model : statistics["ModelCatalog"])
      {
        if (model["Name"] == "isolated")
        {
          CHECK(model["PatchCount"].get<int>() == 12);
        }
      }
    }
  }

  SECTION("a point not located for any other reason still fails closed, naming the patch")
  {
    // The notch is not seen by the placed sections (their knots are not in the cavity), so
    // the record is the same four cells; the mortar's ring samples of the left edge's
    // applied cell cross the cavity and the construction aborts.
    const json manifest = Preflight(true);
    const auto rows = ReadDryRun(patches_path);
    std::set<std::size_t> expected;
    for (const auto &row : rows)
    {
      if (Expected(row, true) > 0)
      {
        expected.insert(row.patch);
      }
    }
    REQUIRE(expected.size() == 4);
    const auto &diagnostics =
        manifest["Identification"]["Diagnostics"]["DomainBoundaryExclusions"];
    REQUIRE(diagnostics["Count"].get<int>() == 4);
    std::vector<std::unique_ptr<Mesh>> meshes;
    meshes.push_back(std::make_unique<Mesh>(MakeObliqueCutLeadMesh(shear, true)));
    LaplaceOperator laplace(iodata, meshes);
    CHECK_THROWS_WITH(SurfaceResponseOperator(iodata, laplace),
                      ContainsSubstring("could not be located") &&
                          ContainsSubstring("patch") && ContainsSubstring("isolated"));
  }

  SECTION("a coupon placed off the mesh fails closed as misplaced, naming the patch")
  {
    // Decision 260 (MAJOR-1 of the block-D review): the exclusion is for a cut THROUGH a
    // placed coupon; a patch whose metal-edge reference (the first conductor reference at
    // the origin section) is not located, or none of whose tested points is, is a
    // misplaced or mis-scaled coupon and aborts, and the exclusion may never leave no
    // applied patch. Patches in the lead's frame (u from the metal into the gap, v the
    // plane normal, w along the edge) on the sheared mesh, tested directly through the
    // collective containment test (identical on 1 and 2 ranks).
    using ResponsePatchData = config::ElectrostaticSolverData::ResponseCorrectionPatchData;
    auto mesh = MakeObliqueCutLeadMesh(shear, false);
    const std::vector<std::array<double, 3>> wide_basis_points = {
        {-2.0, -2.0, 0.0}, {2.0, -2.0, 0.0}, {2.0, 2.0, 0.0}, {-2.0, 2.0, 0.0}};
    auto BasisPoints = [&](int model_idx) -> const std::vector<std::array<double, 3>> *
    { return model_idx == 1 ? &wide_basis_points : &basis_points; };
    auto Spatial = [](int) { return false; };
    auto Name = [](int model_idx) { return std::string(model_idx == 1 ? "wide" : "edge"); };
    auto Patch = [&](const std::array<double, 3> &origin, int model)
    {
      ResponsePatchData patch;
      patch.model = model;
      patch.origin = origin;
      patch.axis_u = {-1.0, 0.0, 0.0};
      patch.axis_v = {0.0, 1.0, 0.0};
      patch.axis_w = {0.0, 0.0, 1.0};
      patch.longitudinal_cell = {-0.03, 0.03};
      return patch;
    };
    auto Find = [&](std::vector<ResponsePatchData> &patches)
    {
      return FindDomainBoundaryExclusions(*mesh, patches, BasisPoints, Spatial, Name, 1.0,
                                          R, {});
    };
    // A well-placed cell on the left long edge (x = 0.25, y = 0.5) far from the cut, and
    // one at z = 0.02 whose section sticks out of the tilted cut (the S1p pattern) while
    // its reference is inside: one exclusion, one applied patch, no abort.
    {
      std::vector<ResponsePatchData> patches = {Patch({0.25, 0.5, 0.3}, 0),
                                                Patch({0.25, 0.5, 0.02}, 0)};
      const auto exclusions = Find(patches);
      REQUIRE(exclusions.tested_patches == 2);
      REQUIRE(exclusions.patches.size() == 1);
      CHECK(exclusions.patches[0].patch == 1);
      CHECK(exclusions.patches[0].outside_points > 0);
      CHECK(exclusions.patches[0].outside_points < exclusions.patches[0].tested_points);
      CHECK(patches[0].weight == 1.0);
      CHECK(patches[1].weight == 0.0);
    }
    // The metal-edge reference off the mesh: a misplaced coupon, named 0-based with its
    // model and point, whatever the other patches.
    {
      std::vector<ResponsePatchData> patches = {Patch({0.25, 0.5, 0.3}, 0),
                                                Patch({2.0, 0.5, 0.3}, 0)};
      CHECK_THROWS_WITH(Find(patches),
                        ContainsSubstring("misplaced or mis-scaled coupon") &&
                            ContainsSubstring("patch 1 (0-based") &&
                            ContainsSubstring("model edge") &&
                            ContainsSubstring("metal-edge reference") &&
                            ContainsSubstring("2.000000000e+00"));
    }
    // Without a conductor reference, every tested point outside aborts the same way.
    {
      std::vector<ResponsePatchData> patches = {Patch({0.25, 0.5, 0.3}, 0),
                                                Patch({2.0, 0.5, 0.3}, 0)};
      patches[1].conductor_references.clear();
      CHECK_THROWS_WITH(Find(patches),
                        ContainsSubstring("misplaced or mis-scaled coupon") &&
                            ContainsSubstring("patch 1 (0-based") &&
                            ContainsSubstring("every one of its tested points"));
    }
    // A mis-scaled coupon (basis points 2.0 beyond a 1.0 domain) with its reference inside
    // is excluded; alone, it leaves no applied patch and the exclusion fails closed.
    {
      std::vector<ResponsePatchData> patches = {Patch({0.5, 0.5, 0.5}, 1)};
      CHECK_THROWS_WITH(Find(patches),
                        ContainsSubstring("leaves no applied surface-response patch"));
    }
  }
#endif
}

// F-DB-a (decisions 442 / 454; DESIGN 2.1): a DomainBoundary-excluded patch is NEVER
// dropped. Its raw claim — here the own-edge interval of each excluded cell, the two
// cut-touching cells [0, 0.0625] and [0.0625, 0.125] of both long edges of the sheared lead
// — keeps the device's raw within-R energy in the corrected interface energies (fixed trace
// / fixed flux on the raw field, self-consistent on the corrected field), reported per type
// in surface-response-domain-boundary-energy.csv as the DomainBoundary share. On the
// out-of-plane-sheared cut (a non-vertical truncation face: the Unsupported class of the
// mirror rule, so F-DB-a is the whole correction) with the isolated-edge library alone
// (corners Missing -> uncovered) the within-R raw energy is partitioned EXACTLY by the
// nearest perimeter foot into the DomainBoundary claims, the uncovered corner windows and
// the applied cells (decision 366 Region entries on the same quadrature: the long edges'
// [0.125, 0.3] and the cap's [0.45, 0.55]): ft = E_out + models + uncovered +
// domainboundary to 1e-10, DB + uncovered = within-R - applied cells to 1e-9, the
// DomainBoundary and uncovered portions pairwise disjoint; with the full library (corners
// matched, nothing uncovered) ft = E_out + models + domainboundary and the raw outputs are
// byte-identical.
TEST_CASE_METHOD(test::SurfaceResponseFiles,
                 "Electrostatic DomainBoundary claims keep their raw energy (F-DB-a)",
                 "[electrostaticsolver][surfaceresponseoperator][domainboundary][3d]"
                 "[Serial][Parallel]")
{
#if !defined(MFEM_USE_GSLIB)
  SKIP("SurfaceResponseOperator requires MFEM_USE_GSLIB");
#else
  constexpr double shear = 0.45;
  constexpr double R = 0.2;
  const fs::path mesh_path = temp.temp_dir / "domain-boundary-lead.mesh";
  const auto basis_path = temp.temp_dir / "domain-boundary-basis-points.csv";
  const auto full_library_path = temp.temp_dir / "fabrication-process-db-full.json";
  const auto edge_library_path = temp.temp_dir / "fabrication-process-db-edges.json";
  if (Mpi::Root(Mpi::World()))
  {
    {
      mfem::Mesh serial = MakeSerialObliqueCutLeadMesh(shear, false);
      std::ofstream output(mesh_path);
      serial.Print(output);
    }
    {
      std::ofstream output(basis_path);
      output << "x,y,z\n"
             << "-0.2,-0.2,0.0\n"
             << "0.2,-0.2,0.0\n"
             << "0.2,0.2,0.0\n"
             << "-0.2,0.2,0.0\n";
    }
    std::ifstream input(convex_library_3d_path);
    REQUIRE(input);
    json library = json::parse(input);
    library["Name"] = "unit-test-process-db-full";
    for (auto &model : library["Models"])
    {
      model["Interfaces"] = {{{"Type", "SA"}, {"Coupon", 1}}};
      if (model["Topology"] == "IsolatedEdge")
      {
        model["BasisPoints"] = basis_path.string();
      }
    }
    {
      std::ofstream output(full_library_path);
      output << library.dump(2) << "\n";
    }
    json edges = library;
    edges["Name"] = "unit-test-process-db-edges";
    edges["Models"] = json::array();
    for (const auto &model : library["Models"])
    {
      if (model["Topology"] == "IsolatedEdge")
      {
        edges["Models"].push_back(model);
      }
    }
    REQUIRE(edges["Models"].size() == 1);
    std::ofstream output(edge_library_path);
    output << edges.dump(2) << "\n";
  }
  Mpi::Barrier(Mpi::World());

  json config = IslandConfig();
  config["Model"]["Mesh"] = mesh_path.string();
  config["Model"]["L0"] = 1.0;
  config["Boundaries"]["Ground"]["Attributes"] = {4};
  // The applied cells as decision-366 Region entries (half-open along, transverse <= R):
  // the long edges beyond the two excluded cells, [0.125, 0.3] (the corner window begins
  // at 0.5 - R), and the cap between the corner windows.
  auto &dielectric = config["Boundaries"]["Postprocessing"]["Dielectric"];
  const json target = dielectric[0];
  int index = 5;
  for (const std::array<double, 6> &segment :
       {std::array<double, 6>{0.25, 0.5, 0.125, 0.25, 0.5, 0.5 - R},
        std::array<double, 6>{0.25 + R, 0.5, 0.5, 0.75 - R, 0.5, 0.5},
        std::array<double, 6>{0.75, 0.5, 0.125, 0.75, 0.5, 0.5 - R}})
  {
    json cell = target;
    cell["Index"] = index++;
    cell.erase("AutomaticEdges");
    cell.erase("EdgeDistances");
    cell.erase("EdgeFrameNormal");
    cell["Region"] = {
        {"Segments", {segment}}, {"Distance", R}, {"Normal", {0.0, 1.0, 0.0}}};
    dielectric.push_back(cell);
  }
  auto &correction = config["Solver"]["Electrostatic"]["ResponseCorrection"];
  correction.erase("PatchConstruction");
  correction["TraceCoupling"] = "SurfaceMortar";
  correction["MortarOversampling"] = 2;
  correction["CorrectionMode"] = "Both";
  config["Solver"]["Linear"] = {{"Tol", 1.0e-12}, {"MaxIts", 400}};
  const fs::path full_dir = temp.temp_dir / "db-full";
  const fs::path edges_dir = temp.temp_dir / "db-edges";
  config["Problem"]["Output"] = full_dir.string();
  correction["Library"] = full_library_path.string();
  correction["UnmatchedPolicy"] = "Error";
  test::RunElectrostatic(config);
  config["Problem"]["Output"] = edges_dir.string();
  correction["Library"] = edge_library_path.string();
  correction["UnmatchedPolicy"] = "Warn";
  test::RunElectrostatic(config);
  Mpi::Barrier(Mpi::World());
  if (!Mpi::Root(Mpi::World()))
  {
    return;
  }

  auto ModelEnergySum = [&](const fs::path &dir, int evaluation)
  {
    const Table models = test::LoadCsv(dir / "surface-response-model-energy.csv");
    const Column &eval = test::ColumnByHeader(models, "evaluation");
    const Column &energy = test::ColumnByHeader(models, "fabricated surface energy[4] (J)");
    double sum = 0.0;
    for (std::size_t i = 0; i < eval.data.size(); i++)
    {
      if (static_cast<int>(std::lround(eval.data[i])) == evaluation)
      {
        sum += energy.data[i];
      }
    }
    return sum;
  };
  auto OutsideEnergy = [&](const fs::path &dir)
  {
    const Table edge = test::LoadCsv(dir / "surface-Q-edge.csv");
    const Column &outside = test::ColumnByHeader(edge, "E_out (J)");
    REQUIRE(outside.data.size() == 1);
    return outside.data[0];
  };
  auto RawEnergy = [&](const fs::path &dir, int interface_index)
  {
    const Table surface = test::LoadCsv(dir / "surface-Q.csv");
    const Table domain = test::LoadCsv(dir / "domain-E.csv");
    const Column &participation =
        test::ColumnByHeader(surface, fmt::format("p_surf[{}]", interface_index));
    const Column &energy = test::ColumnByHeader(domain, "E_elec (J)");
    REQUIRE(participation.data.size() == 1);
    REQUIRE(energy.data.size() == 1);
    return participation.data[0] * energy.data[0];
  };
  // The per-type energy table of a run keyed by (evaluation, type).
  struct TypeTable
  {
    std::map<std::pair<int, std::string>, std::vector<std::string>> rows;
    std::size_t portions = 0, length = 0, energy = 0, share = 0, extra = 0;
  };
  auto ReadTypeTable =
      [&](const fs::path &path, const std::string &label, const std::string &extra)
  {
    const test::TextCsv csv = test::ReadTextCsv(path);
    TypeTable table;
    const std::size_t type = csv.Column("type"), evaluation = csv.Column("evaluation");
    table.portions = csv.Column("portions");
    table.length = csv.Column("length (m)");
    table.extra = csv.Column(extra);
    table.energy = csv.Column(label + " raw energy[4] (J)");
    table.share = csv.Column(label + " share[4]");
    for (const auto &row : csv.rows)
    {
      table.rows[{std::stoi(row[evaluation]), row[type]}] = row;
    }
    return table;
  };
  auto CorrectedEnergies = [&](const fs::path &dir)
  {
    const Table corrected = test::LoadCsv(dir / "surface-Q-corrected.csv");
    REQUIRE(corrected.n_rows() == 1);
    return std::array<double, 3>{
        test::ColumnByHeader(corrected, "E_surf postprocessed fixed-trace[4] (J)").data[0],
        test::ColumnByHeader(corrected, "E_surf postprocessed fixed-flux[4] (J)").data[0],
        test::ColumnByHeader(corrected, "E_surf corrected[4] (J)").data[0]};
  };
  // The segments of a record's portion list as (P0, P1) pairs.
  using Segment = std::array<std::array<double, 3>, 2>;
  auto Portions = [](const json &list)
  {
    std::vector<Segment> segments;
    for (const auto &entry : list)
    {
      segments.push_back({entry["P0"].get<std::array<double, 3>>(),
                          entry["P1"].get<std::array<double, 3>>()});
    }
    return segments;
  };
  // The overlap length of two collinear segments (0 when not collinear or disjoint).
  auto Overlap = [](const Segment &a, const Segment &b)
  {
    std::array<double, 3> direction{};
    double length = 0.0;
    for (int d = 0; d < 3; d++)
    {
      direction[d] = a[1][d] - a[0][d];
      length += direction[d] * direction[d];
    }
    length = std::sqrt(length);
    for (int d = 0; d < 3; d++)
    {
      direction[d] /= length;
    }
    auto Along = [&](const std::array<double, 3> &p, double &transverse)
    {
      double along = 0.0;
      std::array<double, 3> r{};
      for (int d = 0; d < 3; d++)
      {
        r[d] = p[d] - a[0][d];
        along += r[d] * direction[d];
      }
      transverse = 0.0;
      for (int d = 0; d < 3; d++)
      {
        transverse += std::pow(r[d] - along * direction[d], 2);
      }
      transverse = std::sqrt(transverse);
      return along;
    };
    double t0 = 0.0, t1 = 0.0;
    const double b0 = Along(b[0], t0), b1 = Along(b[1], t1);
    if (t0 > 1.0e-9 || t1 > 1.0e-9)
    {
      return 0.0;
    }
    return std::max(0.0,
                    std::min(length, std::max(b0, b1)) - std::max(0.0, std::min(b0, b1)));
  };

  // Full library: nothing uncovered, four DomainBoundary cells kept; ft = E_out + models +
  // domainboundary (the identity closes only because the hole is filled).
  REQUIRE(fs::is_regular_file(full_dir / "surface-response-domain-boundary-energy.csv"));
  CHECK_FALSE(fs::exists(full_dir / "surface-response-uncovered-energy.csv"));
  {
    const TypeTable db =
        ReadTypeTable(full_dir / "surface-response-domain-boundary-energy.csv",
                      "domainboundary", "geometric cells");
    // One source x 3 evaluations x (the "isolated edge" model topology + Total).
    REQUIRE(db.rows.size() == 6);
    REQUIRE(db.rows.count({0, "isolated edge"}) == 1);
    const auto &row = db.rows.at({0, "isolated edge"});
    CHECK(std::stoi(row[db.portions]) == 4);
    CHECK(std::stoi(row[db.extra]) == 4);
    CHECK_THAT(std::stod(row[db.length]), WithinAbs(4.0 * 0.0625, 1.0e-9));
    const double db_ft = std::stod(db.rows.at({0, "Total"})[db.energy]);
    CHECK(db_ft > 0.0);
    CHECK(std::stod(row[db.energy]) == db_ft);
    CHECK(std::stod(db.rows.at({1, "Total"})[db.energy]) == db_ft);
    const auto energies = CorrectedEnergies(full_dir);
    const double outside = OutsideEnergy(full_dir);
    CHECK_THAT(energies[0],
               WithinRel(outside + ModelEnergySum(full_dir, 0) + db_ft, 1.0e-10));
    CHECK_THAT(energies[1],
               WithinRel(outside + ModelEnergySum(full_dir, 1) + db_ft, 1.0e-10));
    CHECK_THAT(std::stod(db.rows.at({0, "Total"})[db.share]),
               WithinRel(db_ft / energies[0], 1.0e-9));
    const double db_sc = std::stod(db.rows.at({2, "Total"})[db.energy]);
    REQUIRE(std::isfinite(energies[2]));
    // The corrected field of this tiny fixture is weaker near the cut: the same order,
    // not the same value.
    CHECK(db_sc > 0.0);
    CHECK(db_sc != db_ft);
    CHECK(db_sc < 2.0 * db_ft);
    CHECK_THAT(std::stod(db.rows.at({2, "Total"})[db.share]),
               WithinRel(db_sc / energies[2], 1.0e-9));
    std::ifstream metadata_input(full_dir / "palace.json");
    REQUIRE(metadata_input);
    const auto metadata = json::parse(metadata_input);
    const auto &record =
        metadata.at("SurfaceResponse").at("Diagnostics").at("DomainBoundaryExclusions");
    CHECK(record.at("Count").get<int>() == 4);
    const auto &raw = record.at("RawPortions");
    CHECK(raw.at("Count").get<int>() == 4);
    CHECK(raw.at("GeometricCells").get<int>() == 4);
    CHECK(raw.at("DuplicatePatches").get<int>() == 0);
    CHECK_THAT(raw.at("Length").get<double>(), WithinAbs(0.25, 1.0e-9));
    CHECK(raw.at("ByType").at("isolated edge").at("Portions").get<int>() == 4);
    // Every portion lies on a long edge's line (x = 0.25 or 0.75, y = 0.5) within
    // [0, 0.125] of the cut, each of length 0.0625, pairwise disjoint.
    const auto portions = Portions(raw.at("Portions"));
    for (const auto &portion : portions)
    {
      for (const auto &end : portion)
      {
        CHECK((std::abs(end[0] - 0.25) < 1.0e-9 || std::abs(end[0] - 0.75) < 1.0e-9));
        CHECK_THAT(end[1], WithinAbs(0.5, 1.0e-9));
        CHECK(end[2] >= -1.0e-9);
        CHECK(end[2] <= 0.125 + 1.0e-9);
      }
      CHECK_THAT(std::abs(portion[1][2] - portion[0][2]), WithinAbs(0.0625, 1.0e-9));
    }
    for (std::size_t i = 0; i < portions.size(); i++)
    {
      for (std::size_t j = i + 1; j < portions.size(); j++)
      {
        CHECK(Overlap(portions[i], portions[j]) <= 1.0e-9);
      }
    }
  }

  // Edge-only library: the raw outputs are byte-identical to the full run (the same raw
  // solve); the corner windows are uncovered, the four cells DomainBoundary, and the
  // nearest-foot partition of the within-R raw energy is exact: within - applied cells =
  // uncovered + domainboundary, the two portion sets disjoint.
  for (const char *file : {"terminal-C.csv", "terminal-V.csv", "domain-E.csv",
                           "surface-Q.csv", "surface-Q-edge.csv"})
  {
    INFO(file);
    CHECK(test::ReadFile(full_dir / file) == test::ReadFile(edges_dir / file));
  }
  REQUIRE(fs::is_regular_file(edges_dir / "surface-response-uncovered-energy.csv"));
  REQUIRE(fs::is_regular_file(edges_dir / "surface-response-domain-boundary-energy.csv"));
  {
    const TypeTable uncovered =
        ReadTypeTable(edges_dir / "surface-response-uncovered-energy.csv", "uncovered",
                      "clipped by spatial support (m)");
    const TypeTable db =
        ReadTypeTable(edges_dir / "surface-response-domain-boundary-energy.csv",
                      "domainboundary", "geometric cells");
    const double uncovered_ft =
        std::stod(uncovered.rows.at({0, "Total"})[uncovered.energy]);
    const double db_ft = std::stod(db.rows.at({0, "Total"})[db.energy]);
    CHECK(uncovered_ft > 0.0);
    CHECK(db_ft > 0.0);
    CHECK_THAT(std::stod(db.rows.at({0, "Total"})[db.length]), WithinAbs(0.25, 1.0e-9));
    const auto energies = CorrectedEnergies(edges_dir);
    const double outside = OutsideEnergy(edges_dir);
    CHECK_THAT(
        energies[0],
        WithinRel(outside + ModelEnergySum(edges_dir, 0) + uncovered_ft + db_ft, 1.0e-10));
    CHECK_THAT(
        energies[1],
        WithinRel(outside + ModelEnergySum(edges_dir, 1) + uncovered_ft + db_ft, 1.0e-10));
    const double within = RawEnergy(edges_dir, 4) - outside;
    const double cells =
        RawEnergy(edges_dir, 5) + RawEnergy(edges_dir, 6) + RawEnergy(edges_dir, 7);
    CHECK(within > cells + uncovered_ft);
    CHECK_THAT(uncovered_ft + db_ft, WithinRel(within - cells, 1.0e-9));
    std::ifstream metadata_input(edges_dir / "palace.json");
    REQUIRE(metadata_input);
    const auto metadata = json::parse(metadata_input);
    const auto &diagnostics = metadata.at("SurfaceResponse").at("Diagnostics");
    const auto db_portions = Portions(
        diagnostics.at("DomainBoundaryExclusions").at("RawPortions").at("Portions"));
    const auto uncovered_portions = Portions(diagnostics.at("Uncovered").at("Portions"));
    REQUIRE(db_portions.size() == 4);
    REQUIRE(uncovered_portions.size() == 8);
    for (const auto &a : db_portions)
    {
      for (const auto &b : uncovered_portions)
      {
        CHECK(Overlap(a, b) <= 1.0e-9);
      }
    }
  }
#endif
}

// The SYMMETRY TEST (boundary-cut DESIGN 4.1, M1): for the mirror-symmetric pad solved on
// the FULL domain (no cut through its features) and on the HALF (the wall as a natural
// face), the half's corrected energies under F-DB-c equal HALF the full's: raw, E_out,
// models, uncovered and the per-type terms. S-45 (s = 1): the short edges meet the wall at
// 45 degrees, the virtual corners are 90 degrees (convex at the top, concave at the bottom:
// exact fixture models), s = R so s_half = R and the second-order term is exactly zero (a
// 1e-9-exact case, decision 454 MINOR-4); the pad's real corners at x = 0.5 are 45 / 135
// degrees (Missing in the fixture library: uncovered on both sides). S-90 (s = 0): the
// perpendicular meeting forms no feature (Continued) and the half is bitwise the straight
// continuation. The control with Mirror "Off": every cut-crossing cell is a DomainBoundary
// exclusion whose raw claim is kept (F-DB-a), the ft identity still closes, and the half
// under-reads the full by the raw-vs-model difference of those cells (reported).
TEST_CASE_METHOD(test::SurfaceResponseFiles,
                 "Electrostatic mirror symmetry: half with F-DB-c = full / 2",
                 "[electrostaticsolver][surfaceresponseoperator][domainboundary][mirror]"
                 "[symmetry][3d][Serial][Parallel]")
{
#if !defined(MFEM_USE_GSLIB)
  SKIP("SurfaceResponseOperator requires MFEM_USE_GSLIB");
#else
  constexpr double R = 0.2;
  const auto basis_path = temp.temp_dir / "mirror-pad-basis-points.csv";
  const auto library_path = temp.temp_dir / "fabrication-process-mirror-pad.json";
  const auto stack_library_path =
      temp.temp_dir / "fabrication-process-mirror-pad-stack.json";
  const auto curved_library_path =
      temp.temp_dir / "fabrication-process-mirror-pad-curved.json";
  // A FOURTH library, written by the decision-557 wedge section from the half mesh's own
  // identification: the curved library plus a Signature-keyed SpatialEdgeCluster model per
  // mirror-formed wedge key (the real-half placement by the contract).
  const auto wedge_library_path =
      temp.temp_dir / "fabrication-process-mirror-pad-wedge.json";
  if (Mpi::Root(Mpi::World()))
  {
    {
      std::ofstream output(basis_path);
      output << "x,y,z\n"
             << "-0.2,-0.2,0.0\n"
             << "0.2,-0.2,0.0\n"
             << "0.2,0.2,0.0\n"
             << "-0.2,0.2,0.0\n";
    }
    // The isolated edge, the convex 90 and the concave 90 corner (the concave fixture's
    // corner relabelled: the test is the symmetry identity, not the coupon physics).
    std::ifstream input(convex_library_3d_path);
    REQUIRE(input);
    json library = json::parse(input);
    library["Name"] = "unit-test-process-mirror-pad";
    // The identity test needs symmetric PLACEMENT, not coupon physics: the fixture's
    // isolated / convex-90 models are relabelled into the concave 90, the convex / concave
    // 120 and 152 corners (the S-60 / S-76 virtual corners: s_half > R, the trim path) and
    // the same-conductor gap / strip 0.25 (the parallel cases at d = h = 0.125 < R).
    json convex, isolated;
    for (auto &model : library["Models"])
    {
      model["Interfaces"] = {{{"Type", "SA"}, {"Coupon", 1}}};
      if (model["Topology"] == "IsolatedEdge")
      {
        model["BasisPoints"] = basis_path.string();
        isolated = model;
      }
      if (model["Topology"] == "ConvexCorner")
      {
        convex = model;
      }
    }
    REQUIRE(!convex.is_null());
    REQUIRE(!isolated.is_null());
    for (const auto &[topology, angle] :
         std::vector<std::pair<std::string, double>>{{"ConcaveCorner", 90.0},
                                                     {"ConvexCorner", 120.0},
                                                     {"ConcaveCorner", 120.0},
                                                     {"ConvexCorner", 152.0},
                                                     {"ConcaveCorner", 152.0}})
    {
      json corner = convex;
      corner["Name"] = (topology == "ConvexCorner" ? "convex-corner-" : "concave-corner-") +
                       std::to_string(static_cast<int>(angle));
      corner["Topology"] = topology;
      corner["Angle"] = angle;
      library["Models"].push_back(corner);
    }
    for (const char *topology : {"SameConductorGap", "SameConductorStrip"})
    {
      json pair = isolated;
      pair["Name"] = std::string(topology) + "-0.25";
      pair["Topology"] = topology;
      pair["Separation"] = 0.25;
      pair["SeparationTolerance"] = 1.0e-8;
      pair["Reference"] = {0.0, 0.0, 0.0};
      library["Models"].push_back(pair);
    }
    std::ofstream output(library_path);
    output << library.dump(2) << "\n";
    // A SECOND library for the decision-512 chain sections: the same models plus a
    // curvature family of relabelled isolated models (every node the anchor's matrices AT
    // THE ANCHOR'S COUPON DEPTH, so that every per-length node equals the anchor and any
    // Lagrange blend - whose weights sum to one - IS the anchor, positive semidefinite), so
    // that the bent chains carry the first-order curvature term and their CurvedEdge
    // sections a model - the symmetric PLACEMENT is the test, not the coupon physics. The
    // M1 sections of record keep the first library byte for byte.
    library["Name"] = "unit-test-process-mirror-pad-curved";
    const double anchor_depth =
        isolated.value("CouponDepth", library.at("CouponDepth").get<double>());
    for (const char *convexity : {"Convex", "Concave"})
    {
      for (const double kappa : {0.1, 0.25, 0.5, 0.8})
      {
        json curved = isolated;
        curved["Name"] =
            std::string("curved-edge-") + convexity + "-" + std::to_string(kappa);
        curved["Topology"] = "CurvedEdge";
        curved["Kappa"] = kappa;
        curved["Convexity"] = convexity;
        curved["CouponDepth"] = anchor_depth;
        library["Models"].push_back(curved);
      }
    }
    std::ofstream curved_output(curved_library_path);
    curved_output << library.dump(2) << "\n";
    // A THIRD library for the decision-536 / 537 stack section: the first library plus a
    // relabelled four-edge ParallelEdgeCluster (two strips 0.125 wide separated by a
    // 0.125 gap: edges at lateral offsets 0 / 0.125 / 0.25 / 0.375 from the first side,
    // gap directions alternating, two conductors) on the isolated model's matrices: the
    // symmetric PLACEMENT and the F-DB-a raw-claim mapping are the test, not the coupon
    // physics.
    json stack_library = library;
    stack_library["Name"] = "unit-test-process-mirror-pad-stack";
    {
      json stack = isolated;
      stack["Name"] = "stack-4edge-0.125";
      stack["Topology"] = "ParallelEdgeCluster";
      // 2e-3 (1e-2 R): the bent stack of the decision-545 section tilts its chords by up
      // to 0.13 rad, so the sides' perpendicular offsets read 0.125 cos(theta) = 0.1240.
      stack["EdgeOffsetTolerance"] = 2.0e-3;
      stack["Edges"] = {{{"Offset", 0.0}, {"GapDirection", -1}, {"Conductor", 1}},
                        {{"Offset", 0.125}, {"GapDirection", 1}, {"Conductor", 1}},
                        {{"Offset", 0.25}, {"GapDirection", -1}, {"Conductor", 2}},
                        {{"Offset", 0.375}, {"GapDirection", 1}, {"Conductor", 2}}};
      stack["ConductorReferences"] = {{0.0625, 0.0, 0.0}, {0.3125, 0.0, 0.0}};
      // Its response matrices: the four basis points + one conductor state (two
      // conductors), a diagonal (PSD) matrix in the fixtures' CSV forms (domain: basis_i,
      // basis_j, Q_ij; surface: per interface 1 and edge, Q_ij and the whole-box
      // Q_total_ij).
      const auto stack_domain = temp.temp_dir / "mirror-pad-stack-domain.csv";
      const auto stack_surface = temp.temp_dir / "mirror-pad-stack-surface.csv";
      {
        std::ofstream domain(stack_domain), surface(stack_surface);
        domain << "basis_i,basis_j,Q_ij (J)\n";
        surface << "interface,edge,R (m),basis_i,basis_j,Q_ij (J),Q_total_ij (J)\n";
        for (int i = 1; i <= 5; i++)
        {
          for (int j = i; j <= 5; j++)
          {
            const double value = i == j ? 1.0e-12 : 0.0;
            domain << i << "," << j << "," << value << "\n";
            for (int edge = 1; edge <= 4; edge++)
            {
              surface << "1," << edge << "," << R << "," << i << "," << j << ","
                      << 0.25 * value << "," << 0.25 * value << "\n";
            }
          }
        }
      }
      for (const char *key : {"FabricatedMatrix", "ThinMatrix"})
      {
        stack[key] = stack_domain.string();
      }
      for (const char *key : {"FabricatedSurfaceMatrix", "ThinSurfaceMatrix"})
      {
        stack[key] = stack_surface.string();
      }
      stack_library["Models"].push_back(stack);
    }
    std::ofstream stack_output(stack_library_path);
    stack_output << stack_library.dump(2) << "\n";
  }
  Mpi::Barrier(Mpi::World());

  struct Energies
  {
    double raw = 0.0, outside = 0.0;
    std::array<double, 3> corrected{};  // ft, ff, sc
    std::array<double, 2> models{};     // evaluation 0, 1
    double uncovered_ft = 0.0, domain_boundary_ft = 0.0;
    // Per applied patch (0-based index): its ft SA energy and its cell length (the first
    // source, evaluation 0), from the PatchEnergy export.
    std::map<int, std::pair<double, double>> patch_energy;
    // Per model name: its ft SA energy (evaluation 0), from
    // surface-response-model-energy.csv and the record's ModelCatalog.
    std::map<std::string, double> model_energy;
    json diagnostics;
    json manifest_summary;
    // The identification manifest's Features and Segments (the dry run of the case).
    json manifest_features;
    json manifest_segments;
  };
  auto Run = [&](const std::string &name, double s, bool full, const std::string &mirror,
                 const std::vector<PadRectangle> &pads = kDefaultPads, int nz = 8,
                 const PadBend &bend = {}, int library_variant = 0)
  {
    const fs::path mesh_path = temp.temp_dir / (name + ".mesh");
    if (Mpi::Root(Mpi::World()))
    {
      mfem::Mesh serial = MakeSerialMirrorPadMesh(s, full, pads, nz, bend);
      std::ofstream output(mesh_path);
      // Full precision: the sheared coordinates (s = 1 / tan theta) are irrational, and a
      // wall written at 6 digits is no longer one plane for the truncation-plane fit.
      output.precision(17);
      serial.Print(output);
    }
    Mpi::Barrier(Mpi::World());
    json config = IslandConfig();
    config["Model"]["Mesh"] = mesh_path.string();
    config["Model"]["L0"] = 1.0;
    // The ground is the far face y = 1 (and its image); every other wall is natural. The
    // pad's short edges meet the walls x = 0 and x = 1 (and, in the full mesh, the image of
    // x = 0); the walls z = 0 / z = 2 are 0.75 > 3 R from the pad (outside the band).
    config["Boundaries"]["Ground"]["Attributes"] = {4};
    auto &correction = config["Solver"]["Electrostatic"]["ResponseCorrection"];
    correction.erase("PatchConstruction");
    const fs::path &variant_library_path =
        library_variant == 3 ? wedge_library_path
        : library_variant == 2
            ? stack_library_path
            : (library_variant == 1 ? curved_library_path : library_path);
    correction["Library"] = variant_library_path.string();
    if (library_variant == 3)
    {
      // The wedge library's radius (the decision-557 section).
      config["Boundaries"]["Postprocessing"]["Dielectric"][0]["EdgeDistances"] = {0.1};
    }
    correction["UnmatchedPolicy"] = "Warn";
    correction["TraceCoupling"] = "SurfaceMortar";
    correction["MortarOversampling"] = 2;
    correction["CorrectionMode"] = "Both";
    correction["PatchEnergy"] = true;
    correction["DomainBoundary"] = {{"Mirror", mirror}};
    config["Solver"]["Linear"] = {{"Tol", 1.0e-13}, {"MaxIts", 600}};
    const fs::path dir = temp.temp_dir / name;
    config["Problem"]["Output"] = dir.string();
    // The identification manifest of the case (the dry run; its Features are read by the
    // decision-512 sections, which run on the curved library, and exported as evidence),
    // beside the output directory.
    const fs::path manifest_path = temp.temp_dir / (name + "-dryrun") / "requirements.json";
    if (library_variant != 0 || std::getenv("PALACE_SYMMETRY_DEBUG_DIR"))
    {
      fs::create_directories(manifest_path.parent_path());
      Mpi::Barrier(Mpi::World());
      IoData iodata(config, false);
      std::vector<std::unique_ptr<Mesh>> meshes;
      mfem::Mesh serial_manifest = MakeSerialMirrorPadMesh(s, full, pads, nz, bend);
      meshes.push_back(std::make_unique<Mesh>(
          std::make_unique<mfem::ParMesh>(Mpi::World(), serial_manifest)));
      WriteSurfaceResponseRequirements(iodata, *meshes.front(), manifest_path.string());
      if (const char *debug_dir = std::getenv("PALACE_SYMMETRY_DEBUG_DIR");
          debug_dir && Mpi::Root(Mpi::World()))
      {
        // The dry run before the solve (a failing solve still leaves its record).
        fs::create_directories(fs::path(debug_dir) / name);
        fs::copy_file(manifest_path, fs::path(debug_dir) / name / "requirements.json",
                      fs::copy_options::overwrite_existing);
        fs::copy_file(manifest_path.parent_path() / "surface-response-patches.csv",
                      fs::path(debug_dir) / name / "dryrun-patches.csv",
                      fs::copy_options::overwrite_existing);
      }
    }
    if (const char *debug_dir = std::getenv("PALACE_SYMMETRY_DEBUG_DIR"))
    {
      // Evidence export (the lane record): the dry run, mesh, library and config of this
      // case into the given directory, and its postpro after the run.
      const fs::path debug = fs::path(debug_dir) / name;
      fs::create_directories(debug);
      if (Mpi::Root(Mpi::World()))
      {
        fs::copy_file(manifest_path, debug / "requirements.json",
                      fs::copy_options::overwrite_existing);
        fs::copy_file(mesh_path, debug / "mesh.mesh", fs::copy_options::overwrite_existing);
        fs::copy_file(variant_library_path, debug / "library.json",
                      fs::copy_options::overwrite_existing);
        fs::copy_file(basis_path, debug / "basis.csv",
                      fs::copy_options::overwrite_existing);
        json debug_config = config;
        debug_config["Model"]["Mesh"] = "mesh.mesh";
        debug_config["Solver"]["Electrostatic"]["ResponseCorrection"]["Library"] =
            "library.json";
        debug_config["Problem"]["Output"] = "postpro";
        std::ofstream output(debug / "config.json");
        output << debug_config.dump(2) << "\n";
      }
    }
    test::RunElectrostatic(config);
    Mpi::Barrier(Mpi::World());
    Energies e;
    if (!Mpi::Root(Mpi::World()))
    {
      return e;
    }
    if (const char *debug_dir = std::getenv("PALACE_SYMMETRY_DEBUG_DIR"))
    {
      fs::copy(dir, fs::path(debug_dir) / name / "postpro",
               fs::copy_options::recursive | fs::copy_options::overwrite_existing);
    }
    {
      const Table surface = test::LoadCsv(dir / "surface-Q.csv");
      const Table domain = test::LoadCsv(dir / "domain-E.csv");
      e.raw = test::ColumnByHeader(surface, "p_surf[4]").data[0] *
              test::ColumnByHeader(domain, "E_elec (J)").data[0];
      const Table edge = test::LoadCsv(dir / "surface-Q-edge.csv");
      e.outside = test::ColumnByHeader(edge, "E_out (J)").data[0];
      const Table corrected = test::LoadCsv(dir / "surface-Q-corrected.csv");
      e.corrected = {
          test::ColumnByHeader(corrected, "E_surf postprocessed fixed-trace[4] (J)")
              .data[0],
          test::ColumnByHeader(corrected, "E_surf postprocessed fixed-flux[4] (J)").data[0],
          test::ColumnByHeader(corrected, "E_surf corrected[4] (J)").data[0]};
      const Table models = test::LoadCsv(dir / "surface-response-model-energy.csv");
      const Column &eval = test::ColumnByHeader(models, "evaluation");
      const Column &energy =
          test::ColumnByHeader(models, "fabricated surface energy[4] (J)");
      for (std::size_t i = 0; i < eval.data.size(); i++)
      {
        const int evaluation = static_cast<int>(std::lround(eval.data[i]));
        if (evaluation < 2)
        {
          e.models[evaluation] += energy.data[i];
        }
      }
      auto TotalRow = [&](const fs::path &path, const std::string &label)
      {
        if (!fs::is_regular_file(path))
        {
          return 0.0;
        }
        const test::TextCsv csv = test::ReadTextCsv(path);
        const std::size_t type = csv.Column("type"), evaluation = csv.Column("evaluation");
        const std::size_t energy_column = csv.Column(label + " raw energy[4] (J)");
        for (const auto &row : csv.rows)
        {
          if (row[type] == "Total" && std::stoi(row[evaluation]) == 0)
          {
            return std::stod(row[energy_column]);
          }
        }
        return 0.0;
      };
      e.uncovered_ft = TotalRow(dir / "surface-response-uncovered-energy.csv", "uncovered");
      e.domain_boundary_ft =
          TotalRow(dir / "surface-response-domain-boundary-energy.csv", "domainboundary");
      {
        const Table patches = test::LoadCsv(dir / "surface-response-patch-energy.csv");
        const Column &eval = test::ColumnByHeader(patches, "evaluation");
        const Column &index = test::ColumnByHeader(patches, "patch");
        const Column &begin = test::ColumnByHeader(patches, "cell begin (m)");
        const Column &end = test::ColumnByHeader(patches, "cell end (m)");
        const Column &energy =
            test::ColumnByHeader(patches, "fabricated surface energy[4] (J)");
        for (std::size_t i = 0; i < eval.data.size(); i++)
        {
          if (static_cast<int>(std::lround(eval.data[i])) == 0)
          {
            e.patch_energy[static_cast<int>(std::lround(index.data[i])) - 1] = {
                energy.data[i], end.data[i] - begin.data[i]};
          }
        }
      }
      std::ifstream metadata_input(dir / "palace.json");
      REQUIRE(metadata_input);
      const json metadata = json::parse(metadata_input);
      e.diagnostics = metadata.at("SurfaceResponse").at("Diagnostics");
      if (fs::is_regular_file(manifest_path))
      {
        std::ifstream manifest_input(manifest_path);
        const json manifest = json::parse(manifest_input);
        e.manifest_features = manifest.at("Identification").at("Features");
        e.manifest_segments = manifest.at("Identification").at("Segments");
      }
      {
        std::map<int, std::string> names;
        for (const auto &model : metadata.at("SurfaceResponse").at("ModelCatalog"))
        {
          names[model.at("Index").get<int>()] = model.at("Name").get<std::string>();
        }
        const Table models = test::LoadCsv(dir / "surface-response-model-energy.csv");
        const Column &eval = test::ColumnByHeader(models, "evaluation");
        const Column &index = test::ColumnByHeader(models, "model");
        const Column &energy =
            test::ColumnByHeader(models, "fabricated surface energy[4] (J)");
        for (std::size_t i = 0; i < eval.data.size(); i++)
        {
          if (static_cast<int>(std::lround(eval.data[i])) == 0)
          {
            e.model_energy[names.at(static_cast<int>(std::lround(index.data[i])))] +=
                energy.data[i];
          }
        }
      }
    }
    return e;
  };
  auto CheckIdentity = [&](const Energies &e)
  {
    CHECK_THAT(e.corrected[0],
               WithinRel(e.outside + e.models[0] + e.uncovered_ft + e.domain_boundary_ft,
                         1.0e-10));
    CHECK_THAT(e.corrected[1],
               WithinRel(e.outside + e.models[1] + e.uncovered_ft + e.domain_boundary_ft,
                         1.0e-10));
  };

  SECTION("S-45: the virtual 90-degree corners (s_half = R, exact to 1e-9)")
  {
    const Energies full = Run("mirror-pad-45-full", 1.0, true, "Natural");
    const Energies half = Run("mirror-pad-45-half", 1.0, false, "Natural");
    // The checks run on the root; every rank stays in the section (a return inside a
    // SECTION would end the test body early on that rank and desynchronise Catch2's
    // section discovery across the ranks: the next pass then hangs in a collective).
    if (Mpi::Root(Mpi::World()))
    {
      CheckIdentity(full);
      CheckIdentity(half);
      // The raw interface energies and E_out: the half is exactly half the full (the
      // discrete half problem is the restriction of the symmetric full one).
      CHECK_THAT(2.0 * half.raw, WithinRel(full.raw, 1.0e-9));
      CHECK_THAT(2.0 * half.outside, WithinRel(full.outside, 1.0e-9));
      // Nothing is a DomainBoundary exclusion in the half: its cut-crossing cells are
      // Mirrored, the two virtual corners carry weight 1 / 2 and a mirror arm trim.
      const auto &exclusions = half.diagnostics.at("DomainBoundaryExclusions");
      CHECK(exclusions.at("Count").get<int>() == 0);
      CHECK(exclusions.at("Mirrored").at("Count").get<int>() > 0);
      CHECK(exclusions.at("Mirrored").at("Points").get<int>() > 0);
      // Four virtual corners (two per wall: convex at the top edge, concave at the bottom).
      CHECK(exclusions.at("Mirrored").at("HalfVertices").get<int>() == 4);
      CHECK(exclusions.at("Mirrored").at("ArmTrims").get<int>() == 4);
      const auto &band = half.diagnostics.at("MirrorBand");
      INFO(band.dump());
      REQUIRE(band.at("MirrorFormedFeatures").size() == 4);
      std::map<std::string, int> formed;
      for (const auto &entry : band.at("MirrorFormedFeatures"))
      {
        formed[entry.at("Type").get<std::string>()]++;
        CHECK(entry.at("Status") == "Modelled");
        CHECK_THAT(entry.at("Key").get<std::string>(),
                   ContainsSubstring("\"AngleDegrees\":90.0"));
      }
      CHECK(formed["ConvexCorner"] == 2);
      CHECK(formed["ConcaveCorner"] == 2);
      CHECK(half.domain_boundary_ft == 0.0);
      // The full has no cut through the pad: nothing mirrored, nothing excluded.
      // The full has no cut through the pad at the (former) wall x = 1; its ends at x = 0
      // and the image wall are virtual corners on both sides alike.
      CHECK(full.diagnostics.at("DomainBoundaryExclusions").at("Count").get<int>() == 0);
      CHECK(full.diagnostics.at("MirrorBand").at("MirrorFormedFeatures").size() == 4);
      // THE IDENTITY: the half's model energies, uncovered energy and corrected energies
      // are half the full's (the half-energy identity of DESIGN 2.2.4; the s_half term is
      // zero at 90 degrees).
      CHECK_THAT(2.0 * half.models[0], WithinRel(full.models[0], 1.0e-9));
      CHECK_THAT(2.0 * half.models[1], WithinRel(full.models[1], 1.0e-9));
      CHECK_THAT(2.0 * half.uncovered_ft, WithinRel(full.uncovered_ft, 1.0e-9));
      CHECK_THAT(2.0 * half.corrected[0], WithinRel(full.corrected[0], 1.0e-9));
      CHECK_THAT(2.0 * half.corrected[1], WithinRel(full.corrected[1], 1.0e-9));
      // Per model (isolated, convex-corner-90, concave-corner-90): the half's two
      // half-weight virtual corners at the cut + two at the outer wall against the full's
      // two real corners + four virtual ones, class by class (the aggregates above could
      // hide a transfer between classes; patch indices differ between the meshes, so the
      // comparison is per model, not per patch).
      REQUIRE(full.model_energy.size() == 3);
      REQUIRE(half.model_energy.size() == 3);
      for (const auto &[name, energy] : full.model_energy)
      {
        INFO("model " << name);
        REQUIRE(half.model_energy.count(name) == 1);
        CHECK_THAT(2.0 * half.model_energy.at(name), WithinRel(energy, 1.0e-9));
      }
      // The self-consistent column: the symmetric sc solution's restriction (accepted on
      // both or unavailable on both). 1e-5, not 1e-9: the sc field is the converged iterate
      // of the self-consistent solve (its stopping tolerance, not roundoff, sets the two
      // meshes' agreement); the ft / ff columns above are direct evaluations.
      if (std::isfinite(full.corrected[2]) && std::isfinite(half.corrected[2]))
      {
        CHECK_THAT(2.0 * half.corrected[2], WithinRel(full.corrected[2], 1.0e-5));
      }
    }
  }

  SECTION("S-22.5-wedge (decision 557 (4)): the mirror-formed 2-edge wedge clusters are "
          "emitted with their contract and placed on their real half, half = full / 2")
  {
    // The pad's short edges meet the walls at theta = 22.5 degrees (s = 1 / tan theta):
    // each edge and its image form a 45-degree wedge - a virtual corner (R of each arm at
    // the apex, formed, Missing here: no 45-degree model) plus the two arms beyond the
    // corner window, within 2 R of each other: a SpatialEdgeCluster that no mirror
    // placement merges (Unmerged; the O1 f17 class). The top edge's wedge is a 2-edge
    // cluster (the concave wedge, RealPortions one of two); the bottom edge's wedge joins
    // the pad's left-end corners into a 4-edge cluster (two real portions of four). Both
    // are EMITTED as Missing requirements with their contract, then PLACED on their real
    // half once the library carries a Signature-keyed model of their key: the coupon's own
    // energy on the half is exactly half the full's (the full's real wedge clusters carry
    // the same keys), the touched real cells leave the DomainBoundary raw term (owned by
    // the coupon).
    const double s = 1.0 / std::tan(22.5 * M_PI / 180.0);
    // The pad [0.5, 1] x [0.75, 1.25] at R = 0.1 (the wedge library's radius): the corner
    // windows R and the wedge clusters' arms 2.8 R from the apex leave an isolated stretch
    // before the real corners of the left end; the left end's perpendicular distance 0.5 /
    // sqrt(1 + s^2) = 0.19 from the sheared wall x = 0 is beyond R and its corners' images
    // beyond 2 R: that wall forms nothing.
    const std::vector<PadRectangle> wedge_pads = {{0.5, 0.75, 1.0, 1.25}};
    constexpr double R_wedge = 0.1;
    auto WriteWedgeLibrary = [&](const std::map<std::string, json> &wedge_keys)
    {
      std::ifstream input(curved_library_path);
      REQUIRE(input);
      json library = json::parse(input);
      library["Name"] = "unit-test-process-mirror-pad-wedge";
      library["MatchingRadius"] = R_wedge;
      library["CouponDepth"] = R_wedge;
      json corner;
      // The fixtures' surface matrices carry their within-R rows at R = 0.2 m: the same
      // rows re-keyed at R_wedge (the symmetric placement is the test, not the physics).
      auto RekeyedSurfaceMatrix = [&](const std::string &source)
      {
        const fs::path target =
            temp.temp_dir / (fs::path(source).stem().string() + "-R0.1.csv");
        std::ifstream matrix(source);
        REQUIRE(matrix);
        std::ofstream output(target);
        std::string line;
        std::getline(matrix, line);
        output << line << "\n";
        while (std::getline(matrix, line))
        {
          // interface,edge,R (m),... : the third field.
          std::size_t first = line.find(','), second = line.find(',', first + 1),
                      third = line.find(',', second + 1);
          REQUIRE(third != std::string::npos);
          output << line.substr(0, second + 1) << R_wedge << line.substr(third) << "\n";
        }
        return target.string();
      };
      for (auto &model : library["Models"])
      {
        if (model.contains("CouponDepth"))
        {
          model["CouponDepth"] = R_wedge;
        }
        for (const char *key : {"FabricatedSurfaceMatrix", "ThinSurfaceMatrix"})
        {
          if (model.contains(key))
          {
            model[key] = RekeyedSurfaceMatrix(model[key].get<std::string>());
          }
        }
        if (model["Name"] == "convex-corner-90")
        {
          corner = model;
        }
      }
      REQUIRE(!corner.is_null());
      // The 45-degree wedge corners (convex at the top edge, concave at the bottom):
      // relabelled corner coupons, so that the apex virtual corners are Modelled
      // (HalfByMirror weight 1 / 2) and their ownership by the wedge coupons (rule B4,
      // decision 584 (2) MINOR-7) is exercised: on the half the virtual corner, on the full
      // the real one, owned alike.
      for (const auto &[topology, name] : std::vector<std::pair<std::string, std::string>>{
               {"ConvexCorner", "convex-corner-45"},
               {"ConcaveCorner", "concave-corner-45"}})
      {
        json wedge_corner = corner;
        wedge_corner["Name"] = name;
        wedge_corner["Topology"] = topology;
        wedge_corner["Angle"] = 45.0;
        library["Models"].push_back(wedge_corner);
      }
      // One Signature-keyed model per wedge key: its Edges the Signature's portions in the
      // canonical frame (ChordedSignatureEdges), its matching volume the Signature's Box (R
      // above and below the plane), on the corner coupon's matrices with its 12 knots
      // (three rings of four) on the box's corners at z = -R, 0, R - the matching surface
      // spans the box as a built coupon's does, so the continuation ownership (the box =
      // the bbox of the basis points) and the domain-boundary test (the image half's
      // corners beyond the plane) read the box. No ENTRY stamp (an ordinary coupon of the
      // key, as the O1 model b0b764b21b95 is): the per-Edge weights come from the
      // identification's contract.
      for (const auto &[hash, signature] : wedge_keys)
      {
        json wedge = corner;
        wedge["Name"] = "wedge-cluster-" + hash.substr(0, 12);
        wedge["Topology"] = "SpatialEdgeCluster";
        wedge.erase("Angle");
        wedge.erase("AngleTolerance");
        wedge["Signature"] = signature;
        wedge["Interfaces"] = SignatureInterfaces(signature);
        wedge["Edges"] = ChordedSignatureEdges(signature, R_wedge);
        wedge["EdgePositionTolerance"] = 1.0e-6;
        wedge["EdgeAngleTolerance"] = 1.0e-6;
        const auto box = signature.at("Box").get<std::array<double, 4>>();
        json support = json::array();
        for (const double x : {box[0], box[2]})
        {
          for (const double y : {box[1], box[3]})
          {
            for (const double z : {-1.0, 1.0})
            {
              support.push_back({x * R_wedge, y * R_wedge, z * R_wedge});
            }
          }
        }
        wedge["SupportPoints"] = support;
        const fs::path basis =
            temp.temp_dir / (wedge["Name"].get<std::string>() + "-basis.csv");
        {
          std::ofstream points(basis);
          points << "x,y,z\n";
          for (const double z : {-1.0, 0.0, 1.0})
          {
            for (const auto &[x, y] :
                 std::vector<std::pair<double, double>>{{box[0], box[1]},
                                                        {box[2], box[1]},
                                                        {box[2], box[3]},
                                                        {box[0], box[3]}})
            {
              points << x * R_wedge << "," << y * R_wedge << "," << z * R_wedge << "\n";
            }
          }
        }
        wedge["BasisPoints"] = basis.string();
        library["Models"].push_back(wedge);
      }
      std::ofstream output(wedge_library_path);
      output << library.dump(2) << "\n";
    };
    auto IsConfiguration = [](const json &feature)
    {
      return feature.contains("Mirror") && feature.at("Mirror").is_object() &&
             feature.at("Mirror").value("Status", "") == "Unmerged" &&
             feature.at("Mirror").contains("ExtendedFeature");
    };
    // Step 1: the half without a wedge model: the configurations are Missing requirements
    // WITH their contract; their own cells are DomainBoundary (UnmergedTopology, decision
    // 481), the raw term positive.
    if (Mpi::Root(Mpi::World()))
    {
      WriteWedgeLibrary({});
    }
    Mpi::Barrier(Mpi::World());
    const Energies missing =
        Run("mirror-pad-wedge-half-missing", s, false, "Natural", wedge_pads, 8, {}, 3);
    std::map<std::string, json> keys;  // hash -> the configuration's Signature
    std::set<int> touched_ids;         // the real features the configurations touch
    if (Mpi::Root(Mpi::World()))
    {
      CheckIdentity(missing);
      INFO(missing.diagnostics.at("MirrorBand").dump());
      std::map<int, int> edge_counts;
      for (const auto &feature : missing.manifest_features)
      {
        if (!IsConfiguration(feature))
        {
          if (feature.contains("Mirror") && feature.at("Mirror").is_object() &&
              feature.at("Mirror").value("Status", "") == "Unmerged" &&
              feature.at("Type") == "IsolatedEdge")
          {
            // The touched real features with cells (a touched real cluster, Missing, has
            // uncovered portions clipped by the coupon's box instead).
            touched_ids.insert(feature.at("Id").get<int>());
          }
          continue;
        }
        CHECK(feature.at("Type") == "SpatialEdgeCluster");
        CHECK(feature.at("Match").at("Status") == "Missing");
        REQUIRE(feature.at("Mirror").contains("Contract"));
        const auto &contract = feature.at("Mirror").at("Contract");
        const int edge_count = feature.at("Signature").at("EdgeCount").get<int>();
        edge_counts[edge_count] = static_cast<int>(contract.at("RealPortions").size());
        CHECK(contract.at("RealPortions").size() > 0);
        CHECK(contract.at("RealPortions").size() <
              feature.at("Signature").at("Portions").size());
        {
          nlohmann::json frame = contract.at("Frame");
          // the recorded frame's handedness (= Features[].Chirality for a chiral key)
          CHECK((frame.at("Chirality").get<int>() == 1 ||
                 frame.at("Chirality").get<int>() == -1));
          CHECK((feature.at("Chirality").get<int>() == 0 ||
                 frame.at("Chirality").get<int>() == feature.at("Chirality").get<int>()));
          frame.erase("Chirality");
          CHECK(frame == feature.at("Frame"));
        }
        CHECK_THAT(contract.at("RealLengthOverR").get<double>(),
                   WithinRel(contract.at("ImageLengthOverR").get<double>(), 1.0e-6));
        keys[feature.at("Hash").get<std::string>()] = feature.at("Signature");
      }
      CHECK(keys.size() == 2);
      CHECK(edge_counts == std::map<int, int>{{2, 1}, {4, 2}});
      // The apex virtual 45-degree corners are Modelled at weight 1 / 2 (HalfByMirror) and
      // owned by nothing while the configurations are Missing.
      int half_corners = 0;
      for (const auto &feature : missing.manifest_features)
      {
        if ((feature.at("Type") == "ConvexCorner" ||
             feature.at("Type") == "ConcaveCorner") &&
            feature.contains("Mirror") && feature.at("Mirror").is_object() &&
            feature.at("Mirror").value("Status", "") == "Modelled")
        {
          half_corners++;
          CHECK_THAT(feature.at("Signature").at("AngleDegrees").get<double>(),
                     WithinAbs(45.0, 1.0e-6));
        }
      }
      // ONE apex corner feature: the top edge's convex 45-degree wedge; the bottom edge's
      // concave wedge apex is absorbed into the 4-edge cluster (its corner sites joined the
      // cluster: no corner feature of its own).
      CHECK(half_corners == 1);
      CHECK(missing.diagnostics.at("ContinuationOwnership")
                .at("Vertices")
                .at("Count")
                .get<int>() == 0);
      CHECK(!touched_ids.empty());
      CHECK(missing.domain_boundary_ft > 0.0);
      const auto &exclusions = missing.diagnostics.at("DomainBoundaryExclusions");
      CHECK(exclusions.at("Reasons").contains("UnmergedTopology"));
      for (const auto &patch : exclusions.at("Patches"))
      {
        CHECK(patch.at("Reason") == "UnmergedTopology");
        CHECK(touched_ids.count(patch.at("Feature").get<int>()) == 1);
      }
      CHECK(!missing.diagnostics.contains("MirrorFormedPlacement"));
      WriteWedgeLibrary(keys);
    }
    Mpi::Barrier(Mpi::World());
    // Step 2: the half and the full with the wedge library.
    const Energies full =
        Run("mirror-pad-wedge-full", s, true, "Natural", wedge_pads, 8, {}, 3);
    const Energies half =
        Run("mirror-pad-wedge-half", s, false, "Natural", wedge_pads, 8, {}, 3);
    if (Mpi::Root(Mpi::World()))
    {
      CheckIdentity(full);
      CheckIdentity(half);
      INFO(half.diagnostics.at("DomainBoundaryExclusions").dump());
      CHECK_THAT(2.0 * half.raw, WithinRel(full.raw, 1.0e-9));
      CHECK_THAT(2.0 * half.outside, WithinRel(full.outside, 1.0e-9));
      // The configurations are Matched (the note names the real-half placement); the full's
      // real wedge clusters carry the same keys.
      std::set<std::string> half_models, full_models;
      for (const auto &feature : half.manifest_features)
      {
        if (IsConfiguration(feature))
        {
          CHECK(feature.at("Match").at("Status") == "Matched");
          CHECK_THAT(feature.at("Match").at("Note").get<std::string>(),
                     ContainsSubstring("REAL half"));
          half_models.insert(feature.at("Match").at("Model").get<std::string>());
        }
      }
      for (const auto &feature : full.manifest_features)
      {
        if (feature.at("Type") == "SpatialEdgeCluster")
        {
          CHECK(feature.at("Match").at("Status") == "Matched");
          CHECK(!IsConfiguration(feature));
          full_models.insert(feature.at("Match").at("Model").get<std::string>());
        }
      }
      CHECK(half_models.size() == 2);
      CHECK(half_models == full_models);
      // The placement record: one coupon per configuration at weight 1 / 2 (the real length
      // fraction of a symmetric key), its edges split real / image, its claims real only.
      REQUIRE(half.diagnostics.contains("MirrorFormedPlacement"));
      const auto &placement = half.diagnostics.at("MirrorFormedPlacement");
      REQUIRE(placement.at("Count").get<int>() == 2);
      std::set<std::size_t> coupon_patches;
      for (const auto &entry : placement.at("Patches"))
      {
        coupon_patches.insert(entry.at("Patch").get<std::size_t>());
        CHECK_THAT(entry.at("Weight").get<double>(), WithinAbs(0.5, 1.0e-6));
        CHECK(entry.at("RealEdges").get<int>() > 0);
        CHECK(entry.at("ImageEdges").get<int>() > 0);
        CHECK(entry.at("RealEdges").get<int>() == entry.at("ImageEdges").get<int>());
        CHECK(entry.at("RealClaims").size() >=
              static_cast<std::size_t>(entry.at("RealEdges").get<int>()));
        CHECK(entry.at("Reach").get<double>() > 3.0 * R_wedge);
        CHECK(half_models.count(entry.at("Model").get<std::string>()) == 1);
      }
      CHECK(!full.diagnostics.contains("MirrorFormedPlacement"));
      // Nothing is DomainBoundary any more: the coupons are Mirrored (their image half
      // evaluated by even extension within their reach), the touched real features' cells
      // inside the box are owned by the coupon (weight 0), the raw term is zero.
      const auto &exclusions = half.diagnostics.at("DomainBoundaryExclusions");
      CHECK(exclusions.at("Count").get<int>() == 0);
      CHECK(half.domain_boundary_ft == 0.0);
      CHECK(full.domain_boundary_ft == 0.0);
      std::set<std::size_t> mirrored;
      for (const auto &entry : exclusions.at("Mirrored").at("Patches"))
      {
        mirrored.insert(entry.at("Patch").get<std::size_t>());
      }
      for (const std::size_t patch : coupon_patches)
      {
        CHECK(mirrored.count(patch) == 1);
      }
      std::set<int> owned_features;
      for (const auto &cell : half.diagnostics.at("ContinuationOwnership").at("OwnedCells"))
      {
        owned_features.insert(cell.at("Feature").get<int>());
        REQUIRE(cell.at("Owners").size() == 1);
        CHECK(coupon_patches.count(
                  cell.at("Owners")[0].at("SpatialPatch").get<std::size_t>()) == 1);
      }
      for (const int id : touched_ids)
      {
        CHECK(owned_features.count(id) == 1);
      }
      // MINOR-7 (decision 584 (2)): the apex virtual corners (HalfByMirror 0.5) are owned
      // by the wedge coupons under rule B4 (their vertex is the coupon's chain-piece end on
      // the plane inside its box; the full symmetric signature models the apex) - weight 0
      // in the dry run - exactly as the FULL owns its two real 45-degree corners by its
      // real wedge clusters: the identity half = full / 2 holds with the corners owned
      // alike (their model energy 0 on both). Predicted on O1 f35 (the virtual 22.5-degree
      // apex corner owned by the b0b764b21b95 coupon, as the real instance's apex corner
      // f16 is by the real coupon).
      auto OwnedVertices = [&](const Energies &e)
      {
        std::map<std::string, std::set<std::size_t>>
            owners;  // corner model -> owner patches
        for (const auto &record :
             e.diagnostics.at("ContinuationOwnership").at("Vertices").at("Records"))
        {
          CHECK(record.at("Kind") == "Vertex");
          for (const auto &owner : record.at("Owners"))
          {
            owners[record.at("Model").get<std::string>()].insert(
                owner.at("SpatialPatch").get<std::size_t>());
          }
        }
        return owners;
      };
      const auto half_owned = OwnedVertices(half), full_owned = OwnedVertices(full);
      // ONE apex corner feature (the top edge's convex wedge; the bottom edge's apex is
      // inside the 4-edge cluster): owned on the half by the configuration's coupon, on the
      // full by the real wedge cluster's, its model energy 0 on both.
      CHECK(half.diagnostics.at("ContinuationOwnership")
                .at("Vertices")
                .at("Count")
                .get<int>() == 1);
      CHECK(full.diagnostics.at("ContinuationOwnership")
                .at("Vertices")
                .at("Count")
                .get<int>() == 1);
      REQUIRE(half_owned.count("convex-corner-45") == 1);
      REQUIRE(full_owned.count("convex-corner-45") == 1);
      for (const std::size_t owner : half_owned.at("convex-corner-45"))
      {
        CHECK(coupon_patches.count(owner) == 1);
      }
      REQUIRE(half.model_energy.count("convex-corner-45") == 1);
      REQUIRE(full.model_energy.count("convex-corner-45") == 1);
      CHECK(half.model_energy.at("convex-corner-45") == 0.0);
      CHECK(full.model_energy.at("convex-corner-45") == 0.0);
      // THE IDENTITY on the coupons: each wedge model's energy on the half is half the
      // full's (the full's real cluster at weight 1 against the half's configuration at
      // weight 1 / 2 with its image half sampled by even extension). The isolated model's
      // energy is not compared: the 4-edge key's box (grown on the device plan in its
      // canonical frame) is not symmetric about the plane, so the full itself owns the two
      // mirrored left edges' cells unequally - the full's own treatment, not the half's.
      for (const auto &name : half_models)
      {
        INFO("model " << name);
        REQUIRE(half.model_energy.count(name) == 1);
        REQUIRE(full.model_energy.count(name) == 1);
        CHECK(half.model_energy.at(name) > 0.0);
        // The half's coupon carries the real length fraction f (exactly 1 / 2 for an
        // exactly symmetric key; the 4-edge key's serialised lengths differ by a few
        // quanta): the half's energy = f x the full's. 1e-8: the model energies are read
        // back from the energy CSV's printed digits.
        double fraction = 0.0;
        for (const auto &entry : placement.at("Patches"))
        {
          if (entry.at("Model") == name)
          {
            fraction = entry.at("Weight").get<double>();
          }
        }
        CHECK_THAT(fraction, WithinAbs(0.5, 1.0e-5));
        CHECK_THAT(half.model_energy.at(name),
                   WithinRel(fraction * full.model_energy.at(name), 1.0e-8));
      }
      // The patch dry run carries the DomainWeight column (decision 559 MINOR-3 / O-14).
      std::ifstream dry_run(temp.temp_dir / "mirror-pad-wedge-half-dryrun" /
                            "surface-response-patches.csv");
      REQUIRE(dry_run);
      std::string header;
      std::getline(dry_run, header);
      CHECK_THAT(header, EndsWith(",StripBegin,StripEnd,DomainWeight"));
      std::string row;
      int rows = 0;
      while (std::getline(dry_run, row))
      {
        rows++;
        CHECK_THAT(row, EndsWith(",1"));
      }
      CHECK(rows > 0);
    }
    Mpi::Barrier(Mpi::World());
    // Step 3 (the ENTRY stamp, impl-B5 CONTRACT.md section 3; the dry run alone): the wedge
    // models stamped with the consumer's MirrorFormed record and per-Edge Weights. A stamp
    // that agrees with the identification's contract places exactly as the unstamped
    // coupon; a stamp whose RealPortions are the IMAGE portions (self-consistent with its
    // own Weights, so it loads) fails closed at the placement by name; an Edge Weight
    // disagreeing with the stamp's RealPortions fails closed at the library load.
    {
      // The contracts of the half's configurations by key (from step 1's manifest: the
      // identification is library-independent).
      std::map<std::string, json> contracts;
      if (Mpi::Root(Mpi::World()))
      {
        for (const auto &feature : missing.manifest_features)
        {
          if (IsConfiguration(feature))
          {
            contracts[feature.at("Hash").get<std::string>()] =
                feature.at("Mirror").at("Contract");
          }
        }
        REQUIRE(contracts.size() == 2);
      }
      auto StampedLibrary = [&](const std::string &variant)
      {
        // The path on every rank (the config below is built on every rank); written on
        // the root.
        const fs::path path =
            temp.temp_dir / ("fabrication-process-mirror-pad-wedge-" + variant + ".json");
        if (!Mpi::Root(Mpi::World()))
        {
          return path;
        }
        std::ifstream input(wedge_library_path);
        REQUIRE(input);
        json library = json::parse(input);
        for (auto &model : library["Models"])
        {
          if (model["Topology"] != "SpatialEdgeCluster")
          {
            continue;
          }
          const std::string hash =
              SignatureKeyAndHash(model["Signature"], "SpatialEdgeCluster").second;
          REQUIRE(contracts.count(hash) == 1);
          const auto &contract = contracts.at(hash);
          std::vector<int> real_portions =
              contract.at("RealPortions").get<std::vector<int>>();
          const std::size_t portions = model["Signature"]["Portions"].size();
          if (variant == "image-portions")
          {
            // The complement: the stamp claims the image portions as real.
            std::vector<int> complement;
            for (std::size_t i = 0; i < portions; i++)
            {
              if (std::find(real_portions.begin(), real_portions.end(),
                            static_cast<int>(i)) == real_portions.end())
              {
                complement.push_back(static_cast<int>(i));
              }
            }
            real_portions = complement;
          }
          // ChordedSignatureEdges: one Edge per straight portion, in portion order.
          REQUIRE(model["Edges"].size() == portions);
          json edge_portions = json::array();
          for (std::size_t e = 0; e < portions; e++)
          {
            const bool real = std::find(real_portions.begin(), real_portions.end(),
                                        static_cast<int>(e)) != real_portions.end();
            model["Edges"][e]["Weight"] = real ? 1.0 : 0.0;
            edge_portions.push_back(e);
          }
          if (variant == "flipped-weight")
          {
            model["Edges"][0]["Weight"] = 1.0 - model["Edges"][0]["Weight"].get<double>();
          }
          const double real_over_R = contract.at("RealLengthOverR").get<double>();
          const double image_over_R = contract.at("ImageLengthOverR").get<double>();
          model["MirrorFormed"] = {
              {"Version", 1},
              {"RealPortions", real_portions},
              {"EdgePortions", edge_portions},
              {"Planes", contract.at("Planes")},
              {"RealLengthOverR", real_over_R},
              {"ImageLengthOverR", image_over_R},
              {"RealLengthFraction", real_over_R / (real_over_R + image_over_R)},
              {"Rule", "unit test: the consumer's ENTRY stamp (impl-B5 CONTRACT.md s3)"}};
        }
        std::ofstream output(path);
        output << library.dump(2) << "\n";
        return path;
      };
      auto Preflight = [&](const fs::path &library, const std::string &name)
      {
        json config = IslandConfig();
        const fs::path mesh_path = temp.temp_dir / "mirror-pad-wedge-half.mesh";
        config["Model"]["Mesh"] = mesh_path.string();
        config["Model"]["L0"] = 1.0;
        config["Boundaries"]["Ground"]["Attributes"] = {4};
        config["Boundaries"]["Postprocessing"]["Dielectric"][0]["EdgeDistances"] = {
            R_wedge};
        auto &correction = config["Solver"]["Electrostatic"]["ResponseCorrection"];
        correction.erase("PatchConstruction");
        correction["Library"] = library.string();
        correction["UnmatchedPolicy"] = "Warn";
        correction["TraceCoupling"] = "SurfaceMortar";
        correction["DomainBoundary"] = {{"Mirror", "Natural"}};
        config["Problem"]["Output"] = (temp.temp_dir / name).string();
        const fs::path manifest_path = temp.temp_dir / name / "requirements.json";
        fs::create_directories(manifest_path.parent_path());
        IoData iodata(config, false);
        mfem::Mesh serial = MakeSerialMirrorPadMesh(s, false, wedge_pads, 8, {});
        Mesh mesh(std::make_unique<mfem::ParMesh>(Mpi::World(), serial));
        WriteSurfaceResponseRequirements(iodata, mesh, manifest_path.string());
        Mpi::Barrier(Mpi::World());
        // The manifest is written on the root; every rank returns (a non-root rank reading
        // a file the root writes would desynchronise the ranks).
        if (!Mpi::Root(Mpi::World()))
        {
          return json();
        }
        std::ifstream manifest_input(manifest_path);
        return json::parse(manifest_input);
      };
      const fs::path consistent = StampedLibrary("consistent");
      const fs::path image_portions = StampedLibrary("image-portions");
      const fs::path flipped = StampedLibrary("flipped-weight");
      Mpi::Barrier(Mpi::World());
      {
        const json manifest = Preflight(consistent, "mirror-pad-wedge-stamped-consistent");
        if (Mpi::Root(Mpi::World()))
        {
          const auto &placement =
              manifest.at("Identification").at("Diagnostics").at("MirrorFormedPlacement");
          CHECK(placement.at("Count").get<int>() == 2);
          for (const auto &entry : placement.at("Patches"))
          {
            CHECK_THAT(entry.at("Weight").get<double>(), WithinAbs(0.5, 1.0e-6));
          }
          int contract_requirements = 0;
          for (const auto &requirement : manifest.at("Requirements"))
          {
            if (!requirement.contains("MirrorFormedContract"))
            {
              continue;
            }
            contract_requirements++;
            CHECK(requirement.at("Status") == "Exact");
            CHECK(requirement.at("MirrorFormed").get<bool>());
            CHECK(requirement.at("Topology") == "SpatialEdgeCluster");
            CHECK(requirement.at("MirrorFormedContract").at("Version").get<int>() == 1);
          }
          CHECK(contract_requirements == 2);
        }
      }
      CHECK_THROWS_WITH(Preflight(image_portions, "mirror-pad-wedge-stamped-image"),
                        ContainsSubstring("the identification's contract reads the portion "
                                          "as an IMAGE") &&
                            ContainsSubstring("inconsistent MirrorFormed record"));
      CHECK_THROWS_WITH(Preflight(flipped, "mirror-pad-wedge-stamped-flipped"),
                        ContainsSubstring("carries a MirrorFormed ENTRY record but edge"));
    }
  }

  SECTION("S-90: a perpendicular meeting forms nothing (Continued) and is half the full")
  {
    const Energies full = Run("mirror-pad-90-full", 0.0, true, "Natural");
    const Energies half = Run("mirror-pad-90-half", 0.0, false, "Natural");
    // The checks run on the root; every rank stays in the section (a return inside a
    // SECTION would end the test body early on that rank and desynchronise Catch2's
    // section discovery across the ranks: the next pass then hangs in a collective).
    if (Mpi::Root(Mpi::World()))
    {
      CheckIdentity(full);
      CheckIdentity(half);
      const auto &band = half.diagnostics.at("MirrorBand");
      CHECK(band.at("MirrorFormedFeatures").empty());
      CHECK(band.at("ContinuedFeatures").get<int>() == 2);  // the two short edges
      CHECK(half.diagnostics.at("DomainBoundaryExclusions").at("Count").get<int>() == 0);
      CHECK_THAT(2.0 * half.raw, WithinRel(full.raw, 1.0e-9));
      CHECK_THAT(2.0 * half.outside, WithinRel(full.outside, 1.0e-9));
      CHECK_THAT(2.0 * half.models[0], WithinRel(full.models[0], 1.0e-9));
      CHECK_THAT(2.0 * half.uncovered_ft, WithinRel(full.uncovered_ft, 1.0e-9));
      CHECK_THAT(2.0 * half.corrected[0], WithinRel(full.corrected[0], 1.0e-9));
      CHECK_THAT(2.0 * half.corrected[1], WithinRel(full.corrected[1], 1.0e-9));
    }
  }

  // The oblique cases with s_half > R (decision 473 (2): the trim path and the second-order
  // term): theta = atan(1 / s) exactly 60 / 76 degrees, virtual corners 120 / 152 degrees
  // (relabelled exact models), s = R / max(|cos|, |sin|) = 1.1547 R / 1.1326 R, s_half =
  // 1.0774 R / 1.0663 R. The identity 2 x half - full = sum over the virtual corners of
  // (E[s_half, s] - E[s_half, R]) on the real arm (DESIGN 2.2.3), bounded by the arm's
  // energy density at s_half x (s - R) / 2: measured on the half's first kept arm cell (the
  // MirrorArmTrim record names the trimmed cells) and recorded with the measured residual.
  auto ObliqueCase = [&](const std::string &label, double theta_degrees)
  {
    const double s = 1.0 / std::tan(theta_degrees * M_PI / 180.0);
    const Energies full = Run("mirror-pad-" + label + "-full", s, true, "Natural");
    const Energies half = Run("mirror-pad-" + label + "-half", s, false, "Natural");
    if (Mpi::Root(Mpi::World()))
    {
      CheckIdentity(full);
      CheckIdentity(half);
      CHECK_THAT(2.0 * half.raw, WithinRel(full.raw, 1.0e-9));
      CHECK_THAT(2.0 * half.outside, WithinRel(full.outside, 1.0e-9));
      const auto &exclusions = half.diagnostics.at("DomainBoundaryExclusions");
      INFO(half.diagnostics.at("MirrorBand").dump());
      INFO(half.diagnostics.at("MirrorArmTrim").dump());
      CHECK(exclusions.at("Count").get<int>() == 0);
      CHECK(exclusions.at("Mirrored").at("HalfVertices").get<int>() == 4);
      const auto &trims = half.diagnostics.at("MirrorArmTrim");
      REQUIRE(trims.at("Count").get<int>() == 4);
      const double angle = 2.0 * theta_degrees;
      const double exit = R / std::max(std::abs(std::cos(angle * M_PI / 180.0)),
                                       std::abs(std::sin(angle * M_PI / 180.0)));
      const double s_half = 0.5 * (R + exit);
      double bound = 0.0;
      int trimmed_cells = 0;
      for (const auto &corner : trims.at("Corners"))
      {
        CHECK_THAT(corner.at("AngleDegrees").get<double>(), WithinAbs(angle, 1.0e-6));
        CHECK_THAT(corner.at("HalfStartOverR").get<double>(),
                   WithinAbs(s_half / R, 1.0e-9));
        CHECK_THAT(corner.at("TrimmedLength").get<double>(), WithinAbs(s_half - R, 1.0e-9));
        CHECK(corner.at("Weight").get<double>() == 0.5);
        // The real arm's first cell was clipped to begin at s_half (not wholly removed at
        // this cell size); its energy density bounds the second-order term.
        REQUIRE(!corner.at("Patches").empty());
        double density = 0.0;
        for (const auto &patch : corner.at("Patches"))
        {
          const auto it = half.patch_energy.find(patch.get<int>());
          REQUIRE(it != half.patch_energy.end());
          CHECK(it->second.second > 0.0);
          density = std::max(density, it->second.first / it->second.second);
          trimmed_cells++;
        }
        bound += density * (exit - R) / 2.0;
      }
      CHECK(trimmed_cells >= 4);
      const double residual = 2.0 * half.models[0] - full.models[0];
      INFO("residual 2 x half - full (models, ft) = " << residual << " J; bound " << bound
                                                      << " J; full " << full.models[0]);
      if (const char *debug_dir = std::getenv("PALACE_SYMMETRY_DEBUG_DIR"))
      {
        // The lane record: the measured second-order residual beside its bound.
        std::ofstream record(fs::path(debug_dir) / ("symmetry-" + label + ".json"));
        record << json{{"Case", label},
                       {"ThetaDegrees", theta_degrees},
                       {"CornerAngleDegrees", angle},
                       {"ExitOverR", exit / R},
                       {"HalfStartOverR", s_half / R},
                       {"FullModels", full.models[0]},
                       {"TwiceHalfModels", 2.0 * half.models[0]},
                       {"Residual", residual},
                       {"ResidualRelative", residual / full.models[0]},
                       {"Bound", bound},
                       {"BoundRelative", bound / full.models[0]},
                       {"FullFt", full.corrected[0]},
                       {"TwiceHalfFt", 2.0 * half.corrected[0]},
                       {"FtResidualRelative",
                        (2.0 * half.corrected[0] - full.corrected[0]) / full.corrected[0]}}
                      .dump(2)
               << "\n";
      }
      CHECK(std::abs(residual) <= bound + 1.0e-9 * full.models[0]);
      CHECK_THAT(2.0 * half.uncovered_ft, WithinRel(full.uncovered_ft, 1.0e-9));
      CHECK(std::abs(2.0 * half.corrected[0] - full.corrected[0]) <=
            bound + 1.0e-9 * full.corrected[0]);
      CHECK(std::abs(2.0 * half.corrected[1] - full.corrected[1]) <=
            bound + 1.0e-9 * full.corrected[1]);
    }
  };

  SECTION("S-60: virtual 120-degree corners, s_half = 1.0774 R (the trim path)")
  {
    ObliqueCase("60", 60.0);
  }

  SECTION("S-76: virtual 152-degree corners, s_half = 1.0663 R")
  {
    ObliqueCase("76", 76.0);
  }

  // The parallel cases (s = 0, the wall x = 1): an edge parallel to the wall at d = h =
  // 0.125
  // (< R) with the gap toward the wall (the slot's left edge: same-conductor GAP 2 d = 0.25
  // with the image) or with metal toward the wall (the L's strip along the wall: same-
  // conductor STRIP 0.25). The pair is placed on its real side only (side factor 1 / 2; the
  // image side frames the midline); the full has the real pair. Exact: <= 1e-9.
  // The pair must be longer than the cluster reach at its ends (2 R each side): a slot /
  // strip of 1.25 - 1.375 (> 6 R) on the nz = 16 grid (h_z = 0.125), its middle a clean
  // pair.
  auto ParallelCase = [&](const std::string &label, const std::vector<PadRectangle> &pads,
                          const std::string &type)
  {
    const Energies full =
        Run("mirror-pad-" + label + "-full", 0.0, true, "Natural", pads, 16);
    const Energies half =
        Run("mirror-pad-" + label + "-half", 0.0, false, "Natural", pads, 16);
    if (Mpi::Root(Mpi::World()))
    {
      CheckIdentity(full);
      CheckIdentity(half);
      const auto &band = half.diagnostics.at("MirrorBand");
      INFO(band.dump());
      bool formed = false;
      for (const auto &entry : band.at("MirrorFormedFeatures"))
      {
        if (entry.at("Type") == type)
        {
          formed = true;
          CHECK(entry.at("Status") == "Modelled");
        }
      }
      CHECK(formed);
      CHECK_THAT(2.0 * half.raw, WithinRel(full.raw, 1.0e-9));
      CHECK_THAT(2.0 * half.outside, WithinRel(full.outside, 1.0e-9));
      // THE PAIR IDENTITY: the half's real side at its side factor = half the full's pair.
      const std::string model = type + "-0.25";
      REQUIRE(half.model_energy.count(model) == 1);
      REQUIRE(full.model_energy.count(model) == 1);
      CHECK_THAT(2.0 * half.model_energy.at(model),
                 WithinRel(full.model_energy.at(model), 1.0e-9));
      // MISSING-EQUIVALENCE (decision 481): the unmerged mirror-formed cluster at the
      // pair's end is read as a Missing feature - its OWN cells raw (every DomainBoundary
      // patch overlaps one of its real portions; none is a pair cell), its NEIGHBOURS
      // applied: the pair cells whose support reaches the cluster (UnmergedSupportReach,
      // distance
      // <= R) keep their model, exactly as next to the full's real Missing cluster - which
      // is why the pair identity above is exact on both parallel cases.
      const auto &half_exclusions = half.diagnostics.at("DomainBoundaryExclusions");
      REQUIRE(!half.diagnostics.at("MirrorBand").at("UnmergedFeatures").empty());
      std::set<int> excluded_patches;
      for (const auto &patch : half_exclusions.at("Patches"))
      {
        excluded_patches.insert(patch.at("Patch").get<int>());
        CHECK(patch.at("Reason") == "UnmergedTopology");
        CHECK(patch.at("UnmergedTopology").at("OwnFootprint").get<bool>());
        CHECK(patch.at("Topology") != "same-conductor gap");
        CHECK(patch.at("Topology") != "same-conductor strip");
      }
      int pair_cells_reaching = 0;
      for (const auto &reach : half_exclusions.at("UnmergedSupportReach").at("Patches"))
      {
        CHECK(reach.at("Distance").get<double>() <= R + 1.0e-9);
        if (!excluded_patches.count(reach.at("Patch").get<int>()))
        {
          pair_cells_reaching++;  // applied although its support reaches the cluster
        }
      }
      CHECK(pair_cells_reaching >= 2);  // one pair cell at each end of the pair
      const bool clean_ends = half_exclusions.at("Count").get<int>() == 0;
      if (clean_ends)
      {
        // No blocked patch at the pair's ends (the slot: its real clusters at the slot's
        // mouths are Missing on both meshes - raw uncovered - and the Unmerged
        // mirror-formed clusters touch them alone): the whole window is exact.
        CHECK_THAT(2.0 * half.models[0], WithinRel(full.models[0], 1.0e-9));
        CHECK_THAT(2.0 * half.models[1], WithinRel(full.models[1], 1.0e-9));
        CHECK_THAT(2.0 * half.uncovered_ft, WithinRel(full.uncovered_ft, 1.0e-9));
        CHECK_THAT(2.0 * half.corrected[0], WithinRel(full.corrected[0], 1.0e-9));
        CHECK_THAT(2.0 * half.corrected[1], WithinRel(full.corrected[1], 1.0e-9));
      }
      else
      {
        // The strip's ends are no symmetric-identity region: the L junction (the half's
        // real cluster at the bridge; the full's two-corner cluster of the finger root) and
        // the strip's top against the wall (the top corner 0.125 from the wall with its
        // image: a two-vertex configuration, Unmerged -> DomainBoundary on the half) are
        // MISSING configurations whose raw claims differ between the two meshes by
        // construction - the F-DB-a BRACKET of DESIGN 2.2.5, not the identity. Recorded
        // (the whole-window residual beside the pair identity above); the Reasons record
        // names the UnmergedTopology route. Note also the sub-quadrature cells: a
        // DomainBoundary cell shorter than the quadrature spacing of the interface integral
        // reads 0 raw energy (the F2 nearest-foot filter is quadrature-resolved; irrelevant
        // at production cell sizes, recorded here).
        const double raw_kept = 2.0 * (half.uncovered_ft + half.domain_boundary_ft) +
                                full.uncovered_ft + full.domain_boundary_ft;
        const double residual = 2.0 * half.corrected[0] - full.corrected[0];
        INFO("unmerged ends: residual " << residual << " J; raw kept " << raw_kept);
        CHECK(half.diagnostics.at("DomainBoundaryExclusions")
                  .at("Reasons")
                  .contains("UnmergedTopology"));
        CHECK(std::abs(residual) < 0.05 * full.corrected[0]);
        if (const char *debug_dir = std::getenv("PALACE_SYMMETRY_DEBUG_DIR"))
        {
          std::ofstream record(fs::path(debug_dir) / ("symmetry-" + label + ".json"));
          record << json{{"Case", label},
                         {"PairModel", model},
                         {"FullPair", full.model_energy.at(model)},
                         {"TwiceHalfPair", 2.0 * half.model_energy.at(model)},
                         {"FullFt", full.corrected[0]},
                         {"TwiceHalfFt", 2.0 * half.corrected[0]},
                         {"Residual", residual},
                         {"RawKept", raw_kept},
                         {"HalfDomainBoundary", half.domain_boundary_ft},
                         {"FullUncovered", full.uncovered_ft}}
                        .dump(2)
                 << "\n";
        }
      }
    }
  };

  SECTION(
      "S-par-gap: a slot edge parallel to the wall at d = 0.125 (gap 0.25 with the image)")
  {
    // The pad z in [0.25, 1.75] with a slot x in [0.875, 1], z in [0.375, 1.625]: three
    // rectangles.
    ParallelCase(
        "par-gap",
        {{0.0, 0.25, 0.875, 1.75}, {0.875, 0.25, 1.0, 0.375}, {0.875, 1.625, 1.0, 1.75}},
        "SameConductorGap");
  }

  SECTION("S-par-strip: a strip along the wall, d = 0.125 (strip 0.25 with the image)")
  {
    // An L: the bridge z in [0.25, 0.375] from x = 0 to the wall and the strip x in [0.875,
    // 1] up to z = 1.75.
    ParallelCase("par-strip", {{0.0, 0.25, 1.0, 0.375}, {0.875, 0.375, 1.0, 1.75}},
                 "SameConductorStrip");
  }

  SECTION("S-corner-box: a pair along one wall meeting the other wall (two planes)")
  {
    // The pad x in [0.25, 1], z in [1.0, 1.875] (nz = 16): its top edge runs 0.125 below
    // the z = 2 wall (a same-conductor GAP 0.25 with the z image, Modelled) and meets the x
    // = 1 wall perpendicularly (the straight joint continues the pair through the x image);
    // its top-left corner at (0.25, 1.875) with its z image is a two-vertex configuration
    // (Unmerged -> DomainBoundary, raw kept) 0.75 > 2 R from the x wall, so it is the SAME
    // configuration on both sides. Samples of the pair's cells near the box corner lie
    // beyond the z plane (and are reflected); the composition through both planes is
    // exercised by the pure-function test. The pair, raw, E_out and DomainBoundary terms
    // are exact (1e-9); the window to 1e-6 (measured 2.7e-7, recorded).
    const std::vector<PadRectangle> pads = {{0.25, 1.0, 1.0, 1.875}};
    const Energies full = Run("mirror-pad-box-full", 0.0, true, "Natural", pads, 16);
    const Energies half = Run("mirror-pad-box-half", 0.0, false, "Natural", pads, 16);
    if (Mpi::Root(Mpi::World()))
    {
      CheckIdentity(full);
      CheckIdentity(half);
      INFO(half.diagnostics.at("MirrorBand").dump());
      INFO(half.diagnostics.at("DomainBoundaryExclusions").dump());
      const auto &band = half.diagnostics.at("MirrorBand");
      bool gap = false;
      for (const auto &entry : band.at("MirrorFormedFeatures"))
      {
        gap = gap ||
              (entry.at("Type") == "SameConductorGap" && entry.at("Status") == "Modelled");
      }
      CHECK(gap);
      CHECK(!band.at("UnmergedFeatures").empty());
      const auto &exclusions = half.diagnostics.at("DomainBoundaryExclusions");
      CHECK(exclusions.at("Count").get<int>() > 0);
      CHECK(exclusions.at("Reasons").contains("UnmergedTopology"));
      CHECK(exclusions.at("Mirrored").at("Points").get<int>() > 0);
      CHECK_THAT(2.0 * half.raw, WithinRel(full.raw, 1.0e-9));
      CHECK_THAT(2.0 * half.outside, WithinRel(full.outside, 1.0e-9));
      REQUIRE(half.model_energy.count("SameConductorGap-0.25") == 1);
      REQUIRE(full.model_energy.count("SameConductorGap-0.25") == 1);
      CHECK_THAT(2.0 * half.model_energy.at("SameConductorGap-0.25"),
                 WithinRel(full.model_energy.at("SameConductorGap-0.25"), 1.0e-9));
      CHECK_THAT(2.0 * half.domain_boundary_ft, WithinRel(full.domain_boundary_ft, 1.0e-9));
      // Measured 2.7e-7 relative (the gap pair, the raw / E_out terms and the
      // DomainBoundary term are exact to 1e-9; the residual sits in the isolated-edge /
      // corner cells Mirrored through the z plane on both meshes): recorded, bounded at
      // 1e-6.
      CHECK_THAT(
          2.0 * (half.models[0] + half.uncovered_ft + half.domain_boundary_ft),
          WithinRel(full.models[0] + full.uncovered_ft + full.domain_boundary_ft, 1.0e-6));
      CHECK_THAT(2.0 * half.corrected[0], WithinRel(full.corrected[0], 1.0e-6));
      CHECK_THAT(2.0 * half.corrected[1], WithinRel(full.corrected[1], 1.0e-6));
      if (const char *debug_dir = std::getenv("PALACE_SYMMETRY_DEBUG_DIR"))
      {
        std::ofstream record(fs::path(debug_dir) / "symmetry-box.json");
        record << json{{"Case", "box"},
                       {"FullFt", full.corrected[0]},
                       {"TwiceHalfFt", 2.0 * half.corrected[0]},
                       {"ResidualRelative",
                        (2.0 * half.corrected[0] - full.corrected[0]) / full.corrected[0]},
                       {"HalfDomainBoundary", half.domain_boundary_ft},
                       {"FullDomainBoundary", full.domain_boundary_ft},
                       {"HalfReasons", exclusions.at("Reasons")}}
                      .dump(2)
               << "\n";
      }
    }
  }

  // Decision 512 (the boundary-cut CHAINS lane; DESIGN ERRATA-7): the rerun-2 J1 finding -
  // a real chain with curved sub-portions (BendRadiusOverR) straight at one plane and
  // clipped by a virtual corner at another was refused as Continued (its extended portions
  // differed from the unextended ones) and became an Unmerged IsolatedEdge whose own cells
  // (rule 481) went DomainBoundary raw, the WHOLE chain (S5 1163.5 um). The fixture: the
  // default pad with its edges along x bent in z by t(x) (PadBend: a slope per x cell).
  SECTION("S-90-bent: a bent chain straight at the symmetry plane and oblique at the far "
          "wall is Continued; half = full / 2 exact")
  {
    // t(x): slope 1 / tan 76 on [0, 0.25], 0.6 of it on [0.25, 0.5], 0.2 of it on [0.5,
    // 0.75], 0 beyond: the edges meet the far wall x = 0 at theta = 76 degrees (virtual
    // 152-degree corners, Modelled: the S-76 class on BOTH meshes alike, so they cancel in
    // 2 x half - full) and turn by ~0.1 rad at x = 0.25 / 0.5 / 0.75 (sub-noise joints:
    // implied sagitta 0.003 < 0.05 R = 0.01 on 0.25-long pieces; the turn density of a
    // joint is spread over its two half-runs: 0.38-0.39 on [0.125, 0.625], so the windowed
    // bend radius stays above 10 R everywhere (~12.8 R): ONE IsolatedEdge per edge WITH
    // curved sub-portions, BendRadiusOverR set and PortionTurns non-zero - the S5 class,
    // no CurvedEdge section), then run straight into the plane x = 1 (a straight joint).
    // The image continues the last run collinearly, so the extended run spreads the x =
    // 0.75 joint's turn over a longer half-run (0.13 instead of 0.2 per unit): a turn
    // difference far inside the joint noise rule, tolerated by the Continued test.
    const double slope = 1.0 / std::tan(76.0 * M_PI / 180.0);
    const PadBend bend = {slope,       slope,       0.6 * slope, 0.6 * slope,
                          0.2 * slope, 0.2 * slope, 0.0,         0.0};
    const Energies full =
        Run("mirror-pad-bent-full", 0.0, true, "Natural", kDefaultPads, 8, bend, 1);
    const Energies half =
        Run("mirror-pad-bent-half", 0.0, false, "Natural", kDefaultPads, 8, bend, 1);
    if (Mpi::Root(Mpi::World()))
    {
      CheckIdentity(full);
      CheckIdentity(half);
      const auto &band = half.diagnostics.at("MirrorBand");
      INFO(band.dump());
      INFO(half.diagnostics.at("DomainBoundaryExclusions").dump());
      // The two chains (the pad's bottom and top edges) are Continued at x = 1 although
      // their far ends are clipped by the virtual corners at x = 0 (the unextended
      // IsolatedEdge runs from the wall, the extended one from the corner's claim: the
      // exact-portions test of record refused this); nothing is Unmerged, no real feature
      // is touched, no cell is DomainBoundary. The chains are ONE IsolatedEdge each with
      // curved sub-portions (no real corner, no CurvedEdge section).
      CHECK(band.at("ContinuedFeatures").get<int>() == 2);
      {
        int isolated = 0, curved = 0, real_corners = 0;
        for (const auto &feature : half.manifest_features)
        {
          const std::string type = feature.at("Type").get<std::string>();
          const bool mirror_formed = feature.contains("Mirror") &&
                                     !feature.at("Mirror").is_null() &&
                                     feature.at("Mirror").at("Status") != "Continued";
          isolated += type == "IsolatedEdge" ? 1 : 0;
          curved += type == "CurvedEdge" ? 1 : 0;
          // A real corner would be a kink read as a corner instead of a sub-noise joint.
          real_corners +=
              ((type == "ConvexCorner" || type == "ConcaveCorner") && !mirror_formed) ? 1
                                                                                      : 0;
          if (type == "IsolatedEdge")
          {
            REQUIRE(!feature.at("BendRadiusOverR").is_null());
            CHECK(feature.at("BendRadiusOverR").get<double>() >= 10.0);
            CHECK(feature.at("BendRadiusOverR").get<double>() < 20.0);
            double total_turn = 0.0;
            for (const auto &turn : feature.at("PortionTurns"))
            {
              total_turn += std::abs(turn.get<double>());
            }
            CHECK(total_turn > 0.1);
            REQUIRE(!feature.at("Mirror").is_null());
            CHECK(feature.at("Mirror").at("Status") == "Continued");
          }
        }
        CHECK(isolated == 2);
        CHECK(curved == 0);
        CHECK(real_corners == 0);
      }
      CHECK(band.at("UnmergedFeatures").empty());
      CHECK(band.at("TouchedRealFeatures").empty());
      std::map<std::string, int> formed;
      for (const auto &entry : band.at("MirrorFormedFeatures"))
      {
        formed[entry.at("Type").get<std::string>()]++;
        CHECK(entry.at("Status") == "Modelled");
        CHECK_THAT(entry.at("Key").get<std::string>(),
                   ContainsSubstring("\"AngleDegrees\":152.0"));
      }
      CHECK(formed["ConvexCorner"] == 1);
      CHECK(formed["ConcaveCorner"] == 1);
      const auto &exclusions = half.diagnostics.at("DomainBoundaryExclusions");
      CHECK(exclusions.at("Count").get<int>() == 0);
      CHECK(exclusions.at("Mirrored").at("Count").get<int>() > 0);
      CHECK(exclusions.at("Mirrored").at("HalfVertices").get<int>() == 2);
      CHECK(half.domain_boundary_ft == 0.0);
      // THE IDENTITY, exact: raw, E_out, models, uncovered, ft, ff, and per model class.
      CHECK_THAT(2.0 * half.raw, WithinRel(full.raw, 1.0e-9));
      CHECK_THAT(2.0 * half.outside, WithinRel(full.outside, 1.0e-9));
      CHECK_THAT(2.0 * half.models[0], WithinRel(full.models[0], 1.0e-9));
      CHECK_THAT(2.0 * half.models[1], WithinRel(full.models[1], 1.0e-9));
      CHECK_THAT(2.0 * half.uncovered_ft, WithinRel(full.uncovered_ft, 1.0e-9));
      CHECK_THAT(2.0 * half.corrected[0], WithinRel(full.corrected[0], 1.0e-9));
      CHECK_THAT(2.0 * half.corrected[1], WithinRel(full.corrected[1], 1.0e-9));
      // Per class: the virtual corners to 1e-9 of the total (they hold ~1e-2 of the models,
      // at roundoff level between the half's and the reflected full's coordinates); the
      // straight-like chain's anchor ("isolated") and first-order curvature node ("...kappa
      // 0.1...") TOGETHER to 1e-9 - the Continued chain keeps the UNEXTENDED turns
      // (ERRATA-4 MINOR-2) while the full reads the joint at x = 0.75 spread over its run
      // continued through the plane: the first-order weights of the last run's cells
      // differ, a reallocation between the anchor and the node that this fixture's family
      // (every node the anchor) leaves exactly invariant in total.
      REQUIRE(half.model_energy.size() == full.model_energy.size());
      double half_chain = 0.0, full_chain = 0.0;
      for (const auto &[name, energy] : full.model_energy)
      {
        INFO("model " << name);
        REQUIRE(half.model_energy.count(name) == 1);
        if (name.find("corner") != std::string::npos)
        {
          CHECK_THAT(2.0 * half.model_energy.at(name),
                     WithinAbs(energy, 1.0e-9 * full.models[0]));
        }
        else
        {
          half_chain += half.model_energy.at(name);
          full_chain += energy;
        }
      }
      CHECK(full_chain > 0.0);
      CHECK_THAT(2.0 * half_chain, WithinRel(full_chain, 1.0e-9));
      if (const char *debug_dir = std::getenv("PALACE_SYMMETRY_DEBUG_DIR"))
      {
        std::ofstream record(fs::path(debug_dir) / "symmetry-bent.json");
        record << json{{"Case", "bent"},
                       {"FullFt", full.corrected[0]},
                       {"TwiceHalfFt", 2.0 * half.corrected[0]},
                       {"ResidualRelative",
                        (2.0 * half.corrected[0] - full.corrected[0]) / full.corrected[0]},
                       {"ContinuedFeatures", band.at("ContinuedFeatures")},
                       {"HalfDomainBoundary", half.domain_boundary_ft}}
                      .dump(2)
               << "\n";
      }
    }
  }

  SECTION("S-90-curved-end: a chain curved up to the plane - the unmerged joined bend owns "
          "the re-classified sliver only (decision 512 (c)), the accounting exact")
  {
    // The edges along x: slope -tan 30 on [0, 0.5] (theta = 60 degrees at the far wall x =
    // 0: virtual 120-degree corners, Modelled, on BOTH meshes alike), then a CIRCULAR arc
    // of radius rho = 1 = 5 R centred ON the plane x = 1 (the chord vertices at x = 0.5 ...
    // 1.0 lie exactly on the circle, tangent-continuous with the straight run at 0.5 and
    // perpendicular to the plane at 1.0: joints of ~0.13 rad, all of one sense and all
    // sub-noise; the arc rule fits the circle). The image of the arc is the SAME circle, so
    // the extended run reads one bend arc through the plane. Unextended run: the arc's
    // density ends at its last joint (x = 0.875) and the window is clipped at the chain
    // end, so the windowed bend radius leaves the curved regime before the plane: a
    // CurvedEdge (kappa = R / rho = 0.2, on the curved family) and an IsolatedEdge sliver
    // [~0.875, 1.0] before it. Extended run: the joined arc continues through the plane
    // (the same convexity, the same radius), so the CurvedEdge extends over the sliver -
    // ONE unmerged CurvedEdge per edge (no mirror placement) whose real portions the
    // unextended run identifies identically except that sliver. Rule 481 + decision 512
    // (c): the sliver's cells are the configuration's own (DomainBoundary, raw kept); the
    // CurvedEdge's cells keep their curved-family model (they would have been the
    // configuration's own under the whole-portion reading of ERRATA-5).
    PadBend bend(8, -std::tan(30.0 * M_PI / 180.0));
    for (int i = 4; i < 8; i++)
    {
      auto Circle = [](double x) { return -std::sqrt(1.0 - (1.0 - x) * (1.0 - x)); };
      bend[i] = (Circle((i + 1) * 0.125) - Circle(i * 0.125)) / 0.125;
    }
    const Energies full =
        Run("mirror-pad-curved-end-full", 0.0, true, "Natural", kDefaultPads, 8, bend, 1);
    const Energies half =
        Run("mirror-pad-curved-end-half", 0.0, false, "Natural", kDefaultPads, 8, bend, 1);
    if (Mpi::Root(Mpi::World()))
    {
      CheckIdentity(full);
      CheckIdentity(half);
      const auto &band = half.diagnostics.at("MirrorBand");
      const auto &exclusions = half.diagnostics.at("DomainBoundaryExclusions");
      INFO(band.dump());
      INFO(exclusions.dump());
      // The unextended reading: per edge one IsolatedEdge (the run from the far wall and
      // the sliver) and one CurvedEdge on the curved family; the virtual 120-degree
      // corners.
      std::set<int> isolated_ids, curved_ids, configuration_ids;
      for (const auto &feature : half.manifest_features)
      {
        const std::string type = feature.at("Type").get<std::string>();
        const bool configuration = feature.contains("Mirror") &&
                                   feature.at("Mirror").is_object() &&
                                   feature.at("Mirror").value("Status", "") == "Unmerged" &&
                                   feature.at("Mirror").contains("ExtendedFeature");
        if (configuration)
        {
          // The unmerged joined bend emitted as a feature with its contract (decision 557
          // (1)): a CurvedEdge key carries no Portions (RealPortions empty), no mirror
          // placement in this version (Missing, the note names the 2D family's consumer),
          // its image portions listed apart.
          configuration_ids.insert(feature.at("Id").get<int>());
          CHECK(type == "CurvedEdge");
          CHECK(feature.at("Match").at("Status") == "Missing");
          CHECK_THAT(feature.at("Match").at("Note").get<std::string>(),
                     ContainsSubstring("no mirror placement in this version"));
          const auto &contract = feature.at("Mirror").at("Contract");
          CHECK(contract.at("Version").get<int>() == 1);
          CHECK(contract.at("RealPortions").empty());
          CHECK(contract.at("RealLengthOverR").get<double>() > 0.0);
          CHECK(contract.at("ImageLengthOverR").get<double>() > 0.0);
          {
            nlohmann::json frame = contract.at("Frame");
            CHECK((frame.at("Chirality").get<int>() == 1 ||
                   frame.at("Chirality").get<int>() == -1));
            CHECK((feature.at("Chirality").get<int>() == 0 ||
                   frame.at("Chirality").get<int>() == feature.at("Chirality").get<int>()));
            frame.erase("Chirality");
            CHECK(frame == feature.at("Frame"));
          }
          CHECK(feature.contains("ImagePortions"));
          continue;
        }
        if (type == "IsolatedEdge")
        {
          isolated_ids.insert(feature.at("Id").get<int>());
        }
        else if (type == "CurvedEdge")
        {
          curved_ids.insert(feature.at("Id").get<int>());
          CHECK(feature.at("Match").at("Status") == "Matched");
          CHECK_THAT(feature.at("Match").at("Model").get<std::string>(),
                     ContainsSubstring("kappa"));
          // Not touched: identified identically by both runs.
          CHECK((!feature.contains("Mirror") || feature.at("Mirror").is_null()));
        }
      }
      CHECK(isolated_ids.size() == 2);
      CHECK(curved_ids.size() == 2);
      CHECK(configuration_ids.size() == 2);
      std::map<std::string, int> formed;
      for (const auto &entry : band.at("MirrorFormedFeatures"))
      {
        if (entry.at("Status") == "Modelled")
        {
          formed[entry.at("Type").get<std::string>()]++;
          CHECK_THAT(entry.at("Key").get<std::string>(),
                     ContainsSubstring("\"AngleDegrees\":120.0"));
        }
      }
      CHECK(formed["ConvexCorner"] == 1);
      CHECK(formed["ConcaveCorner"] == 1);
      // Two unmerged joined bends (the bottom and the top edge), each with identical real
      // portions (the unextended CurvedEdge's three chords [0.5, 0.875]) and a differing
      // one (the sliver [0.875, 1.0], read as the IsolatedEdge's by the unextended run):
      // 3 : 1 in length (the record's lengths are in the operator's mesh units).
      REQUIRE(band.at("UnmergedFeatures").size() == 2);
      double identical = 0.0, differing = 0.0;
      for (const auto &entry : band.at("UnmergedFeatures"))
      {
        CHECK(entry.at("Type") == "CurvedEdge");
        const double entry_identical = entry.at("IdenticalRealLength").get<double>();
        const double entry_differing = entry.at("DifferingRealLength").get<double>();
        CHECK(entry_differing > 0.0);
        CHECK_THAT(entry_identical / (entry_identical + entry_differing),
                   WithinAbs(0.76, 0.06));
        identical += entry_identical;
        differing += entry_differing;
        REQUIRE(entry.at("RealFeatures").size() == 1);
        CHECK(isolated_ids.count(entry.at("RealFeatures")[0].get<int>()) == 1);
        for (const auto &portion : entry.at("Portions"))
        {
          if (portion.at("Image").get<bool>())
          {
            continue;
          }
          if (portion.at("Identical").get<bool>())
          {
            CHECK(curved_ids.count(portion.at("RealFeature").get<int>()) == 1);
          }
          else
          {
            CHECK(portion.at("Differs") == "Type");
            CHECK(isolated_ids.count(portion.at("RealFeature").get<int>()) == 1);
          }
        }
      }
      std::set<int> touched;
      for (const auto &id : band.at("TouchedRealFeatures"))
      {
        touched.insert(id.get<int>());
      }
      CHECK(touched == isolated_ids);
      // The own cells: on the differing slivers only (every DomainBoundary patch is a cell
      // of a touched IsolatedEdge - its isolated-edge patch and the co-located first-order
      // "curved edge" blend patch of the same cell), never the CurvedEdge features' cells,
      // which keep their curved-family model (applied / Mirrored). The excluded own-edge
      // length (distinct geometric cells, in the config's units) is the two slivers' one
      // mesh segment each (0.125 along x, tilted by the last chord's slope 0.063).
      CHECK(exclusions.at("Count").get<int>() > 0);
      CHECK(exclusions.at("Reasons").contains("UnmergedTopology"));
      std::map<std::tuple<int, int, double, double>, double> own_cells;
      for (const auto &patch : exclusions.at("Patches"))
      {
        CHECK(patch.at("Reason") == "UnmergedTopology");
        CHECK(patch.at("UnmergedTopology").at("OwnFootprint").get<bool>());
        CHECK(isolated_ids.count(patch.at("Feature").get<int>()) == 1);
        const auto cell = patch.at("Cell").get<std::array<double, 2>>();
        own_cells[{patch.at("Feature").get<int>(), patch.at("Segment").get<int>(),
                   std::round(cell[0] * 1.0e9), std::round(cell[1] * 1.0e9)}] =
            cell[1] - cell[0];
      }
      double own_length = 0.0;
      for (const auto &[key, length] : own_cells)
      {
        own_length += length;
      }
      INFO("own cells " << own_cells.size() << ", length " << own_length << ", differing "
                        << differing);
      CHECK_THAT(own_length, WithinAbs(2.0 * 0.125 * std::hypot(1.0, bend[7]), 1.0e-6));
      bool curved_model = false;
      for (const auto &[name, energy] : half.model_energy)
      {
        if (name.find("kappa0.2") != std::string::npos)
        {
          curved_model = true;
          CHECK(energy > 0.0);
        }
      }
      CHECK(curved_model);
      CHECK(half.domain_boundary_ft > 0.0);
      // The symmetric terms are exact: raw, E_out. The window differs from half the full by
      // the slivers alone: raw on the half (the DomainBoundary term), modelled on the full
      // by its single CurvedEdge over the plane - 2 x half - full = 2 x DB_half - (the
      // full's curved-family energy on the slivers' mirror pair), the latter read at the
      // full's mean curved energy per length (the fixture's coupon response dwarfs the raw
      // field energy, so the residual is the whole model energy of the slivers; within 15 %
      // for the non-uniform density along the arc).
      CHECK_THAT(2.0 * half.raw, WithinRel(full.raw, 1.0e-9));
      CHECK_THAT(2.0 * half.outside, WithinRel(full.outside, 1.0e-9));
      const double residual = 2.0 * half.corrected[0] - full.corrected[0];
      double full_curved_energy = 0.0, full_curved_length = 0.0;
      for (const auto &[name, energy] : full.model_energy)
      {
        if (name.find("kappa0.2") != std::string::npos)
        {
          full_curved_energy += energy;
        }
      }
      for (const auto &feature : full.manifest_features)
      {
        if (feature.at("Type") == "CurvedEdge")
        {
          full_curved_length += feature.at("Length").get<double>();
        }
      }
      // The unmerged configurations are listed with their contract and a MergedFeature.
      for (const auto &entry : band.at("UnmergedFeatures"))
      {
        CHECK(configuration_ids.count(entry.at("MergedFeature").get<int>()) == 1);
        CHECK(!entry.contains("ContractRefused"));
      }
      REQUIRE(full_curved_length > 0.0);
      const double expected_residual =
          2.0 * half.domain_boundary_ft -
          full_curved_energy * (2.0 * own_length / full_curved_length);
      INFO("residual " << residual << " J, expected " << expected_residual
                       << " J; half DomainBoundary " << half.domain_boundary_ft);
      CHECK_THAT(residual, WithinRel(expected_residual, 0.15));
      if (const char *debug_dir = std::getenv("PALACE_SYMMETRY_DEBUG_DIR"))
      {
        std::ofstream record(fs::path(debug_dir) / "symmetry-curved-end.json");
        record << json{{"Case", "curved-end"},
                       {"FullFt", full.corrected[0]},
                       {"TwiceHalfFt", 2.0 * half.corrected[0]},
                       {"Residual", residual},
                       {"ExpectedResidual", expected_residual},
                       {"ResidualRelative", residual / full.corrected[0]},
                       {"HalfDomainBoundary", half.domain_boundary_ft},
                       {"IdenticalRealLength", identical},
                       {"DifferingRealLength", differing},
                       {"OwnCellLength", own_length},
                       {"ContinuedFeatures", band.at("ContinuedFeatures")},
                       {"Reasons", exclusions.at("Reasons")}}
                      .dump(2)
               << "\n";
      }
    }
  }

  // Decisions 536 / 537 (the boundary-cut STACKS lane; DESIGN ERRATA-8): the rerun-2 J2
  // solves of C2p / C2 / S1b / C3 aborted in the first postprocessing - "an uncovered
  // perimeter portion of feature N (parallel-edge cluster) lies within the tolerance of a
  // perimeter segment at its midpoint but not at its ends". A DomainBoundary translational
  // cell's raw claim (F-DB-a) was rebuilt from the patch FRAME, origin + EdgeOffset AxisU +
  // cell AxisW; for a pair / stack AxisW follows the partner / last side's chord (the far
  // foot clamped at a side's end next to an oblique plane), so the interval left the member
  // edge by cell x sin(angle) - 8-75 nm against the postprocessor's 1e-3 R. The raw claim
  // is now the cell's recorded pre-image on its OWN segment (Provenance::own_cell).
  SECTION("S-45-stack: a four-edge stack meeting the plane obliquely - its Unmerged "
          "cells are raw-kept on their member edges, the accounting exact")
  {
    // Two strips along x (z in [1.5, 1.625] and [1.75, 1.875]: edges 0.125 = 0.625 R apart,
    // a four-edge two-conductor ParallelEdgeCluster on the stack-4edge-0.125 model) meeting
    // both walls at 45 degrees: with its image the stack is a bent stack (no mirror
    // placement) -> an Unmerged SpatialEdgeCluster at each wall whose own cells (the stack
    // cells on the differing pieces) are DomainBoundary, raw kept; the same configuration
    // on both meshes at x = 0.
    // The default pad lowered to z in [0.5, 1.0] (0.5 = 2.5 R below the stack: outside the
    // stack rule's reach; the S-45 configuration: its cells and virtual corners are the
    // applied patches the fail-closed exclusion test requires) and the stack z in [1.5,
    // 1.625] / [1.75, 1.875] (its top edge 0.125 below the z = 2 wall: the z image joins
    // the configuration, Unmerged too). The stack, 5 R long between two oblique walls, is
    // wholly within its clusters' reach: every stack cell is a DomainBoundary own cell.
    const std::vector<PadRectangle> pads = {
        {0.0, 0.5, 1.0, 1.0}, {0.0, 1.5, 1.0, 1.625}, {0.0, 1.75, 1.0, 1.875}};
    const Energies half =
        Run("mirror-pad-stack-half", 1.0, false, "Natural", pads, 16, {}, 2);
    if (Mpi::Root(Mpi::World()))
    {
      // The accounting identity: ft = E_out + models + uncovered + domainboundary (the
      // DomainBoundary term integrated by the postprocessor over the raw claims, which the
      // verification at surfacepostoperator.cpp UncoveredEdgeIntervals accepted: every
      // claim end on a perimeter segment).
      CheckIdentity(half);
      const auto &band = half.diagnostics.at("MirrorBand");
      const auto &exclusions = half.diagnostics.at("DomainBoundaryExclusions");
      INFO(band.dump());
      INFO(exclusions.dump());
      // The stack feature, matched, and the Unmerged configurations touching it.
      // The MATCHED four-edge stack (at the oblique walls the sides end staggered, so the
      // stack rule also reads short three-edge / pair pieces there: Missing, raw
      // uncovered).
      int stack_id = -1;
      for (const auto &feature : half.manifest_features)
      {
        if (feature.at("Type") == "ParallelEdgeCluster" &&
            feature.at("Match").at("Status") == "Matched")
        {
          REQUIRE(stack_id < 0);
          stack_id = feature.at("Id").get<int>();
          CHECK(feature.at("Match").at("Model") == "stack-4edge-0.125");
          CHECK(feature.at("Signature").at("Edges").size() == 4);
          CHECK(feature.at("Length").get<double>() > 2.0);
        }
      }
      REQUIRE(stack_id >= 0);
      bool stack_unmerged = false;
      for (const auto &entry : band.at("UnmergedFeatures"))
      {
        for (const auto &id : entry.at("RealFeatures"))
        {
          stack_unmerged = stack_unmerged || id.get<int>() == stack_id;
        }
      }
      CHECK(stack_unmerged);
      // Its own cells are DomainBoundary (Reason UnmergedTopology), raw kept: the raw
      // portions of the parallel-edge cluster type exist and EVERY raw portion's ends lie
      // on its own identification segment (the pre-image), none dropped or skipped.
      CHECK(exclusions.at("Count").get<int>() > 0);
      CHECK(exclusions.at("Reasons").contains("UnmergedTopology"));
      const auto &raw = exclusions.at("RawPortions");
      CHECK(raw.at("ByType").contains("parallel-edge cluster"));
      double stack_raw_length = 0.0;
      int stack_raw_portions = 0;
      for (const auto &portion : raw.at("Portions"))
      {
        const int segment = portion.at("Segment").get<int>();
        REQUIRE(segment >= 0);
        const auto &key =
            half.manifest_segments.at(static_cast<std::size_t>(segment)).at("Key");
        auto Off = [&](const json &p)
        {
          std::array<double, 3> a = key[0], b = key[1], q = p;
          double ab2 = 0.0, t = 0.0;
          for (int d = 0; d < 3; d++)
          {
            ab2 += (b[d] - a[d]) * (b[d] - a[d]);
            t += (q[d] - a[d]) * (b[d] - a[d]);
          }
          t = std::clamp(t / ab2, 0.0, 1.0);
          double off2 = 0.0;
          for (int d = 0; d < 3; d++)
          {
            const double r = q[d] - (a[d] + t * (b[d] - a[d]));
            off2 += r * r;
          }
          return std::sqrt(off2);
        };
        CHECK(Off(portion.at("P0")) <= 1.0e-9);
        CHECK(Off(portion.at("P1")) <= 1.0e-9);
        if (portion.at("Type") == "parallel-edge cluster")
        {
          stack_raw_length += portion.at("Length").get<double>();
          stack_raw_portions++;
        }
      }
      CHECK(stack_raw_portions > 0);
      INFO("stack raw portions " << stack_raw_portions << ", length " << stack_raw_length);
      // On this fixture the stack's member edges are exactly parallel and its sides end
      // together, so the frame reconstruction of record coincides with the pre-image (the
      // defect needs a skewed frame: the S-bent-stack section next reproduces it end to
      // end, the synthetic case "Domain-boundary raw claim of a pair / stack cell" below at
      // the production skew); here the dry run's first-side cells confirm the two readings
      // agree on a parallel stack.
      const auto rows = ReadDryRun(temp.temp_dir / "mirror-pad-stack-half-dryrun" /
                                   "surface-response-patches.csv");
      double worst_frame_offset = 0.0;
      for (const auto &patch : exclusions.at("Patches"))
      {
        if (patch.at("Topology") != "parallel-edge cluster")
        {
          continue;
        }
        const std::size_t index = patch.at("Patch").get<std::size_t>();
        const auto row = std::find_if(rows.begin(), rows.end(),
                                      [&](const DryRunRow &r) { return r.patch == index; });
        REQUIRE(row != rows.end());
        const int segment = patch.at("Segment").get<int>();
        const auto &key =
            half.manifest_segments.at(static_cast<std::size_t>(segment)).at("Key");
        std::array<double, 3> a = key[0], b = key[1];
        double origin_off = 0.0, ab2 = 0.0, t = 0.0;
        for (int d = 0; d < 3; d++)
        {
          ab2 += (b[d] - a[d]) * (b[d] - a[d]);
          t += (row->origin[d] - a[d]) * (b[d] - a[d]);
        }
        t = std::clamp(t / ab2, 0.0, 1.0);
        for (int d = 0; d < 3; d++)
        {
          const double r = row->origin[d] - (a[d] + t * (b[d] - a[d]));
          origin_off += r * r;
        }
        if (std::sqrt(origin_off) > 1.0e-9)
        {
          continue;  // another side (EdgeOffset != 0): the frame formula needs it
        }
        for (const double c : row->strip)
        {
          double off2 = 0.0, tt = 0.0;
          std::array<double, 3> q{};
          for (int d = 0; d < 3; d++)
          {
            q[d] = row->origin[d] + c * row->axis_w[d];
            tt += (q[d] - a[d]) * (b[d] - a[d]);
          }
          tt = std::clamp(tt / ab2, 0.0, 1.0);
          for (int d = 0; d < 3; d++)
          {
            const double r = q[d] - (a[d] + tt * (b[d] - a[d]));
            off2 += r * r;
          }
          worst_frame_offset = std::max(worst_frame_offset, std::sqrt(off2));
        }
      }
      INFO("worst frame-reconstruction offset of a first-side stack cell: "
           << worst_frame_offset);
      CHECK(worst_frame_offset <= 1.0e-9);
      CHECK(half.domain_boundary_ft > 0.0);
      if (const char *debug_dir = std::getenv("PALACE_SYMMETRY_DEBUG_DIR"))
      {
        std::ofstream record(fs::path(debug_dir) / "symmetry-stack.json");
        record << json{{"Case", "stack"},
                       {"StackFeature", stack_id},
                       {"StackRawPortions", stack_raw_portions},
                       {"StackRawLength", stack_raw_length},
                       {"WorstFrameOffset", worst_frame_offset},
                       {"HalfDomainBoundary", half.domain_boundary_ft},
                       {"Ft", half.corrected[0]},
                       {"Reasons", exclusions.at("Reasons")}}
                      .dump(2)
               << "\n";
      }
    }
  }

  // Decision 545 MAJOR-2: the defect of record end to end. A BENT four-edge stack (the C2 /
  // C3 / S1b class: a CPW stack meandering through sub-noise joints, its sides' chords
  // staggered): t(x) with the slopes 0.13 / 0.10 / 0.07 / 0 on the quarter runs (kinks of
  // 0.03 rad at x = 0.25 and 0.5, 0.07 at x = 0.75; one chain per edge, windowed bend
  // radius 17.9 R). With its images at both walls the stack is a bent stack (no mirror
  // placement: Unmerged at both planes), so EVERY cell of the matched stack (x in [0.275,
  // 1]; the steeper piece next to x = 0 reads as Missing bent clusters, raw uncovered) is
  // the configuration's own, DomainBoundary raw. The stack bends TOWARD its first side
  // (slopes decreasing), so a first-side sample just past a joint has its far foot on the
  // LAST side's PREVIOUS chord (the perpendicular foot onto it is nearer than onto the own
  // chord's partner): AxisW follows that chord, skewed from the own segment by the kink.
  // The frame reconstruction of record (origin + c AxisW, EdgeOffset 0 on the first side)
  // then puts the cell's midpoint ON its segment (|c_mid| ~ 0.005 x 0.03 = 1.5e-4 < 1e-3
  // R = 2e-4) and its ends OFF it (0.036 x 0.03 = 1.1e-3 > 2e-4): exactly the acceptance
  // UncoveredEdgeIntervals (surfacepostoperator.cpp) refuses with "lies within the
  // tolerance of a perimeter segment at its midpoint but not at its ends" - the J2 abort.
  // (The first Gauss point after the x = 0.5 joint: the foot's shift 0.375 x 0.10 = 0.037
  // exceeds its 0.026 from the joint; at x = 0.75 the previous slope 0.07 shifts the foot
  // by 0.026, clamped at the joint, no stagger.) The solve passes with the pre-image
  // claims; the OLD claims (the placement's frames from the dry run of the SAME case, the
  // formula of record) are then handed to the real postprocessor
  // (SurfacePostOperator::GetInterfaceUncoveredEdgeEnergies on the fixture mesh) and
  // refused with that message; the new claims (the solve's RawPortions) are accepted.
  SECTION("S-bent-stack (decision 545): a bent four-edge stack oblique at the plane - the "
          "frame claims of record abort the postprocessor, the pre-image claims pass")
  {
    const std::vector<PadRectangle> pads = {
        {0.0, 0.5, 1.0, 1.0}, {0.0, 1.5, 1.0, 1.625}, {0.0, 1.75, 1.0, 1.875}};
    const PadBend bend = {0.13, 0.13, 0.10, 0.10, 0.07, 0.07, 0.0, 0.0};
    const std::string name = "mirror-pad-bent-stack-half";
    const Energies half = Run(name, 0.0, false, "Natural", pads, 16, bend, 2);
    // The claims, read from the case's records on every rank (the postprocessor call below
    // is collective): the dry run's placement (frames and cells) and the solve's manifest
    // (the identification segments, the DomainBoundary patches and raw portions).
    std::vector<SurfacePostOperator::UncoveredPerimeterPortion> frame_claims, record_claims;
    double worst_frame_offset = 0.0, worst_frame_midpoint = 0.0, worst_record_offset = 0.0;
    int stack_raw_portions = 0, skewed_cells = 0;
    {
      std::ifstream metadata_input(temp.temp_dir / name / "palace.json");
      REQUIRE(metadata_input);
      const json metadata = json::parse(metadata_input);
      const json &exclusions =
          metadata.at("SurfaceResponse").at("Diagnostics").at("DomainBoundaryExclusions");
      std::ifstream manifest_input(temp.temp_dir / (name + "-dryrun") /
                                   "requirements.json");
      REQUIRE(manifest_input);
      const json manifest = json::parse(manifest_input);
      const json &segments = manifest.at("Identification").at("Segments");
      const auto rows =
          ReadDryRun(temp.temp_dir / (name + "-dryrun") / "surface-response-patches.csv");
      auto Off = [&](int segment, const std::array<double, 3> &q)
      {
        const auto &key = segments.at(static_cast<std::size_t>(segment)).at("Key");
        std::array<double, 3> a = key[0], b = key[1];
        double ab2 = 0.0, t = 0.0;
        for (int d = 0; d < 3; d++)
        {
          ab2 += (b[d] - a[d]) * (b[d] - a[d]);
          t += (q[d] - a[d]) * (b[d] - a[d]);
        }
        t = std::clamp(t / ab2, 0.0, 1.0);
        double off2 = 0.0;
        for (int d = 0; d < 3; d++)
        {
          const double r = q[d] - (a[d] + t * (b[d] - a[d]));
          off2 += r * r;
        }
        return std::sqrt(off2);
      };
      for (const auto &portion : exclusions.at("RawPortions").at("Portions"))
      {
        if (portion.at("Type") != "parallel-edge cluster")
        {
          continue;
        }
        const int segment = portion.at("Segment").get<int>();
        REQUIRE(segment >= 0);
        const std::array<double, 3> p0 = portion.at("P0"), p1 = portion.at("P1");
        worst_record_offset =
            std::max({worst_record_offset, Off(segment, p0), Off(segment, p1)});
        record_claims.push_back({p0, p1, "parallel-edge cluster", portion.at("Feature")});
        stack_raw_portions++;
      }
      for (const auto &patch : exclusions.at("Patches"))
      {
        if (patch.at("Topology") != "parallel-edge cluster")
        {
          continue;
        }
        const std::size_t index = patch.at("Patch").get<std::size_t>();
        const auto row = std::find_if(rows.begin(), rows.end(),
                                      [&](const DryRunRow &r) { return r.patch == index; });
        REQUIRE(row != rows.end());
        const int segment = patch.at("Segment").get<int>();
        if (Off(segment, row->origin) > 1.0e-9)
        {
          continue;  // another side (EdgeOffset != 0): its frame claim needs the offset
        }
        std::array<double, 3> p0{}, p1{}, midpoint{};
        for (int d = 0; d < 3; d++)
        {
          p0[d] = row->origin[d] + row->strip[0] * row->axis_w[d];
          p1[d] = row->origin[d] + row->strip[1] * row->axis_w[d];
          midpoint[d] = 0.5 * (p0[d] + p1[d]);
        }
        const double off = std::max(Off(segment, p0), Off(segment, p1));
        worst_frame_offset = std::max(worst_frame_offset, off);
        if (off > 1.0e-3 * R)
        {
          skewed_cells++;
          worst_frame_midpoint = std::max(worst_frame_midpoint, Off(segment, midpoint));
        }
        frame_claims.push_back({p0, p1, "parallel-edge cluster", patch.at("Feature")});
      }
    }
    INFO("stack raw portions " << stack_raw_portions << ", first-side frame claims "
                               << frame_claims.size() << " of which skewed beyond 1e-3 R: "
                               << skewed_cells << " (worst ends " << worst_frame_offset
                               << ", their midpoints at most " << worst_frame_midpoint
                               << "); record claims at most " << worst_record_offset
                               << " off their segments");
    // The solve passed the postprocessor's acceptance with the pre-image claims (the
    // identity holds on its raw-kept DomainBoundary term); every record claim on its
    // segment; the frame claims of the first side skewed beyond the tolerance with their
    // midpoints inside it (the J2 configuration).
    if (Mpi::Root(Mpi::World()))
    {
      CheckIdentity(half);
      CHECK(half.domain_boundary_ft > 0.0);
    }
    CHECK(stack_raw_portions > 0);
    CHECK(worst_record_offset <= 1.0e-9);
    CHECK(skewed_cells > 0);
    CHECK(worst_frame_offset > 1.0e-3 * R);
    CHECK(worst_frame_midpoint <= 1.0e-3 * R);
    // The real postprocessor on the fixture mesh (the pad's interface 4 with its edge
    // distance R, a smooth nonzero field): the frame claims of record are REFUSED by
    // UncoveredEdgeIntervals, the pre-image claims accepted.
    {
      mfem::Mesh serial = MakeSerialMirrorPadMesh(0.0, false, pads, 16, bend);
      Mesh mesh(std::make_unique<mfem::ParMesh>(Mpi::World(), serial));
      mfem::H1_FECollection h1(1, 3);
      mfem::ND_FECollection nd(1, 3);
      FiniteElementSpace h1_space(mesh, &h1), nd_space(mesh, &nd);
      config::MaterialData material;
      material.attributes = {1};
      config::PeriodicBoundaryData periodic;
      MaterialOperator materials({material}, periodic, ProblemType::ELECTROSTATIC, mesh);
      GridFunction field(nd_space, false);
      mfem::VectorFunctionCoefficient coefficient(3,
                                                  [](const mfem::Vector &x, mfem::Vector &E)
                                                  {
                                                    E.SetSize(3);
                                                    E = 0.0;
                                                    E(1) = 1.0 + 0.5 * x(0) + 0.25 * x(2);
                                                  });
      field.Real().ProjectCoefficient(coefficient);
      config::InterfaceDielectricData data;
      data.attributes = {9};
      data.type = InterfaceDielectric::MA;
      data.t = 0.002;
      data.epsilon_r = 4.0;
      data.edge_attributes = {9};
      data.edge_distances = {R};
      config::BoundaryPostData postpro;
      postpro.dielectric.emplace(4, data);
      SurfacePostOperator surf_post_op(postpro, ProblemType::ELECTROSTATIC, materials,
                                       h1_space, nd_space);
      CHECK_THROWS_WITH(
          surf_post_op.GetInterfaceUncoveredEdgeEnergies({4}, field, nullptr, frame_claims),
          Catch::Matchers::ContainsSubstring(
              "lies within the tolerance of a perimeter segment at its midpoint but not "
              "at its ends"));
      const auto accepted = surf_post_op.GetInterfaceUncoveredEdgeEnergies(
          {4}, field, nullptr, record_claims);
      REQUIRE(accepted.size() == 1);
      CHECK(accepted.at(4).energy > 0.0);
      CHECK(accepted.at(4).by_type.count("parallel-edge cluster") == 1);
    }
    if (const char *debug_dir = std::getenv("PALACE_SYMMETRY_DEBUG_DIR");
        debug_dir && Mpi::Root(Mpi::World()))
    {
      std::ofstream record(fs::path(debug_dir) / "symmetry-bent-stack.json");
      record << json{{"Case", "bent-stack"},
                     {"StackRawPortions", stack_raw_portions},
                     {"FirstSideFrameClaims", frame_claims.size()},
                     {"SkewedCells", skewed_cells},
                     {"WorstFrameOffset", worst_frame_offset},
                     {"WorstFrameMidpoint", worst_frame_midpoint},
                     {"WorstRecordOffset", worst_record_offset},
                     {"Tolerance", 1.0e-3 * R},
                     {"HalfDomainBoundary", half.domain_boundary_ft},
                     {"Ft", half.corrected[0]},
                     {"Reasons",
                      half.diagnostics.at("DomainBoundaryExclusions").at("Reasons")}}
                    .dump(2)
             << "\n";
    }
  }

  SECTION("S-off (control): Mirror Off keeps the raw claims and the identity, not the half")
  {
    const Energies full = Run("mirror-pad-45-full-off", 1.0, true, "Natural");
    const Energies half = Run("mirror-pad-45-half-off", 1.0, false, "Off");
    // The checks run on the root; every rank stays in the section (a return inside a
    // SECTION would end the test body early on that rank and desynchronise Catch2's
    // section discovery across the ranks: the next pass then hangs in a collective).
    if (Mpi::Root(Mpi::World()))
    {
      CheckIdentity(half);
      const auto &exclusions = half.diagnostics.at("DomainBoundaryExclusions");
      CHECK(exclusions.at("Count").get<int>() > 0);
      CHECK(exclusions.at("Mirrored").at("Count").get<int>() == 0);
      CHECK(half.domain_boundary_ft > 0.0);
      CHECK(half.diagnostics.at("MirrorBand").at("MirrorFormedFeatures").empty());
      // The DomainBoundary term is the raw energy of the cut cells; the half's ft differs
      // from half the full's by the raw-vs-model difference of those cells (not zero).
      CHECK_THAT(2.0 * half.raw, WithinRel(full.raw, 1.0e-9));
      CHECK(std::abs(2.0 * half.corrected[0] - full.corrected[0]) >
            1.0e-6 * full.corrected[0]);
    }
  }
#endif
}

// Decision 473 MINOR-9: the truncation planes' constants are the quantised grid's (the key
// every rank groups the faces on), so they are identical for any rank count and partition;
// on the sheared pad box (s = 1) they equal the analytic planes to the grid (1e-9 x the
// coordinate scale 2): x - z = 0 (attribute 3, Natural), -x + z = 1 (5, Natural), z = 0 (1)
// and z = 2 (6) (vertical to the process normal y: Natural), y = 0 (2) Unsupported. Run at
// 1 and 2 ranks ([Serial][Parallel]): the same assertions on the same constants.
// Decisions 536 / 537 (DESIGN ERRATA-8): the F-DB-a raw claim of a DomainBoundary
// translational cell is its recorded pre-image on its own segment, not the frame
// reconstruction origin + EdgeOffset AxisU + cell AxisW. A pair's / stack's frame follows
// the PARTNER side (AxisU to the sample's foot on the partner, AxisW = AxisU x AxisV):
// where the foot is clamped at the partner's end or the chords differ, AxisW is skewed off
// the own segment and the frame ends leave it by cell x sin(angle) - the rerun-2 J2 abort
// (8-75 nm against the postprocessor's 1e-3 R tolerance on C2p / C2 / S1b / C3). Reproduced
// here on a synthetic stack cell: the frame reconstruction (a legacy patch without the
// record) fails the postprocessor's acceptance, the pre-image lies on the segment exactly.
TEST_CASE("Domain-boundary raw claim of a pair / stack cell lies on its own segment",
          "[surfaceresponseoperator][domainboundary][stacks][Serial]")
{
  using config::ElectrostaticSolverData;
  constexpr double R = 1.9;                 // um
  constexpr double tolerance = 1.0e-3 * R;  // the postprocessor's acceptance
  // The own segment: along x at z = 0 (the metal plane), the sample at s = 5.0 on it.
  const std::array<double, 3> p0 = {100.0, 20.0, 0.0}, tangent = {1.0, 0.0, 0.0};
  const double sample_s = 5.0;
  // The stack frame of a second-side cell: AxisU points from the first side's foot
  // (EdgeOffset -3.9 away along it) to the far side's foot, clamped 0.4 along the edge at
  // the far side's end next to an oblique plane: AxisU = (0.4, 3.9 x 2, 0) normalised,
  // AxisV the process normal, AxisW = AxisU x AxisV (skewed off the tangent by atan(0.4
  // / 7.8)).
  const double skew = std::atan2(0.4, 7.8);
  const std::array<double, 3> axis_u = {std::sin(skew), std::cos(skew), 0.0};
  const std::array<double, 3> axis_v = {0.0, 0.0, 1.0};
  const std::array<double, 3> axis_w = {axis_u[1] * axis_v[2] - axis_u[2] * axis_v[1],
                                        axis_u[2] * axis_v[0] - axis_u[0] * axis_v[2],
                                        axis_u[0] * axis_v[1] - axis_u[1] * axis_v[0]};
  const double edge_offset = 3.9;  // the own edge from the origin along AxisU
  const std::array<double, 2> cell = {-0.9, 0.9};  // the cell about the sample along AxisW
  ElectrostaticSolverData::ResponseCorrectionData config;
  ElectrostaticSolverData::ResponseCorrectionModelData model;
  model.idx = 1;
  model.topology = "parallel-edge cluster";
  config.models.push_back(model);
  ElectrostaticSolverData::ResponseCorrectionPatchData patch;
  patch.model = 1;
  // The origin: the sample's foot on the first side, EdgeOffset away from the sample.
  const std::array<double, 3> sample = {p0[0] + sample_s * tangent[0],
                                        p0[1] + sample_s * tangent[1],
                                        p0[2] + sample_s * tangent[2]};
  for (int d = 0; d < 3; d++)
  {
    patch.origin[d] = sample[d] - edge_offset * axis_u[d];
  }
  patch.axis_u = axis_u;
  patch.axis_v = axis_v;
  patch.axis_w = axis_w;
  patch.weight = 1.0;
  patch.longitudinal_cell = cell;
  patch.provenance.feature = 2;
  patch.provenance.segment = 7;
  patch.provenance.s0 = 0.0;
  patch.provenance.s1 = 12.0;
  patch.provenance.edge_offset = edge_offset;
  patch.provenance.coupon_depth = 1.9;
  const double projection = tangent[0] * axis_w[0] + tangent[1] * axis_w[1];
  // The pre-image of the cell on the own segment: the cell offsets along AxisW are the
  // arc offsets projected by Dot(tangent, AxisW) (LongitudinalCellOffsets).
  for (int k = 0; k < 2; k++)
  {
    const double s = sample_s + cell[k] / projection;
    for (int d = 0; d < 3; d++)
    {
      patch.provenance.own_cell[k][d] = p0[d] + s * tangent[d];
    }
  }
  auto Off = [&](const std::array<double, 3> &q)
  {
    double t = 0.0;
    for (int d = 0; d < 3; d++)
    {
      t += (q[d] - p0[d]) * tangent[d];
    }
    double off2 = 0.0;
    for (int d = 0; d < 3; d++)
    {
      const double r = q[d] - (p0[d] + t * tangent[d]);
      off2 += r * r;
    }
    return std::sqrt(off2);
  };
  DomainBoundaryExclusions exclusions;
  DomainBoundaryExclusion exclusion;
  exclusion.patch = 0;
  exclusion.reason = "UnmergedTopology";
  exclusions.patches.push_back(exclusion);
  ContinuationOwnership ownership;
  UncoveredSpatialSupportClipping clipping;

  SECTION("the frame reconstruction (a legacy patch without the record) leaves the edge")
  {
    patch.provenance.has_own_cell = false;
    const auto portions =
        CollectDomainBoundaryPortions(exclusions, {patch}, config, ownership, clipping, R);
    REQUIRE(portions.portions.size() == 1);
    const auto &portion = portions.portions.front();
    CHECK(portion.feature == 2);
    CHECK(portion.topology == "parallel-edge cluster");
    // The midpoint is on the edge (the own point), the ends 0.9 x sin(skew) = 46 nm off it:
    // the postprocessor's UncoveredEdgeIntervals accepts the midpoint and refuses the ends
    // ("lies within the tolerance of a perimeter segment at its midpoint but not at its
    // ends") - the J2 abort.
    std::array<double, 3> midpoint{};
    for (int d = 0; d < 3; d++)
    {
      midpoint[d] = 0.5 * (portion.p0[d] + portion.p1[d]);
    }
    CHECK(Off(midpoint) <= tolerance);
    CHECK(Off(portion.p0) > tolerance);
    CHECK(Off(portion.p1) > tolerance);
    CHECK_THAT(Off(portion.p0), WithinRel(0.9 * std::sin(skew), 1.0e-9));
  }

  SECTION("the recorded pre-image lies on the own segment: the whole cell, once")
  {
    patch.provenance.has_own_cell = true;
    const auto portions =
        CollectDomainBoundaryPortions(exclusions, {patch}, config, ownership, clipping, R);
    REQUIRE(portions.portions.size() == 1);
    const auto &portion = portions.portions.front();
    CHECK(portion.feature == 2);
    CHECK(portion.segment == 7);
    CHECK(Off(portion.p0) <= 1.0e-12);
    CHECK(Off(portion.p1) <= 1.0e-12);
    // The interval is the cell's arc length on the segment (the frame cell projected by
    // Dot(tangent, AxisW)), about the sample.
    double length = 0.0;
    for (int d = 0; d < 3; d++)
    {
      length += (portion.p1[d] - portion.p0[d]) * (portion.p1[d] - portion.p0[d]);
    }
    CHECK_THAT(std::sqrt(length), WithinRel((cell[1] - cell[0]) / projection, 1.0e-12));
    CHECK_THAT(portions.length, WithinRel((cell[1] - cell[0]) / projection, 1.0e-12));
    CHECK(portions.geometric_cells == 1);
    // A co-located first-order split patch of the same cell maps to the same interval,
    // counted once.
    auto split = patch;
    split.weight = 0.3;
    exclusions.patches.push_back(exclusion);
    exclusions.patches.back().patch = 1;
    const auto both = CollectDomainBoundaryPortions(exclusions, {patch, split}, config,
                                                    ownership, clipping, R);
    CHECK(both.portions.size() == 1);
    CHECK(both.duplicate_patches == 1);
    CHECK_THAT(both.length, WithinRel(portions.length, 1.0e-12));
  }

  SECTION(
      "the part of a cell attributed to a spatial coupon maps like the cell (decision 545 "
      "MAJOR-1): OwnEdgePointAt on the record vs the frame")
  {
    // ApplyContinuationOwnership attributes the owned part [removed_lo, removed_hi] of a
    // cell to its owning spatial coupon through OwnEdgePointAt; a DomainBoundary coupon
    // keeps those parts raw (CollectDomainBoundaryPortions). Without the record the frame
    // reconstruction leaves the segment by the skew; with it the points lie on the segment
    // and agree with the clip of the same offsets.
    const double removed_lo = cell[0], removed_hi = cell[0] + 0.5;
    patch.provenance.has_own_cell = false;
    const auto frame_lo = OwnEdgePointAt(patch, removed_lo);
    const auto frame_hi = OwnEdgePointAt(patch, removed_hi);
    CHECK(Off(frame_lo) > tolerance);
    CHECK_THAT(Off(frame_lo), WithinRel(0.9 * std::sin(skew), 1.0e-9));
    CHECK(Off(frame_hi) > tolerance);
    patch.provenance.has_own_cell = true;
    const auto record_lo = OwnEdgePointAt(patch, removed_lo);
    const auto record_hi = OwnEdgePointAt(patch, removed_hi);
    CHECK(Off(record_lo) <= 1.0e-12);
    CHECK(Off(record_hi) <= 1.0e-12);
    // The attributed part's arc length is the offsets' length over the projection.
    double length2 = 0.0;
    for (int d = 0; d < 3; d++)
    {
      length2 += (record_hi[d] - record_lo[d]) * (record_hi[d] - record_lo[d]);
    }
    CHECK_THAT(std::sqrt(length2), WithinRel(0.5 / projection, 1.0e-12));
    // The same affine rule as the clip: clipping the cell to [removed_lo, removed_hi]
    // records exactly these two points.
    auto clipped = patch;
    ClipLongitudinalCell(clipped, removed_lo, removed_hi);
    for (int d = 0; d < 3; d++)
    {
      CHECK_THAT(clipped.provenance.own_cell[0][d], WithinAbs(record_lo[d], 1.0e-12));
      CHECK_THAT(clipped.provenance.own_cell[1][d], WithinAbs(record_hi[d], 1.0e-12));
    }
  }

  SECTION("a cell whose segment runs against AxisW: the record follows the cell's offsets")
  {
    // The placement orders own_cell as the cell offsets (own_cell[k] = the pre-image of
    // longitudinal_cell[k]); with the segment running against AxisW (projection < 0) the
    // far-arc end is own_cell[0]. A clip of the cell's LOW offsets must then move that
    // end, not the other (the first gate run's 78-nm shifts on trimmed reversed cells).
    patch.provenance.has_own_cell = true;
    for (int d = 0; d < 3; d++)
    {
      patch.axis_w[d] = -patch.axis_w[d];
      patch.axis_u[d] = -patch.axis_u[d];  // keep AxisW = AxisU x AxisV
    }
    const double reversed_projection = -projection;
    for (int k = 0; k < 2; k++)
    {
      const double s = sample_s + cell[k] / reversed_projection;
      for (int d = 0; d < 3; d++)
      {
        patch.provenance.own_cell[k][d] = p0[d] + s * tangent[d];
      }
    }
    for (int d = 0; d < 3; d++)
    {
      patch.origin[d] = sample[d] + edge_offset * axis_u[d];  // EdgeOffset along -AxisU
    }
    const auto before = patch.provenance.own_cell;
    ClipLongitudinalCell(patch, cell[0] + 0.3, cell[1]);
    const auto portions =
        CollectDomainBoundaryPortions(exclusions, {patch}, config, ownership, clipping, R);
    REQUIRE(portions.portions.size() == 1);
    const auto &portion = portions.portions.front();
    CHECK(Off(portion.p0) <= 1.0e-12);
    CHECK(Off(portion.p1) <= 1.0e-12);
    CHECK_THAT(portions.length, WithinRel((cell[1] - cell[0] - 0.3) / projection, 1.0e-12));
    // own_cell[0] (the low offset's pre-image, the far-arc end here) moved by 0.3 /
    // |projection|; own_cell[1] stayed.
    double moved0 = 0.0, moved1 = 0.0;
    for (int d = 0; d < 3; d++)
    {
      moved0 += (portion.p0[d] - before[0][d]) * (portion.p0[d] - before[0][d]);
      moved1 += (portion.p1[d] - before[1][d]) * (portion.p1[d] - before[1][d]);
    }
    CHECK_THAT(std::sqrt(moved0), WithinRel(0.3 / projection, 1.0e-12));
    CHECK(std::sqrt(moved1) <= 1.0e-12);
  }

  SECTION("a clipped cell (an arm trim, an ownership clip) keeps the kept part's pre-image")
  {
    // The clip keeps [c0 + 0.3, c1] of the cell: the record follows (ClipOwnCell through
    // ClipLongitudinalCell), so the raw claim is the kept part on the segment, not the
    // whole cell nor the frame reconstruction about the moved origin.
    patch.provenance.has_own_cell = true;
    const auto before = patch.provenance.own_cell;
    ClipLongitudinalCell(patch, cell[0] + 0.3, cell[1]);
    CHECK_THAT(patch.longitudinal_cell[1] - patch.longitudinal_cell[0],
               WithinRel(cell[1] - cell[0] - 0.3, 1.0e-12));
    const auto portions =
        CollectDomainBoundaryPortions(exclusions, {patch}, config, ownership, clipping, R);
    REQUIRE(portions.portions.size() == 1);
    const auto &portion = portions.portions.front();
    CHECK(Off(portion.p0) <= 1.0e-12);
    CHECK(Off(portion.p1) <= 1.0e-12);
    // The kept end is the recorded far end; the clipped end moved 0.3 / projection along
    // the segment.
    CHECK_THAT(portions.length, WithinRel((cell[1] - cell[0] - 0.3) / projection, 1.0e-12));
    double far_end_moved = 0.0, near_end_moved = 0.0;
    for (int d = 0; d < 3; d++)
    {
      far_end_moved += (portion.p1[d] - before[1][d]) * (portion.p1[d] - before[1][d]);
      near_end_moved += (portion.p0[d] - before[0][d]) * (portion.p0[d] - before[0][d]);
    }
    CHECK(std::sqrt(far_end_moved) <= 1.0e-12);
    CHECK_THAT(std::sqrt(near_end_moved), WithinRel(0.3 / projection, 1.0e-12));
  }
}

TEST_CASE("Truncation planes are deterministic across rank counts",
          "[surfaceresponsemirror][mirror][mirrorplanes][3d][Serial][Parallel]")
{
  mfem::Mesh serial = MakeSerialMirrorPadMesh(1.0, false);
  mfem::ParMesh mesh(Mpi::World(), serial);
  const auto planes = FitTruncationPlanes(mesh, {1, 2, 3, 5, 6}, {0.0, 1.0, 0.0}, 0.2);
  REQUIRE(planes.size() == 5);
  const double r = 1.0 / std::sqrt(2.0);
  // The offset quantum: 1e-9 x the coordinate scale (the largest |coordinate| of the
  // bounding box: 2 here, at z = 2 and x = 1 + s (z - 1) = 2).
  const double grid = 1.0e-9 * 2.0;
  std::map<int, const MirrorPlane *> by_attribute;
  for (const auto &plane : planes)
  {
    by_attribute.emplace(plane.attribute, &plane);
  }
  auto Expect = [&](int attribute, std::array<double, 3> normal, double offset,
                    const std::string &status)
  {
    REQUIRE(by_attribute.count(attribute) == 1);
    const auto &plane = *by_attribute.at(attribute);
    INFO("attribute " << attribute);
    CHECK(plane.status == status);
    for (int d = 0; d < 3; d++)
    {
      CHECK_THAT(plane.normal[d], WithinAbs(normal[d], 1.0e-9));
    }
    CHECK_THAT(plane.offset, WithinAbs(offset, grid));
    CHECK((plane.faces == 64 || plane.faces == 128));
  };
  Expect(3, {r, 0.0, -r}, 0.0, "Natural");
  Expect(5, {-r, 0.0, r}, r, "Natural");
  Expect(1, {0.0, 0.0, -1.0}, 0.0, "Natural");
  Expect(6, {0.0, 0.0, 1.0}, 2.0, "Natural");
  Expect(2, {0.0, -1.0, 0.0}, 0.0, "Unsupported");
  // The constants are exact multiples of the grid: a point ON the analytic plane is at
  // most half a grid step from the fitted one, the same on every rank.
  for (const auto &plane : planes)
  {
    CHECK_THAT(plane.offset, WithinAbs(std::round(plane.offset / grid) * grid, 1.0e-15));
  }
}

}  // namespace palace
