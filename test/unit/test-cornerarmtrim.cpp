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
#include "driver.hpp"
#include "fem/fespace.hpp"
#include "fem/gridfunction.hpp"
#include "fem/integrator.hpp"
#include "fem/mesh.hpp"
#include "models/laplaceoperator.hpp"
#include "models/materialoperator.hpp"
#include "models/surfacepostoperator.hpp"
#include "models/surfaceresponseidentification.hpp"
#include "models/surfaceresponseoperator.hpp"
#include "utils/communication.hpp"
#include "utils/iodata.hpp"
#include "utils/omp.hpp"
#include "utils/outputdir.hpp"
#include "utils/tablecsv.hpp"
#include "utils/timer.hpp"

namespace palace
{

namespace fs = std::filesystem;
using json = nlohmann::json;
using namespace Catch::Matchers;
using test::ColumnByHeader;
using test::LoadCsv;
using test::ReadFile;
using test::ReadTextCsv;
using test::RunElectrostatic;
using test::TextCsv;

namespace
{

// A PEC lead (attribute 9) on the plane y = 0.5 of the unit cube meshed n x n x n (h = 1
// / n): x in [x_min, x_max] from the bottom face z = 0 (where the lead is cut by the
// domain) to its cap at z = z_cap (0.5 unless stated). With x_max < 1 the cap ends in two
// 90-degree convex corners; with x_max = 1 the lead reaches the wall x = 1 (a truncation
// cut) and has one corner at (x_min, 0.5, z_cap). The mesh is then sheared IN THE METAL
// PLANE, x -> x + shear (z - 0.5): the cap (along x) is unchanged while the long edges
// (along z) tilt to the direction (shear, 0, 1), so the corner angle between the arms
// becomes 90 + atan(shear) degrees (shear = 1 / sqrt(3): 120 degrees) and the lead stays
// planar on y = 0.5. With fillet > 0 (a multiple of h; unsheared) both cap corners are
// ROUNDED as the island fixture rounds its corners (SurfaceResponseFiles::MakeIslandMesh):
// the lead-edge vertices within the fillet square of each corner are projected onto the
// quarter circle of radius fillet tangent to both arms (fillet = 2 h: tangent points, 22.5
// / 45 / 67.5-degree joints, a 4-chord polyline whose vertices lie on the arc), so the
// identification fits one rounded corner of radius fillet per cap corner.
mfem::Mesh MakeSerialLeadMesh(int n, double x_min, double x_max, double shear,
                              double fillet, double z_cap = 0.5)
{
  const double h = 1.0 / n;
  auto Vertex = [n](int i, int j, int k) { return i + (n + 1) * (j + (n + 1) * k); };
  mfem::Mesh serial(3, (n + 1) * (n + 1) * (n + 1), n * n * n, 0, 3);
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
  // z).
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
  // The lead: the interior faces on y = 0.5 with x in [x_min, x_max], z in [0, 0.5].
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
    double xmin = 1.0, xmax = 0.0, zmax = 0.0;
    for (const int vertex : vertices)
    {
      const double *point = serial.GetVertex(vertex);
      on_plane = on_plane && std::abs(point[1] - 0.5) < 1.0e-12;
      xmin = std::min(xmin, point[0]);
      xmax = std::max(xmax, point[0]);
      zmax = std::max(zmax, point[2]);
    }
    if (on_plane && xmin >= x_min - 1.0e-12 && xmax <= x_max + 1.0e-12 &&
        zmax <= z_cap + 1.0e-12)
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
    point[0] += shear * (point[2] - 0.5);
  }
  if (fillet > 0.0)
  {
    constexpr double tolerance = 1.0e-12;
    constexpr double pi = 3.14159265358979323846;
    for (int vertex = 0; vertex < serial.GetNV(); vertex++)
    {
      double *point = serial.GetVertex(vertex);
      if (std::abs(point[1] - 0.5) > tolerance)
      {
        continue;
      }
      for (const double sign_x : {-1.0, 1.0})
      {
        // The cap corner at (corner_x, 0.5, 0.5); its fillet centre fillet inside along
        // both arms; a vertex on the cap line (z = 0.5) or on the long edge (x = corner_x)
        // within the fillet square moves onto the arc at the angle of its arm position.
        const double corner_x = sign_x < 0.0 ? x_min : x_max, corner_z = z_cap;
        const double center_x = corner_x - sign_x * fillet, center_z = corner_z - fillet;
        const double local_x = sign_x * (point[0] - center_x);
        const double local_z = point[2] - center_z;
        if (local_x < -tolerance || local_x > fillet + tolerance || local_z < -tolerance ||
            local_z > fillet + tolerance)
        {
          continue;
        }
        double angle;
        if (std::abs(point[2] - corner_z) <= tolerance)
        {
          angle = 0.5 * pi - 0.25 * pi * local_x / fillet;
        }
        else if (std::abs(point[0] - corner_x) <= tolerance)
        {
          angle = 0.25 * pi * local_z / fillet;
        }
        else
        {
          continue;
        }
        point[0] = center_x + sign_x * fillet * std::cos(angle);
        point[2] = center_z + fillet * std::sin(angle);
        break;
      }
    }
  }
  return serial;
}

mfem::Mesh MakeSerialLeadMesh(double shear, double x_max)
{
  return MakeSerialLeadMesh(8, 0.25, x_max, shear, 0.0);
}

std::unique_ptr<mfem::ParMesh> MakeLeadMesh(double shear, double x_max)
{
  mfem::Mesh serial = MakeSerialLeadMesh(shear, x_max);
  return std::make_unique<mfem::ParMesh>(Mpi::World(), serial);
}

// The rounded lead of the corner-arm extension tests: 16 x 16 x 16 (h = 0.0625), fillets of
// 2 h = 0.125 (0.625 R at R = 0.2; the rounded_library_3d model's radius) on both cap
// corners, cap from x_min to x_max at z_cap.
std::unique_ptr<mfem::ParMesh> MakeRoundedLeadMesh(double x_min, double x_max,
                                                   double z_cap = 0.5)
{
  mfem::Mesh serial = MakeSerialLeadMesh(16, x_min, x_max, 0.0, 0.125, z_cap);
  return std::make_unique<mfem::ParMesh>(Mpi::World(), serial);
}

// Move the mesh vertex at `from` (within 1e-9) to `to`; exactly one vertex must match.
void MoveVertex(mfem::Mesh &mesh, const std::array<double, 3> &from,
                const std::array<double, 3> &to)
{
  int found = -1;
  for (int vertex = 0; vertex < mesh.GetNV(); vertex++)
  {
    const double *point = mesh.GetVertex(vertex);
    if (std::abs(point[0] - from[0]) < 1.0e-9 && std::abs(point[1] - from[1]) < 1.0e-9 &&
        std::abs(point[2] - from[2]) < 1.0e-9)
    {
      REQUIRE(found < 0);
      found = vertex;
    }
  }
  REQUIRE(found >= 0);
  for (int d = 0; d < 3; d++)
  {
    mesh.GetVertex(found)[d] = to[d];
  }
}

// The sheared ROUNDED lead of decision 520 MINOR-1: the lead x in [0.125, 1] (the right
// edge on the wall), sheared so that the left long edge tilts to (-shear, 0, -1) from the
// cap corner C = (0.125, 0.5, 0.5): a convex corner of angle theta = 90 + atan(shear)
// degrees (120 at shear 1 / sqrt(3), 135 at shear 1), given a TRUE circular fillet of
// radius 0.125 in physical coordinates: the cap vertex at x = C_x + h moves to the tangent
// point T_A = C + t_d a (t_d = r / tan(theta / 2)), the first edge vertex below C to T_B =
// C + t_d b, the corner vertex onto the arc midpoint; a 2-chord polyline on the circle.
std::unique_ptr<mfem::ParMesh> MakeShearedRoundedLeadMesh(double shear)
{
  constexpr double r = 0.125, h = 1.0 / 16, x_min = 0.125;
  mfem::Mesh serial = MakeSerialLeadMesh(16, x_min, 1.0, shear, 0.0);
  const double n = std::sqrt(1.0 + shear * shear);
  const std::array<double, 3> C = {x_min, 0.5, 0.5};
  const std::array<double, 3> a = {1.0, 0.0, 0.0};
  const std::array<double, 3> b = {-shear / n, 0.0, -1.0 / n};
  const double theta = std::acos(-shear / n);  // between a and b
  const double t_d = r / std::tan(0.5 * theta);
  REQUIRE(t_d < 2.0 * h);
  REQUIRE(t_d <= n * h + 1.0e-12);
  const std::array<double, 3> center = {C[0] + t_d, 0.5, 0.5 - r};
  const double mid = 0.5 * 3.14159265358979323846 + 0.5 * (3.14159265358979323846 - theta);
  MoveVertex(serial, {x_min + h, 0.5, 0.5}, {C[0] + t_d * a[0], 0.5, 0.5});
  MoveVertex(serial, {x_min - shear * h, 0.5, 0.5 - h},
             {C[0] + t_d * b[0], 0.5, C[2] + t_d * b[2]});
  MoveVertex(serial, {x_min, 0.5, 0.5},
             {center[0] + r * std::cos(mid), 0.5, center[2] + r * std::sin(mid)});
  return std::make_unique<mfem::ParMesh>(Mpi::World(), serial);
}

// The rounded VIRTUAL corner of decision 520 MAJOR-1: the 16 x 16 x 16 lead x in [0.375, 1]
// (the right edge on the wall, the cap at z = 0.5) sheared uniformly in its metal plane by
// tan(30 deg) (affine hexes; the F1 fixture's construction): the right wall tilts by 30
// degrees to x = 1 + (z - 0.5) tan(30 deg), so that the cap (along x) meets it at 60
// degrees, the left long edge tilts to (-tan 30, 0, -1) (the cap's left corner is a sharp
// 120-degree convex corner, F1's), and that edge meets the bottom cut at 60 degrees (a
// sharp 120-degree virtual corner with its image, placed HalfByMirror: the no-stretch
// control). The cap's right end is a HALF fillet of radius 0.125 (the model's) whose
// virtual corner C = (1, 0.5, 0.5) lies ON the wall and whose arc meets the wall
// perpendicularly at its midpoint: the mirror image completes a 120-degree ROUNDED corner
// (the cap and its image; t_d = r / tan(60 deg)) placed HalfByMirror. Built in physical
// (sheared) coordinates: the cap vertices x = 14 / 16 and 15 / 16 move onto the arc at the
// tangent point T_A = C - t_d x and at 75 degrees, the corner vertex onto the arc midpoint
// (60 degrees; on the wall).
constexpr double kVirtualLeadShear = 0.57735026918962576;  // tan(30 deg)
constexpr double kVirtualLeadLeft = 0.375;

std::unique_ptr<mfem::ParMesh> MakeVirtualRoundedLeadMesh()
{
  constexpr double r = 0.125, pi = 3.14159265358979323846;
  mfem::Mesh serial = MakeSerialLeadMesh(16, kVirtualLeadLeft, 1.0, kVirtualLeadShear, 0.0);
  const double t_d = r / std::tan(60.0 * pi / 180.0);
  const std::array<double, 3> center = {1.0 - t_d, 0.5, 0.5 - r};
  auto OnArc = [&](double degrees)
  {
    const double a = degrees * pi / 180.0;
    return std::array<double, 3>{center[0] + r * std::cos(a), 0.5,
                                 center[2] + r * std::sin(a)};
  };
  MoveVertex(serial, {0.875, 0.5, 0.5}, OnArc(90.0));
  MoveVertex(serial, {0.9375, 0.5, 0.5}, OnArc(75.0));
  MoveVertex(serial, {1.0, 0.5, 0.5}, OnArc(60.0));
  return std::make_unique<mfem::ParMesh>(Mpi::World(), serial);
}

struct DryRunRow
{
  std::size_t patch = 0;
  int feature = -1;
  int segment = -1;
  std::string topology;
  double weight = 0.0, model_weight = 0.0, quadrature_weight = 0.0;
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
                    segment = Column("Segment"), topology = Column("Topology"),
                    weight = Column("Weight"), model_weight = Column("ModelWeight"),
                    quadrature_weight = Column("QuadratureWeight"), s0 = Column("S0"),
                    s1 = Column("S1"), origin = Column("OriginX"),
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
    row.segment = std::stoi(fields[segment]);
    row.topology = fields[topology];
    row.weight = std::stod(fields[weight]);
    row.model_weight = std::stod(fields[model_weight]);
    row.quadrature_weight = std::stod(fields[quadrature_weight]);
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

// Distance from the corner vertex to the cell end of a translational dry-run row nearest
// to it (the own edge of an isolated edge is the cell line itself).
double NearestCellEndDistance(const DryRunRow &row, const std::array<double, 3> &vertex)
{
  double nearest = std::numeric_limits<double>::infinity();
  for (const double offset : row.strip)
  {
    double distance2 = 0.0;
    for (int d = 0; d < 3; d++)
    {
      const double delta = row.origin[d] + offset * row.axis_w[d] - vertex[d];
      distance2 += delta * delta;
    }
    nearest = std::min(nearest, std::sqrt(distance2));
  }
  return nearest;
}

// The sum of quadrature x model weights per (feature, segment, portion) of the isolated
// edges: 1 on an untouched portion, 1 - removed / portion length on a trimmed one.
std::map<std::tuple<int, int, double, double>, double>
PortionQuadratureSums(const std::vector<DryRunRow> &rows)
{
  std::map<std::tuple<int, int, double, double>, double> sums;
  for (const auto &row : rows)
  {
    if (row.segment < 0)
    {
      continue;
    }
    sums[{row.feature, row.segment, std::round(row.s0 * 1.0e9) / 1.0e9,
          std::round(row.s1 * 1.0e9) / 1.0e9}] += row.quadrature_weight * row.model_weight;
  }
  return sums;
}

}  // namespace

// Decision 394 F1: a matched corner's coupon is calibrated on the matching square |u|, |v|
// <= R of its canonical frame (u = the first arm), which contains the second arm up to s =
// R / max(|cos theta|, |sin theta|) from the vertex while the vertex window claims R along
// each arm, so the second arm's straight cells begin at s (their part before s removed)
// and the stretch [R, s) is no longer modelled twice; at 90 degrees s = R and nothing
// changes (the legacy layout). The lead sheared to a 120-degree convex corner (s = 1.1547
// R) against the unsheared 90-degree lead: the dry run, the manifest record and the
// operator record (1 and 2 ranks).
TEST_CASE_METHOD(test::SurfaceResponseFiles, "SurfaceResponseOperator corner-arm trim",
                 "[surfaceresponseoperator][cornerarmtrim][3d][Serial][Parallel]")
{
#if !defined(MFEM_USE_GSLIB)
  SKIP("SurfaceResponseOperator requires MFEM_USE_GSLIB");
#else
  constexpr double R = 0.2;
  const double shear = 1.0 / std::sqrt(3.0);  // 90 + atan(shear) = 120 degrees
  const auto basis_path = temp.temp_dir / "corner-arm-trim-basis-points.csv";
  const auto library_path = temp.temp_dir / "fabrication-process-corner-arm-trim-3d.json";
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
    library["Name"] = "unit-test-process-corner-arm-trim-3d";
    json corner_120;
    for (auto &model : library["Models"])
    {
      model["Interfaces"] = {{{"Type", "SA"}, {"Coupon", 1}}};
      if (model["Topology"] == "IsolatedEdge")
      {
        model["BasisPoints"] = basis_path.string();
      }
      if (model["Topology"] == "ConvexCorner")
      {
        corner_120 = model;
        corner_120["Name"] = "convex-corner-120";
        corner_120["Angle"] = 120.0;
      }
    }
    REQUIRE(!corner_120.is_null());
    library["Models"].push_back(corner_120);
    std::ofstream output(library_path);
    output << library.dump(2) << "\n";
  }
  Mpi::Barrier(Mpi::World());
  json config = IslandConfig();
  config["Boundaries"]["Ground"]["Attributes"] = {4};
  auto &correction = config["Solver"]["Electrostatic"]["ResponseCorrection"];
  correction.erase("PatchConstruction");
  correction["Library"] = library_path.string();
  correction["TraceCoupling"] = "SurfaceMortar";
  correction["MortarOversampling"] = 2;
  // The lead reaches the wall x = 1 obliquely: under the mirror rule (boundary-cut DESIGN
  // 2.2) that meeting is a virtual corner; this test is about the F1 trim of the REAL
  // corner alone, so the mirror is off (its own tests: [mirror], test-domainboundary.cpp).
  correction["DomainBoundary"] = {{"Mirror", "Off"}};
  IoData iodata(config, false);
  iodata.boundaries.cracked_attributes.insert(9);
  const auto manifest_path = temp.temp_dir / "surface-response-requirements-trim.json";
  const auto patches_path = temp.temp_dir / "surface-response-patches.csv";
  auto Preflight = [&](double mesh_shear, double x_max)
  {
    Mesh dry_run_mesh(MakeLeadMesh(mesh_shear, x_max));
    WriteSurfaceResponseRequirements(iodata, dry_run_mesh, manifest_path.string());
    Mpi::Barrier(Mpi::World());
    std::ifstream input(manifest_path);
    REQUIRE(input);
    return json::parse(input);
  };
  const std::array<double, 3> vertex = {0.25, 0.5, 0.5};

  SECTION("a 120-degree corner's second-arm cells begin where the arm exits the square")
  {
    const json manifest = Preflight(shear, 1.0);
    const auto rows = ReadDryRun(patches_path);
    const auto &record = manifest["Identification"]["Diagnostics"]["CornerArmTrim"];
    const double theta = 120.0 * std::acos(-1.0) / 180.0;
    const double exit = R / std::max(std::abs(std::cos(theta)), std::abs(std::sin(theta)));
    REQUIRE(record["Count"].get<int>() == 1);
    const auto &corner = record["Corners"][0];
    CHECK(corner["Topology"] == "ConvexCorner");
    CHECK_THAT(corner["AngleDegrees"].get<double>(), WithinAbs(120.0, 1.0e-6));
    CHECK_THAT(corner["ExitDistanceOverR"].get<double>(),
               WithinAbs(2.0 / std::sqrt(3.0), 1.0e-9));
    CHECK_THAT(corner["TrimmedLength"].get<double>(), WithinAbs(exit - R, 1.0e-9));
    // The second arm (the tilted long edge; the cap along +x is the canonical first arm,
    // the second being counterclockwise about the plane normal) begins at R from the
    // vertex, so exactly the geometric trim is removed from its first cell.
    CHECK_THAT(corner["RemovedCellLength"].get<double>(), WithinAbs(exit - R, 1.0e-9));
    CHECK(corner["RemovedUncoveredLength"].get<double>() == 0.0);
    REQUIRE(corner["Patches"].size() == 1);
    REQUIRE(record["Cells"].size() == 1);
    const std::size_t clipped = corner["Patches"][0].get<std::size_t>();
    CHECK(record["Cells"][0]["Patch"].get<std::size_t>() == clipped);
    CHECK_THAT(record["Cells"][0]["RemovedLength"].get<double>(),
               WithinAbs(exit - R, 1.0e-9));
    CHECK_THAT(manifest["Summary"]["CornerArmTrim"]["TrimmedLength"].get<double>(),
               WithinAbs(exit - R, 1.0e-9));
    CHECK(manifest["Summary"]["CornerArmTrim"]["Corners"].get<int>() == 1);
    // The corner's canonical first arm (Frame.Axes[0]; the second is counterclockwise
    // about the plane normal): its cells begin at R as before (one exactly there); every
    // cell of the second arm lies at least s from the vertex, the clipped one exactly
    // there, with the dry run's weight formula holding on the kept cell.
    std::array<double, 3> first_arm{};
    int corner_features = 0;
    for (const auto &feature : manifest["Identification"]["Features"])
    {
      if (feature["Type"] == "ConvexCorner")
      {
        corner_features++;
        CHECK(feature["Match"]["Status"] == "Matched");
        CHECK(feature["Match"]["Model"] == "convex-corner-120");
        first_arm = feature["Frame"]["Axes"][0].get<std::array<double, 3>>();
      }
    }
    REQUIRE(corner_features == 1);
    int corner_rows = 0, first_arm_cells = 0, second_arm_cells = 0, at_radius = 0;
    for (const auto &row : rows)
    {
      if (row.topology == "convex corner")
      {
        corner_rows++;
        continue;
      }
      REQUIRE(row.topology == "isolated edge");
      const double distance = NearestCellEndDistance(row, vertex);
      const double along_first =
          std::abs(row.axis_w[0] * first_arm[0] + row.axis_w[1] * first_arm[1] +
                   row.axis_w[2] * first_arm[2]);
      if (along_first > 0.9)
      {
        first_arm_cells++;
        CHECK(distance >= R - 1.0e-9);
        at_radius += std::abs(distance - R) < 1.0e-9 ? 1 : 0;
        CHECK(row.patch != clipped);
      }
      else
      {
        second_arm_cells++;
        CHECK(distance >= exit - 1.0e-9);
        if (row.patch == clipped)
        {
          CHECK_THAT(distance, WithinAbs(exit, 1.0e-9));
          CHECK_THAT(row.quadrature_weight * (row.s1 - row.s0),
                     WithinAbs(row.strip[1] - row.strip[0], 1.0e-9));
        }
      }
    }
    CHECK(corner_rows == 1);
    CHECK(first_arm_cells > 0);
    CHECK(second_arm_cells > 0);
    CHECK(at_radius == 1);
    // The portion quadrature sums: 1 everywhere except the trimmed portion, which reads
    // 1 - removed / portion length (the A7 audit's rule with the trim record).
    const auto sums = PortionQuadratureSums(rows);
    int trimmed_portions = 0;
    for (const auto &[key, sum] : sums)
    {
      const auto &[feature, segment, s0, s1] = key;
      if (std::abs(sum - 1.0) > 1.0e-9)
      {
        trimmed_portions++;
        // S0 / S1 are written on the manifest's length grid (1e-10 R).
        CHECK_THAT(sum, WithinAbs(1.0 - (exit - R) / (s1 - s0), 1.0e-7));
      }
    }
    CHECK(trimmed_portions == 1);
    // The operator applies the same trimmed cells and carries the same record.
    std::vector<std::unique_ptr<Mesh>> meshes;
    meshes.push_back(std::make_unique<Mesh>(MakeLeadMesh(shear, 1.0)));
    LaplaceOperator laplace(iodata, meshes);
    SurfaceResponseOperator response(iodata, laplace);
    const auto statistics = response.GetStatistics();
    const auto &operator_record = statistics["Diagnostics"]["CornerArmTrim"];
    REQUIRE(operator_record["Count"].get<int>() == 1);
    CHECK_THAT(operator_record["TrimmedLength"].get<double>(), WithinAbs(exit - R, 1.0e-9));
    CHECK_THAT(operator_record["RemovedCellLength"].get<double>(),
               WithinAbs(exit - R, 1.0e-9));
    CHECK(operator_record["Cells"][0]["Patch"].get<std::size_t>() == clipped);
    CHECK(statistics["Diagnostics"]["Uncovered"]["Count"].get<int>() == 0);
  }

  SECTION("a 90-degree corner keeps the legacy layout")
  {
    const json manifest = Preflight(0.0, 0.75);
    const auto rows = ReadDryRun(patches_path);
    const auto &record = manifest["Identification"]["Diagnostics"]["CornerArmTrim"];
    CHECK(record["Count"].get<int>() == 0);
    // A sharp corner's claim ends at R = the square exit: no extension either (decision
    // 511; the golden byte identity below covers the patches).
    CHECK(manifest["Identification"]["Diagnostics"]["CornerArmExtension"]["Count"]
              .get<int>() == 0);
    CHECK(manifest["Summary"]["CornerArmExtension"]["Corners"].get<int>() == 0);
    CHECK(record["TrimmedLength"].get<double>() == 0.0);
    CHECK(record["RemovedCellLength"].get<double>() == 0.0);
    CHECK(record["Corners"].empty());
    CHECK(record["Cells"].empty());
    CHECK(manifest["Summary"]["CornerArmTrim"]["Corners"].get<int>() == 0);
    CHECK(manifest["Summary"]["Uncovered"]["Features"].get<int>() == 0);
    // Both corners' arm cells begin at exactly R and every portion is tiled by its cells
    // with unit quadrature sum (nothing clipped).
    const std::array<double, 3> other_vertex = {0.75, 0.5, 0.5};
    int corner_rows = 0, at_radius = 0;
    for (const auto &row : rows)
    {
      if (row.topology == "convex corner")
      {
        corner_rows++;
        continue;
      }
      const double distance = std::min(NearestCellEndDistance(row, vertex),
                                       NearestCellEndDistance(row, other_vertex));
      CHECK(distance >= R - 1.0e-9);
      at_radius += std::abs(distance - R) < 1.0e-9 ? 1 : 0;
    }
    CHECK(corner_rows == 2);
    CHECK(at_radius == 4);  // the two arms of each corner
    for (const auto &[key, sum] : PortionQuadratureSums(rows))
    {
      CHECK_THAT(sum, WithinAbs(1.0, 1.0e-12));
    }
    CHECK(rows.size() == 18);
    // Byte identity with the trim-disabled dry run (decision 399 MINOR-6): the golden
    // test/data/surfaceresponse/corner-arm-trim-90deg-trim-disabled-patches.csv was written
    // by this very dry run with the ApplyCornerArmTrim call removed (serial and 2 ranks
    // identical; its README records the provenance), so a 90-degree corner's patches are
    // the legacy layout to the byte, not only in the asserted columns.
    if (Mpi::Root(Mpi::World()))
    {
      const fs::path golden =
          fs::path(PALACE_TEST_DATA_DIR) /
          "surfaceresponse/corner-arm-trim-90deg-trim-disabled-patches.csv";
      REQUIRE(fs::is_regular_file(golden));
      CHECK(ReadFile(patches_path) == ReadFile(golden));
    }
  }
#endif
}

// Decision 511 O2 (i) (fillet-basis design 2026-10-07 section 7.1): a matched ROUNDED
// corner's identification claim runs R along each arm from the TANGENT point, i.e. to t_d +
// R from the virtual corner (t_d = r / tan(theta / 2); 90 deg: t_d = r), while its coupon
// is calibrated on the square |u|, |v| <= R about the virtual corner, so the stretch [R,
// t_d
// + R) of each arm was claimed by the corner (no straight cell, not uncovered) yet outside
// the coupon: its within-R energy was modelled by nothing. The F1 rule completed in the
// other direction: the arm's straight cells begin where the arm exits the square, the cell
// beginning at the claim end extended back to R. The rounded lead (fillets r = 0.125 =
// 0.625 R on both cap corners, matched to the rounded_library_3d model of that radius):
// the dry run, the manifest record and the operator record; the short-cap lead for the
// stretch without a host (the neighbouring claim abuts: recorded as unhosted, not silent).
TEST_CASE_METHOD(test::SurfaceResponseFiles, "SurfaceResponseOperator corner-arm extension",
                 "[surfaceresponseoperator][cornerarmtrim][3d][Serial][Parallel]")
{
#if !defined(MFEM_USE_GSLIB)
  SKIP("SurfaceResponseOperator requires MFEM_USE_GSLIB");
#else
  constexpr double R = 0.2, r = 0.125;
  const auto basis_path = temp.temp_dir / "corner-arm-extension-basis-points.csv";
  const auto library_path =
      temp.temp_dir / "fabrication-process-corner-arm-extension-3d.json";
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
    std::ifstream input(rounded_library_3d_path);
    REQUIRE(input);
    json library = json::parse(input);
    library["Name"] = "unit-test-process-corner-arm-extension-3d";
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
  config["Boundaries"]["Ground"]["Attributes"] = {4};
  auto &correction = config["Solver"]["Electrostatic"]["ResponseCorrection"];
  correction.erase("PatchConstruction");
  correction["Library"] = library_path.string();
  correction["TraceCoupling"] = "SurfaceMortar";
  correction["MortarOversampling"] = 2;
  correction["DomainBoundary"] = {{"Mirror", "Off"}};
  IoData iodata(config, false);
  iodata.boundaries.cracked_attributes.insert(9);
  const auto manifest_path = temp.temp_dir / "surface-response-requirements-extension.json";
  const auto patches_path = temp.temp_dir / "surface-response-patches.csv";
  auto Preflight = [&](double x_min, double x_max, double z_cap = 0.5)
  {
    Mesh dry_run_mesh(MakeRoundedLeadMesh(x_min, x_max, z_cap));
    WriteSurfaceResponseRequirements(iodata, dry_run_mesh, manifest_path.string());
    Mpi::Barrier(Mpi::World());
    std::ifstream input(manifest_path);
    REQUIRE(input);
    return json::parse(input);
  };
  auto Distance = [](const std::array<double, 3> &a, const std::array<double, 3> &b)
  {
    double d2 = 0.0;
    for (int d = 0; d < 3; d++)
    {
      d2 += (a[d] - b[d]) * (a[d] - b[d]);
    }
    return std::sqrt(d2);
  };

  SECTION("both arms of each rounded corner extended back to the square exit")
  {
    // Virtual corners at (0.125, 0.5, 0.5) and (0.875, 0.5, 0.5); tangent points t_d = r
    // along each arm; claims to t_d + R = 0.325 from the virtual corner; square exit R.
    const double x_min = 0.125, x_max = 0.875;
    const json manifest = Preflight(x_min, x_max);
    const auto rows = ReadDryRun(patches_path);
    const std::array<std::array<double, 3>, 2> virtual_corners = {
        std::array<double, 3>{x_min, 0.5, 0.5}, std::array<double, 3>{x_max, 0.5, 0.5}};
    int corner_features = 0;
    for (const auto &feature : manifest["Identification"]["Features"])
    {
      if (feature["Type"] == "ConvexCorner")
      {
        corner_features++;
        CHECK(feature["Match"]["Status"] == "Matched");
        CHECK(feature["Match"]["Model"] == "convex-corner-90-r0.125");
        CHECK_THAT(feature["Signature"]["CornerRadiusOverR"].get<double>(),
                   WithinAbs(r / R, 1.0e-9));
        const auto origin = feature["Frame"]["Origin"].get<std::array<double, 3>>();
        CHECK(std::min(Distance(origin, virtual_corners[0]),
                       Distance(origin, virtual_corners[1])) < 1.0e-6);
      }
    }
    REQUIRE(corner_features == 2);
    const auto &record = manifest["Identification"]["Diagnostics"]["CornerArmExtension"];
    REQUIRE(record["Count"].get<int>() == 2);
    CHECK(record["ArmsExtended"].get<int>() == 4);
    CHECK(record["ArmsUnhosted"].get<int>() == 0);
    CHECK_THAT(record["StretchLength"].get<double>(), WithinAbs(4.0 * r, 1.0e-9));
    CHECK_THAT(record["ExtendedCellLength"].get<double>(), WithinAbs(4.0 * r, 1.0e-9));
    CHECK(record["ExtendedUncoveredLength"].get<double>() == 0.0);
    CHECK(record["UnhostedLength"].get<double>() == 0.0);
    REQUIRE(record["Cells"].size() == 4);
    std::set<std::size_t> extended;
    for (const auto &corner : record["Corners"])
    {
      CHECK(corner["Topology"] == "ConvexCorner");
      CHECK_THAT(corner["AngleDegrees"].get<double>(), WithinAbs(90.0, 1.0e-6));
      CHECK_THAT(corner["CornerRadiusOverR"].get<double>(), WithinAbs(r / R, 1.0e-9));
      REQUIRE(corner["Arms"].size() == 2);
      std::set<int> arms;
      for (const auto &arm : corner["Arms"])
      {
        arms.insert(arm["Arm"].get<int>());
        CHECK(arm["Hosted"].get<bool>());
        CHECK_THAT(arm["ExitDistanceOverR"].get<double>(), WithinAbs(1.0, 1.0e-9));
        CHECK_THAT(arm["ClaimEndOverR"].get<double>(), WithinAbs((r + R) / R, 1.0e-7));
        CHECK_THAT(arm["StretchLength"].get<double>(), WithinAbs(r, 1.0e-7));
        CHECK_THAT(arm["ExtendedCellLength"].get<double>(), WithinAbs(r, 1.0e-7));
        CHECK(arm["ExtendedUncoveredLength"].get<double>() == 0.0);
      }
      CHECK(arms == std::set<int>{0, 1});
      REQUIRE(corner["Patches"].size() == 2);
      for (const auto &patch : corner["Patches"])
      {
        extended.insert(patch.get<std::size_t>());
      }
    }
    CHECK(extended.size() == 4);
    for (const auto &cell : record["Cells"])
    {
      CHECK(extended.count(cell["Patch"].get<std::size_t>()) == 1);
      CHECK_THAT(cell["ExtendedLength"].get<double>(), WithinAbs(r, 1.0e-7));
      CHECK(cell["Model"] == "isolated");
    }
    CHECK(manifest["Summary"]["CornerArmExtension"]["Corners"].get<int>() == 2);
    CHECK_THAT(
        manifest["Summary"]["CornerArmExtension"]["ExtendedCellLength"].get<double>(),
        WithinAbs(4.0 * r, 1.0e-9));
    CHECK(manifest["Summary"]["CornerArmExtension"]["UnhostedLength"].get<double>() == 0.0);
    // Perpendicular arms: no F1 trim; nothing uncovered.
    CHECK(manifest["Identification"]["Diagnostics"]["CornerArmTrim"]["Count"].get<int>() ==
          0);
    CHECK(manifest["Summary"]["Uncovered"]["Features"].get<int>() == 0);
    // Every straight cell lies at least R from both virtual corners; exactly one cell end
    // per arm (4) sits at R: the extended cells, whose dry-run weight formula holds.
    int corner_rows = 0, at_exit = 0;
    for (const auto &row : rows)
    {
      if (row.topology == "convex corner")
      {
        corner_rows++;
        continue;
      }
      REQUIRE(row.topology == "isolated edge");
      const double distance = std::min(NearestCellEndDistance(row, virtual_corners[0]),
                                       NearestCellEndDistance(row, virtual_corners[1]));
      CHECK(distance >= R - 1.0e-7);
      if (std::abs(distance - R) < 1.0e-7)
      {
        at_exit++;
        CHECK(extended.count(row.patch) == 1);
        CHECK_THAT(row.quadrature_weight * (row.s1 - row.s0),
                   WithinAbs(row.strip[1] - row.strip[0], 1.0e-9));
      }
      else
      {
        CHECK(extended.count(row.patch) == 0);
      }
    }
    CHECK(corner_rows == 2);
    CHECK(at_exit == 4);
    // The portion quadrature sums: 1 + extended / portion length on every portion adjacent
    // to a claim end (the A7 rule read with the extension record), 1 elsewhere; the four
    // stretches add up to 4 t_d.
    int extended_portions = 0;
    double gained_total = 0.0;
    for (const auto &[key, sum] : PortionQuadratureSums(rows))
    {
      const auto &[feature, segment, s0, s1] = key;
      if (std::abs(sum - 1.0) > 1.0e-9)
      {
        extended_portions++;
        const double gained = (sum - 1.0) * (s1 - s0);
        CHECK(gained > 0.0);
        gained_total += gained;
      }
    }
    CHECK(extended_portions == 4);
    CHECK_THAT(gained_total, WithinAbs(4.0 * r, 1.0e-7));
    // The operator applies the same extended cells and carries the same record; the
    // geometry cache it writes carries their own-segment pre-image.
    const auto cache_path = temp.temp_dir / "response-geometry-corner-arm-extension.json";
    test::GeometryCacheEnvGuard cache_env(cache_path.string(), true);
    std::vector<std::unique_ptr<Mesh>> meshes;
    meshes.push_back(std::make_unique<Mesh>(MakeRoundedLeadMesh(x_min, x_max)));
    LaplaceOperator laplace(iodata, meshes);
    SurfaceResponseOperator response(iodata, laplace);
    const auto statistics = response.GetStatistics();
    const auto &operator_record = statistics["Diagnostics"]["CornerArmExtension"];
    REQUIRE(operator_record["Count"].get<int>() == 2);
    CHECK(operator_record["ArmsExtended"].get<int>() == 4);
    CHECK_THAT(operator_record["ExtendedCellLength"].get<double>(),
               WithinAbs(4.0 * r, 1.0e-9));
    REQUIRE(operator_record["Cells"].size() == 4);
    for (const auto &cell : operator_record["Cells"])
    {
      CHECK(extended.count(cell["Patch"].get<std::size_t>()) == 1);
    }
    CHECK(statistics["Diagnostics"]["Uncovered"]["Count"].get<int>() == 0);
    // The own-segment pre-image (decision 537) of every extended cell follows the
    // extension's clip (ClipOwnCell extrapolates the recorded pre-image affinely to the
    // kept offset below the old cell): the cached OwnCell begins exactly R from the cell's
    // virtual corner on the arm line, its far end is the cell's (R + the cell length, at
    // or beyond the claim end), and both ends equal the frame reconstruction (origin +
    // EdgeOffset AxisU + c AxisW), the rule the trims locate cells by on a straight arm.
    Mpi::Barrier(Mpi::World());
    std::ifstream cache_input(cache_path);
    REQUIRE(cache_input);
    const json cache = json::parse(cache_input);
    CHECK(cache["Version"] == 15);
    const auto &cached_patches = cache["Patches"];
    for (const std::size_t p : extended)
    {
      REQUIRE(p < cached_patches.size());
      const auto &cached = cached_patches[p];
      REQUIRE(!cached["OwnCell"].is_null());
      const auto own_cell = cached["OwnCell"].get<std::array<std::array<double, 3>, 2>>();
      const auto cell = cached["LongitudinalCell"].get<std::array<double, 2>>();
      const auto origin = cached["Origin"].get<std::array<double, 3>>();
      const auto axis_u = cached["AxisU"].get<std::array<double, 3>>();
      const auto axis_w = cached["AxisW"].get<std::array<double, 3>>();
      const double edge_offset = cached["EdgeOffset"].get<double>();
      for (int e = 0; e < 2; e++)
      {
        for (int d = 0; d < 3; d++)
        {
          CHECK_THAT(
              own_cell[e][d],
              WithinAbs(origin[d] + edge_offset * axis_u[d] + cell[e] * axis_w[d], 1.0e-9));
        }
      }
      std::array<double, 2> distance{};
      std::array<std::size_t, 2> corner{};
      for (int e = 0; e < 2; e++)
      {
        const double d0 = Distance(own_cell[e], virtual_corners[0]);
        const double d1 = Distance(own_cell[e], virtual_corners[1]);
        corner[e] = d0 < d1 ? 0 : 1;
        distance[e] = std::min(d0, d1);
      }
      const int near = distance[0] < distance[1] ? 0 : 1;
      CHECK_THAT(distance[near], WithinAbs(R, 1.0e-7));
      CHECK_THAT(Distance(own_cell[1 - near], virtual_corners[corner[near]]),
                 WithinAbs(R + (cell[1] - cell[0]), 1.0e-7));
      CHECK(Distance(own_cell[1 - near], virtual_corners[corner[near]]) >= r + R - 1.0e-7);
    }
  }

  SECTION("short arms: the stretch without a host is recorded, never silent")
  {
    // The cap at z = 0.25: each long arm runs 0.125 < R from its tangent point to the
    // domain cut z = 0, so the corner claims the whole arm (claim end t_d + 0.125 = 0.25
    // from the virtual corner, 0.05 beyond the square exit R) and no straight cell nor
    // uncovered portion begins at the claim end: the stretch is unhosted and recorded
    // (Hosted false, UnhostedLength), with the cap arms extended as above. (A short cap
    // between the two fillets is a different case: the identification joins the corners
    // into one spatial cluster, which owns its fillets as arc portions; design 7.2.)
    const json manifest = Preflight(0.125, 0.875, 0.25);
    const auto &record = manifest["Identification"]["Diagnostics"]["CornerArmExtension"];
    REQUIRE(record["Count"].get<int>() == 2);
    CHECK(record["ArmsExtended"].get<int>() == 2);
    CHECK(record["ArmsUnhosted"].get<int>() == 2);
    CHECK_THAT(record["UnhostedLength"].get<double>(), WithinAbs(2.0 * 0.05, 1.0e-7));
    CHECK_THAT(record["ExtendedCellLength"].get<double>(), WithinAbs(2.0 * r, 1.0e-7));
    CHECK_THAT(record["StretchLength"].get<double>(),
               WithinAbs(2.0 * r + 2.0 * 0.05, 1.0e-7));
    for (const auto &corner : record["Corners"])
    {
      REQUIRE(corner["Arms"].size() == 2);
      int hosted = 0, unhosted = 0;
      for (const auto &arm : corner["Arms"])
      {
        if (arm["Hosted"].get<bool>())
        {
          hosted++;
          CHECK_THAT(arm["ExtendedCellLength"].get<double>(), WithinAbs(r, 1.0e-7));
          CHECK_THAT(arm["StretchLength"].get<double>(), WithinAbs(r, 1.0e-7));
        }
        else
        {
          unhosted++;
          CHECK(arm["ExtendedCellLength"].get<double>() == 0.0);
          CHECK(arm["ExtendedUncoveredLength"].get<double>() == 0.0);
          CHECK_THAT(arm["ClaimEndOverR"].get<double>(), WithinAbs(0.25 / R, 1.0e-7));
          CHECK_THAT(arm["StretchLength"].get<double>(), WithinAbs(0.05, 1.0e-7));
        }
      }
      CHECK(hosted == 1);
      CHECK(unhosted == 1);
    }
    CHECK_THAT(manifest["Summary"]["CornerArmExtension"]["UnhostedLength"].get<double>(),
               WithinAbs(2.0 * 0.05, 1.0e-7));
    CHECK(manifest["Summary"]["Uncovered"]["Features"].get<int>() == 0);
  }
#endif
}

// Decision 520 MAJOR-1: a matched ROUNDED virtual (mirror-formed) corner, placed
// HalfByMirror with its real arm's cells from s_half = (R + s) / 2 (the mirror arm trim),
// has its real-arm stretch [s_half, t_d + R) extended like a whole corner's (the 5ecabc3aae
// extension skipped every virtual corner silently). The ramp-sheared lead whose cap ends in
// a half fillet on the 30-degree-tilted wall: the rounded 120-degree virtual corner (the
// cap + its image) is Matched to the r = 0.125 model and HalfByMirror, its one real arm
// (the cap) extended from the claim end t_d + R = 1.361 R back to s_half = 1.077 R and
// recorded HalfByMirror; the cap's cells end exactly at s_half from the virtual corner
// (the stretch once, in the extended cell); the vertex coupon's raw claims stop at s_half
// (the operator record and the geometry cache); the mirror arm trim removes nothing from
// the extended cell (the cap's cells begin beyond s_half). The lead's left end (the sharp
// 120-degree corner, the tilted edge reaching the bottom cut beside the right wall) forms
// an unmergeable mirror configuration whose touched real features are Status Unmerged — the
// decision-512 lane's topology, incidental here and not asserted: the assertions below read
// the virtual corner's record, the cap cells' geometry and the raw claims only.
TEST_CASE_METHOD(test::SurfaceResponseFiles,
                 "SurfaceResponseOperator corner-arm extension of a virtual rounded corner",
                 "[surfaceresponseoperator][cornerarmtrim][mirror][3d][Serial][Parallel]")
{
#if !defined(MFEM_USE_GSLIB)
  SKIP("SurfaceResponseOperator requires MFEM_USE_GSLIB");
#else
  constexpr double R = 0.2, r = 0.125;
  const auto basis_path = temp.temp_dir / "virtual-rounded-basis-points.csv";
  const auto library_path = temp.temp_dir / "fabrication-process-virtual-rounded-3d.json";
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
    std::ifstream rounded_input(rounded_library_3d_path),
        sharp_input(convex_library_3d_path);
    REQUIRE(rounded_input);
    REQUIRE(sharp_input);
    json library = json::parse(rounded_input);
    const json sharp = json::parse(sharp_input);
    library["Name"] = "unit-test-process-virtual-rounded-3d";
    // The sharp 120-degree model (the cap's left corner and the virtual corner at the
    // bottom cut; F1's) and the rounded 120-degree model (the virtual rounded corner)
    // beside the rounded 90.
    for (const auto &model : sharp["Models"])
    {
      if (model["Topology"] == "ConvexCorner")
      {
        json sharp_120 = model;
        sharp_120["Name"] = "convex-corner-120";
        sharp_120["Angle"] = 120.0;
        library["Models"].push_back(sharp_120);
      }
    }
    for (const auto &model : json(library["Models"]))
    {
      if (model["Name"] == "convex-corner-90-r0.125")
      {
        json rounded_120 = model;
        rounded_120["Name"] = "convex-corner-120-r0.125";
        rounded_120["Angle"] = 120.0;
        library["Models"].push_back(rounded_120);
      }
    }
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
  config["Boundaries"]["Ground"]["Attributes"] = {4};
  auto &correction = config["Solver"]["Electrostatic"]["ResponseCorrection"];
  correction.erase("PatchConstruction");
  correction["Library"] = library_path.string();
  correction["TraceCoupling"] = "SurfaceMortar";
  correction["MortarOversampling"] = 2;
  IoData iodata(config, false);
  iodata.boundaries.cracked_attributes.insert(9);
  const auto manifest_path = temp.temp_dir / "surface-response-requirements-virtual.json";
  const auto patches_path = temp.temp_dir / "surface-response-patches.csv";
  const std::array<double, 3> rounded_corner = {1.0, 0.5, 0.5};
  const std::array<double, 3> whole_corner = {kVirtualLeadLeft, 0.5, 0.5};
  const double theta = 120.0 * M_PI / 180.0, t_d = r / std::tan(0.5 * theta);
  const double s_1 = R / std::max(std::abs(std::cos(theta)), std::abs(std::sin(theta)));
  const double s_half = 0.5 * (R + s_1), claim_end = t_d + R;
  auto Distance = [](const std::array<double, 3> &a, const std::array<double, 3> &b)
  {
    double d2 = 0.0;
    for (int d = 0; d < 3; d++)
    {
      d2 += (a[d] - b[d]) * (a[d] - b[d]);
    }
    return std::sqrt(d2);
  };
  json manifest;
  {
    Mesh dry_run_mesh(MakeVirtualRoundedLeadMesh());
    WriteSurfaceResponseRequirements(iodata, dry_run_mesh, manifest_path.string());
    Mpi::Barrier(Mpi::World());
    std::ifstream input(manifest_path);
    REQUIRE(input);
    manifest = json::parse(input);
  }
  const auto rows = ReadDryRun(patches_path);
  int rounded_id = -1, whole_id = -1;
  for (const auto &feature : manifest["Identification"]["Features"])
  {
    if (feature["Type"] != "ConvexCorner")
    {
      continue;
    }
    REQUIRE(feature["Match"]["Status"] == "Matched");
    const auto origin = feature["Frame"]["Origin"].get<std::array<double, 3>>();
    CHECK_THAT(feature["Signature"]["AngleDegrees"].get<double>(),
               WithinAbs(120.0, 1.0e-6));
    if (feature["Match"]["Model"] == "convex-corner-120-r0.125")
    {
      rounded_id = feature["Id"];
      CHECK(feature["Mirror"]["Status"] == "Modelled");
      CHECK(Distance(origin, rounded_corner) < 1.0e-6);
      CHECK_THAT(feature["Signature"]["CornerRadiusOverR"].get<double>(),
                 WithinAbs(r / R, 1.0e-6));
    }
    else if (Distance(origin, whole_corner) < 1.0e-6)
    {
      CHECK(feature["Match"]["Model"] == "convex-corner-120");
      whole_id = feature["Id"];
    }
  }
  REQUIRE(rounded_id >= 0);
  REQUIRE(whole_id >= 0);
  // The extension record: the rounded virtual corner alone, HalfByMirror, its one real arm
  // (the cap) from s_half with the stretch claim end - s_half hosted by the cap's cell.
  const auto &record = manifest["Identification"]["Diagnostics"]["CornerArmExtension"];
  REQUIRE(record["Count"].get<int>() == 1);
  CHECK(record["ArmsExtended"].get<int>() == 1);
  CHECK(record["ArmsUnhosted"].get<int>() == 0);
  CHECK(record["UnhostedLength"].get<double>() == 0.0);
  const double stretch = claim_end - s_half;
  CHECK_THAT(record["ExtendedCellLength"].get<double>(), WithinAbs(stretch, 1.0e-6));
  REQUIRE(record["Cells"].size() == 1);
  const json *corner = &record["Corners"][0];
  CHECK((*corner)["Feature"].get<int>() == rounded_id);
  CHECK((*corner)["HalfByMirror"].get<bool>());
  CHECK_THAT((*corner)["AngleDegrees"].get<double>(), WithinAbs(120.0, 1.0e-6));
  REQUIRE((*corner)["Arms"].size() == 1);
  const auto &arm = (*corner)["Arms"][0];
  CHECK(arm["Hosted"].get<bool>());
  CHECK_THAT(arm["ExitDistanceOverR"].get<double>(), WithinAbs(s_half / R, 1.0e-9));
  CHECK_THAT(arm["ClaimEndOverR"].get<double>(), WithinAbs(claim_end / R, 1.0e-6));
  CHECK_THAT(arm["StretchLength"].get<double>(), WithinAbs(stretch, 1.0e-6));
  CHECK_THAT(arm["ExtendedCellLength"].get<double>(), WithinAbs(stretch, 1.0e-6));
  const auto direction = arm["Direction"].get<std::array<double, 3>>();
  CHECK_THAT(direction[0], WithinAbs(-1.0, 1.0e-6));  // the cap, away from the wall
  REQUIRE((*corner)["Patches"].size() == 1);
  const std::size_t extended = (*corner)["Patches"][0].get<std::size_t>();
  // The virtual rounded corner's mirror arm trim record removes nothing from the extended
  // cell (the cap's cells begin beyond s_half); nothing uncovered.
  const auto &mirror_record = manifest["Identification"]["Diagnostics"]["MirrorArmTrim"];
  bool rounded_mirror_record = false;
  for (const auto &entry : mirror_record["Corners"])
  {
    if (entry["Feature"].get<int>() == rounded_id)
    {
      rounded_mirror_record = true;
      CHECK_THAT(entry["HalfStartOverR"].get<double>(), WithinAbs(s_half / R, 1.0e-9));
      CHECK(entry["Patches"].empty());
      // The record's real-arm direction is the cap (-x), not a chord of the real half arc
      // (decision 533 MINOR-3: ApplyMirrorArmTrim reads the arm from the first real
      // portion; pinned here).
      const auto mirror_arm = entry["Arm"].get<std::array<double, 3>>();
      CHECK_THAT(mirror_arm[0], WithinAbs(-1.0, 1.0e-6));
      CHECK_THAT(mirror_arm[2], WithinAbs(0.0, 1.0e-6));
    }
  }
  CHECK(rounded_mirror_record);
  for (const auto &cell : mirror_record["Cells"])
  {
    CHECK(cell["Patch"].get<std::size_t>() != extended);
  }
  CHECK(manifest["Summary"]["Uncovered"]["Features"].get<int>() == 0);
  // The cap's cells (AxisW along x): every cell end at least s_half from the virtual
  // corner, exactly one at s_half (the extended cell, with the dry-run weight formula
  // holding on it: the stretch [s_half, claim end) counted once, in that cell), the cell
  // ends of the extended cell at s_half and at its old far end; the virtual corner's
  // patch at weight 1 / 2 (HalfByMirror).
  int at_rounded_exit = 0, rounded_rows = 0;
  for (const auto &row : rows)
  {
    if (row.topology == "convex corner")
    {
      if (row.feature == rounded_id)
      {
        rounded_rows++;
        CHECK_THAT(row.weight, WithinAbs(0.5, 1.0e-12));
      }
      continue;
    }
    REQUIRE(row.topology == "isolated edge");
    if (std::abs(std::abs(row.axis_w[0]) - 1.0) > 1.0e-6)
    {
      continue;  // the tilted left edge's cells
    }
    const double to_rounded = NearestCellEndDistance(row, rounded_corner);
    CHECK(to_rounded >= s_half - 1.0e-7);
    if (std::abs(to_rounded - s_half) < 1.0e-7)
    {
      at_rounded_exit++;
      CHECK(row.patch == extended);
      CHECK_THAT(row.quadrature_weight * (row.s1 - row.s0),
                 WithinAbs(row.strip[1] - row.strip[0], 1.0e-9));
      CHECK(row.strip[1] - row.strip[0] > stretch);
    }
  }
  CHECK(rounded_rows == 1);
  CHECK(at_rounded_exit == 1);
  // The operator: the same record, and (through the geometry cache) the rounded vertex
  // coupon's raw claims end at the exit: nothing of the stretch is claimed raw.
  const auto cache_path = temp.temp_dir / "virtual-rounded-geometry-cache.json";
  {
    test::GeometryCacheEnvGuard cache_env(cache_path.string(), true);
    std::vector<std::unique_ptr<Mesh>> meshes;
    meshes.push_back(std::make_unique<Mesh>(MakeVirtualRoundedLeadMesh()));
    LaplaceOperator laplace(iodata, meshes);
    SurfaceResponseOperator response(iodata, laplace);
    const auto statistics = response.GetStatistics();
    const auto &operator_record = statistics["Diagnostics"]["CornerArmExtension"];
    REQUIRE(operator_record["Count"].get<int>() == 1);
    CHECK(operator_record["Corners"][0]["HalfByMirror"].get<bool>());
    CHECK_THAT(operator_record["ExtendedCellLength"].get<double>(),
               WithinAbs(stretch, 1.0e-6));
    Mpi::Barrier(Mpi::World());
  }
  if (Mpi::Root(Mpi::World()))
  {
    std::ifstream input(cache_path);
    REQUIRE(input);
    const json cache = json::parse(input);
    int rounded_patches = 0;
    for (const auto &patch : cache["Patches"])
    {
      if (patch["Feature"].get<int>() != rounded_id || patch["Segment"].get<int>() >= 0)
      {
        continue;
      }
      rounded_patches++;
      CHECK_THAT(patch["Weight"].get<double>(), WithinAbs(0.5, 1.0e-12));
      REQUIRE(!patch["RawClaims"].empty());
      double claimed = 0.0;
      for (const auto &claim : patch["RawClaims"])
      {
        const auto p0 = claim["P0"].get<std::array<double, 3>>();
        const auto p1 = claim["P1"].get<std::array<double, 3>>();
        CHECK(claim["Segment"].get<int>() >= 0);  // no synthetic stretch: none trimmed
        CHECK(Distance(p0, rounded_corner) <= s_half + 1.0e-7);
        CHECK(Distance(p1, rounded_corner) <= s_half + 1.0e-7);
        claimed += Distance(p0, p1);
      }
      // The real half arc (two chords of 15 degrees) + the cap from the tangent point to
      // the exit, s_half - t_d.
      CHECK_THAT(
          claimed,
          WithinAbs(2.0 * 2.0 * r * std::sin(7.5 * M_PI / 180.0) + (s_half - t_d), 1.0e-6));
    }
    CHECK(rounded_patches == 1);
  }
#endif
}

// Decision 520 MINOR-1: the non-90 rounded branches on one fixture, the sheared rounded
// lead with a true circular fillet (r = 0.125) at a 120- or 135-degree convex corner. At
// 120 degrees (s_1 = 1.1547 R) the claim end t_d + R = 1.361 R lies beyond both exits: both
// arms are extended (arm 0 by t_d, arm 1 by t_d + R - s_1) and the F1 trim record removes
// nothing. At 135 degrees (s_1 = 1.4142 R) the claim end 1.259 R lies beyond s_0 = R but
// before s_1: arm 0 is extended by t_d, arm 1 is F1-trimmed from the claim end to s_1 and
// the vertex coupon's raw stretch starts at the CLAIM END, not at R (FillVertexRawClaims).
TEST_CASE_METHOD(test::SurfaceResponseFiles,
                 "SurfaceResponseOperator corner-arm extension of a non-90 rounded corner",
                 "[surfaceresponseoperator][cornerarmtrim][3d][Serial][Parallel]")
{
#if !defined(MFEM_USE_GSLIB)
  SKIP("SurfaceResponseOperator requires MFEM_USE_GSLIB");
#else
  constexpr double R = 0.2, r = 0.125;
  const auto basis_path = temp.temp_dir / "non90-rounded-basis-points.csv";
  const auto library_path = temp.temp_dir / "fabrication-process-non90-rounded-3d.json";
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
    std::ifstream input(rounded_library_3d_path);
    REQUIRE(input);
    json library = json::parse(input);
    library["Name"] = "unit-test-process-non90-rounded-3d";
    std::vector<json> obtuse;
    for (auto &model : library["Models"])
    {
      model["Interfaces"] = {{{"Type", "SA"}, {"Coupon", 1}}};
      if (model["Topology"] == "IsolatedEdge")
      {
        model["BasisPoints"] = basis_path.string();
      }
      if (model["Topology"] == "ConvexCorner")
      {
        for (const double angle : {120.0, 135.0})
        {
          json rounded = model;
          rounded["Name"] = fmt::format("convex-corner-{:g}-r0.125", angle);
          rounded["Angle"] = angle;
          obtuse.push_back(rounded);
        }
      }
    }
    for (auto &model : obtuse)
    {
      library["Models"].push_back(std::move(model));
    }
    std::ofstream output(library_path);
    output << library.dump(2) << "\n";
  }
  Mpi::Barrier(Mpi::World());
  json config = IslandConfig();
  config["Boundaries"]["Ground"]["Attributes"] = {4};
  auto &correction = config["Solver"]["Electrostatic"]["ResponseCorrection"];
  correction.erase("PatchConstruction");
  correction["Library"] = library_path.string();
  correction["TraceCoupling"] = "SurfaceMortar";
  correction["MortarOversampling"] = 2;
  // The lead reaches the wall x = 1: the mirror is off (its own tests above and in
  // test-domainboundary.cpp); this test is about the whole rounded corner at x_min.
  correction["DomainBoundary"] = {{"Mirror", "Off"}};
  IoData iodata(config, false);
  iodata.boundaries.cracked_attributes.insert(9);
  const auto manifest_path = temp.temp_dir / "surface-response-requirements-non90.json";
  const std::array<double, 3> C = {0.125, 0.5, 0.5};
  auto Distance = [](const std::array<double, 3> &a, const std::array<double, 3> &b)
  {
    double d2 = 0.0;
    for (int d = 0; d < 3; d++)
    {
      d2 += (a[d] - b[d]) * (a[d] - b[d]);
    }
    return std::sqrt(d2);
  };
  auto Run = [&](double shear, double angle_degrees)
  {
    const double theta = angle_degrees * M_PI / 180.0;
    const double t_d = r / std::tan(0.5 * theta);
    const double s_1 = R / std::max(std::abs(std::cos(theta)), std::abs(std::sin(theta)));
    const double claim_end = t_d + R;
    json manifest;
    {
      Mesh dry_run_mesh(MakeShearedRoundedLeadMesh(shear));
      WriteSurfaceResponseRequirements(iodata, dry_run_mesh, manifest_path.string());
      Mpi::Barrier(Mpi::World());
      std::ifstream input(manifest_path);
      REQUIRE(input);
      manifest = json::parse(input);
    }
    int corner_id = -1;
    for (const auto &feature : manifest["Identification"]["Features"])
    {
      if (feature["Type"] == "ConvexCorner")
      {
        REQUIRE(corner_id < 0);
        corner_id = feature["Id"];
        CHECK(feature["Match"]["Status"] == "Matched");
        CHECK(feature["Match"]["Model"] ==
              fmt::format("convex-corner-{:g}-r0.125", angle_degrees));
        CHECK_THAT(feature["Signature"]["AngleDegrees"].get<double>(),
                   WithinAbs(angle_degrees, 1.0e-6));
        CHECK(Distance(feature["Frame"]["Origin"].get<std::array<double, 3>>(), C) <
              1.0e-6);
      }
    }
    REQUIRE(corner_id >= 0);
    const auto &record = manifest["Identification"]["Diagnostics"]["CornerArmExtension"];
    const auto &trim = manifest["Identification"]["Diagnostics"]["CornerArmTrim"];
    REQUIRE(record["Count"].get<int>() == 1);
    CHECK_FALSE(record["Corners"][0]["HalfByMirror"].get<bool>());
    CHECK(record["ArmsUnhosted"].get<int>() == 0);
    REQUIRE(trim["Count"].get<int>() == 1);  // the F1 record of every non-90 matched corner
    CHECK_THAT(trim["TrimmedLength"].get<double>(), WithinAbs(s_1 - R, 1.0e-7));
    std::map<int, json> arms;
    for (const auto &arm : record["Corners"][0]["Arms"])
    {
      arms[arm["Arm"].get<int>()] = arm;
    }
    // Arm 0 (the first arm, exit R): extended by t_d in both cases.
    REQUIRE(arms.count(0) == 1);
    CHECK(arms.at(0)["Hosted"].get<bool>());
    CHECK_THAT(arms.at(0)["ExitDistanceOverR"].get<double>(), WithinAbs(1.0, 1.0e-9));
    CHECK_THAT(arms.at(0)["ClaimEndOverR"].get<double>(), WithinAbs(claim_end / R, 1.0e-6));
    CHECK_THAT(arms.at(0)["ExtendedCellLength"].get<double>(), WithinAbs(t_d, 1.0e-6));
    // The vertex coupon's raw claims through the geometry cache.
    const auto cache_path =
        temp.temp_dir / fmt::format("non90-cache-{:g}.json", angle_degrees);
    {
      test::GeometryCacheEnvGuard cache_env(cache_path.string(), true);
      std::vector<std::unique_ptr<Mesh>> meshes;
      meshes.push_back(std::make_unique<Mesh>(MakeShearedRoundedLeadMesh(shear)));
      LaplaceOperator laplace(iodata, meshes);
      SurfaceResponseOperator response(iodata, laplace);
      const auto statistics = response.GetStatistics();
      CHECK(statistics["Diagnostics"]["CornerArmExtension"]["Count"].get<int>() == 1);
      CHECK(statistics["Diagnostics"]["CornerArmTrim"]["Count"].get<int>() == 1);
      Mpi::Barrier(Mpi::World());
    }
    std::vector<std::pair<double, double>> synthetic;  // (near, far) distances from C
    double farthest_portion_claim = 0.0;
    if (Mpi::Root(Mpi::World()))
    {
      std::ifstream input(cache_path);
      REQUIRE(input);
      const json cache = json::parse(input);
      int vertex_patches = 0;
      for (const auto &patch : cache["Patches"])
      {
        if (patch["Feature"].get<int>() != corner_id || patch["Segment"].get<int>() >= 0)
        {
          continue;
        }
        vertex_patches++;
        for (const auto &claim : patch["RawClaims"])
        {
          const double d0 = Distance(claim["P0"].get<std::array<double, 3>>(), C);
          const double d1 = Distance(claim["P1"].get<std::array<double, 3>>(), C);
          if (claim["Segment"].get<int>() < 0)
          {
            synthetic.emplace_back(std::min(d0, d1), std::max(d0, d1));
          }
          else
          {
            farthest_portion_claim = std::max({farthest_portion_claim, d0, d1});
          }
        }
      }
      CHECK(vertex_patches == 1);
    }
    if (claim_end > s_1 + 1.0e-9)
    {
      // 120 degrees: arm 1 extended too (exit s_1 > R), the trim removed nothing, no raw
      // stretch, every portion claim within s_1 of the corner.
      REQUIRE(arms.count(1) == 1);
      CHECK(arms.at(1)["Hosted"].get<bool>());
      CHECK_THAT(arms.at(1)["ExitDistanceOverR"].get<double>(), WithinAbs(s_1 / R, 1.0e-9));
      CHECK_THAT(arms.at(1)["ExtendedCellLength"].get<double>(),
                 WithinAbs(claim_end - s_1, 1.0e-6));
      CHECK(record["ArmsExtended"].get<int>() == 2);
      CHECK_THAT(record["ExtendedCellLength"].get<double>(),
                 WithinAbs(t_d + claim_end - s_1, 1.0e-6));
      CHECK(trim["RemovedCellLength"].get<double>() == 0.0);
      if (Mpi::Root(Mpi::World()))
      {
        CHECK(synthetic.empty());
        CHECK(farthest_portion_claim <= s_1 + 1.0e-7);
      }
    }
    else
    {
      // 135 degrees: arm 1 is F1's: its cells cut from the claim end to s_1 and the raw
      // stretch [claim end, s_1) (not [R, s_1)) claimed by the vertex coupon.
      CHECK(arms.count(1) == 0);
      CHECK(record["ArmsExtended"].get<int>() == 1);
      CHECK_THAT(trim["RemovedCellLength"].get<double>(),
                 WithinAbs(s_1 - claim_end, 1.0e-6));
      if (Mpi::Root(Mpi::World()))
      {
        REQUIRE(synthetic.size() == 1);
        CHECK_THAT(synthetic[0].first, WithinAbs(claim_end, 1.0e-6));
        CHECK_THAT(synthetic[0].second, WithinAbs(s_1, 1.0e-6));
        CHECK(farthest_portion_claim <= claim_end + 1.0e-7);
      }
    }
  };
  SECTION("120 degrees: both arms extended, the second from s_1 > R")
  {
    Run(1.0 / std::sqrt(3.0), 120.0);
  }
  SECTION("135 degrees: the first arm extended, the second F1-trimmed from the claim end")
  {
    Run(1.0, 135.0);
  }
#endif
}

// Decision 511 O2 (i), the accounting (design review MINOR-4): the stretch's raw within-R
// energy is present EXACTLY ONCE in the corrected interface energy after the fix. On the
// rounded lead solved with a library carrying the rounded corner model only (the straight
// runs Missing under UnmatchedPolicy Warn), the uncovered portions begin at the claim ends
// and are extended back to the square exits, so the uncovered raw energy equals the raw
// within-R energy of the straight runs INCLUDING the four stretches (decision 366 Region
// entries on the same quadrature: an exact partition), each stretch's energy positive and
// counted once; before the fix the uncovered portions stopped at the claim ends and the
// stretches' energy was in no term. With the full library (corner + isolated edge, nothing
// uncovered) the fixed-trace energy is outside + models exactly with the extended cells in
// the operator record.
TEST_CASE_METHOD(test::SurfaceResponseFiles,
                 "Electrostatic corner-arm extension keeps the stretch energy once",
                 "[electrostaticsolver][surfaceresponseoperator][cornerarmtrim][3d][Serial]"
                 "[Parallel]")
{
#if !defined(MFEM_USE_GSLIB)
  SKIP("SurfaceResponseOperator requires MFEM_USE_GSLIB");
#else
  constexpr double R = 0.2, r = 0.125, x_min = 0.125, x_max = 0.875;
  const fs::path mesh_path = temp.temp_dir / "extension-lead.mesh";
  const auto basis_path = temp.temp_dir / "extension-basis-points.csv";
  const auto full_library_path = temp.temp_dir / "fabrication-process-extension-full.json";
  const auto corner_library_path =
      temp.temp_dir / "fabrication-process-extension-corner.json";
  if (Mpi::Root(Mpi::World()))
  {
    {
      mfem::Mesh serial = MakeSerialLeadMesh(16, x_min, x_max, 0.0, r);
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
    std::ifstream input(rounded_library_3d_path);
    REQUIRE(input);
    json library = json::parse(input);
    library["Name"] = "unit-test-process-extension-full";
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
    json corner = library;
    corner["Name"] = "unit-test-process-extension-corner";
    corner["Models"] = json::array();
    for (const auto &model : library["Models"])
    {
      if (model["Topology"] == "ConvexCorner")
      {
        corner["Models"].push_back(model);
      }
    }
    REQUIRE(corner["Models"].size() == 1);  // the rounded 90-degree model alone
    std::ofstream output(corner_library_path);
    output << corner.dump(2) << "\n";
  }
  Mpi::Barrier(Mpi::World());

  json config = IslandConfig();
  config["Model"]["Mesh"] = mesh_path.string();
  config["Model"]["L0"] = 1.0;
  config["Boundaries"]["Ground"]["Attributes"] = {4};
  // The decision-366 regions on the same quadrature: the three straight runs as the
  // extended uncovered portions tile them (from the square exit R of each virtual corner),
  // their parts up to the claim ends (t_d + R from the virtual corner), and the four
  // stretches [R, t_d + R) themselves.
  const double exit_z = 0.5 - R, claim_z = 0.5 - r - R;
  const double exit_lo = x_min + R, claim_lo = x_min + r + R;
  const double exit_hi = x_max - R, claim_hi = x_max - r - R;
  const std::vector<std::array<double, 6>> regions = {
      {x_min, 0.5, 0.0, x_min, 0.5, exit_z},      // 5: left arm to the exit
      {x_max, 0.5, 0.0, x_max, 0.5, exit_z},      // 6: right arm to the exit
      {exit_lo, 0.5, 0.5, exit_hi, 0.5, 0.5},     // 7: cap between the exits
      {x_min, 0.5, 0.0, x_min, 0.5, claim_z},     // 8: left arm to the claim end
      {x_max, 0.5, 0.0, x_max, 0.5, claim_z},     // 9: right arm to the claim end
      {claim_lo, 0.5, 0.5, claim_hi, 0.5, 0.5},   // 10: cap between the claim ends
      {x_min, 0.5, claim_z, x_min, 0.5, exit_z},  // 11: left stretch
      {x_max, 0.5, claim_z, x_max, 0.5, exit_z},  // 12: right stretch
      {exit_lo, 0.5, 0.5, claim_lo, 0.5, 0.5},    // 13: cap stretch, left corner
      {claim_hi, 0.5, 0.5, exit_hi, 0.5, 0.5}};   // 14: cap stretch, right corner
  auto &dielectric = config["Boundaries"]["Postprocessing"]["Dielectric"];
  const json target = dielectric[0];
  int index = 5;
  for (const auto &segment : regions)
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
  correction["DomainBoundary"] = {{"Mirror", "Off"}};
  config["Solver"]["Linear"] = {{"Tol", 1.0e-12}, {"MaxIts", 400}};
  const fs::path full_dir = temp.temp_dir / "extension-full";
  const fs::path corner_dir = temp.temp_dir / "extension-corner";
  config["Problem"]["Output"] = full_dir.string();
  correction["Library"] = full_library_path.string();
  correction["UnmatchedPolicy"] = "Error";
  RunElectrostatic(config);
  config["Problem"]["Output"] = corner_dir.string();
  correction["Library"] = corner_library_path.string();
  correction["UnmatchedPolicy"] = "Warn";
  RunElectrostatic(config);
  Mpi::Barrier(Mpi::World());
  if (!Mpi::Root(Mpi::World()))
  {
    return;
  }

  auto ModelEnergySum = [&](const fs::path &dir, int evaluation)
  {
    const Table models = LoadCsv(dir / "surface-response-model-energy.csv");
    const Column &eval = ColumnByHeader(models, "evaluation");
    const Column &energy = ColumnByHeader(models, "fabricated surface energy[4] (J)");
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
    const Table edge = LoadCsv(dir / "surface-Q-edge.csv");
    const Column &interface = ColumnByHeader(edge, "interface");
    const Column &outside = ColumnByHeader(edge, "E_out (J)");
    REQUIRE(interface.data.size() == 1);
    CHECK(static_cast<int>(std::lround(interface.data[0])) == 4);
    return outside.data[0];
  };
  auto RawEnergy = [&](const fs::path &dir, int interface_index)
  {
    const Table surface = LoadCsv(dir / "surface-Q.csv");
    const Table domain = LoadCsv(dir / "domain-E.csv");
    const Column &participation =
        ColumnByHeader(surface, fmt::format("p_surf[{}]", interface_index));
    const Column &energy = ColumnByHeader(domain, "E_elec (J)");
    REQUIRE(participation.data.size() == 1);
    REQUIRE(energy.data.size() == 1);
    return participation.data[0] * energy.data[0];
  };
  auto ExtensionRecord = [&](const fs::path &dir)
  {
    std::ifstream metadata_input(dir / "palace.json");
    REQUIRE(metadata_input);
    const auto metadata = json::parse(metadata_input);
    return metadata.at("SurfaceResponse").at("Diagnostics");
  };
  // The Total row of a portion-energy table (uncovered / DomainBoundary raw claims) on the
  // raw field (evaluation 0); 0 when the table was not written (no such portions).
  auto PortionTableTotal = [&](const fs::path &path)
  {
    if (!fs::exists(path))
    {
      return 0.0;
    }
    const TextCsv table = ReadTextCsv(path);
    const std::size_t type_column = table.Column("type");
    const std::size_t evaluation_column = table.Column("evaluation");
    std::size_t energy_column = 0;
    for (std::size_t c = 0; c < table.header.size(); c++)
    {
      if (table.header[c].find("raw energy[4]") != std::string::npos)
      {
        energy_column = c;
      }
    }
    for (const auto &row : table.rows)
    {
      if (std::stoi(row[evaluation_column]) == 0 && row[type_column] == "Total")
      {
        return std::stod(row[energy_column]);
      }
    }
    FAIL("no Total row in " << path.string());
    return 0.0;
  };

  // The full library: nothing uncovered, the extended cells in the record, the fixed-trace
  // energy outside + models exactly.
  CHECK_FALSE(fs::exists(full_dir / "surface-response-uncovered-energy.csv"));
  {
    const Table corrected = LoadCsv(full_dir / "surface-Q-corrected.csv");
    REQUIRE(corrected.n_rows() == 1);
    const double ft =
        ColumnByHeader(corrected, "E_surf postprocessed fixed-trace[4] (J)").data[0];
    // Outside + models (+ the DomainBoundary raw claims of the cut lead end, if any: the
    // F-DB-a term of the same identity).
    const double domain_boundary =
        PortionTableTotal(full_dir / "surface-response-domain-boundary-energy.csv");
    CHECK_THAT(ft, WithinRel(OutsideEnergy(full_dir) + ModelEnergySum(full_dir, 0) +
                                 domain_boundary,
                             1.0e-10));
    const auto diagnostics = ExtensionRecord(full_dir);
    CHECK(diagnostics.at("Uncovered").at("Count").get<int>() == 0);
    const auto &extension = diagnostics.at("CornerArmExtension");
    CHECK(extension.at("Count").get<int>() == 2);
    CHECK(extension.at("ArmsExtended").get<int>() == 4);
    CHECK(extension.at("Cells").size() == 4);
    CHECK_THAT(extension.at("ExtendedCellLength").get<double>(),
               WithinAbs(4.0 * r, 1.0e-9));
    CHECK(extension.at("ExtendedUncoveredLength").get<double>() == 0.0);
  }

  // The corner-only library: the raw outputs are byte-identical to the full run (the same
  // raw solve); the three straight runs are uncovered, their portions extended to the
  // square exits (4 x t_d), and the uncovered raw energy is exactly the regions' raw
  // within-R energy including every stretch, each positive and counted once.
  for (const char *file : {"terminal-C.csv", "terminal-V.csv", "domain-E.csv",
                           "surface-Q.csv", "surface-Q-edge.csv"})
  {
    INFO(file);
    CHECK(ReadFile(full_dir / file) == ReadFile(corner_dir / file));
  }
  {
    const auto diagnostics = ExtensionRecord(corner_dir);
    const auto &extension = diagnostics.at("CornerArmExtension");
    CHECK(extension.at("Count").get<int>() == 2);
    CHECK(extension.at("ArmsExtended").get<int>() == 4);
    CHECK(extension.at("Cells").empty());
    CHECK(extension.at("ExtendedCellLength").get<double>() == 0.0);
    CHECK_THAT(extension.at("ExtendedUncoveredLength").get<double>(),
               WithinAbs(4.0 * r, 1.0e-9));
    const auto &uncovered = diagnostics.at("Uncovered");
    CHECK(uncovered.at("Features").get<int>() == 3);
    // The straight runs from the square exits: 2 x (0.5 - R) + (x_max - x_min - 2 R).
    CHECK_THAT(uncovered.at("Length").get<double>(),
               WithinAbs(2.0 * (0.5 - R) + (x_max - x_min - 2.0 * R), 1.0e-9));
  }
  REQUIRE(fs::is_regular_file(corner_dir / "surface-response-uncovered-energy.csv"));
  const TextCsv uncovered =
      ReadTextCsv(corner_dir / "surface-response-uncovered-energy.csv");
  const std::size_t type_column = uncovered.Column("type");
  const std::size_t evaluation_column = uncovered.Column("evaluation");
  const std::size_t energy_column = uncovered.Column("uncovered raw energy[4] (J)");
  double uncovered_ft = 0.0;
  bool found = false;
  for (const auto &row : uncovered.rows)
  {
    if (std::stoi(row[evaluation_column]) == 0 && row[type_column] == "Total")
    {
      uncovered_ft = std::stod(row[energy_column]);
      found = true;
    }
  }
  REQUIRE(found);
  CHECK(uncovered_ft > 0.0);
  const Table corrected = LoadCsv(corner_dir / "surface-Q-corrected.csv");
  REQUIRE(corrected.n_rows() == 1);
  const double ft =
      ColumnByHeader(corrected, "E_surf postprocessed fixed-trace[4] (J)").data[0];
  CHECK_THAT(
      ft,
      WithinRel(
          OutsideEnergy(corner_dir) + ModelEnergySum(corner_dir, 0) + uncovered_ft +
              PortionTableTotal(corner_dir / "surface-response-domain-boundary-energy.csv"),
          1.0e-10));
  // The uncovered energy = the raw within-R energy of the straight runs from the exits
  // (regions 5-7), = their parts up to the claim ends (8-10) + the four stretches (11-14).
  const double to_exits =
      RawEnergy(corner_dir, 5) + RawEnergy(corner_dir, 6) + RawEnergy(corner_dir, 7);
  const double to_claims =
      RawEnergy(corner_dir, 8) + RawEnergy(corner_dir, 9) + RawEnergy(corner_dir, 10);
  double stretches = 0.0;
  for (int region = 11; region <= 14; region++)
  {
    const double stretch = RawEnergy(corner_dir, region);
    CHECK(stretch > 0.0);
    stretches += stretch;
  }
  CHECK_THAT(uncovered_ft, WithinRel(to_exits, 1.0e-9));
  CHECK_THAT(to_exits, WithinRel(to_claims + stretches, 1.0e-9));
  CHECK(stretches > 0.01 * to_exits);  // the stretches carry a visible share
  // Before the fix the uncovered portions ended at the claim ends: the energy they would
  // have kept is to_claims, short of the stretches now counted once.
  CHECK_THAT(uncovered_ft - to_claims, WithinRel(stretches, 1.0e-6));
#endif
}

// Decision 394 F2: an uncovered requirement (a feature the library has no model for) keeps
// its RAW within-R surface energy in the corrected interface energies instead of losing it
// with the modelled perimeter's within-R energy, and the share is reported per feature
// type. The two-corner lead solved with a library that has the corner model (nothing
// uncovered: the fixed-trace energy is the outside energy plus the model energies exactly
// and no uncovered table is written) and with a library without it (both corners Missing
// under UnmatchedPolicy Warn): the fixed-trace energy is outside + models + the uncovered
// energy, the uncovered energy equals the raw within-R energy left by the isolated-edge
// portions' translational cells (decision 366 Region entries on the same quadrature, an
// exact partition at a 90-degree convex corner), and the table carries the ConvexCorner and
// Total rows for the raw-field (0, 1) and corrected-field (2) evaluations.
TEST_CASE_METHOD(test::SurfaceResponseFiles,
                 "Electrostatic uncovered requirements keep their raw energy",
                 "[electrostaticsolver][surfaceresponseoperator][cornerarmtrim][3d][Serial]"
                 "[Parallel]")
{
#if !defined(MFEM_USE_GSLIB)
  SKIP("SurfaceResponseOperator requires MFEM_USE_GSLIB");
#else
  constexpr double R = 0.2;
  const fs::path mesh_path = temp.temp_dir / "uncovered-lead.mesh";
  const auto basis_path = temp.temp_dir / "uncovered-basis-points.csv";
  const auto full_library_path = temp.temp_dir / "fabrication-process-uncovered-full.json";
  const auto edge_library_path = temp.temp_dir / "fabrication-process-uncovered-edges.json";
  {
    if (Mpi::Root(Mpi::World()))
    {
      {
        mfem::Mesh serial = MakeSerialLeadMesh(0.0, 0.75);
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
      library["Name"] = "unit-test-process-uncovered-full";
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
      edges["Name"] = "unit-test-process-uncovered-edges";
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
  }
  Mpi::Barrier(Mpi::World());

  json config = IslandConfig();
  config["Model"]["Mesh"] = mesh_path.string();
  // The mesh unit is the metre of the fixture libraries' surface-matrix rows (R = 0.2 m).
  config["Model"]["L0"] = 1.0;
  config["Boundaries"]["Ground"]["Attributes"] = {4};
  // The raw within-R partition of decision 366 on the same quadrature: the translational
  // cells (half-open along, transverse <= R) of the three isolated-edge portions — the long
  // edges from the cut z = 0 to the corner windows, the cap between the windows — as
  // non-target entries; the corners' windows (the nearest perimeter point within R of a
  // vertex along either arm, the exterior sectors included) are their complement within R.
  auto &dielectric = config["Boundaries"]["Postprocessing"]["Dielectric"];
  const json target = dielectric[0];
  int index = 5;
  for (const std::array<double, 6> &segment :
       {std::array<double, 6>{0.25, 0.5, 0.0, 0.25, 0.5, 0.5 - R},
        std::array<double, 6>{0.25 + R, 0.5, 0.5, 0.75 - R, 0.5, 0.5},
        std::array<double, 6>{0.75, 0.5, 0.0, 0.75, 0.5, 0.5 - R}})
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
  const fs::path full_dir = temp.temp_dir / "uncovered-full";
  const fs::path edges_dir = temp.temp_dir / "uncovered-edges";
  config["Problem"]["Output"] = full_dir.string();
  correction["Library"] = full_library_path.string();
  correction["UnmatchedPolicy"] = "Error";
  RunElectrostatic(config);
  config["Problem"]["Output"] = edges_dir.string();
  correction["Library"] = edge_library_path.string();
  correction["UnmatchedPolicy"] = "Warn";
  RunElectrostatic(config);
  Mpi::Barrier(Mpi::World());
  if (!Mpi::Root(Mpi::World()))
  {
    return;
  }

  auto ModelEnergySum = [&](const fs::path &dir, int evaluation)
  {
    const Table models = LoadCsv(dir / "surface-response-model-energy.csv");
    const Column &eval = ColumnByHeader(models, "evaluation");
    const Column &energy = ColumnByHeader(models, "fabricated surface energy[4] (J)");
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
    const Table edge = LoadCsv(dir / "surface-Q-edge.csv");
    const Column &interface = ColumnByHeader(edge, "interface");
    const Column &outside = ColumnByHeader(edge, "E_out (J)");
    REQUIRE(interface.data.size() == 1);
    CHECK(static_cast<int>(std::lround(interface.data[0])) == 4);
    return outside.data[0];
  };
  // The raw interface energy of an entry: its participation (surface-Q.csv) times the
  // domain energy (domain-E.csv), one source.
  auto RawEnergy = [&](const fs::path &dir, int interface_index)
  {
    const Table surface = LoadCsv(dir / "surface-Q.csv");
    const Table domain = LoadCsv(dir / "domain-E.csv");
    const Column &participation =
        ColumnByHeader(surface, fmt::format("p_surf[{}]", interface_index));
    const Column &energy = ColumnByHeader(domain, "E_elec (J)");
    REQUIRE(participation.data.size() == 1);
    REQUIRE(energy.data.size() == 1);
    return participation.data[0] * energy.data[0];
  };

  // Nothing uncovered: no uncovered table, the record counts 0, and the fixed-trace
  // interface energy is outside + models exactly (the pre-existing identity).
  CHECK_FALSE(fs::exists(full_dir / "surface-response-uncovered-energy.csv"));
  const Table full_corrected = LoadCsv(full_dir / "surface-Q-corrected.csv");
  REQUIRE(full_corrected.n_rows() == 1);
  const double full_ft =
      ColumnByHeader(full_corrected, "E_surf postprocessed fixed-trace[4] (J)").data[0];
  CHECK_THAT(full_ft,
             WithinRel(OutsideEnergy(full_dir) + ModelEnergySum(full_dir, 0), 1.0e-10));
  {
    std::ifstream metadata_input(full_dir / "palace.json");
    REQUIRE(metadata_input);
    const auto metadata = json::parse(metadata_input);
    const auto &diagnostics = metadata.at("SurfaceResponse").at("Diagnostics");
    CHECK(diagnostics.at("Uncovered").at("Count").get<int>() == 0);
    CHECK(diagnostics.at("CornerArmTrim").at("Count").get<int>() == 0);
  }

  // Both corners uncovered: the raw outputs are byte-identical to the fully covered run
  // (the same raw solve), the uncovered table exists, and the fixed-trace energy carries
  // the uncovered raw within-R energy of the two corner windows.
  for (const char *file : {"terminal-C.csv", "terminal-V.csv", "domain-E.csv",
                           "surface-Q.csv", "surface-Q-edge.csv"})
  {
    INFO(file);
    CHECK(ReadFile(full_dir / file) == ReadFile(edges_dir / file));
  }
  REQUIRE(fs::is_regular_file(edges_dir / "surface-response-uncovered-energy.csv"));
  const TextCsv uncovered =
      ReadTextCsv(edges_dir / "surface-response-uncovered-energy.csv");
  const std::size_t type_column = uncovered.Column("type");
  const std::size_t evaluation_column = uncovered.Column("evaluation");
  const std::size_t portions_column = uncovered.Column("portions");
  const std::size_t length_column = uncovered.Column("length (m)");
  const std::size_t energy_column = uncovered.Column("uncovered raw energy[4] (J)");
  const std::size_t share_column = uncovered.Column("uncovered share[4]");
  // One source x 3 evaluations x (ConvexCorner + Total).
  REQUIRE(uncovered.rows.size() == 6);
  std::map<std::pair<int, std::string>, std::vector<std::string>> by_key;
  for (const auto &row : uncovered.rows)
  {
    by_key[{std::stoi(row[evaluation_column]), row[type_column]}] = row;
  }
  REQUIRE(by_key.count({0, "ConvexCorner"}) == 1);
  REQUIRE(by_key.count({0, "Total"}) == 1);
  REQUIRE(by_key.count({1, "Total"}) == 1);
  REQUIRE(by_key.count({2, "Total"}) == 1);
  // The two corners' windows: R along each arm (two mesh segments of h = 0.125 each), 0.8
  // in total.
  CHECK(std::stoi(by_key.at({0, "ConvexCorner"})[portions_column]) == 8);
  CHECK_THAT(std::stod(by_key.at({0, "ConvexCorner"})[length_column]),
             WithinAbs(4.0 * R, 1.0e-9));
  CHECK_THAT(std::stod(by_key.at({0, "Total"})[length_column]), WithinAbs(4.0 * R, 1.0e-9));
  const double uncovered_ft = std::stod(by_key.at({0, "Total"})[energy_column]);
  CHECK(uncovered_ft > 0.0);
  CHECK(std::stod(by_key.at({0, "ConvexCorner"})[energy_column]) == uncovered_ft);
  CHECK(std::stod(by_key.at({1, "Total"})[energy_column]) == uncovered_ft);
  const Table edges_corrected = LoadCsv(edges_dir / "surface-Q-corrected.csv");
  REQUIRE(edges_corrected.n_rows() == 1);
  const double edges_ft =
      ColumnByHeader(edges_corrected, "E_surf postprocessed fixed-trace[4] (J)").data[0];
  const double edges_ff =
      ColumnByHeader(edges_corrected, "E_surf postprocessed fixed-flux[4] (J)").data[0];
  const double edges_sc =
      ColumnByHeader(edges_corrected, "E_surf corrected[4] (J)").data[0];
  const double outside = OutsideEnergy(edges_dir);
  CHECK_THAT(edges_ft,
             WithinRel(outside + ModelEnergySum(edges_dir, 0) + uncovered_ft, 1.0e-10));
  CHECK_THAT(edges_ff,
             WithinRel(outside + ModelEnergySum(edges_dir, 1) + uncovered_ft, 1.0e-10));
  CHECK_THAT(std::stod(by_key.at({0, "Total"})[share_column]),
             WithinRel(uncovered_ft / edges_ft, 1.0e-9));
  CHECK_THAT(std::stod(by_key.at({1, "Total"})[share_column]),
             WithinRel(uncovered_ft / edges_ff, 1.0e-9));
  // The uncovered energy is the raw within-R energy of the corner windows: the raw
  // interface energy minus the outside energy minus the three isolated-edge portions'
  // cells (the decision-366 regions, indices 5-7).
  const double within = RawEnergy(edges_dir, 4) - outside;
  const double cells =
      RawEnergy(edges_dir, 5) + RawEnergy(edges_dir, 6) + RawEnergy(edges_dir, 7);
  CHECK(within > cells);
  CHECK_THAT(uncovered_ft, WithinRel(within - cells, 1.0e-9));
  // The self-consistent column (accepted here) carries the corrected field's uncovered
  // energy: finite, positive, of the same order, and its share is against E_sc.
  const double uncovered_sc = std::stod(by_key.at({2, "Total"})[energy_column]);
  REQUIRE(std::isfinite(edges_sc));
  CHECK(uncovered_sc > 0.0);
  CHECK(uncovered_sc != uncovered_ft);
  CHECK_THAT(uncovered_sc, WithinRel(uncovered_ft, 0.5));
  CHECK_THAT(std::stod(by_key.at({2, "Total"})[share_column]),
             WithinRel(uncovered_sc / edges_sc, 1.0e-9));
  {
    std::ifstream metadata_input(edges_dir / "palace.json");
    REQUIRE(metadata_input);
    const auto metadata = json::parse(metadata_input);
    const auto &record = metadata.at("SurfaceResponse").at("Diagnostics").at("Uncovered");
    CHECK(record.at("Count").get<int>() == 8);
    CHECK(record.at("Features").get<int>() == 2);
    CHECK_THAT(record.at("Length").get<double>(), WithinAbs(4.0 * R, 1.0e-9));
    CHECK(record.at("ByType").at("ConvexCorner").at("Portions").get<int>() == 8);
  }
#endif
}

// Decision 399 MAJOR-1: an uncovered portion inside a matched spatial cluster's support box
// would be counted twice (the cluster coupon models its whole box; F2 keeps the raw
// within-R energy of the portion), so the placement clips the uncovered portions by every
// matched cluster's box exactly as the continuation ownership clips the translational cells
// (the same bounds, the same strict-interior test). On the top face of the unit cube (an MA
// interface under an exactly represented field; the perimeter's tree segments are the mesh
// edges) a portion of the edge y = 0 crossing a box face is clipped exactly at the face,
// a portion straddling a box is cut into two pieces, a portion wholly inside is removed, a
// portion touching a face is untouched, and the uncovered energies of the kept pieces and
// the removed piece partition the portion's energy to 1e-10; without a cluster box (no
// support, or a vertex coupon's support without claims) the portions are untouched bitwise.
TEST_CASE("SurfaceResponseOperator uncovered portions clipped by spatial supports",
          "[surfaceresponseoperator][cornerarmtrim][Serial][Parallel]")
{
  using Portion =
      config::ElectrostaticSolverData::ResponseCorrectionData::UncoveredPortionData;
  using Claim =
      config::ElectrostaticSolverData::ResponseCorrectionPatchData::Provenance::Claim;
  auto Same = [](const Portion &a, const Portion &b)
  {
    return a.feature == b.feature && a.topology == b.topology && a.segment == b.segment &&
           a.p0 == b.p0 && a.p1 == b.p1;
  };
  // The portion of the edge y = 0, z = 1 (a tree segment [0, 0.5] of the 2 x 2 x 1 mesh
  // below) from x = 0.05 to x = 0.45, feature 7.
  const Portion portion{7, "SpatialEdgeCluster", 3, {0.05, 0.0, 1.0}, {0.45, 0.0, 1.0}};
  // A matched cluster's box (claims non-empty) whose x face at 0.2137 the portion crosses,
  // R above and below the plane; a vertex coupon's box (no claims) containing the portion.
  SpatialSupportBounds cluster;
  cluster.patch = 2;
  cluster.min = {0.2137, -0.3, 0.75};
  cluster.max = {0.9, 0.3, 1.25};
  cluster.claims = {Claim{3, {0.5, 0.0, 1.0}, {0.9, 0.0, 1.0}}};
  SpatialSupportBounds corner = cluster;
  corner.patch = 5;
  corner.min = {-0.1, -0.3, 0.75};
  corner.max = {0.6, 0.3, 1.25};
  corner.claims.clear();

  SECTION("no cluster box: bitwise untouched")
  {
    std::vector<Portion> portions = {portion, portion};
    portions[1].feature = 8;
    const Portion second = portions[1];
    const auto none = ClipUncoveredPortionsBySpatialSupport(portions, {}, 3);
    CHECK(none.clips.empty());
    CHECK(none.clipped_portions == 0);
    CHECK(none.removed_length == 0.0);
    REQUIRE(portions.size() == 2);
    CHECK(Same(portions[0], portion));
    CHECK(Same(portions[1], second));
    const auto vertex_only = ClipUncoveredPortionsBySpatialSupport(portions, {corner}, 3);
    CHECK(vertex_only.clips.empty());
    CHECK(vertex_only.clipped_portions == 0);
    REQUIRE(portions.size() == 2);
    CHECK(Same(portions[0], portion));
    CHECK(Same(portions[1], second));
  }

  SECTION("a portion crossing a box face is clipped exactly at the face")
  {
    std::vector<Portion> portions = {portion};
    const auto clipping =
        ClipUncoveredPortionsBySpatialSupport(portions, {corner, cluster}, 3);
    REQUIRE(clipping.clips.size() == 1);
    CHECK(clipping.clips[0].feature == 7);
    CHECK(clipping.clips[0].topology == "SpatialEdgeCluster");
    CHECK(clipping.clips[0].segment == 3);
    CHECK(clipping.clips[0].spatial_patch == 2);
    CHECK_THAT(clipping.clips[0].length, WithinAbs(0.45 - 0.2137, 1.0e-12));
    CHECK(clipping.clipped_portions == 1);
    CHECK(clipping.removed_portions == 0);
    CHECK(clipping.split_portions == 0);
    CHECK_THAT(clipping.removed_length, WithinAbs(0.45 - 0.2137, 1.0e-12));
    REQUIRE(clipping.removed_by_feature.size() == 1);
    CHECK_THAT(clipping.removed_by_feature.at(7), WithinAbs(0.45 - 0.2137, 1.0e-12));
    REQUIRE(portions.size() == 1);
    CHECK(portions[0].feature == 7);
    CHECK(portions[0].segment == 3);
    CHECK(portions[0].p0 == portion.p0);  // the untouched end keeps its coordinates
    CHECK_THAT(portions[0].p1[0], WithinAbs(0.2137, 1.0e-12));
    CHECK(portions[0].p1[1] == 0.0);
    CHECK(portions[0].p1[2] == 1.0);
    // The same portion reversed is clipped at the same face.
    std::vector<Portion> reversed = {
        Portion{7, "SpatialEdgeCluster", 3, portion.p1, portion.p0}};
    ClipUncoveredPortionsBySpatialSupport(reversed, {cluster}, 3);
    REQUIRE(reversed.size() == 1);
    CHECK_THAT(reversed[0].p0[0], WithinAbs(0.2137, 1.0e-12));
    CHECK(reversed[0].p1 == portion.p0);
  }

  SECTION("a box inside the portion splits it, a portion inside a box is removed, a face "
          "touched is not inside")
  {
    SpatialSupportBounds middle = cluster;
    middle.min[0] = 0.15;
    middle.max[0] = 0.3;
    SpatialSupportBounds whole = cluster;
    whole.patch = 9;
    whole.min[0] = -0.5;
    whole.max[0] = 0.6;
    SpatialSupportBounds touching = cluster;
    touching.patch = 11;
    touching.min[0] = 0.45;
    touching.max[0] = 0.8;
    std::vector<Portion> portions = {portion, portion, portion};
    portions[1].feature = 8;
    portions[2].feature = 9;
    portions[2].p0 = {0.85, 0.0, 1.0};  // beyond both boxes: wholly outside
    portions[2].p1 = {0.95, 0.0, 1.0};
    const auto clipping =
        ClipUncoveredPortionsBySpatialSupport(portions, {middle, touching}, 3);
    CHECK(clipping.clipped_portions == 2);
    CHECK(clipping.split_portions == 2);
    CHECK(clipping.removed_portions == 0);
    CHECK_THAT(clipping.removed_length, WithinAbs(2.0 * 0.15, 1.0e-12));
    REQUIRE(portions.size() == 5);
    CHECK(portions[0].feature == 7);
    CHECK(portions[0].p0 == portion.p0);
    CHECK_THAT(portions[0].p1[0], WithinAbs(0.15, 1.0e-12));
    CHECK(portions[1].feature == 7);
    CHECK_THAT(portions[1].p0[0], WithinAbs(0.3, 1.0e-12));
    CHECK(portions[1].p1 == portion.p1);
    CHECK(portions[2].feature == 8);
    CHECK(portions[3].feature == 8);
    CHECK(portions[4].feature == 9);
    CHECK(portions[4].p0 == std::array<double, 3>{0.85, 0.0, 1.0});
    std::vector<Portion> removed = {portion};
    const auto wholly = ClipUncoveredPortionsBySpatialSupport(removed, {whole}, 3);
    CHECK(removed.empty());
    CHECK(wholly.removed_portions == 1);
    CHECK(wholly.clipped_portions == 1);
    CHECK_THAT(wholly.removed_length, WithinAbs(0.4, 1.0e-12));
    // Two boxes sharing the portion: the union is removed once, each box records its part.
    std::vector<Portion> shared = {portion};
    const auto two = ClipUncoveredPortionsBySpatialSupport(shared, {middle, cluster}, 3);
    REQUIRE(two.clips.size() == 2);
    CHECK_THAT(two.removed_length, WithinAbs(0.45 - 0.15, 1.0e-12));
    REQUIRE(shared.size() == 1);
    CHECK_THAT(shared[0].p1[0], WithinAbs(0.15, 1.0e-12));
  }

  SECTION("the uncovered energies of the kept and removed pieces partition the portion's")
  {
    auto serial = mfem::Mesh::MakeCartesian3D(2, 2, 1, mfem::Element::TETRAHEDRON);
    int top = -1;
    mfem::Vector center(3);
    for (int be = 0; be < serial.GetNBE(); be++)
    {
      auto *T = serial.GetBdrElementTransformation(be);
      T->Transform(mfem::Geometries.GetCenter(T->GetGeometryType()), center);
      if (std::abs(center(2) - 1.0) < 1e-12)
      {
        top = serial.GetBdrAttribute(be);
        break;
      }
    }
    REQUIRE(top > 0);
    Mesh mesh(std::make_unique<mfem::ParMesh>(Mpi::World(), serial));
    mfem::H1_FECollection h1(2, 3);
    mfem::ND_FECollection nd(2, 3);
    FiniteElementSpace h1_space(mesh, &h1), nd_space(mesh, &nd);
    fem::DefaultIntegrationOrder::p_trial = 2;
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
                                                  E(2) = 1.0 + x(0) + 0.5 * x(1);
                                                });
    field.Real().ProjectCoefficient(coefficient);
    config::InterfaceDielectricData data;
    data.attributes = {top};
    data.type = InterfaceDielectric::MA;
    data.t = 0.2;
    data.epsilon_r = 2.0;
    data.edge_attributes = {top};
    data.edge_distances = {0.25};
    config::BoundaryPostData postpro;
    postpro.dielectric.emplace(1, data);
    SurfacePostOperator surf_post_op(postpro, ProblemType::ELECTROSTATIC, materials,
                                     h1_space, nd_space);
    auto Energy = [&](const std::vector<Portion> &portions)
    {
      std::vector<SurfacePostOperator::UncoveredPerimeterPortion> perimeter;
      for (const auto &entry : portions)
      {
        perimeter.push_back({entry.p0, entry.p1, entry.topology, entry.feature});
      }
      const auto energies =
          surf_post_op.GetInterfaceUncoveredEdgeEnergies({1}, field, nullptr, perimeter);
      REQUIRE(energies.size() == 1);
      return energies.at(1).energy;
    };
    const double whole = Energy({portion});
    REQUIRE(whole > 0.0);
    // Crossing: kept [0.05, 0.2137] + removed [0.2137, 0.45].
    std::vector<Portion> kept = {portion};
    ClipUncoveredPortionsBySpatialSupport(kept, {cluster}, 3);
    REQUIRE(kept.size() == 1);
    const Portion removed{7, "SpatialEdgeCluster", 3, kept[0].p1, portion.p1};
    const double kept_energy = Energy(kept), removed_energy = Energy({removed});
    CHECK(kept_energy > 0.0);
    CHECK(removed_energy > 0.0);
    CHECK_THAT(kept_energy + removed_energy, WithinRel(whole, 1.0e-10));
    // Straddling: two kept pieces + the removed middle.
    SpatialSupportBounds middle = cluster;
    middle.min[0] = 0.15;
    middle.max[0] = 0.3;
    std::vector<Portion> pieces = {portion};
    ClipUncoveredPortionsBySpatialSupport(pieces, {middle}, 3);
    REQUIRE(pieces.size() == 2);
    const Portion middle_piece{7, "SpatialEdgeCluster", 3, pieces[0].p1, pieces[1].p0};
    CHECK_THAT(Energy(pieces) + Energy({middle_piece}), WithinRel(whole, 1.0e-10));
    CHECK(Energy(pieces) < whole);
  }
}

TEST_CASE("SurfaceResponseOperator mirror arm direction of a rounded virtual corner",
          "[surfaceresponseoperator][cornerarmtrim][mirror][Serial][Parallel]")
{
  // A rounded 90-degree virtual corner at the vertex (0, 0) with the arms +x and +y and the
  // fillet radius 0.1 (arc centre (0.1, 0.1)): the real half is the chord of arc 0 between
  // the arc points at 22.5 and 45 degrees (whose far end from the vertex, (0.0617, 0.0076),
  // points 7 degrees off the arm) and the straight arm along +x from the tangent point
  // (0.1, 0) (segment 1); the image arm (segment 2, at or beyond real_segments) along +y.
  // The portions list the chord FIRST (the perimeter order from the vertex): the arm must
  // still be +x (decision 533 MINOR-3), the chord's direction never taken while a straight
  // real portion exists.
  IdentificationResult identification;
  IdentifiedSegment chord, straight, image;
  const double phi_a = 22.5 * M_PI / 180.0, phi_b = 45.0 * M_PI / 180.0;
  const std::array<double, 3> arc_a = {0.1 - 0.1 * std::sin(phi_a),
                                       0.1 - 0.1 * std::cos(phi_a), 0.0};
  const std::array<double, 3> arc_b = {0.1 - 0.1 * std::sin(phi_b),
                                       0.1 - 0.1 * std::cos(phi_b), 0.0};
  chord.key = {arc_b, arc_a};  // lexicographic: arc_b.x < arc_a.x
  chord.arc = 0;
  straight.key = {{{0.1, 0.0, 0.0}, {0.9, 0.0, 0.0}}};
  image.key = {{{0.0, 0.1, 0.0}, {0.0, 0.9, 0.0}}};
  identification.segments = {chord, straight, image};
  identification.real_segments = 2;
  IdentifiedFeature feature;
  feature.origin = {0.0, 0.0, 0.0};
  feature.axes = {{{1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, {0.0, 0.0, 1.0}}};
  IdentifiedPortion on_chord, on_straight, on_image;
  on_chord.segment = 0;
  on_straight.segment = 1;
  on_image.segment = 2;
  feature.portions = {on_chord, on_straight, on_image};
  SECTION("the chord listed first: the straight real portion gives the arm")
  {
    const auto arm =
        MirrorArmDirection(identification, feature, identification.real_segments);
    REQUIRE(arm);
    CHECK_THAT((*arm)[0], Catch::Matchers::WithinAbs(1.0, 1.0e-12));
    CHECK_THAT((*arm)[1], Catch::Matchers::WithinAbs(0.0, 1.0e-12));
    CHECK_THAT((*arm)[2], Catch::Matchers::WithinAbs(0.0, 1.0e-12));
  }
  SECTION("the straight listed first: the same arm")
  {
    feature.portions = {on_straight, on_chord, on_image};
    const auto arm =
        MirrorArmDirection(identification, feature, identification.real_segments);
    REQUIRE(arm);
    CHECK_THAT((*arm)[0], Catch::Matchers::WithinAbs(1.0, 1.0e-12));
    CHECK_THAT((*arm)[1], Catch::Matchers::WithinAbs(0.0, 1.0e-12));
  }
  SECTION("only the chord is real: the chord's far end gives the direction (the pre-fillet "
          "reading)")
  {
    feature.portions = {on_chord, on_image};
    const auto arm =
        MirrorArmDirection(identification, feature, identification.real_segments);
    REQUIRE(arm);
    const double norm = std::hypot(arc_a[0], arc_a[1]);
    CHECK_THAT((*arm)[0], Catch::Matchers::WithinAbs(arc_a[0] / norm, 1.0e-12));
    CHECK_THAT((*arm)[1], Catch::Matchers::WithinAbs(arc_a[1] / norm, 1.0e-12));
  }
  SECTION("no real portion: empty")
  {
    feature.portions = {on_image};
    CHECK_FALSE(MirrorArmDirection(identification, feature, identification.real_segments));
  }
}

}  // namespace palace
