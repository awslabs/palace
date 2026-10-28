// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "fixtures.hpp"
#include "surfaceresponse-fixtures.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iterator>
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
#include "driver.hpp"
#include "fem/fespace.hpp"
#include "fem/gridfunction.hpp"
#include "fem/integrator.hpp"
#include "fem/mesh.hpp"
#include "linalg/vector.hpp"
#include "models/laplaceoperator.hpp"
#include "models/materialoperator.hpp"
#include "models/postoperator.hpp"
#include "models/spaceoperator.hpp"
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
    CHECK(cache["Version"] == 19);
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
  // The Total row of a portion-energy table (uncovered / DomainBoundary /
  // ConductorConsistency raw claims) on the raw field (evaluation 0); 0 when the table was
  // not written (no such portions).
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
    // F-DB-a term of the same identity; + the raw claims of the corner coupons the
    // conductor-consistency gate excludes on this lead: the S7diag R1 term, decision 758).
    const double domain_boundary =
        PortionTableTotal(full_dir / "surface-response-domain-boundary-energy.csv");
    const double conductor_consistency =
        PortionTableTotal(full_dir / "surface-response-conductor-consistency-energy.csv");
    CHECK_THAT(ft, WithinRel(OutsideEnergy(full_dir) + ModelEnergySum(full_dir, 0) +
                                 domain_boundary + conductor_consistency,
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
      WithinRel(OutsideEnergy(corner_dir) + ModelEnergySum(corner_dir, 0) + uncovered_ft +
                    PortionTableTotal(corner_dir /
                                      "surface-response-domain-boundary-energy.csv") +
                    PortionTableTotal(corner_dir /
                                      "surface-response-conductor-consistency-energy.csv"),
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

// Decisions 751 / 752 (the maxwellop diagnosis): the Maxwell path's corrected interface
// energy was E_out + E_fab alone (postoperator.cpp MeasureSurfaceResponseCorrection), so
// the raw within-R energy of every UNMATCHED feature was dropped, while the electrostatic
// path keeps it (decision 394 F2; the case above). The Maxwell analogue of that case on the
// same lead, library and decision-366 Region partition, through a directly constructed
// eigenmode PostOperator fed E = -grad V (any curl-free field exercises the same assembly;
// the corrected field is the raw field, so the self-consistent column equals the fixed
// trace): (a) the Maxwell-configured operator's uncovered portions are the electrostatic
// operator's; (b) the identity: the raw within-R energy of the corner windows (raw minus
// outside minus the three isolated-edge cells) is the uncovered energy, and the corrected
// energy with the model term removed and the covered cells put back is the raw energy
// (nothing within R is dropped); (c) fixed trace / fixed flux / self-consistent = outside +
// model + uncovered; (d) the Maxwell surface-response-uncovered-energy.csv carries the
// ConvexCorner and Total rows for the three evaluations, and the operator's statistics
// record the portions.
TEST_CASE_METHOD(test::SurfaceResponseFiles,
                 "Maxwell uncovered requirements keep their raw energy",
                 "[postoperator][surfaceresponseoperator][cornerarmtrim][3d][maxwell]"
                 "[Serial][Parallel]")
{
#if !defined(MFEM_USE_GSLIB)
  SKIP("SurfaceResponseOperator requires MFEM_USE_GSLIB");
#else
  constexpr double R = 0.2;
  const fs::path mesh_path = temp.temp_dir / "maxwell-uncovered-lead.mesh";
  const auto basis_path = temp.temp_dir / "maxwell-uncovered-basis-points.csv";
  const auto edge_library_path =
      temp.temp_dir / "fabrication-process-maxwell-uncovered-edges.json";
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
    library["Name"] = "unit-test-process-maxwell-uncovered-edges";
    json models = json::array();
    for (auto &model : library["Models"])
    {
      if (model["Topology"] == "IsolatedEdge")
      {
        model["Interfaces"] = {{{"Type", "SA"}, {"Coupon", 1}}};
        model["BasisPoints"] = basis_path.string();
        models.push_back(model);
      }
    }
    REQUIRE(models.size() == 1);
    library["Models"] = std::move(models);
    std::ofstream output(edge_library_path);
    output << library.dump(2) << "\n";
  }
  Mpi::Barrier(Mpi::World());

  // The electrostatic form of the case above (the same mesh, library, Region cells and
  // UnmatchedPolicy) with the DomainBoundary mirror OFF (the Maxwell path runs no mirror:
  // the like-for-like electrostatic comparator is the Mirror-Off operator, decision 743's
  // identity comparator; here the mirror changes nothing, the two Missing corners are far
  // from the cut), then its Maxwell form: an eigenmode problem with the same metal set, the
  // lead and the cut face PEC (every other face natural: a PEC box would make the box edges
  // metal perimeter and the lead's cut a non-manifold junction with Endpoint features).
  json config = IslandConfig();
  config["Model"]["Mesh"] = mesh_path.string();
  config["Model"]["L0"] = 1.0;
  config["Boundaries"]["Ground"]["Attributes"] = {4};
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
  correction["Library"] = edge_library_path.string();
  correction["UnmatchedPolicy"] = "Warn";
  correction["DomainBoundary"] = {{"Mirror", "Off"}};
  IoData electrostatic_iodata(config, false);
  electrostatic_iodata.boundaries.cracked_attributes.insert(9);
  std::vector<std::unique_ptr<Mesh>> electrostatic_meshes;
  electrostatic_meshes.push_back(std::make_unique<Mesh>(MakeLeadMesh(0.0, 0.75)));
  LaplaceOperator laplace_op(electrostatic_iodata, electrostatic_meshes);
  SurfaceResponseOperator electrostatic_response(electrostatic_iodata, laplace_op);

  const fs::path output_dir = temp.temp_dir / "maxwell-uncovered";
  json maxwell_config = config;
  maxwell_config["Problem"]["Type"] = "Eigenmode";
  maxwell_config["Problem"]["Output"] = output_dir.string();
  maxwell_config["Boundaries"]["Ground"]["Attributes"] = {4, 9};
  maxwell_config["Boundaries"].erase("Terminal");
  maxwell_config["Solver"] = {{"Order", 1},
                              {"Eigenmode", {{"Target", 1.0}}},
                              {"SurfaceResponseCorrection",
                               {{"Library", edge_library_path.string()},
                                {"TargetInterfaces", {4}},
                                {"UnmatchedPolicy", "Warn"},
                                {"CorrectionMode", "Both"}}}};
  IoData iodata(maxwell_config, false);
  iodata.boundaries.cracked_attributes.insert(9);
  MPI_Comm comm = Mpi::World();
  MakeOutputFolder(iodata, comm);
  std::vector<std::unique_ptr<Mesh>> meshes;
  meshes.push_back(std::make_unique<Mesh>(MakeLeadMesh(0.0, 0.75)));
  SpaceOperator space_op(iodata, meshes);
  PostOperator<ProblemType::EIGENMODE> post_op(iodata, space_op);
  auto *response = post_op.GetSurfaceResponseOperator();
  REQUIRE(response);
  REQUIRE(response->GetTargetInterfaces() == std::set<int>{4});

  // (a) The Maxwell operator lists the two Missing corners' portions as the electrostatic
  // operator does: 2 corners x 2 arms x 2 mesh segments (h = 0.125) of R along each arm.
  const auto &maxwell_portions = response->GetUncoveredPortions();
  const auto &electrostatic_portions = electrostatic_response.GetUncoveredPortions();
  REQUIRE(electrostatic_portions.size() == 8);
  CHECK(maxwell_portions.size() == electrostatic_portions.size());
  double maxwell_length = 0.0;
  std::map<std::string, int> maxwell_portions_by_type;
  for (std::size_t i = 0;
       i < std::min(maxwell_portions.size(), electrostatic_portions.size()); i++)
  {
    const auto &maxwell = maxwell_portions[i];
    const auto &electrostatic = electrostatic_portions[i];
    maxwell_portions_by_type[maxwell.topology]++;
    CHECK(maxwell.topology == electrostatic.topology);
    CHECK(maxwell.feature == electrostatic.feature);
    CHECK(maxwell.segment == electrostatic.segment);
    double length2 = 0.0;
    for (int d = 0; d < 3; d++)
    {
      CHECK_THAT(maxwell.p0[d], WithinAbs(electrostatic.p0[d], 1.0e-12));
      CHECK_THAT(maxwell.p1[d], WithinAbs(electrostatic.p1[d], 1.0e-12));
      length2 += (maxwell.p1[d] - maxwell.p0[d]) * (maxwell.p1[d] - maxwell.p0[d]);
    }
    maxwell_length += std::sqrt(length2);
  }
  CHECK_THAT(maxwell_length, WithinAbs(4.0 * R, 1.0e-9));
  CHECK(maxwell_portions_by_type == std::map<std::string, int>{{"ConvexCorner", 8}});
  const json statistics = response->GetStatistics();
  const bool records_uncovered =
      statistics.contains("Diagnostics") && statistics["Diagnostics"].contains("Uncovered");
  CHECK(records_uncovered);
  if (records_uncovered)
  {
    CHECK(statistics["Diagnostics"]["Uncovered"]["Count"].get<int>() == 8);
    CHECK(statistics["Diagnostics"]["Uncovered"]["Features"].get<int>() == 2);
    CHECK_THAT(statistics["Diagnostics"]["Uncovered"]["Length"].get<double>(),
               WithinAbs(4.0 * R, 1.0e-9));
  }

  // The raw field E = -grad V of a smooth potential, measured and printed as mode 1; the
  // corrected field is the same field.
  mfem::ParGridFunction potential(&space_op.GetH1Space().Get());
  mfem::FunctionCoefficient potential_coefficient(
      [](const mfem::Vector &x)
      { return std::sin(1.3 * x[0] + 0.4) * std::cosh(0.7 * x[1]) + 0.3 * x[2] * x[2]; });
  potential.ProjectCoefficient(potential_coefficient);
  mfem::ParDiscreteLinearOperator gradient(&space_op.GetH1Space().Get(),
                                           &space_op.GetNDSpace().Get());
  gradient.AddDomainInterpolator(new mfem::GradientInterpolator());
  gradient.Assemble();
  GridFunction field(space_op.GetNDSpace(), true);
  gradient.Mult(potential, field.Real());
  field.Real() *= -1.0;
  ComplexVector e(space_op.GetNDSpace().GetTrueVSize());
  field.Real().GetTrueDofs(e.Real());
  e.Imag() = 0.0;
  field.Real().SetFromTrueDofs(e.Real());
  field.Imag() = 0.0;
  field.Real().ExchangeFaceNbrData();
  field.Imag().ExchangeFaceNbrData();
  ComplexVector b(space_op.GetRTSpace().GetTrueVSize());
  b = 0.0;
  const std::complex<double> omega(1.0, 0.0);
  post_op.SetSurfaceResponseCorrectedField(e, nullptr, omega, 1.0);
  post_op.MeasureAndPrintAll(0, e, b, omega, 0.0, 0.0, 1);

  // The terms of the assembly, each from its own evaluation of the same field.
  const auto &surf_post_op = post_op.GetSurfacePostOp();
  auto RawEnergy = [&](int interface)
  { return surf_post_op.GetInterfaceElectricFieldEnergy(interface, field, nullptr); };
  const auto edge_energies = surf_post_op.GetInterfaceEdgeElectricFieldEnergies(4, field);
  REQUIRE(edge_energies.size() == 1);
  CHECK_THAT(edge_energies[0].distance, WithinAbs(R, 1.0e-12));
  const double outside = edge_energies[0].energy_outside;
  // The uncovered energy of the corner windows from the electrostatic operator's portions
  // (the reference set of (a)); the Maxwell operator's portions give the same energy.
  auto UncoveredEnergy = [&](const auto &portions)
  {
    std::vector<SurfacePostOperator::UncoveredPerimeterPortion> perimeter;
    for (const auto &portion : portions)
    {
      perimeter.push_back({portion.p0, portion.p1, portion.topology, portion.feature});
    }
    const auto energies =
        surf_post_op.GetInterfaceUncoveredEdgeEnergies({4}, field, nullptr, perimeter);
    REQUIRE(energies.size() == 1);
    return energies.at(4).energy;
  };
  const double uncovered = UncoveredEnergy(electrostatic_portions);
  CHECK(uncovered > 0.0);
  CHECK_THAT(UncoveredEnergy(maxwell_portions), WithinRel(uncovered, 1.0e-12));
  const auto maxwell_response = response->GetMaxwellResponse(field, omega);
  const double model_ft = maxwell_response.fabricated_surface_energy.at(4);
  const double model_ff = maxwell_response.fabricated_surface_energy_fixed_flux.at(4);
  const double raw = RawEnergy(4);
  const double cells = RawEnergy(5) + RawEnergy(6) + RawEnergy(7);
  CHECK(raw - outside > cells);
  // (b) The decision-366 partition: the corner windows' raw within-R energy is the
  // uncovered energy.
  CHECK_THAT(uncovered, WithinRel(raw - outside - cells, 1.0e-9));

  Mpi::Barrier(comm);
  if (!Mpi::Root(comm))
  {
    return;
  }
  auto Joule = [&](double energy)
  { return iodata.units.Dimensionalize<Units::ValueType::ENERGY>(energy); };
  auto Trim = [](std::string text)
  {
    const auto begin = text.find_first_not_of(" \t\r");
    const auto end = text.find_last_not_of(" \t\r");
    return begin == std::string::npos ? std::string() : text.substr(begin, end - begin + 1);
  };
  auto ReadTrimmedCsv = [&](const fs::path &path)
  {
    TextCsv csv = ReadTextCsv(path);
    for (auto &name : csv.header)
    {
      name = Trim(name);
    }
    for (auto &row : csv.rows)
    {
      for (auto &cell : row)
      {
        cell = Trim(cell);
      }
    }
    return csv;
  };
  const TextCsv corrected_table = ReadTrimmedCsv(output_dir / "surface-Q-corrected.csv");
  const std::size_t interface_column = corrected_table.Column("interface");
  const auto corrected_row =
      std::find_if(corrected_table.rows.begin(), corrected_table.rows.end(),
                   [&](const auto &row) { return std::stoi(row[interface_column]) == 4; });
  REQUIRE(corrected_row != corrected_table.rows.end());
  const double csv_raw =
      std::stod((*corrected_row)[corrected_table.Column("E_surf raw (J)")]);
  const double csv_ft =
      std::stod((*corrected_row)[corrected_table.Column("E_surf corrected (J)")]);
  const double csv_ff = std::stod(
      (*corrected_row)[corrected_table.Column("E_surf corrected fixed-flux (J)")]);
  const double csv_sc =
      std::stod((*corrected_row)[corrected_table.Column("E_surf self-consistent (J)")]);
  CHECK_THAT(csv_raw, WithinRel(Joule(raw), 1.0e-10));
  // (c) Fixed trace / fixed flux / self-consistent: outside + model + uncovered.
  CHECK_THAT(csv_ft, WithinRel(Joule(outside + model_ft + uncovered), 1.0e-10));
  CHECK_THAT(csv_ff, WithinRel(Joule(outside + model_ff + uncovered), 1.0e-10));
  REQUIRE(std::isfinite(csv_sc));
  CHECK_THAT(csv_sc, WithinRel(csv_ft, 1.0e-10));
  // (b) The identity: the corrected energy with the model term removed and the covered
  // cells put back is the raw energy.
  CHECK_THAT(csv_ft - Joule(model_ft) + Joule(cells), WithinRel(Joule(raw), 1.0e-10));

  // (d) The Maxwell uncovered table: (ConvexCorner, Total) x (fixed trace, fixed flux,
  // self-consistent) for mode 1, excitation 0.
  REQUIRE(fs::is_regular_file(output_dir / "surface-response-uncovered-energy.csv"));
  const TextCsv uncovered_table =
      ReadTrimmedCsv(output_dir / "surface-response-uncovered-energy.csv");
  const std::size_t mode_column = uncovered_table.Column("m");
  const std::size_t excitation_column = uncovered_table.Column("exc");
  const std::size_t evaluation_column = uncovered_table.Column("evaluation");
  const std::size_t type_column = uncovered_table.Column("type");
  const std::size_t portions_column = uncovered_table.Column("portions");
  const std::size_t length_column = uncovered_table.Column("length (m)");
  const std::size_t clipped_column =
      uncovered_table.Column("clipped by spatial support (m)");
  const std::size_t energy_column = uncovered_table.Column("uncovered raw energy[4] (J)");
  const std::size_t share_column = uncovered_table.Column("uncovered share[4]");
  REQUIRE(uncovered_table.rows.size() == 6);
  std::map<std::pair<int, std::string>, std::vector<std::string>> by_key;
  for (const auto &row : uncovered_table.rows)
  {
    CHECK(std::stoi(row[mode_column]) == 1);
    CHECK(std::stoi(row[excitation_column]) == 0);
    by_key[{std::stoi(row[evaluation_column]), row[type_column]}] = row;
  }
  REQUIRE(by_key.size() == 6);
  const std::array<double, 3> corrected_by_evaluation = {csv_ft, csv_ff, csv_sc};
  for (int evaluation = 0; evaluation < 3; evaluation++)
  {
    for (const char *type : {"ConvexCorner", "Total"})
    {
      REQUIRE(by_key.count({evaluation, type}) == 1);
      const auto &row = by_key.at({evaluation, type});
      CHECK(std::stoi(row[portions_column]) == 8);
      CHECK_THAT(std::stod(row[length_column]), WithinAbs(4.0 * R, 1.0e-9));
      CHECK(std::stod(row[clipped_column]) == 0.0);
      CHECK_THAT(std::stod(row[energy_column]), WithinRel(Joule(uncovered), 1.0e-10));
      CHECK_THAT(std::stod(row[share_column]),
                 WithinRel(Joule(uncovered) / corrected_by_evaluation[evaluation], 1.0e-9));
    }
  }
#endif
}

// M-5 (decisions 845 O-C / 846): the Maxwell path writes surface-response-model-energy.csv
// with the electrostatic driver's column set exactly and in order as a strict prefix, then
// the Maxwell keys "m" / "exc". On the lead fixture above (one isolated-edge model, three
// Region cells, two Missing corners; the ES comparator = the Mirror-Off operator of the
// same mesh and library) the eigenmode PostOperator fed E = -grad V writes, for mode 1, the
// rows (ev 0 fixed trace, ev 1 fixed flux, ev 2 self-consistent = the same field) whose
// values are the electrostatic operator's ModelContribution of V to 1e-10 relative (the
// line integrals of the exact gradient of a p1 potential are its differences), with
// `source` = the mode index, 0 for the three non-target interfaces, and the ev-0 row's
// model term equal to surface-Q-corrected.csv's corrected minus outside minus uncovered
// energy.
TEST_CASE_METHOD(test::SurfaceResponseFiles,
                 "Maxwell model energies are written as the electrostatic rows",
                 "[postoperator][surfaceresponseoperator][cornerarmtrim][3d][maxwell]"
                 "[Serial][Parallel]")
{
#if !defined(MFEM_USE_GSLIB)
  SKIP("SurfaceResponseOperator requires MFEM_USE_GSLIB");
#else
  constexpr double R = 0.2;
  const fs::path mesh_path = temp.temp_dir / "maxwell-model-energy-lead.mesh";
  const auto basis_path = temp.temp_dir / "maxwell-model-energy-basis-points.csv";
  const auto edge_library_path =
      temp.temp_dir / "fabrication-process-maxwell-model-energy-edges.json";
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
    library["Name"] = "unit-test-process-maxwell-model-energy-edges";
    json models = json::array();
    for (auto &model : library["Models"])
    {
      if (model["Topology"] == "IsolatedEdge")
      {
        model["Interfaces"] = {{{"Type", "SA"}, {"Coupon", 1}}};
        model["BasisPoints"] = basis_path.string();
        models.push_back(model);
      }
    }
    REQUIRE(models.size() == 1);
    library["Models"] = std::move(models);
    std::ofstream output(edge_library_path);
    output << library.dump(2) << "\n";
  }
  Mpi::Barrier(Mpi::World());

  json config = IslandConfig();
  config["Model"]["Mesh"] = mesh_path.string();
  config["Model"]["L0"] = 1.0;
  config["Boundaries"]["Ground"]["Attributes"] = {4};
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
  correction["Library"] = edge_library_path.string();
  correction["UnmatchedPolicy"] = "Warn";
  correction["DomainBoundary"] = {{"Mirror", "Off"}};
  IoData electrostatic_iodata(config, false);
  electrostatic_iodata.boundaries.cracked_attributes.insert(9);
  std::vector<std::unique_ptr<Mesh>> electrostatic_meshes;
  electrostatic_meshes.push_back(std::make_unique<Mesh>(MakeLeadMesh(0.0, 0.75)));
  LaplaceOperator laplace_op(electrostatic_iodata, electrostatic_meshes);
  SurfaceResponseOperator electrostatic_response(electrostatic_iodata, laplace_op);

  const fs::path output_dir = temp.temp_dir / "maxwell-model-energy";
  json maxwell_config = config;
  maxwell_config["Problem"]["Type"] = "Eigenmode";
  maxwell_config["Problem"]["Output"] = output_dir.string();
  maxwell_config["Boundaries"]["Ground"]["Attributes"] = {4, 9};
  maxwell_config["Boundaries"].erase("Terminal");
  maxwell_config["Solver"] = {{"Order", 1},
                              {"Eigenmode", {{"Target", 1.0}}},
                              {"SurfaceResponseCorrection",
                               {{"Library", edge_library_path.string()},
                                {"TargetInterfaces", {4}},
                                {"UnmatchedPolicy", "Warn"},
                                {"CorrectionMode", "Both"}}}};
  IoData iodata(maxwell_config, false);
  iodata.boundaries.cracked_attributes.insert(9);
  MPI_Comm comm = Mpi::World();
  MakeOutputFolder(iodata, comm);
  std::vector<std::unique_ptr<Mesh>> meshes;
  meshes.push_back(std::make_unique<Mesh>(MakeLeadMesh(0.0, 0.75)));
  SpaceOperator space_op(iodata, meshes);
  PostOperator<ProblemType::EIGENMODE> post_op(iodata, space_op);
  auto *response = post_op.GetSurfaceResponseOperator();
  REQUIRE(response);
  REQUIRE(response->GetPatchCount() == electrostatic_response.GetPatchCount());

  // The same smooth potential on both meshes (the same serial mesh and partition): the
  // electrostatic operator reads V, the Maxwell operator E = -grad V. V vanishes on the
  // lead's plane y = 0.5 and its field is normal to it (a PEC carries no tangential field;
  // the Maxwell trace gauges the conductor at 0 where the electrostatic trace samples V on
  // it), so the two traces agree coefficient by coefficient.
  mfem::FunctionCoefficient potential_coefficient(
      [](const mfem::Vector &x) { return std::sinh(0.7 * (x[1] - 0.5)); });
  mfem::ParGridFunction electrostatic_potential(&laplace_op.GetH1Space().Get());
  electrostatic_potential.ProjectCoefficient(potential_coefficient);
  Vector electrostatic_potential_true;
  electrostatic_potential.GetTrueDofs(electrostatic_potential_true);
  const auto electrostatic =
      electrostatic_response.GetElectrostaticResponse(electrostatic_potential_true);
  REQUIRE(electrostatic.model_contributions.size() == 1);
  const auto &reference = electrostatic.model_contributions.front();
  CHECK(reference.patch_count > 0.0);
  CHECK(reference.fabricated_surface_energy.at(4) > 0.0);

  mfem::ParGridFunction potential(&space_op.GetH1Space().Get());
  potential.ProjectCoefficient(potential_coefficient);
  mfem::ParDiscreteLinearOperator gradient(&space_op.GetH1Space().Get(),
                                           &space_op.GetNDSpace().Get());
  gradient.AddDomainInterpolator(new mfem::GradientInterpolator());
  gradient.Assemble();
  GridFunction field(space_op.GetNDSpace(), true);
  gradient.Mult(potential, field.Real());
  field.Real() *= -1.0;
  ComplexVector e(space_op.GetNDSpace().GetTrueVSize());
  field.Real().GetTrueDofs(e.Real());
  e.Imag() = 0.0;
  field.Real().SetFromTrueDofs(e.Real());
  field.Imag() = 0.0;
  field.Real().ExchangeFaceNbrData();
  field.Imag().ExchangeFaceNbrData();
  ComplexVector b(space_op.GetRTSpace().GetTrueVSize());
  b = 0.0;
  const std::complex<double> omega(1.0, 0.0);
  post_op.SetSurfaceResponseCorrectedField(e, nullptr, omega, 1.0);
  post_op.MeasureAndPrintAll(0, e, b, omega, 0.0, 0.0, 1);

  // The operator's own per-model entry of the same field equals the electrostatic one
  // (every rank; the entries are replicated).
  const auto maxwell = response->GetMaxwellResponse(field, omega);
  REQUIRE(maxwell.model_contributions.size() == 1);
  const auto &entry = maxwell.model_contributions.front();
  CHECK(entry.model == reference.model);
  CHECK(entry.patch_count == reference.patch_count);
  CHECK_THAT(entry.patch_weight, WithinRel(reference.patch_weight, 1.0e-12));
  CHECK_THAT(entry.domain_correction, WithinRel(reference.domain_correction, 1.0e-10));
  CHECK_THAT(entry.domain_correction_fixed_flux,
             WithinRel(reference.domain_correction_fixed_flux, 1.0e-10));
  CHECK_THAT(entry.fabricated_surface_energy.at(4),
             WithinRel(reference.fabricated_surface_energy.at(4), 1.0e-10));
  CHECK_THAT(entry.fabricated_surface_energy_fixed_flux.at(4),
             WithinRel(reference.fabricated_surface_energy_fixed_flux.at(4), 1.0e-10));
  const auto &surf_post_op = post_op.GetSurfacePostOp();
  const auto edge_energies = surf_post_op.GetInterfaceEdgeElectricFieldEnergies(4, field);
  REQUIRE(edge_energies.size() == 1);
  const double outside = edge_energies[0].energy_outside;
  std::vector<SurfacePostOperator::UncoveredPerimeterPortion> perimeter;
  for (const auto &portion : response->GetUncoveredPortions())
  {
    perimeter.push_back({portion.p0, portion.p1, portion.topology, portion.feature});
  }
  const auto uncovered_energies =
      surf_post_op.GetInterfaceUncoveredEdgeEnergies({4}, field, nullptr, perimeter);
  REQUIRE(uncovered_energies.size() == 1);
  const double uncovered = uncovered_energies.at(4).energy;

  Mpi::Barrier(comm);
  if (!Mpi::Root(comm))
  {
    return;
  }
  auto Joule = [&](double energy)
  { return iodata.units.Dimensionalize<Units::ValueType::ENERGY>(energy); };
  auto Trim = [](std::string text)
  {
    const auto begin = text.find_first_not_of(" \t\r");
    const auto end = text.find_last_not_of(" \t\r");
    return begin == std::string::npos ? std::string() : text.substr(begin, end - begin + 1);
  };
  auto ReadTrimmedCsv = [&](const fs::path &path)
  {
    TextCsv csv = ReadTextCsv(path);
    for (auto &name : csv.header)
    {
      name = Trim(name);
    }
    for (auto &row : csv.rows)
    {
      for (auto &cell : row)
      {
        cell = Trim(cell);
      }
    }
    return csv;
  };
  REQUIRE(fs::is_regular_file(output_dir / "surface-response-model-energy.csv"));
  const TextCsv table = ReadTrimmedCsv(output_dir / "surface-response-model-energy.csv");
  // The electrostatic driver's column set, in order, as a strict prefix (every
  // postprocessing interface: the target 4 and the three Region cells 5-7), then the keys.
  const std::vector<std::string> expected_header = {"source",
                                                    "evaluation",
                                                    "model",
                                                    "patch count",
                                                    "patch weight",
                                                    "domain correction (J)",
                                                    "fabricated surface energy[4] (J)",
                                                    "fabricated surface energy[5] (J)",
                                                    "fabricated surface energy[6] (J)",
                                                    "fabricated surface energy[7] (J)",
                                                    "m",
                                                    "exc"};
  CHECK(table.header == expected_header);
  REQUIRE(table.rows.size() == 3);
  const TextCsv corrected_table = ReadTrimmedCsv(output_dir / "surface-Q-corrected.csv");
  const std::size_t interface_column = corrected_table.Column("interface");
  const auto corrected_row =
      std::find_if(corrected_table.rows.begin(), corrected_table.rows.end(),
                   [&](const auto &row) { return std::stoi(row[interface_column]) == 4; });
  REQUIRE(corrected_row != corrected_table.rows.end());
  const double csv_ft =
      std::stod((*corrected_row)[corrected_table.Column("E_surf corrected (J)")]);
  for (int evaluation = 0; evaluation < 3; evaluation++)
  {
    const auto &row = table.rows[evaluation];
    INFO("evaluation " << evaluation);
    CHECK(std::stoi(row[table.Column("evaluation")]) == evaluation);
    CHECK(std::stoi(row[table.Column("source")]) == 1);
    CHECK(std::stoi(row[table.Column("m")]) == 1);
    CHECK(std::stoi(row[table.Column("exc")]) == 0);
    CHECK(std::stoi(row[table.Column("model")]) == reference.model);
    CHECK(std::stod(row[table.Column("patch count")]) == reference.patch_count);
    CHECK_THAT(std::stod(row[table.Column("patch weight")]),
               WithinRel(reference.patch_weight, 1.0e-10));
    const bool fixed_flux = evaluation == 1;
    CHECK_THAT(std::stod(row[table.Column("domain correction (J)")]),
               WithinRel(Joule(fixed_flux ? reference.domain_correction_fixed_flux
                                          : reference.domain_correction),
                         1.0e-10));
    CHECK_THAT(
        std::stod(row[table.Column("fabricated surface energy[4] (J)")]),
        WithinRel(Joule(fixed_flux ? reference.fabricated_surface_energy_fixed_flux.at(4)
                                   : reference.fabricated_surface_energy.at(4)),
                  1.0e-10));
    for (int cell = 5; cell <= 7; cell++)
    {
      CHECK(
          std::stod(
              row[table.Column(fmt::format("fabricated surface energy[{}] (J)", cell))]) ==
          0.0);
    }
  }
  // The ev-0 model term is the coupon term of surface-Q-corrected.csv: corrected minus the
  // outside and the uncovered raw energy.
  CHECK_THAT(std::stod(table.rows[0][table.Column("fabricated surface energy[4] (J)")]),
             WithinRel(csv_ft - Joule(outside + uncovered), 1.0e-9));
#endif
}

// PatchEnergy on the Maxwell path (decisions 913 / 920 / 922 B; the 851 (C) per-model RAW
// partition formalised). On the lead fixture of the case above with "PatchEnergy": true:
// (1) the applied patches' raw claims (GetModelRawClaimPortions) carry one label "m<idx>"
// for the isolated-edge model on both the electrostatic and the Maxwell operator, with the
// model's patch count; the surface postoperator partitions the within-R energy of any
// field EXACTLY over the labels plus the unclaimed remainder: Σ labels + unclaimed = E_raw
// − E_out to 1e-9 on both postoperators (the two operators' portions give the same
// partition); the unclaimed remainder equals the uncovered requirements' energy (the two
// Missing corners' portions); the Maxwell path writes surface-response-model-raw-energy.csv
// / -groups.csv with those values. (2) The Maxwell per-patch export: the patch table at
// construction, surface-response-patch-energy.csv (the electrostatic columns as a strict
// prefix, then "m" / "exc" / "max loop residual") whose per-model sums are the
// surface-response-model-energy.csv rows, and surface-response-patch-loops.csv (920) whose
// contour-loop maximum is the confidence table's "loop residual" (the 0.05-gate quantity)
// and whose response-weighted RMS, recomputed from the per-path rows and the per-patch
// response energies, is its "response-weighted loop residual". (3) Without PatchEnergy
// (the default) none of the files exists, the claims are empty and
// surface-response-model-energy.csv is byte-identical. On main the Maxwell configuration
// refuses "PatchEnergy" (configfile.cpp) and the operator has no per-patch records.
TEST_CASE_METHOD(test::SurfaceResponseFiles,
                 "Maxwell per-patch export and per-model raw partition (decision 922 B)",
                 "[postoperator][surfaceresponseoperator][cornerarmtrim][3d][maxwell]"
                 "[patchenergy][Serial][Parallel]")
{
#if !defined(MFEM_USE_GSLIB)
  SKIP("SurfaceResponseOperator requires MFEM_USE_GSLIB");
#else
  constexpr double R = 0.2;
  const fs::path mesh_path = temp.temp_dir / "model-raw-lead.mesh";
  const auto basis_path = temp.temp_dir / "model-raw-basis-points.csv";
  const auto edge_library_path = temp.temp_dir / "fabrication-process-model-raw-edges.json";
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
    library["Name"] = "unit-test-process-model-raw-edges";
    json models = json::array();
    for (auto &model : library["Models"])
    {
      if (model["Topology"] == "IsolatedEdge")
      {
        model["Interfaces"] = {{{"Type", "SA"}, {"Coupon", 1}}};
        model["BasisPoints"] = basis_path.string();
        models.push_back(model);
      }
    }
    REQUIRE(models.size() == 1);
    library["Models"] = std::move(models);
    std::ofstream output(edge_library_path);
    output << library.dump(2) << "\n";
  }
  Mpi::Barrier(Mpi::World());

  json config = IslandConfig();
  config["Model"]["Mesh"] = mesh_path.string();
  config["Model"]["L0"] = 1.0;
  config["Boundaries"]["Ground"]["Attributes"] = {4};
  auto &correction = config["Solver"]["Electrostatic"]["ResponseCorrection"];
  correction.erase("PatchConstruction");
  correction["Library"] = edge_library_path.string();
  correction["UnmatchedPolicy"] = "Warn";
  correction["DomainBoundary"] = {{"Mirror", "Off"}};
  correction["PatchEnergy"] = true;
  const fs::path electrostatic_output_dir = temp.temp_dir / "model-raw-electrostatic";
  config["Problem"]["Output"] = electrostatic_output_dir.string();
  IoData electrostatic_iodata(config, false);
  electrostatic_iodata.boundaries.cracked_attributes.insert(9);
  MPI_Comm comm = Mpi::World();
  MakeOutputFolder(electrostatic_iodata, comm);
  std::vector<std::unique_ptr<Mesh>> electrostatic_meshes;
  electrostatic_meshes.push_back(std::make_unique<Mesh>(MakeLeadMesh(0.0, 0.75)));
  LaplaceOperator laplace_op(electrostatic_iodata, electrostatic_meshes);
  SurfaceResponseOperator electrostatic_response(electrostatic_iodata, laplace_op);
  PostOperator<ProblemType::ELECTROSTATIC> electrostatic_post_op(electrostatic_iodata,
                                                                 laplace_op);

  const fs::path output_dir = temp.temp_dir / "model-raw-maxwell";
  json maxwell_config = config;
  maxwell_config["Problem"]["Type"] = "Eigenmode";
  maxwell_config["Problem"]["Output"] = output_dir.string();
  maxwell_config["Boundaries"]["Ground"]["Attributes"] = {4, 9};
  maxwell_config["Boundaries"].erase("Terminal");
  maxwell_config["Solver"] = {{"Order", 1},
                              {"Eigenmode", {{"Target", 1.0}}},
                              {"SurfaceResponseCorrection",
                               {{"Library", edge_library_path.string()},
                                {"TargetInterfaces", {4}},
                                {"UnmatchedPolicy", "Warn"},
                                {"CorrectionMode", "PostprocessOnly"},
                                {"PatchEnergy", true}}}};
  IoData iodata(maxwell_config, false);
  iodata.boundaries.cracked_attributes.insert(9);
  MakeOutputFolder(iodata, comm);
  std::vector<std::unique_ptr<Mesh>> meshes;
  meshes.push_back(std::make_unique<Mesh>(MakeLeadMesh(0.0, 0.75)));
  SpaceOperator space_op(iodata, meshes);
  PostOperator<ProblemType::EIGENMODE> post_op(iodata, space_op);
  auto *response = post_op.GetSurfaceResponseOperator();
  REQUIRE(response);
  REQUIRE(response->GetPatchCount() == electrostatic_response.GetPatchCount());

  // The claims: one label for the one applied model on both operators, every applied patch
  // of the model on it, the same portion count (the same placed patches).
  const auto &electrostatic_claims = electrostatic_response.GetModelRawClaimPortions();
  const auto &maxwell_claims = response->GetModelRawClaimPortions();
  REQUIRE(!electrostatic_claims.portions.empty());
  REQUIRE(electrostatic_claims.portions.size() == maxwell_claims.portions.size());
  REQUIRE(electrostatic_claims.members.size() == 1);
  REQUIRE(maxwell_claims.members.size() == 1);
  const auto &member = electrostatic_claims.members.front();
  const std::string label = fmt::format("m{}", member.model);
  CHECK(member.label == label);
  CHECK(maxwell_claims.members.front().label == label);
  CHECK(maxwell_claims.members.front().model == member.model);
  CHECK(maxwell_claims.members.front().patches == member.patches);
  CHECK_THAT(maxwell_claims.members.front().weight_length,
             WithinRel(member.weight_length, 1.0e-12));
  CHECK_THAT(maxwell_claims.length, WithinRel(electrostatic_claims.length, 1.0e-12));
  for (const auto &portion : electrostatic_claims.portions)
  {
    CHECK(portion.topology == label);
  }
  auto Perimeter = [](const auto &portions)
  {
    std::vector<SurfacePostOperator::UncoveredPerimeterPortion> perimeter;
    for (const auto &portion : portions)
    {
      perimeter.push_back({portion.p0, portion.p1, portion.topology, portion.feature});
    }
    return perimeter;
  };

  mfem::FunctionCoefficient potential_coefficient(
      [](const mfem::Vector &x) { return std::sinh(0.7 * (x[1] - 0.5)); });
  auto Gradient = [&](auto &h1_fespace, auto &nd_fespace)
  {
    mfem::ParGridFunction potential(&h1_fespace.Get());
    potential.ProjectCoefficient(potential_coefficient);
    mfem::ParDiscreteLinearOperator gradient(&h1_fespace.Get(), &nd_fespace.Get());
    gradient.AddDomainInterpolator(new mfem::GradientInterpolator());
    gradient.Assemble();
    auto field = std::make_unique<GridFunction>(nd_fespace, false);
    gradient.Mult(potential, field->Real());
    field->Real() *= -1.0;
    Vector e(nd_fespace.GetTrueVSize());
    field->Real().GetTrueDofs(e);
    field->Real().SetFromTrueDofs(e);
    field->Real().ExchangeFaceNbrData();
    return field;
  };
  // The exact partition identity on both postoperators with both operators' claims:
  // Σ labels + unclaimed = E_raw − E_out to 1e-9; the two claim sets partition alike.
  auto CheckPartition = [&](const SurfacePostOperator &surf_post_op,
                            const GridFunction &field, const char *path)
  {
    INFO(path);
    const double raw = surf_post_op.GetInterfaceElectricFieldEnergy(4, field);
    const auto edge_energies = surf_post_op.GetInterfaceEdgeElectricFieldEnergies(4, field);
    REQUIRE(edge_energies.size() == 1);
    const double within = raw - edge_energies[0].energy_outside;
    REQUIRE(within > 0.0);
    const auto uncovered = surf_post_op.GetInterfaceUncoveredEdgeEnergies(
        {4}, field, nullptr, Perimeter(response->GetUncoveredPortions()));
    std::map<std::string, double> reference;
    for (const auto *claims : {&electrostatic_claims, &maxwell_claims})
    {
      const auto partition = surf_post_op.GetInterfacePartitionedEdgeEnergies(
          {4}, field, nullptr, Perimeter(claims->portions));
      REQUIRE(partition.size() == 1);
      const auto &by_label = partition.at(4);
      REQUIRE(by_label.size() == 2);
      REQUIRE(by_label.count(label) == 1);
      REQUIRE(by_label.count("unclaimed") == 1);
      CHECK(by_label.at(label) > 0.0);
      double sum = 0.0;
      for (const auto &[key, energy] : by_label)
      {
        (void)key;
        sum += energy;
      }
      CHECK_THAT(sum, WithinRel(within, 1.0e-9));
      // The unclaimed remainder is the uncovered requirements' energy (the claims and the
      // uncovered portions tile the perimeter within R).
      CHECK_THAT(by_label.at("unclaimed"), WithinRel(uncovered.at(4).energy, 1.0e-9));
      if (reference.empty())
      {
        reference = by_label;
      }
      else
      {
        for (const auto &[key, energy] : by_label)
        {
          CHECK_THAT(energy, WithinRel(reference.at(key), 1.0e-12));
        }
      }
    }
    return reference.at(label);
  };
  const auto electrostatic_field =
      Gradient(laplace_op.GetH1Space(), laplace_op.GetNDSpace());
  CheckPartition(electrostatic_post_op.GetSurfacePostOp(), *electrostatic_field,
                 "electrostatic");
  const auto maxwell_field = Gradient(space_op.GetH1Space(), space_op.GetNDSpace());
  const double model_raw =
      CheckPartition(post_op.GetSurfacePostOp(), *maxwell_field, "maxwell");

  // The Maxwell path's tables of the same field.
  GridFunction field(space_op.GetNDSpace(), true);
  field.Real() = maxwell_field->Real();
  field.Imag() = 0.0;
  ComplexVector e(space_op.GetNDSpace().GetTrueVSize());
  field.Real().GetTrueDofs(e.Real());
  e.Imag() = 0.0;
  ComplexVector b(space_op.GetRTSpace().GetTrueVSize());
  b = 0.0;
  const std::complex<double> omega(1.0, 0.0);
  post_op.MeasureAndPrintAll(0, e, b, omega, 0.0, 0.0, 1);
  // A second print (mode 2) with a ROTATIONAL field: the line-integral trace lift of a
  // non-conservative field closes its loops with O(1) residuals, so that the per-path
  // records (decision 920) are read on a non-trivial case (the gradient field above
  // closes every loop to roundoff).
  {
    mfem::VectorFunctionCoefficient rotational(3,
                                               [](const mfem::Vector &x, mfem::Vector &v)
                                               {
                                                 v.SetSize(3);
                                                 v(0) = 0.4 * (x[1] - 0.5) + 0.1 * x[2];
                                                 v(1) = -0.4 * (x[0] - 0.5) + 0.3;
                                                 v(2) = 0.2 * x[0] * (x[1] - 0.5);
                                               });
    GridFunction rotational_field(space_op.GetNDSpace(), true);
    rotational_field.Real().ProjectCoefficient(rotational);
    rotational_field.Imag() = 0.0;
    ComplexVector e2(space_op.GetNDSpace().GetTrueVSize());
    rotational_field.Real().GetTrueDofs(e2.Real());
    e2.Imag() = 0.0;
    post_op.MeasureAndPrintAll(1, e2, b, omega, 0.0, 0.0, 2);
  }
  Mpi::Barrier(comm);
  if (!Mpi::Root(comm))
  {
    return;
  }
  auto Joule = [&](double energy)
  { return iodata.units.Dimensionalize<Units::ValueType::ENERGY>(energy); };
  auto Trim = [](std::string text)
  {
    const auto begin = text.find_first_not_of(" \t\r");
    const auto end = text.find_last_not_of(" \t\r");
    return begin == std::string::npos ? std::string() : text.substr(begin, end - begin + 1);
  };
  // The integer columns print in the tables' numeric form ("1.60e+01"): read as numbers.
  auto Int = [](const std::string &field)
  { return static_cast<int>(std::lround(std::stod(field))); };
  auto ReadTrimmedCsv = [&](const fs::path &path)
  {
    TextCsv csv = ReadTextCsv(path);
    for (auto &name : csv.header)
    {
      name = Trim(name);
    }
    for (auto &row : csv.rows)
    {
      for (auto &cell : row)
      {
        cell = Trim(cell);
      }
    }
    return csv;
  };
  REQUIRE(fs::is_regular_file(output_dir / "surface-response-model-raw-energy.csv"));
  const TextCsv table =
      ReadTrimmedCsv(output_dir / "surface-response-model-raw-energy.csv");
  const std::vector<std::string> expected_header = {
      "source", "evaluation", "label", "raw energy[4] (J)", "m", "exc"};
  CHECK(table.header == expected_header);
  REQUIRE(table.rows.size() == 4);  // 2 labels x 2 modes
  std::map<std::string, std::vector<std::string>> by_label;
  for (const auto &row : table.rows)
  {
    CHECK(Int(row[table.Column("evaluation")]) == 0);
    CHECK(Int(row[table.Column("exc")]) == 0);
    CHECK(Int(row[table.Column("source")]) == Int(row[table.Column("m")]));
    if (Int(row[table.Column("m")]) == 1)
    {
      by_label[row[table.Column("label")]] = row;
    }
  }
  REQUIRE(by_label.count(label) == 1);
  REQUIRE(by_label.count("unclaimed") == 1);
  CHECK_THAT(std::stod(by_label.at(label)[table.Column("raw energy[4] (J)")]),
             WithinRel(Joule(model_raw), 1.0e-10));
  REQUIRE(fs::is_regular_file(output_dir / "surface-response-model-raw-groups.csv"));
  const TextCsv groups =
      ReadTrimmedCsv(output_dir / "surface-response-model-raw-groups.csv");
  CHECK(groups.header ==
        std::vector<std::string>{"label", "model", "patches", "weight length (m)"});
  REQUIRE(groups.rows.size() == 1);
  CHECK(groups.rows[0][groups.Column("label")] == label);
  CHECK(Int(groups.rows[0][groups.Column("model")]) == member.model);
  CHECK(Int(groups.rows[0][groups.Column("patches")]) == member.patches);

  // (2) The per-patch export of the same print.
  REQUIRE(fs::is_regular_file(output_dir / "surface-response-patches.csv"));
  const TextCsv patch_table = ReadTrimmedCsv(output_dir / "surface-response-patches.csv");
  CHECK(patch_table.header ==
        std::vector<std::string>{"model", "origin x (m)", "origin y (m)", "origin z (m)",
                                 "axis u x", "axis u y", "axis u z", "axis v x", "axis v y",
                                 "axis v z", "axis w x", "axis w y", "axis w z", "weight"});
  CHECK(patch_table.rows.size() == static_cast<std::size_t>(response->GetPatchCount()));
  REQUIRE(fs::is_regular_file(output_dir / "surface-response-patch-energy.csv"));
  const TextCsv patch_energy =
      ReadTrimmedCsv(output_dir / "surface-response-patch-energy.csv");
  CHECK(patch_energy.header ==
        std::vector<std::string>{
            "source", "evaluation", "patch", "model", "feature", "weight", "cell begin (m)",
            "cell end (m)", "origin x (m)", "origin y (m)", "origin z (m)",
            "domain correction (J)", "fabricated surface energy[4] (J)", "m", "exc",
            "max loop residual", "loop share", "loop fallback", "loop weight (J)"});
  // Evaluation 0 and 1 for every applied patch and mode (PostprocessOnly: no evaluation 2).
  REQUIRE(patch_energy.rows.size() ==
          4 * static_cast<std::size_t>(response->GetPatchCount()));
  const TextCsv model_energy =
      ReadTrimmedCsv(output_dir / "surface-response-model-energy.csv");
  std::map<std::tuple<int, int, int>, std::array<double, 4>> model_sums;  // (m, model, ev)
  double max_patch_loop = 0.0;                                            // mode 2
  std::map<int, double> patch_response_energy;  // mode 2: the weight of the loop RMS
  for (const auto &row : patch_energy.rows)
  {
    const int mode = Int(row[patch_energy.Column("m")]);
    CHECK(Int(row[patch_energy.Column("source")]) == mode);
    CHECK(Int(row[patch_energy.Column("exc")]) == 0);
    const int ev = Int(row[patch_energy.Column("evaluation")]);
    const int model = Int(row[patch_energy.Column("model")]);
    const int patch = Int(row[patch_energy.Column("patch")]);
    auto &sums = model_sums[{mode, model, ev}];
    sums[0] += 1.0;
    sums[1] += std::stod(row[patch_energy.Column("weight")]);
    sums[2] += std::stod(row[patch_energy.Column("domain correction (J)")]);
    const double energy =
        std::stod(row[patch_energy.Column("fabricated surface energy[4] (J)")]);
    sums[3] += energy;
    if (mode == 2)
    {
      // The loop-RMS weight is the patch's response energy on the CHAIN trace (decision 955
      // D3: the detector is independent of the least-squares lift), exported as "loop
      // weight (J)"; on the rotational field it differs from the fitted energy above.
      patch_response_energy[patch] = std::stod(row[patch_energy.Column("loop weight (J)")]);
      max_patch_loop = std::max(max_patch_loop,
                                std::stod(row[patch_energy.Column("max loop residual")]));
    }
  }
  REQUIRE(model_energy.rows.size() == 4);
  for (const auto &row : model_energy.rows)
  {
    const int mode = Int(row[model_energy.Column("m")]);
    const int ev = Int(row[model_energy.Column("evaluation")]);
    const int model = Int(row[model_energy.Column("model")]);
    REQUIRE(model_sums.count({mode, model, ev}) == 1);
    const auto &sums = model_sums.at({mode, model, ev});
    CHECK_THAT(sums[0],
               WithinRel(std::stod(row[model_energy.Column("patch count")]), 1.0e-12));
    CHECK_THAT(sums[1],
               WithinRel(std::stod(row[model_energy.Column("patch weight")]), 1.0e-9));
    CHECK_THAT(
        sums[2],
        WithinRel(std::stod(row[model_energy.Column("domain correction (J)")]), 1.0e-9));
    CHECK_THAT(
        sums[3],
        WithinRel(std::stod(row[model_energy.Column("fabricated surface energy[4] (J)")]),
                  1.0e-9));
  }
  // No spatial patch on this fixture: no traces file.
  CHECK(!fs::exists(output_dir / "surface-response-traces.csv"));
  // (920) The per-path loop residuals: every row a contour loop or a hat of an applied
  // patch; the contour-loop maximum = the gate quantity; the response-weighted RMS
  // recomputed with the operator's weights (a patch's response energy split over its paths
  // by the squared loop scale) = the confidence table's.
  REQUIRE(fs::is_regular_file(output_dir / "surface-response-patch-loops.csv"));
  const TextCsv loops = ReadTrimmedCsv(output_dir / "surface-response-patch-loops.csv");
  CHECK(loops.header == std::vector<std::string>{"source", "patch", "model", "path", "kind",
                                                 "loop residual", "loop scale (V)", "m",
                                                 "exc"});
  REQUIRE(!loops.rows.empty());
  const TextCsv confidence = ReadTrimmedCsv(output_dir / "surface-response-confidence.csv");
  REQUIRE(confidence.rows.size() == 2);
  REQUIRE(Int(confidence.rows[1][confidence.Column("m")]) == 2);
  const double gate_loop_residual =
      std::stod(confidence.rows[1][confidence.Column("loop residual")]);
  const double gate_weighted =
      std::stod(confidence.rows[1][confidence.Column("response-weighted loop residual")]);
  REQUIRE(gate_loop_residual > 1.0e-3);  // the rotational field: a real residual
  double max_loop = 0.0;
  std::map<int, std::vector<std::pair<double, double>>> paths_by_patch;
  for (const auto &row : loops.rows)
  {
    CHECK(Int(row[loops.Column("source")]) == Int(row[loops.Column("m")]));
    if (Int(row[loops.Column("m")]) != 2)
    {
      continue;
    }
    const int patch = Int(row[loops.Column("patch")]);
    const int kind = Int(row[loops.Column("kind")]);
    CHECK((kind == 0 || kind == 1));
    CHECK(patch_response_energy.count(patch) == 1);
    const double residual = std::stod(row[loops.Column("loop residual")]);
    const double scale = std::stod(row[loops.Column("loop scale (V)")]);
    CHECK(scale > 0.0);
    if (kind == 0)
    {
      max_loop = std::max(max_loop, residual);
      paths_by_patch[patch].emplace_back(residual, scale * scale);
    }
  }
  CHECK(max_loop > 0.0);
  CHECK_THAT(max_loop, WithinRel(gate_loop_residual, 1.0e-9));
  CHECK_THAT(max_patch_loop, WithinRel(gate_loop_residual, 1.0e-9));
  double weighted_squared = 0.0, weight_sum = 0.0;
  for (const auto &[patch, paths] : paths_by_patch)
  {
    double scale_squared = 0.0;
    for (const auto &[residual, s2] : paths)
    {
      (void)residual;
      scale_squared += s2;
    }
    for (const auto &[residual, s2] : paths)
    {
      const double weight = patch_response_energy.at(patch) * s2 / scale_squared;
      weighted_squared += weight * residual * residual;
      weight_sum += weight;
    }
  }
  REQUIRE(weight_sum > 0.0);
  {
    double a2 = 0.0, an = 0.0, b2 = 0.0, bn = 0.0;
    for (const auto &[patch, paths] : paths_by_patch)
    {
      for (const auto &[residual, s2] : paths)
      {
        a2 += residual * residual;
        an += 1.0;
        b2 += s2 * residual * residual;
        bn += s2;
      }
    }
    UNSCOPED_INFO("unweighted RMS " << std::sqrt(a2 / an) << " scale2-weighted RMS "
                                    << std::sqrt(b2 / bn) << " patches "
                                    << paths_by_patch.size() << " energies "
                                    << patch_response_energy.size());
  }
  CHECK_THAT(std::sqrt(weighted_squared / weight_sum), WithinRel(gate_weighted, 1.0e-7));
#endif
}

// The default (no "PatchEnergy"): the same fixture writes none of the per-patch files, the
// operator collects no claims, and surface-response-model-energy.csv is byte-identical to
// the PatchEnergy run's (the export reads values the totals compute anyway).
TEST_CASE_METHOD(test::SurfaceResponseFiles,
                 "Maxwell per-patch export is off by default (decision 922 B)",
                 "[postoperator][surfaceresponseoperator][cornerarmtrim][3d][maxwell]"
                 "[patchenergy][Serial][Parallel]")
{
#if !defined(MFEM_USE_GSLIB)
  SKIP("SurfaceResponseOperator requires MFEM_USE_GSLIB");
#else
  const fs::path mesh_path = temp.temp_dir / "patch-off-lead.mesh";
  const auto basis_path = temp.temp_dir / "patch-off-basis-points.csv";
  const auto edge_library_path = temp.temp_dir / "fabrication-process-patch-off-edges.json";
  MPI_Comm comm = Mpi::World();
  if (Mpi::Root(comm))
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
    library["Name"] = "unit-test-process-patch-off-edges";
    json models = json::array();
    for (auto &model : library["Models"])
    {
      if (model["Topology"] == "IsolatedEdge")
      {
        model["Interfaces"] = {{{"Type", "SA"}, {"Coupon", 1}}};
        model["BasisPoints"] = basis_path.string();
        models.push_back(model);
      }
    }
    library["Models"] = std::move(models);
    std::ofstream output(edge_library_path);
    output << library.dump(2) << "\n";
  }
  Mpi::Barrier(comm);
  json config = IslandConfig();
  config["Model"]["Mesh"] = mesh_path.string();
  config["Model"]["L0"] = 1.0;
  config["Problem"]["Type"] = "Eigenmode";
  config["Boundaries"]["Ground"]["Attributes"] = {4, 9};
  config["Boundaries"].erase("Terminal");
  auto Run = [&](bool patch_energy, const char *stem)
  {
    json run = config;
    const fs::path output_dir = temp.temp_dir / stem;
    run["Problem"]["Output"] = output_dir.string();
    run["Solver"] = {{"Order", 1},
                     {"Eigenmode", {{"Target", 1.0}}},
                     {"SurfaceResponseCorrection",
                      {{"Library", edge_library_path.string()},
                       {"TargetInterfaces", {4}},
                       {"UnmatchedPolicy", "Warn"},
                       {"CorrectionMode", "PostprocessOnly"}}}};
    if (patch_energy)
    {
      run["Solver"]["SurfaceResponseCorrection"]["PatchEnergy"] = true;
    }
    IoData iodata(run, false);
    iodata.boundaries.cracked_attributes.insert(9);
    MakeOutputFolder(iodata, comm);
    std::vector<std::unique_ptr<Mesh>> meshes;
    meshes.push_back(std::make_unique<Mesh>(MakeLeadMesh(0.0, 0.75)));
    SpaceOperator space_op(iodata, meshes);
    PostOperator<ProblemType::EIGENMODE> post_op(iodata, space_op);
    auto *response = post_op.GetSurfaceResponseOperator();
    REQUIRE(response);
    CHECK(response->GetModelRawClaimPortions().portions.empty() == !patch_energy);
    mfem::FunctionCoefficient potential_coefficient(
        [](const mfem::Vector &x) { return std::sinh(0.7 * (x[1] - 0.5)); });
    mfem::ParGridFunction potential(&space_op.GetH1Space().Get());
    potential.ProjectCoefficient(potential_coefficient);
    mfem::ParDiscreteLinearOperator gradient(&space_op.GetH1Space().Get(),
                                             &space_op.GetNDSpace().Get());
    gradient.AddDomainInterpolator(new mfem::GradientInterpolator());
    gradient.Assemble();
    GridFunction field(space_op.GetNDSpace(), true);
    gradient.Mult(potential, field.Real());
    field.Real() *= -1.0;
    ComplexVector e(space_op.GetNDSpace().GetTrueVSize());
    field.Real().GetTrueDofs(e.Real());
    e.Imag() = 0.0;
    ComplexVector b(space_op.GetRTSpace().GetTrueVSize());
    b = 0.0;
    post_op.MeasureAndPrintAll(0, e, b, std::complex<double>(1.0, 0.0), 0.0, 0.0, 1);
    Mpi::Barrier(comm);
    return output_dir;
  };
  const fs::path on_dir = Run(true, "patch-on");
  const fs::path off_dir = Run(false, "patch-off");
  if (!Mpi::Root(comm))
  {
    return;
  }
  for (const char *name :
       {"surface-response-patches.csv", "surface-response-patch-energy.csv",
        "surface-response-patch-loops.csv", "surface-response-model-raw-energy.csv",
        "surface-response-model-raw-groups.csv"})
  {
    CHECK(fs::is_regular_file(on_dir / name));
    CHECK(!fs::exists(off_dir / name));
  }
  for (const char *name : {"surface-response-model-energy.csv", "surface-Q-corrected.csv",
                           "surface-response-confidence.csv"})
  {
    INFO(name);
    CHECK(ReadFile(on_dir / name) == ReadFile(off_dir / name));
  }
#endif
}

// OI-MX-LOOPGATE (decision 955) fixture: the 922-B lead (MakeLeadMesh(0.0, 0.75)) with the
// one isolated-edge model of the convex library on a 4-knot ring (the square of side 2 R =
// 0.4 shifted R / 2 along v: the path's start knot, the nearest to the PEC anchor (-R, 0),
// is then unique, while the rotation about the square's centre stays uniform along every
// side), under an Eigenmode PostOperator with "PatchEnergy" and the one target interface 4.
// Returns the Maxwell configuration; writes the mesh / basis / library files on root. The
// field E = -a + c z x (x - o) (a in the lead plane, o between the two long edges) is exact
// in the lowest-order ND hex space.
struct LoopFitLeadFixture
{
  static constexpr double R = 0.2;
  static constexpr std::array<std::array<double, 2>, 4> knots = {
      {{-R, -0.5 * R}, {R, -0.5 * R}, {R, 1.5 * R}, {-R, 1.5 * R}}};
  static constexpr std::array<double, 3> a = {0.3, 0.7, 0.0};
  static constexpr std::array<double, 3> o = {0.5, 0.5, 0.0};
  static std::array<double, 3> Field(const std::array<double, 3> &x, double c)
  {
    return {-a[0] - c * (x[1] - o[1]), -a[1] + c * (x[0] - o[0]), -a[2]};
  }
  static double Dot(const std::array<double, 3> &x, const std::array<double, 3> &y)
  {
    return x[0] * y[0] + x[1] * y[1] + x[2] * y[2];
  }
  static json Config(const test::SurfaceResponseFiles &files, const fs::path &output_dir,
                     const std::string &tag)
  {
    const fs::path mesh_path = files.temp.temp_dir / (tag + "-lead.mesh");
    const auto basis_path = files.temp.temp_dir / (tag + "-basis-points.csv");
    const auto edge_library_path =
        files.temp.temp_dir / ("fabrication-process-" + tag + "-edges.json");
    if (Mpi::Root(Mpi::World()))
    {
      {
        mfem::Mesh serial = MakeSerialLeadMesh(0.0, 0.75);
        std::ofstream output(mesh_path);
        serial.Print(output);
      }
      {
        std::ofstream output(basis_path);
        output << "x,y,z\n";
        for (const auto &knot : knots)
        {
          output << knot[0] << "," << knot[1] << ",0.0\n";
        }
      }
      std::ifstream input(files.convex_library_3d_path);
      REQUIRE(input);
      json library = json::parse(input);
      library["Name"] = "unit-test-process-" + tag + "-edges";
      json models = json::array();
      for (auto &model : library["Models"])
      {
        if (model["Topology"] == "IsolatedEdge")
        {
          model["Interfaces"] = {{{"Type", "SA"}, {"Coupon", 1}}};
          model["BasisPoints"] = basis_path.string();
          models.push_back(model);
        }
      }
      REQUIRE(models.size() == 1);
      library["Models"] = std::move(models);
      std::ofstream output(edge_library_path);
      output << library.dump(2) << "\n";
    }
    Mpi::Barrier(Mpi::World());
    json config = files.IslandConfig();
    config["Model"]["Mesh"] = mesh_path.string();
    config["Model"]["L0"] = 1.0;
    config["Boundaries"]["Ground"]["Attributes"] = {4, 9};
    config["Boundaries"].erase("Terminal");
    config["Problem"]["Type"] = "Eigenmode";
    config["Problem"]["Output"] = output_dir.string();
    config["Solver"] = {{"Order", 1},
                        {"Eigenmode", {{"Target", 1.0}}},
                        {"SurfaceResponseCorrection",
                         {{"Library", edge_library_path.string()},
                          {"TargetInterfaces", {4}},
                          {"UnmatchedPolicy", "Warn"},
                          {"CorrectionMode", "PostprocessOnly"},
                          {"PatchEnergy", true}}}};
    return config;
  }
  // The field as a complex ND grid function (imaginary part 0) and its real true dofs.
  static std::unique_ptr<GridFunction> MakeField(SpaceOperator &space_op, double c,
                                                 Vector &e_real)
  {
    mfem::VectorFunctionCoefficient coefficient(3,
                                                [c](const mfem::Vector &x, mfem::Vector &v)
                                                {
                                                  const auto value =
                                                      Field({x[0], x[1], x[2]}, c);
                                                  v.SetSize(3);
                                                  v(0) = value[0];
                                                  v(1) = value[1];
                                                  v(2) = value[2];
                                                });
    auto field = std::make_unique<GridFunction>(space_op.GetNDSpace(), true);
    field->Real().ProjectCoefficient(coefficient);
    field->Imag() = 0.0;
    e_real.SetSize(space_op.GetNDSpace().GetTrueVSize());
    field->Real().GetTrueDofs(e_real);
    field->Real().SetFromTrueDofs(e_real);
    field->Real().ExchangeFaceNbrData();
    return field;
  }
  // The fitted (gradient-part) knot potentials of a patch: -a + c z x (q - o) (q = the
  // square's centre) integrated from the PEC anchor A = origin - R u, minus the anchor
  // line's own circulating integral (a tree step in no cycle: the chain's value stands).
  // Also the knot positions and the anchor for the caller's line integrals.
  struct Geometry
  {
    std::array<std::array<double, 3>, 4> device{};
    std::array<double, 3> anchor{};
    std::array<double, 3> centre{};
    std::size_t start = 0;
    bool long_edge = false;
    std::array<double, 4> fit{};
  };
  static Geometry Expect(const SurfaceResponseOperator::PatchAssignment &assignment,
                         double c)
  {
    Geometry g;
    g.long_edge = std::abs(assignment.axis_w[2]) > 0.5;
    for (std::size_t k = 0; k < 4; k++)
    {
      for (int d = 0; d < 3; d++)
      {
        g.device[k][d] = assignment.origin[d] + knots[k][0] * assignment.axis_u[d] +
                         knots[k][1] * assignment.axis_v[d];
        g.centre[d] += 0.25 * g.device[k][d];
      }
    }
    for (int d = 0; d < 3; d++)
    {
      g.anchor[d] = assignment.origin[d] - R * assignment.axis_u[d];
    }
    double start_distance = mfem::infinity();
    for (std::size_t k = 0; k < 4; k++)
    {
      std::array<double, 3> delta;
      for (int d = 0; d < 3; d++)
      {
        delta[d] = g.device[k][d] - g.anchor[d];
      }
      const double distance = std::sqrt(Dot(delta, delta));
      if (distance < start_distance - 1.0e-6)
      {
        g.start = k;
        start_distance = distance;
      }
    }
    REQUIRE(g.start == 0);
    const std::array<double, 3> gradient = {a[0] + c * (g.centre[1] - o[1]),
                                            a[1] - c * (g.centre[0] - o[0]), a[2]};
    std::array<double, 3> anchor_midpoint, anchor_delta;
    for (int d = 0; d < 3; d++)
    {
      anchor_midpoint[d] = 0.5 * (g.anchor[d] + g.device[g.start][d]);
      anchor_delta[d] = g.device[g.start][d] - g.anchor[d];
    }
    const std::array<double, 3> rotation_at_midpoint = {
        -c * (anchor_midpoint[1] - g.centre[1]), c * (anchor_midpoint[0] - g.centre[0]),
        0.0};
    const double anchor_circulation = Dot(rotation_at_midpoint, anchor_delta);
    for (std::size_t k = 0; k < 4; k++)
    {
      std::array<double, 3> delta;
      for (int d = 0; d < 3; d++)
      {
        delta[d] = g.device[k][d] - g.anchor[d];
      }
      g.fit[k] = Dot(gradient, delta) - anchor_circulation;
    }
    return g;
  }
};

// OI-MX-LOOPGATE (decision 955), the FAILING-FIRST case (main's interfaces only): on the
// fixture above with the circulating field (c = 0.5) the Maxwell trace of GetPatchTraces
// equals the gradient part's knot potentials on every long-edge patch (the chain lift of
// main carries the ring's loop EMF into the knots after the start knot in route order and
// FAILS), and the print writes surface-response-loop-fit.csv and the "loop share" column of
// surface-response-patch-energy.csv (absent on main).
TEST_CASE_METHOD(test::SurfaceResponseFiles,
                 "Maxwell loop fit recovers the gradient potentials of a circulating field "
                 "(decision 955, failing-first)",
                 "[postoperator][surfaceresponseoperator][cornerarmtrim][3d][maxwell]"
                 "[loopfit][Serial][Parallel]")
{
#if !defined(MFEM_USE_GSLIB)
  SKIP("SurfaceResponseOperator requires MFEM_USE_GSLIB");
#else
  const fs::path output_dir = temp.temp_dir / "loop-fit-failing-first";
  const json config = LoopFitLeadFixture::Config(*this, output_dir, "loop-fit-ff");
  IoData iodata(config, false);
  iodata.boundaries.cracked_attributes.insert(9);
  MPI_Comm comm = Mpi::World();
  MakeOutputFolder(iodata, comm);
  std::vector<std::unique_ptr<Mesh>> meshes;
  meshes.push_back(std::make_unique<Mesh>(MakeLeadMesh(0.0, 0.75)));
  SpaceOperator space_op(iodata, meshes);
  PostOperator<ProblemType::EIGENMODE> post_op(iodata, space_op);
  auto *response = post_op.GetSurfaceResponseOperator();
  REQUIRE(response);
  Vector e_real;
  const auto field = LoopFitLeadFixture::MakeField(space_op, 0.5, e_real);
  const auto traces = response->GetPatchTraces(e_real);
  REQUIRE(traces.size() == static_cast<std::size_t>(response->GetPatchCount()));
  if (Mpi::Root(comm))
  {
    const auto &assignments = response->GetPatchAssignments();
    REQUIRE(assignments.size() == traces.size());
    int long_edges = 0;
    for (std::size_t p = 0; p < assignments.size(); p++)
    {
      const auto geometry = LoopFitLeadFixture::Expect(assignments[p], 0.5);
      if (!geometry.long_edge)
      {
        continue;
      }
      long_edges++;
      INFO("patch " << assignments[p].global_index);
      REQUIRE(traces[p].patch == assignments[p].global_index);
      REQUIRE(traces[p].coefficients.size() == 4);
      for (std::size_t k = 0; k < 4; k++)
      {
        CHECK_THAT(traces[p].coefficients[k], WithinAbs(geometry.fit[k], 1.0e-10));
      }
    }
    REQUIRE(long_edges > 0);
  }
  ComplexVector e(space_op.GetNDSpace().GetTrueVSize()),
      b(space_op.GetRTSpace().GetTrueVSize());
  e.Real() = e_real;
  e.Imag() = 0.0;
  b = 0.0;
  post_op.MeasureAndPrintAll(0, e, b, std::complex<double>(1.0, 0.0), 0.0, 0.0, 1);
  Mpi::Barrier(comm);
  if (!Mpi::Root(comm))
  {
    return;
  }
  CHECK(fs::is_regular_file(output_dir / "surface-response-loop-fit.csv"));
  CHECK(fs::is_regular_file(output_dir / "surface-response-loop-remainder-energy.csv"));
  const TextCsv patch_energy =
      ReadTextCsv(output_dir / "surface-response-patch-energy.csv");
  bool share_column = false;
  for (const auto &name : patch_energy.header)
  {
    share_column = share_column || name.find("loop share") != std::string::npos;
  }
  CHECK(share_column);
#endif
}

// OI-MX-LOOPGATE (decisions 950 / 952 / 954 / 955): the weighted least-squares gradient /
// circulating separation of the Maxwell line-integral trace lift, on the 922-B lead fixture
// (a 4-knot square ring of side 2 R = 0.4 around the PEC lead's edges, shifted R / 2 off
// the edge along v, one isolated-edge model, the PEC anchor R inside the edge). The field E
// = -a + c z x (x - o) is exact in the lowest-order ND hex space, so every line integral is
// the midpoint rule on its straight legs and every expectation below is analytic per patch
// from the patch table: on a long edge (axis +-z) the rotation about the square's centre q
// is the uniform tangential field c R along the ring (plus a constant field c z x (q - o),
// a gradient), so the ring's loop misclosure is c R x its perimeter while the fit recovers
// the gradient potentials exactly (the anchor line, in no cycle, keeps its own circulating
// integral: a gauge term named in the lane record); on the cap (axis +-x) the rotation is
// normal to the ring plane and the loops close. (1) c = 0 (curl-free): the fit reproduces
// the chain = the analytic potentials, every share is 0 to roundoff, no remainder, the
// loop-fit file is written and the remainder file carries the Total rows with 0 (passes on
// main except for the files); (2) c = 0.5 (a remainder on every long-edge patch, no
// fallback): the fitted trace = the gradient potentials; per patch the circulating norm r^T
// G^-1 r = (sum of the ring integrals)^2 / perimeter and the route norm = sum L_j^2 / l_j;
// the chain-based loop residual = the analytic one; the remainder energy = sum_p s_p x the
// raw within-R energy on the patch's own claims (an independent partitioned integral with
// per-patch labels); the corrected E_surf = E_out + fabricated + uncovered + remainder on
// both closures; (3) c = 10 (a fallback on every long-edge patch): those patches carry 0
// energies and the fallback flag, rho = 1 on their claims, the same identities. (4) The
// fit's transpose: <Mult p, q> = <p, Mult q> on random vectors. On main (2) and (3) FAIL:
// the chain's knot values carry the loop EMF in route order and none of the records exists.
TEST_CASE_METHOD(test::SurfaceResponseFiles,
                 "Maxwell loop fit separates the circulating part of the trace lift "
                 "(decision 955)",
                 "[postoperator][surfaceresponseoperator][cornerarmtrim][3d][maxwell]"
                 "[loopfit][Serial][Parallel]")
{
#if !defined(MFEM_USE_GSLIB)
  SKIP("SurfaceResponseOperator requires MFEM_USE_GSLIB");
#else
  const fs::path output_dir = temp.temp_dir / "loop-fit-maxwell";
  const json config = LoopFitLeadFixture::Config(*this, output_dir, "loop-fit");
  IoData iodata(config, false);
  iodata.boundaries.cracked_attributes.insert(9);
  MPI_Comm comm = Mpi::World();
  MakeOutputFolder(iodata, comm);
  std::vector<std::unique_ptr<Mesh>> meshes;
  meshes.push_back(std::make_unique<Mesh>(MakeLeadMesh(0.0, 0.75)));
  SpaceOperator space_op(iodata, meshes);
  PostOperator<ProblemType::EIGENMODE> post_op(iodata, space_op);
  auto *response = post_op.GetSurfaceResponseOperator();
  REQUIRE(response);
  const int patch_count = response->GetPatchCount();
  REQUIRE(patch_count > 0);
  // Diagnostics.LoopFit (A1.4; decision 962): the remainder is measured on the 3D patches'
  // raw claims (the two-dimensional path reports false and measures none).
  {
    const auto loop_fit_diagnostics = response->GetStatistics()["Diagnostics"]["LoopFit"];
    CHECK(loop_fit_diagnostics["RemainderMeasured"] == true);
    CHECK(loop_fit_diagnostics["Cycles"].get<long long int>() == patch_count);
    CHECK(loop_fit_diagnostics["UnconstrainedKnots"].get<long long int>() == 0);
  }

  // Per patch (root: the patch table), the analytic line integrals of the five legs the
  // path evaluates (the anchor line A -> the start knot and the four ring sides in cyclic
  // order; a linear field is exactly integrated by the midpoint rule on a straight leg).
  struct PatchExpectation
  {
    bool long_edge = false;
    std::array<double, 4> fit{};  // the fitted knot potentials
    double circulation = 0.0;     // (sum of the ring integrals)^2 / perimeter
    double route_norm = 0.0;      // sum_j L_j^2 / l_j
    double share = 0.0;
    double residual = 0.0;  // |sum of the ring integrals| / sum |ring integrals|
  };
  auto Expect = [&](const SurfaceResponseOperator::PatchAssignment &assignment, double c)
  {
    const auto geometry = LoopFitLeadFixture::Expect(assignment, c);
    PatchExpectation expectation;
    expectation.long_edge = geometry.long_edge;
    expectation.fit = geometry.fit;
    auto Leg = [&](const std::array<double, 3> &p, const std::array<double, 3> &q)
    {
      std::array<double, 3> midpoint, delta;
      for (int d = 0; d < 3; d++)
      {
        midpoint[d] = 0.5 * (p[d] + q[d]);
        delta[d] = q[d] - p[d];
      }
      return std::make_pair(
          LoopFitLeadFixture::Dot(LoopFitLeadFixture::Field(midpoint, c), delta),
          std::sqrt(LoopFitLeadFixture::Dot(delta, delta)));
    };
    const auto [anchor_integral, anchor_length] =
        Leg(geometry.anchor, geometry.device[geometry.start]);
    double ring_sum = 0.0, ring_abs = 0.0, perimeter = 0.0;
    expectation.route_norm = anchor_integral * anchor_integral / anchor_length;
    for (std::size_t k = 0; k < 4; k++)
    {
      const auto [integral, length] = Leg(geometry.device[k], geometry.device[(k + 1) % 4]);
      ring_sum += integral;
      ring_abs += std::abs(integral);
      perimeter += length;
      expectation.route_norm += integral * integral / length;
    }
    expectation.circulation = ring_sum * ring_sum / perimeter;
    expectation.share = expectation.circulation / expectation.route_norm;
    expectation.residual = ring_abs > 0.0 ? std::abs(ring_sum) / ring_abs : 0.0;
    return expectation;
  };

  ComplexVector b(space_op.GetRTSpace().GetTrueVSize());
  b = 0.0;
  const std::complex<double> omega(1.0, 0.0);
  struct Case
  {
    double c;
    int mode;
    bool remainder;
    bool fallback;
  };
  const std::array<Case, 3> cases = {
      {{0.0, 1, false, false}, {0.5, 2, true, false}, {10.0, 3, true, true}}};
  std::map<int, double> expected_remainder;  // mode -> sum_p rho_p E_p (interface 4)
  for (const auto &test_case : cases)
  {
    INFO("c = " << test_case.c);
    Vector e_real;
    const auto field = LoopFitLeadFixture::MakeField(space_op, test_case.c, e_real);
    const auto measured = response->GetMaxwellResponse(*field, omega, true);
    const auto traces = response->GetPatchTraces(e_real);
    REQUIRE(traces.size() == static_cast<std::size_t>(patch_count));
    REQUIRE(measured.patch_loop_shares.size() == static_cast<std::size_t>(patch_count));
    std::map<int, double> rho;  // the remainder patches' rho_p
    for (const auto &remainder : measured.loop_remainders)
    {
      rho[remainder.patch] = remainder.fallback ? 1.0 : remainder.share;
    }
    if (Mpi::Root(comm))
    {
      const auto &assignments = response->GetPatchAssignments();
      REQUIRE(assignments.size() == static_cast<std::size_t>(patch_count));
      std::map<int, double> loop_residuals;
      for (const auto &loop : measured.path_loop_residuals)
      {
        if (loop.kind == 0)
        {
          loop_residuals[loop.patch] = std::max(loop_residuals[loop.patch], loop.residual);
        }
      }
      int long_edges = 0;
      for (std::size_t p = 0; p < assignments.size(); p++)
      {
        const auto &assignment = assignments[p];
        const auto &trace = traces[p];
        const auto &share = measured.patch_loop_shares[p];
        REQUIRE(trace.patch == assignment.global_index);
        REQUIRE(share.patch == assignment.global_index);
        REQUIRE(trace.coefficients.size() == 4);
        const auto expectation = Expect(assignment, test_case.c);
        INFO("patch " << assignment.global_index
                      << (expectation.long_edge ? " (long edge)" : " (cap)"));
        long_edges += expectation.long_edge ? 1 : 0;
        // (1) / (2) / (3) (i): the fitted trace = the gradient potentials.
        for (std::size_t k = 0; k < 4; k++)
        {
          CHECK_THAT(trace.coefficients[k], WithinAbs(expectation.fit[k], 1.0e-10));
        }
        // (ii) the share, and the remainder record's norms.
        if (test_case.c == 0.0 || !expectation.long_edge)
        {
          CHECK(expectation.share <= 1.0e-24);
          CHECK(share.share <= 1.0e-14);
          CHECK(share.fallback == 0);
          CHECK(rho.count(assignment.global_index) == 0);
        }
        else
        {
          // The constants were chosen so that the long-edge patches are remainder patches
          // (c = 0.5) or fallback patches (c = 10): asserted on the analytic share.
          REQUIRE(expectation.share > SurfaceResponseOperator::kLoopFitShareFloor);
          REQUIRE((expectation.share > SurfaceResponseOperator::kLoopFitFallbackShare) ==
                  test_case.fallback);
          CHECK_THAT(share.share, WithinRel(expectation.share, 1.0e-9));
          CHECK(share.fallback == (test_case.fallback ? 1 : 0));
          const auto remainder = std::find_if(
              measured.loop_remainders.begin(), measured.loop_remainders.end(),
              [&](const auto &entry) { return entry.patch == assignment.global_index; });
          REQUIRE(remainder != measured.loop_remainders.end());
          CHECK_THAT(remainder->circulation, WithinRel(expectation.circulation, 1.0e-9));
          CHECK_THAT(remainder->total, WithinRel(expectation.route_norm, 1.0e-9));
          CHECK(remainder->fallback == (test_case.fallback ? 1 : 0));
          // (iii) the chain-based loop residual export = the analytic ring residual.
          REQUIRE(loop_residuals.count(assignment.global_index) == 1);
          CHECK(expectation.residual > 1.0e-3);
          CHECK_THAT(loop_residuals.at(assignment.global_index),
                     WithinRel(expectation.residual, 1.0e-9));
        }
      }
      REQUIRE(long_edges > 0);
      CHECK(measured.loop_remainder_patch_count ==
            (test_case.remainder && !test_case.fallback ? long_edges : 0));
      CHECK(measured.loop_fallback_patch_count == (test_case.fallback ? long_edges : 0));
      CHECK(static_cast<int>(measured.loop_remainders.size()) ==
            (test_case.remainder ? long_edges : 0));
    }
    // (iv) the remainder energy = sum_p rho_p x the raw within-R energy on the patch's own
    // claims: the claims per patch through the public portion accessor, labelled per patch,
    // partitioned by the (unweighted) exact partition.
    if (!measured.loop_remainders.empty())
    {
      std::vector<SurfacePostOperator::UncoveredPerimeterPortion> perimeter;
      for (const auto &remainder : measured.loop_remainders)
      {
        const auto portions = response->GetLoopRemainderPortions({remainder});
        REQUIRE(!portions.portions.empty());
        for (const auto &portion : portions.portions)
        {
          perimeter.push_back({portion.p0, portion.p1, fmt::format("p{}", remainder.patch),
                               portion.feature});
        }
      }
      const auto partition = post_op.GetSurfacePostOp().GetInterfacePartitionedEdgeEnergies(
          {4}, *field, nullptr, perimeter);
      double expected = 0.0;
      int positive = 0;
      for (const auto &remainder : measured.loop_remainders)
      {
        const auto label = fmt::format("p{}", remainder.patch);
        positive += partition.at(4).at(label) > 0.0 ? 1 : 0;
        REQUIRE(partition.at(4).count(label) == 1);
        CHECK(partition.at(4).at(label) >= 0.0);
        expected += rho.at(remainder.patch) * partition.at(4).at(label);
      }
      REQUIRE(positive > 0);
      REQUIRE(expected > 0.0);
      expected_remainder[test_case.mode] = expected;
    }
    else
    {
      expected_remainder[test_case.mode] = 0.0;
    }
    ComplexVector e(space_op.GetNDSpace().GetTrueVSize());
    e.Real() = e_real;
    e.Imag() = 0.0;
    post_op.MeasureAndPrintAll(test_case.mode - 1, e, b, omega, 0.0, 0.0, test_case.mode);
  }

  // (4) The transpose of the fit through the operator's Mult.
  {
    Vector p(space_op.GetNDSpace().GetTrueVSize()), q(p.Size()), Mp, Mq;
    auto *p_data = p.HostWrite();
    auto *q_data = q.HostWrite();
    for (int i = 0; i < p.Size(); i++)
    {
      p_data[i] = std::sin(0.37 * (i + 1 + 11 * Mpi::Rank(comm)));
      q_data[i] = std::cos(0.53 * (i + 3 + 7 * Mpi::Rank(comm)));
    }
    p.SetSubVector(space_op.GetNDDbcTDofLists().back(), 0.0);
    q.SetSubVector(space_op.GetNDDbcTDofLists().back(), 0.0);
    response->Mult(p, Mp);
    response->Mult(q, Mq);
    const double pMq = linalg::Dot(comm, p, Mq), qMp = linalg::Dot(comm, q, Mp);
    REQUIRE(std::abs(pMq) > 0.0);
    CHECK_THAT(pMq, WithinRel(qMp, 1.0e-10));
  }

  Mpi::Barrier(comm);
  if (!Mpi::Root(comm))
  {
    return;
  }
  auto Joule = [&](double energy)
  { return iodata.units.Dimensionalize<Units::ValueType::ENERGY>(energy); };
  auto Trim = [](std::string text)
  {
    const auto begin = text.find_first_not_of(" \t\r");
    const auto end = text.find_last_not_of(" \t\r");
    return begin == std::string::npos ? std::string() : text.substr(begin, end - begin + 1);
  };
  auto Int = [](const std::string &field)
  { return static_cast<int>(std::lround(std::stod(field))); };
  auto ReadTrimmedCsv = [&](const fs::path &path)
  {
    TextCsv csv = ReadTextCsv(path);
    for (auto &name : csv.header)
    {
      name = Trim(name);
    }
    for (auto &row : csv.rows)
    {
      for (auto &cell : row)
      {
        cell = Trim(cell);
      }
    }
    return csv;
  };
  // The loop-fit file: one row per (mode, evaluation 0 / 1) with the counts.
  REQUIRE(fs::is_regular_file(output_dir / "surface-response-loop-fit.csv"));
  const TextCsv fit = ReadTrimmedCsv(output_dir / "surface-response-loop-fit.csv");
  CHECK(fit.header ==
        std::vector<std::string>{"source", "evaluation", "m", "exc", "max loop share",
                                 "response-weighted loop share", "remainder patches",
                                 "fallback patches", "weighted claim length (m)",
                                 "withheld domain correction (J)"});
  REQUIRE(fit.rows.size() == 6);
  std::map<int, std::array<int, 2>> fit_counts;  // mode -> (remainder, fallback)
  std::map<int, double> fit_max_share;
  for (const auto &row : fit.rows)
  {
    const int mode = Int(row[fit.Column("m")]);
    CHECK(Int(row[fit.Column("source")]) == mode);
    CHECK((Int(row[fit.Column("evaluation")]) == 0 ||
           Int(row[fit.Column("evaluation")]) == 1));
    fit_counts[mode] = {Int(row[fit.Column("remainder patches")]),
                        Int(row[fit.Column("fallback patches")])};
    fit_max_share[mode] = std::stod(row[fit.Column("max loop share")]);
  }
  CHECK(fit_max_share.at(1) <= 1.0e-14);
  CHECK(fit_counts.at(1) == std::array<int, 2>{0, 0});
  CHECK(fit_counts.at(2)[0] > 0);
  CHECK(fit_counts.at(2)[1] == 0);
  CHECK(fit_counts.at(3)[0] == 0);
  CHECK(fit_counts.at(3)[1] == fit_counts.at(2)[0]);
  CHECK(fit_max_share.at(3) > SurfaceResponseOperator::kLoopFitFallbackShare);
  // The remainder file: per mode and evaluation a Total row (0 on mode 1) and the model row
  // of the one isolated-edge model where a remainder exists.
  REQUIRE(fs::is_regular_file(output_dir / "surface-response-loop-remainder-energy.csv"));
  const TextCsv remainder =
      ReadTrimmedCsv(output_dir / "surface-response-loop-remainder-energy.csv");
  CHECK(remainder.header == std::vector<std::string>{"source", "evaluation", "model",
                                                     "patches", "weighted claim length (m)",
                                                     "remainder energy[4] (J)", "m",
                                                     "exc"});
  std::map<std::pair<int, int>, double> remainder_total;  // (mode, ev) -> Total energy
  std::map<int, int> remainder_rows;                      // mode -> rows per evaluation 0
  for (const auto &row : remainder.rows)
  {
    const int mode = Int(row[remainder.Column("m")]);
    const int ev = Int(row[remainder.Column("evaluation")]);
    if (row[remainder.Column("model")] == "Total")
    {
      remainder_total[{mode, ev}] =
          std::stod(row[remainder.Column("remainder energy[4] (J)")]);
      CHECK(Int(row[remainder.Column("patches")]) ==
            fit_counts.at(mode)[0] + fit_counts.at(mode)[1]);
    }
    else if (ev == 0)
    {
      remainder_rows[mode]++;
      CHECK(row[remainder.Column("model")].rfind("m", 0) == 0);
    }
  }
  CHECK(remainder_rows.count(1) == 0);
  CHECK(remainder_rows.at(2) == 1);
  CHECK(remainder_rows.at(3) == 1);
  for (const auto &test_case : cases)
  {
    INFO("mode " << test_case.mode);
    for (const int ev : {0, 1})
    {
      REQUIRE(remainder_total.count({test_case.mode, ev}) == 1);
      if (expected_remainder.at(test_case.mode) == 0.0)
      {
        CHECK(remainder_total.at({test_case.mode, ev}) == 0.0);
      }
      else
      {
        CHECK_THAT(remainder_total.at({test_case.mode, ev}),
                   WithinRel(Joule(expected_remainder.at(test_case.mode)), 1.0e-9));
      }
    }
  }
  // The corrected E_surf = E_out + fabricated + uncovered + remainder on both closures.
  const TextCsv edge = ReadTrimmedCsv(output_dir / "surface-Q-edge.csv");
  const TextCsv corrected = ReadTrimmedCsv(output_dir / "surface-Q-corrected.csv");
  const TextCsv model_energy =
      ReadTrimmedCsv(output_dir / "surface-response-model-energy.csv");
  const TextCsv uncovered =
      ReadTrimmedCsv(output_dir / "surface-response-uncovered-energy.csv");
  const TextCsv patch_energy =
      ReadTrimmedCsv(output_dir / "surface-response-patch-energy.csv");
  CHECK(patch_energy.header ==
        std::vector<std::string>{
            "source", "evaluation", "patch", "model", "feature", "weight", "cell begin (m)",
            "cell end (m)", "origin x (m)", "origin y (m)", "origin z (m)",
            "domain correction (J)", "fabricated surface energy[4] (J)", "m", "exc",
            "max loop residual", "loop share", "loop fallback", "loop weight (J)"});
  for (const auto &test_case : cases)
  {
    INFO("mode " << test_case.mode);
    double e_out = std::numeric_limits<double>::quiet_NaN();
    for (const auto &row : edge.rows)
    {
      if (Int(row[edge.Column("m")]) == test_case.mode &&
          Int(row[edge.Column("interface")]) == 4)
      {
        e_out = std::stod(row[edge.Column("E_out (J)")]);
      }
    }
    REQUIRE(std::isfinite(e_out));
    std::array<double, 2> fabricated{}, kept{}, surf{};
    for (const auto &row : model_energy.rows)
    {
      if (Int(row[model_energy.Column("m")]) == test_case.mode)
      {
        fabricated[Int(row[model_energy.Column("evaluation")])] +=
            std::stod(row[model_energy.Column("fabricated surface energy[4] (J)")]);
      }
    }
    for (const auto &row : uncovered.rows)
    {
      if (Int(row[uncovered.Column("m")]) == test_case.mode &&
          row[uncovered.Column("type")] == "Total" &&
          Int(row[uncovered.Column("evaluation")]) < 2)
      {
        kept[Int(row[uncovered.Column("evaluation")])] +=
            std::stod(row[uncovered.Column("uncovered raw energy[4] (J)")]);
      }
    }
    for (const auto &row : corrected.rows)
    {
      if (Int(row[corrected.Column("m")]) == test_case.mode &&
          Int(row[corrected.Column("interface")]) == 4)
      {
        surf[0] = std::stod(row[corrected.Column("E_surf corrected (J)")]);
        surf[1] = std::stod(row[corrected.Column("E_surf corrected fixed-flux (J)")]);
      }
    }
    for (const int ev : {0, 1})
    {
      CHECK_THAT(surf[ev], WithinRel(e_out + fabricated[ev] + kept[ev] +
                                         remainder_total.at({test_case.mode, ev}),
                                     1.0e-10));
    }
    // The per-patch columns: a fallback patch carries 0 energies and the flag.
    int fallback_rows = 0;
    for (const auto &row : patch_energy.rows)
    {
      if (Int(row[patch_energy.Column("m")]) != test_case.mode)
      {
        continue;
      }
      const int flag = Int(row[patch_energy.Column("loop fallback")]);
      const double share = std::stod(row[patch_energy.Column("loop share")]);
      CHECK(std::stod(row[patch_energy.Column("loop weight (J)")]) >= 0.0);
      if (flag == 1)
      {
        fallback_rows++;
        CHECK(share > SurfaceResponseOperator::kLoopFitFallbackShare);
        CHECK(std::stod(row[patch_energy.Column("fabricated surface energy[4] (J)")]) ==
              0.0);
        CHECK(std::stod(row[patch_energy.Column("domain correction (J)")]) == 0.0);
      }
      else
      {
        CHECK(share <= SurfaceResponseOperator::kLoopFitFallbackShare);
      }
    }
    CHECK(fallback_rows == 2 * fit_counts.at(test_case.mode)[1]);  // evaluations 0 and 1
  }
#endif
}

// PatchEnergy's spatial-trace export on the Maxwell path (decision 922 B-3; the d919 review
// MAJOR-3): the Legacy spatial-cluster fixture (the offset corner pair: one spatial
// edge-cluster patch, six corner patches, the isolated-edge cells) under an Eigenmode
// PostOperator with "PatchEnergy" on a complex field: surface-response-traces.csv carries
// one row per coefficient of every spatial trace of GetMaxwellResponse(E, omega, true): the
// 1-based placed patch index, the model, the 1-based coefficient, the conductor-state flag
// (a coefficient past the contour size), the real and imaginary values in volts, the mode
// and excitation keys, and nothing else; the patch-energy table carries the patch set twice
// (evaluations 0 / 1).
TEST_CASE_METHOD(test::SurfaceResponseFiles,
                 "Maxwell spatial-trace export of a spatial cluster (decision 922 B-3)",
                 "[postoperator][surfaceresponseoperator][cornerarmtrim][3d][maxwell]"
                 "[patchenergy][spatial][Serial][Parallel]")
{
#if !defined(MFEM_USE_GSLIB)
  SKIP("SurfaceResponseOperator requires MFEM_USE_GSLIB");
#else
  MPI_Comm comm = Mpi::World();
  const fs::path output_dir = temp.temp_dir / "spatial-traces-maxwell";
  json config = ConvexMaxwellIslandConfig();
  config["Model"]["L0"] = 1.0;
  config["Problem"]["Output"] = output_dir.string();
  config["Boundaries"]["Ground"]["Attributes"] = {1, 2, 3, 4, 5, 6, 9, 10};
  config["Boundaries"]["Postprocessing"]["Dielectric"][0]["Attributes"] = {9, 10};
  auto &correction = config["Solver"]["SurfaceResponseCorrection"];
  correction["Library"] = spatial_cluster_library_3d_path.string();
  correction["CorrectionMode"] = "PostprocessOnly";
  correction["PatchEnergy"] = true;
  IoData iodata(config, false);
  iodata.boundaries.cracked_attributes.insert(9);
  iodata.boundaries.cracked_attributes.insert(10);
  MakeOutputFolder(iodata, comm);
  std::vector<std::unique_ptr<Mesh>> meshes;
  meshes.push_back(std::make_unique<Mesh>(MakeOffsetCornerPairMesh()));
  SpaceOperator space_op(iodata, meshes);
  PostOperator<ProblemType::EIGENMODE> post_op(iodata, space_op);
  auto *response = post_op.GetSurfaceResponseOperator();
  REQUIRE(response);
  REQUIRE(response->GetPatchCount() > 0);

  GridFunction field(space_op.GetNDSpace(), true);
  mfem::VectorConstantCoefficient real_coefficient = ConstantFieldCoefficient();
  mfem::Vector imaginary_value(3);
  imaginary_value[0] = -0.3;
  imaginary_value[1] = 0.5;
  imaginary_value[2] = 0.2;
  mfem::VectorConstantCoefficient imaginary_coefficient(imaginary_value);
  field.Real().ProjectCoefficient(real_coefficient);
  field.Imag().ProjectCoefficient(imaginary_coefficient);
  const std::complex<double> omega(1.0, 0.0);
  const auto result = response->GetMaxwellResponse(field, omega, true);
  REQUIRE(!result.spatial_traces.empty());
  ComplexVector e(space_op.GetNDSpace().GetTrueVSize());
  field.Real().GetTrueDofs(e.Real());
  field.Imag().GetTrueDofs(e.Imag());
  ComplexVector b(space_op.GetRTSpace().GetTrueVSize());
  b = 0.0;
  post_op.MeasureAndPrintAll(0, e, b, omega, 0.0, 0.0, 1);
  Mpi::Barrier(comm);
  if (!Mpi::Root(comm))
  {
    return;
  }
  auto Trim = [](std::string text)
  {
    const auto begin = text.find_first_not_of(" \t\r");
    const auto end = text.find_last_not_of(" \t\r");
    return begin == std::string::npos ? std::string() : text.substr(begin, end - begin + 1);
  };
  auto Int = [](const std::string &text)
  { return static_cast<int>(std::lround(std::stod(text))); };
  auto ReadTrimmedCsv = [&](const fs::path &path)
  {
    TextCsv csv = ReadTextCsv(path);
    for (auto &name : csv.header)
    {
      name = Trim(name);
    }
    for (auto &row : csv.rows)
    {
      for (auto &cell : row)
      {
        cell = Trim(cell);
      }
    }
    return csv;
  };
  auto Volt = [&](double value)
  { return iodata.units.Dimensionalize<Units::ValueType::VOLTAGE>(value); };
  REQUIRE(fs::is_regular_file(output_dir / "surface-response-traces.csv"));
  const TextCsv traces = ReadTrimmedCsv(output_dir / "surface-response-traces.csv");
  CHECK(traces.header == std::vector<std::string>{"i", "patch", "model", "coefficient",
                                                  "conductor state", "value (V)",
                                                  "value imag (V)", "m", "exc"});
  std::size_t expected_rows = 0;
  for (const auto &trace : result.spatial_traces)
  {
    expected_rows += trace.coefficients.size();
  }
  REQUIRE(traces.rows.size() == expected_rows);
  std::size_t row = 0;
  for (const auto &trace : result.spatial_traces)
  {
    for (std::size_t k = 0; k < trace.coefficients.size(); k++, row++)
    {
      const auto &cells = traces.rows[row];
      CHECK(Int(cells[traces.Column("i")]) == 1);
      CHECK(Int(cells[traces.Column("patch")]) == trace.patch + 1);
      CHECK(Int(cells[traces.Column("model")]) == trace.model);
      CHECK(Int(cells[traces.Column("coefficient")]) == static_cast<int>(k) + 1);
      CHECK(Int(cells[traces.Column("conductor state")]) ==
            (static_cast<int>(k) >= trace.contour_size ? 1 : 0));
      CHECK_THAT(
          std::stod(cells[traces.Column("value (V)")]),
          WithinAbs(Volt(trace.coefficients[k].real()),
                    1.0e-10 * std::max(1.0, std::abs(Volt(trace.coefficients[k].real())))));
      CHECK_THAT(
          std::stod(cells[traces.Column("value imag (V)")]),
          WithinAbs(Volt(trace.coefficients[k].imag()),
                    1.0e-10 * std::max(1.0, std::abs(Volt(trace.coefficients[k].imag())))));
      CHECK(Int(cells[traces.Column("m")]) == 1);
      CHECK(Int(cells[traces.Column("exc")]) == 0);
    }
  }
  // A spatial trace per spatial patch, fewer than the applied patches (the translational
  // cells carry none), and the imaginary part is a genuine column (the field is complex).
  CHECK(result.spatial_traces.size() < static_cast<std::size_t>(response->GetPatchCount()));
  bool imaginary_nonzero = false;
  for (const auto &trace : result.spatial_traces)
  {
    for (const auto &coefficient : trace.coefficients)
    {
      imaginary_nonzero = imaginary_nonzero || std::abs(coefficient.imag()) > 1.0e-12;
    }
  }
  CHECK(imaginary_nonzero);
  REQUIRE(fs::is_regular_file(output_dir / "surface-response-patch-energy.csv"));
  const TextCsv patch_energy =
      ReadTrimmedCsv(output_dir / "surface-response-patch-energy.csv");
  CHECK(patch_energy.rows.size() ==
        2 * static_cast<std::size_t>(response->GetPatchCount()));
#endif
}

// OI-ZERO-IFACE-NAN (decisions 891 / 922 C): a postprocessing Dielectric interface that is
// NOT a response target (S6's L2 ground attr 126, no EdgeDistances) read NaN in every
// corrected column of the Maxwell surface-Q-corrected.csv on main (the measurement
// initialised every interface to NaN and filled the targets), where the electrostatic
// driver carries such an interface's RAW energy forward. With the fix the untargeted row
// reads corrected = fixed-flux = raw, spread 0, matched / corner fractions 0 / 0 and the
// "response target" flag 0 (the column exists only when such an interface exists); the
// target row's values are those of the same fixture without the extra interface, column by
// column by name (the self-consistent columns, NaN when the mode is PostprocessOnly on both
// paths, are the one "not computed" class left as is).
TEST_CASE_METHOD(test::SurfaceResponseFiles,
                 "Maxwell untargeted postprocessing interface carries its raw energy "
                 "(decision 922 C)",
                 "[postoperator][surfaceresponseoperator][cornerarmtrim][3d][maxwell]"
                 "[untargeted][Serial][Parallel]")
{
#if !defined(MFEM_USE_GSLIB)
  SKIP("SurfaceResponseOperator requires MFEM_USE_GSLIB");
#else
  const fs::path mesh_path = temp.temp_dir / "untargeted-lead.mesh";
  const auto basis_path = temp.temp_dir / "untargeted-basis-points.csv";
  const auto edge_library_path =
      temp.temp_dir / "fabrication-process-untargeted-edges.json";
  MPI_Comm comm = Mpi::World();
  if (Mpi::Root(comm))
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
    library["Name"] = "unit-test-process-untargeted-edges";
    json models = json::array();
    for (auto &model : library["Models"])
    {
      if (model["Topology"] == "IsolatedEdge")
      {
        model["Interfaces"] = {{{"Type", "SA"}, {"Coupon", 1}}};
        model["BasisPoints"] = basis_path.string();
        models.push_back(model);
      }
    }
    library["Models"] = std::move(models);
    std::ofstream output(edge_library_path);
    output << library.dump(2) << "\n";
  }
  Mpi::Barrier(comm);
  json config = IslandConfig();
  config["Model"]["Mesh"] = mesh_path.string();
  config["Model"]["L0"] = 1.0;
  config["Problem"]["Type"] = "Eigenmode";
  config["Boundaries"]["Ground"]["Attributes"] = {4, 9};
  config["Boundaries"].erase("Terminal");
  config["Solver"] = {{"Order", 1},
                      {"Eigenmode", {{"Target", 1.0}}},
                      {"SurfaceResponseCorrection",
                       {{"Library", edge_library_path.string()},
                        {"TargetInterfaces", {4}},
                        {"UnmatchedPolicy", "Warn"},
                        {"CorrectionMode", "PostprocessOnly"}}}};
  auto Run = [&](bool extra_interface, const char *stem)
  {
    json run = config;
    const fs::path output_dir = temp.temp_dir / stem;
    run["Problem"]["Output"] = output_dir.string();
    if (extra_interface)
    {
      // The sheet (attribute 9) once more as an untargeted SA interface (S6's [2] / [3]
      // share their attributes the same way): no EdgeDistances, not in TargetInterfaces, a
      // non-zero raw energy (an MS entry of this fixture's sheet reads 0: the zero-energy
      // class, handled the same way).
      run["Boundaries"]["Postprocessing"]["Dielectric"].push_back({{"Index", 5},
                                                                   {"Attributes", {9}},
                                                                   {"Type", "SA"},
                                                                   {"Thickness", 0.002},
                                                                   {"Permittivity", 4.0},
                                                                   {"LossTan", 0.002}});
    }
    IoData iodata(run, false);
    iodata.boundaries.cracked_attributes.insert(9);
    MakeOutputFolder(iodata, comm);
    std::vector<std::unique_ptr<Mesh>> meshes;
    meshes.push_back(std::make_unique<Mesh>(MakeLeadMesh(0.0, 0.75)));
    SpaceOperator space_op(iodata, meshes);
    PostOperator<ProblemType::EIGENMODE> post_op(iodata, space_op);
    REQUIRE(post_op.GetSurfaceResponseOperator());
    mfem::VectorFunctionCoefficient rotational(3,
                                               [](const mfem::Vector &x, mfem::Vector &v)
                                               {
                                                 v.SetSize(3);
                                                 v(0) = 0.4 * (x[1] - 0.5) + 0.1 * x[2];
                                                 v(1) = -0.4 * (x[0] - 0.5) + 0.3;
                                                 v(2) = 0.2 * x[0] * (x[1] - 0.5);
                                               });
    GridFunction field(space_op.GetNDSpace(), true);
    field.Real().ProjectCoefficient(rotational);
    field.Imag() = 0.0;
    ComplexVector e(space_op.GetNDSpace().GetTrueVSize());
    field.Real().GetTrueDofs(e.Real());
    e.Imag() = 0.0;
    ComplexVector b(space_op.GetRTSpace().GetTrueVSize());
    b = 0.0;
    post_op.MeasureAndPrintAll(0, e, b, std::complex<double>(1.0, 0.0), 0.0, 0.0, 1);
    Mpi::Barrier(comm);
    return output_dir;
  };
  const fs::path with_dir = Run(true, "untargeted-with");
  const fs::path without_dir = Run(false, "untargeted-without");
  if (!Mpi::Root(comm))
  {
    return;
  }
  auto Trim = [](std::string text)
  {
    const auto begin = text.find_first_not_of(" \t\r");
    const auto end = text.find_last_not_of(" \t\r");
    return begin == std::string::npos ? std::string() : text.substr(begin, end - begin + 1);
  };
  auto ReadTrimmedCsv = [&](const fs::path &path)
  {
    TextCsv csv = ReadTextCsv(path);
    for (auto &name : csv.header)
    {
      name = Trim(name);
    }
    for (auto &row : csv.rows)
    {
      for (auto &cell : row)
      {
        cell = Trim(cell);
      }
    }
    return csv;
  };
  const TextCsv with = ReadTrimmedCsv(with_dir / "surface-Q-corrected.csv");
  const TextCsv without = ReadTrimmedCsv(without_dir / "surface-Q-corrected.csv");
  // The flag column exists only with the untargeted interface; the other columns are the
  // same set in the same order.
  REQUIRE(with.header.size() == without.header.size() + 1);
  CHECK(with.header.back() == "response target");
  CHECK(std::vector<std::string>(with.header.begin(), with.header.end() - 1) ==
        without.header);
  REQUIRE(without.rows.size() == 1);
  REQUIRE(with.rows.size() == 2);
  auto Row = [&](const TextCsv &csv, int interface) -> const std::vector<std::string> &
  {
    for (const auto &row : csv.rows)
    {
      if (std::lround(std::stod(row[csv.Column("interface")])) == interface)
      {
        return row;
      }
    }
    FAIL("interface " << interface << " row missing");
    return csv.rows.front();
  };
  const auto &target_with = Row(with, 4);
  const auto &target_without = Row(without, 4);
  const auto &untargeted = Row(with, 5);
  // The target row: every shared column equal as text (the extra interface changes nothing
  // of the target's measurement); flag 1.
  for (const auto &name : without.header)
  {
    INFO(name);
    CHECK(target_with[with.Column(name)] == target_without[without.Column(name)]);
  }
  CHECK(std::lround(std::stod(target_with[with.Column("response target")])) == 1);
  // The untargeted row: raw carried into the corrected columns, spread 0, fractions 0 / 0,
  // flag 0, no NaN outside the self-consistent columns (PostprocessOnly).
  CHECK(std::lround(std::stod(untargeted[with.Column("response target")])) == 0);
  const double raw = std::stod(untargeted[with.Column("E_surf raw (J)")]);
  CHECK(raw > 0.0);
  CHECK_THAT(std::stod(untargeted[with.Column("E_surf corrected (J)")]),
             WithinRel(raw, 1.0e-12));
  CHECK_THAT(std::stod(untargeted[with.Column("E_surf corrected fixed-flux (J)")]),
             WithinRel(raw, 1.0e-12));
  CHECK(std::stod(untargeted[with.Column("trace closure spread")]) == 0.0);
  CHECK(std::stod(untargeted[with.Column("matched edge fraction")]) == 0.0);
  CHECK(std::stod(untargeted[with.Column("unmodeled corner neighborhood fraction")]) ==
        0.0);
  // p_surf corrected = raw / E_norm corrected (the corrected normalisation of the run).
  CHECK_THAT(
      std::stod(untargeted[with.Column("p_surf corrected")]),
      WithinRel(raw / std::stod(untargeted[with.Column("E_norm corrected (J)")]), 1.0e-9));
  for (const auto &name : with.header)
  {
    if (name.find("self-consistent") != std::string::npos)
    {
      continue;
    }
    INFO(name);
    CHECK(!std::isnan(std::stod(untargeted[with.Column(name)])));
  }
#endif
}

// OI-SCP-HOLE through the operator, the preflight and both postoperators (DESIGN A.4 /
// A-U2; decisions 922 / 934 / 939; the d919 review MAJOR-2). The offset-corner-pair mesh in
// the identification-based construction (the devices' form; a Legacy corner patch carries
// no feature and no raw footprint) with a Signature-keyed copy of the spatial-cluster
// library built from the mesh's own identification (one model per feature hash: the
// fixture model of that topology with the feature's Signature). The SpatialClusterPriority
// rule tests the placed BASIS-POINT boxes (CollectSpatialSupports), so the straddling
// fixture moves ONE of the synthetic cluster coupon's four basis knots to the device point
// (1.12, 0.52, 1.70): the signature, edges, frame, claims and weight are the control's
// (asserted), the cluster's basis box now reaches the basis box of island 2's top-left
// convex corner (1.125, 0.5, 1.75), and the rule drops that corner and no other (decision
// 939: the deviation from "the same basis" is named; the SCP arithmetic does not depend on
// the response values). Its window = the two arms' R windows (0.4): the +x arm (z = 1.75,
// above the box) lies wholly outside, the -z arm outside on [1.70, 1.75]: OutsideLength
// 0.25 = R + 0.05, the kept portions tiling exactly those pieces. Asserted: the operator's
// Diagnostics record and portions (Maxwell and electrostatic constructors alike), the
// preflight's record equal to the operator's with the dry-run patch table unchanged by the
// pass (form (i)), the kept pieces' energies (kept = the direct pieces, kept + the inside
// piece = the whole window: the uncovered machinery's exact partition), the Maxwell
// PostOperator's add-back (corrected = outside + models + the SpatialClusterPriority raw
// energy on both closures; the 4th table written with the operator's energy), the
// electrostatic driver's add-back (the same identity per source on surface-Q-corrected.csv;
// the palace.json record); a two-knot variant whose box holds the whole window: Count 1,
// OutsideLength 0, no portion, no table; the control (the fixture's knots): no record, no
// table — there the B4 vertex rule owns the free corners whole (this small device lies
// inside the cluster's Signature box; the known B4 hole), where the SpatialClusterPriority
// pass, which runs first, takes precedence in the straddling fixture.
TEST_CASE_METHOD(
    test::SurfaceResponseFiles,
    "SpatialClusterPriority raw remainder through the operator / the preflight "
    "/ both postoperators (decisions 922 / 934)",
    "[postoperator][surfaceresponseoperator][cornerarmtrim][3d][scphole]"
    "[Serial][Parallel]")
{
#if !defined(MFEM_USE_GSLIB)
  SKIP("SurfaceResponseOperator requires MFEM_USE_GSLIB");
#else
  constexpr double R = 0.2;
  // The identification of this mesh clusters the two facing corners WITH the two short
  // edges' far corners (six edges), so the free convex corners are the two left corners of
  // island 1 and the two top corners of island 2. The dropped corner is island 2's top-left
  // corner (1.125, 0.5, 1.75) in (x, y, z): its +x arm along the top edge z = 1.75 and its
  // -z arm along the left edge x = 1.125. The cluster's first basis knot is moved to the
  // device point (1.12, 0.52, 1.70): the cluster's basis-point box then spans x in [1.12,
  // ~1.215], z up to 1.70, inside the corner's basis-point box (x in [1.145, 1.205], z in
  // [1.67, 1.73]: the points lie 0.02-0.08 into the metal quadrant); the +x arm (z = 1.75,
  // above the box) is wholly outside, the -z arm outside on [1.70, 1.75].
  const std::array<double, 3> corner = {1.125, 0.5, 1.75};
  const std::array<double, 3> first_arm_end = {1.325, 0.5, 1.75};
  const std::array<double, 3> second_arm_end = {1.125, 0.5, 1.55};
  const std::vector<std::array<double, 3>> straddling_knots = {{1.12, 0.52, 1.70}};
  // Two knots whose box contains the whole window (x in [1.125, 1.325], z in [1.55, 1.75]).
  const std::vector<std::array<double, 3>> inside_knots = {{1.12, 0.52, 1.76},
                                                           {1.33, 0.52, 1.76}};
  const double box_face = straddling_knots[0][2];
  const std::array<double, 3> box_face_point = {corner[0], corner[1], box_face};
  const double outside_length = R + (corner[2] - box_face);
  MPI_Comm comm = Mpi::World();
  const fs::path mesh_path = temp.temp_dir / "scp-offset-corner-pair.mesh";
  const fs::path control_library = temp.temp_dir / "surface-process-scp-control.json";
  const fs::path straddling_library = temp.temp_dir / "surface-process-scp-straddling.json";
  const fs::path inside_library = temp.temp_dir / "surface-process-scp-inside.json";
  // The (ii-a) reference (decisions 958 (b) / 969): the same Signature-keyed library
  // WITHOUT the cluster model, so that the dropped corner's coupon is APPLIED (weight 1)
  // and its PatchEnergy row is the fabricated energy the diagnostic column must reproduce.
  const fs::path corner_only_library =
      temp.temp_dir / "surface-process-scp-corner-only.json";
  std::vector<std::unique_ptr<Mesh>> meshes;
  meshes.push_back(std::make_unique<Mesh>(MakeOffsetCornerPairMesh()));
  // The Maxwell configuration of the spatial-cluster case in the identification-based
  // construction (the devices' form: a Legacy corner patch carries no feature and no raw
  // footprint).
  auto MaxwellConfig = [&](const fs::path &library, const fs::path &output)
  {
    json config = ConvexMaxwellIslandConfig();
    config["Model"]["Mesh"] = mesh_path.string();
    config["Model"]["L0"] = 1.0;
    config["Problem"]["Output"] = output.string();
    config["Boundaries"]["Ground"]["Attributes"] = {1, 2, 3, 4, 5, 6, 9, 10};
    config["Boundaries"]["Postprocessing"]["Dielectric"][0]["Attributes"] = {9, 10};
    auto &correction = config["Solver"]["SurfaceResponseCorrection"];
    correction["Library"] = library.string();
    correction["CorrectionMode"] = "PostprocessOnly";
    correction.erase("PatchConstruction");
    return config;
  };
  auto MaxwellIoData = [&](json config)
  {
    IoData iodata(std::move(config), false);
    iodata.boundaries.cracked_attributes.insert(9);
    iodata.boundaries.cracked_attributes.insert(10);
    return iodata;
  };
  // The identification of the mesh (every feature Missing against the fixture's Edges-keyed
  // library): the Signature-keyed libraries below are built from its features.
  const fs::path identification_dir = temp.temp_dir / "scp-identification";
  {
    auto config = MaxwellConfig(spatial_cluster_library_3d_path, identification_dir);
    config["Solver"]["SurfaceResponseCorrection"]["UnmatchedPolicy"] = "Warn";
    IoData iodata = MaxwellIoData(config);
    MakeOutputFolder(iodata, comm);
    WriteSurfaceResponseRequirements(
        iodata, *meshes.back(),
        (identification_dir / "surface-response-requirements.json").string());
  }
  Mpi::Barrier(comm);
  if (Mpi::Root(comm))
  {
    {
      mfem::Mesh serial = MakeSerialOffsetCornerPairMesh();
      std::ofstream output(mesh_path);
      serial.Print(output);
    }
    std::ifstream manifest_input(identification_dir / "surface-response-requirements.json");
    REQUIRE(manifest_input);
    const json manifest = json::parse(manifest_input);
    const auto &features = manifest["Identification"]["Features"];
    std::ifstream library_input(spatial_cluster_library_3d_path);
    REQUIRE(library_input);
    const json fixture_library = json::parse(library_input);
    REQUIRE(fixture_library["Models"].size() == 3);
    REQUIRE(fixture_library["Models"][0]["Topology"] == "IsolatedEdge");
    REQUIRE(fixture_library["Models"][1]["Topology"] == "ConvexCorner");
    REQUIRE(fixture_library["Models"][2]["Name"] == "offset-corner-pair");
    // The cluster feature's canonical frame: a Signature-keyed model's BasisPoints are
    // local coordinates of this frame (origin + sum_k local[k] axes[k]).
    const auto cluster_feature =
        std::find_if(features.begin(), features.end(), [](const auto &feature)
                     { return feature["Type"] == "SpatialEdgeCluster"; });
    REQUIRE(cluster_feature != features.end());
    const auto frame_origin =
        (*cluster_feature)["Frame"]["Origin"].get<std::array<double, 3>>();
    const auto frame_axes =
        (*cluster_feature)["Frame"]["Axes"].get<std::array<std::array<double, 3>, 3>>();
    auto Local = [&](const std::array<double, 3> &device)
    {
      std::array<double, 3> local{};
      for (int k = 0; k < 3; k++)
      {
        for (int d = 0; d < 3; d++)
        {
          local[k] += (device[d] - frame_origin[d]) * frame_axes[k][d];
        }
      }
      return local;
    };
    // A Signature-keyed copy of the fixture library: one model per distinct feature hash,
    // the fixture model of the feature's topology with the feature's Signature (the
    // corner's angle from it); the cluster model keyed by its Signature instead of its
    // Edges, its first basis knots replaced by `knots` (device coordinates; empty = the
    // fixture's).
    auto WriteLibrary = [&](const std::vector<std::array<double, 3>> &knots,
                            const fs::path &path, const char *suffix,
                            bool with_cluster = true)
    {
      json library = fixture_library;
      library["Name"] = std::string("unit-test-process-spatial-cluster-scp-") + suffix;
      library["Models"] = json::array();
      std::set<std::string> hashes;
      for (const auto &feature : features)
      {
        if (!hashes.insert(feature["Hash"].get<std::string>()).second)
        {
          continue;
        }
        const std::string type = feature["Type"];
        REQUIRE((type == "IsolatedEdge" || type == "ConvexCorner" ||
                 type == "SpatialEdgeCluster"));
        if (type == "SpatialEdgeCluster" && !with_cluster)
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
          model["Name"] = std::string("offset-corner-pair-scp-") + suffix;
          if (!knots.empty())
          {
            const std::vector<std::array<double, 3>> fixture_knots = {
                {0.025, -0.025, 0.02},
                {0.050, -0.050, 0.02},
                {0.075, -0.075, -0.02},
                {0.100, -0.100, -0.02}};
            REQUIRE(knots.size() < fixture_knots.size());
            const fs::path points_path =
                temp.temp_dir /
                (std::string("spatial-cluster-scp-") + suffix + "-points.csv");
            std::ofstream output(points_path);
            output << "x,y,z\n";
            for (std::size_t k = 0; k < fixture_knots.size(); k++)
            {
              const auto local = k < knots.size() ? Local(knots[k]) : fixture_knots[k];
              output << std::setprecision(17) << local[0] << "," << local[1] << ","
                     << local[2] << "\n";
            }
            model["BasisPoints"] = points_path.string();
          }
          // The SupportPoints keep containing every knot (they are not the box the
          // ownership rules read: that is the placed basis points' box).
          json support = json::array();
          for (const double x : {-1.2, 1.2})
          {
            for (const double y : {-1.2, 1.2})
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
    WriteLibrary({}, control_library, "control");
    WriteLibrary(straddling_knots, straddling_library, "straddling");
    WriteLibrary(inside_knots, inside_library, "inside");
    WriteLibrary({}, corner_only_library, "corner-only", false);
  }
  Mpi::Barrier(comm);

  auto Near = [](const std::array<double, 3> &a, const std::array<double, 3> &b)
  {
    return std::abs(a[0] - b[0]) < 1.0e-9 && std::abs(a[1] - b[1]) < 1.0e-9 &&
           std::abs(a[2] - b[2]) < 1.0e-9;
  };
  auto Perimeter = [](const auto &portions)
  {
    std::vector<SurfacePostOperator::UncoveredPerimeterPortion> perimeter;
    for (const auto &portion : portions)
    {
      perimeter.push_back({portion.p0, portion.p1, portion.topology, portion.feature});
    }
    return perimeter;
  };
  // The record of a placement (the operator's Diagnostics or the preflight's manifest
  // entry) against the analytic exclusion.
  auto CheckRecord = [&](const json &record, double expected_outside)
  {
    REQUIRE(record["Count"].get<int>() == 1);
    CHECK_THAT(record["WindowLength"].get<double>(), WithinAbs(2.0 * R, 1.0e-9));
    CHECK_THAT(record["OutsideLength"].get<double>(), WithinAbs(expected_outside, 1.0e-9));
    REQUIRE(record["Patches"].size() == 1);
    const auto &entry = record["Patches"][0];
    CHECK(entry["Model"] == "convex-corner-90");
    CHECK(entry["Topology"] == "convex corner");
    CHECK_THAT(entry["Origin"][0].get<double>(), WithinAbs(corner[0], 1.0e-9));
    CHECK_THAT(entry["Origin"][1].get<double>(), WithinAbs(corner[1], 1.0e-9));
    CHECK_THAT(entry["Origin"][2].get<double>(), WithinAbs(corner[2], 1.0e-9));
    CHECK_THAT(entry["OutsideLength"].get<double>(), WithinAbs(expected_outside, 1.0e-9));
    CHECK(entry["Owner"]["Model"].get<std::string>().rfind("offset-corner-pair-scp-", 0) ==
          0);
    if (expected_outside > 0.0)
    {
      CHECK(record["RawPortions"]["Count"].get<int>() >= 2);
      CHECK_THAT(record["RawPortions"]["Length"].get<double>(),
                 WithinAbs(expected_outside, 1.0e-9));
    }
    else
    {
      CHECK(record["RawPortions"]["Count"].get<int>() == 0);
    }
  };
  // The operator's portions tile the kept pieces exactly: the whole +x arm (as the
  // identification cut it, over its perimeter segments) and the -z arm from the vertex to
  // the box face; nothing else.
  auto OnArm = [&](const std::array<double, 3> &p, const std::array<double, 3> &end)
  {
    for (int d = 0; d < 3; d++)
    {
      const double lo = std::min(corner[d], end[d]) - 1.0e-9;
      const double hi = std::max(corner[d], end[d]) + 1.0e-9;
      if (p[d] < lo || p[d] > hi)
      {
        return false;
      }
    }
    return true;
  };
  auto CheckPortions = [&](const SurfaceResponseOperator &response)
  {
    const auto &portions = response.GetSpatialClusterPriorityPortions();
    REQUIRE(portions.size() >= 2);
    double first_arm_length = 0.0, second_arm_length = 0.0;
    for (const auto &portion : portions)
    {
      CHECK(portion.topology == "convex corner");
      CHECK(portion.feature == portions.front().feature);
      double length = 0.0;
      for (int d = 0; d < 3; d++)
      {
        length += (portion.p1[d] - portion.p0[d]) * (portion.p1[d] - portion.p0[d]);
      }
      length = std::sqrt(length);
      if (OnArm(portion.p0, first_arm_end) && OnArm(portion.p1, first_arm_end))
      {
        first_arm_length += length;
      }
      else if (OnArm(portion.p0, box_face_point) && OnArm(portion.p1, box_face_point))
      {
        second_arm_length += length;
      }
      else
      {
        FAIL("a SpatialClusterPriority portion off the kept pieces");
      }
    }
    CHECK_THAT(first_arm_length, WithinAbs(R, 1.0e-9));
    CHECK_THAT(second_arm_length, WithinAbs(corner[2] - box_face, 1.0e-9));
    const auto &exclusions = response.GetSpatialClusterPriorityExclusions();
    REQUIRE(exclusions.entries.size() == 1);
    CHECK_THAT(exclusions.entries[0].window_length, WithinAbs(2.0 * R, 1.0e-9));
    CHECK_THAT(exclusions.entries[0].outside_length, WithinAbs(outside_length, 1.0e-9));
    CHECK_THAT(exclusions.portions.length, WithinAbs(outside_length, 1.0e-9));
    CHECK(exclusions.portions.portions.size() == portions.size());
  };
  auto AppliedAtCorner = [&](const SurfaceResponseOperator &response)
  {
    int count = 0;
    for (const auto &assignment : response.GetPatchAssignments())
    {
      count += Near(assignment.origin, corner) ? 1 : 0;
    }
    Mpi::GlobalSum(1, &count, Mpi::World());
    return count;
  };

  // (1) The Maxwell constructor and its postoperator.
  const fs::path straddling_dir = temp.temp_dir / "scp-maxwell-straddling";
  const fs::path fixture_dir = temp.temp_dir / "scp-maxwell-control";
  const fs::path inside_dir = temp.temp_dir / "scp-maxwell-inside";
  IoData straddling_iodata =
      MaxwellIoData(MaxwellConfig(straddling_library, straddling_dir));
  IoData fixture_iodata = MaxwellIoData(MaxwellConfig(control_library, fixture_dir));
  IoData inside_iodata = MaxwellIoData(MaxwellConfig(inside_library, inside_dir));
  MakeOutputFolder(straddling_iodata, comm);
  MakeOutputFolder(fixture_iodata, comm);
  MakeOutputFolder(inside_iodata, comm);
  SpaceOperator space_op(straddling_iodata, meshes);
  PostOperator<ProblemType::EIGENMODE> post_op(straddling_iodata, space_op);
  auto *response = post_op.GetSurfaceResponseOperator();
  REQUIRE(response);
  const json statistics = response->GetStatistics();
  REQUIRE(statistics["Diagnostics"].contains("SpatialClusterPriority"));
  CheckRecord(statistics["Diagnostics"]["SpatialClusterPriority"], outside_length);
  CheckPortions(*response);
  CHECK(AppliedAtCorner(*response) == 0);
  const auto &cluster_record =
      statistics["Diagnostics"]["SpatialClusterPriority"]["Patches"][0]["Owner"]["Box"];
  CHECK_THAT(cluster_record["Max"][2].get<double>(), WithinAbs(box_face, 1.0e-9));
  CHECK_THAT(cluster_record["Min"][0].get<double>(),
             WithinAbs(straddling_knots[0][0], 1.0e-9));
  CHECK(cluster_record["Max"][0].get<double>() > corner[0]);

  std::vector<std::unique_ptr<Mesh>> fixture_meshes;
  fixture_meshes.push_back(std::make_unique<Mesh>(MakeOffsetCornerPairMesh()));
  SpaceOperator fixture_space_op(fixture_iodata, fixture_meshes);
  PostOperator<ProblemType::EIGENMODE> fixture_post_op(fixture_iodata, fixture_space_op);
  auto *fixture_response = fixture_post_op.GetSurfaceResponseOperator();
  REQUIRE(fixture_response);
  const json fixture_statistics = fixture_response->GetStatistics();
  CHECK(!fixture_statistics["Diagnostics"].contains("SpatialClusterPriority"));
  CHECK(fixture_response->GetSpatialClusterPriorityPortions().empty());
  // In the control the corner is not applied either: this small device lies inside the
  // cluster's Signature box, so the continuation ownership's vertex rule (B4) owns every
  // free corner whole (no raw remainder: the owner coupon models its window). The
  // SpatialClusterPriority pass runs before it and takes precedence where its basis-point
  // boxes overlap (the straddling case above: the corner is in the SpatialClusterPriority
  // record, its outside window kept raw).
  CHECK(AppliedAtCorner(*fixture_response) == 0);
  {
    const auto &vertices =
        fixture_statistics["Diagnostics"]["ContinuationOwnership"]["Vertices"];
    int owned_corner = 0;
    for (const auto &vertex : vertices["Records"])
    {
      const auto origin = vertex["Origin"].get<std::array<double, 3>>();
      owned_corner += Near(origin, corner) ? 1 : 0;
    }
    CHECK(owned_corner == 1);
  }
  // The same cluster patch (one spatial edge-cluster, the same frame and weight) in both:
  // the owning rank's entry, replicated by the sum (zero elsewhere).
  auto ClusterAssignment = [&](const SurfaceResponseOperator &r, const std::string &model)
  {
    const auto names = r.GetModelNames();
    std::array<double, 10> values{};
    int count = 0;
    for (const auto &assignment : r.GetPatchAssignments())
    {
      if (names.at(assignment.model).rfind(model, 0) == 0)
      {
        for (int d = 0; d < 3; d++)
        {
          values[d] = assignment.origin[d];
          values[3 + d] = assignment.axis_u[d];
          values[6 + d] = assignment.axis_v[d];
        }
        values[9] = assignment.weight;
        count++;
      }
    }
    Mpi::GlobalSum(1, &count, Mpi::World());
    Mpi::GlobalSum(static_cast<int>(values.size()), values.data(), Mpi::World());
    REQUIRE(count == 1);
    return values;
  };
  {
    const auto wide = ClusterAssignment(*response, "offset-corner-pair-scp-straddling");
    const auto fixture =
        ClusterAssignment(*fixture_response, "offset-corner-pair-scp-control");
    for (std::size_t i = 0; i < wide.size(); i++)
    {
      CHECK_THAT(wide[i], WithinAbs(fixture[i], 1.0e-12));
    }
  }

  std::vector<std::unique_ptr<Mesh>> inside_meshes;
  inside_meshes.push_back(std::make_unique<Mesh>(MakeOffsetCornerPairMesh()));
  SpaceOperator inside_space_op(inside_iodata, inside_meshes);
  PostOperator<ProblemType::EIGENMODE> inside_post_op(inside_iodata, inside_space_op);
  auto *inside_response = inside_post_op.GetSurfaceResponseOperator();
  REQUIRE(inside_response);
  const json inside_statistics = inside_response->GetStatistics();
  REQUIRE(inside_statistics["Diagnostics"].contains("SpatialClusterPriority"));
  CheckRecord(inside_statistics["Diagnostics"]["SpatialClusterPriority"], 0.0);
  CHECK(inside_response->GetSpatialClusterPriorityPortions().empty());
  CHECK(AppliedAtCorner(*inside_response) == 0);

  // (2) The preflight (A-U2): the same record on the dry run; the dry-run patch table is
  // not changed by the pass (form (i)): it lists the dropped corner as it lists every built
  // patch, its row identical to the control's (where the B4 vertex rule owns it).
  auto DryRunCornerRow = [&](const fs::path &dir)
  {
    const TextCsv table = ReadTextCsv(dir / "surface-response-patches.csv");
    std::vector<std::string> corner_row;
    int corner_rows = 0;
    for (const auto &row : table.rows)
    {
      const std::array<double, 3> origin = {std::stod(row[table.Column("OriginX")]),
                                            std::stod(row[table.Column("OriginY")]),
                                            std::stod(row[table.Column("OriginZ")])};
      if (Near(origin, corner) && row[table.Column("Topology")] == "convex corner")
      {
        corner_rows++;
        corner_row = row;
      }
    }
    CHECK(corner_rows == 1);
    return std::make_pair(table.rows.size(), corner_row);
  };
  {
    const fs::path manifest_path = straddling_dir / "surface-response-requirements.json";
    WriteSurfaceResponseRequirements(straddling_iodata, *meshes.back(),
                                     manifest_path.string());
    const fs::path fixture_manifest_path =
        fixture_dir / "surface-response-requirements.json";
    WriteSurfaceResponseRequirements(fixture_iodata, *fixture_meshes.back(),
                                     fixture_manifest_path.string());
    Mpi::Barrier(comm);
    if (Mpi::Root(comm))
    {
      std::ifstream input(manifest_path);
      REQUIRE(input);
      const json manifest = json::parse(input);
      REQUIRE(manifest["Identification"]["Diagnostics"].contains("SpatialClusterPriority"));
      const auto &record =
          manifest["Identification"]["Diagnostics"]["SpatialClusterPriority"];
      CheckRecord(record, outside_length);
      const auto &operator_record = statistics["Diagnostics"]["SpatialClusterPriority"];
      CHECK_THAT(record["OutsideLength"].get<double>(),
                 WithinAbs(operator_record["OutsideLength"].get<double>(), 1.0e-12));
      CHECK_THAT(record["WindowLength"].get<double>(),
                 WithinAbs(operator_record["WindowLength"].get<double>(), 1.0e-12));
      CHECK(record["Patches"][0]["Feature"] == operator_record["Patches"][0]["Feature"]);
      CHECK(record["Patches"][0]["Owner"]["Box"] ==
            operator_record["Patches"][0]["Owner"]["Box"]);
      CHECK(record["RawPortions"]["Count"] == operator_record["RawPortions"]["Count"]);
      CHECK_THAT(
          manifest["Summary"]["SpatialClusterPriority"]["OutsideLength"].get<double>(),
          WithinAbs(outside_length, 1.0e-9));
      CHECK(manifest["Summary"]["SpatialClusterPriority"]["Patches"].get<int>() == 1);
      std::ifstream fixture_input(fixture_manifest_path);
      REQUIRE(fixture_input);
      const json fixture_manifest = json::parse(fixture_input);
      CHECK(!fixture_manifest["Identification"]["Diagnostics"].contains(
          "SpatialClusterPriority"));
      CHECK(!fixture_manifest["Summary"].contains("SpatialClusterPriority"));
      const auto [straddling_rows, straddling_corner] = DryRunCornerRow(straddling_dir);
      const auto [control_rows, control_corner] = DryRunCornerRow(fixture_dir);
      CHECK(straddling_rows == control_rows);
      CHECK(straddling_corner == control_corner);
    }
  }

  // (3) The energies of a prescribed gradient field: the kept pieces' raw within-R energy
  // (the uncovered machinery), its exact partition and the Maxwell PostOperator's add-back.
  mfem::FunctionCoefficient potential_coefficient(
      [](const mfem::Vector &x) { return std::sinh(0.7 * (x[0] - 1.0)) + 0.3 * x[2]; });
  auto Gradient = [&](auto &h1_fespace, auto &nd_fespace)
  {
    mfem::ParGridFunction potential(&h1_fespace.Get());
    potential.ProjectCoefficient(potential_coefficient);
    mfem::ParDiscreteLinearOperator gradient(&h1_fespace.Get(), &nd_fespace.Get());
    gradient.AddDomainInterpolator(new mfem::GradientInterpolator());
    gradient.Assemble();
    auto field = std::make_unique<GridFunction>(nd_fespace, false);
    gradient.Mult(potential, field->Real());
    field->Real() *= -1.0;
    Vector e(nd_fespace.GetTrueVSize());
    field->Real().GetTrueDofs(e);
    field->Real().SetFromTrueDofs(e);
    field->Real().ExchangeFaceNbrData();
    return field;
  };
  const auto field = Gradient(space_op.GetH1Space(), space_op.GetNDSpace());
  const auto &surf_post_op = post_op.GetSurfacePostOp();
  auto PieceEnergy =
      [&](const std::vector<SurfacePostOperator::UncoveredPerimeterPortion> &perimeter)
  {
    return surf_post_op.GetInterfaceUncoveredEdgeEnergies({4}, *field, nullptr, perimeter)
        .at(4)
        .energy;
  };
  const double kept = PieceEnergy(Perimeter(response->GetSpatialClusterPriorityPortions()));
  // The direct pieces, one per perimeter segment (the mesh edges of 0.125 along the island
  // edges; the uncovered machinery assigns a portion to the tree segment of its midpoint):
  // the +x arm [1.125, 1.25] + [1.25, 1.325] at z = 1.75, the -z arm's outside piece
  // [1.70, 1.75] and its inside pieces [1.625, 1.70] + [1.55, 1.625] at x = 1.125.
  const std::array<double, 3> first_arm_node = {1.25, 0.5, 1.75};
  const std::array<double, 3> second_arm_node = {1.125, 0.5, 1.625};
  auto Pieces =
      [&](std::initializer_list<std::pair<std::array<double, 3>, std::array<double, 3>>>
              ends)
  {
    std::vector<SurfacePostOperator::UncoveredPerimeterPortion> perimeter;
    for (const auto &[p0, p1] : ends)
    {
      perimeter.push_back({p0, p1, "convex corner", 0});
    }
    return perimeter;
  };
  const double first_arm_a = PieceEnergy(Pieces({{corner, first_arm_node}}));
  const double first_arm_b = PieceEnergy(Pieces({{first_arm_node, first_arm_end}}));
  const double second_arm_outside = PieceEnergy(Pieces({{corner, box_face_point}}));
  const double second_arm_inside = PieceEnergy(
      Pieces({{box_face_point, second_arm_node}, {second_arm_node, second_arm_end}}));
  const double whole_window = PieceEnergy(Pieces({{corner, first_arm_node},
                                                  {first_arm_node, first_arm_end},
                                                  {corner, second_arm_node},
                                                  {second_arm_node, second_arm_end}}));
  CHECK(kept > 0.0);
  CHECK(first_arm_a > 0.0);
  CHECK(first_arm_b > 0.0);
  CHECK(second_arm_inside > 0.0);
  CHECK(second_arm_outside >= 0.0);
  INFO("kept " << kept << " = " << first_arm_a << " + " << first_arm_b << " + "
               << second_arm_outside << "; inside " << second_arm_inside << "; window "
               << whole_window);
  CHECK_THAT(kept, WithinRel(first_arm_a + first_arm_b + second_arm_outside, 1.0e-9));
  CHECK_THAT(kept + second_arm_inside, WithinRel(whole_window, 1.0e-9));
  // The electrostatic postoperator's reading of the same pieces (the same surface
  // quadrature on the same mesh) equals the Maxwell one.
  {
    json config = IslandConfig();
    config["Model"]["Mesh"] = mesh_path.string();
    config["Model"]["L0"] = 1.0;
    config["Problem"]["Output"] = (temp.temp_dir / "scp-es-postop").string();
    config["Boundaries"]["Terminal"] = {{{"Index", 1}, {"Attributes", {9}}},
                                        {{"Index", 2}, {"Attributes", {10}}}};
    config["Boundaries"]["Postprocessing"]["Dielectric"][0]["Attributes"] = {9, 10};
    config["Solver"]["Electrostatic"]["ResponseCorrection"]["Library"] =
        straddling_library.string();
    config["Solver"]["Electrostatic"]["ResponseCorrection"].erase("PatchConstruction");
    IoData electrostatic_iodata(config, false);
    electrostatic_iodata.boundaries.cracked_attributes.insert(9);
    electrostatic_iodata.boundaries.cracked_attributes.insert(10);
    MakeOutputFolder(electrostatic_iodata, comm);
    std::vector<std::unique_ptr<Mesh>> electrostatic_meshes;
    electrostatic_meshes.push_back(std::make_unique<Mesh>(MakeOffsetCornerPairMesh()));
    LaplaceOperator laplace_op(electrostatic_iodata, electrostatic_meshes);
    SurfaceResponseOperator electrostatic_response(electrostatic_iodata, laplace_op);
    PostOperator<ProblemType::ELECTROSTATIC> electrostatic_post_op(electrostatic_iodata,
                                                                   laplace_op);
    const json electrostatic_statistics = electrostatic_response.GetStatistics();
    REQUIRE(electrostatic_statistics["Diagnostics"].contains("SpatialClusterPriority"));
    CheckRecord(electrostatic_statistics["Diagnostics"]["SpatialClusterPriority"],
                outside_length);
    CheckPortions(electrostatic_response);
    CHECK(AppliedAtCorner(electrostatic_response) == 0);
    const auto electrostatic_field =
        Gradient(laplace_op.GetH1Space(), laplace_op.GetNDSpace());
    const double electrostatic_kept =
        electrostatic_post_op.GetSurfacePostOp()
            .GetInterfaceUncoveredEdgeEnergies(
                {4}, *electrostatic_field, nullptr,
                Perimeter(electrostatic_response.GetSpatialClusterPriorityPortions()))
            .at(4)
            .energy;
    CHECK_THAT(electrostatic_kept, WithinRel(kept, 1.0e-9));
  }

  // The Maxwell PostOperator's tables on the same field: corrected = outside + models + the
  // SpatialClusterPriority raw energy on both closures (the fixture: outside + models, no
  // table).
  auto Measure = [&](PostOperator<ProblemType::EIGENMODE> &p, SpaceOperator &s,
                     const GridFunction &real_field)
  {
    GridFunction complex_field(s.GetNDSpace(), true);
    complex_field.Real() = real_field.Real();
    complex_field.Imag() = 0.0;
    ComplexVector e(s.GetNDSpace().GetTrueVSize());
    complex_field.Real().GetTrueDofs(e.Real());
    e.Imag() = 0.0;
    ComplexVector b(s.GetRTSpace().GetTrueVSize());
    b = 0.0;
    p.MeasureAndPrintAll(0, e, b, std::complex<double>(1.0, 0.0), 0.0, 0.0, 1);
  };
  Measure(post_op, space_op, *field);
  const auto fixture_field =
      Gradient(fixture_space_op.GetH1Space(), fixture_space_op.GetNDSpace());
  Measure(fixture_post_op, fixture_space_op, *fixture_field);
  const auto inside_field =
      Gradient(inside_space_op.GetH1Space(), inside_space_op.GetNDSpace());
  Measure(inside_post_op, inside_space_op, *inside_field);
  const double kept_joule =
      straddling_iodata.units.Dimensionalize<Units::ValueType::ENERGY>(kept);
  Mpi::Barrier(comm);

  // (4) The electrostatic driver on the same mesh and libraries (its add-back and record).
  const fs::path es_straddling_dir = temp.temp_dir / "scp-es-straddling";
  const fs::path es_fixture_dir = temp.temp_dir / "scp-es-control";
  const fs::path es_corner_only_dir = temp.temp_dir / "scp-es-corner-only";
  {
    json config = IslandConfig();
    config["Model"]["Mesh"] = mesh_path.string();
    config["Model"]["L0"] = 1.0;
    config["Boundaries"]["Terminal"] = {{{"Index", 1}, {"Attributes", {9}}},
                                        {{"Index", 2}, {"Attributes", {10}}}};
    config["Boundaries"]["Postprocessing"]["Dielectric"][0]["Attributes"] = {9, 10};
    auto &correction = config["Solver"]["Electrostatic"]["ResponseCorrection"];
    correction["CorrectionMode"] = "PostprocessOnly";
    correction.erase("PatchConstruction");
    config["Solver"]["Linear"] = {{"Tol", 1.0e-12}, {"MaxIts", 400}};
    config["Problem"]["Output"] = es_straddling_dir.string();
    correction["Library"] = straddling_library.string();
    RunElectrostatic(config);
    config["Problem"]["Output"] = es_fixture_dir.string();
    correction["Library"] = control_library.string();
    RunElectrostatic(config);
    // The (ii-a) reference: the corner APPLIED (no cluster model: the cluster is Missing,
    // UnmatchedPolicy Warn) with PatchEnergy, on the same raw field.
    config["Problem"]["Output"] = es_corner_only_dir.string();
    correction["Library"] = corner_only_library.string();
    correction["UnmatchedPolicy"] = "Warn";
    correction["PatchEnergy"] = true;
    RunElectrostatic(config);
  }
  Mpi::Barrier(comm);
  if (!Mpi::Root(comm))
  {
    return;
  }

  auto Trim = [](std::string text)
  {
    const auto begin = text.find_first_not_of(" \t\r");
    const auto end = text.find_last_not_of(" \t\r");
    return begin == std::string::npos ? std::string() : text.substr(begin, end - begin + 1);
  };
  auto ReadTrimmedCsv = [&](const fs::path &path)
  {
    TextCsv csv = ReadTextCsv(path);
    for (auto &name : csv.header)
    {
      name = Trim(name);
    }
    for (auto &row : csv.rows)
    {
      for (auto &cell : row)
      {
        cell = Trim(cell);
      }
    }
    return csv;
  };
  // The one row of a one-interface, one-print table.
  auto Single = [&](const fs::path &path, const std::string &column)
  {
    const TextCsv table = ReadTrimmedCsv(path);
    REQUIRE(table.rows.size() == 1);
    return std::stod(table.rows[0][table.Column(column)]);
  };
  auto ModelEnergySum = [&](const fs::path &dir, int evaluation)
  {
    const TextCsv models = ReadTrimmedCsv(dir / "surface-response-model-energy.csv");
    double sum = 0.0;
    for (const auto &row : models.rows)
    {
      if (std::lround(std::stod(row[models.Column("evaluation")])) == evaluation)
      {
        sum += std::stod(row[models.Column("fabricated surface energy[4] (J)")]);
      }
    }
    return sum;
  };
  const char *scp_table = "surface-response-spatial-cluster-priority-energy.csv";
  auto PortionTableTotal = [&](const fs::path &path, int evaluation)
  {
    const TextCsv table = ReadTrimmedCsv(path);
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
      if (std::lround(std::stod(row[table.Column("evaluation")])) == evaluation &&
          row[table.Column("type")] == "Total")
      {
        return std::stod(row[energy_column]);
      }
    }
    FAIL("no Total row in " << path.string());
    return 0.0;
  };
  {
    REQUIRE(fs::is_regular_file(straddling_dir / scp_table));
    const double scp_ft = PortionTableTotal(straddling_dir / scp_table, 0);
    const double scp_ff = PortionTableTotal(straddling_dir / scp_table, 1);
    CHECK_THAT(scp_ft, WithinRel(kept_joule, 1.0e-10));
    CHECK_THAT(scp_ff, WithinRel(kept_joule, 1.0e-10));
    const TextCsv table = ReadTrimmedCsv(straddling_dir / scp_table);
    int corner_rows = 0;
    for (const auto &row : table.rows)
    {
      if (row[table.Column("type")] == "convex corner" &&
          std::lround(std::stod(row[table.Column("evaluation")])) == 0)
      {
        corner_rows++;
        CHECK(std::lround(std::stod(row[table.Column("portions")])) ==
              static_cast<long>(response->GetSpatialClusterPriorityPortions().size()));
        CHECK_THAT(std::stod(row[table.Column("length (m)")]),
                   WithinAbs(outside_length, 1.0e-9));
      }
    }
    CHECK(corner_rows == 1);
    const double outside = Single(straddling_dir / "surface-Q-edge.csv", "E_out (J)");
    const double corrected =
        Single(straddling_dir / "surface-Q-corrected.csv", "E_surf corrected (J)");
    const double corrected_fixed_flux = Single(straddling_dir / "surface-Q-corrected.csv",
                                               "E_surf corrected fixed-flux (J)");
    CHECK_THAT(corrected,
               WithinRel(outside + ModelEnergySum(straddling_dir, 0) + scp_ft, 1.0e-10));
    CHECK_THAT(corrected_fixed_flux,
               WithinRel(outside + ModelEnergySum(straddling_dir, 1) + scp_ff, 1.0e-10));
    CHECK(!fs::exists(fixture_dir / scp_table));
    CHECK_THAT(Single(fixture_dir / "surface-Q-corrected.csv", "E_surf corrected (J)"),
               WithinRel(Single(fixture_dir / "surface-Q-edge.csv", "E_out (J)") +
                             ModelEnergySum(fixture_dir, 0),
                         1.0e-10));
    CHECK(!fs::exists(inside_dir / scp_table));
    CHECK_THAT(Single(inside_dir / "surface-Q-corrected.csv", "E_surf corrected (J)"),
               WithinRel(Single(inside_dir / "surface-Q-edge.csv", "E_out (J)") +
                             ModelEnergySum(inside_dir, 0),
                         1.0e-10));
  }
  // The electrostatic driver: the same identity on its fixed-trace / fixed-flux columns for
  // both terminal sources and the record in palace.json; the control library writes no
  // table and no record.
  // The value of the row of `source` in a per-source table (the "i" / "source" key).
  auto SourceValue =
      [&](const fs::path &path, const char *key, int source, const std::string &column)
  {
    const TextCsv table = ReadTrimmedCsv(path);
    for (const auto &row : table.rows)
    {
      if (std::lround(std::stod(row[table.Column(key)])) == source)
      {
        return std::stod(row[table.Column(column)]);
      }
    }
    FAIL("no row of source " << source << " in " << path.string());
    return 0.0;
  };
  auto SourceModelEnergySum = [&](const fs::path &dir, int source, int evaluation)
  {
    const TextCsv models = ReadTrimmedCsv(dir / "surface-response-model-energy.csv");
    double sum = 0.0;
    for (const auto &row : models.rows)
    {
      if (std::lround(std::stod(row[models.Column("source")])) == source &&
          std::lround(std::stod(row[models.Column("evaluation")])) == evaluation)
      {
        sum += std::stod(row[models.Column("fabricated surface energy[4] (J)")]);
      }
    }
    return sum;
  };
  auto SourcePortionTotal = [&](const fs::path &path, int source, int evaluation)
  {
    if (!fs::exists(path))
    {
      return 0.0;
    }
    const TextCsv table = ReadTrimmedCsv(path);
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
      if (std::lround(std::stod(row[table.Column("source")])) == source &&
          std::lround(std::stod(row[table.Column("evaluation")])) == evaluation &&
          row[table.Column("type")] == "Total")
      {
        return std::stod(row[energy_column]);
      }
    }
    FAIL("no Total row of source " << source << " in " << path.string());
    return 0.0;
  };
  {
    REQUIRE(fs::is_regular_file(es_straddling_dir / scp_table));
    for (const int source : {1, 2})
    {
      INFO("source " << source);
      const double scp_ft = SourcePortionTotal(es_straddling_dir / scp_table, source, 0);
      const double scp_ff = SourcePortionTotal(es_straddling_dir / scp_table, source, 1);
      CHECK(scp_ft > 0.0);
      const double outside =
          SourceValue(es_straddling_dir / "surface-Q-edge.csv", "i", source, "E_out (J)");
      // The other kept-raw classes of the identity, if the run wrote them.
      const double other_kept =
          SourcePortionTotal(es_straddling_dir /
                                 "surface-response-domain-boundary-energy.csv",
                             source, 0) +
          SourcePortionTotal(es_straddling_dir / "surface-response-uncovered-energy.csv",
                             source, 0);
      CHECK_THAT(SourceValue(es_straddling_dir / "surface-Q-corrected.csv", "i", source,
                             "E_surf postprocessed fixed-trace[4] (J)"),
                 WithinRel(outside + SourceModelEnergySum(es_straddling_dir, source, 0) +
                               other_kept + scp_ft,
                           1.0e-10));
      CHECK_THAT(SourceValue(es_straddling_dir / "surface-Q-corrected.csv", "i", source,
                             "E_surf postprocessed fixed-flux[4] (J)"),
                 WithinRel(outside + SourceModelEnergySum(es_straddling_dir, source, 1) +
                               other_kept + scp_ff,
                           1.0e-10));
      CHECK_THAT(
          SourceValue(es_fixture_dir / "surface-Q-corrected.csv", "i", source,
                      "E_surf postprocessed fixed-trace[4] (J)"),
          WithinRel(
              SourceValue(es_fixture_dir / "surface-Q-edge.csv", "i", source, "E_out (J)") +
                  SourceModelEnergySum(es_fixture_dir, source, 0) +
                  SourcePortionTotal(es_fixture_dir /
                                         "surface-response-domain-boundary-energy.csv",
                                     source, 0) +
                  SourcePortionTotal(
                      es_fixture_dir / "surface-response-uncovered-energy.csv", source, 0),
              1.0e-10));
    }
    std::ifstream metadata_input(es_straddling_dir / "palace.json");
    REQUIRE(metadata_input);
    const json metadata = json::parse(metadata_input);
    const auto &diagnostics = metadata.at("SurfaceResponse").at("Diagnostics");
    REQUIRE(diagnostics.contains("SpatialClusterPriority"));
    CheckRecord(diagnostics.at("SpatialClusterPriority"), outside_length);
    CHECK(!fs::exists(es_fixture_dir / scp_table));
    std::ifstream fixture_metadata_input(es_fixture_dir / "palace.json");
    REQUIRE(fixture_metadata_input);
    const json fixture_metadata = json::parse(fixture_metadata_input);
    CHECK(!fixture_metadata.at("SurfaceResponse")
               .at("Diagnostics")
               .contains("SpatialClusterPriority"));
    // The raw solves are the same (the libraries differ only in the cluster's support).
    for (const char *file :
         {"terminal-C.csv", "domain-E.csv", "surface-Q.csv", "surface-Q-edge.csv"})
    {
      INFO(file);
      CHECK(ReadFile(es_straddling_dir / file) == ReadFile(es_fixture_dir / file));
      CHECK(ReadFile(es_straddling_dir / file) == ReadFile(es_corner_only_dir / file));
    }
  }

  // (5) The (ii-a) clipped-correction DIAGNOSTIC (OI-SCP-BOXGROW; decisions 958 (b) / 960 /
  //     969): the straddling run writes the companion table with ONE row per (source,
  //     evaluation 0 / 1) for the dropped corner; per interface 4: RawOut = the SCP table's
  //     add-back (the one dropped patch IS the 'convex corner' type), RawWhole > RawOut > 0
  //     (the whole window carries the inside piece too), Fab = the fabricated energy the
  //     SAME corner coupon yields when APPLIED on the same raw field (the corner-only run's
  //     PatchEnergy row: the diagnostic-only patch samples the device trace at the coupon's
  //     own points exactly as the applied patch does), Ratio = Fab / RawWhole and
  //     ClippedCorrection = (Ratio - 1) x RawOut (weight 1); the record in palace.json; no
  //     table on the control run; every existing table of the straddling run untouched by
  //     the diagnostic (asserted through the add-back identity above and the raw-solve
  //     identity).
  {
    const fs::path clipped_table =
        "surface-response-spatial-cluster-priority-clipped-correction.csv";
    REQUIRE(fs::is_regular_file(es_straddling_dir / clipped_table));
    CHECK(!fs::exists(es_fixture_dir / clipped_table));
    CHECK(!fs::exists(es_corner_only_dir / clipped_table));
    // The table carries a '#' comment line naming the parents before its header.
    std::vector<std::string> lines;
    {
      std::ifstream input(es_straddling_dir / clipped_table);
      std::string line;
      while (std::getline(input, line))
      {
        lines.push_back(line);
      }
    }
    REQUIRE(lines.size() == 2 + 4);  // the comment, the header, 2 sources x 2 evaluations
    CHECK_THAT(lines[0], ContainsSubstring("DERIVED"));
    CHECK_THAT(lines[0], ContainsSubstring("parents"));
    const fs::path trimmed = temp.temp_dir / "scp-es-clipped-correction-trimmed.csv";
    {
      std::ofstream output(trimmed);
      for (std::size_t i = 1; i < lines.size(); i++)
      {
        output << lines[i] << '\n';
      }
    }
    const TextCsv clipped = ReadTrimmedCsv(trimmed);
    REQUIRE(clipped.rows.size() == 4);
    std::ifstream metadata_input(es_straddling_dir / "palace.json");
    REQUIRE(metadata_input);
    const json metadata = json::parse(metadata_input);
    const auto &scp_record =
        metadata.at("SurfaceResponse").at("Diagnostics").at("SpatialClusterPriority");
    const int dropped_patch = scp_record.at("Patches").at(0).at("Patch").get<int>();
    REQUIRE(scp_record.contains("ClippedCorrection"));
    CHECK(scp_record.at("ClippedCorrection").at("Patches") == 1);
    CHECK(scp_record.at("ClippedCorrection").at("File") == clipped_table.string());
    CHECK(scp_record.at("ClippedCorrection").contains("Parents"));
    CHECK(metadata.at("SurfaceResponse").at("Interpolation").at("DiagnosticPatches") == 1);
    const TextCsv patch_energy =
        ReadTrimmedCsv(es_corner_only_dir / "surface-response-patch-energy.csv");
    auto AppliedCornerEnergy = [&](int source, int evaluation)
    {
      for (const auto &row : patch_energy.rows)
      {
        if (std::lround(std::stod(row[patch_energy.Column("source")])) == source &&
            std::lround(std::stod(row[patch_energy.Column("evaluation")])) == evaluation &&
            std::abs(std::stod(row[patch_energy.Column("origin x (m)")]) - corner[0]) <
                1.0e-9 &&
            std::abs(std::stod(row[patch_energy.Column("origin y (m)")]) - corner[1]) <
                1.0e-9 &&
            std::abs(std::stod(row[patch_energy.Column("origin z (m)")]) - corner[2]) <
                1.0e-9)
        {
          CHECK_THAT(std::stod(row[patch_energy.Column("weight")]),
                     WithinAbs(1.0, 1.0e-12));
          return std::stod(row[patch_energy.Column("fabricated surface energy[4] (J)")]);
        }
      }
      FAIL("no applied corner patch row of source " << source << " / evaluation "
                                                    << evaluation);
      return 0.0;
    };
    int rows_seen = 0;
    for (const int source : {1, 2})
    {
      for (const int evaluation : {0, 1})
      {
        INFO("source " << source << " evaluation " << evaluation);
        const auto row_it = std::find_if(
            clipped.rows.begin(), clipped.rows.end(),
            [&](const std::vector<std::string> &row)
            {
              return std::lround(std::stod(row[clipped.Column("source")])) == source &&
                     std::lround(std::stod(row[clipped.Column("evaluation")])) ==
                         evaluation;
            });
        REQUIRE(row_it != clipped.rows.end());
        rows_seen++;
        const auto &row = *row_it;
        CHECK(std::lround(std::stod(row[clipped.Column("patch")])) == dropped_patch + 1);
        CHECK_THAT(std::stod(row[clipped.Column("weight")]), WithinAbs(1.0, 1.0e-12));
        const double fab = std::stod(row[clipped.Column("fab energy[4] (J)")]);
        const double whole = std::stod(row[clipped.Column("raw whole energy[4] (J)")]);
        const double outside = std::stod(row[clipped.Column("raw outside energy[4] (J)")]);
        const double ratio = std::stod(row[clipped.Column("ratio[4]")]);
        const double correction =
            std::stod(row[clipped.Column("clipped correction[4] (J)")]);
        CHECK(outside > 0.0);
        CHECK(whole > outside);
        CHECK_THAT(outside, WithinRel(SourcePortionTotal(es_straddling_dir / scp_table,
                                                         source, evaluation),
                                      1.0e-12));
        CHECK_THAT(fab, WithinRel(AppliedCornerEnergy(source, evaluation), 1.0e-10));
        CHECK_THAT(ratio, WithinRel(fab / whole, 1.0e-12));
        CHECK_THAT(correction, WithinRel((ratio - 1.0) * outside, 1.0e-12));
      }
    }
    CHECK(rows_seen == 4);
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
  cluster.SetMeshAlignedBox({0.2137, -0.3, 0.75}, {0.9, 0.3, 1.25});
  cluster.claims = {Claim{3, {0.5, 0.0, 1.0}, {0.9, 0.0, 1.0}}};
  SpatialSupportBounds corner = cluster;
  corner.patch = 5;
  corner.SetMeshAlignedBox({-0.1, -0.3, 0.75}, {0.6, 0.3, 1.25});
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
    middle.SetMeshAlignedBox({0.15, -0.3, 0.75}, {0.3, 0.3, 1.25});
    SpatialSupportBounds whole = cluster;
    whole.patch = 9;
    whole.SetMeshAlignedBox({-0.5, -0.3, 0.75}, {0.6, 0.3, 1.25});
    SpatialSupportBounds touching = cluster;
    touching.patch = 11;
    touching.SetMeshAlignedBox({0.45, -0.3, 0.75}, {0.8, 0.3, 1.25});
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
    middle.SetMeshAlignedBox({0.15, -0.3, 0.75}, {0.3, 0.3, 1.25});
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
