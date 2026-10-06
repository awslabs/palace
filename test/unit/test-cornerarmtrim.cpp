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
#include "fem/mesh.hpp"
#include "models/laplaceoperator.hpp"
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

namespace
{

// A PEC lead (attribute 9) on the plane y = 0.5 of the unit cube meshed 8 x 8 x 8 (h =
// 0.125): x in [0.25, x_max] from the bottom face z = 0 (where the lead is cut by the
// domain) to its cap at z = 0.5. With x_max = 0.75 the cap ends in two 90-degree convex
// corners; with x_max = 1 the lead reaches the wall x = 1 (a truncation cut) and has one
// corner at (0.25, 0.5, 0.5). The mesh is then sheared IN THE METAL PLANE, x -> x + shear
// (z - 0.5): the cap (along x) is unchanged while the long edges (along z) tilt to the
// direction (shear, 0, 1), so the corner angle between the arms becomes 90 + atan(shear)
// degrees (shear = 1 / sqrt(3): 120 degrees) and the lead stays planar on y = 0.5.
mfem::Mesh MakeSerialLeadMesh(double shear, double x_max)
{
  constexpr int n = 8;
  constexpr double h = 1.0 / n;
  auto Vertex = [](int i, int j, int k) { return i + (n + 1) * (j + (n + 1) * k); };
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
  // The lead: the interior faces on y = 0.5 with x in [0.25, x_max], z in [0, 0.5].
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
    if (on_plane && xmin >= 0.25 - 1.0e-12 && xmax <= x_max + 1.0e-12 &&
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
    point[0] += shear * (point[2] - 0.5);
  }
  return serial;
}

std::unique_ptr<mfem::ParMesh> MakeLeadMesh(double shear, double x_max)
{
  mfem::Mesh serial = MakeSerialLeadMesh(shear, x_max);
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

std::string ReadFile(const fs::path &path)
{
  std::ifstream input(path, std::ios::binary);
  REQUIRE(input);
  return std::string(std::istreambuf_iterator<char>(input), {});
}

Table LoadCsv(const fs::path &path)
{
  TableWithCSVFile wrapped(path.string(), /*load_existing_file=*/true);
  return std::move(wrapped.table);
}

const Column &ColumnByHeader(const Table &table, const std::string &header)
{
  for (auto it = table.cbegin(); it != table.cend(); ++it)
  {
    if (it->header_text == header)
    {
      return *it;
    }
  }
  FAIL("No column \"" << header << "\"");
  return *table.cbegin();
}

// The plain text table of surface-response-uncovered-energy.csv: header fields and rows.
struct TextCsv
{
  std::vector<std::string> header;
  std::vector<std::vector<std::string>> rows;
  std::size_t Column(const std::string &name) const
  {
    const auto it = std::find(header.begin(), header.end(), name);
    REQUIRE(it != header.end());
    return static_cast<std::size_t>(it - header.begin());
  }
};

TextCsv ReadTextCsv(const fs::path &path)
{
  std::ifstream input(path);
  REQUIRE(input);
  TextCsv csv;
  std::string line;
  std::getline(input, line);
  std::stringstream header(line);
  std::string field;
  while (std::getline(header, field, ','))
  {
    csv.header.push_back(field);
  }
  while (std::getline(input, line))
  {
    if (line.empty())
    {
      continue;
    }
    std::vector<std::string> fields;
    std::stringstream stream(line);
    while (std::getline(stream, field, ','))
    {
      fields.push_back(field);
    }
    REQUIRE(fields.size() == csv.header.size());
    csv.rows.push_back(fields);
  }
  return csv;
}

void RunElectrostatic(json config)
{
  IoData iodata(std::move(config), /*print=*/false);
  MPI_Comm comm = Mpi::World();
  MakeOutputFolder(iodata, comm);
  const int omp_threads = utils::ConfigureOmp();
  BlockTimer::Reset();
  palace::Run(iodata, comm, omp_threads, /*git_tag=*/nullptr);
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
  }
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

}  // namespace palace
