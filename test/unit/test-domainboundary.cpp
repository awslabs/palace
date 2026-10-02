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
#include <memory>
#include <set>
#include <sstream>
#include <string>
#include <vector>
#include <mfem.hpp>
#include <catch2/catch_test_macros.hpp>
#include <nlohmann/json.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>
#include "fem/mesh.hpp"
#include "models/laplaceoperator.hpp"
#include "models/surfaceresponseidentification.hpp"
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
std::unique_ptr<mfem::ParMesh> MakeObliqueCutLeadMesh(double shear, bool notch)
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
  return std::make_unique<mfem::ParMesh>(Mpi::World(), serial);
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
// containment test whose decision is independent of the rank count ([Parallel]) and of
// the locator path (PALACE_RESPONSE_USE_GSLIB_POINTS). Here the lead's long edges meet the
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

    // The operator: the same decision (14 applied patches, the same record), with the
    // default locator path and with FindPointsGSLIB forced.
    auto Construct = [&](const char *gslib)
    {
      if (gslib)
      {
        setenv("PALACE_RESPONSE_USE_GSLIB_POINTS", gslib, 1);
      }
      std::vector<std::unique_ptr<Mesh>> meshes;
      meshes.push_back(std::make_unique<Mesh>(MakeObliqueCutLeadMesh(shear, false)));
      LaplaceOperator laplace(iodata, meshes);
      SurfaceResponseOperator response(iodata, laplace);
      if (gslib)
      {
        unsetenv("PALACE_RESPONSE_USE_GSLIB_POINTS");
      }
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
    };
    Construct(nullptr);
    Construct("1");
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

}  // namespace palace
