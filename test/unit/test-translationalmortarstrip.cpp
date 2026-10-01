// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "fixtures.hpp"
#include "surfaceresponse-fixtures.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <functional>
#include <map>
#include <memory>
#include <sstream>
#include <string>
#include <tuple>
#include <vector>
#include <mfem.hpp>
#include <catch2/catch_test_macros.hpp>
#include <nlohmann/json.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include "fem/mesh.hpp"
#include "linalg/vector.hpp"
#include "models/laplaceoperator.hpp"
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

// One row of the surface-response patch dry run (surface-response-patches.csv).
struct DryRunPatch
{
  std::string topology;
  double weight = 0.0;
  std::array<double, 3> origin{};
  std::array<double, 3> axis_w{};
  std::array<double, 2> strip{};
};

std::vector<DryRunPatch> ReadDryRunPatches(const fs::path &path)
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
  const std::size_t topology = Column("Topology"), weight = Column("Weight"),
                    origin_x = Column("OriginX"), axis_wx = Column("AxisWX"),
                    strip_begin = Column("StripBegin"), strip_end = Column("StripEnd");
  std::vector<DryRunPatch> patches;
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
    DryRunPatch patch;
    patch.topology = fields[topology];
    patch.weight = std::stod(fields[weight]);
    for (int d = 0; d < 3; d++)
    {
      patch.origin[d] = std::stod(fields[origin_x + d]);
      patch.axis_w[d] = std::stod(fields[axis_wx + d]);
    }
    patch.strip = {std::stod(fields[strip_begin]), std::stod(fields[strip_end])};
    patches.push_back(patch);
  }
  return patches;
}

// The electrostatic response of the island for a projected potential.
SurfaceResponseOperator::ElectrostaticResponse
FabricatedResponse(const SurfaceResponseOperator &response, LaplaceOperator &laplace,
                   const std::function<double(const mfem::Vector &)> &potential)
{
  mfem::ParGridFunction potential_gf(&laplace.GetH1Space().Get());
  mfem::FunctionCoefficient coefficient(potential);
  potential_gf.ProjectCoefficient(coefficient);
  Vector potential_true;
  potential_gf.GetTrueDofs(potential_true);
  return response.GetElectrostaticResponse(potential_true);
}

}  // namespace

// The translational surface mortar (TraceCoupling SurfaceMortar) projects the device trace
// of every translational patch over the patch's longitudinal cell of its edge portion — the
// strip of the portion length that the patch's quadrature weight measures — and not in the
// single cross-section of its quadrature point nor over a strip of another length (decision
// 152: before the fix the dimensionless Patch.weight = l / CouponDepth was used as the
// strip length in the nondimensional mesh coordinates, a strip of l x Lc / CouponDepth
// centred on the Gauss point: five cells here, 3.8 on the transmon, one cross-section where
// that length is below the mortar resolution). For a potential that is linear along the
// edge the strip average is the value at the cell midpoint, which differs from the value at
// the Gauss point: with the isolated-edge quadrature (two Gauss points per portion, cells =
// the halves) the fabricated surface energy ratio of the longitudinally varying potential
// to the longitudinally constant one must be the weighted mean of the squared potential
// factor at the cell midpoints. The Features construction (the device path) is checked
// with the recorded strips of the patch dry run; the legacy 3D construction through the
// patch assignments (its dry run carries no quadrature patches).
TEST_CASE_METHOD(test::SurfaceResponseFiles,
                 "SurfaceResponseOperator translational mortar strip",
                 "[surfaceresponseoperator][3d][mortar][Serial][Parallel]")
{
#if !defined(MFEM_USE_GSLIB)
  SKIP("SurfaceResponseOperator requires MFEM_USE_GSLIB");
#else
  // The PEC island x, z in [0.25, 0.75] on y = 0.5 (MakeIslandMesh): four 0.5-long edges;
  // CouponDepth 0.2, R 0.2 for both libraries.
  constexpr double coupon_depth = 0.2;
  // The potential factor along the perimeter: a function of z alone, so that it is
  // constant across the cross-sections of the two edges x = const (their trace scales
  // with it exactly) and vanishes on the cross-sections (half-width 0.08 < 0.125) of the
  // two edges z = const, which then contribute the same to both potentials' energies; the
  // Q1 hat of the mesh plane z = 0.5 (kinks on the mesh planes 0.375 / 0.5 / 0.625, so the
  // Q1 hexahedra represent the potential exactly), linear on every quadrature portion and
  // cell (the portions of both constructions end at the mesh vertex z = 0.5), so that its
  // strip average is its value at the cell midpoint (the slices sample the cell midpoints).
  // The transverse profile y - 0.5 vanishes on the metal plane (the conductor reference).
  auto Factor = [](double, double z)
  { return std::max(0.0, 1.0 - std::abs(z - 0.5) / 0.125); };
  auto ConstantPotential = [](const mfem::Vector &x) { return x[1] - 0.5; };
  auto VaryingPotential = [&](const mfem::Vector &x)
  { return (x[1] - 0.5) * Factor(x[0], x[2]); };
  const auto &rule = mfem::IntRules.Get(mfem::Geometry::SEGMENT, 2);
  REQUIRE(rule.GetNPoints() == 2);
  const double gauss_offset = 0.5 - rule.IntPoint(0).x;  // 1 / (2 sqrt(3))

  // Group the isolated-edge patches (Gauss points) by island edge, with their longitudinal
  // coordinate; the sorted coordinates of one edge pair up per portion as
  // s_mid -+ L (1/2 - x_0), L the portion length; the cells are the halves of the portion.
  struct EdgePatch
  {
    double s;
    std::size_t index;
  };
  auto GroupByEdge = [&](const auto &origins)
  {
    std::map<std::string, std::vector<EdgePatch>> edges;
    for (std::size_t i = 0; i < origins.size(); i++)
    {
      const auto &origin = origins[i];
      CHECK_THAT(origin[1], WithinAbs(0.5, 1.0e-12));
      const bool on_x_edge =
          std::abs(origin[0] - 0.25) < 1.0e-9 || std::abs(origin[0] - 0.75) < 1.0e-9;
      const bool on_z_edge =
          std::abs(origin[2] - 0.25) < 1.0e-9 || std::abs(origin[2] - 0.75) < 1.0e-9;
      REQUIRE(on_x_edge != on_z_edge);
      edges[on_x_edge ? "x=" + std::to_string(origin[0]) : "z=" + std::to_string(origin[2])]
          .push_back({on_x_edge ? origin[2] : origin[0], i});
    }
    REQUIRE(edges.size() == 4);
    for (auto &[name, edge_patches] : edges)
    {
      std::sort(edge_patches.begin(), edge_patches.end(),
                [](const EdgePatch &a, const EdgePatch &b) { return a.s < b.s; });
      REQUIRE(edge_patches.size() % 2 == 0);
    }
    return edges;
  };
  struct Expectation
  {
    double total_weight = 0.0, strip = 0.0, point = 0.0;
  };
  // Per portion (pair of patches): weights, cell midpoints, the expected strip and point
  // energy factors; `Cell(name, lower, upper, cell_lower, half_length)` checks the strips.
  auto Expect = [&](const std::map<std::string, std::vector<EdgePatch>> &edges,
                    const auto &origins, const auto &weights, const auto &Cell,
                    double discrimination)
  {
    Expectation expectation;
    for (const auto &[name, edge_patches] : edges)
    {
      const bool along_z = name[0] == 'x';
      for (std::size_t k = 0; k < edge_patches.size(); k += 2)
      {
        const std::size_t lower = edge_patches[k].index, upper = edge_patches[k + 1].index;
        const double s_lower = edge_patches[k].s, s_upper = edge_patches[k + 1].s;
        const double length = (s_upper - s_lower) / (2.0 * gauss_offset);
        const double s_mid = 0.5 * (s_lower + s_upper);
        CHECK_THAT(weights[lower], WithinRel(0.5 * length / coupon_depth, 1.0e-9));
        CHECK_THAT(weights[upper], WithinRel(0.5 * length / coupon_depth, 1.0e-9));
        auto FactorAt = [&](std::size_t index, double s)
        { return along_z ? Factor(origins[index][0], s) : Factor(s, origins[index][2]); };
        expectation.total_weight += weights[lower] + weights[upper];
        expectation.strip +=
            weights[lower] * std::pow(FactorAt(lower, s_mid - 0.25 * length), 2) +
            weights[upper] * std::pow(FactorAt(upper, s_mid + 0.25 * length), 2);
        expectation.point += weights[lower] * std::pow(FactorAt(lower, s_lower), 2) +
                             weights[upper] * std::pow(FactorAt(upper, s_upper), 2);
        Cell(along_z, lower, s_mid - 0.5 * length, 0.5 * length);
        Cell(along_z, upper, s_mid, 0.5 * length);
      }
    }
    expectation.strip /= expectation.total_weight;
    expectation.point /= expectation.total_weight;
    // The two semantics differ by a^2 L^2 / 48 per unit weight (slope a, portion length
    // L): on the edges x = const 0.5 % of the mean for the 0.05 portions adjacent to the
    // edge midpoint (Features), 6.7 % for the legacy mesh-edge portions; the check
    // discriminates them.
    REQUIRE(std::abs(expectation.strip - expectation.point) >
            discrimination * expectation.strip);
    return expectation;
  };

  SECTION("Features construction")
  {
    // The convex library (isolated edge + convex-corner-90) keyed on SA only, so that the
    // island's SA-only features match: 4 corner windows of R on each arm, 4 isolated
    // portions of 0.5 - 2 R = 0.1 centred on the edge midpoints, each framed in two
    // portions at the mesh vertex of the midpoint (16 quadrature patches).
    const auto library_path = temp.temp_dir / "fabrication-process-strip-features-3d.json";
    if (Mpi::Root(Mpi::World()))
    {
      std::ifstream input(convex_library_3d_path);
      REQUIRE(input);
      json library = json::parse(input);
      library["Name"] = "unit-test-process-strip-features-3d";
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
    IoData iodata(config, false);
    iodata.boundaries.cracked_attributes.insert(9);

    const auto manifest_path = temp.temp_dir / "surface-response-requirements-strip.json";
    const auto patches_path = temp.temp_dir / "surface-response-patches.csv";
    {
      Mesh dry_run_mesh(MakeIslandMesh());
      WriteSurfaceResponseRequirements(iodata, dry_run_mesh, manifest_path.string());
    }
    Mpi::Barrier(Mpi::World());
    const auto dry_run = ReadDryRunPatches(patches_path);
    std::vector<DryRunPatch> isolated;
    int corner_patches = 0;
    for (const auto &patch : dry_run)
    {
      if (patch.topology == "isolated edge")
      {
        isolated.push_back(patch);
      }
      else
      {
        REQUIRE(patch.topology == "convex corner");
        corner_patches++;
      }
    }
    REQUIRE(corner_patches == 4);
    REQUIRE(isolated.size() == 16);

    std::vector<std::unique_ptr<Mesh>> meshes;
    meshes.push_back(std::make_unique<Mesh>(MakeIslandMesh()));
    LaplaceOperator laplace(iodata, meshes);
    SurfaceResponseOperator response(iodata, laplace);
    REQUIRE(response.GetPatchCount() == static_cast<int>(dry_run.size()));
    auto IsolatedEnergy = [&](const auto &potential)
    {
      const auto result = FabricatedResponse(response, laplace, potential);
      const auto &names = response.GetModelNames();
      double energy = 0.0;
      int isolated_models = 0;
      for (const auto &contribution : result.model_contributions)
      {
        if (names.at(contribution.model) == "isolated")
        {
          energy += contribution.fabricated_surface_energy.at(4);
          isolated_models++;
        }
      }
      REQUIRE(isolated_models == 1);
      return energy;
    };
    const double constant_energy = IsolatedEnergy(ConstantPotential);
    REQUIRE(constant_energy > 0.0);
    const double varying_energy = IsolatedEnergy(VaryingPotential);

    std::vector<std::array<double, 3>> origins;
    std::vector<double> weights;
    for (const auto &patch : isolated)
    {
      origins.push_back(patch.origin);
      weights.push_back(patch.weight);
    }
    const auto edges = GroupByEdge(origins);
    // The recorded strips: the cells of the portion as offsets along AxisW from the Gauss
    // points (AxisW is +-the edge direction), containing their point and tiling it.
    auto Cell = [&](bool along_z, std::size_t index, double cell_lower, double half_length)
    {
      const auto &patch = isolated[index];
      const double direction = along_z ? patch.axis_w[2] : patch.axis_w[0];
      CHECK_THAT(std::abs(direction), WithinAbs(1.0, 1.0e-12));
      CHECK(patch.strip[0] <= 0.0);
      CHECK(patch.strip[1] >= 0.0);
      const double s_origin = along_z ? patch.origin[2] : patch.origin[0];
      const double begin = s_origin + direction * patch.strip[0];
      const double end = s_origin + direction * patch.strip[1];
      CHECK_THAT(std::min(begin, end), WithinAbs(cell_lower, 1.0e-12));
      CHECK_THAT(std::max(begin, end), WithinAbs(cell_lower + half_length, 1.0e-12));
      CHECK_THAT(half_length, WithinAbs(0.025, 1.0e-12));
    };
    const auto expectation = Expect(edges, origins, weights, Cell, 0.004);
    CHECK_THAT(expectation.total_weight, WithinRel(4 * 0.1 / coupon_depth, 1.0e-12));
    CHECK_THAT(varying_energy / constant_energy, WithinRel(expectation.strip, 1.0e-8));
  }

  SECTION("Legacy construction")
  {
    // The concave library under the legacy construction: every edge is isolated-edge
    // patches over the full perimeter (the island corners case asserts the count and the
    // weight perimeter / CouponDepth), one portion per mesh edge (0.125), cells 0.0625. The
    // mortar resolution is the element size (0.125) / MortarOversampling: oversampling 2
    // samples every cell in ONE slice (its midpoint), oversampling 4 in TWO slices at the
    // cell's quarter points, whose uniform 1 / n average is again the midpoint value of the
    // linear potential factor, so both must give the same energy ratio while the mortar
    // point count of the second more than doubles (the contour samples double and every
    // cell adds a slice: with one slice per cell it could at most double).
    long long int single_slice_points = 0;
    for (const int oversampling : {2, 4})
    {
      INFO("MortarOversampling " << oversampling);
      json config = IslandConfig();
      auto &correction = config["Solver"]["Electrostatic"]["ResponseCorrection"];
      correction["TraceCoupling"] = "SurfaceMortar";
      correction["MortarOversampling"] = oversampling;
      IoData iodata(config, false);
      iodata.boundaries.cracked_attributes.insert(9);
      std::vector<std::unique_ptr<Mesh>> meshes;
      meshes.push_back(std::make_unique<Mesh>(MakeIslandMesh()));
      LaplaceOperator laplace(iodata, meshes);
      SurfaceResponseOperator response(iodata, laplace);
      const double constant_energy =
          FabricatedResponse(response, laplace, ConstantPotential)
              .fabricated_surface_energy.at(4);
      REQUIRE(constant_energy > 0.0);
      const double varying_energy = FabricatedResponse(response, laplace, VaryingPotential)
                                        .fabricated_surface_energy.at(4);
      const long long int points =
          response.GetStatistics()["Interpolation"]["PointQueries"]["Total"]
              .get<long long int>();
      if (oversampling == 2)
      {
        single_slice_points = points;
      }
      else
      {
        CHECK(points > 2 * single_slice_points);
        CHECK(points <= 4 * single_slice_points);
      }
      // The patch assignments are gathered on the root rank.
      if (Mpi::Root(Mpi::World()))
      {
        const auto &assignments = response.GetPatchAssignments();
        REQUIRE(static_cast<int>(assignments.size()) == response.GetPatchCount());
        std::vector<std::array<double, 3>> origins;
        std::vector<double> weights;
        for (const auto &patch : assignments)
        {
          origins.push_back(patch.origin);
          weights.push_back(patch.weight);
        }
        const auto edges = GroupByEdge(origins);
        const auto expectation =
            Expect(edges, origins, weights, [](bool, std::size_t, double, double) {}, 0.05);
        CHECK_THAT(expectation.total_weight, WithinRel(2.0 / coupon_depth, 1.0e-12));
        CHECK_THAT(varying_energy / constant_energy, WithinRel(expectation.strip, 1.0e-8));
      }
    }
  }
#endif
}

namespace
{

// A PEC strip x in [0.375, 0.625], z in [0.125, 0.875] on the plane y = 0.5 of the unit
// box (8 x 4 x 8 hexahedra, boundary attribute 9 like MakeIslandMesh), then displaced in x
// by `shear(x, z)` (continuous, small: the hexahedra stay valid) so that its two long sides
// are no longer parallel straight lines while remaining a same-conductor strip of width
// 0.25 < 2R = 0.4 for the identification.
std::unique_ptr<mfem::ParMesh>
MakeShearedStripMesh(const std::function<double(double, double)> &shear)
{
  mfem::Mesh serial =
      mfem::Mesh::MakeCartesian3D(8, 4, 8, mfem::Element::HEXAHEDRON, 1.0, 1.0, 1.0);
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
    const bool inside_strip = xmin >= 0.375 - 1.0e-12 && xmax <= 0.625 + 1.0e-12 &&
                              zmin >= 0.125 - 1.0e-12 && zmax <= 0.875 + 1.0e-12;
    if (on_plane && inside_strip)
    {
      serial.AddBdrElement(serial.GetFace(face)->Duplicate(&serial));
      serial.SetBdrAttribute(serial.GetNBE() - 1, 9);
    }
  }
  serial.FinalizeTopology();
  serial.Finalize();
  serial.Transform(
      [&](const mfem::Vector &x, mfem::Vector &y)
      {
        y = x;
        y(0) += shear(x(0), x(2));
      });
  return std::make_unique<mfem::ParMesh>(Mpi::World(), serial);
}

// One pair patch of the dry run: its portion on the sample's own segment and the frame.
struct PairPatch
{
  int segment = -1;
  double s0 = 0.0, s1 = 0.0, quadrature_weight = 0.0;
  std::array<double, 3> origin{};
  std::array<double, 3> axis_w{};
  std::array<double, 2> strip{};
};

std::vector<PairPatch> ReadPairPatches(const fs::path &path, const std::string &topology)
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
  const std::size_t topology_column = Column("Topology"), segment = Column("Segment"),
                    s0 = Column("S0"), s1 = Column("S1"),
                    quadrature_weight = Column("QuadratureWeight"),
                    origin_x = Column("OriginX"), axis_wx = Column("AxisWX"),
                    strip_begin = Column("StripBegin"), strip_end = Column("StripEnd");
  std::vector<PairPatch> patches;
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
    if (fields[topology_column] != topology)
    {
      continue;
    }
    PairPatch patch;
    patch.segment = std::stoi(fields[segment]);
    patch.s0 = std::stod(fields[s0]);
    patch.s1 = std::stod(fields[s1]);
    patch.quadrature_weight = std::stod(fields[quadrature_weight]);
    for (int d = 0; d < 3; d++)
    {
      patch.origin[d] = std::stod(fields[origin_x + d]);
      patch.axis_w[d] = std::stod(fields[axis_wx + d]);
    }
    patch.strip = {std::stod(fields[strip_begin]), std::stod(fields[strip_end])};
    patches.push_back(patch);
  }
  return patches;
}

}  // namespace

// A pair patch's AxisW follows the partner side (AxisU points at the sample's closest foot
// on the other side), so on the slow tapers and sub-noise polyline bends the identification
// classifies as pairs it is NOT parallel to the sample's own segment: 1 - |cos| = 2.0e-4 on
// every patch of the taper below (each side tilted by atan(0.01) off the strip axis, the
// sides 0.02 rad apart, all 16 feet interior) and 5.4e-4 on one patch of the bend (an 8 deg
// joint of both sides at the mesh plane z = 0.5, noise for the joint rule: implied sagitta
// 0.1875 tan(2 deg) = 0.0066 < 0.05 R = 0.01; the outer side's Gauss point 0.125 / (2 sqrt
// 3) = 0.0264 above the joint lies within 0.25 sin(8 deg) = 0.035 of it, so its foot is
// clamped at the partner's joint vertex and its frame is tilted by atan(0.0264 cos(8 deg) /
// (0.25 - 0.0264 sin(8 deg))) = 6.1 deg off the strip axis, 1.9 deg off its own segment;
// the other 15 feet are interior, |cos| = 1). Before this fix the strip construction
// required |cos| within 1e-8 of 1 and aborted on both (review of decisions 150-152,
// MAJOR-1). The cell is projected onto AxisW by the cosine and oriented by the sign of the
// dot product, so the cells of every portion, mapped back to the arc through that cosine,
// tile the portion.
TEST_CASE_METHOD(test::SurfaceResponseFiles,
                 "SurfaceResponseOperator translational mortar strip of tapered and bent "
                 "pairs",
                 "[surfaceresponseoperator][3d][mortar][Serial][Parallel]")
{
#if !defined(MFEM_USE_GSLIB)
  SKIP("SurfaceResponseOperator requires MFEM_USE_GSLIB");
#else
  // The strip library (isolated edge, concave corner, same-conductor strip) keyed on SA
  // only, its strip coupon at the separation the identification reads for the shape (the
  // library builder's workflow: the requirements manifest of a first dry run names the
  // separation; the signature match admits 1e-3 R).
  const auto library_path = temp.temp_dir / "fabrication-process-strip-pairs-3d.json";
  // One strip model per identified separation (the bend shape reads two exact strips).
  auto WriteLibrary = [&](const std::vector<double> &separations)
  {
    if (Mpi::Root(Mpi::World()))
    {
      std::ifstream input(strip_library_3d_path);
      REQUIRE(input);
      json library = json::parse(input);
      library["Name"] = "unit-test-process-strip-pairs-3d";
      json models = json::array();
      for (auto &model : library["Models"])
      {
        model["Interfaces"] = {{{"Type", "SA"}, {"Coupon", 1}}};
        if (model["Topology"] == "SameConductorStrip")
        {
          for (std::size_t i = 0; i < separations.size(); i++)
          {
            json copy = model;
            copy["Separation"] = separations[i];
            if (i > 0)
            {
              copy["Name"] = model["Name"].get<std::string>() + "-" + std::to_string(i);
            }
            models.push_back(copy);
          }
        }
        else
        {
          models.push_back(model);
        }
      }
      library["Models"] = models;
      std::ofstream output(library_path);
      output << library.dump(2) << "\n";
    }
    Mpi::Barrier(Mpi::World());
  };
  json config = IslandConfig();
  auto &correction = config["Solver"]["Electrostatic"]["ResponseCorrection"];
  correction.erase("PatchConstruction");
  correction["Library"] = library_path.string();
  // The strip's end corners across the 0.25 gap are cluster events without a model.
  correction["UnmatchedPolicy"] = "Warn";
  correction["TraceCoupling"] = "SurfaceMortar";
  correction["MortarOversampling"] = 2;
  const auto &rule = mfem::IntRules.Get(mfem::Geometry::SEGMENT, 2);
  REQUIRE(rule.GetNPoints() == 2);

  struct Shape
  {
    std::string name;
    std::function<double(double, double)> shear;
    // The smallest |cos| a pair patch's frame may have to its own segment.
    double minimum_cosine;
  };
  constexpr double taper_slope = 0.01;
  const double bend_slope = std::tan(8.0 * M_PI / 180.0);
  // The transverse profile of the displacement: +-1 on the strip's sides (|x - 0.5| =
  // 0.125; a uniform stretch about the strip axis inside), falling linearly to 0 at the box
  // walls x = 0 / 1 (the walls stay planar and perpendicular to the metal plane: the
  // identification's reference process normal is the metal plane's).
  auto Profile = [](double x)
  {
    const double u = x - 0.5;
    return std::abs(u) <= 0.125 ? u / 0.125 : std::copysign((0.5 - std::abs(u)) / 0.375, u);
  };
  const std::vector<Shape> shapes = {
      // Each side tilted by atan(taper_slope) away from the strip's axis (the width 0.25 +
      // 2 taper_slope (z - 0.5)): a slow taper, a pair for the identification.
      {"taper", [=](double x, double z) { return taper_slope * (z - 0.5) * Profile(x); },
       std::cos(2.0 * std::atan(taper_slope))},
      // Both sides turn by 8 deg at the mesh plane z = 0.5 (the same displacement on both
      // sides: the width is preserved).
      {"bend", [=](double x, double z)
       { return bend_slope * std::max(0.0, z - 0.5) * std::abs(Profile(x)); },
       std::cos(std::atan(bend_slope))}};
  for (const auto &shape : shapes)
  {
    DYNAMIC_SECTION(shape.name)
    {
      const auto manifest_path =
          temp.temp_dir / ("surface-response-requirements-" + shape.name + ".json");
      const auto patches_path = temp.temp_dir / "surface-response-patches.csv";
      // First dry run: the identified strip(s) and their separations. The taper is one
      // strip at its mean width (the nominal 0.25 shifted by the taper's asymmetric
      // portion, beyond the 1e-3 R match); the bend is TWO exact strips (USER decision 184
      // (3), 2026-10-01: a piece's separation is its own geometry, exact groups within
      // 1e-3 R): the straight half at 0.25 and the sheared half at the perpendicular width
      // 0.25 cos(8 deg) of a pure shear, 1 % apart.
      WriteLibrary({0.25});
      IoData requirements_iodata(config, false);
      requirements_iodata.boundaries.cracked_attributes.insert(9);
      {
        Mesh dry_run_mesh(MakeShearedStripMesh(shape.shear));
        WriteSurfaceResponseRequirements(requirements_iodata, dry_run_mesh,
                                         manifest_path.string());
      }
      Mpi::Barrier(Mpi::World());
      std::vector<double> separations;
      {
        std::ifstream input(manifest_path);
        REQUIRE(input);
        const json manifest = json::parse(input);
        for (const auto &requirement : manifest["Requirements"])
        {
          if (requirement["Topology"] == "SameConductorStrip")
          {
            separations.push_back(requirement["Geometry"]["Separation"].get<double>());
          }
        }
        REQUIRE(separations.size() == (shape.name == "taper" ? 1 : 2));
      }
      for (const double separation : separations)
      {
        CHECK_THAT(separation, WithinAbs(0.25, 0.005));
      }
      WriteLibrary(separations);
      IoData iodata(config, false);
      iodata.boundaries.cracked_attributes.insert(9);
      {
        Mesh dry_run_mesh(MakeShearedStripMesh(shape.shear));
        WriteSurfaceResponseRequirements(iodata, dry_run_mesh, manifest_path.string());
      }
      Mpi::Barrier(Mpi::World());
      const auto pairs = ReadPairPatches(patches_path, "same-conductor strip");
      // Both sides of the pair are patched: at least the portions clear of the end clusters
      // (0.75 - 2 x 0.2 of each 0.75 side, framed at the mesh vertices: >= 3 x 2 Gauss
      // points per side).
      REQUIRE(pairs.size() >= 12);

      // Per portion (segment, s0, s1): its two Gauss patches in position order.
      std::map<std::tuple<int, double, double>, std::vector<std::size_t>> portions;
      for (std::size_t i = 0; i < pairs.size(); i++)
      {
        portions[{pairs[i].segment, pairs[i].s0, pairs[i].s1}].push_back(i);
      }
      int tilted_frames = 0;
      for (const auto &[key, indices] : portions)
      {
        REQUIRE(indices.size() == 2);
        const auto &first = pairs[indices[0]];
        const auto &second = pairs[indices[1]];
        const double length = first.s1 - first.s0;
        REQUIRE(length > 0.0);
        // The sign of Dot(own tangent, AxisW): the origins (gap midpoints) of the two Gauss
        // patches move along the side with the arc parameter.
        std::array<double, 3> displacement{};
        double along_w = 0.0, displacement_norm = 0.0;
        for (int d = 0; d < 3; d++)
        {
          displacement[d] = second.origin[d] - first.origin[d];
          along_w += displacement[d] * first.axis_w[d];
          displacement_norm += displacement[d] * displacement[d];
        }
        displacement_norm = std::sqrt(displacement_norm);
        REQUIRE(displacement_norm > 0.0);
        REQUIRE(std::abs(along_w) > 0.9 * displacement_norm);
        const std::array<double, 2> unit_cells[2] = {{0.0, rule.IntPoint(0).weight},
                                                     {rule.IntPoint(0).weight, 1.0}};
        for (int q = 0; q < 2; q++)
        {
          const auto &patch = pairs[indices[q]];
          CHECK_THAT(patch.quadrature_weight, WithinRel(rule.IntPoint(q).weight, 1.0e-12));
          CHECK(patch.strip[0] <= 0.0);
          CHECK(patch.strip[1] >= 0.0);
          // The projection cosine implied by the cell length; the frame must be within the
          // regime the fix admits and the cell shorter than the arc cell, never longer.
          const double cosine =
              (patch.strip[1] - patch.strip[0]) / (length * patch.quadrature_weight);
          // (S0 / S1 are the manifest's quantized portion parameters: 1e-10 of slack.)
          CHECK(cosine <= 1.0 + 1.0e-9);
          CHECK(cosine >= shape.minimum_cosine - 1.0e-9);
          if (1.0 - cosine > 1.0e-8)
          {
            tilted_frames++;
          }
          // Mapped back to the arc through the cosine and the sign, the cell is the Gauss
          // cell of the portion: [s0, mid] and [mid, s1].
          const double sign = along_w > 0.0 ? 1.0 : -1.0;
          const double t_q = patch.s0 + length * rule.IntPoint(q).x;
          const double begin = t_q + sign * patch.strip[0] / cosine;
          const double end = t_q + sign * patch.strip[1] / cosine;
          CHECK_THAT(std::min(begin, end),
                     WithinAbs(patch.s0 + length * unit_cells[q][0], 1.0e-9));
          CHECK_THAT(std::max(begin, end),
                     WithinAbs(patch.s0 + length * unit_cells[q][1], 1.0e-9));
        }
      }
      // The frames the parallel-to-1e-8 rule aborted on: every patch of the taper, at the
      // cosine of the 0.02 rad between the sides; the one clamped-foot patch of the bend.
      if (shape.name == "taper")
      {
        CHECK(tilted_frames == static_cast<int>(pairs.size()));
        for (const auto &patch : pairs)
        {
          CHECK_THAT((patch.strip[1] - patch.strip[0]) /
                         ((patch.s1 - patch.s0) * patch.quadrature_weight),
                     WithinAbs(shape.minimum_cosine, 1.0e-9));
        }
      }
      else
      {
        CHECK(tilted_frames == 1);
      }

      // The operator itself builds and evaluates on the same mesh (the mortar samples the
      // projected strips).
      std::vector<std::unique_ptr<Mesh>> meshes;
      meshes.push_back(std::make_unique<Mesh>(MakeShearedStripMesh(shape.shear)));
      LaplaceOperator laplace(iodata, meshes);
      SurfaceResponseOperator response(iodata, laplace);
      CHECK(response.GetPatchCount() >= static_cast<int>(pairs.size()));
      const auto result = FabricatedResponse(response, laplace, [](const mfem::Vector &x)
                                             { return x[1] - 0.5; });
      CHECK(std::isfinite(result.fabricated_surface_energy.at(4)));
      CHECK(result.fabricated_surface_energy.at(4) > 0.0);
    }
  }
#endif
}

}  // namespace palace
