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
// single cross-section of its quadrature point (decision 152: before the fix the
// dimensionless Patch.weight = l / CouponDepth was used as the strip length, so every
// patch was one cross-section). For a potential that is linear along the edge the strip
// average is the value at the cell midpoint, which differs from the value at the Gauss
// point: with the isolated-edge quadrature (two Gauss points per portion, cells = the
// halves) the fabricated surface energy ratio of the longitudinally varying potential to
// the longitudinally constant one must be the weighted mean of the squared potential
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
    // weight perimeter / CouponDepth), one portion per mesh edge (0.125).
    json config = IslandConfig();
    auto &correction = config["Solver"]["Electrostatic"]["ResponseCorrection"];
    correction["TraceCoupling"] = "SurfaceMortar";
    correction["MortarOversampling"] = 2;
    IoData iodata(config, false);
    iodata.boundaries.cracked_attributes.insert(9);
    std::vector<std::unique_ptr<Mesh>> meshes;
    meshes.push_back(std::make_unique<Mesh>(MakeIslandMesh()));
    LaplaceOperator laplace(iodata, meshes);
    SurfaceResponseOperator response(iodata, laplace);
    const double constant_energy = FabricatedResponse(response, laplace, ConstantPotential)
                                       .fabricated_surface_energy.at(4);
    REQUIRE(constant_energy > 0.0);
    const double varying_energy = FabricatedResponse(response, laplace, VaryingPotential)
                                      .fabricated_surface_energy.at(4);
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
#endif
}

}  // namespace palace
