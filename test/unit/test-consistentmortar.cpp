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
#include <random>
#include <sstream>
#include <string>
#include <vector>
#include <fmt/format.h>
#include <mfem.hpp>
#include <catch2/catch_test_macros.hpp>
#include <nlohmann/json.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>
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

// The straight coupon of the test: a square matching contour of half-size R traversed
// counterclockwise from (-R, -R) with knots every `spacing` (generate_edge_response
// write_bases), the metal band [0, t] on the left face x = -R (edge at the origin, gap
// toward +x) between the band knots (-R, t) and (-R, 0). The generator constrains the band
// and publishes the free knots only.
struct BandContour
{
  std::vector<std::array<double, 3>> knots;  // every knot, contour order
  std::vector<bool> band;                    // the two band knots
  std::vector<int> free_indices;             // knots published by a legacy library
};

BandContour MakeBandContour(double R, double t, double spacing)
{
  BandContour contour;
  const int per_side = static_cast<int>(std::lround(2.0 * R / spacing));
  auto Point = [&](double distance) -> std::array<double, 3>
  {
    const double side = 2.0 * R;
    if (distance < side)
    {
      return {-R + distance, -R, 0.0};
    }
    distance -= side;
    if (distance < side)
    {
      return {R, -R + distance, 0.0};
    }
    distance -= side;
    if (distance < side)
    {
      return {R - distance, R, 0.0};
    }
    distance -= side;
    return {-R, R - distance, 0.0};
  };
  for (int k = 0; k < 4 * per_side; k++)
  {
    const auto point = Point(k * spacing);
    const bool on_band = std::abs(point[0] + R) < 1.0e-12 &&
                         (std::abs(point[1]) < 1.0e-12 || std::abs(point[1] - t) < 1.0e-12);
    contour.knots.push_back(point);
    contour.band.push_back(on_band);
    if (!on_band)
    {
      contour.free_indices.push_back(k);
    }
  }
  return contour;
}

void WriteKnots(const fs::path &path, const BandContour &contour, bool all)
{
  std::ofstream output(path);
  output << "x,y,z\n";
  for (int k = 0; k < static_cast<int>(contour.knots.size()); k++)
  {
    if (all || !contour.band[k])
    {
      output << fmt::format("{:.16e},{:.16e},{:.16e}\n", contour.knots[k][0],
                            contour.knots[k][1], contour.knots[k][2]);
    }
  }
}

// Synthetic SPD response matrices on n knots, zero rows / columns at the `zero` knots:
// Q_ij = scale (i == j ? 2 : 0.2), the thin matrices `thin_factor` times the fabricated.
void WriteMatrices(const fs::path &fabricated, const fs::path &thin,
                   const fs::path &fabricated_surface, const fs::path &thin_surface, int n,
                   const std::vector<bool> &zero, double scale, double thin_factor,
                   double radius_m)
{
  auto Entry = [&](int i, int j)
  { return zero[i] || zero[j] ? 0.0 : scale * (i == j ? 2.0 : 0.2); };
  for (const auto &[path, factor] :
       {std::make_pair(fabricated, 1.0), std::make_pair(thin, thin_factor)})
  {
    std::ofstream output(path);
    output << "basis_i,basis_j,Q_ij (J)\n";
    for (int i = 0; i < n; i++)
    {
      for (int j = i; j < n; j++)
      {
        output << i + 1 << "," << j + 1 << ","
               << fmt::format("{:.16e}", factor * Entry(i, j)) << "\n";
      }
    }
  }
  for (const auto &[path, factor] :
       {std::make_pair(fabricated_surface, 1.0), std::make_pair(thin_surface, thin_factor)})
  {
    std::ofstream output(path);
    output << "interface,edge,R (m),basis_i,basis_j,Q_ij (J),Q_total_ij (J)\n";
    for (int edge = 1; edge <= 2; edge++)
    {
      for (int i = 0; i < n; i++)
      {
        for (int j = i; j < n; j++)
        {
          const double value = 0.5 * factor * Entry(i, j);
          output << "1," << edge << "," << radius_m << "," << i + 1 << "," << j + 1 << ","
                 << fmt::format("{:.16e},{:.16e}", value, value) << "\n";
        }
      }
    }
  }
}

struct DryRunRow
{
  std::size_t patch = 0;
  std::string topology;
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
  const std::size_t patch = Column("Patch"), topology = Column("Topology"),
                    origin = Column("OriginX"), axis_u = Column("AxisUX"),
                    axis_v = Column("AxisVX"), axis_w = Column("AxisWX"),
                    strip = Column("StripBegin");
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
    row.topology = fields[topology];
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

Vector ProjectPotential(LaplaceOperator &laplace,
                        const std::function<double(const mfem::Vector &)> &potential)
{
  mfem::ParGridFunction potential_gf(&laplace.GetH1Space().Get());
  mfem::FunctionCoefficient coefficient(potential);
  potential_gf.ProjectCoefficient(coefficient);
  Vector potential_true;
  potential_gf.GetTrueDofs(potential_true);
  return potential_true;
}

}  // namespace

// Decision 404 D1 (the consistent mortar). The straight coupon generators constrain the
// metal band [0, t] where the fabricated metal meets the matching contour and publish the
// free knots only, so the coupon's hats adjacent to the band ramp to zero at (-R, t) and
// (-R, 0); the runtime's translational surface mortar must project onto that basis — with
// the band vertices inserted at library load from the topology and
// Fabrication.MetalThickness (a legacy library), or taken from ZeroTraceIndices (the
// library contract) — and not onto hats spanning the band (the straight-class MS lift of
// decisions 400 / 403). On the PEC island (x, z in [0.25, 0.75] on y = 0.5) with R = 3 h_y,
// t = h_y and knots every h_y, a device potential that vanishes on the band, is linear
// between the knots across the normal and linear along the edges is a Q1 function of the
// mesh AND lies in the coupon's hat space, so the consistent projection must return its
// value at every knot (at the cell midpoint along the strip: the nodal strip average) to
// rounding, while the legacy mortar (one hat segment across the band) misreads the
// band-adjacent knots. The contract route (ZeroTraceIndices on the full knot list) is the
// same projection; a library whose free knots lie on the band is refused; the record names
// the rule; the transpose is the adjoint (the self-consistent operator stays symmetric) and
// the correction is positive semidefinite for a positive defect.
TEST_CASE_METHOD(
    test::SurfaceResponseFiles, "SurfaceResponseOperator consistent translational mortar",
    "[surfaceresponseoperator][consistentmortar][3d][mortar][Serial][Parallel]")
{
#if !defined(MFEM_USE_GSLIB)
  SKIP("SurfaceResponseOperator requires MFEM_USE_GSLIB");
#else
  constexpr int normal_elements = 16;
  constexpr double h = 1.0 / normal_elements;  // 0.0625
  constexpr double R = 3.0 * h, t = h, spacing = h;
  const BandContour contour = MakeBandContour(R, t, spacing);
  const int n_all = static_cast<int>(contour.knots.size());
  const int n_free = static_cast<int>(contour.free_indices.size());
  REQUIRE(n_all == 24);
  REQUIRE(n_free == 22);
  std::vector<int> band_indices;
  for (int k = 0; k < n_all; k++)
  {
    if (contour.band[k])
    {
      band_indices.push_back(k);
    }
  }
  REQUIRE(band_indices.size() == 2);

  const fs::path dir = temp.temp_dir / "consistent-mortar";
  const fs::path free_points = dir / "basis-points-free.csv";
  const fs::path all_points = dir / "basis-points-all.csv";
  std::map<std::string, fs::path> libraries;
  for (const std::string name : {"consistent", "legacy", "contract", "disagreeing"})
  {
    libraries[name] = dir / ("library-" + name + ".json");
  }
  if (Mpi::Root(Mpi::World()))
  {
    fs::create_directories(dir);
    WriteKnots(free_points, contour, false);
    WriteKnots(all_points, contour, true);
    const std::vector<bool> no_zero(n_free, false);
    WriteMatrices(dir / "free-fabricated.csv", dir / "free-thin.csv",
                  dir / "free-fabricated-surface.csv", dir / "free-thin-surface.csv",
                  n_free, no_zero, 1.0e-12, 0.5, R);
    WriteMatrices(dir / "all-fabricated.csv", dir / "all-thin.csv",
                  dir / "all-fabricated-surface.csv", dir / "all-thin-surface.csv", n_all,
                  contour.band, 1.0e-12, 0.5, R);
    auto Library = [&](const std::string &prefix, const fs::path &points, bool thickness)
    {
      json library = {
          {"Version", 3},
          {"TraceLiftVersion", 2},
          {"Name", "unit-test-consistent-mortar-" + prefix},
          {"MatchingRadius", R},
          {"Fabrication",
           {{"InterfaceLayers", {{"SA", {{"Thickness", 0.002}, {"Permittivity", 4.0}}}}}}},
          {"Models",
           {{{"Name", "isolated"},
             {"Topology", "IsolatedEdge"},
             {"CouponDepth", 0.2},
             {"FabricatedMatrix", (dir / (prefix + "-fabricated.csv")).string()},
             {"ThinMatrix", (dir / (prefix + "-thin.csv")).string()},
             {"FabricatedSurfaceMatrix",
              (dir / (prefix + "-fabricated-surface.csv")).string()},
             {"ThinSurfaceMatrix", (dir / (prefix + "-thin-surface.csv")).string()},
             {"BasisPoints", points.string()},
             {"Interfaces", {{{"Type", "SA"}, {"Coupon", 1}}}}}}}};
      if (thickness)
      {
        library["Fabrication"]["MetalThickness"] = t;
      }
      return library;
    };
    auto Write = [&](const std::string &name, const json &library)
    {
      std::ofstream output(libraries.at(name));
      output << library.dump(2) << "\n";
    };
    Write("consistent", Library("free", free_points, true));
    Write("legacy", Library("free", free_points, false));
    json contract = Library("all", all_points, true);
    contract["Models"][0]["ZeroTraceIndices"] = {band_indices[0] + 1, band_indices[1] + 1};
    Write("contract", contract);
    // Every knot listed, none constrained: two free knots on the band.
    Write("disagreeing", Library("all", all_points, true));
  }
  Mpi::Barrier(Mpi::World());

  auto Config = [&](const std::string &library)
  {
    json config = IslandConfig();
    config["Boundaries"]["Postprocessing"]["Dielectric"][0]["EdgeDistances"] = {R};
    auto &correction = config["Solver"]["Electrostatic"]["ResponseCorrection"];
    correction.erase("PatchConstruction");
    correction["Library"] = libraries.at(library).string();
    correction["TraceCoupling"] = "SurfaceMortar";
    correction["MortarOversampling"] = 2;
    correction["UnmatchedPolicy"] = "Warn";
    return config;
  };
  auto MakeMesh = [&]()
  { return MakeIslandMesh(false, false, false, false, false, false, normal_elements); };

  // The device potential: g(v) h(z) with v = y - 0.5 the height above the metal plane,
  // g = 0 on the band [0, t], 3 (v - t) above it, 2 v below it (kinks on mesh planes), and
  // h linear along z (the strip average of a linear function is its midpoint value).
  auto Potential = [&](const mfem::Vector &x)
  {
    const double v = x[1] - 0.5;
    const double g = v > t ? 3.0 * (v - t) : (v < 0.0 ? 2.0 * v : 0.0);
    return g * (1.0 + 0.5 * (x[2] - 0.5));
  };
  auto ExpectedAtKnot = [&](const DryRunRow &row, const std::array<double, 3> &knot)
  {
    const double s = 0.5 * (row.strip[0] + row.strip[1]);
    mfem::Vector position(3);
    for (int d = 0; d < 3; d++)
    {
      position[d] = row.origin[d] + knot[0] * row.axis_u[d] + knot[1] * row.axis_v[d] +
                    s * row.axis_w[d];
    }
    return Potential(position);
  };

  struct Run
  {
    std::vector<DryRunRow> rows;
    std::vector<SurfaceResponseOperator::PatchTrace> traces;
    SurfaceResponseOperator::ElectrostaticResponse response;
    json diagnostics;
  };
  // The dry run (patch frames and strips) and the operator on the same mesh: the traces of
  // every translational patch and the response of the projected potential.
  auto Evaluate = [&](const std::string &library, auto &&extra)
  {
    const json config = Config(library);
    IoData iodata(config, false);
    iodata.boundaries.cracked_attributes.insert(9);
    const auto manifest_path = dir / ("requirements-" + library + ".json");
    {
      Mesh dry_run_mesh(MakeMesh());
      WriteSurfaceResponseRequirements(iodata, dry_run_mesh, manifest_path.string());
    }
    Mpi::Barrier(Mpi::World());
    Run run;
    run.rows = ReadDryRun(dir / "surface-response-patches.csv");
    std::vector<std::unique_ptr<Mesh>> meshes;
    meshes.push_back(std::make_unique<Mesh>(MakeMesh()));
    LaplaceOperator laplace(iodata, meshes);
    SurfaceResponseOperator response(iodata, laplace);
    REQUIRE(response.GetPatchCount() == static_cast<int>(run.rows.size()));
    const Vector V = ProjectPotential(laplace, Potential);
    run.traces = response.GetPatchTraces(V);
    run.response = response.GetElectrostaticResponse(V);
    run.diagnostics = response.GetStatistics()["Diagnostics"]["ConsistentMortar"];
    extra(response, laplace, V);
    return run;
  };
  auto NoExtra = [](const SurfaceResponseOperator &, LaplaceOperator &, const Vector &) {};

  // Coefficient errors per knot of the free basis against the potential at the knot.
  struct KnotErrors
  {
    double band_adjacent = 0.0, elsewhere = 0.0, scale = 0.0;
    int patches = 0;
  };
  auto Errors = [&](const Run &run, const std::vector<int> &knot_of_coefficient)
  {
    KnotErrors errors;
    std::map<std::size_t, const DryRunRow *> rows;
    for (const auto &row : run.rows)
    {
      rows.emplace(row.patch, &row);
    }
    for (const auto &trace : run.traces)
    {
      const auto row = rows.find(static_cast<std::size_t>(trace.patch));
      REQUIRE(row != rows.end());
      REQUIRE(row->second->topology == "isolated edge");
      REQUIRE(trace.contour_size == static_cast<int>(knot_of_coefficient.size()));
      errors.patches++;
      for (int i = 0; i < trace.contour_size; i++)
      {
        const int k = knot_of_coefficient[i];
        const double expected = ExpectedAtKnot(*row->second, contour.knots[k]);
        const double error = std::abs(trace.coefficients[i] - expected);
        errors.scale = std::max(errors.scale, std::abs(expected));
        // The two free knots next to the band on the left face: (-R, 2 t) and (-R, -t).
        const bool adjacent = std::abs(contour.knots[k][0] + R) < 1.0e-12 &&
                              (std::abs(contour.knots[k][1] - 2.0 * t) < 1.0e-12 ||
                               std::abs(contour.knots[k][1] + t) < 1.0e-12);
        (adjacent ? errors.band_adjacent : errors.elsewhere) =
            std::max(adjacent ? errors.band_adjacent : errors.elsewhere, error);
      }
    }
    return errors;
  };
  std::vector<int> all_knots(n_all);
  for (int k = 0; k < n_all; k++)
  {
    all_knots[k] = k;
  }

  SECTION("the consistent mortar reads the nodal strip average; the legacy mortar does not")
  {
    const Run consistent = Evaluate("consistent", NoExtra);
    REQUIRE(consistent.traces.size() >= 8);
    const auto consistent_errors = Errors(consistent, contour.free_indices);
    REQUIRE(consistent_errors.scale > 0.1);
    CHECK(consistent_errors.band_adjacent <= 1.0e-9 * consistent_errors.scale);
    CHECK(consistent_errors.elsewhere <= 1.0e-9 * consistent_errors.scale);
    CHECK(consistent.diagnostics["TranslationalModels"].get<int>() == 1);
    CHECK(consistent.diagnostics["WithInsertedBandVertices"].get<int>() == 1);
    CHECK(consistent.diagnostics["WithZeroTraceIndices"].get<int>() == 0);
    REQUIRE(consistent.diagnostics["Models"].size() == 1);
    const auto &record = consistent.diagnostics["Models"][0];
    CHECK(record["Source"] == "RuntimeRule");
    CHECK_THAT(record["Rule"].get<std::string>(), ContainsSubstring("inserted"));
    REQUIRE(record["BandVertices"].size() == 2);
    // Contour order on the left face is downward: (-R, t) then (-R, 0).
    CHECK_THAT(record["BandVertices"][0][0].get<double>(), WithinAbs(-R, 1.0e-14));
    CHECK_THAT(record["BandVertices"][0][1].get<double>(), WithinAbs(t, 1.0e-14));
    CHECK_THAT(record["BandVertices"][1][0].get<double>(), WithinAbs(-R, 1.0e-14));
    CHECK_THAT(record["BandVertices"][1][1].get<double>(), WithinAbs(0.0, 1.0e-14));

    const Run legacy = Evaluate("legacy", NoExtra);
    REQUIRE(legacy.traces.size() == consistent.traces.size());
    const auto legacy_errors = Errors(legacy, contour.free_indices);
    // One hat segment from (-R, 2 t) to (-R, -t) across the band: the projection of the
    // trace (0 on [0, t], linear beside it) misreads the band-adjacent knots by a finite
    // fraction of the trace amplitude.
    CHECK(legacy_errors.band_adjacent > 1.0e-2 * legacy_errors.scale);
    CHECK(legacy.diagnostics["WithInsertedBandVertices"].get<int>() == 0);
    CHECK(legacy.diagnostics["WithoutBand"].get<int>() == 1);
    CHECK(legacy.diagnostics["Models"][0]["Source"] == "None");
    CHECK_THAT(legacy.diagnostics["Models"][0]["Rule"].get<std::string>(),
               ContainsSubstring("MetalThickness"));
    // The two libraries have the same free basis and matrices: only the projection differs.
    CHECK(legacy.response.fabricated_surface_energy.at(4) !=
          consistent.response.fabricated_surface_energy.at(4));
  }

  SECTION("the library contract (ZeroTraceIndices) is the same projection")
  {
    const Run consistent = Evaluate("consistent", NoExtra);
    const Run contract = Evaluate("contract", NoExtra);
    REQUIRE(contract.traces.size() == consistent.traces.size());
    const auto contract_errors = Errors(contract, all_knots);
    CHECK(contract_errors.band_adjacent <= 1.0e-9 * contract_errors.scale);
    CHECK(contract_errors.elsewhere <= 1.0e-9 * contract_errors.scale);
    for (std::size_t p = 0; p < contract.traces.size(); p++)
    {
      const auto &a = contract.traces[p];
      const auto &b = consistent.traces[p];
      REQUIRE(a.patch == b.patch);
      REQUIRE(a.contour_size == n_all);
      REQUIRE(b.contour_size == n_free);
      for (int i = 0, j = 0; i < n_all; i++)
      {
        if (contour.band[i])
        {
          CHECK(a.coefficients[i] == 0.0);
          continue;
        }
        CHECK_THAT(a.coefficients[i], WithinAbs(b.coefficients[j++], 1.0e-11));
      }
    }
    CHECK_THAT(contract.response.fabricated_surface_energy.at(4),
               WithinRel(consistent.response.fabricated_surface_energy.at(4), 1.0e-10));
    CHECK_THAT(contract.response.domain_correction,
               WithinRel(consistent.response.domain_correction, 1.0e-10));
    CHECK(contract.diagnostics["WithZeroTraceIndices"].get<int>() == 1);
    CHECK(contract.diagnostics["Models"][0]["Source"] == "ZeroTraceIndices");
    CHECK(contract.diagnostics["Models"][0]["BandVertices"].empty());
  }

  SECTION("a library whose free knots lie on the band is refused")
  {
    const json config = Config("disagreeing");
    IoData iodata(config, false);
    iodata.boundaries.cracked_attributes.insert(9);
    Mesh dry_run_mesh(MakeMesh());
    CHECK_THROWS_WITH(
        WriteSurfaceResponseRequirements(iodata, dry_run_mesh,
                                         (dir / "requirements-disagreeing.json").string()),
        ContainsSubstring("constrained metal band"));
  }

  SECTION("the transpose is the adjoint and the correction of a positive defect is PSD")
  {
    std::mt19937 generator(404);
    std::uniform_real_distribution<double> uniform(-1.0, 1.0);
    auto Adjoint = [&](const SurfaceResponseOperator &response, LaplaceOperator &laplace,
                       const Vector &V)
    {
      const int n = laplace.GetH1Space().GetTrueVSize();
      Vector x1(n), x2(n), y1(n), y2(n);
      for (int i = 0; i < n; i++)
      {
        x1(i) = uniform(generator);
        x2(i) = uniform(generator);
      }
      response.FixedTraceDomainDefectMult(x1, y1);
      response.FixedTraceDomainDefectMult(x2, y2);
      const double a12 = linalg::Dot(laplace.GetComm(), x1, y2);
      const double a21 = linalg::Dot(laplace.GetComm(), x2, y1);
      const double a11 = linalg::Dot(laplace.GetComm(), x1, y1);
      const double a22 = linalg::Dot(laplace.GetComm(), x2, y2);
      REQUIRE(std::abs(a11) > 0.0);
      // Symmetry of Pᵀ W D P (the transpose of the consistent mortar is its adjoint).
      CHECK_THAT(a12, WithinRel(a21, 1.0e-10));
      // Q_fab = 2 Q_thin here: D is positive definite on the trace space, so the form is
      // positive and equals twice the fixed-trace domain correction.
      CHECK(a11 > 0.0);
      CHECK(a22 > 0.0);
      Vector yV(n);
      response.FixedTraceDomainDefectMult(V, yV);
      CHECK_THAT(
          0.5 * linalg::Dot(laplace.GetComm(), V, yV),
          WithinRel(response.GetElectrostaticResponse(V).domain_correction, 1.0e-10));
    };
    Evaluate("consistent", Adjoint);
    Evaluate("contract", Adjoint);
  }
#endif
}

}  // namespace palace
