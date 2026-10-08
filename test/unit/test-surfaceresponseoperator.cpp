// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "fixtures.hpp"
#include "surfaceresponse-fixtures.hpp"

#include <array>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <limits>
#include <map>
#include <memory>
#include <optional>
#include <set>
#include <sstream>
#include <string_view>
#include <tuple>
#include <vector>
#include <mfem.hpp>
#include <catch2/catch_test_macros.hpp>
#include <nlohmann/json.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>
#include "fem/gridfunction.hpp"
#include "fem/mesh.hpp"
#include "linalg/vector.hpp"
#include "models/boundarymodeoperator.hpp"
#include "models/cornertracebasis.hpp"
#include "models/laplaceoperator.hpp"
#include "models/spaceoperator.hpp"
#include "models/surfaceresponseidentification.hpp"
#include "models/surfaceresponseoperator.hpp"
#include "utils/communication.hpp"
#include "utils/edgedistance.hpp"
#include "utils/iodata.hpp"
#include "utils/metaledge.hpp"

namespace palace
{

namespace fs = std::filesystem;

using json = nlohmann::json;
using namespace Catch::Matchers;

namespace
{

// Corner-coupon basis files of a unit-test corner family. The rule layout (the corner
// family's trace basis rule, MakeCornerBoxSeed + BuildCornerTraceBasis; with a segment
// connectivity angle in degrees the band merge is keyed by the rule's layout at that angle
// and the TraceBasis record carries ConnectivityAngleDegrees, corner-qualification block
// 2026-09-29; without one a legacy node, exact matches only) or the lane-2
// angle-independent 8-knot layout (every ring the fixed knots; the zero set = the knots on
// the metal footprint, the corner-family review's defective basis) at one angle; written by
// the calling rank. Returns {basis-points, trace-vertices, trace-triangles, zero set
// (1-based), contour groups}.
struct CornerBasisFiles
{
  fs::path points, vertices, triangles;
  std::vector<int> zero_trace_indices;
  std::vector<int> contour_groups;
  json trace_basis;
};

CornerBasisFiles
WriteCornerBasisFiles(const fs::path &directory, const std::string &tag,
                      double angle_degrees, bool convex, double radius,
                      double metal_thickness, double overetch_depth, bool rule_layout,
                      std::optional<double> connectivity_degrees = std::nullopt)
{
  const CornerTraceBasisRule rule;
  const auto seed =
      MakeCornerBoxSeed(radius, metal_thickness, overetch_depth, convex, rule);
  const double angle = angle_degrees * M_PI / 180.0;
  MFEM_VERIFY(rule_layout || !connectivity_degrees,
              "A segment connectivity angle needs the rule layout!");
  std::optional<double> connectivity;
  if (connectivity_degrees)
  {
    connectivity = *connectivity_degrees * M_PI / 180.0;
  }
  CornerBasisFiles files;
  files.points = directory / ("corner-basis-" + tag + "-points.csv");
  files.vertices = directory / ("corner-basis-" + tag + "-vertices.csv");
  files.triangles = directory / ("corner-basis-" + tag + "-triangles.csv");
  files.contour_groups = seed.contour_groups;
  std::vector<std::array<double, 3>> knots;
  std::vector<int> zero;
  std::vector<ConstructedCornerTraceBasis::Vertex> vertices;
  std::vector<std::array<int, 3>> triangles;
  if (rule_layout)
  {
    const auto basis =
        BuildCornerTraceBasis(seed.points, seed.contour_groups, seed.zero_trace_indices,
                              angle, convex, rule, connectivity);
    knots = basis.knots;
    for (int k = 0; k < static_cast<int>(knots.size()); k++)
    {
      if (basis.zero[k])
      {
        zero.push_back(k);
      }
    }
    vertices = basis.vertices;
    triangles = basis.triangles;
    files.trace_basis = {{"RingSize", rule.ring_size},
                         {"MetalInteriorKnots", rule.metal_interior_knots},
                         {"FreeKnots", rule.free_knots},
                         {"Fractions", rule.fractions}};
    if (connectivity_degrees)
    {
      files.trace_basis["ConnectivityAngleDegrees"] = *connectivity_degrees;
    }
  }
  else
  {
    // The lane-2 layout: the seed's fixed rings at every level; PEC = the knots on the
    // closed metal footprint (the sector 0 .. angle for a convex corner) at z = 0 and
    // z = metal_thickness.
    knots = seed.points;
    int offset = 0;
    for (const int size : seed.contour_groups)
    {
      const double z = knots[offset][2];
      if (std::abs(z) <= 1.0e-12 || std::abs(z - metal_thickness) <= 1.0e-12)
      {
        for (int i = 0; i < size; i++)
        {
          const auto &point = knots[offset + i];
          const double first = point[1];
          const double second = point[0] * std::sin(angle) - point[1] * std::cos(angle);
          const bool in_wedge = first >= -1.0e-12 && second >= -1.0e-12;
          if (convex ? in_wedge : !in_wedge)
          {
            zero.push_back(offset + i);
          }
        }
      }
      offset += size;
    }
    for (int k = 0; k < static_cast<int>(knots.size()); k++)
    {
      ConstructedCornerTraceBasis::Vertex vertex;
      vertex.point = knots[k];
      vertex.basis = k;
      vertices.push_back(vertex);
    }
    // The lane-2 connect_rings / cap fans (equal rings).
    const int ring_size = seed.contour_groups.front();
    const int ring_count = static_cast<int>(seed.contour_groups.size());
    const int outer = ring_count - 2;
    auto Connect = [&](int first, int second)
    {
      for (int i = 0; i < ring_size; i++)
      {
        const int next = (i + 1) % ring_size;
        triangles.push_back({first + i, first + next, second + next});
        triangles.push_back({first + i, second + next, second + i});
      }
    };
    for (int r = 0; r + 1 < outer; r++)
    {
      Connect(r * ring_size, (r + 1) * ring_size);
    }
    Connect((outer - 1) * ring_size, outer * ring_size);
    Connect((outer + 1) * ring_size, 0);
    for (int i = 1; i + 1 < ring_size; i++)
    {
      triangles.push_back(
          {outer * ring_size, outer * ring_size + i, outer * ring_size + i + 1});
      triangles.push_back({(outer + 1) * ring_size + i + 1, (outer + 1) * ring_size + i,
                           (outer + 1) * ring_size});
    }
  }
  const std::set<int> zero_set(zero.begin(), zero.end());
  for (const int index : zero)
  {
    files.zero_trace_indices.push_back(index + 1);
  }
  {
    std::ofstream output(files.points);
    output << std::setprecision(17) << "x,y,z\n";
    for (const auto &knot : knots)
    {
      output << knot[0] << "," << knot[1] << "," << knot[2] << "\n";
    }
  }
  {
    std::ofstream output(files.vertices);
    output << std::setprecision(17)
           << "vertex,x,y,z,basis,conductor,parent_a,parent_b,weight_a\n";
    for (std::size_t v = 0; v < vertices.size(); v++)
    {
      const auto &vertex = vertices[v];
      output << v + 1 << "," << vertex.point[0] << "," << vertex.point[1] << ","
             << vertex.point[2] << ",";
      if (vertex.basis >= 0)
      {
        output << vertex.basis + 1 << "," << (zero_set.count(vertex.basis) ? 1 : 0)
               << ",0,0,0\n";
      }
      else
      {
        output << "0,0," << vertex.parent_a + 1 << "," << vertex.parent_b + 1 << ","
               << vertex.weight_a << "\n";
      }
    }
  }
  {
    std::ofstream output(files.triangles);
    output << "triangle,vertex_i,vertex_j,vertex_k\n";
    for (std::size_t t = 0; t < triangles.size(); t++)
    {
      output << t + 1 << "," << triangles[t][0] + 1 << "," << triangles[t][1] + 1 << ","
             << triangles[t][2] + 1 << "\n";
    }
  }
  return files;
}

// Synthetic N x N response matrices (domain and one-interface surface, with the within-R
// column) of a unit-test corner coupon: diagonal `diagonal`, off-diagonal couplings
// decaying with the index distance; `bump` is added to the diagonal entry of basis index
// `bump_index` (0-based; none when negative).
void WriteCornerMatrices(const fs::path &domain_path, const fs::path &surface_path,
                         int size, double diagonal, double coupling_scale, double radius_m,
                         int bump_index = -1, double bump = 0.0)
{
  std::ofstream domain(domain_path);
  domain << "basis_i,basis_j,Q_ij (J)\n";
  std::ofstream surface(surface_path);
  surface << "interface,edge,R (m),basis_i,basis_j,Q_ij (J),Q_total_ij (J)\n";
  for (int i = 0; i < size; i++)
  {
    for (int j = 0; j < size; j++)
    {
      const double coupling = 1.0 / (1.0 + std::abs(i - j));
      const double value =
          (i == j ? diagonal + (i == bump_index ? bump : 0.0) : coupling_scale * coupling) *
          1.0e-12;
      if (j >= i)
      {
        domain << i + 1 << "," << j + 1 << "," << value << "\n";
        surface << "1,1," << radius_m << "," << i + 1 << "," << j + 1 << "," << value << ","
                << value << "\n";
      }
    }
  }
}

// A hexahedral box mesh whose metal island (cracked boundary attribute 9 on the plane
// y = 0.5) is an exact star-shaped polygon: the radial map of every concentric square of
// the (x, z) grid onto the polygon (its vertices must lie on grid rays so the outline is
// the exact polygon at every level), blended to the identity between the levels 2 and 3 so
// the PEC box walls stay on the bounding box.
std::unique_ptr<mfem::ParMesh>
MakePolygonIslandMesh(const std::vector<std::array<double, 2>> &polygon, double extent,
                      double h)
{
  const int n = static_cast<int>(std::lround(extent / h));
  mfem::Mesh serial =
      mfem::Mesh::MakeCartesian3D(n, 4, n, mfem::Element::HEXAHEDRON, extent, 1.0, extent);
  const double c = 0.5 * extent;
  auto Level = [&](const double *point)
  { return std::max(std::abs(point[0] - c), std::abs(point[2] - c)); };
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
    double level_max = 0.0;
    for (const int vertex : vertices)
    {
      const double *point = serial.GetVertex(vertex);
      on_plane = on_plane && std::abs(point[1] - 0.5) < 1.0e-12;
      level_max = std::max(level_max, Level(point));
    }
    if (on_plane && level_max <= 1.0 + 1.0e-9)
    {
      serial.AddBdrElement(serial.GetFace(face)->Duplicate(&serial));
      serial.SetBdrAttribute(serial.GetNBE() - 1, 9);
    }
  }
  auto PolygonRadius = [&](double ux, double uz)
  {
    double r = mfem::infinity();
    for (std::size_t k = 0; k < polygon.size(); k++)
    {
      const auto &a = polygon[k], &b = polygon[(k + 1) % polygon.size()];
      const double nx = b[1] - a[1], nz = a[0] - b[0];  // outward normal of side a -> b
      const double denominator = nx * ux + nz * uz;
      if (denominator > 1.0e-14)
      {
        r = std::min(r, (nx * a[0] + nz * a[1]) / denominator);
      }
    }
    return r;
  };
  for (int vertex = 0; vertex < serial.GetNV(); vertex++)
  {
    double *point = serial.GetVertex(vertex);
    const double lx = point[0] - c, lz = point[2] - c;
    const double norm = std::hypot(lx, lz);
    const double level = std::max(std::abs(lx), std::abs(lz));
    if (norm > 1.0e-12 && level < 3.0 - 1.0e-9)
    {
      const double ux = lx / norm, uz = lz / norm;
      const double square = 1.0 / std::max(std::abs(ux), std::abs(uz));
      const double polygon_scale = PolygonRadius(ux, uz) / square;
      const double blend = level <= 2.0 ? 1.0 : (3.0 - level);
      const double scale = 1.0 + blend * (polygon_scale - 1.0);
      point[0] = c + scale * lx;
      point[2] = c + scale * lz;
    }
  }
  serial.FinalizeTopology();
  serial.Finalize();
  return std::make_unique<mfem::ParMesh>(Mpi::World(), serial);
}

}  // namespace

// The SurfaceResponseOperator cases (TEST_CASE_METHOD on test::SurfaceResponseFiles; every
// case writes the fixture files it reads and builds its own meshes and operators). The 2D
// cases run on 4 x 4 .. 10 x 4 triangle meshes, the 3D CPW cases on the cpw3d-surface-nc
// test mesh (one mesh read per boundary configuration), the island cases on synthetic
// hexahedral boxes. The quadratic (high-order) spatial cluster round trips are [Long].

TEST_CASE_METHOD(test::SurfaceResponseFiles, "SurfaceResponseOperator explicit patches 2D",
                 "[surfaceresponseoperator][2d][patch][Serial][Parallel]")
{
#if !defined(MFEM_USE_GSLIB)
  SKIP("SurfaceResponseOperator requires MFEM_USE_GSLIB");
#else
  json config = {
      {"Problem", {{"Type", "Electrostatic"}, {"Output", temp.temp_dir.string()}}},
      {"Model", {{"Mesh", "unused.msh"}}},
      {"Domains", {{"Materials", {{{"Attributes", {1}}}}}}},
      {"Boundaries",
       {{"Ground", {{"Attributes", {1, 3, 4}}}},
        {"Terminal", {{{"Index", 1}, {"Attributes", {2}}}}}}},
      {"Solver",
       {{"Order", 1},
        {"Electrostatic",
         {{"ResponseCorrection",
           {{"Models",
             {{{"Index", 1},
               {"FabricatedMatrix", fabricated_path.string()},
               {"ThinMatrix", thin_path.string()},
               {"FabricatedSurfaceMatrix", fabricated_surface_path.string()},
               {"ThinSurfaceMatrix", thin_surface_path.string()},
               {"BasisPoints", points_path.string()},
               {"Interfaces", {{{"Target", 4}, {"Coupon", 1}}}}}}},
            {"Patches",
             {{{"Model", 1},
               {"Origin", {0.87, 0.35, 0.0}},
               {"AxisU", {1.0, 0.0, 0.0}},
               {"AxisV", {0.0, 1.0, 0.0}},
               {"Reference", {0.0, 0.0, 0.0}}}}}}}}}}}};
  IoData iodata(config, false);

  mfem::Mesh serial_mesh =
      mfem::Mesh::MakeCartesian2D(4, 4, mfem::Element::TRIANGLE, false, 1.0, 1.0);
  while (serial_mesh.GetNE() < Mpi::Size(Mpi::World()))
  {
    serial_mesh.UniformRefinement();
  }

  auto parallel_mesh = std::make_unique<mfem::ParMesh>(Mpi::World(), serial_mesh);
  std::vector<std::unique_ptr<Mesh>> meshes;
  meshes.push_back(std::make_unique<Mesh>(std::move(parallel_mesh)));

  LaplaceOperator laplace_op(iodata, meshes);
  auto prescribed_config = config;
  prescribed_config["Boundaries"] = {{"Ground", {{"Attributes", {1, 4}}}},
                                     {"PrescribedPotential",
                                      {{{"Index", 1},
                                        {"Attributes", {2}},
                                        {"TerminalAttributes", {3}},
                                        {"DataFile", zero_trace_path.string()}}}}};
  prescribed_config["Solver"]["Electrostatic"] = json::object();
  IoData prescribed_iodata(prescribed_config, false);
  LaplaceOperator prescribed_laplace(prescribed_iodata, meshes);
  auto prescribed_stiffness = prescribed_laplace.GetStiffnessMatrix();
  Vector prescribed_excitation, prescribed_rhs;
  prescribed_laplace.GetExcitationVector(1, *prescribed_stiffness, prescribed_excitation,
                                         prescribed_rhs);
  double prescribed_max = prescribed_excitation.Normlinf();
  Mpi::GlobalMax(1, &prescribed_max, Mpi::World());
  CHECK_THAT(
      prescribed_max,
      WithinRel(1.0 / prescribed_iodata.units.GetScaleFactor<Units::ValueType::VOLTAGE>(),
                1.0e-12));

  // A trace surface and conductor surface can share boundary dofs. The lift for one
  // combined prescribed-potential source must equal the sum of independently generated
  // trace and conductor-state lifts; response matrices rely on this exact superposition.
  auto split_source_config = prescribed_config;
  split_source_config["Boundaries"] = {
      {"Ground", {{"Attributes", {1, 4}}}},
      {"PrescribedPotential",
       {{{"Index", 1},
         {"Attributes", {2}},
         {"DataFile", shared_boundary_trace_path.string()}},
        {{"Index", 2},
         {"Attributes", {2}},
         {"TerminalAttributes", {3}},
         {"DataFile", zero_trace_path.string()}}}}};
  IoData split_source_iodata(split_source_config, false);
  LaplaceOperator split_source_laplace(split_source_iodata, meshes);
  auto split_source_stiffness = split_source_laplace.GetStiffnessMatrix();
  Vector trace_excitation, trace_rhs, conductor_excitation, conductor_rhs;
  split_source_laplace.GetExcitationVector(1, *split_source_stiffness, trace_excitation,
                                           trace_rhs);
  split_source_laplace.GetExcitationVector(2, *split_source_stiffness, conductor_excitation,
                                           conductor_rhs);

  auto combined_source_config = prescribed_config;
  combined_source_config["Boundaries"] = {
      {"Ground", {{"Attributes", {1, 4}}}},
      {"PrescribedPotential",
       {{{"Index", 1},
         {"Attributes", {2}},
         {"TerminalAttributes", {3}},
         {"DataFile", shared_boundary_trace_path.string()}}}}};
  IoData combined_source_iodata(combined_source_config, false);
  LaplaceOperator combined_source_laplace(combined_source_iodata, meshes);
  auto combined_source_stiffness = combined_source_laplace.GetStiffnessMatrix();
  Vector combined_excitation, combined_rhs;
  combined_source_laplace.GetExcitationVector(1, *combined_source_stiffness,
                                              combined_excitation, combined_rhs);
  trace_excitation += conductor_excitation;
  trace_excitation -= combined_excitation;
  trace_rhs += conductor_rhs;
  trace_rhs -= combined_rhs;
  CHECK(linalg::Norml2(Mpi::World(), trace_excitation) == 0.0);
  CHECK(linalg::Norml2(Mpi::World(), trace_rhs) == 0.0);

  SurfaceResponseOperator response(iodata, laplace_op);
  REQUIRE(response.GetBasisSize() == 4);
  REQUIRE(response.GetPatchCount() == 1);
  REQUIRE(response.HasSurfaceResponse());
  auto compact_config = config;
  auto &compact_model =
      compact_config["Solver"]["Electrostatic"]["ResponseCorrection"]["Models"][0];
  compact_model["FabricatedSurfaceMatrix"] = compact_fabricated_surface_path.string();
  compact_model["ThinSurfaceMatrix"] = compact_thin_surface_path.string();
  IoData compact_iodata(compact_config, false);
  LaplaceOperator compact_laplace_op(compact_iodata, meshes);
  SurfaceResponseOperator compact_response(compact_iodata, compact_laplace_op);

  const int size = laplace_op.GetH1Space().GetTrueVSize();
  Vector x(size), y(size), Cx, Cy, Ctx, compact_Cx;
  for (int i = 0; i < size; i++)
  {
    x(i) = 0.17 + 0.03 * (i + 1) * (Mpi::Rank(Mpi::World()) + 1);
    y(i) = -0.11 + 0.02 * (i + 2) * (Mpi::Rank(Mpi::World()) + 2);
  }
  response.Mult(x, Cx);
  response.Mult(y, Cy);
  response.MultTranspose(x, Ctx);
  compact_response.Mult(x, compact_Cx);

  double lhs = x * Cy;
  double rhs = Cx * y;
  double norm = Cx * Cx;
  Mpi::GlobalSum(1, &lhs, Mpi::World());
  Mpi::GlobalSum(1, &rhs, Mpi::World());
  Mpi::GlobalSum(1, &norm, Mpi::World());
  CHECK_THAT(lhs, WithinRel(rhs, 1.0e-12));
  CHECK(norm > 0.0);

  Ctx -= Cx;
  double transpose_error = Ctx * Ctx;
  Mpi::GlobalSum(1, &transpose_error, Mpi::World());
  CHECK(transpose_error == 0.0);
  compact_Cx -= Cx;
  double compact_error = compact_Cx * compact_Cx;
  Mpi::GlobalSum(1, &compact_error, Mpi::World());
  CHECK(compact_error == 0.0);

  Vector essential_values;
  Cx.GetSubVector(laplace_op.GetDbcTDofList(), essential_values);
  CHECK(essential_values.Normlinf() == 0.0);

  Vector prescribed(size), eliminated_rhs(size);
  prescribed = 0.0;
  const auto &essential = laplace_op.GetDbcTDofList();
  for (int i = 0; i < essential.Size(); i++)
  {
    prescribed(essential[i]) = 0.2 + 0.01 * (i + 1);
  }
  eliminated_rhs = 0.0;
  response.EliminateRHS(prescribed, eliminated_rhs);
  eliminated_rhs.GetSubVector(essential, essential_values);
  CHECK(essential_values.Normlinf() == 0.0);
  double rhs_norm = eliminated_rhs * eliminated_rhs;
  Mpi::GlobalSum(1, &rhs_norm, Mpi::World());
  CHECK(rhs_norm > 0.0);

  const auto energy = response.GetEnergyCorrection(x);
  REQUIRE(energy.interfaces.size() == 1);
  CHECK_THAT(energy.interfaces.at(4), WithinRel(energy.domain, 1.0e-12));
  const auto fabricated_energy = response.GetFabricatedSurfaceEnergy(x);
  REQUIRE(fabricated_energy.size() == 1);
  CHECK(fabricated_energy.at(4) > energy.interfaces.at(4));
  const auto compact_fabricated_energy = compact_response.GetFabricatedSurfaceEnergy(x);
  CHECK_THAT(compact_fabricated_energy.at(4), WithinRel(fabricated_energy.at(4), 1.0e-12));
  const auto local_electrostatic_response = response.GetElectrostaticResponse(x);
  CHECK_THAT(local_electrostatic_response.domain_correction,
             WithinRel(energy.domain, 1.0e-12));

  // The fixed-trace domain defect as a bilinear form (the corrected capacitance,
  // PostprocessCorrectedTerminals): on full vectors (x carries essential values), its
  // quadratic form is twice the fixed-trace domain correction, it is symmetric, its
  // off-diagonal value is the polarization of the quadratic form, and with the essential
  // rows zeroed it is the correction EliminateRHS subtracts (the same apply path).
  {
    Vector Dx, Dy;
    response.FixedTraceDomainDefectMult(x, Dx);
    response.FixedTraceDomainDefectMult(y, Dy);
    const double xDx = linalg::Dot<Vector>(Mpi::World(), x, Dx);
    const double yDy = linalg::Dot<Vector>(Mpi::World(), y, Dy);
    const double xDy = linalg::Dot<Vector>(Mpi::World(), x, Dy);
    const double yDx = linalg::Dot<Vector>(Mpi::World(), y, Dx);
    const auto response_y = response.GetElectrostaticResponse(y);
    CHECK(xDx != 0.0);
    CHECK_THAT(0.5 * xDx,
               WithinRel(local_electrostatic_response.domain_correction, 1.0e-12));
    CHECK_THAT(0.5 * yDy, WithinRel(response_y.domain_correction, 1.0e-12));
    CHECK(xDy != 0.0);
    CHECK_THAT(xDy, WithinRel(yDx, 1.0e-12));
    Vector sum(x);
    sum += y;
    const double polarization = response.GetElectrostaticResponse(sum).domain_correction -
                                local_electrostatic_response.domain_correction -
                                response_y.domain_correction;
    CHECK_THAT(xDy, WithinRel(polarization, 1.0e-9));
    Vector eliminated(size);
    eliminated = 0.0;
    response.EliminateRHS(x, eliminated);
    eliminated += Dx;
    eliminated.SetSubVector(laplace_op.GetDbcTDofList(), 0.0);
    double eliminated_error = eliminated * eliminated;
    Mpi::GlobalSum(1, &eliminated_error, Mpi::World());
    CHECK(eliminated_error == 0.0);
    // The masked operator drops the essential part of the trace: it is not the form.
    Dx.GetSubVector(laplace_op.GetDbcTDofList(), essential_values);
    double essential_norm = essential_values.Size() > 0 ? essential_values.Normlinf() : 0.0;
    Mpi::GlobalMax(1, &essential_norm, Mpi::World());
    CHECK(essential_norm > 0.0);
  }
  CHECK_THAT(local_electrostatic_response.fabricated_surface_energy.at(4),
             WithinRel(fabricated_energy.at(4), 1.0e-12));
  CHECK(std::isfinite(local_electrostatic_response.domain_correction_fixed_flux));
  CHECK(local_electrostatic_response.fabricated_surface_energy_fixed_flux.at(4) > 0.0);
  CHECK(local_electrostatic_response.maximum_trace_closure_spread > 0.0);
  // One patch and one mapped interface have no averaging: the response-weighted local
  // spread is exactly the interface aggregate, and all response weight either passes or
  // fails the local 5% limit.
  CHECK_THAT(local_electrostatic_response.response_weighted_trace_closure_spread,
             WithinRel(local_electrostatic_response.maximum_trace_closure_spread, 1.0e-12));
  CHECK_THAT(local_electrostatic_response.trace_closure_response_failure_fraction,
             WithinAbs(local_electrostatic_response.maximum_trace_closure_spread > 0.05
                           ? 1.0
                           : 0.0,
                       1.0e-12));
  CHECK(local_electrostatic_response.confident ==
        (local_electrostatic_response.maximum_trace_closure_spread <= 0.05));
  REQUIRE(local_electrostatic_response.model_contributions.size() == 1);
  const auto &model_contribution = local_electrostatic_response.model_contributions.front();
  CHECK(model_contribution.model == 1);
  CHECK_THAT(model_contribution.patch_count, WithinAbs(1.0, 1.0e-12));
  CHECK_THAT(model_contribution.domain_correction,
             WithinRel(local_electrostatic_response.domain_correction, 1.0e-12));
  CHECK_THAT(
      model_contribution.fabricated_surface_energy.at(4),
      WithinRel(local_electrostatic_response.fabricated_surface_energy.at(4), 1.0e-12));
  const auto statistics = response.GetStatistics();
  REQUIRE(statistics["ModelCatalog"].size() == 1);
  CHECK(statistics["ModelCatalog"][0]["Name"] == "model-1");
  CHECK(statistics["ModelCatalog"][0]["Topology"] == "Explicit");
  CHECK(statistics["ModelCatalog"][0]["PatchCount"] == 1);

  auto separated_config = config;
  separated_config["Solver"]["Electrostatic"]["ResponseCorrection"]["Patches"].push_back(
      {{"Model", 1},
       {"Origin", {0.35, 0.70, 0.0}},
       {"AxisU", {1.0, 0.0, 0.0}},
       {"AxisV", {0.0, 1.0, 0.0}},
       {"Reference", {0.0, 0.0, 0.0}}});
  IoData separated_iodata(separated_config, false);
  SurfaceResponseOperator separated_response(separated_iodata, laplace_op);
  CHECK(separated_response.GetPatchCount() == 2);
  CHECK(separated_response.GetBasisSize() == 8);

  auto overlapping_config = config;
  overlapping_config["Solver"]["Electrostatic"]["ResponseCorrection"]["Patches"].push_back(
      {{"Model", 1},
       {"Origin", {0.90, 0.35, 0.0}},
       {"AxisU", {1.0, 0.0, 0.0}},
       {"AxisV", {0.0, 1.0, 0.0}},
       {"Reference", {0.0, 0.0, 0.0}}});
  IoData overlapping_iodata(overlapping_config, false);
  CHECK_THROWS_WITH(SurfaceResponseOperator(overlapping_iodata, laplace_op),
                    Catch::Matchers::ContainsSubstring("coupled multi-edge coupon model"));

#endif
}

TEST_CASE_METHOD(test::SurfaceResponseFiles, "SurfaceResponseOperator automatic 2D library",
                 "[surfaceresponseoperator][2d][automatic][cache][Serial][Parallel]")
{
#if !defined(MFEM_USE_GSLIB)
  SKIP("SurfaceResponseOperator requires MFEM_USE_GSLIB");
#else
  json automatic_config = AutomaticConfig2D();
  IoData automatic_iodata(automatic_config, false);
  automatic_iodata.boundaries.cracked_attributes.insert(9);
  automatic_iodata.boundaries.cracked_attributes.insert(10);

  mfem::Mesh automatic_serial = MakeAutomatic2DMesh();
  auto automatic_parallel = std::make_unique<mfem::ParMesh>(Mpi::World(), automatic_serial);
  std::vector<std::unique_ptr<Mesh>> automatic_meshes;
  automatic_meshes.push_back(std::make_unique<Mesh>(std::move(automatic_parallel)));
  LaplaceOperator automatic_laplace(automatic_iodata, automatic_meshes);
  std::shared_ptr<const SurfaceResponseGeometry> automatic_geometry;
  SurfaceResponseOperator automatic_response(automatic_iodata, automatic_laplace,
                                             &automatic_geometry);
  REQUIRE(automatic_geometry);
  CHECK(automatic_response.GetPatchCount() == 2);
  CHECK(automatic_response.GetBasisSize() == 8);
  CHECK(automatic_response.HasSurfaceResponse());
  SurfaceResponseOperator cached_automatic_response(automatic_iodata, automatic_laplace,
                                                    &automatic_geometry);
  CHECK(cached_automatic_response.GetPatchCount() == automatic_response.GetPatchCount());
  CHECK(cached_automatic_response.GetBasisSize() == automatic_response.GetBasisSize());
  CHECK(cached_automatic_response.HasSurfaceResponse());
  const auto response_statistics = automatic_response.GetStatistics();
  const auto cached_response_statistics = cached_automatic_response.GetStatistics();
  CHECK(response_statistics["Version"] == 1);
  CHECK(response_statistics["Correction"]["Models"] == 1);
  CHECK(response_statistics["Correction"]["Patches"] == 2);
  CHECK(response_statistics["Correction"]["TraceCoefficients"] == 8);
  CHECK(response_statistics["Interpolation"]["PointQueries"]["Total"] ==
        response_statistics["Interpolation"]["StencilRows"]["Total"]);
  CHECK(response_statistics["Interpolation"]["StencilNonzeros"]["Total"] > 0);
  CHECK(response_statistics["Communication"]["PointSendItems"]["Total"] ==
        response_statistics["Communication"]["PointReceiveItems"]["Total"]);
  CHECK(cached_response_statistics["Geometry"] == response_statistics["Geometry"]);
  CHECK(cached_response_statistics["Matching"] == response_statistics["Matching"]);

  const auto requirements_path = temp.temp_dir / "surface-response-requirements.json";
  WriteSurfaceResponseRequirements(automatic_iodata, *automatic_meshes.back(),
                                   requirements_path.string());
  std::ifstream requirements_input(requirements_path);
  REQUIRE(requirements_input);
  const json requirements = json::parse(requirements_input);
  CHECK(requirements["Version"] == 1);
  CHECK(requirements["Complete"]);
  CHECK(requirements["MeshDimension"] == 2);
  CHECK_FALSE(requirements["Maxwell"]);
  CHECK(requirements["Summary"]["Counts"]["Exact"] == 2);
  CHECK(requirements["Summary"]["Counts"]["Missing"] == 0);
  REQUIRE(requirements.contains("Statistics"));
  CHECK(requirements["Statistics"]["Version"] == 1);
  CHECK(requirements["Statistics"]["Geometry"]["TargetGroups"] == 1);
  CHECK(requirements["Statistics"]["Geometry"]["EdgeSites2D"] == 2);
  REQUIRE(requirements["Requirements"].size() == 1);
  CHECK(requirements["Requirements"][0]["Topology"] == "IsolatedEdge");
  CHECK(requirements["Requirements"][0]["Status"] == "Exact");
  CHECK(requirements["Requirements"][0]["Count"] == 2);
  CHECK(requirements["Requirements"][0]["Interfaces"][0]["Target"] == 4);
  CHECK(requirements["Requirements"][0]["BoundaryCondition"] == json({{"Type", "PEC"}}));

  std::ifstream empty_library_input(library_path);
  REQUIRE(empty_library_input);
  auto empty_library = json::parse(empty_library_input);
  empty_library["Name"] = "unit-test-empty-process";
  empty_library["Models"] = json::array();
  const auto empty_library_path = temp.temp_dir / "surface-process-empty.json";
  std::ofstream empty_library_output(empty_library_path);
  empty_library_output << empty_library.dump(2) << "\n";
  empty_library_output.close();
  auto empty_library_config = automatic_config;
  empty_library_config["Solver"]["Electrostatic"]["ResponseCorrection"]["Library"] =
      empty_library_path.string();
  IoData empty_library_iodata(empty_library_config, false);
  empty_library_iodata.boundaries.cracked_attributes.insert(9);
  empty_library_iodata.boundaries.cracked_attributes.insert(10);
  const auto empty_requirements_path =
      temp.temp_dir / "surface-response-requirements-empty.json";
  WriteSurfaceResponseRequirements(empty_library_iodata, *automatic_meshes.back(),
                                   empty_requirements_path.string());
  std::ifstream empty_requirements_input(empty_requirements_path);
  REQUIRE(empty_requirements_input);
  const json empty_requirements = json::parse(empty_requirements_input);
  CHECK_FALSE(empty_requirements["Complete"]);
  CHECK(empty_requirements["Summary"]["Counts"]["Exact"] == 0);
  CHECK(empty_requirements["Summary"]["Counts"]["Missing"] == 2);
  REQUIRE(empty_requirements["Requirements"].size() == 1);
  CHECK(empty_requirements["Requirements"][0]["Topology"] == "IsolatedEdge");
  CHECK(empty_requirements["Requirements"][0]["Status"] == "Missing");
  CHECK(empty_requirements["Requirements"][0]["Count"] == 2);
  CHECK_THROWS(SurfaceResponseOperator(empty_library_iodata, automatic_laplace));

  auto thickness_mismatch_config = automatic_config;
  thickness_mismatch_config["Boundaries"]["Postprocessing"]["Dielectric"][0]["Thickness"] =
      0.003;
  IoData thickness_mismatch_iodata(thickness_mismatch_config, false);
  thickness_mismatch_iodata.boundaries.cracked_attributes.insert(9);
  thickness_mismatch_iodata.boundaries.cracked_attributes.insert(10);
  CHECK_THROWS_WITH(
      SurfaceResponseOperator(thickness_mismatch_iodata, automatic_laplace),
      Catch::Matchers::ContainsSubstring("does not match fabrication-process response "
                                         "library \"unit-test-process\" thickness"));

  auto permittivity_mismatch_config = automatic_config;
  permittivity_mismatch_config["Boundaries"]["Postprocessing"]["Dielectric"][0]
                              ["Permittivity"] = 4.1;
  IoData permittivity_mismatch_iodata(permittivity_mismatch_config, false);
  permittivity_mismatch_iodata.boundaries.cracked_attributes.insert(9);
  permittivity_mismatch_iodata.boundaries.cracked_attributes.insert(10);
  CHECK_THROWS_WITH(
      SurfaceResponseOperator(permittivity_mismatch_iodata, automatic_laplace),
      Catch::Matchers::ContainsSubstring("does not match fabrication-process response "
                                         "library \"unit-test-process\" permittivity"));

  auto missing_layer_config = automatic_config;
  missing_layer_config["Solver"]["Electrostatic"]["ResponseCorrection"]["Library"] =
      missing_layer_library_path.string();
  IoData missing_layer_iodata(missing_layer_config, false);
  missing_layer_iodata.boundaries.cracked_attributes.insert(9);
  missing_layer_iodata.boundaries.cracked_attributes.insert(10);
  CHECK_THROWS_WITH(
      SurfaceResponseOperator(missing_layer_iodata, automatic_laplace),
      Catch::Matchers::ContainsSubstring(
          "has a SA surface response but no matching Fabrication.InterfaceLayers entry"));

  auto legacy_config = automatic_config;
  legacy_config["Solver"]["Electrostatic"]["ResponseCorrection"]["Library"] =
      legacy_library_path.string();
  IoData legacy_iodata(legacy_config, false);
  legacy_iodata.boundaries.cracked_attributes.insert(9);
  legacy_iodata.boundaries.cracked_attributes.insert(10);
  SurfaceResponseOperator legacy_response(legacy_iodata, automatic_laplace);
  CHECK(legacy_response.GetPatchCount() == automatic_response.GetPatchCount());

  json exact_pair_config_2d = ExactPairConfig2D();
  IoData exact_pair_iodata_2d(exact_pair_config_2d, false);
  exact_pair_iodata_2d.boundaries.cracked_attributes.insert(9);
  exact_pair_iodata_2d.boundaries.cracked_attributes.insert(10);
  SurfaceResponseOperator exact_pair_response_2d(exact_pair_iodata_2d, automatic_laplace);
  CHECK(exact_pair_response_2d.GetPatchCount() == 1);
  CHECK(exact_pair_response_2d.GetBasisSize() == 4);

  auto interpolated_pair_config_2d = exact_pair_config_2d;
  interpolated_pair_config_2d["Solver"]["Electrostatic"]["ResponseCorrection"]["Library"] =
      interpolated_pair_library_2d_path.string();
  IoData interpolated_pair_iodata_2d(interpolated_pair_config_2d, false);
  interpolated_pair_iodata_2d.boundaries.cracked_attributes.insert(9);
  interpolated_pair_iodata_2d.boundaries.cracked_attributes.insert(10);
  SurfaceResponseOperator interpolated_pair_response_2d(interpolated_pair_iodata_2d,
                                                        automatic_laplace);
  CHECK(interpolated_pair_response_2d.GetPatchCount() == 2);
  CHECK(interpolated_pair_response_2d.GetBasisSize() == 8);
  CHECK_THAT(interpolated_pair_response_2d.GetPatchWeight(),
             WithinRel(exact_pair_response_2d.GetPatchWeight(), 1.0e-12));

  Vector pair_probe(automatic_laplace.GetH1Space().GetTrueVSize());
  for (int i = 0; i < pair_probe.Size(); i++)
  {
    pair_probe(i) = std::sin(0.11 * (i + 1 + Mpi::Rank(Mpi::World())));
  }
  Vector exact_pair_correction, interpolated_pair_correction;
  exact_pair_response_2d.Mult(pair_probe, exact_pair_correction);
  interpolated_pair_response_2d.Mult(pair_probe, interpolated_pair_correction);
  interpolated_pair_correction.Add(-1.0, exact_pair_correction);
  CHECK(linalg::Norml2(Mpi::World(), interpolated_pair_correction) <=
        1.0e-12 * std::max(linalg::Norml2(Mpi::World(), exact_pair_correction), 1.0e-300));

  auto invalid_depth_config = automatic_config;
  invalid_depth_config["Solver"]["Electrostatic"]["ResponseCorrection"]["Library"] =
      invalid_library_path.string();
  IoData invalid_depth_iodata(invalid_depth_config, false);
  invalid_depth_iodata.boundaries.cracked_attributes.insert(9);
  invalid_depth_iodata.boundaries.cracked_attributes.insert(10);
  CHECK_THROWS_WITH(SurfaceResponseOperator(invalid_depth_iodata, automatic_laplace),
                    Catch::Matchers::ContainsSubstring(
                        "Fabrication-process response-model CouponDepth must be positive"));

  auto disconnected_cluster_config = automatic_config;
  disconnected_cluster_config["Solver"]["Electrostatic"]["ResponseCorrection"]["Library"] =
      disconnected_cluster_library_3d_path.string();
  IoData disconnected_cluster_iodata(disconnected_cluster_config, false);
  disconnected_cluster_iodata.boundaries.cracked_attributes.insert(9);
  disconnected_cluster_iodata.boundaries.cracked_attributes.insert(10);
  CHECK_THROWS_WITH(SurfaceResponseOperator(disconnected_cluster_iodata, automatic_laplace),
                    Catch::Matchers::ContainsSubstring(
                        "OpenContourPaths must connect every conductor reference"));

#endif
}

TEST_CASE_METHOD(test::SurfaceResponseFiles, "SurfaceResponseOperator boundary-mode 2D",
                 "[surfaceresponseoperator][2d][boundarymode][impedance][Serial][Parallel]")
{
#if !defined(MFEM_USE_GSLIB)
  SKIP("SurfaceResponseOperator requires MFEM_USE_GSLIB");
#else
  json automatic_config = AutomaticConfig2D();
  IoData automatic_iodata(automatic_config, false);
  automatic_iodata.boundaries.cracked_attributes.insert(9);
  automatic_iodata.boundaries.cracked_attributes.insert(10);
  mfem::Mesh automatic_serial = MakeAutomatic2DMesh();
  auto automatic_parallel = std::make_unique<mfem::ParMesh>(Mpi::World(), automatic_serial);
  std::vector<std::unique_ptr<Mesh>> automatic_meshes;
  automatic_meshes.push_back(std::make_unique<Mesh>(std::move(automatic_parallel)));
  LaplaceOperator automatic_laplace(automatic_iodata, automatic_meshes);
  SurfaceResponseOperator automatic_response(automatic_iodata, automatic_laplace);

  json boundary_mode_config = BoundaryModeConfig2D();
  IoData boundary_mode_iodata(boundary_mode_config, false);
  boundary_mode_iodata.boundaries.cracked_attributes.insert(9);
  boundary_mode_iodata.boundaries.cracked_attributes.insert(10);
  MaterialOperator boundary_mode_material(boundary_mode_iodata, *automatic_meshes.back());
  BoundaryModeOperator boundary_mode_op(boundary_mode_iodata, automatic_meshes,
                                        boundary_mode_material);
  SurfaceResponseOperator boundary_mode_response(boundary_mode_iodata, boundary_mode_op);
  CHECK(boundary_mode_response.GetPatchCount() == automatic_response.GetPatchCount());
  CHECK(boundary_mode_response.GetBasisSize() == automatic_response.GetBasisSize());
  CHECK(boundary_mode_response.HasSurfaceResponse());
  CHECK(boundary_mode_response.GetTargetInterfaces() == std::set<int>{4});

  GridFunction boundary_mode_field(boundary_mode_op.GetNDSpace(), true);
  mfem::Vector boundary_mode_field_value(2);
  boundary_mode_field_value[0] = 0.7;
  boundary_mode_field_value[1] = -0.4;
  mfem::VectorConstantCoefficient boundary_mode_field_coefficient(
      boundary_mode_field_value);
  boundary_mode_field.Real().ProjectCoefficient(boundary_mode_field_coefficient);
  boundary_mode_field.Imag() = 0.0;
  const auto boundary_mode_result =
      boundary_mode_response.GetMaxwellResponse(boundary_mode_field, 0.0);
  CHECK(boundary_mode_result.fabricated_surface_energy.at(4) > 0.0);
  CHECK(boundary_mode_result.fabricated_surface_energy_fixed_flux.at(4) > 0.0);
  CHECK(boundary_mode_result.loop_residual < 1.0e-10);
  CHECK_THAT(boundary_mode_result.matched_length_fraction, WithinAbs(1.0, 1.0e-12));

  // Finite-impedance boundary-mode edges use the local metal surface as their voltage
  // reference. For a conservative field normal to the sheet this is gauge-equivalent to
  // the PEC reference displaced into the metal.
  auto impedance_boundary_mode_config = boundary_mode_config;
  impedance_boundary_mode_config["Boundaries"].erase("PEC");
  impedance_boundary_mode_config["Boundaries"]["Impedance"] = {
      {{"Attributes", {9, 10}}, {"Ls", 1.0e-13}}};
  impedance_boundary_mode_config["Solver"]["SurfaceResponseCorrection"]["Library"] =
      impedance_library_path.string();
  IoData impedance_boundary_mode_iodata(impedance_boundary_mode_config, false);
  impedance_boundary_mode_iodata.boundaries.cracked_attributes.insert(9);
  impedance_boundary_mode_iodata.boundaries.cracked_attributes.insert(10);
  MaterialOperator impedance_boundary_mode_material(impedance_boundary_mode_iodata,
                                                    *automatic_meshes.back());
  BoundaryModeOperator impedance_boundary_mode_op(
      impedance_boundary_mode_iodata, automatic_meshes, impedance_boundary_mode_material);
  SurfaceResponseOperator impedance_boundary_mode_response(impedance_boundary_mode_iodata,
                                                           impedance_boundary_mode_op);
  CHECK(impedance_boundary_mode_response.GetPatchCount() ==
        boundary_mode_response.GetPatchCount());
  GridFunction boundary_mode_normal_field(boundary_mode_op.GetNDSpace(), true);
  mfem::Vector boundary_mode_normal_value(2);
  boundary_mode_normal_value[0] = 0.0;
  boundary_mode_normal_value[1] = -1.0;
  mfem::VectorConstantCoefficient boundary_mode_normal_coefficient(
      boundary_mode_normal_value);
  boundary_mode_normal_field.Real().ProjectCoefficient(boundary_mode_normal_coefficient);
  boundary_mode_normal_field.Imag() = 0.0;
  GridFunction impedance_boundary_mode_normal_field(impedance_boundary_mode_op.GetNDSpace(),
                                                    true);
  impedance_boundary_mode_normal_field.Real().ProjectCoefficient(
      boundary_mode_normal_coefficient);
  impedance_boundary_mode_normal_field.Imag() = 0.0;
  const auto boundary_mode_normal_result =
      boundary_mode_response.GetMaxwellResponse(boundary_mode_normal_field, 0.0);
  const auto impedance_boundary_mode_result =
      impedance_boundary_mode_response.GetMaxwellResponse(
          impedance_boundary_mode_normal_field, 0.0);
  CHECK_THAT(impedance_boundary_mode_result.domain_correction,
             WithinRel(boundary_mode_normal_result.domain_correction, 1.0e-10));
  CHECK_THAT(
      impedance_boundary_mode_result.fabricated_surface_energy.at(4),
      WithinRel(boundary_mode_normal_result.fabricated_surface_energy.at(4), 1.0e-10));
  CHECK(impedance_boundary_mode_result.loop_residual < 1.0e-10);
  CHECK(impedance_boundary_mode_result.boundary_law_verified);

  auto nondimensionalized_impedance_config = impedance_boundary_mode_config;
  nondimensionalized_impedance_config["Model"]["L0"] = 1.0e-6;
  nondimensionalized_impedance_config["Model"]["Lc"] = 1.0;
  IoData nondimensionalized_impedance_iodata(nondimensionalized_impedance_config, false);
  nondimensionalized_impedance_iodata.boundaries.cracked_attributes.insert(9);
  nondimensionalized_impedance_iodata.boundaries.cracked_attributes.insert(10);
  auto nondimensionalized_serial = std::make_unique<mfem::Mesh>(automatic_serial);
  nondimensionalized_impedance_iodata.NondimensionalizeInputs(nondimensionalized_serial);
  auto nondimensionalized_parallel =
      std::make_unique<mfem::ParMesh>(Mpi::World(), *nondimensionalized_serial);
  std::vector<std::unique_ptr<Mesh>> nondimensionalized_meshes;
  nondimensionalized_meshes.push_back(
      std::make_unique<Mesh>(std::move(nondimensionalized_parallel)));
  MaterialOperator nondimensionalized_impedance_material(
      nondimensionalized_impedance_iodata, *nondimensionalized_meshes.back());
  BoundaryModeOperator nondimensionalized_impedance_op(
      nondimensionalized_impedance_iodata, nondimensionalized_meshes,
      nondimensionalized_impedance_material);
  SurfaceResponseOperator nondimensionalized_impedance_response(
      nondimensionalized_impedance_iodata, nondimensionalized_impedance_op);
  GridFunction nondimensionalized_impedance_field(
      nondimensionalized_impedance_op.GetNDSpace(), true);
  nondimensionalized_impedance_field.Real().ProjectCoefficient(
      boundary_mode_normal_coefficient);
  nondimensionalized_impedance_field.Imag() = 0.0;
  const auto nondimensionalized_impedance_result =
      nondimensionalized_impedance_response.GetMaxwellResponse(
          nondimensionalized_impedance_field, 0.0);
  CHECK(nondimensionalized_impedance_result.boundary_law_verified);
  CHECK(nondimensionalized_impedance_result.loop_residual < 1.0e-10);

  const auto impedance_requirements_path =
      temp.temp_dir / "surface-response-requirements-impedance.json";
  WriteSurfaceResponseRequirements(nondimensionalized_impedance_iodata,
                                   *nondimensionalized_meshes.back(),
                                   impedance_requirements_path.string());
  std::ifstream impedance_requirements_input(impedance_requirements_path);
  REQUIRE(impedance_requirements_input);
  const json impedance_requirements = json::parse(impedance_requirements_input);
  REQUIRE(impedance_requirements["Requirements"].size() == 1);
  const auto impedance_requirement_law =
      impedance_requirements["Requirements"][0]["BoundaryCondition"];
  CHECK(impedance_requirement_law["Type"] == "Impedance");
  CHECK(impedance_requirement_law["Rs"] == 0.0);
  CHECK_THAT(impedance_requirement_law["Ls"].get<double>(), WithinRel(1.0e-13, 1.0e-12));
  CHECK(impedance_requirement_law["Cs"] == 0.0);
  CHECK_FALSE(impedance_requirement_law.contains("Parameters"));
  CHECK_FALSE(impedance_requirement_law.contains("ParametersVerified"));

  const auto preflight_impedance_library_path =
      temp.temp_dir / "fabrication-process-impedance-preflight.json";
  if (Mpi::Root(Mpi::World()))
  {
    std::ifstream input(library_path);
    REQUIRE(input);
    auto preflight_library = json::parse(input);
    for (auto &model : preflight_library["Models"])
    {
      model["BoundaryCondition"] = impedance_requirement_law;
      model["BoundaryLawQualification"] = {{"Version", 1},
                                           {"Status", "Unqualified"},
                                           {"Calibration", "QuasiElectrostatic"},
                                           {"FrequencyUniversal", false}};
    }
    std::ofstream output(preflight_impedance_library_path);
    output << preflight_library.dump(2) << "\n";
  }
  Mpi::Barrier(Mpi::World());
  auto preflight_impedance_config = nondimensionalized_impedance_config;
  preflight_impedance_config["Solver"]["SurfaceResponseCorrection"]["Library"] =
      preflight_impedance_library_path.string();
  IoData preflight_impedance_iodata(preflight_impedance_config, false);
  preflight_impedance_iodata.boundaries.cracked_attributes.insert(9);
  preflight_impedance_iodata.boundaries.cracked_attributes.insert(10);
  std::unique_ptr<mfem::Mesh> no_impedance_mesh;
  preflight_impedance_iodata.NondimensionalizeInputs(no_impedance_mesh);
  SurfaceResponseOperator preflight_impedance_response(preflight_impedance_iodata,
                                                       nondimensionalized_impedance_op);
  CHECK(preflight_impedance_response.GetPatchCount() ==
        nondimensionalized_impedance_response.GetPatchCount());
  const auto preflight_impedance_result = preflight_impedance_response.GetMaxwellResponse(
      nondimensionalized_impedance_field, 0.0);
  CHECK_FALSE(preflight_impedance_result.boundary_law_verified);
  CHECK_FALSE(preflight_impedance_result.confident);

  auto impedance_mismatch_config = impedance_boundary_mode_config;
  impedance_mismatch_config["Boundaries"]["Impedance"][0]["Ls"] = 2.0e-13;
  IoData impedance_mismatch_iodata(impedance_mismatch_config, false);
  impedance_mismatch_iodata.boundaries.cracked_attributes.insert(9);
  impedance_mismatch_iodata.boundaries.cracked_attributes.insert(10);
  MaterialOperator impedance_mismatch_material(impedance_mismatch_iodata,
                                               *automatic_meshes.back());
  BoundaryModeOperator impedance_mismatch_op(impedance_mismatch_iodata, automatic_meshes,
                                             impedance_mismatch_material);
  CHECK_THROWS_WITH(
      SurfaceResponseOperator(impedance_mismatch_iodata, impedance_mismatch_op),
      Catch::Matchers::ContainsSubstring(
          "Automatic fabrication-process response matching failed"));

  auto invalid_boundary_law_config = impedance_boundary_mode_config;
  invalid_boundary_law_config["Solver"]["SurfaceResponseCorrection"]["Library"] =
      invalid_boundary_law_library_path.string();
  IoData invalid_boundary_law_iodata(invalid_boundary_law_config, false);
  CHECK_THROWS_WITH(
      SurfaceResponseOperator(invalid_boundary_law_iodata, impedance_boundary_mode_op),
      Catch::Matchers::ContainsSubstring(
          "Unknown fabrication-process response BoundaryCondition key \"Inductance\""));

  auto legacy_impedance_config = impedance_boundary_mode_config;
  legacy_impedance_config["Solver"]["SurfaceResponseCorrection"]["Library"] =
      legacy_impedance_library_path.string();
  IoData legacy_impedance_iodata(legacy_impedance_config, false);
  legacy_impedance_iodata.boundaries.cracked_attributes.insert(9);
  legacy_impedance_iodata.boundaries.cracked_attributes.insert(10);
  SurfaceResponseOperator legacy_impedance_response(legacy_impedance_iodata,
                                                    impedance_boundary_mode_op);
  const auto legacy_impedance_result = legacy_impedance_response.GetMaxwellResponse(
      impedance_boundary_mode_normal_field, 0.0);
  CHECK_FALSE(legacy_impedance_result.boundary_law_verified);
  CHECK_FALSE(legacy_impedance_result.closure_independent_confident);
  CHECK_FALSE(legacy_impedance_result.confident);

  auto conductivity_boundary_mode_config = boundary_mode_config;
  conductivity_boundary_mode_config["Boundaries"].erase("PEC");
  conductivity_boundary_mode_config["Boundaries"]["Conductivity"] = {
      {{"Attributes", {9, 10}},
       {"Conductivity", 5.8e7},
       {"Permeability", 1.2},
       {"Thickness", 1.0e-7},
       {"External", true}}};
  conductivity_boundary_mode_config["Solver"]["SurfaceResponseCorrection"]["Library"] =
      conductivity_library_path.string();
  IoData conductivity_boundary_mode_iodata(conductivity_boundary_mode_config, false);
  conductivity_boundary_mode_iodata.boundaries.cracked_attributes.insert(9);
  conductivity_boundary_mode_iodata.boundaries.cracked_attributes.insert(10);
  MaterialOperator conductivity_boundary_mode_material(conductivity_boundary_mode_iodata,
                                                       *automatic_meshes.back());
  BoundaryModeOperator conductivity_boundary_mode_op(conductivity_boundary_mode_iodata,
                                                     automatic_meshes,
                                                     conductivity_boundary_mode_material);
  SurfaceResponseOperator conductivity_boundary_mode_response(
      conductivity_boundary_mode_iodata, conductivity_boundary_mode_op);
  GridFunction conductivity_boundary_mode_field(conductivity_boundary_mode_op.GetNDSpace(),
                                                true);
  conductivity_boundary_mode_field.Real().ProjectCoefficient(
      boundary_mode_normal_coefficient);
  conductivity_boundary_mode_field.Imag() = 0.0;
  const auto conductivity_boundary_mode_result =
      conductivity_boundary_mode_response.GetMaxwellResponse(
          conductivity_boundary_mode_field, 0.0);
  CHECK(conductivity_boundary_mode_result.boundary_law_verified);
  CHECK(conductivity_boundary_mode_result.loop_residual < 1.0e-10);

  const auto conductivity_requirements_path =
      temp.temp_dir / "surface-response-requirements-conductivity.json";
  WriteSurfaceResponseRequirements(conductivity_boundary_mode_iodata,
                                   *automatic_meshes.back(),
                                   conductivity_requirements_path.string());
  std::ifstream conductivity_requirements_input(conductivity_requirements_path);
  REQUIRE(conductivity_requirements_input);
  const json conductivity_requirements = json::parse(conductivity_requirements_input);
  REQUIRE(conductivity_requirements["Requirements"].size() == 1);
  const auto conductivity_requirement_law =
      conductivity_requirements["Requirements"][0]["BoundaryCondition"];
  CHECK(conductivity_requirement_law["Type"] == "Conductivity");
  CHECK_THAT(conductivity_requirement_law["Conductivity"].get<double>(),
             WithinRel(5.8e7, 1.0e-12));
  CHECK_THAT(conductivity_requirement_law["Permeability"].get<double>(),
             WithinRel(1.2, 1.0e-12));
  CHECK_THAT(conductivity_requirement_law["Thickness"].get<double>(),
             WithinRel(2.0e-7, 1.0e-12));
  CHECK_FALSE(conductivity_requirement_law["External"].get<bool>());

  const auto preflight_conductivity_library_path =
      temp.temp_dir / "fabrication-process-conductivity-preflight.json";
  if (Mpi::Root(Mpi::World()))
  {
    std::ifstream input(library_path);
    REQUIRE(input);
    auto preflight_library = json::parse(input);
    for (auto &model : preflight_library["Models"])
    {
      model["BoundaryCondition"] = conductivity_requirement_law;
    }
    std::ofstream output(preflight_conductivity_library_path);
    output << preflight_library.dump(2) << "\n";
  }
  Mpi::Barrier(Mpi::World());
  auto preflight_conductivity_config = conductivity_boundary_mode_config;
  preflight_conductivity_config["Solver"]["SurfaceResponseCorrection"]["Library"] =
      preflight_conductivity_library_path.string();
  IoData preflight_conductivity_iodata(preflight_conductivity_config, false);
  preflight_conductivity_iodata.boundaries.cracked_attributes.insert(9);
  preflight_conductivity_iodata.boundaries.cracked_attributes.insert(10);
  SurfaceResponseOperator preflight_conductivity_response(preflight_conductivity_iodata,
                                                          conductivity_boundary_mode_op);
  CHECK(preflight_conductivity_response.GetPatchCount() ==
        conductivity_boundary_mode_response.GetPatchCount());

  auto rational_boundary_mode_config = boundary_mode_config;
  rational_boundary_mode_config["Boundaries"].erase("PEC");
  rational_boundary_mode_config["Boundaries"]["RationalImpedance"] = {
      {{"Attributes", {9, 10}},
       {"Numerator", rational_impedance_law["Numerator"]},
       {"Denominator", rational_impedance_law["Denominator"]}}};
  rational_boundary_mode_config["Solver"]["SurfaceResponseCorrection"]["Library"] =
      rational_impedance_library_path.string();
  IoData rational_boundary_mode_iodata(rational_boundary_mode_config, false);
  rational_boundary_mode_iodata.boundaries.cracked_attributes.insert(9);
  rational_boundary_mode_iodata.boundaries.cracked_attributes.insert(10);
  MaterialOperator rational_boundary_mode_material(rational_boundary_mode_iodata,
                                                   *automatic_meshes.back());
  BoundaryModeOperator rational_boundary_mode_op(
      rational_boundary_mode_iodata, automatic_meshes, rational_boundary_mode_material);
  SurfaceResponseOperator rational_boundary_mode_response(rational_boundary_mode_iodata,
                                                          rational_boundary_mode_op);
  GridFunction rational_boundary_mode_field(rational_boundary_mode_op.GetNDSpace(), true);
  rational_boundary_mode_field.Real().ProjectCoefficient(boundary_mode_normal_coefficient);
  rational_boundary_mode_field.Imag() = 0.0;
  const auto rational_boundary_mode_result =
      rational_boundary_mode_response.GetMaxwellResponse(rational_boundary_mode_field, 0.0);
  CHECK(rational_boundary_mode_result.boundary_law_verified);
  CHECK(rational_boundary_mode_result.loop_residual < 1.0e-10);

  auto nondimensionalized_rational_config = rational_boundary_mode_config;
  nondimensionalized_rational_config["Model"]["L0"] = 1.0e-6;
  nondimensionalized_rational_config["Model"]["Lc"] = 1.0;
  IoData nondimensionalized_rational_iodata(nondimensionalized_rational_config, false);
  nondimensionalized_rational_iodata.boundaries.cracked_attributes.insert(9);
  nondimensionalized_rational_iodata.boundaries.cracked_attributes.insert(10);
  auto nondimensionalized_rational_serial = std::make_unique<mfem::Mesh>(automatic_serial);
  nondimensionalized_rational_iodata.NondimensionalizeInputs(
      nondimensionalized_rational_serial);
  auto nondimensionalized_rational_parallel =
      std::make_unique<mfem::ParMesh>(Mpi::World(), *nondimensionalized_rational_serial);
  std::vector<std::unique_ptr<Mesh>> nondimensionalized_rational_meshes;
  nondimensionalized_rational_meshes.push_back(
      std::make_unique<Mesh>(std::move(nondimensionalized_rational_parallel)));
  MaterialOperator nondimensionalized_rational_material(
      nondimensionalized_rational_iodata, *nondimensionalized_rational_meshes.back());
  BoundaryModeOperator nondimensionalized_rational_op(nondimensionalized_rational_iodata,
                                                      nondimensionalized_rational_meshes,
                                                      nondimensionalized_rational_material);
  SurfaceResponseOperator nondimensionalized_rational_response(
      nondimensionalized_rational_iodata, nondimensionalized_rational_op);

  const auto rational_requirements_path =
      temp.temp_dir / "surface-response-requirements-rational-impedance.json";
  WriteSurfaceResponseRequirements(nondimensionalized_rational_iodata,
                                   *nondimensionalized_rational_meshes.back(),
                                   rational_requirements_path.string());
  std::ifstream rational_requirements_input(rational_requirements_path);
  REQUIRE(rational_requirements_input);
  const json rational_requirements = json::parse(rational_requirements_input);
  REQUIRE(rational_requirements["Requirements"].size() == 1);
  const auto rational_requirement_law =
      rational_requirements["Requirements"][0]["BoundaryCondition"];
  CHECK(rational_requirement_law["Type"] == "RationalImpedance");
  REQUIRE(rational_requirement_law["Numerator"].size() == 2);
  REQUIRE(rational_requirement_law["Denominator"].size() == 3);
  CHECK_FALSE(rational_requirement_law.contains("ParametersVerified"));

  const auto preflight_rational_library_path =
      temp.temp_dir / "fabrication-process-rational-impedance-preflight.json";
  if (Mpi::Root(Mpi::World()))
  {
    std::ifstream input(library_path);
    REQUIRE(input);
    auto preflight_library = json::parse(input);
    for (auto &model : preflight_library["Models"])
    {
      model["BoundaryCondition"] = rational_requirement_law;
    }
    std::ofstream output(preflight_rational_library_path);
    output << preflight_library.dump(2) << "\n";
  }
  Mpi::Barrier(Mpi::World());
  auto preflight_rational_config = nondimensionalized_rational_config;
  preflight_rational_config["Solver"]["SurfaceResponseCorrection"]["Library"] =
      preflight_rational_library_path.string();
  IoData preflight_rational_iodata(preflight_rational_config, false);
  preflight_rational_iodata.boundaries.cracked_attributes.insert(9);
  preflight_rational_iodata.boundaries.cracked_attributes.insert(10);
  std::unique_ptr<mfem::Mesh> no_rational_mesh;
  preflight_rational_iodata.NondimensionalizeInputs(no_rational_mesh);
  SurfaceResponseOperator preflight_rational_response(preflight_rational_iodata,
                                                      nondimensionalized_rational_op);
  CHECK(preflight_rational_response.GetPatchCount() ==
        nondimensionalized_rational_response.GetPatchCount());

#endif
}

TEST_CASE_METHOD(
    test::SurfaceResponseFiles, "SurfaceResponseOperator 2D pairs and clusters",
    "[surfaceresponseoperator][2d][cluster][boundarymode][cache][Serial][Parallel]")
{
#if !defined(MFEM_USE_GSLIB)
  SKIP("SurfaceResponseOperator requires MFEM_USE_GSLIB");
#else
  json automatic_config = AutomaticConfig2D();
  IoData automatic_iodata(automatic_config, false);
  automatic_iodata.boundaries.cracked_attributes.insert(9);
  automatic_iodata.boundaries.cracked_attributes.insert(10);
  mfem::Mesh automatic_serial = MakeAutomatic2DMesh();
  auto automatic_parallel = std::make_unique<mfem::ParMesh>(Mpi::World(), automatic_serial);
  std::vector<std::unique_ptr<Mesh>> automatic_meshes;
  automatic_meshes.push_back(std::make_unique<Mesh>(std::move(automatic_parallel)));
  LaplaceOperator automatic_laplace(automatic_iodata, automatic_meshes);

  json boundary_mode_config = BoundaryModeConfig2D();
  auto different_pair_config = boundary_mode_config;
  different_pair_config["Boundaries"]["Postprocessing"]["Dielectric"][0]["EdgeAttributes"] =
      {9, 10};
  different_pair_config["Boundaries"]["Postprocessing"]["Dielectric"][0]["EdgeDistances"] =
      {0.25};
  different_pair_config["Solver"]["SurfaceResponseCorrection"]["Library"] =
      different_pair_library_2d_path.string();
  IoData different_pair_iodata(different_pair_config, false);
  different_pair_iodata.boundaries.cracked_attributes.insert(9);
  different_pair_iodata.boundaries.cracked_attributes.insert(10);

  mfem::Mesh different_pair_serial =
      mfem::Mesh::MakeCartesian2D(10, 4, mfem::Element::TRIANGLE, false, 1.0, 1.0);
  for (int face = 0; face < different_pair_serial.GetNumFaces(); face++)
  {
    int element1, element2;
    different_pair_serial.GetFaceElements(face, &element1, &element2);
    if (element1 < 0 || element2 < 0)
    {
      continue;
    }
    mfem::Array<int> vertices;
    different_pair_serial.GetFaceVertices(face, vertices);
    if (vertices.Size() != 2)
    {
      continue;
    }
    const double *p0 = different_pair_serial.GetVertex(vertices[0]);
    const double *p1 = different_pair_serial.GetVertex(vertices[1]);
    const double xmin = std::min(p0[0], p1[0]);
    const double xmax = std::max(p0[0], p1[0]);
    if (std::abs(p0[1] - 0.5) < 1.0e-12 && std::abs(p1[1] - 0.5) < 1.0e-12 &&
        (xmax <= 0.3 + 1.0e-12 || xmin >= 0.7 - 1.0e-12))
    {
      different_pair_serial.AddBdrElement(
          different_pair_serial.GetFace(face)->Duplicate(&different_pair_serial));
      different_pair_serial.SetBdrAttribute(different_pair_serial.GetNBE() - 1,
                                            xmax <= 0.3 + 1.0e-12 ? 9 : 10);
    }
  }
  different_pair_serial.FinalizeTopology();
  different_pair_serial.Finalize();
  while (different_pair_serial.GetNE() < Mpi::Size(Mpi::World()))
  {
    different_pair_serial.UniformRefinement();
  }
  auto different_pair_parallel =
      std::make_unique<mfem::ParMesh>(Mpi::World(), different_pair_serial);
  std::vector<std::unique_ptr<Mesh>> different_pair_meshes;
  different_pair_meshes.push_back(
      std::make_unique<Mesh>(std::move(different_pair_parallel)));
  MaterialOperator different_pair_material(different_pair_iodata,
                                           *different_pair_meshes.back());
  BoundaryModeOperator different_pair_mode(different_pair_iodata, different_pair_meshes,
                                           different_pair_material);
  SurfaceResponseOperator different_pair_response(different_pair_iodata,
                                                  different_pair_mode);
  CHECK(different_pair_response.GetPatchCount() == 1);
  CHECK(different_pair_response.GetBasisSize() == 4);

  auto parallel_cluster_config_2d = automatic_config;
  parallel_cluster_config_2d["Boundaries"]["Ground"]["Attributes"] = {1, 3, 4, 9};
  parallel_cluster_config_2d["Boundaries"]["Terminal"][0]["Attributes"] = {2, 10};
  parallel_cluster_config_2d["Boundaries"]["Postprocessing"]["Dielectric"][0]
                            ["EdgeAttributes"] = {9, 10};
  parallel_cluster_config_2d["Boundaries"]["Postprocessing"]["Dielectric"][0]
                            ["EdgeDistances"] = {0.25};
  parallel_cluster_config_2d["Solver"]["Electrostatic"]["ResponseCorrection"]["Library"] =
      parallel_cluster_library_2d_path.string();
  IoData parallel_cluster_iodata_2d(parallel_cluster_config_2d, false);
  parallel_cluster_iodata_2d.boundaries.cracked_attributes.insert(9);
  parallel_cluster_iodata_2d.boundaries.cracked_attributes.insert(10);

  mfem::Mesh parallel_cluster_serial_2d =
      mfem::Mesh::MakeCartesian2D(10, 4, mfem::Element::TRIANGLE, false, 1.0, 1.0);
  for (int face = 0; face < parallel_cluster_serial_2d.GetNumFaces(); face++)
  {
    int element1, element2;
    parallel_cluster_serial_2d.GetFaceElements(face, &element1, &element2);
    if (element1 < 0 || element2 < 0)
    {
      continue;
    }
    mfem::Array<int> vertices;
    parallel_cluster_serial_2d.GetFaceVertices(face, vertices);
    if (vertices.Size() != 2)
    {
      continue;
    }
    const double *p0 = parallel_cluster_serial_2d.GetVertex(vertices[0]);
    const double *p1 = parallel_cluster_serial_2d.GetVertex(vertices[1]);
    const double xmin = std::min(p0[0], p1[0]);
    const double xmax = std::max(p0[0], p1[0]);
    const bool left_ground = xmax <= 0.2 + 1.0e-12;
    const bool trace = xmin >= 0.4 - 1.0e-12 && xmax <= 0.6 + 1.0e-12;
    const bool right_ground = xmin >= 0.8 - 1.0e-12;
    if (std::abs(p0[1] - 0.5) < 1.0e-12 && std::abs(p1[1] - 0.5) < 1.0e-12 &&
        (left_ground || trace || right_ground))
    {
      parallel_cluster_serial_2d.AddBdrElement(
          parallel_cluster_serial_2d.GetFace(face)->Duplicate(&parallel_cluster_serial_2d));
      parallel_cluster_serial_2d.SetBdrAttribute(parallel_cluster_serial_2d.GetNBE() - 1,
                                                 trace ? 10 : 9);
    }
  }
  parallel_cluster_serial_2d.FinalizeTopology();
  parallel_cluster_serial_2d.Finalize();
  while (parallel_cluster_serial_2d.GetNE() < Mpi::Size(Mpi::World()))
  {
    parallel_cluster_serial_2d.UniformRefinement();
  }
  auto parallel_cluster_parallel_2d =
      std::make_unique<mfem::ParMesh>(Mpi::World(), parallel_cluster_serial_2d);
  std::vector<std::unique_ptr<Mesh>> parallel_cluster_meshes_2d;
  parallel_cluster_meshes_2d.push_back(
      std::make_unique<Mesh>(std::move(parallel_cluster_parallel_2d)));
  LaplaceOperator parallel_cluster_laplace_2d(parallel_cluster_iodata_2d,
                                              parallel_cluster_meshes_2d);
  SurfaceResponseOperator parallel_cluster_response_2d(parallel_cluster_iodata_2d,
                                                       parallel_cluster_laplace_2d);
  CHECK(parallel_cluster_response_2d.GetPatchCount() == 1);
  CHECK(parallel_cluster_response_2d.GetBasisSize() == 5);

  auto mortar_cluster_config_2d = parallel_cluster_config_2d;
  mortar_cluster_config_2d["Solver"]["Electrostatic"]["ResponseCorrection"]
                          ["TraceCoupling"] = "SurfaceMortar";
  IoData mortar_cluster_iodata_2d(mortar_cluster_config_2d, false);
  mortar_cluster_iodata_2d.boundaries.cracked_attributes.insert(9);
  mortar_cluster_iodata_2d.boundaries.cracked_attributes.insert(10);
  LaplaceOperator mortar_cluster_laplace_2d(mortar_cluster_iodata_2d,
                                            parallel_cluster_meshes_2d);
  SurfaceResponseOperator mortar_cluster_response_2d(mortar_cluster_iodata_2d,
                                                     mortar_cluster_laplace_2d);
  mfem::FunctionCoefficient linear_trace_coefficient([](const mfem::Vector &x)
                                                     { return x[0] + 2.0 * x[1]; });
  mfem::ParGridFunction linear_trace_field(&parallel_cluster_laplace_2d.GetH1Space().Get());
  linear_trace_field.ProjectCoefficient(linear_trace_coefficient);
  Vector linear_trace_true;
  linear_trace_field.GetTrueDofs(linear_trace_true);
  const auto collocated_linear_response =
      parallel_cluster_response_2d.GetElectrostaticResponse(linear_trace_true);
  const auto mortar_linear_response =
      mortar_cluster_response_2d.GetElectrostaticResponse(linear_trace_true);
  CHECK_THAT(mortar_linear_response.domain_correction,
             WithinAbs(collocated_linear_response.domain_correction, 1.0e-12));
  CHECK_THAT(mortar_linear_response.domain_correction_fixed_flux,
             WithinAbs(collocated_linear_response.domain_correction_fixed_flux, 1.0e-12));
  CHECK_THAT(
      mortar_linear_response.fabricated_surface_energy.at(4),
      WithinAbs(collocated_linear_response.fabricated_surface_energy.at(4), 1.0e-12));

  auto parallel_cluster_boundary_config_2d = parallel_cluster_config_2d;
  parallel_cluster_boundary_config_2d["Problem"]["Type"] = "BoundaryMode";
  parallel_cluster_boundary_config_2d["Boundaries"].erase("Ground");
  parallel_cluster_boundary_config_2d["Boundaries"].erase("Terminal");
  parallel_cluster_boundary_config_2d["Boundaries"]["PEC"] = {{"Attributes", {9, 10}}};
  parallel_cluster_boundary_config_2d["Solver"] = {
      {"Order", 1},
      {"BoundaryMode", {{"Freq", 5.0}}},
      {"SurfaceResponseCorrection",
       {{"Library", parallel_cluster_library_2d_path.string()},
        {"UnmatchedPolicy", "Error"}}}};
  IoData parallel_cluster_boundary_iodata_2d(parallel_cluster_boundary_config_2d, false);
  parallel_cluster_boundary_iodata_2d.boundaries.cracked_attributes.insert(9);
  parallel_cluster_boundary_iodata_2d.boundaries.cracked_attributes.insert(10);
  MaterialOperator parallel_cluster_boundary_material_2d(
      parallel_cluster_boundary_iodata_2d, *parallel_cluster_meshes_2d.back());
  BoundaryModeOperator parallel_cluster_boundary_op_2d(
      parallel_cluster_boundary_iodata_2d, parallel_cluster_meshes_2d,
      parallel_cluster_boundary_material_2d);
  SurfaceResponseOperator parallel_cluster_boundary_response_2d(
      parallel_cluster_boundary_iodata_2d, parallel_cluster_boundary_op_2d);
  CHECK(parallel_cluster_boundary_response_2d.GetPatchCount() == 1);
  CHECK(parallel_cluster_boundary_response_2d.GetBasisSize() == 6);

  // Finite-impedance conductors use the same connected-component ownership as PEC.
  // In particular, the two disconnected ground strips share attribute 9 but remain
  // distinct conductors, selecting the three-conductor cluster response.
  auto impedance_cluster_boundary_config_2d = parallel_cluster_boundary_config_2d;
  impedance_cluster_boundary_config_2d["Boundaries"].erase("PEC");
  impedance_cluster_boundary_config_2d["Boundaries"]["Impedance"] = {
      {{"Attributes", {9, 10}}, {"Ls", 1.0e-13}}};
  impedance_cluster_boundary_config_2d["Solver"]["SurfaceResponseCorrection"]["Library"] =
      impedance_parallel_cluster_library_2d_path.string();
  IoData impedance_cluster_boundary_iodata_2d(impedance_cluster_boundary_config_2d, false);
  impedance_cluster_boundary_iodata_2d.boundaries.cracked_attributes.insert(9);
  impedance_cluster_boundary_iodata_2d.boundaries.cracked_attributes.insert(10);
  MaterialOperator impedance_cluster_boundary_material_2d(
      impedance_cluster_boundary_iodata_2d, *parallel_cluster_meshes_2d.back());
  BoundaryModeOperator impedance_cluster_boundary_op_2d(
      impedance_cluster_boundary_iodata_2d, parallel_cluster_meshes_2d,
      impedance_cluster_boundary_material_2d);
  SurfaceResponseOperator impedance_cluster_boundary_response_2d(
      impedance_cluster_boundary_iodata_2d, impedance_cluster_boundary_op_2d);
  CHECK(impedance_cluster_boundary_response_2d.GetPatchCount() == 1);
  CHECK(impedance_cluster_boundary_response_2d.GetBasisSize() == 6);

  mfem::FunctionCoefficient cluster_potential_coefficient(
      [](const mfem::Vector &x)
      {
        const double distance = std::abs(x[0] - 0.5);
        if (distance <= 0.1)
        {
          return 1.0;
        }
        if (distance >= 0.3)
        {
          return 0.0;
        }
        const double t = (distance - 0.1) / 0.2;
        return 1.0 - 3.0 * t * t + 2.0 * t * t * t;
      });
  mfem::ParGridFunction parallel_cluster_potential_2d(
      &parallel_cluster_laplace_2d.GetH1Space().Get());
  parallel_cluster_potential_2d.ProjectCoefficient(cluster_potential_coefficient);
  Vector parallel_cluster_potential_true_2d;
  parallel_cluster_potential_2d.GetTrueDofs(parallel_cluster_potential_true_2d);
  const auto parallel_cluster_electrostatic_result_2d =
      parallel_cluster_response_2d.GetElectrostaticResponse(
          parallel_cluster_potential_true_2d);

  // A disabled translational domain correction must remove the cluster from both the
  // self-consistent operator and corrected-domain accounting without disabling its
  // fabricated surface-energy evaluation.
  auto disabled_translational_config_2d = parallel_cluster_config_2d;
  disabled_translational_config_2d["Solver"]["Electrostatic"]["ResponseCorrection"]
                                  ["TranslationalDomainCorrection"] = "Disabled";
  IoData disabled_translational_iodata_2d(disabled_translational_config_2d, false);
  disabled_translational_iodata_2d.boundaries.cracked_attributes.insert(9);
  disabled_translational_iodata_2d.boundaries.cracked_attributes.insert(10);
  LaplaceOperator disabled_translational_laplace_2d(disabled_translational_iodata_2d,
                                                    parallel_cluster_meshes_2d);
  SurfaceResponseOperator disabled_translational_response_2d(
      disabled_translational_iodata_2d, disabled_translational_laplace_2d);
  Vector disabled_translational_action(parallel_cluster_potential_true_2d.Size());
  disabled_translational_response_2d.Mult(parallel_cluster_potential_true_2d,
                                          disabled_translational_action);
  CHECK(disabled_translational_action.Norml2() == 0.0);
  const auto disabled_translational_result_2d =
      disabled_translational_response_2d.GetElectrostaticResponse(
          parallel_cluster_potential_true_2d, false);
  CHECK(disabled_translational_result_2d.domain_correction == 0.0);
  CHECK_THAT(
      disabled_translational_result_2d.fabricated_surface_energy.at(4),
      WithinRel(parallel_cluster_electrostatic_result_2d.fabricated_surface_energy.at(4),
                1.0e-12));

  auto fixed_flux_translational_config_2d = parallel_cluster_config_2d;
  fixed_flux_translational_config_2d["Solver"]["Electrostatic"]["ResponseCorrection"]
                                    ["TranslationalDomainCorrection"] = "FixedFlux";
  IoData fixed_flux_translational_iodata_2d(fixed_flux_translational_config_2d, false);
  fixed_flux_translational_iodata_2d.boundaries.cracked_attributes.insert(9);
  fixed_flux_translational_iodata_2d.boundaries.cracked_attributes.insert(10);
  LaplaceOperator fixed_flux_translational_laplace_2d(fixed_flux_translational_iodata_2d,
                                                      parallel_cluster_meshes_2d);
  SurfaceResponseOperator fixed_flux_translational_response_2d(
      fixed_flux_translational_iodata_2d, fixed_flux_translational_laplace_2d);
  Vector fixed_flux_translational_action(parallel_cluster_potential_true_2d.Size());
  fixed_flux_translational_response_2d.Mult(parallel_cluster_potential_true_2d,
                                            fixed_flux_translational_action);
  const auto fixed_flux_translational_result_2d =
      fixed_flux_translational_response_2d.GetElectrostaticResponse(
          parallel_cluster_potential_true_2d, false);
  CHECK_THAT(
      fixed_flux_translational_result_2d.domain_correction,
      WithinRel(parallel_cluster_electrostatic_result_2d.domain_correction_fixed_flux,
                1.0e-12));
  const auto fixed_flux_translational_energy_2d =
      fixed_flux_translational_response_2d.GetEnergyCorrection(
          parallel_cluster_potential_true_2d);
  CHECK_THAT(fixed_flux_translational_energy_2d.domain,
             WithinRel(fixed_flux_translational_result_2d.domain_correction, 1.0e-12));
  // The corrected capacitance's bilinear defect is the FIXED-TRACE form irrespective of
  // the model's domain-coupling mode: half its quadratic form is the fixed-trace domain
  // correction (and not the fixed-flux one the operator applies).
  {
    Vector fixed_trace_action;
    fixed_flux_translational_response_2d.FixedTraceDomainDefectMult(
        parallel_cluster_potential_true_2d, fixed_trace_action);
    const double quadratic_form = linalg::Dot<Vector>(
        Mpi::World(), parallel_cluster_potential_true_2d, fixed_trace_action);
    CHECK_THAT(
        0.5 * quadratic_form,
        WithinRel(parallel_cluster_electrostatic_result_2d.domain_correction, 1.0e-12));
    CHECK_THAT(0.5 * quadratic_form,
               !WithinRel(fixed_flux_translational_result_2d.domain_correction, 1.0e-6));
  }

  GridFunction parallel_cluster_boundary_field_2d(
      parallel_cluster_boundary_op_2d.GetNDSpace(), true);
  mfem::ParGridFunction parallel_cluster_boundary_potential_2d(
      &parallel_cluster_boundary_op_2d.GetH1Space().Get());
  parallel_cluster_boundary_potential_2d.ProjectCoefficient(cluster_potential_coefficient);
  mfem::ParDiscreteLinearOperator parallel_cluster_gradient_2d(
      &parallel_cluster_boundary_op_2d.GetH1Space().Get(),
      &parallel_cluster_boundary_op_2d.GetNDSpace().Get());
  parallel_cluster_gradient_2d.AddDomainInterpolator(new mfem::GradientInterpolator());
  parallel_cluster_gradient_2d.Assemble();
  parallel_cluster_gradient_2d.Mult(parallel_cluster_boundary_potential_2d,
                                    parallel_cluster_boundary_field_2d.Real());
  parallel_cluster_boundary_field_2d.Real() *= -1.0;
  parallel_cluster_boundary_field_2d.Imag() = 0.0;
  const auto parallel_cluster_boundary_result_2d =
      parallel_cluster_boundary_response_2d.GetMaxwellResponse(
          parallel_cluster_boundary_field_2d, 0.0);
  CHECK(parallel_cluster_boundary_result_2d.loop_residual < 1.0e-10);
  CHECK_THAT(
      parallel_cluster_boundary_result_2d.domain_correction,
      WithinRel(parallel_cluster_electrostatic_result_2d.domain_correction, 1.0e-10));
  CHECK_THAT(
      parallel_cluster_boundary_result_2d.domain_correction_fixed_flux,
      WithinRel(parallel_cluster_electrostatic_result_2d.domain_correction_fixed_flux,
                1.0e-10));
  CHECK_THAT(
      parallel_cluster_boundary_result_2d.fabricated_surface_energy.at(4),
      WithinRel(parallel_cluster_electrostatic_result_2d.fabricated_surface_energy.at(4),
                1.0e-10));

  json exact_pair_config_2d = ExactPairConfig2D();
  // Regression (2D pair conductor references, 2026-09-26): the 1- / 2-site patches of
  // BuildAutomaticResponseData2D carry the selected model's conductor references — the
  // former "if empty" guard never fired because a ResponsePatchData starts with the
  // configuration default {{0, 0, 0}} (a same-conductor gap pair read its reference in the
  // gap, a different-conductor gap pair aborted on the reference count). Checked through
  // the response-geometry cache, which serialises every patch.
  {
    const auto cache_path = temp.temp_dir / "reference-pairs-2d-cache.json";
    auto PatchReferences = [&](const fs::path &path)
    {
      std::vector<std::vector<std::array<double, 3>>> references;
      if (Mpi::Root(Mpi::World()))
      {
        std::ifstream input(path);
        REQUIRE(input);
        const json cache = json::parse(input);
        for (const auto &patch : cache["Patches"])
        {
          references.push_back(
              patch["ConductorReferences"].get<std::vector<std::array<double, 3>>>());
        }
      }
      return references;
    };
    auto CheckReferences = [&](const std::vector<std::vector<std::array<double, 3>>> &found,
                               const std::vector<std::array<double, 3>> &expected)
    {
      if (Mpi::Root(Mpi::World()))
      {
        REQUIRE(!found.empty());
        for (const auto &references : found)
        {
          REQUIRE(references.size() == expected.size());
          for (std::size_t i = 0; i < expected.size(); i++)
          {
            for (int d = 0; d < 3; d++)
            {
              CHECK_THAT(references[i][d], WithinAbs(expected[i][d], 1.0e-12));
            }
          }
        }
      }
    };
    test::GeometryCacheEnvGuard cache_env(cache_path.string(), true);
    // Same-conductor strip pair (two sites) with an explicit model reference.
    auto reference_strip_config_2d = exact_pair_config_2d;
    reference_strip_config_2d["Solver"]["Electrostatic"]["ResponseCorrection"]["Library"] =
        reference_strip_library_2d_path.string();
    IoData reference_strip_iodata_2d(reference_strip_config_2d, false);
    reference_strip_iodata_2d.boundaries.cracked_attributes.insert(9);
    reference_strip_iodata_2d.boundaries.cracked_attributes.insert(10);
    SurfaceResponseOperator reference_strip_response_2d(reference_strip_iodata_2d,
                                                        automatic_laplace);
    CHECK(reference_strip_response_2d.GetPatchCount() == 1);
    Mpi::Barrier(Mpi::World());
    CheckReferences(PatchReferences(cache_path), {{0.05, -0.1, 0.0}});
    // Different-conductor gap pair (two sites of different conductors) with two references.
    auto reference_gap_config_2d = automatic_config;
    reference_gap_config_2d["Boundaries"]["Ground"]["Attributes"] = {1, 3, 4, 9};
    reference_gap_config_2d["Boundaries"]["Terminal"][0]["Attributes"] = {2, 10};
    reference_gap_config_2d["Boundaries"]["Postprocessing"]["Dielectric"][0]
                           ["EdgeAttributes"] = {9, 10};
    reference_gap_config_2d["Boundaries"]["Postprocessing"]["Dielectric"][0]
                           ["EdgeDistances"] = {0.25};
    reference_gap_config_2d["Solver"]["Electrostatic"]["ResponseCorrection"]["Library"] =
        reference_gap_library_2d_path.string();
    IoData reference_gap_iodata_2d(reference_gap_config_2d, false);
    reference_gap_iodata_2d.boundaries.cracked_attributes.insert(9);
    reference_gap_iodata_2d.boundaries.cracked_attributes.insert(10);
    mfem::Mesh reference_gap_serial =
        mfem::Mesh::MakeCartesian2D(10, 4, mfem::Element::TRIANGLE, false, 1.0, 1.0);
    for (int face = 0; face < reference_gap_serial.GetNumFaces(); face++)
    {
      int element1, element2;
      reference_gap_serial.GetFaceElements(face, &element1, &element2);
      if (element1 < 0 || element2 < 0)
      {
        continue;
      }
      mfem::Array<int> vertices;
      reference_gap_serial.GetFaceVertices(face, vertices);
      if (vertices.Size() != 2)
      {
        continue;
      }
      const double *p0 = reference_gap_serial.GetVertex(vertices[0]);
      const double *p1 = reference_gap_serial.GetVertex(vertices[1]);
      const double xmin = std::min(p0[0], p1[0]);
      const double xmax = std::max(p0[0], p1[0]);
      if (std::abs(p0[1] - 0.5) < 1.0e-12 && std::abs(p1[1] - 0.5) < 1.0e-12 &&
          (xmax <= 0.3 + 1.0e-12 || xmin >= 0.7 - 1.0e-12))
      {
        reference_gap_serial.AddBdrElement(
            reference_gap_serial.GetFace(face)->Duplicate(&reference_gap_serial));
        reference_gap_serial.SetBdrAttribute(reference_gap_serial.GetNBE() - 1,
                                             xmax <= 0.3 + 1.0e-12 ? 9 : 10);
      }
    }
    reference_gap_serial.FinalizeTopology();
    reference_gap_serial.Finalize();
    while (reference_gap_serial.GetNE() < Mpi::Size(Mpi::World()))
    {
      reference_gap_serial.UniformRefinement();
    }
    auto reference_gap_parallel =
        std::make_unique<mfem::ParMesh>(Mpi::World(), reference_gap_serial);
    std::vector<std::unique_ptr<Mesh>> reference_gap_meshes;
    reference_gap_meshes.push_back(
        std::make_unique<Mesh>(std::move(reference_gap_parallel)));
    LaplaceOperator reference_gap_laplace(reference_gap_iodata_2d, reference_gap_meshes);
    SurfaceResponseOperator reference_gap_response_2d(reference_gap_iodata_2d,
                                                      reference_gap_laplace);
    CHECK(reference_gap_response_2d.GetPatchCount() == 1);
    Mpi::Barrier(Mpi::World());
    CheckReferences(PatchReferences(cache_path), {{-0.2, 0.0, 0.0}, {0.2, 0.0, 0.0}});
  }

#endif
}

TEST_CASE_METHOD(
    test::SurfaceResponseFiles, "SurfaceResponseOperator axisymmetric curvature family",
    "[surfaceresponseoperator][2d][axisymmetric][curvature][cache][Serial][Parallel]")
{
#if !defined(MFEM_USE_GSLIB)
  SKIP("SurfaceResponseOperator requires MFEM_USE_GSLIB");
#else
  json automatic_config = AutomaticConfig2D();
  IoData automatic_iodata(automatic_config, false);
  automatic_iodata.boundaries.cracked_attributes.insert(9);
  automatic_iodata.boundaries.cracked_attributes.insert(10);
  mfem::Mesh automatic_serial = MakeAutomatic2DMesh();
  auto automatic_parallel = std::make_unique<mfem::ParMesh>(Mpi::World(), automatic_serial);
  std::vector<std::unique_ptr<Mesh>> automatic_meshes;
  automatic_meshes.push_back(std::make_unique<Mesh>(std::move(automatic_parallel)));
  LaplaceOperator automatic_laplace(automatic_iodata, automatic_meshes);

  // Curvature family on an axisymmetric (r, z) mesh (decision 92): the metal line from
  // r = 0.25 to 0.75 at R = 0.1 has a concave edge (gap toward -r) at kappa = 0.4 and a
  // convex one at kappa = 0.1333; each is corrected by a runtime model interpolated in
  // kappa between the straight anchor and the family's coupons (cubic Lagrange on the four
  // nearest nodes) with patch weight 2 pi r / CouponDepth; a missing convexity is never
  // straight.
  {
    const auto curved_library_path = temp.temp_dir / "fabrication-process-curved-2d.json";
    const auto convex_only_library_path =
        temp.temp_dir / "fabrication-process-curved-convex-only-2d.json";
    if (Mpi::Root(Mpi::World()))
    {
      std::ifstream input(library_path);
      json curved_library = json::parse(input);
      curved_library["Name"] = "unit-test-curved-2d";
      auto anchor = curved_library["Models"][0];
      anchor["CouponDepth"] = 1055.0;
      curved_library["Models"] = {anchor};
      // Every kappa node carries its own matrices (the anchor's scaled by 1 + kappa), so
      // a blend differs from the anchor and the cache round trip below is discriminating.
      auto ScaledMatrixFile =
          [&](const std::string &source, const std::string &name, double factor)
      {
        std::ifstream input(source);
        const auto target = temp.temp_dir / name;
        std::ofstream output(target);
        std::string line;
        std::getline(input, line);
        output << line << "\n";
        while (std::getline(input, line))
        {
          const auto comma = line.rfind(',');
          output << line.substr(0, comma + 1) << std::setprecision(17)
                 << factor * std::stod(line.substr(comma + 1)) << "\n";
        }
        return target.string();
      };
      for (const char *convexity : {"Convex", "Concave"})
      {
        for (const double kappa : {0.1, 0.25, 0.5, 0.8})
        {
          auto model = anchor;
          model["Name"] = std::string("curved-") + convexity + "-" + std::to_string(kappa);
          model["Topology"] = "CurvedEdge";
          model["Kappa"] = kappa;
          model["Convexity"] = convexity;
          model["CouponDepth"] = 2.0 * M_PI * 0.1 / kappa;
          const std::string suffix = model["Name"].get<std::string>();
          for (const char *key : {"FabricatedMatrix", "ThinMatrix",
                                  "FabricatedSurfaceMatrix", "ThinSurfaceMatrix"})
          {
            if (model.contains(key))
            {
              model[key] = ScaledMatrixFile(model[key].get<std::string>(),
                                            suffix + "-" + key + ".csv", 1.0 + kappa);
            }
          }
          curved_library["Models"].push_back(model);
        }
      }
      std::ofstream output(curved_library_path);
      output << curved_library.dump(2) << "\n";
      auto convex_only = curved_library;
      convex_only["Models"].erase(convex_only["Models"].begin() + 5,
                                  convex_only["Models"].end());
      std::ofstream convex_only_output(convex_only_library_path);
      convex_only_output << convex_only.dump(2) << "\n";
    }
    Mpi::Barrier(Mpi::World());
    auto curved_config = automatic_config;
    curved_config["Model"]["Axisymmetric"] = true;
    curved_config["Solver"]["Electrostatic"]["ResponseCorrection"]["Library"] =
        curved_library_path.string();
    IoData curved_iodata(curved_config, false);
    curved_iodata.boundaries.cracked_attributes.insert(9);
    curved_iodata.boundaries.cracked_attributes.insert(10);
    automatic_meshes.back()->SetAxisymmetric(true);
    SurfaceResponseOperator curved_response(curved_iodata, automatic_laplace);
    CHECK(curved_response.GetPatchCount() == 2);
    CHECK(curved_response.GetBasisSize() == 8);
    const auto curved_statistics = curved_response.GetStatistics();
    REQUIRE(curved_statistics["ModelCatalog"].size() == 2);
    std::set<std::string> catalog;
    for (const auto &entry : curved_statistics["ModelCatalog"])
    {
      catalog.insert(entry["Name"].get<std::string>());
      CHECK(entry["Topology"] == "curved edge");
    }
    CHECK(catalog.count("isolated@concave-kappa0.4-cubic") == 1);
    CHECK(catalog.count("isolated@convex-kappa0.133333333-cubic") == 1);

    // The response-geometry cache carries the blend of every runtime model (version 2):
    // a run reloading the cache reproduces the blended response exactly, and a cache
    // without the Blend entries is refused (never the anchor's straight matrices).
    {
      const auto cache_path = temp.temp_dir / "response-geometry-curved.json";
      test::GeometryCacheEnvGuard cache_env(cache_path.string(), true);
      SurfaceResponseOperator written_response(curved_iodata, automatic_laplace);
      Mpi::Barrier(Mpi::World());
      cache_env.DisableWrite();
      SurfaceResponseOperator loaded_response(curved_iodata, automatic_laplace);
      CHECK(loaded_response.GetPatchCount() == curved_response.GetPatchCount());
      CHECK(loaded_response.GetBasisSize() == curved_response.GetBasisSize());
      mfem::FunctionCoefficient curved_trace_coefficient(
          [](const mfem::Vector &x) { return x[0] + 0.5 * x[1] * x[1]; });
      mfem::ParGridFunction curved_trace_field(&automatic_laplace.GetH1Space().Get());
      curved_trace_field.ProjectCoefficient(curved_trace_coefficient);
      Vector curved_trace_true;
      curved_trace_field.GetTrueDofs(curved_trace_true);
      const auto fresh = curved_response.GetElectrostaticResponse(curved_trace_true);
      const auto reloaded = loaded_response.GetElectrostaticResponse(curved_trace_true);
      CHECK(fresh.domain_correction != 0.0);
      CHECK_THAT(reloaded.domain_correction, WithinRel(fresh.domain_correction, 1.0e-12));
      for (const auto &[interface, energy] : fresh.fabricated_surface_energy)
      {
        CHECK_THAT(reloaded.fabricated_surface_energy.at(interface),
                   WithinRel(energy, 1.0e-12));
      }
      std::ifstream cache_input(cache_path);
      REQUIRE(cache_input);
      json cache = json::parse(cache_input);
      CHECK(cache["Version"] == 14);
      REQUIRE(cache["Models"].size() == 2);
      for (auto &model : cache["Models"])
      {
        // Four sources per cubic blend; the weights carry the coupon-depth rescaling
        // (Lagrange weight x D_anchor / D_node), so they do not sum to one.
        REQUIRE(model["Blend"].size() == 4);
        for (const auto &source : model["Blend"])
        {
          CHECK(std::isfinite(source["Weight"].get<double>()));
          CHECK(!source["FabricatedMatrix"].get<std::string>().empty());
        }
        model.erase("Blend");
      }
      const auto stripped_path = temp.temp_dir / "response-geometry-curved-stripped.json";
      if (Mpi::Root(Mpi::World()))
      {
        std::ofstream stripped(stripped_path);
        stripped << cache.dump(2) << "\n";
      }
      Mpi::Barrier(Mpi::World());
      test::GeometryCacheEnvGuard stripped_cache_env(stripped_path.string(), false);
      CHECK_THROWS(SurfaceResponseOperator(curved_iodata, automatic_laplace));
    }

    const auto curved_requirements_path =
        temp.temp_dir / "surface-response-requirements-curved.json";
    WriteSurfaceResponseRequirements(curved_iodata, *automatic_meshes.back(),
                                     curved_requirements_path.string());
    std::ifstream curved_input(curved_requirements_path);
    REQUIRE(curved_input);
    const json curved_requirements = json::parse(curved_input);
    auto Lagrange = [](const std::vector<double> &nodes, double x)
    {
      std::vector<double> weights;
      for (std::size_t i = 0; i < nodes.size(); i++)
      {
        double w = 1.0;
        for (std::size_t j = 0; j < nodes.size(); j++)
        {
          if (j != i)
          {
            w *= (x - nodes[j]) / (nodes[i] - nodes[j]);
          }
        }
        weights.push_back(w);
      }
      return weights;
    };
    int curved_records = 0;
    for (const auto &record : curved_requirements["Requirements"])
    {
      if (record["Topology"] != "CurvedEdge")
      {
        continue;
      }
      curved_records++;
      CHECK(record["Status"] == "Interpolated");
      CHECK(record["Geometry"]["InterpolationRule"] == "cubic");
      CHECK_THAT(record["Geometry"]["KappaMax"].get<double>(), WithinAbs(0.8, 1.0e-12));
      const double kappa = record["Geometry"]["Kappa"].get<double>();
      const bool convex = record["Geometry"]["Convexity"] == "Convex";
      CHECK_THAT(kappa, WithinRel(convex ? 0.1 / 0.75 : 0.1 / 0.25, 1.0e-9));
      const std::vector<double> nodes = convex ? std::vector<double>{0.0, 0.1, 0.25, 0.5}
                                               : std::vector<double>{0.1, 0.25, 0.5, 0.8};
      const auto weights = Lagrange(nodes, kappa);
      REQUIRE(record["SelectedModels"].size() == 4);
      double sum = 0.0;
      for (std::size_t i = 0; i < 4; i++)
      {
        const auto &model = record["SelectedModels"][i];
        CHECK_THAT(model["Weight"].get<double>(), WithinAbs(weights[i], 1.0e-6));
        sum += model["Weight"].get<double>();
        if (nodes[i] == 0.0)
        {
          CHECK(model["Name"] == "isolated");
        }
        else
        {
          CHECK(model["Topology"] == "CurvedEdge");
        }
      }
      CHECK_THAT(sum, WithinAbs(1.0, 1.0e-6));
    }
    CHECK(curved_records == 2);

    // Exact-node and first-order rules (the same metal line at other matching radii; the
    // rule is keyed to 1 / StraightBendRadiusOverR = 0.1, not to the library's smallest
    // node). R = 0.075: the convex edge at r = 0.75 sits on the kappa = 0.1 node (exact),
    // the concave one at r = 0.25 is cubic at kappa 0.3. R = 0.025: the convex edge at
    // kappa = 1/30 is linear between the anchor and the 0.1 node (weights 2/3, 1/3), the
    // concave one exact at 0.1; without the 0.1 node the convex edge is refused (never a
    // linear rule on the next node).
    auto RescaledCurvedRun =
        [&](double R, const std::vector<std::string> &drop_models, const std::string &tag)
    {
      const auto path = temp.temp_dir / ("fabrication-process-curved-" + tag + ".json");
      if (Mpi::Root(Mpi::World()))
      {
        std::ifstream input(curved_library_path);
        json library = json::parse(input);
        library["MatchingRadius"] = R;
        json models = json::array();
        for (const auto &model : library["Models"])
        {
          if (std::find(drop_models.begin(), drop_models.end(),
                        model["Name"].get<std::string>()) == drop_models.end())
          {
            models.push_back(model);
          }
        }
        library["Models"] = models;
        std::ofstream output(path);
        output << library.dump(2) << "\n";
      }
      Mpi::Barrier(Mpi::World());
      auto config = curved_config;
      config["Boundaries"]["Postprocessing"]["Dielectric"][0]["EdgeDistances"] = {R};
      config["Solver"]["Electrostatic"]["ResponseCorrection"]["Library"] = path.string();
      return config;
    };
    auto CurvedRecords = [&](const json &config, const std::string &tag)
    {
      IoData iodata(config, false);
      iodata.boundaries.cracked_attributes.insert(9);
      iodata.boundaries.cracked_attributes.insert(10);
      const auto path =
          temp.temp_dir / ("surface-response-requirements-curved-" + tag + ".json");
      WriteSurfaceResponseRequirements(iodata, *automatic_meshes.back(), path.string());
      std::ifstream input(path);
      REQUIRE(input);
      const json requirements = json::parse(input);
      std::map<std::string, json> records;  // by convexity
      for (const auto &record : requirements["Requirements"])
      {
        if (record["Topology"] == "CurvedEdge")
        {
          records[record["Geometry"]["Convexity"].get<std::string>()] = record;
        }
      }
      REQUIRE(records.size() == 2);
      return records;
    };
    {
      const auto records = CurvedRecords(RescaledCurvedRun(0.075, {}, "exact"), "exact");
      const auto &convex = records.at("Convex");
      CHECK(convex["Status"] == "Exact");
      CHECK(convex["Geometry"]["InterpolationRule"] == "exact");
      CHECK_THAT(convex["Geometry"]["Kappa"].get<double>(), WithinAbs(0.1, 1.0e-12));
      CHECK_THAT(convex["Geometry"]["FirstOrderKappa"].get<double>(),
                 WithinAbs(0.1, 1.0e-12));
      REQUIRE(convex["SelectedModels"].size() == 1);
      CHECK(convex["SelectedModels"][0]["Name"] == "curved-Convex-0.100000");
      CHECK_THAT(convex["SelectedModels"][0]["Weight"].get<double>(),
                 WithinAbs(1.0, 1.0e-12));
      const auto &concave = records.at("Concave");
      CHECK(concave["Geometry"]["InterpolationRule"] == "cubic");
      CHECK_THAT(concave["Geometry"]["Kappa"].get<double>(), WithinAbs(0.3, 1.0e-12));
      CHECK(concave["SelectedModels"].size() == 4);
    }
    {
      const auto records = CurvedRecords(RescaledCurvedRun(0.025, {}, "linear"), "linear");
      const auto &convex = records.at("Convex");
      CHECK(convex["Status"] == "Interpolated");
      CHECK(convex["Geometry"]["InterpolationRule"] == "linear");
      CHECK_THAT(convex["Geometry"]["Kappa"].get<double>(), WithinAbs(1.0 / 30.0, 1.0e-12));
      REQUIRE(convex["SelectedModels"].size() == 2);
      CHECK(convex["SelectedModels"][0]["Name"] == "isolated");
      CHECK_THAT(convex["SelectedModels"][0]["Weight"].get<double>(),
                 WithinAbs(2.0 / 3.0, 1.0e-9));
      CHECK(convex["SelectedModels"][1]["Name"] == "curved-Convex-0.100000");
      CHECK_THAT(convex["SelectedModels"][1]["Weight"].get<double>(),
                 WithinAbs(1.0 / 3.0, 1.0e-9));
      const auto &concave = records.at("Concave");
      CHECK(concave["Geometry"]["InterpolationRule"] == "exact");
      CHECK(concave["SelectedModels"][0]["Name"] == "curved-Concave-0.100000");
    }
    {
      // No convex 0.1 node: kappa = 1/30 is in the first-order regime and is refused
      // (the 0.25 node is not used for a linear rule).
      const auto config =
          RescaledCurvedRun(0.025, {"curved-Convex-0.100000"}, "no-first-order-node");
      IoData iodata(config, false);
      iodata.boundaries.cracked_attributes.insert(9);
      iodata.boundaries.cracked_attributes.insert(10);
      CHECK_THROWS_WITH(SurfaceResponseOperator(iodata, automatic_laplace),
                        Catch::Matchers::ContainsSubstring("response matching failed"));
    }

    // Without concave coupons the concave edge is unmatched (reported, not straight).
    auto convex_only_config = curved_config;
    convex_only_config["Solver"]["Electrostatic"]["ResponseCorrection"]["Library"] =
        convex_only_library_path.string();
    IoData convex_only_iodata(convex_only_config, false);
    convex_only_iodata.boundaries.cracked_attributes.insert(9);
    convex_only_iodata.boundaries.cracked_attributes.insert(10);
    CHECK_THROWS_WITH(SurfaceResponseOperator(convex_only_iodata, automatic_laplace),
                      Catch::Matchers::ContainsSubstring("response matching failed"));
    automatic_meshes.back()->SetAxisymmetric(false);
  }

#endif
}

TEST_CASE_METHOD(test::SurfaceResponseFiles,
                 "SurfaceResponseOperator 3D CPW legacy electrostatic",
                 "[surfaceresponseoperator][3d][legacy][cpw][Serial][Parallel]")
{
#if !defined(MFEM_USE_GSLIB)
  SKIP("SurfaceResponseOperator requires MFEM_USE_GSLIB");
#else
  json config_3d = Cpw3dLegacyConfig();
  IoData iodata_3d(config_3d, false);
  auto mesh_3d = mesh::ReadMesh(iodata_3d, Mpi::World());
  const auto geometry_3d = ExtractMetalEdgeGeometry(
      *mesh_3d, iodata_3d.boundaries, JointNoiseExtractionFor(iodata_3d.boundaries));
  const auto segment_indices =
      GetInterfaceMetalEdgeSegmentIndices(geometry_3d, 1, InterfaceDielectric::SA);
  double physical_edge_length = 0.0;
  std::set<std::size_t> physical_edge_vertices;
  for (const std::size_t segment_index : segment_indices)
  {
    const auto &segment = geometry_3d.segments[segment_index];
    physical_edge_vertices.insert(segment.vertices.begin(), segment.vertices.end());
    const auto &p0 = geometry_3d.vertices[segment.vertices[0]].coordinate;
    const auto &p1 = geometry_3d.vertices[segment.vertices[1]].coordinate;
    double length_squared = 0.0;
    for (int d = 0; d < 3; d++)
    {
      length_squared += (p1[d] - p0[d]) * (p1[d] - p0[d]);
    }
    physical_edge_length += std::sqrt(length_squared);
  }
  const int physical_endpoint_count = static_cast<int>(std::count_if(
      physical_edge_vertices.begin(), physical_edge_vertices.end(),
      [&](std::size_t vertex)
      {
        return geometry_3d.vertices[vertex].physical_type == MetalEdgeVertexType::ENDPOINT;
      }));
  REQUIRE(physical_endpoint_count > 0);
  CHECK(static_cast<int>(
            std::count_if(physical_edge_vertices.begin(), physical_edge_vertices.end(),
                          [&](std::size_t vertex)
                          {
                            return geometry_3d.vertices[vertex].physical_type ==
                                       MetalEdgeVertexType::ENDPOINT &&
                                   geometry_3d.vertices[vertex].on_truncation_boundary;
                          })) == physical_endpoint_count);
  std::vector<std::unique_ptr<Mesh>> meshes_3d;
  meshes_3d.push_back(std::make_unique<Mesh>(std::move(mesh_3d)));
  LaplaceOperator laplace_3d(iodata_3d, meshes_3d);
  SurfaceResponseOperator response_3d(iodata_3d, laplace_3d);
  const auto &line_rule =
      mfem::IntRules.Get(mfem::Geometry::SEGMENT, 2 * iodata_3d.solver.order);
  CHECK(response_3d.GetPatchCount() ==
        static_cast<int>(segment_indices.size()) * line_rule.GetNPoints());
  CHECK(response_3d.GetBasisSize() == 4 * response_3d.GetPatchCount());
  CHECK_THAT(response_3d.GetPatchWeight(), WithinRel(0.5 * physical_edge_length, 1.0e-12));
  CHECK(response_3d.HasSurfaceResponse());

  const auto endpoint_requirements_path =
      temp.temp_dir / "surface-response-requirements-endpoint.json";
  WriteSurfaceResponseRequirements(iodata_3d, *meshes_3d.back(),
                                   endpoint_requirements_path.string());
  std::ifstream endpoint_requirements_input(endpoint_requirements_path);
  REQUIRE(endpoint_requirements_input);
  const json endpoint_requirements = json::parse(endpoint_requirements_input);
  const auto endpoint_requirement = std::find_if(
      endpoint_requirements["Requirements"].begin(),
      endpoint_requirements["Requirements"].end(),
      [](const auto &requirement)
      {
        return requirement["Topology"] == "Endpoint" && requirement["Status"] == "Missing";
      });
  CHECK(endpoint_requirement == endpoint_requirements["Requirements"].end());

  auto endpoint_config_3d = config_3d;
  endpoint_config_3d["Solver"]["Electrostatic"]["ResponseCorrection"]["Library"] =
      endpoint_library_3d_path.string();
  IoData endpoint_iodata_3d(endpoint_config_3d, false);
  ShareCrackedAttributes(endpoint_iodata_3d, iodata_3d);
  LaplaceOperator endpoint_laplace_3d(endpoint_iodata_3d, meshes_3d);
  SurfaceResponseOperator endpoint_response_3d(endpoint_iodata_3d, endpoint_laplace_3d);
  CHECK(endpoint_response_3d.GetPatchCount() == response_3d.GetPatchCount());
  CHECK(endpoint_response_3d.GetBasisSize() == response_3d.GetBasisSize());
  CHECK_THAT(endpoint_response_3d.GetPatchWeight(),
             WithinRel(response_3d.GetPatchWeight(), 1.0e-12));

  auto coupled_config_3d = config_3d;
  for (auto &interface : coupled_config_3d["Boundaries"]["Postprocessing"]["Dielectric"])
  {
    interface["EdgeDistances"] = {0.2, 2.0, 7.0};
  }
  coupled_config_3d["Solver"]["Electrostatic"]["ResponseCorrection"]["Library"] =
      coupled_library_3d_path.string();
  IoData coupled_iodata_3d(coupled_config_3d, false);
  coupled_iodata_3d.boundaries.cracked_attributes.insert(1);
  coupled_iodata_3d.boundaries.cracked_attributes.insert(2);
  ShareCrackedAttributes(coupled_iodata_3d, iodata_3d);
  LaplaceOperator coupled_laplace_3d(coupled_iodata_3d, meshes_3d);
  SurfaceResponseOperator coupled_response_3d(coupled_iodata_3d, coupled_laplace_3d);
  CHECK(coupled_response_3d.GetPatchCount() ==
        static_cast<int>(segment_indices.size() / 2) * line_rule.GetNPoints());
  CHECK(coupled_response_3d.GetBasisSize() == 5 * coupled_response_3d.GetPatchCount());
  CHECK_THAT(coupled_response_3d.GetPatchWeight(),
             WithinRel(0.25 * physical_edge_length, 1.0e-12));
  CHECK(coupled_response_3d.HasSurfaceResponse());

  auto missing_pair_config_3d = coupled_config_3d;
  missing_pair_config_3d["Solver"]["Electrostatic"]["ResponseCorrection"]["Library"] =
      missing_pair_library_3d_path.string();
  missing_pair_config_3d["Solver"]["Electrostatic"]["ResponseCorrection"]
                        ["UnmatchedPolicy"] = "Warn";
  IoData missing_pair_iodata_3d(missing_pair_config_3d, false);
  missing_pair_iodata_3d.boundaries.cracked_attributes.insert(1);
  missing_pair_iodata_3d.boundaries.cracked_attributes.insert(2);
  ShareCrackedAttributes(missing_pair_iodata_3d, iodata_3d);
  const auto missing_pair_requirements_path =
      temp.temp_dir / "surface-response-requirements-missing-pair.json";
  WriteSurfaceResponseRequirements(missing_pair_iodata_3d, *meshes_3d.back(),
                                   missing_pair_requirements_path.string());
  std::ifstream missing_pair_requirements_input(missing_pair_requirements_path);
  REQUIRE(missing_pair_requirements_input);
  const json missing_pair_requirements = json::parse(missing_pair_requirements_input);
  CHECK_FALSE(missing_pair_requirements["Complete"]);
  CHECK(missing_pair_requirements["Summary"]["Counts"]["Missing"].get<int>() > 0);
  const auto missing_pair_requirement =
      std::find_if(missing_pair_requirements["Requirements"].begin(),
                   missing_pair_requirements["Requirements"].end(),
                   [](const auto &requirement)
                   {
                     return requirement["Topology"] == "DifferentConductorGap" &&
                            requirement["Status"] == "Missing";
                   });
  REQUIRE(missing_pair_requirement != missing_pair_requirements["Requirements"].end());
  CHECK((*missing_pair_requirement)["Geometry"]["EdgeCount"] == 2);
  // The version-1 Separation is derived from the signature's SeparationOverR on the
  // recorded 1e-6 R signature grid (R = 7), hence the 7e-6 tolerance.
  CHECK_THAT((*missing_pair_requirement)["Geometry"]["Separation"].get<double>(),
             WithinAbs(12.0, 7.0e-6));

  auto interpolated_coupled_config_3d = coupled_config_3d;
  interpolated_coupled_config_3d["Solver"]["Electrostatic"]["ResponseCorrection"]
                                ["Library"] = interpolated_coupled_library_3d_path.string();
  IoData interpolated_coupled_iodata_3d(interpolated_coupled_config_3d, false);
  interpolated_coupled_iodata_3d.boundaries.cracked_attributes.insert(1);
  interpolated_coupled_iodata_3d.boundaries.cracked_attributes.insert(2);
  ShareCrackedAttributes(interpolated_coupled_iodata_3d, iodata_3d);
  LaplaceOperator interpolated_coupled_laplace_3d(interpolated_coupled_iodata_3d,
                                                  meshes_3d);
  SurfaceResponseOperator interpolated_coupled_response_3d(interpolated_coupled_iodata_3d,
                                                           interpolated_coupled_laplace_3d);
  CHECK(interpolated_coupled_response_3d.GetPatchCount() ==
        2 * coupled_response_3d.GetPatchCount());
  CHECK(interpolated_coupled_response_3d.GetBasisSize() ==
        2 * coupled_response_3d.GetBasisSize());
  CHECK_THAT(interpolated_coupled_response_3d.GetPatchWeight(),
             WithinRel(coupled_response_3d.GetPatchWeight(), 1.0e-12));

  auto parallel_cluster_config_3d = coupled_config_3d;
  for (auto &interface :
       parallel_cluster_config_3d["Boundaries"]["Postprocessing"]["Dielectric"])
  {
    interface["EdgeDistances"] = {0.2, 2.0, 11.0};
  }
  parallel_cluster_config_3d["Solver"]["Electrostatic"]["ResponseCorrection"]["Library"] =
      parallel_cluster_library_3d_path.string();
  IoData parallel_cluster_iodata_3d(parallel_cluster_config_3d, false);
  parallel_cluster_iodata_3d.boundaries.cracked_attributes.insert(1);
  parallel_cluster_iodata_3d.boundaries.cracked_attributes.insert(2);
  ShareCrackedAttributes(parallel_cluster_iodata_3d, iodata_3d);
  LaplaceOperator parallel_cluster_laplace_3d(parallel_cluster_iodata_3d, meshes_3d);
  SurfaceResponseOperator parallel_cluster_response_3d(parallel_cluster_iodata_3d,
                                                       parallel_cluster_laplace_3d);
  CHECK(parallel_cluster_response_3d.GetPatchCount() ==
        static_cast<int>(segment_indices.size() / 4) * line_rule.GetNPoints());
  // Conductor identity is the edge-connected metal component (phase 3): the two ground
  // planes cut by the simulation box are distinct conductors, so the CPW cross-section is
  // the three-conductor four-edge cluster (6 basis functions per patch, not 5).
  CHECK(parallel_cluster_response_3d.GetBasisSize() ==
        6 * parallel_cluster_response_3d.GetPatchCount());
  CHECK_THAT(parallel_cluster_response_3d.GetPatchWeight(),
             WithinRel(0.125 * physical_edge_length, 1.0e-12));

  const auto parallel_cluster_requirements_path =
      temp.temp_dir / "surface-response-requirements-parallel-cluster.json";
  WriteSurfaceResponseRequirements(parallel_cluster_iodata_3d, *meshes_3d.back(),
                                   parallel_cluster_requirements_path.string());
  std::ifstream parallel_cluster_requirements_input(parallel_cluster_requirements_path);
  REQUIRE(parallel_cluster_requirements_input);
  const json parallel_cluster_requirements =
      json::parse(parallel_cluster_requirements_input);
  const auto parallel_cluster_requirement =
      std::find_if(parallel_cluster_requirements["Requirements"].begin(),
                   parallel_cluster_requirements["Requirements"].end(),
                   [](const auto &requirement)
                   {
                     return requirement["Topology"] == "ParallelEdgeCluster" &&
                            requirement["Status"] == "Exact";
                   });
  REQUIRE(parallel_cluster_requirement !=
          parallel_cluster_requirements["Requirements"].end());
  CHECK((*parallel_cluster_requirement)["Geometry"]["EdgeCount"] == 4);
  CHECK((*parallel_cluster_requirement)["Geometry"]["Edges"].size() == 4);
  std::set<int> parallel_cluster_conductors;
  for (const auto &edge : (*parallel_cluster_requirement)["Geometry"]["Edges"])
  {
    const int conductor = edge["Conductor"];
    CHECK(conductor > 0);
    parallel_cluster_conductors.insert(conductor);
  }
  // Three geometric conductors (left ground, centre strip, right ground), see above.
  CHECK(parallel_cluster_conductors == std::set<int>{1, 2, 3});
  CHECK((*parallel_cluster_requirement)["TotalEdgeLength"].get<double>() > 0.0);
  // Version-1 record semantics (review fix-4 m-F): Count is the number of mesh segments
  // carrying the coupon (per_segment classes), Instances the number of FEATURE instances
  // grouped into the coupon and DistinctSignatures the number of distinct signatures among
  // them (the library builder's coupon Instances). One 4-edge stack along the CPW: one
  // instance with one signature over several segments.
  CHECK((*parallel_cluster_requirement)["Instances"] == 1);
  CHECK((*parallel_cluster_requirement)["DistinctSignatures"] == 1);
  CHECK((*parallel_cluster_requirement)["Count"].get<int>() >
        (*parallel_cluster_requirement)["Instances"].get<int>());
  CHECK((*parallel_cluster_requirement)["ParameterSpread"] == 0.0);

  // An exact multi-edge coupon is self-contained. It must not require redundant
  // two-edge models for every pair in the active cluster.
  auto parallel_cluster_only_config_3d = parallel_cluster_config_3d;
  parallel_cluster_only_config_3d["Solver"]["Electrostatic"]["ResponseCorrection"]
                                 ["Library"] =
                                     parallel_cluster_only_library_3d_path.string();
  IoData parallel_cluster_only_iodata_3d(parallel_cluster_only_config_3d, false);
  parallel_cluster_only_iodata_3d.boundaries.cracked_attributes.insert(1);
  parallel_cluster_only_iodata_3d.boundaries.cracked_attributes.insert(2);
  ShareCrackedAttributes(parallel_cluster_only_iodata_3d, iodata_3d);
  LaplaceOperator parallel_cluster_only_laplace_3d(parallel_cluster_only_iodata_3d,
                                                   meshes_3d);
  SurfaceResponseOperator parallel_cluster_only_response_3d(
      parallel_cluster_only_iodata_3d, parallel_cluster_only_laplace_3d);
  CHECK(parallel_cluster_only_response_3d.GetPatchCount() ==
        parallel_cluster_response_3d.GetPatchCount());
  CHECK(parallel_cluster_only_response_3d.GetBasisSize() ==
        parallel_cluster_response_3d.GetBasisSize());

#endif
}

TEST_CASE_METHOD(test::SurfaceResponseFiles, "SurfaceResponseOperator 3D CPW Maxwell",
                 "[surfaceresponseoperator][3d][maxwell][cpw][Serial][Parallel]")
{
#if !defined(MFEM_USE_GSLIB)
  SKIP("SurfaceResponseOperator requires MFEM_USE_GSLIB");
#else
  // The electrostatic (Laplace, order 3) counterparts of the Maxwell operators below, built
  // on one mesh read (their own checks are the legacy electrostatic case's).
  json config_3d = Cpw3dLegacyConfig();
  IoData iodata_3d(config_3d, false);
  auto mesh_3d = mesh::ReadMesh(iodata_3d, Mpi::World());
  const auto geometry_3d = ExtractMetalEdgeGeometry(
      *mesh_3d, iodata_3d.boundaries, JointNoiseExtractionFor(iodata_3d.boundaries));
  const auto segment_indices =
      GetInterfaceMetalEdgeSegmentIndices(geometry_3d, 1, InterfaceDielectric::SA);
  std::vector<std::unique_ptr<Mesh>> meshes_3d;
  meshes_3d.push_back(std::make_unique<Mesh>(std::move(mesh_3d)));
  LaplaceOperator laplace_3d(iodata_3d, meshes_3d);
  SurfaceResponseOperator response_3d(iodata_3d, laplace_3d);

  auto coupled_config_3d = config_3d;
  for (auto &interface : coupled_config_3d["Boundaries"]["Postprocessing"]["Dielectric"])
  {
    interface["EdgeDistances"] = {0.2, 2.0, 7.0};
  }
  coupled_config_3d["Solver"]["Electrostatic"]["ResponseCorrection"]["Library"] =
      coupled_library_3d_path.string();
  IoData coupled_iodata_3d(coupled_config_3d, false);
  coupled_iodata_3d.boundaries.cracked_attributes.insert(1);
  coupled_iodata_3d.boundaries.cracked_attributes.insert(2);
  ShareCrackedAttributes(coupled_iodata_3d, iodata_3d);
  LaplaceOperator coupled_laplace_3d(coupled_iodata_3d, meshes_3d);
  SurfaceResponseOperator coupled_response_3d(coupled_iodata_3d, coupled_laplace_3d);

  auto interpolated_coupled_config_3d = coupled_config_3d;
  interpolated_coupled_config_3d["Solver"]["Electrostatic"]["ResponseCorrection"]
                                ["Library"] = interpolated_coupled_library_3d_path.string();
  IoData interpolated_coupled_iodata_3d(interpolated_coupled_config_3d, false);
  interpolated_coupled_iodata_3d.boundaries.cracked_attributes.insert(1);
  interpolated_coupled_iodata_3d.boundaries.cracked_attributes.insert(2);
  ShareCrackedAttributes(interpolated_coupled_iodata_3d, iodata_3d);
  LaplaceOperator interpolated_coupled_laplace_3d(interpolated_coupled_iodata_3d,
                                                  meshes_3d);
  SurfaceResponseOperator interpolated_coupled_response_3d(interpolated_coupled_iodata_3d,
                                                           interpolated_coupled_laplace_3d);

  auto parallel_cluster_config_3d = coupled_config_3d;
  for (auto &interface :
       parallel_cluster_config_3d["Boundaries"]["Postprocessing"]["Dielectric"])
  {
    interface["EdgeDistances"] = {0.2, 2.0, 11.0};
  }
  parallel_cluster_config_3d["Solver"]["Electrostatic"]["ResponseCorrection"]["Library"] =
      parallel_cluster_library_3d_path.string();
  IoData parallel_cluster_iodata_3d(parallel_cluster_config_3d, false);
  parallel_cluster_iodata_3d.boundaries.cracked_attributes.insert(1);
  parallel_cluster_iodata_3d.boundaries.cracked_attributes.insert(2);
  ShareCrackedAttributes(parallel_cluster_iodata_3d, iodata_3d);
  LaplaceOperator parallel_cluster_laplace_3d(parallel_cluster_iodata_3d, meshes_3d);
  SurfaceResponseOperator parallel_cluster_response_3d(parallel_cluster_iodata_3d,
                                                       parallel_cluster_laplace_3d);

  // The Maxwell (SpaceOperator, order 1) problems: one mesh read for the PEC metal
  // configurations, one for the finite-impedance metal.
  auto maxwell_config_3d = config_3d;
  maxwell_config_3d["Problem"]["Type"] = "Eigenmode";
  maxwell_config_3d["Boundaries"]["Ground"]["Attributes"] = {1, 2};
  maxwell_config_3d["Boundaries"].erase("Terminal");
  maxwell_config_3d["Solver"] = {{"Order", 1},
                                 {"Eigenmode", {{"Target", 1.0}}},
                                 {"SurfaceResponseCorrection",
                                  {{"Library", library_3d_path.string()},
                                   {"TargetInterfaces", {1, 2, 3}},
                                   {"UnmatchedPolicy", "Error"},
                                   {"PatchConstruction", "Legacy"}}}};
  IoData maxwell_iodata_3d(maxwell_config_3d, false);
  auto maxwell_mesh_3d = mesh::ReadMesh(maxwell_iodata_3d, Mpi::World());
  std::vector<std::unique_ptr<Mesh>> maxwell_meshes_3d;
  maxwell_meshes_3d.push_back(std::make_unique<Mesh>(std::move(maxwell_mesh_3d)));
  SpaceOperator maxwell_space_3d(maxwell_iodata_3d, maxwell_meshes_3d);
  SurfaceResponseOperator maxwell_response_3d(maxwell_iodata_3d, maxwell_space_3d);
  CHECK(maxwell_response_3d.GetPatchCount() > 0);
  CHECK(maxwell_response_3d.HasSurfaceResponse());
  CHECK(maxwell_response_3d.GetTargetInterfaces() == std::set<int>{1, 2, 3});

  auto endpoint_maxwell_config_3d = maxwell_config_3d;
  endpoint_maxwell_config_3d["Solver"]["SurfaceResponseCorrection"]["Library"] =
      endpoint_library_3d_path.string();
  IoData endpoint_maxwell_iodata_3d(endpoint_maxwell_config_3d, false);
  ShareCrackedAttributes(endpoint_maxwell_iodata_3d, maxwell_iodata_3d);
  SpaceOperator endpoint_maxwell_space_3d(endpoint_maxwell_iodata_3d, maxwell_meshes_3d);
  SurfaceResponseOperator endpoint_maxwell_response_3d(endpoint_maxwell_iodata_3d,
                                                       endpoint_maxwell_space_3d);
  CHECK(endpoint_maxwell_response_3d.GetPatchCount() ==
        maxwell_response_3d.GetPatchCount());

  // A finite-impedance metal does not provide an interior equipotential anchor. The
  // automatic Maxwell trace instead references the local metal edge point and must still
  // reproduce a conservative quasi-electrostatic field.
  auto impedance_maxwell_config_3d = maxwell_config_3d;
  impedance_maxwell_config_3d["Boundaries"].erase("Ground");
  impedance_maxwell_config_3d["Boundaries"]["Impedance"] = {
      {{"Attributes", {1, 2}}, {"Ls", 1.0e-13}}};
  impedance_maxwell_config_3d["Solver"]["SurfaceResponseCorrection"]["Library"] =
      impedance_library_3d_path.string();
  IoData impedance_maxwell_iodata_3d(impedance_maxwell_config_3d, false);
  auto impedance_maxwell_mesh_3d =
      mesh::ReadMesh(impedance_maxwell_iodata_3d, Mpi::World());
  std::vector<std::unique_ptr<Mesh>> impedance_maxwell_meshes_3d;
  impedance_maxwell_meshes_3d.push_back(
      std::make_unique<Mesh>(std::move(impedance_maxwell_mesh_3d)));
  SpaceOperator impedance_maxwell_space_3d(impedance_maxwell_iodata_3d,
                                           impedance_maxwell_meshes_3d);
  SurfaceResponseOperator impedance_maxwell_response_3d(impedance_maxwell_iodata_3d,
                                                        impedance_maxwell_space_3d);
  CHECK(impedance_maxwell_response_3d.GetPatchCount() ==
        maxwell_response_3d.GetPatchCount());
  GridFunction impedance_maxwell_field(impedance_maxwell_space_3d.GetNDSpace(), true);
  mfem::Vector impedance_gradient(3);
  impedance_gradient = 0.0;
  impedance_gradient[1] = -1.0;
  mfem::VectorConstantCoefficient impedance_gradient_coefficient(impedance_gradient);
  impedance_maxwell_field.Real().ProjectCoefficient(impedance_gradient_coefficient);
  impedance_maxwell_field.Imag() = 0.0;
  const auto impedance_gradient_response =
      impedance_maxwell_response_3d.GetMaxwellResponse(impedance_maxwell_field, 0.0);
  CHECK(std::abs(impedance_gradient_response.domain_correction) > 0.0);
  CHECK(impedance_gradient_response.loop_residual < 1.0e-10);
  CHECK(impedance_gradient_response.boundary_law_verified);

  GridFunction maxwell_field(maxwell_space_3d.GetNDSpace(), true);
  mfem::Vector constant_field(3);
  constant_field[0] = 0.7;
  constant_field[1] = -0.4;
  constant_field[2] = 0.2;
  mfem::VectorConstantCoefficient field_coefficient(constant_field);
  maxwell_field.Real().ProjectCoefficient(field_coefficient);
  maxwell_field.Imag() = 0.0;
  GridFunction endpoint_maxwell_field(endpoint_maxwell_space_3d.GetNDSpace(), true);
  endpoint_maxwell_field.Real().ProjectCoefficient(field_coefficient);
  endpoint_maxwell_field.Imag() = 0.0;
  const auto endpoint_maxwell_result =
      endpoint_maxwell_response_3d.GetMaxwellResponse(endpoint_maxwell_field, 0.0);
  CHECK(endpoint_maxwell_result.loop_residual < 1.0e-10);
  const auto real_response = maxwell_response_3d.GetMaxwellResponse(maxwell_field, 0.0);
  CHECK(std::abs(real_response.domain_correction) > 0.0);
  CHECK(std::abs(real_response.domain_correction_fixed_flux) > 0.0);
  REQUIRE(real_response.fabricated_surface_energy.size() == 3);
  REQUIRE(real_response.fabricated_surface_energy_fixed_flux.size() == 3);
  for (const auto &[interface, energy] : real_response.fabricated_surface_energy_fixed_flux)
  {
    (void)interface;
    CHECK(energy > 0.0);
  }
  CHECK(real_response.loop_residual < 1.0e-10);
  CHECK(real_response.response_weighted_loop_residual < 1.0e-10);
  CHECK(real_response.loop_response_failure_fraction == 0.0);
  CHECK(real_response.kR == 0.0);
  CHECK(real_response.maximum_trace_closure_spread > 0.0);
  CHECK(real_response.response_weighted_trace_closure_spread > 0.0);
  CHECK(real_response.trace_closure_response_failure_fraction >= 0.0);
  CHECK(real_response.trace_closure_response_failure_fraction <= 1.0);
  CHECK_THAT(real_response.matched_length_fraction, WithinAbs(1.0, 1.0e-12));
  REQUIRE(real_response.matched_length_fraction_by_interface.size() == 3);
  for (const auto &[interface, fraction] :
       real_response.matched_length_fraction_by_interface)
  {
    (void)interface;
    CHECK_THAT(fraction, WithinAbs(1.0, 1.0e-12));
  }

  // The contour reconstruction is also the trace operator used by the self-consistent
  // Maxwell mass correction. Verify its energy identity and transpose symmetry.
  Vector field_true, correction_true;
  maxwell_field.Real().GetTrueDofs(field_true);
  field_true.SetSubVector(maxwell_space_3d.GetNDDbcTDofLists().back(), 0.0);
  maxwell_field.Real().SetFromTrueDofs(field_true);
  const auto operator_response = maxwell_response_3d.GetMaxwellResponse(maxwell_field, 0.0);
  maxwell_response_3d.Mult(field_true, correction_true);
  CHECK_THAT(0.5 * linalg::Dot(Mpi::World(), field_true, correction_true),
             WithinRel(operator_response.domain_correction, 1.0e-10));

  Vector probe(field_true.Size()), correction_probe;
  auto *probe_data = probe.HostWrite();
  for (int i = 0; i < probe.Size(); i++)
  {
    probe_data[i] = std::sin(0.37 * (i + 1 + 11 * Mpi::Rank(Mpi::World())));
  }
  probe.SetSubVector(maxwell_space_3d.GetNDDbcTDofLists().back(), 0.0);
  maxwell_response_3d.Mult(probe, correction_probe);
  CHECK_THAT(linalg::Dot(Mpi::World(), probe, correction_true),
             WithinRel(linalg::Dot(Mpi::World(), field_true, correction_probe), 1.0e-10));

  // A Maxwell trace reconstructed from E = -grad(V) must reproduce the H1 coupon trace
  // relative to the PEC. This catches a missing contour-voltage gauge even when the
  // contour-loop residual and complex-field scaling tests pass.
  mfem::ParGridFunction potential(&laplace_3d.GetH1Space().Get());
  mfem::FunctionCoefficient potential_coefficient([](const mfem::Vector &x)
                                                  { return x[1]; });
  potential.ProjectCoefficient(potential_coefficient);
  Vector potential_true;
  potential.GetTrueDofs(potential_true);
  mfem::Vector normal_field(3);
  normal_field = 0.0;
  normal_field[1] = -1.0;
  mfem::VectorConstantCoefficient normal_field_coefficient(normal_field);
  maxwell_field.Real().ProjectCoefficient(normal_field_coefficient);
  maxwell_field.Imag() = 0.0;

  const auto electrostatic_correction = response_3d.GetEnergyCorrection(potential_true);
  const auto electrostatic_surfaces =
      response_3d.GetFabricatedSurfaceEnergy(potential_true);
  const auto electrostatic_response = response_3d.GetElectrostaticResponse(potential_true);
  const auto gradient_response = maxwell_response_3d.GetMaxwellResponse(maxwell_field, 0.0);
  CHECK_THAT(gradient_response.domain_correction,
             WithinRel(electrostatic_correction.domain, 1.0e-10));
  REQUIRE(gradient_response.fabricated_surface_energy.size() ==
          electrostatic_surfaces.size());
  for (const auto &[interface, energy] : electrostatic_surfaces)
  {
    CHECK_THAT(gradient_response.fabricated_surface_energy.at(interface),
               WithinRel(energy, 1.0e-10));
    CHECK_THAT(
        gradient_response.fabricated_surface_energy_fixed_flux.at(interface),
        WithinRel(electrostatic_response.fabricated_surface_energy_fixed_flux.at(interface),
                  1.0e-10));
  }
  CHECK_THAT(gradient_response.domain_correction_fixed_flux,
             WithinRel(electrostatic_response.domain_correction_fixed_flux, 1.0e-10));
  CHECK_THAT(
      gradient_response.response_weighted_trace_closure_spread,
      WithinRel(electrostatic_response.response_weighted_trace_closure_spread, 1.0e-10));
  CHECK_THAT(
      gradient_response.trace_closure_response_failure_fraction,
      WithinAbs(electrostatic_response.trace_closure_response_failure_fraction, 1.0e-12));
  CHECK(gradient_response.loop_residual < 1.0e-10);

  maxwell_field.Real().ProjectCoefficient(field_coefficient);
  maxwell_field.Imag() = maxwell_field.Real();
  const auto complex_response = maxwell_response_3d.GetMaxwellResponse(maxwell_field, 0.0);
  CHECK_THAT(complex_response.domain_correction,
             WithinRel(2.0 * real_response.domain_correction, 1.0e-10));
  CHECK_THAT(complex_response.domain_correction_fixed_flux,
             WithinRel(2.0 * real_response.domain_correction_fixed_flux, 1.0e-10));
  for (const auto &[interface, energy] : real_response.fabricated_surface_energy)
  {
    CHECK_THAT(complex_response.fabricated_surface_energy.at(interface),
               WithinRel(2.0 * energy, 1.0e-10));
  }
  for (const auto &[interface, energy] : real_response.fabricated_surface_energy_fixed_flux)
  {
    CHECK_THAT(complex_response.fabricated_surface_energy_fixed_flux.at(interface),
               WithinRel(2.0 * energy, 1.0e-10));
  }
  CHECK_THAT(complex_response.response_weighted_trace_closure_spread,
             WithinRel(real_response.response_weighted_trace_closure_spread, 1.0e-10));
  CHECK_THAT(complex_response.trace_closure_response_failure_fraction,
             WithinAbs(real_response.trace_closure_response_failure_fraction, 1.0e-12));

  mfem::VectorFunctionCoefficient rotational_field_coefficient(
      3,
      [](const mfem::Vector &x, mfem::Vector &field)
      {
        field[0] = -x[1];
        field[1] = 0.0;
        field[2] = 0.0;
      });
  maxwell_field.Real().ProjectCoefficient(rotational_field_coefficient);
  maxwell_field.Imag() = 0.0;
  const auto rotational_response =
      maxwell_response_3d.GetMaxwellResponse(maxwell_field, 0.0);
  CHECK(rotational_response.loop_residual > 0.05);
  CHECK(rotational_response.response_weighted_loop_residual > 0.0);
  CHECK(rotational_response.response_weighted_loop_residual <=
        rotational_response.loop_residual);
  CHECK(rotational_response.loop_response_failure_fraction > 0.0);
  CHECK(rotational_response.loop_response_failure_fraction <= 1.0);

  auto coupled_maxwell_config_3d = maxwell_config_3d;
  for (auto &interface :
       coupled_maxwell_config_3d["Boundaries"]["Postprocessing"]["Dielectric"])
  {
    interface["EdgeDistances"] = {0.2, 2.0, 7.0};
  }
  coupled_maxwell_config_3d["Solver"]["SurfaceResponseCorrection"]["Library"] =
      coupled_library_3d_path.string();
  IoData coupled_maxwell_iodata_3d(coupled_maxwell_config_3d, false);
  ShareCrackedAttributes(coupled_maxwell_iodata_3d, maxwell_iodata_3d);
  SpaceOperator coupled_maxwell_space_3d(coupled_maxwell_iodata_3d, maxwell_meshes_3d);
  SurfaceResponseOperator coupled_maxwell_response_3d(coupled_maxwell_iodata_3d,
                                                      coupled_maxwell_space_3d);
  const auto &maxwell_line_rule = mfem::IntRules.Get(
      mfem::Geometry::SEGMENT, 2 * coupled_maxwell_iodata_3d.solver.order);
  CHECK(coupled_maxwell_response_3d.GetPatchCount() ==
        static_cast<int>(segment_indices.size() / 2) * maxwell_line_rule.GetNPoints());

  auto interpolated_coupled_maxwell_config_3d = coupled_maxwell_config_3d;
  interpolated_coupled_maxwell_config_3d["Solver"]["SurfaceResponseCorrection"]["Library"] =
      interpolated_coupled_library_3d_path.string();
  IoData interpolated_coupled_maxwell_iodata_3d(interpolated_coupled_maxwell_config_3d,
                                                false);
  ShareCrackedAttributes(interpolated_coupled_maxwell_iodata_3d, maxwell_iodata_3d);
  SpaceOperator interpolated_coupled_maxwell_space_3d(
      interpolated_coupled_maxwell_iodata_3d, maxwell_meshes_3d);
  SurfaceResponseOperator interpolated_coupled_maxwell_response_3d(
      interpolated_coupled_maxwell_iodata_3d, interpolated_coupled_maxwell_space_3d);
  CHECK(interpolated_coupled_maxwell_response_3d.GetPatchCount() ==
        2 * coupled_maxwell_response_3d.GetPatchCount());
  CHECK(interpolated_coupled_maxwell_response_3d.GetBasisSize() ==
        2 * coupled_maxwell_response_3d.GetBasisSize());

  GridFunction coupled_maxwell_field(coupled_maxwell_space_3d.GetNDSpace(), true);
  auto cpw_potential = [](const mfem::Vector &x)
  {
    const double distance = std::abs(x[0] - 62.0);
    if (distance <= 10.0)
    {
      return 1.0;
    }
    if (distance >= 22.0)
    {
      return 0.0;
    }
    const double t = (distance - 10.0) / 12.0;
    return 1.0 - 3.0 * t * t + 2.0 * t * t * t;
  };
  mfem::FunctionCoefficient transverse_potential_coefficient(cpw_potential);
  mfem::ParGridFunction maxwell_potential(&coupled_maxwell_space_3d.GetH1Space().Get());
  maxwell_potential.ProjectCoefficient(transverse_potential_coefficient);
  mfem::ParDiscreteLinearOperator gradient(&coupled_maxwell_space_3d.GetH1Space().Get(),
                                           &coupled_maxwell_space_3d.GetNDSpace().Get());
  gradient.AddDomainInterpolator(new mfem::GradientInterpolator());
  gradient.Assemble();
  gradient.Mult(maxwell_potential, coupled_maxwell_field.Real());
  coupled_maxwell_field.Real() *= -1.0;
  coupled_maxwell_field.Imag() = 0.0;
  const auto coupled_maxwell_result =
      coupled_maxwell_response_3d.GetMaxwellResponse(coupled_maxwell_field, 0.0);
  // These matrices differ only in the appended V_B - V_A coefficient, so a nonzero
  // domain correction proves that the independent conductor state is active.
  CHECK(std::abs(coupled_maxwell_result.domain_correction) > 0.0);
  CHECK(coupled_maxwell_result.fabricated_surface_energy.size() == 3);
  CHECK(coupled_maxwell_result.loop_residual < 1.0e-10);

  GridFunction interpolated_coupled_maxwell_field(
      interpolated_coupled_maxwell_space_3d.GetNDSpace(), true);
  mfem::ParGridFunction interpolated_maxwell_potential(
      &interpolated_coupled_maxwell_space_3d.GetH1Space().Get());
  interpolated_maxwell_potential.ProjectCoefficient(transverse_potential_coefficient);
  mfem::ParDiscreteLinearOperator interpolated_gradient(
      &interpolated_coupled_maxwell_space_3d.GetH1Space().Get(),
      &interpolated_coupled_maxwell_space_3d.GetNDSpace().Get());
  interpolated_gradient.AddDomainInterpolator(new mfem::GradientInterpolator());
  interpolated_gradient.Assemble();
  interpolated_gradient.Mult(interpolated_maxwell_potential,
                             interpolated_coupled_maxwell_field.Real());
  interpolated_coupled_maxwell_field.Real() *= -1.0;
  interpolated_coupled_maxwell_field.Imag() = 0.0;
  const auto interpolated_coupled_maxwell_result =
      interpolated_coupled_maxwell_response_3d.GetMaxwellResponse(
          interpolated_coupled_maxwell_field, 0.0);
  CHECK_THAT(interpolated_coupled_maxwell_result.maximum_library_distance,
             WithinRel(4.0 / 7.0, 1.0e-12));
  CHECK(interpolated_coupled_maxwell_result.loop_residual < 1.0e-10);

  mfem::ParGridFunction transverse_potential(&coupled_laplace_3d.GetH1Space().Get());
  transverse_potential.ProjectCoefficient(transverse_potential_coefficient);
  Vector transverse_potential_true;
  transverse_potential.GetTrueDofs(transverse_potential_true);
  const auto coupled_electrostatic_result =
      coupled_response_3d.GetElectrostaticResponse(transverse_potential_true);
  mfem::ParGridFunction interpolated_transverse_potential(
      &interpolated_coupled_laplace_3d.GetH1Space().Get());
  interpolated_transverse_potential.ProjectCoefficient(transverse_potential_coefficient);
  Vector interpolated_transverse_potential_true;
  interpolated_transverse_potential.GetTrueDofs(interpolated_transverse_potential_true);
  const auto interpolated_coupled_electrostatic_result =
      interpolated_coupled_response_3d.GetElectrostaticResponse(
          interpolated_transverse_potential_true);
  CHECK_THAT(interpolated_coupled_electrostatic_result.domain_correction,
             WithinRel(coupled_electrostatic_result.domain_correction, 1.0e-12));
  CHECK_THAT(interpolated_coupled_electrostatic_result.domain_correction_fixed_flux,
             WithinRel(coupled_electrostatic_result.domain_correction_fixed_flux, 1.0e-12));
  CHECK_THAT(interpolated_coupled_maxwell_result.domain_correction,
             WithinRel(coupled_maxwell_result.domain_correction, 1.0e-12));
  CHECK_THAT(interpolated_coupled_maxwell_result.domain_correction_fixed_flux,
             WithinRel(coupled_maxwell_result.domain_correction_fixed_flux, 1.0e-12));
  for (const auto &[interface, energy] :
       coupled_electrostatic_result.fabricated_surface_energy)
  {
    CHECK_THAT(
        interpolated_coupled_electrostatic_result.fabricated_surface_energy.at(interface),
        WithinRel(energy, 1.0e-12));
    CHECK_THAT(
        interpolated_coupled_electrostatic_result.fabricated_surface_energy_fixed_flux.at(
            interface),
        WithinRel(
            coupled_electrostatic_result.fabricated_surface_energy_fixed_flux.at(interface),
            1.0e-12));
    CHECK_THAT(
        interpolated_coupled_maxwell_result.fabricated_surface_energy.at(interface),
        WithinRel(coupled_maxwell_result.fabricated_surface_energy.at(interface), 1.0e-12));
    CHECK_THAT(
        interpolated_coupled_maxwell_result.fabricated_surface_energy_fixed_flux.at(
            interface),
        WithinRel(coupled_maxwell_result.fabricated_surface_energy_fixed_flux.at(interface),
                  1.0e-12));
  }
  CHECK_THAT(coupled_maxwell_result.domain_correction,
             WithinRel(coupled_electrostatic_result.domain_correction, 1.0e-10));
  CHECK_THAT(coupled_maxwell_result.domain_correction_fixed_flux,
             WithinRel(coupled_electrostatic_result.domain_correction_fixed_flux, 1.0e-10));
  for (const auto &[interface, energy] :
       coupled_electrostatic_result.fabricated_surface_energy)
  {
    CHECK_THAT(coupled_maxwell_result.fabricated_surface_energy.at(interface),
               WithinRel(energy, 1.0e-5));
    CHECK_THAT(
        coupled_maxwell_result.fabricated_surface_energy_fixed_flux.at(interface),
        WithinRel(
            coupled_electrostatic_result.fabricated_surface_energy_fixed_flux.at(interface),
            1.0e-5));
  }

  auto parallel_cluster_maxwell_config_3d = coupled_maxwell_config_3d;
  for (auto &interface :
       parallel_cluster_maxwell_config_3d["Boundaries"]["Postprocessing"]["Dielectric"])
  {
    interface["EdgeDistances"] = {0.2, 2.0, 11.0};
  }
  parallel_cluster_maxwell_config_3d["Solver"]["SurfaceResponseCorrection"]["Library"] =
      parallel_cluster_library_3d_path.string();
  IoData parallel_cluster_maxwell_iodata_3d(parallel_cluster_maxwell_config_3d, false);
  ShareCrackedAttributes(parallel_cluster_maxwell_iodata_3d, maxwell_iodata_3d);
  SpaceOperator parallel_cluster_maxwell_space_3d(parallel_cluster_maxwell_iodata_3d,
                                                  maxwell_meshes_3d);
  SurfaceResponseOperator parallel_cluster_maxwell_response_3d(
      parallel_cluster_maxwell_iodata_3d, parallel_cluster_maxwell_space_3d);
  CHECK(parallel_cluster_maxwell_response_3d.GetPatchCount() ==
        static_cast<int>(segment_indices.size() / 4) * maxwell_line_rule.GetNPoints());
  CHECK(parallel_cluster_maxwell_response_3d.GetBasisSize() ==
        6 * parallel_cluster_maxwell_response_3d.GetPatchCount());

  auto parallel_cluster_only_maxwell_config_3d = parallel_cluster_maxwell_config_3d;
  parallel_cluster_only_maxwell_config_3d["Solver"]["SurfaceResponseCorrection"]
                                         ["Library"] =
                                             parallel_cluster_only_library_3d_path.string();
  IoData parallel_cluster_only_maxwell_iodata_3d(parallel_cluster_only_maxwell_config_3d,
                                                 false);
  ShareCrackedAttributes(parallel_cluster_only_maxwell_iodata_3d, maxwell_iodata_3d);
  SpaceOperator parallel_cluster_only_maxwell_space_3d(
      parallel_cluster_only_maxwell_iodata_3d, maxwell_meshes_3d);
  SurfaceResponseOperator parallel_cluster_only_maxwell_response_3d(
      parallel_cluster_only_maxwell_iodata_3d, parallel_cluster_only_maxwell_space_3d);
  CHECK(parallel_cluster_only_maxwell_response_3d.GetPatchCount() ==
        parallel_cluster_maxwell_response_3d.GetPatchCount());
  CHECK(parallel_cluster_only_maxwell_response_3d.GetBasisSize() ==
        parallel_cluster_maxwell_response_3d.GetBasisSize());

  GridFunction parallel_cluster_maxwell_field(
      parallel_cluster_maxwell_space_3d.GetNDSpace(), true);
  mfem::ParGridFunction parallel_cluster_maxwell_potential(
      &parallel_cluster_maxwell_space_3d.GetH1Space().Get());
  parallel_cluster_maxwell_potential.ProjectCoefficient(transverse_potential_coefficient);
  mfem::ParDiscreteLinearOperator parallel_cluster_gradient(
      &parallel_cluster_maxwell_space_3d.GetH1Space().Get(),
      &parallel_cluster_maxwell_space_3d.GetNDSpace().Get());
  parallel_cluster_gradient.AddDomainInterpolator(new mfem::GradientInterpolator());
  parallel_cluster_gradient.Assemble();
  parallel_cluster_gradient.Mult(parallel_cluster_maxwell_potential,
                                 parallel_cluster_maxwell_field.Real());
  parallel_cluster_maxwell_field.Real() *= -1.0;
  parallel_cluster_maxwell_field.Imag() = 0.0;
  const auto parallel_cluster_maxwell_result =
      parallel_cluster_maxwell_response_3d.GetMaxwellResponse(
          parallel_cluster_maxwell_field, 0.0);

  mfem::ParGridFunction parallel_cluster_electrostatic_potential(
      &parallel_cluster_laplace_3d.GetH1Space().Get());
  parallel_cluster_electrostatic_potential.ProjectCoefficient(
      transverse_potential_coefficient);
  Vector parallel_cluster_electrostatic_true;
  parallel_cluster_electrostatic_potential.GetTrueDofs(parallel_cluster_electrostatic_true);
  const auto parallel_cluster_electrostatic_result =
      parallel_cluster_response_3d.GetElectrostaticResponse(
          parallel_cluster_electrostatic_true);
  CHECK(std::abs(parallel_cluster_maxwell_result.domain_correction) > 0.0);
  CHECK(parallel_cluster_maxwell_result.loop_residual < 1.0e-10);
  CHECK_THAT(parallel_cluster_maxwell_result.domain_correction,
             WithinRel(parallel_cluster_electrostatic_result.domain_correction, 1.0e-10));
  CHECK_THAT(parallel_cluster_maxwell_result.domain_correction_fixed_flux,
             WithinRel(parallel_cluster_electrostatic_result.domain_correction_fixed_flux,
                       1.0e-10));
  for (const auto &[interface, energy] :
       parallel_cluster_electrostatic_result.fabricated_surface_energy)
  {
    CHECK_THAT(parallel_cluster_maxwell_result.fabricated_surface_energy.at(interface),
               WithinRel(energy, 1.0e-5));
    CHECK_THAT(
        parallel_cluster_maxwell_result.fabricated_surface_energy_fixed_flux.at(interface),
        WithinRel(
            parallel_cluster_electrostatic_result.fabricated_surface_energy_fixed_flux.at(
                interface),
            1.0e-4));
  }

  Vector cluster_state(parallel_cluster_maxwell_space_3d.GetNDSpace().GetTrueVSize());
  Vector cluster_probe(cluster_state.Size());
  auto *cluster_state_data = cluster_state.HostWrite();
  auto *cluster_probe_data = cluster_probe.HostWrite();
  for (int i = 0; i < cluster_state.Size(); i++)
  {
    cluster_state_data[i] = std::cos(0.19 * (i + 1 + 5 * Mpi::Rank(Mpi::World())));
    cluster_probe_data[i] = std::sin(0.27 * (i + 1 + 7 * Mpi::Rank(Mpi::World())));
  }
  cluster_state.SetSubVector(parallel_cluster_maxwell_space_3d.GetNDDbcTDofLists().back(),
                             0.0);
  cluster_probe.SetSubVector(parallel_cluster_maxwell_space_3d.GetNDDbcTDofLists().back(),
                             0.0);
  Vector cluster_correction, cluster_probe_correction;
  parallel_cluster_maxwell_response_3d.Mult(cluster_state, cluster_correction);
  parallel_cluster_maxwell_response_3d.Mult(cluster_probe, cluster_probe_correction);
  CHECK_THAT(linalg::Dot(Mpi::World(), cluster_probe, cluster_correction),
             WithinRel(linalg::Dot(Mpi::World(), cluster_state, cluster_probe_correction),
                       1.0e-10));

  // Exercise line reconstruction with a nontrivial discrete gradient. A globally
  // constant field does not detect integration errors when contour lines cross elements.
  mfem::FunctionCoefficient curved_potential_coefficient(
      [cpw_potential](const mfem::Vector &x)
      { return cpw_potential(x) + x[1] * (1.0 + 0.02 * x[0] + 0.03 * x[2]); });
  maxwell_potential.ProjectCoefficient(curved_potential_coefficient);
  gradient.Mult(maxwell_potential, coupled_maxwell_field.Real());
  coupled_maxwell_field.Real() *= -1.0;
  const auto curved_maxwell_result =
      coupled_maxwell_response_3d.GetMaxwellResponse(coupled_maxwell_field, 0.0);

  mfem::ParGridFunction electrostatic_potential(&coupled_laplace_3d.GetH1Space().Get());
  electrostatic_potential.ProjectCoefficient(curved_potential_coefficient);
  Vector electrostatic_potential_true;
  electrostatic_potential.GetTrueDofs(electrostatic_potential_true);
  const auto curved_electrostatic_result =
      coupled_response_3d.GetElectrostaticResponse(electrostatic_potential_true);
  CHECK_THAT(curved_maxwell_result.domain_correction,
             WithinRel(curved_electrostatic_result.domain_correction, 1.0e-10));
  CHECK_THAT(curved_maxwell_result.domain_correction_fixed_flux,
             WithinRel(curved_electrostatic_result.domain_correction_fixed_flux, 1.0e-10));
  for (const auto &[interface, energy] :
       curved_electrostatic_result.fabricated_surface_energy)
  {
    CHECK_THAT(curved_maxwell_result.fabricated_surface_energy.at(interface),
               WithinRel(energy, 1.0e-2));
    CHECK_THAT(
        curved_maxwell_result.fabricated_surface_energy_fixed_flux.at(interface),
        WithinRel(
            curved_electrostatic_result.fabricated_surface_energy_fixed_flux.at(interface),
            1.0e-2));
  }

#endif
}

TEST_CASE_METHOD(test::SurfaceResponseFiles, "SurfaceResponseOperator island corners",
                 "[surfaceresponseoperator][3d][corner][legacy][Serial][Parallel]")
{
#if !defined(MFEM_USE_GSLIB)
  SKIP("SurfaceResponseOperator requires MFEM_USE_GSLIB");
#else
  json island_config = IslandConfig();
  IoData concave_island_iodata(island_config, false);
  concave_island_iodata.boundaries.cracked_attributes.insert(9);
  auto island_geometry_mesh = MakeIslandMesh();
  const auto island_geometry =
      ExtractMetalEdgeGeometry(*island_geometry_mesh, concave_island_iodata.boundaries,
                               JointNoiseExtractionFor(concave_island_iodata.boundaries));
  const auto island_segments =
      GetInterfaceMetalEdgeSegmentIndices(island_geometry, 4, InterfaceDielectric::SA);
  std::set<std::size_t> island_vertices;
  double island_perimeter = 0.0;
  for (const std::size_t segment_index : island_segments)
  {
    const auto &segment = island_geometry.segments[segment_index];
    island_vertices.insert(segment.vertices.begin(), segment.vertices.end());
    const auto &p0 = island_geometry.vertices[segment.vertices[0]].coordinate;
    const auto &p1 = island_geometry.vertices[segment.vertices[1]].coordinate;
    double length_squared = 0.0;
    for (int d = 0; d < 3; d++)
    {
      length_squared += (p1[d] - p0[d]) * (p1[d] - p0[d]);
    }
    island_perimeter += std::sqrt(length_squared);
  }
  const int island_corners = static_cast<int>(
      std::count_if(island_vertices.begin(), island_vertices.end(),
                    [&](std::size_t vertex)
                    {
                      return island_geometry.vertices[vertex].physical_type ==
                             MetalEdgeVertexType::CORNER;
                    }));
  REQUIRE(island_corners == 4);
  CHECK_THAT(island_perimeter, WithinAbs(2.0, 1.0e-12));

  std::vector<std::unique_ptr<Mesh>> concave_island_meshes;
  concave_island_meshes.push_back(std::make_unique<Mesh>(MakeIslandMesh()));
  LaplaceOperator concave_island_laplace(concave_island_iodata, concave_island_meshes);
  SurfaceResponseOperator concave_island_response(concave_island_iodata,
                                                  concave_island_laplace);
  const auto &island_line_rule = mfem::IntRules.Get(mfem::Geometry::SEGMENT, 2);
  CHECK(concave_island_response.GetPatchCount() ==
        static_cast<int>(island_segments.size()) * island_line_rule.GetNPoints());
  CHECK(concave_island_response.GetBasisSize() ==
        4 * concave_island_response.GetPatchCount());
  CHECK_THAT(concave_island_response.GetPatchWeight(),
             WithinRel(island_perimeter / 0.2, 1.0e-12));

  auto convex_island_config = island_config;
  convex_island_config["Solver"]["Electrostatic"]["ResponseCorrection"]["Library"] =
      convex_library_3d_path.string();
  IoData convex_island_iodata(convex_island_config, false);
  convex_island_iodata.boundaries.cracked_attributes.insert(9);
  std::vector<std::unique_ptr<Mesh>> convex_island_meshes;
  convex_island_meshes.push_back(std::make_unique<Mesh>(MakeIslandMesh()));
  LaplaceOperator convex_island_laplace(convex_island_iodata, convex_island_meshes);
  SurfaceResponseOperator convex_island_response(convex_island_iodata,
                                                 convex_island_laplace);
  const int removed_straight_patches = 2 * island_corners * island_line_rule.GetNPoints();
  CHECK(convex_island_response.GetPatchCount() == concave_island_response.GetPatchCount() -
                                                      removed_straight_patches +
                                                      island_corners);
  CHECK(convex_island_response.GetBasisSize() ==
        4 * (convex_island_response.GetPatchCount() - island_corners) +
            12 * island_corners);
  const double expected_convex_weight =
      (island_perimeter - 2.0 * island_corners * 0.2) / 0.2 + island_corners;
  CHECK_THAT(convex_island_response.GetPatchWeight(),
             WithinRel(expected_convex_weight, 1.0e-12));

  // Regression (lane J double count, decision 112(a)): a spatial (3D box) coupon adds the
  // energy within R of its edges (`Q_ij (J)` at the matching radius), not the whole-box
  // `Q_total_ij (J)` — the device keeps its own raw energy beyond R inside the box. The
  // island's convex corners are 3D box coupons: inflating their Q_total 3x with Q_ij
  // unchanged must leave the fabricated surface energy unchanged (the old assembly read
  // Q_total and grew), scaling Q_ij as well must change it (the corners contribute), and a
  // corner file without the within-R column or with its rows at another radius is refused.
  {
    mfem::ParGridFunction island_potential(&convex_island_laplace.GetH1Space().Get());
    mfem::FunctionCoefficient island_potential_coefficient(
        [](const mfem::Vector &x) { return x[1] * (1.0 + 0.3 * x[0] - 0.2 * x[2]); });
    island_potential.ProjectCoefficient(island_potential_coefficient);
    Vector island_potential_true;
    island_potential.GetTrueDofs(island_potential_true);
    const auto convex_island_result =
        convex_island_response.GetElectrostaticResponse(island_potential_true);
    REQUIRE(convex_island_result.fabricated_surface_energy.at(4) > 0.0);
    auto VariantConfig = [&](const auto &library_path)
    {
      auto config = convex_island_config;
      config["Solver"]["Electrostatic"]["ResponseCorrection"]["Library"] =
          library_path.string();
      return config;
    };
    {
      IoData inflated_iodata(VariantConfig(inflated_box_convex_library_3d_path), false);
      inflated_iodata.boundaries.cracked_attributes.insert(9);
      SurfaceResponseOperator inflated_response(inflated_iodata, convex_island_laplace);
      const auto inflated_result =
          inflated_response.GetElectrostaticResponse(island_potential_true);
      CHECK_THAT(inflated_result.fabricated_surface_energy.at(4),
                 WithinRel(convex_island_result.fabricated_surface_energy.at(4), 1.0e-12));
      CHECK_THAT(inflated_result.fabricated_surface_energy_fixed_flux.at(4),
                 WithinRel(convex_island_result.fabricated_surface_energy_fixed_flux.at(4),
                           1.0e-12));
    }
    {
      IoData scaled_iodata(VariantConfig(scaled_convex_library_3d_path), false);
      scaled_iodata.boundaries.cracked_attributes.insert(9);
      SurfaceResponseOperator scaled_response(scaled_iodata, convex_island_laplace);
      const auto scaled_result =
          scaled_response.GetElectrostaticResponse(island_potential_true);
      CHECK(scaled_result.fabricated_surface_energy.at(4) >
            1.01 * convex_island_result.fabricated_surface_energy.at(4));
    }
    {
      IoData legacy_iodata(VariantConfig(legacy_compact_convex_library_3d_path), false);
      legacy_iodata.boundaries.cracked_attributes.insert(9);
      CHECK_THROWS(SurfaceResponseOperator(legacy_iodata, convex_island_laplace));
    }
    {
      IoData other_radius_iodata(VariantConfig(other_radius_convex_library_3d_path), false);
      other_radius_iodata.boundaries.cracked_attributes.insert(9);
      CHECK_THROWS(SurfaceResponseOperator(other_radius_iodata, convex_island_laplace));
    }
  }

  auto touching_geometry_mesh = MakeTouchingIslandMesh();
  const auto touching_geometry =
      ExtractMetalEdgeGeometry(*touching_geometry_mesh, convex_island_iodata.boundaries,
                               JointNoiseExtractionFor(convex_island_iodata.boundaries));
  const auto touching_segments =
      GetInterfaceMetalEdgeSegmentIndices(touching_geometry, 4, InterfaceDielectric::SA);
  std::set<std::size_t> touching_vertices;
  for (const std::size_t segment_index : touching_segments)
  {
    const auto &segment = touching_geometry.segments[segment_index];
    touching_vertices.insert(segment.vertices.begin(), segment.vertices.end());
  }
  const int touching_junctions = static_cast<int>(
      std::count_if(touching_vertices.begin(), touching_vertices.end(),
                    [&](std::size_t vertex)
                    {
                      return touching_geometry.vertices[vertex].physical_type ==
                             MetalEdgeVertexType::JUNCTION;
                    }));
  REQUIRE(touching_junctions == 1);

  std::vector<std::unique_ptr<Mesh>> touching_island_meshes;
  touching_island_meshes.push_back(std::make_unique<Mesh>(MakeTouchingIslandMesh()));
  LaplaceOperator touching_island_laplace(convex_island_iodata, touching_island_meshes);
  SurfaceResponseOperator touching_island_response(convex_island_iodata,
                                                   touching_island_laplace);
  auto junction_island_config = convex_island_config;
  junction_island_config["Solver"]["Electrostatic"]["ResponseCorrection"]["Library"] =
      junction_library_3d_path.string();
  IoData junction_island_iodata(junction_island_config, false);
  junction_island_iodata.boundaries.cracked_attributes.insert(9);
  std::vector<std::unique_ptr<Mesh>> junction_island_meshes;
  junction_island_meshes.push_back(std::make_unique<Mesh>(MakeTouchingIslandMesh()));
  LaplaceOperator junction_island_laplace(junction_island_iodata, junction_island_meshes);
  SurfaceResponseOperator junction_island_response(junction_island_iodata,
                                                   junction_island_laplace);
  // Phase 3 of the identification fix: conductor identity is the edge-connected metal
  // component, so the two squares touching at one vertex are two conductors (a point
  // contact carries no galvanic connection); the legacy junction patch requires one
  // conductor on every arm and is not built. The vertex is reported as a PointContact.
  CHECK(junction_island_response.GetPatchCount() ==
        touching_island_response.GetPatchCount());
  CHECK(junction_island_response.GetBasisSize() == touching_island_response.GetBasisSize());
  CHECK_THAT(junction_island_response.GetPatchWeight(),
             WithinRel(touching_island_response.GetPatchWeight(), 1.0e-12));

#endif
}

TEST_CASE_METHOD(test::SurfaceResponseFiles, "SurfaceResponseOperator cap-interior hats",
                 "[surfaceresponseoperator][3d][spatial][cache][Serial][Parallel]")
{
#if !defined(MFEM_USE_GSLIB)
  SKIP("SurfaceResponseOperator requires MFEM_USE_GSLIB");
#else
  json island_config = IslandConfig();
  // Cap-interior hats (decision 112(b), library key InteriorTraceCount): the trailing basis
  // points of a spatial model lie on no contour and are ordinary trace coefficients. On the
  // offset corner pair (electrostatic, Collocated): the cap-hat library whose open paths
  // partition BasisPoints - InteriorTraceCount loads and, with zero hat rows, reproduces
  // the ring-only model's response exactly; a hat energy enters the fabricated surface
  // energy (and not the domain defect, equal in both coupons); the geometry cache carries
  // the key (round trip exact) and a cache without it fails the partition check; refused:
  // open paths summing to BasisPoints, the key without a TraceMesh, the key on a
  // translational model.
  {
    auto cap_hat_config = island_config;
    cap_hat_config["Boundaries"]["Terminal"] = {{{"Index", 1}, {"Attributes", {9}}},
                                                {{"Index", 2}, {"Attributes", {10}}}};
    cap_hat_config["Boundaries"]["Postprocessing"]["Dielectric"][0]["Attributes"] = {9, 10};
    auto CapHatIoData = [&](const auto &library_path)
    {
      auto config = cap_hat_config;
      config["Solver"]["Electrostatic"]["ResponseCorrection"]["Library"] =
          library_path.string();
      IoData iodata(config, false);
      iodata.boundaries.cracked_attributes.insert(9);
      iodata.boundaries.cracked_attributes.insert(10);
      return iodata;
    };
    std::vector<std::unique_ptr<Mesh>> cap_hat_meshes;
    cap_hat_meshes.push_back(std::make_unique<Mesh>(MakeOffsetCornerPairMesh()));
    IoData ring_only_iodata = CapHatIoData(spatial_cluster_library_3d_path);
    LaplaceOperator cap_hat_laplace(ring_only_iodata, cap_hat_meshes);
    SurfaceResponseOperator ring_only_response(ring_only_iodata, cap_hat_laplace);
    // One spatial edge-cluster patch (the two facing corners) among the isolated-edge and
    // convex-corner patches of the two islands: the hats add two coefficients.
    const int ring_only_patches = ring_only_response.GetPatchCount();
    const int ring_only_basis = ring_only_response.GetBasisSize();
    REQUIRE(ring_only_patches > 0);
    mfem::ParGridFunction cap_hat_potential(&cap_hat_laplace.GetH1Space().Get());
    mfem::FunctionCoefficient cap_hat_potential_coefficient(
        [](const mfem::Vector &x) { return x[1] * (1.0 + 0.3 * x[0] - 0.2 * x[2]); });
    cap_hat_potential.ProjectCoefficient(cap_hat_potential_coefficient);
    Vector cap_hat_potential_true;
    cap_hat_potential.GetTrueDofs(cap_hat_potential_true);
    const auto ring_only_result =
        ring_only_response.GetElectrostaticResponse(cap_hat_potential_true);
    REQUIRE(ring_only_result.fabricated_surface_energy.at(4) > 0.0);
    REQUIRE(ring_only_result.domain_correction != 0.0);

    IoData cap_hat_iodata = CapHatIoData(cap_hat_spatial_cluster_library_3d_path);
    SurfaceResponseOperator cap_hat_response(cap_hat_iodata, cap_hat_laplace);
    CHECK(cap_hat_response.GetPatchCount() == ring_only_patches);
    CHECK(cap_hat_response.GetBasisSize() == ring_only_basis + 2);
    const auto cap_hat_result =
        cap_hat_response.GetElectrostaticResponse(cap_hat_potential_true);
    CHECK_THAT(cap_hat_result.domain_correction,
               WithinRel(ring_only_result.domain_correction, 1.0e-12));
    CHECK_THAT(cap_hat_result.domain_correction_fixed_flux,
               WithinRel(ring_only_result.domain_correction_fixed_flux, 1.0e-9));
    CHECK_THAT(cap_hat_result.fabricated_surface_energy.at(4),
               WithinRel(ring_only_result.fabricated_surface_energy.at(4), 1.0e-12));
    CHECK_THAT(
        cap_hat_result.fabricated_surface_energy_fixed_flux.at(4),
        WithinRel(ring_only_result.fabricated_surface_energy_fixed_flux.at(4), 1.0e-9));

    IoData cap_hat_loaded_iodata =
        CapHatIoData(cap_hat_loaded_spatial_cluster_library_3d_path);
    SurfaceResponseOperator cap_hat_loaded_response(cap_hat_loaded_iodata, cap_hat_laplace);
    CHECK(cap_hat_loaded_response.GetBasisSize() == ring_only_basis + 2);
    const auto cap_hat_loaded_result =
        cap_hat_loaded_response.GetElectrostaticResponse(cap_hat_potential_true);
    CHECK_THAT(cap_hat_loaded_result.domain_correction,
               WithinRel(ring_only_result.domain_correction, 1.0e-12));
    // The one cluster patch among ~60 patches of the two islands: its two hats at a
    // 1e-12 J diagonal raise the interface's fabricated surface energy by ~0.26 %.
    CHECK(cap_hat_loaded_result.fabricated_surface_energy.at(4) >
          1.001 * ring_only_result.fabricated_surface_energy.at(4));

    {
      const auto cache_path = temp.temp_dir / "response-geometry-cap-hats.json";
      test::GeometryCacheEnvGuard cache_env(cache_path.string(), true);
      SurfaceResponseOperator written_response(cap_hat_loaded_iodata, cap_hat_laplace);
      Mpi::Barrier(Mpi::World());
      cache_env.DisableWrite();
      SurfaceResponseOperator reloaded_response(cap_hat_loaded_iodata, cap_hat_laplace);
      CHECK(reloaded_response.GetPatchCount() == ring_only_patches);
      CHECK(reloaded_response.GetBasisSize() == ring_only_basis + 2);
      const auto reloaded_result =
          reloaded_response.GetElectrostaticResponse(cap_hat_potential_true);
      CHECK_THAT(reloaded_result.domain_correction,
                 WithinRel(cap_hat_loaded_result.domain_correction, 1.0e-12));
      CHECK_THAT(reloaded_result.fabricated_surface_energy.at(4),
                 WithinRel(cap_hat_loaded_result.fabricated_surface_energy.at(4), 1.0e-12));
      std::ifstream cache_input(cache_path);
      REQUIRE(cache_input);
      json cache = json::parse(cache_input);
      cache_input.close();
      CHECK(cache["Version"] == 14);
      int cap_hat_models = 0;
      for (auto &model : cache["Models"])
      {
        if (model["Name"] == "offset-corner-pair-cap-hats")
        {
          cap_hat_models++;
          CHECK(model["InteriorTraceCount"] == 2);
          CHECK(model["SpatialBasis"] == true);
          model["InteriorTraceCount"] = 0;
        }
      }
      CHECK(cap_hat_models == 1);
      // Without the key the open paths partition four of six BasisPoints: refused.
      const auto stale_path = temp.temp_dir / "response-geometry-cap-hats-stale.json";
      if (Mpi::Root(Mpi::World()))
      {
        std::ofstream stale(stale_path);
        stale << cache.dump(2) << "\n";
      }
      Mpi::Barrier(Mpi::World());
      test::GeometryCacheEnvGuard stale_cache_env(stale_path.string(), false);
      CHECK_THROWS_WITH(SurfaceResponseOperator(cap_hat_loaded_iodata, cap_hat_laplace),
                        Catch::Matchers::ContainsSubstring("do not partition the contour"));
    }

    IoData cap_hat_full_partition_iodata =
        CapHatIoData(cap_hat_full_partition_library_3d_path);
    CHECK_THROWS_WITH(
        SurfaceResponseOperator(cap_hat_full_partition_iodata, cap_hat_laplace),
        Catch::Matchers::ContainsSubstring("OpenContourPaths contain an invalid"));
    IoData cap_hat_without_trace_mesh_iodata =
        CapHatIoData(cap_hat_without_trace_mesh_library_3d_path);
    CHECK_THROWS_WITH(
        SurfaceResponseOperator(cap_hat_without_trace_mesh_iodata, cap_hat_laplace),
        Catch::Matchers::ContainsSubstring("InteriorTraceCount requires"));
    IoData cap_hat_translational_iodata =
        CapHatIoData(cap_hat_translational_library_3d_path);
    CHECK_THROWS_WITH(
        SurfaceResponseOperator(cap_hat_translational_iodata, cap_hat_laplace),
        Catch::Matchers::ContainsSubstring("InteriorTraceCount requires"));
  }

#endif
}

TEST_CASE_METHOD(test::SurfaceResponseFiles,
                 "SurfaceResponseOperator features patch construction",
                 "[surfaceresponseoperator][3d][features][placement][Serial][Parallel]")
{
#if !defined(MFEM_USE_GSLIB)
  SKIP("SurfaceResponseOperator requires MFEM_USE_GSLIB");
#else
  json island_config = IslandConfig();
  // The island's perimeter and corner count (asserted by the island corners case).
  IoData concave_island_iodata(island_config, false);
  concave_island_iodata.boundaries.cracked_attributes.insert(9);
  const auto island_edges =
      SummarizeInterfaceEdges(*MakeIslandMesh(), concave_island_iodata.boundaries, 4);
  const int island_corners = island_edges.corners;
  const double island_perimeter = island_edges.length;
  // Features-driven patch construction (the default; SURFACE-RESPONSE-IDENTIFICATION.md
  // (e)): the preflight is the patch dry run. With the legacy convex library every feature
  // is Missing (its models map three interface types, the island's features carry SA only:
  // key-based matching, no parametric tolerance), so nothing is patched and the solve path
  // aborts under UnmatchedPolicy = Error. A signature-keyed library built from the
  // manifest's own features patches every feature exactly once: the corner windows as one
  // patch each in the feature frame, the isolated edges as one quadrature per portion.
  {
    auto features_island_config = island_config;
    auto &features_correction =
        features_island_config["Solver"]["Electrostatic"]["ResponseCorrection"];
    features_correction.erase("PatchConstruction");
    features_correction["Library"] = convex_library_3d_path.string();
    IoData features_island_iodata(features_island_config, false);
    features_island_iodata.boundaries.cracked_attributes.insert(9);
    REQUIRE(features_island_iodata.solver.electrostatic.response_correction
                ->patch_construction ==
            config::ElectrostaticSolverData::ResponseCorrectionData::PatchConstruction::
                FEATURES);
    auto features_mesh = MakeIslandMesh();
    Mesh features_island_mesh(std::move(features_mesh));
    const auto features_manifest_path =
        temp.temp_dir / "surface-response-requirements-features-island.json";
    const auto features_patches_path = temp.temp_dir / "surface-response-patches.csv";
    fs::remove(features_patches_path);
    WriteSurfaceResponseRequirements(features_island_iodata, features_island_mesh,
                                     features_manifest_path.string());
    Mpi::Barrier(Mpi::World());
    std::ifstream features_manifest_input(features_manifest_path);
    REQUIRE(features_manifest_input);
    const json features_manifest = json::parse(features_manifest_input);
    const auto &features = features_manifest["Identification"]["Features"];
    REQUIRE(features.size() == 8);  // 4 convex corners + 4 isolated edges
    for (const auto &feature : features)
    {
      CHECK(feature["Match"]["Status"] == "Missing");
      if (feature["Type"] == "ConvexCorner")
      {
        // The corner frame: x = the first arm, y = the second (counterclockwise about the
        // process normal), z = the process normal (design (b) 4).
        const auto axes =
            feature["Frame"]["Axes"].get<std::array<std::array<double, 3>, 3>>();
        CHECK_THAT(std::hypot(axes[0][0], axes[0][1], axes[0][2]), WithinAbs(1.0, 1.0e-12));
        CHECK_THAT(std::hypot(axes[1][0], axes[1][1], axes[1][2]), WithinAbs(1.0, 1.0e-12));
        CHECK_THAT(axes[2][1], WithinAbs(1.0, 1.0e-12));
        CHECK_THAT(axes[0][0] * axes[1][0] + axes[0][1] * axes[1][1] +
                       axes[0][2] * axes[1][2],
                   WithinAbs(0.0, 1.0e-12));
        // Right-handed: x cross y = z.
        CHECK_THAT(axes[0][2] * axes[1][0] - axes[0][0] * axes[1][2],
                   WithinAbs(1.0, 1.0e-12));
      }
    }
    auto ReadPatches = [](const fs::path &path)
    {
      std::ifstream input(path);
      REQUIRE(input);
      std::string line;
      std::getline(input, line);  // header
      std::vector<std::vector<std::string>> rows;
      while (std::getline(input, line))
      {
        std::vector<std::string> fields;
        std::stringstream stream(line);
        std::string field;
        while (std::getline(stream, field, ','))
        {
          fields.push_back(field);
        }
        rows.push_back(std::move(fields));
      }
      return rows;
    };
    CHECK(ReadPatches(features_patches_path).empty());
    {
      std::vector<std::unique_ptr<Mesh>> meshes;
      meshes.push_back(std::make_unique<Mesh>(MakeIslandMesh()));
      LaplaceOperator laplace(features_island_iodata, meshes);
      CHECK_THROWS_WITH(SurfaceResponseOperator(features_island_iodata, laplace),
                        ContainsSubstring("8 identified feature(s) have no library model"));
    }

    // Signature-keyed library from the manifest: the legacy corner matrices (12 knots, 3
    // contour groups) for the corners, the straight-edge matrices for the isolated edges.
    const auto signature_library_path =
        temp.temp_dir / "fabrication-process-signature-3d.json";
    if (Mpi::Root(Mpi::World()))
    {
      json signature_library = {{"Version", 2},
                                {"Name", "unit-test-signature-3d"},
                                {"MatchingRadius", 0.2},
                                {"CouponDepth", 0.2},
                                {"Models", json::array()}};
      std::set<std::string> hashes;
      for (const auto &feature : features)
      {
        if (!hashes.insert(feature["Hash"].get<std::string>()).second)
        {
          continue;
        }
        json model = {{"Name", feature["Type"].get<std::string>() + "-" +
                                   feature["Hash"].get<std::string>().substr(0, 12)},
                      {"Topology", feature["Type"]},
                      {"Signature", feature["Signature"]}};
        if (feature["Type"] == "ConvexCorner")
        {
          model["Angle"] = feature["Signature"]["AngleDegrees"];
          model["FabricatedMatrix"] = corner_fabricated_path.string();
          model["ThinMatrix"] = corner_thin_path.string();
          model["BasisPoints"] = corner_points_path.string();
          model["ContourGroups"] = {4, 4, 4};
        }
        else
        {
          model["FabricatedMatrix"] = fabricated_path.string();
          model["ThinMatrix"] = thin_path.string();
          model["BasisPoints"] = points_path.string();
        }
        signature_library["Models"].push_back(model);
      }
      std::ofstream output(signature_library_path);
      output << signature_library.dump(2) << "\n";
    }
    Mpi::Barrier(Mpi::World());
    features_correction["Library"] = signature_library_path.string();
    IoData signature_island_iodata(features_island_config, false);
    signature_island_iodata.boundaries.cracked_attributes.insert(9);
    const auto signature_manifest_path =
        temp.temp_dir / "surface-response-requirements-signature-island.json";
    WriteSurfaceResponseRequirements(signature_island_iodata, features_island_mesh,
                                     signature_manifest_path.string());
    Mpi::Barrier(Mpi::World());
    std::ifstream signature_manifest_input(signature_manifest_path);
    REQUIRE(signature_manifest_input);
    const json signature_manifest = json::parse(signature_manifest_input);
    const auto &signature_identification = signature_manifest["Identification"];
    std::map<int, json> features_by_id;
    for (const auto &feature : signature_identification["Features"])
    {
      CHECK(feature["Match"]["Status"] == "Matched");
      features_by_id.emplace(feature["Id"].get<int>(), feature);
    }
    const auto rows = ReadPatches(features_patches_path);
    // Columns: Patch, Feature, Topology, Model, ModelIndex, Weight, ModelWeight,
    // QuadratureWeight, SideFactor, CouponDepth, Segment, S0, S1, Origin(3), AxisU(3),
    // AxisV(3), AxisW(3), StripBegin, StripEnd.
    std::set<int> patched_features;
    std::map<int, std::set<std::tuple<int, double, double>>> intervals_by_feature;
    std::map<std::tuple<int, int, double, double>, double> quadrature_sums;
    int corner_patches = 0;
    double total_weight = 0.0, isolated_length = 0.0;
    std::set<std::array<double, 3>> corner_origins;
    for (const auto &row : rows)
    {
      REQUIRE(row.size() == 27);
      const int feature = std::stoi(row[1]);
      patched_features.insert(feature);
      const double weight = std::stod(row[5]);
      total_weight += weight;
      const int segment = std::stoi(row[10]);
      const double s0 = std::stod(row[11]), s1 = std::stod(row[12]);
      if (features_by_id.at(feature)["Type"] == "ConvexCorner")
      {
        CHECK(segment == -1);
        CHECK_THAT(weight, WithinAbs(1.0, 1.0e-12));
        corner_patches++;
        corner_origins.insert({std::stod(row[13]), std::stod(row[14]), std::stod(row[15])});
        // The patch frame is the feature frame: a right-handed frame with z = the process
        // normal (0, 1, 0), x and y along the island edges.
        CHECK_THAT(std::stod(row[23]), WithinAbs(1.0, 1.0e-12));
        CHECK_THAT(std::abs(std::stod(row[16])) + std::abs(std::stod(row[18])),
                   WithinAbs(1.0, 1.0e-12));
        CHECK_THAT(std::stod(row[16]) * std::stod(row[21]) -
                       std::stod(row[18]) * std::stod(row[19]),
                   WithinAbs(-1.0, 1.0e-12));
      }
      else
      {
        CHECK(segment >= 0);
        CHECK(std::stod(row[8]) == 1.0);  // side factor of a single edge
        CHECK_THAT(std::stod(row[9]), WithinRel(0.2, 1.0e-12));  // coupon depth
        CHECK_THAT(weight, WithinRel(std::stod(row[7]) * (s1 - s0) / 0.2, 1.0e-12));
        if (intervals_by_feature[feature].insert({segment, s0, s1}).second)
        {
          isolated_length += s1 - s0;
        }
        quadrature_sums[{feature, segment, s0, s1}] += std::stod(row[7]);
      }
    }
    CHECK(patched_features.size() == features_by_id.size());
    CHECK(corner_patches == island_corners);
    CHECK(corner_origins.size() == static_cast<std::size_t>(island_corners));
    for (const auto &[key, sum] : quadrature_sums)
    {
      CHECK_THAT(sum, WithinAbs(1.0, 1.0e-12));
    }
    // Every portion of every isolated-edge feature is exactly one quadrature interval.
    for (const auto &[id, feature] : features_by_id)
    {
      if (feature["Type"] != "IsolatedEdge")
      {
        continue;
      }
      std::set<std::tuple<int, double, double>> portions;
      for (const auto &portion : feature["Portions"])
      {
        portions.insert(
            {portion[0].get<int>(), portion[1].get<double>(), portion[2].get<double>()});
      }
      CHECK(portions == intervals_by_feature.at(id));
    }
    // The whole perimeter minus the corner windows (R along each of the 8 arms) is
    // isolated.
    CHECK_THAT(isolated_length,
               WithinAbs(island_perimeter - 2.0 * island_corners * 0.2, 1.0e-9));
    CHECK_THAT(
        isolated_length,
        WithinAbs(signature_identification["Totals"]["AssignedLength"].get<double>() -
                      2.0 * island_corners * 0.2,
                  1.0e-9));
    CHECK_THAT(total_weight, WithinRel(isolated_length / 0.2 + island_corners, 1.0e-12));
    // The solve path builds the same patches (no field solve: the operator only assembles
    // the local coupon responses).
    std::vector<std::unique_ptr<Mesh>> signature_meshes;
    signature_meshes.push_back(std::make_unique<Mesh>(MakeIslandMesh()));
    LaplaceOperator signature_laplace(signature_island_iodata, signature_meshes);
    SurfaceResponseOperator signature_response(signature_island_iodata, signature_laplace);
    CHECK(signature_response.GetPatchCount() == static_cast<int>(rows.size()));
    CHECK_THAT(signature_response.GetPatchWeight(), WithinRel(total_weight, 1.0e-12));
    CHECK(signature_response.GetBasisSize() ==
          4 * (static_cast<int>(rows.size()) - island_corners) + 12 * island_corners);
  }

#endif
}

TEST_CASE_METHOD(test::SurfaceResponseFiles, "SurfaceResponseOperator rounded islands",
                 "[surfaceresponseoperator][3d][curvature][rounded][Serial][Parallel]")
{
#if !defined(MFEM_USE_GSLIB)
  SKIP("SurfaceResponseOperator requires MFEM_USE_GSLIB");
#else
  json island_config = IslandConfig();
  const auto &island_line_rule = mfem::IntRules.Get(mfem::Geometry::SEGMENT, 2);
  json rounded_island_config = RoundedIslandConfig();
  IoData rounded_island_iodata(rounded_island_config, false);
  rounded_island_iodata.boundaries.cracked_attributes.insert(9);
  auto rounded_geometry_mesh = MakeIslandMesh(true);
  const auto rounded_geometry =
      ExtractMetalEdgeGeometry(*rounded_geometry_mesh, rounded_island_iodata.boundaries,
                               JointNoiseExtractionFor(rounded_island_iodata.boundaries));
  const auto rounded_segments =
      GetInterfaceMetalEdgeSegmentIndices(rounded_geometry, 4, InterfaceDielectric::SA);
  std::set<std::size_t> rounded_vertices;
  for (const std::size_t segment_index : rounded_segments)
  {
    const auto &segment = rounded_geometry.segments[segment_index];
    rounded_vertices.insert(segment.vertices.begin(), segment.vertices.end());
  }
  // Under the geometric joint noise rule (kJointNoiseSagittaOverRadius = 0.05, USER
  // decision 121 (B)) the sampled fillet joints (8-45 deg per chord on chords far below R:
  // implied sagitta (c / 2) tan(t / 4) below 0.05 R) are REGULAR joints of the perimeter
  // extraction (the 1 deg angular threshold of 117(4) made them corners); the
  // identification reads the non-collinear joints regardless and its arc rule fits the
  // fillets. The legacy per-group classifier below (PatchConstruction "Legacy", comparison
  // only) re-extracts the perimeter at its own 30 deg corner class
  // (kLegacyCornerTurnToleranceDegrees) so that its rounded-run rule reads the fillets as
  // REGULAR runs.
  CHECK(std::count_if(rounded_vertices.begin(), rounded_vertices.end(),
                      [&](std::size_t vertex)
                      {
                        return rounded_geometry.vertices[vertex].physical_type ==
                               MetalEdgeVertexType::CORNER;
                      }) == 0);

  std::vector<std::unique_ptr<Mesh>> rounded_island_meshes;
  rounded_island_meshes.push_back(std::make_unique<Mesh>(MakeIslandMesh(true)));
  LaplaceOperator rounded_island_laplace(rounded_island_iodata, rounded_island_meshes);
  SurfaceResponseOperator rounded_island_response(rounded_island_iodata,
                                                  rounded_island_laplace);
  constexpr int rounded_corner_count = 4;
  constexpr int remaining_straight_intervals = 8;
  CHECK(rounded_island_response.GetPatchCount() ==
        remaining_straight_intervals * island_line_rule.GetNPoints() +
            rounded_corner_count);
  CHECK(rounded_island_response.GetBasisSize() ==
        4 * remaining_straight_intervals * island_line_rule.GetNPoints() +
            12 * rounded_corner_count);
  CHECK_THAT(rounded_island_response.GetPatchWeight(), WithinRel(6.0, 1.0e-12));

  auto rounded_concave_island_config = island_config;
  rounded_concave_island_config["Solver"]["Electrostatic"]["ResponseCorrection"]
                               ["Library"] = rounded_concave_library_3d_path.string();
  rounded_concave_island_config["Boundaries"]["Ground"]["Attributes"].push_back(9);
  rounded_concave_island_config["Boundaries"].erase("Terminal");
  rounded_concave_island_config["Boundaries"]["Postprocessing"]["Dielectric"][0]
                               ["EdgeExcludeAttributes"] = {1, 2, 3, 4, 5, 6};
  IoData rounded_concave_island_iodata(rounded_concave_island_config, false);
  rounded_concave_island_iodata.boundaries.cracked_attributes.insert(9);
  std::vector<std::unique_ptr<Mesh>> rounded_concave_island_meshes;
  rounded_concave_island_meshes.push_back(
      std::make_unique<Mesh>(MakeIslandMesh(true, false, true)));
  LaplaceOperator rounded_concave_island_laplace(rounded_concave_island_iodata,
                                                 rounded_concave_island_meshes);
  SurfaceResponseOperator rounded_concave_island_response(rounded_concave_island_iodata,
                                                          rounded_concave_island_laplace);
  CHECK(rounded_concave_island_response.GetPatchCount() ==
        remaining_straight_intervals * island_line_rule.GetNPoints() +
            rounded_corner_count);
  CHECK(rounded_concave_island_response.GetBasisSize() ==
        4 * remaining_straight_intervals * island_line_rule.GetNPoints() +
            12 * rounded_corner_count);
  CHECK_THAT(rounded_concave_island_response.GetPatchWeight(), WithinRel(6.0, 1.0e-12));

  auto interpolated_rounded_island_config = rounded_island_config;
  interpolated_rounded_island_config["Solver"]["Electrostatic"]["ResponseCorrection"]
                                    ["Library"] =
                                        interpolated_rounded_library_3d_path.string();
  IoData interpolated_rounded_island_iodata(interpolated_rounded_island_config, false);
  interpolated_rounded_island_iodata.boundaries.cracked_attributes.insert(9);
  std::vector<std::unique_ptr<Mesh>> interpolated_rounded_island_meshes;
  interpolated_rounded_island_meshes.push_back(
      std::make_unique<Mesh>(MakeIslandMesh(true)));
  LaplaceOperator interpolated_rounded_island_laplace(interpolated_rounded_island_iodata,
                                                      interpolated_rounded_island_meshes);
  SurfaceResponseOperator interpolated_rounded_island_response(
      interpolated_rounded_island_iodata, interpolated_rounded_island_laplace);
  CHECK(interpolated_rounded_island_response.GetPatchCount() ==
        remaining_straight_intervals * island_line_rule.GetNPoints() +
            2 * rounded_corner_count);
  CHECK(interpolated_rounded_island_response.GetBasisSize() ==
        4 * remaining_straight_intervals * island_line_rule.GetNPoints() +
            2 * 12 * rounded_corner_count);
  CHECK_THAT(interpolated_rounded_island_response.GetPatchWeight(),
             WithinRel(6.0, 1.0e-12));
  auto unqualified_interpolated_rounded_island_config = rounded_island_config;
  unqualified_interpolated_rounded_island_config
      ["Solver"]["Electrostatic"]["ResponseCorrection"]["Library"] =
          unqualified_interpolated_rounded_library_3d_path.string();
  IoData unqualified_interpolated_rounded_island_iodata(
      unqualified_interpolated_rounded_island_config, false);
  unqualified_interpolated_rounded_island_iodata.boundaries.cracked_attributes.insert(9);
  std::vector<std::unique_ptr<Mesh>> unqualified_interpolated_rounded_island_meshes;
  unqualified_interpolated_rounded_island_meshes.push_back(
      std::make_unique<Mesh>(MakeIslandMesh(true)));
  LaplaceOperator unqualified_interpolated_rounded_island_laplace(
      unqualified_interpolated_rounded_island_iodata,
      unqualified_interpolated_rounded_island_meshes);
  CHECK_THROWS_WITH(
      SurfaceResponseOperator(unqualified_interpolated_rounded_island_iodata,
                              unqualified_interpolated_rounded_island_laplace),
      Catch::Matchers::ContainsSubstring(
          "Automatic fabrication-process response matching failed"));
  Vector interpolation_probe(rounded_island_response.Height());
  auto *interpolation_probe_data = interpolation_probe.HostWrite();
  for (int i = 0; i < interpolation_probe.Size(); i++)
  {
    interpolation_probe_data[i] = std::cos(0.17 * (i + 1 + 3 * Mpi::Rank(Mpi::World())));
  }
  Vector rounded_correction, interpolated_rounded_correction;
  rounded_island_response.Mult(interpolation_probe, rounded_correction);
  interpolated_rounded_island_response.Mult(interpolation_probe,
                                            interpolated_rounded_correction);
  interpolated_rounded_correction.Add(-1.0, rounded_correction);
  CHECK(linalg::Norml2(Mpi::World(), interpolated_rounded_correction) <=
        1.0e-12 * std::max(linalg::Norml2(Mpi::World(), rounded_correction), 1.0e-300));

#endif
}

TEST_CASE_METHOD(test::SurfaceResponseFiles,
                 "SurfaceResponseOperator high-order rounded island",
                 "[surfaceresponseoperator][3d][curvature][highorder][Serial][Parallel]")
{
#if !defined(MFEM_USE_GSLIB)
  SKIP("SurfaceResponseOperator requires MFEM_USE_GSLIB");
#else
  json island_config = IslandConfig();
  json rounded_island_config = RoundedIslandConfig();
  IoData rounded_island_iodata(rounded_island_config, false);
  rounded_island_iodata.boundaries.cracked_attributes.insert(9);
  // The linear rounded island (the rounded islands case asserts its patches).
  std::vector<std::unique_ptr<Mesh>> rounded_island_meshes;
  rounded_island_meshes.push_back(std::make_unique<Mesh>(MakeIslandMesh(true)));
  LaplaceOperator rounded_island_laplace(rounded_island_iodata, rounded_island_meshes);
  SurfaceResponseOperator rounded_island_response(rounded_island_iodata,
                                                  rounded_island_laplace);

  // The same fillet represented by coarse quadratic edges must select the same four
  // rounded-corner coupons and leave the same straight response intervals.
  std::vector<std::unique_ptr<Mesh>> high_order_rounded_island_meshes;
  high_order_rounded_island_meshes.push_back(
      std::make_unique<Mesh>(MakeIslandMesh(true, false, false, false, false, true)));
  LaplaceOperator high_order_rounded_island_laplace(rounded_island_iodata,
                                                    high_order_rounded_island_meshes);
  SurfaceResponseOperator high_order_rounded_island_response(
      rounded_island_iodata, high_order_rounded_island_laplace);
  CHECK(high_order_rounded_island_response.GetPatchCount() ==
        rounded_island_response.GetPatchCount());
  CHECK(high_order_rounded_island_response.GetBasisSize() ==
        rounded_island_response.GetBasisSize());
  CHECK_THAT(high_order_rounded_island_response.GetPatchWeight(),
             WithinRel(rounded_island_response.GetPatchWeight(), 1.0e-12));

#endif
}

TEST_CASE_METHOD(test::SurfaceResponseFiles,
                 "SurfaceResponseOperator high-order spatial cluster round trips",
                 "[surfaceresponseoperator][3d][highorder][placement][signature][Long]["
                 "Serial][Parallel]")
{
#if !defined(MFEM_USE_GSLIB)
  SKIP("SurfaceResponseOperator requires MFEM_USE_GSLIB");
#else
  json island_config = IslandConfig();
  // Nearby curved conductors require a spatial cluster rather than independent corner
  // responses. Preflight must retain sampled face loops without an eight-vertex cap, and
  // its exact mask must round-trip through the process library on the same geometry.
  json high_order_spatial_config = HighOrderSpatialConfig();
  IoData high_order_spatial_iodata(high_order_spatial_config, false);
  high_order_spatial_iodata.boundaries.cracked_attributes.insert(9);
  high_order_spatial_iodata.boundaries.cracked_attributes.insert(10);
  auto high_order_spatial_mesh = MakeIslandMesh(true, false, false, true, false, true);
  const auto high_order_spatial_requirements_path =
      temp.temp_dir / "surface-response-requirements-high-order-spatial.json";
  WriteSurfaceResponseRequirements(high_order_spatial_iodata, *high_order_spatial_mesh,
                                   high_order_spatial_requirements_path.string());
  std::ifstream high_order_spatial_requirements_input(high_order_spatial_requirements_path);
  REQUIRE(high_order_spatial_requirements_input);
  const auto high_order_spatial_requirements =
      json::parse(high_order_spatial_requirements_input);
  // The plan-view mask round trip is a legacy-classifier contract: the version-2
  // Requirements are derived from the identification features (matched by signature), so
  // the mask lives in the LegacyRequirements comparison table.
  const auto high_order_spatial_requirement =
      std::find_if(high_order_spatial_requirements["LegacyRequirements"].begin(),
                   high_order_spatial_requirements["LegacyRequirements"].end(),
                   [](const auto &requirement)
                   {
                     return requirement["Topology"] == "SpatialEdgeCluster" &&
                            requirement["Geometry"].contains("PlanViewFacets") &&
                            requirement["Geometry"].contains("PlanViewBoundary");
                   });
  REQUIRE(high_order_spatial_requirement !=
          high_order_spatial_requirements["LegacyRequirements"].end());
  const auto &high_order_spatial_facets =
      (*high_order_spatial_requirement)["Geometry"]["PlanViewFacets"];
  CHECK(std::any_of(high_order_spatial_facets.begin(), high_order_spatial_facets.end(),
                    [](const auto &facet) { return facet["Points"].size() > 8; }));

  std::ifstream exact_high_order_spatial_library_input(spatial_cluster_library_3d_path);
  REQUIRE(exact_high_order_spatial_library_input);
  auto exact_high_order_spatial_library =
      json::parse(exact_high_order_spatial_library_input);
  auto &exact_high_order_spatial_model = exact_high_order_spatial_library["Models"].back();
  exact_high_order_spatial_model["Name"] = "high-order-curved-spatial-exact-mask";
  exact_high_order_spatial_model["Edges"] =
      (*high_order_spatial_requirement)["Geometry"]["Edges"];
  for (auto &edge : exact_high_order_spatial_model["Edges"])
  {
    edge["BoundaryCondition"] = "PEC";
  }
  exact_high_order_spatial_model["PlanViewBoundary"] =
      (*high_order_spatial_requirement)["Geometry"]["PlanViewBoundary"];
  const auto exact_high_order_spatial_library_path =
      temp.temp_dir / "fabrication-process-high-order-spatial-exact-mask-3d.json";
  std::ofstream exact_high_order_spatial_library_output(
      exact_high_order_spatial_library_path);
  exact_high_order_spatial_library_output << exact_high_order_spatial_library.dump(2)
                                          << "\n";
  exact_high_order_spatial_library_output.close();

  high_order_spatial_config["Solver"]["Electrostatic"]["ResponseCorrection"]["Library"] =
      exact_high_order_spatial_library_path.string();
  IoData exact_high_order_spatial_iodata(high_order_spatial_config, false);
  exact_high_order_spatial_iodata.boundaries.cracked_attributes.insert(9);
  exact_high_order_spatial_iodata.boundaries.cracked_attributes.insert(10);
  const auto exact_high_order_spatial_requirements_path =
      temp.temp_dir / "surface-response-requirements-high-order-spatial-exact.json";
  WriteSurfaceResponseRequirements(exact_high_order_spatial_iodata,
                                   *high_order_spatial_mesh,
                                   exact_high_order_spatial_requirements_path.string());
  std::ifstream exact_high_order_spatial_requirements_input(
      exact_high_order_spatial_requirements_path);
  REQUIRE(exact_high_order_spatial_requirements_input);
  const auto exact_high_order_spatial_requirements =
      json::parse(exact_high_order_spatial_requirements_input);
  // Legacy Edges/PlanViewBoundary model round trip (comparison table only; a version-2
  // model matches through its Signature, tested below).
  const auto matched_high_order_spatial =
      std::find_if(exact_high_order_spatial_requirements["LegacyRequirements"].begin(),
                   exact_high_order_spatial_requirements["LegacyRequirements"].end(),
                   [](const auto &requirement)
                   {
                     return requirement["Topology"] == "SpatialEdgeCluster" &&
                            requirement["Status"] == "Exact" &&
                            requirement["SelectedModels"][0]["Name"] ==
                                "high-order-curved-spatial-exact-mask";
                   });
  CHECK(matched_high_order_spatial !=
        exact_high_order_spatial_requirements["LegacyRequirements"].end());

  // Version-2 round trip: a model carrying a feature's canonical Signature matches that
  // feature by key (and only that feature), independent of the mesh ordering.
  {
    const auto &identification = high_order_spatial_requirements["Identification"];
    REQUIRE(identification["Version"] == 2);
    const auto cluster_feature = std::find_if(
        identification["Features"].begin(), identification["Features"].end(),
        [](const auto &feature) { return feature["Type"] == "SpatialEdgeCluster"; });
    REQUIRE(cluster_feature != identification["Features"].end());
    CHECK((*cluster_feature)["Match"]["Status"] == "Missing");
    auto signature_library = exact_high_order_spatial_library;
    auto &signature_model = signature_library["Models"].back();
    signature_model["Name"] = "v2-signature-cluster";
    signature_model["Signature"] = (*cluster_feature)["Signature"];
    const auto signature_library_path =
        temp.temp_dir / "fabrication-process-v2-signature-3d.json";
    high_order_spatial_config["Solver"]["Electrostatic"]["ResponseCorrection"]["Library"] =
        signature_library_path.string();
    signature_model.erase("Edges");  // keyed and built by the Signature alone
    std::ofstream signature_library_output(signature_library_path);
    signature_library_output << signature_library.dump(2) << "\n";
    signature_library_output.close();
    IoData signature_iodata(high_order_spatial_config, false);
    signature_iodata.boundaries.cracked_attributes.insert(9);
    signature_iodata.boundaries.cracked_attributes.insert(10);
    const auto signature_requirements_path =
        temp.temp_dir / "surface-response-requirements-v2-signature.json";
    WriteSurfaceResponseRequirements(signature_iodata, *high_order_spatial_mesh,
                                     signature_requirements_path.string());
    std::ifstream signature_requirements_input(signature_requirements_path);
    REQUIRE(signature_requirements_input);
    const auto signature_requirements = json::parse(signature_requirements_input);
    const auto &signature_identification = signature_requirements["Identification"];
    CHECK(signature_identification["GeometryDigest"] == identification["GeometryDigest"]);
    int matched_clusters = 0;
    for (const auto &feature : signature_identification["Features"])
    {
      if (feature["Match"]["Status"] == "Matched")
      {
        CHECK(feature["Type"] == "SpatialEdgeCluster");
        CHECK(feature["Hash"] == (*cluster_feature)["Hash"]);
        CHECK(feature["Match"]["Model"] == "v2-signature-cluster");
        matched_clusters++;
      }
    }
    CHECK(matched_clusters >= 1);
    CHECK(std::any_of(signature_requirements["Requirements"].begin(),
                      signature_requirements["Requirements"].end(),
                      [](const auto &requirement)
                      {
                        return requirement["Topology"] == "SpatialEdgeCluster" &&
                               requirement["Status"] == "Exact";
                      }));

    // A version-2 model that also stores its Edges in the canonical frame (the library
    // builder's contract: Point = P x R, Interval along gap x normal, the Signature's own
    // portions, an arc portion chorded at 5 deg / 0.25 R as signature_library
    // .cluster_plan_view_edges does) is placed with the identity map, so the patch dry run
    // maps every model edge endpoint onto the feature's portions in the mesh: a straight
    // edge endpoint onto a portion endpoint, a chord endpoint onto the fitted arc within
    // the claimed angular range. (A re-canonicalisation of the stored Edges read the arc
    // chords / the tangent of the opposite sign in another frame and placed the coupon
    // elsewhere on the device: the 2394fdb0c failure class.)
    {
      auto edges_library = signature_library;
      auto &edges_model = edges_library["Models"].back();
      edges_model["Name"] = "v2-signature-cluster-edges";
      const double radius = edges_library["MatchingRadius"].get<double>();
      edges_model["Interfaces"] = SignatureInterfaces((*cluster_feature)["Signature"]);
      edges_model["Edges"] = ChordedSignatureEdges((*cluster_feature)["Signature"], radius);
      CHECK(edges_model["Edges"].size() >
            (*cluster_feature)["Signature"]["Portions"].size());
      const auto edges_library_path =
          temp.temp_dir / "fabrication-process-v2-signature-edges-3d.json";
      std::ofstream edges_library_output(edges_library_path);
      edges_library_output << edges_library.dump(2) << "\n";
      edges_library_output.close();
      auto edges_config = high_order_spatial_config;
      auto &edges_correction =
          edges_config["Solver"]["Electrostatic"]["ResponseCorrection"];
      edges_correction["Library"] = edges_library_path.string();
      edges_correction.erase("PatchConstruction");  // Features (the default)
      IoData edges_iodata(edges_config, false);
      edges_iodata.boundaries.cracked_attributes.insert(9);
      edges_iodata.boundaries.cracked_attributes.insert(10);
      const auto edges_requirements_path =
          temp.temp_dir / "surface-response-requirements-v2-signature-edges.json";
      const auto edges_patches_path = temp.temp_dir / "surface-response-patches.csv";
      WriteSurfaceResponseRequirements(edges_iodata, *high_order_spatial_mesh,
                                       edges_requirements_path.string());
      std::ifstream edges_requirements_input(edges_requirements_path);
      REQUIRE(edges_requirements_input);
      const auto edges_requirements = json::parse(edges_requirements_input);
      const auto &edges_identification = edges_requirements["Identification"];
      const auto patch_rows = ReadPatchRows(edges_patches_path);
      int checked_clusters = 0;
      for (const auto &feature : edges_identification["Features"])
      {
        if (feature["Type"] != "SpatialEdgeCluster")
        {
          continue;
        }
        REQUIRE(feature["Match"]["Status"] == "Matched");
        REQUIRE(feature["Match"]["Model"] == "v2-signature-cluster-edges");
        CheckPlacedModelEdges(feature, edges_identification, patch_rows,
                              edges_model["Edges"], 1.0e-6 * radius);
        checked_clusters++;
      }
      CHECK(checked_clusters >= 1);
    }
  }

#endif
}

TEST_CASE_METHOD(
    test::SurfaceResponseFiles,
    "SurfaceResponseOperator spatial cluster contracts on a linear rounded mesh",
    "[surfaceresponseoperator][3d][placement][signature][Serial][Parallel]")
{
#if !defined(MFEM_USE_GSLIB)
  SKIP("SurfaceResponseOperator requires MFEM_USE_GSLIB");
#else
  json island_config = IslandConfig();
  // The version-2 library contracts of the high-order spatial cluster case on the same two
  // rounded islands with linear fillets (the arcs are the fitted linear fillets): the
  // identification costs 1e-2 of the quadratic mesh's.
  json linear_spatial_config = HighOrderSpatialConfig();
  IoData linear_spatial_iodata(linear_spatial_config, false);
  linear_spatial_iodata.boundaries.cracked_attributes.insert(9);
  linear_spatial_iodata.boundaries.cracked_attributes.insert(10);
  auto linear_spatial_mesh = MakeIslandMesh(true, false, false, true);
  const auto linear_spatial_requirements_path =
      temp.temp_dir / "surface-response-requirements-linear-spatial.json";
  WriteSurfaceResponseRequirements(linear_spatial_iodata, *linear_spatial_mesh,
                                   linear_spatial_requirements_path.string());
  std::ifstream linear_spatial_requirements_input(linear_spatial_requirements_path);
  REQUIRE(linear_spatial_requirements_input);
  const auto linear_spatial_requirements = json::parse(linear_spatial_requirements_input);
  // The plan-view mask round trip is a legacy-classifier contract: the version-2
  // Requirements are derived from the identification features (matched by signature), so
  // the mask lives in the LegacyRequirements comparison table.
  const auto linear_spatial_requirement =
      std::find_if(linear_spatial_requirements["LegacyRequirements"].begin(),
                   linear_spatial_requirements["LegacyRequirements"].end(),
                   [](const auto &requirement)
                   {
                     return requirement["Topology"] == "SpatialEdgeCluster" &&
                            requirement["Geometry"].contains("PlanViewFacets") &&
                            requirement["Geometry"].contains("PlanViewBoundary");
                   });
  REQUIRE(linear_spatial_requirement !=
          linear_spatial_requirements["LegacyRequirements"].end());
  std::ifstream exact_linear_spatial_library_input(spatial_cluster_library_3d_path);
  REQUIRE(exact_linear_spatial_library_input);
  auto exact_linear_spatial_library = json::parse(exact_linear_spatial_library_input);
  auto &exact_linear_spatial_model = exact_linear_spatial_library["Models"].back();
  exact_linear_spatial_model["Name"] = "linear-curved-spatial-exact-mask";
  exact_linear_spatial_model["Edges"] = (*linear_spatial_requirement)["Geometry"]["Edges"];
  for (auto &edge : exact_linear_spatial_model["Edges"])
  {
    edge["BoundaryCondition"] = "PEC";
  }
  exact_linear_spatial_model["PlanViewBoundary"] =
      (*linear_spatial_requirement)["Geometry"]["PlanViewBoundary"];
  const auto exact_linear_spatial_library_path =
      temp.temp_dir / "fabrication-process-linear-spatial-exact-mask-3d.json";
  std::ofstream exact_linear_spatial_library_output(exact_linear_spatial_library_path);
  exact_linear_spatial_library_output << exact_linear_spatial_library.dump(2) << "\n";
  exact_linear_spatial_library_output.close();

  // Version-2 round trip: a model carrying a feature's canonical Signature matches that
  // feature by key (and only that feature), independent of the mesh ordering.
  {
    const auto &identification = linear_spatial_requirements["Identification"];
    REQUIRE(identification["Version"] == 2);
    const auto cluster_feature = std::find_if(
        identification["Features"].begin(), identification["Features"].end(),
        [](const auto &feature) { return feature["Type"] == "SpatialEdgeCluster"; });
    REQUIRE(cluster_feature != identification["Features"].end());
    CHECK((*cluster_feature)["Match"]["Status"] == "Missing");
    auto signature_library = exact_linear_spatial_library;
    auto &signature_model = signature_library["Models"].back();
    signature_model["Name"] = "v2-signature-cluster";
    signature_model["Signature"] = (*cluster_feature)["Signature"];
    const auto signature_library_path =
        temp.temp_dir / "fabrication-process-v2-signature-3d.json";
    linear_spatial_config["Solver"]["Electrostatic"]["ResponseCorrection"]["Library"] =
        signature_library_path.string();
    // A Signature model's Edges must be the Signature's portions in its canonical frame:
    // the legacy Edges inherited here (the classifier's mesh-frame description, process
    // normal +y) are refused at library load (fail closed), never re-canonicalised.
    {
      std::ofstream output(signature_library_path);
      output << signature_library.dump(2) << "\n";
    }
    {
      IoData legacy_edges_iodata(linear_spatial_config, false);
      legacy_edges_iodata.boundaries.cracked_attributes.insert(9);
      legacy_edges_iodata.boundaries.cracked_attributes.insert(10);
      CHECK_THROWS_WITH(
          WriteSurfaceResponseRequirements(
              legacy_edges_iodata, *linear_spatial_mesh,
              (temp.temp_dir / "surface-response-requirements-v2-legacy-edges.json")
                  .string()),
          Catch::Matchers::ContainsSubstring("canonical frame"));
    }
    // The load-time contract is not geometric only: a cluster whose geometry is symmetric
    // under x -> -x but whose labels are not (two facing edges of different conductors, the
    // left one with the full interface set, the right one with SA alone) has mirrored
    // Edges that land on the portions exactly, with every Conductor / InterfaceSlot
    // attached to the wrong portion; the correct Edges load, the mirrored Edges and the
    // Edges with the interface slots exchanged are refused.
    {
      auto labelled_library = signature_library;
      auto &labelled_model = labelled_library["Models"].back();
      labelled_model["Name"] = "v2-signature-cluster-asymmetric-labels";
      labelled_model["Signature"] = {{"Type", "SpatialEdgeCluster"},
                                     {"EdgeCount", 2},
                                     {"Portions",
                                      {{{"Conductor", 1},
                                        {"Gap", {1.0, 0.0}},
                                        {"Interfaces", {"MA", "MS", "SA"}},
                                        {"Law", "{\"Type\":\"PEC\"}"},
                                        {"P", {-1.0, -1.0, -1.0, 1.0}}},
                                       {{"Conductor", 2},
                                        {"Gap", {-1.0, 0.0}},
                                        {"Interfaces", {"SA"}},
                                        {"Law", "{\"Type\":\"PEC\"}"},
                                        {"P", {1.0, -1.0, 1.0, 1.0}}}}},
                                     {"Vertices", json::array()}};
      const double radius = labelled_library["MatchingRadius"].get<double>();
      labelled_model["Interfaces"] = SignatureInterfaces(labelled_model["Signature"]);
      labelled_model["Edges"] = ChordedSignatureEdges(labelled_model["Signature"], radius);
      REQUIRE(labelled_model["Edges"].size() == 2);
      REQUIRE(labelled_model["Edges"][0]["InterfaceSlot"] == 0);  // [MA, MS, SA]
      REQUIRE(labelled_model["Edges"][1]["InterfaceSlot"] == 1);  // [SA]
      const auto labelled_library_path =
          temp.temp_dir / "fabrication-process-v2-signature-asymmetric-labels-3d.json";
      auto labelled_config = linear_spatial_config;
      labelled_config["Solver"]["Electrostatic"]["ResponseCorrection"]["Library"] =
          labelled_library_path.string();
      auto WriteAndIdentify = [&](const json &library, const std::string &suffix)
      {
        {
          std::ofstream output(labelled_library_path);
          output << library.dump(2) << "\n";
        }
        IoData iodata(labelled_config, false);
        iodata.boundaries.cracked_attributes.insert(9);
        iodata.boundaries.cracked_attributes.insert(10);
        WriteSurfaceResponseRequirements(
            iodata, *linear_spatial_mesh,
            (temp.temp_dir /
             ("surface-response-requirements-v2-labels-" + suffix + ".json"))
                .string());
      };
      CHECK_NOTHROW(WriteAndIdentify(labelled_library, "correct"));
      // Mirror x -> -x: Point.x and GapDirection.x change sign, the Interval (along
      // gap x normal, which flips) is reversed; every endpoint lands on the other portion.
      auto mirrored_library = labelled_library;
      for (auto &edge : mirrored_library["Models"].back()["Edges"])
      {
        edge["Point"][0] = -edge["Point"][0].get<double>();
        edge["GapDirection"][0] = -edge["GapDirection"][0].get<double>();
        const auto interval = edge["Interval"].get<std::array<double, 2>>();
        edge["Interval"] = {-interval[1], -interval[0]};
      }
      CHECK_THROWS_WITH(WriteAndIdentify(mirrored_library, "mirrored"),
                        Catch::Matchers::ContainsSubstring("relabelled Conductor") &&
                            Catch::Matchers::ContainsSubstring(
                                "(Conductor 1, InterfaceSlot 0 = [MA, MS, SA]) lies on "
                                "Signature portion 1 (Conductor 2, Interfaces [SA])"));
      // The right geometry and conductors with the interface slots exchanged.
      auto swapped_slots_library = labelled_library;
      swapped_slots_library["Models"].back()["Edges"][0]["InterfaceSlot"] = 1;
      swapped_slots_library["Models"].back()["Edges"][1]["InterfaceSlot"] = 0;
      CHECK_THROWS_WITH(WriteAndIdentify(swapped_slots_library, "swapped-slots"),
                        Catch::Matchers::ContainsSubstring("relabelled Conductor") &&
                            Catch::Matchers::ContainsSubstring(
                                "(Conductor 1, InterfaceSlot 1 = [SA]) lies on Signature "
                                "portion 0 (Conductor 1, Interfaces [MA, MS, SA])"));
    }
    signature_model.erase("Edges");  // keyed and built by the Signature alone
    std::ofstream signature_library_output(signature_library_path);
    signature_library_output << signature_library.dump(2) << "\n";
    signature_library_output.close();
    IoData signature_iodata(linear_spatial_config, false);
    signature_iodata.boundaries.cracked_attributes.insert(9);
    signature_iodata.boundaries.cracked_attributes.insert(10);
    const auto signature_requirements_path =
        temp.temp_dir / "surface-response-requirements-v2-signature.json";
    WriteSurfaceResponseRequirements(signature_iodata, *linear_spatial_mesh,
                                     signature_requirements_path.string());
    std::ifstream signature_requirements_input(signature_requirements_path);
    REQUIRE(signature_requirements_input);
    const auto signature_requirements = json::parse(signature_requirements_input);
    const auto &signature_identification = signature_requirements["Identification"];
    CHECK(signature_identification["GeometryDigest"] == identification["GeometryDigest"]);
    int matched_clusters = 0;
    for (const auto &feature : signature_identification["Features"])
    {
      if (feature["Match"]["Status"] == "Matched")
      {
        CHECK(feature["Type"] == "SpatialEdgeCluster");
        CHECK(feature["Hash"] == (*cluster_feature)["Hash"]);
        CHECK(feature["Match"]["Model"] == "v2-signature-cluster");
        matched_clusters++;
      }
    }
    CHECK(matched_clusters >= 1);
    CHECK(std::any_of(signature_requirements["Requirements"].begin(),
                      signature_requirements["Requirements"].end(),
                      [](const auto &requirement)
                      {
                        return requirement["Topology"] == "SpatialEdgeCluster" &&
                               requirement["Status"] == "Exact";
                      }));
    CHECK(signature_requirements["Summary"]["Counts"]["QuantumNearMatched"] == 0);

    // Quantum near-match (block (b) DESIGN section 4, decision 303): a model whose
    // Signature differs from the feature's by ONE signature quantum in one number (a 1e-6 R
    // shift of a portion end: the same geometry at the grid) matches the feature, Status
    // Matched / Exact, with Match.QuantumNearMatch naming both keys and the differing
    // number, counted in Summary.QuantumNearMatched; exactly 4 on-grid quanta still match
    // (the half-quantum inclusive rule, decision 317 MINOR-1); a 5-quantum shift stays
    // Missing; two library models within 8 quanta of each other (8 on-grid included) are
    // refused, 9 apart are admitted.
    {
      auto Perturbed = [&](int quanta, const std::string &name)
      {
        auto library = signature_library;
        auto &model = library["Models"].back();
        model["Name"] = name;
        model["Signature"]["Portions"][0]["P"][0] =
            model["Signature"]["Portions"][0]["P"][0].get<double>() + quanta * 1.0e-6;
        return library;
      };
      auto Preflight = [&](const json &library, const std::string &tag)
      {
        const auto library_path =
            temp.temp_dir / ("fabrication-process-near-match-" + tag + "-3d.json");
        std::ofstream output(library_path);
        output << library.dump(2) << "\n";
        output.close();
        auto config = linear_spatial_config;
        config["Solver"]["Electrostatic"]["ResponseCorrection"]["Library"] =
            library_path.string();
        IoData iodata(config, false);
        iodata.boundaries.cracked_attributes.insert(9);
        iodata.boundaries.cracked_attributes.insert(10);
        const auto requirements_path =
            temp.temp_dir / ("surface-response-requirements-near-match-" + tag + ".json");
        WriteSurfaceResponseRequirements(iodata, *linear_spatial_mesh,
                                         requirements_path.string());
        std::ifstream requirements_input(requirements_path);
        REQUIRE(requirements_input);
        return json::parse(requirements_input);
      };
      const json near = Preflight(Perturbed(1, "near-one-quantum"), "one");
      const auto [model_key, model_hash] = SignatureKeyAndHash(
          Perturbed(1, "near-one-quantum")["Models"].back()["Signature"],
          "SpatialEdgeCluster");
      (void)model_key;
      int near_matched = 0;
      for (const auto &feature : near["Identification"]["Features"])
      {
        if (feature["Type"] != "SpatialEdgeCluster" ||
            feature["Hash"] != (*cluster_feature)["Hash"])
        {
          continue;
        }
        REQUIRE(feature["Match"]["Status"] == "Matched");
        CHECK(feature["Match"]["Model"] == "near-one-quantum");
        CHECK_THAT(feature["Match"]["Deviation"].get<double>(),
                   WithinAbs(1.0 / 4.5, 1.0e-6));
        REQUIRE(feature["Match"].contains("QuantumNearMatch"));
        const auto &record = feature["Match"]["QuantumNearMatch"];
        CHECK(record["ModelKey"] == model_hash);
        CHECK(record["FeatureKey"] == feature["Hash"]);
        CHECK_THAT(record["MaxDeltaQuanta"].get<double>(), WithinAbs(1.0, 1.0e-6));
        CHECK(record["DifferingNumbers"]["Count"] == 1);
        CHECK(record["DifferingNumbers"]["Paths"] == json({"Portions[0].P[0]"}));
        near_matched++;
      }
      CHECK(near_matched >= 1);
      CHECK(near["Summary"]["Counts"]["QuantumNearMatched"] == near_matched);
      CHECK(near["Summary"]["QuantumNearMatch"]["Keys"].size() == 1);
      CHECK(near["Summary"]["QuantumNearMatch"]["Keys"][0]["ModelKey"] == model_hash);
      CHECK(near["Summary"]["QuantumNearMatch"]["Keys"][0]["FeatureKey"] ==
            (*cluster_feature)["Hash"]);
      CHECK(std::any_of(near["Requirements"].begin(), near["Requirements"].end(),
                        [](const auto &requirement)
                        {
                          return requirement["Topology"] == "SpatialEdgeCluster" &&
                                 requirement["Status"] == "Exact" &&
                                 requirement["SelectedModels"][0]["Name"] ==
                                     "near-one-quantum";
                        }));
      const json four = Preflight(Perturbed(4, "near-four-quanta"), "four");
      int four_matched = 0;
      for (const auto &feature : four["Identification"]["Features"])
      {
        if (feature["Type"] == "SpatialEdgeCluster" &&
            feature["Hash"] == (*cluster_feature)["Hash"])
        {
          REQUIRE(feature["Match"]["Status"] == "Matched");
          CHECK(feature["Match"]["Model"] == "near-four-quanta");
          CHECK(feature["Match"]["Deviation"].get<double>() <= 1.0);
          CHECK_THAT(feature["Match"]["QuantumNearMatch"]["MaxDeltaQuanta"].get<double>(),
                     WithinAbs(4.0, 1.0e-6));
          four_matched++;
        }
      }
      CHECK(four_matched >= 1);
      CHECK(four["Summary"]["Counts"]["QuantumNearMatched"] == four_matched);
      const json far = Preflight(Perturbed(5, "far-five-quanta"), "five");
      for (const auto &feature : far["Identification"]["Features"])
      {
        if (feature["Type"] == "SpatialEdgeCluster" &&
            feature["Hash"] == (*cluster_feature)["Hash"])
        {
          CHECK(feature["Match"]["Status"] == "Missing");
          CHECK(!feature["Match"].contains("QuantumNearMatch"));
        }
      }
      CHECK(far["Summary"]["Counts"]["QuantumNearMatched"] == 0);
      // Two models one geometry: the exact model and the one-quantum model together.
      auto duplicate = signature_library;
      duplicate["Models"].push_back(Perturbed(1, "near-one-quantum")["Models"].back());
      CHECK_THROWS_WITH(Preflight(duplicate, "duplicate"),
                        Catch::Matchers::ContainsSubstring("two models, one geometry"));
      auto duplicate_eight = signature_library;
      duplicate_eight["Models"].push_back(
          Perturbed(8, "near-eight-quanta")["Models"].back());
      CHECK_THROWS_WITH(Preflight(duplicate_eight, "duplicate-eight"),
                        Catch::Matchers::ContainsSubstring("two models, one geometry"));
      auto distinct_nine = signature_library;
      distinct_nine["Models"].push_back(Perturbed(9, "far-nine-quanta")["Models"].back());
      const json nine = Preflight(distinct_nine, "distinct-nine");
      CHECK(nine["Summary"]["Counts"]["QuantumNearMatched"] == 0);  // the exact model wins
    }

    // A version-2 model that also stores its Edges in the canonical frame (the library
    // builder's contract: Point = P x R, Interval along gap x normal, the Signature's own
    // portions, an arc portion chorded at 5 deg / 0.25 R as signature_library
    // .cluster_plan_view_edges does) is placed with the identity map, so the patch dry run
    // maps every model edge endpoint onto the feature's portions in the mesh: a straight
    // edge endpoint onto a portion endpoint, a chord endpoint onto the fitted arc within
    // the claimed angular range. (A re-canonicalisation of the stored Edges read the arc
    // chords / the tangent of the opposite sign in another frame and placed the coupon
    // elsewhere on the device: the 2394fdb0c failure class.)
    {
      auto edges_library = signature_library;
      auto &edges_model = edges_library["Models"].back();
      edges_model["Name"] = "v2-signature-cluster-edges";
      const double radius = edges_library["MatchingRadius"].get<double>();
      edges_model["Interfaces"] = SignatureInterfaces((*cluster_feature)["Signature"]);
      edges_model["Edges"] = ChordedSignatureEdges((*cluster_feature)["Signature"], radius);
      CHECK(edges_model["Edges"].size() >
            (*cluster_feature)["Signature"]["Portions"].size());
      const auto edges_library_path =
          temp.temp_dir / "fabrication-process-v2-signature-edges-3d.json";
      std::ofstream edges_library_output(edges_library_path);
      edges_library_output << edges_library.dump(2) << "\n";
      edges_library_output.close();
      auto edges_config = linear_spatial_config;
      auto &edges_correction =
          edges_config["Solver"]["Electrostatic"]["ResponseCorrection"];
      edges_correction["Library"] = edges_library_path.string();
      edges_correction.erase("PatchConstruction");  // Features (the default)
      IoData edges_iodata(edges_config, false);
      edges_iodata.boundaries.cracked_attributes.insert(9);
      edges_iodata.boundaries.cracked_attributes.insert(10);
      const auto edges_requirements_path =
          temp.temp_dir / "surface-response-requirements-v2-signature-edges.json";
      const auto edges_patches_path = temp.temp_dir / "surface-response-patches.csv";
      WriteSurfaceResponseRequirements(edges_iodata, *linear_spatial_mesh,
                                       edges_requirements_path.string());
      std::ifstream edges_requirements_input(edges_requirements_path);
      REQUIRE(edges_requirements_input);
      const auto edges_requirements = json::parse(edges_requirements_input);
      const auto &edges_identification = edges_requirements["Identification"];
      const auto patch_rows = ReadPatchRows(edges_patches_path);
      int checked_clusters = 0;
      for (const auto &feature : edges_identification["Features"])
      {
        if (feature["Type"] != "SpatialEdgeCluster")
        {
          continue;
        }
        REQUIRE(feature["Match"]["Status"] == "Matched");
        REQUIRE(feature["Match"]["Model"] == "v2-signature-cluster-edges");
        CheckPlacedModelEdges(feature, edges_identification, patch_rows,
                              edges_model["Edges"], 1.0e-6 * radius);
        checked_clusters++;
      }
      CHECK(checked_clusters >= 1);
    }
  }

#endif
}

TEST_CASE_METHOD(test::SurfaceResponseFiles,
                 "SurfaceResponseOperator filleted finger arc portions",
                 "[surfaceresponseoperator][3d][curvature][placement][Serial][Parallel]")
{
#if !defined(MFEM_USE_GSLIB)
  SKIP("SurfaceResponseOperator requires MFEM_USE_GSLIB");
#else
  json island_config = IslandConfig();
  // Version-2 round trip with ARC portions (option A): a 1.5 R x 6 R finger with 0.5 R
  // fillets is one SpatialEdgeCluster per end (two fillet arcs, the end edge, R along both
  // sides). The library builder chords each arc portion of the Signature at 5 deg / 0.25 R
  // and stores the chords as the model's Edges in the canonical frame
  // (signature_library.cluster_plan_view_edges); the patch construction must place such a
  // model with the identity (its chords cannot be re-canonicalised into arc portions), so
  // every chord endpoint lands on the feature's arc in the mesh and every straight edge
  // endpoint on a portion endpoint, for both finger ends (rotated frames) and independent
  // of the mesh's chord count (4 / 8 chords per fillet).
  {
    constexpr double finger_R = 0.2;  // the spatial library's MatchingRadius
    constexpr double fillet_radius = 0.5 * finger_R;
    constexpr double finger_x0 = 0.4, finger_x1 = 1.6, finger_z0 = 0.35, finger_z1 = 0.65;
    auto MakeFilletedFingerMesh = [&](int elements_per_fillet)
    {
      const double h = fillet_radius / elements_per_fillet;
      const int nx = static_cast<int>(std::lround(2.0 / h)),
                nz = static_cast<int>(std::lround(1.0 / h));
      mfem::Mesh serial =
          mfem::Mesh::MakeCartesian3D(nx, 4, nz, mfem::Element::HEXAHEDRON, 2.0, 1.0, 1.0);
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
        bool on_plane = true, inside = true;
        for (const int vertex : vertices)
        {
          const double *point = serial.GetVertex(vertex);
          on_plane = on_plane && std::abs(point[1] - 0.5) < 1.0e-12;
          inside = inside && point[0] >= finger_x0 - 1.0e-12 &&
                   point[0] <= finger_x1 + 1.0e-12 && point[2] >= finger_z0 - 1.0e-12 &&
                   point[2] <= finger_z1 + 1.0e-12;
        }
        if (on_plane && inside)
        {
          serial.AddBdrElement(serial.GetFace(face)->Duplicate(&serial));
          serial.SetBdrAttribute(serial.GetNBE() - 1, 9);
        }
      }
      // Each corner square [0, r]^2 of the finger (local coordinates from the fillet
      // centre outward) is mapped square-by-square onto concentric quarter circles: the
      // level-m square boundary onto the arc of radius m, so the outline becomes a
      // polyline of 2 elements_per_fillet chords on the design circle and no element is
      // inverted.
      for (int vertex = 0; vertex < serial.GetNV(); vertex++)
      {
        double *point = serial.GetVertex(vertex);
        if (std::abs(point[1] - 0.5) > 1.0e-12)
        {
          continue;
        }
        for (const double sign_x : {-1.0, 1.0})
        {
          for (const double sign_z : {-1.0, 1.0})
          {
            const double corner_x = sign_x > 0 ? finger_x1 : finger_x0;
            const double corner_z = sign_z > 0 ? finger_z1 : finger_z0;
            const double center_x = corner_x - sign_x * fillet_radius;
            const double center_z = corner_z - sign_z * fillet_radius;
            const double local_x = sign_x * (point[0] - center_x);
            const double local_z = sign_z * (point[2] - center_z);
            if (local_x < -1.0e-12 || local_x > fillet_radius + 1.0e-12 ||
                local_z < -1.0e-12 || local_z > fillet_radius + 1.0e-12)
            {
              continue;
            }
            const double m = std::max(local_x, local_z);
            if (m <= 1.0e-12)
            {
              break;
            }
            const double angle =
                std::abs(local_z - m) <= 1.0e-12
                    ? 0.5 * std::acos(-1.0) - 0.25 * std::acos(-1.0) * local_x / m
                    : 0.25 * std::acos(-1.0) * local_z / m;
            point[0] = center_x + sign_x * m * std::cos(angle);
            point[2] = center_z + sign_z * m * std::sin(angle);
            break;
          }
        }
      }
      serial.FinalizeTopology();
      serial.Finalize();
      return std::make_unique<mfem::ParMesh>(Mpi::World(), serial);
    };
    auto finger_config = island_config;
    finger_config["Solver"]["Electrostatic"]["ResponseCorrection"]["Library"] =
        spatial_cluster_library_3d_path.string();
    finger_config["Solver"]["Electrostatic"]["ResponseCorrection"]["UnmatchedPolicy"] =
        "Warn";
    finger_config["Solver"]["Electrostatic"]["ResponseCorrection"].erase(
        "PatchConstruction");
    // The feature's Signature from the coarse mesh (4 chords per fillet).
    std::optional<json> finger_signature;
    std::optional<std::string> finger_hash;
    const auto finger_library_path =
        temp.temp_dir / "fabrication-process-v2-arc-cluster.json";
    for (const int elements_per_fillet : {2, 4})
    {
      auto finger_mesh = MakeFilletedFingerMesh(elements_per_fillet);
      const auto finger_requirements_path =
          temp.temp_dir / "surface-response-requirements-arc-finger.json";
      const auto finger_patches_path = temp.temp_dir / "surface-response-patches.csv";
      if (!finger_signature)
      {
        IoData finger_iodata(finger_config, false);
        finger_iodata.boundaries.cracked_attributes.insert(9);
        WriteSurfaceResponseRequirements(finger_iodata, *finger_mesh,
                                         finger_requirements_path.string());
        std::ifstream input(finger_requirements_path);
        REQUIRE(input);
        const auto requirements = json::parse(input);
        for (const auto &feature : requirements["Identification"]["Features"])
        {
          if (feature["Type"] != "SpatialEdgeCluster")
          {
            continue;
          }
          int arcs = 0;
          for (const auto &portion : feature["Signature"]["Portions"])
          {
            arcs += portion.contains("Arc") ? 1 : 0;
          }
          REQUIRE(arcs == 2);
          REQUIRE(feature["Signature"]["EdgeCount"].get<int>() == 5);
          if (!finger_signature)
          {
            finger_signature = feature["Signature"];
            finger_hash = feature["Hash"].get<std::string>();
          }
          else
          {
            CHECK(feature["Hash"].get<std::string>() == *finger_hash);  // congruent ends
          }
        }
        REQUIRE(finger_signature);
        // The builder's model: verbatim Signature, Edges = straight portions and the
        // chorded arcs (5 deg / 0.25 R) in the canonical frame, gap x normal intervals.
        // The spatial cluster library's last model (BasisPoints, tolerances and the
        // coupled matrices) carries the finger's Signature and Edges.
        std::ifstream finger_library_input(spatial_cluster_library_3d_path);
        auto finger_library = json::parse(finger_library_input);
        auto &finger_model = finger_library["Models"].back();
        finger_model["Name"] = "v2-signature-arc-cluster-edges";
        finger_model["Signature"] = *finger_signature;
        for (const char *key : {"ConductorReferences", "OpenContourPaths", "SupportPoints",
                                "PlanViewBoundary"})
        {
          finger_model.erase(key);
        }
        finger_model["Interfaces"] = SignatureInterfaces(*finger_signature);
        finger_model["Edges"] = ChordedSignatureEdges(*finger_signature, finger_R);
        CHECK(finger_model["Edges"].size() >= 3 + 2 * 18);  // 18 chords per quarter circle
        std::ofstream output(finger_library_path);
        output << finger_library.dump(2) << "\n";
        output.close();
        finger_config["Solver"]["Electrostatic"]["ResponseCorrection"]["Library"] =
            finger_library_path.string();
      }
      IoData finger_iodata(finger_config, false);
      finger_iodata.boundaries.cracked_attributes.insert(9);
      WriteSurfaceResponseRequirements(finger_iodata, *finger_mesh,
                                       finger_requirements_path.string());
      std::ifstream input(finger_requirements_path);
      REQUIRE(input);
      const auto requirements = json::parse(input);
      const auto &identification = requirements["Identification"];
      const auto patch_rows = ReadPatchRows(finger_patches_path);
      std::ifstream library_input(finger_library_path);
      const json model_edges = json::parse(library_input)["Models"].back()["Edges"];
      int checked_clusters = 0;
      for (const auto &feature : identification["Features"])
      {
        if (feature["Type"] != "SpatialEdgeCluster")
        {
          continue;
        }
        INFO("elements per fillet " << elements_per_fillet << ", feature "
                                    << feature["Id"]);
        CHECK(feature["Hash"].get<std::string>() == *finger_hash);  // chord-count free
        REQUIRE(feature["Match"]["Status"] == "Matched");
        REQUIRE(feature["Match"]["Model"] == "v2-signature-arc-cluster-edges");
        CHECK(CheckPlacedModelEdges(feature, identification, patch_rows, model_edges,
                                    1.0e-6 * finger_R) == 2);  // both fillets claimed
        checked_clusters++;
      }
      CHECK(checked_clusters == 2);  // both finger ends
    }
  }

#endif
}

TEST_CASE_METHOD(test::SurfaceResponseFiles,
                 "SurfaceResponseOperator C4 curvature family rings",
                 "[surfaceresponseoperator][3d][curvature][features][Serial][Parallel]")
{
#if !defined(MFEM_USE_GSLIB)
  SKIP("SurfaceResponseOperator requires MFEM_USE_GSLIB");
#else
  json island_config = IslandConfig();
  // Curvature family on the three-dimensional FEATURES path (C4, decision 108): two thin
  // annuli of width 3 R on one plane. Annulus A (3 R .. 6 R) has two CurvedEdge features
  // (inner edge concave at kappa 1/3, outer convex at 1/6), matched to the family's cubic
  // blend at their kappa (a runtime model named <anchor>@<convexity>-kappa<k>-cubic with
  // Lagrange weights on the four nodes); annulus B (20 R .. 23 R) has two straight-like
  // isolated edges whose bend annotation and per-portion turns evaluate the first-order
  // term: every quadrature point splits between the anchor (model weight 1 - a) and the
  // family node at kappa 0.1 of the turn's convexity (a = kappa / 0.1 = 0.5 inner, 20 / 23
  // x 0.5 outer). A library without the concave nodes leaves the concave CurvedEdge
  // unmatched with the reason (never straight) and the concave first-order term at the
  // anchor alone (recorded).
  {
    constexpr double ring_R = 0.2;
    auto MakeRingMesh =
        [&](double extent, double h, const std::vector<std::array<double, 2>> &rings)
    {
      const int n = static_cast<int>(std::lround(extent / h));
      mfem::Mesh serial = mfem::Mesh::MakeCartesian3D(n, 4, n, mfem::Element::HEXAHEDRON,
                                                      extent, 1.0, extent);
      const double c = 0.5 * extent;
      auto Level = [&](const double *point)
      { return std::max(std::abs(point[0] - c), std::abs(point[2] - c)); };
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
        double level_min = mfem::infinity(), level_max = 0.0;
        for (const int vertex : vertices)
        {
          const double *point = serial.GetVertex(vertex);
          on_plane = on_plane && std::abs(point[1] - 0.5) < 1.0e-12;
          level_min = std::min(level_min, Level(point));
          level_max = std::max(level_max, Level(point));
        }
        if (!on_plane)
        {
          continue;
        }
        for (const auto &[r_in, r_out] : rings)
        {
          if (level_min >= r_in - 1.0e-9 && level_max <= r_out + 1.0e-9)
          {
            serial.AddBdrElement(serial.GetFace(face)->Duplicate(&serial));
            serial.SetBdrAttribute(serial.GetNBE() - 1, 9);
          }
        }
      }
      // Concentric squares of the (x, z) grid onto concentric circles (radial projection
      // of every vertex onto the circle of its square's half-side): the ring outlines
      // become inscribed polylines on the design circles at the grid's chord count.
      for (int vertex = 0; vertex < serial.GetNV(); vertex++)
      {
        double *point = serial.GetVertex(vertex);
        const double lx = point[0] - c, lz = point[2] - c;
        const double m = std::max(std::abs(lx), std::abs(lz));
        const double norm = std::hypot(lx, lz);
        if (norm > 1.0e-12)
        {
          point[0] = c + m * lx / norm;
          point[2] = c + m * lz / norm;
        }
      }
      serial.FinalizeTopology();
      serial.Finalize();
      return std::make_unique<mfem::ParMesh>(Mpi::World(), serial);
    };
    // Extent 12: the outer ring (radius 4.6) stays more than 2 R from the PEC box walls (a
    // planar edge within 2 R of a wall is a CrossLayer exclusion).
    auto ring_mesh = MakeRingMesh(
        12.0, 0.1, {{3.0 * ring_R, 6.0 * ring_R}, {20.0 * ring_R, 23.0 * ring_R}});
    // Library: the SA-only isolated anchor at R 0.2 (matched by key) and the CurvedEdge
    // nodes of both convexities at kappa 0.1 / 0.25 / 0.5 / 0.8 (the anchor's matrices; the
    // weights are checked, not the values).
    const auto ring_library_path = temp.temp_dir / "fabrication-process-curved-3d.json";
    const auto ring_convex_only_path =
        temp.temp_dir / "fabrication-process-curved-convex-only-3d.json";
    if (Mpi::Root(Mpi::World()))
    {
      json ring_library = {
          {"Version", 3},
          {"TraceLiftVersion", 2},
          {"Name", "unit-test-curved-3d"},
          {"MatchingRadius", ring_R},
          {"CouponDepth", 0.2},
          {"Fabrication",
           {{"InterfaceLayers", {{"SA", {{"Thickness", 0.002}, {"Permittivity", 4.0}}}}}}},
          {"Models", json::array()}};
      json anchor = {{"Name", "isolated"},
                     {"Topology", "IsolatedEdge"},
                     {"CouponDepth", 0.2},
                     {"FabricatedMatrix", fabricated_path.string()},
                     {"ThinMatrix", thin_path.string()},
                     {"FabricatedSurfaceMatrix", fabricated_surface_path.string()},
                     {"ThinSurfaceMatrix", thin_surface_path.string()},
                     {"BasisPoints", points_path.string()},
                     {"Interfaces", {{{"Type", "SA"}, {"Coupon", 1}}}}};
      ring_library["Models"].push_back(anchor);
      for (const char *convexity : {"Convex", "Concave"})
      {
        for (const double kappa : {0.1, 0.25, 0.5, 0.8})
        {
          auto model = anchor;
          model["Name"] = std::string("curved-") + convexity + "-" + std::to_string(kappa);
          model["Topology"] = "CurvedEdge";
          model["Kappa"] = kappa;
          model["Convexity"] = convexity;
          model["CouponDepth"] = 2.0 * M_PI * ring_R / kappa;
          ring_library["Models"].push_back(model);
        }
      }
      std::ofstream output(ring_library_path);
      output << ring_library.dump(2) << "\n";
      auto convex_only = ring_library;
      convex_only["Models"].erase(convex_only["Models"].begin() + 5,
                                  convex_only["Models"].end());
      std::ofstream convex_only_output(ring_convex_only_path);
      convex_only_output << convex_only.dump(2) << "\n";
    }
    Mpi::Barrier(Mpi::World());
    auto ring_config = island_config;
    auto &ring_correction = ring_config["Solver"]["Electrostatic"]["ResponseCorrection"];
    ring_correction.erase("PatchConstruction");
    ring_correction["UnmatchedPolicy"] = "Warn";
    const auto ring_manifest_path =
        temp.temp_dir / "surface-response-requirements-rings.json";
    const auto ring_patches_path = temp.temp_dir / "surface-response-patches.csv";
    auto RunRings = [&](const fs::path &library)
    {
      ring_correction["Library"] = library.string();
      IoData ring_iodata(ring_config, false);
      ring_iodata.boundaries.cracked_attributes.insert(9);
      fs::remove(ring_patches_path);
      WriteSurfaceResponseRequirements(ring_iodata, *ring_mesh,
                                       ring_manifest_path.string());
      Mpi::Barrier(Mpi::World());
      std::ifstream input(ring_manifest_path);
      REQUIRE(input);
      return json::parse(input);
    };
    // Per feature: the sum over its patch rows of ModelWeight x QuadratureWeight per
    // (Segment, S0, S1) interval and the set of models with their ModelWeight range.
    auto PatchSummary = [&](int feature_id)
    {
      std::map<std::tuple<int, double, double>, double> interval_sums;
      std::map<std::string, std::pair<double, double>> model_weights;  // min, max
      for (const auto &row : ReadPatchRows(ring_patches_path))
      {
        if (std::stoi(row[1]) != feature_id)
        {
          continue;
        }
        const double model_weight = std::stod(row[6]), quadrature = std::stod(row[7]);
        interval_sums[{std::stoi(row[10]), std::stod(row[11]), std::stod(row[12])}] +=
            model_weight * quadrature;
        auto [it, inserted] =
            model_weights.emplace(row[3], std::make_pair(model_weight, model_weight));
        if (!inserted)
        {
          it->second.first = std::min(it->second.first, model_weight);
          it->second.second = std::max(it->second.second, model_weight);
        }
      }
      return std::make_pair(interval_sums, model_weights);
    };
    {
      const json manifest = RunRings(ring_library_path);
      const auto &features = manifest["Identification"]["Features"];
      int curved = 0, straight_like = 0;
      for (const auto &feature : features)
      {
        const std::string type = feature["Type"].get<std::string>();
        if (type != "CurvedEdge" && type != "IsolatedEdge")
        {
          continue;
        }
        REQUIRE(feature["Match"]["Status"] == "Matched");
        const auto [interval_sums, model_weights] = PatchSummary(feature["Id"].get<int>());
        REQUIRE(!interval_sums.empty());
        for (const auto &[interval, sum] : interval_sums)
        {
          (void)interval;
          CHECK_THAT(sum, WithinAbs(1.0, 1.0e-9));
        }
        if (type == "CurvedEdge")
        {
          curved++;
          const double radius = feature["Signature"]["RadiusOverR"].get<double>();
          const std::string convexity =
              feature["Signature"]["Convexity"].get<std::string>();
          const bool inner = radius < 4.5;
          CHECK_THAT(radius, WithinAbs(inner ? 3.0 : 6.0, 1.0e-3));
          CHECK(convexity == (inner ? "Concave" : "Convex"));
          const std::string model = feature["Match"]["Model"].get<std::string>();
          CHECK(model.rfind("isolated@" + std::string(inner ? "concave" : "convex") +
                                "-kappa",
                            0) == 0);
          CHECK(model.find("-cubic") != std::string::npos);
          CHECK(feature["Match"]["Note"].get<std::string>().find("cubic") !=
                std::string::npos);
          REQUIRE(model_weights.size() == 1);
          CHECK(model_weights.begin()->first == model);
          CHECK_THAT(model_weights.begin()->second.second, WithinAbs(1.0, 1.0e-12));
          CHECK_THAT(feature["TurnTowardMetal"].get<double>(),
                     WithinRel((inner ? -1.0 : 1.0) * 2.0 * M_PI, 0.02));
        }
        else
        {
          REQUIRE(!feature["BendRadiusOverR"].is_null());
          const double bend = feature["BendRadiusOverR"].get<double>();
          if (bend > 100.0)
          {
            continue;  // no bend of note
          }
          straight_like++;
          const bool inner = bend < 21.5;
          CHECK_THAT(bend, WithinAbs(inner ? 20.0 : 23.0, 1.0e-3));
          CHECK(feature["Match"]["Model"] == "isolated");
          CHECK_THAT(feature["TurnTowardMetal"].get<double>(),
                     WithinRel((inner ? -1.0 : 1.0) * 2.0 * M_PI, 0.02));
          // Anchor + the node of the turn's convexity at a = (R / rho) / 0.1.
          const double a = (inner ? 1.0 / 20.0 : 1.0 / 23.0) / 0.1;
          REQUIRE(model_weights.size() == 2);
          const std::string node = std::string("curved-") + (inner ? "Concave" : "Convex") +
                                   "-" + std::to_string(0.1);
          REQUIRE(model_weights.count("isolated") == 1);
          REQUIRE(model_weights.count(node) == 1);
          CHECK_THAT(model_weights.at(node).first, WithinAbs(a, 0.02));
          CHECK_THAT(model_weights.at(node).second, WithinAbs(a, 0.02));
          CHECK_THAT(model_weights.at("isolated").first, WithinAbs(1.0 - a, 0.02));
        }
      }
      CHECK(curved == 2);
      CHECK(straight_like == 2);
      // The version-1 records: the curved groups Interpolated with their family.
      int interpolated = 0;
      for (const auto &record : manifest["Requirements"])
      {
        if (record["Topology"] != "CurvedEdge")
        {
          continue;
        }
        CHECK(record["Status"] == "Interpolated");
        REQUIRE(record.contains("CurvatureFamily"));
        CHECK(record["CurvatureFamily"]["InterpolationRule"] == "cubic");
        double weights = 0.0;
        for (const auto &node : record["CurvatureFamily"]["Nodes"])
        {
          weights += node["Weight"].get<double>();
        }
        CHECK_THAT(weights, WithinAbs(1.0, 1.0e-9));
        interpolated++;
      }
      CHECK(interpolated == 2);
    }
    {
      // Without concave coupons: the inner CurvedEdge is unmatched with the reason and the
      // inner straight-like edge keeps the anchor alone.
      const json manifest = RunRings(ring_convex_only_path);
      int unmatched = 0, anchor_only = 0;
      for (const auto &feature : manifest["Identification"]["Features"])
      {
        if (feature["Type"] == "CurvedEdge" &&
            feature["Signature"]["Convexity"] == "Concave")
        {
          CHECK(feature["Match"]["Status"] == "Missing");
          CHECK(feature["Match"]["Note"].get<std::string>().find("no concave") !=
                std::string::npos);
          unmatched++;
        }
        if (feature["Type"] == "IsolatedEdge" && !feature["BendRadiusOverR"].is_null() &&
            feature["BendRadiusOverR"].get<double>() < 21.5)
        {
          const auto [interval_sums, model_weights] =
              PatchSummary(feature["Id"].get<int>());
          REQUIRE(model_weights.size() == 1);
          CHECK(model_weights.begin()->first == "isolated");
          CHECK_THAT(model_weights.begin()->second.first, WithinAbs(1.0, 1.0e-12));
          anchor_only++;
        }
      }
      CHECK(unmatched == 1);
      CHECK(anchor_only == 1);
    }
  }

#endif
}

TEST_CASE_METHOD(test::SurfaceResponseFiles, "SurfaceResponseOperator corner angle family",
                 "[surfaceresponseoperator][3d][corner][features][Serial][Parallel]")
{
#if !defined(MFEM_USE_GSLIB)
  SKIP("SurfaceResponseOperator requires MFEM_USE_GSLIB");
#else
  json island_config = IslandConfig();
  // Angle-interpolated corner family on the FEATURES path (USER decision 121 (C), the
  // kink-aware stencils of the corner-qualification block 2026-09-29, decisions 136 / 137):
  // a hexagonal island with vertices (+-L, 0), (+-0.8 L, +-L) on one plane (L = 1 = 5 R)
  // has two convex corners of 157.38 deg (turn 22.62 deg, at (+-L, 0)) and four of
  // 101.31 deg (turn 78.69 deg). The synthetic family is built on one box basis with
  // SEGMENT connectivity records (TraceBasis ConnectivityAngleDegrees): the segment
  // [90, 135] keyed 112.5 (nodes 90 / 105 / 120 / 135) and the segment [153.435, 180] keyed
  // 166.7 (nodes 153.434948822922 = the knot-corner passage / 160 / 165 / the straight
  // anchor 180), so both corner classes are cubic Lagrange on the four nodes of their
  // segment (weights summing to one, never across a knot-corner passage of the trace
  // basis); the runtime models are named <base>@corner-angle<deg>-cubic on the nearest
  // node of the stencil (105 for 101.31, 160 for 157.38); every corner is patched once
  // (weight one) in its frame.
  // Without the 90 deg coupon the 101 deg corners are sharper than the sharpest coupon:
  // unmatched with the reason (never the nearest node); without the anchor the 157 deg
  // corners are quadratic on the three remaining nodes of their segment (the highest order
  // the segment supports, never across the passage); without the 160 / 165 / 180 nodes they
  // are wider than the widest coupon: unmatched with the reason (no extrapolation) while
  // the 101 deg corners keep their segment.
  {
    constexpr double hex_R = 0.2;
    constexpr double hex_L = 1.0;
    constexpr double hex_g = 0.8;
    constexpr double hex_passage = 153.434948822922;  // convex knot-corner passage (free 2)
    auto MakeHexagonMesh = [&](double extent, double h)
    {
      const int n = static_cast<int>(std::lround(extent / h));
      mfem::Mesh serial = mfem::Mesh::MakeCartesian3D(n, 4, n, mfem::Element::HEXAHEDRON,
                                                      extent, 1.0, extent);
      const double c = 0.5 * extent;
      auto Level = [&](const double *point)
      { return std::max(std::abs(point[0] - c), std::abs(point[2] - c)); };
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
        double level_max = 0.0;
        for (const int vertex : vertices)
        {
          const double *point = serial.GetVertex(vertex);
          on_plane = on_plane && std::abs(point[1] - 0.5) < 1.0e-12;
          level_max = std::max(level_max, Level(point));
        }
        if (on_plane && level_max <= hex_L + 1.0e-9)
        {
          serial.AddBdrElement(serial.GetFace(face)->Duplicate(&serial));
          serial.SetBdrAttribute(serial.GetNBE() - 1, 9);
        }
      }
      // Radial map of every concentric square of the (x, z) grid onto the hexagon with
      // vertices (+-1, 0), (+-g, +-1) at the same level: the square's boundary points land
      // on the hexagon's sides (collinear between its vertices, which are grid points), so
      // the island outline is the exact hexagon at every level.
      const std::vector<std::array<double, 2>> hexagon = {{1.0, 0.0},     {hex_g, 1.0},
                                                          {-hex_g, 1.0},  {-1.0, 0.0},
                                                          {-hex_g, -1.0}, {hex_g, -1.0}};
      auto PolygonRadius = [&](double ux, double uz)
      {
        double r = mfem::infinity();
        for (std::size_t k = 0; k < hexagon.size(); k++)
        {
          const auto &a = hexagon[k], &b = hexagon[(k + 1) % hexagon.size()];
          const double nx = b[1] - a[1], nz = a[0] - b[0];  // outward normal of side a -> b
          const double denominator = nx * ux + nz * uz;
          if (denominator > 1.0e-14)
          {
            r = std::min(r, (nx * a[0] + nz * a[1]) / denominator);
          }
        }
        return r;
      };
      // The map is the hexagon map up to level 2 L, blends linearly to the identity at
      // level 3 L and leaves the outer grid (the PEC box walls stay on the bounding box, so
      // that the island's plane remains the process plane) unchanged.
      for (int vertex = 0; vertex < serial.GetNV(); vertex++)
      {
        double *point = serial.GetVertex(vertex);
        const double lx = point[0] - c, lz = point[2] - c;
        const double norm = std::hypot(lx, lz);
        const double level = std::max(std::abs(lx), std::abs(lz));
        if (norm > 1.0e-12 && level < 3.0 * hex_L - 1.0e-9)
        {
          const double ux = lx / norm, uz = lz / norm;
          const double square = 1.0 / std::max(std::abs(ux), std::abs(uz));
          const double hexagon_scale = PolygonRadius(ux, uz) / square;
          const double blend = level <= 2.0 * hex_L ? 1.0 : (3.0 * hex_L - level) / hex_L;
          const double scale = 1.0 + blend * (hexagon_scale - 1.0);
          point[0] = c + scale * lx;
          point[2] = c + scale * lz;
        }
      }
      serial.FinalizeTopology();
      serial.Finalize();
      return std::make_unique<mfem::ParMesh>(Mpi::World(), serial);
    };
    // Extent 8: the hexagon (level 1) stays more than 2 R from the PEC box walls.
    auto hexagon_mesh = MakeHexagonMesh(8.0, 0.1);
    const auto corner_family_path =
        temp.temp_dir / "fabrication-process-corner-family.json";
    const auto corner_family_no90_path =
        temp.temp_dir / "fabrication-process-corner-family-no90.json";
    const auto corner_family_no_anchor_path =
        temp.temp_dir / "fabrication-process-corner-family-no-anchor.json";
    const auto corner_family_no_wide_path =
        temp.temp_dir / "fabrication-process-corner-family-no-wide.json";
    if (Mpi::Root(Mpi::World()))
    {
      json family_library = {
          {"Version", 3},
          {"TraceLiftVersion", 2},
          {"Name", "unit-test-corner-family-3d"},
          {"MatchingRadius", hex_R},
          {"Fabrication",
           {{"InterfaceLayers", {{"SA", {{"Thickness", 0.002}, {"Permittivity", 4.0}}}}}}},
          {"Models", json::array()}};
      json isolated = {{"Name", "isolated"},
                       {"Topology", "IsolatedEdge"},
                       {"CouponDepth", 0.2},
                       {"FabricatedMatrix", fabricated_path.string()},
                       {"ThinMatrix", thin_path.string()},
                       {"FabricatedSurfaceMatrix", fabricated_surface_path.string()},
                       {"ThinSurfaceMatrix", thin_surface_path.string()},
                       {"BasisPoints", points_path.string()},
                       {"Interfaces", {{{"Type", "SA"}, {"Coupon", 1}}}}};
      family_library["Models"].push_back(isolated);
      // The family's coupons are built on the trace basis rule (corner-family review
      // 2026-09-29): one knot semantics for every node, 72 knots, the same zero set; each
      // with its segment's connectivity (one triangulation per segment).
      const auto family_fabricated_path = temp.temp_dir / "corner-family-fabricated.csv";
      const auto family_thin_path = temp.temp_dir / "corner-family-thin.csv";
      const auto family_fabricated_surface_path =
          temp.temp_dir / "corner-family-fabricated-surface.csv";
      const auto family_thin_surface_path =
          temp.temp_dir / "corner-family-thin-surface.csv";
      WriteCornerMatrices(family_fabricated_path, family_fabricated_surface_path, 72, 3.0,
                          0.05, hex_R);
      WriteCornerMatrices(family_thin_path, family_thin_surface_path, 72, 1.0, 0.01, hex_R);
      const std::vector<std::pair<double, double>> family_nodes = {
          {90.0, 112.5},        {105.0, 112.5}, {120.0, 112.5}, {135.0, 112.5},
          {hex_passage, 166.7}, {160.0, 166.7}, {165.0, 166.7}, {180.0, 166.7}};
      for (const auto &[angle, connectivity] : family_nodes)
      {
        std::ostringstream tag;
        tag << std::setprecision(12) << angle;
        const auto files =
            WriteCornerBasisFiles(temp.temp_dir, "family-" + tag.str(), angle, true, hex_R,
                                  0.05 * hex_R, 0.025 * hex_R, true, connectivity);
        json model = {{"Name", "convex-corner-" + tag.str()},
                      {"Topology", "ConvexCorner"},
                      {"Angle", angle},
                      {"AngleDegrees", angle},
                      {"Convexity", "Convex"},
                      {"AngleTolerance", 1.0e-6},
                      {"CornerRadius", 0.0},
                      {"CornerRadiusTolerance", 0.0},
                      {"FabricatedMatrix", family_fabricated_path.string()},
                      {"ThinMatrix", family_thin_path.string()},
                      {"FabricatedSurfaceMatrix", family_fabricated_surface_path.string()},
                      {"ThinSurfaceMatrix", family_thin_surface_path.string()},
                      {"BasisPoints", files.points.string()},
                      {"TraceMesh",
                       {{"Vertices", files.vertices.string()},
                        {"Triangles", files.triangles.string()}}},
                      {"ContourGroups", files.contour_groups},
                      {"ZeroTraceIndices", files.zero_trace_indices},
                      {"TraceBasis", files.trace_basis},
                      {"Interfaces", {{{"Type", "SA"}, {"Coupon", 1}}}}};
        family_library["Models"].push_back(model);
      }
      std::ofstream output(corner_family_path);
      output << family_library.dump(2) << "\n";
      auto no90 = family_library;
      no90["Name"] = "unit-test-corner-family-no90-3d";
      no90["Models"].erase(no90["Models"].begin() + 1);
      std::ofstream no90_output(corner_family_no90_path);
      no90_output << no90.dump(2) << "\n";
      auto no_anchor = family_library;
      no_anchor["Name"] = "unit-test-corner-family-no-anchor-3d";
      no_anchor["Models"].erase(no_anchor["Models"].end() - 1);
      std::ofstream no_anchor_output(corner_family_no_anchor_path);
      no_anchor_output << no_anchor.dump(2) << "\n";
      // Without the 160 / 165 / 180 nodes the wide segment is the passage node alone.
      auto no_wide = family_library;
      no_wide["Name"] = "unit-test-corner-family-no-wide-3d";
      no_wide["Models"].erase(no_wide["Models"].end() - 3, no_wide["Models"].end());
      std::ofstream no_wide_output(corner_family_no_wide_path);
      no_wide_output << no_wide.dump(2) << "\n";
    }
    Mpi::Barrier(Mpi::World());
    auto hexagon_config = island_config;
    auto &hexagon_correction =
        hexagon_config["Solver"]["Electrostatic"]["ResponseCorrection"];
    hexagon_correction.erase("PatchConstruction");
    hexagon_correction["UnmatchedPolicy"] = "Warn";
    const auto hexagon_manifest_path =
        temp.temp_dir / "surface-response-requirements-hexagon.json";
    const auto hexagon_patches_path = temp.temp_dir / "surface-response-patches.csv";
    auto RunHexagon = [&](const fs::path &library)
    {
      hexagon_correction["Library"] = library.string();
      IoData hexagon_iodata(hexagon_config, false);
      hexagon_iodata.boundaries.cracked_attributes.insert(9);
      fs::remove(hexagon_patches_path);
      WriteSurfaceResponseRequirements(hexagon_iodata, *hexagon_mesh,
                                       hexagon_manifest_path.string());
      Mpi::Barrier(Mpi::World());
      std::ifstream input(hexagon_manifest_path);
      REQUIRE(input);
      return json::parse(input);
    };
    const double wide_angle =
        180.0 - 2.0 * std::atan(1.0 - hex_g) * 180.0 / M_PI;                   // 157.38
    const double narrow_angle = 90.0 + std::atan(1.0 - hex_g) * 180.0 / M_PI;  // 101.31
    {
      const json manifest = RunHexagon(corner_family_path);
      int wide = 0, narrow = 0;
      for (const auto &feature : manifest["Identification"]["Features"])
      {
        if (feature["Type"] != "ConvexCorner")
        {
          continue;
        }
        const double angle = feature["Signature"]["AngleDegrees"].get<double>();
        CHECK(feature["Signature"]["CornerRadiusOverR"].get<double>() == 0.0);
        CHECK(feature["Match"]["Status"] == "Matched");
        const std::string model = feature["Match"]["Model"].get<std::string>();
        if (std::abs(angle - wide_angle) < 1.0e-6)
        {
          // Cubic on the segment [153.435, 180] keyed 166.7; the base is the nearest node
          // (160: 2.62 deg away, the passage node 3.95).
          CHECK(model.find("convex-corner-160@corner-angle") != std::string::npos);
          CHECK(model.find("-cubic") != std::string::npos);
          wide++;
        }
        else
        {
          // Cubic on the segment [90, 135] keyed 112.5; the base is the nearest node (105).
          CHECK_THAT(angle, WithinAbs(narrow_angle, 1.0e-6));
          CHECK(model.find("convex-corner-105@corner-angle") != std::string::npos);
          CHECK(model.find("-cubic") != std::string::npos);
          narrow++;
        }
        // One patch of weight one per corner.
        int patches = 0;
        for (const auto &row : ReadPatchRows(hexagon_patches_path))
        {
          if (std::stoi(row[1]) == feature["Id"].get<int>())
          {
            CHECK(row[3] == model);
            CHECK_THAT(std::stod(row[6]), WithinAbs(1.0, 1.0e-12));
            patches++;
          }
        }
        CHECK(patches == 1);
      }
      CHECK(wide == 2);
      CHECK(narrow == 4);
      // The version-1 records: Interpolated with the family selection and its weights.
      int records = 0;
      for (const auto &record : manifest["Requirements"])
      {
        if (record["Topology"] != "ConvexCorner")
        {
          continue;
        }
        CHECK(record["Status"] == "Interpolated");
        REQUIRE(record.contains("CornerFamily"));
        const auto &family = record["CornerFamily"];
        const double angle = family["AngleDegrees"].get<double>();
        CHECK(family["Convexity"] == "Convex");
        CHECK_THAT(family["MaxTurnDegrees"].get<double>(), WithinAbs(90.0, 1.0e-9));
        // The smallest node turn is the 165 deg node's (the first-order record field; the
        // stencil itself is the segment's cubic).
        CHECK_THAT(family["FirstOrderTurnDegrees"].get<double>(), WithinAbs(15.0, 1.0e-9));
        // (The node angles are the models' Angle in radians back in degrees, the group's
        // AngleDegrees the representative on the 1e-6 deg signature grid: keys rounded,
        // weights compared at 1e-6.)
        std::map<double, double> weights;
        for (const auto &node : family["Nodes"])
        {
          weights[std::round(node["AngleDegrees"].get<double>() * 1.0e6) * 1.0e-6] =
              node["Weight"].get<double>();
        }
        double sum = 0.0;
        for (const auto &[node_angle, weight] : weights)
        {
          (void)node_angle;
          sum += weight;
        }
        CHECK_THAT(sum, WithinAbs(1.0, 1.0e-9));
        // Both classes: cubic Lagrange on the four nodes of the angle's segment (the
        // segment's connectivity angle recorded), weights the Lagrange polynomials in the
        // turn.
        CHECK(family["InterpolationRule"] == "cubic");
        REQUIRE(weights.size() == 4);
        const bool wide_record = std::abs(angle - wide_angle) < 1.0e-6;
        CHECK_THAT(family["ConnectivityAngleDegrees"].get<double>(),
                   WithinAbs(wide_record ? 166.7 : 112.5, 1.0e-9));
        const std::vector<double> expected_nodes =
            wide_record ? std::vector<double>{hex_passage, 160.0, 165.0, 180.0}
                        : std::vector<double>{90.0, 105.0, 120.0, 135.0};
        for (const double node_angle : expected_nodes)
        {
          CHECK(weights.count(std::round(node_angle * 1.0e6) * 1.0e-6) == 1);
        }
        {
          const double t = 180.0 - (wide_record ? wide_angle : narrow_angle);
          for (const auto &[node_angle, weight] : weights)
          {
            double expected = 1.0;
            const double ti = 180.0 - node_angle;
            for (const auto &[other_angle, other_weight] : weights)
            {
              (void)other_weight;
              const double tj = 180.0 - other_angle;
              if (other_angle != node_angle)
              {
                expected *= (t - tj) / (ti - tj);
              }
            }
            CHECK_THAT(weight, WithinAbs(expected, 1.0e-6));
          }
        }
        records++;
      }
      CHECK(records == 2);
    }
    {
      // Without the 90 deg coupon the 101 deg corners are sharper than the sharpest coupon
      // (105 deg): unmatched with the reason; the 157 deg corners keep their cubic segment.
      const json manifest = RunHexagon(corner_family_no90_path);
      int refused = 0, cubic = 0;
      for (const auto &feature : manifest["Identification"]["Features"])
      {
        if (feature["Type"] != "ConvexCorner")
        {
          continue;
        }
        if (feature["Match"]["Status"] == "Missing")
        {
          CHECK(feature["Match"]["Note"].get<std::string>().find("sharper than") !=
                std::string::npos);
          refused++;
        }
        else
        {
          CHECK(feature["Match"]["Model"].get<std::string>().find(
                    "convex-corner-160@corner-angle") != std::string::npos);
          CHECK(feature["Match"]["Model"].get<std::string>().find("-cubic") !=
                std::string::npos);
          cubic++;
        }
      }
      CHECK(refused == 4);
      CHECK(cubic == 2);
    }
    {
      // Without the straight anchor the 157 deg corners are quadratic on the three
      // remaining nodes of their segment (153.435 / 160 / 165 keyed 166.7, base 160: the
      // highest order the segment supports); the 101 deg corners keep their cubic segment
      // [90, 135].
      const json manifest = RunHexagon(corner_family_no_anchor_path);
      int quadratic = 0, cubic = 0;
      for (const auto &feature : manifest["Identification"]["Features"])
      {
        if (feature["Type"] != "ConvexCorner")
        {
          continue;
        }
        REQUIRE(feature["Match"]["Status"] == "Matched");
        const std::string model = feature["Match"]["Model"].get<std::string>();
        if (model.find("-quadratic") != std::string::npos)
        {
          CHECK(model.find("convex-corner-160@corner-angle") != std::string::npos);
          quadratic++;
        }
        else
        {
          CHECK(model.find("convex-corner-105@corner-angle") != std::string::npos);
          CHECK(model.find("-cubic") != std::string::npos);
          cubic++;
        }
      }
      CHECK(quadratic == 2);
      CHECK(cubic == 4);
    }
    {
      // Without the 160 / 165 / 180 nodes the 157 deg corners are wider than the widest
      // coupon (the passage node 153.435): unmatched with the reason (no extrapolation);
      // the 101 deg corners keep their cubic segment [90, 135].
      const json manifest = RunHexagon(corner_family_no_wide_path);
      int refused = 0, cubic = 0;
      for (const auto &feature : manifest["Identification"]["Features"])
      {
        if (feature["Type"] != "ConvexCorner")
        {
          continue;
        }
        if (feature["Match"]["Status"] == "Missing")
        {
          CHECK(feature["Match"]["Note"].get<std::string>().find("wider than") !=
                std::string::npos);
          refused++;
        }
        else
        {
          CHECK(feature["Match"]["Model"].get<std::string>().find(
                    "convex-corner-105@corner-angle") != std::string::npos);
          CHECK(feature["Match"]["Model"].get<std::string>().find("-cubic") !=
                std::string::npos);
          cubic++;
        }
      }
      CHECK(refused == 2);
      CHECK(cubic == 4);
    }
  }

#endif
}

TEST_CASE_METHOD(test::SurfaceResponseFiles, "SurfaceResponseOperator Maxwell islands",
                 "[surfaceresponseoperator][3d][maxwell][corner][Serial][Parallel]")
{
#if !defined(MFEM_USE_GSLIB)
  SKIP("SurfaceResponseOperator requires MFEM_USE_GSLIB");
#else
  json island_config = IslandConfig();
  json convex_island_config = ConvexIslandConfig();
  IoData convex_island_iodata(convex_island_config, false);
  convex_island_iodata.boundaries.cracked_attributes.insert(9);
  const auto &island_line_rule = mfem::IntRules.Get(mfem::Geometry::SEGMENT, 2);
  // The island's straight segments and corners (asserted by the island corners case) and
  // the electrostatic patch counts the Maxwell operators reproduce.
  const auto island_edges =
      SummarizeInterfaceEdges(*MakeIslandMesh(), convex_island_iodata.boundaries, 4);
  const auto &island_segments = island_edges.segments;
  const int island_corners = island_edges.corners;
  const int removed_straight_patches = 2 * island_corners * island_line_rule.GetNPoints();
  const auto touching_segments =
      SummarizeInterfaceEdges(*MakeTouchingIslandMesh(), convex_island_iodata.boundaries, 4)
          .segments;
  auto junction_island_config = convex_island_config;
  junction_island_config["Solver"]["Electrostatic"]["ResponseCorrection"]["Library"] =
      junction_library_3d_path.string();
  IoData junction_island_iodata(junction_island_config, false);
  junction_island_iodata.boundaries.cracked_attributes.insert(9);
  std::vector<std::unique_ptr<Mesh>> junction_island_meshes;
  junction_island_meshes.push_back(std::make_unique<Mesh>(MakeTouchingIslandMesh()));
  LaplaceOperator junction_island_laplace(junction_island_iodata, junction_island_meshes);
  SurfaceResponseOperator junction_island_response(junction_island_iodata,
                                                   junction_island_laplace);
  json rounded_island_config = RoundedIslandConfig();
  constexpr int rounded_corner_count = 4;
  mfem::VectorConstantCoefficient field_coefficient = ConstantFieldCoefficient();

  json convex_maxwell_island_config = ConvexMaxwellIslandConfig();
  IoData convex_maxwell_island_iodata(convex_maxwell_island_config, false);
  convex_maxwell_island_iodata.boundaries.cracked_attributes.insert(9);
  std::vector<std::unique_ptr<Mesh>> convex_maxwell_island_meshes;
  convex_maxwell_island_meshes.push_back(std::make_unique<Mesh>(MakeIslandMesh()));
  SpaceOperator convex_maxwell_island_space(convex_maxwell_island_iodata,
                                            convex_maxwell_island_meshes);
  SurfaceResponseOperator convex_maxwell_island_response(convex_maxwell_island_iodata,
                                                         convex_maxwell_island_space);
  CHECK(convex_maxwell_island_response.GetPatchCount() ==
        static_cast<int>(island_segments.size()) * island_line_rule.GetNPoints() -
            removed_straight_patches + island_corners);

  GridFunction island_field(convex_maxwell_island_space.GetNDSpace(), true);
  island_field.Real().ProjectCoefficient(field_coefficient);
  island_field.Imag() = 0.0;
  const auto constant_island_response =
      convex_maxwell_island_response.GetMaxwellResponse(island_field, 0.0);
  CHECK(constant_island_response.loop_residual < 1.0e-10);
  CHECK(constant_island_response.corner_neighborhood_fraction == 0.0);
  CHECK(constant_island_response.closure_independent_confident);
  CHECK(constant_island_response.maximum_trace_closure_spread > 0.05);
  CHECK(constant_island_response.response_weighted_trace_closure_spread > 0.0);
  CHECK(constant_island_response.trace_closure_response_failure_fraction > 0.0);
  CHECK_FALSE(constant_island_response.confident);

  auto impedance_convex_maxwell_config = convex_maxwell_island_config;
  impedance_convex_maxwell_config["Boundaries"]["Ground"]["Attributes"] = {1, 2, 3,
                                                                           4, 5, 6};
  impedance_convex_maxwell_config["Boundaries"]["Impedance"] = {
      {{"Attributes", {9}}, {"Ls", 1.0e-13}}};
  impedance_convex_maxwell_config["Solver"]["SurfaceResponseCorrection"]["Library"] =
      finite_impedance_convex_library_3d_path.string();
  IoData impedance_convex_maxwell_iodata(impedance_convex_maxwell_config, false);
  impedance_convex_maxwell_iodata.boundaries.cracked_attributes.insert(9);
  std::vector<std::unique_ptr<Mesh>> impedance_convex_maxwell_meshes;
  impedance_convex_maxwell_meshes.push_back(std::make_unique<Mesh>(MakeIslandMesh()));
  SpaceOperator impedance_convex_maxwell_space(impedance_convex_maxwell_iodata,
                                               impedance_convex_maxwell_meshes);
  SurfaceResponseOperator impedance_convex_maxwell_response(impedance_convex_maxwell_iodata,
                                                            impedance_convex_maxwell_space);
  CHECK(impedance_convex_maxwell_response.GetPatchCount() ==
        convex_maxwell_island_response.GetPatchCount());
  GridFunction impedance_convex_maxwell_field(impedance_convex_maxwell_space.GetNDSpace(),
                                              true);
  impedance_convex_maxwell_field.Real().ProjectCoefficient(field_coefficient);
  impedance_convex_maxwell_field.Imag() = 0.0;
  const auto impedance_convex_maxwell_result =
      impedance_convex_maxwell_response.GetMaxwellResponse(impedance_convex_maxwell_field,
                                                           0.0);
  CHECK(impedance_convex_maxwell_result.loop_residual < 1.0e-10);
  CHECK_THAT(impedance_convex_maxwell_result.matched_length_fraction,
             WithinAbs(1.0, 1.0e-12));
  CHECK(impedance_convex_maxwell_result.corner_neighborhood_fraction == 0.0);
  CHECK(impedance_convex_maxwell_result.boundary_law_verified);

  auto junction_maxwell_config = convex_maxwell_island_config;
  junction_maxwell_config["Solver"]["SurfaceResponseCorrection"]["Library"] =
      junction_library_3d_path.string();
  IoData junction_maxwell_iodata(junction_maxwell_config, false);
  junction_maxwell_iodata.boundaries.cracked_attributes.insert(9);
  std::vector<std::unique_ptr<Mesh>> junction_maxwell_meshes;
  junction_maxwell_meshes.push_back(std::make_unique<Mesh>(MakeTouchingIslandMesh()));
  SpaceOperator junction_maxwell_space(junction_maxwell_iodata, junction_maxwell_meshes);
  SurfaceResponseOperator junction_maxwell_response(junction_maxwell_iodata,
                                                    junction_maxwell_space);
  CHECK(junction_maxwell_response.GetPatchCount() ==
        junction_island_response.GetPatchCount());
  GridFunction junction_maxwell_field(junction_maxwell_space.GetNDSpace(), true);
  junction_maxwell_field.Real().ProjectCoefficient(field_coefficient);
  junction_maxwell_field.Imag() = 0.0;
  const auto junction_maxwell_result =
      junction_maxwell_response.GetMaxwellResponse(junction_maxwell_field, 0.0);
  CHECK(junction_maxwell_result.loop_residual < 1.0e-10);
  CHECK(junction_maxwell_result.corner_neighborhood_fraction == 0.0);

  const auto junction_requirements_path =
      temp.temp_dir / "surface-response-requirements-junction.json";
  WriteSurfaceResponseRequirements(junction_maxwell_iodata, *junction_maxwell_meshes.back(),
                                   junction_requirements_path.string());
  std::ifstream junction_requirements_input(junction_requirements_path);
  REQUIRE(junction_requirements_input);
  const json junction_requirements = json::parse(junction_requirements_input);
  // Phase 3 of the identification fix: the two squares touch at one vertex, a point
  // contact between two edge-connected metal components. The legacy classifier builds no
  // junction there (its arms must share one conductor); the version-2 identification
  // reports a Junction feature whose signature carries the arm conductors (two labels) and
  // flags the vertex as a PointContact, so the degenerate geometry is never silent.
  CHECK(std::none_of(junction_requirements["LegacyRequirements"].begin(),
                     junction_requirements["LegacyRequirements"].end(),
                     [](const auto &requirement)
                     { return requirement["Topology"] == "Junction"; }));
  const auto &junction_identification = junction_requirements["Identification"];
  const auto junction_feature =
      std::find_if(junction_identification["Features"].begin(),
                   junction_identification["Features"].end(),
                   [](const auto &feature) { return feature["Type"] == "Junction"; });
  REQUIRE(junction_feature != junction_identification["Features"].end());
  CHECK((*junction_feature)["Match"]["Status"] == "Missing");
  const auto &arm_conductors = (*junction_feature)["Signature"]["ArmConductors"];
  REQUIRE(arm_conductors.size() == 4);
  CHECK(std::set<int>(arm_conductors.begin(), arm_conductors.end()) == std::set<int>{1, 2});
  CHECK((*junction_feature)["Signature"]["ArmAnglesDegrees"].size() == 4);
  CHECK(std::count_if(junction_identification["Vertices"].begin(),
                      junction_identification["Vertices"].end(), [](const auto &vertex)
                      { return vertex.value("PointContact", false); }) == 1);

  auto impedance_junction_config = junction_maxwell_config;
  impedance_junction_config["Boundaries"]["Ground"]["Attributes"] = {1, 2, 3, 4, 5, 6};
  impedance_junction_config["Boundaries"]["Impedance"] = {
      {{"Attributes", {9}}, {"Ls", 1.0e-13}}};
  IoData impedance_junction_iodata(impedance_junction_config, false);
  impedance_junction_iodata.boundaries.cracked_attributes.insert(9);
  std::vector<std::unique_ptr<Mesh>> impedance_junction_meshes;
  impedance_junction_meshes.push_back(std::make_unique<Mesh>(MakeTouchingIslandMesh()));
  SpaceOperator impedance_junction_space(impedance_junction_iodata,
                                         impedance_junction_meshes);
  SurfaceResponseOperator impedance_junction_response(impedance_junction_iodata,
                                                      impedance_junction_space);
  // No legacy junction patch at the two-conductor point contact (see above).
  CHECK(impedance_junction_response.GetPatchCount() ==
        static_cast<int>(touching_segments.size()) * island_line_rule.GetNPoints());
  CHECK(impedance_junction_response.GetBasisSize() ==
        4 * static_cast<int>(touching_segments.size()) * island_line_rule.GetNPoints());
  GridFunction impedance_junction_field(impedance_junction_space.GetNDSpace(), true);
  impedance_junction_field.Real().ProjectCoefficient(field_coefficient);
  impedance_junction_field.Imag() = 0.0;
  const auto impedance_junction_result =
      impedance_junction_response.GetMaxwellResponse(impedance_junction_field, 0.0);
  CHECK(impedance_junction_result.loop_residual < 1.0e-10);

  // Separated fabrication planes are independent placements of the same local process.
  // They should both match unless their radius-R neighborhoods actually interact.
  auto multilayer_maxwell_config = convex_maxwell_island_config;
  multilayer_maxwell_config["Boundaries"]["Ground"]["Attributes"] = {1, 2, 3, 4,
                                                                     5, 6, 9, 10};
  auto second_interface =
      multilayer_maxwell_config["Boundaries"]["Postprocessing"]["Dielectric"][0];
  second_interface["Index"] = 5;
  second_interface["Attributes"] = {10};
  multilayer_maxwell_config["Boundaries"]["Postprocessing"]["Dielectric"].push_back(
      second_interface);
  multilayer_maxwell_config["Solver"]["SurfaceResponseCorrection"]["TargetInterfaces"] = {
      4, 5};
  IoData multilayer_maxwell_iodata(multilayer_maxwell_config, false);
  multilayer_maxwell_iodata.boundaries.cracked_attributes.insert(9);
  multilayer_maxwell_iodata.boundaries.cracked_attributes.insert(10);
  std::vector<std::unique_ptr<Mesh>> multilayer_maxwell_meshes;
  multilayer_maxwell_meshes.push_back(
      std::make_unique<Mesh>(MakeIslandMesh(false, false, false, false, true)));
  SpaceOperator multilayer_maxwell_space(multilayer_maxwell_iodata,
                                         multilayer_maxwell_meshes);
  SurfaceResponseOperator multilayer_maxwell_response(multilayer_maxwell_iodata,
                                                      multilayer_maxwell_space);
  CHECK(multilayer_maxwell_response.GetPatchCount() ==
        2 * convex_maxwell_island_response.GetPatchCount());
  CHECK(multilayer_maxwell_response.GetTargetInterfaces() == std::set<int>{4, 5});
  GridFunction multilayer_field(multilayer_maxwell_space.GetNDSpace(), true);
  multilayer_field.Real().ProjectCoefficient(field_coefficient);
  multilayer_field.Imag() = 0.0;
  const auto multilayer_result =
      multilayer_maxwell_response.GetMaxwellResponse(multilayer_field, 0.0);
  CHECK_THAT(multilayer_result.matched_length_fraction, WithinAbs(1.0, 1.0e-12));
  CHECK(multilayer_result.fabricated_surface_energy.count(4) == 1);
  CHECK(multilayer_result.fabricated_surface_energy.count(5) == 1);

  Vector island_true, island_correction, island_probe, island_probe_correction;
  island_field.Real().GetTrueDofs(island_true);
  auto *island_data = island_true.HostWrite();
  for (int i = 0; i < island_true.Size(); i++)
  {
    island_data[i] = std::cos(0.23 * (i + 1 + 7 * Mpi::Rank(Mpi::World())));
  }
  island_true.SetSubVector(convex_maxwell_island_space.GetNDDbcTDofLists().back(), 0.0);
  island_field.Real().SetFromTrueDofs(island_true);
  const auto random_island_response =
      convex_maxwell_island_response.GetMaxwellResponse(island_field, 0.0);
  convex_maxwell_island_response.Mult(island_true, island_correction);
  CHECK_THAT(0.5 * linalg::Dot(Mpi::World(), island_true, island_correction),
             WithinRel(random_island_response.domain_correction, 1.0e-10));
  island_probe.SetSize(island_true.Size());
  auto *island_probe_data = island_probe.HostWrite();
  for (int i = 0; i < island_probe.Size(); i++)
  {
    island_probe_data[i] = std::sin(0.31 * (i + 1 + 5 * Mpi::Rank(Mpi::World())));
  }
  island_probe.SetSubVector(convex_maxwell_island_space.GetNDDbcTDofLists().back(), 0.0);
  convex_maxwell_island_response.Mult(island_probe, island_probe_correction);
  CHECK_THAT(
      linalg::Dot(Mpi::World(), island_probe, island_correction),
      WithinRel(linalg::Dot(Mpi::World(), island_true, island_probe_correction), 1.0e-10));

  // A rounded-corner reference lies inside the PEC footprint. On a tetrahedral mesh,
  // exact line integration rejects an anchor path through that internal boundary. The
  // process-plane contour instead starts at a library-declared zero-trace knot.
  auto rounded_maxwell_island_config = rounded_island_config;
  rounded_maxwell_island_config["Problem"]["Type"] = "Eigenmode";
  rounded_maxwell_island_config["Boundaries"]["Ground"]["Attributes"] = {1, 2, 3, 4,
                                                                         5, 6, 9};
  rounded_maxwell_island_config["Boundaries"].erase("Terminal");
  rounded_maxwell_island_config["Solver"] = {
      {"Order", 1},
      {"Eigenmode", {{"Target", 1.0}}},
      {"SurfaceResponseCorrection",
       {{"Library", rounded_library_3d_path.string()},
        {"TargetInterfaces", {4}},
        {"UnmatchedPolicy", "Error"},
        {"PatchConstruction", "Legacy"}}}};
  IoData rounded_maxwell_island_iodata(rounded_maxwell_island_config, false);
  rounded_maxwell_island_iodata.boundaries.cracked_attributes.insert(9);
  std::vector<std::unique_ptr<Mesh>> rounded_maxwell_island_meshes;
  rounded_maxwell_island_meshes.push_back(
      std::make_unique<Mesh>(MakeIslandMesh(true, true)));
  SpaceOperator rounded_maxwell_island_space(rounded_maxwell_island_iodata,
                                             rounded_maxwell_island_meshes);
  SurfaceResponseOperator rounded_maxwell_island_response(rounded_maxwell_island_iodata,
                                                          rounded_maxwell_island_space);
  CHECK(rounded_maxwell_island_response.GetPatchCount() > rounded_corner_count);

  GridFunction rounded_island_field(rounded_maxwell_island_space.GetNDSpace(), true);
  rounded_island_field.Real().ProjectCoefficient(field_coefficient);
  rounded_island_field.Imag() = 0.0;
  const auto rounded_island_result =
      rounded_maxwell_island_response.GetMaxwellResponse(rounded_island_field, 0.0);
  CHECK(rounded_island_result.loop_residual < 1.0e-10);

  // Fixed-flux closure acts only on the free trace subspace. Matrix rows associated
  // with exact PEC trace knots are calibration artifacts and must not change either
  // closure when their free-free blocks are unchanged.
  auto constrained_perturbed_rounded_config = rounded_maxwell_island_config;
  constrained_perturbed_rounded_config["Solver"]["SurfaceResponseCorrection"]["Library"] =
      constrained_perturbed_rounded_library_3d_path.string();
  IoData constrained_perturbed_rounded_iodata(constrained_perturbed_rounded_config, false);
  constrained_perturbed_rounded_iodata.boundaries.cracked_attributes.insert(9);
  SurfaceResponseOperator constrained_perturbed_rounded_response(
      constrained_perturbed_rounded_iodata, rounded_maxwell_island_space);
  const auto constrained_perturbed_rounded_result =
      constrained_perturbed_rounded_response.GetMaxwellResponse(rounded_island_field, 0.0);
  CHECK_THAT(constrained_perturbed_rounded_result.domain_correction,
             WithinRel(rounded_island_result.domain_correction, 1.0e-12));
  CHECK_THAT(constrained_perturbed_rounded_result.domain_correction_fixed_flux,
             WithinRel(rounded_island_result.domain_correction_fixed_flux, 1.0e-12));
  CHECK_THAT(constrained_perturbed_rounded_result.fabricated_surface_energy.at(4),
             WithinRel(rounded_island_result.fabricated_surface_energy.at(4), 1.0e-12));
  CHECK_THAT(
      constrained_perturbed_rounded_result.fabricated_surface_energy_fixed_flux.at(4),
      WithinRel(rounded_island_result.fabricated_surface_energy_fixed_flux.at(4), 1.0e-12));
  CHECK_THAT(
      constrained_perturbed_rounded_result.response_weighted_trace_closure_spread,
      WithinRel(rounded_island_result.response_weighted_trace_closure_spread, 1.0e-12));
  CHECK_THAT(
      constrained_perturbed_rounded_result.trace_closure_response_failure_fraction,
      WithinAbs(rounded_island_result.trace_closure_response_failure_fraction, 1.0e-12));

  auto impedance_rounded_maxwell_config = rounded_maxwell_island_config;
  impedance_rounded_maxwell_config["Boundaries"]["Ground"]["Attributes"] = {1, 2, 3,
                                                                            4, 5, 6};
  impedance_rounded_maxwell_config["Boundaries"]["Impedance"] = {
      {{"Attributes", {9}}, {"Ls", 1.0e-13}}};
  impedance_rounded_maxwell_config["Solver"]["SurfaceResponseCorrection"]["Library"] =
      finite_impedance_rounded_library_3d_path.string();
  IoData impedance_rounded_maxwell_iodata(impedance_rounded_maxwell_config, false);
  impedance_rounded_maxwell_iodata.boundaries.cracked_attributes.insert(9);
  std::vector<std::unique_ptr<Mesh>> impedance_rounded_maxwell_meshes;
  impedance_rounded_maxwell_meshes.push_back(
      std::make_unique<Mesh>(MakeIslandMesh(true, true)));
  SpaceOperator impedance_rounded_maxwell_space(impedance_rounded_maxwell_iodata,
                                                impedance_rounded_maxwell_meshes);
  SurfaceResponseOperator impedance_rounded_maxwell_response(
      impedance_rounded_maxwell_iodata, impedance_rounded_maxwell_space);
  CHECK(impedance_rounded_maxwell_response.GetPatchCount() ==
        rounded_maxwell_island_response.GetPatchCount());
  GridFunction impedance_rounded_maxwell_field(impedance_rounded_maxwell_space.GetNDSpace(),
                                               true);
  impedance_rounded_maxwell_field.Real().ProjectCoefficient(field_coefficient);
  impedance_rounded_maxwell_field.Imag() = 0.0;
  const auto impedance_rounded_maxwell_result =
      impedance_rounded_maxwell_response.GetMaxwellResponse(impedance_rounded_maxwell_field,
                                                            0.0);
  CHECK(impedance_rounded_maxwell_result.loop_residual < 1.0e-10);
  CHECK_THAT(impedance_rounded_maxwell_result.matched_length_fraction,
             WithinAbs(1.0, 1.0e-12));
  CHECK(impedance_rounded_maxwell_result.corner_neighborhood_fraction == 0.0);
  CHECK_FALSE(impedance_rounded_maxwell_result.boundary_law_verified);
  CHECK_FALSE(impedance_rounded_maxwell_result.closure_independent_confident);

  auto rounded_concave_maxwell_config = rounded_maxwell_island_config;
  rounded_concave_maxwell_config["Solver"]["SurfaceResponseCorrection"]["Library"] =
      rounded_concave_library_3d_path.string();
  rounded_concave_maxwell_config["Boundaries"]["Postprocessing"]["Dielectric"][0]
                                ["EdgeExcludeAttributes"] = {1, 2, 3, 4, 5, 6};
  IoData rounded_concave_maxwell_iodata(rounded_concave_maxwell_config, false);
  rounded_concave_maxwell_iodata.boundaries.cracked_attributes.insert(9);
  std::vector<std::unique_ptr<Mesh>> rounded_concave_maxwell_meshes;
  rounded_concave_maxwell_meshes.push_back(
      std::make_unique<Mesh>(MakeIslandMesh(true, true, true)));
  SpaceOperator rounded_concave_maxwell_space(rounded_concave_maxwell_iodata,
                                              rounded_concave_maxwell_meshes);
  SurfaceResponseOperator rounded_concave_maxwell_response(rounded_concave_maxwell_iodata,
                                                           rounded_concave_maxwell_space);
  CHECK(rounded_concave_maxwell_response.GetPatchCount() ==
        rounded_maxwell_island_response.GetPatchCount());

  GridFunction rounded_concave_maxwell_field(rounded_concave_maxwell_space.GetNDSpace(),
                                             true);
  rounded_concave_maxwell_field.Real().ProjectCoefficient(field_coefficient);
  rounded_concave_maxwell_field.Imag() = 0.0;
  const auto rounded_concave_maxwell_result =
      rounded_concave_maxwell_response.GetMaxwellResponse(rounded_concave_maxwell_field,
                                                          0.0);
  CHECK(rounded_concave_maxwell_result.loop_residual < 1.0e-10);
  CHECK(rounded_concave_maxwell_result.corner_neighborhood_fraction == 0.0);

  auto interpolated_rounded_maxwell_island_config = rounded_maxwell_island_config;
  interpolated_rounded_maxwell_island_config["Solver"]["SurfaceResponseCorrection"]
                                            ["Library"] =
                                                interpolated_rounded_library_3d_path
                                                    .string();
  IoData interpolated_rounded_maxwell_island_iodata(
      interpolated_rounded_maxwell_island_config, false);
  interpolated_rounded_maxwell_island_iodata.boundaries.cracked_attributes.insert(9);
  std::vector<std::unique_ptr<Mesh>> interpolated_rounded_maxwell_island_meshes;
  interpolated_rounded_maxwell_island_meshes.push_back(
      std::make_unique<Mesh>(MakeIslandMesh(true, true)));
  SpaceOperator interpolated_rounded_maxwell_island_space(
      interpolated_rounded_maxwell_island_iodata,
      interpolated_rounded_maxwell_island_meshes);
  SurfaceResponseOperator interpolated_rounded_maxwell_island_response(
      interpolated_rounded_maxwell_island_iodata,
      interpolated_rounded_maxwell_island_space);
  GridFunction interpolated_rounded_island_field(
      interpolated_rounded_maxwell_island_space.GetNDSpace(), true);
  interpolated_rounded_island_field.Real().ProjectCoefficient(field_coefficient);
  interpolated_rounded_island_field.Imag() = 0.0;
  const auto interpolated_rounded_island_result =
      interpolated_rounded_maxwell_island_response.GetMaxwellResponse(
          interpolated_rounded_island_field, 0.0);
  CHECK(interpolated_rounded_island_result.loop_residual < 1.0e-10);
  CHECK_THAT(interpolated_rounded_island_result.maximum_library_distance,
             WithinRel(0.25, 1.0e-12));
  CHECK_THAT(interpolated_rounded_island_result.domain_correction,
             WithinRel(rounded_island_result.domain_correction, 1.0e-12));
  for (const auto &[interface, energy] : rounded_island_result.fabricated_surface_energy)
  {
    CHECK_THAT(interpolated_rounded_island_result.fabricated_surface_energy.at(interface),
               WithinRel(energy, 1.0e-12));
    CHECK_THAT(
        interpolated_rounded_island_result.fabricated_surface_energy_fixed_flux.at(
            interface),
        WithinRel(rounded_island_result.fabricated_surface_energy_fixed_flux.at(interface),
                  1.0e-12));
  }

  // A local interface-signature conflict must not invalidate every segment carrying that
  // signature. The two islands are separated by less than 2R only along their facing
  // sides, so Warn omits those local segments while Error remains strict.
  auto mixed_signature_config = convex_maxwell_island_config;
  mixed_signature_config["Boundaries"]["Ground"]["Attributes"] = {1, 2, 3, 4, 5, 6, 9, 10};
  mixed_signature_config["Boundaries"]["Postprocessing"]["Dielectric"] = {
      {{"Index", 4},
       {"Attributes", {9}},
       {"Type", "SA"},
       {"Thickness", 0.002},
       {"Permittivity", 4.0},
       {"AutomaticEdges", true},
       {"EdgeDistances", {0.2}},
       {"EdgeFrameNormal", {0.0, 1.0, 0.0}}},
      {{"Index", 5},
       {"Attributes", {10}},
       {"Type", "MS"},
       {"Thickness", 0.002},
       {"Permittivity", 11.47},
       {"AutomaticEdges", true},
       {"EdgeDistances", {0.2}},
       {"EdgeFrameNormal", {0.0, 1.0, 0.0}}}};
  mixed_signature_config["Solver"]["SurfaceResponseCorrection"]["TargetInterfaces"] = {4,
                                                                                       5};
  mixed_signature_config["Solver"]["SurfaceResponseCorrection"]["UnmatchedPolicy"] = "Warn";
  IoData mixed_signature_iodata(mixed_signature_config, false);
  mixed_signature_iodata.boundaries.cracked_attributes.insert(9);
  mixed_signature_iodata.boundaries.cracked_attributes.insert(10);
  std::vector<std::unique_ptr<Mesh>> mixed_signature_meshes;
  mixed_signature_meshes.push_back(
      std::make_unique<Mesh>(MakeIslandMesh(false, false, false, true)));
  SpaceOperator mixed_signature_space(mixed_signature_iodata, mixed_signature_meshes);
  SurfaceResponseOperator mixed_signature_response(mixed_signature_iodata,
                                                   mixed_signature_space);
  CHECK(mixed_signature_response.GetPatchCount() > 0);

  GridFunction mixed_signature_field(mixed_signature_space.GetNDSpace(), true);
  mixed_signature_field.Real().ProjectCoefficient(field_coefficient);
  mixed_signature_field.Imag() = 0.0;
  const auto mixed_signature_result =
      mixed_signature_response.GetMaxwellResponse(mixed_signature_field, 0.0);
  CHECK(mixed_signature_result.matched_length_fraction > 0.0);
  CHECK(mixed_signature_result.matched_length_fraction < 1.0);
  CHECK(mixed_signature_result.fabricated_surface_energy.count(4) == 1);
  CHECK(mixed_signature_result.fabricated_surface_energy.count(5) == 1);

  auto local_interaction_config = mixed_signature_config;
  local_interaction_config["Boundaries"]["Postprocessing"]["Dielectric"] = {
      {{"Index", 4},
       {"Attributes", {9, 10}},
       {"Type", "SA"},
       {"Thickness", 0.002},
       {"Permittivity", 4.0},
       {"AutomaticEdges", true},
       {"EdgeDistances", {0.2}},
       {"EdgeFrameNormal", {0.0, 1.0, 0.0}}}};
  local_interaction_config["Solver"]["SurfaceResponseCorrection"]["TargetInterfaces"] = {4};
  IoData local_interaction_iodata(local_interaction_config, false);
  local_interaction_iodata.boundaries.cracked_attributes.insert(9);
  local_interaction_iodata.boundaries.cracked_attributes.insert(10);
  std::vector<std::unique_ptr<Mesh>> local_interaction_meshes;
  local_interaction_meshes.push_back(
      std::make_unique<Mesh>(MakeIslandMesh(false, false, false, true)));
  SpaceOperator local_interaction_space(local_interaction_iodata, local_interaction_meshes);
  SurfaceResponseOperator local_interaction_response(local_interaction_iodata,
                                                     local_interaction_space);
  CHECK(local_interaction_response.GetPatchCount() > 0);

  GridFunction local_interaction_field(local_interaction_space.GetNDSpace(), true);
  local_interaction_field.Real().ProjectCoefficient(field_coefficient);
  local_interaction_field.Imag() = 0.0;
  const auto local_interaction_result =
      local_interaction_response.GetMaxwellResponse(local_interaction_field, 0.0);
  CHECK(local_interaction_result.matched_length_fraction > 0.0);
  CHECK(local_interaction_result.matched_length_fraction < 1.0);
  CHECK(local_interaction_result.fabricated_surface_energy.count(4) == 1);

  auto strict_mixed_signature_config = mixed_signature_config;
  strict_mixed_signature_config["Solver"]["SurfaceResponseCorrection"]["UnmatchedPolicy"] =
      "Error";
  IoData strict_mixed_signature_iodata(strict_mixed_signature_config, false);
  strict_mixed_signature_iodata.boundaries.cracked_attributes.insert(9);
  strict_mixed_signature_iodata.boundaries.cracked_attributes.insert(10);
  std::vector<std::unique_ptr<Mesh>> strict_mixed_signature_meshes;
  strict_mixed_signature_meshes.push_back(
      std::make_unique<Mesh>(MakeIslandMesh(false, false, false, true)));
  SpaceOperator strict_mixed_signature_space(strict_mixed_signature_iodata,
                                             strict_mixed_signature_meshes);
  CHECK_THROWS_WITH(
      SurfaceResponseOperator(strict_mixed_signature_iodata, strict_mixed_signature_space),
      Catch::Matchers::ContainsSubstring("different interface mapping"));
#endif
}

TEST_CASE_METHOD(
    test::SurfaceResponseFiles, "SurfaceResponseOperator spatial cluster Maxwell",
    "[surfaceresponseoperator][3d][maxwell][spatial][placement][Serial][Parallel]")
{
#if !defined(MFEM_USE_GSLIB)
  SKIP("SurfaceResponseOperator requires MFEM_USE_GSLIB");
#else
  json convex_maxwell_island_config = ConvexMaxwellIslandConfig();
  mfem::VectorConstantCoefficient field_coefficient = ConstantFieldCoefficient();

  // Two disconnected island corners separated diagonally by less than 2R form one
  // localized four-edge neighborhood. It contains perpendicular and endpoint-adjacent
  // interactions, so no longitudinal pair or parallel-cluster coupon can represent it.
  auto spatial_cluster_config = convex_maxwell_island_config;
  spatial_cluster_config["Boundaries"]["Ground"]["Attributes"] = {1, 2, 3, 4, 5, 6, 9, 10};
  spatial_cluster_config["Boundaries"]["Postprocessing"]["Dielectric"][0]["Attributes"] = {
      9, 10};
  spatial_cluster_config["Solver"]["SurfaceResponseCorrection"]["Library"] =
      spatial_cluster_library_3d_path.string();
  IoData spatial_cluster_iodata(spatial_cluster_config, false);
  spatial_cluster_iodata.boundaries.cracked_attributes.insert(9);
  spatial_cluster_iodata.boundaries.cracked_attributes.insert(10);
  std::vector<std::unique_ptr<Mesh>> spatial_cluster_meshes;
  spatial_cluster_meshes.push_back(std::make_unique<Mesh>(MakeOffsetCornerPairMesh()));
  SpaceOperator spatial_cluster_space(spatial_cluster_iodata, spatial_cluster_meshes);
  SurfaceResponseOperator spatial_cluster_response(spatial_cluster_iodata,
                                                   spatial_cluster_space);
  CHECK(spatial_cluster_response.GetPatchCount() > 0);
  CHECK(spatial_cluster_response.GetTargetInterfaces() == std::set<int>{4});
  CHECK_THAT(spatial_cluster_response.GetPatchWeight(), WithinRel(17.25, 1.0e-12));

  GridFunction spatial_cluster_field(spatial_cluster_space.GetNDSpace(), true);
  spatial_cluster_field.Real().ProjectCoefficient(field_coefficient);
  spatial_cluster_field.Imag() = 0.0;
  const auto spatial_cluster_result =
      spatial_cluster_response.GetMaxwellResponse(spatial_cluster_field, 0.0);
  CHECK(std::abs(spatial_cluster_result.domain_correction) > 0.0);
  CHECK(spatial_cluster_result.loop_residual < 1.0e-10);
  CHECK_THAT(spatial_cluster_result.matched_length_fraction, WithinAbs(1.0, 1.0e-12));
  CHECK(spatial_cluster_result.corner_neighborhood_fraction == 0.0);

  // Maxwell response correction refuses cap-interior trace coefficients
  // (InteriorTraceCount): every Maxwell coefficient must lie on a contour.
  {
    auto cap_hat_maxwell_config = spatial_cluster_config;
    cap_hat_maxwell_config["Solver"]["SurfaceResponseCorrection"]["Library"] =
        cap_hat_spatial_cluster_library_3d_path.string();
    IoData cap_hat_maxwell_iodata(cap_hat_maxwell_config, false);
    cap_hat_maxwell_iodata.boundaries.cracked_attributes.insert(9);
    cap_hat_maxwell_iodata.boundaries.cracked_attributes.insert(10);
    CHECK_THROWS_WITH(
        SurfaceResponseOperator(cap_hat_maxwell_iodata, spatial_cluster_space),
        Catch::Matchers::ContainsSubstring(
            "Maxwell response correction does not support cap-interior"));
  }

  auto missing_spatial_cluster_config = spatial_cluster_config;
  missing_spatial_cluster_config["Solver"]["SurfaceResponseCorrection"]["Library"] =
      convex_library_3d_path.string();
  missing_spatial_cluster_config["Solver"]["SurfaceResponseCorrection"]["UnmatchedPolicy"] =
      "Warn";
  IoData missing_spatial_cluster_iodata(missing_spatial_cluster_config, false);
  missing_spatial_cluster_iodata.boundaries.cracked_attributes.insert(9);
  missing_spatial_cluster_iodata.boundaries.cracked_attributes.insert(10);
  const auto spatial_requirements_path =
      temp.temp_dir / "surface-response-requirements-spatial.json";
  WriteSurfaceResponseRequirements(missing_spatial_cluster_iodata,
                                   *spatial_cluster_meshes.back(),
                                   spatial_requirements_path.string());
  std::ifstream spatial_requirements_input(spatial_requirements_path);
  REQUIRE(spatial_requirements_input);
  const json spatial_requirements = json::parse(spatial_requirements_input);
  const auto spatial_requirement =
      std::find_if(spatial_requirements["LegacyRequirements"].begin(),
                   spatial_requirements["LegacyRequirements"].end(),
                   [](const auto &requirement)
                   {
                     return requirement["Topology"] == "SpatialEdgeCluster" &&
                            requirement["Status"] == "Missing" &&
                            requirement["Geometry"].contains("Edges");
                   });
  REQUIRE(spatial_requirement != spatial_requirements["LegacyRequirements"].end());
  const auto &spatial_edges = (*spatial_requirement)["Geometry"]["Edges"];
  REQUIRE(spatial_edges.size() >= 2);
  std::set<int> spatial_conductors;
  for (const auto &edge : spatial_edges)
  {
    const int conductor = edge["Conductor"];
    CHECK(conductor > 0);
    spatial_conductors.insert(conductor);
    CHECK(edge["Point"].size() == 3);
    CHECK(edge["GapDirection"].size() == 3);
    CHECK(edge["ProcessNormal"].size() == 3);
    REQUIRE(edge["Interval"].size() == 2);
    CHECK(std::isfinite(edge["Interval"][0].get<double>()));
    CHECK(std::isfinite(edge["Interval"][1].get<double>()));
    CHECK(edge["Interval"][0].get<double>() <= 0.0);
    CHECK(edge["Interval"][1].get<double>() >= 0.0);
    CHECK(edge["Interval"][1].get<double>() > edge["Interval"][0].get<double>());
    CHECK(edge["InterfaceSlot"].get<int>() >= 0);
    CHECK(edge["BoundaryCondition"]["Type"] == "PEC");
  }
  CHECK(*spatial_conductors.begin() == 1);
  REQUIRE((*spatial_requirement)["Geometry"].contains("PlanViewFacets"));
  const auto &spatial_facets = (*spatial_requirement)["Geometry"]["PlanViewFacets"];
  REQUIRE(!spatial_facets.empty());
  std::set<int> facet_conductors;
  for (const auto &facet : spatial_facets)
  {
    facet_conductors.insert(facet["Conductor"].get<int>());
    REQUIRE(facet["Points"].size() >= 3);
    for (const auto &point : facet["Points"])
    {
      REQUIRE(point.size() == 3);
      CHECK(std::isfinite(point[0].get<double>()));
      CHECK(std::isfinite(point[1].get<double>()));
      CHECK(std::isfinite(point[2].get<double>()));
    }
  }
  CHECK(facet_conductors == spatial_conductors);
  REQUIRE((*spatial_requirement)["Geometry"].contains("PlanViewBoundary"));
  const auto &plan_view_boundary = (*spatial_requirement)["Geometry"]["PlanViewBoundary"];
  REQUIRE(plan_view_boundary.is_array());
  REQUIRE(!plan_view_boundary.empty());

  std::ifstream empty_spatial_library_input(convex_library_3d_path);
  REQUIRE(empty_spatial_library_input);
  auto empty_spatial_library = json::parse(empty_spatial_library_input);
  empty_spatial_library["Name"] = "unit-test-empty-spatial-process";
  empty_spatial_library["Models"] = json::array();
  const auto empty_spatial_library_path =
      temp.temp_dir / "surface-process-empty-spatial-3d.json";
  std::ofstream empty_spatial_library_output(empty_spatial_library_path);
  empty_spatial_library_output << empty_spatial_library.dump(2) << "\n";
  empty_spatial_library_output.close();
  auto empty_spatial_config = missing_spatial_cluster_config;
  empty_spatial_config["Solver"]["SurfaceResponseCorrection"]["Library"] =
      empty_spatial_library_path.string();
  IoData empty_spatial_iodata(empty_spatial_config, false);
  empty_spatial_iodata.boundaries.cracked_attributes.insert(9);
  empty_spatial_iodata.boundaries.cracked_attributes.insert(10);
  const auto empty_spatial_requirements_path =
      temp.temp_dir / "surface-response-requirements-empty-spatial.json";
  WriteSurfaceResponseRequirements(empty_spatial_iodata, *spatial_cluster_meshes.back(),
                                   empty_spatial_requirements_path.string());
  std::ifstream empty_spatial_requirements_input(empty_spatial_requirements_path);
  REQUIRE(empty_spatial_requirements_input);
  const json empty_spatial_requirements = json::parse(empty_spatial_requirements_input);
  CHECK_FALSE(empty_spatial_requirements["Complete"]);
  CHECK(empty_spatial_requirements["Summary"]["Counts"]["Exact"] == 0);
  const auto empty_spatial_requirement =
      std::find_if(empty_spatial_requirements["LegacyRequirements"].begin(),
                   empty_spatial_requirements["LegacyRequirements"].end(),
                   [](const auto &requirement)
                   {
                     return requirement["Topology"] == "SpatialEdgeCluster" &&
                            requirement["Status"] == "Missing" &&
                            requirement["Geometry"].contains("PlanViewBoundary");
                   });
  REQUIRE(empty_spatial_requirement !=
          empty_spatial_requirements["LegacyRequirements"].end());
  CHECK((*empty_spatial_requirement)["Geometry"].contains("PlanViewFacets"));

  std::ifstream spatial_cluster_library_input(spatial_cluster_library_3d_path);
  REQUIRE(spatial_cluster_library_input);
  auto exact_mask_library = json::parse(spatial_cluster_library_input);
  exact_mask_library["Name"] = "unit-test-process-spatial-cluster-exact-mask-3d";
  auto &exact_mask_model = exact_mask_library["Models"].back();
  exact_mask_model["Name"] = "offset-corner-pair-exact-mask";
  exact_mask_model["Edges"] = spatial_edges;
  std::map<int, json> spatial_references;
  for (auto &edge : exact_mask_model["Edges"])
  {
    edge["BoundaryCondition"] = "PEC";
    spatial_references.try_emplace(edge["Conductor"].get<int>(), edge["Point"]);
  }
  if (spatial_references.size() == 1)
  {
    exact_mask_model.erase("ConductorReferences");
    exact_mask_model.erase("OpenContourPaths");
    exact_mask_model["ContourGroups"] = {4};
    exact_mask_model["Reference"] = spatial_references.begin()->second;
  }
  else
  {
    exact_mask_model.erase("Reference");
    exact_mask_model["ConductorReferences"] = json::array();
    for (const auto &[conductor, reference] : spatial_references)
    {
      (void)conductor;
      exact_mask_model["ConductorReferences"].push_back(reference);
    }
  }
  exact_mask_model["PlanViewBoundary"] = plan_view_boundary;
  const auto exact_mask_library_path =
      temp.temp_dir / "surface-process-spatial-cluster-exact-mask-3d.json";
  std::ofstream exact_mask_output(exact_mask_library_path);
  exact_mask_output << exact_mask_library.dump(2) << "\n";
  exact_mask_output.close();

  auto exact_mask_config = spatial_cluster_config;
  exact_mask_config["Solver"]["SurfaceResponseCorrection"]["Library"] =
      exact_mask_library_path.string();
  IoData exact_mask_iodata(exact_mask_config, false);
  exact_mask_iodata.boundaries.cracked_attributes.insert(9);
  exact_mask_iodata.boundaries.cracked_attributes.insert(10);
  SurfaceResponseOperator exact_mask_response(exact_mask_iodata, spatial_cluster_space);
  CHECK(exact_mask_response.GetPatchCount() == spatial_cluster_response.GetPatchCount());
  CHECK_THAT(exact_mask_response.GetPatchWeight(),
             WithinRel(spatial_cluster_response.GetPatchWeight(), 1.0e-12));

  auto CheckPlanViewMaskStatus = [&](const fs::path &library_path,
                                     std::string_view expected_status,
                                     std::string_view suffix)
  {
    auto config = spatial_cluster_config;
    config["Solver"]["SurfaceResponseCorrection"]["Library"] = library_path.string();
    config["Solver"]["SurfaceResponseCorrection"]["UnmatchedPolicy"] = "Warn";
    IoData iodata(config, false);
    iodata.boundaries.cracked_attributes.insert(9);
    iodata.boundaries.cracked_attributes.insert(10);
    const auto path = temp.temp_dir / ("surface-response-requirements-mask-" +
                                       std::string(suffix) + ".json");
    WriteSurfaceResponseRequirements(iodata, *spatial_cluster_meshes.back(), path.string());
    std::ifstream input(path);
    REQUIRE(input);
    const json manifest = json::parse(input);
    const auto requirement = std::find_if(
        manifest["LegacyRequirements"].begin(), manifest["LegacyRequirements"].end(),
        [&](const auto &entry)
        {
          return entry["Topology"] == "SpatialEdgeCluster" &&
                 entry["Geometry"].contains("Edges") &&
                 entry["Geometry"]["Edges"].size() == spatial_edges.size();
        });
    REQUIRE(requirement != manifest["LegacyRequirements"].end());
    CHECK((*requirement)["Status"] == expected_status);
  };
  CheckPlanViewMaskStatus(exact_mask_library_path, "Exact", "exact");

  auto mismatched_mask_library = exact_mask_library;
  mismatched_mask_library["Name"] = "unit-test-process-spatial-cluster-mismatched-mask-3d";
  auto &mismatched_coordinate =
      mismatched_mask_library["Models"].back()["PlanViewBoundary"][0]["Segments"][0][0][0];
  mismatched_coordinate = mismatched_coordinate.get<long long int>() + 1;
  const auto mismatched_mask_library_path =
      temp.temp_dir / "surface-process-spatial-cluster-mismatched-mask-3d.json";
  std::ofstream mismatched_mask_output(mismatched_mask_library_path);
  mismatched_mask_output << mismatched_mask_library.dump(2) << "\n";
  mismatched_mask_output.close();
  CheckPlanViewMaskStatus(mismatched_mask_library_path, "Missing", "mismatched");

  // Spatial clusters match each physical edge's complete boundary law. The second
  // island is finite impedance while the first remains PEC; its straight edges,
  // corners, and coupled four-edge neighborhood all have exact object-form models.
  auto mixed_impedance_spatial_cluster_config = spatial_cluster_config;
  mixed_impedance_spatial_cluster_config["Boundaries"]["Ground"]["Attributes"] = {
      1, 2, 3, 4, 5, 6, 9};
  mixed_impedance_spatial_cluster_config["Boundaries"]["Impedance"] = {
      {{"Attributes", {10}}, {"Ls", 2.0e-13}}};
  mixed_impedance_spatial_cluster_config["Solver"]["SurfaceResponseCorrection"]["Library"] =
      spatial_cluster_mixed_impedance_library_3d_path.string();
  IoData mixed_impedance_spatial_cluster_iodata(mixed_impedance_spatial_cluster_config,
                                                false);
  mixed_impedance_spatial_cluster_iodata.boundaries.cracked_attributes.insert(9);
  mixed_impedance_spatial_cluster_iodata.boundaries.cracked_attributes.insert(10);
  std::vector<std::unique_ptr<Mesh>> mixed_impedance_spatial_cluster_meshes;
  mixed_impedance_spatial_cluster_meshes.push_back(
      std::make_unique<Mesh>(MakeOffsetCornerPairMesh()));
  SpaceOperator mixed_impedance_spatial_cluster_space(
      mixed_impedance_spatial_cluster_iodata, mixed_impedance_spatial_cluster_meshes);
  SurfaceResponseOperator mixed_impedance_spatial_cluster_response(
      mixed_impedance_spatial_cluster_iodata, mixed_impedance_spatial_cluster_space);
  GridFunction mixed_impedance_spatial_cluster_field(
      mixed_impedance_spatial_cluster_space.GetNDSpace(), true);
  mixed_impedance_spatial_cluster_field.Real().ProjectCoefficient(field_coefficient);
  mixed_impedance_spatial_cluster_field.Imag() = 0.0;
  const auto mixed_impedance_spatial_cluster_result =
      mixed_impedance_spatial_cluster_response.GetMaxwellResponse(
          mixed_impedance_spatial_cluster_field, 0.0);
  CHECK_THAT(mixed_impedance_spatial_cluster_result.matched_length_fraction,
             WithinAbs(1.0, 1.0e-12));
  CHECK(mixed_impedance_spatial_cluster_result.corner_neighborhood_fraction == 0.0);
  CHECK(mixed_impedance_spatial_cluster_result.boundary_law_verified);

  auto parameter_mismatch_spatial_cluster_config = mixed_impedance_spatial_cluster_config;
  parameter_mismatch_spatial_cluster_config
      ["Solver"]["SurfaceResponseCorrection"]["Library"] =
          spatial_cluster_parameter_mismatch_library_3d_path.string();
  parameter_mismatch_spatial_cluster_config["Solver"]["SurfaceResponseCorrection"]
                                           ["UnmatchedPolicy"] = "Warn";
  IoData parameter_mismatch_spatial_cluster_iodata(
      parameter_mismatch_spatial_cluster_config, false);
  parameter_mismatch_spatial_cluster_iodata.boundaries.cracked_attributes.insert(9);
  parameter_mismatch_spatial_cluster_iodata.boundaries.cracked_attributes.insert(10);
  SurfaceResponseOperator parameter_mismatch_spatial_cluster_response(
      parameter_mismatch_spatial_cluster_iodata, mixed_impedance_spatial_cluster_space);
  CHECK(parameter_mismatch_spatial_cluster_response.GetPatchCount() > 0);
  CHECK(parameter_mismatch_spatial_cluster_response.GetPatchCount() <
        mixed_impedance_spatial_cluster_response.GetPatchCount());
  const auto parameter_mismatch_spatial_cluster_result =
      parameter_mismatch_spatial_cluster_response.GetMaxwellResponse(
          mixed_impedance_spatial_cluster_field, 0.0);
  CHECK(parameter_mismatch_spatial_cluster_result.matched_length_fraction > 0.0);
  CHECK(parameter_mismatch_spatial_cluster_result.matched_length_fraction < 1.0);

  auto cross_layer_spatial_cluster_config = spatial_cluster_config;
  cross_layer_spatial_cluster_config["Boundaries"]["Postprocessing"]["Dielectric"] = {
      {{"Index", 4},
       {"Attributes", {9}},
       {"Type", "SA"},
       {"Thickness", 0.002},
       {"Permittivity", 4.0},
       {"AutomaticEdges", true},
       {"EdgeDistances", {0.2}},
       {"EdgeFrameNormal", {0.0, 1.0, 0.0}}},
      {{"Index", 5},
       {"Attributes", {10}},
       {"Type", "SA"},
       {"Thickness", 0.002},
       {"Permittivity", 4.0},
       {"AutomaticEdges", true},
       {"EdgeDistances", {0.2}},
       {"EdgeFrameNormal", {0.0, 1.0, 0.0}}}};
  cross_layer_spatial_cluster_config["Solver"]["SurfaceResponseCorrection"]["Library"] =
      cross_layer_spatial_cluster_library_3d_path.string();
  cross_layer_spatial_cluster_config["Solver"]["SurfaceResponseCorrection"]
                                    ["TargetInterfaces"] = {4, 5};
  IoData cross_layer_spatial_cluster_iodata(cross_layer_spatial_cluster_config, false);
  cross_layer_spatial_cluster_iodata.boundaries.cracked_attributes.insert(9);
  cross_layer_spatial_cluster_iodata.boundaries.cracked_attributes.insert(10);
  std::vector<std::unique_ptr<Mesh>> cross_layer_spatial_cluster_meshes;
  cross_layer_spatial_cluster_meshes.push_back(
      std::make_unique<Mesh>(MakeOffsetCornerPairMesh()));
  SpaceOperator cross_layer_spatial_cluster_space(cross_layer_spatial_cluster_iodata,
                                                  cross_layer_spatial_cluster_meshes);
  SurfaceResponseOperator cross_layer_spatial_cluster_response(
      cross_layer_spatial_cluster_iodata, cross_layer_spatial_cluster_space);
  CHECK(cross_layer_spatial_cluster_response.GetPatchCount() > 0);
  CHECK(cross_layer_spatial_cluster_response.GetTargetInterfaces() == std::set<int>{4, 5});

  GridFunction cross_layer_spatial_cluster_field(
      cross_layer_spatial_cluster_space.GetNDSpace(), true);
  cross_layer_spatial_cluster_field.Real().ProjectCoefficient(field_coefficient);
  cross_layer_spatial_cluster_field.Imag() = 0.0;
  const auto cross_layer_spatial_cluster_result =
      cross_layer_spatial_cluster_response.GetMaxwellResponse(
          cross_layer_spatial_cluster_field, 0.0);
  CHECK(cross_layer_spatial_cluster_result.loop_residual < 1.0e-10);
  CHECK_THAT(cross_layer_spatial_cluster_result.matched_length_fraction,
             WithinAbs(1.0, 1.0e-12));
  CHECK(cross_layer_spatial_cluster_result.corner_neighborhood_fraction == 0.0);
  CHECK(cross_layer_spatial_cluster_result.fabricated_surface_energy.count(4) == 1);
  CHECK(cross_layer_spatial_cluster_result.fabricated_surface_energy.count(5) == 1);

  const auto cross_layer_requirements_path =
      temp.temp_dir / "surface-response-requirements-cross-layer.json";
  WriteSurfaceResponseRequirements(cross_layer_spatial_cluster_iodata,
                                   *cross_layer_spatial_cluster_meshes.back(),
                                   cross_layer_requirements_path.string());
  std::ifstream cross_layer_requirements_input(cross_layer_requirements_path);
  REQUIRE(cross_layer_requirements_input);
  const json cross_layer_requirements = json::parse(cross_layer_requirements_input);
  const auto cross_layer_requirement =
      std::find_if(cross_layer_requirements["LegacyRequirements"].begin(),
                   cross_layer_requirements["LegacyRequirements"].end(),
                   [](const auto &requirement)
                   {
                     if (requirement["Topology"] != "SpatialEdgeCluster" ||
                         !requirement["Geometry"].contains("Edges"))
                     {
                       return false;
                     }
                     std::set<int> slots;
                     for (const auto &edge : requirement["Geometry"]["Edges"])
                     {
                       slots.insert(edge["InterfaceSlot"].template get<int>());
                     }
                     return slots.size() == 2;
                   });
  REQUIRE(cross_layer_requirement != cross_layer_requirements["LegacyRequirements"].end());
  std::set<int> cross_layer_slots;
  for (const auto &edge : (*cross_layer_requirement)["Geometry"]["Edges"])
  {
    cross_layer_slots.insert(edge["InterfaceSlot"].get<int>());
    CHECK(edge["Conductor"].get<int>() > 0);
    CHECK(edge["BoundaryCondition"]["Type"] == "PEC");
  }
  CHECK(cross_layer_slots.size() == 2);

  // Cross-layer replacement is all-or-nothing. A library with no multi-slot model, or
  // with a multi-slot model that omits one physical target type, must leave the
  // interaction neighborhood unmatched instead of applying independent local coupons.
  auto CheckCrossLayerSpatialClusterMismatch = [&](const fs::path &library)
  {
    CAPTURE(library.string());
    auto mismatch_config = cross_layer_spatial_cluster_config;
    mismatch_config["Solver"]["SurfaceResponseCorrection"]["Library"] = library.string();
    mismatch_config["Solver"]["SurfaceResponseCorrection"]["UnmatchedPolicy"] = "Warn";
    IoData mismatch_iodata(mismatch_config, false);
    mismatch_iodata.boundaries.cracked_attributes.insert(9);
    mismatch_iodata.boundaries.cracked_attributes.insert(10);
    SurfaceResponseOperator mismatch_response(mismatch_iodata,
                                              cross_layer_spatial_cluster_space);
    CHECK(mismatch_response.GetPatchCount() > 0);
    CHECK(mismatch_response.GetPatchCount() <
          cross_layer_spatial_cluster_response.GetPatchCount());
    const auto result =
        mismatch_response.GetMaxwellResponse(cross_layer_spatial_cluster_field, 0.0);
    CHECK(result.matched_length_fraction > 0.0);
    CHECK(result.matched_length_fraction < 1.0);
  };
  CheckCrossLayerSpatialClusterMismatch(spatial_cluster_library_3d_path);
  CheckCrossLayerSpatialClusterMismatch(
      incomplete_cross_layer_spatial_cluster_library_3d_path);

  // A spatial cluster is an exact local topology signature. Perturbing one edge's
  // position or orientation, changing its metal boundary condition, or presenting one
  // more physical edge than the coupon describes must omit the interaction neighborhood.
  auto CheckSpatialClusterMismatch = [&](const fs::path &library)
  {
    CAPTURE(library.string());
    auto mismatch_config = spatial_cluster_config;
    mismatch_config["Solver"]["SurfaceResponseCorrection"]["Library"] = library.string();
    mismatch_config["Solver"]["SurfaceResponseCorrection"]["UnmatchedPolicy"] = "Warn";
    IoData mismatch_iodata(mismatch_config, false);
    mismatch_iodata.boundaries.cracked_attributes.insert(9);
    mismatch_iodata.boundaries.cracked_attributes.insert(10);
    SurfaceResponseOperator mismatch_response(mismatch_iodata, spatial_cluster_space);
    CHECK(mismatch_response.GetPatchCount() > 0);
    CHECK(mismatch_response.GetPatchCount() < spatial_cluster_response.GetPatchCount());
    const auto result = mismatch_response.GetMaxwellResponse(spatial_cluster_field, 0.0);
    CHECK(result.matched_length_fraction > 0.0);
    CHECK(result.matched_length_fraction < 1.0);
  };
  CheckSpatialClusterMismatch(spatial_cluster_position_mismatch_library_3d_path);
  CheckSpatialClusterMismatch(spatial_cluster_orientation_mismatch_library_3d_path);
  CheckSpatialClusterMismatch(spatial_cluster_interval_mismatch_library_3d_path);
  CheckSpatialClusterMismatch(spatial_cluster_extra_edge_library_3d_path);
  CheckSpatialClusterMismatch(spatial_cluster_impedance_mismatch_library_3d_path);

#endif
}

TEST_CASE_METHOD(test::SurfaceResponseFiles, "SurfaceResponseOperator paired aperture",
                 "[surfaceresponseoperator][3d][maxwell][aperture][Serial][Parallel]")
{
#if !defined(MFEM_USE_GSLIB)
  SKIP("SurfaceResponseOperator requires MFEM_USE_GSLIB");
#else
  json convex_maxwell_island_config = ConvexMaxwellIslandConfig();
  mfem::VectorConstantCoefficient field_coefficient = ConstantFieldCoefficient();

  // Two aperture perimeters are disconnected edge-graph components even though the
  // metal bridge between them is one physical PEC strip. Its outward-facing edges must
  // select a strip coupon from local geometry instead of being rejected as different
  // conductors.
  auto MakePairedApertureMesh = []()
  {
    constexpr double extent = 3.0;
    mfem::Mesh serial = mfem::Mesh::MakeCartesian3D(24, 4, 24, mfem::Element::HEXAHEDRON,
                                                    extent, 1.0, extent);
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
      double xmin = extent, xmax = 0.0;
      double zmin = extent, zmax = 0.0;
      for (const int vertex : vertices)
      {
        const double *point = serial.GetVertex(vertex);
        on_plane = on_plane && std::abs(point[1] - 0.5) < 1.0e-12;
        xmin = std::min(xmin, point[0]);
        xmax = std::max(xmax, point[0]);
        zmin = std::min(zmin, point[2]);
        zmax = std::max(zmax, point[2]);
      }
      const bool left_aperture = xmin >= 0.75 - 1.0e-12 && xmax <= 1.25 + 1.0e-12 &&
                                 zmin >= 0.5 - 1.0e-12 && zmax <= 2.5 + 1.0e-12;
      const bool right_aperture = xmin >= 1.5 - 1.0e-12 && xmax <= 2.0 + 1.0e-12 &&
                                  zmin >= 0.5 - 1.0e-12 && zmax <= 2.5 + 1.0e-12;
      if (on_plane && !left_aperture && !right_aperture)
      {
        serial.AddBdrElement(serial.GetFace(face)->Duplicate(&serial));
        serial.SetBdrAttribute(serial.GetNBE() - 1, 9);
      }
    }
    serial.FinalizeTopology();
    serial.Finalize();
    return std::make_unique<mfem::ParMesh>(Mpi::World(), serial);
  };

  auto strip_aperture_config = convex_maxwell_island_config;
  strip_aperture_config["Solver"]["SurfaceResponseCorrection"]["Library"] =
      strip_library_3d_path.string();
  strip_aperture_config["Solver"]["SurfaceResponseCorrection"]["UnmatchedPolicy"] = "Warn";
  strip_aperture_config["Boundaries"]["Postprocessing"]["Dielectric"][0]
                       ["EdgeExcludeAttributes"] = {1, 2, 3, 4, 5, 6};
  IoData strip_aperture_iodata(strip_aperture_config, false);
  strip_aperture_iodata.boundaries.cracked_attributes.insert(9);
  auto strip_geometry_mesh = MakePairedApertureMesh();
  const auto strip_geometry =
      ExtractMetalEdgeGeometry(*strip_geometry_mesh, strip_aperture_iodata.boundaries,
                               JointNoiseExtractionFor(strip_aperture_iodata.boundaries));
  auto strip_segments =
      GetInterfaceMetalEdgeSegmentIndices(strip_geometry, 4, InterfaceDielectric::SA);
  ExcludeMetalEdgeSegmentIndices(*strip_geometry_mesh, strip_geometry, {1, 2, 3, 4, 5, 6},
                                 strip_segments);
  std::set<int> strip_perimeter_components;
  for (const std::size_t segment : strip_segments)
  {
    strip_perimeter_components.insert(strip_geometry.segments[segment].component);
  }
  CHECK(strip_perimeter_components.size() == 2);

  std::vector<std::unique_ptr<Mesh>> strip_aperture_meshes;
  strip_aperture_meshes.push_back(std::make_unique<Mesh>(MakePairedApertureMesh()));
  SpaceOperator strip_aperture_space(strip_aperture_iodata, strip_aperture_meshes);
  SurfaceResponseOperator strip_aperture_response(strip_aperture_iodata,
                                                  strip_aperture_space);
  CHECK(strip_aperture_response.GetPatchCount() > 0);

  auto no_strip_aperture_config = strip_aperture_config;
  no_strip_aperture_config["Solver"]["SurfaceResponseCorrection"]["Library"] =
      concave_library_3d_path.string();
  IoData no_strip_aperture_iodata(no_strip_aperture_config, false);
  no_strip_aperture_iodata.boundaries.cracked_attributes.insert(9);
  std::vector<std::unique_ptr<Mesh>> no_strip_aperture_meshes;
  no_strip_aperture_meshes.push_back(std::make_unique<Mesh>(MakePairedApertureMesh()));
  SpaceOperator no_strip_aperture_space(no_strip_aperture_iodata, no_strip_aperture_meshes);
  SurfaceResponseOperator no_strip_aperture_response(no_strip_aperture_iodata,
                                                     no_strip_aperture_space);

  GridFunction strip_aperture_field(strip_aperture_space.GetNDSpace(), true);
  strip_aperture_field.Real().ProjectCoefficient(field_coefficient);
  strip_aperture_field.Imag() = 0.0;
  const auto strip_aperture_result =
      strip_aperture_response.GetMaxwellResponse(strip_aperture_field, 0.0);
  GridFunction no_strip_aperture_field(no_strip_aperture_space.GetNDSpace(), true);
  no_strip_aperture_field.Real().ProjectCoefficient(field_coefficient);
  no_strip_aperture_field.Imag() = 0.0;
  const auto no_strip_aperture_result =
      no_strip_aperture_response.GetMaxwellResponse(no_strip_aperture_field, 0.0);
  CHECK(strip_aperture_result.matched_length_fraction >
        no_strip_aperture_result.matched_length_fraction);

#endif
}

namespace
{

// Shared setup of the SurfaceResponseOperatorCornerTraceBasis cases (the corner family's
// trace basis rule, corner-family review 2026-09-29): the libraries (the isolated edge plus
// one convex corner coupon at 90 / 120 / 135 / 165 / 180 degrees on the lane-2 layout or
// the rule's layout, and the family library with the rule's nodes: the legacy tie coupons
// 90 / 135 / 150 / 165 / 180 (exact matches only) and the segment [90, 135] keyed 112.5
// with the per-side 90+ / 135- coupons and the 105 / 120 nodes carrying TraceBasis
// ConnectivityAngleDegrees, corner-qualification block 2026-09-29), the device islands (a
// house (90, 90, 120, 120, 120 degrees), a gable (90, 90, 165, 105, 105, 165) and a steep
// house with two 112.5-degree corners and a 135-degree apex, vertices on grid rays) and
// the runtime helpers: the family IoData with the default Collocated lift or the protocol's
// SurfaceMortar lift (MortarOversampling 2) and the fabricated surface energy per corner
// patch of every runtime model.
struct CornerTraceBasisFixture
{
  test::SharedTempDir temp;
  static constexpr double R = 0.2, t = 0.01, oe = 0.005;
  const CornerTraceBasisRule rule;
  const fs::path points_path = temp.temp_dir / "isolated-points.csv";
  const fs::path isolated_domain_path = temp.temp_dir / "isolated-domain.csv";
  const fs::path isolated_surface_path = temp.temp_dir / "isolated-surface.csv";
  const fs::path corner_fabricated_path = temp.temp_dir / "corner-fabricated.csv";
  const fs::path corner_thin_path = temp.temp_dir / "corner-thin.csv";
  const fs::path corner_fabricated_surface_path =
      temp.temp_dir / "corner-fabricated-surface.csv";
  const fs::path corner_thin_surface_path = temp.temp_dir / "corner-thin-surface.csv";
  json base_library;
  std::map<std::string, fs::path> libraries;
  // The family's nodes: angle, segment connectivity (legacy when absent) and model name
  // (the per-side coupons beside a legacy tie coupon carry the -c<connectivity> suffix).
  struct FamilyNode
  {
    double angle;
    std::optional<double> connectivity;
    std::string name;
  };
  const std::vector<FamilyNode> family_nodes = {{90.0, std::nullopt, "convex-corner-90"},
                                                {90.0, 112.5, "convex-corner-90-c112.5"},
                                                {105.0, 112.5, "convex-corner-105"},
                                                {120.0, 112.5, "convex-corner-120"},
                                                {135.0, 112.5, "convex-corner-135-c112.5"},
                                                {135.0, std::nullopt, "convex-corner-135"},
                                                {150.0, std::nullopt, "convex-corner-150"},
                                                {165.0, std::nullopt, "convex-corner-165"},
                                                {180.0, std::nullopt, "convex-corner-180"}};
  const std::vector<std::array<double, 2>> house = {
      {-1.0, -1.0},
      {1.0, -1.0},
      {1.0, 0.2},
      {0.0, 0.2 + std::tan(30.0 * M_PI / 180.0)},
      {-1.0, 0.2}};
  const double s15 = std::sin(15.0 * M_PI / 180.0), c15 = std::cos(15.0 * M_PI / 180.0);
  const double gable_s = (1.25 - 0.2) / (c15 + s15 / 0.8);
  const double gable_r = (1.0 - gable_s * s15) / 0.8;
  const std::vector<std::array<double, 2>> gable = {{-1.0, -1.0},
                                                    {1.0, -1.0},
                                                    {1.0, 0.2},
                                                    {0.8 * gable_r, gable_r},
                                                    {-0.8 * gable_r, gable_r},
                                                    {-1.0, 0.2}};
  json config;
  std::unique_ptr<mfem::ParMesh> house_mesh, gable_mesh, steep_mesh;
  const std::vector<std::array<double, 2>> steep_house = {
      {-1.0, -1.0},
      {1.0, -1.0},
      {1.0, 0.2},
      {0.0, 0.2 + std::tan(22.5 * M_PI / 180.0)},
      {-1.0, 0.2}};

  CornerTraceBasisFixture();
  json CornerModel(double angle, const CornerBasisFiles &files, bool with_rule) const;
  json ConfigFor(const fs::path &library) const;
  json Requirements(const fs::path &library, mfem::ParMesh &mesh,
                    const std::string &tag) const;
  IoData FamilyIoData(bool mortar) const;
  static Vector ProjectedPotential(LaplaceOperator &laplace);
  // Per runtime model: fabricated surface energy per patch (every corner is one patch of
  // weight one).
  static std::map<std::string, double>
  PerPatchEnergies(const SurfaceResponseOperator &response,
                   const SurfaceResponseOperator::ElectrostaticResponse &result);
  std::pair<std::map<std::string, int>, std::map<std::string, double>>
  CornerEnergies(mfem::ParMesh &mesh, const std::string &tag, bool mortar) const;
};

CornerTraceBasisFixture::CornerTraceBasisFixture()
{
  base_library = {
      {"Version", 3},
      {"TraceLiftVersion", 2},
      {"Name", "unit-test-corner-trace-basis"},
      {"MatchingRadius", R},
      {"Fabrication",
       {{"InterfaceLayers", {{"SA", {{"Thickness", 0.002}, {"Permittivity", 4.0}}}}}}},
      {"Models",
       {{{"Name", "isolated"},
         {"Topology", "IsolatedEdge"},
         {"CouponDepth", R},
         {"FabricatedMatrix", isolated_domain_path.string()},
         {"ThinMatrix", isolated_domain_path.string()},
         {"FabricatedSurfaceMatrix", isolated_surface_path.string()},
         {"ThinSurfaceMatrix", isolated_surface_path.string()},
         {"BasisPoints", points_path.string()},
         {"Interfaces", {{{"Type", "SA"}, {"Coupon", 1}}}}}}}};
  if (Mpi::Root(Mpi::World()))
  {
    {
      std::ofstream output(points_path);
      output << "x,y,z\n-0.16,-0.12,0.0\n0.16,-0.12,0.0\n0.16,0.12,0.0\n-0.16,0.12,0.0\n";
      std::ofstream domain(isolated_domain_path);
      domain << "basis_i,basis_j,Q_ij (J)\n";
      std::ofstream surface(isolated_surface_path);
      surface << "interface,edge,basis_i,basis_j,Q_total_ij (J)\n";
      for (int i = 1; i <= 4; i++)
      {
        for (int j = 1; j <= 4; j++)
        {
          const double value = (i == j ? 2.0 : 0.2) * 1.0e-12;
          if (j >= i)
          {
            domain << i << "," << j << "," << value << "\n";
            surface << "1,1," << i << "," << j << "," << value << "\n";
          }
        }
      }
    }
    WriteCornerMatrices(corner_fabricated_path, corner_fabricated_surface_path, 72, 3.0,
                        0.05, R);
    WriteCornerMatrices(corner_thin_path, corner_thin_surface_path, 72, 1.0, 0.01, R);
    for (const bool rule_layout : {false, true})
    {
      for (const double angle : {90.0, 120.0, 135.0, 165.0, 180.0})
      {
        const std::string tag =
            (rule_layout ? "rule-" : "lane2-") + std::to_string(static_cast<int>(angle));
        const auto files =
            WriteCornerBasisFiles(temp.temp_dir, tag, angle, true, R, t, oe, rule_layout);
        auto library = base_library;
        library["Name"] = "unit-test-corner-trace-basis-" + tag;
        library["Models"].push_back(CornerModel(angle, files, rule_layout));
        libraries[tag] = temp.temp_dir / ("library-" + tag + ".json");
        std::ofstream output(libraries[tag]);
        output << library.dump(2) << "\n";
      }
    }
    auto family = base_library;
    family["Name"] = "unit-test-corner-trace-basis-family";
    for (const auto &node : family_nodes)
    {
      const auto files =
          WriteCornerBasisFiles(temp.temp_dir, "family-" + node.name, node.angle, true, R,
                                t, oe, true, node.connectivity);
      json model = CornerModel(node.angle, files, true);
      model["Name"] = node.name;
      family["Models"].push_back(model);
    }
    libraries["family"] = temp.temp_dir / "library-family.json";
    std::ofstream output(libraries["family"]);
    output << family.dump(2) << "\n";
  }
  else
  {
    for (const bool rule_layout : {false, true})
    {
      for (const double angle : {90.0, 120.0, 135.0, 165.0, 180.0})
      {
        const std::string tag =
            (rule_layout ? "rule-" : "lane2-") + std::to_string(static_cast<int>(angle));
        libraries[tag] = temp.temp_dir / ("library-" + tag + ".json");
      }
    }
    libraries["family"] = temp.temp_dir / "library-family.json";
  }
  Mpi::Barrier(Mpi::World());

  config = {
      {"Problem", {{"Type", "Electrostatic"}, {"Output", temp.temp_dir.string()}}},
      {"Model", {{"Mesh", "unused.msh"}}},
      {"Domains", {{"Materials", {{{"Attributes", {1}}}}}}},
      {"Boundaries",
       {{"Ground", {{"Attributes", {1, 2, 3, 4, 5, 6}}}},
        {"Terminal", {{{"Index", 1}, {"Attributes", {9}}}}},
        {"Postprocessing",
         {{"Dielectric",
           {{{"Index", 4},
             {"Attributes", {9}},
             {"Type", "SA"},
             {"Thickness", 0.002},
             {"Permittivity", 4.0},
             {"AutomaticEdges", true},
             {"EdgeDistances", {R}},
             {"EdgeFrameNormal", {0.0, 1.0, 0.0}}}}}}}}},
      {"Solver",
       {{"Order", 1},
        {"Electrostatic",
         {{"ResponseCorrection",
           {{"Library", ""}, {"TargetInterfaces", {4}}, {"UnmatchedPolicy", "Error"}}}}}}}};
  house_mesh = MakePolygonIslandMesh(house, 8.0, 0.1);
  gable_mesh = MakePolygonIslandMesh(gable, 8.0, 0.1);
  steep_mesh = MakePolygonIslandMesh(steep_house, 8.0, 0.1);
}

json CornerTraceBasisFixture::CornerModel(double angle, const CornerBasisFiles &files,
                                          bool with_rule) const
{
  json model = {
      {"Name", "convex-corner-" + std::to_string(static_cast<int>(angle))},
      {"Topology", "ConvexCorner"},
      {"Angle", angle},
      {"AngleDegrees", angle},
      {"Convexity", "Convex"},
      {"AngleTolerance", 1.0e-6},
      {"CornerRadius", 0.0},
      {"CornerRadiusTolerance", 0.0},
      {"FabricatedMatrix", corner_fabricated_path.string()},
      {"ThinMatrix", corner_thin_path.string()},
      {"FabricatedSurfaceMatrix", corner_fabricated_surface_path.string()},
      {"ThinSurfaceMatrix", corner_thin_surface_path.string()},
      {"BasisPoints", files.points.string()},
      {"TraceMesh",
       {{"Vertices", files.vertices.string()}, {"Triangles", files.triangles.string()}}},
      {"ContourGroups", files.contour_groups},
      {"ZeroTraceIndices", files.zero_trace_indices},
      {"Interfaces", {{{"Type", "SA"}, {"Coupon", 1}}}}};
  if (with_rule)
  {
    model["TraceBasis"] = files.trace_basis;
  }
  return model;
}

json CornerTraceBasisFixture::ConfigFor(const fs::path &library) const
{
  auto result = config;
  result["Solver"]["Electrostatic"]["ResponseCorrection"]["Library"] = library.string();
  return result;
}

json CornerTraceBasisFixture::Requirements(const fs::path &library, mfem::ParMesh &mesh,
                                           const std::string &tag) const
{
  IoData iodata(ConfigFor(library), false);
  iodata.boundaries.cracked_attributes.insert(9);
  const auto manifest_path = temp.temp_dir / ("requirements-" + tag + ".json");
  WriteSurfaceResponseRequirements(iodata, mesh, manifest_path.string());
  Mpi::Barrier(Mpi::World());
  std::ifstream input(manifest_path);
  REQUIRE(input);
  return json::parse(input);
}

IoData CornerTraceBasisFixture::FamilyIoData(bool mortar) const
{
  auto family_config = ConfigFor(libraries.at("family"));
  if (mortar)
  {
    auto &correction = family_config["Solver"]["Electrostatic"]["ResponseCorrection"];
    correction["TraceCoupling"] = "SurfaceMortar";
    correction["MortarOversampling"] = 2;
  }
  IoData iodata(family_config, false);
  iodata.boundaries.cracked_attributes.insert(9);
  return iodata;
}

Vector CornerTraceBasisFixture::ProjectedPotential(LaplaceOperator &laplace)
{
  mfem::ParGridFunction potential(&laplace.GetH1Space().Get());
  mfem::FunctionCoefficient potential_coefficient(
      [](const mfem::Vector &x)
      { return (x[1] - 0.5) * (1.0 + 0.1 * (x[0] - 4.0) - 0.05 * (x[2] - 4.0)); });
  potential.ProjectCoefficient(potential_coefficient);
  Vector potential_true;
  potential.GetTrueDofs(potential_true);
  return potential_true;
}

std::map<std::string, double> CornerTraceBasisFixture::PerPatchEnergies(
    const SurfaceResponseOperator &response,
    const SurfaceResponseOperator::ElectrostaticResponse &result)
{
  std::map<std::string, double> per_patch;
  const auto &names = response.GetModelNames();
  for (const auto &contribution : result.model_contributions)
  {
    const auto &name = names.at(contribution.model);
    if (name.find("corner") != std::string::npos)
    {
      REQUIRE(contribution.patch_count > 0.0);
      per_patch[name] =
          contribution.fabricated_surface_energy.at(4) / contribution.patch_count;
      CHECK(std::isfinite(per_patch[name]));
      CHECK(per_patch[name] > 0.0);
    }
  }
  return per_patch;
}

std::pair<std::map<std::string, int>, std::map<std::string, double>>
CornerTraceBasisFixture::CornerEnergies(mfem::ParMesh &mesh, const std::string &tag,
                                        bool mortar) const
{
  const json manifest = Requirements(libraries.at("family"), mesh, tag);
  std::map<std::string, int> matched;  // model name -> corners
  for (const auto &feature : manifest["Identification"]["Features"])
  {
    if (feature["Type"] == "ConvexCorner")
    {
      REQUIRE(feature["Match"]["Status"] == "Matched");
      matched[feature["Match"]["Model"].get<std::string>()]++;
    }
  }
  IoData iodata = FamilyIoData(mortar);
  std::vector<std::unique_ptr<Mesh>> meshes;
  meshes.push_back(std::make_unique<Mesh>(std::make_unique<mfem::ParMesh>(mesh)));
  LaplaceOperator laplace(iodata, meshes);
  SurfaceResponseOperator response(iodata, laplace);
  const auto result = response.GetElectrostaticResponse(ProjectedPotential(laplace));
  return std::make_pair(matched, PerPatchEnergies(response, result));
}

}  // namespace

// The corner family's trace basis (corner-family review 2026-09-29, root cause of the
// non-90-degree MS / MA over-correction, and the supervisor's rule): (1) the rule's layout
// — the 90-degree convex node is the lane-2 fixed layout, every node has the same zero set
// and 72 knots, the second arm's crossing is a PEC knot, the box corners are slave vertices
// with a partition of unity; (2) the fail-closed library gate — a corner coupon whose metal
// arm crosses a box ring between a free knot and a PEC knot (the lane-2 angle-independent
// layout at 120 / 165 degrees) is refused at library load, the lane-2 layout at 90 / 135 /
// 180 degrees (arms on knot rays) and the rule's layout at every angle load; (3) the
// runtime — an island with exact 120-degree (house) and exact 165 / 105-degree corners
// matched to exact nodes of a rule-built family (the legacy tie coupons ahead of the
// per-side 90+ / 135- coupons) and an interpolated angle (112.5 degrees: cubic on the
// segment [90, 135] keyed 112.5, corner-qualification block 2026-09-29) whose runtime
// basis is constructed by the rule at the device angle with the segment's connectivity,
// run with the default
// Collocated lift (the trace sampled at the knots; the trace mesh is not read) and with the
// protocol's SurfaceMortar lift (MortarOversampling 2: the angle-specific trace meshes with
// their slave box-corner vertices enter the mortar mass and the lifts through
// MortarVertex::ForEachBasis) — under both lifts every 105 / 120 / 165 / 112.5 corner patch
// has a fabricated surface energy within a factor two of the 90-degree patches on the same
// synthetic matrices (the physics of the fix — the fabricated MS of a corner patch within
// 50 % of the isolated 2R scale — is the library gate and the device check on the rebuilt
// coupons); (4) the constructed basis of the interpolated corner (points, slave vertices,
// triangles) round-trips through the response-geometry cache: a mortar run reloading the
// cache written by the previous run reproduces every model contribution.
TEST_CASE_METHOD(CornerTraceBasisFixture, "SurfaceResponseOperatorCornerTraceBasis",
                 "[surfaceresponseoperator][corner][tracebasis][Serial][Parallel]")
{
#if !defined(MFEM_USE_GSLIB)
  SKIP("SurfaceResponseOperator requires MFEM_USE_GSLIB");
#else
  // (1) The rule's layout.
  {
    const auto seed = MakeCornerBoxSeed(R, t, oe, true, rule);
    REQUIRE(seed.points.size() == 72);
    REQUIRE(seed.contour_groups == std::vector<int>(9, 8));
    CHECK(seed.zero_trace_indices == std::vector<int>{28, 29, 30, 36, 37, 38});
    const auto ninety = BuildCornerTraceBasis(
        seed.points, seed.contour_groups, seed.zero_trace_indices, M_PI / 2.0, true, rule);
    for (std::size_t k = 0; k < seed.points.size(); k++)
    {
      for (int d = 0; d < 3; d++)
      {
        CHECK_THAT(ninety.knots[k][d], WithinAbs(seed.points[k][d], 1.0e-14));
      }
    }
    CHECK(ninety.vertices.size() == 72);  // no slave: every corner is a knot
    CHECK(ninety.triangles.size() == 140);
    for (const double angle : {75.0, 105.0, 120.0, 135.0, 150.0, 165.0, 180.0})
    {
      const auto basis =
          BuildCornerTraceBasis(seed.points, seed.contour_groups, seed.zero_trace_indices,
                                angle * M_PI / 180.0, true, rule);
      std::vector<int> zero;
      for (int k = 0; k < 72; k++)
      {
        if (basis.zero[k])
        {
          zero.push_back(k);
        }
      }
      CHECK(zero == seed.zero_trace_indices);
      // The second arm's crossing of the z = 0 ring is the PEC knot at slot 6 (1-based 31).
      const auto crossing = SquarePerimeterPoint(
          R, 0.0, ArmCrossingFractions(R, angle * M_PI / 180.0).second);
      CHECK_THAT(basis.knots[30][0], WithinAbs(crossing[0], 1.0e-12));
      CHECK_THAT(basis.knots[30][1], WithinAbs(crossing[1], 1.0e-12));
      if (angle == 165.0)
      {
        CHECK_THAT(crossing[0], WithinAbs(-R, 1.0e-12));
        CHECK_THAT(crossing[1], WithinAbs(R * std::tan(15.0 * M_PI / 180.0), 1.0e-12));
      }
      // Slave vertices: the box corners that are no knot, with a partition of unity.
      int slaves = 0;
      for (const auto &vertex : basis.vertices)
      {
        if (vertex.basis < 0)
        {
          slaves++;
          CHECK(std::abs(std::abs(vertex.point[0]) - R) < 1.0e-12);
          CHECK(std::abs(std::abs(vertex.point[1]) - R) < 1.0e-12);
          CHECK(vertex.weight_a >= 0.0);
          CHECK(vertex.weight_a <= 1.0);
        }
      }
      CHECK(slaves == (angle == 135.0 ? 6 : 8));
      CHECK(basis.triangles.size() == 140 + 2 * slaves);
      // Every free knot is outside the metal footprint and every ring of the fixed
      // layout is untouched.
      for (int r : {0, 1, 2, 5, 6, 7, 8})
      {
        for (int i = 0; i < 8; i++)
        {
          CHECK_THAT(basis.knots[8 * r + i][0],
                     WithinAbs(seed.points[8 * r + i][0], 1.0e-14));
          CHECK_THAT(basis.knots[8 * r + i][1],
                     WithinAbs(seed.points[8 * r + i][1], 1.0e-14));
        }
      }
      CHECK(CheckCornerBasisCrossings(basis.knots, basis.contour_groups, zero,
                                      angle * M_PI / 180.0, true, 1.0e-9 * R)
                .empty());
    }
    // The concave family: crossings and metal interior at the slots 0 / 6 / 7.
    const auto concave_seed = MakeCornerBoxSeed(R, t, oe, false, rule);
    CHECK(concave_seed.zero_trace_indices == std::vector<int>{24, 30, 31, 32, 38, 39});
    const auto concave = BuildCornerTraceBasis(
        concave_seed.points, concave_seed.contour_groups, concave_seed.zero_trace_indices,
        120.0 * M_PI / 180.0, false, rule);
    CHECK(CheckCornerBasisCrossings(concave.knots, concave.contour_groups,
                                    concave_seed.zero_trace_indices, 120.0 * M_PI / 180.0,
                                    false, 1.0e-9 * R)
              .empty());
    // The fixed layout's eight points by fraction (SquarePerimeterPoint) are the
    // generator's square_ring coordinates (review m6: the C++ places every knot by its
    // fraction, the generator reuses square_ring's coordinates at the fixed fractions; the
    // Python test_ring_points_agree_with_square_perimeter_point pins the same table).
    const std::vector<std::array<double, 3>> fixed = {
        {-R, 0.0, t}, {-R, -R, t}, {0.0, -R, t}, {R, -R, t},
        {R, 0.0, t},  {R, R, t},   {0.0, R, t},  {-R, R, t}};
    for (int k = 0; k < 8; k++)
    {
      const auto point = SquarePerimeterPoint(R, t, k / 8.0);
      for (int d = 0; d < 3; d++)
      {
        CHECK_THAT(point[d], WithinAbs(fixed[k][d], 1.0e-14 * R));  // rounding residues
      }
    }
    // Knot coincidence (review m3): a free knot within 1e-6 of a box corner's fraction
    // snaps onto the corner (no slave there, the knot exactly at the corner); one outside
    // the band keeps its position and the corner its slave, at least 8e-6 R away. Convex
    // free 2 (slot 0) sits at the corner (-R, -R) (fraction 1/8) when the second arm
    // crosses the left side at fraction 15/16, i.e. at (-R, R / 2): 153.435 degrees; the
    // arm's crossing moves by 8 R per unit fraction, free 2 by two thirds of that.
    for (const double offset : {7.5e-7, 3.0e-6})
    {
      const double crossing_y = 0.5 * R - 8.0 * R * offset;
      const double angle = std::atan2(crossing_y, -R);
      const auto basis = BuildCornerTraceBasis(seed.points, seed.contour_groups,
                                               seed.zero_trace_indices, angle, true, rule);
      const bool snapped = offset * 2.0 / 3.0 <= 1.0e-6;
      int slaves = 0;
      double nearest_slave = std::numeric_limits<double>::infinity();
      for (const auto &vertex : basis.vertices)
      {
        if (vertex.basis < 0)
        {
          slaves++;
          nearest_slave =
              std::min(nearest_slave, std::hypot(vertex.point[0] - basis.knots[24][0],
                                                 vertex.point[1] - basis.knots[24][1]));
        }
      }
      CHECK(slaves == (snapped ? 6 : 8));  // both metal rings
      for (const int knot : {24, 32})      // free 2 of the z = 0 and z = t rings
      {
        // The box radius is read off the seed's points (rounding residues of 1e-16 R).
        CHECK_THAT(basis.knots[knot][1], WithinAbs(-R, 1.0e-14));
        if (snapped)
        {
          CHECK_THAT(basis.knots[knot][0], WithinAbs(-R, 1.0e-14));
        }
        else
        {
          // Past the corner on the bottom side, by two thirds of the arm's 8 R offset.
          CHECK_THAT(basis.knots[knot][0],
                     WithinAbs(-R + 8.0 * R * offset * 2.0 / 3.0, 1.0e-12));
        }
      }
      if (!snapped)
      {
        CHECK(nearest_slave >= 8.0e-6 * R);
      }
      CHECK(CheckCornerBasisCrossings(basis.knots, basis.contour_groups,
                                      seed.zero_trace_indices, angle, true, 1.0e-9 * R)
                .empty());
    }
    // The footprint test of the gate is tolerant on the arms (review m4): a FREE knot of a
    // concave straight anchor at (R, -4.4e-16) lies on the first arm (the metal boundary),
    // hence on the closed footprint, and is refused; the same ring with that knot at
    // (R, R / 2) passes.
    {
      const double z = 0.0;
      for (const double y : {-4.4e-16, 0.5 * R})
      {
        const std::vector<std::array<double, 3>> ring = {
            {R, 0.0, z}, {R, y, z},    {R, R, z},    {0.0, R, z},
            {-R, R, z},  {-R, 0.0, z}, {0.0, -R, z}, {-R, -R, z}};
        // PEC: both crossings of the arms with the ring ((R, 0) and (-R, 0)) and the knots
        // on the lower half-plane (the concave anchor's metal).
        const std::vector<int> zero = {0, 5, 6, 7};
        const std::string reason =
            CheckCornerBasisCrossings(ring, {8}, zero, M_PI, false, 1.0e-9 * R);
        if (y < 0.0)
        {
          CHECK(reason.find("free trace knot 2") != std::string::npos);
          CHECK(reason.find("metal footprint") != std::string::npos);
        }
        else
        {
          CHECK(reason.empty());
        }
      }
    }
  }

  // (2) The library gate, fail closed: the lane-2 layout at 120 / 165 degrees is refused
  // (the second arm crosses the ring between a free knot and a PEC knot), at 90 / 135 / 180
  // it loads (arms on knot rays); the rule's layout loads at every angle.
  for (const double angle : {120.0, 165.0})
  {
    IoData iodata(
        ConfigFor(libraries.at("lane2-" + std::to_string(static_cast<int>(angle)))), false);
    iodata.boundaries.cracked_attributes.insert(9);
    const auto manifest_path = temp.temp_dir / "requirements-refused.json";
    CHECK_THROWS_WITH(
        WriteSurfaceResponseRequirements(iodata, *house_mesh, manifest_path.string()),
        ContainsSubstring("trace basis gate") &&
            ContainsSubstring("no PEC (ZeroTraceIndices) knot"));
  }
  // A single-coupon library leaves the house's other corners unmatched (Warn): a lane-2
  // coupon forms no family (its free hats would cross the metal at another angle: the
  // reason names the missing TraceBasis rule), a rule-built 120-degree coupon is an exact
  // node of the house's 120-degree corners and refuses the sharper 90-degree ones.
  config["Solver"]["Electrostatic"]["ResponseCorrection"]["UnmatchedPolicy"] = "Warn";
  for (const double angle : {90.0, 135.0, 180.0})
  {
    const json manifest =
        Requirements(libraries.at("lane2-" + std::to_string(static_cast<int>(angle))),
                     *house_mesh, "lane2-" + std::to_string(static_cast<int>(angle)));
    int matched = 0, missing = 0;
    for (const auto &feature : manifest["Identification"]["Features"])
    {
      if (feature["Type"] != "ConvexCorner")
      {
        continue;
      }
      if (feature["Match"]["Status"] == "Matched")
      {
        matched++;
      }
      else
      {
        CHECK(feature["Match"]["Note"].get<std::string>().find("TraceBasis") !=
              std::string::npos);
        missing++;
      }
    }
    CHECK(matched == (angle == 90.0 ? 2 : 0));
    CHECK(missing == (angle == 90.0 ? 3 : 5));
  }
  for (const double angle : {90.0, 120.0, 135.0, 165.0, 180.0})
  {
    const json manifest =
        Requirements(libraries.at("rule-" + std::to_string(static_cast<int>(angle))),
                     *house_mesh, "rule-" + std::to_string(static_cast<int>(angle)));
    int matched = 0;
    for (const auto &feature : manifest["Identification"]["Features"])
    {
      if (feature["Type"] == "ConvexCorner" && feature["Match"]["Status"] == "Matched")
      {
        matched++;
      }
    }
    CHECK(matched == (angle == 90.0 ? 2 : (angle == 120.0 ? 3 : 0)));
  }
  config["Solver"]["Electrostatic"]["ResponseCorrection"]["UnmatchedPolicy"] = "Error";

  // (3) The runtime on the rule-built family with the default Collocated lift (the
  // trace sampled at the knots; the trace mesh is not read): exact nodes on the house
  // (120) and the gable (165 / 105), the interpolated 112.5-degree angle of the steep
  // house (cubic on the segment [90, 135] keyed 112.5, the runtime basis constructed by the
  // rule with the segment's connectivity).
  {
    const bool mortar = false;
    const std::string lift = "-collocated";
    // Exact nodes are the library models themselves (the signature match within the
    // AngleTolerance; the family is consulted only for angles without a coupon).
    {
      const auto [matched, per_patch] = CornerEnergies(*house_mesh, "house" + lift, mortar);
      REQUIRE(matched.size() == 2);
      CHECK(matched.at("convex-corner-90") == 2);
      CHECK(matched.at("convex-corner-120") == 3);
      const double ninety = per_patch.at("convex-corner-90");
      const double one_twenty = per_patch.at("convex-corner-120");
      CHECK(one_twenty > 0.5 * ninety);
      CHECK(one_twenty < 2.0 * ninety);
    }
    {
      const auto [matched, per_patch] = CornerEnergies(*gable_mesh, "gable" + lift, mortar);
      REQUIRE(matched.size() == 3);
      CHECK(matched.at("convex-corner-90") == 2);
      CHECK(matched.at("convex-corner-165") == 2);
      CHECK(matched.at("convex-corner-105") == 2);
      const double ninety = per_patch.at("convex-corner-90");
      for (const auto &name : {"convex-corner-165", "convex-corner-105"})
      {
        CHECK(per_patch.at(name) > 0.5 * ninety);
        CHECK(per_patch.at(name) < 2.0 * ninety);
      }
    }
    {
      const auto [matched, per_patch] =
          CornerEnergies(*steep_mesh, "steep-house" + lift, mortar);
      REQUIRE(matched.size() == 3);
      CHECK(matched.at("convex-corner-90") == 2);
      CHECK(matched.at("convex-corner-135") == 1);
      std::string interpolated;
      for (const auto &[name, count] : matched)
      {
        if (name.find("@corner-angle112.5-cubic") != std::string::npos)
        {
          interpolated = name;
          CHECK(count == 2);
        }
      }
      REQUIRE(!interpolated.empty());
      // The base is the nearest node of the segment's stencil; 112.5 is equidistant from
      // 105 and 120, so either base is a correct pick (the recorded runs take 120: the
      // models' Angle round trip through radians leaves 120 nearer by a rounding residue;
      // the assertion does not depend on that residue). The exact 90 / 135 corners take
      // the legacy tie coupons.
      CHECK((interpolated.rfind("convex-corner-120@", 0) == 0 ||
             interpolated.rfind("convex-corner-105@", 0) == 0));
      const double ninety = per_patch.at("convex-corner-90");
      CHECK(per_patch.at(interpolated) > 0.5 * ninety);
      CHECK(per_patch.at(interpolated) < 2.0 * ninety);
    }
  }
#endif
}

TEST_CASE_METHOD(
    CornerTraceBasisFixture, "SurfaceResponseOperatorCornerTraceBasisMortar",
    "[surfaceresponseoperator][corner][tracebasis][mortar][cache][Serial][Parallel]")
{
#if !defined(MFEM_USE_GSLIB)
  SKIP("SurfaceResponseOperator requires MFEM_USE_GSLIB");
#else
  // (3) The runtime on the rule-built family with the protocol's SurfaceMortar lift
  // (MortarOversampling 2), which reads the trace meshes (slave box-corner vertices).
  {
    const bool mortar = true;
    const std::string lift = "-mortar";
    // Exact nodes are the library models themselves (the signature match within the
    // AngleTolerance; the family is consulted only for angles without a coupon).
    {
      const auto [matched, per_patch] = CornerEnergies(*house_mesh, "house" + lift, mortar);
      REQUIRE(matched.size() == 2);
      CHECK(matched.at("convex-corner-90") == 2);
      CHECK(matched.at("convex-corner-120") == 3);
      const double ninety = per_patch.at("convex-corner-90");
      const double one_twenty = per_patch.at("convex-corner-120");
      CHECK(one_twenty > 0.5 * ninety);
      CHECK(one_twenty < 2.0 * ninety);
    }
    {
      const auto [matched, per_patch] = CornerEnergies(*gable_mesh, "gable" + lift, mortar);
      REQUIRE(matched.size() == 3);
      CHECK(matched.at("convex-corner-90") == 2);
      CHECK(matched.at("convex-corner-165") == 2);
      CHECK(matched.at("convex-corner-105") == 2);
      const double ninety = per_patch.at("convex-corner-90");
      for (const auto &name : {"convex-corner-165", "convex-corner-105"})
      {
        CHECK(per_patch.at(name) > 0.5 * ninety);
        CHECK(per_patch.at(name) < 2.0 * ninety);
      }
    }
    {
      const auto [matched, per_patch] =
          CornerEnergies(*steep_mesh, "steep-house" + lift, mortar);
      REQUIRE(matched.size() == 3);
      CHECK(matched.at("convex-corner-90") == 2);
      CHECK(matched.at("convex-corner-135") == 1);
      std::string interpolated;
      for (const auto &[name, count] : matched)
      {
        if (name.find("@corner-angle112.5-cubic") != std::string::npos)
        {
          interpolated = name;
          CHECK(count == 2);
        }
      }
      REQUIRE(!interpolated.empty());
      // The base is the nearest node of the segment's stencil; 112.5 is equidistant from
      // 105 and 120, so either base is a correct pick (the recorded runs take 120: the
      // models' Angle round trip through radians leaves 120 nearer by a rounding residue;
      // the assertion does not depend on that residue). The exact 90 / 135 corners take
      // the legacy tie coupons.
      CHECK((interpolated.rfind("convex-corner-120@", 0) == 0 ||
             interpolated.rfind("convex-corner-105@", 0) == 0));
      const double ninety = per_patch.at("convex-corner-90");
      CHECK(per_patch.at(interpolated) > 0.5 * ninety);
      CHECK(per_patch.at(interpolated) < 2.0 * ninety);
    }
  }

  // (4) The response-geometry cache carries the constructed basis of the interpolated
  // 112.5-degree corner (points, slave vertices, triangles): a SurfaceMortar run reloading
  // the cache written by the previous run has identical model contributions.
  {
    IoData iodata = FamilyIoData(true);
    std::vector<std::unique_ptr<Mesh>> meshes;
    meshes.push_back(std::make_unique<Mesh>(std::make_unique<mfem::ParMesh>(*steep_mesh)));
    LaplaceOperator laplace(iodata, meshes);
    const auto cache_path = temp.temp_dir / "response-geometry-steep-house.json";
    test::GeometryCacheEnvGuard cache_env(cache_path.string(), true);
    SurfaceResponseOperator written(iodata, laplace);
    Mpi::Barrier(Mpi::World());
    cache_env.DisableWrite();
    SurfaceResponseOperator loaded(iodata, laplace);
    const Vector potential_true = ProjectedPotential(laplace);
    const auto fresh = written.GetElectrostaticResponse(potential_true);
    const auto reloaded = loaded.GetElectrostaticResponse(potential_true);
    CHECK(loaded.GetModelNames() == written.GetModelNames());
    REQUIRE(reloaded.model_contributions.size() == fresh.model_contributions.size());
    REQUIRE(fresh.model_contributions.size() >= 3);  // 90, 135 and the constructed 112.5
    for (std::size_t m = 0; m < fresh.model_contributions.size(); m++)
    {
      const auto &a = fresh.model_contributions[m];
      const auto &b = reloaded.model_contributions[m];
      CHECK(b.model == a.model);
      CHECK(b.patch_count == a.patch_count);
      CHECK(b.patch_weight == a.patch_weight);
      CHECK(b.domain_correction == a.domain_correction);
      CHECK(b.domain_correction_fixed_flux == a.domain_correction_fixed_flux);
      CHECK(b.fabricated_surface_energy == a.fabricated_surface_energy);
      CHECK(b.fabricated_surface_energy_fixed_flux ==
            a.fabricated_surface_energy_fixed_flux);
    }
    const auto per_patch = PerPatchEnergies(loaded, reloaded);
    bool constructed = false;
    for (const auto &[name, energy] : per_patch)
    {
      constructed =
          constructed || name.find("@corner-angle112.5-cubic") != std::string::npos;
    }
    CHECK(constructed);

    // (5) Decision 352 follow-up (1): the per-patch contributions of every model kind on
    // this island (isolated edges, exact and interpolated corners) sum per model to the
    // ModelContribution (energies, patch count, weight) to roundoff, are ordered by patch
    // index without repetition, and leave the model contributions themselves unchanged.
    const auto with_patches = loaded.GetElectrostaticResponse(potential_true, true, true);
    REQUIRE(with_patches.model_contributions.size() == reloaded.model_contributions.size());
    REQUIRE(!with_patches.patch_contributions.empty());
    CHECK(reloaded.patch_contributions.empty());
    for (std::size_t m = 0; m < reloaded.model_contributions.size(); m++)
    {
      const auto &expected = reloaded.model_contributions[m];
      const auto &a = with_patches.model_contributions[m];
      CHECK(a.domain_correction == expected.domain_correction);
      CHECK(a.fabricated_surface_energy == expected.fabricated_surface_energy);
      double count = 0.0, weight = 0.0, domain = 0.0, domain_fixed_flux = 0.0;
      std::map<int, double> surface, surface_fixed_flux;
      for (const auto &patch : with_patches.patch_contributions)
      {
        if (patch.model != expected.model)
        {
          continue;
        }
        count += 1.0;
        weight += patch.weight;
        domain += patch.domain_correction;
        domain_fixed_flux += patch.domain_correction_fixed_flux;
        for (const auto &[interface, energy] : patch.fabricated_surface_energy)
        {
          surface[interface] += energy;
          surface_fixed_flux[interface] +=
              patch.fabricated_surface_energy_fixed_flux.at(interface);
        }
      }
      CHECK(count == expected.patch_count);
      CHECK_THAT(weight, WithinRel(expected.patch_weight, 1.0e-12));
      CHECK_THAT(domain, WithinRel(expected.domain_correction, 1.0e-10));
      CHECK_THAT(domain_fixed_flux,
                 WithinRel(expected.domain_correction_fixed_flux, 1.0e-10));
      REQUIRE(surface.size() == expected.fabricated_surface_energy.size());
      for (const auto &[interface, energy] : expected.fabricated_surface_energy)
      {
        CHECK_THAT(surface.at(interface), WithinRel(energy, 1.0e-10));
        CHECK_THAT(surface_fixed_flux.at(interface),
                   WithinRel(expected.fabricated_surface_energy_fixed_flux.at(interface),
                             1.0e-10));
      }
    }
    for (std::size_t p = 1; p < with_patches.patch_contributions.size(); p++)
    {
      CHECK(with_patches.patch_contributions[p - 1].patch <
            with_patches.patch_contributions[p].patch);
    }
    // The corrected-field form (no fixed flux) carries the mode's domain correction only.
    const auto corrected_form =
        loaded.GetElectrostaticResponse(potential_true, false, true);
    REQUIRE(corrected_form.patch_contributions.size() ==
            with_patches.patch_contributions.size());
    for (const auto &patch : corrected_form.patch_contributions)
    {
      CHECK(patch.domain_correction_fixed_flux == 0.0);
      CHECK(patch.fabricated_surface_energy_fixed_flux.at(4) == 0.0);
    }
  }
#endif
}

// The PSD fallback of an angle-interpolated corner (decision 374 (B), scope decision 376):
// the stencil's blended fabricated and thin DOMAIN matrices on the free knots must be PSD
// beyond roundoff by the operator's own criterion (min eigenvalue >= -1e-9 x max |eig|,
// PositiveSemidefiniteInverseProduct); otherwise the convex LINEAR blend of the two
// bracketing nodes replaces the stencil for every matrix of the feature, recorded per
// feature (Match.InterpolationFallback: stencil, Lagrange weights, min eigenvalues, linear
// nodes and weights) and in the requirement record (CornerFamily), the base node and the
// constructed basis unchanged; every PSD blend is left untouched (recorded under
// Match.BlendEigenvalues only). Synthetic family on the NON-UNIFORM nodes 90 / 105 / 110 /
// 135 (one segment keyed 112.5): the house's 120-degree corners are cubic on all four with
// the Lagrange weights 1/6, -2, 2.7, 2/15 (sum one). (a) Node 105 carries a diagonal bump
// on the free knot 0 in its fabricated and thin matrices (3 -> 5 and 1 -> 2, x 1e-12): the
// cubic blend's entry is 3 - 2 x 2 = -1 (thin 1 - 2 x 1 = -1), not PSD -> the fallback to
// 110 / 135 with weights 0.6 / 0.4 (base 110), the operator constructs (without the
// fallback PositiveSemidefiniteInverseProduct aborts on the fabricated matrix). (b) The
// same bump on the thin matrix only: the fallback is taken (the operator never tested the
// thin matrix). (c) The bump on a PEC knot (zero-trace index 28) of both matrices: the free
// block is untouched, the cubic stays. (d) The uniform family (identical matrices on every
// node, the fixture's `family` library) on the steep house: the 112.5-degree cubic stays
// with its Lagrange weights, no fallback recorded.
TEST_CASE_METHOD(CornerTraceBasisFixture, "SurfaceResponseOperatorCornerBlendPositivity",
                 "[surfaceresponseoperator][corner][tracebasis][psd][Serial][Parallel]")
{
#if !defined(MFEM_USE_GSLIB)
  SKIP("SurfaceResponseOperator requires MFEM_USE_GSLIB");
#else
  constexpr double connectivity = 112.5;
  const std::vector<double> node_angles = {90.0, 105.0, 110.0, 135.0};
  struct Variant
  {
    std::string tag;
    int bump_index;
    bool bump_fabricated, bump_thin;
  };
  const std::vector<Variant> variants = {{"bump-both", 0, true, true},
                                         {"bump-thin", 0, false, true},
                                         {"bump-pec", 28, true, true}};
  std::map<std::string, fs::path> variant_libraries;
  for (const auto &variant : variants)
  {
    variant_libraries[variant.tag] =
        temp.temp_dir / ("library-psd-" + variant.tag + ".json");
  }
  if (Mpi::Root(Mpi::World()))
  {
    for (const auto &variant : variants)
    {
      const auto fabricated_path =
          temp.temp_dir / ("psd-" + variant.tag + "-fabricated.csv");
      const auto thin_path = temp.temp_dir / ("psd-" + variant.tag + "-thin.csv");
      const auto fabricated_surface_path =
          temp.temp_dir / ("psd-" + variant.tag + "-fabricated-surface.csv");
      const auto thin_surface_path =
          temp.temp_dir / ("psd-" + variant.tag + "-thin-surface.csv");
      WriteCornerMatrices(fabricated_path, fabricated_surface_path, 72, 3.0, 0.05, R,
                          variant.bump_fabricated ? variant.bump_index : -1, 2.0);
      WriteCornerMatrices(thin_path, thin_surface_path, 72, 1.0, 0.01, R,
                          variant.bump_thin ? variant.bump_index : -1, 1.0);
      auto library = base_library;
      library["Name"] = "unit-test-corner-blend-positivity-" + variant.tag;
      for (const double angle : node_angles)
      {
        const auto files = WriteCornerBasisFiles(
            temp.temp_dir,
            "psd-" + variant.tag + "-" + std::to_string(static_cast<int>(angle)), angle,
            true, R, t, oe, true, connectivity);
        json model = CornerModel(angle, files, true);
        if (angle == 105.0)
        {
          model["FabricatedMatrix"] = fabricated_path.string();
          model["ThinMatrix"] = thin_path.string();
          model["FabricatedSurfaceMatrix"] = fabricated_surface_path.string();
          model["ThinSurfaceMatrix"] = thin_surface_path.string();
        }
        library["Models"].push_back(model);
      }
      std::ofstream output(variant_libraries.at(variant.tag));
      output << library.dump(2) << "\n";
    }
  }
  Mpi::Barrier(Mpi::World());

  // The cubic Lagrange weights of 120 on the nodes (sum one; -2 on 105).
  std::vector<double> cubic_weights;
  for (const double node : node_angles)
  {
    double weight = 1.0;
    for (const double other : node_angles)
    {
      if (other != node)
      {
        weight *= (120.0 - other) / (node - other);
      }
    }
    cubic_weights.push_back(weight);
  }
  CHECK_THAT(cubic_weights[1], WithinAbs(-2.0, 1.0e-12));
  CHECK_THAT(cubic_weights[2], WithinAbs(2.7, 1.0e-12));

  auto CornerFeatures = [](const json &manifest, double angle)
  {
    std::vector<json> features;
    for (const auto &feature : manifest["Identification"]["Features"])
    {
      if (feature["Type"] == "ConvexCorner" &&
          std::abs(feature["Signature"]["AngleDegrees"].get<double>() - angle) < 1.0e-6)
      {
        features.push_back(feature);
      }
    }
    return features;
  };
  auto CornerRecord = [](const json &manifest, double angle)
  {
    for (const auto &record : manifest["Requirements"])
    {
      if (record["Topology"] == "ConvexCorner" &&
          std::abs(record["Geometry"]["AngleDegrees"].get<double>() - angle) < 1.0e-6)
      {
        return record;
      }
    }
    FAIL("no ConvexCorner requirement record at " << angle << " deg");
    return json();
  };
  const double tolerance = 1.0e-9;

  // (a) + (b): the fallback is taken, recorded and applied; the operator constructs.
  for (const std::string tag : {"bump-both", "bump-thin"})
  {
    const bool fabricated_bumped = tag == "bump-both";
    const json manifest =
        Requirements(variant_libraries.at(tag), *house_mesh, "psd-" + tag);
    const auto corners = CornerFeatures(manifest, 120.0);
    REQUIRE(corners.size() == 3);
    for (const auto &feature : corners)
    {
      REQUIRE(feature["Match"]["Status"] == "Matched");
      REQUIRE(feature["Match"]["Model"] == "convex-corner-110@corner-angle120-linear");
      CHECK(feature["Match"]["Note"].get<std::string>().find("PSD fallback") !=
            std::string::npos);
      REQUIRE(feature["Match"].contains("InterpolationFallback"));
      const auto &fallback = feature["Match"]["InterpolationFallback"];
      CHECK(fallback["StencilRule"] == "cubic");
      REQUIRE(fallback["Stencil"].size() == 4);
      REQUIRE(fallback["StencilWeights"].size() == 4);
      for (std::size_t k = 0; k < node_angles.size(); k++)
      {
        CHECK(fallback["Stencil"][k]["Name"] ==
              "convex-corner-" + std::to_string(static_cast<int>(node_angles[k])));
        CHECK_THAT(fallback["Stencil"][k]["AngleDegrees"].get<double>(),
                   WithinAbs(node_angles[k], 1.0e-9));
        CHECK_THAT(fallback["StencilWeights"][k].get<double>(),
                   WithinAbs(cubic_weights[k], 1.0e-9));
      }
      const auto &cubic = fallback["MinEigenvalueRelative"];
      CHECK(cubic["DomainPositiveSemidefinite"] == false);
      // The bumped matrix reads a negative mode far beyond roundoff (the free block's
      // entry -1 x 1e-12 against a largest eigenvalue of a few 1e-12); the unbumped
      // fabricated matrix of (b) is the nodes' common matrix (the weights sum to one): PSD.
      CHECK(cubic["ThinMatrix"].get<double>() < -1.0e-2);
      if (fabricated_bumped)
      {
        CHECK(cubic["FabricatedMatrix"].get<double>() < -1.0e-2);
      }
      else
      {
        CHECK(cubic["FabricatedMatrix"].get<double>() >= -tolerance);
      }
      CHECK(cubic["FabricatedSurfaceMatrix"].contains("1"));
      CHECK(cubic["ThinSurfaceMatrix"].contains("1"));
      CHECK_THAT(fallback["NegativeToleranceRelative"].get<double>(),
                 WithinAbs(tolerance, 0.0));
      REQUIRE(fallback["LinearNodes"].size() == 2);
      CHECK(fallback["LinearNodes"][0] == "convex-corner-110");
      CHECK(fallback["LinearNodes"][1] == "convex-corner-135");
      REQUIRE(fallback["LinearWeights"].size() == 2);
      CHECK_THAT(fallback["LinearWeights"][0].get<double>(), WithinAbs(0.6, 1.0e-12));
      CHECK_THAT(fallback["LinearWeights"][1].get<double>(), WithinAbs(0.4, 1.0e-12));
      const auto &linear = fallback["LinearMinEigenvalueRelative"];
      CHECK(linear["DomainPositiveSemidefinite"] == true);
      CHECK(linear["FabricatedMatrix"].get<double>() >= -tolerance);
      CHECK(linear["ThinMatrix"].get<double>() >= -tolerance);
      // The applied blend's record is the linear one.
      REQUIRE(feature["Match"].contains("BlendEigenvalues"));
      CHECK(feature["Match"]["BlendEigenvalues"] == linear);
    }
    // The exact 90-degree corners carry no blend record.
    const auto exact = CornerFeatures(manifest, 90.0);
    REQUIRE(exact.size() == 2);
    for (const auto &feature : exact)
    {
      CHECK(feature["Match"]["Model"] == "convex-corner-90");
      CHECK(!feature["Match"].contains("BlendEigenvalues"));
      CHECK(!feature["Match"].contains("InterpolationFallback"));
    }
    // The requirement record: Interpolated, the LINEAR selection with the fallback record.
    const json record = CornerRecord(manifest, 120.0);
    CHECK(record["Status"] == "Interpolated");
    const auto &family = record["CornerFamily"];
    CHECK(family["InterpolationRule"] == "linear");
    CHECK(family["Base"] == "convex-corner-110");
    REQUIRE(family["Nodes"].size() == 2);
    CHECK(family["Nodes"][0]["Name"] == "convex-corner-110");
    CHECK_THAT(family["Nodes"][0]["Weight"].get<double>(), WithinAbs(0.6, 1.0e-12));
    CHECK(family["Nodes"][1]["Name"] == "convex-corner-135");
    CHECK_THAT(family["Nodes"][1]["Weight"].get<double>(), WithinAbs(0.4, 1.0e-12));
    CHECK_THAT(family["ConnectivityAngleDegrees"].get<double>(),
               WithinAbs(connectivity, 1.0e-9));
    REQUIRE(family.contains("InterpolationFallback"));
    CHECK(family["InterpolationFallback"] ==
          corners.front()["Match"]["InterpolationFallback"]);
    REQUIRE(family.contains("BlendEigenvalues"));
    CHECK(family["BlendEigenvalues"]["DomainPositiveSemidefinite"] == true);

    // The operator on the fallback library: the linear runtime model (0.6 x the 110 node +
    // 0.4 x the 135 node = the nodes' common matrices) constructs and its three patches
    // read a fabricated surface energy within a factor two of the exact 90-degree patches.
    {
      IoData iodata(ConfigFor(variant_libraries.at(tag)), false);
      iodata.boundaries.cracked_attributes.insert(9);
      std::vector<std::unique_ptr<Mesh>> meshes;
      meshes.push_back(
          std::make_unique<Mesh>(std::make_unique<mfem::ParMesh>(*house_mesh)));
      LaplaceOperator laplace(iodata, meshes);
      SurfaceResponseOperator response(iodata, laplace);
      const auto result = response.GetElectrostaticResponse(ProjectedPotential(laplace));
      const auto per_patch = PerPatchEnergies(response, result);
      REQUIRE(per_patch.count("convex-corner-110@corner-angle120-linear") == 1);
      REQUIRE(per_patch.count("convex-corner-90") == 1);
      const double ninety = per_patch.at("convex-corner-90");
      const double interpolated = per_patch.at("convex-corner-110@corner-angle120-linear");
      CHECK(interpolated > 0.5 * ninety);
      CHECK(interpolated < 2.0 * ninety);
      for (const auto &contribution : result.model_contributions)
      {
        if (response.GetModelNames().at(contribution.model) ==
            "convex-corner-110@corner-angle120-linear")
        {
          CHECK_THAT(contribution.patch_count, WithinAbs(3.0, 1.0e-12));
        }
      }
    }
  }

  // (c) The bump on a PEC knot leaves the free block untouched: the cubic stays, recorded
  // PSD, no fallback.
  {
    const json manifest =
        Requirements(variant_libraries.at("bump-pec"), *house_mesh, "psd-bump-pec");
    const auto corners = CornerFeatures(manifest, 120.0);
    REQUIRE(corners.size() == 3);
    for (const auto &feature : corners)
    {
      REQUIRE(feature["Match"]["Status"] == "Matched");
      CHECK(feature["Match"]["Model"] == "convex-corner-110@corner-angle120-cubic");
      CHECK(!feature["Match"].contains("InterpolationFallback"));
      REQUIRE(feature["Match"].contains("BlendEigenvalues"));
      const auto &eigenvalues = feature["Match"]["BlendEigenvalues"];
      CHECK(eigenvalues["DomainPositiveSemidefinite"] == true);
      CHECK(eigenvalues["FabricatedMatrix"].get<double>() >= -tolerance);
      CHECK(eigenvalues["ThinMatrix"].get<double>() >= -tolerance);
    }
    const json record = CornerRecord(manifest, 120.0);
    const auto &family = record["CornerFamily"];
    CHECK(family["InterpolationRule"] == "cubic");
    CHECK(!family.contains("InterpolationFallback"));
    REQUIRE(family["Nodes"].size() == 4);
    for (std::size_t k = 0; k < node_angles.size(); k++)
    {
      CHECK_THAT(family["Nodes"][k]["Weight"].get<double>(),
                 WithinAbs(cubic_weights[k], 1.0e-12));
    }
  }

  // (d) The uniform family (the fixture's: identical matrices on every node) on the steep
  // house: the 112.5-degree cubic on 90 / 105 / 120 / 135 is unchanged (its Lagrange
  // weights to 1e-12), recorded PSD, no fallback.
  {
    const json manifest = Requirements(libraries.at("family"), *steep_mesh, "psd-uniform");
    const auto corners = CornerFeatures(manifest, 112.5);
    REQUIRE(corners.size() == 2);
    for (const auto &feature : corners)
    {
      REQUIRE(feature["Match"]["Status"] == "Matched");
      const std::string model = feature["Match"]["Model"].get<std::string>();
      CHECK(model.find("@corner-angle112.5-cubic") != std::string::npos);
      CHECK(!feature["Match"].contains("InterpolationFallback"));
      REQUIRE(feature["Match"].contains("BlendEigenvalues"));
      CHECK(feature["Match"]["BlendEigenvalues"]["DomainPositiveSemidefinite"] == true);
      CHECK(feature["Match"]["Note"].get<std::string>().find("PSD fallback") ==
            std::string::npos);
    }
    const json record = CornerRecord(manifest, 112.5);
    const auto &family = record["CornerFamily"];
    CHECK(family["InterpolationRule"] == "cubic");
    CHECK(!family.contains("InterpolationFallback"));
    REQUIRE(family["Nodes"].size() == 4);
    const std::vector<double> uniform_nodes = {90.0, 105.0, 120.0, 135.0};
    for (const auto &node : family["Nodes"])
    {
      const double node_angle = node["AngleDegrees"].get<double>();
      double expected = 1.0;
      for (const double other : uniform_nodes)
      {
        if (std::abs(other - node_angle) > 1.0e-6)
        {
          expected *= (112.5 - other) / (node_angle - other);
        }
      }
      CHECK_THAT(node["Weight"].get<double>(), WithinAbs(expected, 1.0e-12));
    }
  }
#endif
}

TEST_CASE("SurfaceResponseOperatorTranslationalStretchOwnership",
          "[surfaceresponseoperator][Serial]")
{
  // Decisions 224 / 236 / 242 / 252: a translational STRETCH (every longitudinal cell of
  // one feature stretch) strictly inside one spatial support's box is RECORDED with its
  // class — a Continuation (it continues a claim of the cluster through the claim cut by
  // the side's OWN edge: a cell on the claim's mesh segment, or parallel, a cell end on the
  // side's own edge (the cell end shifted by the provenance edge offset along AxisU)
  // abutting a claim end along the chain and transversely within the tolerance, AND
  // extending beyond that claim end; the double count of the coupon's straight
  // continuation) or Foreign (a model mismatch, not a double count) — never an
  // abort; a stack end adjacent to a cluster, whose first cells lie inside the box while
  // the stretch continues outside, is the ordinary stack-end configuration and records
  // nothing. Cells of a 3-cell portion along +x from the origin x = 0: [0, 1], [1, 2.5],
  // [2.5, 4] (mesh units); a pair's cells sit on its midline y with its own edge at
  // y + edge_offset (AxisU = +y).
  using Patch = config::ElectrostaticSolverData::ResponseCorrectionPatchData;
  auto Cell = [](int feature, int stretch, double begin, double end, double y = 0.0,
                 double edge_offset = 0.0)
  {
    Patch patch;
    patch.origin = {0.5 * (begin + end), y, 0.0};
    patch.axis_u = {0.0, 1.0, 0.0};
    patch.axis_v = {0.0, 0.0, 1.0};
    patch.axis_w = {1.0, 0.0, 0.0};
    patch.longitudinal_cell = {begin - patch.origin[0], end - patch.origin[0]};
    patch.provenance.feature = feature;
    patch.provenance.stretch = stretch;
    patch.provenance.segment = 7;
    patch.provenance.edge_offset = edge_offset;
    return patch;
  };
  // A spatial support over x in [-3, 1.5] (its claims end at x = 0, the box reaches R = 1.5
  // beyond them), |y| <= 3, |z| <= 2; the cluster claims the edge y = 0 from x = -2.5 to 0
  // (on mesh segment 3) and the edge x = -2.5 from y = 0 to 1 (segment 5); the cells above
  // lie on mesh segment 7.
  const double R = 1.5;
  const double tolerance = kSignatureParameterToleranceOverRadius * R;
  using Claim = Patch::Provenance::Claim;
  SpatialSupportBounds box;
  box.patch = 0;
  box.min = {-3.0, -3.0, -2.0};
  box.max = {1.5, 3.0, 2.0};
  box.claims = {Claim{3, {-2.5, 0.0, 0.0}, {0.0, 0.0, 0.0}},
                Claim{5, {-2.5, 0.0, 0.0}, {-2.5, 1.0, 0.0}}};
  SECTION("a stack end adjacent to the cluster records nothing")
  {
    const std::vector<Patch> patches = {Cell(4, 0, 0.0, 1.0), Cell(4, 0, 1.0, 2.5),
                                        Cell(4, 0, 2.5, 4.0)};
    CHECK(
        FindTranslationalStretchInsideSpatialSupport(patches, {box}, 3, tolerance).empty());
    // Two sides of a pair on one stretch index are judged by their own cells: the far side
    // outside the box keeps its stretch outside too.
    std::vector<Patch> sides = patches;
    sides.push_back(Cell(4, 0, 0.0, 1.0, 2.0));
    CHECK(FindTranslationalStretchInsideSpatialSupport(sides, {box}, 3, tolerance).empty());
  }
  SECTION("a stretch wholly inside the box on a claim's segment is its continuation")
  {
    // The claim boundary cut mesh segment 7 (a claim on it from x = -2.5 to 0): the stretch
    // beyond the cut is recorded as a Continuation whatever its lateral offset (the cells
    // of a pair sit on the pair's midline, here y = 0.5).
    std::vector<Patch> patches = {Cell(4, 0, 0.0, 1.0), Cell(4, 0, 1.0, 2.5),
                                  Cell(4, 0, 2.5, 4.0), Cell(4, 1, 0.3, 0.8, 0.5),
                                  Cell(4, 1, 0.8, 1.2, 0.5)};
    SpatialSupportBounds cut = box;
    cut.claims.push_back(Claim{7, {-2.5, 0.0, 0.0}, {0.0, 0.0, 0.0}});
    const auto records =
        FindTranslationalStretchInsideSpatialSupport(patches, {cut}, 3, tolerance);
    REQUIRE(records.size() == 1);
    const auto &record = records.front();
    CHECK(record.feature == 4);
    CHECK(record.stretch == 1);
    CHECK(record.first_patch == 3);
    CHECK(record.patch_count == 2);
    CHECK(record.spatial_patch == 0);
    CHECK_THAT(record.length, WithinAbs(0.9, 1.0e-12));
    CHECK_THAT(record.lo[0], WithinAbs(0.3, 1.0e-12));
    CHECK_THAT(record.hi[0], WithinAbs(1.2, 1.0e-12));
    CHECK(record.continuation);
    // Without the cut segment among the claims the same stretch (0.3 past the claim end,
    // beyond the tolerance) is foreign.
    const auto foreign =
        FindTranslationalStretchInsideSpatialSupport(patches, {box}, 3, tolerance);
    REQUIRE(foreign.size() == 1);
    CHECK(!foreign.front().continuation);
  }
  SECTION("a stretch abutting a claim end along the chain is its continuation")
  {
    // The claim boundary snapped onto the mesh vertex at x = 0 (segment 3 ends there, the
    // stretch starts on segment 7): parallel, the side's own edge (a pair midline at
    // y = 0.5 whose own edge is y = 0: edge offset -0.5) abutting within 1e-3 R along and
    // across, every cell beyond the claim end -> Continuation; the pair's OTHER side (the
    // same midline cells, own edge y = 1, segment 8) is foreign (decision 252: never on
    // midline proximity alone), as is a perpendicular stretch starting there, a parallel
    // one 0.3 away along the chain, or one abutting the claim x = -2.5 along y but 3.4 away
    // from it transversely. A parallel stretch whose own edge is the claim x = -2.5 (cells
    // at x = -1.8, edge offset -0.7) lying BESIDE it over the claim's own range y in [0, 1]
    // (one end aligned with each claim end) is foreign: it does not extend through the
    // claim cut; the same stretch beyond the claim end (y in [1, 2]) is its continuation.
    const std::vector<Patch> abutting = {Cell(4, 1, 1.0e-4 * R, 0.5, 0.5, -0.5),
                                         Cell(4, 1, 0.5, 1.2, 0.5, -0.5)};
    const auto records =
        FindTranslationalStretchInsideSpatialSupport(abutting, {box}, 3, tolerance);
    REQUIRE(records.size() == 1);
    CHECK(records.front().continuation);
    CHECK_THAT(records.front().length, WithinAbs(1.2 - 1.0e-4 * R, 1.0e-12));
    std::vector<Patch> far_side = {Cell(4, 2, 1.0e-4 * R, 0.5, 0.5, 0.5),
                                   Cell(4, 2, 0.5, 1.2, 0.5, 0.5)};
    for (auto &patch : far_side)
    {
      patch.provenance.segment = 8;
    }
    const auto unclaimed =
        FindTranslationalStretchInsideSpatialSupport(far_side, {box}, 3, tolerance);
    REQUIRE(unclaimed.size() == 1);
    CHECK(!unclaimed.front().continuation);
    // The pre-252 criterion (the midline cells within R of the claim end, no edge identity)
    // would have called the far side a continuation too: with no edge offset the midline
    // itself is the own edge, 0.5 off the claim's line, foreign.
    const auto midline = FindTranslationalStretchInsideSpatialSupport(
        {Cell(4, 1, 1.0e-4 * R, 0.5, 0.5), Cell(4, 1, 0.5, 1.2, 0.5)}, {box}, 3, tolerance);
    REQUIRE(midline.size() == 1);
    CHECK(!midline.front().continuation);
    Patch perpendicular = Cell(4, 1, 0.0, 1.0);
    perpendicular.origin = {0.9, 0.5, 0.0};
    perpendicular.axis_w = {0.0, 1.0, 0.0};
    perpendicular.axis_u = {1.0, 0.0, 0.0};
    const auto turned =
        FindTranslationalStretchInsideSpatialSupport({perpendicular}, {box}, 3, tolerance);
    REQUIRE(turned.size() == 1);
    CHECK(!turned.front().continuation);
    perpendicular.origin = {-1.8, 0.5,
                            0.0};  // y in [0, 1] at x = -1.8: alongside the claim
    perpendicular.provenance.edge_offset = -0.7;  // the own edge is the claim x = -2.5
    const auto beside =
        FindTranslationalStretchInsideSpatialSupport({perpendicular}, {box}, 3, tolerance);
    REQUIRE(beside.size() == 1);
    CHECK(!beside.front().continuation);
    Patch beyond_cell = perpendicular;  // y in [1, 2] at x = -1.8: past the claim end y = 1
    beyond_cell.origin = {-1.8, 1.5, 0.0};
    const auto beyond =
        FindTranslationalStretchInsideSpatialSupport({beyond_cell}, {box}, 3, tolerance);
    REQUIRE(beyond.size() == 1);
    CHECK(beyond.front().continuation);
    CHECK_THAT(beyond.front().length, WithinAbs(1.0, 1.0e-12));
    // The same cells as the pair's other side (own edge x = -1.1, unclaimed) are foreign.
    Patch beyond_far_side = beyond_cell;
    beyond_far_side.provenance.edge_offset = 0.7;
    beyond_far_side.provenance.segment = 8;
    CHECK(!FindTranslationalStretchInsideSpatialSupport({beyond_far_side}, {box}, 3,
                                                        tolerance)
               .front()
               .continuation);
    // Starting within the tolerance before the claim end still counts as beyond it; a
    // stretch straddling the claim end (cells y in [0.7, 1] and [1, 1.7], one cell end
    // exactly on the claim end) reaches back alongside the claim and does not.
    beyond_cell.origin = {-1.8, 1.5 - 0.5e-3 * R, 0.0};
    CHECK(FindTranslationalStretchInsideSpatialSupport({beyond_cell}, {box}, 3, tolerance)
              .front()
              .continuation);
    Patch straddle_before = perpendicular, straddle_after = perpendicular;
    straddle_before.origin = {-1.8, 0.85, 0.0};
    straddle_before.longitudinal_cell = {-0.15, 0.15};
    straddle_after.origin = {-1.8, 1.35, 0.0};
    straddle_after.longitudinal_cell = {-0.35, 0.35};
    const auto straddle = FindTranslationalStretchInsideSpatialSupport(
        {straddle_before, straddle_after}, {box}, 3, tolerance);
    REQUIRE(straddle.size() == 1);
    CHECK_THAT(straddle.front().length, WithinAbs(1.0, 1.0e-12));
    CHECK(!straddle.front().continuation);
    // The direction test also applies to the x = 0 claim end: cells over x in [-1, 0] at
    // y = 0.5 (own edge y = 0) abut it but lie alongside the claim.
    const auto back = FindTranslationalStretchInsideSpatialSupport(
        {Cell(4, 1, -1.0, -0.4, 0.5, -0.5), Cell(4, 1, -0.4, 0.0, 0.5, -0.5)}, {box}, 3,
        tolerance);
    REQUIRE(back.size() == 1);
    CHECK(!back.front().continuation);
    const auto apart = FindTranslationalStretchInsideSpatialSupport(
        {Cell(4, 1, 0.3, 0.8, 0.5, -0.5), Cell(4, 1, 0.8, 1.2, 0.5, -0.5)}, {box}, 3,
        tolerance);
    REQUIRE(apart.size() == 1);
    CHECK(!apart.front().continuation);
  }
  SECTION("a stretch wholly inside the box off every claim is foreign")
  {
    const std::vector<Patch> patches = {Cell(4, 1, -2.0, -1.5, 0.5),
                                        Cell(4, 1, -1.5, -0.8, 0.5)};
    const auto records =
        FindTranslationalStretchInsideSpatialSupport(patches, {box}, 3, tolerance);
    REQUIRE(records.size() == 1);
    CHECK(records.front().feature == 4);
    CHECK(records.front().stretch == 1);
    CHECK_THAT(records.front().length, WithinAbs(1.2, 1.0e-12));
    CHECK(!records.front().continuation);
  }
  SECTION("every stretch inside is recorded, in (feature, stretch) order")
  {
    const std::vector<Patch> patches = {Cell(9, 0, -2.0, -1.0, 0.5), Cell(4, 2, 0.0, 1.2),
                                        Cell(4, 1, -1.0, -0.5, 1.5)};
    const auto records =
        FindTranslationalStretchInsideSpatialSupport(patches, {box}, 3, tolerance);
    REQUIRE(records.size() == 3);
    CHECK(records[0].feature == 4);
    CHECK(records[0].stretch == 1);
    CHECK(!records[0].continuation);
    CHECK(records[1].feature == 4);
    CHECK(records[1].stretch == 2);
    CHECK(records[1].continuation);
    CHECK(records[2].feature == 9);
    CHECK(records[2].stretch == 0);
    CHECK(!records[2].continuation);
  }
  SECTION("touching the box face is not inside")
  {
    const std::vector<Patch> patches = {Cell(4, 1, -2.0, -1.5), Cell(4, 1, -1.5, 1.5)};
    CHECK(
        FindTranslationalStretchInsideSpatialSupport(patches, {box}, 3, tolerance).empty());
  }
  SECTION("point-in-z and unattributed patches are never judged")
  {
    Patch spatial = Cell(5, 0, -1.0, -1.0);  // {0, 0} cell
    Patch explicit_patch = Cell(-1, -1, -2.0, -1.0);
    CHECK(FindTranslationalStretchInsideSpatialSupport({spatial, explicit_patch}, {box}, 3,
                                                       tolerance)
              .empty());
  }
}

TEST_CASE("SurfaceResponseOperatorContinuationOwnership",
          "[surfaceresponseoperator][Serial]")
{
  // Decisions 236 (2) / 244: the cells of a translational stretch that continues a
  // cluster's claim are owned by that coupon inside its box and clipped exactly at the box
  // face; foreign cells, cells outside the box and the stack-end cells of a stretch that
  // continues no claim are untouched; a cell on the continuations of two coupons is removed
  // once and attributed by the midpoint between the two continued claim ends. A pair's or
  // stack's side is owned only where its OWN edge continues the claim (decision 252). Cells
  // along +x on mesh segment 7 at y = 0 carry weight 0.1 x length (the weight is linear in
  // the cell length) and a Maxwell anchor at the origin; a pair's cells sit on its midline
  // y with the side's own edge at y + edge_offset (AxisU = +y).
  using Patch = config::ElectrostaticSolverData::ResponseCorrectionPatchData;
  auto Cell = [](int feature, int stretch, double begin, double end, double y = 0.0,
                 double edge_offset = 0.0)
  {
    Patch patch;
    patch.origin = {0.5 * (begin + end), y, 0.0};
    patch.axis_u = {0.0, 1.0, 0.0};
    patch.axis_v = {0.0, 0.0, 1.0};
    patch.axis_w = {1.0, 0.0, 0.0};
    patch.longitudinal_cell = {begin - patch.origin[0], end - patch.origin[0]};
    patch.weight = 0.1 * (end - begin);
    patch.provenance.feature = feature;
    patch.provenance.stretch = stretch;
    patch.provenance.segment = 7;
    patch.provenance.edge_offset = edge_offset;
    patch.provenance.s0 = 0.0;
    patch.provenance.s1 = 4.0;
    patch.provenance.quadrature_weight = (end - begin) / 4.0;
    patch.maxwell_conductor_anchors = {patch.origin};
    return patch;
  };
  auto Same = [](const Patch &a, const Patch &b)
  {
    for (int d = 0; d < 3; d++)
    {
      CHECK_THAT(a.origin[d], WithinAbs(b.origin[d], 1.0e-12));
      CHECK_THAT(a.maxwell_conductor_anchors.front()[d],
                 WithinAbs(b.maxwell_conductor_anchors.front()[d], 1.0e-12));
    }
    CHECK_THAT(a.longitudinal_cell[0], WithinAbs(b.longitudinal_cell[0], 1.0e-12));
    CHECK_THAT(a.longitudinal_cell[1], WithinAbs(b.longitudinal_cell[1], 1.0e-12));
    CHECK_THAT(a.weight, WithinAbs(b.weight, 1.0e-12));
    CHECK_THAT(a.provenance.quadrature_weight,
               WithinAbs(b.provenance.quadrature_weight, 1.0e-12));
  };
  // The spatial support of the stretch-record case: box x in [-3, 1.5], |y| <= 3, |z| <= 2;
  // the cluster claims the edge y = 0 from x = -2.5 to 0 (segment 3) and the edge x = -2.5
  // from y = 0 to 1 (segment 5); its claim cut at x = 0 continues along +x to the face
  // x = 1.5 (R = 1.5).
  const double R = 1.5;
  const double tolerance = kSignatureParameterToleranceOverRadius * R;
  using Claim = Patch::Provenance::Claim;
  SpatialSupportBounds box;
  box.patch = 0;
  box.min = {-3.0, -3.0, -2.0};
  box.max = {1.5, 3.0, 2.0};
  box.claims = {Claim{3, {-2.5, 0.0, 0.0}, {0.0, 0.0, 0.0}},
                Claim{5, {-2.5, 0.0, 0.0}, {-2.5, 1.0, 0.0}}};
  SECTION("a continuation cell wholly inside the box keeps weight 0")
  {
    std::vector<Patch> patches = {Cell(4, 0, 0.0, 1.0), Cell(4, 0, 1.0, 1.4)};
    const auto ownership = ApplyContinuationOwnership(patches, {box}, 3, tolerance, R);
    REQUIRE(ownership.cells.size() == 2);
    CHECK(ownership.wholly_owned_cells == 2);
    CHECK(ownership.clipped_cells == 0);
    CHECK(ownership.shared_cells == 0);
    CHECK_THAT(ownership.owned_length, WithinAbs(1.4, 1.0e-12));
    CHECK_THAT(ownership.owned_by_support.at(0), WithinAbs(1.4, 1.0e-12));
    CHECK_THAT(ownership.owned_by_stretch.at(std::make_tuple(4, 0, std::size_t{0})),
               WithinAbs(1.4, 1.0e-12));
    for (const auto &patch : patches)
    {
      CHECK(patch.weight == 0.0);
      CHECK(patch.provenance.quadrature_weight == 0.0);
      CHECK(patch.longitudinal_cell == std::array<double, 2>{0.0, 0.0});
    }
    CHECK(ownership.cells[0].owners == std::vector<std::size_t>{0});
    CHECK_THAT(ownership.cells[0].attributed.front(), WithinAbs(1.0, 1.0e-12));
    CHECK_THAT(ownership.cells[1].owned_length, WithinAbs(0.4, 1.0e-12));
  }
  SECTION("a straddling cell keeps exactly its fraction outside, as a cell of the kept "
          "interval")
  {
    std::vector<Patch> patches = {Cell(4, 0, 0.0, 1.0), Cell(4, 0, 1.0, 2.5),
                                  Cell(4, 0, 2.5, 4.0)};
    const auto before = patches;
    const auto ownership = ApplyContinuationOwnership(patches, {box}, 3, tolerance, R);
    REQUIRE(ownership.cells.size() == 2);
    CHECK(ownership.wholly_owned_cells == 1);
    CHECK(ownership.clipped_cells == 1);
    CHECK_THAT(ownership.owned_length, WithinAbs(1.5, 1.0e-12));
    CHECK(ownership.cells[1].patch == 1);
    CHECK_THAT(ownership.cells[1].owned_length, WithinAbs(0.5, 1.0e-12));
    CHECK_THAT(ownership.cells[1].cell_length, WithinAbs(1.5, 1.0e-12));
    // The clipped cell [1, 2.5] keeps [1.5, 2.5]: weight x 1 / 1.5, origin at x = 2, the
    // cell symmetric about it, the quadrature weight scaled alike, the anchor moved with
    // the origin — identical to a cell built directly on [1.5, 2.5].
    Same(patches[1], Cell(4, 0, 1.5, 2.5));
    CHECK_THAT(patches[1].weight, WithinAbs(before[1].weight / 1.5, 1.0e-12));
    CHECK_THAT(patches[1].provenance.quadrature_weight,
               WithinAbs(before[1].provenance.quadrature_weight / 1.5, 1.0e-12));
    CHECK(patches[1].provenance.s0 == before[1].provenance.s0);
    CHECK(patches[1].provenance.s1 == before[1].provenance.s1);
    // (iv) The cell beyond the box is untouched.
    Same(patches[2], before[2]);
    // Idempotent: the placed cells are owned no further.
    auto again = patches;
    const auto repeat = ApplyContinuationOwnership(again, {box}, 3, tolerance, R);
    CHECK(repeat.owned_length == 0.0);
    Same(again[1], patches[1]);
    // (v) Energy consistency: the correction is additive in the weights, so the removed
    // weight is exactly the owned fraction of every owned cell, once.
    double removed = 0.0, expected = 0.0;
    for (std::size_t i = 0; i < patches.size(); i++)
    {
      removed += before[i].weight - patches[i].weight;
    }
    for (const auto &cell : ownership.cells)
    {
      expected += before[cell.patch].weight * cell.owned_length / cell.cell_length;
    }
    CHECK_THAT(removed, WithinAbs(expected, 1.0e-12));
    CHECK_THAT(removed, WithinAbs(0.1 * 1.5, 1.0e-12));
  }
  SECTION("a foreign cell inside the box is unchanged")
  {
    // Off every claim (y = 2, x in [-2, -0.8]) and alongside the claim x = -2.5 over its
    // own range (x = -1.8, y in [0, 1]): foreign, untouched.
    Patch beside = Cell(4, 2, 0.0, 1.0);
    beside.origin = {-1.8, 0.5, 0.0};
    beside.axis_w = {0.0, 1.0, 0.0};
    beside.axis_u = {1.0, 0.0, 0.0};
    beside.maxwell_conductor_anchors = {beside.origin};
    std::vector<Patch> patches = {Cell(4, 1, -2.0, -1.5, 2.0), Cell(4, 1, -1.5, -0.8, 2.0),
                                  beside};
    const auto before = patches;
    const auto ownership = ApplyContinuationOwnership(patches, {box}, 3, tolerance, R);
    CHECK(ownership.cells.empty());
    CHECK(ownership.owned_length == 0.0);
    for (std::size_t i = 0; i < patches.size(); i++)
    {
      Same(patches[i], before[i]);
    }
  }
  SECTION("a stack end adjacent to the cluster is owned up to the face, unchanged beyond")
  {
    // The stack's cells sit on its first side y = 0.5; the side whose own edge is the
    // claimed edge y = 0 (edge offset -0.5) continues the claim through its cut (abutting
    // x = 0): owned inside the box; the cells past the face keep their patches. A
    // perpendicular stretch starting at the claim end continues nothing and keeps its
    // first cell inside the box.
    Patch perpendicular = Cell(5, 0, 0.0, 1.0);
    perpendicular.origin = {0.0, 0.5, 0.0};
    perpendicular.axis_w = {0.0, 1.0, 0.0};
    perpendicular.axis_u = {1.0, 0.0, 0.0};
    perpendicular.maxwell_conductor_anchors = {perpendicular.origin};
    std::vector<Patch> patches = {
        Cell(4, 0, 0.0, 1.0, 0.5, -0.5), Cell(4, 0, 1.0, 2.5, 0.5, -0.5),
        Cell(4, 0, 2.5, 4.0, 0.5, -0.5), Cell(4, 0, 4.0, 5.0, 0.5, -0.5), perpendicular};
    const auto before = patches;
    const auto ownership = ApplyContinuationOwnership(patches, {box}, 3, tolerance, R);
    REQUIRE(ownership.cells.size() == 2);
    CHECK_THAT(ownership.owned_length, WithinAbs(1.5, 1.0e-12));
    CHECK(patches[0].weight == 0.0);
    Same(patches[1], Cell(4, 0, 1.5, 2.5, 0.5, -0.5));
    Same(patches[2], before[2]);
    Same(patches[3], before[3]);
    Same(patches[4], before[4]);
  }
  SECTION("a pair with one edge claimed: the claimed side owned, the other side untouched")
  {
    // A pair of separation 1 whose both sides' cells sit on the midline y = 0.5 beyond the
    // claim end x = 0: side 0's own edge is the claimed edge y = 0 (segment 7 continues
    // segment 3 at the mesh vertex x = 0; edge offset -0.5), side 1's own edge y = 1
    // (segment 8, edge offset +0.5) is not claimed. Decision 252: side 0 is owned inside
    // the box (its first cell wholly, its second clipped at the face x = 1.5), side 1's
    // cells stay whole although they abut the same claim end within R — the coupon's
    // twins do not carry that edge. The removed weight is side 0's owned fraction only.
    std::vector<Patch> patches = {
        Cell(4, 0, 0.0, 1.0, 0.5, -0.5), Cell(4, 0, 1.0, 2.5, 0.5, -0.5),
        Cell(4, 1, 0.0, 1.0, 0.5, 0.5), Cell(4, 1, 1.0, 2.5, 0.5, 0.5)};
    patches[2].provenance.segment = patches[3].provenance.segment = 8;
    const auto before = patches;
    const auto ownership = ApplyContinuationOwnership(patches, {box}, 3, tolerance, R);
    REQUIRE(ownership.cells.size() == 2);
    CHECK(ownership.cells[0].patch == 0);
    CHECK(ownership.cells[1].patch == 1);
    CHECK(ownership.wholly_owned_cells == 1);
    CHECK(ownership.clipped_cells == 1);
    CHECK_THAT(ownership.owned_length, WithinAbs(1.5, 1.0e-12));
    CHECK_THAT(ownership.owned_by_stretch.at(std::make_tuple(4, 0, std::size_t{0})),
               WithinAbs(1.5, 1.0e-12));
    CHECK(ownership.owned_by_stretch.count(std::make_tuple(4, 1, std::size_t{0})) == 0);
    CHECK(patches[0].weight == 0.0);
    Same(patches[1], Cell(4, 0, 1.5, 2.5, 0.5, -0.5));
    Same(patches[2], before[2]);
    Same(patches[3], before[3]);
    double removed = 0.0;
    for (std::size_t i = 0; i < patches.size(); i++)
    {
      removed += before[i].weight - patches[i].weight;
    }
    CHECK_THAT(removed, WithinAbs(0.1 * 1.5, 1.0e-12));
    // The record agrees (one classifier): side 0 wholly inside is a Continuation, side 1
    // Foreign.
    const std::vector<Patch> inside = {Cell(4, 0, 0.0, 1.0, 0.5, -0.5),
                                       Cell(4, 1, 0.0, 1.0, 0.5, 0.5)};
    const auto records =
        FindTranslationalStretchInsideSpatialSupport(inside, {box}, 3, tolerance);
    REQUIRE(records.size() == 2);
    CHECK(records[0].stretch == 0);
    CHECK(records[0].continuation);
    CHECK(records[1].stretch == 1);
    CHECK(!records[1].continuation);
    // The segment branch is the own edge too: side 1 on the claimed segment 3 (a claim cut
    // mid-segment) is owned whatever its offset.
    std::vector<Patch> on_claimed_segment = {Cell(4, 1, 0.0, 1.0, 0.5, 0.5)};
    on_claimed_segment[0].provenance.segment = 3;
    CHECK_THAT(
        ApplyContinuationOwnership(on_claimed_segment, {box}, 3, tolerance, R).owned_length,
        WithinAbs(1.0, 1.0e-12));
  }
  SECTION("fresh and cache-round-tripped patches classify identically")
  {
    // The geometry cache (version 10) carries the mesh segment and the own-edge offset of
    // every patch with the cluster patch's claims, support box and chain: the ownership on
    // the cached patches is
    // the ownership on the fresh ones (a cache without the segment would send every
    // stretch through the abutment branch; without the offset every pair side would be
    // judged on its midline).
    test::SharedTempDir temp;
    config::ElectrostaticSolverData::ResponseCorrectionData data;
    data.matching_radius = R;
    data.models.emplace_back().idx = 1;
    data.models.back().name = "gap";
    data.models.back().fabricated_matrix = "fabricated.csv";
    data.models.back().thin_matrix = "thin.csv";
    data.models.back().basis_points = "points.csv";
    data.models.emplace_back().idx = 2;
    data.models.back().name = "cluster";
    data.models.back().fabricated_matrix = "fabricated.csv";
    data.models.back().thin_matrix = "thin.csv";
    data.models.back().basis_points = "points.csv";
    data.models.back().spatial_basis = true;
    data.patches = {Cell(4, 0, 0.0, 1.0, 0.5, -0.5), Cell(4, 0, 1.0, 2.5, 0.5, -0.5),
                    Cell(4, 1, 0.0, 1.0, 0.5, 0.5),  Cell(4, 1, 1.0, 2.5, 0.5, 0.5),
                    Cell(4, 2, 0.3, 0.8, 0.5, 0.5),  Patch{}};
    for (auto &patch : data.patches)
    {
      patch.model = 1;
    }
    data.patches[2].provenance.segment = data.patches[3].provenance.segment = 8;
    data.patches[4].provenance.segment = 3;
    data.patches.back().model = 2;
    data.patches.back().provenance.feature = 8;
    data.patches.back().provenance.claims = box.claims;
    // A contract-3 model's support box and chain (rule B4) ride along in the cache (v7).
    data.patches.back().provenance.has_support_box = true;
    data.patches.back().provenance.support_box = {-2.0, -2.0, 1.0, 2.0};
    data.patches.back().provenance.chain = {{0.0, 0.0, 1.0, 0.0}, {-1.5, 0.0, -1.5, 1.0}};
    // A legacy-contract alias record (USER decision 283) rides along in the cache.
    data.legacy_contract = {{"legacy-model",
                             std::string(64, 'a'),
                             std::string(64, 'b'),
                             "unit test: USER decision 283",
                             {8, 12}}};
    const auto cache_path = temp.temp_dir / "response-geometry-ownership.json";
    WriteResponseGeometryCache(cache_path, data);
    const auto cached = ReadResponseGeometryCache(cache_path, data);
    REQUIRE(cached.patches.size() == data.patches.size());
    REQUIRE(cached.legacy_contract.size() == 1);
    CHECK(cached.legacy_contract.front().model == "legacy-model");
    CHECK(cached.legacy_contract.front().key == std::string(64, 'a'));
    CHECK(cached.legacy_contract.front().context_digest == std::string(64, 'b'));
    CHECK(cached.legacy_contract.front().reason == "unit test: USER decision 283");
    CHECK(cached.legacy_contract.front().features == std::vector<int>{8, 12});
    for (std::size_t i = 0; i < data.patches.size(); i++)
    {
      CHECK(cached.patches[i].provenance.feature == data.patches[i].provenance.feature);
      CHECK(cached.patches[i].provenance.segment == data.patches[i].provenance.segment);
      CHECK(cached.patches[i].provenance.stretch == data.patches[i].provenance.stretch);
      CHECK(cached.patches[i].provenance.edge_offset ==
            data.patches[i].provenance.edge_offset);
      CHECK(cached.patches[i].provenance.claims.size() ==
            data.patches[i].provenance.claims.size());
      CHECK(cached.patches[i].provenance.has_support_box ==
            data.patches[i].provenance.has_support_box);
      CHECK(cached.patches[i].provenance.support_box ==
            data.patches[i].provenance.support_box);
      CHECK(cached.patches[i].provenance.chain == data.patches[i].provenance.chain);
    }
    const auto fresh_records =
        FindTranslationalStretchInsideSpatialSupport(data.patches, {box}, 3, tolerance);
    const auto cached_records =
        FindTranslationalStretchInsideSpatialSupport(cached.patches, {box}, 3, tolerance);
    REQUIRE(fresh_records.size() == 1);  // stretch 2 wholly inside, on the claimed segment
    CHECK(fresh_records.front().continuation);
    REQUIRE(cached_records.size() == fresh_records.size());
    CHECK(cached_records.front().continuation == fresh_records.front().continuation);
    auto fresh = data.patches, reloaded = cached.patches;
    const auto fresh_ownership = ApplyContinuationOwnership(fresh, {box}, 3, tolerance, R);
    const auto cached_ownership =
        ApplyContinuationOwnership(reloaded, {box}, 3, tolerance, R);
    CHECK_THAT(fresh_ownership.owned_length, WithinAbs(1.5 + 0.5, 1.0e-12));
    CHECK_THAT(cached_ownership.owned_length,
               WithinAbs(fresh_ownership.owned_length, 1.0e-12));
    REQUIRE(cached_ownership.cells.size() == fresh_ownership.cells.size());
    for (std::size_t i = 0; i < fresh_ownership.cells.size(); i++)
    {
      CHECK(cached_ownership.cells[i].patch == fresh_ownership.cells[i].patch);
      CHECK_THAT(cached_ownership.cells[i].owned_length,
                 WithinAbs(fresh_ownership.cells[i].owned_length, 1.0e-12));
    }
    for (std::size_t i = 0; i < fresh.size(); i++)
    {
      CHECK_THAT(reloaded[i].weight, WithinAbs(fresh[i].weight, 1.0e-12));
    }
    // A version-6 cache (no support box, no chain) is refused.
    std::ifstream input(cache_path);
    nlohmann::json stale = nlohmann::json::parse(input);
    input.close();
    CHECK(stale["Version"] == 14);
    stale["Version"] = 6;
    const auto stale_path = temp.temp_dir / "response-geometry-ownership-stale.json";
    {
      std::ofstream output(stale_path);
      output << stale.dump(2) << "\n";
    }
    CHECK_THROWS_WITH(ReadResponseGeometryCache(stale_path, data),
                      Catch::Matchers::ContainsSubstring("cache version 6"));
  }
  SECTION("vertex ownership (rule B4): a corner patch on a chain piece end inside a "
          "contract-3 box is owned once; interior chain points, other planes, other boxes "
          "and legacy models are not")
  {
    // The contract-3 support in its local frame (units of R = 1.5): box x in [-2, 1],
    // y in [-2, 2]; the chain runs from the claim cut (0, 0) along +x to the device corner
    // (0.6, 0) and up to the face (0.6, 2). The spatial patch frame is the identity about
    // the origin; a second support (patch 1) has its frame origin 0.9 R further along x
    // and its own chain reaching the SAME device corner (local (-0.3, 0)), so the corner
    // lies on a chain end inside both boxes.
    auto Spatial = [&](double shift_x, std::vector<std::array<double, 4>> chain)
    {
      Patch patch;
      patch.origin = {shift_x * R, 0.0, 0.0};
      patch.axis_u = {1.0, 0.0, 0.0};
      patch.axis_v = {0.0, 1.0, 0.0};
      patch.axis_w = {0.0, 0.0, 1.0};
      patch.weight = 1.0;
      patch.provenance.feature = 8;
      patch.provenance.claims = box.claims;
      patch.provenance.has_support_box = true;
      patch.provenance.support_box = {-2.0, -2.0, 1.0, 2.0};
      patch.provenance.chain = std::move(chain);
      return patch;
    };
    const std::vector<std::array<double, 4>> chain_a = {{0.0, 0.0, 0.6, 0.0},
                                                        {0.6, 0.0, 0.6, 2.0}};
    const std::vector<std::array<double, 4>> chain_b = {{-0.9, 0.0, -0.3, 0.0},
                                                        {-0.3, 0.0, -0.3, 2.0}};
    auto Corner = [&](int feature, double x, double y, double z = 0.0)
    {
      Patch patch;
      patch.origin = {x * R, y * R, z * R};
      patch.axis_u = {1.0, 0.0, 0.0};
      patch.axis_v = {0.0, 1.0, 0.0};
      patch.axis_w = {0.0, 0.0, 1.0};
      patch.weight = 1.0;
      patch.provenance.feature = feature;
      return patch;
    };
    std::vector<Patch> patches = {Spatial(0.0, chain_a),     Spatial(0.9, chain_b),
                                  Corner(20, 0.6, 0.0),       // the chain corner (both)
                                  Corner(21, 0.3, 0.0),       // on a chain piece interior
                                  Corner(22, 0.6, 0.0, 0.5),  // another plane
                                  Corner(23, 0.6, 2.4),       // outside both boxes
                                  Corner(24, 1.5, 0.0)};      // inside box B, off its chain
    std::vector<SpatialSupportBounds> supports;
    for (std::size_t i = 0; i < 2; i++)
    {
      SpatialSupportBounds support = box;
      support.patch = i;
      support.has_support_box = true;
      support.support_box = patches[i].provenance.support_box;
      support.chain = patches[i].provenance.chain;
      supports.push_back(support);
    }
    const auto ownership = ApplyContinuationOwnership(patches, supports, 3, tolerance, R);
    REQUIRE(ownership.vertices.size() == 1);
    CHECK(ownership.shared_vertices == 1);
    const auto &shared = ownership.vertices[0];
    CHECK(shared.patch == 2);
    CHECK(shared.feature == 20);
    CHECK(shared.owners == std::vector<std::size_t>{0, 1});
    CHECK_THAT(shared.face_distance_over_r, WithinAbs(0.4, 1.0e-12));  // to x = 1
    CHECK(shared.arm_outside_box);
    CHECK_THAT(shared.lost_arm_length_over_r, WithinAbs(0.6, 1.0e-12));  // R - 0.4 R
    CHECK_THAT(patches[2].weight, WithinAbs(0.0, 1.0e-15));
    for (const std::size_t untouched : {3, 4, 5, 6})
    {
      CHECK_THAT(patches[untouched].weight, WithinAbs(1.0, 1.0e-15));
    }
    // Idempotent: a second pass owns nothing more (weight-0 patches are skipped).
    const auto repeat = ApplyContinuationOwnership(patches, supports, 3, tolerance, R);
    CHECK(repeat.vertices.empty());
    // A legacy (claims-only) support owns no vertex.
    std::vector<Patch> legacy = {Spatial(0.0, chain_a), Corner(20, 0.6, 0.0)};
    legacy[0].provenance.has_support_box = false;
    legacy[0].provenance.chain.clear();
    std::vector<SpatialSupportBounds> legacy_supports = {box};
    CHECK(ApplyContinuationOwnership(legacy, legacy_supports, 3, tolerance, R)
              .vertices.empty());
    CHECK_THAT(legacy[1].weight, WithinAbs(1.0, 1.0e-15));
    // The Diagnostics entry lists the records with Kind Vertex.
    config::ElectrostaticSolverData::ResponseCorrectionData data;
    data.matching_radius = R;
    data.models.emplace_back().idx = 0;
    data.models.back().name = "cluster";
    data.models.back().spatial_basis = true;
    data.models.emplace_back().idx = 1;
    data.models.back().name = "corner";
    for (auto &patch : patches)
    {
      patch.model = patch.provenance.has_support_box ? 0 : 1;
    }
    data.patches = patches;
    const auto record = DescribeContinuationOwnership(ownership, supports, data, 1.0);
    CHECK(record["Vertices"]["Count"] == 1);
    CHECK(record["Vertices"]["Shared"] == 1);
    CHECK(record["Vertices"]["Records"][0]["Kind"] == "Vertex");
    CHECK(record["Vertices"]["Records"][0]["Owners"].size() == 2);
    CHECK(record["Vertices"]["Records"][0]["ArmOutsideBox"] == true);
    CHECK_THAT(record["Vertices"]["Records"][0]["LostArmLengthOverR"].get<double>(),
               WithinAbs(0.6, 1.0e-12));
  }
  SECTION("a cell on the continuations of two coupons is removed once, attributed by the "
          "midpoint between the claim ends")
  {
    // A second cluster whose box spans x in [0.5, 6] claims the edge y = 0 from x = 2 to 4
    // (segment 9): the bridge x in [0, 2] continues both claims (A's end x = 0 and B's end
    // x = 2); the cell [1, 2] lies inside both boxes, the cell [0, 1] inside A wholly and
    // inside B over [0.5, 1]. Both are removed wholly (the complement of every owning box),
    // once; the attribution splits at the midpoint x = 1 between the claim ends: A owns
    // [0, 1], B owns [1, 2]. The cell [4, 5] beyond B's claim is B's continuation on the
    // far side, inside B's box only.
    SpatialSupportBounds other;
    other.patch = 9;
    other.min = {0.5, -3.0, -2.0};
    other.max = {6.0, 3.0, 2.0};
    other.claims = {Claim{9, {2.0, 0.0, 0.0}, {4.0, 0.0, 0.0}},
                    Claim{11, {4.0, 0.0, 0.0}, {4.0, 1.0, 0.0}}};
    std::vector<Patch> patches = {Cell(4, 0, 0.0, 1.0), Cell(4, 0, 1.0, 2.0),
                                  Cell(6, 0, 4.0, 5.0), Cell(6, 0, 5.0, 7.0)};
    const auto before = patches;
    const auto ownership =
        ApplyContinuationOwnership(patches, {box, other}, 3, tolerance, R);
    REQUIRE(ownership.cells.size() == 4);
    CHECK(ownership.shared_cells == 2);
    CHECK_THAT(ownership.shared_length, WithinAbs(2.0, 1.0e-12));
    CHECK_THAT(ownership.owned_length, WithinAbs(2.0 + 1.0 + 1.0, 1.0e-12));
    CHECK(patches[0].weight == 0.0);
    CHECK(patches[1].weight == 0.0);
    CHECK(patches[2].weight == 0.0);
    Same(patches[3], Cell(6, 0, 6.0, 7.0));
    CHECK(ownership.cells[0].owners == std::vector<std::size_t>{0, 9});
    CHECK_THAT(ownership.cells[0].attributed[0], WithinAbs(1.0, 1.0e-12));
    CHECK_THAT(ownership.cells[0].attributed[1], WithinAbs(0.0, 1.0e-12));
    CHECK(ownership.cells[1].owners == std::vector<std::size_t>{0, 9});
    CHECK_THAT(ownership.cells[1].attributed[0], WithinAbs(0.0, 1.0e-12));
    CHECK_THAT(ownership.cells[1].attributed[1], WithinAbs(1.0, 1.0e-12));
    CHECK_THAT(ownership.owned_by_support.at(0), WithinAbs(1.0, 1.0e-12));
    CHECK_THAT(ownership.owned_by_support.at(9), WithinAbs(1.0 + 2.0, 1.0e-12));
    // The removed weight equals the owned fraction of every owned cell, once.
    double removed = 0.0, expected = 0.0;
    for (std::size_t i = 0; i < patches.size(); i++)
    {
      removed += before[i].weight - patches[i].weight;
    }
    for (const auto &cell : ownership.cells)
    {
      expected += before[cell.patch].weight * cell.owned_length / cell.cell_length;
    }
    CHECK_THAT(removed, WithinAbs(expected, 1.0e-12));
    CHECK_THAT(removed, WithinAbs(0.1 * 4.0, 1.0e-12));
  }
  SECTION("the stretch record carries the owned length per record and in total")
  {
    // The 0.9-long stretch 1 lies wholly inside the box on the claim's own segment 3 (a
    // claim cut mid-segment); stretch 0 straddles the face. The spatial patch is patch 4.
    config::ElectrostaticSolverData::ResponseCorrectionData data;
    data.models.emplace_back().idx = 1;
    data.models.back().name = "isolated-edge";
    data.models.emplace_back().idx = 2;
    data.models.back().name = "cluster";
    data.models.back().spatial_basis = true;
    data.patches = {Cell(4, 0, 0.0, 1.0), Cell(4, 0, 1.0, 2.5), Cell(4, 1, 0.3, 0.8, 0.5),
                    Cell(4, 1, 0.8, 1.2, 0.5), Patch{}};
    for (auto &patch : data.patches)
    {
      patch.model = 1;
    }
    data.patches[2].provenance.segment = data.patches[3].provenance.segment = 3;
    data.patches.back().model = 2;
    data.patches.back().provenance.feature = 8;
    SpatialSupportBounds spatial = box;
    spatial.patch = 4;
    const auto records =
        FindTranslationalStretchInsideSpatialSupport(data.patches, {spatial}, 3, tolerance);
    REQUIRE(records.size() == 1);
    CHECK(records.front().continuation);
    const auto ownership =
        ApplyContinuationOwnership(data.patches, {spatial}, 3, tolerance, R);
    const auto description = DescribeTranslationalOwnershipRecords(records, {spatial}, data,
                                                                   2.0, {}, &ownership);
    CHECK_THAT(description["OwnedLength"].get<double>(), WithinAbs(2.0 * 0.9, 1.0e-12));
    CHECK_THAT(description["OwnedLengthTotal"].get<double>(),
               WithinAbs(2.0 * (0.9 + 1.5), 1.0e-12));
    CHECK_THAT(description["Records"][0]["OwnedLength"].get<double>(),
               WithinAbs(2.0 * 0.9, 1.0e-12));
    const auto summary = DescribeContinuationOwnership(ownership, {spatial}, data, 2.0);
    CHECK(summary["Cells"] == 4);
    CHECK(summary["WhollyOwnedCells"] == 3);
    CHECK(summary["ClippedCells"] == 1);
    CHECK_THAT(summary["OwnedLength"].get<double>(), WithinAbs(2.0 * 2.4, 1.0e-12));
    CHECK(summary["BySupport"][0]["SpatialPatch"] == 4);
    CHECK(summary["BySupport"][0]["SpatialModel"] == "cluster");
    CHECK_THAT(summary["BySupport"][0]["OwnedLength"].get<double>(),
               WithinAbs(2.0 * 2.4, 1.0e-12));
    CHECK(summary["OwnedCells"][1]["Owners"][0]["SpatialFeature"] == 8);
    CHECK(DescribeContinuationOwnershipSummary(summary).find("4 translational cell(s)") !=
          std::string::npos);
  }
}

TEST_CASE("SurfaceResponseOperatorSpatialSupportMarginOverlaps",
          "[surfaceresponseoperator][Serial]")
{
  // Decision 244: the two 19-edge coupons of the S1p re-preflight (R = 1.9, boxes
  // x in [568.55, 592.75] and [583.05, 607.25], the same y and z): A's claim on
  // y = -108.325 ends at its cut x = 587.05 and continues to A's face 592.75 over 2.80 of
  // B's claim [589.95, 593.0]; B's cut at 589.95 continues to B's face 583.05 over 3.80 of
  // A's claims [583.25, 587.05]; the two continuations share the bridge [587.05, 589.95]
  // (2.90): 9.50 of edge corrected by both coupons, recorded (no claim of either inside
  // the other's claims hull). A claim of B reaching into A's hull is a true overlap.
  using Patch = config::ElectrostaticSolverData::ResponseCorrectionPatchData;
  using Claim = Patch::Provenance::Claim;
  const double R = 1.9;
  const double tolerance = kSignatureParameterToleranceOverRadius * R;
  const double y = -108.325;
  SpatialSupportBounds a, b;
  a.patch = 204;
  a.min = {568.55, -128.2, 2.85};
  a.max = {592.75, -102.625, 6.8};
  a.claims = {Claim{494, {583.25, y, 4.8}, {583.55, y, 4.8}},
              Claim{514, {583.55, y, 4.8}, {587.05, y, 4.8}},
              Claim{520, {583.25, y, 4.8}, {583.25, -120.0, 4.8}}};
  b.patch = 205;
  b.min = {583.05, -128.2, 2.85};
  b.max = {607.25, -102.625, 6.8};
  b.claims = {Claim{646, {589.95, y, 4.8}, {593.0, y, 4.8}},
              Claim{650, {593.0, y, 4.8}, {593.0, -120.0, 4.8}}};
  SECTION("a margins-only overlap is recorded with the S1p pair's 9.50 um")
  {
    const auto overlaps = FindSpatialSupportMarginOverlaps({a, b}, 3, tolerance);
    REQUIRE(overlaps.size() == 1);
    const auto &overlap = overlaps.front();
    CHECK(overlap.first_patch == 204);
    CHECK(overlap.second_patch == 205);
    CHECK(!overlap.claim_in_hull);
    CHECK_THAT(overlap.overlap_min[0], WithinAbs(583.05, 1.0e-12));
    CHECK_THAT(overlap.overlap_max[0], WithinAbs(592.75, 1.0e-12));
    CHECK_THAT(overlap.first_margin_over_second_claims, WithinAbs(2.80, 1.0e-9));
    CHECK_THAT(overlap.second_margin_over_first_claims, WithinAbs(3.80, 1.0e-9));
    CHECK_THAT(overlap.margin_over_margin, WithinAbs(2.90, 1.0e-9));
    config::ElectrostaticSolverData::ResponseCorrectionData data;
    data.models.emplace_back().idx = 1;
    data.models.back().name = "spatialedgecluster_edgecount-19";
    data.models.back().spatial_basis = true;
    data.patches.resize(206);
    for (auto &patch : data.patches)
    {
      patch.model = 1;
    }
    const auto description = DescribeSpatialSupportMarginOverlaps(overlaps, data, 1.0);
    CHECK(description["Count"] == 1);
    CHECK(description["ClaimInHull"] == 0);
    CHECK_THAT(description["DoubleCountedLength"].get<double>(), WithinAbs(9.50, 1.0e-9));
    CHECK_THAT(description["MarginOverClaimsLength"].get<double>(),
               WithinAbs(6.60, 1.0e-9));
    CHECK_THAT(description["MarginOverMarginLength"].get<double>(),
               WithinAbs(2.90, 1.0e-9));
    CHECK(description["Pairs"][0]["ClaimInHull"] == false);
    CHECK(DescribeSpatialSupportMarginOverlapWarning(description).find("9.500000e+00") !=
          std::string::npos);
    // The bridge's cells (feature 26: 0.8056 + 1.2889 + 0.8056) continue both claims: owned
    // once, 1.45 attributed to each coupon by the midpoint x = 588.5 of the claim ends.
    auto Cell = [&](double begin, double end)
    {
      Patch patch;
      patch.origin = {0.5 * (begin + end), y, 4.8};
      patch.axis_u = {0.0, 1.0, 0.0};
      patch.axis_v = {0.0, 0.0, 1.0};
      patch.axis_w = {1.0, 0.0, 0.0};
      patch.longitudinal_cell = {begin - patch.origin[0], end - patch.origin[0]};
      patch.weight = end - begin;
      patch.provenance.feature = 26;
      patch.provenance.stretch = 0;
      patch.provenance.segment = 646;
      return patch;
    };
    std::vector<Patch> bridge = {Cell(587.05, 587.8556), Cell(587.8556, 589.1444),
                                 Cell(589.1444, 589.95)};
    const auto ownership = ApplyContinuationOwnership(bridge, {a, b}, 3, tolerance, R);
    CHECK(ownership.shared_cells == 3);
    CHECK_THAT(ownership.owned_length, WithinAbs(2.90, 1.0e-9));
    CHECK_THAT(ownership.owned_by_support.at(204), WithinAbs(1.45, 1.0e-9));
    CHECK_THAT(ownership.owned_by_support.at(205), WithinAbs(1.45, 1.0e-9));
    for (const auto &patch : bridge)
    {
      CHECK(patch.weight == 0.0);
    }
  }
  SECTION("a claim inside the other's claims hull is a true overlap")
  {
    // B's second claim turns into A's hull: x = 586 from y to -115 (its far end strictly
    // inside A's claims' bounding box x in [583.25, 587.05], y in [-120, -108.325]).
    SpatialSupportBounds reaching = b;
    reaching.claims = {Claim{646, {586.0, y, 4.8}, {593.0, y, 4.8}},
                       Claim{650, {586.0, y, 4.8}, {586.0, -115.0, 4.8}}};
    const auto overlaps = FindSpatialSupportMarginOverlaps({a, reaching}, 3, tolerance);
    REQUIRE(overlaps.size() == 1);
    CHECK(overlaps.front().claim_in_hull);
    // The claim SEGMENT is tested, not only its ends (decision 252): B's claim y = -115
    // from x = 582 to 588 crosses A's hull (x in [583.25, 587.05]) with both ends outside
    // it, while B's own hull (x in [582, 593], y in [-115, -108.325]) holds no end of A's
    // claims strictly inside (A's claims lie on its face y = -108.325).
    SpatialSupportBounds crossing = b;
    crossing.claims = {Claim{646, {589.95, y, 4.8}, {593.0, y, 4.8}},
                       Claim{650, {582.0, -115.0, 4.8}, {588.0, -115.0, 4.8}}};
    const auto crossed = FindSpatialSupportMarginOverlaps({a, crossing}, 3, tolerance);
    REQUIRE(crossed.size() == 1);
    CHECK(crossed.front().claim_in_hull);
    // A claim lying on A's hull face x = 587.05 (A's claims' cut) over y in [-110, -115]
    // touches the hull within the tolerance and is not inside it.
    SpatialSupportBounds grazing = b;
    grazing.claims = {Claim{646, {589.95, y, 4.8}, {593.0, y, 4.8}},
                      Claim{650, {587.05, -110.0, 4.8}, {587.05, -115.0, 4.8}}};
    const auto grazed = FindSpatialSupportMarginOverlaps({a, grazing}, 3, tolerance);
    REQUIRE(grazed.size() == 1);
    CHECK(!grazed.front().claim_in_hull);
  }
  SECTION("separate boxes and claim-less (corner) supports are not pairs")
  {
    SpatialSupportBounds apart = b;
    apart.min[0] = 592.75;
    CHECK(FindSpatialSupportMarginOverlaps({a, apart}, 3, tolerance).empty());
    SpatialSupportBounds corner = b;
    corner.claims.clear();
    CHECK(FindSpatialSupportMarginOverlaps({a, corner}, 3, tolerance).empty());
  }
}

}  // namespace palace
