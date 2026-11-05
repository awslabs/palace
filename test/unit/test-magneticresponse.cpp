// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "fixtures.hpp"
#include "surfaceresponse-fixtures.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <functional>
#include <memory>
#include <optional>
#include <string>
#include <utility>
#include <vector>
#include <fmt/format.h>
#include <mfem.hpp>
#include <catch2/catch_test_macros.hpp>
#include <nlohmann/json.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>
#include "fem/gridfunction.hpp"
#include "fem/mesh.hpp"
#include "linalg/operator.hpp"
#include "linalg/vector.hpp"
#include "models/postoperator.hpp"
#include "models/romoperator.hpp"
#include "models/spaceoperator.hpp"
#include "models/surfaceresponseoperator.hpp"
#include "utils/communication.hpp"
#include "utils/geodata.hpp"
#include "utils/iodata.hpp"
#include "utils/units.hpp"

namespace palace
{

namespace fs = std::filesystem;

using json = nlohmann::json;
using Catch::Matchers::WithinAbs;
using Catch::Matchers::WithinRel;

namespace
{

// The magnetic (inductive) leaves of the mx-sc-magnetic lane (DESIGN-v2 section 3; USER
// 1221 D3; decision 1229) on the legacy 3D IsolatedEdge library of the fixtures: the vacuum
// twin's domain matrices (4 x 4 on the 4-knot basis, the electric CSV form), the
// MagneticSource provenance and the library's Fabrication block the r2 rule compares it
// against. `defect` = S_m,fab - S_m,thin as written (S_thin = the fixture's thin matrix).
struct MagneticLibraryFiles
{
  fs::path library;
  fs::path fabricated;
  fs::path thin;
  std::array<std::array<double, 4>, 4> defect{};
};

constexpr std::array<std::array<double, 4>, 4> kThinVacuum = {
    {{{1.0e-12, 0.1e-12, 0.05e-12, 0.02e-12}},
     {{0.1e-12, 0.5e-12, 0.08e-12, 0.04e-12}},
     {{0.05e-12, 0.08e-12, 0.7e-12, 0.09e-12}},
     {{0.02e-12, 0.04e-12, 0.09e-12, 0.4e-12}}}};

// A PSD defect (diagonally dominant), the default leaf.
constexpr std::array<std::array<double, 4>, 4> kPsdDefect = {
    {{{2.0e-12, 0.4e-12, 0.1e-12, 0.1e-12}},
     {{0.4e-12, 1.5e-12, 0.2e-12, 0.1e-12}},
     {{0.1e-12, 0.2e-12, 1.8e-12, 0.3e-12}},
     {{0.1e-12, 0.1e-12, 0.3e-12, 1.1e-12}}}};

json ProcessOfRecord()
{
  return {{"MetalThickness", 0.01},      {"OveretchDepth", 0.005},
          {"SidewallAngle", 90.0},       {"TopRoundingRadius", 0.0},
          {"BottomRoundingRadius", 0.0}, {"SubstratePermittivity", 11.45}};
}

json DefaultMagneticSource()
{
  return {{"Version", 1},
          {"Generator", "unit-test"},
          {"VacuumPermittivity", 1.0},
          {"Fabrication", ProcessOfRecord()},
          {"ElectricRegenerationMaxRelativeDeviation", 0.0},
          {"BasisPointCount", 4},
          {"ConductorStateCount", 0},
          {"BasisPointsIdentical", true}};
}

// Writes the magnetic twin of `base_library` (the fixture's 3D library) with the given
// defect and MagneticSource into `directory` under `tag`; `edit` may alter the library
// before it is written (root rank only; every rank returns the paths).
MagneticLibraryFiles WriteMagneticLibrary(
    const fs::path &base_library, const fs::path &directory, const std::string &tag,
    const std::array<std::array<double, 4>, 4> &defect,
    json source = DefaultMagneticSource(), const std::function<void(json &)> &edit = {})
{
  MagneticLibraryFiles files;
  files.library = directory / ("fabrication-process-magnetic-" + tag + ".json");
  files.fabricated = directory / ("magnetic-fabricated-" + tag + ".csv");
  files.thin = directory / ("magnetic-thin-" + tag + ".csv");
  files.defect = defect;
  if (Mpi::Root(Mpi::World()))
  {
    auto write = [](const fs::path &path, const std::array<std::array<double, 4>, 4> &m)
    {
      std::ofstream output(path);
      output.precision(17);
      output << "basis_i,basis_j,Q_ij (J)\n";
      for (std::size_t i = 0; i < 4; i++)
      {
        for (std::size_t j = i; j < 4; j++)
        {
          output << i + 1 << "," << j + 1 << "," << m[i][j] << "\n";
        }
      }
    };
    std::array<std::array<double, 4>, 4> fabricated = kThinVacuum;
    for (std::size_t i = 0; i < 4; i++)
    {
      for (std::size_t j = 0; j < 4; j++)
      {
        fabricated[i][j] += defect[i][j];
      }
    }
    write(files.fabricated, fabricated);
    write(files.thin, kThinVacuum);
    std::ifstream input(base_library);
    json library = json::parse(input);
    library["Name"] = "unit-test-process-magnetic-" + tag;
    const json process = ProcessOfRecord();
    for (const auto &[key, value] : process.items())
    {
      library["Fabrication"][key] = value;
    }
    for (auto &model : library["Models"])
    {
      model["MagneticFabricatedMatrix"] = files.fabricated.string();
      model["MagneticThinMatrix"] = files.thin.string();
      model["MagneticSource"] = source;
    }
    if (edit)
    {
      edit(library);
    }
    std::ofstream output(files.library);
    output << library.dump(2) << "\n";
  }
  Mpi::Barrier(Mpi::World());
  return files;
}

// The 3D CPW Maxwell problem of the fixtures (the test mesh cpw3d-surface-nc.msh, order 1,
// the legacy translational IsolatedEdge library) on the given library.
json MaxwellCpw3dConfig(const test::SurfaceResponseFiles &files, const fs::path &library)
{
  json config = files.Cpw3dLegacyConfig();
  config["Problem"]["Type"] = "Eigenmode";
  config["Boundaries"]["Ground"]["Attributes"] = {1, 2};
  config["Boundaries"].erase("Terminal");
  config["Solver"] = {{"Order", 1},
                      {"Eigenmode", {{"Target", 1.0}}},
                      {"SurfaceResponseCorrection",
                       {{"Library", library.string()},
                        {"TargetInterfaces", {1, 2, 3}},
                        {"UnmatchedPolicy", "Error"},
                        {"PatchConstruction", "Legacy"}}}};
  return config;
}

// E = 1/2 B x r: a linear field (exactly in the lowest-order Nedelec space) with the
// constant curl B, so that the flux trace along every contour is exact and every closed
// loop of B . n vanishes.
mfem::VectorFunctionCoefficient ConstantCurlField(const std::array<double, 3> &B)
{
  return mfem::VectorFunctionCoefficient(3,
                                         [B](const mfem::Vector &x, mfem::Vector &E)
                                         {
                                           E.SetSize(3);
                                           E[0] = 0.5 * (B[1] * x[2] - B[2] * x[1]);
                                           E[1] = 0.5 * (B[2] * x[0] - B[0] * x[2]);
                                           E[2] = 0.5 * (B[0] * x[1] - B[1] * x[0]);
                                         });
}

// Phase 3 (a) (DESIGN-v2 a.2 / a.3 / a.5; decisions 1284 Q-1 / 1292): the SHORTED CPW STUB
// of the R0-a tests (test::MakeSerialShortedCpwMesh; the one-label 4-edge PEC-cluster
// model, offsets 0 / 0.125 / 0.375 / 0.5, R 0.2, the 4-knot ring with its top knots at v =
// 0.1) WITH magnetic leaves (diagonal 1e-12 thin / PSD defect), the PEC-cluster anchor =
// the ground's edge point (MAJOR-1's case). `piece_form` = Edges[].LocalPiece 1 / 2 / 2 / 3
// + three ConductorReferences (-0.05 / 0.25 / 0.55: ground-left / strip / ground-right) +
// OpenContourPaths (knot 1: 1 -> 2; knot 2: 2 -> 3; knots 4, 3: 1 -> 3), so the gate
// compares a probe with its PIECE's flux state (the a.3 control's substrate); `rules` are
// written as the library's MagneticRules. Written by the root; returns the library path.
fs::path WriteStubMagneticLibrary(const fs::path &directory,
                                  const fs::path &template_library,
                                  const std::string &suffix, bool piece_form,
                                  const json &rules = json::object())
{
  const int basis_size = piece_form ? 6 : 4;
  const auto library_path =
      directory / fmt::format("fabrication-process-magnetic-stub-{}.json", suffix);
  const auto points_path = directory / fmt::format("magstub-basis-points-{}.csv", suffix);
  const auto domain_path = directory / fmt::format("magstub-domain-{}.csv", suffix);
  const auto surface_path = directory / fmt::format("magstub-surface-{}.csv", suffix);
  const auto thin_path = directory / fmt::format("magstub-magnetic-thin-{}.csv", suffix);
  const auto fab_path = directory / fmt::format("magstub-magnetic-fab-{}.csv", suffix);
  if (Mpi::Root(Mpi::World()))
  {
    {
      std::ofstream output(points_path);
      output << "x,y,z\n0.0625,0.1,0.0\n0.4375,0.1,0.0\n0.4375,-0.1,0.0\n0.0625,-0.1,0.0\n";
    }
    auto write_diagonal = [&](const fs::path &path, double diagonal)
    {
      std::ofstream output(path);
      output.precision(17);
      output << "basis_i,basis_j,Q_ij (J)\n";
      for (int i = 1; i <= basis_size; i++)
      {
        for (int j = i; j <= basis_size; j++)
        {
          output << i << "," << j << "," << (i == j ? diagonal : 0.0) << "\n";
        }
      }
    };
    write_diagonal(domain_path, 1.0e-12);
    write_diagonal(thin_path, 1.0e-12);
    write_diagonal(fab_path, 3.0e-12);  // dS_m = 2e-12 I: PSD, Applied
    {
      std::ofstream output(surface_path);
      output << "interface,edge,R (m),basis_i,basis_j,Q_ij (J),Q_total_ij (J)\n";
      for (int edge = 1; edge <= 4; edge++)
      {
        for (int i = 1; i <= basis_size; i++)
        {
          for (int j = i; j <= basis_size; j++)
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
    library["Name"] = "unit-test-process-magnetic-stub-" + suffix;
    const json process = ProcessOfRecord();
    for (const auto &[key, value] : process.items())
    {
      library["Fabrication"][key] = value;
    }
    json model = library["Models"][0];
    model["Name"] =
        piece_form ? "cpw-one-label-four-edge-pieces" : "cpw-one-label-four-edge";
    model["Topology"] = "ParallelEdgeCluster";
    model["Interfaces"] = {{{"Type", "SA"}, {"Coupon", 1}}};
    model["BasisPoints"] = points_path.string();
    model["FabricatedMatrix"] = domain_path.string();
    model["ThinMatrix"] = domain_path.string();
    model["FabricatedSurfaceMatrix"] = surface_path.string();
    model["ThinSurfaceMatrix"] = surface_path.string();
    model["EdgeOffsetTolerance"] = 1.0e-6;
    model["Edges"] = {{{"Offset", 0.0}, {"GapDirection", 1}, {"Conductor", 1}},
                      {{"Offset", 0.125}, {"GapDirection", -1}, {"Conductor", 1}},
                      {{"Offset", 0.375}, {"GapDirection", 1}, {"Conductor", 1}},
                      {{"Offset", 0.5}, {"GapDirection", -1}, {"Conductor", 1}}};
    for (const char *key : {"Separation", "SeparationTolerance", "Angle", "AngleTolerance",
                            "TraceMesh", "TraceMeshTriangles", "ConductorReferences",
                            "Reference", "OpenContourPaths", "ContourGroups"})
    {
      model.erase(key);
    }
    if (piece_form)
    {
      const int pieces[4] = {1, 2, 2, 3};
      for (int k = 0; k < 4; k++)
      {
        model["Edges"][k]["LocalPiece"] = pieces[k];
      }
      model["ConductorReferences"] = {
          {-0.05, 0.0, 0.0}, {0.25, 0.0, 0.0}, {0.55, 0.0, 0.0}};
      model["OpenContourPaths"] = {
          {{"Indices", {1}}, {"StartConductor", 1}, {"EndConductor", 2}},
          {{"Indices", {2}}, {"StartConductor", 2}, {"EndConductor", 3}},
          {{"Indices", {4, 3}}, {"StartConductor", 1}, {"EndConductor", 3}}};
    }
    else
    {
      model["Reference"] = {0.25, 0.0, 0.0};
    }
    json source = DefaultMagneticSource();
    source["ConductorStateCount"] = piece_form ? 2 : 0;
    model["MagneticFabricatedMatrix"] = fab_path.string();
    model["MagneticThinMatrix"] = thin_path.string();
    model["MagneticSource"] = source;
    library["Models"] = {model};
    if (!rules.empty())
    {
      library["MagneticRules"] = rules;
    }
    std::ofstream output(library_path);
    output << library.dump(2) << "\n";
  }
  Mpi::Barrier(Mpi::World());
  return library_path;
}

// The stub's Maxwell config on `library` (the R0-a Maxwell form: the automatic patches,
// collocated trace coupling, the short's corner cluster unmatched = raw).
json StubMaxwellConfig(const test::SurfaceResponseFiles &files, const fs::path &library)
{
  json config = files.ConvexMaxwellIslandConfig();
  config["Boundaries"]["Ground"]["Attributes"] = {1, 2, 3, 4, 5, 6, 9};
  auto &correction = config["Solver"]["SurfaceResponseCorrection"];
  correction.erase("PatchConstruction");
  correction["Library"] = library.string();
  correction["UnmatchedPolicy"] = "Warn";
  correction["TraceCoupling"] = "Collocated";
  return config;
}

// The stub's GAP-FLUX field: E = f(z) x_hat (the LONGITUDINAL component) with f = 0 on the
// grounds, rising linearly across slot A (z in [0.75, 0.875]) to V on the strip and falling
// across slot B (z in [1.125, 1.25]); the kinks lie on element faces (h = 0.125), so the
// field is exact in ND1 on the hexahedra and B = curl E = f'(z) y_hat is the gap flux
// (normal to the sheet, in the slots only; zero elsewhere). int (w x B) . dl = f(end) -
// f(start): a probe on the strip reads V against the ground anchor; the piece-form strip
// state reads V too.
mfem::VectorFunctionCoefficient GapFluxField(double voltage)
{
  return mfem::VectorFunctionCoefficient(3,
                                         [voltage](const mfem::Vector &x, mfem::Vector &E)
                                         {
                                           E.SetSize(3);
                                           E = 0.0;
                                           const double z = x[2];
                                           double f = 0.0;
                                           if (z > 0.75 && z < 0.875)
                                           {
                                             f = voltage * (z - 0.75) / 0.125;
                                           }
                                           else if (z >= 0.875 && z <= 1.125)
                                           {
                                             f = voltage;
                                           }
                                           else if (z > 1.125 && z < 1.25)
                                           {
                                             f = voltage * (1.25 - z) / 0.125;
                                           }
                                           E[0] = f;
                                         });
}

ComplexVector ProjectField(SpaceOperator &space,
                           mfem::VectorFunctionCoefficient &coefficient)
{
  GridFunction field(space.GetNDSpace(), true);
  field.Real().ProjectCoefficient(coefficient);
  field.Imag() = 0.0;
  ComplexVector E(field.Real().ParFESpace()->GetTrueVSize());
  field.Real().GetTrueDofs(E.Real());
  field.Imag().GetTrueDofs(E.Imag());
  return E;
}

}  // namespace

TEST_CASE_METHOD(test::SurfaceResponseFiles,
                 "SurfaceResponseOperator magnetic leaves: flux trace, dK, dual loop gate",
                 "[surfaceresponseoperator][3d][maxwell][magnetic][Serial][Parallel]")
{
#if !defined(MFEM_USE_GSLIB)
  SKIP("SurfaceResponseOperator requires MFEM_USE_GSLIB");
#else
  const auto files =
      WriteMagneticLibrary(library_3d_path, temp.temp_dir, "psd", kPsdDefect);
  IoData iodata(MaxwellCpw3dConfig(*this, files.library), false);
  auto device_mesh = mesh::ReadMesh(iodata, Mpi::World());
  std::vector<std::unique_ptr<Mesh>> meshes;
  meshes.push_back(std::make_unique<Mesh>(std::move(device_mesh)));
  SpaceOperator space(iodata, meshes);
  SurfaceResponseOperator response(iodata, space);
  REQUIRE(response.GetPatchCount() > 0);
  REQUIRE(response.HasMagneticModels());
  for (const auto &record : response.GetMagneticModelRecords())
  {
    INFO(record.name << ": " << record.status << " " << record.rule << " "
                     << record.reason);
    CHECK(record.status == "Applied");
    CHECK(record.psd_ratio < 1.0e-12);
    CHECK(record.negative_eigenvalue_count == 0);
    CHECK(record.same_mesh_deviation == 0.0);
  }

  // The constant-curl field: the flux trace of a knot relative to the next is
  // -(w x B) . (q - p) (w = u x v the plane normal of the patch), as the electric trace of
  // the constant field E0 = w x B would be.
  // The exact linear field reads the exact flux trace (the knot differences below); the
  // operator and the postprocessing read the field with its PEC (sheet) dofs zeroed, as dK
  // does internally, so the quadratic forms are compared on that zeroed field.
  const std::array<double, 3> B = {0.3, -0.5, 0.8};
  auto coefficient = ConstantCurlField(B);
  GridFunction field(space.GetNDSpace(), true);
  field.Real().ProjectCoefficient(coefficient);
  field.Imag() = 0.0;
  Vector x_exact, x;
  field.Real().GetTrueDofs(x_exact);
  x = x_exact;
  x.SetSubVector(space.GetNDDbcTDofLists().back(), 0.0);
  field.Real().SetFromTrueDofs(x);
  const auto exact_traces = response.GetMagneticPatchTraces(x_exact);
  const auto traces = response.GetMagneticPatchTraces(x);
  const auto &assignments = response.GetPatchAssignments();
  const std::array<std::array<double, 2>, 4> knots = {
      {{-0.08, -0.06}, {0.08, -0.06}, {0.08, 0.06}, {-0.08, 0.06}}};
  // The scale of the loaded defect (the electric CSV convention, 1/2 a^T C a = the
  // nondimensional energy): C = 2 V_scale^2 / E_scale x Q.
  const double scale =
      2.0 * std::pow(iodata.units.GetScaleFactor<Units::ValueType::VOLTAGE>(), 2) /
      iodata.units.GetScaleFactor<Units::ValueType::ENERGY>();
  double expected_form = 0.0;
  double max_difference_error = 0.0;
  if (Mpi::Root(Mpi::World()))
  {
    REQUIRE(traces.size() == assignments.size());
    REQUIRE(exact_traces.size() == traces.size());
    for (std::size_t t = 0; t < traces.size(); t++)
    {
      const auto &entry = traces[t];
      const auto &exact = exact_traces[t];
      REQUIRE(exact.patch == entry.patch);
      const auto assignment =
          std::find_if(assignments.begin(), assignments.end(),
                       [&](const auto &a) { return a.global_index == entry.patch; });
      REQUIRE(assignment != assignments.end());
      REQUIRE(entry.coefficients.size() == 4);
      const auto &u = assignment->axis_u;
      const auto &v = assignment->axis_v;
      const std::array<double, 3> w = {u[1] * v[2] - u[2] * v[1], u[2] * v[0] - u[0] * v[2],
                                       u[0] * v[1] - u[1] * v[0]};
      const std::array<double, 3> g = {w[1] * B[2] - w[2] * B[1], w[2] * B[0] - w[0] * B[2],
                                       w[0] * B[1] - w[1] * B[0]};
      for (int i = 0; i < 4; i++)
      {
        const int next = (i + 1) % 4;
        double expected = 0.0;
        for (int d = 0; d < 3; d++)
        {
          expected -= g[d] * ((knots[next][0] - knots[i][0]) * u[d] +
                              (knots[next][1] - knots[i][1]) * v[d]);
        }
        max_difference_error =
            std::max(max_difference_error,
                     std::abs(exact.coefficients[next] - exact.coefficients[i] - expected));
      }
      for (int i = 0; i < 4; i++)
      {
        for (int j = 0; j < 4; j++)
        {
          expected_form += assignment->weight * scale * entry.coefficients[i] *
                           files.defect[i][j] * entry.coefficients[j];
        }
      }
    }
    CHECK(max_difference_error < 1.0e-10);
  }
  Mpi::Broadcast(1, &expected_form, 0, Mpi::World());

  // dK = P_m^T W dS_m P_m: the quadratic form on the flux trace, symmetric, PSD, and zero
  // on a curl-free field.
  SurfaceResponseStiffnessOperator dK(response);
  Vector y(x.Size());
  dK.Mult(x, y);
  const double form = linalg::Dot(space.GetComm(), x, y);
  CHECK(form > 0.0);
  CHECK_THAT(form, WithinRel(expected_form, 1.0e-9));
  Vector z(x.Size()), dKz(x.Size());
  z.Randomize(7);
  z.SetSubVector(space.GetNDDbcTDofLists().back(), 0.0);
  dK.Mult(z, dKz);
  CHECK_THAT(linalg::Dot(space.GetComm(), z, y),
             WithinRel(linalg::Dot(space.GetComm(), x, dKz), 1.0e-9));
  CHECK(linalg::Dot(space.GetComm(), z, dKz) >= -1.0e-12 * std::abs(form));
  // A gradient normal to every plane y = const (grad(y^2 / 2): zero on the sheet's PEC
  // dofs, so the zeroing leaves it curl-free): dK annihilates it.
  mfem::VectorFunctionCoefficient normal_gradient(3,
                                                  [](const mfem::Vector &p, mfem::Vector &E)
                                                  {
                                                    E.SetSize(3);
                                                    E = 0.0;
                                                    E[1] = p[1];
                                                  });
  GridFunction gradient(space.GetNDSpace(), true);
  gradient.Real().ProjectCoefficient(normal_gradient);
  Vector grad_true, dK_grad(x.Size());
  gradient.Real().GetTrueDofs(grad_true);
  dK.Mult(grad_true, dK_grad);
  CHECK(linalg::Norml2(space.GetComm(), dK_grad) <=
        1.0e-10 * std::abs(form) / std::max(linalg::Norml2(space.GetComm(), x), 1.0e-300));

  // The postprocessing: the magnetic domain correction 1/2 w a_m^H dS_m a_m / |omega|^2
  // equals 1/2 <x, dK x> at |omega| = 1, every closed loop of B . n vanishes (the dual loop
  // gate passes) and every patch carries the term.
  const auto result = response.GetMaxwellResponse(field, 1.0);
  CHECK(result.has_magnetic);
  CHECK(result.magnetic_loop_residual < 1.0e-10);
  CHECK(result.magnetic_weighted_loop_residual < 1.0e-10);
  CHECK(result.magnetic_loop_response_failure_fraction == 0.0);
  CHECK(result.magnetic_confident);
  CHECK(result.magnetic_patch_count == response.GetPatchCount());
  CHECK(result.magnetic_excluded_patch_count == 0);
  // No electric loop-fit fallback on a closed-loop field: nothing withheld (decision 1235).
  CHECK(result.magnetic_fallback_patch_count == 0);
  CHECK(result.loop_fallback_magnetic_domain_correction == 0.0);
  CHECK_THAT(result.magnetic_domain_correction, WithinRel(0.5 * form, 1.0e-9));
  // |omega| = 2 halves nothing of E but the trace is B = curl E / omega: a quarter.
  const auto result_2 = response.GetMaxwellResponse(field, 2.0);
  CHECK_THAT(result_2.magnetic_domain_correction, WithinRel(0.125 * form, 1.0e-9));
  REQUIRE(result.model_magnetic_domain_corrections.size() == 1);
  CHECK_THAT(result.model_magnetic_domain_corrections[0],
             WithinRel(result.magnetic_domain_correction, 1.0e-12));
#endif
}

TEST_CASE_METHOD(
    test::SurfaceResponseFiles,
    "SurfaceResponseOperator magnetic leaves: refusal rules r1-r4 and the switch",
    "[surfaceresponseoperator][3d][maxwell][magnetic][Serial][Parallel]")
{
#if !defined(MFEM_USE_GSLIB)
  SKIP("SurfaceResponseOperator requires MFEM_USE_GSLIB");
#else
  // One mesh read; every library below is loaded on it.
  IoData reader(MaxwellCpw3dConfig(*this, library_3d_path), false);
  auto device_mesh = mesh::ReadMesh(reader, Mpi::World());
  std::vector<std::unique_ptr<Mesh>> meshes;
  meshes.push_back(std::make_unique<Mesh>(std::move(device_mesh)));
  auto Load = [&](const fs::path &library, bool allow_fallback, bool magnetic = true)
  {
    json config = MaxwellCpw3dConfig(*this, library);
    if (allow_fallback)
    {
      config["Solver"]["SurfaceResponseCorrection"]["AllowElectricOnlyFallback"] = true;
    }
    if (!magnetic)
    {
      config["Solver"]["SurfaceResponseCorrection"]["Magnetic"] = false;
    }
    IoData iodata(config, false);
    ShareCrackedAttributes(iodata, reader);
    SpaceOperator space(iodata, meshes);
    return std::make_unique<SurfaceResponseOperator>(iodata, space);
  };
  auto Status = [](const SurfaceResponseOperator &response)
  {
    const auto records = response.GetMagneticModelRecords();
    REQUIRE(records.size() == 1);
    return records.front();
  };

  // r1: the twin's basis count differs from the model's.
  auto r1_source = DefaultMagneticSource();
  r1_source["BasisPointCount"] = 5;
  const auto r1 =
      WriteMagneticLibrary(library_3d_path, temp.temp_dir, "r1", kPsdDefect, r1_source);
  CHECK_THROWS_WITH(Load(r1.library, false),
                    Catch::Matchers::ContainsSubstring("REFUSED by rule r1"));
  {
    const auto response = Load(r1.library, true);
    CHECK_FALSE(response->HasMagneticModels());
    const auto status = Status(*response);
    CHECK(status.status == "Refused");
    CHECK(status.rule == "r1");
  }
  // r2: the twin's recorded process differs from the library's Fabrication block.
  auto r2_source = DefaultMagneticSource();
  r2_source["Fabrication"]["OveretchDepth"] = 0.05;
  const auto r2 =
      WriteMagneticLibrary(library_3d_path, temp.temp_dir, "r2", kPsdDefect, r2_source);
  CHECK_THROWS_WITH(Load(r2.library, false),
                    Catch::Matchers::ContainsSubstring("REFUSED by rule r2"));
  {
    const auto status = Status(*Load(r2.library, true));
    CHECK(status.status == "Refused");
    CHECK(status.rule == "r2");
    CHECK(status.reason.find("OveretchDepth") != std::string::npos);
  }
  // r3: the same-mesh proof exceeds the library's SameMeshTolerance (1e-12 of record).
  auto r3_source = DefaultMagneticSource();
  r3_source["ElectricRegenerationMaxRelativeDeviation"] = 8.1e-5;
  const auto r3 =
      WriteMagneticLibrary(library_3d_path, temp.temp_dir, "r3", kPsdDefect, r3_source);
  CHECK_THROWS_WITH(Load(r3.library, false),
                    Catch::Matchers::ContainsSubstring("REFUSED by rule r3"));
  {
    const auto status = Status(*Load(r3.library, true));
    CHECK(status.status == "Refused");
    CHECK(status.rule == "r3");
    CHECK_THAT(status.same_mesh_deviation, WithinAbs(-1.0, 0.0));
  }
  // A looser SameMeshTolerance in the library's MagneticRules admits the same leaf.
  const auto r3_loose = WriteMagneticLibrary(
      library_3d_path, temp.temp_dir, "r3-loose", kPsdDefect, r3_source,
      [](json &library) { library["MagneticRules"] = {{"SameMeshTolerance", 1.0e-4}}; });
  {
    const auto response = Load(r3_loose.library, false);
    CHECK(response->HasMagneticModels());
    const auto status = Status(*response);
    CHECK(status.status == "Applied");
    CHECK_THAT(status.same_mesh_deviation, WithinAbs(8.1e-5, 1.0e-20));
  }
  // r3, TWO-TIER (decision 1242): MagneticRules.SameMeshToleranceWithMeshIdentity applies
  // only to a twin whose MagneticSource.MeshOfRecordSha256Match proves BOTH meshes
  // identical to the model's meshes of record; absent (the default), SameMeshTolerance
  // applies everywhere.
  auto r3_rounding_source = DefaultMagneticSource();
  r3_rounding_source["ElectricRegenerationMaxRelativeDeviation"] = 1.0e-10;
  const auto two_tier = [](json &library)
  { library["MagneticRules"] = {{"SameMeshToleranceWithMeshIdentity", 1.0e-9}}; };
  // (a) 1e-10 WITHOUT the identity: refused by r3 at the 1e-12 tier.
  const auto r3_no_identity =
      WriteMagneticLibrary(library_3d_path, temp.temp_dir, "r3-two-tier-no-identity",
                           kPsdDefect, r3_rounding_source, two_tier);
  {
    const auto status = Status(*Load(r3_no_identity.library, true));
    CHECK(status.status == "Refused");
    CHECK(status.rule == "r3");
    CHECK(status.reason.find("mesh identity not proven") != std::string::npos);
  }
  // (b) 1e-10 WITH the identity: admitted at the 1e-9 tier.
  auto r3_identity_source = r3_rounding_source;
  r3_identity_source["MeshOfRecordSha256Match"] = {{"Thin", true}, {"Fabricated", true}};
  const auto r3_identity =
      WriteMagneticLibrary(library_3d_path, temp.temp_dir, "r3-two-tier-identity",
                           kPsdDefect, r3_identity_source, two_tier);
  {
    const auto response = Load(r3_identity.library, false);
    CHECK(response->HasMagneticModels());
    const auto status = Status(*response);
    CHECK(status.status == "Applied");
    CHECK_THAT(status.same_mesh_deviation, WithinAbs(1.0e-10, 1.0e-25));
  }
  // (c) 1e-8 WITH the identity: refused (above the 1e-9 tier).
  auto r3_identity_large = r3_identity_source;
  r3_identity_large["ElectricRegenerationMaxRelativeDeviation"] = 1.0e-8;
  const auto r3_identity_refused =
      WriteMagneticLibrary(library_3d_path, temp.temp_dir, "r3-two-tier-identity-refused",
                           kPsdDefect, r3_identity_large, two_tier);
  CHECK_THROWS_WITH(Load(r3_identity_refused.library, false),
                    Catch::Matchers::ContainsSubstring("REFUSED by rule r3"));
  {
    const auto status = Status(*Load(r3_identity_refused.library, true));
    CHECK(status.status == "Refused");
    CHECK(status.rule == "r3");
    CHECK(status.reason.find("SameMeshToleranceWithMeshIdentity") != std::string::npos);
  }
  // (d) the identity WITHOUT the leaf: SameMeshTolerance (1e-12) applies everywhere.
  const auto r3_identity_default =
      WriteMagneticLibrary(library_3d_path, temp.temp_dir, "r3-identity-default-rules",
                           kPsdDefect, r3_identity_source);
  {
    const auto status = Status(*Load(r3_identity_default.library, true));
    CHECK(status.status == "Refused");
    CHECK(status.rule == "r3");
  }
  // (e) a leaf below SameMeshTolerance is malformed.
  const auto r3_bad_rules = WriteMagneticLibrary(
      library_3d_path, temp.temp_dir, "r3-two-tier-bad", kPsdDefect, r3_identity_source,
      [](json &library)
      {
        library["MagneticRules"] = {{"SameMeshTolerance", 1.0e-9},
                                    {"SameMeshToleranceWithMeshIdentity", 1.0e-10}};
      });
  CHECK_THROWS_WITH(
      Load(r3_bad_rules.library, true),
      Catch::Matchers::ContainsSubstring("SameMeshToleranceWithMeshIdentity"));
  // r4, two-level (ERRATUM 1227): a defect with one negative eigenvalue; its ratio
  // -lambda_min / ||dS_m|| against PsdFindingRatio (1e-2) and PsdRefusalRatio (0.1).
  auto finding_defect = kPsdDefect;
  finding_defect[3][3] = -0.1e-12;  // ~ -0.05 of the largest eigenvalue: a FINDING
  const auto r4_finding =
      WriteMagneticLibrary(library_3d_path, temp.temp_dir, "r4-finding", finding_defect);
  {
    const auto response = Load(r4_finding.library, false);
    CHECK(response->HasMagneticModels());
    const auto status = Status(*response);
    CHECK(status.status == "Finding");
    CHECK(status.rule == "r4");
    CHECK(status.psd_ratio > 1.0e-2);
    CHECK(status.psd_ratio <= 0.1);
    CHECK(status.negative_eigenvalue_count == 1);
  }
  auto refusal_defect = kPsdDefect;
  refusal_defect[3][3] = -1.5e-12;  // ~ -0.6 of the largest eigenvalue: a REFUSAL
  const auto r4_refusal =
      WriteMagneticLibrary(library_3d_path, temp.temp_dir, "r4-refusal", refusal_defect);
  {
    // Never fatal: electric-only by name, without the fallback switch.
    const auto response = Load(r4_refusal.library, false);
    CHECK_FALSE(response->HasMagneticModels());
    const auto status = Status(*response);
    CHECK(status.status == "Refused");
    CHECK(status.rule == "r4");
    CHECK(status.psd_ratio > 0.1);
  }
  // The library-level PsdRefusalRatio moves the line (a library parameter, not a constant).
  const auto r4_moved = WriteMagneticLibrary(
      library_3d_path, temp.temp_dir, "r4-moved", refusal_defect, DefaultMagneticSource(),
      [](json &library) { library["MagneticRules"] = {{"PsdRefusalRatio", 0.9}}; });
  {
    const auto response = Load(r4_moved.library, false);
    CHECK(response->HasMagneticModels());
    CHECK(Status(*response).status == "Finding");
    // Decision 1251: without MagneticRules.ConductorConsistencyTolerance the dual gate
    // reads the electric gate's constant (bitwise as before).
    CHECK_THAT(response->MagneticConductorConsistencyTolerance(),
               WithinAbs(SurfaceResponseOperator::kConductorConsistencyTolerance, 0.0));
  }
  // Decision 1251 (INTERIM): the dual gate's tolerance as a library parameter; a negative
  // value in the file is malformed.
  const auto gate_moved = WriteMagneticLibrary(
      library_3d_path, temp.temp_dir, "gate-moved", finding_defect, DefaultMagneticSource(),
      [](json &library)
      { library["MagneticRules"] = {{"ConductorConsistencyTolerance", 0.25}}; });
  {
    const auto response = Load(gate_moved.library, false);
    CHECK(response->HasMagneticModels());
    CHECK(Status(*response).status == "Finding");
    CHECK_THAT(response->MagneticConductorConsistencyTolerance(), WithinAbs(0.25, 0.0));
  }
  const auto gate_bad = WriteMagneticLibrary(
      library_3d_path, temp.temp_dir, "gate-bad", finding_defect, DefaultMagneticSource(),
      [](json &library)
      { library["MagneticRules"] = {{"ConductorConsistencyTolerance", -0.5}}; });
  CHECK_THROWS_WITH(Load(gate_bad.library, true),
                    Catch::Matchers::ContainsSubstring("ConductorConsistencyTolerance"));
  // The switch (decision 1229 Q1 (i)): Magnetic = false reads no leaf at all.
  {
    const auto response = Load(r4_finding.library, false, false);
    CHECK_FALSE(response->HasMagneticModels());
    CHECK(Status(*response).status == "Absent");
  }
  // A leaf without its MagneticSource is malformed at the catalog.
  const auto incomplete = WriteMagneticLibrary(
      library_3d_path, temp.temp_dir, "incomplete", kPsdDefect, DefaultMagneticSource(),
      [](json &library) { library["Models"][0].erase("MagneticSource"); });
  CHECK_THROWS_WITH(Load(incomplete.library, true),
                    Catch::Matchers::ContainsSubstring("MagneticSource together"));
#endif
}

// D6 (USER 1221; decision 1229), the corrected frequency of record, against the Phase-0
// measurement (measure/REPORT.md RESULTS; prereg-gates.md G-d FORM; review d1227 MAJOR-2):
// the fixed-trace shift (U_m,defect - U_e,defect) / (2 U_raw) on the measured within-R
// defects. In the quasi-TEM reading U_e,defect / U = dC'/C' (N1) at fixed voltage and
// U_m,defect / U = dC'_0/C'_0 (N2) at fixed flux, so the shift is -1/2 (N1 - N2) = the
// REPORT's MAGNETIC (TRUE) df/f: +3.257 % (3 / 3 / 3, N1 -2.434 %, N2 +4.080 %) and
// +4.584 % (2 / 2 / 2, N1 -3.378 %, N2 +5.790 %), to the REPORT's printed 1e-3 pt and the
// gate's 1e-4 pt on the derived value; the gates withhold it (NaN), never the trace-closure
// spread; the electric half alone gives the REPORT's ELECTRIC-ONLY +1.217 % / +1.689 %.
// Phase 3 (a) t1 / t3 / t4 / t5 (DESIGN-v2 a.5; AMENDMENT 1275 MAJOR-A): the LIFTED
// magnetic probe legs on the shorted CPW stub with a magnetic leaf.
TEST_CASE_METHOD(
    test::SurfaceResponseFiles,
    "SurfaceResponseOperator magnetic probe form: LIFTED legs, FLAT, BAND and "
    "the setbacks (Phase 3 (a))",
    "[surfaceresponseoperator][maxwell][magnetic][probeform][Serial][Parallel]")
{
#if !defined(MFEM_USE_GSLIB)
  SKIP("SurfaceResponseOperator requires MFEM_USE_GSLIB");
#else
  auto Build = [&](const json &run_config, std::vector<std::unique_ptr<Mesh>> &meshes)
  {
    IoData iodata(run_config, false);
    iodata.boundaries.cracked_attributes.insert(9);
    mfem::Mesh serial = test::MakeSerialShortedCpwMesh();
    meshes.push_back(
        std::make_unique<Mesh>(std::make_unique<mfem::ParMesh>(Mpi::World(), serial)));
    return iodata;
  };
  constexpr double R = 0.2;
  const std::array<double, 3> B = {0.3, -0.5, 0.8};
  const double B_norm = std::sqrt(B[0] * B[0] + B[1] * B[1] + B[2] * B[2]);
  // One operator per probe form on the SAME mesh / config (the library differs by its
  // MagneticRules only): FLAT (0 / 0), LIFTED without setback (0.5 / 0), LIFTED (the
  // defaults 0.5 / 0.5), LIFTED with a setback above the strip's width cap (0.5 / 1.0) and
  // BAND (3 Gauss heights, no setback).
  struct Form
  {
    std::string name;
    json rules;
    std::vector<std::unique_ptr<Mesh>> meshes;
    std::unique_ptr<IoData> iodata;
    std::unique_ptr<SpaceOperator> space;
    std::unique_ptr<SurfaceResponseOperator> response;
    std::vector<std::vector<std::complex<double>>> uniform_readings;
  };
  std::vector<Form> forms;
  for (const auto &[name, rules] : std::vector<std::pair<std::string, json>>{
           {"flat", {{"ProbeLiftOverR", 0.0}, {"ProbeSetbackOverR", 0.0}}},
           {"lift", {{"ProbeLiftOverR", 0.5}, {"ProbeSetbackOverR", 0.0}}},
           {"lifted", json::object()},
           {"capped", {{"ProbeLiftOverR", 0.5}, {"ProbeSetbackOverR", 1.0}}},
           {"band",
            {{"ProbeLiftOverR", 0.5},
             {"ProbeSetbackOverR", 0.0},
             {"ProbeBandHeights", 3}}}})
  {
    Form form;
    form.name = name;
    form.rules = rules;
    const auto library =
        WriteStubMagneticLibrary(temp.temp_dir, convex_library_3d_path, name, false, rules);
    form.iodata =
        std::make_unique<IoData>(Build(StubMaxwellConfig(*this, library), form.meshes));
    form.space = std::make_unique<SpaceOperator>(*form.iodata, form.meshes);
    form.response = std::make_unique<SurfaceResponseOperator>(*form.iodata, *form.space);
    REQUIRE(form.response->GetPatchCount() > 0);
    REQUIRE(form.response->HasMagneticModels());
    REQUIRE(form.response->HasConductorConsistencyProbes());
    auto uniform = ConstantCurlField(B);
    form.uniform_readings =
        form.response->EvaluateProbeReadings(ProjectField(*form.space, uniform), true);
    forms.push_back(std::move(form));
  }
  auto &flat = forms[0];
  auto &lift = forms[1];
  auto &lifted = forms[2];
  auto &capped = forms[3];
  auto &band = forms[4];
  CHECK(flat.response->MagneticProbeForm() == "FLAT");
  CHECK(lift.response->MagneticProbeForm() == "LIFTED");
  CHECK(lifted.response->MagneticProbeForm() == "LIFTED");
  CHECK(band.response->MagneticProbeForm() == "BAND");
  const auto &flat_probes = flat.response->GetMaxwellProbes();
  const std::size_t local_patches = flat_probes.size();
  for (const auto &form : forms)
  {
    REQUIRE(form.response->GetMaxwellProbes().size() == local_patches);
    REQUIRE(form.uniform_readings.size() == local_patches);
  }
  const double tolerance = 1.0e-10 * B_norm * R;
  int compared = 0;
  for (std::size_t patch = 0; patch < local_patches; patch++)
  {
    REQUIRE(flat_probes[patch].size() == 8);  // 4 edges x 2 end faces
    for (std::size_t q = 0; q < flat_probes[patch].size(); q++)
    {
      const auto &flat_probe = flat_probes[patch][q];
      // t3: FLAT's magnetic legs ARE the electric legs (the a92686289e gate verbatim).
      CHECK(flat_probe.magnetic_legs == flat_probe.legs);
      CHECK(flat_probe.magnetic_path.empty());
      // t1: path independence on the uniform B (closed loop flat -> lifted = zero flux):
      // the lifted reading without setback equals the flat reading; BAND (every height
      // reads the same value, the weights sum to 1) equals it too.
      const auto flat_reading = flat.uniform_readings[patch][q];
      CHECK(std::abs(lift.uniform_readings[patch][q] - flat_reading) < tolerance);
      CHECK(std::abs(band.uniform_readings[patch][q] - flat_reading) < tolerance);
      // The lifted paths: anchor* -> A -> B -> probe* (4 points for one height; the BAND
      // path climbs through its 3 heights on both sides: 8 points).
      const auto &lift_probe = lift.response->GetMaxwellProbes()[patch][q];
      const auto &band_probe = band.response->GetMaxwellProbes()[patch][q];
      REQUIRE(lift_probe.magnetic_path.size() == 4);
      REQUIRE(band_probe.magnetic_path.size() == 8);
      // The reading of a lifted path on the uniform B is -(probe* - anchor*) . (w x B)
      // with w = u x n (the patch frame); u and n read from the path itself (the vertical
      // leg gives n, the horizontal leg u) wherever the horizontal leg exists.
      auto Expected = [&](const auto &probe)
      {
        const auto &path = probe.magnetic_path;
        const auto &a0 = path.front();
        const auto &a1 = path[1];
        const auto &b1 = path[path.size() - 2];
        std::array<double, 3> n{}, u{}, delta{};
        double n_norm = 0.0, u_norm = 0.0;
        for (int d = 0; d < 3; d++)
        {
          n[d] = a1[d] - a0[d];
          u[d] = b1[d] - a1[d];
          delta[d] = path.back()[d] - a0[d];
          n_norm += n[d] * n[d];
          u_norm += u[d] * u[d];
        }
        n_norm = std::sqrt(n_norm);
        u_norm = std::sqrt(u_norm);
        if (u_norm < 1.0e-12)
        {
          return std::optional<double>();  // anchor* and probe* coincide: a zero path
        }
        for (int d = 0; d < 3; d++)
        {
          n[d] /= n_norm;
          u[d] /= u_norm;
        }
        const std::array<double, 3> w = {u[1] * n[2] - u[2] * n[1],
                                         u[2] * n[0] - u[0] * n[2],
                                         u[0] * n[1] - u[1] * n[0]};
        const std::array<double, 3> wxB = {w[1] * B[2] - w[2] * B[1],
                                           w[2] * B[0] - w[0] * B[2],
                                           w[0] * B[1] - w[1] * B[0]};
        return std::optional<double>(
            -(delta[0] * wxB[0] + delta[1] * wxB[1] + delta[2] * wxB[2]));
      };
      for (const auto *form : {&lift, &lifted, &capped})
      {
        const auto &probe = form->response->GetMaxwellProbes()[patch][q];
        const auto expected = Expected(probe);
        const auto reading = form->uniform_readings[patch][q];
        CHECK(std::abs(reading.imag()) < tolerance);
        if (expected)
        {
          INFO(form->name << " patch " << patch << " probe " << q);
          CHECK_THAT(reading.real(), WithinAbs(*expected, tolerance));
          compared++;
        }
        else
        {
          CHECK(std::abs(reading) < tolerance);
        }
      }
      // t4 / t5 geometry (the defaults 0.5 / 0.5 / 0.45 / 1e-6): every path point sits off
      // the sheet (y = 0.5) by eps = 1e-6 R at its ends and by h_m = 0.5 R on the
      // horizontal leg; probe* is set back s = 0.5 R = 0.1 from its edge point on every
      // edge (the strip's cap 0.45 x 0.25 = 0.1125 does not bind); under 1.0 R the cap
      // binds on the strip (0.1125) and not on the grounds (0.2).
      const auto &lifted_probe = lifted.response->GetMaxwellProbes()[patch][q];
      const auto &capped_probe = capped.response->GetMaxwellProbes()[patch][q];
      const double eps = 1.0e-6 * R, h_m = 0.5 * R;
      for (const auto *probe : {&lifted_probe, &capped_probe})
      {
        const auto &path = probe->magnetic_path;
        CHECK_THAT(std::abs(path[0][1] - 0.5), WithinAbs(eps, 1.0e-12));
        CHECK_THAT(std::abs(path[1][1] - 0.5), WithinAbs(eps + h_m, 1.0e-12));
        CHECK_THAT(std::abs(path[2][1] - 0.5), WithinAbs(eps + h_m, 1.0e-12));
        CHECK_THAT(std::abs(path[3][1] - 0.5), WithinAbs(eps, 1.0e-12));
        // The vertical legs are vertical (x and z unchanged).
        CHECK_THAT(path[1][0], WithinAbs(path[0][0], 1.0e-12));
        CHECK_THAT(path[1][2], WithinAbs(path[0][2], 1.0e-12));
        CHECK_THAT(path[2][0], WithinAbs(path[3][0], 1.0e-12));
        CHECK_THAT(path[2][2], WithinAbs(path[3][2], 1.0e-12));
      }
      auto LateralSetback = [](const auto &probe)
      {
        // probe* - the probe point along the lateral axis z (the stub's cross-section).
        return std::abs(probe.magnetic_path.back()[2] - probe.point[2]);
      };
      CHECK_THAT(LateralSetback(lifted_probe), WithinAbs(0.1, 1.0e-12));
      const bool strip_edge = lifted_probe.edge == 1 || lifted_probe.edge == 2;
      CHECK_THAT(LateralSetback(capped_probe),
                 WithinAbs(strip_edge ? 0.1125 : 0.2, 1.0e-12));
      // t5 (MAJOR-1): the anchor is the ground's edge point (offset 0, piece 1); anchor* is
      // set back INTO the ground (toward -u, the unbounded side: no cap): the edge-0
      // probe's probe* coincides laterally with anchor* (both 0.1 into the ground).
      if (lifted_probe.edge == 0)
      {
        CHECK_THAT(lifted_probe.magnetic_path.back()[2],
                   WithinAbs(lifted_probe.magnetic_path.front()[2], 1.0e-12));
        CHECK_THAT(capped_probe.magnetic_path.back()[2],
                   WithinAbs(capped_probe.magnetic_path.front()[2], 1.0e-12));
      }
      else
      {
        // anchor* lies 0.1 (0.2) outside the ground's edge, i.e. 0.1 (0.2) further from the
        // strip than the edge-0 point: |anchor* - probe*_edge0| = 0 and the anchor* of this
        // probe equals that of the edge-0 probe of the same end face.
        for (const auto &other : lifted.response->GetMaxwellProbes()[patch])
        {
          if (other.edge == 0 && other.end_face == lifted_probe.end_face)
          {
            CHECK_THAT(lifted_probe.magnetic_path.front()[2],
                       WithinAbs(other.magnetic_path.front()[2], 1.0e-12));
            CHECK_THAT(std::abs(other.magnetic_path.back()[2] - other.point[2]),
                       WithinAbs(0.1, 1.0e-12));
          }
        }
      }
    }
  }
  int global_compared = compared;
  Mpi::GlobalSum(1, &global_compared, Mpi::World());
  CHECK(global_compared > 0);
  // The statistics name the form and count the magnetic legs (FLAT adds none).
  CHECK(flat.response->GetStatistics()["Correction"]["MagneticProbeLines"]["Total"]
            .get<long long int>() == 0);
  CHECK(lifted.response->GetStatistics()["Correction"]["MagneticProbeLines"]["Total"]
            .get<long long int>() > 0);
  CHECK(lifted.response->GetStatistics()["Correction"]["MagneticProbeForm"] == "LIFTED");
#endif
}

// Phase 3 (a) t2 = the a.3 DISCRIMINATION CONTROL (decisions 1256 / 1292 (B)): on the
// piece-form stub (the gate compares a probe with its PIECE's flux state) under the
// gap-flux field, the correctly mapped probes read ratio ~0 at FLAT and LIFTED, and
// ConductorConsistency.ControlFlipProbePiece (RecordOnly) relabelling the strip's right
// edge (edge 2) as the right ground (piece 3) makes every patch read ratio ~1 (the strip's
// flux V against the ground's ~0) at both forms: the gate catches a one-label mis-mapping
// by >= 2x its tolerance. The control is refused without RecordOnly and when it applies
// nowhere.
TEST_CASE_METHOD(test::SurfaceResponseFiles,
                 "SurfaceResponseOperator magnetic gate: the flipped-label control "
                 "(ControlFlipProbePiece, Phase 3 (a))",
                 "[surfaceresponseoperator][maxwell][magnetic][probeform][control][Serial]"
                 "[Parallel]")
{
#if !defined(MFEM_USE_GSLIB)
  SKIP("SurfaceResponseOperator requires MFEM_USE_GSLIB");
#else
  auto Build = [&](const json &run_config, std::vector<std::unique_ptr<Mesh>> &meshes)
  {
    IoData iodata(run_config, false);
    iodata.boundaries.cracked_attributes.insert(9);
    mfem::Mesh serial = test::MakeSerialShortedCpwMesh();
    meshes.push_back(
        std::make_unique<Mesh>(std::make_unique<mfem::ParMesh>(Mpi::World(), serial)));
    return iodata;
  };
  const double tolerance = SurfaceResponseOperator::kConductorConsistencyTolerance;
  const json flat_rules = {{"ProbeLiftOverR", 0.0}, {"ProbeSetbackOverR", 0.0}};
  for (const auto &[form, rules] :
       std::vector<std::pair<std::string, json>>{{"FLAT", flat_rules}, {"LIFTED", {}}})
  {
    const auto library = WriteStubMagneticLibrary(temp.temp_dir, convex_library_3d_path,
                                                  "pieces-" + form, true, rules);
    json config = StubMaxwellConfig(*this, library);
    // The correctly mapped piece form: every probe within the tolerance, nothing excluded.
    {
      std::vector<std::unique_ptr<Mesh>> meshes;
      IoData iodata = Build(config, meshes);
      SpaceOperator space(iodata, meshes);
      SurfaceResponseOperator response(iodata, space);
      REQUIRE(response.GetPatchCount() > 0);
      REQUIRE(response.HasMagneticModels());
      CHECK(response.MagneticProbeForm() == form);
      CHECK(response.GetControlFlippedProbeRecords().empty());
      auto field = GapFluxField(1.0);
      response.ApplyMagneticConductorConsistencyGate(ProjectField(space, field), 0, 5.0);
      const auto &records = response.GetMagneticGateRecords();
      REQUIRE(records.size() == static_cast<std::size_t>(response.GetPatchCount()));
      for (const auto &record : records)
      {
        INFO(form << " patch " << record.patch << " ratio " << record.max_ratio);
        CHECK(record.max_ratio < 0.5 * tolerance);
        CHECK_FALSE(record.excluded);
      }
    }
    // The CONTROL: edge 2 (the strip's right edge, piece 2) relabelled piece 3.
    {
      json control = config;
      control["Solver"]["SurfaceResponseCorrection"]["ConductorConsistency"] = {
          {"RecordOnly", true},
          {"ControlFlipProbePiece",
           {{{"Model", "cpw-one-label-four-edge-pieces"}, {"Edge", 2}, {"Piece", 3}}}}};
      std::vector<std::unique_ptr<Mesh>> meshes;
      IoData iodata = Build(control, meshes);
      SpaceOperator space(iodata, meshes);
      SurfaceResponseOperator response(iodata, space);
      // One flip per end face of every local patch, recorded by name.
      const auto &flips = response.GetControlFlippedProbeRecords();
      CHECK(flips.size() == 2 * response.GetMaxwellProbes().size());
      for (const auto &flip : flips)
      {
        CHECK(flip == "cpw-one-label-four-edge-pieces:2:2->3");
      }
      const int flipped_pieces[4] = {1, 2, 3, 3};  // edges 0..3; edge 2's own piece is 2
      for (const auto &probes : response.GetMaxwellProbes())
      {
        for (const auto &probe : probes)
        {
          CHECK(probe.piece == flipped_pieces[probe.edge]);
        }
      }
      auto field = GapFluxField(1.0);
      response.ApplyMagneticConductorConsistencyGate(ProjectField(space, field), 0, 5.0);
      const auto &records = response.GetMagneticGateRecords();
      REQUIRE(records.size() == static_cast<std::size_t>(response.GetPatchCount()));
      for (const auto &record : records)
      {
        INFO(form << " patch " << record.patch << " ratio " << record.max_ratio);
        // The flipped probe reads the strip's flux V against the right ground's state (~0):
        // ratio ~1 = V / the trace amplitude (the strip state), >= 2x the tolerance.
        CHECK(record.max_ratio > 2.0 * tolerance);
        CHECK_THAT(record.max_ratio, WithinAbs(1.0, 0.2));
        CHECK(record.excluded);  // the verdict is recorded; RecordOnly excludes nothing
      }
      const auto statistics = response.GetStatistics();
      int global_flips = static_cast<int>(flips.size());
      Mpi::GlobalSum(1, &global_flips, Mpi::World());
      CHECK(statistics["Correction"]["ControlFlipProbePiece"]["Count"].get<int>() ==
            global_flips);
      CHECK(statistics["Correction"]["ControlFlipProbePiece"]["Forms"].size() == 1);
      // The electric gate's twin reading (information): the flipped probe departs there
      // too.
      const auto electric =
          response.ApplyConductorConsistencyGate(ProjectField(space, field), 0, 5.0, false);
      CHECK_FALSE(electric.empty());
    }
    // REFUSED: the control without RecordOnly; a Model / Edge that flips nothing.
    {
      json bad = config;
      bad["Solver"]["SurfaceResponseCorrection"]["ConductorConsistency"] = {
          {"ControlFlipProbePiece",
           {{{"Model", "cpw-one-label-four-edge-pieces"}, {"Edge", 2}, {"Piece", 3}}}}};
      CHECK_THROWS_WITH(IoData(bad, false),
                        Catch::Matchers::ContainsSubstring("RecordOnly"));
      json dangling = config;
      dangling["Solver"]["SurfaceResponseCorrection"]["ConductorConsistency"] = {
          {"RecordOnly", true},
          {"ControlFlipProbePiece",
           {{{"Model", "cpw-one-label-four-edge-pieces"}, {"Edge", 7}, {"Piece", 3}}}}};
      std::vector<std::unique_ptr<Mesh>> meshes;
      IoData iodata = Build(dangling, meshes);
      SpaceOperator space(iodata, meshes);
      CHECK_THROWS_WITH(SurfaceResponseOperator(iodata, space),
                        Catch::Matchers::ContainsSubstring("flipped no probe"));
    }
  }
#endif
}

// The hexahedral stub's blind line quadrature at the slot edge (the H9a limitation the R0-a
// tests record): a twin anchor line that crosses the edge where f'(z) jumps reads ~2 % of
// the gap flux (0.0197 on the chain, 0.0202 with the magnetic loop fit); the bound the
// re-anchored gate readings are held to on this fixture (the (ii) tetrahedral meshes split
// the coverage exactly; the device readings are the measurement).
constexpr double kHexSlotEdgeQuadratureBound = 0.025;

// Decision 1302: the ControlStateAnchorSetbackOverR DIAGNOSTIC CONTROL re-anchors the
// MAGNETIC chain's anchor / end / return lines at set-back conductor anchors (RecordOnly +
// PostprocessOnly only). On the piece-form stub under the gap-flux field the re-anchored
// states equal the edge-anchored ones in the continuum (path independence is exact for this
// field); on the HEXAHEDRAL fixture the twin anchor lines cross the slot edge where the
// integrand f'(z) jumps and the line quadrature samples a hexahedral mesh blindly (no face
// split: the exact segment intersections exist on simplices only, the H9a site), so the
// twins read the jump with a ~2 % quadrature error (0.0197 on the chain; 0.0202 once the
// magnetic loop fit, R-MAG-17, redistributes the twin loop's misclosure): every gate ratio
// stays at that fixture bound and uniform (the same re-anchored geometry on every cell).
// The twin lines are counted as contour lines; the control is recorded; it is refused
// without RecordOnly, on a self-consistent run and out of range.
TEST_CASE_METHOD(test::SurfaceResponseFiles,
                 "SurfaceResponseOperator magnetic gate: the state-anchor setback control "
                 "(ControlStateAnchorSetbackOverR, decision 1302)",
                 "[surfaceresponseoperator][maxwell][magnetic][probeform][control][Serial]"
                 "[Parallel]")
{
#if !defined(MFEM_USE_GSLIB)
  SKIP("SurfaceResponseOperator requires MFEM_USE_GSLIB");
#else
  auto Build = [&](const json &run_config, std::vector<std::unique_ptr<Mesh>> &meshes)
  {
    IoData iodata(run_config, false);
    iodata.boundaries.cracked_attributes.insert(9);
    mfem::Mesh serial = test::MakeSerialShortedCpwMesh();
    meshes.push_back(
        std::make_unique<Mesh>(std::make_unique<mfem::ParMesh>(Mpi::World(), serial)));
    return iodata;
  };
  const auto library = WriteStubMagneticLibrary(temp.temp_dir, convex_library_3d_path,
                                                "pieces-stateanchor", true);
  json config = StubMaxwellConfig(*this, library);
  config["Solver"]["SurfaceResponseCorrection"]["CorrectionMode"] = "PostprocessOnly";
  long long electric_lines = 0;
  {
    std::vector<std::unique_ptr<Mesh>> meshes;
    IoData iodata = Build(config, meshes);
    SpaceOperator space(iodata, meshes);
    SurfaceResponseOperator response(iodata, space);
    electric_lines =
        response.GetStatistics()["Correction"]["ContourLines"]["Total"].get<long long>();
    CHECK_FALSE(
        response.GetStatistics()["Correction"].contains("ControlStateAnchorSetbackOverR"));
  }
  {
    json control = config;
    control["Solver"]["SurfaceResponseCorrection"]["ConductorConsistency"] = {
        {"RecordOnly", true}, {"ControlStateAnchorSetbackOverR", 0.5}};
    std::vector<std::unique_ptr<Mesh>> meshes;
    IoData iodata = Build(control, meshes);
    SpaceOperator space(iodata, meshes);
    SurfaceResponseOperator response(iodata, space);
    const auto statistics = response.GetStatistics();
    REQUIRE(statistics["Correction"].contains("ControlStateAnchorSetbackOverR"));
    CHECK(statistics["Correction"]["ControlStateAnchorSetbackOverR"]["Patches"]["Total"]
              .get<long long>() == response.GetPatchCount());
    // The twin anchor / end / return lines are contour lines beside the electric ones.
    CHECK(statistics["Correction"]["ContourLines"]["Total"].get<long long>() >
          electric_lines);
    auto field = GapFluxField(1.0);
    response.ApplyMagneticConductorConsistencyGate(ProjectField(space, field), 0, 5.0);
    const auto &records = response.GetMagneticGateRecords();
    REQUIRE(records.size() == static_cast<std::size_t>(response.GetPatchCount()));
    double lowest = mfem::infinity(), highest = 0.0;
    for (const auto &record : records)
    {
      INFO("patch " << record.patch << " ratio " << record.max_ratio);
      CHECK(record.max_ratio < kHexSlotEdgeQuadratureBound);
      lowest = std::min(lowest, record.max_ratio);
      highest = std::max(highest, record.max_ratio);
    }
    CHECK(highest - lowest < 1.0e-6);  // the same re-anchored geometry on every cell
  }
  {
    json bad = config;
    bad["Solver"]["SurfaceResponseCorrection"]["ConductorConsistency"] = {
        {"ControlStateAnchorSetbackOverR", 0.5}};
    CHECK_THROWS_WITH(IoData(bad, false), Catch::Matchers::ContainsSubstring("RecordOnly"));
    json self_consistent = config;
    self_consistent["Solver"]["SurfaceResponseCorrection"]["CorrectionMode"] = "Both";
    self_consistent["Solver"]["SurfaceResponseCorrection"]["ConductorConsistency"] = {
        {"RecordOnly", true}, {"ControlStateAnchorSetbackOverR", 0.5}};
    CHECK_THROWS_WITH(IoData(self_consistent, false),
                      Catch::Matchers::ContainsSubstring("PostprocessOnly"));
    json out_of_range = config;
    out_of_range["Solver"]["SurfaceResponseCorrection"]["ConductorConsistency"] = {
        {"RecordOnly", true}, {"ControlStateAnchorSetbackOverR", 1.5}};
    CHECK_THROWS_WITH(IoData(out_of_range, false),
                      Catch::Matchers::ContainsSubstring("(0, 1]"));
  }
#endif
}

// R-MAG-17 (USER 1308 U-2, decision 1309): the library leaf
// MagneticRules.StateAnchorSetbackOverR makes the set-back state anchors the magnetic
// TERM's path: the magnetic chain, its own loop fit, the dual gate and dK read the twin
// lines. On the piece-form stub (R 0.2, s_a = 0.1): every conductor anchor is an EDGE point
// of its piece (the ground's u = 0, the strip's 0.125, the right ground's 0.5: the six
// anchor forms of AMENDMENT 1275 MAJOR-A), so each moves s_a = 0.1 into its piece (-0.1 u
// for the ground, +0.1 u for the strip and the right ground). Under the constant-curl field
// the flux trace is path-independent, so every ring knot's state shifts by the one constant
// c0 = -(w x B) . (anchor_0 - anchor_0*) = -0.1 (w x B) . u and both conductor states by 2
// c0 (their own anchor's move adds another c0); the postprocessed correction equals 1/2 <x,
// dK x> (the transpose on the twins is the adjoint) and dK stays symmetric. The gate under
// the leaf reads the gap-flux field as the 1302 control did (the hexahedral fixture's
// bound). Absent = 0 = bitwise with the legacy form; the control is refused by name on a
// library whose leaf already re-anchors.
TEST_CASE_METHOD(test::SurfaceResponseFiles,
                 "SurfaceResponseOperator magnetic term: the state-anchor setback leaf "
                 "(MagneticRules.StateAnchorSetbackOverR, R-MAG-17)",
                 "[surfaceresponseoperator][maxwell][magnetic][probeform][stateanchor]"
                 "[Serial][Parallel]")
{
#if !defined(MFEM_USE_GSLIB)
  SKIP("SurfaceResponseOperator requires MFEM_USE_GSLIB");
#else
  auto Build = [&](const json &run_config, std::vector<std::unique_ptr<Mesh>> &meshes)
  {
    IoData iodata(run_config, false);
    iodata.boundaries.cracked_attributes.insert(9);
    mfem::Mesh serial = test::MakeSerialShortedCpwMesh();
    meshes.push_back(
        std::make_unique<Mesh>(std::make_unique<mfem::ParMesh>(Mpi::World(), serial)));
    return iodata;
  };
  // No sheet nudge: the set-back anchors are exactly anchor - s_a u (the shifts below are
  // exact); the gate check runs with the nudge of record.
  const json leaf_rules = {{"StateAnchorSetbackOverR", 0.5}, {"ProbeSheetNudgeOverR", 0.0}};
  const json legacy_rules = {{"ProbeSheetNudgeOverR", 0.0}};
  const auto leaf_library = WriteStubMagneticLibrary(
      temp.temp_dir, convex_library_3d_path, "pieces-stateanchor-leaf", true, leaf_rules);
  const auto legacy_library =
      WriteStubMagneticLibrary(temp.temp_dir, convex_library_3d_path,
                               "pieces-stateanchor-legacy", true, legacy_rules);
  const std::array<double, 3> B = {0.3, -0.5, 0.8};
  struct Operator
  {
    std::vector<std::unique_ptr<Mesh>> meshes;
    std::unique_ptr<IoData> iodata;
    std::unique_ptr<SpaceOperator> space;
    std::unique_ptr<SurfaceResponseOperator> response;
  };
  auto Make = [&](const fs::path &library)
  {
    Operator op;
    op.iodata =
        std::make_unique<IoData>(Build(StubMaxwellConfig(*this, library), op.meshes));
    op.space = std::make_unique<SpaceOperator>(*op.iodata, op.meshes);
    op.response = std::make_unique<SurfaceResponseOperator>(*op.iodata, *op.space);
    return op;
  };
  auto legacy = Make(legacy_library);
  auto leaf = Make(leaf_library);
  REQUIRE(leaf.response->GetPatchCount() == legacy.response->GetPatchCount());
  REQUIRE(leaf.response->GetPatchCount() > 0);
  {
    const auto statistics = leaf.response->GetStatistics();
    const auto legacy_statistics = legacy.response->GetStatistics();
    CHECK(statistics["Correction"]["MagneticStateAnchorSetbackOverR"].get<double>() == 0.5);
    CHECK(
        statistics["Correction"]["MagneticStateAnchorPatches"]["Total"].get<long long>() ==
        leaf.response->GetPatchCount());
    CHECK_FALSE(statistics["Correction"].contains("ControlStateAnchorSetbackOverR"));
    CHECK(
        legacy_statistics["Correction"]["MagneticStateAnchorSetbackOverR"].get<double>() ==
        0.0);
    CHECK(legacy_statistics["Correction"]["MagneticStateAnchorPatches"]["Total"]
              .get<long long>() == 0);
    // The twin anchor / end / return lines are contour lines beside the electric ones.
    CHECK(statistics["Correction"]["ContourLines"]["Total"].get<long long>() >
          legacy_statistics["Correction"]["ContourLines"]["Total"].get<long long>());
    CHECK(leaf.response->GetMagneticProbeRules().state_anchor_setback_over_r == 0.5);
  }
  // The constant-curl field: the exact shifts of the states.
  auto coefficient = ConstantCurlField(B);
  GridFunction field(leaf.space->GetNDSpace(), true);
  field.Real().ProjectCoefficient(coefficient);
  field.Imag() = 0.0;
  Vector x_exact, x;
  field.Real().GetTrueDofs(x_exact);
  x = x_exact;
  x.SetSubVector(leaf.space->GetNDDbcTDofLists().back(), 0.0);
  field.Real().SetFromTrueDofs(x);
  const auto leaf_traces = leaf.response->GetMagneticPatchTraces(x_exact);
  const auto legacy_traces = legacy.response->GetMagneticPatchTraces(x_exact);
  const auto &assignments = leaf.response->GetPatchAssignments();
  if (Mpi::Root(Mpi::World()))
  {
    REQUIRE(leaf_traces.size() == legacy_traces.size());
    REQUIRE(leaf_traces.size() == assignments.size());
    for (std::size_t t = 0; t < leaf_traces.size(); t++)
    {
      const auto &shifted = leaf_traces[t];
      const auto &reference = legacy_traces[t];
      REQUIRE(shifted.patch == reference.patch);
      REQUIRE(shifted.coefficients.size() == 6);
      REQUIRE(reference.coefficients.size() == 6);
      const auto assignment =
          std::find_if(assignments.begin(), assignments.end(),
                       [&](const auto &a) { return a.global_index == shifted.patch; });
      REQUIRE(assignment != assignments.end());
      const auto &u = assignment->axis_u;
      const auto &v = assignment->axis_v;
      const std::array<double, 3> w = {u[1] * v[2] - u[2] * v[1], u[2] * v[0] - u[0] * v[2],
                                       u[0] * v[1] - u[1] * v[0]};
      const std::array<double, 3> g = {w[1] * B[2] - w[2] * B[1], w[2] * B[0] - w[0] * B[2],
                                       w[0] * B[1] - w[1] * B[0]};
      // anchor - anchor* = +0.1 u for the ground anchor (s_a = 0.5 R = 0.1 into the
      // unbounded ground).
      const double c0 = -0.1 * (g[0] * u[0] + g[1] * u[1] + g[2] * u[2]);
      INFO("patch " << shifted.patch << " c0 " << c0);
      REQUIRE(std::abs(c0) > 1.0e-3);
      for (int i = 0; i < 4; i++)
      {
        CHECK_THAT(shifted.coefficients[i] - reference.coefficients[i],
                   WithinAbs(c0, 1.0e-10));
      }
      // The strip's and the right ground's anchors (their edge points 0.125 / 0.5) move
      // +0.1 u into their piece: -(w x B) . (0.1 u) = c0 more on each conductor state.
      CHECK_THAT(shifted.coefficients[4] - reference.coefficients[4],
                 WithinAbs(2.0 * c0, 1.0e-10));
      CHECK_THAT(shifted.coefficients[5] - reference.coefficients[5],
                 WithinAbs(2.0 * c0, 1.0e-10));
    }
  }
  // dK under the leaf: symmetric, and the postprocessed correction (the chain + the
  // magnetic loop fit, no transpose) = 1/2 <x, dK x> (the forward + the transpose on the
  // twin lines): the transpose is the adjoint.
  {
    SurfaceResponseStiffnessOperator dK(*leaf.response);
    Vector y(x.Size()), z(x.Size()), dKz(x.Size());
    dK.Mult(x, y);
    const double form = linalg::Dot(leaf.space->GetComm(), x, y);
    CHECK(form > 0.0);
    z.Randomize(11);
    z.SetSubVector(leaf.space->GetNDDbcTDofLists().back(), 0.0);
    dK.Mult(z, dKz);
    CHECK_THAT(linalg::Dot(leaf.space->GetComm(), z, y),
               WithinRel(linalg::Dot(leaf.space->GetComm(), x, dKz), 1.0e-9));
    CHECK(linalg::Dot(leaf.space->GetComm(), z, dKz) >= -1.0e-12 * std::abs(form));
    const auto result = leaf.response->GetMaxwellResponse(field, 1.0);
    REQUIRE(result.has_magnetic);
    CHECK(result.magnetic_fallback_patch_count == 0);
    CHECK(result.magnetic_excluded_patch_count == 0);
    CHECK_THAT(result.magnetic_domain_correction, WithinRel(0.5 * form, 1.0e-9));
    // On the EXACT linear field (the sheet dofs not zeroed: the zeroing bends the discrete
    // curl in the sheet's elements, which the twin lines starting on the metal cross) every
    // twin loop closes: the magnetic loop residual at roundoff under the leaf as in the
    // legacy form.
    GridFunction exact_field(leaf.space->GetNDSpace(), true);
    exact_field.Real().SetFromTrueDofs(x_exact);
    exact_field.Imag() = 0.0;
    CHECK(leaf.response->GetMaxwellResponse(exact_field, 1.0).magnetic_loop_residual <
          1.0e-10);
    CHECK(legacy.response->GetMaxwellResponse(exact_field, 1.0).magnetic_loop_residual <
          1.0e-10);
    // The legacy form reads a different correction (the states moved by c0).
    SurfaceResponseStiffnessOperator dK_legacy(*legacy.response);
    Vector y_legacy(x.Size());
    dK_legacy.Mult(x, y_legacy);
    CHECK(std::abs(linalg::Dot(legacy.space->GetComm(), x, y_legacy) - form) >
          1.0e-6 * std::abs(form));
  }
  // The gate under the leaf on the gap-flux field (the nudge of record): every ratio at the
  // hexahedral fixture's blind-quadrature bound and uniform across the cells (as the 1302
  // control read it; the magnetic loop fit now applied).
  {
    const auto gated_library = WriteStubMagneticLibrary(
        temp.temp_dir, convex_library_3d_path, "pieces-stateanchor-gate", true,
        {{"StateAnchorSetbackOverR", 0.5}});
    auto gated = Make(gated_library);
    auto gap = GapFluxField(1.0);
    gated.response->ApplyMagneticConductorConsistencyGate(ProjectField(*gated.space, gap),
                                                          0, 5.0);
    const auto &records = gated.response->GetMagneticGateRecords();
    REQUIRE(records.size() == static_cast<std::size_t>(gated.response->GetPatchCount()));
    double lowest = mfem::infinity(), highest = 0.0;
    for (const auto &record : records)
    {
      INFO("patch " << record.patch << " ratio " << record.max_ratio);
      CHECK(record.max_ratio < kHexSlotEdgeQuadratureBound);
      lowest = std::min(lowest, record.max_ratio);
      highest = std::max(highest, record.max_ratio);
    }
    CHECK(highest - lowest < 1.0e-6);
    // The 1302 CONTROL on a library whose leaf already re-anchors: refused by name.
    json control = StubMaxwellConfig(*this, gated_library);
    control["Solver"]["SurfaceResponseCorrection"]["CorrectionMode"] = "PostprocessOnly";
    control["Solver"]["SurfaceResponseCorrection"]["ConductorConsistency"] = {
        {"RecordOnly", true}, {"ControlStateAnchorSetbackOverR", 0.5}};
    std::vector<std::unique_ptr<Mesh>> meshes;
    IoData iodata = Build(control, meshes);
    SpaceOperator space(iodata, meshes);
    CHECK_THROWS_WITH(SurfaceResponseOperator(iodata, space),
                      Catch::Matchers::ContainsSubstring("StateAnchorSetbackOverR"));
  }
  // The leaf's bounds: outside [0, 1] = malformed, refused by name.
  {
    const auto bad = WriteStubMagneticLibrary(temp.temp_dir, convex_library_3d_path,
                                              "pieces-stateanchor-bad", true,
                                              {{"StateAnchorSetbackOverR", 1.5}});
    std::vector<std::unique_ptr<Mesh>> meshes;
    IoData iodata = Build(StubMaxwellConfig(*this, bad), meshes);
    SpaceOperator space(iodata, meshes);
    CHECK_THROWS_WITH(SurfaceResponseOperator(iodata, space),
                      Catch::Matchers::ContainsSubstring("StateAnchorSetbackOverR"));
  }
#endif
}

// Decision 1323 (the (b) W06 line): the state anchor's EDGE SNAP. A conductor anchor is the
// device's edge point while the patch's conductor_probe_edges carry the device sides' feet;
// an anchor within StateAnchorEdgeToleranceOverR x R of an edge is ON that edge (its piece,
// the setback measured from the edge offset), beyond it the anchor must lie in metal
// (refused by name); the tolerance must stay below half the narrowest model gap / width
// (CONDITION 1, asserted at load); the on-edge / snapped counts and the farthest snap are
// recorded (CONDITION 2). On the stub every anchor sits exactly on its edge: all on-edge,
// none snapped from a nonzero distance, max 0 = bitwise to the R-MAG-17 form.
TEST_CASE("SnapStateAnchorToEdge: the nearest edge within the tolerance, metal beyond it",
          "[surfaceresponseoperator][magnetic][probeform][stateanchor][Serial]")
{
  // The stub's cross-section: edges 0 / 0.125 / 0.375 / 0.5 (ground | slot A | strip |
  // slot B | ground), pieces 1 / 2 / 2 / 3.
  const std::vector<double> offsets = {0.0, 0.125, 0.375, 0.5};
  const std::vector<int> gap_directions = {1, -1, 1, -1};
  const std::vector<int> pieces = {1, 2, 2, 3};
  // An anchor a hair into slot A from the ground edge (the near-matched DCG case): snapped
  // onto edge 0 (piece 1) at its distance.
  auto snap = SnapStateAnchorToEdge(offsets, gap_directions, pieces, 0.004, 0.02);
  CHECK(snap.snapped);
  CHECK(snap.piece == 1);
  CHECK(snap.u_edge == 0.0);
  CHECK_THAT(snap.distance, WithinAbs(0.004, 1.0e-15));
  // Into slot A from the strip edge: snapped onto edge 1 (piece 2).
  snap = SnapStateAnchorToEdge(offsets, gap_directions, pieces, 0.12, 0.02);
  CHECK((snap.snapped && snap.piece == 2 && snap.u_edge == 0.125));
  // Exactly on an edge: on the edge at distance 0.
  snap = SnapStateAnchorToEdge(offsets, gap_directions, pieces, 0.375, 0.02);
  CHECK((snap.snapped && snap.piece == 2 && snap.distance == 0.0));
  // Beyond the tolerance inside the slot: no piece (the caller refuses).
  snap = SnapStateAnchorToEdge(offsets, gap_directions, pieces, 0.06, 0.02);
  CHECK_FALSE(snap.snapped);
  CHECK(snap.piece == 0);
  CHECK_THAT(snap.distance, WithinAbs(0.06, 1.0e-15));
  // Beyond the tolerance inside metal: the metal piece, unsnapped, u_edge = u0.
  snap = SnapStateAnchorToEdge(offsets, gap_directions, pieces, 0.25, 0.02);
  CHECK((!snap.snapped && snap.piece == 2 && snap.u_edge == 0.25));
  snap = SnapStateAnchorToEdge(offsets, gap_directions, pieces, -0.3, 0.02);
  CHECK((!snap.snapped && snap.piece == 1));
  // The setback measured from the edge: edge 0's piece (the unbounded ground) moves the
  // anchor to -s from the EDGE, not from u0.
  snap = SnapStateAnchorToEdge(offsets, gap_directions, pieces, 0.004, 0.02);
  const double move =
      (snap.u_edge - 0.004) + LateralSetbackMove(offsets, gap_directions, pieces,
                                                 snap.u_edge, snap.piece, 0.1, 0.45);
  CHECK_THAT(0.004 + move, WithinAbs(-0.1, 1.0e-15));
  // Decision 1334 (MAJOR-2): SetBackStateAnchor = the ONE rule shared by the magnetic
  // term's state anchors and the LIFTED probe form's anchor end, the expression above
  // verbatim: an anchor a hair into slot A from the ground edge (the near-matched DCG
  // case) is set back s = 0.5 R = 0.1 from the EDGE into the unbounded ground.
  const StateAnchorSetback ground =
      SetBackStateAnchor(offsets, gap_directions, pieces, 0.004, 0.02, 0.1, 0.45);
  CHECK((ground.snap.snapped && ground.snap.piece == 1 && ground.snap.u_edge == 0.0));
  CHECK(ground.move == move);
  CHECK_THAT(0.004 + ground.move, WithinAbs(-0.1, 1.0e-15));
  // A hair into slot A from the strip edge and a hair into the strip metal: both snap onto
  // the strip edge 0.125 (piece 2) and land at the SAME set-back point 0.225 (the strip's
  // width 0.25 caps the setback at 0.45 x 0.25 = 0.1125 >= 0.1: s applies).
  const StateAnchorSetback gap_side =
      SetBackStateAnchor(offsets, gap_directions, pieces, 0.121, 0.02, 0.1, 0.45);
  const StateAnchorSetback metal_side =
      SetBackStateAnchor(offsets, gap_directions, pieces, 0.129, 0.02, 0.1, 0.45);
  CHECK((gap_side.snap.snapped && gap_side.snap.piece == 2));
  CHECK((metal_side.snap.snapped && metal_side.snap.piece == 2));
  CHECK_THAT(0.121 + gap_side.move, WithinAbs(0.225, 1.0e-15));
  CHECK_THAT(0.129 + metal_side.move, WithinAbs(0.225, 1.0e-15));
  // Beyond the tolerance inside the slot: no piece and no move (the caller refuses).
  const StateAnchorSetback in_gap =
      SetBackStateAnchor(offsets, gap_directions, pieces, 0.06, 0.02, 0.1, 0.45);
  CHECK((!in_gap.snap.snapped && in_gap.snap.piece == 0 && in_gap.move == 0.0));
  // BITWISE with the legacy probe-anchor expression where no anchor snaps from a nonzero
  // distance: exactly on an edge (u_edge == u0) and well inside metal (unsnapped), the
  // move IS LateralSetbackMove of the anchor's own piece ((u0 - u0) + x == x).
  for (const double u0 : {0.0, 0.125, 0.375, 0.5, 0.25, -0.3, 0.7})
  {
    const StateAnchorSetback exact =
        SetBackStateAnchor(offsets, gap_directions, pieces, u0, 0.02, 0.1, 0.45);
    const int piece = LocalPieceOfOffset(offsets, gap_directions, pieces, u0);
    REQUIRE(piece > 0);
    CHECK(exact.snap.piece == piece);
    CHECK(exact.move ==
          LateralSetbackMove(offsets, gap_directions, pieces, u0, piece, 0.1, 0.45));
  }
  // ProbeSetbackOverR 0 (the lift-only form) moves nothing on a metal anchor.
  CHECK(SetBackStateAnchor(offsets, gap_directions, pieces, 0.375, 0.02, 0.0, 0.45).move ==
        0.0);
}

TEST_CASE_METHOD(test::SurfaceResponseFiles,
                 "SurfaceResponseOperator magnetic term: the state-anchor edge tolerance "
                 "(MagneticRules.StateAnchorEdgeToleranceOverR, decision 1323)",
                 "[surfaceresponseoperator][maxwell][magnetic][probeform][stateanchor]"
                 "[Serial][Parallel]")
{
#if !defined(MFEM_USE_GSLIB)
  SKIP("SurfaceResponseOperator requires MFEM_USE_GSLIB");
#else
  auto Build = [&](const json &run_config, std::vector<std::unique_ptr<Mesh>> &meshes)
  {
    IoData iodata(run_config, false);
    iodata.boundaries.cracked_attributes.insert(9);
    mfem::Mesh serial = test::MakeSerialShortedCpwMesh();
    meshes.push_back(
        std::make_unique<Mesh>(std::make_unique<mfem::ParMesh>(Mpi::World(), serial)));
    return iodata;
  };
  // The stub's narrowest gap 0.125 at R 0.2 = 0.625 R: the default tolerance 0.1 R is
  // admitted; on the stub every conductor anchor is exactly on its edge.
  {
    const auto library = WriteStubMagneticLibrary(temp.temp_dir, convex_library_3d_path,
                                                  "pieces-edgetol-default", true,
                                                  {{"StateAnchorSetbackOverR", 0.5}});
    std::vector<std::unique_ptr<Mesh>> meshes;
    IoData iodata = Build(StubMaxwellConfig(*this, library), meshes);
    SpaceOperator space(iodata, meshes);
    SurfaceResponseOperator response(iodata, space);
    REQUIRE(response.GetPatchCount() > 0);
    const auto &rules = response.GetMagneticProbeRules();
    CHECK(rules.state_anchor_edge_tolerance_over_r == 0.1);
    CHECK_THAT(rules.state_anchor_min_gap_over_r, WithinRel(0.625, 1.0e-12));
    const auto statistics = response.GetStatistics();
    const auto &snapped = statistics["Correction"]["MagneticStateAnchorSnapped"];
    // 3 conductor-state anchors per patch (the ground, the strip, the right ground), each
    // on its edge; none moved by the snap.
    CHECK(snapped["OnEdge"].get<int>() ==
          3 * statistics["Correction"]["MagneticStateAnchorPatches"]["Total"].get<int>());
    CHECK(snapped["Count"].get<int>() == 0);
    CHECK(snapped["MaxOverR"].get<double>() == 0.0);
    CHECK(snapped["ToleranceOverR"].get<double>() == 0.1);
  }
  // CONDITION 1: a tolerance at or above half the narrowest gap (0.3125 R) is refused at
  // load by name when EITHER setback is in effect (decision 1334: the term's
  // StateAnchorSetbackOverR or the probe form's ProbeSetbackOverR, whose default 0.5
  // alone triggers the check); a library with BOTH setbacks 0 is not checked.
  {
    const auto bad = WriteStubMagneticLibrary(
        temp.temp_dir, convex_library_3d_path, "pieces-edgetol-wide", true,
        {{"StateAnchorSetbackOverR", 0.5}, {"StateAnchorEdgeToleranceOverR", 0.4}});
    std::vector<std::unique_ptr<Mesh>> meshes;
    IoData iodata = Build(StubMaxwellConfig(*this, bad), meshes);
    SpaceOperator space(iodata, meshes);
    CHECK_THROWS_WITH(SurfaceResponseOperator(iodata, space),
                      Catch::Matchers::ContainsSubstring("StateAnchorEdgeToleranceOverR"));
    const auto bad_probe = WriteStubMagneticLibrary(
        temp.temp_dir, convex_library_3d_path, "pieces-edgetol-wide-probe", true,
        {{"ProbeSetbackOverR", 0.5}, {"StateAnchorEdgeToleranceOverR", 0.4}});
    std::vector<std::unique_ptr<Mesh>> probe_meshes;
    IoData probe_iodata = Build(StubMaxwellConfig(*this, bad_probe), probe_meshes);
    SpaceOperator probe_space(probe_iodata, probe_meshes);
    CHECK_THROWS_WITH(SurfaceResponseOperator(probe_iodata, probe_space),
                      Catch::Matchers::ContainsSubstring("StateAnchorEdgeToleranceOverR"));
    const auto legacy = WriteStubMagneticLibrary(
        temp.temp_dir, convex_library_3d_path, "pieces-edgetol-legacy", true,
        {{"StateAnchorEdgeToleranceOverR", 0.4}, {"ProbeSetbackOverR", 0.0}});
    std::vector<std::unique_ptr<Mesh>> legacy_meshes;
    IoData legacy_iodata = Build(StubMaxwellConfig(*this, legacy), legacy_meshes);
    SpaceOperator legacy_space(legacy_iodata, legacy_meshes);
    SurfaceResponseOperator legacy_response(legacy_iodata, legacy_space);
    CHECK(legacy_response.GetPatchCount() > 0);
    CHECK(legacy_response.GetStatistics()["Correction"]["MagneticStateAnchorSetbackOverR"]
              .get<double>() == 0.0);
  }
  // The leaf's bounds: outside [0, 1] = malformed, refused by name.
  {
    const auto bad = WriteStubMagneticLibrary(
        temp.temp_dir, convex_library_3d_path, "pieces-edgetol-bad", true,
        {{"StateAnchorSetbackOverR", 0.5}, {"StateAnchorEdgeToleranceOverR", 1.5}});
    std::vector<std::unique_ptr<Mesh>> meshes;
    IoData iodata = Build(StubMaxwellConfig(*this, bad), meshes);
    SpaceOperator space(iodata, meshes);
    CHECK_THROWS_WITH(SurfaceResponseOperator(iodata, space),
                      Catch::Matchers::ContainsSubstring("StateAnchorEdgeToleranceOverR"));
  }
#endif
}

// Phase 3 (a) t6: the probe-form leaves (R-MAG-10 / 12): absent = the defaults (LIFTED
// 0.5 / 0.5 / 0.45 / 1e-6 / 1), 0 / 0 = FLAT, a Gauss-height count >= 2 = BAND, and a value
// outside its bounds is MALFORMED (refused at load by name).
TEST_CASE_METHOD(
    test::SurfaceResponseFiles,
    "SurfaceResponseOperator magnetic probe-form leaves: defaults, FLAT, BAND, "
    "malformed (Phase 3 (a))",
    "[surfaceresponseoperator][3d][maxwell][magnetic][probeform][Serial]")
{
#if !defined(MFEM_USE_GSLIB)
  SKIP("SurfaceResponseOperator requires MFEM_USE_GSLIB");
#else
  auto Load = [&](const std::string &tag, const json &rules)
  {
    const auto files = WriteMagneticLibrary(library_3d_path, temp.temp_dir, tag, kPsdDefect,
                                            DefaultMagneticSource(),
                                            [&](json &library)
                                            {
                                              if (!rules.empty())
                                              {
                                                library["MagneticRules"] = rules;
                                              }
                                            });
    IoData iodata(MaxwellCpw3dConfig(*this, files.library), false);
    auto device_mesh = mesh::ReadMesh(iodata, Mpi::World());
    std::vector<std::unique_ptr<Mesh>> meshes;
    meshes.push_back(std::make_unique<Mesh>(std::move(device_mesh)));
    SpaceOperator space(iodata, meshes);
    SurfaceResponseOperator response(iodata, space);
    return std::make_pair(response.MagneticProbeForm(), response.GetMagneticProbeRules());
  };
  {
    const auto [form, rules] = Load("defaults", json::object());
    CHECK(form == "LIFTED");
    CHECK(rules.probe_lift_over_r == 0.5);
    CHECK(rules.probe_setback_over_r == 0.5);
    CHECK(rules.probe_setback_width_fraction == 0.45);
    CHECK(rules.probe_sheet_nudge_over_r == 1.0e-6);
    CHECK(rules.probe_band_heights == 1);
  }
  {
    const auto [form, rules] =
        Load("flat", {{"ProbeLiftOverR", 0.0}, {"ProbeSetbackOverR", 0.0}});
    CHECK(form == "FLAT");
    CHECK(rules.probe_lift_over_r == 0.0);
  }
  {
    const auto [form, rules] = Load("band", {{"ProbeBandHeights", 3},
                                             {"ProbeSetbackWidthFraction", 0.3},
                                             {"ProbeSheetNudgeOverR", 1.0e-5}});
    CHECK(form == "BAND");
    CHECK(rules.probe_band_heights == 3);
    CHECK(rules.probe_setback_width_fraction == 0.3);
    CHECK(rules.probe_sheet_nudge_over_r == 1.0e-5);
  }
  for (const auto &[tag, rules] : std::vector<std::pair<std::string, json>>{
           {"lift-high", {{"ProbeLiftOverR", 1.5}}},
           {"setback-negative", {{"ProbeSetbackOverR", -0.1}}},
           {"fraction-half", {{"ProbeSetbackWidthFraction", 0.5}}},
           {"nudge-negative", {{"ProbeSheetNudgeOverR", -1.0e-6}}},
           {"heights-zero", {{"ProbeBandHeights", 0}}},
           {"unknown-key", {{"ProbeLift", 0.5}}}})
  {
    INFO(tag);
    CHECK_THROWS_WITH(Load(tag, rules),
                      Catch::Matchers::ContainsSubstring(rules.begin().key().substr(0, 9)));
  }
#endif
}

// Decision 1309 Q3 (USER 1308 U-1; OI-MAGNETIC-DUAL-GATE-3D-ENCLOSED-FLUX): the dual gate's
// ConductorConsistencyTolerance is a scalar (one value at every order, as before) or a
// PER-ORDER object {"<order>": tolerance}: the solve's Order (the 3D CPW config runs at
// Order 1) selects its entry; an order absent from the object is refused by name; a
// malformed key or value is refused at load. The value in effect is recorded.
TEST_CASE_METHOD(test::SurfaceResponseFiles,
                 "SurfaceResponseOperator magnetic gate tolerance: scalar or per-order "
                 "(ConductorConsistencyTolerance, decision 1309)",
                 "[surfaceresponseoperator][3d][maxwell][magnetic][probeform][Serial]")
{
#if !defined(MFEM_USE_GSLIB)
  SKIP("SurfaceResponseOperator requires MFEM_USE_GSLIB");
#else
  auto Load = [&](const std::string &tag, const json &rules)
  {
    const auto files = WriteMagneticLibrary(library_3d_path, temp.temp_dir, tag, kPsdDefect,
                                            DefaultMagneticSource(), [&](json &library)
                                            { library["MagneticRules"] = rules; });
    IoData iodata(MaxwellCpw3dConfig(*this, files.library), false);
    REQUIRE(iodata.solver.order == 1);
    auto device_mesh = mesh::ReadMesh(iodata, Mpi::World());
    std::vector<std::unique_ptr<Mesh>> meshes;
    meshes.push_back(std::make_unique<Mesh>(std::move(device_mesh)));
    SpaceOperator space(iodata, meshes);
    SurfaceResponseOperator response(iodata, space);
    const auto statistics = response.GetStatistics();
    CHECK(statistics["Correction"]["MagneticGateTolerance"].get<double>() ==
          response.MagneticConductorConsistencyTolerance());
    return std::make_pair(response.MagneticConductorConsistencyTolerance(),
                          response.GetMagneticProbeRules());
  };
  {
    const auto [tolerance, rules] = Load("tol-absent", json::object());
    CHECK(tolerance == SurfaceResponseOperator::kConductorConsistencyTolerance);
    CHECK(rules.conductor_consistency_tolerance_by_order.empty());
  }
  {
    const auto [tolerance, rules] =
        Load("tol-scalar", {{"ConductorConsistencyTolerance", 0.25}});
    CHECK(tolerance == 0.25);
    CHECK(rules.conductor_consistency_tolerance == 0.25);
    CHECK(rules.conductor_consistency_tolerance_by_order.empty());
  }
  {
    const auto [tolerance, rules] =
        Load("tol-per-order",
             {{"ConductorConsistencyTolerance", {{"1", 0.07}, {"3", 0.13}, {"4", 0.08}}}});
    CHECK(tolerance == 0.07);
    CHECK(rules.conductor_consistency_tolerance < 0.0);
    REQUIRE(rules.conductor_consistency_tolerance_by_order.size() == 3);
    CHECK(rules.conductor_consistency_tolerance_by_order.at(3) == 0.13);
  }
  CHECK_THROWS_WITH(Load("tol-missing-order",
                         {{"ConductorConsistencyTolerance", {{"3", 0.13}, {"4", 0.08}}}}),
                    Catch::Matchers::ContainsSubstring("no entry for the solve's Order 1"));
  for (const auto &[tag, rules] : std::vector<std::pair<std::string, json>>{
           {"tol-bad-key", {{"ConductorConsistencyTolerance", {{"p3", 0.1}}}}},
           {"tol-zero-key", {{"ConductorConsistencyTolerance", {{"0", 0.1}}}}},
           {"tol-negative-entry", {{"ConductorConsistencyTolerance", {{"1", -0.1}}}}},
           {"tol-empty", {{"ConductorConsistencyTolerance", json::object()}}},
           {"tol-negative", {{"ConductorConsistencyTolerance", -0.1}}},
           {"tol-string", {{"ConductorConsistencyTolerance", "0.1"}}}})
  {
    INFO(tag);
    CHECK_THROWS_WITH(Load(tag, rules),
                      Catch::Matchers::ContainsSubstring("ConductorConsistencyTolerance"));
  }
#endif
}

// Phase 3 (b) COL-3 (R-MAG-15): the within-R VOLUME energies per patch box on the stub (the
// one-label 4-edge PEC cluster: the box u in [0 - R, 0.5 + R], v in [-R, R], w = the cell).
// Uniform E and B: the densities are exact at every quadrature point (energy / resolved
// volume = the density to roundoff), the resolved volume is the exact box volume to the
// quadrature resolution (the boundary layer of one element, h = 0.125, around the box).
TEST_CASE_METHOD(test::SurfaceResponseFiles,
                 "SurfaceResponseOperator within-R volume energies per patch box (COL-3, "
                 "Phase 3 (b))",
                 "[surfaceresponseoperator][maxwell][magnetic][withinr][Serial][Parallel]")
{
#if !defined(MFEM_USE_GSLIB)
  SKIP("SurfaceResponseOperator requires MFEM_USE_GSLIB");
#else
  const auto library =
      WriteStubMagneticLibrary(temp.temp_dir, convex_library_3d_path, "withinr", false);
  json config = StubMaxwellConfig(*this, library);
  // The leaf without PatchEnergy is refused; with it, read.
  {
    json bad = config;
    bad["Solver"]["SurfaceResponseCorrection"]["Diagnostics"] = {
        {"WithinRVolumeEnergy", true}};
    CHECK_THROWS_WITH(IoData(bad, false),
                      Catch::Matchers::ContainsSubstring("PatchEnergy"));
  }
  config["Solver"]["SurfaceResponseCorrection"]["PatchEnergy"] = true;
  config["Solver"]["SurfaceResponseCorrection"]["Diagnostics"] = {
      {"WithinRVolumeEnergy", true}};
  IoData iodata(config, false);
  iodata.boundaries.cracked_attributes.insert(9);
  std::vector<std::unique_ptr<Mesh>> meshes;
  mfem::Mesh serial = test::MakeSerialShortedCpwMesh();
  meshes.push_back(
      std::make_unique<Mesh>(std::make_unique<mfem::ParMesh>(Mpi::World(), serial)));
  SpaceOperator space(iodata, meshes);
  SurfaceResponseOperator response(iodata, space);
  REQUIRE(response.GetPatchCount() > 0);
  CHECK(response.WithinRVolumeEnergyRequested());
  const std::array<double, 3> E0 = {0.7, -0.2, 0.4}, B0 = {0.3, 0.9, -0.5};
  mfem::VectorConstantCoefficient e_coefficient(
      mfem::Vector(const_cast<double *>(E0.data()), 3));
  mfem::VectorConstantCoefficient b_coefficient(
      mfem::Vector(const_cast<double *>(B0.data()), 3));
  GridFunction E(space.GetNDSpace(), true), B(space.GetRTSpace(), true);
  E.Real().ProjectCoefficient(e_coefficient);
  E.Imag() = 0.0;
  B.Real().ProjectCoefficient(b_coefficient);
  B.Imag() = 0.0;
  const auto &mat_op = space.GetMaterialOp();
  const mfem::DenseMatrix epsilon = mat_op.GetPermittivityReal(1);
  const mfem::DenseMatrix muinv = mat_op.GetInvPermeability(1);
  mfem::Vector e(const_cast<double *>(E0.data()), 3), b(const_cast<double *>(B0.data()), 3),
      tmp(3);
  epsilon.Mult(e, tmp);
  const double electric_density = 0.5 * (e * tmp);
  muinv.Mult(b, tmp);
  const double magnetic_density = 0.5 * (b * tmp);
  const auto energies = response.GetWithinRVolumeEnergies(E, B, mat_op);
  REQUIRE(energies.size() == static_cast<std::size_t>(response.GetPatchCount()));
  // The local patches' boxes give the exact volumes (every rank checks its own patches
  // against the gathered, sorted global list).
  const auto &probes = response.GetMaxwellProbes();
  constexpr double R = 0.2, h = 0.125;
  double total_volume = 0.0, total_exact = 0.0;
  int resolved = 0, unresolved = 0;
  for (const auto &entry : energies)
  {
    REQUIRE(entry.defined);
    if (entry.volume == 0.0)
    {
      // A cell thinner than the quadrature spacing resolves no point (the quadrature-
      // resolution limitation of the estimate, named in the schema): its energies are 0.
      CHECK(entry.electric == 0.0);
      CHECK(entry.magnetic == 0.0);
      unresolved++;
      continue;
    }
    resolved++;
    CHECK_THAT(entry.electric / entry.volume, WithinRel(electric_density, 1.0e-10));
    CHECK_THAT(entry.magnetic / entry.volume, WithinRel(magnetic_density, 1.0e-10));
    total_volume += entry.volume;
  }
  CHECK(resolved > unresolved);
  // The exact box volume of every LOCAL patch: (0.5 + 2R) x 2R x the cell length (the two
  // end-face probes of edge 0 give the cell along x).
  for (const auto &patch_probes : probes)
  {
    double x0 = mfem::infinity(), x1 = -mfem::infinity();
    for (const auto &probe : patch_probes)
    {
      x0 = std::min(x0, probe.point[0]);
      x1 = std::max(x1, probe.point[0]);
    }
    const double cell = x1 - x0;
    REQUIRE(cell > 0.0);
    total_exact += (0.5 + 2.0 * R) * (2.0 * R) * cell;
  }
  Mpi::GlobalSum(1, &total_exact, Mpi::World());
  // The quadrature resolves the box faces to one element: the boundary layer bound h x the
  // boxes' total surface (the cross-section A = (0.5 + 2R) x 2R, the placed length L).
  const double A = (0.5 + 2.0 * R) * (2.0 * R);
  const double L = total_exact / A;
  const double surface = 2.0 * A + 2.0 * ((0.5 + 2.0 * R) + 2.0 * R) * L;
  CHECK(std::abs(total_volume - total_exact) <= h * surface);
  CHECK(total_volume > 0.5 * total_exact);
#endif
}

TEST_CASE("FixedTraceFrequencyShift on the Phase-0 numbers", "[magnetic][Serial]")
{
  struct Cell
  {
    double dC_over_C;      // N1 = a^T dQ a / W_full
    double dC0_over_C0;    // N2 = a_m^T dS_m a_m / W_full,vac
    double expected;       // the REPORT's MAGNETIC (TRUE) df/f
    double electric_only;  // the REPORT's ELECTRIC-ONLY df/f
  };
  const std::array<Cell, 2> cells = {{{-2.434e-2, 4.080e-2, 3.257e-2, 1.217e-2},
                                      {-3.378e-2, 5.790e-2, 4.584e-2, 1.689e-2}}};
  // An arbitrary raw energy (nondimensional): the shift is scale-free.
  const double raw_energy = 2.0811e-13 * 7.3;
  for (const auto &cell : cells)
  {
    SurfaceResponseOperator::MaxwellResponse response;
    response.has_magnetic = true;
    response.closure_independent_confident = true;
    response.magnetic_confident = true;
    response.confident = false;  // the trace-closure spread fails: NOT a gate of D6
    response.maximum_trace_closure_spread = 0.085;
    response.domain_correction = cell.dC_over_C * raw_energy;
    response.magnetic_domain_correction = cell.dC0_over_C0 * raw_energy;
    const double shift = FixedTraceFrequencyShift(response, raw_energy);
    CHECK_THAT(shift, WithinAbs(-0.5 * (cell.dC_over_C - cell.dC0_over_C0), 1.0e-6));
    CHECK_THAT(shift, WithinAbs(cell.expected, 1.0e-5));
    // The electric half alone (no magnetic leaf applied) is the electric-only reading.
    SurfaceResponseOperator::MaxwellResponse electric = response;
    electric.magnetic_domain_correction = 0.0;
    CHECK_THAT(FixedTraceFrequencyShift(electric, raw_energy),
               WithinAbs(cell.electric_only, 1.0e-5));
    // The gates: magnetic_confident and has_magnetic withhold the shift; decision 1250:
    // closure_independent_confident does NOT (a device with ends: the (ii) stub's short
    // and open end are 0.95 % unmatched, so the closure-independent gate's coverage
    // criterion >= 1 - 1e-10 fails while the shift is read unchanged, the flag beside it).
    SurfaceResponseOperator::MaxwellResponse gated = response;
    gated.closure_independent_confident = false;
    gated.matched_length_fraction = 0.990517;
    CHECK_THAT(FixedTraceFrequencyShift(gated, raw_energy), WithinAbs(shift, 0.0));
    gated.confident = false;
    gated.maximum_trace_closure_spread = 0.0834;
    CHECK_THAT(FixedTraceFrequencyShift(gated, raw_energy), WithinAbs(shift, 0.0));
    gated = response;
    gated.magnetic_confident = false;
    CHECK(std::isnan(FixedTraceFrequencyShift(gated, raw_energy)));
    gated = response;
    gated.has_magnetic = false;
    CHECK(std::isnan(FixedTraceFrequencyShift(gated, raw_energy)));
  }
}

// Decision 1250: the matched energy share the fixed-trace shift rests on = 1 - (the raw
// within-R energy kept raw by the electric correction) / (the raw within-R energy of the
// target interfaces): 1 at full coverage (bitwise), the (ii) stub's 0.95 % unmatched ends
// read as their energy share, NaN without within-R energy.
TEST_CASE("FixedTraceShiftCoverage", "[magnetic][Serial]")
{
  const double within = 2.0811e-13;
  CHECK(FixedTraceShiftCoverage(within, 0.0) == 1.0);
  CHECK_THAT(FixedTraceShiftCoverage(within, 0.0095 * within), WithinAbs(0.9905, 1.0e-12));
  CHECK_THAT(FixedTraceShiftCoverage(1.0, 0.25), WithinAbs(0.75, 0.0));
  CHECK(std::isnan(FixedTraceShiftCoverage(0.0, 0.0)));
  CHECK(std::isnan(FixedTraceShiftCoverage(-1.0, 0.0)));
}

// The Phase-1 gate G-d loader read (prereg-gates.md): the leaves of the library named by
// PALACE_MAGNETIC_LIBRARY as the operator loads them (r1-r4 with the library's
// MagneticRules, the PSD statistic, the scaled defect's norms), written to
// PALACE_MAGNETIC_LIBRARY_RECORD for the lane's Python reader to cross-check against its
// own evaluation of the same files.
TEST_CASE("DescribeMagneticLeaves of PALACE_MAGNETIC_LIBRARY", "[magnetic][Serial]")
{
  const char *library = std::getenv("PALACE_MAGNETIC_LIBRARY");
  const char *record = std::getenv("PALACE_MAGNETIC_LIBRARY_RECORD");
  if (!library || !record)
  {
    SKIP("PALACE_MAGNETIC_LIBRARY / PALACE_MAGNETIC_LIBRARY_RECORD are not set");
  }
  const auto description = DescribeMagneticLeaves(library, Units(1.0, 1.0), false);
  REQUIRE(description.contains("Models"));
  CHECK(!description["Models"].empty());
  for (const auto &model : description["Models"])
  {
    INFO(model.dump());
    CHECK(model["Status"].get<std::string>() != "Refused");
  }
  if (Mpi::Root(Mpi::World()))
  {
    std::ofstream output(record);
    output << description.dump(1) << "\n";
  }
}

}  // namespace palace
