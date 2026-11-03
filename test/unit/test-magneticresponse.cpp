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
#include <string>
#include <vector>
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
  }
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
    // The gates: closure_independent_confident, magnetic_confident, has_magnetic.
    SurfaceResponseOperator::MaxwellResponse gated = response;
    gated.closure_independent_confident = false;
    CHECK(std::isnan(FixedTraceFrequencyShift(gated, raw_energy)));
    gated = response;
    gated.magnetic_confident = false;
    CHECK(std::isnan(FixedTraceFrequencyShift(gated, raw_energy)));
    gated = response;
    gated.has_magnetic = false;
    CHECK(std::isnan(FixedTraceFrequencyShift(gated, raw_energy)));
  }
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
