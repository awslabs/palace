// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include <array>
#include <cmath>
#include <complex>
#include <filesystem>
#include <fstream>
#include <memory>
#include <string>
#include <utility>
#include <vector>
#include <mfem.hpp>
#include <catch2/catch_test_macros.hpp>
#include <nlohmann/json.hpp>
#include <catch2/generators/catch_generators.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include "driver.hpp"
#include "fem/mesh.hpp"
#include "linalg/vector.hpp"
#include "models/postoperator.hpp"
#include "models/spaceoperator.hpp"
#include "surfaceresponse-fixtures.hpp"
#include "utils/communication.hpp"
#include "utils/geodata.hpp"
#include "utils/iodata.hpp"
#include "utils/omp.hpp"
#include "utils/outputdir.hpp"
#include "utils/tablecsv.hpp"
#include "utils/timer.hpp"

namespace palace
{

using json = nlohmann::json;
using Catch::Matchers::WithinRel;
namespace fs = std::filesystem;

// P4 S3-0 (review r2084 F-6; decision 2090): the self-consistent driven path through the
// last frequency. The convex PEC island of the fixtures (written to a mesh file) driven by
// a lumped port on the exterior face z = 0 (attribute 1), the island sheet (attribute 9) in
// Ground, a uniform sweep of four frequencies below the box's lowest mode (about 0.1 GHz),
// CorrectionMode Both: the corrected solve runs at every frequency, and after the last one
// the driver releases its solvers (the SY1 resolution of upstream #1087 in
// DrivenSolver::SweepUniform). The references are our main's values: the same problem run
// by the reference build P of main 81425c2d40 (palace-cd-81425c2-43a431a7ed16, np 1;
// tools/p4_scdrv_fixture_config_v1.py + p4_scdrv_runs_v1.sh of lane p4), as written in its
// surface-Q-corrected.csv.
namespace
{

constexpr int kSweepFrequencies = 4;
// Legacy form against P: MEASURED <= 4.0e-13 relative (np 1 and np 2, this build and P's
// own np 2 run); the bound is 250 x that.
constexpr double kLegacyVersusMainRelativeTolerance = 1.0e-10;
// Default form (the island sheet uncracked) against the legacy references: MEASURED
// 3.2e-5 (np 1) and 1.3e-4 (np 2) relative; the bound is the owner's line for a difference
// attributable to cracking on a corrected participation (decision 1566 (5)).
constexpr double kDefaultVersusLegacyRelativeTolerance = 1.0e-3;

// The references hold at 1 and 2 ranks only. On more ranks the corrected values of this
// fixture depend on the MPI partition while the raw field does not (MEASURED against np 1,
// on our pre-merge main as well: np 3 1.1e-2, np 4 7.8e-3, np 6 2.0e-2 relative; review
// P4-batch2-r2109 M-1 / M-2; open item OI-CORRECTED-VALUES-PARTITION-DEPENDENT-C3-FIXTURE).
constexpr int kMaxRanksWithMeasuredReferences = 2;

}  // namespace

TEST_CASE_METHOD(
    test::SurfaceResponseFiles,
    "Self-consistent driven island runs through its last frequency",
    "[drivensolver][surfaceresponseoperator][correctedsolve][Serial][Parallel]")
{
#if !defined(MFEM_USE_GSLIB)
  SKIP("SurfaceResponseOperator requires MFEM_USE_GSLIB");
#else
  MPI_Comm comm = Mpi::World();
  if (Mpi::Size(comm) > kMaxRanksWithMeasuredReferences)
  {
    SKIP(
        "References measured at 1-2 ranks only: on more ranks the corrected values of this "
        "fixture depend on the MPI partition (up to 2e-2 relative at 3-6 ranks; the raw "
        "field does not); open item OI-CORRECTED-VALUES-PARTITION-DEPENDENT-C3-FIXTURE");
  }
  const auto crack_policy = GENERATE(test::CrackPolicy::Legacy, test::CrackPolicy::Default);
  INFO("crack policy: " << test::CrackPolicyName(crack_policy));
  const fs::path mesh_path = temp.temp_dir / "self-consistent-driven-island.mesh";
  if (Mpi::Root(comm))
  {
    mfem::Mesh serial = MakeSerialIslandMesh();
    std::ofstream output(mesh_path);
    serial.Print(output);
  }
  Mpi::Barrier(comm);
  json config = ConvexMaxwellIslandConfig(crack_policy);
  const fs::path output = temp.temp_dir / "self-consistent-driven";
  config["Problem"] = {{"Type", "Driven"}, {"Verbose", 0}, {"Output", output.string()}};
  config["Model"]["Mesh"] = mesh_path.string();
  // Mesh units are metres, unscaled, as the library's R (m).
  config["Model"]["L0"] = 1.0;
  config["Model"]["Lc"] = 1.0;
  config["Boundaries"]["Ground"]["Attributes"] = {2, 3, 4, 5, 6, 9};
  config["Boundaries"]["LumpedPort"] = {{{"Index", 1},
                                         {"Attributes", {1}},
                                         {"R", 50.0},
                                         {"Direction", "+X"},
                                         {"Excitation", 1}}};
  config["Solver"]["Driven"] = {
      {"Samples",
       {{{"Type", "Linear"}, {"MinFreq", 0.02}, {"MaxFreq", 0.08}, {"FreqStep", 0.02}}}}};
  config["Solver"]["Linear"] = {
      {"Type", "Default"}, {"KSPType", "GMRES"}, {"Tol", 1.0e-12}, {"MaxIts", 400}};
  auto &correction = config["Solver"]["SurfaceResponseCorrection"];
  correction["Developer"] = true;
  correction["CorrectionMode"] = "Both";
  correction["SolveTol"] = 1.0e-12;
  {
    IoData iodata(config, false);
    MakeOutputFolder(iodata, comm);
    BlockTimer::Reset();
    palace::Run(iodata, comm, utils::ConfigureOmp(), /*git_tag=*/nullptr);
    Mpi::Barrier(comm);
  }
  if (!Mpi::Root(comm))
  {
    return;
  }
  // Every frequency ran its corrected solve (the totals are saved when the driver releases
  // the corrected solver after the last one).
  std::ifstream metadata(output / "palace.json");
  REQUIRE(metadata);
  const json run_metadata = json::parse(metadata);
  CHECK(run_metadata.at("SurfaceResponse")
            .at("CorrectedSolver")
            .at("Driven")
            .at("LinearSolves") == kSweepFrequencies);
  const Table table = test::LoadCsv(output / "surface-Q-corrected.csv");
  const auto &frequencies = test::ColumnByHeader(table, "f (GHz)").data;
  REQUIRE(frequencies.size() == kSweepFrequencies);
  CHECK_THAT(frequencies.back(), WithinRel(0.08, 1.0e-12));
  const std::vector<std::pair<std::string, std::array<double, kSweepFrequencies>>>
      references = {
          {"p_surf corrected",
           {0.03251911452622, 0.034946641213, 0.039589317411509999, 0.047429695007490003}},
          {"p_surf corrected fixed-flux",
           {0.0025214304518989998, 0.0027147782796970001, 0.0030874838283299999,
            0.0037250888080689998}},
          {"p_surf self-consistent",
           {0.025449240048319999, 0.0275616578931, 0.031652251354559999,
            0.038698485224579998}},
          {"E_norm corrected (J)",
           {6.3069056867340005e-11, 1.11496725789e-10, 1.4046383290580001e-10,
            1.7435415302449999e-10}},
          {"E_norm corrected fixed-flux (J)",
           {6.1172348804359997e-11, 1.0789319581579999e-10, 1.353207663103e-10,
            1.6670558180939999e-10}},
          {"E_norm self-consistent (J)",
           {6.2962955859730003e-11, 1.11673557536e-10, 1.4171814555580001e-10,
            1.7887717361100001e-10}},
          {"E_surf corrected (J)",
           {2.0509498833290001e-12, 3.8964360725720004e-12, 5.5608672657449996e-12,
            8.2695643012440003e-12}},
          {"E_surf corrected fixed-flux (J)",
           {1.5424182308950001e-13, 2.9290610452790002e-13, 4.178006776203e-13,
            6.20993097041e-13}},
          {"E_surf self-consistent (J)",
           {1.6023593778260001e-12, 3.0779083885120002e-12, 4.485698364634e-12,
            6.9222756600010003e-12}}};
  const double tolerance = (crack_policy == test::CrackPolicy::Legacy)
                               ? kLegacyVersusMainRelativeTolerance
                               : kDefaultVersusLegacyRelativeTolerance;
  for (const auto &[header, reference] : references)
  {
    const auto &values = test::ColumnByHeader(table, header).data;
    REQUIRE(values.size() == kSweepFrequencies);
    for (int i = 0; i < kSweepFrequencies; i++)
    {
      INFO(header << " at f = " << frequencies[i] << " GHz");
      CHECK_THAT(values[i], WithinRel(reference[i], tolerance));
    }
  }
#endif
}

// The capacitor energy of the self-consistent normalization (review of the sync with
// Palace 0.20.0, physics MAJOR-1): E_norm self-consistent = E_elec + the time-averaged
// stored energy 1/4 |C| |V|^2 of every lumped capacitor at the corrected field's port
// voltage + the domain correction, as the raw normalization uses Palace's E_cap. A
// complex field, given as both the raw and the corrected field, on the island of the
// self-consistent driven case with a lumped port with C > 0: the self-consistent and the
// postprocessed normalizations of the same field are equal, and their capacitor part is
// Palace's E_cap = 1/4 |C| |V|^2 of the printed port voltage (the former literal
// 1/2 |C| |V|^2 counted it twice).
TEST_CASE_METHOD(
    test::SurfaceResponseFiles,
    "Self-consistent normalization carries the time-averaged capacitor energy",
    "[postoperator][surfaceresponseoperator][correctedsolve][Serial][Parallel]")
{
#if !defined(MFEM_USE_GSLIB)
  SKIP("SurfaceResponseOperator requires MFEM_USE_GSLIB");
#else
  MPI_Comm comm = Mpi::World();
  const fs::path mesh_path = temp.temp_dir / "self-consistent-capacitor-island.mesh";
  if (Mpi::Root(comm))
  {
    mfem::Mesh serial = MakeSerialIslandMesh();
    std::ofstream output(mesh_path);
    serial.Print(output);
  }
  Mpi::Barrier(comm);
  constexpr double port_C = 1.0e-11;  // F
  json config = ConvexMaxwellIslandConfig(test::CrackPolicy::Default);
  const fs::path output = temp.temp_dir / "self-consistent-capacitor";
  config["Problem"] = {{"Type", "Driven"}, {"Verbose", 0}, {"Output", output.string()}};
  config["Model"]["Mesh"] = mesh_path.string();
  config["Model"]["L0"] = 1.0;
  config["Model"]["Lc"] = 1.0;
  config["Boundaries"]["Ground"]["Attributes"] = {2, 3, 4, 5, 6, 9};
  config["Boundaries"]["LumpedPort"] = {{{"Index", 1},
                                         {"Attributes", {1}},
                                         {"R", 50.0},
                                         {"C", port_C},
                                         {"Direction", "+X"},
                                         {"Excitation", 1}}};
  config["Solver"]["Driven"] = {
      {"Samples",
       {{{"Type", "Linear"}, {"MinFreq", 0.02}, {"MaxFreq", 0.08}, {"FreqStep", 0.02}}}}};
  auto &correction = config["Solver"]["SurfaceResponseCorrection"];
  correction["Developer"] = true;
  correction["CorrectionMode"] = "Both";
  {
    IoData iodata(config, false);
    std::vector<std::unique_ptr<Mesh>> mesh;
    {
      auto smesh = mesh::Load(iodata, comm);
      int ne = smesh ? smesh->GetNE() : 0;
      Mpi::GlobalMax(1, &ne, comm);
      if (Mpi::Size(comm) > ne)
      {
        SKIP("Mesh has fewer elements than processes");
      }
      iodata.NondimensionalizeInputs(smesh);
      std::vector<std::unique_ptr<mfem::ParMesh>> mfem_mesh;
      mfem_mesh.push_back(mesh::Partition(iodata, std::move(smesh), comm));
      mesh::RefineMesh(iodata, mfem_mesh);
      for (auto &m : mfem_mesh)
      {
        mesh.push_back(std::make_unique<Mesh>(std::move(m)));
      }
    }
    MakeOutputFolder(iodata, comm);
    SpaceOperator space_op(iodata, mesh);
    PostOperator<ProblemType::DRIVEN> post_op(iodata, space_op);
    // A complex field: the gradients of two smooth potentials as its real and imaginary
    // parts, zero on the perfect conductors (their tangential true dofs).
    auto GradientField = [&](double (*potential)(const mfem::Vector &), Vector &out)
    {
      mfem::ParGridFunction phi(&space_op.GetH1Space().Get());
      mfem::FunctionCoefficient phi_coefficient(potential);
      phi.ProjectCoefficient(phi_coefficient);
      mfem::ParDiscreteLinearOperator gradient(&space_op.GetH1Space().Get(),
                                               &space_op.GetNDSpace().Get());
      gradient.AddDomainInterpolator(new mfem::GradientInterpolator());
      gradient.Assemble();
      mfem::ParGridFunction field(&space_op.GetNDSpace().Get());
      gradient.Mult(phi, field);
      field.GetTrueDofs(out);
      out.SetSubVector(space_op.GetNDDbcTDofLists().back(), 0.0);
    };
    ComplexVector e(space_op.GetNDSpace().GetTrueVSize()),
        b(space_op.GetRTSpace().GetTrueVSize());
    GradientField([](const mfem::Vector &x)
                  { return std::sin(1.3 * x[0] + 0.4) * std::cosh(0.7 * x[1]); }, e.Real());
    GradientField([](const mfem::Vector &x)
                  { return 0.6 * std::cos(0.9 * x[0] + 0.2 * x[2]); }, e.Imag());
    b = 0.0;
    const std::complex<double> omega(iodata.solver.driven.sample_f.back(), 0.0);
    post_op.SetSurfaceResponseCorrectedField(e, nullptr, omega, 1.0);
    post_op.MeasureAndPrintAll(1, 0, e, b, omega);
  }
  Mpi::Barrier(comm);
  if (!Mpi::Root(comm))
  {
    return;
  }
  const Table surface_q = test::LoadCsv(output / "surface-Q-corrected.csv");
  const Table domain_e = test::LoadCsv(output / "domain-E.csv");
  const Table port_v = test::LoadCsv(output / "port-V.csv");
  const double E_norm_postprocessed =
      test::ColumnByHeader(surface_q, "E_norm corrected (J)").data.at(0);
  const double E_norm_self_consistent =
      test::ColumnByHeader(surface_q, "E_norm self-consistent (J)").data.at(0);
  const double E_cap = test::ColumnByHeader(domain_e, "E_cap (J)").data.at(0);
  const std::complex<double> V(test::ColumnByHeader(port_v, "Re{V[1]} (V)").data.at(0),
                               test::ColumnByHeader(port_v, "Im{V[1]} (V)").data.at(0));
  INFO("E_norm corrected " << E_norm_postprocessed << " J, E_norm self-consistent "
                           << E_norm_self_consistent << " J, E_cap " << E_cap
                           << " J, V = " << V.real() << " + " << V.imag() << "i V");
  // The case has power: the field is complex and the capacitor carries a part of the
  // normalization well above the print precision of the tables.
  REQUIRE(std::abs(V.real()) > 0.0);
  REQUIRE(std::abs(V.imag()) > 0.0);
  REQUIRE(E_cap > 1.0e-2 * E_norm_postprocessed);
  // Palace's E_cap is the time average 1/4 |C| |V|^2 of the peak-phasor port voltage.
  CHECK_THAT(E_cap, WithinRel(0.25 * port_C * std::norm(V), 1.0e-9));
  // The capacitor part of the self-consistent normalization is Palace's E_cap: the two
  // normalizations of the same field are equal.
  const double capacitor_part = E_norm_self_consistent - (E_norm_postprocessed - E_cap);
  CHECK_THAT(capacitor_part, WithinRel(E_cap, 1.0e-8));
  CHECK_THAT(E_norm_self_consistent, WithinRel(E_norm_postprocessed, 1.0e-10));
#endif
}

}  // namespace palace
