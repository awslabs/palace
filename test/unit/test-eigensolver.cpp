// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include <cmath>
#include <complex>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <iterator>
#include <optional>
#include <string>
#include <utility>
#include <vector>
#include <unistd.h>
#include <fmt/format.h>
#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>
#include "driver.hpp"
#include "drivers/eigensolver.hpp"
#include "surfaceresponse-fixtures.hpp"
#include "utils/communication.hpp"
#include "utils/filesystem.hpp"
#include "utils/iodata.hpp"
#include "utils/omp.hpp"
#include "utils/outputdir.hpp"
#include "utils/tablecsv.hpp"
#include "utils/timer.hpp"

using palace::BuildComplexSumOperator;
using palace::ComplexWrapperOperator;
using palace::EigenSolverBackend;
using palace::EigenvalueSolver;
using palace::NonlinearEigenSolver;
using palace::internal::BuildResponseCorrectedMass;
using palace::internal::BuildResponseCorrectedStiffnesses;
using palace::internal::EigenvalueToAngularFrequency;
using palace::internal::GetEigenvalueTarget;
using palace::internal::MaximumWeightModeAssignment;
using palace::internal::VerifySelfConsistentMagneticNonlinearSolver;

TEST_CASE("Eigenmode response-corrected mass", "[eigensolver][surface-response][Serial]")
{
  mfem::DenseMatrix mass_real(2), mass_imag(2), response(2);
  mass_real = 0.0;
  mass_imag = 0.0;
  response = 0.0;
  mass_real(0, 0) = 2.0;
  mass_real(1, 1) = 3.0;
  mass_imag(0, 1) = 0.5;
  mass_imag(1, 0) = 0.5;
  response(0, 0) = -0.25;
  response(1, 1) = 0.75;

  ComplexWrapperOperator mass(&mass_real, &mass_imag);
  auto corrected = BuildResponseCorrectedMass(mass, response);
  REQUIRE(corrected.real);
  REQUIRE(corrected.op);
  CHECK(corrected.op->Real() == corrected.real.get());
  CHECK(corrected.op->Imag() == &mass_imag);

  mfem::Vector x({2.0, -1.0}), y(2);
  corrected.op->Real()->Mult(x, y);
  CHECK(y[0] == Catch::Approx(3.5));
  CHECK(y[1] == Catch::Approx(-3.75));

  auto sum = BuildComplexSumOperator({{2.0, corrected.op.get()}});
  palace::ComplexVector xc(2), yc(2);
  xc.Real() = x;
  xc.Imag() = 0.0;
  sum->Mult(xc, yc);
  CHECK(yc.Real()[0] == Catch::Approx(7.0));
  CHECK(yc.Real()[1] == Catch::Approx(-7.5));
  CHECK(yc.Imag()[0] == Catch::Approx(-1.0));
  CHECK(yc.Imag()[1] == Catch::Approx(2.0));
}

// OI-UPSYNC-NULL-KP (P3-11): with magnetic leaves the corrected eigenproblem builds K + dK
// and, only when the HYBRID polynomial seed Kp = K + A2_0 exists, Kp + dK. On the SLP path
// (an A2 term, no seed) Kp is null: the corrected Kp was built from it.
TEST_CASE("Eigenmode response-corrected stiffness without a polynomial seed",
          "[eigensolver][surface-response][Serial]")
{
  mfem::DenseMatrix stiffness_real(2), seed_real(2), seed_imag(2), response(2);
  stiffness_real = 0.0;
  seed_real = 0.0;
  seed_imag = 0.0;
  response = 0.0;
  stiffness_real(0, 0) = 2.0;
  stiffness_real(1, 1) = 3.0;
  seed_real(0, 0) = 5.0;
  seed_real(1, 1) = 7.0;
  seed_imag(0, 1) = 0.5;
  seed_imag(1, 0) = 0.5;
  response(0, 0) = -0.25;
  response(1, 1) = 0.75;
  ComplexWrapperOperator stiffness(&stiffness_real, nullptr);
  ComplexWrapperOperator seed(&seed_real, &seed_imag);
  mfem::Vector x({2.0, -1.0}), y(2);

  // SLP: no seed, so only K + dK.
  auto slp = BuildResponseCorrectedStiffnesses(stiffness, nullptr, response);
  REQUIRE(slp.K.op);
  CHECK_FALSE(slp.Kp.op);
  CHECK_FALSE(slp.Kp.real);
  slp.K.op->Real()->Mult(x, y);
  CHECK(y[0] == Catch::Approx(3.5));
  CHECK(y[1] == Catch::Approx(-3.75));

  // HYBRID: both twins; Kp + dK keeps the seed's imaginary part.
  auto hybrid = BuildResponseCorrectedStiffnesses(stiffness, &seed, response);
  REQUIRE(hybrid.K.op);
  REQUIRE(hybrid.Kp.op);
  CHECK(hybrid.Kp.op->Imag() == &seed_imag);
  hybrid.Kp.op->Real()->Mult(x, y);
  CHECK(y[0] == Catch::Approx(9.5));
  CHECK(y[1] == Catch::Approx(-7.75));
}

TEST_CASE("Eigenmode self-consistent magnetic defect refused under SLP",
          "[eigensolver][surface-response][Serial]")
{
  // Electric-only libraries and the HYBRID seed keep the self-consistent corrected solve.
  CHECK_NOTHROW(
      VerifySelfConsistentMagneticNonlinearSolver(false, NonlinearEigenSolver::SLP));
  CHECK_NOTHROW(
      VerifySelfConsistentMagneticNonlinearSolver(false, NonlinearEigenSolver::HYBRID));
  CHECK_NOTHROW(
      VerifySelfConsistentMagneticNonlinearSolver(true, NonlinearEigenSolver::HYBRID));

  // K + dK beside the exact A2 of the SLP solver is refused by name.
  CHECK_THROWS_WITH(
      VerifySelfConsistentMagneticNonlinearSolver(true, NonlinearEigenSolver::SLP),
      Catch::Matchers::ContainsSubstring("NonlinearType \"SLP\"") &&
          Catch::Matchers::ContainsSubstring("magnetic (inductive) defect") &&
          Catch::Matchers::ContainsSubstring("OI-EIGEN-SLP-MAGNETIC-SELF-CONSISTENT"));
}

TEST_CASE("Eigenmode spectral conventions", "[eigensolver][surface-response][Serial]")
{
  using Which = EigenvalueSolver::WhichType;

  SECTION("Linear generalized eigenproblem")
  {
    auto target = GetEigenvalueTarget(3.0, false, EigenSolverBackend::SLEPC,
                                      NonlinearEigenSolver::HYBRID);
    CHECK(target.shift == std::complex<double>(9.0, 0.0));
    CHECK(target.which == Which::TARGET_REAL);

    target = GetEigenvalueTarget(3.0, false, EigenSolverBackend::ARPACK,
                                 NonlinearEigenSolver::HYBRID);
    CHECK(target.which == Which::LARGEST_REAL);
    CHECK(EigenvalueToAngularFrequency({9.0, 0.0}, false) ==
          std::complex<double>(3.0, 0.0));
  }

  SECTION("Polynomial and nonlinear eigenproblems")
  {
    auto target = GetEigenvalueTarget(3.0, true, EigenSolverBackend::SLEPC,
                                      NonlinearEigenSolver::HYBRID);
    CHECK(target.shift == std::complex<double>(0.0, 3.0));
    CHECK(target.which == Which::TARGET_IMAGINARY);

    target = GetEigenvalueTarget(3.0, true, EigenSolverBackend::SLEPC,
                                 NonlinearEigenSolver::SLP);
    CHECK(target.which == Which::TARGET_MAGNITUDE);

    target = GetEigenvalueTarget(3.0, true, EigenSolverBackend::ARPACK,
                                 NonlinearEigenSolver::HYBRID);
    CHECK(target.which == Which::SMALLEST_IMAGINARY);

    const auto omega = EigenvalueToAngularFrequency({-0.25, 3.0}, true);
    CHECK(omega.real() == Catch::Approx(3.0));
    CHECK(omega.imag() == Catch::Approx(0.25));
  }
}

TEST_CASE("Eigenmode maximum-weight assignment", "[eigensolver][surface-response][Serial]")
{
  SECTION("Empty")
  {
    CHECK(MaximumWeightModeAssignment({}).empty());
  }

  SECTION("Global optimum differs from row-greedy matching")
  {
    const auto assignment = MaximumWeightModeAssignment({{0.90, 0.80}, {0.85, 0.10}});
    REQUIRE(assignment.size() == 2);
    CHECK(assignment[0] == 1);
    CHECK(assignment[1] == 0);
  }

  SECTION("Rectangular")
  {
    const auto assignment =
        MaximumWeightModeAssignment({{0.20, 0.95, 0.10}, {0.80, 0.10, 0.70}});
    REQUIRE(assignment.size() == 2);
    CHECK(assignment[0] == 1);
    CHECK(assignment[1] == 0);
  }

  SECTION("Invalid")
  {
    CHECK_THROWS(MaximumWeightModeAssignment({{1.0}, {0.5}}));
    CHECK_THROWS(MaximumWeightModeAssignment({{1.0, 0.0}, {0.5}}));
    CHECK_THROWS(MaximumWeightModeAssignment({{1.0, -0.1}}));
  }
}

namespace palace
{

// SD-507 (eigen #8) / SD-503 (P3-09): with field output the eigen driver stages the raw
// and the assigned self-consistent corrected eigenvectors in host memory and releases both
// eigensolvers before postprocessing (upstream #825); without field output the
// postprocessing reads them from the solvers. Both forms must write bitwise the same
// tables, and the first written mode must carry the corrected field (registered on the
// lazily created ParaView collection). Two modes, CorrectionMode Both, on the convex PEC
// island of the fixtures (written to a mesh file), in both crack forms.
TEST_CASE_METHOD(test::SurfaceResponseFiles,
                 "Eigenmode corrected postprocessing staged = unstaged",
                 "[eigensolver][surface-response][Serial][Parallel]")
{
#if !defined(MFEM_USE_GSLIB)
  SKIP("SurfaceResponseOperator requires MFEM_USE_GSLIB");
#else
  // Legacy: the island sheet (attribute 9) is cracked; default (the key unset): it is not.
  const auto crack_policy = GENERATE(test::CrackPolicy::Legacy, test::CrackPolicy::Default);
  const bool legacy = (crack_policy == test::CrackPolicy::Legacy);
  INFO("crack policy: " << test::CrackPolicyName(crack_policy));
  MPI_Comm comm = Mpi::World();
  const fs::path mesh_path = temp.temp_dir / "staged-island.mesh";
  if (Mpi::Root(comm))
  {
    mfem::Mesh serial = MakeSerialIslandMesh();
    std::ofstream output(mesh_path);
    serial.Print(output);
  }
  Mpi::Barrier(comm);
  // The builder writes the crack key of the policy (SD-3301: legacy = true, default =
  // unset).
  test::json config = ConvexMaxwellIslandConfig(crack_policy);
  config["Model"]["Mesh"] = mesh_path.string();
  // Mesh units are metres, unscaled, as the library's R (m).
  config["Model"]["L0"] = 1.0;
  config["Model"]["Lc"] = 1.0;
  config["Solver"]["SurfaceResponseCorrection"]["CorrectionMode"] = "Both";
  // The two lowest modes of the 2 m PEC box (about 0.1 GHz).
  config["Solver"]["Eigenmode"] = {{"Target", 0.05}, {"N", 2}};
  auto Run = [&](const std::string &stem, int save)
  {
    test::json run = config;
    const fs::path output = temp.temp_dir / (stem + (legacy ? "-legacy" : "-default"));
    run["Problem"]["Output"] = output.string();
    run["Solver"]["Eigenmode"]["Save"] = save;
    IoData iodata(run, false);
    MakeOutputFolder(iodata, comm);
    BlockTimer::Reset();
    palace::Run(iodata, comm, utils::ConfigureOmp(), /*git_tag=*/nullptr);
    Mpi::Barrier(comm);
    return output;
  };
  const fs::path unstaged = Run("unstaged", 0);
  const fs::path staged = Run("staged", 2);
  if (!Mpi::Root(comm))
  {
    return;
  }

  // Both modes carry a corrected and a self-consistent surface participation.
  const Table table = test::LoadCsv(unstaged / "surface-Q-corrected.csv");
  const auto &modes = test::ColumnByHeader(table, "m").data;
  const auto &corrected = test::ColumnByHeader(table, "p_surf corrected").data;
  const auto &self_consistent = test::ColumnByHeader(table, "p_surf self-consistent").data;
  REQUIRE(modes.size() == 2);
  for (std::size_t i = 0; i < modes.size(); i++)
  {
    INFO("mode " << modes[i]);
    CHECK(std::isfinite(corrected[i]));
    CHECK(corrected[i] > 0.0);
    CHECK(std::isfinite(self_consistent[i]));
    CHECK(self_consistent[i] > 0.0);
  }

  // Every table is bitwise the same with and without the staging.
  std::vector<std::string> names;
  for (const auto &entry : fs::directory_iterator(unstaged))
  {
    if (entry.is_regular_file() && entry.path().extension() == ".csv")
    {
      names.push_back(entry.path().filename().string());
    }
  }
  std::size_t staged_count = 0;
  for (const auto &entry : fs::directory_iterator(staged))
  {
    staged_count += (entry.is_regular_file() && entry.path().extension() == ".csv");
  }
  REQUIRE(names.size() > 1);
  CHECK(staged_count == names.size());
  for (const auto &name : names)
  {
    INFO("table " << name);
    REQUIRE(fs::is_regular_file(staged / name));
    CHECK(test::ReadFile(staged / name) == test::ReadFile(unstaged / name));
  }

  // The staged run wrote fields, and its first written mode holds the corrected field.
  CHECK_FALSE(fs::exists(unstaged / "paraview"));
  REQUIRE(fs::is_directory(staged / "paraview"));
  fs::path first_cycle;
  for (const auto &entry : fs::recursive_directory_iterator(staged / "paraview"))
  {
    const auto name = entry.path().filename().string();
    const auto collection = entry.path().parent_path().filename().string();
    const bool boundary = collection.size() >= 9 &&
                          collection.compare(collection.size() - 9, 9, "_boundary") == 0;
    if (entry.is_directory() && name.rfind("Cycle", 0) == 0 && !boundary &&
        (first_cycle.empty() || entry.path() < first_cycle))
    {
      first_cycle = entry.path();
    }
  }
  REQUIRE_FALSE(first_cycle.empty());
  INFO("first cycle " << first_cycle.string());
  bool has_corrected_field = false;
  for (const auto &entry : fs::directory_iterator(first_cycle))
  {
    has_corrected_field =
        has_corrected_field ||
        (entry.is_regular_file() &&
         test::ReadFile(entry.path()).find("E_corrected_real") != std::string::npos);
  }
  CHECK(has_corrected_field);
#endif
}

// Decision 2159 (OI-EIGEN-SC-CORRECTED-LINEAR-SOLVE-NONCONVERGED-RESULT-UNFLAGGED): the
// corrected linear solves of the self-consistent eigenproblem that did not reach their
// tolerance are counted in palace.json (SurfaceResponse.CorrectedSolver.Eigenmode.
// NonConvergedLinearSolves, the total over the solve: its shift-invert solves build one
// subspace for all modes) and named by one line at the end of the log when there are any.
// The convex PEC island (CorrectionMode Both, two modes): with the default linear solver
// (a sparse direct preconditioner of the raw operator) every corrected solve converges;
// with Solver.Linear.MaxIts 2 the corrected solves (the raw preconditioner on the corrected
// operator) stop at the limit.
TEST_CASE_METHOD(test::SurfaceResponseFiles,
                 "Eigenmode self-consistent corrected solves that did not converge are "
                 "counted",
                 "[eigensolver][surface-response][Serial][Parallel]")
{
#if !defined(MFEM_USE_GSLIB)
  SKIP("SurfaceResponseOperator requires MFEM_USE_GSLIB");
#else
  MPI_Comm comm = Mpi::World();
  const fs::path mesh_path = temp.temp_dir / "unconverged-island.mesh";
  if (Mpi::Root(comm))
  {
    mfem::Mesh serial = MakeSerialIslandMesh();
    std::ofstream output(mesh_path);
    serial.Print(output);
  }
  Mpi::Barrier(comm);
  test::json config = ConvexMaxwellIslandConfig(test::CrackPolicy::Legacy);
  config["Model"]["Mesh"] = mesh_path.string();
  // Mesh units are metres, unscaled, as the library's R (m).
  config["Model"]["L0"] = 1.0;
  config["Model"]["Lc"] = 1.0;
  config["Solver"]["SurfaceResponseCorrection"]["CorrectionMode"] = "Both";
  // The two lowest modes of the 2 m PEC box (about 0.1 GHz).
  config["Solver"]["Eigenmode"] = {{"Target", 0.05}, {"N", 2}, {"Save", 0}};
  constexpr const char *kUnconvergedLine =
      "corrected linear solves did not reach their tolerance";
  // The run's palace.json and, on the root, its log.
  auto Run = [&](const std::string &stem, std::optional<int> max_its)
  {
    test::json run = config;
    const fs::path output = temp.temp_dir / stem;
    run["Problem"]["Output"] = output.string();
    if (max_its)
    {
      run["Solver"]["Linear"]["MaxIts"] = *max_its;
    }
    IoData iodata(run, false);
    MakeOutputFolder(iodata, comm);
    BlockTimer::Reset();
    std::fflush(stdout);
    const int saved = dup(fileno(stdout));
    std::FILE *log = std::fopen((temp.temp_dir / (stem + ".log")).c_str(), "w");
    // Every rank's flag reduced before the collective run (review r2157 m-3): a rank-local
    // REQUIRE would leave the other ranks alone in palace::Run.
    bool log_opened = (log != nullptr);
    Mpi::GlobalAnd(1, &log_opened, comm);
    REQUIRE(log_opened);
    dup2(fileno(log), fileno(stdout));
    palace::Run(iodata, comm, utils::ConfigureOmp(), /*git_tag=*/nullptr);
    std::fflush(stdout);
    dup2(saved, fileno(stdout));
    close(saved);
    std::fclose(log);
    Mpi::Barrier(comm);
    std::ifstream text(temp.temp_dir / (stem + ".log"));
    std::ifstream metadata(output / "palace.json");
    bool metadata_opened = static_cast<bool>(metadata);
    Mpi::GlobalAnd(1, &metadata_opened, comm);
    REQUIRE(metadata_opened);
    return std::pair{test::json::parse(metadata),
                     std::string(std::istreambuf_iterator<char>(text),
                                 std::istreambuf_iterator<char>())};
  };
  const auto [converged_metadata, converged_log] = Run("converged", std::nullopt);
  const auto [limited_metadata, limited_log] = Run("limited", 2);
  if (!Mpi::Root(comm))
  {
    return;
  }
  const auto &converged =
      converged_metadata.at("SurfaceResponse").at("CorrectedSolver").at("Eigenmode");
  const auto &limited =
      limited_metadata.at("SurfaceResponse").at("CorrectedSolver").at("Eigenmode");
  INFO("default: " << converged.dump() << "; Solver.Linear.MaxIts 2: " << limited.dump());
  CHECK(converged.at("LinearSolves").get<long long>() > 0);
  CHECK(converged.at("NonConvergedLinearSolves").get<long long>() == 0);
  CHECK(converged_log.find(kUnconvergedLine) == std::string::npos);
  const auto limited_count = limited.at("NonConvergedLinearSolves").get<long long>();
  CHECK(limited_count > 0);
  CHECK(limited_count <= limited.at("LinearSolves").get<long long>());
  const auto line =
      fmt::format("{:d} of {:d} {}", limited_count,
                  limited.at("LinearSolves").get<long long>(), kUnconvergedLine);
  INFO("expected log line: " << line);
  CHECK(limited_log.find(line) != std::string::npos);
  CHECK(limited_log.find("Solver.Linear.MaxIts") != std::string::npos);
#endif
}

}  // namespace palace
