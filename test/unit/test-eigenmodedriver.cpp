// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

// End-to-end tests of the eigenmode driver with frequency-dependent materials on small box
// meshes, run in process (unlike the regression cases, so that coverage and sanitizer
// builds include the nonlinear eigensolver paths).

#include <algorithm>
#include <cmath>
#include <complex>
#include <filesystem>
#include <fstream>
#include <limits>
#include <string>
#include <vector>
#include <mfem.hpp>
#include <catch2/catch_test_macros.hpp>
#include <nlohmann/json.hpp>
#include "driver.hpp"
#include "fixtures.hpp"
#include "utils/communication.hpp"
#include "utils/iodata.hpp"
#include "utils/omp.hpp"
#include "utils/outputdir.hpp"
#include "utils/tablecsv.hpp"
#include "utils/timer.hpp"

namespace palace
{

namespace
{

using json = nlohmann::json;
namespace fs = std::filesystem;

// A 1 cm cube of 2 x 2 x 2 hexahedra (L0 = 1 cm). Boundary attributes of
// mfem::Mesh::MakeCartesian3D: 1 (z = 0), 2 (y = 0), 3 (x = 1), 4 (y = 1), 5 (x = 0),
// 6 (z = 1).
fs::path WriteBoxMesh(const fs::path &dir)
{
  const auto path = dir / "box.mesh";
  if (Mpi::Root(Mpi::World()))
  {
    auto mesh =
        mfem::Mesh::MakeCartesian3D(2, 2, 2, mfem::Element::HEXAHEDRON, 1.0, 1.0, 1.0);
    std::ofstream out(path);
    mesh.Print(out);
  }
  Mpi::Barrier(Mpi::World());
  return path;
}

// The complex eigenfrequencies f = Re{f} + i Im{f} (GHz) and backward errors of an
// eigenmode solve of the given configuration, with its output in dir / name.
struct Eigenvalues
{
  std::vector<std::complex<double>> f;
  std::vector<double> backward_error;
};

Eigenvalues RunEigenmode(json config, const fs::path &dir, const std::string &name)
{
  MPI_Comm comm = Mpi::World();
  const auto output = dir / name;
  config["Problem"]["Output"] = output.string();
  IoData iodata(config, false);
  MakeOutputFolder(iodata, comm);
  BlockTimer::Reset();
  Run(iodata, comm, utils::ConfigureOmp());
  Mpi::Barrier(comm);
  TableWithCSVFile eig((output / "eig.csv").string(), /*load_existing_file=*/true);
  // Columns: m, Re{f} (GHz), Im{f} (GHz), Q, Error (Bkwd.), Error (Abs.).
  Eigenvalues out;
  for (std::size_t i = 0; i < eig.table[1].data.size(); i++)
  {
    out.f.emplace_back(eig.table[1].data[i], eig.table[2].data[i]);
    out.backward_error.push_back(eig.table[4].data[i]);
  }
  return out;
}

json BoxConfig(const fs::path &mesh, const json &permittivity, const json &boundaries,
               int n, double target, const std::string &nonlinear_type)
{
  json eigenmode = {{"N", n}, {"Tol", 1.0e-9}, {"Target", target}, {"Type", "SLEPc"}};
  if (!nonlinear_type.empty())
  {
    eigenmode["NonlinearType"] = nonlinear_type;
  }
  return {
      {"Problem", {{"Type", "Eigenmode"}, {"Verbose", 0}}},
      {"Model", {{"Mesh", mesh.string()}, {"L0", 1.0e-2}}},
      {"Domains", {{"Materials", {{{"Attributes", {1}}, {"Permittivity", permittivity}}}}}},
      {"Boundaries", boundaries},
      {"Solver",
       {{"Order", 2},
        {"Device", "CPU"},
        {"Linear", {{"Type", "Default"}}},
        {"Eigenmode", eigenmode}}}};
}

// Lorentz term: χ = Δε f₀² / (f₀² − f² + i γ f), with s = i 2π f and the frequencies in GHz
// (the factors 2π cancel).
constexpr double eps_inf = 2.0, delta_eps = 0.1, f0 = 14.0, damping = 0.5;

std::complex<double> LorentzPermittivity(std::complex<double> f)
{
  const std::complex<double> i(0.0, 1.0);
  return eps_inf + delta_eps * f0 * f0 / (f0 * f0 - f * f + i * damping * f);
}

json LorentzPermittivityConfig()
{
  return {{"HighFrequency", eps_inf},
          {"Terms",
           {{{"Type", "Lorentz"},
             {"DeltaPermittivity", delta_eps},
             {"ResonanceFrequency", f0},
             {"DampingFrequency", damping}}}}};
}

}  // namespace

TEST_CASE_METHOD(test::SharedTempDir,
                 "Eigenmode driver: Lorentz-filled PEC box from the constant modes",
                 "[eigenmodedriver][Serial][Parallel]")
{
#if !defined(PALACE_WITH_SLEPC)
  SKIP("Requires the SLEPc eigensolver");
#endif
  // In a PEC box filled with one material, the discrete modes of ε(ω) on the same mesh are
  // exactly those of a constant ε∞ rescaled by the material: f² ε(f) = f_h² ε∞, with f_h a
  // discrete eigenfrequency of the constant ε∞. Each Lorentz mode found by the nonlinear
  // eigensolver (Quasi-Newton after the Hybrid seed) must satisfy this for some f_h,
  // including the modes below the resonance f₀, where ε(f) differs most from ε∞.
  const auto mesh = WriteBoxMesh(temp_dir);
  const json pec = {{"PEC", {{"Attributes", {1, 2, 3, 4, 5, 6}}}}};
  const auto constant =
      RunEigenmode(BoxConfig(mesh, eps_inf, pec, 16, 8.0, ""), temp_dir, "constant");
  const auto lorentz = RunEigenmode(
      BoxConfig(mesh, LorentzPermittivityConfig(), pec, 6, 8.0, ""), temp_dir, "lorentz");
  REQUIRE(constant.f.size() >= 16);
  REQUIRE(lorentz.f.size() >= 6);
  bool below = false;
  for (std::size_t i = 0; i < lorentz.f.size(); i++)
  {
    const auto f = lorentz.f[i];
    CAPTURE(f, lorentz.backward_error[i]);
    CHECK(lorentz.backward_error[i] <= 1.0e-8);
    CHECK(f.imag() > 0.0);  // Damped
    double rel = std::numeric_limits<double>::infinity();
    for (const auto f_h : constant.f)
    {
      const double lhs = std::norm(f_h) * eps_inf;  // f_h is real for a lossless box
      rel = std::min(rel, std::abs(f * f * LorentzPermittivity(f) - lhs) / lhs);
    }
    CHECK(rel <= 1.0e-7);
    below = below || f.real() < f0;
  }
  CHECK(below);
}

TEST_CASE_METHOD(
    test::SharedTempDir,
    "Eigenmode driver: second-order absorbing boundaries next to a Lorentz medium",
    "[eigenmodedriver][Serial][Parallel]")
{
#if !defined(PALACE_WITH_SLEPC)
  SKIP("Requires the SLEPc eigensolver");
#endif
  // The faces x = 0 and x = 1 are second-order absorbing boundaries next to the Lorentz
  // medium, and the others PEC. The Hybrid nonlinear eigensolver, whose polynomial seed
  // freezes the absorbing boundary pole terms (SpaceOperator::GetFarfieldPoleTerms), must
  // find the requested modes, and each of them must be a mode found independently by the
  // SLP solver, which does not use the seed.
  const auto mesh = WriteBoxMesh(temp_dir);
  const json boundaries = {{"PEC", {{"Attributes", {1, 2, 4, 6}}}},
                           {"Absorbing", {{"Attributes", {3, 5}}, {"Order", 2}}}};
  constexpr int n = 4;
  const auto hybrid = RunEigenmode(
      BoxConfig(mesh, LorentzPermittivityConfig(), boundaries, n, 8.0, "Hybrid"), temp_dir,
      "hybrid");
  const auto slp = RunEigenmode(
      BoxConfig(mesh, LorentzPermittivityConfig(), boundaries, 2 * n, 8.0, "SLP"), temp_dir,
      "slp");
  REQUIRE(hybrid.f.size() >= n);
  REQUIRE(slp.f.size() >= n);
  for (std::size_t i = 0; i < n; i++)
  {
    const auto f = hybrid.f[i];
    CAPTURE(f, hybrid.backward_error[i]);
    CHECK(hybrid.backward_error[i] <= 1.0e-8);
    CHECK(f.imag() > 0.0);  // Radiating through the absorbing faces
    double rel = std::numeric_limits<double>::infinity();
    for (const auto g : slp.f)
    {
      rel = std::min(rel, std::abs(f - g) / std::abs(g));
    }
    CHECK(rel <= 1.0e-6);
  }
}

}  // namespace palace
