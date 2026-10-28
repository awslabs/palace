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
#include <iterator>
#include <limits>
#include <string>
#include <vector>
#include <Eigen/Dense>
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

// A 1 x 0.8 x 0.6 cm box of 2 x 2 x 2 hexahedra (L0 = 1 cm), without degenerate modes.
// Boundary attributes of mfem::Mesh::MakeCartesian3D: 1 (z = 0), 2 (y = 0), 3 (x = 1), 4 (y
// = 1), 5 (x = 0), 6 (z = 1).
fs::path WriteBoxMesh(const fs::path &dir)
{
  const auto path = dir / "box.mesh";
  if (Mpi::Root(Mpi::World()))
  {
    auto mesh =
        mfem::Mesh::MakeCartesian3D(2, 2, 2, mfem::Element::HEXAHEDRON, 1.0, 0.8, 0.6);
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
        {"Linear", {{"Type", "Default"}, {"Tol", 1.0e-12}, {"MaxIts", 1000}}},
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
  // eigensolver (Quasi-Newton after the Hybrid seed) must satisfy this for some f_h.
  const auto mesh = WriteBoxMesh(temp_dir);
  const json pec = {{"PEC", {{"Attributes", {1, 2, 3, 4, 5, 6}}}}};
  const auto constant =
      RunEigenmode(BoxConfig(mesh, eps_inf, pec, 16, 8.0, ""), temp_dir, "constant");
  const auto lorentz = RunEigenmode(
      BoxConfig(mesh, LorentzPermittivityConfig(), pec, 6, 8.0, ""), temp_dir, "lorentz");
  REQUIRE(constant.f.size() >= 16);
  REQUIRE(lorentz.f.size() >= 6);
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
  }

  // The roots above the target are those of the quartic f² (ε∞ D(f) + Δε f₀²) =
  // f_h² ε∞ D(f), D(f) = f₀² − f² + i γ f, for every f_h. With the target below the
  // resonance and the constant modes above it, the requested modes are roots between the
  // target and the pole, on the other side of the pole from the constant modes, distinct,
  // and include the root nearest the target. The roots accumulate at the pole, where the
  // polynomial seed of the Hybrid method resolves them only approximately, so the modes
  // found need not be the next ones in order after the first.
  std::vector<std::complex<double>> expected;
  const std::complex<double> i(0.0, 1.0);
  for (const auto f_h : constant.f)
  {
    const double h2 = std::norm(f_h) * eps_inf;
    // Monic quartic f⁴ + c₃ f³ + c₂ f² + c₁ f + c₀ (divided by −ε∞), as a companion matrix.
    const std::complex<double> c3 = -i * damping,
                               c2 = -(eps_inf * f0 * f0 + delta_eps * f0 * f0 + h2) /
                                    eps_inf,
                               c1 = i * damping * h2 / eps_inf, c0 = h2 * f0 * f0 / eps_inf;
    Eigen::Matrix4cd companion = Eigen::Matrix4cd::Zero();
    companion(0, 0) = -c3;
    companion(0, 1) = -c2;
    companion(0, 2) = -c1;
    companion(0, 3) = -c0;
    companion(1, 0) = companion(2, 1) = companion(3, 2) = 1.0;
    // A named solver: the eigenvalues are a reference into it, and in C++20 a temporary in
    // the range of a range-based for loop would be destroyed before the loop body.
    const Eigen::ComplexEigenSolver<Eigen::Matrix4cd> roots(companion);
    REQUIRE(roots.info() == Eigen::Success);
    for (const auto f : roots.eigenvalues())
    {
      if (f.real() >= 8.0 && f.real() < f0)
      {
        expected.push_back(f);
      }
    }
  }
  auto by_frequency = [](auto a, auto b) { return a.real() < b.real(); };
  std::ranges::sort(expected, by_frequency);
  REQUIRE(expected.size() >= 6);
  std::vector<std::complex<double>> found;
  std::ranges::copy_if(lorentz.f, std::back_inserter(found),
                       [](auto f) { return f.real() >= 8.0; });
  std::ranges::sort(found, by_frequency);
  REQUIRE(found.size() >= 6);
  auto Match = [&](std::complex<double> f)
  {
    for (std::size_t k = 0; k < expected.size(); k++)
    {
      if (std::abs(f - expected[k]) <= 1.0e-7 * std::abs(expected[k]))
      {
        return static_cast<int>(k);
      }
    }
    return -1;
  };
  std::vector<int> matched;
  for (const auto f : found)
  {
    CAPTURE(f);
    const int k = Match(f);
    CHECK(k >= 0);  // A root on the target's side of the resonance
    CHECK(std::ranges::find(matched, k) == matched.end());
    matched.push_back(k);
  }
  CHECK(matched.front() == 0);
}

TEST_CASE_METHOD(
    test::SharedTempDir,
    "Eigenmode driver: second-order absorbing boundaries next to a Debye medium",
    "[eigenmodedriver][Serial][Parallel]")
{
#if !defined(PALACE_WITH_SLEPC)
  SKIP("Requires the SLEPc eigensolver");
#endif
  // The faces x = 0 and x = 1 are second-order absorbing boundaries next to a Debye medium
  // (ε∞ = 2, Δε = 0.5, τ = 0.01 ns), and the others PEC. The Hybrid nonlinear eigensolver
  // freezes the pole terms of the absorbing boundaries in its polynomial seed
  // (SpaceOperator::GetFarfieldPoleTerms), including those next to the frequency-dependent
  // material. Each mode found must be an eigenpair (backward error), and they must be the
  // four modes of lowest frequency above the target, two of them strongly damped. These
  // reference values were computed with this solver and confirmed by their backward
  // errors; without the frozen terms next to the material, the seed missed the two damped
  // modes and Quasi-Newton returned two others instead, (23.035 + 2.534i) and
  // (23.291 + 2.215i) GHz. SLP, as an independent solver, did not finish in minutes here.
  const auto mesh = WriteBoxMesh(temp_dir);
  const json boundaries = {{"PEC", {{"Attributes", {1, 2, 4, 6}}}},
                           {"Absorbing", {{"Attributes", {3, 5}}, {"Order", 2}}}};
  const json debye = {
      {"HighFrequency", eps_inf},
      {"Terms",
       {{{"Type", "Debye"}, {"DeltaPermittivity", 0.5}, {"RelaxationTime", 0.01}}}}};
  constexpr int n = 4;
  const auto hybrid = RunEigenmode(BoxConfig(mesh, debye, boundaries, n, 18.0, "Hybrid"),
                                   temp_dir, "hybrid");
  REQUIRE(hybrid.f.size() >= n);
  const std::vector<std::complex<double>> expected = {{19.108678, 2.399168},
                                                      {22.125867, 26.837326},
                                                      {22.287635, 6.409018},
                                                      {22.567106, 25.977182}};
  auto found = hybrid.f;
  std::ranges::sort(found, [](auto a, auto b) { return a.real() < b.real(); });
  for (std::size_t i = 0; i < hybrid.f.size(); i++)
  {
    CAPTURE(hybrid.f[i]);
    CHECK(hybrid.backward_error[i] <= 1.0e-8);
    CHECK(hybrid.f[i].imag() > 0.0);  // Radiating through the absorbing faces
  }
  for (std::size_t k = 0; k < n; k++)
  {
    CAPTURE(k, found[k], expected[k]);
    CHECK(std::abs(found[k] - expected[k]) <= 1.0e-6 * std::abs(expected[k]));
  }
}

}  // namespace palace
