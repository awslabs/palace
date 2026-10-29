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

// A 1 x 0.8 x 0.6 cm box of 2 x 2 x 2 hexahedra (L0 = 1 cm).
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

  // The requested modes are the N lowest above the target, counted with multiplicity: the
  // roots of the quartic f² (ε∞ D(f) + Δε f₀²) = f_h² ε∞ D(f), D(f) = f₀² − f² + i γ f, for
  // every constant mode f_h, including both modes of a degenerate f_h. With the target
  // below the resonance and the constant modes above it, the six lowest roots lie between
  // the target and the pole, on the other side of the pole from the constant modes: one of
  // them skipped or replaced by a root beyond the pole fails. (This box has degenerate
  // constant modes, at 24.586, 30.767, and 33.750 GHz.)
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
      if (f.real() >= 8.0)
      {
        expected.push_back(f);
      }
    }
  }
  auto by_frequency = [](auto a, auto b) { return a.real() < b.real(); };
  std::ranges::sort(expected, by_frequency);
  REQUIRE(expected.size() >= 6);
  REQUIRE(expected[5].real() < f0);  // All on the target's side of the resonance
  std::vector<std::complex<double>> found;
  std::ranges::copy_if(lorentz.f, std::back_inserter(found),
                       [](auto f) { return f.real() >= 8.0; });
  std::ranges::sort(found, by_frequency);
  REQUIRE(found.size() >= 6);
  for (std::size_t k = 0; k < 6; k++)
  {
    CAPTURE(k, found[k], expected[k]);
    CHECK(std::abs(found[k] - expected[k]) <= 1.0e-7 * std::abs(expected[k]));
  }
}

TEST_CASE_METHOD(test::SharedTempDir,
                 "Eigenmode driver: Hybrid seed of a Lorentz-filled PEC box",
                 "[eigenmodedriver][Serial][Parallel]")
{
#if !defined(PALACE_WITH_SLEPC)
  SKIP("Requires the SLEPc eigensolver");
#endif
  // The Hybrid seed freezes the permittivity at the target in its λ² block. In a PEC box
  // filled with one material, its eigenvalues are then exactly the discrete constant modes
  // f_h of ε∞ rescaled to ε(target): f = f_h √(ε∞ / ε(target)). With RefineNonlinear =
  // false these seeds are the output, before Quasi-Newton refines them. A permittivity
  // interpolated across the resonance (here at 14 GHz, between the nodes 8, 16, and 24 GHz)
  // gives other seeds, about 18 GHz for the lowest.
  const auto mesh = WriteBoxMesh(temp_dir);
  const json pec = {{"PEC", {{"Attributes", {1, 2, 3, 4, 5, 6}}}}};
  const auto constant =
      RunEigenmode(BoxConfig(mesh, eps_inf, pec, 16, 8.0, ""), temp_dir, "constant");
  auto config = BoxConfig(mesh, LorentzPermittivityConfig(), pec, 6, 8.0, "");
  config["Solver"]["Eigenmode"]["RefineNonlinear"] = false;
  config["Solver"]["Eigenmode"]["LinearTol"] = 1.0e-10;
  const auto seed = RunEigenmode(config, temp_dir, "seed");
  REQUIRE(seed.f.size() >= 6);
  const auto scale = std::sqrt(eps_inf / LorentzPermittivity(8.0));
  for (const auto f : seed.f)
  {
    CAPTURE(f);
    double rel = std::numeric_limits<double>::infinity();
    for (const auto f_h : constant.f)
    {
      rel = std::min(rel, std::abs(f - f_h * scale) / std::abs(f_h * scale));
    }
    CHECK(rel <= 1.0e-7);
  }
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
  // material. Each mode found must be an eigenpair (backward error) of a damped mode, and
  // the mode of lowest frequency above the target, 19.108678 + 2.399168i GHz (a reference
  // value with a backward error below 1e-8), must be among them. Which of the strongly
  // damped modes next to it are found depends on the seed: with one process the next
  // three were 22.126 + 26.837i, 22.288 + 6.409i, and 22.567 + 25.977i GHz, with two
  // processes 22.288 + 6.409i, 22.567 + 25.977i, and 23.035 + 2.534i GHz, and without the
  // frozen terms next to the material 22.288 + 6.409i, 23.035 + 2.534i, and
  // 23.291 + 2.215i GHz. This test therefore does not detect a missing freeze.
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
  for (std::size_t i = 0; i < hybrid.f.size(); i++)
  {
    CAPTURE(hybrid.f[i]);
    CHECK(hybrid.backward_error[i] <= 1.0e-8);
    CHECK(hybrid.f[i].imag() > 0.0);  // Radiating through the absorbing faces
  }
  const std::complex<double> lowest(19.108678, 2.399168);
  CHECK(std::ranges::any_of(hybrid.f, [&](auto f)
                            { return std::abs(f - lowest) <= 1.0e-6 * std::abs(lowest); }));
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
  // medium, and the others PEC, with the target below the resonance. The Hybrid seed
  // freezes the permittivity and the first-order absorbing boundary coefficient at the
  // target: interpolated across the resonance, they gave seeds that were boundary gradients
  // of the absorbing faces, from which Quasi-Newton found no mode. The requested modes must
  // be the lowest found independently by the SLP solver, which does not use the seed, and
  // match the reference values (backward errors below 1e-10).
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
  auto by_frequency = [](auto a, auto b) { return a.real() < b.real(); };
  auto found = hybrid.f, lowest = slp.f;
  std::ranges::sort(found, by_frequency);
  std::ranges::sort(lowest, by_frequency);
  const std::vector<std::complex<double>> reference = {{13.37300, 0.5602253},
                                                       {13.70413, 0.2946123},
                                                       {13.82612, 0.2657585},
                                                       {13.83002, 0.2586597}};
  for (std::size_t i = 0; i < n; i++)
  {
    CAPTURE(i, found[i], lowest[i], reference[i]);
    CHECK(std::abs(found[i] - lowest[i]) <= 1.0e-6 * std::abs(lowest[i]));
    CHECK(std::abs(found[i] - reference[i]) <= 1.0e-6 * std::abs(reference[i]));
  }
  for (std::size_t i = 0; i < hybrid.f.size(); i++)
  {
    CAPTURE(hybrid.f[i]);
    CHECK(hybrid.backward_error[i] <= 1.0e-8);
  }
}

}  // namespace palace
