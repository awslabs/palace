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
#include <variant>
#include <vector>
#include <Eigen/Dense>
#include <mfem.hpp>
#include <catch2/catch_test_macros.hpp>
#include <nlohmann/json.hpp>
#include "driver.hpp"
#include "fem/mesh.hpp"
#include "fixtures.hpp"
#include "linalg/vector.hpp"
#include "models/spaceoperator.hpp"
#include "test-helpers.hpp"
#include "utils/communication.hpp"
#include "utils/geodata.hpp"
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
  // Each seed is a rescaled constant mode, and the output is the eigenpairs of one seed
  // solve, in order of frequency: no constant mode appears more often than its
  // multiplicity, as it could in the pool of two seed solves used for Quasi-Newton.
  std::vector<int> count(constant.f.size(), 0);
  for (std::size_t i = 0; i < seed.f.size(); i++)
  {
    const auto f = seed.f[i];
    CAPTURE(i, f);
    if (i > 0)
    {
      CHECK(f.real() >= seed.f[i - 1].real());
    }
    double rel = std::numeric_limits<double>::infinity();
    for (std::size_t h = 0; h < constant.f.size(); h++)
    {
      const auto ref = constant.f[h] * scale;
      rel = std::min(rel, std::abs(f - ref) / std::abs(ref));
      if (std::abs(f - ref) <= 1.0e-7 * std::abs(ref))
      {
        count[h]++;
      }
    }
    CHECK(rel <= 1.0e-7);
  }
  for (std::size_t h = 0; h < constant.f.size(); h++)
  {
    // The multiplicity of f_h among the constant modes (whose list holds each copy).
    int multiplicity = 0;
    for (const auto g : constant.f)
    {
      multiplicity += (std::abs(g - constant.f[h]) <= 1.0e-7 * std::abs(constant.f[h]));
    }
    CAPTURE(constant.f[h], multiplicity);
    CHECK(count[h] <= multiplicity);
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
  // (ε∞ = 2, Δε = 0.5, τ = 0.01 ns), and the others PEC. Each reported mode must satisfy
  // the full operator's backward-error check and have damping. The lowest mode above
  // 18 GHz was derived from the complete constant-material ABC2 polynomial spectrum,
  // mapped through the Debye law and checked on the current full operator. The separate
  // [debyeoracle] case compares all four lowest modes by that independent construction.
  // The old 19.108678 + 2.399168i reference used the opposite ABC2 boundary sign.
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
  const std::complex<double> lowest(18.78468849275, 10.05621002136);
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
  // medium, and the others PEC, with the target below the resonance. The references use
  // the current ABC2 sign: the complete constant-material spectrum is mapped through all
  // Lorentz quartic branches and checked on the full operator in [lorentzoracle]. That
  // independent check retains distinct nonlinear roots which share the same field; it
  // does not rely on an SLP reference solve with physical-field deflation.
  const auto mesh = WriteBoxMesh(temp_dir);
  const json boundaries = {{"PEC", {{"Attributes", {1, 2, 4, 6}}}},
                           {"Absorbing", {{"Attributes", {3, 5}}, {"Order", 2}}}};
  constexpr int n = 4;
  const auto hybrid = RunEigenmode(
      BoxConfig(mesh, LorentzPermittivityConfig(), boundaries, n, 8.0, "Hybrid"), temp_dir,
      "hybrid");
  REQUIRE(hybrid.f.size() >= n);
  auto found = hybrid.f;
  std::ranges::sort(found, [](auto a, auto b) { return a.real() < b.real(); });
  const std::vector<std::complex<double>> reference = {{12.46381247035, 1.451488651821},
                                                       {13.61113256125, 0.3868815645070},
                                                       {13.80584778635, 0.2800203047645},
                                                       {13.82612379780, 0.2657585107049}};
  for (std::size_t i = 0; i < n; i++)
  {
    CAPTURE(i, found[i], reference[i]);
    CHECK(std::abs(found[i] - reference[i]) <= 1.0e-6 * std::abs(reference[i]));
  }
  for (std::size_t i = 0; i < hybrid.f.size(); i++)
  {
    CAPTURE(hybrid.f[i]);
    CHECK(hybrid.backward_error[i] <= 1.0e-8);
  }
}

// Frequency-dependent terms which are not rational in the frequency.

TEST_CASE_METHOD(
    test::SharedTempDir,
    "Eigenmode driver: Djordjevic-Sarkar-filled PEC box from the constant modes",
    "[eigenmodedriver][Serial][Parallel]")
{
#if !defined(PALACE_WITH_SLEPC)
  SKIP("Requires the SLEPc eigensolver");
#endif
  // As the Lorentz-filled box: in a PEC box filled with one material, the discrete modes of
  // ε(f) are those of the constant ε∞ rescaled, f² ε(f) = f_h² ε∞, with the
  // Djordjevic-Sarkar ε(f) = ε∞ + S log((i f + f_u) / (i f + f_l)) (the factors 2π cancel).
  // The oracle is one scalar Newton solve per constant mode f_h (independent of the
  // eigensolver); the requested modes are the N lowest of these above the target.
  constexpr double strength = 0.05, f_lo = 0.1, f_hi = 100.0;
  const std::complex<double> i(0.0, 1.0);
  auto Eps = [&](std::complex<double> f)
  { return eps_inf + strength * std::log((i * f + f_hi) / (i * f + f_lo)); };
  auto dEps = [&](std::complex<double> f)
  { return strength * (i / (i * f + f_hi) - i / (i * f + f_lo)); };
  const auto mesh = WriteBoxMesh(temp_dir);
  const json pec = {{"PEC", {{"Attributes", {1, 2, 3, 4, 5, 6}}}}};
  const json ds = {{"HighFrequency", eps_inf},
                   {"Terms",
                    {{{"Type", "DjordjevicSarkar"},
                      {"Strength", strength},
                      {"LowerFrequency", f_lo},
                      {"UpperFrequency", f_hi}}}}};
  constexpr int n = 6;
  const auto constant =
      RunEigenmode(BoxConfig(mesh, eps_inf, pec, 16, 8.0, ""), temp_dir, "constant");
  const auto found_all = RunEigenmode(BoxConfig(mesh, ds, pec, n, 8.0, ""), temp_dir, "ds");
  REQUIRE(constant.f.size() >= 16);
  REQUIRE(found_all.f.size() >= n);
  std::vector<std::complex<double>> expected;
  for (const auto f_h : constant.f)
  {
    const double h2 = std::norm(f_h) * eps_inf;  // f_h is real for a lossless box
    std::complex<double> f = std::abs(f_h) * std::sqrt(eps_inf / Eps(std::abs(f_h)));
    for (int it = 0; it < 100; it++)
    {
      const auto g = f * f * Eps(f) - h2, dg = 2.0 * f * Eps(f) + f * f * dEps(f);
      const auto df = g / dg;
      f -= df;
      if (std::abs(df) <= 1.0e-15 * std::abs(f))
      {
        break;
      }
    }
    REQUIRE(std::abs(f * f * Eps(f) - h2) <= 1.0e-12 * h2);
    if (f.real() >= 8.0)
    {
      expected.push_back(f);
    }
  }
  auto by_frequency = [](auto a, auto b) { return a.real() < b.real(); };
  std::ranges::sort(expected, by_frequency);
  REQUIRE(expected.size() >= n);
  std::vector<std::complex<double>> found;
  std::ranges::copy_if(found_all.f, std::back_inserter(found),
                       [](auto f) { return f.real() >= 8.0; });
  std::ranges::sort(found, by_frequency);
  REQUIRE(found.size() >= n);
  for (std::size_t k = 0; k < n; k++)
  {
    CAPTURE(k, found[k], expected[k]);
    CHECK(std::abs(found[k] - expected[k]) <= 1.0e-7 * std::abs(expected[k]));
  }
  for (std::size_t k = 0; k < found_all.f.size(); k++)
  {
    CAPTURE(found_all.f[k], found_all.backward_error[k]);
    CHECK(found_all.backward_error[k] <= 1.0e-8);
    CHECK(found_all.f[k].imag() > 0.0);  // Lossy
  }
}

TEST_CASE_METHOD(test::SharedTempDir,
                 "Eigenmode driver: surface conductivity on a face of a PEC box",
                 "[eigenmodedriver][Serial][Parallel]")
{
#if !defined(PALACE_WITH_SLEPC)
  SKIP("Requires the SLEPc eigensolver");
#endif
  // The face z = 1 has a finite surface conductivity (a principal √ in ω, not rational),
  // the others PEC, in the constant ε∞ box. For zero thickness, z = √λ gives an exact
  // quartic operator. [conductivityoracle] computes its complete spectrum and checks the
  // four references below on the full native T, independently of Hybrid and SLP.
  // Hybrid must return those lowest four with backward errors below 1e-8 and damping.
  const auto mesh = WriteBoxMesh(temp_dir);
  const json boundaries = {
      {"PEC", {{"Attributes", {1, 2, 3, 4, 5}}}},
      {"Conductivity", {{{"Attributes", {6}}, {"Conductivity", 1.0e5}}}}};
  constexpr int n = 4;
  const auto hybrid = RunEigenmode(BoxConfig(mesh, eps_inf, boundaries, n, 8.0, "Hybrid"),
                                   temp_dir, "hybrid");
  REQUIRE(hybrid.f.size() >= n);
  auto found = hybrid.f;
  std::ranges::sort(found, [](auto a, auto b) { return a.real() < b.real(); });
  const std::vector<std::complex<double>> lowest = {{17.02214338926, 0.008699223002311},
                                                    {20.66424920490, 0.01433407222388},
                                                    {22.15176051954, 0.01291535754073},
                                                    {24.56454562049, 0.02128372813413}};
  for (std::size_t k = 0; k < n; k++)
  {
    CAPTURE(k, found[k], lowest[k]);
    CHECK(std::abs(found[k] - lowest[k]) <= 1.0e-6 * std::abs(lowest[k]));
  }
  for (std::size_t k = 0; k < hybrid.f.size(); k++)
  {
    CAPTURE(hybrid.f[k], hybrid.backward_error[k]);
    CHECK(hybrid.backward_error[k] <= 1.0e-8);
    CHECK(hybrid.f[k].imag() > 0.0);
  }
}

TEST_CASE_METHOD(test::SharedTempDir,
                 "Debye ABC2: current full operator at the historical reference",
                 "[debyereference][Serial]")
{
#if !defined(PALACE_WITH_SLEPC)
  SKIP("Requires the SLEPc eigensolver");
#endif
  using namespace std::complex_literals;
  REQUIRE(Mpi::Size(Mpi::World()) == 1);
  test::IntegrationSettingsGuard settings_guard;
  const auto mesh_path = WriteBoxMesh(temp_dir);
  const json boundaries = {{"PEC", {{"Attributes", {1, 2, 4, 6}}}},
                           {"Absorbing", {{"Attributes", {3, 5}}, {"Order", 2}}}};
  const json debye = {
      {"HighFrequency", eps_inf},
      {"Terms",
       {{{"Type", "Debye"}, {"DeltaPermittivity", 0.5}, {"RelaxationTime", 0.01}}}}};
  IoData iodata(BoxConfig(mesh_path, debye, boundaries, 4, 18.0, "Hybrid"), false);

  // Use the driver's mesh preparation and nondimensionalization, not a recreated matrix.
  auto serial_mesh = mesh::Load(iodata, Mpi::World());
  if (iodata.model.Lc <= 0.0)
  {
    iodata.model.Lc = mesh::ComputeReferenceLength(serial_mesh, Mpi::World());
  }
  iodata.NondimensionalizeInputs(serial_mesh);
  std::vector<std::unique_ptr<mfem::ParMesh>> par_mesh;
  par_mesh.push_back(mesh::Partition(iodata, std::move(serial_mesh), Mpi::World()));
  mesh::RefineMesh(iodata, par_mesh);
  std::vector<std::unique_ptr<Mesh>> mesh;
  for (auto &m : par_mesh)
  {
    mesh.push_back(std::make_unique<Mesh>(std::move(m)));
  }

  SpaceOperator space_op(iodata, mesh);
  auto K = space_op.GetStiffnessMatrix<ComplexOperator>(Operator::DIAG_ONE);
  auto C = space_op.GetDampingMatrix<ComplexOperator>(Operator::DIAG_ZERO);
  auto M = space_op.GetMassMatrix<ComplexOperator>(Operator::DIAG_ZERO);
  REQUIRE(K);
  REQUIRE(M);
  const std::complex<double> frequency(19.108678, 2.399168);
  const auto omega =
      2.0 * M_PI * iodata.units.Nondimensionalize<Units::ValueType::FREQUENCY>(frequency);
  const auto lambda = 1i * omega;
  auto A2 = space_op.GetExtraSystemOperator(omega, Operator::DIAG_ZERO);
  REQUIRE(A2);
  auto T = space_op.GetSystemMatrix(1.0 + 0.0i, lambda, lambda * lambda, K.get(), C.get(),
                                    M.get(), A2.get());
  REQUIRE(T);
  REQUIRE(T->Height() == T->Width());
  const int total = T->Width();
  REQUIRE(total == 300);
  const auto &essential = space_op.GetNDDbcTDofLists().back();
  REQUIRE(essential.Size() == 144);
  std::vector<bool> constrained(total, false);
  for (int dof : essential)
  {
    REQUIRE(dof >= 0);
    REQUIRE(dof < total);
    constrained[dof] = true;
  }
  std::vector<int> free;
  for (int dof = 0; dof < total; dof++)
  {
    if (!constrained[dof])
    {
      free.push_back(dof);
    }
  }
  const int dim = static_cast<int>(free.size());
  REQUIRE(dim == 156);

  // Every free Euclidean unit vector is included; gradient DOFs are not projected out.
  // Remove both essential rows and columns so eliminated PEC rows cannot fake a nullspace.
  Eigen::MatrixXcd physical(dim, dim);
  ComplexVector x(total), y(total);
  x.UseDevice(true);
  y.UseDevice(true);
  std::vector<std::complex<double>> column(total);
  for (int j = 0; j < dim; j++)
  {
    x = 0.0;
    x.Real()(free[j]) = 1.0;
    T->Mult(x, y);
    y.Get(column.data(), total, false);
    for (int i = 0; i < dim; i++)
    {
      physical(i, j) = column[free[i]];
    }
  }
  REQUIRE(physical.allFinite());
  const Eigen::JacobiSVD<Eigen::MatrixXcd> svd(physical);
  REQUIRE(svd.info() == Eigen::Success);
  REQUIRE(svd.singularValues().allFinite());
  const double sigma_max = svd.singularValues().maxCoeff();
  const double sigma_min = svd.singularValues().minCoeff();
  REQUIRE(sigma_max > 0.0);
  const double ratio = sigma_min / sigma_max;
  REQUIRE(std::isfinite(ratio));
  // Raw conditioning measurement only: no root-acceptance threshold or disk certificate.
  Mpi::Print("Debye ABC2 full-T point: total {:d}, PEC {:d}, free {:d}, "
             "f = {:.12e}{:+.12e}i GHz, lambda = {:.12e}{:+.12e}i, "
             "sigma_min = {:.12e}, sigma_max = {:.12e}, ratio = {:.12e}\n",
             total, essential.Size(), dim, frequency.real(), frequency.imag(),
             lambda.real(), lambda.imag(), sigma_min, sigma_max, ratio);
}

namespace
{

enum class FixedBoxSpectrumFixture
{
  Debye,
  Lorentz,
  Conductivity
};

// Only these three fixed box laws share this extraction and validation path.
void CheckFixedBoxSpectrumLowestFour(const fs::path &temp_dir,
                                     FixedBoxSpectrumFixture fixture)
{
  using namespace std::complex_literals;
  const auto comm = Mpi::World();
  REQUIRE(Mpi::Size(comm) == 1);
  test::IntegrationSettingsGuard settings_guard;
  const bool is_debye = fixture == FixedBoxSpectrumFixture::Debye;
  const bool is_conductivity = fixture == FixedBoxSpectrumFixture::Conductivity;
  const char *model = is_conductivity ? "Conductivity" : (is_debye ? "Debye" : "Lorentz");
  const double target_ghz = is_debye ? 18.0 : 8.0;
  const auto mesh_path = WriteBoxMesh(temp_dir);
  const json boundaries =
      is_conductivity
          ? json{{"PEC", {{"Attributes", {1, 2, 3, 4, 5}}}},
                 {"Conductivity", {{{"Attributes", {6}}, {"Conductivity", 1.0e5}}}}}
          : json{{"PEC", {{"Attributes", {1, 2, 4, 6}}}},
                 {"Absorbing", {{"Attributes", {3, 5}}, {"Order", 2}}}};
  const json debye = {
      {"HighFrequency", eps_inf},
      {"Terms",
       {{{"Type", "Debye"}, {"DeltaPermittivity", 0.5}, {"RelaxationTime", 0.01}}}}};
  const json permittivity =
      is_conductivity ? json(eps_inf) : (is_debye ? debye : LorentzPermittivityConfig());
  constexpr int n = 4;
  IoData iodata(BoxConfig(mesh_path, permittivity, boundaries, n, target_ghz, "Hybrid"),
                false);
  auto serial_mesh = mesh::Load(iodata, comm);
  if (iodata.model.Lc <= 0.0)
  {
    iodata.model.Lc = mesh::ComputeReferenceLength(serial_mesh, comm);
  }
  iodata.NondimensionalizeInputs(serial_mesh);
  std::vector<std::unique_ptr<mfem::ParMesh>> par_mesh;
  par_mesh.push_back(mesh::Partition(iodata, std::move(serial_mesh), comm));
  mesh::RefineMesh(iodata, par_mesh);
  std::vector<std::unique_ptr<Mesh>> mesh;
  for (auto &m : par_mesh)
  {
    mesh.push_back(std::make_unique<Mesh>(std::move(m)));
  }

  SpaceOperator space_op(iodata, mesh);
  auto K = space_op.GetStiffnessMatrix<ComplexOperator>(Operator::DIAG_ONE);
  auto C = space_op.GetDampingMatrix<ComplexOperator>(Operator::DIAG_ZERO);
  auto M = space_op.GetMassMatrix<ComplexOperator>(Operator::DIAG_ZERO);
  REQUIRE(K);
  REQUIRE_FALSE(C);  // All frequency-dependent boundary damping is in A2.
  REQUIRE(M);
  std::unique_ptr<ComplexOperator> damping, curlcurl, surface_mass;
  double kappa = 0.0;
  if (is_conductivity)
  {
    REQUIRE(space_op.GetSurfaceConductivityOp().Size() == 1);
    const auto &law = iodata.boundaries.conductivity.front();
    REQUIRE(law.sigma > 0.0);
    REQUIRE(law.mu_r == 1.0);
    REQUIRE(law.h == 0.0);
    REQUIRE_FALSE(law.external);
    // sigma is already sigma_hat = sigma_SI Z0 Lc. Face 6 is exterior (weight 1).
    kappa = std::sqrt(law.sigma / law.mu_r);
    REQUIRE(std::isfinite(kappa));
    surface_mass = space_op.GetSurfaceConductivityBoundaryMatrix<ComplexOperator>(
        0, Operator::DIAG_ZERO);
    REQUIRE(surface_mass);
  }
  else
  {
    REQUIRE(space_op.GetFarfieldOp().NumDispersiveMaterials() == 1);
    damping = space_op.GetFarfieldDispersiveDampingMatrix(0, Operator::DIAG_ZERO);
    curlcurl = space_op.GetFarfieldDispersiveCurlCurlMatrix(0, Operator::DIAG_ZERO);
    REQUIRE(damping);
    REQUIRE(curlcurl);
  }
  const int total = space_op.GetNDSpace().GetTrueVSize();
  REQUIRE(total == 300);
  const auto &essential = space_op.GetNDDbcTDofLists().back();
  REQUIRE(essential.Size() == (is_conductivity ? 168 : 144));
  std::vector<bool> constrained(total, false);
  for (int dof : essential)
  {
    REQUIRE(dof >= 0);
    REQUIRE(dof < total);
    constrained[dof] = true;
  }
  std::vector<int> free;
  for (int dof = 0; dof < total; dof++)
  {
    if (!constrained[dof])
    {
      free.push_back(dof);
    }
  }
  const int dim = static_cast<int>(free.size());
  REQUIRE(dim == (is_conductivity ? 132 : 156));
  ComplexVector x(total), y(total);
  x.UseDevice(true);
  y.UseDevice(true);
  std::vector<std::complex<double>> column(total);
  // Full free Euclidean basis, including gradients; remove both PEC rows and columns.
  auto Extract = [&](const ComplexOperator &op) -> Eigen::MatrixXcd
  {
    REQUIRE(op.Height() == total);
    REQUIRE(op.Width() == total);
    Eigen::MatrixXcd dense(dim, dim);
    for (int j = 0; j < dim; j++)
    {
      x = 0.0;
      x.Real().HostReadWrite()[free[j]] = 1.0;
      op.Mult(x, y);
      y.Get(column.data(), total, false);
      for (int i = 0; i < dim; i++)
      {
        dense(i, j) = column[free[i]];
      }
    }
    REQUIRE(dense.allFinite());
    return dense;
  };
  const Eigen::MatrixXcd Kinf = Extract(*K), Minf = Extract(*M);
  const Eigen::LLT<Eigen::MatrixXcd> mass(Minf);
  REQUIRE(mass.info() == Eigen::Success);
  Eigen::MatrixXcd companion;
  if (is_conductivity)
  {
    // The owner returns iB. For h=0, delta=sqrt(2/(mu sigma omega)),
    // Z=(1+i)/(sigma delta), so s(lambda)=kappa sqrt(lambda) in Im{lambda}>0.
    // P(z)=Minf z^4+kappa B z+Kinf, lambda=z^2, state (u,z u,z^2 u,z^3 u).
    // Minf is SPD: no clearing, poles or infinity. Zero/gradient copies are retained.
    const Eigen::MatrixXcd B = Extract(*surface_mass) / 1i;
    companion = Eigen::MatrixXcd::Zero(4 * dim, 4 * dim);
    companion.block(0, dim, dim, dim).setIdentity();
    companion.block(dim, 2 * dim, dim, dim).setIdentity();
    companion.block(2 * dim, 3 * dim, dim, dim).setIdentity();
    companion.block(3 * dim, 0, dim, dim) = -mass.solve(Kinf);
    companion.block(3 * dim, dim, dim, dim) = -kappa * mass.solve(B);
  }
  else
  {
    // These owner matrices carry i, including the current negative sign in Mff.
    const Eigen::MatrixXcd Cinf = Extract(*damping) / 1i, Mff = Extract(*curlcurl) / 1i;
    // Tinf(mu)=K+mu Cinf+mu^2 Minf-0.5/mu Mff. Clear its pole, then use
    // the strong cubic companion with state (u,mu u,mu^2 u). No inverse or subspace.
    companion = Eigen::MatrixXcd::Zero(3 * dim, 3 * dim);
    companion.block(0, dim, dim, dim).setIdentity();
    companion.block(dim, 2 * dim, dim, dim).setIdentity();
    companion.block(2 * dim, 0, dim, dim) = 0.5 * mass.solve(Mff);
    companion.block(2 * dim, dim, dim, dim) = -mass.solve(Kinf);
    companion.block(2 * dim, 2 * dim, dim, dim) = -mass.solve(Cinf);
  }
  REQUIRE(companion.allFinite());
  const Eigen::ComplexEigenSolver<Eigen::MatrixXcd> spectrum(companion);
  REQUIRE(spectrum.info() == Eigen::Success);
  REQUIRE(spectrum.eigenvalues().allFinite());
  REQUIRE(spectrum.eigenvectors().allFinite());

  const double a = is_debye ? 0.25 : 0.05;
  double tau = 0.0, omega0_squared = 0.0, gamma = 0.0;
  if (!is_conductivity)
  {
    const auto &term = iodata.domains.materials.front().permittivity_terms.front();
    if (is_debye)
    {
      const auto &law = std::get<config::DebyeData>(term);
      tau = law.tau;  // Already nondimensionalized, not the input ns value.
      REQUIRE(law.delta == 0.5);
      REQUIRE(tau > 0.0);
    }
    else
    {
      const auto &law = std::get<config::LorentzData>(term);
      omega0_squared = law.omega_0 * law.omega_0;
      gamma = law.gamma;  // Both rates already use the owner's nondimensional units.
      REQUIRE(law.delta == delta_eps);
      REQUIRE(omega0_squared > 0.0);
      REQUIRE(gamma > 0.0);
    }
  }
  struct Candidate
  {
    std::complex<double> lambda;
    int mode;
  };
  std::vector<Candidate> candidates;
  int mapped_values = 0;
  for (int mode = 0; mode < spectrum.eigenvalues().size(); mode++)
  {
    const auto mu = spectrum.eigenvalues()(mode);  // z for the conductivity fixture.
    Eigen::VectorXcd values;
    if (is_conductivity)
    {
      values.resize(1);
      values(0) = mu * mu;  // Every z maps to lambda=z^2; no extra root search.
    }
    else if (is_debye)
    {
      // mu = lambda r, r^2 = 1 + a/(1 + tau lambda): all three scalar roots.
      Eigen::Matrix3cd scalar = Eigen::Matrix3cd::Zero();
      scalar(0, 1) = scalar(1, 2) = 1.0;
      scalar(2, 0) = mu * mu / tau;
      scalar(2, 1) = mu * mu;
      scalar(2, 2) = -(1.0 + a) / tau;
      REQUIRE(scalar.allFinite());
      const Eigen::ComplexEigenSolver<Eigen::Matrix3cd> roots(scalar, false);
      REQUIRE(roots.info() == Eigen::Success);
      values = roots.eigenvalues();
      // Cleared mu=0 roots map to lambda=0 or -(1+a)/tau, both below this target.
      // The Debye pole is not a scalar root. Stationary points of lambda r have
      // Re{f} = sqrt(8a-a^2)/(8 pi tau_ns), about 5.54 GHz, also below 18 GHz.
    }
    else
    {
      // D(lambda)(lambda^2-mu^2) + a omega0^2 lambda^2 = 0,
      // D = lambda^2 + gamma lambda + omega0^2. The quartic is monic.
      Eigen::Matrix4cd scalar = Eigen::Matrix4cd::Zero();
      scalar(0, 1) = scalar(1, 2) = scalar(2, 3) = 1.0;
      scalar(3, 0) = omega0_squared * mu * mu;
      scalar(3, 1) = gamma * mu * mu;
      scalar(3, 2) = mu * mu - (1.0 + a) * omega0_squared;
      scalar(3, 3) = -gamma;
      REQUIRE(scalar.allFinite());
      const Eigen::ComplexEigenSolver<Eigen::Matrix4cd> roots(scalar, false);
      REQUIRE(roots.info() == Eigen::Success);
      values = roots.eigenvalues();
      // D=0 is not a root: the remaining a omega0^2 lambda^2 is nonzero there.
      // Cleared mu=0 gives lambda=0 twice and the epsilon zeros near
      // +/-14.34355+0.25i GHz. The positive zero is in band: do not delete it.
      // The owner sqrt cut is at Im{f}=0.25 between the positive-frequency pole
      // near 13.99777 and zero near 14.34355 GHz. GetDispersiveScale selects its
      // passive lossy limit; no separate analytic sqrt or stationary-point filter.
    }
    REQUIRE(values.allFinite());
    mapped_values += static_cast<int>(values.size());
    for (const auto lambda : values)
    {
      if (lambda.imag() <= iodata.solver.eigenmode.target)
      {
        continue;
      }
      // The target is tested first: Im{lambda}>0 gives Re{omega}>0 and fixes the
      // owner's principal skin-depth branch. Zero/defective gradient copies get no filter.
      const auto owner_value =
          is_conductivity
              ? space_op.GetSurfaceConductivityOp().EvaluateScalar(0, lambda / 1i)
              : lambda * space_op.GetFarfieldOp().GetDispersiveScale(0, lambda / 1i);
      const auto model_value = is_conductivity ? kappa * mu : mu;
      // Choose the analytic sheet's sign, not a distance cutoff.
      const double sheet = std::real(std::conj(model_value) * owner_value);
      REQUIRE(std::isfinite(sheet));
      if (sheet < 0.0)
      {
        continue;
      }
      if (is_debye)
      {
        REQUIRE(sheet > 0.0);
      }
      // Lorentz and conductivity sheet=0 remains a structural candidate. If selected,
      // it must pass the unsquared owner relation and original full T below.
      candidates.push_back({lambda, mode});
    }
  }
  std::stable_sort(candidates.begin(), candidates.end(),
                   [](const auto &lhs, const auto &rhs)
                   { return lhs.lambda.imag() < rhs.lambda.imag(); });
  Mpi::Print("{} fixed-box spectrum: free {:d}, companion {:d}, mapped values {:d}, "
             "outgoing or structural candidates above target {:d}\n",
             model, dim, static_cast<int>(spectrum.eigenvalues().size()), mapped_values,
             candidates.size());
  REQUIRE(candidates.size() >= n);

  // Native backward normalization: ||T(lambda)x|| / (||x|| (||K||+|lambda|^2||M||)).
  // C is absent; no A2 norm. The eliminated K block is identity, hence max(1, ||Kfree||).
  const Eigen::JacobiSVD<Eigen::MatrixXcd> k_norm(Kinf), m_norm(Minf);
  REQUIRE(k_norm.info() == Eigen::Success);
  REQUIRE(m_norm.info() == Eigen::Success);
  REQUIRE(k_norm.singularValues().allFinite());
  REQUIRE(m_norm.singularValues().allFinite());
  const double normK = std::max(1.0, k_norm.singularValues().maxCoeff());
  const double normM = m_norm.singularValues().maxCoeff();
  REQUIRE(normM > 0.0);
  std::vector<std::complex<double>> expected;
  for (int k = 0; k < n; k++)
  {
    // Select the lowest four BEFORE validation; a failed lower candidate stops the case.
    // Every companion copy and every mapped branch is retained: no field/frequency dedup.
    const auto &candidate = candidates[k];
    const auto lambda = candidate.lambda;
    const auto mu = spectrum.eigenvalues()(candidate.mode);  // z for conductivity.
    CAPTURE(k, candidate.mode, lambda, mu);
    INFO("The lowest-four prefix is fixed, including any structural zero/sheet-zero copy. "
         "A failed owner relation or full-T validation leaves the oracle inconclusive; "
         "no higher candidate substitutes.");
    if (is_conductivity)
    {
      Mpi::Print("Conductivity selected {:d}: source {:d}, z = {:.12e}{:+.12e}i, "
                 "lambda = {:.12e}{:+.12e}i\n",
                 k + 1, candidate.mode, mu.real(), mu.imag(), lambda.real(), lambda.imag());
      const auto owner_s =
          space_op.GetSurfaceConductivityOp().EvaluateScalar(0, lambda / 1i);
      REQUIRE(std::abs(owner_s - kappa * mu) <=
              iodata.solver.eigenmode.tol * kappa * std::abs(mu));
    }
    else
    {
      const auto scale = space_op.GetFarfieldOp().GetDispersiveScale(0, lambda / 1i);
      const auto mapped = lambda * scale;
      if (!is_debye)
      {
        Mpi::Print("Lorentz ABC2 selected {:d}: source {:d}, mu = {:.12e}{:+.12e}i, "
                   "lambda = {:.12e}{:+.12e}i\n",
                   k + 1, candidate.mode, mu.real(), mu.imag(), lambda.real(),
                   lambda.imag());
      }
      // No division by |mu|, so an exact clearing root is not discarded or rescaled.
      REQUIRE(std::abs(mapped - mu) <= iodata.solver.eigenmode.tol * std::abs(mu));
      if (!is_debye)
      {
        INFO("An exact epsilon-zero scale makes the original ABC2 operator undefined; "
             "a selected such candidate leaves the oracle inconclusive.");
        REQUIRE(scale != std::complex<double>{});
      }
    }
    Eigen::VectorXcd field = spectrum.eigenvectors().col(candidate.mode).head(dim);
    const double field_norm = field.norm();
    REQUIRE(std::isfinite(field_norm));
    REQUIRE(field_norm > 0.0);  // Structural, no numerical rank cutoff.
    field /= field_norm;
    x = 0.0;
    auto *xr = x.Real().HostReadWrite();
    auto *xi = x.Imag().HostReadWrite();
    for (int i = 0; i < dim; i++)
    {
      xr[free[i]] = field(i).real();
      xi[free[i]] = field(i).imag();
    }
    auto A2 = space_op.GetExtraSystemOperator(lambda / 1i, Operator::DIAG_ZERO);
    REQUIRE(A2);
    auto T = space_op.GetSystemMatrix(1.0 + 0.0i, lambda, lambda * lambda, K.get(), C.get(),
                                      M.get(), A2.get());
    REQUIRE(T);
    T->Mult(x, y);
    const double action = linalg::Norml2(comm, y) / linalg::Norml2(comm, x);
    const double backward = action / (normK + std::norm(lambda) * normM);
    const auto frequency =
        iodata.units.Dimensionalize<Units::ValueType::FREQUENCY>(lambda / 1i) /
        (2.0 * M_PI);
    Mpi::Print("{} fixed-box oracle mode {:d}: f = {:.12e}{:+.12e}i GHz, "
               "action = {:.12e}, backward = {:.12e}\n",
               model, k + 1, frequency.real(), frequency.imag(), action, backward);
    REQUIRE(std::isfinite(action));
    REQUIRE(std::isfinite(backward));
    // Use the fixture's backward-error check, not a new absolute-action gate.
    REQUIRE(backward <= 1.0e-8);
    CHECK(frequency.imag() > 0.0);  // Damping is a check, never a selection filter.
    expected.push_back(frequency);
  }

  const auto hybrid = RunEigenmode(
      BoxConfig(mesh_path, permittivity, boundaries, n, target_ghz, "Hybrid"), temp_dir,
      is_conductivity ? "hybrid-conductivity-oracle"
                      : (is_debye ? "hybrid-oracle" : "hybrid-lorentz-oracle"));
  REQUIRE(hybrid.f.size() >= n);
  auto found = hybrid.f;
  std::ranges::sort(found, [](auto lhs, auto rhs) { return lhs.real() < rhs.real(); });
  for (std::size_t k = 0; k < hybrid.f.size(); k++)
  {
    CAPTURE(k, hybrid.f[k], hybrid.backward_error[k]);
    CHECK(hybrid.backward_error[k] <= 1.0e-8);
    CHECK(hybrid.f[k].imag() > 0.0);
  }
  for (int k = 0; k < n; k++)
  {
    CAPTURE(k, found[k], expected[k]);
    Mpi::Print("{} fixed-box Hybrid comparison {:d}: found {:.12e}{:+.12e}i GHz, "
               "oracle {:.12e}{:+.12e}i GHz\n",
               model, k + 1, found[k].real(), found[k].imag(), expected[k].real(),
               expected[k].imag());
    CHECK(std::abs(found[k] - expected[k]) <= 1.0e-6 * std::abs(expected[k]));
  }
}

}  // namespace

TEST_CASE_METHOD(test::SharedTempDir,
                 "Debye ABC2: lowest four from the complete constant-material spectrum",
                 "[debyeoracle][Serial]")
{
#if !defined(PALACE_WITH_SLEPC)
  SKIP("Requires the SLEPc eigensolver");
#endif
  CheckFixedBoxSpectrumLowestFour(temp_dir, FixedBoxSpectrumFixture::Debye);
}

TEST_CASE_METHOD(test::SharedTempDir,
                 "Lorentz ABC2: lowest four from the complete constant-material spectrum",
                 "[lorentzoracle][Serial]")
{
#if !defined(PALACE_WITH_SLEPC)
  SKIP("Requires the SLEPc eigensolver");
#endif
  CheckFixedBoxSpectrumLowestFour(temp_dir, FixedBoxSpectrumFixture::Lorentz);
}

TEST_CASE_METHOD(test::SharedTempDir,
                 "Surface conductivity: lowest four from the complete skin-effect spectrum",
                 "[conductivityoracle][Serial]")
{
#if !defined(PALACE_WITH_SLEPC)
  SKIP("Requires the SLEPc eigensolver");
#endif
  CheckFixedBoxSpectrumLowestFour(temp_dir, FixedBoxSpectrumFixture::Conductivity);
}

}  // namespace palace
