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
#include <numbers>
#include <numeric>
#include <string>
#include <tuple>
#include <vector>
#include <Eigen/Dense>
#include <mfem.hpp>
#include <catch2/catch_test_macros.hpp>
#include <nlohmann/json.hpp>
#include <catch2/generators/catch_generators.hpp>
#include <catch2/generators/catch_generators_range.hpp>
#include "driver.hpp"
#include "fem/mesh.hpp"
#include "fixtures.hpp"
#include "linalg/operator.hpp"
#include "linalg/vector.hpp"
#include "models/spaceoperator.hpp"
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
               int n, double target, const std::string &nonlinear_type,
               const std::string &backend = "SLEPc")
{
  json eigenmode = {{"N", n}, {"Tol", 1.0e-9}, {"Target", target}, {"Type", backend}};
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

// The partially filled cavity of the review of the exact pole terms: a 3 x 2 x 1 cm PEC box
// with a Lorentz slab below z = 0.3 cm and vacuum above, here on 6 x 4 x 4 hexahedra of
// order 1 (the review used 12 x 8 x 10 of order 2), the slab one layer of elements. The
// vertices are moved in x and y, keeping the faces of the box and the plane of the slab, so
// that no mode is degenerate by a symmetry of the mesh (an exactly degenerate pair at the
// end of the requested modes can lose one member in the Krylov solve).
constexpr double slab_eps_inf = 3.8, slab_delta_eps = 0.2, slab_f0 = 10.0,
                 slab_damping = 0.1;

void MoveSlabVertex(const mfem::Vector &x, mfem::Vector &y)
{
  y = x;
  y(0) = x(0) + 0.15 * std::sin(std::numbers::pi * x(0) / 3.0) * x(1) / 2.0;
  y(1) = x(1) + 0.1 * std::sin(std::numbers::pi * x(1) / 2.0) * x(0) / 3.0;
  y(2) = (x(2) <= 0.25) ? 1.2 * x(2) : 0.3 + (x(2) - 0.25) * 0.7 / 0.75;
}

fs::path WriteSlabMesh(const fs::path &dir)
{
  const auto path = dir / "slab.mesh";
  if (Mpi::Root(Mpi::World()))
  {
    auto mesh =
        mfem::Mesh::MakeCartesian3D(6, 4, 4, mfem::Element::HEXAHEDRON, 3.0, 2.0, 1.0);
    mesh.Transform(MoveSlabVertex);
    for (int e = 0; e < mesh.GetNE(); e++)
    {
      mfem::Vector center;
      mesh.GetElementCenter(e, center);
      mesh.SetAttribute(e, (center(2) < 0.3) ? 2 : 1);
    }
    mesh.SetAttributes();
    std::ofstream out(path);
    mesh.Print(out);
  }
  Mpi::Barrier(Mpi::World());
  return path;
}

json SlabConfig(const fs::path &mesh, int n, double target, const std::string &backend)
{
  const json pec = {{"PEC", {{"Attributes", {1, 2, 3, 4, 5, 6}}}}};
  auto config = BoxConfig(mesh, 1.0, pec, n, target, "", backend);
  config["Domains"]["Materials"].push_back({{"Attributes", {2}},
                                            {"Permittivity",
                                             {{"HighFrequency", slab_eps_inf},
                                              {"Terms",
                                               {{{"Type", "Lorentz"},
                                                 {"DeltaPermittivity", slab_delta_eps},
                                                 {"ResonanceFrequency", slab_f0},
                                                 {"DampingFrequency", slab_damping}}}}}}});
  config["Solver"]["Order"] = 1;
  return config;
}

// The matrix of an operator on the true dofs (numbered by process), on every process.
Eigen::MatrixXcd DenseMatrix(const ComplexOperator &A, MPI_Comm comm)
{
  const int n_local = A.Height();
  std::vector<int> counts(Mpi::Size(comm)), offsets(Mpi::Size(comm), 0);
  MPI_Allgather(&n_local, 1, MPI_INT, counts.data(), 1, MPI_INT, comm);
  std::partial_sum(counts.begin(), counts.end() - 1, offsets.begin() + 1);
  const int n = offsets.back() + counts.back(), offset = offsets[Mpi::Rank(comm)];
  Eigen::MatrixXcd D(n, n);
  ComplexVector x(n_local), y(n_local);
  std::vector<double> re(n), im(n);
  for (int j = 0; j < n; j++)
  {
    x = 0.0;
    if (j >= offset && j < offset + n_local)
    {
      x.Real().HostReadWrite()[j - offset] = 1.0;
    }
    A.Mult(x, y);
    MPI_Allgatherv(y.Real().HostRead(), n_local, MPI_DOUBLE, re.data(), counts.data(),
                   offsets.data(), MPI_DOUBLE, comm);
    MPI_Allgatherv(y.Imag().HostRead(), n_local, MPI_DOUBLE, im.data(), counts.data(),
                   offsets.data(), MPI_DOUBLE, comm);
    for (int i = 0; i < n; i++)
    {
      D(i, j) = {re[i], im[i]};
    }
  }
  return D;
}

// The discrete eigenvalue problem T(λ) x = (K + λ² M + λ² χ(f) B) x = 0 of a PEC cavity
// with one Lorentz material, λ = i κ f (f in GHz), with K, M, and the unit mass B of the
// material restricted to the dofs which are not PEC, and χ(f) its closed form from the
// configuration. With the eigenpairs K V = M V Λ, Vᵀ M V = I, and B = Pᵀ L Lᵀ P on the
// support of B, det T(λ) = det M det D det(I + g(λ) Uᵀ D⁻¹ U), D = Λ + λ², g = λ² χ, U = (P
// V)ᵀ L, and its logarithmic derivative needs only the small matrices on the material.
struct DenseLorentzCavity
{
  double kappa, eps_delta, f_res, gamma;
  Eigen::VectorXd Lambda;
  Eigen::MatrixXd U;

  // χ(f) and dχ/df, the frequencies in GHz.
  std::pair<std::complex<double>, std::complex<double>> Chi(std::complex<double> f) const
  {
    const std::complex<double> i(0.0, 1.0), d = f_res * f_res - f * f + i * gamma * f;
    return {eps_delta * f_res * f_res / d,
            eps_delta * f_res * f_res * (2.0 * f - i * gamma) / (d * d)};
  }

  // (d/df) log det T(λ(f)).
  std::complex<double> LogDerivative(std::complex<double> f) const
  {
    const std::complex<double> i(0.0, 1.0), l = i * kappa * f;
    const auto [chi, dchi] = Chi(f);
    const std::complex<double> g = l * l * chi,
                               dg = 2.0 * l * chi + l * l * dchi / (i * kappa);
    const Eigen::VectorXcd D = (Lambda.array() + l * l).matrix();
    const Eigen::MatrixXcd S =
        D.cwiseInverse().asDiagonal() * U.cast<std::complex<double>>();
    const Eigen::MatrixXcd W = U.transpose() * S, dW = -2.0 * l * S.transpose() * S;
    const Eigen::MatrixXcd A = Eigen::MatrixXcd::Identity(W.rows(), W.cols()) + g * W;
    return i * kappa *
           ((2.0 * l * D.cwiseInverse()).sum() +
            A.partialPivLu().solve(dg * W + g * dW).trace());
  }

  // The number of eigenvalues inside the ellipse with the axis [a, b] of the real frequency
  // and the half height h about Im{f} = c: (1/2πi) ∮ (d/df) log det T df, by the
  // trapezoidal rule (it is not an integer if the nodes do not resolve the integrand).
  std::complex<double> Count(double a, double b, double h, double c, int nodes) const
  {
    const std::complex<double> i(0.0, 1.0);
    std::complex<double> sum = 0.0;
    for (int k = 0; k < nodes; k++)
    {
      const double t = 2.0 * std::numbers::pi * k / nodes;
      const std::complex<double> f = 0.5 * (a + b) + 0.5 * (b - a) * std::cos(t) +
                                     i * (c + h * std::sin(t)),
                                 df = -0.5 * (b - a) * std::sin(t) + i * h * std::cos(t);
      sum += LogDerivative(f) * df;
    }
    return sum / (i * static_cast<double>(nodes));
  }
};

}  // namespace

TEST_CASE_METHOD(test::SharedTempDir,
                 "Eigenmode driver: Lorentz-filled PEC box from the constant modes",
                 "[eigenmodedriver][Serial][Parallel]")
{
  // The eigenvalue problem with the pole terms of the Lorentz medium is linearized exactly
  // by each eigensolver backend of the build.
  std::vector<std::string> backends;
#if defined(PALACE_WITH_SLEPC)
  backends.push_back("SLEPc");
#endif
#if defined(PALACE_WITH_ARPACK)
  backends.push_back("ARPACK");
#endif
  const std::string backend = GENERATE_REF(from_range(backends));
  CAPTURE(backend);
  // In a PEC box filled with one material, the discrete modes of ε(ω) on the same mesh are
  // exactly those of a constant ε∞ rescaled by the material: f² ε(f) = f_h² ε∞, with f_h a
  // discrete eigenfrequency of the constant ε∞. Each Lorentz mode found must satisfy this
  // for some f_h.
  const auto mesh = WriteBoxMesh(temp_dir);
  const json pec = {{"PEC", {{"Attributes", {1, 2, 3, 4, 5, 6}}}}};
  const auto constant = RunEigenmode(BoxConfig(mesh, eps_inf, pec, 16, 8.0, "", backend),
                                     temp_dir, "constant");
  // The linearization is exact with and without the nonlinear refinement, which has
  // nothing to refine.
  const bool refine = GENERATE(true, false);
  CAPTURE(refine);
  auto lorentz_config =
      BoxConfig(mesh, LorentzPermittivityConfig(), pec, 6, 8.0, "", backend);
  lorentz_config["Solver"]["Eigenmode"]["RefineNonlinear"] = refine;
  const auto lorentz = RunEigenmode(lorentz_config, temp_dir, "lorentz");
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

  // The requested modes are the lowest above the target: the roots of the quartic
  // f² (ε∞ D(f) + Δε f₀²) = f_h² ε∞ D(f), D(f) = f₀² − f² + i γ f, for every constant mode
  // f_h. With the target below the resonance and the constant modes above it, the lowest
  // roots lie between the target and the pole, on the other side of the pole from the
  // constant modes: a root skipped, or replaced by a root beyond the pole, fails. The
  // multiplicity of a root is not checked. This box has degenerate constant modes (at
  // 24.586, 30.767, and 33.750 GHz), and the eigenvalue solver can return one mode of an
  // exactly degenerate pair where the pair ends the requested modes.
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
  // Every root up to the sixth mode found is among the modes found.
  for (const auto e : expected)
  {
    if (e.real() > found[5].real())
    {
      break;
    }
    CAPTURE(e, found[5]);
    CHECK(std::ranges::any_of(found, [e](auto f)
                              { return std::abs(f - e) <= 1.0e-7 * std::abs(e); }));
  }
}

TEST_CASE_METHOD(test::SharedTempDir,
                 "Eigenmode driver: partially filled Lorentz cavity, the modes nearest "
                 "above the target",
                 "[eigenmodedriver][Serial][Parallel]")
{
  // Below the resonance of the slab (the pole at 9.99987 + 0.05i GHz) its modes accumulate
  // at it. Above it there are modes on the surface of the slab, where Re{ε} < 0, which
  // accumulate where ε = -1, and modes above the zero of ε at 10.26 GHz.
  //
  // Independent of the solver: the number of eigenvalues of the discrete problem in each of
  // the windows of frequency [5, 9.7], [9.7, 9.988] (up to the pole), [10.1, 10.25], and
  // [10.5, 13.62] GHz, by the argument principle on det T(λ) along an ellipse about the
  // window, with dense matrices. This count depends on the discretization (the matrices K,
  // M, and B of the slab assembled by Palace, checked against A2(λ)) and on the closed form
  // of the Lorentz ε(f), not on the pole terms of the material model, the linearization,
  // the eigenvalue solver, its shift, or its linear solves.
  //
  // From the solver: one long list of the 60 modes nearest above 5 GHz, which must hold as
  // many modes in each window as are counted, so that it is complete in the windows. The
  // solves for other targets and mode counts, below the pole, close below it, just above
  // it, and above the zero of ε, each from the start of a window, must then return the
  // modes nearest above the target in that list.
  const auto mesh = WriteSlabMesh(temp_dir);
  MPI_Comm comm = Mpi::World();
  const std::complex<double> i(0.0, 1.0);
  // Each window [a, b] with the half height h of its ellipse about Im{f} = c (the modes
  // have 0 < Im{f} < 0.05 GHz).
  constexpr double c = 0.025;
  const std::vector<std::tuple<double, double, double>> windows = {
      {5.0, 9.7, 0.5}, {9.7, 9.988, 0.1}, {10.1, 10.25, 0.1}, {10.5, 13.62, 0.5}};
  std::vector<long> counts;
  {
    IoData iodata(SlabConfig(mesh, 1, 5.0, "Default"), false);
    auto smesh = mesh::Load(iodata, comm);
    iodata.model.Lc = mesh::ComputeReferenceLength(smesh, comm);  // As the driver
    iodata.NondimensionalizeInputs(smesh);
    std::vector<std::unique_ptr<Mesh>> par_mesh;
    par_mesh.push_back(
        std::make_unique<Mesh>(mesh::Partition(iodata, std::move(smesh), comm)));
    SpaceOperator space_op(iodata, par_mesh);
    const auto K = DenseMatrix(
        *space_op.GetStiffnessMatrix<ComplexOperator>(Operator::DIAG_ONE), comm);
    const auto M =
        DenseMatrix(*space_op.GetMassMatrix<ComplexOperator>(Operator::DIAG_ZERO), comm);
    const auto B_op =
        space_op.GetFrequencyDependentPermittivityMassMatrix(1, Operator::DIAG_ZERO);
    REQUIRE(B_op);
    const Eigen::MatrixXcd B = DenseMatrix(*B_op, comm) / i;  // On the imaginary slot
    DenseLorentzCavity cavity{
        2.0 * std::numbers::pi *
            iodata.units.Nondimensionalize<Units::ValueType::FREQUENCY>(1.0),
        slab_delta_eps, slab_f0, slab_damping};
    const std::complex<double> f_check(7.3, 0.4), l_check = i * cavity.kappa * f_check;
    const auto A2 =
        DenseMatrix(*space_op.GetExtraSystemMatrix(l_check / i, Operator::DIAG_ZERO), comm);
    CHECK((A2 - l_check * l_check * cavity.Chi(f_check).first * B).norm() <=
          1.0e-12 * A2.norm());
    CHECK(K.imag().norm() == 0.0);
    CHECK(M.imag().norm() == 0.0);

    // The dofs which are not PEC, where M has a nonzero diagonal, and among them the
    // support of B.
    std::vector<int> dofs, support;
    for (int j = 0; j < M.rows(); j++)
    {
      if (M(j, j) != 0.0)
      {
        if (B(j, j) != 0.0)
        {
          support.push_back(static_cast<int>(dofs.size()));
        }
        dofs.push_back(j);
      }
    }
    const Eigen::MatrixXd K_r = K(dofs, dofs).real(), M_r = M(dofs, dofs).real(),
                          B_s = B(dofs, dofs).real()(support, support);
    CHECK(std::abs(B_s.norm() - B.norm()) <= 1.0e-14 * B.norm());  // Zero outside
    const Eigen::GeneralizedSelfAdjointEigenSolver<Eigen::MatrixXd> modes(K_r, M_r);
    const Eigen::LLT<Eigen::MatrixXd> llt(B_s);
    REQUIRE(modes.info() == Eigen::Success);
    REQUIRE(llt.info() == Eigen::Success);
    const Eigen::MatrixXd PV = modes.eigenvectors()(support, Eigen::seqN(0, K_r.cols()));
    cavity.Lambda = modes.eigenvalues();
    cavity.U = PV.transpose() * Eigen::MatrixXd(llt.matrixL());
    for (const auto &[a, b, h] : windows)
    {
      const auto count = cavity.Count(a, b, h, c, 256);
      CAPTURE(a, b, count);
      CHECK(std::abs(count - std::round(count.real())) <= 1.0e-5);
      counts.push_back(std::lround(count.real()));
    }
  }

  std::vector<std::string> backends;
#if defined(PALACE_WITH_SLEPC)
  backends.push_back("SLEPc");
#endif
#if defined(PALACE_WITH_ARPACK)
  backends.push_back("ARPACK");
#endif
  // The first n modes of a solve (those nearest above its target), by frequency.
  auto Modes = [](const Eigenvalues &e, std::size_t n)
  {
    REQUIRE(e.f.size() >= n);
    for (std::size_t k = 0; k < n; k++)
    {
      CAPTURE(e.f[k]);
      CHECK(e.backward_error[k] <= 1.0e-8);
    }
    std::vector<std::complex<double>> f(e.f.begin(), e.f.begin() + n);
    std::ranges::sort(f, [](auto a, auto b) { return a.real() < b.real(); });
    return f;
  };
  for (const auto &backend : backends)
  {
    CAPTURE(backend);
    const auto list =
        Modes(RunEigenmode(SlabConfig(mesh, 60, 5.0, backend), temp_dir, backend), 60);
    for (std::size_t w = 0; w < windows.size(); w++)
    {
      const auto [a, b, h] = windows[w];
      CAPTURE(a, b);
      REQUIRE(list.back().real() > b);
      auto inside = [a, b, h](auto f)
      {
        return std::norm((f.real() - 0.5 * (a + b)) / (0.5 * (b - a))) +
                   std::norm((f.imag() - c) / h) <
               1.0;
      };
      CHECK(std::ranges::count_if(list, inside) == counts[w]);
    }
    // Each target and number of modes, with the end of the windows from the target in which
    // the list is complete.
    for (const auto &[target, n, end] : std::vector<std::tuple<double, int, double>>{
             {5.0, 6, 9.988}, {9.7, 8, 9.988}, {10.1, 3, 10.25}, {10.5, 3, 13.62}})
    {
      CAPTURE(target, n);
      const auto found = Modes(RunEigenmode(SlabConfig(mesh, n, target, backend), temp_dir,
                                            backend + "-" + std::to_string(n)),
                               n);
      std::vector<std::complex<double>> nearest;
      std::ranges::copy_if(list, std::back_inserter(nearest),
                           [target](auto f) { return f.real() >= target; });
      REQUIRE(nearest.size() >= static_cast<std::size_t>(n));
      REQUIRE(nearest[n - 1].real() < end);
      for (int k = 0; k < n; k++)
      {
        CAPTURE(k, found[k], nearest[k]);
        CHECK(std::abs(found[k] - nearest[k]) <= 1.0e-7 * std::abs(nearest[k]));
      }
    }
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
  // eigensolver). The check is that every root up to the last mode found is among the modes
  // found, without the multiplicity of a root: as for the Lorentz-filled box, the linear
  // eigensolver of the seed can return one mode only of an exactly degenerate pair.
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
  for (const auto e : expected)
  {
    if (e.real() > found[n - 1].real())
    {
      break;
    }
    CAPTURE(e, found[n - 1]);
    CHECK(std::ranges::any_of(found, [e](auto f)
                              { return std::abs(f - e) <= 1.0e-7 * std::abs(e); }));
  }
  for (std::size_t k = 0; k < found_all.f.size(); k++)
  {
    CAPTURE(found_all.f[k], found_all.backward_error[k]);
    CHECK(found_all.backward_error[k] <= 1.0e-8);
    CHECK(found_all.f[k].imag() > 0.0);  // Lossy
  }
}

}  // namespace palace
