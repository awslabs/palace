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
#include <variant>
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

// The eigensolver backends of the build.
std::vector<std::string> Backends()
{
  std::vector<std::string> backends;
#if defined(PALACE_WITH_SLEPC)
  backends.push_back("SLEPc");
#endif
#if defined(PALACE_WITH_ARPACK)
  backends.push_back("ARPACK");
#endif
  return backends;
}

// Boxes with a fixed spectrum: second-order absorbing boundaries on the faces x = 0 and
// x = 1 next to a constant (ε = ε∞), a Debye (ε∞ = 2, Δε = 0.5, τ = 0.01 ns), or the
// Lorentz medium, the others PEC; or a finite conductivity on the face z = 1, the others
// PEC, next to the constant medium. The targets are 8 GHz, and 18 GHz for Debye.
enum class FixedBox
{
  Constant,
  Debye,
  Lorentz,
  Conductivity
};

json FixedBoxBoundaries(FixedBox fixture)
{
  if (fixture == FixedBox::Conductivity)
  {
    return {{"PEC", {{"Attributes", {1, 2, 3, 4, 5}}}},
            {"Conductivity", {{{"Attributes", {6}}, {"Conductivity", 1.0e5}}}}};
  }
  return {{"PEC", {{"Attributes", {1, 2, 4, 6}}}},
          {"Absorbing", {{"Attributes", {3, 5}}, {"Order", 2}}}};
}

json FixedBoxPermittivity(FixedBox fixture)
{
  switch (fixture)
  {
    case FixedBox::Debye:
      return {
          {"HighFrequency", eps_inf},
          {"Terms",
           {{{"Type", "Debye"}, {"DeltaPermittivity", 0.5}, {"RelaxationTime", 0.01}}}}};
    case FixedBox::Lorentz:
      return LorentzPermittivityConfig();
    default:
      return eps_inf;
  }
}

double FixedBoxTarget(FixedBox fixture)
{
  return (fixture == FixedBox::Debye) ? 18.0 : 8.0;
}

// The n lowest eigenfrequencies (GHz) above the target of a fixed box, from the complete
// spectrum of a dense companion problem on one process, independent of the eigensolvers.
// Each is checked on the complete operator T(λ) of SpaceOperator. For the dispersive
// media, each eigenvalue μ of the constant-medium companion problem with the absorbing
// boundary scalar of λ r(λ), r = √(ε(λ) / ε∞), gives the roots λ of the scalar relation
// λ r(λ) = μ on the sheet of the operator. For the conductivity, the companion problem is
// in z = √λ, and λ = z².
std::vector<std::complex<double>> FixedBoxSpectrum(const fs::path &temp_dir,
                                                   FixedBox fixture, int n)
{
  using namespace std::complex_literals;
  const auto comm = Mpi::World();
  REQUIRE(Mpi::Size(comm) == 1);
  test::IntegrationSettingsGuard settings_guard;
  const bool is_debye = fixture == FixedBox::Debye;
  const bool is_lorentz = fixture == FixedBox::Lorentz;
  const bool is_conductivity = fixture == FixedBox::Conductivity;
  const bool is_constant = fixture == FixedBox::Constant;
  const auto mesh_path = WriteBoxMesh(temp_dir);
  // The configuration names a backend of the build, which this function does not use.
  IoData iodata(BoxConfig(mesh_path, FixedBoxPermittivity(fixture),
                          FixedBoxBoundaries(fixture), n, FixedBoxTarget(fixture), "",
                          Backends().front()),
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
  REQUIRE(M);
  // The first-order absorbing boundary term is in C next to the constant medium, and in A2
  // next to the dispersive media.
  REQUIRE(static_cast<bool>(C) == is_constant);
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
    surface_mass = space_op.GetSurfaceConductivityBoundaryMatrix<ComplexOperator>(
        0, Operator::DIAG_ZERO);
    REQUIRE(surface_mass);
  }
  else if (is_constant)
  {
    curlcurl =
        space_op.GetFarfieldBoundaryCurlCurlMatrix<ComplexOperator>(Operator::DIAG_ZERO);
    REQUIRE(curlcurl);
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
  const Eigen::MatrixXcd Kd = Extract(*K), Md = Extract(*M);
  const Eigen::MatrixXcd Cd = C ? Extract(*C) : Eigen::MatrixXcd::Zero(dim, dim);
  const Eigen::LLT<Eigen::MatrixXcd> mass(Md);
  REQUIRE(mass.info() == Eigen::Success);
  Eigen::MatrixXcd companion;
  if (is_conductivity)
  {
    // The owner returns iB. For h = 0, δ = √(2 / (μ σ ω)), Z = (1 + i) / (σ δ), so the
    // scalar is κ √λ in Im λ > 0. P(z) = M z⁴ + κ B z + K, λ = z², state (u, z u, z² u,
    // z³ u). M is SPD: no clearing, poles or infinity.
    const Eigen::MatrixXcd B = Extract(*surface_mass) / 1i;
    companion = Eigen::MatrixXcd::Zero(4 * dim, 4 * dim);
    companion.block(0, dim, dim, dim).setIdentity();
    companion.block(dim, 2 * dim, dim, dim).setIdentity();
    companion.block(2 * dim, 3 * dim, dim, dim).setIdentity();
    companion.block(3 * dim, 0, dim, dim) = -mass.solve(Kd);
    companion.block(3 * dim, dim, dim, dim) = -kappa * mass.solve(B);
  }
  else
  {
    // T(μ) = K + μ C + μ² M - 0.5 / μ M_ff, with the damping matrix C of the first-order
    // term (the owner matrices of the dispersive media carry i). Clear its pole, then use
    // the strong cubic companion with state (u, μ u, μ² u). No inverse or subspace.
    const Eigen::MatrixXcd Cinf =
        is_constant ? Cd : Eigen::MatrixXcd(Extract(*damping) / 1i);
    const Eigen::MatrixXcd Mff =
        is_constant ? Extract(*curlcurl) : Eigen::MatrixXcd(Extract(*curlcurl) / 1i);
    companion = Eigen::MatrixXcd::Zero(3 * dim, 3 * dim);
    companion.block(0, dim, dim, dim).setIdentity();
    companion.block(dim, 2 * dim, dim, dim).setIdentity();
    companion.block(2 * dim, 0, dim, dim) = 0.5 * mass.solve(Mff);
    companion.block(2 * dim, dim, dim, dim) = -mass.solve(Kd);
    companion.block(2 * dim, 2 * dim, dim, dim) = -mass.solve(Cinf);
  }
  REQUIRE(companion.allFinite());
  const Eigen::ComplexEigenSolver<Eigen::MatrixXcd> spectrum(companion);
  REQUIRE(spectrum.info() == Eigen::Success);
  REQUIRE(spectrum.eigenvalues().allFinite());
  REQUIRE(spectrum.eigenvectors().allFinite());

  // The roots λ of each companion eigenvalue μ (z for the conductivity) above the target,
  // on the sheet of the operator, lowest Im λ (Re f) first. Every companion copy and every
  // mapped branch is retained: no field or frequency deduplication. For the constant
  // medium, there is none below the target other than the gradients at λ = 0 (a multiple
  // root, which round-off scatters to |λ| of about 1e-5 times the target here).
  const double a = is_debye ? 0.25 : 0.05;
  double tau = 0.0, omega0_squared = 0.0, gamma = 0.0;
  if (is_debye || is_lorentz)
  {
    const auto &term = iodata.domains.materials.front().permittivity_terms.front();
    if (is_debye)
    {
      const auto &law = std::get<config::DebyeData>(term);
      tau = law.tau;  // Already nondimensionalized, not the input ns value.
      REQUIRE(law.delta == 0.5);
    }
    else
    {
      const auto &law = std::get<config::LorentzData>(term);
      omega0_squared = law.omega_0 * law.omega_0;
      gamma = law.gamma;
      REQUIRE(law.delta == delta_eps);
    }
  }
  struct Candidate
  {
    std::complex<double> lambda;
    int mode;
  };
  std::vector<Candidate> candidates;
  int below = 0;
  for (int mode = 0; mode < spectrum.eigenvalues().size(); mode++)
  {
    const auto mu = spectrum.eigenvalues()(mode);
    Eigen::VectorXcd values;
    if (is_constant)
    {
      values = Eigen::VectorXcd::Constant(1, mu);
    }
    else if (is_conductivity)
    {
      values = Eigen::VectorXcd::Constant(1, mu * mu);  // λ = z²
    }
    else if (is_debye)
    {
      // μ = λ r, r² = 1 + a / (1 + τ λ): all three scalar roots.
      Eigen::Matrix3cd scalar = Eigen::Matrix3cd::Zero();
      scalar(0, 1) = scalar(1, 2) = 1.0;
      scalar(2, 0) = mu * mu / tau;
      scalar(2, 1) = mu * mu;
      scalar(2, 2) = -(1.0 + a) / tau;
      const Eigen::ComplexEigenSolver<Eigen::Matrix3cd> roots(scalar, false);
      REQUIRE(roots.info() == Eigen::Success);
      values = roots.eigenvalues();
    }
    else
    {
      // D(λ) (λ² - μ²) + a ω₀² λ² = 0, D = λ² + γ λ + ω₀². The quartic is monic.
      Eigen::Matrix4cd scalar = Eigen::Matrix4cd::Zero();
      scalar(0, 1) = scalar(1, 2) = scalar(2, 3) = 1.0;
      scalar(3, 0) = omega0_squared * mu * mu;
      scalar(3, 1) = gamma * mu * mu;
      scalar(3, 2) = mu * mu - (1.0 + a) * omega0_squared;
      scalar(3, 3) = -gamma;
      const Eigen::ComplexEigenSolver<Eigen::Matrix4cd> roots(scalar, false);
      REQUIRE(roots.info() == Eigen::Success);
      values = roots.eigenvalues();
    }
    REQUIRE(values.allFinite());
    for (const auto lambda : values)
    {
      if (lambda.imag() <= iodata.solver.eigenmode.target)
      {
        below += (is_constant && lambda.imag() > 1.0e-2 * iodata.solver.eigenmode.target);
        continue;
      }
      if (!is_constant)
      {
        // The analytic sheet of the operator's scalar, not a distance cutoff.
        const auto owner_value =
            is_conductivity
                ? space_op.GetSurfaceConductivityOp().EvaluateScalar(0, lambda / 1i)
                : lambda * space_op.GetFarfieldOp().GetDispersiveScale(0, lambda / 1i);
        const auto model_value = is_conductivity ? kappa * mu : mu;
        if (std::real(std::conj(model_value) * owner_value) < 0.0)
        {
          continue;
        }
      }
      candidates.push_back({lambda, mode});
    }
  }
  CHECK(below == 0);
  std::ranges::stable_sort(candidates, [](const auto &l, const auto &r)
                           { return l.lambda.imag() < r.lambda.imag(); });
  REQUIRE(candidates.size() >= static_cast<std::size_t>(n));

  // Native backward normalization: ||T(λ) x|| / (||x|| (||K|| + |λ| ||C|| + |λ|² ||M||)).
  // The eliminated K block is identity, hence max(1, ||K_free||).
  auto Norm = [](const Eigen::MatrixXcd &A)
  {
    const Eigen::JacobiSVD<Eigen::MatrixXcd> svd(A);
    REQUIRE(svd.info() == Eigen::Success);
    return svd.singularValues().maxCoeff();
  };
  const double normK = std::max(1.0, Norm(Kd)), normC = Norm(Cd), normM = Norm(Md);
  std::vector<std::complex<double>> expected;
  for (int k = 0; k < n; k++)
  {
    // The lowest n are selected before the validation: a failed candidate stops the case,
    // no higher candidate substitutes.
    const auto [lambda, mode] = candidates[k];
    CAPTURE(k, mode, lambda);
    Eigen::VectorXcd field = spectrum.eigenvectors().col(mode).head(dim);
    const double field_norm = field.norm();
    REQUIRE(std::isfinite(field_norm));
    REQUIRE(field_norm > 0.0);
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
    const double backward =
        action / (normK + std::abs(lambda) * normC + std::norm(lambda) * normM);
    const auto frequency =
        iodata.units.Dimensionalize<Units::ValueType::FREQUENCY>(lambda / 1i) /
        (2.0 * M_PI);
    Mpi::Print(
        "Fixed-box spectrum mode {:d}: f = {:.12e}{:+.12e}i GHz, backward = {:.3e}\n",
        k + 1, frequency.real(), frequency.imag(), backward);
    REQUIRE(backward <= 1.0e-8);
    expected.push_back(frequency);
  }
  return expected;
}

// The lowest modes above the target of the fixed boxes (GHz), from FixedBoxSpectrum.
const std::vector<std::complex<double>> constant_abc2_modes = {
    {13.33373469207, 2.199734381129}, {14.22474477916, 9.077464578857},
    {17.19541109405, 11.72289377094}, {18.36358844556, 1.956428758461},
    {20.36005578510, 8.643387453640}, {20.84330195020, 11.46050049783},
    {23.07797431387, 1.617380220753}, {24.05954758472, 10.13988796594},
    {24.08380007773, 1.199616515463}, {24.94974214409, 11.40373997771},
    {25.66438125839, 9.413930030888}, {26.59336710685, 10.64146610520},
    {27.25680784376, 6.628162295503}, {27.69726264615, 1.315861288834},
    {29.12393495954, 5.066430628881}, {29.34965173756, 11.44100635637},
    {30.41173891639, 11.78279664698}, {31.55163281795, 4.876212150301},
    {33.07768200068, 1.037107381349}, {33.55149993009, 0.6518134200846}};

const std::vector<std::complex<double>> debye_abc2_modes = {
    {18.78468849275, 10.05621002136}, {19.15947802623, 13.11310664957},
    {21.99149984087, 2.914592596579}, {22.64959793205, 11.80948177085},
    {23.03491379891, 2.533598861728}, {23.56024034526, 13.17004899462},
    {24.3610299302, 11.08994663097},  {25.30774245855, 12.39958516034},
    {26.09860497143, 8.230931140732}, {26.71499449196, 2.78497533712},
    {28.06820103421, 6.669399914204}, {28.1772488368, 13.27729093113},
    {29.27861665967, 13.6439073093},  {30.56652584218, 6.529912233014},
    {32.20594353367, 2.6393389868},   {32.69602397393, 2.256953880317},
    {33.47395494536, 10.4834215293},  {33.73151430188, 13.44239894606},
    {34.21908439059, 14.28300183515}, {35.7350850681, 5.442231528656}};
const std::vector<std::complex<double>> lorentz_abc2_modes = {
    {12.46381247035, 1.451488651821},  {13.61113256125, 0.386881564507},
    {13.80584778635, 0.2800203047645}, {13.8261237978, 0.2657585107049},
    {13.88060535601, 0.259063120424},  {13.90480581917, 0.3091485605215},
    {13.90682894278, 0.2871822102945}, {13.91629411739, 0.4145898481345},
    {13.92054511244, 0.2763882939216}, {13.92230478494, 0.2524306892286},
    {13.92390065601, 0.3348785506493}, {13.92473971357, 0.2502573340685},
    {13.93049773057, 0.3526688139333}, {13.93624274171, 0.328131721448},
    {13.93827833162, 0.2508939527882}, {13.94014257405, 0.261295375264},
    {13.94285904993, 0.3425550514265}, {13.9438587734, 0.2817881307931},
    {13.94453680484, 0.2585689822111}, {13.94485496478, 0.310042540685}};
const std::vector<std::complex<double>> conductivity_modes = {
    {17.02214338926, 0.008699223002311}, {20.6642492049, 0.01433407222388},
    {22.15176051954, 0.01291535754073},  {24.56454562049, 0.02128372813413},
    {24.5747616102, 0.01105397630204},   {25.13254121752, 0.01062400027702},
    {27.73476928962, 0.009220599844565}, {28.70493787589, 0.01138231435152},
    {30.74286211977, 0.02382730150408},  {30.75876804998, 0.007895542114848},
    {32.02135897455, 0.007433465205159}, {33.72463234975, 0.02496482159314},
    {33.74269154441, 0.006872295052889}, {34.14537826432, 0.0124649416367},
    {37.09635117434, 0.02437944684332},  {37.94507460581, 0.02356658003522},
    {38.45919258548, 0.02667483233399},  {38.48017616205, 0.005643742636361},
    {39.40365887558, 0.02737660987959},  {39.40877679638, 0.02226670221516}};

const std::vector<std::complex<double>> &FixedBoxModes(FixedBox fixture)
{
  switch (fixture)
  {
    case FixedBox::Debye:
      return debye_abc2_modes;
    case FixedBox::Lorentz:
      return lorentz_abc2_modes;
    case FixedBox::Conductivity:
      return conductivity_modes;
    default:
      return constant_abc2_modes;
  }
}

// Solve a fixed box for n modes with an eigensolver backend, and check that it returns the
// requested modes without gaps. They are those above the target which SLEPc and ARPACK
// select for TARGET_IMAGINARY: for SLEPc the n lowest in Re f, and for ARPACK those with
// the largest -Im θ of the shift-and-invert eigenvalue θ = 1 / (λ - σ), which prefers modes
// with less damping.
void CheckFixedBoxModes(const fs::path &temp_dir, FixedBox fixture,
                        const std::string &backend, int n,
                        const json &linear = json::object())
{
  const auto mesh = WriteBoxMesh(temp_dir);
  const double target = FixedBoxTarget(fixture);
  auto config = BoxConfig(mesh, FixedBoxPermittivity(fixture), FixedBoxBoundaries(fixture),
                          n, target, "", backend);
  config["Solver"]["Linear"].update(linear);
  const auto found_all = RunEigenmode(config, temp_dir, "fixed-box");
  REQUIRE(found_all.f.size() >= static_cast<std::size_t>(n));
  const auto &modes = FixedBoxModes(fixture);
  auto expected = modes;
  REQUIRE(expected.size() >= static_cast<std::size_t>(n));
  if (backend == "ARPACK")
  {
    // -Im θ = b / (a² + b²) for λ - σ ∝ a + i b, a = -Im f, b = Re f - target. No mode
    // beyond the list (with a larger b) can be among the requested: -Im θ ≤ 1 / b.
    auto key = [target](std::complex<double> f)
    { return -std::imag(1.0 / std::complex<double>(-f.imag(), f.real() - target)); };
    std::ranges::stable_sort(expected, [&](auto a, auto b) { return key(a) > key(b); });
    REQUIRE(1.0 / (modes.back().real() - target) < key(expected[n - 1]));
  }
  expected.resize(n);
  auto by_frequency = [](auto a, auto b) { return a.real() < b.real(); };
  std::ranges::sort(expected, by_frequency);
  auto found = found_all.f;
  std::ranges::sort(found, by_frequency);
  for (int k = 0; k < n; k++)
  {
    CAPTURE(k, found[k], expected[k]);
    CHECK(std::abs(found[k] - expected[k]) <= 1.0e-7 * std::abs(expected[k]));
  }
  for (std::size_t k = 0; k < found_all.f.size(); k++)
  {
    CAPTURE(found_all.f[k], found_all.backward_error[k]);
    CHECK(found_all.backward_error[k] <= 1.0e-8);
    CHECK(found_all.f[k].imag() > 0.0);  // Lossy or radiating
  }
}

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
                 "Eigenmode driver: warning of the longitudinal modes at a zero of the "
                 "permittivity",
                 "[eigenmodedriver][Serial][Parallel]")
{
  // At the zero of ε(f) of the Lorentz medium which fills the PEC box, the root of
  // f² − iγf − f₀² (1 + Δε/ε∞), every gradient field which vanishes on the walls solves the
  // equations without a magnetic field. From a target just below the zero, the solve
  // returns some of these fields first (the eigenvalue solver finds only some vectors of a
  // degenerate eigenvalue), and the warning lists exactly the modes at the zero. From 8
  // GHz, below the resonance, it returns electromagnetic modes and does not warn.
  const std::complex<double> i(0.0, 1.0);
  const std::complex<double> f_zero =
      0.5 * (i * damping +
             std::sqrt(4.0 * f0 * f0 * (1.0 + delta_eps / eps_inf) - damping * damping));
  const auto mesh = WriteBoxMesh(temp_dir);
  const json pec = {{"PEC", {{"Attributes", {1, 2, 3, 4, 5, 6}}}}};
  for (const double target : {8.0, 14.3})
  {
    CAPTURE(target, f_zero);
    REQUIRE(target < f_zero.real());
    Eigenvalues modes;
    const auto out = test::CaptureStdout(
        [&]
        {
          modes = RunEigenmode(
              BoxConfig(mesh, LorentzPermittivityConfig(), pec, 4, target, "", "Default"),
              temp_dir, "target_" + std::to_string(static_cast<int>(target)));
        });
    const bool longitudinal = (target > f0);
    REQUIRE(modes.f.size() >= 4);
    std::vector<std::size_t> listed;
    for (std::size_t k = 0; k < modes.f.size(); k++)
    {
      CAPTURE(k, modes.f[k]);
      const bool at_zero = std::abs(modes.f[k] - f_zero) <= 1.0e-6 * std::abs(f_zero);
      if (k == 0 || !longitudinal)  // The first mode, or any from 8 GHz
      {
        CHECK(at_zero == longitudinal);
      }
      if (at_zero)
      {
        listed.push_back(k + 1);
      }
    }
    if (Mpi::Root(Mpi::World()))
    {
      std::string warning = (listed.size() > 1) ? "Eigenmodes " : "Eigenmode ";
      for (std::size_t k = 0; k < listed.size(); k++)
      {
        warning += (k > 0 ? ", " : "") + std::to_string(listed[k]);
      }
      warning += (listed.size() > 1) ? " have" : " has";
      warning += " a magnetic energy below 1e-08 of the electric energy!";
      CAPTURE(out, warning);
      CHECK((out.find("magnetic energy below") != std::string::npos) == longitudinal);
      CHECK((out.find(warning) != std::string::npos) == longitudinal);
    }
  }
}

TEST_CASE_METHOD(test::SharedTempDir,
                 "Eigenmode driver: warning of a resonance near the search window of the "
                 "nonlinear eigensolver",
                 "[eigenmodedriver][Serial][Parallel]")
{
  // A weak Lorentz term in the PEC box, from 8 to 20 GHz (Target to TargetUpper). With a
  // Djordjevic-Sarkar term in the same material, the problem goes to the nonlinear
  // eigensolver, which warns of a resonance within 3 half-widths of the window: at 21 GHz,
  // 2 half-widths above it, and not at 22 GHz, 4 half-widths above it. Without the
  // Djordjevic-Sarkar term the pole terms are linearized exactly, and a resonance inside
  // the window is not warned of.
  constexpr double gamma = 1.0;  // GHz, the half-width is 0.5 GHz
  const auto mesh = WriteBoxMesh(temp_dir);
  const json pec = {{"PEC", {{"Attributes", {1, 2, 3, 4, 5, 6}}}}};
  const json ds = {{"Type", "DjordjevicSarkar"},
                   {"Strength", 0.05},
                   {"LowerFrequency", 0.1},
                   {"UpperFrequency", 100.0}};
  for (const auto &[f_res, with_ds, warned] : std::vector<std::tuple<double, bool, bool>>{
           {21.0, true, true}, {22.0, true, false}, {14.0, false, false}})
  {
    CAPTURE(f_res, with_ds);
    // The pole is at Im p = √(f₀² − γ²/4).
    json terms = {{{"Type", "Lorentz"},
                   {"DeltaPermittivity", 1.0e-3},
                   {"ResonanceFrequency", std::sqrt(f_res * f_res + 0.25 * gamma * gamma)},
                   {"DampingFrequency", gamma}}};
    if (with_ds)
    {
      terms.push_back(ds);
    }
    auto config = BoxConfig(mesh, {{"HighFrequency", eps_inf}, {"Terms", terms}}, pec, 2,
                            8.0, "", "Default");
    config["Solver"]["Eigenmode"]["TargetUpper"] = 20.0;
    // A name, not the structured binding, in the lambda (as in the slab cavity test).
    const auto name = "resonance_" + std::to_string(static_cast<int>(f_res));
    Eigenvalues modes;
    const auto out =
        test::CaptureStdout([&] { modes = RunEigenmode(config, temp_dir, name); });
    CHECK(modes.f.size() >= 2);
    if (Mpi::Root(Mpi::World()))
    {
      CAPTURE(out);
      CHECK((out.find("linearized exactly") == std::string::npos) == with_ds);
      CHECK((out.find("has a resonance at") != std::string::npos) == warned);
      CHECK((out.find("The frequency-dependent Permittivity of material 1 has a resonance "
                      "at 21 GHz, within 3 half-widths of the search window from Target "
                      "to TargetUpper, 8 to 20 GHz!") != std::string::npos) == warned);
    }
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
      // Copies, as clang cannot capture a structured binding in a lambda with OpenMP.
      auto inside = [a = a, b = b, h = h](auto f)
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
                           [target = target](auto f) { return f.real() >= target; });
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

// Frequency-dependent boundary terms: rational in the frequency (the second-order
// absorbing boundary next to a constant medium), or fitted with rational functions.

TEST_CASE_METHOD(test::SharedTempDir,
                 "Eigenmode driver: second-order absorbing boundaries next to a constant "
                 "medium",
                 "[eigenmodedriver][Serial][Parallel]")
{
  // The term -0.5 / λ M_ff of the second-order absorbing boundary is rational, and each
  // eigensolver backend of the build linearizes it exactly: it must return the requested
  // modes without gaps, for each number of modes. The Hybrid seed and Quasi-Newton
  // refinement skipped some of these, differently for each N.
  const std::string backend = GENERATE_REF(from_range(Backends()));
  const int n = GENERATE(4, 6, 8);
  CAPTURE(backend, n);
  CheckFixedBoxModes(temp_dir, FixedBox::Constant, backend, n);
}

TEST_CASE_METHOD(test::SharedTempDir,
                 "Eigenmode driver: second-order absorbing boundaries next to a Debye "
                 "medium",
                 "[eigenmodedriver][Serial][Parallel]")
{
  // The scalars √(ε / ε∞) and its inverse of the absorbing boundary next to the medium
  // are fitted; the cut of the square root is outside the window of the search. The
  // linear solver is exact (ComplexCoarseSolve, one level): the default solver with a
  // SuperLU coarse solve returns wrong solutions for this box, while reporting
  // convergence.
  const std::string backend = GENERATE_REF(from_range(Backends()));
  CAPTURE(backend);
  CheckFixedBoxModes(temp_dir, FixedBox::Debye, backend, 4,
                     {{"ComplexCoarseSolve", true}, {"MGMaxLevels", 1}});
}

TEST_CASE_METHOD(test::SharedTempDir,
                 "Eigenmode driver: second-order absorbing boundaries next to a Lorentz "
                 "medium",
                 "[eigenmodedriver][Serial][Parallel]")
{
  // The target is below the resonance. The branch cut of √(ε / ε∞) runs between the
  // material pole and the permittivity zero, at 14.00 and 14.34 GHz, in the window of the
  // search; the four lowest modes are below it.
  const std::string backend = GENERATE_REF(from_range(Backends()));
  CAPTURE(backend);
  CheckFixedBoxModes(temp_dir, FixedBox::Lorentz, backend, 4);
}

TEST_CASE_METHOD(test::SharedTempDir,
                 "Eigenmode driver: finite conductivity on a face of a PEC box",
                 "[eigenmodedriver][Serial][Parallel]")
{
  // The skin-effect scalar √λ of the face z = 1 is fitted; its branch point λ = 0 is
  // outside the window of the search.
  const std::string backend = GENERATE_REF(from_range(Backends()));
  CAPTURE(backend);
  CheckFixedBoxModes(temp_dir, FixedBox::Conductivity, backend, 4);
}

TEST_CASE_METHOD(test::SharedTempDir,
                 "Eigenmode driver: modes next to the branch cut of a fitted boundary term",
                 "[eigenmodedriver][Serial][Parallel]")
{
  // In the Lorentz box, the rational fit of √(ε / ε∞) is not valid in a strip around its
  // branch cut, from the material pole to the permittivity zero, where the modes
  // accumulate. With six modes, the fifth and sixth lowest are in the strip: the solve
  // must name the strip in a warning, and each mode found in it, whose error is that of the
  // complete operator.
  const std::string backend = GENERATE_REF(from_range(Backends()));
  CAPTURE(backend);
  const auto mesh = WriteBoxMesh(temp_dir);
  const auto fixture = FixedBox::Lorentz;
  const auto config =
      BoxConfig(mesh, FixedBoxPermittivity(fixture), FixedBoxBoundaries(fixture), 6,
                FixedBoxTarget(fixture), "", backend);
  Eigenvalues found;
  const auto out =
      test::CaptureStdout([&] { found = RunEigenmode(config, temp_dir, "strip"); });
  if (!Mpi::Root(Mpi::World()))
  {
    return;  // Only the root process prints
  }
  // The strip around the cut, a segment at Im f = γ / 2 (in GHz, the damping frequency is
  // γ / 2π) between the pole and the zero of ε.
  const std::string strip_text = "Its rational fit is not valid within ";
  const auto pos = out.find(strip_text);
  REQUIRE(pos != std::string::npos);
  const double strip = std::stod(out.substr(pos + strip_text.size()));
  CAPTURE(strip);
  REQUIRE(strip > 0.0);
  const double im = 0.5 * damping;
  const double pole = std::sqrt(f0 * f0 - im * im),
               zero = std::sqrt(f0 * f0 * (1.0 + delta_eps / eps_inf) - im * im);
  int inside = 0;
  for (std::size_t k = 0; k < found.f.size(); k++)
  {
    const auto f = found.f[k];
    const double d =
        std::abs(f - std::complex<double>(std::clamp(f.real(), pole, zero), im));
    const bool warned =
        out.find(fmt::format("Eigenmode {:d} at ", k + 1)) != std::string::npos;
    CAPTURE(k, f, d, warned);
    if (d < 0.9 * strip)
    {
      inside++;
      CHECK(warned);
    }
    else if (d > 1.1 * strip)
    {
      CHECK_FALSE(warned);
    }
  }
  CHECK(inside >= 1);
}

TEST_CASE_METHOD(
    test::SharedTempDir,
    "Eigenmode driver: spectra of the boxes with frequency-dependent boundaries",
    "[eigenmodedriver][Serial]")
{
  // The reference modes of the driver tests above, from the complete spectrum of a dense
  // companion problem (independent of the eigensolvers).
  const auto fixture = GENERATE(FixedBox::Constant, FixedBox::Debye, FixedBox::Lorentz,
                                FixedBox::Conductivity);
  const auto &modes = FixedBoxModes(fixture);
  CAPTURE(static_cast<int>(fixture));
  const auto expected = FixedBoxSpectrum(temp_dir, fixture, 20);
  REQUIRE(modes.size() <= expected.size());
  for (std::size_t k = 0; k < modes.size(); k++)
  {
    CAPTURE(k, expected[k], modes[k]);
    CHECK(std::abs(expected[k] - modes[k]) <= 1.0e-10 * std::abs(expected[k]));
  }
}

}  // namespace palace
