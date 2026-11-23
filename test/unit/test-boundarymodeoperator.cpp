// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include <algorithm>
#include <array>
#include <cmath>
#include <complex>
#include <functional>
#include <limits>
#include <memory>
#include <numbers>
#include <string>
#include <tuple>
#include <utility>
#include <vector>
#include <mfem.hpp>
#include <catch2/catch_test_macros.hpp>
#include <nlohmann/json.hpp>
#include <catch2/generators/catch_generators.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>
#include "drivers/boundarymodesolver.hpp"

#include "fem/fespace.hpp"
#include "fem/mesh.hpp"
#include "fixtures.hpp"
#include "linalg/vector.hpp"
#include "models/boundarymodeoperator.hpp"
#include "models/materialoperator.hpp"
#include "models/modeeigensolver.hpp"
#include "models/modeoperatorassembly.hpp"
#include "models/surfaceterms.hpp"
#include "models/waveportoperator.hpp"
#include "test-helpers.hpp"
#include "utils/communication.hpp"
#include "utils/constants.hpp"
#include "utils/filesystem.hpp"
#include "utils/geodata.hpp"
#include "utils/iodata.hpp"
#include "utils/units.hpp"

namespace palace
{
using namespace Catch::Matchers;

namespace
{

// Solve for modes of a 2D rectangular waveguide cross-section using
// ModeEigenSolver. Returns eigenvalues as complex kn values.
// MakeCartesian2D boundary attributes: bottom=1, right=2, top=3, left=4.
struct ModeResult
{
  std::vector<std::complex<double>> kn;
  std::vector<std::complex<double>> final_kn;
  std::vector<double> reduced_backward_errors;
  std::vector<double> final_reduced_backward_errors;
  int num_converged;
  ModeEigenSolver::ReducedModelStats reduced_stats;
  ModeEigenSolver::ReducedModelStats first_target_stats;
  std::size_t reduced_basis_size = 0;
  std::size_t first_target_basis_size = 0;
  double reduced_tol = 0.0;
  int complex_exact_converged = -1;
  // Shift-and-invert target of the in-band solve that produced kn.
  double kn_target = 0.0;
};

ModeResult SolveRectangularModes(double width, double height, double freq_ghz,
                                 double epsilon_r, int order, int num_modes,
                                 const std::function<void(IoData &)> &configure_bcs,
                                 bool exercise_reduced_model = false,
                                 int reduced_evaluations = 1,
                                 std::size_t reduced_training_capacity = 16,
                                 std::vector<double> reduced_training_factors = {0.9, 1.1},
                                 double adaptive_tol = 1.0e-3, double eig_tol = 1.0e-8,
                                 double reduced_basis_limit_factor = 0.0)
{
  MPI_Comm comm = Mpi::World();
  Units units(1.0, 1.0);
  IoData iodata(units);
  iodata.model.Lc = 1.0;

  auto &material = iodata.domains.materials.emplace_back();
  material.attributes = {1};
  material.epsilon_r.s = {epsilon_r, epsilon_r, epsilon_r};

  // Default: PEC on all boundaries.
  iodata.boundaries.pec.attributes = {1, 2, 3, 4};

  // Let the caller configure specific BCs.
  configure_bcs(iodata);

  iodata.solver.order = order;
  iodata.solver.boundary_mode.freq = freq_ghz;
  iodata.solver.boundary_mode.n = num_modes;
  iodata.solver.boundary_mode.tol = eig_tol;
  iodata.solver.linear.tol = eig_tol;
  iodata.solver.linear.max_it = 200;

  auto serial_mesh = std::make_unique<mfem::Mesh>(
      mfem::Mesh::MakeCartesian2D(10, 5, mfem::Element::TRIANGLE, false, width, height));
  iodata.NondimensionalizeInputs(serial_mesh);
  auto par_mesh = std::make_unique<mfem::ParMesh>(comm, *serial_mesh);
  iodata.CheckConfiguration();
  Mesh palace_mesh(std::move(par_mesh));

  auto nd_fec = std::make_unique<mfem::ND_FECollection>(order, palace_mesh.Dimension());
  auto h1_fec = std::make_unique<mfem::H1_FECollection>(order, palace_mesh.Dimension());
  FiniteElementSpace nd_fespace(palace_mesh, nd_fec.get());
  FiniteElementSpace h1_fespace(palace_mesh, h1_fec.get());
  MaterialOperator mat_op(iodata, palace_mesh);

  SurfaceTerms surface_terms(iodata, mat_op, palace_mesh.Get());

  mfem::Array<int> nd_dbc_tdof_list, h1_dbc_tdof_list;
  {
    const auto &pmesh = palace_mesh.Get();
    int bdr_attr_max = pmesh.bdr_attributes.Size() ? pmesh.bdr_attributes.Max() : 0;
    auto dbc_marker = mesh::AttrToMarker(bdr_attr_max, iodata.boundaries.pec.attributes);
    nd_fespace.Get().GetEssentialTrueDofs(dbc_marker, nd_dbc_tdof_list);
    h1_fespace.Get().GetEssentialTrueDofs(dbc_marker, h1_dbc_tdof_list);
  }

  int nd_size = nd_fespace.GetTrueVSize();
  mfem::Array<int> dbc_tdof_list;
  dbc_tdof_list.Append(nd_dbc_tdof_list);
  for (int i = 0; i < h1_dbc_tdof_list.Size(); i++)
  {
    dbc_tdof_list.Append(nd_size + h1_dbc_tdof_list[i]);
  }

  double omega = 2.0 * std::numbers::pi *
                 iodata.units.Nondimensionalize<Units::ValueType::FREQUENCY>(freq_ghz);

  // ModeEigenSolver requires a positive Krylov subspace size (num_vec). Mirror the
  // formula used by IoData::CheckConfiguration for eigenmode.max_size.
  const int num_vec = std::max(2 * num_modes, num_modes + 15);
  ModeEigenSolver mode_solver(
      mat_op, nullptr, surface_terms, nd_fespace, h1_fespace, dbc_tdof_list, num_modes,
      num_vec, eig_tol, EigenvalueSolver::WhichType::LARGEST_REAL, iodata.solver.linear,
      iodata.solver.boundary_mode.type, 0, nd_fespace.GetComm());

  auto target_at = [&](std::complex<double> w)
  { return w.real() * std::sqrt(1.1 * mat_op.GetMaxMuEpsilon()); };
  auto solve_at = [&](std::complex<double> w)
  {
    const double kn_target = target_at(w);
    const double sigma = -kn_target * kn_target;
    return std::make_pair(mode_solver.Solve(w, sigma), sigma);
  };

  if (exercise_reduced_model)
  {
    mode_solver.ConfigureReducedModelTraining(reduced_training_capacity);
    for (double factor : reduced_training_factors)
    {
      solve_at(factor * omega);
    }
    mode_solver.EnableReducedModel(adaptive_tol, reduced_basis_limit_factor);
  }

  auto result = solve_at(omega).first;
  ModeResult out;
  out.num_converged = result.num_converged;
  out.kn_target = target_at(omega);
  for (int i = 0; i < result.num_converged; i++)
  {
    // Capture the first in-band evaluation, which is the reduced result under test.
    out.kn.push_back(mode_solver.GetPropagationConstant(i));
    if (exercise_reduced_model)
    {
      out.reduced_backward_errors.push_back(
          mode_solver.GetError(i, EigenvalueSolver::ErrorType::BACKWARD));
    }
  }

  if (exercise_reduced_model)
  {
    out.first_target_stats = mode_solver.GetReducedModelStats();
    out.first_target_basis_size = mode_solver.GetReducedBasisSize();
    // The final evaluation may converge a different number of modes than the first one,
    // for example a reduced solve after an exact fallback.
    int final_converged = result.num_converged;
    for (int i = 1; i < reduced_evaluations; i++)
    {
      final_converged = solve_at(omega).first.num_converged;
    }
    for (int i = 0; i < final_converged; i++)
    {
      out.final_kn.push_back(mode_solver.GetPropagationConstant(i));
      out.final_reduced_backward_errors.push_back(
          mode_solver.GetError(i, EigenvalueSolver::ErrorType::BACKWARD));
    }
    // Complex-frequency queries must bypass the real-axis reduced model. Issue one complex
    // query solely to verify stats after retaining the reduced real-frequency result above.
    out.complex_exact_converged =
        solve_at(std::complex<double>(omega, 1.0e-3 * omega)).first.num_converged;
  }
  out.reduced_stats = mode_solver.GetReducedModelStats();
  out.reduced_basis_size = mode_solver.GetReducedBasisSize();
  out.reduced_tol = mode_solver.GetReducedTolerance();
  return out;
}

// Solve for the fundamental mode of the rectangular waveguide cross-section (at 100 GHz,
// where it is the only propagating mode) through the 2D BoundaryModeOperator path used by
// the BoundaryMode driver, which uses p-multigrid preconditioning when mg_max_levels > 1.
std::vector<std::complex<double>>
SolveRectangularModesMultigrid(int mg_max_levels,
                               const std::function<void(IoData &)> &configure_bcs)
{
  constexpr double width = 1000.0, height = 500.0, freq_ghz = 100.0, eig_tol = 1.0e-8;
  constexpr int order = 2, num_modes = 1;
  MPI_Comm comm = Mpi::World();
  Units units(1.0, 1.0);
  IoData iodata(units);
  iodata.problem.type = ProblemType::BOUNDARYMODE;
  iodata.model.Lc = 1.0;

  auto &material = iodata.domains.materials.emplace_back();
  material.attributes = {1};
  material.epsilon_r.s = {4.0, 4.0, 4.0};
  iodata.boundaries.pec.attributes = {1, 2, 3, 4};
  configure_bcs(iodata);

  iodata.solver.order = order;
  iodata.solver.boundary_mode.freq = freq_ghz;
  iodata.solver.boundary_mode.n = num_modes;
  iodata.solver.boundary_mode.tol = eig_tol;
  iodata.solver.linear.tol = 1.0e-10;
  iodata.solver.linear.max_it = 200;
  iodata.solver.linear.mg_max_levels = mg_max_levels;

  auto serial_mesh = std::make_unique<mfem::Mesh>(
      mfem::Mesh::MakeCartesian2D(10, 5, mfem::Element::TRIANGLE, false, width, height));
  iodata.NondimensionalizeInputs(serial_mesh);
  std::vector<std::unique_ptr<Mesh>> mesh;
  mesh.push_back(
      std::make_unique<Mesh>(std::make_unique<mfem::ParMesh>(comm, *serial_mesh)));
  iodata.CheckConfiguration();

  MaterialOperator mat_op(iodata, *mesh.back());
  BoundaryModeOperator mode_op(iodata, mesh, mat_op);
  REQUIRE(mode_op.GetNDSpaceHierarchy().GetNumLevels() ==
          static_cast<std::size_t>(mg_max_levels));

  const int nd_size = mode_op.GetNDTrueVSize();
  mfem::Array<int> dbc_tdof_list;
  dbc_tdof_list.Append(mode_op.GetNDDbcTDofLists().back());
  for (auto tdof : mode_op.GetH1DbcTDofLists().back())
  {
    dbc_tdof_list.Append(nd_size + tdof);
  }

  const int num_vec = std::max(2 * num_modes, num_modes + 15);
  ModeEigenSolver mode_solver(mode_op, dbc_tdof_list, num_modes, num_vec, eig_tol,
                              EigenvalueSolver::WhichType::LARGEST_REAL,
                              iodata.solver.linear, iodata.solver.boundary_mode.type, 0);
  const double omega =
      2.0 * M_PI * iodata.units.Nondimensionalize<Units::ValueType::FREQUENCY>(freq_ghz);
  const double kn_target = omega * std::sqrt(1.1 * mat_op.GetMaxMuEpsilon());
  auto result = mode_solver.Solve(omega, -kn_target * kn_target);

  std::vector<std::complex<double>> kn;
  for (int i = 0; i < result.num_converged; i++)
  {
    kn.push_back(mode_solver.GetPropagationConstant(i));
  }
  return kn;
}

// Components of the 2D mode eigenproblem matrix for a W × H rectangle (MakeCartesian2D
// boundary attributes: bottom = 1, x = W side = 2, top = 3, x = 0 side = 4), with elements
// below y = H/4 of attribute 1 (ε = 4) and the rest of attribute 2 (vacuum). The elements
// are distributed over the processes in vertical strips, so that all processes but the
// last have no elements on the side x = W. Boundaries listed as PEC are essential.
class RectangleModeModel
{
public:
  static constexpr double W = 2.0, H = 1.0;

private:
  IoData iodata;
  std::unique_ptr<Mesh> mesh;
  std::unique_ptr<mfem::ND_FECollection> nd_fec;
  std::unique_ptr<mfem::H1_FECollection> h1_fec;
  std::unique_ptr<FiniteElementSpace> nd_fespace, h1_fespace;
  std::unique_ptr<MaterialOperator> mat_op;
  std::unique_ptr<SurfaceTerms> surface_terms;
  std::unique_ptr<mfem::HypreParMatrix> Atnr, Atni, Btnr, Bttr;
  mfem::Array<int> dbc_tdof_list;
  std::unique_ptr<mode_assembly::ModeOperatorModel> model;

public:
  RectangleModeModel(const std::function<void(IoData &)> &configure_bcs)
    : iodata(Units(1.0, 1.0))
  {
    MPI_Comm comm = Mpi::World();
    iodata.problem.type = ProblemType::BOUNDARYMODE;
    iodata.model.Lc = 1.0;
    auto &substrate = iodata.domains.materials.emplace_back();
    substrate.attributes = {1};
    substrate.epsilon_r.s = {4.0, 4.0, 4.0};
    auto &vacuum = iodata.domains.materials.emplace_back();
    vacuum.attributes = {2};
    configure_bcs(iodata);
    iodata.solver.order = 2;
    iodata.solver.boundary_mode.freq = 1.0;
    iodata.solver.boundary_mode.n = 1;

    constexpr int nx = 8, ny = 4;
    auto serial_mesh = std::make_unique<mfem::Mesh>(
        mfem::Mesh::MakeCartesian2D(nx, ny, mfem::Element::TRIANGLE, false, W, H));
    const int num_procs = Mpi::Size(comm);
    REQUIRE(num_procs <= nx);
    std::vector<int> partitioning(serial_mesh->GetNE());
    for (int i = 0; i < serial_mesh->GetNE(); i++)
    {
      mfem::Vector center;
      serial_mesh->GetElementCenter(i, center);
      serial_mesh->SetAttribute(i, (center(1) < 0.25 * H) ? 1 : 2);
      partitioning[i] =
          std::min(static_cast<int>(center(0) / W * num_procs), num_procs - 1);
    }
    serial_mesh->SetAttributes();
    iodata.NondimensionalizeInputs(serial_mesh);
    mesh = std::make_unique<Mesh>(
        std::make_unique<mfem::ParMesh>(comm, *serial_mesh, partitioning.data()));
    iodata.CheckConfiguration();

    nd_fec = std::make_unique<mfem::ND_FECollection>(iodata.solver.order, 2);
    h1_fec = std::make_unique<mfem::H1_FECollection>(iodata.solver.order, 2);
    nd_fespace = std::make_unique<FiniteElementSpace>(*mesh, nd_fec.get());
    h1_fespace = std::make_unique<FiniteElementSpace>(*mesh, h1_fec.get());
    mat_op = std::make_unique<MaterialOperator>(iodata, *mesh);
    surface_terms = std::make_unique<SurfaceTerms>(iodata, *mat_op, mesh->Get());

    std::tie(Atnr, Atni) = mode_assembly::AssembleAtn(*nd_fespace, *h1_fespace, *mat_op);
    Btnr.reset(Atnr->Transpose());
    *Btnr *= -1.0;
    Bttr = std::get<0>(mode_assembly::AssembleBtt(*nd_fespace, *mat_op));

    // Essential (PEC) true DOFs of the block system: ND, then H1 offset by the ND size.
    {
      const auto &pmesh = mesh->Get();
      auto dbc_marker =
          mesh::AttrToMarker(pmesh.bdr_attributes.Max(), iodata.boundaries.pec.attributes);
      mfem::Array<int> nd_dbc, h1_dbc;
      nd_fespace->Get().GetEssentialTrueDofs(dbc_marker, nd_dbc);
      h1_fespace->Get().GetEssentialTrueDofs(dbc_marker, h1_dbc);
      dbc_tdof_list = nd_dbc;
      for (auto tdof : h1_dbc)
      {
        dbc_tdof_list.Append(nd_fespace->GetTrueVSize() + tdof);
      }
    }
    model = std::make_unique<mode_assembly::ModeOperatorModel>(
        *nd_fespace, *h1_fespace, *mat_op, nullptr, *surface_terms, *Bttr, Atnr.get(),
        Atni.get(), Btnr.get(), dbc_tdof_list);
  }

  const auto &GetModel() const { return *model; }
  const auto &GetMaterialOp() const { return *mat_op; }
  const auto &GetDbcTDofList() const { return dbc_tdof_list; }
  const auto &GetUnits() const { return iodata.units; }
  const auto &GetSurfaceTerms() const { return *surface_terms; }

  const mode_assembly::OperatorComponent &
  GetComponent(mode_assembly::CoefficientType type) const
  {
    const auto &components = model->GetComponents();
    auto it =
        std::find_if(components.begin(), components.end(),
                     [type](const auto &component) { return component.type == type; });
    REQUIRE(it != components.end());
    return *it;
  }

  // Pairing xᴴ A x of a component A with x = [Eₜ; Eₙ] for the in-plane field Eₜ = eₜ ŷ,
  // tangential to the sides x = 0 and x = W, and a constant out-of-plane field Eₙ = eₙ.
  // Returns the pairing with the real and the imaginary part of the component.
  std::complex<double> Pairing(const mode_assembly::OperatorComponent &component, double et,
                               double en) const
  {
    MPI_Comm comm = Mpi::World();
    const int nd_size = nd_fespace->GetTrueVSize(), h1_size = h1_fespace->GetTrueVSize();
    ComplexVector x(nd_size + h1_size), y(nd_size + h1_size);
    x = 0.0;
    {
      mfem::ParGridFunction E(&nd_fespace->Get());
      mfem::Vector field(2);
      field(0) = 0.0;
      field(1) = et;
      mfem::VectorConstantCoefficient coeff(field);
      E.ProjectCoefficient(coeff);
      Vector t(nd_size);
      E.GetTrueDofs(t);
      x.Real().SetVector(t, 0);
    }
    auto *X = x.Real().HostReadWrite();
    for (int i = nd_size; i < nd_size + h1_size; i++)
    {
      X[i] = en;
    }
    component.op->Mult(x, y);
    return {linalg::Dot(comm, x.Real(), y.Real()), linalg::Dot(comm, x.Real(), y.Imag())};
  }
};
// Exact mode index of the TM0 mode of a parallel plate guide with a PEC plate at y = 0, a
// vacuum gap of height h, and a London slab of thickness d and penetration depth lambda,
// whose back face is PEC or free (PMC), at free-space wavenumber k0 (all lengths in um).
// With p^2 = k0^2 (n^2 - 1) and q^2 = 1/lambda^2 + p^2, H_x = cosh(p y) in the gap and cosh
// or sinh of q (h + d - y) in the slab, and continuity of H_x and of (1/eps) dH_x/dy at y =
// h, where eps = 1 - 1/(k0 lambda)^2 is the effective permittivity of the slab, give p
// tanh(p h) = -(q/eps) coth(q d) (free back face) or -(q/eps) tanh(q d) (PEC back face). In
// the quasi-static limit, n^2 = 1 + (lambda/h) coth(d/lambda) or tanh(d/lambda) (Swihart).
double ExactLondonSlabIndex(double h, double d, double lambda, double k0, bool pec_back)
{
  const double eps = 1.0 - 1.0 / (k0 * k0 * lambda * lambda);
  auto F = [&](double n)
  {
    const double p = k0 * std::sqrt(n * n - 1.0);
    const double q = std::sqrt(1.0 / (lambda * lambda) + p * p);
    const double g = pec_back ? std::tanh(q * d) : 1.0 / std::tanh(q * d);
    return p * std::tanh(p * h) + q * g / eps;
  };
  double lo = 1.0 + 1.0e-12, hi = 10.0;
  for (int it = 0; it < 200; it++)
  {
    const double mid = 0.5 * (lo + hi);
    ((F(mid) < 0.0) ? lo : hi) = mid;
  }
  return 0.5 * (lo + hi);
}

// Propagation constant normalized by the free-space wavenumber (n_eff) of a mode of a 2D
// cross-section of width w, in um: a PEC plate at y = 0 (attribute 1), a vacuum gap of
// height h, and a London superconductor slab of thickness d with penetration depth lambda
// above it, whose back face y = h + d (attribute 2) is PEC or left free (PMC). The sides x
// = 0, w (attribute 3) are PEC or free.
double SolveLondonSlabMode(double w, double h, double d, double lambda, bool pec_back,
                           bool pec_sides, double freq_ghz, double n_target, int order,
                           int ny_gap = 8)
{
  MPI_Comm comm = Mpi::World();
  Units units(1.0e-6, 1.0e-6);
  IoData iodata(units);
  iodata.model.Lc = 1.0;
  auto &vacuum = iodata.domains.materials.emplace_back();
  vacuum.attributes = {1};
  auto &london = iodata.domains.materials.emplace_back();
  london.attributes = {2};
  london.lambda_L = lambda;
  iodata.boundaries.pec.attributes = {1};
  if (pec_back)
  {
    iodata.boundaries.pec.attributes.push_back(2);
  }
  if (pec_sides)
  {
    iodata.boundaries.pec.attributes.push_back(3);
  }
  iodata.solver.order = order;
  iodata.solver.boundary_mode.freq = freq_ghz;
  iodata.solver.boundary_mode.n = 1;
  iodata.solver.boundary_mode.tol = 1.0e-12;
  iodata.solver.linear.tol = 1.0e-12;
  iodata.solver.linear.max_it = 200;

  // Structured quadrilaterals: 2 across the width, ny_gap across the gap, and enough across
  // the slab to resolve the penetration depth.
  constexpr int nx = 2;
  const int ny_slab = std::max(8, static_cast<int>(std::ceil(8.0 * d / lambda)));
  std::vector<double> y;
  for (int j = 0; j <= ny_gap; j++)
  {
    y.push_back(h * j / ny_gap);
  }
  for (int j = 1; j <= ny_slab; j++)
  {
    y.push_back(h + d * j / ny_slab);
  }
  const int ny = static_cast<int>(y.size()) - 1;
  auto serial_mesh =
      std::make_unique<mfem::Mesh>(2, (nx + 1) * (ny + 1), nx * ny, 2 * (nx + ny));
  for (int j = 0; j <= ny; j++)
  {
    for (int i = 0; i <= nx; i++)
    {
      serial_mesh->AddVertex(w * i / nx, y[j]);
    }
  }
  auto v = [&](int i, int j) { return i + (nx + 1) * j; };
  for (int j = 0; j < ny; j++)
  {
    for (int i = 0; i < nx; i++)
    {
      serial_mesh->AddQuad(v(i, j), v(i + 1, j), v(i + 1, j + 1), v(i, j + 1),
                           (j < ny_gap) ? 1 : 2);
    }
  }
  for (int i = 0; i < nx; i++)
  {
    serial_mesh->AddBdrSegment(v(i, 0), v(i + 1, 0), 1);
    serial_mesh->AddBdrSegment(v(i + 1, ny), v(i, ny), 2);
  }
  for (int j = 0; j < ny; j++)
  {
    serial_mesh->AddBdrSegment(v(0, j + 1), v(0, j), 3);
    serial_mesh->AddBdrSegment(v(nx, j), v(nx, j + 1), 3);
  }
  serial_mesh->FinalizeQuadMesh(1, 1, true);

  iodata.NondimensionalizeInputs(serial_mesh);
  auto par_mesh = std::make_unique<mfem::ParMesh>(comm, *serial_mesh);
  iodata.CheckConfiguration();
  Mesh palace_mesh(std::move(par_mesh));

  auto nd_fec = std::make_unique<mfem::ND_FECollection>(order, palace_mesh.Dimension());
  auto h1_fec = std::make_unique<mfem::H1_FECollection>(order, palace_mesh.Dimension());
  FiniteElementSpace nd_fespace(palace_mesh, nd_fec.get());
  FiniteElementSpace h1_fespace(palace_mesh, h1_fec.get());
  MaterialOperator mat_op(iodata, palace_mesh);
  SurfaceTerms surface_terms(iodata, mat_op, palace_mesh.Get());

  mfem::Array<int> nd_dbc_tdof_list, h1_dbc_tdof_list;
  {
    const auto &pmesh = palace_mesh.Get();
    int bdr_attr_max = pmesh.bdr_attributes.Size() ? pmesh.bdr_attributes.Max() : 0;
    auto dbc_marker = mesh::AttrToMarker(bdr_attr_max, iodata.boundaries.pec.attributes);
    nd_fespace.Get().GetEssentialTrueDofs(dbc_marker, nd_dbc_tdof_list);
    h1_fespace.Get().GetEssentialTrueDofs(dbc_marker, h1_dbc_tdof_list);
  }
  const int nd_size = nd_fespace.GetTrueVSize();
  mfem::Array<int> dbc_tdof_list;
  dbc_tdof_list.Append(nd_dbc_tdof_list);
  for (int i = 0; i < h1_dbc_tdof_list.Size(); i++)
  {
    dbc_tdof_list.Append(nd_size + h1_dbc_tdof_list[i]);
  }

  const double omega =
      2.0 * std::numbers::pi *
      iodata.units.Nondimensionalize<Units::ValueType::FREQUENCY>(freq_ghz);
  // Shift-and-invert about the target, as in BoundaryModeSolver with a target n_eff.
  const int num_vec = 16;
  ModeEigenSolver mode_solver(
      mat_op, nullptr, surface_terms, nd_fespace, h1_fespace, dbc_tdof_list, 1, num_vec,
      1.0e-12, EigenvalueSolver::WhichType::LARGEST_MAGNITUDE, iodata.solver.linear,
      iodata.solver.boundary_mode.type, 0, nd_fespace.GetComm());
  const double kn_target = omega * n_target;
  auto result = mode_solver.Solve(omega, -kn_target * kn_target);
  REQUIRE(result.num_converged >= 1);
  const auto kn = mode_solver.GetPropagationConstant(0);
  CHECK(std::abs(kn.imag()) < 1.0e-6 * std::abs(kn.real()));
  return kn.real() / omega;
}

// Effective index n_eff = kn / k0 (Im{n_eff} < 0) of the TM0 mode of a parallel plate line
// with PEC planes at y = 0 and y = t + d, a conductor of conductivity sigma (and relative
// permittivity 1) in 0 < y < t, and a dielectric of relative permittivity eps_d in t < y <
// t + d (lengths in um). With H_x = cos(kc y) in the conductor and cos(kd (t + d - y)) in
// the dielectric, continuity of H_x and of (1/eps) dH_x/dy at y = t gives (kc/eps_c)
// tan(kc t) + (kd/eps_d) tan(kd d) = 0, where kj^2 = k0^2 (eps_j - n_eff^2) and eps_c = 1 -
// i sigma / (omega eps0). The left-hand side is even in kc and kd, so the branches of the
// square roots do not matter.
std::complex<double> ExactParallelPlateConductorIndex(double t, double d, double sigma,
                                                      double eps_d, double freq_ghz)
{
  const double omega = 2.0 * std::numbers::pi * freq_ghz * 1.0e9;
  const double k0 = omega / electromagnetics::c0_ * 1.0e-6;
  const std::complex<double> eps_c(1.0, -sigma / (omega * electromagnetics::epsilon0_));
  auto F = [&](std::complex<double> n)
  {
    const std::complex<double> kc = k0 * std::sqrt(eps_c - n * n);
    const std::complex<double> kd = k0 * std::sqrt(eps_d - n * n);
    return kc / eps_c * std::tan(kc * t) + kd / eps_d * std::tan(kd * d);
  };
  // Secant iteration from the mode of the line with a perfect conductor.
  std::complex<double> n0 = std::sqrt(eps_d), n1 = n0 * std::complex<double>(1.0, -0.01);
  std::complex<double> F0 = F(n0), F1 = F(n1);
  for (int it = 0; it < 100 && std::abs(n1 - n0) > 1.0e-15 * std::abs(n1); it++)
  {
    const std::complex<double> n2 = n1 - F1 * (n1 - n0) / (F1 - F0);
    n0 = n1;
    F0 = F1;
    n1 = n2;
    F1 = F(n1);
  }
  return n1;
}

struct ParallelPlateModeResult
{
  std::complex<double> n_eff;
  int num_solves, num_iterations;
};

// Effective index of the TM0 mode of a parallel plate line of width 2 um, with layers of
// attribute 1 in 0 < y < t and 2 in t < y < t + d (configured by the caller, along with
// the boundaries: bottom 1, right 2, top 3, left 4), solved through the 2D
// BoundaryModeOperator path of the BoundaryMode driver, with a sparse direct or p-multigrid
// (mg_max_levels > 1) preconditioner of the given type, and the numbers of linear solves
// and iterations. The elements are 0.25 um.
ParallelPlateModeResult
SolveParallelPlateMode(double t, double d, double freq_ghz, double n_target,
                       int mg_max_levels, bool pc_mat_real, bool complex_coarse_solve,
                       const std::function<void(IoData &)> &configure)
{
  constexpr double width = 2.0, h = 0.25, tol = 1.0e-10;
  constexpr int order = 2;
  MPI_Comm comm = Mpi::World();
  Units units(1.0e-6, 1.0e-6);
  IoData iodata(units);
  iodata.problem.type = ProblemType::BOUNDARYMODE;
  iodata.model.Lc = 1.0;
  configure(iodata);
  iodata.solver.order = order;
  iodata.solver.boundary_mode.freq = freq_ghz;
  iodata.solver.boundary_mode.n = 1;
  iodata.solver.boundary_mode.tol = tol;
  iodata.solver.linear.tol = tol;
  iodata.solver.linear.max_it = 200;
  iodata.solver.linear.mg_max_levels = mg_max_levels;
  iodata.solver.linear.pc_mat_real = pc_mat_real;
  iodata.solver.linear.complex_coarse_solve = complex_coarse_solve;

  auto serial_mesh = std::make_unique<mfem::Mesh>(mfem::Mesh::MakeCartesian2D(
      static_cast<int>(std::round(width / h)), static_cast<int>(std::round((t + d) / h)),
      mfem::Element::TRIANGLE, false, width, t + d));
  for (int i = 0; i < serial_mesh->GetNE(); i++)
  {
    mfem::Vector center;
    serial_mesh->GetElementCenter(i, center);
    serial_mesh->SetAttribute(i, (center(1) < t) ? 1 : 2);
  }
  serial_mesh->SetAttributes();
  iodata.NondimensionalizeInputs(serial_mesh);
  std::vector<std::unique_ptr<Mesh>> mesh;
  mesh.push_back(
      std::make_unique<Mesh>(std::make_unique<mfem::ParMesh>(comm, *serial_mesh)));
  iodata.CheckConfiguration();

  MaterialOperator mat_op(iodata, *mesh.back());
  BoundaryModeOperator mode_op(iodata, mesh, mat_op);
  const int nd_size = mode_op.GetNDTrueVSize();
  mfem::Array<int> dbc_tdof_list;
  dbc_tdof_list.Append(mode_op.GetNDDbcTDofLists().back());
  for (auto tdof : mode_op.GetH1DbcTDofLists().back())
  {
    dbc_tdof_list.Append(nd_size + tdof);
  }
  ModeEigenSolver mode_solver(mode_op, dbc_tdof_list, 1, 16, tol,
                              EigenvalueSolver::WhichType::LARGEST_MAGNITUDE,
                              iodata.solver.linear, iodata.solver.boundary_mode.type, 0);
  const double omega =
      2.0 * std::numbers::pi *
      iodata.units.Nondimensionalize<Units::ValueType::FREQUENCY>(freq_ghz);
  const double kn_target = omega * n_target;
  auto result = mode_solver.Solve(omega, -kn_target * kn_target);
  REQUIRE(result.num_converged >= 1);
  const auto *ksp = mode_solver.GetLinearSolver();
  REQUIRE(ksp);
  return {mode_solver.GetPropagationConstant(0) / omega, ksp->NumTotalMult(),
          ksp->NumTotalMultIterations()};
}

// The two modes nearest to the shift target kₙ = omega n_target of the cross-section of
// test::ExactAbsorbingParallelPlateKn for W = 1 m and H = 0.5 m, at 0.5 GHz (the free-space
// wavenumber k = omega, in 1/m), with a mesh of n x n / 2 quadrilaterals of the given
// order, solved through the 2D BoundaryModeOperator path of the BoundaryMode driver with a
// sparse direct or p-multigrid (mg_max_levels > 1) preconditioner: their propagation
// constants, and the numbers of linear solves and iterations. MakeCartesian2D boundary
// attributes: bottom = 1, right = 2, top = 3, left = 4.
constexpr double absorbing_parallel_plate_freq_ghz = 0.5;
struct AbsorbingParallelPlateModes
{
  std::vector<std::complex<double>> kn;
  int num_solves, num_iterations;
};
AbsorbingParallelPlateModes SolveAbsorbingParallelPlateModes(int abc_order, double n_target,
                                                             int order, int n,
                                                             int mg_max_levels,
                                                             bool pc_mat_real)
{
  constexpr double W = 1.0, H = 0.5, freq_ghz = absorbing_parallel_plate_freq_ghz,
                   tol = 1.0e-10;
  constexpr int num_modes = 2;
  MPI_Comm comm = Mpi::World();
  IoData iodata(Units(1.0, 1.0));
  iodata.problem.type = ProblemType::BOUNDARYMODE;
  iodata.model.L0 = 1.0;
  iodata.model.Lc = 1.0;
  iodata.domains.materials.emplace_back().attributes = {1};
  iodata.boundaries.pec.attributes = {1, 3};
  iodata.boundaries.farfield.attributes = {2, 4};
  iodata.boundaries.farfield.order = abc_order;
  iodata.solver.order = order;
  iodata.solver.boundary_mode.freq = freq_ghz;
  iodata.solver.boundary_mode.n = num_modes;
  iodata.solver.boundary_mode.tol = tol;
  iodata.solver.linear.tol = tol;
  iodata.solver.linear.max_it = 200;
  iodata.solver.linear.mg_max_levels = mg_max_levels;
  iodata.solver.linear.pc_mat_real = pc_mat_real;

  auto serial_mesh = std::make_unique<mfem::Mesh>(
      mfem::Mesh::MakeCartesian2D(n, n / 2, mfem::Element::QUADRILATERAL, false, W, H));
  iodata.NondimensionalizeInputs(serial_mesh);
  std::vector<std::unique_ptr<Mesh>> mesh;
  mesh.push_back(
      std::make_unique<Mesh>(std::make_unique<mfem::ParMesh>(comm, *serial_mesh)));
  iodata.CheckConfiguration();

  MaterialOperator mat_op(iodata, *mesh.back());
  BoundaryModeOperator mode_op(iodata, mesh, mat_op);
  const int nd_size = mode_op.GetNDTrueVSize();
  mfem::Array<int> dbc_tdof_list;
  dbc_tdof_list.Append(mode_op.GetNDDbcTDofLists().back());
  for (auto tdof : mode_op.GetH1DbcTDofLists().back())
  {
    dbc_tdof_list.Append(nd_size + tdof);
  }
  ModeEigenSolver mode_solver(mode_op, dbc_tdof_list, num_modes, 16, tol,
                              EigenvalueSolver::WhichType::LARGEST_MAGNITUDE,
                              iodata.solver.linear, iodata.solver.boundary_mode.type, 0);
  const double omega =
      2.0 * std::numbers::pi *
      iodata.units.Nondimensionalize<Units::ValueType::FREQUENCY>(freq_ghz);
  const double kn_target = omega * n_target;
  auto result = mode_solver.Solve(omega, -kn_target * kn_target);
  REQUIRE(result.num_converged >= num_modes);
  AbsorbingParallelPlateModes modes;
  for (int i = 0; i < num_modes; i++)
  {
    modes.kn.push_back(mode_solver.GetPropagationConstant(i));
  }
  const auto *ksp = mode_solver.GetLinearSolver();
  REQUIRE(ksp);
  modes.num_solves = ksp->NumTotalMult();
  modes.num_iterations = ksp->NumTotalMultIterations();
  return modes;
}

// Whether a propagation constant is within the relative tolerance of the expected one.
bool HasMode(const std::vector<std::complex<double>> &kn, std::complex<double> kn_expected,
             double rtol)
{
  return std::ranges::any_of(
      kn,
      [&](auto k) { return std::abs(k - kn_expected) <= rtol * std::abs(kn_expected); });
}

}  // namespace

TEST_CASE("ModeEigenSolver PEC", "[boundarymodeoperator][Serial]")
{
  // Rectangular waveguide: 1000×500 μm (L0=1e-6), ε=4, f=500 GHz.
  // Analytical kn for TE10 mode:
  //   kc = π / (a * L0) = π / 1e-3 ≈ 3141.6 1/m
  //   ω = 2π * 500e9 ≈ 3.1416e12 rad/s
  //   kn = sqrt(ω²ε/c² - kc²) = sqrt(4*(π*1e12/c)² - (π/1e-3)²)
  //   In nondimensional units (Lc = a*L0 = 1e-3 m):
  //     kn_nd = kn * Lc
  auto result = SolveRectangularModes(1000.0, 500.0, 500.0, 4.0, 2, 3, [](IoData &) {});

  REQUIRE(result.num_converged >= 1);

  double kn_real = result.kn[0].real();
  double kn_imag = result.kn[0].imag();
  CAPTURE(kn_real, kn_imag, result.num_converged);

  // First mode should be propagating (real kn, negligible imaginary part).
  CHECK(kn_real > 0.0);
  CHECK(std::abs(kn_imag) < 1.0e-6 * std::abs(kn_real));

  // Analytical kn for TE10 mode of rectangular waveguide with PEC walls:
  //   a = 1000 μm = 1e-3 m, ε_r = 4, f = 500 GHz
  //   kc = π / a = π / 1e-3 m
  //   kn = sqrt(ω²ε_r/c² - kc²)
  //   kn ≈ 20708 1/m → nondimensional (×Lc where Lc = 1e-6 m) ≈ 0.02071
  // Allow 5% tolerance for the coarse 10×5 mesh at order 2.
  CHECK_THAT(kn_real, WithinRel(0.02071, 0.05));
}

TEST_CASE("ModeEigenSolver dispersive material in an oblique cross-section frame",
          "[boundarymodeoperator][Serial][Parallel]")
{
  // A port on an oblique plane rotates the materials into the frame of its cross-section
  // (RotateMaterialDefinitions). A scalar dispersive substrate gives the same mode as in
  // the axis-aligned frame, to the tolerance of the two eigenvalue solves (the default of
  // SolveRectangularModes); an anisotropic one is still rejected.
  constexpr double eig_tol = 1.0e-8;
  auto configure = [](const std::array<mfem::Vector, 3> *frame, double eps_zz)
  {
    return [frame, eps_zz](IoData &iodata)
    {
      const double rate = electromagnetics::c0_ / iodata.model.L0;
      auto &material = iodata.domains.materials[0];
      material.epsilon_r.s = {4.0, 4.0, eps_zz};
      material.permittivity_terms.push_back(config::DebyeData{0.5, 1.0 / rate});
      if (frame)
      {
        RotateMaterialDefinitions(iodata.domains.materials, (*frame)[0], (*frame)[1],
                                  (*frame)[2]);
      }
    };
  };
  const auto ref =
      SolveRectangularModes(1000.0, 500.0, 500.0, 4.0, 2, 1, configure(nullptr, 4.0));
  REQUIRE(ref.num_converged >= 1);
  for (const auto &frame : test::ObliqueFrames(2))
  {
    const auto rotated =
        SolveRectangularModes(1000.0, 500.0, 500.0, 4.0, 2, 1, configure(&frame, 4.0));
    REQUIRE(rotated.num_converged >= 1);
    CHECK(std::abs(rotated.kn[0] - ref.kn[0]) <= eig_tol * std::abs(ref.kn[0]));
    CHECK_THROWS(
        SolveRectangularModes(1000.0, 500.0, 500.0, 4.0, 2, 1, configure(&frame, 4.5)));
  }
}

TEST_CASE("ModeEigenSolver guarded reduced real-frequency solve",
          "[boundarymodeoperator][Serial][Parallel]")
{
  constexpr int num_modes = 3;
  auto exact =
      SolveRectangularModes(1000.0, 500.0, 500.0, 4.0, 2, num_modes, [](IoData &) {});
  REQUIRE(exact.num_converged >= num_modes);
  CHECK(exact.reduced_basis_size == 0);
  CHECK(exact.reduced_stats.reduced_solves == 0);
  CHECK(exact.reduced_stats.exact_solves == 1);

  // The mode shapes of the homogeneous waveguide do not depend on frequency, so the
  // eigenvectors of the second training solve only differ from the basis by their errors,
  // which are added to the basis if they exceed its orthogonalization tolerance. With the
  // default eigenvalue tolerance, these errors are close to that tolerance, so whether they
  // are added, which can make the first reduced solve fall back to an exact solve, depends
  // on roundoff. Solve the eigenproblems accurately to keep the errors well below it.
  auto reduced = SolveRectangularModes(
      1000.0, 500.0, 500.0, 4.0, 2, num_modes, [](IoData &) {}, true, 21, 16, {0.9, 1.1},
      1.0e-3, 1.0e-12);
  REQUIRE(reduced.num_converged >= num_modes);
  REQUIRE(reduced.reduced_basis_size >= num_modes);
  CHECK(reduced.reduced_stats.reduced_solves == 21);
  CHECK(reduced.reduced_stats.exact_solves == 2);
  CHECK(reduced.reduced_stats.worst_residual <= reduced.reduced_tol);
  REQUIRE(reduced.reduced_backward_errors.size() == num_modes);
  for (double error : reduced.reduced_backward_errors)
  {
    CHECK(error <= reduced.reduced_tol);
  }
  CHECK(reduced.complex_exact_converged >= num_modes);
  CHECK(reduced.reduced_stats.offline_basis_rank >= num_modes);
  CHECK(reduced.reduced_stats.online_basis_cap >=
        reduced.reduced_stats.offline_basis_rank + 4 * num_modes);
  for (int i = 0; i < num_modes; i++)
  {
    CHECK_THAT(reduced.kn[i].real(), WithinRel(exact.kn[i].real(), 1.0e-6));
    CHECK_THAT(reduced.kn[i].imag(), WithinAbs(exact.kn[i].imag(), 1.0e-8));
  }
}

// With a basis size limit, the reduced model is only used while its basis is within the
// limit, which also caps its online enrichment, and otherwise all solves are exact without
// enriching the basis. The decision only depends on the sizes of the basis and of the
// problem, so it is the same on all processes and in all runs.
TEST_CASE("ModeEigenSolver reduced model basis size limit",
          "[boundarymodeoperator][Serial][Parallel]")
{
  constexpr int num_modes = 3, num_evaluations = 4;
  auto exact =
      SolveRectangularModes(1000.0, 500.0, 500.0, 4.0, 2, num_modes, [](IoData &) {});
  REQUIRE(exact.num_converged >= num_modes);

  // The basis of the two training solves is within the limit for wave ports, of 2√N vectors
  // (55 vectors for the N = 761 unknowns of this cross-section). The eigenproblems are
  // solved accurately, as for the guarded reduced solve above.
  auto within = SolveRectangularModes(
      1000.0, 500.0, 500.0, 4.0, 2, num_modes, [](IoData &) {}, true, num_evaluations, 16,
      {0.9, 1.1}, 1.0e-3, 1.0e-12, 2.0);
  REQUIRE(within.num_converged >= num_modes);
  const auto &within_stats = within.reduced_stats;
  CHECK_FALSE(within_stats.over_basis_limit);
  CHECK(within_stats.basis_limit == 55);
  CHECK(within_stats.basis_limit >= within_stats.offline_basis_rank);
  CHECK(within_stats.online_basis_cap <= within_stats.basis_limit);
  CHECK(within_stats.reduced_solves == num_evaluations);
  CHECK(within_stats.exact_solves == 2);

  // The same basis is above a limit of a single vector.
  auto over = SolveRectangularModes(
      1000.0, 500.0, 500.0, 4.0, 2, num_modes, [](IoData &) {}, true, num_evaluations, 16,
      {0.9, 1.1}, 1.0e-3, 1.0e-12, 1.0e-3);
  REQUIRE(over.num_converged >= num_modes);
  const auto &over_stats = over.reduced_stats;
  CHECK(over_stats.over_basis_limit);
  CHECK(over_stats.basis_limit == 1);
  CHECK(over_stats.reduced_solves == 0);
  CHECK(over_stats.exact_solves == 2 + num_evaluations);
  CHECK(over.reduced_basis_size == over_stats.offline_basis_rank);

  for (const auto *result : {&within, &over})
  {
    REQUIRE(result->kn.size() >= num_modes);
    REQUIRE(result->final_kn.size() >= num_modes);
    for (int i = 0; i < num_modes; i++)
    {
      CHECK_THAT(result->kn[i].real(), WithinRel(exact.kn[i].real(), 1.0e-6));
      CHECK_THAT(result->final_kn[i].real(), WithinRel(exact.kn[i].real(), 1.0e-6));
    }
  }
}

TEST_CASE("ModeEigenSolver reduced basis capacity lifecycle",
          "[boundarymodeoperator][Serial]")
{
  auto result =
      SolveRectangularModes(1000.0, 500.0, 500.0, 4.0, 2, 1, [](IoData &) {}, true, 1, 1);
  REQUIRE(result.num_converged >= 1);
  CHECK(result.reduced_stats.offline_basis_rank == 1);
  CHECK(result.reduced_stats.online_basis_cap == 5);
  CHECK(result.reduced_stats.reduced_solves == 1);
}

TEST_CASE("ModeEigenSolver Impedance shifts kn", "[boundarymodeoperator][Serial]")
{
  auto pec_result = SolveRectangularModes(1000.0, 500.0, 500.0, 4.0, 2, 3, [](IoData &) {});

  // Use a large enough inductance so the impedance shift is well above numerical
  // noise across different BLAS/LAPACK implementations.
  auto imp_result = SolveRectangularModes(1000.0, 500.0, 500.0, 4.0, 2, 3,
                                          [](IoData &iodata)
                                          {
                                            iodata.boundaries.pec.attributes = {1, 3, 4};
                                            auto &imp =
                                                iodata.boundaries.impedance.emplace_back();
                                            imp.attributes = {2};
                                            imp.Ls = 1.0e-8;
                                          });

  REQUIRE(pec_result.num_converged >= 1);
  REQUIRE(imp_result.num_converged >= 1);

  CHECK(imp_result.kn[0].real() > pec_result.kn[0].real());
}

TEST_CASE("ModeEigenSolver rational impedance affine component",
          "[boundarymodeoperator][Serial]")
{
  auto configure_rational = [](IoData &iodata)
  {
    iodata.boundaries.pec.attributes = {1, 3, 4};
    auto &rz = iodata.boundaries.rational_impedance.emplace_back();
    rz.attributes = {2};
    // Series RL impedance Z(s) = R + sL gives the genuinely rational Robin coefficient
    // g(s) = s/(R+sL), rather than merely duplicating the linear resistive component.
    rz.num = {1.0e-12, 50.0};
    rz.den = {1.0};
  };
  auto exact = SolveRectangularModes(1000.0, 500.0, 500.0, 4.0, 2, 1, configure_rational);
  auto reduced =
      SolveRectangularModes(1000.0, 500.0, 500.0, 4.0, 2, 1, configure_rational, true);

  REQUIRE(exact.num_converged >= 1);
  REQUIRE(reduced.num_converged >= 1);
  CHECK(reduced.reduced_stats.reduced_solves == 1);
  CHECK(reduced.reduced_stats.worst_residual <= reduced.reduced_tol);
  CHECK_THAT(reduced.kn[0].real(), WithinRel(exact.kn[0].real(), 1.0e-6));
  CHECK_THAT(reduced.kn[0].imag(), WithinAbs(exact.kn[0].imag(), 1.0e-8));
}

TEST_CASE("ModeEigenSolver reduced rejection fallback enrichment",
          "[boundarymodeoperator][Serial][Parallel]")
{
  auto configure_rational = [](IoData &iodata)
  {
    iodata.boundaries.pec.attributes = {1, 3};
    auto &rz = iodata.boundaries.rational_impedance.emplace_back();
    rz.attributes = {2, 4};
    rz.num = {1.0e-12, 1.0};
    rz.den = {1.0};
  };
  auto result = SolveRectangularModes(1000.0, 500.0, 500.0, 4.0, 2, 1, configure_rational,
                                      true, 2, 16, {0.1}, 1.0e-4, 1.0e-10);

  REQUIRE(result.num_converged >= 1);
  CHECK(result.first_target_stats.offline_basis_rank == 1);
  CHECK(result.first_target_stats.exact_solves == 2);
  CHECK(result.first_target_stats.fallbacks == 1);
  CHECK(result.first_target_stats.reduced_solves == 0);
  CHECK(result.first_target_basis_size > result.first_target_stats.offline_basis_rank);

  CHECK(result.reduced_stats.exact_solves == 2);
  CHECK(result.reduced_stats.fallbacks == 1);
  CHECK(result.reduced_stats.reduced_solves == 1);
  CHECK(result.reduced_basis_size == result.first_target_basis_size);
  CHECK(result.reduced_stats.worst_residual <= result.reduced_tol);
  REQUIRE(result.final_reduced_backward_errors.size() == 1);
  CHECK(result.final_reduced_backward_errors[0] <= result.reduced_tol);
  REQUIRE(result.final_kn.size() == 1);
  CHECK_THAT(result.final_kn[0].real(), WithinRel(result.kn[0].real(), 1.0e-6));
  CHECK_THAT(result.final_kn[0].imag(), WithinAbs(result.kn[0].imag(), 1.0e-8));
}

TEST_CASE("ModeEigenSolver Conductivity adds loss", "[boundarymodeoperator][Serial]")
{
  auto pec_result = SolveRectangularModes(1000.0, 500.0, 500.0, 4.0, 2, 3, [](IoData &) {});

  auto configure_conductivity = [](IoData &iodata)
  {
    iodata.boundaries.pec.attributes = {1, 3, 4};
    auto &cond = iodata.boundaries.conductivity.emplace_back();
    cond.attributes = {2};
    cond.sigma = 5.0e7;
    cond.h = 0.001;
  };
  auto cond_result =
      SolveRectangularModes(1000.0, 500.0, 500.0, 4.0, 2, 3, configure_conductivity);
  auto cond_reduced =
      SolveRectangularModes(1000.0, 500.0, 500.0, 4.0, 2, 1, configure_conductivity, true);

  REQUIRE(pec_result.num_converged >= 1);
  REQUIRE(cond_result.num_converged >= 1);
  REQUIRE(cond_reduced.num_converged >= 1);

  CHECK(std::abs(cond_result.kn[0].imag()) > std::abs(pec_result.kn[0].imag()));
  CHECK(cond_reduced.reduced_stats.reduced_solves == 1);
  CHECK(cond_reduced.reduced_stats.worst_residual <= cond_reduced.reduced_tol);
  CHECK_THAT(cond_reduced.kn[0].real(), WithinRel(cond_result.kn[0].real(), 1.0e-6));
  CHECK_THAT(cond_reduced.kn[0].imag(), WithinAbs(cond_result.kn[0].imag(), 1.0e-8));
}

TEST_CASE("ModeEigenSolver reduced solve with second-order absorbing boundaries",
          "[boundarymodeoperator][Serial][Parallel]")
{
  // The second-order absorbing boundary terms make B depend on frequency, and the reduced
  // model projects its frequency-independent part and the components of the terms
  // separately, as for A. The first reduced solve after training falls back to the exact
  // solve, whose mode enriches the basis, so that the second one is accepted, with the mode
  // of the exact solve.
  auto configure_farfield = [](IoData &iodata)
  {
    iodata.boundaries.pec.attributes = {1, 3};
    iodata.boundaries.farfield.attributes = {2, 4};
    iodata.boundaries.farfield.order = 2;
  };
  auto exact = SolveRectangularModes(1000.0, 500.0, 500.0, 4.0, 2, 1, configure_farfield);
  auto reduced =
      SolveRectangularModes(1000.0, 500.0, 500.0, 4.0, 2, 1, configure_farfield, true, 2);
  REQUIRE(exact.num_converged >= 1);
  REQUIRE(reduced.final_kn.size() == 1);
  CAPTURE(exact.kn[0], reduced.final_kn[0], reduced.reduced_stats.exact_solves,
          reduced.reduced_stats.fallbacks, reduced.reduced_stats.worst_residual,
          reduced.reduced_tol);
  CHECK(exact.kn[0].imag() < -1.0e-6 * exact.kn[0].real());
  CHECK(reduced.reduced_stats.fallbacks == 1);
  CHECK(reduced.reduced_stats.reduced_solves == 1);
  CHECK(reduced.final_reduced_backward_errors[0] <= reduced.reduced_tol);
  CHECK_THAT(reduced.final_kn[0].real(), WithinRel(exact.kn[0].real(), 1.0e-6));
  CHECK_THAT(reduced.final_kn[0].imag(), WithinRel(exact.kn[0].imag(), 1.0e-4));
}

TEST_CASE("ModeEigenSolver p-multigrid preconditioning", "[boundarymodeoperator][Serial]")
{
  // The H1 block of the multigrid preconditioner is negative definite (the diffusion term
  // keeps its sign from the integration by parts), and essential DOFs are eliminated with a
  // unit diagonal, so its Chebyshev smoothers see a diagonal of mixed sign (with essential
  // boundaries) or a negative diagonal (without). The multigrid preconditioned solve must
  // reproduce the mode of the sparse direct one.
  auto check = [](const std::function<void(IoData &)> &configure_bcs)
  {
    const auto direct = SolveRectangularModesMultigrid(1, configure_bcs);
    const auto multigrid = SolveRectangularModesMultigrid(2, configure_bcs);
    REQUIRE(direct.size() >= 1);
    REQUIRE(multigrid.size() >= 1);
    CAPTURE(direct[0], multigrid[0]);
    CHECK_THAT(multigrid[0].real(), WithinRel(direct[0].real(), 1.0e-6));
    CHECK_THAT(multigrid[0].imag(),
               WithinAbs(direct[0].imag(), 1.0e-6 * std::abs(direct[0])));
  };

  SECTION("PEC walls")
  {
    check([](IoData &) {});
  }

  SECTION("Impedance walls")
  {
    check(
        [](IoData &iodata)
        {
          iodata.boundaries.pec.attributes.clear();
          auto &imp = iodata.boundaries.impedance.emplace_back();
          imp.attributes = {1, 2, 3, 4};
          imp.Ls = 1.0e-8;
        });
  }
}

TEST_CASE("ModeEigenSolver preconditioners of lossy modes",
          "[boundarymodeoperator][Serial][Parallel]")
{
  // The linear solves of the mode eigenproblem with each preconditioner: sparse direct (the
  // exact complex matrix, or its real-valued approximation) or p-multigrid (complex or
  // real-valued), for loss terms which dominate the system matrix where they apply. A
  // preconditioner without them stagnates the Krylov solver, or with left preconditioning,
  // can meet the tolerance with an inaccurate solution.
  struct Preconditioner
  {
    int mg_max_levels;
    bool pc_mat_real, complex_coarse_solve;
  };
  const Preconditioner preconditioners[] = {{1, false, true},  {1, false, false},
                                            {1, true, false},  {2, false, true},
                                            {2, false, false}, {2, true, false}};
  // The shift target is that of the BoundaryMode driver, above the index of the dielectric:
  // a target equal to it would cancel the mass term of the shifted system matrix for the
  // in-plane field in the dielectric, which makes it singular or nearly so.
  constexpr double t = 2.0, d = 3.0, freq_ghz = 10.0, eps_d = 3.9;
  const double n_target = std::sqrt(1.1 * eps_d);
  auto check = [&](const std::function<void(IoData &)> &configure,
                   std::complex<double> n_ref, double atol)
  {
    for (const auto &pc : preconditioners)
    {
      const auto result =
          SolveParallelPlateMode(t, d, freq_ghz, n_target, pc.mg_max_levels, pc.pc_mat_real,
                                 pc.complex_coarse_solve, configure);
      CAPTURE(pc.mg_max_levels, pc.pc_mat_real, pc.complex_coarse_solve, result.n_eff,
              n_ref, result.num_solves, result.num_iterations);
      REQUIRE(result.num_solves > 0);
      // All solves converge, in at most a fourth of the maximum number of iterations on
      // average (here at most 34 for multigrid and 6 for the sparse direct preconditioner).
      CHECK(result.num_iterations <= 50 * result.num_solves);
      CHECK(std::abs(result.n_eff - n_ref) <= atol);
    }
  };

  SECTION("Volume conductor")
  {
    // Aluminum (sigma = 3.77e7 S/m, skin depth 0.82 um at 10 GHz) 2 um thick on a PEC plane
    // and oxide 3 um thick up to a second PEC plane: the current along the line flows in
    // the volume of the aluminum, where the conductivity term of the system matrix exceeds
    // the permittivity term by sigma / (omega eps0) = 7e7, compared to the exact solution
    // (the discretization error is 3e-5 of the imaginary part of the effective index).
    constexpr double sigma = 3.77e7;
    const auto n_exact = ExactParallelPlateConductorIndex(t, d, sigma, eps_d, freq_ghz);
    check(
        [&](IoData &iodata)
        {
          auto &conductor = iodata.domains.materials.emplace_back();
          conductor.attributes = {1};
          conductor.sigma.s = {sigma, sigma, sigma};
          auto &oxide = iodata.domains.materials.emplace_back();
          oxide.attributes = {2};
          oxide.epsilon_r.s = {eps_d, eps_d, eps_d};
          iodata.boundaries.pec.attributes = {1, 3};
          iodata.boundaries.pmc.attributes = {2, 4};
        },
        n_exact, 2.0e-4 * std::abs(n_exact.imag()));
  }

  SECTION("Resistive sheet")
  {
    // Oxide between a PEC plane and a resistive sheet of 1 mOhm per square, whose damping
    // term dominates the boundary terms of the system matrix, compared to the solution
    // with the exact complex preconditioner.
    auto configure = [&](IoData &iodata)
    {
      auto &oxide = iodata.domains.materials.emplace_back();
      oxide.attributes = {1, 2};
      oxide.epsilon_r.s = {eps_d, eps_d, eps_d};
      iodata.boundaries.pec.attributes = {1};
      iodata.boundaries.pmc.attributes = {2, 4};
      auto &sheet = iodata.boundaries.impedance.emplace_back();
      sheet.attributes = {3};
      sheet.Rs = 1.0e-3;
    };
    const auto n_ref =
        SolveParallelPlateMode(t, d, freq_ghz, n_target, 1, false, true, configure).n_eff;
    CHECK(n_ref.imag() < -1.0e-4 * n_ref.real());
    check(configure, n_ref, 1.0e-3 * std::abs(n_ref.imag()));
  }
}

TEST_CASE("ModeEigenSolver multigrid preconditioners of impedance sheets",
          "[boundarymodeoperator][Serial][Parallel]")
{
  // The p-multigrid preconditioners of the in-plane (ND) and out-of-plane (H1) blocks have
  // the terms of the impedance boundaries of the system matrix. Here oxide between two
  // sheets of 1 mOhm per square in parallel with 10 fH per square (impedance boundaries, or
  // rational impedance boundaries of the same parallel RL admittance), whose terms dominate
  // the boundary terms of the system matrix. A preconditioner without them stagnates the
  // Krylov solver (in-plane block) or doubles its number of iterations (out-of-plane block;
  // here about 17 per solve with the complex and 24 with the real-valued preconditioner).
  // The same holds for finite conductivity boundaries of the same surface resistance. The
  // solutions are compared to the solution with the exact complex sparse direct
  // preconditioner.
  constexpr double t = 2.0, d = 3.0, freq_ghz = 10.0, eps_d = 3.9, Rs = 1.0e-3,
                   Ls = 1.0e-14;
  const double n_target = std::sqrt(1.1 * eps_d);
  auto check = [&](const std::function<void(IoData &)> &configure_sheets)
  {
    auto configure = [&](IoData &iodata)
    {
      auto &oxide = iodata.domains.materials.emplace_back();
      oxide.attributes = {1, 2};
      oxide.epsilon_r.s = {eps_d, eps_d, eps_d};
      iodata.boundaries.pmc.attributes = {2, 4};
      configure_sheets(iodata);
    };
    const auto ref =
        SolveParallelPlateMode(t, d, freq_ghz, n_target, 1, false, true, configure);
    for (bool pc_mat_real : {false, true})
    {
      const auto result = SolveParallelPlateMode(t, d, freq_ghz, n_target, 2, pc_mat_real,
                                                 !pc_mat_real, configure);
      CAPTURE(pc_mat_real, result.n_eff, ref.n_eff, result.num_solves,
              result.num_iterations);
      REQUIRE(result.num_solves > 0);
      CHECK(result.num_iterations <= 35 * result.num_solves);
      CHECK(std::abs(result.n_eff - ref.n_eff) <= 1.0e-6 * std::abs(ref.n_eff));
    }
  };

  SECTION("Impedance")
  {
    check(
        [&](IoData &iodata)
        {
          auto &sheets = iodata.boundaries.impedance.emplace_back();
          sheets.attributes = {1, 3};
          sheets.Rs = Rs;
          sheets.Ls = Ls;
        });
  }

  SECTION("Rational impedance")
  {
    // Zs(s) = s Ls Rs / (s Ls + Rs).
    check(
        [&](IoData &iodata)
        {
          auto &sheets = iodata.boundaries.rational_impedance.emplace_back();
          sheets.attributes = {1, 3};
          sheets.num = {Ls * Rs, 0.0};
          sheets.den = {Ls, Rs};
        });
  }

  SECTION("Surface conductivity")
  {
    // Finite conductivity boundaries of the same surface resistance, 1 mOhm per square at
    // 10 GHz (sqrt(mu0 omega / (2 sigma)) for sigma = 4e10 S/m).
    check(
        [&](IoData &iodata)
        {
          auto &sheets = iodata.boundaries.conductivity.emplace_back();
          sheets.attributes = {1, 3};
          sheets.sigma = 4.0e10;
        });
  }
}

TEST_CASE("ModeOperatorModel farfield damping uses the neighboring material",
          "[boundarymodeoperator][Serial][Parallel]")
{
  // The first-order absorbing boundary condition adds iω/Z₀ times the boundary mass of the
  // in-plane (ND, tangential trace) and out-of-plane (H1, with the negative sign of the
  // H1 block) field components, where 1/Z₀ = √(ε/μ) is taken from the domain material
  // adjacent to each boundary element. The absorbing side x = W spans both materials, with
  // 1/Z₀ = 2 along y < H/4 and 1/Z₀ = 1 above, so a unit field tangential to it gives
  // ±(2 H/4 + 3 H/4) (the frequency scalar iω is applied separately).
  using RMM = RectangleModeModel;
  RMM rect([](IoData &iodata) { iodata.boundaries.farfield.attributes = {2}; });
  const auto &component = rect.GetComponent(mode_assembly::CoefficientType::OMEGA);
  const double expected = 2.0 * 0.25 * RMM::H + 1.0 * 0.75 * RMM::H;
  const auto pair_t = rect.Pairing(component, 1.0, 0.0);
  const auto pair_n = rect.Pairing(component, 0.0, 1.0);
  CHECK(pair_t.real() == 0.0);
  CHECK_THAT(pair_t.imag(), WithinRel(expected, 1.0e-12));
  CHECK(pair_n.real() == 0.0);
  CHECK_THAT(pair_n.imag(), WithinRel(-expected, 1.0e-12));
}

TEST_CASE("ModeEigenSolver absorbing boundaries of first and second order",
          "[boundarymodeoperator][Serial][Parallel]")
{
  // For m = 0 and m = 1, the modes Eₓ = 0 of test::ExactAbsorbingParallelPlateKn.
  // For m = 0, Eₙ = 0, and for m = 1, Eₙ varies along the absorbing boundaries, so
  // that the boundary terms of the second-order condition couple it to the in-plane
  // field. For m = 1, the mode Hₓ = 0 of the second-order condition is that of the
  // first-order condition. The two conditions give modes 1e-2 apart, and the
  // discretization error of these modes is less than 1e-7.
  const int abc_order = GENERATE(1, 2);
  const double omega = 2.0 * std::numbers::pi *
                       Units(1.0, 1.0).Nondimensionalize<Units::ValueType::FREQUENCY>(
                           absorbing_parallel_plate_freq_ghz);
  for (int m : {0, 1})
  {
    // The target is the real part of the exact mode.
    const auto kn_exact =
        test::ExactAbsorbingParallelPlateKn(1.0, 0.5, omega, m, abc_order);
    const auto kn_other =
        test::ExactAbsorbingParallelPlateKn(1.0, 0.5, omega, m, 3 - abc_order);
    const auto modes = SolveAbsorbingParallelPlateModes(abc_order, kn_exact.real() / omega,
                                                        3, 16, 1, false);
    CAPTURE(abc_order, m, omega, modes.kn, kn_exact, kn_other);
    CHECK(HasMode(modes.kn, kn_exact, 1.0e-6));
    CHECK(std::abs(kn_other - kn_exact) > 1.0e-2 * std::abs(kn_exact));
    if (m > 0 && abc_order > 1)
    {
      CHECK(HasMode(modes.kn, kn_other, 1.0e-6));
    }
  }
}

TEST_CASE("ModeEigenSolver multigrid preconditioners of second-order absorbing boundaries",
          "[boundarymodeoperator][Serial][Parallel]")
{
  // The p-multigrid preconditioners of the in-plane (ND) and out-of-plane (H1) blocks, and
  // their coupling, have the terms of the second-order absorbing boundaries of the system
  // matrix (or their magnitudes, for the real-valued preconditioner), for the modes of the
  // previous test with m = 1 on a coarser mesh, compared to the solution with the sparse
  // direct preconditioner (here 70 to 80 iterations per solve, and 60 to 70 with
  // first-order absorbing boundaries).
  const double omega = 2.0 * std::numbers::pi *
                       Units(1.0, 1.0).Nondimensionalize<Units::ValueType::FREQUENCY>(
                           absorbing_parallel_plate_freq_ghz);
  const double n_target =
      test::ExactAbsorbingParallelPlateKn(1.0, 0.5, omega, 1, 2).real() / omega;
  const auto ref = SolveAbsorbingParallelPlateModes(2, n_target, 2, 8, 1, false);
  for (bool pc_mat_real : {false, true})
  {
    const auto modes = SolveAbsorbingParallelPlateModes(2, n_target, 2, 8, 2, pc_mat_real);
    CAPTURE(pc_mat_real, modes.kn, ref.kn, modes.num_solves, modes.num_iterations);
    REQUIRE(modes.num_solves > 0);
    CHECK(modes.num_iterations <= 100 * modes.num_solves);
    for (const auto &kn : ref.kn)
    {
      CHECK(HasMode(modes.kn, kn, 1.0e-8));
    }
  }
}

TEST_CASE("ModeOperatorModel impedance damping with an anisotropic absorbing boundary",
          "[boundarymodeoperator][Serial][Parallel]")
{
  // The absorbing side y = 0 lies along the substrate, here of permittivity diag(4, 9, 1),
  // so its in-plane inverse impedance diag(2, 3) has the first entry of the conductance 2
  // of the impedance side x = W. The impedance side must keep its own coefficient 2 I: a
  // unit field ŷ, tangential to it and normal to the absorbing side, gives 2 H in the
  // in-plane block, and 3 H with the tensor of the absorbing side.
  using RMM = RectangleModeModel;
  RMM rect(
      [](IoData &iodata)
      {
        iodata.domains.materials[0].epsilon_r.s = {4.0, 9.0, 1.0};
        iodata.boundaries.farfield.attributes = {1};
        auto &imp = iodata.boundaries.impedance.emplace_back();
        imp.attributes = {2};
        imp.Rs = 0.5 * electromagnetics::Z0_;
      });
  const auto &component = rect.GetComponent(mode_assembly::CoefficientType::OMEGA);
  const auto pair_t = rect.Pairing(component, 1.0, 0.0);
  CHECK(pair_t.real() == 0.0);
  CHECK_THAT(pair_t.imag(), WithinRel(2.0 * RMM::H, 1.0e-12));
}

TEST_CASE("ModeOperatorModel surface capacitance of impedance boundaries",
          "[boundarymodeoperator][Serial][Parallel]")
{
  // The surface capacitance Cs of an impedance boundary adds -ω² Cs times the boundary mass
  // of the in-plane (ND, tangential trace) field component and +ω² Cs times that of the
  // out-of-plane (H1) component, with the sign of the H1 block, so for a unit field
  // tangential to the impedance side x = W, the quadratic component differs by ∓Cs H from
  // that of the cross-section without it (the frequency scalar ω² is applied separately).
  using RMM = RectangleModeModel;
  const double Cs = 2.0e-6 * electromagnetics::epsilon0_;  // F/sq
  RMM bare([](IoData &) {});
  RMM capacitive(
      [&](IoData &iodata)
      {
        auto &imp = iodata.boundaries.impedance.emplace_back();
        imp.attributes = {2};
        imp.Cs = Cs;
      });
  const double expected =
      capacitive.GetUnits().Nondimensionalize<Units::ValueType::CAPACITANCE>(Cs) * RMM::H;
  const auto &bare_component =
      bare.GetComponent(mode_assembly::CoefficientType::OMEGA_SQUARED);
  const auto &component =
      capacitive.GetComponent(mode_assembly::CoefficientType::OMEGA_SQUARED);
  const auto pair_t =
      capacitive.Pairing(component, 1.0, 0.0) - bare.Pairing(bare_component, 1.0, 0.0);
  const auto pair_n =
      capacitive.Pairing(component, 0.0, 1.0) - bare.Pairing(bare_component, 0.0, 1.0);
  CHECK_THAT(pair_t.real(), WithinRel(-expected, 1.0e-10));
  CHECK(pair_t.imag() == 0.0);
  CHECK_THAT(pair_n.real(), WithinRel(expected, 1.0e-10));
  CHECK(pair_n.imag() == 0.0);
}

TEST_CASE("ModeOperatorModel dispersive permittivity matches its value at the frequency",
          "[boundarymodeoperator][Serial][Parallel]")
{
  // At any one frequency, a dispersive substrate is the constant material with ε = ε∞ +
  // χ(iω), the imaginary part entering as a loss tangent, so the two assemble the same mode
  // operator, at real and at complex ω. Rates are entered in rad/s and nondimensionalized
  // by t_c = L₀/c₀: Debye, complex pole pair, and Djordjevic-Sarkar terms.
  auto make_dispersive = [](IoData &iodata)
  {
    const double rate = electromagnetics::c0_ / iodata.model.L0;
    iodata.problem.type = ProblemType::DRIVEN;
    auto &substrate = iodata.domains.materials[0];
    substrate.epsilon_r.s = {4.0, 4.0, 4.0};
    substrate.permittivity_terms.push_back(
        config::PermittivityPoleData{{-1.0 * rate, 0.0}, {1.5 * rate, 0.0}});
    substrate.permittivity_terms.push_back(
        config::PermittivityPoleData{{-0.1 * rate, 2.5 * rate}, {0.3 * rate, -0.2 * rate}});
    substrate.permittivity_terms.push_back(
        config::DjordjevicSarkarData{0.2, 0.01 * rate, 100.0 * rate});
  };
  const RectangleModeModel dispersive(make_dispersive);
  REQUIRE(dispersive.GetMaterialOp().HasFrequencyDependentPermittivityA2(0));

  for (const std::complex<double> omega :
       {std::complex<double>(2.0, 0.0), std::complex<double>(2.0, 0.1)})
  {
    CAPTURE(omega);
    const std::complex<double> s = std::complex<double>(0.0, 1.0) * omega;
    const std::complex<double> chi =
        dispersive.GetMaterialOp().EvaluateFrequencyDependentPermittivityA2(0, s) / (s * s);
    const double eps = 4.0 + chi.real(), tandelta = -chi.imag() / eps;
    const RectangleModeModel constant(
        [eps, tandelta](IoData &iodata)
        {
          iodata.problem.type = ProblemType::DRIVEN;
          auto &substrate = iodata.domains.materials[0];
          substrate.epsilon_r.s = {eps, eps, eps};
          substrate.tandelta.s = {tandelta, tandelta, tandelta};
        });

    const double sigma = -2.0e-1;
    const auto A = dispersive.GetModel().Assemble(omega, sigma);
    const auto B = constant.GetModel().Assemble(omega, sigma);
    ComplexVector x(A->Width()), y(A->Height()), z(A->Height());
    x.UseDevice(true);
    y.UseDevice(true);
    z.UseDevice(true);
    linalg::SetRandom(Mpi::World(), x);
    A->Mult(x, y);
    B->Mult(x, z);
    const double norm = linalg::Norml2(Mpi::World(), z);
    z -= y;
    CHECK(linalg::Norml2(Mpi::World(), z) <= 1.0e-12 * norm);
  }
}

TEST_CASE("ModeOperatorModel dispersive permittivity at absorbing port edges",
          "[boundarymodeoperator][Serial][Parallel]")
{
  // Absorbing edges at x = 0 and x = W neighbor the substrate and vacuum. A lossless
  // Lorentz substrate at real ω is the constant material ε(ω), whose absorbing condition
  // uses √(ε(ω)/μ) (and, for the second-order condition, the light speed of ε(ω)), so the
  // two assemble the same mode operator A, and the same B.
  const int order = GENERATE(1, 2);
  const double eps_inf = 4.0, w0 = 3.0, delta = 0.7;
  const double omega = 2.0;
  const double eps = eps_inf + delta * w0 * w0 / (w0 * w0 - omega * omega);
  auto make = [&](bool dispersive)
  {
    return [=](IoData &iodata)
    {
      const double rate = electromagnetics::c0_ / iodata.model.L0;  // Rates in rad/s
      iodata.problem.type = ProblemType::DRIVEN;
      iodata.boundaries.farfield.attributes = {2, 4};
      iodata.boundaries.farfield.order = order;
      auto &substrate = iodata.domains.materials[0];
      if (dispersive)
      {
        substrate.epsilon_r.s = {eps_inf, eps_inf, eps_inf};
        substrate.permittivity_terms.push_back(config::LorentzData{delta, w0 * rate, 0.0});
      }
      else
      {
        substrate.epsilon_r.s = {eps, eps, eps};
      }
    };
  };
  const RectangleModeModel dispersive(make(true)), constant(make(false));
  REQUIRE_THAT(dispersive.GetMaterialOp()
                   .EvaluateFrequencyDependentPermittivity(0, {0.0, omega})
                   .real(),
               WithinRel(eps, 1.0e-12));
  CHECK(dispersive.GetModel().HasFrequencyDependentB() == (order > 1));
  const double sigma = -2.0e-1;
  const auto A = dispersive.GetModel().Assemble(omega, sigma);
  const auto B = constant.GetModel().Assemble(omega, sigma);
  ComplexVector x(A->Width()), y(A->Height()), z(A->Height());
  x.UseDevice(true);
  y.UseDevice(true);
  z.UseDevice(true);
  linalg::SetRandom(Mpi::World(), x);
  A->Mult(x, y);
  B->Mult(x, z);
  const double norm = linalg::Norml2(Mpi::World(), z);
  z -= y;
  CHECK(linalg::Norml2(Mpi::World(), z) <= 1.0e-12 * norm);

  // The frequency-dependent part of B, the sum of the components scaled by their
  // coefficients in B.
  auto ApplyB = [&](const RectangleModeModel &rect, ComplexVector &bx)
  {
    ComplexVector t(bx.Size());
    t.UseDevice(true);
    bx = 0.0;
    for (const auto &component : rect.GetModel().GetComponents())
    {
      const auto coeff = rect.GetModel().EvaluateBCoefficient(component, omega);
      if (coeff != 0.0)
      {
        component.op->Mult(x, t);
        bx.Add(coeff, t);
      }
    }
  };
  ApplyB(dispersive, y);
  ApplyB(constant, z);
  const double norm_b = linalg::Norml2(Mpi::World(), z);
  CHECK((norm_b > 0.0) == (order > 1));
  z -= y;
  CHECK(linalg::Norml2(Mpi::World(), z) <= 1.0e-12 * norm_b);
}

TEST_CASE("ModeOperatorModel decouples essential DOFs in every component",
          "[boundarymodeoperator][Serial][Parallel]")
{
  // Each component, and so the assembled matrix, must decouple the essential DOFs from the
  // rest of the system, with a unit diagonal from the constant component. This includes
  // the ND/H1 coupling blocks, which appear in a component without their transposed
  // counterpart (Atn in the constant and -Btn in the shift component, and those of the
  // second-order absorbing boundary, on the side x = W, which ends on a PEC side).
  RectangleModeModel rect(
      [](IoData &iodata)
      {
        iodata.boundaries.pec.attributes = {1, 4};
        iodata.boundaries.farfield.attributes = {2};
        iodata.boundaries.farfield.order = 2;
      });
  MPI_Comm comm = Mpi::World();
  const auto &dbc_tdof_list = rect.GetDbcTDofList();
  int num_dbc = dbc_tdof_list.Size();
  Mpi::GlobalSum(1, &num_dbc, comm);
  REQUIRE(num_dbc > 0);

  // Apply A to x supported on the essential or on the other DOFs, and return the norms of
  // y = A x - c x on the essential and on the other DOFs (with c the expected diagonal of
  // the essential DOFs for x on these).
  auto Residuals = [&](const ComplexOperator &A, bool essential, double c)
  {
    const int n = A.Height();
    std::vector<bool> is_dbc(n, false);
    for (auto tdof : dbc_tdof_list)
    {
      is_dbc[tdof] = true;
    }
    ComplexVector x(n), y(n);
    x = 0.0;
    for (int i = 0; i < n; i++)
    {
      if (is_dbc[i] == essential)
      {
        x.Real()[i] = 1.0 + 0.1 * (i % 7);
        x.Imag()[i] = 0.5 - 0.1 * (i % 5);
      }
    }
    A.Mult(x, y);
    if (essential)
    {
      y.AXPY(-c, x);
    }
    double norms[2] = {0.0, 0.0};  // Essential, other DOFs
    for (int i = 0; i < n; i++)
    {
      norms[is_dbc[i] ? 0 : 1] += std::norm(std::complex<double>(y.Real()[i], y.Imag()[i]));
    }
    Mpi::GlobalSum(2, norms, comm);
    return std::array<double, 2>{std::sqrt(norms[0]), std::sqrt(norms[1])};
  };

  for (const auto &component : rect.GetModel().GetComponents())
  {
    CAPTURE(static_cast<int>(component.type));
    const double c =
        (component.type == mode_assembly::CoefficientType::CONSTANT) ? 1.0 : 0.0;
    const auto from_dbc = Residuals(*component.op, true, c);
    CHECK(from_dbc[0] == 0.0);
    CHECK(from_dbc[1] == 0.0);
    const auto from_other = Residuals(*component.op, false, 0.0);
    CHECK(from_other[0] == 0.0);
  }
  const auto A = rect.GetModel().Assemble(1.0e-2, -2.0e-4);
  const auto from_dbc = Residuals(*A, true, 1.0);
  CHECK(from_dbc[0] == 0.0);
  CHECK(from_dbc[1] == 0.0);
  CHECK(Residuals(*A, false, 0.0)[0] == 0.0);
}

TEST_CASE("ModeOperatorModel boundary components on processes without the boundary",
          "[boundarymodeoperator][Serial][Parallel]")
{
  // Surface conductivity and rational impedance boundaries contribute a unit boundary mass
  // component each (the frequency dependence is applied separately), which is assembled
  // collectively also by processes without elements on the boundary (all but the last
  // here). The conductivity term comes first, scaled by its surface admittance iω/Z(ω),
  // then the rational impedance term, scaled by its Robin coefficient iω/Zs(iω).
  using RMM = RectangleModeModel;
  RMM rect(
      [](IoData &iodata)
      {
        auto &cond = iodata.boundaries.conductivity.emplace_back();
        cond.attributes = {2};
        cond.sigma = 5.0e7;
        auto &rz = iodata.boundaries.rational_impedance.emplace_back();
        rz.attributes = {2};
        rz.num = {1.0e-12, 50.0};
        rz.den = {1.0};
      });
  const std::complex<double> omega = 1.0e-2;
  const auto &terms = rect.GetSurfaceTerms();
  const std::complex<double> scales[] = {
      terms.GetSurfaceConductivityOp().EvaluateScalar(0, omega),
      terms.GetRationalImpedanceOp().EvalRobinCoefficient(
          0, std::complex<double>(0.0, 1.0) * omega)};
  REQUIRE(std::abs(scales[0] - scales[1]) > 1.0e-3 * std::abs(scales[0]));
  std::vector<int> indices;
  for (const auto &component : rect.GetModel().GetComponents())
  {
    if (component.type != mode_assembly::CoefficientType::SURFACE_TERM)
    {
      continue;
    }
    CAPTURE(component.index);
    indices.push_back(component.index);
    REQUIRE((component.index == 0 || component.index == 1));
    CHECK(rect.GetModel().EvaluateCoefficient(component, omega, 0.0) ==
          scales[component.index]);
    const auto pair_t = rect.Pairing(component, 1.0, 0.0);
    const auto pair_n = rect.Pairing(component, 0.0, 1.0);
    CHECK_THAT(pair_t.real(), WithinRel(RMM::H, 1.0e-12));
    CHECK(pair_t.imag() == 0.0);
    CHECK_THAT(pair_n.real(), WithinRel(-RMM::H, 1.0e-12));
    CHECK(pair_n.imag() == 0.0);
  }
  CHECK(indices == std::vector<int>{0, 1});
}

TEST_CASE("ModeOperatorModel second-order absorbing boundary components",
          "[boundarymodeoperator][Serial][Parallel]")
{
  // The second-order absorbing boundary on the side x = W, which only the last process has
  // elements on, adds the components -[0, G; 0, D] and -[M, 0; Gᵀ, 0], with coefficients β
  // and σβ in A, and -β in B for the latter, where β = 1/(2iω). With the coefficient
  // μ⁻¹c₀ = 1/2 along the substrate (y < H/4) and 1 above, a unit field ŷ, tangential to
  // the side, gives -(H/8 + 3H/4) for M, and a constant out-of-plane field gives zero (no
  // tangential derivative along the side).
  using RMM = RectangleModeModel;
  RMM rect(
      [](IoData &iodata)
      {
        iodata.boundaries.farfield.attributes = {2};
        iodata.boundaries.farfield.order = 2;
      });
  const std::complex<double> omega = 1.0e-2,
                             beta = 0.5 / (std::complex<double>(0.0, 1.0) * omega);
  constexpr double sigma = -2.0e-4;
  CHECK(rect.GetModel().HasFrequencyDependentB());
  int num_components = 0;
  for (const auto &component : rect.GetModel().GetComponents())
  {
    const bool shift =
        (component.type == mode_assembly::CoefficientType::SECOND_ORDER_FARFIELD_SHIFT);
    if (!shift && component.type != mode_assembly::CoefficientType::SECOND_ORDER_FARFIELD)
    {
      CHECK(rect.GetModel().EvaluateBCoefficient(component, omega) == 0.0);
      continue;
    }
    CAPTURE(shift);
    num_components++;
    CHECK(component.index == 0);
    CHECK(rect.GetModel().EvaluateCoefficient(component, omega, sigma) ==
          (shift ? sigma * beta : beta));
    CHECK(rect.GetModel().EvaluateBCoefficient(component, omega) == (shift ? -beta : 0.0));
    const auto pair_t = rect.Pairing(component, 1.0, 0.0);
    const auto pair_n = rect.Pairing(component, 0.0, 1.0);
    if (shift)
    {
      CHECK_THAT(pair_t.real(), WithinRel(-(0.5 * 0.25 + 0.75) * RMM::H, 1.0e-12));
    }
    else
    {
      CHECK_THAT(pair_t.real(), WithinAbs(0.0, 1.0e-12));
    }
    CHECK(pair_t.imag() == 0.0);
    CHECK_THAT(pair_n.real(), WithinAbs(0.0, 1.0e-12));
    CHECK(pair_n.imag() == 0.0);
  }
  CHECK(num_components == 2);
}

TEST_CASE("ModeEigenSolver London slab", "[boundarymodeoperator][Serial][Parallel]")
{
  // London superconductor slabs of thickness d, compared to exact solutions.
  constexpr double lambda = 0.1;
  SECTION("Out-of-plane current")
  {
    // TM0 mode of a parallel plate guide with gap h, with the slab carrying the current
    // along the propagation direction, compared to the exact dispersion relation. The
    // frequency is high enough for the propagation constant to be well conditioned.
    constexpr double h = 0.5, freq_ghz = 500.0;
    const double k0 =
        2.0 * std::numbers::pi * freq_ghz * 1.0e9 / electromagnetics::c0_ * 1.0e-6;
    for (double d : {0.01, 0.1, 1.0})
    {
      for (bool pec_back : {false, true})
      {
        const double n_exact = ExactLondonSlabIndex(h, d, lambda, k0, pec_back);
        const double n_eff =
            SolveLondonSlabMode(1.0, h, d, lambda, pec_back, false, freq_ghz, n_exact, 3);
        CAPTURE(d, pec_back, n_exact, n_eff);
        CHECK_THAT(n_eff, WithinRel(n_exact, 1.0e-8));
      }
    }
  }
  SECTION("In-plane current")
  {
    // TE1 mode between the PEC plate and the slab (PEC sides), with E_x = sin(k_y y) in
    // the gap and the current along the slab, perpendicular to the propagation direction.
    // In the slab, q^2 = 1/lambda^2 - k_y^2 and k_y cot(k_y h) = -q tanh(q d) (free back
    // face).
    constexpr double h = 500.0, d = 0.5, freq_ghz = 500.0;
    const double k0 =
        2.0 * std::numbers::pi * freq_ghz * 1.0e9 / electromagnetics::c0_ * 1.0e-6;
    auto F = [&](double ky)
    {
      const double q = std::sqrt(1.0 / (lambda * lambda) - ky * ky);
      return ky / std::tan(ky * h) + q * std::tanh(q * d);
    };
    double lo = 0.9 * std::numbers::pi / h, hi = (1.0 - 1.0e-12) * std::numbers::pi / h;
    for (int it = 0; it < 200; it++)
    {
      const double mid = 0.5 * (lo + hi);
      ((F(lo) * F(mid) <= 0.0) ? hi : lo) = mid;
    }
    const double ky = 0.5 * (lo + hi);
    const double n_exact = std::sqrt(1.0 - (ky / k0) * (ky / k0));
    const double n_pec =
        std::sqrt(1.0 - (std::numbers::pi / h / k0) * (std::numbers::pi / h / k0));
    const double n_eff =
        SolveLondonSlabMode(50.0, h, d, lambda, false, true, freq_ghz, n_exact, 3, 32);
    CAPTURE(n_exact, n_pec, n_eff);
    CHECK_THAT(n_eff - n_pec, WithinRel(n_exact - n_pec, 1.0e-5));
  }
}

TEST_CASE("BoundaryModeOperator interior PEC on a nonconforming mesh",
          "[boundarymodeoperator][Serial][Parallel]")
{
  // A PEC line through the interior of the cross-section at x = 4, and an interior
  // boundary without boundary condition at y = 2 for x > 4, with the elements of the
  // quadrant x > 4, y > 2 refined. The edges of both lines on the unrefined side are
  // master edges without boundary elements: those of the PEC line must be essential, and
  // those of the other line not.
  const auto type = GENERATE(mfem::Element::TRIANGLE, mfem::Element::QUADRILATERAL);
  const int order = GENERATE(1, 2);
  MPI_Comm comm = Mpi::World();
  constexpr double width = 8.0, height = 4.0, x_pec = 4.0, y_other = 2.0, eps = 1.0e-9;
  constexpr int pec_attr = 5, other_attr = 6;

  auto serial_mesh = std::make_unique<mfem::Mesh>(
      mfem::Mesh::MakeCartesian2D(8, 4, type, false, width, height));
  {
    auto &smesh = *serial_mesh;
    mfem::Array<int> v;
    for (int f = 0; f < smesh.GetNumFaces(); f++)
    {
      int e1, e2;
      smesh.GetFaceElements(f, &e1, &e2);
      if (e1 < 0 || e2 < 0)
      {
        continue;
      }
      smesh.GetFaceVertices(f, v);
      const double *a = smesh.GetVertex(v[0]), *b = smesh.GetVertex(v[1]);
      if (std::abs(a[0] - x_pec) < eps && std::abs(b[0] - x_pec) < eps)
      {
        smesh.AddBdrSegment(v[0], v[1], pec_attr);
      }
      else if (std::abs(a[1] - y_other) < eps && std::abs(b[1] - y_other) < eps &&
               std::min(a[0], b[0]) > x_pec - eps)
      {
        smesh.AddBdrSegment(v[0], v[1], other_attr);
      }
    }
    smesh.FinalizeTopology();
    smesh.Finalize(true, true);
    smesh.EnsureNCMesh(true);
    mfem::Array<int> marked;
    mfem::Vector c;
    for (int e = 0; e < smesh.GetNE(); e++)
    {
      smesh.GetElementCenter(e, c);
      if (c(0) > x_pec && c(1) > y_other)
      {
        marked.Append(e);
      }
    }
    smesh.GeneralRefinement(marked, 1, 0);
  }

  Units units(1.0, 1.0);
  IoData iodata(units);
  iodata.model.Lc = 1.0;
  auto &material = iodata.domains.materials.emplace_back();
  material.attributes = {1};
  material.epsilon_r.s = {1.0, 1.0, 1.0};
  iodata.boundaries.pec.attributes = {pec_attr};
  iodata.solver.order = order;
  iodata.solver.boundary_mode.freq = 1.0;
  iodata.solver.boundary_mode.n = 1;
  iodata.NondimensionalizeInputs(serial_mesh);
  auto par_mesh = std::make_unique<mfem::ParMesh>(comm, *serial_mesh);
  iodata.CheckConfiguration();
  std::vector<std::unique_ptr<Mesh>> mesh;
  mesh.push_back(std::make_unique<Mesh>(std::move(par_mesh)));
  MaterialOperator mat_op(iodata, *mesh.back());
  BoundaryModeOperator mode_op(iodata, mesh, mat_op);

  auto CheckDbcTDofs = [&](const FiniteElementSpace &fespace,
                           const mfem::Array<int> &dbc_tdof_list, bool edge_dofs)
  {
    const auto &fes = fespace.Get();
    const auto &pmesh = *fes.GetParMesh();
    std::vector<bool> is_dbc(fes.GetTrueVSize(), false);
    for (auto t : dbc_tdof_list)
    {
      is_dbc[t] = true;
    }
    // Numbers of checked true DOFs on the PEC line, on the other line, and on the master
    // edges of the PEC line (which have DOFs if the space has edge DOFs).
    int counts[3] = {0, 0, 0};
    auto Check = [&](const mfem::Array<int> &dofs, bool pec, bool master)
    {
      for (auto d : dofs)
      {
        const int t = fes.GetLocalTDofNumber((d >= 0) ? d : -1 - d);
        if (t < 0)
        {
          continue;
        }
        CHECK(is_dbc[t] == pec);
        counts[pec ? 0 : 1]++;
        counts[2] += (pec && master);
      }
    };
    auto OnPEC = [&](const double *x) { return std::abs(x[0] - x_pec) < eps; };
    auto OnOther = [&](const double *x)
    { return std::abs(x[1] - y_other) < eps && x[0] > x_pec + eps && x[0] < width - eps; };
    const auto &edge_list = pmesh.ncmesh->GetEdgeList();
    mfem::Array<int> v, dofs;
    for (int e = 0; e < pmesh.GetNEdges(); e++)
    {
      pmesh.GetEdgeVertices(e, v);
      const double *a = pmesh.GetVertex(v[0]), *b = pmesh.GetVertex(v[1]);
      const bool pec = OnPEC(a) && OnPEC(b);
      const bool other = std::abs(a[1] - y_other) < eps && std::abs(b[1] - y_other) < eps &&
                         std::min(a[0], b[0]) > x_pec - eps;
      if (pec || other)
      {
        fes.GetEdgeInteriorDofs(e, dofs);
        Check(dofs, pec,
              edge_list.GetMeshIdAndType(e).type ==
                  mfem::NCMesh::NCList::MeshIdType::MASTER);
      }
    }
    for (int i = 0; i < pmesh.GetNV(); i++)
    {
      const double *x = pmesh.GetVertex(i);
      if (OnPEC(x) || OnOther(x))
      {
        fes.GetVertexDofs(i, dofs);
        Check(dofs, OnPEC(x), false);
      }
    }
    Mpi::GlobalSum(3, counts, comm);
    CHECK(counts[0] > 0);
    CHECK(counts[1] > 0);
    CHECK((counts[2] > 0) == edge_dofs);
  };
  CheckDbcTDofs(mode_op.GetNDSpace(), mode_op.GetNDDbcTDofLists().back(), true);
  CheckDbcTDofs(mode_op.GetH1Space(), mode_op.GetH1DbcTDofLists().back(), order > 1);
}

// BoundaryMode cross-section extracted from a nonconforming 3D mesh: the wave port of
// cpw_wave_2dmode on the uncracked mesh, with the elements next to the port refined on
// one side of the traces (as an adapted mesh can be), with straight or curved elements.
// The 2D mesh must keep its boundary elements through partitioning and uniform refinement,
// the sides of the traces must be identified for error estimation, and the DoFs on the PEC
// edges (master edges included) must be essential, and no others.
TEST_CASE_METHOD(palace::test::SharedTempDir,
                 "BoundaryMode cross-section of a nonconforming 3D mesh",
                 "[boundarymodeoperator][Serial][Parallel]")
{
  const bool curved = GENERATE(false, true);
  MPI_Comm comm = Mpi::World();
  const auto cpw_dir = fs::path(PALACE_TEST_DATA_DIR) / "regression" / "input" / "cpw";
  auto buffer = PreprocessFile((cpw_dir / "cpw_wave_2dmode.json").string().c_str());
  auto setup = nlohmann::json::parse(buffer);
  setup["Model"]["Mesh"] = (cpw_dir / setup["Model"]["Mesh"].get<std::string>()).string();
  setup["Model"]["CrackInternalBoundaryElements"] = false;
  setup["Model"]["Refinement"] = {{"UniformLevels", 1}};
  setup["Problem"]["Output"] = temp_dir.string();
  IoData iodata(setup, /*print=*/false);
  const std::vector<int> pec = iodata.boundaries.pec.attributes;

  auto smesh = mesh::Load(iodata, comm);
  if (smesh)
  {
    // Refine the elements next to the cross-section on one side of the traces, the
    // interior PEC boundary elements.
    double z_trace = mfem::infinity();
    for (int be = 0; be < smesh->GetNBE(); be++)
    {
      int e1, e2;
      smesh->GetFaceElements(smesh->GetBdrElementFaceIndex(be), &e1, &e2);
      if (e2 >= 0 &&
          std::find(pec.begin(), pec.end(), smesh->GetBdrAttribute(be)) != pec.end())
      {
        mfem::Array<int> v;
        smesh->GetBdrElementVertices(be, v);
        z_trace = std::min(z_trace, smesh->GetVertex(v[0])[2]);
      }
    }
    REQUIRE(z_trace < mfem::infinity());
    if (curved)
    {
      // Curve the elements in the plane of the cross-section with a smooth displacement
      // along y and z, which keeps the planes x = const and the exterior boundaries. Its
      // direction varies, such that the Jacobian determinant of the elements is quadratic.
      mfem::Vector min, max;
      smesh->GetBoundingBox(min, max);
      smesh->SetCurvature(2);
      auto &nodes = *smesh->GetNodes();
      const auto &fes = *nodes.FESpace();
      for (int i = 0; i < fes.GetNDofs(); i++)
      {
        const int iy = fes.DofToVDof(i, 1), iz = fes.DofToVDof(i, 2);
        const double sy = (nodes(iy) - min(1)) / (max(1) - min(1)),
                     sz = (nodes(iz) - min(2)) / (max(2) - min(2));
        const double sin_z = std::sin(std::numbers::pi * sz);
        nodes(iy) += 0.02 * (max(1) - min(1)) * std::sin(std::numbers::pi * sy) * sin_z;
        nodes(iz) +=
            0.02 * (max(2) - min(2)) * std::sin(2.0 * std::numbers::pi * sy) * sin_z;
      }
    }
    const auto &attrs = iodata.solver.boundary_mode.attributes;
    smesh->EnsureNCMesh(true);
    mfem::Array<int> marked;
    mfem::Vector c;
    for (int be = 0; be < smesh->GetNBE(); be++)
    {
      if (std::find(attrs.begin(), attrs.end(), smesh->GetBdrAttribute(be)) == attrs.end())
      {
        continue;
      }
      int e1, e2;
      smesh->GetFaceElements(smesh->GetBdrElementFaceIndex(be), &e1, &e2);
      smesh->GetElementCenter(e1, c);
      if (c(2) > z_trace)
      {
        marked.Append(e1);
      }
    }
    marked.Sort();
    marked.Unique();
    REQUIRE(marked.Size() > 0);
    smesh->GeneralRefinement(marked, 1, 0);
  }
  BoundaryModeSolver solver(iodata, Mpi::Root(comm), Mpi::Size(comm), 1, nullptr);
  solver.Preprocess(iodata, smesh, comm);
  std::vector<std::unique_ptr<mfem::ParMesh>> mfem_mesh;
  mfem_mesh.push_back(mesh::Partition(iodata, std::move(smesh), comm));

  // Segments of the PEC boundary elements of the 2D mesh before its uniform refinement,
  // gathered from all ranks (the other wave ports are relabelled as PEC), and the numbers
  // of PEC boundary elements and of those in the interior.
  auto CountPEC = [&](const mfem::ParMesh &pmesh, std::vector<double> *segs)
  {
    std::array<int, 2> n = {0, 0};
    for (int be = 0; be < pmesh.GetNBE(); be++)
    {
      if (std::find(pec.begin(), pec.end(), pmesh.GetBdrAttribute(be)) == pec.end())
      {
        continue;
      }
      n[0]++;
      n[1] += pmesh.FaceIsInterior(pmesh.GetBdrElementFaceIndex(be));
      mfem::Array<int> v;
      pmesh.GetBdrElementVertices(be, v);
      for (int i : v)
      {
        double x[2];
        pmesh.GetNode(i, x);
        if (segs)
        {
          segs->insert(segs->end(), x, x + 2);
        }
      }
    }
    Mpi::GlobalSum(2, n.data(), comm);
    return n;
  };
  std::vector<double> segs;
  const auto ne = mfem_mesh.back()->GetGlobalNE();
  const auto nbe_pec = CountPEC(*mfem_mesh.back(), &segs);
  CHECK(nbe_pec[0] > 0);
  CHECK(nbe_pec[1] > 0);
  mesh::RefineMesh(iodata, mfem_mesh);
  std::vector<std::unique_ptr<Mesh>> mesh;
  for (auto &m : mfem_mesh)
  {
    mesh.push_back(std::make_unique<Mesh>(std::move(m)));
  }
  const auto &pmesh = mesh.back()->Get();
  REQUIRE(pmesh.Nonconforming());
  REQUIRE(pmesh.GetGlobalNE() == 4 * ne);
  CHECK(CountPEC(pmesh, nullptr) == std::array<int, 2>{2 * nbe_pec[0], 2 * nbe_pec[1]});
  bool any_crack = mesh.back()->GetCrackSides(iodata.boundaries.GetSheetAttributes()).Any();
  Mpi::GlobalOr(1, &any_crack, comm);
  CHECK(any_crack);
  MaterialOperator mat_op(iodata, *mesh.back());
  BoundaryModeOperator mode_op(iodata, mesh, mat_op);

  int nloc = static_cast<int>(segs.size()), nproc = Mpi::Size(comm);
  std::vector<int> sizes(nproc), displs(nproc, 0);
  MPI_Allgather(&nloc, 1, MPI_INT, sizes.data(), 1, MPI_INT, comm);
  for (int r = 1; r < nproc; r++)
  {
    displs[r] = displs[r - 1] + sizes[r - 1];
  }
  std::vector<double> all_segs(displs.back() + sizes.back());
  MPI_Allgatherv(segs.data(), nloc, MPI_DOUBLE, all_segs.data(), sizes.data(),
                 displs.data(), MPI_DOUBLE, comm);
  auto OnPEC = [&](const double *x)
  {
    for (std::size_t k = 0; k < all_segs.size(); k += 4)
    {
      const double *a = &all_segs[k], *b = &all_segs[k + 2];
      const double ab2 = (b[0] - a[0]) * (b[0] - a[0]) + (b[1] - a[1]) * (b[1] - a[1]);
      const double ax_ab = (x[0] - a[0]) * (b[0] - a[0]) + (x[1] - a[1]) * (b[1] - a[1]);
      const double ax2 = (x[0] - a[0]) * (x[0] - a[0]) + (x[1] - a[1]) * (x[1] - a[1]);
      const double t = ax_ab / ab2;
      if (t > -1.0e-9 && t < 1.0 + 1.0e-9 && ax2 - t * t * ab2 < 1.0e-12 * ab2)
      {
        return true;
      }
    }
    return false;
  };

  // Numbers of checked true DoFs on the PEC edges and off them, and on master edges.
  int counts[3] = {0, 0, 0};
  auto CheckDbcTDofs =
      [&](const FiniteElementSpace &fespace, const mfem::Array<int> &dbc_tdof_list)
  {
    const auto &fes = fespace.Get();
    std::vector<bool> is_dbc(fes.GetTrueVSize(), false);
    for (auto t : dbc_tdof_list)
    {
      is_dbc[t] = true;
    }
    auto Check = [&](const mfem::Array<int> &dofs, bool on_pec, bool master)
    {
      for (auto d : dofs)
      {
        const int t = fes.GetLocalTDofNumber((d >= 0) ? d : -1 - d);
        if (t < 0)
        {
          continue;
        }
        CHECK(is_dbc[t] == on_pec);
        counts[on_pec ? 0 : 1]++;
        counts[2] += (on_pec && master);
      }
    };
    const auto &edge_list = pmesh.ncmesh->GetEdgeList();
    mfem::Array<int> v, dofs;
    for (int e = 0; e < pmesh.GetNEdges(); e++)
    {
      double x0[2], x1[2];
      pmesh.GetEdgeVertices(e, v);
      pmesh.GetNode(v[0], x0);
      pmesh.GetNode(v[1], x1);
      const double mid[2] = {0.5 * (x0[0] + x1[0]), 0.5 * (x0[1] + x1[1])};
      fes.GetEdgeInteriorDofs(e, dofs);
      Check(dofs, OnPEC(mid),
            edge_list.GetMeshIdAndType(e).type == mfem::NCMesh::NCList::MeshIdType::MASTER);
    }
    for (int i = 0; i < pmesh.GetNV(); i++)
    {
      double x[2];
      pmesh.GetNode(i, x);
      fes.GetVertexDofs(i, dofs);
      Check(dofs, OnPEC(x), false);
    }
  };
  if (!curved)
  {
    // The PEC edges are straight segments.
    CheckDbcTDofs(mode_op.GetNDSpace(), mode_op.GetNDDbcTDofLists().back());
    CheckDbcTDofs(mode_op.GetH1Space(), mode_op.GetH1DbcTDofLists().back());
    Mpi::GlobalSum(3, counts, comm);
    CHECK(counts[0] > 0);
    CHECK(counts[1] > 0);
    CHECK(counts[2] > 0);
  }
}

// BoundaryMode cross-section of a 3D mesh crossed by an impedance, conductivity, or
// rational impedance sheet, or by a sheet which is PEC only for the 2D solve (WavePortPEC):
// the uncracked metal trace of cpw_wave_2dmode crossing its wave port face. Where the sheet
// meets the cross-section, the 2D mesh has boundary elements with the attribute of the
// sheet, the same as for a PEC trace, on which the 2D solve applies the boundary condition
// of the sheet.
TEST_CASE_METHOD(palace::test::SharedTempDir,
                 "BoundaryMode cross-section edges of impedance sheets",
                 "[boundarymodeoperator][Serial][Parallel]")
{
  MPI_Comm comm = Mpi::World();
  constexpr int trace_attr = 11;
  const auto cpw_dir = fs::path(PALACE_TEST_DATA_DIR) / "regression" / "input" / "cpw";
  const std::string config =
      PreprocessFile((cpw_dir / "cpw_wave_2dmode.json").string().c_str()).str();

  // Numbers of boundary elements of the 2D mesh with the attribute of the trace, and of
  // those in its interior, with the given boundary condition on the trace.
  auto CountTraceEdges = [&](const std::string &type, const nlohmann::json &trace_bc)
  {
    auto setup = nlohmann::json::parse(config);
    setup["Model"]["Mesh"] = (cpw_dir / setup["Model"]["Mesh"].get<std::string>()).string();
    setup["Model"]["CrackInternalBoundaryElements"] = false;
    setup["Problem"]["Output"] = temp_dir.string();
    auto &boundaries = setup["Boundaries"];
    REQUIRE(boundaries["PEC"]["Attributes"] == nlohmann::json({8, 9, trace_attr}));
    if (type != "PEC")
    {
      boundaries["PEC"]["Attributes"] = {8, 9};
      nlohmann::json bc = trace_bc;
      bc["Attributes"] = nlohmann::json::array({trace_attr});
      // WavePortPEC is one boundary, the other types are arrays of boundaries.
      boundaries[type] = (type == "WavePortPEC") ? bc : nlohmann::json::array({bc});
    }
    IoData iodata(setup, /*print=*/false);
    auto smesh = mesh::Load(iodata, comm);
    BoundaryModeSolver solver(iodata, Mpi::Root(comm), Mpi::Size(comm), 1, nullptr);
    solver.Preprocess(iodata, smesh, comm);
    auto pmesh = mesh::Partition(iodata, std::move(smesh), comm);
    std::array<int, 2> n = {0, 0};
    for (int be = 0; be < pmesh->GetNBE(); be++)
    {
      if (pmesh->GetBdrAttribute(be) == trace_attr)
      {
        n[0]++;
        n[1] += pmesh->FaceIsInterior(pmesh->GetBdrElementFaceIndex(be));
      }
    }
    Mpi::GlobalSum(2, n.data(), comm);
    return n;
  };
  const auto pec = CountTraceEdges("PEC", {});
  CHECK(pec[1] > 0);
  const std::pair<std::string, nlohmann::json> sheets[] = {
      {"Impedance", {{"Rs", 1.0e-3}, {"Ls", 1.0e-12}}},
      {"Conductivity", {{"Conductivity", 5.0e7}}},
      {"RationalImpedance", {{"Numerator", {1.0e-12, 1.0e-3}}, {"Denominator", {1.0}}}},
      {"WavePortPEC", {}}};
  for (const auto &[type, bc] : sheets)
  {
    CAPTURE(type);
    CHECK(CountTraceEdges(type, bc) == pec);
  }
}

TEST_CASE("Mode target distance", "[boundarymodeoperator][Serial]")
{
  using mode_assembly::TargetDistance;

  // Reported wave-port selections of a strongly evanescent mode of a lossy cross-section
  // whose real part lies closer to the target than that of the guided mode, so ranking by
  // |Re{kn} - kn_target| picked it: a CPW with absorbing and conductivity boundaries at
  // 1 GHz (awslabs/palace#920, in 1/m) and examples/cpw/cpw_wave_uniform.json at 0.1 GHz
  // after a sweep to 50 GHz (awslabs/palace#996, nondimensional).
  struct ReportedSelection
  {
    int issue;
    double kn_target;
    std::complex<double> kn_guided, kn_evanescent;
  };
  const ReportedSelection reported[] = {
      {920, 60.2, {38.07, -4.445}, {51.03, 5677.0}},
      {996, 2.981703e-2, {1.988221e-2, -5.136255e-7}, {3.341085e-2, 1.574046e1}}};
  for (const auto &r : reported)
  {
    CAPTURE(r.issue);
    REQUIRE(std::abs(r.kn_evanescent.real() - r.kn_target) <
            std::abs(r.kn_guided.real() - r.kn_target));
    CHECK(TargetDistance(r.kn_guided, r.kn_target) <
          TargetDistance(r.kn_evanescent, r.kn_target));
  }
  const double kn_target = reported[1].kn_target;
  const std::complex<double> kn_cpw = reported[1].kn_guided;
  const std::complex<double> kn_evanescent = reported[1].kn_evanescent;

  // Lossless propagating modes below the target rank by descending kn, ahead of all
  // lossless evanescent modes, which rank by ascending attenuation.
  const double kt = 1.0;
  CHECK(TargetDistance(0.9 * kt, kt) < TargetDistance(0.5 * kt, kt));
  CHECK(TargetDistance(0.5 * kt, kt) < TargetDistance(1.0e-3 * kt, kt));
  CHECK(TargetDistance(1.0e-3 * kt, kt) <
        TargetDistance(std::complex<double>(0.0, 1.0e-3 * kt), kt));
  CHECK(TargetDistance(std::complex<double>(0.0, 1.0e-3 * kt), kt) <
        TargetDistance(std::complex<double>(0.0, 10.0 * kt), kt));

  // Invariant under conjugation, so the same for the two roots ±iα of a lossless mode below
  // cutoff.
  for (const auto kn : {kn_cpw, kn_evanescent, std::complex<double>(0.3, -0.7)})
  {
    CHECK(TargetDistance(std::conj(kn), kn_target) == TargetDistance(kn, kn_target));
  }

  // Non-finite values rank last, after any finite candidate.
  const double inf = std::numeric_limits<double>::infinity();
  const double nan = std::numeric_limits<double>::quiet_NaN();
  CHECK(TargetDistance(std::complex<double>(nan, 0.0), kt) == inf);
  CHECK(TargetDistance(std::complex<double>(0.0, nan), kt) == inf);
  CHECK(TargetDistance(std::complex<double>(nan, inf), kt) == inf);
  CHECK(TargetDistance(std::complex<double>(inf, 0.0), kt) == inf);
  CHECK(TargetDistance(kn_evanescent, kt) < inf);
}

TEST_CASE("Mode propagation constant branch", "[boundarymodeoperator][Serial]")
{
  using mode_assembly::PropagationConstant;

  // Lossless modes, with the roundoff of the mode solve in Im{kn^2} of either sign: the
  // root which propagates in the +n direction above cutoff, and kn = -i alpha, which decays
  // in the +n direction, below cutoff.
  for (const double roundoff : {-1.0e-14, 0.0, 1.0e-14})
  {
    CAPTURE(roundoff);
    const auto kn_above = PropagationConstant({4.0, roundoff});
    const auto kn_below = PropagationConstant({-4.0, roundoff});
    CHECK(std::abs(kn_above - 2.0) <= 1.0e-14);
    CHECK(std::abs(kn_below - std::complex<double>(0.0, -2.0)) <= 1.0e-14);
    CHECK(kn_below.imag() < 0.0);
  }

  // At a complex frequency, the analytic continuation of the root at the real frequency
  // Re{omega}: kn = n omega for a TEM mode, with kn^2 = n^2 omega^2, also for |Im{omega}| >
  // Re{omega}, where Re{kn^2} < 0, and, for a mode with kn^2 = n^2 omega^2 - kc^2, the
  // continuation of the root decaying below cutoff, or of the propagating one above.
  const double n = 1.5, kc = 3.0;
  for (const auto omega : {std::complex<double>(1.0, 0.1), std::complex<double>(1.0, 2.0),
                           std::complex<double>(1.0, 10.0)})
  {
    CAPTURE(omega);
    CHECK(std::abs(PropagationConstant(n * n * omega * omega, omega) - n * omega) <=
          1.0e-14 * std::abs(n * omega));
  }
  for (const double omega_r : {1.0, 1.9, 2.1, 3.0})
  {
    CAPTURE(omega_r);
    auto kn2 = [&](std::complex<double> omega) { return n * n * omega * omega - kc * kc; };
    const auto kn_real = PropagationConstant(kn2(omega_r), omega_r);
    CHECK(((omega_r < kc / n) ? kn_real.imag() < 0.0 : kn_real.real() > 0.0));
    for (const double omega_i : {1.0e-3, 1.0e-1, 1.0})
    {
      CAPTURE(omega_i);
      // Continuous along omega_r + i t, for t from 0 to omega_i.
      std::complex<double> kn_prev = kn_real;
      for (int k = 1; k <= 1000; k++)
      {
        const std::complex<double> omega(omega_r, omega_i * k / 1000.0);
        const auto kn = PropagationConstant(kn2(omega), omega);
        CHECK(std::abs(kn - kn_prev) <= 0.1 * std::abs(kn_prev));
        kn_prev = kn;
      }
    }
  }

  // Passive lossy modes, with Im{kn^2} < 0, above, at, and below cutoff: the principal
  // root, with Re{kn} >= 0 and Im{kn} <= 0, continuous through cutoff.
  for (const auto kn2 :
       {std::complex<double>(1.0, -0.5), std::complex<double>(1.0e-20, -0.5),
        std::complex<double>(0.0, -0.5), std::complex<double>(-1.0e-20, -0.5),
        std::complex<double>(-1.0, -0.5)})
  {
    CAPTURE(kn2);
    const auto kn = PropagationConstant(kn2);
    CHECK(std::abs(kn - std::sqrt(kn2)) <= 1.0e-14 * std::abs(kn));
    CHECK(kn.real() >= 0.0);
    CHECK(kn.imag() <= 0.0);
    CHECK(std::abs(kn * kn - kn2) <= 1.0e-14 * std::abs(kn2));
  }
}

TEST_CASE("ModeEigenSolver takes the decaying root below cutoff",
          "[boundarymodeoperator][Serial][Parallel]")
{
  // Below the TE10 cutoff (75 GHz for this 1000 x 500 um guide with eps = 4) every mode of
  // the lossless guide is evanescent, with kn^2 < 0 up to the roundoff of the mode solve,
  // whose sign selected the root returned by the principal square root: kn = -i alpha,
  // which decays in the +n direction, or the growing kn = +i alpha.
  for (const double freq_ghz : {20.0, 35.0, 50.0, 65.0})
  {
    auto result =
        SolveRectangularModes(1000.0, 500.0, freq_ghz, 4.0, 2, 3, [](IoData &) {});
    REQUIRE(result.num_converged >= 3);
    for (const auto &kn : result.kn)
    {
      CAPTURE(freq_ghz, kn);
      CHECK(kn.imag() < 0.0);
      CHECK(std::abs(kn.real()) <= 1.0e-6 * std::abs(kn));
    }
  }
}

TEST_CASE("ModeEigenSolver ranks evanescent modes below cutoff by attenuation",
          "[boundarymodeoperator][Serial][Parallel]")
{
  // Below the TE10 cutoff (75 GHz for this 1000 x 500 um guide with eps = 4) every mode of
  // a lossless guide is evanescent: the real parts of the computed kn are round-off and
  // carry no ranking information, and ranking on them returned an arbitrary evanescent mode
  // as mode 1 (e.g. |S11| > 0 dB for a lossless shorted guide driven below cutoff). Both
  // the exact and the reduced solve must rank the modes least-attenuated first, so mode 1
  // is the continuation of TE10 through cutoff.
  constexpr int num_modes = 3;
  constexpr double freq_ghz = 60.0, epsilon_r = 4.0;
  auto exact = SolveRectangularModes(1000.0, 500.0, freq_ghz, epsilon_r, 2, num_modes,
                                     [](IoData &) {});
  REQUIRE(exact.num_converged >= num_modes);
  for (int i = 0; i < exact.num_converged; i++)
  {
    CAPTURE(i, exact.kn[i]);
    REQUIRE(std::abs(exact.kn[i].real()) <= 1.0e-6 * std::abs(exact.kn[i].imag()));
    CHECK(exact.kn[i].imag() < 0.0);
    if (i > 0)
    {
      CHECK(std::abs(exact.kn[i - 1].imag()) <=
            (1.0 + 1.0e-6) * std::abs(exact.kn[i].imag()));
    }
  }

  // TE10 attenuation constant sqrt((pi/a)^2 - eps k0^2), nondimensionalized by 1 um.
  const double kc = std::numbers::pi / 1000.0;
  const double k = std::sqrt(epsilon_r) * 2.0 * std::numbers::pi * freq_ghz * 1.0e9 /
                   electromagnetics::c0_ * 1.0e-6;
  CHECK_THAT(std::abs(exact.kn[0].imag()), WithinRel(std::sqrt(kc * kc - k * k), 1.0e-3));

  // Train the reduced model below cutoff as well (0.9 and 1.1 times the frequency). The
  // reduced solve takes the same decaying root as the exact one.
  auto reduced = SolveRectangularModes(
      1000.0, 500.0, freq_ghz, epsilon_r, 2, num_modes, [](IoData &) {}, true);
  REQUIRE(reduced.reduced_stats.reduced_solves == 1);
  REQUIRE(reduced.kn.size() == num_modes);
  for (int i = 0; i < num_modes; i++)
  {
    CAPTURE(i, reduced.kn[i], exact.kn[i]);
    CHECK(std::abs(reduced.kn[i] - exact.kn[i]) <= 1.0e-6 * std::abs(exact.kn[i]));
  }
}

TEST_CASE("ModeEigenSolver ranks lossy modes by complex target distance",
          "[boundarymodeoperator][Serial][Parallel]")
{
  // A resistive top wall gives the modes attenuation comparable to the spacing of their
  // phase constants, so ranking by |Re{kn} - kn_target| and by |kn - kn_target| disagree.
  // Exact and reduced solves must both return the modes in ascending complex distance.
  auto configure_resistive_wall = [](IoData &iodata)
  {
    iodata.boundaries.pec.attributes = {1, 2, 4};
    auto &imp = iodata.boundaries.impedance.emplace_back();
    imp.attributes = {3};
    imp.Rs = 377.0;
  };
  constexpr int num_modes = 3;
  auto exact = SolveRectangularModes(1000.0, 500.0, 200.0, 4.0, 2, num_modes,
                                     configure_resistive_wall);
  REQUIRE(exact.num_converged >= num_modes);
  const double kn_target = exact.kn_target;
  for (int i = 1; i < exact.num_converged; i++)
  {
    CAPTURE(i, exact.kn[i - 1], exact.kn[i]);
    CHECK(mode_assembly::TargetDistance(exact.kn[i - 1], kn_target) <=
          mode_assembly::TargetDistance(exact.kn[i], kn_target));
  }

  // Guard the premise: some returned mode has a real part closer to the target than the
  // first-ranked mode, so a real-part ranking would select differently.
  bool real_part_ranking_differs = false;
  for (int i = 1; i < exact.num_converged; i++)
  {
    real_part_ranking_differs =
        real_part_ranking_differs ||
        std::abs(exact.kn[i].real() - kn_target) < std::abs(exact.kn[0].real() - kn_target);
  }
  REQUIRE(real_part_ranking_differs);

  // The reduced wave-port solve, trained with all ranked modes around the target, must
  // select the same modes in the same order.
  auto reduced = SolveRectangularModes(1000.0, 500.0, 200.0, 4.0, 2, num_modes,
                                       configure_resistive_wall, true, 1, 16,
                                       {0.95, 0.98, 1.02, 1.05}, 1.0e-2);
  REQUIRE(reduced.reduced_stats.reduced_solves == 1);
  REQUIRE(reduced.kn.size() == num_modes);
  for (int i = 0; i < num_modes; i++)
  {
    CAPTURE(i, reduced.kn[i], exact.kn[i]);
    CHECK(std::abs(reduced.kn[i] - exact.kn[i]) <= 1.0e-5 * std::abs(exact.kn[i]));
  }
}

// The Poynting power of a boundary mode is the time average of its peak phasors,
// P = 1/2 ∫ (Eₜ x Hₜ⋆) ⋅ ẑ with Hₜ = μ⁻¹ kₙ ẑ x Eₜ / ω for Eₙ = 0: a uniform field Eₜ = e x̂
// in vacuum on the unit square carries P = 1/2 (kₙ⋆ / ω) |e|² (nondimensional, μ = 1). The
// boundary mode solver normalizes its modes to |P| = 1, so this weight makes them carry 1
// W.
TEST_CASE("Boundary mode Poynting power is the time average",
          "[boundarymodeoperator][Serial][Parallel]")
{
  const auto comm = Mpi::World();
  IoData iodata(Units(0.013, 0.43));
  iodata.problem.type = ProblemType::BOUNDARYMODE;
  iodata.solver.order = 1;
  iodata.solver.boundary_mode.freq = 2.0;
  iodata.solver.boundary_mode.n = 1;
  iodata.domains.materials.emplace_back().attributes = {1};
  iodata.CheckConfiguration();
  auto serial = mfem::Mesh::MakeCartesian2D(2, 2, mfem::Element::QUADRILATERAL);
  serial.EnsureNodes();
  std::vector<std::unique_ptr<Mesh>> meshes;
  meshes.push_back(std::make_unique<Mesh>(std::make_unique<mfem::ParMesh>(comm, serial)));
  MaterialOperator mat_op(iodata, *meshes.back());
  BoundaryModeOperator mode_op(iodata, meshes, mat_op);
  constexpr double e = 0.7, omega = 2.0;
  const std::complex<double> kn(1.5, -0.2);
  mfem::Vector field(2);
  field = 0.0;
  field(0) = e;
  mfem::VectorConstantCoefficient coeff(field);
  GridFunction E(mode_op.GetNDSpace(), true);
  E.Real().ProjectCoefficient(coeff);
  E.Imag() = 0.0;
  ComplexVector et(mode_op.GetNDTrueVSize()), en(mode_op.GetH1TrueVSize());
  E.Real().GetTrueDofs(et.Real());
  et.Imag() = 0.0;
  en = 0.0;
  const std::complex<double> P = mode_op.ComputePoyntingPower(omega, kn, et, en);
  const std::complex<double> expected = 0.5 * std::conj(kn) / omega * e * e;
  CHECK_THAT(P.real(), WithinRel(expected.real(), 1.0e-12));
  CHECK_THAT(P.imag(), WithinRel(expected.imag(), 1.0e-12));
}

TEST_CASE("BoundaryModeOperator rejects multiple boundary conditions on an attribute",
          "[boundarymodeoperator][Serial][Parallel]")
{
  // As in the 3D operators, no boundary attribute may have more than one of the boundary
  // conditions of the 2D mode solve, whose terms would otherwise all be added on it.
  // MakeCartesian2D boundary attributes: bottom = 1, right = 2, top = 3, left = 4.
  auto Construct = [](const std::function<void(IoData &)> &configure_bcs)
  {
    MPI_Comm comm = Mpi::World();
    IoData iodata(Units(1.0, 1.0));
    iodata.problem.type = ProblemType::BOUNDARYMODE;
    iodata.model.Lc = 1.0;
    iodata.domains.materials.emplace_back().attributes = {1};
    iodata.boundaries.pec.attributes = {1, 3};
    configure_bcs(iodata);
    iodata.solver.order = 1;
    iodata.solver.boundary_mode.freq = 1.0;
    iodata.solver.boundary_mode.n = 1;
    iodata.CheckConfiguration();
    auto serial = mfem::Mesh::MakeCartesian2D(2, 2, mfem::Element::QUADRILATERAL);
    std::vector<std::unique_ptr<Mesh>> mesh;
    mesh.push_back(std::make_unique<Mesh>(std::make_unique<mfem::ParMesh>(comm, serial)));
    MaterialOperator mat_op(iodata, *mesh.back());
    BoundaryModeOperator mode_op(iodata, mesh, mat_op);
  };
  auto Impedance = [](IoData &iodata, std::vector<int> attributes)
  {
    auto &data = iodata.boundaries.impedance.emplace_back();
    data.attributes = std::move(attributes);
    data.Rs = 50.0;
  };
  auto Conductivity = [](IoData &iodata, std::vector<int> attributes)
  {
    auto &data = iodata.boundaries.conductivity.emplace_back();
    data.attributes = std::move(attributes);
    data.sigma = 5.8e7;
  };
  auto RationalImpedance = [](IoData &iodata, std::vector<int> attributes)
  {
    auto &data = iodata.boundaries.rational_impedance.emplace_back();
    data.attributes = std::move(attributes);
    data.num = {1.0e-12, 50.0};
    data.den = {1.0};
  };
  const auto multiple_bc = Catch::Matchers::ContainsSubstring(
      "Boundary attribute 2 should not be specified with multiple BC!");

  // One boundary condition on each attribute, and the same PEC attribute given as "PEC"
  // and "WavePortPEC".
  CHECK_NOTHROW(Construct(
      [&](IoData &iodata)
      {
        iodata.boundaries.auxpec.attributes = {1};
        Impedance(iodata, {2});
        iodata.boundaries.farfield.attributes = {4};
      }));
  CHECK_NOTHROW(Construct(
      [&](IoData &iodata)
      {
        Conductivity(iodata, {2});
        RationalImpedance(iodata, {4});
      }));

  // The same attribute in two boundary conditions.
  const std::vector<std::function<void(IoData &)>> configurations = {
      [&](IoData &iodata)
      {
        Impedance(iodata, {2});
        iodata.boundaries.farfield.attributes = {2};
      },
      [&](IoData &iodata)
      {
        iodata.boundaries.auxpec.attributes = {2, 4};
        Impedance(iodata, {2});
      },
      [&](IoData &iodata)
      {
        iodata.boundaries.pec.attributes = {2};
        Conductivity(iodata, {2, 4});
      },
      [&](IoData &iodata)
      {
        Conductivity(iodata, {2});
        RationalImpedance(iodata, {2});
      },
      [&](IoData &iodata)
      {
        RationalImpedance(iodata, {2});
        iodata.boundaries.farfield.attributes = {2, 4};
      }};
  for (const auto &configure_bcs : configurations)
  {
    CHECK_THROWS_WITH(Construct(configure_bcs), multiple_bc);
  }
}

}  // namespace palace
