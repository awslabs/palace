// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include <algorithm>
#include <cmath>
#include <complex>
#include <limits>
#include <vector>
#include <catch2/catch_test_macros.hpp>
#include <nlohmann/json.hpp>
#include <catch2/generators/catch_generators.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>

#include "fem/bilinearform.hpp"
#include "fem/integrator.hpp"
#include "fem/libceed/operator.hpp"
#include "fem/mesh.hpp"
#include "fem/qfunctions/coeff/coeff_qf.h"
#include "fixtures.hpp"
#include "linalg/hypre.hpp"
#include "linalg/rap.hpp"
#include "models/materialoperator.hpp"
#include "models/pml.hpp"
#include "models/spaceoperator.hpp"
#include "utils/communication.hpp"
#include "utils/configfile.hpp"
#include "utils/filesystem.hpp"
#include "utils/geodata.hpp"
#include "utils/iodata.hpp"
#include "utils/labels.hpp"
#include "utils/units.hpp"

namespace palace
{
namespace
{

struct IntegrationSettingsGuard
{
  int pa_order_threshold = BilinearForm::pa_order_threshold;
  int p_trial = fem::DefaultIntegrationOrder::p_trial;
  bool q_order_jac = fem::DefaultIntegrationOrder::q_order_jac;
  int q_order_extra_pk = fem::DefaultIntegrationOrder::q_order_extra_pk;
  int q_order_extra_qk = fem::DefaultIntegrationOrder::q_order_extra_qk;

  ~IntegrationSettingsGuard()
  {
    BilinearForm::pa_order_threshold = pa_order_threshold;
    fem::DefaultIntegrationOrder::p_trial = p_trial;
    fem::DefaultIntegrationOrder::q_order_jac = q_order_jac;
    fem::DefaultIntegrationOrder::q_order_extra_pk = q_order_extra_pk;
    fem::DefaultIntegrationOrder::q_order_extra_qk = q_order_extra_qk;
  }
};

struct ImaginaryHierarchyData
{
  std::vector<HYPRE_Int> coarse_row_offsets;
  std::vector<HYPRE_Int> coarse_columns;
  std::vector<double> coarse_values;
  CeedInt fine_suboperators = 0;
};

ImaginaryHierarchyData InspectImaginaryHierarchy(const ComplexOperator &op)
{
  const auto *mg_op = dynamic_cast<const ComplexMultigridOperator *>(&op);
  REQUIRE(mg_op);
  REQUIRE(mg_op->GetNumLevels() == 2);

  const auto *coarse_op =
      dynamic_cast<const ComplexParOperator *>(&mg_op->GetOperatorAtLevel(0));
  REQUIRE(coarse_op);
  const auto *coarse_imag =
      dynamic_cast<const hypre::HypreCSRMatrix *>(coarse_op->LocalOperator().Imag());
  REQUIRE(coarse_imag);
  hypre_CSRMatrixMigrate(*coarse_imag, HYPRE_MEMORY_HOST);

  ImaginaryHierarchyData data;
  data.coarse_row_offsets.assign(coarse_imag->GetI(),
                                 coarse_imag->GetI() + coarse_imag->Height() + 1);
  data.coarse_columns.assign(coarse_imag->GetJ(), coarse_imag->GetJ() + coarse_imag->NNZ());
  data.coarse_values.assign(coarse_imag->GetData(),
                            coarse_imag->GetData() + coarse_imag->NNZ());

  const auto *fine_op =
      dynamic_cast<const ComplexParOperator *>(&mg_op->GetOperatorAtLevel(1));
  REQUIRE(fine_op);
  const auto *fine_imag =
      dynamic_cast<const ceed::Operator *>(fine_op->LocalOperator().Imag());
  REQUIRE(fine_imag);
  for (std::size_t i = 0; i < fine_imag->Size(); i++)
  {
    CeedInt num_suboperators = 0;
    REQUIRE(CeedOperatorCompositeGetNumSub((*fine_imag)[i], &num_suboperators) == 0);
    data.fine_suboperators += num_suboperators;
  }
  return data;
}

// Load and nondimensionalize the mesh of a configuration, as main.cpp does.
std::vector<std::unique_ptr<Mesh>> LoadMesh(IoData &iodata, MPI_Comm comm)
{
  auto smesh = mesh::Load(iodata, comm);
  if (iodata.model.Lc <= 0.0)
  {
    iodata.model.Lc = mesh::ComputeReferenceLength(smesh, comm);
  }
  iodata.NondimensionalizeInputs(smesh);
  std::vector<std::unique_ptr<mfem::ParMesh>> mfem_mesh;
  mfem_mesh.push_back(mesh::Partition(iodata, std::move(smesh), comm));
  mesh::RefineMesh(iodata, mfem_mesh);
  std::vector<std::unique_ptr<Mesh>> mesh;
  for (auto &m : mfem_mesh)
  {
    mesh.push_back(std::make_unique<Mesh>(std::move(m)));
  }
  return mesh;
}

// Relative difference |y - y_ref| / |y_ref|.
double RelativeDifference(MPI_Comm comm, const ComplexVector &y, const ComplexVector &y_ref)
{
  ComplexVector diff(y.Size());
  diff.UseDevice(true);
  diff = y;
  diff -= y_ref;
  const double norm = linalg::Norml2(comm, y_ref);
  REQUIRE(norm > 0.0);
  return linalg::Norml2(comm, diff) / norm;
}

}  // namespace

TEST_CASE("SpaceOperator retains coarse support while omitting fine exact zeros",
          "[spaceoperator][Serial][Parallel]")
{
  using namespace std::complex_literals;

  MPI_Comm comm = Mpi::World();
  mfem::Mesh serial_mesh = mfem::Mesh::MakeCartesian3D(1, 1, 1, mfem::Element::TETRAHEDRON);
  while (serial_mesh.GetNE() < Mpi::Size(comm))
  {
    serial_mesh.UniformRefinement();
  }
  std::vector<std::unique_ptr<Mesh>> mesh;
  mesh.push_back(std::make_unique<Mesh>(comm, serial_mesh));

  IntegrationSettingsGuard settings_guard;
  config::SolverData solver;
  solver.order = 2;
  solver.pa_order_threshold = 2;
  solver.linear.mg_max_levels = 2;
  solver.linear.mg_coarsening = MultigridCoarsening::LINEAR;
  solver.linear.pc_mat_real = false;
  solver.linear.pc_mat_shifted = 0;
  BilinearForm::pa_order_threshold = solver.pa_order_threshold;
  fem::DefaultIntegrationOrder::p_trial = solver.order;
  fem::DefaultIntegrationOrder::q_order_jac = solver.q_order_jac;
  fem::DefaultIntegrationOrder::q_order_extra_pk = solver.q_order_extra;
  fem::DefaultIntegrationOrder::q_order_extra_qk = solver.q_order_extra;

  config::MaterialData material;
  material.attributes = {1};
  config::DomainData domains;
  domains.attributes = {1};
  domains.materials = {material};
  config::BoundaryData boundaries;
  Units units(1.0, 1.0);
  SpaceOperator space_op(solver, domains, boundaries, ProblemType::EIGENMODE, units, mesh);

  constexpr double omega = 2.0;
  const std::complex<double> lambda_zero = 1i * omega;
  const double epsilon = std::numeric_limits<double>::epsilon();
  const std::complex<double> lambda_tiny = epsilon + 1i * omega;
  auto Assemble = [&space_op](std::complex<double> lambda)
  {
    return space_op.GetPreconditionerMatrix<ComplexOperator>(1.0 + 0.0i, lambda,
                                                             lambda * lambda, lambda / 1i);
  };

  auto zero_pc = Assemble(lambda_zero);
  auto zero_data = InspectImaginaryHierarchy(*zero_pc);
  auto tiny_pc = Assemble(lambda_tiny);
  auto tiny_data = InspectImaginaryHierarchy(*tiny_pc);

  CHECK(zero_data.coarse_row_offsets == tiny_data.coarse_row_offsets);
  CHECK(zero_data.coarse_columns == tiny_data.coarse_columns);
  CHECK(std::ranges::all_of(zero_data.coarse_values,
                            [](double value) { return value == 0.0; }));
  CHECK(std::ranges::any_of(tiny_data.coarse_values,
                            [](double value) { return value != 0.0; }));
  CHECK(zero_data.fine_suboperators == 0);
  CHECK(tiny_data.fine_suboperators > 0);
}

TEST_CASE("CPU preconditioner quadrature data preserves the multigrid operators",
          "[spaceoperator][Serial][Parallel]")
{
  const auto element = GENERATE(mfem::Element::TETRAHEDRON, mfem::Element::HEXAHEDRON);
  const auto comm = Mpi::World();
  auto serial_mesh = mfem::Mesh::MakeCartesian3D(2, 2, 2, element);
  serial_mesh.SetCurvature(2);
  serial_mesh.Transform(
      [](const mfem::Vector &x, mfem::Vector &y)
      {
        y = x;
        y(0) += 0.05 * x(1) * x(2);
        y(1) += 0.03 * x(0) * x(0);
      });
  std::vector<std::unique_ptr<Mesh>> mesh;
  mesh.push_back(std::make_unique<Mesh>(comm, serial_mesh));

  IntegrationSettingsGuard settings_guard;
  config::SolverData solver;
  solver.order = 3;
  solver.pa_order_threshold = 2;
  solver.linear.mg_max_levels = 3;
  solver.linear.mg_coarsening = MultigridCoarsening::LINEAR;
  solver.linear.pc_mat_shifted = 0;
  BilinearForm::pa_order_threshold = solver.pa_order_threshold;
  fem::DefaultIntegrationOrder::p_trial = solver.order;
  fem::DefaultIntegrationOrder::q_order_jac = solver.q_order_jac;
  fem::DefaultIntegrationOrder::q_order_extra_pk = solver.q_order_extra;
  fem::DefaultIntegrationOrder::q_order_extra_qk = solver.q_order_extra;

  config::MaterialData material;
  material.attributes = {1};
  material.mu_r.s = {1.0, 2.0, 3.0};
  material.epsilon_r.s = {2.0, 4.0, 5.0};
  config::DomainData domains;
  domains.attributes = {1};
  domains.materials = {material};
  config::BoundaryData boundaries;
  Units units(1.0, 1.0);
  SpaceOperator space_op(solver, domains, boundaries, ProblemType::EIGENMODE, units, mesh);
  auto pc = space_op.GetPreconditionerMatrix<Operator>(1.5, 0.0, 2.0, 0.0);
  const auto *mg = dynamic_cast<const MultigridOperator *>(pc.get());
  REQUIRE(mg);
  REQUIRE(mg->GetNumLevels() == 3);

  const auto &mat = space_op.GetMaterialOp();
  MaterialPropertyCoefficient curl(mat.GetAttributeToMaterial(),
                                   mat.GetCurlCurlInvPermeability(), 1.5);
  MaterialPropertyCoefficient mass(mat.GetAttributeToMaterial(), mat.GetPermittivityAbs(),
                                   2.0);
  for (bool auxiliary : {false, true})
  {
    const auto &spaces = auxiliary ? space_op.GetH1Spaces() : space_op.GetNDSpaces();
    for (std::size_t l = 0; l < mg->GetNumLevels(); l++)
    {
      const auto &space = spaces.GetFESpaceAtLevel(l);
      // Independently assemble each level without cached quadrature data. This also
      // checks reuse of the cached tensors between polynomial levels.
      BilinearForm form(space);
      if (auxiliary)
      {
        form.AddDomainIntegrator<DiffusionIntegrator>(mass);
      }
      else
      {
        form.AddDomainIntegrator<CurlCurlMassIntegrator>(curl, mass);
      }
      ParOperator reference(form.Assemble(false), space);
      const auto &actual =
          auxiliary ? mg->GetAuxiliaryOperatorAtLevel(l) : mg->GetOperatorAtLevel(l);
      Vector x(space.GetTrueVSize()), expected(x.Size()), result(x.Size());
      for (int seed : {42, 314159})
      {
        x.Randomize(seed + Mpi::Rank(comm));
        reference.Mult(x, expected);
        actual.Mult(x, result);
        result.Add(-1.0, expected);
        CHECK_THAT(linalg::Norml2(comm, result) / linalg::Norml2(comm, expected),
                   Catch::Matchers::WithinAbs(0.0, 2.0e-13));
      }
    }
  }
}

TEST_CASE("SpaceOperator frequency-dependent PML at complex frequency",
          "[spaceoperator][pml][Serial][Parallel][GPU]")
{
  using namespace std::complex_literals;

  // Unit cube with a PML layer z ∈ [0.75, 1] (attribute 2).
  MPI_Comm comm = Mpi::World();
  mfem::Mesh serial_mesh = mfem::Mesh::MakeCartesian3D(4, 4, 4, mfem::Element::HEXAHEDRON);
  for (int i = 0; i < serial_mesh.GetNE(); i++)
  {
    mfem::Vector center;
    serial_mesh.GetElementCenter(i, center);
    serial_mesh.SetAttribute(i, (center(2) > 0.75) ? 2 : 1);
  }
  serial_mesh.SetAttributes();
  std::vector<std::unique_ptr<Mesh>> mesh;
  mesh.push_back(std::make_unique<Mesh>(comm, serial_mesh));

  IntegrationSettingsGuard settings_guard;
  config::SolverData solver;
  solver.order = 2;
  solver.linear.mg_max_levels = 1;
  solver.linear.pc_mat_real = false;
  solver.linear.pc_mat_shifted = 0;
  fem::DefaultIntegrationOrder::p_trial = solver.order;

  constexpr double omega0 = 2.0;
  auto MakeSpaceOperator = [&](bool frequency_dependent)
  {
    config::MaterialData vacuum, pml;
    vacuum.attributes = {1};
    pml.attributes = {2};
    pml.epsilon_r.s = {2.0, 2.0, 2.0};
    config::PMLData pml_pml;
    pml_pml.autodetect_geometry = false;
    pml_pml.directions = {false, false, false, false, false, true};
    pml_pml.thickness = {0.0, 0.0, 0.0, 0.0, 0.0, 0.25};
    pml_pml.kappa_max = {1.0, 1.0, 1.5};
    pml_pml.alpha_max = {0.0, 0.0, 0.2};
    pml_pml.frequency_dependent = frequency_dependent;
    pml_pml.reference_frequency = omega0;
    config::DomainData domains;
    domains.attributes = {1, 2};
    domains.materials = {vacuum, pml};
    domains.pml = {pml_pml};
    domains.pml[0].attributes = pml.attributes;
    config::BoundaryData boundaries;
    Units units(1.0, 1.0);
    return std::make_unique<SpaceOperator>(solver, domains, boundaries,
                                           ProblemType::EIGENMODE, units, mesh);
  };
  auto fd_op = MakeSpaceOperator(true);
  const ComplexOperator *C0 = nullptr;
  REQUIRE(fd_op->GetMaterialOp().HasPML());
  REQUIRE(fd_op->GetMaterialOp().HasFrequencyDependentPML());

  ComplexVector x(fd_op->GetNDSpace().GetTrueVSize()), y1(x.Size()), y2(x.Size());
  x.UseDevice(true);
  y1.UseDevice(true);
  y2.UseDevice(true);
  x.Real().Randomize(42 + Mpi::Rank(comm));
  x.Imag().Randomize(314159 + Mpi::Rank(comm));
  auto RelErr = [&](ComplexVector &a, const ComplexVector &b)
  {
    a.Add(-1.0, b);
    return linalg::Norml2(comm, a) / linalg::Norml2(comm, b);
  };
  auto A2 = [&](std::complex<double> omega)
  {
    auto op = fd_op->GetExtraSystemMatrix(omega, Operator::DIAG_ZERO);
    REQUIRE(op);
    return op;
  };

  SECTION("Complex frequency overload matches the real frequency one for real ω")
  {
    auto A2r = fd_op->GetExtraSystemMatrix<ComplexOperator>(omega0, Operator::DIAG_ZERO);
    REQUIRE(A2r);
    A2r->Mult(x, y1);
    A2(omega0)->Mult(x, y2);
    CHECK_THAT(RelErr(y2, y1), Catch::Matchers::WithinAbs(0.0, 1.0e-14));
  }

  SECTION("System matrix at real ω matches the static PML with ω₀ = ω")
  {
    auto static_op = MakeSpaceOperator(false);
    REQUIRE(!static_op->GetMaterialOp().HasFrequencyDependentPML());
    auto Ks = static_op->GetStiffnessMatrix<ComplexOperator>(Operator::DIAG_ZERO);
    auto Ms = static_op->GetMassMatrix<ComplexOperator>(Operator::DIAG_ZERO);
    auto As = static_op->GetSystemMatrix(1.0 + 0.0i, 0.0 + 0.0i, -omega0 * omega0 + 0.0i,
                                         Ks.get(), C0, Ms.get());
    auto Kf = fd_op->GetStiffnessMatrix<ComplexOperator>(Operator::DIAG_ZERO);
    auto Mf = fd_op->GetMassMatrix<ComplexOperator>(Operator::DIAG_ZERO);
    auto A2f = A2(omega0);
    auto Af = fd_op->GetSystemMatrix(1.0 + 0.0i, 0.0 + 0.0i, -omega0 * omega0 + 0.0i,
                                     Kf.get(), C0, Mf.get(), A2f.get());
    As->Mult(x, y1);
    Af->Mult(x, y2);
    CHECK_THAT(RelErr(y2, y1), Catch::Matchers::WithinAbs(0.0, 1.0e-13));
  }

  SECTION("Frozen-stretch PML matrices are the static PML matrices")
  {
    // The frozen-stretch matrices at ω₀ reproduce A2(ω₀) and the static PML terms of the
    // stiffness and mass matrices for the reference frequency ω₀.
    auto P = fd_op->GetFrozenPMLMatrices(omega0, Operator::DIAG_ZERO);
    REQUIRE(P[0]);
    REQUIRE(!P[1]);
    REQUIRE(P[2]);
    ComplexVector t(x.Size());
    t.UseDevice(true);
    P[0]->Mult(x, y1);
    P[2]->Mult(x, t);
    y1.Add(-omega0 * omega0, t);
    A2(omega0)->Mult(x, y2);
    CHECK_THAT(RelErr(y1, y2), Catch::Matchers::WithinAbs(0.0, 1.0e-14));

    auto static_op = MakeSpaceOperator(false);
    auto Ks = static_op->GetStiffnessMatrix<ComplexOperator>(Operator::DIAG_ZERO);
    auto Kf = fd_op->GetStiffnessMatrix<ComplexOperator>(Operator::DIAG_ZERO);
    Ks->Mult(x, y1);
    Kf->Mult(x, y2);
    P[0]->AddMult(x, y2);
    CHECK_THAT(RelErr(y2, y1), Catch::Matchers::WithinAbs(0.0, 1.0e-14));
    auto Ms = static_op->GetMassMatrix<ComplexOperator>(Operator::DIAG_ZERO);
    auto Mf = fd_op->GetMassMatrix<ComplexOperator>(Operator::DIAG_ZERO);
    Ms->Mult(x, y1);
    Mf->Mult(x, y2);
    P[2]->AddMult(x, y2);
    CHECK_THAT(RelErr(y2, y1), Catch::Matchers::WithinAbs(0.0, 1.0e-14));

    // Without the PML terms, A2 is empty.
    CHECK(!fd_op->GetExtraSystemMatrix(std::complex<double>(omega0, 0.3),
                                       Operator::DIAG_ZERO, false));
  }

  SECTION("A2(ω) is the analytic continuation in ω")
  {
    // Cauchy-Riemann equations for the operator action, with central differences.
    const std::complex<double> omega = {omega0, 0.3};
    const double h = 1.0e-5;
    ComplexVector t(x.Size());
    t.UseDevice(true);
    A2(omega + h)->Mult(x, y1);
    A2(omega - h)->Mult(x, t);
    y1.Add(-1.0, t);  // ∂A2/∂Re{ω} x (2 h)
    A2(omega + 1i * h)->Mult(x, y2);
    A2(omega - 1i * h)->Mult(x, t);
    y2.Add(-1.0, t);
    y2 *= -1i;  // -i ∂A2/∂Im{ω} x (2 h) = ∂A2/∂Re{ω} x (2 h)
    CHECK_THAT(RelErr(y2, y1), Catch::Matchers::WithinAbs(0.0, 1.0e-8));

    // The stretch is evaluated at the complex frequency, not its real part.
    A2(omega)->Mult(x, y1);
    A2(omega.real())->Mult(x, y2);
    CHECK(RelErr(y2, y1) > 1.0e-3);
  }

  SECTION("Preconditioner matches the system matrix at complex ω")
  {
    const std::complex<double> omega = {omega0, 0.3}, lambda = 1i * omega;
    auto K = fd_op->GetStiffnessMatrix<ComplexOperator>(Operator::DIAG_ZERO);
    auto M = fd_op->GetMassMatrix<ComplexOperator>(Operator::DIAG_ZERO);
    auto A2c = A2(omega);
    auto A = fd_op->GetSystemMatrix(1.0 + 0.0i, lambda, lambda * lambda, K.get(), C0,
                                    M.get(), A2c.get());
    auto P = fd_op->GetPreconditionerMatrix<ComplexOperator>(1.0 + 0.0i, lambda,
                                                             lambda * lambda, omega);
    const auto *mg = dynamic_cast<const ComplexMultigridOperator *>(P.get());
    REQUIRE(mg);
    A->Mult(x, y1);
    mg->GetFinestOperator().Mult(x, y2);
    CHECK_THAT(RelErr(y2, y1), Catch::Matchers::WithinAbs(0.0, 1.0e-12));
  }
}

TEST_CASE("SpaceOperator real-valued coarse solve with PML regions",
          "[spaceoperator][pml][Serial][Parallel][GPU]")
{
  using namespace std::complex_literals;

  // Unit cube (vacuum, lossless) with a strongly absorbing PML layer z ∈ [0.75, 1]
  // (attribute 2), with two multigrid levels.
  MPI_Comm comm = Mpi::World();
  mfem::Mesh serial_mesh = mfem::Mesh::MakeCartesian3D(4, 4, 4, mfem::Element::HEXAHEDRON);
  for (int i = 0; i < serial_mesh.GetNE(); i++)
  {
    mfem::Vector center;
    serial_mesh.GetElementCenter(i, center);
    serial_mesh.SetAttribute(i, (center(2) > 0.75) ? 2 : 1);
  }
  serial_mesh.SetAttributes();
  std::vector<std::unique_ptr<Mesh>> mesh;
  mesh.push_back(std::make_unique<Mesh>(comm, serial_mesh));

  IntegrationSettingsGuard settings_guard;
  constexpr double omega = 2.0;
  auto MakeSpaceOperator = [&](bool pc_mat_real, bool complex_coarse_solve)
  {
    config::SolverData solver;
    solver.order = 2;
    solver.linear.mg_max_levels = 2;
    solver.linear.pc_mat_real = pc_mat_real;
    solver.linear.pc_mat_shifted = 0;
    solver.linear.complex_coarse_solve = complex_coarse_solve;
    fem::DefaultIntegrationOrder::p_trial = solver.order;
    config::MaterialData vacuum;
    vacuum.attributes = {1, 2};
    config::PMLData pml;
    pml.attributes = {2};
    pml.reference_frequency = omega;
    pml.sigma_max = {0.0, 0.0, 100.0};
    config::DomainData domains;
    domains.attributes = {1, 2};
    domains.materials = {vacuum};
    domains.pml = {pml};
    config::BoundaryData boundaries;
    Units units(1.0, 1.0);
    return std::make_unique<SpaceOperator>(solver, domains, boundaries, ProblemType::DRIVEN,
                                           units, mesh);
  };
  auto Coarse = [](const std::unique_ptr<ComplexOperator> &P) -> const ComplexOperator &
  {
    const auto *mg = dynamic_cast<const ComplexMultigridOperator *>(P.get());
    REQUIRE(mg);
    REQUIRE(mg->GetNumLevels() == 2);
    return mg->GetOperatorAtLevel(0);
  };

  // Real-valued approximation of the preconditioner (PCMatReal).
  auto real_op = MakeSpaceOperator(true, false);
  auto P_real = real_op->GetPreconditionerMatrix<ComplexOperator>(
      1.0 + 0.0i, 0.0 + 0.0i, -omega * omega + 0.0i, omega);
  const ComplexOperator &B_real = Coarse(P_real);
  REQUIRE(B_real.Real());
  REQUIRE(!B_real.Imag());
  Vector x(B_real.Width()), y_real(B_real.Height());
  x.UseDevice(true);
  y_real.UseDevice(true);
  x.Randomize(42 + Mpi::Rank(comm));
  B_real.Real()->Mult(x, y_real);

  // With a real-valued coarse solve, the coarse-level operator of the complex
  // preconditioner has the real-valued approximation Re{c} |T| of the PML terms as its real
  // part, and no PML terms in its imaginary part (zero for this lossless problem): the sum
  // of the two parts factored by the coarse solver is the coarse-level operator of the
  // real-valued approximation.
  for (bool complex_coarse_solve : {false, true})
  {
    auto space_op = MakeSpaceOperator(false, complex_coarse_solve);
    auto P = space_op->GetPreconditionerMatrix<ComplexOperator>(
        1.0 + 0.0i, 0.0 + 0.0i, -omega * omega + 0.0i, omega);
    const ComplexOperator &B = Coarse(P);
    REQUIRE(B.Real());
    Vector y(B.Height()), t(B.Height());
    y.UseDevice(true);
    t.UseDevice(true);
    B.Real()->Mult(x, y);
    if (B.Imag())
    {
      B.Imag()->Mult(x, t);
      y += t;
    }
    t = y;
    t -= y_real;
    const double err = linalg::Norml2(comm, t) / linalg::Norml2(comm, y_real);
    INFO("ComplexCoarseSolve = " << complex_coarse_solve);
    if (complex_coarse_solve)
    {
      // The complex coarse-level operator has the exact PML terms.
      CHECK(err > 1.0e-2);
    }
    else
    {
      CHECK_THAT(err, Catch::Matchers::WithinAbs(0.0, 1.0e-12));
    }
  }
}

TEST_CASE("SpaceOperator PML blocks with static and frequency-dependent stretches",
          "[spaceoperator][pml][Serial][Parallel]")
{
  using namespace std::complex_literals;

  // Unit cube with PML layers z ∈ [0, 0.25] (attribute 2) and z ∈ [0.75, 1] (attribute 3).
  MPI_Comm comm = Mpi::World();
  mfem::Mesh serial_mesh = mfem::Mesh::MakeCartesian3D(4, 4, 4, mfem::Element::HEXAHEDRON);
  for (int i = 0; i < serial_mesh.GetNE(); i++)
  {
    mfem::Vector center;
    serial_mesh.GetElementCenter(i, center);
    serial_mesh.SetAttribute(i, (center(2) < 0.25) ? 2 : ((center(2) > 0.75) ? 3 : 1));
  }
  serial_mesh.SetAttributes();
  std::vector<std::unique_ptr<Mesh>> mesh;
  mesh.push_back(std::make_unique<Mesh>(comm, serial_mesh));

  IntegrationSettingsGuard settings_guard;
  config::SolverData solver;
  solver.order = 2;
  solver.linear.mg_max_levels = 1;
  solver.linear.pc_mat_real = false;
  solver.linear.pc_mat_shifted = 0;
  fem::DefaultIntegrationOrder::p_trial = solver.order;

  constexpr double omega0 = 2.0;
  auto MakeBlock = [](std::vector<int> attributes, int order, bool frequency_dependent)
  {
    config::PMLData pml;
    pml.attributes = std::move(attributes);
    pml.order = order;
    pml.kappa_max = {1.0, 1.0, 1.5};
    pml.frequency_dependent = frequency_dependent;
    pml.reference_frequency = omega0;
    return pml;
  };
  auto MakeSpaceOperator = [&](const std::vector<config::PMLData> &pml)
  {
    config::MaterialData vacuum, dielectric;
    vacuum.attributes = {1, 3};
    dielectric.attributes = {2};
    dielectric.epsilon_r.s = {2.0, 2.0, 2.0};
    config::DomainData domains;
    domains.attributes = {1, 2, 3};
    domains.materials = {vacuum, dielectric};
    domains.pml = pml;
    config::BoundaryData boundaries;
    Units units(1.0, 1.0);
    return std::make_unique<SpaceOperator>(solver, domains, boundaries, ProblemType::DRIVEN,
                                           units, mesh);
  };

  ComplexVector x, y1, y2;
  auto RelErr = [&](ComplexVector &a, const ComplexVector &b)
  {
    a.Add(-1.0, b);
    return linalg::Norml2(comm, a) / linalg::Norml2(comm, b);
  };
  auto Apply = [&](SpaceOperator &space_op, std::complex<double> omega, ComplexVector &y)
  {
    const ComplexOperator *C0 = nullptr;
    auto K = space_op.GetStiffnessMatrix<ComplexOperator>(Operator::DIAG_ZERO);
    auto M = space_op.GetMassMatrix<ComplexOperator>(Operator::DIAG_ZERO);
    auto A2 = space_op.GetExtraSystemMatrix(omega, Operator::DIAG_ZERO);
    auto A = space_op.GetSystemMatrix(1.0 + 0.0i, 0.0 + 0.0i, -omega * omega, K.get(), C0,
                                      M.get(), A2.get());
    if (x.Size() == 0)
    {
      x.SetSize(space_op.GetNDSpace().GetTrueVSize());
      x.UseDevice(true);
      x.Real().Randomize(42 + Mpi::Rank(comm));
      x.Imag().Randomize(314159 + Mpi::Rank(comm));
    }
    y.SetSize(x.Size());
    y.UseDevice(true);
    A->Mult(x, y);
    return A2 != nullptr;
  };

  // A static block with a different grading and a frequency-dependent block: at ω = ω₀,
  // the system matrix is the one of two static blocks.
  auto mixed = MakeSpaceOperator({MakeBlock({2}, 2, false), MakeBlock({3}, 3, true)});
  REQUIRE(mixed->GetMaterialOp().GetPMLLayers().size() == 2);
  REQUIRE(mixed->GetMaterialOp().HasFrequencyDependentPML());
  auto both_static =
      MakeSpaceOperator({MakeBlock({2}, 2, false), MakeBlock({3}, 3, false)});
  REQUIRE(!both_static->GetMaterialOp().HasFrequencyDependentPML());
  CHECK(Apply(*mixed, omega0, y1));
  CHECK(!Apply(*both_static, omega0, y2));
  CHECK_THAT(RelErr(y1, y2), Catch::Matchers::WithinAbs(0.0, 1.0e-13));

  // At another frequency, only the frequency-dependent block changes.
  Apply(*mixed, 1.5 * omega0, y1);
  Apply(*both_static, 1.5 * omega0, y2);
  CHECK(RelErr(y1, y2) > 1.0e-3);

  // Blocks with the same parameters are a single block (σ_max is given, since by default it
  // depends on the materials of each block).
  auto Block = [&](std::vector<int> attributes)
  {
    auto pml = MakeBlock(std::move(attributes), 3, false);
    pml.sigma_max = {0.0, 0.0, 20.0};
    return pml;
  };
  auto single = MakeSpaceOperator({Block({2, 3})});
  auto split = MakeSpaceOperator({Block({2}), Block({3})});
  REQUIRE(single->GetMaterialOp().GetPMLLayers().size() == 1);
  Apply(*single, omega0, y1);
  Apply(*split, omega0, y2);
  CHECK_THAT(RelErr(y1, y2), Catch::Matchers::WithinAbs(0.0, 1.0e-14));
}

TEST_CASE("SpaceOperator PML subdomain true DOFs are independent of the partitioning",
          "[spaceoperator][pml][Serial][Parallel]")
{
  using namespace std::complex_literals;

  // Box [0, 1] x [0, 1] x [0, 2] with a PML layer z > 1.5 (attribute 2), with the PML
  // elements on the first or last process and the others on the remaining processes, so
  // that the true DOFs on the interface with the PML layer are owned by a process with
  // only PML elements or only physical elements.
  MPI_Comm comm = Mpi::World();
  const int np = Mpi::Size(comm);
  constexpr int order = 2, nz = 8;
  if (np > nz - 1)
  {
    SKIP("Test requires at most " << nz - 1 << " processes");
  }
  const bool pml_first = GENERATE(false, true);
  mfem::Mesh serial_mesh =
      mfem::Mesh::MakeCartesian3D(2, 2, nz, mfem::Element::HEXAHEDRON, 1.0, 1.0, 2.0);
  std::vector<int> partitioning(serial_mesh.GetNE());
  for (int i = 0; i < serial_mesh.GetNE(); i++)
  {
    mfem::Vector center;
    serial_mesh.GetElementCenter(i, center);
    const bool is_pml = (center(2) > 1.5);
    serial_mesh.SetAttribute(i, is_pml ? 2 : 1);
    const int layer = static_cast<int>(center(2) / (2.0 / nz));
    const int pml_rank = pml_first ? 0 : np - 1;
    const int physical_rank = (layer * (np - 1)) / (nz - 2) + (pml_first ? 1 : 0);
    partitioning[i] = (np == 1) ? 0 : (is_pml ? pml_rank : physical_rank);
  }
  serial_mesh.SetAttributes();

  // Expected number of PML subdomain unknowns: the DOFs of the PML elements which are not
  // DOFs of physical elements.
  mfem::ND_FECollection fec(order, serial_mesh.Dimension());
  mfem::FiniteElementSpace serial_fespace(&serial_mesh, &fec);
  std::vector<char> marker(serial_fespace.GetVSize(), 0);
  mfem::Array<int> vdofs;
  for (int i = 0; i < serial_mesh.GetNE(); i++)
  {
    serial_fespace.GetElementVDofs(i, vdofs);
    for (const int vdof : vdofs)
    {
      marker[(vdof >= 0) ? vdof : -1 - vdof] |= (serial_mesh.GetAttribute(i) == 2) ? 1 : 2;
    }
  }
  const auto n_pml_expected = std::ranges::count(marker, 1);
  REQUIRE(n_pml_expected > 0);

  std::vector<std::unique_ptr<Mesh>> mesh;
  mesh.push_back(std::make_unique<Mesh>(comm, serial_mesh, partitioning.data()));

  IntegrationSettingsGuard settings_guard;
  config::SolverData solver;
  solver.order = order;
  solver.linear.mg_max_levels = 2;
  solver.linear.mg_coarsening = MultigridCoarsening::LINEAR;
  solver.linear.pml_subdomain_solver = PMLSubdomainSolver::DIRECT;
  solver.linear.pc_mat_real = false;
  solver.linear.pc_mat_shifted = 0;
  fem::DefaultIntegrationOrder::p_trial = solver.order;

  config::MaterialData vacuum, pml;
  vacuum.attributes = {1};
  pml.attributes = {2};
  config::PMLData pml_pml;
  pml_pml.autodetect_geometry = false;
  pml_pml.directions = {false, false, false, false, false, true};
  pml_pml.thickness = {0.0, 0.0, 0.0, 0.0, 0.0, 0.5};
  pml_pml.reference_frequency = 2.0;
  config::DomainData domains;
  domains.attributes = {1, 2};
  domains.materials = {vacuum, pml};
  domains.pml = {pml_pml};
  domains.pml[0].attributes = pml.attributes;
  config::BoundaryData boundaries;
  Units units(1.0, 1.0);
  SpaceOperator space_op(solver, domains, boundaries, ProblemType::DRIVEN, units, mesh);

  constexpr double omega = 2.0;
  auto P = space_op.GetPreconditionerMatrix<ComplexOperator>(1.0 + 0.0i, 1i * omega,
                                                             -omega * omega + 0.0i, omega);
  const auto *mg = dynamic_cast<const ComplexMultigridOperator *>(P.get());
  REQUIRE(mg);
  REQUIRE(mg->GetNumLevels() > 1);
  const auto &tdofs = mg->GetPMLTrueDofs();
  CHECK(std::ranges::is_sorted(tdofs));
  CHECK(std::ranges::adjacent_find(tdofs) == tdofs.end());
  long long n_pml = static_cast<long long>(tdofs.size());
  Mpi::GlobalSum(1, &n_pml, comm);
  CHECK(n_pml == n_pml_expected);
}

TEST_CASE("SpaceOperator rejects PML regions in 2D simulations",
          "[spaceoperator][pml][Serial][Parallel]")
{
  // Rectangle [0, 2] x [0, 1] with a PML layer x > 1.5 (attribute 2).
  MPI_Comm comm = Mpi::World();
  mfem::Mesh serial_mesh =
      mfem::Mesh::MakeCartesian2D(8, 4, mfem::Element::QUADRILATERAL, true, 2.0, 1.0);
  for (int i = 0; i < serial_mesh.GetNE(); i++)
  {
    mfem::Vector center;
    serial_mesh.GetElementCenter(i, center);
    serial_mesh.SetAttribute(i, (center(0) > 1.5) ? 2 : 1);
  }
  serial_mesh.SetAttributes();
  std::vector<std::unique_ptr<Mesh>> mesh;
  mesh.push_back(std::make_unique<Mesh>(comm, serial_mesh));

  IntegrationSettingsGuard settings_guard;
  config::SolverData solver;
  solver.order = 1;
  solver.linear.mg_max_levels = 1;
  fem::DefaultIntegrationOrder::p_trial = solver.order;

  config::MaterialData vacuum, pml;
  vacuum.attributes = {1};
  pml.attributes = {2};
  config::PMLData pml_pml;
  pml_pml.reference_frequency = 2.0;
  config::DomainData domains;
  domains.attributes = {1, 2};
  domains.materials = {vacuum, pml};
  domains.pml = {pml_pml};
  domains.pml[0].attributes = pml.attributes;
  config::BoundaryData boundaries;
  Units units(1.0, 1.0);
  for (auto problem_type : {ProblemType::DRIVEN, ProblemType::EIGENMODE})
  {
    CHECK_THROWS_WITH(SpaceOperator(solver, domains, boundaries, problem_type, units, mesh),
                      Catch::Matchers::ContainsSubstring(
                          "PML regions are only supported for 3D simulations"));
  }

  // Other simulation types treat the PML regions as regular materials.
  SpaceOperator transient_op(solver, domains, boundaries, ProblemType::TRANSIENT, units,
                             mesh);
  CHECK(!transient_op.GetMaterialOp().HasPML());
}

TEST_CASE("SpaceOperator PML with unit stretch reproduces the bulk operators",
          "[spaceoperator][pml][Serial][Parallel][GPU]")
{
  using namespace std::complex_literals;

  // With σ = 0, κ = 1, the PML tensors equal the background material tensors, so all PML
  // terms (including the Floquet terms) must reproduce the standard bulk integrators.
  MPI_Comm comm = Mpi::World();
  mfem::Mesh serial_mesh = mfem::Mesh::MakeCartesian3D(3, 3, 3, mfem::Element::HEXAHEDRON);
  for (int i = 0; i < serial_mesh.GetNE(); i++)
  {
    mfem::Vector center;
    serial_mesh.GetElementCenter(i, center);
    serial_mesh.SetAttribute(i, (center(0) > 2.0 / 3.0 || center(2) > 2.0 / 3.0) ? 2 : 1);
  }
  serial_mesh.SetAttributes();
  std::vector<std::unique_ptr<Mesh>> mesh;
  mesh.push_back(std::make_unique<Mesh>(comm, serial_mesh));

  IntegrationSettingsGuard settings_guard;
  config::SolverData solver;
  solver.order = 2;
  solver.pa_order_threshold = 2;
  solver.linear.mg_max_levels = 2;
  solver.linear.mg_coarsening = MultigridCoarsening::LINEAR;
  solver.linear.pc_mat_shifted = 0;
  BilinearForm::pa_order_threshold = solver.pa_order_threshold;
  fem::DefaultIntegrationOrder::p_trial = solver.order;
  fem::DefaultIntegrationOrder::q_order_jac = solver.q_order_jac;
  fem::DefaultIntegrationOrder::q_order_extra_pk = solver.q_order_extra;
  fem::DefaultIntegrationOrder::q_order_extra_qk = solver.q_order_extra;

  const bool isotropic = GENERATE(false, true);
  const bool scaled = GENERATE(false, true);
  auto MakeSpaceOperator = [&](bool with_pml)
  {
    config::MaterialData vacuum, layer;
    vacuum.attributes = {1};
    layer.attributes = {2};
    if (isotropic)
    {
      layer.mu_r.s = {1.5, 1.5, 1.5};
      layer.epsilon_r.s = {3.0, 3.0, 3.0};
      layer.tandelta.s = {2.0e-2, 2.0e-2, 2.0e-2};
    }
    else
    {
      // Rotated material axes, so the background tensors are full.
      const double c = std::cos(0.3), s = std::sin(0.3);
      const std::array<std::array<double, 3>, 3> axes = {
          {{c, s, 0.0}, {-s * c, c * c, s}, {s * s, -c * s, c}}};
      layer.mu_r.s = {1.0, 2.0, 3.0};
      layer.mu_r.v = axes;
      layer.epsilon_r.s = {2.0, 4.0, 5.0};
      layer.epsilon_r.v = axes;
      layer.tandelta.s = {1.0e-2, 2.0e-2, 3.0e-2};
      layer.tandelta.v = axes;
    }
    config::DomainData domains;
    domains.attributes = {1, 2};
    domains.materials = {vacuum, layer};
    if (with_pml)
    {
      auto &pml = domains.pml.emplace_back();
      pml.attributes = layer.attributes;
      pml.autodetect_geometry = false;
      pml.directions = {false, true, false, false, false, true};
      pml.thickness = {0.0, 1.0 / 3.0, 0.0, 0.0, 0.0, 1.0 / 3.0};
      pml.sigma_max = {0.0, 0.0, 0.0};
      pml.reference_frequency = 2.0;
    }
    config::BoundaryData boundaries;
    boundaries.periodic.wave_vector = {0.3, -0.2, 0.5};
    boundaries.periodic.floquet_reference_freq = scaled ? 1.5 : 0.0;
    Units units(1.0, 1.0);
    return std::make_unique<SpaceOperator>(
        solver, domains, boundaries, scaled ? ProblemType::DRIVEN : ProblemType::EIGENMODE,
        units, mesh);
  };
  auto pml_op = MakeSpaceOperator(true), ref_op = MakeSpaceOperator(false);
  REQUIRE(pml_op->GetMaterialOp().HasPML());
  REQUIRE(pml_op->GetMaterialOp().HasWaveVector());
  REQUIRE(pml_op->GetMaterialOp().HasFloquetFrequencyScaling() == scaled);
  REQUIRE(!ref_op->GetMaterialOp().HasPML());

  auto CheckEqual = [&](const auto *A, const auto *B, auto &x, auto &y1, auto &y2)
  {
    REQUIRE(A);
    REQUIRE(B);
    for (int seed : {42, 314159})
    {
      if constexpr (std::is_same_v<std::decay_t<decltype(x)>, ComplexVector>)
      {
        x.Real().Randomize(seed + Mpi::Rank(comm));
        x.Imag().Randomize(2 * seed + 1 + Mpi::Rank(comm));
      }
      else
      {
        x.Randomize(seed + Mpi::Rank(comm));
      }
      A->Mult(x, y1);
      B->Mult(x, y2);
      const double norm = linalg::Norml2(comm, y2);
      REQUIRE(norm > 0.0);
      y1.Add(-1.0, y2);
      CHECK_THAT(linalg::Norml2(comm, y1) / norm, Catch::Matchers::WithinAbs(0.0, 1.0e-13));
    }
  };
  auto CheckComplex = [&](const ComplexOperator *A, const ComplexOperator *B)
  {
    ComplexVector x(A->Width()), y1(A->Height()), y2(A->Height());
    x.UseDevice(true);
    y1.UseDevice(true);
    y2.UseDevice(true);
    CheckEqual(A, B, x, y1, y2);
  };

  SECTION("System matrices")
  {
    constexpr auto diag = Operator::DIAG_ONE;
    CheckComplex(pml_op->GetStiffnessMatrix<ComplexOperator>(diag).get(),
                 ref_op->GetStiffnessMatrix<ComplexOperator>(diag).get());
    CheckComplex(pml_op->GetMassMatrix<ComplexOperator>(Operator::DIAG_ZERO).get(),
                 ref_op->GetMassMatrix<ComplexOperator>(Operator::DIAG_ZERO).get());
    if (scaled)
    {
      CheckComplex(pml_op->GetDampingMatrix<ComplexOperator>(Operator::DIAG_ZERO).get(),
                   ref_op->GetDampingMatrix<ComplexOperator>(Operator::DIAG_ZERO).get());
    }
    CHECK(
        !pml_op->GetExtraSystemMatrix(std::complex<double>(2.0, 0.1), Operator::DIAG_ZERO));
  }

  SECTION("Preconditioner")
  {
    // The real-valued approximation uses the entrywise magnitudes of the PML tensors, which
    // matches the bulk approximation for isotropic materials.
    for (bool pc_mat_real : {false, true})
    {
      if (pc_mat_real && !isotropic)
      {
        continue;
      }
      solver.linear.pc_mat_real = pc_mat_real;
      auto pml_pc_op = MakeSpaceOperator(true), ref_pc_op = MakeSpaceOperator(false);
      const std::complex<double> omega = {2.0, 0.1}, lambda = 1i * omega;
      auto P1 = pml_pc_op->GetPreconditionerMatrix<ComplexOperator>(1.0 + 0.0i, lambda,
                                                                    lambda * lambda, omega);
      auto P2 = ref_pc_op->GetPreconditionerMatrix<ComplexOperator>(1.0 + 0.0i, lambda,
                                                                    lambda * lambda, omega);
      const auto *mg1 = dynamic_cast<const ComplexMultigridOperator *>(P1.get());
      const auto *mg2 = dynamic_cast<const ComplexMultigridOperator *>(P2.get());
      REQUIRE(mg1);
      REQUIRE(mg2);
      REQUIRE(mg1->GetNumLevels() == 2);
      for (std::size_t l = 0; l < mg1->GetNumLevels(); l++)
      {
        CheckComplex(&mg1->GetOperatorAtLevel(l), &mg2->GetOperatorAtLevel(l));
        CheckComplex(&mg1->GetAuxiliaryOperatorAtLevel(l),
                     &mg2->GetAuxiliaryOperatorAtLevel(l));
      }
    }
  }
}

TEST_CASE("SpaceOperator PML operators match an MFEM coefficient reference",
          "[spaceoperator][pml][Serial][Parallel][GPU]")
{
  using namespace std::complex_literals;

  // Unit cube with a PML layer on the +x and +z faces (including the edge region), on
  // tetrahedral and hexahedral meshes. The PML curl-curl and mass operators are compared
  // to MFEM integrators with the PML tensors evaluated at the quadrature points by an
  // independent implementation of the stretch.
  const auto element = GENERATE(mfem::Element::TETRAHEDRON, mfem::Element::HEXAHEDRON);
  MPI_Comm comm = Mpi::World();
  auto serial_mesh = mfem::Mesh::MakeCartesian3D(3, 3, 3, element);
  for (int i = 0; i < serial_mesh.GetNE(); i++)
  {
    mfem::Vector center;
    serial_mesh.GetElementCenter(i, center);
    serial_mesh.SetAttribute(i, (center(0) > 2.0 / 3.0 || center(2) > 2.0 / 3.0) ? 2 : 1);
  }
  serial_mesh.SetAttributes();
  std::vector<std::unique_ptr<Mesh>> mesh;
  mesh.push_back(std::make_unique<Mesh>(comm, serial_mesh));

  IntegrationSettingsGuard settings_guard;
  config::SolverData solver;
  solver.order = 2;
  solver.linear.mg_max_levels = 1;
  fem::DefaultIntegrationOrder::p_trial = solver.order;

  config::MaterialData vacuum, layer;
  vacuum.attributes = {1};
  layer.attributes = {2};
  layer.mu_r.s = {1.0, 2.0, 1.5};
  layer.epsilon_r.s = {2.0, 3.0, 4.0};
  layer.tandelta.s = {1.0e-2, 2.0e-2, 3.0e-2};
  config::PMLData layer_pml;
  layer_pml.autodetect_geometry = false;
  layer_pml.directions = {false, true, false, false, false, true};
  layer_pml.thickness = {0.0, 1.0 / 3.0, 0.0, 0.0, 0.0, 1.0 / 3.0};
  layer_pml.sigma_max = {3.0, 0.0, 5.0};
  layer_pml.kappa_max = {1.5, 1.0, 2.0};
  layer_pml.alpha_max = {0.2, 0.0, 0.1};
  layer_pml.reference_frequency = 2.0;
  config::DomainData domains;
  domains.attributes = {1, 2};
  domains.materials = {vacuum, layer};
  domains.pml = {layer_pml};
  domains.pml[0].attributes = layer.attributes;
  config::BoundaryData boundaries;
  Units units(1.0, 1.0);
  SpaceOperator space_op(solver, domains, boundaries, ProblemType::DRIVEN, units, mesh);
  REQUIRE(space_op.GetMaterialOp().HasPML());
  const auto &p = space_op.GetMaterialOp().GetPMLLayers().at(0).GetStretch();

  // Diagonal background material of the layer.
  std::array<double, 3> b_mu_inv, b_eps_re, b_eps_im;
  for (int i = 0; i < 3; i++)
  {
    b_mu_inv[i] = 1.0 / layer.mu_r.s[i];
    b_eps_re[i] = layer.epsilon_r.s[i];
    b_eps_im[i] = -layer.epsilon_r.s[i] * layer.tandelta.s[i];
  }

  // Reference stretch factors and tensors at a physical point.
  auto Stretch = [&p](const mfem::Vector &x)
  {
    std::array<std::complex<double>, 3> s;
    const std::array<double, 3> xp = {x(0), x(1), x(2)};
    for (int a = 0; a < 3; a++)
    {
      const auto &g = p.geometry;
      double r = 0.0;
      if (g.thickness[2 * a] > 0.0 && xp[a] < g.inner[2 * a])
      {
        r = (g.inner[2 * a] - xp[a]) / g.thickness[2 * a];
      }
      else if (g.thickness[2 * a + 1] > 0.0 && xp[a] > g.inner[2 * a + 1])
      {
        r = (xp[a] - g.inner[2 * a + 1]) / g.thickness[2 * a + 1];
      }
      const double shape = std::pow(std::min(r, 1.0), p.order);
      const bool pos = (xp[a] > p.geometry.inner[2 * a + 1]);
      const double sigma = shape * p.sigma_max[2 * a + (pos ? 1 : 0)];
      s[a] = 1.0 + (p.kappa_max[a] - 1.0) * shape +
             sigma / (p.alpha_max[a] * shape + 1i * p.reference_frequency);
    }
    return s;
  };
  auto MakeCoefficient = [&](bool muinv, bool imag)
  {
    return mfem::MatrixFunctionCoefficient(
        3,
        [=](const mfem::Vector &x, mfem::DenseMatrix &T)
        {
          // Bulk material values are evaluated by the PWConstCoefficient below.
          const auto s = Stretch(x);
          const auto det = s[0] * s[1] * s[2];
          T = 0.0;
          for (int i = 0; i < 3; i++)
          {
            const std::complex<double> b =
                muinv ? b_mu_inv[i] : b_eps_re[i] + 1i * b_eps_im[i];
            const std::complex<double> t =
                muinv ? b * s[i] * s[i] / det : b * det / (s[i] * s[i]);
            T(i, i) = imag ? t.imag() : t.real();
          }
        });
  };

  const auto &pfes = space_op.GetNDSpace().Get();
  const int q_order = fem::DefaultIntegrationOrder::Get(
      mesh[0]->Get(), mesh[0]->Get().GetElementGeometry(0));
  auto AssembleReference = [&](bool muinv, bool imag)
  {
    mfem::ParBilinearForm a(const_cast<mfem::ParFiniteElementSpace *>(&pfes));
    mfem::Array<int> pml_marker({0, 1}), bulk_marker({1, 0});
    auto *pml_coef = new mfem::MatrixFunctionCoefficient(MakeCoefficient(muinv, imag));
    auto *bulk_coef = new mfem::ConstantCoefficient(muinv ? 1.0 : 1.0);
    const auto &ir = mfem::IntRules.Get(mesh[0]->Get().GetElementGeometry(0), q_order);
    mfem::BilinearFormIntegrator *integ_pml, *integ_bulk = nullptr;
    if (muinv)
    {
      integ_pml = new mfem::CurlCurlIntegrator(*pml_coef);
      if (!imag)
      {
        integ_bulk = new mfem::CurlCurlIntegrator(*bulk_coef);
      }
    }
    else
    {
      integ_pml = new mfem::VectorFEMassIntegrator(*pml_coef);
      if (!imag)
      {
        integ_bulk = new mfem::VectorFEMassIntegrator(*bulk_coef);
      }
    }
    integ_pml->SetIntRule(&ir);
    a.AddDomainIntegrator(integ_pml, pml_marker);
    if (integ_bulk)
    {
      integ_bulk->SetIntRule(&ir);
      a.AddDomainIntegrator(integ_bulk, bulk_marker);
    }
    a.Assemble();
    a.Finalize();
    std::unique_ptr<mfem::HypreParMatrix> A(a.ParallelAssemble());
    delete pml_coef;
    delete bulk_coef;
    return A;
  };

  for (bool muinv : {true, false})
  {
    auto Ar = AssembleReference(muinv, false), Ai = AssembleReference(muinv, true);
    auto A = muinv ? space_op.GetStiffnessMatrix<ComplexOperator>(Operator::DIAG_ZERO)
                   : space_op.GetMassMatrix<ComplexOperator>(Operator::DIAG_ZERO);
    REQUIRE(A);
    ComplexVector x(A->Width()), y(A->Height()), y_ref(A->Height());
    x.Real().Randomize(42 + Mpi::Rank(comm));
    x.Imag().Randomize(4242 + Mpi::Rank(comm));
    A->Mult(x, y);
    Vector t(A->Height());
    Ar->Mult(x.Real(), y_ref.Real());
    Ai->Mult(x.Imag(), t);
    y_ref.Real().Add(-1.0, t);
    Ar->Mult(x.Imag(), y_ref.Imag());
    Ai->Mult(x.Real(), t);
    y_ref.Imag().Add(1.0, t);
    const double norm = linalg::Norml2(comm, y_ref);
    y.Add(-1.0, y_ref);
    INFO("muinv = " << muinv << ", element = " << static_cast<int>(element));
    CHECK_THAT(linalg::Norml2(comm, y) / norm, Catch::Matchers::WithinAbs(0.0, 1.0e-12));
  }
}

namespace
{

// Unit cube of tetrahedra with a lossy dielectric (attribute 1, z < 1/3) and vacuum
// (attribute 3), under a Cartesian PML layer on +z (attribute 2, z > 2/3) with the lossy
// dielectric as background, static or frequency-dependent with α > 0.
std::unique_ptr<SpaceOperator>
MakePMLSpaceOperator(std::vector<std::unique_ptr<Mesh>> &mesh, bool frequency_dependent,
                     int order = 2, int mg_max_levels = 1)
{
  auto serial_mesh = mfem::Mesh::MakeCartesian3D(3, 3, 3, mfem::Element::TETRAHEDRON);
  for (int i = 0; i < serial_mesh.GetNE(); i++)
  {
    mfem::Vector center;
    serial_mesh.GetElementCenter(i, center);
    serial_mesh.SetAttribute(
        i, (center(2) > 2.0 / 3.0) ? 2 : ((center(2) < 1.0 / 3.0) ? 1 : 3));
  }
  serial_mesh.SetAttributes();
  config::PMLData pml;
  pml.attributes = {2};
  pml.autodetect_geometry = false;
  pml.directions = {false, false, false, false, false, true};
  pml.thickness = {0.0, 0.0, 0.0, 0.0, 0.0, 1.0 / 3.0};
  pml.kappa_max.fill(1.5);
  pml.frequency_dependent = frequency_dependent;
  pml.alpha_max.fill(frequency_dependent ? 0.3 : 0.0);
  pml.reference_frequency = 2.0;
  mesh.clear();
  mesh.push_back(std::make_unique<Mesh>(Mpi::World(), serial_mesh));

  config::SolverData solver;
  solver.order = order;
  solver.linear.mg_max_levels = mg_max_levels;
  solver.linear.mg_coarsening = MultigridCoarsening::LINEAR;
  fem::DefaultIntegrationOrder::p_trial = solver.order;
  config::MaterialData dielectric, vacuum;
  dielectric.attributes = {1, 2};
  dielectric.mu_r.s = {1.0, 2.0, 1.5};
  dielectric.epsilon_r.s = {2.0, 3.0, 4.0};
  dielectric.tandelta.s = {1.0e-2, 2.0e-2, 3.0e-2};
  vacuum.attributes = {3};
  config::DomainData domains;
  domains.attributes = {1, 2, 3};
  domains.materials = {dielectric, vacuum};
  domains.pml = {pml};
  config::BoundaryData boundaries;
  Units units(1.0, 1.0);
  auto space_op = std::make_unique<SpaceOperator>(solver, domains, boundaries,
                                                  ProblemType::DRIVEN, units, mesh);
  REQUIRE(space_op->GetMaterialOp().GetPMLLayers().size() == 1);
  return space_op;
}

// Relative difference over all processes of local vectors, which may be empty.
double RelativeError(MPI_Comm comm, const Vector &y, const Vector &y_ref)
{
  double norm[2] = {0.0, 0.0};
  if (y_ref.Size() > 0)
  {
    Vector d(y);
    d -= y_ref;
    norm[0] = std::pow(d.Norml2(), 2);
    norm[1] = std::pow(y_ref.Norml2(), 2);
  }
  Mpi::GlobalSum(2, norm, comm);
  REQUIRE(norm[1] > 0.0);
  return std::sqrt(norm[0] / norm[1]);
}

// The PML integrator T with its default stored quadrature data gives the action of the
// evaluation of the tensors at every application, on the local elements of the layer.
template <typename T>
void CheckPMLQuadratureData(const FiniteElementSpace &space, const pml::Layer &layer,
                            const pml::ContextHeader &header)
{
  const auto ctx = layer.PackContext(header);
  Vector y, y_ref;
  if (!ctx.empty())
  {
    auto Assemble = [&](bool assemble_q_data)
    {
      BilinearForm a(space);
      a.AddDomainIntegratorOnAttributes<T>(layer.GetAttributes(), ctx.data(),
                                           ctx.size() * sizeof(CeedIntScalar),
                                           assemble_q_data);
      return a.PartialAssemble();
    };
    const auto A = Assemble(true), A_ref = Assemble(false);
    Vector x(A->Width());
    y.SetSize(A->Height());
    y_ref.SetSize(A->Height());
    x.UseDevice(true);
    y.UseDevice(true);
    y_ref.UseDevice(true);
    x.Randomize(42);
    A->Mult(x, y);
    A_ref->Mult(x, y_ref);
  }
  CHECK_THAT(RelativeError(space.GetComm(), y, y_ref),
             Catch::Matchers::WithinAbs(0.0, 1.0e-12));
}

}  // namespace

TEST_CASE("SpaceOperator PML quadrature data matches the evaluation of the tensors",
          "[spaceoperator][pml][Serial][Parallel][GPU]")
{
  // Every PML term and output part with the stored tensors, against their evaluation at
  // every application, for a static stretch and a frequency-dependent one with α > 0.
  const bool frequency_dependent = GENERATE(false, true);
  IntegrationSettingsGuard settings_guard;
  std::vector<std::unique_ptr<Mesh>> mesh;
  auto space_op = MakePMLSpaceOperator(mesh, frequency_dependent);
  const auto &layer = space_op->GetMaterialOp().GetPMLLayers()[0];
  const auto &nd = space_op->GetNDSpace(), &h1 = space_op->GetH1Space();
  pml::ContextHeader header;
  header.c_muinv = {0.7, -0.4};
  header.c_eps = {1.3, 0.2};
  header.omega = 1.1;
  const std::array<double, 3> k = {0.3, -0.2, 0.5};
  header.wave_vector_cross = {0.0, k[2], -k[1], -k[2], 0.0, k[0], k[1], -k[0], 0.0};
  for (auto part :
       {pml::TensorPart::REAL, pml::TensorPart::IMAG, pml::TensorPart::ABS,
        pml::TensorPart::STRETCH_MAGNITUDE_REAL, pml::TensorPart::STRETCH_MAGNITUDE_IMAG})
  {
    header.part = part;
    INFO("part = " << static_cast<int>(part));
    CheckPMLQuadratureData<CurlCurlPMLIntegrator>(nd, layer, header);
    CheckPMLQuadratureData<VectorFEMassPMLIntegrator>(nd, layer, header);
    CheckPMLQuadratureData<CurlCurlMassPMLIntegrator>(nd, layer, header);
    CheckPMLQuadratureData<FloquetMassPMLIntegrator>(nd, layer, header);
    CheckPMLQuadratureData<FloquetCrossPMLIntegrator>(nd, layer, header);
    CheckPMLQuadratureData<DiffusionPMLIntegrator>(h1, layer, header);
    CheckPMLQuadratureData<FloquetDiffusionPMLIntegrator>(h1, layer, header);
  }
}

// With a frequency dependent Floquet wave vector k = (ω / ω_ref) k_ref
// ("FloquetReferenceFrequency"), the Floquet terms of the domain and of the second-order
// absorbing boundary are split by their power of ω between the stiffness, damping, mass,
// and extra system matrices. At ω = ω_ref, the system matrix must be that of the fixed wave
// vector k_ref, whose boundary Floquet terms are in the extra system matrix i·(0.5/ω)·M_ff.
// The comparison holds for any vector, also where the periodic field is not uniform on the
// boundary, which the mixed boundary curl terms need to contribute.
TEST_CASE_METHOD(palace::test::SharedTempDir,
                 "SpaceOperator Floquet terms of the second-order absorbing boundary",
                 "[spaceoperator][Serial][Parallel]")
{
  using namespace std::complex_literals;
  using json = nlohmann::json;
  MPI_Comm comm = Mpi::World();

  // Box of the floquet_absorbing regression case, periodic in x and y, with a lumped port
  // sheet at z = 0 (attribute 5, the excitation, with a frequency independent term) and the
  // absorbing boundary at z = L (attribute 6). The wave vector (1/cm) has two tangential
  // components and is in the first Brillouin zone, so that it is not wrapped for the fixed
  // wave vector.
  json config;
  config["Problem"] = {{"Type", "Driven"}, {"Verbose", 0}, {"Output", temp_dir}};
  config["Model"] = {
      {"Mesh", fs::path(PALACE_TEST_DATA_DIR) /
                   "regression/input/floquet_absorbing/mesh/floquet_absorbing.msh"},
      {"L0", 1.0e-2}};
  config["Domains"] = {
      {"Materials", json::array({json::object({{"Attributes", json::array({1})},
                                               {"Permeability", 2.0},
                                               {"Permittivity", 3.0}})})}};
  config["Boundaries"] = {
      {"Periodic",
       {{"FloquetWaveVector", json::array({0.3, 0.4, 0.0})},
        {"BoundaryPairs",
         json::array({json::object({{"DonorAttributes", json::array({1})},
                                    {"ReceiverAttributes", json::array({2})}}),
                      json::object({{"DonorAttributes", json::array({3})},
                                    {"ReceiverAttributes", json::array({4})}})})}}},
      {"LumpedPort", json::array({json::object({{"Index", 1},
                                                {"Attributes", json::array({5})},
                                                {"Excitation", uint(1)},
                                                {"R", 50.0},
                                                {"Direction", "+X"}})})},
      {"Absorbing", {{"Attributes", json::array({6})}, {"Order", 2}}}};
  config["Solver"] = {{"Order", 2UL},
                      {"Device", "CPU"},
                      {"Driven", {{"MinFreq", 2.0}, {"MaxFreq", 4.0}, {"FreqStep", 1.0}}}};
  auto config_scaled = config;
  config_scaled["Boundaries"]["Periodic"]["FloquetReferenceFrequency"] = 3.0;

  IoData iodata(config, false), iodata_scaled(config_scaled, false);
  auto mesh = LoadMesh(iodata, comm);
  auto mesh_scaled = LoadMesh(iodata_scaled, comm);
  SpaceOperator space_op(iodata, mesh), space_op_scaled(iodata_scaled, mesh_scaled);
  REQUIRE(space_op.GetMaterialOp().HasWaveVector());
  REQUIRE_FALSE(space_op.GetMaterialOp().HasFloquetFrequencyScaling());
  REQUIRE(space_op_scaled.GetMaterialOp().HasFloquetFrequencyScaling());
  const double omega_ref = space_op_scaled.GetMaterialOp().GetFloquetOmegaRef();

  ComplexVector x(space_op.GetNDSpace().GetTrueVSize()), y(x.Size()), y_scaled(x.Size());
  REQUIRE(space_op_scaled.GetNDSpace().GetTrueVSize() == x.Size());
  x.UseDevice(true);
  y.UseDevice(true);
  y_scaled.UseDevice(true);
  x.Real().Randomize(1 + Mpi::Rank(comm));
  x.Imag().Randomize(7 + Mpi::Rank(comm));

  // A = K + iω C - ω² M + A2(ω), as assembled by the driven solver.
  auto ApplySystemMatrix = [&x](SpaceOperator &op, double omega, ComplexVector &y)
  {
    auto K = op.GetStiffnessMatrix<ComplexOperator>(Operator::DIAG_ONE);
    auto C = op.GetDampingMatrix<ComplexOperator>(Operator::DIAG_ZERO);
    auto M = op.GetMassMatrix<ComplexOperator>(Operator::DIAG_ZERO);
    auto A2 = op.GetExtraSystemOperator(omega, Operator::DIAG_ZERO);
    auto A = op.GetSystemMatrix(1.0 + 0.0i, 1i * omega, -omega * omega + 0.0i, K.get(),
                                C.get(), M.get(), A2.get());
    A->Mult(x, y);
  };

  SECTION("Frequency dependent wave vector at the reference frequency")
  {
    ApplySystemMatrix(space_op, omega_ref, y);
    ApplySystemMatrix(space_op_scaled, omega_ref, y_scaled);
    CHECK(RelativeDifference(comm, y_scaled, y) < 1.0e-12);

    // Away from the reference frequency, the wave vectors differ.
    ApplySystemMatrix(space_op, 1.5 * omega_ref, y);
    ApplySystemMatrix(space_op_scaled, 1.5 * omega_ref, y_scaled);
    CHECK(RelativeDifference(comm, y_scaled, y) > 1.0e-3);
  }

  SECTION("Frequency independent boundary matrix")
  {
    // The second-order absorbing boundary is the only frequency dependent term, so that
    // A2(ω) = i·(0.5/ω)·M_ff. With the fixed wave vector, M_ff includes the boundary
    // Floquet terms and is complex.
    for (auto *op : {&space_op, &space_op_scaled})
    {
      const double omega = 1.5 * omega_ref;
      auto A2 = op->GetExtraSystemOperator(omega, Operator::DIAG_ZERO);
      auto M_ff =
          op->GetFarfieldBoundaryCurlCurlMatrix<ComplexOperator>(Operator::DIAG_ZERO);
      REQUIRE(A2);
      REQUIRE(M_ff);
      A2->Mult(x, y);
      M_ff->Mult(x, y_scaled);
      y_scaled *= 1i * 0.5 / omega;
      CHECK(RelativeDifference(comm, y_scaled, y) < 1.0e-12);
    }
  }
}

// The driven solvers apply the system operator through the finest level of the
// preconditioner when that level is the sparse part of A = K + iω C - ω² (Mr + i Mi) +
// A2(ω) (see SpaceOperator::GetSystemOperator). The wave-port modal correction, which is
// the matrix-free part of A2, is covered in test-waveportoperator.cpp.
TEST_CASE("Driven system operator from the preconditioner",
          "[spaceoperator][Serial][Parallel]")
{
  using namespace std::complex_literals;

  const auto comm = Mpi::World();
  auto serial_mesh = mfem::Mesh::MakeCartesian3D(2, 2, 2, mfem::Element::TETRAHEDRON);
  // A material loss tangent and an electrical conductivity cannot be combined on the same
  // attribute, so split the domain: the first material gives the mass term an imaginary
  // part, the second gives a damping term, and both are needed for every part of A to be
  // nonzero.
  for (int i = 0; i < serial_mesh.GetNE(); i++)
  {
    serial_mesh.SetAttribute(i, 1 + (i % 2));
  }
  serial_mesh.SetAttributes();
  std::vector<std::unique_ptr<Mesh>> mesh;
  mesh.push_back(std::make_unique<Mesh>(comm, serial_mesh));

  IntegrationSettingsGuard settings_guard;
  config::SolverData solver;
  solver.order = 2;
  solver.pa_order_threshold = 2;
  solver.linear.mg_max_levels = 2;
  solver.linear.mg_coarsening = MultigridCoarsening::LINEAR;
  solver.linear.pc_mat_real = false;
  solver.linear.pc_mat_shifted = 0;
  BilinearForm::pa_order_threshold = solver.pa_order_threshold;
  fem::DefaultIntegrationOrder::p_trial = solver.order;
  fem::DefaultIntegrationOrder::q_order_jac = solver.q_order_jac;
  fem::DefaultIntegrationOrder::q_order_extra_pk = solver.q_order_extra;
  fem::DefaultIntegrationOrder::q_order_extra_qk = solver.q_order_extra;

  config::MaterialData lossy;
  lossy.attributes = {1};
  lossy.mu_r.s = {1.0, 2.0, 3.0};
  lossy.epsilon_r.s = {2.0, 4.0, 5.0};
  lossy.tandelta.s = {0.1, 0.1, 0.1};
  config::MaterialData conductive;
  conductive.attributes = {2};
  conductive.epsilon_r.s = {3.0, 3.0, 3.0};
  conductive.sigma.s = {0.5, 0.5, 0.5};
  config::DomainData domains;
  domains.attributes = {1, 2};
  domains.materials = {lossy, conductive};

  // Essential dofs, so the diagonal policies of the two operators are exercised, and a
  // surface impedance boundary contributing to the boundary stiffness, damping and mass.
  // It has two attributes: repeated terms on one boundary group must update only the listed
  // attributes even when equal material properties initially share coefficient storage.
  config::ImpedanceData impedance;
  impedance.Rs = 2.0;
  impedance.Ls = 0.5;
  impedance.Cs = 0.25;
  impedance.attributes = {4, 5};
  config::BoundaryData boundaries;
  boundaries.pec.attributes = {1, 2, 3};
  boundaries.impedance = {impedance};
  Units units(1.0, 1.0);

  const double omega = 2.0;
  auto GetPreconditioner = [&](SpaceOperator &space_op)
  {
    return space_op.GetPreconditionerMatrix<ComplexOperator>(1.0 + 0.0i, 1i * omega,
                                                             -omega * omega + 0.0i, omega);
  };

  // The system operator from the preconditioner and its transpose must match the ones
  // assembled from K, C, M and GetExtraSystemOperator, as the driven solvers do otherwise.
  auto CheckSystemOperator = [&](const config::BoundaryData &test_boundaries)
  {
    SpaceOperator space_op(solver, domains, test_boundaries, ProblemType::DRIVEN, units,
                           mesh);
    auto K = space_op.GetStiffnessMatrix<ComplexOperator>(Operator::DIAG_ONE);
    auto C = space_op.GetDampingMatrix<ComplexOperator>(Operator::DIAG_ZERO);
    auto M = space_op.GetMassMatrix<ComplexOperator>(Operator::DIAG_ZERO);
    auto A2 = space_op.GetExtraSystemOperator(omega, Operator::DIAG_ZERO);
    auto A_ref = space_op.GetSystemMatrix(1.0 + 0.0i, 1i * omega, -omega * omega + 0.0i,
                                          K.get(), C.get(), M.get(), A2.get());
    auto P = GetPreconditioner(space_op);
    auto A = space_op.GetSystemOperator(omega, *P);
    REQUIRE(A);

    ComplexVector x(A->Width()), y_ref(A->Height()), y(A->Height());
    x.UseDevice(true);
    y_ref.UseDevice(true);
    y.UseDevice(true);
    linalg::SetRandom(comm, x);
    for (bool transpose : {false, true})
    {
      CAPTURE(transpose);
      if (transpose)
      {
        A_ref->MultTranspose(x, y_ref);
        A->MultTranspose(x, y);
      }
      else
      {
        A_ref->Mult(x, y_ref);
        A->Mult(x, y);
      }
      y -= y_ref;
      CHECK(linalg::Norml2(comm, y) <= 1.0e-12 * linalg::Norml2(comm, y_ref));
    }
  };

  // The preconditioner is another matrix, so there is no system operator from it.
  auto CheckNoSystemOperator = [&](const config::SolverData &test_solver,
                                   const config::BoundaryData &test_boundaries)
  {
    SpaceOperator space_op(test_solver, domains, test_boundaries, ProblemType::DRIVEN,
                           units, mesh);
    auto P = GetPreconditioner(space_op);
    CHECK(!space_op.GetSystemOperator(omega, *P));
  };

  auto farfield_boundaries = boundaries;
  farfield_boundaries.impedance = {};
  farfield_boundaries.farfield.order = 2;
  farfield_boundaries.farfield.attributes = {4, 5, 6};

  SECTION("Surface impedance boundary")
  {
    CheckSystemOperator(boundaries);
  }

  SECTION("Second-order farfield boundary")
  {
    CheckSystemOperator(farfield_boundaries);
  }

  SECTION("Surface conductivity boundary")
  {
    auto conductivity_boundaries = boundaries;
    auto &conductivity = conductivity_boundaries.conductivity.emplace_back();
    conductivity.sigma = 0.5;
    conductivity.attributes = {6};
    CheckSystemOperator(conductivity_boundaries);
  }

  SECTION("Rational impedance boundary")
  {
    // Series RL impedance Z(s) = 0.5 s + 2.
    auto rational_boundaries = boundaries;
    auto &rational = rational_boundaries.rational_impedance.emplace_back();
    rational.num = {0.5, 2.0};
    rational.den = {1.0};
    rational.attributes = {6};
    CheckSystemOperator(rational_boundaries);
  }

  SECTION("Real-valued or shifted preconditioner")
  {
    auto real_solver = solver;
    real_solver.linear.pc_mat_real = true;
    CheckNoSystemOperator(real_solver, boundaries);
    auto shifted_solver = solver;
    shifted_solver.linear.pc_mat_shifted = 1;
    CheckNoSystemOperator(shifted_solver, boundaries);
  }

  SECTION("AMS")
  {
    // AMS and BoomerAMG flip the sign of the second-order absorbing boundary term in the
    // preconditioner.
    auto ams_solver = solver;
    ams_solver.linear.type = LinearSolver::AMS;
    CheckNoSystemOperator(ams_solver, farfield_boundaries);
  }

  SECTION("Single multigrid level")
  {
    // Its only level is the assembled coarse matrix; keep the system matrix-free.
    auto single_solver = solver;
    single_solver.linear.mg_max_levels = 1;
    CheckNoSystemOperator(single_solver, boundaries);
  }

  SECTION("Floquet wave vector")
  {
    // The domain Floquet terms are in the stiffness matrix for a fixed wave vector, and in
    // the damping and mass matrices when it scales with frequency. The Floquet terms of the
    // second-order farfield boundary are in the extra system matrix.
    auto floquet_boundaries = boundaries;
    floquet_boundaries.periodic.wave_vector = {0.1, 0.2, 0.0};
    floquet_boundaries.periodic.floquet_reference_freq = GENERATE(0.0, 0.3);
    floquet_boundaries.farfield.order = 2;
    floquet_boundaries.farfield.attributes = {6};
    CheckSystemOperator(floquet_boundaries);
  }
}

}  // namespace palace
