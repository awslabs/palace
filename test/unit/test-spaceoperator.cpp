// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include <algorithm>
#include <complex>
#include <limits>
#include <vector>
#include <catch2/catch_test_macros.hpp>
#include <nlohmann/json.hpp>
#include <catch2/generators/catch_generators.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "fem/bilinearform.hpp"
#include "fem/integrator.hpp"
#include "fem/libceed/operator.hpp"
#include "fem/mesh.hpp"
#include "fixtures.hpp"
#include "linalg/hypre.hpp"
#include "linalg/rap.hpp"
#include "models/materialoperator.hpp"
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
    // Periodic wave-vector terms are outside the reuse contract.
    auto floquet_boundaries = boundaries;
    floquet_boundaries.periodic.wave_vector = {0.1, 0.2, 0.0};
    CheckNoSystemOperator(solver, floquet_boundaries);
  }
}

}  // namespace palace
