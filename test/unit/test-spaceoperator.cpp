// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include <algorithm>
#include <complex>
#include <limits>
#include <map>
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
#include "linalg/solver.hpp"
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

struct SpaceOperatorExposed : public SpaceOperator
{
  using SpaceOperator::pc_term_coarse;
  using SpaceOperator::pc_terms_built;
  using SpaceOperator::SpaceOperator;
};

class RecordingSolver : public mfem::Solver
{
public:
  int operator_height = 0;
  int operator_width = 0;

  void SetOperator(const mfem::Operator &op) override
  {
    operator_height = op.Height();
    operator_width = op.Width();
    height = operator_height;
    width = operator_width;
  }

  void Mult(const Vector &, Vector &y) const override { y = 0.0; }
};

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

// The local rows of an assembled parallel matrix as maps from the global column index to
// the value, which is independent of the order of the entries within a row and of the local
// numbering of the off-diagonal columns.
using AssembledRows = std::vector<std::map<HYPRE_BigInt, double>>;

AssembledRows GetAssembledRows(const mfem::HypreParMatrix &A)
{
  hypre_CSRMatrix *csr = hypre_MergeDiagAndOffd((hypre_ParCSRMatrix *)A);
  hypre_CSRMatrixMigrate(csr, HYPRE_MEMORY_HOST);
  const HYPRE_Int *I = hypre_CSRMatrixI(csr);
  const HYPRE_BigInt *J = hypre_CSRMatrixBigJ(csr);
  const double *data = hypre_CSRMatrixData(csr);
  AssembledRows rows(hypre_CSRMatrixNumRows(csr));
  for (std::size_t i = 0; i < rows.size(); i++)
  {
    for (HYPRE_Int j = I[i]; j < I[i + 1]; j++)
    {
      rows[i][J[j]] = data[j];
    }
  }
  hypre_CSRMatrixDestroy(csr);
  return rows;
}

// The sparsity pattern of the local rows, as the global column indices of each of them.
std::vector<std::vector<HYPRE_BigInt>> GetPattern(const AssembledRows &rows)
{
  std::vector<std::vector<HYPRE_BigInt>> pattern(rows.size());
  for (std::size_t i = 0; i < rows.size(); i++)
  {
    for (const auto &[column, value] : rows[i])
    {
      pattern[i].push_back(column);
    }
  }
  return pattern;
}

// The assembled parallel matrix of one part of a complex parallel operator.
const mfem::HypreParMatrix &GetAssembledPart(const ComplexParOperator &op, bool imag)
{
  const auto *part = dynamic_cast<const ParOperator *>(imag ? op.Imag() : op.Real());
  REQUIRE(part);
  return part->ParallelAssemble();
}

struct ImaginaryHierarchyData
{
  AssembledRows coarse_rows;
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

  // The coarsest level is either assembled here or combined from the cached
  // frequency-independent term matrices, so compare it after the parallel assembly.
  ImaginaryHierarchyData data;
  constexpr bool imag = true;
  data.coarse_rows = GetAssembledRows(GetAssembledPart(*coarse_op, imag));

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

  auto AllValues = [](const AssembledRows &rows)
  {
    std::vector<double> values;
    for (const auto &row : rows)
    {
      for (const auto &[column, value] : row)
      {
        values.push_back(value);
      }
    }
    return values;
  };
  const auto zero_values = AllValues(zero_data.coarse_rows);
  const auto tiny_values = AllValues(tiny_data.coarse_rows);
  CHECK(GetPattern(zero_data.coarse_rows) == GetPattern(tiny_data.coarse_rows));
  CHECK(std::ranges::all_of(zero_values, [](double value) { return value == 0.0; }));
  CHECK(std::ranges::any_of(tiny_values, [](double value) { return value != 0.0; }));
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

// The Chebyshev smoothers ask every preconditioner level above the coarsest for its
// diagonal at each frequency, which the level operators assemble with libCEED. The diagonal
// is linear in the material property coefficients, so it is instead combined from the
// cached diagonals of the fixed spatial terms (stiffness, damping, real and imaginary mass,
// and Floquet mass and curl) using the same scalings as the
// coefficients themselves (SpaceOperator::GetPreconditionerScalars). A configuration whose
// coefficients are not those terms scaled by real numbers keeps the operator assembly.
TEST_CASE("Complex preconditioner diagonal is combined from cached terms",
          "[spaceoperator][Serial]")
{
  using namespace std::complex_literals;

  const auto comm = Mpi::World();
  auto serial_mesh = mfem::Mesh::MakeCartesian3D(1, 1, 1, mfem::Element::TETRAHEDRON);
  // A material loss tangent and an electrical conductivity cannot be combined on the same
  // attribute, so split the domain: both the damping and the imaginary mass term are then
  // nonzero, on different attributes.
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

  // Essential dofs, so the eliminated diagonal entries are exercised, and a surface
  // impedance boundary contributing to the boundary stiffness, damping and mass terms.
  config::ImpedanceData impedance;
  impedance.Rs = 2.0;
  impedance.Ls = 0.5;
  impedance.Cs = 0.25;
  impedance.attributes = {4};
  config::BoundaryData boundaries;
  boundaries.pec.attributes = {1, 2, 3};
  boundaries.impedance = {impedance};
  Units units(1.0, 1.0);

  // A second-order farfield with no attributes contributes no matrix but selects the
  // independent per-frequency assembly used as the numerical reference below.
  auto reference_boundaries = boundaries;
  reference_boundaries.farfield.order = 2;

  constexpr double omega = 2.0;
  // Compare the diagonals returned by cached and per-frequency operators on every level
  // above the coarsest of both hierarchies. The type of a3 selects the same real- or
  // complex-frequency boundary stamping path as the solver.
  auto DiagonalDifference = [&](SpaceOperator &space_op, SpaceOperator &reference_op,
                                std::complex<double> lambda, auto a3)
  {
    auto pc = space_op.GetPreconditionerMatrix<ComplexOperator>(1.0 + 0.0i, lambda,
                                                                lambda * lambda, a3);
    auto pc_ref = reference_op.GetPreconditionerMatrix<ComplexOperator>(
        1.0 + 0.0i, lambda, lambda * lambda, a3);
    const auto *mg = dynamic_cast<const ComplexMultigridOperator *>(pc.get());
    const auto *mg_ref = dynamic_cast<const ComplexMultigridOperator *>(pc_ref.get());
    REQUIRE(mg);
    REQUIRE(mg_ref);
    REQUIRE(mg->GetNumLevels() == 2);
    REQUIRE(mg_ref->GetNumLevels() == mg->GetNumLevels());
    double max_rel_diff = 0.0;
    for (bool aux : {false, true})
    {
      const auto &dbc_tdof_lists =
          aux ? space_op.GetH1DbcTDofLists() : space_op.GetNDDbcTDofLists();
      for (std::size_t l = 1; l < mg->GetNumLevels(); l++)
      {
        const auto &level_op =
            aux ? mg->GetAuxiliaryOperatorAtLevel(l) : mg->GetOperatorAtLevel(l);
        const auto &level_ref =
            aux ? mg_ref->GetAuxiliaryOperatorAtLevel(l) : mg_ref->GetOperatorAtLevel(l);
        ComplexVector diag(level_op.Height()), diag_ref(level_ref.Height());
        diag.UseDevice(true);
        diag_ref.UseDevice(true);
        level_op.AssembleDiagonal(diag);
        level_ref.AssembleDiagonal(diag_ref);
        for (bool imag : {false, true})
        {
          const Vector &ref = imag ? diag_ref.Imag() : diag_ref.Real();
          Vector diff(imag ? diag.Imag() : diag.Real());
          diff -= ref;
          double norm[2] = {diff.Normlinf(), ref.Normlinf()};
          Mpi::GlobalMax(2, norm, comm);
          CHECK(norm[1] > 0.0);
          max_rel_diff = std::max(max_rel_diff, norm[0] / norm[1]);
        }

        // The essential dofs carry the real scalar selected by DIAG_ONE: 1 + 0i.
        const auto &dbc_tdof_list = dbc_tdof_lists[l];
        int num_dbc = dbc_tdof_list.Size();
        Mpi::GlobalSum(1, &num_dbc, comm);
        REQUIRE(num_dbc > 0);
        const auto *dbc = dbc_tdof_list.HostRead();
        const auto *diag_real = diag.Real().HostRead();
        const auto *diag_imag = diag.Imag().HostRead();
        bool eliminated = true;
        for (int i = 0; i < dbc_tdof_list.Size(); i++)
        {
          eliminated = eliminated && diag_real[dbc[i]] == 1.0 && diag_imag[dbc[i]] == 0.0;
        }
        CHECK(eliminated);
      }
    }
    return max_rel_diff;
  };

  SECTION("Cached term diagonals are reused")
  {
    SpaceOperatorExposed space_op(solver, domains, boundaries, ProblemType::DRIVEN, units,
                                  mesh);
    SpaceOperator reference_op(solver, domains, reference_boundaries, ProblemType::DRIVEN,
                               units, mesh);
    CHECK(!space_op.pc_terms_built);
    CHECK(DiagonalDifference(space_op, reference_op, 1i * omega, omega) <= 1.0e-13);
    REQUIRE(space_op.pc_terms_built);
    const auto *stiffness = space_op.pc_term_coarse.stiffness.get();
    REQUIRE(stiffness);
    CHECK(DiagonalDifference(space_op, reference_op, 1.5i * omega, 1.5 * omega) <= 1.0e-13);
    CHECK(space_op.pc_term_coarse.stiffness.get() == stiffness);
  }

  SECTION("Complex frequency")
  {
    const std::complex<double> lambda = 0.1 + 1i * omega;
    SpaceOperatorExposed space_op(solver, domains, boundaries, ProblemType::EIGENMODE,
                                  units, mesh);
    SpaceOperator reference_op(solver, domains, reference_boundaries,
                               ProblemType::EIGENMODE, units, mesh);
    CHECK(DiagonalDifference(space_op, reference_op, lambda, lambda / 1i) <= 1.0e-13);
    CHECK(space_op.pc_terms_built);
  }

  SECTION("Shifted preconditioner mass coefficient")
  {
    auto shifted_solver = solver;
    shifted_solver.linear.pc_mat_shifted = 1;
    SpaceOperatorExposed space_op(shifted_solver, domains, boundaries, ProblemType::DRIVEN,
                                  units, mesh);
    SpaceOperator reference_op(shifted_solver, domains, reference_boundaries,
                               ProblemType::DRIVEN, units, mesh);
    CHECK(DiagonalDifference(space_op, reference_op, 1i * omega, omega) <= 1.0e-13);
    CHECK(space_op.pc_terms_built);
  }

  for (bool frequency_scaled : {false, true})
  {
    DYNAMIC_SECTION("Floquet wave vector, frequency scaled = " << frequency_scaled)
    {
      auto floquet_boundaries = boundaries;
      floquet_boundaries.periodic.wave_vector = {0.1, 0.0, 0.0};
      floquet_boundaries.periodic.floquet_reference_freq = frequency_scaled ? 1.0 : 0.0;
      auto floquet_reference_boundaries = floquet_boundaries;
      floquet_reference_boundaries.farfield.order = 2;
      SpaceOperatorExposed space_op(solver, domains, floquet_boundaries,
                                    ProblemType::DRIVEN, units, mesh);
      SpaceOperator reference_op(solver, domains, floquet_reference_boundaries,
                                 ProblemType::DRIVEN, units, mesh);
      REQUIRE(space_op.GetMaterialOp().HasWaveVector());
      CHECK(DiagonalDifference(space_op, reference_op, 1i * omega, omega) <= 1.0e-13);
      CHECK(space_op.pc_terms_built);
    }
  }

  SECTION("Frequency-dependent boundary keeps the per-frequency diagonal")
  {
    auto farfield_boundaries = boundaries;
    farfield_boundaries.farfield.order = 2;
    farfield_boundaries.farfield.attributes = {5, 6};
    SpaceOperatorExposed space_op(solver, domains, farfield_boundaries, ProblemType::DRIVEN,
                                  units, mesh);
    auto pc = space_op.GetPreconditionerMatrix<ComplexOperator>(
        1.0 + 0.0i, 1i * omega, -omega * omega + 0.0i, omega);
    REQUIRE(pc);
    CHECK(!space_op.pc_terms_built);
  }
}

// The coarsest level of the complex preconditioner hierarchy is a sparse matrix which is
// linear in the frequency scalings of the fixed preconditioner terms (stiffness, damping,
// real and imaginary mass, and Floquet mass and curl), so it is combined from the
// cached parallel matrices of those terms instead of being assembled and RAP'd at every
// frequency (SpaceOperator::CombinePreconditionerTermMatrices). The combination must have
// the same sparsity pattern as the per-frequency assembly, including the exact zeros of the
// terms which the current frequency does not scale in, the same values up to rounding, and
// the same essential true dof elimination.
TEST_CASE("Complex preconditioner coarsest level is combined from cached terms",
          "[spaceoperator][Serial][Parallel]")
{
  using namespace std::complex_literals;

  const auto comm = Mpi::World();
  auto serial_mesh = mfem::Mesh::MakeCartesian3D(1, 1, 1, mfem::Element::TETRAHEDRON);
  // A material loss tangent and an electrical conductivity cannot be combined on the same
  // attribute, so split the domain: both the damping and the imaginary mass term are then
  // nonzero, on different attributes.
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

  // Essential dofs, so the eliminated rows are exercised, and a surface impedance boundary
  // contributing to the boundary stiffness, damping and mass terms.
  config::ImpedanceData impedance;
  impedance.Rs = 2.0;
  impedance.Ls = 0.5;
  impedance.Cs = 0.25;
  impedance.attributes = {4};
  config::BoundaryData boundaries;
  boundaries.pec.attributes = {1, 2, 3};
  boundaries.impedance = {impedance};
  Units units(1.0, 1.0);

  // The reference: a second-order farfield boundary on no attributes contributes to no
  // coefficient, but excludes the cached combination, so this operator assembles and RAP's
  // its coarsest level at every frequency, which is what the combination replaces.
  auto reference_boundaries = boundaries;
  reference_boundaries.farfield.order = 2;

  auto CoarseOperator = [](const ComplexOperator &pc)
  {
    const auto *mg = dynamic_cast<const ComplexMultigridOperator *>(&pc);
    REQUIRE(mg);
    REQUIRE(mg->GetNumLevels() == 2);
    const auto *coarse_op =
        dynamic_cast<const ComplexParOperator *>(&mg->GetOperatorAtLevel(0));
    REQUIRE(coarse_op);
    return coarse_op;
  };

  // Compare both parts of the combined coarsest level with the per-frequency assembly of
  // the reference, at the driven coefficients of one frequency.
  auto CheckCoarseLevel =
      [&](SpaceOperator &space_op, SpaceOperator &reference_op, double omega)
  {
    const std::complex<double> a1 = 1i * omega;
    auto pc =
        space_op.GetPreconditionerMatrix<ComplexOperator>(1.0 + 0.0i, a1, a1 * a1, omega);
    auto pc_ref = reference_op.GetPreconditionerMatrix<ComplexOperator>(1.0 + 0.0i, a1,
                                                                        a1 * a1, omega);
    const auto *coarse_op = CoarseOperator(*pc);
    const auto *coarse_ref_op = CoarseOperator(*pc_ref);

    const auto &dbc_tdof_list = space_op.GetNDDbcTDofLists()[0];
    int num_dbc = dbc_tdof_list.Size();
    Mpi::GlobalSum(1, &num_dbc, comm);
    REQUIRE(num_dbc > 0);
    double max_rel_diff = 0.0;
    for (bool imag : {false, true})
    {
      const auto &A = GetAssembledPart(*coarse_op, imag);
      const auto &A_ref = GetAssembledPart(*coarse_ref_op, imag);
      const auto rows = GetAssembledRows(A);
      const auto rows_ref = GetAssembledRows(A_ref);
      REQUIRE(rows.size() == rows_ref.size());
      CHECK(GetPattern(rows) == GetPattern(rows_ref));

      // The entries agree to rounding, relative to the largest entry of the row.
      for (std::size_t i = 0; i < rows.size(); i++)
      {
        double row_scale = 0.0, row_diff = 0.0;
        for (const auto &[column, value] : rows_ref[i])
        {
          row_scale = std::max(row_scale, std::abs(value));
          const auto it = rows[i].find(column);
          row_diff = std::max(row_diff,
                              std::abs((it != rows[i].end() ? it->second : 0.0) - value));
        }
        max_rel_diff =
            std::max(max_rel_diff, (row_scale > 0.0) ? row_diff / row_scale : row_diff);
      }

      // The essential dofs carry the eliminated values of the policy of the part exactly:
      // the real part uses DIAG_ONE and the imaginary part DIAG_ZERO.
      const auto *dbc = dbc_tdof_list.HostRead();
      bool eliminated = true;
      for (int i = 0; i < dbc_tdof_list.Size(); i++)
      {
        for (const auto &[column, value] : rows[dbc[i]])
        {
          const bool diagonal = (column == A.RowPart()[0] + dbc[i]);
          eliminated = eliminated && value == (diagonal && !imag ? 1.0 : 0.0);
        }
      }
      CHECK(eliminated);

      // The first entry in each row of the diagonal block is the diagonal one, as the
      // per-frequency assembly leaves it (ParOperator::ParallelAssemble).
      hypre_CSRMatrix *diag = hypre_ParCSRMatrixDiag((hypre_ParCSRMatrix *)A);
      hypre_CSRMatrixMigrate(diag, HYPRE_MEMORY_HOST);
      const HYPRE_Int *I = hypre_CSRMatrixI(diag);
      const HYPRE_Int *J = hypre_CSRMatrixJ(diag);
      bool diagonal_first = true;
      for (HYPRE_Int i = 0; i < hypre_CSRMatrixNumRows(diag); i++)
      {
        diagonal_first = diagonal_first && (I[i] == I[i + 1] || J[I[i]] == i);
      }
      CHECK(diagonal_first);
    }
    Mpi::GlobalMax(1, &max_rel_diff, comm);
    // Rounding level: the combination applies the term scalings after the assembly of the
    // term matrices instead of before (the largest measured here is 1.2e-15, a few ULP).
    CHECK(max_rel_diff <= 1.0e-13);

    // The combined level applies as the complex operator of its two assembled parts, which
    // is what the reference level applies matrix-free.
    ComplexVector x(coarse_op->Width()), y(coarse_op->Height()), y_ref(coarse_op->Height());
    x.UseDevice(true);
    y.UseDevice(true);
    y_ref.UseDevice(true);
    linalg::SetRandom(comm, x);
    coarse_op->Mult(x, y);
    coarse_ref_op->Mult(x, y_ref);
    y.AXPY(-1.0, y_ref);
    CHECK_THAT(linalg::Norml2(comm, y) / linalg::Norml2(comm, y_ref),
               Catch::Matchers::WithinAbs(0.0, 1.0e-12));
  };

  // The parallel registration exercises the collective cache and matrix assembly once. The
  // scalar variants below are rank-independent and run only in the serial registration.
  if (Mpi::Size(comm) > 1)
  {
    SpaceOperatorExposed space_op(solver, domains, boundaries, ProblemType::DRIVEN, units,
                                  mesh);
    SpaceOperator reference_op(solver, domains, reference_boundaries, ProblemType::DRIVEN,
                               units, mesh);
    CheckCoarseLevel(space_op, reference_op, 2.0);
    CHECK(space_op.pc_terms_built);
    return;
  }

  SECTION("Lossy materials")
  {
    SpaceOperatorExposed space_op(solver, domains, boundaries, ProblemType::DRIVEN, units,
                                  mesh);
    SpaceOperator reference_op(solver, domains, reference_boundaries, ProblemType::DRIVEN,
                               units, mesh);
    CheckCoarseLevel(space_op, reference_op, 2.0);
    REQUIRE(space_op.pc_terms_built);
    const auto *stiffness = space_op.pc_term_coarse.stiffness.get();
    REQUIRE(stiffness);
    CheckCoarseLevel(space_op, reference_op, 3.5);
    CHECK(space_op.pc_term_coarse.stiffness.get() == stiffness);
  }

  SECTION("Coarse solver handoff")
  {
    SpaceOperator space_op(solver, domains, boundaries, ProblemType::DRIVEN, units, mesh);
    const std::complex<double> a1 = 2.0i;
    auto pc =
        space_op.GetPreconditionerMatrix<ComplexOperator>(1.0 + 0.0i, a1, a1 * a1, 2.0);
    const auto *coarse_op = CoarseOperator(*pc);
    auto recording = std::make_unique<RecordingSolver>();
    auto *recording_ptr = recording.get();
    MfemWrapperSolver<ComplexOperator> coarse_solver(std::move(recording), false, true,
                                                     false);
    coarse_solver.SetOperator(*coarse_op);
    CHECK(coarse_solver.Height() == coarse_op->Height());
    CHECK(coarse_solver.Width() == coarse_op->Width());
    CHECK(recording_ptr->operator_height == 2 * coarse_op->Height());
    CHECK(recording_ptr->operator_width == 2 * coarse_op->Width());
    ComplexVector x(coarse_op->Width()), y(coarse_op->Height());
    x.UseDevice(true);
    y.UseDevice(true);
    linalg::SetRandom(comm, x);
    coarse_solver.Mult(x, y);
    CHECK(linalg::Norml2(comm, y) == 0.0);
  }

  SECTION("Lossless materials")
  {
    // Without a loss tangent or a conductivity the imaginary mass term has an exactly zero
    // coefficient and no cached matrix, while the per-frequency assembly of the level keeps
    // its configured, exactly zero block: the pattern must be the same either way.
    auto lossless = lossy;
    lossless.attributes = {1, 2};
    lossless.tandelta.s = {0.0, 0.0, 0.0};
    auto lossless_domains = domains;
    lossless_domains.materials = {lossless};
    SpaceOperator space_op(solver, lossless_domains, boundaries, ProblemType::DRIVEN, units,
                           mesh);
    SpaceOperator reference_op(solver, lossless_domains, reference_boundaries,
                               ProblemType::DRIVEN, units, mesh);
    REQUIRE(!space_op.GetMaterialOp().HasLossTangent());
    REQUIRE(!space_op.GetMaterialOp().HasConductivity());
    CheckCoarseLevel(space_op, reference_op, 2.0);
  }

  for (bool frequency_scaled : {false, true})
  {
    DYNAMIC_SECTION("Floquet wave vector, frequency scaled = " << frequency_scaled)
    {
      auto floquet_boundaries = boundaries;
      floquet_boundaries.periodic.wave_vector = {0.1, 0.0, 0.0};
      floquet_boundaries.periodic.floquet_reference_freq = frequency_scaled ? 1.0 : 0.0;
      auto floquet_reference_boundaries = floquet_boundaries;
      floquet_reference_boundaries.farfield.order = 2;
      SpaceOperatorExposed space_op(solver, domains, floquet_boundaries,
                                    ProblemType::DRIVEN, units, mesh);
      SpaceOperator reference_op(solver, domains, floquet_reference_boundaries,
                                 ProblemType::DRIVEN, units, mesh);
      REQUIRE(space_op.GetMaterialOp().HasWaveVector());
      CheckCoarseLevel(space_op, reference_op, 2.0);
      CHECK(space_op.pc_terms_built);
    }
  }
}

}  // namespace palace
