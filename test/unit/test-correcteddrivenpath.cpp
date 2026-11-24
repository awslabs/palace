// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include <complex>
#include <memory>
#include <vector>
#include <catch2/catch_test_macros.hpp>
#include <nlohmann/json.hpp>
#include <catch2/generators/catch_generators.hpp>
#include "drivers/drivensolver.hpp"
#include "fem/mesh.hpp"
#include "fixtures.hpp"
#include "linalg/ksp.hpp"
#include "linalg/operator.hpp"
#include "linalg/vector.hpp"
#include "models/spaceoperator.hpp"
#include "utils/communication.hpp"
#include "utils/filesystem.hpp"
#include "utils/geodata.hpp"
#include "utils/iodata.hpp"

// ERRATUM E-3 (b) (decision 1582 (1)): the response-corrected driven solve must see the
// same system on upstream's matrix-free system-operator path (SpaceOperator::
// GetSystemOperator, applied through the finest level of the preconditioner) and on the
// assembled fallback. The case drives the uniform sweep's own construction
// (internal::BuildDrivenSystemOperator, called by DrivenSolver::SweepUniform) on both
// paths, so it fails if the response terms are dropped on either one.

namespace palace
{

using json = nlohmann::json;
using namespace std::complex_literals;

namespace
{

// The unit cube (6 tetrahedra) with a lumped port on attribute 100 (z = 0) and a PEC face
// on attribute 101 (z = 1), at Order 2. The cube is 1 cm (L0), so at 4 GHz omega^2 dM is
// comparable to the stiffness and the response terms are well above round-off.
// MGMaxLevels 2 builds a p-multigrid preconditioner whose finest level is the system
// operator (the matrix-free path); MGMaxLevels 1 leaves GetSystemOperator empty (the
// assembled path). The cube has no interior boundary, so the crack policy does not act;
// the key is stated as in the RomOperator response cases.
json DrivenPathConfig(const fs::path &output, int mg_max_levels)
{
  json config;
  config["Problem"] = {{"Type", "Driven"}, {"Verbose", 0}, {"Output", output}};
  config["Model"] = {
      {"Mesh", fs::path(PALACE_TEST_DATA_DIR) / "lumpedport_mesh/cube_mesh_1_1_1_tet.msh"},
      {"L0", 1.0e-2},
      {"Refinement", json::object()},
      {"CrackInternalBoundaryElements", false}};
  config["Domains"] = {
      {"Materials", json::array({json::object({{"Attributes", json::array({1})},
                                               {"Permeability", 1.0},
                                               {"Permittivity", 1.0},
                                               {"LossTan", 0.0}})})}};
  config["Boundaries"] = {
      {"PEC", {{"Attributes", json::array({101})}}},
      {"LumpedPort", json::array({json::object({{"Index", 1},
                                                {"R", 50.0},
                                                {"Excitation", uint(1)},
                                                {"Attributes", json::array({100})},
                                                {"Direction", "+X"}})})}};
  config["Solver"] = {{"Order", 2},
                      {"Device", "CPU"},
                      {"Driven",
                       {{"Samples", json::array({json::object({{"Type", "Linear"},
                                                               {"MinFreq", 2.0},
                                                               {"MaxFreq", 4.0},
                                                               {"FreqStep", 2.0}})})}}},
                      {"Linear",
                       {{"Type", "Default"},
                        {"KSPType", "GMRES"},
                        {"MaxIts", 200},
                        {"Tol", 1.0e-10},
                        {"MGMaxLevels", mg_max_levels},
                        {"PCMatReal", false},
                        {"PCMatShifted", false}}},
                      {"SurfaceResponseCorrection",
                       {{"Library", (output / "unused-process-library.json").string()},
                        {"Developer", true},
                        {"CorrectionMode", "Both"},
                        {"SolveTol", 1.0e-10}}}};
  return config;
}

std::vector<std::unique_ptr<Mesh>> LoadDrivenPathMesh(IoData &iodata, MPI_Comm comm)
{
  std::vector<std::unique_ptr<Mesh>> mesh;
  auto smesh = mesh::Load(iodata, comm);
  int ne = smesh ? smesh->GetNE() : 0;
  Mpi::GlobalMax(1, &ne, comm);
  if (Mpi::Size(comm) > ne)
  {
    SKIP("Mesh has fewer elements than processes");
  }
  if (iodata.model.Lc <= 0.0)
  {
    iodata.model.Lc = mesh::ComputeReferenceLength(smesh, comm);
  }
  iodata.NondimensionalizeInputs(smesh);
  std::vector<std::unique_ptr<mfem::ParMesh>> mfem_mesh;
  mfem_mesh.push_back(mesh::Partition(iodata, std::move(smesh), comm));
  mesh::RefineMesh(iodata, mfem_mesh);
  for (auto &m : mfem_mesh)
  {
    mesh.push_back(std::make_unique<Mesh>(std::move(m)));
  }
  return mesh;
}

double RelativeDifference(MPI_Comm comm, const ComplexVector &a, const ComplexVector &b)
{
  ComplexVector difference(a);
  linalg::AXPY(-1.0, b, difference);
  const double norm = linalg::Norml2(comm, b);
  REQUIRE(norm > 0.0);
  return linalg::Norml2(comm, difference) / norm;
}

}  // namespace

TEST_CASE_METHOD(palace::test::SharedTempDir,
                 "Response-corrected driven solve: matrix-free equals assembled",
                 "[drivensolver][correctedsolve][Serial][Parallel]")
{
  // The response mass dM alone, and dM with the magnetic stiffness dK beside it. The
  // stand-ins are the (DIAG_ZERO) mass matrix, as in the RomOperator response cases.
  const bool with_stiffness = GENERATE(false, true);
  MPI_Comm comm = Mpi::World();

  IoData iodata_mf(DrivenPathConfig(temp_dir / "matrix-free", 2), false);
  IoData iodata_asm(DrivenPathConfig(temp_dir / "assembled", 1), false);
  auto mesh_mf = LoadDrivenPathMesh(iodata_mf, comm);
  auto mesh_asm = LoadDrivenPathMesh(iodata_asm, comm);
  SpaceOperator space_op_mf(iodata_mf, mesh_mf);
  SpaceOperator space_op_asm(iodata_asm, mesh_asm);
  REQUIRE(space_op_mf.GetNDSpaces().GetNumLevels() == 2);
  REQUIRE(space_op_asm.GetNDSpaces().GetNumLevels() == 1);
  const int n = space_op_mf.GetNDSpace().GetTrueVSize();
  REQUIRE(space_op_asm.GetNDSpace().GetTrueVSize() == n);
  const double omega = iodata_mf.solver.driven.sample_f.back();
  REQUIRE(iodata_asm.solver.driven.sample_f.back() == omega);

  auto dM_mf = space_op_mf.GetMassMatrix<ComplexOperator>(Operator::DIAG_ZERO);
  auto dM_asm = space_op_asm.GetMassMatrix<ComplexOperator>(Operator::DIAG_ZERO);
  REQUIRE(dM_mf);
  REQUIRE(dM_asm);
  std::unique_ptr<ComplexOperator> dK_mf, dK_asm;
  if (with_stiffness)
  {
    dK_mf = space_op_mf.GetMassMatrix<ComplexOperator>(Operator::DIAG_ZERO);
    dK_asm = space_op_asm.GetMassMatrix<ComplexOperator>(Operator::DIAG_ZERO);
  }

  // The preconditioner matrices of the sweep; the matrix-free path is really taken on the
  // multigrid configuration and not on the other.
  auto P_mf = space_op_mf.GetPreconditionerMatrix<ComplexOperator>(
      1.0 + 0.0i, 1i * omega, -omega * omega + 0.0i, omega);
  auto P_asm = space_op_asm.GetPreconditionerMatrix<ComplexOperator>(
      1.0 + 0.0i, 1i * omega, -omega * omega + 0.0i, omega);
  REQUIRE(space_op_mf.GetSystemOperator(omega, *P_mf));
  REQUIRE_FALSE(space_op_asm.GetSystemOperator(omega, *P_asm));

  std::unique_ptr<ComplexOperator> K_mf, C_mf, M_mf, K_asm, C_asm, M_asm;
  auto system_mf = internal::BuildDrivenSystemOperator(
      space_op_mf, omega, *P_mf, K_mf, C_mf, M_mf, dM_mf.get(), dK_mf.get());
  auto system_asm = internal::BuildDrivenSystemOperator(
      space_op_asm, omega, *P_asm, K_asm, C_asm, M_asm, dM_asm.get(), dK_asm.get());
  REQUIRE(system_mf.A);
  REQUIRE(system_asm.A);
  CHECK_FALSE(system_mf.assembled);
  CHECK_FALSE(K_mf);  // The matrix-free path assembles no K, C, M.
  CHECK(system_asm.assembled);
  CHECK(K_asm);
  REQUIRE(system_mf.corrected);
  REQUIRE(system_asm.corrected);

  ComplexVector x(n), y_mf(n), y_asm(n);
  x.UseDevice(true);
  y_mf.UseDevice(true);
  y_asm.UseDevice(true);
  linalg::SetRandom(comm, x);

  // The raw system operators agree (this also shows that both configurations number the
  // true dofs alike, so the vectors below compare entry by entry).
  system_mf.A->Mult(x, y_mf);
  system_asm.A->Mult(x, y_asm);
  CHECK(RelativeDifference(comm, y_mf, y_asm) < 1.0e-12);

  // On each path the corrected operator is A - omega^2 dM (+ dK): the response terms are
  // present whichever operator A is.
  auto ResponseTermError = [&](const internal::DrivenSystemOperator &system,
                               const ComplexOperator &dM, const ComplexOperator *dK)
  {
    ComplexVector raw(n), corrected(n), response(n);
    raw.UseDevice(true);
    corrected.UseDevice(true);
    response.UseDevice(true);
    system.A->Mult(x, raw);
    system.corrected->Mult(x, corrected);
    linalg::AXPY(-1.0, raw, corrected);
    dM.Mult(x, response);
    response *= std::complex<double>(-omega * omega, 0.0);
    if (dK)
    {
      dK->AddMult(x, response, 1.0);
    }
    return RelativeDifference(comm, corrected, response);
  };
  CHECK(ResponseTermError(system_mf, *dM_mf, dK_mf.get()) < 1.0e-10);
  CHECK(ResponseTermError(system_asm, *dM_asm, dK_asm.get()) < 1.0e-10);

  // The corrected operators agree to round-off.
  system_mf.corrected->Mult(x, y_mf);
  system_asm.corrected->Mult(x, y_asm);
  CHECK(RelativeDifference(comm, y_mf, y_asm) < 1.0e-12);

  // The solves in the driver's form (the corrected solver: initial guess on, the
  // correction's SolveTol).
  auto Solve = [&](IoData &iodata, SpaceOperator &space_op, const ComplexOperator &A,
                   const ComplexOperator &P, bool corrected, ComplexVector &E)
  {
    ComplexKspSolver ksp(iodata, space_op.GetNDSpaces(), &space_op.GetH1Spaces());
    if (corrected)
    {
      ksp.SetInitialGuess(true);
      ksp.SetRelTol(iodata.solver.surface_response_correction->solve_tol);
    }
    ksp.SetOperators(A, P);
    ComplexVector RHS(n);
    RHS.UseDevice(true);
    space_op.GetExcitationVector(1, omega, RHS);
    E.SetSize(n);
    E.UseDevice(true);
    E = 0.0;
    ksp.Mult(RHS, E);
    return RHS;
  };
  ComplexVector E_raw_mf, E_mf, E_asm;
  Solve(iodata_mf, space_op_mf, *system_mf.A, *P_mf, false, E_raw_mf);
  Solve(iodata_mf, space_op_mf, *system_mf.corrected, *P_mf, true, E_mf);
  const auto RHS_asm =
      Solve(iodata_asm, space_op_asm, *system_asm.corrected, *P_asm, true, E_asm);

  // The corrected solutions agree within the solve tolerance, the matrix-free one solves
  // the assembled corrected system, and the correction moves the solution.
  CHECK(RelativeDifference(comm, E_mf, E_asm) < 1.0e-8);
  ComplexVector residual(n);
  residual.UseDevice(true);
  system_asm.corrected->Mult(E_mf, residual);
  CHECK(RelativeDifference(comm, residual, RHS_asm) < 1.0e-8);
  CHECK(RelativeDifference(comm, E_raw_mf, E_mf) > 1.0e-6);
}

}  // namespace palace
