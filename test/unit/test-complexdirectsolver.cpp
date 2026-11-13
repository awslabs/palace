// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include <algorithm>
#include <cmath>
#include <complex>
#include <memory>
#include <string>
#include <utility>
#include <vector>
#include <mfem.hpp>
#include <catch2/catch_test_macros.hpp>
#include <nlohmann/json.hpp>
#include <catch2/generators/catch_generators.hpp>
#include <catch2/generators/catch_generators_range.hpp>
#include "fem/mesh.hpp"
#include "linalg/ksp.hpp"
#include "linalg/mumps.hpp"
#include "linalg/mumpscomplex.hpp"
#include "linalg/operator.hpp"
#include "linalg/solver.hpp"
#include "linalg/strumpack.hpp"
#include "linalg/strumpackcomplex.hpp"
#include "linalg/superlu.hpp"
#include "linalg/superlucomplex.hpp"
#include "linalg/vector.hpp"
#include "models/spaceoperator.hpp"
#include "utils/communication.hpp"
#include "utils/configfile.hpp"
#include "utils/geodata.hpp"
#include "utils/iodata.hpp"
#include "utils/units.hpp"

#if defined(MFEM_USE_STRUMPACK) || defined(MFEM_USE_SUPERLU) || defined(MFEM_USE_MUMPS)

namespace palace
{

namespace
{

// Complex-valued indefinite matrix K - k² (1 - i τ(x)) M (+ i C with a convection matrix C
// for an unsymmetric matrix) of an H1 space on a hexahedral mesh of the unit cube, with a
// loss τ(x) in part of the domain, as a ComplexOperator with HypreParMatrix real and
// imaginary parts.
std::unique_ptr<ComplexOperator> AssembleTestMatrix(mfem::ParFiniteElementSpace &fespace,
                                                    double imag_scale = 1.0,
                                                    bool unsymmetric = false)
{
  constexpr double k2 = 20.0;
  mfem::ConstantCoefficient one(1.0), mk2(-k2);
  mfem::FunctionCoefficient tau([](const mfem::Vector &x)
                                { return (x(2) > 0.5) ? 5.0 * (x(2) - 0.5) : 0.0; });
  mfem::ProductCoefficient k2tau(k2 * imag_scale, tau);
  mfem::Vector b(3);
  b(0) = 1.0;
  b(1) = 2.0;
  b(2) = 3.0;
  mfem::VectorConstantCoefficient beta(b);
  mfem::ParBilinearForm ar(&fespace), ai(&fespace);
  ar.AddDomainIntegrator(new mfem::DiffusionIntegrator(one));
  ar.AddDomainIntegrator(new mfem::MassIntegrator(mk2));
  ai.AddDomainIntegrator(new mfem::MassIntegrator(k2tau));
  if (unsymmetric)
  {
    ai.AddDomainIntegrator(new mfem::ConvectionIntegrator(beta));
  }
  ar.Assemble();
  ar.Finalize();
  ai.Assemble();
  ai.Finalize();
  return std::make_unique<ComplexWrapperOperator>(
      std::unique_ptr<Operator>(ar.ParallelAssemble()),
      std::unique_ptr<Operator>(ai.ParallelAssemble()));
}

// Diagonal matrix (2 + i) M_L with the lumped mass matrix M_L, which has a different
// sparsity pattern than the matrices of AssembleTestMatrix.
std::unique_ptr<ComplexOperator>
AssembleDiagonalTestMatrix(mfem::ParFiniteElementSpace &fespace)
{
  mfem::ConstantCoefficient two(2.0), one(1.0);
  mfem::ParBilinearForm ar(&fespace), ai(&fespace);
  ar.AddDomainIntegrator(new mfem::LumpedIntegrator(new mfem::MassIntegrator(two)));
  ai.AddDomainIntegrator(new mfem::LumpedIntegrator(new mfem::MassIntegrator(one)));
  ar.Assemble();
  ar.Finalize();
  ai.Assemble();
  ai.Finalize();
  return std::make_unique<ComplexWrapperOperator>(
      std::unique_ptr<Operator>(ar.ParallelAssemble()),
      std::unique_ptr<Operator>(ai.ParallelAssemble()));
}

// The solution of the complex arithmetic solver matches the solution of the solver of the
// equivalent real-valued system, also for repeated factorizations of matrices with the same
// sparsity pattern and with a different sparsity pattern.
void TestComplexDirectSolver(Solver<ComplexOperator> &solver, Solver<ComplexOperator> &ref,
                             double tol, bool unsymmetric = false)
{
  MPI_Comm comm = Mpi::World();
  mfem::Mesh serial_mesh =
      mfem::Mesh::MakeCartesian3D(6, 6, 6, mfem::Element::HEXAHEDRON, 1.0, 1.0, 1.0);
  mfem::ParMesh mesh(comm, serial_mesh);
  mfem::H1_FECollection fec(2, 3);
  mfem::ParFiniteElementSpace fespace(&mesh, &fec);
  auto A = AssembleTestMatrix(fespace, 1.0, unsymmetric);

  ComplexVector b(A->Height()), x(A->Height()), x_ref(A->Height()), r(A->Height());
  b.Real().Randomize(1);
  b.Imag().Randomize(2);
  const double norm_b = linalg::Norml2(comm, b);
  auto Residual = [&](const ComplexOperator &M)
  {
    M.Mult(x, r);
    r -= b;
    return linalg::Norml2(comm, r) / norm_b;
  };

  solver.SetOperator(*A);
  solver.Mult(b, x);
  ref.SetOperator(*A);
  ref.Mult(b, x_ref);
  CHECK(Residual(*A) < tol);
  x -= x_ref;
  CHECK(linalg::Norml2(comm, x) / linalg::Norml2(comm, x_ref) < 10.0 * tol);

  // Same sparsity pattern, different values (reuse of the reordering).
  auto A2 = AssembleTestMatrix(fespace, 2.0, unsymmetric);
  solver.SetOperator(*A2);
  solver.Mult(b, x);
  CHECK(Residual(*A2) < tol);

  // Different sparsity pattern, and back.
  auto D = AssembleDiagonalTestMatrix(fespace);
  solver.SetOperator(*D);
  solver.Mult(b, x);
  CHECK(Residual(*D) < tol);
  solver.SetOperator(*A);
  solver.Mult(b, x);
  CHECK(Residual(*A) < tol);
}

}  // namespace

#if defined(MFEM_USE_STRUMPACK)
TEST_CASE("STRUMPACK complex arithmetic solver", "[strumpack][Serial][Parallel]")
{
  // Complex symmetric and unsymmetric matrices, with and without the reuse of the
  // reordering.
  MPI_Comm comm = Mpi::World();
  const bool unsymmetric = GENERATE(false, true);
  const bool reorder_reuse = GENERATE(true, false);
  StrumpackComplexSolver solver(comm, SymbolicFactorization::DEFAULT, reorder_reuse, 0);
  MfemWrapperSolver<ComplexOperator> ref(
      std::make_unique<StrumpackSolver>(comm, SymbolicFactorization::DEFAULT,
                                        SparseCompression::NONE, 0.0, 0, 0, true, 0),
      true, true, false, true);
  TestComplexDirectSolver(solver, ref, 1.0e-11, unsymmetric);
}

// As for StrumpackSolver (see test-strumpack.cpp), the parallel matrix products of hypre
// can store explicit zeros without their transposes, so the sparsity pattern given to
// STRUMPACK is not declared symmetric: a five-point Laplacian with a loss, plus explicit
// zeros which couple distant unknowns in one direction only.
TEST_CASE("STRUMPACK complex arithmetic solver unsymmetric sparsity pattern",
          "[strumpack][Serial][Parallel]")
{
  MPI_Comm comm = Mpi::World();
  const int P = Mpi::Size(comm), rank = Mpi::Rank(comm);
  const int m = 24;
  const HYPRE_BigInt N = m * m;
  const HYPRE_BigInt first = (N * rank) / P, last = (N * (rank + 1)) / P;
  const int n_loc = static_cast<int>(last - first);
  std::vector<int> Ir(n_loc + 1, 0), Ii(n_loc + 1, 0);
  std::vector<HYPRE_BigInt> Jr, Ji;
  std::vector<double> Vr, Vi;
  for (HYPRE_BigInt i = first; i < last; i++)
  {
    const HYPRE_BigInt ix = i % m, iy = i / m;
    std::vector<std::pair<HYPRE_BigInt, double>> row = {{i, 4.0}};
    if (ix > 0)
    {
      row.emplace_back(i - 1, -1.0);
    }
    if (ix < m - 1)
    {
      row.emplace_back(i + 1, -1.0);
    }
    if (iy > 0)
    {
      row.emplace_back(i - m, -1.0);
    }
    if (iy < m - 1)
    {
      row.emplace_back(i + m, -1.0);
    }
    if (i < N / 2 && i % 5 == 0)
    {
      row.emplace_back(N - 1 - i, 0.0);
    }
    if (i >= N / 2 && i % 7 == 3 && (i - N / 2) % 5 != 0)
    {
      row.emplace_back(i - N / 2, 0.0);
    }
    std::sort(row.begin(), row.end());
    for (const auto &[j, v] : row)
    {
      Jr.push_back(j);
      Vr.push_back(v);
    }
    Ir[i - first + 1] = static_cast<int>(Jr.size());
    Ji.push_back(i);
    Vi.push_back(0.5);
    Ii[i - first + 1] = static_cast<int>(Ji.size());
  }
  HYPRE_BigInt row_starts[2] = {first, last};
  auto *Ar = new mfem::HypreParMatrix(comm, n_loc, N, N, Ir.data(), Jr.data(), Vr.data(),
                                      row_starts, row_starts);
  auto *Ai = new mfem::HypreParMatrix(comm, n_loc, N, N, Ii.data(), Ji.data(), Vi.data(),
                                      row_starts, row_starts);
  ComplexWrapperOperator A{std::unique_ptr<Operator>(Ar), std::unique_ptr<Operator>(Ai)};

  ComplexVector x_ex(n_loc), b(n_loc), x(n_loc), r(n_loc);
  for (int i = 0; i < n_loc; i++)
  {
    const double t = static_cast<double>(first + i) / static_cast<double>(N);
    x_ex.Real()[i] = std::sin(7.0 * t) + t;
    x_ex.Imag()[i] = std::cos(3.0 * t);
  }
  A.Mult(x_ex, b);
  StrumpackComplexSolver solver(comm, SymbolicFactorization::DEFAULT, true, 0);
  for (int k = 0; k < 2; k++)
  {
    // Factorization, and refactorization reusing the reordering.
    solver.SetOperator(A);
    x = 0.0;
    solver.Mult(b, x);
    r = x;
    r -= x_ex;
    CHECK(linalg::Norml2(comm, r) / linalg::Norml2(comm, x_ex) < 1.0e-12);
  }
}
#endif

#if defined(MFEM_USE_SUPERLU)
TEST_CASE("SuperLU complex arithmetic solver", "[superlu][Serial][Parallel]")
{
  // Complex symmetric and unsymmetric matrices, with and without the reuse of the row and
  // column permutations.
  MPI_Comm comm = Mpi::World();
  const bool unsymmetric = GENERATE(false, true);
  const bool reorder_reuse = GENERATE(true, false);
  SuperLUComplexSolver solver(comm, SymbolicFactorization::DEFAULT, reorder_reuse, 0);
  MfemWrapperSolver<ComplexOperator> ref(
      std::make_unique<SuperLUSolver>(comm, SymbolicFactorization::DEFAULT, false, true, 0),
      true, true, false, true);
  TestComplexDirectSolver(solver, ref, 1.0e-10, unsymmetric);
}
#endif

#if defined(MFEM_USE_MUMPS)
TEST_CASE("MUMPS complex arithmetic solver", "[mumps][Serial][Parallel]")
{
  // Complex symmetric and unsymmetric matrices.
  MPI_Comm comm = Mpi::World();
  const auto sym = GENERATE(MatrixSymmetry::SYMMETRIC, MatrixSymmetry::UNSYMMETRIC);
  MumpsComplexSolver solver(comm, sym, SymbolicFactorization::DEFAULT, 0.0, true, 0);
  MfemWrapperSolver<ComplexOperator> ref(
      std::make_unique<MumpsSolver>(comm, sym, SymbolicFactorization::DEFAULT, 0.0, true,
                                    0),
      true, true, false, true);
  TestComplexDirectSolver(solver, ref, 1.0e-10, sym == MatrixSymmetry::UNSYMMETRIC);
}
#endif

namespace
{

using json = nlohmann::json;

// The sparse direct solver types of this build ("Type" of config["Solver"]["Linear"]).
std::vector<std::string> SparseDirectSolverTypes()
{
  std::vector<std::string> types;
#if defined(MFEM_USE_STRUMPACK)
  types.push_back("STRUMPACK");
#endif
#if defined(MFEM_USE_SUPERLU)
  types.push_back("SuperLU");
#endif
#if defined(MFEM_USE_MUMPS)
  types.push_back("MUMPS");
#endif
  return types;
}

// Whether the solver is the complex arithmetic solver of the sparse direct solver type.
bool IsComplexArithmeticSolver(const std::string &type, Solver<ComplexOperator> *pc)
{
#if defined(MFEM_USE_STRUMPACK)
  if (type == "STRUMPACK")
  {
    return dynamic_cast<StrumpackComplexSolver *>(pc) != nullptr;
  }
#endif
#if defined(MFEM_USE_SUPERLU)
  if (type == "SuperLU")
  {
    return dynamic_cast<SuperLUComplexSolver *>(pc) != nullptr;
  }
#endif
#if defined(MFEM_USE_MUMPS)
  if (type == "MUMPS")
  {
    return dynamic_cast<MumpsComplexSolver *>(pc) != nullptr;
  }
#endif
  return false;
}

// ComplexKspSolver whose preconditioner can be inspected.
class InspectableComplexKspSolver : public ComplexKspSolver
{
public:
  using ComplexKspSolver::ComplexKspSolver;

  Solver<ComplexOperator> *GetPreconditioner() { return pc.get(); }
};

// Driven problem on a cube of lossy dielectric with PEC walls at z = 0 and z = 1 (essential
// boundary conditions) and a lumped port on x = 0, without multigrid: the preconditioner
// matrix is a ComplexParOperator with real and imaginary parts and eliminated essential
// degrees of freedom, built and solved as in a simulation.
struct DrivenTestProblem
{
  IoData iodata{Units(1.0e-3, 1.0e-3)};
  std::vector<std::unique_ptr<Mesh>> mesh;
  std::unique_ptr<SpaceOperator> space_op;

  DrivenTestProblem(const json &linear)
  {
    iodata.model.L0 = 1.0e-3;
    iodata.model.Lc = 1.0;
    iodata.problem.type = ProblemType::DRIVEN;
    iodata.solver.order = 2;
    iodata.solver.driven.sample_f = {1.0};
    iodata.solver.linear = config::LinearSolverData(linear);
    iodata.domains = config::DomainData(
        json{{"Materials", json::array({json::object({{"Attributes", json::array({1})},
                                                      {"Permittivity", 4.0},
                                                      {"LossTan", 0.3}})})}});
    // MakeCartesian3D boundary attributes: 1 (z = 0), 6 (z = sz), 5 (x = 0).
    iodata.boundaries = config::BoundaryData(
        json{{"PEC", json::object({{"Attributes", json::array({1, 6})}})},
             {"LumpedPort", json::array({json::object({{"Index", 1},
                                                       {"R", 50.0},
                                                       {"Excitation", 1},
                                                       {"Attributes", json::array({5})},
                                                       {"Direction", "+Z"}})})}});
    iodata.CheckConfiguration();
    auto smesh = std::make_unique<mfem::Mesh>(
        mfem::Mesh::MakeCartesian3D(4, 4, 4, mfem::Element::HEXAHEDRON, 1.0, 1.0, 1.0));
    iodata.NondimensionalizeInputs(smesh);
    mesh.push_back(
        std::make_unique<Mesh>(mesh::Partition(iodata, std::move(smesh), Mpi::World())));
    space_op = std::make_unique<SpaceOperator>(iodata, mesh);
  }

  std::unique_ptr<InspectableComplexKspSolver> MakeKspSolver()
  {
    return std::make_unique<InspectableComplexKspSolver>(iodata, space_op->GetNDSpaces(),
                                                         &space_op->GetH1Spaces());
  }

  // Solve the driven problem at the given multiple of the sample frequency, from a zero
  // initial guess. Returns the number of linear solver iterations.
  int Solve(InspectableComplexKspSolver &ksp, double f_scale, ComplexVector &E)
  {
    using namespace std::complex_literals;
    const double omega = f_scale * iodata.solver.driven.sample_f.front();
    auto K = space_op->GetStiffnessMatrix<ComplexOperator>(Operator::DIAG_ONE);
    auto C = space_op->GetDampingMatrix<ComplexOperator>(Operator::DIAG_ZERO);
    auto M = space_op->GetMassMatrix<ComplexOperator>(Operator::DIAG_ZERO);
    auto A2 = space_op->GetExtraSystemOperator(omega, Operator::DIAG_ZERO);
    auto A = space_op->GetSystemMatrix(1.0 + 0.0i, 1i * omega, -omega * omega + 0.0i,
                                       K.get(), C.get(), M.get(), A2.get());
    auto P = space_op->GetPreconditionerMatrix<ComplexOperator>(
        1.0 + 0.0i, 1i * omega, -omega * omega + 0.0i, omega);
    ksp.SetOperators(*A, *P);
    ComplexVector RHS;
    RHS.UseDevice(true);
    REQUIRE(space_op->GetExcitationVector(1, omega, RHS));
    E.SetSize(RHS.Size());
    E.UseDevice(true);
    E = 0.0;
    const int it = ksp.NumTotalMultIterations();
    ksp.Mult(RHS, E);
    CHECK(ksp.GetConverged());
    return ksp.NumTotalMultIterations() - it;
  }
};

json LinearOptions(const std::string &type, json options)
{
  options["Type"] = type;
  options["KSPType"] = "GMRES";
  options["Tol"] = 1.0e-10;
  options["MaxIts"] = 200;
  options["MGMaxLevels"] = 1;
  return options;
}

}  // namespace

TEST_CASE("ComplexKspSolver sparse direct solver with ComplexCoarseSolve",
          "[strumpack][superlu][mumps][Serial][Parallel]")
{
  // ConfigurePreconditionerSolver uses the complex arithmetic solvers where they apply,
  // and they solve the ComplexParOperator of a simulation (assembled again for a new
  // operator after the previous assembly was released, essential boundary conditions
  // eliminated) as the equivalent real-valued system of twice the size does.
  MPI_Comm comm = Mpi::World();
  const auto type = GENERATE(from_range(SparseDirectSolverTypes()));
  INFO("type " << type);

  // Reference: the equivalent real-valued system of twice the size (DropSmallEntries keeps
  // it with ComplexCoarseSolve).
  ComplexVector E_ref;
  int n = 0;
  {
    DrivenTestProblem problem(LinearOptions(
        type, json{{"ComplexCoarseSolve", true}, {"DropSmallEntries", true}}));
    auto ksp = problem.MakeKspSolver();
    auto *pc = dynamic_cast<MfemWrapperSolver<ComplexOperator> *>(ksp->GetPreconditioner());
    REQUIRE(pc);
    CHECK(problem.Solve(*ksp, 1.0, E_ref) <= 2);
    n = E_ref.Size();
    CHECK(pc->GetSolver().Height() == 2 * n);
  }

  // Complex arithmetic, also with the ParMETIS column ordering.
  for (const std::string ordering : {"Default", "ParMETIS"})
  {
    INFO("ordering " << ordering);
    DrivenTestProblem problem(LinearOptions(
        type, json{{"ComplexCoarseSolve", true}, {"ColumnOrdering", ordering}}));
    auto ksp = problem.MakeKspSolver();
    CHECK(IsComplexArithmeticSolver(type, ksp->GetPreconditioner()));
    ComplexVector E;
    CHECK(problem.Solve(*ksp, 1.0, E) <= 2);
    E -= E_ref;
    CHECK(linalg::Norml2(comm, E) / linalg::Norml2(comm, E_ref) < 1.0e-8);
    CHECK(problem.Solve(*ksp, 1.2, E) <= 2);
  }

  // The cases which keep the equivalent real-valued system of twice the size.
  std::vector<json> block_options = {
      json{{"ComplexCoarseSolve", true}, {"PCMatReal", true}},
      json{{"ComplexCoarseSolve", true}, {"DropSmallEntries", true}}};
  if (type == "SuperLU")
  {
    block_options.push_back(
        json{{"ComplexCoarseSolve", true}, {"SuperLU3DCommunicator", true}});
  }
  if (type == "STRUMPACK")
  {
    block_options.push_back(
        json{{"ComplexCoarseSolve", true}, {"STRUMPACKCompressionType", "BLR"}});
  }
  for (const auto &options : block_options)
  {
    INFO("options " << options.dump());
    DrivenTestProblem problem(LinearOptions(type, options));
    auto ksp = problem.MakeKspSolver();
    CHECK(dynamic_cast<MfemWrapperSolver<ComplexOperator> *>(ksp->GetPreconditioner()));
  }

  // The real-valued approximation.
  {
    DrivenTestProblem problem(LinearOptions(type, json{{"ComplexCoarseSolve", false}}));
    auto ksp = problem.MakeKspSolver();
    auto *pc = dynamic_cast<MfemWrapperSolver<ComplexOperator> *>(ksp->GetPreconditioner());
    REQUIRE(pc);
    ComplexVector E;
    problem.Solve(*ksp, 1.0, E);
    CHECK(pc->GetSolver().Height() == n);
  }
}

}  // namespace palace

#endif
