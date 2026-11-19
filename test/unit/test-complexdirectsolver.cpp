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
#include "linalg/complexcoarse.hpp"
#include "linalg/complexcsr.hpp"
#include "linalg/iterative.hpp"
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
#include "test-helpers.hpp"
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
  // reordering. The solvers store the permutation once per node, when STRUMPACK is built
  // with that option: with STRUMPACK_SHARED_PERM_GROUP=1, every process acts as a node, so
  // that on several processes the root broadcasts the permutation to the other node leaders
  // as on several nodes. STRUMPACK reads the variable on every process at each reordering,
  // so it is set before the solvers are constructed and restored after they are destroyed.
  MPI_Comm comm = Mpi::World();
  const bool unsymmetric = GENERATE(false, true);
  const bool reorder_reuse = GENERATE(true, false);
  const bool process_as_node = GENERATE(false, true);
  CAPTURE(unsymmetric, reorder_reuse, process_as_node);
  ScopedEnvironmentVariable group("STRUMPACK_SHARED_PERM_GROUP",
                                  process_as_node ? "1" : nullptr);
  StrumpackComplexSolver solver(comm, SymbolicFactorization::DEFAULT,
                                SparseCompression::NONE, 0.0, 0, 0, false, reorder_reuse,
                                0);
  MfemWrapperSolver<ComplexOperator> ref(
      std::make_unique<StrumpackSolver>(comm, SymbolicFactorization::DEFAULT,
                                        SparseCompression::NONE, 0.0, 0, 0, true, 0),
      true, true, false, true);
  TestComplexDirectSolver(solver, ref, 1.0e-11, unsymmetric);
}

TEST_CASE("STRUMPACK complex arithmetic solver with compression",
          "[strumpack][Serial][Parallel]")
{
  // Every compression type of StrumpackSolver: the types with a compression tolerance on a
  // matrix whose largest separators exceed the minimum size for BLR compression, and the
  // ZFP types, which compress small fronts too, on a smaller one. STRUMPACK factors the
  // matrix without compression for the types which need larger separators (HSS, HODLR) or
  // a library it was not built with (ButterflyPACK for the HODLR types, ZFP for the ZFP
  // types). The compressed factorization is checked as a solver, within the accuracy of
  // the compression (which STRUMPACK fixes for the lossy ZFP compression), and as the
  // preconditioner of GMRES, also when it is computed again for a matrix with the same
  // sparsity pattern and with another one.
  MPI_Comm comm = Mpi::World();

  // Compression type, lossy precision (0 is lossless), bound for the residual of one solve
  // (with a compression tolerance of 1e-6), and for the GMRES iterations.
  struct Case
  {
    SparseCompression type;
    int lossy_prec;
    double tol;
    int max_it;
  };
  auto Test = [&](int mesh_n, const std::vector<Case> &cases)
  {
    mfem::Mesh serial_mesh = mfem::Mesh::MakeCartesian3D(
        mesh_n, mesh_n, mesh_n, mfem::Element::HEXAHEDRON, 1.0, 1.0, 1.0);
    mfem::ParMesh mesh(comm, serial_mesh);
    mfem::H1_FECollection fec(2, 3);
    mfem::ParFiniteElementSpace fespace(&mesh, &fec);
    const auto A = AssembleTestMatrix(fespace, 1.0);
    const auto A2 = AssembleTestMatrix(fespace, 2.0);
    const auto D = AssembleDiagonalTestMatrix(fespace);

    ComplexVector b(A->Height()), x(A->Height()), x_ref(A->Height()), r(A->Height());
    b.Real().Randomize(1);
    b.Imag().Randomize(2);
    const double norm_b = linalg::Norml2(comm, b);
    {
      StrumpackComplexSolver ref(comm, SymbolicFactorization::DEFAULT,
                                 SparseCompression::NONE, 0.0, 0, 0, false, true, 0);
      ref.SetOperator(*A);
      ref.Mult(b, x_ref);
    }
    for (const auto &[type, lossy_prec, tol, max_it] : cases)
    {
      INFO("compression " << static_cast<int>(type) << ", lossy precision " << lossy_prec);
      StrumpackComplexSolver solver(comm, SymbolicFactorization::DEFAULT, type, 1.0e-6, 1,
                                    lossy_prec, false, true, 0);
      GmresSolver<ComplexOperator> gmres(comm, 0);
      gmres.SetRelTol(1.0e-10);
      gmres.SetMaxIter(max_it);
      gmres.SetPreconditioner(solver);
      for (const auto *M : {A.get(), A2.get(), D.get()})
      {
        solver.SetOperator(*M);
        x = 0.0;
        solver.Mult(b, x);
        M->Mult(x, r);
        r -= b;
        CHECK(linalg::Norml2(comm, r) / norm_b < tol);
        gmres.SetOperator(*M);
        x = 0.0;
        gmres.Mult(b, x);
        CHECK(gmres.GetConverged());
        if (M == A.get())
        {
          x -= x_ref;
          CHECK(linalg::Norml2(comm, x) / linalg::Norml2(comm, x_ref) < 1.0e-8);
        }
      }
    }
  };
  Test(12, {{SparseCompression::HSS, 16, 1.0e-4, 5},
            {SparseCompression::BLR, 16, 1.0e-4, 5},
            {SparseCompression::HODLR, 16, 1.0e-4, 5},
            {SparseCompression::BLR_HODLR, 16, 1.0e-4, 5}});
  Test(6, {{SparseCompression::ZFP, 0, 1.0e-10, 2},
           {SparseCompression::ZFP, 16, 0.5, 30},
           {SparseCompression::ZFP_BLR_HODLR, 16, 0.5, 30}});
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
  StrumpackComplexSolver solver(comm, SymbolicFactorization::DEFAULT,
                                SparseCompression::NONE, 0.0, 0, 0, false, true, 0);
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
  // column permutations, on a 2D and on a 3D process grid (which has a depth above one
  // from 8 processes).
  MPI_Comm comm = Mpi::World();
  const bool unsymmetric = GENERATE(false, true);
  const bool reorder_reuse = GENERATE(true, false);
  const bool use_3d = GENERATE(false, true);
  SuperLUComplexSolver solver(comm, SymbolicFactorization::DEFAULT, use_3d, false,
                              reorder_reuse, 0);
  MfemWrapperSolver<ComplexOperator> ref(
      std::make_unique<SuperLUSolver>(comm, SymbolicFactorization::DEFAULT, use_3d, true,
                                      0),
      true, true, false, true);
  TestComplexDirectSolver(solver, ref, 1.0e-10, unsymmetric);
}
#endif

#if defined(MFEM_USE_MUMPS)
TEST_CASE("MUMPS complex arithmetic solver", "[mumps][Serial][Parallel]")
{
  // Complex symmetric and unsymmetric matrices, without and with block low-rank
  // compression (which the fronts of this matrix are too small for: the solution is the
  // same).
  MPI_Comm comm = Mpi::World();
  const auto sym = GENERATE(MatrixSymmetry::SYMMETRIC, MatrixSymmetry::UNSYMMETRIC);
  const double blr_tol = GENERATE(0.0, 1.0e-8);
  MumpsComplexSolver solver(comm, sym, SymbolicFactorization::DEFAULT, blr_tol, false, true,
                            0);
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
// degrees of freedom, built and solved as in a simulation. Without losses (a lossless
// dielectric and an inductive port), its imaginary part is assembled as zeros.
struct DrivenTestProblem
{
  IoData iodata{Units(1.0e-3, 1.0e-3)};
  std::vector<std::unique_ptr<Mesh>> mesh;
  std::unique_ptr<SpaceOperator> space_op;

  DrivenTestProblem(const json &linear, bool lossless = false)
  {
    iodata.model.L0 = 1.0e-3;
    iodata.model.Lc = 1.0;
    iodata.problem.type = ProblemType::DRIVEN;
    iodata.solver.order = 2;
    iodata.solver.driven.sample_f = {1.0};
    iodata.solver.linear = config::LinearSolverData(linear);
    iodata.domains = config::DomainData(json{
        {"Materials", json::array({json::object({{"Attributes", json::array({1})},
                                                 {"Permittivity", 4.0},
                                                 {"LossTan", lossless ? 0.0 : 0.3}})})}});
    // MakeCartesian3D boundary attributes: 1 (z = 0), 6 (z = sz), 5 (x = 0).
    iodata.boundaries = config::BoundaryData(
        json{{"PEC", json::object({{"Attributes", json::array({1, 6})}})},
             {"LumpedPort",
              json::array({json::object({{"Index", 1},
                                         {lossless ? "L" : "R", lossless ? 1.0e-9 : 50.0},
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

  // Reference: the solution with the real-valued approximation of the preconditioner
  // matrix.
  ComplexVector E_ref;
  int n = 0;
  {
    DrivenTestProblem problem(LinearOptions(type, json{{"ComplexCoarseSolve", false}}));
    auto ksp = problem.MakeKspSolver();
    problem.Solve(*ksp, 1.0, E_ref);
    n = E_ref.Size();
  }

  // Complex arithmetic, also with the ParMETIS column ordering, with small entries dropped
  // (the eliminated entries of the essential boundary conditions), for STRUMPACK with
  // compression (which the fronts of this matrix are too small for), and for SuperLU_DIST
  // with a 3D process grid.
  std::vector<json> complex_options = {
      json{{"ComplexCoarseSolve", true}},
      json{{"ComplexCoarseSolve", true}, {"ColumnOrdering", "ParMETIS"}},
      json{{"ComplexCoarseSolve", true}, {"DropSmallEntries", true}}};
  if (type == "STRUMPACK")
  {
    complex_options.push_back(
        json{{"ComplexCoarseSolve", true}, {"STRUMPACKCompressionType", "BLR"}});
  }
  if (type == "SuperLU")
  {
    complex_options.push_back(
        json{{"ComplexCoarseSolve", true}, {"SuperLU3DCommunicator", true}});
  }
  for (const auto &options : complex_options)
  {
    INFO("options " << options.dump());
    DrivenTestProblem problem(LinearOptions(type, options));
    auto ksp = problem.MakeKspSolver();
    auto *pc = dynamic_cast<ComplexCoarseSolver *>(ksp->GetPreconditioner());
    REQUIRE(pc);
    ComplexVector E;
    CHECK(problem.Solve(*ksp, 1.0, E) <= 2);
    CHECK(IsComplexArithmeticSolver(type, pc->GetComplexSolver()));
    CHECK(!pc->GetRealSolver());
    E -= E_ref;
    CHECK(linalg::Norml2(comm, E) / linalg::Norml2(comm, E_ref) < 1.0e-8);
    CHECK(problem.Solve(*ksp, 1.2, E) <= 2);
  }

  // A model without losses: the imaginary part of the preconditioner matrix is assembled
  // as zeros, and the real part alone is factored in real arithmetic (by the real-valued
  // solver, which drops the small entries itself), for the same solution as without
  // ComplexCoarseSolve.
  {
    ComplexVector E_ref0;
    {
      DrivenTestProblem problem(LinearOptions(type, json{{"ComplexCoarseSolve", false}}),
                                true);
      auto ksp = problem.MakeKspSolver();
      CHECK(problem.Solve(*ksp, 1.0, E_ref0) <= 2);
    }
    for (const auto &options :
         {json{{"ComplexCoarseSolve", true}},
          json{{"ComplexCoarseSolve", true}, {"DropSmallEntries", true}}})
    {
      INFO("options " << options.dump());
      DrivenTestProblem problem(LinearOptions(type, options), true);
      auto ksp = problem.MakeKspSolver();
      auto *pc = dynamic_cast<ComplexCoarseSolver *>(ksp->GetPreconditioner());
      REQUIRE(pc);
      ComplexVector E;
      CHECK(problem.Solve(*ksp, 1.0, E) <= 2);
      CHECK(!pc->GetComplexSolver());
      auto *real_pc =
          dynamic_cast<MfemWrapperSolver<ComplexOperator> *>(pc->GetRealSolver());
      REQUIRE(real_pc);
      CHECK(real_pc->GetSolver().Height() == n);
      E -= E_ref0;
      CHECK(linalg::Norml2(comm, E) / linalg::Norml2(comm, E_ref0) < 1.0e-8);
    }
  }

  // The mixed precision solver of STRUMPACK keeps the real-valued solver, for the
  // equivalent real-valued system of twice the size.
  if (type == "STRUMPACK")
  {
    DrivenTestProblem problem(
        LinearOptions("STRUMPACK-MP", json{{"ComplexCoarseSolve", true}}));
    auto ksp = problem.MakeKspSolver();
    auto *pc = dynamic_cast<MfemWrapperSolver<ComplexOperator> *>(ksp->GetPreconditioner());
    REQUIRE(pc);
    ComplexVector E;
    problem.Solve(*ksp, 1.0, E);
    CHECK(pc->GetSolver().Height() == 2 * n);
  }

  // The real-valued approximation, and the real-valued preconditioner matrix of PCMatReal
  // (also with ComplexCoarseSolve): the real-valued solver, for a matrix of the same size.
  for (const auto &options : {json{{"ComplexCoarseSolve", false}},
                              json{{"ComplexCoarseSolve", true}, {"PCMatReal", true}}})
  {
    INFO("options " << options.dump());
    DrivenTestProblem problem(LinearOptions(type, options));
    auto ksp = problem.MakeKspSolver();
    auto *pc = dynamic_cast<MfemWrapperSolver<ComplexOperator> *>(ksp->GetPreconditioner());
    REQUIRE(pc);
    ComplexVector E;
    problem.Solve(*ksp, 1.0, E);
    CHECK(pc->GetSolver().Height() == n);
  }
}

TEST_CASE("ComplexCoarseSolver zero imaginary part",
          "[strumpack][superlu][mumps][Serial][Parallel]")
{
  // The solver which ConfigurePreconditionerSolver builds for ComplexCoarseSolve factors
  // the real part alone, in real arithmetic, when the imaginary part is exactly zero, the
  // complex-valued matrix in complex arithmetic otherwise, and follows the operators it is
  // given, with one solver at a time.
  MPI_Comm comm = Mpi::World();
  const auto type = GENERATE(from_range(SparseDirectSolverTypes()));
  INFO("type " << type);
  DrivenTestProblem problem(LinearOptions(type, json{{"ComplexCoarseSolve", true}}));
  auto ksp = problem.MakeKspSolver();
  auto *pc = dynamic_cast<ComplexCoarseSolver *>(ksp->GetPreconditioner());
  REQUIRE(pc);
  CHECK((!pc->GetRealSolver() && !pc->GetComplexSolver()));

  // Matrices of TestComplexDirectSolver: with a loss, with an imaginary part of explicit
  // zeros, with an imaginary part which is tiny (its square underflows) and zero on part of
  // the domain, and without imaginary part.
  mfem::Mesh serial_mesh =
      mfem::Mesh::MakeCartesian3D(6, 6, 6, mfem::Element::HEXAHEDRON, 1.0, 1.0, 1.0);
  mfem::ParMesh mesh(comm, serial_mesh);
  mfem::H1_FECollection fec(2, 3);
  mfem::ParFiniteElementSpace fespace(&mesh, &fec);
  const auto A = AssembleTestMatrix(fespace, 1.0);
  const auto A_zero = AssembleTestMatrix(fespace, 0.0);
  const auto A_tiny = AssembleTestMatrix(fespace, 1.0e-300);
  const ComplexWrapperOperator A_real(A_zero->Real(), nullptr);
  REQUIRE(A_zero->Imag());

  const int n = A->Height();
  ComplexVector b(n), x(n), x_real(n), r(n);
  b.Real().Randomize(1);
  b.Imag().Randomize(2);
  const double norm_b = linalg::Norml2(comm, b);
  constexpr double tol = 1.0e-10;

  // Factors and solves, checks which solver is used (the real-valued one for a matrix of
  // size n, not for the real-valued system of twice the size), and returns the residual.
  auto Solve = [&](const ComplexOperator &M, bool real)
  {
    pc->SetOperator(M);
    if (real)
    {
      CHECK(!pc->GetComplexSolver());
      auto *real_pc =
          dynamic_cast<MfemWrapperSolver<ComplexOperator> *>(pc->GetRealSolver());
      REQUIRE(real_pc);
      CHECK(real_pc->GetSolver().Height() == n);
    }
    else
    {
      CHECK(!pc->GetRealSolver());
      CHECK(IsComplexArithmeticSolver(type, pc->GetComplexSolver()));
    }
    x = 0.0;
    pc->Mult(b, x);
    M.Mult(x, r);
    r -= b;
    return linalg::Norml2(comm, r) / norm_b;
  };
  auto DifferenceFromReal = [&]()
  {
    r = x;
    r -= x_real;
    return linalg::Norml2(comm, r) / linalg::Norml2(comm, x_real);
  };

  CHECK(Solve(*A_zero, true) < tol);
  x_real = x;
  CHECK(Solve(*A_tiny, false) < tol);
  CHECK(DifferenceFromReal() < 10.0 * tol);
  CHECK(Solve(*A, false) < tol);
  CHECK(Solve(A_real, true) < tol);
  CHECK(DifferenceFromReal() < 10.0 * tol);
  CHECK(Solve(*A, false) < tol);
  CHECK(Solve(*A_zero, true) < tol);
  CHECK(DifferenceFromReal() < 10.0 * tol);
}

namespace
{

// Value far below the tolerance of DropSmallEntries for a row with entries of order one.
constexpr double tiny_entry = 1.0e-40;

// The matrix of AssembleTestMatrix with the unknowns of the boundary eliminated from the
// rows and columns of its real and imaginary parts (unit diagonal), as for essential
// boundary conditions. The eliminated entries stay in the sparsity pattern: as zeros in
// the imaginary part and, in the real part, as tiny entries, so that the solution of a
// system with a zero right-hand side for the eliminated unknowns shows whether they were
// dropped (it is then exactly zero for these unknowns).
std::pair<std::unique_ptr<mfem::HypreParMatrix>, std::unique_ptr<mfem::HypreParMatrix>>
AssembleEliminatedTestMatrix(mfem::ParFiniteElementSpace &fespace,
                             const mfem::Array<int> &ess_tdofs)
{
  constexpr double k2 = 20.0;
  mfem::ConstantCoefficient one(1.0), mk2(-k2);
  mfem::FunctionCoefficient tau([](const mfem::Vector &x)
                                { return (x(2) > 0.5) ? 5.0 * (x(2) - 0.5) : 0.0; });
  mfem::ProductCoefficient k2tau(k2, tau);
  mfem::ParBilinearForm ar(&fespace), ai(&fespace);
  ar.AddDomainIntegrator(new mfem::DiffusionIntegrator(one));
  ar.AddDomainIntegrator(new mfem::MassIntegrator(mk2));
  ai.AddDomainIntegrator(new mfem::MassIntegrator(k2tau));
  ar.Assemble();
  ar.Finalize();
  ai.Assemble();
  ai.Finalize();
  std::unique_ptr<mfem::HypreParMatrix> Ar(ar.ParallelAssemble()),
      Ai(ai.ParallelAssemble());
  delete Ar->EliminateRowsCols(ess_tdofs);
  delete Ai->EliminateRowsCols(ess_tdofs);
  Ar->HostReadWrite();
  auto *A = (hypre_ParCSRMatrix *)*Ar;
  for (hypre_CSRMatrix *block : {hypre_ParCSRMatrixDiag(A), hypre_ParCSRMatrixOffd(A)})
  {
    auto *data = hypre_CSRMatrixData(block);
    std::replace(data, data + hypre_CSRMatrixNumNonzeros(block), 0.0, tiny_entry);
  }
  Ar->HypreReadWrite();
  return {std::move(Ar), std::move(Ai)};
}

// Unknowns of the whole boundary of the mesh.
mfem::Array<int> BoundaryTrueDofs(mfem::ParFiniteElementSpace &fespace)
{
  mfem::Array<int> ess_bdr(fespace.GetParMesh()->bdr_attributes.Max()), ess_tdofs;
  ess_bdr = 1;
  fespace.GetEssentialTrueDofs(ess_bdr, ess_tdofs);
  return ess_tdofs;
}

}  // namespace

TEST_CASE("ComplexLocalCSR drops small entries",
          "[strumpack][superlu][mumps][Serial][Parallel]")
{
  // The complex-valued matrix which the solvers in complex arithmetic factor with
  // DropSmallEntries. For a real-valued matrix, it has the entries which
  // mfem::HypreParMatrix::DropSmallEntries leaves (what MfemWrapperSolver factors). For a
  // complex-valued matrix whose imaginary part has the sparsity pattern of the real part,
  // the same entries are dropped: the eliminated ones, whose imaginary part is zero.
  MPI_Comm comm = Mpi::World();
  mfem::Mesh serial_mesh =
      mfem::Mesh::MakeCartesian3D(6, 6, 6, mfem::Element::HEXAHEDRON, 1.0, 1.0, 1.0);
  mfem::ParMesh mesh(comm, serial_mesh);
  mfem::H1_FECollection fec(2, 3);
  mfem::ParFiniteElementSpace fespace(&mesh, &fec);
  const auto ess_tdofs = BoundaryTrueDofs(fespace);
  const auto [Ar, Ai] = AssembleEliminatedTestMatrix(fespace, ess_tdofs);
  auto GlobalNNZ = [&](const ComplexLocalCSR &csr)
  {
    HYPRE_BigInt nnz = csr.J.size();
    Mpi::GlobalSum(1, &nnz, comm);
    return nnz;
  };

  // Real part alone, also with a small diagonal entry, which is kept, in the first row of
  // every process which is not eliminated.
  {
    mfem::HypreParMatrix Ar_diag(*Ar);
    Ar_diag.HostReadWrite();
    hypre_CSRMatrix *diag = hypre_ParCSRMatrixDiag((hypre_ParCSRMatrix *)Ar_diag);
    for (int i = 0; i < Ar_diag.Height(); i++)
    {
      if (ess_tdofs.Find(i) < 0)
      {
        for (HYPRE_Int k = hypre_CSRMatrixI(diag)[i]; k < hypre_CSRMatrixI(diag)[i + 1];
             k++)
        {
          if (hypre_CSRMatrixJ(diag)[k] == i)
          {
            hypre_CSRMatrixData(diag)[k] = tiny_entry;
          }
        }
        break;
      }
    }
    Ar_diag.HypreReadWrite();
    mfem::HypreParMatrix Ar_dropped(Ar_diag);
    Ar_dropped.DropSmallEntries(drop_small_entries_tol);
    REQUIRE(Ar_dropped.NNZ() < Ar_diag.NNZ());
    const ComplexLocalCSR ref(ComplexWrapperOperator(&Ar_dropped, nullptr));
    const ComplexLocalCSR all(ComplexWrapperOperator(&Ar_diag, nullptr));
    const ComplexLocalCSR csr(ComplexWrapperOperator(&Ar_diag, nullptr), true);
    CHECK(GlobalNNZ(all) == Ar_diag.NNZ());
    CHECK(GlobalNNZ(csr) == Ar_dropped.NNZ());
    bool same = (csr.I == ref.I && csr.J == ref.J && csr.V == ref.V);
    Mpi::GlobalAnd(1, &same, comm);
    CHECK(same);
  }

  // Real and imaginary parts: same sparsity pattern as for the real part alone, and every
  // eliminated row is left with its diagonal entry.
  const ComplexLocalCSR csr_real(ComplexWrapperOperator(Ar.get(), nullptr), true);
  const ComplexLocalCSR csr(ComplexWrapperOperator(Ar.get(), Ai.get()), true);
  bool same = (csr.I == csr_real.I && csr.J == csr_real.J);
  for (const int i : ess_tdofs)
  {
    same = same && (csr.I[i + 1] - csr.I[i] == 1) &&
           (csr.J[csr.I[i]] == csr.first_row + i) &&
           (csr.V[csr.I[i]] == std::complex<double>(1.0, 1.0));
  }
  Mpi::GlobalAnd(1, &same, comm);
  CHECK(same);
}

TEST_CASE("Complex arithmetic solvers drop small entries",
          "[strumpack][superlu][mumps][Serial][Parallel]")
{
  // With DropSmallEntries, each solver in complex arithmetic factors the matrix without
  // its small entries, for the same solution as without it: the solution is exactly zero
  // for the eliminated unknowns, which the tiny entries would otherwise couple to the
  // others. STRUMPACK itself ignores the values of entries below the machine precision,
  // but not their positions: its factors have fewer entries without them. The sparsity
  // pattern of the factored matrix then changes with the values of the matrix, also when
  // the reordering is reused.
  MPI_Comm comm = Mpi::World();
  const auto type = GENERATE(from_range(SparseDirectSolverTypes()));
  INFO("type " << type);
  auto MakeSolver = [&](bool drop_small_entries) -> std::unique_ptr<Solver<ComplexOperator>>
  {
#if defined(MFEM_USE_STRUMPACK)
    if (type == "STRUMPACK")
    {
      return std::make_unique<StrumpackComplexSolver>(comm, SymbolicFactorization::DEFAULT,
                                                      SparseCompression::NONE, 0.0, 0, 0,
                                                      drop_small_entries, true, 0);
    }
#endif
#if defined(MFEM_USE_SUPERLU)
    if (type == "SuperLU")
    {
      return std::make_unique<SuperLUComplexSolver>(comm, SymbolicFactorization::DEFAULT,
                                                    false, drop_small_entries, true, 0);
    }
#endif
#if defined(MFEM_USE_MUMPS)
    if (type == "MUMPS")
    {
      return std::make_unique<MumpsComplexSolver>(comm, MatrixSymmetry::SYMMETRIC,
                                                  SymbolicFactorization::DEFAULT, 0.0,
                                                  drop_small_entries, true, 0);
    }
#endif
    return nullptr;
  };

  mfem::Mesh serial_mesh =
      mfem::Mesh::MakeCartesian3D(6, 6, 6, mfem::Element::HEXAHEDRON, 1.0, 1.0, 1.0);
  mfem::ParMesh mesh(comm, serial_mesh);
  mfem::H1_FECollection fec(2, 3);
  mfem::ParFiniteElementSpace fespace(&mesh, &fec);
  const auto ess_tdofs = BoundaryTrueDofs(fespace);
  const auto [Ar, Ai] = AssembleEliminatedTestMatrix(fespace, ess_tdofs);
  const ComplexWrapperOperator A(Ar.get(), Ai.get());
  const auto A_full = AssembleTestMatrix(fespace, 1.0);

  const int n = A.Height();
  ComplexVector b(n), x(n), x_all(n), r(n);
  b.Real().Randomize(1);
  b.Imag().Randomize(2);
  for (const int i : ess_tdofs)
  {
    b.Real()[i] = b.Imag()[i] = 0.0;
  }
  const double norm_b = linalg::Norml2(comm, b);
  constexpr double tol = 1.0e-10;
  auto Residual = [&](const ComplexOperator &M)
  {
    M.Mult(x, r);
    r -= b;
    return linalg::Norml2(comm, r) / norm_b;
  };
  auto MaxEliminated = [&](const ComplexVector &v)
  {
    double max = 0.0;
    for (const int i : ess_tdofs)
    {
      max = std::max(max, std::hypot(v.Real()[i], v.Imag()[i]));
    }
    Mpi::GlobalMax(1, &max, comm);
    return max;
  };

  // Without DropSmallEntries, the tiny entries are factored.
  std::size_t factor_nnz_all = 0;
  {
    auto solver = MakeSolver(false);
    solver->SetOperator(A);
    solver->Mult(b, x_all);
#if defined(MFEM_USE_STRUMPACK)
    if (auto *strumpack = dynamic_cast<StrumpackComplexSolver *>(solver.get()))
    {
      factor_nnz_all = strumpack->GetSolver().factor_nonzeros();
    }
    else
#endif
    {
      CHECK(MaxEliminated(x_all) > 0.0);
    }
  }

  // With DropSmallEntries: the matrix, then a matrix with the same stored entries of which
  // none is small, and the first one again.
  auto solver = MakeSolver(true);
  for (const ComplexOperator *M : {static_cast<const ComplexOperator *>(&A),
                                   static_cast<const ComplexOperator *>(A_full.get()),
                                   static_cast<const ComplexOperator *>(&A)})
  {
    solver->SetOperator(*M);
    x = 0.0;
    solver->Mult(b, x);
    CHECK(Residual(*M) < tol);
    if (M == &A)
    {
      CHECK(MaxEliminated(x) == 0.0);
#if defined(MFEM_USE_STRUMPACK)
      if (auto *strumpack = dynamic_cast<StrumpackComplexSolver *>(solver.get()))
      {
        CHECK(strumpack->GetSolver().factor_nonzeros() < factor_nnz_all);
      }
#endif
      x -= x_all;
      CHECK(linalg::Norml2(comm, x) / linalg::Norml2(comm, x_all) < tol);
    }
  }
}

}  // namespace palace

#endif
