// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include <cmath>
#include <complex>
#include <memory>
#include <type_traits>
#include <vector>
#include <mfem.hpp>
#include <catch2/catch_test_macros.hpp>
#include "fem/mesh.hpp"
#include "linalg/distrelaxation.hpp"
#include "linalg/gmg.hpp"
#include "linalg/ksp.hpp"
#include "linalg/operator.hpp"
#include "linalg/solver.hpp"
#include "linalg/vector.hpp"
#include "models/laplaceoperator.hpp"
#include "models/spaceoperator.hpp"
#include "utils/communication.hpp"
#include "utils/iodata.hpp"
#include "utils/units.hpp"

namespace palace
{

namespace
{

// Linear solver whose preconditioner can be applied on its own.
template <OperatorType OperType>
class InspectableKspSolver : public BaseKspSolver<OperType>
{
public:
  using BaseKspSolver<OperType>::BaseKspSolver;

  Solver<OperType> &GetPreconditioner() { return *this->pc; }
};

// An electrostatic problem on the unit cube at order 3, so with p-multigrid, solved with
// GMRES.
struct ElectrostaticProblem
{
  Units units;
  IoData iodata;
  std::vector<std::unique_ptr<Mesh>> mesh;
  std::unique_ptr<LaplaceOperator> laplace_op;
  std::unique_ptr<Operator> K;

  ElectrostaticProblem(MultigridCycle cycle, bool cheby_4th = true)
    : units(1.0, 1.0), iodata(units)
  {
    auto smesh = std::make_unique<mfem::Mesh>(
        mfem::Mesh::MakeCartesian3D(4, 4, 4, mfem::Element::HEXAHEDRON, 1.0, 1.0, 1.0));
    iodata.model.Lc = 1.0;
    iodata.problem.type = ProblemType::ELECTROSTATIC;
    iodata.domains.materials.emplace_back().attributes = {1};
    for (int i = 0; i < smesh->bdr_attributes.Size(); i++)
    {
      iodata.boundaries.pec.attributes.push_back(smesh->bdr_attributes[i]);
    }
    iodata.solver.order = 3;
    iodata.solver.linear.krylov_solver = KrylovSolver::GMRES;
    iodata.solver.linear.mg_cycle = cycle;
    iodata.solver.linear.mg_smooth_cheby_4th = cheby_4th;
    iodata.solver.linear.tol = 1.0e-10;
    iodata.NondimensionalizeInputs(smesh);
    iodata.CheckConfiguration();
    mesh.push_back(
        std::make_unique<Mesh>(std::make_unique<mfem::ParMesh>(Mpi::World(), *smesh)));
    laplace_op = std::make_unique<LaplaceOperator>(iodata, mesh);
    K = laplace_op->GetStiffnessMatrix();
  }
};

// The shifted operator K - σ M of an eigenmode problem in a PEC cube at order 2, so with
// p-multigrid and distributive relaxation smoothing, below the first eigenvalue 2π².
struct ShiftedMaxwellProblem
{
  Units units;
  IoData iodata;
  std::vector<std::unique_ptr<Mesh>> mesh;
  std::unique_ptr<SpaceOperator> space_op;
  std::unique_ptr<ComplexOperator> K, C, M, A, P;

  ShiftedMaxwellProblem(MultigridCycle cycle, int smooth_order = -1, int smooth_it = 1)
    : units(1.0, 1.0), iodata(units)
  {
    using namespace std::complex_literals;
    auto smesh = std::make_unique<mfem::Mesh>(
        mfem::Mesh::MakeCartesian3D(4, 4, 4, mfem::Element::HEXAHEDRON, 1.0, 1.0, 1.0));
    iodata.model.Lc = 1.0;
    iodata.problem.type = ProblemType::EIGENMODE;
    iodata.domains.materials.emplace_back().attributes = {1};
    for (int i = 0; i < smesh->bdr_attributes.Size(); i++)
    {
      iodata.boundaries.pec.attributes.push_back(smesh->bdr_attributes[i]);
    }
    iodata.solver.order = 2;
    iodata.solver.linear.krylov_solver = KrylovSolver::GMRES;
    iodata.solver.linear.mg_cycle = cycle;
    iodata.solver.linear.mg_smooth_order = smooth_order;
    iodata.solver.linear.mg_smooth_it = smooth_it;
    iodata.solver.linear.tol = 1.0e-10;
    iodata.NondimensionalizeInputs(smesh);
    iodata.CheckConfiguration();
    mesh.push_back(
        std::make_unique<Mesh>(std::make_unique<mfem::ParMesh>(Mpi::World(), *smesh)));
    space_op = std::make_unique<SpaceOperator>(iodata, mesh);
    constexpr double target = 3.0;
    K = space_op->GetStiffnessMatrix<ComplexOperator>(Operator::DIAG_ONE);
    C = space_op->GetDampingMatrix<ComplexOperator>(Operator::DIAG_ZERO);
    M = space_op->GetMassMatrix<ComplexOperator>(Operator::DIAG_ZERO);
    A = space_op->GetSystemMatrix(1.0 + 0.0i, 1i * target, -target * target + 0.0i, K.get(),
                                  C.get(), M.get());
    P = space_op->GetPreconditionerMatrix<ComplexOperator>(
        1.0 + 0.0i, 1i * target, -target * target + 0.0i, target + 0.0i);
  }
};

// Relative difference of two vectors.
template <typename VecType>
double RelativeDifference(MPI_Comm comm, const VecType &x, const VecType &y)
{
  VecType d(x);
  d -= y;
  return linalg::Norml2(comm, d) / linalg::Norml2(comm, y);
}

}  // namespace

TEST_CASE("Geometric multigrid one-sided V-cycle smoothing orders", "[gmg][Serial]")
{
  // With the default smoothing order of symmetric V-cycles, max(2 p, 4) for the finest order
  // p, one-sided V-cycles smooth each level at twice its own default order, max(4 p_l, 8).
  using V = std::vector<int>;
  CHECK(GetOneSidedSmoothOrders({1, 2, 4}, 8) == V{8, 8, 16});
  CHECK(GetOneSidedSmoothOrders({1, 2, 3}, 6) == V{8, 8, 12});
  CHECK(GetOneSidedSmoothOrders({1, 1, 2}, 4) == V{8, 8, 8});
  CHECK(GetOneSidedSmoothOrders({1, 1}, 4) == V{8, 8});

  // With another smoothing order, twice it on every level.
  CHECK(GetOneSidedSmoothOrders({1, 2, 4}, 6) == V{12, 12, 12});
  CHECK(GetOneSidedSmoothOrders({1, 2, 4}, 16) == V{32, 32, 32});
}

TEST_CASE("Geometric multigrid one-sided V-cycles", "[gmg][Serial][Parallel][GPU]")
{
  // One-sided V-cycles, with twice the smoothing order, make a nonsymmetric preconditioner
  // for the electrostatic problem, with which GMRES converges as with symmetric V-cycles.
  MPI_Comm comm = Mpi::World();
  for (bool one_sided : {false, true})
  {
    CAPTURE(one_sided);
    ElectrostaticProblem problem(one_sided ? MultigridCycle::ONE_SIDED
                                           : MultigridCycle::SYMMETRIC);
    REQUIRE(problem.iodata.solver.linear.mg_smooth_order == 6);
    REQUIRE(problem.laplace_op->GetH1Spaces().GetNumLevels() > 1);
    const auto &K = *problem.K;
    InspectableKspSolver<Operator> ksp(problem.iodata, problem.laplace_op->GetH1Spaces());
    ksp.SetOperators(K, K);

    // (y, B x) = (x, B y) only for symmetric V-cycles. The strong smoothers of this small
    // problem leave little asymmetry to the one-sided ones.
    const int n = K.Height();
    Vector x(n), y(n), Bx(n), By(n), r(n);
    for (auto *v : {&x, &y, &Bx, &By, &r})
    {
      v->UseDevice(true);
    }
    linalg::SetRandom(comm, x, 1);
    linalg::SetRandom(comm, y, 2);
    const auto &pc = ksp.GetPreconditioner();
    pc.Mult(x, Bx);
    pc.Mult(y, By);
    const double asymmetry = std::abs(linalg::Dot(comm, y, Bx) - linalg::Dot(comm, x, By)) /
                             (linalg::Norml2(comm, y) * linalg::Norml2(comm, Bx));
    if (one_sided)
    {
      CHECK(asymmetry > 1.0e-9);
    }
    else
    {
      CHECK(asymmetry < 1.0e-12);
    }

    // GMRES converges to the solution of K u = x.
    Vector u(n);
    u.UseDevice(true);
    u = 0.0;
    ksp.Mult(x, u);
    CHECK(ksp.GetConverged());
    K.Mult(u, r);
    r -= x;
    CHECK(linalg::Norml2(comm, r) <= 1.0e-8 * linalg::Norml2(comm, x));
  }
}

TEST_CASE("Geometric multigrid one-sided V-cycles with distributive relaxation",
          "[gmg][Serial][Parallel][GPU]")
{
  // The smoother of the shifted Maxwell problem splits into the auxiliary space correction
  // of its last iteration and the rest, which one-sided V-cycles apply after and before the
  // coarse-grid correction. GMRES converges with one-sided V-cycles with one or two
  // smoothing iterations.
  MPI_Comm comm = Mpi::World();
  for (int smooth_it : {1, 2})
  {
    CAPTURE(smooth_it);
    ShiftedMaxwellProblem problem(MultigridCycle::ONE_SIDED, -1, smooth_it);
    REQUIRE(problem.iodata.solver.linear.mg_smooth_aux == 1);
    const auto &space_op = *problem.space_op;
    REQUIRE(space_op.GetNDSpaces().GetNumLevels() > 1);
    const auto &A = *problem.A;

    // Mult2 is MultExceptLastAuxiliary2 followed by MultAuxiliary2.
    const auto &mg_op = dynamic_cast<const ComplexMultigridOperator &>(*problem.P);
    DistRelaxationSmoother<ComplexOperator> B(
        comm, space_op.GetGradMatrix(), smooth_it, 1,
        2 * problem.iodata.solver.linear.mg_smooth_order, 1.0, 0.0, true);
    B.SetOperators(mg_op.GetFinestOperator(), mg_op.GetFinestAuxiliaryOperator());
    B.SetInitialGuess(false);
    const int n = A.Height();
    ComplexVector x(n), y(n), z(n), r(n);
    for (auto *v : {&x, &y, &z, &r})
    {
      v->UseDevice(true);
    }
    linalg::SetRandom(comm, x, 1);
    B.Mult2(x, y, r);
    B.MultExceptLastAuxiliary2(x, z, r);
    z -= y;
    CHECK(linalg::Norml2(comm, z) > 1.0e-3 * linalg::Norml2(comm, y));
    z += y;
    B.MultAuxiliary2(x, z, r);
    z -= y;
    CHECK(linalg::Norml2(comm, z) <= 1.0e-12 * linalg::Norml2(comm, y));

    // GMRES converges to the solution of A u = x, in 9 iterations with one smoothing
    // iteration (50 without the auxiliary space correction).
    ComplexKspSolver ksp(problem.iodata, problem.space_op->GetNDSpaces(),
                         &problem.space_op->GetH1Spaces());
    ksp.SetOperators(A, *problem.P);
    ComplexVector u(n);
    u.UseDevice(true);
    u = 0.0;
    ksp.Mult(x, u);
    CHECK(ksp.GetConverged());
    CHECK(ksp.NumTotalMultIterations() <= 20);
    A.Mult(u, r);
    r -= x;
    CHECK(linalg::Norml2(comm, r) <= 1.0e-8 * linalg::Norml2(comm, x));
  }
}

TEST_CASE("Geometric multigrid V-cycle switch", "[gmg][Serial][Parallel][GPU]")
{
  // Switching between one-sided and symmetric V-cycles gives the preconditioner constructed
  // with the other V-cycles: for the 4th-kind and 1st-kind Chebyshev smoothers of an
  // electrostatic problem, and for the distributive relaxation of a shifted Maxwell
  // problem.
  MPI_Comm comm = Mpi::World();
  auto CheckSwitch = [comm](auto &one_sided, auto &symmetric, const auto &x)
  {
    using VecType = std::decay_t<decltype(x)>;
    VecType y_one_sided(x.Size()), y_symmetric(x.Size()), y(x.Size());
    for (auto *v : {&y_one_sided, &y_symmetric, &y})
    {
      v->UseDevice(true);
    }
    REQUIRE(one_sided.IsOneSided());
    REQUIRE_FALSE(symmetric.IsOneSided());
    one_sided.Mult(x, y_one_sided);
    symmetric.Mult(x, y_symmetric);
    CHECK(RelativeDifference(comm, y_one_sided, y_symmetric) > 1.0e-6);
    symmetric.SetOneSided(true);
    symmetric.Mult(x, y);
    CHECK(RelativeDifference(comm, y, y_one_sided) <= 1.0e-12);
    symmetric.SetOneSided(false);
    symmetric.Mult(x, y);
    CHECK(RelativeDifference(comm, y, y_symmetric) <= 1.0e-12);
    one_sided.SetOneSided(false);
    one_sided.Mult(x, y);
    CHECK(RelativeDifference(comm, y, y_symmetric) <= 1.0e-12);
  };
  for (bool cheby_4th : {true, false})
  {
    CAPTURE(cheby_4th);
    ElectrostaticProblem problem(MultigridCycle::ONE_SIDED, cheby_4th);
    auto &spaces = problem.laplace_op->GetH1Spaces();
    InspectableKspSolver<Operator> ksp_one_sided(problem.iodata, spaces);
    problem.iodata.solver.linear.mg_cycle = MultigridCycle::SYMMETRIC;
    InspectableKspSolver<Operator> ksp_symmetric(problem.iodata, spaces);
    ksp_one_sided.SetOperators(*problem.K, *problem.K);
    ksp_symmetric.SetOperators(*problem.K, *problem.K);
    Vector x(problem.K->Height());
    x.UseDevice(true);
    linalg::SetRandom(comm, x, 1);
    CheckSwitch(dynamic_cast<GeometricMultigridSolver<Operator> &>(
                    ksp_one_sided.GetPreconditioner()),
                dynamic_cast<GeometricMultigridSolver<Operator> &>(
                    ksp_symmetric.GetPreconditioner()),
                x);
  }
  {
    ShiftedMaxwellProblem problem(MultigridCycle::ONE_SIDED);
    auto &nd_spaces = problem.space_op->GetNDSpaces();
    auto &h1_spaces = problem.space_op->GetH1Spaces();
    InspectableKspSolver<ComplexOperator> ksp_one_sided(problem.iodata, nd_spaces,
                                                        &h1_spaces);
    problem.iodata.solver.linear.mg_cycle = MultigridCycle::SYMMETRIC;
    InspectableKspSolver<ComplexOperator> ksp_symmetric(problem.iodata, nd_spaces,
                                                        &h1_spaces);
    ksp_one_sided.SetOperators(*problem.A, *problem.P);
    ksp_symmetric.SetOperators(*problem.A, *problem.P);
    ComplexVector x(problem.A->Height());
    x.UseDevice(true);
    linalg::SetRandom(comm, x, 1);
    CheckSwitch(dynamic_cast<GeometricMultigridSolver<ComplexOperator> &>(
                    ksp_one_sided.GetPreconditioner()),
                dynamic_cast<GeometricMultigridSolver<ComplexOperator> &>(
                    ksp_symmetric.GetPreconditioner()),
                x);
  }
}

TEST_CASE("Geometric multigrid V-cycles chosen on the first solve",
          "[gmg][Serial][Parallel]")
{
  // With "Auto", the first solve keeps the one-sided V-cycles if they need at most as many
  // GMRES iterations as symmetric ones, otherwise the symmetric ones, with their solution,
  // and the next solves use them. On the shifted Maxwell problem, symmetric V-cycles need
  // fewer iterations with smoothing order 2, and one-sided ones with order 1.
  MPI_Comm comm = Mpi::World();
  for (int smooth_order : {1, 2})
  {
    CAPTURE(smooth_order);
    ShiftedMaxwellProblem problem(MultigridCycle::AUTO, smooth_order);
    auto &nd_spaces = problem.space_op->GetNDSpaces();
    auto &h1_spaces = problem.space_op->GetH1Spaces();
    const int n = problem.A->Height();
    ComplexVector x1(n), x2(n);
    x1.UseDevice(true);
    x2.UseDevice(true);
    linalg::SetRandom(comm, x1, 1);
    linalg::SetRandom(comm, x2, 2);
    auto Solve = [&](InspectableKspSolver<ComplexOperator> &ksp, const ComplexVector &x,
                     ComplexVector &u)
    {
      u.SetSize(n);
      u.UseDevice(true);
      u = 0.0;
      const int its = ksp.NumTotalMultIterations();
      ksp.Mult(x, u);
      CHECK(ksp.GetConverged());
      return ksp.NumTotalMultIterations() - its;
    };

    // The same solves with each V-cycle.
    ComplexVector u_one_sided, u_symmetric, u, u2;
    problem.iodata.solver.linear.mg_cycle = MultigridCycle::ONE_SIDED;
    InspectableKspSolver<ComplexOperator> ksp_one_sided(problem.iodata, nd_spaces,
                                                        &h1_spaces);
    problem.iodata.solver.linear.mg_cycle = MultigridCycle::SYMMETRIC;
    InspectableKspSolver<ComplexOperator> ksp_symmetric(problem.iodata, nd_spaces,
                                                        &h1_spaces);
    problem.iodata.solver.linear.mg_cycle = MultigridCycle::AUTO;
    InspectableKspSolver<ComplexOperator> ksp(problem.iodata, nd_spaces, &h1_spaces);
    for (auto *ksp_c : {&ksp_one_sided, &ksp_symmetric, &ksp})
    {
      ksp_c->SetOperators(*problem.A, *problem.P);
    }
    const int its_one_sided = Solve(ksp_one_sided, x1, u_one_sided);
    const int its_symmetric = Solve(ksp_symmetric, x1, u_symmetric);
    const bool one_sided = (its_one_sided <= its_symmetric);
    CHECK(one_sided == (smooth_order == 1));

    // The first solve tries both, the second uses the kept V-cycles.
    CHECK(Solve(ksp, x1, u) == its_symmetric + (one_sided ? its_one_sided : its_symmetric));
    CHECK(RelativeDifference(comm, u, one_sided ? u_one_sided : u_symmetric) <= 1.0e-12);
    const auto &gmg =
        dynamic_cast<GeometricMultigridSolver<ComplexOperator> &>(ksp.GetPreconditioner());
    CHECK(gmg.IsOneSided() == one_sided);
    CHECK(Solve(ksp, x2, u) == Solve(one_sided ? ksp_one_sided : ksp_symmetric, x2, u2));
    CHECK(RelativeDifference(comm, u, u2) <= 1.0e-12);
  }
}

}  // namespace palace
