// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include <cmath>
#include <memory>
#include <vector>
#include <mfem.hpp>
#include <catch2/catch_test_macros.hpp>
#include "fem/mesh.hpp"
#include "linalg/ksp.hpp"
#include "linalg/operator.hpp"
#include "linalg/solver.hpp"
#include "linalg/vector.hpp"
#include "models/laplaceoperator.hpp"
#include "utils/communication.hpp"
#include "utils/iodata.hpp"
#include "utils/units.hpp"

namespace palace
{

namespace
{

// Linear solver whose preconditioner can be applied on its own.
class InspectableKspSolver : public KspSolver
{
public:
  using KspSolver::KspSolver;

  const Solver<Operator> &GetPreconditioner() const { return *pc; }
};

}  // namespace

TEST_CASE("Geometric multigrid one-sided V-cycles", "[gmg][Serial][Parallel][GPU]")
{
  // An electrostatic problem on the unit cube at order 3, so with p-multigrid, solved with
  // GMRES. One-sided V-cycles, with twice the default smoothing order, make a nonsymmetric
  // preconditioner, with which GMRES converges as with symmetric V-cycles.
  MPI_Comm comm = Mpi::World();
  for (bool one_sided : {false, true})
  {
    CAPTURE(one_sided);
    auto smesh = std::make_unique<mfem::Mesh>(
        mfem::Mesh::MakeCartesian3D(4, 4, 4, mfem::Element::HEXAHEDRON, 1.0, 1.0, 1.0));
    Units units(1.0, 1.0);
    IoData iodata(units);
    iodata.model.Lc = 1.0;
    iodata.problem.type = ProblemType::ELECTROSTATIC;
    iodata.domains.materials.emplace_back().attributes = {1};
    for (int i = 0; i < smesh->bdr_attributes.Size(); i++)
    {
      iodata.boundaries.pec.attributes.push_back(smesh->bdr_attributes[i]);
    }
    iodata.solver.order = 3;
    iodata.solver.linear.krylov_solver = KrylovSolver::GMRES;
    iodata.solver.linear.mg_one_sided = one_sided;
    iodata.solver.linear.tol = 1.0e-10;
    iodata.NondimensionalizeInputs(smesh);
    iodata.CheckConfiguration();
    REQUIRE(iodata.solver.linear.mg_one_sided == one_sided);
    REQUIRE(iodata.solver.linear.mg_smooth_order == (one_sided ? 12 : 6));
    std::vector<std::unique_ptr<Mesh>> mesh;
    mesh.push_back(std::make_unique<Mesh>(std::make_unique<mfem::ParMesh>(comm, *smesh)));
    LaplaceOperator laplace_op(iodata, mesh);
    REQUIRE(laplace_op.GetH1Spaces().GetNumLevels() > 1);
    auto K = laplace_op.GetStiffnessMatrix();
    InspectableKspSolver ksp(iodata, laplace_op.GetH1Spaces());
    ksp.SetOperators(*K, *K);

    // (y, B x) = (x, B y) only for symmetric V-cycles. The strong smoothers of this small
    // problem leave little asymmetry to the one-sided ones.
    const int n = K->Height();
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
    K->Mult(u, r);
    r -= x;
    CHECK(linalg::Norml2(comm, r) <= 1.0e-8 * linalg::Norml2(comm, x));
  }
}

}  // namespace palace
