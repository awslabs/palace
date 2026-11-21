// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include <algorithm>
#include <cmath>
#include <complex>
#include <memory>
#include <numeric>
#include <set>
#include <vector>
#include <mfem.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators_all.hpp>
#include "fem/fieldtransfer.hpp"
#include "linalg/eps.hpp"
#include "linalg/iterative.hpp"
#include "linalg/jacobi.hpp"
#include "linalg/ksp.hpp"
#include "linalg/operator.hpp"
#include "linalg/slepc.hpp"
#include "linalg/vector.hpp"
#include "utils/communication.hpp"

namespace palace
{

TEST_CASE("Vector field transfer through refinement and rebalancing",
          "[eigenmode][Serial][Parallel][GPU]")
{
  // H(curl) fields of a space of order p, linear in the coordinates, are kept exactly
  // through nonconforming refinement and rebalancing of their mesh.
  const auto type = GENERATE(mfem::Element::TETRAHEDRON, mfem::Element::HEXAHEDRON);
  const int order = GENERATE(1, 2);
  MPI_Comm comm = Mpi::World();
  auto smesh = mfem::Mesh::MakeCartesian3D(3, 3, 3, type, 1.0, 1.0, 1.0);
  smesh.EnsureNCMesh(true);
  mfem::ParMesh pmesh(comm, smesh);
  mfem::ND_FECollection nd_fec(order, pmesh.Dimension());

  // Two fields: a constant one, exact for every order, and a linear one for order 2.
  auto Field = [order](const mfem::Vector &x, int k, mfem::Vector &E)
  {
    E.SetSize(3);
    E(0) = 1.0 + ((order > 1) ? 2.0 * x(1) - x(2) : 0.0);
    E(1) = -0.5 + ((order > 1) ? 3.0 * x(0) + x(2) : 0.0);
    E(2) = 2.0 + ((order > 1) ? x(0) - 4.0 * x(1) : 0.0);
    if (k == 1)
    {
      E *= -0.25;
      E(0) += 1.0;
    }
  };
  auto TrueDofs = [&](mfem::ParFiniteElementSpace &fespace, int k)
  {
    mfem::VectorFunctionCoefficient coeff(
        3, [&Field, k](const mfem::Vector &x, mfem::Vector &E) { Field(x, k, E); });
    mfem::ParGridFunction gf(&fespace);
    gf.ProjectCoefficient(coeff);
    Vector tv(fespace.GetTrueVSize());
    tv.UseDevice(true);
    gf.GetTrueDofs(tv);
    return tv;
  };

  VectorFieldTransfer transfer(pmesh, order);
  {
    mfem::ParFiniteElementSpace nd_fespace(&pmesh, &nd_fec);
    for (int k = 0; k < 2; k++)
    {
      transfer.AddField(nd_fespace, TrueDofs(nd_fespace, k));
    }
  }
  CHECK(transfer.Size() == 2u);
  CHECK(transfer.IsCurrent(pmesh));

  // Refine the elements in a corner, twice, then rebalance.
  for (int it = 0; it < 2; it++)
  {
    mfem::Array<int> marked;
    for (int e = 0; e < pmesh.GetNE(); e++)
    {
      mfem::Vector c;
      pmesh.GetElementCenter(e, c);
      if (c(0) < 0.5 && c(1) < 0.5)
      {
        marked.Append(e);
      }
    }
    pmesh.GeneralRefinement(marked, 1, 0);
    CHECK(!transfer.IsCurrent(pmesh));
    transfer.Update();
    CHECK(transfer.IsCurrent(pmesh));
  }
  pmesh.Rebalance();
  transfer.Update();
  CHECK(transfer.IsCurrent(pmesh));

  mfem::ParFiniteElementSpace nd_fespace(&pmesh, &nd_fec);
  for (int k = 0; k < 2; k++)
  {
    const auto ref = TrueDofs(nd_fespace, k);
    Vector tv;
    tv.UseDevice(true);
    transfer.GetField(k, nd_fespace, tv);
    REQUIRE(tv.Size() == ref.Size());
    tv -= ref;
    const double err = linalg::Norml2(comm, tv), norm = linalg::Norml2(comm, ref);
    CHECK(err <= 1.0e-12 * norm);
  }
}

TEST_CASE("Vector field transfer through conformal refinement, new marks and a move",
          "[eigenmode][Serial][Parallel]")
{
  // H(curl) fields of a space of order p, linear in the coordinates, are kept exactly
  // through a conformal refinement of tetrahedra, a new marking of the mesh, which reorders
  // the vertices of elements, and a parallel move of the elements. Without putting each
  // element's values back through its new vertex order, they are not. Also on triangles,
  // and on tetrahedra with curved (order 2) nodes, where the values at the element centers
  // are compared before and after the marking, as the space does not hold linear fields.
  const int order = GENERATE(1, 2);
  const bool remark = GENERATE(true, false);
  const int kind = GENERATE(0, 1, 2);  // Straight tetrahedra, triangles, curved tetrahedra
  CAPTURE(order, remark, kind);
  const int dim = (kind == 1) ? 2 : 3;
  MPI_Comm comm = Mpi::World();
  const int size = Mpi::Size(comm);
  auto smesh =
      (dim == 3)
          ? mfem::Mesh::MakeCartesian3D(3, 3, 3, mfem::Element::TETRAHEDRON, 1.0, 1.0, 1.0)
          : mfem::Mesh::MakeCartesian2D(4, 4, mfem::Element::TRIANGLE, true, 1.0, 1.0);
  if (kind == 2)
  {
    smesh.SetCurvature(2);
  }
  smesh.Transform(  // Unequal edge lengths, so that new marks reorder vertices
      [dim](const mfem::Vector &x, mfem::Vector &y)
      {
        y = x;
        const double z = (dim == 3) ? x(2) : 0.0;
        // Irregular, so that the new marks give many vertex orders.
        y(0) += 0.08 * std::sin(17.0 * x(1) + 13.0 * z + 5.0 * x(0));
        y(1) += 0.08 * std::sin(15.0 * z + 12.0 * x(0) + 3.0 * x(1));
        if (dim == 3)
        {
          y(2) += 0.08 * std::sin(13.0 * x(0) + 14.0 * x(1) + 7.0 * x(2));
        }
        else
        {
          y(1) *= 1.7;  // Stretched, so that the new marks of the bisected triangles differ
        }
      });
  smesh.Finalize(true);
  std::vector<int> part(smesh.GetNE(), 0);
  for (int e = 0; e < smesh.GetNE(); e++)
  {
    part[e] = (e % 5 == 0) ? (e / 5) % size : 0;
  }
  mfem::ParMesh pmesh(comm, smesh, part.data());
  mfem::ND_FECollection nd_fec(order, pmesh.Dimension());
  auto Field = [order, dim](const mfem::Vector &x, int k, mfem::Vector &E)
  {
    E.SetSize(dim);
    const double z = (dim == 3) ? x(2) : 0.0;
    // A rotation, which the lowest order holds, and at order 2 any linear field.
    E(0) = 1.0 + 0.5 * z - 0.3 * x(1) + ((order > 1) ? 2.0 * x(1) - z : 0.0);
    E(1) = -0.5 + 0.3 * x(0) - 0.2 * z + ((order > 1) ? 3.0 * x(0) + z : 0.0);
    if (dim == 3)
    {
      E(2) = 2.0 + 0.2 * x(1) - 0.5 * x(0) + ((order > 1) ? x(0) - 4.0 * x(1) : 0.0);
    }
    if (k == 1)
    {
      E *= -0.25;
      E(0) += 1.0;
    }
  };
  auto TrueDofs = [&](mfem::ParFiniteElementSpace &fespace, int k)
  {
    mfem::VectorFunctionCoefficient coeff(
        dim, [&Field, k](const mfem::Vector &x, mfem::Vector &E) { Field(x, k, E); });
    mfem::ParGridFunction gf(&fespace);
    gf.ProjectCoefficient(coeff);
    Vector tv(fespace.GetTrueVSize());
    gf.GetTrueDofs(tv);
    return tv;
  };
  VectorFieldTransfer transfer(pmesh, order);
  {
    mfem::ParFiniteElementSpace nd_fespace(&pmesh, &nd_fec);
    for (int k = 0; k < 2; k++)
    {
      transfer.AddField(nd_fespace, TrueDofs(nd_fespace, k));
    }
  }

  // Refine the elements in a corner, conformally; mark the mesh again; move the elements to
  // a balanced partition in the order of their global numbers.
  mfem::Array<int> marked;
  for (int e = 0; e < pmesh.GetNE(); e++)
  {
    mfem::Vector c;
    pmesh.GetElementCenter(e, c);
    if (c(0) < 0.5 && c(1) < 0.5)
    {
      marked.Append(e);
    }
  }
  pmesh.GeneralRefinement(marked, 0, 0);
  REQUIRE(pmesh.Conforming());
  transfer.Update();
  auto VertexLists = [&pmesh]()
  {
    std::vector<int> v;
    for (int e = 0; e < pmesh.GetNE(); e++)
    {
      const auto *el = pmesh.GetElement(e);
      v.insert(v.end(), el->GetVertices(), el->GetVertices() + el->GetNVertices());
    }
    return v;
  };
  // The values of the fields at the element centers, which a marking does not move.
  auto CenterValues = [&]()
  {
    mfem::ParFiniteElementSpace nd_fespace(&pmesh, &nd_fec);
    std::vector<double> values;
    for (int k = 0; k < 2; k++)
    {
      Vector tv;
      transfer.GetField(k, nd_fespace, tv);
      mfem::ParGridFunction gf(&nd_fespace);
      gf.SetFromTrueDofs(tv);
      gf.HostRead();
      mfem::Vector v;
      for (int e = 0; e < pmesh.GetNE(); e++)
      {
        gf.GetVectorValue(e, mfem::Geometries.GetCenter(pmesh.GetElementBaseGeometry(e)),
                          v);
        values.insert(values.end(), v.begin(), v.end());
      }
    }
    return values;
  };
  // The L2 fields of the transfer at a point of each element given by barycentric
  // coordinates on its vertices in the order of `verts`, which is the same physical point
  // before and after the marking (the weights follow the vertices, not their order).
  const double bary[4] = {0.1, 0.2, 0.3, 0.4};
  auto L2Values = [&](const std::vector<int> &verts, const std::vector<int> &ref)
  {
    const int nv = dim + 1;
    std::vector<double> values;
    for (int k = 0; k < 2; k++)
    {
      const auto &gf = transfer.GetL2Field(k);
      gf.HostRead();
      mfem::Vector v;
      for (int e = 0; e < pmesh.GetNE(); e++)
      {
        double w[4] = {0.0, 0.0, 0.0, 0.0};
        for (int j = 0; j < nv; j++)
        {
          for (int i = 0; i < nv; i++)
          {
            w[j] += (verts[nv * e + j] == ref[nv * e + i])
                        ? bary[i] / (dim == 3 ? 1.0 : 0.6)
                        : 0.0;
          }
        }
        mfem::IntegrationPoint ip;
        ip.Set3(w[1], w[2], (dim == 3) ? w[3] : 0.0);
        gf.GetVectorValue(e, ip, v);
        values.insert(values.end(), v.begin(), v.end());
      }
    }
    return values;
  };
  const auto old_vertices = VertexLists();
  const auto old_values = CenterValues();
  const auto old_l2 = L2Values(old_vertices, old_vertices);
  pmesh.Finalize(true);
  const auto new_vertices = VertexLists();
  int reordered = 0;
  for (std::size_t i = 0; i < old_vertices.size(); i++)
  {
    reordered += (old_vertices[i] != new_vertices[i]);
  }
  Mpi::GlobalSum(1, &reordered, comm);
  CHECK(reordered > 0);
  // The orders of the vertices which the marking gave the elements.
  {
    const int nv = dim + 1;
    std::vector<int> seen(256, 0);  // Over all processes, by the code of the permutation
    for (int e = 0; e < pmesh.GetNE(); e++)
    {
      int code = 0;
      for (int j = 0; j < nv; j++)
      {
        code = 4 * code + static_cast<int>(std::find(old_vertices.begin() + nv * e,
                                                     old_vertices.begin() + nv * (e + 1),
                                                     new_vertices[nv * e + j]) -
                                           (old_vertices.begin() + nv * e));
      }
      seen[code] = 1;
    }
    Mpi::GlobalMax(static_cast<int>(seen.size()), seen.data(), comm);
    const int orders = std::accumulate(seen.begin(), seen.end(), 0);
    CHECK(orders == ((dim == 3) ? 12 : 3));  // Every order a marking can give
  }
  transfer.Remark(remark ? old_vertices : new_vertices, 1);
  CHECK(transfer.IsCurrent(pmesh, 1));
  CHECK_FALSE(transfer.IsCurrent(pmesh, 0));
  {
    const auto new_l2 = L2Values(new_vertices, old_vertices);
    double diff = 0.0, norm = 0.0;
    for (std::size_t i = 0; i < new_l2.size(); i++)
    {
      diff = std::max(diff, std::abs(new_l2[i] - old_l2[i]));
      norm = std::max(norm, std::abs(old_l2[i]));
    }
    Mpi::GlobalMax(1, &diff, comm);
    Mpi::GlobalMax(1, &norm, comm);
    CHECK((remark ? (diff <= 1.0e-12 * norm) : (diff > 1.0e-3 * norm)));
  }
  {
    const auto new_values = CenterValues();
    REQUIRE(new_values.size() == old_values.size());
    double diff = 0.0, norm = 0.0;
    for (std::size_t i = 0; i < new_values.size(); i++)
    {
      diff = std::max(diff, std::abs(new_values[i] - old_values[i]));
      norm = std::max(norm, std::abs(old_values[i]));
    }
    Mpi::GlobalMax(1, &diff, comm);
    Mpi::GlobalMax(1, &norm, comm);
    if (remark)
    {
      CHECK(diff <= 1.0e-12 * norm);
    }
    else if (order > 1 && dim == 3)  // At the center of a triangle, the values do not tell
    {
      CHECK(diff > 1.0e-3 * norm);
    }
  }
  {
    const long long ne = pmesh.GetGlobalNE(), offset = pmesh.GetGlobalElementNum(0);
    mfem::Array<int> partition(pmesh.GetNE());
    for (int e = 0; e < pmesh.GetNE(); e++)
    {
      partition[e] = static_cast<int>(((offset + e) * size) / ne);
    }
    pmesh.Rebalance(partition);
  }
  transfer.Update();
  CHECK(transfer.IsCurrent(pmesh));

  mfem::ParFiniteElementSpace nd_fespace(&pmesh, &nd_fec);
  for (int k = 0; k < 2; k++)
  {
    const auto ref = TrueDofs(nd_fespace, k);
    Vector tv;
    transfer.GetField(k, nd_fespace, tv);
    REQUIRE(tv.Size() == ref.Size());
    tv -= ref;
    const double err = linalg::Norml2(comm, tv), norm = linalg::Norml2(comm, ref);
    if (remark && kind != 2)  // The curved space does not hold the linear fields
    {
      CHECK(err <= 1.0e-12 * norm);
    }
    else if (!remark)
    {
      CHECK(err > 1.0e-3 * norm);  // The new marks move the values of the linear fields
    }
  }
}

#if defined(PALACE_WITH_SLEPC)

namespace
{

// The symmetric tridiagonal stencil [-1, 2 + i / n, -1] minus s times the identity (real,
// as a complex operator) on a single process: distinct eigenvalues in (0, 5).
class TridiagonalEigenTestOperator : public ComplexOperator
{
private:
  const double s;

public:
  TridiagonalEigenTestOperator(int n, double s = 0.0) : ComplexOperator(n), s(s) {}

  bool IsReal() const override { return true; }

  void AssembleDiagonal(ComplexVector &diag) const override
  {
    diag.SetSize(height);
    diag = 0.0;
    auto *D = diag.Real().HostWrite();
    for (int i = 0; i < height; i++)
    {
      D[i] = 2.0 + static_cast<double>(i) / height - s;
    }
  }

  void Mult(const ComplexVector &x, ComplexVector &y) const override
  {
    y.SetSize(height);
    const auto *XR = x.Real().HostRead();
    const auto *XI = x.Imag().HostRead();
    auto *YR = y.Real().HostWrite();
    auto *YI = y.Imag().HostWrite();
    for (int i = 0; i < height; i++)
    {
      const double d = 2.0 + static_cast<double>(i) / height - s;
      YR[i] =
          d * XR[i] - ((i > 0) ? XR[i - 1] : 0.0) - ((i < height - 1) ? XR[i + 1] : 0.0);
      YI[i] =
          d * XI[i] - ((i > 0) ? XI[i - 1] : 0.0) - ((i < height - 1) ? XI[i + 1] : 0.0);
    }
  }

  void AddMult(const ComplexVector &x, ComplexVector &y,
               const std::complex<double> a = 1.0) const override
  {
    ComplexVector t;
    Mult(x, t);
    y.Add(a, t);
  }
};

class IdentityEigenTestOperator : public ComplexOperator
{
public:
  IdentityEigenTestOperator(int n) : ComplexOperator(n) {}

  bool IsReal() const override { return true; }

  void Mult(const ComplexVector &x, ComplexVector &y) const override { y = x; }

  void AddMult(const ComplexVector &x, ComplexVector &y,
               const std::complex<double> a = 1.0) const override
  {
    y.Add(a, x);
  }
};

struct TridiagonalEigenTestResult
{
  std::vector<std::complex<double>> lambda;
  std::vector<ComplexVector> x;
  int its;      // Number of restarts
  int lin_its;  // Number of linear solver iterations
};

// Shift-and-invert Krylov-Schur solve for the nev eigenvalues of the tridiagonal operator
// closest to the shift s, with GMRES and Jacobi preconditioning for the linear solves, from
// the start vector v0, and with the linear solver tolerance first_tol for the first Krylov
// cycle if positive.
TridiagonalEigenTestResult SolveTridiagonalEigenTest(int n, double s, int nev, int ncv,
                                                     double tol, double lin_tol,
                                                     const ComplexVector &v0,
                                                     double first_tol = 0.0)
{
  MPI_Comm comm = Mpi::World();
  TridiagonalEigenTestOperator K(n), A(n, s);
  IdentityEigenTestOperator M(n);
  auto gmres = std::make_unique<GmresSolver<ComplexOperator>>(comm, 0);
  gmres->SetRelTol(lin_tol);
  gmres->SetMaxIter(10 * n);
  gmres->SetRestartDim(n);
  ComplexKspSolver ksp(std::move(gmres),
                       std::make_unique<JacobiSmoother<ComplexOperator>>(comm));
  ksp.SetOperators(A, A);
  slepc::SlepcEPSSolver eps(comm, 0);
  eps.SetType(slepc::SlepcEigenvalueSolver::Type::KRYLOVSCHUR);
  eps.SetProblemType(slepc::SlepcEigenvalueSolver::ProblemType::GEN_NON_HERMITIAN);
  eps.SetOperators(K, M, EigenvalueSolver::ScaleType::NONE);
  eps.SetNumModes(nev, ncv);
  eps.SetTol(tol);
  eps.SetWhichEigenpairs(EigenvalueSolver::WhichType::TARGET_REAL);
  eps.SetLinearSolver(ksp);
  eps.SetShiftInvert(s);
  eps.SetInitialSpace(v0);
  eps.SetFirstCycleLinearSolverTol(first_tol);
  const int num_conv = eps.Solve();
  REQUIRE(num_conv >= nev);
  TridiagonalEigenTestResult result;
  result.its = eps.GetNumIterations();
  result.lin_its = ksp.NumTotalMultIterations();
  for (int i = 0; i < nev; i++)
  {
    result.lambda.push_back(eps.GetEigenvalue(i));
    result.x.emplace_back(n);
    eps.GetEigenvector(i, result.x.back());
  }
  return result;
}

}  // namespace

TEST_CASE("Eigenmode warm start", "[eigenmode][Serial]")
{
  // From a start vector combining the wanted eigenvectors, the eigenvalues are only as
  // accurate as the linear solves of the first Krylov cycle, whose error the error estimate
  // of the eigenvalue solver does not see: a smaller tolerance for that cycle gives the
  // accuracy of exact solves. The shift is below the spectrum, so that GMRES converges for
  // any tolerance.
  constexpr int n = 120, nev = 2, ncv = 8;
  const double s = 0.0, tol = 1.0e-10, exact_tol = 1.0e-12, lin_tol = 1.0e-4;
  ComplexVector v0(n);
  linalg::SetRandom(Mpi::World(), v0);
  const auto exact = SolveTridiagonalEigenTest(n, s, nev, ncv, tol, exact_tol, v0);
  CHECK(exact.its > 1);

  ComplexVector v(n), r(n);
  linalg::SetRandom(Mpi::World(), r, 1);
  v = exact.x[0];
  v.Add(1.0, exact.x[1]);
  v.Add(1.0e-12 / std::sqrt(static_cast<double>(n)), r);
  const auto warm = SolveTridiagonalEigenTest(n, s, nev, ncv, tol, lin_tol, v);
  const auto warm_tight =
      SolveTridiagonalEigenTest(n, s, nev, ncv, tol, lin_tol, v, exact_tol);
  CHECK(warm_tight.its == 1);
  double err = 0.0, err_tight = 0.0;
  for (int i = 0; i < nev; i++)
  {
    const double a = std::abs(exact.lambda[i]);
    err = std::max(err, std::abs(warm.lambda[i] - exact.lambda[i]) / a);
    err_tight = std::max(err_tight, std::abs(warm_tight.lambda[i] - exact.lambda[i]) / a);
  }
  CHECK(err > 100.0 * tol);
  CHECK(err_tight <= tol);

  // The tolerance is restored after the first cycle: from a random vector, the solve takes
  // more linear solver iterations than without the smaller tolerance, and fewer than with
  // it for all cycles.
  const auto cold = SolveTridiagonalEigenTest(n, s, nev, ncv, tol, lin_tol, v0);
  const auto cold_tight =
      SolveTridiagonalEigenTest(n, s, nev, ncv, tol, lin_tol, v0, exact_tol);
  CHECK(cold_tight.lin_its > cold.lin_its);
  CHECK(cold_tight.lin_its < exact.lin_its);
}

#endif

}  // namespace palace
