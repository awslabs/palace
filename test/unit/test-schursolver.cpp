// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include <algorithm>
#include <cmath>
#include <memory>
#include <string>
#include <vector>
#include <mfem.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>
#include <catch2/generators/catch_generators_range.hpp>
#include "linalg/mumpsschur.hpp"
#include "linalg/schursolver.hpp"
#include "linalg/strumpackschur.hpp"
#include "utils/communication.hpp"

namespace palace
{

namespace
{

// The dense matrix of A, on every rank.
mfem::DenseMatrix Dense(const mfem::HypreParMatrix &A)
{
  const int n = static_cast<int>(A.GetGlobalNumRows()), n_loc = A.Height();
  const HYPRE_BigInt row0 = A.GetRowStarts()[0];
  mfem::DenseMatrix D(n);
  D = 0.0;
  mfem::Vector e(n_loc), col(n_loc);
  for (int j = 0; j < n; j++)
  {
    e = 0.0;
    if (j >= row0 && j < row0 + n_loc)
    {
      e(static_cast<int>(j - row0)) = 1.0;
    }
    A.Mult(e, col);
    for (int i = 0; i < n_loc; i++)
    {
      D(static_cast<int>(row0) + i, j) = col(i);
    }
  }
  MPI_Allreduce(MPI_IN_PLACE, D.Data(), n * n, MPI_DOUBLE, MPI_SUM, A.GetComm());
  return D;
}

// The global vector of the local parts x, on every rank.
mfem::Vector Gather(const mfem::Vector &x, int n, HYPRE_BigInt row0, MPI_Comm comm)
{
  mfem::Vector g(n);
  g = 0.0;
  for (int i = 0; i < x.Size(); i++)
  {
    g(static_cast<int>(row0) + i) = x(i);
  }
  MPI_Allreduce(MPI_IN_PLACE, g.GetData(), n, MPI_DOUBLE, MPI_SUM, comm);
  return g;
}

}  // namespace

TEST_CASE("Schur solvers give the Schur complement and the internal solves",
          "[schursolver][Serial][Parallel]")
{
  MPI_Comm comm = Mpi::World();
  mfem::Mesh serial = mfem::Mesh::MakeCartesian2D(6, 6, mfem::Element::QUADRILATERAL);
  mfem::ParMesh mesh(comm, serial);
  serial.Clear();
  mfem::H1_FECollection fec(2, 2);
  mfem::ParFiniteElementSpace fes(&mesh, &fec);
  mfem::ConstantCoefficient one(1.0);
  mfem::ParBilinearForm a(&fes);
  a.AddDomainIntegrator(new mfem::DiffusionIntegrator(one));
  a.AddDomainIntegrator(new mfem::MassIntegrator(one));
  a.Assemble();
  a.Finalize();
  std::unique_ptr<mfem::HypreParMatrix> A(a.ParallelAssemble());
  const int n = static_cast<int>(A->GetGlobalNumRows()), n_loc = A->Height();
  const HYPRE_BigInt row0 = A->GetRowStarts()[0];

  // Scattered Schur variables, in an order other than that of their rows, and rows of
  // neither kind (as Dirichlet rows), eliminated for the partial factorizations.
  std::vector<HYPRE_BigInt> schur_vars;
  std::vector<int> internal;
  std::vector<int> pos(n, -1);
  std::vector<char> other(n, 0);
  for (int g = n - 1; g >= 0; g--)
  {
    if (g % 5 == 2)
    {
      pos[g] = static_cast<int>(schur_vars.size());
      schur_vars.push_back(g);
    }
  }
  for (int g = 0; g < n; g++)
  {
    other[g] = (pos[g] < 0 && g % 7 == 3);
    if (pos[g] < 0 && !other[g])
    {
      internal.push_back(g);
    }
  }
  const int ns = static_cast<int>(schur_vars.size()),
            ni = static_cast<int>(internal.size());
  mfem::Array<int> other_rows;
  for (int i = 0; i < n_loc; i++)
  {
    if (other[row0 + i])
    {
      other_rows.Append(i);
    }
  }
  mfem::HypreParMatrix A_elim(*A);
  A_elim.EliminateBC(other_rows, mfem::Operator::DiagonalPolicy::DIAG_ONE);

  // The dense references: A_II^-1 and S = A_SS - A_SI A_II^-1 A_IS.
  const mfem::DenseMatrix D = Dense(*A);
  mfem::DenseMatrix A_II(ni), A_IS(ni, ns), S_ref(ns);
  for (int i = 0; i < ni; i++)
  {
    for (int j = 0; j < ni; j++)
    {
      A_II(i, j) = D(internal[i], internal[j]);
    }
    for (int j = 0; j < ns; j++)
    {
      A_IS(i, j) = D(internal[i], static_cast<int>(schur_vars[j]));
    }
  }
  for (int i = 0; i < ns; i++)
  {
    for (int j = 0; j < ns; j++)
    {
      S_ref(i, j) = D(static_cast<int>(schur_vars[i]), static_cast<int>(schur_vars[j]));
    }
  }
  mfem::DenseMatrixInverse A_II_inv(A_II);
  mfem::DenseMatrix Y(ni, ns), AtY(ns);
  A_II_inv.Mult(A_IS, Y);
  mfem::MultAtB(A_IS, Y, AtY);
  S_ref.Add(-1.0, AtY);

  // An exact internal solver for the back-solves.
  auto internal_solve =
      [&](const std::vector<const mfem::Vector *> &X, const std::vector<mfem::Vector *> &Yv)
  {
    for (std::size_t k = 0; k < X.size(); k++)
    {
      const mfem::Vector x = Gather(*X[k], n, row0, comm);
      mfem::Vector b(ni), y(ni);
      for (int i = 0; i < ni; i++)
      {
        b(i) = x(internal[i]);
      }
      A_II_inv.Mult(b, y);
      mfem::Vector full(n);
      full = 0.0;
      for (int i = 0; i < ni; i++)
      {
        full(internal[i]) = y(i);
      }
      Yv[k]->SetSize(n_loc);
      for (int i = 0; i < n_loc; i++)
      {
        (*Yv[k])(i) = full(static_cast<int>(row0) + i);
      }
    }
  };

  std::vector<std::string> backends = {"back-solves"};
#if defined(MFEM_USE_MUMPS)
  backends.push_back("MUMPS");
#endif
#if defined(PALACE_WITH_STRUMPACK_SCHUR)
  backends.push_back("STRUMPACK");
#endif
  const std::string backend = GENERATE_COPY(from_range(backends));
  CAPTURE(backend, ns);
  std::unique_ptr<SchurSolver<double>> solver;
  if (backend == "back-solves")
  {
    // Batches of 4: the last one is partial.
    REQUIRE(ns % 4 != 0);
    solver = std::make_unique<BackSolveSchurSolver>(*A, schur_vars, internal_solve, 4);
  }
#if defined(MFEM_USE_MUMPS)
  if (backend == "MUMPS")
  {
    solver = std::make_unique<MumpsSchurSolver>(A_elim, schur_vars);
  }
#endif
#if defined(PALACE_WITH_STRUMPACK_SCHUR)
  if (backend == "STRUMPACK")
  {
    solver = std::make_unique<StrumpackSchurSolver>(A_elim, schur_vars);
  }
#endif

  // The Schur complement, on rank 0.
  const std::vector<double> S = solver->ReleaseSchur();
  if (Mpi::Root(comm))
  {
    REQUIRE(S.size() == static_cast<std::size_t>(ns) * ns);
    double err = 0.0;
    for (int j = 0; j < ns; j++)
    {
      for (int i = 0; i < ns; i++)
      {
        err =
            std::max(err, std::abs(S[i + static_cast<std::size_t>(ns) * j] - S_ref(i, j)));
      }
    }
    CHECK(err <= 1.0e-12 * S_ref.MaxMaxNorm());
  }
  else
  {
    CHECK(S.empty());
  }

  // Two internal solves in place: A_II^-1 b on the internal variables, 0 on the Schur ones
  // (on the other rows, 0 for the back-solves and b for the partial factorizations).
  std::vector<mfem::Vector> x(2, mfem::Vector(n_loc));
  for (int i = 0; i < n_loc; i++)
  {
    x[0](i) = std::sin(1.0 + static_cast<double>(row0 + i));
    x[1](i) = std::cos(2.0 * static_cast<double>(row0 + i));
  }
  std::vector<mfem::Vector> ref(2);
  for (int k = 0; k < 2; k++)
  {
    const mfem::Vector g = Gather(x[k], n, row0, comm);
    mfem::Vector b(ni), y(ni);
    for (int i = 0; i < ni; i++)
    {
      b(i) = g(internal[i]);
    }
    A_II_inv.Mult(b, y);
    ref[k].SetSize(n);
    ref[k] = 0.0;
    for (int i = 0; i < ni; i++)
    {
      ref[k](internal[i]) = y(i);
    }
  }
  solver->SolveInternal({&x[0], &x[1]}, {&x[0], &x[1]});
  for (int k = 0; k < 2; k++)
  {
    double err = 0.0;
    for (int i = 0; i < n_loc; i++)
    {
      if (!other[row0 + i])
      {
        err = std::max(err, std::abs(x[k](i) - ref[k](static_cast<int>(row0) + i)));
      }
    }
    MPI_Allreduce(MPI_IN_PLACE, &err, 1, MPI_DOUBLE, MPI_MAX, comm);
    CHECK(err <= 1.0e-12 * ref[k].Normlinf());
  }
}

}  // namespace palace
