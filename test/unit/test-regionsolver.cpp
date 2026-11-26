// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include <cmath>
#include <complex>
#include <memory>
#include <string>
#include <type_traits>
#include <vector>
#include <mfem.hpp>
#include <catch2/catch_template_test_macros.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>
#include <catch2/generators/catch_generators_range.hpp>
#include "linalg/blockcyclic.hpp"
#include "linalg/lowertriangle.hpp"
#include "linalg/vector.hpp"
#include "models/regionsolver.hpp"
#include "utils/communication.hpp"

namespace palace
{

namespace
{

using cplx = std::complex<double>;

// A pseudo-random value in [-1, 1] x [-1, 1] i from (i, j).
cplx Value(int i, int j)
{
  const double a = std::sin(12.9898 * i + 78.233 * j + 1.0) * 43758.5453,
               b = std::sin(39.346 * i + 11.135 * j + 2.0) * 23421.631;
  return {2.0 * (a - std::floor(a)) - 1.0, 2.0 * (b - std::floor(b)) - 1.0};
}

template <typename T>
T Entry(cplx v)
{
  if constexpr (std::is_same_v<T, double>)
  {
    return v.real();
  }
  else
  {
    return v;
  }
}

// Entry i of a real vector, or of a complex one.
cplx Get(const mfem::Vector &x, int i)
{
  return x(i);
}
cplx Get(const ComplexVector &x, int i)
{
  return {x.Real()(i), x.Imag()(i)};
}
void Set(mfem::Vector &x, int i, cplx v)
{
  x(i) = v.real();
}
void Set(ComplexVector &x, int i, cplx v)
{
  x.Real()(i) = v.real();
  x.Imag()(i) = v.imag();
}

}  // namespace

TEMPLATE_TEST_CASE("RegionSolver solves a system with a dense interface term",
                   "[substructure][Serial][Parallel]", double, cplx)
{
  // The system K + M (real, positive definite) or K + i M (complex symmetric) of H1 on a
  // square, with the DOFs g % 5 == 2 as Γ (in reverse order) and those g % 7 == 3 off it
  // pinned (outside the system), and a dense term X on Γ, by each sparse direct solver on
  // grids of one process and of all, in small blocks: the solutions of two right-hand
  // sides, in place, have residuals at rounding and are 0 on the pinned DOFs. Without Γ,
  // it solves the system.
  using T = TestType;
  using VecType = typename RegionSolver<T>::VecType;
  constexpr bool real = std::is_same_v<T, double>;
  MPI_Comm comm = Mpi::World();
  mfem::Mesh serial = mfem::Mesh::MakeCartesian2D(6, 6, mfem::Element::QUADRILATERAL);
  mfem::ParMesh pmesh(comm, serial);
  mfem::H1_FECollection fec(2, 2);
  mfem::ParFiniteElementSpace fes(&pmesh, &fec);
  mfem::ConstantCoefficient one(1.0);
  auto assemble = [&](mfem::BilinearFormIntegrator *integ)
  {
    mfem::ParBilinearForm a(&fes);
    a.AddDomainIntegrator(integ);
    a.Assemble();
    a.Finalize();
    return std::unique_ptr<mfem::HypreParMatrix>(a.ParallelAssemble());
  };
  const auto K = assemble(new mfem::DiffusionIntegrator(one));
  const auto M = assemble(new mfem::MassIntegrator(one));
  const int nt = K->Height(), n = static_cast<int>(K->GetGlobalNumRows());
  const HYPRE_BigInt row0 = K->GetRowStarts()[0];

  const bool with_gamma = GENERATE(true, false);
  auto on_gamma = [&](HYPRE_BigInt g) { return with_gamma && g % 5 == 2; };
  std::vector<HYPRE_BigInt> gamma;
  for (HYPRE_BigInt g = n - 1; g >= 0; g--)
  {
    if (on_gamma(g))
    {
      gamma.push_back(g);
    }
  }
  std::vector<char> keep(nt);
  for (int i = 0; i < nt; i++)
  {
    keep[i] = on_gamma(row0 + i) || (row0 + i) % 7 != 3;
  }
  const int nG = static_cast<int>(gamma.size());
  // The term X, symmetric, each entry computed once.
  std::vector<cplx> Xs(static_cast<std::size_t>(nG) * nG);
  for (int a = 0; a < nG; a++)
  {
    for (int b = 0; b <= a; b++)
    {
      Xs[static_cast<std::size_t>(a) * nG + b] = Xs[static_cast<std::size_t>(b) * nG + a] =
          0.1 * Value(a, b) + ((a == b) ? cplx(1.0 + 0.1 * nG, 0.0) : cplx(0.0));
    }
  }
  auto X = [&](int a, int b) { return Xs[static_cast<std::size_t>(a) * nG + b]; };

  std::vector<LinearSolver> types;
#if defined(MFEM_USE_MUMPS)
  types.push_back(LinearSolver::MUMPS);
#endif
#if defined(MFEM_USE_STRUMPACK)
  types.push_back(LinearSolver::STRUMPACK);
#endif
#if defined(MFEM_USE_SUPERLU)
  types.push_back(LinearSolver::SUPERLU);
#endif
  if (types.empty())
  {
    return;
  }
  const LinearSolver type = GENERATE_COPY(from_range(types));
  const int nproc = GENERATE_COPY(1, Mpi::Size(comm));
  CAPTURE(static_cast<int>(type), nproc, with_gamma);
  BlockCyclicGrid grid(comm, nG, nproc, 4);
  SideSystem side({K.get(), M.get()}, keep);
  auto val = side.Values<T>({T(1.0), Entry<T>(real ? cplx(1.0) : cplx(0.0, 1.0))});
  RegionSolver<T> solver(grid, type, SymbolicFactorization::DEFAULT, side.n_glob,
                         static_cast<int>(side.rows.size()),
                         {std::move(side.irn), std::move(side.jcn), std::move(val)},
                         side.Indices(gamma), side.rows,
                         [&](std::vector<T> &C)
                         {
                           for (int jl = 0; jl < grid.LocalCols(); jl++)
                           {
                             for (int il = 0; il < grid.LocalRows(); il++)
                             {
                               C[il + static_cast<std::size_t>(grid.LeadingDim()) * jl] +=
                                   Entry<T>(X(grid.GlobalRow(il), grid.GlobalCol(jl)));
                             }
                           }
                         });

  std::vector<VecType> b(2, VecType(nt)), y(2, VecType(nt));
  for (int k = 0; k < 2; k++)
  {
    for (int i = 0; i < nt; i++)
    {
      Set(b[k], i, Value(static_cast<int>(row0) + i, k));
    }
    y[k] = b[k];
  }
  solver.Mult({&y[0], &y[1]}, {&y[0], &y[1]});

  // The residual b - (K + c M) y - E X E^T y on the system, and y on the pinned DOFs.
  const cplx c = real ? cplx(1.0) : cplx(0.0, 1.0);
  for (int k = 0; k < 2; k++)
  {
    mfem::Vector yr(nt), yi(nt), t(nt);
    for (int i = 0; i < nt; i++)
    {
      yr(i) = Get(y[k], i).real();
      yi(i) = Get(y[k], i).imag();
    }
    std::vector<cplx> Ay(nt, 0.0);
    auto add = [&](const mfem::HypreParMatrix &A, const mfem::Vector &v, cplx a)
    {
      A.Mult(v, t);
      for (int i = 0; i < nt; i++)
      {
        Ay[i] += a * t(i);
      }
    };
    add(*K, yr, 1.0);
    add(*K, yi, cplx(0.0, 1.0));
    add(*M, yr, c);
    add(*M, yi, c * cplx(0.0, 1.0));
    std::vector<cplx> yg(nG, 0.0);
    for (int a = 0; a < nG; a++)
    {
      if (gamma[a] >= row0 && gamma[a] < row0 + nt)
      {
        yg[a] = Get(y[k], static_cast<int>(gamma[a] - row0));
      }
    }
    MPI_Allreduce(MPI_IN_PLACE, yg.data(), 2 * nG, MPI_DOUBLE, MPI_SUM, comm);
    for (int a = 0; a < nG; a++)
    {
      if (gamma[a] >= row0 && gamma[a] < row0 + nt)
      {
        for (int q = 0; q < nG; q++)
        {
          Ay[gamma[a] - row0] += Entry<T>(X(a, q)) * yg[q];
        }
      }
    }
    double res[3] = {0.0, 0.0, 0.0};  // the residual, the right-hand side, y pinned
    for (int i = 0; i < nt; i++)
    {
      if (keep[i])
      {
        res[0] += std::norm(Get(b[k], i) - Ay[i]);
        res[1] += std::norm(Get(b[k], i));
      }
      else
      {
        res[2] += std::norm(Get(y[k], i));
      }
    }
    Mpi::GlobalSum(3, res, comm);
    CHECK(std::sqrt(res[0] / res[1]) < 1.0e-10);
    CHECK(res[2] == 0.0);
  }
}

}  // namespace palace
