// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "superlucomplex.hpp"

#if defined(MFEM_USE_SUPERLU)

#include <algorithm>
#include <cmath>
#include <complex>
#include <utility>
#include <vector>
#include "linalg/complexcsr.hpp"
#include "utils/communication.hpp"

// SuperLU_DIST complex double precision header.
#include "superlu_zdefs.h"

namespace palace
{

struct SuperLUComplexSolver::Data
{
  gridinfo_t grid;
  superlu_dist_options_t options;
  zScalePermstruct_t ScalePermstruct;
  zLUstruct_t LUstruct;
  zSOLVEstruct_t SOLVEstruct;
  SuperMatrix A;
  bool has_A = false, initialized = false;
  int_t n_glob = 0, n_loc = 0;

  // Sparsity pattern of the factored matrix, for the reuse of the row and column
  // permutations.
  std::vector<int_t> I, J;

  // Right-hand side and solution (overwritten in place).
  std::vector<doublecomplex> B;
  double berr = 0.0;

  void DestroyMatrix()
  {
    if (has_A)
    {
      Destroy_CompRowLoc_Matrix_dist(&A);
      has_A = false;
    }
  }

  // Create the distributed matrix from the local rows (SuperLU_DIST takes ownership of
  // copies of the arrays).
  void CreateMatrix(int_t first_row, const std::vector<std::complex<double>> &V)
  {
    DestroyMatrix();
    const int_t nnz_loc = static_cast<int_t>(J.size());
    auto *nzval = doublecomplexMalloc_dist(nnz_loc);
    auto *colind = intMalloc_dist(nnz_loc);
    auto *rowptr = intMalloc_dist(n_loc + 1);
    MFEM_VERIFY(nzval && colind && rowptr, "SuperLU_DIST matrix allocation failed!");
    std::copy(J.begin(), J.end(), colind);
    std::copy(I.begin(), I.end(), rowptr);
    std::copy(V.begin(), V.end(), reinterpret_cast<std::complex<double> *>(nzval));
    zCreate_CompRowLoc_Matrix_dist(&A, n_glob, n_glob, nnz_loc, n_loc, first_row, nzval,
                                   colind, rowptr, SLU_NR_loc, SLU_Z, SLU_GE);
    has_A = true;
  }

  void DestroyFactorization()
  {
    if (initialized)
    {
      zDestroy_LU(n_glob, &grid, &LUstruct);
      if (options.SolveInitialized == YES)
      {
        zSolveFinalize(&options, &SOLVEstruct);
      }
      zScalePermstructFree(&ScalePermstruct);
      zLUstructFree(&LUstruct);
      initialized = false;
    }
  }

  // Factor the matrix (nrhs = 0) or solve in place for the right-hand side B with the
  // factorization (nrhs = 1).
  void Call(int nrhs)
  {
    SuperLUStat_t stat;
    PStatInit(&stat);
    int info = -1;
    pzgssvx(&options, &A, &ScalePermstruct, B.data(), static_cast<int>(n_loc), nrhs, &grid,
            &LUstruct, &SOLVEstruct, &berr, &stat, &info);
    PStatFree(&stat);
    MFEM_VERIFY(info == 0, "SuperLU_DIST complex-valued "
                               << ((nrhs > 0) ? "solve" : "factorization")
                               << " failed (info = " << info << ")!");
    options.Fact = FACTORED;
  }
};

SuperLUComplexSolver::SuperLUComplexSolver(MPI_Comm comm, SymbolicFactorization reorder,
                                           bool reorder_reuse, int print)
  : Solver<ComplexOperator>(), comm(comm), reorder_reuse(reorder_reuse),
    data(std::make_unique<Data>())
{
  // 2D process grid as close to square as possible, as for SuperLUSolver.
  const int np = Mpi::Size(comm);
  int nprow = static_cast<int>(std::sqrt(static_cast<double>(np)));
  while (nprow > 1 && np % nprow != 0)
  {
    nprow--;
  }
  superlu_gridinit(comm, nprow, np / nprow, &data->grid);
  if (print > 1)
  {
    Mpi::Print(comm, " SuperLUComplexSolver: Using 2D processor grid {:d} x {:d}\n", nprow,
               np / nprow);
  }

  // Options as for SuperLUSolver.
  auto &options = data->options;
  set_default_options_dist(&options);
  options.PrintStat = (print > 1) ? YES : NO;
  options.Equil = NO;
  options.ReplaceTinyPivot = NO;
  options.IterRefine = NOREFINE;
  options.SymPattern = YES;
  switch (reorder)
  {
    case SymbolicFactorization::METIS:
      options.ColPerm = METIS_AT_PLUS_A;
      break;
    case SymbolicFactorization::PARMETIS:
      options.ColPerm = PARMETIS;
      options.ParSymbFact = YES;
      break;
    case SymbolicFactorization::AMD:
    case SymbolicFactorization::RCM:
      options.ColPerm = MMD_AT_PLUS_A;
      break;
    case SymbolicFactorization::SCOTCH:
    case SymbolicFactorization::PTSCOTCH:
    case SymbolicFactorization::PORD:
    case SymbolicFactorization::DEFAULT:
      // Should have good default.
      break;
  }

    // Relaxed supernodes of up to 200 columns, as for SuperLUSolver.
#if SUPERLU_DIST_MAJOR_VERSION >= 8
  options.superlu_relax = 200;
#endif
}

SuperLUComplexSolver::~SuperLUComplexSolver()
{
  data->DestroyFactorization();
  data->DestroyMatrix();
  superlu_gridexit(&data->grid);
}

void SuperLUComplexSolver::SetOperator(const ComplexOperator &op)
{
  // Local rows of the complex-valued matrix.
  ComplexLocalCSR csr(op);
  const int_t n_glob = static_cast<int_t>(csr.n_glob);
  const int_t n_loc = static_cast<int_t>(csr.n_loc);
  std::vector<int_t> I(csr.I.begin(), csr.I.end()), J(csr.J.begin(), csr.J.end());

  // Reuse the row and column permutations for a matrix with the same sparsity pattern.
  auto &d = *data;
  int same = d.initialized && reorder_reuse && n_glob == d.n_glob && I == d.I && J == d.J;
  Mpi::GlobalMin(1, &same, comm);
  if (same)
  {
    d.options.Fact = SamePattern_SameRowPerm;
    zZeroLblocks(d.grid.iam, n_glob, &d.grid, &d.LUstruct);
    zZeroUblocks(d.grid.iam, n_glob, &d.grid, &d.LUstruct);
  }
  else
  {
    d.DestroyFactorization();
    zScalePermstructInit(n_glob, n_glob, &d.ScalePermstruct);
    zLUstructInit(n_glob, &d.LUstruct);
    d.options.Fact = DOFACT;
    d.options.SolveInitialized = NO;
    d.options.RefineInitialized = NO;
    d.initialized = true;
  }
  d.n_glob = n_glob;
  d.n_loc = n_loc;
  d.I = std::move(I);
  d.J = std::move(J);
  d.CreateMatrix(static_cast<int_t>(csr.first_row), csr.V);
  d.B.resize(n_loc);
  d.Call(0);
  this->height = op.Height();
  this->width = op.Width();
}

void SuperLUComplexSolver::Mult(const ComplexVector &x, ComplexVector &y) const
{
  auto &d = *data;
  const int n = static_cast<int>(d.n_loc);
  const auto *xr = x.Real().HostRead(), *xi = x.Imag().HostRead();
  for (int i = 0; i < n; i++)
  {
    d.B[i] = {xr[i], xi[i]};
  }
  d.Call(1);
  auto *yr = y.Real().HostWrite(), *yi = y.Imag().HostWrite();
  for (int i = 0; i < n; i++)
  {
    yr[i] = d.B[i].r;
    yi[i] = d.B[i].i;
  }
}

}  // namespace palace

#endif
