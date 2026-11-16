// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "superlucomplex.hpp"

#if defined(MFEM_USE_SUPERLU)

#include <algorithm>
#include <cmath>
#include <complex>
#include <cstddef>
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

  // Hash of the sparsity pattern of the factored matrix, for the reuse of the row and
  // column permutations.
  std::size_t pattern = 0;

  // Right-hand side and solution (overwritten in place).
  std::vector<std::complex<double>> B;
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
  void CreateMatrix(const ComplexLocalCSR &csr)
  {
    DestroyMatrix();
    const int_t nnz_loc = static_cast<int_t>(csr.J.size());
    auto *nzval = doublecomplexMalloc_dist(nnz_loc);
    auto *colind = intMalloc_dist(nnz_loc);
    auto *rowptr = intMalloc_dist(n_loc + 1);
    MFEM_VERIFY(nzval && colind && rowptr, "SuperLU_DIST matrix allocation failed!");
    std::copy(csr.J.begin(), csr.J.end(), colind);
    std::copy(csr.I.begin(), csr.I.end(), rowptr);
    std::copy(csr.V.begin(), csr.V.end(), reinterpret_cast<std::complex<double> *>(nzval));
    zCreate_CompRowLoc_Matrix_dist(&A, n_glob, n_glob, nnz_loc, n_loc,
                                   static_cast<int_t>(csr.first_row), nzval, colind, rowptr,
                                   SLU_NR_loc, SLU_Z, SLU_GE);
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
    pzgssvx(&options, &A, &ScalePermstruct, reinterpret_cast<doublecomplex *>(B.data()),
            static_cast<int>(n_loc), nrhs, &grid, &LUstruct, &SOLVEstruct, &berr, &stat,
            &info);
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
  auto &d = *data;
  {
    // Local rows of the complex-valued matrix (freed before the factorization).
    ComplexLocalCSR csr(op);

    // SuperLU_DIST gets a copy of the matrix: release the parallel assembled parts (as
    // MfemWrapperSolver does for the sparse direct solvers).
    ReleaseParallelAssemble(op);

    // Reuse the row and column permutations for a matrix with the same sparsity pattern on
    // all processes.
    const int_t n_glob = static_cast<int_t>(csr.n_glob);
    if (csr.SamePattern(comm, d.pattern) && d.initialized && reorder_reuse)
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
    d.n_loc = static_cast<int_t>(csr.n_loc);
    d.CreateMatrix(csr);
    d.B.resize(csr.n_loc);
  }
  d.Call(0);
  this->height = op.Height();
  this->width = op.Width();
}

void SuperLUComplexSolver::Mult(const ComplexVector &x, ComplexVector &y) const
{
  auto &d = *data;
  x.Get(d.B.data(), static_cast<int>(d.B.size()), false);
  d.Call(1);
  y.Set(d.B.data(), static_cast<int>(d.B.size()), false);
}

}  // namespace palace

#endif
