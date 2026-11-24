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
#include "linalg/superlu.hpp"
#include "linalg/vector.hpp"
#include "utils/communication.hpp"

// SuperLU_DIST complex double precision header.
#include "superlu_zdefs.h"

namespace palace
{

namespace
{

// Relative residual above which a factorization with robust pivoting is considered
// inaccurate (as for SuperLUSolver).
constexpr double robust_pivoting_tol = 1.0e-10;

}  // namespace

struct SuperLUComplexSolver::Data
{
  // Process grid: 2D, or 3D with a depth above one.
  gridinfo_t grid;
  gridinfo3d_t grid3d;
  int npdep = 1;
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

  // The 2D process grid, or the one of the process layer with a 3D process grid.
  gridinfo_t &Grid2D() { return (npdep > 1) ? grid3d.grid2d : grid; }

  // With a 3D process grid, all the process layers hold L and U factors (those of their
  // part of the elimination tree) and the solve structures, and the matrix gathered on the
  // process layers is kept between the calls of the driver (freed in the order of
  // mfem::SuperLUSolver).
  void DestroyFactorization()
  {
    if (initialized)
    {
      if (npdep > 1)
      {
        if (options.SolveInitialized == YES)
        {
          zSolveFinalize(&options, &SOLVEstruct);
        }
        zDestroy_LU(n_glob, &grid3d.grid2d, &LUstruct);
        zDestroy_A3d_gathered_on_2d(&SOLVEstruct, &grid3d);
        SOLVEstruct.A3d = nullptr;
      }
      else
      {
        zDestroy_LU(n_glob, &grid, &LUstruct);
        if (options.SolveInitialized == YES)
        {
          zSolveFinalize(&options, &SOLVEstruct);
        }
      }
      zScalePermstructFree(&ScalePermstruct);
      zLUstructFree(&LUstruct);
      initialized = false;
    }
  }

  // Factor the matrix (nrhs = 0, or nrhs = 1 to also solve) or solve with the
  // factorization (nrhs = 1), in place for the right-hand side B.
  void Call(int nrhs)
  {
    const bool factor = (options.Fact != FACTORED);
    SuperLUStat_t stat;
    PStatInit(&stat);
    int info = -1;
    auto *b = reinterpret_cast<doublecomplex *>(B.data());
    if (npdep > 1)
    {
      pzgssvx3d(&options, &A, &ScalePermstruct, b, static_cast<int>(n_loc), nrhs, &grid3d,
                &LUstruct, &SOLVEstruct, &berr, &stat, &info);
    }
    else
    {
      pzgssvx(&options, &A, &ScalePermstruct, b, static_cast<int>(n_loc), nrhs, &grid,
              &LUstruct, &SOLVEstruct, &berr, &stat, &info);
    }
    PStatFree(&stat);
    MFEM_VERIFY(info == 0, "SuperLU_DIST complex-valued "
                               << (factor ? "factorization" : "solve")
                               << " failed (info = " << info << ")!");
    options.Fact = FACTORED;
  }
};

SuperLUComplexSolver::SuperLUComplexSolver(MPI_Comm comm, SymbolicFactorization reorder,
                                           bool use_3d, bool drop_small_entries,
                                           bool reorder_reuse, int print,
                                           bool robust_pivoting)
  : Solver<ComplexOperator>(), comm(comm), drop_small_entries(drop_small_entries),
    reorder_reuse(reorder_reuse), robust_pivoting(robust_pivoting), refine(false),
    print(print), data(std::make_unique<Data>())
{
  // Process grid as for SuperLUSolver: the depth of a 3D process grid, and process layers
  // as close to square as possible.
  const int np = Mpi::Size(comm);
  const int npdep = data->npdep = GetSuperLUGridDepth(np, use_3d);
  MFEM_VERIFY(!robust_pivoting || npdep == 1,
              "SuperLUComplexSolver supports robust pivoting on 2D process grids only!");
#if SUPERLU_DIST_MAJOR_VERSION < 9
  MFEM_VERIFY(npdep == 1, "SuperLUComplexSolver requires SuperLU_DIST 9 or later for a "
                          "3D process grid!");
#endif
  int nprow = static_cast<int>(std::sqrt(static_cast<double>(np / npdep)));
  while (nprow > 1 && (np / npdep) % nprow != 0)
  {
    nprow--;
  }
  const int npcol = np / (npdep * nprow);
  if (npdep > 1)
  {
    superlu_gridinit3d(comm, nprow, npcol, npdep, &data->grid3d);
  }
  else
  {
    superlu_gridinit(comm, nprow, npcol, &data->grid);
  }
  if (print > 1)
  {
    if (npdep > 1)
    {
      Mpi::Print(comm,
                 " SuperLUComplexSolver: Using 3D processor grid {:d} x {:d} x {:d}\n",
                 nprow, npcol, npdep);
    }
    else
    {
      Mpi::Print(comm, " SuperLUComplexSolver: Using 2D processor grid {:d} x {:d}\n",
                 nprow, npcol);
    }
  }

  // Options as for SuperLUSolver.
  auto &options = data->options;
  set_default_options_dist(&options);
  if (npdep > 1)
  {
    options.Algo3d = YES;
  }
  options.PrintStat = (print > 1) ? YES : NO;
  options.Equil = NO;
  options.ReplaceTinyPivot = NO;
  options.IterRefine = NOREFINE;
  options.SymPattern = robust_pivoting ? NO : YES;
  switch (reorder)
  {
    case SymbolicFactorization::METIS:
      options.ColPerm = METIS_AT_PLUS_A;
      break;
    case SymbolicFactorization::PARMETIS:
      // SuperLU_DIST does not support ParMETIS (parallel symbolic factorization) with a 3D
      // process grid and returns an error for it.
      if (npdep > 1)
      {
        Mpi::Warning(comm, "SuperLU_DIST does not support ParMETIS ordering with a 3D "
                           "process grid, using METIS instead!\n");
        options.ColPerm = METIS_AT_PLUS_A;
      }
      else
      {
        options.ColPerm = PARMETIS;
        options.ParSymbFact = YES;
      }
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
  if (data->npdep > 1)
  {
    superlu_gridexit3d(&data->grid3d);
  }
  else
  {
    superlu_gridexit(&data->grid);
  }
}

void SuperLUComplexSolver::SetOperator(const ComplexOperator &op)
{
  auto &d = *data;
  {
    // Local rows of the complex-valued matrix (freed before the factorization).
    ComplexLocalCSR csr(op, drop_small_entries);

    // SuperLU_DIST gets a copy of the matrix: release the parallel assembled parts (as
    // MfemWrapperSolver does for the sparse direct solvers).
    ReleaseParallelAssemble(op);

    // Reuse the row and column permutations for a matrix with the same sparsity pattern on
    // all processes (when small entries are dropped, the sparsity pattern can change with
    // the values). With a 3D process grid, the factors are zeroed on process layer 0 (the
    // driver fills and zeroes those of the other layers), as in mfem::SuperLUSolver.
    const int_t n_glob = static_cast<int_t>(csr.n_glob);
    if (csr.SamePattern(comm, d.pattern) && d.initialized && reorder_reuse)
    {
      d.options.Fact = SamePattern_SameRowPerm;
      if (d.npdep == 1 || d.grid3d.zscp.Iam == 0)
      {
        auto &grid = d.Grid2D();
        zZeroLblocks(grid.iam, n_glob, &grid, &d.LUstruct);
        zZeroUblocks(grid.iam, n_glob, &grid, &d.LUstruct);
      }
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
    d.B.assign(csr.n_loc, 0.0);
  }

  this->height = op.Height();
  this->width = op.Width();
  if (!robust_pivoting)
  {
    // With a 3D process grid, the factorization is followed by a solve (of a zero
    // right-hand side) in the same call, as in mfem::SuperLUSolver, whose first solve
    // factors the matrix. This works around a defect of SuperLU_DIST (9.3.1): the driver
    // for 3D process grids allocates part of its solve structures in the first solve, and
    // zDestroy_A3d_gathered_on_2d frees them whether or not they were allocated, which is
    // an invalid free when a factorization is destroyed or replaced before any solve.
    d.Call((d.npdep > 1) ? 1 : 0);
    return;
  }

  // Check the accuracy of the factorization, and factor again with tiny pivot replacement
  // (which bounds the element growth of the elimination) and iterative refinement (which
  // recovers the accuracy lost to the perturbed pivots) if it is inaccurate (see
  // SuperLUSolver). The first solve of the check factors the matrix: SuperLU_DIST sets up
  // the iterative refinement for a new matrix only in a call which factors it, and a later
  // solve would refine with the setup of the previous matrix.
  const double res = FactorizationResidual(op);
  if (res > robust_pivoting_tol && !refine)
  {
    if (print > 0)
    {
      Mpi::Print(comm,
                 " SuperLUComplexSolver: Inaccurate factorization (relative residual "
                 "{:.3e}), using tiny pivot replacement and iterative refinement\n",
                 res);
    }
    refine = true;
    d.options.ReplaceTinyPivot = YES;
    d.options.IterRefine = SLU_DOUBLE;
    SetOperator(op);
    return;
  }
  if (res > robust_pivoting_tol)
  {
    Mpi::Warning(comm,
                 "SuperLU_DIST factorization is inaccurate (relative residual {:.3e})!\n"
                 "Consider using STRUMPACK or MUMPS instead.\n",
                 res);
  }
  else if (print > 1)
  {
    Mpi::Print(comm, " SuperLUComplexSolver: Factorization relative residual {:.3e}\n",
               res);
  }
}

double SuperLUComplexSolver::FactorizationResidual(const ComplexOperator &op) const
{
  ComplexVector x(op.Width()), b(op.Height()), r(op.Height());
  x.UseDevice(true);
  b.UseDevice(true);
  r.UseDevice(true);
  x.Real().Randomize(1);
  x.Imag().Randomize(2);
  op.Mult(x, b);
  Mult(b, x);
  op.Mult(x, r);
  r -= b;
  return linalg::Norml2(comm, r) / linalg::Norml2(comm, b);
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
