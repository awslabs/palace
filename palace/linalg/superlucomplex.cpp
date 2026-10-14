// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "superlucomplex.hpp"

#if defined(MFEM_USE_SUPERLU)

#include <algorithm>
#include <cmath>
#include <complex>
#include <utility>
#include <vector>
#include "linalg/rap.hpp"
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

const mfem::HypreParMatrix *GetAssembled(const Operator *op)
{
  if (!op)
  {
    return nullptr;
  }
  if (const auto *hA = dynamic_cast<const mfem::HypreParMatrix *>(op))
  {
    return hA;
  }
  const auto *PtAP = dynamic_cast<const ParOperator *>(op);
  MFEM_VERIFY(PtAP, "SuperLUComplexSolver requires HypreParMatrix or ParOperator parts!");
  return &PtAP->ParallelAssemble();
}

}  // namespace

struct SuperLUComplexSolver::Data
{
  gridinfo_t grid;
  superlu_dist_options_t options;
  zScalePermstruct_t ScalePermstruct;
  zLUstruct_t LUstruct;
  zSOLVEstruct_t SOLVEstruct;
  SuperMatrix A;
  bool has_A = false, initialized = false;
  int_t n_glob = 0, n_loc = 0, first_row = 0;
  std::vector<int_t> I, J;
  std::vector<std::complex<double>> V;
  mutable std::vector<doublecomplex> B;
  mutable std::vector<double> berr;

  void DestroyMatrix()
  {
    if (has_A)
    {
      Destroy_CompRowLoc_Matrix_dist(&A);
      has_A = false;
    }
  }

  // Create the distributed matrix from the local rows (I, J, V). The solve with a
  // factorization overwrites the matrix (scaling and permutation), so that each
  // factorization requires a new matrix.
  void CreateMatrix()
  {
    DestroyMatrix();
    const int_t nnz_loc = static_cast<int_t>(J.size());
    auto *nzval = doublecomplexMalloc_dist(nnz_loc);
    auto *colind = intMalloc_dist(nnz_loc);
    auto *rowptr = intMalloc_dist(n_loc + 1);
    std::copy(J.begin(), J.end(), colind);
    std::copy(I.begin(), I.end(), rowptr);
    std::copy(V.begin(), V.end(), reinterpret_cast<std::complex<double> *>(nzval));
    zCreate_CompRowLoc_Matrix_dist(&A, n_glob, n_glob, nnz_loc, n_loc, first_row, nzval,
                                   colind, rowptr, SLU_NR_loc, SLU_Z, SLU_GE);
    has_A = true;
  }

  // Refactor with the same sparsity pattern and row permutation.
  void ResetFactors()
  {
    options.Fact = SamePattern_SameRowPerm;
    zZeroLblocks(grid.iam, n_glob, &grid, &LUstruct);
    zZeroUblocks(grid.iam, n_glob, &grid, &LUstruct);
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

  // Solve in place for the right-hand side B (with the factorization for Fact != FACTORED).
  void Solve() const
  {
    auto &self = const_cast<Data &>(*this);
    SuperLUStat_t stat;
    PStatInit(&stat);
    int info = -1;
    pzgssvx(&self.options, &self.A, &self.ScalePermstruct, self.B.data(),
            static_cast<int>(n_loc), 1, &self.grid, &self.LUstruct, &self.SOLVEstruct,
            self.berr.data(), &stat, &info);
    PStatFree(&stat);
    MFEM_VERIFY(info == 0, "SuperLU_DIST complex solve failed (info = " << info << ")!");
    self.options.Fact = FACTORED;
  }
};

SuperLUComplexSolver::SuperLUComplexSolver(MPI_Comm comm, SymbolicFactorization reorder,
                                           bool reorder_reuse, int print,
                                           bool robust_pivoting)
  : Solver<ComplexOperator>(), comm(comm), reorder(reorder), reorder_reuse(reorder_reuse),
    robust_pivoting(robust_pivoting), refine(false), print(print),
    data(std::make_unique<Data>())
{
  // 2D process grid as close to square as possible.
  const int np = Mpi::Size(comm);
  int nprow = static_cast<int>(std::sqrt(static_cast<double>(np)));
  while (nprow > 1 && np % nprow != 0)
  {
    nprow--;
  }
  superlu_gridinit(comm, nprow, np / nprow, &data->grid);

  // Options as for SuperLUSolver.
  auto &options = data->options;
  set_default_options_dist(&options);
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
      break;
  }
}

SuperLUComplexSolver::~SuperLUComplexSolver()
{
  data->DestroyFactorization();
  data->DestroyMatrix();
  superlu_gridexit(&data->grid);
}

void SuperLUComplexSolver::SetOperator(const ComplexOperator &op)
{
  const auto *hAr = GetAssembled(op.Real());
  const auto *hAi = GetAssembled(op.Imag());
  MFEM_VERIFY(hAr || hAi, "Empty ComplexOperator for SuperLUComplexSolver!");
  const auto &hA = hAr ? *hAr : *hAi;
  const int_t n_glob = static_cast<int_t>(hA.GetGlobalNumRows());
  const int_t n_loc = static_cast<int_t>(hA.Height());
  const int_t first_row = static_cast<int_t>(hA.GetRowStarts()[0]);

  // Local rows of the real and imaginary parts with global column indices, merged into a
  // complex-valued matrix (SuperLU_DIST takes ownership of these arrays).
  auto Merge = [](const mfem::HypreParMatrix *A) -> hypre_CSRMatrix *
  {
    if (!A)
    {
      return nullptr;
    }
    hypre_CSRMatrix *csr = hypre_MergeDiagAndOffd(
        (hypre_ParCSRMatrix *)const_cast<mfem::HypreParMatrix &>(*A));
    hypre_CSRMatrixMigrate(csr, HYPRE_MEMORY_HOST);
    return csr;
  };
  hypre_CSRMatrix *csr[2] = {Merge(hAr), Merge(hAi)};
  std::vector<int_t> I(n_loc + 1, 0), J;
  std::vector<std::complex<double>> V;
  std::vector<std::pair<HYPRE_BigInt, std::complex<double>>> row;
  for (int_t i = 0; i < n_loc; i++)
  {
    row.clear();
    for (int part = 0; part < 2; part++)
    {
      if (!csr[part])
      {
        continue;
      }
      const std::complex<double> s = (part == 0) ? 1.0 : std::complex<double>(0.0, 1.0);
      const HYPRE_Int *Ip = hypre_CSRMatrixI(csr[part]);
      const HYPRE_BigInt *Jp = hypre_CSRMatrixBigJ(csr[part]);
      const double *Dp = hypre_CSRMatrixData(csr[part]);
      for (HYPRE_Int k = Ip[i]; k < Ip[i + 1]; k++)
      {
        row.emplace_back(Jp[k], s * Dp[k]);
      }
    }
    std::sort(row.begin(), row.end(),
              [](const auto &a, const auto &b) { return a.first < b.first; });
    for (const auto &[j, v] : row)
    {
      if (static_cast<int_t>(J.size()) > I[i] && J.back() == static_cast<int_t>(j))
      {
        V.back() += v;
      }
      else
      {
        J.push_back(static_cast<int_t>(j));
        V.push_back(v);
      }
    }
    I[i + 1] = static_cast<int_t>(J.size());
  }
  for (auto *c : csr)
  {
    if (c)
    {
      hypre_CSRMatrixDestroy(c);
    }
  }
  for (const auto *part : {op.Real(), op.Imag()})
  {
    if (const auto *PtAP = dynamic_cast<const ParOperator *>(part))
    {
      PtAP->StealParallelAssemble();
    }
  }

  // Reuse the row and column permutations for a matrix with the same sparsity pattern.
  int same = data->initialized && reorder_reuse && n_glob == data->n_glob && I == data->I &&
             J == data->J;
  Mpi::GlobalMin(1, &same, comm);
  auto &options = data->options;
  if (same)
  {
    data->ResetFactors();
  }
  else
  {
    data->DestroyFactorization();
    zScalePermstructInit(n_glob, n_glob, &data->ScalePermstruct);
    zLUstructInit(n_glob, &data->LUstruct);
    options.Fact = DOFACT;
    options.SolveInitialized = NO;
    options.RefineInitialized = NO;
    data->initialized = true;
  }
  data->n_glob = n_glob;
  data->n_loc = n_loc;
  data->first_row = first_row;
  data->I = std::move(I);
  data->J = std::move(J);
  data->V = std::move(V);
  data->CreateMatrix();
  data->B.resize(n_loc);
  data->berr.resize(1);
  this->height = op.Height();
  this->width = op.Width();

  // Factor and check the accuracy of the factorization (robust pivoting).
  double res = Factor(op);
  if (robust_pivoting && res > robust_pivoting_tol && !refine)
  {
    if (print > 0)
    {
      Mpi::Print(comm,
                 " SuperLUComplexSolver: Inaccurate factorization (relative residual "
                 "{:.3e}), using tiny pivot replacement and iterative refinement\n",
                 res);
    }
    refine = true;
    options.ReplaceTinyPivot = YES;
    options.IterRefine = SLU_DOUBLE;
    data->ResetFactors();
    data->CreateMatrix();
    res = Factor(op);
  }
  if (robust_pivoting && res > robust_pivoting_tol)
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

double SuperLUComplexSolver::Factor(const ComplexOperator &op)
{
  // The first solve with Fact != FACTORED performs the factorization.
  const int n = static_cast<int>(data->n_loc);
  ComplexVector x(n), b(n), r(n);
  x.Real().Randomize(1);
  x.Imag().Randomize(2);
  op.Mult(x, b);
  const auto *br = b.Real().HostRead(), *bi = b.Imag().HostRead();
  for (int i = 0; i < n; i++)
  {
    data->B[i] = {br[i], bi[i]};
  }
  data->Solve();
  auto *xr = x.Real().HostWrite(), *xi = x.Imag().HostWrite();
  for (int i = 0; i < n; i++)
  {
    xr[i] = data->B[i].r;
    xi[i] = data->B[i].i;
  }
  op.Mult(x, r);
  r -= b;
  return linalg::Norml2(comm, r) / linalg::Norml2(comm, b);
}

void SuperLUComplexSolver::Mult(const ComplexVector &x, ComplexVector &y) const
{
  const int n = static_cast<int>(data->n_loc);
  const auto *xr = x.Real().HostRead(), *xi = x.Imag().HostRead();
  for (int i = 0; i < n; i++)
  {
    data->B[i] = {xr[i], xi[i]};
  }
  data->Solve();
  auto *yr = y.Real().HostWrite(), *yi = y.Imag().HostWrite();
  for (int i = 0; i < n; i++)
  {
    yr[i] = data->B[i].r;
    yi[i] = data->B[i].i;
  }
}

}  // namespace palace

#endif
