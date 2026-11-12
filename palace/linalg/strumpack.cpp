// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "strumpack.hpp"

#include <cstdio>
#include <cstdlib>
#include <type_traits>
#include "utils/communication.hpp"

#if defined(MFEM_USE_STRUMPACK)

namespace palace
{

namespace
{

strumpack::CompressionType GetCompressionType(SparseCompression type)
{
  switch (type)
  {
    case SparseCompression::HSS:
      return strumpack::CompressionType::HSS;
    case SparseCompression::BLR:
      return strumpack::CompressionType::BLR;
    case SparseCompression::HODLR:
      return strumpack::CompressionType::HODLR;
    case SparseCompression::ZFP:
      return strumpack::CompressionType::LOSSY;
    case SparseCompression::BLR_HODLR:
      return strumpack::CompressionType::BLR_HODLR;
      break;
    case SparseCompression::ZFP_BLR_HODLR:
      return strumpack::CompressionType::ZFP_BLR_HODLR;
      break;
    case SparseCompression::NONE:
      return strumpack::CompressionType::NONE;
  }
  return strumpack::CompressionType::NONE;  // For compiler warning
}

}  // namespace

template <typename StrumpackSolverType>
StrumpackSolverBase<StrumpackSolverType>::StrumpackSolverBase(
    MPI_Comm comm, SymbolicFactorization reorder, SparseCompression compression,
    double lr_tol, int butterfly_l, int lossy_prec, bool reorder_reuse, int print)
  : StrumpackSolverType(comm), comm(comm)
{
  // Configure the solver.
  this->SetPrintFactorStatistics(print > 1);
  this->SetPrintSolveStatistics(print > 1);
  this->SetKrylovSolver(strumpack::KrylovSolver::DIRECT);  // Always as a preconditioner or
                                                           // direct solver
  this->SetMatching(strumpack::MatchingJob::NONE);
  switch (reorder)
  {
    case SymbolicFactorization::METIS:
      this->SetReorderingStrategy(strumpack::ReorderingStrategy::METIS);
      // this->SetReorderingStrategy(strumpack::ReorderingStrategy::AND);
      break;
    case SymbolicFactorization::PARMETIS:
      this->SetReorderingStrategy(strumpack::ReorderingStrategy::PARMETIS);
      break;
    case SymbolicFactorization::SCOTCH:
      this->SetReorderingStrategy(strumpack::ReorderingStrategy::SCOTCH);
      break;
    case SymbolicFactorization::PTSCOTCH:
      this->SetReorderingStrategy(strumpack::ReorderingStrategy::PTSCOTCH);
      break;
    case SymbolicFactorization::AMD:
      this->SetReorderingStrategy(strumpack::ReorderingStrategy::AMD);
      // this->SetReorderingStrategy(strumpack::ReorderingStrategy::MMD);
      break;
    case SymbolicFactorization::RCM:
      this->SetReorderingStrategy(strumpack::ReorderingStrategy::RCM);
    case SymbolicFactorization::PORD:
    case SymbolicFactorization::DEFAULT:
      // Should have good default.
      break;
  }
  this->SetReorderingReuse(
      reorder_reuse);  // If true repeated calls use same sparsity pattern

  // Configure compression.
  this->SetCompression(GetCompressionType(compression));
  switch (compression)
  {
    case SparseCompression::ZFP:
      if (lossy_prec <= 0)
      {
        this->SetCompression(strumpack::CompressionType::LOSSLESS);
      }
      else
      {
        this->SetCompressionLossyPrecision(lossy_prec);
      }
      break;
    case SparseCompression::ZFP_BLR_HODLR:
      this->SetCompressionLossyPrecision(lossy_prec);
    case SparseCompression::HODLR:
    case SparseCompression::BLR_HODLR:
      this->SetCompressionButterflyLevels(butterfly_l);
    case SparseCompression::HSS:
    case SparseCompression::BLR:
      this->SetCompressionRelTol(lr_tol);
      break;
    case SparseCompression::NONE:
      break;
  }
  // if (mfem::Device::Allows(mfem::Backend::DEVICE_MASK))
  // {
  //   this->EnableGPU();  // XX TODO: GPU support disabled for now
  // }
  // else
  {
    this->DisableGPU();
  }

  // LAB (throwaway A/B branch): Solver.Linear.STRUMPACKButterflyLevels = 1000 + flags
  // selects candidate changes without compression: 1 FACTOR_MEMORY proportional mapping,
  // 2 cached solve tree dimensions, 4 front cost 1000 in the FACTOR_MEMORY mapping, 8
  // min-max process split, 16 small-P dynamic programming split, 32 selective inversion,
  // 64 small right-hand side solve kernels, 128 selective inversion of the triangular
  // factors (implies 32), 256 narrow TRSM for any block size, 512 / 1024 triangular
  // inversion condition limit 1e4 / unlimited, 2048 BLACS contexts on demand, 4096 memory
  // parts A-C (interval table, sparse equilibration, no permutation copy), 8192 node-shared
  // permutation.
  if (compression == SparseCompression::NONE && butterfly_l >= 1000)
  {
    const int flags = butterfly_l - 1000;
    if (flags & 2)
    {
      setenv("STRUMPACK_LAB_DIMS_CACHE", "1", 1);
    }
    if (flags & 4)
    {
      setenv("STRUMPACK_PM_FRONT_COST", "1000", 1);
    }
    if (flags & 8)
    {
      setenv("STRUMPACK_PM_MINMAX", "1", 1);
    }
    if (flags & 16)
    {
      setenv("STRUMPACK_PM_DP", "16", 1);
    }
    if (flags & 64)
    {
      setenv("STRUMPACK_LAB_KERNELS", "1", 1);
    }
    if (flags & 128)
    {
      setenv("STRUMPACK_SI_TRI", "1", 1);
    }
    if (flags & 256)
    {
      setenv("STRUMPACK_LAB_NARROW_ANY_NB", "1", 1);
    }
    if (flags & 512)
    {
      setenv("STRUMPACK_SI_TRI_MAX_COND", "1e4", 1);
    }
    if (flags & 1024)
    {
      setenv("STRUMPACK_SI_TRI_MAX_COND", "1e300", 1);
    }
    if (flags & 2048)
    {
      setenv("STRUMPACK_LAB_BLACS_LAZY", "1", 1);
    }
    if (flags & 4096)
    {
      setenv("STRUMPACK_LAB_MEM_INTERVALS", "1", 1);
      setenv("STRUMPACK_LAB_MEM_EQUIL", "1", 1);
      setenv("STRUMPACK_LAB_MEM_NOPERMCOPY", "1", 1);
    }
    if (flags & 8192)
    {
      setenv("STRUMPACK_LAB_MEM_SHARED_PERM", "1", 1);
    }
    if constexpr (std::is_same_v<StrumpackSolverType, mfem::STRUMPACKSolver>)
    {
      if (flags & 1)
      {
        this->solver_->options().set_proportional_mapping(
            strumpack::ProportionalMapping::FACTOR_MEMORY);
      }
#if defined(STRUMPACK_HAS_SELECTIVE_INVERSION)
      if (flags & (32 | 128))
      {
        this->solver_->options().enable_selective_inversion();
        setenv("STRUMPACK_SI_VERBOSE", "1", 1);
      }
#endif
    }
    if (Mpi::Root(comm))
    {
      std::printf("STRUMPACK LAB flags %d\n", flags);
    }
  }
}

template <typename StrumpackSolverType>
void StrumpackSolverBase<StrumpackSolverType>::SetOperator(const Operator &op)
{
  // Convert the input operator to a distributed STRUMPACK matrix. This is very similar to
  // the MFEM's STRUMPACKRowLocMatrix from a HypreParMatrix but avoids using the
  // communicator from the Hypre matrix in the case that the solver is constructed on a
  // different communicator. The sparsity pattern is not declared symmetric: hypre's
  // parallel matrix products can store explicit zeros without their transposes on small
  // partitions, and STRUMPACK then orders and builds its elimination tree from an
  // unsymmetric graph, which fails or hangs. STRUMPACK symmetrizes the pattern instead,
  // which adds no entries when it is already symmetric.
  const auto *hA = dynamic_cast<const mfem::HypreParMatrix *>(&op);
  MFEM_VERIFY(hA && hA->GetGlobalNumRows() == hA->GetGlobalNumCols(),
              "StrumpackSolver requires a square HypreParMatrix operator!");
  auto *parcsr = (hypre_ParCSRMatrix *)const_cast<mfem::HypreParMatrix &>(*hA);
  hypre_CSRMatrix *csr = hypre_MergeDiagAndOffd(parcsr);
  hypre_CSRMatrixMigrate(csr, HYPRE_MEMORY_HOST);

  // Create the STRUMPACKRowLocMatrix by taking the internal data from a hypre_CSRMatrix.
  HYPRE_BigInt glob_n = hypre_ParCSRMatrixGlobalNumRows(parcsr);
  HYPRE_BigInt first_row = hypre_ParCSRMatrixFirstRowIndex(parcsr);
  HYPRE_Int n_loc = hypre_CSRMatrixNumRows(csr);
  HYPRE_Int *I = hypre_CSRMatrixI(csr);
  HYPRE_BigInt *J = hypre_CSRMatrixBigJ(csr);
  double *data = hypre_CSRMatrixData(csr);

  // Safe to delete the matrix since STRUMPACK copies it on input. Also clean up the Hypre
  // data structure once we are done with it.
#if !defined(HYPRE_BIGINT)
  mfem::STRUMPACKRowLocMatrix A(comm, n_loc, first_row, glob_n, glob_n, I, J, data, false);
#else
  int n_loc_int = static_cast<int>(n_loc);
  MFEM_ASSERT(n_loc == (HYPRE_Int)n_loc_int,
              "Overflow error for local sparse matrix size!");
  mfem::Array<int> II(n_loc_int + 1);
  for (int i = 0; i <= n_loc_int; i++)
  {
    II[i] = static_cast<int>(I[i]);
    MFEM_ASSERT(I[i] == (HYPRE_Int)II[i], "Overflow error for local sparse matrix index!");
  }
  mfem::STRUMPACKRowLocMatrix A(comm, n_loc_int, first_row, glob_n, glob_n, II.HostRead(),
                                J, data, false);
#endif
  StrumpackSolverType::SetOperator(A);
  hypre_CSRMatrixDestroy(csr);
}

template <typename StrumpackSolverType>
void StrumpackSolverBase<StrumpackSolverType>::SetReorderReuse(bool reorder_reuse)
{
  StrumpackSolverType::SetReorderingReuse(
      reorder_reuse);  // If true repeated calls use same sparsity pattern
}

template class StrumpackSolverBase<mfem::STRUMPACKSolver>;
template class StrumpackSolverBase<mfem::STRUMPACKMixedPrecisionSolver>;

}  // namespace palace

#endif
