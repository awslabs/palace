// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "strumpack.hpp"

#include <type_traits>

#if defined(MFEM_USE_STRUMPACK)

#include <StrumpackConfig.hpp>
#include "utils/communication.hpp"

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
    double lr_tol, int butterfly_l, int lossy_prec, bool reorder_reuse,
    bool selective_inversion, int print)
  : StrumpackSolverType(comm), comm(comm),
    default_reorder(reorder == SymbolicFactorization::DEFAULT)
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
      break;
    case SymbolicFactorization::DEFAULT:
      // Chosen from the matrix size in SetOperator.
      break;
    case SymbolicFactorization::PORD:
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

  if constexpr (std::is_same_v<StrumpackSolverType, mfem::STRUMPACKSolver>)
  {
    // The solve phase is applied once per Krylov iteration and is bound by the slowest
    // process: map the elimination tree to processes by factor memory (what the solve
    // reads) rather than by factorization flops.
    this->solver_->options().set_proportional_mapping(
        strumpack::ProportionalMapping::FACTOR_MEMORY);
#if defined(STRUMPACK_HAS_SELECTIVE_INVERSION)
    // Store the inverses of the distributed separator blocks so that the distributed
    // triangular solves become matrix products (latency bound at scale); ill-conditioned
    // blocks keep their triangular solves.
    if (selective_inversion)
    {
      this->solver_->options().enable_selective_inversion();
    }
#endif
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
#if defined(STRUMPACK_USE_PTSCOTCH)
  if (default_reorder)
  {
    // Parallel nested dissection (PT-Scotch distributed separators on at most 256
    // processes, METIS for the subgraph of each process; see
    // extern/patch/STRUMPACK/patch_ptscotch_ordering.diff) instead of METIS on a single
    // process, whose time grows with the matrix size while the other processes wait. Small
    // matrices on few processes, such as those of the wave port mode problems, which are
    // factorized at each exact port solve, are ordered faster by METIS.
    this->SetReorderingStrategy((Mpi::Size(comm) > 1 && glob_n >= 50000)
                                    ? strumpack::ReorderingStrategy::PTSCOTCH
                                    : strumpack::ReorderingStrategy::METIS);
  }
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
