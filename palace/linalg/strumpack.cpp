// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "strumpack.hpp"

#if defined(MFEM_USE_STRUMPACK)

#include <complex>
#include <type_traits>

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

template <typename ScalarType>
void ConfigureStrumpackCompression(strumpack::SPOptions<ScalarType> &opts,
                                   SparseCompression compression, double lr_tol,
                                   int butterfly_l, int lossy_prec)
{
  opts.set_compression(GetCompressionType(compression));
  switch (compression)
  {
    case SparseCompression::ZFP:
      if (lossy_prec <= 0)
      {
        opts.set_compression(strumpack::CompressionType::LOSSLESS);
      }
      else
      {
        opts.set_lossy_precision(lossy_prec);
      }
      break;
    case SparseCompression::ZFP_BLR_HODLR:
      opts.set_lossy_precision(lossy_prec);
    case SparseCompression::HODLR:
    case SparseCompression::BLR_HODLR:
      opts.HODLR_options().set_butterfly_levels(butterfly_l);
    case SparseCompression::HSS:
    case SparseCompression::BLR:
      opts.set_compression_rel_tol(lr_tol);
      break;
    case SparseCompression::NONE:
      break;
  }
}

template void ConfigureStrumpackCompression(strumpack::SPOptions<double> &,
                                            SparseCompression, double, int, int);
template void ConfigureStrumpackCompression(strumpack::SPOptions<std::complex<double>> &,
                                            SparseCompression, double, int, int);

template <typename ScalarType>
void ConfigureStrumpackReordering(strumpack::SPOptions<ScalarType> &opts,
                                  SymbolicFactorization reorder)
{
  switch (reorder)
  {
    case SymbolicFactorization::METIS:
      opts.set_reordering_method(strumpack::ReorderingStrategy::METIS);
      break;
    case SymbolicFactorization::PARMETIS:
      opts.set_reordering_method(strumpack::ReorderingStrategy::PARMETIS);
      break;
    case SymbolicFactorization::SCOTCH:
      opts.set_reordering_method(strumpack::ReorderingStrategy::SCOTCH);
      break;
    case SymbolicFactorization::PTSCOTCH:
      opts.set_reordering_method(strumpack::ReorderingStrategy::PTSCOTCH);
      break;
    case SymbolicFactorization::AMD:
      opts.set_reordering_method(strumpack::ReorderingStrategy::AMD);
      break;
    case SymbolicFactorization::RCM:
      opts.set_reordering_method(strumpack::ReorderingStrategy::RCM);
      break;
    case SymbolicFactorization::PORD:
    case SymbolicFactorization::DEFAULT:
      opts.set_reordering_method(GetStrumpackDefaultReordering());
      break;
  }
}

template void ConfigureStrumpackReordering(strumpack::SPOptions<double> &,
                                           SymbolicFactorization);
template void ConfigureStrumpackReordering(strumpack::SPOptions<std::complex<double>> &,
                                           SymbolicFactorization);

template <typename StrumpackSolverType>
StrumpackSolverBase<StrumpackSolverType>::StrumpackSolverBase(
    MPI_Comm comm, SymbolicFactorization reorder, SparseCompression compression,
    double lr_tol, int butterfly_l, int lossy_prec, bool reorder_reuse,
    bool selective_inversion, int print, bool factor_size_mapping)
  : StrumpackSolverType(comm), comm(comm)
{
  // Configure the solver.
  this->SetPrintFactorStatistics(print > 1);
  this->SetPrintSolveStatistics(print > 1);
  this->SetKrylovSolver(strumpack::KrylovSolver::DIRECT);  // Always as a preconditioner or
                                                           // direct solver
  this->SetMatching(strumpack::MatchingJob::NONE);
  this->SetReorderingReuse(
      reorder_reuse);  // If true repeated calls use same sparsity pattern

  // Configure compression.
  ConfigureStrumpackCompression(this->solver_->options(), compression, lr_tol, butterfly_l,
                                lossy_prec);

  // Set the column ordering and store the permutation once per node. The mixed precision
  // solver computes the ordering in its single precision solver, which has its own
  // options: the column ordering is not set there, so that it is always the default of
  // STRUMPACK, METIS.
  if constexpr (std::is_same_v<StrumpackSolverType, mfem::STRUMPACKMixedPrecisionSolver>)
  {
    EnableStrumpackSharedPermutation(this->solver_->solver().options());
  }
  else
  {
    ConfigureStrumpackReordering(this->solver_->options(), reorder);
    EnableStrumpackSharedPermutation(this->solver_->options());

    // Map the elimination tree by factor size (the mixed precision solver keeps the
    // STRUMPACK default).
    if (factor_size_mapping)
    {
      SetStrumpackFactorSizeMapping(this->solver_->options());
    }
  }

  // Selective inversion of the distributed fronts. Never for the mixed precision solver:
  // the inverted factors of its single precision factorization would not be accurate
  // enough with the condition number limit which Palace sets.
  if constexpr (std::is_same_v<StrumpackSolverType, mfem::STRUMPACKSolver>)
  {
    if (selective_inversion)
    {
      EnableStrumpackSelectiveInversion(this->solver_->options());
    }
  }
  // if (mfem::Device::Allows(mfem::Backend::DEVICE_MASK))
  // {
  //   this->EnableGPU();  // XX TODO: GPU support disabled for now
  // }
  // else
  {
    this->DisableGPU();
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

template <typename StrumpackSolverType>
strumpack::ReorderingStrategy
StrumpackSolverBase<StrumpackSolverType>::GetReorderingStrategy() const
{
  if constexpr (std::is_same_v<StrumpackSolverType, mfem::STRUMPACKMixedPrecisionSolver>)
  {
    return this->solver_->solver().options().reordering_method();
  }
  else
  {
    return this->solver_->options().reordering_method();
  }
}

template class StrumpackSolverBase<mfem::STRUMPACKSolver>;
template class StrumpackSolverBase<mfem::STRUMPACKMixedPrecisionSolver>;

}  // namespace palace

#endif
