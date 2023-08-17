// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "strumpack.hpp"

#if defined(MFEM_USE_STRUMPACK)

#include "linalg/rap.hpp"

namespace palace
{

namespace
{

strumpack::CompressionType
GetCompressionType(config::LinearSolverData::CompressionType type)
{
  switch (type)
  {
    case config::LinearSolverData::CompressionType::HSS:
      return strumpack::CompressionType::HSS;
    case config::LinearSolverData::CompressionType::BLR:
      return strumpack::CompressionType::BLR;
    case config::LinearSolverData::CompressionType::HODLR:
      return strumpack::CompressionType::HODLR;
    case config::LinearSolverData::CompressionType::ZFP:
      return strumpack::CompressionType::LOSSY;
    case config::LinearSolverData::CompressionType::BLR_HODLR:
      return strumpack::CompressionType::BLR_HODLR;
      break;
    case config::LinearSolverData::CompressionType::ZFP_BLR_HODLR:
      return strumpack::CompressionType::ZFP_BLR_HODLR;
      break;
    case config::LinearSolverData::CompressionType::NONE:
      return strumpack::CompressionType::NONE;
  }
  return strumpack::CompressionType::NONE;  // For compiler warning
}

}  // namespace

template <typename StrumpackSolverType>
StrumpackSolverBase<StrumpackSolverType>::StrumpackSolverBase(
    MPI_Comm comm, config::LinearSolverData::SymFactType reorder,
    config::LinearSolverData::CompressionType compression, double lr_tol, int butterfly_l,
    int lossy_prec, int print)
  : StrumpackSolverType(comm), comm(comm)
{
  // Configure the solver.
  this->SetPrintFactorStatistics(print > 1);
  this->SetPrintSolveStatistics(print > 1);
  this->SetKrylovSolver(strumpack::KrylovSolver::DIRECT);  // Always as a preconditioner or
                                                           // direct solver
  this->SetMatching(strumpack::MatchingJob::NONE);
  if (reorder == config::LinearSolverData::SymFactType::METIS)
  {
    this->SetReorderingStrategy(strumpack::ReorderingStrategy::METIS);
  }
  else if (reorder == config::LinearSolverData::SymFactType::PARMETIS)
  {
    this->SetReorderingStrategy(strumpack::ReorderingStrategy::PARMETIS);
  }
  else if (reorder == config::LinearSolverData::SymFactType::SCOTCH)
  {
    this->SetReorderingStrategy(strumpack::ReorderingStrategy::SCOTCH);
  }
  else if (reorder == config::LinearSolverData::SymFactType::PTSCOTCH)
  {
    this->SetReorderingStrategy(strumpack::ReorderingStrategy::PTSCOTCH);
  }
  else
  {
    // Use default
  }
  this->SetReorderingReuse(true);  // Repeated calls use same sparsity pattern

  // Configure compression.
  this->SetCompression(GetCompressionType(compression));
  switch (compression)
  {
    case config::LinearSolverData::CompressionType::ZFP:
      if (lossy_prec <= 0)
      {
        this->SetCompression(strumpack::CompressionType::LOSSLESS);
      }
      else
      {
        this->SetCompressionLossyPrecision(lossy_prec);
      }
      break;
    case config::LinearSolverData::CompressionType::ZFP_BLR_HODLR:
      this->SetCompressionLossyPrecision(lossy_prec);
    case config::LinearSolverData::CompressionType::HODLR:
    case config::LinearSolverData::CompressionType::BLR_HODLR:
      this->SetCompressionButterflyLevels(butterfly_l);
    case config::LinearSolverData::CompressionType::HSS:
    case config::LinearSolverData::CompressionType::BLR:
      this->SetCompressionRelTol(lr_tol);
      break;
    case config::LinearSolverData::CompressionType::NONE:
      break;
  }
}

template <typename StrumpackSolverType>
void StrumpackSolverBase<StrumpackSolverType>::SetOperator(const Operator &op)
{
  // Convert the input operator to a distributed STRUMPACK matrix (always assume a symmetric
  // sparsity pattern). This is very similar to the MFEM STRUMPACKRowLocMatrix from a
  // HypreParMatrix but avoids using the communicator from the Hypre matrix in the case that
  // the solver is constructed on a different communicator.
  const mfem::HypreParMatrix *hypA;
  const auto *PtAP = dynamic_cast<const ParOperator *>(&op);
  if (PtAP)
  {
    hypA = &PtAP->ParallelAssemble();
  }
  else
  {
    hypA = dynamic_cast<const mfem::HypreParMatrix *>(&op);
    MFEM_VERIFY(hypA, "StrumpackSolver requires a HypreParMatrix operator!");
  }
  auto *parcsr = (hypre_ParCSRMatrix *)const_cast<mfem::HypreParMatrix &>(*hypA);
  hypA->HostRead();
  hypre_CSRMatrix *csr = hypre_MergeDiagAndOffd(parcsr);
  hypA->HypreRead();

  // Create the STRUMPACKRowLocMatrix by taking the internal data from a hypre_CSRMatrix.
  HYPRE_Int n_loc = csr->num_rows;
  HYPRE_BigInt first_row = parcsr->first_row_index;
  HYPRE_Int *I = csr->i;
  HYPRE_BigInt *J = csr->big_j;
  double *data = csr->data;

  // Safe to delete the matrix since STRUMPACK copies it on input. Also clean up the Hypre
  // data structure once we are done with it.
#if !defined(HYPRE_BIGINT)
  mfem::STRUMPACKRowLocMatrix A(comm, n_loc, first_row, hypA->GetGlobalNumRows(),
                                hypA->GetGlobalNumCols(), I, J, data, true);
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
  mfem::STRUMPACKRowLocMatrix A(comm, n_loc_int, first_row, hypA->GetGlobalNumRows(),
                                hypA->GetGlobalNumCols(), II, J, data, true);
#endif
  StrumpackSolverType::SetOperator(A);
  hypre_CSRMatrixDestroy(csr);
}

template class StrumpackSolverBase<mfem::STRUMPACKSolver>;
template class StrumpackSolverBase<mfem::STRUMPACKMixedPrecisionSolver>;

}  // namespace palace

#endif
